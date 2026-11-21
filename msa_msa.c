#include "msa_msa.h"
#include "msa_traceback.h"

/**
 * Performs multiple sequnce alignment on an inputted array of sequences
 * returns the final alignment length for all sequences
 */
int msa(Gene_types** genes, int* gene_sizes, int num_genes, Gene_types* aligned_genes){
    //Calloc for here instead of malloc, to initialize with 0s.  Diags aren't filled
    //This really only matters when looking at it with print, but best to not get confused
    float *dist_matrix = calloc((size_t)num_genes * (size_t)num_genes, sizeof(float));
    if (dist_matrix == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    //calculate total size of all genes and find the largest size out of the genes
    int total_size = 0;
    int max_gene_size = 0;
    for(int i=0; i<num_genes; i++){
        total_size += gene_sizes[i];
        if (gene_sizes[i] > max_gene_size){
            max_gene_size = gene_sizes[i];
        }
    }

    //put genes in an 1D array to use on GPU
    int cur_offset = 0;
    Gene_types *genes_1D = malloc(sizeof(Gene_types) * (size_t)total_size);
    int offsets[num_genes];
    if (genes_1D == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }
    for (int i = 0; i < num_genes; i++){
        memcpy(genes_1D + cur_offset, genes[i], sizeof(Gene_types) * gene_sizes[i]);
        offsets[i] = cur_offset;
        cur_offset += gene_sizes[i];
    }

    //create score matrices for alignment
    int total_combinations = (num_genes * (num_genes-1))/2;
    int max_matrix_size = (max_gene_size+1)*(max_gene_size+1);
    int* all_scores = malloc(sizeof(int) * total_combinations * max_matrix_size);

    //lookup table for i and j pairs
    int* i_k = malloc(sizeof(int) * total_combinations);
    int* j_k = malloc(sizeof(int) * total_combinations);

    int k = 0;

    for(int i = 0; i<num_genes; i++){
        for(int j = i+1; j<num_genes; j++){
            i_k[k] = i;
            j_k[k] = j;

            k++;
        }
    }

    //do pairwise alignments and fill in the distance matrix
    //this can be done in parallel since each pairwise alignment is independent of the others
    #ifdef OPEN_ACC_ON
    #pragma acc data copy(dist_matrix[0:(size_t)num_genes * num_genes]) \
                copyin(genes_1D[0:total_size], gene_sizes[0:num_genes], offsets[0:num_genes], \
                i_k[0:total_combinations], j_k[0:total_combinations]) \
                create(all_scores[0:(size_t)total_combinations*max_matrix_size])
    #endif
    {
        //Create a gang for every combination of genes (total_combinations)
        //vector_length sets the vectors in the gangs to 128
        #ifdef OPEN_ACC_ON
        #pragma acc parallel loop gang num_gangs(total_combinations) vector_length(128)
        #endif
        //iterate through all combinations (each one has it's own gang)
        for(int k = 0; k < total_combinations; k++){
            int i = i_k[k];
            int j = j_k[k];

            //Set up local variables for this specific combination
            const Gene_types *gene_A = &genes_1D[offsets[i]];
            const Gene_types *gene_B = &genes_1D[offsets[j]];
            const int A_size = gene_sizes[i];
            const int B_size = gene_sizes[j];
            const int num_rows = A_size + 1;
            const int num_cols = B_size + 1;
            int *scores = &all_scores[k * max_matrix_size];

            const int total_waves = num_rows + num_cols - 1;

            //Fill up rows with negative values based on gap penalty
            //We can do it parallel because there are no dependencies on last loop iteration!
            #ifdef OPEN_ACC_ON
            #pragma acc loop vector
            #endif
            for (int row = 0; row < num_rows; row++){
                scores[row * num_cols] = row * GAP_PENALTY_VAL;
            }

            //Fill up rows with negative values based on gap penalty
            //We can do it parallel!
            #ifdef OPEN_ACC_ON
            #pragma acc loop vector
            #endif
            for (int col = 0; col < num_cols; col++){
                scores[col] = col * GAP_PENALTY_VAL;
            }

            //Unlike PSA, waves are confined to a gang and run across the 128 vector langs
            for (int wave = 2; wave < total_waves; wave++) {

                // Determine bounds for the current diagonal wave
                int i_start = (wave < num_cols) ? 1 : wave - num_cols + 1;
                int i_end = (wave < num_rows) ? wave - 1 : num_rows - 1;

                #ifdef OPEN_ACC_ON
                #pragma acc loop vector
                #endif
                for (int row = i_start; row <= i_end; row++) {
                    score_matrix_cell(gene_A, gene_B, num_cols, row, wave - row, scores);
                }
            }
        }

        // Parallel traceback, the traceback itself is not parallel, but parallel across all combinations
        #ifdef OPEN_ACC_ON
        #pragma acc parallel loop
        #endif
        for(int k = 0; k < total_combinations; k++){
            int i = i_k[k];
            int j = j_k[k];

            float identity = traceback_for_identity(&genes_1D[offsets[i]], &genes_1D[offsets[j]], gene_sizes[i], gene_sizes[j], &all_scores[k*max_matrix_size]);

            float dist = 1 - identity;
            dist_matrix[i*num_genes + j] = dist;
            dist_matrix[j*num_genes + i] = dist; // distance matrix is symmetric
        }
    }

    free(all_scores);
    free(i_k);
    free(j_k);

    // print distance matrix
    #ifdef PRINT_VAR_OUTPUT
        printf("Distance matrix:\n");
        for (int i = 0; i < num_genes; i++){
            for(int j = 0; j < num_genes; j++){
                printf("%f ", dist_matrix[i * num_genes + j]);
            }
            printf("\n");
        }
    #endif

    //neighbor joining to find order to align the sequences
    int joins[num_genes-1][2]; // store the pairs of joined nodes
    neighbor_joining(num_genes, dist_matrix, joins);

    //align sequences in order determined bt neighbor joining (0-num_genes-1 is individual sequences, anything above is joined node corresponding
    // to the joins array order (eg 3 sequences, 0, 1, and 2 are individual sequences, 3 is the first element of joins (found by joins[3-(num_genes)])))

    //initialize alignments for each sequence

    Gene_types *aligned_sequences = malloc(sizeof(Gene_types) * (size_t)num_genes * (size_t)total_size);
    int sequence_lengths[num_genes]; // store the lengths of the aligned sequences

    // initialize storage for aligned sequences of each join
    int *joined_seq_ids = malloc(sizeof(int) * (size_t)num_genes * (size_t)num_genes); // store the indices of the sequences that were joined
    int joined_num_genes[num_genes-1]; // store the number of sequences in each joined node

    if (aligned_sequences == NULL || joined_seq_ids == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    // go through the joins array and align sequence in that order, storing the aligned sequences in the aligned_sequences array
    for(int i=0; i<(num_genes-1); i++){
        #ifdef PRINT_VAR_OUTPUT
            printf("Join %d\n", i+1);
        #endif
        //counter for number of sequences in joined node
        int num_seq_A = 0;
        int num_seq_B = 0;

        int seqA_index = joins[i][0];
        int seqB_index = joins[i][1];

        int seqA_size = 0;
        int seqB_size = 0;

        Gene_types *seqA = malloc(sizeof(Gene_types) * total_size * num_genes);
        Gene_types *seqB = malloc(sizeof(Gene_types) * total_size * num_genes);

        //gets sequence A and its information (size, number of sequences, and sequence ids) based on whether it is an individual sequence or a joined node
        if (seqA_index < num_genes) {
            // seqA is an individual sequence
            num_seq_A = 1;
            seqA_size = gene_sizes[seqA_index];
            memcpy(seqA, genes[seqA_index], sizeof(Gene_types) * seqA_size);
        } else {
            // seqA is a joined node
            int join_index_A = seqA_index - num_genes;
            num_seq_A = joined_num_genes[join_index_A];

            //combines all of the joined sequences of seqA into a single sequence for alignment
            for (int j = 0; j < num_seq_A; j++) {
                int seq_id = joined_seq_ids[join_index_A * num_genes + j];
                seqA_size = sequence_lengths[seq_id];

                memcpy(seqA + j*seqA_size, aligned_sequences + seq_id * total_size, sizeof(Gene_types) * seqA_size);
            }
        }
        //print sequence a upon loading
        // printf("Alignment A:\n");
        // print_aligned_genes(seqA, num_seq_A, seqA_size);

        //gets sequence B and its information (size, number of sequences, and sequence ids) based on whether it is an individual sequence or a joined node
        if (seqB_index < num_genes) {
            // seqB is an individual sequence
            num_seq_B = 1;
            seqB_size = gene_sizes[seqB_index];
            memcpy(seqB, genes[seqB_index], sizeof(Gene_types) * seqB_size);
        } else {
            // seqB is a joined node
            int join_index_B = seqB_index - num_genes;
            num_seq_B = joined_num_genes[join_index_B];

            //combines all of the joined sequences of seqB into a single sequence for alignment
            for (int j = 0; j < num_seq_B; j++) {
                int seq_id = joined_seq_ids[join_index_B * num_genes + j];
                seqB_size = sequence_lengths[seq_id];

                memcpy(seqB + j*seqB_size, aligned_sequences + seq_id * total_size, sizeof(Gene_types) * seqB_size);
            }
        }
        //print sequence a upon loading
        // printf("Alignment B:\n");
        // print_aligned_genes(seqB, num_seq_B, seqB_size);

        //align the two sequences
        float *score_matrix = malloc(sizeof(float) * (size_t)(seqA_size + 1) * (size_t)(seqB_size + 1));
        // will be stored backwards (last element is at index 0)
        Gene_types *aligned_A = malloc(sizeof(Gene_types) * (size_t)num_seq_A * (size_t)(seqA_size + seqB_size));
        Gene_types *aligned_B = malloc(sizeof(Gene_types) * (size_t)num_seq_B * (size_t)(seqA_size + seqB_size));

        if (score_matrix == NULL || aligned_A == NULL || aligned_B == NULL) {
            printf("Memory allocation failed!\n");
            abort();
        }

        altered_score_matrix_for_msa(seqA, seqB, seqA_size, seqB_size, num_seq_A, num_seq_B, score_matrix);

        int alignment_length = 0;

        float identity_score = altered_traceback_for_msa(seqA, seqB, seqA_size, seqB_size, num_seq_A, num_seq_B, score_matrix, aligned_A, aligned_B, &alignment_length);

        //free seqA and seqB
        free(seqA);
        free(seqB);

        //store the aligned sequences in the aligned_sequences array
        //also store the sequence ids for the join location in the joined_seq_ids array

        int join_index = i; // index for the joined node
        joined_num_genes[join_index] = num_seq_A + num_seq_B; // store the number of sequences in the joined node

        int gene_counter = 0; // counter for the number of sequences in the joined node

        //seq A genes
        if (seqA_index < num_genes) {
            // seqA is an individual sequence
            //reverse aligned_A and store in aligned_sequences
            flip_sequence(aligned_A, alignment_length);
            
            memcpy(aligned_sequences + seqA_index * total_size, aligned_A, sizeof(Gene_types) * alignment_length);
            //store new length of sequence in sequence lengths
            sequence_lengths[seqA_index] = alignment_length;

            // store the sequence id in the joined_seq_ids array
            joined_seq_ids[join_index * num_genes + gene_counter] = seqA_index;
            gene_counter++;

        } else {
            // seqA is a joined node
            int join_index_A = seqA_index - num_genes;

            for (int j = 0; j < num_seq_A; j++) {
                int seq_id = joined_seq_ids[join_index_A * num_genes + j];
                //reverse the sequence
                flip_sequence(aligned_A + (j*(seqA_size+seqB_size)), alignment_length);

                memcpy(aligned_sequences + seq_id * total_size, aligned_A + (j*(seqA_size+seqB_size)), sizeof(Gene_types) * alignment_length);
                // store new length of sequence in sequence lengths
                sequence_lengths[seq_id] = alignment_length;

                // store the sequence id in the joined_seq_ids array
                joined_seq_ids[join_index * num_genes + gene_counter] = seq_id;
                gene_counter++;
            }
        }
        
        //seq B genes
        if (seqB_index < num_genes) {
            // seqB is an individual sequence
            //reverse aligned_B and store in aligned_sequences
            flip_sequence(aligned_B, alignment_length);
            
            memcpy(aligned_sequences + seqB_index * total_size, aligned_B, sizeof(Gene_types) * alignment_length);
            //store new length of sequence in sequence lengths
            sequence_lengths[seqB_index] = alignment_length;

            // store the sequence id in the joined_seq_ids array
            joined_seq_ids[join_index * num_genes + gene_counter] = seqB_index;
            gene_counter++;
        } else {
            // seqB is a joined node
            int join_index_B = seqB_index - num_genes;

            for (int j = 0; j < num_seq_B; j++) {
                int seq_id = joined_seq_ids[join_index_B * num_genes + j];
                //reverse the sequence
                flip_sequence(aligned_B + (j*(seqA_size+seqB_size)), alignment_length);

                memcpy(aligned_sequences + seq_id * total_size, aligned_B + (j*(seqA_size+seqB_size)), sizeof(Gene_types) * alignment_length);
                // store new length of sequence in sequence lengths
                sequence_lengths[seq_id] = alignment_length;

                // store the sequence id in the joined_seq_ids array
                joined_seq_ids[join_index * num_genes + gene_counter] = seq_id;
                gene_counter++;
            }
        }

        free(score_matrix);
        free(aligned_A);
        free(aligned_B);
    }

    //done with msa, copy the aligned sequences to the output array
    int alignment_length = sequence_lengths[0]; // all sequences should have the same length after MSA

    for(int i=0; i<num_genes; i++){
        memcpy(aligned_genes + i * alignment_length, aligned_sequences + i * total_size, sizeof(Gene_types) * alignment_length);
    }

    free(dist_matrix);
    free(genes_1D);
    free(aligned_sequences);
    free(joined_seq_ids);

    return alignment_length;
}

/**
 * Performs neighbor joining to construct a phylogenetic tree from a distance matrix
 */
void neighbor_joining(int num_genes, float *distance_matrix, int joins[num_genes-1][2]) {
    // initiate values
    int L = num_genes;

    // store the indices of the current nodes on the distance matrix
    int *indices = malloc(sizeof(int) * (size_t)num_genes);
    float *net_divergence = malloc(sizeof(float) * (size_t)num_genes);
    float *adjusted_d = malloc(sizeof(float) * (size_t)num_genes * (size_t)num_genes);

    if (indices == NULL || net_divergence == NULL || adjusted_d == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    for (int i = 0; i < L; i++){
        indices[i] = i;
    }

    int new_node_index = num_genes; // index for the new node created by joining two nodes

    int num_joins = 0; // counter for the number of joins made

    //loop until all nodes have been joined into a single tree
    while (L > 2){

        //calculate net divergence
        for (int i = 0; i < L; i++){
            net_divergence[i] = 0.0;
            for (int j = 0; j < L; j++){
                if (i != j){
                    net_divergence[i] += distance_matrix[i * num_genes + j];
                }
            }
            net_divergence[i] /= (L - 2);
            // printf("%f\n", net_divergence[i]);
        }

        //adjust distance matrix and find pair of nodes with smallest distance
        int min_i = 0;
        int min_j = 1;
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                if (i != j){
                    //adjust distance matrix using net divergences
                    adjusted_d[i * L + j] = distance_matrix[i*num_genes + j] - (net_divergence[i] + net_divergence[j]);
                    adjusted_d[j * L + i] = adjusted_d[i * L + j];

                    //find pair with smallest adjusted distance
                    if (adjusted_d[i * L + j] < adjusted_d[min_i * L + min_j]){
                    min_i = i;
                    min_j = j;
                    }
                }
            }
        }

        // store joined pair
        joins[num_joins][0] = indices[min_i];
        joins[num_joins][1] = indices[min_j];

        num_joins++;

        //calculate distances of nodes to new internal node
        //this is not used unless an actual tree is to be built, since we only need ordering, this is commented out
        // float dist_iz = (distance_matrix[min_i * num_genes + min_j] + net_divergence[min_i] - net_divergence[min_j]) / 2.0;
        // float dist_jz = (distance_matrix[min_i * num_genes + min_j] + net_divergence[min_j] - net_divergence[min_i]) / 2.0;

        // create new node and update distance matrix
        //replace min_i location on the distance matrix with distances from the new node instead
        indices[min_i] = new_node_index;
        new_node_index++;
        for (int k = 0; k < L; k++){
            if (k != min_i && k != min_j){
                //calculates new distance of other nodes from the combined node
                distance_matrix[min_i * num_genes + k] = (distance_matrix[min_i * num_genes + k] + distance_matrix[min_j * num_genes + k] - distance_matrix[min_i * num_genes + min_j]) / 2.0;
                distance_matrix[k * num_genes + min_i] = distance_matrix[min_i * num_genes + k];
            }
        }
        

        //replace min_j with last element of matrix to reduce the size to reflect the combining of node i and j
        indices[min_j] = indices[L-1];
        for (int k = 0; k < L; k++){
            if((L-1) != min_j){
                distance_matrix[min_j * num_genes + k] = distance_matrix[(L-1) * num_genes + k];
                distance_matrix[k * num_genes + min_j] = distance_matrix[k * num_genes + (L-1)];
            }
        }

        //decrement L to reflect the reduced size of the distance matrix
        L--;
    }

    //termination
    //store the last two nodes to be combined in the last location of the joins array
    joins[num_joins][0] = indices[0];
    joins[num_joins][1] = indices[1];

    free(indices);
    free(net_divergence);
    free(adjusted_d);
}

/**
 * Creates a score matrix that allows for comparison of multiple sequences at once by averaging the match scores for each pair of sequences in the two groups being compared.
 * This function is similar to the score_matrix function, but it takes in two groups of sequences as opposed to just one sequence in each group.
 * additional variable num_seq_A and num_seq_B tell function how many different sequences are in each group
 * also score matrix is a float to allow for average values across each sequence comparison between groups
 */
void altered_score_matrix_for_msa(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, int num_seq_A, int num_seq_B, float matrix[(A_size+1) * (B_size+1)]) {
    int num_rows = A_size + 1;
    int num_cols = B_size + 1;

    //Fill start of rows with 0, -1, -2, ect.
    for(int i = 0; i < (num_rows); i++){
        matrix[i * num_cols] = (float)(i * GAP_PENALTY_VAL);
    }

    //Fill top of columns with 0, -1, -2, ect.
    for (int j = 0; j < (num_cols); j++){
        matrix[j] = (float)(j * GAP_PENALTY_VAL);
    }

    //Little double check on the array
    // for (int i = 0; i < (num_rows); i++){
    //     for(int j = 0; j < (num_cols); j++){
    //         printf("%f ", matrix[i * num_cols + j]);
    //     }
    //     printf("\n");
    // }

    // Total diagonal waves = (GRID_SIZE + GRID_SIZE - 1)
    // Wave index starts at 2 because (0,0), (0,1), and (1,0) are already boundary-initialized
    int total_waves = num_rows + num_cols - 1;

    //Copy the matrix into location easily accessable to GPU
#ifdef OPEN_ACC_ON
#pragma acc data copy(matrix[0:(size_t)num_rows*num_cols]) \
                 copyin(gene_A[0:(num_seq_A*A_size)], gene_B[0:(num_seq_B*B_size)])
#endif
    {
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < total_waves; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < num_cols) ? 1 : wave - num_cols + 1;
            int i_end = (wave < num_rows) ? wave - 1 : num_rows - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #ifdef OPEN_ACC_ON
            #pragma acc parallel loop gang vector present(matrix, gene_A, gene_B)
            #endif
            for (int i = i_start; i <= i_end; i++) {
                int j = wave - i;

                // this grabs the values from cells top, left and top-left diag.  cell
                //cell(i-1, j), cell(i, j-1), cell(i-1, j-1)
                float top_matrix  = matrix[(i - 1) * num_cols + j];
                float left_matrix = matrix[i * num_cols + (j - 1)];
                float diag_matrix = matrix[(i - 1) * num_cols + (j - 1)];

                float match = 0.0;
                
                #ifdef OPEN_ACC_ON
                #pragma acc loop seq
                #endif
                //scores match for each individial pairing of sequences in the two groups and sums them
                for (int a = 0; a < num_seq_A; a++) {
                    #ifdef OPEN_ACC_ON
                    #pragma acc loop seq
                    #endif
                    for (int b = 0; b < num_seq_B; b++) {
                        match += (gene_A[a*A_size + i - 1] == gene_B[b*B_size + j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;
                    }
                }

                //gets average score for match score
                match /= (num_seq_A * num_seq_B);

                // proceeds with scoring as normal, but with the average match score for the two groups of sequences
                float top_score = top_matrix + GAP_PENALTY_VAL;
                float left_score = left_matrix + GAP_PENALTY_VAL;
                float diag_score = diag_matrix + match;

                if ((top_score >= left_score) && (top_score >= diag_score)){
                    matrix[i * num_cols + j] = top_score;
                }
                else if((left_score >= top_score) && (left_score >= diag_score)){
                    matrix[i * num_cols + j] = left_score;
                }
                else{
                    matrix[i * num_cols + j] = diag_score;
                }
            }
        }
    }

    // //Little double check on the array
    // for (int i = 0; i < num_rows; i++){
    //     for(int j = 0; j < num_cols; j++){
    //         printf("%f ", matrix[i * num_cols + j]);
    //     }
    //     printf("\n");
    // }
}

/**
 * traceback function for msa alignment
 * similar to the traceback algorithm for pairwise, except it allows for two groups of sequences to be aligned, as opposed to just two individual sequences
 * additional variable num_seq_A and num_seq_B tell function how many different sequences are in each group
 */
float altered_traceback_for_msa(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, 
    int num_seq_A, int num_seq_B, float score_matrix[(A_size+1) * (B_size+1)], Gene_types* aligned_A, Gene_types* aligned_B, int *alignment_length){
    // indices for the score matrix
    int i = A_size;
    int j = B_size;

    int num_cols = B_size + 1;

    // counter for aligned sequences
    int align_index = 0;

    int identity_count = 0; // counter for identical matches
    int residues_compared = 0; // counter for total residues compared (excluding gaps)

    // main traceback loop
    while(i > 0 && j > 0){
        int score = score_matrix[i * num_cols + j];
        int score_diag = score_matrix[(i-1) * num_cols + (j-1)];
        int score_left = score_matrix[i * num_cols + (j-1)];
        int score_up = score_matrix[(i-1) * num_cols + j];

        float match = 0.0;
                
        //scores match for each individial pairing of sequences in the two groups and sums them
        for (int a = 0; a < num_seq_A; a++) {
            for (int b = 0; b < num_seq_B; b++) {
                match += (gene_A[a*A_size + i - 1] == gene_B[b*B_size + j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;
            }
        }
        //gets average score for match score
        match /= (num_seq_A * num_seq_B);

        // printf("score: %d, score_diag: %d, score_left: %d, score_up: %d\n", score, score_diag, score_left, score_up);
        // score from match/mismatch
        if (score == (score_diag + match)){
            //aligning all sequences in A
            for (int a = 0; a < num_seq_A; a++) {
                aligned_A[align_index + a*(A_size+B_size)] = gene_A[a*A_size + i - 1];
            }
            //aligning all sequences in B
            for (int b = 0; b < num_seq_B; b++) {
                aligned_B[align_index + b*(A_size+B_size)] = gene_B[b*B_size + j - 1];
            }
            i--;
            j--;

            // count identical matches and total residues compared
            if (match == MATCH_SCORE_VAL) {
                identity_count++;
            }
            residues_compared++;
            // printf("match: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene A
        else if (score == (score_left + GAP_PENALTY_VAL)){
            // adding gap to all seq in A
            for (int a = 0; a < num_seq_A; a++) {
                aligned_A[align_index + a*(A_size+B_size)] = GAP; // GAP represents a gap
            }
            // adding the corresponding nucleotide from gene B to all seq in B
            for (int b = 0; b < num_seq_B; b++) {
                aligned_B[align_index + b*(A_size+B_size)] = gene_B[b*B_size + j-1];
            }
            j--;
            // printf("gap in gene A: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene B
        else {
            // adding the corresponding nucleotide from gene A to all seq in A
            for (int a = 0; a < num_seq_A; a++) {
                aligned_A[align_index + a*(A_size+B_size)] = gene_A[a*A_size + i-1];
            }
            // adding gap to all seq in B
            for (int b = 0; b < num_seq_B; b++) {
                aligned_B[align_index + b*(A_size+B_size)] = GAP; // GAP represents a gap
            }
            i--;
            // printf("gap in gene B: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // update align index
        align_index++;
    }

    // fill in remaining gaps
    while(i > 0){
        for (int a = 0; a < num_seq_A; a++) {
            aligned_A[align_index + a*(A_size+B_size)] = gene_A[a*A_size + i-1];
        }
        for (int b = 0; b < num_seq_B; b++) {
            aligned_B[align_index + b*(A_size+B_size)] = GAP; // GAP represents a gap
        }
        align_index++;
        i--;
    }

    while (j > 0){
        for (int a = 0; a < num_seq_A; a++) {
            aligned_A[align_index + a*(A_size+B_size)] = GAP; // GAP represents a gap
        }
        for (int b = 0; b < num_seq_B; b++) {
            aligned_B[align_index + b*(A_size+B_size)] = gene_B[b*B_size + j-1];
        }
        align_index++;
        j--;
    }

    *alignment_length = align_index;

    float identity_score = (residues_compared > 0) ? (float)identity_count / (float)residues_compared : 0;

    return identity_score;
}


void print_aligned_genes(Gene_types *alignment, int num_genes, int alignment_length){
    for (int i = 0; i < num_genes; i++){
        for(int j = 0; j < alignment_length; j++){
            printf("%c", num_to_nucleotide(alignment[i*alignment_length + j]));
        }
        printf("\n");
    }
}