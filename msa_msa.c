#include "msa_msa.h"

/**
 * Performs multiple sequnce alignment on an inputted array of sequences
 * returns the final alignment length for all sequences
 */
int msa(int** genes, int* gene_sizes, int num_genes, int gap_penalty, int mismatch_penalty, int match_score, int* aligned_genes){
    float dist_matrix[num_genes * num_genes];

    //do pairwise alignments and fill in the distance matrix
    //this can be done in parallel since each pairwise alignment is independent of the others (not implemented yet)
    for(int i=0; i<num_genes; i++){
        for(int j=i+1; j<num_genes; j++){
            int aligned_A[gene_sizes[i] + gene_sizes[j]];
            int aligned_B[gene_sizes[i] + gene_sizes[j]];
            int alignment_length = 0;
            float dist = needleman_wunsch(genes[i], genes[j], gap_penalty, mismatch_penalty, match_score, gene_sizes[i], gene_sizes[j], aligned_A, aligned_B, alignment_length);
            dist_matrix[i*num_genes + j] = dist;
            dist_matrix[j*num_genes + i] = dist; // distance matrix is symmetric
        }
    }

    //neighbor joining to find order to align the sequences
    int joins[num_genes-1][2]; // store the pairs of joined nodes
    neighbor_joining(num_genes, dist_matrix, joins);

    //align sequences in order determined bt neighbor joining (0-num_genes-1 is individual sequences, anything above is joined node corresponding
    // to the joins array order (eg 3 sequences, 0, 1, and 2 are individual sequences, 3 is the first element of joins (found by joins[3-(num_genes)])))

    //initialize alignments for each sequence
    
    int total_size = 0;
    for(int i=0; i<num_genes; i++){
        total_size += gene_sizes[i];
    }

    int aligned_sequences[num_genes*total_size];
    int sequence_lengths[num_genes]; // store the lengths of the aligned sequences

    // initialize storage for aligned sequences of each join
    int joined_seq_ids[num_genes-1][num_genes]; // store the indices of the sequences that were joined
    int joined_num_genes[num_genes-1]; // store the number of sequences in each joined node

    // go through the joins array and align sequence in that order, storing the aligned sequences in the aligned_sequences array
    for(int i=0; i<(num_genes-1); i++){
        //counter for number of sequences in joined node
        int num_seq_A = 0;
        int num_seq_B = 0;

        int seqA_index = joins[i][0];
        int seqB_index = joins[i][1];

        int seqA_size = 0;
        int seqB_size = 0;

        int *seqA = malloc(sizeof(int) * total_size);
        int *seqB = malloc(sizeof(int) * total_size);

        //gets sequence A and its information (size, number of sequences, and sequence ids) based on whether it is an individual sequence or a joined node
        if (seqA_index < num_genes) {
            // seqA is an individual sequence
            num_seq_A = 1;
            seqA_size = gene_sizes[seqA_index];
            memcpy(seqA, genes[seqA_index], sizeof(int) * seqA_size);
        } else {
            // seqA is a joined node
            int join_index_A = seqA_index - num_genes;
            num_seq_A = joined_num_genes[join_index_A];

            //combines all of the joined sequences of seqA into a single sequence for alignment
            for (int j = 0; j < num_seq_A; j++) {
                int seq_id = joined_seq_ids[join_index_A][j];
                seqA_size = sequence_lengths[seq_id];

                memcpy(seqA + j*seqA_size, aligned_sequences + seq_id * total_size, sizeof(int) * seqA_size);
            }
        }

        //gets sequence B and its information (size, number of sequences, and sequence ids) based on whether it is an individual sequence or a joined node
        if (seqB_index < num_genes) {
            // seqB is an individual sequence
            num_seq_B = 1;
            seqB_size = gene_sizes[seqB_index];
            memcpy(seqB, genes[seqB_index], sizeof(int) * seqB_size);
        } else {
            // seqB is a joined node
            int join_index_B = seqB_index - num_genes;
            num_seq_B = joined_num_genes[join_index_B];

            //combines all of the joined sequences of seqB into a single sequence for alignment
            for (int j = 0; j < num_seq_B; j++) {
                int seq_id = joined_seq_ids[join_index_B][j];
                seqB_size = sequence_lengths[seq_id];

                memcpy(seqB + j*seqB_size, aligned_sequences + seq_id * total_size, sizeof(int) * seqB_size);
            }
        }

        //align the two sequences
        float score_matrix[(seqA_size+1) * (seqB_size+1)];
        altered_score_matrix_for_msa(seqA, seqB, gap_penalty, mismatch_penalty, match_score, seqA_size, seqB_size, num_seq_A, num_seq_B, score_matrix);

        // will be stored backwards (last elements is at index 0)
        int aligned_A[num_seq_A * (seqA_size + seqB_size)];
        int aligned_B[num_seq_B * (seqA_size + seqB_size)];

        int alignment_length = 0;

        float identity_score = altered_traceback_for_msa(seqA, seqB, gap_penalty, mismatch_penalty, match_score, seqA_size, seqB_size, num_seq_A, num_seq_B, score_matrix, aligned_A, aligned_B, &alignment_length);

        //store the aligned sequences in the aligned_sequences array
        //also store the sequence ids for the join location in the joined_seq_ids array

        int join_index = i; // index for the joined node
        joined_num_genes[join_index] = num_seq_A + num_seq_B; // store the number of sequences in the joined node

        int gene_counter = 0; // counter for the number of sequences in the joined node

        //seq A genes
        if (seqA_index < num_genes) {
            // seqA is an individual sequence
            //reverse aligned_A and store in aligned_sequences
            int reversed_aligned_A[alignment_length];
            for (int k = 0; k < alignment_length; k++) {
                reversed_aligned_A[k] = aligned_A[alignment_length - 1 - k];
            }
            
            memcpy(aligned_sequences + seqA_index * total_size, reversed_aligned_A, sizeof(int) * alignment_length);
            //store new length of sequence in sequence lengths
            sequence_lengths[seqA_index] = alignment_length;

            // store the sequence id in the joined_seq_ids array
            joined_seq_ids[join_index][gene_counter] = seqA_index;
            gene_counter++;

        } else {
            // seqA is a joined node
            int join_index_A = seqA_index - num_genes;

            for (int j = 0; j < num_seq_A; j++) {
                int seq_id = joined_seq_ids[join_index_A][j];
                //reverse the sequence
                int reversed_aligned_A[alignment_length];
                for (int k = 0; k < alignment_length; k++) {
                    reversed_aligned_A[k] = aligned_A[j*(seqA_size+seqB_size) + alignment_length - 1 - k];
                }
                memcpy(aligned_sequences + seq_id * total_size, reversed_aligned_A, sizeof(int) * alignment_length);
                // store new length of sequence in sequence lengths
                sequence_lengths[seq_id] = alignment_length;

                // store the sequence id in the joined_seq_ids array
                joined_seq_ids[join_index][gene_counter] = seq_id;
                gene_counter++;
            }
        }
        
        //seq B genes
        if (seqB_index < num_genes) {
            // seqB is an individual sequence
            //reverse aligned_B and store in aligned_sequences
            int reversed_aligned_B[alignment_length];
            for (int k = 0; k < alignment_length; k++) {
                reversed_aligned_B[k] = aligned_B[alignment_length - 1 - k];
            }
            
            memcpy(aligned_sequences + seqB_index * total_size, reversed_aligned_B, sizeof(int) * alignment_length);
            //store new length of sequence in sequence lengths
            sequence_lengths[seqB_index] = alignment_length;

            // store the sequence id in the joined_seq_ids array
            joined_seq_ids[join_index][gene_counter] = seqB_index;
            gene_counter++;
        } else {
            // seqB is a joined node
            int join_index_B = seqB_index - num_genes;

            for (int j = 0; j < num_seq_B; j++) {
                int seq_id = joined_seq_ids[join_index_B][j];
                //reverse the sequence
                int reversed_aligned_B[alignment_length];
                for (int k = 0; k < alignment_length; k++) {
                    reversed_aligned_B[k] = aligned_B[j*(seqA_size+seqB_size) + alignment_length - 1 - k];
                }
                memcpy(aligned_sequences + seq_id * total_size, reversed_aligned_B, sizeof(int) * alignment_length);
                // store new length of sequence in sequence lengths
                sequence_lengths[seq_id] = alignment_length;

                // store the sequence id in the joined_seq_ids array
                joined_seq_ids[join_index][gene_counter] = seq_id;
                gene_counter++;
            }
        }
    }

    //done with msa, copy the aligned sequences to the output array
    int alignment_length = sequence_lengths[0]; // all sequences should have the same length after MSA

    for(int i=0; i<num_genes; i++){
        memcpy(aligned_genes + i * alignment_length, aligned_sequences + i * total_size, sizeof(int) * alignment_length);
    }

    return alignment_length;
}

/**
 * Performs neighbor joining to construct a phylogenetic tree from a distance matrix
 */
void neighbor_joining(int num_genes, float *distance_matrix, int **joins) {
    // initiate values
    int L = num_genes;
    int indices[L*L]; // store the indices of the current nodes
    for (int i = 0; i < L; i++){
        indices[i] = i;
    }

    int new_node_index = num_genes; // index for the new node created by joining two nodes

    int num_joins = 0; // counter for the number of joins made

    //loop until all nodes have been joined into a single tree
    while (L > 2){

        //calculate net divergence
        float net_divergence[L];

        for (int i = 0; i < L; i++){
            net_divergence[i] = 0.0;
            for (int j = 0; j < L; j++){
                if (i != j){
                    net_divergence[i] += distance_matrix[i * num_genes + j];
                }
            }
            net_divergence[i] /= (L - 2);
        }

        //adjust distance matrix
        int min_i = 0;
        int min_j = 1;
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                if (i != j){
                    distance_matrix[i * num_genes + j] -= (net_divergence[i] + net_divergence[j]);

                    //find pair with smallest distance
                    if (distance_matrix[i * num_genes + j] < distance_matrix[min_i * num_genes + min_j]){
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
        float dist_iz = (distance_matrix[min_i * num_genes + min_j] + net_divergence[min_i] - net_divergence[min_j]) / 2.0;
        float dist_jz = (distance_matrix[min_i * num_genes + min_j] + net_divergence[min_j] - net_divergence[min_i]) / 2.0;

        // create new node and update distance matrix
        //replace min_i with new node
        indices[min_i] = new_node_index;
        new_node_index++;
        for (int k = 0; k < L; k++){
            if (k != min_i && k != min_j){
                distance_matrix[min_i * num_genes + k] = (distance_matrix[min_i * num_genes + k] + distance_matrix[min_j * num_genes + k] - distance_matrix[min_i * num_genes + min_j]) / 2.0;
                distance_matrix[k * num_genes + min_i] = distance_matrix[min_i * num_genes + k];
            }
        }

        //replace min_j with last element of matrix to reduce the size
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
    joins[num_joins][0] = indices[0];
    joins[num_joins][1] = indices[1];
}

/**
 * Creates a score matrix that allows for comparison of multiple sequences at once by averaging the match scores for each pair of sequences in the two groups being compared.
 * This function is similar to the score_matrix function, but it takes in two groups of sequences as opposed to just one sequence in each group.
 * additional variable num_seq_A and num_seq_B tell function how many different sequences are in each group
 * also score matrix is a float to allow for average values across each sequence comparison between groups
 */
void altered_score_matrix_for_msa(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int num_seq_A, int num_seq_B, float matrix[(A_size+1) * (B_size+1)]) {
    int num_rows = A_size + 1;
    int num_cols = B_size + 1;

    //Fill start of rows with 0, -1, -2, ect.
    for(int i = 0; i < (num_rows); i++){
        matrix[i * num_cols] = (float)(i * gap_penalty);
    }

    //Fill top of columns with 0, -1, -2, ect.
    for (int j = 0; j < (num_cols); j++){
        matrix[j] = (float)(j * gap_penalty);
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
    #pragma acc data copy(matrix[0:(num_rows*num_cols)]), copyin(gene_A[0:(A_size)], gene_B[0:(B_size)])
    {
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < total_waves; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < num_cols) ? 1 : wave - num_cols + 1;
            int i_end = (wave < num_rows) ? wave - 1 : num_rows - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #pragma acc parallel loop
            for (int i = i_start; i <= i_end; i++) {
                int j = wave - i;

                // this grabs the values from cells top, left and top-left diag.  cell
                //cell(i-1, j), cell(i, j-1), cell(i-1, j-1)
                float top_matrix  = matrix[(i - 1) * num_cols + j];
                float left_matrix = matrix[i * num_cols + (j - 1)];
                float diag_matrix = matrix[(i - 1) * num_cols + (j - 1)];

                float match = 0.0;
                
                //scores match for each individial pairing of sequences in the two groups and sums them
                for (int a = 0; a < num_seq_A; a++) {
                    for (int b = 0; b < num_seq_B; b++) {

                        match += (gene_A[a*A_size + i - 1] == gene_B[b*B_size + j - 1]) ? match_score : mismatch_penalty;
                    }
                }

                //gets average score for match score
                match /= (num_seq_A * num_seq_B);

                // proceeds with scoring as normal, but with the average match score for the two groups of sequences
                float top_score = top_matrix + gap_penalty;
                float left_score = left_matrix + gap_penalty;
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

        //Little double check on the array
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
float altered_traceback_for_msa(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, 
    int num_seq_A, int num_seq_B, float score_matrix[(A_size+1) * (B_size+1)], int* aligned_A, int* aligned_B, int *alignment_length){
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

        int match = (gene_A[i - 1] == gene_B[j - 1]) ? match_score : mismatch_penalty;
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
            if (match == match_score) {
                identity_count++;
            }
            residues_compared++;
            // printf("match: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene A
        else if (score == (score_left + gap_penalty)){
            // adding gap to all seq in A
            for (int a = 0; a < num_seq_A; a++) {
                aligned_A[align_index + a*(A_size+B_size)] = -1; //-1 represents a gap
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
                aligned_B[align_index + b*(A_size+B_size)] = -1; //-1 represents a gap
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
            aligned_B[align_index + b*(A_size+B_size)] = -1; //-1 represents a gap
        }
        align_index++;
        i--;
    }

    while (j > 0){
        for (int a = 0; a < num_seq_A; a++) {
            aligned_A[align_index + a*(A_size+B_size)] = -1; //-1 represents a gap
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