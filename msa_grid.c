/******************************************************************************
 * File:        msa_grid.c
 * Author:      Kailey MacFadyen, Stefanie Kozera
 * Date:        2026-11-22
 * Description: Functions related to grid (score matrix) fillout
 ******************************************************************************/

#include "msa_grid.h"

/**
 * Randomly generates a gene of size gene_size
 */
void InitRandomGene(Gene_types *gene_arr, int gene_size)
{
    for(int i = 0; i < gene_size; ++i){
        gene_arr[i] = (Gene_types)((rand() % (GENE_OPT - 1)) + 1);
    }
}

/**
 * Completes first row and column, fills with decending values based on the GAP_PENALTY_VAL
 */
void InitGapPenalty(int *matrix)
{
    #ifdef OPEN_ACC_ON
    #pragma acc loop vector
    #endif
    for(int i = 0; i < GRID_SIZE_Y; ++i)
    {
        //-1, -2, -3, etc, depends on GAP_PENALTY_VAL
        matrix[(size_t)i * GRID_SIZE_X] = i * GAP_PENALTY_VAL;
    }

    #ifdef OPEN_ACC_ON
    #pragma acc loop vector
    #endif
    for(int i = 0; i < GRID_SIZE_X; ++i)
    {
        //-1, -2, -3, etc, depends on GAP_PENALTY_VAL
        matrix[(size_t)i] = i * GAP_PENALTY_VAL;
    }
}

/**
 * Computes one cell of needleman-wunsch.  Reliant on top, left, right neighbor and X and Y genes
 */
int score_matrix_cell(const Gene_types *gene_A, const Gene_types *gene_B, int num_cols, int i, int j, int *matrix)
{
    //Collect values from top, left, and top left
    int top_matrix  = matrix[(size_t)(i - 1) * num_cols + j];
    int left_matrix = matrix[(size_t)i * num_cols + (j - 1)];
    int diag_matrix = matrix[(size_t)(i - 1) * num_cols + (j - 1)];

    //Comparison of genes, if they are the same use MATCH_SCORE_VAL, if not, MISMATCH_SCORE_VAL
    int match = (gene_A[i - 1] == gene_B[j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;

    //Create scores based on gap penalty and the results of the gene match
    int top_score  = top_matrix + GAP_PENALTY_VAL;
    int left_score = left_matrix + GAP_PENALTY_VAL;
    int diag_score = diag_matrix + match;

    //Of the three options we take the highest and use that in the cell we are calculating
    int best;
    if ((top_score >= left_score) && (top_score >= diag_score)){
        best = top_score;
    }
    else if((left_score >= top_score) && (left_score >= diag_score)){
        best = left_score;
    }
    else{
        best = diag_score;
    }

    matrix[(size_t)i * num_cols + j] = best;
    return best;
}

/**
 * loads a FASTA file and converts it to a sequence of nucleotide numbers
 * Returns the size of the sequence loaded into the gene array
 */
int load_FASTA(const char* filename, Gene_types* gene, int max_size) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[max_size];
    int gene_size = 0;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') {
            continue; // Skip header lines
        }
        for (int i = 0; line[i] != '\0' && line[i] != '\n'; i++) {
            // gets enum values for each nucleotide in fasta file
            switch (line[i]) {
                case 'A':
                    gene[gene_size] = A;
                    break;
                case 'C':
                    gene[gene_size] = C;
                    break;
                case 'G':
                    gene[gene_size] = G;
                    break;
                case 'T':
                    gene[gene_size] = T;
                    break;
                default:
                    fprintf(stderr, "Invalid nucleotide: %c\n", line[i]);
                    exit(EXIT_FAILURE);
            }
            gene_size++;
        }
    }

    fclose(file);

    return gene_size;
}

/**
 * OpenACC implementation - wavefront parallelization
 */
void grid_fillout_ACC(PSA_Components *PSA_variable)
{
    // Extract pointers to local variables to avoid host struct pointer dereferencing on GPU
    int *PSA_matrix = PSA_variable->PSA_matrix;
    const Gene_types *gene_X = PSA_variable->gene_X; // adjust type if not char
    const Gene_types *gene_Y = PSA_variable->gene_Y;

    //Copy the matrix into location easily accessable to GPU
    #ifdef OPEN_ACC_ON
    #pragma acc data copy(PSA_matrix[0:(size_t)GRID_SIZE_X*GRID_SIZE_Y]) \
                     copyin(gene_X[0:(size_t)GENE_SIZE_X], gene_Y[0:(size_t)GENE_SIZE_Y])
    #endif
    {
    
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < TOTAL_WAVES; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < GRID_SIZE_X) ? 1 : wave - GRID_SIZE_X + 1;
            int i_end = (wave < GRID_SIZE_Y) ? wave - 1 : GRID_SIZE_Y - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #ifdef OPEN_ACC_ON
            #pragma acc parallel loop gang vector present(PSA_matrix, gene_X, gene_Y)
            #endif
            for (int i = i_start; i <= i_end; i++) {
                int j = wave - i;

                score_matrix_cell(gene_Y, gene_X, GRID_SIZE_X, i, j, PSA_matrix);

            }
        }
    }
}

/**
 * Variable printout for PSA
 */
void var_printout_psa(PSA_Components *PSA_variable)
{

    printf("GENE_X:\n");
    for (int i = 0; i < GENE_SIZE_X; ++i){
            printf("%c", num_to_nucleotide(PSA_variable->gene_X[i]));
    }
    printf("\n");

    printf("GENE_Y:\n");
    for (int i = 0; i < GENE_SIZE_Y; ++i){
        printf("%c", num_to_nucleotide(PSA_variable->gene_Y[i]));
    }
    printf("\n");

    printf("Matrix Filled:\n");
    for(int i = 0; i < GRID_SIZE_Y; ++i){
        for(int j = 0; j < GRID_SIZE_X; ++j){
            printf("%d ", PSA_variable->PSA_matrix[(size_t)i * GRID_SIZE_X + j]);
        }
        printf("\n");
    }

    printf("GENE_X aligned:\n");
    for(int i = 0; i < (PSA_variable->alignment_length); ++i){
        printf("%c", num_to_nucleotide(PSA_variable->gene_X_aligned[i]));
    }
    printf("\n");

    printf("GENE_Y aligned:\n");
    for (int i = 0; i < (PSA_variable->alignment_length); ++i){
        printf("%c", num_to_nucleotide(PSA_variable->gene_Y_aligned[i]));
    }
    printf("\n");

    printf("Alignment Length:\n");
    printf("%d \n", PSA_variable->alignment_length);

    printf("Identity Score:\n");
    printf("%f \n", PSA_variable->identity_score);

    printf("Alignment Distance (small score = more closely aligned):\n");
    printf("%f \n", PSA_variable->alignment_dist);
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
                
                //scores match for each individial pairing of sequences in the two groups and sums them
                // done sequentially to not create inconsistencies when adding to same data space for match
                #ifdef OPEN_ACC_ON
                #pragma acc loop seq
                #endif
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
}