#include "msa_grid.h"

void InitRandomGene(Gene_types *gene_arr, int gene_size)
{
    for(int i = 0; i < gene_size; ++i){
        gene_arr[i] = (Gene_types)((rand() % (GENE_OPT - 1)) + 1);
    }
}

void InitGapPenalty(int *matrix)
{
    #ifdef OPEN_ACC_ON
    #pragma acc loop vector
    #endif
    for(int i = 0; i < GRID_SIZE_Y; ++i)
    {
        matrix[i * GRID_SIZE_X] = i * GAP_PENALTY_VAL;
    }

    #ifdef OPEN_ACC_ON
    #pragma acc loop vector
    #endif
    for(int i = 0; i < GRID_SIZE_X; ++i)
    {
        matrix[i] = i * GAP_PENALTY_VAL;
    }
}

/**
 * Computes one cell of needleman-wunsch.  Reliant on top, left, right neighbor and X and Y genes
 */
int score_matrix_cell(const Gene_types *gene_A, const Gene_types *gene_B, int num_cols, int i, int j, int *matrix)
{
    int top_matrix  = matrix[(i - 1) * num_cols + j];
    int left_matrix = matrix[i * num_cols + (j - 1)];
    int diag_matrix = matrix[(i - 1) * num_cols + (j - 1)];

    int match = (gene_A[i - 1] == gene_B[j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;

    int top_score  = top_matrix + GAP_PENALTY_VAL;
    int left_score = left_matrix + GAP_PENALTY_VAL;
    int diag_score = diag_matrix + match;

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

    matrix[i * num_cols + j] = best;
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
            printf("%d ", PSA_variable->PSA_matrix[i * GRID_SIZE_X + j]);
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
