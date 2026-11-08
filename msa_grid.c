#include "msa_grid.h"

void InitRandomGene(Gene_types *gene_X_arr, Gene_types *gene_Y_arr)
{
    for(int i = 0; i < GENE_SIZE_X; ++i){
        gene_X_arr[i] = (Gene_types)(rand() % GENE_COUNT);
    }

    for(int i = 0; i < GENE_SIZE_Y; ++i){
        gene_Y_arr[i] = (Gene_types)(rand() % GENE_COUNT);
    }
}

void InitGapPenalty(int *matrix)
{
    for(int i = 0; i < GRID_SIZE_Y; ++i)
    {
        ACCESS_2D(matrix, i, 0, GRID_SIZE_X) = i * GAP_PENALTY_VAL;
    }

    for(int i = 0; i < GRID_SIZE_X; ++i)
    {
        ACCESS_2D(matrix, 0, i, GRID_SIZE_X) = i * GAP_PENALTY_VAL;
    }
}

/**
 * loads a FASTA file and converts it to a sequence of nucleotide numbers
 * Returns the size of the sequence loaded into the gene array
 */
int load_FASTA(const char* filename, int* gene, int max_size) {
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

void grid_fillout_ACC(NM_Components *NM_variable)
{

    // Extract pointers to local variables to avoid host struct pointer dereferencing on GPU
    int * restrict NM_matrix = NM_variable->NM_matrix;
    const Gene_types * restrict gene_X = NM_variable->gene_X; // adjust type if not char
    const Gene_types * restrict gene_Y = NM_variable->gene_Y;

    //Copy the matrix into location easily accessable to GPU
    #pragma acc data copy(NM_matrix[0:(size_t)GRID_SIZE_X*GRID_SIZE_Y]) \
                     copyin(gene_X[0:(size_t)GENE_SIZE_X], gene_Y[0:(size_t)GENE_SIZE_Y])
    {
    
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < TOTAL_WAVES; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < GRID_SIZE_X) ? 1 : wave - GRID_SIZE_X + 1;
            int i_end = (wave < GRID_SIZE_Y) ? wave - 1 : GRID_SIZE_Y - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #pragma acc parallel loop gang vector present(NM_matrix, gene_X, gene_Y) async(1)
            for (int i = i_start; i <= i_end; i++) {
                int j = wave - i;

                // this grabs the values from cells top, left and top-left diag.  cell
                //cell(i-1, j), cell(i, j-1), cell(i-1, j-1)
                int top_matrix  = NM_matrix[(i - 1) * GRID_SIZE_X + j];
                int left_matrix = NM_matrix[i * GRID_SIZE_X + (j - 1)];
                int diag_matrix = NM_matrix[(i - 1) * GRID_SIZE_X + (j - 1)];

                int match = (gene_Y[i - 1] == gene_X[j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;

                int top_score = top_matrix + GAP_PENALTY_VAL;
                int left_score = left_matrix + GAP_PENALTY_VAL;
                int diag_score = diag_matrix + match;

                if ((top_score >= left_score) && (top_score >= diag_score)){
                    NM_matrix[i * GRID_SIZE_X + j] = top_score;
                }
                else if((left_score >= top_score) && (left_score >= diag_score)){
                    NM_matrix[i * GRID_SIZE_X + j] = left_score;
                }
                else{
                    NM_matrix[i * GRID_SIZE_X + j] = diag_score;
                }
            }
            #pragma acc wait
        }
    }
}

void var_printout(NM_Components *NM_variable)
{
    printf("GENE_X:\n");
    for (int i = 0; i < GENE_SIZE_X; ++i){
            printf("%d ", NM_variable->gene_X[i]);
    }
    printf("\n");

    printf("GENE_Y:\n");
    for (int i = 0; i < GENE_SIZE_Y; ++i){
        printf("%d ", NM_variable->gene_Y[i]);
    }
    printf("\n");

    printf("Matrix Filled:\n");
    for (int i = 0; i < GRID_SIZE_Y; ++i){
        for(int j = 0; j < GRID_SIZE_X; ++j){
            printf("%d ", ACCESS_2D(NM_variable->NM_matrix, i, j, GRID_SIZE_X));
        }
        printf("\n");
    }
}
