#include "msa.h"
#include "msa_grid.h"
#include "msa_traceback.h"
#include "msa_msa.h"

/*
        {Rows, Width, X}
{Columns [x, x, x, x,
Length    x, x, x, x,
Y}        x, x, x, x]

**** 2D_ACCESS(matrix, rows, columns, GRID_SIZE_X) ****
*/
int main() 
{
    struct timeval start_psa, end_psa;
    long t_us;

    srand(time(NULL));

    gettimeofday(&start_psa,NULL);

    acc_init(acc_device_default);

    gettimeofday(&end_psa,NULL);

    t_us = (end_psa.tv_sec - start_psa.tv_sec)*1000000 + end_psa.tv_usec - start_psa.tv_usec;
    printf("OpenACC init time: %ld\n", t_us);

    #ifdef RUN_PSA
        PSA_Run();
    #endif

    #ifdef RUN_MSA_TEST_VALUES
        MSA_Run_Test_Values();
    #elif defined(RUN_MSA)
        MSA_Run();
    #endif
}

void PSA_Run(void)
{
    struct timeval start_psa, end_psa;
    long t_us;

    PSA_Components *PSA_var = malloc(sizeof(PSA_Components));
    if(PSA_var == NULL){
        printf("Structure allocation failed!\n");
        abort();
    }

    init_PSA_Components(PSA_var);

    gettimeofday(&start_psa,NULL);

    grid_fillout_ACC(PSA_var);

    traceback(PSA_var);

    gettimeofday(&end_psa,NULL);

    var_printout_psa(PSA_var);

    t_us = (end_psa.tv_sec - start_psa.tv_sec)*1000000 + end_psa.tv_usec - start_psa.tv_usec;
    printf("PSA time to complete: %ld\n", t_us);

    deinit_PSA_Components(PSA_var);
}

void MSA_Run(void)
{
    struct timeval start_msa, end_msa;
    long t_us;

    int **genes = malloc(MSA_GENE_GEN_NUMBER * sizeof(int *));
    if (genes == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    for(int i = 0; i < MSA_GENE_GEN_NUMBER; ++i)
    {
        genes[i] = malloc((size_t)MSA_GENE_GEN_LENGTH * sizeof(int));
        if (genes[i] == NULL) {
            printf("Memory allocation failed!\n");
            abort();
        }

        InitRandomGene_Int(genes[i], MSA_GENE_GEN_LENGTH);
    }

    int gene_sizes[MSA_GENE_GEN_NUMBER];

    for(int i = 0; i < MSA_GENE_GEN_NUMBER; ++i)
    {
        gene_sizes[i] = MSA_GENE_GEN_LENGTH;
    }

    int aligned_genes[MSA_GENE_GEN_NUMBER * MSA_GENE_GEN_LENGTH];

    for(int i = 0; i < MSA_GENE_GEN_NUMBER; ++i)
    {
        printf("GENE_%d:\n", i);
        for(int j = 0; j < MSA_GENE_GEN_LENGTH; ++j){
            printf("%c", num_to_nucleotide(genes[i][j]));
        }
        printf("\n");
    }

    gettimeofday(&start_msa,NULL);

    int align_length = msa(genes, gene_sizes, MSA_GENE_GEN_NUMBER, GAP_PENALTY_VAL, MISMATCH_SCORE_VAL, MATCH_SCORE_VAL, aligned_genes);

    gettimeofday(&end_msa,NULL);

    print_aligned_genes(aligned_genes, MSA_GENE_GEN_NUMBER, align_length);

    t_us = (end_msa.tv_sec - start_msa.tv_sec)*1000000 + end_msa.tv_usec - start_msa.tv_usec;

    printf("MSA time to complete: %ld\n", t_us);

    //Cleanup
    for (int i = 0; i < MSA_GENE_GEN_NUMBER; ++i) {
        free(genes[i]);
    }
    free(genes);
}

void MSA_Run_Test_Values(void)
{
    struct timeval start_msa, end_msa;
    long t_us;

    int gene_A[4];
    int gene_B[4];
    int gene_C[5];

    gene_A[0] = A;
    gene_A[1] = T;
    gene_A[2] = G;
    gene_A[3] = C;

    int a_size = 4;

    gene_B[0] = A;
    gene_B[1] = G;
    gene_B[2] = C;
    gene_B[3] = T;

    int b_size = 4;

    gene_C[0] = A;
    gene_C[1] = T;
    gene_C[2] = C;
    gene_C[3] = C;
    gene_C[4] = T;

    int c_size = 5;
    int num_genes = 3;

    int **genes = (int**)malloc(num_genes * sizeof(int *));
    if (genes == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    genes[0] = gene_A;
    genes[1] = gene_B;
    genes[2] = gene_C;

    int gene_sizes[3];
    gene_sizes[0] = a_size;
    gene_sizes[1] = b_size;
    gene_sizes[2] = c_size;

    int aligned_genes[a_size+b_size+c_size];

    gettimeofday(&start_msa,NULL);

    int align_length = msa(genes, gene_sizes, num_genes, GAP_PENALTY_VAL, MISMATCH_SCORE_VAL, MATCH_SCORE_VAL, aligned_genes);

    gettimeofday(&end_msa,NULL);

    print_aligned_genes(aligned_genes, num_genes, align_length);

    t_us = (end_msa.tv_sec - start_msa.tv_sec)*1000000 + end_msa.tv_usec - start_msa.tv_usec;

    printf("MSA time to complete: %ld\n", t_us);
    
}


void init_PSA_Components(PSA_Components *PSA_variable)
{
    //Dynamic allocation
    PSA_variable->gene_X = malloc(GENE_SIZE_X * sizeof(Gene_types));
    if(PSA_variable->gene_X == NULL){
        printf("Array allocation failed!\n");
        abort();
    }

    PSA_variable->gene_Y = malloc(GENE_SIZE_Y * sizeof(Gene_types));
    if(PSA_variable->gene_Y == NULL){
        printf("Array allocation failed!\n");
        abort();
    }

    PSA_variable->gene_X_aligned = malloc((GENE_SIZE_X + GENE_SIZE_Y) * sizeof(Gene_types));
    if(PSA_variable->gene_X_aligned == NULL){
        printf("Array allocation failed!\n");
        abort();
    }

    PSA_variable->gene_Y_aligned = malloc((GENE_SIZE_X + GENE_SIZE_Y) * sizeof(Gene_types));
    if(PSA_variable->gene_Y_aligned == NULL){
        printf("Array allocation failed!\n");
        abort();
    }

    //Needleman-Wunsch is a 2d array hiding in a single pointer.  We're doing this because the performance might be modestly improved (less cache miss probably)
    //See the macro in msa.h "2D_ACCESS"
    PSA_variable->PSA_matrix = malloc(GRID_SIZE_X * GRID_SIZE_Y * sizeof(int));
    if(PSA_variable->PSA_matrix == NULL){
        printf("Array allocation failed!\n");
        abort();
    }

    InitRandomGene(PSA_variable->gene_X, GENE_SIZE_X);
    InitRandomGene(PSA_variable->gene_Y, GENE_SIZE_Y);
    InitGapPenalty(PSA_variable->PSA_matrix);
}

void deinit_PSA_Components(PSA_Components *PSA_variable)
{
    free(PSA_variable->gene_X);
    free(PSA_variable->gene_Y);
    free(PSA_variable->PSA_matrix);
    free(PSA_variable->gene_X_aligned);
    free(PSA_variable->gene_Y_aligned);
    free(PSA_variable);
}

/**
 * converts a nucleotide number to its character representation (A, C, G, T)
 * returns '-' for gaps or invalid numbers
 */
char num_to_nucleotide(Gene_types nucleotide){
    switch(nucleotide){
        case G:
            return 'g';
        case A:
            return 'a';
        case T:
            return 't';
        case C:
            return 'c';
        case GAP:
            return '-';
        default:
            return '-';
    }
}



