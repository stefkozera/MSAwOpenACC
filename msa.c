/******************************************************************************
 * File:        msa.c
 * Author:      Kailey MacFadyen, Stefanie Kozera
 * Date:        2026-11-22
 * Description: Main file of MSA program.  Main function contains calls to PSA
 *              and MSA run.  Collects timings.  Misc. functionality
 ******************************************************************************/

#include "msa.h"
#include "msa_grid.h"
#include "msa_traceback.h"
#include "msa_msa.h"

/**
 * Main, calls PSA_Run() and MSA_Run()
 * OpenACC initiatlization as well
 */
int main() 
{
    //Timer variables
    struct timeval start_psa, end_psa;
    long t_us;

    srand(time(NULL));

    gettimeofday(&start_psa,NULL);

    //Open acc initialization function call.  Separating out openACC init to get more clean readings of PSA and MSA
    #ifdef OPEN_ACC_ON
    acc_init(acc_device_default);
    #endif

    gettimeofday(&end_psa,NULL);

    t_us = (end_psa.tv_sec - start_psa.tv_sec)*1000000 + end_psa.tv_usec - start_psa.tv_usec;
    printf("OpenACC init time: %ld\n", t_us);

    //PSA run
    #ifdef RUN_PSA
        PSA_Run();
    #endif

    //MSA run
    #ifdef RUN_MSA_TEST_VALUES
        MSA_Run_Test_Values();
    #elif defined(RUN_MSA)
        MSA_Run();
    #endif
}

/**
 * Pairwise sequence alignment.  Calls functions to initalize PSA structure, fills out score matrix and performs traceback.
 * Collects timings
 */
void PSA_Run(void)
{
    //Timing variables
    struct timeval start_psa, end_psa;
    long t_us;

    //Malloc allocation of PSA_var, structure of PSA variables
    PSA_Components *PSA_var = malloc(sizeof(PSA_Components));
    if(PSA_var == NULL){
        printf("Structure allocation failed!\n");
        abort();
    }

    //Malloc allocations of interior structure components, 
    init_PSA_Components(PSA_var);

    gettimeofday(&start_psa,NULL);

    //Gap penalty fill out of score matrix
    InitGapPenalty(PSA_var->PSA_matrix);

    //Complete the score matrix - parallelization here
    grid_fillout_ACC(PSA_var);

    //Traceback across score matrix
    traceback(PSA_var);

    gettimeofday(&end_psa,NULL);

    //Variable printout if needed
    #ifdef PRINT_VAR_OUTPUT
        var_printout_psa(PSA_var);
    #endif

    t_us = (end_psa.tv_sec - start_psa.tv_sec)*1000000 + end_psa.tv_usec - start_psa.tv_usec;
    printf("PSA time to complete: %ld\n", t_us);

    //Free
    deinit_PSA_Components(PSA_var);
}

/**
 * Multiple sequence alignment.  Variable initialization, random generation and timings
 * Collects timings
 */
void MSA_Run(void)
{
    //timing variables
    struct timeval start_msa, end_msa;
    long t_us;

    //Allocation of genes - 2d array, row for each gene
    Gene_types **genes = malloc(MSA_GENE_GEN_NUMBER * sizeof(Gene_types *));
    if (genes == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    //allocation of space per gene
    for(int i = 0; i < MSA_GENE_GEN_NUMBER; ++i)
    {
        genes[i] = malloc((size_t)MSA_GENE_GEN_LENGTH * sizeof(Gene_types));
        if (genes[i] == NULL) {
            printf("Memory allocation failed!\n");
            abort();
        }

        //Randomize genes for MSA
        InitRandomGene(genes[i], MSA_GENE_GEN_LENGTH);
    }

    //Array of gene sizes, input into MSA.  These can be variable, was tested and works, but in this implementation they are all the same length
    int gene_sizes[MSA_GENE_GEN_NUMBER];

    for(int i = 0; i < MSA_GENE_GEN_NUMBER; ++i)
    {
        gene_sizes[i] = MSA_GENE_GEN_LENGTH;
    }

    //Output aligned_genes
    Gene_types *aligned_genes = malloc(sizeof(Gene_types) * (size_t)MSA_GENE_GEN_NUMBER * MSA_GENE_GEN_NUMBER * MSA_GENE_GEN_LENGTH);
    if (aligned_genes == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    //Variable printout
    #ifdef PRINT_VAR_OUTPUT
        for(int i = 0; i < MSA_GENE_GEN_NUMBER; ++i)
        {
            printf("GENE_%d:\n", i);
            for(int j = 0; j < MSA_GENE_GEN_LENGTH; ++j){
                printf("%c", num_to_nucleotide(genes[i][j]));
            }
            printf("\n");
        }
    #endif

    gettimeofday(&start_msa,NULL);

    //Multiple Sequence Alignment
    int align_length = msa(genes, gene_sizes, MSA_GENE_GEN_NUMBER, aligned_genes);

    gettimeofday(&end_msa,NULL);

    #ifdef PRINT_VAR_OUTPUT
        print_aligned_genes(aligned_genes, MSA_GENE_GEN_NUMBER, align_length);
    #endif

    t_us = (end_msa.tv_sec - start_msa.tv_sec)*1000000 + end_msa.tv_usec - start_msa.tv_usec;

    printf("MSA time to complete: %ld\n", t_us);

    //Cleanup
    for (int i = 0; i < MSA_GENE_GEN_NUMBER; ++i) {
        free(genes[i]);
    }
    free(genes);
    free(aligned_genes);
}

/**
 * Multiple sequence alignment with consistent values for testing. Showcases variable gene length functionality as well
 * Collects timings
 */
void MSA_Run_Test_Values(void)
{
    struct timeval start_msa, end_msa;
    long t_us;

    // arrays for each gene
    Gene_types gene_A[4];
    Gene_types gene_B[4];
    Gene_types gene_C[5];

    //populate gene A
    gene_A[0] = A;
    gene_A[1] = T;
    gene_A[2] = G;
    gene_A[3] = C;

    int a_size = 4;

    //populate gene B
    gene_B[0] = A;
    gene_B[1] = G;
    gene_B[2] = C;
    gene_B[3] = T;

    int b_size = 4;

    //populate gene C
    gene_C[0] = A;
    gene_C[1] = T;
    gene_C[2] = C;
    gene_C[3] = C;
    gene_C[4] = T;

    int c_size = 5;
    int num_genes = 3; //total number of genes

    //location to store all genes
    Gene_types **genes = malloc(num_genes * sizeof(Gene_types *));
    if (genes == NULL) {
        printf("Memory allocation failed!\n");
        abort();
    }

    //place genes in array
    genes[0] = gene_A;
    genes[1] = gene_B;
    genes[2] = gene_C;

    //all gene sizes in array
    int gene_sizes[3];
    gene_sizes[0] = a_size;
    gene_sizes[1] = b_size;
    gene_sizes[2] = c_size;

    Gene_types aligned_genes[a_size+b_size+c_size]; //space for all aligned genes

    gettimeofday(&start_msa,NULL);

    int align_length = msa(genes, gene_sizes, num_genes, aligned_genes); //performs msa

    gettimeofday(&end_msa,NULL);

    print_aligned_genes(aligned_genes, num_genes, align_length); //prints result

    t_us = (end_msa.tv_sec - start_msa.tv_sec)*1000000 + end_msa.tv_usec - start_msa.tv_usec;

    printf("MSA time to complete: %ld\n", t_us);
    
}

/**
 * PSA_Components structure element memory allocation
 * And random gene generation for Gene_X and Gene_Y
 */
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
    PSA_variable->PSA_matrix = malloc((size_t)GRID_SIZE_X * GRID_SIZE_Y * sizeof(int));
    if(PSA_variable->PSA_matrix == NULL){
        printf("Array allocation failed!\n");
        abort();
    }

    //Random gene creation
    InitRandomGene(PSA_variable->gene_X, GENE_SIZE_X);
    InitRandomGene(PSA_variable->gene_Y, GENE_SIZE_Y);
}

/**
 * PSA_Components deinitialization
 */
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



