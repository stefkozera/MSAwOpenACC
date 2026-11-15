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
    PSA_Run();
    MSA_Run();
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

    

    gettimeofday(&start_msa,NULL);



    gettimeofday(&end_msa,NULL);

    t_us = (end_msa.tv_sec - start_msa.tv_sec)*1000000 + end_msa.tv_usec - start_msa.tv_usec;

    printf("MSA time to complete: %ld\n", t_us);
}


void init_PSA_Components(PSA_Components *PSA_variable)
{
    //Dynamic allocation
    PSA_variable->gene_X = (Gene_types *)malloc(GENE_SIZE_X * sizeof(Gene_types));
    if(PSA_variable->gene_X == NULL){
        printf("Array allocation failed!\n");
        free(PSA_variable);
        abort();
    }

    PSA_variable->gene_Y = (Gene_types *)malloc(GENE_SIZE_Y * sizeof(Gene_types));
    if(PSA_variable->gene_Y == NULL){
        printf("Array allocation failed!\n");
        free(PSA_variable);
        abort();
    }

    PSA_variable->gene_X_aligned = (Gene_types *)malloc((GENE_SIZE_X + GENE_SIZE_Y) * sizeof(Gene_types));
    if(PSA_variable->gene_X_aligned == NULL){
        printf("Array allocation failed!\n");
        free(PSA_variable);
        abort();
    }

    PSA_variable->gene_Y_aligned = (Gene_types *)malloc((GENE_SIZE_X + GENE_SIZE_Y) * sizeof(Gene_types));
    if(PSA_variable->gene_Y_aligned == NULL){
        printf("Array allocation failed!\n");
        free(PSA_variable);
        abort();
    }

    //Needleman-Wunsch is a 2d array hiding in a single pointer.  We're doing this because the performance might be modestly improved (less cache miss probably)
    //See the macro in msa.h "2D_ACCESS"
    PSA_variable->PSA_matrix = (int *)malloc(GRID_SIZE_X * GRID_SIZE_Y * sizeof(int));
    if(PSA_variable->PSA_matrix == NULL){
        printf("Array allocation failed!\n");
        free(PSA_variable);
        abort();
    }

    InitRandomGene(PSA_variable->gene_X, PSA_variable->gene_Y);
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



