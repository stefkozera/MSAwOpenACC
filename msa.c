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
void init_NM_Components(NM_Components *NM_variable)
{
    //Dynamic allocation
    NM_variable->gene_X = (Gene_types *)malloc(GENE_SIZE_X * sizeof(Gene_types));
    if(NM_variable->gene_X == NULL){
        printf("Array allocation failed!\n");
        free(NM_variable);
        abort();
    }

    NM_variable->gene_Y = (Gene_types *)malloc(GENE_SIZE_Y * sizeof(Gene_types));
    if(NM_variable->gene_Y == NULL){
        printf("Array allocation failed!\n");
        free(NM_variable);
        abort();
    }

    //Needleman-Wunsch is a 2d array hiding in a single pointer.  We're doing this because the performance might be modestly improved (less cache miss probably)
    //See the macro in msa.h "2D_ACCESS"
    NM_variable->NM_matrix = (int *)malloc(GRID_SIZE_X * GRID_SIZE_Y * sizeof(int));
    if(NM_variable->NM_matrix == NULL){
        printf("Array allocation failed!\n");
        free(NM_variable);
        abort();
    }

    InitRandomGene(NM_variable->gene_X, NM_variable->gene_Y);
    InitGapPenalty(NM_variable->NM_matrix);
}

void deinit_NM_Components(NM_Components *NM_variable)
{
    free(NM_variable->gene_X);
    free(NM_variable->gene_Y);
    free(NM_variable->NM_matrix);
    free(NM_variable);
}

//#pragma acc routine gang
int main() {

    //Variables for time calculations
    struct timeval start, end;
    long t_us;

    NM_Components *NM_var = malloc(sizeof(NM_Components));
    if(NM_var == NULL){
        printf("Structure allocation failed!\n");
        abort();
    }

    init_NM_Components(NM_var);

    //Seed random number generator
    srand(time(NULL)); //seed random number generator
    gettimeofday(&start,NULL);

    grid_fillout_ACC(NM_var);

    gettimeofday(&end,NULL);

    var_printout(NM_var);

    t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;

    printf("Success! Time was %ld\n", t_us);

    deinit_NM_Components(NM_var);
}




