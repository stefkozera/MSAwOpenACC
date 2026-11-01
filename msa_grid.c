#include "msa_grid.h"

void InitRandomGene(Gene_types *gene_X_arr, Gene_types *gene_Y_arr)
{
    gene_X_arr[0] = C;
    gene_X_arr[1] = C;
    gene_X_arr[2] = G;
    gene_X_arr[3] = T;

    gene_Y_arr[0] = G;
    gene_Y_arr[1] = C;
    gene_Y_arr[2] = A;
    gene_Y_arr[3] = T;

    /*for(int i = 0; i < GENE_SIZE_X; ++i){
        gene_X_arr[i] = (Gene_types)(rand() % GENE_COUNT);
    }

    for(int i = 0; i < GENE_SIZE_Y; ++i){
        gene_Y_arr[i] = (Gene_types)(rand() % GENE_COUNT);
    }*/
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
