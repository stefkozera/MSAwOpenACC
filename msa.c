#include "msa.h"
#include "msa_grid.h"
#include "msa_traceback.h"

/*
        {Rows, Width, X}
{Columns [x, x, x, x,
Length    x, x, x, x,
Y}        x, x, x, x]

**** 2D_ACCESS(matrix, rows, columns, GRID_SIZE_X) ****
*/

#pragma acc routine gang
int main() {

    //Variables for time calculations
    struct timeval start, end;
    long t_us;

    //Dynamic matrix allocation
    Gene_types *gene_X = (Gene_types *)malloc(GENE_SIZE_X * sizeof(Gene_types));
    Gene_types *gene_Y = (Gene_types *)malloc(GENE_SIZE_Y * sizeof(Gene_types));

    //Needleman-Wunsch is a 2d array hiding in a single pointer.  We're doing this because the performance might be modestly improved (less cache miss probably)
    //See the macro in msa.h "2D_ACCESS"
    int *NM_matrix = (int *)malloc(GRID_SIZE_X * GRID_SIZE_Y * sizeof(int));

    //Seed random number generator
    srand(time(NULL)); //seed random number generator

    InitRandomGene(gene_X, gene_Y);
    InitGapPenalty(NM_matrix);

    gettimeofday(&start,NULL);

    //Copy the matrix into location easily accessable to GPU
    #pragma acc data copy(NM_matrix[0:(size_t)GRID_SIZE_X*GRID_SIZE_Y]) \
                     copyin(gene_X[0:(size_t)GENE_SIZE_X], gene_Y[0:(size_t)GENE_SIZE_Y])
    /*{
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < TOTAL_WAVES; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < GRID_SIZE) ? 1 : wave - GRID_SIZE + 1;
            int i_end = (wave < GRID_SIZE) ? wave - 1 : GRID_SIZE - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #pragma acc parallel loop gang independent
            for (int i = i_start; i <= i_end; i++) {
                int j = wave - i;

                // this grabs the values from cells top, left and top-left diag.  cell
                //cell(i-1, j), cell(i, j-1), cell(i-1, j-1)
                int top_matrix  = matrix[(i - 1) * GRID_SIZE + j];
                int left_matrix = matrix[i * GRID_SIZE + (j - 1)];
                int diag_matrix = matrix[(i - 1) * GRID_SIZE + (j - 1)];

                int match = (gene_A[i - 1] == gene_B[j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;

                int top_score = top_matrix + GAP_PENALTY_VAL;
                int left_score = left_matrix + GAP_PENALTY_VAL;
                int diag_score = diag_matrix + match;

                if ((top_score >= left_score) && (top_score >= diag_score)){
                    matrix[i * GRID_SIZE + j] = top_score;
                }
                else if((left_score >= top_score) && (left_score >= diag_score)){
                    matrix[i * GRID_SIZE + j] = left_score;
                }
                else{
                    matrix[i * GRID_SIZE + j] = diag_score;
                }
            }
        }
    }*/

        //Little double check on the array
    /*for (int i = 0; i < GRID_SIZE; i++){
        for(int j = 0; j < GRID_SIZE; j++){
            printf("%d ", matrix[i*GRID_SIZE + j]);
        }
        printf("\n");
    }*/

   gettimeofday(&end,NULL);
   t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;

   printf("Success! Time was %ld\n", t_us);
}


