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
    
    #ifdef ACC_ON
        #pragma acc data copy(NM_matrix[0:GRID_SIZE_X*GRID_SIZE_Y])
        {
    #endif
        
            for (int wave = 2; wave < TOTAL_WAVES; ++wave)
            {
                int i_start = (wave < GRID_SIZE_X) ? 1 : wave - GRID_SIZE_X + 1;
                int i_end = (wave < GRID_SIZE_Y) ? wave - 1 : GRID_SIZE_Y - 1;

                #ifdef ACC_ON
                    #pragma acc parallel loop independent
                #endif

                for (int i = i_start; i <= i_end; i++) 
                {
                    int j = wave - 1;
                    
                    // this grabs the values from cells top, left and top-left diag.  cell
                    //cell(i-1, j), cell(i, j-1), cell(i-1, j-1)
                    int top_matrix = ACCESS_2D(NM_matrix, i-1, j, GRID_SIZE_X);
                    int left_matrix = ACCESS_2D(NM_matrix, i, j-1, GRID_SIZE_X);
                    int diag_matrix = ACCESS_2D(NM_matrix, i-1, j-1, GRID_SIZE_X);

                    int match = (gene_Y[i-1] == gene_X[j-1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;

                    int top_score = top_matrix + GAP_PENALTY_VAL;
                    int left_score = left_matrix + GAP_PENALTY_VAL;
                    int diag_score = diag_matrix + match;


                    if ((top_score >= left_score) && (top_score >= diag_score)){
                        ACCESS_2D(NM_matrix, i, j, GRID_SIZE_X) = top_score;
                    }
                    else if((left_score >= top_score) && (left_score >= diag_score)){
                        ACCESS_2D(NM_matrix, i, j, GRID_SIZE_X) = left_score;
                    }
                    else{
                        ACCESS_2D(NM_matrix, i, j, GRID_SIZE_X) = diag_score;
                    }
                    
                }
            }

    #ifdef ACC_ON
        }
    #endif

    

    #ifdef PRINT_CHECKS //Grid with the gap penalties included
        printf("GENE_X:\n");
        for (int i = 0; i < GENE_SIZE_X; ++i){
            printf("%d ", gene_X[i]);
        }
        printf("\n");

        printf("GENE_Y:\n");
        for (int i = 0; i < GENE_SIZE_Y; ++i){
            printf("%d ", gene_Y[i]);
        }
        printf("\n");

        printf("Matrix Filled:\n");
        for (int i = 0; i < GRID_SIZE_Y; ++i){
            for(int j = 0; j < GRID_SIZE_X; ++j){
                printf("%d ", ACCESS_2D(NM_matrix, i, j, GRID_SIZE_X));
            }
            printf("\n");
        }
    #endif


   /*

   gettimeofday(&start,NULL);
   gettimeofday(&end,NULL);
   t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;

   printf("Success! Time was %ld\n", t_us);
   */

   free(NM_matrix);
   free(gene_X);
   free(gene_Y);
   
}


