#include "msa_fun.h"
#define GENE_SIZE 10
#define GRID_SIZE (GENE_SIZE+1)

#define GAP_PENALTY_VAL -1
#define MATCH_SCORE_VAL 1
#define MISMATCH_SCORE_VAL -1

struct timeval start, end;
long t_us;

int gene_A[GENE_SIZE];
int gene_B[GENE_SIZE];

int main() {

    int *matrix = (int *)malloc(GRID_SIZE * GRID_SIZE * sizeof(int));

    srand(time(NULL)); //seed random number generator

    //Randomly generate a gene out of the enum Gene_types
    for (int i = 0; i < GENE_SIZE; ++i)
    {
        gene_A[i] = (Gene_types)(rand() % GENE_COUNT);
        gene_B[i] = (Gene_types)(rand() % GENE_COUNT);

    }

    //Fill start of rows with 0, -1, -2, ect.
    for(int i = 0; i < GRID_SIZE; i++){
        matrix[i * GRID_SIZE] = i * GAP_PENALTY_VAL;
    }

    //Fill top of columns with 0, -1, -2, ect.
    for (int j = 0; j < GRID_SIZE; j++){
        matrix[j] = j * GAP_PENALTY_VAL;
    }

    //Little double check on the array
    for (int i = 0; i < GRID_SIZE; i++){
        for(int j = 0; j < GRID_SIZE; j++){
            printf("%d ", matrix[i*GRID_SIZE + j]);
        }
        printf("\n");
    }

    // Total diagonal waves = (GRID_SIZE + GRID_SIZE - 1)
    // Wave index starts at 2 because (0,0), (0,1), and (1,0) are already boundary-initialized
    int total_waves = GRID_SIZE + GRID_SIZE - 1;

    //Copy the matrix into location easily accessable to GPU
    #pragma acc data copy(matrix[0:GRID_SIZE*GRID_SIZE])
    {
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < total_waves; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < GRID_SIZE) ? 1 : wave - GRID_SIZE + 1;
            int i_end = (wave < GRID_SIZE) ? wave - 1 : GRID_SIZE - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #pragma acc parallel loop independent
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
    }

        //Little double check on the array
    for (int i = 0; i < GRID_SIZE; i++){
        for(int j = 0; j < GRID_SIZE; j++){
            printf("%d ", matrix[i*GRID_SIZE + j]);
        }
        printf("\n");
    }


   gettimeofday(&start,NULL);
   gettimeofday(&end,NULL);
   t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;

   printf("Success! Time was %ld\n", t_us);
}


