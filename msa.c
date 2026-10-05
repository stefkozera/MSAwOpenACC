#include "msa_fun.h"
#define GENE_SIZE 1000

struct timeval start, end;
long t_us;

int gene_A[GENE_SIZE];
int gene_B[GENE_SIZE];

int main() {

    srand(time(NULL)); //seed random number generator

    for (int i = 0; i < GENE_SIZE; ++i)
    {
        gene_A[i] = (Gene_types)(rand() % GENE_COUNT);
        gene_B[i] = (Gene_types)(rand() % GENE_COUNT);

    }

    

    gettimeofday(&start,NULL);





   gettimeofday(&end,NULL);
   t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;

   printf("Success! Time was %ld\n", t_us);
}

//#pragma acc parallel loop copy(array[0:N])
//  for(int i = 0; i < N; i++) {
//     array[i] = 3.0;
//  }