#include "msa_traceback.h"
#include "msa_msa.h"
#define GENE_SIZE 600
#define GRID_SIZE (GENE_SIZE+1)

#define GAP_PENALTY_VAL -1
#define MATCH_SCORE_VAL 1
#define MISMATCH_SCORE_VAL -1

struct timeval start, end;
long t_us;

int gene_A[GENE_SIZE];
int gene_B[GENE_SIZE];
int gene_C[GENE_SIZE];

int main() {

    srand(time(NULL)); //seed random number generator

    //Randomly generate a gene out of the enum Gene_types
    // for (int i = 0; i < GENE_SIZE; ++i)
    // {
    //     gene_A[i] = (Gene_types)(rand() % GENE_COUNT);
    //     gene_B[i] = (Gene_types)(rand() % GENE_COUNT);

    // }

    // int a_size = load_FASTA("GenBank_Sequences/Hemoglobin_alpha1/homo_sapiens_HBA1.fasta", gene_A, 600);

    // int b_size = load_FASTA("GenBank_Sequences/Hemoglobin_alpha1/mus_musculus_Hba1.fasta", gene_B, 600);

    // int c_size = load_FASTA("GenBank_Sequences/Hemoglobin_alpha1/macaca_mulatta_HBA1.fasta", gene_C, 600);

    // printf("\n");
    
    // for(int i = 0; i < a_size; i++){
    //     printf("%c", num_to_nucleotide(gene_A[i]));
    // }
    // printf("\n");

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
        return 1;
    }

    genes[0] = gene_A;
    genes[1] = gene_B;
    genes[2] = gene_C;

    int gene_sizes[3];
    gene_sizes[0] = a_size;
    gene_sizes[1] = b_size;
    gene_sizes[2] = c_size;

    int aligned_A[a_size+b_size];
    int aligned_B[a_size+b_size];
    int alignment_length = 0;

    int aligned_genes[a_size+b_size+c_size];

    gettimeofday(&start,NULL);

//    float distance = needleman_wunsch(gene_B, gene_C, GAP_PENALTY_VAL, MISMATCH_SCORE_VAL, MATCH_SCORE_VAL, b_size, c_size, aligned_A, aligned_B, &alignment_length);

//    print_alignment(aligned_A, aligned_B, alignment_length);

//    printf("%f\n", distance);
    int align_length = msa(genes, gene_sizes, num_genes, GAP_PENALTY_VAL, MISMATCH_SCORE_VAL, MATCH_SCORE_VAL, aligned_genes);

    gettimeofday(&end,NULL);
    t_us = (end.tv_sec - start.tv_sec)*1000000 + end.tv_usec - start.tv_usec;
    print_aligned_genes(aligned_genes, num_genes, align_length);

    free(genes);

    printf("Success! Time was %ld\n", t_us);
}
