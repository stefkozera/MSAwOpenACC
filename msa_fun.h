#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

typedef enum {
    G,
    A,
    T,
    C,
    GENE_COUNT
} Gene_types;

void score_matrix(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int score_matrix[(A_size+1) * (B_size+1)]);

float traceback(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int score_matrix[(A_size+1) * (B_size+1)], int* aligned_A, int* aligned_B, int *alignment_length);

char num_to_nucleotide(Gene_types nucleotide);

float needleman_wunsch(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int aligned_A[A_size + B_size], int aligned_B[A_size + B_size], int* alignment_length);

int load_FASTA(const char* filename, int* gene, int max_size);

int msa(int** genes, int* gene_sizes, int num_genes, int gap_penalty, int mismatch_penalty, int match_score, int* aligned_genes);

