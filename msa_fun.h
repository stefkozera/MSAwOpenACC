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

int traceback(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int score_matrix[(A_size+1) * (B_size+1)], int* aligned_A, int* aligned_B);

char num_to_nucleotide(Gene_types nucleotide);

void needleman_wunsch(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size);

int load_FASTA(const char* filename, int* gene, int max_size);