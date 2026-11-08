#ifndef MSA_TRACEBACK_H
#define MSA_TRACEBACK_H

#include "msa.h"

void score_matrix(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int matrix[(A_size+1) * (B_size+1)]);
float traceback(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int score_matrix[(A_size+1) * (B_size+1)], int* aligned_A, int* aligned_B, int *alignment_length);
float needleman_wunsch(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int aligned_A[A_size + B_size], int aligned_B[A_size + B_size], int* alignment_length);
char num_to_nucleotide(Gene_types nucleotide);

#endif //MSA_TRACEBACK_H