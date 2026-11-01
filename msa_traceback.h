#ifndef MSA_TRACEBACK_H
#define MSA_TRACEBACK_H

#include "msa.h"

void score_matrix(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int A_size, int B_size, int score_matrix[A_size+1][B_size+1]);

int traceback(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int A_size, int B_size, int score_matrix[A_size+1][B_size+1], int* aligned_A, int* aligned_B);

char num_to_nucleotide(Gene_types nucleotide);

void needleman_wunsch(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int A_size, int B_size, int scores[A_size+1][B_size+1]);

#endif //MSA_TRACEBACK_H