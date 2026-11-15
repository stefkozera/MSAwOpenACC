#ifndef MSA_TRACEBACK_H
#define MSA_TRACEBACK_H

#include "msa.h"

void traceback(NM_Components *NM_variable);
void flip_sequence(Gene_types *gene, int gene_size);

float needleman_wunsch_for_dists(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int* scores);
float traceback_for_identity(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int* score_matrix);
void score_matrix(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int *matrix);
#endif //MSA_TRACEBACK_H