#ifndef MSA_GRID_H
#define MSA_GRID_H

#include "msa.h"

void InitRandomGene(Gene_types *gene_arr, int gene_size);
void InitGapPenalty(int *matrix);

int score_matrix_cell(const Gene_types *gene_A, const Gene_types *gene_B, int num_cols, int i, int j, int *matrix);
int  load_FASTA(const char* filename, Gene_types* gene, int max_size);
void grid_fillout_ACC(PSA_Components *PSA_variable);
void var_printout_psa(PSA_Components *PSA_variable);

#endif //MSA_GRID_H