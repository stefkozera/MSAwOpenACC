#ifndef MSA_GRID_H
#define MSA_GRID_H

#include "msa.h"

void InitRandomGene(Gene_types *gene_arr, int gene_size);
void InitRandomGene_Int(int *gene_arr, int gene_size);
void InitGapPenalty(int *matrix);
int  load_FASTA(const char* filename, int* gene, int max_size);
void grid_fillout_ACC(PSA_Components *PSA_variable);
void var_printout_psa(PSA_Components *PSA_variable);

#endif //MSA_GRID_H