#ifndef MSA_GRID_H
#define MSA_GRID_H

#include "msa.h"

void InitRandomGene(Gene_types *gene_X_arr, Gene_types *gen_Y_arry);
void InitGapPenalty(int *matrix);
int  load_FASTA(const char* filename, int* gene, int max_size);
void grid_fillout_ACC(NM_Components *NM_variable);
void var_printout(NM_Components *NM_variable);

#endif //MSA_GRID_H