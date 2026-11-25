/******************************************************************************
 * File:        msa_grid.h
 * Author:      Kailey MacFadyen, Stefanie Kozera
 * Date:        2026-11-22
 * Description: h file for msa_grid.c
 ******************************************************************************/

#ifndef MSA_GRID_H
#define MSA_GRID_H

#include "msa.h"

void InitRandomGene(Gene_types *gene_arr, int gene_size);
void InitGapPenalty(int *matrix);

int score_matrix_cell(const Gene_types *gene_A, const Gene_types *gene_B, int num_cols, int i, int j, int *matrix);
int  load_FASTA(const char* filename, Gene_types* gene, int max_size);
void grid_fillout_ACC(PSA_Components *PSA_variable);
void var_printout_psa(PSA_Components *PSA_variable);
void altered_score_matrix_for_msa(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, int num_seq_A, int num_seq_B, 
    float matrix[(A_size+1) * (B_size+1)]);

#endif //MSA_GRID_H