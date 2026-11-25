/******************************************************************************
 * File:        msa_msa.h
 * Author:      Kailey MacFadyen, Stefanie Kozera
 * Date:        2026-11-22
 * Description: h file for msa_msa.c
 ******************************************************************************/

#ifndef MSA_MSA_H
#define MSA_MSA_H

#include "msa.h"
#include "msa_grid.h"

int msa(Gene_types** genes, int* gene_sizes, int num_genes, Gene_types* aligned_genes);
void neighbor_joining(int num_genes, float *distance_matrix, int joins[num_genes-1][2]);
void altered_score_matrix_for_msa(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, int num_seq_A, int num_seq_B, float matrix[(A_size+1) * (B_size+1)]);
void print_aligned_genes(Gene_types *alignment, int num_genes, int alignment_length);

#endif //MSA_MSA_H
