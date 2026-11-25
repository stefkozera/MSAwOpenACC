/******************************************************************************
 * File:        msa_traceback.h
 * Author:      Kailey MacFadyen, Stefanie Kozera
 * Date:        2026-11-22
 * Description: h file for msa_traceback.c
 ******************************************************************************/

#ifndef MSA_TRACEBACK_H
#define MSA_TRACEBACK_H

#include "msa.h"

#ifdef OPEN_ACC_ON
#pragma acc routine seq
#endif

void traceback(PSA_Components *PSA_variable);
void flip_sequence(Gene_types *gene, int gene_size);
float traceback_for_identity(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, int* score_matrix);
float altered_traceback_for_msa(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, 
int num_seq_A, int num_seq_B, float score_matrix[(A_size+1) * (B_size+1)], Gene_types* aligned_A, Gene_types* aligned_B, int *alignment_length);
#endif //MSA_TRACEBACK_H