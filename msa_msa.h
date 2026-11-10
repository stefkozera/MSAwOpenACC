#ifndef MSA_MSA_H
#define MSA_MSA_H

#include "msa.h"

int msa(int** genes, int* gene_sizes, int num_genes, int gap_penalty, int mismatch_penalty, int match_score, int* aligned_genes);
void neighbor_joining(int num_genes, float *distance_matrix, int joins[num_genes-1][2]);
void altered_score_matrix_for_msa(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int num_seq_A, int num_seq_B, float matrix[(A_size+1) * (B_size+1)]);
float altered_traceback_for_msa(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, 
    int num_seq_A, int num_seq_B, float score_matrix[(A_size+1) * (B_size+1)], int* aligned_A, int* aligned_B, int *alignment_length);
void print_aligned_genes(int *alignment, int num_genes, int alignment_length);

#endif //MSA_MSA_H
