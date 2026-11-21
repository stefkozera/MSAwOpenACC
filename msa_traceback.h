#ifndef MSA_TRACEBACK_H
#define MSA_TRACEBACK_H

#include "msa.h"

/**
 * Computes a single cell of the Needleman-Wunsch score matrix. The cell at
 * (i, j) only ever reads its three already-finalised neighbours:
 *   top (i-1, j), left (i, j-1) and diag (i-1, j-1)
 * so it is the smallest unit of work that is safe to run in parallel, provided
 * the caller has ordered the cells into anti-diagonal waves.
 *
 * Declared here, defined in msa_traceback.c. Both callers use it: the host path
 * (score_matrix) and the OpenACC wavefront in msa(). 'acc routine seq' has to
 * precede the declaration, otherwise a caller that includes this header cannot
 * compile the call as device code.
 *
 * Being a normal (external) function rather than a static inline one in this
 * header means it is not inlined across translation units -- nvc emits it as a
 * separate device routine -- which costs roughly 17% on the pairwise phase,
 * where it is called once per score-matrix cell.
 *
 * Scoring uses the GAP_PENALTY_VAL / MISMATCH_SCORE_VAL / MATCH_SCORE_VAL values
 * from msa.h; they are compile-time constants of the alignment, not per-call
 * tuning, so they are not passed in.
 */
#ifdef OPEN_ACC_ON
#pragma acc routine seq
#endif

void traceback(PSA_Components *PSA_variable);
void flip_sequence(Gene_types *gene, int gene_size);

float traceback_for_identity(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, int* score_matrix);
#endif //MSA_TRACEBACK_H