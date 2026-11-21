#include "msa_traceback.h"

/**
 * traceback function to construct the alignment of the sequences based on the scoring matrix.
 * max size of alignment is A_size + B_size (in case of all gaps)
 * returns the length of the aligned sequences.
 * gap penalty and mismatch penalty are negative numbers, match score is just the inverse of mismatch penalty
 */
void traceback(PSA_Components *PSA_variable){
    // indices for the score matrix
    int i = GENE_SIZE_Y;
    int j = GENE_SIZE_X;

    size_t num_cols = (size_t)GENE_SIZE_X + 1;

    // counter for aligned sequences
    int align_index = 0;

    int identity_count = 0; // counter for identical matches
    int residues_compared = 0; // counter for total residues compared (excluding gaps)

    // main traceback loop
    while(i > 0 && j > 0){
        int score = PSA_variable->PSA_matrix[(size_t)i * num_cols + j];
        int score_diag = PSA_variable->PSA_matrix[(size_t)(i-1) * num_cols + (j-1)];
        int score_left = PSA_variable->PSA_matrix[(size_t)i * num_cols + (j-1)];

        int match = (PSA_variable->gene_Y[i - 1] == PSA_variable->gene_X[j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;
        // printf("score: %d, score_diag: %d, score_left: %d, score_up: %d\n", score, score_diag, score_left, score_up);
        // score from match/mismatch
        if (score == (score_diag + match)){
            PSA_variable->gene_Y_aligned[align_index] = PSA_variable->gene_Y[i-1];
            PSA_variable->gene_X_aligned[align_index] = PSA_variable->gene_X[j-1];
            i--;
            j--;

            // count identical matches and total residues compared
            if (match == MATCH_SCORE_VAL) {
                identity_count++;
            }
            residues_compared++;
            // printf("match: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene A
        else if (score == (score_left + GAP_PENALTY_VAL)){
            PSA_variable->gene_Y_aligned[align_index] = GAP; //255 represents a gap
            PSA_variable->gene_X_aligned[align_index] = PSA_variable->gene_X[j-1];
            j--;
            // printf("gap in gene A: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene B
        else {
            PSA_variable->gene_Y_aligned[align_index] = PSA_variable->gene_Y[i-1];
            PSA_variable->gene_X_aligned[align_index] = GAP; //255 represents a gap
            i--;
            // printf("gap in gene B: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // update align index
        align_index++;
    }

    // fill in remaining gaps
    while(i > 0){
        PSA_variable->gene_Y_aligned[align_index] = PSA_variable->gene_Y[i-1];
        PSA_variable->gene_X_aligned[align_index] = GAP; //255 represents a gap
        align_index++;
        i--;
    }

    while (j > 0){
        PSA_variable->gene_Y_aligned[align_index] = GAP; //255 represents a gap
        PSA_variable->gene_X_aligned[align_index] = PSA_variable->gene_X[j-1];
        align_index++;
        j--;
    }

    //flip sequence alignments so they are not backward
    flip_sequence(PSA_variable->gene_Y_aligned, align_index);
    flip_sequence(PSA_variable->gene_X_aligned, align_index);

    PSA_variable->alignment_length = align_index;

    PSA_variable->identity_score = (residues_compared > 0) ? (float)identity_count / (float)residues_compared : 0;
    PSA_variable->alignment_dist = 1 - PSA_variable->identity_score;
}

/**
 * Reverses sequence order
 */
void flip_sequence(Gene_types *gene, int gene_size){
    for (int k = 0; k < gene_size/2; k++) {
        Gene_types temp = gene[k];
        gene[k] = gene[gene_size - 1 - k];
        gene[gene_size - 1 - k] = temp;
    }
}

/**
 * traceback for calculating the identity between two sequences
 * does not build aligned arrays
 */
#ifdef OPEN_ACC_ON
#pragma acc routine seq
#endif
float traceback_for_identity(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, int* score_matrix){
    // indices for the score matrix
    int i = A_size;
    int j = B_size;

    size_t num_cols = (size_t)B_size + 1;

    int identity_count = 0; // counter for identical matches
    int residues_compared = 0; // counter for total residues compared (excluding gaps)

    // main traceback loop
    while(i > 0 && j > 0){
        int score = score_matrix[(size_t)i * num_cols + j];
        int score_diag = score_matrix[(size_t)(i-1) * num_cols + (j-1)];
        int score_left = score_matrix[(size_t)i * num_cols + (j-1)];

        int match = (gene_A[i - 1] == gene_B[j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;
        // score from match/mismatch
        if (score == (score_diag + match)){
            i--;
            j--;

            // count identical matches and total residues compared
            if (match == MATCH_SCORE_VAL) {
                identity_count++;
            }
            residues_compared++;
        }
        // score from a gap in gene A
        else if (score == (score_left + GAP_PENALTY_VAL)){
            j--;
        }
        // score from a gap in gene B
        else {
            i--;
        }
    }

    // fill in remaining gaps
    while(i > 0){
        i--;
    }

    while (j > 0){
        j--;
    }

    float identity_score = (residues_compared > 0) ? (float)identity_count / (float)residues_compared : 0;

    return identity_score;
}