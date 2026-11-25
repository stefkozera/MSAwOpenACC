/******************************************************************************
 * File:        msa_traceback.c
 * Author:      Kailey MacFadyen, Stefanie Kozera
 * Date:        2026-11-22
 * Description: Function related to traceback functionality for both PSA and
 *              MSA
 ******************************************************************************/
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

/**
 * Creates a score matrix that allows for comparison of multiple sequences at once by averaging the match scores for each pair of sequences in the two groups being compared.
 * This function is similar to the score_matrix function, but it takes in two groups of sequences as opposed to just one sequence in each group.
 * additional variable num_seq_A and num_seq_B tell function how many different sequences are in each group
 * also score matrix is a float to allow for average values across each sequence comparison between groups
 */
void altered_score_matrix_for_msa(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, int num_seq_A, int num_seq_B, float matrix[(A_size+1) * (B_size+1)]) {
    int num_rows = A_size + 1;
    int num_cols = B_size + 1;

    //Fill start of rows with 0, -1, -2, ect.
    for(int i = 0; i < (num_rows); i++){
        matrix[i * num_cols] = (float)(i * GAP_PENALTY_VAL);
    }

    //Fill top of columns with 0, -1, -2, ect.
    for (int j = 0; j < (num_cols); j++){
        matrix[j] = (float)(j * GAP_PENALTY_VAL);
    }

    //Little double check on the array
    // for (int i = 0; i < (num_rows); i++){
    //     for(int j = 0; j < (num_cols); j++){
    //         printf("%f ", matrix[i * num_cols + j]);
    //     }
    //     printf("\n");
    // }

    // Total diagonal waves = (GRID_SIZE + GRID_SIZE - 1)
    // Wave index starts at 2 because (0,0), (0,1), and (1,0) are already boundary-initialized
    int total_waves = num_rows + num_cols - 1;

    //Copy the matrix into location easily accessable to GPU
#ifdef OPEN_ACC_ON
#pragma acc data copy(matrix[0:(size_t)num_rows*num_cols]) \
                 copyin(gene_A[0:(num_seq_A*A_size)], gene_B[0:(num_seq_B*B_size)])
#endif
    {
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < total_waves; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < num_cols) ? 1 : wave - num_cols + 1;
            int i_end = (wave < num_rows) ? wave - 1 : num_rows - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #ifdef OPEN_ACC_ON
            #pragma acc parallel loop gang vector present(matrix, gene_A, gene_B)
            #endif
            for (int i = i_start; i <= i_end; i++) {
                int j = wave - i;

                // this grabs the values from cells top, left and top-left diag.  cell
                //cell(i-1, j), cell(i, j-1), cell(i-1, j-1)
                float top_matrix  = matrix[(i - 1) * num_cols + j];
                float left_matrix = matrix[i * num_cols + (j - 1)];
                float diag_matrix = matrix[(i - 1) * num_cols + (j - 1)];

                float match = 0.0;
                
                #ifdef OPEN_ACC_ON
                #pragma acc loop seq
                #endif
                //scores match for each individial pairing of sequences in the two groups and sums them
                for (int a = 0; a < num_seq_A; a++) {
                    #ifdef OPEN_ACC_ON
                    #pragma acc loop seq
                    #endif
                    for (int b = 0; b < num_seq_B; b++) {
                        match += (gene_A[a*A_size + i - 1] == gene_B[b*B_size + j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;
                    }
                }

                //gets average score for match score
                match /= (num_seq_A * num_seq_B);

                // proceeds with scoring as normal, but with the average match score for the two groups of sequences
                float top_score = top_matrix + GAP_PENALTY_VAL;
                float left_score = left_matrix + GAP_PENALTY_VAL;
                float diag_score = diag_matrix + match;

                if ((top_score >= left_score) && (top_score >= diag_score)){
                    matrix[i * num_cols + j] = top_score;
                }
                else if((left_score >= top_score) && (left_score >= diag_score)){
                    matrix[i * num_cols + j] = left_score;
                }
                else{
                    matrix[i * num_cols + j] = diag_score;
                }
            }
        }
    }
}

/**
 * traceback function for msa alignment
 * similar to the traceback algorithm for pairwise, except it allows for two groups of sequences to be aligned, as opposed to just two individual sequences
 * additional variable num_seq_A and num_seq_B tell function how many different sequences are in each group
 */
float altered_traceback_for_msa(Gene_types* gene_A, Gene_types* gene_B, int A_size, int B_size, 
    int num_seq_A, int num_seq_B, float score_matrix[(A_size+1) * (B_size+1)], Gene_types* aligned_A, Gene_types* aligned_B, int *alignment_length){
    // indices for the score matrix
    int i = A_size;
    int j = B_size;

    int num_cols = B_size + 1;

    // counter for aligned sequences
    int align_index = 0;

    int identity_count = 0; // counter for identical matches
    int residues_compared = 0; // counter for total residues compared (excluding gaps)

    // main traceback loop
    while(i > 0 && j > 0){
        int score = score_matrix[i * num_cols + j];
        int score_diag = score_matrix[(i-1) * num_cols + (j-1)];
        int score_left = score_matrix[i * num_cols + (j-1)];
        int score_up = score_matrix[(i-1) * num_cols + j];

        float match = 0.0;
                
        //scores match for each individial pairing of sequences in the two groups and sums them
        for (int a = 0; a < num_seq_A; a++) {
            for (int b = 0; b < num_seq_B; b++) {
                match += (gene_A[a*A_size + i - 1] == gene_B[b*B_size + j - 1]) ? MATCH_SCORE_VAL : MISMATCH_SCORE_VAL;
            }
        }
        //gets average score for match score
        match /= (num_seq_A * num_seq_B);

        // printf("score: %d, score_diag: %d, score_left: %d, score_up: %d\n", score, score_diag, score_left, score_up);
        // score from match/mismatch
        if (score == (score_diag + match)){
            //aligning all sequences in A
            for (int a = 0; a < num_seq_A; a++) {
                aligned_A[align_index + a*(A_size+B_size)] = gene_A[a*A_size + i - 1];
            }
            //aligning all sequences in B
            for (int b = 0; b < num_seq_B; b++) {
                aligned_B[align_index + b*(A_size+B_size)] = gene_B[b*B_size + j - 1];
            }
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
            // adding gap to all seq in A
            for (int a = 0; a < num_seq_A; a++) {
                aligned_A[align_index + a*(A_size+B_size)] = GAP; // GAP represents a gap
            }
            // adding the corresponding nucleotide from gene B to all seq in B
            for (int b = 0; b < num_seq_B; b++) {
                aligned_B[align_index + b*(A_size+B_size)] = gene_B[b*B_size + j-1];
            }
            j--;
            // printf("gap in gene A: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene B
        else {
            // adding the corresponding nucleotide from gene A to all seq in A
            for (int a = 0; a < num_seq_A; a++) {
                aligned_A[align_index + a*(A_size+B_size)] = gene_A[a*A_size + i-1];
            }
            // adding gap to all seq in B
            for (int b = 0; b < num_seq_B; b++) {
                aligned_B[align_index + b*(A_size+B_size)] = GAP; // GAP represents a gap
            }
            i--;
            // printf("gap in gene B: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // update align index
        align_index++;
    }

    // fill in remaining gaps
    while(i > 0){
        for (int a = 0; a < num_seq_A; a++) {
            aligned_A[align_index + a*(A_size+B_size)] = gene_A[a*A_size + i-1];
        }
        for (int b = 0; b < num_seq_B; b++) {
            aligned_B[align_index + b*(A_size+B_size)] = GAP; // GAP represents a gap
        }
        align_index++;
        i--;
    }

    while (j > 0){
        for (int a = 0; a < num_seq_A; a++) {
            aligned_A[align_index + a*(A_size+B_size)] = GAP; // GAP represents a gap
        }
        for (int b = 0; b < num_seq_B; b++) {
            aligned_B[align_index + b*(A_size+B_size)] = gene_B[b*B_size + j-1];
        }
        align_index++;
        j--;
    }

    *alignment_length = align_index;

    float identity_score = (residues_compared > 0) ? (float)identity_count / (float)residues_compared : 0;

    return identity_score;
}