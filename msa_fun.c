#include "msa_fun.h"

void score_matrix(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int A_size, int B_size, int score_matrix[A_size+1][B_size+1]) {
    // generate first row and column of score matrix (gap penalties all the way down)
    for(int i = 0; i <= A_size; i++){
        score_matrix[i][0] = i * gap_penalty;
    }
    for(int j = 0; j <=B_size; j++){
        score_matrix[0][j] = j * gap_penalty;
    }

    // fill in scoring matrix (wavefront)
    for(int i = 1; i <= (A_size + B_size -1); i++){

    }
}

/**
 * traceback function to construct the alignment of the sequences based on the scoring matrix.
 * max size of alignment is A_size + B_size (in case of all gaps)
 * returns the length of the aligned sequences.
 * gap penalty and mismatch penalty are negative numbers, match score is just the inverse of mismatch penalty
 */
int traceback(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int A_size, int B_size, int score_matrix[A_size+1][B_size+1], int* aligned_A, int* aligned_B){
    // indices for the score matrix
    int i = A_size;
    int j = B_size;

    // counter for aligned sequences
    int align_index = 0;

    // main traceback loop
    while(i > 0 && j > 0){
        int score = score_matrix[i][j];
        int score_diag = score_matrix[i-1][j-1];
        int score_left = score_matrix[i][j-1];
        int score_up = score_matrix[i-1][j];
        //printf("score: %d, score_diag: %d, score_left: %d, score_up: %d\n", score, score_diag, score_left, score_up);
        // score from match
        if (score == (score_diag - mismatch_penalty) && gene_A[i-1] == gene_B[j-1]){
            aligned_A[align_index] = gene_A[i-1];
            aligned_B[align_index] = gene_B[j-1];
            i--;
            j--;
            //printf("match: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from mismatch
        else if (score == (score_diag + mismatch_penalty)){
            aligned_A[align_index] = gene_A[i-1];
            aligned_B[align_index] = gene_B[j-1];
            i--;
            j--;
            //printf("mismatch: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene A
        else if (score == (score_left + gap_penalty)){
            aligned_A[align_index] = -1; //-1 represents a gap
            aligned_B[align_index] = gene_B[j-1];
            j--;
            //printf("gap in gene A: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // score from a gap in gene B
        else {
            aligned_A[align_index] = gene_A[i-1];
            aligned_B[align_index] = -1; //-1 represents a gap
            i--;
            //printf("gap in gene B: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
        }
        // update align index
        align_index++;
    }

    // fill in remaining gaps
    while(i > 0){
        aligned_A[align_index] = gene_A[i-1];
        aligned_B[align_index] = -1; //-1 represents a gap
        align_index++;
        i--;
    }

    while (j > 0){
        aligned_A[align_index] = -1; //-1 represents a gap
        aligned_B[align_index] = gene_B[j-1];
        align_index++;
        j--;
    }

    return align_index;
}

/**
 * aligns two sequences using entire needleman-wunsch algorithm (score matrix and traceback)
 * prints the aligned sequences
 * can optionally pass a precomputed score matrix
 */
void needleman_wunsch(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int A_size, int B_size, int scores[A_size+1][B_size+1]) {
    // create scoring matrix
    if (scores == NULL)
    {
        int scores[A_size+1][B_size+1];
        score_matrix(gene_A, gene_B, gap_penalty, mismatch_penalty, A_size, B_size, scores);
    }

    // traceback to get aligned sequences
    int aligned_A[A_size + B_size];
    int aligned_B[A_size + B_size];
    int align_len = traceback(gene_A, gene_B, gap_penalty, mismatch_penalty, A_size, B_size, scores, aligned_A, aligned_B);

    //print alignment (starts at align_len-1 and works backwards)
    // print aligned A
    for(int i = align_len-1; i>=0; i--){
        printf("%c", num_to_nucleotide(aligned_A[i]));
    }
    printf("\n");

    //print aligned B
    for(int i = align_len-1; i>=0; i--){
        printf("%c", num_to_nucleotide(aligned_B[i]));
    }
    printf("\n");
}

/**
 * converts a nucleotide number to its character representation (A, C, G, T)
 * returns '-' for gaps or invalid numbers
 */
char num_to_nucleotide(Gene_types nucleotide){
    switch(nucleotide){
        case G:
            return 'G';
        case A:
            return 'A';
        case T:
            return 'T';
        case C:
            return 'C';
        default:
            return '-';
    }
}