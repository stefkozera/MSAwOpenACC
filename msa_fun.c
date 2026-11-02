#include "msa_fun.h"

void score_matrix(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int matrix[A_size+1][B_size+1]) {
    int num_rows = A_size + 1;
    int num_cols = B_size + 1;

    //Fill start of rows with 0, -1, -2, ect.
    for(int i = 0; i < (num_rows); i++){
        matrix[i][0] = i * gap_penalty;
    }

    //Fill top of columns with 0, -1, -2, ect.
    for (int j = 0; j < (num_cols); j++){
        matrix[0][j] = j * gap_penalty;
    }

    //Little double check on the array
    for (int i = 0; i < (num_rows); i++){
        for(int j = 0; j < (num_cols); j++){
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }

    // Total diagonal waves = (GRID_SIZE + GRID_SIZE - 1)
    // Wave index starts at 2 because (0,0), (0,1), and (1,0) are already boundary-initialized
    int total_waves = num_rows + num_cols - 1;

    //Copy the matrix into location easily accessable to GPU
    #pragma acc data copy(matrix[0:(num_rows)][0:(num_cols)])
    {
        // Outer loop stays on the host, not parallel due to dependencies
        for (int wave = 2; wave < total_waves; wave++) {
            
            // Determine bounds for the current diagonal wave
            int i_start = (wave < num_cols) ? 1 : wave - num_cols + 1;
            int i_end = (wave < num_rows) ? wave - 1 : num_rows - 1;

            // Inner loop runs along the diagonal waves, all values calculated independently, so can be done by GPU
            #pragma acc parallel loop independent
            for (int i = i_start; i <= i_end; i++) {
                int j = wave - i;

                // this grabs the values from cells top, left and top-left diag.  cell
                //cell(i-1, j), cell(i, j-1), cell(i-1, j-1)
                int top_matrix  = matrix[i - 1][j];
                int left_matrix = matrix[i][j - 1];
                int diag_matrix = matrix[i - 1][j - 1];

                int match = (gene_A[i - 1] == gene_B[j - 1]) ? match_score : mismatch_penalty;

                int top_score = top_matrix + gap_penalty;
                int left_score = left_matrix + gap_penalty;
                int diag_score = diag_matrix + match;

                if ((top_score >= left_score) && (top_score >= diag_score)){
                    matrix[i][j] = top_score;
                }
                else if((left_score >= top_score) && (left_score >= diag_score)){
                    matrix[i][j] = left_score;
                }
                else{
                    matrix[i][j] = diag_score;
                }
            }
        }
    }

        //Little double check on the array
    for (int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_cols; j++){
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

/**
 * traceback function to construct the alignment of the sequences based on the scoring matrix.
 * max size of alignment is A_size + B_size (in case of all gaps)
 * returns the length of the aligned sequences.
 * gap penalty and mismatch penalty are negative numbers, match score is just the inverse of mismatch penalty
 */
int traceback(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int score_matrix[A_size+1][B_size+1], int* aligned_A, int* aligned_B){
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

        int match = (gene_A[i - 1] == gene_B[j - 1]) ? match_score : mismatch_penalty;
        //printf("score: %d, score_diag: %d, score_left: %d, score_up: %d\n", score, score_diag, score_left, score_up);
        // score from match/mismatch
        if (score == (score_diag + match)){
            aligned_A[align_index] = gene_A[i-1];
            aligned_B[align_index] = gene_B[j-1];
            i--;
            j--;
            //printf("match: aligned_A[%d] = %d, aligned_B[%d] = %d\n", align_index, aligned_A[align_index], align_index, aligned_B[align_index]);
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
void needleman_wunsch(int* gene_A, int* gene_B, int gap_penalty, int mismatch_penalty, int match_score, int A_size, int B_size, int scores[A_size+1][B_size+1]) {
    // create scoring matrix
    if (scores == NULL)
    {
        int scores[A_size+1][B_size+1];
        score_matrix(gene_A, gene_B, gap_penalty, mismatch_penalty, match_score, A_size, B_size, scores);
    }

    // traceback to get aligned sequences
    int aligned_A[A_size + B_size];
    int aligned_B[A_size + B_size];
    int align_len = traceback(gene_A, gene_B, gap_penalty, mismatch_penalty, match_score, A_size, B_size, scores, aligned_A, aligned_B);

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

/**
 * loads a FASTA file and converts it to a sequence of nucleotide numbers
 * Returns the size of the sequence loaded into the gene array
 */
int load_FASTA(const char* filename, int* gene, int max_size) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[max_size];
    int gene_size = 0;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') {
            continue; // Skip header lines
        }
        for (int i = 0; line[i] != '\0' && line[i] != '\n'; i++) {
            // gets enum values for each nucleotide in fasta file
            switch (line[i]) {
                case 'A':
                    gene[gene_size] = A;
                    break;
                case 'C':
                    gene[gene_size] = C;
                    break;
                case 'G':
                    gene[gene_size] = G;
                    break;
                case 'T':
                    gene[gene_size] = T;
                    break;
                default:
                    fprintf(stderr, "Invalid nucleotide: %c\n", line[i]);
                    exit(EXIT_FAILURE);
            }
            gene_size++;
        }
    }

    fclose(file);

    return gene_size;
}