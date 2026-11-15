#ifndef MSA_H
#define MSA_H

//H file stratgy
//msa.h file to include #defines, functions and types to be accessable to all c files
//msa_x.h files to include only #defines, functions and variables that are needed by a specific msa_x.c file

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>

//toggle defines
#define OPEN_ACC_ON
#define PRINT_VAR_OUTPUT
#define RUN_PSA
#define RUN_MSA

//Size of gene, compared during needleman-wunsch 
#define GENE_SIZE_X 6
#define GENE_SIZE_Y 6

//Grid size is one more than gene size, (mismatch scores are the additional row/column)
#define GRID_SIZE_X (GENE_SIZE_X + 1)
#define GRID_SIZE_Y (GENE_SIZE_Y + 1)

//Waves for wavefront parallel processing, calculation done along the anti-diagonals.  This is the number of anti-diagonals
#define TOTAL_WAVES (GRID_SIZE_X + GRID_SIZE_Y - 1)

//Valuations, used during needleman-wunsch, can be adjusted for different outcomes
#define GAP_PENALTY_VAL -1
#define MATCH_SCORE_VAL 1
#define MISMATCH_SCORE_VAL -1

//MACROS
#define ACCESS_2D(array, row, col, width) ((array)[((row) * (width)) + (col)])

typedef enum {
    GAP,
    G,
    A,
    T,
    C,
    GENE_OPT,
} Gene_types;

typedef struct {
    int *PSA_matrix;
    Gene_types *gene_X;
    Gene_types *gene_Y;
    Gene_types *gene_X_aligned;
    Gene_types *gene_Y_aligned;
    int alignment_length;
    float identity_score;
    float alignment_dist;
} PSA_Components;

void PSA_Run(void);
void MSA_Run(void);
char num_to_nucleotide(Gene_types nucleotide);
void init_PSA_Components(PSA_Components *PSA_variable);
void deinit_PSA_Components(PSA_Components *PSA_variable);

#endif //MSA_H