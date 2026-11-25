/******************************************************************************
 * File:        msa.h
 * Author:      Kailey MacFadyen, Stefanie Kozera
 * Date:        2026-11-22
 * Description: Main h file of MSA program.  Contains modifiable defines to change
 *              program functionality
 ******************************************************************************/

#ifndef MSA_H
#define MSA_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <openacc.h>

//toggle defines
//When defined, run with OpenACC in parallel on GPU, when not define run sequentially on CPU
#define OPEN_ACC_ON

//When defined, print output variables in terminal.  Warning, don't use when gene lengths are long
#define PRINT_VAR_OUTPUT

//When defined, run PSA
#define RUN_PSA

//When defined, run MSA
#define RUN_MSA

//When defined, run MSA with test values
#undef RUN_MSA_TEST_VALUES

//Program Wide Defines
//Valuations, used during needleman-wunsch, can be adjusted for different outcomes
#define GAP_PENALTY_VAL -1
#define MATCH_SCORE_VAL 1
#define MISMATCH_SCORE_VAL -1

//MSA Defines
//Length and amount of generated Genes for MSA
#define MSA_GENE_GEN_NUMBER 10
#define MSA_GENE_GEN_LENGTH 10

//PSA Defines
//Size of gene, compared during needleman-wunsch 
#define GENE_SIZE_X 10
#define GENE_SIZE_Y 10

//Grid size is one more than gene size, (mismatch scores are the additional row/column)
#define GRID_SIZE_X (GENE_SIZE_X + 1)
#define GRID_SIZE_Y (GENE_SIZE_Y + 1)

//Waves for wavefront parallel processing, calculation done along the anti-diagonals.  This is the number of anti-diagonals
#define TOTAL_WAVES (GRID_SIZE_X + GRID_SIZE_Y - 1)

//Gene enumeration
typedef enum {
    GAP,
    G,
    A,
    T,
    C,
    GENE_OPT,
} Gene_types;

//PSA structure
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
void MSA_Run_Test_Values(void);
char num_to_nucleotide(Gene_types nucleotide);
void init_PSA_Components(PSA_Components *PSA_variable);
void deinit_PSA_Components(PSA_Components *PSA_variable);

#endif //MSA_H