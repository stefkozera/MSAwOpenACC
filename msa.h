#ifndef MSA_H
#define MSA_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

//Size of gene, compared during needleman-wunsch 
#define GENE_SIZE_X 4
#define GENE_SIZE_Y 4

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
    G,
    A,
    T,
    C,
    GENE_COUNT
} Gene_types;


#undef ACC_ON
#define PRINT_CHECKS

#endif //MSA_H
