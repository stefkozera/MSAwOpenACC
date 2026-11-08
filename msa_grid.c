#include "msa_grid.h"

void InitRandomGene(Gene_types *gene_X_arr, Gene_types *gene_Y_arr)
{
    for(int i = 0; i < GENE_SIZE_X; ++i){
        gene_X_arr[i] = (Gene_types)(rand() % GENE_COUNT);
    }

    for(int i = 0; i < GENE_SIZE_Y; ++i){
        gene_Y_arr[i] = (Gene_types)(rand() % GENE_COUNT);
    }
}

void InitGapPenalty(int *matrix)
{
    for(int i = 0; i < GRID_SIZE_Y; ++i)
    {
        ACCESS_2D(matrix, i, 0, GRID_SIZE_X) = i * GAP_PENALTY_VAL;
    }

    for(int i = 0; i < GRID_SIZE_X; ++i)
    {
        ACCESS_2D(matrix, 0, i, GRID_SIZE_X) = i * GAP_PENALTY_VAL;
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