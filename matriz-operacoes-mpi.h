#include "toolsv3.h"
#include "matrizv3.h"
#include "matriz-operacoesv3.h"

void multiplicarMPI(int *A, int *B, int *C, int MPI_size, int MPI_rank, int chunk_lines_A, int size);
void multiplicaBlocoMPI(int *A, int *B, int *C, int MPI_size, int MPI_rank, int chunk_lines_A, int size, int N);