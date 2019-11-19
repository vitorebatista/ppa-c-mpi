#include "matriz-operacoes-mpi.h"

#include <mpi.h>
#define MASTER 0

/* Variáveis globais*/
mymatriz mat_a, mat_b;
mymatriz *mmult_MATRIZ_SeqC;
mymatriz *mmult_MATRIZ_SeqBlC;
mymatriz *mmult_MATRIZ_MPIC;
mymatriz *mmult_MATRIZ_MPIBlC;
int *A, *B, *C;
int *A1, *B1, *C1;

int main(int argc, char *argv[])
{

    int MPI_rank, MPI_size;

	MPI_Init(&argc, &argv);
    /* Número de processos */
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

    /* Captura o número do processo em execução */
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);

    // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
    // DECLARAÇÃO de VARIÁVEIS

    // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
    // DECLARAÇÃO DE VARIÁVEIS


    //variaveis para manipulacao de arquivos no disco
    char filename[100];
    FILE *fmat;
    int nr_line;
    int *vet_line = NULL;
    int N, M, La, Lb;

    //variaveis para controle de blocos
    matriz_bloco_t **Vsubmat_a = NULL;
    matriz_bloco_t **Vsubmat_b = NULL;
    matriz_bloco_t **Vsubmat_c = NULL;
    
    //For para executar calculo da média
    int n_threads = 4;
    int nro_submatrizes = n_threads;
    int count_for = 10; //numero de repeticoes para média de runtime

    //variaveis para controle de tempo (runtime)
    double start_time, end_time;
    double tempo_MATRIZ_SeqC = 0;
    double tempo_MATRIZ_SeqBlC = 0;
    double tempo_MATRIZ_MPIC = 0;
    double tempo_MATRIZ_MPIBlC = 0;
    double speedup_seqC;
    double speedup_BlC;
    // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%

    if (argc < 3)
    {
        printf("ERRO: Numero de parametros %s <matriz_a> <matriz_b> <threads>\n", argv[0]);
        exit(1);
    }

    if (argv[3] != NULL){
        nro_submatrizes = atoi(argv[3]);
        n_threads = atoi(argv[3]);
    }

    // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
    //Leitura da Matriz A (arquivo)
    fmat = fopen(argv[1], "r");
    if (fmat == NULL)
    {
        printf("Error: Na abertura dos arquivos.");
        exit(1);
    }
    extrai_parametros_matriz(fmat, &N, &La, &vet_line, &nr_line);
    mat_a.matriz = NULL;
    mat_a.lin = N;
    mat_a.col = La;
    if (malocar(&mat_a))
    {
        printf("ERROR: Out of memory\n");
    }
    filein_matriz(mat_a.matriz, N, La, fmat, vet_line, nr_line);
    free(vet_line);
    fclose(fmat);
    // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%

    // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
    //               Leitura da Matriz B (arquivo)
    fmat = fopen(argv[2], "r");
    if (fmat == NULL)
    {
        printf("Error: Na abertura dos arquivos.");
        exit(1);
    }
    extrai_parametros_matriz(fmat, &Lb, &M, &vet_line, &nr_line);
    mat_b.matriz = NULL;
    mat_b.lin = Lb;
    mat_b.col = M;
    if (malocar(&mat_b))
    {
        printf("ERROR: Out of memory\n");
    }
    filein_matriz(mat_b.matriz, Lb, M, fmat, vet_line, nr_line);
    free(vet_line);
    fclose(fmat);



    // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
    // Multiplicação Sequencial
    if(MPI_rank == MASTER){
        mmult_MATRIZ_SeqC = (mymatriz *)malloc(sizeof(mymatriz));
        for (int count = 0; count < count_for; count++)
        {   
            start_time = wtime();
            printf("\rMultiplicação Sequencial, teste %d...             ", count+1);
            fflush(stdout);

            mmult_MATRIZ_SeqC = mmultiplicar(&mat_a, &mat_b, 3);  //1=mais rápido (2.04), 5=mais lento (5.94)
            end_time = wtime();
            tempo_MATRIZ_SeqC += end_time - start_time;
        }
        sprintf(filename, "MATRIZ_SeqC.result");
        fmat = fopen(filename, "w");
        fileout_matriz(mmult_MATRIZ_SeqC, fmat);
        fclose(fmat);
        // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%

        // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
        // Multiplicação Sequencial em Bloco
        printf("\n");
        mmult_MATRIZ_SeqBlC = (mymatriz *)malloc(sizeof(mymatriz));
        for (int count = 0; count < count_for; count++)
        {
            start_time = wtime();
            printf("\rMultiplicação Sequencial em Bloco, teste %d...             ", count+1);
            fflush(stdout);

            Vsubmat_a = particionar_matriz(mat_a.matriz, N, La, 1, nro_submatrizes);
            Vsubmat_a = particionar_matriz(mat_a.matriz, N, La, 1, nro_submatrizes);
            Vsubmat_b = particionar_matriz(mat_b.matriz, Lb, M, 0, nro_submatrizes);
            Vsubmat_c = csubmatrizv2(N, M, nro_submatrizes);
            
            //multiplicacao de blocos
            for (int i = 0; i < nro_submatrizes; i++){
                multiplicar_submatriz (Vsubmat_a[i], Vsubmat_b[i], Vsubmat_c[i]);
            }

            //soma os blocos separados
            mmult_MATRIZ_SeqBlC = msomar(Vsubmat_c[0]->matriz,Vsubmat_c[1]->matriz, 1);
            mmult_MATRIZ_SeqBlC = msomar(Vsubmat_c[0]->matriz,Vsubmat_c[1]->matriz, 1);
            for (int i = 2; i < nro_submatrizes; i++){
                mmult_MATRIZ_SeqBlC = msomar(mmult_MATRIZ_SeqBlC,Vsubmat_c[i]->matriz, 1);	
            }

            end_time = wtime();
            tempo_MATRIZ_SeqBlC += end_time - start_time;
        }
        sprintf(filename, "MATRIZ_SeqBlC.result");
        fmat = fopen(filename, "w");
        fileout_matriz(mmult_MATRIZ_SeqBlC, fmat);
        fclose(fmat);
    }
    // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%

    // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
    // Multiplicação MPI
    int size = N;

	size = size % MPI_size == 0 ? size : size + (MPI_size - size % MPI_size);

	int chunk_lines_A = (int) (size / MPI_size);

	if(MPI_rank == MASTER){
		A = calloc(size * size, sizeof(int));
		B = calloc(size * size, sizeof(int));
		C = calloc(size * size, sizeof(int));

		for(int i = 0; i < size; i++){
			for(int j = 0; j < size; j++){
					A[i * size + j] = mat_a.matriz[i][j];
					B[i * size + j] = mat_b.matriz[i][j];
				
			}
		}
        start_time = wtime();
	}else{
		A = malloc(chunk_lines_A * size * sizeof(int));
		B = malloc(size * size * sizeof(int));
		C = calloc(chunk_lines_A * size, sizeof(int));
	}

    printf("\n");
    mmult_MATRIZ_MPIC = (mymatriz *)malloc(sizeof(mymatriz));
    mmult_MATRIZ_MPIC = malloc(sizeof(mymatriz));
    mmult_MATRIZ_MPIC->matriz = NULL;
    mmult_MATRIZ_MPIC->lin = mat_a.lin;
    mmult_MATRIZ_MPIC->col = mat_b.col;

    //realiza a alocação de memória para matriz resultado
    if (malocar(mmult_MATRIZ_MPIC)) {
        printf("ERROR: Out of memory\n");
        exit(1);
    }else{
        mzerar(mmult_MATRIZ_MPIC);
    }

    printf("\rMultiplicação MPI, rank: %d...             ", MPI_rank);
    fflush(stdout);

    //mzerar(mmult_MATRIZ_MPIC);
    
    multiplicarMPI (A, B, C, MPI_size, MPI_rank, chunk_lines_A, size);	
    
        

    if(MPI_rank == MASTER){
		for(int i = 0; i < size; i++){
			for(int j = 0; j < size; j++){
				mmult_MATRIZ_MPIC->matriz[i][j] = C[i * size + j];
			}
		}
        end_time = wtime();
        tempo_MATRIZ_MPIC += end_time - start_time;
        sprintf(filename, "MATRIZ_MPIC.result");
        fmat = fopen(filename, "w");
        fileout_matriz(mmult_MATRIZ_MPIC, fmat);
        fclose(fmat);
    }
    // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%

    // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
    // Multiplicação MultiThreads em Bloco
    printf("\n");
    if(MPI_rank == MASTER){
		A1 = calloc(size * size, sizeof(int));
		B1 = calloc(size * size, sizeof(int));
		C1 = calloc(size * size, sizeof(int));

		for(int i = 0; i < size; i++){
			for(int j = 0; j < size; j++){
					A1[i * size + j] = mat_a.matriz[i][j];
					B1[i * size + j] = mat_b.matriz[i][j];
				
			}
		}
        start_time = wtime();
	}else{
		A1 = malloc(chunk_lines_A * size * sizeof(int));
		B1 = malloc(size * size * sizeof(int));
		C1 = calloc(chunk_lines_A * size, sizeof(int));
	}

    mmult_MATRIZ_MPIBlC = (mymatriz *)malloc(sizeof(mymatriz));
    mmult_MATRIZ_MPIBlC = malloc(sizeof(mymatriz));
    mmult_MATRIZ_MPIBlC->matriz = NULL;
    mmult_MATRIZ_MPIBlC->lin = mat_a.lin;
    mmult_MATRIZ_MPIBlC->col = mat_b.col;

    //realiza a alocação de memória para matriz resultado
    if (malocar(mmult_MATRIZ_MPIBlC)) {
        printf("ERROR: Out of memory\n");
        exit(1);
    }else{
        mzerar(mmult_MATRIZ_MPIBlC);
    }
    printf("\rMultiplicação MPI em bloco, rank: %d...             ",MPI_rank);
    fflush(stdout);

    //Vsubmat_a = particionar_matriz(mat_a.matriz, N, La, 1, nro_submatrizes);
    //Vsubmat_b = particionar_matriz(mat_b.matriz, Lb, M, 0, nro_submatrizes);
    //Vsubmat_c = csubmatrizv2(N, M, nro_submatrizes);

    
    multiplicaBlocoMPI (A1,B1,C1, MPI_size, 2, MPI_rank, chunk_lines_A, size, N);
    

    if(MPI_rank == MASTER){
		for(int i = 0; i < size; i++){
			for(int j = 0; j < size; j++){
				mmult_MATRIZ_MPIBlC->matriz[i][j] = C1[i * size + j];
                //printf("\nvalor: %d", C1[i * size + j]);
			}
		}		
        end_time = wtime();
    
        tempo_MATRIZ_MPIBlC += end_time - start_time;
        sprintf(filename, "MATRIZ_MPIBlC.result");
        fmat = fopen(filename, "w");
        fileout_matriz(mmult_MATRIZ_MPIBlC, fmat);
        fclose(fmat);
    
    // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%
    
        // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
        // Impressao dos resultados de tempo
        printf("\n\n\tCOMPARAR MATRIZ_SeqC c/ MATRIZ_SeqBlC\n\t");
        mcomparar(mmult_MATRIZ_SeqC, mmult_MATRIZ_SeqBlC);

        printf("\n\tCOMPARAR MATRIZ_SeqC c/ MATRIZ_MPIC\n\t");
        mcomparar(mmult_MATRIZ_SeqC, mmult_MATRIZ_MPIC);
        
        printf("\n\tCOMPARAR MATRIZ_SeqC c/ MATRIZ_MPIBlC\n\t");
        mcomparar(mmult_MATRIZ_SeqC, mmult_MATRIZ_MPIBlC);

        printf("\n\tTempo Médio MATRIZ_SeqC:\t%.6f sec \n", tempo_MATRIZ_SeqC / count_for);
        printf("\tTempo Médio MATRIZ_SeqBlC:\t%.6f sec\n", tempo_MATRIZ_SeqBlC / count_for );
        printf("\tTempo Médio MATRIZ_MPIC:\t%.6f sec \n", tempo_MATRIZ_MPIC / count_for);
        printf("\tTempo Médio MATRIZ_MPIBlC:\t%.6f sec \n", tempo_MATRIZ_MPIBlC / count_for);

        speedup_seqC = (tempo_MATRIZ_SeqC / count_for) / (tempo_MATRIZ_MPIC / count_for);
        speedup_BlC = (tempo_MATRIZ_SeqBlC / count_for) / (tempo_MATRIZ_MPIBlC / count_for);
        printf("\n\tSPEEDUP (MATRIZ_C): \t%.3f (%.2f %c)", speedup_seqC, speedup_seqC*100, 37 );
        printf("\n\tSPEEDUP (MATRIZ_BLC): \t%.3f (%.2f %c)\n\n", speedup_BlC, speedup_BlC*100, 37 );

        // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%

        // %%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%
        //Liberação de memória
        mliberar(mmult_MATRIZ_SeqC);
        mliberar(mmult_MATRIZ_SeqBlC);
        mliberar(mmult_MATRIZ_MPIC);
        mliberar(mmult_MATRIZ_MPIBlC);

        free(mmult_MATRIZ_SeqC);
        free(mmult_MATRIZ_SeqBlC);
        free(mmult_MATRIZ_MPIC);
        free(mmult_MATRIZ_MPIBlC);

        mliberar(&mat_a);
        mliberar(&mat_b);
    }
    // %%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%
    
    MPI_Finalize();	  
    return 0;
}