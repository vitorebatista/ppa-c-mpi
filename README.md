# Matriz MPI

Este programa em sua execução realiza operações básicas de multiplição de estruturas matriciais de quatro diferentes maneiras: sequencial, em blocos, mpi e mpi em blocos. 

Estão presentes nos códigos funções para alocação de memória, geração e preenchimento de todas as posições das matrizes, liberação da memória, etc.

## Compilação

Estando dentro do diretório correto basta executar o *Makefile*.

> $ make

Após execução deste comando os códigos serão compilados e gerados os executáveis e binarios.

## Excecução

Deve-se gerar os arquivos com as matrizes para realizar os testes:

> ./gmat 1000 1000

Em seguida chamar o programa principal que realizará os devidos testes e apresentará os resultados obtidos. Haverá dois parâmetros, o primeiro e segundo serão do arquivo de matriz.

> ./mpirun -np 2 main_mpi 1000x1000-mat.map 1000x1000-mat.map

Será apresentado no terminal um retorno semelhante a este:

```

        COMPARAR MATRIZ_SeqC c/ MATRIZ_SeqBlC
        Matrizes são idênticas!! :) 

        COMPARAR MATRIZ_SeqC c/ MATRIZ_MPIC
        Matrizes são idênticas!! :) 

        COMPARAR MATRIZ_SeqC c/ MATRIZ_MPIBlC
        Matrizes são idênticas!! :) 

        Tempo Médio MATRIZ_SeqC:        0.569333 sec 
        Tempo Médio MATRIZ_SeqBlC:      3.339507 sec
        Tempo Médio MATRIZ_MPIC:        2.420058 sec 
        Tempo Médio MATRIZ_MPIBlC:      1.956970 sec 

        SPEEDUP (MATRIZ_C):     0.235 (23.53 %)
        SPEEDUP (MATRIZ_BLC):   1.706 (170.65 %)

```