# GNU Makefile
# Revisão: ago/2019

CC = mpicc
CCFLAGS = -Wall -O3 
LDFLAGS = 
TARGET = main_mpi gmat help

all: $(TARGET)

%.o: %.c
	$(CC) $(CCFLAGS) -c $<

%: %.o
	$(CC) $(LDFLAGS) $^ -o $@

test: test.c 
main_mpi: main_mpi.c matrizv3.o toolsv3.o matriz-operacoesv3.o matriz-operacoes-mpi.o
			$(CC) $(CCFLAGS) matriz-operacoesv3.o matrizv3.o toolsv3.o matriz-operacoes-mpi.o main_mpi.c -o $@ $(LDFLAGS)

gmat: matrizv3.o toolsv3.o gera_matrizv3.c
		$(CC) $(CCFLAGS) matrizv3.o toolsv3.o gera_matrizv3.c -o $@ $(LDFLAGS)

help:
	@echo
	@echo
	@echo "####### Exemplo de Execução #######"
	@echo "./gmat 1000 1000"
	@echo "mpirun -np 2 main_mpi 1000x1000-mat.map 1000x1000-mat.map"

clean:
	rm -f *.o *~ $(TARGET) *.map *.result
