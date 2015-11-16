CC=mpicxx
CEX=mpiexec
RM=rm
target=convolution
main=main.cpp
CFLAGS= -fopenmp -std=c++11 -o
output= convolution.log

all: compile

compile:
	$(CC) $(CFLAGS) $(target).exe $(main)

parppm:
	export OMP_NUM_THREADS=4
	$(CEX) -n 4 ./$(target).exe ./img/flower.ppm ppm> $(output)

parpgm:
	export OMP_NUM_THREADS=4
	$(CEX) -n 4 ./$(target).exe ./img/lena.pgm pgm> $(output)

seqppm:
	export OMP_NUM_THREADS=4:
	$(CEX) -n 1 ./$(target).exe ./img/football.ppm ppm > $(output)

seqpgm:
	export OMP_NUM_THREADS=4:
	$(CEX) -n 1 ./$(target).exe ./img/dragon.pgm pgm > $(output)


clean:
	$(RM) *.exe
	$(RM) output.ppm
	$(RM) output.pgm
	$(RM) *.log


