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

par:
	export OMP_NUM_THREADS=4
	$(CEX) -n 8 ./$(target).exe> $(output)

seq:
	$(CEX) -n 1 ./$(target).exe> $(output)

parppm:
	export OMP_NUM_THREADS=4
	$(CEX) -n 8 ./$(target).exe ./img/football.ppm ppm> $(output)

parpgm:
	export OMP_NUM_THREADS=4
	$(CEX) -n 10 ./$(target).exe ./img/cake.pgm pgm> $(output)

seqppm:
	export OMP_NUM_THREADS=4:
	$(CEX) -n 1 ./$(target).exe ./img/football.ppm ppm > $(output)

seqpgm:
	export OMP_NUM_THREADS=4:
	$(CEX) -n 1 ./$(target).exe ./img/cake.pgm pgm > $(output)


clean:
	$(RM) *.exe
	$(RM) output.ppm
	$(RM) output.pgm
	$(RM) *.log

