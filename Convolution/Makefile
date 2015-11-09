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
	$(CEX) -n 4 ./$(target).exe > $(output)
seq:
	export OMP_NUM_THREADS=4:
	$(CEX) -n 1 ./$(target).exe > $(output)
clean:
	$(RM) *.exe
	$(RM) output.ppm
	$(RM) output.pgm
	$(RM) *.log


