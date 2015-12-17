all: compress exec

compress:
	gcc  -o distanceM.exe distanceM.c

exec:
	export OMP_NUM_THREADS=4
	./distanceM.exe NC_0*.fna

clean:
	rm *.zip
	rm *.exe
