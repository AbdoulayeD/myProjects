CC=gcc
CFLAGS=-Wall -O3

LDFLAGS=-lm 

OBJS=driver.o rdtsc.o fluid_mihps.o

test: $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

fluid_mihps.o: fluid_mihps.c
	$(CC) -c $< $(CFLAGS)

clean:
	rm -rf $(OBJS) test
