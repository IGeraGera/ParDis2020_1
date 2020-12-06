CC=gcc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3 -std=c99

default: all


v3Cilk: v3Cilk.o mmio.o
	$(CILKCC) $(CFLAGS) -o $@ v3Cilk.o mmio.o -fcilkplus -lpthread

v3Omp: v3Omp.o mmio.o
	$(CC) $(CFLAGS) -o $@ v3Omp.o mmio.o -fopenmp

v3OmpMultiple: v3OmpMultiple.o mmio.o
	$(CC) $(CFLAGS) -o $@ v3OmpMultiple.o mmio.o -fopenmp

v3: v3.o mmio.o
	$(CC) $(CFLAGS) -o $@ v3.o mmio.o   

v3.o:	v3.c mmio.h
	$(CC) $(CFLAGS) -c v3.c 

v3Omp.o: v3Omp.c mmio.h
	$(CC) $(CFLAGS) -c v3Omp.c -fopenmp

v3OmpMultiple.o: v3OmpMultiple.c mmio.h
	$(CC) $(CFLAGS) -c v3OmpMultiple.c -fopenmp

v3Cilk.o: v3Cilk.c mmio.h
	$(CC) $(CFLAGS) -c v3Cilk.c -fcilkplus -pthread

mmio.o: mmio.c
	$(CC) $(CFLAGS) -c mmio.c 

.PHONY: clean

all:	v3 v3Omp v3Cilk v3OmpMultiple

test: 
	./v3 dataset/dblp-2010/dblp-2010.mtx
	./v3Cilk dataset/dblp-2010/dblp-2010.mtx
	./v3Omp dataset/dblp-2010/dblp-2010.mtx

clean:
	-rm -f v3Cilk
	-rm -f v3Cilk.o	
	-rm -f v3Omp 
	-rm -f v3Omp.o
	-rm -f v3OmpMultiple 
	-rm -f v3OmpMultiple.o
	-rm -f v3
	-rm -f v3.o
	-rm -f mmio.o
