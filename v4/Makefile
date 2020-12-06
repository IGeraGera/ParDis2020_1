CC=gcc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
MPICC=mpicc
CFLAGS=-O3 -std=c99

all: v4 v4Cilk v4Omp v4pthread

v4: v4.o mmio.o
	$(CC) $(CFLAGS) -o $@ v4.o mmio.o

v4Cilk: v4Cilk.o mmio.o
	$(CILKCC) $(CFLAGS) -o $@ v4Cilk.o mmio.o -fcilkplus -pthread

v4Omp: v4Omp.o mmio.o
	$(CC) $(CFLAGS) -o $@ v4Omp.o mmio.o -fopenmp

v4OmpMultiple: v4OmpMultiple.o mmio.o
	$(CC) $(CFLAGS) -o $@ v4OmpMultiple.o mmio.o -fopenmp

v4pthread: v4pthread.o mmio.o
	$(CC) $(CFLAGS) -o $@ v4pthread.o mmio.o -pthread -lm

v4.o: v4.c mmio.h
	$(CC) $(CFLAGS) -c v4.c 

v4Cilk.o: v4Cilk.c mmio.h
	$(CILKCC) $(CFLAGS) -c  v4Cilk.c -fcilkplus -pthread

v4Omp.o: v4Omp.c mmio.h
	$(CC) $(CFLAGS) -c v4Omp.c -fopenmp

v4OmpMultiple.o: v4OmpMultiple.c mmio.h
	$(CC) $(CFLAGS) -c v4OmpMultiple.c -fopenmp

v4pthread.o: v4pthread.c mmio.h
	$(CC)  $(CFLAGS) -c v4pthread.c -pthread -lm

mmio.o: mmio.c
	$(CC) $(CFLAGS) -c mmio.c

DATAPATH=dataset/dblp-2010/dblp-2010.mtx

test: 
	./v4 $(DATAPATH)
	./v4Cilk $(DATAPATH)
	./v4Omp $(DATAPATH)
	./v4pthread $(DATAPATH) 10
.PHONY: clean

clean: 
	-rm -f mmio.o
	-rm -f v4.o
	-rm -f v4
	-rm -f v4Cilk.o
	-rm -f v4Cilk
	-rm -f v4Omp.o
	-rm -f v4Omp
	-rm -f v4OmpMultiple.o
	-rm -f v4OmpMultiple
	-rm -f v4pthread.o
	-rm -f v4pthread
