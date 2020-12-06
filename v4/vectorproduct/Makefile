CC=gcc
CFLAGS=-O2

vectorproduct:
	$(CC) $(CFLAGS) -o vectorproduct vectorproductFun.c mmio.c

test:
	python3 testGenerator.py
	./vectorproduct test.mtx

.PHONY: clean

clean: 
	-rm test.mtx 
	-rm vectorproduct
