# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7
MPICC = mpicc
MPIRUN = mpirun -np 4

test_sequential:
	tar -xvzf code.tar.gz
	cd knnring; make lib; cd ..
	cd knnring; cp lib/*.a inc/knnring.h ../; cd ..
	$(CC) tester.c knnring_sequential.a -o $@ -lm -lblas
	./test_sequential


test_mpi_sync:
	tar -xvzf code.tar.gz
	cd knnring; make lib; cd ..
	cd knnring; cp lib/*.a inc/knnring.h ../; cd ..
	$(MPICC) tester_mpi.c knnring_synchronous.a -o $@ -lm -lblas
	$(MPIRUN) ./test_mpi_sync

test_mpi_async:
	tar -xvzf code.tar.gz
	cd knnring; make lib; cd ..
	cd knnring; cp lib/*.a inc/knnring.h ../; cd ..
	$(MPICC) tester_mpi.c knnring_asynchronous.a -o $@ -lm -lblas
	$(MPIRUN) ./test_mpi_async

