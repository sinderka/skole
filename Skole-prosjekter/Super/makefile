all: mpi single openmp superkjapp openmpi

mpi: mpi.c
	mpicc -g -O3 mpi.c -o mpi

single: single.c
	gcc -o single single.c

openmp: openmp.c
	gcc -o openmp openmp.c

superkjapp: superkjapp.c
	mpicc -g -O3 superkjapp.c -o superkjapp

openmpi: openmpi.c
	mpicc -g -O3 openmpi.c -o openmpi
