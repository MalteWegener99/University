#!/bin/bash
gfortran -Wall -g -c -o fortran_array_routine.o  fortran_array_routine.f
gcc      -Wall -g -c -o fortran_arrays.o         fortran_arrays.c 
gcc      -Wall -g    -o fortran_arrays           fortran_arrays.o fortran_array_routine.o -lm -lgfortran
