#!/bin/bash
gfortran -Wall -g -c -o fortran_hypot.o  hypot.f
echo "Build Fortran"
gcc      -Wall -g -c -o hypot.o         hypot_example.c 
echo "Build C"
gcc      -Wall -g    -o hypot           hypot.o fortran_hypot.o -lm -lgfortran
echo "Linked"
