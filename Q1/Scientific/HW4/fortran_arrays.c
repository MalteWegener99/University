// This C program calls a routine from Fortran
// To make this possible 3 things are important :
// 
// 1: the Fortran function name should end with "_" (if using GNU Fortran)
// 2: the arguments should be passed by reference (using a pointer)
// 3: the return value should be also passed as an argument

// Steps to compile this :
//   
//   gfortran -Wall -g -c -o fortran_array_routine.o  fortran_array_routine.f
//   gcc      -Wall -g -c -o fortran_arrays.o         fortran_arrays.c 
//   gcc      -Wall -g    -o fortran_arrays           fortran_arrays.o fortran_array_routine.o -lm -lgfortran


#include <stdio.h>
#include <stdlib.h>

// Declare the calling method for the Fortran routine :
extern void multiply_array_(int *a, int *b, int *c, int *n);

#define ASIZE 10

int main()
{
   int a[ASIZE], b[ASIZE], c[ASIZE];
   int n = ASIZE;
   int x;
   
   // Assign some arbitrary values
   for(x=0;x<n;x++)
   {
      a[x] = x+2;
      b[x] = x*2;
   }
   
   // Call Fortran routine :
   // Add the correct arguments here : don't forget to pass the size
   multiply_array_(a, b, c, &n);
   
   for (x = 0;x < n;x++)
   {
      printf("%d * %d = %d\n",a[x], b[x], c[x]);
   }
   
   exit(0);
}
