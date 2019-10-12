#include <stdio.h>
#include <stdlib.h>

#define DIMX 4
#define DIMY 5

/*
 * Dynamic matrices in C with dimension 2 (or higher)
 * 
 * Method 2 : Here we malloc an array of pointers with size Y. After this
 * we can malloc separate blocks of size X and assign their addresses to
 * the pointers in our first array by using a simple loop.
 * 
 * Kees Lemmens; TU Delft 1996
 */ 

void show(double **array)
{  int x,y;
   
   for(y=0;y<DIMY;y++)
   {
      for(x=0;x<DIMX;x++)
         printf("%f ",array[y][x]);
      printf("\n");
   }
}

int main()
{  
   double **array;
   int x,y;
   
   array = (double **)malloc(DIMY*sizeof(double *));
   if(array == NULL)
      exit(1);
   
   for(y=0;y<DIMY;y++)
   {
      array[y] = (double *)malloc(DIMX*sizeof(double));
      if(array[y] == NULL)
         exit(1);
   }
   
   // Fill with some data that makes checking easy (value = sum of the indexes) :
   for(y=0;y<DIMY;y++)
      for(x=0;x<DIMX;x++)
         array[y][x]=x+y;
   
   show(array);
   free(array);
   exit(0);
}
