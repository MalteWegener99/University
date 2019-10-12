/* Example for Scientific Programming : how to allocate a 
 * 3D matrix (method 2 : non-contiguous memory regions).
 * 
 * Kees Lemmens, May 2013, March 2018
 */

#include <stdio.h>
#include <stdlib.h>

void show(double ***array, int m, int n, int p)
{
   int i,j,k;
   
   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
         for(k=0;k<p;k++)
            printf("array[%d][%d][%d] = %lf\n",i,j,k,array[i][j][k]);
}

int main()
{
   double*** mat3;
   int m,n,p;
   int i,j,k;
   
   m=2; n=3; p=4; /* unequal dimensions should work fine as well */
   
   mat3 = (double***) malloc(m*sizeof(double**));
   if (mat3 == NULL)
    exit(1);
   
   for (i=0; i<m; i++)
   {
      mat3[i] = (double**) malloc(n*sizeof(double*));;
      if (mat3[i] == NULL)
          exit(1);
      
      for (j=0; j<n; j++)
      {
          mat3[i][j] = (double*) malloc(p*sizeof(double));
          if (mat3[i][j] == NULL)
             exit(1);
      }
   }
   
   // Fill with some data that makes checking easy (value = sum of the indexes) :
   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
         for(k=0;k<p;k++)
            mat3[i][j][k] = i + 10*j + 100*k;
   
   show(mat3,m,n,p);
   for (i=0; i<m; i++)
   {  
      for (j=0; j<n; j++)
      {
         free(mat3[i][j]);
      }
   }

   for (i=0; i<m; i++)
   {  
      free(mat3[i]);
   }
   free(mat3);
   
   exit(0);
}
