//
// Allocate 3D matrix - efficient method
// Example as given by Niels Looye (omission in coursebook Wi4260, p. 95)
 
#include <stdio.h>
#include <stdlib.h>

int main()
{
   double*** mat3;
   int m,n,p;
   int i,j,k;
   
   m=2; n=3; p=4;

   printf("This program causes a coredump, unless you fix the error in it !\n");
   
   mat3 = (double***) malloc(m*sizeof(double**));
   if (mat3 == NULL)
      exit(1);
   
   mat3[0] = (double**) malloc(m*n*sizeof(double*));
   if (mat3[0] == NULL)
      exit(1);
   
   mat3[0][0] = (double*) malloc(m*n*p*sizeof(double));
   if (mat3[0][0] == NULL)
      exit(1);
   
   for (i=0; i<m; i++)
   {
      mat3[i] = mat3[0] + i*n;
      // fix is here
      mat3[i][0] = mat3[0][0] + i*n*p;
      for(j=0; j<n; j++)
      {
         mat3[i][j] = mat3[i][0] + j*p; // dit gaat fout
      }
   }
   
   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
         for(k=0;k<p;k++)
             mat3[i][j][k] = i + 10*j + 100*k;
   
   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
         for(k=0;k<p;k++)
             printf("mat3[%d][%d][%d] = %lf\n",i,j,k,mat3[i][j][k]);
   
   exit(0);
}
