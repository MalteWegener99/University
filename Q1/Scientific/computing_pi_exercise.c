// Exercise for loss of accuracy problems in floating point numbers (catastrophic cancellation).
// This uses the problematic formula as used in the example on the Pi Page by Jim Loy.
// 
// Scientific Programming, Kees Lemmens, Feb 2014, Sep 2018

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Note that the recurrence relation :
 * 
 *                (sqrt(1 + t_n^2) - 1) 
 *  t_{n+1} =     ( ----------------  )  (with typo)
 *                (        t_n        )
 * 
 * with t0 = 1/2 also yields an approximation for Pi : 
 * 
 * pi = 6 * 2^n * t_n
 * 
 * The higher the value of n, the better the approximation.
 * 
 */

int main(int argc, char **argv)
{
   double t1,t2,pow2;
   double pi1,pi2;
   int x;
   int n = 10;
   
   if(argc > 1)
     n = atoi(argv[1]);
   else
     printf("You can provide the number of iterations as argument, for now using n=%d !\n\n",n);
   
   t1 = t2 = 1./2; // modify
   pow2 = 1; // modify

   for(x=0;x<n;x++)
   {
      pi1 = 6 * pow2 * t1;
      pi2 = 6 * pow2 * t2;
      pow2 *= 2;
      printf("Step %d : result for method3 = %.16lf  ,  method4 = %.16lf \n",x, pi1, pi2);
      t1  = (sqrt(1+t1*t1)-1)/t1;
      t2  = t2/(sqrt(t2*t2+1)+1);
   }

   exit(0);
}
