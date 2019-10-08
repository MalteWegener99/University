// Example for overflow problems in floating point numbers
// Scientific Programming, Kees Lemmens, Feb 2013

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Compare e.g.
 * ./hypot_example-float 3e154 4e154 (normal overflows, safer doesn't)
 * ./hypot_example-float 3e307 4e307 (normal overflows, safer doesn't)
 * ./hypot_example-float 3e153 4e153 (both work fine)
 * ./hypot_example-float 3e308 4e308 (both fail now)
 * 
*/

void hypot_normal(double x, double y)
{
   double result;
   
   result = sqrt( x*x + y*y );

   printf("Result for normal       method (x=%lg, y = %lg) = %lg\n",x,y,result);
}

void hypot_safer_fabs(double x, double y)
{
   double result;

   if(x > y)
     result = fabs(x) * sqrt( 1.0 + (y/x)*(y/x));
   if(x < y)
     result = fabs(y) * sqrt( 1.0 + (x/y)*(x/y));
   if(x == y)
     result=sqrt(2) * fabs(x);
   
   printf("Result for safe  (fabs) method (x=%lg, y = %lg) = %lg\n",x,y,result);
}

int main(int argc, char **argv)
{
   double x,y;
   
   if(argc > 2)
   {
      x = atof(argv[1]);
      y = atof(argv[2]);
   }
   else
   {
      x=3e10;
      y=4e10;
      printf("You can provide 2 arguments, for now using x=%lg, x=%lg !\n",x,y);
   }
   
   hypot_normal(x,y);
   hypot_safer_fabs(x,y);
   
   exit(0);
}
