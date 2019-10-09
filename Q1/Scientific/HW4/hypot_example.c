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

extern void hypot_normal_(double* x, double *y, double* r);

extern void hypot_safer_fabs_(double* x, double *y, double* r);

int main(int argc, char **argv)
{
   double x,y,r1,r2;
   
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
   
   hypot_normal_(&x, &y, &r1);
   hypot_safer_fabs_(&x, &y, &r2);

   printf("Normal Method: %f\n", r1);
   printf("Safer Method: %f\n", r2);

   exit(0);
}
