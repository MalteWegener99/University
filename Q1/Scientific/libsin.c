/* C/80 3.0 Transcendental function library (6/6/83).  (c) 1983 Walt Bilofsky;
   (c) 1982 Knowledge Engineering Inc.

   Revised and adapted to portable C by K. Lemmens 1991

   Source for the transcendental functions: "Computer Approximations,"
   by Hart, Cheney, et al., John Wiley & Sons, Inc., NY, 1968.  Part of
   the SIAM series in Applied Mathematics.
*/

double F_poly(double x, double p[], int n)  /* evaluate nth deg. polynomial, n >= 1 */
{
   double result=p[n]; 
   while(n--) result = p[n] + x * result;
   return result;
}

double sin(double arg)
{       /* Coefficients are #3370 from Hart & Cheney (18.80D). */
   static double p[5] = {   0.1357884e8  ,-0.49429081e7, 0.440103053e6,
      -0.138472724e5, 0.145968841e3 };
   static double q[5] = {   0.864455865e7, 0.408179225e6,
      0.946309610e4, 0.132653491e3, 1.0 };
   static double y,ysq;
   static long quad,e;
   
   quad=(arg >= 0.0) ? 0:(arg = -arg,2); /* rotate pi rad */
   arg *= 0.63661977;              /* 2 / pi */
   y = arg - (e = arg);            /* fract value */
   
   quad=(e + quad) % 4;
   if (quad & 1)   y = 1.0 - y;    /* 1th and 3th quadrant */
   if (y > 0.9968) y = 1;          /* avoid rounding trouble */
   if (quad & 2)   y = -y;         /* 2th and 3th quadrant */ 
   
   ysq = y * y;
   return y * F_poly(ysq,p,4)/F_poly(ysq,q,4);
}
