/*
 * file for ep.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

typedef int boolean;
typedef struct {
   double real;
   double imag;
} dcomplex;
void wtime_(double *t);

void wtime_(double *t);
double elapsed_time(void);


double start[64], elapsed[64];
void timer_clear(int n);

void timer_start(int n);

void timer_stop(int n);

double timer_read(int n);


double randlc(double *x, double a);

void vranlc(int n, double *x_seed, double a, double y[]);

void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand);


#include "npbparams.h"
static double x[2*(1<<16)];
static double q[10];

/*--------------------------------------------------------------------
      program EMBAR
c-------------------------------------------------------------------*/
/*
c   This is the serial version of the APP Benchmark 1,
c   the "embarassingly parallel" benchmark.
c
c   M is the Log_2 of the number of complex pairs of uniform (0, 1) random
c   numbers.  MK is the Log_2 of the size of each batch of uniform random
c   numbers.  MK can be set for convenience on a given system, since it does
c   not affect the results.
*/
int main(int argc, char **argv);
void wtime_(double *t)
{
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, (void *) 0);
   if (sec<0)
      sec = tv.tv_sec;
   *t = tv.tv_sec-sec+1.0e-6*tv.tv_usec;
}
double elapsed_time(void)
{
   double t;

   wtime_(&t);
   return t;
}
void timer_clear(int n)
{
   elapsed[n] = 0.0;
}
void timer_start(int n)
{
   start[n] = elapsed_time();
}
void timer_stop(int n)
{
   double t, now;

   now = elapsed_time();
   t = now-start[n];
   elapsed[n] += t;
}
double timer_read(int n)
{
   return elapsed[n];
}
double randlc(double *x, double a)
{
   double t1, t2, t3, t4, a1, a2, x1, x2, z;
   t1 = 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*a;
   a1 = (int) t1;
   a2 = a-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*a1;
   t1 = 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5**x;
   x1 = (int) t1;
   x2 = *x-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*x1;
   t1 = a1*x2+a2*x1;
   t2 = (int) (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*t1);
   z = t1-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*t2;
   t3 = 2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*z+a2*x2;
   t4 = (int) (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*t3);
   *x = t3-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*t4;
   return 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5**x;
}
void vranlc(int n, double *x_seed, double a, double y[])
{
   int i;
   double x, t1, t2, t3, t4, a1, a2, x1, x2, z;
   t1 = 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*a;
   a1 = (int) t1;
   a2 = a-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*a1;
   x = *x_seed;
   for(i = 1; i <= n; i += 1) {
      t1 = 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*x;
      x1 = (int) t1;
      x2 = x-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*x1;
      t1 = a1*x2+a2*x1;
      t2 = (int) (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*t1);
      z = t1-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*t2;
      t3 = 2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*z+a2*x2;
      t4 = (int) (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*t3);
      x = t3-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*t4;
      y[i] = 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*x;
   }
   *x_seed = x;
}
void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand)
{
   char *evalue = "1000";

   printf("\n\n %s Benchmark Completed\n", name);

   printf(" Class           =                        %c\n", class);

   if (n2==0&&n3==0)
      printf(" Size            =             %12d\n", n1);
   else
      printf(" Size            =              %3dx%3dx%3d\n", n1, n2, n3);

   printf(" Iterations      =             %12d\n", niter);

   printf(" Threads         =             %12d\n", nthreads);

   printf(" Time in seconds =             %12.2f\n", t);

   printf(" Mop/s total     =             %12.2f\n", mops);

   printf(" Operation type  = %24s\n", optype);

   if (passed_verification)
      printf(" Verification    =               SUCCESSFUL\n");
   else
      printf(" Verification    =             UNSUCCESSFUL\n");

   printf(" Version         =           %12s\n", npbversion);

   printf(" Compile date    =             %12s\n", compiletime);

   printf("\n Compile options:\n");

   printf("    CC           = %s\n", cc);

   printf("    CLINK        = %s\n", clink);

   printf("    C_LIB        = %s\n", c_lib);

   printf("    C_INC        = %s\n", c_inc);

   printf("    CFLAGS       = %s\n", cflags);

   printf("    CLINKFLAGS   = %s\n", clinkflags);

   printf("    RAND         = %s\n", rand);
}
int main(int argc, char **argv)
{

   double Mops, t1, t2, t3, t4, x1, x2, sx, sy, tm, an, tt, gc;
   double dum[3] = {1.0, 1.0, 1.0};

   int np, ierr, node, no_nodes, i, ik, kk, l, k, nit, ierrcode, no_large_nodes, np_add, k_offset, j;
   int nthreads = 1;
   boolean verified;
   char size[13+1];
   /* character*13 */
   
   /*
   c   Because the size of the problem is too large to store in a 32-bit
   c   integer for some classes, we put it into a string (for printing).
   c   Have to strip off the decimal point put in there by the floating
   c   point print statement (internal file)
   */
   
   
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"" - EP Benchmark\n");
   sprintf(size, "%12.0f", pow(2.0, 30+1));
#pragma omp parallel for
   for(j = 13; j >= 1; j += -1)
      if (size[j]=='.')
         size[j] = ' ';
   printf(" Number of random numbers generated: %13s\n", size);

   verified = 0;
   
   /*
   c   Compute the number of "batches" of random number pairs generated 
   c   per processor. Adjust if the number of processors does not evenly 
   c   divide the total number
   */
   np = 1<<30-16;
   
   /*
   c   Call the random number generator functions and initialize
   c   the x-array to reduce the effects of paging on the timings.
   c   Also, call all mathematical functions that are used. Make
   c   sure these initializations cannot be eliminated as dead code.
   */
   vranlc(0, &dum[0], dum[1], &dum[2]);
   dum[0] = randlc(&dum[1], dum[2]);

#pragma omp parallel for
   for(i = 0; i <= 2*(1<<16)-1; i += 1)
      x[i] = -1.0e99;
   Mops = log(sqrt(fabs(1.0>1.0?1.0:1.0)));

   timer_clear(1);
   timer_clear(2);
   timer_clear(3);
   timer_start(1);

   vranlc(0, &t1, 1220703125.0, x);
   
   /*   Compute AN = A ^ (2 * NK) (mod 2^46). */
   
   t1 = 1220703125.0;

   for(i = 1; i <= 16+1; i += 1)
      t2 = randlc(&t1, t1);

   an = t1;
   tt = 271828183.0;
   gc = 0.0;
   sx = 0.0;
   sy = 0.0;

#pragma omp parallel for
   for(i = 0; i <= 10-1; i += 1)
      q[i] = 0.0;
   
   /*
   c   Each instance of this loop may be performed independently. We compute
   c   the k offsets separately to take into account the fact that some nodes
   c   have more numbers to generate than others
   */
   k_offset = -1;
   {
      double t1, t2, t3, t4, x1, x2;
      int kk, i, ik, l;
      double qq[10];
      /* private copy of q[0:NQ-1] */
      
#pragma omp parallel for
      for(i = 0; i <= 9; i += 1)
         qq[i] = 0.0;
      
      //#pragma omp parallel for copyin(x, qq) private(x1, x2, t1, t2, t3, t4, ik, kk, i, l) reduction(+:sx) reduction(+:sy) 
      for(k = 1; k <= np; k += 1) {
         kk = k_offset+k;
         t1 = 271828183.0;
         t2 = an;
         i = 1;
l99999:         ;

/*      Find starting seed t1 for this kk. */

         if (!(i<=100)) goto _break_7;
         ik = kk/2;
         if (2*ik!=kk)
            t3 = randlc(&t1, t2);
         if (ik==0) goto _break_7;
         t3 = randlc(&t2, t2);
         kk = ik;
         i++;
         goto l99999;
_break_7:         ;
         
         /*      Compute uniform pseudorandom numbers. */
         
         if (0==1)
            timer_start(3);
         vranlc(2*(1<<16), &t1, 1220703125.0, x-1);
         if (0==1)
            timer_stop(3);
         
         /*
         c       Compute Gaussian deviates by acceptance-rejection method and 
         c       tally counts in concentric square annuli.  This loop is not 
         c       vectorizable.
         */
         if (0==1)
            timer_start(2);

         for(i = 0; i <= (1<<16)-1; i += 1) {
            x1 = 2.0*x[2*i]-1.0;
            x2 = 2.0*x[2*i+1]-1.0;
            t1 = x1*x1+x2*x2;
            if (t1<=1.0) {
               t2 = sqrt(-2.0*log(t1)/t1);
               t3 = x1*t2;
               /* Xi */
               t4 = x2*t2;
               /* Yi */
               l = fabs(t3)>fabs(t4)?fabs(t3):fabs(t4);
               qq[l] += 1.0;
               /* counts */
               sx = sx+t3;
               /* sum of Xi */
               sy = sy+t4;
            }
         }
         if (0==1)
            timer_stop(2);
      }
#pragma omp parallel for
      for(i = 0; i <= 10-1; i += 1)
         q[i] += qq[i];
   }
   /* end of parallel region */
   
#pragma omp parallel for reduction(+:gc)
   for(i = 0; i <= 10-1; i += 1)
      gc = gc+q[i];

   timer_stop(1);
   tm = timer_read(1);

   nit = 0;
   if (30==24) {
      if (fabs((sx-(-3.247834652034740e3))/sx)<=1.0e-8&&fabs((sy-(-6.958407078382297e3))/sy)<=1.0e-8)
         verified = 1;
   }
   else if (30==25)
      if (fabs((sx-(-2.863319731645753e3))/sx)<=1.0e-8&&fabs((sy-(-6.320053679109499e3))/sy)<=1.0e-8)
         verified = 1;
   else if (30==28)
      if (fabs((sx-(-4.295875165629892e3))/sx)<=1.0e-8&&fabs((sy-(-1.580732573678431e4))/sy)<=1.0e-8)
         verified = 1;
   else if (30==30)
      if (fabs((sx-4.033815542441498e4)/sx)<=1.0e-8&&fabs((sy-(-2.660669192809235e4))/sy)<=1.0e-8)
         verified = 1;
   else if (30==32)
      if (fabs((sx-4.764367927995374e4)/sx)<=1.0e-8&&fabs((sy-(-8.084072988043731e4))/sy)<=1.0e-8)
         verified = 1;

   Mops = pow(2.0, 30+1)/tm/1000000.0;
   
   
   
   
   
   
   
   printf("EP Benchmark Results: \n""CPU Time = %10.4f\n""N = 2^%5d\n""No. Gaussian Pairs = %15.0f\n""Sums = %25.15e %25.15e\n""Counts:\n", tm, 30, gc, sx, sy);
   for(i = 0; i <= 10-1; i += 1)
      printf("%3d %15.0f\n", i, q[i]);
   
   
   
   
   
   c_print_results("EP", 'B', 30+1, 0, 0, nit, nthreads, tm, Mops, "Random numbers generated", verified, "3.0 structured", "01 Dec 2023", "(none)", "(none)", "-lm", "(none)", "(none)", "(none)", "randdp");

   if (0==1) {
      printf("Total time:     %f", timer_read(1));
      printf("Gaussian pairs: %f", timer_read(2));
      printf("Random numbers: %f", timer_read(3));
   }
}
