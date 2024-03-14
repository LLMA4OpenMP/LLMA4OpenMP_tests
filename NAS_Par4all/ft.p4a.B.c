/*
 * file for ft.c
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



/* global variables */
#include "global.h"

/* function declarations */

static void evolve(dcomplex u0[256][256][512], dcomplex u1[256][256][512], int t, int indexmap[256][256][512], int d[3]);
static void compute_initial_conditions(dcomplex u0[256][256][512], int d[3]);
static void ipow46(double a, int exponent, double *result);
static void setup(void);
static void compute_indexmap(int indexmap[256][256][512], int d[3]);
static void print_timers(void);
static void fft(int dir, dcomplex x1[256][256][512], dcomplex x2[256][256][512]);



static void cffts1(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18]);



static void cffts2(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18]);



static void cffts3(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18]);
static void fft_init(int n);

static void cfftz(int is, int m, int n, dcomplex x[512][18], dcomplex y[512][18]);


static void fftz2(int is, int l, int m, int n, int ny, int ny1, dcomplex u[512], dcomplex x[512][18], dcomplex y[512][18]);
static int ilog2(int n);
static void checksum(int i, dcomplex u1[256][256][512], int d[3]);

static void verify(int d1, int d2, int d3, int nt, boolean *verified, char *class);

/*--------------------------------------------------------------------
c FT benchmark
c-------------------------------------------------------------------*/

int main(int argc, char **argv);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void evolve(dcomplex u0[256][256][512], dcomplex u1[256][256][512], int t, int indexmap[256][256][512], int d[3]);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void compute_initial_conditions(dcomplex u0[256][256][512], int d[3]);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void ipow46(double a, int exponent, double *result);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void setup(void);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void compute_indexmap(int indexmap[256][256][512], int d[3]);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void print_timers(void);



/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void fft(int dir, dcomplex x1[256][256][512], dcomplex x2[256][256][512]);



/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void cffts1(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18]);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void cffts2(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18]);

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void cffts3(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18]);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void fft_init(int n);



/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void cfftz(int is, int m, int n, dcomplex x[512][18], dcomplex y[512][18]);



/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void fftz2(int is, int l, int m, int n, int ny, int ny1, dcomplex u[512], dcomplex x[512][18], dcomplex y[512][18]);



/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static int ilog2(int n);



/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void checksum(int i, dcomplex u1[256][256][512], int d[3]);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void verify(int d1, int d2, int d3, int nt, boolean *verified, char *class);
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
   
   /*c-------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int i, ierr;
   
   /*------------------------------------------------------------------
   c u0, u1, u2 are the main arrays in the problem. 
   c Depending on the decomposition, these arrays will have different 
   c dimensions. To accomodate all possibilities, we allocate them as 
   c one-dimensional arrays and pass them to subroutines for different 
   c views
   c  - u0 contains the initial (transformed) initial condition
   c  - u1 and u2 are working arrays
   c  - indexmap maps i,j,k of u0 to the correct i^2+j^2+k^2 for the
   c    time evolution operator. 
   c-----------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c Large arrays are in common so that they are allocated on the
   c heap rather than the stack. This common block is not
   c referenced directly anywhere else. Padding is to avoid accidental 
   c cache problems, since all array sizes are powers of two.
   c-------------------------------------------------------------------*/
   static dcomplex u0[256][256][512];
   static dcomplex pad1[3];
   static dcomplex u1[256][256][512];
   static dcomplex pad2[3];
   static dcomplex u2[256][256][512];
   static dcomplex pad3[3];
   static int indexmap[256][256][512];

   int iter;
   int nthreads = 1;
   double total_time, mflops;
   boolean verified;
   char class;
   
   /*--------------------------------------------------------------------
   c Run the entire problem once to make sure all data is touched. 
   c This reduces variable startup costs, which is important for such a 
   c short benchmark. The other NPB 2 implementations are similar. 
   c-------------------------------------------------------------------*/
#pragma omp parallel for
   for(i = 0; i <= 6; i += 1)
      timer_clear(i);
   setup();

   compute_indexmap(indexmap, dims[2]);

   compute_initial_conditions(u1, dims[0]);
   fft_init(dims[0][0]);
   
   
   fft(1, u1, u0);
   
   
   
   /*--------------------------------------------------------------------
   c Start over from the beginning. Note that all operations must
   c be timed, in contrast to other benchmarks. 
   c-------------------------------------------------------------------*/
#pragma omp parallel for
   for(i = 0; i <= 6; i += 1)
      timer_clear(i);

   timer_start(0);
   if (0==1)
      timer_start(1);

   compute_indexmap(indexmap, dims[2]);

   compute_initial_conditions(u1, dims[0]);
   fft_init(dims[0][0]);
   
   
   if (0==1)
      timer_stop(1);
   if (0==1)
      timer_start(2);

   fft(1, u1, u0);
   if (0==1)
      timer_stop(2);

   for(iter = 1; iter <= niter; iter += 1) {
      if (0==1)
         timer_start(3);

      evolve(u0, u1, iter, indexmap, dims[0]);

      if (0==1)
         timer_stop(3);
      if (0==1)
         timer_start(2);

      fft(-1, u1, u2);

      if (0==1)
         timer_stop(2);
      if (0==1)
         timer_start(4);
      checksum(iter, u2, dims[0]);

      if (0==1)
         timer_stop(4);
   }

   verify(512, 256, 256, niter, &verified, &class);
   /* end parallel */
   
   timer_stop(0);
   total_time = timer_read(0);

   if (total_time!=0.0)
      mflops = 1.0e-6*(double) 33554432*(14.8157+7.19641*log((double) 33554432)+(5.23518+7.21113*log((double) 33554432))*niter)/total_time;
   else
      mflops = 0.0;
   
   
   
   c_print_results("FT", class, 512, 256, 256, niter, nthreads, total_time, mflops, "          floating point", verified, "3.0 structured", "01 Dec 2023", "(none)", "(none)", "-lm", "(none)", "(none)", "(none)", "randdp");
   if (0==1)
      print_timers();
}
static void evolve(dcomplex u0[256][256][512], dcomplex u1[256][256][512], int t, int indexmap[256][256][512], int d[3])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c evolve u0 -> u1 (t time steps) in fourier space
   c-------------------------------------------------------------------*/
   
   int i, j, k;

#pragma omp parallel for private(j, i)
   for(k = 0; k <= d[2]-1; k += 1)
      for(j = 0; j <= d[1]-1; j += 1)
         for(i = 0; i <= d[0]-1; i += 1)
            u1[k][j][i].real = u0[k][j][i].real*ex[t*indexmap[k][j][i]], u1[k][j][i].imag = u0[k][j][i].imag*ex[t*indexmap[k][j][i]];
}
static void compute_initial_conditions(dcomplex u0[256][256][512], int d[3])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c Fill in array u0 with initial conditions from 
   c random number generator 
   c-------------------------------------------------------------------*/
   
   int k;
   double x0, start, an, dummy;
   static double tmp[512*2*512+1];
   int i, j, t;

   start = 314159265.0;
   /*--------------------------------------------------------------------
   c Jump to the starting element for our first plane.
   c-------------------------------------------------------------------*/
   ipow46(1220703125.0, (zstart[0]-1)*2*512*256+(ystart[0]-1)*2*512, &an);
   dummy = randlc(&start, an);
   ipow46(1220703125.0, 2*512*256, &an);
   
   /*--------------------------------------------------------------------
   c Go through by z planes filling in one square at a time.
   c-------------------------------------------------------------------*/
   for(k = 0; k <= dims[0][2]-1; k += 1) {
      x0 = start;
      vranlc(2*512*dims[0][1], &x0, 1220703125.0, tmp);

      t = 1;
      for(j = 0; j <= dims[0][1]-1; j += 1)
         for(i = 0; i <= 511; i += 1) {
            u0[k][j][i].real = tmp[t++];
            u0[k][j][i].imag = tmp[t++];
         }

      if (k!=dims[0][2])
         dummy = randlc(&start, an);
   }
}
static void ipow46(double a, int exponent, double *result)
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c compute a^exponent mod 2^46
   c-------------------------------------------------------------------*/
   
   double dummy, q, r;
   int n, n2;
   
   /*--------------------------------------------------------------------
   c Use
   c   a^n = a^(n/2)*a^(n/2) if n even else
   c   a^n = a*a^(n-1)       if n odd
   c-------------------------------------------------------------------*/
   *result = 1;
   if (exponent==0)
      return;
   q = a;
   r = 1;
   n = exponent;

   while (n>1) {
      n2 = n/2;
      if (n2*2==n) {
         dummy = randlc(&q, q);
         n = n2;
      }
      else {
         dummy = randlc(&r, q);
         n = n-1;
      }
   }
   dummy = randlc(&r, q);
   *result = r;
}
static void setup(void)
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int ierr, i, j, fstatus;
   
   
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"" - FT Benchmark\n\n");

   niter = 20;

   printf(" Size                : %3dx%3dx%3d\n", 512, 256, 256);
   printf(" Iterations          :     %7d\n", niter);
   
   /* 1004 format(' Number of processes :     ', i7)
    1005 format(' Processor array     :     ', i3, 'x', i3)
    1006 format(' WARNING: compiled for ', i5, ' processes. ',
        >       ' Will not verify. ')*/
   
#pragma omp parallel for
   for(i = 0; i <= 2; i += 1) {
      dims[i][0] = 512;
      dims[i][1] = 256;
      dims[i][2] = 256;
   }
   
   
#pragma omp parallel for
   for(i = 0; i <= 2; i += 1) {
      xstart[i] = 1;
      xend[i] = 512;
      ystart[i] = 1;
      yend[i] = 256;
      zstart[i] = 1;
      zend[i] = 256;
   }
   
   /*--------------------------------------------------------------------
   c Set up info for blocking of ffts and transposes.  This improves
   c performance on cache-based systems. Blocking involves
   c working on a chunk of the problem at a time, taking chunks
   c along the first, second, or third dimension. 
   c
   c - In cffts1 blocking is on 2nd dimension (with fft on 1st dim)
   c - In cffts2/3 blocking is on 1st dimension (with fft on 2nd and 3rd dims)
   
   c Since 1st dim is always in processor, we'll assume it's long enough 
   c (default blocking factor is 16 so min size for 1st dim is 16)
   c The only case we have to worry about is cffts1 in a 2d decomposition. 
   c so the blocking factor should not be larger than the 2nd dimension. 
   c-------------------------------------------------------------------*/
   
   fftblock = 16;
   fftblockpad = 18;

   if (fftblock!=16)
      fftblockpad = fftblock+3;
}
static void compute_indexmap(int indexmap[256][256][512], int d[3])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c compute function from local (i,j,k) to ibar^2+jbar^2+kbar^2 
   c for time evolution exponent. 
   c-------------------------------------------------------------------*/
   
   int i, j, k, ii, ii2, jj, ij2, kk;
   double ap;
   
   /*--------------------------------------------------------------------
   c basically we want to convert the fortran indices 
   c   1 2 3 4 5 6 7 8 
   c to 
   c   0 1 2 3 -4 -3 -2 -1
   c The following magic formula does the trick:
   c mod(i-1+n/2, n) - n/2
   c-------------------------------------------------------------------*/
   
#pragma omp parallel for private(j, ii2, ii, k, jj, ij2, kk)
   for(i = 0; i <= dims[2][0]-1; i += 1) {
      ii = (i+1+xstart[2]-2+512/2)%512-512/2;
      ii2 = ii*ii;
      for(j = 0; j <= dims[2][1]-1; j += 1) {
         jj = (j+1+ystart[2]-2+256/2)%256-256/2;
         ij2 = jj*jj+ii2;
         for(k = 0; k <= dims[2][2]-1; k += 1) {
            kk = (k+1+zstart[2]-2+256/2)%256-256/2;
            indexmap[k][j][i] = kk*kk+ij2;
         }
      }
   }
   
   /*--------------------------------------------------------------------
   c compute array of exponentials for time evolution. 
   c-------------------------------------------------------------------*/
   ap = -4.0*1.0e-6*3.141592653589793238*3.141592653589793238;

   ex[0] = 1.0;
   ex[1] = exp(ap);
   for(i = 2; i <= 20*(512*512/4+256*256/4+256*256/4); i += 1)
      ex[i] = ex[i-1]*ex[1];
}
static void print_timers(void)
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int i;
   
   
   
   
   
   
   char *tstrings[] = {"          total ", "          setup ", "            fft ", "         evolve ", "       checksum ", "         fftlow ", "        fftcopy "};

   for(i = 0; i <= 6; i += 1)
      if (timer_read(i)!=0.0)
         printf("timer %2d(%16s( :%10.6f\n", i, tstrings[i], timer_read(i));
}
static void fft(int dir, dcomplex x1[256][256][512], dcomplex x2[256][256][512])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   dcomplex y0[512][18];
   dcomplex y1[512][18];
   
   /*--------------------------------------------------------------------
   c note: args x1, x2 must be different arrays
   c note: args for cfftsx are (direction, layout, xin, xout, scratch)
   c       xin/xout may be the same and it can be somewhat faster
   c       if they are
   c-------------------------------------------------------------------*/
   
   if (dir==1) {
      cffts1(1, dims[0], x1, x1, y0, y1);
      /* x1 -> x1 */
      cffts2(1, dims[1], x1, x1, y0, y1);
      /* x1 -> x1 */
      cffts3(1, dims[2], x1, x2, y0, y1);
   }
   else {
      cffts3(-1, dims[2], x1, x1, y0, y1);
      /* x1 -> x1 */
      cffts2(-1, dims[1], x1, x1, y0, y1);
      /* x1 -> x1 */
      cffts1(-1, dims[0], x1, x2, y0, y1);
   }
}
static void cffts1(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int logd[3];
   int i, j, k, jj;

#pragma omp parallel for
   for(i = 0; i <= 2; i += 1)
      logd[i] = ilog2(d[i]);
   {
      dcomplex y0[512][18];
      dcomplex y1[512][18];

      for(k = 0; k <= d[2]-1; k += 1) {
         jj = 0;
         while (jj<=d[1]-fftblock) {
            /*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
#pragma omp parallel for private(i)
            for(j = 0; j <= fftblock-1; j += 1)
               for(i = 0; i <= d[0]-1; i += 1) {

                  y0[i][j].real = x[k][j+jj][i].real;

                  y0[i][j].imag = x[k][j+jj][i].imag;
               }
            /*          if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTCOPY); */
            
            /*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTLOW); */
            
            cfftz(is, logd[0], d[0], y0, y1);
            
            /*          if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTLOW); */
            /*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
#pragma omp parallel for private(i)
            for(j = 0; j <= fftblock-1; j += 1)
               for(i = 0; i <= d[0]-1; i += 1) {
                  xout[k][j+jj][i].real = y0[i][j].real;
                  xout[k][j+jj][i].imag = y0[i][j].imag;
               }
            jj += fftblock;
         }
      }
   }
}
static void cffts2(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int logd[3];
   int i, j, k, ii;

#pragma omp parallel for
   for(i = 0; i <= 2; i += 1)
      logd[i] = ilog2(d[i]);
   {
      dcomplex y0[512][18];
      dcomplex y1[512][18];

      for(k = 0; k <= d[2]-1; k += 1) {
         ii = 0;
         while (ii<=d[0]-fftblock) {
            /*	    if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
#pragma omp parallel for private(i)
            for(j = 0; j <= d[1]-1; j += 1)
               for(i = 0; i <= fftblock-1; i += 1) {
                  y0[j][i].real = x[k][j][i+ii].real;
                  y0[j][i].imag = x[k][j][i+ii].imag;
               }
            /*	    if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTCOPY); */
            /*	    if (TIMERS_ENABLED == TRUE) timer_start(T_FFTLOW); */
            
            cfftz(is, logd[1], d[1], y0, y1);
            
            /*          if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTLOW); */
            /*          if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
#pragma omp parallel for private(i)
            for(j = 0; j <= d[1]-1; j += 1)
               for(i = 0; i <= fftblock-1; i += 1) {
                  xout[k][j][i+ii].real = y0[j][i].real;
                  xout[k][j][i+ii].imag = y0[j][i].imag;
               }
            ii += fftblock;
         }
      }
   }
}
static void cffts3(int is, int d[3], dcomplex x[256][256][512], dcomplex xout[256][256][512], dcomplex y0[512][18], dcomplex y1[512][18])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int logd[3];
   int i, j, k, ii;

#pragma omp parallel for
   for(i = 0; i <= 2; i += 1)
      logd[i] = ilog2(d[i]);
   {
      dcomplex y0[512][18];
      dcomplex y1[512][18];

      for(j = 0; j <= d[1]-1; j += 1) {
         ii = 0;
         while (ii<=d[0]-fftblock) {
            /*	    if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
#pragma omp parallel for private(i)
            for(k = 0; k <= d[2]-1; k += 1)
               for(i = 0; i <= fftblock-1; i += 1) {
                  y0[k][i].real = x[k][j][i+ii].real;
                  y0[k][i].imag = x[k][j][i+ii].imag;
               }
            
            /*           if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTCOPY); */
            /*           if (TIMERS_ENABLED == TRUE) timer_start(T_FFTLOW); */
            
            cfftz(is, logd[2], d[2], y0, y1);
            /*           if (TIMERS_ENABLED == TRUE) timer_stop(T_FFTLOW); */
            /*           if (TIMERS_ENABLED == TRUE) timer_start(T_FFTCOPY); */
#pragma omp parallel for private(i)
            for(k = 0; k <= d[2]-1; k += 1)
               for(i = 0; i <= fftblock-1; i += 1) {
                  xout[k][j][i+ii].real = y0[k][i].real;
                  xout[k][j][i+ii].imag = y0[k][i].imag;
               }
            ii += fftblock;
         }
      }
   }
}
static void fft_init(int n)
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c compute the roots-of-unity array that will be used for subsequent FFTs. 
   c-------------------------------------------------------------------*/
   
   int m, nu, ku, i, j, ln;
   double t, ti;
   
   
   /*--------------------------------------------------------------------
   c   Initialize the U array with sines and cosines in a manner that permits
   c   stride one access at each FFT iteration.
   c-------------------------------------------------------------------*/
   nu = n;
   m = ilog2(n);
   u[0].real = (double) m;
   u[0].imag = 0.0;
   ku = 1;
   ln = 1;

   for(j = 1; j <= m; j += 1) {
      t = 3.141592653589793238/ln;

#pragma omp parallel for private(ti)
      for(i = 0; i <= ln-1; i += 1) {
         ti = i*t;
         u[i+ku].real = cos(ti);
         u[i+ku].imag = sin(ti);
      }

      ku = ku+ln;
      ln = 2*ln;
   }
}
static void cfftz(int is, int m, int n, dcomplex x[512][18], dcomplex y[512][18])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c   Computes NY N-point complex-to-complex FFTs of X using an algorithm due
   c   to Swarztrauber.  X is both the input and the output array, while Y is a 
   c   scratch array.  It is assumed that N = 2^M.  Before calling CFFTZ to 
   c   perform FFTs, the array U must be initialized by calling CFFTZ with IS 
   c   set to 0 and M set to MX, where MX is the maximum value of M for any 
   c   subsequent call.
   c-------------------------------------------------------------------*/
   
   int i, j, l, mx;
   
   /*--------------------------------------------------------------------
   c   Check if input parameters are invalid.
   c-------------------------------------------------------------------*/
   mx = (int) u[0].real;
   if (is!=1&&is!=-1||m<1||m>mx) {
      
      
      printf("CFFTZ: Either U has not been initialized, or else\n""one of the input parameters is invalid%5d%5d%5d\n", is, m, mx);
      exit(1);
   }
   l = 1;
l99999:   ;

/*--------------------------------------------------------------------
c   Perform one variant of the Stockham FFT.
c-------------------------------------------------------------------*/
   if (!(l<=m)) goto _break_1;
   fftz2(is, l, m, n, fftblock, fftblockpad, u, x, y);
   if (l==m) goto _break_1;
   fftz2(is, l+1, m, n, fftblock, fftblockpad, u, y, x);
   l += 2;
   goto l99999;
_break_1:   ;
   
   /*--------------------------------------------------------------------
   c   Copy Y to X.
   c-------------------------------------------------------------------*/
   if (m%2==1)
#pragma omp parallel for private(i)
      for(j = 0; j <= n-1; j += 1)
         for(i = 0; i <= fftblock-1; i += 1) {
            x[j][i].real = y[j][i].real;
            x[j][i].imag = y[j][i].imag;
         }
}
static void fftz2(int is, int l, int m, int n, int ny, int ny1, dcomplex u[512], dcomplex x[512][18], dcomplex y[512][18])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c   Performs the L-th iteration of the second variant of the Stockham FFT.
   c-------------------------------------------------------------------*/
   
   int k, n1, li, lj, lk, ku, i, j, i11, i12, i21, i22;
   dcomplex u1, x11, x21;
   
   /*--------------------------------------------------------------------
   c   Set initial parameters.
   c-------------------------------------------------------------------*/
   
   n1 = n/2;
   if (l-1==0)
      lk = 1;
   else
      lk = 2<<l-1-1;
   if (m-l==0)
      li = 1;
   else
      li = 2<<m-l-1;
   lj = 2*lk;
   ku = li;

   for(i = 0; i <= li-1; i += 1) {

      i11 = i*lk;
      i12 = i11+n1;
      i21 = i*lj;
      i22 = i21+lk;
      if (is>=1) {
         u1.real = u[ku+i].real;
         u1.imag = u[ku+i].imag;
      }
      else {
         u1.real = u[ku+i].real;
         u1.imag = -u[ku+i].imag;
      }
      
      /*--------------------------------------------------------------------
      c   This loop is vectorizable.
      c-------------------------------------------------------------------*/
      for(k = 0; k <= lk-1; k += 1)
#pragma omp parallel for
         for(j = 0; j <= ny-1; j += 1) {
            double x11real, x11imag;
            double x21real, x21imag;
            x11real = x[i11+k][j].real;
            x11imag = x[i11+k][j].imag;
            x21real = x[i12+k][j].real;
            x21imag = x[i12+k][j].imag;
            y[i21+k][j].real = x11real+x21real;
            y[i21+k][j].imag = x11imag+x21imag;
            y[i22+k][j].real = u1.real*(x11real-x21real)-u1.imag*(x11imag-x21imag);
            y[i22+k][j].imag = u1.real*(x11imag-x21imag)+u1.imag*(x11real-x21real);
         }
   }
}
static int ilog2(int n)
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int nn, lg;

   if (n==1)
      return 0;
   lg = 1;
   nn = 2;
   while (nn<n) {
      nn = nn<<1;
      lg++;
   }

   return lg;
}
static void checksum(int i, dcomplex u1[256][256][512], int d[3])
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int j, q, r, s, ierr;
   dcomplex chk, allchk;

   chk.real = 0.0;
   chk.imag = 0.0;
   
   
   for(j = 1; j <= 1024; j += 1) {
      q = j%512+1;
      if (q>=xstart[0]&&q<=xend[0]) {
         r = 3*j%256+1;
         if (r>=ystart[0]&&r<=yend[0]) {
            s = 5*j%256+1;
            if (s>=zstart[0]&&s<=zend[0])
               chk.real = chk.real+u1[s-zstart[0]][r-ystart[0]][q-xstart[0]].real, chk.imag = chk.imag+u1[s-zstart[0]][r-ystart[0]][q-xstart[0]].imag;
         }
      }
   }
   sums[i].real += chk.real;
   sums[i].imag += chk.imag;
   /* complex % real */
   sums[i].real = sums[i].real/((double) 33554432);
   sums[i].imag = sums[i].imag/((double) 33554432);
   
   
   printf("T = %5d     Checksum = %22.12e %22.12e\n", i, sums[i].real, sums[i].imag);
}
static void verify(int d1, int d2, int d3, int nt, boolean *verified, char *class)
{
   
   /*--------------------------------------------------------------------
   c-------------------------------------------------------------------*/
   
   int ierr, size, i;
   double err, epsilon;
   
   /*--------------------------------------------------------------------
   c   Sample size reference checksums
   c-------------------------------------------------------------------*/
   
   /*--------------------------------------------------------------------
   c   Class S size reference checksums
   c-------------------------------------------------------------------*/
   
   
   
   
   
   
   double vdata_real_s[6+1] = {0.0, 5.546087004964e+02, 5.546385409189e+02, 5.546148406171e+02, 5.545423607415e+02, 5.544255039624e+02, 5.542683411902e+02};
   
   
   
   
   
   
   double vdata_imag_s[6+1] = {0.0, 4.845363331978e+02, 4.865304269511e+02, 4.883910722336e+02, 4.901273169046e+02, 4.917475857993e+02, 4.932597244941e+02};
   /*--------------------------------------------------------------------
   c   Class W size reference checksums
   c-------------------------------------------------------------------*/
   
   
   
   
   
   
   double vdata_real_w[6+1] = {0.0, 5.673612178944e+02, 5.631436885271e+02, 5.594024089970e+02, 5.560698047020e+02, 5.530898991250e+02, 5.504159734538e+02};
   
   
   
   
   
   
   double vdata_imag_w[6+1] = {0.0, 5.293246849175e+02, 5.282149986629e+02, 5.270996558037e+02, 5.260027904925e+02, 5.249400845633e+02, 5.239212247086e+02};
   /*--------------------------------------------------------------------
   c   Class A size reference checksums
   c-------------------------------------------------------------------*/
   
   
   
   
   
   
   double vdata_real_a[6+1] = {0.0, 5.046735008193e+02, 5.059412319734e+02, 5.069376896287e+02, 5.077892868474e+02, 5.085233095391e+02, 5.091487099959e+02};
   
   
   
   
   
   
   double vdata_imag_a[6+1] = {0.0, 5.114047905510e+02, 5.098809666433e+02, 5.098144042213e+02, 5.101336130759e+02, 5.104914655194e+02, 5.107917842803e+02};
   /*--------------------------------------------------------------------
   c   Class B size reference checksums
   c-------------------------------------------------------------------*/
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   double vdata_real_b[20+1] = {0.0, 5.177643571579e+02, 5.154521291263e+02, 5.146409228649e+02, 5.142378756213e+02, 5.139626667737e+02, 5.137423460082e+02, 5.135547056878e+02, 5.133910925466e+02, 5.132470705390e+02, 5.131197729984e+02, 5.130070319283e+02, 5.129070537032e+02, 5.128182883502e+02, 5.127393733383e+02, 5.126691062020e+02, 5.126064276004e+02, 5.125504076570e+02, 5.125002331720e+02, 5.124551951846e+02, 5.124146770029e+02};
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   double vdata_imag_b[20+1] = {0.0, 5.077803458597e+02, 5.088249431599e+02, 5.096208912659e+02, 5.101023387619e+02, 5.103976610617e+02, 5.105948019802e+02, 5.107404165783e+02, 5.108576573661e+02, 5.109577278523e+02, 5.110460304483e+02, 5.111252433800e+02, 5.111968077718e+02, 5.112616233064e+02, 5.113203605551e+02, 5.113735928093e+02, 5.114218460548e+02, 5.114656139760e+02, 5.115053595966e+02, 5.115415130407e+02, 5.115744692211e+02};
   /*--------------------------------------------------------------------
   c   Class C size reference checksums
   c-------------------------------------------------------------------*/
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   double vdata_real_c[20+1] = {0.0, 5.195078707457e+02, 5.155422171134e+02, 5.144678022222e+02, 5.140150594328e+02, 5.137550426810e+02, 5.135811056728e+02, 5.134569343165e+02, 5.133651975661e+02, 5.132955192805e+02, 5.132410471738e+02, 5.131971141679e+02, 5.131605205716e+02, 5.131290734194e+02, 5.131012720314e+02, 5.130760908195e+02, 5.130528295923e+02, 5.130310107773e+02, 5.130103090133e+02, 5.129905029333e+02, 5.129714421109e+02};
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   double vdata_imag_c[20+1] = {0.0, 5.149019699238e+02, 5.127578201997e+02, 5.122251847514e+02, 5.121090289018e+02, 5.121143685824e+02, 5.121496764568e+02, 5.121870921893e+02, 5.122193250322e+02, 5.122454735794e+02, 5.122663649603e+02, 5.122830879827e+02, 5.122965869718e+02, 5.123075927445e+02, 5.123166486553e+02, 5.123241541685e+02, 5.123304037599e+02, 5.123356167976e+02, 5.123399592211e+02, 5.123435588985e+02, 5.123465164008e+02};

   epsilon = 1.0e-12;
   *verified = 1;
   *class = 'U';
   
   
   
   
   if (d1==64&&d2==64&&d3==64&&nt==6) {
      *class = 'S';
      i = 1;
l99997:      ;
      if (!(i<=nt)) goto _break_1;
      err = (sums[i].real-vdata_real_s[i])/vdata_real_s[i];
      if (fabs(err)>epsilon) goto l99999;
      err = (sums[i].imag-vdata_imag_s[i])/vdata_imag_s[i];
      if (fabs(err)>epsilon) goto l99998;
      i++;
      goto l99997;
l99998:      ;
      *verified = 0;
      goto _break_1;
l99999:      ;
      *verified = 0;
_break_1:      ;
   }
   else
      
      
      if (d1==128&&d2==128&&d3==32&&nt==6) {
         *class = 'W';
         i = 1;
l99994:         ;
         if (!(i<=nt)) goto _break_2;
         err = (sums[i].real-vdata_real_w[i])/vdata_real_w[i];
         if (fabs(err)>epsilon) goto l99996;
         err = (sums[i].imag-vdata_imag_w[i])/vdata_imag_w[i];
         if (fabs(err)>epsilon) goto l99995;
         i++;
         goto l99994;
l99995:         ;
         *verified = 0;
         goto _break_2;
l99996:         ;
         *verified = 0;
_break_2:         ;
      }
      else
         
         
         if (d1==256&&d2==256&&d3==128&&nt==6) {
            *class = 'A';
            i = 1;
l99991:            ;
            if (!(i<=nt)) goto _break_3;
            err = (sums[i].real-vdata_real_a[i])/vdata_real_a[i];
            if (fabs(err)>epsilon) goto l99993;
            err = (sums[i].imag-vdata_imag_a[i])/vdata_imag_a[i];
            if (fabs(err)>epsilon) goto l99992;
            i++;
            goto l99991;
l99992:            ;
            *verified = 0;
            goto _break_3;
l99993:            ;
            *verified = 0;
_break_3:            ;
         }
         else
            
            
            if (d1==512&&d2==256&&d3==256&&nt==20) {
               *class = 'B';
               i = 1;
l99988:               ;
               if (!(i<=nt)) goto _break_4;
               err = (sums[i].real-vdata_real_b[i])/vdata_real_b[i];
               if (fabs(err)>epsilon) goto l99990;
               err = (sums[i].imag-vdata_imag_b[i])/vdata_imag_b[i];
               if (fabs(err)>epsilon) goto l99989;
               i++;
               goto l99988;
l99989:               ;
               *verified = 0;
               goto _break_4;
l99990:               ;
               *verified = 0;
_break_4:               ;
            }
            else
               
               
               if (d1==512&&d2==512&&d3==512&&nt==20) {
                  *class = 'C';
                  i = 1;
l99985:                  ;
                  if (!(i<=nt)) goto _break_5;
                  err = (sums[i].real-vdata_real_c[i])/vdata_real_c[i];
                  if (fabs(err)>epsilon) goto l99987;
                  err = (sums[i].imag-vdata_imag_c[i])/vdata_imag_c[i];
                  if (fabs(err)>epsilon) goto l99986;
                  i++;
                  goto l99985;
l99986:                  ;
                  *verified = 0;
                  goto _break_5;
l99987:                  ;
                  *verified = 0;
_break_5:                  ;
               }

   if (*class!='U')
      printf("Result verification successful\n");
   else
      printf("Result verification failed\n");
   printf("class = %1c\n", *class);
}
