/*
 * file for cg.c
 */
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>




typedef int boolean;
typedef struct {
   double real;
   double imag;
} dcomplex;
double randlc(double *x, double a);
void vranlc(int n, double *x_seed, double a, double y[]);
void timer_clear(int n);
void timer_start(int n);
void timer_stop(int n);
double timer_read(int n);






void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand);



#include <sys/time.h>

void wtime_(double *t);


/*
*/
/*c---------------------------------------------------------------------
c---------------------------------------------------------------------*/

double randlc(double *x, double a);

/*c---------------------------------------------------------------------
c---------------------------------------------------------------------*/

void vranlc(int n, double *x_seed, double a, double y[]);


#define wtime wtime_
#include <stdlib.h>
/*  Prototype  */
void wtime_(double *t);
/*****************************************************************/
/******         E  L  A  P  S  E  D  _  T  I  M  E          ******/
/*****************************************************************/

double elapsed_time();

double start[64];
double elapsed[64];
/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/

void timer_clear(int n);

/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/

void timer_start(int n);

/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/

void timer_stop(int n);

/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/

double timer_read(int n);


/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/
#include <stdlib.h>
#include <stdio.h>

void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand);

//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
#include "npbparams.h"






static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;
static int colidx[14000*(11+1)*(11+1)+14000*(11+2)+1];
static int rowstr[14000+1+1];
static int iv[2*14000+1+1];
static int arow[14000*(11+1)*(11+1)+14000*(11+2)+1];
static int acol[14000*(11+1)*(11+1)+14000*(11+2)+1];


static double v[14000+1+1];
static double aelt[14000*(11+1)*(11+1)+14000*(11+2)+1];
static double a[14000*(11+1)*(11+1)+14000*(11+2)+1];
static double x[14000+2+1];
static double z[14000+2+1];
static double p[14000+2+1];
static double q[14000+2+1];
static double r[14000+2+1];


static double amult;
static double tran;




static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm);



static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift);



static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza);

static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);



int main(int argc, char **argv);



static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm);



static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift);



static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza);



static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]);



static int icnvrt(double x, int ipwr2);



static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);
void wtime_(double *t)
{
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, (void *) 0);
   if (sec<0)
      sec = tv.tv_sec;
   *t = tv.tv_sec-sec+1.0e-6*tv.tv_usec;
}
double randlc(double *x, double a)
{
   /*c---------------------------------------------------------------------
   c---------------------------------------------------------------------*/
   /*c---------------------------------------------------------------------
   c
   c   This routine returns a uniform pseudorandom double precision number in the
   c   range (0, 1) by using the linear congruential generator
   c
   c   x_{k+1} = a x_k  (mod 2^46)
   c
   c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
   c   before repeating.  The argument A is the same as 'a' in the above formula,
   c   and X is the same as x_0.  A and X must be odd double precision integers
   c   in the range (1, 2^46).  The returned value RANDLC is normalized to be
   c   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
   c   the new seed x_1, so that subsequent calls to RANDLC using the same
   c   arguments will generate a continuous sequence.
   c
   c   This routine should produce the same results on any computer with at least
   c   48 mantissa bits in double precision floating point data.  On 64 bit
   c   systems, double precision should be disabled.
   c
   c   David H. Bailey     October 26, 1990
   c
   c---------------------------------------------------------------------*/
   double t1;
   double t2;
   double t3;
   double t4;
   double a1;
   double a2;
   double x1;
   double x2;
   double z;
   /*c---------------------------------------------------------------------
   c   Break A into two parts such that A = 2^23 * A1 + A2.
   c---------------------------------------------------------------------*/
   t1 = 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*a;
   a1 = (int) t1;
   a2 = a-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*a1;
   /*c---------------------------------------------------------------------
   c   Break X into two parts such that X = 2^23 * X1 + X2, compute
   c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
   c   X = 2^23 * Z + A2 * X2  (mod 2^46).
   c---------------------------------------------------------------------*/
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
   /*c---------------------------------------------------------------------
   c---------------------------------------------------------------------*/
   /*c---------------------------------------------------------------------
   c
   c   This routine generates N uniform pseudorandom double precision numbers in
   c   the range (0, 1) by using the linear congruential generator
   c
   c   x_{k+1} = a x_k  (mod 2^46)
   c
   c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
   c   before repeating.  The argument A is the same as 'a' in the above formula,
   c   and X is the same as x_0.  A and X must be odd double precision integers
   c   in the range (1, 2^46).  The N results are placed in Y and are normalized
   c   to be between 0 and 1.  X is updated to contain the new seed, so that
   c   subsequent calls to VRANLC using the same arguments will generate a
   c   continuous sequence.  If N is zero, only initialization is performed, and
   c   the variables X, A and Y are ignored.
   c
   c   This routine is the standard version designed for scalar or RISC systems.
   c   However, it should produce the same results on any single processor
   c   computer with at least 48 mantissa bits in double precision floating point
   c   data.  On 64 bit systems, double precision should be disabled.
   c
   c---------------------------------------------------------------------*/
   int i;
   double x;
   double t1;
   double t2;
   double t3;
   double t4;
   double a1;
   double a2;
   double x1;
   double x2;
   double z;
   /*c---------------------------------------------------------------------
   c   Break A into two parts such that A = 2^23 * A1 + A2.
   c---------------------------------------------------------------------*/
   t1 = 0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*a;
   a1 = (int) t1;
   a2 = a-2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*a1;
   x = *x_seed;
   /*c---------------------------------------------------------------------
   c   Generate N results.   This loop is not vectorizable.
   c---------------------------------------------------------------------*/
   for(i = 1; i <= n; i += 1) {
      /*c---------------------------------------------------------------------
      c   Break X into two parts such that X = 2^23 * X1 + X2, compute
      c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
      c   X = 2^23 * Z + A2 * X2  (mod 2^46).
      c---------------------------------------------------------------------*/
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
double elapsed_time()
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
   double t;
   double now;
   now = elapsed_time();
   t = now-start[n];
   elapsed[n] += t;
}
double timer_read(int n)
{
   return elapsed[n];
}
void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand)
{
   char *evalue = "1000";
   printf("\n\n %s Benchmark Completed\n", name);
   printf(" Class           =                        %c\n", class);
   if (n2==0&&n3==0)
      /* as in IS */
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
   int i, j, k, it;
   int nthreads = 1;
   double zeta;
   double rnorm;
   double norm_temp11;
   double norm_temp12;
   double t, mflops;
   char class;
   boolean verified;
   double zeta_verify_value, epsilon;

   firstrow = 1;
   lastrow = 14000;
   firstcol = 1;
   lastcol = 14000;

   if (14000==1400&&11==7&&15==15&&20.0==10.0) {
      class = 'S';
      zeta_verify_value = 8.5971775078648;
   }
   else if (14000==7000&&11==8&&15==15&&20.0==12.0) {
      class = 'W';
      zeta_verify_value = 10.362595087124;
   }
   else if (14000==14000&&11==11&&15==15&&20.0==20.0) {
      class = 'A';
      zeta_verify_value = 17.130235054029;
   }
   else if (14000==75000&&11==13&&15==75&&20.0==60.0) {
      class = 'B';
      zeta_verify_value = 22.712745482631;
   }
   else if (14000==150000&&11==15&&15==75&&20.0==110.0) {
      class = 'C';
      zeta_verify_value = 28.973605592845;
   }
   else
      class = 'U';
   
   
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"" - CG Benchmark\n");
   printf(" Size: %10d\n", 14000);
   printf(" Iterations: %5d\n", 15);

   naa = 14000;
   nzz = 14000*(11+1)*(11+1)+14000*(11+2);
   
   
   tran = 314159265.0;
   amult = 1220703125.0;
   zeta = randlc(&tran, amult);
   
   
   
   
   makea(naa, nzz, a, colidx, rowstr, 11, firstrow, lastrow, firstcol, lastcol, 1.0e-1, arow, acol, aelt, v, iv, 20.0);
   for(j = 1; j <= lastrow-firstrow+1; j += 1)
#pragma omp parallel for
      for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
         colidx[k] = colidx[k]-firstcol+1;
   
   
#pragma omp parallel for
   for(i = 1; i <= 14000+1; i += 1)
      x[i] = 1.0;
#pragma omp parallel for
   for(j = 1; j <= lastcol-firstcol+1; j += 1) {
      q[j] = 0.0;
      z[j] = 0.0;
      r[j] = 0.0;
      p[j] = 0.0;
   }
   zeta = 0.0;
   
   
   
   for(it = 1; it <= 1; it += 1) {
      
      
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      
      
      norm_temp11 = 0.0;
      norm_temp12 = 0.0;
#pragma omp parallel for reduction(+:norm_temp12) reduction(+:norm_temp11)
      for(j = 1; j <= lastcol-firstcol+1; j += 1) {
         norm_temp11 = norm_temp11+x[j]*z[j];
         norm_temp12 = norm_temp12+z[j]*z[j];
      }
      norm_temp12 = 1.0/sqrt(norm_temp12);
      
      
#pragma omp parallel for
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         x[j] = norm_temp12*z[j];
   }
   
   
#pragma omp parallel for
   for(i = 1; i <= 14000+1; i += 1)
      x[i] = 1.0;
   zeta = 0.0;
   
   
   timer_clear(1);
   timer_start(1);
   
   
   
   for(it = 1; it <= 15; it += 1) {
      
      
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      
      
      norm_temp11 = 0.0;
      norm_temp12 = 0.0;

#pragma omp parallel for reduction(+:norm_temp12) reduction(+:norm_temp11)
      for(j = 1; j <= lastcol-firstcol+1; j += 1) {
         norm_temp11 = norm_temp11+x[j]*z[j];
         norm_temp12 = norm_temp12+z[j]*z[j];
      }

      norm_temp12 = 1.0/sqrt(norm_temp12);

      zeta = 20.0+1.0/norm_temp11;

      if (it==1)
         printf("   iteration           ||r||                 zeta\n");
      printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);
      
      
#pragma omp parallel for
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         x[j] = norm_temp12*z[j];
   }

   timer_stop(1);
   
   
   
   t = timer_read(1);

   printf(" Benchmark completed\n");

   epsilon = 1.0e-10;
   if (class!='U')
      if (fabs(zeta-zeta_verify_value)<=epsilon) {
         verified = 1;
         printf(" VERIFICATION SUCCESSFUL\n");
         printf(" Zeta is    %20.12e\n", zeta);
         printf(" Error is   %20.12e\n", zeta-zeta_verify_value);
      }
      else {
         verified = 0;
         printf(" VERIFICATION FAILED\n");
         printf(" Zeta                %20.12e\n", zeta);
         printf(" The correct zeta is %20.12e\n", zeta_verify_value);
      }
   else {
      verified = 0;
      printf(" Problem size unknown\n");
      printf(" NO VERIFICATION PERFORMED\n");
   }

   if (t!=0.0)
      mflops = 2.0*15*14000*(3.0+11*(11+1)+25.0*(5.0+11*(11+1))+3.0)/t/1000000.0;
   else
      mflops = 0.0;
   
   
   
   
   c_print_results("CG", class, 14000, 0, 0, 15, nthreads, t, mflops, "          floating point", verified, "3.0 structured", "26 Nov 2023", "gcc-12", "gcc-12", "-fopenmp -lm", "-I../common -fopenmp", "-lm", "(none)", "randdp");
}
static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm)
{
   static int callcount = 0;
   double d, sum, rho, rho0, alpha, beta;
   int i, j, k;
   int cgit, cgitmax = 25;

   rho = 0.0;
#pragma omp parallel for
   for(j = 1; j <= naa+1; j += 1) {
      q[j] = 0.0;
      z[j] = 0.0;
      r[j] = x[j];
      p[j] = r[j];
   }
   
   
#pragma omp parallel for reduction(+:rho)
   for(j = 1; j <= lastcol-firstcol+1; j += 1)
      rho = rho+r[j]*r[j];

   for(cgit = 1; cgit <= cgitmax; cgit += 1) {
      rho0 = rho;
      d = 0.0;
      rho = 0.0;
      
      
      
      
#pragma omp parallel for private(sum, k)
      for(j = 1; j <= lastrow-firstrow+1; j += 1) {
         sum = 0.0;
         for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
            sum = sum+a[k]*p[colidx[k]];
         q[j] = sum;
      }
      
      
      
      
      
      
      
#pragma omp parallel for reduction(+:d)
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         d = d+p[j]*q[j];

      alpha = rho0/d;
      
      
      
      
      
      for(j = 1; j <= lastcol-firstcol+1; j += 1) {
         z[j] = z[j]+alpha*p[j];
         r[j] = r[j]-alpha*q[j];
         
         
         
         rho = rho+r[j]*r[j];
      }
      
      
      beta = rho/rho0;
      
      
#pragma omp parallel for
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         p[j] = r[j]+beta*p[j];
      callcount++;
   }
   
   
   sum = 0.0;
#pragma omp parallel for private(k, d)
   for(j = 1; j <= lastrow-firstrow+1; j += 1) {
      d = 0.0;
      for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
         d = d+a[k]*z[colidx[k]];
      r[j] = d;
   }
   
   
#pragma omp parallel for private(d) reduction(+:sum)
   for(j = 1; j <= lastcol-firstcol+1; j += 1) {
      d = x[j]-r[j];
      sum = sum+d*d;
   }
   *rnorm = sqrt(sum);
}
static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift)
{
   int i, nnza, iouter, ivelt, ivelt1, irow, nzv;
   
   
   
   double size, ratio, scale;
   int jcol;

   size = 1.0;
   ratio = pow(rcond, 1.0/((double) n));
   nnza = 0;
   
   
#pragma omp parallel for
   for(i = 1; i <= n; i += 1)
      colidx[n+i] = 0;
   for(iouter = 1; iouter <= n; iouter += 1) {
      nzv = nonzer;
      sprnvc(n, nzv, v, iv, &colidx[0], &colidx[n]);
      vecset(n, v, iv, &nzv, iouter, 0.5);
      for(ivelt = 1; ivelt <= nzv; ivelt += 1) {
         jcol = iv[ivelt];
         if (jcol>=firstcol&&jcol<=lastcol) {
            scale = size*v[ivelt];
            for(ivelt1 = 1; ivelt1 <= nzv; ivelt1 += 1) {
               irow = iv[ivelt1];
               if (irow>=firstrow&&irow<=lastrow) {
                  nnza = nnza+1;
                  if (nnza>nz) {

                     printf("Space for matrix elements exceeded in"" makea\n");
                     printf("nnza, nzmax = %d, %d\n", nnza, nz);
                     printf("iouter = %d\n", iouter);
                     exit(1);
                  }
                  acol[nnza] = jcol;
                  arow[nnza] = irow;
                  aelt[nnza] = v[ivelt1]*scale;
               }
            }
         }
      }
      size = size*ratio;
   }
   
   
   for(i = firstrow; i <= lastrow; i += 1)
      if (i>=firstcol&&i<=lastcol) {
         iouter = n+i;
         nnza = nnza+1;
         if (nnza>nz) {
            printf("Space for matrix elements exceeded in makea\n");
            printf("nnza, nzmax = %d, %d\n", nnza, nz);
            printf("iouter = %d\n", iouter);
            exit(1);
         }
         acol[nnza] = i;
         arow[nnza] = i;
         aelt[nnza] = rcond-shift;
      }
   sparse(a, colidx, rowstr, n, arow, acol, aelt, firstrow, lastrow, v, &iv[0], &iv[n], nnza);
}
static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza)
{
   int nrows;
   int i, j, jajp1, nza, k, nzrow;
   double xi;
   
   
   nrows = lastrow-firstrow+1;
   
   
#pragma omp parallel for
   for(j = 1; j <= n; j += 1) {
      rowstr[j] = 0;
      mark[j] = 0;
      x[j] = 0.0;
      mark[j] = 0;
   }
   rowstr[n+1] = 0;

   for(nza = 1; nza <= nnza; nza += 1) {
      j = arow[nza]-firstrow+1+1;
      rowstr[j] = rowstr[j]+1;
   }

   rowstr[1] = 1;
   for(j = 2; j <= nrows+1; j += 1)
      rowstr[j] = rowstr[j]+rowstr[j-1];
   
   
   
   
   for(j = 0; j <= nrows-1; j += 1)
#pragma omp parallel for
      for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
         a[k] = 0.0;

   for(nza = 1; nza <= nnza; nza += 1) {
      j = arow[nza]-firstrow+1;
      k = rowstr[j];
      a[k] = aelt[nza];
      colidx[k] = acol[nza];
      rowstr[j] = rowstr[j]+1;
   }
   
   
   for(j = nrows; j >= 1; j += -1)
      rowstr[j+1] = rowstr[j];
   rowstr[1] = 1;
   
   
   nza = 0;

   jajp1 = rowstr[1];
   for(j = 1; j <= nrows; j += 1) {
      nzrow = 0;
      
      
      for(k = jajp1; k <= rowstr[j+1]-1; k += 1) {
         i = colidx[k];
         x[i] = x[i]+a[k];
         if (mark[i]==0&&x[i]!=0.0) {
            mark[i] = 1;
            nzrow = nzrow+1;
            nzloc[nzrow] = i;
         }
      }
      
      
      for(k = 1; k <= nzrow; k += 1) {
         i = nzloc[k];
         mark[i] = 0;
         xi = x[i];
         x[i] = 0.0;
         if (xi!=0.0) {
            nza = nza+1;
            a[nza] = xi;
            colidx[nza] = i;
         }
      }
      jajp1 = rowstr[j+1];
      rowstr[j+1] = nza+rowstr[1];
   }
}
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[])
{
   int nn1;
   int nzrow, nzv, ii, i;
   double vecelt, vecloc;

   nzv = 0;
   nzrow = 0;
   nn1 = 1;
   do {
      nn1 = 2*nn1;
   }
   while (nn1<n);
   
   
   
   while (nzv<nz) {
      vecelt = randlc(&tran, amult);
      
      
      vecloc = randlc(&tran, amult);
      i = icnvrt(vecloc, nn1)+1;
      if (i>n)
         ;
      else
         
         
         if (mark[i]==0) {
            mark[i] = 1;
            nzrow = nzrow+1;
            nzloc[nzrow] = i;
            nzv = nzv+1;
            v[nzv] = vecelt;
            iv[nzv] = i;
         }
_loop_end_2:      ;
   }

   for(ii = 1; ii <= nzrow; ii += 1) {
      i = nzloc[ii];
      mark[i] = 0;
   }
}
static int icnvrt(double x, int ipwr2)
{
   return (int) (ipwr2*x);
}
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val)
{
   int k;
   boolean set;

   set = 0;
   for(k = 1; k <= *nzv; k += 1)
      if (iv[k]==i) {
         v[k] = val;
         set = 1;
      }
   if (set==0) {
      *nzv = *nzv+1;
      v[*nzv] = val;
      iv[*nzv] = i;
   }
}
