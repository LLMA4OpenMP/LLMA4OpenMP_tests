/*
 * file for cg.c
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



/* global variables */

/* common /partit_size/ */
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;
static int colidx[150000*(15+1)*(15+1)+150000*(15+2)+1];
/* colidx[1:NZ] */
static int rowstr[150000+1+1];
/* rowstr[1:NA+1] */
static int iv[2*150000+1+1];
static int arow[150000*(15+1)*(15+1)+150000*(15+2)+1];
static int acol[150000*(15+1)*(15+1)+150000*(15+2)+1];
/* acol[1:NZ] */

/* common /main_flt_mem/ */
static double v[150000+1+1];
static double aelt[150000*(15+1)*(15+1)+150000*(15+2)+1];
static double a[150000*(15+1)*(15+1)+150000*(15+2)+1];
/* a[1:NZ] */
static double x[150000+2+1];
/* x[1:NA+2] */
static double z[150000+2+1];
/* z[1:NA+2] */
static double p[150000+2+1];
/* p[1:NA+2] */
static double q[150000+2+1];
/* q[1:NA+2] */
static double r[150000+2+1];
/* r[1:NA+2] */
//static double w[NA+2+1];	/* w[1:NA+2] */

/* common /urando/ */
static double amult;
static double tran;

/* function declarations */


//double w[],
static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm);



static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift);



static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza);

static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);

/*--------------------------------------------------------------------
      program cg
--------------------------------------------------------------------*/

int main(int argc, char **argv);


/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/
static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm);


/*---------------------------------------------------------------------
c       generate the test problem for benchmark 6
c       makea generates a sparse matrix with a
c       prescribed sparsity distribution
c
c       parameter    type        usage
c
c       input
c
c       n            i           number of cols/rows of matrix
c       nz           i           nonzeros as declared array size
c       rcond        r*8         condition number
c       shift        r*8         main diagonal shift
c
c       output
c
c       a            r*8         array for nonzeros
c       colidx       i           col indices
c       rowstr       i           row pointers
c
c       workspace
c
c       iv, arow, acol i
c       v, aelt        r*8
c---------------------------------------------------------------------*/
static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift);


/*---------------------------------------------------
c       generate a sparse matrix from a list of
c       [col, row, element] tri
c---------------------------------------------------*/
static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza);


/*---------------------------------------------------------------------
c       generate a sparse n-vector (v, iv)
c       having nzv nonzeros
c
c       mark(i) is set to 1 if position i is nonzero.
c       mark is all zero on entry and is reset to all zero before exit
c       this corrects a performance bug found by John G. Lewis, caused by
c       reinitialization of mark on every one of the n calls to sprnvc
---------------------------------------------------------------------*/
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]);


/*---------------------------------------------------------------------
* scale a double precision number x in (0,1) by a power of 2 and chop it
*---------------------------------------------------------------------*/
static int icnvrt(double x, int ipwr2);


/*--------------------------------------------------------------------
c       set ith element of sparse vector (v, iv) with
c       nzv nonzeros to val
c-------------------------------------------------------------------*/
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
   lastrow = 150000;
   firstcol = 1;
   lastcol = 150000;

   if (150000==1400&&15==7&&75==15&&110.0==10.0) {
      class = 'S';
      zeta_verify_value = 8.5971775078648;
   }
   else if (150000==7000&&15==8&&75==15&&110.0==12.0) {
      class = 'W';
      zeta_verify_value = 10.362595087124;
   }
   else if (150000==14000&&15==11&&75==15&&110.0==20.0) {
      class = 'A';
      zeta_verify_value = 17.130235054029;
   }
   else if (150000==75000&&15==13&&75==75&&110.0==60.0) {
      class = 'B';
      zeta_verify_value = 22.712745482631;
   }
   else if (150000==150000&&15==15&&75==75&&110.0==110.0) {
      class = 'C';
      zeta_verify_value = 28.973605592845;
   }
   else
      class = 'U';
   
   
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"" - CG Benchmark\n");
   printf(" Size: %10d\n", 150000);
   printf(" Iterations: %5d\n", 75);

   naa = 150000;
   nzz = 150000*(15+1)*(15+1)+150000*(15+2);
   
   /*--------------------------------------------------------------------
   c  Initialize random number generator
   c-------------------------------------------------------------------*/
   tran = 314159265.0;
   amult = 1220703125.0;
   zeta = randlc(&tran, amult);
   
   /*--------------------------------------------------------------------
   c  
   c-------------------------------------------------------------------*/
   
   
   makea(naa, nzz, a, colidx, rowstr, 15, firstrow, lastrow, firstcol, lastcol, 1.0e-1, arow, acol, aelt, v, iv, 110.0);
   for(j = 1; j <= lastrow-firstrow+1; j += 1)
#pragma omp parallel for
      for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
         colidx[k] = colidx[k]-firstcol+1;
   
   /*--------------------------------------------------------------------
   c  set starting vector to (1, 1, .... 1)
   c-------------------------------------------------------------------*/
#pragma omp parallel for
   for(i = 1; i <= 150000+1; i += 1)
      x[i] = 1.0;
#pragma omp parallel for
   for(j = 1; j <= lastcol-firstcol+1; j += 1) {
      q[j] = 0.0;
      z[j] = 0.0;
      r[j] = 0.0;
      p[j] = 0.0;
   }
   // end omp parallel
   zeta = 0.0;
   
   /*-------------------------------------------------------------------
   c---->
   c  Do one iteration untimed to init all code and data page tables
   c---->                    (then reinit, start timing, to niter its)
   c-------------------------------------------------------------------*/
   
   for(it = 1; it <= 1; it += 1) {
      
      /*--------------------------------------------------------------------
      c  The call to the conjugate gradient routine:
      c-------------------------------------------------------------------*/
      /* w,*/
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      
      /*--------------------------------------------------------------------
      c  zeta = shift + 1/(x.z)
      c  So, first: (x.z)
      c  Also, find norm of z
      c  So, first: (z.z)
      c-------------------------------------------------------------------*/
      norm_temp11 = 0.0;
      norm_temp12 = 0.0;
#pragma omp parallel for reduction(+:norm_temp12) reduction(+:norm_temp11)
      for(j = 1; j <= lastcol-firstcol+1; j += 1) {
         norm_temp11 = norm_temp11+x[j]*z[j];
         norm_temp12 = norm_temp12+z[j]*z[j];
      }
      norm_temp12 = 1.0/sqrt(norm_temp12);
      
      /*--------------------------------------------------------------------
      c  Normalize z to obtain x
      c-------------------------------------------------------------------*/
#pragma omp parallel for
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         x[j] = norm_temp12*z[j];
   }
   /* end of do one iteration untimed */
   
   /*--------------------------------------------------------------------
   c  set starting vector to (1, 1, .... 1)
   c-------------------------------------------------------------------*/
#pragma omp parallel for
   for(i = 1; i <= 150000+1; i += 1)
      x[i] = 1.0;
   zeta = 0.0;
   
   
   timer_clear(1);
   timer_start(1);
   
   /*--------------------------------------------------------------------
   c---->
   c  Main Iteration for inverse power method
   c---->
   c-------------------------------------------------------------------*/
   
   for(it = 1; it <= 75; it += 1) {
      
      /*--------------------------------------------------------------------
      c  The call to the conjugate gradient routine:
      c-------------------------------------------------------------------*/
      /*, w*/
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      
      /*--------------------------------------------------------------------
      c  zeta = shift + 1/(x.z)
      c  So, first: (x.z)
      c  Also, find norm of z
      c  So, first: (z.z)
      c-------------------------------------------------------------------*/
      norm_temp11 = 0.0;
      norm_temp12 = 0.0;

#pragma omp parallel for reduction(+:norm_temp12) reduction(+:norm_temp11)
      for(j = 1; j <= lastcol-firstcol+1; j += 1) {
         norm_temp11 = norm_temp11+x[j]*z[j];
         norm_temp12 = norm_temp12+z[j]*z[j];
      }

      norm_temp12 = 1.0/sqrt(norm_temp12);

      zeta = 110.0+1.0/norm_temp11;

      if (it==1)
         printf("   iteration           ||r||                 zeta\n");
      printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);
      
      /*--------------------------------------------------------------------
      c  Normalize z to obtain x
      c-------------------------------------------------------------------*/
#pragma omp parallel for
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         x[j] = norm_temp12*z[j];
   }
   /* end parallel */
   
   timer_stop(1);
   
   /*--------------------------------------------------------------------
   c  End of timed section
   c-------------------------------------------------------------------*/
   
   t = timer_read(1);

   printf(" Benchmark completed\n");

   epsilon = 1.0e-10;
   //epsilon = 1.0e-2;
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
      mflops = 2.0*75*150000*(3.0+15*(15+1)+25.0*(5.0+15*(15+1))+3.0)/t/1000000.0;
   else
      mflops = 0.0;
   
   
   
   
   c_print_results("CG", class, 150000, 0, 0, 75, nthreads, t, mflops, "          floating point", verified, "3.0 structured", "01 Dec 2023", "(none)", "(none)", "-lm", "(none)", "(none)", "(none)", "randdp");
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
   
   /*--------------------------------------------------------------------
   c  rho = r.r
   c  Now, obtain the norm of r: First, sum squares of r elements locally...
   c-------------------------------------------------------------------*/
#pragma omp parallel for reduction(+:rho)
   for(j = 1; j <= lastcol-firstcol+1; j += 1)
      rho = rho+r[j]*r[j];
   /* end omp parallel */
   /*--------------------------------------------------------------------
   c---->
   c  The conj grad iteration loop
   c---->
   c-------------------------------------------------------------------*/
   for(cgit = 1; cgit <= cgitmax; cgit += 1) {
      rho0 = rho;
      d = 0.0;
      rho = 0.0;
      
      /*--------------------------------------------------------------------
      c  q = A.p
      c  The partition submatrix-vector multiply: use workspace w
      c---------------------------------------------------------------------
      C
      C  NOTE: this version of the multiply is actually (slightly: maybe %5) 
      C        faster on the sp2 on 16 nodes than is the unrolled-by-2 version 
      C        below.   On the Cray t3d, the reverse is true, i.e., the 
      C        unrolled-by-two version is some 10% faster.  
      C        The unrolled-by-8 version below is significantly faster
      C        on the Cray t3d - overall speed of code is 1.5 times faster.
      */
      
      /* rolled version */
#pragma omp parallel for private(sum, k)
      for(j = 1; j <= lastrow-firstrow+1; j += 1) {
         sum = 0.0;
         for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
            sum = sum+a[k]*p[colidx[k]];
         //w[j] = sum;
         q[j] = sum;
      }
      
      /* unrolled-by-two version
              for (j = 1; j <= lastrow-firstrow+1; j++) {
      	    int iresidue;
      	    double sum1, sum2;
      	    i = rowstr[j]; 
                  iresidue = (rowstr[j+1]-i) % 2;
                  sum1 = 0.0;
                  sum2 = 0.0;
                  if (iresidue == 1) sum1 = sum1 + a[i]*p[colidx[i]];
      	    for (k = i+iresidue; k <= rowstr[j+1]-2; k += 2) {
      		sum1 = sum1 + a[k]   * p[colidx[k]];
      		sum2 = sum2 + a[k+1] * p[colidx[k+1]];
      	    }
                  w[j] = sum1 + sum2;
              }
      */
      /* unrolled-by-8 version
              for (j = 1; j <= lastrow-firstrow+1; j++) {
      	    int iresidue;
                  i = rowstr[j]; 
                  iresidue = (rowstr[j+1]-i) % 8;
                  sum = 0.0;
                  for (k = i; k <= i+iresidue-1; k++) {
                      sum = sum +  a[k] * p[colidx[k]];
                  }
                  for (k = i+iresidue; k <= rowstr[j+1]-8; k += 8) {
                      sum = sum + a[k  ] * p[colidx[k  ]]
                                + a[k+1] * p[colidx[k+1]]
                                + a[k+2] * p[colidx[k+2]]
                                + a[k+3] * p[colidx[k+3]]
                                + a[k+4] * p[colidx[k+4]]
                                + a[k+5] * p[colidx[k+5]]
                                + a[k+6] * p[colidx[k+6]]
                                + a[k+7] * p[colidx[k+7]];
                  }
                  w[j] = sum;
              }
      */
      /*	
      	for (j = 1; j <= lastcol-firstcol+1; j++) {
                  q[j] = w[j];
      	}
      */
      /*--------------------------------------------------------------------
      c  Clear w for reuse...
      c-------------------------------------------------------------------*/
      /*
      	for (j = 1; j <= lastcol-firstcol+1; j++) {
                  w[j] = 0.0;
      	}
      */
      /*--------------------------------------------------------------------
      c  Obtain p.q
      c-------------------------------------------------------------------*/
#pragma omp parallel for reduction(+:d)
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         d = d+p[j]*q[j];
      /*--------------------------------------------------------------------
      c  Obtain alpha = rho / (p.q)
      c-------------------------------------------------------------------*/
      alpha = rho0/d;
      
      /*--------------------------------------------------------------------
      c  Save a temporary of rho
      c-------------------------------------------------------------------*/
      /*	rho0 = rho;*/
      
      /*---------------------------------------------------------------------
      c  Obtain z = z + alpha*p
      c  and    r = r - alpha*q
      c---------------------------------------------------------------------*/
      for(j = 1; j <= lastcol-firstcol+1; j += 1) {
         z[j] = z[j]+alpha*p[j];
         r[j] = r[j]-alpha*q[j];
         //	}
         
         /*---------------------------------------------------------------------
         c  rho = r.r
         c  Now, obtain the norm of r: First, sum squares of r elements locally...
         c---------------------------------------------------------------------*/
         /*
         	for (j = 1; j <= lastcol-firstcol+1; j++) {*/
         rho = rho+r[j]*r[j];
      }
      
      /*--------------------------------------------------------------------
      c  Obtain beta:
      c-------------------------------------------------------------------*/
      beta = rho/rho0;
      
      /*--------------------------------------------------------------------
      c  p = r + beta*p
      c-------------------------------------------------------------------*/
#pragma omp parallel for
      for(j = 1; j <= lastcol-firstcol+1; j += 1)
         p[j] = r[j]+beta*p[j];
      callcount++;
   }
   /* end of do cgit=1,cgitmax */
   
   /*---------------------------------------------------------------------
   c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
   c  First, form A.z
   c  The partition submatrix-vector multiply
   c---------------------------------------------------------------------*/
   sum = 0.0;
#pragma omp parallel for private(k, d)
   for(j = 1; j <= lastrow-firstrow+1; j += 1) {
      d = 0.0;
      for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
         d = d+a[k]*z[colidx[k]];
      r[j] = d;
   }
   
   /*--------------------------------------------------------------------
   c  At this point, r contains A.z
   c-------------------------------------------------------------------*/
#pragma omp parallel for private(d) reduction(+:sum)
   for(j = 1; j <= lastcol-firstcol+1; j += 1) {
      d = x[j]-r[j];
      sum = sum+d*d;
   }
   //end omp parallel
   *rnorm = sqrt(sum);
}
static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift)
{
   int i, nnza, iouter, ivelt, ivelt1, irow, nzv;
   
   /*--------------------------------------------------------------------
   c      nonzer is approximately  (int(sqrt(nnza /n)));
   c-------------------------------------------------------------------*/
   
   double size, ratio, scale;
   int jcol;

   size = 1.0;
   ratio = pow(rcond, 1.0/((double) n));
   nnza = 0;
   
   /*---------------------------------------------------------------------
   c  Initialize colidx(n+1 .. 2n) to zero.
   c  Used by sprnvc to mark nonzero positions
   c---------------------------------------------------------------------*/
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
   
   /*---------------------------------------------------------------------
   c       ... add the identity * rcond to the generated matrix to bound
   c           the smallest eigenvalue from below by rcond
   c---------------------------------------------------------------------*/
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
   
   /*---------------------------------------------------------------------
   c       ... make the sparse matrix from list of elements with duplicates
   c           (v and iv are used as  workspace)
   c---------------------------------------------------------------------*/
   
   sparse(a, colidx, rowstr, n, arow, acol, aelt, firstrow, lastrow, v, &iv[0], &iv[n], nnza);
}
static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza)
{
   int nrows;
   int i, j, jajp1, nza, k, nzrow;
   double xi;
   
   /*--------------------------------------------------------------------
   c    how many rows of result
   c-------------------------------------------------------------------*/
   nrows = lastrow-firstrow+1;
   
   /*--------------------------------------------------------------------
   c     ...count the number of triples in each row
   c-------------------------------------------------------------------*/
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
   
   /*---------------------------------------------------------------------
   c     ... rowstr(j) now is the location of the first nonzero
   c           of row j of a
   c---------------------------------------------------------------------*/
   
   /*---------------------------------------------------------------------
   c     ... preload data pages
   c---------------------------------------------------------------------*/
   for(j = 0; j <= nrows-1; j += 1)
#pragma omp parallel for
      for(k = rowstr[j]; k <= rowstr[j+1]-1; k += 1)
         a[k] = 0.0;
   /*--------------------------------------------------------------------
   c     ... do a bucket sort of the triples on the row index
   c-------------------------------------------------------------------*/
   for(nza = 1; nza <= nnza; nza += 1) {
      j = arow[nza]-firstrow+1;
      k = rowstr[j];
      a[k] = aelt[nza];
      colidx[k] = acol[nza];
      rowstr[j] = rowstr[j]+1;
   }
   
   /*--------------------------------------------------------------------
   c       ... rowstr(j) now points to the first element of row j+1
   c-------------------------------------------------------------------*/
   for(j = nrows; j >= 1; j += -1)
      rowstr[j+1] = rowstr[j];
   rowstr[1] = 1;
   
   /*--------------------------------------------------------------------
   c       ... generate the actual output rows by adding elements
   c-------------------------------------------------------------------*/
   nza = 0;

   jajp1 = rowstr[1];
   for(j = 1; j <= nrows; j += 1) {
      nzrow = 0;
      
      /*--------------------------------------------------------------------
      c          ...loop over the jth row of a
      c-------------------------------------------------------------------*/
      for(k = jajp1; k <= rowstr[j+1]-1; k += 1) {
         i = colidx[k];
         x[i] = x[i]+a[k];
         if (mark[i]==0&&x[i]!=0.0) {
            mark[i] = 1;
            nzrow = nzrow+1;
            nzloc[nzrow] = i;
         }
      }
      
      /*--------------------------------------------------------------------
      c          ... extract the nonzeros of this row
      c-------------------------------------------------------------------*/
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
   
   /*--------------------------------------------------------------------
   c    nn1 is the smallest power of two not less than n
   c-------------------------------------------------------------------*/
   
   while (nzv<nz) {
      vecelt = randlc(&tran, amult);
      
      /*--------------------------------------------------------------------
      c   generate an integer between 1 and n in a portable manner
      c-------------------------------------------------------------------*/
      vecloc = randlc(&tran, amult);
      i = icnvrt(vecloc, nn1)+1;
      if (i>n)
         ;
      else
         
         /*--------------------------------------------------------------------
         c  was this integer generated already?
         c-------------------------------------------------------------------*/
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
