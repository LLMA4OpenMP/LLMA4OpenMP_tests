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

void wtime_(double *t) {
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, (void *) 0);
   if(sec < 0) sec = tv.tv_sec;
   *t = (tv.tv_sec - sec) + 1.0e-6 * tv.tv_usec;
}

void wtime_(double *);
double elapsed_time() {
   double t;
   wtime_(&t);
   
   return (t);
}

double start[64];
double elapsed[64];
void timer_clear(int n) {
   elapsed[n] = 0.0;
}

void timer_start(int n) {
   start[n] = elapsed_time();
}

void timer_stop(int n) {
   double t, now;
   now = elapsed_time();
   t = now - start[n];
   elapsed[n] += t;
}

double timer_read(int n) {
   
   return (elapsed[n]);
}

double randlc(double *x, double a) {
   double t1, t2, t3, t4, a1, a2, x1, x2, z;
   t1 = (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * a;
   a1 = (int) t1;
   a2 = a - (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * a1;
   t1 = (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * (*x);
   x1 = (int) t1;
   x2 = (*x) - (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * x1;
   t1 = a1 * x2 + a2 * x1;
   t2 = (int) ((0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * t1);
   z = t1 - (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * t2;
   t3 = (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * z + a2 * x2;
   t4 = (int) (((0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5)) * t3);
   (*x) = t3 - ((2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0)) * t4;
   
   return (((0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5)) * (*x));
}

void vranlc(int n, double *x_seed, double a, double y[]) {
   int i;
   double x, t1, t2, t3, t4, a1, a2, x1, x2, z;
   t1 = (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * a;
   a1 = (int) t1;
   a2 = a - (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * a1;
   x = *x_seed;
   /*************** Clava msgError **************
   		Variable x could not be categorized into any OpenMP Variable Scopeuse : RWR
   ****************************************/
   for(i = 1; i <= n; i++) {
      t1 = (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * x;
      x1 = (int) t1;
      x2 = x - (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * x1;
      t1 = a1 * x2 + a2 * x1;
      t2 = (int) ((0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * t1);
      z = t1 - (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * t2;
      t3 = (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * z + a2 * x2;
      t4 = (int) (((0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5)) * t3);
      x = t3 - ((2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0)) * t4;
      y[i] = ((0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5)) * x;
   }
   *x_seed = x;
}

void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand) {
   char *evalue = "1000";
   printf("\n\n %s Benchmark Completed\n", name);
   printf(" Class           =                        %c\n", class);
   if(n2 == 0 && n3 == 0) printf(" Size            =             %12d\n", n1);
   else printf(" Size            =              %3dx%3dx%3d\n", n1, n2, n3);
   printf(" Iterations      =             %12d\n", niter);
   printf(" Threads         =             %12d\n", nthreads);
   printf(" Time in seconds =             %12.2f\n", t);
   printf(" Mop/s total     =             %12.2f\n", mops);
   printf(" Operation type  = %24s\n", optype);
   if(passed_verification) printf(" Verification    =               SUCCESSFUL\n");
   else printf(" Verification    =             UNSUCCESSFUL\n");
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

/*CLASS = A*/

/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class of the NPB
c  in this directory. Do not modify it by hand.
*/

static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;
static int colidx[2198001];
static int rowstr[14002];
static int iv[28002];
static int arow[2198001];
static int acol[2198001];
static double v[14002];
static double aelt[2198001];
static double a[2198001];
static double x[14003];
static double z[14003];
static double p[14003];
static double q[14003];
static double r[14003];
static double amult;
static double tran;
static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm);
static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift);
static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza);
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);
int main(int argc, char **argv) {
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
   if(14000 == 1400 && 11 == 7 && 15 == 15 && 20.0 == 10.0) {
      class = 'S';
      zeta_verify_value = 8.5971775078648;
   }
   else if(14000 == 7000 && 11 == 8 && 15 == 15 && 20.0 == 12.0) {
         class = 'W';
         zeta_verify_value = 10.362595087124;
      }
      else if(14000 == 14000 && 11 == 11 && 15 == 15 && 20.0 == 20.0) {
            class = 'A';
            zeta_verify_value = 17.130235054029;
         }
         else if(14000 == 75000 && 11 == 13 && 15 == 75 && 20.0 == 60.0) {
               class = 'B';
               zeta_verify_value = 22.712745482631;
            }
            else if(14000 == 150000 && 11 == 15 && 15 == 75 && 20.0 == 110.0) {
                  class = 'C';
                  zeta_verify_value = 28.973605592845;
               }
               else {
                  class = 'U';
               }
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - CG Benchmark\n");
   printf(" Size: %10d\n", 14000);
   printf(" Iterations: %5d\n", 15);
   naa = 14000;
   nzz = 14000 * (11 + 1) * (11 + 1) + 14000 * (11 + 2);
   tran = 314159265.0;
   amult = 1220703125.0;
   zeta = randlc(&tran, amult);
   makea(naa, nzz, a, colidx, rowstr, 11, firstrow, lastrow, firstcol, lastcol, 1.0e-1, arow, acol, aelt, v, iv, 20.0);
   {
      #pragma omp parallel for default(shared) private(j, k) firstprivate(lastrow, firstrow, firstcol, rowstr) reduction(- : colidx[:2198001])
      for(j = 1; j <= lastrow - firstrow + 1; j++) {
         #pragma omp parallel for default(shared) private(k) firstprivate(j, firstcol, rowstr)
         for(k = rowstr[j]; k < rowstr[j + 1]; k++) {
            colidx[k] = colidx[k] - firstcol + 1;
         }
      }
      #pragma omp parallel for default(shared) private(i)
      for(i = 1; i <= 14000 + 1; i++) {
         x[i] = 1.0;
      }
      #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol)
      for(j = 1; j <= lastcol - firstcol + 1; j++) {
         q[j] = 0.0;
         z[j] = 0.0;
         r[j] = 0.0;
         p[j] = 0.0;
      }
   }
   zeta = 0.0;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(it = 1; it <= 1; it++) {
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      norm_temp11 = 0.0;
      norm_temp12 = 0.0;
      #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, x, z) reduction(+ : norm_temp11) reduction(+ : norm_temp12)
      for(j = 1; j <= lastcol - firstcol + 1; j++) {
         norm_temp11 = norm_temp11 + x[j] * z[j];
         norm_temp12 = norm_temp12 + z[j] * z[j];
      }
      norm_temp12 = 1.0 / sqrt(norm_temp12);
      #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, norm_temp12, z)
      for(j = 1; j <= lastcol - firstcol + 1; j++) {
         x[j] = norm_temp12 * z[j];
      }
   }
   #pragma omp parallel for default(shared) private(i)
   for(i = 1; i <= 14000 + 1; i++) {
      x[i] = 1.0;
   }
   zeta = 0.0;
   timer_clear(1);
   timer_start(1);
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(it = 1; it <= 15; it++) {
      conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);
      norm_temp11 = 0.0;
      norm_temp12 = 0.0;
      #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, x, z) reduction(+ : norm_temp11) reduction(+ : norm_temp12)
      for(j = 1; j <= lastcol - firstcol + 1; j++) {
         norm_temp11 = norm_temp11 + x[j] * z[j];
         norm_temp12 = norm_temp12 + z[j] * z[j];
      }
      norm_temp12 = 1.0 / sqrt(norm_temp12);
      zeta = 20.0 + 1.0 / norm_temp11;
      if(it == 1) {
         printf("   iteration           ||r||                 zeta\n");
      }
      printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);
      #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, norm_temp12, z)
      for(j = 1; j <= lastcol - firstcol + 1; j++) {
         x[j] = norm_temp12 * z[j];
      }
   }
   {
   }
   timer_stop(1);
   t = timer_read(1);
   printf(" Benchmark completed\n");
   epsilon = 1.0e-10;
   if(class != 'U') {
      if(fabs(zeta - zeta_verify_value) <= epsilon) {
         verified = 1;
         printf(" VERIFICATION SUCCESSFUL\n");
         printf(" Zeta is    %20.12e\n", zeta);
         printf(" Error is   %20.12e\n", zeta - zeta_verify_value);
      }
      else {
         verified = 0;
         printf(" VERIFICATION FAILED\n");
         printf(" Zeta                %20.12e\n", zeta);
         printf(" The correct zeta is %20.12e\n", zeta_verify_value);
      }
   }
   else {
      verified = 0;
      printf(" Problem size unknown\n");
      printf(" NO VERIFICATION PERFORMED\n");
   }
   if(t != 0.0) {
      mflops = (2.0 * 15 * 14000) * (3.0 + (11 * (11 + 1)) + 25.0 * (5.0 + (11 * (11 + 1))) + 3.0) / t / 1000000.0;
   }
   else {
      mflops = 0.0;
   }
   c_print_results("CG", class, 14000, 0, 0, 15, nthreads, t, mflops, "          floating point", verified, "3.0 structured", "25 Nov 2023", "gcc-12", "gcc-12", "-fopenmp -lm", "-I../common -fopenmp", "-lm", "(none)", "randdp");
}

static void conj_grad(int colidx[], int rowstr[], double x[], double z[], double a[], double p[], double q[], double r[], double *rnorm) {
   static int callcount = 0;
   double d, sum, rho, rho0, alpha, beta;
   int i, j, k;
   int cgit, cgitmax = 25;
   rho = 0.0;
   {
      #pragma omp parallel for default(shared) private(j) firstprivate(naa, x)
      for(j = 1; j <= naa + 1; j++) {
         q[j] = 0.0;
         z[j] = 0.0;
         r[j] = x[j];
         p[j] = r[j];
      }
      #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, r) reduction(+ : rho)
      for(j = 1; j <= lastcol - firstcol + 1; j++) {
         rho = rho + r[j] * r[j];
      }
   }
   /*************** Clava msgError **************
   		Variable rho could not be categorized into any OpenMP Variable Scopeuse : RWR
   ****************************************/
   for(cgit = 1; cgit <= cgitmax; cgit++) {
      rho0 = rho;
      d = 0.0;
      rho = 0.0;
      {
         #pragma omp parallel for default(shared) private(j, k, sum) firstprivate(lastrow, firstrow, rowstr, a, colidx, p)
         for(j = 1; j <= lastrow - firstrow + 1; j++) {
            sum = 0.0;
            #pragma omp parallel for default(shared) private(k) firstprivate(j, rowstr, a, colidx, p) reduction(+ : sum)
            for(k = rowstr[j]; k < rowstr[j + 1]; k++) {
               sum = sum + a[k] * p[colidx[k]];
            }
            q[j] = sum;
         }
         #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, p, q) reduction(+ : d)
         for(j = 1; j <= lastcol - firstcol + 1; j++) {
            d = d + p[j] * q[j];
         }
         alpha = rho0 / d;
         #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, alpha, p, q) reduction(+ : rho)
         for(j = 1; j <= lastcol - firstcol + 1; j++) {
            z[j] = z[j] + alpha * p[j];
            r[j] = r[j] - alpha * q[j];
            rho = rho + r[j] * r[j];
         }
         beta = rho / rho0;
         #pragma omp parallel for default(shared) private(j) firstprivate(lastcol, firstcol, beta, r)
         for(j = 1; j <= lastcol - firstcol + 1; j++) {
            p[j] = r[j] + beta * p[j];
         }
         callcount++;
      }
   }
   sum = 0.0;
   {
      #pragma omp parallel for default(shared) private(j, k, d) firstprivate(lastrow, firstrow, rowstr, a, colidx, z)
      for(j = 1; j <= lastrow - firstrow + 1; j++) {
         d = 0.0;
         #pragma omp parallel for default(shared) private(k) firstprivate(j, rowstr, a, colidx, z) reduction(+ : d)
         for(k = rowstr[j]; k <= rowstr[j + 1] - 1; k++) {
            d = d + a[k] * z[colidx[k]];
         }
         r[j] = d;
      }
      #pragma omp parallel for default(shared) private(j, d) firstprivate(lastcol, firstcol, x, r) reduction(+ : sum)
      for(j = 1; j <= lastcol - firstcol + 1; j++) {
         d = x[j] - r[j];
         sum = sum + d * d;
      }
   }
   (*rnorm) = sqrt(sum);
}

static void makea(int n, int nz, double a[], int colidx[], int rowstr[], int nonzer, int firstrow, int lastrow, int firstcol, int lastcol, double rcond, int arow[], int acol[], double aelt[], double v[], int iv[], double shift) {
   int i, nnza, iouter, ivelt, ivelt1, irow, nzv;
   double size, ratio, scale;
   int jcol;
   size = 1.0;
   ratio = pow(rcond, (1.0 / (double) n));
   nnza = 0;
   #pragma omp parallel for default(shared) private(i) firstprivate(n)
   for(i = 1; i <= n; i++) {
      colidx[n + i] = 0;
   }
   /*************** Clava msgError **************
   		Loop contains Invalid Statement -> exit#467
   ****************************************/
   for(iouter = 1; iouter <= n; iouter++) {
      nzv = nonzer;
      sprnvc(n, nzv, v, iv, &(colidx[0]), &(colidx[n]));
      vecset(n, v, iv, &nzv, iouter, 0.5);
      /*************** Clava msgError **************
      		Loop contains Invalid Statement -> exit#467
      ****************************************/
      for(ivelt = 1; ivelt <= nzv; ivelt++) {
         jcol = iv[ivelt];
         if(jcol >= firstcol && jcol <= lastcol) {
            scale = size * v[ivelt];
            /*************** Clava msgError **************
            		Loop contains Invalid Statement -> exit#467
            ****************************************/
            for(ivelt1 = 1; ivelt1 <= nzv; ivelt1++) {
               irow = iv[ivelt1];
               if(irow >= firstrow && irow <= lastrow) {
                  nnza = nnza + 1;
                  if(nnza > nz) {
                     printf("Space for matrix elements exceeded in makea\n");
                     printf("nnza, nzmax = %d, %d\n", nnza, nz);
                     printf("iouter = %d\n", iouter);
                     exit(1);
                  }
                  acol[nnza] = jcol;
                  arow[nnza] = irow;
                  aelt[nnza] = v[ivelt1] * scale;
               }
            }
         }
      }
      size = size * ratio;
   }
   /*************** Clava msgError **************
   		Loop contains Invalid Statement -> exit#488
   ****************************************/
   for(i = firstrow; i <= lastrow; i++) {
      if(i >= firstcol && i <= lastcol) {
         iouter = n + i;
         nnza = nnza + 1;
         if(nnza > nz) {
            printf("Space for matrix elements exceeded in makea\n");
            printf("nnza, nzmax = %d, %d\n", nnza, nz);
            printf("iouter = %d\n", iouter);
            exit(1);
         }
         acol[nnza] = i;
         arow[nnza] = i;
         aelt[nnza] = rcond - shift;
      }
   }
   sparse(a, colidx, rowstr, n, arow, acol, aelt, firstrow, lastrow, v, &(iv[0]), &(iv[n]), nnza);
}

static void sparse(double a[], int colidx[], int rowstr[], int n, int arow[], int acol[], double aelt[], int firstrow, int lastrow, double x[], boolean mark[], int nzloc[], int nnza) {
   int nrows;
   int i, j, jajp1, nza, k, nzrow;
   double xi;
   nrows = lastrow - firstrow + 1;
   #pragma omp parallel for default(shared) private(j) firstprivate(n)
   for(j = 1; j <= n; j++) {
      rowstr[j] = 0;
      mark[j] = 0;
   }
   rowstr[n + 1] = 0;
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess rowstr	 use : RW
   ****************************************/
   for(nza = 1; nza <= nnza; nza++) {
      j = (arow[nza] - firstrow + 1) + 1;
      rowstr[j] = rowstr[j] + 1;
   }
   rowstr[1] = 1;
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess rowstr	 use : RW
   ****************************************/
   for(j = 2; j <= nrows + 1; j++) {
      rowstr[j] = rowstr[j] + rowstr[j - 1];
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess a	 use : W
   ****************************************/
   for(j = 0; j <= nrows - 1; j++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(j, rowstr)
      for(k = rowstr[j]; k <= rowstr[j + 1] - 1; k++)
         a[k] = 0.0;
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess rowstr	 use : RW
   		unsolved dependency for arrayAccess a	 use : W
   		unsolved dependency for arrayAccess colidx	 use : W
   ****************************************/
   for(nza = 1; nza <= nnza; nza++) {
      j = arow[nza] - firstrow + 1;
      k = rowstr[j];
      a[k] = aelt[nza];
      colidx[k] = acol[nza];
      rowstr[j] = rowstr[j] + 1;
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess rowstr	 use : RW
   ****************************************/
   for(j = nrows; j >= 1; j--) {
      rowstr[j + 1] = rowstr[j];
   }
   rowstr[1] = 1;
   nza = 0;
   #pragma omp parallel for default(shared) private(i) firstprivate(n)
   for(i = 1; i <= n; i++) {
      x[i] = 0.0;
      mark[i] = 0;
   }
   jajp1 = rowstr[1];
   /*************** Clava msgError **************
   		Variable jajp1 could not be categorized into any OpenMP Variable Scopeuse : RW
   		Variable nza could not be categorized into any OpenMP Variable Scopeuse : RWR
   ****************************************/
   for(j = 1; j <= nrows; j++) {
      nzrow = 0;
      /*************** Clava msgError **************
      		Variable nzrow could not be categorized into any OpenMP Variable Scopeuse : RWR
      ****************************************/
      for(k = jajp1; k < rowstr[j + 1]; k++) {
         i = colidx[k];
         x[i] = x[i] + a[k];
         if(mark[i] == 0 && x[i] != 0.0) {
            mark[i] = 1;
            nzrow = nzrow + 1;
            nzloc[nzrow] = i;
         }
      }
      /*************** Clava msgError **************
      		Variable nza could not be categorized into any OpenMP Variable Scopeuse : RWR
      ****************************************/
      for(k = 1; k <= nzrow; k++) {
         i = nzloc[k];
         mark[i] = 0;
         xi = x[i];
         x[i] = 0.0;
         if(xi != 0.0) {
            nza = nza + 1;
            a[nza] = xi;
            colidx[nza] = i;
         }
      }
      jajp1 = rowstr[j + 1];
      rowstr[j + 1] = nza + rowstr[1];
   }
}

static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]) {
   int nn1;
   int nzrow, nzv, ii, i;
   double vecelt, vecloc;
   nzv = 0;
   nzrow = 0;
   nn1 = 1;
   do  {
      nn1 = 2 * nn1;
   }
   while (nn1 < n);
   while(nzv < nz) {
      vecelt = randlc(&tran, amult);
      vecloc = randlc(&tran, amult);
      i = icnvrt(vecloc, nn1) + 1;
      if(i > n) continue;
      if(mark[i] == 0) {
         mark[i] = 1;
         nzrow = nzrow + 1;
         nzloc[nzrow] = i;
         nzv = nzv + 1;
         v[nzv] = vecelt;
         iv[nzv] = i;
      }
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess mark	 use : W
   ****************************************/
   for(ii = 1; ii <= nzrow; ii++) {
      i = nzloc[ii];
      mark[i] = 0;
   }
}

static int icnvrt(double x, int ipwr2) {
   
   return ((int) (ipwr2 * x));
}

static void vecset(int n, double v[], int iv[], int *nzv, int i, double val) {
   int k;
   boolean set;
   set = 0;
   /*************** Clava msgError **************
   		 Variable Access set is changed inside  of ifstmt
   ****************************************/
   for(k = 1; k <= *nzv; k++) {
      if(iv[k] == i) {
         v[k] = val;
         set = 1;
      }
   }
   if(set == 0) {
      *nzv = *nzv + 1;
      v[*nzv] = val;
      iv[*nzv] = i;
   }
}
