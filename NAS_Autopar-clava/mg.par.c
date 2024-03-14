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

/*parameters*/

/*actual dimension including ghost cells for communications*/

/*size of rhs array*/

/*size of residual array*/

/*size of communication buffer*/

/*maximum number of levels*/

/*---------------------------------------------------------------------*/

/*common /mg3/*/

static int nx[12];
static int ny[12];
static int nz[12];
/*common /ClassType/*/

static char Class;
/*common /my_debug/*/

static int debug_vec[8];
/*common /fap/*/

/*static int ir[MAXLEVEL], m1[MAXLEVEL], m2[MAXLEVEL], m3[MAXLEVEL];*/

static int m1[12];
static int m2[12];
static int m3[12];
static int lt;
static int lb;
/*c---------------------------------------------------------------------
c  Set at m=1024, can handle cases up to 1024^3 case
c---------------------------------------------------------------------*/

static int is1;
static int is2;
static int is3;
static int ie1;
static int ie2;
static int ie3;
static void setup(int *n1, int *n2, int *n3, int lt);
static void mg3P(double ****u, double ***v, double ****r, double a[4], double c[4], int n1, int n2, int n3, int k);
static void psinv(double ***r, double ***u, int n1, int n2, int n3, double c[4], int k);
static void resid(double ***u, double ***v, double ***r, int n1, int n2, int n3, double a[4], int k);
static void rprj3(double ***r, int m1k, int m2k, int m3k, double ***s, int m1j, int m2j, int m3j, int k);
static void interp(double ***z, int mm1, int mm2, int mm3, double ***u, int n1, int n2, int n3, int k);
static void norm2u3(double ***r, int n1, int n2, int n3, double *rnm2, double *rnmu, int nx, int ny, int nz);
static void rep_nrm(double ***u, int n1, int n2, int n3, char *title, int kk);
static void comm3(double ***u, int n1, int n2, int n3, int kk);
static void zran3(double ***z, int n1, int n2, int n3, int nx, int ny, int k);
static void showall(double ***z, int n1, int n2, int n3);
static double power(double a, int n);
static void bubble(double ten[1037][2], int j1[1037][2], int j2[1037][2], int j3[1037][2], int m, int ind);
static void zero3(double ***z, int n1, int n2, int n3);
static void nonzero(double ***z, int n1, int n2, int n3);
int main(int argc, char *argv[]) {
   int k, it;
   double t, tinit, mflops;
   int nthreads = 1;
   double ****u;
   double ***v;
   double ****r;
   double a[4];
   double c[4];
   double rnm2, rnmu;
   double epsilon = 1.0e-8;
   int n1, n2, n3, nit;
   double verify_value;
   boolean verified;
   int i, j, l;
   FILE *fp;
   timer_clear(1);
   timer_clear(2);
   timer_start(2);
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - MG Benchmark\n\n");
   fp = fopen("mg.input", "r");
   if(fp != ((void *) 0)) {
      printf(" Reading from input file mg.input\n");
      fscanf(fp, "%d", &lt);
      while(fgetc(fp) != '\n');
      fscanf(fp, "%d%d%d", &nx[lt], &ny[lt], &nz[lt]);
      while(fgetc(fp) != '\n');
      fscanf(fp, "%d", &nit);
      while(fgetc(fp) != '\n');
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(i = 0; i <= 7; i++) {
         fscanf(fp, "%d", &debug_vec[i]);
      }
      fclose(fp);
   }
   else {
      printf(" No input file. Using compiled defaults\n");
      lt = 8;
      nit = 4;
      nx[lt] = 256;
      ny[lt] = 256;
      nz[lt] = 256;
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(i = 0; i <= 7; i++) {
         debug_vec[i] = 0;
      }
   }
   if((nx[lt] != ny[lt]) || (nx[lt] != nz[lt])) {
      Class = 'U';
   }
   else if(nx[lt] == 32 && nit == 4) {
         Class = 'S';
      }
      else if(nx[lt] == 64 && nit == 40) {
            Class = 'W';
         }
         else if(nx[lt] == 256 && nit == 20) {
               Class = 'B';
            }
            else if(nx[lt] == 512 && nit == 20) {
                  Class = 'C';
               }
               else if(nx[lt] == 256 && nit == 4) {
                     Class = 'A';
                  }
                  else {
                     Class = 'U';
                  }
   a[0] = -8.0 / 3.0;
   a[1] = 0.0;
   a[2] = 1.0 / 6.0;
   a[3] = 1.0 / 12.0;
   if(Class == 'A' || Class == 'S' || Class == 'W') {
      c[0] = -3.0 / 8.0;
      c[1] = 1.0 / 32.0;
      c[2] = -1.0 / 64.0;
      c[3] = 0.0;
   }
   else {
      c[0] = -3.0 / 17.0;
      c[1] = 1.0 / 33.0;
      c[2] = -1.0 / 61.0;
      c[3] = 0.0;
   }
   lb = 1;
   setup(&n1, &n2, &n3, lt);
   u = (double ****) malloc((lt + 1) * sizeof(double ***));
   #pragma omp parallel for default(shared) private(l, k, j) firstprivate(lt, m3, m2, m1)
   for(l = lt; l >= 1; l--) {
      u[l] = (double ***) malloc(m3[l] * sizeof(double **));
      #pragma omp parallel for default(shared) private(k, j) firstprivate(l, m3, m2, m1)
      for(k = 0; k < m3[l]; k++) {
         u[l][k] = (double **) malloc(m2[l] * sizeof(double *));
         #pragma omp parallel for default(shared) private(j) firstprivate(l, k, m2, m1)
         for(j = 0; j < m2[l]; j++) {
            u[l][k][j] = (double *) malloc(m1[l] * sizeof(double));
         }
      }
   }
   v = (double ***) malloc(m3[lt] * sizeof(double **));
   #pragma omp parallel for default(shared) private(k, j) firstprivate(lt, m3, m2, m1)
   for(k = 0; k < m3[lt]; k++) {
      v[k] = (double **) malloc(m2[lt] * sizeof(double *));
      #pragma omp parallel for default(shared) private(j) firstprivate(lt, k, m2, m1)
      for(j = 0; j < m2[lt]; j++) {
         v[k][j] = (double *) malloc(m1[lt] * sizeof(double));
      }
   }
   r = (double ****) malloc((lt + 1) * sizeof(double ***));
   #pragma omp parallel for default(shared) private(l, k, j) firstprivate(lt, m3, m2, m1)
   for(l = lt; l >= 1; l--) {
      r[l] = (double ***) malloc(m3[l] * sizeof(double **));
      #pragma omp parallel for default(shared) private(k, j) firstprivate(l, m3, m2, m1)
      for(k = 0; k < m3[l]; k++) {
         r[l][k] = (double **) malloc(m2[l] * sizeof(double *));
         #pragma omp parallel for default(shared) private(j) firstprivate(l, k, m2, m1)
         for(j = 0; j < m2[l]; j++) {
            r[l][k][j] = (double *) malloc(m1[l] * sizeof(double));
         }
      }
   }
   zero3(u[lt], n1, n2, n3);
   zran3(v, n1, n2, n3, nx[lt], ny[lt], lt);
   norm2u3(v, n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt], nz[lt]);
   printf(" Size: %3dx%3dx%3d (class %1c)\n", nx[lt], ny[lt], nz[lt], Class);
   printf(" Iterations: %3d\n", nit);
   resid(u[lt], v, r[lt], n1, n2, n3, a, lt);
   norm2u3(r[lt], n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt], nz[lt]);
   mg3P(u, v, r, a, c, n1, n2, n3, lt);
   resid(u[lt], v, r[lt], n1, n2, n3, a, lt);
   setup(&n1, &n2, &n3, lt);
   zero3(u[lt], n1, n2, n3);
   zran3(v, n1, n2, n3, nx[lt], ny[lt], lt);
   timer_stop(2);
   timer_start(1);
   resid(u[lt], v, r[lt], n1, n2, n3, a, lt);
   norm2u3(r[lt], n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt], nz[lt]);
   /*************** Clava msgError **************
   		Variables Access as passed arguments Can not be traced inside of function calls : 
   mg3P#362{mg3P(u, v, r, a, c, n1, n2, n3, lt)}
   resid#363{resid(u[lt], v, r[lt], n1, n2, n3, a, lt)}
   ****************************************/
   for(it = 1; it <= nit; it++) {
      mg3P(u, v, r, a, c, n1, n2, n3, lt);
      resid(u[lt], v, r[lt], n1, n2, n3, a, lt);
   }
   norm2u3(r[lt], n1, n2, n3, &rnm2, &rnmu, nx[lt], ny[lt], nz[lt]);
   {
   }
   timer_stop(1);
   t = timer_read(1);
   tinit = timer_read(2);
   verified = 0;
   verify_value = 0.0;
   printf(" Initialization time: %15.3f seconds\n", tinit);
   printf(" Benchmark completed\n");
   if(Class != 'U') {
      if(Class == 'S') {
         verify_value = 0.530770700573e-04;
      }
      else if(Class == 'W') {
            verify_value = 0.250391406439e-17;
         }
         else if(Class == 'A') {
               verify_value = 0.2433365309e-5;
            }
            else if(Class == 'B') {
                  verify_value = 0.180056440132e-5;
               }
               else if(Class == 'C') {
                     verify_value = 0.570674826298e-06;
                  }
      if(fabs(rnm2 - verify_value) <= epsilon) {
         verified = 1;
         printf(" VERIFICATION SUCCESSFUL\n");
         printf(" L2 Norm is %20.12e\n", rnm2);
         printf(" Error is   %20.12e\n", rnm2 - verify_value);
      }
      else {
         verified = 0;
         printf(" VERIFICATION FAILED\n");
         printf(" L2 Norm is             %20.12e\n", rnm2);
         printf(" The correct L2 Norm is %20.12e\n", verify_value);
      }
   }
   else {
      verified = 0;
      printf(" Problem size unknown\n");
      printf(" NO VERIFICATION PERFORMED\n");
   }
   if(t != 0.0) {
      int nn = nx[lt] * ny[lt] * nz[lt];
      mflops = 58. * nit * nn * 1.0e-6 / t;
   }
   else {
      mflops = 0.0;
   }
   c_print_results("MG", Class, nx[lt], ny[lt], nz[lt], nit, nthreads, t, mflops, "          floating point", verified, "3.0 structured", "25 Nov 2023", "gcc-12", "gcc-12", "-fopenmp -lm", "-I../common -fopenmp", "-lm", "(none)", "randdp");
}

static void setup(int *n1, int *n2, int *n3, int lt) {
   int k;
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess nx	 use : RW
   		unsolved dependency for arrayAccess ny	 use : RW
   		unsolved dependency for arrayAccess nz	 use : RW
   ****************************************/
   for(k = lt - 1; k >= 1; k--) {
      nx[k] = nx[k + 1] / 2;
      ny[k] = ny[k + 1] / 2;
      nz[k] = nz[k + 1] / 2;
   }
   #pragma omp parallel for default(shared) private(k) firstprivate(lt, nx, nz, ny)
   for(k = 1; k <= lt; k++) {
      m1[k] = nx[k] + 2;
      m2[k] = nz[k] + 2;
      m3[k] = ny[k] + 2;
   }
   is1 = 1;
   ie1 = nx[lt];
   *n1 = nx[lt] + 2;
   is2 = 1;
   ie2 = ny[lt];
   *n2 = ny[lt] + 2;
   is3 = 1;
   ie3 = nz[lt];
   *n3 = nz[lt] + 2;
   if(debug_vec[1] >= 1) {
      printf(" in setup, \n");
      printf("  lt  nx  ny  nz  n1  n2  n3 is1 is2 is3 ie1 ie2 ie3\n");
      printf("%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n", lt, nx[lt], ny[lt], nz[lt], *n1, *n2, *n3, is1, is2, is3, ie1, ie2, ie3);
   }
}

static void mg3P(double ****u, double ***v, double ****r, double a[4], double c[4], int n1, int n2, int n3, int k) {
   int j;
   /*************** Clava msgError **************
   		Variables Access as passed arguments Can not be traced inside of function calls : 
   rprj3#465{rprj3(r[k], m1[k], m2[k], m3[k], r[k - 1], m1[k - 1], m2[k - 1], m3[k - 1], k)}
   ****************************************/
   for(k = lt; k >= lb + 1; k--) {
      j = k - 1;
      rprj3(r[k], m1[k], m2[k], m3[k], r[j], m1[j], m2[j], m3[j], k);
   }
   k = lb;
   zero3(u[k], m1[k], m2[k], m3[k]);
   psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k);
   /*************** Clava msgError **************
   		Variables Access as passed arguments Can not be traced inside of function calls : 
   interp#489{interp(u[k - 1], m1[k - 1], m2[k - 1], m3[k - 1], u[k], m1[k], m2[k], m3[k], k)}
   resid#490{resid(u[k], r[k], r[k], m1[k], m2[k], m3[k], a, k)}
   psinv#491{psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k)}
   ****************************************/
   for(k = lb + 1; k <= lt - 1; k++) {
      j = k - 1;
      zero3(u[k], m1[k], m2[k], m3[k]);
      interp(u[j], m1[j], m2[j], m3[j], u[k], m1[k], m2[k], m3[k], k);
      resid(u[k], r[k], r[k], m1[k], m2[k], m3[k], a, k);
      psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k);
   }
   j = lt - 1;
   k = lt;
   interp(u[j], m1[j], m2[j], m3[j], u[lt], n1, n2, n3, k);
   resid(u[lt], v, r[lt], n1, n2, n3, a, k);
   psinv(r[lt], u[lt], n1, n2, n3, c, k);
}

static void psinv(double ***r, double ***u, int n1, int n2, int n3, double c[4], int k) {
   int i3, i2, i1;
   double r1[1037];
   double r2[1037];
   #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(n3, n2, n1, r, c, r1, r2)
   for(i3 = 1; i3 < n3 - 1; i3++) {
      #pragma omp parallel for default(shared) private(i2, i1) firstprivate(n2, n1, i3, r, c, r1, r2)
      for(i2 = 1; i2 < n2 - 1; i2++) {
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, i2, i3, r)
         for(i1 = 0; i1 < n1; i1++) {
            r1[i1] = r[i3][i2 - 1][i1] + r[i3][i2 + 1][i1] + r[i3 - 1][i2][i1] + r[i3 + 1][i2][i1];
            r2[i1] = r[i3 - 1][i2 - 1][i1] + r[i3 - 1][i2 + 1][i1] + r[i3 + 1][i2 - 1][i1] + r[i3 + 1][i2 + 1][i1];
         }
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, i3, i2, r, r1, r2, c)
         for(i1 = 1; i1 < n1 - 1; i1++) {
            u[i3][i2][i1] = u[i3][i2][i1] + c[0] * r[i3][i2][i1] + c[1] * (r[i3][i2][i1 - 1] + r[i3][i2][i1 + 1] + r1[i1]) + c[2] * (r2[i1] + r1[i1 - 1] + r1[i1 + 1]);
         }
      }
   }
   comm3(u, n1, n2, n3, k);
   if(debug_vec[0] >= 1) {
      rep_nrm(u, n1, n2, n3, "   psinv", k);
   }
   if(debug_vec[3] >= k) {
      showall(u, n1, n2, n3);
   }
}

static void resid(double ***u, double ***v, double ***r, int n1, int n2, int n3, double a[4], int k) {
   int i3, i2, i1;
   double u1[1037];
   double u2[1037];
   #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(n3, n2, n1, u, a, v, u1, u2)
   for(i3 = 1; i3 < n3 - 1; i3++) {
      #pragma omp parallel for default(shared) private(i2, i1) firstprivate(n2, n1, i3, u, a, v, u1, u2)
      for(i2 = 1; i2 < n2 - 1; i2++) {
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, i2, i3, u)
         for(i1 = 0; i1 < n1; i1++) {
            u1[i1] = u[i3][i2 - 1][i1] + u[i3][i2 + 1][i1] + u[i3 - 1][i2][i1] + u[i3 + 1][i2][i1];
            u2[i1] = u[i3 - 1][i2 - 1][i1] + u[i3 - 1][i2 + 1][i1] + u[i3 + 1][i2 - 1][i1] + u[i3 + 1][i2 + 1][i1];
         }
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, i3, i2, u2, u1, a, u, v)
         for(i1 = 1; i1 < n1 - 1; i1++) {
            r[i3][i2][i1] = v[i3][i2][i1] - a[0] * u[i3][i2][i1] - a[2] * (u2[i1] + u1[i1 - 1] + u1[i1 + 1]) - a[3] * (u2[i1 - 1] + u2[i1 + 1]);
         }
      }
   }
   comm3(r, n1, n2, n3, k);
   if(debug_vec[0] >= 1) {
      rep_nrm(r, n1, n2, n3, "   resid", k);
   }
   if(debug_vec[2] >= k) {
      showall(r, n1, n2, n3);
   }
}

static void rprj3(double ***r, int m1k, int m2k, int m3k, double ***s, int m1j, int m2j, int m3j, int k) {
   int j3, j2, j1, i3, i2, i1, d1, d2, d3;
   double x1[1037];
   double y1[1037];
   double x2;
   double y2;
   if(m1k == 3) {
      d1 = 2;
   }
   else {
      d1 = 1;
   }
   if(m2k == 3) {
      d2 = 2;
   }
   else {
      d2 = 1;
   }
   if(m3k == 3) {
      d3 = 2;
   }
   else {
      d3 = 1;
   }
   #pragma omp parallel for default(shared) private(j3, j2, j1, i3, i2, i1, y2, x2) firstprivate(m3j, d3, m2j, d2, m1j, d1, r, x1, y1)
   for(j3 = 1; j3 < m3j - 1; j3++) {
      i3 = 2 * j3 - d3;
      #pragma omp parallel for default(shared) private(j2, j1, i2, i1, y2, x2) firstprivate(m2j, d2, m1j, d1, i3, j3, d3, r, x1, y1)
      for(j2 = 1; j2 < m2j - 1; j2++) {
         i2 = 2 * j2 - d2;
         #pragma omp parallel for default(shared) private(j1, i1) firstprivate(m1j, d1, j2, i3, d2, i2, j3, d3, r)
         for(j1 = 1; j1 < m1j; j1++) {
            i1 = 2 * j1 - d1;
            x1[i1] = r[i3 + 1][i2][i1] + r[i3 + 1][i2 + 2][i1] + r[i3][i2 + 1][i1] + r[i3 + 2][i2 + 1][i1];
            y1[i1] = r[i3][i2][i1] + r[i3 + 2][i2][i1] + r[i3][i2 + 2][i1] + r[i3 + 2][i2 + 2][i1];
         }
         #pragma omp parallel for default(shared) private(j1, i1, y2, x2) firstprivate(m1j, d1, j3, j2, d3, d2, i3, i2, r, x1, y1)
         for(j1 = 1; j1 < m1j - 1; j1++) {
            i1 = 2 * j1 - d1;
            y2 = r[i3][i2][i1 + 1] + r[i3 + 2][i2][i1 + 1] + r[i3][i2 + 2][i1 + 1] + r[i3 + 2][i2 + 2][i1 + 1];
            x2 = r[i3 + 1][i2][i1 + 1] + r[i3 + 1][i2 + 2][i1 + 1] + r[i3][i2 + 1][i1 + 1] + r[i3 + 2][i2 + 1][i1 + 1];
            s[j3][j2][j1] = 0.5 * r[i3 + 1][i2 + 1][i1 + 1] + 0.25 * (r[i3 + 1][i2 + 1][i1] + r[i3 + 1][i2 + 1][i1 + 2] + x2) + 0.125 * (x1[i1] + x1[i1 + 2] + y2) + 0.0625 * (y1[i1] + y1[i1 + 2]);
         }
      }
   }
   comm3(s, m1j, m2j, m3j, k - 1);
   if(debug_vec[0] >= 1) {
      rep_nrm(s, m1j, m2j, m3j, "   rprj3", k - 1);
   }
   if(debug_vec[4] >= k) {
      showall(s, m1j, m2j, m3j);
   }
}

static void interp(double ***z, int mm1, int mm2, int mm3, double ***u, int n1, int n2, int n3, int k) {
   int i3, i2, i1, d1, d2, d3, t1, t2, t3;
   double z1[1037];
   double z2[1037];
   double z3[1037];
   if(n1 != 3 && n2 != 3 && n3 != 3) {
      #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(mm3, mm2, mm1, z, z1, z2, z3)
      for(i3 = 0; i3 < mm3 - 1; i3++) {
         #pragma omp parallel for default(shared) private(i2, i1) firstprivate(mm2, mm1, i3, z, z1, z2, z3)
         for(i2 = 0; i2 < mm2 - 1; i2++) {
            #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i2, i3, z)
            for(i1 = 0; i1 < mm1; i1++) {
               z1[i1] = z[i3][i2 + 1][i1] + z[i3][i2][i1];
               z2[i1] = z[i3 + 1][i2][i1] + z[i3][i2][i1];
               z3[i1] = z[i3 + 1][i2 + 1][i1] + z[i3 + 1][i2][i1] + z1[i1];
            }
            #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i3, i2, z)
            for(i1 = 0; i1 < mm1 - 1; i1++) {
               u[2 * i3][2 * i2][2 * i1] = u[2 * i3][2 * i2][2 * i1] + z[i3][i2][i1];
               u[2 * i3][2 * i2][2 * i1 + 1] = u[2 * i3][2 * i2][2 * i1 + 1] + 0.5 * (z[i3][i2][i1 + 1] + z[i3][i2][i1]);
            }
            #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i2, i3, z1)
            for(i1 = 0; i1 < mm1 - 1; i1++) {
               u[2 * i3][2 * i2 + 1][2 * i1] = u[2 * i3][2 * i2 + 1][2 * i1] + 0.5 * z1[i1];
               u[2 * i3][2 * i2 + 1][2 * i1 + 1] = u[2 * i3][2 * i2 + 1][2 * i1 + 1] + 0.25 * (z1[i1] + z1[i1 + 1]);
            }
            #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i3, i2, z2)
            for(i1 = 0; i1 < mm1 - 1; i1++) {
               u[2 * i3 + 1][2 * i2][2 * i1] = u[2 * i3 + 1][2 * i2][2 * i1] + 0.5 * z2[i1];
               u[2 * i3 + 1][2 * i2][2 * i1 + 1] = u[2 * i3 + 1][2 * i2][2 * i1 + 1] + 0.25 * (z2[i1] + z2[i1 + 1]);
            }
            #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i3, i2, z3)
            for(i1 = 0; i1 < mm1 - 1; i1++) {
               u[2 * i3 + 1][2 * i2 + 1][2 * i1] = u[2 * i3 + 1][2 * i2 + 1][2 * i1] + 0.25 * z3[i1];
               u[2 * i3 + 1][2 * i2 + 1][2 * i1 + 1] = u[2 * i3 + 1][2 * i2 + 1][2 * i1 + 1] + 0.125 * (z3[i1] + z3[i1 + 1]);
            }
         }
      }
   }
   else {
      if(n1 == 3) {
         d1 = 2;
         t1 = 1;
      }
      else {
         d1 = 1;
         t1 = 0;
      }
      if(n2 == 3) {
         d2 = 2;
         t2 = 1;
      }
      else {
         d2 = 1;
         t2 = 0;
      }
      if(n3 == 3) {
         d3 = 2;
         t3 = 1;
      }
      else {
         d3 = 1;
         t3 = 0;
      }
      {
         #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(d3, mm3, d2, mm2, d1, mm1, t1, t2, z)
         for(i3 = d3; i3 <= mm3 - 1; i3++) {
            #pragma omp parallel for default(shared) private(i2, i1) firstprivate(d2, mm2, d1, mm1, i3, d3, t1, z)
            for(i2 = d2; i2 <= mm2 - 1; i2++) {
               #pragma omp parallel for default(shared) private(i1) firstprivate(d1, mm1, i3, i2, d3, d2, z)
               for(i1 = d1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - d3 - 1][2 * i2 - d2 - 1][2 * i1 - d1 - 1] = u[2 * i3 - d3 - 1][2 * i2 - d2 - 1][2 * i1 - d1 - 1] + z[i3 - 1][i2 - 1][i1 - 1];
               }
               #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i3, i2, d3, d2, t1, z)
               for(i1 = 1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - d3 - 1][2 * i2 - d2 - 1][2 * i1 - t1 - 1] = u[2 * i3 - d3 - 1][2 * i2 - d2 - 1][2 * i1 - t1 - 1] + 0.5 * (z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2 - 1][i1 - 1]);
               }
            }
            #pragma omp parallel for default(shared) private(i2, i1) firstprivate(mm2, d1, mm1, i3, d3, t2, t1, z)
            for(i2 = 1; i2 <= mm2 - 1; i2++) {
               #pragma omp parallel for default(shared) private(i1) firstprivate(d1, mm1, i3, i2, d3, t2, z)
               for(i1 = d1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - d3 - 1][2 * i2 - t2 - 1][2 * i1 - d1 - 1] = u[2 * i3 - d3 - 1][2 * i2 - t2 - 1][2 * i1 - d1 - 1] + 0.5 * (z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
               }
               #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i3, i2, d3, t2, t1, z)
               for(i1 = 1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - d3 - 1][2 * i2 - t2 - 1][2 * i1 - t1 - 1] = u[2 * i3 - d3 - 1][2 * i2 - t2 - 1][2 * i1 - t1 - 1] + 0.25 * (z[i3 - 1][i2][i1] + z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
               }
            }
         }
         #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(mm3, d2, mm2, d1, mm1, t3, t1, t2, z)
         for(i3 = 1; i3 <= mm3 - 1; i3++) {
            #pragma omp parallel for default(shared) private(i2, i1) firstprivate(d2, mm2, d1, mm1, i3, t3, t1, z)
            for(i2 = d2; i2 <= mm2 - 1; i2++) {
               #pragma omp parallel for default(shared) private(i1) firstprivate(d1, mm1, i2, i3, t3, d2, z)
               for(i1 = d1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - t3 - 1][2 * i2 - d2 - 1][2 * i1 - d1 - 1] = u[2 * i3 - t3 - 1][2 * i2 - d2 - 1][2 * i1 - d1 - 1] + 0.5 * (z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
               }
               #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i2, i3, t3, d2, t1, z)
               for(i1 = 1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - t3 - 1][2 * i2 - d2 - 1][2 * i1 - t1 - 1] = u[2 * i3 - t3 - 1][2 * i2 - d2 - 1][2 * i1 - t1 - 1] + 0.25 * (z[i3][i2 - 1][i1] + z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2 - 1][i1 - 1]);
               }
            }
            #pragma omp parallel for default(shared) private(i2, i1) firstprivate(mm2, d1, mm1, i3, t3, t2, t1, z)
            for(i2 = 1; i2 <= mm2 - 1; i2++) {
               #pragma omp parallel for default(shared) private(i1) firstprivate(d1, mm1, i2, i3, t3, t2, z)
               for(i1 = d1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - t3 - 1][2 * i2 - t2 - 1][2 * i1 - d1 - 1] = u[2 * i3 - t3 - 1][2 * i2 - t2 - 1][2 * i1 - d1 - 1] + 0.25 * (z[i3][i2][i1 - 1] + z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
               }
               #pragma omp parallel for default(shared) private(i1) firstprivate(mm1, i2, i3, t3, t2, t1, z)
               for(i1 = 1; i1 <= mm1 - 1; i1++) {
                  u[2 * i3 - t3 - 1][2 * i2 - t2 - 1][2 * i1 - t1 - 1] = u[2 * i3 - t3 - 1][2 * i2 - t2 - 1][2 * i1 - t1 - 1] + 0.125 * (z[i3][i2][i1] + z[i3][i2 - 1][i1] + z[i3][i2][i1 - 1] + z[i3][i2 - 1][i1 - 1] + z[i3 - 1][i2][i1] + z[i3 - 1][i2 - 1][i1] + z[i3 - 1][i2][i1 - 1] + z[i3 - 1][i2 - 1][i1 - 1]);
               }
            }
         }
      }
   }
   if(debug_vec[0] >= 1) {
      rep_nrm(z, mm1, mm2, mm3, "z: inter", k - 1);
      rep_nrm(u, n1, n2, n3, "u: inter", k);
   }
   if(debug_vec[5] >= k) {
      showall(z, mm1, mm2, mm3);
      showall(u, n1, n2, n3);
   }
}

static void norm2u3(double ***r, int n1, int n2, int n3, double *rnm2, double *rnmu, int nx, int ny, int nz) {
   double s = 0.0;
   int i3, i2, i1, n;
   double a = 0.0, tmp = 0.0;
   n = nx * ny * nz;
   /*************** Clava msgError **************
   		Variable tmp could not be categorized into any OpenMP Variable Scopeuse : RW
   ****************************************/
   for(i3 = 1; i3 < n3 - 1; i3++) {
      /*************** Clava msgError **************
      		Variable tmp could not be categorized into any OpenMP Variable Scopeuse : RW
      ****************************************/
      for(i2 = 1; i2 < n2 - 1; i2++) {
         /*************** Clava msgError **************
         		Variable tmp could not be categorized into any OpenMP Variable Scopeuse : RW
         ****************************************/
         for(i1 = 1; i1 < n1 - 1; i1++) {
            s = s + r[i3][i2][i1] * r[i3][i2][i1];
            a = fabs(r[i3][i2][i1]);
            if(a > tmp) tmp = a;
         }
      }
   }
   *rnmu = tmp;
   *rnm2 = sqrt(s / (double) n);
}

static void rep_nrm(double ***u, int n1, int n2, int n3, char *title, int kk) {
   double rnm2, rnmu;
   norm2u3(u, n1, n2, n3, &rnm2, &rnmu, nx[kk], ny[kk], nz[kk]);
   printf(" Level%2d in %8s: norms =%21.14e%21.14e\n", kk, title, rnm2, rnmu);
}

static void comm3(double ***u, int n1, int n2, int n3, int kk) {
   int i1, i2, i3;
   {
      #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(n3, n2, n1)
      for(i3 = 1; i3 < n3 - 1; i3++) {
         #pragma omp parallel for default(shared) private(i2) firstprivate(n2, n1, i3)
         for(i2 = 1; i2 < n2 - 1; i2++) {
            u[i3][i2][n1 - 1] = u[i3][i2][1];
            u[i3][i2][0] = u[i3][i2][n1 - 2];
         }
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, n2, i3)
         for(i1 = 0; i1 < n1; i1++) {
            u[i3][n2 - 1][i1] = u[i3][1][i1];
            u[i3][0][i1] = u[i3][n2 - 2][i1];
         }
      }
      #pragma omp parallel for default(shared) private(i2, i1) firstprivate(n2, n1, n3)
      for(i2 = 0; i2 < n2; i2++) {
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, n3, i2)
         for(i1 = 0; i1 < n1; i1++) {
            u[n3 - 1][i2][i1] = u[1][i2][i1];
            u[0][i2][i1] = u[n3 - 2][i2][i1];
         }
      }
   }
}

static void zran3(double ***z, int n1, int n2, int n3, int nx, int ny, int k) {
   int i0, m0, m1;
   int i1, i2, i3, d1, e1, e2, e3;
   double xx, x0, x1, a1, a2, ai;
   double ten[10][2];
   double best;
   int i;
   int j1[10][2];
   int j2[10][2];
   int j3[10][2];
   int jg[4][10][2];
   double rdummy;
   a1 = power(pow(5.0, 13), nx);
   a2 = power(pow(5.0, 13), nx * ny);
   zero3(z, n1, n2, n3);
   i = is1 - 1 + nx * (is2 - 1 + ny * (is3 - 1));
   ai = power(pow(5.0, 13), i);
   d1 = ie1 - is1 + 1;
   e1 = ie1 - is1 + 2;
   e2 = ie2 - is2 + 2;
   e3 = ie3 - is3 + 2;
   x0 = 314159265.e0;
   rdummy = randlc(&x0, ai);
   /*************** Clava msgError **************
   		Variables Access as passed arguments Can not be traced inside of function calls : 
   vranlc#865{vranlc(d1, &xx, pow(5.0, 13), &(z[i3][i2][0]))}
   ****************************************/
   for(i3 = 1; i3 < e3; i3++) {
      x1 = x0;
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      vranlc#865{vranlc(d1, &xx, pow(5.0, 13), &(z[i3][i2][0]))}
      ****************************************/
      for(i2 = 1; i2 < e2; i2++) {
         xx = x1;
         vranlc(d1, &xx, pow(5.0, 13), &(z[i3][i2][0]));
         rdummy = randlc(&x1, a1);
      }
      rdummy = randlc(&x0, a2);
   }
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(i = 0; i < 10; i++) {
      ten[i][1] = 0.0;
      j1[i][1] = 0;
      j2[i][1] = 0;
      j3[i][1] = 0;
      ten[i][0] = 1.0;
      j1[i][0] = 0;
      j2[i][0] = 0;
      j3[i][0] = 0;
   }
   /*************** Clava msgError **************
   		Variables Access as passed arguments Can not be traced inside of function calls : 
   bubble#926{bubble(ten, j1, j2, j3, 10, 1)}
   bubble#933{bubble(ten, j1, j2, j3, 10, 0)}
   ****************************************/
   for(i3 = 1; i3 < n3 - 1; i3++) {
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      bubble#926{bubble(ten, j1, j2, j3, 10, 1)}
      bubble#933{bubble(ten, j1, j2, j3, 10, 0)}
      ****************************************/
      for(i2 = 1; i2 < n2 - 1; i2++) {
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         bubble#926{bubble(ten, j1, j2, j3, 10, 1)}
         bubble#933{bubble(ten, j1, j2, j3, 10, 0)}
         ****************************************/
         for(i1 = 1; i1 < n1 - 1; i1++) {
            if(z[i3][i2][i1] > ten[0][1]) {
               ten[0][1] = z[i3][i2][i1];
               j1[0][1] = i1;
               j2[0][1] = i2;
               j3[0][1] = i3;
               bubble(ten, j1, j2, j3, 10, 1);
            }
            if(z[i3][i2][i1] < ten[0][0]) {
               ten[0][0] = z[i3][i2][i1];
               j1[0][0] = i1;
               j2[0][0] = i2;
               j3[0][0] = i3;
               bubble(ten, j1, j2, j3, 10, 0);
            }
         }
      }
   }
   i1 = 10 - 1;
   i0 = 10 - 1;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(i = 10 - 1; i >= 0; i--) {
      best = z[j3[i1][1]][j2[i1][1]][j1[i1][1]];
      if(best == z[j3[i1][1]][j2[i1][1]][j1[i1][1]]) {
         jg[0][i][1] = 0;
         jg[1][i][1] = is1 - 1 + j1[i1][1];
         jg[2][i][1] = is2 - 1 + j2[i1][1];
         jg[3][i][1] = is3 - 1 + j3[i1][1];
         i1 = i1 - 1;
      }
      else {
         jg[0][i][1] = 0;
         jg[1][i][1] = 0;
         jg[2][i][1] = 0;
         jg[3][i][1] = 0;
      }
      ten[i][1] = best;
      best = z[j3[i0][0]][j2[i0][0]][j1[i0][0]];
      if(best == z[j3[i0][0]][j2[i0][0]][j1[i0][0]]) {
         jg[0][i][0] = 0;
         jg[1][i][0] = is1 - 1 + j1[i0][0];
         jg[2][i][0] = is2 - 1 + j2[i0][0];
         jg[3][i][0] = is3 - 1 + j3[i0][0];
         i0 = i0 - 1;
      }
      else {
         jg[0][i][0] = 0;
         jg[1][i][0] = 0;
         jg[2][i][0] = 0;
         jg[3][i][0] = 0;
      }
      ten[i][0] = best;
   }
   m1 = i1 + 1;
   m0 = i0 + 1;
   #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(n3, n2, n1)
   for(i3 = 0; i3 < n3; i3++) {
      #pragma omp parallel for default(shared) private(i2, i1) firstprivate(n2, n1, i3)
      for(i2 = 0; i2 < n2; i2++) {
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, i3, i2)
         for(i1 = 0; i1 < n1; i1++) {
            z[i3][i2][i1] = 0.0;
         }
      }
   }
   /*************** Clava msgError **************
   		 Array access z[j3[i][0]][j2[i][0]][j1[i][0]] which is used for writing has subscript of arrayType j3[i][0]
   ****************************************/
   for(i = 10 - 1; i >= m0; i--) {
      z[j3[i][0]][j2[i][0]][j1[i][0]] = -1.0;
   }
   /*************** Clava msgError **************
   		 Array access z[j3[i][1]][j2[i][1]][j1[i][1]] which is used for writing has subscript of arrayType j3[i][1]
   ****************************************/
   for(i = 10 - 1; i >= m1; i--) {
      z[j3[i][1]][j2[i][1]][j1[i][1]] = 1.0;
   }
   comm3(z, n1, n2, n3, k);
}

static void showall(double ***z, int n1, int n2, int n3) {
   int i1, i2, i3;
   int m1, m2, m3;
   m1 = (((n1) < (18)) ? (n1) : (18));
   m2 = (((n2) < (14)) ? (n2) : (14));
   m3 = (((n3) < (18)) ? (n3) : (18));
   printf("\n");
   /*************** Clava msgError **************
   		Variables Access as passed arguments Can not be traced inside of function calls : 
   printf#1018{printf("%6.3f", z[i3][i2][i1])}
   printf#1020{printf("\n")}
   printf#1022{printf(" - - - - - - - \n")}
   ****************************************/
   for(i3 = 0; i3 < m3; i3++) {
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      printf#1018{printf("%6.3f", z[i3][i2][i1])}
      printf#1020{printf("\n")}
      ****************************************/
      for(i1 = 0; i1 < m1; i1++) {
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         printf#1018{printf("%6.3f", z[i3][i2][i1])}
         ****************************************/
         for(i2 = 0; i2 < m2; i2++) {
            printf("%6.3f", z[i3][i2][i1]);
         }
         printf("\n");
      }
      printf(" - - - - - - - \n");
   }
   printf("\n");
}

static double power(double a, int n) {
   double aj;
   int nj;
   double rdummy;
   double power;
   power = 1.0;
   nj = n;
   aj = a;
   while(nj != 0) {
      if((nj % 2) == 1) rdummy = randlc(&power, aj);
      rdummy = randlc(&aj, aj);
      nj = nj / 2;
   }
   
   return (power);
}

static void bubble(double ten[1037][2], int j1[1037][2], int j2[1037][2], int j3[1037][2], int m, int ind) {
   double temp;
   int i, j_temp;
   if(ind == 1) {
      /*************** Clava msgError **************
      		Loop contains Invalid Statement -> ReturnStmt#1069
      ****************************************/
      for(i = 0; i < m - 1; i++) {
         if(ten[i][ind] > ten[i + 1][ind]) {
            temp = ten[i + 1][ind];
            ten[i + 1][ind] = ten[i][ind];
            ten[i][ind] = temp;
            j_temp = j1[i + 1][ind];
            j1[i + 1][ind] = j1[i][ind];
            j1[i][ind] = j_temp;
            j_temp = j2[i + 1][ind];
            j2[i + 1][ind] = j2[i][ind];
            j2[i][ind] = j_temp;
            j_temp = j3[i + 1][ind];
            j3[i + 1][ind] = j3[i][ind];
            j3[i][ind] = j_temp;
         }
         else {
            
            return;
         }
      }
   }
   else {
      /*************** Clava msgError **************
      		Loop contains Invalid Statement -> ReturnStmt#1093
      ****************************************/
      for(i = 0; i < m - 1; i++) {
         if(ten[i][ind] < ten[i + 1][ind]) {
            temp = ten[i + 1][ind];
            ten[i + 1][ind] = ten[i][ind];
            ten[i][ind] = temp;
            j_temp = j1[i + 1][ind];
            j1[i + 1][ind] = j1[i][ind];
            j1[i][ind] = j_temp;
            j_temp = j2[i + 1][ind];
            j2[i + 1][ind] = j2[i][ind];
            j2[i][ind] = j_temp;
            j_temp = j3[i + 1][ind];
            j3[i + 1][ind] = j3[i][ind];
            j3[i][ind] = j_temp;
         }
         else {
            
            return;
         }
      }
   }
}

static void zero3(double ***z, int n1, int n2, int n3) {
   int i1, i2, i3;
   #pragma omp parallel for default(shared) private(i3, i2, i1) firstprivate(n3, n2, n1)
   for(i3 = 0; i3 < n3; i3++) {
      #pragma omp parallel for default(shared) private(i2, i1) firstprivate(n2, n1, i3)
      for(i2 = 0; i2 < n2; i2++) {
         #pragma omp parallel for default(shared) private(i1) firstprivate(n1, i3, i2)
         for(i1 = 0; i1 < n1; i1++) {
            z[i3][i2][i1] = 0.0;
         }
      }
   }
}
