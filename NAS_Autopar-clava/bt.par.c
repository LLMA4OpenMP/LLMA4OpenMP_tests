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

/*COMMON block: global*/

static int grid_points[3];
/*grid_ponts(1:3)*/

/*COMMON block: constants*/

static double tx1;
static double tx2;
static double tx3;
static double ty1;
static double ty2;
static double ty3;
static double tz1;
static double tz2;
static double tz3;
static double dx1;
static double dx2;
static double dx3;
static double dx4;
static double dx5;
static double dy1;
static double dy2;
static double dy3;
static double dy4;
static double dy5;
static double dz1;
static double dz2;
static double dz3;
static double dz4;
static double dz5;
static double dssp;
static double dt;
static double ce[5][13];
/*ce(5,13)*/

static double dxmax;
static double dymax;
static double dzmax;
static double xxcon1;
static double xxcon2;
static double xxcon3;
static double xxcon4;
static double xxcon5;
static double dx1tx1;
static double dx2tx1;
static double dx3tx1;
static double dx4tx1;
static double dx5tx1;
static double yycon1;
static double yycon2;
static double yycon3;
static double yycon4;
static double yycon5;
static double dy1ty1;
static double dy2ty1;
static double dy3ty1;
static double dy4ty1;
static double dy5ty1;
static double zzcon1;
static double zzcon2;
static double zzcon3;
static double zzcon4;
static double zzcon5;
static double dz1tz1;
static double dz2tz1;
static double dz3tz1;
static double dz4tz1;
static double dz5tz1;
static double dnxm1;
static double dnym1;
static double dnzm1;
static double c1c2;
static double c1c5;
static double c3c4;
static double c1345;
static double conz1;
static double c1;
static double c2;
static double c3;
static double c4;
static double c5;
static double c4dssp;
static double c5dssp;
static double dtdssp;
static double dttx1;
static double dttx2;
static double dtty1;
static double dtty2;
static double dttz1;
static double dttz2;
static double c2dttx1;
static double c2dtty1;
static double c2dttz1;
static double comz1;
static double comz4;
static double comz5;
static double comz6;
static double c3c4tx3;
static double c3c4ty3;
static double c3c4tz3;
static double c2iv;
static double con43;
static double con16;
/*
c   to improve cache performance, grid dimensions padded by 1
c   for even number sizes only.
*/

/*COMMON block: fields*/

static double us[65][65][65];
static double vs[65][65][65];
static double ws[65][65][65];
static double qs[65][65][65];
static double rho_i[65][65][65];
static double square[65][65][65];
static double forcing[65][65][65][6];
static double u[65][65][65][5];
static double rhs[65][65][65][5];
static double lhs[65][65][65][3][5][5];
/*COMMON block: work_1d*/

static double cuf[64];
static double q[64];
static double ue[64][5];
static double buf[64][5];
/*
c   to improve cache performance, grid dimensions (first two for these
c   to arrays) padded by 1 for even number sizes only.
*/

/*COMMON block: work_lhs*/

static double fjac[65][65][64][5][5];
/*fjac(5, 5, 0:IMAX/2*2, 0:JMAX/2*2, 0:KMAX-1)*/

static double njac[65][65][64][5][5];
/*njac(5, 5, 0:IMAX/2*2, 0:JMAX/2*2, 0:KMAX-1)*/

static double tmp1;
static double tmp2;
static double tmp3;
static void add();
static void adi();
static void error_norm(double rms[5]);
static void rhs_norm(double rms[5]);
static void exact_rhs();
static void exact_solution(double xi, double eta, double zeta, double dtemp[5]);
static void initialize();
static void lhsinit();
static void lhsx();
static void lhsy();
static void lhsz();
static void compute_rhs();
static void set_constants();
static void verify(int no_time_steps, char *class, boolean *verified);
static void x_solve();
static void x_backsubstitute();
static void x_solve_cell();
static void matvec_sub(double ablock[5][5], double avec[5], double bvec[5]);
static void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5]);
static void binvcrhs(double lhs[5][5], double c[5][5], double r[5]);
static void binvrhs(double lhs[5][5], double r[5]);
static void y_solve();
static void y_backsubstitute();
static void y_solve_cell();
static void z_solve();
static void z_backsubstitute();
static void z_solve_cell();
int main(int argc, char **argv) {
   int niter, step, n3;
   int nthreads = 1;
   double navg, mflops;
   double tmax;
   boolean verified;
   char class;
   FILE *fp;
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - BT Benchmark\n\n");
   fp = fopen("inputbt.data", "r");
   if(fp != ((void *) 0)) {
      printf(" Reading from input file inputbt.data");
      fscanf(fp, "%d", &niter);
      while(fgetc(fp) != '\n');
      fscanf(fp, "%lg", &dt);
      while(fgetc(fp) != '\n');
      fscanf(fp, "%d%d%d", &grid_points[0], &grid_points[1], &grid_points[2]);
      fclose(fp);
   }
   else {
      printf(" No input file inputbt.data. Using compiled defaults\n");
      niter = 200;
      dt = 0.0008;
      grid_points[0] = 64;
      grid_points[1] = 64;
      grid_points[2] = 64;
   }
   printf(" Size: %3dx%3dx%3d\n", grid_points[0], grid_points[1], grid_points[2]);
   printf(" Iterations: %3d   dt: %10.6f\n", niter, dt);
   if(grid_points[0] > 64 || grid_points[1] > 64 || grid_points[2] > 64) {
      printf(" %dx%dx%d\n", grid_points[0], grid_points[1], grid_points[2]);
      printf(" Problem size too big for compiled array sizes\n");
      exit(1);
   }
   set_constants();
   initialize();
   lhsinit();
   exact_rhs();
   adi();
   initialize();
   timer_clear(1);
   timer_start(1);
   /*************** Clava msgError **************
   		Variables Access as passed arguments Can not be traced inside of function calls : 
   printf#390{printf(" Time step %4d\n", step)}
   compute_rhs#393{compute_rhs()}
   x_solve#394{x_solve()}
   y_solve#395{y_solve()}
   z_solve#396{z_solve()}
   add#397{add()}
   ****************************************/
   for(step = 1; step <= niter; step++) {
      if(step % 20 == 0 || step == 1) {
         printf(" Time step %4d\n", step);
      }
      adi();
   }
   {
   }
   timer_stop(1);
   tmax = timer_read(1);
   verify(niter, &class, &verified);
   n3 = grid_points[0] * grid_points[1] * grid_points[2];
   navg = (grid_points[0] + grid_points[1] + grid_points[2]) / 3.0;
   if(tmax != 0.0) {
      mflops = 1.0e-6 * (double) niter * (3478.8 * (double) n3 - 17655.7 * ((navg) * (navg)) + 28023.7 * navg) / tmax;
   }
   else {
      mflops = 0.0;
   }
   c_print_results("BT", class, grid_points[0], grid_points[1], grid_points[2], niter, nthreads, tmax, mflops, "          floating point", verified, "3.0 structured", "25 Nov 2023", "gcc-12", "gcc-12", "-fopenmp -lm", "-I../common -fopenmp", "-lm", "(none)", "(none)");
}

static void add() {
   int i, j, k, m;
   #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(grid_points, rhs)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, grid_points, rhs)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, grid_points, rhs)
         for(k = 1; k < grid_points[2] - 1; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               u[i][j][k][m] = u[i][j][k][m] + rhs[i][j][k][m];
            }
         }
      }
   }
}

static void adi() {
   compute_rhs();
   x_solve();
   y_solve();
   z_solve();
   add();
}

static void error_norm(double rms[5]) {
   int i, j, k, m, d;
   double xi;
   double eta;
   double zeta;
   double u_exact[5];
   double add;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      rms[m] = 0.0;
   }
   #pragma omp parallel for default(shared) private(i, j, k, m, xi, eta, zeta, add) firstprivate(dnxm1, dnym1, dnzm1, grid_points, ce, u, u_exact) reduction(+ : rms[:5])
   for(i = 0; i < grid_points[0]; i++) {
      xi = (double) i * dnxm1;
      #pragma omp parallel for default(shared) private(j, k, m, eta, zeta, add) firstprivate(dnym1, dnzm1, xi, i, grid_points, ce, u, u_exact) reduction(+ : rms[:5])
      for(j = 0; j < grid_points[1]; j++) {
         eta = (double) j * dnym1;
         #pragma omp parallel for default(shared) private(k, m, zeta, add) firstprivate(dnzm1, eta, xi, i, j, grid_points, ce, u, u_exact) reduction(+ : rms[:5])
         for(k = 0; k < grid_points[2]; k++) {
            zeta = (double) k * dnzm1;
            exact_solution(xi, eta, zeta, u_exact);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               add = u[i][j][k][m] - u_exact[m];
               rms[m] = rms[m] + add * add;
            }
         }
      }
   }
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(d = 0; d <= 2; d++) {
         rms[m] = rms[m] / (double) (grid_points[d] - 2);
      }
      rms[m] = sqrt(rms[m]);
   }
}

static void rhs_norm(double rms[5]) {
   int i, j, k, d, m;
   double add;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      rms[m] = 0.0;
   }
   #pragma omp parallel for default(shared) private(i, j, k, m, add) firstprivate(grid_points, rhs) reduction(+ : rms[:5])
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, m, add) firstprivate(i, grid_points, rhs) reduction(+ : rms[:5])
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, m, add) firstprivate(i, j, grid_points, rhs) reduction(+ : rms[:5])
         for(k = 1; k < grid_points[2] - 1; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               add = rhs[i][j][k][m];
               rms[m] = rms[m] + add * add;
            }
         }
      }
   }
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(d = 0; d <= 2; d++) {
         rms[m] = rms[m] / (double) (grid_points[d] - 2);
      }
      rms[m] = sqrt(rms[m]);
   }
}

static void exact_rhs() {
   {
      double dtemp[5];
      double xi;
      double eta;
      double zeta;
      double dtpp;
      int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(grid_points)
      for(i = 0; i < grid_points[0]; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(grid_points)
         for(j = 0; j < grid_points[1]; j++) {
            #pragma omp parallel for default(shared) private(k, m) firstprivate(grid_points)
            for(k = 0; k < grid_points[2]; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  forcing[i][j][k][m] = 0.0;
               }
            }
         }
      }
      #pragma omp parallel for default(shared) private(j, k, i, m) firstprivate(dnym1, dnzm1, dnxm1, tx2, dx1tx1, c2, xxcon1, dx2tx1, xxcon2, dx3tx1, dx4tx1, c1, xxcon3, xxcon4, xxcon5, dx5tx1, dssp, grid_points, ce, dtemp, ue, buf, cuf, q)
      for(j = 1; j < grid_points[1] - 1; j++) {
         eta = (double) j * dnym1;
         #pragma omp parallel for default(shared) private(k, i, m) firstprivate(dnzm1, dnxm1, tx2, dx1tx1, c2, xxcon1, dx2tx1, xxcon2, dx3tx1, dx4tx1, c1, xxcon3, xxcon4, xxcon5, dx5tx1, dssp, grid_points, ce, dtemp, ue, buf, cuf, q)
         for(k = 1; k < grid_points[2] - 1; k++) {
            zeta = (double) k * dnzm1;
            #pragma omp parallel for default(shared) private(i, m) firstprivate(dnxm1, grid_points, ce, dtemp)
            for(i = 0; i < grid_points[0]; i++) {
               xi = (double) i * dnxm1;
               exact_solution(xi, eta, zeta, dtemp);
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  ue[i][m] = dtemp[m];
               }
               dtpp = 1.0 / dtemp[0];
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 1; m <= 4; m++) {
                  buf[i][m] = dtpp * dtemp[m];
               }
               cuf[i] = buf[i][1] * buf[i][1];
               buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] + buf[i][3] * buf[i][3];
               q[i] = 0.5 * (buf[i][1] * ue[i][1] + buf[i][2] * ue[i][2] + buf[i][3] * ue[i][3]);
            }
            #pragma omp parallel for default(shared) private(i) firstprivate(tx2, dx1tx1, c2, xxcon1, dx2tx1, xxcon2, dx3tx1, dx4tx1, c1, xxcon3, xxcon4, xxcon5, dx5tx1, grid_points, ue, q, buf, cuf)
            for(i = 1; i < grid_points[0] - 1; i++) {
               im1 = i - 1;
               ip1 = i + 1;
               forcing[i][j][k][0] = forcing[i][j][k][0] - tx2 * (ue[ip1][1] - ue[im1][1]) + dx1tx1 * (ue[ip1][0] - 2.0 * ue[i][0] + ue[im1][0]);
               forcing[i][j][k][1] = forcing[i][j][k][1] - tx2 * ((ue[ip1][1] * buf[ip1][1] + c2 * (ue[ip1][4] - q[ip1])) - (ue[im1][1] * buf[im1][1] + c2 * (ue[im1][4] - q[im1]))) + xxcon1 * (buf[ip1][1] - 2.0 * buf[i][1] + buf[im1][1]) + dx2tx1 * (ue[ip1][1] - 2.0 * ue[i][1] + ue[im1][1]);
               forcing[i][j][k][2] = forcing[i][j][k][2] - tx2 * (ue[ip1][2] * buf[ip1][1] - ue[im1][2] * buf[im1][1]) + xxcon2 * (buf[ip1][2] - 2.0 * buf[i][2] + buf[im1][2]) + dx3tx1 * (ue[ip1][2] - 2.0 * ue[i][2] + ue[im1][2]);
               forcing[i][j][k][3] = forcing[i][j][k][3] - tx2 * (ue[ip1][3] * buf[ip1][1] - ue[im1][3] * buf[im1][1]) + xxcon2 * (buf[ip1][3] - 2.0 * buf[i][3] + buf[im1][3]) + dx4tx1 * (ue[ip1][3] - 2.0 * ue[i][3] + ue[im1][3]);
               forcing[i][j][k][4] = forcing[i][j][k][4] - tx2 * (buf[ip1][1] * (c1 * ue[ip1][4] - c2 * q[ip1]) - buf[im1][1] * (c1 * ue[im1][4] - c2 * q[im1])) + 0.5 * xxcon3 * (buf[ip1][0] - 2.0 * buf[i][0] + buf[im1][0]) + xxcon4 * (cuf[ip1] - 2.0 * cuf[i] + cuf[im1]) + xxcon5 * (buf[ip1][4] - 2.0 * buf[i][4] + buf[im1][4]) + dx5tx1 * (ue[ip1][4] - 2.0 * ue[i][4] + ue[im1][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               i = 1;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (5.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
               i = 2;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (-4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               #pragma omp parallel for default(shared) private(i) firstprivate(dssp, grid_points, ue)
               for(i = 1 * 3; i <= grid_points[0] - 3 * 1 - 1; i++) {
                  forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m] + ue[i + 2][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               i = grid_points[0] - 3;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 6.0 * ue[i][m] - 4.0 * ue[i + 1][m]);
               i = grid_points[0] - 2;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[i - 2][m] - 4.0 * ue[i - 1][m] + 5.0 * ue[i][m]);
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, k, j, m) firstprivate(dnxm1, dnzm1, dnym1, ty2, dy1ty1, yycon2, dy2ty1, c2, yycon1, dy3ty1, dy4ty1, c1, yycon3, yycon4, yycon5, dy5ty1, dssp, grid_points, ce, dtemp, ue, buf, cuf, q)
      for(i = 1; i < grid_points[0] - 1; i++) {
         xi = (double) i * dnxm1;
         #pragma omp parallel for default(shared) private(k, j, m) firstprivate(dnzm1, dnym1, ty2, dy1ty1, yycon2, dy2ty1, c2, yycon1, dy3ty1, dy4ty1, c1, yycon3, yycon4, yycon5, dy5ty1, dssp, grid_points, ce, dtemp, ue, buf, cuf, q)
         for(k = 1; k < grid_points[2] - 1; k++) {
            zeta = (double) k * dnzm1;
            #pragma omp parallel for default(shared) private(j, m) firstprivate(dnym1, grid_points, ce, dtemp)
            for(j = 0; j < grid_points[1]; j++) {
               eta = (double) j * dnym1;
               exact_solution(xi, eta, zeta, dtemp);
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  ue[j][m] = dtemp[m];
               }
               dtpp = 1.0 / dtemp[0];
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 1; m <= 4; m++) {
                  buf[j][m] = dtpp * dtemp[m];
               }
               cuf[j] = buf[j][2] * buf[j][2];
               buf[j][0] = cuf[j] + buf[j][1] * buf[j][1] + buf[j][3] * buf[j][3];
               q[j] = 0.5 * (buf[j][1] * ue[j][1] + buf[j][2] * ue[j][2] + buf[j][3] * ue[j][3]);
            }
            #pragma omp parallel for default(shared) private(j) firstprivate(ty2, dy1ty1, yycon2, dy2ty1, c2, yycon1, dy3ty1, dy4ty1, c1, yycon3, yycon4, yycon5, dy5ty1, grid_points, ue, buf, q, cuf)
            for(j = 1; j < grid_points[1] - 1; j++) {
               jm1 = j - 1;
               jp1 = j + 1;
               forcing[i][j][k][0] = forcing[i][j][k][0] - ty2 * (ue[jp1][2] - ue[jm1][2]) + dy1ty1 * (ue[jp1][0] - 2.0 * ue[j][0] + ue[jm1][0]);
               forcing[i][j][k][1] = forcing[i][j][k][1] - ty2 * (ue[jp1][1] * buf[jp1][2] - ue[jm1][1] * buf[jm1][2]) + yycon2 * (buf[jp1][1] - 2.0 * buf[j][1] + buf[jm1][1]) + dy2ty1 * (ue[jp1][1] - 2.0 * ue[j][1] + ue[jm1][1]);
               forcing[i][j][k][2] = forcing[i][j][k][2] - ty2 * ((ue[jp1][2] * buf[jp1][2] + c2 * (ue[jp1][4] - q[jp1])) - (ue[jm1][2] * buf[jm1][2] + c2 * (ue[jm1][4] - q[jm1]))) + yycon1 * (buf[jp1][2] - 2.0 * buf[j][2] + buf[jm1][2]) + dy3ty1 * (ue[jp1][2] - 2.0 * ue[j][2] + ue[jm1][2]);
               forcing[i][j][k][3] = forcing[i][j][k][3] - ty2 * (ue[jp1][3] * buf[jp1][2] - ue[jm1][3] * buf[jm1][2]) + yycon2 * (buf[jp1][3] - 2.0 * buf[j][3] + buf[jm1][3]) + dy4ty1 * (ue[jp1][3] - 2.0 * ue[j][3] + ue[jm1][3]);
               forcing[i][j][k][4] = forcing[i][j][k][4] - ty2 * (buf[jp1][2] * (c1 * ue[jp1][4] - c2 * q[jp1]) - buf[jm1][2] * (c1 * ue[jm1][4] - c2 * q[jm1])) + 0.5 * yycon3 * (buf[jp1][0] - 2.0 * buf[j][0] + buf[jm1][0]) + yycon4 * (cuf[jp1] - 2.0 * cuf[j] + cuf[jm1]) + yycon5 * (buf[jp1][4] - 2.0 * buf[j][4] + buf[jm1][4]) + dy5ty1 * (ue[jp1][4] - 2.0 * ue[j][4] + ue[jm1][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               j = 1;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (5.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
               j = 2;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (-4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               #pragma omp parallel for default(shared) private(j) firstprivate(dssp, grid_points, ue)
               for(j = 1 * 3; j <= grid_points[1] - 3 * 1 - 1; j++) {
                  forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               j = grid_points[1] - 3;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m]);
               j = grid_points[1] - 2;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[j - 2][m] - 4.0 * ue[j - 1][m] + 5.0 * ue[j][m]);
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(dnxm1, dnym1, dnzm1, tz2, dz1tz1, zzcon2, dz2tz1, dz3tz1, c2, zzcon1, dz4tz1, c1, zzcon3, zzcon4, zzcon5, dz5tz1, dssp, grid_points, ce, dtemp, ue, buf, cuf, q)
      for(i = 1; i < grid_points[0] - 1; i++) {
         xi = (double) i * dnxm1;
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(dnym1, dnzm1, tz2, dz1tz1, zzcon2, dz2tz1, dz3tz1, c2, zzcon1, dz4tz1, c1, zzcon3, zzcon4, zzcon5, dz5tz1, dssp, grid_points, ce, dtemp, ue, buf, cuf, q)
         for(j = 1; j < grid_points[1] - 1; j++) {
            eta = (double) j * dnym1;
            #pragma omp parallel for default(shared) private(k, m) firstprivate(dnzm1, grid_points, ce, dtemp)
            for(k = 0; k < grid_points[2]; k++) {
               zeta = (double) k * dnzm1;
               exact_solution(xi, eta, zeta, dtemp);
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  ue[k][m] = dtemp[m];
               }
               dtpp = 1.0 / dtemp[0];
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 1; m <= 4; m++) {
                  buf[k][m] = dtpp * dtemp[m];
               }
               cuf[k] = buf[k][3] * buf[k][3];
               buf[k][0] = cuf[k] + buf[k][1] * buf[k][1] + buf[k][2] * buf[k][2];
               q[k] = 0.5 * (buf[k][1] * ue[k][1] + buf[k][2] * ue[k][2] + buf[k][3] * ue[k][3]);
            }
            #pragma omp parallel for default(shared) private(k) firstprivate(tz2, dz1tz1, zzcon2, dz2tz1, dz3tz1, c2, zzcon1, dz4tz1, c1, zzcon3, zzcon4, zzcon5, dz5tz1, grid_points, ue, buf, q, cuf)
            for(k = 1; k < grid_points[2] - 1; k++) {
               km1 = k - 1;
               kp1 = k + 1;
               forcing[i][j][k][0] = forcing[i][j][k][0] - tz2 * (ue[kp1][3] - ue[km1][3]) + dz1tz1 * (ue[kp1][0] - 2.0 * ue[k][0] + ue[km1][0]);
               forcing[i][j][k][1] = forcing[i][j][k][1] - tz2 * (ue[kp1][1] * buf[kp1][3] - ue[km1][1] * buf[km1][3]) + zzcon2 * (buf[kp1][1] - 2.0 * buf[k][1] + buf[km1][1]) + dz2tz1 * (ue[kp1][1] - 2.0 * ue[k][1] + ue[km1][1]);
               forcing[i][j][k][2] = forcing[i][j][k][2] - tz2 * (ue[kp1][2] * buf[kp1][3] - ue[km1][2] * buf[km1][3]) + zzcon2 * (buf[kp1][2] - 2.0 * buf[k][2] + buf[km1][2]) + dz3tz1 * (ue[kp1][2] - 2.0 * ue[k][2] + ue[km1][2]);
               forcing[i][j][k][3] = forcing[i][j][k][3] - tz2 * ((ue[kp1][3] * buf[kp1][3] + c2 * (ue[kp1][4] - q[kp1])) - (ue[km1][3] * buf[km1][3] + c2 * (ue[km1][4] - q[km1]))) + zzcon1 * (buf[kp1][3] - 2.0 * buf[k][3] + buf[km1][3]) + dz4tz1 * (ue[kp1][3] - 2.0 * ue[k][3] + ue[km1][3]);
               forcing[i][j][k][4] = forcing[i][j][k][4] - tz2 * (buf[kp1][3] * (c1 * ue[kp1][4] - c2 * q[kp1]) - buf[km1][3] * (c1 * ue[km1][4] - c2 * q[km1])) + 0.5 * zzcon3 * (buf[kp1][0] - 2.0 * buf[k][0] + buf[km1][0]) + zzcon4 * (cuf[kp1] - 2.0 * cuf[k] + cuf[km1]) + zzcon5 * (buf[kp1][4] - 2.0 * buf[k][4] + buf[km1][4]) + dz5tz1 * (ue[kp1][4] - 2.0 * ue[k][4] + ue[km1][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               k = 1;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (5.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
               k = 2;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (-4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               #pragma omp parallel for default(shared) private(k) firstprivate(dssp, grid_points, ue)
               for(k = 1 * 3; k <= grid_points[2] - 3 * 1 - 1; k++) {
                  forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m] + ue[k + 2][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               k = grid_points[2] - 3;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 6.0 * ue[k][m] - 4.0 * ue[k + 1][m]);
               k = grid_points[2] - 2;
               forcing[i][j][k][m] = forcing[i][j][k][m] - dssp * (ue[k - 2][m] - 4.0 * ue[k - 1][m] + 5.0 * ue[k][m]);
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(grid_points)
      for(i = 1; i < grid_points[0] - 1; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(grid_points)
         for(j = 1; j < grid_points[1] - 1; j++) {
            #pragma omp parallel for default(shared) private(k, m) firstprivate(grid_points)
            for(k = 1; k < grid_points[2] - 1; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  forcing[i][j][k][m] = -1.0 * forcing[i][j][k][m];
               }
            }
         }
      }
   }
}

static void exact_solution(double xi, double eta, double zeta, double dtemp[5]) {
   int m;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      dtemp[m] = ce[m][0] + xi * (ce[m][1] + xi * (ce[m][4] + xi * (ce[m][7] + xi * ce[m][10]))) + eta * (ce[m][2] + eta * (ce[m][5] + eta * (ce[m][8] + eta * ce[m][11]))) + zeta * (ce[m][3] + zeta * (ce[m][6] + zeta * (ce[m][9] + zeta * ce[m][12])));
   }
}

static void initialize() {
   {
      int i, j, k, m, ix, iy, iz;
      double xi;
      double eta;
      double zeta;
      double Pface[2][3][5];
      double Pxi;
      double Peta;
      double Pzeta;
      double temp[5];
      #pragma omp parallel for default(shared) private(i, j, k, m)
      for(i = 0; i < 64; i++) {
         #pragma omp parallel for default(shared) private(j, k, m)
         for(j = 0; j < 64; j++) {
            #pragma omp parallel for default(shared) private(k, m)
            for(k = 0; k < 64; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  u[i][j][k][m] = 1.0;
               }
            }
         }
      }
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      exact_solution#855{exact_solution((double) ix, eta, zeta, &(Pface[ix][0][0]))}
      exact_solution#860{exact_solution(xi, (double) iy, zeta, &Pface[iy][1][0])}
      exact_solution#865{exact_solution(xi, eta, (double) iz, &Pface[iz][2][0])}
      ****************************************/
      for(i = 0; i < grid_points[0]; i++) {
         xi = (double) i * dnxm1;
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         exact_solution#855{exact_solution((double) ix, eta, zeta, &(Pface[ix][0][0]))}
         exact_solution#860{exact_solution(xi, (double) iy, zeta, &Pface[iy][1][0])}
         exact_solution#865{exact_solution(xi, eta, (double) iz, &Pface[iz][2][0])}
         ****************************************/
         for(j = 0; j < grid_points[1]; j++) {
            eta = (double) j * dnym1;
            /*************** Clava msgError **************
            		Variables Access as passed arguments Can not be traced inside of function calls : 
            exact_solution#855{exact_solution((double) ix, eta, zeta, &(Pface[ix][0][0]))}
            exact_solution#860{exact_solution(xi, (double) iy, zeta, &Pface[iy][1][0])}
            exact_solution#865{exact_solution(xi, eta, (double) iz, &Pface[iz][2][0])}
            ****************************************/
            for(k = 0; k < grid_points[2]; k++) {
               zeta = (double) k * dnzm1;
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(ix = 0; ix < 2; ix++) {
                  exact_solution((double) ix, eta, zeta, &(Pface[ix][0][0]));
               }
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(iy = 0; iy < 2; iy++) {
                  exact_solution(xi, (double) iy, zeta, &Pface[iy][1][0]);
               }
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(iz = 0; iz < 2; iz++) {
                  exact_solution(xi, eta, (double) iz, &Pface[iz][2][0]);
               }
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  Pxi = xi * Pface[1][0][m] + (1.0 - xi) * Pface[0][0][m];
                  Peta = eta * Pface[1][1][m] + (1.0 - eta) * Pface[0][1][m];
                  Pzeta = zeta * Pface[1][2][m] + (1.0 - zeta) * Pface[0][2][m];
                  u[i][j][k][m] = Pxi + Peta + Pzeta - Pxi * Peta - Pxi * Pzeta - Peta * Pzeta + Pxi * Peta * Pzeta;
               }
            }
         }
      }
      i = 0;
      xi = 0.0;
      #pragma omp parallel for default(shared) private(j, k, m) firstprivate(dnym1, dnzm1, grid_points, ce, temp)
      for(j = 0; j < grid_points[1]; j++) {
         eta = (double) j * dnym1;
         #pragma omp parallel for default(shared) private(k, m) firstprivate(dnzm1, grid_points, ce, temp)
         for(k = 0; k < grid_points[2]; k++) {
            zeta = (double) k * dnzm1;
            exact_solution(xi, eta, zeta, temp);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               u[i][j][k][m] = temp[m];
            }
         }
      }
      i = grid_points[0] - 1;
      xi = 1.0;
      #pragma omp parallel for default(shared) private(j, k, m) firstprivate(dnym1, dnzm1, grid_points, ce, temp)
      for(j = 0; j < grid_points[1]; j++) {
         eta = (double) j * dnym1;
         #pragma omp parallel for default(shared) private(k, m) firstprivate(dnzm1, grid_points, ce, temp)
         for(k = 0; k < grid_points[2]; k++) {
            zeta = (double) k * dnzm1;
            exact_solution(xi, eta, zeta, temp);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               u[i][j][k][m] = temp[m];
            }
         }
      }
      j = 0;
      eta = 0.0;
      #pragma omp parallel for default(shared) private(i, k, m) firstprivate(dnxm1, dnzm1, grid_points, ce, temp)
      for(i = 0; i < grid_points[0]; i++) {
         xi = (double) i * dnxm1;
         #pragma omp parallel for default(shared) private(k, m) firstprivate(dnzm1, grid_points, ce, temp)
         for(k = 0; k < grid_points[2]; k++) {
            zeta = (double) k * dnzm1;
            exact_solution(xi, eta, zeta, temp);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               u[i][j][k][m] = temp[m];
            }
         }
      }
      j = grid_points[1] - 1;
      eta = 1.0;
      #pragma omp parallel for default(shared) private(i, k, m) firstprivate(dnxm1, dnzm1, grid_points, ce, temp)
      for(i = 0; i < grid_points[0]; i++) {
         xi = (double) i * dnxm1;
         #pragma omp parallel for default(shared) private(k, m) firstprivate(dnzm1, grid_points, ce, temp)
         for(k = 0; k < grid_points[2]; k++) {
            zeta = (double) k * dnzm1;
            exact_solution(xi, eta, zeta, temp);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               u[i][j][k][m] = temp[m];
            }
         }
      }
      k = 0;
      zeta = 0.0;
      #pragma omp parallel for default(shared) private(i, j, m) firstprivate(dnxm1, dnym1, grid_points, ce, temp)
      for(i = 0; i < grid_points[0]; i++) {
         xi = (double) i * dnxm1;
         #pragma omp parallel for default(shared) private(j, m) firstprivate(dnym1, grid_points, ce, temp)
         for(j = 0; j < grid_points[1]; j++) {
            eta = (double) j * dnym1;
            exact_solution(xi, eta, zeta, temp);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               u[i][j][k][m] = temp[m];
            }
         }
      }
      k = grid_points[2] - 1;
      zeta = 1.0;
      #pragma omp parallel for default(shared) private(i, j, m) firstprivate(dnxm1, dnym1, grid_points, ce, temp)
      for(i = 0; i < grid_points[0]; i++) {
         xi = (double) i * dnxm1;
         #pragma omp parallel for default(shared) private(j, m) firstprivate(dnym1, grid_points, ce, temp)
         for(j = 0; j < grid_points[1]; j++) {
            eta = (double) j * dnym1;
            exact_solution(xi, eta, zeta, temp);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               u[i][j][k][m] = temp[m];
            }
         }
      }
   }
}

static void lhsinit() {
   {
      int i, j, k, m, n;
      #pragma omp parallel for default(shared) private(i, j, k, m, n) firstprivate(grid_points)
      for(i = 0; i < grid_points[0]; i++) {
         #pragma omp parallel for default(shared) private(j, k, m, n) firstprivate(grid_points)
         for(j = 0; j < grid_points[1]; j++) {
            #pragma omp parallel for default(shared) private(k, m, n) firstprivate(grid_points)
            for(k = 0; k < grid_points[2]; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  /*************** Clava msgError **************
                  		 Loop Iteration number is too low
                  ****************************************/
                  for(n = 0; n < 5; n++) {
                     lhs[i][j][k][0][m][n] = 0.0;
                     lhs[i][j][k][1][m][n] = 0.0;
                     lhs[i][j][k][2][m][n] = 0.0;
                  }
               }
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(grid_points)
      for(i = 0; i < grid_points[0]; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(grid_points)
         for(j = 0; j < grid_points[1]; j++) {
            #pragma omp parallel for default(shared) private(k, m) firstprivate(grid_points)
            for(k = 0; k < grid_points[2]; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  lhs[i][j][k][1][m][m] = 1.0;
               }
            }
         }
      }
   }
}

static void lhsx() {
   int i, j, k;
   #pragma omp parallel for default(shared) private(j, k, i, tmp1, tmp2, tmp3) firstprivate(c2, c1, con43, c3c4, c1345, dt, tx1, tx2, dx1, dx2, dx3, dx4, dx5, grid_points, u)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k, i, tmp1, tmp2, tmp3) firstprivate(j, c2, c1, con43, c3c4, c1345, dt, tx1, tx2, dx1, dx2, dx3, dx4, dx5, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         #pragma omp parallel for default(shared) private(i, tmp1, tmp2, tmp3) firstprivate(j, k, c2, c1, con43, c3c4, c1345, grid_points, u)
         for(i = 0; i < grid_points[0]; i++) {
            tmp1 = 1.0 / u[i][j][k][0];
            tmp2 = tmp1 * tmp1;
            tmp3 = tmp1 * tmp2;
            fjac[i][j][k][0][0] = 0.0;
            fjac[i][j][k][0][1] = 1.0;
            fjac[i][j][k][0][2] = 0.0;
            fjac[i][j][k][0][3] = 0.0;
            fjac[i][j][k][0][4] = 0.0;
            fjac[i][j][k][1][0] = -(u[i][j][k][1] * tmp2 * u[i][j][k][1]) + c2 * 0.50 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][1][1] = (2.0 - c2) * (u[i][j][k][1] / u[i][j][k][0]);
            fjac[i][j][k][1][2] = -c2 * (u[i][j][k][2] * tmp1);
            fjac[i][j][k][1][3] = -c2 * (u[i][j][k][3] * tmp1);
            fjac[i][j][k][1][4] = c2;
            fjac[i][j][k][2][0] = -(u[i][j][k][1] * u[i][j][k][2]) * tmp2;
            fjac[i][j][k][2][1] = u[i][j][k][2] * tmp1;
            fjac[i][j][k][2][2] = u[i][j][k][1] * tmp1;
            fjac[i][j][k][2][3] = 0.0;
            fjac[i][j][k][2][4] = 0.0;
            fjac[i][j][k][3][0] = -(u[i][j][k][1] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][3][1] = u[i][j][k][3] * tmp1;
            fjac[i][j][k][3][2] = 0.0;
            fjac[i][j][k][3][3] = u[i][j][k][1] * tmp1;
            fjac[i][j][k][3][4] = 0.0;
            fjac[i][j][k][4][0] = (c2 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2 - c1 * (u[i][j][k][4] * tmp1)) * (u[i][j][k][1] * tmp1);
            fjac[i][j][k][4][1] = c1 * u[i][j][k][4] * tmp1 - 0.50 * c2 * (3.0 * u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][4][2] = -c2 * (u[i][j][k][2] * u[i][j][k][1]) * tmp2;
            fjac[i][j][k][4][3] = -c2 * (u[i][j][k][3] * u[i][j][k][1]) * tmp2;
            fjac[i][j][k][4][4] = c1 * (u[i][j][k][1] * tmp1);
            njac[i][j][k][0][0] = 0.0;
            njac[i][j][k][0][1] = 0.0;
            njac[i][j][k][0][2] = 0.0;
            njac[i][j][k][0][3] = 0.0;
            njac[i][j][k][0][4] = 0.0;
            njac[i][j][k][1][0] = -con43 * c3c4 * tmp2 * u[i][j][k][1];
            njac[i][j][k][1][1] = con43 * c3c4 * tmp1;
            njac[i][j][k][1][2] = 0.0;
            njac[i][j][k][1][3] = 0.0;
            njac[i][j][k][1][4] = 0.0;
            njac[i][j][k][2][0] = -c3c4 * tmp2 * u[i][j][k][2];
            njac[i][j][k][2][1] = 0.0;
            njac[i][j][k][2][2] = c3c4 * tmp1;
            njac[i][j][k][2][3] = 0.0;
            njac[i][j][k][2][4] = 0.0;
            njac[i][j][k][3][0] = -c3c4 * tmp2 * u[i][j][k][3];
            njac[i][j][k][3][1] = 0.0;
            njac[i][j][k][3][2] = 0.0;
            njac[i][j][k][3][3] = c3c4 * tmp1;
            njac[i][j][k][3][4] = 0.0;
            njac[i][j][k][4][0] = -(con43 * c3c4 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (c3c4 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (c3c4 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - c1345 * tmp2 * u[i][j][k][4];
            njac[i][j][k][4][1] = (con43 * c3c4 - c1345) * tmp2 * u[i][j][k][1];
            njac[i][j][k][4][2] = (c3c4 - c1345) * tmp2 * u[i][j][k][2];
            njac[i][j][k][4][3] = (c3c4 - c1345) * tmp2 * u[i][j][k][3];
            njac[i][j][k][4][4] = (c1345) * tmp1;
         }
         #pragma omp parallel for default(shared) private(i, tmp1, tmp2) firstprivate(dt, tx1, tx2, j, k, dx1, dx2, dx3, dx4, dx5, grid_points, fjac, njac)
         for(i = 1; i < grid_points[0] - 1; i++) {
            tmp1 = dt * tx1;
            tmp2 = dt * tx2;
            lhs[i][j][k][0][0][0] = -tmp2 * fjac[i - 1][j][k][0][0] - tmp1 * njac[i - 1][j][k][0][0] - tmp1 * dx1;
            lhs[i][j][k][0][0][1] = -tmp2 * fjac[i - 1][j][k][0][1] - tmp1 * njac[i - 1][j][k][0][1];
            lhs[i][j][k][0][0][2] = -tmp2 * fjac[i - 1][j][k][0][2] - tmp1 * njac[i - 1][j][k][0][2];
            lhs[i][j][k][0][0][3] = -tmp2 * fjac[i - 1][j][k][0][3] - tmp1 * njac[i - 1][j][k][0][3];
            lhs[i][j][k][0][0][4] = -tmp2 * fjac[i - 1][j][k][0][4] - tmp1 * njac[i - 1][j][k][0][4];
            lhs[i][j][k][0][1][0] = -tmp2 * fjac[i - 1][j][k][1][0] - tmp1 * njac[i - 1][j][k][1][0];
            lhs[i][j][k][0][1][1] = -tmp2 * fjac[i - 1][j][k][1][1] - tmp1 * njac[i - 1][j][k][1][1] - tmp1 * dx2;
            lhs[i][j][k][0][1][2] = -tmp2 * fjac[i - 1][j][k][1][2] - tmp1 * njac[i - 1][j][k][1][2];
            lhs[i][j][k][0][1][3] = -tmp2 * fjac[i - 1][j][k][1][3] - tmp1 * njac[i - 1][j][k][1][3];
            lhs[i][j][k][0][1][4] = -tmp2 * fjac[i - 1][j][k][1][4] - tmp1 * njac[i - 1][j][k][1][4];
            lhs[i][j][k][0][2][0] = -tmp2 * fjac[i - 1][j][k][2][0] - tmp1 * njac[i - 1][j][k][2][0];
            lhs[i][j][k][0][2][1] = -tmp2 * fjac[i - 1][j][k][2][1] - tmp1 * njac[i - 1][j][k][2][1];
            lhs[i][j][k][0][2][2] = -tmp2 * fjac[i - 1][j][k][2][2] - tmp1 * njac[i - 1][j][k][2][2] - tmp1 * dx3;
            lhs[i][j][k][0][2][3] = -tmp2 * fjac[i - 1][j][k][2][3] - tmp1 * njac[i - 1][j][k][2][3];
            lhs[i][j][k][0][2][4] = -tmp2 * fjac[i - 1][j][k][2][4] - tmp1 * njac[i - 1][j][k][2][4];
            lhs[i][j][k][0][3][0] = -tmp2 * fjac[i - 1][j][k][3][0] - tmp1 * njac[i - 1][j][k][3][0];
            lhs[i][j][k][0][3][1] = -tmp2 * fjac[i - 1][j][k][3][1] - tmp1 * njac[i - 1][j][k][3][1];
            lhs[i][j][k][0][3][2] = -tmp2 * fjac[i - 1][j][k][3][2] - tmp1 * njac[i - 1][j][k][3][2];
            lhs[i][j][k][0][3][3] = -tmp2 * fjac[i - 1][j][k][3][3] - tmp1 * njac[i - 1][j][k][3][3] - tmp1 * dx4;
            lhs[i][j][k][0][3][4] = -tmp2 * fjac[i - 1][j][k][3][4] - tmp1 * njac[i - 1][j][k][3][4];
            lhs[i][j][k][0][4][0] = -tmp2 * fjac[i - 1][j][k][4][0] - tmp1 * njac[i - 1][j][k][4][0];
            lhs[i][j][k][0][4][1] = -tmp2 * fjac[i - 1][j][k][4][1] - tmp1 * njac[i - 1][j][k][4][1];
            lhs[i][j][k][0][4][2] = -tmp2 * fjac[i - 1][j][k][4][2] - tmp1 * njac[i - 1][j][k][4][2];
            lhs[i][j][k][0][4][3] = -tmp2 * fjac[i - 1][j][k][4][3] - tmp1 * njac[i - 1][j][k][4][3];
            lhs[i][j][k][0][4][4] = -tmp2 * fjac[i - 1][j][k][4][4] - tmp1 * njac[i - 1][j][k][4][4] - tmp1 * dx5;
            lhs[i][j][k][1][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dx1;
            lhs[i][j][k][1][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
            lhs[i][j][k][1][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
            lhs[i][j][k][1][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
            lhs[i][j][k][1][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
            lhs[i][j][k][1][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
            lhs[i][j][k][1][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dx2;
            lhs[i][j][k][1][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
            lhs[i][j][k][1][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
            lhs[i][j][k][1][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
            lhs[i][j][k][1][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
            lhs[i][j][k][1][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
            lhs[i][j][k][1][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dx3;
            lhs[i][j][k][1][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
            lhs[i][j][k][1][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
            lhs[i][j][k][1][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
            lhs[i][j][k][1][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
            lhs[i][j][k][1][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
            lhs[i][j][k][1][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dx4;
            lhs[i][j][k][1][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
            lhs[i][j][k][1][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
            lhs[i][j][k][1][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
            lhs[i][j][k][1][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
            lhs[i][j][k][1][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
            lhs[i][j][k][1][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dx5;
            lhs[i][j][k][2][0][0] = tmp2 * fjac[i + 1][j][k][0][0] - tmp1 * njac[i + 1][j][k][0][0] - tmp1 * dx1;
            lhs[i][j][k][2][0][1] = tmp2 * fjac[i + 1][j][k][0][1] - tmp1 * njac[i + 1][j][k][0][1];
            lhs[i][j][k][2][0][2] = tmp2 * fjac[i + 1][j][k][0][2] - tmp1 * njac[i + 1][j][k][0][2];
            lhs[i][j][k][2][0][3] = tmp2 * fjac[i + 1][j][k][0][3] - tmp1 * njac[i + 1][j][k][0][3];
            lhs[i][j][k][2][0][4] = tmp2 * fjac[i + 1][j][k][0][4] - tmp1 * njac[i + 1][j][k][0][4];
            lhs[i][j][k][2][1][0] = tmp2 * fjac[i + 1][j][k][1][0] - tmp1 * njac[i + 1][j][k][1][0];
            lhs[i][j][k][2][1][1] = tmp2 * fjac[i + 1][j][k][1][1] - tmp1 * njac[i + 1][j][k][1][1] - tmp1 * dx2;
            lhs[i][j][k][2][1][2] = tmp2 * fjac[i + 1][j][k][1][2] - tmp1 * njac[i + 1][j][k][1][2];
            lhs[i][j][k][2][1][3] = tmp2 * fjac[i + 1][j][k][1][3] - tmp1 * njac[i + 1][j][k][1][3];
            lhs[i][j][k][2][1][4] = tmp2 * fjac[i + 1][j][k][1][4] - tmp1 * njac[i + 1][j][k][1][4];
            lhs[i][j][k][2][2][0] = tmp2 * fjac[i + 1][j][k][2][0] - tmp1 * njac[i + 1][j][k][2][0];
            lhs[i][j][k][2][2][1] = tmp2 * fjac[i + 1][j][k][2][1] - tmp1 * njac[i + 1][j][k][2][1];
            lhs[i][j][k][2][2][2] = tmp2 * fjac[i + 1][j][k][2][2] - tmp1 * njac[i + 1][j][k][2][2] - tmp1 * dx3;
            lhs[i][j][k][2][2][3] = tmp2 * fjac[i + 1][j][k][2][3] - tmp1 * njac[i + 1][j][k][2][3];
            lhs[i][j][k][2][2][4] = tmp2 * fjac[i + 1][j][k][2][4] - tmp1 * njac[i + 1][j][k][2][4];
            lhs[i][j][k][2][3][0] = tmp2 * fjac[i + 1][j][k][3][0] - tmp1 * njac[i + 1][j][k][3][0];
            lhs[i][j][k][2][3][1] = tmp2 * fjac[i + 1][j][k][3][1] - tmp1 * njac[i + 1][j][k][3][1];
            lhs[i][j][k][2][3][2] = tmp2 * fjac[i + 1][j][k][3][2] - tmp1 * njac[i + 1][j][k][3][2];
            lhs[i][j][k][2][3][3] = tmp2 * fjac[i + 1][j][k][3][3] - tmp1 * njac[i + 1][j][k][3][3] - tmp1 * dx4;
            lhs[i][j][k][2][3][4] = tmp2 * fjac[i + 1][j][k][3][4] - tmp1 * njac[i + 1][j][k][3][4];
            lhs[i][j][k][2][4][0] = tmp2 * fjac[i + 1][j][k][4][0] - tmp1 * njac[i + 1][j][k][4][0];
            lhs[i][j][k][2][4][1] = tmp2 * fjac[i + 1][j][k][4][1] - tmp1 * njac[i + 1][j][k][4][1];
            lhs[i][j][k][2][4][2] = tmp2 * fjac[i + 1][j][k][4][2] - tmp1 * njac[i + 1][j][k][4][2];
            lhs[i][j][k][2][4][3] = tmp2 * fjac[i + 1][j][k][4][3] - tmp1 * njac[i + 1][j][k][4][3];
            lhs[i][j][k][2][4][4] = tmp2 * fjac[i + 1][j][k][4][4] - tmp1 * njac[i + 1][j][k][4][4] - tmp1 * dx5;
         }
      }
   }
}

static void lhsy() {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k, tmp1, tmp2, tmp3) firstprivate(c2, c1, c3c4, con43, c1345, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, tmp1, tmp2, tmp3) firstprivate(i, c2, c1, c3c4, con43, c1345, grid_points, u)
      for(j = 0; j < grid_points[1]; j++) {
         #pragma omp parallel for default(shared) private(k, tmp1, tmp2, tmp3) firstprivate(i, j, c2, c1, c3c4, con43, c1345, grid_points, u)
         for(k = 1; k < grid_points[2] - 1; k++) {
            tmp1 = 1.0 / u[i][j][k][0];
            tmp2 = tmp1 * tmp1;
            tmp3 = tmp1 * tmp2;
            fjac[i][j][k][0][0] = 0.0;
            fjac[i][j][k][0][1] = 0.0;
            fjac[i][j][k][0][2] = 1.0;
            fjac[i][j][k][0][3] = 0.0;
            fjac[i][j][k][0][4] = 0.0;
            fjac[i][j][k][1][0] = -(u[i][j][k][1] * u[i][j][k][2]) * tmp2;
            fjac[i][j][k][1][1] = u[i][j][k][2] * tmp1;
            fjac[i][j][k][1][2] = u[i][j][k][1] * tmp1;
            fjac[i][j][k][1][3] = 0.0;
            fjac[i][j][k][1][4] = 0.0;
            fjac[i][j][k][2][0] = -(u[i][j][k][2] * u[i][j][k][2] * tmp2) + 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2);
            fjac[i][j][k][2][1] = -c2 * u[i][j][k][1] * tmp1;
            fjac[i][j][k][2][2] = (2.0 - c2) * u[i][j][k][2] * tmp1;
            fjac[i][j][k][2][3] = -c2 * u[i][j][k][3] * tmp1;
            fjac[i][j][k][2][4] = c2;
            fjac[i][j][k][3][0] = -(u[i][j][k][2] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][3][1] = 0.0;
            fjac[i][j][k][3][2] = u[i][j][k][3] * tmp1;
            fjac[i][j][k][3][3] = u[i][j][k][2] * tmp1;
            fjac[i][j][k][3][4] = 0.0;
            fjac[i][j][k][4][0] = (c2 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2 - c1 * u[i][j][k][4] * tmp1) * u[i][j][k][2] * tmp1;
            fjac[i][j][k][4][1] = -c2 * u[i][j][k][1] * u[i][j][k][2] * tmp2;
            fjac[i][j][k][4][2] = c1 * u[i][j][k][4] * tmp1 - 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + 3.0 * u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2);
            fjac[i][j][k][4][3] = -c2 * (u[i][j][k][2] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][4][4] = c1 * u[i][j][k][2] * tmp1;
            njac[i][j][k][0][0] = 0.0;
            njac[i][j][k][0][1] = 0.0;
            njac[i][j][k][0][2] = 0.0;
            njac[i][j][k][0][3] = 0.0;
            njac[i][j][k][0][4] = 0.0;
            njac[i][j][k][1][0] = -c3c4 * tmp2 * u[i][j][k][1];
            njac[i][j][k][1][1] = c3c4 * tmp1;
            njac[i][j][k][1][2] = 0.0;
            njac[i][j][k][1][3] = 0.0;
            njac[i][j][k][1][4] = 0.0;
            njac[i][j][k][2][0] = -con43 * c3c4 * tmp2 * u[i][j][k][2];
            njac[i][j][k][2][1] = 0.0;
            njac[i][j][k][2][2] = con43 * c3c4 * tmp1;
            njac[i][j][k][2][3] = 0.0;
            njac[i][j][k][2][4] = 0.0;
            njac[i][j][k][3][0] = -c3c4 * tmp2 * u[i][j][k][3];
            njac[i][j][k][3][1] = 0.0;
            njac[i][j][k][3][2] = 0.0;
            njac[i][j][k][3][3] = c3c4 * tmp1;
            njac[i][j][k][3][4] = 0.0;
            njac[i][j][k][4][0] = -(c3c4 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (con43 * c3c4 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (c3c4 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - c1345 * tmp2 * u[i][j][k][4];
            njac[i][j][k][4][1] = (c3c4 - c1345) * tmp2 * u[i][j][k][1];
            njac[i][j][k][4][2] = (con43 * c3c4 - c1345) * tmp2 * u[i][j][k][2];
            njac[i][j][k][4][3] = (c3c4 - c1345) * tmp2 * u[i][j][k][3];
            njac[i][j][k][4][4] = (c1345) * tmp1;
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, tmp1, tmp2) firstprivate(dt, ty1, ty2, dy1, dy2, dy3, dy4, dy5, grid_points, fjac, njac)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, tmp1, tmp2) firstprivate(dt, ty1, ty2, i, dy1, dy2, dy3, dy4, dy5, grid_points, fjac, njac)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, tmp1, tmp2) firstprivate(dt, ty1, ty2, j, i, dy1, dy2, dy3, dy4, dy5, grid_points, fjac, njac)
         for(k = 1; k < grid_points[2] - 1; k++) {
            tmp1 = dt * ty1;
            tmp2 = dt * ty2;
            lhs[i][j][k][0][0][0] = -tmp2 * fjac[i][j - 1][k][0][0] - tmp1 * njac[i][j - 1][k][0][0] - tmp1 * dy1;
            lhs[i][j][k][0][0][1] = -tmp2 * fjac[i][j - 1][k][0][1] - tmp1 * njac[i][j - 1][k][0][1];
            lhs[i][j][k][0][0][2] = -tmp2 * fjac[i][j - 1][k][0][2] - tmp1 * njac[i][j - 1][k][0][2];
            lhs[i][j][k][0][0][3] = -tmp2 * fjac[i][j - 1][k][0][3] - tmp1 * njac[i][j - 1][k][0][3];
            lhs[i][j][k][0][0][4] = -tmp2 * fjac[i][j - 1][k][0][4] - tmp1 * njac[i][j - 1][k][0][4];
            lhs[i][j][k][0][1][0] = -tmp2 * fjac[i][j - 1][k][1][0] - tmp1 * njac[i][j - 1][k][1][0];
            lhs[i][j][k][0][1][1] = -tmp2 * fjac[i][j - 1][k][1][1] - tmp1 * njac[i][j - 1][k][1][1] - tmp1 * dy2;
            lhs[i][j][k][0][1][2] = -tmp2 * fjac[i][j - 1][k][1][2] - tmp1 * njac[i][j - 1][k][1][2];
            lhs[i][j][k][0][1][3] = -tmp2 * fjac[i][j - 1][k][1][3] - tmp1 * njac[i][j - 1][k][1][3];
            lhs[i][j][k][0][1][4] = -tmp2 * fjac[i][j - 1][k][1][4] - tmp1 * njac[i][j - 1][k][1][4];
            lhs[i][j][k][0][2][0] = -tmp2 * fjac[i][j - 1][k][2][0] - tmp1 * njac[i][j - 1][k][2][0];
            lhs[i][j][k][0][2][1] = -tmp2 * fjac[i][j - 1][k][2][1] - tmp1 * njac[i][j - 1][k][2][1];
            lhs[i][j][k][0][2][2] = -tmp2 * fjac[i][j - 1][k][2][2] - tmp1 * njac[i][j - 1][k][2][2] - tmp1 * dy3;
            lhs[i][j][k][0][2][3] = -tmp2 * fjac[i][j - 1][k][2][3] - tmp1 * njac[i][j - 1][k][2][3];
            lhs[i][j][k][0][2][4] = -tmp2 * fjac[i][j - 1][k][2][4] - tmp1 * njac[i][j - 1][k][2][4];
            lhs[i][j][k][0][3][0] = -tmp2 * fjac[i][j - 1][k][3][0] - tmp1 * njac[i][j - 1][k][3][0];
            lhs[i][j][k][0][3][1] = -tmp2 * fjac[i][j - 1][k][3][1] - tmp1 * njac[i][j - 1][k][3][1];
            lhs[i][j][k][0][3][2] = -tmp2 * fjac[i][j - 1][k][3][2] - tmp1 * njac[i][j - 1][k][3][2];
            lhs[i][j][k][0][3][3] = -tmp2 * fjac[i][j - 1][k][3][3] - tmp1 * njac[i][j - 1][k][3][3] - tmp1 * dy4;
            lhs[i][j][k][0][3][4] = -tmp2 * fjac[i][j - 1][k][3][4] - tmp1 * njac[i][j - 1][k][3][4];
            lhs[i][j][k][0][4][0] = -tmp2 * fjac[i][j - 1][k][4][0] - tmp1 * njac[i][j - 1][k][4][0];
            lhs[i][j][k][0][4][1] = -tmp2 * fjac[i][j - 1][k][4][1] - tmp1 * njac[i][j - 1][k][4][1];
            lhs[i][j][k][0][4][2] = -tmp2 * fjac[i][j - 1][k][4][2] - tmp1 * njac[i][j - 1][k][4][2];
            lhs[i][j][k][0][4][3] = -tmp2 * fjac[i][j - 1][k][4][3] - tmp1 * njac[i][j - 1][k][4][3];
            lhs[i][j][k][0][4][4] = -tmp2 * fjac[i][j - 1][k][4][4] - tmp1 * njac[i][j - 1][k][4][4] - tmp1 * dy5;
            lhs[i][j][k][1][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dy1;
            lhs[i][j][k][1][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
            lhs[i][j][k][1][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
            lhs[i][j][k][1][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
            lhs[i][j][k][1][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
            lhs[i][j][k][1][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
            lhs[i][j][k][1][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dy2;
            lhs[i][j][k][1][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
            lhs[i][j][k][1][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
            lhs[i][j][k][1][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
            lhs[i][j][k][1][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
            lhs[i][j][k][1][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
            lhs[i][j][k][1][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dy3;
            lhs[i][j][k][1][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
            lhs[i][j][k][1][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
            lhs[i][j][k][1][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
            lhs[i][j][k][1][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
            lhs[i][j][k][1][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
            lhs[i][j][k][1][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dy4;
            lhs[i][j][k][1][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
            lhs[i][j][k][1][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
            lhs[i][j][k][1][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
            lhs[i][j][k][1][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
            lhs[i][j][k][1][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
            lhs[i][j][k][1][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dy5;
            lhs[i][j][k][2][0][0] = tmp2 * fjac[i][j + 1][k][0][0] - tmp1 * njac[i][j + 1][k][0][0] - tmp1 * dy1;
            lhs[i][j][k][2][0][1] = tmp2 * fjac[i][j + 1][k][0][1] - tmp1 * njac[i][j + 1][k][0][1];
            lhs[i][j][k][2][0][2] = tmp2 * fjac[i][j + 1][k][0][2] - tmp1 * njac[i][j + 1][k][0][2];
            lhs[i][j][k][2][0][3] = tmp2 * fjac[i][j + 1][k][0][3] - tmp1 * njac[i][j + 1][k][0][3];
            lhs[i][j][k][2][0][4] = tmp2 * fjac[i][j + 1][k][0][4] - tmp1 * njac[i][j + 1][k][0][4];
            lhs[i][j][k][2][1][0] = tmp2 * fjac[i][j + 1][k][1][0] - tmp1 * njac[i][j + 1][k][1][0];
            lhs[i][j][k][2][1][1] = tmp2 * fjac[i][j + 1][k][1][1] - tmp1 * njac[i][j + 1][k][1][1] - tmp1 * dy2;
            lhs[i][j][k][2][1][2] = tmp2 * fjac[i][j + 1][k][1][2] - tmp1 * njac[i][j + 1][k][1][2];
            lhs[i][j][k][2][1][3] = tmp2 * fjac[i][j + 1][k][1][3] - tmp1 * njac[i][j + 1][k][1][3];
            lhs[i][j][k][2][1][4] = tmp2 * fjac[i][j + 1][k][1][4] - tmp1 * njac[i][j + 1][k][1][4];
            lhs[i][j][k][2][2][0] = tmp2 * fjac[i][j + 1][k][2][0] - tmp1 * njac[i][j + 1][k][2][0];
            lhs[i][j][k][2][2][1] = tmp2 * fjac[i][j + 1][k][2][1] - tmp1 * njac[i][j + 1][k][2][1];
            lhs[i][j][k][2][2][2] = tmp2 * fjac[i][j + 1][k][2][2] - tmp1 * njac[i][j + 1][k][2][2] - tmp1 * dy3;
            lhs[i][j][k][2][2][3] = tmp2 * fjac[i][j + 1][k][2][3] - tmp1 * njac[i][j + 1][k][2][3];
            lhs[i][j][k][2][2][4] = tmp2 * fjac[i][j + 1][k][2][4] - tmp1 * njac[i][j + 1][k][2][4];
            lhs[i][j][k][2][3][0] = tmp2 * fjac[i][j + 1][k][3][0] - tmp1 * njac[i][j + 1][k][3][0];
            lhs[i][j][k][2][3][1] = tmp2 * fjac[i][j + 1][k][3][1] - tmp1 * njac[i][j + 1][k][3][1];
            lhs[i][j][k][2][3][2] = tmp2 * fjac[i][j + 1][k][3][2] - tmp1 * njac[i][j + 1][k][3][2];
            lhs[i][j][k][2][3][3] = tmp2 * fjac[i][j + 1][k][3][3] - tmp1 * njac[i][j + 1][k][3][3] - tmp1 * dy4;
            lhs[i][j][k][2][3][4] = tmp2 * fjac[i][j + 1][k][3][4] - tmp1 * njac[i][j + 1][k][3][4];
            lhs[i][j][k][2][4][0] = tmp2 * fjac[i][j + 1][k][4][0] - tmp1 * njac[i][j + 1][k][4][0];
            lhs[i][j][k][2][4][1] = tmp2 * fjac[i][j + 1][k][4][1] - tmp1 * njac[i][j + 1][k][4][1];
            lhs[i][j][k][2][4][2] = tmp2 * fjac[i][j + 1][k][4][2] - tmp1 * njac[i][j + 1][k][4][2];
            lhs[i][j][k][2][4][3] = tmp2 * fjac[i][j + 1][k][4][3] - tmp1 * njac[i][j + 1][k][4][3];
            lhs[i][j][k][2][4][4] = tmp2 * fjac[i][j + 1][k][4][4] - tmp1 * njac[i][j + 1][k][4][4] - tmp1 * dy5;
         }
      }
   }
}

static void lhsz() {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k, tmp1, tmp2, tmp3) firstprivate(c2, c1, c3c4, con43, c3, c4, c1345, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, tmp1, tmp2, tmp3) firstprivate(i, c2, c1, c3c4, con43, c3, c4, c1345, grid_points, u)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, tmp1, tmp2, tmp3) firstprivate(i, j, c2, c1, c3c4, con43, c3, c4, c1345, grid_points, u)
         for(k = 0; k < grid_points[2]; k++) {
            tmp1 = 1.0 / u[i][j][k][0];
            tmp2 = tmp1 * tmp1;
            tmp3 = tmp1 * tmp2;
            fjac[i][j][k][0][0] = 0.0;
            fjac[i][j][k][0][1] = 0.0;
            fjac[i][j][k][0][2] = 0.0;
            fjac[i][j][k][0][3] = 1.0;
            fjac[i][j][k][0][4] = 0.0;
            fjac[i][j][k][1][0] = -(u[i][j][k][1] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][1][1] = u[i][j][k][3] * tmp1;
            fjac[i][j][k][1][2] = 0.0;
            fjac[i][j][k][1][3] = u[i][j][k][1] * tmp1;
            fjac[i][j][k][1][4] = 0.0;
            fjac[i][j][k][2][0] = -(u[i][j][k][2] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][2][1] = 0.0;
            fjac[i][j][k][2][2] = u[i][j][k][3] * tmp1;
            fjac[i][j][k][2][3] = u[i][j][k][2] * tmp1;
            fjac[i][j][k][2][4] = 0.0;
            fjac[i][j][k][3][0] = -(u[i][j][k][3] * u[i][j][k][3] * tmp2) + 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2);
            fjac[i][j][k][3][1] = -c2 * u[i][j][k][1] * tmp1;
            fjac[i][j][k][3][2] = -c2 * u[i][j][k][2] * tmp1;
            fjac[i][j][k][3][3] = (2.0 - c2) * u[i][j][k][3] * tmp1;
            fjac[i][j][k][3][4] = c2;
            fjac[i][j][k][4][0] = (c2 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * tmp2 - c1 * (u[i][j][k][4] * tmp1)) * (u[i][j][k][3] * tmp1);
            fjac[i][j][k][4][1] = -c2 * (u[i][j][k][1] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][4][2] = -c2 * (u[i][j][k][2] * u[i][j][k][3]) * tmp2;
            fjac[i][j][k][4][3] = c1 * (u[i][j][k][4] * tmp1) - 0.50 * c2 * ((u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + 3.0 * u[i][j][k][3] * u[i][j][k][3]) * tmp2);
            fjac[i][j][k][4][4] = c1 * u[i][j][k][3] * tmp1;
            njac[i][j][k][0][0] = 0.0;
            njac[i][j][k][0][1] = 0.0;
            njac[i][j][k][0][2] = 0.0;
            njac[i][j][k][0][3] = 0.0;
            njac[i][j][k][0][4] = 0.0;
            njac[i][j][k][1][0] = -c3c4 * tmp2 * u[i][j][k][1];
            njac[i][j][k][1][1] = c3c4 * tmp1;
            njac[i][j][k][1][2] = 0.0;
            njac[i][j][k][1][3] = 0.0;
            njac[i][j][k][1][4] = 0.0;
            njac[i][j][k][2][0] = -c3c4 * tmp2 * u[i][j][k][2];
            njac[i][j][k][2][1] = 0.0;
            njac[i][j][k][2][2] = c3c4 * tmp1;
            njac[i][j][k][2][3] = 0.0;
            njac[i][j][k][2][4] = 0.0;
            njac[i][j][k][3][0] = -con43 * c3c4 * tmp2 * u[i][j][k][3];
            njac[i][j][k][3][1] = 0.0;
            njac[i][j][k][3][2] = 0.0;
            njac[i][j][k][3][3] = con43 * c3 * c4 * tmp1;
            njac[i][j][k][3][4] = 0.0;
            njac[i][j][k][4][0] = -(c3c4 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (c3c4 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (con43 * c3c4 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - c1345 * tmp2 * u[i][j][k][4];
            njac[i][j][k][4][1] = (c3c4 - c1345) * tmp2 * u[i][j][k][1];
            njac[i][j][k][4][2] = (c3c4 - c1345) * tmp2 * u[i][j][k][2];
            njac[i][j][k][4][3] = (con43 * c3c4 - c1345) * tmp2 * u[i][j][k][3];
            njac[i][j][k][4][4] = (c1345) * tmp1;
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, tmp1, tmp2) firstprivate(dt, tz1, tz2, dz1, dz2, dz3, dz4, dz5, grid_points, fjac, njac)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, tmp1, tmp2) firstprivate(dt, tz1, tz2, i, dz1, dz2, dz3, dz4, dz5, grid_points, fjac, njac)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, tmp1, tmp2) firstprivate(dt, tz1, tz2, i, j, dz1, dz2, dz3, dz4, dz5, grid_points, fjac, njac)
         for(k = 1; k < grid_points[2] - 1; k++) {
            tmp1 = dt * tz1;
            tmp2 = dt * tz2;
            lhs[i][j][k][0][0][0] = -tmp2 * fjac[i][j][k - 1][0][0] - tmp1 * njac[i][j][k - 1][0][0] - tmp1 * dz1;
            lhs[i][j][k][0][0][1] = -tmp2 * fjac[i][j][k - 1][0][1] - tmp1 * njac[i][j][k - 1][0][1];
            lhs[i][j][k][0][0][2] = -tmp2 * fjac[i][j][k - 1][0][2] - tmp1 * njac[i][j][k - 1][0][2];
            lhs[i][j][k][0][0][3] = -tmp2 * fjac[i][j][k - 1][0][3] - tmp1 * njac[i][j][k - 1][0][3];
            lhs[i][j][k][0][0][4] = -tmp2 * fjac[i][j][k - 1][0][4] - tmp1 * njac[i][j][k - 1][0][4];
            lhs[i][j][k][0][1][0] = -tmp2 * fjac[i][j][k - 1][1][0] - tmp1 * njac[i][j][k - 1][1][0];
            lhs[i][j][k][0][1][1] = -tmp2 * fjac[i][j][k - 1][1][1] - tmp1 * njac[i][j][k - 1][1][1] - tmp1 * dz2;
            lhs[i][j][k][0][1][2] = -tmp2 * fjac[i][j][k - 1][1][2] - tmp1 * njac[i][j][k - 1][1][2];
            lhs[i][j][k][0][1][3] = -tmp2 * fjac[i][j][k - 1][1][3] - tmp1 * njac[i][j][k - 1][1][3];
            lhs[i][j][k][0][1][4] = -tmp2 * fjac[i][j][k - 1][1][4] - tmp1 * njac[i][j][k - 1][1][4];
            lhs[i][j][k][0][2][0] = -tmp2 * fjac[i][j][k - 1][2][0] - tmp1 * njac[i][j][k - 1][2][0];
            lhs[i][j][k][0][2][1] = -tmp2 * fjac[i][j][k - 1][2][1] - tmp1 * njac[i][j][k - 1][2][1];
            lhs[i][j][k][0][2][2] = -tmp2 * fjac[i][j][k - 1][2][2] - tmp1 * njac[i][j][k - 1][2][2] - tmp1 * dz3;
            lhs[i][j][k][0][2][3] = -tmp2 * fjac[i][j][k - 1][2][3] - tmp1 * njac[i][j][k - 1][2][3];
            lhs[i][j][k][0][2][4] = -tmp2 * fjac[i][j][k - 1][2][4] - tmp1 * njac[i][j][k - 1][2][4];
            lhs[i][j][k][0][3][0] = -tmp2 * fjac[i][j][k - 1][3][0] - tmp1 * njac[i][j][k - 1][3][0];
            lhs[i][j][k][0][3][1] = -tmp2 * fjac[i][j][k - 1][3][1] - tmp1 * njac[i][j][k - 1][3][1];
            lhs[i][j][k][0][3][2] = -tmp2 * fjac[i][j][k - 1][3][2] - tmp1 * njac[i][j][k - 1][3][2];
            lhs[i][j][k][0][3][3] = -tmp2 * fjac[i][j][k - 1][3][3] - tmp1 * njac[i][j][k - 1][3][3] - tmp1 * dz4;
            lhs[i][j][k][0][3][4] = -tmp2 * fjac[i][j][k - 1][3][4] - tmp1 * njac[i][j][k - 1][3][4];
            lhs[i][j][k][0][4][0] = -tmp2 * fjac[i][j][k - 1][4][0] - tmp1 * njac[i][j][k - 1][4][0];
            lhs[i][j][k][0][4][1] = -tmp2 * fjac[i][j][k - 1][4][1] - tmp1 * njac[i][j][k - 1][4][1];
            lhs[i][j][k][0][4][2] = -tmp2 * fjac[i][j][k - 1][4][2] - tmp1 * njac[i][j][k - 1][4][2];
            lhs[i][j][k][0][4][3] = -tmp2 * fjac[i][j][k - 1][4][3] - tmp1 * njac[i][j][k - 1][4][3];
            lhs[i][j][k][0][4][4] = -tmp2 * fjac[i][j][k - 1][4][4] - tmp1 * njac[i][j][k - 1][4][4] - tmp1 * dz5;
            lhs[i][j][k][1][0][0] = 1.0 + tmp1 * 2.0 * njac[i][j][k][0][0] + tmp1 * 2.0 * dz1;
            lhs[i][j][k][1][0][1] = tmp1 * 2.0 * njac[i][j][k][0][1];
            lhs[i][j][k][1][0][2] = tmp1 * 2.0 * njac[i][j][k][0][2];
            lhs[i][j][k][1][0][3] = tmp1 * 2.0 * njac[i][j][k][0][3];
            lhs[i][j][k][1][0][4] = tmp1 * 2.0 * njac[i][j][k][0][4];
            lhs[i][j][k][1][1][0] = tmp1 * 2.0 * njac[i][j][k][1][0];
            lhs[i][j][k][1][1][1] = 1.0 + tmp1 * 2.0 * njac[i][j][k][1][1] + tmp1 * 2.0 * dz2;
            lhs[i][j][k][1][1][2] = tmp1 * 2.0 * njac[i][j][k][1][2];
            lhs[i][j][k][1][1][3] = tmp1 * 2.0 * njac[i][j][k][1][3];
            lhs[i][j][k][1][1][4] = tmp1 * 2.0 * njac[i][j][k][1][4];
            lhs[i][j][k][1][2][0] = tmp1 * 2.0 * njac[i][j][k][2][0];
            lhs[i][j][k][1][2][1] = tmp1 * 2.0 * njac[i][j][k][2][1];
            lhs[i][j][k][1][2][2] = 1.0 + tmp1 * 2.0 * njac[i][j][k][2][2] + tmp1 * 2.0 * dz3;
            lhs[i][j][k][1][2][3] = tmp1 * 2.0 * njac[i][j][k][2][3];
            lhs[i][j][k][1][2][4] = tmp1 * 2.0 * njac[i][j][k][2][4];
            lhs[i][j][k][1][3][0] = tmp1 * 2.0 * njac[i][j][k][3][0];
            lhs[i][j][k][1][3][1] = tmp1 * 2.0 * njac[i][j][k][3][1];
            lhs[i][j][k][1][3][2] = tmp1 * 2.0 * njac[i][j][k][3][2];
            lhs[i][j][k][1][3][3] = 1.0 + tmp1 * 2.0 * njac[i][j][k][3][3] + tmp1 * 2.0 * dz4;
            lhs[i][j][k][1][3][4] = tmp1 * 2.0 * njac[i][j][k][3][4];
            lhs[i][j][k][1][4][0] = tmp1 * 2.0 * njac[i][j][k][4][0];
            lhs[i][j][k][1][4][1] = tmp1 * 2.0 * njac[i][j][k][4][1];
            lhs[i][j][k][1][4][2] = tmp1 * 2.0 * njac[i][j][k][4][2];
            lhs[i][j][k][1][4][3] = tmp1 * 2.0 * njac[i][j][k][4][3];
            lhs[i][j][k][1][4][4] = 1.0 + tmp1 * 2.0 * njac[i][j][k][4][4] + tmp1 * 2.0 * dz5;
            lhs[i][j][k][2][0][0] = tmp2 * fjac[i][j][k + 1][0][0] - tmp1 * njac[i][j][k + 1][0][0] - tmp1 * dz1;
            lhs[i][j][k][2][0][1] = tmp2 * fjac[i][j][k + 1][0][1] - tmp1 * njac[i][j][k + 1][0][1];
            lhs[i][j][k][2][0][2] = tmp2 * fjac[i][j][k + 1][0][2] - tmp1 * njac[i][j][k + 1][0][2];
            lhs[i][j][k][2][0][3] = tmp2 * fjac[i][j][k + 1][0][3] - tmp1 * njac[i][j][k + 1][0][3];
            lhs[i][j][k][2][0][4] = tmp2 * fjac[i][j][k + 1][0][4] - tmp1 * njac[i][j][k + 1][0][4];
            lhs[i][j][k][2][1][0] = tmp2 * fjac[i][j][k + 1][1][0] - tmp1 * njac[i][j][k + 1][1][0];
            lhs[i][j][k][2][1][1] = tmp2 * fjac[i][j][k + 1][1][1] - tmp1 * njac[i][j][k + 1][1][1] - tmp1 * dz2;
            lhs[i][j][k][2][1][2] = tmp2 * fjac[i][j][k + 1][1][2] - tmp1 * njac[i][j][k + 1][1][2];
            lhs[i][j][k][2][1][3] = tmp2 * fjac[i][j][k + 1][1][3] - tmp1 * njac[i][j][k + 1][1][3];
            lhs[i][j][k][2][1][4] = tmp2 * fjac[i][j][k + 1][1][4] - tmp1 * njac[i][j][k + 1][1][4];
            lhs[i][j][k][2][2][0] = tmp2 * fjac[i][j][k + 1][2][0] - tmp1 * njac[i][j][k + 1][2][0];
            lhs[i][j][k][2][2][1] = tmp2 * fjac[i][j][k + 1][2][1] - tmp1 * njac[i][j][k + 1][2][1];
            lhs[i][j][k][2][2][2] = tmp2 * fjac[i][j][k + 1][2][2] - tmp1 * njac[i][j][k + 1][2][2] - tmp1 * dz3;
            lhs[i][j][k][2][2][3] = tmp2 * fjac[i][j][k + 1][2][3] - tmp1 * njac[i][j][k + 1][2][3];
            lhs[i][j][k][2][2][4] = tmp2 * fjac[i][j][k + 1][2][4] - tmp1 * njac[i][j][k + 1][2][4];
            lhs[i][j][k][2][3][0] = tmp2 * fjac[i][j][k + 1][3][0] - tmp1 * njac[i][j][k + 1][3][0];
            lhs[i][j][k][2][3][1] = tmp2 * fjac[i][j][k + 1][3][1] - tmp1 * njac[i][j][k + 1][3][1];
            lhs[i][j][k][2][3][2] = tmp2 * fjac[i][j][k + 1][3][2] - tmp1 * njac[i][j][k + 1][3][2];
            lhs[i][j][k][2][3][3] = tmp2 * fjac[i][j][k + 1][3][3] - tmp1 * njac[i][j][k + 1][3][3] - tmp1 * dz4;
            lhs[i][j][k][2][3][4] = tmp2 * fjac[i][j][k + 1][3][4] - tmp1 * njac[i][j][k + 1][3][4];
            lhs[i][j][k][2][4][0] = tmp2 * fjac[i][j][k + 1][4][0] - tmp1 * njac[i][j][k + 1][4][0];
            lhs[i][j][k][2][4][1] = tmp2 * fjac[i][j][k + 1][4][1] - tmp1 * njac[i][j][k + 1][4][1];
            lhs[i][j][k][2][4][2] = tmp2 * fjac[i][j][k + 1][4][2] - tmp1 * njac[i][j][k + 1][4][2];
            lhs[i][j][k][2][4][3] = tmp2 * fjac[i][j][k + 1][4][3] - tmp1 * njac[i][j][k + 1][4][3];
            lhs[i][j][k][2][4][4] = tmp2 * fjac[i][j][k + 1][4][4] - tmp1 * njac[i][j][k + 1][4][4] - tmp1 * dz5;
         }
      }
   }
}

static void compute_rhs() {
   int i, j, k, m;
   double rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;
   #pragma omp parallel for default(shared) private(i, j, k, rho_inv) firstprivate(grid_points, u)
   for(i = 0; i < grid_points[0]; i++) {
      #pragma omp parallel for default(shared) private(j, k, rho_inv) firstprivate(i, grid_points, u)
      for(j = 0; j < grid_points[1]; j++) {
         #pragma omp parallel for default(shared) private(k, rho_inv) firstprivate(i, j, grid_points, u)
         for(k = 0; k < grid_points[2]; k++) {
            rho_inv = 1.0 / u[i][j][k][0];
            rho_i[i][j][k] = rho_inv;
            us[i][j][k] = u[i][j][k][1] * rho_inv;
            vs[i][j][k] = u[i][j][k][2] * rho_inv;
            ws[i][j][k] = u[i][j][k][3] * rho_inv;
            square[i][j][k] = 0.5 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) * rho_inv;
            qs[i][j][k] = square[i][j][k] * rho_inv;
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(grid_points, forcing)
   for(i = 0; i < grid_points[0]; i++) {
      #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, grid_points, forcing)
      for(j = 0; j < grid_points[1]; j++) {
         #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, grid_points, forcing)
         for(k = 0; k < grid_points[2]; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rhs[i][j][k][m] = forcing[i][j][k][m];
            }
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, uijk, up1, um1) firstprivate(dx1tx1, tx2, c2, dx2tx1, xxcon2, con43, dx3tx1, dx4tx1, c1, dx5tx1, xxcon3, xxcon4, xxcon5, grid_points, us, u, square, vs, ws, qs, rho_i)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, uijk, up1, um1) firstprivate(i, dx1tx1, tx2, c2, dx2tx1, xxcon2, con43, dx3tx1, dx4tx1, c1, dx5tx1, xxcon3, xxcon4, xxcon5, grid_points, us, u, square, vs, ws, qs, rho_i)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, uijk, up1, um1) firstprivate(i, j, dx1tx1, tx2, c2, dx2tx1, xxcon2, con43, dx3tx1, dx4tx1, c1, dx5tx1, xxcon3, xxcon4, xxcon5, grid_points, us, u, square, vs, ws, qs, rho_i)
         for(k = 1; k < grid_points[2] - 1; k++) {
            uijk = us[i][j][k];
            up1 = us[i + 1][j][k];
            um1 = us[i - 1][j][k];
            rhs[i][j][k][0] = rhs[i][j][k][0] + dx1tx1 * (u[i + 1][j][k][0] - 2.0 * u[i][j][k][0] + u[i - 1][j][k][0]) - tx2 * (u[i + 1][j][k][1] - u[i - 1][j][k][1]);
            rhs[i][j][k][1] = rhs[i][j][k][1] + dx2tx1 * (u[i + 1][j][k][1] - 2.0 * u[i][j][k][1] + u[i - 1][j][k][1]) + xxcon2 * con43 * (up1 - 2.0 * uijk + um1) - tx2 * (u[i + 1][j][k][1] * up1 - u[i - 1][j][k][1] * um1 + (u[i + 1][j][k][4] - square[i + 1][j][k] - u[i - 1][j][k][4] + square[i - 1][j][k]) * c2);
            rhs[i][j][k][2] = rhs[i][j][k][2] + dx3tx1 * (u[i + 1][j][k][2] - 2.0 * u[i][j][k][2] + u[i - 1][j][k][2]) + xxcon2 * (vs[i + 1][j][k] - 2.0 * vs[i][j][k] + vs[i - 1][j][k]) - tx2 * (u[i + 1][j][k][2] * up1 - u[i - 1][j][k][2] * um1);
            rhs[i][j][k][3] = rhs[i][j][k][3] + dx4tx1 * (u[i + 1][j][k][3] - 2.0 * u[i][j][k][3] + u[i - 1][j][k][3]) + xxcon2 * (ws[i + 1][j][k] - 2.0 * ws[i][j][k] + ws[i - 1][j][k]) - tx2 * (u[i + 1][j][k][3] * up1 - u[i - 1][j][k][3] * um1);
            rhs[i][j][k][4] = rhs[i][j][k][4] + dx5tx1 * (u[i + 1][j][k][4] - 2.0 * u[i][j][k][4] + u[i - 1][j][k][4]) + xxcon3 * (qs[i + 1][j][k] - 2.0 * qs[i][j][k] + qs[i - 1][j][k]) + xxcon4 * (up1 * up1 - 2.0 * uijk * uijk + um1 * um1) + xxcon5 * (u[i + 1][j][k][4] * rho_i[i + 1][j][k] - 2.0 * u[i][j][k][4] * rho_i[i][j][k] + u[i - 1][j][k][4] * rho_i[i - 1][j][k]) - tx2 * ((c1 * u[i + 1][j][k][4] - c2 * square[i + 1][j][k]) * up1 - (c1 * u[i - 1][j][k][4] - c2 * square[i - 1][j][k]) * um1);
         }
      }
   }
   i = 1;
   #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, dssp, grid_points, u)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (5.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]);
         }
      }
   }
   i = 2;
   #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, dssp, grid_points, u)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (-4.0 * u[i - 1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(dssp, grid_points, u)
   for(i = 3; i < grid_points[0] - 3; i++) {
      #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, dssp, grid_points, u)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, dssp, grid_points, u)
         for(k = 1; k < grid_points[2] - 1; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i - 2][j][k][m] - 4.0 * u[i - 1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]);
            }
         }
      }
   }
   i = grid_points[0] - 3;
   #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, dssp, grid_points, u)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i - 2][j][k][m] - 4.0 * u[i - 1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m]);
         }
      }
   }
   i = grid_points[0] - 2;
   #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, dssp, grid_points, u)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i - 2][j][k][m] - 4. * u[i - 1][j][k][m] + 5.0 * u[i][j][k][m]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, vijk, vp1, vm1) firstprivate(dy1ty1, ty2, dy2ty1, yycon2, c2, dy3ty1, con43, dy4ty1, c1, dy5ty1, yycon3, yycon4, yycon5, grid_points, vs, u, us, square, ws, qs, rho_i)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, vijk, vp1, vm1) firstprivate(i, dy1ty1, ty2, dy2ty1, yycon2, c2, dy3ty1, con43, dy4ty1, c1, dy5ty1, yycon3, yycon4, yycon5, grid_points, vs, u, us, square, ws, qs, rho_i)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, vijk, vp1, vm1) firstprivate(i, j, dy1ty1, ty2, dy2ty1, yycon2, c2, dy3ty1, con43, dy4ty1, c1, dy5ty1, yycon3, yycon4, yycon5, grid_points, vs, u, us, square, ws, qs, rho_i)
         for(k = 1; k < grid_points[2] - 1; k++) {
            vijk = vs[i][j][k];
            vp1 = vs[i][j + 1][k];
            vm1 = vs[i][j - 1][k];
            rhs[i][j][k][0] = rhs[i][j][k][0] + dy1ty1 * (u[i][j + 1][k][0] - 2.0 * u[i][j][k][0] + u[i][j - 1][k][0]) - ty2 * (u[i][j + 1][k][2] - u[i][j - 1][k][2]);
            rhs[i][j][k][1] = rhs[i][j][k][1] + dy2ty1 * (u[i][j + 1][k][1] - 2.0 * u[i][j][k][1] + u[i][j - 1][k][1]) + yycon2 * (us[i][j + 1][k] - 2.0 * us[i][j][k] + us[i][j - 1][k]) - ty2 * (u[i][j + 1][k][1] * vp1 - u[i][j - 1][k][1] * vm1);
            rhs[i][j][k][2] = rhs[i][j][k][2] + dy3ty1 * (u[i][j + 1][k][2] - 2.0 * u[i][j][k][2] + u[i][j - 1][k][2]) + yycon2 * con43 * (vp1 - 2.0 * vijk + vm1) - ty2 * (u[i][j + 1][k][2] * vp1 - u[i][j - 1][k][2] * vm1 + (u[i][j + 1][k][4] - square[i][j + 1][k] - u[i][j - 1][k][4] + square[i][j - 1][k]) * c2);
            rhs[i][j][k][3] = rhs[i][j][k][3] + dy4ty1 * (u[i][j + 1][k][3] - 2.0 * u[i][j][k][3] + u[i][j - 1][k][3]) + yycon2 * (ws[i][j + 1][k] - 2.0 * ws[i][j][k] + ws[i][j - 1][k]) - ty2 * (u[i][j + 1][k][3] * vp1 - u[i][j - 1][k][3] * vm1);
            rhs[i][j][k][4] = rhs[i][j][k][4] + dy5ty1 * (u[i][j + 1][k][4] - 2.0 * u[i][j][k][4] + u[i][j - 1][k][4]) + yycon3 * (qs[i][j + 1][k] - 2.0 * qs[i][j][k] + qs[i][j - 1][k]) + yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk + vm1 * vm1) + yycon5 * (u[i][j + 1][k][4] * rho_i[i][j + 1][k] - 2.0 * u[i][j][k][4] * rho_i[i][j][k] + u[i][j - 1][k][4] * rho_i[i][j - 1][k]) - ty2 * ((c1 * u[i][j + 1][k][4] - c2 * square[i][j + 1][k]) * vp1 - (c1 * u[i][j - 1][k][4] - c2 * square[i][j - 1][k]) * vm1);
         }
      }
   }
   j = 1;
   #pragma omp parallel for default(shared) private(i, k, m) firstprivate(j, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(j, i, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (5.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m] + u[i][j + 2][k][m]);
         }
      }
   }
   j = 2;
   #pragma omp parallel for default(shared) private(i, k, m) firstprivate(j, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(j, i, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (-4.0 * u[i][j - 1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m] + u[i][j + 2][k][m]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, dssp, grid_points, u)
      for(j = 3; j < grid_points[1] - 3; j++) {
         #pragma omp parallel for default(shared) private(k, m) firstprivate(j, i, dssp, grid_points, u)
         for(k = 1; k < grid_points[2] - 1; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j - 2][k][m] - 4.0 * u[i][j - 1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m] + u[i][j + 2][k][m]);
            }
         }
      }
   }
   j = grid_points[1] - 3;
   #pragma omp parallel for default(shared) private(i, k, m) firstprivate(j, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(j, i, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j - 2][k][m] - 4.0 * u[i][j - 1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m]);
         }
      }
   }
   j = grid_points[1] - 2;
   #pragma omp parallel for default(shared) private(i, k, m) firstprivate(j, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(k, m) firstprivate(j, i, dssp, grid_points, u)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j - 2][k][m] - 4. * u[i][j - 1][k][m] + 5. * u[i][j][k][m]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, wijk, wp1, wm1) firstprivate(dz1tz1, tz2, dz2tz1, zzcon2, dz3tz1, c2, dz4tz1, con43, c1, dz5tz1, zzcon3, zzcon4, zzcon5, grid_points, ws, u, us, vs, square, qs, rho_i)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, wijk, wp1, wm1) firstprivate(i, dz1tz1, tz2, dz2tz1, zzcon2, dz3tz1, c2, dz4tz1, con43, c1, dz5tz1, zzcon3, zzcon4, zzcon5, grid_points, ws, u, us, vs, square, qs, rho_i)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, wijk, wp1, wm1) firstprivate(i, j, dz1tz1, tz2, dz2tz1, zzcon2, dz3tz1, c2, dz4tz1, con43, c1, dz5tz1, zzcon3, zzcon4, zzcon5, grid_points, ws, u, us, vs, square, qs, rho_i)
         for(k = 1; k < grid_points[2] - 1; k++) {
            wijk = ws[i][j][k];
            wp1 = ws[i][j][k + 1];
            wm1 = ws[i][j][k - 1];
            rhs[i][j][k][0] = rhs[i][j][k][0] + dz1tz1 * (u[i][j][k + 1][0] - 2.0 * u[i][j][k][0] + u[i][j][k - 1][0]) - tz2 * (u[i][j][k + 1][3] - u[i][j][k - 1][3]);
            rhs[i][j][k][1] = rhs[i][j][k][1] + dz2tz1 * (u[i][j][k + 1][1] - 2.0 * u[i][j][k][1] + u[i][j][k - 1][1]) + zzcon2 * (us[i][j][k + 1] - 2.0 * us[i][j][k] + us[i][j][k - 1]) - tz2 * (u[i][j][k + 1][1] * wp1 - u[i][j][k - 1][1] * wm1);
            rhs[i][j][k][2] = rhs[i][j][k][2] + dz3tz1 * (u[i][j][k + 1][2] - 2.0 * u[i][j][k][2] + u[i][j][k - 1][2]) + zzcon2 * (vs[i][j][k + 1] - 2.0 * vs[i][j][k] + vs[i][j][k - 1]) - tz2 * (u[i][j][k + 1][2] * wp1 - u[i][j][k - 1][2] * wm1);
            rhs[i][j][k][3] = rhs[i][j][k][3] + dz4tz1 * (u[i][j][k + 1][3] - 2.0 * u[i][j][k][3] + u[i][j][k - 1][3]) + zzcon2 * con43 * (wp1 - 2.0 * wijk + wm1) - tz2 * (u[i][j][k + 1][3] * wp1 - u[i][j][k - 1][3] * wm1 + (u[i][j][k + 1][4] - square[i][j][k + 1] - u[i][j][k - 1][4] + square[i][j][k - 1]) * c2);
            rhs[i][j][k][4] = rhs[i][j][k][4] + dz5tz1 * (u[i][j][k + 1][4] - 2.0 * u[i][j][k][4] + u[i][j][k - 1][4]) + zzcon3 * (qs[i][j][k + 1] - 2.0 * qs[i][j][k] + qs[i][j][k - 1]) + zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk + wm1 * wm1) + zzcon5 * (u[i][j][k + 1][4] * rho_i[i][j][k + 1] - 2.0 * u[i][j][k][4] * rho_i[i][j][k] + u[i][j][k - 1][4] * rho_i[i][j][k - 1]) - tz2 * ((c1 * u[i][j][k + 1][4] - c2 * square[i][j][k + 1]) * wp1 - (c1 * u[i][j][k - 1][4] - c2 * square[i][j][k - 1]) * wm1);
         }
      }
   }
   k = 1;
   #pragma omp parallel for default(shared) private(i, j, m) firstprivate(k, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, m) firstprivate(k, i, dssp, grid_points, u)
      for(j = 1; j < grid_points[1] - 1; j++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (5.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m] + u[i][j][k + 2][m]);
         }
      }
   }
   k = 2;
   #pragma omp parallel for default(shared) private(i, j, m) firstprivate(k, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, m) firstprivate(k, i, dssp, grid_points, u)
      for(j = 1; j < grid_points[1] - 1; j++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (-4.0 * u[i][j][k - 1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m] + u[i][j][k + 2][m]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i, dssp, grid_points, u)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, m) firstprivate(i, j, dssp, grid_points, u)
         for(k = 3; k < grid_points[2] - 3; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j][k - 2][m] - 4.0 * u[i][j][k - 1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m] + u[i][j][k + 2][m]);
            }
         }
      }
   }
   k = grid_points[2] - 3;
   #pragma omp parallel for default(shared) private(i, j, m) firstprivate(k, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, m) firstprivate(k, i, dssp, grid_points, u)
      for(j = 1; j < grid_points[1] - 1; j++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j][k - 2][m] - 4.0 * u[i][j][k - 1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m]);
         }
      }
   }
   k = grid_points[2] - 2;
   #pragma omp parallel for default(shared) private(i, j, m) firstprivate(k, dssp, grid_points, u)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, m) firstprivate(k, i, dssp, grid_points, u)
      for(j = 1; j < grid_points[1] - 1; j++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            rhs[i][j][k][m] = rhs[i][j][k][m] - dssp * (u[i][j][k - 2][m] - 4.0 * u[i][j][k - 1][m] + 5.0 * u[i][j][k][m]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(j, k, m, i) firstprivate(dt, grid_points)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k, m, i) firstprivate(j, dt, grid_points)
      for(k = 1; k < grid_points[2] - 1; k++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            #pragma omp parallel for default(shared) private(i) firstprivate(j, k, m, dt, grid_points)
            for(i = 1; i < grid_points[0] - 1; i++) {
               rhs[i][j][k][m] = rhs[i][j][k][m] * dt;
            }
         }
      }
   }
}

static void set_constants() {
   ce[0][0] = 2.0;
   ce[0][1] = 0.0;
   ce[0][2] = 0.0;
   ce[0][3] = 4.0;
   ce[0][4] = 5.0;
   ce[0][5] = 3.0;
   ce[0][6] = 0.5;
   ce[0][7] = 0.02;
   ce[0][8] = 0.01;
   ce[0][9] = 0.03;
   ce[0][10] = 0.5;
   ce[0][11] = 0.4;
   ce[0][12] = 0.3;
   ce[1][0] = 1.0;
   ce[1][1] = 0.0;
   ce[1][2] = 0.0;
   ce[1][3] = 0.0;
   ce[1][4] = 1.0;
   ce[1][5] = 2.0;
   ce[1][6] = 3.0;
   ce[1][7] = 0.01;
   ce[1][8] = 0.03;
   ce[1][9] = 0.02;
   ce[1][10] = 0.4;
   ce[1][11] = 0.3;
   ce[1][12] = 0.5;
   ce[2][0] = 2.0;
   ce[2][1] = 2.0;
   ce[2][2] = 0.0;
   ce[2][3] = 0.0;
   ce[2][4] = 0.0;
   ce[2][5] = 2.0;
   ce[2][6] = 3.0;
   ce[2][7] = 0.04;
   ce[2][8] = 0.03;
   ce[2][9] = 0.05;
   ce[2][10] = 0.3;
   ce[2][11] = 0.5;
   ce[2][12] = 0.4;
   ce[3][0] = 2.0;
   ce[3][1] = 2.0;
   ce[3][2] = 0.0;
   ce[3][3] = 0.0;
   ce[3][4] = 0.0;
   ce[3][5] = 2.0;
   ce[3][6] = 3.0;
   ce[3][7] = 0.03;
   ce[3][8] = 0.05;
   ce[3][9] = 0.04;
   ce[3][10] = 0.2;
   ce[3][11] = 0.1;
   ce[3][12] = 0.3;
   ce[4][0] = 5.0;
   ce[4][1] = 4.0;
   ce[4][2] = 3.0;
   ce[4][3] = 2.0;
   ce[4][4] = 0.1;
   ce[4][5] = 0.4;
   ce[4][6] = 0.3;
   ce[4][7] = 0.05;
   ce[4][8] = 0.04;
   ce[4][9] = 0.03;
   ce[4][10] = 0.1;
   ce[4][11] = 0.3;
   ce[4][12] = 0.2;
   c1 = 1.4;
   c2 = 0.4;
   c3 = 0.1;
   c4 = 1.0;
   c5 = 1.4;
   dnxm1 = 1.0 / (double) (grid_points[0] - 1);
   dnym1 = 1.0 / (double) (grid_points[1] - 1);
   dnzm1 = 1.0 / (double) (grid_points[2] - 1);
   c1c2 = c1 * c2;
   c1c5 = c1 * c5;
   c3c4 = c3 * c4;
   c1345 = c1c5 * c3c4;
   conz1 = (1.0 - c1c5);
   tx1 = 1.0 / (dnxm1 * dnxm1);
   tx2 = 1.0 / (2.0 * dnxm1);
   tx3 = 1.0 / dnxm1;
   ty1 = 1.0 / (dnym1 * dnym1);
   ty2 = 1.0 / (2.0 * dnym1);
   ty3 = 1.0 / dnym1;
   tz1 = 1.0 / (dnzm1 * dnzm1);
   tz2 = 1.0 / (2.0 * dnzm1);
   tz3 = 1.0 / dnzm1;
   dx1 = 0.75;
   dx2 = 0.75;
   dx3 = 0.75;
   dx4 = 0.75;
   dx5 = 0.75;
   dy1 = 0.75;
   dy2 = 0.75;
   dy3 = 0.75;
   dy4 = 0.75;
   dy5 = 0.75;
   dz1 = 1.0;
   dz2 = 1.0;
   dz3 = 1.0;
   dz4 = 1.0;
   dz5 = 1.0;
   dxmax = (((dx3) > (dx4)) ? (dx3) : (dx4));
   dymax = (((dy2) > (dy4)) ? (dy2) : (dy4));
   dzmax = (((dz2) > (dz3)) ? (dz2) : (dz3));
   dssp = 0.25 * (((dx1) > ((((dy1) > (dz1)) ? (dy1) : (dz1)))) ? (dx1) : ((((dy1) > (dz1)) ? (dy1) : (dz1))));
   c4dssp = 4.0 * dssp;
   c5dssp = 5.0 * dssp;
   dttx1 = dt * tx1;
   dttx2 = dt * tx2;
   dtty1 = dt * ty1;
   dtty2 = dt * ty2;
   dttz1 = dt * tz1;
   dttz2 = dt * tz2;
   c2dttx1 = 2.0 * dttx1;
   c2dtty1 = 2.0 * dtty1;
   c2dttz1 = 2.0 * dttz1;
   dtdssp = dt * dssp;
   comz1 = dtdssp;
   comz4 = 4.0 * dtdssp;
   comz5 = 5.0 * dtdssp;
   comz6 = 6.0 * dtdssp;
   c3c4tx3 = c3c4 * tx3;
   c3c4ty3 = c3c4 * ty3;
   c3c4tz3 = c3c4 * tz3;
   dx1tx1 = dx1 * tx1;
   dx2tx1 = dx2 * tx1;
   dx3tx1 = dx3 * tx1;
   dx4tx1 = dx4 * tx1;
   dx5tx1 = dx5 * tx1;
   dy1ty1 = dy1 * ty1;
   dy2ty1 = dy2 * ty1;
   dy3ty1 = dy3 * ty1;
   dy4ty1 = dy4 * ty1;
   dy5ty1 = dy5 * ty1;
   dz1tz1 = dz1 * tz1;
   dz2tz1 = dz2 * tz1;
   dz3tz1 = dz3 * tz1;
   dz4tz1 = dz4 * tz1;
   dz5tz1 = dz5 * tz1;
   c2iv = 2.5;
   con43 = 4.0 / 3.0;
   con16 = 1.0 / 6.0;
   xxcon1 = c3c4tx3 * con43 * tx3;
   xxcon2 = c3c4tx3 * tx3;
   xxcon3 = c3c4tx3 * conz1 * tx3;
   xxcon4 = c3c4tx3 * con16 * tx3;
   xxcon5 = c3c4tx3 * c1c5 * tx3;
   yycon1 = c3c4ty3 * con43 * ty3;
   yycon2 = c3c4ty3 * ty3;
   yycon3 = c3c4ty3 * conz1 * ty3;
   yycon4 = c3c4ty3 * con16 * ty3;
   yycon5 = c3c4ty3 * c1c5 * ty3;
   zzcon1 = c3c4tz3 * con43 * tz3;
   zzcon2 = c3c4tz3 * tz3;
   zzcon3 = c3c4tz3 * conz1 * tz3;
   zzcon4 = c3c4tz3 * con16 * tz3;
   zzcon5 = c3c4tz3 * c1c5 * tz3;
}

static void verify(int no_time_steps, char *class, boolean *verified) {
   double xcrref[5];
   double xceref[5];
   double xcrdif[5];
   double xcedif[5];
   double epsilon;
   double xce[5];
   double xcr[5];
   double dtref;
   int m;
   epsilon = 1.0e-08;
   error_norm(xce);
   compute_rhs();
   rhs_norm(xcr);
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      xcr[m] = xcr[m] / dt;
   }
   *class = 'U';
   *verified = 1;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      xcrref[m] = 1.0;
      xceref[m] = 1.0;
   }
   if(grid_points[0] == 12 && grid_points[1] == 12 && grid_points[2] == 12 && no_time_steps == 60) {
      *class = 'S';
      dtref = 1.0e-2;
      xcrref[0] = 1.7034283709541311e-01;
      xcrref[1] = 1.2975252070034097e-02;
      xcrref[2] = 3.2527926989486055e-02;
      xcrref[3] = 2.6436421275166801e-02;
      xcrref[4] = 1.9211784131744430e-01;
      xceref[0] = 4.9976913345811579e-04;
      xceref[1] = 4.5195666782961927e-05;
      xceref[2] = 7.3973765172921357e-05;
      xceref[3] = 7.3821238632439731e-05;
      xceref[4] = 8.9269630987491446e-04;
   }
   else if(grid_points[0] == 24 && grid_points[1] == 24 && grid_points[2] == 24 && no_time_steps == 200) {
         *class = 'W';
         dtref = 0.8e-3;
         xcrref[0] = 0.1125590409344e+03;
         xcrref[1] = 0.1180007595731e+02;
         xcrref[2] = 0.2710329767846e+02;
         xcrref[3] = 0.2469174937669e+02;
         xcrref[4] = 0.2638427874317e+03;
         xceref[0] = 0.4419655736008e+01;
         xceref[1] = 0.4638531260002e+00;
         xceref[2] = 0.1011551749967e+01;
         xceref[3] = 0.9235878729944e+00;
         xceref[4] = 0.1018045837718e+02;
      }
      else if(grid_points[0] == 64 && grid_points[1] == 64 && grid_points[2] == 64 && no_time_steps == 200) {
            *class = 'A';
            dtref = 0.8e-3;
            xcrref[0] = 1.0806346714637264e+02;
            xcrref[1] = 1.1319730901220813e+01;
            xcrref[2] = 2.5974354511582465e+01;
            xcrref[3] = 2.3665622544678910e+01;
            xcrref[4] = 2.5278963211748344e+02;
            xceref[0] = 4.2348416040525025e+00;
            xceref[1] = 4.4390282496995698e-01;
            xceref[2] = 9.6692480136345650e-01;
            xceref[3] = 8.8302063039765474e-01;
            xceref[4] = 9.7379901770829278e+00;
         }
         else if(grid_points[0] == 102 && grid_points[1] == 102 && grid_points[2] == 102 && no_time_steps == 200) {
               *class = 'B';
               dtref = 3.0e-4;
               xcrref[0] = 1.4233597229287254e+03;
               xcrref[1] = 9.9330522590150238e+01;
               xcrref[2] = 3.5646025644535285e+02;
               xcrref[3] = 3.2485447959084092e+02;
               xcrref[4] = 3.2707541254659363e+03;
               xceref[0] = 5.2969847140936856e+01;
               xceref[1] = 4.4632896115670668e+00;
               xceref[2] = 1.3122573342210174e+01;
               xceref[3] = 1.2006925323559144e+01;
               xceref[4] = 1.2459576151035986e+02;
            }
            else if(grid_points[0] == 162 && grid_points[1] == 162 && grid_points[2] == 162 && no_time_steps == 200) {
                  *class = 'C';
                  dtref = 1.0e-4;
                  xcrref[0] = 0.62398116551764615e+04;
                  xcrref[1] = 0.50793239190423964e+03;
                  xcrref[2] = 0.15423530093013596e+04;
                  xcrref[3] = 0.13302387929291190e+04;
                  xcrref[4] = 0.11604087428436455e+05;
                  xceref[0] = 0.16462008369091265e+03;
                  xceref[1] = 0.11497107903824313e+02;
                  xceref[2] = 0.41207446207461508e+02;
                  xceref[3] = 0.37087651059694167e+02;
                  xceref[4] = 0.36211053051841265e+03;
               }
               else {
                  *verified = 0;
               }
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      xcrdif[m] = fabs((xcr[m] - xcrref[m]) / xcrref[m]);
      xcedif[m] = fabs((xce[m] - xceref[m]) / xceref[m]);
   }
   if(*class != 'U') {
      printf(" Verification being performed for class %1c\n", *class);
      printf(" accuracy setting for epsilon = %20.13e\n", epsilon);
      if(fabs(dt - dtref) > epsilon) {
         *verified = 0;
         *class = 'U';
         printf(" DT does not match the reference value of %15.8e\n", dtref);
      }
   }
   else {
      printf(" Unknown class\n");
   }
   if(*class != 'U') {
      printf(" Comparison of RMS-norms of residual\n");
   }
   else {
      printf(" RMS-norms of residual\n");
   }
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      if(*class == 'U') {
         printf("          %2d%20.13e\n", m, xcr[m]);
      }
      else if(xcrdif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
         }
         else {
            printf("          %2d%20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
         }
   }
   if(*class != 'U') {
      printf(" Comparison of RMS-norms of solution error\n");
   }
   else {
      printf(" RMS-norms of solution error\n");
   }
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      if(*class == 'U') {
         printf("          %2d%20.13e\n", m, xce[m]);
      }
      else if(xcedif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
         }
         else {
            printf("          %2d%20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
         }
   }
   if(*class == 'U') {
      printf(" No reference values provided\n");
      printf(" No verification performed\n");
   }
   else if(*verified == 1) {
         printf(" Verification Successful\n");
      }
      else {
         printf(" Verification failed\n");
      }
}

static void x_solve() {
   lhsx();
   x_solve_cell();
   x_backsubstitute();
}

static void x_backsubstitute() {
   int i, j, k, m, n;
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess rhs	 use : RW
   ****************************************/
   for(i = grid_points[0] - 2; i >= 0; i--) {
      #pragma omp parallel for default(shared) private(j, k, m, n) firstprivate(i, grid_points, lhs)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k, m, n) firstprivate(i, j, grid_points, lhs)
         for(k = 1; k < grid_points[2] - 1; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(n = 0; n < 5; n++) {
                  rhs[i][j][k][m] = rhs[i][j][k][m] - lhs[i][j][k][2][m][n] * rhs[i + 1][j][k][n];
               }
            }
         }
      }
   }
}

static void x_solve_cell() {
   int i, j, k, isize;
   isize = grid_points[0] - 1;
   #pragma omp parallel for default(shared) private(j, k) firstprivate(grid_points)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(j, grid_points)
      for(k = 1; k < grid_points[2] - 1; k++) {
         binvcrhs(lhs[0][j][k][1], lhs[0][j][k][2], rhs[0][j][k]);
      }
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess lhs	 use : RW
   		unsolved dependency for arrayAccess rhs	 use : RW
   ****************************************/
   for(i = 1; i < isize; i++) {
      #pragma omp parallel for default(shared) private(j, k) firstprivate(i, grid_points)
      for(j = 1; j < grid_points[1] - 1; j++) {
         #pragma omp parallel for default(shared) private(k) firstprivate(i, j, grid_points)
         for(k = 1; k < grid_points[2] - 1; k++) {
            matvec_sub(lhs[i][j][k][0], rhs[i - 1][j][k], rhs[i][j][k]);
            matmul_sub(lhs[i][j][k][0], lhs[i - 1][j][k][2], lhs[i][j][k][1]);
            binvcrhs(lhs[i][j][k][1], lhs[i][j][k][2], rhs[i][j][k]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(j, k) firstprivate(isize, i, grid_points)
   for(j = 1; j < grid_points[1] - 1; j++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(isize, j, i, grid_points)
      for(k = 1; k < grid_points[2] - 1; k++) {
         matvec_sub(lhs[isize][j][k][0], rhs[isize - 1][j][k], rhs[isize][j][k]);
         matmul_sub(lhs[isize][j][k][0], lhs[isize - 1][j][k][2], lhs[isize][j][k][1]);
         binvrhs(lhs[i][j][k][1], rhs[i][j][k]);
      }
   }
}

static void matvec_sub(double ablock[5][5], double avec[5], double bvec[5]) {
   int i;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(i = 0; i < 5; i++) {
      bvec[i] = bvec[i] - ablock[i][0] * avec[0] - ablock[i][1] * avec[1] - ablock[i][2] * avec[2] - ablock[i][3] * avec[3] - ablock[i][4] * avec[4];
   }
}

static void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5]) {
   int j;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(j = 0; j < 5; j++) {
      cblock[0][j] = cblock[0][j] - ablock[0][0] * bblock[0][j] - ablock[0][1] * bblock[1][j] - ablock[0][2] * bblock[2][j] - ablock[0][3] * bblock[3][j] - ablock[0][4] * bblock[4][j];
      cblock[1][j] = cblock[1][j] - ablock[1][0] * bblock[0][j] - ablock[1][1] * bblock[1][j] - ablock[1][2] * bblock[2][j] - ablock[1][3] * bblock[3][j] - ablock[1][4] * bblock[4][j];
      cblock[2][j] = cblock[2][j] - ablock[2][0] * bblock[0][j] - ablock[2][1] * bblock[1][j] - ablock[2][2] * bblock[2][j] - ablock[2][3] * bblock[3][j] - ablock[2][4] * bblock[4][j];
      cblock[3][j] = cblock[3][j] - ablock[3][0] * bblock[0][j] - ablock[3][1] * bblock[1][j] - ablock[3][2] * bblock[2][j] - ablock[3][3] * bblock[3][j] - ablock[3][4] * bblock[4][j];
      cblock[4][j] = cblock[4][j] - ablock[4][0] * bblock[0][j] - ablock[4][1] * bblock[1][j] - ablock[4][2] * bblock[2][j] - ablock[4][3] * bblock[3][j] - ablock[4][4] * bblock[4][j];
   }
}

static void binvcrhs(double lhs[5][5], double c[5][5], double r[5]) {
   double pivot, coeff;
   pivot = 1.00 / lhs[0][0];
   lhs[0][1] = lhs[0][1] * pivot;
   lhs[0][2] = lhs[0][2] * pivot;
   lhs[0][3] = lhs[0][3] * pivot;
   lhs[0][4] = lhs[0][4] * pivot;
   c[0][0] = c[0][0] * pivot;
   c[0][1] = c[0][1] * pivot;
   c[0][2] = c[0][2] * pivot;
   c[0][3] = c[0][3] * pivot;
   c[0][4] = c[0][4] * pivot;
   r[0] = r[0] * pivot;
   coeff = lhs[1][0];
   lhs[1][1] = lhs[1][1] - coeff * lhs[0][1];
   lhs[1][2] = lhs[1][2] - coeff * lhs[0][2];
   lhs[1][3] = lhs[1][3] - coeff * lhs[0][3];
   lhs[1][4] = lhs[1][4] - coeff * lhs[0][4];
   c[1][0] = c[1][0] - coeff * c[0][0];
   c[1][1] = c[1][1] - coeff * c[0][1];
   c[1][2] = c[1][2] - coeff * c[0][2];
   c[1][3] = c[1][3] - coeff * c[0][3];
   c[1][4] = c[1][4] - coeff * c[0][4];
   r[1] = r[1] - coeff * r[0];
   coeff = lhs[2][0];
   lhs[2][1] = lhs[2][1] - coeff * lhs[0][1];
   lhs[2][2] = lhs[2][2] - coeff * lhs[0][2];
   lhs[2][3] = lhs[2][3] - coeff * lhs[0][3];
   lhs[2][4] = lhs[2][4] - coeff * lhs[0][4];
   c[2][0] = c[2][0] - coeff * c[0][0];
   c[2][1] = c[2][1] - coeff * c[0][1];
   c[2][2] = c[2][2] - coeff * c[0][2];
   c[2][3] = c[2][3] - coeff * c[0][3];
   c[2][4] = c[2][4] - coeff * c[0][4];
   r[2] = r[2] - coeff * r[0];
   coeff = lhs[3][0];
   lhs[3][1] = lhs[3][1] - coeff * lhs[0][1];
   lhs[3][2] = lhs[3][2] - coeff * lhs[0][2];
   lhs[3][3] = lhs[3][3] - coeff * lhs[0][3];
   lhs[3][4] = lhs[3][4] - coeff * lhs[0][4];
   c[3][0] = c[3][0] - coeff * c[0][0];
   c[3][1] = c[3][1] - coeff * c[0][1];
   c[3][2] = c[3][2] - coeff * c[0][2];
   c[3][3] = c[3][3] - coeff * c[0][3];
   c[3][4] = c[3][4] - coeff * c[0][4];
   r[3] = r[3] - coeff * r[0];
   coeff = lhs[4][0];
   lhs[4][1] = lhs[4][1] - coeff * lhs[0][1];
   lhs[4][2] = lhs[4][2] - coeff * lhs[0][2];
   lhs[4][3] = lhs[4][3] - coeff * lhs[0][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[0][4];
   c[4][0] = c[4][0] - coeff * c[0][0];
   c[4][1] = c[4][1] - coeff * c[0][1];
   c[4][2] = c[4][2] - coeff * c[0][2];
   c[4][3] = c[4][3] - coeff * c[0][3];
   c[4][4] = c[4][4] - coeff * c[0][4];
   r[4] = r[4] - coeff * r[0];
   pivot = 1.00 / lhs[1][1];
   lhs[1][2] = lhs[1][2] * pivot;
   lhs[1][3] = lhs[1][3] * pivot;
   lhs[1][4] = lhs[1][4] * pivot;
   c[1][0] = c[1][0] * pivot;
   c[1][1] = c[1][1] * pivot;
   c[1][2] = c[1][2] * pivot;
   c[1][3] = c[1][3] * pivot;
   c[1][4] = c[1][4] * pivot;
   r[1] = r[1] * pivot;
   coeff = lhs[0][1];
   lhs[0][2] = lhs[0][2] - coeff * lhs[1][2];
   lhs[0][3] = lhs[0][3] - coeff * lhs[1][3];
   lhs[0][4] = lhs[0][4] - coeff * lhs[1][4];
   c[0][0] = c[0][0] - coeff * c[1][0];
   c[0][1] = c[0][1] - coeff * c[1][1];
   c[0][2] = c[0][2] - coeff * c[1][2];
   c[0][3] = c[0][3] - coeff * c[1][3];
   c[0][4] = c[0][4] - coeff * c[1][4];
   r[0] = r[0] - coeff * r[1];
   coeff = lhs[2][1];
   lhs[2][2] = lhs[2][2] - coeff * lhs[1][2];
   lhs[2][3] = lhs[2][3] - coeff * lhs[1][3];
   lhs[2][4] = lhs[2][4] - coeff * lhs[1][4];
   c[2][0] = c[2][0] - coeff * c[1][0];
   c[2][1] = c[2][1] - coeff * c[1][1];
   c[2][2] = c[2][2] - coeff * c[1][2];
   c[2][3] = c[2][3] - coeff * c[1][3];
   c[2][4] = c[2][4] - coeff * c[1][4];
   r[2] = r[2] - coeff * r[1];
   coeff = lhs[3][1];
   lhs[3][2] = lhs[3][2] - coeff * lhs[1][2];
   lhs[3][3] = lhs[3][3] - coeff * lhs[1][3];
   lhs[3][4] = lhs[3][4] - coeff * lhs[1][4];
   c[3][0] = c[3][0] - coeff * c[1][0];
   c[3][1] = c[3][1] - coeff * c[1][1];
   c[3][2] = c[3][2] - coeff * c[1][2];
   c[3][3] = c[3][3] - coeff * c[1][3];
   c[3][4] = c[3][4] - coeff * c[1][4];
   r[3] = r[3] - coeff * r[1];
   coeff = lhs[4][1];
   lhs[4][2] = lhs[4][2] - coeff * lhs[1][2];
   lhs[4][3] = lhs[4][3] - coeff * lhs[1][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[1][4];
   c[4][0] = c[4][0] - coeff * c[1][0];
   c[4][1] = c[4][1] - coeff * c[1][1];
   c[4][2] = c[4][2] - coeff * c[1][2];
   c[4][3] = c[4][3] - coeff * c[1][3];
   c[4][4] = c[4][4] - coeff * c[1][4];
   r[4] = r[4] - coeff * r[1];
   pivot = 1.00 / lhs[2][2];
   lhs[2][3] = lhs[2][3] * pivot;
   lhs[2][4] = lhs[2][4] * pivot;
   c[2][0] = c[2][0] * pivot;
   c[2][1] = c[2][1] * pivot;
   c[2][2] = c[2][2] * pivot;
   c[2][3] = c[2][3] * pivot;
   c[2][4] = c[2][4] * pivot;
   r[2] = r[2] * pivot;
   coeff = lhs[0][2];
   lhs[0][3] = lhs[0][3] - coeff * lhs[2][3];
   lhs[0][4] = lhs[0][4] - coeff * lhs[2][4];
   c[0][0] = c[0][0] - coeff * c[2][0];
   c[0][1] = c[0][1] - coeff * c[2][1];
   c[0][2] = c[0][2] - coeff * c[2][2];
   c[0][3] = c[0][3] - coeff * c[2][3];
   c[0][4] = c[0][4] - coeff * c[2][4];
   r[0] = r[0] - coeff * r[2];
   coeff = lhs[1][2];
   lhs[1][3] = lhs[1][3] - coeff * lhs[2][3];
   lhs[1][4] = lhs[1][4] - coeff * lhs[2][4];
   c[1][0] = c[1][0] - coeff * c[2][0];
   c[1][1] = c[1][1] - coeff * c[2][1];
   c[1][2] = c[1][2] - coeff * c[2][2];
   c[1][3] = c[1][3] - coeff * c[2][3];
   c[1][4] = c[1][4] - coeff * c[2][4];
   r[1] = r[1] - coeff * r[2];
   coeff = lhs[3][2];
   lhs[3][3] = lhs[3][3] - coeff * lhs[2][3];
   lhs[3][4] = lhs[3][4] - coeff * lhs[2][4];
   c[3][0] = c[3][0] - coeff * c[2][0];
   c[3][1] = c[3][1] - coeff * c[2][1];
   c[3][2] = c[3][2] - coeff * c[2][2];
   c[3][3] = c[3][3] - coeff * c[2][3];
   c[3][4] = c[3][4] - coeff * c[2][4];
   r[3] = r[3] - coeff * r[2];
   coeff = lhs[4][2];
   lhs[4][3] = lhs[4][3] - coeff * lhs[2][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[2][4];
   c[4][0] = c[4][0] - coeff * c[2][0];
   c[4][1] = c[4][1] - coeff * c[2][1];
   c[4][2] = c[4][2] - coeff * c[2][2];
   c[4][3] = c[4][3] - coeff * c[2][3];
   c[4][4] = c[4][4] - coeff * c[2][4];
   r[4] = r[4] - coeff * r[2];
   pivot = 1.00 / lhs[3][3];
   lhs[3][4] = lhs[3][4] * pivot;
   c[3][0] = c[3][0] * pivot;
   c[3][1] = c[3][1] * pivot;
   c[3][2] = c[3][2] * pivot;
   c[3][3] = c[3][3] * pivot;
   c[3][4] = c[3][4] * pivot;
   r[3] = r[3] * pivot;
   coeff = lhs[0][3];
   lhs[0][4] = lhs[0][4] - coeff * lhs[3][4];
   c[0][0] = c[0][0] - coeff * c[3][0];
   c[0][1] = c[0][1] - coeff * c[3][1];
   c[0][2] = c[0][2] - coeff * c[3][2];
   c[0][3] = c[0][3] - coeff * c[3][3];
   c[0][4] = c[0][4] - coeff * c[3][4];
   r[0] = r[0] - coeff * r[3];
   coeff = lhs[1][3];
   lhs[1][4] = lhs[1][4] - coeff * lhs[3][4];
   c[1][0] = c[1][0] - coeff * c[3][0];
   c[1][1] = c[1][1] - coeff * c[3][1];
   c[1][2] = c[1][2] - coeff * c[3][2];
   c[1][3] = c[1][3] - coeff * c[3][3];
   c[1][4] = c[1][4] - coeff * c[3][4];
   r[1] = r[1] - coeff * r[3];
   coeff = lhs[2][3];
   lhs[2][4] = lhs[2][4] - coeff * lhs[3][4];
   c[2][0] = c[2][0] - coeff * c[3][0];
   c[2][1] = c[2][1] - coeff * c[3][1];
   c[2][2] = c[2][2] - coeff * c[3][2];
   c[2][3] = c[2][3] - coeff * c[3][3];
   c[2][4] = c[2][4] - coeff * c[3][4];
   r[2] = r[2] - coeff * r[3];
   coeff = lhs[4][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[3][4];
   c[4][0] = c[4][0] - coeff * c[3][0];
   c[4][1] = c[4][1] - coeff * c[3][1];
   c[4][2] = c[4][2] - coeff * c[3][2];
   c[4][3] = c[4][3] - coeff * c[3][3];
   c[4][4] = c[4][4] - coeff * c[3][4];
   r[4] = r[4] - coeff * r[3];
   pivot = 1.00 / lhs[4][4];
   c[4][0] = c[4][0] * pivot;
   c[4][1] = c[4][1] * pivot;
   c[4][2] = c[4][2] * pivot;
   c[4][3] = c[4][3] * pivot;
   c[4][4] = c[4][4] * pivot;
   r[4] = r[4] * pivot;
   coeff = lhs[0][4];
   c[0][0] = c[0][0] - coeff * c[4][0];
   c[0][1] = c[0][1] - coeff * c[4][1];
   c[0][2] = c[0][2] - coeff * c[4][2];
   c[0][3] = c[0][3] - coeff * c[4][3];
   c[0][4] = c[0][4] - coeff * c[4][4];
   r[0] = r[0] - coeff * r[4];
   coeff = lhs[1][4];
   c[1][0] = c[1][0] - coeff * c[4][0];
   c[1][1] = c[1][1] - coeff * c[4][1];
   c[1][2] = c[1][2] - coeff * c[4][2];
   c[1][3] = c[1][3] - coeff * c[4][3];
   c[1][4] = c[1][4] - coeff * c[4][4];
   r[1] = r[1] - coeff * r[4];
   coeff = lhs[2][4];
   c[2][0] = c[2][0] - coeff * c[4][0];
   c[2][1] = c[2][1] - coeff * c[4][1];
   c[2][2] = c[2][2] - coeff * c[4][2];
   c[2][3] = c[2][3] - coeff * c[4][3];
   c[2][4] = c[2][4] - coeff * c[4][4];
   r[2] = r[2] - coeff * r[4];
   coeff = lhs[3][4];
   c[3][0] = c[3][0] - coeff * c[4][0];
   c[3][1] = c[3][1] - coeff * c[4][1];
   c[3][2] = c[3][2] - coeff * c[4][2];
   c[3][3] = c[3][3] - coeff * c[4][3];
   c[3][4] = c[3][4] - coeff * c[4][4];
   r[3] = r[3] - coeff * r[4];
}

static void binvrhs(double lhs[5][5], double r[5]) {
   double pivot, coeff;
   pivot = 1.00 / lhs[0][0];
   lhs[0][1] = lhs[0][1] * pivot;
   lhs[0][2] = lhs[0][2] * pivot;
   lhs[0][3] = lhs[0][3] * pivot;
   lhs[0][4] = lhs[0][4] * pivot;
   r[0] = r[0] * pivot;
   coeff = lhs[1][0];
   lhs[1][1] = lhs[1][1] - coeff * lhs[0][1];
   lhs[1][2] = lhs[1][2] - coeff * lhs[0][2];
   lhs[1][3] = lhs[1][3] - coeff * lhs[0][3];
   lhs[1][4] = lhs[1][4] - coeff * lhs[0][4];
   r[1] = r[1] - coeff * r[0];
   coeff = lhs[2][0];
   lhs[2][1] = lhs[2][1] - coeff * lhs[0][1];
   lhs[2][2] = lhs[2][2] - coeff * lhs[0][2];
   lhs[2][3] = lhs[2][3] - coeff * lhs[0][3];
   lhs[2][4] = lhs[2][4] - coeff * lhs[0][4];
   r[2] = r[2] - coeff * r[0];
   coeff = lhs[3][0];
   lhs[3][1] = lhs[3][1] - coeff * lhs[0][1];
   lhs[3][2] = lhs[3][2] - coeff * lhs[0][2];
   lhs[3][3] = lhs[3][3] - coeff * lhs[0][3];
   lhs[3][4] = lhs[3][4] - coeff * lhs[0][4];
   r[3] = r[3] - coeff * r[0];
   coeff = lhs[4][0];
   lhs[4][1] = lhs[4][1] - coeff * lhs[0][1];
   lhs[4][2] = lhs[4][2] - coeff * lhs[0][2];
   lhs[4][3] = lhs[4][3] - coeff * lhs[0][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[0][4];
   r[4] = r[4] - coeff * r[0];
   pivot = 1.00 / lhs[1][1];
   lhs[1][2] = lhs[1][2] * pivot;
   lhs[1][3] = lhs[1][3] * pivot;
   lhs[1][4] = lhs[1][4] * pivot;
   r[1] = r[1] * pivot;
   coeff = lhs[0][1];
   lhs[0][2] = lhs[0][2] - coeff * lhs[1][2];
   lhs[0][3] = lhs[0][3] - coeff * lhs[1][3];
   lhs[0][4] = lhs[0][4] - coeff * lhs[1][4];
   r[0] = r[0] - coeff * r[1];
   coeff = lhs[2][1];
   lhs[2][2] = lhs[2][2] - coeff * lhs[1][2];
   lhs[2][3] = lhs[2][3] - coeff * lhs[1][3];
   lhs[2][4] = lhs[2][4] - coeff * lhs[1][4];
   r[2] = r[2] - coeff * r[1];
   coeff = lhs[3][1];
   lhs[3][2] = lhs[3][2] - coeff * lhs[1][2];
   lhs[3][3] = lhs[3][3] - coeff * lhs[1][3];
   lhs[3][4] = lhs[3][4] - coeff * lhs[1][4];
   r[3] = r[3] - coeff * r[1];
   coeff = lhs[4][1];
   lhs[4][2] = lhs[4][2] - coeff * lhs[1][2];
   lhs[4][3] = lhs[4][3] - coeff * lhs[1][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[1][4];
   r[4] = r[4] - coeff * r[1];
   pivot = 1.00 / lhs[2][2];
   lhs[2][3] = lhs[2][3] * pivot;
   lhs[2][4] = lhs[2][4] * pivot;
   r[2] = r[2] * pivot;
   coeff = lhs[0][2];
   lhs[0][3] = lhs[0][3] - coeff * lhs[2][3];
   lhs[0][4] = lhs[0][4] - coeff * lhs[2][4];
   r[0] = r[0] - coeff * r[2];
   coeff = lhs[1][2];
   lhs[1][3] = lhs[1][3] - coeff * lhs[2][3];
   lhs[1][4] = lhs[1][4] - coeff * lhs[2][4];
   r[1] = r[1] - coeff * r[2];
   coeff = lhs[3][2];
   lhs[3][3] = lhs[3][3] - coeff * lhs[2][3];
   lhs[3][4] = lhs[3][4] - coeff * lhs[2][4];
   r[3] = r[3] - coeff * r[2];
   coeff = lhs[4][2];
   lhs[4][3] = lhs[4][3] - coeff * lhs[2][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[2][4];
   r[4] = r[4] - coeff * r[2];
   pivot = 1.00 / lhs[3][3];
   lhs[3][4] = lhs[3][4] * pivot;
   r[3] = r[3] * pivot;
   coeff = lhs[0][3];
   lhs[0][4] = lhs[0][4] - coeff * lhs[3][4];
   r[0] = r[0] - coeff * r[3];
   coeff = lhs[1][3];
   lhs[1][4] = lhs[1][4] - coeff * lhs[3][4];
   r[1] = r[1] - coeff * r[3];
   coeff = lhs[2][3];
   lhs[2][4] = lhs[2][4] - coeff * lhs[3][4];
   r[2] = r[2] - coeff * r[3];
   coeff = lhs[4][3];
   lhs[4][4] = lhs[4][4] - coeff * lhs[3][4];
   r[4] = r[4] - coeff * r[3];
   pivot = 1.00 / lhs[4][4];
   r[4] = r[4] * pivot;
   coeff = lhs[0][4];
   r[0] = r[0] - coeff * r[4];
   coeff = lhs[1][4];
   r[1] = r[1] - coeff * r[4];
   coeff = lhs[2][4];
   r[2] = r[2] - coeff * r[4];
   coeff = lhs[3][4];
   r[3] = r[3] - coeff * r[4];
}

static void y_solve() {
   lhsy();
   y_solve_cell();
   y_backsubstitute();
}

static void y_backsubstitute() {
   int i, j, k, m, n;
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess rhs	 use : RW
   ****************************************/
   for(j = grid_points[1] - 2; j >= 0; j--) {
      #pragma omp parallel for default(shared) private(i, k, m, n) firstprivate(j, grid_points, lhs)
      for(i = 1; i < grid_points[0] - 1; i++) {
         #pragma omp parallel for default(shared) private(k, m, n) firstprivate(j, i, grid_points, lhs)
         for(k = 1; k < grid_points[2] - 1; k++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(n = 0; n < 5; n++) {
                  rhs[i][j][k][m] = rhs[i][j][k][m] - lhs[i][j][k][2][m][n] * rhs[i][j + 1][k][n];
               }
            }
         }
      }
   }
}

static void y_solve_cell() {
   int i, j, k, jsize;
   jsize = grid_points[1] - 1;
   #pragma omp parallel for default(shared) private(i, k) firstprivate(grid_points)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(i, grid_points)
      for(k = 1; k < grid_points[2] - 1; k++) {
         binvcrhs(lhs[i][0][k][1], lhs[i][0][k][2], rhs[i][0][k]);
      }
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess lhs	 use : RW
   		unsolved dependency for arrayAccess rhs	 use : RW
   ****************************************/
   for(j = 1; j < jsize; j++) {
      #pragma omp parallel for default(shared) private(i, k) firstprivate(j, grid_points)
      for(i = 1; i < grid_points[0] - 1; i++) {
         #pragma omp parallel for default(shared) private(k) firstprivate(j, i, grid_points)
         for(k = 1; k < grid_points[2] - 1; k++) {
            matvec_sub(lhs[i][j][k][0], rhs[i][j - 1][k], rhs[i][j][k]);
            matmul_sub(lhs[i][j][k][0], lhs[i][j - 1][k][2], lhs[i][j][k][1]);
            binvcrhs(lhs[i][j][k][1], lhs[i][j][k][2], rhs[i][j][k]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, k) firstprivate(jsize, grid_points)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(jsize, i, grid_points)
      for(k = 1; k < grid_points[2] - 1; k++) {
         matvec_sub(lhs[i][jsize][k][0], rhs[i][jsize - 1][k], rhs[i][jsize][k]);
         matmul_sub(lhs[i][jsize][k][0], lhs[i][jsize - 1][k][2], lhs[i][jsize][k][1]);
         binvrhs(lhs[i][jsize][k][1], rhs[i][jsize][k]);
      }
   }
}

static void z_solve() {
   lhsz();
   z_solve_cell();
   z_backsubstitute();
}

static void z_backsubstitute() {
   int i, j, k, m, n;
   #pragma omp parallel for default(shared) private(i, j, k, m, n) firstprivate(grid_points, lhs)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j, k, m, n) firstprivate(i, grid_points, lhs)
      for(j = 1; j < grid_points[1] - 1; j++) {
         /*************** Clava msgError **************
         		unsolved dependency for arrayAccess rhs	 use : RW
         ****************************************/
         for(k = grid_points[2] - 2; k >= 0; k--) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(n = 0; n < 5; n++) {
                  rhs[i][j][k][m] = rhs[i][j][k][m] - lhs[i][j][k][2][m][n] * rhs[i][j][k + 1][n];
               }
            }
         }
      }
   }
}

static void z_solve_cell() {
   int i, j, k, ksize;
   ksize = grid_points[2] - 1;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(grid_points)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j) firstprivate(i, grid_points)
      for(j = 1; j < grid_points[1] - 1; j++) {
         binvcrhs(lhs[i][j][0][1], lhs[i][j][0][2], rhs[i][j][0]);
      }
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess lhs	 use : RW
   		unsolved dependency for arrayAccess rhs	 use : RW
   ****************************************/
   for(k = 1; k < ksize; k++) {
      #pragma omp parallel for default(shared) private(i, j) firstprivate(k, grid_points)
      for(i = 1; i < grid_points[0] - 1; i++) {
         #pragma omp parallel for default(shared) private(j) firstprivate(k, i, grid_points)
         for(j = 1; j < grid_points[1] - 1; j++) {
            matvec_sub(lhs[i][j][k][0], rhs[i][j][k - 1], rhs[i][j][k]);
            matmul_sub(lhs[i][j][k][0], lhs[i][j][k - 1][2], lhs[i][j][k][1]);
            binvcrhs(lhs[i][j][k][1], lhs[i][j][k][2], rhs[i][j][k]);
         }
      }
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(ksize, grid_points)
   for(i = 1; i < grid_points[0] - 1; i++) {
      #pragma omp parallel for default(shared) private(j) firstprivate(ksize, i, grid_points)
      for(j = 1; j < grid_points[1] - 1; j++) {
         matvec_sub(lhs[i][j][ksize][0], rhs[i][j][ksize - 1], rhs[i][j][ksize]);
         matmul_sub(lhs[i][j][ksize][0], lhs[i][j][ksize - 1][2], lhs[i][j][ksize][1]);
         binvrhs(lhs[i][j][ksize][1], rhs[i][j][ksize]);
      }
   }
}
