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

/*full problem size*/

/*number of iterations and how often to print the norm*/

static int nx;
static int ny;
static int nz;
static int nx0;
static int ny0;
static int nz0;
static int ist;
static int iend;
static int jst;
static int jend;
static int ii1;
static int ii2;
static int ji1;
static int ji2;
static int ki1;
static int ki2;
static double dxi;
static double deta;
static double dzeta;
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
static double u[64][65][65][5];
static double rsd[64][65][65][5];
static double frct[64][65][65][5];
static double flux[64][65][65][5];
static int ipr;
static int inorm;
static int itmax;
static int invert;
static double dt;
static double omega;
static double tolrsd[5];
static double rsdnm[5];
static double errnm[5];
static double frc;
static double ttotal;
/*common /cjac/*/

static double a[64][64][5][5];
static double b[64][64][5][5];
static double c[64][64][5][5];
static double d[64][64][5][5];
static double ce[5][13];
static double maxtime;
static void blts(int nx, int ny, int nz, int k, double omega, double v[64][65][65][5], double ldz[64][64][5][5], double ldy[64][64][5][5], double ldx[64][64][5][5], double d[64][64][5][5], int ist, int iend, int jst, int jend, int nx0, int ny0);
static void buts(int nx, int ny, int nz, int k, double omega, double v[64][65][65][5], double tv[64][64][5], double d[64][64][5][5], double udx[64][64][5][5], double udy[64][64][5][5], double udz[64][64][5][5], int ist, int iend, int jst, int jend, int nx0, int ny0);
static void domain();
static void erhs();
static void error();
static void exact(int i, int j, int k, double u000ijk[5]);
static void jacld(int k);
static void jacu(int k);
static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend, double v[64][65][65][5], double sum[5]);
static void pintgr();
static void read_input();
static void rhs();
static void setbv();
static void setcoeff();
static void setiv();
static void ssor();
static void verify(double xcr[5], double xce[5], double xci, char *class, boolean *verified);
int main(int argc, char **argv) {
   char class;
   boolean verified;
   double mflops;
   int nthreads = 1;
   read_input();
   domain();
   setcoeff();
   setbv();
   setiv();
   erhs();
   {
   }
   ssor();
   error();
   pintgr();
   verify(rsdnm, errnm, frc, &class, &verified);
   mflops = (double) itmax * (1984.77 * (double) nx0 * (double) ny0 * (double) nz0 - 10923.3 * (((double) (nx0 + ny0 + nz0) / 3.0) * ((double) (nx0 + ny0 + nz0) / 3.0)) + 27770.9 * (double) (nx0 + ny0 + nz0) / 3.0 - 144010.0) / (maxtime * 1000000.0);
   c_print_results("LU", class, nx0, ny0, nz0, itmax, nthreads, maxtime, mflops, "          floating point", verified, "3.0 structured", "25 Nov 2023", "gcc-12", "gcc-12", "-fopenmp -lm", "-I../common -fopenmp", "-lm", "(none)", "(none)");
}

static void blts(int nx, int ny, int nz, int k, double omega, double v[64][65][65][5], double ldz[64][64][5][5], double ldy[64][64][5][5], double ldx[64][64][5][5], double d[64][64][5][5], int ist, int iend, int jst, int jend, int nx0, int ny0) {
   int i, j, m;
   double tmp, tmp1;
   double tmat[5][5];
   #pragma omp parallel for default(shared) private(i, j, m) firstprivate(ist, iend, jst, jend, k, omega, ldz)
   for(i = ist; i <= iend; i++) {
      #pragma omp parallel for default(shared) private(j, m) firstprivate(jst, jend, k, i, omega, ldz)
      for(j = jst; j <= jend; j++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            v[i][j][k][m] = v[i][j][k][m] - omega * (ldz[i][j][m][0] * v[i][j][k - 1][0] + ldz[i][j][m][1] * v[i][j][k - 1][1] + ldz[i][j][m][2] * v[i][j][k - 1][2] + ldz[i][j][m][3] * v[i][j][k - 1][3] + ldz[i][j][m][4] * v[i][j][k - 1][4]);
         }
      }
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess v	 use : RW
   ****************************************/
   for(i = ist; i <= iend; i++) {
      /*************** Clava msgError **************
      		unsolved dependency for arrayAccess v	 use : RW
      ****************************************/
      for(j = jst; j <= jend; j++) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            v[i][j][k][m] = v[i][j][k][m] - omega * (ldy[i][j][m][0] * v[i][j - 1][k][0] + ldx[i][j][m][0] * v[i - 1][j][k][0] + ldy[i][j][m][1] * v[i][j - 1][k][1] + ldx[i][j][m][1] * v[i - 1][j][k][1] + ldy[i][j][m][2] * v[i][j - 1][k][2] + ldx[i][j][m][2] * v[i - 1][j][k][2] + ldy[i][j][m][3] * v[i][j - 1][k][3] + ldx[i][j][m][3] * v[i - 1][j][k][3] + ldy[i][j][m][4] * v[i][j - 1][k][4] + ldx[i][j][m][4] * v[i - 1][j][k][4]);
         }
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            tmat[m][0] = d[i][j][m][0];
            tmat[m][1] = d[i][j][m][1];
            tmat[m][2] = d[i][j][m][2];
            tmat[m][3] = d[i][j][m][3];
            tmat[m][4] = d[i][j][m][4];
         }
         tmp1 = 1.0 / tmat[0][0];
         tmp = tmp1 * tmat[1][0];
         tmat[1][1] = tmat[1][1] - tmp * tmat[0][1];
         tmat[1][2] = tmat[1][2] - tmp * tmat[0][2];
         tmat[1][3] = tmat[1][3] - tmp * tmat[0][3];
         tmat[1][4] = tmat[1][4] - tmp * tmat[0][4];
         v[i][j][k][1] = v[i][j][k][1] - v[i][j][k][0] * tmp;
         tmp = tmp1 * tmat[2][0];
         tmat[2][1] = tmat[2][1] - tmp * tmat[0][1];
         tmat[2][2] = tmat[2][2] - tmp * tmat[0][2];
         tmat[2][3] = tmat[2][3] - tmp * tmat[0][3];
         tmat[2][4] = tmat[2][4] - tmp * tmat[0][4];
         v[i][j][k][2] = v[i][j][k][2] - v[i][j][k][0] * tmp;
         tmp = tmp1 * tmat[3][0];
         tmat[3][1] = tmat[3][1] - tmp * tmat[0][1];
         tmat[3][2] = tmat[3][2] - tmp * tmat[0][2];
         tmat[3][3] = tmat[3][3] - tmp * tmat[0][3];
         tmat[3][4] = tmat[3][4] - tmp * tmat[0][4];
         v[i][j][k][3] = v[i][j][k][3] - v[i][j][k][0] * tmp;
         tmp = tmp1 * tmat[4][0];
         tmat[4][1] = tmat[4][1] - tmp * tmat[0][1];
         tmat[4][2] = tmat[4][2] - tmp * tmat[0][2];
         tmat[4][3] = tmat[4][3] - tmp * tmat[0][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[0][4];
         v[i][j][k][4] = v[i][j][k][4] - v[i][j][k][0] * tmp;
         tmp1 = 1.0 / tmat[1][1];
         tmp = tmp1 * tmat[2][1];
         tmat[2][2] = tmat[2][2] - tmp * tmat[1][2];
         tmat[2][3] = tmat[2][3] - tmp * tmat[1][3];
         tmat[2][4] = tmat[2][4] - tmp * tmat[1][4];
         v[i][j][k][2] = v[i][j][k][2] - v[i][j][k][1] * tmp;
         tmp = tmp1 * tmat[3][1];
         tmat[3][2] = tmat[3][2] - tmp * tmat[1][2];
         tmat[3][3] = tmat[3][3] - tmp * tmat[1][3];
         tmat[3][4] = tmat[3][4] - tmp * tmat[1][4];
         v[i][j][k][3] = v[i][j][k][3] - v[i][j][k][1] * tmp;
         tmp = tmp1 * tmat[4][1];
         tmat[4][2] = tmat[4][2] - tmp * tmat[1][2];
         tmat[4][3] = tmat[4][3] - tmp * tmat[1][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[1][4];
         v[i][j][k][4] = v[i][j][k][4] - v[i][j][k][1] * tmp;
         tmp1 = 1.0 / tmat[2][2];
         tmp = tmp1 * tmat[3][2];
         tmat[3][3] = tmat[3][3] - tmp * tmat[2][3];
         tmat[3][4] = tmat[3][4] - tmp * tmat[2][4];
         v[i][j][k][3] = v[i][j][k][3] - v[i][j][k][2] * tmp;
         tmp = tmp1 * tmat[4][2];
         tmat[4][3] = tmat[4][3] - tmp * tmat[2][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[2][4];
         v[i][j][k][4] = v[i][j][k][4] - v[i][j][k][2] * tmp;
         tmp1 = 1.0 / tmat[3][3];
         tmp = tmp1 * tmat[4][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[3][4];
         v[i][j][k][4] = v[i][j][k][4] - v[i][j][k][3] * tmp;
         v[i][j][k][4] = v[i][j][k][4] / tmat[4][4];
         v[i][j][k][3] = v[i][j][k][3] - tmat[3][4] * v[i][j][k][4];
         v[i][j][k][3] = v[i][j][k][3] / tmat[3][3];
         v[i][j][k][2] = v[i][j][k][2] - tmat[2][3] * v[i][j][k][3] - tmat[2][4] * v[i][j][k][4];
         v[i][j][k][2] = v[i][j][k][2] / tmat[2][2];
         v[i][j][k][1] = v[i][j][k][1] - tmat[1][2] * v[i][j][k][2] - tmat[1][3] * v[i][j][k][3] - tmat[1][4] * v[i][j][k][4];
         v[i][j][k][1] = v[i][j][k][1] / tmat[1][1];
         v[i][j][k][0] = v[i][j][k][0] - tmat[0][1] * v[i][j][k][1] - tmat[0][2] * v[i][j][k][2] - tmat[0][3] * v[i][j][k][3] - tmat[0][4] * v[i][j][k][4];
         v[i][j][k][0] = v[i][j][k][0] / tmat[0][0];
      }
   }
}

static void buts(int nx, int ny, int nz, int k, double omega, double v[64][65][65][5], double tv[64][64][5], double d[64][64][5][5], double udx[64][64][5][5], double udy[64][64][5][5], double udz[64][64][5][5], int ist, int iend, int jst, int jend, int nx0, int ny0) {
   int i, j, m;
   double tmp, tmp1;
   double tmat[5][5];
   #pragma omp parallel for default(shared) private(i, j, m) firstprivate(iend, ist, jend, jst, k, omega, udz, v)
   for(i = iend; i >= ist; i--) {
      #pragma omp parallel for default(shared) private(j, m) firstprivate(jend, jst, k, i, omega, udz, v)
      for(j = jend; j >= jst; j--) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            tv[i][j][m] = omega * (udz[i][j][m][0] * v[i][j][k + 1][0] + udz[i][j][m][1] * v[i][j][k + 1][1] + udz[i][j][m][2] * v[i][j][k + 1][2] + udz[i][j][m][3] * v[i][j][k + 1][3] + udz[i][j][m][4] * v[i][j][k + 1][4]);
         }
      }
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess v	 use : RW
   ****************************************/
   for(i = iend; i >= ist; i--) {
      /*************** Clava msgError **************
      		unsolved dependency for arrayAccess v	 use : RW
      ****************************************/
      for(j = jend; j >= jst; j--) {
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            tv[i][j][m] = tv[i][j][m] + omega * (udy[i][j][m][0] * v[i][j + 1][k][0] + udx[i][j][m][0] * v[i + 1][j][k][0] + udy[i][j][m][1] * v[i][j + 1][k][1] + udx[i][j][m][1] * v[i + 1][j][k][1] + udy[i][j][m][2] * v[i][j + 1][k][2] + udx[i][j][m][2] * v[i + 1][j][k][2] + udy[i][j][m][3] * v[i][j + 1][k][3] + udx[i][j][m][3] * v[i + 1][j][k][3] + udy[i][j][m][4] * v[i][j + 1][k][4] + udx[i][j][m][4] * v[i + 1][j][k][4]);
         }
         /*************** Clava msgError **************
         		 Loop Iteration number is too low
         ****************************************/
         for(m = 0; m < 5; m++) {
            tmat[m][0] = d[i][j][m][0];
            tmat[m][1] = d[i][j][m][1];
            tmat[m][2] = d[i][j][m][2];
            tmat[m][3] = d[i][j][m][3];
            tmat[m][4] = d[i][j][m][4];
         }
         tmp1 = 1.0 / tmat[0][0];
         tmp = tmp1 * tmat[1][0];
         tmat[1][1] = tmat[1][1] - tmp * tmat[0][1];
         tmat[1][2] = tmat[1][2] - tmp * tmat[0][2];
         tmat[1][3] = tmat[1][3] - tmp * tmat[0][3];
         tmat[1][4] = tmat[1][4] - tmp * tmat[0][4];
         tv[i][j][1] = tv[i][j][1] - tv[i][j][0] * tmp;
         tmp = tmp1 * tmat[2][0];
         tmat[2][1] = tmat[2][1] - tmp * tmat[0][1];
         tmat[2][2] = tmat[2][2] - tmp * tmat[0][2];
         tmat[2][3] = tmat[2][3] - tmp * tmat[0][3];
         tmat[2][4] = tmat[2][4] - tmp * tmat[0][4];
         tv[i][j][2] = tv[i][j][2] - tv[i][j][0] * tmp;
         tmp = tmp1 * tmat[3][0];
         tmat[3][1] = tmat[3][1] - tmp * tmat[0][1];
         tmat[3][2] = tmat[3][2] - tmp * tmat[0][2];
         tmat[3][3] = tmat[3][3] - tmp * tmat[0][3];
         tmat[3][4] = tmat[3][4] - tmp * tmat[0][4];
         tv[i][j][3] = tv[i][j][3] - tv[i][j][0] * tmp;
         tmp = tmp1 * tmat[4][0];
         tmat[4][1] = tmat[4][1] - tmp * tmat[0][1];
         tmat[4][2] = tmat[4][2] - tmp * tmat[0][2];
         tmat[4][3] = tmat[4][3] - tmp * tmat[0][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[0][4];
         tv[i][j][4] = tv[i][j][4] - tv[i][j][0] * tmp;
         tmp1 = 1.0 / tmat[1][1];
         tmp = tmp1 * tmat[2][1];
         tmat[2][2] = tmat[2][2] - tmp * tmat[1][2];
         tmat[2][3] = tmat[2][3] - tmp * tmat[1][3];
         tmat[2][4] = tmat[2][4] - tmp * tmat[1][4];
         tv[i][j][2] = tv[i][j][2] - tv[i][j][1] * tmp;
         tmp = tmp1 * tmat[3][1];
         tmat[3][2] = tmat[3][2] - tmp * tmat[1][2];
         tmat[3][3] = tmat[3][3] - tmp * tmat[1][3];
         tmat[3][4] = tmat[3][4] - tmp * tmat[1][4];
         tv[i][j][3] = tv[i][j][3] - tv[i][j][1] * tmp;
         tmp = tmp1 * tmat[4][1];
         tmat[4][2] = tmat[4][2] - tmp * tmat[1][2];
         tmat[4][3] = tmat[4][3] - tmp * tmat[1][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[1][4];
         tv[i][j][4] = tv[i][j][4] - tv[i][j][1] * tmp;
         tmp1 = 1.0 / tmat[2][2];
         tmp = tmp1 * tmat[3][2];
         tmat[3][3] = tmat[3][3] - tmp * tmat[2][3];
         tmat[3][4] = tmat[3][4] - tmp * tmat[2][4];
         tv[i][j][3] = tv[i][j][3] - tv[i][j][2] * tmp;
         tmp = tmp1 * tmat[4][2];
         tmat[4][3] = tmat[4][3] - tmp * tmat[2][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[2][4];
         tv[i][j][4] = tv[i][j][4] - tv[i][j][2] * tmp;
         tmp1 = 1.0 / tmat[3][3];
         tmp = tmp1 * tmat[4][3];
         tmat[4][4] = tmat[4][4] - tmp * tmat[3][4];
         tv[i][j][4] = tv[i][j][4] - tv[i][j][3] * tmp;
         tv[i][j][4] = tv[i][j][4] / tmat[4][4];
         tv[i][j][3] = tv[i][j][3] - tmat[3][4] * tv[i][j][4];
         tv[i][j][3] = tv[i][j][3] / tmat[3][3];
         tv[i][j][2] = tv[i][j][2] - tmat[2][3] * tv[i][j][3] - tmat[2][4] * tv[i][j][4];
         tv[i][j][2] = tv[i][j][2] / tmat[2][2];
         tv[i][j][1] = tv[i][j][1] - tmat[1][2] * tv[i][j][2] - tmat[1][3] * tv[i][j][3] - tmat[1][4] * tv[i][j][4];
         tv[i][j][1] = tv[i][j][1] / tmat[1][1];
         tv[i][j][0] = tv[i][j][0] - tmat[0][1] * tv[i][j][1] - tmat[0][2] * tv[i][j][2] - tmat[0][3] * tv[i][j][3] - tmat[0][4] * tv[i][j][4];
         tv[i][j][0] = tv[i][j][0] / tmat[0][0];
         v[i][j][k][0] = v[i][j][k][0] - tv[i][j][0];
         v[i][j][k][1] = v[i][j][k][1] - tv[i][j][1];
         v[i][j][k][2] = v[i][j][k][2] - tv[i][j][2];
         v[i][j][k][3] = v[i][j][k][3] - tv[i][j][3];
         v[i][j][k][4] = v[i][j][k][4] - tv[i][j][4];
      }
   }
}

static void domain() {
   nx = nx0;
   ny = ny0;
   nz = nz0;
   if(nx < 4 || ny < 4 || nz < 4) {
      printf("     SUBDOMAIN SIZE IS TOO SMALL - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n     TO 4 THEY ARE CURRENTLY%3d%3d%3d\n", nx, ny, nz);
      exit(1);
   }
   if(nx > 64 || ny > 64 || nz > 64) {
      printf("     SUBDOMAIN SIZE IS TOO LARGE - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n     CURRENTLY%4d%4d%4d\n", nx, ny, nz);
      exit(1);
   }
   ist = 1;
   iend = nx - 2;
   jst = 1;
   jend = ny - 2;
}

static void erhs() {
   {
      int i, j, k, m;
      int iglob, jglob;
      int L1, L2;
      int ist1, iend1;
      int jst1, jend1;
      double dsspm;
      double xi, eta, zeta;
      double q;
      double u21, u31, u41;
      double tmp;
      double u21i, u31i, u41i, u51i;
      double u21j, u31j, u41j, u51j;
      double u21k, u31k, u41k, u51k;
      double u21im1, u31im1, u41im1, u51im1;
      double u21jm1, u31jm1, u41jm1, u51jm1;
      double u21km1, u31km1, u41km1, u51km1;
      dsspm = dssp;
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(nx, ny, nz)
      for(i = 0; i < nx; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(ny, nz)
         for(j = 0; j < ny; j++) {
            #pragma omp parallel for default(shared) private(k, m) firstprivate(nz)
            for(k = 0; k < nz; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  frct[i][j][k][m] = 0.0;
               }
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(nx, nx0, ny, ny0, nz, ce)
      for(i = 0; i < nx; i++) {
         iglob = i;
         xi = ((double) (iglob)) / (nx0 - 1);
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(ny, ny0, nz, ce)
         for(j = 0; j < ny; j++) {
            jglob = j;
            eta = ((double) (jglob)) / (ny0 - 1);
            #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, ce)
            for(k = 0; k < nz; k++) {
               zeta = ((double) (k)) / (nz - 1);
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = ce[m][0] + ce[m][1] * xi + ce[m][2] * eta + ce[m][3] * zeta + ce[m][4] * xi * xi + ce[m][5] * eta * eta + ce[m][6] * zeta * zeta + ce[m][7] * xi * xi * xi + ce[m][8] * eta * eta * eta + ce[m][9] * zeta * zeta * zeta + ce[m][10] * xi * xi * xi * xi + ce[m][11] * eta * eta * eta * eta + ce[m][12] * zeta * zeta * zeta * zeta;
               }
            }
         }
      }
      L1 = 0;
      L2 = nx - 1;
      #pragma omp parallel for default(shared) private(i, j, k) firstprivate(jst, jend, nz, rsd)
      for(i = L1; i <= L2; i++) {
         #pragma omp parallel for default(shared) private(j, k) firstprivate(jst, jend, nz, rsd)
         for(j = jst; j <= jend; j++) {
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, rsd)
            for(k = 1; k < nz - 1; k++) {
               flux[i][j][k][0] = rsd[i][j][k][1];
               u21 = rsd[i][j][k][1] / rsd[i][j][k][0];
               q = 0.50 * (rsd[i][j][k][1] * rsd[i][j][k][1] + rsd[i][j][k][2] * rsd[i][j][k][2] + rsd[i][j][k][3] * rsd[i][j][k][3]) / rsd[i][j][k][0];
               flux[i][j][k][1] = rsd[i][j][k][1] * u21 + 0.40e+00 * (rsd[i][j][k][4] - q);
               flux[i][j][k][2] = rsd[i][j][k][2] * u21;
               flux[i][j][k][3] = rsd[i][j][k][3] * u21;
               flux[i][j][k][4] = (1.40e+00 * rsd[i][j][k][4] - 0.40e+00 * q) * u21;
            }
         }
      }
      #pragma omp parallel for default(shared) private(j, k, i, m) firstprivate(jst, jend, nz, ist, iend, tx2, tx3, dx1, tx1, dx2, dx3, dx4, dx5, nx, rsd)
      for(j = jst; j <= jend; j++) {
         #pragma omp parallel for default(shared) private(k, i, m) firstprivate(nz, ist, iend, tx2, tx3, dx1, tx1, dx2, dx3, dx4, dx5, nx, rsd)
         for(k = 1; k <= nz - 2; k++) {
            #pragma omp parallel for default(shared) private(i, m) firstprivate(ist, iend, tx2, flux)
            for(i = ist; i <= iend; i++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  frct[i][j][k][m] = frct[i][j][k][m] - tx2 * (flux[i + 1][j][k][m] - flux[i - 1][j][k][m]);
               }
            }
            #pragma omp parallel for default(shared) private(i) firstprivate(ist, tx3, rsd)
            for(i = ist; i <= L2; i++) {
               tmp = 1.0 / rsd[i][j][k][0];
               u21i = tmp * rsd[i][j][k][1];
               u31i = tmp * rsd[i][j][k][2];
               u41i = tmp * rsd[i][j][k][3];
               u51i = tmp * rsd[i][j][k][4];
               tmp = 1.0 / rsd[i - 1][j][k][0];
               u21im1 = tmp * rsd[i - 1][j][k][1];
               u31im1 = tmp * rsd[i - 1][j][k][2];
               u41im1 = tmp * rsd[i - 1][j][k][3];
               u51im1 = tmp * rsd[i - 1][j][k][4];
               flux[i][j][k][1] = (4.0 / 3.0) * tx3 * (u21i - u21im1);
               flux[i][j][k][2] = tx3 * (u31i - u31im1);
               flux[i][j][k][3] = tx3 * (u41i - u41im1);
               flux[i][j][k][4] = 0.50 * (1.0 - 1.40e+00 * 1.40e+00) * tx3 * ((u21i * u21i + u31i * u31i + u41i * u41i) - (u21im1 * u21im1 + u31im1 * u31im1 + u41im1 * u41im1)) + (1.0 / 6.0) * tx3 * (u21i * u21i - u21im1 * u21im1) + 1.40e+00 * 1.40e+00 * tx3 * (u51i - u51im1);
            }
            #pragma omp parallel for default(shared) private(i) firstprivate(ist, iend, dx1, tx1, tx3, dx2, dx3, dx4, dx5, rsd, flux)
            for(i = ist; i <= iend; i++) {
               frct[i][j][k][0] = frct[i][j][k][0] + dx1 * tx1 * (rsd[i - 1][j][k][0] - 2.0 * rsd[i][j][k][0] + rsd[i + 1][j][k][0]);
               frct[i][j][k][1] = frct[i][j][k][1] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][1] - flux[i][j][k][1]) + dx2 * tx1 * (rsd[i - 1][j][k][1] - 2.0 * rsd[i][j][k][1] + rsd[i + 1][j][k][1]);
               frct[i][j][k][2] = frct[i][j][k][2] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][2] - flux[i][j][k][2]) + dx3 * tx1 * (rsd[i - 1][j][k][2] - 2.0 * rsd[i][j][k][2] + rsd[i + 1][j][k][2]);
               frct[i][j][k][3] = frct[i][j][k][3] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][3] - flux[i][j][k][3]) + dx4 * tx1 * (rsd[i - 1][j][k][3] - 2.0 * rsd[i][j][k][3] + rsd[i + 1][j][k][3]);
               frct[i][j][k][4] = frct[i][j][k][4] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][4] - flux[i][j][k][4]) + dx5 * tx1 * (rsd[i - 1][j][k][4] - 2.0 * rsd[i][j][k][4] + rsd[i + 1][j][k][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               frct[1][j][k][m] = frct[1][j][k][m] - dsspm * (+5.0 * rsd[1][j][k][m] - 4.0 * rsd[2][j][k][m] + rsd[3][j][k][m]);
               frct[2][j][k][m] = frct[2][j][k][m] - dsspm * (-4.0 * rsd[1][j][k][m] + 6.0 * rsd[2][j][k][m] - 4.0 * rsd[3][j][k][m] + rsd[4][j][k][m]);
            }
            ist1 = 3;
            iend1 = nx - 4;
            #pragma omp parallel for default(shared) private(i, m) firstprivate(rsd)
            for(i = ist1; i <= iend1; i++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  frct[i][j][k][m] = frct[i][j][k][m] - dsspm * (rsd[i - 2][j][k][m] - 4.0 * rsd[i - 1][j][k][m] + 6.0 * rsd[i][j][k][m] - 4.0 * rsd[i + 1][j][k][m] + rsd[i + 2][j][k][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               frct[nx - 3][j][k][m] = frct[nx - 3][j][k][m] - dsspm * (rsd[nx - 5][j][k][m] - 4.0 * rsd[nx - 4][j][k][m] + 6.0 * rsd[nx - 3][j][k][m] - 4.0 * rsd[nx - 2][j][k][m]);
               frct[nx - 2][j][k][m] = frct[nx - 2][j][k][m] - dsspm * (rsd[nx - 4][j][k][m] - 4.0 * rsd[nx - 3][j][k][m] + 5.0 * rsd[nx - 2][j][k][m]);
            }
         }
      }
      L1 = 0;
      L2 = ny - 1;
      #pragma omp parallel for default(shared) private(i, j, k) firstprivate(ist, iend, nz, rsd)
      for(i = ist; i <= iend; i++) {
         #pragma omp parallel for default(shared) private(j, k) firstprivate(nz, rsd)
         for(j = L1; j <= L2; j++) {
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, rsd)
            for(k = 1; k <= nz - 2; k++) {
               flux[i][j][k][0] = rsd[i][j][k][2];
               u31 = rsd[i][j][k][2] / rsd[i][j][k][0];
               q = 0.50 * (rsd[i][j][k][1] * rsd[i][j][k][1] + rsd[i][j][k][2] * rsd[i][j][k][2] + rsd[i][j][k][3] * rsd[i][j][k][3]) / rsd[i][j][k][0];
               flux[i][j][k][1] = rsd[i][j][k][1] * u31;
               flux[i][j][k][2] = rsd[i][j][k][2] * u31 + 0.40e+00 * (rsd[i][j][k][4] - q);
               flux[i][j][k][3] = rsd[i][j][k][3] * u31;
               flux[i][j][k][4] = (1.40e+00 * rsd[i][j][k][4] - 0.40e+00 * q) * u31;
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, k, j, m) firstprivate(ist, iend, nz, jst, jend, ty2, ty3, dy1, ty1, dy2, dy3, dy4, dy5, ny, rsd)
      for(i = ist; i <= iend; i++) {
         #pragma omp parallel for default(shared) private(k, j, m) firstprivate(nz, jst, jend, ty2, ty3, dy1, ty1, dy2, dy3, dy4, dy5, ny, rsd)
         for(k = 1; k <= nz - 2; k++) {
            #pragma omp parallel for default(shared) private(j, m) firstprivate(jst, jend, ty2, flux)
            for(j = jst; j <= jend; j++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  frct[i][j][k][m] = frct[i][j][k][m] - ty2 * (flux[i][j + 1][k][m] - flux[i][j - 1][k][m]);
               }
            }
            #pragma omp parallel for default(shared) private(j) firstprivate(jst, ty3, rsd)
            for(j = jst; j <= L2; j++) {
               tmp = 1.0 / rsd[i][j][k][0];
               u21j = tmp * rsd[i][j][k][1];
               u31j = tmp * rsd[i][j][k][2];
               u41j = tmp * rsd[i][j][k][3];
               u51j = tmp * rsd[i][j][k][4];
               tmp = 1.0 / rsd[i][j - 1][k][0];
               u21jm1 = tmp * rsd[i][j - 1][k][1];
               u31jm1 = tmp * rsd[i][j - 1][k][2];
               u41jm1 = tmp * rsd[i][j - 1][k][3];
               u51jm1 = tmp * rsd[i][j - 1][k][4];
               flux[i][j][k][1] = ty3 * (u21j - u21jm1);
               flux[i][j][k][2] = (4.0 / 3.0) * ty3 * (u31j - u31jm1);
               flux[i][j][k][3] = ty3 * (u41j - u41jm1);
               flux[i][j][k][4] = 0.50 * (1.0 - 1.40e+00 * 1.40e+00) * ty3 * ((u21j * u21j + u31j * u31j + u41j * u41j) - (u21jm1 * u21jm1 + u31jm1 * u31jm1 + u41jm1 * u41jm1)) + (1.0 / 6.0) * ty3 * (u31j * u31j - u31jm1 * u31jm1) + 1.40e+00 * 1.40e+00 * ty3 * (u51j - u51jm1);
            }
            #pragma omp parallel for default(shared) private(j) firstprivate(jst, jend, dy1, ty1, ty3, dy2, dy3, dy4, dy5, rsd, flux)
            for(j = jst; j <= jend; j++) {
               frct[i][j][k][0] = frct[i][j][k][0] + dy1 * ty1 * (rsd[i][j - 1][k][0] - 2.0 * rsd[i][j][k][0] + rsd[i][j + 1][k][0]);
               frct[i][j][k][1] = frct[i][j][k][1] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][1] - flux[i][j][k][1]) + dy2 * ty1 * (rsd[i][j - 1][k][1] - 2.0 * rsd[i][j][k][1] + rsd[i][j + 1][k][1]);
               frct[i][j][k][2] = frct[i][j][k][2] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][2] - flux[i][j][k][2]) + dy3 * ty1 * (rsd[i][j - 1][k][2] - 2.0 * rsd[i][j][k][2] + rsd[i][j + 1][k][2]);
               frct[i][j][k][3] = frct[i][j][k][3] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][3] - flux[i][j][k][3]) + dy4 * ty1 * (rsd[i][j - 1][k][3] - 2.0 * rsd[i][j][k][3] + rsd[i][j + 1][k][3]);
               frct[i][j][k][4] = frct[i][j][k][4] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][4] - flux[i][j][k][4]) + dy5 * ty1 * (rsd[i][j - 1][k][4] - 2.0 * rsd[i][j][k][4] + rsd[i][j + 1][k][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               frct[i][1][k][m] = frct[i][1][k][m] - dsspm * (+5.0 * rsd[i][1][k][m] - 4.0 * rsd[i][2][k][m] + rsd[i][3][k][m]);
               frct[i][2][k][m] = frct[i][2][k][m] - dsspm * (-4.0 * rsd[i][1][k][m] + 6.0 * rsd[i][2][k][m] - 4.0 * rsd[i][3][k][m] + rsd[i][4][k][m]);
            }
            jst1 = 3;
            jend1 = ny - 4;
            #pragma omp parallel for default(shared) private(j, m) firstprivate(rsd)
            for(j = jst1; j <= jend1; j++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  frct[i][j][k][m] = frct[i][j][k][m] - dsspm * (rsd[i][j - 2][k][m] - 4.0 * rsd[i][j - 1][k][m] + 6.0 * rsd[i][j][k][m] - 4.0 * rsd[i][j + 1][k][m] + rsd[i][j + 2][k][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               frct[i][ny - 3][k][m] = frct[i][ny - 3][k][m] - dsspm * (rsd[i][ny - 5][k][m] - 4.0 * rsd[i][ny - 4][k][m] + 6.0 * rsd[i][ny - 3][k][m] - 4.0 * rsd[i][ny - 2][k][m]);
               frct[i][ny - 2][k][m] = frct[i][ny - 2][k][m] - dsspm * (rsd[i][ny - 4][k][m] - 4.0 * rsd[i][ny - 3][k][m] + 5.0 * rsd[i][ny - 2][k][m]);
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(ist, iend, jst, jend, nz, tz2, tz3, dz1, tz1, dz2, dz3, dz4, dz5, rsd)
      for(i = ist; i <= iend; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(jst, jend, nz, tz2, tz3, dz1, tz1, dz2, dz3, dz4, dz5, rsd)
         for(j = jst; j <= jend; j++) {
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, rsd)
            for(k = 0; k <= nz - 1; k++) {
               flux[i][j][k][0] = rsd[i][j][k][3];
               u41 = rsd[i][j][k][3] / rsd[i][j][k][0];
               q = 0.50 * (rsd[i][j][k][1] * rsd[i][j][k][1] + rsd[i][j][k][2] * rsd[i][j][k][2] + rsd[i][j][k][3] * rsd[i][j][k][3]) / rsd[i][j][k][0];
               flux[i][j][k][1] = rsd[i][j][k][1] * u41;
               flux[i][j][k][2] = rsd[i][j][k][2] * u41;
               flux[i][j][k][3] = rsd[i][j][k][3] * u41 + 0.40e+00 * (rsd[i][j][k][4] - q);
               flux[i][j][k][4] = (1.40e+00 * rsd[i][j][k][4] - 0.40e+00 * q) * u41;
            }
            #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, tz2, flux)
            for(k = 1; k <= nz - 2; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  frct[i][j][k][m] = frct[i][j][k][m] - tz2 * (flux[i][j][k + 1][m] - flux[i][j][k - 1][m]);
               }
            }
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, tz3, rsd)
            for(k = 1; k <= nz - 1; k++) {
               tmp = 1.0 / rsd[i][j][k][0];
               u21k = tmp * rsd[i][j][k][1];
               u31k = tmp * rsd[i][j][k][2];
               u41k = tmp * rsd[i][j][k][3];
               u51k = tmp * rsd[i][j][k][4];
               tmp = 1.0 / rsd[i][j][k - 1][0];
               u21km1 = tmp * rsd[i][j][k - 1][1];
               u31km1 = tmp * rsd[i][j][k - 1][2];
               u41km1 = tmp * rsd[i][j][k - 1][3];
               u51km1 = tmp * rsd[i][j][k - 1][4];
               flux[i][j][k][1] = tz3 * (u21k - u21km1);
               flux[i][j][k][2] = tz3 * (u31k - u31km1);
               flux[i][j][k][3] = (4.0 / 3.0) * tz3 * (u41k - u41km1);
               flux[i][j][k][4] = 0.50 * (1.0 - 1.40e+00 * 1.40e+00) * tz3 * ((u21k * u21k + u31k * u31k + u41k * u41k) - (u21km1 * u21km1 + u31km1 * u31km1 + u41km1 * u41km1)) + (1.0 / 6.0) * tz3 * (u41k * u41k - u41km1 * u41km1) + 1.40e+00 * 1.40e+00 * tz3 * (u51k - u51km1);
            }
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, dz1, tz1, tz3, dz2, dz3, dz4, dz5, rsd, flux)
            for(k = 1; k <= nz - 2; k++) {
               frct[i][j][k][0] = frct[i][j][k][0] + dz1 * tz1 * (rsd[i][j][k + 1][0] - 2.0 * rsd[i][j][k][0] + rsd[i][j][k - 1][0]);
               frct[i][j][k][1] = frct[i][j][k][1] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][1] - flux[i][j][k][1]) + dz2 * tz1 * (rsd[i][j][k + 1][1] - 2.0 * rsd[i][j][k][1] + rsd[i][j][k - 1][1]);
               frct[i][j][k][2] = frct[i][j][k][2] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][2] - flux[i][j][k][2]) + dz3 * tz1 * (rsd[i][j][k + 1][2] - 2.0 * rsd[i][j][k][2] + rsd[i][j][k - 1][2]);
               frct[i][j][k][3] = frct[i][j][k][3] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][3] - flux[i][j][k][3]) + dz4 * tz1 * (rsd[i][j][k + 1][3] - 2.0 * rsd[i][j][k][3] + rsd[i][j][k - 1][3]);
               frct[i][j][k][4] = frct[i][j][k][4] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][4] - flux[i][j][k][4]) + dz5 * tz1 * (rsd[i][j][k + 1][4] - 2.0 * rsd[i][j][k][4] + rsd[i][j][k - 1][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               frct[i][j][1][m] = frct[i][j][1][m] - dsspm * (+5.0 * rsd[i][j][1][m] - 4.0 * rsd[i][j][2][m] + rsd[i][j][3][m]);
               frct[i][j][2][m] = frct[i][j][2][m] - dsspm * (-4.0 * rsd[i][j][1][m] + 6.0 * rsd[i][j][2][m] - 4.0 * rsd[i][j][3][m] + rsd[i][j][4][m]);
            }
            #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, rsd)
            for(k = 3; k <= nz - 4; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  frct[i][j][k][m] = frct[i][j][k][m] - dsspm * (rsd[i][j][k - 2][m] - 4.0 * rsd[i][j][k - 1][m] + 6.0 * rsd[i][j][k][m] - 4.0 * rsd[i][j][k + 1][m] + rsd[i][j][k + 2][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               frct[i][j][nz - 3][m] = frct[i][j][nz - 3][m] - dsspm * (rsd[i][j][nz - 5][m] - 4.0 * rsd[i][j][nz - 4][m] + 6.0 * rsd[i][j][nz - 3][m] - 4.0 * rsd[i][j][nz - 2][m]);
               frct[i][j][nz - 2][m] = frct[i][j][nz - 2][m] - dsspm * (rsd[i][j][nz - 4][m] - 4.0 * rsd[i][j][nz - 3][m] + 5.0 * rsd[i][j][nz - 2][m]);
            }
         }
      }
   }
}

static void error() {
   int i, j, k, m;
   int iglob, jglob;
   double tmp;
   double u000ijk[5];
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      errnm[m] = 0.0;
   }
   #pragma omp parallel for default(shared) private(i, j, k, m, iglob, jglob, tmp) firstprivate(ist, iend, jst, jend, nz, nx0, ny0, ce, u, u000ijk) reduction(+ : errnm[:5])
   for(i = ist; i <= iend; i++) {
      iglob = i;
      #pragma omp parallel for default(shared) private(j, k, m, jglob, tmp) firstprivate(jst, jend, nz, iglob, nx0, ny0, i, ce, u, u000ijk) reduction(+ : errnm[:5])
      for(j = jst; j <= jend; j++) {
         jglob = j;
         #pragma omp parallel for default(shared) private(k, m, tmp) firstprivate(nz, jglob, iglob, nx0, ny0, i, j, ce, u, u000ijk) reduction(+ : errnm[:5])
         for(k = 1; k <= nz - 2; k++) {
            exact(iglob, jglob, k, u000ijk);
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               tmp = (u000ijk[m] - u[i][j][k][m]);
               errnm[m] = errnm[m] + tmp * tmp;
            }
         }
      }
   }
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      errnm[m] = sqrt(errnm[m] / ((nx0 - 2) * (ny0 - 2) * (nz0 - 2)));
   }
}

static void exact(int i, int j, int k, double u000ijk[5]) {
   int m;
   double xi, eta, zeta;
   xi = ((double) i) / (nx0 - 1);
   eta = ((double) j) / (ny0 - 1);
   zeta = ((double) k) / (nz - 1);
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      u000ijk[m] = ce[m][0] + ce[m][1] * xi + ce[m][2] * eta + ce[m][3] * zeta + ce[m][4] * xi * xi + ce[m][5] * eta * eta + ce[m][6] * zeta * zeta + ce[m][7] * xi * xi * xi + ce[m][8] * eta * eta * eta + ce[m][9] * zeta * zeta * zeta + ce[m][10] * xi * xi * xi * xi + ce[m][11] * eta * eta * eta * eta + ce[m][12] * zeta * zeta * zeta * zeta;
   }
}

static void jacld(int k) {
   int i, j;
   double r43;
   double c1345;
   double c34;
   double tmp1, tmp2, tmp3;
   r43 = (4.0 / 3.0);
   c1345 = 1.40e+00 * 1.00e-01 * 1.00e+00 * 1.40e+00;
   c34 = 1.00e-01 * 1.00e+00;
   #pragma omp parallel for default(shared) private(i, j, tmp1, tmp2, tmp3) firstprivate(ist, iend, jst, jend, k, tx1, dx1, ty1, dy1, tz1, dz1, dt, r43, c34, dx2, dy2, dz2, dx3, dy3, dz3, dx4, dy4, dz4, c1345, dx5, dy5, dz5, tz2, ty2, tx2, u)
   for(i = ist; i <= iend; i++) {
      #pragma omp parallel for default(shared) private(j, tmp1, tmp2, tmp3) firstprivate(jst, jend, i, k, tx1, dx1, ty1, dy1, tz1, dz1, dt, r43, c34, dx2, dy2, dz2, dx3, dy3, dz3, dx4, dy4, dz4, c1345, dx5, dy5, dz5, tz2, ty2, tx2, u)
      for(j = jst; j <= jend; j++) {
         tmp1 = 1.0 / u[i][j][k][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         d[i][j][0][0] = 1.0 + dt * 2.0 * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
         d[i][j][0][1] = 0.0;
         d[i][j][0][2] = 0.0;
         d[i][j][0][3] = 0.0;
         d[i][j][0][4] = 0.0;
         d[i][j][1][0] = dt * 2.0 * (tx1 * (-r43 * c34 * tmp2 * u[i][j][k][1]) + ty1 * (-c34 * tmp2 * u[i][j][k][1]) + tz1 * (-c34 * tmp2 * u[i][j][k][1]));
         d[i][j][1][1] = 1.0 + dt * 2.0 * (tx1 * r43 * c34 * tmp1 + ty1 * c34 * tmp1 + tz1 * c34 * tmp1) + dt * 2.0 * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
         d[i][j][1][2] = 0.0;
         d[i][j][1][3] = 0.0;
         d[i][j][1][4] = 0.0;
         d[i][j][2][0] = dt * 2.0 * (tx1 * (-c34 * tmp2 * u[i][j][k][2]) + ty1 * (-r43 * c34 * tmp2 * u[i][j][k][2]) + tz1 * (-c34 * tmp2 * u[i][j][k][2]));
         d[i][j][2][1] = 0.0;
         d[i][j][2][2] = 1.0 + dt * 2.0 * (tx1 * c34 * tmp1 + ty1 * r43 * c34 * tmp1 + tz1 * c34 * tmp1) + dt * 2.0 * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
         d[i][j][2][3] = 0.0;
         d[i][j][2][4] = 0.0;
         d[i][j][3][0] = dt * 2.0 * (tx1 * (-c34 * tmp2 * u[i][j][k][3]) + ty1 * (-c34 * tmp2 * u[i][j][k][3]) + tz1 * (-r43 * c34 * tmp2 * u[i][j][k][3]));
         d[i][j][3][1] = 0.0;
         d[i][j][3][2] = 0.0;
         d[i][j][3][3] = 1.0 + dt * 2.0 * (tx1 * c34 * tmp1 + ty1 * c34 * tmp1 + tz1 * r43 * c34 * tmp1) + dt * 2.0 * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
         d[i][j][3][4] = 0.0;
         d[i][j][4][0] = dt * 2.0 * (tx1 * (-(r43 * c34 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - (c1345) * tmp2 * u[i][j][k][4]) + ty1 * (-(c34 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (r43 * c34 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - (c1345) * tmp2 * u[i][j][k][4]) + tz1 * (-(c34 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (r43 * c34 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - (c1345) * tmp2 * u[i][j][k][4]));
         d[i][j][4][1] = dt * 2.0 * (tx1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k][1] + ty1 * (c34 - c1345) * tmp2 * u[i][j][k][1] + tz1 * (c34 - c1345) * tmp2 * u[i][j][k][1]);
         d[i][j][4][2] = dt * 2.0 * (tx1 * (c34 - c1345) * tmp2 * u[i][j][k][2] + ty1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k][2] + tz1 * (c34 - c1345) * tmp2 * u[i][j][k][2]);
         d[i][j][4][3] = dt * 2.0 * (tx1 * (c34 - c1345) * tmp2 * u[i][j][k][3] + ty1 * (c34 - c1345) * tmp2 * u[i][j][k][3] + tz1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k][3]);
         d[i][j][4][4] = 1.0 + dt * 2.0 * (tx1 * c1345 * tmp1 + ty1 * c1345 * tmp1 + tz1 * c1345 * tmp1) + dt * 2.0 * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
         tmp1 = 1.0 / u[i][j][k - 1][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         a[i][j][0][0] = -dt * tz1 * dz1;
         a[i][j][0][1] = 0.0;
         a[i][j][0][2] = 0.0;
         a[i][j][0][3] = -dt * tz2;
         a[i][j][0][4] = 0.0;
         a[i][j][1][0] = -dt * tz2 * (-(u[i][j][k - 1][1] * u[i][j][k - 1][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[i][j][k - 1][1]);
         a[i][j][1][1] = -dt * tz2 * (u[i][j][k - 1][3] * tmp1) - dt * tz1 * c34 * tmp1 - dt * tz1 * dz2;
         a[i][j][1][2] = 0.0;
         a[i][j][1][3] = -dt * tz2 * (u[i][j][k - 1][1] * tmp1);
         a[i][j][1][4] = 0.0;
         a[i][j][2][0] = -dt * tz2 * (-(u[i][j][k - 1][2] * u[i][j][k - 1][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[i][j][k - 1][2]);
         a[i][j][2][1] = 0.0;
         a[i][j][2][2] = -dt * tz2 * (u[i][j][k - 1][3] * tmp1) - dt * tz1 * (c34 * tmp1) - dt * tz1 * dz3;
         a[i][j][2][3] = -dt * tz2 * (u[i][j][k - 1][2] * tmp1);
         a[i][j][2][4] = 0.0;
         a[i][j][3][0] = -dt * tz2 * (-(u[i][j][k - 1][3] * tmp1) * (u[i][j][k - 1][3] * tmp1) + 0.50 * 0.40e+00 * ((u[i][j][k - 1][1] * u[i][j][k - 1][1] + u[i][j][k - 1][2] * u[i][j][k - 1][2] + u[i][j][k - 1][3] * u[i][j][k - 1][3]) * tmp2)) - dt * tz1 * (-r43 * c34 * tmp2 * u[i][j][k - 1][3]);
         a[i][j][3][1] = -dt * tz2 * (-0.40e+00 * (u[i][j][k - 1][1] * tmp1));
         a[i][j][3][2] = -dt * tz2 * (-0.40e+00 * (u[i][j][k - 1][2] * tmp1));
         a[i][j][3][3] = -dt * tz2 * (2.0 - 0.40e+00) * (u[i][j][k - 1][3] * tmp1) - dt * tz1 * (r43 * c34 * tmp1) - dt * tz1 * dz4;
         a[i][j][3][4] = -dt * tz2 * 0.40e+00;
         a[i][j][4][0] = -dt * tz2 * ((0.40e+00 * (u[i][j][k - 1][1] * u[i][j][k - 1][1] + u[i][j][k - 1][2] * u[i][j][k - 1][2] + u[i][j][k - 1][3] * u[i][j][k - 1][3]) * tmp2 - 1.40e+00 * (u[i][j][k - 1][4] * tmp1)) * (u[i][j][k - 1][3] * tmp1)) - dt * tz1 * (-(c34 - c1345) * tmp3 * (u[i][j][k - 1][1] * u[i][j][k - 1][1]) - (c34 - c1345) * tmp3 * (u[i][j][k - 1][2] * u[i][j][k - 1][2]) - (r43 * c34 - c1345) * tmp3 * (u[i][j][k - 1][3] * u[i][j][k - 1][3]) - c1345 * tmp2 * u[i][j][k - 1][4]);
         a[i][j][4][1] = -dt * tz2 * (-0.40e+00 * (u[i][j][k - 1][1] * u[i][j][k - 1][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[i][j][k - 1][1];
         a[i][j][4][2] = -dt * tz2 * (-0.40e+00 * (u[i][j][k - 1][2] * u[i][j][k - 1][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[i][j][k - 1][2];
         a[i][j][4][3] = -dt * tz2 * (1.40e+00 * (u[i][j][k - 1][4] * tmp1) - 0.50 * 0.40e+00 * ((u[i][j][k - 1][1] * u[i][j][k - 1][1] + u[i][j][k - 1][2] * u[i][j][k - 1][2] + 3.0 * u[i][j][k - 1][3] * u[i][j][k - 1][3]) * tmp2)) - dt * tz1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k - 1][3];
         a[i][j][4][4] = -dt * tz2 * (1.40e+00 * (u[i][j][k - 1][3] * tmp1)) - dt * tz1 * c1345 * tmp1 - dt * tz1 * dz5;
         tmp1 = 1.0 / u[i][j - 1][k][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         b[i][j][0][0] = -dt * ty1 * dy1;
         b[i][j][0][1] = 0.0;
         b[i][j][0][2] = -dt * ty2;
         b[i][j][0][3] = 0.0;
         b[i][j][0][4] = 0.0;
         b[i][j][1][0] = -dt * ty2 * (-(u[i][j - 1][k][1] * u[i][j - 1][k][2]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[i][j - 1][k][1]);
         b[i][j][1][1] = -dt * ty2 * (u[i][j - 1][k][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy2;
         b[i][j][1][2] = -dt * ty2 * (u[i][j - 1][k][1] * tmp1);
         b[i][j][1][3] = 0.0;
         b[i][j][1][4] = 0.0;
         b[i][j][2][0] = -dt * ty2 * (-(u[i][j - 1][k][2] * tmp1) * (u[i][j - 1][k][2] * tmp1) + 0.50 * 0.40e+00 * ((u[i][j - 1][k][1] * u[i][j - 1][k][1] + u[i][j - 1][k][2] * u[i][j - 1][k][2] + u[i][j - 1][k][3] * u[i][j - 1][k][3]) * tmp2)) - dt * ty1 * (-r43 * c34 * tmp2 * u[i][j - 1][k][2]);
         b[i][j][2][1] = -dt * ty2 * (-0.40e+00 * (u[i][j - 1][k][1] * tmp1));
         b[i][j][2][2] = -dt * ty2 * ((2.0 - 0.40e+00) * (u[i][j - 1][k][2] * tmp1)) - dt * ty1 * (r43 * c34 * tmp1) - dt * ty1 * dy3;
         b[i][j][2][3] = -dt * ty2 * (-0.40e+00 * (u[i][j - 1][k][3] * tmp1));
         b[i][j][2][4] = -dt * ty2 * 0.40e+00;
         b[i][j][3][0] = -dt * ty2 * (-(u[i][j - 1][k][2] * u[i][j - 1][k][3]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[i][j - 1][k][3]);
         b[i][j][3][1] = 0.0;
         b[i][j][3][2] = -dt * ty2 * (u[i][j - 1][k][3] * tmp1);
         b[i][j][3][3] = -dt * ty2 * (u[i][j - 1][k][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy4;
         b[i][j][3][4] = 0.0;
         b[i][j][4][0] = -dt * ty2 * ((0.40e+00 * (u[i][j - 1][k][1] * u[i][j - 1][k][1] + u[i][j - 1][k][2] * u[i][j - 1][k][2] + u[i][j - 1][k][3] * u[i][j - 1][k][3]) * tmp2 - 1.40e+00 * (u[i][j - 1][k][4] * tmp1)) * (u[i][j - 1][k][2] * tmp1)) - dt * ty1 * (-(c34 - c1345) * tmp3 * (((u[i][j - 1][k][1]) * (u[i][j - 1][k][1]))) - (r43 * c34 - c1345) * tmp3 * (((u[i][j - 1][k][2]) * (u[i][j - 1][k][2]))) - (c34 - c1345) * tmp3 * (((u[i][j - 1][k][3]) * (u[i][j - 1][k][3]))) - c1345 * tmp2 * u[i][j - 1][k][4]);
         b[i][j][4][1] = -dt * ty2 * (-0.40e+00 * (u[i][j - 1][k][1] * u[i][j - 1][k][2]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[i][j - 1][k][1];
         b[i][j][4][2] = -dt * ty2 * (1.40e+00 * (u[i][j - 1][k][4] * tmp1) - 0.50 * 0.40e+00 * ((u[i][j - 1][k][1] * u[i][j - 1][k][1] + 3.0 * u[i][j - 1][k][2] * u[i][j - 1][k][2] + u[i][j - 1][k][3] * u[i][j - 1][k][3]) * tmp2)) - dt * ty1 * (r43 * c34 - c1345) * tmp2 * u[i][j - 1][k][2];
         b[i][j][4][3] = -dt * ty2 * (-0.40e+00 * (u[i][j - 1][k][2] * u[i][j - 1][k][3]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[i][j - 1][k][3];
         b[i][j][4][4] = -dt * ty2 * (1.40e+00 * (u[i][j - 1][k][2] * tmp1)) - dt * ty1 * c1345 * tmp1 - dt * ty1 * dy5;
         tmp1 = 1.0 / u[i - 1][j][k][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         c[i][j][0][0] = -dt * tx1 * dx1;
         c[i][j][0][1] = -dt * tx2;
         c[i][j][0][2] = 0.0;
         c[i][j][0][3] = 0.0;
         c[i][j][0][4] = 0.0;
         c[i][j][1][0] = -dt * tx2 * (-(u[i - 1][j][k][1] * tmp1) * (u[i - 1][j][k][1] * tmp1) + 0.40e+00 * 0.50 * (u[i - 1][j][k][1] * u[i - 1][j][k][1] + u[i - 1][j][k][2] * u[i - 1][j][k][2] + u[i - 1][j][k][3] * u[i - 1][j][k][3]) * tmp2) - dt * tx1 * (-r43 * c34 * tmp2 * u[i - 1][j][k][1]);
         c[i][j][1][1] = -dt * tx2 * ((2.0 - 0.40e+00) * (u[i - 1][j][k][1] * tmp1)) - dt * tx1 * (r43 * c34 * tmp1) - dt * tx1 * dx2;
         c[i][j][1][2] = -dt * tx2 * (-0.40e+00 * (u[i - 1][j][k][2] * tmp1));
         c[i][j][1][3] = -dt * tx2 * (-0.40e+00 * (u[i - 1][j][k][3] * tmp1));
         c[i][j][1][4] = -dt * tx2 * 0.40e+00;
         c[i][j][2][0] = -dt * tx2 * (-(u[i - 1][j][k][1] * u[i - 1][j][k][2]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[i - 1][j][k][2]);
         c[i][j][2][1] = -dt * tx2 * (u[i - 1][j][k][2] * tmp1);
         c[i][j][2][2] = -dt * tx2 * (u[i - 1][j][k][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx3;
         c[i][j][2][3] = 0.0;
         c[i][j][2][4] = 0.0;
         c[i][j][3][0] = -dt * tx2 * (-(u[i - 1][j][k][1] * u[i - 1][j][k][3]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[i - 1][j][k][3]);
         c[i][j][3][1] = -dt * tx2 * (u[i - 1][j][k][3] * tmp1);
         c[i][j][3][2] = 0.0;
         c[i][j][3][3] = -dt * tx2 * (u[i - 1][j][k][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx4;
         c[i][j][3][4] = 0.0;
         c[i][j][4][0] = -dt * tx2 * ((0.40e+00 * (u[i - 1][j][k][1] * u[i - 1][j][k][1] + u[i - 1][j][k][2] * u[i - 1][j][k][2] + u[i - 1][j][k][3] * u[i - 1][j][k][3]) * tmp2 - 1.40e+00 * (u[i - 1][j][k][4] * tmp1)) * (u[i - 1][j][k][1] * tmp1)) - dt * tx1 * (-(r43 * c34 - c1345) * tmp3 * (((u[i - 1][j][k][1]) * (u[i - 1][j][k][1]))) - (c34 - c1345) * tmp3 * (((u[i - 1][j][k][2]) * (u[i - 1][j][k][2]))) - (c34 - c1345) * tmp3 * (((u[i - 1][j][k][3]) * (u[i - 1][j][k][3]))) - c1345 * tmp2 * u[i - 1][j][k][4]);
         c[i][j][4][1] = -dt * tx2 * (1.40e+00 * (u[i - 1][j][k][4] * tmp1) - 0.50 * 0.40e+00 * ((3.0 * u[i - 1][j][k][1] * u[i - 1][j][k][1] + u[i - 1][j][k][2] * u[i - 1][j][k][2] + u[i - 1][j][k][3] * u[i - 1][j][k][3]) * tmp2)) - dt * tx1 * (r43 * c34 - c1345) * tmp2 * u[i - 1][j][k][1];
         c[i][j][4][2] = -dt * tx2 * (-0.40e+00 * (u[i - 1][j][k][2] * u[i - 1][j][k][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[i - 1][j][k][2];
         c[i][j][4][3] = -dt * tx2 * (-0.40e+00 * (u[i - 1][j][k][3] * u[i - 1][j][k][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[i - 1][j][k][3];
         c[i][j][4][4] = -dt * tx2 * (1.40e+00 * (u[i - 1][j][k][1] * tmp1)) - dt * tx1 * c1345 * tmp1 - dt * tx1 * dx5;
      }
   }
}

static void jacu(int k) {
   int i, j;
   double r43;
   double c1345;
   double c34;
   double tmp1, tmp2, tmp3;
   r43 = (4.0 / 3.0);
   c1345 = 1.40e+00 * 1.00e-01 * 1.00e+00 * 1.40e+00;
   c34 = 1.00e-01 * 1.00e+00;
   #pragma omp parallel for default(shared) private(i, j, tmp1, tmp2, tmp3) firstprivate(iend, ist, jend, jst, k, tx1, dx1, ty1, dy1, tz1, dz1, dt, r43, c34, dx2, dy2, dz2, dx3, dy3, dz3, dx4, dy4, dz4, c1345, dx5, dy5, dz5, tx2, ty2, tz2, u)
   for(i = iend; i >= ist; i--) {
      #pragma omp parallel for default(shared) private(j, tmp1, tmp2, tmp3) firstprivate(jend, jst, i, k, tx1, dx1, ty1, dy1, tz1, dz1, dt, r43, c34, dx2, dy2, dz2, dx3, dy3, dz3, dx4, dy4, dz4, c1345, dx5, dy5, dz5, tx2, ty2, tz2, u)
      for(j = jend; j >= jst; j--) {
         tmp1 = 1.0 / u[i][j][k][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         d[i][j][0][0] = 1.0 + dt * 2.0 * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
         d[i][j][0][1] = 0.0;
         d[i][j][0][2] = 0.0;
         d[i][j][0][3] = 0.0;
         d[i][j][0][4] = 0.0;
         d[i][j][1][0] = dt * 2.0 * (tx1 * (-r43 * c34 * tmp2 * u[i][j][k][1]) + ty1 * (-c34 * tmp2 * u[i][j][k][1]) + tz1 * (-c34 * tmp2 * u[i][j][k][1]));
         d[i][j][1][1] = 1.0 + dt * 2.0 * (tx1 * r43 * c34 * tmp1 + ty1 * c34 * tmp1 + tz1 * c34 * tmp1) + dt * 2.0 * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
         d[i][j][1][2] = 0.0;
         d[i][j][1][3] = 0.0;
         d[i][j][1][4] = 0.0;
         d[i][j][2][0] = dt * 2.0 * (tx1 * (-c34 * tmp2 * u[i][j][k][2]) + ty1 * (-r43 * c34 * tmp2 * u[i][j][k][2]) + tz1 * (-c34 * tmp2 * u[i][j][k][2]));
         d[i][j][2][1] = 0.0;
         d[i][j][2][2] = 1.0 + dt * 2.0 * (tx1 * c34 * tmp1 + ty1 * r43 * c34 * tmp1 + tz1 * c34 * tmp1) + dt * 2.0 * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
         d[i][j][2][3] = 0.0;
         d[i][j][2][4] = 0.0;
         d[i][j][3][0] = dt * 2.0 * (tx1 * (-c34 * tmp2 * u[i][j][k][3]) + ty1 * (-c34 * tmp2 * u[i][j][k][3]) + tz1 * (-r43 * c34 * tmp2 * u[i][j][k][3]));
         d[i][j][3][1] = 0.0;
         d[i][j][3][2] = 0.0;
         d[i][j][3][3] = 1.0 + dt * 2.0 * (tx1 * c34 * tmp1 + ty1 * c34 * tmp1 + tz1 * r43 * c34 * tmp1) + dt * 2.0 * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
         d[i][j][3][4] = 0.0;
         d[i][j][4][0] = dt * 2.0 * (tx1 * (-(r43 * c34 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - (c1345) * tmp2 * u[i][j][k][4]) + ty1 * (-(c34 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (r43 * c34 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - (c1345) * tmp2 * u[i][j][k][4]) + tz1 * (-(c34 - c1345) * tmp3 * (((u[i][j][k][1]) * (u[i][j][k][1]))) - (c34 - c1345) * tmp3 * (((u[i][j][k][2]) * (u[i][j][k][2]))) - (r43 * c34 - c1345) * tmp3 * (((u[i][j][k][3]) * (u[i][j][k][3]))) - (c1345) * tmp2 * u[i][j][k][4]));
         d[i][j][4][1] = dt * 2.0 * (tx1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k][1] + ty1 * (c34 - c1345) * tmp2 * u[i][j][k][1] + tz1 * (c34 - c1345) * tmp2 * u[i][j][k][1]);
         d[i][j][4][2] = dt * 2.0 * (tx1 * (c34 - c1345) * tmp2 * u[i][j][k][2] + ty1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k][2] + tz1 * (c34 - c1345) * tmp2 * u[i][j][k][2]);
         d[i][j][4][3] = dt * 2.0 * (tx1 * (c34 - c1345) * tmp2 * u[i][j][k][3] + ty1 * (c34 - c1345) * tmp2 * u[i][j][k][3] + tz1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k][3]);
         d[i][j][4][4] = 1.0 + dt * 2.0 * (tx1 * c1345 * tmp1 + ty1 * c1345 * tmp1 + tz1 * c1345 * tmp1) + dt * 2.0 * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
         tmp1 = 1.0 / u[i + 1][j][k][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         a[i][j][0][0] = -dt * tx1 * dx1;
         a[i][j][0][1] = dt * tx2;
         a[i][j][0][2] = 0.0;
         a[i][j][0][3] = 0.0;
         a[i][j][0][4] = 0.0;
         a[i][j][1][0] = dt * tx2 * (-(u[i + 1][j][k][1] * tmp1) * (u[i + 1][j][k][1] * tmp1) + 0.40e+00 * 0.50 * (u[i + 1][j][k][1] * u[i + 1][j][k][1] + u[i + 1][j][k][2] * u[i + 1][j][k][2] + u[i + 1][j][k][3] * u[i + 1][j][k][3]) * tmp2) - dt * tx1 * (-r43 * c34 * tmp2 * u[i + 1][j][k][1]);
         a[i][j][1][1] = dt * tx2 * ((2.0 - 0.40e+00) * (u[i + 1][j][k][1] * tmp1)) - dt * tx1 * (r43 * c34 * tmp1) - dt * tx1 * dx2;
         a[i][j][1][2] = dt * tx2 * (-0.40e+00 * (u[i + 1][j][k][2] * tmp1));
         a[i][j][1][3] = dt * tx2 * (-0.40e+00 * (u[i + 1][j][k][3] * tmp1));
         a[i][j][1][4] = dt * tx2 * 0.40e+00;
         a[i][j][2][0] = dt * tx2 * (-(u[i + 1][j][k][1] * u[i + 1][j][k][2]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[i + 1][j][k][2]);
         a[i][j][2][1] = dt * tx2 * (u[i + 1][j][k][2] * tmp1);
         a[i][j][2][2] = dt * tx2 * (u[i + 1][j][k][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx3;
         a[i][j][2][3] = 0.0;
         a[i][j][2][4] = 0.0;
         a[i][j][3][0] = dt * tx2 * (-(u[i + 1][j][k][1] * u[i + 1][j][k][3]) * tmp2) - dt * tx1 * (-c34 * tmp2 * u[i + 1][j][k][3]);
         a[i][j][3][1] = dt * tx2 * (u[i + 1][j][k][3] * tmp1);
         a[i][j][3][2] = 0.0;
         a[i][j][3][3] = dt * tx2 * (u[i + 1][j][k][1] * tmp1) - dt * tx1 * (c34 * tmp1) - dt * tx1 * dx4;
         a[i][j][3][4] = 0.0;
         a[i][j][4][0] = dt * tx2 * ((0.40e+00 * (u[i + 1][j][k][1] * u[i + 1][j][k][1] + u[i + 1][j][k][2] * u[i + 1][j][k][2] + u[i + 1][j][k][3] * u[i + 1][j][k][3]) * tmp2 - 1.40e+00 * (u[i + 1][j][k][4] * tmp1)) * (u[i + 1][j][k][1] * tmp1)) - dt * tx1 * (-(r43 * c34 - c1345) * tmp3 * (((u[i + 1][j][k][1]) * (u[i + 1][j][k][1]))) - (c34 - c1345) * tmp3 * (((u[i + 1][j][k][2]) * (u[i + 1][j][k][2]))) - (c34 - c1345) * tmp3 * (((u[i + 1][j][k][3]) * (u[i + 1][j][k][3]))) - c1345 * tmp2 * u[i + 1][j][k][4]);
         a[i][j][4][1] = dt * tx2 * (1.40e+00 * (u[i + 1][j][k][4] * tmp1) - 0.50 * 0.40e+00 * ((3.0 * u[i + 1][j][k][1] * u[i + 1][j][k][1] + u[i + 1][j][k][2] * u[i + 1][j][k][2] + u[i + 1][j][k][3] * u[i + 1][j][k][3]) * tmp2)) - dt * tx1 * (r43 * c34 - c1345) * tmp2 * u[i + 1][j][k][1];
         a[i][j][4][2] = dt * tx2 * (-0.40e+00 * (u[i + 1][j][k][2] * u[i + 1][j][k][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[i + 1][j][k][2];
         a[i][j][4][3] = dt * tx2 * (-0.40e+00 * (u[i + 1][j][k][3] * u[i + 1][j][k][1]) * tmp2) - dt * tx1 * (c34 - c1345) * tmp2 * u[i + 1][j][k][3];
         a[i][j][4][4] = dt * tx2 * (1.40e+00 * (u[i + 1][j][k][1] * tmp1)) - dt * tx1 * c1345 * tmp1 - dt * tx1 * dx5;
         tmp1 = 1.0 / u[i][j + 1][k][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         b[i][j][0][0] = -dt * ty1 * dy1;
         b[i][j][0][1] = 0.0;
         b[i][j][0][2] = dt * ty2;
         b[i][j][0][3] = 0.0;
         b[i][j][0][4] = 0.0;
         b[i][j][1][0] = dt * ty2 * (-(u[i][j + 1][k][1] * u[i][j + 1][k][2]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[i][j + 1][k][1]);
         b[i][j][1][1] = dt * ty2 * (u[i][j + 1][k][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy2;
         b[i][j][1][2] = dt * ty2 * (u[i][j + 1][k][1] * tmp1);
         b[i][j][1][3] = 0.0;
         b[i][j][1][4] = 0.0;
         b[i][j][2][0] = dt * ty2 * (-(u[i][j + 1][k][2] * tmp1) * (u[i][j + 1][k][2] * tmp1) + 0.50 * 0.40e+00 * ((u[i][j + 1][k][1] * u[i][j + 1][k][1] + u[i][j + 1][k][2] * u[i][j + 1][k][2] + u[i][j + 1][k][3] * u[i][j + 1][k][3]) * tmp2)) - dt * ty1 * (-r43 * c34 * tmp2 * u[i][j + 1][k][2]);
         b[i][j][2][1] = dt * ty2 * (-0.40e+00 * (u[i][j + 1][k][1] * tmp1));
         b[i][j][2][2] = dt * ty2 * ((2.0 - 0.40e+00) * (u[i][j + 1][k][2] * tmp1)) - dt * ty1 * (r43 * c34 * tmp1) - dt * ty1 * dy3;
         b[i][j][2][3] = dt * ty2 * (-0.40e+00 * (u[i][j + 1][k][3] * tmp1));
         b[i][j][2][4] = dt * ty2 * 0.40e+00;
         b[i][j][3][0] = dt * ty2 * (-(u[i][j + 1][k][2] * u[i][j + 1][k][3]) * tmp2) - dt * ty1 * (-c34 * tmp2 * u[i][j + 1][k][3]);
         b[i][j][3][1] = 0.0;
         b[i][j][3][2] = dt * ty2 * (u[i][j + 1][k][3] * tmp1);
         b[i][j][3][3] = dt * ty2 * (u[i][j + 1][k][2] * tmp1) - dt * ty1 * (c34 * tmp1) - dt * ty1 * dy4;
         b[i][j][3][4] = 0.0;
         b[i][j][4][0] = dt * ty2 * ((0.40e+00 * (u[i][j + 1][k][1] * u[i][j + 1][k][1] + u[i][j + 1][k][2] * u[i][j + 1][k][2] + u[i][j + 1][k][3] * u[i][j + 1][k][3]) * tmp2 - 1.40e+00 * (u[i][j + 1][k][4] * tmp1)) * (u[i][j + 1][k][2] * tmp1)) - dt * ty1 * (-(c34 - c1345) * tmp3 * (((u[i][j + 1][k][1]) * (u[i][j + 1][k][1]))) - (r43 * c34 - c1345) * tmp3 * (((u[i][j + 1][k][2]) * (u[i][j + 1][k][2]))) - (c34 - c1345) * tmp3 * (((u[i][j + 1][k][3]) * (u[i][j + 1][k][3]))) - c1345 * tmp2 * u[i][j + 1][k][4]);
         b[i][j][4][1] = dt * ty2 * (-0.40e+00 * (u[i][j + 1][k][1] * u[i][j + 1][k][2]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[i][j + 1][k][1];
         b[i][j][4][2] = dt * ty2 * (1.40e+00 * (u[i][j + 1][k][4] * tmp1) - 0.50 * 0.40e+00 * ((u[i][j + 1][k][1] * u[i][j + 1][k][1] + 3.0 * u[i][j + 1][k][2] * u[i][j + 1][k][2] + u[i][j + 1][k][3] * u[i][j + 1][k][3]) * tmp2)) - dt * ty1 * (r43 * c34 - c1345) * tmp2 * u[i][j + 1][k][2];
         b[i][j][4][3] = dt * ty2 * (-0.40e+00 * (u[i][j + 1][k][2] * u[i][j + 1][k][3]) * tmp2) - dt * ty1 * (c34 - c1345) * tmp2 * u[i][j + 1][k][3];
         b[i][j][4][4] = dt * ty2 * (1.40e+00 * (u[i][j + 1][k][2] * tmp1)) - dt * ty1 * c1345 * tmp1 - dt * ty1 * dy5;
         tmp1 = 1.0 / u[i][j][k + 1][0];
         tmp2 = tmp1 * tmp1;
         tmp3 = tmp1 * tmp2;
         c[i][j][0][0] = -dt * tz1 * dz1;
         c[i][j][0][1] = 0.0;
         c[i][j][0][2] = 0.0;
         c[i][j][0][3] = dt * tz2;
         c[i][j][0][4] = 0.0;
         c[i][j][1][0] = dt * tz2 * (-(u[i][j][k + 1][1] * u[i][j][k + 1][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[i][j][k + 1][1]);
         c[i][j][1][1] = dt * tz2 * (u[i][j][k + 1][3] * tmp1) - dt * tz1 * c34 * tmp1 - dt * tz1 * dz2;
         c[i][j][1][2] = 0.0;
         c[i][j][1][3] = dt * tz2 * (u[i][j][k + 1][1] * tmp1);
         c[i][j][1][4] = 0.0;
         c[i][j][2][0] = dt * tz2 * (-(u[i][j][k + 1][2] * u[i][j][k + 1][3]) * tmp2) - dt * tz1 * (-c34 * tmp2 * u[i][j][k + 1][2]);
         c[i][j][2][1] = 0.0;
         c[i][j][2][2] = dt * tz2 * (u[i][j][k + 1][3] * tmp1) - dt * tz1 * (c34 * tmp1) - dt * tz1 * dz3;
         c[i][j][2][3] = dt * tz2 * (u[i][j][k + 1][2] * tmp1);
         c[i][j][2][4] = 0.0;
         c[i][j][3][0] = dt * tz2 * (-(u[i][j][k + 1][3] * tmp1) * (u[i][j][k + 1][3] * tmp1) + 0.50 * 0.40e+00 * ((u[i][j][k + 1][1] * u[i][j][k + 1][1] + u[i][j][k + 1][2] * u[i][j][k + 1][2] + u[i][j][k + 1][3] * u[i][j][k + 1][3]) * tmp2)) - dt * tz1 * (-r43 * c34 * tmp2 * u[i][j][k + 1][3]);
         c[i][j][3][1] = dt * tz2 * (-0.40e+00 * (u[i][j][k + 1][1] * tmp1));
         c[i][j][3][2] = dt * tz2 * (-0.40e+00 * (u[i][j][k + 1][2] * tmp1));
         c[i][j][3][3] = dt * tz2 * (2.0 - 0.40e+00) * (u[i][j][k + 1][3] * tmp1) - dt * tz1 * (r43 * c34 * tmp1) - dt * tz1 * dz4;
         c[i][j][3][4] = dt * tz2 * 0.40e+00;
         c[i][j][4][0] = dt * tz2 * ((0.40e+00 * (u[i][j][k + 1][1] * u[i][j][k + 1][1] + u[i][j][k + 1][2] * u[i][j][k + 1][2] + u[i][j][k + 1][3] * u[i][j][k + 1][3]) * tmp2 - 1.40e+00 * (u[i][j][k + 1][4] * tmp1)) * (u[i][j][k + 1][3] * tmp1)) - dt * tz1 * (-(c34 - c1345) * tmp3 * (((u[i][j][k + 1][1]) * (u[i][j][k + 1][1]))) - (c34 - c1345) * tmp3 * (((u[i][j][k + 1][2]) * (u[i][j][k + 1][2]))) - (r43 * c34 - c1345) * tmp3 * (((u[i][j][k + 1][3]) * (u[i][j][k + 1][3]))) - c1345 * tmp2 * u[i][j][k + 1][4]);
         c[i][j][4][1] = dt * tz2 * (-0.40e+00 * (u[i][j][k + 1][1] * u[i][j][k + 1][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[i][j][k + 1][1];
         c[i][j][4][2] = dt * tz2 * (-0.40e+00 * (u[i][j][k + 1][2] * u[i][j][k + 1][3]) * tmp2) - dt * tz1 * (c34 - c1345) * tmp2 * u[i][j][k + 1][2];
         c[i][j][4][3] = dt * tz2 * (1.40e+00 * (u[i][j][k + 1][4] * tmp1) - 0.50 * 0.40e+00 * ((u[i][j][k + 1][1] * u[i][j][k + 1][1] + u[i][j][k + 1][2] * u[i][j][k + 1][2] + 3.0 * u[i][j][k + 1][3] * u[i][j][k + 1][3]) * tmp2)) - dt * tz1 * (r43 * c34 - c1345) * tmp2 * u[i][j][k + 1][3];
         c[i][j][4][4] = dt * tz2 * (1.40e+00 * (u[i][j][k + 1][3] * tmp1)) - dt * tz1 * c1345 * tmp1 - dt * tz1 * dz5;
      }
   }
}

static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend, double v[64][65][65][5], double sum[5]) {
   {
      int i, j, k, m;
      double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(m = 0; m < 5; m++) {
         sum[m] = 0.0;
      }
      #pragma omp parallel for default(shared) private(i, j, k) firstprivate(ist, iend, jst, jend, nz0, v)
      for(i = ist; i <= iend; i++) {
         #pragma omp parallel for default(shared) private(j, k) firstprivate(jst, jend, nz0, v)
         for(j = jst; j <= jend; j++) {
            #pragma omp parallel for default(shared) private(k) firstprivate(nz0, v)
            for(k = 1; k <= nz0 - 2; k++) {
               sum0 = sum0 + v[i][j][k][0] * v[i][j][k][0];
               sum1 = sum1 + v[i][j][k][1] * v[i][j][k][1];
               sum2 = sum2 + v[i][j][k][2] * v[i][j][k][2];
               sum3 = sum3 + v[i][j][k][3] * v[i][j][k][3];
               sum4 = sum4 + v[i][j][k][4] * v[i][j][k][4];
            }
         }
      }
      {
         sum[0] += sum0;
         sum[1] += sum1;
         sum[2] += sum2;
         sum[3] += sum3;
         sum[4] += sum4;
      }
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(m = 0; m < 5; m++) {
         sum[m] = sqrt(sum[m] / ((nx0 - 2) * (ny0 - 2) * (nz0 - 2)));
      }
   }
}

static void pintgr() {
   int i, j, k;
   int ibeg, ifin, ifin1;
   int jbeg, jfin, jfin1;
   int iglob, iglob1, iglob2;
   int jglob, jglob1, jglob2;
   double phi1[66][66];
   double phi2[66][66];
   double frc1, frc2, frc3;
   ibeg = nx;
   ifin = 0;
   iglob1 = -1;
   iglob2 = nx - 1;
   if(iglob1 >= ii1 && iglob2 < ii2 + nx) ibeg = 0;
   if(iglob1 >= ii1 - nx && iglob2 <= ii2) ifin = nx;
   if(ii1 >= iglob1 && ii1 <= iglob2) ibeg = ii1;
   if(ii2 >= iglob1 && ii2 <= iglob2) ifin = ii2;
   jbeg = ny;
   jfin = -1;
   jglob1 = 0;
   jglob2 = ny - 1;
   if(jglob1 >= ji1 && jglob2 < ji2 + ny) jbeg = 0;
   if(jglob1 > ji1 - ny && jglob2 <= ji2) jfin = ny;
   if(ji1 >= jglob1 && ji1 <= jglob2) jbeg = ji1;
   if(ji2 >= jglob1 && ji2 <= jglob2) jfin = ji2;
   ifin1 = ifin;
   jfin1 = jfin;
   if(ifin1 == ii2) ifin1 = ifin - 1;
   if(jfin1 == ji2) jfin1 = jfin - 1;
   #pragma omp parallel for default(shared) private(i, k)
   for(i = 0; i <= 64 + 1; i++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(i)
      for(k = 0; k <= 64 + 1; k++) {
         phi1[i][k] = 0.0;
         phi2[i][k] = 0.0;
      }
   }
   #pragma omp parallel for default(shared) private(i, j, iglob, jglob, k) firstprivate(ibeg, ifin, jbeg, jfin, ki1, ki2, u)
   for(i = ibeg; i <= ifin; i++) {
      iglob = i;
      #pragma omp parallel for default(shared) private(j, jglob, k) firstprivate(jbeg, jfin, ki1, i, ki2, u)
      for(j = jbeg; j <= jfin; j++) {
         jglob = j;
         k = ki1;
         phi1[i][j] = 0.40e+00 * (u[i][j][k][4] - 0.50 * (((u[i][j][k][1]) * (u[i][j][k][1])) + ((u[i][j][k][2]) * (u[i][j][k][2])) + ((u[i][j][k][3]) * (u[i][j][k][3]))) / u[i][j][k][0]);
         k = ki2;
         phi2[i][j] = 0.40e+00 * (u[i][j][k][4] - 0.50 * (((u[i][j][k][1]) * (u[i][j][k][1])) + ((u[i][j][k][2]) * (u[i][j][k][2])) + ((u[i][j][k][3]) * (u[i][j][k][3]))) / u[i][j][k][0]);
      }
   }
   frc1 = 0.0;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(ibeg, ifin1, jbeg, jfin1, phi1, phi2) reduction(+ : frc1)
   for(i = ibeg; i <= ifin1; i++) {
      #pragma omp parallel for default(shared) private(j) firstprivate(jbeg, jfin1, i, phi1, phi2) reduction(+ : frc1)
      for(j = jbeg; j <= jfin1; j++) {
         frc1 = frc1 + (phi1[i][j] + phi1[i + 1][j] + phi1[i][j + 1] + phi1[i + 1][j + 1] + phi2[i][j] + phi2[i + 1][j] + phi2[i][j + 1] + phi2[i + 1][j + 1]);
      }
   }
   frc1 = dxi * deta * frc1;
   #pragma omp parallel for default(shared) private(i, k)
   for(i = 0; i <= 64 + 1; i++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(i)
      for(k = 0; k <= 64 + 1; k++) {
         phi1[i][k] = 0.0;
         phi2[i][k] = 0.0;
      }
   }
   jglob = jbeg;
   if(jglob == ji1) {
      #pragma omp parallel for default(shared) private(i, k, iglob) firstprivate(ibeg, ifin, ki1, ki2, jbeg, u)
      for(i = ibeg; i <= ifin; i++) {
         iglob = i;
         #pragma omp parallel for default(shared) private(k) firstprivate(ki1, ki2, i, jbeg, u)
         for(k = ki1; k <= ki2; k++) {
            phi1[i][k] = 0.40e+00 * (u[i][jbeg][k][4] - 0.50 * (((u[i][jbeg][k][1]) * (u[i][jbeg][k][1])) + ((u[i][jbeg][k][2]) * (u[i][jbeg][k][2])) + ((u[i][jbeg][k][3]) * (u[i][jbeg][k][3]))) / u[i][jbeg][k][0]);
         }
      }
   }
   jglob = jfin;
   if(jglob == ji2) {
      #pragma omp parallel for default(shared) private(i, k, iglob) firstprivate(ibeg, ifin, ki1, ki2, jfin, u)
      for(i = ibeg; i <= ifin; i++) {
         iglob = i;
         #pragma omp parallel for default(shared) private(k) firstprivate(ki1, ki2, i, jfin, u)
         for(k = ki1; k <= ki2; k++) {
            phi2[i][k] = 0.40e+00 * (u[i][jfin][k][4] - 0.50 * (((u[i][jfin][k][1]) * (u[i][jfin][k][1])) + ((u[i][jfin][k][2]) * (u[i][jfin][k][2])) + ((u[i][jfin][k][3]) * (u[i][jfin][k][3]))) / u[i][jfin][k][0]);
         }
      }
   }
   frc2 = 0.0;
   #pragma omp parallel for default(shared) private(i, k) firstprivate(ibeg, ifin1, ki1, ki2, phi1, phi2) reduction(+ : frc2)
   for(i = ibeg; i <= ifin1; i++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(ki1, ki2, i, phi1, phi2) reduction(+ : frc2)
      for(k = ki1; k <= ki2 - 1; k++) {
         frc2 = frc2 + (phi1[i][k] + phi1[i + 1][k] + phi1[i][k + 1] + phi1[i + 1][k + 1] + phi2[i][k] + phi2[i + 1][k] + phi2[i][k + 1] + phi2[i + 1][k + 1]);
      }
   }
   frc2 = dxi * dzeta * frc2;
   #pragma omp parallel for default(shared) private(i, k)
   for(i = 0; i <= 64 + 1; i++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(i)
      for(k = 0; k <= 64 + 1; k++) {
         phi1[i][k] = 0.0;
         phi2[i][k] = 0.0;
      }
   }
   iglob = ibeg;
   if(iglob == ii1) {
      #pragma omp parallel for default(shared) private(j, k, jglob) firstprivate(jbeg, jfin, ki1, ki2, ibeg, u)
      for(j = jbeg; j <= jfin; j++) {
         jglob = j;
         #pragma omp parallel for default(shared) private(k) firstprivate(ki1, ki2, ibeg, j, u)
         for(k = ki1; k <= ki2; k++) {
            phi1[j][k] = 0.40e+00 * (u[ibeg][j][k][4] - 0.50 * (((u[ibeg][j][k][1]) * (u[ibeg][j][k][1])) + ((u[ibeg][j][k][2]) * (u[ibeg][j][k][2])) + ((u[ibeg][j][k][3]) * (u[ibeg][j][k][3]))) / u[ibeg][j][k][0]);
         }
      }
   }
   iglob = ifin;
   if(iglob == ii2) {
      #pragma omp parallel for default(shared) private(j, k, jglob) firstprivate(jbeg, jfin, ki1, ki2, ifin, u)
      for(j = jbeg; j <= jfin; j++) {
         jglob = j;
         #pragma omp parallel for default(shared) private(k) firstprivate(ki1, ki2, ifin, j, u)
         for(k = ki1; k <= ki2; k++) {
            phi2[j][k] = 0.40e+00 * (u[ifin][j][k][4] - 0.50 * (((u[ifin][j][k][1]) * (u[ifin][j][k][1])) + ((u[ifin][j][k][2]) * (u[ifin][j][k][2])) + ((u[ifin][j][k][3]) * (u[ifin][j][k][3]))) / u[ifin][j][k][0]);
         }
      }
   }
   frc3 = 0.0;
   #pragma omp parallel for default(shared) private(j, k) firstprivate(jbeg, jfin1, ki1, ki2, phi1, phi2) reduction(+ : frc3)
   for(j = jbeg; j <= jfin1; j++) {
      #pragma omp parallel for default(shared) private(k) firstprivate(ki1, ki2, j, phi1, phi2) reduction(+ : frc3)
      for(k = ki1; k <= ki2 - 1; k++) {
         frc3 = frc3 + (phi1[j][k] + phi1[j + 1][k] + phi1[j][k + 1] + phi1[j + 1][k + 1] + phi2[j][k] + phi2[j + 1][k] + phi2[j][k + 1] + phi2[j + 1][k + 1]);
      }
   }
   frc3 = deta * dzeta * frc3;
   frc = 0.25 * (frc1 + frc2 + frc3);
}

static void read_input() {
   FILE *fp;
   printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - LU Benchmark\n\n");
   fp = fopen("inputlu.data", "r");
   if(fp != ((void *) 0)) {
      printf(" Reading from input file inputlu.data\n");
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      fscanf(fp, "%d%d", &ipr, &inorm);
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      fscanf(fp, "%d", &itmax);
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      fscanf(fp, "%lf", &dt);
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      fscanf(fp, "%lf", &omega);
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      fscanf(fp, "%lf%lf%lf%lf%lf", &tolrsd[0], &tolrsd[1], &tolrsd[2], &tolrsd[3], &tolrsd[4]);
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      while(fgetc(fp) != '\n');
      fscanf(fp, "%d%d%d", &nx0, &ny0, &nz0);
      while(fgetc(fp) != '\n');
      fclose(fp);
   }
   else {
      ipr = 1;
      inorm = 250;
      itmax = 250;
      dt = 2.0;
      omega = 1.2;
      tolrsd[0] = 1.0e-8;
      tolrsd[1] = 1.0e-8;
      tolrsd[2] = 1.0e-8;
      tolrsd[3] = 1.0e-8;
      tolrsd[4] = 1.0e-8;
      nx0 = 64;
      ny0 = 64;
      nz0 = 64;
   }
   if(nx0 < 4 || ny0 < 4 || nz0 < 4) {
      printf("     PROBLEM SIZE IS TOO SMALL - \n     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n");
      exit(1);
   }
   if(nx0 > 64 || ny0 > 64 || nz0 > 64) {
      printf("     PROBLEM SIZE IS TOO LARGE - \n     NX, NY AND NZ SHOULD BE EQUAL TO \n     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY\n");
      exit(1);
   }
   printf(" Size: %3dx%3dx%3d\n", nx0, ny0, nz0);
   printf(" Iterations: %3d\n", itmax);
}

static void rhs() {
   {
      int i, j, k, m;
      int L1, L2;
      int ist1, iend1;
      int jst1, jend1;
      double q;
      double u21, u31, u41;
      double tmp;
      double u21i, u31i, u41i, u51i;
      double u21j, u31j, u41j, u51j;
      double u21k, u31k, u41k, u51k;
      double u21im1, u31im1, u41im1, u51im1;
      double u21jm1, u31jm1, u41jm1, u51jm1;
      double u21km1, u31km1, u41km1, u51km1;
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(nx, ny, nz, frct)
      for(i = 0; i <= nx - 1; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(ny, nz, frct)
         for(j = 0; j <= ny - 1; j++) {
            #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, frct)
            for(k = 0; k <= nz - 1; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = -frct[i][j][k][m];
               }
            }
         }
      }
      L1 = 0;
      L2 = nx - 1;
      #pragma omp parallel for default(shared) private(i, j, k) firstprivate(jst, jend, nz, u)
      for(i = L1; i <= L2; i++) {
         #pragma omp parallel for default(shared) private(j, k) firstprivate(jst, jend, nz, u)
         for(j = jst; j <= jend; j++) {
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, u)
            for(k = 1; k <= nz - 2; k++) {
               flux[i][j][k][0] = u[i][j][k][1];
               u21 = u[i][j][k][1] / u[i][j][k][0];
               q = 0.50 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) / u[i][j][k][0];
               flux[i][j][k][1] = u[i][j][k][1] * u21 + 0.40e+00 * (u[i][j][k][4] - q);
               flux[i][j][k][2] = u[i][j][k][2] * u21;
               flux[i][j][k][3] = u[i][j][k][3] * u21;
               flux[i][j][k][4] = (1.40e+00 * u[i][j][k][4] - 0.40e+00 * q) * u21;
            }
         }
      }
      #pragma omp parallel for default(shared) private(j, k, i, m) firstprivate(jst, jend, nz, ist, iend, tx2, nx, tx3, dx1, tx1, dx2, dx3, dx4, dx5, dssp, u)
      for(j = jst; j <= jend; j++) {
         #pragma omp parallel for default(shared) private(k, i, m) firstprivate(nz, ist, iend, tx2, nx, tx3, dx1, tx1, dx2, dx3, dx4, dx5, dssp, u)
         for(k = 1; k <= nz - 2; k++) {
            #pragma omp parallel for default(shared) private(i, m) firstprivate(ist, iend, tx2, flux)
            for(i = ist; i <= iend; i++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = rsd[i][j][k][m] - tx2 * (flux[i + 1][j][k][m] - flux[i - 1][j][k][m]);
               }
            }
            L2 = nx - 1;
            #pragma omp parallel for default(shared) private(i) firstprivate(ist, tx3, u)
            for(i = ist; i <= L2; i++) {
               tmp = 1.0 / u[i][j][k][0];
               u21i = tmp * u[i][j][k][1];
               u31i = tmp * u[i][j][k][2];
               u41i = tmp * u[i][j][k][3];
               u51i = tmp * u[i][j][k][4];
               tmp = 1.0 / u[i - 1][j][k][0];
               u21im1 = tmp * u[i - 1][j][k][1];
               u31im1 = tmp * u[i - 1][j][k][2];
               u41im1 = tmp * u[i - 1][j][k][3];
               u51im1 = tmp * u[i - 1][j][k][4];
               flux[i][j][k][1] = (4.0 / 3.0) * tx3 * (u21i - u21im1);
               flux[i][j][k][2] = tx3 * (u31i - u31im1);
               flux[i][j][k][3] = tx3 * (u41i - u41im1);
               flux[i][j][k][4] = 0.50 * (1.0 - 1.40e+00 * 1.40e+00) * tx3 * ((((u21i) * (u21i)) + ((u31i) * (u31i)) + ((u41i) * (u41i))) - (((u21im1) * (u21im1)) + ((u31im1) * (u31im1)) + ((u41im1) * (u41im1)))) + (1.0 / 6.0) * tx3 * (((u21i) * (u21i)) - ((u21im1) * (u21im1))) + 1.40e+00 * 1.40e+00 * tx3 * (u51i - u51im1);
            }
            #pragma omp parallel for default(shared) private(i) firstprivate(ist, iend, dx1, tx1, tx3, dx2, dx3, dx4, dx5, u, flux)
            for(i = ist; i <= iend; i++) {
               rsd[i][j][k][0] = rsd[i][j][k][0] + dx1 * tx1 * (u[i - 1][j][k][0] - 2.0 * u[i][j][k][0] + u[i + 1][j][k][0]);
               rsd[i][j][k][1] = rsd[i][j][k][1] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][1] - flux[i][j][k][1]) + dx2 * tx1 * (u[i - 1][j][k][1] - 2.0 * u[i][j][k][1] + u[i + 1][j][k][1]);
               rsd[i][j][k][2] = rsd[i][j][k][2] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][2] - flux[i][j][k][2]) + dx3 * tx1 * (u[i - 1][j][k][2] - 2.0 * u[i][j][k][2] + u[i + 1][j][k][2]);
               rsd[i][j][k][3] = rsd[i][j][k][3] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][3] - flux[i][j][k][3]) + dx4 * tx1 * (u[i - 1][j][k][3] - 2.0 * u[i][j][k][3] + u[i + 1][j][k][3]);
               rsd[i][j][k][4] = rsd[i][j][k][4] + tx3 * 1.00e-01 * 1.00e+00 * (flux[i + 1][j][k][4] - flux[i][j][k][4]) + dx5 * tx1 * (u[i - 1][j][k][4] - 2.0 * u[i][j][k][4] + u[i + 1][j][k][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rsd[1][j][k][m] = rsd[1][j][k][m] - dssp * (+5.0 * u[1][j][k][m] - 4.0 * u[2][j][k][m] + u[3][j][k][m]);
               rsd[2][j][k][m] = rsd[2][j][k][m] - dssp * (-4.0 * u[1][j][k][m] + 6.0 * u[2][j][k][m] - 4.0 * u[3][j][k][m] + u[4][j][k][m]);
            }
            ist1 = 3;
            iend1 = nx - 4;
            #pragma omp parallel for default(shared) private(i, m) firstprivate(dssp, u)
            for(i = ist1; i <= iend1; i++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = rsd[i][j][k][m] - dssp * (u[i - 2][j][k][m] - 4.0 * u[i - 1][j][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i + 1][j][k][m] + u[i + 2][j][k][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rsd[nx - 3][j][k][m] = rsd[nx - 3][j][k][m] - dssp * (u[nx - 5][j][k][m] - 4.0 * u[nx - 4][j][k][m] + 6.0 * u[nx - 3][j][k][m] - 4.0 * u[nx - 2][j][k][m]);
               rsd[nx - 2][j][k][m] = rsd[nx - 2][j][k][m] - dssp * (u[nx - 4][j][k][m] - 4.0 * u[nx - 3][j][k][m] + 5.0 * u[nx - 2][j][k][m]);
            }
         }
      }
      L1 = 0;
      L2 = ny - 1;
      #pragma omp parallel for default(shared) private(i, j, k) firstprivate(ist, iend, nz, u)
      for(i = ist; i <= iend; i++) {
         #pragma omp parallel for default(shared) private(j, k) firstprivate(nz, u)
         for(j = L1; j <= L2; j++) {
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, u)
            for(k = 1; k <= nz - 2; k++) {
               flux[i][j][k][0] = u[i][j][k][2];
               u31 = u[i][j][k][2] / u[i][j][k][0];
               q = 0.50 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) / u[i][j][k][0];
               flux[i][j][k][1] = u[i][j][k][1] * u31;
               flux[i][j][k][2] = u[i][j][k][2] * u31 + 0.40e+00 * (u[i][j][k][4] - q);
               flux[i][j][k][3] = u[i][j][k][3] * u31;
               flux[i][j][k][4] = (1.40e+00 * u[i][j][k][4] - 0.40e+00 * q) * u31;
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, k, j, m) firstprivate(ist, iend, nz, jst, jend, ty2, ny, ty3, dy1, ty1, dy2, dy3, dy4, dy5, dssp, u)
      for(i = ist; i <= iend; i++) {
         #pragma omp parallel for default(shared) private(k, j, m) firstprivate(nz, jst, jend, ty2, ny, ty3, dy1, ty1, dy2, dy3, dy4, dy5, dssp, u)
         for(k = 1; k <= nz - 2; k++) {
            #pragma omp parallel for default(shared) private(j, m) firstprivate(jst, jend, ty2, flux)
            for(j = jst; j <= jend; j++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = rsd[i][j][k][m] - ty2 * (flux[i][j + 1][k][m] - flux[i][j - 1][k][m]);
               }
            }
            L2 = ny - 1;
            #pragma omp parallel for default(shared) private(j) firstprivate(jst, ty3, u)
            for(j = jst; j <= L2; j++) {
               tmp = 1.0 / u[i][j][k][0];
               u21j = tmp * u[i][j][k][1];
               u31j = tmp * u[i][j][k][2];
               u41j = tmp * u[i][j][k][3];
               u51j = tmp * u[i][j][k][4];
               tmp = 1.0 / u[i][j - 1][k][0];
               u21jm1 = tmp * u[i][j - 1][k][1];
               u31jm1 = tmp * u[i][j - 1][k][2];
               u41jm1 = tmp * u[i][j - 1][k][3];
               u51jm1 = tmp * u[i][j - 1][k][4];
               flux[i][j][k][1] = ty3 * (u21j - u21jm1);
               flux[i][j][k][2] = (4.0 / 3.0) * ty3 * (u31j - u31jm1);
               flux[i][j][k][3] = ty3 * (u41j - u41jm1);
               flux[i][j][k][4] = 0.50 * (1.0 - 1.40e+00 * 1.40e+00) * ty3 * ((((u21j) * (u21j)) + ((u31j) * (u31j)) + ((u41j) * (u41j))) - (((u21jm1) * (u21jm1)) + ((u31jm1) * (u31jm1)) + ((u41jm1) * (u41jm1)))) + (1.0 / 6.0) * ty3 * (((u31j) * (u31j)) - ((u31jm1) * (u31jm1))) + 1.40e+00 * 1.40e+00 * ty3 * (u51j - u51jm1);
            }
            #pragma omp parallel for default(shared) private(j) firstprivate(jst, jend, dy1, ty1, ty3, dy2, dy3, dy4, dy5, u, flux)
            for(j = jst; j <= jend; j++) {
               rsd[i][j][k][0] = rsd[i][j][k][0] + dy1 * ty1 * (u[i][j - 1][k][0] - 2.0 * u[i][j][k][0] + u[i][j + 1][k][0]);
               rsd[i][j][k][1] = rsd[i][j][k][1] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][1] - flux[i][j][k][1]) + dy2 * ty1 * (u[i][j - 1][k][1] - 2.0 * u[i][j][k][1] + u[i][j + 1][k][1]);
               rsd[i][j][k][2] = rsd[i][j][k][2] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][2] - flux[i][j][k][2]) + dy3 * ty1 * (u[i][j - 1][k][2] - 2.0 * u[i][j][k][2] + u[i][j + 1][k][2]);
               rsd[i][j][k][3] = rsd[i][j][k][3] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][3] - flux[i][j][k][3]) + dy4 * ty1 * (u[i][j - 1][k][3] - 2.0 * u[i][j][k][3] + u[i][j + 1][k][3]);
               rsd[i][j][k][4] = rsd[i][j][k][4] + ty3 * 1.00e-01 * 1.00e+00 * (flux[i][j + 1][k][4] - flux[i][j][k][4]) + dy5 * ty1 * (u[i][j - 1][k][4] - 2.0 * u[i][j][k][4] + u[i][j + 1][k][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rsd[i][1][k][m] = rsd[i][1][k][m] - dssp * (+5.0 * u[i][1][k][m] - 4.0 * u[i][2][k][m] + u[i][3][k][m]);
               rsd[i][2][k][m] = rsd[i][2][k][m] - dssp * (-4.0 * u[i][1][k][m] + 6.0 * u[i][2][k][m] - 4.0 * u[i][3][k][m] + u[i][4][k][m]);
            }
            jst1 = 3;
            jend1 = ny - 4;
            #pragma omp parallel for default(shared) private(j, m) firstprivate(dssp, u)
            for(j = jst1; j <= jend1; j++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = rsd[i][j][k][m] - dssp * (u[i][j - 2][k][m] - 4.0 * u[i][j - 1][k][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j + 1][k][m] + u[i][j + 2][k][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rsd[i][ny - 3][k][m] = rsd[i][ny - 3][k][m] - dssp * (u[i][ny - 5][k][m] - 4.0 * u[i][ny - 4][k][m] + 6.0 * u[i][ny - 3][k][m] - 4.0 * u[i][ny - 2][k][m]);
               rsd[i][ny - 2][k][m] = rsd[i][ny - 2][k][m] - dssp * (u[i][ny - 4][k][m] - 4.0 * u[i][ny - 3][k][m] + 5.0 * u[i][ny - 2][k][m]);
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(ist, iend, jst, jend, nz, tz2, tz3, dz1, tz1, dz2, dz3, dz4, dz5, dssp, u)
      for(i = ist; i <= iend; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(jst, jend, nz, tz2, tz3, dz1, tz1, dz2, dz3, dz4, dz5, dssp, u)
         for(j = jst; j <= jend; j++) {
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, u)
            for(k = 0; k <= nz - 1; k++) {
               flux[i][j][k][0] = u[i][j][k][3];
               u41 = u[i][j][k][3] / u[i][j][k][0];
               q = 0.50 * (u[i][j][k][1] * u[i][j][k][1] + u[i][j][k][2] * u[i][j][k][2] + u[i][j][k][3] * u[i][j][k][3]) / u[i][j][k][0];
               flux[i][j][k][1] = u[i][j][k][1] * u41;
               flux[i][j][k][2] = u[i][j][k][2] * u41;
               flux[i][j][k][3] = u[i][j][k][3] * u41 + 0.40e+00 * (u[i][j][k][4] - q);
               flux[i][j][k][4] = (1.40e+00 * u[i][j][k][4] - 0.40e+00 * q) * u41;
            }
            #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, tz2, flux)
            for(k = 1; k <= nz - 2; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = rsd[i][j][k][m] - tz2 * (flux[i][j][k + 1][m] - flux[i][j][k - 1][m]);
               }
            }
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, tz3, u)
            for(k = 1; k <= nz - 1; k++) {
               tmp = 1.0 / u[i][j][k][0];
               u21k = tmp * u[i][j][k][1];
               u31k = tmp * u[i][j][k][2];
               u41k = tmp * u[i][j][k][3];
               u51k = tmp * u[i][j][k][4];
               tmp = 1.0 / u[i][j][k - 1][0];
               u21km1 = tmp * u[i][j][k - 1][1];
               u31km1 = tmp * u[i][j][k - 1][2];
               u41km1 = tmp * u[i][j][k - 1][3];
               u51km1 = tmp * u[i][j][k - 1][4];
               flux[i][j][k][1] = tz3 * (u21k - u21km1);
               flux[i][j][k][2] = tz3 * (u31k - u31km1);
               flux[i][j][k][3] = (4.0 / 3.0) * tz3 * (u41k - u41km1);
               flux[i][j][k][4] = 0.50 * (1.0 - 1.40e+00 * 1.40e+00) * tz3 * ((((u21k) * (u21k)) + ((u31k) * (u31k)) + ((u41k) * (u41k))) - (((u21km1) * (u21km1)) + ((u31km1) * (u31km1)) + ((u41km1) * (u41km1)))) + (1.0 / 6.0) * tz3 * (((u41k) * (u41k)) - ((u41km1) * (u41km1))) + 1.40e+00 * 1.40e+00 * tz3 * (u51k - u51km1);
            }
            #pragma omp parallel for default(shared) private(k) firstprivate(nz, dz1, tz1, tz3, dz2, dz3, dz4, dz5, u, flux)
            for(k = 1; k <= nz - 2; k++) {
               rsd[i][j][k][0] = rsd[i][j][k][0] + dz1 * tz1 * (u[i][j][k - 1][0] - 2.0 * u[i][j][k][0] + u[i][j][k + 1][0]);
               rsd[i][j][k][1] = rsd[i][j][k][1] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][1] - flux[i][j][k][1]) + dz2 * tz1 * (u[i][j][k - 1][1] - 2.0 * u[i][j][k][1] + u[i][j][k + 1][1]);
               rsd[i][j][k][2] = rsd[i][j][k][2] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][2] - flux[i][j][k][2]) + dz3 * tz1 * (u[i][j][k - 1][2] - 2.0 * u[i][j][k][2] + u[i][j][k + 1][2]);
               rsd[i][j][k][3] = rsd[i][j][k][3] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][3] - flux[i][j][k][3]) + dz4 * tz1 * (u[i][j][k - 1][3] - 2.0 * u[i][j][k][3] + u[i][j][k + 1][3]);
               rsd[i][j][k][4] = rsd[i][j][k][4] + tz3 * 1.00e-01 * 1.00e+00 * (flux[i][j][k + 1][4] - flux[i][j][k][4]) + dz5 * tz1 * (u[i][j][k - 1][4] - 2.0 * u[i][j][k][4] + u[i][j][k + 1][4]);
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rsd[i][j][1][m] = rsd[i][j][1][m] - dssp * (+5.0 * u[i][j][1][m] - 4.0 * u[i][j][2][m] + u[i][j][3][m]);
               rsd[i][j][2][m] = rsd[i][j][2][m] - dssp * (-4.0 * u[i][j][1][m] + 6.0 * u[i][j][2][m] - 4.0 * u[i][j][3][m] + u[i][j][4][m]);
            }
            #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, dssp, u)
            for(k = 3; k <= nz - 4; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  rsd[i][j][k][m] = rsd[i][j][k][m] - dssp * (u[i][j][k - 2][m] - 4.0 * u[i][j][k - 1][m] + 6.0 * u[i][j][k][m] - 4.0 * u[i][j][k + 1][m] + u[i][j][k + 2][m]);
               }
            }
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(m = 0; m < 5; m++) {
               rsd[i][j][nz - 3][m] = rsd[i][j][nz - 3][m] - dssp * (u[i][j][nz - 5][m] - 4.0 * u[i][j][nz - 4][m] + 6.0 * u[i][j][nz - 3][m] - 4.0 * u[i][j][nz - 2][m]);
               rsd[i][j][nz - 2][m] = rsd[i][j][nz - 2][m] - dssp * (u[i][j][nz - 4][m] - 4.0 * u[i][j][nz - 3][m] + 5.0 * u[i][j][nz - 2][m]);
            }
         }
      }
   }
}

static void setbv() {
   {
      int i, j, k;
      int iglob, jglob;
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      exact#1735{exact(iglob, jglob, 0, &u[i][j][0][0])}
      exact#1736{exact(iglob, jglob, nz - 1, &u[i][j][nz - 1][0])}
      ****************************************/
      for(i = 0; i < nx; i++) {
         iglob = i;
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         exact#1735{exact(iglob, jglob, 0, &u[i][j][0][0])}
         exact#1736{exact(iglob, jglob, nz - 1, &u[i][j][nz - 1][0])}
         ****************************************/
         for(j = 0; j < ny; j++) {
            jglob = j;
            exact(iglob, jglob, 0, &u[i][j][0][0]);
            exact(iglob, jglob, nz - 1, &u[i][j][nz - 1][0]);
         }
      }
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      exact#1746{exact(iglob, 0, k, &u[i][0][k][0])}
      ****************************************/
      for(i = 0; i < nx; i++) {
         iglob = i;
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         exact#1746{exact(iglob, 0, k, &u[i][0][k][0])}
         ****************************************/
         for(k = 0; k < nz; k++) {
            exact(iglob, 0, k, &u[i][0][k][0]);
         }
      }
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      exact#1756{exact(iglob, ny0 - 1, k, &u[i][ny - 1][k][0])}
      ****************************************/
      for(i = 0; i < nx; i++) {
         iglob = i;
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         exact#1756{exact(iglob, ny0 - 1, k, &u[i][ny - 1][k][0])}
         ****************************************/
         for(k = 0; k < nz; k++) {
            exact(iglob, ny0 - 1, k, &u[i][ny - 1][k][0]);
         }
      }
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      exact#1766{exact(0, jglob, k, &u[0][j][k][0])}
      ****************************************/
      for(j = 0; j < ny; j++) {
         jglob = j;
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         exact#1766{exact(0, jglob, k, &u[0][j][k][0])}
         ****************************************/
         for(k = 0; k < nz; k++) {
            exact(0, jglob, k, &u[0][j][k][0]);
         }
      }
      /*************** Clava msgError **************
      		Variables Access as passed arguments Can not be traced inside of function calls : 
      exact#1776{exact(nx0 - 1, jglob, k, &u[nx - 1][j][k][0])}
      ****************************************/
      for(j = 0; j < ny; j++) {
         jglob = j;
         /*************** Clava msgError **************
         		Variables Access as passed arguments Can not be traced inside of function calls : 
         exact#1776{exact(nx0 - 1, jglob, k, &u[nx - 1][j][k][0])}
         ****************************************/
         for(k = 0; k < nz; k++) {
            exact(nx0 - 1, jglob, k, &u[nx - 1][j][k][0]);
         }
      }
   }
}

static void setcoeff() {
   dxi = 1.0 / (nx0 - 1);
   deta = 1.0 / (ny0 - 1);
   dzeta = 1.0 / (nz0 - 1);
   tx1 = 1.0 / (dxi * dxi);
   tx2 = 1.0 / (2.0 * dxi);
   tx3 = 1.0 / dxi;
   ty1 = 1.0 / (deta * deta);
   ty2 = 1.0 / (2.0 * deta);
   ty3 = 1.0 / deta;
   tz1 = 1.0 / (dzeta * dzeta);
   tz2 = 1.0 / (2.0 * dzeta);
   tz3 = 1.0 / dzeta;
   ii1 = 1;
   ii2 = nx0 - 2;
   ji1 = 1;
   ji2 = ny0 - 3;
   ki1 = 2;
   ki2 = nz0 - 2;
   dx1 = 0.75;
   dx2 = dx1;
   dx3 = dx1;
   dx4 = dx1;
   dx5 = dx1;
   dy1 = 0.75;
   dy2 = dy1;
   dy3 = dy1;
   dy4 = dy1;
   dy5 = dy1;
   dz1 = 1.00;
   dz2 = dz1;
   dz3 = dz1;
   dz4 = dz1;
   dz5 = dz1;
   dssp = ((((dx1) > ((((dy1) > (dz1)) ? (dy1) : (dz1)))) ? (dx1) : ((((dy1) > (dz1)) ? (dy1) : (dz1))))) / 4.0;
   ce[0][0] = 2.0;
   ce[0][1] = 0.0;
   ce[0][2] = 0.0;
   ce[0][3] = 4.0;
   ce[0][4] = 5.0;
   ce[0][5] = 3.0;
   ce[0][6] = 5.0e-01;
   ce[0][7] = 2.0e-02;
   ce[0][8] = 1.0e-02;
   ce[0][9] = 3.0e-02;
   ce[0][10] = 5.0e-01;
   ce[0][11] = 4.0e-01;
   ce[0][12] = 3.0e-01;
   ce[1][0] = 1.0;
   ce[1][1] = 0.0;
   ce[1][2] = 0.0;
   ce[1][3] = 0.0;
   ce[1][4] = 1.0;
   ce[1][5] = 2.0;
   ce[1][6] = 3.0;
   ce[1][7] = 1.0e-02;
   ce[1][8] = 3.0e-02;
   ce[1][9] = 2.0e-02;
   ce[1][10] = 4.0e-01;
   ce[1][11] = 3.0e-01;
   ce[1][12] = 5.0e-01;
   ce[2][0] = 2.0;
   ce[2][1] = 2.0;
   ce[2][2] = 0.0;
   ce[2][3] = 0.0;
   ce[2][4] = 0.0;
   ce[2][5] = 2.0;
   ce[2][6] = 3.0;
   ce[2][7] = 4.0e-02;
   ce[2][8] = 3.0e-02;
   ce[2][9] = 5.0e-02;
   ce[2][10] = 3.0e-01;
   ce[2][11] = 5.0e-01;
   ce[2][12] = 4.0e-01;
   ce[3][0] = 2.0;
   ce[3][1] = 2.0;
   ce[3][2] = 0.0;
   ce[3][3] = 0.0;
   ce[3][4] = 0.0;
   ce[3][5] = 2.0;
   ce[3][6] = 3.0;
   ce[3][7] = 3.0e-02;
   ce[3][8] = 5.0e-02;
   ce[3][9] = 4.0e-02;
   ce[3][10] = 2.0e-01;
   ce[3][11] = 1.0e-01;
   ce[3][12] = 3.0e-01;
   ce[4][0] = 5.0;
   ce[4][1] = 4.0;
   ce[4][2] = 3.0;
   ce[4][3] = 2.0;
   ce[4][4] = 1.0e-01;
   ce[4][5] = 4.0e-01;
   ce[4][6] = 3.0e-01;
   ce[4][7] = 5.0e-02;
   ce[4][8] = 4.0e-02;
   ce[4][9] = 3.0e-02;
   ce[4][10] = 1.0e-01;
   ce[4][11] = 3.0e-01;
   ce[4][12] = 2.0e-01;
}

static void setiv() {
   {
      int i, j, k, m;
      int iglob, jglob;
      double xi, eta, zeta;
      double pxi, peta, pzeta;
      double ue_1jk[5];
      double ue_nx0jk[5];
      double ue_i1k[5];
      double ue_iny0k[5];
      double ue_ij1[5];
      double ue_ijnz[5];
      #pragma omp parallel for default(shared) private(j, k, i, m) firstprivate(ny, nz, ny0, nx, nx0, ce, ue_1jk, ue_nx0jk, ue_i1k, ue_iny0k, ue_ij1, ue_ijnz)
      for(j = 0; j < ny; j++) {
         jglob = j;
         #pragma omp parallel for default(shared) private(k, i, m) firstprivate(nz, ny0, nx, nx0, ce, ue_1jk, ue_nx0jk, ue_i1k, ue_iny0k, ue_ij1, ue_ijnz)
         for(k = 1; k < nz - 1; k++) {
            zeta = ((double) k) / (nz - 1);
            if(jglob != 0 && jglob != ny0 - 1) {
               eta = ((double) (jglob)) / (ny0 - 1);
               #pragma omp parallel for default(shared) private(i, m) firstprivate(nx, nx0, ny0, nz, ce, ue_1jk, ue_nx0jk, ue_i1k, ue_iny0k, ue_ij1, ue_ijnz)
               for(i = 0; i < nx; i++) {
                  iglob = i;
                  if(iglob != 0 && iglob != nx0 - 1) {
                     xi = ((double) (iglob)) / (nx0 - 1);
                     exact(0, jglob, k, ue_1jk);
                     exact(nx0 - 1, jglob, k, ue_nx0jk);
                     exact(iglob, 0, k, ue_i1k);
                     exact(iglob, ny0 - 1, k, ue_iny0k);
                     exact(iglob, jglob, 0, ue_ij1);
                     exact(iglob, jglob, nz - 1, ue_ijnz);
                     /*************** Clava msgError **************
                     		 Loop Iteration number is too low
                     ****************************************/
                     for(m = 0; m < 5; m++) {
                        pxi = (1.0 - xi) * ue_1jk[m] + xi * ue_nx0jk[m];
                        peta = (1.0 - eta) * ue_i1k[m] + eta * ue_iny0k[m];
                        pzeta = (1.0 - zeta) * ue_ij1[m] + zeta * ue_ijnz[m];
                        u[i][j][k][m] = pxi + peta + pzeta - pxi * peta - peta * pzeta - pzeta * pxi + pxi * peta * pzeta;
                     }
                  }
               }
            }
         }
      }
   }
}

static void ssor() {
   int i, j, k, m;
   int istep;
   double tmp;
   double delunm[5];
   double tv[64][64][5];
   tmp = 1.0 / (omega * (2.0 - omega));
   {
      #pragma omp parallel for default(shared) private(i, j, k, m)
      for(i = 0; i < 64; i++) {
         #pragma omp parallel for default(shared) private(j, k, m) firstprivate(i)
         for(j = 0; j < 64; j++) {
            /*************** Clava msgError **************
            		 Loop Iteration number is too low
            ****************************************/
            for(k = 0; k < 5; k++) {
               /*************** Clava msgError **************
               		 Loop Iteration number is too low
               ****************************************/
               for(m = 0; m < 5; m++) {
                  a[i][j][k][m] = 0.0;
                  b[i][j][k][m] = 0.0;
                  c[i][j][k][m] = 0.0;
                  d[i][j][k][m] = 0.0;
               }
            }
         }
      }
   }
   rhs();
   l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
   timer_clear(1);
   timer_start(1);
   /*************** Clava msgError **************
   		Loop contains Invalid Statement -> exit#2895
   ****************************************/
   for(istep = 1; istep <= itmax; istep++) {
      if(istep % 20 == 0 || istep == itmax || istep == 1) {
         printf(" Time step %4d\n", istep);
      }
      {
         #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(ist, iend, jst, jend, nz, dt)
         for(i = ist; i <= iend; i++) {
            #pragma omp parallel for default(shared) private(j, k, m) firstprivate(jst, jend, nz, dt, i)
            for(j = jst; j <= jend; j++) {
               #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, dt, i, j)
               for(k = 1; k <= nz - 2; k++) {
                  /*************** Clava msgError **************
                  		 Loop Iteration number is too low
                  ****************************************/
                  for(m = 0; m < 5; m++) {
                     rsd[i][j][k][m] = dt * rsd[i][j][k][m];
                  }
               }
            }
         }
         /*************** Clava msgError **************
         		unsolved dependency for arrayAccess rsd	 use : RW
         ****************************************/
         for(k = 1; k <= nz - 2; k++) {
            jacld(k);
            blts(nx, ny, nz, k, omega, rsd, a, b, c, d, ist, iend, jst, jend, nx0, ny0);
         }
         /*************** Clava msgError **************
         		unsolved dependency for arrayAccess rsd	 use : RW
         ****************************************/
         for(k = nz - 2; k >= 1; k--) {
            jacu(k);
            buts(nx, ny, nz, k, omega, rsd, tv, d, a, b, c, ist, iend, jst, jend, nx0, ny0);
         }
         #pragma omp parallel for default(shared) private(i, j, k, m) firstprivate(ist, iend, jst, jend, nz, tmp, rsd)
         for(i = ist; i <= iend; i++) {
            #pragma omp parallel for default(shared) private(j, k, m) firstprivate(jst, jend, nz, tmp, i, rsd)
            for(j = jst; j <= jend; j++) {
               #pragma omp parallel for default(shared) private(k, m) firstprivate(nz, tmp, i, j, rsd)
               for(k = 1; k <= nz - 2; k++) {
                  /*************** Clava msgError **************
                  		 Loop Iteration number is too low
                  ****************************************/
                  for(m = 0; m < 5; m++) {
                     u[i][j][k][m] = u[i][j][k][m] + tmp * rsd[i][j][k][m];
                  }
               }
            }
         }
      }
      if(istep % inorm == 0) {
         l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsd, delunm);
      }
      rhs();
      if((istep % inorm == 0) || (istep == itmax)) {
         l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
      }
      if((rsdnm[0] < tolrsd[0]) && (rsdnm[1] < tolrsd[1]) && (rsdnm[2] < tolrsd[2]) && (rsdnm[3] < tolrsd[3]) && (rsdnm[4] < tolrsd[4])) {
         exit(1);
      }
   }
   timer_stop(1);
   maxtime = timer_read(1);
}

static void verify(double xcr[5], double xce[5], double xci, char *class, boolean *verified) {
   double xcrref[5];
   double xceref[5];
   double xciref;
   double xcrdif[5];
   double xcedif[5];
   double xcidif;
   double epsilon;
   double dtref;
   int m;
   epsilon = 1.0e-08;
   *class = 'U';
   *verified = 1;
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(m = 0; m < 5; m++) {
      xcrref[m] = 1.0;
      xceref[m] = 1.0;
   }
   xciref = 1.0;
   if(nx0 == 12 && ny0 == 12 && nz0 == 12 && itmax == 50) {
      *class = 'S';
      dtref = 5.0e-1;
      xcrref[0] = 1.6196343210976702e-02;
      xcrref[1] = 2.1976745164821318e-03;
      xcrref[2] = 1.5179927653399185e-03;
      xcrref[3] = 1.5029584435994323e-03;
      xcrref[4] = 3.4264073155896461e-02;
      xceref[0] = 6.4223319957960924e-04;
      xceref[1] = 8.4144342047347926e-05;
      xceref[2] = 5.8588269616485186e-05;
      xceref[3] = 5.8474222595157350e-05;
      xceref[4] = 1.3103347914111294e-03;
      xciref = 7.8418928865937083;
   }
   else if(nx0 == 33 && ny0 == 33 && nz0 == 33 && itmax == 300) {
         *class = 'W';
         dtref = 1.5e-3;
         xcrref[0] = 0.1236511638192e+02;
         xcrref[1] = 0.1317228477799e+01;
         xcrref[2] = 0.2550120713095e+01;
         xcrref[3] = 0.2326187750252e+01;
         xcrref[4] = 0.2826799444189e+02;
         xceref[0] = 0.4867877144216;
         xceref[1] = 0.5064652880982e-01;
         xceref[2] = 0.9281818101960e-01;
         xceref[3] = 0.8570126542733e-01;
         xceref[4] = 0.1084277417792e+01;
         xciref = 0.1161399311023e+02;
      }
      else if(nx0 == 64 && ny0 == 64 && nz0 == 64 && itmax == 250) {
            *class = 'A';
            dtref = 2.0e+0;
            xcrref[0] = 7.7902107606689367e+02;
            xcrref[1] = 6.3402765259692870e+01;
            xcrref[2] = 1.9499249727292479e+02;
            xcrref[3] = 1.7845301160418537e+02;
            xcrref[4] = 1.8384760349464247e+03;
            xceref[0] = 2.9964085685471943e+01;
            xceref[1] = 2.8194576365003349;
            xceref[2] = 7.3473412698774742;
            xceref[3] = 6.7139225687777051;
            xceref[4] = 7.0715315688392578e+01;
            xciref = 2.6030925604886277e+01;
         }
         else if(nx0 == 102 && ny0 == 102 && nz0 == 102 && itmax == 250) {
               *class = 'B';
               dtref = 2.0e+0;
               xcrref[0] = 3.5532672969982736e+03;
               xcrref[1] = 2.6214750795310692e+02;
               xcrref[2] = 8.8333721850952190e+02;
               xcrref[3] = 7.7812774739425265e+02;
               xcrref[4] = 7.3087969592545314e+03;
               xceref[0] = 1.1401176380212709e+02;
               xceref[1] = 8.1098963655421574;
               xceref[2] = 2.8480597317698308e+01;
               xceref[3] = 2.5905394567832939e+01;
               xceref[4] = 2.6054907504857413e+02;
               xciref = 4.7887162703308227e+01;
            }
            else if(nx0 == 162 && ny0 == 162 && nz0 == 162 && itmax == 250) {
                  *class = 'C';
                  dtref = 2.0e+0;
                  xcrref[0] = 1.03766980323537846e+04;
                  xcrref[1] = 8.92212458801008552e+02;
                  xcrref[2] = 2.56238814582660871e+03;
                  xcrref[3] = 2.19194343857831427e+03;
                  xcrref[4] = 1.78078057261061185e+04;
                  xceref[0] = 2.15986399716949279e+02;
                  xceref[1] = 1.55789559239863600e+01;
                  xceref[2] = 5.41318863077207766e+01;
                  xceref[3] = 4.82262643154045421e+01;
                  xceref[4] = 4.55902910043250358e+02;
                  xciref = 6.66404553572181300e+01;
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
   xcidif = fabs((xci - xciref) / xciref);
   if(*class != 'U') {
      printf("\n Verification being performed for class %1c\n", *class);
      printf(" Accuracy setting for epsilon = %20.13e\n", epsilon);
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
         printf("          %2d  %20.13e\n", m, xcr[m]);
      }
      else if(xcrdif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
         }
         else {
            printf("          %2d  %20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
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
         printf("          %2d  %20.13e\n", m, xce[m]);
      }
      else if(xcedif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
         }
         else {
            printf("          %2d  %20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
         }
   }
   if(*class != 'U') {
      printf(" Comparison of surface integral\n");
   }
   else {
      printf(" Surface integral\n");
   }
   if(*class == 'U') {
      printf("              %20.13e\n", xci);
   }
   else if(xcidif > epsilon) {
         *verified = 0;
         printf(" FAILURE:     %20.13e%20.13e%20.13e\n", xci, xciref, xcidif);
      }
      else {
         printf("              %20.13e%20.13e%20.13e\n", xci, xciref, xcidif);
      }
   if(*class == 'U') {
      printf(" No reference values provided\n");
      printf(" No verification performed\n");
   }
   else if(*verified) {
         printf(" Verification Successful\n");
      }
      else {
         printf(" Verification failed\n");
      }
}
