﻿#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

typedef int boolean;
typedef struct { double real; double imag; } dcomplex;

#define TRUE	1
#define FALSE	0

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define	pow2(a) ((a)*(a))

#define get_real(c) c.real
#define get_imag(c) c.imag
#define cadd(c,a,b) (c.real = a.real + b.real, c.imag = a.imag + b.imag)
#define csub(c,a,b) (c.real = a.real - b.real, c.imag = a.imag - b.imag)
#define cmul(c,a,b) (c.real = a.real * b.real - a.imag * b.imag, \
                     c.imag = a.real * b.imag + a.imag * b.real)
#define crmul(c,a,b) (c.real = a.real * b, c.imag = a.imag * b)

#if defined(USE_POW)
#define r23 pow(0.5, 23.0)
#define r46 (r23*r23)
#define t23 pow(2.0, 23.0)
#define t46 (t23*t23)
#else
#define r23 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define r46 (r23*r23)
#define t23 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)
#define t46 (t23*t23)
#endif

#define wtime wtime_

void wtime(double *t)
{
  static int sec = -1;
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  if (sec < 0) sec = tv.tv_sec;
  *t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}
void wtime( double * );
double elapsed_time( void )
{
    double t;

    wtime( &t );
    return( t );
}

double start[64], elapsed[64];
void timer_clear( int n )
{
    elapsed[n] = 0.0;
}
void timer_start( int n )
{
    start[n] = elapsed_time();
}
void timer_stop( int n )
{
    double t, now;

    now = elapsed_time();
    t = now - start[n];
    elapsed[n] += t;

}
double timer_read( int n )
{
    return( elapsed[n] );
}

double randlc (double *x, double a) {
    double t1,t2,t3,t4,a1,a2,x1,x2,z;
    t1 = r23 * a;
    a1 = (int)t1;
    a2 = a - t23 * a1;
    t1 = r23 * (*x);
    x1 = (int)t1;
    x2 = (*x) - t23 * x1;
    t1 = a1 * x2 + a2 * x1;
    t2 = (int)(r23 * t1);
    z = t1 - t23 * t2;
    t3 = t23 * z + a2 * x2;
    t4 = (int)(r46 * t3);
    (*x) = t3 - t46 * t4;
    return (r46 * (*x));
}
void vranlc (int n, double *x_seed, double a, double y[]) {
    int i;
    double x,t1,t2,t3,t4,a1,a2,x1,x2,z;
    t1 = r23 * a;
    a1 = (int)t1;
    a2 = a - t23 * a1;
    x = *x_seed;
    for (i = 1; i <= n; i++) {
        t1 = r23 * x;
        x1 = (int)t1;
        x2 = x - t23 * x1;
        t1 = a1 * x2 + a2 * x1;
        t2 = (int)(r23 * t1);
        z = t1 - t23 * t2;
        t3 = t23 * z + a2 * x2;
        t4 = (int)(r46 * t3);
        x = t3 - t46 * t4;
        y[i] = r46 * x;
    }
    *x_seed = x;
}
void c_print_results( char   *name,
                      char   class,
                      int    n1, 
                      int    n2,
                      int    n3,
                      int    niter,
		      int    nthreads,
                      double t,
                      double mops,
		      char   *optype,
                      int    passed_verification,
                      char   *npbversion,
                      char   *compiletime,
                      char   *cc,
                      char   *clink,
                      char   *c_lib,
                      char   *c_inc,
                      char   *cflags,
                      char   *clinkflags,
		      char   *rand)
{
    char *evalue="1000";

    printf( "\n\n %s Benchmark Completed\n", name ); 

    printf( " Class           =                        %c\n", class );

    if( n2 == 0 && n3 == 0 )
        printf( " Size            =             %12d\n", n1 );   /* as in IS */
    else
        printf( " Size            =              %3dx%3dx%3d\n", n1,n2,n3 );

    printf( " Iterations      =             %12d\n", niter );
    
    printf( " Threads         =             %12d\n", nthreads );
 
    printf( " Time in seconds =             %12.2f\n", t );

    printf( " Mop/s total     =             %12.2f\n", mops );

    printf( " Operation type  = %24s\n", optype);

    if( passed_verification )
        printf( " Verification    =               SUCCESSFUL\n" );
    else
        printf( " Verification    =             UNSUCCESSFUL\n" );

    printf( " Version         =           %12s\n", npbversion );

    printf( " Compile date    =             %12s\n", compiletime );

    printf( "\n Compile options:\n" );

    printf( "    CC           = %s\n", cc );

    printf( "    CLINK        = %s\n", clink );

    printf( "    C_LIB        = %s\n", c_lib );

    printf( "    C_INC        = %s\n", c_inc );

    printf( "    CFLAGS       = %s\n", cflags );

    printf( "    CLINKFLAGS   = %s\n", clinkflags );

    printf( "    RAND         = %s\n", rand );
}
/* CLASS = A */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class of the NPB
c  in this directory. Do not modify it by hand.
*/
#define	NA	14000
#define	NONZER	11
#define	NITER	15
#define	SHIFT	20.0
#define	RCOND	1.0e-1
#define	CONVERTDOUBLE	FALSE
#define COMPILETIME "25 Nov 2023"
#define NPBVERSION "3.0 structured"
#define CS1 "gcc-12"
#define CS2 "gcc-12"
#define CS3 "-fopenmp -lm"
#define CS4 "-I../common -fopenmp"
#define CS5 "-lm"
#define CS6 "(none)"
#define CS7 "randdp"

#define	NZ	NA*(NONZER+1)*(NONZER+1)+NA*(NONZER+2)




static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;


static int colidx[NZ+1];	
static int rowstr[NA+1+1];	
static int iv[2*NA+1+1];	
static int arow[NZ+1];		
static int acol[NZ+1];		


static double v[NA+1+1];	
static double aelt[NZ+1];	
static double a[NZ+1];		
static double x[NA+2+1];	
static double z[NA+2+1];	
static double p[NA+2+1];	
static double q[NA+2+1];	
static double r[NA+2+1];	


static double amult;
static double tran;


static void conj_grad (int colidx[], int rowstr[], double x[], double z[],
		       double a[], double p[], double q[], double r[],
		       		       double *rnorm);
static void makea(int n, int nz, double a[], int colidx[], int rowstr[],
		  int nonzer, int firstrow, int lastrow, int firstcol,
		  int lastcol, double rcond, int arow[], int acol[],
		  double aelt[], double v[], int iv[], double shift );
static void sparse(double a[], int colidx[], int rowstr[], int n,
		   int arow[], int acol[], double aelt[],
		   int firstrow, int lastrow,
		   double x[], boolean mark[], int nzloc[], int nnza);
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[],
		   int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);



int main(int argc, char **argv) {
    int	i, j, k, it;
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
    lastrow  = NA;
    firstcol = 1;
    lastcol  = NA;

    if (NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10.0) {
	class = 'S';
	zeta_verify_value = 8.5971775078648;
    } else if (NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12.0) {
	class = 'W';
	zeta_verify_value = 10.362595087124;
    } else if (NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20.0) {
	class = 'A';
	zeta_verify_value = 17.130235054029;
    } else if (NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60.0) {
	class = 'B';
	zeta_verify_value = 22.712745482631;
    } else if (NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110.0) {
	class = 'C';
	zeta_verify_value = 28.973605592845;
    } else {
	class = 'U';
    }

    printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
	   " - CG Benchmark\n");
    printf(" Size: %10d\n", NA);
    printf(" Iterations: %5d\n", NITER);

    naa = NA;
    nzz = NZ;


    tran    = 314159265.0;
    amult   = 1220703125.0;
    zeta    = randlc( &tran, amult );


    makea(naa, nzz, a, colidx, rowstr, NONZER,
	  firstrow, lastrow, firstcol, lastcol, 
	  RCOND, arow, acol, aelt, v, iv, SHIFT);
    

{	
    for (j = 1; j <= lastrow - firstrow + 1; j++) {
	for (k = rowstr[j]; k < rowstr[j+1]; k++) {
            colidx[k] = colidx[k] - firstcol + 1;
	}
    }


    for (i = 1; i <= NA+1; i++) {
	x[i] = 1.0;
    }
      for (j = 1; j <= lastcol-firstcol+1; j++) {
         q[j] = 0.0;
         z[j] = 0.0;
         r[j] = 0.0;
         p[j] = 0.0;
      }
}    zeta  = 0.0;



    for (it = 1; it <= 1; it++) {


	conj_grad (colidx, rowstr, x, z, a, p, q, r, &rnorm);


	norm_temp11 = 0.0;
	norm_temp12 = 0.0;
	for (j = 1; j <= lastcol-firstcol+1; j++) {
            norm_temp11 = norm_temp11 + x[j]*z[j];
            norm_temp12 = norm_temp12 + z[j]*z[j];
	}
	norm_temp12 = 1.0 / sqrt( norm_temp12 );


	for (j = 1; j <= lastcol-firstcol+1; j++) {
            x[j] = norm_temp12*z[j];
	}
	
    } 


    for (i = 1; i <= NA+1; i++) {
         x[i] = 1.0;
    }  
    zeta  = 0.0;


    timer_clear( 1 );
    timer_start( 1 );



    for (it = 1; it <= NITER; it++) {


	conj_grad(colidx, rowstr, x, z, a, p, q, r, &rnorm);


	norm_temp11 = 0.0;
	norm_temp12 = 0.0;

	for (j = 1; j <= lastcol-firstcol+1; j++) {
            norm_temp11 = norm_temp11 + x[j]*z[j];
            norm_temp12 = norm_temp12 + z[j]*z[j];
	}

	norm_temp12 = 1.0 / sqrt( norm_temp12 );

	zeta = SHIFT + 1.0 / norm_temp11;

	if( it == 1 ) {
	  printf("   iteration           ||r||                 zeta\n");
	}
	printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);


	for (j = 1; j <= lastcol-firstcol+1; j++) {
            x[j] = norm_temp12*z[j];
	}
    } 
{
} 

    timer_stop( 1 );



    t = timer_read( 1 );

    printf(" Benchmark completed\n");

    epsilon = 1.0e-10;
        if (class != 'U') {
	if (fabs(zeta - zeta_verify_value) <= epsilon) {
            verified = TRUE;
	    printf(" VERIFICATION SUCCESSFUL\n");
	    printf(" Zeta is    %20.12e\n", zeta);
	    printf(" Error is   %20.12e\n", zeta - zeta_verify_value);
	} else {
            verified = FALSE;
	    printf(" VERIFICATION FAILED\n");
	    printf(" Zeta                %20.12e\n", zeta);
	    printf(" The correct zeta is %20.12e\n", zeta_verify_value);
	}
    } else {
	verified = FALSE;
	printf(" Problem size unknown\n");
	printf(" NO VERIFICATION PERFORMED\n");
    }

    if ( t != 0.0 ) {
	mflops = (2.0*NITER*NA)
	    * (3.0+(NONZER*(NONZER+1)) + 25.0*(5.0+(NONZER*(NONZER+1))) + 3.0 )
	    / t / 1000000.0;
    } else {
	mflops = 0.0;
    }

    c_print_results("CG", class, NA, 0, 0, NITER, nthreads, t, 
		    mflops, "          floating point", 
		    verified, NPBVERSION, COMPILETIME,
		    CS1, CS2, CS3, CS4, CS5, CS6, CS7);
}


static void conj_grad (
    int colidx[],	
    int rowstr[],	
    double x[],		
    double z[],		
    double a[],		
    double p[],		
    double q[],		
    double r[],		
        double *rnorm )

    

{
    static int callcount = 0;
    double d, sum, rho, rho0, alpha, beta;
    int i, j, k;
    int cgit, cgitmax = 25;

    rho = 0.0;
    

{
    for (j = 1; j <= naa+1; j++) {
	q[j] = 0.0;
	z[j] = 0.0;
	r[j] = x[j];
	p[j] = r[j];
	    }


    for (j = 1; j <= lastcol-firstcol+1; j++) {
	rho = rho + r[j]*r[j];
    }
}

    for (cgit = 1; cgit <= cgitmax; cgit++) {
      rho0 = rho;
      d = 0.0;
      rho = 0.0;
{
      


    
	for (j = 1; j <= lastrow-firstrow+1; j++) {
            sum = 0.0;
	    for (k = rowstr[j]; k < rowstr[j+1]; k++) {
		sum = sum + a[k]*p[colidx[k]];
	    }
                        q[j] = sum;
	}
	






	for (j = 1; j <= lastcol-firstcol+1; j++) {
            d = d + p[j]*q[j];
	}

	alpha = rho0 / d;


	


	for (j = 1; j <= lastcol-firstcol+1; j++) {
            z[j] = z[j] + alpha*p[j];
            r[j] = r[j] - alpha*q[j];
            


            rho = rho + r[j]*r[j];
	}


	beta = rho / rho0;


	for (j = 1; j <= lastcol-firstcol+1; j++) {
            p[j] = r[j] + beta*p[j];
	}
    callcount++;
    } 
    } 


    sum = 0.0;
    
{
    for (j = 1; j <= lastrow-firstrow+1; j++) {
	d = 0.0;
	for (k = rowstr[j]; k <= rowstr[j+1]-1; k++) {
            d = d + a[k]*z[colidx[k]];
	}
	r[j] = d;
    }


    for (j = 1; j <= lastcol-firstcol+1; j++) {
	d = x[j] - r[j];
	sum = sum + d*d;
    }
}     (*rnorm) = sqrt(sum);
}


static void makea(
    int n,
    int nz,
    double a[],		
    int colidx[],	
    int rowstr[],	
    int nonzer,
    int firstrow,
    int lastrow,
    int firstcol,
    int lastcol,
    double rcond,
    int arow[],		
    int acol[],		
    double aelt[],	
    double v[],		
    int iv[],		
    double shift )
{
    int i, nnza, iouter, ivelt, ivelt1, irow, nzv;



    double size, ratio, scale;
    int jcol;

    size = 1.0;
    ratio = pow(rcond, (1.0 / (double)n));
    nnza = 0;


    for (i = 1; i <= n; i++) {
	colidx[n+i] = 0;
    }
    for (iouter = 1; iouter <= n; iouter++) {
	nzv = nonzer;
	sprnvc(n, nzv, v, iv, &(colidx[0]), &(colidx[n]));
	vecset(n, v, iv, &nzv, iouter, 0.5);
	for (ivelt = 1; ivelt <= nzv; ivelt++) {
	    jcol = iv[ivelt];
	    if (jcol >= firstcol && jcol <= lastcol) {
		scale = size * v[ivelt];
		for (ivelt1 = 1; ivelt1 <= nzv; ivelt1++) {
	            irow = iv[ivelt1];
                    if (irow >= firstrow && irow <= lastrow) {
			nnza = nnza + 1;
			if (nnza > nz) {
			    printf("Space for matrix elements exceeded in"
				   " makea\n");
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


    for (i = firstrow; i <= lastrow; i++) {
	if (i >= firstcol && i <= lastcol) {
	    iouter = n + i;
	    nnza = nnza + 1;
	    if (nnza > nz) {
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


    sparse(a, colidx, rowstr, n, arow, acol, aelt,
	   firstrow, lastrow, v, &(iv[0]), &(iv[n]), nnza);
}


static void sparse(
    double a[],		
    int colidx[],	
    int rowstr[],	
    int n,
    int arow[],		
    int acol[],		
    double aelt[],	
    int firstrow,
    int lastrow,
    double x[],		
    boolean mark[],	
    int nzloc[],	
    int nnza)

{
    int nrows;
    int i, j, jajp1, nza, k, nzrow;
    double xi;


    nrows = lastrow - firstrow + 1;


    for (j = 1; j <= n; j++) {
	rowstr[j] = 0;
	mark[j] = FALSE;
    }
    rowstr[n+1] = 0;
    
    for (nza = 1; nza <= nnza; nza++) {
	j = (arow[nza] - firstrow + 1) + 1;
	rowstr[j] = rowstr[j] + 1;
    }

    rowstr[1] = 1;
    for (j = 2; j <= nrows+1; j++) {
	rowstr[j] = rowstr[j] + rowstr[j-1];
    }


    

      for(j = 0;j <= nrows-1;j++) {
         for(k = rowstr[j];k <= rowstr[j+1]-1;k++)
	       a[k] = 0.0;
      }

    for (nza = 1; nza <= nnza; nza++) {
	j = arow[nza] - firstrow + 1;
	k = rowstr[j];
	a[k] = aelt[nza];
	colidx[k] = acol[nza];
	rowstr[j] = rowstr[j] + 1;
    }


    for (j = nrows; j >= 1; j--) {
	rowstr[j+1] = rowstr[j];
    }
    rowstr[1] = 1;


    nza = 0;
    for (i = 1; i <= n; i++) {
	x[i] = 0.0;
	mark[i] = FALSE;
    }

    jajp1 = rowstr[1];
    for (j = 1; j <= nrows; j++) {
	nzrow = 0;
	

	for (k = jajp1; k < rowstr[j+1]; k++) {
            i = colidx[k];
            x[i] = x[i] + a[k];
            if ( mark[i] == FALSE && x[i] != 0.0) {
		mark[i] = TRUE;
		nzrow = nzrow + 1;
		nzloc[nzrow] = i;
	    }
	}


	for (k = 1; k <= nzrow; k++) {
            i = nzloc[k];
            mark[i] = FALSE;
            xi = x[i];
            x[i] = 0.0;
            if (xi != 0.0) {
		nza = nza + 1;
		a[nza] = xi;
		colidx[nza] = i;
	    }
	}
	jajp1 = rowstr[j+1];
	rowstr[j+1] = nza + rowstr[1];
    }
}


static void sprnvc(
    int n,
    int nz,
    double v[],		
    int iv[],		
    int nzloc[],	
    int mark[] ) 	
{
    int nn1;
    int nzrow, nzv, ii, i;
    double vecelt, vecloc;

    nzv = 0;
    nzrow = 0;
    nn1 = 1;
    do {
	nn1 = 2 * nn1;
    } while (nn1 < n);



    while (nzv < nz) {
	vecelt = randlc(&tran, amult);


	vecloc = randlc(&tran, amult);
	i = icnvrt(vecloc, nn1) + 1;
	if (i > n) continue;


	if (mark[i] == 0) {
	    mark[i] = 1;
	    nzrow = nzrow + 1;
	    nzloc[nzrow] = i;
	    nzv = nzv + 1;
	    v[nzv] = vecelt;
	    iv[nzv] = i;
	}
    }

    for (ii = 1; ii <= nzrow; ii++) {
	i = nzloc[ii];
	mark[i] = 0;
    }
}


static int icnvrt(double x, int ipwr2) {
    return ((int)(ipwr2 * x));
}


static void vecset(
    int n,
    double v[],	
    int iv[],	
    int *nzv,
    int i,
    double val)
{
    int k;
    boolean set;

    set = FALSE;
    for (k = 1; k <= *nzv; k++) {
	if (iv[k] == i) {
            v[k] = val;
            set  = TRUE;
	}
    }
    if (set == FALSE) {
	*nzv = *nzv + 1;
	v[*nzv] = val;
	iv[*nzv] = i;
    }
}
