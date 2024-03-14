/*
 * file for is.c
 */





void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand);

void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, int nthreads, double t, double mops, char *optype, int passed_verification, char *npbversion, char *compiletime, char *cc, char *clink, char *c_lib, char *c_inc, char *cflags, char *clinkflags, char *rand);

void wtime_(double *t);
double elapsed_time();

double start[64];
double elapsed[64];

void timer_clear(int n);


void timer_start(int n);


void timer_stop(int n);


double timer_read(int n);


#include <sys/time.h>

void wtime_(double *t);


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

#include "npbparams.h"
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
typedef int INT_TYPE;





INT_TYPE *key_buff_ptr_global;


int passed_verification;
INT_TYPE key_array[1<<23], key_buff1[1<<23], key_buff2[1<<23], partial_verify_vals[5];


























INT_TYPE test_index_array[5], test_rank_array[5], S_test_index_array[5] = {48427, 17148, 23627, 62548, 4431}, S_test_rank_array[5] = {0, 18, 346, 64917, 65463}, W_test_index_array[5] = {357773, 934767, 875723, 898999, 404505}, W_test_rank_array[5] = {1249, 11698, 1039987, 1043896, 1048018}, A_test_index_array[5] = {2112377, 662041, 5336171, 3642833, 4250760}, A_test_rank_array[5] = {104, 17523, 123928, 8288932, 8388264}, B_test_index_array[5] = {41869, 812306, 5102857, 18232239, 26860214}, B_test_rank_array[5] = {33422937, 10244, 59149, 33135281, 99}, C_test_index_array[5] = {44172927, 72999161, 74326391, 129606274, 21736814}, C_test_rank_array[5] = {61147, 882988, 266290, 133997595, 133525895};






double randlc(double *X, double *A);

void full_verify();
double randlc(double *X, double *A);

void create_seq(double seed, double a);

void full_verify();

void rank(int iteration);







 main(int argc, char **argv);
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
void wtime_(double *t)
{
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, (void *) 0);
   if (sec<0)
      sec = tv.tv_sec;
   *t = tv.tv_sec-sec+1.0e-6*tv.tv_usec;
}
//static __uint64_t __bswap_64(__uint64_t __bsx)
//{
  // return (__bsx&0xff00000000000000ull)>>56|(__bsx&0x00ff000000000000ull)>>40|(__bsx&0x0000ff0000000000ull)>>24|(__bsx&0x000000ff00000000ull)>>8|(__bsx&0x00000000ff000000ull)<<8|(__bsx&0x0000000000ff0000ull)<<24|(__bsx&0x000000000000ff00ull)<<40|(__bsx&0x00000000000000ffull)<<56;
//}
double randlc(double *X, double *A)
{
   static int KS = 0;
   static double R23, R46, T23, T46;
   double T1, T2, T3, T4;
   double A1;
   double A2;
   double X1;
   double X2;
   double Z;
   int i, j;

   if (KS==0) {
      R23 = 1.0;
      R46 = 1.0;
      T23 = 1.0;
      T46 = 1.0;

#pragma omp parallel for reduction(*:T23) reduction(*:R23)
      for(i = 1; i <= 23; i += 1) {
         R23 = 0.50*R23;
         T23 = 2.0*T23;
      }
#pragma omp parallel for reduction(*:T46) reduction(*:R46)
      for(i = 1; i <= 46; i += 1) {
         R46 = 0.50*R46;
         T46 = 2.0*T46;
      }
      KS = 1;
   }
   
   
   
   T1 = R23**A;
   j = T1;
   A1 = j;
   A2 = *A-T23*A1;
   
   
   
   T1 = R23**X;
   j = T1;
   X1 = j;
   X2 = *X-T23*X1;
   T1 = A1*X2+A2*X1;

   j = R23*T1;
   T2 = j;
   Z = T1-T23*T2;
   T3 = T23*Z+A2*X2;
   j = R46*T3;
   T4 = j;
   *X = T3-T46*T4;
   return R46**X;
}
void create_seq(double seed, double a)
{
   double x;
   int i, j, k;
   k = (1<<19)/4;

   for(i = 0; i <= (1<<23)-1; i += 1) {
      x = randlc(&seed, &a);
      x += randlc(&seed, &a);
      x += randlc(&seed, &a);
      x += randlc(&seed, &a);

      key_array[i] = k*x;
   }
}
void full_verify()
{
   INT_TYPE i, j;
   INT_TYPE k;
   INT_TYPE m, unique_keys;
   
   
   
   
   for(i = 0; i <= (1<<23)-1; i += 1)
      key_array[--key_buff_ptr_global[key_buff2[i]]] = key_buff2[i];
   
   
   
   
   j = 0;
#pragma omp parallel for reduction(+:j)
   for(i = 1; i <= (1<<23)-1; i += 1)
      if (key_array[i-1]>key_array[i])
         j++;
   
   
   if (j!=0)

      printf("Full_verify: number of keys out of sort: %d\n", j);
   else
      passed_verification++;
}
void rank(int iteration)
{

   INT_TYPE i, j, k;
   INT_TYPE l, m;

   INT_TYPE shift = 19-10;
   INT_TYPE key;
   INT_TYPE min_key_val, max_key_val;
   INT_TYPE prv_buff1[1<<19];
   key_array[iteration] = iteration;
   key_array[iteration+10] = (1<<19)-iteration;
   
   
   
#pragma omp parallel for
   for(i = 0; i <= 4; i += 1)
      partial_verify_vals[i] = key_array[test_index_array[i]];

#pragma omp parallel for
   for(i = 0; i <= (1<<19)-1; i += 1)
      prv_buff1[i] = 0;
   
   
   for(i = 0; i <= (1<<23)-1; i += 1) {
      key_buff2[i] = key_array[i];
      
      
      
      
      
      prv_buff1[key_buff2[i]]++;
   }

   for(i = 0; i <= (1<<19)-1-1; i += 1)
      prv_buff1[i+1] += prv_buff1[i];
   
   
#pragma omp parallel for
   for(i = 0; i <= (1<<19)-1; i += 1) {
      key_buff1[i] = 0;
      key_buff1[i] += prv_buff1[i];
   }
   
   
   
   
   for(i = 0; i <= 4; i += 1) {
      k = partial_verify_vals[i];
      if (0<=k&&k<=(1<<23)-1) {
         if ('A'=='S') {
_switch_8_case_S:            ;
            if (i<=2)
               if (key_buff1[k-1]!=test_rank_array[i]+iteration)
                  
                  
                  printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
               else
                  passed_verification++;
            else if (key_buff1[k-1]!=test_rank_array[i]-iteration)
               
               
               printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
            else
               passed_verification++;
         }
         else if ('A'=='W') {
_switch_8_case_W:            ;
            if (i<2)
               if (key_buff1[k-1]!=test_rank_array[i]+(iteration-2))
                  
                  
                  printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
               else
                  passed_verification++;
            else if (key_buff1[k-1]!=test_rank_array[i]-iteration)
               
               
               printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
            else
               passed_verification++;
         }
         else if ('A'=='A') {
_switch_8_case_A:            ;
            if (i<=2)
               if (key_buff1[k-1]!=test_rank_array[i]+(iteration-1))
                  
                  
                  printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
               else
                  passed_verification++;
            else if (key_buff1[k-1]!=test_rank_array[i]-(iteration-1))
               
               
               printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
            else
               passed_verification++;
         }
         else if ('A'=='B') {
_switch_8_case_B:            ;
            if (i==1||i==2||i==4)
               if (key_buff1[k-1]!=test_rank_array[i]+iteration)
                  
                  
                  printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
               else
                  passed_verification++;
            else if (key_buff1[k-1]!=test_rank_array[i]-iteration)
               
               
               printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
            else
               passed_verification++;
         }
         else if ('A'=='C') {
_switch_8_case_C:            ;
            if (i<=2)
               if (key_buff1[k-1]!=test_rank_array[i]+iteration)
                  
                  
                  printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
               else
                  passed_verification++;
            else if (key_buff1[k-1]!=test_rank_array[i]-iteration)
               
               
               printf("Failed partial verification: ""iteration %d, test key %d\n", iteration, i);
            else
               passed_verification++;
         }
_break_8:         ;
      }
   }
   
   
   
   
   
   
   if (iteration==10)
      key_buff_ptr_global = key_buff1;
}
int main(int argc, char **argv)
{
   double t1, t2;
   
   int i, iteration, itemp;
   int nthreads = 1;
   double timecounter, maxtime;
   
   
   
   t1 = omp_get_wtime();
#pragma omp parallel for
   for(i = 0; i <= 4; i += 1) {
      if ('A'=='S') {
_switch_2_case_S:         ;
         test_index_array[i] = S_test_index_array[i];
         test_rank_array[i] = S_test_rank_array[i];
      }
      else if ('A'=='A') {
_switch_2_case_A:         ;
         test_index_array[i] = A_test_index_array[i];
         test_rank_array[i] = A_test_rank_array[i];
      }
      else if ('A'=='W') {
_switch_2_case_W:         ;
         test_index_array[i] = W_test_index_array[i];
         test_rank_array[i] = W_test_rank_array[i];
      }
      else if ('A'=='B') {
_switch_2_case_B:         ;
         test_index_array[i] = B_test_index_array[i];
         test_rank_array[i] = B_test_rank_array[i];
      }
      else if ('A'=='C') {
_switch_2_case_C:         ;
         test_index_array[i] = C_test_index_array[i];
         test_rank_array[i] = C_test_rank_array[i];
      }
_break_2:      ;
   }
   
   
   
   
   
   printf("\n\n NAS Parallel Benchmarks 2.3 OpenMP C version"" - IS Benchmark\n\n");
   printf(" Size:  %d  (class %c)\n", 1<<23, 'A');
   printf(" Iterations:   %d\n", 10);
   
   
   timer_clear(0);
   
   
   
   create_seq(314159265.00, 1220703125.00);
   
   
   
   rank(1);
   
   
   passed_verification = 0;

   if ('A'!='S')
      printf("\n   iteration\n");
   
   
   timer_start(0);
   
   
   
   
   for(iteration = 1; iteration <= 10; iteration += 1) {
      if ('A'!='S')
         printf("        %d\n", iteration);

      rank(iteration);
   }
   
   t2 = omp_get_wtime();
   timer_stop(0);
   timecounter = timer_read(0);
   
   
   
   full_verify();
   
   
   
   
   if (passed_verification!=5*10+1)
      passed_verification = 0;
   
   
   
   
   
   
   
   
   
   
   
   c_print_results("IS", 'A', 1<<23, 0, 0, 10, nthreads, timecounter, (double) (10*(1<<23))/timecounter/1000000., "keys ranked", passed_verification, "3.0 structured", "17 Jan 2020", "(none)", "(none)", "-lm", "(none)", "(none)", "(none)", "randlc");

   printf("\n%lf\n", t2 - t1);
}
