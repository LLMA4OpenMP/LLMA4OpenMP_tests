#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
/*
This file is generated automatically by the setparams utility.
It sets the number of processors and the class of the NPB
in this directory. Do not modify it by hand.*/

typedef int INT_TYPE;
INT_TYPE *key_buff_ptr_global;
int passed_verification;
INT_TYPE key_array[8388608];
INT_TYPE key_buff1[8388608];
INT_TYPE key_buff2[8388608];
INT_TYPE partial_verify_vals[5];
INT_TYPE test_index_array[5];
INT_TYPE test_rank_array[5];
INT_TYPE S_test_index_array[5] = {48427, 17148, 23627, 62548, 4431};
INT_TYPE S_test_rank_array[5] = {0, 18, 346, 64917, 65463};
INT_TYPE W_test_index_array[5] = {357773, 934767, 875723, 898999, 404505};
INT_TYPE W_test_rank_array[5] = {1249, 11698, 1039987, 1043896, 1048018};
INT_TYPE A_test_index_array[5] = {2112377, 662041, 5336171, 3642833, 4250760};
INT_TYPE A_test_rank_array[5] = {104, 17523, 123928, 8288932, 8388264};
INT_TYPE B_test_index_array[5] = {41869, 812306, 5102857, 18232239, 26860214};
INT_TYPE B_test_rank_array[5] = {33422937, 10244, 59149, 33135281, 99};
INT_TYPE C_test_index_array[5] = {44172927, 72999161, 74326391, 129606274, 21736814};
INT_TYPE C_test_rank_array[5] = {61147, 882988, 266290, 133997595, 133525895};
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

double randlc(double *X, double *A);
void full_verify();
double randlc(double *X, double *A) {
   static int KS = 0;
   static double R23, R46, T23, T46;
   double T1, T2, T3, T4;
   double A1;
   double A2;
   double X1;
   double X2;
   double Z;
   int i, j;
   if(KS == 0) {
      R23 = 1.0;
      R46 = 1.0;
      T23 = 1.0;
      T46 = 1.0;
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(i = 1; i <= 23; i++) {
         R23 = 0.50 * R23;
         T23 = 2.0 * T23;
      }
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(i = 1; i <= 46; i++) {
         R46 = 0.50 * R46;
         T46 = 2.0 * T46;
      }
      KS = 1;
   }
   T1 = R23 * *A;
   j = T1;
   A1 = j;
   A2 = *A - T23 * A1;
   T1 = R23 * *X;
   j = T1;
   X1 = j;
   X2 = *X - T23 * X1;
   T1 = A1 * X2 + A2 * X1;
   j = R23 * T1;
   T2 = j;
   Z = T1 - T23 * T2;
   T3 = T23 * Z + A2 * X2;
   j = R46 * T3;
   T4 = j;
   *X = T3 - T46 * T4;
   
   return (R46 * *X);
}

void create_seq(double seed, double a) {
   double x;
   int i, j, k;
   k = (1 << 19) / 4;
   #pragma omp parallel for default(shared) private(i, x) firstprivate(seed, a, k)
   for(i = 0; i < (1 << 23); i++) {
      x = randlc(&seed, &a);
      x += randlc(&seed, &a);
      x += randlc(&seed, &a);
      x += randlc(&seed, &a);
      key_array[i] = k * x;
   }
}

void full_verify() {
   INT_TYPE i, j;
   INT_TYPE k;
   INT_TYPE m, unique_keys;
   /*************** Clava msgError **************
   		 Array access key_array[--key_buff_ptr_global[key_buff2[i]]] which is used for writing has subscript of arrayType --key_buff_ptr_global[key_buff2[i]]
   ****************************************/
   for(i = 0; i < (1 << 23); i++)
      key_array[--key_buff_ptr_global[key_buff2[i]]] = key_buff2[i];
   j = 0;
   #pragma omp parallel for default(shared) private(i) firstprivate(key_array) reduction(+ : j)
   for(i = 1; i < (1 << 23); i++)
      if(key_array[i - 1] > key_array[i]) j++;
   if(j != 0) {
      printf("Full_verify: number of keys out of sort: %d\n", j);
   }
   else passed_verification++;
}

void rank(int iteration) {
   INT_TYPE i, j, k;
   INT_TYPE l, m;
   INT_TYPE shift = 19 - 10;
   INT_TYPE key;
   INT_TYPE min_key_val, max_key_val;
   INT_TYPE prv_buff1[524288];
   {
      key_array[iteration] = iteration;
      key_array[iteration + 10] = (1 << 19) - iteration;
      /*************** Clava msgError **************
      		 Loop Iteration number is too low
      ****************************************/
      for(i = 0; i < 5; i++)
         partial_verify_vals[i] = key_array[test_index_array[i]];
      #pragma omp parallel for default(shared) private(i)
      for(i = 0; i < (1 << 19); i++)
         key_buff1[i] = 0;
   }
   #pragma omp parallel for default(shared) private(i)
   for(i = 0; i < (1 << 19); i++)
      prv_buff1[i] = 0;
   /*************** Clava msgError **************
   		 Array access prv_buff1[key_buff2[i]] which is used for writing has subscript of arrayType key_buff2[i]
   ****************************************/
   for(i = 0; i < (1 << 23); i++) {
      key_buff2[i] = key_array[i];
      prv_buff1[key_buff2[i]]++;
   }
   /*************** Clava msgError **************
   		unsolved dependency for arrayAccess prv_buff1	 use : RWR
   ****************************************/
   for(i = 0; i < (1 << 19) - 1; i++)
      prv_buff1[i + 1] += prv_buff1[i];
   {
      #pragma omp parallel for default(shared) private(i) firstprivate(prv_buff1)
      for(i = 0; i < (1 << 19); i++)
         key_buff1[i] += prv_buff1[i];
   }
   {
      /*************** Clava msgError **************
      		Loop contains Invalid Statement -> BreakStmt#221
      ****************************************/
      for(i = 0; i < 5; i++) {
         k = partial_verify_vals[i];
         if(0 <= k && k <= (1 << 23) - 1) switch ('A') {    case 'S':    if(i <= 2) {       if(key_buff1[k - 1] != test_rank_array[i] + iteration) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    else {       if(key_buff1[k - 1] != test_rank_array[i] - iteration) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    break;    case 'W':    if(i < 2) {       if(key_buff1[k - 1] != test_rank_array[i] + (iteration - 2)) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    else {       if(key_buff1[k - 1] != test_rank_array[i] - iteration) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    break;    case 'A':    if(i <= 2) {       if(key_buff1[k - 1] != test_rank_array[i] + (iteration - 1)) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    else {       if(key_buff1[k - 1] != test_rank_array[i] - (iteration - 1)) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    break;    case 'B':    if(i == 1 || i == 2 || i == 4) {       if(key_buff1[k - 1] != test_rank_array[i] + iteration) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    else {       if(key_buff1[k - 1] != test_rank_array[i] - iteration) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    break;    case 'C':    if(i <= 2) {       if(key_buff1[k - 1] != test_rank_array[i] + iteration) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    else {       if(key_buff1[k - 1] != test_rank_array[i] - iteration) {          printf("Failed partial verification: iteration %d, test key %d\n", iteration, i);       }       else passed_verification++;    }    break; }
      }
      if(iteration == 10) key_buff_ptr_global = key_buff1;
   }
}

int main(int argc, char **argv) {
   double t1, t2;
   t1 = omp_get_wtime();
   int i, iteration, itemp;
   int nthreads = 1;
   double timecounter, maxtime;
   /*************** Clava msgError **************
   		Loop contains Invalid Statement -> BreakStmt#316
   ****************************************/
   for(i = 0; i < 5; i++)
      switch ('A') {
         case 'S':
         test_index_array[i] = S_test_index_array[i];
         test_rank_array[i] = S_test_rank_array[i];
         break;
         case 'A':
         test_index_array[i] = A_test_index_array[i];
         test_rank_array[i] = A_test_rank_array[i];
         break;
         case 'W':
         test_index_array[i] = W_test_index_array[i];
         test_rank_array[i] = W_test_rank_array[i];
         break;
         case 'B':
         test_index_array[i] = B_test_index_array[i];
         test_rank_array[i] = B_test_rank_array[i];
         break;
         case 'C':
         test_index_array[i] = C_test_index_array[i];
         test_rank_array[i] = C_test_rank_array[i];
         break;
      }
   ;
   printf("\n\n NAS Parallel Benchmarks 2.3 OpenMP C version - IS Benchmark\n\n");
   printf(" Size:  %d  (class %c)\n", (1 << 23), 'A');
   printf(" Iterations:   %d\n", 10);
   create_seq(314159265.00, 1220703125.00);
   rank(1);
   passed_verification = 0;
   if('A' != 'S') printf("\n   iteration\n");
   /*************** Clava msgError **************
   		 Loop Iteration number is too low
   ****************************************/
   for(iteration = 1; iteration <= 10; iteration++) {
      if('A' != 'S') printf("        %d\n", iteration);
      rank(iteration);
   }
   t2 = omp_get_wtime();
   full_verify();
   if(passed_verification != 5 * 10 + 1) passed_verification = 0;
   c_print_results("IS", 'A', (1 << 23), 0, 0, 10, nthreads, timecounter, ((double) (10 * (1 << 23))) / timecounter / 1000000., "keys ranked", passed_verification, "3.0 structured", "25 Nov 2023", "gcc-12", "gcc-12", "-fopenmp -lm", "-I../common -fopenmp", "-lm", "(none)", "randlc");
   printf("\n%lf\n", t2 - t1);
}
