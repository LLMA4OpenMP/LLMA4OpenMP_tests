/*
Copyright (c) 2017, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory
Written by Chunhua Liao, Pei-Hung Lin, Joshua Asplund,
Markus Schordan, and Ian Karlin
(email: liao6@llnl.gov, lin32@llnl.gov, asplund1@llnl.gov,
schordan1@llnl.gov, karlin1@llnl.gov)
LLNL-CODE-732144
All rights reserved.
This file is part of DataRaceBench. For details, see
https://github.com/LLNL/dataracebench. Please also see the LICENSE file
for our additional BSD notice.
Redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the disclaimer below.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the disclaimer (as noted below)
  in the documentation and/or other materials provided with the
  distribution.
* Neither the name of the LLNS/LLNL nor the names of its contributors
  may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/
/* 
Only the outmost loop can be parallelized. 
  
*/
#include <omp.h> 
int n = 100;
int m = 100;
double b[100][100];

int init()
{
  int i;
  int j;
  int k;
  
    #pragma omp parallel for private(j)
  for (i = 0; i <= n - 1; i += 1) {
      for (j = 0; j <= m - 1; j += 1) {
        b[i][j] = (i * j);
      }
    }
  return 0;
}

void foo()
{
  int i;
  int j;
  
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    //#pragma omp parallel for private(j)
    #pragma omp parallel for private(j)
    //in 38 and 39 ChatGPT was wrong more than half of the times
    //in this case the loop is just slightly different than 38 & 39, but ChatGPT get it right all the time
  for (i = 0; i <= n - 1; i += 1) {
        for (j = 0; j <= m - 1 - 1; j += 1) {
          b[i][j] = b[i][j + 1];
        }
      }
}

int print()
{
  int i;
  int j;
  int k;
  for (i = 0; i <= n - 1; i += 1) {
    for (j = 0; j <= m - 1; j += 1) {
      printf("%lf\n",b[i][j]);
    }
  }
  return 0;
}

int main()
{
  init();
  foo();
  print();
  return 0;
}