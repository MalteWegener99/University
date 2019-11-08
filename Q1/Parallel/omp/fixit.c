﻿/******************************************************************************
* FILE: fixit.c
* 
*   This very simple program contains errors. Find them and fix.
* 
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N 500

int main (int argc, char *argv[]) 
{
int nthreads, tid, i, j;
double a[N][N];
/*
/* Fork a team of threads */
#pragma omp parallel shared(nthreads) private(i,j,tid,a)
  {

  /* Obtain/print thread info */
  tid = omp_get_thread_num();
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);

  /* Each thread works on its own private copy of the array */

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j] = tid + i + j;
  /* For confirmation */
  printf("Thread %d done. Last element= %lf\n",tid,a[N-1][N-1]);
     /* %d, %lf - print a decimal number and a long floating (double) number repectively */

  }  /* All threads join master thread and disband */

}



