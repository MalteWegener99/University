/*
 * hello.c
 */

#include <stdio.h>
#include <omp.h>


int main()
{
    int tid, num_threads;
    
#pragma omp parallel private(tid, num_threads)
    {
        tid=omp_get_thread_num();
		num_threads = omp_get_num_threads();
        printf("Hello from thread No %d from %d total threads.\n",tid, num_threads);
    }
        return 0;
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
