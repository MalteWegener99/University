#include <stdio.h>
#include <omp.h>

int main(){
    printf("Hello World\n");

    int num_threads;
    int my_num;

    #pragma omp parallel private(num_threads, my_num)
    {
        my_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        
        printf("I am %2d of %d\n", my_num+1, num_threads);
    }
}