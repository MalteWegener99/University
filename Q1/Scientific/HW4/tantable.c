#include <stdio.h>
#include <time.h>

double tan(double f);
#define    pi2 1.57079632679489661923

int main(int argc, char* argv[]){
    double x = pi2/10;
    for(int i = 0; i <= 10; i++){
        printf("tan(%1.9f) = %4e \n", x*i, tan(x*i));
    }
    double tmp;
    clock_t before = clock();
    for(int i = 0; i <= 10000; i++){
        tmp = tan(x*i);
    }
    clock_t difference = clock() - before;
    int msec = difference * 1000 / CLOCKS_PER_SEC;
    printf("Homemade took %i ms\n", msec);
}
