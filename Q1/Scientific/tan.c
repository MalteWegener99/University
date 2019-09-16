#include <math.h>
#include <stdio.h>

#ifndef DOUBLE_PREC
    #define TYPE float
    #define FN tanf
#else
    #define TYPE double
    #define FN tan
#endif

int main(int argc, char* argv[]){
    TYPE x = M_PI/20;
    for(int i = 0; i <= 10; i++){
        printf("tan(%1.9f) = %4e \n", x*i, FN(x*i));
    }
}
