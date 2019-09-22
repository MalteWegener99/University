#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#define CONTROL printf("Get here\n");

int main(){
    printf("Table of 23^45");
    mpz_t a, r, z;
    mpz_init(a);
    mpz_init(r);
    mpz_init(z);
    mpz_set_str(a, "23", 10);
    mpz_set_str(r, "0", 10);
    mpz_pow_ui(a,a,45);
    
    for(int i = 1; i < 11; i++){
        mpz_mul_si(r, a, i);
        mpz_set_si(z, i);
        gmp_printf("%Zd * %Zd = %Zd\n", a, z, r);
    }

    mpz_clear(a);
    mpz_clear(r);
    mpz_clear(z);
    
    printf("Table of tan with 40 digit precision\n");
    // 133 = 40 * log(10)/log(2)
    mpfr_set_default_prec(133);
    mpfr_t pi2, r1, r2;
    mpfr_init(pi2);
    mpfr_init(r1);
    mpfr_init(r2);
    mpfr_set_str(pi2, "3.14159265358979323846264338327950288419716939937510", 10, GMP_RNDD);
    mpfr_div_ui(pi2, pi2, 20u, GMP_RNDD);

    for(unsigned int i = 1; i < 11; i++){
        mpfr_mul_ui(r1, pi2, i, GMP_RNDD);
        mpfr_tan(r2, r1, GMP_RNDD);
        mpfr_printf("tan(%.40RNf) = %.40RNf \n", r1, r2);
    }

    mpfr_clear(pi2);
    mpfr_clear(r1);
    mpfr_clear(r2);
}