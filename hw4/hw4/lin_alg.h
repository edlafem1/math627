#ifndef LIN_ALG
#define LIN_ALG
#include "utilities.h"
#include "memory.h"

/* Stuff */

void naive_inner_product(double *A, double *B, double *D, int m, int k, int n);

void blas1_inner_product(double *A, double *B, double *C, int m, int k, int n);

void blas2_inner_product(double *A, double *B, double *C, int m, int k, int n);

void blas3_inner_product(double *A, double *B, double *C, int m, int k, int n);

double frobenius_check(double *known, double *computed, int m, int n, int id, int np);



#endif