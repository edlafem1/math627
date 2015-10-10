#ifndef LIN_ALG
#define LIN_ALG
#include "utilities.h"
#include "memory.h"

/* Stuff */
void naive_matrix_mul(double *A, double *B, double *D, int m, int k, int n);

double frobenius_norm(double *known, double *computed, int m, int n, int id, int np);

#endif