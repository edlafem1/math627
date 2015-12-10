#ifndef LIN_ALG
#define LIN_ALG

#include "memory.h"
#include <string.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef BLAS
#include <mkl.h>
#endif

void identity(double *X, int m, int n, int zeroed);

void printMatrix(double *B, int m, int n);

double dot_product(double *x, double *y, int length);

double euclidean_norm(double *x, int m);

#endif