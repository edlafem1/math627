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

/*
Computes dot product of two n-length column vectors and returns result.
Note, l_x and l_y represent only the portion of the vectors that this function will process
*/
double parallel_dot_product(double *l_x, double *l_y, int n, int id, int np);

/*
Calculates the Eucildean Norm of a column vector of length n in parallel.
Relies on the dot_product function.
Note, l_x represents only the portion of the vector that this function will process
*/
double parallel_euclidean_norm(double *l_x, int n, int id, int np);

/*
Performs the multiplication y = Ax where A is a matrix and x is a column vector.
Assumes dimensions of A are n x n and dimensions of x and y are n x 1.
*/
void matrix_vector_mult(double *l_y, double *l_A, double *l_x, double *temp_y, double *y, int n, int id, int np);

void naive_inner_product(double *A, double *B, double *D, int m, int k, int n);

#ifdef BLAS
void blas1_inner_product(double *A, double *B, double *C, int m, int k, int n);

void blas2_inner_product(double *A, double *B, double *C, int m, int k, int n);

void blas3_inner_product(double *A, double *B, double *C, int m, int k, int n);

void parallel_blas3_product(double *A, double *B, double *C, int m, int k, int n, int id, int np);
#endif


#endif