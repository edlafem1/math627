#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef BLAS
#include <mkl.h>
#endif

double serial_dot (double *x, double *y, int n);

/*
Computes dot product of two n-length column vectors and returns result.
Note, l_x and l_y represent only the portion of the vectors that this function will process
*/
double dot_product_parallel(double *l_x, double *l_y, int n, int id, int np);

/*
Calculates the Eucildean Norm of a column vector of length n in parallel.
Relies on the dot_product_parallel function.
Note, l_x represents only the portion of the vector that this function will process
*/
double euclidean_norm_parallel(double *l_x, int n, int id, int np);

/*
Performs the multiplication y = Ax where A is a matrix and x is a column vector.
Assumes dimensions of A are n x n and dimensions of x and y are n x 1.
*/
void matrix_vector_mult_parallel(double *l_y, double *l_A, double *l_x, double *temp_y, double *y, int n, int id, int np);
#endif

