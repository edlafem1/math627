#ifndef UTILITIES_H
	#define UTILITIES_H
	#ifndef MAIN_H
		#include "main.h"
	#endif

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

/*
Uses the Power method to approximate the largest eigenvalue of the matrix A in parallel. l_y should represent y = A*x.
Returns the number of iterations taken, stores the eigenvale approximation in lambda and this process's chunk of the eigenvector approximation in l_x and the error in err.
*/
int eigenvalue_approximation_parallel(double *lambda, double *err, double *l_x, double *l_A, double *l_y, double *temp_nvector, double *y, double tol, int itmax, int n, int id, int np);

/*
Prints the square matrix A with n rows and n columns distributed over np processes.
*/
void print_Square_Matrix(double *l_A, int id, int n, int np);

/*
Every process will send a message to process 0 which prints a message containing the operation and value that individual processes currently have.
Used to verify that all processes have the same result, or that they have a value unique to their id.
The double printed does not have all the digits required for scientific precision, only enough to verify sameness or differentess and still be easy to read.
*/
void print_result_every_process(char *operation, double value, int id, int np);
#endif