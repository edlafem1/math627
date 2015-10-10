#include "utilities.h"

/* 09/12/02-10/10/02, updated 02/21/08 and 02/25/09 by Matthias K. Gobbert */

double serial_dot (double *x, double *y, int n)
{
  double dp;
#ifndef BLAS
  int i;
#endif

#ifdef BLAS
  dp = cblas_ddot(n, x, 1, y, 1);
#else
  dp = 0.0;
  for(i = 0; i < n; i++)
    dp += x[i] * y[i];
#endif

  return dp;
}


#include "utilities.h"

/*
Computes dot product of two n-length column vectors and returns result.
Note, l_x and l_y represent only the portion of the vector that this function will process
*/
double dot_product_parallel(double *l_x, double *l_y, int n, int id, int np) {
	int l_n = n / np; // how many products each process will compute
	double l_sum = 0; // the local sum of the products each process computes

	for (int l_i = 0; l_i < l_n; l_i++) {
		// id*l_n is the starting element for this process, (id+1)*l_n is the starting element for the next process
		l_sum += (l_x[l_i] * l_y[l_i]);
	}

	double dot_product = 0;
	MPI_Allreduce(&l_sum, &dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return dot_product;
}

/*
Calculates the Eucildean Norm of a column vector of length n in parallel.
Relies on the dot_product_parallel function.
*/
double euclidean_norm_parallel(double *l_x, int n, int id, int np) {
	return sqrt(dot_product_parallel(l_x, l_x, n, id, np));
}

void matrix_vector_mult_parallel(double *l_y, double *l_A, double *l_x, double *temp_y, double *y, int n, int id, int np) {
	int l_n = n / np;
	for (int i = 0; i < l_n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == 0)
				temp_y[j] = 0; // just in case this temp variable has been modified by anything else
			temp_y[j] += (l_A[j + i*n] * l_x[i]);
		}
	}

	MPI_Reduce(temp_y, y, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Scatter(y, l_n, MPI_DOUBLE, l_y, l_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}