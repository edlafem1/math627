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

// rewrite assuming l_A is a bunch of columns of A, l_x is a bunch of rows, and l_y n x 1. Use an MPI_Reduce call. NO ALLGATHER of x. Use a scatter to put result on all processes.
void matrix_vector_mult_parallel(double *l_y, double *l_A, double *l_x, double *temp_y, double *y, int n, int id, int np) {
	int l_n = n / np;
	for (int i = 0; i < l_n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == 0)
				temp_y[0] = 0;
			temp_y[j] += (l_A[j + i*n] * l_x[i]);
			printf("id=%i, i=%i, j=%i, l_A[j+i*n]=%f, l_x[i]=%f\n");
		}
	}
	//for (int j = 0; j < n; j++)
	//	printf("id: %i, j=%i, partial=%f\n", id, j, temp_y[j]);

	MPI_Reduce(temp_y, y, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Scatter(y, l_n, MPI_DOUBLE, l_y, l_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

int eigenvalue_approximation_parallel(double *lambda, double *err, double *l_x, double *l_A, double *l_y, double *temp_nvector, double *y, double tol, int itmax, int n, int id, int np) {
	int l_n = n / np;
	double norm_x;
	double lambdaold = 0;
	int iterations = 0;
	*err = 1.0 + tol; // to ensure 1 pass through loop
	*lambda = 0;
	matrix_vector_mult_parallel(l_y, l_A, l_x, temp_nvector, y, n, id, np); // just to make sure, in case we change the calling code:

	while ((*err > tol) && (iterations < itmax)) {
		++iterations; // increment iteration counter
		lambdaold = *lambda;

		//norm_y = euclidean_norm_parallel(l_y, n, id, np); // Euclidean vector norm of vector y
		memcpy(l_x, l_y, l_n * sizeof(double)); // takes place of x=y/normy but without scaling, so really just means x = y
		matrix_vector_mult_parallel(l_y, l_A, l_x, temp_nvector, y, n, id, np); // y_new = A * x = A * y_old
		*lambda = dot_product_parallel(l_x, l_y, n, id, np) / dot_product_parallel(l_x, l_x, n, id, np); // eigenvalue approximation using Rayleigh quotient with eigenvector x
		*err = fabs((*lambda - lambdaold) / *lambda);
	}

	if (iterations == itmax && id == 0) {
		printf("Max number of iterations reached. Answer given is the most recent approximation.\n");
	}

	// compute scaled eigenvector x
	norm_x = euclidean_norm_parallel(l_x, n, id, np);
	double norm_y = euclidean_norm_parallel(l_y, n, id, np);
	for (int row = 0; row < l_n; row++) {
		l_x[row] = l_x[row] / norm_x;
		l_y[row] = l_y[row] / norm_y;
	}
	return iterations;
}

void print_Square_Matrix(double *l_A, int id, int n, int np) {
	// create a matrix A for use on only process 0 to gather all local pieces into one place
	double *A;
	int destination = 0;
	if (id == 0) {
		A = (double *)calloc(n*n, sizeof(double));
	}

	int l_n = n / np;
	// PRINT n x n array
	// gather all l_A which have been set in the setup_example function into A for printing by process 0
	MPI_Gather(l_A, n*l_n, MPI_DOUBLE, A, n*l_n, MPI_DOUBLE, destination, MPI_COMM_WORLD);
	if (id == 0) {
		printf("Matrix A:\n");
		int row, col;
		for (row = 0; row < n; row++) {
			for (col = 0; col < n; col++) {
				printf("% -24.16e   ", A[row + col * n]);
			}
			printf("\n");
		}
		free(A); // only allocated on process 0
	}
}


/*
Every process will send a message to process 0 which prints a message containing the operation and value that individual processes currently have.
Used to verify that all processes have the same result, or that they have a value unique to their id.
The double printed does not have all the digits required for scientific precision, only enough to verify sameness or differentess and still be easy to read.
*/
void print_result_every_process(char *operation, double value, int id, int np) {
	char *message = (char *)malloc(100 * sizeof(char));
	sprintf(message, "%s % -24.16e", operation, value);
	if (id == 0) {
		MPI_Status status;
		printf("From 0: %s\n", message);
		for (int i = 1; i < np; i++) {
			MPI_Recv(message, 100, MPI_CHAR, i, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);
			printf("From %i: %s\n", status.MPI_SOURCE, message);
		}
	}
	else {
		MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}
	free(message);
}