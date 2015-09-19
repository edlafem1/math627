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
	double norm = sqrt(dot_product_parallel(l_x, l_x, n, id, np));

	return norm;
}

void matrix_vector_mult_parallel(double *l_y, double *l_A, double *l_x, int n, int id, int np) {
	int l_n = n / np;
	
}


/*
Every process will send a message to process 0 which prints a message containing the operation and value that individual processes currently have.
Used to verify that all processes have the same result, or that they have a value unique to their id.
The double printed does not have all the digits required for scientific precision, only enough to verify sameness or differentess and still be easy to read.
*/
void print_result_every_process(char *operation, double value, int id, int np) {
	char *message;
	asprintf(message, "Process %i: %s % 09.9f", id, operation, value);
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
}