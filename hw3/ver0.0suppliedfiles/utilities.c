#include "utilities.h"

/*
Computes dot product of two n-length column vectors and returns result
*/
double dot_product_parallel(double *x, double *y, int n, int id, int np) {
	int l_n = n / np; // how many products each process will compute
	double l_sum = 0; // the local sum of the products each process computes

	for (int i = id*l_n; i < (id + 1)*l_n; i++) {
		// id*l_n is the starting element for this process, (id+1)*l_n is the starting element for the next process
		l_sum += (x[i] * y[i]);
	}

	double dot_product = 0;
	MPI_Allreduce(&l_sum, &dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// to test that all processes have same value for dot product:

	char *message;
	sprintf(message, "Process %i has dot product % 09.9f", dot_product);
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

	return dot_product;
}