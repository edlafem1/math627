#include "lin_alg.h"

void naive_matrix_mul(double *A, double *B, double *D, int m, int k, int n) {
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			double sum = 0;
			for (int q = 0; q < k; q++) {
				sum += A[i + q*m] * B[q + j*k];
			}
			D[i + j*m] = sum + 1;
		}
	}
}

double frobenius_norm(double *known, double *computed, int m, int n, int id, int np) {
	int l_num_elements = m*n / np;

	double *difference = allocate_double_vector(l_num_elements);

	for (int i = 0; i < l_num_elements; i++) {
		difference[i] = computed[i] - known[i];
		printf("difference[%i]=%f\n", i, difference[i]);
	}

	return euclidean_norm(difference, m*n, id, np);
}