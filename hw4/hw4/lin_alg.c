#include "lin_alg.h"

void naive_inner_product(double *A, double *B, double *D, int m, int k, int n) {
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			double sum = 0;
			for (int q = 0; q < k; q++) {
				sum += A[i + q*m] * B[q + j*k];
			}
			D[i + j*m] = sum;
		}
	}
}

void blas1_inner_product(double *A, double *B, double *C, int m, int k, int n) {
    // uses cblas_ddot function
    //rows of A(m by k) dot columns of B (k by n)
    // incA = m
    for (int res_col = 0; res_col < n; ++res_col) {
        for (int res_row = 0; res_row < m; ++res_row) {
            C[res_row + n*res_col] = cblas_ddot(k, &(A[res_row]), m, &(B[res_col*k]), 1);
        }
    }
}

double frobenius_norm(double *known, double *computed, int m, int n, int id, int np) {
	int l_num_elements = m*n / np;

	double *difference = allocate_double_vector(l_num_elements);

	for (int i = 0; i < l_num_elements; i++) {
		difference[i] = computed[i] - known[i];
	}

	return euclidean_norm(difference, m*n, id, np);
}