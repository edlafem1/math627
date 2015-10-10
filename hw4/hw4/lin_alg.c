#include "lin_alg.h"

void naive_matrix_mul(double *A, double *B, double *D, int m, int k, int n) {
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