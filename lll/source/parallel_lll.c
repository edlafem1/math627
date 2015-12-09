#include "parallel_lll.h"

void parallel_LLL(double *B, double *D, double *U, double *M, double w, int m, int n, int id, int np) {
    int k = 1;
    int gamma;

    while (k < n) {
        gamma = closest_integer(U[k*n + (k - 1)]);
        if (D[k] < (w -
            (U[k*n + (k - 1)] - gamma)
            *(U[k*n + (k - 1)] - gamma)
            )*D[k - 1]) {
            reduceSwapRestore(k, gamma, B, D, U, M, m, n);
            k = max(k - 1, 1);
        }
        else {
            ++k;
        }
    }

    // Begin parallel
    // Even k
    for (int k = 0 + id * (n / np); k < n; k += 2) {

    }
    // Odd k
    for (int k = 1 + id * (n / np); k < n; k += 2) {

    }

    // End parallel

    for (k = 1; k < n; ++k) {
        for (int i = k - 1; i >= 0; --i) {
            if (fabs(U[k*n + i])>0.5 + NUM_ERR) {
                reduce(U, B, M, i, k, m, n);
            }
        }
    }
}