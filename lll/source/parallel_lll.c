#include "parallel_lll.h"

void parallel_LLL(double *B, double *D, double *U, double *M, double w, int m, int n, int id, int np) {
#ifdef MPI_INCLUDE
    int k = 1;
    int gamma;
    int f = 0;

    /*
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
    */
    while (f == 0) {
        f = 1;
        // Even math k (odd in computer terms). Math says start at 2, which is 1 in computer terms
        // Begin Parallel
        for (int k = 1 + id * (n / np); k < 0 + (id + 1) * (n / np); k += 2) {
            gamma = closest_integer(U[k*n + (k - 1)]);
            if (D[k] < (w -
                (U[k*n + (k - 1)] - gamma)
                *(U[k*n + (k - 1)] - gamma)
                )*D[k - 1]) {
                f = 0;
                reduceSwapRestore(k, gamma, B, D, U, M, m, n);
            }
        }
        // End Parallel
        MPI_Barrier(MPI_COMM_WORLD);
        // Odd math k (even in computer terms). Math says start at 3, which is 2 in computer terms
        // Begin Parallel
        for (int k = 2 + id * (n / np); k < 1 + (id + 1) * (n / np); k += 2) {
            gamma = closest_integer(U[k*n + (k - 1)]);
            if (D[k] < (w -
                (U[k*n + (k - 1)] - gamma)
                *(U[k*n + (k - 1)] - gamma)
                )*D[k - 1]) {
                f = 0;
                reduceSwapRestore(k, gamma, B, D, U, M, m, n);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // End Parallel
    // NEED TO GET MATRICES ON ALL PARTS NOW

    for (k = 1; k < n; ++k) {

        for (int i = k - 1; i >= 0; --i) {
            if (fabs(U[k*n + i])>0.5 + NUM_ERR) {
                reduce(U, B, M, i, k, m, n);
            }
        }


    }
#else
    return;
#endif
}
