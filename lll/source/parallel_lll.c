#include "parallel_lll.h"

void parallel_LLL(double *B, double *D, double *U, double *M, double w, int m, int n, int id, int np) {

    int k = 1;
    int gamma;
    int f = 0;
    while (f == 0) {
        f = 1;
        // Even math k (odd in computer terms). Math says start at 2, which is 1 in computer terms
        // Begin Parallel
        for (int k = 2 + id * (n / np) - 1; k < 1  + (id + 1) * (n / np) - 1 && k < n; k += 2) {
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
#ifdef MPI_INCLUDE
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        // Odd math k (even in computer terms). Math says start at 3, which is 2 in computer terms
        // Begin Parallel
        for (int k = 3 + id * (n / np) - 1; k < 3 + (id + 1) * (n / np) - 1 && k < n; k += 2) {
            gamma = closest_integer(U[k*n + (k - 1)]);
            if (D[k] < (w -
                (U[k*n + (k - 1)] - gamma)
                *(U[k*n + (k - 1)] - gamma)
                )*D[k - 1]) {
                f = 0;
                reduceSwapRestore(k, gamma, B, D, U, M, m, n);
            }
        }
#ifdef MPI_INCLUDE
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }


    // End Parallel
    // NEED TO GET MATRICES ON ALL PARTS NOW
    int i, j, start;
    for (k = 2 * n - 3; k >= 1; k--) {
        if (k <= n - 1) {
            start = 1;
        }
        else {
            start = k - n + 2;
        }
        for (i = start; i < (k + 3) / 2; i++) {
            j = k + 2 - i;
#ifdef DEBUG_LLL
            printf("U[%i][%i]=%lf\n", i, j, U[(j - 0)*n + (i - 0)]);
#endif
            if (fabs(U[(j-1)*n+(i-1)]) > 0.5+NUM_ERR) {
                reduce(U, B, M, i-1, j-1, m, n);
            }
        }
    }
}
