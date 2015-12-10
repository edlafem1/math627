#include "lin_alg.h"

/**
Embeds the idendity matrix of dimension min(m,n) inside the matrix X
of dimensions m x n.
If zeroed != 0, we assume X is filled with zeroes and set each
element X_i,i to 1 for 0 <= i < min(m,n). Otherwise, we
iterate over each element of X

*/
void identity(double *X, int m, int n, int zeroed) {
    int min_dim = (m < n) ? m : n;
    if (!zeroed)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (i == j)
                    X[i*m + j] = 1;
                else
                    X[i*m + j] = 0;
            }
        }
    else
        for (int i = 0; i < min_dim; i++) {
            X[i*m + i] = 1;
        }
}

void printMatrix(double *B, int m, int n) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            printf("% 8.4f ", B[i*m + j]);
        }
        printf("\n");
    }
}

double dot_product(double *x, double *y, int length) {
    double product = 0;
#ifdef BLAS
    product = cblas_ddot(length, x, 1, y, 1);
#else
    for (int i = 0; i < length; i++) {
        product += (x[i] * y[i]);
    }
#endif
    return product;
}

double euclidean_norm(double *x, int m) {
    return sqrt(dot_product(x, x, m));
}