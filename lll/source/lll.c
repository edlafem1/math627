#include "lll.h"
/*
B is dimension m x n, m>=n
D is dimension n x n but since it is a diagonal, we just represent it as in dimension n
U is dimension n x n
X is dimension n x n
M is dimension n x n
*/


int size_reduced(double *U, int m, int n) {
    for (int i = 0; i < n-1; i++) {
        for (int j = i + 1; j <= n-1; j++) {
            if (fabs(U[i*m + j]) > 0.5)
                return -1;
        }
    }
    return 1;
}

int LLL_reduced(double *D, double *U, double w, int m, int n) {
    for (int i = 1; i <= n-1; i++) {
        if (w*D[i - 1] > D[i] + (U[i*m + (i - 1)])*(U[i*m + (i - 1)])*D[i - 1]) {
            return -1;
        }
    }
    return 1;
}

int closest_integer(double x) {
    return (int)(x + 0.5);
}

/*
B is the vector E in Driver. i.e. It is the Gram Schmidt orthognolized vectors
*/
void reduce(double *U, double *B, double *M, int i, int j, int m, int n) {
    double gamma = (double)closest_integer(U[j*m + i]);
    /*
    M_ij = I_n - gamma*e_i*e_j'
    For any matrix A:
    A*M_ij = A where gamma*element i of every row is subtracted from every element of column j
    i.e. a_{r,j}-gamma*a_{r,i} where r is the row

    In Mathematica:
    U = U[[All,j]]-gamma*U[[All,i]]
    */

    for (int r = 0; r < m; r++) {
        U[j*m + r] -= gamma*U[i*m + r];
        B[j*m + r] -= gamma*B[i*m + r];
        M[j*m + r] -= gamma*M[i*m + r];
    }
}

/**
Precondition: X is the identity matrix I_n
Postcondition: X is the identity matrix I_n
*/
void swap_restore(double *U, double *B, double *D, double *M, double *X, int i, int m, int n) {
    double u = U[i*m + (i - 1)];
    double d_hat_m1 = D[i] + u*u*D[i - 1];
    D[i] = (D[i - 1] * D[i]) / d_hat_m1;
    double epsilon = (u*D[i - 1]) / d_hat_m1;
    D[i - 1] = d_hat_m1;
    U[i*m + (i - 1)] = epsilon;

    double temp;

    for (int j = 0; j < m; j++) {
        temp = U[i*m + j];
        U[i*m + j] = U[(i - 1)*m + j];
        U[(i - 1)*m + j] = temp;
    }
    for (int j = 0; j < m; j++) {
        temp = B[i*m + j];
        B[i*m + j] = B[(i - 1)*m + j];
        B[(i - 1)*m + j] = temp;
    }
    for (int j = 0; j < m; j++) {
        temp = M[i*m + j];
        M[i*m + j] = M[(i - 1)*m + j];
        M[(i - 1)*m + j] = temp;
    }

    X[(i - 2)*m + (i - 2)] = epsilon;
    X[(i - 2)*m + (i - 1)] = 1;
    X[(i - 1)*m + (i - 2)] = 1 - epsilon*u;
    X[(i - 1)*m + (i - 1)] = -u;

    naive_inner_product(X, U, U, n, n, n);

    // reset back to Identity matrix
    X[(i - 2)*m + (i - 2)] = 1;
    X[(i - 2)*m + (i - 1)] = 0;
    X[(i - 1)*m + (i - 2)] = 0;
    X[(i - 1)*m + (i - 1)] = 1;
}