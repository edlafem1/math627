#include "lll.h"
/*
B is dimension m x n, m>=n
D is dimension n x n but since it is a diagonal, we just represent it as in dimension n
U is dimension n x n
M is dimension n x n
*/


int size_reduced(double *U, int m, int n) {
    // iterating over upper triangular matrix, excluding main diagonal of all 1's
    for (int i = 1; i < n-1; i++) {
        for (int j = 0; j < i; j++) {
            if (fabs(U[i*n + j]) > 0.5 + NUM_ERR) {
                printf("U[%i][%i]=%f\n", j, i, U[i*n+j]);
                return -1;
            }
        }
    }
    return 1;
}

int LLL_reduced(double *D, double *U, double w, int m, int n) {
    for (int i = 1; i <= n-1; i++) {
        if (w*D[i - 1] > D[i] + (U[i*n + (i - 1)])*(U[i*n + (i - 1)])*D[i - 1]) {
            return -1;
        }
    }
    return 1;
}

int closest_integer(double x) {
    if (x < 0) {
        return (int)(x - 0.5);
    }
    else {
        return (int)(x + 0.5);
    }
}

/*
B is the Gram Schmidt orthognolized vectors
*/
void reduce(double *U, double *B, double *M, int i, int j, int m, int n) {
    double gamma = (double)closest_integer(U[j*n + i]);
    /*
    M_ij = I_n - gamma*e_i*e_j'
    For any matrix A:
    A*M_ij = A where gamma*element i of every row is subtracted from every element of column j
    i.e. a_{r,j}-gamma*a_{r,i} where r is the row

    In Mathematica:
    U = U[[All,j]]-gamma*U[[All,i]]
    */
    /*
    for (int r = 0; r < m; r++) {
        if (r < n) {
            U[j*n + r] -= gamma*U[i*n + r];
            M[j*n + r] -= gamma*M[i*n + r];
        }
        B[j*m + r] -= gamma*B[i*m + r];
    }
    */
    for (int k = 0; k < i; ++k) {
        U[j*n + k] -= gamma*U[i*n + k];
    }
    U[j*n + i] -= gamma;

    for (int q = 0; q < n; q++) {
        M[j*n + q] -= gamma*M[i*n + q];
    }

    for (int r = 0; r < m; ++r) {
        B[j*m + r] -= gamma*B[i*m + r];
    }
}

/**
Precondition: X is the identity matrix I_n
Postcondition: X is the identity matrix I_n
Math condition: 2<=i<=n. Computer condition: 1<=i<=n-1
*/
void swap_restore(double *U, double *B, double *D, double *M, int i, int m, int n) {
    double u = U[i*n + (i - 1)];
    double d_hat_m1 = D[i] + u*u*D[i - 1];
    D[i] = (D[i - 1] * D[i]) / d_hat_m1;
    double epsilon = (u*D[i - 1]) / d_hat_m1;
    D[i - 1] = d_hat_m1;

    // swap columns i, i-1
    double temp;
    for (int j = 0; j <= i-2; j++) { // these conditions swap the elements in the upper triangular that are not 0
        temp = U[i*n + j];
        U[i*n + j] = U[(i - 1)*n + j];
        U[(i - 1)*n + j] = temp;
    }
    for (int j = 0; j < m; j++) {
        temp = B[i*m + j];
        B[i*m + j] = B[(i - 1)*m + j];
        B[(i - 1)*m + j] = temp;
    }
    for (int j = 0; j < n; j++) {
        temp = M[i*n + j];
        M[i*n + j] = M[(i - 1)*n + j];
        M[(i - 1)*n + j] = temp;
    }

    /*
    Doing operation U=XU
    */
    for (int j = i + 1; j < n; j++) {
        double u1 = U[j*n + (i - 1)];
        double u2 = U[j*n + i];
        U[j*n + (i - 1)] = u1*epsilon + (1 - epsilon*u)*u2;
        U[j*n + i] = u1 - u*u2;
    }
    U[i*n + (i - 1)] = epsilon;
}

void LLL(double *B, double *D, double *U, double *M, double w, int m, int n) {
#ifdef DEBUG_LLL
    printf("-----------------------------\n\n");
#endif
//    identity(M, n, n, 1); // putting this in main

    int k = 1; //math: k=2
    while (k < n) { //math: k <= n
        if (fabs(U[k*n + (k - 1)]) > 0.5+ NUM_ERR) { //Need to add NUM_ERR to account for machine error
#ifdef DEBUG_LLL
            printf("Reduce %i, %i\n", k - 1, k);
#endif
            reduce(U, B, M, k - 1, k, m, n);
#ifdef DEBUG_LLL
            printf("U: \n");
            printMatrix(U, n, n);
            printf("B: \n");
            printMatrix(B, m, n);
#endif
        }
        if (D[k] < (w - (U[k*n + (k - 1)])*(U[k*n + (k - 1)]))*D[k - 1]) {
#ifdef DEBUG_LLL
            printf("SwapRestore %i\n", k);
#endif
            swap_restore(U, B, D, M, k, m, n);
#ifdef DEBUG_LLL
            printf("U: \n");
            printMatrix(U, n, n);
            printf("B: \n");
            printMatrix(B, m, n);
            printf("D: \n");
            printMatrix(D, n, 1);
#endif
            k = max(k - 1, 1); //math: k=max(k-1,2)
        }
        else {
            for (int i = k - 2; i >= 0; i--) { //math: i=k-2 down to 1
                if (fabs(U[k*n + i])>0.5+ NUM_ERR) { // need this check for machine error
#ifdef DEBUG_LLL
                    printf("reduce %i, %i\n", i, k);
#endif
                    reduce(U, B, M, i, k, m, n);
#ifdef DEBUG_LLL
                    printf("U: \n");
                    printMatrix(U, n, n);
                    printf("B: \n");
                    printMatrix(B, m, n);
#endif
                }
            }
            k++;
        }
    }
#ifdef DEBUG_LLL
    printf("-----------------------------\n\n");
#endif
}