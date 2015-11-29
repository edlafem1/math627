#include "lll.h"
/*
B is dimension m x n, m>=n
D is dimension n x n but since it is a diagonal, we just represent it as in dimension n
U is dimension n x n
X is dimension n x n
M is dimension n x n
*/


int size_reduced(double *U, int m, int n) {
    // iterating over upper triangular matrix, excluding main diagonal of all 1's
    for (int i = 1; i < n-1; i++) {
        for (int j = 0; j < i; j++) {
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
    //printf("Questionable num %f\n", x);
    return (int)(x + 0.5);
}

/*
B is the vector E in Driver. i.e. It is the Gram Schmidt orthognolized vectors
*/
void reduce(double *U, double *B, double *M, int i, int j, int m, int n) {
    double gamma = (double)closest_integer(U[j*m + i]);
    //printf("gamma=%f\n", gamma);
    /*
    M_ij = I_n - gamma*e_i*e_j'
    For any matrix A:
    A*M_ij = A where gamma*element i of every row is subtracted from every element of column j
    i.e. a_{r,j}-gamma*a_{r,i} where r is the row

    In Mathematica:
    U = U[[All,j]]-gamma*U[[All,i]]
    */

    for (int r = 0; r < m; r++) {
        //printf("%f - %f = %f\n", U[j*m + r], gamma*U[i*m + r], U[j*m + r] - gamma*U[i*m + r]);
        U[j*m + r] -= gamma*U[i*m + r];
        B[j*m + r] -= gamma*B[i*m + r];
        if (r < m)
            M[j*n + r] -= gamma*M[i*n + r];
    }
}

/**
Precondition: X is the identity matrix I_n
Postcondition: X is the identity matrix I_n
Math condition: 2<=i<=n. Computer condition: 1<=i<=n-1
*/
void swap_restore(double *U, double *B, double *D, double *M, int i, int m, int n) {
    double u = U[i*m + (i - 1)];
    double d_hat_m1 = D[i] + u*u*D[i - 1];
    D[i] = (D[i - 1] * D[i]) / d_hat_m1;
    double epsilon = (u*D[i - 1]) / d_hat_m1;
    D[i - 1] = d_hat_m1;

    // swap columns i, i-1
    double temp;
    for (int j = 0; j <= i-2; j++) { // these conditions swap the elements in the upper triangular that are not 0
        temp = U[i*m + j];
        U[i*m + j] = U[(i - 1)*m + j];
        U[(i - 1)*m + j] = temp;
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
        double u1 = U[j*m + (i - 1)];
        double u2 = U[j*m + i];
        U[j*m + (i - 1)] = u1*epsilon + (1 - epsilon*u)*u2;
        U[j*m + i] = u1 - u*u2;
    }
    U[i*m + (i - 1)] = epsilon;
}

void LLL(double *B, double *D, double *U, double *M, double w, int m, int n) {
    printf("-----------------------------\n\n");
    
    identity(M, m, n, 1);

    int k = 1; //math: k=2
    while (k < n) { //math: k <= n
        //printf("[%i]: Top of while\n", k);
        if (fabs(U[k*m + (k - 1)]) > 0.5+1e-14) { //Need to add 1e-14 to account for machine error
            /*
            printf("[%i]: Reduce(%i, %i)\n", k, k-1, k);
            printf("B\n");
            printMatrix(B, m, n);
            */
            reduce(U, B, M, k - 1, k, m, n);
            /*
            printf("B\n");
            printMatrix(B, m, n);
            printf("D\n");
            printMatrix(D, m, 1);
            printf("U\n");
            printMatrix(U, m, n);
            */
        }
        //printf("[%i]: After first reduce if\n", k);
        if (D[k] < (w - (U[k*m + (k - 1)])*(U[k*m + (k - 1)]))*D[k - 1]) {
            /*
            printf("[%i]: SwapRestore(%i)\n", k, k);
            printf("B\n");
            printMatrix(B, m, n);
            */
            swap_restore(U, B, D, M, k, m, n);
            /*
            printf("B\n");
            printMatrix(B, m, n);
            printf("D\n");
            printMatrix(D, m, 1);
            printf("U\n");
            printMatrix(U, m, n);
            */
            k = max(k - 1, 1); //math: k=max(k-1,2)
            //printf("[%i]: After first swap restore\n", k);
        }
        else {
            //printf("[%i]: before else for loop\n", k);
            for (int i = k - 2; i >= 0; i--) { //math: i=k-2 down to 1
                if (fabs(U[k*m + i])>0.5+1e-14) { // need this check for machine error
                    /*
                    printf("[%i, %i]: reduce(%i,%i)\n", k, i, i, k);
                    printf("B\n");
                    printMatrix(B, m, n);
                    */
                    reduce(U, B, M, i, k, m, n);
                    /*
                    printf("B\n");
                    printMatrix(B, m, n);
                    printf("D\n");
                    printMatrix(D, m, 1);
                    printf("U\n");
                    printMatrix(U, m, n);
                    */
                }
                //printf("[%i, %i]: after second reduce if\n", k, i);
            }
            //printf("[%i]: After else for loop\n", k);
            k++;
        }
        //printf("[%i]: Bottom of while\n", k);
    }
}