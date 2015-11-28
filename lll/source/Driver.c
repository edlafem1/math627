#include "gram_schmidt.h"
#include <stdlib.h>
#include <stdio.h>

void printMatrix(double *B, int m, int n) {
    printf("Matrix:\n");
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            printf("%f   ", B[i*m + j]);
        }
        printf("\n");
    }
}

int main() {
    int m = 2;
    int n = 2;

    double *A = (double *) calloc(m*n, sizeof(double));
    double *E = (double *)calloc(m*n, sizeof(double));
    double *U = (double *)calloc(m*n, sizeof(double));
    double *D = (double *)calloc(m, sizeof(double)); // diagonal matrix
    /*
    // j is row, i is column
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            B[j + i*m] = i ^ j;
        }
    }
    */
    A[0 + 0 * m] = 1;
    A[0 + 1 * m] = 1;
    A[1 + 0 * m] = 2;
    A[1 + 1 * m] = 0;

    printMatrix(A, m, n);

    gramschmidt_process(A, E, m, n);
    printf("A:\n");
    printMatrix(A, m, n);

    printf("Q:\n");
    printMatrix(E, m, n);

    qdu_decomposition(A, E, D, U, m, n);

    printf("D:\n");
    printMatrix(D, m, 1);

    printf("U:\n");
    printMatrix(U, m, n);

    return 0;
}