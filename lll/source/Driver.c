#include <stdlib.h>
#include <stdio.h>

#include "lin_alg.h"
#include "gram_schmidt.h"
#include "lll.h"

int main() {
    int m = 3;
    int n = 3;

    double *A = (double *) calloc(m*n, sizeof(double));
    double *Q = (double *)calloc(m*n, sizeof(double));
    double *D = (double *)calloc(m, sizeof(double)); // diagonal matrix
    double *U = (double *)calloc(m*n, sizeof(double));
    double *M = (double *)calloc(m*n, sizeof(double));

    A[0 + 0 * m] = 1;
    A[0 + 1 * m] = -1;
    A[0 + 2 * m] = 3;
    A[1 + 0 * m] = 1;
    A[1 + 1 * m] = 0;
    A[1 + 2 * m] = 5;
    A[2 + 0 * m] = 1;
    A[2 + 1 * m] = 2;
    A[2 + 2 * m] = 6;

    printf("A:\n");
    printMatrix(A, m, n);

    gramschmidt_process(A, Q, m, n);
    /*
    printf("A:\n");
    printMatrix(A, m, n);
    */

    qdu_decomposition(A, Q, D, U, m, n);

    printf("Q:\n");
    printMatrix(Q, m, n);

    printf("D:\n");
    printMatrix(D, m, 1);

    printf("U:\n");
    printMatrix(U, m, n);

    double w = .75;
    LLL(Q, D, U, M, w, m, n);
    printf("M:\n");
    printMatrix(M, m, n);

    printf("D:\n");
    printMatrix(D, m, 1);

    printf("U:\n");
    printMatrix(U, m, n);

    printf("Q:\n");
    printMatrix(Q, m, n);

    printf("Is size reduced? %s\n", (size_reduced(U, m, n)==1) ? "yes" : "no");

    printf("Is LLL reduced? %s\n", (LLL_reduced(D, U, w, m, n)==1) ? "yes" : "no");

    free(A);
    free(Q);
    free(D);
    free(U);
    free(M);
    return 0;
}