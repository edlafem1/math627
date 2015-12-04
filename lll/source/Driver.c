#include <stdlib.h>
#include <stdio.h>

#include "lin_alg.h"
#include "gram_schmidt.h"
#include "lll.h"

int main() {
    int m = 3;
    int n = 3;

    double *B = (double *) calloc(m*n, sizeof(double));
    double *Q = (double *)calloc(m*n, sizeof(double));
    double *D = (double *)calloc(n, sizeof(double)); // diagonal matrix
    double *U = (double *)calloc(n*n, sizeof(double));
    double *M = (double *)calloc(m*n, sizeof(double));

    B[0 + 0 * m] = 1;
    B[0 + 1 * m] = -1;
    B[0 + 2 * m] = 3;
    B[1 + 0 * m] = 1;
    B[1 + 1 * m] = 0;
    B[1 + 2 * m] = 5;
    B[2 + 0 * m] = 1;
    B[2 + 1 * m] = 2;
    B[2 + 2 * m] = 6;
    
    printf("Initial Basis:\n");
    printMatrix(B, m, n);

    gramschmidt_process(B, Q, m, n);
    /*
    printf("B:\n");
    printMatrix(B, m, n);
    */

    qdu_decomposition(B, Q, D, U, m, n);

    /*
    printf("Q:\n");
    printMatrix(Q, m, n);
    */
    free(Q);
    /*
    printf("D:\n");
    printMatrix(D, m, 1);

    printf("U:\n");
    printMatrix(U, m, n);
    */

    double w = .75;
    LLL(B, D, U, M, w, m, n);
    /*
    printf("M:\n");
    printMatrix(M, m, n);

    printf("D:\n");
    printMatrix(D, m, 1);

    printf("U:\n");
    printMatrix(U, m, n);
    */

    printf("Final Basis:\n");
    printMatrix(B, m, n);

    printf("Is size reduced? %s\n", (size_reduced(U, m, n)==1) ? "yes" : "no");

    printf("Is LLL reduced? %s\n", (LLL_reduced(D, U, w, m, n)==1) ? "yes" : "no");

    free(B);
    free(D);
    free(U);
    free(M);
    return 0;
}