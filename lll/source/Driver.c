#include <stdlib.h>
#include <stdio.h>

#include "lin_alg.h"
#include "gram_schmidt.h"
#include "lll.h"
#include "delayed_lll.h"

void get_B_values(double *B, int m, int n, char *filename) {
    /*
    Mathematica and MuPad print something like this:
    ( a  b  c )
    ( d  e  f )
    ( g  h  i )
    Which is interpreted as list of vectors.
    We think of it like this:
    [ a  d  g ]
    [ b  d  h ]
    [ c  f  i ]
    Which needs to be in the input file looking like the first representation
    */
    FILE *file = fopen(filename, "r");
    int i = 0;
    int j;
    while (i < m*n) {
        fscanf(file, "%i", &(j));
        //printf("%i\n", j);
        B[i] = (double)j;
        ++i;
    }
    fclose(file);
}

int main() {
    /*
    Memory Estimates:
    4n^2 + n
    With 64 GB, n=46341 max
    */


    int m = 8;
    int n = 8;

    double *B = (double *) calloc(m*n, sizeof(double));
    double *Q = (double *)calloc(m*n, sizeof(double));
    double *D = (double *)calloc(n, sizeof(double)); // diagonal matrix
    double *U = (double *)calloc(n*n, sizeof(double));
    double *M = (double *)calloc(n*n, sizeof(double));

    /*
    To make vectors in Matlab:
    M = randi(r,m,n)
    Creates m rows of linearly independent n-vectors with values x between 1<=x<=r.
    To write to file:
    dlmwrite('filename',M) will be a comma deliminated file of the values
    */
    /*
    B[0 + 0 * m] = 1;
    B[0 + 1 * m] = -1;
    B[0 + 2 * m] = 3;
    B[1 + 0 * m] = 1;
    B[1 + 1 * m] = 0;
    B[1 + 2 * m] = 5;
    B[2 + 0 * m] = 1;
    B[2 + 1 * m] = 2;
    B[2 + 2 * m] = 6;

    B[3 + 0 * m] = 1;
    B[3 + 1 * m] = 0;
    B[3 + 2 * m] = 0;
    */
    char *filename = "input.txt";
    get_B_values(B, m, n, filename);
    
    printf("Initial Basis:\n");
    printMatrix(B, m, n);

    gramschmidt_process(B, Q, m, n);
    
    printf("Q:\n");
    printMatrix(Q, m, n);
    

    qdu_decomposition(B, Q, D, U, m, n);

    
    printf("D:\n");
    printMatrix(D, n, 1);
    
    free(Q);
    
    printf("U:\n");
    printMatrix(U, n, n);
    

    double w = .75;
    delayed_LLL(B, D, U, M, w, m, n);
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
    printf("M: \n");
    printMatrix(M, n, n);

    printf("U:\n");
    printMatrix(U, n, n);
    printf("Is size reduced? %s\n", (size_reduced(U, m, n)==1) ? "yes" : "no");

    printf("Is LLL reduced? %s\n", (LLL_reduced(D, U, w, m, n)==1) ? "yes" : "no");

    free(B);
    free(D);
    free(U);
    free(M);
    return 0;
}