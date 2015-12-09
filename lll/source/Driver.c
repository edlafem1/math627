#include <stdlib.h>
#include <stdio.h>

#include "lin_alg.h"
#include "gram_schmidt.h"
#include "lll.h"
#include "delayed_lll.h"

#define DATA_FOLDER "~/student_user/bases/"

/**
Reads in a column oriented matrix of dimension m x n from file referenced by filename into B.
The data in the file must be arranged such that each group of m successive integers is one 
column of B. Note: every element MUST be an integer.

Returns 0 on success, -1 on error.
*/
int get_Initial_Basis(double *B, int m, int n, char *filename) {
    /*
    To make vectors in Matlab:
    M = randi(r,m,n)
    Creates m rows of linearly independent n-vectors with values x between 1<=x<=r.
    To write to file:
    dlmwrite('filename',M) will be a comma deliminated file of the values

    Mathematica and MuPad print something like this:
    ( a  b  c )
    ( d  e  f )
    ( g  h  i )
    Which is interpreted as list of vectors.
    We think of it like this:
    [ a  d  g ]
    [ b  d  h ]
    [ c  f  i ]
    Which needs to be in the input file looking like the first representation.
    */
    printf("filename: %s\n", filename);
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\nQuiting.\n");
        return -1;
    }
    double j;

    for (int i = 0; i < m*n; i++) {
        fscanf(file, "%lf", &(j));
        B[i] = j;
    }
    fclose(file);
    return 0;
}

/**
Memory Estimates assuming m == n:
4n^2 + n
With 64 GB, n=46341 max
*/
int main(int argc, char *argv[]) {
    /**
    m x n is the dimension of the lattice basis.
    */
    int m, n;

    /**
    A constant used as a goodness level.
    0.25 < w <= 1, however polynomial time guaranteed for only 0.25 < w < 1.
    The closer w gets to 1, the better the resulting basis.
    By default, we have this value at 0.75 because it is a common value and that used by
    MuPad (part of MATLAB).
    */
    double w = 0.75;

    char filename[256];
    if (argc >= 3) {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        if (m <= 0 || n <= 0 || n > m) {
            fprintf(stderr, "m and n must satisfy 0 < n <= m.\nQuiting\n");
            exit(0);
        }
        sprintf(filename, "%s%ix%i.dat", DATA_FOLDER, m, n);
        if (argc == 4) {
            w = atof(argv[3]);
            if (w <= .25 || w >= 1) {
                fprintf(stderr, "w must be > 0.25 and < 1.\nQuiting.\n");
                exit(0);
            }
        }
    }
    else {
        m = 8;
        n = 8;
        sprintf(filename, "input.txt");
    }
    /**
    B is the initial basis. 
    Dimensions m x n.
    B = [b_1, b_2, ..., b_n] where b_i are m-length vectors.
    */
    double *B = (double *) calloc(m*n, sizeof(double));
    /**
    Q is the gram-schmidt orthogonalized basis. 
    Dimensions m x n.
    Q = [b_1*, b_2*, ..., b_n*] where b_i* are orthogonal m-length vectors.
    */
    double *Q = (double *)calloc(m*n, sizeof(double));
    /**
    D is a diagonal matrix with the L2 norm of the gram-schmidt vectors on the main diagonal, zeros elsewhere. 
    Dimension is n x n, but we can represent it with just a vector of length n to save memory.
    Dimensions n x 1.
    D=diag(d_i), d_i = ||b_i*||^2
    */
    double *D = (double *)calloc(n, sizeof(double));
    /**
    U is an upper-triangular matrix with ones on the main diagonal. 
    Dimensions n x n.
    */
    double *U = (double *)calloc(n*n, sizeof(double));
    /**
    M is a unimodular matrix that relates two bases for the same lattice by C=BM where
    B is the original basis and C is the new basis. Initially, the LLL algorithm forces this to be
    the identity matrix I_n but relies on the assumption that it starts out as filled with zeros.
    Dimensions n x n.
    */
    double *M = (double *)calloc(n*n, sizeof(double));
    if (get_Initial_Basis(B, m, n, filename) != 0) {
        exit(0);
    }
    
#ifdef DEBUG_LLL
    printf("Initial Basis:\n");
    printMatrix(B, m, n);
#endif

    gramschmidt_process(B, Q, m, n);

#ifdef DEBUG_LLL
    printf("Q:\n");
    printMatrix(Q, m, n);
#endif   

    qdu_decomposition(B, Q, D, U, m, n);

#ifdef DEBUG_LLL  
    printf("D:\n");
    printMatrix(D, n, 1);
#endif

    free(Q);
    
#ifdef DEBUG_LLL
    printf("U:\n");
    printMatrix(U, n, n);
#endif

    delayed_LLL(B, D, U, M, w, m, n);

    printf("Final Basis:\n");
    printMatrix(B, m, n);

#ifdef DEBUG_LLL

    printf("D:\n");
    printMatrix(D, m, 1);

    printf("M: \n");
    printMatrix(M, n, n);

    printf("U:\n");
    printMatrix(U, n, n);
#endif

    printf("Is size reduced? %s\n", (size_reduced(U, m, n)==1) ? "yes" : "no");

    printf("Is LLL reduced? %s\n", (LLL_reduced(D, U, w, m, n)==1) ? "yes" : "no");

    free(B);
    free(D);
    free(U);
    free(M);
    return 0;
}