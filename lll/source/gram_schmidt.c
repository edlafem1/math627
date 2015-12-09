#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
Computes the projection of m-vectors uj on ui.
*/
void vj_projection(double *ui, double *uj, int m) {
    double numerator, denominator, quotient;
    numerator = dot_product(ui, uj, m);
    denominator = dot_product(ui, ui, m);
    quotient = numerator / denominator;

    for (int i = 0; i < m; i++) {
        uj[i] -= (quotient * ui[i]);
    }
}

/**
Normalizes a vector of length m with the Euclidean norm.
*/
void normalize(double *v, int m) {
    double norm = euclidean_norm(v, m);
    for (int i = 0; i < m; i++) {
        v[i] /= norm;
    }
}
/*
B is a column major matrix. Dimensions m x n.
Stores orthogonalized vectors in Q. Dimensions m x n.
*/
void gramschmidt_process(double *B, double *Q, int m, int n) {
    double *vi, *vj;

    memcpy(Q, B, sizeof(double)*m*n);

    for (int i = 0; i < n; i++) {
        vi = &(Q[i*m]);

        normalize(vi, m);

        for (int j = i + 1; j < n; j++) {
            vj = &(Q[j*m]);
            vj_projection(vi, vj, m);
        }

    }
}

/*
Takes Gram-Schmidt orthogonalized vectors and computes D and U such that:
B=Q(D^1/2)U.

B is the initial basis.
Dimensions m x n.
B = [b_1, b_2, ..., b_n] where b_i are m-length vectors.

Q is the gram-schmidt orthogonalized basis.
Dimensions m x n.
Q = [b_1*, b_2*, ..., b_n*] where b_i* are orthogonal m-length vectors.

D is a diagonal matrix with the L2 norm of the gram-schmidt vectors on the main diagonal, zeros elsewhere.
Dimension is n x n, but we can represent it with just a vector of length n to save memory.
Dimensions n x 1.
D=diag(d_i), d_i = ||b_i*||^2.

U is an upper-triangular matrix with ones on the main diagonal.
Dimensions n x n.
*/
void qdu_decomposition(double *B, double *Q, double *D, double *U, int m, int n) {
    double denominator;
    // see if maybe switching order of the loops and doing condition i <= j will be faster
    for (int i = 0; i < n; i++) {
        denominator = dot_product(&(Q[i*m]), &(B[i*m]), m); // this is R[i,i]
        for (int j = i; j < n; j++) {
            if (i == j) {
                // we use this indexing because j is changing most frequently and we reduce cache misses. Note, we can switch i and j here because i==j
                U[i*n + j] = 1;
                D[i] = denominator*denominator; // (D[i,j] = R[i,i], but D is only diagonal, so we represent it with just a vector)^2 for LLL algo
            }
            else {                
                U[j*n + i] = dot_product(&(Q[i*m]), &(B[j*m]), m) / denominator; // U[i,j] = R[i,j] / R[i,i]
            }

        }
    }
}