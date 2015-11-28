#include <math.h>
#include <stdlib.h>
#include <string.h>

double dot_product(double *x, double *y, int length) {
    double product = 0;

    for (int i = 0; i < length; i++) {
        product += (x[i] * y[i]);
    }
    return product;
}

double euclidean_norm(double *x, int m) {
    return sqrt(dot_product(x, x, m));
}

void vj_projection(double *ui, double *uj, int m) {
    double numerator, denominator, quotient;
    numerator = dot_product(ui, uj, m);
    denominator = dot_product(ui, ui, m);
    quotient = numerator / denominator;

    for (int i = 0; i < m; i++) {
        uj[i] -= (quotient * ui[i]);
    }
}
void normalize(double *v, int m) {
    double norm = euclidean_norm(v, m);
    for (int i = 0; i < m; i++) {
        v[i] /= norm;
    }
}
/*
A is column major matrix of dimension m x n
Note: what if B in m x n is such that m > n? Does the bottom part of E mean anything? We want E n x n, but for now it is m x n because we have not answered this question.
*/
void gramschmidt_process(double *A, double *E, int m, int n) {
    double *vi, *vj;

    memcpy(E, A, sizeof(double)*m*n);

    for (int i = 0; i < n; i++) {
        vi = &(E[i*m]);

        normalize(vi, m);

        for (int j = i + 1; j < n; j++) {
            vj = &(E[j*m]);
            vj_projection(vi, vj, m); //comment says remove component in direction vi
        }

    }
}

/*
See above considerations about dimension of other matrices
*/
void qdu_decomposition(double *A, double *E, double *D, double *U, int m, int n) {
    double denominator;
    // see if maybe switching order of the loops and doing condition i <= j will be faster
    for (int i = 0; i < m; i++) {
        denominator = dot_product(&(E[i*m]), &(A[i*m]), m); // this is R[i,i]
        for (int j = i; j < m; j++) {
            if (i == j) {
                // we use this indexing because j is changing most frequently and we reduce cache misses. Note, we can switch i and j here because i==j
                U[i*m + j] = 1;
                D[i] = denominator; // D[i,j] = R[i,i], but D is only diagonal, so we represent it with just a vector
            }
            else {                
                U[j*m + i] = dot_product(&(E[i*m]), &(A[j*m]), m) / denominator; // U[i,j] = R[i,j] / R[i,i]
            }

        }
    }
}