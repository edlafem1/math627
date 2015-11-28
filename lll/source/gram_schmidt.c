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

/*
A is column major matrix of dimension m x n
*/
void gramschmidt_process(double *A, double *E, int m, int n) {
    double norm;

    double *vi, *vj;

    for (int i = 0; i < n; i++) {
        vi = &(A[i*m]);

        norm = euclidean_norm(vi, m);
        for (int iindex = 0; iindex < m; iindex++) {
            E[i*m + iindex] = vi[iindex] / norm;
        }

        for (int j = i + 1; j < n; j++) {
            vj = &(A[j*m]);
            vj_projection(&(E[i*m]), vj, m); //comment says remove component in direction vi
        }

    }
}

void qr_decomposition(double *A, double *E, int m, int n) {

}