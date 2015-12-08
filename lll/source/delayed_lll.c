#include "delayed_lll.h"

void reduceSwapRestore(int i, int gamma, double *B, double *D, double *U, double *M, int m, int n) {
    double u = U[i*n + (i - 1)];
    double d_hat_m = D[i] + (u - gamma)*(u - gamma)*D[i - 1];
    D[i] = (D[i - 1] * D[i]) / d_hat_m;

    double epsilon = ((u - gamma)*D[i - 1]) / d_hat_m;
    D[i - 1] = d_hat_m;

    // Update i-1 and i columns of B
    double tempB;
    for (int k = 0; k < m; ++k) {
        tempB = B[(i - 1)*m + k];
        B[(i - 1)*m + k] = B[i*m + k] - gamma*tempB;
        B[i*m + k] = tempB;
    }

    //Update i-1 and i columns of M
    double tempM;
    for (int k = 0; k < n; ++k) {
        tempM = M[(i - 1)*n + k];
        M[(i - 1)*n + k] = M[i*m + k] - gamma*tempM;
        M[i*m + k] = tempM;
    }

    // Update i-1 and i columsn of U
    double tempU;
    for (int k = 0; k <= i - 2; ++k) {
        tempU = U[(i - 1)*n + k];
        U[(i - 1)*n + k] = U[i*n + k] - gamma*tempU;
        U[i*n + k] = tempU;
    }

    // Do U=X^-1 * U
    u = U[i*n + (i - 1)];
    double u1, u2;
    for (int k = i + 1; k < n; ++k) {
        u1 = U[k*n + (i - 1)];
        u2 = U[k*n + i];
        U[k*n + (i - 1)] = u1*epsilon + (1 - epsilon*u + gamma*epsilon)*u2;
        U[k*n + i] = u1 + (gamma - u)*u2;
    }
    U[i*n + (i - 1)] = epsilon;
}

void delayed_LLL(double *B, double *D, double *U, double *M, double w, int m, int n) {
    identity(M, n, n, 1);

    int k = 1;
    double gamma;
    
    while (k < n) {
        gamma = closest_integer(U[k*n + (k - 1)]);
        if (D[k] < (w - 
                        (U[k*n + (k - 1)] - gamma)
                       *(U[k*n + (k - 1)] - gamma)
                   )*D[k - 1]) {
            reduceSwapRestore(k, gamma, B, D, U, M, m, n);
            k = max(k - 1, 1);
        }
        else {
            ++k;
        }
    }
    
    for (k = 1; k < n; ++k) {
        for (int i = k - 1; i >= 0; --i) {
            if (fabs(U[k*n + i])>0.5 + NUM_ERR) {
                reduce(U, B, M, i, k, m, n);
            }
        }
    }
}