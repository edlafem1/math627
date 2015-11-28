int size_reduced(double *U, int m, int n) {
    for (int i = 0; i < n-1; i++) {
        for (int j = i + 1; j <= n-1; j++) {
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
    return (int)(x + 0.5);
}

/*
B is the vector E in Driver. i.e. It is the Gram Schmidt orthognolized vectors
*/
void reduce(double *U, double *B, double *M, int i, int j, int m, int n) {
    double gamma = (double)closest_integer(U[j*m + i]);
    /*
    M_ij = I_n - gamma*e_i*e_j'
    For any matrix A:
    A*M_ij = A where gamma*element i of every row is subtracted from every element of column j
    i.e. a_{r,j}-gamma*a_{r,i} where r is the row
    */

    for (int r = 0; r < m; r++) {
        U[j*m + r] -= gamma*U[i*m + r];
        B[j*m + r] -= gamma*B[i*m + r];
        M[j*m + r] -= gamma*M[i*m + r];
    }
}