int dot_product(int *x, int *y, int length) {
    int product = 0;

    for (int i = 0; i < length; i++) {
        product += (x[i] * y[i]);
    }
    return product;
}

/*
B is column major matrix of dimension m x n
*/
void gramschmidt_process(int **B, int m, int n) {
    int *u_k, *u_j; // represents vectors
    int numerator, denominator;

    for (int k = 0; k < n; k++) {
        u_k = B[k*m];
        
        for (int j = 0; j < k - 1; j++) {
            u_j = B[j*m];
            numerator = dot_product(u_k, u_j, m);
            denominator = dot_product(u_j, u_j, m);
            for (int l = 0; l < m; l++) {
                u_k[l] -= ((numerator / denominator)*u_j[l]);
            }
        }
    }
}