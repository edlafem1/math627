#include "utilities.h"

/*
Computes dot product of two n-length column vectors and returns result.
Note, l_x and l_y represent only the portion of the vector that this function will process
*/
double dot_product(double *l_x, double *l_y, int n, int id, int np) {
    double dot_product = 0;
#ifdef BLASXXX
#ifndef PARALLEL
    dot_product = cblas_ddot(n, l_x, 1, l_y, 1); // increment by 1 because we are in serial
#endif
#else
    int l_n = n / np; // how many products each process will compute
    double l_sum = 0; // the local sum of the products each process computes

    for (int l_i = 0; l_i < l_n; l_i++) {
        // id*l_n is the starting element for this process, (id+1)*l_n is the starting element for the next process
        l_sum += (l_x[l_i] * l_y[l_i]);
    }

#ifdef PARALLEL
    if (np > 1)
        MPI_Allreduce(&l_sum, &dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    else {
        dot_product = l_sum;
    }
#else
    dot_product = l_sum;
#endif

#endif
    return dot_product;
}

/*
Calculates the Eucildean Norm of a column vector of length n.
Relies on the dot_product function.
*/
double euclidean_norm(double *l_x, int n, int id, int np) {
#ifdef BLASXXX
    return cblas_dnrm2(n, l_x, 1); // euclidean norm BLAS level 1
#else
    return sqrt(dot_product(l_x, l_x, n, id, np));
#endif
}

void matrix_vector_mult(double *l_y, double *l_A, double *l_x, double *temp_y, double *y, int n, int id, int np) {
    int l_n = n / np;
    for (int i = 0; i < l_n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0)
                temp_y[j] = 0; // just in case this temp variable has been modified by anything else
            temp_y[j] += (l_A[j + i*n] * l_x[i]);
        }
    }
#ifdef PARALLEL
    MPI_Reduce(temp_y, y, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Scatter(y, l_n, MPI_DOUBLE, l_y, l_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

// prints an mxn matrix
void print_Matrix(double *l_matrix, int m, int n, int id, int np) {
    // create a matrix A for use on only process 0 to gather all local pieces into one place
    double *A;
#ifdef PARALLEL
    int destination = 0;
    if (np > 1) {
        if (id == 0) {
            printf("trying to allocate for print matrix\n");
            A = allocate_double_vector(m*n);
            printf("allocated for print matrix\n");
        }
        int l_n = n / np;

        // gather all l_A which have been set in the setup_example function into A for printing by process 0
        MPI_Gather(l_matrix, m*l_n, MPI_DOUBLE, A, m*l_n, MPI_DOUBLE, destination, MPI_COMM_WORLD);
    }
    else {
        A = l_matrix;
    }
#else
    A = l_matrix;
#endif
#ifdef PARALLEL
    if (id == destination) {
#endif
        int row, col;
        for (row = 0; row < m; row++) {
            for (col = 0; col < n; col++) {
                printf("% -24.16e   ", A[row + col * m]);
            }
            printf("\n");
        }
#ifdef PARALLEL
        if (np > 1)
            free(A); // only allocated on process destination=0
    }
#endif
}