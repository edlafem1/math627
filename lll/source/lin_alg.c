#include "lin_alg.h"

/**
Embeds the idendity matrix of dimension min(m,n) inside the matrix X
of dimensions m x n.
If zeroed != 0, we assume X is filled with zeroes and set each
element X_i,i to 1 for 0 <= i < min(m,n). Otherwise, we
iterate over each element of X

*/
void identity(double *X, int m, int n, int zeroed) {
    int min_dim = (m < n) ? m : n;
    if (!zeroed)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (i == j)
                    X[i*m + j] = 1;
                else
                    X[i*m + j] = 0;
            }
        }
    else
        for (int i = 0; i < min_dim; i++) {
            X[i*m + i] = 1;
        }
}

void printMatrix(double *B, int m, int n) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            printf("% 8.4f ", B[i*m + j]);
        }
        printf("\n");
    }
}

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





















#ifdef PARALLEL
/*
Computes dot product of two n-length column vectors and returns result.
Note, l_x and l_y represent only the portion of the vector that this function will process
*/
double parallel_dot_product(double *l_x, double *l_y, int n, int id, int np) {
    double dot_product = 0;
#ifdef BLAS
#ifndef PARALLEL
    dot_product = cblas_ddot(n, l_x, 1, l_y, 1); // increment by 1 because we are in serial only here
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
double parallel_euclidean_norm(double *l_x, int n, int id, int np) {
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
#endif


void naive_inner_product(double *A, double *B, double *D, int m, int k, int n) {
    for (int j = 0; j < n; j++) {
        for (int q = 0; q < k; q++) {
            for (int i = 0; i < m; i++) {
                D[i + j*m] += A[i + q*m] * B[q + j*k];
            }
        }
    }
}

#ifdef BLAS

void blas1_inner_product(double *A, double *B, double *C, int m, int k, int n) {
    // uses cblas_ddot function
    //rows of A(m by k) dot columns of B (k by n)
    // incA = m
    for (int res_col = 0; res_col < n; ++res_col) {
        for (int res_row = 0; res_row < m; ++res_row) {
            C[res_row + m*res_col] = cblas_ddot(k, &(A[res_row]), m, &(B[res_col*k]), 1);
        }
    }
}

void blas2_inner_product(double *A, double *B, double *C, int m, int k, int n) {
    //http://www.cs.utexas.edu/users/flame/pubs/SUMMA2d3dTOMS.pdf
    //https://software.intel.com/en-us/node/520751#94156EDE-4ADD-4830-940E-1CA5688ABE88

	for (int row = 0; row < k; ++row) {
		cblas_dger(CblasColMajor, m, n, 1, &(A[row*m]), 1, &(B[row]), k, C, m);
	}

}

void blas3_inner_product(double *A, double *B, double *C, int m, int k, int n) {
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A, m, B, k, 0, C, m);
}



void parallel_blas3_product(double *A, double *B, double *C, int m, int k, int n, int id, int np) {
    if (k % np != 0) {
        if (id == 0)
            fprintf(stderr, "k is not divisible by np.\n");

        MPI_Abort(MPI_COMM_WORLD, 1);
        
    }
    MPI_Status status;
    int l_k = k / np;
    double *l_A = allocate_double_vector(l_k * m);
    double *l_B = allocate_double_vector(l_k * k);
    MPI_Datatype block_col_t;
    MPI_Datatype block_row_t;
    
    // for blocks in B = k x n
    MPI_Type_vector(
        n,          // count = number of blocks, i.e. length of column * l_k(num rows)
        l_k,              // blocklen = number of things in each block
        k,              // stride = difference between start of blocks
        MPI_DOUBLE,     // old datatype
        &block_row_t    // new datatype
        );
    MPI_Type_commit(&block_row_t);

    // for column of A= m x k
    MPI_Type_contiguous(
        m * l_k,        // count = number of items
        MPI_DOUBLE,     // old_type = type of items
        &block_col_t    // new_mpi_type = the new datatype
        );
    MPI_Type_commit(&block_col_t);

    if (id == 0) {
        // copy correct elements from A to l_A
        memcpy(l_A, A, sizeof(double) * l_k * m);
        for (int i = 1; i < np; ++i) {
            MPI_Send(&(A[0 + m*(i*l_k)]), 1, block_col_t, i, 0, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(l_A, (m*l_k), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }


    if (id == 0) {
        // copy numbers from B to l_B
        
        for (int col = 0; col < n; ++col) {
            for (int row = 0; row < l_k; ++row) {
                l_B[row + l_k*col] = B[row + k*col];
            }
        }

        for (int i = 1; i < np; ++i) {
            MPI_Send(&(B[i*l_k]), 1, block_row_t, i, 0, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(l_B, (l_k*n), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }

    /*
    //debugging only
    for (int i = 0; i < l_k*n; ++i) {
        printf("[%i]: Row l_B[%i]=%f\n", id, i, l_B[i]);
    }

    for (int i = 0; i < l_k*m; ++i) {
        printf("[%i]: Col l_A[%i]=%f\n", id, i, l_A[i]);
    }
    */




    // C only matters on process 0 and should be allocated outside this function
    double *local_C = allocate_double_vector(m*n);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, l_k, 1, l_A, m, l_B, l_k, 0, local_C, m);
    MPI_Reduce(local_C, C, m*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    free(local_C);
    free(l_A);
    free(l_B);
}
#endif