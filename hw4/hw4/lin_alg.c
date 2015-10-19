#include "lin_alg.h"

void naive_inner_product(double *A, double *B, double *D, int m, int k, int n) {
    for (int j = 0; j < n; j++) {
        for (int q = 0; q < k; q++) {
            for (int i = 0; i < m; i++) {
                D[i + j*m] += A[i + q*m] * B[q + j*k];
            }
        }
    }
}

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

double frobenius_check(double *known, double *computed, int m, int n, int id, int np) {
	int l_num_elements = m*n / np;

    if (computed == NULL) {
        return euclidean_norm(known, m*n, id, np);
    }

	double *difference = allocate_double_vector(l_num_elements);

	for (int i = 0; i < l_num_elements; i++) {
		difference[i] = computed[i] - known[i];
	}
    double e_norm = euclidean_norm(difference, m*n, id, np);
    free(difference);

    return e_norm;
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
/////////////////////////////////////////////////////////////////////////////
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


    for (int i = 0; i < l_k*n; ++i) {
        printf("[%i]: Row l_B[%i]=%f\n", id, i, l_B[i]);
    }

    for (int i = 0; i < l_k*m; ++i) {
        printf("[%i]: Col l_A[%i]=%f\n", id, i, l_A[i]);
    }





    // C only matters on process 0 and should be allocated outside this function
    double *local_C = allocate_double_vector(m*n);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, l_k, 1, l_A, m, l_B, k, 0, local_C, m);
    MPI_Reduce(local_C, C, m*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    printf("[%i] local_C[1,1]=%f\n", id, local_C[1 + m * 1]);
    free(local_C);
    free(l_A);
    free(l_B);
}