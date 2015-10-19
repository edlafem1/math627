#include "main.h"

int main(int argc, char *argv[])
{
    int id, np, processor_name_len;
#ifdef PARALLEL
    char processor_name[MPI_MAX_PROCESSOR_NAME];
#endif
    double *l_x, *l_y, *temp_nvector, *y;

    double *l_A, *l_B, *l_C, *l_D;

#ifdef PARALLEL
    MPI_Init(&argc, &argv);
    /* Check processes: */
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Get_processor_name(processor_name, &processor_name_len);
#else
    np = 1;
    id = 0;
#endif
    /* Do stuff here */
    int m, k, n;
    if (argc != 4) {
        printf("Usage: %s [m] [k] [n]. Using default values 3, 4, 3.\n", argv[0]);
        m = 3;
        k = 4;
        n = 3;
    }
    else {
        m = (int)atof(argv[1]);
        k = (int)atof(argv[2]);
        n = (int)atof(argv[3]);
    }
    if (id == 0) {
        double start_time, end_time;

        double *A = allocate_double_vector(m*k);
        double *B = allocate_double_vector(k*n);
        double *C = allocate_double_vector(m*n);
        setABC_example(A, B, C, m, k, n);
        printf("np=%i\n", np);
        /*
        printf("Matrix A:\n");
        print_Matrix(A, m, k, id, np);
        printf("Matrix B:\n");
        print_Matrix(B, k, n, id, np);
        printf("Matrix C:\n");
        print_Matrix(C, m, n, id, np);
    */        
    
        // begin compute
        printf("\n");
        double *D = allocate_double_vector(m*n);
        /*
        start_time = MPI_Wtime();
        blas3_inner_product(A, B, D, m, k, n);
        end_time = MPI_Wtime();
        printf("BLAS3:\n");
        //print_Matrix(D, m, n, id, np);
        printf("Frobenius Norm: %f\n", frobenius_check(D, C, m, n, id, np));
        printf("Elapsed Time: %f\n", end_time - start_time);
        printf("\n");

        free(D);
        D = allocate_double_vector(m*n);
        start_time =  MPI_Wtime();
        blas2_inner_product(A, B, D, m, k, n);
        end_time =  MPI_Wtime();
        printf("BLAS2:\n");
        //print_Matrix(D, m, n, id, np);
        printf("Frobenius Norm: %f\n", frobenius_check(D, C, m, n, id, np));
        printf("Elapsed Time: %f\n", end_time - start_time);
        printf("\n");

        free(D);
        D = allocate_double_vector(m*n);
        start_time =  MPI_Wtime();
        blas1_inner_product(A, B, D, m, k, n);
        end_time =  MPI_Wtime();
        printf("BLAS1:\n");
        //print_Matrix(D, m, n, id, np);
        printf("Frobenius Norm: %f\n", frobenius_check(D, C, m, n, id, np));
        printf("Elapsed Time: %f\n", end_time - start_time);
        printf("\n");
        */
        free(D);
        D = allocate_double_vector(m*n);
        start_time =  MPI_Wtime();
        naive_inner_product(A, B, D, m, k, n);
        end_time =  MPI_Wtime();
        printf("Naive:\n");
        //print_Matrix(D, m, n, id, np);
        printf("Frobenius Norm: %f\n", frobenius_check(D, C, m, n, id, np));
        printf("Elapsed Time: %f\n", end_time - start_time);
        printf("\n");




        free(A);
        free(B);
        free(C);
        free(D);
    }



#ifdef PARALLEL
    MPI_Finalize();
#endif
}
