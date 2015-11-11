#include "main.h"

int main(int argc, char **argv) {
    int i, j, id, np, processor_name_len;
    int maxit, idleft, idright, flag, iter;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int dimensions, N, l_N, n, l_n;
    double *l_r, *l_x, *l_p, *l_q;
    double *x, *y, *gl, *gr;
    double h, tol, relres;

    MPI_Init(&argc, &argv);

    /* Check processes: */
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Get_processor_name(processor_name, &processor_name_len);

    /* process command-line inputs: */
    if (argc != 3)
    {
        if (id == 0)
        {
            printf("Usage: \"./poisson N dim\" \n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    N = atoi(argv[1]);
    dimensions = atoi(argv[2]);

    if (dimensions == 2) {
        n = N*N;
    }
    else if (dimensions == 3) {
        n = N*N*N;
    }
    else {
        printf("Error: dimensions must equal 2 or 3");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* number of processes np must divide n: */
    if ((n % np) != 0)
    {
        if (id == 0)
        {
            printf("Error: np must divide n!\n");
            printf("  n = %d, np = %d, n%%np = %d\n", n, np, (n%np));
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* calculate size of local blocks: */
    l_N = N / np;
    l_n = n / np;

    if (id == 0)
    {
        printf("n = %d, np = %d, l_n = %d\n", n, np, l_n);
        printf("\n");
        fflush(stdout);
    }
    /*
    vectors l_x, l_r, l_p, l_q should be of
    length l_n and gl and gr should be of length l_n / l_N

    */
    x = allocate_double_vector(N);          /* x is a N length vector*/
    y = allocate_double_vector(N);          /* y is a N length vector*/
    gl = allocate_double_vector(l_n / l_N); 
    gr = allocate_double_vector(l_n / l_N);
    
    l_x = allocate_double_vector(l_n);                          
    l_r = allocate_double_vector(l_n);
    l_p = allocate_double_vector(l_n);
    l_q = allocate_double_vector(l_n);

    double start_time, end_time;

    /*Beginning of cg method: follows matlab code*/
    h = 1.0 / (N + 1.0);
    for (i = 1; i <= N; i++) {
        x[i - 1] = i*h;
        y[i - 1] = i*h;
    }

    /*Setup B*/
    setupB(l_r, x, y, l_N, N, h, id);

    tol = 1.0e-6;
    maxit = 99999;

    if (id>0) {
        idleft = id - 1;
    }
    else {
        idleft = MPI_PROC_NULL;
    }

    if (id<np - 1) {
        idright = id + 1;
    }
    else {
        idright = MPI_PROC_NULL;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    

    cg(l_x, &flag, &relres, &iter,  /*output*/
        l_r, tol, maxit,            /*input*/
        l_p, l_q, l_n, l_N, N, id, idleft, idright, np, MPI_COMM_WORLD, gl, gr);

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    
    double *full = allocate_double_vector(n);
    MPI_Gather(l_p, l_n, MPI_DOUBLE, full, l_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double diff_norm = fd_norm(full, h, N);
    /*
    if (id == 0) {
        for (i = 0; i < n; i++) {
            printf("%i: % -24.16e\n", id, full[i]);
        }
    }
    */
    free_vector(full);
    

    if (id == 0) {
        printf("%15s %15s %22s %15s %22s\n", "N", "DOF", "relres", "iter", "time");
        printf("%15d %15.f %22.16e %15d %22.16e\n", N, (double)N*N, relres, iter, (end_time - start_time));
        printf("||u-u_h||=%22.16e\n", diff_norm);
        /*
        printf("N:            %d\n", N);
        printf("DOF:          %d\n", N*(double)N);
        printf("relres:       %22.16e\n", relres);
        printf("iter:         %d\n", iter);
        printf("elapsed time: %22.16e\n", (end_time - start_time));
        */
    }


    free_vector(x);
    free_vector(y);
    free_vector(l_r);
    free_vector(gl);
    free_vector(gr);
    free_vector(l_x);
    free_vector(l_p);
    free_vector(l_q);
    MPI_Finalize();

    return 0;
}