#include "utilities.h"
#define M_PI 3.14159265358979323846

void setupB(double *l_r, double *x, double *y, int l_N, int N, double h, int id) {
    // copy code from Gobbert notes #18
    // l_r is where we store B
    int l_j, j, i;
    for(l_j = 0; l_j < l_N; l_j++) {
        j = l_j + id*l_N;
        for (i = 0; i < N; i++) {
            l_r[i + N*l_j] = (h*h)*f(x[i], y[i]);
        }
    }

}

double f(double X, double Y) {
    return ((-2.0*(M_PI*M_PI))*(cos(2.0*M_PI*X)*((sin(M_PI*Y))*(sin(M_PI*Y))) + ((sin(M_PI*X))*(sin(M_PI*X)))*cos(2.0*M_PI*Y))); // << THIS IS MATLAB CODE
}

/*
Computes dot product of two n-length column vectors and returns result.
Note, l_x and l_y represent only the portion of the vector that this function will process
*/
double parallel_dot(double *l_x, double *l_y, int l_n, MPI_Comm comm) {
    int id, np;
    double dot_product = 0;
    double l_sum = 0; // the local sum of the products each process computes

    for (int l_i = 0; l_i < l_n; l_i++) {
        // id*l_n is the starting element for this process, (id+1)*l_n is the starting element for the next process
        l_sum += (l_x[l_i] * l_y[l_i]);
    }

    if (np > 1)
        MPI_Allreduce(&l_sum, &dot_product, 1, MPI_DOUBLE, MPI_SUM, comm);
    else {
        dot_product = l_sum;
    }
    return dot_product;
}

double *allocate_double_vector(int n)
{
    double *x;

    x = (double*)calloc(n, sizeof(double));

    if (x == NULL)
    {
        fprintf(stderr, "Problem allocating memory for vector\n");
#ifdef PARALLEL
        MPI_Abort(MPI_COMM_WORLD, 1);
#else
        exit(1);
#endif
    }

    return x;
}

void free_vector(void *x)
{
    free(x);
}



/*
Need to make my parallel_dot function work with his cg.c
Make Ax function do blocks A, B, D(C is a waitAll)
    communication can be done after block B


In main:
    In 2-dimensions, N=n*n
    l_N = N/np;
    l_n = n/np
    allocate variables
    h=1/(N+1)
    setup X and Y
    setupB

    tolerance stuff
    if(id>0) idleft=id-1;
    else     idleft=MPI_PROC_NULL;
    if(id<np-1) idright=id+1;
    else        idright=MPI_PROC_NULL;

    call cg(); //solution vector is l_x

    print results
*/

