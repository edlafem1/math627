#include "utilities.h"
#define M_PI 3.14159265358979323846

void setupB(double *l_r, double *x, double *y, int l_N, int N, double h, int id) {
    // copy code from Gobbert notes #18
    // l_r is where we store B
    int l_j, j, i;
    for(l_j = 0; l_j < l_N; l_j++) {
        j = l_j + id*l_N; // <<<<<<<<<<<<<<< WHY IS THIS NOT USED? Ahh, see below
        for (i = 0; i < N; i++) {
            l_r[i + N*l_j] = (h*h)*f(x[i], y[j]); // y[i] changed to y[j];
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
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return x;
}

void free_vector(void *x)
{
    free(x);
}

double F(double x, double y) {
    // sin^2(pi*x)sin^2(pi*y)
    return (sin(M_PI*x)*sin(M_PI*x))*(sin(M_PI*y)*sin(M_PI*y));
}

double fd_norm(double *l_x, int h, int N) {
    double diff, max = 0;
    for (int i = 1; i <= N; ++i) { // col(i.e. x axis)
        for (int j = 1; j <= N; ++j) { // row(i.e. y axis)
            diff = F(i*h, j*h) - l_x[N*i + j];
            if (diff < 0) diff *= -1;
            if (diff > max) max = diff;
        }
    }

    return max;
}

