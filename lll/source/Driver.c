#include "Driver.h"

#define DATA_FOLDER "/home/edlafem1/student_user/bases/"

/**
Reads in a column oriented matrix of dimension m x n from file referenced by filename into B.
The data in the file must be arranged such that each group of m successive integers is one 
column of B. Note: every element MUST be an integer.

Returns 0 on success, -1 on error.
*/
int get_Initial_Basis(double *B, int m, int n, char *filename) {
    /*
    To make vectors in Matlab:
    M = randi(r,m,n)
    Creates m rows of linearly independent n-vectors with values x between 1<=x<=r.
    To write to file:
    dlmwrite('filename',M) will be a comma deliminated file of the values

    Mathematica and MuPad print something like this:
    ( a  b  c )
    ( d  e  f )
    ( g  h  i )
    Which is interpreted as list of vectors.
    We think of it like this:
    [ a  d  g ]
    [ b  d  h ]
    [ c  f  i ]
    Which needs to be in the input file looking like the first representation.
    */
    fprintf(stdout, "filename: %s\n", filename);
    fprintf(stdout, "Dimensions are %i x %i\n", m, n);
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file.\nQuiting.\n");
        return -1;
    }
    double j;

    for (int i = 0; i < m*n; i++) {
        fscanf(file, "%lf", &(j));
        //printf("%lf\n", j);
        B[i] = j;
    }
    fclose(file);
    return 0;
}

/**
Memory Estimates assuming m == n:
4n^2 + n
With 64 GB, n=46341 max
*/
int main(int argc, char *argv[]) {
    fprintf(stdout, "Starting point\n");
#ifdef MPI_INCLUDE
    fprintf(stdout, "MPI Included\n");
    MPI_Init(&argc, &argv);
#endif

#ifdef DELAYED_LLL
    fprintf(stdout, "Delayed LLL\n");
#else
    fprintf(stdout, "Normal LLL\n");
#endif

    /**
    m x n is the dimension of the lattice basis.
    */
    int m, n;

    /**
    A constant used as a goodness level.
    0.25 < w <= 1, however polynomial time guaranteed for only 0.25 < w < 1.
    The closer w gets to 1, the better the resulting basis.
    By default, we have this value at 0.75 because it is a common value and that used by
    MuPad (part of MATLAB).
    */
    double w = 0.75;

    char filename[512];
    if (argc >= 3) {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        if (m <= 0 || n <= 0 || n > m) {
            fprintf(stderr, "m and n must satisfy 0 < n <= m.\nQuiting\n");
#ifdef MPI_INCLUDE
            MPI_Abort(MPI_COMM_WORLD, 1);
#else
            exit(1);
#endif
        }
        sprintf(filename, "%s%ix%i.dat", DATA_FOLDER, m, n);
        if (argc == 4) {
            w = atof(argv[3]);
            if (w <= .25 || w >= 1) {
                fprintf(stderr, "w must be > 0.25 and < 1.\nQuiting.\n");
#ifdef MPI_INCLUDE
                MPI_Abort(MPI_COMM_WORLD, 1);
#else
                exit(1);
#endif
            }
        }
    }
    else {
        m = 8;
        n = 8;
        sprintf(filename, "input.txt");
    }
    /**
    B is the initial basis. 
    Dimensions m x n.
    B = [b_1, b_2, ..., b_n] where b_i are m-length vectors.
    */
    double *B = allocate_double_vector(m*n);
    fprintf(stdout, "Alloc B\n");
    /**
    Q is the gram-schmidt orthogonalized basis. 
    Dimensions m x n.
    Q = [b_1*, b_2*, ..., b_n*] where b_i* are orthogonal m-length vectors.
    */
    double *Q = allocate_double_vector(m*n);
    fprintf(stdout, "Alloc Q\n");
    /**
    D is a diagonal matrix with the L2 norm of the gram-schmidt vectors on the main diagonal, zeros elsewhere. 
    Dimension is n x n, but we can represent it with just a vector of length n to save memory.
    Dimensions n x 1.
    D=diag(d_i), d_i = ||b_i*||^2
    */
    double *D = allocate_double_vector(n);
    fprintf(stdout, "Alloc D\n");
    /**
    U is an upper-triangular matrix with ones on the main diagonal. 
    Dimensions n x n.
    */
    double *U = allocate_double_vector(n*n);
    fprintf(stdout, "Alloc U\n");
    /**
    M is a unimodular matrix that relates two bases for the same lattice by C=BM where
    B is the original basis and C is the new basis. Initially, the LLL algorithm forces this to be
    the identity matrix I_n but relies on the assumption that it starts out as filled with zeros.
    Dimensions n x n.
    */
    double *M = allocate_double_vector(n*n);
    fprintf(stdout, "Alloc M\n");


    int basis_initialized = 0;
    //basis_initialized = get_Initial_Basis(B, m, n, filename);
    if (basis_initialized != 0) {
#ifdef MPI_INCLUDE
        MPI_Abort(MPI_COMM_WORLD, 1);
#else
        exit(1);
#endif
    }
    fprintf(stdout, "Got basis\n");

#ifdef DEBUG_LLL
    fprintf(stdout, "Initial Basis:\n");
    printMatrix(B, m, n);
#endif

    //gramschmidt_process(B, Q, m, n);
    fprintf(stdout, "Done GS\n");
#ifdef DEBUG_LLL
    fprintf(stdout, "Q:\n");
    printMatrix(Q, m, n);
#endif   

    //qdu_decomposition(B, Q, D, U, m, n);
    fprintf(stdout, "Done QDR\n");
#ifdef DEBUG_LLL  
    fprintf(stdout, "D:\n");
    printMatrix(D, n, 1);
#endif

    free(Q);
    
#ifdef DEBUG_LLL
    fprintf(stdout, "U:\n");
    printMatrix(U, n, n);
#endif

    //identity(M, n, n, 1);
    fprintf(stdout, "Done identity\n");
#ifdef MPI_INCLUDE
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    fprintf(stdout, "Time started\n");
#endif
#ifdef DELAYED_LLL
   delayed_LLL(B, D, U, M, w, m, n);
#else
    LLL(B, D, U, M, w, m, n);
#endif
    fprintf(stdout, "Done with LLL stuff\n");
#ifdef MPI_INCLUDE
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    fprintf(stdout, "Time stopped\n");
#endif

    if (m * n < 128 * 128) {
        fprintf(stdout, "Final Basis:\n");
        printMatrix(B, m, n);
    }
    else {
        fprintf(stdout, "Final basis too large to print out.\n");
    }
#ifdef DEBUG_LLL

    fprintf(stdout, "D:\n");
    printMatrix(D, m, 1);

    fprintf(stdout, "M: \n");
    printMatrix(M, n, n);

    fprintf(stdout, "U:\n");
    printMatrix(U, n, n);
#endif

    //fprintf(stdout, "Is size reduced? %s\n", (size_reduced(U, m, n)==1) ? "yes" : "no");

    //fprintf(stdout, "Is LLL reduced? %s\n", (LLL_reduced(D, U, w, m, n)==1) ? "yes" : "no");

#ifdef MPI_INCLUDE
    fprintf(stdout, "Elapsed time for LLL algorithm only: %lf\n", (end_time - start_time));
#endif

    free(B);
    free(D);
    free(U);
    free(M);

#ifdef MPI_INCLUDE
    MPI_Finalize();
#endif

    return 0;
}