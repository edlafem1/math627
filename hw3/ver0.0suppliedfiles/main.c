#include "main.h"

/* 09/12/02-10/10/02, updated 02/07/08 by Matthias K. Gobbert */

void setup_example (double *l_A, int n, int l_n, int id, int np)
{
  int i, j, l_j;

  /* example in global notation:
   *   A(i,j) = 1 / (1 + i + j)
   * for 0<=i<n and 0<=j<n   ==> A is n x n matrix
   */
  for (l_j = 0; l_j < l_n; l_j++)
    for (i = 0; i < n; i++)
    {
      j = l_j + id*l_n;
      l_A[i+l_j*n] = 1.0 / ((double) (1 + i + j));
    }
  /* Note 1: loops ordered to access memory of l_A contiguously */
  /* Note 2: systematic choice of variable names and uses:
   * 0<=l_j<l_n is the index into the local l_A,
   * j is the mathematical index into the global A (not actually set up),
   * hence, j = l_j + id*l_n transforms l_j to j such that we get
   * id*l_n<=j<(id+1)*l_n on Process id.
   */
}

int main (int argc, char *argv[])
{
  int id, np, processor_name_len;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int n;
  double double_n;
  int l_n, l_i;
  double *l_A, *l_x, *l_y;
  double tol;
  int itmax;

  MPI_Init(&argc, &argv);

  /* Check processes: */
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Get_processor_name(processor_name, &processor_name_len);
  /*
  sprintf(message, "Hello from %3.3d of %3.3d on %s\n", id, np, processor_name);
  if (id == 0) {
    printf("%s\n", message);
    for (i = 1; i < np; i++) {
      MPI_Recv(message, 100, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      printf("%s\n", message);
    }
  } else {
    comm = MPI_COMM_WORLD;
    MPI_Send(message, 1+strlen(message), MPI_CHAR, 0, 0, comm);
  }
  */

  /* test output of command-line arguments: */
  /**/
  if (id == 0) {
    printf("argc = %d\n", argc);
    for (int i = 0; i < argc; i++) {
      printf("%s\n", argv[i]);
    }
  }
  /**/

  /* process command-line inputs: */
  if (argc != 4)
  {
    if (id == 0)
    {
      printf("Usage: \"./power n tol itmax\" \n");
      printf("  with int n, double tol, and int itmax\n");
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  n     =      atoi(argv[1]);
  tol   =      atof(argv[2]);
  itmax = (int)atof(argv[3]);

  /* number of processes np must divide n: */
  if ( (n % np) != 0)
  {
    if (id == 0)
    {
      printf("Error: np must divide n!\n");
      printf("  n = %d, np = %d, n%%np = %d\n", n, np, (n%np));
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* compute size of local blocks: */
  l_n = n / np;
  if (id == 0)
  {
    printf("n = %d, np = %d, l_n = %d\n", n, np, l_n);
    printf("\n");
    fflush(stdout);
  }

  /* allocate storage: */
  l_A = allocate_double_vector((n*l_n)); /* l_A is n x l_n matrix */
  l_x = allocate_double_vector((l_n));   /* l_x is l_n-vector */
  l_y = allocate_double_vector((l_n));   /* l_y is l_n-vector */

  /* setup example: */
  setup_example(l_A, n, l_n, id, np);
  for (l_i = 0; l_i < l_n; l_i++) {
    l_x[l_i] = 1.0 / sqrt(((double)(n)));
  }

  double *A;
  int destination = 0;
  if (id == 0) {
	A = (double *) calloc(n*n, sizeof(double));
  }
  
  MPI_Gather(l_A, n*l_n, MPI_DOUBLE, A, n*l_n, MPI_DOUBLE, destination, MPI_COMM_WORLD);
  
  /* PRINT STUFF HERE */
  if (id == 0) {
	  int i, j;
	  for (i = 0; i < n; i++) {
		  for (j = 0; j < n; j++) {
			  printf("%- 09.9f   ", A[i*n + j]);
		  }
		  printf("\n");
	  }
  }






  free_vector(l_y);
  free_vector(l_x);
  free_vector(l_A);
  free(A);
  MPI_Finalize();

  return 0;
}
