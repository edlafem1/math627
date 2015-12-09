#include "memory.h"

/* 09/12/02-10/10/02, updated 02/07/08 by Matthias K. Gobbert */

int *allocate_int_vector (int n)
{
  int *x;

  x = (int*) calloc (n, sizeof(int));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef MPI_INCLUDE
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}

double *allocate_double_vector (int n)
{
  double *x;

  x = (double*) calloc (n, sizeof(double));

  if (x == NULL)
  {
    fprintf (stderr, "Problem allocating memory for vector\n");
#ifdef MPI_INCLUDE
    MPI_Abort (MPI_COMM_WORLD, 1);
#else
    exit (1);
#endif
  }

  return x;
}

void free_vector (void *x)
{
  free (x);
}

