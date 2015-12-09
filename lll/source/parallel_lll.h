#ifndef PARALLEL_LLL_H
#define PARALLEL_LLL_H

#ifdef MPI_INCLUDE
#include <mpi.h>
#endif

#include "lin_alg.h"
#include "delayed_lll.h"

void parallel_LLL(double *B, double *D, double *U, double *M, double w, int m, int n, int id, int np);

#endif