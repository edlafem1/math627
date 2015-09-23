#ifndef MEMORY_H
#define MEMORY_H

#include <stdio.h>
#include <stdlib.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

int *allocate_int_vector (int n);
double *allocate_double_vector (int n);
void free_vector (void *x);

#endif

