#ifndef UTILITIES_H
#define UTILITIES_H

#include <math.h>

#include "main.h"

void setupB(double *l_r, double *x, double *y, int l_N, int N, double h, int id);

double f(double X, double Y);

double parallel_dot(double *l_x, double *l_y, int l_n, MPI_Comm comm);

double *allocate_double_vector(int n);

void free_vector(void *x);

#endif
