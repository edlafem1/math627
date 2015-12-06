#ifndef LLL_H
#define LLL_H

#include <math.h>
#include "lin_alg.h"

#define NUM_ERR 1e-14
/**
Returns 1 if |u[i,j]| <= 1/2 for all 0 <= i < j <= n-1
*/
int size_reduced(double *U, int m, int n);

int LLL_reduced(double *D, double *U, double w, int m, int n);

int closest_integer(double x);

void reduce(double *U, double *B, double *M, int i, int j, int m, int n);

void LLL(double *B, double *D, double *U, double *M, double w, int m, int n);

#endif