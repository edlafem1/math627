#ifndef LLL_H
#define LLL_H

#include <math.h>
#include "lin_alg.h"

/**
Returns 1 if |u[i,j]| <= 1/2 for all 0 <= i < j <= n-1
*/
int size_reduced(double *U, int m, int n);

int LLL_reduced(double *D, double *U, double w, int m, int n);

void LLL(double *B, double *D, double *U, double *M, double w, int m, int n);

#endif