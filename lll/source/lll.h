#ifndef LLL_H
#define LLL_H

#include <math.h>

/**
Returns 1 if |u[i,j]| <= 1/2 for all 0 <= i < j <= n-1
*/
int size_reduced(double *U, int m, int n);

#endif