#ifndef DELAYED_LLL_H
#define DELAYED_LLL_H

#include "lll.h"
#include <math.h>
#include "lin_alg.h"

void delayed_LLL(double *B, double *D, double *U, double *M, double w, int m, int n);

void reduceSwapRestore(int i, int gamma, double *B, double *D, double *U, double *M, int m, int n);

#endif