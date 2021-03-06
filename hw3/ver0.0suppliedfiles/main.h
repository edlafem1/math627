#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utilities.h"
#ifdef PARALLEL
#include <mpi.h>
#endif

#include "memory.h"

void setup_example (double *l_A, int n, int l_n, int id, int np);

#endif

