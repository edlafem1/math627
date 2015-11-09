#ifndef AX_H
#define AX_H

#include "main.h"

void Ax(double *l_v, double *l_u, int l_n, int l_N, int N,
    int id, int idleft, int idright, int np, MPI_Comm comm,
    double *gl, double *gr);

#endif
