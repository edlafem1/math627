void setABC_example (double *A, double *B, double *C, int m, int k, int n)
{
  int i, q, j;
  double ri, rj, rk;

  /* set up example:
   * in mathematical counting:
   * i = 1, 2, ..., m, q = 1, 2, ..., k, j = 1, 2, ..., n
   * A = (A(i,q)) m-by-k matrix, A(i,q) = i,
   * B = (B(q,j)) k-by-n matrix, B(q,j) = 10*q + j,
   * then results
   * C = A B = (C(i,j)) m-by-n matrix, C(i,j) = (5*k*k + (5+j)*k) * i
   */
  rk = (double) k;
  for (q = 0; q < k; q++)
    for (i = 0; i < m; i++)
      A[i+q*m] = (double) (i+1);
  for (j = 0; j < n; j++)
    for (q = 0; q < k; q++)
      B[q+j*k] = (double) (10*(q+1) + j+1);
  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
    {
      ri = (double) i;
      rj = (double) j;
      C[i+j*m] = ((5.0*rk*rk + (5.0+rj+1.0)*k) * (ri+1.0));
    }
}

