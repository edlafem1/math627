#ifndef GRAM_SCHMIDT_H
#define GRAM_SCHMIDT_H

void gramschmidt_process(double *A, double *E, int m, int n);

void qdu_decomposition(double *A, double *E, double *D, double *U, int m, int n);

#endif
