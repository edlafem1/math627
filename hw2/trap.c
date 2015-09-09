/* trap.c -- Parallel Trapezoidal Rule, first version
 *
 * Input: None.
 * Output:  Estimate of the integral from a to b of f(x)
 *    using the trapezoidal rule and n trapezoids.
 *
 * Algorithm:
 *    1.  Each process calculates "its" interval of
 *        integration.
 *    2.  Each process estimates the integral of f(x)
 *        over its interval using the trapezoidal rule.
 *    3a. Each process != 0 sends its integral to 0.
 *    3b. Process 0 sums the calculations received from
 *        the individual processes and prints the result.
 *
 * Notes:  
 *    1.  f(x), a, b, and n are all hardwired.
 *    2.  The number of processes (p) should evenly divide
 *        the number of trapezoids (n = 1024)
 *
 * See Chap. 4, pp. 56 & ff. in PPMPI.
 */
#include <stdio.h>

/* We'll be using MPI routines, definitions, etc. */
#include "mpi.h"


main(int argc, char** argv) {
    int         id;   /* My process rank           */
    int         np;         /* The number of processes   */
    double       a = 0.0;   /* Left endpoint             */
    double       b = 1.0;   /* Right endpoint            */
    int         n = 1024;  /* Number of trapezoids      */
    double       h;         /* Trapezoid base length     */
    double       local_a;   /* Left endpoint my process  */
    double       local_b;   /* Right endpoint my process */
    int         local_n;   /* Number of trapezoids for  */
                           /* my calculation            */
    double       approximation;  /* approximation of integral over my interval */
    double       total;     /* Total approximation            */
    int         source;    /* Process sending approximation  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;

		
    double Trap(double local_a, double local_b, int local_n,
              double h);    /* Calculate local approximation of integral  */

    /* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &np);
	/*
	if (argc < 4 && id == 0) {
		printf("Not enough arguments supplied.\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	} else {
		a = atoi(argv[1]);
		b = atoi(argv[2]);
		n = (int)atof(argv[3]);
	}
	*/
    h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/np;  /* So is the number of trapezoids */

    /* Length of each process' interval of
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + id*local_n*h;
    local_b = local_a + local_n*h;
    approximation = Trap(local_a, local_b, local_n, h);

    /* Add up the approximations calculated by each process */
    if (id == 0) {
        total = approximation;
        for (source = 1; source < np; source++) {
            MPI_Recv(&approximation, 1, MPI_DOUBLE, source, tag,
                MPI_COMM_WORLD, &status);
            total = total + approximation;
        }
    } else {  
        MPI_Send(&approximation, 1, MPI_DOUBLE, dest,
            tag, MPI_COMM_WORLD);
    }

    /* Print the result */
    if (id == 0) {
        printf("With:\n");
		printf("\tNumber of processes np = %i\n", np);
		printf("\tn = %i trapezoids,\n", n);
		printf("\th = %24.16e\n", h);
		printf("\th^2 = %24.16e\n", h*h);
		printf("\tTrue value = %24.16e\n", (1.0/3.0)*(b*b*b-a*a*a));
		printf("\tTrue error = %24.16e\n", total-(1.0/3.0)*(b*b*b-a*a*a));
		printf("\tApproximation = %24.16e\n", total);
    }

    /* Shut down MPI */
    MPI_Finalize();
} /*  main  */


double Trap(
          double  local_a   /* in */,
          double  local_b   /* in */,
          int    local_n   /* in */,
          double  h         /* in */) {

    double approximation;   /* Store result in approximation  */
    double x;
    int i;

    double f(double x); /* function we're integrating */

    approximation = (f(local_a) + f(local_b))/2.0;
    x = local_a;
    for (i = 1; i <= local_n-1; i++) {
        x = x + h;
        approximation = approximation + f(x);
    }
    approximation = approximation*h;
    return approximation;
} /*  Trap  */


double f(double x) {
    double return_val;
    /* Calculate f(x). */
    /* Store calculation in return_val. */
    return_val = x*x;
    return return_val;
} /* f */


