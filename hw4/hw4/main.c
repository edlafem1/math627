#include "main.h"

int main(int argc, char *argv[])
{
	int id, np, processor_name_len;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double *l_x, *l_y, *temp_nvector, *y;

	double *l_A, *l_B, *l_C, *l_D;


	MPI_Init(&argc, &argv);
	/* Check processes: */
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Get_processor_name(processor_name, &processor_name_len);

	/* Do stuff here */
	int m, k, n;
	if (argc != 4) {
		printf("Usage: %s [m] [k] [n]. Using default values 3, 4, 3.\n", argv[0]);
		m = 3;
		k = 4;
		n = 3;
	}
	else {
		m = (int)atof(argv[1]);
		k = (int)atof(argv[2]);
		n = (int)atof(argv[3]);
	}
	if (id == 0) {
		double *A = allocate_double_vector(m*k);
		double *B = allocate_double_vector(k*n);
		double *C = allocate_double_vector(m*n);
		setABC_example(A, B, C, m, k, n);

		printf("Matrix A:\n");
		print_Matrix(A, m, k, id, np);
		printf("Matrix B:\n");
		print_Matrix(B, k, n, id, np);
		printf("Matrix C:\n");
		print_Matrix(C, m, n, id, np);

	}


	MPI_Finalize();
}