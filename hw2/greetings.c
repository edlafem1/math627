/* greetings.c -- greetings program
 *
 * Send a message from all processes with rank != 0 to process 0.
 * Process 0 prints the messages received.
 *
 * Input: none.
 * Output: contents of messages received by process 0.
 *
 * See Chapter 3, pp. 41 & ff in PPMPI.
 */
#include <stdio.h> 
#include <string.h> 
#include "mpi.h"

main(int argc, char* argv[]) {
    int id; /* rank of process */
    int np; /* number of processes */
    int source; /* rank of sender */
    int dest; /* rank of receiver */
    int tag = 0; /* tag for messages */
    char message[100]; /* storage for message */
    MPI_Status status; /* return status for */
                               /* receive */
    /* Start up MPI */
    MPI_Init(&argc, &argv);
    /* Find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (id == 0) {
        /* Create message */
        sprintf(message, "Greetings from process %d!",
            id);
        dest = (id + 1) % np;
        /* Use strlen+1 so that '\0' gets transmitted */
        MPI_Send(message, strlen(message)+1, MPI_CHAR,
            dest, tag, MPI_COMM_WORLD);
			
		source = (id - 1 + np) % np;
        MPI_Recv(message, 100, MPI_CHAR, source, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
        printf("My ID: %d, From: %i, Message: %s\n", id, status.MPI_SOURCE, message);
    } else {
		source = (id - 1 + np) % np;
        MPI_Recv(message, 100, MPI_CHAR, source, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
        printf("My ID: %d, From: %i, Message: %s\n", id, status.MPI_SOURCE, message);
        
		/* Create message */
        sprintf(message, "Greetings from process %d!",
            id);
        dest = (id + 1) % np;
        /* Use strlen+1 so that '\0' gets transmitted */
        MPI_Send(message, strlen(message)+1, MPI_CHAR,
            dest, tag, MPI_COMM_WORLD);
    }
    /* Shut down MPI */
    MPI_Finalize();
} /* main */
