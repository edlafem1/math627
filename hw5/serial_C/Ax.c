#include "Ax.h"


void Ax(double *l_v, double *l_u, int l_n, int l_N, int N,
        int id, int idleft, int idright, int np, MPI_Comm comm,
        double *gl, double *gr) {

    int i, l_j, j;
    double temp;
    MPI_Status stat;


    if (id % 2 == 0) {
        MPI_Recv(gl, N, MPI_DOUBLE, idleft, 1, comm, &stat);
        MPI_Recv(gr, N, MPI_DOUBLE, idright, 2, comm, &stat);
        MPI_Send(&(l_u[N*l_N - N]), N, MPI_DOUBLE, idright, 1, comm);
        MPI_Send(l_u, N, MPI_DOUBLE, idleft, 2, comm);
    }
    else {
        MPI_Send(&(l_u[N*l_N - N]), N, MPI_DOUBLE, idright, 1, comm);
        MPI_Send(l_u, N, MPI_DOUBLE, idleft, 2, comm);
        MPI_Recv(gl, N, MPI_DOUBLE, idleft, 1, comm, &stat);
        MPI_Recv(gr, N, MPI_DOUBLE, idright, 2, comm, &stat);
    }

    /*------------------Block A------------------*/
    /*for(l_j=0; l_j<l_N;l_j++){
    for(i=0; i<N; i++){
    temp=4.0*l_u[i+N*l_j];

    if(l_j>0)temp -= l_u[i+N*(l_j-1)];
    if(i>0)temp -=l_u[i-1+N*l_j];
    if(i<N-1)temp -= l_u[i+1+N*l_j];
    if(l_j<N-1)temp -= l_u[i+N*(l_j+1)];
    l_v[i+N*l_j]=temp;
    }
    }*/

    /*------------------Block B------------------*/
    for (l_j = 1; l_j<l_N - 1; l_j++) {
        for (i = 0; i<N; i++) {
            temp = 4.0*l_u[i + N*l_j];

            temp -= l_u[i + N*(l_j - 1)];
            if (i>0)temp -= l_u[i - 1 + N*l_j];
            if (i<N - 1)temp -= l_u[i + 1 + N*l_j];
            temp -= l_u[i + N*(l_j + 1)];
            l_v[i + N*l_j] = temp;
        }
    }

    /* if(id%2==0){
    MPI_Recv(gl, N, MPI_DOUBLE, idleft, 1, comm, &stat);
    MPI_Recv(gr, N, MPI_DOUBLE, idright, 2, comm, &stat);
    MPI_Send(&(l_u[N*l_N-N]), N, MPI_DOUBLE, idright, 1, comm);
    MPI_Send(l_u, N, MPI_DOUBLE, idleft, 2, comm);
    }else{
    MPI_Send(&(l_u[N*l_N-N]), N, MPI_DOUBLE, idright, 1, comm);
    MPI_Send(l_u, N, MPI_DOUBLE, idleft, 2, comm);
    MPI_Recv(gl, N, MPI_DOUBLE, idleft, 1, comm, &stat);
    MPI_Recv(gr, N, MPI_DOUBLE, idright, 2, comm, &stat);
    }*/

    /*------------------Block C------------------*/
    //MPI_Waitall for parallel should go here

    /*------------------Block D------------------*/
    l_j = 0;
    j = l_j + id*l_N;

    for (i = 0; i<N; i++) {
        temp = 4.0*l_u[i + N*l_j];

        if (j>0)temp -= gl[i];
        if (i>0)temp -= l_u[i - 1 + N*l_j];
        if (i<N - 1)temp -= l_u[i + 1 + N*l_j];
        temp -= l_u[i + N*(l_j + 1)];
        l_v[i + N*l_j] = temp;
    }

    l_j = l_N - 1;
    j = l_j + id*l_N;

    for (i = 0; i<N; i++) {
        temp = 4.0*l_u[i + N*l_j];

        temp -= l_u[i + N*(l_j - 1)];
        if (i>0)temp -= l_u[i - 1 + N*l_j];
        if (i<N - 1)temp -= l_u[i + 1 + N*l_j];
        if (j<N - 1)temp -= gr[i];
        l_v[i + N*l_j] = temp;
    }































    return;

    ///////////////////////////////////////////////////////////////////////////////////
    MPI_Status statuses[4];
    MPI_Request requests[4];
    /*---------------Communication---------------*/
    // when using non-blocking comms, need to have a waitall after
    if (id % 2 == 0) {
        // 1 3 2 4
        MPI_Recv(gl, N, MPI_DOUBLE, idleft, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&(l_u[N*l_N - N]), N, MPI_DOUBLE, idright, 1, MPI_COMM_WORLD);
        MPI_Recv(gr, N, MPI_DOUBLE, idright, 2, MPI_COMM_WORLD, &status);
        MPI_Send(l_u, N, MPI_DOUBLE, idleft, 2, MPI_COMM_WORLD);
    }
    else {
        // 3 1 4 2
        MPI_Recv(gr, N, MPI_DOUBLE, idright, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(gl, N, MPI_DOUBLE, idleft, 1, MPI_COMM_WORLD, &status);
        MPI_Send(l_u, N, MPI_DOUBLE, idleft, 2, MPI_COMM_WORLD);
        MPI_Send(&(l_u[N*l_N - N]), N, MPI_DOUBLE, idright, 1, MPI_COMM_WORLD);
    }
    /*------------------Block A------------------*/
    // this is serial code
    for (l_j = 0; l_j < l_N; l_j++) {
        for (i = 0; i < N; i++) {
            temp = 4.0 * l_u[i + N*l_j];

            if (l_j > 0)     temp -= l_u[i + N*(l_j-1)];
            if (i > 0)       temp -= l_u[i - 1 + N*l_j];
            if (i < N - 1)   temp -= l_u[i + 1 + N*l_j];
            if (l_j < N - 1) temp -= l_u[i + N*(l_j + 1)];
            
            l_v[i + N*l_j] = temp;
        }
    }

    /*------------------Block B------------------*/
    for (l_j = 1; l_j < l_N-1; l_j++) {
        for (i = 0; i < N; i++) {
            temp = 4.0 * l_u[i + N*l_j];

                             temp -= l_u[i + N*(l_j - 1)];
            if (i > 0)       temp -= l_u[i - 1 + N*l_j];
            if (i < N - 1)   temp -= l_u[i + 1 + N*l_j];
                             temp -= l_u[i + N*(l_j + 1)];

            l_v[i + N*l_j] = temp;
        }
    }
    
    /*------------------Block D------------------*/
    l_j = 0; 
    j = l_j + id*l_N;
    for (i = 0; i < N; i++) {
        temp = 4.0 * l_u[i + N*l_j];

        if (j > 0)       temp -= gl[i];
        if (i > 0)       temp -= l_u[i - 1 + N*l_j];
        if (i < N - 1)   temp -= l_u[i + 1 + N*l_j];
                         temp -= l_u[i + N*(l_j + 1)];

        l_v[i + N*l_j] = temp;
    }
    /*---------------Block D part 2--------------*/
    l_j = l_N - 1;
    j = l_j + id*l_N;
    for (i = 0; i < N; i++) {
        temp = 4.0 * l_u[i + N*l_j];

                            temp -= l_u[i + N*(l_j - 1)];
        if (i > 0)       temp -= l_u[i - 1 + N*l_j];
        if (i < N - 1)   temp -= l_u[i + 1 + N*l_j];
        if (j < N - 1)   temp -= gr[i];

        l_v[i + N*l_j] = temp;
    }
    
}