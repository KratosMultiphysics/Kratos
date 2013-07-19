/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include "mpi.h"
#include <stdio.h>
#include <string.h>

#include "mpe.h"

#define BUFLEN 512

int main(int argc, char *argv[])
{
    int myid, numprocs, next, namelen, idx;
    char buffer[BUFLEN], processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);

    fprintf(stderr,"Process %d on %s\n", myid, processor_name);
    strcpy(buffer,"hello there");
    if (myid == numprocs-1)
	next = 0;
    else
	next = myid+1;

    for (idx = 0; idx < 1000; idx++) {
        if (myid == 0) {
            /*
   	    printf("%d sending '%s' \n",myid,buffer);
            */
	    MPI_Send(buffer, strlen(buffer)+1, MPI_CHAR, next, 99,
                     MPI_COMM_WORLD);
	    MPI_Send(buffer, strlen(buffer)+1, MPI_CHAR, MPI_PROC_NULL, 299,
                     MPI_COMM_WORLD);
            /*
	    printf("%d receiving \n",myid);
            */
	    MPI_Recv(buffer, BUFLEN, MPI_CHAR, MPI_ANY_SOURCE, 99,
                     MPI_COMM_WORLD, &status);
            /*
	    printf("%d received '%s' \n",myid,buffer);
            */
        }
        else {
            /*
 	    printf("%d receiving  \n",myid);
            */
	    MPI_Recv(buffer, BUFLEN, MPI_CHAR, MPI_ANY_SOURCE, 99,
                     MPI_COMM_WORLD, &status);
	    MPI_Recv(buffer, BUFLEN, MPI_CHAR, MPI_PROC_NULL, 299,
                     MPI_COMM_WORLD, &status);
            /*
	    printf("%d received '%s' \n",myid,buffer);
            */
	    MPI_Send(buffer, strlen(buffer)+1, MPI_CHAR, next, 99,
                     MPI_COMM_WORLD);
            /*
	    printf("%d sent '%s' \n",myid,buffer);
            */
        }
        /* MPI_Barrier(MPI_COMM_WORLD); */
        if ( idx % 100 == 0 )
            MPE_Log_sync_clocks();
    }
    MPI_Finalize();
    return (0);
}
