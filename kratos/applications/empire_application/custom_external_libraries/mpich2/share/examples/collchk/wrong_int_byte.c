/*
   (C) 2004 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[] )
{
    int rank, size;
    int ibuff = 0x0;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if ( rank == size-1 )
        /* Create pathological case */
        MPI_Bcast( &ibuff, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD );
    else
        MPI_Bcast( &ibuff, 1, MPI_INT, 0, MPI_COMM_WORLD );

    MPI_Finalize();
    return 0;
}
