/*
   (C) 2004 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[] )
{
    int rank;
    int ibuff = 0x0;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    /* root is changing with rank */
    MPI_Bcast( &ibuff, 1, MPI_INT, rank, MPI_COMM_WORLD );

    MPI_Finalize();
    return 0;
}
