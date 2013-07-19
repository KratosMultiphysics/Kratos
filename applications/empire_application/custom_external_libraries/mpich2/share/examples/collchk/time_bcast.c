/*
   (C) 2004 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define COUNT 131072

int main( int argc, char *argv[] )
{
    double *dbuff;
    double  time_init, time_final;
    int     rank;
    int     count, num_itr, idx;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if ( argc > 1 && argv[1] != NULL )
        count = atoi( argv[1] );
    else
        count = 1;

    if ( argc > 2 && argv[2] != NULL )
        num_itr = atoi( argv[2] );
    else
        num_itr = 1;

    dbuff = (double*) malloc( count * sizeof(double) );
    if ( dbuff == NULL ) {
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );
    time_init   = MPI_Wtime();

    for ( idx = 0; idx < num_itr; idx++ ) {
        MPI_Bcast( dbuff, count, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    }

    /* MPI_Barrier( MPI_COMM_WORLD ); */
    time_final  = MPI_Wtime();

    fprintf( stdout, "time taken by %dX%d MPI_Bcast() at rank %d = %f\n",
                      count, num_itr, rank, time_final - time_init );

    free( dbuff );
    MPI_Finalize();
    return 0;
}
