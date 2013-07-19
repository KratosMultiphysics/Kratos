/*
 *
 *  (C) 2003 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

int main( int argc, char *argv[] )
{
    MPI_Comm  comm;
    double    *buf;
    int       rank, size;
    int       num_itr, count, ii, idx; 
    double    time_init, time_final;
    /* int       errs = 0; */

    MPI_Init( &argc, &argv );

        comm = MPI_COMM_WORLD;

        MPI_Comm_size( comm, &size );
        MPI_Comm_rank( comm, &rank );

        if ( argc > 1 && argv[1] != NULL )
            count = atoi( argv[1] );
        else
            count = 1;

        if ( argc > 2 && argv[2] != NULL )
            num_itr = atoi( argv[2] );
        else
            num_itr = 1;

        /* Contiguous data */
        buf = (double *) malloc( count * sizeof(double) );
        for ( ii = 0; ii < count; ii++ )
            buf[ii] = rank + ii;

        MPI_Barrier( comm );
        MPI_Barrier( comm );
        time_init   = MPI_Wtime();

        for ( idx = 0; idx < num_itr; idx++ ) {
          MPI_Allreduce( MPI_IN_PLACE, buf, count, MPI_DOUBLE, MPI_SUM, comm );
        }

        /* MPI_Barrier( comm ); */
        time_final  = MPI_Wtime();

        fprintf( stdout, "time taken by %dx%d MPI_Allreduce() at rank %d = %f\n",
                         count, num_itr, rank, time_final - time_init );

        /* Check the results */
        /*
        for ( ii = 0; ii < count; ii++ ) {
            int result = ii * size + (size*(size-1))/2;
            if ( buf[ii] != result ) {
                errs ++;
                if (errs < 10) {
                    fprintf( stderr, "buf[%d] = %d expected %d\n",
                                     ii, buf[ii], result );
                }
            }
        }
        */
        free( buf );

    MPI_Finalize();
    return 0;
}
