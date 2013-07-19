/*
   (C) 2004 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[] )
{
    int    rank, size;
    double dbuff = 0x0;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if ( rank != size-1 ) {
        /* create pathological case */
        MPI_Datatype types[2] = { MPI_INT, MPI_FLOAT };
        int          blks[2]  = { 1, 1};
        MPI_Aint     displs[2] = {0, sizeof(float) };
        MPI_Datatype flt_int_type;
        MPI_Type_struct( 2, blks, displs, types, &flt_int_type );
        MPI_Type_commit( &flt_int_type );
        MPI_Bcast( &dbuff, 1, flt_int_type, 0, MPI_COMM_WORLD );
        MPI_Type_free( &flt_int_type );
    }
    else
        MPI_Bcast( &dbuff, 1, MPI_FLOAT_INT, 0, MPI_COMM_WORLD );

    MPI_Finalize();
    return 0;
}
