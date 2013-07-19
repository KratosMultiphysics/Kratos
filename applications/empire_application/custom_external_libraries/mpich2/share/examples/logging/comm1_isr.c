/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main( int argc, char *argv[] )
{
    int          comm_rank, comm_size, comm_neighbor;
    int          world_rank, world_size, world_neighbor;
    int          icolor, namelen, ibuffer;
    char         processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm     splited_comm, duped_comm, inter_comm, *comm_ptr;
    MPI_Request  world_request, comm_request;
    MPI_Status   world_status, comm_status;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
    MPI_Get_processor_name( processor_name, &namelen );

    fprintf( stdout, "world_rank %d on %s\n", world_rank, processor_name );
    fflush( stdout );

    if ( world_rank == world_size - 1 )
        world_neighbor = 0;
    else
        world_neighbor = world_rank + 1;

    MPI_Irecv( &ibuffer, 1, MPI_INT, MPI_ANY_SOURCE,
               99, MPI_COMM_WORLD, &world_request );
    MPI_Send( &world_rank, 1, MPI_INT, world_neighbor,
               99, MPI_COMM_WORLD );
    MPI_Wait( &world_request, &world_status );

    /* Split all processes into 2 separate intracommunicators */
    icolor = world_rank % 2;
    MPI_Comm_split( MPI_COMM_WORLD, icolor, world_rank, &splited_comm );

    /* Put in a Comm_dup so local comm ID is different in 2 splited comm */
    if ( icolor == 0 ) {
        MPI_Comm_dup( splited_comm, &duped_comm );
        comm_ptr  = &duped_comm;
    }
    else
        comm_ptr  = &splited_comm;

    MPI_Comm_size( *comm_ptr, &comm_size );
    MPI_Comm_rank( *comm_ptr, &comm_rank );

    if ( comm_rank == comm_size - 1 )
        comm_neighbor = 0;
    else
        comm_neighbor = comm_rank + 1;

    MPI_Irecv( &ibuffer, 1, MPI_INT, MPI_ANY_SOURCE,
               999, *comm_ptr, &comm_request );
    MPI_Send( &comm_rank, 1, MPI_INT, comm_neighbor, 999, *comm_ptr );
    MPI_Wait( &comm_request, &comm_status );

    /* Form an intercomm between the 2 splited intracomm's */
    if ( icolor == 0 )
        MPI_Intercomm_create( *comm_ptr, 0, MPI_COMM_WORLD, 1,
                              9090, &inter_comm );
    else
        MPI_Intercomm_create( *comm_ptr, 0, MPI_COMM_WORLD, 0,
                              9090, &inter_comm );

    if ( comm_rank == 0 ) {
        MPI_Irecv( &ibuffer, 1, MPI_INT, 0,
                   9999, inter_comm, &comm_request );
        MPI_Send( &comm_rank, 1, MPI_INT, 0, 9999, inter_comm );
        MPI_Wait( &comm_request, &comm_status );
    }

    /* Free all communicators created */
    MPI_Comm_free( &inter_comm );
    if ( icolor == 0 )
        MPI_Comm_free( &duped_comm );
    MPI_Comm_free( &splited_comm );

    MPI_Finalize();
    return( 0 );
}
