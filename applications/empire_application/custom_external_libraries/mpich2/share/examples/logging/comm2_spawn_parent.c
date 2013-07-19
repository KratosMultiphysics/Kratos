#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[] )
{
    MPI_Comm intercomm;
    char     processor_name[MPI_MAX_PROCESSOR_NAME];
    int      err, errcodes[256], rank, num_procs;
    int      namelen;
    char     str[10] = "none";

    MPI_Init( &argc, &argv );

    MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Get_processor_name( processor_name, &namelen );

    err = MPI_Comm_spawn( "comm2_spawn_child", MPI_ARGV_NULL, num_procs,
                          MPI_INFO_NULL, 0, MPI_COMM_WORLD,
                          &intercomm, errcodes );
    if ( err != MPI_SUCCESS )
        printf( "Error in MPI_Comm_spawn\n" );

    MPI_Send( "Hello", 6, MPI_CHAR, rank, 101, intercomm );
    MPI_Recv( str, 4, MPI_CHAR, rank, 102, intercomm, MPI_STATUS_IGNORE );

    printf( "Parent %d on %s received from child: %s.\n",
            rank, processor_name, str );
    fflush( stdout );

    MPI_Finalize();

    return 0;
}
