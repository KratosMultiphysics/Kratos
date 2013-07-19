/*
   (C) 2007 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include "mpi.h"
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#define BUFLEN 512
#define NTIMES 50
#define MAX_THREADS 10

/*
    Concurrent send and recv by multiple threads on each process. 
*/
void *thd_sendrecv( void * );
void *thd_sendrecv( void *comm_ptr )
{
    MPI_Comm     comm;
    int         my_rank, num_procs, next, buffer_size, namelen, idx;
    char        buffer[BUFLEN], processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Status  status;

    comm = *(MPI_Comm *) comm_ptr;

    MPI_Comm_size( comm, &num_procs );
    MPI_Comm_rank( comm, &my_rank );
    MPI_Get_processor_name( processor_name, &namelen );

    fprintf( stdout, "Process %d on %s\n", my_rank, processor_name );
    strcpy( buffer, "hello there" );
    buffer_size = strlen(buffer)+1;

    if ( my_rank == num_procs-1 )
        next = 0;
    else
        next = my_rank+1;

    for ( idx = 0; idx < NTIMES; idx++ ) {
        if (my_rank == 0) {
            /*
            printf("%d sending '%s' \n",my_rank,buffer);
            */
            MPI_Send(buffer, buffer_size, MPI_CHAR, next, 99, comm);
            MPI_Send(buffer, buffer_size, MPI_CHAR, MPI_PROC_NULL, 299, comm);
            /*
            printf("%d receiving \n",my_rank);
            */
            MPI_Recv(buffer, BUFLEN, MPI_CHAR, MPI_ANY_SOURCE, 99,
                     comm, &status);
            /*
            printf("%d received '%s' \n",my_rank,buffer);
            */
        }
        else {
            /*
            printf("%d receiving  \n",my_rank);
            */
            MPI_Recv(buffer, BUFLEN, MPI_CHAR, MPI_ANY_SOURCE, 99,
                     comm, &status);
            MPI_Recv(buffer, BUFLEN, MPI_CHAR, MPI_PROC_NULL, 299,
                     comm, &status);
            /*
            printf("%d received '%s' \n",my_rank,buffer);
            */
            MPI_Send(buffer, buffer_size, MPI_CHAR, next, 99, comm);
            /*
            printf("%d sent '%s' \n",my_rank,buffer);
            */
        }
        /* MPI_Barrier(comm); */
    }

    return 0;
}



int main( int argc,char *argv[] )
{
    MPI_Comm   comm[ MAX_THREADS ];
    pthread_t  thd_id[ MAX_THREADS ];
    int        my_rank, ii, provided;
    int        num_threads;

    MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
    if ( provided != MPI_THREAD_MULTIPLE ) {
        fprintf( stderr, "Aborting, MPI_THREAD_MULTIPLE is needed...\n" );
        fflush( stderr );
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    if ( my_rank == 0 ) {
        if ( argc != 2 ) {
            fprintf( stderr, "Error: %s num_threads\n", argv[0] );
            fflush( stderr );
            MPI_Abort( MPI_COMM_WORLD, 1 );
        }
        num_threads = atoi( argv[1] );
        if ( num_threads < 1 ) {
            fprintf( stderr, "Error: Input num_threads=%d < 1 \n",
                             num_threads );
            fflush( stderr );
            MPI_Abort( MPI_COMM_WORLD, 1 );
        }
        if ( num_threads > MAX_THREADS ) {
            fprintf( stderr, "Error: Input num_threads=%d < %d \n",
                             num_threads, MAX_THREADS );
            fflush( stderr );
            MPI_Abort( MPI_COMM_WORLD, 1 );
        }
        MPI_Bcast( &num_threads, 1, MPI_INT, 0, MPI_COMM_WORLD );
    }
    else
        MPI_Bcast( &num_threads, 1, MPI_INT, 0, MPI_COMM_WORLD );

    MPI_Barrier( MPI_COMM_WORLD );

    for ( ii=0; ii < num_threads; ii++ ) {
        MPI_Comm_dup( MPI_COMM_WORLD, &comm[ii] );
    }

    for ( ii=1; ii < num_threads; ii++ ) {
        pthread_create( &thd_id[ii], NULL, thd_sendrecv, (void *) &comm[ii] );
    }
    thd_sendrecv( (void *) &comm[0] );

    for ( ii=1; ii < num_threads; ii++ )
        pthread_join( thd_id[ii], NULL );

    MPI_Finalize();
    return 0;
}
