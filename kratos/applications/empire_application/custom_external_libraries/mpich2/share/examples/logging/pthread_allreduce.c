/*
   (C) 2007 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include "mpi.h"
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#define SIZE 512
#define NTIMES 50
#define MAX_THREADS 10

/*
    Measures the time taken by concurrent calls to send and receive
    by multiple threads on a node. 
*/
void *thd_allreduce( void * );
void *thd_allreduce( void *comm_ptr )
{
    MPI_Comm  comm;
    int      *inbuf, *outbuf;
    double    stime, etime;
    int       ii;

    inbuf = (int *) malloc( SIZE * sizeof(int) );
    if ( inbuf == NULL ) {
        printf( "Cannot allocate buffer\n" );
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    outbuf = (int *) malloc( SIZE * sizeof(int) );
    if ( outbuf == NULL ) {
        printf( "Cannot allocate buffer\n" );
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    comm = *(MPI_Comm *) comm_ptr;

    stime = MPI_Wtime();
    for ( ii = 0; ii < NTIMES; ii++ ) {
        MPI_Allreduce( inbuf, outbuf, SIZE, MPI_INT, MPI_MAX, comm );
    }
    etime = MPI_Wtime();

    printf( "Time = %f ms\n", ((etime-stime)*1000)/NTIMES );

    free(inbuf);
    free(outbuf);

    pthread_exit( NULL );
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
        pthread_create( &thd_id[ii], NULL, thd_allreduce, (void *) &comm[ii] );
    }

    for ( ii=0; ii < num_threads; ii++ )
        pthread_join( thd_id[ii], NULL );

    MPI_Finalize();
    return 0;
}
