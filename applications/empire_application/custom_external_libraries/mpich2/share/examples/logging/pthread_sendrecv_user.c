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

#include "mpe_log.h"
#include "mpe_log_thread.h"
extern CLOG_CommSet_t  *CLOG_CommSet;

typedef struct {
    int    start_evtID;
    int    final_evtID;
    char  *name;
    char  *color;
} USER_State;

#define  USER_MAX_STATES   2

static USER_State user_states[ USER_MAX_STATES ];

#define USER_SENDRECV_ID  0
#define USER_RECVSEND_ID  1

void USER_Init_log( void );
void USER_Init_log( void )
{
    USER_State *state;
    int         idx;
    int         myrank;

    /* Define each state's legend name and color */
    state = &user_states[ USER_SENDRECV_ID ];
    state->name   = "USER_SendRecv";
    state->color  = "MistyRose";

    state = &user_states[ USER_RECVSEND_ID ];
    state->name   = "USER_RecvSend";
    state->color  = "HotPink";

    PMPI_Comm_rank( MPI_COMM_WORLD, &myrank );

    for ( idx = 0; idx < USER_MAX_STATES; idx++ ) {
        state = &user_states[ idx ];
        MPE_Log_get_state_eventIDs( &(state->start_evtID),
                                    &(state->final_evtID) );
        if ( myrank == 0 ) {
            MPE_Describe_comm_state( MPI_COMM_WORLD,
                                     state->start_evtID, state->final_evtID,
                                     state->name, state->color, NULL );
        }
    }
}

/* Define macro to do communicator-thread user-defined logging */
#define USER_LOG_STATE_DECL          \
          USER_State      *state;

#define USER_LOG_STATE_BEGIN(name,comm) \
    state   = &user_states[ name ]; \
    MPE_Log_comm_event( comm, state->start_evtID, NULL );

#define USER_LOG_STATE_END \
    MPE_Log_comm_event( comm, state->final_evtID, NULL );


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

    /* Declare variables for user-defined communicator-thread logging */
    USER_LOG_STATE_DECL

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
            USER_LOG_STATE_BEGIN( USER_SENDRECV_ID, comm )
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
            USER_LOG_STATE_END
        }
        else {
            USER_LOG_STATE_BEGIN( USER_RECVSEND_ID, comm )
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
            USER_LOG_STATE_END
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

    /* Define/describe MPE user-defined states */
    USER_Init_log();

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
