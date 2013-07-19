/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mpi.h"
#include "mpe.h"

#define ITER_COUNT  5

double f( double );
double f( double a )
{
    return (4.0 / (1.0 + a*a));
}

int main( int argc, char *argv[] )
{
    int  n, myid, numprocs, ii, jj;
    double PI25DT = 3.141592653589793238462643;
    double mypi, pi, h, sum, x;
    double startwtime = 0.0, endwtime;
    int namelen; 
    int event1a, event1b, event2a, event2b,
        event3a, event3b, event4a, event4b;
    char processor_name[ MPI_MAX_PROCESSOR_NAME ];

    MPE_LOG_BYTES  bytebuf;
    int            bytebuf_pos;


    MPI_Init( &argc, &argv );
        
        MPI_Pcontrol( 0 );

    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &myid );

    MPI_Get_processor_name( processor_name, &namelen );
    fprintf( stderr, "Process %d running on %s\n", myid, processor_name );

    /*
        MPE_Init_log() & MPE_Finish_log() are NOT needed when
        liblmpe.a is linked with this program.  In that case,
        MPI_Init() would have called MPE_Init_log() already.
    */
#if defined( NO_MPI_LOGGING )
    MPE_Init_log();
#endif

    /*  Get event ID from MPE, user should NOT assign event ID directly */
    event1a = MPE_Log_get_event_number(); 
    event1b = MPE_Log_get_event_number(); 
    event2a = MPE_Log_get_event_number(); 
    event2b = MPE_Log_get_event_number(); 
    event3a = MPE_Log_get_event_number(); 
    event3b = MPE_Log_get_event_number(); 
    event4a = MPE_Log_get_event_number(); 
    event4b = MPE_Log_get_event_number(); 

    if ( myid == 0 ) {
        MPE_Describe_state( event1a, event1b, "Broadcast", "red" );
        MPE_Describe_info_state( event2a, event2b, "Sync", "orange",
                                 "source = %s()'s line %d." );
        MPE_Describe_info_state( event3a, event3b, "Compute", "blue",
                                 "mypi = %E computed at iteration %d." );
        MPE_Describe_info_state( event4a, event4b, "Reduce", "green",
                                 "final pi = %E at iteration %d." );
    }

    if ( myid == 0 ) {
        n = 1000000;
        startwtime = MPI_Wtime();
    }
    MPI_Barrier( MPI_COMM_WORLD );

    MPI_Pcontrol( 1 );
    /*
    MPE_Start_log();
    */

    for ( jj = 0; jj < ITER_COUNT; jj++ ) {
        MPE_Log_event( event1a, 0, NULL );
        MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD );
        MPE_Log_event( event1b, 0, NULL );
    
        MPE_Log_event( event2a, 0, NULL );
        MPI_Barrier( MPI_COMM_WORLD );
            int line_num;
            bytebuf_pos = 0;
            MPE_Log_pack( bytebuf, &bytebuf_pos, 's',
                          sizeof(__func__)-1, __func__ );
            line_num = __LINE__;
            MPE_Log_pack( bytebuf, &bytebuf_pos, 'd', 1, &line_num );
        MPE_Log_event( event2b, 0, bytebuf );

        MPE_Log_event( event3a, 0, NULL );
        h   = 1.0 / (double) n;
        sum = 0.0;
        for ( ii = myid + 1; ii <= n; ii += numprocs ) {
            x = h * ((double)ii - 0.5);
            sum += f(x);
        }
        mypi = h * sum;
            bytebuf_pos = 0;
            MPE_Log_pack( bytebuf, &bytebuf_pos, 'E', 1, &mypi );
            MPE_Log_pack( bytebuf, &bytebuf_pos, 'd', 1, &jj );
        MPE_Log_event( event3b, 0, bytebuf );

        pi = 0.0;
        MPE_Log_event( event4a, 0, NULL );
        MPI_Reduce( &mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            bytebuf_pos = 0;
            MPE_Log_pack( bytebuf, &bytebuf_pos, 'E', 1, &pi );
            MPE_Log_pack( bytebuf, &bytebuf_pos, 'd', 1, &jj );
        MPE_Log_event( event4b, 0, bytebuf );
    }
#if defined( NO_MPI_LOGGING )
    if ( argv != NULL )
        MPE_Finish_log( argv[0] );
    else
        MPE_Finish_log( "cpilog" );
#endif

    if ( myid == 0 ) {
        endwtime = MPI_Wtime();
        printf( "pi is approximately %.16f, Error is %.16f\n",
                pi, fabs(pi - PI25DT) );
        printf( "wall clock time = %f\n", endwtime-startwtime );
    }

    MPI_Finalize();
    return( 0 );
}
