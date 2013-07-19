/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include "mpi.h"
#include "mpe.h"
#include <math.h>
#include <stdio.h>

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
    int event1, event2, event3;
    char processor_name[ MPI_MAX_PROCESSOR_NAME ];

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

    /*
        user should NOT assign eventIDs directly in MPE_Describe_state()
        Get the eventIDs for user-defined STATES(rectangles) from
        MPE_Log_get_state_eventIDs() instead of the deprecated function
        MPE_Log_get_event_number().
    */
    MPE_Log_get_state_eventIDs( &event1a, &event1b );
    MPE_Log_get_state_eventIDs( &event2a, &event2b );
    MPE_Log_get_state_eventIDs( &event3a, &event3b );
    MPE_Log_get_state_eventIDs( &event4a, &event4b );

    if ( myid == 0 ) {
        MPE_Describe_state( event1a, event1b, "Broadcast", "red" );
        MPE_Describe_state( event2a, event2b, "Sync", "orange" );
        MPE_Describe_state( event3a, event3b, "Compute", "blue" );
        MPE_Describe_state( event4a, event4b, "Reduce", "green" );
    }

    /* Get event ID for Solo-Event(single timestamp object) from MPE */
    MPE_Log_get_solo_eventID( &event1 );
    MPE_Log_get_solo_eventID( &event2 );
    MPE_Log_get_solo_eventID( &event3 );

    if ( myid == 0 ) {
       MPE_Describe_event( event1, "Broadcast Post", "white" );
       MPE_Describe_event( event2, "Compute Start", "purple" );
       MPE_Describe_event( event3, "Compute End", "navy" );
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

    for ( jj = 0; jj < 5; jj++ ) {
        MPE_Log_event( event1a, 0, NULL );
        MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD );
        MPE_Log_event( event1b, 0, NULL );

        MPE_Log_event( event1, 0, NULL );
    
        MPE_Log_event( event2a, 0, NULL );
        MPI_Barrier( MPI_COMM_WORLD );
        MPE_Log_event( event2b, 0, NULL );

        MPE_Log_event( event2, 0, NULL );
        MPE_Log_event( event3a, 0, NULL );
        h   = 1.0 / (double) n;
        sum = 0.0;
        for ( ii = myid + 1; ii <= n; ii += numprocs ) {
            x = h * ((double)ii - 0.5);
            sum += f(x);
        }
        mypi = h * sum;
        MPE_Log_event( event3b, 0, NULL );
        MPE_Log_event( event3, 0, NULL );

        pi = 0.0;
        MPE_Log_event( event4a, 0, NULL );
        MPI_Reduce( &mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MPE_Log_event( event4b, 0, NULL );

        MPE_Log_sync_clocks();
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
