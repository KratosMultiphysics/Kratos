/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mpe_log.h"

#if defined( USEC_TIMING )
#include "rdtsc.h"
#define MYTIME_T  unsigned long long
#endif

#define ITER_COUNT  5

int main( int argc, char *argv[] )
{
    int    event1a, event1b, event2a, event2b;
    int    jj;

    MPE_LOG_BYTES  bytebuf;
    int            bytebuf_pos;

#if defined( USEC_TIMING )
    MYTIME_T       Atimes[ITER_COUNT], Btimes[ITER_COUNT];
#endif

    /*
        MPE_Init_log() & MPE_Finish_log() are NOT needed when
        liblmpe.a is linked with this program.  In that case,
        MPI_Init() would have called MPE_Init_log() already.
    */
    MPE_Init_log();

    /*  Get event ID from MPE, user should NOT assign event ID directly */
    event1a = MPE_Log_get_event_number(); 
    event1b = MPE_Log_get_event_number(); 
    event2a = MPE_Log_get_event_number(); 
    event2b = MPE_Log_get_event_number(); 

    MPE_Describe_state( event1a, event1b, "Blank", "red" );
    MPE_Describe_info_state( event2a, event2b, "GetRank", "orange",
                             "source = %s()'s line %d." );

    /*
    MPE_Start_log();
    */

    for ( jj = 0; jj < ITER_COUNT; jj++ ) {
#if defined( USEC_TIMING )
        TIME_PRE( Atimes[jj] );
#endif
        MPE_Log_event( event1a, 0, NULL );
        MPE_Log_event( event1b, 0, NULL );
#if defined( USEC_TIMING )
        TIME_POST( Atimes[jj] );
#endif
    }
    
    for ( jj = 0; jj < ITER_COUNT; jj++ ) {
#if defined( USEC_TIMING )
        TIME_PRE( Btimes[jj] );
#endif
        MPE_Log_event( event2a, 0, NULL );
            int line_num;
            bytebuf_pos = 0;
            MPE_Log_pack( bytebuf, &bytebuf_pos, 's',
                          sizeof(__func__)-1, __func__ );
            line_num = __LINE__;
            MPE_Log_pack( bytebuf, &bytebuf_pos, 'd', 1, &line_num );
        MPE_Log_event( event2b, 0, bytebuf );
#if defined( USEC_TIMING )
        TIME_POST( Btimes[jj] );
#endif
    }

    MPE_Finish_log( argv[0] );

#if defined( USEC_TIMING )
    for ( jj = 0; jj < ITER_COUNT; jj++ ) {
         printf( "Atimes[%d] = %5.3f\n", jj, USECS( Atimes[jj] ) );
         printf( "Btimes[%d] = %5.3f\n", jj, USECS( Btimes[jj] ) );
    }
#endif

    return( 0 );
}
