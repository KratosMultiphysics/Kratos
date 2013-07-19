/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include "mpi.h"
#include "mpe.h"

#define KBYTES       1024
#define MBYTES       1048576
#define PATH_STRLEN  256

void set_tmpfilename( char*, int );

int main( int argc, char *argv[] )
{
    char    processor_name[ MPI_MAX_PROCESSOR_NAME ];
    char    tmpfilename[ PATH_STRLEN ];
    char    data[ MBYTES ]; 
    double  startwtime, endwtime;
    int     event1a, event1b;
    int     myid, numprocs;
    int     fd, namelen, ii, ierr; 

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &myid );
    MPI_Get_processor_name( processor_name, &namelen );
    fprintf( stderr, "Process %d running on %s\n", myid, processor_name );

    /*  Get event ID from MPE, user should NOT assign event ID  */
    event1a = MPE_Log_get_event_number(); 
    event1b = MPE_Log_get_event_number(); 

    startwtime  = 0.0;
    if ( myid == 0 ) {
        MPE_Describe_state( event1a, event1b, "DiskIO", "purple" );
        startwtime = MPI_Wtime();
    }

    set_tmpfilename( tmpfilename, myid );
    fd = open( tmpfilename, O_RDWR|O_CREAT|O_TRUNC, 0600 );
    if ( fd == -1 ) {
        fprintf( stderr, "open() fails!\n" );
        fflush( stderr );
        MPI_Abort( MPI_COMM_WORLD, -1 );
    }

    MPI_Barrier( MPI_COMM_WORLD );
    for ( ii = 0; ii < 128; ii++ ) {
        MPE_Log_event( event1a, 0, NULL );
        ierr = write( fd, data, 64*KBYTES );
        MPE_Log_event( event1b, 0, NULL );
    }    
    MPI_Barrier( MPI_COMM_WORLD );
    close( fd );

    if ( myid == 0 ) {
        endwtime = MPI_Wtime();
        printf( "wall clock time = %f\n", endwtime-startwtime );
    }
    MPI_Finalize();
    return( 0 );
}

void set_tmpfilename( char *tmp_pathname, int my_rank )
{
    char   *env_tmpdir = NULL;
    char    tmpdirname_ref[ PATH_STRLEN ] = "";
    char    tmpdirname[ PATH_STRLEN ] = "";
    char    tmpfilename[ PATH_STRLEN ] = "";

    /* MPE_TMPDIR takes precedence over TMPDIR */
    env_tmpdir = (char *) getenv( "MPE_TMPDIR" );
    if ( env_tmpdir == NULL )
        env_tmpdir = (char *) getenv( "TMPDIR" );
    if ( env_tmpdir == NULL )
        env_tmpdir = (char *) getenv( "TMP" );
    if ( env_tmpdir == NULL )
        env_tmpdir = (char *) getenv( "TEMP" );

    /*  Set tmpdirname_ref to TMPDIR if available  */
    if ( env_tmpdir != NULL )
        strcat( tmpdirname_ref, env_tmpdir );
    else
#ifdef HAVE_WINDOWS_H
        if ( GetTempPath( PATH_STRLEN, tmpdirname_ref ) == 0 )
            strcat( tmpdirname_ref, "\\");
#else
        strcat( tmpdirname_ref, "/tmp" );
#endif

    if ( env_tmpdir != NULL )
        strcpy( tmpdirname, env_tmpdir );
    else
        strcpy( tmpdirname, tmpdirname_ref );

    if ( strlen( tmpdirname ) <= 0 ) {
        fprintf( stderr, __FILE__":CLOG_Util_get_tmpfilename() - \n"
                         "\t""strlen(tmpdirname) = %d\n",
                         (int)strlen( tmpdirname ) );
        fflush( stderr );
        exit( 1 );
    }

    /*  Set the local tmp filename then tmp_pathname */
    strcpy( tmp_pathname, tmpdirname );
    sprintf( tmpfilename, "/tmp_taskID=%04d_XXXXXX", my_rank );
    strcat( tmp_pathname, tmpfilename );
}
