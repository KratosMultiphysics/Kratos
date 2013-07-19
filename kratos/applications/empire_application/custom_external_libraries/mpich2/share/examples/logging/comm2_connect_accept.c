/*
   (C) 2001 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/
#define _GNU_SOURCE

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"




static int is_server = 0;
static int is_client = 0;

void handle_error( int errcode, char *str );
void handle_error( int errcode, char *str )
{
    char msg[ MPI_MAX_ERROR_STRING ];
    int  resultlen;
    MPI_Error_string( errcode, msg, &resultlen );
    fprintf( stderr, "%s: %s\n", str, msg );
}

/*
   open a port, waiting for a connection from a client,
*/
MPI_Comm server_init( MPI_Comm comm );
MPI_Comm server_init( MPI_Comm comm )
{
    char      port_name[ MPI_MAX_PORT_NAME ];
    MPI_Comm  newcomm;

    MPI_Open_port( MPI_INFO_NULL, port_name );
    fprintf( stdout, "server: port opened at %s\n", port_name );
    MPI_Publish_name( "mpe_port_name", MPI_INFO_NULL, port_name );
    MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, comm, &newcomm );

    return newcomm;
}

/*
   look up the available port, then connect to the server with the port name.
*/
MPI_Comm client_init( MPI_Comm comm );
MPI_Comm client_init( MPI_Comm comm )
{
    MPI_Comm  newcomm;
    int       ret;
    char      port_name[ MPI_MAX_PORT_NAME ];

    ret = MPI_Lookup_name( "mpe_port_name", MPI_INFO_NULL, port_name );
    if ( ret != MPI_SUCCESS ) {
        handle_error( ret, "MPI_Lookup_name" );
        return 0;
    }
    fprintf( stdout, "client: found open port at %s\n", port_name );
    MPI_Comm_connect( port_name, MPI_INFO_NULL, 0, comm, &newcomm );

    return newcomm;
}

void usage( char * name );
void usage( char * name )
{
    fprintf( stderr, "usage: %s [-s|-c]\n", name );
    fprintf( stderr, "      specify one and only one of -s or -c\n" );
    exit( -1 );
}

int parse_args( int argc, char ** argv );
int parse_args( int argc, char ** argv )
{
    int c;
    while ( (c = getopt( argc, argv, "csp:" ) ) != -1 ) {
        switch (c) {
            case 's':
                is_server = 1;
                break;
            case 'c':
                is_client = 1;
                break;
            case '?':
            case ':':
            default:
                usage(argv[0]);
        }
    }
    if ( (is_client == 0 ) && (is_server == 0) ) {
        usage( argv[0] );
    }
    return 0;
}


int main ( int argc, char ** argv )
{
    MPI_Comm     intercomm, intracomm;
    MPI_Request  request;
    MPI_Status   status;
    int          ibuffer;
    int          neighbor_rank, rank, size;
    int         *usize, aflag;
    


    MPI_Init(&argc, &argv);
#if 1
    /* temporary hack for MPICH2: if we inquire about MPI_UNIVERSE_SIZE,
     * MPICH2 will promote our singleton init into a full-fleged MPI
     * environment */
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE, &usize, &aflag);
#endif

    intercomm = MPI_COMM_NULL;

#if defined( SERVER )
    is_server = 1;
    intercomm = server_init( MPI_COMM_WORLD );
#elif defined( CLIENT )
    is_client = 1;
    intercomm = client_init( MPI_COMM_WORLD );
#else
    parse_args( argc, argv );
    if ( is_server ) {
        intercomm = server_init( MPI_COMM_WORLD );
    }
    else if ( is_client ) {
        intercomm = client_init( MPI_COMM_WORLD );
    }
#endif

    if ( intercomm == MPI_COMM_NULL ) {
        if ( is_server ) {
            fprintf( stderr, "Server returns NULL intercommunicator!" );
        }
        else if ( is_client ) {
            fprintf( stderr, "Client returns NULL intercommunicator!" );
        }
        else {
            fprintf( stderr, "Unknown server/client: NULL intercommunicator!" );
        }
        return -1;
    }

    MPI_Comm_rank( intercomm, &rank );

    if ( rank == 0 ) {
        MPI_Irecv( &ibuffer, 1, MPI_INT, 0,
                   9999, intercomm, &request );
        MPI_Send( &rank, 1, MPI_INT, 0, 9999, intercomm );
        MPI_Wait( &request, &status );
    }


    MPI_Intercomm_merge( intercomm, 0, &intracomm );
    MPI_Comm_rank( intracomm, &rank );
    MPI_Comm_size( intracomm, &size );

    fprintf( stdout, "[%d/%d] after Intercomm_merge()\n", rank, size );

    if ( rank == size - 1 )
        neighbor_rank = 0;
    else
        neighbor_rank = rank + 1;

    MPI_Irecv( &ibuffer, 1, MPI_INT, MPI_ANY_SOURCE,
               999, intracomm, &request );
    MPI_Send( &rank, 1, MPI_INT, neighbor_rank, 999, intracomm );
    MPI_Wait( &request, &status );

    MPI_Comm_free( &intracomm );
    MPI_Comm_disconnect( &intercomm );
    MPI_Finalize();
    return 0;
}
