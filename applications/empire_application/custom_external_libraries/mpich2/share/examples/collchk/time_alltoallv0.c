/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

/*
  This program tests MPI_Alltoallv by having processor each process 
  send data to two neighbors only, using counts of 0 for the other processes.
  This idiom is sometimes used for halo exchange operations.

  Because there are separate send and receive types to alltoallv,
  there need to be tests to rearrange data on the fly.  Not done yet.
  
  Currently, the test uses only MPI_INT; this is adequate for testing systems
  that use point-to-point operations
 */

int main( int argc, char **argv )
{

    MPI_Comm comm;
    double   time_init, time_final;
    double   *sbuf, *rbuf, *p;
    int      rank, size;
    int      *sendcounts, *recvcounts, *rdispls, *sdispls;
    int      num_itr, idx, ii, err;
    int      left, right, length;
    
    MPI_Init( &argc, &argv );
    err = 0;
    
      comm = MPI_COMM_WORLD;

      MPI_Comm_size( comm, &size );
      MPI_Comm_rank( comm, &rank );

      if ( argc > 1 && argv[1] != NULL )
          length = atoi( argv[1] );
      else
          length = 1;

      if ( argc > 2 && argv[2] != NULL )
          num_itr = atoi( argv[2] );
      else
          num_itr = 1;
      
      /* Create and load the arguments to alltoallv */
      sendcounts = (int *)malloc( size * sizeof(int) );
      recvcounts = (int *)malloc( size * sizeof(int) );
      rdispls    = (int *)malloc( size * sizeof(int) );
      sdispls    = (int *)malloc( size * sizeof(int) );
      if (!sendcounts || !recvcounts || !rdispls || !sdispls) {
        fprintf( stderr, "Could not allocate arg items!\n" );
        MPI_Abort( comm, 1 );
      }

      /* Get the neighbors */
      left  = (rank - 1 + size) % size;
      right = (rank + 1) % size;

      /* Set the defaults */
      for ( ii = 0; ii < size; ii++) {
          sendcounts[ii] = 0;
          recvcounts[ii] = 0;
          rdispls[ii]    = 0;
          sdispls[ii]    = 0;
      }

      /* Get the buffers */
      sbuf = (double *)malloc( 2 * length * sizeof(double) );
      rbuf = (double *)malloc( 2 * length * sizeof(double) );
      if (!sbuf || !rbuf) {
          fprintf( stderr, "Could not allocate buffers!\n" );
          MPI_Abort( comm, 1 );
      }

      /* Load up the buffers */
      for ( ii = 0; ii < length; ii++) {
          sbuf[ii]        = ii + 100000*rank;
          sbuf[ii+length] = ii + 100000*rank;
          rbuf[ii]        = -ii;
          rbuf[ii+length] = -ii-length;
      }
      sendcounts[left]  = length;
      sendcounts[right] = length;
      recvcounts[left]  = length;
      recvcounts[right] = length;
      rdispls[left]     = 0;
      rdispls[right]    = length;
      sdispls[left]     = 0;
      sdispls[right]    = length;

      MPI_Barrier( comm );
      MPI_Barrier( comm );
      time_init   = MPI_Wtime();

      for ( idx = 0; idx < num_itr; idx++ ) {
          MPI_Alltoallv( sbuf, sendcounts, sdispls, MPI_DOUBLE,
                         rbuf, recvcounts, rdispls, MPI_DOUBLE, comm );
      }

      /* MPI_Barrier( comm ); */
      time_final  = MPI_Wtime();

      fprintf( stdout, "time taken by %dx%d MPI_AlltoAllv() at rank %d = %f\n",
                       length, num_itr, rank, time_final - time_init );

      /* Check rbuf */
      p = rbuf;          /* left */

      for ( ii = 0; ii < length; ii++) {
          if (p[ii] != ii + 100000 * left) {
              if (err < 10) {
                  fprintf( stderr, "[%d from %d] got %f expected %d for %dth\n",
                           rank, left, p[ii], ii + 100000 * left, ii );
              }
              err++;
          }
      }

      p = rbuf + length; /* right */
      for ( ii = 0; ii < length; ii++) {
          if (p[ii] != ii + 100000 * right) {
              if (err < 10) {
                  fprintf( stderr, "[%d from %d] got %f expected %d for %dth\n",
                           rank, right, p[ii], ii + 100000 * right, ii );
              }
              err++;
          }
      }

      free( rbuf );
      free( sbuf );
      free( sdispls );
      free( rdispls );
      free( recvcounts );
      free( sendcounts );

    MPI_Finalize();
    return 0;
}
