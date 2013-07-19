!
!  (C) 2004 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
!

      program  wrong_reals
      implicit none

      include 'mpif.h'

      integer myrank, numprocs, ierr
      real*8  dval
      real*4  fval

      call MPI_INIT( ierr )

      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      print *, "Process ", myrank, " of ", numprocs, " is alive"

      if ( myrank .eq. numprocs-1 ) then
          call MPI_BCAST( fval, 2, MPI_REAL4,
     &                    0, MPI_COMM_WORLD, ierr ) 
      else
          call MPI_BCAST( dval, 1, MPI_REAL8,
     &                    0, MPI_COMM_WORLD, ierr ) 
      endif

      call MPI_Finalize( ierr )

      end
