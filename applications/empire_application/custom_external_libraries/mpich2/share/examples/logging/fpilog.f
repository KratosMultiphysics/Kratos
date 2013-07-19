!
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
!
!**********************************************************************
!   pi.f - compute pi by integrating f(x) = 4/(1 + x**2)     
!     
!   Each node: 
!    1) receives the number of rectangles used in the approximation.
!    2) calculates the areas of it's rectangles.
!    3) Synchronizes for a global summation.
!   Node 0 prints the result.
!
!  Variables:
!
!    pi  the calculated result
!    n   number of points of integration.  
!    x           midpoint of each rectangle's interval
!    f           function to integrate
!    sum,pi      area of rectangles
!    tmp         temporary scratch space for global summation
!    i           do loop index
!****************************************************************************
      double precision function f( a )
      implicit none
      double precision a
          f = 4.d0 / (1.d0 + a*a)
          return
      end
!
!
!
      program main
      implicit none

      include 'mpif.h'
      include 'mpe_logf.h'

      double precision  PI25DT
      parameter        (PI25DT = 3.141592653589793238462643d0)

      double precision  mypi, pi, h, sum, x
      integer n, myid, numprocs, ii, idx
      double precision f
      external f

      integer event1a, event1b, event2a, event2b
      integer event3a, event3b, event4a, event4b
      integer ierr

      call MPI_Init( ierr )

      call MPI_Comm_rank( MPI_COMM_WORLD, myid, ierr )
      call MPI_Comm_size( MPI_COMM_WORLD, numprocs, ierr )
      write(6,*) "Process ", myid, " of ", numprocs, " is alive"

! Demonstrate the use of MPE_Log_state_eventIDs() and MPE_Log_solo_eventID()
! which replace the deprecated function MPE_Log_get_event_number.    
!
      ierr = MPE_Log_get_state_eventIDs( event1a, event1b )
      ierr = MPE_Log_get_state_eventIDs( event2a, event2b )
      ierr = MPE_Log_get_solo_eventID( event3a )
      ierr = MPE_Log_get_solo_eventID( event3b )
      ierr = MPE_Log_get_state_eventIDs( event4a, event4b )

! Demonstrate the use MPE_Describe_event() to describe single-timestamped
! drawable, i.e. event.  Caution: One can use MPE_Describe_state() instead
! of 2 MPE_Dresribe_event() calls.  The difference is that one will see
! one state instead of 2 events.
      if ( myid .eq. 0 ) then
          ierr = MPE_Describe_state( event1a, event1b,
     &                               "User_Broadcast", "red" )
          ierr = MPE_Describe_state( event2a, event2b,
     &                               "User_Barrier", "blue" )
          ierr = MPE_Describe_event( event3a, "User_Compute_Start",
     &                               "orange" )
          ierr = MPE_Describe_event( event3b, "User_Compute_Final",
     &                               "orange" )
          ierr = MPE_Describe_state( event4a, event4b,
     &                               "User_Reduce", "green" )
          write(6,*) "event IDs are ", event1a, event1b, ", ",
     &                                 event2a, event2b, ", ",
     &                                 event3a, event3b, ", ",
     &                                 event4a, event4b
      endif

      if ( myid .eq. 0 ) then
!         write(6,98)
! 98      format('Enter the number of intervals: (0 quits)')
!         read(5,99) n
! 99      format(i10)
          n = 1000000
          write(6,*) 'The number of intervals =', n
!         check for quit signal
!         if ( n .le. 0 ) goto 30
      endif

      call MPI_Barrier( MPI_COMM_WORLD, ierr )

      do idx = 1, 5

          ierr = MPE_Log_event( event1a, 0, '' )
          call MPI_Bcast( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
          ierr = MPE_Log_event( event1b, 0, '' )

          call MPI_Pcontrol( 0, ierr )

          ierr = MPE_Log_event( event2a, 0, '' )
          call MPI_Barrier( MPI_COMM_WORLD, ierr )
          ierr = MPE_Log_event( event2b, 0, '' )

          call MPI_Pcontrol( 1, ierr )

          ierr = MPE_Log_event( event3a, 0, '' )
          h = 1.0d0/n
          sum  = 0.0d0
          do ii = myid+1, n, numprocs
              x = h * (dble(ii) - 0.5d0)
              sum = sum + f(x)
          enddo
          mypi = h * sum
          ierr = MPE_Log_event( event3b, 0, '' )

          ierr = MPE_Log_event( event4a, 0, '' )
          pi = 0.0d0
          call MPI_Reduce( mypi, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM,
     &                     0, MPI_COMM_WORLD, ierr )
          ierr = MPE_Log_event( event4b, 0, '' )

          if ( myid .eq. 0 ) then
              write(6, 97) pi, abs(pi - PI25DT)
 97           format('  pi is approximately: ', F18.16,
     +               '  Error is: ', F18.16)
          endif

      enddo
!     - Only GNU fortran does not flush stdout, so calling flush() is
!       absolutely needed with GNU compiler to get all stdout messages.
!     - XLF needs flush_() instead of flush() otherwise needs -qextname=flush
!     - Pathscale fortran compiler needs -intrinsic=G77{or PGI}.
!     call flush(6)

      call MPI_Finalize( ierr )

      end
