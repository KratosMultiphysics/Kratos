      program test_OpenMP

#ifdef _OPENMP
      use omp_lib
#endif
      implicit none

      include 'mpif.h'
      include 'mpe_logf.h'

      integer world_size, world_rank, ierr
      integer blkA_startevt, blkA_finalevt
      integer blkB_startevt, blkB_finalevt
      integer buffer
      integer max
      parameter (max = 50)
      integer ii

      call MPI_Init( ierr )
      call MPI_Comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      call MPI_Comm_size( MPI_COMM_WORLD, world_size, ierr )

      ierr = MPE_Log_get_state_eventIDs( blkA_startevt, blkA_finalevt )
      ierr = MPE_Log_get_state_eventIDs( blkB_startevt, blkB_finalevt )
      if ( world_rank .eq. 0 ) then
          ierr = MPE_Describe_state( blkA_startevt, blkA_finalevt,
     &                               "Do block A", "orange" )
          ierr = MPE_Describe_state( blkB_startevt, blkB_finalevt,
     &                               "Do block B", "green" )
      endif

      if ( world_rank .eq. 0 ) then
          buffer = 1000
      else
          buffer = 0
      endif
      call MPI_Bcast( buffer, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO
      do ii = 1, max
         ierr = MPE_Log_event( blkA_startevt, 0, '' )
#ifdef _OPENMP
         write(*,*) '[', world_rank, '] ',
     &              'thdID: ', omp_get_thread_num(),
     &              ', ii = ', buffer+ii
#else
         write(*,*) '[', world_rank, '] ', 'ii = ', ii
#endif
         ierr = MPE_Log_event( blkA_finalevt, 0, '' )
      enddo
!$OMP END DO nowait
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii) 
!$OMP DO
      do ii = 1, max
         ierr = MPE_Log_event( blkB_startevt, 0, '' )
#ifdef _OPENMP
         write(*,*) '[', world_rank, '] ',
     &              'thdID: ', omp_get_thread_num(),
     &              ', ii = ', buffer+ii
#else
         write(*,*) '[', world_rank, '] ', 'ii = ', ii
#endif
         ierr = MPE_Log_event( blkB_finalevt, 0, '' )
      enddo
!$OMP END DO nowait
!$OMP END PARALLEL

      call MPI_Finalize( ierr )

      end
