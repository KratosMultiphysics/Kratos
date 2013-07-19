!
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
!

      program main
      implicit none

      include 'mpif.h'

      character*(MPI_MAX_PROCESSOR_NAME)  processor_name
      integer    comm_rank, comm_size, comm_neighbor
      integer    world_rank, world_size, world_neighbor
      integer    icolor, namelen, ibuffer
      integer    splited_comm, duped_comm, inter_comm, comm
      integer    world_request, comm_request
      integer    world_status(MPI_STATUS_SIZE)
      integer    comm_status(MPI_STATUS_SIZE)
      integer    ierr

      call MPI_Init( ierr )
      call MPI_Comm_size( MPI_COMM_WORLD, world_size, ierr )
      call MPI_Comm_rank( MPI_COMM_WORLD, world_rank, ierr )
      call MPI_Get_processor_name( processor_name, namelen, ierr )
      print *, "world_rank ", world_rank, " on ",
     &      processor_name(1:namelen)

      if ( world_rank .eq. world_size - 1 ) then
          world_neighbor = 0
      else
          world_neighbor = world_rank + 1
      endif

      call MPI_Irecv( ibuffer, 1, MPI_INTEGER, MPI_ANY_SOURCE,
     &                99, MPI_COMM_WORLD, world_request, ierr )
      call MPI_Send( world_rank, 1, MPI_INTEGER, world_neighbor,
     &               99, MPI_COMM_WORLD, ierr )
      call MPI_Wait( world_request, world_status, ierr )

!     Split all processes into 2 separate intracommunicators
      icolor  = world_rank - 2 * (world_rank / 2)
      call MPI_Comm_split( MPI_COMM_WORLD, icolor, world_rank,
     &                     splited_comm, ierr )

!     Put in a Comm_dup so local comm ID is different in 2 splited comm
      if ( icolor .eq. 0 ) then
          call MPI_Comm_dup( splited_comm, duped_comm, ierr )
          comm  = duped_comm
      else
          comm  = splited_comm
      endif

      call MPI_Comm_size( comm, comm_size, ierr )
      call MPI_Comm_rank( comm, comm_rank, ierr )

      if ( comm_rank .eq. comm_size - 1 ) then
          comm_neighbor  = 0
      else
          comm_neighbor  = comm_rank + 1
      endif

      call MPI_Irecv( ibuffer, 1, MPI_INTEGER, MPI_ANY_SOURCE,
     &                999, comm, comm_request, ierr )
      call MPI_Send( comm_rank, 1, MPI_INTEGER, comm_neighbor,
     &               999, comm, ierr )
      call MPI_Wait( comm_request, comm_status, ierr )

!     Form an intercomm between the 2 splited intracomm's
      if ( icolor .eq. 0 ) then
          call MPI_Intercomm_create( comm, 0, MPI_COMM_WORLD, 1,
     &                               9090, inter_comm, ierr )
      else
          call MPI_Intercomm_create( comm, 0, MPI_COMM_WORLD, 0,
     &                               9090, inter_comm, ierr )
      endif

      if ( comm_rank .eq. 0 ) then
          call MPI_Irecv( ibuffer, 1, MPI_INTEGER, 0,
     &                    9999, inter_comm, comm_request, ierr )
          call MPI_Send( comm_rank, 1, MPI_INTEGER, 0,
     &                   9999, inter_comm, ierr )
          call MPI_Wait( comm_request, comm_status, ierr )
      endif

!     Free all communicators created
      call MPI_Comm_free( inter_comm, ierr )
      if ( icolor .eq. 0 ) then
          call MPI_Comm_free( duped_comm, ierr )
      endif
      call MPI_Comm_free( splited_comm, ierr )

      call MPI_Finalize( ierr )

      end
