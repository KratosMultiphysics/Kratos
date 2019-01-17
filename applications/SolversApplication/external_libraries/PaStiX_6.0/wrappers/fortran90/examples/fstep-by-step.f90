!
! @file fstep-by-step.f90
!
! Fortran 90 example using a matrix read with the spm driver.
!
! @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.1
! @author Mathieu Faverge
! @date 2018-07-16
!
program fsimple
  use iso_c_binding
  use pastix_enums
  use spmf
  use pastixf
  ! use mpi_f08
  implicit none

  type(pastix_data_t),        pointer                      :: pastix_data
  type(pastix_order_t),       pointer                      :: order => null()
  type(spmatrix_t),           pointer                      :: spm
  type(spmatrix_t),           pointer                      :: spm2
  integer(kind=pastix_int_t), target                       :: iparm(iparm_size)
  real(kind=c_double),        target                       :: dparm(dparm_size)
  integer(c_int)                                           :: info
  integer(kind=pastix_int_t)                               :: nrhs
  real(kind=c_double), dimension(:,:), allocatable, target :: x0, x, b
  type(c_ptr)                                              :: x0_ptr, x_ptr, b_ptr
  integer(kind=pastix_int_t), dimension(:), pointer        :: permtab
  integer                                                  :: i, j, nfact, nsolv

  nfact = 2
  nsolv = 3

  !
  ! Initialize the problem
  !   1- The matrix
  allocate( spm )
  call spmReadDriver( SpmDriverLaplacian, "d:10:10:10:2.", spm, info )

  allocate( spm2 )
  call spmCheckAndCorrect( spm, spm2, info )
  if ( info .ne. 0 ) then
     call spmExit( spm )
     spm = spm2
  end if
  deallocate( spm2 )

  call spmPrintInfo( spm )

  !   2- The right hand side
  nrhs = 10
  allocate(x0(spm%n,nrhs))
  allocate(x( spm%n,nrhs))
  allocate(b( spm%n,nrhs))
  x0_ptr = c_loc(x0)
  x_ptr  = c_loc(x)
  b_ptr  = c_loc(b)

  ! 1- Initialize the parameters and the solver
  call pastixInitParam( iparm, dparm )
  call pastixInit( pastix_data, 0, iparm, dparm )

  ! 2- Perform ordering, symbolic factorization, and analyze steps
  call pastix_subtask_order( pastix_data, spm, order, info )
  call pastix_subtask_symbfact( pastix_data, info )
  call pastix_subtask_reordering( pastix_data, info )
  call pastix_subtask_blend( pastix_data, info )

  ! If needed, get the generated ordering
  call pastixOrderGet( pastix_data, order )

  ! Convert the permtab to Fortran array
  call c_f_pointer( order%permtab, permtab, [order%vertnbr] )
  print *, permtab(1:10)

  ! 3- Factorize nfact times the matrix
  do i=0,nfact
     ! Perform the numerical factorization
     call pastix_subtask_spm2bcsc( pastix_data, spm, info )
     call pastix_subtask_bcsc2ctab( pastix_data, info )
     call pastix_subtask_sopalin( pastix_data, info )

     ! Perform nsolv solve steps
     do j=0,nsolv

        call spmGenRHS( SpmRhsRndX, nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, info )
        x = b

        ! 4- Solve the problem
        call pastix_task_solve( pastix_data, nrhs, x_ptr, spm%n, info )

        ! 5- Refine the solution
        call pastix_task_refine( pastix_data, spm%n, nrhs, b_ptr, spm%n, x_ptr, spm%n, info )

        ! Check the solution
        call spmCheckAxb( dparm(DPARM_EPSILON_REFINEMENT), nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, x_ptr, spm%n, info )

     end do
  end do

  ! 6- Destroy the C data structure
  call pastixFinalize( pastix_data )

  call spmExit( spm )
  deallocate( spm )
  deallocate(x0)
  deallocate(x)
  deallocate(b)

end program fsimple
