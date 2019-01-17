!
! @file spm_driver.f90
!
! Fortran 90 example using a matrix read with the spm driver.
!
! @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
program spm_driver
  use iso_c_binding
  use spmf
  ! use mpi_f08
  implicit none

  type(spmatrix_t),           target                       :: spm
  type(spmatrix_t),           target                       :: spm2
  real(kind=c_double)                                      :: normA
  real(kind=c_double)                                      :: eps = 1.e-15
  integer(c_int)                                           :: info
  integer(kind=spm_int_t)                                  :: nrhs
  real(kind=c_double), dimension(:,:), allocatable, target :: x0, x, b
  type(c_ptr)                                              :: x0_ptr, x_ptr, b_ptr

  !
  ! Initialize the problem
  !   1- The matrix
  call spmReadDriver( SpmDriverLaplacian, "d:10:10:10:4.", spm, info )

  call spmCheckAndCorrect( spm, spm2, info )
  if ( info .ne. 0 ) then
     call spmExit( spm )
     spm = spm2
  end if

  call spmPrintInfo( spm )

  ! Scale A for better stability with low-rank computations
  call spmNorm( SpmFrobeniusNorm, spm, normA )
  call spmScalMatrix( 1. / normA, spm )

  !   2- The right hand side
  nrhs = 10
  allocate(x0(spm%n, nrhs))
  allocate(x( spm%n, nrhs))
  allocate(b( spm%n, nrhs))
  x0_ptr = c_loc(x0)
  x_ptr  = c_loc(x)
  b_ptr  = c_loc(b)

  ! Compute b = A * x, with x random
  call spmGenRHS( SpmRhsRndX, nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, info )

  ! Copy x0 into x
  x = x0

  !
  ! Check the solution
  !
  call spmCheckAxb( eps, nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, x_ptr, spm%n, info )

  call spmExit( spm )
  deallocate(x0)
  deallocate(x)
  deallocate(b)

end program spm_driver
