!
! @file spm_user.f90
!
! Fortran 90 example using a laplacian matrix.
!
! @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
program spm_user
  use iso_c_binding
  use spmf
  implicit none

  integer(kind=spm_int_t), dimension(:), pointer :: rowptr
  integer(kind=spm_int_t), dimension(:), pointer :: colptr
  real(kind=c_double),  dimension(:), pointer    :: values
  real(kind=c_double),  dimension(:,:), allocatable, target  :: x0, x, b
  type(c_ptr)                                                :: x0_ptr, x_ptr, b_ptr
  real(kind=c_double)                                        :: eps = 1.e-15
  type(spmatrix_t),        target                            :: spm
  type(spmatrix_t),        target                            :: spm2
  integer(kind=spm_int_t)                                    :: dim1, dim2, dim3, n, nnz
  integer(kind=spm_int_t)                                    :: i, j, k, l, nrhs
  integer(c_int)                                             :: info

  !
  ! Generate a 10x10x10 complex Laplacian in IJV format
  !
  dim1 = 10
  dim2 = 10
  dim3 = 10
  n    = dim1 * dim2 * dim3
  nnz  = (2*(dim1)-1) * dim2 * dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1)

  !
  ! Create the spm out of the internal data
  !
  call spmInit( spm )
  spm%mtxtype = SpmSymmetric
  spm%flttype = SpmDouble
  spm%fmttype = SpmIJV
  spm%n       = n
  spm%nnz     = nnz
  spm%dof     = 1

  call spmUpdateComputedFields( spm )
  call spmAlloc( spm )

  call c_f_pointer( spm%rowptr, rowptr, [nnz] )
  call c_f_pointer( spm%colptr, colptr, [nnz] )
  call c_f_pointer( spm%values, values, [nnz] )

  l = 1
  do i=1,dim1
     do j=1,dim2
        do k=1,dim3
           rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
           values(l) = 6.

           if (i == 1) then
              values(l) = values(l) - 1.
           end if
           if (i == dim1) then
              values(l) = values(l) - 1.
           end if
           if (j == 1) then
              values(l) = values(l) - 1.
           end if
           if (j == dim2) then
              values(l) = values(l) - 1.
           end if
           if (k == 1) then
              values(l) = values(l) - 1.
           end if
           if (k == dim3) then
              values(l) = values(l) - 1.
           end if

           values(l) = values(l) * 8.
           l = l + 1

           if (i < dim1) then
              rowptr(l) =  i    + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = - 1.
              l = l + 1
           end if
           if (j < dim2) then
              rowptr(l) = (i-1) + dim1 *  j    + dim1 * dim2 * (k-1) + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = - 1.
              l = l + 1
           end if
           if (k < dim3) then
              rowptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 *  k    + 1
              colptr(l) = (i-1) + dim1 * (j-1) + dim1 * dim2 * (k-1) + 1
              values(l) = -1.
              l = l + 1
           end if
        end do
     end do
  end do

  if ( l .ne. nnz+1 ) then
     write(6,*) 'l ', l, " nnz ", nnz
  end if

  call spmCheckAndCorrect( spm, spm2, info )
  if ( info .ne. 0 ) then
     call spmExit( spm )
     spm = spm2
  end if

  call spmPrintInfo( spm )

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

end program spm_user
