
!
! @file spmf.f90
!
! SPM Fortran 90 wrapper
!
! @copyright 2017      Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
! This file has been automatically generated with gen_wrappers.py
!
module spmf
  use iso_c_binding
  use spm_enums
  implicit none

  type, bind(c) :: spmatrix_t
     integer(c_int)          :: mtxtype
     integer(c_int)          :: flttype
     integer(c_int)          :: fmttype
     integer(kind=spm_int_t) :: gN
     integer(kind=spm_int_t) :: n
     integer(kind=spm_int_t) :: gnnz
     integer(kind=spm_int_t) :: nnz
     integer(kind=spm_int_t) :: gNexp
     integer(kind=spm_int_t) :: nexp
     integer(kind=spm_int_t) :: gnnzexp
     integer(kind=spm_int_t) :: nnzexp
     integer(kind=spm_int_t) :: dof
     type(c_ptr)             :: dofs
     integer(c_int)          :: layout
     type(c_ptr)             :: colptr
     type(c_ptr)             :: rowptr
     type(c_ptr)             :: loc2glob
     type(c_ptr)             :: values
  end type spmatrix_t

  interface
     subroutine spmInit_c(spm) &
          bind(c, name='spmInit')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmInit_c
  end interface

  interface
     subroutine spmAlloc_c(spm) &
          bind(c, name='spmAlloc')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmAlloc_c
  end interface

  interface
     subroutine spmExit_c(spm) &
          bind(c, name='spmExit')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmExit_c
  end interface

  interface
     function spmCopy_c(spm) &
          bind(c, name='spmCopy')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr)        :: spmCopy_c
       type(c_ptr), value :: spm
     end function spmCopy_c
  end interface

  interface
     subroutine spmBase_c(spm, baseval) &
          bind(c, name='spmBase')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: baseval
     end subroutine spmBase_c
  end interface

  interface
     function spmFindBase_c(spm) &
          bind(c, name='spmFindBase')
       use iso_c_binding
       import spm_int_t
       import spmatrix_t
       implicit none
       integer(kind=spm_int_t)   :: spmFindBase_c
       type(c_ptr),        value :: spm
     end function spmFindBase_c
  end interface

  interface
     function spmConvert_c(ofmttype, ospm) &
          bind(c, name='spmConvert')
       use iso_c_binding
       import spmatrix_t
       implicit none
       integer(kind=c_int)        :: spmConvert_c
       integer(kind=c_int), value :: ofmttype
       type(c_ptr),         value :: ospm
     end function spmConvert_c
  end interface

  interface
     subroutine spmUpdateComputedFields_c(spm) &
          bind(c, name='spmUpdateComputedFields')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmUpdateComputedFields_c
  end interface

  interface
     subroutine spmGenFakeValues_c(spm) &
          bind(c, name='spmGenFakeValues')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spmGenFakeValues_c
  end interface

  interface
     function spmNorm_c(ntype, spm) &
          bind(c, name='spmNorm')
       use iso_c_binding
       import spmatrix_t
       implicit none
       real(kind=c_double)   :: spmNorm_c
       integer(c_int), value :: ntype
       type(c_ptr),    value :: spm
     end function spmNorm_c
  end interface

  interface
     function spmMatVec_c(trans, alpha, spm, x, beta, y) &
          bind(c, name='spmMatVec')
       use iso_c_binding
       import spmatrix_t
       implicit none
       integer(kind=c_int)        :: spmMatVec_c
       integer(c_int),      value :: trans
       real(kind=c_double), value :: alpha
       type(c_ptr),         value :: spm
       type(c_ptr),         value :: x
       real(kind=c_double), value :: beta
       type(c_ptr),         value :: y
     end function spmMatVec_c
  end interface

  interface
     function spmMatMat_c(trans, n, alpha, A, B, ldb, beta, C, ldc) &
          bind(c, name='spmMatMat')
       use iso_c_binding
       import spm_int_t
       import spmatrix_t
       implicit none
       integer(kind=c_int)            :: spmMatMat_c
       integer(c_int),          value :: trans
       integer(kind=spm_int_t), value :: n
       real(kind=c_double),     value :: alpha
       type(c_ptr),             value :: A
       type(c_ptr),             value :: B
       integer(kind=spm_int_t), value :: ldb
       real(kind=c_double),     value :: beta
       type(c_ptr),             value :: C
       integer(kind=spm_int_t), value :: ldc
     end function spmMatMat_c
  end interface

  interface
     subroutine spmScalMatrix_c(alpha, spm) &
          bind(c, name='spmScalMatrix')
       use iso_c_binding
       import spmatrix_t
       implicit none
       real(kind=c_double), value :: alpha
       type(c_ptr),         value :: spm
     end subroutine spmScalMatrix_c
  end interface

  interface
     subroutine spmScalVector_c(flt, alpha, n, x, incx) &
          bind(c, name='spmScalVector')
       use iso_c_binding
       import spm_int_t
       implicit none
       integer(c_int),          value :: flt
       real(kind=c_double),     value :: alpha
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: incx
     end subroutine spmScalVector_c
  end interface

  interface
     function spmSort_c(spm) &
          bind(c, name='spmSort')
       use iso_c_binding
       import spmatrix_t
       implicit none
       integer(kind=c_int)   :: spmSort_c
       type(c_ptr),    value :: spm
     end function spmSort_c
  end interface

  interface
     function spmMergeDuplicate_c(spm) &
          bind(c, name='spmMergeDuplicate')
       use iso_c_binding
       import spm_int_t
       import spmatrix_t
       implicit none
       integer(kind=spm_int_t)   :: spmMergeDuplicate_c
       type(c_ptr),        value :: spm
     end function spmMergeDuplicate_c
  end interface

  interface
     function spmSymmetrize_c(spm) &
          bind(c, name='spmSymmetrize')
       use iso_c_binding
       import spm_int_t
       import spmatrix_t
       implicit none
       integer(kind=spm_int_t)   :: spmSymmetrize_c
       type(c_ptr),        value :: spm
     end function spmSymmetrize_c
  end interface

  interface
     function spmCheckAndCorrect_c(spm_in, spm_out) &
          bind(c, name='spmCheckAndCorrect')
       use iso_c_binding
       import spmatrix_t
       implicit none
       integer(kind=c_int)   :: spmCheckAndCorrect_c
       type(c_ptr),    value :: spm_in
       type(c_ptr),    value :: spm_out
     end function spmCheckAndCorrect_c
  end interface

  interface
     function spmGenRHS_c(type, nrhs, spm, x, ldx, b, ldb) &
          bind(c, name='spmGenRHS')
       use iso_c_binding
       import spm_int_t
       import spmatrix_t
       implicit none
       integer(kind=c_int)            :: spmGenRHS_c
       integer(c_int),          value :: type
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: ldx
       type(c_ptr),             value :: b
       integer(kind=spm_int_t), value :: ldb
     end function spmGenRHS_c
  end interface

  interface
     function spmCheckAxb_c(eps, nrhs, spm, x0, ldx0, b, ldb, x, ldx) &
          bind(c, name='spmCheckAxb')
       use iso_c_binding
       import spm_int_t
       import spmatrix_t
       implicit none
       integer(kind=c_int)            :: spmCheckAxb_c
       real(kind=c_double),     value :: eps
       integer(kind=spm_int_t), value :: nrhs
       type(c_ptr),             value :: spm
       type(c_ptr),             value :: x0
       integer(kind=spm_int_t), value :: ldx0
       type(c_ptr),             value :: b
       integer(kind=spm_int_t), value :: ldb
       type(c_ptr),             value :: x
       integer(kind=spm_int_t), value :: ldx
     end function spmCheckAxb_c
  end interface

  interface
     function spmIntConvert_c(n, input) &
          bind(c, name='spmIntConvert')
       use iso_c_binding
       import spm_int_t
       implicit none
       type(c_ptr)                    :: spmIntConvert_c
       integer(kind=spm_int_t), value :: n
       type(c_ptr),             value :: input
     end function spmIntConvert_c
  end interface

  interface
     function spmLoad_c(spm, infile) &
          bind(c, name='spmLoad')
       use iso_c_binding
       import spmatrix_t
       implicit none
       integer(kind=c_int)   :: spmLoad_c
       type(c_ptr),    value :: spm
       type(c_ptr),    value :: infile
     end function spmLoad_c
  end interface

  interface
     function spmSave_c(spm, outfile) &
          bind(c, name='spmSave')
       use iso_c_binding
       import spmatrix_t
       implicit none
       integer(kind=c_int)   :: spmSave_c
       type(c_ptr),    value :: spm
       type(c_ptr),    value :: outfile
     end function spmSave_c
  end interface

  interface
     function spmReadDriver_c(driver, filename, spm) &
          bind(c, name='spmReadDriver')
       use iso_c_binding
       import spmatrix_t
       implicit none
       integer(kind=c_int)   :: spmReadDriver_c
       integer(c_int), value :: driver
       type(c_ptr),    value :: filename
       type(c_ptr),    value :: spm
     end function spmReadDriver_c
  end interface

  interface
     function spmParseLaplacianInfo_c(filename, flttype, dim1, dim2, dim3, &
          alpha, beta) &
          bind(c, name='spmParseLaplacianInfo')
       use iso_c_binding
       import spm_int_t
       implicit none
       integer(kind=c_int)   :: spmParseLaplacianInfo_c
       type(c_ptr),    value :: filename
       type(c_ptr),    value :: flttype
       type(c_ptr),    value :: dim1
       type(c_ptr),    value :: dim2
       type(c_ptr),    value :: dim3
       type(c_ptr),    value :: alpha
       type(c_ptr),    value :: beta
     end function spmParseLaplacianInfo_c
  end interface

  interface
     subroutine spm2Dense_c(spm) &
          bind(c, name='spm2Dense')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
     end subroutine spm2Dense_c
  end interface

  interface
     subroutine spmPrint_c(spm, f) &
          bind(c, name='spmPrint')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
       type(c_ptr), value :: f
     end subroutine spmPrint_c
  end interface

  interface
     subroutine spmPrintInfo_c(spm, f) &
          bind(c, name='spmPrintInfo')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr), value :: spm
       type(c_ptr), value :: f
     end subroutine spmPrintInfo_c
  end interface

  interface
     function spmExpand_c(spm) &
          bind(c, name='spmExpand')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr)        :: spmExpand_c
       type(c_ptr), value :: spm
     end function spmExpand_c
  end interface

  interface
     function spmDofExtend_c(spm, type, dof) &
          bind(c, name='spmDofExtend')
       use iso_c_binding
       import spmatrix_t
       implicit none
       type(c_ptr)                :: spmDofExtend_c
       type(c_ptr),         value :: spm
       integer(kind=c_int), value :: type
       integer(kind=c_int), value :: dof
     end function spmDofExtend_c
  end interface

contains

  ! Wrappers of the C functions.
  subroutine spmInit(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(inout), target :: spm

    call spmInit_c(c_loc(spm))
  end subroutine spmInit

  subroutine spmAlloc(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(inout), target :: spm

    call spmAlloc_c(c_loc(spm))
  end subroutine spmAlloc

  subroutine spmExit(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(inout), target :: spm

    call spmExit_c(c_loc(spm))
  end subroutine spmExit

  subroutine spmCopy(spm, spmo)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(in),  target  :: spm
    type(spmatrix_t), intent(out), pointer :: spmo

    call c_f_pointer(spmCopy_c(c_loc(spm)), spmo)
  end subroutine spmCopy

  subroutine spmBase(spm, baseval)
    use iso_c_binding
    implicit none
    type(spmatrix_t),    intent(inout), target :: spm
    integer(kind=c_int), intent(in)            :: baseval

    call spmBase_c(c_loc(spm), baseval)
  end subroutine spmBase

  subroutine spmFindBase(spm, value)
    use iso_c_binding
    implicit none
    type(spmatrix_t),        intent(in), target :: spm
    integer(kind=spm_int_t), intent(out)        :: value

    value = spmFindBase_c(c_loc(spm))
  end subroutine spmFindBase

  subroutine spmConvert(ofmttype, ospm, info)
    use iso_c_binding
    implicit none
    integer(kind=c_int), intent(in)            :: ofmttype
    type(spmatrix_t),    intent(inout), target :: ospm
    integer(kind=c_int), intent(out)           :: info

    info = spmConvert_c(ofmttype, c_loc(ospm))
  end subroutine spmConvert

  subroutine spmUpdateComputedFields(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(inout), target :: spm

    call spmUpdateComputedFields_c(c_loc(spm))
  end subroutine spmUpdateComputedFields

  subroutine spmGenFakeValues(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(inout), target :: spm

    call spmGenFakeValues_c(c_loc(spm))
  end subroutine spmGenFakeValues

  subroutine spmNorm(ntype, spm, value)
    use iso_c_binding
    implicit none
    integer(c_int),      intent(in)         :: ntype
    type(spmatrix_t),    intent(in), target :: spm
    real(kind=c_double), intent(out)        :: value

    value = spmNorm_c(ntype, c_loc(spm))
  end subroutine spmNorm

  subroutine spmMatVec(trans, alpha, spm, x, beta, y, info)
    use iso_c_binding
    implicit none
    integer(c_int),      intent(in)            :: trans
    real(kind=c_double), intent(in)            :: alpha
    type(spmatrix_t),    intent(in),    target :: spm
    type(c_ptr),         intent(in),    target :: x
    real(kind=c_double), intent(in)            :: beta
    type(c_ptr),         intent(inout), target :: y
    integer(kind=c_int), intent(out)           :: info

    info = spmMatVec_c(trans, alpha, c_loc(spm), x, beta, y)
  end subroutine spmMatVec

  subroutine spmMatMat(trans, n, alpha, A, B, ldb, beta, C, ldc, info)
    use iso_c_binding
    implicit none
    integer(c_int),          intent(in)            :: trans
    integer(kind=spm_int_t), intent(in)            :: n
    real(kind=c_double),     intent(in)            :: alpha
    type(spmatrix_t),        intent(in),    target :: A
    type(c_ptr),             intent(in),    target :: B
    integer(kind=spm_int_t), intent(in)            :: ldb
    real(kind=c_double),     intent(in)            :: beta
    type(c_ptr),             intent(inout), target :: C
    integer(kind=spm_int_t), intent(in)            :: ldc
    integer(kind=c_int),     intent(out)           :: info

    info = spmMatMat_c(trans, n, alpha, c_loc(A), B, ldb, beta, C, ldc)
  end subroutine spmMatMat

  subroutine spmScalMatrix(alpha, spm)
    use iso_c_binding
    implicit none
    real(kind=c_double), intent(in)            :: alpha
    type(spmatrix_t),    intent(inout), target :: spm

    call spmScalMatrix_c(alpha, c_loc(spm))
  end subroutine spmScalMatrix

  subroutine spmScalVector(flt, alpha, n, x, incx)
    use iso_c_binding
    implicit none
    integer(c_int),          intent(in)            :: flt
    real(kind=c_double),     intent(in)            :: alpha
    integer(kind=spm_int_t), intent(in)            :: n
    type(c_ptr),             intent(inout), target :: x
    integer(kind=spm_int_t), intent(in)            :: incx

    call spmScalVector_c(flt, alpha, n, x, incx)
  end subroutine spmScalVector

  subroutine spmSort(spm, info)
    use iso_c_binding
    implicit none
    type(spmatrix_t),    intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = spmSort_c(c_loc(spm))
  end subroutine spmSort

  subroutine spmMergeDuplicate(spm, value)
    use iso_c_binding
    implicit none
    type(spmatrix_t),        intent(inout), target :: spm
    integer(kind=spm_int_t), intent(out)           :: value

    value = spmMergeDuplicate_c(c_loc(spm))
  end subroutine spmMergeDuplicate

  subroutine spmSymmetrize(spm, value)
    use iso_c_binding
    implicit none
    type(spmatrix_t),        intent(inout), target :: spm
    integer(kind=spm_int_t), intent(out)           :: value

    value = spmSymmetrize_c(c_loc(spm))
  end subroutine spmSymmetrize

  subroutine spmCheckAndCorrect(spm_in, spm_out, info)
    use iso_c_binding
    implicit none
    type(spmatrix_t),    intent(in),    target :: spm_in
    type(spmatrix_t),    intent(inout), target :: spm_out
    integer(kind=c_int), intent(out)           :: info

    info = spmCheckAndCorrect_c(c_loc(spm_in), c_loc(spm_out))
  end subroutine spmCheckAndCorrect

  subroutine spmGenRHS(type, nrhs, spm, x, ldx, b, ldb, info)
    use iso_c_binding
    implicit none
    integer(c_int),          intent(in)         :: type
    integer(kind=spm_int_t), intent(in)         :: nrhs
    type(spmatrix_t),        intent(in), target :: spm
    type(c_ptr),             intent(in), target :: x
    integer(kind=spm_int_t), intent(in)         :: ldx
    type(c_ptr),             intent(in), target :: b
    integer(kind=spm_int_t), intent(in)         :: ldb
    integer(kind=c_int),     intent(out)        :: info

    info = spmGenRHS_c(type, nrhs, c_loc(spm), x, ldx, b, ldb)
  end subroutine spmGenRHS

  subroutine spmCheckAxb(eps, nrhs, spm, x0, ldx0, b, ldb, x, ldx, info)
    use iso_c_binding
    implicit none
    real(kind=c_double),     intent(in)         :: eps
    integer(kind=spm_int_t), intent(in)         :: nrhs
    type(spmatrix_t),        intent(in), target :: spm
    type(c_ptr),             intent(in), target :: x0
    integer(kind=spm_int_t), intent(in)         :: ldx0
    type(c_ptr),             intent(in), target :: b
    integer(kind=spm_int_t), intent(in)         :: ldb
    type(c_ptr),             intent(in), target :: x
    integer(kind=spm_int_t), intent(in)         :: ldx
    integer(kind=c_int),     intent(out)        :: info

    info = spmCheckAxb_c(eps, nrhs, c_loc(spm), x0, ldx0, b, ldb, x, ldx)
  end subroutine spmCheckAxb

  subroutine spmIntConvert(n, input, value)
    use iso_c_binding
    implicit none
    integer(kind=spm_int_t), intent(in)             :: n
    integer(kind=c_int),     intent(inout), target  :: input
    integer(kind=spm_int_t), intent(out),   pointer :: value

    call c_f_pointer(spmIntConvert_c(n, c_loc(input)), value)
  end subroutine spmIntConvert

  subroutine spmLoad(spm, info)
    use iso_c_binding
    implicit none
    type(spmatrix_t),    intent(inout), target :: spm
    integer(kind=c_int), intent(out)           :: info

    info = spmLoad_c(c_loc(spm), c_null_ptr)
  end subroutine spmLoad

  subroutine spmSave(spm, info)
    use iso_c_binding
    implicit none
    type(spmatrix_t),    intent(in), target :: spm
    integer(kind=c_int), intent(out)        :: info

    info = spmSave_c(c_loc(spm), c_null_ptr)
  end subroutine spmSave

  subroutine spmReadDriver(driver, filename, spm, info)
    use iso_c_binding
    implicit none
    integer(c_int),         intent(in)            :: driver
    character(kind=c_char), intent(in),    target :: filename
    type(spmatrix_t),       intent(inout), target :: spm
    integer(kind=c_int),    intent(out)           :: info

    info = spmReadDriver_c(driver, c_loc(filename), c_loc(spm))
  end subroutine spmReadDriver

  subroutine spmParseLaplacianInfo(filename, flttype, dim1, dim2, dim3, alpha, &
       beta, info)
    use iso_c_binding
    implicit none
    character(kind=c_char),  intent(in),    target :: filename
    integer(c_int),          intent(inout), target :: flttype
    integer(kind=spm_int_t), intent(inout), target :: dim1
    integer(kind=spm_int_t), intent(inout), target :: dim2
    integer(kind=spm_int_t), intent(inout), target :: dim3
    real(kind=c_double),     intent(inout), target :: alpha
    real(kind=c_double),     intent(inout), target :: beta
    integer(kind=c_int),     intent(out)           :: info

    info = spmParseLaplacianInfo_c(c_loc(filename), c_loc(flttype), &
         c_loc(dim1), c_loc(dim2), c_loc(dim3), c_loc(alpha), c_loc(beta))
  end subroutine spmParseLaplacianInfo

  subroutine spm2Dense(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(in), target :: spm

    call spm2Dense_c(c_loc(spm))
  end subroutine spm2Dense

  subroutine spmPrint(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(in), target :: spm

    call spmPrint_c(c_loc(spm), c_null_ptr)
  end subroutine spmPrint

  subroutine spmPrintInfo(spm)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(in), target :: spm

    call spmPrintInfo_c(c_loc(spm), c_null_ptr)
  end subroutine spmPrintInfo

  subroutine spmExpand(spm, spmo)
    use iso_c_binding
    implicit none
    type(spmatrix_t), intent(in),  target  :: spm
    type(spmatrix_t), intent(out), pointer :: spmo

    call c_f_pointer(spmExpand_c(c_loc(spm)), spmo)
  end subroutine spmExpand

  subroutine spmDofExtend(spm, type, dof, spmo)
    use iso_c_binding
    implicit none
    type(spmatrix_t),    intent(in),  target  :: spm
    integer(kind=c_int), intent(in)           :: type
    integer(kind=c_int), intent(in)           :: dof
    type(spmatrix_t),    intent(out), pointer :: spmo

    call c_f_pointer(spmDofExtend_c(c_loc(spm), type, dof), spmo)
  end subroutine spmDofExtend


end module spmf
