!=========================================================================================
!Copyright (c) 2009-2019, The Regents of the University of Massachusetts, Amherst.
!E. Polizzi research lab
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification, 
!are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list of conditions 
!   and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
!   and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
!    products derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
!ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!==========================================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Include Sparse Primitives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!  include 'dzlsprim.f90' !! Sparse primitives



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED SPARSE INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

! List of routines:
!-------------------

!{D,Z}FEAST_{SCSR,HCSR,GCSR}{PEV}


!111111111111!!!!!! EXPERT ROUTINE (CORE) - using MKL-PARDISO
! dfeast_scsrpevx  (wrapper only)
! zfeast_hcsrpevx  (wrapper only)
! dfeast_gcsrpevx  (wrapper only)
! zfeast_gcsrpevx
! zfeast_scsrpevx 

!222222222222!!!!  DEFAULT ROUTINES (Wrappers to expert)
! dfeast_scsrpev
! zfeast_hcsrpev
! dfeast_gcsrpev
! zfeast_gcsrpev
! zfeast_scsrpev



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Remark: simple Wrapper to  zfeast_scsrpevx 
  !
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL SYMMETRIC SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        REAL DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                   
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  double precision,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  !------------------------------------------
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: zsa
  integer ::j,nnza
  logical :: mkl


  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to IFEAST 
     call difeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121147 ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! implement wrapping  
  allocate(zsa(maxval(isa(n+1,1:dmax+1))-1,dmax+1))
  do j=1,dmax+1
     nnza=isa(n+1,j)-1
     call ZLACP2( 'F', nnza, 1,sa(1,j) , nnza, zsa(1,j), nnza )
  enddo

  call zfeast_scsrpevx(UPLO,dmax,N,zsa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(zsa)

#endif

end subroutine dfeast_scsrpevx








subroutine zfeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Remark: simple Wrapper to  zfeast_gcsrpevx
  !
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX HERMITIAN SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        COMPLEX DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) or (2*M0) to store left eig. residuals if fpm(15)=0: Relative Residual of the solution (1-norm)
  !                                                                                  
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  complex(kind=(kind(1.0d0))),dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  !-----------------------------------------------------
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: fsa
  integer,dimension(:,:),allocatable :: fisa,fjsa
  integer ::j,nnza
  logical :: mkl


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to IFEAST
     call zifeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

#ifdef MKL 


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141247 ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((UPLO=='F').or.(UPLO=='f')) then
     call zfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  else
     !!  FORMAT conversion to full csr
     allocate(fisa(n+1,dmax+1))
     nnza=maxval(isa(n+1,1:dmax+1))-1
     nnza=2*nnza ! slightly overestimated
     allocate(fjsa(nnza,dmax+1))
     allocate(fsa(nnza,dmax+1))
     do j=1,dmax+1
        call zhcsr_convert_full(n,UPLO,sa(1,j),isa(1,j),jsa(1,j),fsa(1,j),fisa(1,j),fjsa(1,j))
     enddo

     call zfeast_gcsrpevx(dmax,N,fsa,fisa,fjsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     deallocate(fsa)
     deallocate(fisa)
     deallocate(fjsa)

  end if


#endif




end subroutine zfeast_hcsrpevx






subroutine dfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Remark: simple Wrapper to  zfeast_gcsrpevx
  !
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL  SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        REAL DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) or (2*M0) to store left eig. residuals if fpm(15)=0: Relative Residual of the solution (1-norm)
  !                                                                                  
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  double precision,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  ! !------------------------------------------
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: zsa
  integer ::j,nnza
  logical :: mkl


  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to IFEAST 
     call difeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121347 ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! implement wrapping  
  allocate(zsa(maxval(isa(n+1,1:dmax+1))-1,dmax+1))
  do j=1,dmax+1
     nnza=isa(n+1,j)-1
     call ZLACP2( 'F', nnza, 1,sa(1,j) , nnza, zsa(1,j), nnza )
  enddo

  call zfeast_gcsrpevx(dmax,N,zsa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(zsa)

#endif



end subroutine dfeast_gcsrpevx






subroutine zfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX  SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        COMPLEX DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) or (2*M0) to store left eig. residuals if fpm(15)=0: Relative Residual of the solution (1-norm)
  !                                                                                
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  complex(kind=(kind(1.0d0))),dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k,j
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,Aq,Bq
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  integer,dimension(:),allocatable :: isaz,jsaz,tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! matrix format
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: fsa
  double precision,dimension(:),allocatable :: ddiag
  integer :: opt,nnza,nnz
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:,:),allocatable :: cfsa !single precision copy
  complex,dimension(:,:),allocatable :: cwork,caux
!!!!!for pardiso
  integer(8),dimension(:,:),allocatable :: pt 
  integer,dimension(:,:),allocatable :: iparm
  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to IFEAST 
     call zifeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141347 ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!


  rank=0
  nb_procs=1
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
  !----------------------

  IF ( N<=0 ) THEN
     INFO = -102
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_GCSRPEVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT MATRIX CONVERSION 
!!! A and B already in full format
!!!!!!!!!!!! copy sa (needed for scaling)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! for all matrices A
  nnzA=maxval(isa(n+1,1:dmax+1))-1
  allocate(fsa(nnzA,dmax+1))
  do j=1,dmax+1
     call ZCOPY(isa(n+1,j)-1,sa(1,j),1, fsa(1,j), 1 )
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  endif


!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  isaz(1:n+1)=isa(1:n+1,1) ! initialize to A[1]
  nnz=isaz(n+1)-1 
  allocate(jsaz(nnz))
  jsaz(1:nnz)=jsa(1:nnz,1) ! initialize to A[1]

  allocate(saz(1,1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(n+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:n+1)=isaz(1:n+1)
     tempjsaz(1:nnz)=jsaz(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(N,N,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,sa(1,j),isa(1,j),jsa(1,j),saz,isaz,jsaz)
     deallocate(jsaz)
     nnz=isaz(n+1)-1 
     allocate(jsaz(nnz))
     opt=2 ! get jsaz
     call zaddcsr(N,N,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,sa(1,j),isa(1,j),jsa(1,j),saz,isaz,jsaz)
     deallocate(tempisaz)
     deallocate(tempjsaz)
  end do

  deallocate(tempsaz) ! deallocate dummy
  deallocate(saz) ! deallocate dummy       
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     Ze=Emid
     saz(1:nnz,1)=ZZERO
     do j=1,dmax+1
        call zinc_addcsr(N,N,Ze**(j-1),sa(1,j),isa(1,j),jsa(1,j),saz(1,1),isaz,jsaz)
     enddo

     allocate(ddiag(N))
     ddiag(1:N)=DONE
     ! extract diagonal 
     do i=1,N
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix
           if ((jsaz(k)==i).and.(abs(saz(k,1))/=DZERO)) ddiag(i)=abs(saz(k,1)) 
        enddo
     enddo
     !scale matrices A 
     do j=1,dmax+1
        do i=1,N
           do k=isa(i,j),isa(i+1,j)-1   
              fsa(k,j)=fsa(k,j)/(sqrt(ddiag(i))*sqrt(ddiag(jsa(k,j))))
           enddo
        end do
     end do

     deallocate(saz) ! deallocate dummy
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE RHS if fpm(5)=1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((fpm(41)==1).and.(fpm(5)==1)) then
     do i=1,N
        X(i,1:M0)=X(i,1:M0)*sqrt(ddiag(i))
        if (fpm(15)==0) X(i,M0+1:2*M0)=X(i,M0+1:2*M0)*sqrt(ddiag(i))
     enddo
  end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision   
     allocate(cfsa(nnzA,dmax+1)) ! overestimation
     do j=1,dmax+1 
        nnzA= isa(n+1,j)-1 ! real size
        call ZLAG2C(nnzA,1,fsa(1,j),nnzA,cfsa(1,j),nnzA,infoloc1)
     end do
     allocate(csaz(nnz,nfact)) !! Az matrix with potential for multiple factorizations
  else ! double precision
     allocate(saz(nnz,nfact))  !! Az matrix with potential for multiple factorizations
  end if


  if (infoloc1/=0) then
     info=-1
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!  use same factorization for (normal+transpose solve)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1 
  MNUM=1
  MTYPE=13      ! complex and unsymmetric

  !! check if matrix Hermitian
  if (mod(fpm(30),1000)/100==2) then
     MTYPE=3 ! complex and structurally symmetric
  end if

  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))
  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  end do
  !IPARM(2,:)=3 ! parallel omp nested dissection  (sensible to #threads- no consistency- same runs with different results)
  !IPARM(4)=11 !CGS solver
  IPARM(25,:)=1 ! parallel omp rhs solve
  IPARM(11,:)=0 ! disable scaling (taking care by feast fpm(41))
!!!!!!!!!!!!
  if (fpm(64)==1) then
     do i=1,64
        if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
     enddo
  endif
!!!!!!!!!!!!
  IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used
  MSGLVL=0!0 !0- no output, 1- output
  !PHASE=11 ! symbolic factorization (do it only once)

  if (fpm(42)==1) then ! mixed precision (single precision solver)       
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<
     allocate(caux(N,M0))
     allocate(cwork(N,M0))
  else ! double precision only
     allocate(zaux(N,M0))
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization 
  do while (ijob/=0)   
     call zfeast_grcipevx(ijob,dmax,N,Ze,work,zwork,Aq,Bq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    
     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize P(Ze)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization
        if (fpm(42)==1) then !single precision fact- csaz
csaz(1:nnz,id)=CZERO
           do j=1,dmax+1
              call cinc_addcsr(N,N,cmplx(Ze)**(j-1),cfsa(1,j),isa(1,j),jsa(1,j),csaz(1,id),isaz,jsaz)
           enddo !! get csaz

call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
        else ! full precision

           saz(1:nnz,id)=ZZERO      
           do j=1,dmax+1
              call zinc_addcsr(N,N,Ze**(j-1),fsa(1,j),isa(1,j),jsa(1,j),saz(1,id),isaz,jsaz)
           enddo !! get saz

           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=0 ! normal solve

        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else ! full precision
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     case(21) !!solve the linear system P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
        id=fpm(33)
        PHASE=33 ! solve          
        IPARM(12,id)=1 ! transpose conjugate solve

        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else ! full precision
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsa(1,fpm(57)),isa(1,fpm(57)),jsa(1,fpm(57)),X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(31) !! perform multiplication A^H(fpm(57))*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsa(1,fpm(57)),isa(1,fpm(57)),jsa(1,fpm(57)),X(1,fpm(34)),ZZERO,work(1,fpm(34)))

     end select
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (.not.((fpm(5)==1).and.(loop==0))) then ! at least one loop for fpm(5) (if not we do not enter case(10)-memory factorization- )
PHASE=-1 
  do i=1,nfact
     if (fpm(42)==1) then
        call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,infoloc2)
     else
        call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,infoloc2)
     endif
  end do
  if (infoloc2/=0) then
     info=-2
     return
  end if
  endif

  deallocate(PT)
  deallocate(iparm)

  deallocate(Aq)
  deallocate(Bq)
  deallocate(work)
  deallocate(zwork)

  deallocate(isaz)
  deallocate(jsaz)

  deallocate(fsa)


  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(cfsa)
  else
     deallocate(zaux)
     deallocate(saz)
  end if


  if (fpm(41)==1) then ! scale back the solution
     do i=1,N
        X(i,1:M0)=X(i,1:M0)/sqrt(ddiag(i))
        if (fpm(15)==0) X(i,M0+1:2*M0)=X(i,M0+1:2*M0)/sqrt(ddiag(i))
     enddo
     deallocate(ddiag)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! done with MKL PARDISO
#endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine zfeast_gcsrpevx





subroutine zfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX SYMMETRIC SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        COMPLEX DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                 
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  complex(kind=(kind(1.0d0))),dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k,j
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,Aq,Bq
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  integer,dimension(:),allocatable :: isaz,jsaz,tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! matrix format
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: usa
  double precision,dimension(:),allocatable :: ddiag
  integer,dimension(:,:), allocatable :: uisa,ujsa
  integer :: opt,nnza,nnz
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:,:),allocatable :: susa !single precision copy
  complex,dimension(:,:),allocatable :: cwork,caux
!!!!!for pardiso
  integer(8),dimension(:,:),allocatable :: pt 
  integer,dimension(:,:),allocatable :: iparm
  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to IFEAST 
     call zifeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141147 ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!


  rank=0
  nb_procs=1
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
  !----------------------

  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF ( N<=0 ) THEN
     INFO = -102
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_SCSRPEVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! for all matrices A
  nnzA=maxval(isa(n+1,1:dmax+1))-1
  if ((UPLO=='F').or.(UPLO=='f')) nnzA=nnzA/2+n ! slightly overestimated
  allocate(uisa(n+1,dmax+1))
  allocate(ujsa(nnzA,dmax+1))
  allocate(usa(nnzA,dmax+1))
  do j=1,dmax+1
     call zcsr_convert_upper(n,UPLO,sa(1,j),isa(1,j),jsa(1,j),usa(1,j),uisa(1,j),ujsa(1,j))
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  endif





!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  isaz(1:n+1)=uisa(1:n+1,1) ! initialize to A[1]
  nnz=isaz(n+1)-1 
  allocate(jsaz(nnz))
  jsaz(1:nnz)=ujsa(1:nnz,1) ! initialize to A[1]

  allocate(saz(1,1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(n+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:n+1)=isaz(1:n+1)
     tempjsaz(1:nnz)=jsaz(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(N,N,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,usa(1,j),uisa(1,j),ujsa(1,j),saz,isaz,jsaz)
     deallocate(jsaz)
     nnz=isaz(n+1)-1 
     allocate(jsaz(nnz))
     opt=2 ! get jsaz
     call zaddcsr(N,N,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,usa(1,j),uisa(1,j),ujsa(1,j),saz,isaz,jsaz)
     deallocate(tempisaz)
     deallocate(tempjsaz)
  end do

  deallocate(tempsaz) ! deallocate dummy
  deallocate(saz) ! deallocate dummy       
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     Ze=Emid
     saz(1:nnz,1)=ZZERO
     do j=1,dmax+1
        call zinc_addcsr(N,N,Ze**(j-1),usa(1,j),uisa(1,j),ujsa(1,j),saz(1,1),isaz,jsaz)
     enddo

     allocate(ddiag(N))
     ddiag(1:N)=DONE
     ! extract diagonal 
     do i=1,N
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix
           if ((jsaz(k)==i).and.(abs(saz(k,1))/=DZERO)) ddiag(i)=abs(saz(k,1)) 
        enddo
     enddo
     !scale matrices A 
     do j=1,dmax+1
        do i=1,N
           do k=uisa(i,j),uisa(i+1,j)-1   
              usa(k,j)=usa(k,j)/(sqrt(ddiag(i))*sqrt(ddiag(ujsa(k,j))))
           enddo
        end do
     end do

     deallocate(saz) ! deallocate dummy

  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE RHS if fpm(5)=1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((fpm(41)==1).and.(fpm(5)==1)) then
     do i=1,N
        X(i,1:M0)=X(i,1:M0)*sqrt(ddiag(i)) 
     enddo
  end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision   
     allocate(susa(nnzA,dmax+1)) ! overestimation
     do j=1,dmax+1 
        nnzA= uisa(n+1,j)-1 ! real size
        call ZLAG2C(nnzA,1,usa(1,j),nnzA,susa(1,j),nnzA,infoloc1)
     end do
     allocate(csaz(nnz,nfact)) !! Az matrix with potential for multiple factorizations
  else ! double precision
     allocate(saz(nnz,nfact))  !! Az matrix with potential for multiple factorizations

  end if


  if (infoloc1/=0) then
     info=-1
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  allocate(work(N,M0))
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!  use same factorization for (normal+transpose solve)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1 
  MNUM=1
  MTYPE=6      ! complex and symmetric
  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))
  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  end do
  !IPARM(2,:)=3 ! parallel omp nested dissection (sensible to #threads- no consistency- same runs with different results)
  !IPARM(4)=11 !CGS solver
  IPARM(25,:)=1 ! parallel omp rhs solve
  IPARM(11,:)=0 ! disable scaling (taking care by feast fpm(41))
!!!!!!!!!!!!
  if (fpm(64)==1) then
     do i=1,64
        if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
     enddo
  endif
!!!!!!!!!!!!
  IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used
  MSGLVL=0!0 !0- no output, 1- output
  !PHASE=11 ! symbolic factorization (do it only once)

  if (fpm(42)==1) then ! mixed precision (single precision solver)       
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<
     allocate(caux(N,M0))
     allocate(cwork(N,M0))
  else ! double precision only
     allocate(zaux(N,M0))
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization 
  do while (ijob/=0)   
     call zfeast_srcipevx(ijob,dmax,N,Ze,work,zwork,Aq,Bq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize P(Ze)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization
        if (fpm(42)==1) then !single precision fact- csaz

           csaz(1:nnz,id)=CZERO      
           do j=1,dmax+1
              call cinc_addcsr(N,N,cmplx(Ze)**(j-1),susa(1,j),uisa(1,j),ujsa(1,j),csaz(1,id),isaz,jsaz)
           enddo !! get csaz


           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
        else ! full precision

           saz(1:nnz,id)=ZZERO      
           do j=1,dmax+1
              call zinc_addcsr(N,N,Ze**(j-1),usa(1,j),uisa(1,j),ujsa(1,j),saz(1,id),isaz,jsaz)
           enddo !! get saz

           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve

        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else ! full precision
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('U','N',N,N,fpm(25),ZONE,usa(1,fpm(57)),uisa(1,fpm(57)),ujsa(1,fpm(57)),X(1,fpm(24)),ZZERO,work(1,fpm(24)))

     end select
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (.not.((fpm(5)==1).and.(loop==0))) then ! at least one loop for fpm(5) (if not we do not enter case(10)-memory factorization- )
PHASE=-1 
  do i=1,nfact
     if (fpm(42)==1) then
        call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,infoloc2)
     else
        call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,infoloc2)
     endif
  end do
  if (infoloc2/=0) then
     info=-2
     return
  end if
  endif

  deallocate(PT)
  deallocate(iparm)

  deallocate(Aq)
  deallocate(Bq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

  deallocate(usa)
  deallocate(uisa)
  deallocate(ujsa)


  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(susa)
  else
     deallocate(zaux)
     deallocate(saz)
  end if


  if (fpm(41)==1) then ! scale back the solution
     do i=1,N
        X(i,1:M0)=X(i,1:M0)/sqrt(ddiag(i))
     enddo
     deallocate(ddiag)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! done with MKL PARDISO
#endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine zfeast_scsrpevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine dfeast_scsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL SYMMETRIC SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        REAL DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                   
  !  info       (output)       INTEGER: Error handling (0: successful exit)          
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  double precision,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r 
  double precision :: epsout
  integer :: loop
  integer :: M0
  complex(kind=(kind(1.0d0))),dimension(*):: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_scsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121146
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine dfeast_scsrpev







subroutine zfeast_hcsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX HERMITIAN SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        COMPLEX DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) or (2*M0) to store left eig. residuals if fpm(15)=0: Relative Residual of the solution (1-norm)
  !                                                                                  
  !  info       (output)       INTEGER: Error handling (0: successful exit)          
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  complex(kind=(kind(1.0d0))),dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r 
  double precision :: epsout
  integer :: loop
  integer :: M0
  complex(kind=(kind(1.0d0))),dimension(*):: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_scsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141246
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)


end subroutine zfeast_hcsrpev








subroutine dfeast_gcsrpev(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL  SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        REAL DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) or (2*M0) to store left eig. residuals if fpm(15)=0: Relative Residual of the solution (1-norm)
  !                                                                                
  !  info       (output)       INTEGER: Error handling (0: successful exit)   
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  double precision,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121346
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine dfeast_gcsrpev







subroutine zfeast_gcsrpev(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX  SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        COMPLEX DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) or (2*M0) to store left eig. residuals if fpm(15)=0: Relative Residual of the solution (1-norm)
  !                                                                                  
  !  info       (output)       INTEGER: Error handling (0: successful exit)   
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  complex(kind=(kind(1.0d0))),dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_gcsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141346
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zfeast_gcsrpev


subroutine zfeast_scsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX SYMMETRIC SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  ! 
  !  Arguments
  !  =========
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  sa         (input)        COMPLEX DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid  (input)        COMPLEX DOUBLE PRECISION: search interval
  !   r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                
  !  info       (output)       INTEGER: Error handling (0: successful exit)          
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  complex(kind=(kind(1.0d0))),dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r 
  double precision :: epsout
  integer :: loop
  integer :: M0
  complex(kind=(kind(1.0d0))),dimension(*):: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_scsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141146
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)


end subroutine zfeast_scsrpev





