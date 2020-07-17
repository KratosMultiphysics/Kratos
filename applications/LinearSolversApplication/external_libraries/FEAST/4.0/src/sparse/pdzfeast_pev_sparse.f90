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
!include 'dzlsprim.f90' !! Sparse primitives
!include 'sclsprim.f90' !! for mixed precision


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED SPARSE INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

! List of routines:
!-------------------

!P{D,Z}FEAST_{SCSR,HCSR,GCSR}{PEV}


!111111111111!!!!!! EXPERT ROUTINE (CORE) - using MKL-CLUSTER-PARDISO
! pdfeast_scsrpevx  (wrapper only)
! pzfeast_hcsrpevx  
! pdfeast_gcsrpevx  (wrapper only)
! pzfeast_gcsrpevx
! pzfeast_scsrpevx 

!222222222222!!!!  DEFAULT ROUTINES (Wrappers to expert)
! pdfeast_scsrpev
! pzfeast_hcsrpev
! pdfeast_gcsrpev
! pzfeast_gcsrpev
! pzfeast_scsrpev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine pdfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Remark: simple Wrapper to  pzfeast_scsrpevx 
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
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !  include 'f90_noruntime_interface.fi'
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
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
  complex(kind=(kind(1.0d0))),dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
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

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
      call pdifeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=221147 ! code name
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

  call pzfeast_scsrpevx(UPLO,dmax,N,zsa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(zsa)

#endif

end subroutine pdfeast_scsrpevx







subroutine pzfeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !  sa         (input)        COMPLEX HERMITIAN DOUBLE PRECISION (nnzmax,dmax+1):  Matrices A- CSR format with nnzmax=max of nnz (non-zero elements) among all A sparse matrices
  !  isa        (input)        INTEGER(N+1,dmax+1): CSR row array of Matrices A
  !  jsa        (input)        INTEGER(nnzmax,dmax+1): CSR column array of Matrices A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output)        INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
  !                                                                                 
  !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
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
  complex(kind=(kind(1.0d0))),dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k,j,M02,jj
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,Aq,Bq,zaux,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
   complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  integer,dimension(:),allocatable :: isaz,jsaz,tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa
  integer,dimension(:), allocatable :: fisa,fjsa
  integer :: opt,nnza,nnz
!!!!!for cluster pardiso
  integer(8),dimension(:,:),allocatable :: pt
  !integer,dimension(64) :: iparm
  integer,dimension(:,:),allocatable :: iparm

  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
type ccsr_value
     complex,dimension(:),allocatable :: sa 
  end type ccsr_value
  type(ccsr_value), dimension(:),allocatable :: cmatAj 
  !complex,dimension(:),allocatable :: cfsa !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:,:),allocatable :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(:),allocatable :: matAj
  type(zcsr), dimension(:,:,:),allocatable :: matAjb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
     call pzifeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=241247   ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!
  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1    

  !----------------------------------------------
  L2_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm(49)
  call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
  call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  !----------------------

  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF( N<=0 ) THEN
     INFO = -102
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'PZFEAST_HCSRPEVX', -INFO+100 )
     RETURN
  END IF


  infoloc1=0
  infoloc2=0

!!!!!!!! Left eigenvectors computed?
  M02=M0 ! No
  if (fpm(15)==0) M02=2*M0 ! yes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Find which distribution scheme is used for the input matrix/rhs: global or local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call pfeast_distribution_type(N,isa(1,1),jsa(1,1),fpm(49),distype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
   
     Ntotal=N
!!! FORMAT CONVERSION TO FULL-CSR and DISTRIBUTION
    allocate(Nsize(nb_procs3))
     allocate(matAj(dmax+1))
     allocate(fisa(n+1))
     
 do j=1,dmax+1    
     nnza=isa(n+1,j)-1
     if ((UPLO/='F').and.(UPLO/='f')) nnza=2*nnza ! slightly overestimated
     allocate(fsa(nnza))
     allocate(fjsa(nnza))
     call   zhcsr_convert_full(n,UPLO,sa(1,j),isa(1,j),jsa(1,j),fsa,fisa,fjsa)
!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     call zcsr_distribute_row(Ntotal,fsa,fisa,fjsa,matAj(j),startj,endj,Nlocal,Nsize,fpm(49)) !! A
     deallocate(fsa)
     deallocate(fjsa)
  end do

  deallocate(fisa)
!!! DISTRIBUTE the RHS by rows among all L3 procs

     allocate(Xj(Nlocal,M02))
     Xj(1:Nlocal,1:M02)=X(startj:endj,1:M02) !!<< blas?

  end if !! global distype=0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for distributed user matrices !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (distype==1) then!! matrix A, and B locally known by row

     Nlocal=N ! since N local
     Ntotal=0
     call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code) ! Ntotal

     !! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)

     !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank3
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

!!! COPY the RHS (for name consistency)
     allocate(Xj(Nlocal,M02))
     if (fpm(5)==1) Xj(1:Nlocal,1:M02)=X(1:Nlocal,1:M02) !!<< blas?

!!! FORM the row distributed full matrix for cluster Pardiso   
     ! first simple copy
     ! the case F is then taking care of
     ! the cases L and U use this temp copy and will be modified to
     ! full format once the local blocks csr are defined below

     ! simple copy
      do j=1,dmax+1
        nnzA=isa(n+1,j)-1
        allocate(matAj(j)%isa(N+1))
        allocate(matAj(j)%jsa(nnzA))
        allocate(matAj(j)%sa(nnzA))
        call ZCOPY(nnza,sa(1,j),1, matAj(j)%sa, 1 )
        matAj(j)%isa(1:n+1)= isa(1:n+1,j)
        matAj(j)%jsa(1:nnzA)= jsa(1:nnzA,j)
        matAj(j)%n=N
        matAj(j)%nnz=nnzA
     end do
    
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank2+1,fpm(8),nb_procs2
        nfact=nfact+1
     end do
  endif




!!!!!!! Set up for Az matrix (isaz and jsaz)
  allocate(isaz(nlocal+1))
  isaz(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=isaz(nlocal+1)-1 
  allocate(jsaz(nnz))
  jsaz(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(saz(1,1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=isaz(1:nlocal+1)
     tempjsaz(1:nnz)=jsaz(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)     
     deallocate(jsaz)
     nnz=isaz(nlocal+1)-1 
     allocate(jsaz(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)  
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
        call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa(1),matAj(j)%isa(1),matAj(j)%jsa(1),saz(1,1),isaz,jsaz)
     enddo

     allocate(ddiag(Ntotal))
     ddiag(1:Ntotal)=DZERO
     ! extract diagonal locally
     do i=1,Nlocal
        ddiag(i+startj-1)=DONE
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix   
           if ((jsaz(k)==startj+i-1).and.(abs(saz(k,1))/=DZERO))  ddiag(i+startj-1)=abs(saz(k,1))
        enddo
     enddo
     call MPI_ALLREDUCE(MPI_IN_PLACE,ddiag,Ntotal,MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 

     !scale matrices A
     do j=1,dmax+1
     do i=1,Nlocal
        do k=matAj(j)%isa(i),matAj(j)%isa(i+1)-1   
           matAj(j)%sa(k)=matAj(j)%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj(j)%jsa(k))))
        enddo
     end do
     end do
     deallocate(saz) ! deallocate dummy

  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE RHS if fpm(5)=1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((fpm(41)==1).and.(fpm(5)==1)) then
     do i=1,Nlocal
        Xj(i,1:M02)=Xj(i,1:M02)*sqrt(ddiag(startj+i-1))
     enddo
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! DISTRIBUTE by Row/Column blocks (for mpi-matvec)-- prerequisite: matrix distributed by rows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  UPLO2='F' ! full csr local blocks by default
  allocate(mapA(nb_procs3,nb_procs3,dmax+1))
  allocate(matAjb(nb_procs3,nb_procs3,dmax+1))
   do j=1,dmax+1
  call zget_local_block_csr(Nsize,mapA(1,1,j),matAj(j),matAjb(1,1,j),fpm(49),nb_procs3)
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! HANDLE THE PARTICULAR CASE of user distributed matrices via Lower or upper format -- we need to transform them into full row csr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (distype==1) then
     if ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u') ) then 

        do j=1,dmax+1
 !! transform into full local block csr format        
        call zhlbcsr_uplo_to_csr(Nsize,mapA(1,1,j),matAjb(1,1,j),L3_COMM_WORLD,nb_procs3)
 !! convert into row distributed csr format
        deallocate(matAj(j)%jsa)
        deallocate(matAj(j)%sa)
        deallocate(matAj(j)%isa)
        call zlbcsr_distribute_row(Nsize,mapA(1,1,j),matAjb(1,1,j),matAj(j),L3_COMM_WORLD,nb_procs3)
     end do


!!! we need to reconstruct isaz and jsaz, as full local rows
  !allocate(isaz(nlocal+1))
  isaz(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=isaz(nlocal+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  jsaz(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(saz(1,1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=isaz(1:nlocal+1)
     tempjsaz(1:nnz)=jsaz(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)     
     deallocate(jsaz)
     nnz=isaz(nlocal+1)-1 
     allocate(jsaz(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)  
     deallocate(tempisaz)
     deallocate(tempjsaz)
  end do

  deallocate(tempsaz) ! deallocate dummy
  deallocate(saz) ! deallocate dummy   

end if
end if



     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
  ! copy matrix single precision
     allocate(cmatAj(dmax+1))
 do j=1,dmax+1
        nnza=matAj(j)%isa(Nlocal+1)-1
     allocate(cmatAj(j)%sa(nnza))
     call ZLAG2C(nnza,1,matAj(j)%sa,nnzA,cmatAj(j)%sa,nnzA,infoloc1)
  end do
     allocate(csaz(nnz,nfact)) !! Az matrix with potential for multiple factorizations
  else ! double precision
     allocate(saz(nnz,nfact))  
  end if

  if (infoloc1/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  allocate(work(Nlocal,M02))
  allocate(zwork(Nlocal,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  cluster pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1
  MNUM=1
  MTYPE=13      ! complex and unsymmetric
 

  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))

  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  enddo

  IPARM(11,:)=0 ! disable scaling !done by feast fpm(41)
  ! IPARM(2,:)=10 ! distributed nested (found a bug for 4x4 matrix, 2mpi)
  if (nb_procs3==1) then ! fall into pardiso
     !IPARM(2,:)=3 ! parallel omp nested dissection  (sensible to #threads- no consistency- same runs with different results)
     IPARM(25,:)=1 ! parallel omp rhs solve
  endif


!!!!!!!!!!!!
  if (fpm(64)==1) then
     do i=1,64
        if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
     enddo
  endif
!!!!!!!!!!!!
  IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used
  MSGLVL=0 !0- no output, 1- output
  !PHASE=11 ! symbolic factorization (do it only once)

  IPARM(40,:)=2 ! matrix storage (distributed)
  IPARM(41,:)=startj ! start row for local storage of saz
  IPARM(42,:)=endj ! end row for local storage of saz
!!!!!!!!!!


  if (fpm(42)==1) then ! mixed precision (single precision solver)
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<
     allocate(caux(Nlocal,M0))
     allocate(cwork(Nlocal,M0))

  else   !double precision   
     allocate(zaux(Nlocal,M0))

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0)   
  call zfeast_grcipevx(ijob,dmax,Nlocal,Ze,work,zwork,Aq,Bq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     case(10) !! form and factorize P(Ze)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3
        PHASE=12 !include the symbolic factorization
        if (fpm(42)==1) then !single precision fact- csaz
           
 csaz(1:nnz,id)=CZERO      
           do j=1,dmax+1
              call cinc_addcsr(Nlocal,Ntotal,cmplx(Ze)**(j-1),cmatAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,csaz(1,id),isaz,jsaz)
           enddo !! get csaz
           

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if


        else ! double precision fact- saz
           
 saz(1:nnz,id)=ZZERO      
           do j=1,dmax+1
              call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz(1,id),isaz,jsaz)
           enddo !! get saz
           
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
           end if

        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=0 ! normal solve
        if (fpm(42)==1) then !single precision solve
           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              !print *,'enter'
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
              !print *,'leave'
           end if

           call CLAG2Z(Nlocal, fpm(23), cwork, Nlocal, zwork, Nlocal, infoloc1)

        else ! double precision solve

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)

           end if
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=1 ! transpose conjugate solve
        !print *,21,id,rank3
        if (fpm(42)==1) then !single precision solve
           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           endif


           call CLAG2Z(Nlocal, fpm(23), cwork, Nlocal, zwork, Nlocal, infoloc1)

        else ! double precision solve
           if (nb_procs3>1) then

              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)


           end if


        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if


     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
        call zlbcsrmm(UPLO2,'N',Nsize,mapA(1,1,fpm(57)),fpm(25),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(31) !! perform multiplication A(fpm(57)^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapA(1,1,fpm(57)),fpm(35),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)


     end select
  end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!

if (.not.((fpm(5)==1).and.(loop==0))) then ! at least one loop for fpm(5) (if not we do not enter case(10)-memory factorization- )
PHASE=-1 
  do i=1,nfact
     if (fpm(42)==1) then
        if (nb_procs3>1) then
           call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
        else
           call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,infoloc2)
        end if
     else
        if (nb_procs3>1) then
           call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
        else
           call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,infoloc2)
        end if

     end if
  enddo
  if (infoloc2/=0) then
     info=-2
     return
  end if
  endif

  deallocate(PT)
deallocate(IPARM)
  !!  scale back the solution
  if (fpm(41)==1) then
     do i=1,Nlocal
        Xj(i,1:M02)=Xj(i,1:M02)/sqrt(ddiag(startj+i-1))
     enddo
     deallocate(ddiag)  
  endif
  !! put back result in all processors

  if (distype==0) then
     X(1:Ntotal,1:M02)=0.0d0
     X(startj:endj,1:M02)=Xj(1:Nlocal,1:M02)
     call MPI_ALLREDUCE(MPI_IN_PLACE,X,Ntotal*M02,MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
  elseif (distype==1) then
     X(1:Nlocal,1:M02)=Xj(1:Nlocal,1:M02) !!<< blas?
  endif
  deallocate(Nsize)
  deallocate(Xj)
do j=1,dmax+1
  deallocate(matAj(j)%isa)
  deallocate(matAj(j)%jsa)
  deallocate(matAj(j)%sa)
  end do
  deallocate(matAj)
  deallocate(mapA)
  deallocate(matAjb)

  deallocate(Aq)
  deallocate(Bq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

 ! if (distype==0) then
 !    deallocate(fsa)
 !    deallocate(fsb)
 ! end if
  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
      do j=1,dmax+1
     deallocate(cmatAj(j)%sa)
  end do
  deallocate(cmatAj) 
  else
     deallocate(zaux)
     deallocate(saz)
  end if


#endif
!!!! done with MKL CLUSTER PARDISO!!!! done with MKL CLUSTER PARDISO

end subroutine pzfeast_hcsrpevx









subroutine pdfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !    Remark: simple Wrapper to  pzfeast_gcsrpevx
  !
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL SPARSE FORMAT 
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
  !  fpm (input/output)        INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
  !                                                                                 
  !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
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
  complex(kind=(kind(1.0d0))),dimension(*):: Zne,Wne

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
      call pdifeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=221347 ! code name
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

  call pzfeast_gcsrpevx(dmax,N,zsa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(zsa)

#endif

end subroutine pdfeast_gcsrpevx


  



subroutine pzfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !  fpm (input/output)        INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
  !                                                                                  
  !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
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
  complex(kind=(kind(1.0d0))),dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k,j,M02,jj
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,Aq,Bq,zaux,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
   complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  integer,dimension(:),allocatable :: isaz,jsaz,tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
!!!! matrix format
!  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa
  integer :: opt,nnza,nnz
!!!!!for cluster pardiso
  integer(8),dimension(:,:),allocatable :: pt
  !integer,dimension(64) :: iparm
  integer,dimension(:,:),allocatable :: iparm

  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
type ccsr_value
     complex,dimension(:),allocatable :: sa 
  end type ccsr_value
  type(ccsr_value), dimension(:),allocatable :: cmatAj 
  !complex,dimension(:),allocatable :: cfsa !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:,:),allocatable :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(:),allocatable :: matAj
  type(zcsr), dimension(:,:,:),allocatable :: matAjb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
     call pzifeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=241347   ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!
  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1    

  !----------------------------------------------
  L2_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm(49)
  call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
  call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  !----------------------


  IF( N<=0 ) THEN
     INFO = -102
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'PZFEAST_GCSRPEVX', -INFO+100 )
     RETURN
  END IF


  infoloc1=0
  infoloc2=0

!!!!!!!! Left eigenvectors computed?
  M02=M0 ! No
  if (fpm(15)==0) M02=2*M0 ! yes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Find which distribution scheme is used for the input matrix/rhs: global or local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call pfeast_distribution_type(N,isa(1,1),jsa(1,1),fpm(49),distype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!!! Perform a copy of the matrix values
    ! nnza=isa(n+1)-1
    ! allocate(fsa(nnza))
    ! call ZCOPY(nnza,sa,1, fsa, 1 )
    ! nnzb=isb(n+1)-1
    ! allocate(fsb(nnzb))
    ! call ZCOPY( nnzb,sb,1, fsb, 1 )

!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
      allocate(matAj(dmax+1))
  do j=1,dmax+1
     call zcsr_distribute_row(Ntotal,sa(1,j),isa(1,j),jsa(1,j),matAj(j),startj,endj,Nlocal,Nsize,fpm(49)) !! A
  end do
!!! DISTRIBUTE the RHS by rows among all L3 procs

     allocate(Xj(Nlocal,M02))
     Xj(1:Nlocal,1:M02)=X(startj:endj,1:M02) !!<< blas?

     !print *,rank3,N,Nlocal
     !info=-99
     !return

  end if !! global distype=0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for distributed user matrices !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (distype==1) then!! matrix A, and B locally known by row

     Nlocal=N ! since N local
     Ntotal=0
     call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code) ! Ntotal

     !! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)

     !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank3
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

!!! COPY the RHS (for name consistency)
     allocate(Xj(Nlocal,M02))
     if (fpm(5)==1) Xj(1:Nlocal,1:M02)=X(1:Nlocal,1:M02) !!<< blas?

!!! FORM the row distributed full matrix for cluster Pardiso
     ! simple copy
      do j=1,dmax+1
        nnzA=isa(n+1,j)-1
        allocate(matAj(j)%isa(N+1))
        allocate(matAj(j)%jsa(nnzA))
        allocate(matAj(j)%sa(nnzA))
        call ZCOPY(nnza,sa(1,j),1, matAj(j)%sa, 1 )
        matAj(j)%isa(1:n+1)= isa(1:n+1,j)
        matAj(j)%jsa(1:nnzA)= jsa(1:nnzA,j)
        matAj(j)%n=N
        matAj(j)%nnz=nnzA
     end do
    

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank2+1,fpm(8),nb_procs2
        nfact=nfact+1
     end do
  endif




!!!!!!! Set up for Az matrix (isaz and jsaz)
  allocate(isaz(nlocal+1))
  isaz(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=isaz(nlocal+1)-1 
  allocate(jsaz(nnz))
  jsaz(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(saz(1,1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=isaz(1:nlocal+1)
     tempjsaz(1:nnz)=jsaz(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)     
     deallocate(jsaz)
     nnz=isaz(nlocal+1)-1 
     allocate(jsaz(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)  
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
        call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa(1),matAj(j)%isa(1),matAj(j)%jsa(1),saz(1,1),isaz,jsaz)
     enddo

     allocate(ddiag(Ntotal))
     ddiag(1:Ntotal)=DZERO
     ! extract diagonal locally
     do i=1,Nlocal
        ddiag(i+startj-1)=DONE
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix   
           if ((jsaz(k)==startj+i-1).and.(abs(saz(k,1))/=DZERO))  ddiag(i+startj-1)=abs(saz(k,1))
        enddo
     enddo
     call MPI_ALLREDUCE(MPI_IN_PLACE,ddiag,Ntotal,MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 

     !scale matrices A
     do j=1,dmax+1
     do i=1,Nlocal
        do k=matAj(j)%isa(i),matAj(j)%isa(i+1)-1   
           matAj(j)%sa(k)=matAj(j)%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj(j)%jsa(k))))
        enddo
     end do
     end do
     deallocate(saz) ! deallocate dummy

  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE RHS if fpm(5)=1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((fpm(41)==1).and.(fpm(5)==1)) then
     do i=1,Nlocal
        Xj(i,1:M02)=Xj(i,1:M02)*sqrt(ddiag(startj+i-1))
     enddo
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! DISTRIBUTE by Row/Column blocks (for mpi-matvec)-- prerequisite: matrix distributed by rows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  UPLO2='F' ! full csr local blocks by default
  allocate(mapA(nb_procs3,nb_procs3,dmax+1))
  allocate(matAjb(nb_procs3,nb_procs3,dmax+1))
   do j=1,dmax+1
  call zget_local_block_csr(Nsize,mapA(1,1,j),matAj(j),matAjb(1,1,j),fpm(49),nb_procs3)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision
     allocate(cmatAj(dmax+1))
 do j=1,dmax+1
        nnza=matAj(j)%isa(Nlocal+1)-1
     allocate(cmatAj(j)%sa(nnza))
     call ZLAG2C(nnza,1,matAj(j)%sa,nnzA,cmatAj(j)%sa,nnzA,infoloc1)
  end do
     allocate(csaz(nnz,nfact)) !! Az matrix with potential for multiple factorizations
  else ! double precision
     allocate(saz(nnz,nfact))  
  end if

  if (infoloc1/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  allocate(work(Nlocal,M02))
  allocate(zwork(Nlocal,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  cluster pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1
  MNUM=1
  MTYPE=13      ! complex and unsymmetric

 !! check if matrix Hermitian
!  if (mod(fpm(30),1000)/100==2) then
!     MTYPE=3 ! complex and structurally symmetric
!  end if

  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))

  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  enddo

  IPARM(11,:)=0 ! disable scaling !done by feast fpm(41)
  ! IPARM(2,:)=10 ! distributed nested (found a bug for 4x4 matrix, 2mpi)
  if (nb_procs3==1) then ! fall into pardiso
     !IPARM(2,:)=3 ! parallel omp nested dissection  (sensible to #threads- no consistency- same runs with different results)
     IPARM(25,:)=1 ! parallel omp rhs solve
  endif


!!!!!!!!!!!!
  if (fpm(64)==1) then
     do i=1,64
        if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
     enddo
  endif
!!!!!!!!!!!!
  IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used
  MSGLVL=0 !0- no output, 1- output
  !PHASE=11 ! symbolic factorization (do it only once)

  IPARM(40,:)=2 ! matrix storage (distributed)
  IPARM(41,:)=startj ! start row for local storage of saz
  IPARM(42,:)=endj ! end row for local storage of saz
!!!!!!!!!!


  if (fpm(42)==1) then ! mixed precision (single precision solver)
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<
     allocate(caux(Nlocal,M0))
     allocate(cwork(Nlocal,M0))

  else   !double precision   
     allocate(zaux(Nlocal,M0))

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0)   
  call zfeast_grcipevx(ijob,dmax,Nlocal,Ze,work,zwork,Aq,Bq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     case(10) !! form and factorize P(Ze)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3
        PHASE=12 !include the symbolic factorization
        if (fpm(42)==1) then !single precision fact- csaz
           
 csaz(1:nnz,id)=CZERO      
           do j=1,dmax+1
              call cinc_addcsr(Nlocal,Ntotal,cmplx(Ze)**(j-1),cmatAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,csaz(1,id),isaz,jsaz)
           enddo !! get csaz
           

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if


        else ! double precision fact- saz
           
 saz(1:nnz,id)=ZZERO      
           do j=1,dmax+1
              call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz(1,id),isaz,jsaz)
           enddo !! get saz
           
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
           end if

        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=0 ! normal solve
        if (fpm(42)==1) then !single precision solve
           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              !print *,'enter'
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
              !print *,'leave'
           end if

           call CLAG2Z(Nlocal, fpm(23), cwork, Nlocal, zwork, Nlocal, infoloc1)

        else ! double precision solve

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)

           end if
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=1 ! transpose conjugate solve
        !print *,21,id,rank3
        if (fpm(42)==1) then !single precision solve
           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           endif


           call CLAG2Z(Nlocal, fpm(23), cwork, Nlocal, zwork, Nlocal, infoloc1)

        else ! double precision solve
           if (nb_procs3>1) then

              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)


           end if


        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if


     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
        call zlbcsrmm(UPLO2,'N',Nsize,mapA(1,1,fpm(57)),fpm(25),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(31) !! perform multiplication A(fpm(57)^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapA(1,1,fpm(57)),fpm(35),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)


     end select
  end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!

if (.not.((fpm(5)==1).and.(loop==0))) then ! at least one loop for fpm(5) (if not we do not enter case(10)-memory factorization- )
PHASE=-1 
  do i=1,nfact
     if (fpm(42)==1) then
        if (nb_procs3>1) then
           call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
        else
           call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,infoloc2)
        end if
     else
        if (nb_procs3>1) then
           call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
        else
           call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,infoloc2)
        end if

     end if
  enddo
  if (infoloc2/=0) then
     info=-2
     return
  end if
  endif

  deallocate(PT)
  deallocate(IPARM)

  !!  scale back the solution
  if (fpm(41)==1) then
     do i=1,Nlocal
        Xj(i,1:M02)=Xj(i,1:M02)/sqrt(ddiag(startj+i-1))
     enddo
     deallocate(ddiag)  
  endif

  !! put back result in all processors

  if (distype==0) then
     X(1:Ntotal,1:M02)=0.0d0
     X(startj:endj,1:M02)=Xj(1:Nlocal,1:M02)
     call MPI_ALLREDUCE(MPI_IN_PLACE,X,Ntotal*M02,MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
  elseif (distype==1) then
     X(1:Nlocal,1:M02)=Xj(1:Nlocal,1:M02) !!<< blas?
  endif

  deallocate(Nsize)
  deallocate(Xj)
do j=1,dmax+1
  deallocate(matAj(j)%isa)
  deallocate(matAj(j)%jsa)
  deallocate(matAj(j)%sa)
end do
deallocate(matAj)
  deallocate(mapA)
  deallocate(matAjb)


  deallocate(Aq)
  deallocate(Bq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

 ! if (distype==0) then
 !    deallocate(fsa)
 !    deallocate(fsb)
 ! end if

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
      do j=1,dmax+1
     deallocate(cmatAj(j)%sa)
  end do
  deallocate(cmatAj) 
  else
     deallocate(zaux)
     deallocate(saz)
  end if


#endif
!!!! done with MKL CLUSTER PARDISO!!!! done with MKL CLUSTER PARDISO

end subroutine pzfeast_gcsrpevx











subroutine pzfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !  include 'f90_noruntime_interface.fi'
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
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
  complex(kind=(kind(1.0d0))),dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k,j,jj
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,Aq,Bq,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  double precision, dimension(:),allocatable ::ddiag
  integer,dimension(:),allocatable :: isaz,jsaz,tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: usa
  integer,dimension(:), allocatable :: uisa,ujsa
  integer :: opt,nnza,nnz
!!!!!for cluster pardiso
  integer(8),dimension(:,:),allocatable :: pt 
  integer,dimension(:,:),allocatable :: iparm
  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  type ccsr_value
     complex,dimension(:),allocatable :: sa 
  end type ccsr_value
  !complex,dimension(:,:),allocatable :: susa !single precision copy
  type(ccsr_value), dimension(:),allocatable :: cmatAj 
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:,:),allocatable :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr), dimension(:),allocatable :: matAj
  type(zcsr), dimension(:,:,:),allocatable :: matAjb


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
      call pzifeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=241147  ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!


  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1    

  !----------------------------------------------
  L2_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm(49)
  call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
  call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  !----------------------


  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF ( N<=0 ) THEN
     INFO = -102
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'PZFEAST_SCSRPEVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Find which distribution scheme is used for the input matrix/rhs: global or local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call pfeast_distribution_type(N,isa(1,1),jsa(1,1),fpm(49),distype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!! FORMAT CONVERSION TO CSR-UPPER  and
!!! DISTRIBUTE the (upper) sparse matrix by rows among all L3 procs  
      allocate(Nsize(nb_procs3))
     allocate(matAj(dmax+1))
     allocate(uisa(n+1))
     do j=1,dmax+1
     nnzA=isa(n+1,j)-1
     if ((UPLO=='F').or.(UPLO=='f')) nnzA=nnzA/2+n ! slightly overestimated
  allocate(ujsa(nnzA))
  allocate(usa(nnzA))
 call zcsr_convert_upper(n,UPLO,sa(1,j),isa(1,j),jsa(1,j),usa,uisa,ujsa)
 call zcsr_distribute_row(Ntotal,usa,uisa,ujsa,matAj(j),startj,endj,Nlocal,Nsize,fpm(49)) !! Aj
 deallocate(usa)
 deallocate(ujsa)
  enddo  
deallocate(uisa)
!!! DISTRIBUTE the RHS by rows among all L3 procs
     allocate(Xj(Nlocal,M0))
     Xj(1:Nlocal,1:M0)=X(startj:endj,1:M0) !!<< blas?

     
  end if !! global distype=0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for distributed user matrices !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==1) then!! matrix A, and B locally known by row

     Nlocal=N ! since N local
     Ntotal=0
     call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code) ! Ntotal

     !! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)

     !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank3
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

!!! COPY the RHS (for name consistency)..!!<<try allocatable
     allocate(Xj(Nlocal,M0))
     if (fpm(5)==1) Xj(1:Nlocal,1:M0)=X(1:Nlocal,1:M0) !!<< blas?


     allocate(matAj(dmax+1))
!!! FORM the row distributed upper matrix for Pardiso   
     if ((UPLO=='F').or.(UPLO=='f')) then !matrix need to be reduced to upper part
        !! for A
        do j=1,dmax+1
        nnzA=isa(n+1,j)-1 ! overestimation
        allocate(matAj(j)%isa(N+1))
        allocate(matAj(j)%jsa(nnzA))
        allocate(matAj(j)%sa(nnzA))
        matAj(j)%isa(1)=1
        do i=1,N
           jj=0
           do k=isa(i,j),isa(i+1,j)-1
              if (jsa(k,j)>=i+startj-1) then
                 jj=jj+1
                 matAj(j)%jsa(matAj(j)%isa(i)+jj-1)=jsa(k,j)
                 matAj(j)%sa(matAj(j)%isa(i)+jj-1)=sa(k,j)
              end if
           enddo
           matAj(j)%isa(i+1)=matAj(j)%isa(i)+j 
        enddo
        matAj(j)%n=N
        matAj(j)%nnz=matAj(j)%isa(n+1)-1
     end do
       
     else ! include the following cases
        !! ((UPLO=='U').or.(UPLO=='u')) !! upper-csr already -- simple copy
        !! ((UPLO=='L').or.(UPLO=='l')) !! temp-copy !! transpose difficult..could be done once the local blocks csr are defined below

        !! for all matrices A
 do j=1,dmax+1
        nnzA=isa(n+1,j)-1
        allocate(matAj(j)%isa(N+1))
        allocate(matAj(j)%jsa(nnzA))
        allocate(matAj(j)%sa(nnzA))
        call ZCOPY(nnza,sa(1,j),1, matAj(j)%sa, 1 )
        matAj(j)%isa(1:n+1)= isa(1:n+1,j)
        matAj(j)%jsa(1:nnzA)= jsa(1:nnzA,j)
        matAj(j)%n=N
        matAj(j)%nnz=nnzA
     end do
     
  end if
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank2+1,fpm(8),nb_procs2
        nfact=nfact+1
     end do
  endif


!!!!!!! Set up for Az matrix (isaz and jsaz)
  allocate(isaz(nlocal+1))
  isaz(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=isaz(nlocal+1)-1 
  allocate(jsaz(nnz))
  jsaz(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(saz(1,1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=isaz(1:nlocal+1)
     tempjsaz(1:nnz)=jsaz(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)     
     deallocate(jsaz)
     nnz=isaz(nlocal+1)-1 
     allocate(jsaz(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)  
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
        call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa(1),matAj(j)%isa(1),matAj(j)%jsa(1),saz(1,1),isaz,jsaz)
     enddo

     allocate(ddiag(Ntotal))
     ddiag(1:Ntotal)=DZERO
     ! extract diagonal locally
     do i=1,Nlocal
        ddiag(i+startj-1)=DONE
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix   
           if ((jsaz(k)==startj+i-1).and.(abs(saz(k,1))/=DZERO))  ddiag(i+startj-1)=abs(saz(k,1))
        enddo
     enddo
     call MPI_ALLREDUCE(MPI_IN_PLACE,ddiag,Ntotal,MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 
     !scale matrices A
     do j=1,dmax+1
     do i=1,Nlocal
        do k=matAj(j)%isa(i),matAj(j)%isa(i+1)-1   
           matAj(j)%sa(k)=matAj(j)%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj(j)%jsa(k))))
        enddo
     end do
  end do
     deallocate(saz) ! deallocate dummy

  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE RHS if fpm(5)=1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((fpm(41)==1).and.(fpm(5)==1)) then
     do i=1,Nlocal
        Xj(i,1:M0)=Xj(i,1:M0)*sqrt(ddiag(startj+i-1))
     enddo
  end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! DISTRIBUTE by Row/Column blocks (for mpi-matvec)-- prerequisite: matrix distributed by rows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  UPLO2='U' ! Upper csr local blocks by default
  !  if (distype==1) then
  !     if ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  !  end if
  allocate(mapA(nb_procs3,nb_procs3,dmax+1))
  allocate(matAjb(nb_procs3,nb_procs3,dmax+1))
  do j=1,dmax+1
  call zget_local_block_csr(Nsize,mapA(1,1,j),matAj(j),matAjb(1,1,j),fpm(49),nb_procs3)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! HANDLE THE PARTICULAR CASE of user distributed matrices via Lower format (get it upper for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (distype==1) then
     if ((UPLO=='L').or.(UPLO=='l')) then !! perform the transpose
        !! create local block upper format
do j=1,dmax+1
        call zlbcsr_transpose(Nsize,mapA(1,1,j),matAjb(1,1,j),L3_COMM_WORLD,nb_procs3)
        !! convert into row distributed csr format
        deallocate(matAj(j)%jsa)
        deallocate(matAj(j)%sa)
        deallocate(matAj(j)%isa)
        call zlbcsr_distribute_row(Nsize,mapA(1,1,j),matAjb(1,1,j),matAj(j),L3_COMM_WORLD,nb_procs3)
 end do
        

!!! we need to reconstruct isaz and jsaz, as upper local rows
  !allocate(isaz(nlocal+1))
  isaz(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=isaz(nlocal+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  jsaz(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(saz(1,1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=isaz(1:nlocal+1)
     tempjsaz(1:nnz)=jsaz(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)     
     deallocate(jsaz)
     nnz=isaz(nlocal+1)-1 
     allocate(jsaz(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz,isaz,jsaz)  
     deallocate(tempisaz)
     deallocate(tempjsaz)
  end do

  deallocate(tempsaz) ! deallocate dummy
  deallocate(saz) ! deallocate dummy       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision
     allocate(cmatAj(dmax+1))
 do j=1,dmax+1
        nnza=matAj(j)%isa(Nlocal+1)-1
     allocate(cmatAj(j)%sa(nnza))
     call ZLAG2C(nnza,1,matAj(j)%sa,nnzA,cmatAj(j)%sa,nnzA,infoloc1)
  end do
     allocate(csaz(nnz,nfact)) !! Az matrix with potential for multiple factorizations
  else ! double precision
     allocate(saz(nnz,nfact))  
  end if

  if (infoloc1/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  allocate(work(Nlocal,M0))
  allocate(zwork(Nlocal,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  cluster pardiso initialization
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
  enddo
  IPARM(11,:)=0 ! disable scaling !done by feast fpm(41)
  !IPARM(2,:)=10 !distributed nested dissection! (found a bug for 4x4 matrix, 2mpi or 1mpi)
  if (nb_procs3==1) then ! fall into pardiso
     !IPARM(2,:)=3 ! parallel omp nested dissection  (sensible to #threads- no consistency- same runs with different results)
     IPARM(25,:)=1 ! parallel omp rhs solve
  endif
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

  IPARM(40,:)=2 ! matrix storage (distributed or centralized)
  IPARM(41,:)=startj ! start row for local storage of saz
  IPARM(42,:)=endj ! end row for local storage of saz
!!!!!!!!!!


  if (fpm(42)==1) then ! mixed precision (single precision solver)
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<   
     allocate(caux(Nlocal,M0))
     allocate(cwork(Nlocal,M0))
  else !double precision 

     allocate(zaux(Nlocal,M0))
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
 call zfeast_srcipevx(ijob,dmax,Nlocal,Ze,work,zwork,Aq,Bq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)    
!print *,'ijob',rank2,ijob
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
              call cinc_addcsr(Nlocal,Ntotal,cmplx(Ze)**(j-1),cmatAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,csaz(1,id),isaz,jsaz)
           enddo !! get csaz

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if

        else ! double precision fact- saz
         

  saz(1:nnz,id)=ZZERO      
           do j=1,dmax+1
              call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,saz(1,id),isaz,jsaz)
           enddo !! get saz

           
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
           end if
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


        if (fpm(42)==1) then !single precision solve
           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              !print *,rank2,id
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           endif

           call CLAG2Z(Nlocal, fpm(23), cwork, Nlocal, zwork, Nlocal, infoloc1)

        else ! double precision solve
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
           end if
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

        call zlbcsrmm(UPLO2,'N',Nsize,mapA(1,1,fpm(57)),fpm(25),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     end select
  end do
!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!

if (.not.((fpm(5)==1).and.(loop==0))) then ! at least one loop for fpm(5) (if not we do not enter case(10)-memory factorization- )
PHASE=-1 
  do i=1,nfact
     if (fpm(42)==1) then
        if (nb_procs3>1) then
           call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
        else
           call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,cwork,caux,infoloc2)
        end if
     else
        if (nb_procs3>1) then
           call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
        else
           call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,fpm(23),IPARM(1,i),MSGLVL,zwork,zaux,infoloc2)
        end if
     end if
  enddo
  if (infoloc2/=0) then
     info=-2
     return
  end if
  endif

  deallocate(PT)
  deallocate(IPARM)


  !!  scale back the solution
  if (fpm(41)==1) then
     do i=1,Nlocal
        Xj(i,1:M0)=Xj(i,1:M0)/sqrt(ddiag(startj+i-1))
     enddo
     deallocate(ddiag)  
  endif

  !! put back result in all processors

  if (distype==0) then
     X(1:Ntotal,1:M0)=0.0d0
     X(startj:endj,1:M0)=Xj(1:Nlocal,1:M0)
     call MPI_ALLREDUCE(MPI_IN_PLACE,X,Ntotal*M0,MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
  elseif (distype==1) then
     X(1:Nlocal,1:M0)=Xj(1:Nlocal,1:M0) !!<< blas?
  endif

  deallocate(Nsize)
  deallocate(Xj)
do j=1,dmax+1
  deallocate(matAj(j)%isa)
  deallocate(matAj(j)%jsa)
  deallocate(matAj(j)%sa)
end do
deallocate(matAj)
  deallocate(mapA)
  deallocate(matAjb)

  deallocate(Aq)
  deallocate(Bq)
  deallocate(work)
  deallocate(zwork)

  deallocate(isaz)
  deallocate(jsaz)

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
 do j=1,dmax+1
     deallocate(cmatAj(j)%sa)
  end do
  deallocate(cmatAj)  
  else
     deallocate(zaux)
     deallocate(saz)
  endif



#endif
!!!! done with MKL CLUSTER PARDISO

end subroutine pzfeast_scsrpevx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pdfeast_scsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
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
  !    include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,dmax
  integer,dimension(n+1,*) :: isa
  double precision,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: sa
  integer,dimension(maxval(isa(n+1,1:dmax+1))-1,*) :: jsa
  integer,dimension(*) :: fpm
  double precision :: epsout
   complex(kind=(kind(1.0d0))) :: Emid
  integer :: loop
  double precision :: r
  integer :: M0
    complex(kind=(kind(1.0d0))),dimension(*)  :: E
    complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pdfeast_scsrpevx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. 
#ifdef MKL
  mkl=.true.
#endif

  fpm(30)=221146
   if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pdfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdfeast_scsrpev




subroutine pzfeast_hcsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
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
  !     include 'f90_noruntime_interface.fi'
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
    complex(kind=(kind(1.0d0))),dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pzfeast_hcsrpevx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. 
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=241246
   if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
   allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzfeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzfeast_hcsrpev



subroutine pdfeast_gcsrpev(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output)        INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
  !                                                                                   
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !     include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: pdfeast_gcsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
   logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. 
#ifdef MKL
  mkl=.true.
#endif

  fpm(30)=221346
   if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pdfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdfeast_gcsrpev




subroutine pzfeast_gcsrpev(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !  fpm (input/output)        INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
  !                                                                                  
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !    include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: pzfeast_gcsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. 
#ifdef MKL
  mkl=.true.
#endif

  fpm(30)=241346
   if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzfeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzfeast_gcsrpev




subroutine pzfeast_scsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !  fpm (input/output)        INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input)  Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution 
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                  
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !    include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: pzfeast_scsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. 
#ifdef MKL
  mkl=.true.
#endif

  fpm(30)=241146
   if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzfeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pzfeast_scsrpev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



