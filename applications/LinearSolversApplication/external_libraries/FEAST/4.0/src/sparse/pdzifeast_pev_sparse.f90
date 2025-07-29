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

!P{D,Z}IFEAST_{SCSR,HCSR,GCSR}{PEV}


!111111111111!!!!!! EXPERT ROUTINE (CORE) - using IFEAST/Bicgstab
! pdifeast_scsrpevx  (wrapper only)
! pzifeast_hcsrpevx  
! pdifeast_gcsrpevx  (wrapper only)
! pzifeast_gcsrpevx
! pzifeast_scsrpevx 

!222222222222!!!!  DEFAULT ROUTINES (Wrappers to expert)
! pdifeast_scsrpev
! pzifeast_hcsrpev
! pdifeast_gcsrpev
! pzifeast_gcsrpev
! pzifeast_scsrpev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine pdifeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Remark: simple Wrapper to  pzifeast_scsrpevx 
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


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=222147 ! code name
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

  call pzifeast_scsrpevx(UPLO,dmax,N,zsa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(zsa)

end subroutine pdifeast_scsrpevx







subroutine pzifeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  integer,dimension(:),allocatable :: tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa
  integer,dimension(:), allocatable :: fisa,fjsa
  integer :: opt,nnza,nnz

!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout
  double precision :: DZNRM2

  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:,:),allocatable :: mapA
  integer,dimension(:,:),allocatable :: mapZ,mapC
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matZj
  type(zcsr),dimension(:),allocatable :: matAj
  type(zcsr), dimension(:,:,:),allocatable :: matAjb,matZjb
  !! Az matrix single + copy
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr) :: matCj
  type(ccsr), dimension(:,:,:),allocatable :: matCjb
  type(ccsr), dimension(:,:,:),allocatable :: smatAjb 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=242247   ! code name
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
     CALL XERBLA( 'PZIFEAST_HCSRPEVX', -INFO+100 )
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

!!! FORM the row distributed full matrix 
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
  allocate(matZj%isa(nlocal+1))
  matZj%isa(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=matZj%isa(nlocal+1)-1 
  allocate(matZj%jsa(nnz))
  matZj%jsa(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(matZj%sa(1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=matZj%isa(1:nlocal+1)
     tempjsaz(1:nnz)=matZj%jsa(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)     
     deallocate(matZj%jsa)
     nnz=matZj%isa(nlocal+1)-1 
     allocate(matZj%jsa(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)  
     deallocate(tempisaz)
     deallocate(tempjsaz)
  end do

  deallocate(tempsaz) ! deallocate dummy
  deallocate(matZj%sa)
  allocate(matZj%sa(nnz))
  matZj%sa(1:nnz)= ZZERO ! initialize
  matZj%n=Nlocal
  matZj%nnz=nnz     
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     Ze=Emid

     do j=1,dmax+1
        call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa(1),matAj(j)%isa(1),matAj(j)%jsa(1),matZj%sa,matZj%isa,matZj%jsa)
     enddo

     allocate(ddiag(Ntotal))
     ddiag(1:Ntotal)=DZERO
     ! extract diagonal locally
     do i=1,Nlocal
        ddiag(i+startj-1)=DONE
        do k=matZj%isa(i),matZj%isa(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix   
           if ((matZj%jsa(k)==startj+i-1).and.(abs(matZj%sa(k))/=DZERO))  ddiag(i+startj-1)=abs(matZj%sa(k))
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
        !  allocate(matZj%isa(nlocal+1))
        matZj%isa(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
        nnz=matZj%isa(nlocal+1)-1
        deallocate(matZj%jsa)
        allocate(matZj%jsa(nnz))
        matZj%jsa(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

        deallocate(matZj%sa)
        allocate(matZj%sa(1)) ! dummy
        allocate(tempsaz(1)) ! dummy

        do j=2,dmax+1 ! successives summation of csr matrices
           allocate(tempisaz(nlocal+1)) 
           allocate(tempjsaz(nnz))
           tempisaz(1:nlocal+1)=matZj%isa(1:nlocal+1)
           tempjsaz(1:nnz)=matZj%jsa(1:nnz)    
           opt=1 ! get isaz
           call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)     
           deallocate(matZj%jsa)
           nnz=matZj%isa(nlocal+1)-1 
           allocate(matZj%jsa(nnz))
           opt=2 ! get jsaz
           call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)  
           deallocate(tempisaz)
           deallocate(tempjsaz)
        end do
        deallocate(tempsaz) ! deallocate dummy
        deallocate(matZj%sa)
        allocate(matZj%sa(nnz))
        matZj%sa(1:nnz)= ZZERO ! initialize
        matZj%n=Nlocal
        matZj%nnz=nnz 

     end if
  end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then

     ! copy matrix single precision lbcsr matrices A 
     allocate(smatAjb(nb_procs3,nb_procs3,dmax+1))
     do j=1,dmax+1
        call zlbcsr2c(Nsize,mapA(1,1,j),matAjb(1,1,j),smatAjb(1,1,j),fpm(49),nb_procs3)
     end do
     ! get Az matrix local block csr with single precision (first get new single prec. matrix by rows)
     allocate(matCj%isa(Nlocal+1))
     allocate(matCj%jsa(nnz))
     allocate(matCj%sa(nnz))
     matCj%sa(1:nnz)=CZERO ! initialization (important for cinc_add below)
     matCj%isa(1:nlocal+1)= matZj%isa(1:nlocal+1)
     matCj%jsa(1:nnz)= matZj%jsa(1:nnz)
     matCj%n=Nlocal
     matCj%nnz=nnz

     allocate(mapC(nb_procs3,nb_procs3))
     allocate(matCjb(nb_procs3,nb_procs3,nfact))


     if (fpm(10)==1) then ! store factorization
        id=0
        do i=rank2+1,fpm(8),nb_procs2   !! multiple fact possible
           id=id+1
           call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,id),fpm(49),nb_procs3)
        end do
     else ! only 1 factorization 
        call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,1),fpm(49),nb_procs3)
     end if

  else ! double precision

     matZj%sa(1:nnz)=ZZERO ! initialization (important for zinc_add below)

     ! get Az matrix local block csr 
     allocate(mapZ(nb_procs3,nb_procs3))
     allocate(matZjb(nb_procs3,nb_procs3,nfact))


     if (fpm(10)==1) then ! store factorization
        id=0
        do i=rank2+1,fpm(8),nb_procs2   !! multiple fact possible
           id=id+1
           call zget_local_block_csr(Nsize,mapZ,matZj,matZjb(1,1,id),fpm(49),nb_procs3)
        end do
     else ! only 1 factorization 
        call zget_local_block_csr(Nsize,mapZ,matZj,matZjb(1,1,1),fpm(49),nb_procs3)
     end if

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  allocate(work(Nlocal,M02))
  allocate(zwork(Nlocal,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! iterative method  set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  allocate(nres(M0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(Nlocal,M0))
     allocate(cwork(Nlocal,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(Nlocal,M0))
  endif


  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if ((fpm(1)/=0).and.(rank3==0)) com=.true.
  !#ifdef MPI
  !     if (com.and.(rank/=0)) com=.false. ! comment only in rank 0 
  !#endif
  if (fpm(1)<0) then
     fout=abs(fpm(1))+200!fpm(60) !!file id name
  else
     fout=6 !screen
  endif




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

        if (fpm(42)==1) then !single precision fact- csaz

           do j=1,dmax+1
              call cinc_addlbcsr(Nsize,cmplx(Ze)**(j-1),mapA(1,1,j),smatAjb(1,1,j),mapC,matCjb(1,1,id),fpm(49),nb_procs3)
           enddo

        else ! double precision fact- saz

           do j=1,dmax+1
              call zinc_addlbcsr(Nsize,Ze**(j-1),mapA(1,1,j),matAjb(1,1,j),mapZ,matZjb(1,1,id),fpm(49),nb_procs3)    
           enddo
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)

!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!!!!!!!!!!       

        if (fpm(42)==1) then !single precision solve

           caux(1:Nlocal,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(Nlocal,zwork(1,i),1)**2
           enddo
           call MPI_ALLREDUCE(MPI_IN_PLACE,norm,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 
           do i=1,fpm(23)
              norm(i)=sqrt(norm(i))
              call ZSCAL(Nlocal,ZONE/norm(i),zwork(1,i),1)
           enddo

           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)

           call pcbicgstab(fpm(44),uplo,'N',Nsize,mapC,fpm(23),matCjb(1,1,id),cwork,caux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call CLAG2Z(Nlocal, fpm(23), caux, Nlocal, zwork, Nlocal, infoloc1)

           do i=1,fpm(23)
              call ZSCAL(Nlocal,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! double precision solve

           zaux(1:Nlocal,1:fpm(23))=ZZERO !! initial guess
           call pzbicgstab(fpm(44),uplo,'N',Nsize,mapZ,fpm(23),matZjb(1,1,id),zwork,zaux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call ZCOPY( Nlocal*fpm(23), zaux, 1, zwork, 1 )

        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

        fpm(60)=fpm(60)+linloops
        if (com) then
           if (nb_procs2>1) then
              if (.not.((rank2>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank2'
                 write(fout,'(I4)',advance='no') rank2
              endif
           endif
           write(fout,'(A)',advance='no') '  #it  '
           write(fout,'(I4)',advance='no') linloops
           write(fout,'(A)',advance='no') '; res min='
           write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
           write(fout,'(A)',advance='no') '; res max='
           write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
           write(fout,*)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)

!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!!!!!!!!!!       

        if (fpm(42)==1) then !single precision solve

           caux(1:Nlocal,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(Nlocal,zwork(1,i),1)**2
           enddo
           call MPI_ALLREDUCE(MPI_IN_PLACE,norm,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 
           do i=1,fpm(23)
              norm(i)=sqrt(norm(i))
              call ZSCAL(Nlocal,ZONE/norm(i),zwork(1,i),1)
           enddo

           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)

           call pcbicgstab(fpm(44),uplo2,'C',Nsize,mapC,fpm(23),matCjb(1,1,id),cwork,caux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call CLAG2Z(Nlocal, fpm(23), caux, Nlocal, zwork, Nlocal, infoloc1)

           do i=1,fpm(23)
              call ZSCAL(Nlocal,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! double precision solve

           zaux(1:Nlocal,1:fpm(23))=ZZERO !! initial guess
           call pzbicgstab(fpm(44),uplo2,'C',Nsize,mapZ,fpm(23),matZjb(1,1,id),zwork,zaux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call ZCOPY( Nlocal*fpm(23), zaux, 1, zwork, 1 )

        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

        fpm(60)=fpm(60)+linloops
        if (com) then
           if (nb_procs2>1) then
              if (.not.((rank2>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank2'
                 write(fout,'(I4)',advance='no') rank2
              endif
           endif
           write(fout,'(A)',advance='no') '  #it* '
           write(fout,'(I4)',advance='no') linloops
           write(fout,'(A)',advance='no') '; res min='
           write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
           write(fout,'(A)',advance='no') '; res max='
           write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
           write(fout,*)
        endif



     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
        call zlbcsrmm(UPLO2,'N',Nsize,mapA(1,1,fpm(57)),fpm(25),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(31) !! perform multiplication A(fpm(57)^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapA(1,1,fpm(57)),fpm(35),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)


     end select
  end do


  !loop=fpm(60)  ! new value for loop


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


  deallocate(matZj%isa)
  deallocate(matZj%jsa)
  deallocate(matZj%sa)


  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(matCj%isa)
     deallocate(matCj%jsa)
     deallocate(matCj%sa)
     deallocate(mapC)
     deallocate(matCjb)
     deallocate(smatAjb)
     deallocate(norm)
  else
     deallocate(zaux)
     deallocate(mapZ)
     deallocate(matZjb)
  end if

end subroutine pzifeast_hcsrpevx









subroutine pdifeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !    Remark: simple Wrapper to  pzifeast_gcsrpevx
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

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=222347 ! code name
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

  call pzifeast_gcsrpevx(dmax,N,zsa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(zsa)

end subroutine pdifeast_gcsrpevx






subroutine pzifeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  integer,dimension(:),allocatable :: tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
!!!! matrix format
  !  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa
  integer :: opt,nnza,nnz

!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout
  double precision :: DZNRM2

!!!!!!!!!!!!!!!!!!! for mixed precision
  !  complex,dimension(:,:),allocatable :: csaz
  !type ccsr_value
  !   complex,dimension(:),allocatable :: sa 
  !end type ccsr_value
  !type(ccsr_value), dimension(:),allocatable :: cmatAj 
  !complex,dimension(:),allocatable :: cfsa !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:,:),allocatable :: mapA
  integer,dimension(:,:),allocatable :: mapC,mapZ
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matZj
  type(zcsr),dimension(:),allocatable :: matAj
  type(zcsr), dimension(:,:,:),allocatable :: matAjb,matZjb
  !! Az matrix single
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr) :: matCj
  type(ccsr), dimension(:,:,:),allocatable :: matCjb
  type(ccsr), dimension(:,:,:),allocatable :: smatAjb




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=242347   ! code name
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
     CALL XERBLA( 'PZIFEAST_GCSRPEVX', -INFO+100 )
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
!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     allocate(matAj(dmax+1))
     do j=1,dmax+1
        call zcsr_distribute_row(Ntotal,sa(1,j),isa(1,j),jsa(1,j),matAj(j),startj,endj,Nlocal,Nsize,fpm(49)) !! A
     end do
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

!!! FORM the row distributed full matrix 
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
  allocate(matZj%isa(nlocal+1))
  matZj%isa(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=matZj%isa(nlocal+1)-1 
  allocate(matZj%jsa(nnz))
  matZj%jsa(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(matZj%sa(1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=matZj%isa(1:nlocal+1)
     tempjsaz(1:nnz)=matZj%jsa(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)     
     deallocate(matZj%jsa)
     nnz=matZj%isa(nlocal+1)-1 
     allocate(matZj%jsa(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)  
     deallocate(tempisaz)
     deallocate(tempjsaz)
  end do

  deallocate(tempsaz) ! deallocate dummy
  deallocate(matZj%sa)
  allocate(matZj%sa(nnz))
  matZj%sa(1:nnz)= ZZERO ! initialize
  matZj%n=Nlocal
  matZj%nnz=nnz    
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     Ze=Emid

     do j=1,dmax+1
        call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa(1),matAj(j)%isa(1),matAj(j)%jsa(1),matZj%sa,matZj%isa,matZj%jsa)
     enddo

     allocate(ddiag(Ntotal))
     ddiag(1:Ntotal)=DZERO
     ! extract diagonal locally
     do i=1,Nlocal
        ddiag(i+startj-1)=DONE
        do k=matZj%isa(i),matZj%isa(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix   
           if ((matZj%jsa(k)==startj+i-1).and.(abs(matZj%sa(k))/=DZERO))  ddiag(i+startj-1)=abs(matZj%sa(k))
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
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then

     ! copy matrix single precision lbcsr matrices A 
     allocate(smatAjb(nb_procs3,nb_procs3,dmax+1))
     do j=1,dmax+1
        call zlbcsr2c(Nsize,mapA(1,1,j),matAjb(1,1,j),smatAjb(1,1,j),fpm(49),nb_procs3)
     end do
     ! get Az matrix local block csr with single precision (first get new single prec. matrix by rows)
     allocate(matCj%isa(Nlocal+1))
     allocate(matCj%jsa(nnz))
     allocate(matCj%sa(nnz))
     matCj%sa(1:nnz)=CZERO ! initialization (important for cinc_add below)
     matCj%isa(1:nlocal+1)= matZj%isa(1:nlocal+1)
     matCj%jsa(1:nnz)= matZj%jsa(1:nnz)
     matCj%n=Nlocal
     matCj%nnz=nnz

     allocate(mapC(nb_procs3,nb_procs3))
     allocate(matCjb(nb_procs3,nb_procs3,nfact))


     if (fpm(10)==1) then ! store factorization
        id=0
        do i=rank2+1,fpm(8),nb_procs2   !! multiple fact possible
           id=id+1
           call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,id),fpm(49),nb_procs3)
        end do
     else ! only 1 factorization 
        call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,1),fpm(49),nb_procs3)
     end if

  else ! double precision

     matZj%sa(1:nnz)=ZZERO ! initialization (important for zinc_add below)

     ! get Az matrix local block csr 
     allocate(mapZ(nb_procs3,nb_procs3))
     allocate(matZjb(nb_procs3,nb_procs3,nfact))


     if (fpm(10)==1) then ! store factorization
        id=0
        do i=rank2+1,fpm(8),nb_procs2   !! multiple fact possible
           id=id+1
           call zget_local_block_csr(Nsize,mapZ,matZj,matZjb(1,1,id),fpm(49),nb_procs3)
        end do
     else ! only 1 factorization 
        call zget_local_block_csr(Nsize,mapZ,matZj,matZjb(1,1,1),fpm(49),nb_procs3)
     end if

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  allocate(work(Nlocal,M02))
  allocate(zwork(Nlocal,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! iterative method  set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  allocate(nres(M0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(Nlocal,M0))
     allocate(cwork(Nlocal,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(Nlocal,M0))
  endif


  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if ((fpm(1)/=0).and.(rank3==0)) com=.true.
  !#ifdef MPI
  !     if (com.and.(rank/=0)) com=.false. ! comment only in rank 0 
  !#endif
  if (fpm(1)<0) then
     fout=abs(fpm(1))+200!fpm(60) !!file id name
  else
     fout=6 !screen
  endif


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



        if (fpm(42)==1) then !single precision fact- csaz

           do j=1,dmax+1
              call cinc_addlbcsr(Nsize,cmplx(Ze)**(j-1),mapA(1,1,j),smatAjb(1,1,j),mapC,matCjb(1,1,id),fpm(49),nb_procs3)
           enddo

        else ! double precision fact- saz

           do j=1,dmax+1
              call zinc_addlbcsr(Nsize,Ze**(j-1),mapA(1,1,j),matAjb(1,1,j),mapZ,matZjb(1,1,id),fpm(49),nb_procs3)
           enddo
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)

!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!!!!!!!!!!       

        if (fpm(42)==1) then !single precision solve

           caux(1:Nlocal,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(Nlocal,zwork(1,i),1)**2
           enddo
           call MPI_ALLREDUCE(MPI_IN_PLACE,norm,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 
           do i=1,fpm(23)
              norm(i)=sqrt(norm(i))
              call ZSCAL(Nlocal,ZONE/norm(i),zwork(1,i),1)
           enddo

           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)

           call pcbicgstab(fpm(44),uplo2,'N',Nsize,mapC,fpm(23),matCjb(1,1,id),cwork,caux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call CLAG2Z(Nlocal, fpm(23), caux, Nlocal, zwork, Nlocal, infoloc1)

           do i=1,fpm(23)
              call ZSCAL(Nlocal,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! double precision solve

           zaux(1:Nlocal,1:fpm(23))=ZZERO !! initial guess
           call pzbicgstab(fpm(44),uplo2,'N',Nsize,mapZ,fpm(23),matZjb(1,1,id),zwork,zaux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call ZCOPY( Nlocal*fpm(23), zaux, 1, zwork, 1 )

        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

        fpm(60)=fpm(60)+linloops
        if (com) then
           if (nb_procs2>1) then
              if (.not.((rank2>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank2'
                 write(fout,'(I4)',advance='no') rank2
              endif
           endif
           write(fout,'(A)',advance='no') '  #it  '
           write(fout,'(I4)',advance='no') linloops
           write(fout,'(A)',advance='no') '; res min='
           write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
           write(fout,'(A)',advance='no') '; res max='
           write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
           write(fout,*)
        endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)

!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!!!!!!!!!!       

        if (fpm(42)==1) then !single precision solve

           caux(1:Nlocal,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(Nlocal,zwork(1,i),1)**2
           enddo
           call MPI_ALLREDUCE(MPI_IN_PLACE,norm,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 
           do i=1,fpm(23)
              norm(i)=sqrt(norm(i))
              call ZSCAL(Nlocal,ZONE/norm(i),zwork(1,i),1)
           enddo

           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)

           call pcbicgstab(fpm(44),uplo2,'C',Nsize,mapC,fpm(23),matCjb(1,1,id),cwork,caux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call CLAG2Z(Nlocal, fpm(23), caux, Nlocal, zwork, Nlocal, infoloc1)

           do i=1,fpm(23)
              call ZSCAL(Nlocal,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! double precision solve

           zaux(1:Nlocal,1:fpm(23))=ZZERO !! initial guess
           call pzbicgstab(fpm(44),uplo2,'C',Nsize,mapZ,fpm(23),matZjb(1,1,id),zwork,zaux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call ZCOPY( Nlocal*fpm(23), zaux, 1, zwork, 1 )

        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

        fpm(60)=fpm(60)+linloops
        if (com) then
           if (nb_procs2>1) then
              if (.not.((rank2>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank2'
                 write(fout,'(I4)',advance='no') rank2
              endif
           endif
           write(fout,'(A)',advance='no') '  #it* '
           write(fout,'(I4)',advance='no') linloops
           write(fout,'(A)',advance='no') '; res min='
           write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
           write(fout,'(A)',advance='no') '; res max='
           write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
           write(fout,*)
        endif


     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
        call zlbcsrmm(UPLO2,'N',Nsize,mapA(1,1,fpm(57)),fpm(25),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(31) !! perform multiplication A(fpm(57)^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapA(1,1,fpm(57)),fpm(35),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)


     end select
  end do


  !loop=fpm(60)  ! new value for loop


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

  deallocate(matZj%isa)
  deallocate(matZj%jsa)
  deallocate(matZj%sa)

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(matCj%isa)
     deallocate(matCj%jsa)
     deallocate(matCj%sa)
     deallocate(mapC)
     deallocate(matCjb)
     deallocate(smatAjb)
     deallocate(norm)
  else
     deallocate(zaux)
     deallocate(mapZ)
     deallocate(matZjb)
  end if


end subroutine pzifeast_gcsrpevx











subroutine pzifeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: tempsaz
  double precision, dimension(:),allocatable ::ddiag
  integer,dimension(:),allocatable :: tempisaz,tempjsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
  ! complex(kind=(kind(1.0d0))),dimension(:),allocatable :: usa
  ! integer,dimension(:), allocatable :: uisa,ujsa
  integer :: opt,nnza,nnz

!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout
  double precision :: DZNRM2

!!!!!!!!!!!!!!!!!!! for mixed precision
  !complex,dimension(:,:),allocatable :: csaz
  !type ccsr_value
  !   complex,dimension(:),allocatable :: sa 
  !end type ccsr_value
  !complex,dimension(:,:),allocatable :: susa !single precision copy
  !type(ccsr_value), dimension(:),allocatable :: cmatAj 
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:,:),allocatable :: mapA
  integer,dimension(:,:),allocatable :: mapZ,mapC
  type zcsr
     integer ::n,m,nnz
     complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matZj
  type(zcsr), dimension(:),allocatable :: matAj
  type(zcsr), dimension(:,:,:),allocatable :: matAjb,matZjb
  !! Az matrix single
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr) :: matCj
  type(ccsr), dimension(:,:,:),allocatable :: matCjb
  type(ccsr), dimension(:,:,:),allocatable :: smatAjb

!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=242147  ! code name
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
     CALL XERBLA( 'PZIFEAST_SCSRPEVX', -INFO+100 )
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
!!! DISTRIBUTE the sparse matrix by rows among all L3 procs  
     allocate(Nsize(nb_procs3))
     allocate(matAj(dmax+1))
     do j=1,dmax+1
        nnzA=isa(n+1,j)-1
        call zcsr_distribute_row(Ntotal,sa(1,j),isa(1,j),jsa(1,j),matAj(j),startj,endj,Nlocal,Nsize,fpm(49)) !! Aj
     enddo

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
!!! FORM the row distributed matrix   
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
  allocate(matZj%isa(nlocal+1))
  matZj%isa(1:nlocal+1)=matAj(1)%isa(1:nlocal+1) ! initialize to A[1]
  nnz=matZj%isa(nlocal+1)-1 
  allocate(matZj%jsa(nnz))
  matZj%jsa(1:nnz)=matAj(1)%jsa(1:nnz) ! initialize to A[1]

  allocate(matZj%sa(1)) ! dummy
  allocate(tempsaz(1)) ! dummy

  do j=2,dmax+1 ! successives summation of csr matrices
     allocate(tempisaz(nlocal+1)) 
     allocate(tempjsaz(nnz))
     tempisaz(1:nlocal+1)=matZj%isa(1:nlocal+1)
     tempjsaz(1:nnz)=matZj%jsa(1:nnz)    
     opt=1 ! get isaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)     
     deallocate(matZj%jsa)
     nnz=matZj%isa(nlocal+1)-1 
     allocate(matZj%jsa(nnz))
     opt=2 ! get jsaz
     call zaddcsr(Nlocal,Ntotal,opt,ZONE,tempsaz,tempisaz,tempjsaz,ZONE,matAj(j)%sa,matAj(j)%isa,matAj(j)%jsa,matZj%sa,matZj%isa,matZj%jsa)  
     deallocate(tempisaz)
     deallocate(tempjsaz)
  end do
  deallocate(tempsaz) ! deallocate dummy
  deallocate(matZj%sa)
  allocate(matZj%sa(nnz))
  matZj%sa(1:nnz)= ZZERO ! initialize
  matZj%n=Nlocal
  matZj%nnz=nnz
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     Ze=Emid
     do j=1,dmax+1
        call zinc_addcsr(Nlocal,Ntotal,Ze**(j-1),matAj(j)%sa(1),matAj(j)%isa(1),matAj(j)%jsa(1),matZj%sa,matZj%isa,matZj%jsa)
     enddo

     allocate(ddiag(Ntotal))
     ddiag(1:Ntotal)=DZERO
     ! extract diagonal locally
     do i=1,Nlocal
        ddiag(i+startj-1)=DONE
        do k=matZj%isa(i),matZj%isa(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix   
           if ((matZj%jsa(k)==startj+i-1).and.(abs(matZj%sa(k))/=DZERO))  ddiag(i+startj-1)=abs(matZj%sa(k))
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

  allocate(mapA(nb_procs3,nb_procs3,dmax+1))
  allocate(matAjb(nb_procs3,nb_procs3,dmax+1))
  do j=1,dmax+1
     call zget_local_block_csr(Nsize,mapA(1,1,j),matAj(j),matAjb(1,1,j),fpm(49),nb_procs3)
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then

     ! copy matrix single precision lbcsr matrices A 
     allocate(smatAjb(nb_procs3,nb_procs3,dmax+1))
     do j=1,dmax+1
        call zlbcsr2c(Nsize,mapA(1,1,j),matAjb(1,1,j),smatAjb(1,1,j),fpm(49),nb_procs3)
     end do
     ! get Az matrix local block csr with single precision (first get new single prec. matrix by rows)
     allocate(matCj%isa(Nlocal+1))
     allocate(matCj%jsa(nnz))
     allocate(matCj%sa(nnz))
     matCj%sa(1:nnz)=CZERO ! initialization (important for cinc_add below)
     matCj%isa(1:nlocal+1)= matZj%isa(1:nlocal+1)
     matCj%jsa(1:nnz)= matZj%jsa(1:nnz)
     matCj%n=Nlocal
     matCj%nnz=nnz

     allocate(mapC(nb_procs3,nb_procs3))
     allocate(matCjb(nb_procs3,nb_procs3,nfact))


     if (fpm(10)==1) then ! store factorization
        id=0
        do i=rank2+1,fpm(8),nb_procs2   !! multiple fact possible
           id=id+1
           call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,id),fpm(49),nb_procs3)
        end do
     else ! only 1 factorization 
        call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,1),fpm(49),nb_procs3)
     end if

  else ! double precision

     matZj%sa(1:nnz)=ZZERO ! initialization (important for zinc_add below)

     ! get Az matrix local block csr 
     allocate(mapZ(nb_procs3,nb_procs3))
     allocate(matZjb(nb_procs3,nb_procs3,nfact))


     if (fpm(10)==1) then ! store factorization
        id=0
        do i=rank2+1,fpm(8),nb_procs2   !! multiple fact possible
           id=id+1
           call zget_local_block_csr(Nsize,mapZ,matZj,matZjb(1,1,id),fpm(49),nb_procs3)
        end do
     else ! only 1 factorization 
        call zget_local_block_csr(Nsize,mapZ,matZj,matZjb(1,1,1),fpm(49),nb_procs3)
     end if

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(Aq(M0*dmax,M0*dmax))
  allocate(Bq(M0*dmax,M0*dmax))
  allocate(work(Nlocal,M0))
  allocate(zwork(Nlocal,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! iterative method  set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  allocate(nres(M0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(Nlocal,M0))
     allocate(cwork(Nlocal,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(Nlocal,M0))
  endif


  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if ((fpm(1)/=0).and.(rank3==0)) com=.true.
  !#ifdef MPI
  !     if (com.and.(rank/=0)) com=.false. ! comment only in rank 0 
  !#endif
  if (fpm(1)<0) then
     fout=abs(fpm(1))+200!fpm(60) !!file id name
  else
     fout=6 !screen
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_srcipevx(ijob,dmax,Nlocal,Ze,work,zwork,Aq,Bq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize P(Ze)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        id=fpm(33) !! id of factorization (for fpm(10) flag)

        if (fpm(42)==1) then !single precision fact- csaz

           do j=1,dmax+1
              call cinc_addlbcsr(Nsize,cmplx(Ze)**(j-1),mapA(1,1,j),smatAjb(1,1,j),mapC,matCjb(1,1,id),fpm(49),nb_procs3)
           enddo

        else ! double precision fact- saz 
           do j=1,dmax+1
              call zinc_addlbcsr(Nsize,Ze**(j-1),mapA(1,1,j),matAjb(1,1,j),mapZ,matZjb(1,1,id),fpm(49),nb_procs3)
           enddo
        end if
        !print *,'id-nnz',id,nnz
        !        print *,sum(abs(matZjb(1,1,id)%sa(1:nnz)))
        !return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)

!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!!!!!!!!!!       

        if (fpm(42)==1) then !single precision solve

           caux(1:Nlocal,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(Nlocal,zwork(1,i),1)**2
           enddo
           call MPI_ALLREDUCE(MPI_IN_PLACE,norm,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code) 
           do i=1,fpm(23)
              norm(i)=sqrt(norm(i))
              call ZSCAL(Nlocal,ZONE/norm(i),zwork(1,i),1)
           enddo

           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)

           call pcbicgstab(fpm(44),uplo,'N',Nsize,mapC,fpm(23),matCjb(1,1,id),cwork,caux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call CLAG2Z(Nlocal, fpm(23), caux, Nlocal, zwork, Nlocal, infoloc1)

           do i=1,fpm(23)
              call ZSCAL(Nlocal,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! double precision solve

           zaux(1:Nlocal,1:fpm(23))=ZZERO !! initial guess
           call pzbicgstab(fpm(44),uplo,'N',Nsize,mapZ,fpm(23),matZjb(1,1,id),zwork,zaux,nres,linloops,lintargeterror,comb,fpm(49),nb_procs3,infoloc2)
           call ZCOPY( Nlocal*fpm(23), zaux, 1, zwork, 1 )

        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

        fpm(60)=fpm(60)+linloops
        if (com) then
           if (nb_procs2>1) then
              if (.not.((rank2>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank2'
                 write(fout,'(I4)',advance='no') rank2
              endif
           endif
           write(fout,'(A)',advance='no') '  #it  '
           write(fout,'(I4)',advance='no') linloops
           write(fout,'(A)',advance='no') '; res min='
           write(fout,'(ES25.16)',advance='no') minval(nres(1:fpm(23)))
           write(fout,'(A)',advance='no') '; res max='
           write(fout,'(ES25.16)',advance='no') maxval(nres(1:fpm(23)))
           write(fout,*)
        endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call zlbcsrmm(UPLO,'N',Nsize,mapA(1,1,fpm(57)),fpm(25),ZONE,matAjb(1,1,fpm(57)),Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     end select
  end do


  !loop=fpm(60)  ! new value for loop



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

  deallocate(matZj%isa)
  deallocate(matZj%jsa)
  deallocate(matZj%sa)

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(matCj%isa)
     deallocate(matCj%jsa)
     deallocate(matCj%sa)
     deallocate(mapC)
     deallocate(matCjb)
     deallocate(smatAjb)
     deallocate(norm)
  else
     deallocate(zaux)
     deallocate(mapZ)
     deallocate(matZjb)
  endif


end subroutine pzifeast_scsrpevx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pdifeast_scsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pdifeast_scsrpevx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=222146
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pdifeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdifeast_scsrpev




subroutine pzifeast_hcsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pzifeast_hcsrpevx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=242246
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzifeast_hcsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzifeast_hcsrpev



subroutine pdifeast_gcsrpev(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pdifeast_gcsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=222346
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pdifeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdifeast_gcsrpev




subroutine pzifeast_gcsrpev(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pzifeast_gcsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=242346
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzifeast_gcsrpevx(dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzifeast_gcsrpev




subroutine pzifeast_scsrpev(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pzifeast_scsrpevx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=242146
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzifeast_scsrpevx(UPLO,dmax,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pzifeast_scsrpev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



