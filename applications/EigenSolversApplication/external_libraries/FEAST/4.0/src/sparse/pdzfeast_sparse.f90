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

!P{D,Z}FEAST_{SCSR,HCSR,GCSR}{EV,GV}


!111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE) - using MKL-CLUSTER-PARDISO
! pdfeast_scsrgvx
! pzfeast_hcsrgvx
! pdfeast_gcsrgvx
! pzfeast_gcsrgvx
! pzfeast_scsrgvx

!222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
! pdfeast_scsrgv
! pzfeast_hcsrgv
! pdfeast_gcsrgv
! pzfeast_gcsrgv
! pzfeast_scsrgv

!333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
! pdfeast_scsrevx
! pzfeast_hcsrevx
! pdfeast_gcsrevx
! pzfeast_gcsrevx
! pzfeast_scsrevx

!44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
! pdfeast_scsrev
! pzfeast_hcsrev
! pdfeast_gcsrev
! pzfeast_gcsrev
! pzfeast_scsrev


!5555555555555!!! Routines for Estimating the Search interval Emin,Emax for the real symmetric or Hermitian problem
!pdfeast_scsrev_search
!pzfeast_hcsrev_search
!pdfeast_scsrgv_search
!pzfeast_hcsrgv_search


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pdfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        REAL DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) REAL DOUBLE PRECISION(N,M0) : 
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
  ! James Kestyn 2016-2018
  ! Eric Polizzi 2019 
  !=====================================================================
  implicit none
  !  include 'f90_noruntime_interface.fi'
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  double precision,dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k,j
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  double precision, dimension(:,:),allocatable ::work,Aq,Sq,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
  double precision,dimension(:),allocatable :: usa,usb
  integer,dimension(:), allocatable :: uisa,ujsa,uisb,ujsb
  integer :: opt,nnza,nnzb,nnz
!!!!!for cluster pardiso
  integer(8),dimension(:,:),allocatable :: pt 
  integer,dimension(:,:),allocatable :: iparm
  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl

!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  real,dimension(:),allocatable :: susa,susb !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork

!!!!!!!!!!!!!!!!! case 50
  !  integer :: info_lap,lwork_loc
  !  double precision,dimension(:),allocatable :: work_loc
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: matAj,matBj
  type(dcsr), dimension(:,:),allocatable :: matAjb,matBjb


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif
  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
     call pdifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  


  !  LWORK_LOC=3*M0-1
  !  allocate(WORK_LOC,LWORK_LOC))

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=221145  ! code name
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
     CALL XERBLA( 'PDFEAST_SCSRGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Find which distribution scheme is used for the input matrix/rhs: global or local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call pfeast_distribution_type(N,isa,jsa,fpm(49),distype)
  !fpm(50)=distype


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!! FORMAT CONVERSION TO CSR-UPPER     
     !! for A
     allocate(uisa(n+1))
     nnzA=isa(n+1)-1
     if ((UPLO=='F').or.(UPLO=='f')) nnzA=nnzA/2+n ! slightly overestimated
     allocate(ujsa(nnzA))
     allocate(usa(nnzA))
     call dcsr_convert_upper(n,UPLO,sa,isa,jsa,usa,uisa,ujsa)
     !! for B
     allocate(uisb(n+1))
     nnzB=isb(n+1)-1
     if ((UPLO=='F').or.(UPLO=='f')) nnzB=nnzB/2+n ! slightly overestimated
     allocate(ujsb(nnzB))
     allocate(usb(nnzB))
     call dcsr_convert_upper(n,UPLO,sb,isb,jsb,usb,uisb,ujsb)

!!! DISTRIBUTE the (upper) sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call dcsr_distribute_row(Ntotal,usa,uisa,ujsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call dcsr_distribute_row(Ntotal,usb,uisb,ujsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B


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

!!! FORM the row distributed upper matrix for Pardiso   
     if ((UPLO=='F').or.(UPLO=='f')) then !matrix need to be reduced to upper part
        !! for A
        nnzA=isa(n+1)-1 ! overestimation
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%isa(1)=1
        do i=1,N
           j=0
           do k=isa(i),isa(i+1)-1
              if (jsa(k)>=i+startj-1) then
                 j=j+1
                 matAj%jsa(matAj%isa(i)+j-1)=jsa(k)
                 matAj%sa(matAj%isa(i)+j-1)=sa(k)
              end if
           enddo
           matAj%isa(i+1)=matAj%isa(i)+j 
        enddo
        matAj%n=N
        matAj%nnz=matAj%isa(n+1)-1
        !! for B
        nnzB=isb(n+1)-1 ! overestimation
        allocate(matBj%isa(N+1))
        allocate(matBj%jsa(nnzB))
        allocate(matBj%sa(nnzB))
        matBj%isa(1)=1
        do i=1,N
           j=0
           do k=isb(i),isb(i+1)-1
              if (jsb(k)>=i+startj-1) then
                 j=j+1
                 matBj%jsa(matBj%isa(i)+j-1)=jsb(k)
                 matBj%sa(matBj%isa(i)+j-1)=sb(k)
              end if
           enddo
           matBj%isa(i+1)=matBj%isa(i)+j 
        enddo
        matBj%n=N
        matBj%nnz=matBj%isa(n+1)-1

     else ! include the following cases
        !! ((UPLO=='U').or.(UPLO=='u')) !! upper-csr already -- simple copy
        !! ((UPLO=='L').or.(UPLO=='l')) !! temp-copy !! transpose difficult..could be done once the local blocks csr are defined below

        !! for A
        nnzA=isa(n+1)-1
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%sa(1:nnzA)=sa(1:nnzA)
        matAj%isa(1:n+1)= isa(1:n+1)
        matAj%jsa(1:nnzA)= jsa(1:nnzA)
        matAj%n=N
        matAj%nnz=nnzA
        !! for B
        nnzB=isb(n+1)-1
        allocate(matBj%isa(N+1))
        allocate(matBj%jsa(nnzB))
        allocate(matBj%sa(nnzB))
        matBj%sa(1:nnzB)=sb(1:nnzB)
        matBj%isa(1:n+1)= isb(1:n+1)
        matBj%jsa(1:nnzB)= jsb(1:nnzB)
        matBj%n=N
        matBj%nnz=nnzB    
     end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank2+1,fpm(2),nb_procs2
        nfact=nfact+1
     end do
  endif


!!!!!!! Set up for Az matrix (isaz and jsaz)

  allocate(isaz(Nlocal+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 !! get isaz
  call zdaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  nnz=isaz(Nlocal+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !! get jsaz
  call zdaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zdaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,1),isaz,jsaz)

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
     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo

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
  allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
  call dget_local_block_csr(Nsize,mapA,matAj,matAjb,fpm(49),nb_procs3)
  allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call dget_local_block_csr(Nsize,mapB,matBj,matBjb,fpm(49),nb_procs3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! HANDLE THE PARTICULAR CASE of user distributed matrices via Lower format (get it upper for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (distype==1) then
     if ((UPLO=='L').or.(UPLO=='l')) then !! perform the transpose
        !! create local block upper format

        call dlbcsr_transpose(Nsize,mapA,matAjb,L3_COMM_WORLD,nb_procs3)
        call dlbcsr_transpose(Nsize,mapB,matBjb,L3_COMM_WORLD,nb_procs3)

        !! convert into row distributed csr format
        deallocate(matAj%jsa)
        deallocate(matAj%sa)
        deallocate(matAj%isa)
        call dlbcsr_distribute_row(Nsize,mapA,matAjb,matAj,L3_COMM_WORLD,nb_procs3)
        deallocate(matBj%jsa)
        deallocate(matBj%sa)
        deallocate(matBj%isa)
        call dlbcsr_distribute_row(Nsize,mapB,matBjb,matBj,L3_COMM_WORLD,nb_procs3)

!!$        UPLO2='L'        
!!$        p=rank3+1
!!$
!!$        !! A
!!$deallocate(matAj%jsa)
!!$deallocate(matAj%sa)
!!$allocate(matAj%jsa(1:sum(mapA(p:nb_procs3,p))))
!!$allocate(matAj%sa(1:sum(mapA(p:nb_procs3,p))))
!!$call dget_distribute_row_from_block_csr(p,nb_procs3,Nsize,mapA,matAjb,matAj)  
!!$
!!$ !! B      
!!$deallocate(matBj%jsa)
!!$deallocate(matBj%sa)
!!$allocate(matBj%jsa(1:sum(mapB(p:nb_procs3,p))))
!!$allocate(matBj%sa(1:sum(mapB(p:nb_procs3,p))))
!!$call dget_distribute_row_from_block_csr(p,nb_procs3,Nsize,mapB,matBjb,matBj)  


!!! we need to reconstruct isaz and jsaz, as upper local rows
        allocate(saz(1,1)) ! dummy
        opt=1 !! get isaz
        call zdaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
        nnz=isaz(Nlocal+1)-1
        deallocate(jsaz)
        allocate(jsaz(nnz))
        opt=2 !! get jsaz
        call zdaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
        deallocate(saz) ! dummy
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision   
     allocate(susa(matAj%isa(Nlocal+1)-1))
     allocate(susb(matBj%isa(Nlocal+1)-1))
     call DLAG2S(matAj%isa(Nlocal+1)-1,1,matAj%sa,matAj%isa(Nlocal+1)-1,susa,matAj%isa(Nlocal+1)-1,infoloc1)
     call DLAG2S(matBj%isa(Nlocal+1)-1,1,matBj%sa,matBj%isa(Nlocal+1)-1,susb,matBj%isa(Nlocal+1)-1,infoloc1)     
     !susa(1:matAj%isa(Nlocal+1)-1)=real(matAj%sa(1:matAj%isa(Nlocal+1)-1))
     !susb(1:matBj%isa(Nlocal+1)-1)=real(matBj%sa(1:matBj%isa(Nlocal+1)-1)) 
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
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
     !! Symbolic factorization
!!$       opt=3  !! get csaz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call csaddcsr(Nlocal,Ntotal,opt,-CONE,susa,matAj%isa,matAj%jsa,cmplx(Ze),susb,matBj%isa,matBj%jsa,csaz(1,i),isaz,jsaz)
!!$          if (nb_procs3>1) then
!!$             call  cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
!!$             else
!!$                call  pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,infoloc2)
!!$                endif
!!$       enddo
!!$       
  else !double precision 

     allocate(zaux(Nlocal,M0))
     !! Symbolic factorization
!!$       opt=3 !! get saz 
!!$ do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$    Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$    call zdaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,i),isaz,jsaz)
!!$     if (nb_procs3>1) then
!!$        call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
!!$     else
!!$  call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,infoloc2)
!!$     end if
!!$       enddo
  end if

!!$ 
!!$  if (infoloc2/=0) then
!!$     info=-2
!!$     return
!!$  end if
!!$
!!$  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
     call dfeast_srcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,Xj,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization

        if (fpm(42)==1) then !single precision fact- csaz
           call csaddcsr(Nlocal,Ntotal,opt,-CONE,susa,matAj%isa,matAj%jsa,cmplx(Ze),susb,matBj%isa,matBj%jsa,csaz(1,id),isaz,jsaz) !! get csaz
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if

        else ! double precision fact- saz
           call zdaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,id),isaz,jsaz) !! get saz
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
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call dlbcsrmm(UPLO2,'N',Nsize,mapA,fpm(25),DONE,matAjb,Xj(1,fpm(24)),DZERO,work(1,fpm(24)),fpm(49),nb_procs3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call dlbcsrmm(UPLO2,'N',Nsize,mapB,fpm(25),DONE,matBjb,Xj(1,fpm(24)),DZERO,work(1,fpm(24)),fpm(49),nb_procs3)

!!$     case(50) ! compute the inner product Aq=X^work
!!$        print *,'inner'
!!$        ! Aq(1:M0,1:M0)=matmul(transpose(X(1:N,1:M0)),work(1:N,1:M0))
!!$        call DGEMM('T','N',fpm(23),fpm(25),N,1.0d0,Xj(1,1),N,work(1,fpm(24)),N,0.0d0,Aq(1,fpm(24)),M0) 
!!$
!!$     case(60) !! solve the reduced eigenvalue problem Aqxq=\lambda_qBqxq
!!$        print *,'eigen'
!!$        call DSYGV(1, 'V', 'L', fpm(23),Aq,M0,Sq,M0,E,work_loc,Lwork_loc,INFO_lap)

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
  deallocate(matAj%isa)
  deallocate(matAj%jsa)
  deallocate(matAj%sa)
  deallocate(matBj%isa)
  deallocate(matBj%jsa)
  deallocate(matBj%sa)
  deallocate(mapA)
  deallocate(mapB)
  deallocate(matAjb)
  deallocate(matBjb)


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

  if (distype==0) then
     deallocate(usa)
     deallocate(uisa)
     deallocate(ujsa)

     deallocate(usb)
     deallocate(uisb)
     deallocate(ujsb)
  end if

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(susa)
     deallocate(susb)
  else
     deallocate(zaux)
     deallocate(saz)
  endif



#endif
!!!! done with MKL CLUSTER PARDISO

end subroutine pdfeast_scsrgvx










subroutine pzfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) REAL DOUBLE PRECISION(N,M0) : 
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
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k,j
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,Aq,Sq,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  integer,dimension(:), allocatable :: fisa,fjsa,fisb,fjsb
  integer :: opt,nnza,nnzb,nnz
!!!!!for cluster pardiso
  integer(8),dimension(:,:),allocatable :: pt 
  integer,dimension(:,:),allocatable :: iparm
  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj,matBj
  type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
     call pzifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=241245  ! code name
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
     CALL XERBLA( 'PZFEAST_HCSRGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Find which distribution scheme is used for the input matrix/rhs: global or local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call pfeast_distribution_type(N,isa,jsa,fpm(49),distype)
  !fpm(50)=distype


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!! FORMAT CONVERSION TO FULL-CSR
     !! for A   
     allocate(fisa(n+1))
     nnza=isa(n+1)-1
     if ((UPLO/='F').and.(UPLO/='f')) nnza=2*nnza ! slightly overestimated
     allocate(fsa(nnza))
     allocate(fjsa(nnza))
     call   zhcsr_convert_full(n,UPLO,sa,isa,jsa,fsa,fisa,fjsa)
     !! for B   
     allocate(fisb(n+1))
     nnzb=isa(n+1)-1
     if ((UPLO/='F').and.(UPLO/='f')) nnzb=2*nnzb ! slightly overestimated
     allocate(fsb(nnzb))
     allocate(fjsb(nnzb))
     call   zhcsr_convert_full(n,UPLO,sb,isb,jsb,fsb,fisb,fjsb)

!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,fsa,fisa,fjsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call zcsr_distribute_row(Ntotal,fsb,fisb,fjsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B

!!! DISTRIBUTE the RHS by rows among all L3 procs
     allocate(Xj(Nlocal,M0))
     Xj(1:Nlocal,1:M0)=X(startj:endj,1:M0) !!<< blas?

  end if !! global distype=0

  !print *,'rank',rank3+1,matAj%sa(:)

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
     allocate(Xj(Nlocal,M0))
     if (fpm(5)==1) Xj(1:Nlocal,1:M0)=X(1:Nlocal,1:M0) !!<< blas?


!!! FORM the row distributed full matrix for cluster Pardiso   
     ! first simple copy
     ! the case F is then taking care of
     ! the cases L and U use this temp copy and will be modified to
     ! full format once the local blocks csr are defined below

     nnzA=isa(n+1)-1
     allocate(matAj%isa(N+1))
     allocate(matAj%jsa(nnzA))
     allocate(matAj%sa(nnzA))
     matAj%sa(1:nnzA)=sa(1:nnzA)
     matAj%isa(1:n+1)= isa(1:n+1)
     matAj%jsa(1:nnzA)= jsa(1:nnzA)
     matAj%n=N
     matAj%nnz=nnzA
     !! for B
     nnzB=isb(n+1)-1
     allocate(matBj%isa(N+1))
     allocate(matBj%jsa(nnzB))
     allocate(matBj%sa(nnzB))
     matBj%sa(1:nnzB)=sb(1:nnzB)
     matBj%isa(1:n+1)= isb(1:n+1)
     matBj%jsa(1:nnzB)= jsb(1:nnzB)
     matBj%n=N
     matBj%nnz=nnzB  

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank2+1,fpm(2),nb_procs2
        nfact=nfact+1
     end do
  endif


!!!!!!! Set up for Az matrix (isaz and jsaz)

  allocate(isaz(Nlocal+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 !! get isaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  nnz=isaz(Nlocal+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !! get jsaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,1),isaz,jsaz)

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
     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo

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

  UPLO2='F' ! full csr local blocks by default
  allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapA,matAj,matAjb,fpm(49),nb_procs3)
  allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapB,matBj,matBjb,fpm(49),nb_procs3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! HANDLE THE PARTICULAR CASE of user distributed matrices via Lower or upper format -- we need to transform them into full row csr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (distype==1) then
     if ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u') ) then 


        !! transform into full local block csr format        
        call zhlbcsr_uplo_to_csr(Nsize,mapA,matAjb,L3_COMM_WORLD,nb_procs3)
        call zhlbcsr_uplo_to_csr(Nsize,mapB,matBjb,L3_COMM_WORLD,nb_procs3)  
        !! convert into row distributed csr format
        deallocate(matAj%jsa)
        deallocate(matAj%sa)
        deallocate(matAj%isa)
        call zlbcsr_distribute_row(Nsize,mapA,matAjb,matAj,L3_COMM_WORLD,nb_procs3)
        deallocate(matBj%jsa)
        deallocate(matBj%sa)
        deallocate(matBj%isa)
        call zlbcsr_distribute_row(Nsize,mapB,matBjb,matBj,L3_COMM_WORLD,nb_procs3)


!!! we need to reconstruct isaz and jsaz, as full local rows
        allocate(saz(1,1)) ! dummy
        opt=1 !! get isaz
        call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
        nnz=isaz(Nlocal+1)-1
        deallocate(jsaz)
        allocate(jsaz(nnz))
        opt=2 !! get jsaz
        call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
        deallocate(saz) ! dummy
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision   
     allocate(cfsa(matAj%isa(Nlocal+1)-1))
     allocate(cfsb(matBj%isa(Nlocal+1)-1))
     call ZLAG2C(matAj%isa(Nlocal+1)-1,1,matAj%sa,matAj%isa(Nlocal+1)-1,cfsa,matAj%isa(Nlocal+1)-1,infoloc1)
     call ZLAG2C(matBj%isa(Nlocal+1)-1,1,matBj%sa,matBj%isa(Nlocal+1)-1,cfsb,matBj%isa(Nlocal+1)-1,infoloc1)     
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
  allocate(work(Nlocal,M0))
  allocate(zwork(Nlocal,M0))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  cluster pardiso initialization
!!!!!!!!  use same factorization for (normal+transpose solve)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1
  MNUM=1
  MTYPE=13!3       !Rq: it should be 3 since complex and structurally symmetric but cluster-pardiso fails
  ! for system2 with nL3>=2 (segmentation fault)
  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))
  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  enddo
  IPARM(11,:)=0 ! disable scaling !done by feast fpm(41)
  ! IPARM(2)=10 ! distributed nested (found a bug for 4x4 matrix, 2mpi)
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

  IPARM(40,:)=2 ! matrix storage (distributed)
  IPARM(41,:)=startj ! start row for local storage of saz
  IPARM(42,:)=endj ! end row for local storage of saz
!!!!!!!!!!


  if (fpm(42)==1) then ! mixed precision (single precision solver)
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<   
     allocate(caux(Nlocal,M0))
     allocate(cwork(Nlocal,M0))
!!$  !! Symbolic factorization
!!$       opt=3  !! get csaz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call caddcsr(Nlocal,Ntotal,opt,-CONE,cfsa,matAj%isa,matAj%jsa,cmplx(Ze),cfsb,matBj%isa,matBj%jsa,csaz(1,i),isaz,jsaz)
!!$          if (nb_procs3>1) then
!!$          call  cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
!!$else
!!$ call  pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,infoloc2)
!!$   endif
!!$       enddo
!!$       
  else !double precision 

     allocate(zaux(Nlocal,M0))
     !! Symbolic factorization
!!$       opt=3 !! get saz 
!!$ do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$    Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$    call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,i),isaz,jsaz)
!!$
!!$    if (nb_procs3>1) then
!!$call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
!!$else
!!$call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,infoloc2)
!!$end if
!!$enddo

  end if

!!$ 
!!$  if (infoloc2/=0) then
!!$     info=-2
!!$     return
!!$  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_hrcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,Xj,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization
        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(Nlocal,Ntotal,opt,-CONE,cfsa,matAj%isa,matAj%jsa,cmplx(Ze),cfsb,matBj%isa,matBj%jsa,csaz(1,id),isaz,jsaz) !! get csaz

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if


        else ! double precision fact- saz
           call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,id),isaz,jsaz) !! get saz    

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
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=0 ! normal solve


        if (fpm(42)==1) then !single precision solve
           call ZLAG2C(Nlocal, fpm(23), zwork, Nlocal, cwork, Nlocal, infoloc1)
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if

           call CLAG2Z(Nlocal, fpm(23), cwork, Nlocal, zwork, Nlocal, infoloc1)

        else ! double precision solve
           ! print *,'bef',zwork(1:N,1)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system (ZeB-A)^H x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call zlbcsrmm(UPLO2,'N',Nsize,mapA,fpm(25),ZONE,matAjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call zlbcsrmm(UPLO2,'N',Nsize,mapB,fpm(25),ZONE,matBjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)



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
  deallocate(matAj%isa)
  deallocate(matAj%jsa)
  deallocate(matAj%sa)
  deallocate(matBj%isa)
  deallocate(matBj%jsa)
  deallocate(matBj%sa)
  deallocate(mapA)
  deallocate(mapB)
  deallocate(matAjb)
  deallocate(matBjb)


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

  if (distype==0) then
     deallocate(fsa)
     deallocate(fisa)
     deallocate(fjsa)

     deallocate(fsb)
     deallocate(fisb)
     deallocate(fjsb)
  end if

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(cfsa)
     deallocate(cfsb)
  else
     deallocate(zaux)
     deallocate(saz)
  endif


#endif
!!!! done with MKL CLUSTER PARDISO

end subroutine pzfeast_hcsrgvx





subroutine pdfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        REAL DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
  integer :: N
  double precision,dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
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
  integer :: ijob,infoloc1,infoloc2,i,s,k,j,M02
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,Aq,Sq,zaux,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  integer :: opt,nnza,nnzb,nnz
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
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj,matBj
  type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
   call  pdifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=221345   ! code name
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
     CALL XERBLA( 'PDFEAST_GCSRGVX', -INFO+100 )
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

  call pfeast_distribution_type(N,isa,jsa,fpm(49),distype)
  !fpm(50)=distype



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!!! Perform a copy of the matrix values (transform into complex)
     nnza=isa(n+1)-1
     allocate(fsa(nnza))
     call ZLACP2( 'F', nnza, 1,sa , nnza, fsa, nnza )

     nnzb=isb(n+1)-1
     allocate(fsb(nnzb))
     call ZLACP2( 'F', nnzb, 1,sb , nnzb, fsb, nnzb ) 

!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,fsa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call zcsr_distribute_row(Ntotal,fsb,isb,jsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B

     deallocate(fsa)
     deallocate(fsb)
     
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
     ! simple copy

     nnzA=isa(n+1)-1
     allocate(matAj%isa(N+1))
     allocate(matAj%jsa(nnzA))
     allocate(matAj%sa(nnzA))
     call ZLACP2( 'F', nnza, 1,sa , nnza, matAj%sa, nnza )
     matAj%isa(1:n+1)= isa(1:n+1)
     matAj%jsa(1:nnzA)= jsa(1:nnzA)
     matAj%n=N
     matAj%nnz=nnzA
     !! for B
     nnzB=isb(n+1)-1
     allocate(matBj%isa(N+1))
     allocate(matBj%jsa(nnzB))
     allocate(matBj%sa(nnzB))
     call ZLACP2( 'F', nnzb, 1,sb , nnzb, matBj%sa, nnzb )
     matBj%isa(1:n+1)= isb(1:n+1)
     matBj%jsa(1:nnzB)= jsb(1:nnzB)
     matBj%n=N
     matBj%nnz=nnzB  

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
  allocate(isaz(Nlocal+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 !! get isaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  nnz=isaz(Nlocal+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !! get jsaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid!+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,1),isaz,jsaz)

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

     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo

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
  allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapA,matAj,matAjb,fpm(49),nb_procs3)
  allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapB,matBj,matBjb,fpm(49),nb_procs3)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision   
     allocate(cfsa(matAj%isa(Nlocal+1)-1))
     allocate(cfsb(matBj%isa(Nlocal+1)-1))
     call ZLAG2C(matAj%isa(Nlocal+1)-1,1,matAj%sa,matAj%isa(Nlocal+1)-1,cfsa,matAj%isa(Nlocal+1)-1,infoloc1)
     call ZLAG2C(matBj%isa(Nlocal+1)-1,1,matBj%sa,matBj%isa(Nlocal+1)-1,cfsb,matBj%isa(Nlocal+1)-1,infoloc1)     
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
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

     !! symbolic factorization ==< Does not work with MTYPE=13, PT ends up pointing to the same factorization?, also Ze choice not optimal
!!$  opt=3  !! get csaz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=Emid!+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call caddcsr(Nlocal,Ntotal,opt,-CONE,cfsa,matAj%isa,matAj%jsa,cmplx(Ze),cfsb,matBj%isa,matBj%jsa,csaz(1,i),isaz,jsaz)
!!$          if (nb_procs3>1) then
!!$          call  cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,M0,IPARM(1,i),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
!!$else
!!$ call  pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,i),isaz,jsaz,idum,M0,IPARM(1,i),MSGLVL,cwork,caux,infoloc2)
!!$   endif
!!$       enddo
!!$       


  else   !double precision   
     allocate(zaux(Nlocal,M0))

     !! Symbolic factorization
!!$       opt=3 !! get saz 
!!$ do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$    Ze=Emid!+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$    call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,i),isaz,jsaz)
!!$
!!$    if (nb_procs3>1) then
!!$call cluster_sparse_solver(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,M0,IPARM(1,i),MSGLVL,zwork,zaux,L3_COMM_WORLD,infoloc2)
!!$else
!!$call pardiso(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,saz(1,i),isaz,jsaz,idum,M0,IPARM(1,i),MSGLVL,zwork,zaux,infoloc2)
!!$end if
!!$enddo

  end if

!!$    if (infoloc2/=0) then
!!$       info=-2
!!$       return
!!$    end if
!!$
!!$   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call dfeast_grcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3
        PHASE=12 !include the symbolic factorization

        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(Nlocal,Ntotal,opt,-CONE,cfsa,matAj%isa,matAj%jsa,cmplx(Ze),cfsb,matBj%isa,matBj%jsa,csaz(1,id),isaz,jsaz) !! get csaz

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if


        else ! double precision fact- saz
           call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,id),isaz,jsaz) !! get saz    

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
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=0 ! normal solve

     
        !print *,id,sum(abs(zwork(1:N,1))),sum(abs(saz(1:100,id)))
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


!print *,sum(abs(zwork(1:N,1))),infoloc2!,nb_procs3
        
        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
        call zlbcsrmm(UPLO2,'N',Nsize,mapA,fpm(25),ZONE,matAjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapA,fpm(35),ZONE,matAjb,Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call zlbcsrmm(UPLO2,'N',Nsize,mapB,fpm(25),ZONE,matBjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapB,fpm(35),ZONE,matBjb,Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)



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
  deallocate(matAj%isa)
  deallocate(matAj%jsa)
  deallocate(matAj%sa)
  deallocate(matBj%isa)
  deallocate(matBj%jsa)
  deallocate(matBj%sa)
  deallocate(mapA)
  deallocate(mapB)
  deallocate(matAjb)
  deallocate(matBjb)


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

  !if (distype==0) then
  !   deallocate(fsa)
  !   deallocate(fsb)
  !end if



  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(cfsa)
     deallocate(cfsb)
  else
     deallocate(zaux)
     deallocate(saz)
  end if


#endif
!!!! done with MKL CLUSTER PARDISO!!!! done with MKL CLUSTER PARDISO

end subroutine pdfeast_gcsrgvx





subroutine pzfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
  !=====================================================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
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
  integer :: ijob,infoloc1,infoloc2,i,s,k,j,M02
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,Aq,Sq,zaux,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
!!!! matrix format
  !complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  integer :: opt,nnza,nnzb,nnz
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
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj,matBj
  type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
     call pzifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=241345   ! code name
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
     CALL XERBLA( 'PZFEAST_GCSRGVX', -INFO+100 )
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

  call pfeast_distribution_type(N,isa,jsa,fpm(49),distype)
  !fpm(50)=distype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!!! Perform a copy of the matrix values (transform into complex)
   !  nnza=isa(n+1)-1
   !  allocate(fsa(nnza))
   !  call ZCOPY(nnza,sa,1, fsa, 1 )
   !  nnzb=isb(n+1)-1
   !  allocate(fsb(nnzb))
   !  call ZCOPY( nnzb,sb,1, fsb, 1 )

!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call zcsr_distribute_row(Ntotal,sb,isb,jsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B

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

     nnzA=isa(n+1)-1
     allocate(matAj%isa(N+1))
     allocate(matAj%jsa(nnzA))
     allocate(matAj%sa(nnzA))
     call ZCOPY(nnza,sa,1, matAj%sa, 1 )
     matAj%isa(1:n+1)= isa(1:n+1)
     matAj%jsa(1:nnzA)= jsa(1:nnzA)
     matAj%n=N
     matAj%nnz=nnzA
     !! for B
     nnzB=isb(n+1)-1
     allocate(matBj%isa(N+1))
     allocate(matBj%jsa(nnzB))
     allocate(matBj%sa(nnzB))
     call ZCOPY(nnzb,sb,1, matBj%sa, 1 )
     matBj%isa(1:n+1)= isb(1:n+1)
     matBj%jsa(1:nnzB)= jsb(1:nnzB)
     matBj%n=N
     matBj%nnz=nnzB  

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

  allocate(isaz(Nlocal+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 !! get isaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  nnz=isaz(Nlocal+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !! get jsaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid!+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,1),isaz,jsaz)

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

     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo

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
  allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapA,matAj,matAjb,fpm(49),nb_procs3)
  allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapB,matBj,matBjb,fpm(49),nb_procs3)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision   
     allocate(cfsa(matAj%isa(Nlocal+1)-1))
     allocate(cfsb(matBj%isa(Nlocal+1)-1))
     call ZLAG2C(matAj%isa(Nlocal+1)-1,1,matAj%sa,matAj%isa(Nlocal+1)-1,cfsa,matAj%isa(Nlocal+1)-1,infoloc1)
     call ZLAG2C(matBj%isa(Nlocal+1)-1,1,matBj%sa,matBj%isa(Nlocal+1)-1,cfsb,matBj%isa(Nlocal+1)-1,infoloc1)     
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
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
     call dfeast_grcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3
        PHASE=12 !include the symbolic factorization

        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(Nlocal,Ntotal,opt,-CONE,cfsa,matAj%isa,matAj%jsa,cmplx(Ze),cfsb,matBj%isa,matBj%jsa,csaz(1,id),isaz,jsaz) !! get csaz

           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if


        else ! double precision fact- saz
           call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,id),isaz,jsaz) !! get saz    

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
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
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
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
        call zlbcsrmm(UPLO2,'N',Nsize,mapA,fpm(25),ZONE,matAjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapA,fpm(35),ZONE,matAjb,Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call zlbcsrmm(UPLO2,'N',Nsize,mapB,fpm(25),ZONE,matBjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call zlbcsrmm(UPLO2,'C',Nsize,mapB,fpm(35),ZONE,matBjb,Xj(1,fpm(34)),ZZERO,work(1,fpm(34)),fpm(49),nb_procs3)



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
  deallocate(matAj%isa)
  deallocate(matAj%jsa)
  deallocate(matAj%sa)
  deallocate(matBj%isa)
  deallocate(matBj%jsa)
  deallocate(matBj%sa)
  deallocate(mapA)
  deallocate(mapB)
  deallocate(matAjb)
  deallocate(matBjb)


  deallocate(Aq)
  deallocate(Sq)
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
     deallocate(cfsa)
     deallocate(cfsb)
  else
     deallocate(zaux)
     deallocate(saz)
  end if


#endif
!!!! done with MKL CLUSTER PARDISO!!!! done with MKL CLUSTER PARDISO

end subroutine pzfeast_gcsrgvx











subroutine pzfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B COMPLEX SYMMETRIC:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
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
  integer :: ijob,infoloc1,infoloc2,i,s,k,j
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,Aq,Sq,Xj
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: usa,usb
  integer,dimension(:), allocatable :: uisa,ujsa,uisb,ujsb
  integer :: opt,nnza,nnzb,nnz
!!!!!for cluster pardiso
  integer(8),dimension(:,:),allocatable :: pt 
  integer,dimension(:,:),allocatable :: iparm
  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: susa,susb !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB
  type zcsr
     integer ::n,m,nnz
     complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj,matBj
  type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb


  mkl=.false. ! initialization

#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to PIFEAST
      call pzifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel cluster mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef MKL  


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=241145  ! code name
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
     CALL XERBLA( 'PZFEAST_SCSRGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Find which distribution scheme is used for the input matrix/rhs: global or local
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call pfeast_distribution_type(N,isa,jsa,fpm(49),distype)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!! FORMAT CONVERSION TO CSR-UPPER     
     !! for A
     allocate(uisa(n+1))
     nnzA=isa(n+1)-1
     if ((UPLO=='F').or.(UPLO=='f')) nnzA=nnzA/2+n ! slightly overestimated
     allocate(ujsa(nnzA))
     allocate(usa(nnzA))
     call zcsr_convert_upper(n,UPLO,sa,isa,jsa,usa,uisa,ujsa)
     !! for B
     allocate(uisb(n+1))
     nnzB=isb(n+1)-1
     if ((UPLO=='F').or.(UPLO=='f')) nnzB=nnzB/2+n ! slightly overestimated
     allocate(ujsb(nnzB))
     allocate(usb(nnzB))
     call zcsr_convert_upper(n,UPLO,sb,isb,jsb,usb,uisb,ujsb)

!!! DISTRIBUTE the (upper) sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,usa,uisa,ujsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call zcsr_distribute_row(Ntotal,usb,uisb,ujsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B


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

!!! FORM the row distributed upper matrix for Pardiso   
     if ((UPLO=='F').or.(UPLO=='f')) then !matrix need to be reduced to upper part
        !! for A
        nnzA=isa(n+1)-1 ! overestimation
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%isa(1)=1
        do i=1,N
           j=0
           do k=isa(i),isa(i+1)-1
              if (jsa(k)>=i+startj-1) then
                 j=j+1
                 matAj%jsa(matAj%isa(i)+j-1)=jsa(k)
                 matAj%sa(matAj%isa(i)+j-1)=sa(k)
              end if
           enddo
           matAj%isa(i+1)=matAj%isa(i)+j 
        enddo
        matAj%n=N
        matAj%nnz=matAj%isa(n+1)-1
        !! for B
        nnzB=isb(n+1)-1 ! overestimation
        allocate(matBj%isa(N+1))
        allocate(matBj%jsa(nnzB))
        allocate(matBj%sa(nnzB))
        matBj%isa(1)=1
        do i=1,N
           j=0
           do k=isb(i),isb(i+1)-1
              if (jsb(k)>=i+startj-1) then
                 j=j+1
                 matBj%jsa(matBj%isa(i)+j-1)=jsb(k)
                 matBj%sa(matBj%isa(i)+j-1)=sb(k)
              end if
           enddo
           matBj%isa(i+1)=matBj%isa(i)+j 
        enddo
        matBj%n=N
        matBj%nnz=matBj%isa(n+1)-1

     else ! include the following cases
        !! ((UPLO=='U').or.(UPLO=='u')) !! upper-csr already -- simple copy
        !! ((UPLO=='L').or.(UPLO=='l')) !! temp-copy !! transpose difficult..could be done once the local blocks csr are defined below

        !! for A
        nnzA=isa(n+1)-1
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%sa(1:nnzA)=sa(1:nnzA)
        matAj%isa(1:n+1)= isa(1:n+1)
        matAj%jsa(1:nnzA)= jsa(1:nnzA)
        matAj%n=N
        matAj%nnz=nnzA
        !! for B
        nnzB=isb(n+1)-1
        allocate(matBj%isa(N+1))
        allocate(matBj%jsa(nnzB))
        allocate(matBj%sa(nnzB))
        matBj%sa(1:nnzB)=sb(1:nnzB)
        matBj%isa(1:n+1)= isb(1:n+1)
        matBj%jsa(1:nnzB)= jsb(1:nnzB)
        matBj%n=N
        matBj%nnz=nnzB    
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

  allocate(isaz(Nlocal+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 !! get isaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  nnz=isaz(Nlocal+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !! get jsaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid !ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,1),isaz,jsaz)

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
     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo

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
  allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapA,matAj,matAjb,fpm(49),nb_procs3)
  allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapB,matBj,matBjb,fpm(49),nb_procs3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! HANDLE THE PARTICULAR CASE of user distributed matrices via Lower format (get it upper for pardiso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (distype==1) then
     if ((UPLO=='L').or.(UPLO=='l')) then !! perform the transpose
        !! create local block upper format

        call zlbcsr_transpose(Nsize,mapA,matAjb,L3_COMM_WORLD,nb_procs3)
        call zlbcsr_transpose(Nsize,mapB,matBjb,L3_COMM_WORLD,nb_procs3)

        !! convert into row distributed csr format
        deallocate(matAj%jsa)
        deallocate(matAj%sa)
        deallocate(matAj%isa)
        call zlbcsr_distribute_row(Nsize,mapA,matAjb,matAj,L3_COMM_WORLD,nb_procs3)
        deallocate(matBj%jsa)
        deallocate(matBj%sa)
        deallocate(matBj%isa)
        call zlbcsr_distribute_row(Nsize,mapB,matBjb,matBj,L3_COMM_WORLD,nb_procs3)


!!! we need to reconstruct isaz and jsaz, as upper local rows
        allocate(saz(1,1)) ! dummy
        opt=1 !! get isaz
        call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
        nnz=isaz(Nlocal+1)-1
        deallocate(jsaz)
        allocate(jsaz(nnz))
        opt=2 !! get jsaz
        call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,saz,isaz,jsaz)
        deallocate(saz) ! dummy
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision   
     allocate(susa(matAj%isa(Nlocal+1)-1))
     allocate(susb(matBj%isa(Nlocal+1)-1))
     call ZLAG2C(matAj%isa(Nlocal+1)-1,1,matAj%sa,matAj%isa(Nlocal+1)-1,susa,matAj%isa(Nlocal+1)-1,infoloc1)
     call ZLAG2C(matBj%isa(Nlocal+1)-1,1,matBj%sa,matBj%isa(Nlocal+1)-1,susb,matBj%isa(Nlocal+1)-1,infoloc1)     
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
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
     call zfeast_srcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization

        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(Nlocal,Ntotal,opt,-CONE,susa,matAj%isa,matAj%jsa,cmplx(Ze),susb,matBj%isa,matBj%jsa,csaz(1,id),isaz,jsaz) !! get csaz
           if (nb_procs3>1) then
              call cluster_sparse_solver(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,L3_COMM_WORLD,infoloc2)
           else
              call pardiso(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,Ntotal,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           end if

        else ! double precision fact- saz
           call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,saz(1,id),isaz,jsaz) !! get saz
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
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call zlbcsrmm(UPLO2,'N',Nsize,mapA,fpm(25),ZONE,matAjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call zlbcsrmm(UPLO2,'N',Nsize,mapB,fpm(25),ZONE,matBjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)


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
  deallocate(matAj%isa)
  deallocate(matAj%jsa)
  deallocate(matAj%sa)
  deallocate(matBj%isa)
  deallocate(matBj%jsa)
  deallocate(matBj%sa)
  deallocate(mapA)
  deallocate(mapB)
  deallocate(matAjb)
  deallocate(matBjb)


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

  if (distype==0) then
     deallocate(usa)
     deallocate(uisa)
     deallocate(ujsa)

     deallocate(usb)
     deallocate(uisb)
     deallocate(ujsb)
  end if

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(susa)
     deallocate(susb)
  else
     deallocate(zaux)
     deallocate(saz)
  endif



#endif
!!!! done with MKL CLUSTER PARDISO

end subroutine pzfeast_scsrgvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pdfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        REAL DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) REAL DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! James Kestyn 2016-2018
  ! Eric Polizzi 2019
  !===================================================================== 
  implicit none
  !    include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  double precision,dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pdfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=221144
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return

!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
       call pdfeast_scsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,loop,fpm(1),fpm(9),fpm(40),fpm(49),Emin,Emax,M0,mode,info)
      if (info/=0) then
              info=7
              return  ! failed search
      endif
  end if
!!!!!!!!!!!!!!  
    
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pdfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdfeast_scsrgv




subroutine pzfeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
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
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pzfeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241244
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return

!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
       call pzfeast_hcsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,loop,fpm(1),fpm(9),fpm(40),fpm(49),Emin,Emax,M0,mode,info)
      if (info/=0) then
              info=7
              return  ! failed search
      endif
  end if
!!!!!!!!!!!!!!  
  
  
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))    
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzfeast_hcsrgv



subroutine pdfeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        REAL DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
  !=====================================================================
  implicit none
  !     include 'f90_noruntime_interface.fi'
  integer :: N
  double precision,dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
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
!!!    Wrapper Routine to expert routine: pdfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
 
  fpm(30)=221344
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pdfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdfeast_gcsrgv


!!$

subroutine pzfeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
  !=====================================================================
  implicit none
  !    include 'f90_noruntime_interface.fi'
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
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
!!!    Wrapper Routine to expert routine: pzfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241344
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzfeast_gcsrgv


!!$

subroutine pzfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX SYMMETRIC:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  sb         (input)        COMPLEX DOUBLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
  !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
  !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
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
  !=====================================================================
  implicit none
  !    include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
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
!!!    Wrapper Routine to expert routine: pzfeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241144
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pzfeast_scsrgv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!$

subroutine pdfeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) REAL DOUBLE PRECISION(N,M0) : 
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
  ! James Kestyn 2016-2018
  ! Eric Polizzi 2019
  !============================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  double precision,dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=221143
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pdfeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call feastdefault(fpm,info)
  if (info/=0) return

  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n ! n is global or local
     sb(i)=1.0d0
     jsb(i)=i+startj-1 ! account for local distribution if any 
     isb(i)=i
  enddo
  isb(n+1)=n+1


  call pdfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pdfeast_scsrevx





subroutine pzfeast_hcsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)
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
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  character(len=1) :: UPLO
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pzfeast_hcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241243
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  
  call feastdefault(fpm,info)
  if (info/=0) return

  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if


  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n! n is global or local
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i+startj-1 ! account for local distribution if any 
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call pzfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pzfeast_hcsrevx



subroutine pdfeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
  !=====================================================================
  implicit none
  !    include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  integer :: N
  double precision,dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pdfeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=221343
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))

  do i=1,n ! n is global or local
     sb(i)=1.0d0
     jsb(i)=i+startj-1 ! account for local distribution if any 
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call pdfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)


end subroutine pdfeast_gcsrevx




subroutine pzfeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
  !Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pzfeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241343
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return

  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n  ! n is global or local
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i+startj-1 ! account for local distribution if any 
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call pzfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pzfeast_gcsrevx



subroutine pzfeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  !Eric Polizzi 2019
  !=====================================================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  character(len=1) :: UPLO
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pzfeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241143
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return

  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n ! n is global or local
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i+startj-1 ! account for local distribution if any 
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call pzfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pzfeast_scsrevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine pdfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  X          (input/output) REAL DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! James Kestyn 2016-2018
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !    include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  double precision,dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pdfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. 
#ifdef MKL
  mkl=.true.
#endif

  fpm(30)=221142
   if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return


!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
     call pdfeast_scsrev_search(UPLO,N,sa,isa,jsa,loop,fpm(1),fpm(9),fpm(40),fpm(49),Emin,Emax,M0,mode,info)
      if (info/=0) then
              info=7 
              return  ! stochastic estimate
      endif
  end if
!!!!!!!!!!!!!!  
  
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour

  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if


  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n ! n is global or local
     sb(i)=1.0d0
     jsb(i)=i+startj-1 ! account for local distribution if any 
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call pdfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pdfeast_scsrev








subroutine pzfeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
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
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  character(len=1) :: UPLO
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pzfeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. 
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241242
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return

!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
     call pzfeast_hcsrev_search(UPLO,N,sa,isa,jsa,loop,fpm(1),fpm(9),fpm(40),fpm(49),Emin,Emax,M0,mode,info)
      if (info/=0) then
              info=7
              return  ! failed search
      endif
  end if
!!!!!!!!!!!!!!   
  
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour


  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if



  ! identity B matrix- option for standard eigenvalue problem          
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i+startj-1 ! account for local distribution if any
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call pzfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzfeast_hcsrev




subroutine pdfeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
  !=====================================================================
  implicit none
  !    include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  integer :: N,LDA
  double precision,dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
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
!!!    Wrapper Routine to expert routine: pdfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=221342
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour

  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=1.0d0
     jsb(i)=i+startj-1 ! account for local distribution if any
     isb(i)=i
  enddo
  isb(n+1)=n+1
  call pdfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pdfeast_gcsrev




subroutine pzfeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
  !Eric Polizzi 2019
  !=====================================================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
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
!!!    Wrapper Routine to expert routine: pzfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code    
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241342
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if


  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i+startj-1 ! account for local distribution if any
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call pzfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine pzfeast_gcsrev




subroutine pzfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: SPARSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
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
  !Eric Polizzi 2019
  !=====================================================================
  implicit none
  !   include 'f90_noruntime_interface.fi'
  include 'mpif.h'
  character(len=1) :: UPLO
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
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
!!!    Wrapper Routine to expert routine: pzfeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
 logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization !!<<?? to do for all routines
#ifdef MKL
  mkl=.true.
#endif
  
  fpm(30)=241142
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour

  !is matrix/rhs distribution global or local?
  L3_COMM_WORLD=fpm(49)
  call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  rank3=0 ! default
  nb_procs3=1 ! default
  startj=1 ! default
  if (distype==1) then ! local
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
     ! obtain global Nsize
     allocate(Nsize(nb_procs3))
     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     ! get local startj
     do i=1,rank3
        startj=startj+Nsize(i)
     enddo
     deallocate(Nsize)
  end if


  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i+startj-1 ! account for local distribution if any
     isb(i)=i
  enddo
  isb(n+1)=n+1
  call pzfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine pzfeast_scsrev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!5555555555555555555555555555555555555555555555555555555555555555555555555555555
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pdfeast_scsrev_search(UPLO,N,sa,isa,jsa,itloop,fpm1,fpm9,fpm40,fpm49,Emin,Emax,M0,M,info)
 !  Purpose 
  !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Standard eigenvalue problem version: AX=XE  where A is real symmetric 
  ! 
  !  DOUBLE PRECISION MPI version
  !
  !  Remark: The code starts by implementing few step of Arnoldi. Few options are hardcoded.
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  itloop     (output)       INTEGER- number of stochastic search needed 
  !  fpm1       (input)        INTEGER- Feast parameter fpm(1)
  !  fpm9       (input)        INTEGER- Feast parameter fpm(9)
  !  fpm40      (input)        INTEGER- Feast parameter fpm(40)
  !                                     if fpm(40)=-1 the M0/2 lowest eigenvalues are returned
  !                                     if fpm(40)=1  the M0/2 largest eigenvalues are returned
  !  fpm49      (input)        INTEGER- Feast parameter fpm(49)
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ==================================================================== 
  
  implicit none
    !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40,fpm49
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  double precision,dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: L2_COMM_WORLD,L3_COMM_WORLD,rank2,rank3,code,nb_procs2,nb_procs3,Ntotal
  character(len=3) :: ctemp
  integer(8) :: fout
  !-----
   type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: matAj
   type(dcsr), dimension(:,:),allocatable :: matAjb
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA 
 integer:: Nlocal,startj,endj,distype,nnzA,infofeast
!-----

info=0 

  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1    

  !----------------------------------------------
  L2_COMM_WORLD=fpm9
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm49
  call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
  call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  !----------------------

  !! local or global?
 call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 if (distype==0) then ! global- matrix globally known
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    Ntotal=N
 !!! DISTRIBUTE the  sparse matrix by rows among all procs 
     allocate(Nsize(nb_procs3))
     call dcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,L3_COMM_WORLD) !! A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else ! local! matrix locally known by row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   Nlocal=N
  ! get Ntotal
   Ntotal=0
   call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
 !! obtain global Nsize
     allocate(Nsize(nb_procs3))

     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)

 !! for A (simple copy for consistency with the rest of the code)
        nnzA=isa(n+1)-1
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%sa(1:nnzA)=sa(1:nnzA)
        matAj%isa(1:n+1)= isa(1:n+1)
        matAj%jsa(1:nnzA)= jsa(1:nnzA)
        matAj%n=N
        matAj%nnz=nnzA
end if

!! DISTRIBUTE sparse row into block csr format (needed for matvec with uplo=l,u,f in Arnoldi)
allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
 call dget_local_block_csr(Nsize,mapA,matAj,matAjb,L3_COMM_WORLD,nb_procs3)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call dcheck_feast_srci_input(0.0d0,1.0d0,M0,Ntotal,info)
 if (info/=0) return
 
  ! comments?
  com=.false. ! default
  if ((fpm1/=0).and.(rank2==0).and.(rank3==0)) com=.true. !only rank0 is commenting

  ! file or screen output?
  if (com) then
   if (fpm1<0) then
           write(ctemp,'(I3)') abs(fpm1)
           fout=abs(fpm1)+200 ! file number default
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append')
           write(fout,*)  
        elseif (fpm1>0) then
           fout=6 !screen
   end if
  end if
  ! print header   
    if (com ) then   
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') '***** FEAST Search Estimate for Emin-Emax *****'
     write(fout,'(A)') '*** looking for the M0/2 extremal eigenvalue **' 
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') ' Routine pdfeast_scsrev_search'
     write(fout,'(A)',advance='no') ' Subspace M0='
     write(fout,'(I4)') M0
     write(fout,'(A)',advance='no') ' Looking for '
     write(fout,'(I4)',advance='no') max(M0/2,1)
     if (fpm40<0) then
        write(fout,'(A)') ' lowest eigenvalues '
     else
        write(fout,'(A)') ' largest eigenvalues '
     end if
     write(fout,*)
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Search the edge of the spectrum using few steps of Arnoldi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  eps=1d-2    
  itmax=min(50,Ntotal) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge)
  allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(Nlocal,itmax)) 
  
!!!!! Search the edge of the spectrum using few steps of Arnoldi
  it=itmax
  call pdarnoldi(UPLO,Nsize,mapA,matAjb,dE,dX,dres,eps,it,L3_COMM_WORLD,nb_procs3)
 
!---
if (com) then
   write(fout,'(A)',advance='no') ' Arnoldi results (Extremal eig. estimated after '
   write(fout,'(I4)',advance='no') it
 write(fout,'(A)') ' steps):'
   Do i=1,it
      if ((i==1).or.(i==it)) then !! add comments for printing them all 
      write(fout,'(I4)',advance='no') i
      write(fout,'(ES25.16)',advance='no') dE(i)
      if (i==1) then
         write(fout,'(A)',advance='no') ' (lowest)  res='
         write(fout,'(ES25.16)',advance='no') dres(i)
   elseif (i==it) then
      write(fout,'(A)',advance='no') ' (largest) res='
      write(fout,'(ES25.16)',advance='no') dres(i)
   end if
   write(fout,*)
end if !!
end Do
end if 

if ((it==itmax).and.(dres(1)>eps).and.(dres(it)>eps)) then
!if (it==itmax) then
   info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
          !! search has failed (ill-conditionned matrix?)
   return
end if

itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
Emin=dE(min(3,Ntotal))-(dE(min(3,Ntotal))-dE(1))*1.05d0
Emax=dE(itmax-min(Ntotal-1,2))+(dE(itmax)-dE(itmax-min(Ntotal-1,2)))*1.05d0

if (com) then
write(fout,*)
 write(fout,'(A)') ' Start FEAST stochastic search:'
 write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
 i=0
 write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,Ntotal
end if


!!!!!!! Find a more appropriate Emax or Emin to start with
  if (fpm40 < 0) then ! search lowest- change Emax
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=dE(2)-Emin
     else
        diff=dE(min(4,n-1))-Emin !dE(1) has converged faster - eig more compact close to E1
     end if
     Emax = Emin+diff

  else ! search largest - change Emin
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=Emax-dE(itmax-min(3,n-1))
     else
        diff=Emax-dE(itmax-1) !dE(1) has converged faster - eig more compact close to E1
     end if

     Emin = Emax-diff
  endif


  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!! Start the stochastic search process
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  MW = max(M0/2,1) ! number of wanted eigenvalue
  call feastinit(fpm)
  
 MS=min(10,Ntotal) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS)) ! local 

  
  fpm(9)=fpm9
  fpm(49)=fpm49

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1
  
  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1
    
     call pdfeast_scsrev(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)
   
       if (infofeast/=5) then ! stochastic estimate
     infofeast=71
     return
  end if
     
     if (com)  write(fout,'(I2,X,2ES25.16,X,I7)') itloop,Emin,Emax,M
  
  
  if (init) then 
     if (M <Mw) then ! not enough eigenvalue found (initial Emin-Emax not well set-up)
        if (fpm40 <0) then
           Emax=Emin+(itloop+1)*diff
        else
           Emin=Emax-(itloop+1)*diff
        endif
     else
        init=.false.
   endif     
end if
   
   if ((M > M0).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax - diff
        else
           Emin = Emin + diff
        endif

     else if ((M < MW).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax + diff
        else
           Emin = Emin - diff
        endif
     endif

  enddo

 
  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Start the exact search (using IFEAST)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  MS=min(M0,Ntotal) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS))

  call feastinit(fpm)
  
  fpm(9)=fpm9
  fpm(49)=fpm49
  
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
  call pdfeast_scsrev(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

  if (infofeast/=0) then
     info=72
     return
  end if
     
  
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,M,' (non-stochastic)'

  if (fpm40<0) then
  Emax = (dE(MW)+dE(MW+1))/2.0d0
else
   if (M>MW) then
      Emin=(dE(M-MW)+dE(M-MW+1))/2.0d0
   else
      Emin=(dE(M0-(MW-M))+dE(M0-(MW-M)+1))/2.0d0 
   endif
end if
 
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,MW,' (adjusted)'

  
if ((com).and.(fpm1<0)) close(fout)


 deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)


end subroutine pdfeast_scsrev_search





subroutine pzfeast_hcsrev_search(UPLO,N,sa,isa,jsa,itloop,fpm1,fpm9,fpm40,fpm49,Emin,Emax,M0,M,info)
   !  Purpose 
   !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Standard eigenvalue problem version: AX=XE  where A is complex Hermitian 
  ! 
  !  DOUBLE PRECISION MPI version
  !
  !  Remark: The code starts by implementing few step of Arnoldi. Few options are hardcoded.
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !
  !  N          (input)        INTEGER: Size system
  !  sa         (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  itloop     (output)       INTEGER- number of stochastic search needed 
  !  fpm1       (input)        INTEGER- Feast parameter fpm(1)
  !  fpm9       (input)        INTEGER- Feast parameter fpm(9)
  !  fpm40      (input)        INTEGER- Feast parameter fpm(40)
  !                                     if fpm(40)=-1 the M0/2 lowest eigenvalues are returned
  !                                     if fpm(40)=1  the M0/2 largest eigenvalues are returned
  !  fpm49       (input)       INTEGER- Feast parameter fpm(49)
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ==================================================================== 
  implicit none
    !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40,fpm49
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: L2_COMM_WORLD,L3_COMM_WORLD,rank2,rank3,code,nb_procs2,nb_procs3,Ntotal
  character(len=3) :: ctemp
  integer(8) :: fout
  !-----
 type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj
   type(zcsr), dimension(:,:),allocatable :: matAjb
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA 
 integer:: Nlocal,startj,endj,distype,nnzA,infofeast
!-----

   
info=0 

  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1    

  !----------------------------------------------
  L2_COMM_WORLD=fpm9
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm49
  call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
  call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  !----------------------

  !! local or global?
 call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 if (distype==0) then ! global- matrix globally known
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    Ntotal=N
 !!! DISTRIBUTE the  sparse matrix by rows among all procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,L3_COMM_WORLD) !! A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else ! local! matrix locally known by row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   Nlocal=N
  ! get Ntotal
   Ntotal=0
   call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
 !! obtain global Nsize
     allocate(Nsize(nb_procs3))

     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)

 !! for A (simple copy for consistency with the rest of the code)
        nnzA=isa(n+1)-1
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%sa(1:nnzA)=sa(1:nnzA)
        matAj%isa(1:n+1)= isa(1:n+1)
        matAj%jsa(1:nnzA)= jsa(1:nnzA)
        matAj%n=N
        matAj%nnz=nnzA
end if

!! DISTRIBUTE sparse row into block csr format (needed for matvec with uplo=l,u,f in Arnoldi)
allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
 call zget_local_block_csr(Nsize,mapA,matAj,matAjb,L3_COMM_WORLD,nb_procs3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call dcheck_feast_srci_input(0.0d0,1.0d0,M0,Ntotal,info)
 if (info/=0) return
 
  ! comments?
  com=.false. ! default
  if ((fpm1/=0).and.(rank2==0).and.(rank3==0)) com=.true. !only rank0 is commenting

  ! file or screen output?
  if (com) then
   if (fpm1<0) then
           write(ctemp,'(I3)') abs(fpm1)
           fout=abs(fpm1)+200 ! file number default
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append')
           write(fout,*)  
        elseif (fpm1>0) then
           fout=6 !screen
        end if
     end if
  ! print header   
    if (com ) then   
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') '***** FEAST Search Estimate for Emin-Emax *****'
     write(fout,'(A)') '*** looking for the M0/2 extremal eigenvalue **' 
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') ' Routine pzfeast_hcsrev_search'
     write(fout,'(A)',advance='no') ' Subspace M0='
     write(fout,'(I4)') M0
     write(fout,'(A)',advance='no') ' Looking for '
     write(fout,'(I4)',advance='no') max(M0/2,1)
     if (fpm40<0) then
        write(fout,'(A)') ' lowest eigenvalues '
     else
        write(fout,'(A)') ' largest eigenvalues '
     end if
     write(fout,*)
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Search the edge of the spectrum using few steps of Arnoldi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  eps=1d-2    
  itmax=min(50,Ntotal) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge)
  allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(Nlocal,itmax)) 
  
!!!!! Search the edge of the spectrum using few steps of Arnoldi
  it=itmax
  call pzharnoldi(UPLO,Nsize,mapA,matAjb,dE,dX,dres,eps,it,L3_COMM_WORLD,nb_procs3)
 
!---
if (com) then
   write(fout,'(A)',advance='no') ' Arnoldi results (Extremal eig. estimated after '
   write(fout,'(I4)',advance='no') it
 write(fout,'(A)') ' steps):'
   Do i=1,it
      if ((i==1).or.(i==it)) then !! add comments for printing them all 
      write(fout,'(I4)',advance='no') i
      write(fout,'(ES25.16)',advance='no') dE(i)
      if (i==1) then
         write(fout,'(A)',advance='no') ' (lowest)  res='
         write(fout,'(ES25.16)',advance='no') dres(i)
   elseif (i==it) then
      write(fout,'(A)',advance='no') ' (largest) res='
      write(fout,'(ES25.16)',advance='no') dres(i)
   end if
   write(fout,*)
end if !!
end Do
end if 

if ((it==itmax).and.(dres(1)>eps).and.(dres(it)>eps)) then
!if (it==itmax) then
   info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
          !! search has failed (ill-conditionned matrix?)
   return
end if

itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
Emin=dE(min(3,Ntotal))-(dE(min(3,Ntotal))-dE(1))*1.05d0
Emax=dE(itmax-min(Ntotal-1,2))+(dE(itmax)-dE(itmax-min(Ntotal-1,2)))*1.05d0

if (com) then
write(fout,*)
 write(fout,'(A)') ' Start FEAST stochastic search:'
 write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
 i=0
 write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,Ntotal
end if


!!!!!!! Find a more appropriate Emax or Emin to start with
  if (fpm40 < 0) then ! search lowest- change Emax
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=dE(2)-Emin
     else
        diff=dE(min(4,n-1))-Emin !dE(1) has converged faster - eig more compact close to E1
     end if
     Emax = Emin+diff

  else ! search largest - change Emin
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=Emax-dE(itmax-min(3,n-1))
     else
        diff=Emax-dE(itmax-1) !dE(1) has converged faster - eig more compact close to E1
     end if

     Emin = Emax-diff
  endif


  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!! Start the stochastic search process
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  MW = max(M0/2,1) ! number of wanted eigenvalue
  call feastinit(fpm)

  
  
 MS=min(10,Ntotal) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS)) ! local 
  
  fpm(9)=fpm9
  fpm(49)=fpm49

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1

  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1
     call pzfeast_hcsrev(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

       if (infofeast/=5) then ! stochastic estimate
     infofeast=71
     return
  end if
     
     if (com)  write(fout,'(I2,X,2ES25.16,X,I7)') itloop,Emin,Emax,M
    
  if (init) then 
     if (M <Mw) then ! not enough eigenvalue found (initial Emin-Emax not well set-up)
        if (fpm40 <0) then
           Emax=Emin+(itloop+1)*diff
        else
           Emin=Emax-(itloop+1)*diff
        endif
     else
        init=.false.
   endif     
end if
   
   if ((M > M0).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax - diff
        else
           Emin = Emin + diff
        endif

     else if ((M < MW).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax + diff
        else
           Emin = Emin - diff
        endif
     endif

  enddo

 
  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Start the exact search (using IFEAST)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  MS=min(M0,Ntotal) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS))

  call feastinit(fpm)
  
  fpm(9)=fpm9
  fpm(49)=fpm49
  
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
  call pzfeast_hcsrev(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

  if (infofeast/=0) then
     info=72
     return
  end if
     
  
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,M,' (non-stochastic)'

  if (fpm40<0) then
  Emax = (dE(MW)+dE(MW+1))/2.0d0
else
   if (M>MW) then
      Emin=(dE(M-MW)+dE(M-MW+1))/2.0d0
   else
      Emin=(dE(M0-(MW-M))+dE(M0-(MW-M)+1))/2.0d0 
   endif
end if
 
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,MW,' (adjusted)'

  
if ((com).and.(fpm1<0)) close(fout)


 deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)


end subroutine pzfeast_hcsrev_search






subroutine pdfeast_scsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,itloop,fpm1,fpm9,fpm40,fpm49,Emin,Emax,M0,M,info)
   !  Purpose 
  !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Generalized eigenvalue problem version: AX=BXE  where A is real symmetric, B is spd 
  ! 
  !  DOUBLE PRECISION MPI version
  !
  !  Remark: The code starts by implementing few step of Arnoldi. Few options are hardcoded.
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !
  !  N          (input)        INTEGER: Size system
  !  sa,sb      (input)        REAL DOUBLE PRECISION (isa(N+1)-1):  Matrix A,B- CSR format 
  !  isa,isb    (input)        INTEGER(N+1): CSR row array of Matrix A,B
  !  jsa,jsb    (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A,B
  !  itloop     (output)       INTEGER- number of stochastic search needed 
  !  fpm1       (input)        INTEGER- Feast parameter fpm(1)
  !  fpm9       (input)        INTEGER- Feast parameter fpm(9)
  !  fpm40      (input)        INTEGER- Feast parameter fpm(40)
  !                                     if fpm(40)=-1 the M0/2 lowest eigenvalues are returned
  !                                     if fpm(40)=1  the M0/2 largest eigenvalues are returned
  !  fpm49       (input)        INTEGER- Feast parameter fpm(9)
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================   
  implicit none
    !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40,fpm49
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  double precision,dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: L2_COMM_WORLD,L3_COMM_WORLD,rank2,rank3,code,nb_procs2,nb_procs3,Ntotal
  character(len=3) :: ctemp
  integer(8) :: fout
  !-----
   type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: matAj,matBj
   type(dcsr), dimension(:,:),allocatable :: matAjb,matBjb
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB 
 integer:: Nlocal,startj,endj,distype,nnzA,nnzB,infofeast
!-----

   
info=0 

  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1    

  !----------------------------------------------
  L2_COMM_WORLD=fpm9
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm49
  call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
  call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  !----------------------

  !! local or global?
 call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 if (distype==0) then ! global- matrix globally known
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    Ntotal=N
 !!! DISTRIBUTE the  sparse matrix by rows among all procs 
     allocate(Nsize(nb_procs3))
     call dcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,L3_COMM_WORLD) !! A
      call dcsr_distribute_row(Ntotal,sb,isb,jsb,matBj,startj,endj,Nlocal,Nsize,L3_COMM_WORLD) !! B

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else ! local! matrix locally known by row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   Nlocal=N
  ! get Ntotal
   Ntotal=0
   call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
 !! obtain global Nsize
     allocate(Nsize(nb_procs3))

     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)

 !! for A (simple copy for consistency with the rest of the code)
        nnzA=isa(n+1)-1
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%sa(1:nnzA)=sa(1:nnzA)
        matAj%isa(1:n+1)= isa(1:n+1)
        matAj%jsa(1:nnzA)= jsa(1:nnzA)
        matAj%n=N
        matAj%nnz=nnzA
 !! for B
        nnzB=isb(n+1)-1
        allocate(matBj%isa(N+1))
        allocate(matBj%jsa(nnzB))
        allocate(matBj%sa(nnzB))
        matBj%sa(1:nnzB)=sb(1:nnzB)
        matBj%isa(1:n+1)= isb(1:n+1)
        matBj%jsa(1:nnzB)= jsb(1:nnzB)
        matBj%n=N
        matBj%nnz=nnzB
        
end if

!! DISTRIBUTE sparse row into block csr format (needed for matvec with uplo=l,u,f in Arnoldi)
allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
 call dget_local_block_csr(Nsize,mapA,matAj,matAjb,L3_COMM_WORLD,nb_procs3)
 allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call dget_local_block_csr(Nsize,mapB,matBj,matBjb,L3_COMM_WORLD,nb_procs3)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call dcheck_feast_srci_input(0.0d0,1.0d0,M0,Ntotal,info)
 if (info/=0) return
 
  ! comments?
  com=.false. ! default
  if ((fpm1/=0).and.(rank2==0).and.(rank3==0)) com=.true. !only rank0 is commenting

  ! file or screen output?
  if (com) then
   if (fpm1<0) then
           write(ctemp,'(I3)') abs(fpm1)
           fout=abs(fpm1)+200 ! file number default
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append')
           write(fout,*)  
        elseif (fpm1>0) then
           fout=6 !screen
        end if
     end if
  ! print header   
    if (com ) then   
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') '***** FEAST Search Estimate for Emin-Emax *****'
     write(fout,'(A)') '*** looking for the M0/2 extremal eigenvalue **' 
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') ' Routine pdfeast_scsrgv_search'
     write(fout,'(A)',advance='no') ' Subspace M0='
     write(fout,'(I4)') M0
     write(fout,'(A)',advance='no') ' Looking for '
     write(fout,'(I4)',advance='no') max(M0/2,1)
     if (fpm40<0) then
        write(fout,'(A)') ' lowest eigenvalues '
     else
        write(fout,'(A)') ' largest eigenvalues '
     end if
     write(fout,*)
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Search the edge of the spectrum using few steps of Arnoldi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  eps=1d-2    
  itmax=min(50,Ntotal) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge)
  allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(Nlocal,itmax)) 
  
!!!!! Search the edge of the spectrum using few steps of Arnoldi
  it=itmax

call pdgarnoldi(UPLO,Nsize,mapA,matAjb,mapB,matBjb,dE,dX,dres,eps,it,1d-3,30,L3_COMM_WORLD,nb_procs3) ! 2 constant arguments for bicgstab for B matrix

!---
if (com) then
   write(fout,'(A)',advance='no') ' Arnoldi results (Extremal eig. estimated after '
   write(fout,'(I4)',advance='no') it
 write(fout,'(A)') ' steps):'
   Do i=1,it
      if ((i==1).or.(i==it)) then !! add comments for printing them all 
      write(fout,'(I4)',advance='no') i
      write(fout,'(ES25.16)',advance='no') dE(i)
      if (i==1) then
         write(fout,'(A)',advance='no') ' (lowest)  res='
         write(fout,'(ES25.16)',advance='no') dres(i)
   elseif (i==it) then
      write(fout,'(A)',advance='no') ' (largest) res='
      write(fout,'(ES25.16)',advance='no') dres(i)
   end if
   write(fout,*)
end if !!
end Do
end if 

if ((it==itmax).and.(dres(1)>eps).and.(dres(it)>eps)) then
!if (it==itmax) then
   info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
          !! search has failed (ill-conditionned matrix?)
   return
end if

itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
Emin=dE(min(3,Ntotal))-(dE(min(3,Ntotal))-dE(1))*1.05d0
Emax=dE(itmax-min(Ntotal-1,2))+(dE(itmax)-dE(itmax-min(Ntotal-1,2)))*1.05d0

if (com) then
write(fout,*)
 write(fout,'(A)') ' Start FEAST stochastic search:'
 write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
 i=0
 write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,Ntotal
end if


!!!!!!! Find a more appropriate Emax or Emin to start with
  if (fpm40 < 0) then ! search lowest- change Emax
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=dE(2)-Emin
     else
        diff=dE(min(4,n-1))-Emin !dE(1) has converged faster - eig more compact close to E1
     end if
     Emax = Emin+diff

  else ! search largest - change Emin
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=Emax-dE(itmax-min(3,n-1))
     else
        diff=Emax-dE(itmax-1) !dE(1) has converged faster - eig more compact close to E1
     end if

     Emin = Emax-diff
  endif


  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!! Start the stochastic search process
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  MW = max(M0/2,1) ! number of wanted eigenvalue
  call feastinit(fpm)
  
 MS=min(10,Ntotal) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS)) ! local 

  
  fpm(9)=fpm9
  fpm(49)=fpm49

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1
  
  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1
    
     call pdfeast_scsrgv(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,matBj%sa,matBj%isa,matBj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)
   
       if (infofeast/=5) then ! stochastic estimate
     infofeast=71
     return
  end if
     
     if (com)  write(fout,'(I2,X,2ES25.16,X,I7)') itloop,Emin,Emax,M
  
  
  if (init) then 
     if (M <Mw) then ! not enough eigenvalue found (initial Emin-Emax not well set-up)
        if (fpm40 <0) then
           Emax=Emin+(itloop+1)*diff
        else
           Emin=Emax-(itloop+1)*diff
        endif
     else
        init=.false.
   endif     
end if
   
   if ((M > M0).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax - diff
        else
           Emin = Emin + diff
        endif

     else if ((M < MW).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax + diff
        else
           Emin = Emin - diff
        endif
     endif

  enddo

 
  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Start the exact search (using IFEAST)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  MS=min(M0,Ntotal) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS))

  call feastinit(fpm)
  
  fpm(9)=fpm9
  fpm(49)=fpm49
  
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
  call pdfeast_scsrgv(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,matBj%sa,matBj%isa,matBj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

  if (infofeast/=0) then
     info=72
     return
  end if
     
  
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,M,' (non-stochastic)'

  if (fpm40<0) then
  Emax = (dE(MW)+dE(MW+1))/2.0d0
else
   if (M>MW) then
      Emin=(dE(M-MW)+dE(M-MW+1))/2.0d0
   else
      Emin=(dE(M0-(MW-M))+dE(M0-(MW-M)+1))/2.0d0 
   endif
end if
 
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,MW,' (adjusted)'

  
if ((com).and.(fpm1<0)) close(fout)


 deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)


end subroutine pdfeast_scsrgv_search





subroutine pzfeast_hcsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,itloop,fpm1,fpm9,fpm40,fpm49,Emin,Emax,M0,M,info)
   !  Purpose 
  !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Generalized eigenvalue problem version: AX=BXE  where A is complex Hermitian 
  ! 
  !  DOUBLE PRECISION MPI version
  !
  !  Remark: The code starts by implementing few step of Arnoldi. Few options are hardcoded.
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !
  !  N          (input)        INTEGER: Size system
  !  sa,sb      (input)        COMPLEX DOUBLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
  !  isa,isb    (input)        INTEGER(N+1): CSR row array of Matrix A
  !  jsa,jsb    (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
  !  itloop     (output)       INTEGER- number of stochastic search needed 
  !  fpm1       (input)        INTEGER- Feast parameter fpm(1)
  !  fpm9       (input)        INTEGER- Feast parameter fpm(9)
  !  fpm40      (input)        INTEGER- Feast parameter fpm(40)
  !                                     if fpm(40)=-1 the M0/2 lowest eigenvalues are returned
  !                                     if fpm(40)=1  the M0/2 largest eigenvalues are returned
  !  fpm49       (input)       INTEGER- Feast parameter fpm(49)
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ==================================================================== 
  implicit none
    !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40,fpm49
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: L2_COMM_WORLD,L3_COMM_WORLD,rank2,rank3,code,nb_procs2,nb_procs3,Ntotal
  character(len=3) :: ctemp
  integer(8) :: fout
  !-----
   type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj,matBj
   type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB 
 integer:: Nlocal,startj,endj,distype,nnzA,nnzB,infofeast
!-----

   
info=0 

  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1    

  !----------------------------------------------
  L2_COMM_WORLD=fpm9
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm49
  call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
  call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  !----------------------

  !! local or global?
 call pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 if (distype==0) then ! global- matrix globally known
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    Ntotal=N
 !!! DISTRIBUTE the  sparse matrix by rows among all procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,L3_COMM_WORLD) !! A
      call zcsr_distribute_row(Ntotal,sb,isb,jsb,matBj,startj,endj,Nlocal,Nsize,L3_COMM_WORLD) !! B

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
else ! local! matrix locally known by row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   Nlocal=N
  ! get Ntotal
   Ntotal=0
   call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
 !! obtain global Nsize
     allocate(Nsize(nb_procs3))

     Nsize(1:nb_procs3)=0
     Nsize(rank3+1)=N
     call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)

 !! for A (simple copy for consistency with the rest of the code)
        nnzA=isa(n+1)-1
        allocate(matAj%isa(N+1))
        allocate(matAj%jsa(nnzA))
        allocate(matAj%sa(nnzA))
        matAj%sa(1:nnzA)=sa(1:nnzA)
        matAj%isa(1:n+1)= isa(1:n+1)
        matAj%jsa(1:nnzA)= jsa(1:nnzA)
        matAj%n=N
        matAj%nnz=nnzA
 !! for B
        nnzB=isb(n+1)-1
        allocate(matBj%isa(N+1))
        allocate(matBj%jsa(nnzB))
        allocate(matBj%sa(nnzB))
        matBj%sa(1:nnzB)=sb(1:nnzB)
        matBj%isa(1:n+1)= isb(1:n+1)
        matBj%jsa(1:nnzB)= jsb(1:nnzB)
        matBj%n=N
        matBj%nnz=nnzB
        
end if

!! DISTRIBUTE sparse row into block csr format (needed for matvec with uplo=l,u,f in Arnoldi)
allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
 call zget_local_block_csr(Nsize,mapA,matAj,matAjb,L3_COMM_WORLD,nb_procs3)
 allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapB,matBj,matBjb,L3_COMM_WORLD,nb_procs3)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call dcheck_feast_srci_input(0.0d0,1.0d0,M0,Ntotal,info)
 if (info/=0) return
 
  ! comments?
  com=.false. ! default
  if ((fpm1/=0).and.(rank2==0).and.(rank3==0)) com=.true. !only rank0 is commenting

  ! file or screen output?
  if (com) then
   if (fpm1<0) then
           write(ctemp,'(I3)') abs(fpm1)
           fout=abs(fpm1)+200 ! file number default
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append')
           write(fout,*)  
        elseif (fpm1>0) then
           fout=6 !screen
        end if
     end if
  ! print header   
    if (com ) then   
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') '***** FEAST Search Estimate for Emin-Emax *****'
     write(fout,'(A)') '*** looking for the M0/2 extremal eigenvalue **' 
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') ' Routine pzfeast_hcsrgv_search'
     write(fout,'(A)',advance='no') ' Subspace M0='
     write(fout,'(I4)') M0
     write(fout,'(A)',advance='no') ' Looking for '
     write(fout,'(I4)',advance='no') max(M0/2,1)
     if (fpm40<0) then
        write(fout,'(A)') ' lowest eigenvalues '
     else
        write(fout,'(A)') ' largest eigenvalues '
     end if
     write(fout,*)
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Search the edge of the spectrum using few steps of Arnoldi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  eps=1d-2    
  itmax=min(50,Ntotal) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge)
  allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(Nlocal,itmax)) 
  
!!!!! Search the edge of the spectrum using few steps of Arnoldi
  it=itmax

call pzhgarnoldi(UPLO,Nsize,mapA,matAjb,mapB,matBjb,dE,dX,dres,eps,it,1d-3,30,L3_COMM_WORLD,nb_procs3) ! 2 constant arguments for bicgstab for B matrix

!---
if (com) then
   write(fout,'(A)',advance='no') ' Arnoldi results (Extremal eig. estimated after '
   write(fout,'(I4)',advance='no') it
 write(fout,'(A)') ' steps):'
   Do i=1,it
      if ((i==1).or.(i==it)) then !! add comments for printing them all 
      write(fout,'(I4)',advance='no') i
      write(fout,'(ES25.16)',advance='no') dE(i)
      if (i==1) then
         write(fout,'(A)',advance='no') ' (lowest)  res='
         write(fout,'(ES25.16)',advance='no') dres(i)
   elseif (i==it) then
      write(fout,'(A)',advance='no') ' (largest) res='
      write(fout,'(ES25.16)',advance='no') dres(i)
   end if
   write(fout,*)
end if !!
end Do
end if 


!info=-1
!return


if ((it==itmax).and.(dres(1)>eps).and.(dres(it)>eps)) then
!if (it==itmax) then
   info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
          !! search has failed (ill-conditionned matrix?)
   return
end if

itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
Emin=dE(min(3,Ntotal))-(dE(min(3,Ntotal))-dE(1))*1.05d0
Emax=dE(itmax-min(Ntotal-1,2))+(dE(itmax)-dE(itmax-min(Ntotal-1,2)))*1.05d0

if (com) then
write(fout,*)
 write(fout,'(A)') ' Start FEAST stochastic search:'
 write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
 i=0
 write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,Ntotal
end if


!!!!!!! Find a more appropriate Emax or Emin to start with
  if (fpm40 < 0) then ! search lowest- change Emax
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=dE(2)-Emin
     else
        diff=dE(min(4,n-1))-Emin !dE(1) has converged faster - eig more compact close to E1
     end if
     Emax = Emin+diff

  else ! search largest - change Emin
     if ((abs(dE(itmax))>(abs(dE(1))))) then ! Eitmax largest amplitude (cv faster)- eig more compact close to Eitmax 
        diff=Emax-dE(itmax-min(3,n-1))
     else
        diff=Emax-dE(itmax-1) !dE(1) has converged faster - eig more compact close to E1
     end if

     Emin = Emax-diff
  endif


  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!! Start the stochastic search process
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  MW = max(M0/2,1) ! number of wanted eigenvalue
  call feastinit(fpm)
  
 MS=min(10,Ntotal) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS)) ! local 

  
  fpm(9)=fpm9
  fpm(49)=fpm49

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1
  
  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1
    
     call pzfeast_hcsrgv(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,matBj%sa,matBj%isa,matBj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)
   
       if (infofeast/=5) then ! stochastic estimate
     infofeast=71
     return
  end if
     
     if (com)  write(fout,'(I2,X,2ES25.16,X,I7)') itloop,Emin,Emax,M
  
  
  if (init) then 
     if (M <Mw) then ! not enough eigenvalue found (initial Emin-Emax not well set-up)
        if (fpm40 <0) then
           Emax=Emin+(itloop+1)*diff
        else
           Emin=Emax-(itloop+1)*diff
        endif
     else
        init=.false.
   endif     
end if
   
   if ((M > M0).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax - diff
        else
           Emin = Emin + diff
        endif

     else if ((M < MW).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40 <0) then
           Emax = Emax + diff
        else
           Emin = Emin - diff
        endif
     endif

  enddo

 
  deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Start the exact search (using IFEAST)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  MS=min(M0,Ntotal) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(Nlocal,MS))

  call feastinit(fpm)
  
  fpm(9)=fpm9
  fpm(49)=fpm49
  
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
  call pzfeast_hcsrgv(UPLO,Nlocal,matAj%sa,matAj%isa,matAj%jsa,matBj%sa,matBj%isa,matBj%jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

  if (infofeast/=0) then
     info=72
     return
  end if
     
  
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,M,' (non-stochastic)'

  if (fpm40<0) then
  Emax = (dE(MW)+dE(MW+1))/2.0d0
else
   if (M>MW) then
      Emin=(dE(M-MW)+dE(M-MW+1))/2.0d0
   else
      Emin=(dE(M0-(MW-M))+dE(M0-(MW-M)+1))/2.0d0 
   endif
end if
 
  if (com)  write(fout,'(A,X,2ES25.16,X,I7,A)') '  ',Emin,Emax,MW,' (adjusted)'

  
if ((com).and.(fpm1<0)) close(fout)


 deallocate(dE)     
  deallocate(dres)   
  deallocate(dX)


end subroutine pzfeast_hcsrgv_search


