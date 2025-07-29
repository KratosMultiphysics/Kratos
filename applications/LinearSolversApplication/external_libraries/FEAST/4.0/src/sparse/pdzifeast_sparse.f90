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


!111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE) - using IFEAST/bicgstab
! pdifeast_scsrgvx
! pzifeast_hcsrgvx
! pdifeast_gcsrgvx (wrapper only)
! pzifeast_gcsrgvx
! pzifeast_scsrgvx

!222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
! pdifeast_scsrgv
! pzifeast_hcsrgv
! pdifeast_gcsrgv
! pzifeast_gcsrgv
! pzifeast_scsrgv

!333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
! pdifeast_scsrevx
! pzifeast_hcsrevx
! pdifeast_gcsrevx
! pzifeast_gcsrevx
! pzifeast_scsrevx

!44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
! pdifeast_scsrev
! pzifeast_hcsrev
! pdifeast_gcsrev
! pzifeast_gcsrev
! pzifeast_scsrev



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pdifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  double precision, dimension(:),allocatable ::ddiag
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD

!!!!! matrix format
  integer :: opt,nnza,nnzb,nnz
  complex,dimension(:,:),allocatable :: caux,cwork
 double precision :: DZNRM2
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout
  
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB,mapZ,mapC
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: matAj,matBj
  type(dcsr), dimension(:,:),allocatable :: matAjb,matBjb
  !! single copy
 type scsr
     integer ::n,m,nnz
     real,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type scsr
  type(scsr), dimension(:,:),allocatable :: smatAjb,smatBjb
  !! Az matrix double
  type zcsr
     integer ::n,m,nnz
     complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
 type(zcsr) :: matZj
  type(zcsr), dimension(:,:,:),allocatable :: matZjb
!! Az matrix single
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
 type(ccsr) :: matCj
  type(ccsr), dimension(:,:,:),allocatable :: matCjb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=222145  ! code name
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
     CALL XERBLA( 'PDIFEAST_SCSRGVX', -INFO+100 )
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
!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call dcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call dcsr_distribute_row(Ntotal,sb,isb,jsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B

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

!!! COPY the RHS (for name consistency)
     allocate(Xj(Nlocal,M0))
     if (fpm(5)==1) Xj(1:Nlocal,1:M0)=X(1:Nlocal,1:M0) !!<< blas?

!!! FORM the row distributed matrix (simple copy)
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

  allocate(matZj%isa(Nlocal+1))
  allocate(matZj%sa(1)) ! dummy
  allocate(matZj%jsa(1))! dummy
  opt=1 !! get isaz
  call zdaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
  nnz=matZj%isa(Nlocal+1)-1
  deallocate(matZj%jsa)
  deallocate(matZj%sa)
  allocate(matZj%jsa(nnz))
   allocate(matZj%sa(nnz))
  opt=2 !! get jsaz
  call zdaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
  matZj%sa(1:nnz)= ZZERO
  matZj%n=Nlocal
  matZj%nnz=nnz
 
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zdaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)

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
     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo

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
 
  allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
  call dget_local_block_csr(Nsize,mapA,matAj,matAjb,fpm(49),nb_procs3)
  allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call dget_local_block_csr(Nsize,mapB,matBj,matBjb,fpm(49),nb_procs3)

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision A and B
 allocate(smatAjb(nb_procs3,nb_procs3))
     call dlbcsr2s(Nsize,mapA,matAjb,smatAjb,fpm(49),nb_procs3)
 allocate(smatBjb(nb_procs3,nb_procs3))
     call dlbcsr2s(Nsize,mapB,matBjb,smatBjb,fpm(49),nb_procs3)
     ! get Az matrix local block csr with single precision (first get new single prec. matrix by rows)
        allocate(matCj%isa(Nlocal+1))
        allocate(matCj%jsa(nnz))
        allocate(matCj%sa(nnz))
        matCj%sa(1:nnz)=CZERO
        matCj%isa(1:nlocal+1)= matZj%isa(1:nlocal+1)
        matCj%jsa(1:nnz)= matZj%jsa(1:nnz)
        matCj%n=Nlocal
        matCj%nnz=nnz
  
 allocate(mapC(nb_procs3,nb_procs3))
 allocate(matCjb(nb_procs3,nb_procs3,nfact))

 if (fpm(10)==1) then ! store factorization
    id=0
    do i=rank2+1,fpm(2),nb_procs2   !! multiple fact possible
       id=id+1
 call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,id),fpm(49),nb_procs3)
end do
else ! only 1 factorization 
    call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,1),fpm(49),nb_procs3)
end if

  else ! double precision

  ! get Az matrix local block csr 
 allocate(mapZ(nb_procs3,nb_procs3))
 allocate(matZjb(nb_procs3,nb_procs3,nfact))

 if (fpm(10)==1) then ! store factorization
    id=0
   do i=rank2+1,fpm(2),nb_procs2   !! multiple fact possible
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
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
     call dfeast_srcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,Xj,mode,res,info,Zne,Wne)

     
     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        id=fpm(33) !! id of factorization (for fpm(10) flag)
        !!! form the matrix (version compatible with fpm(10) flag)
       
        if (fpm(42)==1) then !single precision fact- csaz
  call csaddlbcsr(Nsize,-CONE,mapA,smatAjb,cmplx(Ze),mapB,smatBjb,mapC,matCjb(1,1,id),fpm(49),nb_procs3)

else ! double precision fact- saz
            call zdaddlbcsr(Nsize,-ZONE,mapA,matAjb,Ze,mapB,matBjb,mapZ,matZjb(1,1,id),fpm(49),nb_procs3)         
           
end if
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
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
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        call dlbcsrmm(UPLO,'N',Nsize,mapA,fpm(25),DONE,matAjb,Xj(1,fpm(24)),DZERO,work(1,fpm(24)),fpm(49),nb_procs3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call dlbcsrmm(UPLO,'N',Nsize,mapB,fpm(25),DONE,matBjb,Xj(1,fpm(24)),DZERO,work(1,fpm(24)),fpm(49),nb_procs3)

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
      deallocate(smatBjb)
       deallocate(norm)
  else
     deallocate(zaux)
      deallocate(mapZ)
      deallocate(matZjb)
  endif

end subroutine pdifeast_scsrgvx










subroutine pzifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  double precision, dimension(:),allocatable ::ddiag
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD
  
!!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  integer,dimension(:), allocatable :: fisa,fjsa,fisb,fjsb
  integer :: opt,nnza,nnzb,nnz
  complex,dimension(:,:),allocatable :: caux,cwork
  
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout
   double precision :: DZNRM2

!!!!!for cluster pardiso
  !integer(8),dimension(:,:),allocatable :: pt 
  !integer,dimension(:,:),allocatable :: iparm
  !integer :: mtype
  !integer :: MAXFCT,MNUM,PHASE,MSGLVL
  !integer :: idum
  !logical :: mkl
!!!!!!!!!!!!!!!!!!! for mixed precision
!  complex,dimension(:,:),allocatable :: csaz
!  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
 
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB,mapZ,mapC
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj,matBj,matZj
  type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb
  type(zcsr), dimension(:,:,:),allocatable :: matZjb
!! Az matrix single + copy
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
 type(ccsr) :: matCj
  type(ccsr), dimension(:,:,:),allocatable :: matCjb
 type(ccsr), dimension(:,:),allocatable :: smatAjb,smatBjb
 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=242245  ! code name
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
     CALL XERBLA( 'PZIFEAST_HCSRGVX', -INFO+100 )
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
!!! FORMAT CONVERSION TO FULL-CSR!!!! Remark: everything should be full csr since we work with zS-H which is actually unsymmetric 

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

     deallocate(fsa,fisa,fjsa)
     deallocate(fsb,fisb,fjsb)
     
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
     allocate(Xj(Nlocal,M0))
     if (fpm(5)==1) Xj(1:Nlocal,1:M0)=X(1:Nlocal,1:M0) !!<< blas?


!!! FORM the row distributed full matrix 
     ! first simple copy
     ! the case F is then taking care of
     ! the cases L and U use this temp copy and will be modified to
     ! full format once the local blocks csr are defined below

     !!for A
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

  allocate(matZj%isa(Nlocal+1))
  allocate(matZj%sa(1)) ! dummy
  allocate(matZj%jsa(1))! dummy
  opt=1 !! get isaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
  nnz=matZj%isa(Nlocal+1)-1
  deallocate(matZj%jsa)
  deallocate(matZj%sa)
  allocate(matZj%jsa(nnz))
   allocate(matZj%sa(nnz))
  opt=2 !! get jsaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
 matZj%sa(1:nnz)= ZZERO
  matZj%n=Nlocal
  matZj%nnz=nnz
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
  opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)

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
     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo
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
        opt=1 !! get isaz
        call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
  nnz=matZj%isa(Nlocal+1)-1
  deallocate(matZj%jsa)
  deallocate(matZj%sa)
  allocate(matZj%jsa(nnz))
  allocate(matZj%sa(nnz))
  opt=2 !! get jsaz
 call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
 matZj%sa(1:nnz)= ZZERO
  matZj%n=Nlocal
  matZj%nnz=nnz
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then
     
     ! copy matrix single precision A and B
 allocate(smatAjb(nb_procs3,nb_procs3))
 call zlbcsr2c(Nsize,mapA,matAjb,smatAjb,fpm(49),nb_procs3)
 allocate(smatBjb(nb_procs3,nb_procs3))
 call zlbcsr2c(Nsize,mapB,matBjb,smatBjb,fpm(49),nb_procs3)
 
     ! get Az matrix local block csr with single precision (first get new single prec. matrix by rows)
        allocate(matCj%isa(Nlocal+1))
        allocate(matCj%jsa(nnz))
        allocate(matCj%sa(nnz))
        matCj%sa(1:nnz)=CZERO
        matCj%isa(1:nlocal+1)= matZj%isa(1:nlocal+1)
        matCj%jsa(1:nnz)= matZj%jsa(1:nnz)
        matCj%n=Nlocal
        matCj%nnz=nnz
 
 allocate(mapC(nb_procs3,nb_procs3))
 allocate(matCjb(nb_procs3,nb_procs3,nfact))

 if (fpm(10)==1) then ! store factorization
    id=0
    do i=rank2+1,fpm(2),nb_procs2   !! multiple fact possible
       id=id+1
 call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,id),fpm(49),nb_procs3)
end do
else ! only 1 factorization 
    call cget_local_block_csr(Nsize,mapC,matCj,matCjb(1,1,1),fpm(49),nb_procs3)
end if

  else ! double precision
     
  ! get Az matrix local block csr 
 allocate(mapZ(nb_procs3,nb_procs3))
 allocate(matZjb(nb_procs3,nb_procs3,nfact))

 if (fpm(10)==1) then ! store factorization
    id=0
   do i=rank2+1,fpm(2),nb_procs2   !! multiple fact possible
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
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
     call zfeast_hrcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,Xj,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        id=fpm(33) !! id of factorization (for fpm(10) flag)

  if (fpm(42)==1) then !single precision fact- csaz
  call caddlbcsr(Nsize,-CONE,mapA,smatAjb,cmplx(Ze),mapB,smatBjb,mapC,matCjb(1,1,id),fpm(49),nb_procs3)

else ! double precision fact- saz
            call zaddlbcsr(Nsize,-ZONE,mapA,matAjb,Ze,mapB,matBjb,mapZ,matZjb(1,1,id),fpm(49),nb_procs3)         
           
end if
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system (ZeB-A)^H x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      deallocate(smatBjb)
       deallocate(norm)
  else
     deallocate(zaux)
      deallocate(mapZ)
      deallocate(matZjb)
  endif


end subroutine pzifeast_hcsrgvx











subroutine pdifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !    Remark: simple Wrapper to  pzifeast_gcsrpevx
  !
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
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: zsa,zsb
  integer :: nnza,nnzb
!!!!!!!!!!!!!!!!!!!!!!


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=222345   ! code name
     call feastdefault(fpm,info)
  endif
  if (info/=0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! local or global copy 
  nnza=isa(n+1)-1
  allocate(zsa(nnza))
  call ZLACP2( 'F', nnza, 1,sa , nnza, zsa, nnza )

 nnzb=isb(n+1)-1
 allocate(zsb(nnzb))
 call ZLACP2( 'F', nnzb, 1,sb , nnzb, zsb, nnzb )
 
  !! implement wrapping
  call pzifeast_gcsrgvx(N,zsa,isa,jsa,zsb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(zsa)
  deallocate(zsb)
  
end subroutine pdifeast_gcsrgvx





subroutine pzifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !Eric Polizzi 2019
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
!  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  double precision, dimension(:),allocatable ::ddiag
 ! integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
!!!! matrix format
  !complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  integer :: opt,nnza,nnzb,nnz
!!!!!!!!!!!!!!!!!!! for mixed precision
!  complex,dimension(:,:),allocatable :: csaz
!  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
  complex,dimension(:,:),allocatable :: caux,cwork

!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout
   double precision :: DZNRM2
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  character(len=1) :: UPLO2
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB,mapZ,mapC
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj,matBj,matZj
  type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb
 type(zcsr), dimension(:,:,:),allocatable :: matZjb
!! Az matrix single + copy
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
 type(ccsr) :: matCj
  type(ccsr), dimension(:,:,:),allocatable :: matCjb
 type(ccsr), dimension(:,:),allocatable :: smatAjb,smatBjb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=242345   ! code name
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
     CALL XERBLA( 'PZIFEAST_GCSRGVX', -INFO+100 )
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
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!! Distribution schemes for Global matrix known in all rank !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (distype==0) then !! matrix A, and B globally known
     Ntotal=N
!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call zcsr_distribute_row(Ntotal,sb,isb,jsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B

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
  allocate(matZj%isa(Nlocal+1))
  allocate(matZj%sa(1)) ! dummy
  allocate(matZj%jsa(1))! dummy
  opt=1 !! get isaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
  nnz=matZj%isa(Nlocal+1)-1
  deallocate(matZj%jsa)
  deallocate(matZj%sa)
  allocate(matZj%jsa(nnz))
   allocate(matZj%sa(nnz))
  opt=2 !! get jsaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
 matZj%sa(1:nnz)= ZZERO
  matZj%n=Nlocal
  matZj%nnz=nnz
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (fpm(41)==1) then ! if scaling on
  opt=3 !! get saz with given shift to build scaling
     Ze=Emid
     call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)

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
     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo
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
 ! copy matrix single precision A and B
 allocate(smatAjb(nb_procs3,nb_procs3))
 call zlbcsr2c(Nsize,mapA,matAjb,smatAjb,fpm(49),nb_procs3)
 allocate(smatBjb(nb_procs3,nb_procs3))
 call zlbcsr2c(Nsize,mapB,matBjb,smatBjb,fpm(49),nb_procs3)
  ! get Az matrix local block csr with single precision (first get new single prec. matrix by rows)
        allocate(matCj%isa(Nlocal+1))
        allocate(matCj%jsa(nnz))
        allocate(matCj%sa(nnz))
        matCj%sa(1:nnz)=CZERO
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
  allocate(work(Nlocal,M02))
  allocate(zwork(Nlocal,M0))



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
     call dfeast_grcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33) !! id of factorization (for fpm(10) flag)
       
  if (fpm(42)==1) then !single precision fact- csaz
  call caddlbcsr(Nsize,-CONE,mapA,smatAjb,cmplx(Ze),mapB,smatBjb,mapC,matCjb(1,1,id),fpm(49),nb_procs3)

else ! double precision fact- saz
            call zaddlbcsr(Nsize,-ZONE,mapA,matAjb,Ze,mapB,matBjb,mapZ,matZjb(1,1,id),fpm(49),nb_procs3)         
           
end if
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
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
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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
      deallocate(smatBjb)
       deallocate(norm)
  else
     deallocate(zaux)
       deallocate(mapZ)
      deallocate(matZjb)
  end if


end subroutine pzifeast_gcsrgvx




subroutine pzifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B COMPLEX SYMMETRIC :: SPARSE FORMAT 
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
  double precision, dimension(:),allocatable ::ddiag
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  !integer :: rank,code,nb_procs,NEW_COMM_WORLD

!!!!! matrix format
  integer :: opt,nnza,nnzb,nnz
  complex,dimension(:,:),allocatable :: caux,cwork
 
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout
  double precision :: DZNRM2
!!!!!!!!!!!!!!!! pfeast
  integer :: rank2,rank3,nb_procs2,nb_procs3,code
  integer :: L2_COMM_WORLD,L3_COMM_WORLD
  integer :: p,Ntotal,Nlocal,startj,endj,distype
  integer,dimension(:),allocatable :: Nsize
  integer,dimension(:,:),allocatable :: mapA,mapB,mapZ,mapC
  !! Az matrix double
  type zcsr
     integer ::n,m,nnz
     complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
 type(zcsr) ::matAj,matBj, matZj
 type(zcsr), dimension(:,:,:),allocatable :: matZjb
  type(zcsr), dimension(:,:),allocatable :: matAjb,matBjb
!! Az matrix single
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
 type(ccsr) :: matCj
  type(ccsr), dimension(:,:,:),allocatable :: matCjb
  type(ccsr), dimension(:,:),allocatable :: smatAjb,smatBjb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=242145  ! code name
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
     CALL XERBLA( 'PZIFEAST_SCSRGVX', -INFO+100 )
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
!!! DISTRIBUTE the  sparse matrix by rows among all L3 procs 
     allocate(Nsize(nb_procs3))
     call zcsr_distribute_row(Ntotal,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,fpm(49)) !! A
     call zcsr_distribute_row(Ntotal,sb,isb,jsb,matBj,startj,endj,Nlocal,Nsize,fpm(49)) !! B

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

!!! COPY the RHS (for name consistency)
     allocate(Xj(Nlocal,M0))
     if (fpm(5)==1) Xj(1:Nlocal,1:M0)=X(1:Nlocal,1:M0) !!<< blas?

!!! FORM the row distributed matrix (simple copy)
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

  allocate(matZj%isa(Nlocal+1))
  allocate(matZj%sa(1)) ! dummy
  allocate(matZj%jsa(1))! dummy
  opt=1 !! get isaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
  nnz=matZj%isa(Nlocal+1)-1
  deallocate(matZj%jsa)
  deallocate(matZj%sa)
  allocate(matZj%jsa(nnz))
   allocate(matZj%sa(nnz))
  opt=2 !! get jsaz
  call zaddcsr(Nlocal,Ntotal,opt,ZONE,matAj%sa,matAj%isa,matAj%jsa,ZONE,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)
 matZj%sa(1:nnz)= ZZERO
  matZj%n=Nlocal
  matZj%nnz=nnz
 
!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid
     call zaddcsr(Nlocal,Ntotal,opt,-ZONE,matAj%sa,matAj%isa,matAj%jsa,Ze,matBj%sa,matBj%isa,matBj%jsa,matZj%sa,matZj%isa,matZj%jsa)

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
     !scale matrix A and B
     do i=1,Nlocal
        do k=matAj%isa(i),matAj%isa(i+1)-1   
           matAj%sa(k)=matAj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matAj%jsa(k))))
        enddo
        do k=matBj%isa(i),matBj%isa(i+1)-1
           matBj%sa(k)=matBj%sa(k)/(sqrt(ddiag(startj+i-1))*sqrt(ddiag(matBj%jsa(k))))
        enddo

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
 
  allocate(mapA(nb_procs3,nb_procs3))
  allocate(matAjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapA,matAj,matAjb,fpm(49),nb_procs3)
  allocate(mapB(nb_procs3,nb_procs3))
  allocate(matBjb(nb_procs3,nb_procs3))
  call zget_local_block_csr(Nsize,mapB,matBj,matBjb,fpm(49),nb_procs3)

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then 
     ! copy matrix single precision A and B
 allocate(smatAjb(nb_procs3,nb_procs3))
     call zlbcsr2c(Nsize,mapA,matAjb,smatAjb,fpm(49),nb_procs3)
 allocate(smatBjb(nb_procs3,nb_procs3))
     call zlbcsr2c(Nsize,mapB,matBjb,smatBjb,fpm(49),nb_procs3)
     ! get Az matrix local block csr with single precision (first get new single prec. matrix by rows)
        allocate(matCj%isa(Nlocal+1))
        allocate(matCj%jsa(nnz))
        allocate(matCj%sa(nnz))
        matCj%sa(1:nnz)=CZERO
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
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
     call zfeast_srcix(ijob,Nlocal,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,Xj,mode,res,info,Zne,Wne)

     
     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        id=fpm(33) !! id of factorization (for fpm(10) flag)
        !!! form the matrix (version compatible with fpm(10) flag)
       
        if (fpm(42)==1) then !single precision fact- csaz
  call caddlbcsr(Nsize,-CONE,mapA,smatAjb,cmplx(Ze),mapB,smatBjb,mapC,matCjb(1,1,id),fpm(49),nb_procs3)

else ! double precision fact- saz
            call zaddlbcsr(Nsize,-ZONE,mapA,matAjb,Ze,mapB,matBjb,mapZ,matZjb(1,1,id),fpm(49),nb_procs3)         
           
end if
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
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
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        call zlbcsrmm(UPLO,'N',Nsize,mapA,fpm(25),ZONE,matAjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call zlbcsrmm(UPLO,'N',Nsize,mapB,fpm(25),ZONE,matBjb,Xj(1,fpm(24)),ZZERO,work(1,fpm(24)),fpm(49),nb_procs3)

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
      deallocate(smatBjb)
       deallocate(norm)
  else
     deallocate(zaux)
      deallocate(mapZ)
      deallocate(matZjb)
  endif

end subroutine pzifeast_scsrgvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine pdifeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pdifeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=222144
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
  call pdifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdifeast_scsrgv




subroutine pzifeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pzifeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=242244
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
  call pzifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzifeast_hcsrgv



subroutine pdifeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !Eric Polizzi 2019
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
!!!    Wrapper Routine to expert routine: pdifeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=222344
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pdifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pdifeast_gcsrgv




subroutine pzifeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !Eric Polizzi 2019
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
!!!    Wrapper Routine to expert routine: pzifeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=242344
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzifeast_gcsrgv




subroutine pzifeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !Eric Polizzi 2019
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
!!!    Wrapper Routine to expert routine: pzifeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=242144
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call pzifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine pzifeast_scsrgv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine pdifeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  ! Eric Polizzi 2019
  !=====================================================================
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: pdifeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  fpm(30)=222143
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


  call pdifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pdifeast_scsrevx



subroutine pzifeast_hcsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
!!!    Wrapper Routine to expert routine: pzifeast_hcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code


  fpm(30)=242243
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

  call pzifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pzifeast_hcsrevx



subroutine pdifeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !Eric Polizzi 2019
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
!!!    Wrapper Routine to expert routine: pdifeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code


  fpm(30)=222343
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

  call pdifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)


end subroutine pdifeast_gcsrevx




subroutine pzifeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
!!!    Wrapper Routine to expert routine: pzifeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code

  fpm(30)=242343
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

  call pzifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pzifeast_gcsrevx


subroutine pzifeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
!!!    Wrapper Routine to expert routine: pzifeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code


  fpm(30)=242143
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

  call pzifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine pzifeast_scsrevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine pdifeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! Eric Polizzi 2019
  !=====================================================================
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
!!!    Wrapper Routine to expert routine: pdifeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 


  fpm(30)=222142
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

  call pdifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pdifeast_scsrev






subroutine pzifeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pzifeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 

  fpm(30)=242242
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

  call pzifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pzifeast_hcsrev



subroutine pdifeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  !Eric Polizzi 2019
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
!!!    Wrapper Routine to expert routine: pdifeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

  fpm(30)=222342
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
  call pdifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine pdifeast_gcsrev




subroutine pzifeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pzifeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code    
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

  fpm(30)=242342
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

  call pzifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine pzifeast_gcsrev




subroutine pzifeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
!!!    Wrapper Routine to expert routine: pzifeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb,Nsize
  integer :: i,distype,L3_COMM_WORLD,rank3,nb_procs3,startj,code

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

  fpm(30)=242142
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
  call pzifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine pzifeast_scsrev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



