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

!{D,Z}IFEAST_{SCSR,HCSR,GCSR}{EV,GV}



!111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE) - using IFEAST/Bicgstab
! difeast_scsrgvx
! zifeast_hcsrgvx
! difeast_gcsrgvx
! zifeast_gcsrgvx
! zifeast_scsrgvx

!222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
! difeast_scsrgv
! zifeast_hcsrgv
! difeast_gcsrgv
! zifeast_gcsrgv
! zifeast_scsrgv

!333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
! difeast_scsrevx
! zifeast_hcsrevx
! difeast_gcsrevx
! zifeast_gcsrevx
! zifeast_scsrevx

!44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
! difeast_scsrev
! zifeast_hcsrev
! difeast_gcsrev
! zifeast_gcsrev
! zifeast_scsrev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine difeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  !include 'f90_noruntime_interface.fi'
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
  integer :: ijob,infoloc1,infoloc2,i,s,k
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  double precision, dimension(:,:),allocatable ::work,Aq,Sq
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
  !! matrix format
  double precision, dimension(:),allocatable ::dsa,dsb,ddiag
  integer:: nnz,opt
!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  real,dimension(:),allocatable :: ssa,ssb !single precision copy
  complex, dimension(:,:),allocatable ::cwork,caux
  double precision :: DZNRM2
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout 
!!!!!!!!!!!!!!!!

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=122145 ! code name
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
     CALL XERBLA( 'DIFEAST_SCSRGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! copy the sa,sb (needed for scaling)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(dsa(isa(N+1)-1))
  allocate(dsb(isb(N+1)-1))
  call DCOPY(isa(N+1)-1,sa,1,dsa,1)
  call DCOPY(isb(N+1)-1,sb,1,dsb,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store "factorization"
     nfact=0
     !! nfact is local (number of total "factorization" by contour points)
     do i=rank+1,fpm(2),nb_procs
        nfact=nfact+1
     end do
  end if


!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zdaddcsr(N,N,opt,ZONE,sa,isa,jsa,ZONE,sb,isb,jsb,saz,isaz,jsaz) 
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 ! get jsaz
  call zdaddcsr(N,N,opt,ZONE,sa,isa,jsa,ZONE,sb,isb,jsb,saz,isaz,jsaz) 
  deallocate(saz) ! deallocate dummy       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zdaddcsr(N,N,opt,-ZONE,sa,isa,jsa,Ze,sb,isb,jsb,saz,isaz,jsaz)
     allocate(ddiag(N))
     ! ddiag(1:N)=DZERO
     ddiag(1:N)=DONE
     ! extract diagonal 
     do i=1,N
        do k=isaz(i),isaz(i+1)-1
           !if (jsaz(k)==i) ddiag(i)=abs(saz(k,1))
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix
           if ((jsaz(k)==i).and.(abs(saz(k,1))/=DZERO)) ddiag(i)=abs(saz(k,1)) 
        enddo
     enddo
     !scale matrix A and B
     do i=1,N
        do k=isa(i),isa(i+1)-1   
           dsa(k)=dsa(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsa(k))))
        enddo
        do k=isb(i),isb(i+1)-1
           dsb(k)=dsb(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsb(k))))
        enddo
     end do

     deallocate(saz) ! deallocate dummy

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE RHS if fpm(5)=1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((fpm(41)==1).and.(fpm(5)==1)) then
     do i=1,N
        X(i,1:M0)=X(i,1:M0)*sqrt(ddiag(i)) !!<<
     enddo
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then
     ! copy matrix single precision   
     allocate(ssa(isa(N+1)-1))
     allocate(ssb(isb(N+1)-1))
     call DLAG2S(isa(N+1)-1,1,dsa,isa(N+1)-1,ssa,isa(N+1)-1,infoloc1)
     call DLAG2S(isb(N+1)-1,1,dsb,isb(N+1)-1,ssb,isb(N+1)-1,infoloc1)
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
  allocate(work(N,M0))
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! iterative method  set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  allocate(nres(M0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(N,M0))
     allocate(cwork(N,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(N,M0))
  endif


  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if (fpm(1)/=0) com=.true.
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
     call dfeast_srcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
     case(10) !! form (zeB-A) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
!!! form the matrix (version compatible with fpm(10) flag)
        opt=3
        if (fpm(42)==1) then !single precision fact- csaz
           call csaddcsr(N,N,opt,-CONE,ssa,isa,jsa,cmplx(Ze),ssb,isb,jsb,csaz(1,id),isaz,jsaz) !! get csaz
        else
           call zdaddcsr(N,N,opt,-ZONE,dsa,isa,jsa,Ze,dsb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)


!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !if (fpm(48)==0) then!  (max norm within interval)
           !if (loop>0) then
           !   do i=1,fpm(23)
           !      if (res(i)<100.0d0) nres(i)=0.0d0
           !   enddo
           !endif
        !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!!!

        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)

           call cbicgstab(fpm(44),UPLO,'N',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)

           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo


        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess
!print *,'bef',sum(abs(zwork(1:N,1)))
           call zbicgstab(fpm(44),UPLO,'N',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
!           print *,'aft',sum(abs(zwork(1:N,1)))
!           print *,nres(1),linloops,lintargeterror
!return
           
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        call wdcsrmm(UPLO,'N',N,N,fpm(25),DONE,dsa,isa,jsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wdcsrmm(UPLO,'N',N,N,fpm(25),DONE,dsb,isb,jsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))

     end select
  end do

  !loop=fpm(60)  ! new value for loop


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)
  deallocate(nres)

  deallocate(dsa)
  deallocate(dsb)


  deallocate(isaz)
  deallocate(jsaz)

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(ssa)
     deallocate(ssb)
     deallocate(norm)
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


end subroutine difeast_scsrgvx






subroutine zifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  !include 'f90_noruntime_interface.fi'
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,Aq,Sq,zaux
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
  !! matrix format
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: fsa,fsb
  integer,dimension(:), allocatable :: fisa,fjsa,fisb,fjsb
  double precision, dimension(:),allocatable ::ddiag
  integer:: nnz,opt,nnza,nnzb
!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
  complex, dimension(:,:),allocatable ::cwork,caux
  double precision :: DZNRM2
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable ::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8) :: fout =6
!!!!!!!!!!!!!!!!!!!!!!!!!

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=142245 ! ! code name
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
     CALL XERBLA( 'ZIFEAST_HCSRGVX', -INFO+100 )
     RETURN
  END IF


  infoloc1=0
  infoloc2=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR  !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be full csr since we work with zS-H which is actually unsymmetric 



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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store "factorization"
     nfact=0
     !! nfact is local (number of total "factorization" by contour points)
     do i=rank+1,fpm(2),nb_procs
        nfact=nfact+1
     end do
  end if


!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zaddcsr(N,N,opt,ZONE,fsa,fisa,fjsa,ZONE,fsb,fisb,fjsb,saz,isaz,jsaz) 
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 ! get jsaz
  call zaddcsr(N,N,opt,ZONE,fsa,fisa,fjsa,ZONE,fsb,fisb,fjsb,saz,isaz,jsaz) 
  deallocate(saz) ! deallocate dummy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(N,N,opt,-ZONE,fsa,fisa,fjsa,Ze,fsb,fisb,fjsb,saz,isaz,jsaz)
     allocate(ddiag(N))
     ! ddiag(1:N)=DZERO
     ddiag(1:N)=DONE
     ! extract diagonal 
     do i=1,N
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix
           if ((jsaz(k)==i).and.(abs(saz(k,1))/=DZERO)) ddiag(i)=abs(saz(k,1)) 
        enddo
     enddo
     !scale matrix A and B
     do i=1,N
        do k=fisa(i),fisa(i+1)-1   
           fsa(k)=fsa(k)/(sqrt(ddiag(i))*sqrt(ddiag(fjsa(k))))
        enddo
        do k=fisb(i),fisb(i+1)-1
           fsb(k)=fsb(k)/(sqrt(ddiag(i))*sqrt(ddiag(fjsb(k))))
        enddo
     end do
     deallocate(saz) ! deallocate dummy
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE RHS if fpm(5)=1 !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((fpm(41)==1).and.(fpm(5)==1)) then
     do i=1,N
        X(i,1:M0)=X(i,1:M0)*sqrt(ddiag(i)) !!<<
     enddo
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!! Mixed precision set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(42)==1) then
     ! copy matrix single precision   
     allocate(cfsa(fisa(N+1)-1))
     allocate(cfsb(fisb(N+1)-1))
     call ZLAG2C(fisa(N+1)-1,1,fsa,fisa(N+1)-1,cfsa,fisa(N+1)-1,infoloc1)
     call ZLAG2C(fisb(N+1)-1,1,fsb,fisb(N+1)-1,cfsb,fisb(N+1)-1,infoloc1)
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

  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
  allocate(work(N,M0))
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  iterative method set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(nres(m0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(N,M0))
     allocate(cwork(N,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(N,M0))
  endif


  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if (fpm(1)/=0) com=.true.
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
     call zfeast_hrcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     case(10) !!  form Az=(zeB-A) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
!!! form the matrix (version compatible with fpm(10) flag)
        opt=3

        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(N,N,opt,-CONE,cfsa,fisa,fjsa,cmplx(Ze),cfsb,fisb,fjsb,csaz(1,id),isaz,jsaz) !! get csaz
        else
           call zaddcsr(N,N,opt,-ZONE,fsa,fisa,fjsa,Ze,fsb,fisb,fjsb,saz(1,id),isaz,jsaz) !! get saz
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)          
        !          print *,Ze
        !open(10,name='y',status='replace')
        !             do k=1,fpm(23)
        !                do i=1,N
        ! write(10,*) zwork(i,k)
        !             end do
        !          end do
        !close(10)

        !          comb=.true.
!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !if (fpm(48)==0) then!  (max norm within interval)
        !   if (loop>0) then
        !      do i=1,fpm(23)
        !         if (res(i)<100.0d0) nres(i)=0.0d0
        !      enddo
        !   endif
     !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!

        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)

           call cbicgstab(fpm(44),'F','N',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)

           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess

           call zbicgstab(fpm(44),'F','N',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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

        !stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system (ZeB-A)^H x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)


!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !if (fpm(48)==0) then!  (max norm within interval)
        !   if (loop>0) then
        !      do i=1,fpm(23)
        !         if (res(i)<100.0d0) nres(i)=0.0d0
        !      enddo
        !   endif
        !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!
!!!!!!!!

        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)

           call cbicgstab(fpm(44),'F','C',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)

           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess

           call zbicgstab(fpm(44),'F','C',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzhcsrmm('F','N',N,N,fpm(25),ZONE,fsa,fisa,fjsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        call wzhcsrmm('F','N',N,N,fpm(25),ZONE,fsb,fisb,fjsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     end select
  end do


  !loop=fpm(60) ! new value for loop


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)
  deallocate(nres)

  deallocate(isaz)
  deallocate(jsaz)   

  deallocate(fsa)
  deallocate(fisa)
  deallocate(fjsa)

  deallocate(fsb)
  deallocate(fisb)
  deallocate(fjsb)


  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(cfsa)
     deallocate(cfsb)
     deallocate(norm)
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


end subroutine zifeast_hcsrgvx









subroutine difeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  !include 'f90_noruntime_interface.fi'
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
  integer :: ijob,infoloc1,infoloc2,i,s,k
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,Aq,Sq,zaux
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  ! complex(kind=(kind(1.0d0))),dimension(:),allocatable :: saz0,sbz0
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id   
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
  !! matrix format
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: fsa,fsb
  !  integer,dimension(:), allocatable :: fisa,fjsa,fisb,fjsb
  double precision, dimension(:),allocatable ::ddiag
  integer:: nnz,opt,nnza,nnzb
!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
  complex, dimension(:,:),allocatable ::cwork,caux
  double precision :: DZNRM2
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable ::nres,norm
  integer :: linloops,infob
  double precision :: lintargeterror
  integer(8) :: fout =6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=122345 ! code name
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


  IF( N<=0 ) THEN
     INFO = -102
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DIFEAST_GCSRGVX', -INFO+100 )
     RETURN
  END IF



  infoloc1=0
  infoloc2=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! copy sa,sb (needed for scaling)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Perform a copy of the matrix values (transform into complex)
  nnza=isa(n+1)-1
  allocate(fsa(nnza))
  call ZLACP2( 'F', nnza, 1,sa , nnza, fsa, nnza )

  nnzb=isb(n+1)-1
  allocate(fsb(nnzb))
  call ZLACP2( 'F', nnzb, 1,sb , nnzb, fsb, nnzb ) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store "factorization"
     nfact=0
     !! nfact is local (number of total "factorization" by contour points)
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zdaddcsr(N,N,opt,ZONE,sa,isa,jsa,ZONE,sb,isb,jsb,saz,isaz,jsaz) 
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !get jsaz
  call zdaddcsr(N,N,opt,ZONE,sa,isa,jsa,ZONE,sb,isb,jsb,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,1),isaz,jsaz)
     allocate(ddiag(N))
     !ddiag(1:N)=DZERO
     ddiag(1:N)=DONE
     ! extract diagonal 
     do i=1,N
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix
           if ((jsaz(k)==i).and.(abs(saz(k,1))/=DZERO)) ddiag(i)=abs(saz(k,1)) 
        enddo
     enddo
     !scale matrix A and B
     do i=1,N
        do k=isa(i),isa(i+1)-1   
           fsa(k)=fsa(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsa(k))))
        enddo
        do k=isb(i),isb(i+1)-1
           fsb(k)=fsb(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsb(k))))
        enddo

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
     allocate(cfsa(isa(N+1)-1))
     allocate(cfsb(isb(N+1)-1))
     call ZLAG2C(isa(N+1)-1,1,fsa,isa(N+1)-1,cfsa,isa(N+1)-1,infoloc1)
     call ZLAG2C(isb(N+1)-1,1,fsb,isb(N+1)-1,cfsb,isb(N+1)-1,infoloc1)
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

  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  iterative method set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(nres(m0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(N,M0))
     allocate(cwork(N,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(N,M0))
  endif

  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if (fpm(1)/=0) com=.true.
  if (fpm(1)<0) then
     fout=abs(fpm(1))+200!fpm(60) !!file id name
  else
     fout=6 !screen
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!!!!!!!!!!!!!!!! FEAST-RCI  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization
  do while (ijob/=0) 
     call dfeast_grcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(10) !! form (zeB-A) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
!!! form the matrix (version compatible with fpm(10) flag)
        opt=3
        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(N,N,opt,-CONE,cfsa,isa,jsa,cmplx(Ze),cfsb,isb,jsb,csaz(1,id),isaz,jsaz) !! get csaz
        else
           call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)  
!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !if (fpm(48)==0) then!  (max norm within interval)
        !   if (loop>0) then
        !      do i=1,fpm(23)
        !         if (res(i)<100.0d0) nres(i)=0.0d0
        !      enddo
        !   endif
        !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!        
        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call cbicgstab(fpm(44),'F','N',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)           
           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo


        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess
           call zbicgstab(fpm(44),'F','N',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)


!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !if (fpm(48)==0) then!  (max norm within interval)
        !   if (loop>0) then
        !      do i=1,fpm(23)
        !         if (res(i)<100.0d0) nres(i)=0.0d0
        !      enddo
        !   endif
        !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!
!!!!!!!!

        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)

           call cbicgstab(fpm(44),'F','C',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)

           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess

           call zbicgstab(fpm(44),'F','C',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsa,isa,jsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsa,isa,jsa,X(1,fpm(34)),ZZERO,work(1,fpm(34)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsb,isb,jsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsb,isb,jsb,X(1,fpm(34)),ZZERO,work(1,fpm(34)))


     end select
  end do


  !loop=fpm(60) ! new value for loop


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)
  deallocate(nres)


  deallocate(isaz)
  deallocate(jsaz)

  deallocate(fsa)
  deallocate(fsb)



  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(cfsa)
     deallocate(cfsb)
     deallocate(norm)
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



end subroutine difeast_gcsrgvx






subroutine zifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
  complex(kind=(kind(1.0d0))),dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,Aq,Sq,zaux
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  ! complex(kind=(kind(1.0d0))),dimension(:),allocatable :: saz0,sbz0
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id   
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
  !! matrix format
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: fsa,fsb
  !  integer,dimension(:), allocatable :: fisa,fjsa,fisb,fjsb
  double precision, dimension(:),allocatable ::ddiag
  integer:: nnz,opt,nnza,nnzb
!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
  complex, dimension(:,:),allocatable ::cwork,caux
  double precision :: DZNRM2
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable ::nres,norm
  integer :: linloops,infob
  double precision :: lintargeterror
  integer(8) :: fout =6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=142345 ! code name
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


  IF( N<=0 ) THEN
     INFO = -102
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZIFEAST_GCSRGVX', -INFO+100 )
     RETURN
  END IF



  infoloc1=0
  infoloc2=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! copy sa,sb (needed for scaling)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Perform a copy of the matrix values (transform into complex)
  nnza=isa(n+1)-1
  allocate(fsa(nnza))
  call ZCOPY(nnza,sa,1, fsa, 1 )
  nnzb=isb(n+1)-1
  allocate(fsb(nnzb)) 
  call ZCOPY( nnzb,sb,1, fsb, 1 ) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store "factorization"
     nfact=0
     !! nfact is local (number of total "factorization" by contour points)
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zaddcsr(N,N,opt,ZONE,fsa,isa,jsa,ZONE,fsb,isb,jsb,saz,isaz,jsaz) 
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !get jsaz
  call zaddcsr(N,N,opt,ZONE,fsa,isa,jsa,ZONE,fsb,isb,jsb,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,1),isaz,jsaz)
     allocate(ddiag(N))
     !ddiag(1:N)=DZERO
     ddiag(1:N)=DONE
     ! extract diagonal 
     do i=1,N
        do k=isaz(i),isaz(i+1)-1
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix
           if ((jsaz(k)==i).and.(abs(saz(k,1))/=DZERO)) ddiag(i)=abs(saz(k,1)) 
        enddo
     enddo
     !scale matrix A and B
     do i=1,N
        do k=isa(i),isa(i+1)-1   
           fsa(k)=fsa(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsa(k))))
        enddo
        do k=isb(i),isb(i+1)-1
           fsb(k)=fsb(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsb(k))))
        enddo

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
     allocate(cfsa(isa(N+1)-1))
     allocate(cfsb(isb(N+1)-1))
     call ZLAG2C(isa(N+1)-1,1,fsa,isa(N+1)-1,cfsa,isa(N+1)-1,infoloc1)
     call ZLAG2C(isb(N+1)-1,1,fsb,isb(N+1)-1,cfsb,isb(N+1)-1,infoloc1)
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

  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  iterative method set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(nres(M0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(N,M0))
     allocate(cwork(N,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(N,M0))
  endif

  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if (fpm(1)/=0) com=.true.
  if (fpm(1)<0) then
     fout=abs(fpm(1))+200!fpm(60) !!file id name
  else
     fout=6 !screen
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!!!!!!!!!!!!!!!! FEAST-RCI  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_grcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(10) !! form (zeB-A) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
!!! form the matrix (version compatible with fpm(10) flag)
        opt=3
        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(N,N,opt,-CONE,cfsa,isa,jsa,cmplx(Ze),cfsb,isb,jsb,csaz(1,id),isaz,jsaz) !! get csaz
        else
           call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)  
!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !if (fpm(48)==0) then!  (max norm within interval)
        !   if (loop>0) then
        !      do i=1,fpm(23)
        !         if (res(i)<100.0d0) nres(i)=0.0d0
        !      enddo
        !   endif
        !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!        
        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call cbicgstab(fpm(44),'F','N',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)           
           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo


        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess
           call zbicgstab(fpm(44),'F','N',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)


!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !if (fpm(48)==0) then!  (max norm within interval)
        !   if (loop>0) then
        !      do i=1,fpm(23)
        !         if (res(i)<100.0d0) nres(i)=0.0d0
        !      enddo
        !   endif
     !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!
!!!!!!!!

        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)

           call cbicgstab(fpm(44),'F','C',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)

           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo

        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess

           call zbicgstab(fpm(44),'F','C',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsa,isa,jsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsa,isa,jsa,X(1,fpm(34)),ZZERO,work(1,fpm(34)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsb,isb,jsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsb,isb,jsb,X(1,fpm(34)),ZZERO,work(1,fpm(34)))


     end select
  end do


  !loop=fpm(60) ! new value for loop


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)
  deallocate(nres)


  deallocate(isaz)
  deallocate(jsaz)

  deallocate(fsa)
  deallocate(fsb)



  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(cfsa)
     deallocate(cfsb)
     deallocate(norm)
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

end subroutine zifeast_gcsrgvx





subroutine zifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST SPARSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B COMPLEX  SYMMETRIC:: SPARSE FORMAT 
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
  !include 'f90_noruntime_interface.fi'
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
  complex(kind=kind(1.0d0)),dimension(*)  :: E
  complex(kind=kind(1.0d0)),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,k
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,zaux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,Aq,Sq
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: saz
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
  !! matrix format
  complex(kind=(kind(1.0d0))), dimension(:),allocatable ::zsa,zsb
  double precision,dimension(:),allocatable :: ddiag
  integer:: nnz,opt
!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: csa,csb !single precision copy
  complex, dimension(:,:),allocatable ::cwork,caux
  double precision :: DZNRM2
!!!! for bicgstab
  logical :: comb,com
  double precision, dimension(:),allocatable::nres,norm
  integer :: linloops
  double precision :: lintargeterror
  integer(8)  :: fout 
!!!!!!!!!!!!!!!!

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=142145 ! code name
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
     CALL XERBLA( 'ZIFEAST_SCSRGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! copy the sa,sb (needed for scaling)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(zsa(isa(N+1)-1))
  allocate(zsb(isb(N+1)-1))
  call ZCOPY(isa(N+1)-1,sa,1,zsa,1)
  call ZCOPY(isb(N+1)-1,sb,1,zsb,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store "factorization"
     nfact=0
     !! nfact is local (number of total "factorization" by contour points)
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  end if


!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zaddcsr(N,N,opt,ZONE,zsa,isa,jsa,ZONE,zsb,isb,jsb,saz,isaz,jsaz) 
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 ! get jsaz
  call zaddcsr(N,N,opt,ZONE,zsa,isa,jsa,ZONE,zsb,isb,jsb,saz,isaz,jsaz) 
  deallocate(saz) ! deallocate dummy       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid
     call zaddcsr(N,N,opt,-ZONE,zsa,isa,jsa,Ze,zsb,isb,jsb,saz,isaz,jsaz)
     allocate(ddiag(N))
     ! ddiag(1:N)=DZERO
     ddiag(1:N)=DONE
     ! extract diagonal 
     do i=1,N
        do k=isaz(i),isaz(i+1)-1
           !if (jsaz(k)==i) ddiag(i)=abs(saz(k,1))
           ! include the corner case where Emid is the eigenvalue of a diagonal matrix
           if ((jsaz(k)==i).and.(abs(saz(k,1))/=DZERO)) ddiag(i)=abs(saz(k,1)) 
        enddo
     enddo
     !scale matrix A and B
     do i=1,N
        do k=isa(i),isa(i+1)-1   
           zsa(k)=zsa(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsa(k))))
        enddo
        do k=isb(i),isb(i+1)-1
           zsb(k)=zsb(k)/(sqrt(ddiag(i))*sqrt(ddiag(jsb(k))))
        enddo
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
     allocate(csa(isa(N+1)-1))
     allocate(csb(isb(N+1)-1))
     call ZLAG2C(isa(N+1)-1,1,zsa,isa(N+1)-1,csa,isa(N+1)-1,infoloc1)
     call ZLAG2C(isb(N+1)-1,1,zsb,isb(N+1)-1,csb,isb(N+1)-1,infoloc1)
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


  allocate(Aq(M0,M0))
  allocate(Sq(M0,M0))
  allocate(work(N,M0))
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! iterative method  set-up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  allocate(nres(M0))

  if (fpm(42)==1) then ! mixed precision
     allocate(caux(N,M0))
     allocate(cwork(N,M0))
     allocate(norm(M0))
  else ! full precision
     allocate(zaux(N,M0))
  endif


  comb=.false.
  linloops=10000
  !fpm(60)=0 !!! total number of iterations
  ! comments?
  com=.false. !default
  if (fpm(1)/=0) com=.true.
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
     call zfeast_srcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
     case(10) !! form (zeB-A) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
!!! form the matrix (version compatible with fpm(10) flag)
        opt=3
        if (fpm(42)==1) then !single precision fact- csaz
           call caddcsr(N,N,opt,-CONE,csa,isa,jsa,cmplx(Ze),csb,isb,jsb,csaz(1,id),isaz,jsaz) !! get csaz
        else
           call zaddcsr(N,N,opt,-ZONE,zsa,isa,jsa,Ze,zsb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)


!!!!!!! bicgstab preparation
        lintargeterror=10d0**(-fpm(45))
        if (fpm(47)==1) then
           linloops=min(fpm(46),linloops) ! retain the minimum of linloops
        else
           linloops=fpm(46)
        end if

        nres=1.0d0 ! max rhs norm
        !!if (fpm(48)==0) then!  (max norm within interval)
        !   if (loop>0) then
        !      do i=1,fpm(23)
        !         if (res(i)<100.0d0) nres(i)=0.0d0
        !      enddo
        !   endif
        !else
        if (fpm(48)==-1) then ! min norm
           nres=-1.0d0
        end if
!!!!!!!!!!

        if (fpm(42)==1) then ! mixed precision
           caux(1:N,1:fpm(23))=CZERO !! initial guess
           !! Rq: the rhs needs to be scaled. Indeed, it will progressively goes to 0 (with values as small as 1e-20) along the FEAST iterations. Bicgstab requires an inner product calculation that will could lead to values <1e-38 (which will not work in single precision)
           do i=1,fpm(23)
              norm(i)=DZNRM2(N,zwork(1,i),1)
              call ZSCAL(N,ZONE/norm(i),zwork(1,i),1)
           enddo
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)

           call cbicgstab(fpm(44),UPLO,'N',csaz(1,id),isaz,jsaz,N,fpm(23),cwork,caux,nres,linloops,lintargeterror,comb,infoloc2)

           call CLAG2Z(N, fpm(23), caux, N, zwork, N, infoloc1)
           do i=1,fpm(23)
              call ZSCAL(N,ZONE*norm(i),zwork(1,i),1)
           enddo


        else ! full precision

           zaux(1:N,1:fpm(23))=ZZERO !! initial guess

           call zbicgstab(fpm(44),UPLO,'N',saz(1,id),isaz,jsaz,N,fpm(23),zwork,zaux,nres,linloops,lintargeterror,comb,infoloc2)

           call ZCOPY( N*fpm(23), zaux, 1, zwork, 1 )
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
           if (nb_procs>1) then
              if (.not.((rank>0).and.(fout/=6))) then !all mpi prints on screen but not on file (flush problem)
                 write(fout,'(A)',advance='no') ' rank'
                 write(fout,'(I4)',advance='no') rank
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        call wzcsrmm(UPLO,'N',N,N,fpm(25),ZONE,zsa,isa,jsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm(UPLO,'N',N,N,fpm(25),ZONE,zsb,isb,jsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

     end select
  end do

  !loop=fpm(60)  ! new value for loop


  deallocate(Aq)
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)
  deallocate(nres)

  deallocate(zsa)
  deallocate(zsb)


  deallocate(isaz)
  deallocate(jsaz)

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(csa)
     deallocate(csb)
     deallocate(norm)
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


end subroutine zifeast_scsrgvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine difeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: difeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

  fpm(30)=122144
  call feastdefault(fpm,info)
  if (info/=0) return


!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
     call dfeast_scsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,loop,fpm(1),fpm(9),fpm(40),Emin,Emax,M0,mode,info)
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
  call difeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine difeast_scsrgv




subroutine zifeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=142244
  call feastdefault(fpm,info)
  if (info/=0) return


!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
     call zfeast_hcsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,loop,fpm(1),fpm(9),fpm(40),Emin,Emax,M0,mode,info)
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
  call zifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zifeast_hcsrgv



subroutine difeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! ====================================================================
  implicit none
  ! include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: difeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=122344
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call difeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine difeast_gcsrgv


!!$

subroutine zifeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! ====================================================================
  implicit none
  ! include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=142344
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zifeast_gcsrgv


!!$

subroutine zifeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! ====================================================================
  implicit none
  ! include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=142144
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zifeast_scsrgv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!$

subroutine difeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  !include 'f90_noruntime_interface.fi'
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
  integer,dimension(:),allocatable :: jsb
  integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: difeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  fpm(30)=122143
  call feastdefault(fpm,info)
  if (info/=0) return
  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=1.0d0
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call difeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine difeast_scsrevx





subroutine zifeast_hcsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_hcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i

  fpm(30)=142243
  call feastdefault(fpm,info)
  if (info/=0) return

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call zifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine zifeast_hcsrevx



subroutine difeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: difeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i

  fpm(30)=122343
  call feastdefault(fpm,info)
  if (info/=0) return

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))

  do i=1,n
     sb(i)=1.0d0
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call difeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)


end subroutine difeast_gcsrevx




subroutine zifeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i

  fpm(30)=142343
  call feastdefault(fpm,info)
  if (info/=0) return

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call zifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine zifeast_gcsrevx



subroutine zifeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i

  fpm(30)=142143
  call feastdefault(fpm,info)
  if (info/=0) return

  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call zifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine zifeast_scsrevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine difeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: difeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 

  fpm(30)=122142
  call feastdefault(fpm,info)
  if (info/=0) return
!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
     call dfeast_scsrev_search(UPLO,N,sa,isa,jsa,loop,fpm(1),fpm(9),fpm(40),Emin,Emax,M0,mode,info)
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
  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=1.0d0
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call difeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine difeast_scsrev








subroutine zifeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 

  fpm(30)=142242
  call feastdefault(fpm,info)
  if (info/=0) return
!!!!!!!!!!!!!! Procedure to search interval for extremal eigenvalues  
  if (fpm(40)/=0) then ! find interval for (at least) M0/2 extermal eigenvalues
     call zfeast_hcsrev_search(UPLO,N,sa,isa,jsa,loop,fpm(1),fpm(9),fpm(40),Emin,Emax,M0,mode,info)
      if (info/=0) then
              info=7
              return  ! search failed
      endif
  end if
!!!!!!!!!!!!!!  

  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  ! identity B matrix- option for standard eigenvalue problem          
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call zifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zifeast_hcsrev




subroutine difeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: difeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

  fpm(30)=122342
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=1.0d0
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1
  call difeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine difeast_gcsrev




subroutine zifeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

  fpm(30)=142342
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1

  call zifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine zifeast_gcsrev




subroutine zifeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! Eric Polizzi 2019
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
!!!    Wrapper Routine to expert routine: zifeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne

  fpm(30)=142142
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  ! identity B matrix- option for standard eigenvalue problem
  allocate(isb(n+1))
  allocate(jsb(n))
  allocate(sb(n))
  do i=1,n
     sb(i)=(1.0d0,0.0d0)
     jsb(i)=i
     isb(i)=i
  enddo
  isb(n+1)=n+1
  call zifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine zifeast_scsrev








