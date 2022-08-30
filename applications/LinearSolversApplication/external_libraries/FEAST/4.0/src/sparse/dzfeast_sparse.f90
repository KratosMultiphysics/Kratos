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

!{D,Z}FEAST_{SCSR,HCSR,GCSR}{EV,GV}


!111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE) - using MKL-PARDISO
! dfeast_scsrgvx
! zfeast_hcsrgvx
! dfeast_gcsrgvx
! zfeast_gcsrgvx
! zfeast_scsrgvx

!222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
! dfeast_scsrgv
! zfeast_hcsrgv
! dfeast_gcsrgv
! zfeast_gcsrgv
! zfeast_scsrgv

!333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
! dfeast_scsrevx
! zfeast_hcsrevx
! dfeast_gcsrevx
! zfeast_gcsrevx
! zfeast_scsrevx

!44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
! dfeast_scsrev
! zfeast_hcsrev
! dfeast_gcsrev
! zfeast_gcsrev
! zfeast_scsrev


!5555555555555!!! Routines for Estimating the Search interval Emin,Emax for the real symmetric or Hermitian problem
!dfeast_scsrev_search
!zfeast_hcsrev_search
!dfeast_scsrgv_search
!zfeast_hcsrgv_search


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration node
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2009-2019
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
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! matrix format
  double precision,dimension(:),allocatable :: usa,usb,ddiag
  integer,dimension(:), allocatable :: uisa,ujsa,uisb,ujsb
  integer :: opt,nnza,nnzb,nnz
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  real,dimension(:),allocatable :: susa,susb !single precision copy
  complex,dimension(:,:),allocatable :: cwork,caux
!!!!!for pardiso
  integer(8),dimension(:,:),allocatable :: pt 
  integer,dimension(:,:),allocatable :: iparm
  integer :: mtype
  integer :: MAXFCT,MNUM,PHASE,MSGLVL
  integer :: idum
  logical :: mkl
!!!!!!!!!!!!!!!!!
  !    integer :: t1,t2,tim

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif

  if ((fpm(43)==1).or.(.not.mkl)) then ! switch to IFEAST 
     call difeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121145 ! code name
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
     CALL XERBLA( 'DFEAST_SCSRGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Set up Az matrix !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  nfact=1
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank+1,fpm(2),nb_procs
        nfact=nfact+1
     end do
  endif


!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zdaddcsr(N,N,opt,ZONE,usa,uisa,ujsa,ZONE,usb,uisb,ujsb,saz,isaz,jsaz)
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 ! get jsaz
  call zdaddcsr(N,N,opt,ZONE,usa,uisa,ujsa,ZONE,usb,uisb,ujsb,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy       
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zdaddcsr(N,N,opt,-ZONE,usa,uisa,ujsa,Ze,usb,uisb,ujsb,saz(1,1),isaz,jsaz)

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
        do k=uisa(i),uisa(i+1)-1   
           usa(k)=usa(k)/(sqrt(ddiag(i))*sqrt(ddiag(ujsa(k))))
        enddo
        do k=uisb(i),uisb(i+1)-1
           usb(k)=usb(k)/(sqrt(ddiag(i))*sqrt(ddiag(ujsb(k))))
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
     allocate(susa(uisa(N+1)-1))
     allocate(susb(uisb(N+1)-1))
     call DLAG2S(uisa(N+1)-1,1,usa,uisa(N+1)-1,susa,uisa(N+1)-1,infoloc1)
     call DLAG2S(uisb(N+1)-1,1,usb,uisb(N+1)-1,susb,uisb(N+1)-1,infoloc1)
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


  !    call system_clock(t2,tim)
  !    print *,'time 1:',(t2-t1)*1.0d0/tim
  !     call system_clock(t1,tim)

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
     allocate(caux(N,m0))
     allocate(cwork(N,M0))
     !! Symbolic factorization
!!$       opt=3 ! get csaz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call csaddcsr(N,N,opt,-CONE,susa,uisa,ujsa,cmplx(Ze),susb,uisb,ujsb,csaz(1,i),isaz,jsaz)      
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,infoloc2)
!!$       enddo

  else ! double precision only

     allocate(zaux(N,M0))
     !! Symbolic fatorization
!!$       opt=3 ! get saz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call zdaddcsr(N,N,opt,-ZONE,usa,uisa,ujsa,Ze,usb,uisb,ujsb,saz(1,i),isaz,jsaz) 
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,infoloc2)
!!$       enddo

  end if

!!$    if (infoloc2/=0) then
!!$       info=-2
!!$       return
!!$    end if


  !    call system_clock(t2,tim)
  !print *,'time 2:',(t2-t1)*1.0d0/tim
  !   call system_clock(t1,tim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization 
  do while (ijob/=0)
     call dfeast_srcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization
        !print *,'10',id
        if (fpm(42)==1) then !single precision fact- csaz             
           call csaddcsr(N,N,opt,-CONE,susa,uisa,ujsa,cmplx(Ze),susb,uisb,ujsb,csaz(1,id),isaz,jsaz) !! get csaz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
        else ! full precision
           call zdaddcsr(N,N,opt,-ZONE,usa,uisa,ujsa,Ze,usb,uisb,ujsb,saz(1,id),isaz,jsaz) !! get saz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
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
        !           print *,'Ze',Ze
        !print *,'bef',zwork(1:N,1:M0)         

        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else ! full precision
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        !print *,'aft',zwork(1:N,1:M0)
        !stop
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

        call wdcsrmm('U','N',N,N,fpm(25),DONE,usa,uisa,ujsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call wdcsrmm('U','N',N,N,fpm(25),DONE,usb,uisb,ujsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))

     end select
  end do

  ! call system_clock(t2,tim)
  !    print *,'time 2:',(t2-t1)*1.0d0/tim
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
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

  deallocate(usa)
  deallocate(uisa)
  deallocate(ujsa)

  deallocate(usb)
  deallocate(uisb)
  deallocate(ujsb)

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(susa)
     deallocate(susb)
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

end subroutine dfeast_scsrgvx






subroutine zfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  ! Eric Polizzi 2009-2019
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
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  double precision,dimension(:),allocatable :: ddiag
  integer,dimension(:), allocatable :: fisa,fjsa,fisb,fjsb
  integer :: opt,nnza,nnzb,nnz
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
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
     call zifeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141245 ! code name
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
     CALL XERBLA( 'ZFEAST_HCSRGVX', -INFO+100 )
     RETURN
  END IF


  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 

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
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank+1,fpm(2),nb_procs
        nfact=nfact+1
     end do
  endif


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=ZONE*(Emax+Emin)/2.0d0!+ZIMAG*(Emax-Emin)/2.0d0
     call zaddcsr(N,N,opt,-ZONE,fsa,fisa,fjsa,Ze,fsb,fisb,fjsb,saz(1,1),isaz,jsaz)
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
!!!!!!!!  pardiso initialization
!!!!!!!!  use same factorization for (normal+transpose solve)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1
  MNUM=1
  MTYPE=3     ! complex and structurally symmetric
  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))
  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  end do
  !IPARM(2,:)=3 ! parallel omp nested dissection  (sensible to #threads- no consistency- same runs with different results)
  !IPARM(4)=11 !CGS solver
  IPARM(25,:)=1 ! parallel omp rhs solve
  IPARM(11,:)=0 ! disable scaling
!!!!!!!!!!!!
  if (fpm(64)==1) then
     do i=1,64
        if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
     enddo
  endif
!!!!!!!!!!!!
  IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used
  MSGLVL=0!0 !0- no output, 1- output
  !    PHASE=11 ! symbolic factorization (do it only once)

  if (fpm(42)==1) then ! mixed precision (single precision solver)       
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<
     allocate(caux(N,m0))
     allocate(cwork(N,M0))
     !! Symbolic factorization
!!$       opt=3 ! get csaz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call caddcsr(N,N,opt,-CONE,cfsa,fisa,fjsa,cmplx(Ze),cfsb,fisb,fjsb,csaz(1,i),isaz,jsaz)      
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,infoloc2)
!!$       enddo

  else ! double precision only

     allocate(zaux(N,M0))
     !! Symbolic fatorization
!!$       opt=3 ! get saz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call zaddcsr(N,N,opt,-ZONE,fsa,fisa,fjsa,Ze,fsb,fisb,fjsb,saz(1,i),isaz,jsaz) 
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,infoloc2)
!!$       enddo

  end if

!!$    if (infoloc2/=0) then
!!$       info=-2
!!$       return
!!$    end if
!!$  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_hrcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization

        if (fpm(42)==1) then !single precision fact- csaz             
           call caddcsr(N,N,opt,-CONE,cfsa,fisa,fjsa,cmplx(Ze),cfsb,fisb,fjsb,csaz(1,id),isaz,jsaz) !! get csaz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
        else ! full precision
           call zaddcsr(N,N,opt,-ZONE,fsa,fisa,fjsa,Ze,fsb,fisb,fjsb,saz(1,id),isaz,jsaz) !! get saz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
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
        !print *,'bef',zwork(1:N,1:M0) 
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else ! full precision
           !print *,'bef',zwork(1:N,1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
           !print *,'aft',zwork(1:N,1)
        end if

        !          if (loop==1) stop
        ! print *,'aft',ijob,Ze,zwork(1:2,1)

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

        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else ! full precision
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if
        !print *,'aft',zwork(1:N,1:M0)
        !stop
        ! print *,'aft',ijob,Ze,zwork(1:2,1)

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

        call wzhcsrmm('F','N',N,N,fpm(25),ZONE,fsa,fisa,fjsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        call wzhcsrmm('F','N',N,N,fpm(25),ZONE,fsb,fisb,fjsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))


     end select
  end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory pardiso
!!!!!!!!!!!!!!!!!!!!!!!!

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
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)

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

end subroutine zfeast_hcsrgvx








subroutine dfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! James Kestyn, Eric Polizzi 2015 
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
  !complex(kind=(kind(1.0d0))),dimension(:),allocatable :: saz0,sbz0
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  double precision,dimension(:),allocatable :: ddiag
  !integer,dimension(:), allocatable :: fisa,fjsa,fisb,fjsb
  integer :: opt,nnza,nnzb,nnz
  !transpose matrix
  !double precision,dimension(:),allocatable :: tsa,tsb
  !integer,dimension(:), allocatable :: tisa,tjsa,tisb,tjsb
  !integer :: opt,nnz
!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
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
     call difeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    


#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121345 ! code name
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
     CALL XERBLA( 'DFEAST_GCSRGVX', -INFO+100 )
     RETURN
  END IF



  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT MATRIX CONVERSION 
!!! A and B already in full format
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
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  endif


!!!!!!!!!!!!!!!!! Set up for Az matrix
  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zdaddcsr(N,N,opt,ZONE,sa,isa,jsa,ZONE,sb,isb,jsb,saz,isaz,jsaz) 
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !! et jsaz
  call zdaddcsr(N,N,opt,ZONE,sa,isa,jsa,ZONE,sb,isb,jsb,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy  
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid!+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
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
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1
  MNUM=1
  MTYPE=13  ! complex and unsymmetric
  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))
  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  end do
  !IPARM(2,:)=3 ! parallel omp nested dissection  (sensible to #threads- no consistency- same runs with different results)
  !IPARM(4)=11 !CGS solver
  IPARM(25,:)=1 ! parallel omp rhs solve
  IPARM(11,:)=0 ! disable scaling
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

     !! Symbolic factorization
!!$       opt=3 ! get csaz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=Emid+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call caddcsr(N,N,opt,-CONE,cfsa,isa,jsa,cmplx(Ze),cfsb,isb,jsb,csaz(1,i),isaz,jsaz)      
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,infoloc2)
!!$       enddo

  else ! double precision only
     allocate(zaux(N,M0))
!!$ do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=Emid+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,i),isaz,jsaz)    
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,infoloc2)
!!$       enddo     
  end if
!!$
!!$    if (infoloc2/=0) then
!!$       info=-2
!!$       return
!!$    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization
  do while (ijob/=0) 
     call dfeast_grcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization

        if (fpm(42)==1) then !single precision fact- csaz             
           call caddcsr(N,N,opt,-CONE,cfsa,isa,jsa,cmplx(Ze),cfsb,isb,jsb,csaz(1,id),isaz,jsaz) !! get csaz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
        else ! full precision     
           call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        PHASE=33 ! solve
        IPARM(12,id)=0 ! normal solve
    !print *,id,sum(abs(zwork(1:N,1))),sum(abs(saz(1:100,id)))
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else ! full precision
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if
!print *,sum(abs(zwork(1:N,1))),infoloc2
        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsa,isa,jsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsa,isa,jsa,X(1,fpm(34)),ZZERO,work(1,fpm(34)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !print *,'bef',X(1:2,1)
        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsb,isb,jsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))
        !print *,'aft',work(1:2,1)
        !stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsb,isb,jsb,X(1,fpm(34)),ZZERO,work(1,fpm(34)))


     end select
  end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory pardiso
!!!!!!!!!!!!!!!!!!!!!!!!

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
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)

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

end subroutine dfeast_gcsrgvx





subroutine zfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! James Kestyn, Eric Polizzi 2015
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
  !complex(kind=(kind(1.0d0))),dimension(:),allocatable :: saz0,sbz0
  integer,dimension(:),allocatable :: isaz,jsaz
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO),ZIMAG=(DZERO,DONE)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! matrix format
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: fsa,fsb
  double precision,dimension(:),allocatable :: ddiag
  integer :: opt,nnza,nnzb,nnz
!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: cfsa,cfsb !single precision copy
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
     call zifeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    


#ifdef MKL  

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141345 ! code name
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
     CALL XERBLA( 'ZFEAST_GCSRGVX', -INFO+100 )
     RETURN
  END IF



  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT MATRIX CONVERSION 
!!! A and B already in full format
!!!!!!!!!!!! copy sa,sb (needed for scaling)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Perform a copy of the matrix values 
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
  if (fpm(10)==1) then !! store factorization
     nfact=0
     !! nfact is local (number of total factorization by contour points)
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  endif


!!!!!!!!!!!!!!!!! Set up for Az matrix

  allocate(isaz(n+1))
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zaddcsr(N,N,opt,ZONE,fsa,isa,jsa,ZONE,fsb,isb,jsb,saz,isaz,jsaz) 
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 !! et jsaz
  call zaddcsr(N,N,opt,ZONE,fsa,isa,jsa,ZONE,fsb,isb,jsb,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy  
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid!+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
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
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAXFCT=1
  MNUM=1
  MTYPE=13  ! complex and unsymmetric
  allocate(pt(64,nfact))
  allocate(iparm(64,nfact))
  do i=1,nfact !! multiple factorization
     pt(:,i)=0
     call pardisoinit(PT(1,i),MTYPE,IPARM(1,i))
  end do
  !IPARM(2,:)=3 ! parallel omp nested dissection  (sensible to #threads- no consistency- same runs with different results)
  !IPARM(4)=11 !CGS solver
  IPARM(25,:)=1 ! parallel omp rhs solve
  IPARM(11,:)=0 ! disable scaling
!!!!!!!!!!!!
  if (fpm(64)==1) then
     do i=1,64
        if (fpm(64+i)/=-111) iparm(i,:)=fpm(64+i)
     enddo
  endif
!!!!!!!!!!!!
  IPARM(6,:)=1 ! solution and rhs are input/output, attention zaux is always used
  MSGLVL=0!0 !0- no output, 1- output
  PHASE=11 ! symbolic factorization (do it only once)

  if (fpm(42)==1) then ! mixed precision (single precision solver)
     IPARM(28,:)=1 ! pardiso single precision !!!<<<<
     allocate(caux(N,M0))
     allocate(cwork(N,M0))

     !! Symbolic factorization
!!$       opt=3 ! get csaz
!!$       do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=Emid+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call caddcsr(N,N,opt,-CONE,cfsa,isa,jsa,cmplx(Ze),cfsb,isb,jsb,csaz(1,i),isaz,jsaz)      
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,cwork,caux,infoloc2)
!!$       enddo

  else ! double precision only
     allocate(zaux(N,M0))
!!$ do i=1,nfact
!!$          ! initialize Az matrix for pardiso symbolic fatorization
!!$          Ze=Emid+ZIMAG!ZONE*(Emax+Emin)/2.0d0+ZIMAG*(Emax-Emin)/2.0d0
!!$          call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,i),isaz,jsaz)    
!!$          call PARDISO(PT(1,i),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,i),isaz,jsaz,idum,M0,IPARM,MSGLVL,zwork,zaux,infoloc2)
!!$       enddo     
  end if
!!$
!!$    if (infoloc2/=0) then
!!$       info=-2
!!$       return
!!$    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_grcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization
        if (fpm(42)==1) then !single precision fact- csaz             
           call caddcsr(N,N,opt,-CONE,cfsa,isa,jsa,cmplx(Ze),cfsb,isb,jsb,csaz(1,id),isaz,jsaz) !! get csaz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
        else ! full precision     
           call zaddcsr(N,N,opt,-ZONE,fsa,isa,jsa,Ze,fsb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsa,isa,jsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsa,isa,jsa,X(1,fpm(34)),ZZERO,work(1,fpm(34)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !print *,'bef',X(1:2,1)
        call wzcsrmm('F','N',N,N,fpm(25),ZONE,fsb,isb,jsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))
        !print *,'aft',work(1:2,1)
        !stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('F','C',N,N,fpm(35),ZONE,fsb,isb,jsb,X(1,fpm(34)),ZZERO,work(1,fpm(34)))


     end select
  end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory pardiso
!!!!!!!!!!!!!!!!!!!!!!!!

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
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)

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

end subroutine zfeast_gcsrgvx






subroutine zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! Eric Polizzi 2009-2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N
  complex(kind=(kind(1.0d0))),dimension(*)  :: sa,sb
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
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZIMAG=(DZERO,DONE),ZZERO=(DZERO,DZERO)
  real,parameter ::  SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! matrix format
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: usa,usb
  double precision,dimension(:),allocatable :: ddiag
  integer,dimension(:), allocatable :: uisa,ujsa,uisb,ujsb
  integer :: opt,nnza,nnzb,nnz
!!!!!!!!!!!!!!!!!!! for mixed precision
  complex,dimension(:,:),allocatable :: csaz
  complex,dimension(:),allocatable :: susa,susb !single precision copy
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
     call zifeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
     return ! leave
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!  solver using intel mkl pardiso
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
#ifdef MKL

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141145 ! code name
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
     CALL XERBLA( 'ZFEAST_SCSRGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  allocate(saz(1,1)) ! dummy
  allocate(jsaz(1))! dummy
  opt=1 ! get isaz
  call zaddcsr(N,N,opt,ZONE,usa,uisa,ujsa,ZONE,usb,uisb,ujsb,saz,isaz,jsaz)
  nnz=isaz(n+1)-1
  deallocate(jsaz)
  allocate(jsaz(nnz))
  opt=2 ! get jsaz
  call zaddcsr(N,N,opt,ZONE,usa,uisa,ujsa,ZONE,usb,uisb,ujsb,saz,isaz,jsaz)
  deallocate(saz) ! deallocate dummy       
!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  SCALE THE MATRIX !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(41)==1) then ! if scaling on
     allocate(saz(nnz,1)) ! dummy
     opt=3 !! get saz with given shift to build scaling
     Ze=Emid
     call zaddcsr(N,N,opt,-ZONE,usa,uisa,ujsa,Ze,usb,uisb,ujsb,saz(1,1),isaz,jsaz)

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
        do k=uisa(i),uisa(i+1)-1   
           usa(k)=usa(k)/(sqrt(ddiag(i))*sqrt(ddiag(ujsa(k))))
        enddo
        do k=uisb(i),uisb(i+1)-1
           usb(k)=usb(k)/(sqrt(ddiag(i))*sqrt(ddiag(ujsb(k))))
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
     allocate(susa(uisa(N+1)-1))
     allocate(susb(uisb(N+1)-1))
     call ZLAG2C(uisa(N+1)-1,1,usa,uisa(N+1)-1,susa,uisa(N+1)-1,infoloc1)
     call ZLAG2C(uisb(N+1)-1,1,usb,uisb(N+1)-1,susb,uisb(N+1)-1,infoloc1)
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
     allocate(caux(N,m0))
     allocate(cwork(N,M0))
  else ! double precision only
     allocate(zaux(N,M0))
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_srcix(ijob,N,Ze,work,zwork,Aq,Sq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! form and factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        id=fpm(33) !! id of factorization (for fpm(10) flag)
        opt=3 
        PHASE=12 !include the symbolic factorization
        if (fpm(42)==1) then !single precision fact- csaz           
           call caddcsr(N,N,opt,-CONE,susa,uisa,ujsa,cmplx(Ze),susb,uisb,ujsb,csaz(1,id),isaz,jsaz) !! get csaz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,csaz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,cwork,caux,infoloc2)
        else ! full precision
           call zaddcsr(N,N,opt,-ZONE,usa,uisa,ujsa,Ze,usb,uisb,ujsb,saz(1,id),isaz,jsaz) !! get saz
           call PARDISO(PT(1,id),MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM(1,id),MSGLVL,zwork,zaux,infoloc2)
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
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wzcsrmm('U','N',N,N,fpm(25),ZONE,usa,uisa,ujsa,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call wzcsrmm('U','N',N,N,fpm(25),ZONE,usb,uisb,ujsb,X(1,fpm(24)),ZZERO,work(1,fpm(24)))

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
  deallocate(Sq)
  deallocate(work)
  deallocate(zwork)


  deallocate(isaz)
  deallocate(jsaz)

  deallocate(usa)
  deallocate(uisa)
  deallocate(ujsa)

  deallocate(usb)
  deallocate(uisb)
  deallocate(ujsb)

  if (fpm(42)==1) then ! mixed precision
     deallocate(cwork)
     deallocate(caux)
     deallocate(csaz)
     deallocate(susa)
     deallocate(susb)
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

end subroutine zfeast_scsrgvx






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine dfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! Eric Polizzi 2009-2019
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
!!!    Wrapper Routine to expert routine: dfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!


  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121144
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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
  call dfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine dfeast_scsrgv




subroutine zfeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! Eric Polizzi 2009-2019
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
!!!    Wrapper Routine to expert routine: zfeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141244
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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
  call zfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zfeast_hcsrgv



subroutine dfeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: dfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121344
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine dfeast_gcsrgv


!!$

subroutine zfeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! James Kestyn, Eric Polizzi 2015
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
  fpm(30)=141344
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zfeast_gcsrgv


!!$

subroutine zfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: zfeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141144
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_scsrgv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!$

subroutine dfeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  ! Eric Polizzi 2009-2019
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
!!!    Wrapper Routine to expert routine: dfeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121143
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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


  call dfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine dfeast_scsrevx





subroutine zfeast_hcsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
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
  ! Eric Polizzi 2009-2019
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
!!!    Wrapper Routine to expert routine: zfeast_hcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif

  fpm(30)=141243
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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

  call zfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine zfeast_hcsrevx



subroutine dfeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: dfeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121343
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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

  call dfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)


end subroutine dfeast_gcsrevx




subroutine zfeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: zfeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif

  fpm(30)=141343
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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

  call zfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine zfeast_gcsrevx



subroutine zfeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: zfeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141143
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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

  call zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)

end subroutine zfeast_scsrevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine dfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! Eric Polizzi 2009-2015
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
!!!    Wrapper Routine to expert routine: dfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  logical :: mkl
!!!!!!!!!!!!!!!!!
 
  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121142
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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

  call dfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine dfeast_scsrev








subroutine zfeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
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
  ! Eric Polizzi 2009-2015
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
!!!    Wrapper Routine to expert routine: zfeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141242
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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

  call zfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zfeast_hcsrev




subroutine dfeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: dfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision,dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=121342
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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
  call dfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)

end subroutine dfeast_gcsrev




subroutine zfeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: zfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141342
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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

  call zfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine zfeast_gcsrev




subroutine zfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
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
  ! James Kestyn, Eric Polizzi 2015
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
!!!    Wrapper Routine to expert routine: zfeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: sb ! identity
  integer,dimension(:),allocatable :: isb
  integer,dimension(:),allocatable :: jsb
  integer :: i
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne
  logical :: mkl
!!!!!!!!!!!!!!!!!

  mkl=.false. ! initialization
#ifdef MKL
  mkl=.true.
#endif
  fpm(30)=141142
  if ((fpm(43)==1).or.(.not.mkl))  fpm(30)=fpm(30)+1000 ! switch to IFEAST
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
  call zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(sb)
  deallocate(isb)
  deallocate(jsb)
  deallocate(Zne)
  deallocate(Wne)


end subroutine zfeast_scsrev






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!5555555555555555555555555555555555555555555555555555555555555555555555555555555
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine dfeast_scsrev_search(UPLO,N,sa,isa,jsa,itloop,fpm1,fpm9,fpm40,Emin,Emax,M0,M,info)
  !  Purpose 
  !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Standard eigenvalue problem version: AX=XE  where A is real symmetric 
  ! 
  !  DOUBLE PRECISION version
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
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================  
  implicit none
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40

!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  double precision,dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: NEW_COMM_WORLD,rank,code,nb_procs,infofeast
  character(len=3) :: ctemp
  integer(8),save :: fout
 
  info=0 
  rank=0
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm9
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
  !---------------------------------------------

  call dcheck_feast_srci_input(0.0d0,1.0d0,M0,N,info)
  if (info/=0) return

  ! comments?
  com=.false. ! default
  if (fpm1/=0) com=.true.
  if (com.and.(rank/=0)) com=.false. !only rank 0 is commenting

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
     write(fout,'(A)') ' Routine dfeast_scsrev_search'
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
  itmax=min(50,N) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge
  allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(N,itmax)) 

  it=itmax
  call darnoldi(UPLO,N,sa,isa,jsa,dE,dX,dres,eps,it) 
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
     info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
     !! search has failed (ill-conditionned matrix?)
     return
  end if

  itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
  Emin=dE(min(3,N))-(dE(min(3,N))-dE(1))*1.05d0
  Emax=dE(itmax-min(n-1,2))+(dE(itmax)-dE(itmax-min(n-1,2)))*1.05d0

  if (com) then
     write(fout,*)
     write(fout,'(A)') ' Start FEAST stochastic search:'
     write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
     i=0
     write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,N
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

  MS=min(10,N) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1
  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1

     call dfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)
       if (infofeast/=5) then ! stochastic estimate
     info=71
     return
  end if
     

     if (com)  write(fout,'(I2,X,2ES25.16,X,I7)') itloop,Emin,Emax,M

     if (init) then 
        if (M <Mw) then ! not enough eigenvalue found (initial Emin-Emax not well set-up)
           if (fpm40 < 0) then
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

  MS=min(M0,N) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  call feastinit(fpm)
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
  call dfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

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


end subroutine dfeast_scsrev_search




subroutine zfeast_hcsrev_search(UPLO,N,sa,isa,jsa,itloop,fpm1,fpm9,fpm40,Emin,Emax,M0,M,info)
  !  Purpose 
  !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Standard eigenvalue problem version: AX=XE  where A is complex Hermitian 
  ! 
  !  DOUBLE PRECISION version
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
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ==================================================================== 
  implicit none
  character(len=1) :: UPLO
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40

!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: NEW_COMM_WORLD,rank,code,nb_procs,infofeast
  character(len=3) :: ctemp
  integer(8),save :: fout
 
  info=0 
  rank=0
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm9
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
  !---------------------------------------------

  call dcheck_feast_srci_input(0.0d0,1.0d0,M0,N,info)
  if (info/=0) return

  ! comments?
  com=.false. ! default
  if (fpm1/=0) com=.true.
  if (com.and.(rank/=0)) com=.false. !only rank 0 is commenting

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
     write(fout,'(A)') ' Routine zfeast_hcsrev_search'
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
  itmax=min(50,N) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge
  allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(N,itmax)) 

  it=itmax
  call zharnoldi(UPLO,N,sa,isa,jsa,dE,dX,dres,eps,it) 
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
     info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
     !! search has failed (ill-conditionned matrix?)
     return
  end if

  itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
  Emin=dE(min(3,N))-(dE(min(3,N))-dE(1))*1.05d0
  Emax=dE(itmax-min(n-1,2))+(dE(itmax)-dE(itmax-min(n-1,2)))*1.05d0

  if (com) then
     write(fout,*)
     write(fout,'(A)') ' Start FEAST stochastic search:'
     write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
     i=0
     write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,N
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

  MS=min(10,N) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1

  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1

     call zfeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)
       if (infofeast/=5) then ! stochastic estimate
     info=71
     return
  end if
     

     if (com)  write(fout,'(I2,X,2ES25.16,X,I7)') itloop,Emin,Emax,M

     if (init) then 
        if (M <Mw) then ! not enough eigenvalue found (initial Emin-Emax not well set-up)
           if (fpm40 < 0) then
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

  MS=min(M0,N) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  call feastinit(fpm)
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
  call zfeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

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


end subroutine zfeast_hcsrev_search



subroutine dfeast_scsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,itloop,fpm1,fpm9,fpm40,Emin,Emax,M0,M,info)
 !  Purpose 
  !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Generalized eigenvalue problem version: AX=BXE  where A is real symmetric, B is spd 
  ! 
  !  DOUBLE PRECISION version
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
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================   
  implicit none
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  double precision,dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: NEW_COMM_WORLD,rank,code,nb_procs,infofeast
  character(len=3) :: ctemp
  integer(8),save :: fout

  info=0
rank=0
 !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm9
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
  !---------------------------------------------


 call dcheck_feast_srci_input(0.0d0,1.0d0,M0,N,info)
 if (info/=0) return
 
  ! comments?
  com=.false. ! default
  if (fpm1/=0) com=.true.
  if (com.and.(rank/=0)) com=.false. !only rank 0 is commenting

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
  if (com) then   
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') '***** FEAST Search Estimate for Emin-Emax *****'
     write(fout,'(A)') '*** looking for the M0/2 extremal eigenvalue **' 
     write(fout,'(A)') '***********************************************'
     write(fout,'(A)') ' Routine dfeast_scsrgv_search'
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
  itmax=min(50,N) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge
   allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(N,itmax)) 

  it=itmax
  call dgarnoldi(UPLO,N,sa,isa,jsa,sb,isb,jsb,dE,dX,dres,eps,it,1d-3,30) ! last 2 arguments for bicgstab for B matrix
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

!!!!!!!!!!!!!!!!!!

if ((it==itmax).and.(dres(1)>eps).and.(dres(it)>eps)) then
   info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
          !! search has failed (ill-conditionned matrix?)
   return
end if

itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
Emin=dE(min(3,N))-(dE(min(3,N))-dE(1))*1.05d0
Emax=dE(itmax-min(n-1,2))+(dE(itmax)-dE(itmax-min(n-1,2)))*1.05d0

if (com) then
write(fout,*)
 write(fout,'(A)') ' Start FEAST stochastic search:'
 write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
 i=0
 write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,N
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

!!!!!! Initializing the stochastic search
  MW = max(M0/2,1) ! number of wanted eigenvalue
  call feastinit(fpm)

  MS=min(10,N) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1

  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1
    
     call dfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)
      if (infofeast/=5) then ! stochastic estimate
     info=71
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
        if (fpm40<0) then
           Emax = Emax - diff
        else
           Emin = Emin + diff
        endif

     else if ((M < MW).and.(.not.init)) then
        diff=diff/2.0d0
        if (fpm40<0) then
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


  MS=min(M0,N) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  call feastinit(fpm)
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
   call dfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

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


end subroutine dfeast_scsrgv_search





subroutine zfeast_hcsrgv_search(UPLO,N,sa,isa,jsa,sb,isb,jsb,itloop,fpm1,fpm9,fpm40,Emin,Emax,M0,M,info)
   !  Purpose 
  !  =======
  !  Return an estimation of Emin-Emax that contains the M0/2 lowest or largest eigenvalues
  !
  !  Generalized eigenvalue problem version: AX=BXE  where A is complex Hermitian 
  ! 
  !  DOUBLE PRECISION version
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
  !  Emin,Emax  (output)       DOUBLE PRECISION: Search interval (estimation)
  !  M0         (input)        INTEGER: Size subspace- should be two times the # of wanted eigenvalues
  !   M         (output)       INTEGER : # of eigenvalues found in the search interval (should be M0/2)                          
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2019
  ! ==================================================================== 
  implicit none
  character(len=1) :: UPLO
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  integer :: itloop
  double precision :: Emin,Emax
  integer :: M0,M
  integer :: info,fpm1,fpm9,fpm40

!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: MS,MW,it,itmax,loop,i
  double precision,dimension(:),allocatable :: dE,dres
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: dX
  double precision :: diff,depsout,eps
  integer,dimension(64) :: fpm
  logical :: init,com
  integer :: NEW_COMM_WORLD,rank,code,nb_procs,infofeast
  character(len=3) :: ctemp
  integer(8),save :: fout
 
  info=0 
  rank=0
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm9
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
  !---------------------------------------------

  call dcheck_feast_srci_input(0.0d0,1.0d0,M0,N,info)
  if (info/=0) return

  ! comments?
  com=.false. ! default
  if (fpm1/=0) com=.true.
  if (com.and.(rank/=0)) com=.false. !only rank 0 is commenting

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
     write(fout,'(A)') ' Routine zfeast_hcsrgv_search'
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
  itmax=min(50,N) ! 50 iteration should be enough for Arnoldi iteration to converge to 1E-2 (for one extreme value) and provide finer spectrum at the edge
  allocate(dE(itmax))     
  allocate(dres(itmax))   
  allocate(dX(N,itmax)) 

  it=itmax
  call zhgarnoldi(UPLO,N,sa,isa,jsa,sb,isb,jsb,dE,dX,dres,eps,it,1d-3,30) ! last 2 arguments for bicgstab for B matrix
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

!  info=-1
!  return
  

if ((it==itmax).and.(dres(1)>eps).and.(dres(it)>eps)) then  
     info=70 !! unfortunately Arnoldi was unable to reach eps in itmax iteration
     !! search has failed (ill-conditionned matrix?)
     return
  end if

  itmax=it

!!!!!! we want to "increase Emin/Emax by 5%" (we use here a technique to take care of possible E1 or Eitmax close to zero)
  Emin=dE(min(3,N))-(dE(min(3,N))-dE(1))*1.05d0
  Emax=dE(itmax-min(n-1,2))+(dE(itmax)-dE(itmax-min(n-1,2)))*1.05d0

  if (com) then
     write(fout,*)
     write(fout,'(A)') ' Start FEAST stochastic search:'
     write(fout,'(A)') ' It.|          Emin          |          Emax          |  #eig.'
     i=0
     write(fout,'(I2,X,2ES25.16,X,I7)') i,Emin,Emax,N
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

  MS=min(10,N) ! (10 is enough for the stochastic search)
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  fpm(14)=2
  fpm(45)=2
  fpm(43)=1

  M=0  
  itloop=0
  init=.true.
  do while ((M < MW).or.(M > M0))
     itloop=itloop+1

     call zfeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)
       if (infofeast/=5) then ! stochastic estimate
     info=71
     return
  end if
     
     if (com)  write(fout,'(I2,X,2ES25.16,X,I7)') itloop,Emin,Emax,M

     if (init) then 
        if (M <Mw) then ! not enough eigenvalue found (initial Emin-Emax not well set-up)
           if (fpm40 < 0) then
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

  MS=min(M0,N) 
  allocate(dE(MS))     
  allocate(dres(MS))   
  allocate(dX(N,MS))

  call feastinit(fpm)
  fpm(4)=10  ! 10 iterations max
  fpm(3)=2
  fpm(6)=0 ! cv on trace (let the eigenvalue number the time to stabilize)
  fpm(43)=1 ! ifeast
  call zfeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,depsout,loop,Emin,Emax,MS,dE,dX,M,dres,infofeast)

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

end subroutine zfeast_hcsrgv_search












