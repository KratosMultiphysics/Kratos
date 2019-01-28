  !=========================================================================================
  !Copyright (c) 2009-2015, The Regents of the University of Massachusetts, Amherst.
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
  
  include 'sclsprim.f90' !! Sparse primitives





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED SPARSE INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! SINGLE PRECISION VERSION

  ! List of routines:
  !-------------------

  !{D,Z}FEAST_{SCSR,HCSR,GCSR}{EV,GV}


  !111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE)
  ! sfeast_scsrgvx
  ! cfeast_hcsrgvx
  ! sfeast_gcsrgvx
  ! cfeast_gcsrgvx
  ! cfeast_scsrgvx

  !222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
  ! sfeast_scsrgv
  ! cfeast_hcsrgv
  ! sfeast_gcsrgv
  ! cfeast_gcsrgv
  ! cfeast_scsrgv

  !333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
  ! sfeast_scsrevx
  ! cfeast_hcsrevx
  ! sfeast_gcsrevx
  ! cfeast_gcsrevx
  ! cfeast_scsrevx

  !44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
  ! sfeast_scsrev
  ! cfeast_hcsrev
  ! sfeast_gcsrev
  ! cfeast_gcsrev
  ! cfeast_scsrev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex :: Ze
    complex, dimension(:,:),pointer ::workc,caux
    real, dimension(:,:),pointer ::work,Aq,Sq
    complex,dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
    real,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: opt,nnza,nnzb,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------

    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_SCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)

       !!<<<
       call wallocate_1s(ssa,1,infoloc) ! dummy
       call wallocate_1i(sjsa,1,infoloc) !dummy
       call wallocate_1s(ssb,1,infoloc) !dummy
       call wallocate_1i(sjsb,1,infoloc) !dummy
       if (infoloc/=0) then
          info=-1
          return
       end if
       !!>>>
       opt=1
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call scsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 
       !!<<<
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sjsa)
       call wdeallocate_1s(ssb)
       call wdeallocate_1i(sjsb)
       !!>>>

       call wallocate_1s(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1s(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call scsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       !       nnza=isa(n+1)-1
       !       ssa => sa(1:nnza)
       !       sisa => isa(1:n+1)
       !       sjsa => jsa(1:nnza)

       !       nnzb=isb(n+1)-1
       !       ssb =>  sb(1:nnzb)
       !       sisb => isb(1:n+1)
       !       sjsb => jsb(1:nnzb)


    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1s(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1s(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call scsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call scsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2s(Aq,M0,M0,infoloc)
    call wallocate_2s(Sq,M0,M0,infoloc)
    call wallocate_2s(work,N,M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(2),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

    call wallocate_1i(isaz,n+1,infoloc)
    call wallocate_2c(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='U').or.(UPLO=='u')) then
       call csaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    else
       call csaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get isaz
    end if
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2c(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='U').or.(UPLO=='u')) then
       call csaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
    else
       call csaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get jsaz
    end if
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)


    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=nfact ! Rq: same factorization for (normal+transpose)
    MTYPE=6      ! complex and symmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
    if (fpm(64)==1) then
       do i=1,64
          if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
       enddo
    endif
!!!!!!!!!!!!
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0!0 !0- no output, 1- output
    PHASE=11 ! symbolic factorization (do it only once)
!!!!!!!! single precision pardiso (MKL only)
    IPARM(28)=1 ! pardiso single precision
!!!!!!!!!!!!!!!
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


    ijob=-1 ! initialization 
    do while (ijob/=0)
       call sfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

       select case(ijob)
       case(10) !! factorize (zeB-A)


          if (fpm(10)==1) then
             id=id+1
             if (id==nfact+1) then
                id=1
                fact=.false.
             endif
          endif

          if (fact) then
             opt=3
             if ((UPLO=='U').or.(UPLO=='u')) then
                call csaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call csaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
             end if

             if (PHASE==11) then
                PHASE=12
             else
                PHASE=22
             endif
             !             PHASE=22 ! numerical fact only
             MNUM=id
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

          endif ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)        
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='U').or.(UPLO=='u')) then
             call scsrmm('U','N',N,N,fpm(25),SONE,sa,isa,jsa,X(1,fpm(24)),SZERO,work(1,fpm(24)))
          else
             call scsrmm('U','N',N,N,fpm(25),SONE,ssa,sisa,sjsa,X(1,fpm(24)),SZERO,work(1,fpm(24)))
          endif

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='U').or.(UPLO=='u')) then
             call scsrmm('U','N',N,N,fpm(25),SONE,sb,isb,jsb,X(1,fpm(24)),SZERO,work(1,fpm(24)))
          else
             call scsrmm('U','N',N,N,fpm(25),SONE,ssb,sisb,sjsb,X(1,fpm(24)),SZERO,work(1,fpm(24)))
          end if

       end select
    end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    do MNUM=1,nfact
       call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,MNUM),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
       if (infoloc/=0) then
          info=-2
          return
       end if
    end do

    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2s(work)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)


    call wdeallocate_2c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1s(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1s(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif

  end subroutine sfeast_scsrgvx





  subroutine cfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex :: Ze
    complex, dimension(:,:),pointer ::work,workc,zAq,zSq,caux
    complex,dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! full csr format
    complex(kind=kind(1.0e0)),dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: opt,nnza,nnzb,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------

    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_HCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       ! nnza=isa(n+1)-1
       ! ssa => sa(1:nnza)
       ! sisa => isa(1:n+1)
       ! sjsa => jsa(1:nnza)

       ! nnzb=isb(n+1)-1
       ! ssb => sb(1:nnzb)
       ! sisb => isb(1:n+1)
       ! sjsb => jsb(1:nnzb)

    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) ! may be overestimated
       nnzb=2*(isb(n+1)-1)
       call wallocate_1c(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1c(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call chcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
       call chcsr_uplo_to_csr(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work,N,M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(2),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

    call wallocate_1i(isaz,n+1,infoloc)
    call wallocate_2c(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='F').or.(UPLO=='f')) then
       call caddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    else
       call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get isaz
    endif

    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2c(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='F').or.(UPLO=='f')) then
       call caddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
    else
       call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get jsaz
    endif

!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)

    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=nfact ! Rq: same factorization for (normal+transpose)
    MTYPE=3      ! complex and structurally symmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
    if (fpm(64)==1) then
       do i=1,64
          if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
       enddo
    endif
!!!!!!!!!!!!
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0!0 !0- no output, 1- output
    PHASE=11 ! symbolic factorization (do it only once)
!!!!!!!! single precision pardiso (MKL only)
    IPARM(28)=1 ! pardiso single precision
!!!!!!!!!!!!!!!
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ijob=-1 ! initialization
    do while (ijob/=0) 
       call cfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

       select case(ijob)
       case(10) !! factorize (zeB-A)

          if (fpm(10)==1) then
             id=id+1
             if (id==nfact+1) then
                id=1
                fact=.false.
             endif
          endif

          if (fact) then
             opt=3
             if ((UPLO=='F').or.(UPLO=='f')) then
                call caddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
             endif

             if (PHASE==11) then
                PHASE=12
             else
                PHASE=22
             endif
             !             PHASE=22 ! numerical fact only
             MNUM=id
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

          end if ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          IPARM(12)=0 ! normal solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)    

          if (infoloc/=0) then
             info=-2
             return
          end if

       case(21) !!solve the linear system (ZeB-A)^H x=workc(1:N,1:fpm(23)) result in to workc

          IPARM(12)=1 ! transpose conjugate solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)    
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='F').or.(UPLO=='f')) then
             call chcsrmm('F','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call chcsrmm('F','N',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          endif

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='F').or.(UPLO=='f')) then
             call chcsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call chcsrmm('F','N',N,N,fpm(25),ONEC,ssb,sisb,sjsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          endif

       end select
    end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    do MNUM=1,nfact
       call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,MNUM),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
       if (infoloc/=0) then
          info=-2
          return
       end if
    end do



    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)


    call wdeallocate_2c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1c(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1c(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif


  end subroutine cfeast_hcsrgvx



!!$



  subroutine sfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL REAL :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    integer :: N
    real,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex :: Ze
    complex, dimension(:,:),pointer ::work,workc,zAq,zSq,caux
    complex,dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! transpose matrix
    !    real,dimension(:),pointer :: tsa,tsb
    !    integer,dimension(:), pointer :: tisa,tjsa,tisb,tjsb
    integer :: opt,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------


    INFO = 0
    IF( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_GCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT MATRIX CONVERSION 
!!! A and B already in full format
!!! forming tA and tB Transpose for case 31/41 
    !    call wallocate_1s(tsa,isa(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsa,isa(n+1)-1,infoloc)
    !    call wallocate_1i(tisa,n+1,infoloc)
    !    call wallocate_1s(tsb,isb(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsb,isb(n+1)-1,infoloc)
    !    call wallocate_1i(tisb,n+1,infoloc)
    !    if (infoloc/=0) then
    !       info=-1
    !       return
    !    end if
    !
    !    call scsr_transpose(n,sa,isa,jsa,tsa,tisa,tjsa)
    !    call scsr_transpose(n,sb,isb,jsb,tsb,tisb,tjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work,N,2*M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

!!!!!!!!!!!!!!!! find fpm(22) - # total number of contour points
    fpm(22)=fpm(8) ! default
    if (fpm(29)==1) then! default for generated contour -- half-contour capability case 
       if ((aimag(Emid)==0.0e0).and.((2*(fpm(8)/2))==fpm(8))) fpm(22)=fpm(8)/2  
    end if


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(22),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

    call wallocate_1i(isaz,n+1,infoloc)
    call wallocate_2c(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    call csaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2c(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    call csaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)


    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=nfact  ! Rq: same factorization for (normal+transpose)
    MTYPE=13      ! complex and unsymmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
    if (fpm(64)==1) then
       do i=1,64
          if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
       enddo
    endif
!!!!!!!!!!!!
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0!0 !0- no output, 1- output
    PHASE=11 ! symbolic factorization (do it only once)
!!!!!!!! single precision pardiso (MKL only)
    IPARM(28)=1 ! pardiso single precision
!!!!!!!!!!!!!!!
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ijob=-1 ! initialization
    do while (ijob/=0) 
       call sfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

       select case(ijob)

       case(10) !! factorize (zeB-A)

          if (fpm(10)==1) then
             id=id+1
             if (id==nfact+1) then
                id=1
                fact=.false.
             endif
          endif

          if (fact) then
             opt=3
             call csaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz


             if (PHASE==11) then
                PHASE=12
             else
                PHASE=22
             endif

             !            PHASE=22 ! numerical fact only
             MNUM=id
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          IPARM(12)=0 ! normal solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)  
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          IPARM(12)=1 ! transpose conjugate solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)  
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call sccsrmm('F','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))

       case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call sccsrmm('F','C',N,N,fpm(35),ONEC,sa,isa,jsa,X(1,fpm(34)),ZEROC,work(1,fpm(34)))


       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call sccsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))


       case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call sccsrmm('F','C',N,N,fpm(35),ONEC,sb,isb,jsb,X(1,fpm(34)),ZEROC,work(1,fpm(34)))


       end select
    end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    do MNUM=1,nfact
       call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,MNUM),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
       if (infoloc/=0) then
          info=-2
          return
       end if
    end do


    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)

    call wdeallocate_2c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    !    call wdeallocate_1s(tsa)
    !    call wdeallocate_1i(tisa)
    !    call wdeallocate_1i(tjsa)
    !    call wdeallocate_1s(tsb)
    !    call wdeallocate_1i(tisb)
    !    call wdeallocate_1i(tjsb)

  end subroutine sfeast_gcsrgvx




  subroutine cfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    integer :: N
    complex,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex :: Ze
    complex, dimension(:,:),pointer ::work,workc,zAq,zSq,caux
    complex,dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! transpose conjugate matrix
    !    complex,dimension(:),pointer :: tsa,tsb
    !    integer,dimension(:), pointer :: tisa,tjsa,tisb,tjsb
    integer :: opt,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------


    INFO = 0
    IF( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_GCSRGV', -INFO+100 )
       RETURN
    END IF

    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT MATRIX CONVERSION 
!!! A and B already in full format
!!! forming tA and tB Transpose for case 31/41 
    !    call wallocate_1c(tsa,isa(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsa,isa(n+1)-1,infoloc)
    !    call wallocate_1i(tisa,n+1,infoloc)
    !    call wallocate_1c(tsb,isb(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsb,isb(n+1)-1,infoloc)
    !    call wallocate_1i(tisb,n+1,infoloc)
    !    if (infoloc/=0) then
    !       info=-1
    !       return
    !    end if
    !
    !    call ccsr_htranspose(n,sa,isa,jsa,tsa,tisa,tjsa)
    !    call ccsr_htranspose(n,sb,isb,jsb,tsb,tisb,tjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work,N,2*M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!!!!!!!!!!!!!!! find fpm(22) - # total number of contour points
    fpm(22)=fpm(8) ! default


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(22),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if

!!!!!!!!!!!!!!!!! Set up for Az matrix

    call wallocate_1i(isaz,n+1,infoloc)
    call wallocate_2c(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    call caddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2c(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    call caddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)


    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=nfact  ! Rq: same factorization for (normal+transpose)
    MTYPE=13      ! complex and unsymmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
    if (fpm(64)==1) then
       do i=1,64
          if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
       enddo
    endif
!!!!!!!!!!!!
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0!0 !0- no output, 1- output
    PHASE=11 ! symbolic factorization (do it only once)
!!!!!!!! single precision pardiso (MKL only)
    IPARM(28)=1 ! pardiso single precision
!!!!!!!!!!!!!!!
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ijob=-1 ! initialization
    do while (ijob/=0) 
       call cfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

       select case(ijob)
       case(10) !! factorize (zeB-A)

          if (fpm(10)==1) then
             id=id+1
             if (id==nfact+1) then
                id=1
                fact=.false.
             endif
          endif

          if (fact) then
             opt=3
             call caddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz


             if (PHASE==11) then
                PHASE=12
             else
                PHASE=22
             endif

             !            PHASE=22 ! numerical fact only
             MNUM=id
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)


             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          IPARM(12)=0 ! normal solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)  
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          IPARM(12)=1 ! transpose conjugate solve
          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)  
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ccsrmm('F','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))


       case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call ccsrmm('F','C',N,N,fpm(35),ONEC,sa,isa,jsa,X(1,fpm(34)),ZEROC,work(1,fpm(34)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ccsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))


       case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call ccsrmm('F','C',N,N,fpm(35),ONEC,sb,isb,jsb,X(1,fpm(34)),ZEROC,work(1,fpm(34)))


       end select
    end do


!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    do MNUM=1,nfact
       call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,MNUM),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
       if (infoloc/=0) then
          info=-2
          return
       end if
    end do


    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)

    call wdeallocate_2c(caux)

    call wdeallocate_2c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    !    call wdeallocate_1c(tsa)
    !    call wdeallocate_1i(tisa)
    !    call wdeallocate_1i(tjsa)
    !    call wdeallocate_1c(tsb)
    !    call wdeallocate_1i(tisb)
    !    call wdeallocate_1i(tjsb)



  end subroutine cfeast_gcsrgvx







  subroutine cfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B COMPLEX SYMMETRIC :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    complex :: Emid
    real :: r
    integer :: M0
    complex,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,infoloc,i,s
    complex :: Ze
    complex, dimension(:,:),pointer ::workc,work,zAq,zSq,caux
    complex,dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
    complex,dimension(:),pointer :: ssa,ssb
    integer,dimension(:), pointer :: sisa,sjsa,sisb,sjsb
    integer :: opt,nnza,nnzb,nnz
!!!!!for pardiso
    integer(8),dimension(64) :: pt
    integer,dimension(64) :: iparm
    integer :: mtype
    integer :: MAXFCT,MNUM,PHASE,MSGLVL
    integer :: idum


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------

    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_SCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO CSR-UPPER for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)

       !!<<<
       call wallocate_1c(ssa,1,infoloc) ! dummy
       call wallocate_1i(sjsa,1,infoloc) !dummy
       call wallocate_1c(ssb,1,infoloc) !dummy
       call wallocate_1i(sjsb,1,infoloc) !dummy
       if (infoloc/=0) then
          info=-1
          return
       end if
       !!>>>
       opt=1
       call ccsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call ccsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 
       !!<<<
       call wdeallocate_1c(ssa)
       call wdeallocate_1i(sjsa)
       call wdeallocate_1c(ssb)
       call wdeallocate_1i(sjsb)
       !!>>>

       call wallocate_1c(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1c(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call ccsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call ccsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       ! nnza=isa(n+1)-1
       ! ssa => sa(1:nnza)
       ! sisa => isa(1:n+1)
       ! sjsa => jsa(1:nnza)

       ! nnzb=isb(n+1)-1
       ! ssb =>  sb(1:nnzb)
       ! sisb => isb(1:n+1)
       ! sjsb => jsb(1:nnzb)


    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1c(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1c(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call ccsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call ccsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work,N,M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!!!!!!!!!!!!!!! find fpm(22) - # total number of contour points
    fpm(22)=fpm(8) ! default


!!!!!!!!!!!!!!!!!! Factorizations Set-up  
    fact=.true.
    nfact=0
    !! nfact is local (number of total factorization by contour points)
    do i=rank+1,fpm(22),nb_procs
       nfact=nfact+1
    end do

    if (fpm(10)==1) then
       id=0
    else
       id=1
    end if



!!!!!!!!!!!!!!!!! Set up for Az matrix

    call wallocate_1i(isaz,n+1,infoloc)
    call wallocate_2c(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='U').or.(UPLO=='u')) then
       call caddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    else
       call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get isaz
    end if
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2c(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2c(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='U').or.(UPLO=='u')) then
       call caddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
    else
       call caddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get jsaz
    end if
!!!!!!!!!!!!!!!
    call wallocate_2c(caux,N,M0,infoloc)


    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  pardiso initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MAXFCT=nfact ! Rq: same factorization for (normal+transpose)
    MTYPE=6      ! complex and symmetric 
    call pardisoinit(PT,MTYPE,IPARM)
!!!!!!!!!!!!
    if (fpm(64)==1) then
       do i=1,64
          if (fpm(64+i)/=-9999) iparm(i)=fpm(64+i)
       enddo
    endif
!!!!!!!!!!!!
    IPARM(6)=1 ! solution and rhs are input/output, attention caux is always used
    MSGLVL=0!0 !0- no output, 1- output
    PHASE=11 ! symbolic factorization (do it only once)
!!!!!!!! single precision pardiso (MKL only)
    IPARM(28)=1 ! pardiso single precision
!!!!!!!!!!!!!!!
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


    ijob=-1 ! initialization 
    do while (ijob/=0)
       call cfeast_srcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

       select case(ijob)
       case(10) !! factorize (zeB-A)

          if (fpm(10)==1) then
             id=id+1
             if (id==nfact+1) then
                id=1
                fact=.false.
             endif
          endif

          if (fact) then
             opt=3
             if ((UPLO=='U').or.(UPLO=='u')) then
                call caddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call caddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
             end if

             if (PHASE==11) then
                PHASE=12
             else
                PHASE=22
             endif
             !             PHASE=22 ! numerical fact only
             MNUM=id
             call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)  
             if (infoloc/=0) then
                info=-2
                return
             end if
          endif ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          PHASE=33 ! solve
          MNUM=id
          call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,id),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc) 
          if (infoloc/=0) then
             info=-2
             return
          end if

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
          if ((UPLO=='U').or.(UPLO=='u')) then
             call ccsrmm('U','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call ccsrmm('U','N',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          endif

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
          if ((UPLO=='U').or.(UPLO=='u')) then
             call ccsrmm('U','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call ccsrmm('U','N',N,N,fpm(25),ONEC,ssb,sisb,sjsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          end if

       end select
    end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!! release memory
!!!!!!!!!!!!!!!!!!!!!!!!
    PHASE=-1 

    do MNUM=1,nfact
       call PARDISO(PT,MAXFCT,MNUM,MTYPE,PHASE,n,saz(1,MNUM),isaz,jsaz,idum,fpm(23),IPARM,MSGLVL,workc,caux,infoloc)
       if (infoloc/=0) then
          info=-2
          return
       end if
    end do

    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(caux)


    call wdeallocate_2c(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1c(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1c(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif

  end subroutine cfeast_scsrgvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine sfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !===================================================================== 
    implicit none
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call sfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_scsrgv




  subroutine cfeast_hcsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call cfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine cfeast_hcsrgv



  subroutine sfeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL REAL:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        REAL SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    integer :: N
    real,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call sfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_gcsrgv


!!$

  subroutine cfeast_gcsrgv(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    integer :: N
    complex,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call cfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_gcsrgv


!!$

  subroutine cfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX SYMMETRIC:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !                           Rq: It is a dummy parameter here since matrix is general
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
    !  sb         (input)        COMPLEX SINGLE PRECISION (isb(N+1)-1):  Matrix B- CSR format 
    !  isb        (input)        INTEGER(N+1): CSR row array of Matrix B
    !  jsb        (input)        INTEGER(isb(N+1)-1): CSR column array of Matrix B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution 
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa,sb
    integer,dimension(*) :: isa,jsa,isb,jsb
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call cfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_scsrgv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!$

  subroutine sfeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0e0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call sfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_scsrevx





  subroutine cfeast_hcsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_hcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0e0,0.0e0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call cfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_hcsrevx



  subroutine sfeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL REAL :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    integer :: N
    real,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0e0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call sfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine sfeast_gcsrevx




  subroutine cfeast_gcsrevx(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    integer :: N
    complex,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_gcsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0e0,0.0e0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call cfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine cfeast_gcsrevx



  subroutine cfeast_scsrevx(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX SYMMETRIC:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !                           Rq: It is a dummy parameter here since matrix is general
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution 
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX SINGLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
    complex,dimension(*):: Zne,Wne
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0e0,0.0e0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call cfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine cfeast_scsrevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine sfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    ! ====================================================================
    implicit none
    character(len=1) :: UPLO
    integer :: N
    real,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    real,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_scsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 
    real,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0e0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call sfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_scsrev








  subroutine cfeast_hcsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output) INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    real :: epsout 
    integer :: loop
    real :: Emin,Emax
    integer :: M0
    real,dimension(*)  :: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_hcsrgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 
    complex,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0e0,0.0e0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call cfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_hcsrev




  subroutine sfeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL REAL :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        REAL SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    integer :: N,LDA
    real,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 
    real,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0e0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1
    call sfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_gcsrev




  subroutine cfeast_gcsrev(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX :: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    integer :: N
    complex,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_gcsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne
    complex,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i


    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)
    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0e0,0.0e0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call cfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_gcsrev




  subroutine cfeast_scsrev(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST SPARSE INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX SYMMETRIC:: SPARSE FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !                           Rq: It is a dummy parameter here since matrix is general
    !
    !  N          (input)        INTEGER: Size system
    !  sa         (input)        COMPLEX SINGLE PRECISION (isa(N+1)-1):  Matrix A- CSR format 
    !  isa        (input)        INTEGER(N+1): CSR row array of Matrix A
    !  jsa        (input)        INTEGER(isa(N+1)-1): CSR column array of Matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm (input/output)        INTEGER(*) : FEAST parameters
    !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
    !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
    !  Emid       (input/output) COMPLEX SINGLE PRECISION: middle of the search interval
    !  r          (input/output) REAL SINGLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1 
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX SINGLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution 
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    character(len=1) :: UPLO
    integer :: N
    complex,dimension(*) :: sa
    integer,dimension(*) :: isa,jsa
    integer,dimension(*) :: fpm
    complex :: Emid
    real :: r 
    real :: epsout
    integer :: loop
    integer :: M0
    complex,dimension(*):: E
    complex,dimension(N,*):: X
    integer :: mode
    real,dimension(*)    :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: cfeast_scsrgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 
    complex,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i


    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0e0,0.0e0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call cfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_scsrev







