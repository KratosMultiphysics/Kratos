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
  
  include 'dzlsprim.f90' !! Sparse primitives





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED SPARSE INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

  ! List of routines:
  !-------------------

  !{D,Z}FEAST_{SCSR,HCSR,GCSR}{EV,GV}


  !111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE)
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
    !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,caux
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
    double precision,dimension(:),pointer :: ssa,ssb
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
       CALL XERBLA( 'DFEAST_SCSRGV', -INFO+100 )
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
       call wallocate_1d(ssa,1,infoloc) ! dummy
       call wallocate_1i(sjsa,1,infoloc) !dummy
       call wallocate_1d(ssb,1,infoloc) !dummy
       call wallocate_1i(sjsb,1,infoloc) !dummy
       if (infoloc/=0) then
          info=-1
          return
       end if
       !!>>>
       opt=1
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call dcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 
       !!<<<
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sjsa)
       call wdeallocate_1d(ssb)
       call wdeallocate_1i(sjsb)
       !!>>>

       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1d(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

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
       call wallocate_1d(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1d(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call dcsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call dcsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
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
    call wallocate_2z(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='U').or.(UPLO=='u')) then
       call zdaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz,isaz,jsaz) !! get isaz
    else
       call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get isaz
    end if

    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2z(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='U').or.(UPLO=='u')) then
       call zdaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz,isaz,jsaz) !! get jsaz
    else
       call zdaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz,isaz,jsaz) !! get jsaz
    end if
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)


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
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


    ijob=-1 ! initialization 
    do while (ijob/=0)
       call dfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

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
                call zdaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call zdaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
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
             call dcsrmm('U','N',N,N,fpm(25),DONE,sa,isa,jsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          else
             call dcsrmm('U','N',N,N,fpm(25),DONE,ssa,sisa,sjsa,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          end if

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='U').or.(UPLO=='u')) then
             call dcsrmm('U','N',N,N,fpm(25),DONE,sb,isb,jsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))
          else
             call dcsrmm('U','N',N,N,fpm(25),DONE,ssb,sisb,sjsb,X(1,fpm(24)),DZERO,work(1,fpm(24)))
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

    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)


    call wdeallocate_2z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1d(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1d(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq,caux
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! full csr format
    complex(kind=kind(1.0d0)),dimension(:),pointer :: ssa,ssb
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
       CALL XERBLA( 'ZFEAST_HCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FORMAT CONVERSION TO full CSR for PARDISO !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Remark: everything should be csr since we work with zS-H which is actually unsymmetric 


    if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr already

       !       nnza=isa(n+1)-1
       !       ssa => sa(1:nnza)
       !       sisa => isa(1:n+1)
       !       sjsa => jsa(1:nnza)

       !       nnzb=isb(n+1)-1
       !       ssb => sb(1:nnzb)
       !       sisb => isb(1:n+1)
       !       sjsb => jsb(1:nnzb)

    else !! upper-csr or lower-csr to full csr

       nnza=2*(isa(n+1)-1) ! may be overestimated
       nnzb=2*(isb(n+1)-1)
       call wallocate_1z(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1z(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call zhcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
       call zhcsr_uplo_to_csr(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
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
    call wallocate_2z(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='F').or.(UPLO=='f')) then
       call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    else
       call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get isaz
    end if
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2z(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='F').or.(UPLO=='f')) then
       call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
    else
       call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get jsaz
    endif
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)

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
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ijob=-1 ! initialization
    do while (ijob/=0) 
       call zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

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
                call zaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
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
             call zhcsrmm('F','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call zhcsrmm('F','N',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          end if
       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='F').or.(UPLO=='f')) then
             call zhcsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call zhcsrmm('F','N',N,N,fpm(25),ONEC,ssb,sisb,sjsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
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



    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)


    call wdeallocate_2z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    if ((UPLO/='F').and.(UPLO/='f')) then
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif


  end subroutine zfeast_hcsrgvx



!!$



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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq,caux
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! transpose matrix
    !double precision,dimension(:),pointer :: tsa,tsb
    !integer,dimension(:), pointer :: tisa,tjsa,tisb,tjsb
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
       CALL XERBLA( 'DFEAST_GCSRGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT MATRIX CONVERSION 
!!! A and B already in full format
!!! forming tA and tB Transpose for case 31/41 
    !    call wallocate_1d(tsa,isa(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsa,isa(n+1)-1,infoloc)
    !    call wallocate_1i(tisa,n+1,infoloc)
    !    call wallocate_1d(tsb,isb(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsb,isb(n+1)-1,infoloc)
    !    call wallocate_1i(tisb,n+1,infoloc)
    !    if (infoloc/=0) then
    !       info=-1
    !       return
    !    end if
    !    call dcsr_transpose(n,sa,isa,jsa,tsa,tisa,tjsa)
    !    call dcsr_transpose(n,sb,isb,jsb,tsb,tisb,tjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,2*M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

!!!!!!!!!!!!!!!! find fpm(22) - # total number of contour points
    fpm(22)=fpm(8) ! default
    if (fpm(29)==1) then! default for generated contour -- half-contour capability case 
       if ((aimag(Emid)==0.0d0).and.((2*(fpm(8)/2))==fpm(8))) fpm(22)=fpm(8)/2  
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
    call wallocate_2z(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    call zdaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz,isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2z(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    call zdaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz,isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)


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
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ijob=-1 ! initialization
    do while (ijob/=0) 
       call dfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

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
             call zdaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz


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

          call dzcsrmm('F','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))

       case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call dzcsrmm('F','C',N,N,fpm(35),ONEC,sa,isa,jsa,X(1,fpm(34)),ZEROC,work(1,fpm(34)))


       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call dzcsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))


       case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call dzcsrmm('F','C',N,N,fpm(35),ONEC,sb,isb,jsb,X(1,fpm(34)),ZEROC,work(1,fpm(34)))


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


    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)

    call wdeallocate_2z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    !    call wdeallocate_1d(tsa)
    !    call wdeallocate_1i(tisa)
    !    call wdeallocate_1i(tjsa)
    !    call wdeallocate_1d(tsb)
    !    call wdeallocate_1i(tisb)
    !    call wdeallocate_1i(tjsb)

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq,caux
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!! transpose conjugate matrix
    !complex(kind=(kind(1.0d0))),dimension(:),pointer :: tsa,tsb
    !integer,dimension(:), pointer :: tisa,tjsa,tisb,tjsb
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
       CALL XERBLA( 'ZFEAST_GCSRGV', -INFO+100 )
       RETURN
    END IF

    infoloc=0
!!!!!!!!!!!!!!!!!!!!!!!!!! FORMAT MATRIX CONVERSION 
!!! A and B already in full format
!!! forming tA and tB Transpose for case 31/41 
    !    call wallocate_1z(tsa,isa(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsa,isa(n+1)-1,infoloc)
    !    call wallocate_1i(tisa,n+1,infoloc)
    !    call wallocate_1z(tsb,isb(n+1)-1,infoloc) 
    !    call wallocate_1i(tjsb,isb(n+1)-1,infoloc)
    !    call wallocate_1i(tisb,n+1,infoloc)
    !    if (infoloc/=0) then
    !       info=-1
    !       return
    !    end if
    !
    !    call zcsr_htranspose(n,sa,isa,jsa,tsa,tisa,tjsa)
    !    call zcsr_htranspose(n,sb,isb,jsb,tsb,tisb,tjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,2*M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
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
    call wallocate_2z(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2z(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)


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
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ijob=-1 ! initialization
    do while (ijob/=0) 
       call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

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
             call zaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz


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

          call zcsrmm('F','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))


       case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call zcsrmm('F','C',N,N,fpm(35),ONEC,sa,isa,jsa,X(1,fpm(34)),ZEROC,work(1,fpm(34)))

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call zcsrmm('F','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))


       case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call zcsrmm('F','C',N,N,fpm(35),ONEC,sb,isb,jsb,X(1,fpm(34)),ZEROC,work(1,fpm(34)))


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


    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)

    call wdeallocate_2z(caux)

    call wdeallocate_2z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)

    !    call wdeallocate_1z(tsa)
    !    call wdeallocate_1i(tisa)
    !    call wdeallocate_1i(tjsa)
    !    call wdeallocate_1z(tsb)
    !    call wdeallocate_1i(tisb)
    !    call wdeallocate_1i(tjsb)



  end subroutine zfeast_gcsrgvx







  subroutine zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
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
    !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
    !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
    !                            comments: (input)  Emid,r dummy variables
    !                                      (output) Emid,r calculated/estimated from custom contour
    !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
    !                                             (output) Emid,r unchanged if fpm(29)=1
    !  M0         (input/output) INTEGER: Size subspace
    !  E          (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,work,zAq,zSq,caux
    complex(kind=(kind(1.0d0))),dimension(:,:),pointer :: saz
    integer,dimension(:),pointer :: isaz,jsaz
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
!!!!! csr-upper format
    complex(kind=(kind(1.0d0))),dimension(:),pointer :: ssa,ssb
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
       CALL XERBLA( 'ZFEAST_SCSRGV', -INFO+100 )
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
       call wallocate_1z(ssa,1,infoloc) ! dummy
       call wallocate_1i(sjsa,1,infoloc) !dummy
       call wallocate_1z(ssb,1,infoloc) !dummy
       call wallocate_1i(sjsb,1,infoloc) !dummy
       if (infoloc/=0) then
          info=-1
          return
       end if
       !!>>>
       opt=1
       call zcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       nnza=sisa(n+1)-1 
       call zcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)
       nnzb=sisb(n+1)-1 
       !!<<<
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sjsa)
       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sjsb)
       !!>>>

       call wallocate_1z(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1z(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       opt=2
       call zcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
       call zcsr2csr_up(opt,N,sb,isb,jsb,ssb,sisb,sjsb)

!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 

       !       nnza=isa(n+1)-1
       !       ssa => sa(1:nnza)
       !       sisa => isa(1:n+1)
       !       sjsa => jsa(1:nnza)
       !
       !       nnzb=isb(n+1)-1
       !       ssb =>  sb(1:nnzb)
       !       sisb => isb(1:n+1)
       !       sjsb => jsb(1:nnzb)


    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr

       nnza=isa(n+1)-1
       nnzb=isb(n+1)-1
       call wallocate_1z(ssa,nnza,infoloc)
       call wallocate_1i(sjsa,nnza,infoloc)
       call wallocate_1z(ssb,nnzb,infoloc)
       call wallocate_1i(sjsb,nnzb,infoloc)
       call wallocate_1i(sisa,n+1,infoloc)
       call wallocate_1i(sisb,n+1,infoloc)
       if (infoloc/=0) then
          info=-1
          return
       end if

       call zcsr_transpose(N,sa,isa,jsa,ssa,sisa,sjsa)
       call zcsr_transpose(N,sb,isb,jsb,ssb,sisb,sjsb)

    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
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
    call wallocate_2z(saz,1,1,infoloc) ! dummy
    call wallocate_1i(jsaz,1,infoloc)! dummy
    !!>>>
    opt=1
    if ((UPLO=='U').or.(UPLO=='u')) then
       call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get isaz
    else
       call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get isaz
    end if

    nnz=isaz(n+1)-1
    !!<<<
    call wdeallocate_2z(saz)
    call wdeallocate_1i(jsaz)
    !!>>>
    call wallocate_2z(saz,nnz,nfact,infoloc)
    call wallocate_1i(jsaz,nnz,infoloc)

    opt=2
    if ((UPLO=='U').or.(UPLO=='u')) then
       call zaddcsr(N,opt,ONEC,sa,isa,jsa,ONEC,sb,isb,jsb,saz(1,1),isaz,jsaz) !! get jsaz
    else
       call zaddcsr(N,opt,ONEC,ssa,sisa,sjsa,ONEC,ssb,sisb,sjsb,saz(1,1),isaz,jsaz) !! get jsaz
    end if
!!!!!!!!!!!!!!!
    call wallocate_2z(caux,N,M0,infoloc)


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
    !    call PARDISO(PT,MAXFCT,1,MTYPE,PHASE,n,saz(1,1),isaz,jsaz,idum,M0,IPARM,MSGLVL,workc,caux,infoloc)
!!!!!!!!!!


    ijob=-1 ! initialization 
    do while (ijob/=0)
       call zfeast_srcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

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
                call zaddcsr(N,opt,-ONEC,sa,isa,jsa,Ze,sb,isb,jsb,saz(1,id),isaz,jsaz) !! get saz
             else
                call zaddcsr(N,opt,-ONEC,ssa,sisa,sjsa,Ze,ssb,sisb,sjsb,saz(1,id),isaz,jsaz) !! get saz
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
             call zcsrmm('U','N',N,N,fpm(25),ONEC,sa,isa,jsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call zcsrmm('U','N',N,N,fpm(25),ONEC,ssa,sisa,sjsa,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          end if

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if ((UPLO=='U').or.(UPLO=='u')) then
             call zcsrmm('U','N',N,N,fpm(25),ONEC,sb,isb,jsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
          else
             call zcsrmm('U','N',N,N,fpm(25),ONEC,ssb,sisb,sjsb,X(1,fpm(24)),ZEROC,work(1,fpm(24)))
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

    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(caux)


    call wdeallocate_2z(saz)
    call wdeallocate_1i(isaz)
    call wdeallocate_1i(jsaz)


    if ((UPLO/='U').and.(UPLO/='u')) then
       call wdeallocate_1z(ssa)
       call wdeallocate_1i(sisa)
       call wdeallocate_1i(sjsa)

       call wdeallocate_1z(ssb)
       call wdeallocate_1i(sisb)
       call wdeallocate_1i(sjsb)
    endif

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
    !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !===================================================================== 
    implicit none
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

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     call dfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
   
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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call zfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)


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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call dfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call zfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

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
    !                           Rq: It is a dummy parameter here since matrix is general
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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution 
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

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
    !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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
    double precision,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_scsrgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0d0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call dfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
    !           
    !=====================================================================
    ! Eric Polizzi 2009-2015
    !=====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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

    complex(kind=(kind(1.0d0))),dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0d0,0.0d0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call zfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    double precision,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0d0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call dfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))),dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0d0,0.0d0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call zfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


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
    !                           Rq: It is a dummy parameter here since matrix is general
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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution 
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
    !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
    !           
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))), dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0d0,0.0d0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


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
    !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
    double precision,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0d0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call dfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Eigenvectors-solution
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! Eric Polizzi 2009-2015
    ! ====================================================================
    implicit none
    include 'f90_noruntime_interface.fi'
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

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
    complex(kind=(kind(1.0d0))),dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0d0,0.0d0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call zfeast_hcsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
    double precision,dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=1.0d0
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1
    call dfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
    !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne
    complex(kind=(kind(1.0d0))),dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i


    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)
    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0d0,0.0d0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call zfeast_gcsrgvx(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

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
    !                           Rq: It is a dummy parameter here since matrix is general
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
    !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
    !                                                       On entry: subspace initial guess if fpm(5)=1 
    !                                                       On exit : Right Eigenvectors-solution 
    !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
    !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
    !                                                        if option fpm(6)=1 selected                           
    !  info       (output)       INTEGER: Error handling (0: successful exit)
    !=====================================================================
    ! James Kestyn, Eric Polizzi 2015
    ! ====================================================================
    implicit none
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

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
    complex(kind=(kind(1.0d0))), dimension(n) :: sb ! identity
    integer,dimension(n+1) :: isb
    integer,dimension(n) :: jsb
    integer :: i


    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    ! identity B matrix- option for standard eigenvalue problem
    do i=1,n
       sb(i)=(1.0d0,0.0d0)
       jsb(i)=i
       isb(i)=i
    enddo
    isb(n+1)=n+1

    call zfeast_scsrgvx(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine zfeast_scsrev







