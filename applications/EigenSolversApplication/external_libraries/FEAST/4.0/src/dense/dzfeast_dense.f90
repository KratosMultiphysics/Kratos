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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED DENSE INTERFACES !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

! List of routines:
!-------------------

!{D,Z}FEAST_{SY,HE,GE}{EV,GV}{X}


!111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE)
! dfeast_sygvx
! zfeast_hegvx
! dfeast_gegvx
! zfeast_gegvx
! zfeast_sygvx

!222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
! dfeast_sygv
! zfeast_hegv
! dfeast_gegv
! zfeast_gegv
! zfeast_sygv

!333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
! dfeast_syevx
! zfeast_heevx
! dfeast_geevx
! zfeast_geevx
! zfeast_syevx

!44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
! dfeast_syev
! zfeast_heev
! dfeast_geev
! zfeast_geev
! zfeast_syev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  !  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork
  double precision, dimension(:,:),allocatable ::work,Aq,Bq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  real,dimension(:,:),allocatable :: sA,sB ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121125 ! code name
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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_SYGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Format "conversion"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='L'
  IF ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  IF ((UPLO=='U').or.(UPLO=='u')) UPLO2='U'


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



!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(sA(N,N))
     call DLAG2S(N,N,A,LDA,sA,N,infoloc1)
     if (LDB/=-1) then
        allocate(sB(N,N))
        call DLAG2S(N,N,B,LDB,sB,N,infoloc1)
     endif
     allocate(Ac(N,N,nfact))
  else 
     allocate(Az(N,N,nfact))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up LAPACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(ipivloc(N,nfact))  
  if (fpm(42)==1)  allocate(cwork(N,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(Aq(M0,M0))
  allocate(Bq(M0,M0))
  allocate(work(N,M0))
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization 
  do while (ijob/=0)
     call dfeast_srcix(ijob,N,Ze,work,zwork,Aq,Bq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        id=fpm(33) !! id of factorization (for fpm(10) flag) 


        if (fpm(42)==1) then ! single precision
           if (LDB==-1) then !! standard
              Ac(1:N,1:N,id)=-sA(1:N,1:N)*CONE
              do i=1,N
                 Ac(i,i,id)=Ac(i,i,id)+cmplx(Ze)
              enddo
           else        
              Ac(1:N,1:N,id)=cmplx(Ze)*sB(1:N,1:N)-sA(1:N,1:N)*CONE 
           end if

           call CSYTRF(UPLO2,N,Ac(1,1,id),N,IPIVloc(1,id),cwork,N*M0,INFOloc2)

        else ! double precision

           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)*ZONE
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else        
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ZONE 
           end if

           call ZSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),zwork,N*M0,INFOloc2)

        end if


        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
        id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CSYTRS( UPLO2, N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
           call ZSYTRS( UPLO2, N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
        end if
        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call DSYMM ('L', UPLO2, N, fpm(25), DONE, A, LDA, X(1,fpm(24)), N, DZERO,work(1,fpm(24)), N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call DLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call DSYMM ('L', UPLO2, N, fpm(25), DONE, B, LDB, X(1,fpm(24)), N, DZERO,work(1,fpm(24)), N)
        end if

     end select
  end do



  deallocate(Aq)
  deallocate(Bq)
  deallocate(work)
  deallocate(zwork)
  deallocate(ipivloc)
  if (fpm(42)==1) then
     deallocate(cwork)
     deallocate(sA)
      if (LDB/=-1)  deallocate(sB)
     deallocate(Ac)
  else
     deallocate(Az)
  end if


end subroutine dfeast_sygvx





subroutine zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input/output) COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A
  !                            exit: full format if lower or upper 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input/output) COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B
  !                            exit: full format if lower or upper        
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  !  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,zAq,zBq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az 
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:),allocatable :: cA,cB ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
 character(len=1) :: UPLO2

  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141255 ! code name
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

  INFO = 0
  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF ( N<=0 ) THEN
     INFO = -102
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_HEGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Format CONVERSION !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Lower/upper format to full !!!!!!!!!!!!!!!!!
!!!!!!!! original A and B are modified !!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if ((UPLO=='L').or.(UPLO=='l')) then
     do i=1,N-1
        s=N-i
        call ZCOPY(s,A(i+1,i),1,A(i,i+1),LDA)
        call ZLACGV(s, A(i,i+1), LDA )
     enddo

     if (LDB/=-1) then
        do i=1,N-1
           s=N-i
           call ZCOPY(s,B(i+1,i),1,B(i,i+1),LDB)
           call ZLACGV(s, B(i,i+1), LDB )
        enddo
     end if

!!!!
  elseif ((UPLO=='U').or.(UPLO=='u')) then
!!!!
     do i=1,N-1
        s=N-i
        call ZCOPY(s,A(i,i+1),LDA,A(i+1,i),1)
        call ZLACGV(s, A(i+1,i), 1 )
     enddo

     if (LDB/=-1) then
        do i=1,N-1
           s=N-i
           call ZCOPY(s,B(i,i+1),LDB,B(i+1,i),1)
           call ZLACGV(s, B(i+1,i), 1 )
        enddo
     end if

  endif
  
  !! for case (30,40)
 UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='L'
  IF ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  IF ((UPLO=='U').or.(UPLO=='u')) UPLO2='U'


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

!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(cA(N,N))
     call ZLAG2C(N,N,A,LDA,cA,N,infoloc1)
     if (LDB/=-1) then
        allocate(cB(N,N))
        call ZLAG2C(N,N,B,LDB,cB,N,infoloc1)
     endif
     allocate(Ac(N,N,nfact))
  else
     allocate(Az(N,N,nfact))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up LAPACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(ipivloc(N,nfact))  
  if (fpm(42)==1)  allocate(cwork(N,M0))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(zAq(M0,M0))
  allocate(zBq(M0,M0))
  allocate(work(N,M0))
  allocate(zwork(N,M0))





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_hrcix(ijob,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag) 

        if (fpm(42)==1) then ! single precision
           if (LDB==-1) then !! standard
              Ac(1:N,1:N,id)=-cA(1:N,1:N)
              do i=1,N
                 Ac(i,i,id)=Ac(i,i,id)+cmplx(Ze)
              enddo
           else        
              Ac(1:N,1:N,id)=cmplx(Ze)*cB(1:N,1:N)-cA(1:N,1:N)
           end if

           call CGETRF(N,N,Ac(1,1,id),N,IPIVloc(1,id),INFOloc2)

        else ! double precision

           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else        
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N) 
           end if

           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc2)
        end if


        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CGETRS( 'N', N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else

           
           call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
          

        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CGETRS( 'C', N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
            call ZGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
           

        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

        !        if ((UPLO=='F').or.(UPLO=='f')) then
        !           call ZGEMM('N','N',N,fpm(25),N,ZONE,A,LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
        !        else
        call ZHEMM ('L', UPLO2, N, fpm(25), ZONE, A, LDA, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)
        !        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           !  if ((UPLO=='F').or.(UPLO=='f')) then
           !     call ZGEMM('N','N',N,fpm(25),N,ZONE,B,LDB,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
           !  else
           call ZHEMM ('L', UPLO2, N, fpm(25), ZONE, B, LDB, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)
           !  endif
        end if


     end select
  end do

  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)
  deallocate(ipivloc)
  if (fpm(42)==1) then
     deallocate(cwork)
     deallocate(cA)
      if (LDB/=-1)  deallocate(cB)
     deallocate(Ac)
  else
     deallocate(Az)
  end if



end subroutine zfeast_hegvx







subroutine dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  !  include 'f90_noruntime_interface.fi'
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,zAq,zBq,zA,zB
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  real,dimension(:,:),allocatable :: sA,sB ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121355 ! code name
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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_GEGVX', -INFO+100 )
     RETURN
  END IF


  infoloc1=0
  infoloc2=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! :-( If DZGEMM  not available for case 30-40, we need to make a complex copy of the matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
#ifdef MKL
  ! everything fine no copy needed
#else
  ! (extra copy)
  allocate(zA(N,N))
  call ZLACP2( 'F', N, N,A , LDA, zA, N ) 
  if (LDB/=-1) then
     allocate(zB(N,N))
     call ZLACP2( 'F', N, N,B , LDB, zB, N ) 
  end if
#endif       

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


!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(sA(N,N))
     call DLAG2S(N,N,A,LDA,sA,N,infoloc1)
     if (LDB/=-1) then
        allocate(sB(N,N))
        call DLAG2S(N,N,B,LDB,sB,N,infoloc1)
     endif
     allocate(Ac(N,N,nfact))
  else 
     allocate(Az(N,N,nfact))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up LAPACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(ipivloc(N,nfact))  
  if (fpm(42)==1)  allocate(cwork(N,M0))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(zAq(M0,M0))
  allocate(zBq(M0,M0))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call dfeast_grcix(ijob,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag) 


        if (fpm(42)==1) then ! single precision

           if (LDB==-1) then !! standard
              Ac(1:N,1:N,id)=-sA(1:N,1:N)*CONE
              do i=1,N
                 Ac(i,i,id)=Ac(i,i,id)+cmplx(Ze)
              enddo
           else     
              Ac(1:N,1:N,id)=cmplx(Ze)*sB(1:N,1:N)-sA(1:N,1:N)*CONE
           end if

           call CGETRF(N,N,Ac(1,1,id),N,IPIVloc(1,id),INFOloc2)    

        else ! double precision

           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)*ZONE
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else     
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ZONE
           end if
           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc2)     

        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CGETRS( 'N', N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else 
           call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

        id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CGETRS( 'C', N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else 
           call ZGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MKL
        ! With MKL-blas- one should use 
        call DZGEMM('N','N',N,fpm(25),N,ZONE,A(1,1),LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#else
        ! not optimal (extra copy)       
        call ZGEMM('N','N',N,fpm(25),N,ZONE,zA(1,1),N,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#endif       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(31) !! perform multiplication A^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MKL
        ! With MKL-blas- one should use
        call DZGEMM('T','N',N,fpm(35),N,ZONE,A(1,1),LDA,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#else
        ! not optimal (extra copy)
        call ZGEMM('C','N',N,fpm(35),N,ZONE,zA(1,1),N,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#endif      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
#ifdef MKL
           ! With MKL-blas- one should use
           call DZGEMM('N','N',N,fpm(25),N,ZONE,B(1,1),LDB,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#else
           !  ! not optimal (extra copy)
           call ZGEMM('N','N',N,fpm(25),N,ZONE,zB(1,1),N,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#endif            
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(41)!! perform multiplication B^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
        else
#ifdef MKL
           ! With MKL-blas- one should use
           call DZGEMM('T','N',N,fpm(35),N,ZONE,B(1,1),LDB,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#else
           call ZGEMM('C','N',N,fpm(35),N,ZONE,zB(1,1),N,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#endif
        end if

     end select
  end do


  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)
  deallocate(ipivloc)
  if (fpm(42)==1) then
     deallocate(cwork)
     deallocate(sA)
      if (LDB/=-1)  deallocate(sB)
     deallocate(Ac)
  else
     deallocate(Az)
  end if
  
#ifdef MKL
  ! no copy needed
#else
  ! (extra copy)
  deallocate(zA) 
  if (LDB/=-1) deallocate(zB)
#endif       


end subroutine dfeast_gegvx




subroutine zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,zAq,zBq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:),allocatable :: cA,cB ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD



  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141355 ! code name
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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_GEGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
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

!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(cA(N,N))
     call ZLAG2C(N,N,A,LDA,cA,N,infoloc1)
     if (LDB/=-1) then
        allocate(cB(N,N))
        call ZLAG2C(N,N,B,LDB,cB,N,infoloc1)
     endif
     allocate(Ac(N,N,nfact))
  else 
     allocate(Az(N,N,nfact))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up LAPACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(ipivloc(N,nfact))  
  if (fpm(42)==1)  allocate(cwork(N,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(zAq(M0,M0))
  allocate(zBq(M0,M0))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_grcix(ijob,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33) !! id of factorization (for fpm(10) flag) 


        if (fpm(42)==1) then ! single precision

           if (LDB==-1) then !! standard
              Ac(1:N,1:N,id)=-cA(1:N,1:N)
              do i=1,N
                 Ac(i,i,id)=Ac(i,i,id)+cmplx(Ze)
              enddo
           else     
              Ac(1:N,1:N,id)=cmplx(Ze)*cB(1:N,1:N)-cA(1:N,1:N)
           end if

           call CGETRF(N,N,Ac(1,1,id),N,IPIVloc(1,id),INFOloc2)    

        else ! double precision

           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else     
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)
           end if
           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc2)     

        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CGETRS( 'N', N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else 
           call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CGETRS( 'C', N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)

        else 
           call ZGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call ZGEMM('N','N',N,fpm(25),N,ZONE,A(1,1),LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !print *,fpm(34)
        !print *,A(1:N,1:N)
        !print *,X(1:N,fpm(34))
        call ZGEMM('C','N',N,fpm(35),N,ZONE,A(1,1),LDA,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
         !print *,work(1:N,fpm(34))
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
           call ZGEMM('N','N',N,fpm(25),N,ZONE,B(1,1),LDB,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
        else
           call ZGEMM('C','N',N,fpm(35),N,ZONE,B(1,1),LDB,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
        end if

     end select
  end do


  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)
  deallocate(ipivloc)
  if (fpm(42)==1) then
     deallocate(cwork)
     deallocate(cA)
      if (LDB/=-1)  deallocate(cB)
     deallocate(Ac)
  else
     deallocate(Az)
  end if



end subroutine zfeast_gegvx







subroutine zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B COMPLEX SYMMETRIC :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution
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
  ! James Kestyn, Eric Polizzi 2015
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,work,zAq,zBq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:),allocatable :: cA,cB ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO)  
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141125 ! code name
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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_SYGVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Format "conversion"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='L'
  IF ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  IF ((UPLO=='U').or.(UPLO=='u')) UPLO2='U'


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


!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(cA(N,N))
     call ZLAG2C(N,N,A,LDA,cA,N,infoloc1)
     if (LDB/=-1) then
        allocate(cB(N,N))
        call ZLAG2C(N,N,B,LDB,cB,N,infoloc1)
     endif
     allocate(Ac(N,N,nfact))
  else 
     allocate(Az(N,N,nfact))
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up LAPACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(ipivloc(N,nfact))  
  if (fpm(42)==1)  allocate(cwork(N,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Set up FEAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(zAq(M0,M0))
  allocate(zBq(M0,M0))
  allocate(work(N,M0))
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_srcix(ijob,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! factorize (zeB-A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag) 

        if (fpm(42)==1) then ! single precision
           if (LDB==-1) then !! standard
              Ac(1:N,1:N,id)=-cA(1:N,1:N)
              do i=1,N
                 Ac(i,i,id)=Ac(i,i,id)+cmplx(Ze)
              enddo
           else        
              Ac(1:N,1:N,id)=cmplx(Ze)*cB(1:N,1:N)-cA(1:N,1:N) 
           end if

           call CSYTRF(UPLO2,N,Ac(1,1,id),N,IPIVloc(1,id),cwork,N*M0,INFOloc2)

        else ! double precision

           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else        
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)
           end if

           call ZSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),zwork,N*M0,INFOloc2)

        end if


        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     case(11) !!solve the linear system (ZeB-A)x=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CSYTRS( UPLO2, N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
           call ZSYTRS( UPLO2, N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2 )
        end if
        if (infoloc1/=0) then
           info=-1
           return
        end if
        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

        call ZSYMM ('L', UPLO2, N, fpm(25), ZONE, A, LDA, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call ZSYMM ('L', UPLO2, N, fpm(25), ZONE, B, LDB, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)
        end if

     end select
  end do


  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)
  deallocate(ipivloc)
  if (fpm(42)==1) then
     deallocate(cwork)
     deallocate(cA)
      if (LDB/=-1)  deallocate(cB)
     deallocate(Ac)
  else
     deallocate(Az)
  end if

end subroutine zfeast_sygvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine dfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: dfeast_sygvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=121124
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine dfeast_sygv




subroutine zfeast_hegv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: zfeast_hegvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=141254
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

end subroutine zfeast_hegv



subroutine dfeast_gegv(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: dfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  fpm(30)=121354
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine dfeast_gegv




subroutine zfeast_gegv(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: zfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  fpm(30)=141354
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_gegv




subroutine zfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
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
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: zfeast_sygvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  fpm(30)=141124
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_sygv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine dfeast_syevx(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  character(len=1) :: UPLO
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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
  double precision,dimension(1) :: B ! dummy
  integer :: LDB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_sygvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  fpm(30)=121123
  call feastdefault(fpm,info)
  if (info/=0) return
  call dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine dfeast_syevx





subroutine zfeast_heevx(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_hegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  fpm(30)=141253
  call feastdefault(fpm,info)
  if (info/=0) return
  call zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine zfeast_heevx



subroutine dfeast_geevx(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: dfeast_gegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  double precision,dimension(1) :: B ! dummy
  integer :: LDB

  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  fpm(30)=121353
  call feastdefault(fpm,info)
  if (info/=0) return
  call dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine dfeast_geevx




subroutine zfeast_geevx(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_gegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  fpm(30)=141353
  call feastdefault(fpm,info)
  if (info/=0) return
  call zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine zfeast_geevx



subroutine zfeast_syevx(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
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
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_sygvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  fpm(30)=141123
  call feastdefault(fpm,info)
  if (info/=0) return
  call zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine zfeast_syevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine dfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  ! ====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: dfeast_sygvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  double precision,dimension(1) :: B ! dummy
  integer :: LDB

  fpm(30)=121122
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call dfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine dfeast_syev








subroutine zfeast_heev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  ! ====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_hegvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

  fpm(30)=141252
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_heev




subroutine dfeast_geev(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: dfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  double precision,dimension(1) :: B ! dummy
  integer :: LDB

  fpm(30)=121352
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call dfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine dfeast_geev




subroutine zfeast_geev(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

  fpm(30)=141352
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_geev




subroutine zfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
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
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: zfeast_sygvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB

  fpm(30)=141122
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call zfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_syev







