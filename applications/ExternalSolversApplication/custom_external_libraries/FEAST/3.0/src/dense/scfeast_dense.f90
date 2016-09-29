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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED DENSE INTERFACES !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! SINGLE PRECISION VERSION

! List of routines:
!-------------------

!{S,C}FEAST_{SY,HE,GE}{EV,GV}


!111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE)
! sfeast_sygvx
! cfeast_hegvx
! sfeast_gegvx
! cfeast_gegvx
! cfeast_sygvx

!222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
! sfeast_sygv
! cfeast_hegv
! sfeast_gegv
! cfeast_gegv
! cfeast_sygv

!333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
! sfeast_syevx
! cfeast_heevx
! sfeast_geevx
! cfeast_geevx
! cfeast_syevx

!44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
! sfeast_syev
! cfeast_heev
! sfeast_geev
! cfeast_geev
! cfeast_syev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine sfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE  PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL SINGLE  PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL SINGLE  PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL SINGLE  PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  E          (output)       REAL SINGLE  PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) REAL SINGLE  PRECISION(N,M0) : 
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
  integer :: N,LDA,LDB
  real,dimension(LDA,*):: A
  real,dimension(LDB,*):: B
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
  complex, dimension(:,:),pointer ::workc
  real, dimension(:,:),pointer ::work,Aq,Sq
  complex, dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO)
  logical :: fact
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
     CALL XERBLA( 'SFEAST_SYGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

  call wallocate_2s(Aq,M0,M0,infoloc)
  call wallocate_2s(Sq,M0,M0,infoloc)
  call wallocate_2s(work,N,M0,infoloc)
  call wallocate_2c(workc,N,M0,infoloc)

  if (infoloc/=0) then
     info=-1
     return
  end if

!!! Format conversion
  UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='L'
  IF ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  IF ((UPLO=='U').or.(UPLO=='u')) UPLO2='U'

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3c(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)*ONEC
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else        
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ONEc 
           end if

           call CSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),workc,N*M0,INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
        endif ! fact true

     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        call CSYTRS( UPLO2, N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call SSYMM ('L', UPLO2, N, fpm(25), SONE, A, LDA, X(1,fpm(24)), N, SZERO,work(1,fpm(24)), N)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call SLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call SSYMM ('L', UPLO2, N, fpm(25), SONE, B, LDB, X(1,fpm(24)), N, SZERO,work(1,fpm(24)), N)
        end if

     end select
  end do



  call wdeallocate_2s(Aq)
  call wdeallocate_2s(Sq)
  call wdeallocate_2s(work)
  call wdeallocate_2c(workc)
  call wdeallocate_3c(Az)
  call wdeallocate_2i(ipivloc)

end subroutine sfeast_sygvx





subroutine cfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex,dimension(LDA,*):: A
  complex,dimension(LDB,*):: B
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
  complex, dimension(:,:),pointer ::work,workc,zAq,zSq
  complex, dimension(:,:,:),pointer :: Az
  complex, dimension(:),pointer ::ztmp
  integer, dimension(:,:),pointer ::ipivloc
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
  logical :: fact
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
     CALL XERBLA( 'CFEAST_HEGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3c(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  call wallocate_1c(ztmp,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
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
           if (LDB==-1) then !! standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

              if ((UPLO=='L').or.(UPLO=='l')) then
                 do i=1,N-1
                    s=N-i
                    Az(i:N,i,id)=-A(i:N,i)
                    Az(i,i,id)=Az(i,i,id)+Ze
                    call CCOPY(s,Az(i+1,i,id),1,Az(i,i+1,id),N)
                    call CLACGV(s, Az(i,i+1,id), N )
                 enddo
                 Az(N,N,id)=Ze-A(N,N)*ONEC

              elseif ((UPLO=='U').or.(UPLO=='u')) then
                 do i=1,N-1
                    s=N-i
                    Az(i,i:N,id)=-A(i,i:N)
                    Az(i,i,id)=Az(i,i,id)+Ze
                    call CCOPY(s,Az(i,i+1,id),N,Az(i+1,i,id),1)
                    call CLACGV(s, Az(i+1,i,id), 1 )
                 enddo
                 Az(N,N,id)=Ze-A(N,N)*ONEC

              else

                 Az(1:N,1:N,id)=-A(1:N,1:N)
                 do i=1,N
                    Az(i,i,id)=Az(i,i,id)+Ze
                 enddo

              end if


           else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
              if ((UPLO=='L').or.(UPLO=='l')) then
                 do i=1,N-1
                    s=N-i
                    Az(i:N,i,id)=Ze*B(i:N,i)
                    call CCOPY(s,Az(i+1,i,id),1,Az(i,i+1,id),N)
                    call CLACGV(s, Az(i,i+1,id), N )
                    call CSCAL(s,Ze/conjg(Ze),Az(i,i+1,id),N)
                    Az(i:N,i,id)=Az(i:N,i,id)-A(i:N,i)
                 enddo
                 Az(N,N,id)=Ze*B(N,N)-A(N,N)*ONEC
                 do i=1,N-1
                    s=N-i
                    call CCOPY(s,A(i+1,i),1,ztmp(1),1)
                    call CLACGV(s,ztmp(1),1)
                    call CAXPY(s,-ONEC,ztmp(1),1,Az(i,i+1,id),N)
                 enddo


              elseif ((UPLO=='U').or.(UPLO=='u')) then

                 do i=1,N-1
                    s=N-i
                    Az(i,i:N,id)=Ze*B(i,i:N)
                    call CCOPY(s,Az(i,i+1,id),N,Az(i+1,i,id),1)
                    call CLACGV(s, Az(i+1,i,id), 1 )
                    call CSCAL(s,Ze/conjg(Ze),Az(i+1,i,id),1)
                    Az(i,i:N,id)=Az(i,i:N,id)-A(i,i:N)
                 enddo
                 Az(N,N,id)=Ze*B(N,N)-A(N,N)*ONEC
                 do i=1,N-1
                    s=N-i
                    call CCOPY(s,A(i,i+1),N,ztmp(1),1)
                    call CLACGV(s,ztmp(1),1)
                    call CAXPY(s,-ONEC,ztmp(1),1,Az(i+1,i,id),1)
                 enddo
              else ! full 
                 Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)
              end if
           end if

           call CGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
        end if ! fact true

     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        call CGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(21) !!solve the linear system (ZeB-A)^H x=workc(1:N,1:fpm(23)) result in to workc

        call CGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call CGEMM('N','N',N,fpm(25),N,ONEC,A,LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        else
           call CHEMM ('L', UPLO, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call CLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           if ((UPLO=='F').or.(UPLO=='f')) then
              call CGEMM('N','N',N,fpm(25),N,ONEC,B,LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
           else
              call CHEMM ('L', UPLO, N, fpm(25), ONEC, B, LDB, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
           endif
        end if


     end select
  end do

  call wdeallocate_2c(zAq)
  call wdeallocate_2c(zSq)
  call wdeallocate_2c(work)
  call wdeallocate_2c(workc)
  call wdeallocate_3c(Az)
  call wdeallocate_2i(ipivloc)
  call wdeallocate_1c(ztmp)


end subroutine cfeast_hegvx







subroutine sfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  real,dimension(LDA,*):: A
  real,dimension(LDB,*):: B
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
  complex, dimension(:,:),pointer ::work,workc,zAq,zSq
  complex, dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
  logical :: fact
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SFEAST_GEGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

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
!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3c(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
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
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)*ONEC
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else     
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ONEC
           end if


           call CGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
        end if ! fact true


     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
        call CGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        call CGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call CGEMM('N','N',N,fpm(25),N,ONEC,A(1:N,1:N)*ONEC,LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        ! With MKL-blas- one should use 
        !call SCGEMM('N','N',N,fpm(25),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)


     case(31) !! perform multiplication A^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call CGEMM('C','N',N,fpm(35),N,ONEC,A(1:N,1:N)*ONEC,LDA,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
        ! With MKL-blas- one should use
        ! call SCGEMM('T','N',N,fpm(35),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)


     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call CLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
           call CGEMM('N','N',N,fpm(25),N,ONEC,B(1:N,1:N)*ONEC,LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
           ! With MKL-blas- one should use
           ! call SCGEMM('N','N',N,fpm(25),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)


        end if

     case(41)!! perform multiplication B^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        if (LDB==-1) then ! standard
           call CLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
        else
           call CGEMM('C','N',N,fpm(35),N,ONEC,B(1:N,1:N)*ONEC,LDB,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
           ! With MKL-blas- one should use
           !call SCGEMM('T','N',N,fpm(35),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
        end if

     end select
  end do

  call wdeallocate_2c(zAq)
  call wdeallocate_2c(zSq)
  call wdeallocate_2c(work)
  call wdeallocate_2c(workc)
  call wdeallocate_3c(Az)
  call wdeallocate_2i(ipivloc)


end subroutine sfeast_gegvx




subroutine cfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex,dimension(LDA,*):: A
  complex,dimension(LDB,*):: B
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
  complex, dimension(:,:),pointer ::work,workc,zAq,zSq
  complex, dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
  logical :: fact
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'CFEAST_GEGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

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
!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3c(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
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
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else     
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)
           end if


           call CGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
        end if ! fact true


     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc
        call CGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        call CGETRS( 'C', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call CGEMM('N','N',N,fpm(25),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)

     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        call CGEMM('C','N',N,fpm(35),N,ONEC,A(1:N,1:N),LDA,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call CLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
           call CGEMM('N','N',N,fpm(25),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        end if

     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

        if (LDB==-1) then ! standard
           call CLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
        else
           call CGEMM('C','N',N,fpm(35),N,ONEC,B(1:N,1:N),LDB,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
        end if

     end select
  end do

  call wdeallocate_2c(zAq)
  call wdeallocate_2c(zSq)
  call wdeallocate_2c(work)
  call wdeallocate_2c(workc)
  call wdeallocate_3c(Az)
  call wdeallocate_2i(ipivloc)


end subroutine cfeast_gegvx







subroutine cfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B COMPLEX SYMMETRIC :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex,dimension(LDA,*):: A
  complex,dimension(LDB,*):: B
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
  complex, dimension(:,:),pointer ::workc,work,zAq,zSq
  complex, dimension(:,:,:),pointer :: Az
  integer, dimension(:,:),pointer ::ipivloc
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO), ZEROC=(SZERO,SZERO)
  logical :: fact
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

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
     CALL XERBLA( 'CFEAST_SYGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0

  call wallocate_2c(zAq,M0,M0,infoloc)
  call wallocate_2c(zSq,M0,M0,infoloc)
  call wallocate_2c(work,N,M0,infoloc)
  call wallocate_2c(workc,N,M0,infoloc)

  if (infoloc/=0) then
     info=-1
     return
  end if

!!! Format conversion
  UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='L'
  IF ((UPLO=='L').or.(UPLO=='l')) UPLO2='L'
  IF ((UPLO=='U').or.(UPLO=='u')) UPLO2='U'



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
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_3c(Az,N,N,nfact,infoloc)
  call wallocate_2i(ipivloc,N,nfact,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
           if (LDB==-1) then !! standard
              Az(1:N,1:N,id)=-A(1:N,1:N)
              do i=1,N
                 Az(i,i,id)=Az(i,i,id)+Ze
              enddo
           else        
              Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N) 
           end if

           call CSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),workc,N*M0,INFOloc)     
           if (infoloc/=0) then
              info=-2
              return
           end if
        endif ! fact true

     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

        call CSYTRS( UPLO2, N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call CSYMM ('L', UPLO2, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if (LDB==-1) then ! standard
           call CLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call CSYMM ('L', UPLO2, N, fpm(25), ONEC, B, LDB, X(1,fpm(24)), N, ZEROC,work(1,fpm(24)), N)
        end if

     end select
  end do



  call wdeallocate_2c(zAq)
  call wdeallocate_2c(zSq)
  call wdeallocate_2c(work)
  call wdeallocate_2c(workc)
  call wdeallocate_3c(Az)
  call wdeallocate_2i(ipivloc)

end subroutine cfeast_sygvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine sfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  real,dimension(LDA,*):: A
  real,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: sfeast_sygvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex,dimension(fpm(2)) :: Zne,Wne 

  call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call sfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine sfeast_sygv




subroutine cfeast_hegv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex,dimension(LDA,*):: A
  complex,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: cfeast_hegvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex,dimension(fpm(2)) :: Zne,Wne 

  call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call cfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)


end subroutine cfeast_hegv



subroutine sfeast_gegv(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL REAL:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  real,dimension(LDA,*):: A
  real,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: sfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex, dimension(1:fpm(8)) :: Zne,Wne 

  call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call sfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine sfeast_gegv




subroutine cfeast_gegv(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex,dimension(LDA,*):: A
  complex,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: cfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex, dimension(1:fpm(8)) :: Zne,Wne 

  call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call cfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine cfeast_gegv




subroutine cfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A,B GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
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
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
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
  integer :: N,LDA,LDB
  complex,dimension(LDA,*):: A
  complex,dimension(LDB,*):: B
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
!!!    Wrapper Routine to expert routine: cfeast_sygvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex, dimension(1:fpm(8)) :: Zne,Wne 

  call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call cfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine cfeast_sygv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine sfeast_syevx(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  real,dimension(LDA,*):: A
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
  real,dimension(1) :: B ! dummy
  integer :: LDB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_sygvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call sfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine sfeast_syevx





subroutine cfeast_heevx(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: cfeast_hegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex,dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call cfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine cfeast_heevx



subroutine sfeast_geevx(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  real,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: sfeast_gegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real,dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call sfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine sfeast_geevx




subroutine cfeast_geevx(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: cfeast_gegvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex,dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call cfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine cfeast_geevx



subroutine cfeast_syevx(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
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
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: cfeast_sygvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  complex,dimension(1) :: B ! dummy
  integer :: LDB
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call cfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine cfeast_syevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine sfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  real,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: sfeast_sygvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex,dimension(fpm(2)) :: Zne,Wne 
  real,dimension(1) :: B ! dummy
  integer :: LDB

  call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call sfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine sfeast_syev








subroutine cfeast_heev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: cfeast_hegvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex,dimension(fpm(2)) :: Zne,Wne 
  complex,dimension(1) :: B ! dummy
  integer :: LDB

  call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call cfeast_hegvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

end subroutine cfeast_heev




subroutine sfeast_geev(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL REAL :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  real,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: sfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex, dimension(1:fpm(8)) :: Zne,Wne 
  real,dimension(1) :: B ! dummy
  integer :: LDB

  call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call sfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine sfeast_geev




subroutine cfeast_geev(N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  complex,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: cfeast_gegvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex, dimension(1:fpm(8)) :: Zne,Wne 
  complex,dimension(1) :: B ! dummy
  integer :: LDB

  call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call cfeast_gegvx(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine cfeast_geev




subroutine cfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: DENSE FORMAT 
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
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
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
  integer :: N,LDA
  complex,dimension(LDA,*):: A
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
!!!    Wrapper Routine to expert routine: cfeast_sygvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex, dimension(1:fpm(8)) :: Zne,Wne 
  complex,dimension(1) :: B ! dummy
  integer :: LDB

  call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  call cfeast_sygvx(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

end subroutine cfeast_syev







