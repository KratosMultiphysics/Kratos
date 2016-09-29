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
!!!!!!!!! Include SPIKE-SMP and Banded Prmitives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  include './spike-smp/dzlbprim.f90' !! Banded primitives
  include './spike-smp/spike_smp_utilities.f90' !! spike utilities
  include './spike-smp/dspike_smp.f90' !! spike real double precision - not used for feast banded
  include './spike-smp/zspike_smp.f90' !! spike complex double precision


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED BANDED INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

  ! List of routines:
  !-------------------

  !{D,Z}FEAST_{SB,HB,GB}{EV,GV}


  !111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE)
  ! dfeast_sbgvx
  ! zfeast_hbgvx
  ! dfeast_gbgvx
  ! zfeast_gbgvx
  ! zfeast_sbgvx

  !222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
  ! dfeast_sbgv
  ! zfeast_hbgv
  ! dfeast_gbgv
  ! zfeast_gbgv
  ! zfeast_sbgv

  !333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
  ! dfeast_sbevx
  ! zfeast_hbevx
  ! dfeast_gbevx
  ! zfeast_gbevx
  ! zfeast_sbevx

  !44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
  ! dfeast_sbev
  ! zfeast_hbev
  ! dfeast_gbev
  ! zfeast_gbev
  ! zfeast_sbev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
    logical :: fact
    integer :: nfact,id
    character(len=1) :: UPLO2
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: mlda,mldb,klz,band
    complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
    !! spike 
    integer,dimension(64) :: spm
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workspike


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------


    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF



    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda ) THEN
       INFO = -105
    ELSE IF (klb>=N) THEN
       INFO=-106
    ELSE IF((LDB<mldb ).and.(LDB/=-1)) THEN
       INFO = -108
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SBGV', -INFO+100 )
       RETURN
    END IF



    infoloc=0
    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!! Format conversion
    UPLO2=UPLO
    IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='U'
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

    klz=max(kla,klb)
    band=2*klz+1
    call wallocate_3z(Az,band,N,nfact,infoloc)
    call wallocate_1z(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,klz)
    call wallocate_2z(workspike,klz*klz*spm(10),nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (LDB==-1) then !! standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Az(klz+1:band,1:N,id)=-A(1:kla+1,1:N)*ONEC
                   do i=1,N-1
                      s=min(klz,N-i) 
                      call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze


                elseif ((UPLO=='U').or.(UPLO=='u')) then
                   Az(1:klz+1,1:N,id)=-A(1:klz+1,1:N)*ONEC
                   do i=1,N-1
                      s=min(klz,N-i)
                      call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      Az(klz+1,i,id)=-A(klz+1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(klz+1,N)+Ze


                else !UPLO=F

                   do i=1,N
                      Az(1:band,i,id)=-A(1:band,i)*ONEC
                      Az(klz+1,i,id)=Az(klz+1,i,id)+Ze
                   end do

                end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             else !! generalized   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(kla+1:band,1:N,id)=-A(1:kla+1,1:N)*ONEC
                      Az(kla+1:kla+1+klb,1:N,id)=Az(kla+1:kla+1+klb,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:kla+1,1:N,id)=-A(1:kla+1,1:N)*ONEC
                      Az(kla+1-klb:kla+1,1:N,id)=Az(kla+1-klb:kla+1,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=-A(1:band,i)*ONEC
                         ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                         call ZAXPY(2*klb+1,ONEC,ztmp(1),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(klb+1:band,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1:klb+1+kla,1:N,id)=Az(klb+1:klb+1+kla,1:N,id)-ONEC*A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:klb+1,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1-kla:klb+1,1:N,id)=Az(klb+1-kla:klb+1,1:N,id)-ONEC*A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=B(1:band,i)*Ze
                         ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                         call ZAXPY(2*kla+1,-ONEC,ztmp(1),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif

                end if

             end if !!!!!!!!!!! standard or generalized




!!!!! Factorize
             call ZSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),band,workspike(1,id),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

          endif ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call ZSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call DSBMM(UPLO,n,fpm(25),kla,DONE,A(1,1),LDA,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call DLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
          else
             call DSBMM(UPLO,n,fpm(25),klb,DONE,B(1,1),LDB,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)
          end if

       end select
    end do

    call wdeallocate_2z(workspike)



    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2d(work)
    call wdeallocate_2z(workc)
    call wdeallocate_3z(Az)
    call wdeallocate_1z(ztmp)



  end subroutine dfeast_sbgvx





  subroutine zfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq
    complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    character(len=1) :: UPLO2
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: mlda,mldb,klz,band
    complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
    !! spike 
    integer,dimension(64) :: spm
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workspike
    double precision :: norm,nzero


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------

    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF




    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda ) THEN
       INFO = -105
    ELSE IF (klb>=N) THEN
       INFO=-106 
    ELSE IF((LDB<mldb ).and.(LDB/=-1)) THEN
       INFO = -108
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_HBGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


!!! Format conversion
    UPLO2=UPLO
    IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='U'
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


    klz=max(kla,klb)
    band=2*klz+1
    call wallocate_3z(Az,band,N,nfact,infoloc)
    call wallocate_1z(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,klz)
    call wallocate_2z(workspike,klz*klz*spm(10),nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!  

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
             IF (LDB==-1) then !! standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call ZCOPY(s,A(2,i),1,Az(klz+2,i,id),1)
                      call ZSCAL(s,-ONEC,Az(klz+2,i,id),1)
                      call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      call ZLACGV(s, Az(klz,i+1,id), 2*klz )
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze

                elseif ((UPLO=='U').or.(UPLO=='u')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call ZCOPY(s,A(klz,i+1),LDA-1,Az(klz,i+1,id),2*klz)
                      call ZSCAL(s,-ONEC,Az(klz,i+1,id),2*klz)
                      call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      call ZLACGV(s, Az(klz+2,i,id),1)
                      Az(klz+1,i,id)=-A(klz+1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(klz+1,N)+Ze


                else !UPLO=F

                   do i=1,N
                      call ZCOPY(band,A(1,i),1,Az(1,i,id),1)
                      call ZSCAL(band,-ONEC,Az(1,i,id),1)
                      Az(klz+1,i,id)=Az(klz+1,i,id)+Ze
                   end do


                end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             else !! generalized   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(kla+1:band,1:N,id)=-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                         call ZLACGV(s, Az(klz,i+1,id), 2*klz )
                      enddo
                      Az(kla+1:kla+1+klb,1:N,id)=Az(kla+1:kla+1+klb,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klb,N-i)
                         call ZCOPY(s,B(2,i),1,ztmp(1),1)
                         call ZLACGV(s,ztmp(1),1)
                         call ZAXPY(s,Ze,ztmp(1),1,Az(klz,i+1,id),2*klz)
                      enddo


                   elseif ((UPLO=='U').or.(UPLO=='u')) then


                      Az(1:kla+1,1:N,id)=-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                         call ZLACGV(s, Az(klz+2,i,id), 1 )
                      enddo
                      Az(kla+1-klb:kla+1,1:N,id)=Az(kla+1-klb:kla+1,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klb,N-i)
                         call ZCOPY(s,B(klb,i+1),LDB-1,ztmp(1),1)
                         call ZLACGV(s,ztmp(1),1)
                         call ZAXPY(s,Ze,ztmp(1),1,Az(klz+2,i,id),1)
                      enddo



                   else !UPLO=F


                      do i=1,N
                         call ZCOPY(band,A(1,i),1,Az(1,i,id),1)
                         call ZSCAL(band,-ONEC,Az(1,i,id),1)
                         call ZAXPY(2*klb+1,Ze,B(1,i),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                else ! kla<klb!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then


                      Az(klz+1:band,1:N,id)=Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                         call ZLACGV(s, Az(klz,i+1,id), 2*klz )
                         call ZSCAL(s,Ze/conjg(Ze),Az(klz,i+1,id),2*klz)
                      enddo
                      Az(klz+1:klz+1+kla,1:N,id)=Az(klb+1:klb+1+kla,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(kla,N-i)
                         call ZCOPY(s,A(2,i),1,ztmp(1),1)
                         call ZLACGV(s,ztmp(1),1)
                         call ZAXPY(s,-ONEC,ztmp(1),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then


                      Az(1:klz+1,1:N,id)=Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                         call ZLACGV(s, Az(klz+2,i,id), 1 )
                         call ZSCAL(s,Ze/conjg(Ze),Az(klz+2,i,id),1)
                      enddo
                      Az(klz+1-kla:klz+1,1:N,id)=Az(klz+1-kla:klz+1,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(kla,N-i)
                         call ZCOPY(s,A(kla,i+1),LDA-1,ztmp(1),1)
                         call ZLACGV(s,ztmp(1),1)
                         call ZAXPY(s,-ONEC,ztmp(1),1,Az(klz+2,i,id),1)
                      enddo



                   else !UPLO=F


                      do i=1,N
                         call ZCOPY(band,B(1,i),1,Az(1,i,id),1)
                         call ZSCAL(band,Ze,Az(1,i,id),1)
                         call ZAXPY(2*kla+1,-ONEC,A(1,i),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif

                end if

             end IF ! standard or generalized


             call ZSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),band,workspike(1,id),infoloc)   
             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call ZSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)

       case(21) !!solve the linear system (ZeB-A)^H x=workc(1:N,1:fpm(23)) result in to workc

          call ZSPIKE_GBTRS(spm,'C',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)  


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZHBMM(UPLO2,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
          else
             call ZHBMM(UPLO2,n,fpm(25),klb,ONEC,B(1,1),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       end select
    end do


    call wdeallocate_2z(workspike)



    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_3z(Az)
    call wdeallocate_1z(ztmp)


  end subroutine zfeast_hbgvx







  subroutine dfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL REAL :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq
    complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: klz,kuz,band
    !! spike 
    integer,dimension(64) :: spm
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workspike


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
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF (kua>=N) THEN
       INFO=-104
    ELSE IF(LDA<kla+kua+1 ) THEN
       INFO = -106
    ELSE IF (klb>=N) THEN
       INFO=-107
    ELSE IF (kub>=N) THEN
       INFO=-108
    ELSE IF((LDB<klb+kub+1).and.(LDB/=-1)) THEN
       INFO = -110
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_GBGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!

    klz=max(kla,klb)
    kuz=max(kua,kub)
    band=klz+kuz+1

    call wallocate_3z(Az,band,N,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,max(klz,kuz))
    call wallocate_2z(workspike,spm(10)*max(klz,kuz)**2,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!! 


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
             if (LDB==-1) then !! standard
                Az(1:band,N,id)=ZEROC
                do i=1,N
                   Az(1:band,i,id) = -A(1:band,i)*ONEC
                   Az(kua+1,i,id) = Az(kua+1,i,id) + Ze
                end do

             else  
                Az(1:band,N,id)=ZEROC
                Az(kuz-kua+1:kuz+kla+1,1:N,id) = -ONEC*A(1:kla+kua+1,1:N)
                Az(kuz-kub+1:kuz+klb+1,1:N,id) = Az(kuz-kub+1:kuz+klb+1,1:N,id) + Ze*B(1:kub+klb+1,1:N)

             end if

!!!!! Factorize
             call ZSPIKE_GBTRF(spm,N,klz,kuz,Az(1,1,id),band,workspike(1,id),infoloc)    
             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call ZSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 


       case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call ZSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)  


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZGBMM('N','N',n,fpm(25),kla,kua,ONEC,A(1:kla+kua+1,1:N)*ONEC,kla+kua+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          ! With MKL-blas- one should uncoment the following lbprim routine 
!          call DZGBMM('N','N',n,fpm(25),kla,kua,ONEC,A(1:kla+kua+1,1:N),kla+kua+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)


       case(31) !! perform multiplication A^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call ZGBMM('C','N',N,fpm(35),kla,kua,ONEC,A(1:kla+kua+1,1:N)*ONEC,kla+kua+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
          ! With MKL-blas- one should uncoment the following lbprim routine
!          call DZGBMM('T','N',N,fpm(35),kla,kua,ONEC,A(1:kla+kua+1,1:N),kla+kua+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)


       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
          else
             call ZGBMM('N','N',N,fpm(25),klb,kub,ONEC,B(1:klb+kub+1,1:N)*ONEC,klb+kub+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
             ! With MKL-blas- one should uncoment the following lbprim routine
!             call DZGBMM('N','N',N,fpm(25),klb,kub,ONEC,B(1:klb+kub+1,1:N),klb+kub+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       case(41)!! perform multiplication B^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          if (LDB==-1) then ! standard
             call ZLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
          else
             call ZGBMM('C','N',N,fpm(35),klb,kub,ONEC,B(1:klb+kub+1,1:N)*ONEC,klb+kub+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
             ! With MKL-blas- one should uncoment the following lbprim routine
 !            call DZGBMM('T','N',N,fpm(35),klb,kub,ONEC,B(1+klb+kub+1,1:N),klb+kub+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
          end if
       end select
    end do



    call wdeallocate_2z(workspike)


    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_3z(Az)



  end subroutine dfeast_gbgvx




  subroutine zfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work,workc,zAq,zSq
    complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: klz,kuz,band
    !! spike 
    integer,dimension(64) :: spm
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workspike


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
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF (kua>=N) THEN
       INFO=-104
    ELSE IF(LDA<kla+kua+1 ) THEN
       INFO = -106
    ELSE IF (klb>=N) THEN
       INFO=-107
    ELSE IF (kub>=N) THEN
       INFO=-108
    ELSE IF((LDB<klb+kub+1).and.(LDB/=-1)) THEN
       INFO = -110
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_GBGV', -INFO+100 )
       RETURN
    END IF


    infoloc=0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!


    klz=max(kla,klb)
    kuz=max(kua,kub)
    band=klz+kuz+1


    call wallocate_3z(Az,band,N,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,max(klz,kuz))
    call wallocate_2z(workspike,spm(10)*max(klz,kuz)**2,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!! 



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
             if (LDB==-1) then !! standard

                Az(1:band,N,id)=ZEROC
                do i=1,N
                   Az(1:band,i,id) = -A(1:band,i)
                   Az(kua+1,i,id) = Az(kua+1,i,id) + Ze
                end do

             else     

                Az(1:band,N,id)=ZEROC
                Az(kuz-kua+1:kuz+kla+1,1:N,id) = -A(1:kla+kua+1,1:N)
                Az(kuz-kub+1:kuz+klb+1,1:N,id) = Az(kuz-kub+1:kuz+klb+1,1:N,id) + Ze*B(1:kub+klb+1,1:N)

             end if

!!!!! Factorize
             call ZSPIKE_GBTRF(spm,N,klz,kuz,Az(1,1,id),band,workspike(1,id),infoloc)        
             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call ZSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 


       case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call ZSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)  


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZGBMM('N','N',n,fpm(25),kla,kua,ONEC,A(1,1),kla+kua+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
         

       case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call ZGBMM('C','N',N,fpm(35),kla,kua,ONEC,A(1,1),kla+kua+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
          else  
             call ZGBMM('N','N',N,fpm(25),klb,kub,ONEC,B(1,1),klb+kub+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          if (LDB==-1) then ! standard
             call ZLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
          else
             call ZGBMM('C','N',N,fpm(35),klb,kub,ONEC,B(1,1),klb+kub+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)            
          end if

       end select
    end do


    call wdeallocate_2z(workspike)


    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_3z(Az)


  end subroutine zfeast_gbgvx







  subroutine zfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B COMPLEX SYMMETRIC :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,klb
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
    integer :: ijob,infoloc,i,s
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,work,zAq,zSq
    complex(kind=(kind(1.0d0))), dimension(:,:,:),pointer :: Az
    integer, dimension(:,:),pointer ::ipivloc
    double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO), ZEROC=(DZERO,DZERO)
    logical :: fact
    integer :: nfact,id
    character(len=1) :: UPLO2
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: mlda,mldb,klz,band
    complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
    !! spike 
    integer,dimension(64) :: spm
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workspike


    rank=0
    nb_procs=1
    !----------------------------------------------
#ifdef MPI
    NEW_COMM_WORLD=fpm(9)
    call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
    call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
#endif
    !----------------------


    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda ) THEN
       INFO = -105
    ELSE IF (klb>=N) THEN
       INFO=-106
    ELSE IF((LDB<mldb ).and.(LDB/=-1)) THEN
       INFO = -108
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_SBGV', -INFO+100 )
       RETURN
    END IF



    infoloc=0
    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if

!!! Format conversion
    UPLO2=UPLO
    IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='U'
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


    klz=max(kla,klb)
    band=2*klz+1
    call wallocate_3z(Az,band,N,nfact,infoloc)
    call wallocate_1z(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! spike initilization
    call spikeinit(spm,N,klz)
    call wallocate_2z(workspike,klz*klz*spm(10),nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!

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

             if (LDB==-1) then !! standard

!!!!!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Az(klz+1:band,1:N,id)=-A(1:kla+1,1:N)
                   do i=1,N-1
                      s=min(klz,N-i) 
                      call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze


                elseif ((UPLO=='U').or.(UPLO=='u')) then
                   Az(1:klz+1,1:N,id)=-A(1:klz+1,1:N)
                   do i=1,N-1
                      s=min(klz,N-i)
                      call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      Az(klz+1,i,id)=-A(klz+1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(klz+1,N)+Ze


                else !UPLO=F

                   do i=1,N
                      Az(1:band,i,id)=-A(1:band,i)
                      Az(klz+1,i,id)=Az(klz+1,i,id)+Ze
                   end do

                end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             else !! generalized   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(kla+1:band,1:N,id)=-A(1:kla+1,1:N)
                      Az(kla+1:kla+1+klb,1:N,id)=Az(kla+1:kla+1+klb,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:kla+1,1:N,id)=-A(1:kla+1,1:N)
                      Az(kla+1-klb:kla+1,1:N,id)=Az(kla+1-klb:kla+1,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=-A(1:band,i)
                         ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                         call ZAXPY(2*klb+1,ONEC,ztmp(1),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(klb+1:band,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1:klb+1+kla,1:N,id)=Az(klb+1:klb+1+kla,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:klb+1,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1-kla:klb+1,1:N,id)=Az(klb+1-kla:klb+1,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=B(1:band,i)*Ze
                         ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                         call ZAXPY(2*kla+1,-ONEC,ztmp(1),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif

                end if


             end if  !! standard or generalized


!!!!! Factorize
             call ZSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),band,workspike(1,id),infoloc)   
             if (infoloc/=0) then
                info=-2
                return
             end if
          endif ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc


          call ZSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZSBMM(UPLO2,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)


       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
          else
             call ZSBMM(UPLO2,n,fpm(25),klb,ONEC,B(1,1),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       end select
    end do

    call wdeallocate_2z(workspike)


    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work)
    call wdeallocate_2z(workc)
    call wdeallocate_3z(Az)
    call wdeallocate_1z(ztmp)

  end subroutine zfeast_sbgvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine dfeast_sbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
!!!    Wrapper Routine to expert routine: dfeast_sbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call dfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_sbgv




  subroutine zfeast_hbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B
    !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
!!!    Wrapper Routine to expert routine: zfeast_hbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call zfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine zfeast_hbgv



  subroutine dfeast_gbgv(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL REAL:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
!!!    Wrapper Routine to expert routine: dfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call dfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_gbgv




  subroutine zfeast_gbgv(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
!!!    Wrapper Routine to expert routine: zfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call zfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine zfeast_gbgv




  subroutine zfeast_sbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX SYMMETRIC:: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
!!!    Wrapper Routine to expert routine: zfeast_sbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call zfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine zfeast_sbgv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine dfeast_sbevx(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
    integer :: LDB,klb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_sbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call dfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_sbevx





  subroutine zfeast_hbevx(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: zfeast_hbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
    integer :: LDB,klb

    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call zfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine zfeast_hbevx



  subroutine dfeast_gbevx(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL REAL :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: dfeast_gbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    double precision,dimension(1) :: B ! dummy
    integer :: LDB,klb,kub
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call dfeast_gbgvx(N,kla,kub,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_gbevx




  subroutine zfeast_gbevx(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: zfeast_gbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
    integer :: LDB,klb,kub
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call zfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine zfeast_gbevx



  subroutine zfeast_sbevx(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX SYMMETRIC:: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: zfeast_sbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
    integer :: LDB,klb
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call zfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine zfeast_sbevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine dfeast_sbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: dfeast_sbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
    double precision,dimension(1) :: B ! dummy
    integer :: LDB,klb

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call dfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_sbev






  subroutine zfeast_hbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: zfeast_hbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 
    complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
    integer :: LDB,klb

    call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call zfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine zfeast_hbev




  subroutine dfeast_gbev(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL REAL :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: dfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
    double precision,dimension(1) :: B ! dummy
    integer :: LDB,klb,kub

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call dfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine dfeast_gbev




  subroutine zfeast_gbev(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX :: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: zfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
    complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
    integer :: LDB,klb,kub

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call zfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine zfeast_gbev




  subroutine zfeast_sbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX SYMMETRIC:: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: zfeast_sbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 
    complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
    integer :: LDB,klb

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call zfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine zfeast_sbev







