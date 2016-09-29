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
  
  include './spike-smp/sclbprim.f90' !! Banded primitives
  !include './spike-smp/spike_smp_utilities.f90' !! spike utilities  -- already included in dzfeast_banded.f90
  include './spike-smp/sspike_smp.f90' !! spike real single precision - not used for feast banded
  include './spike-smp/cspike_smp.f90' !! spike complex single precision


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED BANDED INTERFACES !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! SINGLE PRECISION VERSION

  ! List of routines:
  !-------------------

  !{S,C}FEAST_{SB,HB,GB}{EV,GV}


  !111111111111!!!!!! EXPERT ROUTINE for GENERALIZED PROBLEM (CORE)
  ! sfeast_sbgvx
  ! cfeast_hbgvx
  ! sfeast_gbgvx
  ! cfeast_gbgvx
  ! cfeast_sbgvx

  !222222222222!!!!  DEFAULT ROUTINES for GENERALIZED PROBLEM (Wrappers to expert generalized)
  ! sfeast_sbgv
  ! cfeast_hbgv
  ! sfeast_gbgv
  ! cfeast_gbgv
  ! cfeast_sbgv

  !333333333333!!!! EXPERT ROUTINE for STANDARD PROBLEM (Wrappers to expert generalized)
  ! sfeast_sbevx
  ! cfeast_hbevx
  ! sfeast_gbevx
  ! cfeast_gbevx
  ! cfeast_sbevx

  !44444444444!!!!! DEFAULT ROUTINES for STANDARD PROBLEM (Wrappers to expert generalized)
  ! sfeast_sbev
  ! cfeast_hbev
  ! sfeast_gbev
  ! cfeast_gbev
  ! cfeast_sbev




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO)
    logical :: fact
    integer :: nfact,id
    character(len=1) :: UPLO2
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: mlda,mldb,klz,band
    complex, dimension(:),pointer ::ztmp
    !! spike 
    integer,dimension(64) :: spm
    complex, dimension(:,:),pointer ::workspike


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
       CALL XERBLA( 'SFEAST_SBGV', -INFO+100 )
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
    call wallocate_3c(Az,band,N,nfact,infoloc)
    call wallocate_1c(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,klz)
    call wallocate_2c(workspike,klz*klz*spm(10),nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (LDB==-1) then !! standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Az(klz+1:band,1:N,id)=-A(1:kla+1,1:N)*ONEC
                   do i=1,N-1
                      s=min(klz,N-i) 
                      call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze


                elseif ((UPLO=='U').or.(UPLO=='u')) then
                   Az(1:klz+1,1:N,id)=-A(1:klz+1,1:N)*ONEC
                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
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
                         call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:kla+1,1:N,id)=-A(1:kla+1,1:N)*ONEC
                      Az(kla+1-klb:kla+1,1:N,id)=Az(kla+1-klb:kla+1,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=-A(1:band,i)*ONEC
                         ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                         call CAXPY(2*klb+1,ONEC,ztmp(1),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(klb+1:band,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1:klb+1+kla,1:N,id)=Az(klb+1:klb+1+kla,1:N,id)-ONEC*A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:klb+1,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1-kla:klb+1,1:N,id)=Az(klb+1-kla:klb+1,1:N,id)-ONEC*A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=B(1:band,i)*Ze
                         ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                         call CAXPY(2*kla+1,-ONEC,ztmp(1),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif

                end if

             end if !!!!!!!!!!! standard or generalized




!!!!! Factorize
             call CSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),band,workspike(1,id),infoloc)
             if (infoloc/=0) then
                info=-2
                return
             end if

          endif ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call CSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 

       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call SSBMM(UPLO,n,fpm(25),kla,SONE,A(1,1),LDA,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call SLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
          else
             call SSBMM(UPLO,n,fpm(25),klb,SONE,B(1,1),LDB,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)
          end if

       end select
    end do

    call wdeallocate_2c(workspike)



    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2s(work)
    call wdeallocate_2c(workc)
    call wdeallocate_3c(Az)
    call wdeallocate_1c(ztmp)



  end subroutine sfeast_sbgvx





  subroutine cfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    logical :: fact
    integer :: nfact,id
    character(len=1) :: UPLO2
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: mlda,mldb,klz,band
    complex, dimension(:),pointer ::ztmp
    !! spike 
    integer,dimension(64) :: spm
    complex, dimension(:,:),pointer ::workspike
    real :: norm,nzero


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
       CALL XERBLA( 'CFEAST_HBGV', -INFO+100 )
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
    call wallocate_3c(Az,band,N,nfact,infoloc)
    call wallocate_1c(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,klz)
    call wallocate_2c(workspike,klz*klz*spm(10),nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!  

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
             IF (LDB==-1) then !! standard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,A(2,i),1,Az(klz+2,i,id),1)
                      call CSCAL(s,-ONEC,Az(klz+2,i,id),1)
                      call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      call CLACGV(s, Az(klz,i+1,id), 2*klz )
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze

                elseif ((UPLO=='U').or.(UPLO=='u')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,A(klz,i+1),LDA-1,Az(klz,i+1,id),2*klz)
                      call CSCAL(s,-ONEC,Az(klz,i+1,id),2*klz)
                      call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      call CLACGV(s, Az(klz+2,i,id),1)
                      Az(klz+1,i,id)=-A(klz+1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(klz+1,N)+Ze


                else !UPLO=F

                   do i=1,N
                      call CCOPY(band,A(1,i),1,Az(1,i,id),1)
                      call CSCAL(band,-ONEC,Az(1,i,id),1)
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
                         call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                         call CLACGV(s, Az(klz,i+1,id), 2*klz )
                      enddo
                      Az(kla+1:kla+1+klb,1:N,id)=Az(kla+1:kla+1+klb,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klb,N-i)
                         call CCOPY(s,B(2,i),1,ztmp(1),1)
                         call CLACGV(s,ztmp(1),1)
                         call CAXPY(s,Ze,ztmp(1),1,Az(klz,i+1,id),2*klz)
                      enddo


                   elseif ((UPLO=='U').or.(UPLO=='u')) then


                      Az(1:kla+1,1:N,id)=-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                         call CLACGV(s, Az(klz+2,i,id), 1 )
                      enddo
                      Az(kla+1-klb:kla+1,1:N,id)=Az(kla+1-klb:kla+1,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klb,N-i)
                         call CCOPY(s,B(klb,i+1),LDB-1,ztmp(1),1)
                         call CLACGV(s,ztmp(1),1)
                         call CAXPY(s,Ze,ztmp(1),1,Az(klz+2,i,id),1)
                      enddo



                   else !UPLO=F


                      do i=1,N
                         call CCOPY(band,A(1,i),1,Az(1,i,id),1)
                         call CSCAL(band,-ONEC,Az(1,i,id),1)
                         call CAXPY(2*klb+1,Ze,B(1,i),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                else ! kla<klb!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then


                      Az(klz+1:band,1:N,id)=Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                         call CLACGV(s, Az(klz,i+1,id), 2*klz )
                         call CSCAL(s,Ze/conjg(Ze),Az(klz,i+1,id),2*klz)
                      enddo
                      Az(klz+1:klz+1+kla,1:N,id)=Az(klb+1:klb+1+kla,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(kla,N-i)
                         call CCOPY(s,A(2,i),1,ztmp(1),1)
                         call CLACGV(s,ztmp(1),1)
                         call CAXPY(s,-ONEC,ztmp(1),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then


                      Az(1:klz+1,1:N,id)=Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                         call CLACGV(s, Az(klz+2,i,id), 1 )
                         call CSCAL(s,Ze/conjg(Ze),Az(klz+2,i,id),1)
                      enddo
                      Az(klz+1-kla:klz+1,1:N,id)=Az(klz+1-kla:klz+1,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(kla,N-i)
                         call CCOPY(s,A(kla,i+1),LDA-1,ztmp(1),1)
                         call CLACGV(s,ztmp(1),1)
                         call CAXPY(s,-ONEC,ztmp(1),1,Az(klz+2,i,id),1)
                      enddo



                   else !UPLO=F


                      do i=1,N
                         call CCOPY(band,B(1,i),1,Az(1,i,id),1)
                         call CSCAL(band,Ze,Az(1,i,id),1)
                         call CAXPY(2*kla+1,-ONEC,A(1,i),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif

                end if

             end IF ! standard or generalized


             call CSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),band,workspike(1,id),infoloc)   
             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call CSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)

       case(21) !!solve the linear system (ZeB-A)^H x=workc(1:N,1:fpm(23)) result in to workc

          call CSPIKE_GBTRS(spm,'C',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)  


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CHBMM(UPLO2,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call CLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
          else
             call CHBMM(UPLO2,n,fpm(25),klb,ONEC,B(1,1),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       end select
    end do


    call wdeallocate_2c(workspike)



    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)
    call wdeallocate_3c(Az)
    call wdeallocate_1c(ztmp)


  end subroutine cfeast_hbgvx







  subroutine sfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL REAL :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: klz,kuz,band
    !! spike 
    integer,dimension(64) :: spm
    complex, dimension(:,:),pointer ::workspike


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
       CALL XERBLA( 'SFEAST_GBGV', -INFO+100 )
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

    klz=max(kla,klb)
    kuz=max(kua,kub)
    band=klz+kuz+1

    call wallocate_3c(Az,band,N,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,max(klz,kuz))
    call wallocate_2c(workspike,spm(10)*max(klz,kuz)**2,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!! 


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
             call CSPIKE_GBTRF(spm,N,klz,kuz,Az(1,1,id),band,workspike(1,id),infoloc)    
             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call CSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 


       case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call CSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)  


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)


          call CGBMM('N','N',n,fpm(25),kla,kua,ONEC,A(1:kla+kua+1,1:N)*ONEC,kla+kua+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          ! With MKL-blas- one should uncoment the following lbprim routine 
          ! call SCGBMM('N','N',n,fpm(25),kla,kua,ONEC,A(1:kla+kua+1,1:N),kla+kua+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)


       case(31) !! perform multiplication A^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)


          call CGBMM('C','N',N,fpm(35),kla,kua,ONEC,A(1:kla+kua+1,1:N)*ONEC,kla+kua+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
          ! With MKL-blas- one should uncoment the following lbprim routine
          !call SCGBMM('T','N',N,fpm(35),kla,kua,ONEC,A(1:kla+kua+1,1:N),kla+kua+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)


       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call CLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
          else
             call CGBMM('N','N',N,fpm(25),klb,kub,ONEC,B(1:klb+kub+1,1:N)*ONEC,klb+kub+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
             ! With MKL-blas- one should uncoment the following lbprim routine
             !call SCGBMM('N','N',N,fpm(25),klb,kub,ONEC,B(1:klb+kub+1,1:N),klb+kub+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       case(41)!! perform multiplication B^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          if (LDB==-1) then ! standard
             call CLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
          else
             call CGBMM('C','N',N,fpm(35),klb,kub,ONEC,B(1:klb+kub+1,1:N)*ONEC,klb+kub+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
             ! With MKL-blas- one should uncoment the following lbprim routine
             !call SCGBMM('T','N',N,fpm(35),klb,kub,ONEC,B(1+klb+kub+1,1:N),klb+kub+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
          end if

       end select
    end do


    call wdeallocate_2c(workspike)



    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)
    call wdeallocate_3c(Az)



  end subroutine sfeast_gbgvx




  subroutine cfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
    real,parameter :: SONE=1.0e0,SZERO=0.0e0
    complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    logical :: fact
    integer :: nfact,id
    integer :: rank,code,nb_procs,NEW_COMM_WORLD
    !! banded conversion
    integer :: klz,kuz,band
    !! spike 
    integer,dimension(64) :: spm
    complex, dimension(:,:),pointer ::workspike


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
       CALL XERBLA( 'CFEAST_GBGV', -INFO+100 )
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


    klz=max(kla,klb)
    kuz=max(kua,kub)
    band=klz+kuz+1


    call wallocate_3c(Az,band,N,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! spike initilization
    call spikeinit(spm,N,max(klz,kuz))
    call wallocate_2c(workspike,spm(10)*max(klz,kuz)**2,nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!! 



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
             call CSPIKE_GBTRF(spm,N,klz,kuz,Az(1,1,id),band,workspike(1,id),infoloc)        
             if (infoloc/=0) then
                info=-2
                return
             end if
          end if ! fact true


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call CSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 


       case(21) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc

          call CSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N)  


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CGBMM('N','N',n,fpm(25),kla,kua,ONEC,A(1,1),kla+kua+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
        
       case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          call CGBMM('C','N',N,fpm(35),kla,kua,ONEC,A(1,1),kla+kua+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
         
       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call CLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
          else
             call CGBMM('N','N',N,fpm(25),klb,kub,ONEC,B(1,1),klb+kub+1,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)

          if (LDB==-1) then ! standard
             call CLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
          else
             call CGBMM('C','N',N,fpm(35),klb,kub,ONEC,B(1,1),klb+kub+1,X(1,fpm(34)),N,ZEROC,work(1,fpm(34)),N)
          end if

       end select
    end do


    call wdeallocate_2c(workspike)


    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)
    call wdeallocate_3c(Az)


  end subroutine cfeast_gbgvx







  subroutine cfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B COMPLEX SYMMETRIC :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,klb
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
    !! banded conversion
    integer :: mlda,mldb,klz,band
    complex, dimension(:),pointer ::ztmp
    !! spike 
    integer,dimension(64) :: spm
    complex, dimension(:,:),pointer ::workspike


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
       CALL XERBLA( 'CFEAST_SBGV', -INFO+100 )
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
    call wallocate_3c(Az,band,N,nfact,infoloc)
    call wallocate_1c(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! spike initilization
    call spikeinit(spm,N,klz)
    call wallocate_2c(workspike,klz*klz*spm(10),nfact,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if
!!!!!

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

!!!!!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Az(klz+1:band,1:N,id)=-A(1:kla+1,1:N)
                   do i=1,N-1
                      s=min(klz,N-i) 
                      call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze


                elseif ((UPLO=='U').or.(UPLO=='u')) then
                   Az(1:klz+1,1:N,id)=-A(1:klz+1,1:N)
                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
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
                         call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:kla+1,1:N,id)=-A(1:kla+1,1:N)
                      Az(kla+1-klb:kla+1,1:N,id)=Az(kla+1-klb:kla+1,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=-A(1:band,i)
                         ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                         call CAXPY(2*klb+1,ONEC,ztmp(1),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(klb+1:band,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1:klb+1+kla,1:N,id)=Az(klb+1:klb+1+kla,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:klb+1,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1-kla:klb+1,1:N,id)=Az(klb+1-kla:klb+1,1:N,id)-A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:band,i,id)=B(1:band,i)*Ze
                         ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                         call CAXPY(2*kla+1,-ONEC,ztmp(1),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif

                end if


             end if  !! standard or generalized


!!!!! Factorize
             call CSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),band,workspike(1,id),infoloc)   
             if (infoloc/=0) then
                info=-2
                return
             end if
          endif ! fact true

       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:fpm(23)) result in to workc


          call CSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),band,workspike(1,id),workc,N) 


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CSBMM(UPLO2,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)


       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          if (LDB==-1) then ! standard
             call CLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
          else
             call CSBMM(UPLO2,n,fpm(25),klb,ONEC,B(1,1),LDB,X(1,fpm(24)),N,ZEROC,work(1,fpm(24)),N)
          end if

       end select
    end do

    call wdeallocate_2c(workspike)


    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work)
    call wdeallocate_2c(workc)
    call wdeallocate_3c(Az)
    call wdeallocate_1c(ztmp)

  end subroutine cfeast_sbgvx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine sfeast_sbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
!!!    Wrapper Routine to expert routine: sfeast_sbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call sfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_sbgv




  subroutine cfeast_hbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B
    !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
!!!    Wrapper Routine to expert routine: cfeast_hbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call cfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine cfeast_hbgv



  subroutine sfeast_gbgv(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL REAL:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
!!!    Wrapper Routine to expert routine: sfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call sfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_gbgv




  subroutine cfeast_gbgv(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
    !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,LDB,kla,kua,klb,kub
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
!!!    Wrapper Routine to expert routine: cfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call cfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_gbgv




  subroutine cfeast_sbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A,B GENERAL COMPLEX SYMMETRIC:: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
    !  LDB        (input)        INTEGER: Leading dimension of matrix B  
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
    integer :: N,LDA,LDB,kla,klb
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
!!!    Wrapper Routine to expert routine: cfeast_sbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    call cfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_sbgv





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!33333333333333333333333333333333333333333333333333333333333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine sfeast_sbevx(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
    integer :: LDB,klb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: sfeast_sbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call sfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_sbevx





  subroutine cfeast_hbevx(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: cfeast_hbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(1) :: B ! dummy
    integer :: LDB,klb

    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call cfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_hbevx



  subroutine sfeast_gbevx(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL REAL :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: sfeast_gbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    real,dimension(1) :: B ! dummy
    integer :: LDB,klb,kub
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call sfeast_gbgvx(N,kla,kub,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_gbevx




  subroutine cfeast_gbevx(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: cfeast_gbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    complex,dimension(1) :: B ! dummy
    integer :: LDB,klb,kub
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call cfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine cfeast_gbevx



  subroutine cfeast_sbevx(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX SYMMETRIC:: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: cfeast_sbgvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    complex,dimension(1) :: B ! dummy
    integer :: LDB,klb
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call cfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


  end subroutine cfeast_sbevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!44444444444444444444444444444444444444444444444444444444444
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine sfeast_sbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: sfeast_sbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 
    real,dimension(1) :: B ! dummy
    integer :: LDB,klb

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call sfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_sbev






  subroutine cfeast_hbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
    !                           triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: cfeast_hbgvx
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex,dimension(fpm(2)) :: Zne,Wne 
    complex,dimension(1) :: B ! dummy
    integer :: LDB,klb

    call cfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call cfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_hbev




  subroutine sfeast_gbev(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL REAL :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: sfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 
    real,dimension(1) :: B ! dummy
    integer :: LDB,klb,kub

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call sfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine sfeast_gbev




  subroutine cfeast_gbev(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX :: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla,kua
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
!!!    Wrapper Routine to expert routine: cfeast_gbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 
    complex,dimension(1) :: B ! dummy
    integer :: LDB,klb,kub

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    kub=0
    call cfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_gbev




  subroutine cfeast_sbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A GENERAL COMPLEX SYMMETRIC:: BANDED FORMAT 
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
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
    !  LDA        (input)        INTEGER: Leading dimension of matrix A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  fpm        (input/output) INTEGER(*) : FEAST parameters
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
    integer :: N,LDA,kla
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
!!!    Wrapper Routine to expert routine: cfeast_sbgvx
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex, dimension(1:fpm(8)) :: Zne,Wne 
    complex,dimension(1) :: B ! dummy
    integer :: LDB,klb

    call cfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    LDB=-1 ! B is a dummy- option for standard eigenvalue problem
    klb=0
    call cfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  end subroutine cfeast_sbev







