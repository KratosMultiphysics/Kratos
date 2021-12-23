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
!!!!!!!!!!!!!!!!! FEAST PREDEFINED BANDED INTERFACES !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

! List of routines:
!-------------------

!{D,Z}FEAST_{SB,HB,GB}{EV,GV}{X}


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


!!! Auxiliary
!mDZGBMM           :: Same than ZGBMM but A matrix is real --> works only with DZGEMM (MKL-BLAS)- "m" in front to avoid duplication with current spike banded primitive routines in spike-v1



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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=2*klb+1)
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork
  double precision, dimension(:,:),allocatable ::work,Aq,Bq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
   complex, dimension(:,:),allocatable ::cwork ! mixed precision
  real,dimension(:,:),allocatable :: sA,sB ! mixed precision
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
  !! banded conversion
  integer :: mlda,mldb,bandz,klz
  complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
  complex, dimension(:),pointer ::ctmp
  !! spike 
  integer,dimension(64) :: spm
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::zworkspike
  complex, dimension(:,:),pointer ::cworkspike


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121135 ! code name
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

    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF

  
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
     CALL XERBLA( 'DFEAST_SBGVX', -INFO+100 )
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
     do i=rank+1,fpm(2),nb_procs
        nfact=nfact+1
     end do
  endif

 klz=max(kla,klb)
 bandz=2*klz+1

!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(sA(LDA,N)) !!<< LDA may be overestimated (can be optimized)
     call DLAG2S(LDA,N,A,LDA,sA,LDA,infoloc1)
     if (LDB/=-1) then
        allocate(sB(LDB,N))  !!<< LDB may be overestimated (can be optimized)
        call DLAG2S(LDB,N,B,LDB,sB,LDB,infoloc1)
     endif
     allocate(Ac(bandz,N,nfact))
     allocate(ctmp(bandz))
  else 
     allocate(Az(bandz,N,nfact))
      allocate(ztmp(bandz))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up SPIKE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call spikeinit(spm,N,klz)
  if (fpm(42)==1) then
     allocate(cworkspike(klz*klz*spm(10),nfact))
      allocate(cwork(N,M0))
     else
        allocate(zworkspike(klz*klz*spm(10),nfact))
     end if


  
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

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (fpm(42)==1) then ! single precision
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           if (LDB==-1) then !! standard 
              !Ac(1:N,1:N,id)=-sA(1:N,1:N)*CONE
              !do i=1,N
              !   Ac(i,i,id)=Ac(i,i,id)+cmplx(Ze)
              !enddo
!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Ac(klz+1:bandz,1:N,id)=-sA(1:kla+1,1:N)*CONE
                   do i=1,N-1
                      s=min(klz,N-i) 
                      call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                      Ac(klz+1,i,id)=-sA(1,i)+cmplx(Ze)
                   enddo
                   Ac(klz+1,N,id)=-sA(1,N)+cmplx(Ze)


                elseif ((UPLO=='U').or.(UPLO=='u')) then
                   Ac(1:klz+1,1:N,id)=-sA(1:klz+1,1:N)*CONE
                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                      Ac(klz+1,i,id)=-sA(klz+1,i)+cmplx(Ze)
                   enddo
                   Ac(klz+1,N,id)=-sA(klz+1,N)+cmplx(Ze)


                else !UPLO=F

                   do i=1,N
                      Ac(1:bandz,i,id)=-sA(1:bandz,i)*CONE
                      Ac(klz+1,i,id)=Ac(klz+1,i,id)+cmplx(Ze)
                   end do

                end if
              
           else  !! generalized       
              !Ac(1:N,1:N,id)=cmplx(Ze)*sB(1:N,1:N)-sA(1:N,1:N)*CONE
!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Ac(kla+1:bandz,1:N,id)=-sA(1:kla+1,1:N)*CONE
                      Ac(kla+1:kla+1+klb,1:N,id)=Ac(kla+1:kla+1+klb,1:N,id)+cmplx(Ze)*sB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                      enddo

                   elseif ((UPLO=='U').or.(UPLO=='u')) then
                      Ac(1:kla+1,1:N,id)=-sA(1:kla+1,1:N)*CONE
                      Ac(kla+1-klb:kla+1,1:N,id)=Ac(kla+1-klb:kla+1,1:N,id)+cmplx(Ze)*sB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Ac(1:bandz,i,id)=-sA(1:bandz,i)*CONE
                         ctmp(1:2*klb+1)=cmplx(Ze)*sB(1:2*klb+1,i)
                         call CAXPY(2*klb+1,CONE,ctmp(1),1,Ac(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Ac(klb+1:bandz,1:N,id)=cmplx(Ze)*sB(1:klb+1,1:N)
                      Ac(klb+1:klb+1+kla,1:N,id)=Ac(klb+1:klb+1+kla,1:N,id)-CONE*sA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then
                      Ac(1:klb+1,1:N,id)=cmplx(Ze)*sB(1:klb+1,1:N)
                      Ac(klb+1-kla:klb+1,1:N,id)=Ac(klb+1-kla:klb+1,1:N,id)-CONE*sA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Ac(1:bandz,i,id)=sB(1:bandz,i)*cmplx(Ze)
                         ctmp(1:2*kla+1)=sA(1:2*kla+1,i)
                         call CAXPY(2*kla+1,-CONE,ctmp(1),1,Ac(klb+1-kla,i,id),1)
                      end do


                   endif

                end if
              
           end if

           

         call CSPIKE_GBTRF(spm,N,klz,klz,Ac(1,1,id),bandz,cworkspike(1,id),infoloc2)    
        

           
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
        else ! double precision
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

           if (LDB==-1) then !! standard
              !Az(1:N,1:N,id)=-A(1:N,1:N)*ZONE
              !do i=1,N
              !   Az(i,i,id)=Az(i,i,id)+Ze
              !enddo
!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Az(klz+1:bandz,1:N,id)=-A(1:kla+1,1:N)*ZONE
                   do i=1,N-1
                      s=min(klz,N-i) 
                      call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze


                elseif ((UPLO=='U').or.(UPLO=='u')) then
                   Az(1:klz+1,1:N,id)=-A(1:klz+1,1:N)*ZONE
                   do i=1,N-1
                      s=min(klz,N-i)
                      call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      Az(klz+1,i,id)=-A(klz+1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(klz+1,N)+Ze


                else !UPLO=F

                   do i=1,N
                      Az(1:bandz,i,id)=-A(1:bandz,i)*ZONE
                      Az(klz+1,i,id)=Az(klz+1,i,id)+Ze
                   end do

                end if

              
           else !! generalized
              !Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ZONE 

!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(kla+1:bandz,1:N,id)=-A(1:kla+1,1:N)*ZONE
                      Az(kla+1:kla+1+klb,1:N,id)=Az(kla+1:kla+1+klb,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:kla+1,1:N,id)=-A(1:kla+1,1:N)*ZONE
                      Az(kla+1-klb:kla+1,1:N,id)=Az(kla+1-klb:kla+1,1:N,id)+Ze*B(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:bandz,i,id)=-A(1:bandz,i)*ZONE
                         ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                         call ZAXPY(2*klb+1,ZONE,ztmp(1),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(klb+1:bandz,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1:klb+1+kla,1:N,id)=Az(klb+1:klb+1+kla,1:N,id)-ZONE*A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then
                      Az(1:klb+1,1:N,id)=Ze*B(1:klb+1,1:N)
                      Az(klb+1-kla:klb+1,1:N,id)=Az(klb+1-kla:klb+1,1:N,id)-ZONE*A(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Az(1:bandz,i,id)=B(1:bandz,i)*Ze
                         ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                         call ZAXPY(2*kla+1,-ZONE,ztmp(1),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif
                end if
           end if


         call ZSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),bandz,zworkspike(1,id),infoloc2)    
        
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
            call CSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
       
 call ZSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N) 

           
        end if
        if (infoloc1/=0) then
           info=-1
           return
        end if
       ! if (infoloc2/=0) then
       !    info=-2
       !    return
       ! end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
 call DSBMM(UPLO,n,fpm(25),kla,DONE,A(1,1),LDA,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call DLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
            call DSBMM(UPLO,n,fpm(25),klb,DONE,B(1,1),LDB,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)
        end if

     end select
  end do



  deallocate(Aq)
  deallocate(Bq)
  deallocate(work)
  deallocate(zwork)
  if (fpm(42)==1) then
     deallocate(cworkspike)
     deallocate(cwork)
     deallocate(ctmp)
     deallocate(sA)
      if (LDB/=-1)  deallocate(sB)
     deallocate(Ac)
  else
      deallocate(zworkspike)
      deallocate(Az)
      deallocate(ztmp)
  end if


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
  !  A          (input/output) COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A
  !                            exit: full format if lower or upper 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
  !  B          (input/output) COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B
  !                            exit: full format if lower or upper        
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=2*klb+1)
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,zAq,zBq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az 
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:),allocatable :: cA,cB ! mixed precision
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
 character(len=1) :: UPLO2
 !! banded conversion
  integer :: mlda,mldb,bandz,klz
  complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
  complex, dimension(:),pointer ::ctmp
  !! spike 
  integer,dimension(64) :: spm
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::zworkspike
  complex, dimension(:,:),pointer ::cworkspike


 
  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141235 ! code name
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
       CALL XERBLA( 'ZFEAST_HBGVX', -INFO+100 )
       RETURN
    END IF

  infoloc1=0
  infoloc2=0
!!!!
  !! for case (30,40)
 UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='U' !L will not be consistant  here with ZHBMM 
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


 klz=max(kla,klb)
 bandz=2*klz+1

  
!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(cA(LDA,N)) !!<< LDA may be overestimated (can be optimized)
     call ZLAG2C(LDA,N,A,LDA,cA,LDA,infoloc1)
     if (LDB/=-1) then
        allocate(cB(LDB,N))!!<< LDB may be overestimated (can be optimized)
        call ZLAG2C(LDB,N,B,LDB,cB,LDB,infoloc1)
     endif
     allocate(Ac(bandz,N,nfact))
      allocate(ctmp(bandz))
  else
     allocate(Az(bandz,N,nfact))
        allocate(ztmp(bandz))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up SPIKE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call spikeinit(spm,N,klz)
  if (fpm(42)==1) then
     allocate(cworkspike(klz*klz*spm(10),nfact))
      allocate(cwork(N,M0))
     else
        allocate(zworkspike(klz*klz*spm(10),nfact))
     end if

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

!!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,cA(2,i),1,Ac(klz+2,i,id),1)
                      call CSCAL(s,-CONE,Ac(klz+2,i,id),1)
                      call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                      call CLACGV(s, Ac(klz,i+1,id), 2*klz )
                      Ac(klz+1,i,id)=-cA(1,i)+cmplx(Ze)
                   enddo
                   Ac(klz+1,N,id)=-cA(1,N)+cmplx(Ze)

                elseif ((UPLO=='U').or.(UPLO=='u')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,cA(klz,i+1),LDA-1,Ac(klz,i+1,id),2*klz)
                      call CSCAL(s,-CONE,Ac(klz,i+1,id),2*klz)
                      call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                      call CLACGV(s, Ac(klz+2,i,id),1)
                      Ac(klz+1,i,id)=-cA(klz+1,i)+cmplx(Ze)
                   enddo
                   Ac(klz+1,N,id)=-cA(klz+1,N)+cmplx(Ze)


                else !UPLO=F

                   do i=1,N
                      call CCOPY(bandz,cA(1,i),1,Ac(1,i,id),1)
                      call CSCAL(bandz,-CONE,Ac(1,i,id),1)
                      Ac(klz+1,i,id)=Ac(klz+1,i,id)+cmplx(Ze)
                   end do


                end if



              
             else !! generalized
                !!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Ac(kla+1:bandz,1:N,id)=-cA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                         call CLACGV(s, Ac(klz,i+1,id), 2*klz )
                      enddo
                      Ac(kla+1:kla+1+klb,1:N,id)=Ac(kla+1:kla+1+klb,1:N,id)+cmplx(Ze)*cB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klb,N-i)
                         call CCOPY(s,cB(2,i),1,ctmp(1),1)
                         call CLACGV(s,ctmp(1),1)
                         call CAXPY(s,cmplx(Ze),ctmp(1),1,Ac(klz,i+1,id),2*klz)
                      enddo


                   elseif ((UPLO=='U').or.(UPLO=='u')) then


                      Ac(1:kla+1,1:N,id)=-cA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                         call CLACGV(s, Ac(klz+2,i,id), 1 )
                      enddo
                      Ac(kla+1-klb:kla+1,1:N,id)=Ac(kla+1-klb:kla+1,1:N,id)+cmplx(Ze)*cB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klb,N-i)
                         call CCOPY(s,cB(klb,i+1),LDB-1,ctmp(1),1)
                         call CLACGV(s,ctmp(1),1)
                         call CAXPY(s,cmplx(Ze),ctmp(1),1,Ac(klz+2,i,id),1)
                      enddo



                   else !UPLO=F


                      do i=1,N
                         call CCOPY(bandz,cA(1,i),1,Ac(1,i,id),1)
                         call CSCAL(bandz,-CONE,Ac(1,i,id),1)
                         call CAXPY(2*klb+1,cmplx(Ze),cB(1,i),1,Ac(kla+1-klb,i,id),1)
                      end do

                   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                else ! kla<klb!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then


                      Ac(klz+1:bandz,1:N,id)=cmplx(Ze)*cB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                         call CLACGV(s, Ac(klz,i+1,id), 2*klz )
                         call CSCAL(s,cmplx(Ze)/conjg(cmplx(Ze)),Ac(klz,i+1,id),2*klz)
                      enddo
                      Ac(klz+1:klz+1+kla,1:N,id)=Ac(klb+1:klb+1+kla,1:N,id)-cA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(kla,N-i)
                         call CCOPY(s,cA(2,i),1,ctmp(1),1)
                         call CLACGV(s,ctmp(1),1)
                         call CAXPY(s,-CONE,ctmp(1),1,Ac(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then


                      Ac(1:klz+1,1:N,id)=cmplx(Ze)*cB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                         call CLACGV(s, Ac(klz+2,i,id), 1 )
                         call CSCAL(s,cmplx(Ze)/conjg(cmplx(Ze)),Ac(klz+2,i,id),1)
                      enddo
                      Ac(klz+1-kla:klz+1,1:N,id)=Ac(klz+1-kla:klz+1,1:N,id)-cA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(kla,N-i)
                         call CCOPY(s,cA(kla,i+1),LDA-1,ctmp(1),1)
                         call CLACGV(s,ctmp(1),1)
                         call CAXPY(s,-CONE,ctmp(1),1,Ac(klz+2,i,id),1)
                      enddo



                   else !UPLO=F


                      do i=1,N
                         call CCOPY(bandz,cB(1,i),1,Ac(1,i,id),1)
                         call CSCAL(bandz,cmplx(Ze),Ac(1,i,id),1)
                         call CAXPY(2*kla+1,-CONE,cA(1,i),1,Ac(klb+1-kla,i,id),1)
                      end do


                   endif

                end if
                



              
           end if

         
              call CSPIKE_GBTRF(spm,N,klz,klz,Ac(1,1,id),bandz,cworkspike(1,id),infoloc2)   


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
        else ! double precision
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

           if (LDB==-1) then !! standard
              !!!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call ZCOPY(s,A(2,i),1,Az(klz+2,i,id),1)
                      call ZSCAL(s,-ZONE,Az(klz+2,i,id),1)
                      call ZCOPY(s,Az(klz+2,i,id),1,Az(klz,i+1,id),2*klz)
                      call ZLACGV(s, Az(klz,i+1,id), 2*klz )
                      Az(klz+1,i,id)=-A(1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(1,N)+Ze

                elseif ((UPLO=='U').or.(UPLO=='u')) then

                   do i=1,N-1
                      s=min(klz,N-i)
                      call ZCOPY(s,A(klz,i+1),LDA-1,Az(klz,i+1,id),2*klz)
                      call ZSCAL(s,-ZONE,Az(klz,i+1,id),2*klz)
                      call ZCOPY(s,Az(klz,i+1,id),2*klz,Az(klz+2,i,id),1)
                      call ZLACGV(s, Az(klz+2,i,id),1)
                      Az(klz+1,i,id)=-A(klz+1,i)+Ze
                   enddo
                   Az(klz+1,N,id)=-A(klz+1,N)+Ze


                else !UPLO=F

                   do i=1,N
                      call ZCOPY(bandz,A(1,i),1,Az(1,i,id),1)
                      call ZSCAL(bandz,-ZONE,Az(1,i,id),1)
                      Az(klz+1,i,id)=Az(klz+1,i,id)+Ze
                   end do


                end if

             
             else !! generalized

!!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(kla+1:bandz,1:N,id)=-A(1:kla+1,1:N)
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
                         call ZCOPY(bandz,A(1,i),1,Az(1,i,id),1)
                         call ZSCAL(bandz,-ZONE,Az(1,i,id),1)
                         call ZAXPY(2*klb+1,Ze,B(1,i),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                else ! kla<klb!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then


                      Az(klz+1:bandz,1:N,id)=Ze*B(1:klb+1,1:N)
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
                         call ZAXPY(s,-ZONE,ztmp(1),1,Az(klz,i+1,id),2*klz)
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
                         call ZAXPY(s,-ZONE,ztmp(1),1,Az(klz+2,i,id),1)
                      enddo



                   else !UPLO=F


                      do i=1,N
                         call ZCOPY(bandz,B(1,i),1,Az(1,i,id),1)
                         call ZSCAL(bandz,Ze,Az(1,i,id),1)
                         call ZAXPY(2*kla+1,-ZONE,A(1,i),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif

                end if
                
                

              
             end if

              call ZSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),bandz,zworkspike(1,id),infoloc2)   


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
            call CSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else

           
call ZSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N)
   
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

       ! if (infoloc2/=0) then
       !    info=-2
       !    return
       ! end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        id=fpm(33)

  if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
            call CSPIKE_GBTRS(spm,'C',N,klz,klz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
    
      call ZSPIKE_GBTRS(spm,'C',N,klz,klz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N)
     
        end if

        if (infoloc1/=0) then
           info=-1
           return
        end if

   !     if (infoloc2/=0) then
   !        info=-2
   !        return
   !     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

       
 call ZHBMM(UPLO2,n,fpm(25),kla,ZONE,A(1,1),LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call ZHBMM(UPLO2,n,fpm(25),klb,ZONE,B(1,1),LDB,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
                 
        end if


     end select
  end do

  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)
  if (fpm(42)==1) then
     deallocate(cworkspike)
     deallocate(cwork)
     deallocate(ctmp)
     deallocate(cA)
      if (LDB/=-1)  deallocate(cB)
     deallocate(Ac)
  else
     deallocate(zworkspike)
     deallocate(Az)
       deallocate(ztmp)
  end if



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
  !
  !  N          (input)        INTEGER: Size system
  !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
  !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=kla+kua+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
  !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=klb+kub+1)
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
  ! Eric Polizzi 2009-2019
  ! ====================================================================
  implicit none
  !  include 'f90_noruntime_interface.fi'
  integer :: N,LDA,LDB,kla,klb,kua,kub
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
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO), CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
 !! banded conversion
  integer :: mlda,mldb,bandz,klz,kuz
  !! spike 
  integer,dimension(64) :: spm
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::zworkspike
  complex, dimension(:,:),pointer ::cworkspike


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121335 ! code name
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
       CALL XERBLA( 'DFEAST_GBGVX', -INFO+100 )
       RETURN
    END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! :-( If DZGEMM (inside DZGBMM)  not available for case 30-40, we need to make a complex copy of the matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
#ifdef MKL
  ! everything fine no copy needed
#else
  ! (extra copy)
  allocate(zA(LDA,N))
  call ZLACP2( 'F', LDA, N,A , LDA, zA, LDA ) 
  if (LDB/=-1) then
     allocate(zB(LDB,N))
     call ZLACP2( 'F', LDB, N,B , LDB, zB, LDB ) 
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


    klz=max(kla,klb)
    kuz=max(kua,kub)
    bandz=klz+kuz+1

!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(sA(LDA,N)) ! LDA not optimal choice
     call DLAG2S(LDA,N,A,LDA,sA,LDA,infoloc1)
     if (LDB/=-1) then
        allocate(sB(LDB,N)) ! LDB not optimal choice
        call DLAG2S(LDB,N,B,LDB,sB,LDB,infoloc1)
     endif
     allocate(Ac(bandz,N,nfact))
  else 
     allocate(Az(bandz,N,nfact))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up SPIKE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call spikeinit(spm,N,max(klz,kuz))
  if (fpm(42)==1) then
     allocate(cworkspike(spm(10)*max(klz,kuz)**2,nfact))
      allocate(cwork(N,M0))
     else
        allocate(zworkspike(spm(10)*max(klz,kuz)**2,nfact))
     end if

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
             Ac(1:bandz,N,id)=CZERO
                do i=1,N
                   Ac(1:bandz,i,id) = -sA(1:bandz,i)*CONE
                   Ac(kua+1,i,id) = Ac(kua+1,i,id) + cmplx(Ze)
                end do
           else !! generalized     
 Ac(1:bandz,N,id)=CZERO
                Ac(kuz-kua+1:kuz+kla+1,1:N,id) = -CONE*sA(1:kla+kua+1,1:N)
                Ac(kuz-kub+1:kuz+klb+1,1:N,id) = Ac(kuz-kub+1:kuz+klb+1,1:N,id) + cmplx(Ze)*sB(1:kub+klb+1,1:N)

              
             end if
             
  call CSPIKE_GBTRF(spm,N,klz,kuz,Ac(1,1,id),bandz,cworkspike(1,id),infoloc2) 
        else ! double precision

           if (LDB==-1) then !! standard
              Az(1:bandz,N,id)=ZZERO
                do i=1,N
                   Az(1:bandz,i,id) = -A(1:bandz,i)*ZONE
                   Az(kua+1,i,id) = Az(kua+1,i,id) + Ze
                end do
           else !! generalized     
              Az(1:bandz,N,id)=ZZERO
                Az(kuz-kua+1:kuz+kla+1,1:N,id) = -ZONE*A(1:kla+kua+1,1:N)
                Az(kuz-kub+1:kuz+klb+1,1:N,id) = Az(kuz-kub+1:kuz+klb+1,1:N,id) + Ze*B(1:kub+klb+1,1:N)

             end if
             
            call ZSPIKE_GBTRF(spm,N,klz,kuz,Az(1,1,id),bandz,zworkspike(1,id),infoloc2)  

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
            call CSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
       
 call ZSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N) 
end if

        if (infoloc1/=0) then
           info=-1
           return
        end if
        
   !     if (infoloc2/=0) then
   !        info=-2
   !        return
   !     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
 id=fpm(33)
        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
            call CSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
       
 call ZSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N) 

           
        end if
        if (infoloc1/=0) then
           info=-1
           return
        end if
        
   !     if (infoloc2/=0) then
   !        info=-2
   !        return
   !     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MKL
        ! With MKL-blas- one should use 
         call mDZGBMM('N','N',n,fpm(25),kla,kua,ZONE,A,kla+kua+1,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#else
         ! not optimal (extra copy)
           call ZGBMM('N','N',n,fpm(25),kla,kua,ZONE,zA,kla+kua+1,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
   
#endif       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(31) !! perform multiplication A^T*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MKL
        ! With MKL-blas- one should use
        call mDZGBMM('T','N',N,fpm(35),kla,kua,ZONE,A,kla+kua+1,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)

#else
        ! not optimal (extra copy)
       call ZGBMM('C','N',N,fpm(35),kla,kua,ZONE,zA,kla+kua+1,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#endif      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
#ifdef MKL
           ! With MKL-blas- one should use
            call mDZGBMM('N','N',N,fpm(25),klb,kub,ZONE,B,klb+kub+1,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#else
            !  ! not optimal (extra copy)
             call ZGBMM('N','N',N,fpm(25),klb,kub,ZONE,zB,klb+kub+1,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
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
            call mDZGBMM('T','N',N,fpm(35),klb,kub,ZONE,B(1,1),klb+kub+1,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#else
           call ZGBMM('C','N',N,fpm(35),klb,kub,ZONE,zB,klb+kub+1,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)   
           
#endif
        end if

     end select
  end do


  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)

  if (fpm(42)==1) then
      deallocate(cworkspike)
     deallocate(cwork)
     deallocate(sA)
      if (LDB/=-1)  deallocate(sB)
     deallocate(Ac)
  else
      deallocate(zworkspike)
     deallocate(Az)
  end if
  
#ifdef MKL
  ! no copy needed
#else
  ! (extra copy)
  deallocate(zA) 
  if (LDB/=-1) deallocate(zB)
#endif       


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
  !
  !  N          (input)        INTEGER: Size system
  !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
  !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=kla+kua+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
  !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=klb+kub+1)
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
  ! Eric Polizzi 2009-2019
  ! ====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
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
  integer :: ijob,infoloc1,infoloc2,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,zAq,zBq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:),allocatable :: cA,cB ! mixed precision
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  integer :: rank,code,nb_procs,NEW_COMM_WORLD
 !! banded conversion
  integer :: mlda,mldb,bandz,klz,kuz
  !! spike 
  integer,dimension(64) :: spm
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::zworkspike
  complex, dimension(:,:),pointer ::cworkspike



  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141335 ! code name
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
       CALL XERBLA( 'ZFEAST_GBGVX', -INFO+100 )
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

    klz=max(kla,klb)
    kuz=max(kua,kub)
    bandz=klz+kuz+1
  
!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(cA(LDA,N)) ! LDA not optimal
     call ZLAG2C(LDA,N,A,LDA,cA,LDA,infoloc1)
     if (LDB/=-1) then
        allocate(cB(LDB,N)) ! LDB not optimal
        call ZLAG2C(LDB,N,B,LDB,cB,LDB,infoloc1)
     endif
     allocate(Ac(bandz,N,nfact))
  else 
     allocate(Az(bandz,N,nfact))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up SPIKE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call spikeinit(spm,N,max(klz,kuz))
  if (fpm(42)==1) then
     allocate(cworkspike(spm(10)*max(klz,kuz)**2,nfact))
      allocate(cwork(N,M0))
     else
        allocate(zworkspike(spm(10)*max(klz,kuz)**2,nfact))
     end if

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
             Ac(1:bandz,N,id)=CZERO
                do i=1,N
                   Ac(1:bandz,i,id) = -cA(1:bandz,i)
                   Ac(kua+1,i,id) = Ac(kua+1,i,id) + cmplx(Ze)
                end do
           else !! generalized     
 Ac(1:bandz,N,id)=CZERO
                Ac(kuz-kua+1:kuz+kla+1,1:N,id) = -cA(1:kla+kua+1,1:N)
                Ac(kuz-kub+1:kuz+klb+1,1:N,id) = Ac(kuz-kub+1:kuz+klb+1,1:N,id) + cmplx(Ze)*cB(1:kub+klb+1,1:N)

              
             end if
             
  call CSPIKE_GBTRF(spm,N,klz,kuz,Ac(1,1,id),bandz,cworkspike(1,id),infoloc2) 
        else ! double precision

           if (LDB==-1) then !! standard
              Az(1:bandz,N,id)=ZZERO
                do i=1,N
                   Az(1:bandz,i,id) = -A(1:bandz,i)
                   Az(kua+1,i,id) = Az(kua+1,i,id) + Ze
                end do
           else !! generalized     
              Az(1:bandz,N,id)=ZZERO
                Az(kuz-kua+1:kuz+kla+1,1:N,id) = -A(1:kla+kua+1,1:N)
                Az(kuz-kub+1:kuz+klb+1,1:N,id) = Az(kuz-kub+1:kuz+klb+1,1:N,id) + Ze*B(1:kub+klb+1,1:N)

             end if
             
            call ZSPIKE_GBTRF(spm,N,klz,kuz,Az(1,1,id),bandz,zworkspike(1,id),infoloc2)  

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
            call CSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
       
 call ZSPIKE_GBTRS(spm,'N',N,klz,kuz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N) 
end if

  if (infoloc1/=0) then
           info=-1
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     case(21) !!solve the linear system (ZeB-A)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id=fpm(33)
          if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
            call CSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
       
 call ZSPIKE_GBTRS(spm,'C',N,klz,kuz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N) 
end if



        if (infoloc1/=0) then
           info=-1
           return
        end if
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call ZGBMM('N','N',n,fpm(25),kla,kua,ZONE,A(1,1),kla+kua+1,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(31) !! perform multiplication A^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call ZGBMM('C','N',N,fpm(35),kla,kua,ZONE,A(1,1),kla+kua+1,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)),N,work(1,fpm(24)),N)
        else
            call ZGBMM('N','N',N,fpm(25),klb,kub,ZONE,B(1,1),klb+kub+1,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(41)!! perform multiplication B^H*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(35),X(1,fpm(34)),N,work(1,fpm(34)),N)
        else
           call ZGBMM('C','N',N,fpm(35),klb,kub,ZONE,B(1,1),klb+kub+1,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)  
        end if

     end select
  end do


  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)

  if (fpm(42)==1) then
      deallocate(cworkspike)
     deallocate(cwork)
     deallocate(cA)
      if (LDB/=-1)  deallocate(cB)
     deallocate(Ac)
  else
      deallocate(zworkspike)
     deallocate(Az)
  end if
 

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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=2*klb+1)
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
  ! Eric Polizzi 2009-2019
  !=====================================================================
  implicit none
  !  include 'f90_noruntime_interface.fi'
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
!! banded conversion
  integer :: mlda,mldb,bandz,klz
  complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
  complex, dimension(:),pointer ::ctmp
  !! spike 
  integer,dimension(64) :: spm
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::zworkspike
  complex, dimension(:,:),pointer ::cworkspike


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141135 ! code name
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


    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF


  
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
       CALL XERBLA( 'ZFEAST_SYGVX', -INFO+100 )
       RETURN
    END IF

  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Format "conversion"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  UPLO2=UPLO
  IF ((UPLO=='F').or.(UPLO=='f')) UPLO2='U' !L will not be consistant  here with ZSBMM 
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

 klz=max(kla,klb)
 bandz=2*klz+1


!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(cA(LDA,N)) !!<< LDA may be overestimated (can be optimized)
     call ZLAG2C(LDA,N,A,LDA,cA,LDA,infoloc1)
     if (LDB/=-1) then
        allocate(cB(LDB,N))!!<< LDB may be overestimated (can be optimized)
        call ZLAG2C(LDB,N,B,LDB,cB,LDB,infoloc1)
     endif
     allocate(Ac(bandz,N,nfact))
      allocate(ctmp(bandz))
  else
     allocate(Az(bandz,N,nfact))
        allocate(ztmp(bandz))
  endif

  if (infoloc1/=0) then
     info=-1
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Set up SPIKE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call spikeinit(spm,N,klz)
  if (fpm(42)==1) then
     allocate(cworkspike(klz*klz*spm(10),nfact))
      allocate(cwork(N,M0))
     else
        allocate(zworkspike(klz*klz*spm(10),nfact))
     end if


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

     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (fpm(42)==1) then ! single precision
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           if (LDB==-1) then !! standard 
!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Ac(klz+1:bandz,1:N,id)=-cA(1:kla+1,1:N)
                   do i=1,N-1
                      s=min(klz,N-i) 
                      call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                      Ac(klz+1,i,id)=-cA(1,i)+cmplx(Ze)
                   enddo
                   Ac(klz+1,N,id)=-cA(1,N)+cmplx(Ze)


                elseif ((UPLO=='U').or.(UPLO=='u')) then
                   Ac(1:klz+1,1:N,id)=-cA(1:klz+1,1:N)
                   do i=1,N-1
                      s=min(klz,N-i)
                      call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                      Ac(klz+1,i,id)=-cA(klz+1,i)+cmplx(Ze)
                   enddo
                   Ac(klz+1,N,id)=-cA(klz+1,N)+cmplx(Ze)


                else !UPLO=F

                   do i=1,N
                      Ac(1:bandz,i,id)=-cA(1:bandz,i)
                      Ac(klz+1,i,id)=Ac(klz+1,i,id)+cmplx(Ze)
                   end do

                end if
              
           else  !! generalized       
!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Ac(kla+1:bandz,1:N,id)=-cA(1:kla+1,1:N)
                      Ac(kla+1:kla+1+klb,1:N,id)=Ac(kla+1:kla+1+klb,1:N,id)+cmplx(Ze)*cB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                      enddo

                   elseif ((UPLO=='U').or.(UPLO=='u')) then
                      Ac(1:kla+1,1:N,id)=-cA(1:kla+1,1:N)
                      Ac(kla+1-klb:kla+1,1:N,id)=Ac(kla+1-klb:kla+1,1:N,id)+cmplx(Ze)*cB(1:klb+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Ac(1:bandz,i,id)=-cA(1:bandz,i)
                         ctmp(1:2*klb+1)=cmplx(Ze)*cB(1:2*klb+1,i)
                         call CAXPY(2*klb+1,CONE,ctmp(1),1,Ac(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Ac(klb+1:bandz,1:N,id)=cmplx(Ze)*cB(1:klb+1,1:N)
                      Ac(klb+1:klb+1+kla,1:N,id)=Ac(klb+1:klb+1+kla,1:N,id)-cA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz+2,i,id),1,Ac(klz,i+1,id),2*klz)
                      enddo

                   elseif  ((UPLO=='U').or.(UPLO=='u')) then
                      Ac(1:klb+1,1:N,id)=cmplx(Ze)*cB(1:klb+1,1:N)
                      Ac(klb+1-kla:klb+1,1:N,id)=Ac(klb+1-kla:klb+1,1:N,id)-cA(1:kla+1,1:N)
                      do i=1,N-1
                         s=min(klz,N-i)
                         call CCOPY(s,Ac(klz,i+1,id),2*klz,Ac(klz+2,i,id),1)
                      enddo

                   else !UPLO=F

                      do i=1,N
                         Ac(1:bandz,i,id)=cB(1:bandz,i)*cmplx(Ze)
                         ctmp(1:2*kla+1)=cA(1:2*kla+1,i)
                         call CAXPY(2*kla+1,-CONE,ctmp(1),1,Ac(klb+1-kla,i,id),1)
                      end do


                   endif

                end if
              
           end if

           

         call CSPIKE_GBTRF(spm,N,klz,klz,Ac(1,1,id),bandz,cworkspike(1,id),infoloc2)    
        

           
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
        else ! double precision
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

           if (LDB==-1) then !! standard
!!!!! Format CONVERSION
                if ((UPLO=='L').or.(UPLO=='l')) then
                   Az(klz+1:bandz,1:N,id)=-A(1:kla+1,1:N)
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
                      Az(1:bandz,i,id)=-A(1:bandz,i)
                      Az(klz+1,i,id)=Az(klz+1,i,id)+Ze
                   end do

                end if

              
           else !! generalized
              !Az(1:N,1:N,id)=Ze*B(1:N,1:N)-A(1:N,1:N)*ZONE 

!!! Format CONVERSION
                if (kla>=klb) then 

                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(kla+1:bandz,1:N,id)=-A(1:kla+1,1:N)
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
                         Az(1:bandz,i,id)=-A(1:bandz,i)
                         ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                         call ZAXPY(2*klb+1,ZONE,ztmp(1),1,Az(kla+1-klb,i,id),1)
                      end do

                   end if


                else ! kla<klb!!!!!!!!!!!!!!!!


                   if ((UPLO=='L').or.(UPLO=='l')) then
                      Az(klb+1:bandz,1:N,id)=Ze*B(1:klb+1,1:N)
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
                         Az(1:bandz,i,id)=B(1:bandz,i)*Ze
                         ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                         call ZAXPY(2*kla+1,-ZONE,ztmp(1),1,Az(klb+1-kla,i,id),1)
                      end do


                   endif
                end if
           end if


         call ZSPIKE_GBTRF(spm,N,klz,klz,Az(1,1,id),bandz,zworkspike(1,id),infoloc2)    
        
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
            call CSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Ac(1,1,id),bandz,cworkspike(1,id),cwork,N) 

           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
       
 call ZSPIKE_GBTRS(spm,'N',N,klz,klz,fpm(23),Az(1,1,id),bandz,zworkspike(1,id),zwork,N)  
        end if
        if (infoloc1/=0) then
           info=-1
           return
        end if
       ! if (infoloc2/=0) then
       !    info=-2
       !    return
       ! end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

      call ZSBMM(UPLO2,n,fpm(25),kla,ZONE,A(1,1),LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (LDB==-1) then ! standard
           call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )
        else
           call DSBMM(UPLO2,n,fpm(25),klb,ZONE,B(1,1),LDB,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
        end if

     end select
  end do


  deallocate(zAq)
  deallocate(zBq)
  deallocate(work)
  deallocate(zwork)
  if (fpm(42)==1) then
     deallocate(cworkspike)
     deallocate(cwork)
     deallocate(ctmp)
     deallocate(cA)
      if (LDB/=-1)  deallocate(cB)
     deallocate(Ac)
  else
      deallocate(zworkspike)
      deallocate(Az)
      deallocate(ztmp)
  end if

  

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
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=121134
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=2*klb+1)
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

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  fpm(30)=141234
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)

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
  ! Eric Polizzi 2009-2019
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

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  fpm(30)=121334
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !
  !  N          (input)        INTEGER: Size system
  !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
  !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=kla+kua+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
  !  kub        (input)        INTEGER: # of superdiagonals within the band of B 
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=klb+kub+1)
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
  ! Eric Polizzi 2009-2019
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

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  fpm(30)=141334
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
  !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=2*klb+1)
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
  ! Eric Polizzi 2009-2019
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

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  fpm(30)=141134
  call feastdefault(fpm,info) 
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
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
  fpm(30)=121133
  call feastdefault(fpm,info)
  if (info/=0) return
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
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
  fpm(30)=141233
  call feastdefault(fpm,info)
  if (info/=0) return
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
  ! Eric Polizzi 2009-2019
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
  fpm(30)=121333
  call feastdefault(fpm,info)
  if (info/=0) return
  call dfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


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
  ! Eric Polizzi 2009-2019
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
  fpm(30)=141333
  call feastdefault(fpm,info)
  if (info/=0) return
  call zfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)


end subroutine zfeast_gbevx



subroutine zfeast_sbevx(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST BANDED INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A GENERAL COMPLEX SYMMETRIC:: BANDEd FORMAT 
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
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
  ! Eric Polizzi 2009-2019
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
  fpm(30)=141133
  call feastdefault(fpm,info)
  if (info/=0) return
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
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

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  double precision,dimension(1) :: B ! dummy
  integer :: LDB,klb

  fpm(30)=121132
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  klb=0
  call dfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
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

  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB,klb

  fpm(30)=141232
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(2)))
  allocate(Wne(fpm(2)))

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  klb=0
  call zfeast_hbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !
  !  N          (input)        INTEGER: Size system
  !  kla        (input)        INTEGER: # of subdiagonals within the band of A
  !  kua        (input)        INTEGER: # of upperdiagonals within the band of A 
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2019
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

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  double precision,dimension(1) :: B ! dummy
  integer :: LDB,klb,kub

  fpm(30)=121332
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  klb=0
  kub=0
  call dfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2019
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

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB,klb,kub

  fpm(30)=141332
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  klb=0
  kub=0
  call zfeast_gbgvx(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
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
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=2*kla+1)
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
  ! Eric Polizzi 2009-2019
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

  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: Zne,Wne 
  complex(kind=(kind(1.0d0))),dimension(1) :: B ! dummy
  integer :: LDB,klb

  fpm(30)=141132
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  LDB=-1 ! B is a dummy- option for standard eigenvalue problem
  klb=0
  call zfeast_sbgvx(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_sbev





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! AUXILIARY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MKL

subroutine mDZGBMM(TRANSA,TRANSB,n,rhs,kl,ku,alpha,A,LDA,B,LDB,beta,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose
  !  =======
  !  Perform the multiplication:
  !                              C=alpha*op(A)*op(B)+beta*C
  !  where alpha and beta are double complex scalars, x and y are double complex vectors and A is an
  !  n-by-n (non-symmetric) double precision band matrix, with kl sub-diagonals and ku super-diagonals,
  !  B and C are complex double precision n-by-rhs matrices.
  !
  !  It is a  Blocked version (BLAS3) of the routine ZHBMV which allow blocking with 
  !  multiple right-hand sides  (with  generalization to structural symetric case.. not necessarely Hermitian).    
  !
  !
  !  Arguments
  !  =========
  !  TRANSA (input) CHARACTER
  !          'N'    A is normal
  !          'T'    A is transposed
  !          'C'    A is complex conjugate tranposed
  !  TRANSB (input) CHARACTER
  !          'N'    B is normal
  !          'T'    B is transposed
  !          'C'    B is complex conjugate tranposed
  !  N      (input) INTEGER
  !          The number of row and columns of the matrix A.  
  !          The number of rows of matrices B and C.   N >= 0.
  !  rhs    (input) INTEGER
  !          The number of columns of matrices B and C
  !  KU      (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  KL      (input) INTEGER 
  !                 bandwidth = kl+ku+1
  !          The number of subdiagonals within the band of A.  KL >= 0.
  ! alpha   (input) COMPLEX DOUBLE PRECISION
  !  A      (input) DOUBLE PRECISION  
  !          array, dimension (LDA,N); the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! B       (input) COMPLEX DOUBLE PRECISION - 2D array
  ! LDB     (input) INTEGER      
  !          The leading dimension of the array B.  LDB >= N.
  ! beta    (input) COMPLEX DOUBLE PRECISION
  ! C       (input/output) COMPLEX DOUBLE PRECISION - 2D array
  !         On exit, contains the solutions 
  ! LDC     (input) INTEGER      
  !          The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! James Kestyn 2015
  !======================================================================
  implicit none
  CHARACTER :: TRANSA, TRANSB
  INTEGER :: N, RHS, KU,KL,LDA,LDB,LDC
  COMPLEX(kind=kind(1.0d0)):: alpha,beta
  DOUBLE PRECISION,dimension(LDA,*) ::  A
  COMPLEX(kind=kind(1.0d0)),dimension(LDB,*) ::  B
  COMPLEX(kind=kind(1.0d0)),dimension(LDC,*) ::  C
!!!!!!!!!!!!!!

  integer :: info_alloc
  logical :: low,up,full

  integer :: bl, i,j, k , kd, tmp, imin,jmin,iminb,iminc,imin2,jmin2
  character(len=1) :: UPLO2,DIAG,SIDE
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
 ! COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO)

  if( TRANSA == 'N' .or. TRANSA == 'n') then

     do j=1,kl !first kl rows of A
        tmp = (ku)+j
        imin= (ku)+j
        jmin= 1
        call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(1,1), LDB, beta, C(j,1), LDC)
     enddo

     do j=1,N-ku-kl !Middle rows of A (no offset)
        tmp = ku+kl+1
        imin = ku+kl+1
        call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,j), LDA-1, B(j,1), LDB, beta, C(kl+j,1), LDC)
     enddo

     do j=1,ku !last ku rows of A
        tmp = (ku+kl+1)-j
        imin= kl+ku+1
        jmin= (N-kl-ku)+j!(ku+2)-i
        call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(N-kl-ku+j,1), LDB, beta, C(N-ku+j,1), LDC)
     enddo

  else

     do j=1,ku !first ku rows of A^C
        tmp = (kl)+j
        imin= (ku+2)-j
        jmin= j
        call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(1,1), LDB, beta, C(jmin,1), LDC)
     enddo


     do j=1,N-ku-kl
        tmp = ku+kl+1
        imin= 1
        jmin= ku+j
        call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(j,1), LDB, beta, C(jmin,1), LDC)
     enddo

     do j=1,kl !last kl rows of A^C
        tmp = (ku+kl+1)-j
        imin= 1
        jmin= (N-kl)+j!(ku+2)-i
        call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(N-ku-kl+j,1), LDB, beta, C(jmin,1), LDC)
     enddo

  endif

end subroutine mDZGBMM


#endif



