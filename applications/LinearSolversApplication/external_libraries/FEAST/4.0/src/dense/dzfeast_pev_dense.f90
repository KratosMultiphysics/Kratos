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

!{D,Z}FEAST_{SY,HE,GE}{PEV}{X}


!111111111111!!!!!! EXPERT ROUTINE  (CORE)
! dfeast_sypevx
! zfeast_hepevx
! dfeast_gepevx
! zfeast_gepevx
! zfeast_sypevx

!222222222222!!!!  DEFAULT ROUTINES (Wrappers to expert)
! dfeast_sypev
! zfeast_hepev
! dfeast_gepev
! zfeast_gepev
! zfeast_sypev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!11111111111111111111111111111111111111111111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dfeast_sypevx(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL SYMMETRIC DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0): 
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
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,dmax
  double precision,dimension(LDA,N,*):: A
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
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,work,zBq,zAq 
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az,zA
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  real,dimension(:,:,:),allocatable :: sA ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=121127 ! code name
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
     !ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     !   INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_SYPEVX', -INFO+100 )
     RETURN
  END IF

  infoloc1=0
  infoloc2=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! :-( If DZGEMM  not available for case 30, we need to make a complex copy of the matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
#ifdef MKL
  ! copy needed if UPLO/='F' (no symmetric equivalent to DZGEMM)
  IF (.not.((UPLO=='F').or.(UPLO=='f'))) then
     allocate(zA(N,N,dmax+1))
     do i=1,dmax+1
        call ZLACP2( 'F', N, N,A(1,1,i) , LDA, zA(1,1,i), N )
     enddo
  end IF
#else
  ! (extra copy)   
  allocate(zA(N,N,dmax+1))
  do i=1,dmax+1
     call ZLACP2( 'F', N, N,A(1,1,i) , LDA, zA(1,1,i), N )
  enddo
#endif       

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
     allocate(sA(N,N,dmax+1))
     do i=1,dmax+1
        call DLAG2S(N,N,A(1,1,i),LDA,sA(1,1,i),N,infoloc1)
     enddo
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

  allocate(zAq(M0*dmax,M0*dmax))
  allocate(zBq(M0*dmax,M0*dmax))
  allocate(work(N,M0))
  allocate(zwork(N,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_srcipevx(ijob,dmax,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)  !!<< do dfeast_...  
     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(10) !! !! form P(Ze) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        id=fpm(33) !! id of factorization (for fpm(10) flag) 

        if (fpm(42)==1) then ! single precision
           Ac(1:N,1:N,id)=sA(1:N,1:N,1)*SONE
           do i=2,dmax+1
              Ac(1:N,1:N,id)=Ac(1:N,1:N,id)+sA(1:N,1:N,i)*cmplx(Ze)**(i-1)
           enddo

           call CSYTRF(UPLO2,N,Ac(1,1,id),N,IPIVloc(1,id),cwork,N*M0,INFOloc2)     
        else ! double precision
           Az(1:N,1:N,id)=A(1:N,1:N,1)*ZONE
           do i=2,dmax+1
              Az(1:N,1:N,id)=Az(1:N,1:N,id)+A(1:N,1:N,i)*Ze**(i-1)
           enddo

           call ZSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),zwork,N*M0,INFOloc2)     
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system  P(Ze)x=workc(1:N,1:fpm(23)) result in to workc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

#ifdef MKL
        ! copy needed if UPLO/='F' (no symmetric equivalent to DZGEMM)
        IF (.not.((UPLO=='F').or.(UPLO=='f'))) then
           call ZSYMM ('L', UPLO2, N, fpm(25), ZONE, zA(1,1,fpm(57)), N, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)
        else ! full format (no extra copy needed)
           call DZGEMM('T','N',N,fpm(25),N,ZONE,A(1,1,fpm(57)),LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
        end IF
#else
        ! (extra copy)   
        call ZSYMM ('L', UPLO2, N, fpm(25), ZONE, zA(1,1,fpm(57)), N, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)
#endif       

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
     deallocate(Ac)
  else
     deallocate(Az)
  end if

#ifdef MKL
  IF (.not.((UPLO=='F').or.(UPLO=='f'))) deallocate(zA)
#else 
  deallocate(zA)
#endif       

end subroutine dfeast_sypevx



subroutine zfeast_hepevx(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX HERMITIAN- DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input/output)        COMPLEX DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A}- On exit: full format if matrices were lower or upper 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
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
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !           
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,dmax
  complex(kind=(kind(1.0d0))),dimension(LDA,N,*):: A
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc1,infoloc2,i,s,j
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,zAq,zBq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az 
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:,:),allocatable :: cA ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO),ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141257 ! code name
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
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_HEPEVX', -INFO+100 )
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
     do j=1,dmax+1
        do i=1,N-1
           s=N-i
           call ZCOPY(s,A(i+1,i,j),1,A(i,i+1,j),LDA)
           call ZLACGV(s, A(i,i+1,j), LDA )
        enddo
     enddo

!!!!
  elseif ((UPLO=='U').or.(UPLO=='u')) then
!!!!
     do j=1,dmax+1
        do i=1,N-1
           s=N-i
           call ZCOPY(s,A(i,i+1,j),LDA,A(i+1,i,j),1)
           call ZLACGV(s, A(i+1,i,j), 1 )
        enddo
     enddo

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
     do i=rank+1,fpm(8),nb_procs
        nfact=nfact+1
     end do
  endif

!!!! mixed precision set-up
  if (fpm(42)==1) then ! copy single precision
     allocate(cA(N,N,dmax+1))
     do i=1,dmax+1
        call ZLAG2C(N,N,A(1,1,i),LDA,cA(1,1,i),N,infoloc1)
     enddo
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


  allocate(zAq(M0*dmax,M0*dmax))
  allocate(zBq(M0*dmax,M0*dmax))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_grcipevx(ijob,dmax,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(10) !! !! form P(Ze) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag) 

        if (fpm(42)==1) then ! single precision
           Ac(1:N,1:N,id)=cA(1:N,1:N,1)
           do i=2,dmax+1
              Ac(1:N,1:N,id)=Ac(1:N,1:N,id)+cA(1:N,1:N,i)*cmplx(Ze)**(i-1)
           enddo

           call CGETRF(N,N,Ac(1,1,id),N,IPIVloc(1,id),INFOloc2)

        else ! double precision
           Az(1:N,1:N,id)=A(1:N,1:N,1)
           do i=2,dmax+1
              Az(1:N,1:N,id)=Az(1:N,1:N,id)+A(1:N,1:N,i)*Ze**(i-1)
           enddo

           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc2) 
        end if


        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system  P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
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
     case(21) !!solve the linear system  P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        call ZHEMM ('L', UPLO2, N, fpm(25), ZONE, A(1,1,fpm(57)), LDA, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(31) !! perform multiplication A^H(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        call ZHEMM ('L', UPLO2, N, fpm(35), ZONE, A(1,1,fpm(57)), LDA, X(1,fpm(34)), N, ZZERO,work(1,fpm(34)), N)

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
     deallocate(Ac)
  else
     deallocate(Az)
  end if


end subroutine zfeast_hepevx





subroutine dfeast_gepevx(dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
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
  !  include 'f90_noruntime_interface.fi'
  integer :: N,LDA,dmax
  double precision,dimension(LDA,N,*):: A
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
  integer :: ijob,infoloc1,infoloc2,i,s,j
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::work,zwork,zAq,zBq
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: zA
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  real,dimension(:,:,:),allocatable :: sA ! mixed precision
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
     fpm(30)=121357 ! code name
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
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_GEPEVX', -INFO+100 )
     RETURN
  END IF


  infoloc1=0
  infoloc2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! :-( If DZGEMM  not available for case 30-40, we need to make a complex copy of all the matrices 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
#ifdef MKL
  ! everything fine no copy needed
#else
  ! (extra copy)
  allocate(zA(N,N,dmax+1))
  do j=1,dmax+1
     call ZLACP2( 'F', N, N,A(1,1,j) , LDA, zA(1,1,j), N )
  end do
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
     allocate(sA(N,N,dmax+1))
     do i=1,dmax+1
        call DLAG2S(N,N,A(1,1,i),LDA,sA(1,1,i),N,infoloc1)
     end do
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

  allocate(zAq(M0*dmax,M0*dmax))
  allocate(zBq(M0*dmax,M0*dmax))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ijob=-1 ! initialization
  do while (ijob/=0) 
     call zfeast_grcipevx(ijob,dmax,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     case(10) !! !! form P(Ze) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33) !! id of factorization (for fpm(10) flag) 

        if (fpm(42)==1) then ! single precision
           Ac(1:N,1:N,id)=sA(1:N,1:N,1)*CONE
           do i=2,dmax+1
              Ac(1:N,1:N,id)=Ac(1:N,1:N,id)+sA(1:N,1:N,i)*cmplx(Ze)**(i-1)
           enddo

           call CGETRF(N,N,Ac(1,1,id),N,IPIVloc(1,id),INFOloc2)

        else ! double precision
           Az(1:N,1:N,id)=A(1:N,1:N,1)*ZONE
           do i=2,dmax+1
              Az(1:N,1:N,id)=Az(1:N,1:N,id)+A(1:N,1:N,i)*Ze**(i-1)
           enddo

           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc2) 
        end if



        if (infoloc2/=0) then
           info=-2
           return
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system  P(Ze)x=zwork(1:N,1:fpm(23)) result in to zwork
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
     case(21) !!solve the linear system  P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
#ifdef MKL
        ! With MKL-blas- one should use 
        call DZGEMM('N','N',N,fpm(25),N,ZONE,A(1,1,fpm(57)),LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#else
        ! not optimal (extra copy)       
        call ZGEMM('N','N',N,fpm(25),N,ZONE,zA(1,1,fpm(57)),N,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)
#endif       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(31) !! perform multiplication A^H(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
#ifdef MKL
        ! With MKL-blas- one should use
        call DZGEMM('T','N',N,fpm(35),N,ZONE,A(1,1,fpm(57)),LDA,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#else
        ! not optimal (extra copy)
        call ZGEMM('C','N',N,fpm(35),N,ZONE,zA(1,1,fpm(57)),N,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)
#endif      

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
     deallocate(Ac)
  else
     deallocate(Az)
  end if

#ifdef MKL
  ! no copy needed
#else
  ! (extra copy)
  deallocate(zA) 
#endif       


end subroutine dfeast_gepevx






subroutine zfeast_gepevx(dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
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
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  integer :: N,LDA,dmax
  complex(kind=(kind(1.0d0))),dimension(LDA,N,*):: A
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
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,work,zBq,zAq 
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:,:),allocatable :: cA ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141357 ! code name
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

  IF ( N<=0 ) THEN
     INFO = -102
  ELSE IF(LDA<N ) THEN
     INFO = -104
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_GEPEVX', -INFO+100 )
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
     allocate(cA(N,N,dmax+1))
     do i=1,dmax+1
        call ZLAG2C(N,N,A(1,1,i),LDA,cA(1,1,i),N,infoloc1)
     enddo
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

  allocate(zAq(M0*dmax,M0*dmax))
  allocate(zBq(M0*dmax,M0*dmax))
  if (fpm(15)==0) then
     allocate(work(N,2*M0))
  else
     allocate(work(N,M0))
  end if
  allocate(zwork(N,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_grcipevx(ijob,dmax,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(10) !! !! form P(Ze) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        id=fpm(33) !! id of factorization (for fpm(10) flag) 

        if (fpm(42)==1) then ! single precision
           Ac(1:N,1:N,id)=cA(1:N,1:N,1)
           do i=2,dmax+1
              Ac(1:N,1:N,id)=Ac(1:N,1:N,id)+cA(1:N,1:N,i)*cmplx(Ze)**(i-1)
           enddo

           call CGETRF(N,N,Ac(1,1,id),N,IPIVloc(1,id),INFOloc2)

        else ! double precision
           Az(1:N,1:N,id)=A(1:N,1:N,1)
           do i=2,dmax+1
              Az(1:N,1:N,id)=Az(1:N,1:N,id)+A(1:N,1:N,i)*Ze**(i-1)
           enddo

           call ZGETRF(N,N,Az(1,1,id),N,IPIVloc(1,id),INFOloc2) 
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system  P(Ze)x=workc(1:N,1:fpm(23)) result in to workc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        id=fpm(33)

        if (fpm(42)==1) then ! single precision solve
           call ZLAG2C(N, fpm(23), zwork, N, cwork, N, infoloc1)
           call CGETRS( 'N', N, fpm(23), Ac(1,1,id), N, IPIVloc(1,id), cwork, N, INFOloc2 )
           call CLAG2Z(N, fpm(23), cwork, N, zwork, N, infoloc1)
        else
           call ZGETRS( 'N', N, fpm(23), Az(1,1,id), N, IPIVloc(1,id), zwork, N, INFOloc2)          
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
     case(21) !!solve the linear system  P(Ze)^Hx=zwork(1:N,1:fpm(23)) result in to zwork
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        call ZGEMM('N','N',N,fpm(25),N,ZONE,A(1,1,fpm(57)),LDA,X(1,fpm(24)),N,ZZERO,work(1,fpm(24)),N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     case(31) !! perform multiplication A^H(fpm(57))*x(1:N,fpm(34):fpm(34)+fpm(35)-1) result in work(1:N,fpm(34)+fpm(35)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        call ZGEMM('C','N',N,fpm(35),N,ZONE,A(1,1,fpm(57)),LDA,X(1,fpm(34)),N,ZZERO,work(1,fpm(34)),N)



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
     deallocate(Ac)
  else
     deallocate(Az)
  end if

end subroutine zfeast_gepevx







subroutine zfeast_sypevx(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX SYMMETRIC DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  !include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,dmax
  complex(kind=(kind(1.0d0))),dimension(LDA,N,*):: A
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
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable ::zwork,work,zBq,zAq 
  complex(kind=(kind(1.0d0))), dimension(:,:,:),allocatable :: Az
  complex, dimension(:,:,:),allocatable :: Ac ! mixed precision
  complex, dimension(:,:),allocatable ::cwork ! mixed precision
  complex,dimension(:,:,:),allocatable :: cA ! mixed precision
  integer, dimension(:,:),allocatable ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
  real,parameter :: SONE=1.0d0,SZERO=0.0d0
  complex,parameter :: CONE=(SONE,SZERO),CZERO=(SZERO,SZERO)
  integer :: nfact,id
  character(len=1) :: UPLO2
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  info=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(30)==-111) then ! this routine is called directly
     fpm(30)=141127 ! code name
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
     !ELSE IF((LDB<N ).and.(LDB/=-1)) THEN
     !   INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_SYPEVX', -INFO+100 )
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
     allocate(cA(N,N,dmax+1))
     do i=1,dmax+1
        call ZLAG2C(N,N,A(1,1,i),LDA,cA(1,1,i),N,infoloc1)
     enddo
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

  allocate(zAq(M0*dmax,M0*dmax))
  allocate(zBq(M0*dmax,M0*dmax))
  allocate(work(N,M0))
  allocate(zwork(N,M0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! FEAST-RCI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ijob=-1 ! initialization 
  do while (ijob/=0)
     call zfeast_srcipevx(ijob,dmax,N,Ze,work,zwork,zAq,zBq,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)    

     select case(ijob)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(10) !! !! form P(Ze) and factorize preconditioner if any
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        id=fpm(33) !! id of factorization (for fpm(10) flag) 

        if (fpm(42)==1) then ! single precision
           Ac(1:N,1:N,id)=cA(1:N,1:N,1)
           do i=2,dmax+1
              Ac(1:N,1:N,id)=Ac(1:N,1:N,id)+cA(1:N,1:N,i)*cmplx(Ze)**(i-1)
           enddo

           call CSYTRF(UPLO2,N,Ac(1,1,id),N,IPIVloc(1,id),cwork,N*M0,INFOloc2)     
        else ! double precision
           Az(1:N,1:N,id)=A(1:N,1:N,1)
           do i=2,dmax+1
              Az(1:N,1:N,id)=Az(1:N,1:N,id)+A(1:N,1:N,i)*Ze**(i-1)
           enddo

           call ZSYTRF(UPLO2,N,Az(1,1,id),N,IPIVloc(1,id),zwork,N*M0,INFOloc2)     
        end if

        if (infoloc2/=0) then
           info=-2
           return
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case(11) !!solve the linear system  P(Ze)x=workc(1:N,1:fpm(23)) result in to workc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
     case(30) !! perform multiplication A(fpm(57))*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        call ZSYMM ('L', UPLO2, N, fpm(25), ZONE, A(1,1,fpm(57)), LDA, X(1,fpm(24)), N, ZZERO,work(1,fpm(24)), N)

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
     deallocate(Ac)
  else
     deallocate(Az)
  end if

end subroutine zfeast_sypevx







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!22222222222222222222222222222222222222222222222222222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dfeast_sypev(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL SYMMETRIC DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0): 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)      
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,LDA,dmax
  double precision,dimension(LDA,N,*):: A
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
  fpm(30)=121126
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_sypevx(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(Zne)
  deallocate(Wne)

end subroutine dfeast_sypev





subroutine zfeast_hepev(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX HERMITIAN- DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input/output)        COMPLEX DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A}- On exit: full format if matrices were lower or upper 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !        
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,LDA,dmax
  complex(kind=(kind(1.0d0))),dimension(LDA,N,*):: A
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

  fpm(30)=141256
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_hepevx(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(Zne)
  deallocate(Wne)

end subroutine zfeast_hepev




subroutine dfeast_gepev(dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::REAL DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
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
  integer :: N,LDA,dmax
  double precision,dimension(LDA,N,*):: A
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
  fpm(30)=121356
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call dfeast_gepevx(dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine dfeast_gepev



subroutine zfeast_gepev(dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !  X          (input/output) COMPLEX DOUBLE PRECISION(N,M0) or (N,2*M0) to store left eigenvectors if fpm(15)=0: 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                                                    
  !  info       (output)       INTEGER: Error handling (0: successful exit)    
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  integer :: N,LDA,dmax
  complex(kind=(kind(1.0d0))),dimension(LDA,N,*):: A
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
  fpm(30)=141356
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_gepevx(dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)
  deallocate(Zne)
  deallocate(Wne)
end subroutine zfeast_gepev





subroutine zfeast_sypev(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the polynomial  P(E)x=0 eigenvalue problem
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !  with {A}_{0..dmax} stored as A[1] to A[dmax+1] ::COMPLEX SYMMETRIC DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N,dmax+1):  Matrices {A} 
  !  LDA        (input)        INTEGER: Leading dimension of any matrix A (LDA>=N)
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
  !=====================================================================
  ! Eric Polizzi 2019
  !=====================================================================
  implicit none
  character(len=1) :: UPLO
  integer :: N,LDA,dmax
  complex(kind=(kind(1.0d0))),dimension(LDA,N,*):: A
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
  fpm(30)=141126
  call feastdefault(fpm,info)
  if (info/=0) return
  allocate(Zne(fpm(8)))
  allocate(Wne(fpm(8)))
  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
  call zfeast_sypevx(UPLO,dmax,N,A,LDA,fpm,epsout,loop,Emid,r,M0,E,X,mode,res,info,Zne,Wne)

  deallocate(Zne)
  deallocate(Wne)

end subroutine zfeast_sypev






