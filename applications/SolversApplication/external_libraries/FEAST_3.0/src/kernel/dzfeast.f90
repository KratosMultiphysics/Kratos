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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST KERNEL - REVERSE COMMUNICATION INTERFACES !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! List of routines:
!-------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Symmetric (real) /Hermitian (complex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>double
! dfeast_srcix   ! Double Real Symmetric - expert
! dfeast_srci    ! Double Real Symmetric
!
! zfeast_hrcix   ! Double Complex Hermitian - expert
! zfeast_hrci    ! Double Complex Hermitian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Non-Symmetric (real) /Non-Hermitian (complex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> double
! zfeast_grcix   ! Double Complex General - expert    ==>Kernel routine<==
! zfeast_grci    ! Double Complex General
! 
! zfeast_srcix   ! Double Complex Symmetric - expert
! zfeast_srci    ! Double Complex Symmetric
!
! dfeast_grcix   ! Double Real General - expert
! dfeast_grci    ! Double Real General


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RCI ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) - Includes option for custom integration nodes/weight
  !  Solve generalized Aq=lambda Bq eigenvalue problems
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE (or B Identity) 
  !  DOUBLE PRECISION version 
  ! 
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) REAL DOUBLE PRECISION (N,M0)   :  Workspace 
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):  Workspace 
  !  Aq,Sq      (input/output) REAL DOUBLE PRECISION (M0,M0)  :  Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                               
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
  !                
  !                            Expert comment: if fpm(29)=1 Zne, Wne have already been generated using default contour   
  !===========================================================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  double precision, dimension(N,*) ::work
  complex(kind=(kind(1.0d0))),dimension(N,*):: workc
  integer,dimension(*) :: fpm
  double precision,dimension(M0,*):: Aq,Sq
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  double precision,dimension(*)  :: lambda
  double precision,dimension(N,*):: q
  integer :: mode
  double precision,dimension(*) :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  !! parameters
  double precision, Parameter :: pi=3.1415926535897932d0
  double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZEROC=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0
  integer(8),parameter :: fout =6
  !! variable for FEAST
  integer :: i,e,j,k
  integer,dimension(4) :: iseed
  double precision :: theta,r,Emid
  complex(kind=(kind(1.0d0))) :: zxe,zwe,aux
  double precision, dimension(:,:),pointer :: Sqo
  integer, dimension(:),pointer :: fpm_default
  logical :: testconv
  double precision :: trace
  !! Lapack variable (reduced system)
  character(len=1) :: JOBZ,UPLO
  double precision, dimension(:),pointer :: work_loc,work_loc2
  integer :: lwork_loc,info_lap,infoloc
  !! MPI compatibility variables
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

  rank=0
  nb_procs=1
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
  if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
#endif
  !---------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then 
     info=0 ! default value
     if (fpm(1)==1) then
        call wwrite_n(fout)
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout)
        call wwrite_s(fout, 'Routine DFEAST_S{}{}')
        if (fpm(29)==0) call wwrite_s(fout,'X') 
        call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
        call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
        call wwrite_n(fout)

        call wallocate_1i(fpm_default,64,infoloc)
        call feastinit(fpm_default)
        do i=1,19
           if (fpm(i)/=fpm_default(i)) then
              call wwrite_s(fout, '   fpm(')
              call wwrite_i(fout, i)
              call wwrite_s(fout, ')=')
              call wwrite_i(fout, fpm(i))
              call wwrite_n(fout)
           endif
        enddo
        call wdeallocate_1i(fpm_default)

        call wwrite_s(fout, 'Search interval [')
        call wwrite_d(fout,Emin)
        call wwrite_s(fout, '; ')
        call wwrite_d(fout,Emax)
        call wwrite_s(fout, ']')
        call wwrite_n(fout)
     end if
     call check_feast_fpm_input(fpm,info)
     call dcheck_feast_srci_input(Emin,Emax,M0,N,info)

     if (info/=0) fpm(21)=100 ! The End
  
!!!!!!!!!!!!!!!!
     IF (info==0) then
        fpm(22)=fpm(2) ! only one half contour necessary
        loop=0
        if (fpm(1)==1) then
           call wwrite_s(fout, 'Size system    ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,N)
           call wwrite_n(fout)
           call wwrite_s(fout, 'Size subspace  ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,M0)
           call wwrite_n(fout)
           call wwrite_s(fout, '#Linear systems')  
           call wwrite_t(fout) 
           call wwrite_i(fout,fpm(22))
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
           call wwrite_s(fout, '#Loop | #Eig  |       Trace           |     Error-Trace       |     Max-Residual')  
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
        endif
        fpm(23)=min(M0,N) ! 'current M0' size (global value)
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors (option 2 in DLARNV)
           iseed=(/56,890,3456,2333/)
           ! copy work into q for multiplication by B matrix, if stochastic estimate is "on" 
           if (fpm(14)==2) then
              call DLARNV(2,iseed,N*fpm(23),q(1,1)) 
           else
              call DLARNV(2,iseed,N*fpm(23),work(1,1)) 
           endif
        end if

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=DZERO 
#endif
           !------------------------------------------
           ijob=40 !! B*q=>work
           return
        end if
     end IF ! info=0
!!!!!!!!!!!!!!
  end if   !ijob=-1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==1) then !! we initialize a new contour integration
     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(5)==1).or.(fpm(14)==2)) then        
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
     endif

     if ((fpm(29)==1).and.(fpm(16)==2)) then ! Zolotarev 
        if ((loop>0).and.(fpm(23)/nb_procs>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !------------------------------------------------------------------------

     if ((fpm(29)==1).and.(fpm(16)==2)) then !! Zolotarev
        call dset_feast_zolotarev(fpm(2),0,zxe,zwe)
        !   r=(Emax-Emin)/2.0d0
        if ((loop==0).and.(fpm(5)==0)) then 
        !   !q(1:N,fpm(23))=-dble(zwe)*r*work(1:N,1:fpm(23)) 
            q(1:N,1:fpm(23))=DZERO
        else
            if (rank==0) then
             q(1:N,1:fpm(23))=-dble(zwe)*q(1:N,1:fpm(23))
            else 
             q(1:N,1:fpm(23))=DZERO
            endif
        end if
     else
        q(1:N,1:fpm(23))=DZERO
     end if

     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     IF (info==0) then !! will end up checking info errors returned by FEAST drivers
        do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the contour 

           if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
              Ze=Zne(e)
              fpm(20)=e-rank
              ijob=10 ! for fact
              if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
           endif

           if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v 
              call ZLACP2( 'F', N, fpm(23),work , N, workc, N )
              ijob=11 ! for solve
              return
           endif

           if (ijob==11) then 
              !! summation              
              aux=2.0d0*Wne(e)            
              !! aux has been multiplied by 2 to account for the 2 half-contours
              q(1:N,1:fpm(23))=q(1:N,1:fpm(23))-dble(aux*workc(1:N,1:fpm(23))) 

              ijob=-2 ! just for identification
           end if
        end do
     end IF  !! info=0

     !------------------------------------------------
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------------  
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(21)=4 
        !------------------------------------------------
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if
  end IF     ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic estimates
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only q is returned with stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     call DLARNV(2,iseed,N*fpm(23),work(1,1)) 
     call DGEMM('T','N',fpm(23),fpm(23),N,-DONE,work,N,q,N,DZERO,Aq,M0) ! projection
     call DGEMM('T','N',fpm(23),fpm(23),N,DONE,work,N,work,N,DZERO,Sq,M0) ! normalization
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+Aq(i,i)/Sq(i,i)
        res(i)=abs(theta*(N/(DONE*i)))
     enddo
     mode=int(real(res(fpm(23))))+1
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Bq xq     with Aq=Q^TAQ Bq=Q^TBQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form  Bq=> Bq=Q^T B Q
  if (fpm(21)==4) then 
#ifdef MPI 
     work(1:N,1:fpm(23))=DZERO
#endif
     fpm(21)=5 ! preparing reenty
     ijob=40 
     return! mat-vec B*q => work
  end if


  if (fpm(21)==5) then
     !----------------------------------------------
#ifdef MPI
     Sq(1:M0,1:fpm(23))=DZERO
#endif
     !-------------------------------------------------

     call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Sq(1,fpm(24)),M0) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     mode=0
     if (loop>0) then
        do i=fpm(24),fpm(24)+fpm(25)-1
           if (res(i)==DONE) then ! indicator for being inside the search interval
              if (abs(sqrt(abs(Sq(i,i)))-lambda(i))/sqrt(abs(Sq(i,i)))<0.1d0) mode=mode+1
           endif
        enddo
     endif
!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,NEW_COMM_WORLD,code)
#endif
!!!!!!!!!!!!!!!!!!
     fpm(21)=6
  end if


!!!!!!!!!form Aq=> Aq=Q^T A Q 
  if (fpm(21)==6) then
     !------------------------------------------------
#ifdef MPI 
     work(1:N,1:fpm(23))=DZERO
#endif
     fpm(21)=7 ! preparing reentry
     ijob=30 
     return  ! mat-vec A*q => work
  endif


  if (fpm(21)==7) then 
     !------------------------------------------------
#ifdef MPI 
     Aq(1:M0,1:fpm(23))=DZERO
#endif
     !-------------------------------------------------
     call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Aq(1,fpm(24)),M0) !q^tAq
     fpm(21)=8
  endif


  if (fpm(21)==8) then

     !---------------------------------------- !(Aq,Sq known to all processors) 
#ifdef MPI
     if (fpm(23)/nb_procs>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,Aq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,Sq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !---------------------------------------

     !---------------------------------------
     if (fpm(12)==1) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=50
        return
     endif
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==8) then
     ! if using FEAST-MPI ==> solve on a single proc
     if (rank==0) then
        JOBZ='V'
        UPLO='L'
        info_lap=1 ! initialization
        i=1
        LWORK_LOC=3*fpm(23)-1 !! for lapack eig reduced system
        call wallocate_1d(WORK_LOC,LWORK_LOC,infoloc)
        if (infoloc/=0) info=-1


        do while ((info_lap/=0).and.(info==0))
           i=i+1
           if (i==10) info=-3 ! arbitrary maximum
           call wallocate_2d(Sqo,fpm(23),fpm(23),infoloc)
           if (infoloc/=0) info=-1
           call DLACPY( 'F', fpm(23), fpm(23),Sq, M0, Sqo, fpm(23) )
           call DSYGV(1, JOBZ, UPLO, fpm(23),Aq,M0,Sqo,fpm(23),lambda,work_loc,Lwork_loc,INFO_lap)

           if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
           if (info_lap>fpm(23)) then !! Sqo is not spd (a posteriori resize subspace)
              fpm(23)=info_lap-fpm(23)-1 
              if (fpm(1)==1) then
                 call wwrite_s(fout, 'Resize subspace')  
                 call wwrite_t(fout) 
                 call wwrite_i(fout,fpm(23))
                 call wwrite_n(fout)
              end if
           end if
           call wdeallocate_2d(Sqo)
        end do
        call wdeallocate_1d(work_loc) 
     end if !(rank 0)
     !-------------------------------- !(info common to all processors)
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(25)=fpm(23)!! current M0 size (by default) -- global
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------!(Aq==> vectors, lambda and fpm(23), known to all processors) 
#ifdef MPI 
        call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(Aq,fpm(23)*M0,MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of current M0
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------
        fpm(21)=9
     end if
  end if !! fpm(21)=8



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors X=Qxq  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     call DLACPY( 'F', N, fpm(23),q , N, work, N ) 
     !! option - non shifted
#ifdef MPI 
     q(1:N,1:fpm(23))=DZERO
#endif
     call DGEMM('N','N',N,fpm(25),fpm(23),DONE,work(1,1),N,Aq(1,fpm(24)),M0,DZERO,q(1,fpm(24)),N)
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     fpm(21)=10 ! preparing reentry
     ijob=30
     return  ! mat-vec A*q => work
  endif

  if (fpm(21)==10) then
     call ZLACP2( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N )
     fpm(21)=11 ! preparing reentry
     ijob=40
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=DZERO  !! work is also needed for outer-loop if any 
#endif
     !----------------------------------------
     return  ! mat-vec S*q => work 
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !-------- Absolute residual
     do i=fpm(24),fpm(24)+fpm(25)-1
           res(i)=sum(abs(dble(workc(1:N,i))-lambda(i)*work(1:N,i)))/sum(abs(max(abs(Emin),abs(Emax))*(work(1:N,i)))) 
     end do
     !----------------------------------------
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (fpm(21)==11) then
     k=0
     trace=DZERO
     theta=DZERO
!!! count how many eigenvalues have converged + trace + max residual
     do i=1,fpm(23)
        if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
           k=k+1 ! number of eigenvalues (could include spurious)
           trace=trace+lambda(i)
           if (res(i)>theta) theta=res(i) ! max residual 
        endif
     enddo

!!!!!!!! remove spurious if any       
     if ((mode==0).and.(k>0)) mode=k ! wait before looking into mode 
     ! Rq: if all eigenvalues k are spurious...FEAST will not converge
     if (mode<k) then
        do j=1,k-mode
           theta=DZERO
           e=1
           do i=1,fpm(23)
              if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                 if (res(i)>theta) then
                    e=i
                    theta=res(i) !max
                 endif
              end if
           enddo
           trace=trace-lambda(e) ! update trace
           res(e)=-DONE !spurious
        end do
!!! max residual
        theta=DZERO
        do i=1,fpm(23)
           if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
              if (res(i)>theta) theta=res(i)
           end if
        enddo
     end if

     if (mode==0) info=1  ! no eigenvalue detected in the interval
     if (loop>1) then ! wait second iteration (spurious related)
        if ((mode==M0).and.(mode/=N)) info=3 ! size subspace too small
     endif
     if (info/=0) fpm(21)=100 ! The End
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  if (fpm(21)==11) then
     testconv=.false. ! initialization

     !! trace convergence
     if (loop==0) epsout=DONE
     if (loop>0) then
        epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))           
        if ((fpm(6)==0).and.(log10(epsout)<(-fpm(3)))) testconv=.true.
     end if

     !! residual convergence
     if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.

     !! Two loops minimum if spurious are found
     if ((loop<=1).and.(k>mode)) testconv=.false.

     if (fpm(1)==1) then
        call wwrite_i(fout,loop)
        call wwrite_t(fout) 
        call wwrite_i(fout,mode)
        call wwrite_t(fout)
        call wwrite_d(fout,trace)
        call wwrite_t(fout) 
        call wwrite_d(fout,epsout)
        call wwrite_t(fout) 
        call wwrite_d(fout,theta)
        call wwrite_n(fout) 
     end if

     if (.not.testconv) then
        epsout=trace
        if (loop==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

  if (fpm(21)==11) then

     if (.not.testconv) then !!! need FEAST refinement loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! rational function for the inner eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        res(1:fpm(23))=DZERO ! temp indicator for being inside the contour
        call wallocate_1d(work_loc,fpm(23),infoloc)
        call dfeast_rationalx(Zne,Wne,fpm(2),lambda,fpm(23),work_loc)
        if ((fpm(29)==1).and.(fpm(16)==2)) then ! zolotarev case
        call dset_feast_zolotarev(fpm(2),0,zxe,zwe)
        work_loc(1:fpm(23))=work_loc(1:fpm(23))+zwe
        endif
        do i=1,fpm(23)
           if ((lambda(i)>Emin).and.(lambda(i)<Emax))  res(i)=1.0d0
           lambda(i)=work_loc(i)
        enddo
        call wdeallocate_1d(work_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Remark: q and work=S*q is known already, using FEAST-MPI work and q are already distributed


        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPI       
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif

!!! reorder (shift) lambda, eigenvector and residual
        call DLACPY( 'F', N, fpm(23),q , N, work, N ) 
        call wallocate_1d(work_loc,fpm(23),infoloc)
        call wallocate_1d(work_loc2,fpm(23),infoloc)
        call DCOPY(fpm(23),lambda , 1, work_loc, 1 )
        call DCOPY(fpm(23),res , 1, work_loc2, 1 )
        q(1:N,1:fpm(23))=DZERO
        e=0
        j=0
        k=0
        do i=1,fpm(23)
           if ((work_loc(i)>Emin).and.(work_loc(i)<Emax)) then ! inside the search interval
              if (work_loc2(i)/=-DONE) then ! spurious 
                 e=e+1
                 q(1:N,e)=work(1:N,i)
                 lambda(e)=work_loc(i)
                 res(e)=work_loc2(i)
              endif
           else
              j=j+1
              q(1:N,mode+j)=work(1:N,i)
              lambda(mode+j)=work_loc(i)
              res(mode+j)=work_loc2(i)
           end if
           if (work_loc2(i)==-DONE) then ! spurious at the end
              k=k+1
              q(1:N,fpm(23)-k+1)=work(1:N,i)
              lambda(fpm(23)-k+1)=work_loc(i)
              res(fpm(23)-k+1)=work_loc2(i)
           endif
        enddo
        call wdeallocate_1d(work_loc)
        call wdeallocate_1d(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        M0=fpm(23)  ! update value of M0 (new subspace)
        fpm(21)=100 ! The End

     end if ! test convergence
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==100) then !! THE END (ijob=0) 
     ijob=0 !! exit FEAST

     if (fpm(1)==1) then !! Print  Information

        if (info>=200) then
           call wwrite_s(fout, 'PROBLEM with input parameters')
           call wwrite_n(fout) 
        end if

        if ((info>100).and.(info<200)) then
           call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
           call wwrite_n(fout) 
        end if

        if (info==-3) then
           call wwrite_s(fout, 'ERROR with reduced system')  
           call wwrite_n(fout) 
        end if

        if (info==-2) then
           call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
           call wwrite_n(fout) 
        end if

        if (info==-1) then
           call wwrite_s(fout, 'ERROR with Internal memory allocation')  
           call wwrite_n(fout) 
        end if

        if (info==1) then
           call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
           call wwrite_n(fout)
        endif

        if (info==3) then
           call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
           call wwrite_n(fout)
        end if

        if (info==4) then
           call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
           call wwrite_n(fout)
        end if

        if (info==5) then
           call wwrite_s(fout, '==>WARNING: Only stochastic estimation of #eigenvalues returned')  
           call wwrite_n(fout)
        end if

        if (info==2) then
           call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
           call wwrite_n(fout)
        end if

        if (info==0) then
           call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
           call wwrite_n(fout) 
        else
           call wwrite_s(fout, '==>INFO code = ') 
           call wwrite_i(fout,info)
           call wwrite_n(fout)
        end if


        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- END*************************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_n(fout)
     endif
  end if

end subroutine dfeast_srcix






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) - 
  !  Solve generalized Aq=lambda Bq eigenvalue problem within a search interval
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE (or B Identity) 
  !  DOUBLE PRECISION version 
  ! 
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) REAL DOUBLE PRECISION (N,M0)   :  Workspace 
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):  Workspace 
  !  Aq,Sq      (input/output) REAL DOUBLE PRECISION (M0,M0)  :  Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION:  search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                               
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
  !=====================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  double precision, dimension(N,*) ::work
  complex(kind=(kind(1.0d0))),dimension(N,*):: workc
  integer,dimension(*) :: fpm
  double precision,dimension(M0,*):: Aq,Sq
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  double precision,dimension(*)  :: lambda
  double precision,dimension(N,*):: q
  integer :: mode
  double precision,dimension(*) :: res
  integer :: info
!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: dfeast_srcix
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
 
  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour

  call dfeast_srcix(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)

end subroutine dfeast_srci





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) - Includes option for custom integration nodes/weight
  !  Solve generalized Aq=lambda Bq eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE (or Identity)  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                          
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
  !
  !                            Expert comment: if fpm(29)=1 Zne, Wne have already been generated using default contour   
  !=====================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  double precision,dimension(*)  :: lambda
  complex(kind=(kind(1.0d0))),dimension(N,*):: q
  integer :: mode
  double precision,dimension(*) :: res
  integer :: info
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  !! parameters
  double precision, Parameter :: pi=3.1415926535897932d0
  double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZEROC=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0,ab=pi/2.0d0
  integer(8),parameter :: fout =6
  !! variable for FEAST
  integer :: i,e,j,k
  integer,dimension(4) :: iseed
  double precision :: theta,r,Emid
  complex(kind=(kind(1.0d0))) :: zxe,zwe,aux
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: zSqo
  integer, dimension(:),pointer :: fpm_default
  logical :: testconv
  double precision :: trace
  !! Lapack variable (reduced system)
  character(len=1) :: JOBZ,UPLO
  double precision, dimension(:),pointer :: work_loc,work_loc2
  complex(kind=(kind(1.0d0))), dimension(:),pointer :: zwork_loc
  integer :: lwork_loc,info_lap,infoloc
  !! MPI compatibility variables
  integer :: rank,code,nb_procs,NEW_COMM_WORLD

  rank=0
  nb_procs=1
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
  if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
#endif
  !---------------------------------------------

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then 
     info=0 ! default value
     if (fpm(1)==1) then
        call wwrite_n(fout)
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout)
        call wwrite_s(fout, 'Routine ZFEAST_H{}{}')
        if (fpm(29)==0) call wwrite_s(fout,'X') 
        call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
        call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
        call wwrite_n(fout)

        call wallocate_1i(fpm_default,64,infoloc)
        call feastinit(fpm_default)
        do i=1,19
           if (fpm(i)/=fpm_default(i)) then
              call wwrite_s(fout, '   fpm(')
              call wwrite_i(fout, i)
              call wwrite_s(fout, ')=')
              call wwrite_i(fout, fpm(i))
              call wwrite_n(fout)
           endif
        enddo
        call wdeallocate_1i(fpm_default)

        call wwrite_s(fout, 'Search interval [')
        call wwrite_d(fout,Emin)
        call wwrite_s(fout, '; ')
        call wwrite_d(fout,Emax)
        call wwrite_s(fout, ']')
        call wwrite_n(fout)
     end if
     call check_feast_fpm_input(fpm,info)
     call dcheck_feast_srci_input(Emin,Emax,M0,N,info)

     if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!
     IF (info==0) then
        fpm(22)=fpm(2) ! only one half contour necessary here
        loop=0
        if (fpm(1)==1) then
           call wwrite_s(fout, 'Size system    ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,N)
           call wwrite_n(fout)
           call wwrite_s(fout, 'Size subspace  ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,M0)
           call wwrite_n(fout)
           call wwrite_s(fout, '#Linear systems')  
           call wwrite_t(fout) 
           call wwrite_i(fout,fpm(22))
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
           call wwrite_s(fout, '#Loop | #Eig  |       Trace           |     Error-Trace       |     Max-Residual')
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
        end if
        fpm(23)=min(M0,N) ! 'current M0' size (global value)
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors (option 2 in ZLARNV)
           iseed=(/56,890,3456,2333/)
           ! copy work into q for multiplication by B matrix, if stochastic estimate  is "on"
           if (fpm(14)==2) then 
              call ZLARNV(2,iseed,N*fpm(23),q(1,1))
           else
              call ZLARNV(2,iseed,N*fpm(23),work(1,1))
           endif
        endif

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=ZEROC 
#endif
           !------------------------------------------
           ijob=40 !! B*q=>work
           return
        end if
     end IF   ! info=0
!!!!!!!!!!!!!!
  end if !ijob=-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  IF (fpm(21)==1) then !! we initialize a new contour integration

     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(5)==1).or.(fpm(14)==2)) then        
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     endif

     if ((fpm(29)==1).and.(fpm(16)==2)) then ! Zolotarev
        if ((loop>0).and.(fpm(23)/nb_procs>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !------------------------------------------------------------------------

     if ((fpm(29)==1).and.(fpm(16)==2)) then !! Zolotarev
        call dset_feast_zolotarev(fpm(2),0,zxe,zwe)
        !r=(Emax-Emin)/2.0d0
        if ((loop==0).and.(fpm(5)==0)) then
           q(1:N,1:fpm(23))=ZEROC
        else
           if (rank==0) then
           q(1:N,1:fpm(23))=-dble(zwe)*q(1:N,1:fpm(23))
           else
           q(1:N,1:fpm(23))=ZEROC
           end if
        endif
     else
        q(1:N,1:fpm(23))=ZEROC
     end if

     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     IF (info==0) then !! will end up checking info errors returned by FEAST drivers
        do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the contour 

           if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
              Ze=Zne(e)
              fpm(20)=e-rank
              ijob=10 ! for fact
              if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor 
           endif

           if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v 
              call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
              ijob=11 ! for solve
              return
           endif

           if (ijob==11) then 
              !! summation 
              aux=-Wne(e) 
              call ZAXPY(N*fpm(23),aux,workc,1,q,1)
              !!Explicit Factorization of the linear system (complex) (zS-A)^T 
              !!needed if driver not capable to exploit Factorization of (zS-A) for solving (zS-A)^Tq=v          
              ijob=20 ! for fact
              if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
           endif

           if (ijob==20) then!!!! Solve the linear system (complex) (zS-A)^Tq=v  
              call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
              ijob=21 ! for solve with transpose
              return
           end if

           if (ijob==21) then
              aux=-conjg(Wne(e)) 
              call ZAXPY(N*fpm(23),aux,workc,1,q,1)

              ijob=-2 ! just for identification
           end if

        end do
     end IF !! info=0

     !------------------------------------------------
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------------  
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(21)=4 
        !------------------------------------------------
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if
  end IF   ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic Estimate
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only q is returned with stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     call ZLARNV(2,iseed,N*fpm(23),work(1,1)) 
     call ZGEMM('C','N',fpm(23),fpm(23),N,-ONEC,work,N,q,N,ZEROC,zAq,M0) ! projection
     call ZGEMM('C','N',fpm(23),fpm(23),N,ONEC,work,N,work,N,ZEROC,zSq,M0) ! normalization
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+dble(zAq(i,i)/zSq(i,i))
        res(i)=abs(theta*(N/(DONE*i)))
     enddo
     mode=int(real(res(fpm(23))))+1
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Sq xq     with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form Bq=> Bq=Q^T B Q 
  if (fpm(21)==4) then
#ifdef MPI 
     work(1:N,1:fpm(23))=ZEROC
#endif
     fpm(21)=5 ! preparing reentry
     ijob=40 
     return  ! mat-vec B*q => work
  endif

  if (fpm(21)==5) then 
     !------------------------------------------------
#ifdef MPI 
     zSq(1:M0,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     call ZGEMM('C','N',fpm(23),fpm(25),N,ONEC,q(1,1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     mode=0
     if (loop>0) then
        do i=fpm(24),fpm(24)+fpm(25)-1
           if (res(i)==DONE) then ! indicator for being inside the search interval
              if (abs(sqrt(abs(zSq(i,i)))-lambda(i))/sqrt(abs(zSq(i,i)))<0.1d0) mode=mode+1
           endif
        enddo
     endif
!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,NEW_COMM_WORLD,code)
#endif
!!!!!!!!!!!!!!!!!!

     fpm(21)=6
  endif

!!!!!!!!!form  Aq=> Aq=Q^T A Q
  if (fpm(21)==6) then
     !------------------------------------------------
#ifdef MPI 
     work(1:N,1:fpm(23))=ZEROC
#endif 
     fpm(21)=7 ! preparing reenty
     ijob=30 
     return! mat-vec A*q => work
  end if


  if (fpm(21)==7) then

     !------------------------------------------------
#ifdef MPI
     zAq(1:M0,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     call ZGEMM('C','N',fpm(23),fpm(25),N,ONEC,q(1,1),N,work(1,fpm(24)),N,ZEROC,zAq(1,fpm(24)),M0)
     fpm(21)=8
  endif


  if (fpm(21)==8) then
     !----------------------------------------!(zAq,zSq known to all processors) 
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zSq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !---------------------------------------
     if (fpm(12)==1) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=50
        return
     endif
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==8) then
     ! if using FEAST-MPI ==> solve on a single proc
     if (rank==0) then
        JOBZ='V'
        UPLO='L'
        info_lap=1 ! initialization
        i=1
        LWORK_LOC=2*fpm(23)-1 !! for lapack eig reduced system
        call wallocate_1z(zWORK_LOC,LWORK_LOC,infoloc)
        call wallocate_1d(WORK_LOC,3*fpm(23)-2,infoloc)
        if (infoloc/=0) info=-1

        do while ((info_lap/=0).and.(info==0))
           i=i+1
           if (i==10) info=-3 ! arbitrary maximum
           call wallocate_2z(zSqo,fpm(23),fpm(23),infoloc)
           if (infoloc/=0) info=-1

           call ZLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )
           call ZHEGV(1, JOBZ, UPLO, fpm(23), zAq, M0, zSqo, fpm(23), lambda, zWORK_loc,Lwork_loc, WORK_loc, INFO_lap)

           if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
           if (info_lap>fpm(23)) then !! zSqo is not spd (a posteriori resize subspace)
              fpm(23)=info_lap-fpm(23)-1
              if (fpm(1)==1) then
                 call wwrite_s(fout, 'Resize subspace')  
                 call wwrite_t(fout) 
                 call wwrite_i(fout,fpm(23))
                 call wwrite_n(fout)
              end if
           end if
           call wdeallocate_2z(zSqo)
        end do
        call wdeallocate_1z(zwork_loc)
        call wdeallocate_1d(work_loc)
     end if !(rank 0)
     !-------------------------------- !(info common to all processors) 
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(25)=fpm(23) !! current M0 size (by default) -- global
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 
#ifdef MPI
        call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(zAq,M0*fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,NEW_COMM_WORLD,code)
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of current M0
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------
        fpm(21)=9
     end if
  end if !! fpm(21)=8



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors X=Qxq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     call ZLACPY( 'F', N, fpm(23),q , N, work, N ) 
     !! option - non shifted
#ifdef MPI 
     q(1:N,1:fpm(23))=ZEROC
#endif
     call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,work(1,1),N,zAq(1,fpm(24)),M0,ZEROC,q(1,fpm(24)),N)

  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     fpm(21)=10 ! preparing reentry
     ijob=30 
     return  ! mat-vec A*q => work
  endif

  if (fpm(21)==10) then
     call ZLACPY( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N )
     fpm(21)=11 ! preparing reentry
     ijob=40 
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZEROC  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     return  ! mat-vec S*q => work 
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !-------- Absolute residual
     do i=fpm(24),fpm(24)+fpm(25)-1
        res(i)=sum(abs(workc(1:N,i)-lambda(i)*work(1:N,i)))/sum(abs(max(abs(Emin),abs(Emax))*(work(1:N,i))))
     end do
     !----------------------------------------
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==11) then
     k=0
     trace=DZERO
     theta=DZERO
!!! count how many eigenvalues have converged + trace + max residual
     do i=1,fpm(23)
        if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
           k=k+1
           trace=trace+lambda(i)
           if (res(i)>theta) theta=res(i) ! max residual
        endif
     enddo

!!!!!!!!! remove spurious if any       
     if ((mode==0).and.(k>0)) mode=k ! wait before looking into mode 
     ! Rq: if all eigenvalues k are spurious...FEAST will not converge
     if (mode<k) then
        do j=1,k-mode
           theta=DZERO
           e=1
           do i=1,fpm(23)
              if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
                 if (res(i)>theta) then
                    e=i
                    theta=res(i) !max
                 endif
              end if
           enddo
           trace=trace-lambda(e) ! update trace
           res(e)=-DONE !spurious
        end do
!!! max residual
        theta=DZERO
        do i=1,fpm(23)
           if ((lambda(i)>Emin).and.(lambda(i)<Emax)) then ! inside the search interval
              if (res(i)>theta) theta=res(i) 
           end if
        enddo
     end if

     if (mode==0) info=1  ! no eigenvalue detected in the interval
     if (loop>1) then! wait second iteration (spurious related)
        if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
     endif
     if (info/=0) fpm(21)=100 ! The End
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  if (fpm(21)==11) then
     testconv=.false. ! initialization

     !! trace convergence
     if (loop==0) epsout=DONE
     if (loop>0) then
        epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))           
        if ((fpm(6)==0).and.(log10(epsout)<(-fpm(3)))) testconv=.true.
     end if

     !! residual convergence
     if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.

     !! Two loops minimum if spurious are found
     if ((loop<=1).and.(k>mode)) testconv=.false.

     if (fpm(1)==1) then
        call wwrite_i(fout,loop)
        call wwrite_t(fout) 
        call wwrite_i(fout,mode)
        call wwrite_t(fout)
        call wwrite_d(fout,trace)
        call wwrite_t(fout) 
        call wwrite_d(fout,epsout)
        call wwrite_t(fout) 
        call wwrite_d(fout,theta)
        call wwrite_n(fout) 
     end if

     if (.not.testconv) then
        epsout=trace
        if (loop==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

  if (fpm(21)==11) then

     if (.not.testconv) then !!! need FEAST refinement loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! rational function for the inner eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        res(1:fpm(23))=DZERO ! temp indicator for being inside the contour
        call wallocate_1d(work_loc,fpm(23),infoloc)
        call dfeast_rationalx(Zne,Wne,fpm(2),lambda,fpm(23),work_loc)
         if ((fpm(29)==1).and.(fpm(16)==2)) then ! zolotarev case
        call dset_feast_zolotarev(fpm(2),0,zxe,zwe)
        work_loc(1:fpm(23))=work_loc(1:fpm(23))+zwe
        endif
        do i=1,fpm(23)
           if ((lambda(i)>Emin).and.(lambda(i)<Emax))  res(i)=DONE
           lambda(i)=work_loc(i)
        enddo
        call wdeallocate_1d(work_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Remark: q and work=S*q is known already,  using FEAST-MPI work and q are already distributed


        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPI       
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
#endif

!!! reorder (shift) lambda, eigenvector and residual

        call ZLACPY( 'F', N, fpm(23),q , N, work, N )
        call wallocate_1d(work_loc,fpm(23),infoloc)
        call wallocate_1d(work_loc2,fpm(23),infoloc)
        call DCOPY(fpm(23),lambda , 1, work_loc, 1 )
        call DCOPY(fpm(23),res , 1, work_loc2, 1 )
        q(1:N,1:fpm(23))=ZEROC
        e=0
        j=0
        k=0
        do i=1,fpm(23)
           if ((work_loc(i)>Emin).and.(work_loc(i)<Emax)) then ! inside the search interval
              if (work_loc2(i)/=-DONE) then ! spurious 
                 e=e+1
                 q(1:N,e)=work(1:N,i)
                 lambda(e)=work_loc(i)
                 res(e)=work_loc2(i)
              endif
           else
              j=j+1
              q(1:N,mode+j)=work(1:N,i)
              lambda(mode+j)=work_loc(i)
              res(mode+j)=work_loc2(i)
           end if
           if (work_loc2(i)==-DONE) then ! spurious at the end
              k=k+1
              q(1:N,fpm(23)-k+1)=work(1:N,i)
              lambda(fpm(23)-k+1)=work_loc(i)
              res(fpm(23)-k+1)=work_loc2(i)
           endif
        enddo
        call wdeallocate_1d(work_loc)
        call wdeallocate_1d(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        M0=fpm(23)  ! update value of M0 (new subspace)
        fpm(21)=100 ! The End

     end if ! test convergence
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==100) then !! THE END (ijob=0) 
     ijob=0 !! exit FEAST

     if (fpm(1)==1) then !! Print  Information

        if (info>=200) then
           call wwrite_s(fout, 'PROBLEM with input parameters')
           call wwrite_n(fout) 
        end if

        if ((info>100).and.(info<200)) then
           call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
           call wwrite_n(fout) 
        end if

        if (info==-3) then
           call wwrite_s(fout, 'ERROR with reduced system')  
           call wwrite_n(fout) 
        end if

        if (info==-2) then
           call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
           call wwrite_n(fout) 
        end if

        if (info==-1) then
           call wwrite_s(fout, 'ERROR with Internal memory allocation')  
           call wwrite_n(fout) 
        end if

        if (info==1) then
           call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
           call wwrite_n(fout)
        endif

        if (info==3) then
           call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
           call wwrite_n(fout)
        end if

        if (info==4) then
           call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
           call wwrite_n(fout)
        end if

        if (info==5) then
           call wwrite_s(fout, '==>WARNING: Only stochastic estimation of #eigenvalues returned')  
           call wwrite_n(fout)
        end if

        if (info==2) then
           call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
           call wwrite_n(fout)
        end if

        if (info==0) then
           call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
           call wwrite_n(fout) 
        else
           call wwrite_s(fout, '==>INFO code = ') 
           call wwrite_i(fout,info)
           call wwrite_n(fout)
        end if


        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- END*************************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_n(fout)
     endif
  end if

end subroutine zfeast_hrcix






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_hrci(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) 
  !  Solve generalized Aq=lambda Bq eigenvalue problems within a given search interval
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE (or B Identity)  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,40)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output) INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION:  search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                          
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
  !=====================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  double precision,dimension(*)  :: lambda
  complex(kind=(kind(1.0d0))),dimension(N,*):: q
  integer :: mode
  double precision,dimension(*) :: res
  integer :: info
!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_hrcix
!!!    Create Zne, Wne arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))),dimension(fpm(2)) :: Zne,Wne 

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne) 

  fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour

  call zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)

end subroutine zfeast_hrci






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      Aqr=lambda Bqr and A^Hql=conjg(lambda) B^Hql                            
  ! 
  !  Expert comments: This is also a Kernel routine (contains several cases depending of fpm(31))
  !                   fpm(31)=0 ! complex general (default)  
  !                   fpm(31)=1 ! Real general (contour symmetric via real axis)           
  !                   fpm(31)=2 ! Complex symmetric          
  !                   fpm(31)=3 ! Real general               
  !
  !  A and B can be NON-HERMITIAN (could be real or complex)  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,2*M0):   Workspace
  !                            expert comment: size of (N,M0) if fpm(31)=2 (complex symmetric case)  
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input) Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !                            expert comments: (input)  Emid,r provided by user if default contour i.e. fpm(29)=1
  !                                             (output) Emid,r unchanged if fpm(29)=1 
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !                            expert comment: size of (N,M0) if fpm(31)=2 (complex symmetric case - only right vectors)
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right (1:mode) and left (M0+1:M0+mode)
  !                          
  !                            expert comment: size of (M0) if fpm(31)=2 (complex symmetric case - only rigth vectors)
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !
  !                            expert comments: if fpm(29)=1, Wne,Zne are generated using default contour
  !                                             important for testing if an eigenvalue is inside an ellipsoid shape vs polygon contour
  !=========================================================================
  !  James Kestyn - Eric Polizzi - 2015
  !=====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
  integer,dimension(*) :: fpm
  double precision :: epsout
  integer :: loop
  complex(kind=(kind(1.0d0))),dimension(*):: lambda,Zne,Wne
  complex(kind=(kind(1.0d0))),dimension(N,*):: Q
  integer :: mode
  double precision,dimension(*) :: res
  integer :: info
  !! parameters
  double precision, Parameter :: pi=3.1415926535897932d0
  double precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZEROC=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0,ab=pi/2.0d0
  integer*8,parameter :: fout =6
  !! variable for FEAST
  integer :: i,e,j,k,jj,m_min
  integer,dimension(4) :: iseed
  double precision :: theta
  complex(kind=(kind(1.0d0))) :: jac,aux
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: zSqo
  integer, dimension(:),pointer :: fpm_default,color
  logical :: testconv
  complex(kind=(kind(1.0d0))) :: trace, zalpha
  !! Lapack variable 
  double precision, dimension(:),pointer :: work_locr,work_loc2
  complex(kind=(kind(1.0d0))), dimension(:),pointer :: zwork_loc
  complex(kind=(kind(1.0d0))), dimension(:),pointer :: LALPHA,LBETA
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer :: VL,VR
  integer :: lwork_loc,info_lap,infoloc
  character(len=1) :: JOBVL
  !! MPI compatibility variables
  integer :: rank,code,nb_procs,NEW_COMM_WORLD


  rank=0
  nb_procs=1
  !----------------------------------------------
#ifdef MPI
  NEW_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)
  if (rank/=0) fpm(1)=0 ! comment only in rank 0 if any
#endif
  !---------------------------------------------

  if ((fpm(29)==1).and.(fpm(31)==3)) then ! default contour
     !Test if Emid is real and fpm(8) even -> factorization on half-contour -- New ID if contour symmetric via real axis
     if ((aimag(Emid)==0.0d0).and.((2*(fpm(8)/2))==fpm(8))) fpm(31)=1 
  end if
  
  if (fpm(14)==2) fpm(31)=2 ! Compute only right contour if estimate is on


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then 
     info=0 ! default value
     if(fpm(29)==0)then ! expert routines Emid,r needs to be calculated
        ! calculate the center of mass and relative maximum length of custom contour
        Emid=ZEROC
        do i=1,fpm(8)
           Emid = Emid + Zne(i)
        enddo
        Emid = Emid / fpm(8)
        r = abs(maxval(abs(Zne(1:fpm(8))))-Emid)
     endif
     if (fpm(1)==1) then
        call wwrite_n(fout)
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout) 
        call wwrite_s(fout, '*********** FEAST- BEGIN **********************')
        call wwrite_n(fout) 
        call wwrite_s(fout, '***********************************************')  
        call wwrite_n(fout)
        if (fpm(31)==3) then
           call wwrite_s(fout, 'Routine DFEAST_G{}{}') !<<<<<<<<
        elseif (fpm(31)==1) then
           call wwrite_s(fout, 'Routine DFEAST_G{}{}') !<<<<<<<<
        elseif (fpm(31)==2) then
           call wwrite_s(fout, 'Routine ZFEAST_S{}{}') !<<<<<<<<
        elseif (fpm(31)==0) then
           call wwrite_s(fout, 'Routine ZFEAST_G{}{}') !<<<<<<<<
        end if
        if (fpm(29)==0) call wwrite_s(fout,'X') 
        call wwrite_n(fout)
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
        call wwrite_s(fout, 'List of input parameters fpm(1:64)-- if different from default')
        call wwrite_n(fout)

        call wallocate_1i(fpm_default,64,infoloc)
        call feastinit(fpm_default)
        do i=1,19
           if (fpm(i)/=fpm_default(i)) then
              call wwrite_s(fout, '   fpm(')
              call wwrite_i(fout, i)
              call wwrite_s(fout, ')=')
              call wwrite_i(fout, fpm(i))
              call wwrite_n(fout)
           endif
        enddo
        call wdeallocate_1i(fpm_default)


        if(fpm(1)==1 .and. rank==0)then
           if (fpm(29)==1) then
           call wwrite_s(fout, 'Search contour [(Emid), r] (')
           else
           call wwrite_s(fout, 'Center of mass for user custom contour (')
           endif
           call wwrite_d(fout,dble(Emid))
           call wwrite_s(fout, ',')
           call wwrite_d(fout,aimag(Emid))
           call wwrite_s(fout, '); ')
           if (fpm(29)==1) call wwrite_d(fout,r)
           call wwrite_n(fout)
        endif
     end if
     call check_feast_fpm_input(fpm,info)
     call dcheck_feast_grci_input(r,M0,N,info)

     if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!
     IF (info==0) then
        fpm(22)=fpm(8) ! default
        if (fpm(29)==1) then! default generated contour 
          if (fpm(31)==1) fpm(22)=fpm(8)/2  ! fpm(31)=1 -- Emid real--fpm(8) even-   will only have half-contour
        endif  
        loop=0
        if (fpm(1)==1) then
           call wwrite_s(fout, 'Size system    ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,N)
           call wwrite_n(fout)
           call wwrite_s(fout, 'Size subspace  ')  
           call wwrite_t(fout) 
           call wwrite_i(fout,M0)
           call wwrite_n(fout)
           call wwrite_s(fout, '#Linear systems')  
           call wwrite_t(fout) 
           call wwrite_i(fout,fpm(22))
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
           call wwrite_s(fout, '#Loop | #Eig  |          |Trace|       |     Error-Trace       |     Max-Residual')
           call wwrite_n(fout)
           call wwrite_s(fout, '-----------------------------------------------------------------------------------')
           call wwrite_n(fout)
        endif
        fpm(23)=min(M0,N) ! 'current M0' size (global value)
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        fpm(35)=fpm(23) !! 'current M0' size (by default)
        fpm(34)=M0+1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
           fpm(35)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(35)=fpm(35)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(34)=M0+1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors
           iseed=(/56,890,3456,2333/)
           if (fpm(14)==2) then
              ! copy work into q for multiplication by B matrix, if stochastic estimate  is "on" 
              call ZLARNV(4,iseed,N*fpm(23),q(1,1))
           else
              call ZLARNV(4,iseed,N*fpm(23),work(1,1))
              iseed=(/6,80,456,5421/)
              if (fpm(31)/=2) call ZLARNV(4,iseed,N*fpm(23),work(1,M0+1))
           endif
        end if

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=ZEROC 
           if (fpm(31) /= 2) work(1:N,M0+1:M0+fpm(23))=ZEROC 
#endif
           !------------------------------------------
           ijob=40 !! B*Q(1,1:fpm(23))=>work(1:N,1:fpm(23))
           if (fpm(31) /=2 ) fpm(21)=-1
           return
        end if
     end IF ! info=0
!!!!!!!!!!!!!!
  end if   !ijob=-1

  if (fpm(21)==-1) then !! we need to initialize the left eigenvectors as well
     ijob=41!! B^T*Q(1,M0+1:M0+fpm(23))=>work(1:N,M0+1:M0+fpm(23))
     fpm(21)=1 ! prepare reentry
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==1) then !! we initialize a new contour integration
     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(5)==1).or.(fpm(14)==2)) then        
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        if ((fpm(31)/=2) .and. (fpm(23)/nb_procs>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !------------------------------------------------------------------------
     Q(1:N,1:fpm(23))=ZEROC
     if(fpm(31)/=2) Q(1:N,M0+1:M0+fpm(23))=ZEROC
     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     IF (info==0) then !! will end up checking info errors returned by FEAST drivers
        do e=fpm(20)+rank,fpm(22),nb_procs !!!! loop over the full contour  

           if (ijob==-2) then !!Factorize the linear system                   
              Ze = Zne(e)   
              fpm(20)=e-rank
              ijob=10 ! for fact
              if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor 
           end if

           if (ijob==10) then !!Solve the linear system (complex) (zS-A)qr=workc 
              call ZLACPY( 'F', N, fpm(23), work, N, workc, N )
              ijob=11 ! for solve
              return
           end if

           if (ijob==11) then  !! summation
              aux = Wne(e) 
              if (fpm(33)==1) then
                  aux=conjg(aux)
                  workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))
               end if
              call ZAXPY(N*fpm(23),aux,workc,1,Q(1,1),1)

              if(fpm(31)==1) then
                 if (fpm(33)==0) then
                    call ZLACPY( 'F', N, fpm(23),work(1,1) , N, workc, N )
                    workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))
                    fpm(33)=1 
                    ijob=11 ! reuse fact for second solve - must conjugate the solution qr
                    return
                 elseif (fpm(33)==1) then
                    fpm(33)=0
                    ijob=20
                    if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
                 end if
              elseif (fpm(31)==2) then
                 ijob=-2 !no solve for ql if complex symmetric
              elseif ((fpm(31)==0).or.(fpm(31)==3)) then 
                 ijob=20 ! for fact
                 if ((loop==0).or.(fpm(22)>nb_procs)) return ! no need to factorize again if one linear system per processor
              end if
           end if !ijob==11

           if (ijob==20) then !!Solve the linear system (complex) (zS-A)^Hql=workc 
              call ZLACPY( 'F', N, fpm(23),work(1,M0+1) , N, workc, N )
              ijob=21 ! for solve
              return 
           end if

           if (ijob==21) then  
              aux = conjg(Wne(e)) 
              if (fpm(33)==1) then
                  aux=conjg(aux)
                  workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))                   
              end if
              call ZAXPY(N*fpm(23),aux,workc,1,Q(1,M0+1),1)

              ijob=-2 ! default identification

              if (fpm(31)==1) then
                 if (fpm(33)==0) then
                    call ZLACPY( 'F', N, fpm(23),work(1,M0+1) , N, workc, N )
                    workc(1:N,1:fpm(23))= conjg(workc(1:N,1:fpm(23)))
                    fpm(33)=1
                    ijob=21 ! reuse fact for second solve - must conjugate the solution qr
                    return
                 elseif (fpm(33)==1) then
                    fpm(33)=0
                 endif
              end if

           end if
        end do
     end IF !! info=0

     !------------------------------------------------
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------------  
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(21)=4 
        !------------------------------------------------
#ifdef MPI
        if(fpm(31)/=2) call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if


  end IF   ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only qr,ql vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic Estimate
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only qr is returned with stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     call ZLARNV(4,iseed,N*fpm(23),work(1,1)) 
     call ZGEMM('C','N',fpm(23),fpm(23),N,-ONEC,work,N,q,N,ZEROC,zAq,M0) ! projection
     call ZGEMM('C','N',fpm(23),fpm(23),N,ONEC,work,N,work,N,ZEROC,zSq,M0) ! normalization
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+dble(zAq(i,i)/zSq(i,i))
        res(i)=abs(theta*(N/(DONE*i)))
     enddo
     mode=int(real(res(fpm(23))))+1
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Sq xq     with Aq=Q^TAQ Sq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!form  Bq=> Bq=Q^T B Q
  if (fpm(21)==4) then 
     !------------------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     fpm(21)=5 ! preparing reenty
     ijob=40 
     return! mat-vec B.qr => workr
  end  if

  if (fpm(21)==5) then

     ! after first loop zSq(1,1) contains old residual-
     if (loop>0) Ze=zSq(1,1) !  dummy variable

     !------------------------------------------------
#ifdef MPI
     zSq(1:M0,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     if(fpm(31)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,ONEC,Q(1,1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     else
        call ZGEMM('C','N',fpm(23),fpm(35),N,ONEC,Q(1,M0+1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     endif
     !---------------------------------------
#ifdef MPI
     if (fpm(23)/nb_procs>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zSq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !       !---------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     mode=0
     if ((fpm(13)==1) .and. (fpm(37)==0)) then
        if (loop>0) then
           do i=fpm(24),fpm(24)+fpm(25)-1
              if (res(i)==DONE) then ! indicator for being inside the search interval
                 if( abs( (sqrt(zSq(i,i)))-lambda(i) ) / sqrt(abs(zSq(i,i))) < 0.1d0) mode=mode+1
              end if
           end do
        endif
     endif
!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,NEW_COMM_WORLD,code)
#endif
!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! B-biorthogonalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     !-------------------------
     ! Diagonalize Bq=Ql^H B Qr
     !-------------------------

     e=0
     call wallocate_2z(VR,fpm(23),fpm(23),infoloc)
     call wallocate_2z(VL,fpm(23),fpm(23),infoloc)
     call wallocate_1z(LALPHA,fpm(23),infoloc)

     if (rank==0) then
        call wallocate_2z(zSqo,fpm(23),fpm(23),infoloc)
        call ZLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )

        LWORK_LOC=8*fpm(23) !! for lapack eig reduced system
        call wallocate_1z(zWORK_LOC,LWORK_LOC,infoloc)
        call wallocate_1d(WORK_LOCR,2*fpm(23),infoloc)
        if (infoloc/=0) info=-1

        JOBVL='V'
        if (fpm(31)==2) JOBVL='N'        
        call ZGEEV(JOBVL,'V',fpm(23),zSqo,fpm(23),LALPHA,VL,fpm(23),VR,fpm(23),zWORK_LOC,LWORK_LOC,WORK_LOCR,info_lap)

        theta = abs(lalpha(1)) !! Find max lalpha
        do j=1,fpm(23)
           if(theta<abs(lalpha(j))) theta=abs(lalpha(j))
        enddo

        e=0
        do j=1,fpm(23)  !! find dimension of NULL(Bq)=e
           if( abs(LALPHA(j))/theta < 1.0d-16 ) e = e+1
        enddo

        do i=1,e  !! Order eigenpairs by |lalpha_i|
           m_min=1
           theta = abs(lalpha(1))
           do j=1,fpm(23) - i + 1 !find min lalpha and put at end
              if( abs(LALPHA(j)) < theta )then
                 theta = abs(lalpha(j))
                 m_min = j
              endif
           enddo
           call ZSWAP(fpm(23),VR(1,m_min),1,VR(1,fpm(23)-i+1),1) 
    if (fpm(31)/=2) call ZSWAP(fpm(23),VL(1,m_min),1,VL(1,fpm(23)-i+1),1) 
           jac = LALPHA(fpm(23)-i+1)
           LALPHA(fpm(23)-i+1) = LALPHA(m_min)
           LALPHA(m_min)=jac
        enddo

        call wdeallocate_1z(zWORK_LOC)
        call wdeallocate_1d(WORK_LOCR)
        call wdeallocate_2z(zSqo)
     endif ! rank==0

     !-------------
#ifdef MPI
     call MPI_BCAST(e,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
     call MPI_BCAST(VR,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
if (fpm(31)/=2)  call MPI_BCAST(VL,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
     call MPI_BCAST(lalpha,fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
#endif

     !--------------
     ! Create Ur and Ul (bi-B-orthogonal)
     !--------------------
     IF (fpm(15)==1) then  

        call ZLACPY('F',N,fpm(23),Q(1,1),N,workc(1,1),N)
        Q(1:N,1:fpm(23))=ZEROC
        call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,workc(1,1),N,VR(1,fpm(24)),fpm(23),ZEROC,Q(1,fpm(24)),N)
        call ZLACPY('F',N,fpm(23),work(1,1),N,workc(1,1),N)
        work(1:N,1:fpm(23))=ZEROC
        call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,workc(1,1),N,VR(1,fpm(24)),fpm(23),ZEROC,work(1,fpm(24)),N)
        if(fpm(31)/=2) then
           call ZLACPY('F',N,fpm(23),Q(1,M0+1),N,workc(1,1),N)
           Q(1:N,M0+1:M0+fpm(23))=ZEROC
           call ZGEMM('N','N',N,fpm(35),fpm(23),ONEC,workc(1,1),N,VL(1,fpm(24)),fpm(23),ZEROC,Q(1,fpm(34)),N)
        endif

        !--------------------
#ifdef MPI
        if (fpm(23)/nb_procs>=1) then
           call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
           call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
           if (fpm(31)/=2) call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        endif
#endif

        !----------------------------------------
        !---------- Normalization of Ur and Ul
        !----------------------------------------          
        call ZGEMM('N','N',fpm(23),fpm(23),fpm(23),ONEC,zSq(1,1),M0,VR(1,1),fpm(23),ZEROC,zAq(1,1),M0)
if (fpm(31)/=2) then
        call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ONEC,VL(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSq(1,1),M0)
else
        call ZGEMM('T','N',fpm(23),fpm(23),fpm(23),ONEC,VR(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSq(1,1),M0)
endif

        do i=1,fpm(23)
           call ZSCAL(N,ONEC/sqrt(zSq(i,i)),Q(1,i),1)
           call ZSCAL(N,ONEC/sqrt(zSq(i,i)),work(1,i),1)
           if(fpm(31)/=2) call ZSCAL(N,ONEC/sqrt(conjg(zSq(i,i))),Q(1,M0+i),1)
        enddo

     end IF ! fpm(15)==1

     call wdeallocate_1z(LALPHA)
     call wdeallocate_2z(VR)
     call wdeallocate_2z(VL)

     !        if ((e/=0) .and. (loop>0))then
     if (e/=0) then
        fpm(23) = fpm(23)-e
        fpm(25)=fpm(23) ! default - single proc
        fpm(35)=fpm(23) ! default - single proc
        if (fpm(1)==1) then 
           call wwrite_s(fout,'Resize Subspace (Bq)')
           call wwrite_t(fout)
           call wwrite_i(fout,fpm(23))
           call wwrite_n(fout)
        endif


        if (fpm(23)/nb_procs>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(24)=1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 

           fpm(35)=fpm(23)/nb_procs ! local size of 'current M0'
           if (rank==nb_procs-1) fpm(35)=fpm(35)+fpm(23)-(fpm(23)/nb_procs)*nb_procs 
           fpm(34)=M0+1+rank*(fpm(23)/nb_procs) ! local origin of first column for vector q for parallel mat-vec 
        end if

     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !--------------------------------
#ifdef MPI
     zSq(1:M0,1:fpm(23))=ZEROC
#endif
     !--------------------------------
     if(fpm(31)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,ONEC,Q(1,1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     else
        call ZGEMM('C','N',fpm(23),fpm(35),N,ONEC,Q(1,M0+1),N,work(1,fpm(24)),N,ZEROC,zSq(1,fpm(24)),M0)
     endif


     fpm(21)=6
  endif


!!!!!!!!!form Aq=> Aq=Q^T A Q 
  if (fpm(21)==6) then
     !------------------------------------------------
#ifdef MPI 
     work(1:N,1:fpm(23))=ZEROC
#endif 
     fpm(21)=7 ! preparing reentry
     ijob=30 
     return  ! mat-vec A.qr => workr
  endif


  if (fpm(21)==7) then 
     !------------------------------------------------
#ifdef MPI 
     zAq(1:M0,1:fpm(23))=ZEROC
#endif
     !-------------------------------------------------
     if(fpm(31)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,ONEC,Q(1,1),N,work(1,fpm(24)),N,ZEROC,zAq(1,fpm(24)),M0)
     else
        call ZGEMM('C','N',fpm(23),fpm(35),N,ONEC,Q(1,M0+1),N,work(1,fpm(24)),N,ZEROC,zAq(1,fpm(24)),M0)
     endif

     fpm(21)=8
  endif



  if (fpm(21)==8) then
     !----------------------------------------!(zAq,zSq known to all processors) 
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zSq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
     end if
#endif
     !---------------------------------------
     if (fpm(12)==1) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=50
        return
     endif
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==8) then
     ! if using FEAST-MPI ==> solve on a single proc
     if (rank==0) then
        info_lap=1 ! initialization
        LWORK_LOC=8*fpm(23) !! for lapack eig reduced system
        call wallocate_1z(zWORK_LOC,LWORK_LOC,infoloc)
        call wallocate_1d(WORK_LOCR,8*fpm(23),infoloc)
        if (infoloc/=0) info=-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call wallocate_2z(zSqo,fpm(23),fpm(23),infoloc)
        call wallocate_2z(VL,fpm(23),fpm(23),infoloc)
        call wallocate_2z(VR,fpm(23),fpm(23),infoloc)

        if (infoloc/=0) info=-1

        call ZLACPY( 'F', fpm(23), fpm(23),zSq , M0, zSqo, fpm(23) )

        if (fpm(39)==1) then  !!! standard eigenvalue

           do i=1,fpm(23)
              zAq(1:fpm(23),i) = zAq(1:fpm(23),i) / sqrt(zSq(i,i))
              zAq(i,1:fpm(23)) = zAq(i,1:fpm(23)) / sqrt(zSq(i,i))
           enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           JOBVL='V'
           if (fpm(31)==2) JOBVL='N'  
           call ZGEEV(JOBVL,'V',fpm(23),zAq,M0,LAMBDA,VL,fpm(23),VR,fpm(23), zWORK_LOC, LWORK_LOC, WORK_LOCR,INFO_lap)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(info_lap==0)then
              do j=1,fpm(23)
                 call ZSCAL(fpm(23),ONEC/sqrt(zSq(j,j)),VR(j,1),fpm(23))
       if (fpm(32)/=2)  call ZSCAL(fpm(23),ONEC/conjg(sqrt(zSq(j,j))),VL(j,1),fpm(23))
              end do
           endif

        else                !!! generalized eigenvalue

           call wallocate_1z(LALPHA,fpm(23),infoloc)
           call wallocate_1z(LBETA,fpm(23),infoloc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           JOBVL='V'
           if (fpm(31)==2) JOBVL='N'  
           call ZGGEV(JOBVL,'V',fpm(23),zAq,M0,zSqo,fpm(23),LALPHA,LBETA,VL,fpm(23),VR,fpm(23), zWORK_LOC, LWORK_LOC, WORK_LOCR,INFO_lap)
           LAMBDA(1:fpm(23))=LALPHA(1:fpm(23))/LBETA(1:fpm(23))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call wdeallocate_1z(LALPHA)
           call wdeallocate_1z(LBETA)
        endif

        if (info_lap/=0) info=-3 ! problem with lapack

        IF (info_lap==0) then

           call ZGEMM('N','N',fpm(23),fpm(23),fpm(23),ONEC,zSq(1,1),M0,VR(1,1),fpm(23),ZEROC,zAq(1,1),M0)
           if(fpm(31)==2) then
              call ZGEMM('T','N',fpm(23),fpm(23),fpm(23),ONEC,VR(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSqo(1,1),fpm(23))
           else
              call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ONEC,VL(1,1),fpm(23),zAq(1,1),M0,ZEROC,zSqo(1,1),fpm(23))
           endif


           fpm(37)=0
           do j=1,fpm(23) !!!!!! Check othogonality of ZGGEV eigenvectors
              zalpha = zSqo(j,j)
              do jj=1,fpm(23)
                 zalpha = zalpha - zSqo(jj,j)
              enddo
              if (log10(abs(zalpha)/(fpm(23)))>-10) then
                 fpm(37) = 1 
              endif
           enddo

           if (fpm(37)==0) then  !!!!!! normalize
              do j=1,fpm(23)
                 call ZSCAL(fpm(23),ONEC/sqrt(zSqo(j,j)),VR(1,j),1)
                 if(fpm(31)/=2) call ZSCAL(fpm(23),ONEC/conjg(sqrt(zSqo(j,j))),VL(1,j),1)
              enddo
           endif

           call ZLACPY( 'F', fpm(23), fpm(23),VR , fpm(23), zAq, M0 )
    if (fpm(32)/=2) call ZLACPY( 'F', fpm(23), fpm(23),VL , fpm(23), zSq, M0 )

        end IF


        call wdeallocate_2z(zSqo)
        call wdeallocate_2z(VL)
        call wdeallocate_2z(VR)

        call wdeallocate_1z(zwork_loc)
        call wdeallocate_1d(work_locr)
     end if !(rank 0)
     !-------------------------------- !(info common to all processors) 
#ifdef MPI
     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end


     if (info==0) then

        !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 
#ifdef MPI
        call MPI_BCAST(zAq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(zSq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        call MPI_BCAST(fpm(37),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
#endif
        !-----------------------------

        fpm(21)=9
     end if
  end if !! fpm(21)=8



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors Xr=Qr*xr  Xl=Ql*xl 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==9 ) then
     call ZCOPY( N*fpm(23), Q(1,1), 1, workc, 1 ) 
     !! option - non shifted
#ifdef MPI 
     Q(1:N,1:fpm(23))=ZEROC
#endif
     call ZGEMM('N','N',N,fpm(25),fpm(23),ONEC,workc(1,1),N,zAq(1,fpm(24)),M0,ZEROC,Q(1,fpm(24)),N)

     if(fpm(31)/=2)  then
        call ZCOPY( N*fpm(23), Q(1,M0+1), 1, workc, 1 ) 
#ifdef MPI 
        Q(1:N,M0+1:M0+fpm(23))=ZEROC
#endif
        call ZGEMM('N','N',N,fpm(35),fpm(23),ONEC,workc(1,1),N,zSq(1,fpm(24)),M0,ZEROC,Q(1,fpm(34)),N)
     endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residuals ||AXr-lambda*BXr|| ||A^H Xl-conjg(lambda)*B^H Xl|| 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     fpm(21)=10 ! preparing reentry
     ijob=30 
     return  ! mat-vec A*qr => workr
  endif

  if (fpm(21)==10) then
     call ZLACPY( 'F', N, fpm(25),work(1,fpm(24)), N, workc(1,fpm(24)), N ) ! workc = Ax
     fpm(21)=11 ! preparing reentry
     ijob=40 
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZEROC  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     return  ! mat-vec S*qr => workr
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !-------- Absolute residual
     do i=fpm(24),fpm(24)+fpm(25)-1
        res(i)=sum(abs(workc(1:N,i)-lambda(i)*work(1:N,i)))/sum(abs((abs(Emid)+r)*work(1:N,i)))
     end  do
     !----------------------------------------
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------
     fpm(21)=-9 
     if(fpm(31)==2) fpm(21)=-12
  end if


  if (fpm(21)==-9) then
     fpm(21)=-10 ! preparing reentry
     ijob=31 
     return  ! mat-vec A^C.ql => work(1,fpm(34))
  endif

  if (fpm(21)==-10) then
     call ZLACPY( 'F', N, fpm(35),work(1,fpm(34)), N, workc(1,fpm(24)), N ) ! workc = A^C.XL
     fpm(21)=-11 ! preparing reentry
     ijob=41 
     !----------------------------------------
#ifdef MPI
     work(1:N,M0+1:M0+fpm(23))=ZEROC  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     return  ! mat-vec S^C.ql => work
  endif

  if (fpm(21)==-11) then
     !----------------------------------------
#ifdef MPI
     res(M0+1:M0+fpm(23))=DZERO
#endif
     !-------- Absolute residual
     do i=fpm(34),fpm(34)+fpm(35)-1
        res(i)=sum(abs(workc(1:N,i-M0)-conjg(lambda(i-M0))*work(1:N,i)))/sum(abs((abs(Emid)+r)*work(1:N,i)))
     end  do
     !----------------------------------------
#ifdef MPI 
     if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,NEW_COMM_WORLD,code)
#endif
     !-----------------------------------------
     fpm(21)=-12
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==-12) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !---color mapping --> check if eigenvalue inside the contour
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call wallocate_1i(color,fpm(23),infoloc)

     if (fpm(29)==0) then ! expert custom contour
        call zfeast_inside_contourx(Zne,fpm(8),lambda,fpm(23),color)
     else
        call zfeast_inside_contour(Emid,r,fpm(18),fpm(19),lambda,fpm(23),color)
     endif
!!!!!!!!!!!!

     k=0
     trace=ZEROC
     theta=DZERO
!!! count how many eigenvalues have converged + trace + max residual
     do i=1,fpm(23)
        if (color(i)==1) then ! inside the search interval
           k=k+1
           trace=trace+lambda(i)
           if (res(i)>theta) theta=res(i) ! max residual (right residual)
           if (fpm(31)/=2) then
              if (res(M0+i)>theta) theta=res(M0+i) ! max residual (left residual)
           end if
        endif
     enddo
!!!!!!!!! remove spurious if any  ---- (work only with the right subspace) ----
     if ((mode==0).and.(k>0)) mode=k ! wait before looking into mode 
     ! Rq: if all eigenvalues k are spurious...FEAST will not converge
     if (mode<k) then
        do j=1,k-mode
           theta=DZERO
           e=1
           do i=1,fpm(23)
              if (color(i)==1) then ! inside the search interval
                 if (res(i)>theta) then
                    e=i
                    theta=res(i) !max
                 endif
              end if
           enddo
           trace=trace-lambda(e) ! update trace
           res(e)=-DONE !spurious
           if (fpm(31)/=2) res(e+M0)=-DONE !spurious
        end do
!!! max residual
        theta=DZERO
        do i=1,fpm(23)
           if (color(i)==1) then ! inside the search interval
              if (res(i)>theta) theta=res(i) ! max residual (right residual)
              if (fpm(31)/=2) then
                 if (res(M0+i)>theta) theta=res(M0+i) ! max residual (left residual)
              end if
           end if
        enddo
     end if

     if (mode==0) info=1  ! no eigenvalue detected in the interval
     if (loop>1) then! wait second iteration (spurious related)
        if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
     endif
     if (info/=0) then
        fpm(21)=100 ! The End
        call wdeallocate_1i(color)
     end if

  end IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check! FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  if (fpm(21)==-12) then
     testconv=.false. ! initialization

     ! here Ze contains the old trace         
     if (loop==0) epsout=DONE
     if (loop>0) then !! trace convergence
        epsout=abs(trace-Ze)/(abs(Emid) + r) 
        if ((fpm(6)==0).and.(log10((epsout))<(-fpm(3)))) testconv=.true.
     end if

     !! residual convergence
     if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.

     !! Two loops minimum if spurious are found
     if ((loop<=1).and.(k>mode)) testconv=.false.


     if (fpm(1)==1) then
        call wwrite_i(fout,loop)
        call wwrite_t(fout)
        call wwrite_i(fout,mode)
        call wwrite_t(fout)
        call wwrite_d(fout,trace)
        call wwrite_t(fout)
        call wwrite_d(fout,epsout)
        call wwrite_t(fout)
        call wwrite_d(fout,theta)
        call wwrite_n(fout)
     end if

     if (.not.testconv) then
        ZSq(1,1)=trace ! dummy
        if (loop==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FEAST exit IF Convergence - FEAST iteration IF NOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

  if (fpm(21)==-12) then

     if (.not.testconv) then !!! need FEAST refinement loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! rational function for the inner eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        res(1:fpm(23))=DZERO ! temp indicator for being inside the contour
        call wallocate_1z(zwork_loc,fpm(23),infoloc)
        call zfeast_grationalx(Zne,Wne,fpm(8),lambda,fpm(23),zwork_loc)
        do i=1,fpm(23)
           if (color(i)==1)  res(i)=DONE
           lambda(i)=zwork_loc(i)
           !if (fpm(31) == 1) lambda(i) = (2.0d0,0.0d0)*dble(zwork_loc(i))
        enddo
        call wdeallocate_1z(zwork_loc)
        call wdeallocate_1i(color)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Remark: q and work=S*q is known already,  using FEAST-MPI work and q are already distributed
        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    ! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        
#ifdef MPI       
        if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        if (fpm(31)/=2) then
         if (fpm(23)/nb_procs>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,NEW_COMM_WORLD,code)
        endif
#endif

!! reorder (shift) lambda, eigenvector and residual
          call ZCOPY( N*fpm(23),Q(1,1), 1, work(1,1), 1 )
          if(fpm(31)/=2) call ZCOPY( N*fpm(23),Q(1,M0+1), 1, work(1,fpm(23)+1), 1 )
          call wallocate_1z(zwork_loc,fpm(23),infoloc)
          if (fpm(31)/=2) then
          call wallocate_1d(work_loc2,2*fpm(23),infoloc) 
          else
          call wallocate_1d(work_loc2,fpm(23),infoloc) 
          end if
          call ZCOPY(fpm(23),lambda , 1, zwork_loc, 1 )
          call DCOPY(fpm(23),res , 1, work_loc2, 1 )
          if(fpm(31)/=2) call DCOPY(fpm(23), res(M0+1), 1, work_loc2(fpm(23)+1), 1 )

          M0=fpm(23)  ! update value of M0 (new subspace)

          Q(1:N,1:fpm(23))=ZEROC
          if(fpm(31)/=2) Q(1:N,M0+1:M0+fpm(23))=ZEROC
          e=0
          j=0
          k=0
          do i=1,fpm(23)  
          if (color(i)==1) then ! inside the search interval
             if (work_loc2(i)/=-DONE) then ! not spurious 
                    e=e+1
                    q(1:N,e)=work(1:N,i)
                    lambda(e)=zwork_loc(i)
                    res(e)=work_loc2(i)
                    if(fpm(31)/=2) then
                    q(1:N,M0+e)=work(1:N,M0+i)
                    res(M0+e)=work_loc2(M0+i)
                 endif
              endif
           else
              j=j+1
              q(1:N,mode+j)=work(1:N,i)
              lambda(mode+j)=zwork_loc(i)
              res(mode+j)=work_loc2(i)
              if(fpm(31)/=2) then
                 q(1:N,M0+mode+j)=work(1:N,M0+i)
                 res(M0+mode+j)=work_loc2(M0+i)
              endif
           end if
           if (work_loc2(i)==-DONE) then ! spurious at the end
              k=k+1
              q(1:N,fpm(23)-k+1)=work(1:N,i)
              lambda(fpm(23)-k+1)=zwork_loc(i)
              res(fpm(23)-k+1)=work_loc2(i)
              if(fpm(31)/=2) then
                 q(1:N,M0+fpm(23)-k+1)=work(1:N,M0+i)
                 res(M0+fpm(23)-k+1)=work_loc2(M0+i)
              end if
           endif
        end do

        call wdeallocate_1z(zwork_loc)
        call wdeallocate_1d(work_loc2)
        call wdeallocate_1i(color)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=100 ! The End

        if ((info==0).and.(fpm(37)==1)) info=6

      end if ! test convergence
    end if 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (fpm(21)==100) then !! THE END (ijob=0) 

       ijob=0 !! exit FEAST

       if (fpm(1)==1) then !! Print  Information

          if (info>=200) then
             call wwrite_s(fout, 'PROBLEM with input parameters')
             call wwrite_n(fout) 
          end if

          if ((info>100).and.(info<200)) then
             call wwrite_s(fout, 'PROBLEM with FEAST array parameters')
             call wwrite_n(fout) 
          end if

          if (info==-3) then
             call wwrite_s(fout, 'ERROR with reduced system')  
             call wwrite_n(fout) 
          end if

          if (info==-2) then
             call wwrite_s(fout, 'ERROR from Inner Linear System Solver in FEAST driver')  
             call wwrite_n(fout) 
          end if

          if (info==-1) then
             call wwrite_s(fout, 'ERROR with Internal memory allocation')  
             call wwrite_n(fout) 
          end if

          if (info==1) then
             call wwrite_s(fout, '==>WARNING: No eigenvalue has been found in the proposed search interval')
             call wwrite_n(fout)
          endif

          if (info==3) then
             call wwrite_s(fout, '==>WARNING: Size subspace M0 too small')  
             call wwrite_n(fout)
          end if

          if (info==4) then
             call wwrite_s(fout, '==>WARNING: Only the subspace has been returned')  
             call wwrite_n(fout)
          end if

           if (info==5) then
              call wwrite_s(fout, '==>WARNING: Only stochastic estimation of #eigenvalues returned')  
              call wwrite_n(fout)
           end if

          if (info==6) then
             call wwrite_s(fout, '==>WARNING: FEAST converge but subspace is not biorthonormal')
             call wwrite_n(fout) 
          end if

          if (info==2) then
             call wwrite_s(fout, '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)')  
             call wwrite_n(fout)
          end if

          if (info==0) then
             call wwrite_s(fout, '==>FEAST has successfully converged (to desired tolerance)')  
             call wwrite_n(fout) 
          else
             call wwrite_s(fout, '==>INFO code = ') 
             call wwrite_i(fout,info)
             call wwrite_n(fout)
          end if


          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_s(fout, '*********** FEAST- END*************************')
          call wwrite_n(fout) 
          call wwrite_s(fout, '***********************************************')  
          call wwrite_n(fout) 
          call wwrite_n(fout)
       endif
    end if

  end subroutine zfeast_grcix




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_grci(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) 
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      Aqr=lambda Bqr and A^Hql=conjg(lambda) B^Hql  within a given search interval                          
  !  
  !  A and B are considered COMPLEX NON-HERMITIAN   
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,2*M0):   Workspace
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input/output) COMPLEX DOUBLE PRECISION:  middle of the search interval
  !  r          (input/output) REAL DOUBLE PRECISION: horizontal radius of the search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,2*M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)
  !                                                                 Left  Eigenvectors-solution (1:N,M0+1:M0+mode)
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right (1:mode) and left (M0+1:M0+mode)
  !                          
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
  !=====================================================================
  !  James Kestyn - Eric Polizzi - 2015
  !=====================================================================

    implicit none
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))) :: Emid
    double precision :: r
    complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
    complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
    integer,dimension(*) :: fpm
    double precision :: epsout
    integer :: loop
    complex(kind=(kind(1.0d0))),dimension(*):: lambda
    complex(kind=(kind(1.0d0))),dimension(N,*):: Q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine to expert routine: zfeast_grcix
!!!    Create Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)
    
    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    fpm(31)=0  ! ID for Complex Non-Symmetric Matrix
    call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)

end subroutine zfeast_grci

  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine zfeast_srcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, and right qr 
  !      Aqr=lambda Bqr (Remark ql=conjg(qr))                            
  ! 
  !  A and B are COMPLEX SYMMETRIC 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (output)       COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (output)       REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input) Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)                                                               
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right residual(1:mode) 
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(8)): Custom Integration weights
  !
  !=====================================================================
  !  James Kestyn - Eric Polizzi - 2015
  !=====================================================================

    implicit none
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))) :: Emid
    double precision :: r
    complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
    complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
    integer,dimension(*) :: fpm
    double precision :: epsout
    integer :: loop
    complex(kind=(kind(1.0d0))),dimension(*):: lambda
    complex(kind=(kind(1.0d0))),dimension(N,*):: Q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info
    complex(kind=(kind(1.0d0))),dimension(*):: Zne,Wne


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine
!!!    Call expert routine zfeast_grcix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    fpm(31)=2 ! ID for Complex Symmetric Matrix
    call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)

end subroutine zfeast_srcix

  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine zfeast_srci(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) - 
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, and right qr 
  !      Aqr=lambda Bqr (Remark ql=conjg(qr))  within a search interval                        
  ! 
  !  A and B are COMPLEX SYMMETRIC 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input)        COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input)        REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input) Emid,r provided by user
  !                                      
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)                                                               
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right residual(1:mode) 
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
  !=====================================================================
  !  James Kestyn - Eric Polizzi - 2015
  !=====================================================================

    implicit none
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))) :: Emid
    double precision :: r
    complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
    complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
    integer,dimension(*) :: fpm
    double precision :: epsout
    integer :: loop
    complex(kind=(kind(1.0d0))),dimension(*):: lambda
    complex(kind=(kind(1.0d0))),dimension(N,*):: qr
    integer :: mode
    double precision,dimension(*) :: resr
    integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine
!!!    Create Zne, Wne arrays
!!!    Call expert routine zfeast_grcix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)
     
    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    fpm(31)=2 ! ID for Complex Symmetric Matrix
    call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info,Zne,Wne)

end subroutine zfeast_srci

  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine dfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, and right qr 
  !      Aqr=lambda Bqr (Remark ql=conjg(qr))                            
  ! 
  !  A and B are REAL NON-SYMMETRIC 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (output)       COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (output)       REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input) Emid,r dummy variables
  !                                      (output) Emid,r calculated/estimated from custom contour
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)                                                               
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right residual(1:mode) 
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of ARGUMENTS FOR EXPERT ROUTINE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Zne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration nodes
  !  Wne        (input)        COMPLEX DOUBLE PRECISION(fpm(2)): Custom Integration weights
  !
  !=====================================================================
  !  James Kestyn - Eric Polizzi - 2015
  !=====================================================================

    implicit none
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))) :: Emid
    double precision :: r
    complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
    complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
    integer,dimension(*) :: fpm
    double precision :: epsout
    integer :: loop
    complex(kind=(kind(1.0d0))),dimension(*):: lambda
    complex(kind=(kind(1.0d0))),dimension(N,*):: Q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info
    complex(kind=(kind(1.0d0))),dimension(*):: Zne,Wne


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine
!!!    Call expert routine zfeast_grcix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fpm(31)=3 !ID for Real Non-Symmetric matrix
    call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)

end subroutine dfeast_grcix

  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine dfeast_grci(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) - 
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, and right qr 
  !      Aqr=lambda Bqr (Remark ql=conjg(qr))  within a search interval                        
  ! 
  !  A and B are REAL NON-SYMMETRIC 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zSq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm     (input/output)    INTEGER(*) : FEAST parameters (see FEAST documentation)
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emid       (input)        COMPLEX DOUBLE PRECISION: middle of the search interval
  !  r          (input)        REAL DOUBLE PRECISION: horizontal radius of the search interval
  !                            comments: (input) Emid,r provided by user
  !                                      
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       COMPLEX DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Right Eigenvectors-solution (1:N,1:mode)                                                               
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right residual(1:mode) 
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
  !=====================================================================
  !  James Kestyn - Eric Polizzi - 2015
  !=====================================================================

    implicit none
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))) :: Emid
    double precision :: r
    complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
    complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zSq
    integer,dimension(*) :: fpm
    double precision :: epsout
    integer :: loop
    complex(kind=(kind(1.0d0))),dimension(*):: lambda
    complex(kind=(kind(1.0d0))),dimension(N,*):: Q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    Wrapper Routine
!!!    Create Zne, Wne arrays
!!!    Call expert routine zfeast_grcix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 

    call zfeast_gcontour(Emid,r,fpm(8),fpm(17),fpm(18),fpm(19),Zne,Wne)

    fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
    fpm(31)=3 ! ID for Real Non-Symmetric matrix 

    call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
    

end subroutine dfeast_grci




