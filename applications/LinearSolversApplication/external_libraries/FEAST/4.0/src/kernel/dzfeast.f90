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
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions l
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
! dfeast_srcix   ! Double Real Symmetric - expert ==>kernel<==
! dfeast_srci    ! Double Real Symmetric
!
! zfeast_hrcix   ! Double Complex Hermitian - expert ==> kernel <==
! zfeast_hrci    ! Double Complex Hermitian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Non-Symmetric (real) /Non-Hermitian (complex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! zfeast_grcix   ! Double Complex General - expert  ==>kernel<==
! zfeast_grci    ! Double Complex General
! 
! zfeast_srcix   ! Double Complex Symmetric - expert
! zfeast_srci    ! Double Complex Symmetric
!
! dfeast_grcix   ! Double Real General - expert
! dfeast_grci    ! Double Real General


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Polynomial (real)/(complex)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! zfeast_grcipevx ! Double Complex General - expert  ==>kernel<==
! zfeast_grcipev  ! Double Complex General
! 
! zfeast_srcipevx ! Double Complex Symmetric - expert
! zfeast_srcipev  ! Double Complex Symmetric








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RCI ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfeast_srcix(ijob,N,Ze,work,workc,Aq,Bq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)
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
  !  Aq,Bq      (input/output) REAL DOUBLE PRECISION (M0,M0)  :  Workspace for Reduced Eigenvalue System
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
  ! Eric Polizzi 2009-2019 -- James Kestyn 2016-2018 (pfeast)
  ! ====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  !include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  double precision, dimension(N,*) ::work
  complex(kind=(kind(1.0d0))),dimension(N,*):: workc
  integer,dimension(*) :: fpm
  double precision,dimension(M0,*):: Aq,Bq
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
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZZERO=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0
  ! integer(8),parameter :: fout =6
  integer(8),save :: fout
  logical,save :: com
  !! variable for FEAST
  character(len=3) :: ctemp 
  character(len=25) :: name
  integer :: i,e,j,k,Ntotal
  integer,dimension(4) :: iseed
  double precision :: theta,r
  complex(kind=(kind(1.0d0))) :: zxe,zwe,aux
  double precision, dimension(:,:),allocatable :: Bqo
  integer, dimension(:),allocatable :: fpm_default
  logical :: testconv,fact
  double precision :: trace
  !! Lapack variable (reduced system)
  character(len=1) :: JOBZ,UPLO
  double precision, dimension(:),allocatable :: work_loc,work_loc2
  integer :: lwork_loc,info_lap,infoloc
  !! MPI compatibility variables
  integer :: rank2,rank3,code,nb_procs2,nb_procs3,L2_COMM_WORLD,L3_COMM_WORLD
  double precision :: DNRM2
  !! timing
  integer,save :: tini,t1,t10,t11,t30,t40,t50,t60
  integer :: tend,tim,t2
  character(len=8) :: date
!  double precision, dimension(:,:),allocatable :: dwork


  !  if (ijob/=-1) then
!!! reopen log file !!<<
  !     if (fpm(1)<0) then
  !      close(fout)  
  !      write(ctemp,'(I3)') abs(fpm(1))
  !      open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append')       
  !   endif
  !endif

!!! record timing
  if (ijob/=-1) then  
     call system_clock(t2,tim)
     select case(ijob)
     case(10)
        t10=t10+(t2-t1)
     case(11)
        t11=t11+(t2-t1)
     case(30)
        t30=t30+(t2-t1)
     case(40)
        t40=t40+(t2-t1)
     case(50)
        t50=t50+(t2-t1)
     case(60)
        t60=t60+(t2-t1)
     end select
  end if

!!! parallel init
  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1
  !----------------------------------------------
#ifdef MPI
  L2_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm(49)
  if (fpm(49)/=0) then
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  end if
#endif
  !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then

     info=0 ! default value
     if (fpm(30)==-111) then
        fpm(30)=121111 ! code name for direct call to this routine   
        call feastdefault(fpm,info) 
     endif
     if (info/=0) then
        ijob=0 ! leave rci
        return
     endif

     call system_clock(tini,tim)
     t10=0
     t11=0
     t30=0
     t40=0
     t50=0
     t60=0

     ! set up factorization id
     fpm(33)=0

     ! comments?
     com=.false. !default
     if (fpm(1)/=0) com=.true.
#ifdef MPI
     if (com.and.((rank2/=0).or.(rank3/=0))) com=.false. ! comment only in rank 0 (only a single mpi is allowed to comment)
#endif
     !file or screen output?
     if (fpm(1)<0) then
        write(ctemp,'(I3)') abs(fpm(1))
        fout=abs(fpm(1))+200!fpm(60) ! file number default
        if (com) then
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append') !! all mpi procs may have access to the file (if needed)
           !write(fout,'(A)',advance='no') '############## NEW FEAST RUN ###################'
           !write(fout,*)
        endif
     elseif (fpm(1)>0) then
        fout=6 !screen
     end if     

!!!!!!!!!!!!!!!! Find Ntotal
     Ntotal=N
#ifdef MPI           
     if (fpm(49)/=0) then
        Ntotal=0
        call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     end if
#endif
!!!!!!!!!!!!!!!!!
     fpm(26)=Ntotal ! save value (used for stochastic estimate)


     call dcheck_feast_srci_input(Emin,Emax,M0,Ntotal,info)

     if (info/=0) fpm(21)=100 ! The End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF (info==0) then


        if (com) then
           write(fout,*)
           write(fout,'(A)',advance='no') '***********************************************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '*********** FEAST v'
           write(fout,'(F3.1)',advance='no') fpm(31)/10.0d0
           write(fout,'(A)',advance='no') ' BEGIN ******************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '***********************************************'  
           write(fout,*)
        if (fpm(1)<0) then ! only in-file
           call date_and_time(DATE=date)
           write(fout,'(A)',advance='no') 'Date: '//trim(date(1:4))//'-'//trim(date(5:6))//'-'//trim(date(7:8))
           write(fout,*)
        end if
           call feast_name(fpm(30),name)
           write(fout,'(A)') 'Routine '//trim(name)
           if ((mod(fpm(30),10)==2).or. (mod(fpm(30),10)==3)) then
                write(fout,'(A)',advance='no') 'Solving AX=eX with A real symmetric'
           elseif    ((mod(fpm(30),10)==4).or. (mod(fpm(30),10)==5)) then 
              write(fout,'(A)',advance='no') 'Solving AX=eBX with A real symmetric and B spd'
           end if
           write(fout,*)
           if (nb_procs2*nb_procs3>1) then
              write(fout,'(A)',advance='no') '#MPI (total=L2*L3) '
              write(fout,'(I4)',advance='no') nb_procs2*nb_procs3
              write(fout,'(A)',advance='no') '='
              write(fout,'(I4)',advance='no') nb_procs2
              write(fout,'(A)',advance='no') '*'
              write(fout,'(I4)',advance='no') nb_procs3
              write(fout,*)
           endif
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
           write(fout,'(A)',advance='no') 'List of input parameters fpm(1:64)-- if different from default'
           write(fout,*)
           allocate(fpm_default(64))
           call feastinit(fpm_default) ! initialize to -111
           fpm_default(1)=0
           fpm_default(30)=fpm(30) ! name
           fpm_default(14)=fpm(14) ! for special conditional cases
           call feastdefault(fpm_default,infoloc) ! get the default values
           do i=1,64   ! compare with possible mismatch with user values
              if ((i<=19).or.((i>=40).and.(i<=50))) then
                 if ((fpm(i)/=fpm_default(i)).and.(i/=9).and.(i/=49)) then
                    write(fout,'(A)',advance='no') '   fpm('
                    write(fout,'(I2)',advance='no') i
                    write(fout,'(A)',advance='no') ')='
                    write(fout,'(I4)',advance='no') fpm(i)
                    write(fout,*)
                 endif
              endif
           enddo
           deallocate(fpm_default)
           fpm(22)=nb_procs3 ! temp copy
           call dfeast_info(fout,fpm,Emin,Emax,Ntotal,M0)

        end if

!!!!!!!!!!!!!!!!!!!!!!!! 


        if (com) then
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST runs        |'
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           if (fpm(14)/=2) then
              write(fout,'(A)',advance='no') '#It |  #Eig  |          Trace            |     Error-Trace          |     Max-Residual' 
           else
              write(fout,'(A)',advance='no') 'Running average for stochastic estimation (1 -> M0)'  
           endif
           write(fout,*)
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
        endif

!!!!!!!!!!! Some initialization
        fpm(22)=fpm(2) ! only one half contour necessary
        loop=0


        fpm(23)=min(M0,fpm(26)) ! 'current M0' size (global value)
        if (fpm(14)==2) fpm(23)=min(fpm(32),fpm(26),M0) ! stochastic estimate
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs2>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(24)=1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=1 ! prepare reentry

        if (fpm(5)==0) then !!! random vectors (option 2 in DLARNV)
           iseed=(/56,890,3456,2333/)
           iseed(1) = iseed(1) * (1+rank3*13)
           iseed(2) = iseed(2) + rank3*7
           iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
           iseed(4) = iseed(4) + rank3*2


           ! copy work into q for multiplication by B matrix, if stochastic estimate  is "on"
           if (fpm(14)==2) then
              ! random distribution must be uniform for stochastic to work
              call DLARNV(2,iseed,N*fpm(23),q(1,1))             
              !          Ntotal=fpm(26)
              ! allocate(dwork(Ntotal,M0))
              ! call DLARNV(2,iseed,Ntotal*fpm(23),dwork(1,1))
              ! if (rank3==0) q(1:N,1:M0)=dwork(1:N,1:M0)
              ! if (rank3==1) q(1:N,1:M0)=dwork(2917:Ntotal,1:M0) 
              ! deallocate(dwork)
           else              
              call DLARNV(2,iseed,N*fpm(23),work(1,1)) !random distribution
              ! Ntotal=fpm(26)
              ! allocate(dwork(Ntotal,M0))
              ! call DLARNV(2,iseed,Ntotal*fpm(23),dwork(1,1))
              ! if (rank3==0) work(1:N,1:M0)=dwork(1:N,1:M0)
             !   if (rank3==1) work(1:N,1:M0)=dwork(5:Ntotal,1:M0)
              ! if (rank3==1) work(1:N,1:M0)=dwork(2917:Ntotal,1:M0) 
              ! deallocate(dwork)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !           startj=1
              !#ifdef MPI                     
              !call wallocate_1i(Nsize,nb_procs3,infoloc)
              !Nsize(1:nb_procs3)=0
              !Nsize(rank3+1)=N !local
              !call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nb_procs3,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)!global
              !startj=1
              !if (rank3>0) startj=1+sum(Nsize(1:rank3))
              !deallocate(Nsize)
              !#endif
              !do i=1,fpm(23)
              !    do k=1,N
              !       q(k,i)=cos(-pi*sqrt(3.0d0)+i*((startj+k-2)*(2.0d0*pi*sqrt(3.0d0)/(Ntotal-1))))+sin(-pi*sqrt(3.0d0)+i*((startj+k-2)*(2.0d0*pi*sqrt(3.0d0)/(Ntotal-1))))
              !    end do
              !enddo    

              !! call DLARNV(2,iseed,fpm(23)*fpm(23),Bq)
              !! do i=1,fpm(23)
              !!    do k=1,N
              !!       q(k,i)=cos(-pi+i*((startj+k-2)*(2.0d0*pi/(Ntotal-1))))
              !!    end do
              !!enddo  
              !!call DGEMM('N','N',N,fpm(23),fpm(23),DONE,q,N,Bq,M0,DZERO,work,N) ! projection

           endif
        end if

        !! Rq: if fpm(14)=2 or fpm(5)=1, q is the initial guess
        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=DZERO 
#endif
           !------------------------------------------
           if (fpm(5)==1) then ! Compute residual vector R=AX-BXE 
              ijob=30 !!A*q==>work
              fpm(21)=10
           else ! only stochastic estimate
              ijob=40 !! B*q=>work
           endif
           call system_clock(t1,tim)
           return
        endif

     end IF ! info=0

!!!!!!!!!!!!!!
  end if   !ijob=-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==1) then !! we initialize a new contour integration
     !------------------------------------------------------------------------
#ifdef MPI

     if ((loop>0).or.(fpm(14)==2)) then !!fpm(5) condition unncessary!!<<
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
     endif

     if ((loop>0).and.(fpm(23)/nb_procs2>=1)) then
        if (.not.((loop==1).and.(fpm(5)==1))) then ! do not include first iteration with initial guess since Q is known to all
           call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
        end if
     end if
#endif
     !------------------------------------------------------------------------

!!!!!!!!!!! RESIDUAL INVERSE ITERATION 
     if (loop==0) then
        !if ((fpm(29)==1).and.(fpm(16)==2).and.(fpm(5)==1)) then !! Zolo and fpm(5)=1
        !   call dset_zolotarev(fpm(2),0,zxe,zwe) ! get zwe
        !   if (rank2==0) then
        !      q(1:N,1:fpm(23))=-dble(zwe)*q(1:N,1:fpm(23))
        !   else 
        !      q(1:N,1:fpm(23))=DZERO
        !   endif
        !else 
        q(1:N,1:fpm(23))=DZERO
        !end if
     else if (loop>0) then ! general case, works for Zolo as well..res(i) includes the offset
        do i=1,fpm(23)
           if (res(i)>100d0) then
              call DSCAL(N,-(res(i)-200d0), q(1,i), 1)
           else
              call DSCAL(N,-res(i), q(1,i), 1)
           end if
        end do
        if (rank2/=0)  q(1:N,1:fpm(23))=DZERO ! for parallelism
     end if
!!!!!!!!!!!!!!

     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     ! IF (info==0) then !! will end up checking info errors returned by FEAST drivers
     do e=fpm(20)+rank2,fpm(22),nb_procs2 !!!! loop over the contour 

        if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
           ! set up factorization id (if fact. storage is used)
           if (fpm(10)==1) then
              fpm(33)=fpm(33)+1
           else
              fpm(33)=1
           endif

           Ze=Zne(e)
           fpm(20)=e-rank2

           ijob=10 ! for fact
           fact=.false. ! action to factorize
           if ((loop==0).or.((loop==1).and.(fpm(5)==1))) then
              fact=.true. ! always factorize in initial loop. Remark: with input initial guess, loop=1 is the initial loop
           else !(not initial loop)
              if ((fpm(22)>nb_procs2).and.(fpm(10)==0)) fact=.true. ! re-factorize if #mpi<#contour points and fact not stored
           end if
           if (fact) then
              call system_clock(t1,tim)
              return 
           end if

!!$           if ((loop==0).or.(fpm(22)>nb_procs2)) then !no refactorization if one linear system per processor
!!$              if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) then !no refactorization if stored in memory
!!$                 if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) then 
!!$                    call system_clock(t1,tim)
!!$                    return 
!!$                 endif
!!$              end if
!!$           end if


        endif

!!!!!!!!!!!!!!       
        if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v
!!!!!!!!!!!!!!!

           call ZLACP2( 'F', N, fpm(23),work , N, workc, N )   
           ijob=11 ! for solve
           call system_clock(t1,tim)
           return
        endif

!!!!!!!!!!!!!!!!!!!           
        if (ijob==11) then
!!!!!!!!!!!!!!!!!!              
           !! summation              
           aux=2.0d0*Wne(e)            
           !! aux has been multiplied by 2 to account for the 2 half-contours             
           if (loop>0) then ! for residual inverse iteration
              do i=1,fpm(23)
                 call ZSCAL(N,ZONE/(ze-lambda(i)*ZONE), workc(1,i), 1)
              enddo
           endif

           q(1:N,1:fpm(23))=q(1:N,1:fpm(23))-dble(aux*workc(1:N,1:fpm(23)))


           ijob=-2 ! just for identification
        end if
     end do
     !     end IF  !! info=0

     !------------------------------------------------
     !#ifdef MPI
     !call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
     !todo: 2 allreduce and reset the value at -2
     !#endif
     !-----------------------------------------------  
     !     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(21)=4
        !print *,rank2,q(1:N,1:M0)

        !------------------------------------------------
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code) 
#endif
        !-----------------------------------------------
     end if
     !print *,rank2,q(1:N,1:M0)
  end IF     ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !print *,'here',sum(abs(q(1:N,1))),sum(abs(q(1:N,2)))
  !stop






  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic estimates
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     iseed(1) = iseed(1) * (1+rank3*13)
     iseed(2) = iseed(2) + rank3*7
     iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
     iseed(4) = iseed(4) + rank3*2


     call DLARNV(2,iseed,N*fpm(23),work(1,1))

     Ntotal=fpm(26)

!!!!!!!!!!!!!!!!!     
     !allocate(dwork(Ntotal,M0))
     !call DLARNV(2,iseed,Ntotal*fpm(23),dwork(1,1))
     !if (rank3==0) work(1:N,1:M0)=dwork(1:N,1:M0)
     !if (rank3==1) work(1:N,1:M0)=dwork(N:Ntotal,1:M0) !system3/1
     !deallocate(dwork)


     call DGEMM('T','N',fpm(23),fpm(23),N,-DONE,work,N,q,N,DZERO,Aq,M0) ! projection
     call DGEMM('T','N',fpm(23),fpm(23),N,DONE,work,N,work,N,DZERO,Bq,M0) ! normalization

     !----------------------------------------
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,Bq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,Aq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !---------------------------------------
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+Aq(i,i)/Bq(i,i)
        res(i)=abs(theta*(Ntotal/(DONE*i)))
        if (com) then
           write(fout,'(I3)',advance='no') i
           write(fout,'(3X)',advance='no')
           write(fout,'(F10.2)',advance='no') res(i)
           write(fout,*)
        endif
     enddo
     mode=int(real(res(fpm(23))))
     if (com) then
        write(fout,'(A)',advance='no') '# estimated '
        write(fout,'(3X)',advance='no')
        write(fout,'(I6)',advance='no') mode
        write(fout,*)
     endif
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!

!!$ if (fpm(21)==4) then 
!!$ fpm(21)=81
!!$if (fpm(13)==-1) then ! Form and solve the reduced system
!!$   ijob=50
!!$   return
!!$   endif
!!$endif


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
     call system_clock(t1,tim)
     return! mat-vec B*q => work
  end if


  if (fpm(21)==5) then
     !----------------------------------------------
#ifdef MPI
     Bq(1:M0,1:fpm(23))=DZERO
#endif
     !-------------------------------------------------
     fpm(21)=51

     !      print *,'work',rank2,work(1:N,1:M0)

     if (fpm(13)>=2) then ! customize inner product
        ijob=50
        call system_clock(t1,tim)
        return
     endif

     call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Bq(1,fpm(24)),M0) 
  end if


!!!!!!!!!!!!!!  
  if (fpm(21)==51) then
     if (fpm(13)>=2) call DCOPY(M0*fpm(25),Aq(1,fpm(24)),1,Bq(1,fpm(24)),1)

     !---------------------------------------- !(Bq known to all processors) 
#ifdef MPI
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,Bq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
     end if
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,Bq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !---------------------------------------
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==51) then
     mode=0
     if (loop>0) then
        do i=fpm(24),fpm(24)+fpm(25)-1
           if (res(i)>100d0) then ! indicator for being inside the search interval
              if (abs(sqrt(abs(Bq(i,i)))-(res(i)-200.0d0))/sqrt(abs(Bq(i,i)))<0.1d0) mode=mode+1
           endif
        enddo
     endif
!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,L2_COMM_WORLD,code)
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
     call system_clock(t1,tim)
     return  ! mat-vec A*q => work
  endif


  if (fpm(21)==7) then 
     !------------------------------------------------
#ifdef MPI 
     Aq(1:M0,1:fpm(23))=DZERO
#endif
     !-------------------------------------------------
     fpm(21)=8
     if (fpm(13)>=2) then ! customize inner product
        ijob=50
        call system_clock(t1,tim)
        return
     endif
     !if (rank2==1) print *,'work',work(1:N,1:M0)

     call DGEMM('T','N',fpm(23),fpm(25),N,DONE,q(1,1),N,work(1,fpm(24)),N,DZERO,Aq(1,fpm(24)),M0) !q^tAq
  endif


  if (fpm(21)==8) then
     !---------------------------------------- !(Aq known to all processors) 
#ifdef MPI
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,Aq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
     end if
     if(fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,Aq(1,1),M0*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !---------------------------------------
     !if (rank2==0) print *,'Bq',Bq(1:N,1:M0)
     !if (rank2==0) print *,'Aq',Aq(1:N,1:M0)



     !---------------------------------------
     if ((fpm(13)==1).or.(fpm(13)==3)) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=60
        call system_clock(t1,tim)
        return
     endif
     fpm(21)=81
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==81) then
     ! if using FEAST-MPI ==> solve on a single proc
     if ((rank2==0).and.(rank3==0)) then
        JOBZ='V'
        UPLO='L'
        info_lap=1 ! initialization
        i=1
        LWORK_LOC=3*fpm(23)-1 !! for lapack eig reduced system
        allocate(WORK_LOC(LWORK_LOC))

        do while ((info_lap/=0).and.(info==0))
           i=i+1
           if (i==10) info=-3 ! arbitrary maximum
           allocate(Bqo(fpm(23),fpm(23)))
           call DLACPY( 'F', fpm(23), fpm(23),Bq, M0, Bqo, fpm(23) )
           call DSYGV(1, JOBZ, UPLO, fpm(23),Aq,M0,Bqo,fpm(23),lambda,work_loc,Lwork_loc,INFO_lap)
           if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
           if (info_lap>fpm(23)) then !! Bqo is not spd (a posteriori resize subspace)
              fpm(23)=info_lap-fpm(23)-1 
              if (com) then
                 write(fout,'(A)',advance='no') 'Resize subspace'  
                 write(fout,'(3X)',advance='no') 
                 write(fout,'(I6)',advance='no') fpm(23)
                 write(fout,*)
              end if
           end if
           deallocate(Bqo)
        end do
        deallocate(work_loc) 
     end if !(rank 0)
     !-------------------------------- !(info common to all processors)
#ifdef MPI
     if (rank3==0) call MPI_BCAST(info,1,MPI_INTEGER,0,L2_COMM_WORLD,code)
     if(fpm(49)/=0) call MPI_BCAST(info,1,MPI_INTEGER,0,L3_COMM_WORLD,code)
#endif
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(25)=fpm(23)!! current M0 size (by default) -- global
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------!(Aq==> vectors, lambda and fpm(23), known to all processors) 
#ifdef MPI 
        if (rank3==0) call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,L2_COMM_WORLD,code)
        if (fpm(49)/=0) call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,L3_COMM_WORLD,code)

        if (rank3==0) then
           call MPI_BCAST(Aq,fpm(23)*M0,MPI_DOUBLE_PRECISION,0,L2_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,L2_COMM_WORLD,code)
        end if
        if (fpm(49)/=0) then
           call MPI_BCAST(Aq,fpm(23)*M0,MPI_DOUBLE_PRECISION,0,L3_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,L3_COMM_WORLD,code)
        end if

        !if (rank2==0) print *,'Aq',Aq(1:N,1:M0) 


        if (fpm(23)/nb_procs2>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs2 ! local size of current M0
           if (rank2==nb_procs2-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(24)=1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------
        fpm(21)=9
     end if
  end if !! fpm(21)=81

  ! print *,'aq',rank2,Aq(1:M0,1:M0)
  !   call MPI_BARRIER(L2_COMM_WORLD)
  !   call MPI_FINALIZE(code)
  !stop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors X=Qxq  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     !! from previous step work=A*Q, before using work memory space
     !! let us construct first workc=work*Aq which is workc=A*X
     !! Remark: workc is complex but we use only the memory space
     !------------------------------------------------------------------------!!!here
#ifdef MPI
     !  if (rank2==1) print *,'work',work(1:N,1:M0)

     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !------------------------------------------------------------------------
     !  if (rank2==1) print *,'workb',work(1:N,1:M0)
     call DGEMM('N','N',N,fpm(25),fpm(23),DONE,work(1,1),N,Aq(1,fpm(24)),M0,DZERO,workc(1,fpm(24)),N)

     !! construct X (q needed for outer-loop as well)  !!!!here
     call DLACPY( 'F', N, fpm(23),q(1,1) , N, work(1,1), N )

#ifdef MPI 
     q(1:N,1:fpm(23))=DZERO
#endif

     !  if (rank2==1) print *,'worka',work(1:N,1:M0)

     call DGEMM('N','N',N,fpm(25),fpm(23),DONE,work(1,1),N,Aq(1,fpm(24)),M0,DZERO,q(1,fpm(24)),N) 

     !if (rank2==1) print *,'q',q(1:N,1:M0)

     fpm(21)=10 !! no need to compute A*q again
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (fpm(21)==10) then
     !handle special case of input initial guess
     !  if ((fpm(5)==1).and.(loop==0)) call DLACPY( 'F', N, fpm(23),work , N, workc, N ) !using memory space of workc
     if ((fpm(5)==1).and.(loop==0))  call DLACPY( 'F', N, fpm(25),work(1,fpm(24)) , N, workc(1,fpm(24)), N )


     fpm(21)=11 ! preparing reentry
     ijob=40
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=DZERO   
#endif
     !----------------------------------------
     call system_clock(t1,tim)
     return  ! mat-vec B*q => work 
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !---------Compute residual vector R=AX-BXE
     !         work<=workc-work*E (needed also for outer-loop)
     do i=fpm(24),fpm(24)+fpm(25)-1
        !res(i)=sum(abs(work(1:N,i))) ! placeholder (denominator)
        res(i)=DNRM2(N,work(1,i),1)**2 ! square because of mpi compatibility sqrt(mpi1)+sqrt(mpi2) /= sqrt(mpi1+mpi2) 
        call DSCAL(N,-lambda(i),work(1,i),1)
     enddo
     call DAXPY(N*fpm(25),DONE,workc(1,fpm(24)),1,work(1,fpm(24)),1)

     !-------- Residual

     !reduce first denominator
#ifdef MPI 
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif

     do i=fpm(24),fpm(24)+fpm(25)-1
        !res(i)=sum(abs(work(1:N,i)))/(max(abs(Emin),abs(Emax))*res(i))
        res(i)=(DNRM2(N,work(1,i),1)**2)/((max(abs(Emin),abs(Emax))**2)*res(i))
     end do

     !----------------------------------------
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !-----------------------------------------
     do i=1,fpm(23)
        res(i)=sqrt(res(i)) ! norm L2
     enddo
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


     if (loop>2) then ! wait third iteration (spurious and ifeast related)
        if (mode==0) info=1  ! no eigenvalue detected in the interval
        if ((mode==M0).and.(mode/=N))  info=3 ! size subspace too small
     end if
     if (info/=0) fpm(21)=100 ! The End
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  if (fpm(21)==11) then
     testconv=.false. ! initialization

     !! trace convergence
     if (loop==0) epsout=DONE
     if (loop>0) then
        epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))
        !if ((fpm(6)==0).and.(log10(epsout)<(-fpm(3)))) testconv=.true. ! epsout could be 0
        if ((fpm(6)==0).and.(epsout<10.0d0**(-fpm(3))))  testconv=.true.
     end if

     !! residual convergence
     !if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.
     if ((fpm(6)/=0).and.(theta<10.0d0**(-fpm(3))))  testconv=.true.
     
     !! Two loops minimum if spurious are found
     if ((loop<=1).and.(k>mode)) testconv=.false.


!!!! L2 barrier to guarantee synchronization in printing
#ifdef MPI
 call MPI_BARRIER(L2_COMM_WORLD,code)
#endif 

     
     if (com) then
        if (fpm(5)==0) then
           write(fout,'(I3)',advance='no') loop
        else ! no integration in first loop
           if (loop==0) then
              write(fout,'(A)',advance='no') 'N/A'
           else
              write(fout,'(I3)',advance='no') loop-1
           end if
        endif
        write(fout,'(3X)',advance='no') 
        write(fout,'(I4)',advance='no') mode
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') trace
        write(fout,'(3X)',advance='no') 
        write(fout,'(ES25.16)',advance='no') epsout
        write(fout,'(3X)',advance='no') 
        write(fout,'(ES25.16)',advance='no') theta
        write(fout,*) 
     end if

     if (.not.testconv) then
        epsout=trace
        e=loop
        if (fpm(5)==1) e=e-1
        if (e==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif

     !! handle a special case (prevent convergence if no eigenvalue found)
     !! for mode=0 FEAST will exit after few iterations (see above)
     if (mode==0) testconv=.false. 

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
        allocate(WORK_LOC(fpm(23)))
        call dfeast_rationalx(Zne,Wne,fpm(2),lambda,fpm(23),work_loc)
        if ((fpm(29)==1).and.(fpm(16)==2)) then ! zolotarev case
           call dset_zolotarev(fpm(2),0,zxe,zwe)
           work_loc(1:fpm(23))=work_loc(1:fpm(23))+dble(zwe) 
        endif
        do i=1,fpm(23)
           res(i)=work_loc(i) !! save them all
           if ((lambda(i)>Emin).and.(lambda(i)<Emax))  res(i)=200d0+work_loc(i)
        enddo
        deallocate(work_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Remark: q=x and work=Ax-eBx is known already, using FEAST-MPI q and work are already distributed

        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration)

        loop=loop+1

        if (fpm(10)==1) fpm(33)=0 ! initialize factorization id 
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef MPI       
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif

!!! reorder (shift) lambda, eigenvector and residual
        call DLACPY( 'F', N, fpm(23),q , N, work, N ) 
        allocate(WORK_LOC(fpm(23)))
        allocate(WORK_LOC2(fpm(23)))
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
        deallocate(work_loc)
        deallocate(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        M0=fpm(23)  ! update value of M0 (new subspace)
        fpm(21)=100 ! The End

     end if ! test convergence
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==100) then !! THE END (ijob=0) 
     ijob=0 !! exit FEAST

if (((info==1).or.(info==3)).and.(fpm(40)/=0)) info=7 !search failed
     
     if (com) then !! Print  Information

        if (info>=200) then
           write(fout,'(A)',advance='no') 'PROBLEM with input parameters'
           write(fout,*) 
        end if


        if (info==-3) then
           write(fout,'(A)',advance='no') 'ERROR with reduced system'  
           write(fout,*) 
        end if

        !        if (info==-2) then
        !           write(fout,'(A)',advance='no') 'ERROR from Inner Linear System Solver in FEAST driver'  
        !           write(fout,*) 
        !        end if

        !        if (info==-1) then
        !           write(fout,'(A)',advance='no') 'ERROR with Internal memory allocation'
        !           write(fout,*) 
        !        end if

        if (info==1) then
           write(fout,'(A)',advance='no') '==>WARNING: No eigenvalue has been found in the proposed search interval'
           write(fout,*)
        endif

        if (info==2) then
           write(fout,'(A)',advance='no') '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed by fpm(4))'  
           write(fout,*)
        end if

        if (info==3) then
           write(fout,'(A)',advance='no') '==>WARNING: Size subspace M0 too small'
           write(fout,*)
        end if

        if (info==4) then
           write(fout,'(A)',advance='no') '==>WARNING: Only the subspace has been returned'  
           write(fout,*)
        end if

        if (info==5) then
           write(fout,'(A)',advance='no') '==>WARNING: Only stochastic estimation of #eigenvalues returned'  
           write(fout,*)
        end if

        if (info==7) then
           write(fout,'(A)',advance='no') '==>WARNING: Unfortunately, the search for extremum eigenvalues has failed- you need to specify you own interval and set fpm(40)=0'  
           write(fout,*)
        end if


        if (info==0) then
            write(fout,*)
            write(fout,'(A)',advance='no') '==>FEAST has successfully converged'
           if (fpm(6)==0) then  
              write(fout,'(A)',advance='no') ' with Trace tolerance <1E-'
           else
               write(fout,'(A)',advance='no') ' with Residual tolerance <1E-'
            end if
             if (fpm(3)<10) then
               write(fout,'(I1)') fpm(3)
            else
               write(fout,'(I2)') fpm(3)
            end if
            write(fout,'(A)',advance='no') '   # FEAST outside it. '
             write(fout,'(I8)') loop
             if (mod(fpm(30),10000)/1000==2) then !! ifeast
                write(fout,'(A)',advance='no') '   # Inner BiCGstab it.'
                write(fout,'(I8)',advance='no') fpm(60)
                if (fpm(30)/100000==2)    write(fout,'(A)',advance='no') ' for rank2=0'  !pfeast
                 write(fout,*)
              end if
           write(fout,'(A)',advance='no') '   # Eigenvalue found  '
           write(fout,'(I8)',advance='no') mode
           write(fout,'(A)',advance='no') ' from'
            write(fout,'(ES25.16)',advance='no') lambda(1)
            write(fout,'(A)',advance='no') ' to'
            write(fout,'(ES25.16)') lambda(mode)
        else
           write(fout,'(A)',advance='no') '==>INFO code = '
           write(fout,'(I4)',advance='no') info
           write(fout,*)
        end if
        
        if (info>=0) then
           call system_clock(tend,tim) !! total time
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           write(fout,*)
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST-RCI timing  |'
           write(fout,'(A)') '--------------------.------------------.'
           write(fout,'(A)',advance='no') '| Fact. cases(10,20)|'
           write(fout,'(F12.4)',advance='no') (t10*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Solve cases(11,12)|'
           write(fout,'(F12.4)',advance='no') (t11*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| A*x   cases(30,31)|'
           write(fout,'(F12.4)',advance='no') (t30*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| B*x   cases(40,41)|'
           write(fout,'(F12.4)',advance='no') (t40*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Misc. time        |'
           write(fout,'(F12.4)',advance='no') (tend-tini-t10-t11-t30-t40)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Total time (s)    |'
           write(fout,'(F12.4)',advance='no') (tend-tini)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)') '--------------------------------------- '

           write(fout,*)
        endif
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,'(A)',advance='no') '*********** FEAST- END*************************'
        write(fout,*) 
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,*)
        !! close file
        if (fpm(1)<0) close(fout)
     endif

#ifdef MPI
     call MPI_BARRIER(L2_COMM_WORLD,code)
#endif

  end if



end subroutine dfeast_srcix






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfeast_srci(ijob,N,Ze,work,workc,Aq,Bq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
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
  !  Aq,Bq      (input/output) REAL DOUBLE PRECISION (M0,M0)  :  Workspace for Reduced Eigenvalue System
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
  ! Eric Polizzi 2009-2019
  ! ====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  double precision, dimension(N,*) ::work
  complex(kind=(kind(1.0d0))),dimension(N,*):: workc
  integer,dimension(*) :: fpm
  double precision,dimension(M0,*):: Aq,Bq
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


  if (ijob==-1) then
     fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     if (fpm(30)==-111) then
        fpm(30)=121110 ! code name for direct call to this routine
        call feastdefault(fpm,info)
end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if
  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne)
  call dfeast_srcix(ijob,N,Ze,work,workc,Aq,Bq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)

end subroutine dfeast_srci





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)
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
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  ! Eric Polizzi 2009-2019
  ! ====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  !include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZZERO=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0,ab=pi/2.0d0
  integer(8),save :: fout
  logical,save :: com
  !! variable for FEAST
  character(len=3) :: ctemp 
  character(len=25) :: name
  integer :: i,e,j,k,Ntotal
  integer,dimension(4) :: iseed
  double precision :: theta,r
  complex(kind=(kind(1.0d0))) :: zxe,zwe,aux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: zBqo
  integer, dimension(:),allocatable :: fpm_default
  logical :: testconv,fact
  double precision :: trace
  !! Lapack variable (reduced system)
  character(len=1) :: JOBZ,UPLO
  double precision, dimension(:),allocatable :: work_loc,work_loc2
  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: zwork_loc
  integer :: lwork_loc,info_lap,infoloc
  !! MPI compatibility variables
  integer :: rank2,rank3,code,nb_procs2,nb_procs3,L2_COMM_WORLD,L3_COMM_WORLD
  double precision :: DZNRM2
  integer,save :: tini,t1,t10,t11,t30,t40,t50,t60
  integer :: tend,tim,t2
  character(len=8) :: date
  !  double precision, dimension(:,:),allocatable :: dwork
  ! complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: zzwork

!!! record timing
  if (ijob/=-1) then  
     call system_clock(t2,tim)
     select case(ijob)
     case(10,20)
        t10=t10+(t2-t1)
     case(11,21)
        t11=t11+(t2-t1)
     case(30)
        t30=t30+(t2-t1)
     case(40)
        t40=t40+(t2-t1)
     case(50)
        t50=t50+(t2-t1)
     case(60)
        t60=t60+(t2-t1)
     end select
  end if

!!! parallel init
  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1
  !----------------------------------------------
#ifdef MPI
  L2_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm(49)
  if (fpm(49)/=0) then
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  end if
#endif
  !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then

     info=0 ! default value  
     if (fpm(30)==-111) then
        fpm(30)=141211 ! code name for direct all to this routine   
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     endif


     call system_clock(tini,tim)
     t10=0
     t11=0
     t30=0
     t40=0
     t50=0
     t60=0

     ! set up factorization id
     fpm(33)=0

     ! comments?
     com=.false. !default
     if (fpm(1)/=0) com=.true.
#ifdef MPI
     if (com.and.((rank2/=0).or.(rank3/=0))) com=.false. ! comment only in rank 0 (only a single mpi is allowed to comment)
#endif
     !file or screen output?
     if (fpm(1)<0) then
        write(ctemp,'(I3)') abs(fpm(1))
        fout=abs(fpm(1))+200!fpm(60) ! file number default
        if (com) then
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append') !! all mpi procs may have access to the file (if needed)
           !write(fout,'(A)',advance='no') '############## NEW FEAST RUN ###################'
           !write(fout,*)
        endif
     elseif (fpm(1)>0) then
        fout=6 !screen
     endif



!!!!!!!!!!!!!!!! Find Ntotal
     Ntotal=N
#ifdef MPI           
     if (fpm(49)/=0) then
        Ntotal=0
        call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     end if
#endif      
!!!!!!!!!!!!!!!!! 
     fpm(26)=Ntotal


     call dcheck_feast_srci_input(Emin,Emax,M0,Ntotal,info)


     if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF (info==0) then


        if (com) then
           write(fout,*)
           write(fout,'(A)',advance='no') '***********************************************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '*********** FEAST v'
           write(fout,'(F3.1)',advance='no') fpm(31)/10.0d0
           write(fout,'(A)',advance='no') ' BEGIN ******************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '***********************************************'  
           write(fout,*)
            if (fpm(1)<0) then ! only in-file
           call date_and_time(DATE=date)
           write(fout,'(A)',advance='no') 'Date: '//trim(date(1:4))//'-'//trim(date(5:6))//'-'//trim(date(7:8))
           write(fout,*)
        end if
           call feast_name(fpm(30),name)
           write(fout,'(A)') 'Routine '//trim(name)
if ((mod(fpm(30),10)==2).or. (mod(fpm(30),10)==3)) then
                write(fout,'(A)',advance='no') 'Solving AX=eX with A complex Hermitian'
           elseif    ((mod(fpm(30),10)==4).or. (mod(fpm(30),10)==5)) then 
              write(fout,'(A)',advance='no') 'Solving AX=eBX with A complex Hermitian and B hpd'
           end if
           write(fout,*)
           if (nb_procs2*nb_procs3>1) then
              write(fout,'(A)',advance='no') '#MPI (total=L2*L3) '
              write(fout,'(I4)',advance='no') nb_procs2*nb_procs3
              write(fout,'(A)',advance='no') '='
              write(fout,'(I4)',advance='no') nb_procs2
              write(fout,'(A)',advance='no') '*'
              write(fout,'(I4)',advance='no') nb_procs3
              write(fout,*)
           endif
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
           write(fout,'(A)',advance='no') 'List of input parameters fpm(1:64)-- if different from default'
           write(fout,*)
           allocate(fpm_default(64))
           call feastinit(fpm_default) ! initialize to -111
           fpm_default(1)=0
           fpm_default(30)=fpm(30) ! name
           fpm_default(14)=fpm(14) ! for special conditional cases
           call feastdefault(fpm_default,infoloc) ! get the default values
           do i=1,64   ! compare with possible mismatch with user values
              if ((i<=19).or.((i>=40).and.(i<=50))) then
                 if ((fpm(i)/=fpm_default(i)).and.(i/=9).and.(i/=49)) then
                    write(fout,'(A)',advance='no') '   fpm('
                    write(fout,'(I2)',advance='no') i
                    write(fout,'(A)',advance='no') ')='
                    write(fout,'(I4)',advance='no') fpm(i)
                    write(fout,*)
                 endif
              endif
           enddo
           deallocate(fpm_default)
           fpm(22)=nb_procs3 ! temp copy
           call dfeast_info(fout,fpm,Emin,Emax,Ntotal,M0)

        end if


!!!!!!!!!!!!!!!!

        if (com) then
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST runs        |'
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           if (fpm(14)/=2) then
              write(fout,'(A)',advance='no') '#It |  #Eig  |          Trace            |     Error-Trace          |     Max-Residual'  
           else
              write(fout,'(A)',advance='no') 'Running average for stochastic estimation (1 -> M0)'  
           endif
           write(fout,*)
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
        endif


!!!!!!!!!!! Some initialization
        fpm(22)=fpm(2) ! only one half contour necessary
        loop=0


        fpm(23)=min(M0,fpm(26)) ! 'current M0' size (global value)
        if (fpm(14)==2) fpm(23)=min(fpm(32),fpm(26),M0) ! stochastic estimate
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs2>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(24)=1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors (option 2 in ZLARNV)
           iseed=(/56,890,3456,2333/)
           iseed(1) = iseed(1) * (1+rank3*13)
           iseed(2) = iseed(2) + rank3*7
           iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
           iseed(4) = iseed(4) + rank3*2

           ! copy work into q for multiplication by B matrix, if stochastic estimate  is "on"
           if (fpm(14)==2) then 
              call ZLARNV(2,iseed,N*fpm(23),q(1,1))
           else
              call ZLARNV(2,iseed,N*fpm(23),work(1,1)) ! random distribution
              !             allocate(dwork(N,M0))
              !call DLARNV(2,iseed,N*fpm(23),dwork(1,1))
              !work(1:N,1:M0)=dwork(1:N,1:M0)*ZONE
              !deallocate(dwork)
              !  allocate(zzwork(Ntotal,M0))
              ! call ZLARNV(4,iseed,Ntotal*fpm(23),zzwork(1,1))
              ! if (rank3==0) work(1:N,1:M0)=zzwork(1:N,1:M0)
              ! if (rank3==1) work(1:N,1:M0)=zzwork(N+1:Ntotal,1:M0)
              ! deallocate(zzwork)




           endif
        endif

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=ZZERO 
#endif
           !-----------------------------------------
           if (fpm(5)==1) then ! Compute residual vector R=AX-BXE 
              ijob=30 !!A*q==>work
              fpm(21)=10
           else ! only stochastic estimate
              ijob=40 !! B*q=>work
           endif
           call system_clock(t1,tim)
           return
        end if
     end IF    ! info=0
!!!!!!!!!!!!!!
  end if  !ijob=-1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  IF (fpm(21)==1) then !! we initialize a new contour integration

     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(14)==2)) then !! condition on fpm(5) unnecessary !!<< 
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     endif

     if ((loop>0).and.(fpm(23)/nb_procs2>=1)) then
        if (.not.((loop==1).and.(fpm(5)==1))) then
           call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        end if
     end if
#endif
     !------------------------------------------------------------------------

!!!!!!!!!!! RESIDUAL INVERSE ITERATION
     if (loop==0) then
        ! if ((fpm(29)==1).and.(fpm(16)==2).and.(fpm(5)==1)) then !! Zolo and fpm(5)=1
        !    call dset_zolotarev(fpm(2),0,zxe,zwe) ! get zwe
        !    if (rank2==0) then
        !       q(1:N,1:fpm(23))=-dble(zwe)*q(1:N,1:fpm(23))
        !    else 
        !       q(1:N,1:fpm(23))=ZZERO
        !    endif
        ! else 
        q(1:N,1:fpm(23))=ZZERO
        ! end if
     else if (loop>0) then ! general case, works for Zolo as well..res(i) includes the offset
        do i=1,fpm(23)
           if (res(i)>100d0) then
              call ZSCAL(N,-ZONE*(res(i)-200d0), q(1,i), 1)
           else
              call ZSCAL(N,-ZONE*res(i), q(1,i), 1)
           end if
        end do
        if (rank2/=0)  q(1:N,1:fpm(23))=ZZERO ! for parallelism
     end if
!!!!!!!!!!!!!!

     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     ! IF (info==0) then !! will end up checking info errors returned by FEAST drivers
     do e=fpm(20)+rank2,fpm(22),nb_procs2 !!!! loop over the contour 

        if (ijob==-2) then !!Factorize the linear system (complex) (zS-A)
           ! set up factorization id (if fact. storage is used)
           if (fpm(10)==1) then
              fpm(33)=fpm(33)+1
           else
              fpm(33)=1
           endif

           Ze=Zne(e)
           fpm(20)=e-rank2

           ijob=10 ! for fact

           fact=.false. ! action to factorize
           if ((loop==0).or.((loop==1).and.(fpm(5)==1))) then
              fact=.true. ! always factorize in initial loop. Remark: with input initial guess, loop=1 is the initial loop
           else !(not initial loop)
              if ((fpm(22)>nb_procs2).and.(fpm(10)==0)) fact=.true. ! re-factorize if #mpi<#contour points and fact not stored
           end if
           if (fact) then
              call system_clock(t1,tim)
              return 
           end if

!!$           
!!$           if ((loop==0).or.(fpm(22)>nb_procs2)) then !no refactorization if one linear system per processor
!!$              if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) then !no refactorization if stored in memory
!!$                 if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) then 
!!$                    call system_clock(t1,tim)
!!$                    return 
!!$                 endif
!!$              end if
!!$           endif


        end if
!!!!!!!!!!!!!!!!           
        if (ijob==10) then !!Solve the linear system (complex) (zS-A)q=v
!!!!!!!!!!!!!!!!!

           call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
           ijob=11 ! for solve
           call system_clock(t1,tim)
           return
        endif

!!!!!!!!!!!!!!!!!!           
        if (ijob==11) then
!!!!!!!!!!!!!!!!!!              
           !! summation 
           aux=-Wne(e)
           if (loop>0) then ! for residual inverse iteration
              do i=1,fpm(23)
                 call ZSCAL(N,ZONE/(ze-lambda(i)*ZONE), workc(1,i), 1)
              enddo
           endif

           call ZAXPY(N*fpm(23),aux,workc,1,q,1)
           !!Explicit Factorization of the linear system (complex) (zS-A)^T 
           !!needed if driver not capable to exploit Factorization of (zS-A) for solving (zS-A)^Tq=v          
           ijob=20 ! for fact

           fact=.false. ! action to factorize
           if ((loop==0).or.((loop==1).and.(fpm(5)==1))) then
              fact=.true. ! always factorize in initial loop. Remark: with input initial guess, loop=1 is the initial loop
           else !(not initial loop)
              if ((fpm(22)>nb_procs2).and.(fpm(10)==0)) fact=.true. ! re-factorize if #mpi<#contour points and fact not stored
           end if
           if (fact) then
              call system_clock(t1,tim)
              return 
           end if


!!$           if ((loop==0).or.(fpm(22)>nb_procs2)) then !no refactorization if one linear system per processor
!!$              if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) then !no refactorization if stored in memory
!!$                 if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) then 
!!$                    call system_clock(t1,tim)
!!$                    return 
!!$                 endif
!!$              end if
!!$           end if
!!$           
        end if


!!!!!!!!!!!!!!!!           
        if (ijob==20) then!!!! Solve the linear system (complex) (zS-A)^Tq=v
!!!!!!!!!!!!!!!!              
           call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
           ijob=21 ! for solve with transpose
           call system_clock(t1,tim)
           return
        end if

!!!!!!!!!!!!!!!!!
        if (ijob==21) then
!!!!!!!!!!!!!!!!!!!!!!!              
           aux=-conjg(Wne(e))
           if (loop>0) then ! for residual inverse iteration
              do i=1,fpm(23)
                 call ZSCAL(N,ZONE/conjg(ze-lambda(i)*ZONE), workc(1,i), 1)
              enddo
           endif
           call ZAXPY(N*fpm(23),aux,workc,1,q,1)

           ijob=-2 ! just for identification
        end if

     end do
     !     end IF !! info=0

     !------------------------------------------------
     !#ifdef MPI
     !     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
     !#endif
     !-----------------------------------------------  
     !    if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(21)=4 
        !------------------------------------------------
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if
     !print *,rank2,q(1:N,1:M0)  
  end IF   ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !call MPI_BARRIER(L2_COMM_WORLD)
  !stop
  !q(1:N,1:M0)=dble(q(1:N,1:M0))*ZONE

  !print *,sum(q(1:N,1))
  !stop


  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only q vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic Estimate
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only q is returned with stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     iseed(1) = iseed(1) * (1+rank3*13)
     iseed(2) = iseed(2) + rank3*7
     iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
     iseed(4) = iseed(4) + rank3*2

     call ZLARNV(2,iseed,N*fpm(23),work(1,1)) 
     call ZGEMM('C','N',fpm(23),fpm(23),N,-ZONE,work,N,q,N,ZZERO,zAq,M0) ! projection
     call ZGEMM('C','N',fpm(23),fpm(23),N,ZONE,work,N,work,N,ZZERO,zBq,M0) ! normalization
     !----------------------------------------  
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !---------------------------------------

     Ntotal=fpm(26)
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+dble(zAq(i,i)/zBq(i,i))
        res(i)=abs(theta*(Ntotal/(DONE*i)))
        if (com) then
           write(fout,'(I3)',advance='no') i
           write(fout,'(3X)',advance='no')
           write(fout,'(F10.2)',advance='no') res(i)
           write(fout,*)
        endif
     enddo
     mode=int(real(res(fpm(23))))
     if (com) then
        write(fout,'(A)',advance='no') '# estimated '
        write(fout,'(3X)',advance='no')
        write(fout,'(I6)',advance='no') mode
        write(fout,*)
     endif
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Bq xq     with Aq=Q^TAQ Bq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!form Bq=> Bq=Q^T B Q 
  if (fpm(21)==4) then
#ifdef MPI 
     work(1:N,1:fpm(23))=ZZERO
#endif
     fpm(21)=5 ! preparing reentry
     ijob=40
     call system_clock(t1,tim)
     return  ! mat-vec B*q => work
  endif

  if (fpm(21)==5) then 
     !------------------------------------------------
#ifdef MPI 
     zBq(1:M0,1:fpm(23))=ZZERO
#endif
     !-------------------------------------------------
     fpm(21)=51


     if (fpm(13)>=2) then ! customize inner product
        ijob=50
        call system_clock(t1,tim)
        return
     endif

     call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,q(1,1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0)


  endif


!!!!!!!!!!!!!!  
  if (fpm(21)==51) then
     if (fpm(13)>=2) call ZCOPY(M0*fpm(25),zAq(1,fpm(24)),1,zBq(1,fpm(24)),1)

     !---------------------------------------- !(zBq known to all processors) 
#ifdef MPI
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     end if
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !---------------------------------------
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==51) then
     mode=0
     if (loop>0) then
        do i=fpm(24),fpm(24)+fpm(25)-1
           if (res(i)>100d0) then ! indicator for being inside the search interval
              if (abs(sqrt(abs(zBq(i,i)))-(res(i)-200.0d0))/sqrt(abs(zBq(i,i)))<0.1d0) mode=mode+1
           endif
        enddo
     endif
!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,L2_COMM_WORLD,code)
#endif
!!!!!!!!!!!!!!!!!!
     fpm(21)=6
  endif


!!!!!!!!!form  Aq=> Aq=Q^T A Q
  if (fpm(21)==6) then
     !------------------------------------------------
#ifdef MPI 
     work(1:N,1:fpm(23))=ZZERO
#endif 
     fpm(21)=7 ! preparing reenty
     ijob=30
     call system_clock(t1,tim)
     return! mat-vec A*q => work
  end if


  if (fpm(21)==7) then
     !------------------------------------------------
#ifdef MPI
     zAq(1:M0,1:fpm(23))=ZZERO
#endif
     !-------------------------------------------------
     fpm(21)=8
     if (fpm(13)>=2) then ! customize inner product
        ijob=50
        call system_clock(t1,tim)
        return
     endif


     call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,q(1,1),N,work(1,fpm(24)),N,ZZERO,zAq(1,fpm(24)),M0)
  endif


  if (fpm(21)==8) then
     !----------------------------------------!(zAq known to all processors) 
#ifdef MPI 
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     end if
     if(fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
#endif

     !---------------------------------------
     if  ((fpm(13)==1).or.(fpm(13)==3)) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=60
        call system_clock(t1,tim)
        return
     endif
     fpm(21)=81
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==81) then
     ! if using FEAST-MPI ==> solve on a single proc
     if ((rank2==0).and.(rank3==0)) then
        JOBZ='V'
        UPLO='L'
        info_lap=1 ! initialization
        i=1
        LWORK_LOC=2*fpm(23)-1 !! for lapack eig reduced system
        allocate(zWORK_LOC(LWORK_LOC))
        allocate(WORK_LOC(3*fpm(23)-2))

        do while ((info_lap/=0).and.(info==0))
           i=i+1
           if (i==10) info=-3 ! arbitrary maximum
           allocate(zBqo(fpm(23),fpm(23)))

           call ZLACPY( 'F', fpm(23), fpm(23),zBq , M0, zBqo, fpm(23) )
           call ZHEGV(1, JOBZ, UPLO, fpm(23), zAq, M0, zBqo, fpm(23), lambda, zWORK_loc,Lwork_loc, WORK_loc, INFO_lap)

           if ((info_lap<=fpm(23)).and.(info_lap/=0)) info=-3
           if (info_lap>fpm(23)) then !! zBqo is not spd (a posteriori resize subspace)
              fpm(23)=info_lap-fpm(23)-1
              if (com) then
                 write(fout,'(A)',advance='no') 'Resize subspace'  
                 write(fout,'(3X)',advance='no') 
                 write(fout,'(I4)',advance='no') fpm(23)
                 write(fout,*)
              end if
           end if
           deallocate(zBqo)
        end do
        deallocate(zwork_loc)
        deallocate(work_loc)
     end if !(rank 0)
     !-------------------------------- !(info common to all processors)
#ifdef MPI
     if (rank3==0) call MPI_BCAST(info,1,MPI_INTEGER,0,L2_COMM_WORLD,code)
     if (fpm(49)/=0) call MPI_BCAST(info,1,MPI_INTEGER,0,L3_COMM_WORLD,code)
#endif     
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end

     if (info==0) then
        fpm(25)=fpm(23) !! current M0 size (by default) -- global
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 



#ifdef MPI
        if (rank3==0) call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,L2_COMM_WORLD,code)
        if (fpm(49)/=0) call MPI_BCAST(fpm(23),1,MPI_INTEGER,0,L3_COMM_WORLD,code)
        if (rank3==0) then
           call MPI_BCAST(zAq,M0*fpm(23),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,L2_COMM_WORLD,code)
        end if
        if (fpm(49)/=0) then
           call MPI_BCAST(zAq,fpm(23)*M0,MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_PRECISION,0,L3_COMM_WORLD,code)
        end if

        !  if (rank2==0) print *,'Aq',zAq(1:N,1:M0)  

        if (fpm(23)/nb_procs2>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs2 ! local size of current M0
           if (rank2==nb_procs2-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(24)=1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------
        fpm(21)=9
     end if
  end if !! fpm(21)=81


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors X=Qxq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==9) then
     !! from precious step work=A*Q, before using work memory space
     !! let us construct first workc=work*Aq which is workc=A*X
     !------------------------------------------------------------------------
#ifdef MPI
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !------------------------------------------------------------------------

     call ZGEMM('N','N',N,fpm(25),fpm(23),ZONE,work(1,1),N,zAq(1,fpm(24)),M0,ZZERO,workc(1,fpm(24)),N)

     !! construct X (q needed for outer-loop as well)
     call ZLACPY( 'F', N, fpm(23),q(1,1) , N, work(1,1), N )

#ifdef MPI 
     q(1:N,1:fpm(23))=ZZERO
#endif

     call ZGEMM('N','N',N,fpm(25),fpm(23),ZONE,work(1,1),N,zAq(1,fpm(24)),M0,ZZERO,q(1,fpm(24)),N)

     fpm(21)=10 !! no need to compute A*q again
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residual ||AX-lambda*BX||
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (fpm(21)==10) then
     !handle special case of input initial guess
     !     if ((fpm(5)==1).and.(loop==0))  call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
     if ((fpm(5)==1).and.(loop==0))  call ZLACPY( 'F', N, fpm(25),work(1,fpm(24)) , N, workc(1,fpm(24)), N )


     fpm(21)=11 ! preparing reentry
     ijob=40 
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZZERO  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     call system_clock(t1,tim)
     return  ! mat-vec B*q => work 
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !---------Compute residual vector R=AX-BXE
     !         work<=workc-work*E (needed also for outer-loop)
     do i=fpm(24),fpm(24)+fpm(25)-1 
        !res(i)=sum(abs(work(1:N,i))) ! placeholder (denominator)
        res(i)=DZNRM2(N,work(1,i),1)**2 ! square because of mpi compatibility sqrt(mpi1)+sqrt(mpi2) /= sqrt(mpi1+mpi2) 
        call ZSCAL(N,-ZONE*lambda(i),work(1,i),1)
     enddo
     call ZAXPY(N*fpm(25),ZONE,workc(1,fpm(24)),1,work(1,fpm(24)),1)

     !-------- Residual

     !reduce first denominator
#ifdef MPI 
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif

     do i=fpm(24),fpm(24)+fpm(25)-1
        !res(i)=sum(abs(work(1:N,i)))/(max(abs(Emin),abs(Emax))*res(i))
        res(i)=(DZNRM2(N,work(1,i),1)**2)/((max(abs(Emin),abs(Emax))**2)*res(i))
     end do
     !----------------------------------------
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !-----------------------------------------
     do i=1,fpm(23)
        res(i)=sqrt(res(i)) ! norm L2
     enddo


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


     if (loop>2) then! wait third iteration (spurious and ifeast related)
        if (mode==0) info=1  ! no eigenvalue detected in the interval
        if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
     endif
     if (info/=0) fpm(21)=100 ! The End
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Check FEAST Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

  if (fpm(21)==11) then
     testconv=.false. ! initialization

     !! trace convergence
     if (loop==0) epsout=DONE
     if (loop>0) then
        epsout=(abs(trace-epsout))/max(abs(Emin),abs(Emax))           
        !if ((fpm(6)==0).and.(log10(epsout)<(-fpm(3)))) testconv=.true.  ! epsout coule be 0
        if ((fpm(6)==0).and.(epsout<10.0d0**(-fpm(3))))  testconv=.true.
     end if

     !! residual convergence
     !if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.
      if ((fpm(6)/=0).and.(theta<10.0d0**(-fpm(3))))  testconv=.true.
    

     !! Two loops minimum if spurious are found
     if ((loop<=1).and.(k>mode)) testconv=.false.


!!!! L2 barrier to guarantee synchronization in printing
#ifdef MPI
 call MPI_BARRIER(L2_COMM_WORLD,code)
#endif 
     

     if (com) then
        if (fpm(5)==0) then
           write(fout,'(I3)',advance='no') loop
        else ! no integration in first loop
           if (loop==0) then
              write(fout,'(A)',advance='no') 'N/A'
           else
              write(fout,'(I3)',advance='no') loop-1
           end if
        endif
        write(fout,'(3X)',advance='no') 
        write(fout,'(I4)',advance='no') mode
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') trace
        write(fout,'(3X)',advance='no') 
        write(fout,'(ES25.16)',advance='no') epsout
        write(fout,'(3X)',advance='no') 
        write(fout,'(ES25.16)',advance='no') theta
        write(fout,*) 
     end if

     if (.not.testconv) then
        epsout=trace
        e=loop
        if (fpm(5)==1) e=e-1
        if (e==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif

     !! handle a special case (prevent convergence if no eigenvalue found)
     !! for mode=0 FEAST will exit after few iterations (see above)
     if (mode==0) testconv=.false. 

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
        allocate(WORK_LOC(fpm(23)))
        call dfeast_rationalx(Zne,Wne,fpm(2),lambda,fpm(23),work_loc)
        if ((fpm(29)==1).and.(fpm(16)==2)) then ! zolotarev case
           call dset_zolotarev(fpm(2),0,zxe,zwe)
           work_loc(1:fpm(23))=work_loc(1:fpm(23))+dble(zwe)
        endif
        do i=1,fpm(23)
           res(i)=work_loc(i)
           if ((lambda(i)>Emin).and.(lambda(i)<Emax))  res(i)=200d0+work_loc(i)
        enddo
        deallocate(work_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Remark: q=x and work=Ax-eBx is known already, using FEAST-MPI q and work are already distributed 

        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        if (fpm(10)==1) fpm(33)=0 ! initialize factorization id 
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    !!!!!!! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef MPI       
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
#endif

!!! reorder (shift) lambda, eigenvector and residual
        call ZLACPY( 'F', N, fpm(23),q , N, work, N )
        allocate(WORK_LOC(fpm(23)))
        allocate(WORK_LOC2(fpm(23)))
        call DCOPY(fpm(23),lambda , 1, work_loc, 1 )
        call DCOPY(fpm(23),res , 1, work_loc2, 1 )
        q(1:N,1:fpm(23))=ZZERO
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
        deallocate(work_loc)
        deallocate(work_loc2)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        M0=fpm(23)  ! update value of M0 (new subspace)
        fpm(21)=100 ! The End

     end if ! test convergence
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==100) then !! THE END (ijob=0)
     ijob=0 !! exit FEAST


if (((info==1).or.(info==3)).and.(fpm(40)/=0)) info=7 !search failed
     
     if (com) then !! Print  Information

        if (info>=200) then
           write(fout,'(A)',advance='no') 'PROBLEM with input parameters'
           write(fout,*) 
        end if

        if (info==-3) then
           write(fout,'(A)',advance='no') 'ERROR with reduced system'
           write(fout,*) 
        end if

        !        if (info==-2) then
        !           write(fout,'(A)',advance='no') 'ERROR from Inner Linear System Solver in FEAST driver'  
        !           write(fout,*) 
        !        end if

        !        if (info==-1) then
        !           write(fout,'(A)',advance='no') 'ERROR with Internal memory allocation'
        !           write(fout,*) 
        !        end if

        if (info==1) then
           write(fout,'(A)',advance='no') '==>WARNING: No eigenvalue has been found in the proposed search interval'
           write(fout,*)
        endif
        
 if (info==2) then
           write(fout,'(A)',advance='no') '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed by fpm(4))'  
           write(fout,*)
        end if

        
        if (info==3) then
           write(fout,'(A)',advance='no') '==>WARNING: Size subspace M0 too small'  
           write(fout,*)
        end if

        if (info==4) then
           write(fout,'(A)',advance='no') '==>WARNING: Only the subspace has been returned'  
           write(fout,*)
        end if

        if (info==5) then
           write(fout,'(A)',advance='no') '==>WARNING: Only stochastic estimation of #eigenvalues returned'  
           write(fout,*)
        end if

      if (info==7) then
           write(fout,'(A)',advance='no') '==>WARNING: Unfortunately, the search for extremum eigenvalues has failed- you need to specify you own interval and set fpm(40)=0'  
           write(fout,*)
        end if
   
        if (info==0) then
            write(fout,*)
            write(fout,'(A)',advance='no') '==>FEAST has successfully converged'
           if (fpm(6)==0) then  
              write(fout,'(A)',advance='no') ' with Trace tolerance <1E-'
           else
               write(fout,'(A)',advance='no') ' with Residual tolerance <1E-'
            end if
             if (fpm(3)<10) then
               write(fout,'(I1)') fpm(3)
            else
               write(fout,'(I2)') fpm(3)
            end if
              write(fout,'(A)',advance='no') '   # FEAST outside it. '
             write(fout,'(I8)') loop
             if (mod(fpm(30),10000)/1000==2) then !! ifeast
                 write(fout,'(A)',advance='no') '   # Inner BiCGstab it.'
                write(fout,'(I8)',advance='no') fpm(60)
                if (fpm(30)/100000==2)    write(fout,'(A)',advance='no') ' for rank2=0'  !pfeast
                 write(fout,*)
              end if
              write(fout,'(A)',advance='no') '   # Eigenvalue found  '
                write(fout,'(I8)',advance='no') mode
           write(fout,'(A)',advance='no') ' from'
            write(fout,'(ES25.16)',advance='no') lambda(1)
            write(fout,'(A)',advance='no') ' to'
            write(fout,'(ES25.16)') lambda(mode)
        else
           write(fout,'(A)',advance='no') '==>INFO code = '
           write(fout,'(I4)',advance='no') info
           write(fout,*)
        end if
        
       
        if (info>=0) then
           call system_clock(tend,tim) !! total time
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           write(fout,*)
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST-RCI timing  |'
           write(fout,'(A)') '--------------------.------------------.'
           write(fout,'(A)',advance='no') '| Fact. cases(10,20)|'
           write(fout,'(F12.4)',advance='no') (t10*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Solve cases(11,12)|'
           write(fout,'(F12.4)',advance='no') (t11*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| A*x   cases(30,31)|'
           write(fout,'(F12.4)',advance='no') (t30*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| B*x   cases(40,41)|'
           write(fout,'(F12.4)',advance='no') (t40*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Misc. time        |'
           write(fout,'(F12.4)',advance='no') (tend-tini-t10-t11-t30-t40)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Total time (s)    |'
           write(fout,'(F12.4)',advance='no') (tend-tini)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)') '--------------------------------------- '

           write(fout,*)
        endif
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,'(A)',advance='no') '*********** FEAST- END*************************'
        write(fout,*) 
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,*)
        !! close file
        if (fpm(1)<0) close(fout)
     endif
#ifdef MPI
     call MPI_BARRIER(L2_COMM_WORLD,code)
#endif

  end if

end subroutine zfeast_hrcix






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_hrci(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info)
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
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  ! Eric Polizzi 2009-2019
  ! ====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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


  if (ijob==-1) then
     fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     if (fpm(30)==-111) then
        fpm(30)=141210 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if

  call zfeast_contour(Emin,Emax,fpm(2),fpm(16),fpm(18),Zne,Wne) 
  call zfeast_hrcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne)

end subroutine zfeast_hrci






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_grcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      Aqr=lambda Bqr and A^Hql=conjg(lambda) B^Hql                            
  ! 
  !  Expert comments: This is also a Kernel routine (contains several cases depending of fpm(15))           
  !
  !                   fpm(15)=0 ! 2 sided algo (default)
  !                   fpm(15)=1 ! 1 sided algo- right eigenvector
  !                   fpm(15)=2 ! 1 sided algo with right/left conjugate- default for complex symmetric
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
  !                            expert comment: size of (N,M0) if fpm(15)>0  
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !                            expert comment: size of (N,M0) if fpm(15)>0 (only right vectors)
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right (1:mode) and left (M0+1:M0+mode)
  !                          
  !                            expert comment: size of (M0) if fpm(15)>0 (only rigth vectors)
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
  !  Eric Polizzi - 2019
  !=====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  !include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZZERO=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0,ab=pi/2.0d0
  integer*8,save :: fout
  logical, save :: com
  !! variable for FEAST
  character(len=3) :: ctemp 
  character(len=25) :: name
  integer :: i,e,j,k,Ntotal,jj,m_min
  integer,dimension(4) :: iseed
  double precision :: theta
  complex(kind=(kind(1.0d0))) :: jac,aux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: zBqo
  integer, dimension(:),allocatable :: fpm_default,color
  logical :: testconv,fact
  complex(kind=(kind(1.0d0))) :: trace, zalpha
  !! Lapack variable 
  double precision, dimension(:),allocatable :: work_locr,work_loc2
  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: zwork_loc!,tau
  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: LALPHA,LBETA
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: VL,VR
  integer :: lwork_loc,info_lap,infoloc
  character(len=1) :: JOBVL
  !! MPI compatibility variables
  integer :: rank2,rank3,code,nb_procs2,nb_procs3,L2_COMM_WORLD,L3_COMM_WORLD
  double precision :: DZNRM2
  integer,save :: tini,t1,t10,t11,t30,t40,t50,t60
  integer :: tend,tim,t2
  character(len=8) :: date
  ! double precision, dimension(:,:),allocatable :: dwork
  !  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: zzwork

!!!!record timing
  if (ijob/=-1) then  
     call system_clock(t2,tim)
     select case(ijob)
     case(10,20)
        t10=t10+(t2-t1)
     case(11,21)
        t11=t11+(t2-t1)
     case(30,31)
        t30=t30+(t2-t1)
     case(40,41)
        t40=t40+(t2-t1)
     case(50)
        t50=t50+(t2-t1)
     case(60)
        t60=t60+(t2-t1)
     end select
  end if

!!! parallel init
  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1
  !----------------------------------------------
#ifdef MPI
  L2_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm(49)
  if (fpm(49)/=0) then
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  end if
#endif
  !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then

     info=0 ! default value  
     if (fpm(30)==-111) then
        fpm(30)=141311 ! code name for direct all to this routine   
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     endif


     call system_clock(tini,tim)
     t10=0
     t11=0
     t30=0
     t40=0
     t50=0
     t60=0
     ! set up factorization id
     fpm(33)=0

     ! expert routines Emid,r needs to be calculated
     if(fpm(29)==0)then 
        ! calculate the center of mass and relative maximum length of custom contour
        Emid=ZZERO
        do i=1,fpm(8)
           Emid = Emid + Zne(i)
        enddo
        Emid = Emid / fpm(8)
        r = abs(maxval(abs(Zne(1:fpm(8))))-Emid) ! for normalization of residual
     end if



     ! comments?
     com=.false. !default
     if (fpm(1)/=0) com=.true.
#ifdef MPI
     if (com.and.((rank2/=0).or.(rank3/=0))) com=.false. ! comment only in rank 0 (only a single mpi is allowed to comment)
#endif
     !file or screen output?
     if (fpm(1)<0) then
        write(ctemp,'(I3)') abs(fpm(1))
        fout=abs(fpm(1))+200!fpm(60) ! file number default
        if (com) then
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append') !! all mpi procs may have access to the file (if needed)
           !write(fout,'(A)',advance='no') '############## NEW FEAST RUN ###################'
           !write(fout,*)
        endif
     elseif (fpm(1)>0) then
        fout=6 !screen
     endif

!!!!!!!!!!!!!!!! Find Ntotal
     Ntotal=N
#ifdef MPI           
     if (fpm(49)/=0) then
        Ntotal=0
        call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     end if
#endif        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
     fpm(26)=Ntotal



     call dcheck_feast_grci_input(r,M0,Ntotal,info)


     if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     IF (info==0) then


        if (com) then
           write(fout,*)
           write(fout,'(A)',advance='no') '***********************************************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '*********** FEAST v'
           write(fout,'(F3.1)',advance='no') fpm(31)/10.0d0
           write(fout,'(A)',advance='no') ' BEGIN ******************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '***********************************************'  
           write(fout,*)
        if (fpm(1)<0) then ! only in-file
           call date_and_time(DATE=date)
           write(fout,'(A)',advance='no') 'Date: '//trim(date(1:4))//'-'//trim(date(5:6))//'-'//trim(date(7:8))
           write(fout,*)
        end if
           call feast_name(fpm(30),name)
           write(fout,'(A)',advance='no') 'Routine '//trim(name)
           write(fout,*)
           if (nb_procs2*nb_procs3>1) then
              write(fout,'(A)',advance='no') '#MPI (total=L2*L3) '
              write(fout,'(I4)',advance='no') nb_procs2*nb_procs3
              write(fout,'(A)',advance='no') '='
              write(fout,'(I4)',advance='no') nb_procs2
              write(fout,'(A)',advance='no') '*'
              write(fout,'(I4)',advance='no') nb_procs3
              write(fout,*)
           endif
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
           write(fout,'(A)',advance='no') 'List of input parameters fpm(1:64)-- if different from default'
           write(fout,*)
           allocate(fpm_default(64))
           call feastinit(fpm_default) ! initialize to -111
           fpm_default(1)=0
           fpm_default(30)=fpm(30) ! name
           fpm_default(14)=fpm(14) ! for special conditional cases
           call feastdefault(fpm_default,infoloc) ! get the default values
           do i=1,64   ! compare with possible mismatch with user values
              if ((i<=19).or.((i>=40).and.(i<=50))) then
                 if ((fpm(i)/=fpm_default(i)).and.(i/=9).and.(i/=49)) then
                    write(fout,'(A)',advance='no') '   fpm('
                    write(fout,'(I2)',advance='no') i
                    write(fout,'(A)',advance='no') ')='
                    write(fout,'(I4)',advance='no') fpm(i)
                    write(fout,*)
                 endif
              endif
           enddo
           deallocate(fpm_default)
           fpm(22)=nb_procs3 ! temp copy
           call zfeast_info(fout,fpm,Emid,r,Ntotal,M0)



        end if



        if (com) then
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST runs        |'
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           if (fpm(14)/=2) then
              write(fout,'(A)',advance='no') '#It |  #Eig  |         |Trace|           |     Error-Trace          |     Max-Residual'    
           else
              write(fout,'(A)',advance='no') 'Running average for stochastic estimation (1 -> M0)'  
           endif
           write(fout,*)
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'

        end if



!!!!!!!!!!! Some initialization
        fpm(22)=fpm(8) ! full contour necessary
        loop=0


        fpm(23)=min(M0,fpm(26)) ! 'current M0' size (global value)
        if (fpm(14)==2) fpm(23)=min(fpm(32),fpm(26),M0) ! stochastic estimate
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        fpm(35)=fpm(23) !! 'current M0' size (by default)
        fpm(34)=M0+1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs2>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(24)=1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
           fpm(35)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(35)=fpm(35)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(34)=M0+1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors
           iseed=(/56,890,3456,2333/)
           iseed(1) = iseed(1) * (1+rank3*13)
           iseed(2) = iseed(2) + rank3*7
           iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
           iseed(4) = iseed(4) + rank3*2


           if (fpm(14)==2) then
              ! copy work into q for multiplication by B matrix, if stochastic estimate  is "on" 
              call ZLARNV(2,iseed,N*fpm(23),q(1,1)) !!<<4?
           else
              call ZLARNV(2,iseed,N*fpm(23),work(1,1)) !!<<4??

              if (fpm(15)==0) then !! 2 contours for left eigenvectors
                 iseed=(/6,80,456,5421/)
                 iseed(1) = iseed(1) * (1+rank3*13)
                 iseed(2) = iseed(2) + rank3*7
                 iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
                 iseed(4) = iseed(4) + rank3*2
                 call ZLARNV(2,iseed,N*fpm(23),work(1,M0+1))

              end if
           end if
        end if

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=ZZERO 
           if (fpm(15)==0) work(1:N,M0+1:M0+fpm(23))=ZZERO ! 2 contours for left eigenvectors
#endif
           !------------------------------------------
           if (fpm(5)==1) then ! Compute residual vector R=AX-BXE 
              ijob=30 !!A*q==>work
              fpm(21)=10
           else ! only stochastic estimate!! B*Q(1,1:fpm(23))=>work(1:N,1:fpm(23))
              ijob=40 !! B*q=>work
           endif
           call system_clock(t1,tim)
           !         if (fpm(15)==0) fpm(21)=-1  
           return
        end if

     end IF  ! info=0
!!!!!!!!!!!!!!
  end if    !ijob=-1

  !  if (fpm(21)==-1) then !! we need to initialize the left eigenvectors as well
  !     ijob=41!! B^T*Q(1,M0+1:M0+fpm(23))=>work(1:N,M0+1:M0+fpm(23))
  !     fpm(21)=1 ! prepare reentry
  !     return
  !  end if





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==1) then !! we initialize a new contour integration
     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(14)==2)) then !! condition on fpm(5) unnecessary        
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        if ((fpm(15)==0) .and. (fpm(23)/nb_procs2>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     end if
     !---------------------------for q --------------------------------------
     if ((loop>0).and.(fpm(23)/nb_procs2>=1)) then
        if (.not.((loop==1).and.(fpm(5)==1))) then ! do not include first iteration with initial guess since Q is known to all
           if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
           if ((fpm(15)==0) .and. (fpm(23)/nb_procs2>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        end if
     end if
     !------------------------------------
#endif



!!!!!!!!!!! RESIDUAL INVERSE ITERATION

     !print *,loop,N,M0,fpm(23)
     if ((loop==0).or.((loop>0).and.(rank2>0))) then
        Q(1:N,1:fpm(23))=ZZERO
        if(fpm(15)==0) Q(1:N,M0+1:M0+fpm(23))=ZZERO
     end if



!!!!!!!!!!!!!!for residual inverse iteration
     if ((loop>0).and.(rank2==0)) then
        allocate(zwork_loc(fpm(23)))
        call zfeast_grationalx(Zne,Wne,fpm(8),lambda,fpm(23),zwork_loc)

        do i=1,fpm(23)
           call ZSCAL(N,zwork_loc(i), q(1,i), 1)
           if (fpm(15)==0) call ZSCAL(N,conjg(zwork_loc(i)), q(1,M0+i), 1)
        end do
        deallocate(zwork_loc)

     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF



!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     !   IF (info==0) then !! will end up checking info errors returned by FEAST drivers
     do e=fpm(20)+rank2,fpm(22),nb_procs2 !!!! loop over the full contour

        if (ijob==-2) then !!Factorize the linear system

           ! set up factorization id (if fact. storage is used)
           if (fpm(10)==1) then
              fpm(33)=fpm(33)+1
           else
              fpm(33)=1
           endif

           Ze = Zne(e)   
           fpm(20)=e-rank2

           ijob=10 ! for fact
           fact=.false. ! action to factorize
           if ((loop==0).or.((loop==1).and.(fpm(5)==1))) then
              fact=.true. ! always factorize in initial loop. Remark: with input initial guess, loop=1 is the initial loop
           else !(not initial loop)
              if ((fpm(22)>nb_procs2).and.(fpm(10)==0)) fact=.true. ! re-factorize if #mpi<#contour points and fact not stored
           end if
           if (fact) then
              call system_clock(t1,tim)
              return 
           end if

!!$           if ((loop==0).or.(fpm(22)>nb_procs2)) then !no refactorization if one linear system per processor
!!$              if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) then !no refactorization if stored in memory
!!$                 if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) then 
!!$                    call system_clock(t1,tim)
!!$                    return 
!!$                 endif
!!$              end if
!!$           end if

        end if
!!!!!!!!!!!!!!!! 1st factorization
        if (ijob==10) then !!Solve the linear system (complex) (zS-A)qr=workc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           call ZLACPY( 'F', N, fpm(23), work, N, workc, N )
           ijob=11 ! for solve
           call system_clock(t1,tim)
           return
        end if

!!!!!!!!!!!!!!!! solve
        if (ijob==11) then 
!!!!!!!!!!!!!!!!
           !! summation
           aux = Wne(e) 
           !----------------
           if (loop>0) then ! for residual inverse iteration
              do i=1,fpm(23)
                 call ZSCAL(N,ZONE/(ze-lambda(i)), workc(1,i), 1)
              enddo
           endif
           !---------------

           call ZAXPY(N*fpm(23),aux,workc,1,Q(1,1),1)

           if ((fpm(15)==1).or.(fpm(15)==2)) then
              ijob=-2 !no solve for ql if complex symmetric or 1 single contour
           elseif (fpm(15)==0) then 
              ijob=20 ! for fact
              fact=.false. ! action to factorize
              if ((loop==0).or.((loop==1).and.(fpm(5)==1))) then
                 fact=.true. ! always factorize in initial loop. Remark: with input initial guess, loop=1 is the initial loop
              else !(not initial loop)
                 if ((fpm(22)>nb_procs2).and.(fpm(10)==0)) fact=.true. ! re-factorize if #mpi<#contour points and fact not stored
              end if
              if (fact) then
                 call system_clock(t1,tim)
                 return 
              end if


!!$              if ((loop==0).or.(fpm(22)>nb_procs2)) then
!!$                 !no refactorization if one linear system per processor
!!$                 if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) then !no refactorization if stored in memory
!!$                    if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) then 
!!$                       call system_clock(t1,tim)
!!$                       return 
!!$                    end if
!!$                 end if
!!$              end if
!!$              
           end if
        end if   !ijob==11

!!!!!!!!!!!! 2nd factorization
        if (ijob==20) then !!Solve the linear system (complex) (zS-A)^Hql=workc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
           call ZLACPY( 'F', N, fpm(23),work(1,M0+1) , N, workc, N )
           call system_clock(t1,tim)
           ijob=21 ! for solve
           return 
        end if

!!!!!!!!!!!! solve with transpose conjugate
        if (ijob==21) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
           aux = conjg(Wne(e)) 
           if (loop>0) then ! for residual inverse iteration
              do i=1,fpm(23)
                 call ZSCAL(N,ZONE/conjg(ze-lambda(i)), workc(1,i), 1)
              enddo
           endif
           !---------------
           call ZAXPY(N*fpm(23),aux,workc,1,Q(1,M0+1),1)

           ijob=-2 ! default identification
        end if
     end do

     !     end IF !! info=0

     !------------------------------------------------
     !#ifdef MPI
     !     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
     !#endif
     !-----------------------------------------------  
     !    if (info/=0) fpm(21)=100 ! the end
     if (info==0) then
        fpm(21)=4 
        !------------------------------------------------
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        if(fpm(15)==0) call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if

  end IF    ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !print *,sum(abs(q(1:N,1)))!,sum(abs(q(1:N,M0+1)))
  !info=-99
  !ijob=0
  !return
  ! stop


  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only qr,ql vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic Estimate
  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only qr is returned with stochastic estimation in res
     iseed=(/56,890,3456,2333/)
     iseed(1) = iseed(1) * (1+rank3*13)
     iseed(2) = iseed(2) + rank3*7
     iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
     iseed(4) = iseed(4) + rank3*2

     call ZLARNV(2,iseed,N*fpm(23),work(1,1)) 
     call ZGEMM('C','N',fpm(23),fpm(23),N,-ZONE,work,N,q,N,ZZERO,zAq,M0) ! projection
     call ZGEMM('C','N',fpm(23),fpm(23),N,ZONE,work,N,work,N,ZZERO,zBq,M0) ! normalization
     !----------------------------------------  
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !---------------------------------------

     Ntotal=fpm(26)
     theta=DZERO
     do i=1,fpm(23)
        theta=theta+dble(zAq(i,i)/zBq(i,i))
        res(i)=abs(theta*(Ntotal/(DONE*i)))
        if (com) then
           write(fout,'(I3)',advance='no') i
           write(fout,'(3X)',advance='no')
           write(fout,'(F10.2)',advance='no') res(i)
           write(fout,*)
        endif
     enddo
     mode=int(real(res(fpm(23))))
     if (com) then
        write(fout,'(A)',advance='no') '# estimated '
        write(fout,'(3X)',advance='no')
        write(fout,'(I6)',advance='no') mode
        write(fout,*)
     endif
     info=5
     if (info/=0) fpm(21)=100 ! The End
  end if


!!!!!!!!!!!!!! Orthogonalize Q  
!!$  if ((fpm(21)==4).and.(fpm(15)==1)) then
!!$ allocate(tau(M0))
!!$     lwork_loc=M0
!!$     allocate(zwork_loc(lwork_loc))
!!$call ZGEQRF( N, fpm(23), Q, N, TAU, zWORK_loc, LWORK_loc, INFO_lap) 
!!$call ZUNGQR( N, fpm(23), fpm(23), Q, N, TAU, zWORK_loc, LWORK_loc, INFO_lap)
!!$     deallocate(tau)
!!$     deallocate(zwork_loc)
!!$end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! Aq xq=eq Bq xq     with Aq=Q^TAQ Bq=Q^TAQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!form  Bq=> Bq=Q^T B Q
  if (fpm(21)==4) then     
     !------------------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZZERO
#endif
     !-------------------------------------------------
     fpm(21)=5 ! preparing reenty
     ijob=40
     call system_clock(t1,tim)
     return! mat-vec B.qr => workr
  end  if

  if (fpm(21)==5) then
     ! after first loop zBq(1,1) contains old residual
     if (loop>0) Ze=zBq(1,1) !  dummy variable
     !------------------------------------------------
#ifdef MPI
     zBq(1:M0,1:fpm(23))=ZZERO
#endif
     !-------------------------------------------------

     fpm(21)=51

     if (fpm(13)>=2) then ! customize inner product
        ijob=50
        call system_clock(t1,tim)
        return
     endif


     if(fpm(15)==1)then 
        call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0)
     elseif(fpm(15)==2)then 
        call ZGEMM('T','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0)
     else ! default fpm(15)=0
        call ZGEMM('C','N',fpm(23),fpm(35),N,ZONE,Q(1,M0+1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0)
     end if

  end if



  if (fpm(21)==51) then
!!!!!!!!!!!!!!  
     if (fpm(13)>=2) call ZCOPY(M0*fpm(25),zAq(1,fpm(24)),1,zBq(1,fpm(24)),1)
!!!!!!!!!!!!!!   

     !---------------------------------------
#ifdef MPI
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code) 
        call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     end if
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !       !---------------------------------------
!!!fpm(21)=6 !! here
  end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! Spurious test 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==51) then 
     mode=0
     if ((fpm(38)==1) .and. (fpm(37)==0)) then
        if (loop>0) then
           allocate(zwork_loc(fpm(23)))
           call zfeast_grationalx(Zne,Wne,fpm(8),lambda,fpm(23),zwork_loc)
           do i=fpm(24),fpm(24)+fpm(25)-1
              if (res(i)>100d0) then ! indicator for being inside the search interval
                 if (abs(sqrt(zBq(i,i))-(zwork_loc(i)))/sqrt(abs(zBq(i,i)))<0.1d0) mode=mode+1
              end if
           end do
           deallocate(zwork_loc)
        endif
     end if
!!!!!!!!!!!!!!!!!!!!
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,mode,1,MPI_INTEGER,MPI_SUM,L2_COMM_WORLD,code)
#endif
!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! B-biorthogonalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !-------------------------
     ! Diagonalize Bq=Ql^H B Qr
     !-------------------------

     e=0
     allocate(VR(fpm(23),fpm(23)))
     allocate(VL(fpm(23),fpm(23)))
     allocate(LALPHA(fpm(23)))

     if ((rank2==0).and.(rank3==0)) then
        allocate(zBqo(fpm(23),fpm(23)))
        call ZLACPY( 'F', fpm(23), fpm(23),zBq , M0, zBqo, fpm(23) )

        LWORK_LOC=8*fpm(23) !! for lapack eig reduced system
        allocate(zWORK_LOC(LWORK_LOC))
        allocate(WORK_LOCR(2*fpm(23)))


        JOBVL='V'
        if ((fpm(15)==1).or.(fpm(15)==2)) JOBVL='N'
        !print *,'ok',zBqo(1:M0,1:M0)
        call ZGEEV(JOBVL,'V',fpm(23),zBqo,fpm(23),LALPHA,VL,fpm(23),VR,fpm(23),zWORK_LOC,LWORK_LOC,WORK_LOCR,info_lap)


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
           if (fpm(15)==0) call ZSWAP(fpm(23),VL(1,m_min),1,VL(1,fpm(23)-i+1),1) 
           jac = LALPHA(fpm(23)-i+1)
           LALPHA(fpm(23)-i+1) = LALPHA(m_min)
           LALPHA(m_min)=jac
        enddo

        deallocate(zWORK_LOC)
        deallocate(WORK_LOCR)
        deallocate(zBqo)
     endif ! rank==0

     !-------------
#ifdef MPI
     if (rank3==0) then
        call MPI_BCAST(e,1,MPI_INTEGER,0,L2_COMM_WORLD,code)
        call MPI_BCAST(VR,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
        if (fpm(15)==0)  call MPI_BCAST(VL,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
        call MPI_BCAST(lalpha,fpm(23),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
     end if

     if (fpm(49)/=0) then
        call MPI_BCAST(e,1,MPI_INTEGER,0,L3_COMM_WORLD,code)
        call MPI_BCAST(VR,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
        if (fpm(15)==0)  call MPI_BCAST(VL,fpm(23)*fpm(23),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
        call MPI_BCAST(lalpha,fpm(23),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
     end if
#endif

     !--------------
     ! Create Ur and Ul (bi-B-orthogonal)
     !--------------------
     IF (fpm(36)==1) then  

        call ZLACPY('F',N,fpm(23),Q(1,1),N,workc(1,1),N)
        Q(1:N,1:fpm(23))=ZZERO
        call ZGEMM('N','N',N,fpm(25),fpm(23),ZONE,workc(1,1),N,VR(1,fpm(24)),fpm(23),ZZERO,Q(1,fpm(24)),N)
        call ZLACPY('F',N,fpm(23),work(1,1),N,workc(1,1),N)
        work(1:N,1:fpm(23))=ZZERO
        call ZGEMM('N','N',N,fpm(25),fpm(23),ZONE,workc(1,1),N,VR(1,fpm(24)),fpm(23),ZZERO,work(1,fpm(24)),N)
        if(fpm(15)==0) then
           call ZLACPY('F',N,fpm(23),Q(1,M0+1),N,workc(1,1),N)
           Q(1:N,M0+1:M0+fpm(23))=ZZERO
           call ZGEMM('N','N',N,fpm(35),fpm(23),ZONE,workc(1,1),N,VL(1,fpm(24)),fpm(23),ZZERO,Q(1,fpm(34)),N)
        endif

        !--------------------
#ifdef MPI
        if (fpm(23)/nb_procs2>=1) then
           call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
           call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
           if (fpm(15)==0) call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        endif
#endif

        !----------------------------------------
        !---------- Normalization of Ur and Ul
        !----------------------------------------          
        call ZGEMM('N','N',fpm(23),fpm(23),fpm(23),ZONE,zBq(1,1),M0,VR(1,1),fpm(23),ZZERO,zAq(1,1),M0)

        if (fpm(15)==1) then
           call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ZONE,VR(1,1),fpm(23),zAq(1,1),M0,ZZERO,zBq(1,1),M0)
        elseif (fpm(15)==2) then
           call ZGEMM('T','N',fpm(23),fpm(23),fpm(23),ZONE,VR(1,1),fpm(23),zAq(1,1),M0,ZZERO,zBq(1,1),M0)
        else ! default fpm(15)=0
           call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ZONE,VL(1,1),fpm(23),zAq(1,1),M0,ZZERO,zBq(1,1),M0)
        endif

        do i=1,fpm(23)
           call ZSCAL(N,ZONE/sqrt(zBq(i,i)),Q(1,i),1)
           call ZSCAL(N,ZONE/sqrt(zBq(i,i)),work(1,i),1)
           if(fpm(15)==0) call ZSCAL(N,ZONE/sqrt(conjg(zBq(i,i))),Q(1,M0+i),1)
        enddo

     end IF ! fpm(36)==1

     deallocate(LALPHA)
     deallocate(VR)
     deallocate(VL)

     if ((e/=0) .and. (loop>0))then
        !     if (e/=0) then
        fpm(23) = fpm(23)-e
        fpm(25)=fpm(23) ! default - single proc
        fpm(35)=fpm(23) ! default - single proc
        if (com) then 
           write(fout,'(A)',advance='no') 'Resize Subspace (Bq)'
           write(fout,'(3X)',advance='no')
           write(fout,'(I4)',advance='no')fpm(23)
           write(fout,*)
        endif


        if (fpm(23)/nb_procs2>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(24)=1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 

           fpm(35)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(35)=fpm(35)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(34)=M0+1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
        end if

     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !--------------------------------
#ifdef MPI
     zBq(1:M0,1:fpm(23))=ZZERO
#endif
     !--------------------------------
     fpm(21)=6

     if (fpm(13)>=2) then ! customize inner product
        ijob=50
        call system_clock(t1,tim)
        return
     endif

     if (fpm(15)==1)then 
        call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0)
     elseif (fpm(15)==2)then 
        call ZGEMM('T','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0)
     else ! default fpm(15)=0
        call ZGEMM('C','N',fpm(23),fpm(35),N,ZONE,Q(1,M0+1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0)
     endif

  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!form Aq=> Aq=Q^T A Q 
  if (fpm(21)==6) then

!!!!!!!!!!!!!!  
     if (fpm(13)>=2) call ZCOPY(M0*fpm(25),zAq(1,fpm(24)),1,zBq(1,fpm(24)),1)
!!!!!!!!!!!!!!   
     !------------------------------------------------
#ifdef MPI 
     work(1:N,1:fpm(23))=ZZERO
#endif 
     fpm(21)=7 ! preparing reentry
     ijob=30
     call system_clock(t1,tim)
     return  ! mat-vec A.qr => workr
  endif


  if (fpm(21)==7) then 
     !------------------------------------------------
#ifdef MPI 
     zAq(1:M0,1:fpm(23))=ZZERO
#endif
     !-------------------------------------------------

     fpm(21)=8
     if (fpm(13)>=2) then ! customize inner product
        call system_clock(t1,tim)
        ijob=50
        return
     endif

     if (fpm(15)==1)then
        call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zAq(1,fpm(24)),M0)
     elseif (fpm(15)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zAq(1,fpm(24)),M0)
     else ! default fpm(15)=0
        call ZGEMM('C','N',fpm(23),fpm(35),N,ZONE,Q(1,M0+1),N,work(1,fpm(24)),N,ZZERO,zAq(1,fpm(24)),M0)
     endif


  endif



  if (fpm(21)==8) then
     !----------------------------------------!(zAq,zBq known to all processors) 
#ifdef MPI 
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     end if
     if (fpm(49)/=0) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
     end if
#endif
     !---------------------------------------
     if  ((fpm(13)==1).or.(fpm(13)==3)) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=60
        call system_clock(t1,tim)
        return
     endif
     fpm(21)=81
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==81) then
     ! if using FEAST-MPI ==> solve on a single proc
     if ((rank2==0).and.(rank3==0)) then
        info_lap=1 ! initialization
        LWORK_LOC=8*fpm(23) !! for lapack eig reduced system
        allocate(zWORK_LOC(LWORK_LOC))
        allocate(WORK_LOCR(8*fpm(23)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(zBqo(fpm(23),fpm(23)))
        allocate(VL(fpm(23),fpm(23)))
        allocate(VR(fpm(23),fpm(23)))


        call ZLACPY( 'F', fpm(23), fpm(23),zBq , M0, zBqo, fpm(23) )

        if (fpm(39)==1) then  !!! standard eigenvalue
           do i=1,fpm(23)
              zAq(1:fpm(23),i) = zAq(1:fpm(23),i) / sqrt(zBq(i,i))
              zAq(i,1:fpm(23)) = zAq(i,1:fpm(23)) / sqrt(zBq(i,i))
           enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           JOBVL='V'
           if ((fpm(15)==1).or.(fpm(15)==2)) JOBVL='N'  
           call ZGEEV(JOBVL,'V',fpm(23),zAq,M0,LAMBDA,VL,fpm(23),VR,fpm(23), zWORK_LOC, LWORK_LOC, WORK_LOCR,INFO_lap)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           if(info_lap==0)then
              do j=1,fpm(23)
                 call ZSCAL(fpm(23),ZONE/sqrt(zBq(j,j)),VR(j,1),fpm(23))
                 if (fpm(15)==0)  call ZSCAL(fpm(23),ZONE/conjg(sqrt(zBq(j,j))),VL(j,1),fpm(23))
              end do
           endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
        else                !!! generalized eigenvalue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           allocate(LALPHA(fpm(23)))
           allocate(LBETA(fpm(23))) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           JOBVL='V'
           if ((fpm(15)==1).or.(fpm(15)==2)) JOBVL='N'  
           call ZGGEV(JOBVL,'V',fpm(23),zAq,M0,zBqo,fpm(23),LALPHA,LBETA,VL,fpm(23),VR,fpm(23), zWORK_LOC, LWORK_LOC, WORK_LOCR,INFO_lap)
           LAMBDA(1:fpm(23))=LALPHA(1:fpm(23))/LBETA(1:fpm(23))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           deallocate(LALPHA)
           deallocate(LBETA)
        endif

        if (info_lap/=0) info=-3 ! problem with lapack
        !print *,info_lap


        IF (info_lap==0) then

           call ZGEMM('N','N',fpm(23),fpm(23),fpm(23),ZONE,zBq(1,1),M0,VR(1,1),fpm(23),ZZERO,zAq(1,1),M0)

           if(fpm(15)==1) then 
              call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ZONE,VR(1,1),fpm(23),zAq(1,1),M0,ZZERO,zBqo(1,1),fpm(23))
           elseif(fpm(15)==2) then 
              call ZGEMM('T','N',fpm(23),fpm(23),fpm(23),ZONE,VR(1,1),fpm(23),zAq(1,1),M0,ZZERO,zBqo(1,1),fpm(23))
           else ! default fpm(15)=0
              call ZGEMM('C','N',fpm(23),fpm(23),fpm(23),ZONE,VL(1,1),fpm(23),zAq(1,1),M0,ZZERO,zBqo(1,1),fpm(23))
           endif


           fpm(37)=0
           do j=1,fpm(23) !!!!!! Check orthogonality of ZGGEV eigenvectors
              zalpha = abs(zBqo(j,j))
              do jj=1,fpm(23)
                 zalpha = zalpha - abs(zBqo(jj,j))
              enddo
              !if (log10(abs(zalpha)/(fpm(23)))>-10) then
              if (abs(zalpha)/(fpm(23))>1d-10) then
                 fpm(37) = 1 
              endif
           enddo

           if (fpm(37)==0) then   
              do j=1,fpm(23)
                 call ZSCAL(fpm(23),ZONE/sqrt(zBqo(j,j)),VR(1,j),1)
                 if(fpm(15)==0) call ZSCAL(fpm(23),ZONE/conjg(sqrt(zBqo(j,j))),VL(1,j),1)
              enddo
           endif

           call ZLACPY( 'F', fpm(23), fpm(23),VR , fpm(23), zAq, M0 )
           if (fpm(15)==0) call ZLACPY( 'F', fpm(23), fpm(23),VL , fpm(23), zBq, M0 )

        end IF

        ! print *,zAq(1:fpm(23),1)
        !  print *,zBq(1:fpm(23),1)
        !        stop 




        deallocate(zBqo)
        deallocate(VL)
        deallocate(VR)

        deallocate(zwork_loc)
        deallocate(work_locr)
     end if !(rank 0)
     !-------------------------------- !(info common to all processors) 
#ifdef MPI
     if (rank3==0) call MPI_BCAST(info,1,MPI_INTEGER,0,L2_COMM_WORLD,code)
     if (fpm(49)/=0) call MPI_BCAST(info,1,MPI_INTEGER,0,L3_COMM_WORLD,code)
#endif
     !--------------------------------
     if (info/=0) fpm(21)=100 ! the end


     if (info==0) then

        !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 
#ifdef MPI
        if (rank3==0) then 
           call MPI_BCAST(zAq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
           call MPI_BCAST(zBq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
           call MPI_BCAST(fpm(37),1,MPI_INTEGER,0,L2_COMM_WORLD,code)
        end if
        if (fpm(49)/=0) then 
           call MPI_BCAST(zAq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
           call MPI_BCAST(zBq,M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
           call MPI_BCAST(fpm(37),1,MPI_INTEGER,0,L3_COMM_WORLD,code)
        end if
#endif
        !-----------------------------

        fpm(21)=9
     end if
  end if !! fpm(21)=81

  !  print *,zAq(1:fpm(23),1)
  !  print *,zBq(1:fpm(23),1)
  !        stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors Xr=Qr*xr 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==9 ) then
     !! from precious step workr=A*Qr, before using workr memory space
     !! let us construct first workc=workr*Aq which is workc=A*Xr
#ifdef MPI
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
#endif

     call ZGEMM('N','N',N,fpm(25),fpm(23),ZONE,work(1,1),N,zAq(1,fpm(24)),M0,ZZERO,workc(1,fpm(24)),N)

     !! construct X (q needed for outer-loop as well)
     !call ZLACPY( 'F', N, fpm(25),q(1,fpm(24)) , N, work(1,fpm(24)), N )
     call ZLACPY( 'F', N, fpm(23),q(1,1) , N, work(1,1), N )

#ifdef MPI 
     Q(1:N,1:fpm(23))=ZZERO
#endif

     call ZGEMM('N','N',N,fpm(25),fpm(23),ZONE,work(1,1),N,zAq(1,fpm(24)),M0,ZZERO,Q(1,fpm(24)),N)

     fpm(21)=10 !! no need to compute A*qr again
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residuals ||AXr-lambda*BXr|| 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==10) then
     !handle special case of input initial guess

     !if ((fpm(5)==1).and.(loop==0))  call ZLACPY( 'F', N, fpm(23),work , N, workc, N )
     if ((fpm(5)==1).and.(loop==0))  call ZLACPY( 'F', N, fpm(25),work(1,fpm(24)) , N, workc(1,fpm(24)), N )

     fpm(21)=11 ! preparing reentry
     ijob=40 
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZZERO  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     call system_clock(t1,tim)
     return  ! mat-vec B*qr => workr
  endif

  if (fpm(21)==11) then
     !----------------------------------------
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif
     !---------Compute residual vector R=AX-BXE
     !         work<=workc-work*E (needed also for outer-loop)
     do i=fpm(24),fpm(24)+fpm(25)-1
        !res(i)=sum(abs(work(1:N,i))) ! placeholder (denominator)
        res(i)=DZNRM2(N,work(1,i),1)**2 ! square because of mpi compatibility sqrt(mpi1)+sqrt(mpi2) /= sqrt(mpi1+mpi2) 
        call ZSCAL(N,-lambda(i),work(1,i),1)
     enddo

     call ZAXPY(N*fpm(25),ZONE,workc(1,fpm(24)),1,work(1,fpm(24)),1)

     !reduce first denominator
#ifdef MPI 
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif


     do i=fpm(24),fpm(24)+fpm(25)-1
        !res(i)=sum(abs(work(1:N,i)))/(abs(abs(Emid)+r)*res(i))
        res(i)=(DZNRM2(N,work(1,i),1)**2)/((abs(abs(Emid)+r)**2)*res(i))
     end do

     !print *,work(1:N,1),res(1)
     !stop
     !do i=1,M0
     !print *,'res',i,res(i)
     !enddo
     !stop

     !----------------------------------------
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !-----------------------------------------
     do i=1,fpm(23)
        res(i)=sqrt(res(i)) ! norm L2
     enddo

     fpm(21)=-9 
     if ((fpm(15)==1).or.(fpm(15)==2)) fpm(21)=-12
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors Xl=Ql*xl
!!!!!  Compute Residuals ||A^H Xl-conjg(lambda)*B^H Xl|| 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==-9 ) then

     if (.not.((fpm(5)==1).and.(loop==0))) then ! no first time initial guess
        call ZLACPY( 'F', N, fpm(23),Q(1,M0+1) , N, workc(1,1), N )
#ifdef MPI 
        Q(1:N,M0+1:M0+fpm(23))=ZZERO
#endif
        call ZGEMM('N','N',N,fpm(35),fpm(23),ZONE,workc(1,1),N,zBq(1,fpm(24)),M0,ZZERO,Q(1,fpm(34)),N) !!<<
     end if

!!!!!!!!!!!!!!!!!
     fpm(21)=-10 ! preparing reentry  
     ijob=31
     !#ifdef MPI
     !     work(1:N,M0+1:M0+fpm(23))=ZZERO   
     !#endif
     call system_clock(t1,tim)  
     return  ! mat-vec A^C.ql => work(1,fpm(34))
  endif

  if (fpm(21)==-10) then    
     !print *,'q',q(1:N,M0+1)
     !     print *,'A*q',work(1:N,M0+1)

     call ZLACPY( 'F', N, fpm(35),work(1,fpm(34)), N, workc(1,fpm(24)), N ) ! workc = A^C.XL
     fpm(21)=-11 ! preparing reentry
     ijob=41 
     !----------------------------------------
#ifdef MPI
     work(1:N,M0+1:M0+fpm(23))=ZZERO  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     call system_clock(t1,tim)
     return  ! mat-vec S^C.ql => work
  endif

  if (fpm(21)==-11) then

     !print *,'B*q',work(1:N,M0+1)

     !----------------------------------------
#ifdef MPI
     res(M0+1:M0+fpm(23))=DZERO
#endif
     !---------Compute residual vector R=AX-BXE
     !         work<=workc-work*E (needed also for outer-loop)
     do i=fpm(34),fpm(34)+fpm(35)-1
        !  res(i)=sum(abs(work(1:N,i))) ! placeholder (denominator)
        res(i)=DZNRM2(N,work(1,i),1)**2 ! square because of mpi compatibility sqrt(mpi1)+sqrt(mpi2) /= sqrt(mpi1+mpi2) 
        call ZSCAL(N,-conjg(lambda(i-M0)),work(1,i),1)
     enddo
     !    print *,'B*q-scaled',work(1:N,M0+1),-conjg(lambda(1))
     !stop 
     call ZAXPY(N*fpm(35),ZONE,workc(1,fpm(34)-M0),1,work(1,fpm(34)),1)
     !reduce first denominator
#ifdef MPI 
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif

     !-------- Residual
     do i=fpm(34),fpm(34)+fpm(35)-1
        !res(i)=sum(abs(work(1:N,i)))/(abs(abs(Emid)+r)*res(i))
        res(i)=(DZNRM2(N,work(1,i),1)**2)/((abs(abs(Emid)+r)**2)*res(i))
     end do


     !print *,work(1:N,M0+1),res(M0+1)
     !stop    

     !----------------------------------------
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !-----------------------------------------
     do i=1,fpm(23)
        res(M0+i)=sqrt(res(M0+i)) ! norm L2
     enddo

     fpm(21)=-12
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==-12) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !---color mapping --> check if eigenvalue inside the contour
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     allocate(color(fpm(23)))

     if (fpm(29)==0) then ! expert custom contour
        call zfeast_inside_contourx(Zne,fpm(8),lambda,fpm(23),color)
     else
        call zfeast_inside_contour(Emid,r,fpm(18),fpm(19),lambda,fpm(23),color)
     endif
!!!!!!!!!!!!

     k=0
     trace=ZZERO
     theta=DZERO
!!! count how many eigenvalues have converged + trace + max residual
     do i=1,fpm(23)
        if (color(i)==1) then ! inside the search interval
           k=k+1
           trace=trace+lambda(i)
           if (res(i)>theta) theta=res(i) ! max residual (right residual)
           if (fpm(15)==0) then
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
           if (fpm(15)==0) res(e+M0)=-DONE !spurious
        end do
!!! max residual
        theta=DZERO
        do i=1,fpm(23)
           if (color(i)==1) then ! inside the search interval
              if (res(i)>theta) theta=res(i) ! max residual (right residual)
              if (fpm(15)==0) then
                 if (res(M0+i)>theta) theta=res(M0+i) ! max residual (left residual)
              end if
           end if
        enddo
     end if

     !Ntotal=fpm(26)
     if (loop>2) then! wait third iteration (spurious related+ifeast)
        if (mode==0) info=1  ! no eigenvalue detected in the interval
        !if (M0/=1) then !!!!!!<<<<<<<<<<
        if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small !!!<<<<<
        !endif        
     endif
     ! print *,'info',rank3,info,N
     if (info/=0) then
        fpm(21)=100 ! The End
        deallocate(color)
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
        !if ((fpm(6)==0).and.(log10((epsout))<(-fpm(3)))) testconv=.true.
           if ((fpm(6)==0).and.(epsout<10.0d0**(-fpm(3))))  testconv=.true.
     end if

     !! residual convergence
     !if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.
 if ((fpm(6)/=0).and.(theta<10.0d0**(-fpm(3))))  testconv=.true.
     
     !! Two loops minimum if spurious are found
     if ((loop<=1).and.(k>mode)) testconv=.false.

 
!!!! L2 barrier to guarantee synchronization in printing
#ifdef MPI
 call MPI_BARRIER(L2_COMM_WORLD,code)
#endif 
    

     if (com) then
        if (fpm(5)==0) then
           write(fout,'(I3)',advance='no') loop
        else ! no integration in first loop
           if (loop==0) then
              write(fout,'(A)',advance='no') 'N/A'
           else
              write(fout,'(I3)',advance='no') loop-1
           end if
        endif
        write(fout,'(3X)',advance='no')
        write(fout,'(I4)',advance='no') mode
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') abs(trace) !!<<<<
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') epsout
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') theta
        write(fout,*)
     end if

     if (.not.testconv) then
        ZBq(1,1)=trace ! dummy
        e=loop
        if (fpm(5)==1) e=e-1
        if (e==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif

     !! handle a special case (prevent convergence if no eigenvalue found)
     !! for mode=0 FEAST will exit after few iterations (see above)
     if (mode==0) testconv=.false. 

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
        allocate(zwork_loc(fpm(23)))
        call zfeast_grationalx(Zne,Wne,fpm(8),lambda,fpm(23),zwork_loc)
        do i=1,fpm(23)
           res(i)=abs(zwork_loc(i)) !! save them all          
           if (color(i)==1)  res(i)=200d0+res(i)  !!could be used for rci commands (e.g. bicgstab)
        enddo
        deallocate(zwork_loc)
        deallocate(color)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Remark: q=x and work=Ax-eBx is known already, using FEAST-MPI q and work are already distributed   (lef,right vector components are known)     
        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        if (fpm(10)==1) fpm(33)=0 ! initialize factorization id 
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    ! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

#ifdef MPI       
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        if (fpm(15)==0) then
           if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        endif
#endif

        !! reorder (shift) lambda, eigenvector and residual
        call ZCOPY( N*fpm(23),Q(1,1), 1, work(1,1), 1 )
        if(fpm(15)==0) call ZCOPY( N*fpm(23),Q(1,M0+1), 1, work(1,fpm(23)+1), 1 )
        allocate(zwork_loc(fpm(23)))
        if (fpm(15)==0) then
           allocate(work_loc2(2*fpm(23))) 
        else
           allocate(WORK_LOC2(fpm(23))) 
        end if
        call ZCOPY(fpm(23),lambda , 1, zwork_loc, 1 )
        call DCOPY(fpm(23),res , 1, work_loc2, 1 )
        if(fpm(15)==0) call DCOPY(fpm(23), res(M0+1), 1, work_loc2(fpm(23)+1), 1 )

        M0=fpm(23)  ! update value of M0 (new subspace)

        Q(1:N,1:fpm(23))=ZZERO
        if(fpm(15)==0) Q(1:N,M0+1:M0+fpm(23))=ZZERO
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
                 if(fpm(15)==0) then
                    q(1:N,M0+e)=work(1:N,M0+i)
                    res(M0+e)=work_loc2(M0+i)
                 endif
              endif
           else
              j=j+1
              q(1:N,mode+j)=work(1:N,i)
              lambda(mode+j)=zwork_loc(i)
              res(mode+j)=work_loc2(i)
              if(fpm(15)==0) then
                 q(1:N,M0+mode+j)=work(1:N,M0+i)
                 res(M0+mode+j)=work_loc2(M0+i)
              endif
           end if
           if (work_loc2(i)==-DONE) then ! spurious at the end
              k=k+1
              q(1:N,fpm(23)-k+1)=work(1:N,i)
              lambda(fpm(23)-k+1)=zwork_loc(i)
              res(fpm(23)-k+1)=work_loc2(i)
              if(fpm(15)==0) then
                 q(1:N,M0+fpm(23)-k+1)=work(1:N,M0+i)
                 res(M0+fpm(23)-k+1)=work_loc2(M0+i)
              end if
           endif
        end do

        deallocate(zwork_loc)
        deallocate(work_loc2)
        deallocate(color)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=100 ! The End

        !if ((info==0).and.(fpm(37)==1).and.(fpm(15)==0)) info=6 ! only if left eigenvectors are computed
        if ((info==0).and.(fpm(37)==1)) info=6 

     end if ! test convergence
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==100) then !! THE END (ijob=0) 

     ijob=0 !! exit FEAST

     if (com) then !! Print  Information

        if (info>=200) then
           write(fout,'(A)',advance='no') 'PROBLEM with input parameters'
           write(fout,*) 
        end if

        if (info==-3) then
           write(fout,'(A)',advance='no') 'ERROR with reduced system'
           write(fout,*) 
        end if

!!$          if (info==-2) then
!!$             write(fout,'(A)',advance='no') 'ERROR from Inner Linear System Solver in FEAST driver'  
!!$             write(fout,*) 
!!$          end if
!!$
!!$          if (info==-1) then
!!$             write(fout,'(A)',advance='no') 'ERROR with Internal memory allocation'  
!!$             write(fout,*) 
!!$          end if

        if (info==1) then
           write(fout,'(A)',advance='no') '==>WARNING: No eigenvalue has been found in the proposed search interval'
           write(fout,*)
        endif

        if (info==3) then
           write(fout,'(A)',advance='no') '==>WARNING: Size subspace M0 too small'  
           write(fout,*)
        end if

        if (info==4) then
           write(fout,'(A)',advance='no') '==>WARNING: Only the subspace has been returned'  
           write(fout,*)
        end if

        if (info==5) then
           write(fout,'(A)',advance='no') '==>WARNING: Only stochastic estimation of #eigenvalues returned'  
           write(fout,*)
        end if

        if (info==6) then
           if (fpm(15)==0) then
              write(fout,'(A)',advance='no') '==>WARNING: FEAST converges but left/right subspaces are not bi-orthonormal'
           else
              write(fout,'(A)',advance='no') '==>WARNING: FEAST converges but right subspace is not orthonormal'
           end if
           write(fout,*) 
        end if

        if (info==2) then
           write(fout,'(A)',advance='no') '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)'  
           write(fout,*)
        end if


        if (info==0) then
            write(fout,*)
            write(fout,'(A)',advance='no') '==>FEAST has successfully converged'
           if (fpm(6)==0) then  
              write(fout,'(A)',advance='no') ' with Trace tolerance <1E-'
           else
               write(fout,'(A)',advance='no') ' with Residual tolerance <1E-'
            end if
             if (fpm(3)<10) then
               write(fout,'(I1)') fpm(3)
            else
               write(fout,'(I2)') fpm(3)
            end if
              write(fout,'(A)',advance='no') '   # FEAST outside it. '
             write(fout,'(I8)') loop
             if (mod(fpm(30),10000)/1000==2) then !! ifeast
                 write(fout,'(A)',advance='no') '   # Inner BiCGstab it.'
                write(fout,'(I8)',advance='no') fpm(60)
                if (fpm(30)/100000==2)    write(fout,'(A)',advance='no') ' for rank2=0'  !pfeast
                 write(fout,*)
              end if
               write(fout,'(A)',advance='no') '   # Eigenvalue found  '
           write(fout,'(I8)') mode
        else
           write(fout,'(A)',advance='no') '==>INFO code = '
           write(fout,'(I4)',advance='no') info
           write(fout,*)
        end if
        

        
        if (info>=0) then             
           call system_clock(tend,tim) !! total time
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           write(fout,*)
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST-RCI timing  |'
           write(fout,'(A)') '--------------------.------------------.'
           write(fout,'(A)',advance='no') '| Fact. cases(10,20)|'
           write(fout,'(F12.4)',advance='no') (t10*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Solve cases(11,12)|'
           write(fout,'(F12.4)',advance='no') (t11*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| A*x   cases(30,31)|'
           write(fout,'(F12.4)',advance='no') (t30*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| B*x   cases(40,41)|'
           write(fout,'(F12.4)',advance='no') (t40*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Misc. time        |'
           write(fout,'(F12.4)',advance='no') (tend-tini-t10-t11-t30-t40)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Total time (s)    |'
           write(fout,'(F12.4)',advance='no') (tend-tini)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)') '--------------------------------------- '
           write(fout,*)
        endif
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,'(A)',advance='no') '*********** FEAST- END*************************'
        write(fout,*) 
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,*)
        !! close file
        if (fpm(1)<0) close(fout)
     endif
#ifdef MPI
     call MPI_BARRIER(L2_COMM_WORLD,code)
     !if (fpm(49)>0)  call MPI_BARRIER(L3_COMM_WORLD,code)
#endif

  end if

end subroutine zfeast_grcix




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_grci(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info)
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
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !  Eric Polizzi- 2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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


  if (ijob==-1) then
     fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     !fpm(15)=0  ! ID for 2-sided feast (default with name)
     if (fpm(30)==-111) then
        fpm(30)=141310 ! code name for direct call to this routine
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)

end subroutine zfeast_grci




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine zfeast_srcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
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
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !  Eric Polizzi - 2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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

  if (ijob==-1) then
     !!fpm(15)=2 ! ID for 1-sided feast for Complex Symmetric Matrix (default with name)
     if (fpm(30)==-111) then
        fpm(30)=141111 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if

  call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)

end subroutine zfeast_srcix




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine zfeast_srci(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info)
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
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !  Eric Polizzi - 2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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


  if (ijob==-1) then
     fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     !!fpm(15)=2 ! ID for 1 sided feast for Complex Symmetric Matrix (default with name)
     if (fpm(30)==-111) then
        fpm(30)=141110 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if


  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info,Zne,Wne)

end subroutine zfeast_srci




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine dfeast_grcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
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
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !  Eric Polizzi 2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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

  if (ijob==-1) then
     !fpm(15)=0 !ID for 2-sided feast (default)
     if (fpm(30)==-111) then
        fpm(30)=121311 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if

  call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)

end subroutine dfeast_grcix





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine dfeast_grci(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info)
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
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !  Eric Polizzi 2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0,*):: zAq,zBq
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


  if (ijob==-1) then
     fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     !fpm(15)=0 ! ID for 2 sided feast (default)
     if (fpm(30)==-111) then
        fpm(30)=121310 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if


  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  call zfeast_grcix(ijob,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)


end subroutine dfeast_grci



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_grcipevx(ijob,dmax,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve polynomial eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      P(lambda)qr=0 and P(lambda)^Hql=0                            
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  ! 
  !  Expert comments: This is also a Kernel routine (contains several cases depending of fpm(15))           
  !
  !                   fpm(15)=0 ! 2 sided algo (default)
  !                   fpm(15)=1 ! 1 sided algo- right eigenvector
  !                   fpm(15)=2 ! 1 sided algo with right/left conjugate- default for complex symmetric
  !
  !  A can be NON-HERMITIAN (could be real or complex)  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,2*M0):   Workspace
  !                            expert comment: size of (N,M0) if fpm(15)>0  
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !                            expert comment: size of (N,M0) if fpm(15)>0 (only right vectors)
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right (1:mode) and left (M0+1:M0+mode)
  !                          
  !                            expert comment: size of (M0) if fpm(15)>0 (only rigth vectors)
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
  !  Eric Polizzi - 2019
  !=====================================================================

  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  !include "f90_noruntime_interface.fi"
  integer :: ijob,N,M0,dmax
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0*dmax,*):: zAq,ZBq
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
  complex(kind=(kind(1.0d0))),parameter :: ZONE=(DONE,DZERO)
  complex(kind=(kind(1.0d0))),parameter :: ZZERO=(DZERO,DZERO)
  double precision, parameter :: ba=-pi/2.0d0,ab=pi/2.0d0
  integer*8,save :: fout
  logical, save :: com
  !! variable for FEAST
  character(len=3) :: ctemp 
  character(len=25) :: name
  integer :: i,e,j,k,Ntotal,jj,m_min
  integer,dimension(4) :: iseed
  double precision :: theta,theta_min
  complex(kind=(kind(1.0d0))) :: jac,aux
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: zBqo
  integer, dimension(:),allocatable :: fpm_default,color,map,map2
  logical :: testconv,test,fact
  complex(kind=(kind(1.0d0))) :: trace, zalpha
  !! Lapack variable 
  double precision, dimension(:),allocatable :: work_locr,work_loc2,dtau
  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: zwork_loc,tau,lambda2
  complex(kind=(kind(1.0d0))), dimension(:),allocatable :: LALPHA,LBETA,lambdar
  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: VL,VR,VRL,Rmat!,bb,zzwork
  integer :: lwork_loc,info_lap,infoloc
  character(len=1) :: JOBVL
  integer,dimension(3000),save :: ipiv !M0????????????
  !! MPI compatibility variables
  integer :: rank,nb_procs,NEW_COMM_WORLD
  !! feast reduced system
  integer,dimension(64) :: fpm2
  integer :: M02,M2,loop2,info2,itemp,k1,k2
  double precision :: epsout2,dtemp
  double precision,dimension(:),allocatable :: res2


  !! MPI compatibility variables
  integer :: rank2,rank3,code,nb_procs2,nb_procs3,L2_COMM_WORLD,L3_COMM_WORLD
  double precision :: DZNRM2
  integer,save :: tini,t1,t10,t11,t30,t31,t50,t60
  integer :: tend,tim,t2
  character(len=8) :: date
  ! double precision, dimension(:,:),allocatable :: dwork
  !  complex(kind=(kind(1.0d0))), dimension(:,:),allocatable :: zzwork

  rank=0 !!<<

!!!!record timing
  if (ijob/=-1) then  
     call system_clock(t2,tim)
     select case(ijob)
     case(10,20)
        t10=t10+(t2-t1)
     case(11,21)
        t11=t11+(t2-t1)
     case(30)
        t30=t30+(t2-t1)
     case(31)
        t31=t31+(t2-t1)
     case(50)
        t50=t50+(t2-t1)
     case(60)
        t60=t60+(t2-t1)
     end select
  end if

!!! parallel init
  rank2=0
  rank3=0
  nb_procs2=1
  nb_procs3=1
  !----------------------------------------------
#ifdef MPI
  L2_COMM_WORLD=fpm(9)
  call MPI_COMM_RANK(L2_COMM_WORLD,rank2,code)
  call MPI_COMM_SIZE(L2_COMM_WORLD,nb_procs2,code)
  L3_COMM_WORLD=fpm(49)
  if (fpm(49)/=0) then
     call MPI_COMM_RANK(L3_COMM_WORLD,rank3,code)
     call MPI_COMM_SIZE(L3_COMM_WORLD,nb_procs3,code)
  end if
#endif
  !---------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Initialization!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ijob==-1) then

     info=0 ! default value  
     if (fpm(30)==-111) then
        fpm(30)=141317 ! code name for direct all to this routine   
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     endif

     call system_clock(tini,tim)
     t10=0
     t11=0
     t30=0
     t31=0
     t50=0
     t60=0
     ! set up factorization id
     fpm(33)=0

     ! expert routines Emid,r needs to be calculated
     if(fpm(29)==0)then 
        ! calculate the center of mass and relative maximum length of custom contour
        Emid=ZZERO
        do i=1,fpm(8)
           Emid = Emid + Zne(i)
        enddo
        Emid = Emid / fpm(8)
        r = abs(maxval(abs(Zne(1:fpm(8))))-Emid) ! for normalization of residual
     end if



     ! comments?
     com=.false. !default
     if (fpm(1)/=0) com=.true.
#ifdef MPI
     if (com.and.((rank2/=0).or.(rank3/=0))) com=.false. ! comment only in rank 0 (only a single mpi is allowed to comment)
#endif
     !file or screen output?
     if (fpm(1)<0) then
        write(ctemp,'(I3)') abs(fpm(1))
        fout=abs(fpm(1))+200!fpm(60) ! file number default
        if (com) then
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append') !! all mpi procs may have access to the file (if needed)
           !write(fout,'(A)',advance='no') '############## NEW FEAST RUN ###################'
           !write(fout,*)
        endif
     elseif (fpm(1)>0) then
        fout=6 !screen
     endif

!!!!!!!!!!!!!!!! Find Ntotal
     Ntotal=N
#ifdef MPI           
     if (fpm(49)/=0) then
        Ntotal=0
        call MPI_ALLREDUCE(N,Ntotal,1,MPI_INTEGER,MPI_SUM,L3_COMM_WORLD,code)
     end if
#endif        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
     fpm(26)=Ntotal


     call dcheck_feast_grci_input(r,M0,Ntotal,info)


     if (info/=0) fpm(21)=100 ! The End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     IF (info==0) then


        if (com) then
           write(fout,*)
           write(fout,'(A)',advance='no') '***********************************************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '*********** FEAST v'
           write(fout,'(F3.1)',advance='no') fpm(31)/10.0d0
           write(fout,'(A)',advance='no') ' BEGIN ******************'
           write(fout,*) 
           write(fout,'(A)',advance='no') '***********************************************'  
           write(fout,*)
         if (fpm(1)<0) then ! only in-file
           call date_and_time(DATE=date)
           write(fout,'(A)',advance='no') 'Date: '//trim(date(1:4))//'-'//trim(date(5:6))//'-'//trim(date(7:8))
           write(fout,*)
        end if
           call feast_name(fpm(30),name)
           write(fout,'(A)') 'Routine '//trim(name)
           write(fout,'(A)',advance='no') 'Solving P(e)X=0 with P polynomial of degree'
           write(fout,'(I2)',advance='no') dmax
           write(fout,*)
           if (nb_procs2*nb_procs3>1) then
              write(fout,'(A)',advance='no') '#MPI (total=L2*L3) '
              write(fout,'(I4)',advance='no') nb_procs2*nb_procs3
              write(fout,'(A)',advance='no') '='
              write(fout,'(I4)',advance='no') nb_procs2
              write(fout,'(A)',advance='no') '*'
              write(fout,'(I4)',advance='no') nb_procs3
              write(fout,*)
           endif
!!!!!!!!!!!! Print the FEAST parameters which has been changed from default
           write(fout,'(A)',advance='no') 'List of input parameters fpm(1:64)-- if different from default'
           write(fout,*)
           allocate(fpm_default(64))
           call feastinit(fpm_default) ! initialize to -111
           fpm_default(1)=0
           fpm_default(30)=fpm(30) ! name
           fpm_default(14)=fpm(14) ! for special conditional cases
           call feastdefault(fpm_default,infoloc) ! get the default values
           do i=1,64   ! compare with possible mismatch with user values
              if ((i<=19).or.((i>=40).and.(i<=50))) then
                 if ((fpm(i)/=fpm_default(i)).and.(i/=9).and.(i/=49)) then
                    write(fout,'(A)',advance='no') '   fpm('
                    write(fout,'(I2)',advance='no') i
                    write(fout,'(A)',advance='no') ')='
                    write(fout,'(I4)',advance='no') fpm(i)
                    write(fout,*)
                 endif
              endif
           enddo
           deallocate(fpm_default)
           fpm(22)=nb_procs3 ! temp copy
           call zfeast_info(fout,fpm,Emid,r,Ntotal,M0)

        end if


        if (com) then
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST runs        |'
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           if (fpm(14)/=2) then
              write(fout,'(A)',advance='no') '#It |  #Eig  |         |Trace|           |     Error-Trace          |     Max-Residual'    
           else
              write(fout,'(A)',advance='no') 'Running average for stochastic estimation (1 -> M0)'  
           endif
           write(fout,*)
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'

        end if



!!!!!!!!!!! Some initialization
        fpm(22)=fpm(8) ! full contour necessary
        loop=0


        fpm(23)=min(M0,fpm(26)) ! 'current M0' size (global value)
        if (fpm(14)==2) fpm(23)=min(fpm(32),fpm(26),M0) ! stochastic estimate
        fpm(25)=fpm(23) !! 'current M0' size (by default)
        fpm(24)=1  !! origin first column number of vector q for parallel mat-vec (default)
        fpm(35)=fpm(23) !! 'current M0' size (by default)
        fpm(34)=M0+1  !! origin first column number of vector q for parallel mat-vec (default)
        !----------------------------------------------
#ifdef MPI
        if (fpm(23)/nb_procs2>=1) then ! criteria for parallelism of mat-vec
           fpm(25)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(25)=fpm(25)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(24)=1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
           fpm(35)=fpm(23)/nb_procs2 ! local size of 'current M0'
           if (rank2==nb_procs2-1) fpm(35)=fpm(35)+fpm(23)-(fpm(23)/nb_procs2)*nb_procs2 
           fpm(34)=M0+1+rank2*(fpm(23)/nb_procs2) ! local origin of first column for vector q for parallel mat-vec 
        end if
#endif
        !-----------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        fpm(21)=1 ! prepare reentry
        if (fpm(5)==0) then !!! random vectors
           iseed=(/56,890,3456,2333/)
           iseed(1) = iseed(1) * (1+rank3*13)
           iseed(2) = iseed(2) + rank3*7
           iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
           iseed(4) = iseed(4) + rank3*2




           if (fpm(14)==2) then
              ! copy work into q for multiplication by B matrix, if stochastic estimate  is "on" 
              call ZLARNV(2,iseed,N*fpm(23),q(1,1)) !!<<4?
           else
              call ZLARNV(2,iseed,N*fpm(23),work(1,1)) !!<<4??

              !  allocate(zzwork(50,M0))
              !  call ZLARNV(2,iseed,50*fpm(23),zzwork(1,1))
              !  if (rank3==0) work(1:N,1:M0)=zzwork(1:N,1:M0)
              !  if (rank3==1) work(1:N,1:M0)=zzwork(N+1:50,1:M0)
              !  deallocate(zzwork)


              if (fpm(15)==0) then !! 2 contours for left eigenvectors
                 iseed=(/6,80,456,5421/)
                 iseed(1) = iseed(1) * (1+rank3*13)
                 iseed(2) = iseed(2) + rank3*7
                 iseed(3) = iseed(3) * (1+rank3*3)+rank3*5
                 iseed(4) = iseed(4) + rank3*2
                 call ZLARNV(2,iseed,N*fpm(23),work(1,M0+1))

              end if
           end if
        end if

        if ((fpm(5)==1).or.(fpm(14)==2)) then !!!!!! q is the initial guess
           !----------------------------------------
#ifdef MPI
           work(1:N,1:fpm(23))=ZZERO 
           if (fpm(15)==0) work(1:N,M0+1:M0+fpm(23))=ZZERO ! 2 contours for left eigenvectors
#endif
           !------------------------------------------
           if (fpm(5)==1) then ! Compute residual vector R=AX-BXE 
              !ijob=30 !!A0*q==>work at first (for compatibility reasons with fpm(5)=0)
              fpm(57)=0
              fpm(21)=9
           else ! only stochastic estimate!! B*Q(1,1:fpm(23))=>work(1:N,1:fpm(23))
              ijob=40 !! B*q=>work
              !endif
              call system_clock(t1,tim)
              !         if (fpm(15)==0) fpm(21)=-1  
              return
           end if
        end if

     end IF  ! info=0
!!!!!!!!!!!!!!
  end if    !ijob=-1

  !  if (fpm(21)==-1) then !! we need to initialize the left eigenvectors as well
  !     ijob=41!! B^T*Q(1,M0+1:M0+fpm(23))=>work(1:N,M0+1:M0+fpm(23))
  !     fpm(21)=1 ! prepare reentry
  !     return
  !  end if





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CONTOUR INTEGRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==1) then !! we initialize a new contour integration
     !------------------------------------------------------------------------
#ifdef MPI
     if ((loop>0).or.(fpm(14)==2)) then !! condition on fpm(5) unnecessary        
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        if ((fpm(15)==0) .and. (fpm(23)/nb_procs2>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,work(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     end if
     !---------------------------for q --------------------------------------
     if ((loop>0).and.(fpm(23)/nb_procs2>=1)) then
        if (.not.((loop==1).and.(fpm(5)==1))) then ! do not include first iteration with initial guess since Q is known to all
           if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
           if ((fpm(15)==0) .and. (fpm(23)/nb_procs2>=1)) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        end if
     end if
     !------------------------------------
#endif



!!!!!!!!!!! RESIDUAL INVERSE ITERATION

     !print *,loop,N,M0,fpm(23)
     if ((loop==0).or.((loop>0).and.(rank2>0))) then
        Q(1:N,1:fpm(23))=ZZERO
        if(fpm(15)==0) Q(1:N,M0+1:M0+fpm(23))=ZZERO
     end if



!!!!!!!!!!!!!!for residual inverse iteration
     if ((loop>0).and.(rank2==0)) then
        allocate(zwork_loc(fpm(23)))
        call zfeast_grationalx(Zne,Wne,fpm(8),lambda,fpm(23),zwork_loc)

        do i=1,fpm(23)
           call ZSCAL(N,zwork_loc(i), q(1,i), 1)
           if (fpm(15)==0) call ZSCAL(N,conjg(zwork_loc(i)), q(1,M0+i), 1)
        end do
        deallocate(zwork_loc)

     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     fpm(20)=1
     fpm(21)=2
     ijob=-2 ! just initialization 
  end IF


  !print *,sum(abs(work(1:N,1))),sum(abs(work(1:N,M0+1)))
  !info=-99
  !ijob=0
  !return
  !stop 


!!!!!!!!!!!!
  IF (fpm(21)==2) then !! we start or pursue the contour integration

     !   IF (info==0) then !! will end up checking info errors returned by FEAST drivers
     do e=fpm(20)+rank2,fpm(22),nb_procs2 !!!! loop over the full contour

        if (ijob==-2) then !!Factorize the linear system

           ! set up factorization id (if fact. storage is used)
           if (fpm(10)==1) then
              fpm(33)=fpm(33)+1
           else
              fpm(33)=1
           endif

           Ze = Zne(e)   
           fpm(20)=e-rank2

           ijob=10 ! for fact

           fact=.false. ! action to factorize
           if ((loop==0).or.((loop==1).and.(fpm(5)==1))) then
              fact=.true. ! always factorize in initial loop. Remark: with input initial guess, loop=1 is the initial loop
           else !(not initial loop)
              if ((fpm(22)>nb_procs2).and.(fpm(10)==0)) fact=.true. ! re-factorize if #mpi<#contour points and fact not stored
              !if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) fact=.true. ! if do not save fact
              !if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) fact=.true. ! if do not save fact with input initial guess
           end if
           if (fact) then
              call system_clock(t1,tim)
              return 
           end if
           !    if ((loop==0).or.(fpm(22)>nb_procs2)) then !no refactorization if one linear system per processor
           !       if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) then !no refactorization if stored in memory
           !          if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) then   
           !             call system_clock(t1,tim)
           !             return 
           !          endif
           !       end if
           !    end if
        end if
!!!!!!!!!!!!!!!! 1st factorization
        if (ijob==10) then !!Solve the linear system (complex) (zS-A)qr=workc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           call ZLACPY( 'F', N, fpm(23), work, N, workc, N )
           ijob=11 ! for solve
           call system_clock(t1,tim)
           return
        end if

!!!!!!!!!!!!!!!! solve
        if (ijob==11) then 
!!!!!!!!!!!!!!!!
           !! summation
           aux = Wne(e) 
           !----------------
           if (loop>0) then ! for residual inverse iteration
              do i=1,fpm(23)
                 call ZSCAL(N,ZONE/(ze-lambda(i)), workc(1,i), 1)
              enddo
           endif
           !---------------

           call ZAXPY(N*fpm(23),aux,workc,1,Q(1,1),1)

           if ((fpm(15)==1).or.(fpm(15)==2)) then
              ijob=-2 !no solve for ql if complex symmetric or 1 single contour
           elseif (fpm(15)==0) then 
              ijob=20 ! for fact

              fact=.false. ! action to factorize
              if ((loop==0).or.((loop==1).and.(fpm(5)==1))) then
                 fact=.true. ! always factorize in initial loop. Remark: with input initial guess, loop=1 is the initial loop
              else !(not initial loop)
                 if ((fpm(22)>nb_procs2).and.(fpm(10)==0)) fact=.true. ! re-factorize if #mpi<#contour points and fact not stored
                 !if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) fact=.true. ! if do not save fact
                 !if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) fact=.true. ! if do not save fact with input initial guess
              end if

              if (fact) then
                 call system_clock(t1,tim)
                 return 
              end if


              !              if ((loop==0).or.(fpm(22)>nb_procs2)) then  !no refactorization if one linear system per processor
              !                 if (.not.((fpm(10)==1).and.(fpm(5)==0).and.(loop>=1))) then !no refactorization if stored in memory
              !                    if (.not.((fpm(10)==1).and.(fpm(5)==1).and.(loop>=2))) then 
              !                       call system_clock(t1,tim)
              !                       return 
              !                    end if
              !                 end if
           end if

        end if     !ijob==11

!!!!!!!!!!!! 2nd factorization
        if (ijob==20) then !!Solve the linear system (complex) (zS-A)^Hql=workc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
           call ZLACPY( 'F', N, fpm(23),work(1,M0+1) , N, workc, N )
           call system_clock(t1,tim)
           ijob=21 ! for solve
           return 
        end if

!!!!!!!!!!!! solve with transpose conjugate
        if (ijob==21) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
           aux = conjg(Wne(e)) 
           if (loop>0) then ! for residual inverse iteration
              do i=1,fpm(23)
                 call ZSCAL(N,ZONE/conjg(ze-lambda(i)), workc(1,i), 1)
              enddo
           endif
           !---------------
           call ZAXPY(N*fpm(23),aux,workc,1,Q(1,M0+1),1)

           ijob=-2 ! default identification
        end if
     end do

     !     end IF !! info=0

     !------------------------------------------------
     !#ifdef MPI
     !     call MPI_BCAST(info,1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
     !#endif
     !-----------------------------------------------  
     !    if (info/=0) fpm(21)=100 ! the end
     if (info==0) then
        fpm(21)=4 
        !------------------------------------------------
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        if(fpm(15)==0) call MPI_ALLREDUCE(MPI_IN_PLACE,Q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
#endif
        !-----------------------------------------------  
     end if

  end IF    ! fpm(21)==1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  print *,'here',rank3,sum(abs(q(1:N,1)))!,sum(abs(q(1:N,M0)))
  !info=-99
  !ijob=0
  !return
  ! stop 


  if ((fpm(21)==4).and.(fpm(14)==1)) then !! only qr,ql vectors has been computed and is returned
     info=4
     if (info/=0) fpm(21)=100 ! The End
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!! Stochastic Estimate
!!$  if ((fpm(21)==4).and.(fpm(14)==2)) then !! only qr is returned with stochastic estimation in res
!!$     iseed=(/56,890,3456,2333/)
!!$     call ZLARNV(4,iseed,N*fpm(23),work(1,1)) 
!!$     call ZGEMM('C','N',fpm(23),fpm(23),N,-ZONE,work,N,q,N,ZZERO,zAq,M0) ! projection
!!$     call ZGEMM('C','N',fpm(23),fpm(23),N,ZONE,work,N,work,N,ZZERO,zBq,M0) ! normalization
!!$     theta=DZERO
!!$     do i=1,fpm(23)
!!$        theta=theta+dble(zAq(i,i)/zBq(i,i))
!!$        res(i)=abs(theta*(N/(DONE*i)))
!!$     enddo
!!$     mode=int(real(res(fpm(23))))+1
!!$     info=5
!!$     if (info/=0) fpm(21)=100 ! The End
!!$  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! USE QR algo to ORTHONORMALIZE Qr and Ql (if any)!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==4) then 

!!!! orthogonalize q
!!!! Modified GS in-place to get QR factorization
     allocate(Rmat(fpm(23),fpm(23))) !! use for R matrix
#ifdef MPI
     if (fpm(49)/=0) then ! L3 parallel
        call pzqrgs(N,fpm(23),Q,Rmat,L3_COMM_WORLD)
        if (fpm(15)==0) call pzqrgs(N,fpm(23),Q(1,M0+1),Rmat,L3_COMM_WORLD) 
     else
        call zqrgs(N,fpm(23),Q,Rmat)
        if (fpm(15)==0) call zqrgs(N,fpm(23),Q(1,M0+1),Rmat)
     end if
#else
     call zqrgs(N,fpm(23),Q,Rmat)
     if (fpm(15)==0) call zqrgs(N,fpm(23),Q(1,M0+1),Rmat) 
#endif
     deallocate(Rmat)



     !  allocate(bb(1:fpm(23),1:fpm(23)))
     !  call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,Q(1,fpm(24)),N,ZZERO,bb(1,fpm(24)),fpm(23))
     !  do i=1,fpm(23)
     !print *,i,bb(i,i),(sum(abs(bb(i,1:fpm(23))))-abs(bb(i,i)))/fpm(23)
     !     enddo
     !deallocate(bb)

!!$
!!$
!!!!! DIRECT QR orthogonalization <--- difficult to do with MPI
!!$     allocate(tau(M0))
!!$     lwork_loc=M0
!!$     allocate(zwork_loc(lwork_loc))
!!$     call ZGEQRF( N, fpm(23), Q, N, TAU, zWORK_loc, LWORK_loc, INFO2 )
!!$     call ZUNGQR( N, fpm(23), fpm(23), Q, N, TAU, zWORK_loc, LWORK_loc, INFO2 )
!!$
!!$ if (fpm(15)==0) then    
!!$ call ZGEQRF( N, fpm(23), Q(1,M0+1), N, TAU, zWORK_loc, LWORK_loc, INFO2 )
!!$ call ZUNGQR( N, fpm(23), fpm(23), Q(1,M0+1), N, TAU, zWORK_loc, LWORK_loc, INFO2 )
!!$ endif     
!!$     deallocate(tau)
!!$     deallocate(zwork_loc)

!!!!!!!!!!!!! Alternative to direct QR, steps:
!!$   !1- form b=Q^H.Q, 2- solver eigenvalue Delta=x^H.b.x, 3- form new Q=Q.x.Delta^(-1/2)
!!$!!!! Attention==>Does not work some Delta values quickly approaching zero...
!!$
!!$
!!$!!! do it for the right eigenvectors   
!!$   allocate(bb(1:fpm(23),1:fpm(23)))
!!$   call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,Q(1,fpm(24)),N,ZZERO,bb(1,fpm(24)),fpm(23))
!!$ !---------------------------------------
!!$#ifdef MPI
!!$     if (fpm(23)/nb_procs2>=1) then
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,bb(1,1),fpm(23)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
!!$     end if
!!$     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,bb(1,1),fpm(23)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
!!$#endif
!!$!---------------------------------------
!!$
!!$     !! Rq: zBq should contain real values, since complex, let us
!!$     !! consider the Hermitian eigenvalue solver
!!$     allocate(dtau(fpm(23)))
!!$     LWORK_LOC=2*fpm(23)-1
!!$     allocate(zWORK_LOC(LWORK_LOC))
!!$     allocate(WORK_LOCr(3*fpm(23)-2))
!!$     call ZHEEV('V','L',fpm(23),bb,fpm(23), dtau, zWORK_loc, LWORK_loc, WORK_locr, INFO_lap)
!!$     !print *,'eig',info_lap,dtau(1:fpm(23))
!!$deallocate(zwork_loc)
!!$deallocate(work_locr)
!!$
!!$ call ZLACPY('F',N,fpm(23),Q(1,1),N,workc(1,1),N)
!!$! Q(1:N,1:fpm(23))=ZZERO
!!$  call ZGEMM('N','N',N,fpm(23),fpm(23),ZONE,workc(1,1),N,bb(1,1),fpm(23),ZZERO,Q(1,1),N)
!!$do i=1,fpm(23)
!!$  call ZSCAL(N,ZONE/sqrt(dtau(i)),Q(1,i),1)
!!$   enddo
!!$deallocate(dtau)
!!$
!!$deallocate(bb)


  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Form the reduced eigenvalue problem
!!!!!!! original problem ap.e^(p)+ap-1.e^(p-1)+...+a1.e+a0=0
!!!!!!! Companion problem
!!!!!!!   
!!!!!!!    [-ap-1 -ap-2 -a0  ]       [ap 0 0]
!!!!!!! a =[I        0   0   ]    b =[0  I 0]
!!!!!!!    [0        I   0   ]       [0  0 I]
!!!!!!!
!!!!!!! We consider ai=Q^H Ai Q (or Q^T or Q^H_left)
!!!!!!! In this example dmax=3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!form  ap=Q^H Ap Q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  if (fpm(21)==4) then 

     !------------------------------------------------
     !#ifdef MPI
     !     work(1:N,1:fpm(23))=ZZERO
     !#endif
     !-------------------------------------------------
     fpm(21)=5 ! preparing reenty
     fpm(57)=dmax+1 ! Ap at first
     ijob=30
     call system_clock(t1,tim)
     return! mat-vec Ap.qr => workr
  end  if



  if (fpm(21)==5) then

     ! after first loop zAq(1,1) contains old residual-
     if (loop>0) Ze=zAq(1,1) !  dummy variable
     !-------------------------------------------------
     zBq(1:M0*dmax,1:fpm(23)*dmax)=ZZERO !common to all processors
     do i=fpm(23)+1,dmax*fpm(23)
        zBq(i,i)=ZONE
     enddo

     !-------------------------------------------------
     if(fpm(15)==1)then 
        call ZGEMM('C','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0*dmax)
     elseif(fpm(15)==2)then 
        call ZGEMM('T','N',fpm(23),fpm(25),N,ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0*dmax)
     else ! default fpm(15)=0
        call ZGEMM('C','N',fpm(23),fpm(35),N,ZONE,Q(1,M0+1),N,work(1,fpm(24)),N,ZZERO,zBq(1,fpm(24)),M0*dmax)
     endif


     !-------!(zBq known to all processors)-
#ifdef MPI
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*dmax*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     end if
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,zBq(1,1),M0*dmax*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
#endif
     !------------------------------------     


     !    print *,'here',sum(abs(zBq(1:M0*dmax,1))),zBq(M0,1)!,sum(abs(q(1:N,M0)))
     !    if (rank3==0) then
     !    do i=1,M0*dmax
     !       print *,i,zBq(i,1)
     !    enddo
     !    endif
     ! info=-99
     ! ijob=0
     ! return   

     zAq(1:M0*dmax,1:fpm(23)*dmax)=ZZERO ! common to all processors


     fpm(21)=6 ! reentry
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! form  ak=Q^H Ak Q for k=p-1 to 0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  if (fpm(21)==6) then

     !------------------------------------------------
#ifdef MPI 
     if (fpm(57)==2) work(1:N,1:fpm(23))=ZZERO !for * with A0
#endif
     !------------------------------------------------
     fpm(21)=7 ! preparing reentry
     fpm(57)=fpm(57)-1 !Ap-1 to A0
     ijob=30
     call system_clock(t1,tim)
     return  ! mat-vec A.qr => workr
  endif



  if (fpm(21)==7) then 

     if (fpm(15)==1)then !notice the -ZONE
        call ZGEMM('C','N',fpm(23),fpm(25),N,-ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zAq(1,(dmax-fpm(57))*fpm(23)+fpm(24)),dmax*M0)
     elseif (fpm(15)==2)then
        call ZGEMM('T','N',fpm(23),fpm(25),N,-ZONE,Q(1,1),N,work(1,fpm(24)),N,ZZERO,zAq(1,(dmax-fpm(57))*fpm(23)+fpm(24)),dmax*M0)
     else ! default fpm(15)=0
        call ZGEMM('C','N',fpm(23),fpm(35),N,-ZONE,Q(1,M0+1),N,work(1,fpm(24)),N,ZZERO,zAq(1,(dmax-fpm(57))*fpm(23)+fpm(24)),dmax*M0)
     endif


     !----------------!(zAq known to all processors)
     !! a bit tricky because of summation
#ifdef MPI 
     if (fpm(23)/nb_procs2>=1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,(dmax-fpm(57))*fpm(23)+1),M0*dmax*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code) 
     end if
     if (fpm(49)/=0) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,zAq(1,(dmax-fpm(57))*fpm(23)+1),M0*dmax*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L3_COMM_WORLD,code)
     end if
#endif
     !---------------------------------------

!!! work on reentry (loop if needed)
     if (fpm(57)>1) then
        fpm(21)=6
        ijob=-2
        return ! perform a loop to fpm(21)=6
     endif

!!!!! add common identity part to zAq
     do k=2,dmax
        do i=1,fpm(23)
           zAq((k-1)*fpm(23)+i,(k-2)*fpm(23)+i)=ZONE
        enddo
     enddo





     if  (fpm(13)==1) then ! customize eigenvalue solver
        fpm(21)=9 ! preparing reentry - could return new value of M0 in fpm(23) if reduced subspace is needed
        ijob=60
        return
     endif
     fpm(21)=81
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Solve the reduced eigenvalue problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==81) then



     ! if using FEAST-MPI ==> solve on a single proc
     if ((rank2==0).and.(rank3==0)) then

!!!!!!!!!!!!!!!!!!!  FEAST

!!$        M02=fpm(23)
!!$        call wallocate_2z(VRL,fpm(23)*dmax,M02*2,infoloc)
!!$        call wallocate_1d(res2,M02*2,infoloc)
!!$        
!!$call feastinit(fpm2)
!!$!fpm2(1)=1
!!$!fpm2(4)=10
!!$call zfeast_gegv(dmax*fpm(23),zAq,M0*dmax,zBq,M0*dmax,fpm2,epsout2,loop2,Emid,r,M02,lambda,VRL,M2,res2,info2)

!!!!!!!!!!!!! LAPACK companion problem
        LWORK_LOC=8*fpm(23)*dmax !! for lapack eig reduced system
        allocate(zWORK_LOC(LWORK_LOC))
        allocate(WORK_LOCR(8*fpm(23)*dmax))
        allocate(LALPHA(fpm(23)*dmax))
        allocate(LBETA(fpm(23)*dmax))
        allocate(VL(dmax*fpm(23),dmax*fpm(23)))
        allocate(VR(dmax*fpm(23),dmax*fpm(23)))

        ! print *,'aq',ZAq(1,1)
        !  print *,'bq',ZBq(1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        JOBVL='V'
        if ((fpm(15)==1).or.(fpm(15)==2)) JOBVL='N' 
        call ZGGEV(JOBVL,'V',fpm(23)*dmax,zAq,M0*dmax,zBq,M0*dmax,LALPHA,LBETA,VL,dmax*fpm(23),VR,dmax*fpm(23), zWORK_LOC, LWORK_LOC, WORK_LOCR,INFO_lap)

        if (info_lap/=0) then
           info=-3 ! problem with lapack
           !print *,'info_lap',info_lap,fpm(23)*dmax
        end if

        !do i=1,fpm(23)
        !          call ZCOPY(fpm(23),VR(fpm(23)*(dmax-1)+1,i), 1, zAq(1,i), 1 )  
        !          if (fpm(15)==0) call ZCOPY(fpm(23),VL(1,i), 1, zBq(1,i), 1 )
        !       enddo


        IF (info_lap==0) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!!! Sort the Eigenvalues of companion problems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
           allocate(lambda2(1:fpm(23)*dmax))
           lambda2(1:fpm(23)*dmax)=lalpha(1:dmax*fpm(23))/lbeta(1:dmax*fpm(23))

           allocate(color(fpm(23)*dmax))
           if (fpm(29)==0) then ! expert custom contour
              call zfeast_inside_contourx(Zne,fpm(8),lambda2,fpm(23)*dmax,color)
           else
              call zfeast_inside_contour(Emid,r,fpm(18),fpm(19),lambda2,fpm(23)*dmax,color)
           endif


           allocate(map(fpm(23)*dmax)) !! for mapping color
           itemp=sum(color)
           k1=0
           k2=itemp
           do i=1,dmax*fpm(23)
              work_locr(i)=abs(lambda2(i)-Emid) ! from center of mass !!!<<< need to do the color first (eigenvalues inside arbitrary contour)
              if (color(i)==1) then
                 k1=k1+1
                 map(k1)=i
              else
                 k2=k2+1
                 map(k2)=i
              endif
           enddo


           allocate(map2(fpm(23)*dmax)) !! for mapping amplitude
           map2(1:fpm(23)*dmax)=map(1:fpm(23)*dmax)
           !! sort (bubble) from itemp+1 to fpm(23)*dmax
           do i=dmax*fpm(23),itemp+2,-1
              !print *,i,itemp
              do j=itemp+1,i-1
                 if (work_locr(map(j))>work_locr(map(j+1))) then
                    dtemp=work_locr(map(j))
                    work_locr(map(j))=work_locr(map(j+1))
                    work_locr(map(j+1))=dtemp

                    k2=map2(j)
                    map2(j)=map2(j+1)
                    map2(j+1)=k2

                 endif
              enddo
           enddo



!!! reorder eigenvalues and eigenvectors

           do i=1,fpm(23)
              call ZCOPY(fpm(23),VR(fpm(23)*(dmax-1)+1,map2(i)), 1, zAq(1,i), 1 )  
              if (fpm(15)==0) call ZCOPY(fpm(23),VL(1,map2(i)), 1, zBq(1,i), 1 )
              lambda(i)=lambda2(map2(i))
           enddo


           deallocate(map)
           deallocate(map2)
           deallocate(color)              
           deallocate(lambda2)
        end IF

        deallocate(LALPHA)
        deallocate(LBETA)
        deallocate(zwork_loc)
        deallocate(work_locr)
        deallocate(VL)
        deallocate(VR)

     end if !(rank 0)

     !-------------------------------- !(info common to all processors) 
#ifdef MPI
     if (rank3==0) call MPI_BCAST(info,1,MPI_INTEGER,0,L2_COMM_WORLD,code)
     if (fpm(49)/=0) call MPI_BCAST(info,1,MPI_INTEGER,0,L3_COMM_WORLD,code)
#endif
     !--------------------------------


     if (info/=0) fpm(21)=100 ! the end


     if (info==0) then

        !---------------------------------------- !(zAq==> vectors, lambda and fpm(23), known to all processors) 
        !#ifdef MPI        
        !        call MPI_BCAST(zAq(1,1),M0*dmax*(fpm(23)),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)  !!! todo
        !       !        call MPI_BCAST(zAq(1,1,2),M0*( fpm(23) ),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        !        call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_COMPLEX,0,NEW_COMM_WORLD,code)
        !        call MPI_BCAST(fpm(37),1,MPI_INTEGER,0,NEW_COMM_WORLD,code)
        !#endif
        !-----------------------------

#ifdef MPI
        if (rank3==0) then 
           call MPI_BCAST(zAq,M0*dmax*(fpm(23)),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
           call MPI_BCAST(zBq,M0*dmax*(fpm(23)),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_COMPLEX,0,L2_COMM_WORLD,code)
        end if
        if (fpm(49)/=0) then 
           call MPI_BCAST(zAq,M0*dmax*(fpm(23)),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
           call MPI_BCAST(zBq,M0*dmax*(fpm(23)),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
           call MPI_BCAST(lambda,fpm(23),MPI_DOUBLE_COMPLEX,0,L3_COMM_WORLD,code)
        end if
#endif


        fpm(21)=9
     end if

  end if !! fpm(21)=81

  ! print *,zAq(1:fpm(23),1)
  !  print *,zBq(1:fpm(23),1)
  !        stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Ritz vectors Xr=Qr*xr 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==9 ) then
     !! from precious step workr=A0*Qr, before using workr memory space
     !! let us construct first workc=-workr*Aq which is workc=-A0*Xr
     !ifdef MPI
     !    if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,work,N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
     !endif

     !    call ZGEMM('N','N',N,fpm(25),fpm(23),-ZONE,work(1,1),N,zAq(1,fpm(24)),M0*dmax,ZZERO,workc(1,fpm(24)),N) 

     if (.not.((fpm(5)==1).and.(loop==0))) then      
        !! construct X (q needed for outer-loop as well)
        call ZLACPY( 'F', N, fpm(23),Q(1,1) , N, work(1,1), N )


#ifdef MPI 
        Q(1:N,1:fpm(23))=ZZERO
#endif

        call ZGEMM('N','N',N,fpm(25),fpm(23),ZONE,work(1,1),N,zAq(1,fpm(24)),M0*dmax,ZZERO,Q(1,fpm(24)),N)
     end if

     !print *,'sum q',sum(q(1:N,1)),sum(q(1,1:fpm(23)))

!!!!!!!!!!!!!!!
     workc(1:N,fpm(24):fpm(24)+fpm(25)-1)=ZZERO  !! initialize the residual temp
#ifdef MPI
     res(1:fpm(23))=DZERO
#endif  
     fpm(21)=10 ! preparing reenty 
     fpm(57)=0
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residuals ||(A2*Xr*lambda^2)+(A1*Xr*lambda)+A0*Xr||
!!!!! (quadratic example)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==10) then


     !if ((fpm(5)==1).and.(loop==0).and.(fpm(57)==1))  call ZLACPY( 'F', N, fpm(23),work , N, workc, N ) !! * with A0  done at the start
     !     if ((fpm(5)==1).and.(loop==0).and.(fpm(57)==1))  then
     !        workc(1:N,fpm(24):fpm(24)+fpm(25)-1)=ZZERO
     !        call ZAXPY(N*fpm(25),-ZONE,work(1,fpm(24)),1,workc(1,fpm(24)),1)!! initialize the residual temp
     !        !call ZLACPY( 'F', N, fpm(25),work(1,fpm(24)) , N, workc(1,fpm(24)), N )
     !     end if

     fpm(21)=11 ! preparing reentry
     fpm(57)=fpm(57)+1 ! from A1 to Ap (A0 already done)
     ijob=30 
     call system_clock(t1,tim)
     return  ! mat-vec A0*qr => workr
  endif


  if (fpm(21)==11) then

     !---------Compute residual vector 
     !  keep adding  work<=workc+work*E^p (needed also for outer-loop)
     do i=fpm(24),fpm(24)+fpm(25)-1
        !if (fpm(57)==dmax+1)  res(i)=(sum(abs(work(1:N,i))))
        if (fpm(57)==dmax+1)  res(i)=DZNRM2(N,work(1,i),1)**2 ! square because of mpi compatibility sqrt(mpi1)+sqrt(mpi2) /= sqrt(mpi1+mpi2) 
        if (fpm(57)>1)  call ZSCAL(N,lambda(i)**(fpm(57)-1),work(1,i),1)
     enddo

     call ZAXPY(N*fpm(25),-ZONE,work(1,fpm(24)),1,workc(1,fpm(24)),1) 

!!! work on rentry (loop if needed)
     if (fpm(57)<dmax+1) then
        fpm(21)=10
        ijob=-2
        return
     endif

     !reduce first denominator
#ifdef MPI 
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif


     !  save residual for reentry 
     !----------------------------------------
#ifdef MPI
     work(1:N,1:fpm(23))=ZZERO  !! work is also needed for outer-loop if any 
#endif
     !------------------------------------------
     !    work(1:N,1:fpm(23))=workc(1:N,1:fpm(23)) !!<<<<?? parallel
     call ZLACPY( 'F', N, fpm(25),workc(1,fpm(24)) , N, work(1,fpm(24)), N )



     !-------- Residual
     do i=fpm(24),fpm(24)+fpm(25)-1
        !      res(i)=sum(abs(work(1:N,i)))/(abs((abs(Emid)+r)**dmax)*res(i))
        res(i)=(DZNRM2(N,work(1,i),1)**2)/(abs((abs(Emid)+r)**(2*dmax))*res(i))
        !        res(i)=DZNRM2(N,work(1,i),1)/(abs(lambda(i)**dmax)*res(i))
     end do
     !----------------------------------------
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res,fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !-----------------------------------------    
     !if (fpm(15)==2) then
     !     do i=1,M0
     !   print *,i,'resr',res(i)
     !enddo
     !if (loop==2) stop
     !endif
     do i=1,fpm(23)
        res(i)=sqrt(res(i)) ! norm L2
     enddo

     fpm(21)=-9
     !fpm(57)=1 ! go back to A0 (preparing for left eigenvectors)
     if ((fpm(15)==1).or.(fpm(15)==2)) fpm(21)=-12
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Compute Ritz vectors Xl=Ql*xl  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (fpm(21)==-9 ) then

     if (.not.((fpm(5)==1).and.(loop==0))) then ! no first time initial guess
        !if (fpm(57)==1) then ! do it only once
        call ZLACPY( 'F', N, fpm(23),Q(1,M0+1) , N, workc(1,1), N )
#ifdef MPI 
        Q(1:N,M0+1:M0+fpm(23))=ZZERO
#endif
        call ZGEMM('N','N',N,fpm(35),fpm(23),ZONE,workc(1,1),N,zBq(1,fpm(24)),M0*dmax,ZZERO,Q(1,fpm(34)),N) 
        !end if
     end if

!!!!!!!!!!!!!!!!!
     !! workc(1:N,1:fpm(23))=ZZERO  !! initialize for residual temp!!<<help
     workc(1:N,fpm(24):fpm(24)+fpm(25)-1)=ZZERO  !! initialize the residual temp
#ifdef MPI
     res(M0+1:M0+fpm(23))=DZERO
!!!work(1:N,M0+1:M0+fpm(23))=ZZERO  !! work is also needed for outer-loop (final residual)
#endif

     fpm(21)=-10 ! preparing reentry
     fpm(57)=0
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Residuals ||(A2^H*Xl*conjg(lambda)^2)+(A1^H*Xl*conjg(lambda))+A0^H*Xl||
!!!!! (quadratic example)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (fpm(21)==-10) then
     fpm(21)=-11 ! preparing reentry
     fpm(57)=fpm(57)+1 ! from A0^H to Ap^H  
     ijob=31    
     call system_clock(t1,tim)  
     return  ! mat-vec A0^H.ql => workl(1,fpm(34))
  endif


  if (fpm(21)==-11) then
     !  if (fpm(57)==1) print *,'q',q(1:N,M0+1)
     !   print *,'R*q',work(1:N,M0+1)


     !---------Compute residual vector 
     ! keep adding  work<=workc+work*E^p (needed also for outer-loop)
     do i=fpm(34),fpm(34)+fpm(35)-1
        !res(i)=sum(abs(work(1:N,i))) ! placeholder (denominator)
        if (fpm(57)==dmax+1)  res(i)=DZNRM2(N,work(1,i),1)**2 ! square because of mpi compatibility sqrt(mpi1)+sqrt(mpi2) /= sqrt(mpi1+mpi2) 
        if (fpm(57)>1) call ZSCAL(N,conjg(lambda(i-M0))**(fpm(57)-1),work(1,i),1)
     enddo

     !      print *,'B*q-scaled',work(1:N,M0+1),conjg(lambda(1))
     !if (fpm(57)==1)       print *,'workc 0?',workc(1:N,1) 

     call ZAXPY(N*fpm(35),-ZONE,work(1,fpm(34)),1,workc(1,fpm(34)-M0),1)
     !if (fpm(57)==2)
     !      print *,'A*q+B*q-scaled',workc(1:N,1)

!!! work on rentry (loop if needed)
     if (fpm(57)<dmax+1) then
        fpm(21)=-10
        ijob=-2
        return
     endif

     !reduce first denominator
#ifdef MPI 
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
#endif

     !  save residual for reentry
#ifdef MPI
     work(1:N,M0+1:M0+fpm(23))=ZZERO  !! work is also needed for outer-loop (final residual)
#endif       
     !  work(1:N,M0+1:M0+fpm(23))=workc(1:N,1:fpm(23)) !!<<<<?? parallel
     call ZLACPY( 'F', N, fpm(25),workc(1,fpm(24)), N, work(1,fpm(34)), N )


     !-------- Residual
     do i=fpm(34),fpm(34)+fpm(35)-1
        !res(i)=sum(abs(work(1:N,i)))/(abs(abs(Emid)+r)*res(i))
        res(i)=(DZNRM2(N,work(1,i),1)**2)/(abs((abs(Emid)+r)**(2*dmax))*res(i))
        ! res(i)=DZNRM2(N,work(1,i),1)/(abs(lambda(i-M0)**dmax)*res(i))
     end do

     !----------------------------------------

     !print *,work(1:N,M0+1),res(M0+1)
     !stop    

     !do i=1,M0
     !   print *,i,'resrl',res(i)!,res(M0+i)
     !enddo
     !if (loop==2) stop
     !----------------------------------------
#ifdef MPI
     if (fpm(49)/=0) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L3_COMM_WORLD,code)
     if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,res(M0+1),fpm(23),MPI_DOUBLE_PRECISION,MPI_SUM,L2_COMM_WORLD,code)
#endif
     !-----------------------------------------
     do i=1,fpm(23)
        res(M0+i)=sqrt(res(M0+i)) ! norm L2
     enddo

     !print *,'res-r-l',maxval(res(1:20)),maxval(res(M0+1:M0+20))



     fpm(21)=-12
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Compute Trace / count eigenvalues (remove spurious if any)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (fpm(21)==-12) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !---color mapping --> check if eigenvalue inside the contour
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     allocate(color(fpm(23)))

     if (fpm(29)==0) then ! expert custom contour
        call zfeast_inside_contourx(Zne,fpm(8),lambda,fpm(23),color)
     else
        call zfeast_inside_contour(Emid,r,fpm(18),fpm(19),lambda,fpm(23),color)
     endif
!!!!!!!!!!!!
     mode=0
     trace=ZZERO
     theta=DZERO
     theta_min=minval(res(1:fpm(23))) ! will be inside the interval
!!! count how many eigenvalues have converged + trace + max residual
     do i=1,fpm(23)
        if (color(i)==1) then ! inside the search interval
           test=.true.
           if ((loop>2).and.(res(i)>(theta_min*1d5))) test=.false. !! remove likely spurious if detected  - simple approach looking at residual
           if (test) then
              mode=mode+1
              trace=trace+lambda(i)
              if (res(i)>theta) theta=res(i) ! max residual (right residual)
              if (fpm(15)==0) then
                 if (res(M0+i)>theta) theta=res(M0+i) ! max residual (left residual)
              end if
           else
              res(i)=-DONE ! spurious
              if (fpm(15)==0) res(i+M0)=-DONE !spurious
           end if
        end if
     end do


     !  if (mode>20) then
     !     print *,'FOUND spurious'
     !     do i=1,mode
     !        print *,i,lambda(i),res(i)
     !     enddo
     !     stop
     !end if


!!!!!!!!! remove spurious if detected  - simple approach looking at residual


     !print *,'trace',trace,mode
!!!!!!!!! remove spurious if any  ---- (work only with the right subspace) ----
!!$     if ((mode==0).and.(k>0)) mode=k ! wait before looking into mode 
!!$     ! Rq: if all eigenvalues k are spurious...FEAST will not converge
!!$     if (mode<k) then
!!$        do j=1,k-mode
!!$           theta=DZERO
!!$           e=1
!!$           do i=1,fpm(23)
!!$              if (color(i)==1) then ! inside the search interval
!!$                 if (res(i)>theta) then
!!$                    e=i
!!$                    theta=res(i) !max
!!$                 endif
!!$              end if
!!$           enddo
!!$           trace=trace-lambda(e) ! update trace
!!$           res(e)=-DONE !spurious
!!$           if (fpm(15)==0) res(e+M0)=-DONE !spurious
!!$        end do
!!$!!! max residual
!!$        theta=DZERO
!!$        do i=1,fpm(23)
!!$           if (color(i)==1) then ! inside the search interval
!!$              if (res(i)>theta) theta=res(i) ! max residual (right residual)
!!$              if (fpm(15)==0) then
!!$                 if (res(M0+i)>theta) theta=res(M0+i) ! max residual (left residual)
!!$              end if
!!$           end if
!!$        enddo
!!$     end if


     if (loop>2) then! wait third iteration (spurious related+ifeast)
        if (mode==0) info=1  ! no eigenvalue detected in the interval
        if ((mode==M0).and.(mode/=N)) info=3 !size subspace too small
     endif
     if (info/=0) then
        fpm(21)=100 ! The End
        deallocate(color)
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
       ! if ((fpm(6)==0).and.(log10((epsout))<(-fpm(3)))) testconv=.true.
         if ((fpm(6)==0).and.(epsout<10.0d0**(-fpm(3))))  testconv=.true.
     end if

     !! residual convergence
     !if ((fpm(6)/=0).and.(log10(theta)<(-fpm(3)))) testconv=.true.
 if ((fpm(6)/=0).and.(theta<10.0d0**(-fpm(3))))  testconv=.true.
     
     !! Two loops minimum if spurious are found
     !     if ((loop<=1).and.(k>mode)) testconv=.false.


!!!! L2 barrier to guarantee synchronization in printing
#ifdef MPI
 call MPI_BARRIER(L2_COMM_WORLD,code)
#endif 

     if (com) then
        if (fpm(5)==0) then
           write(fout,'(I3)',advance='no') loop
        else ! no integration in first loop
           if (loop==0) then
              write(fout,'(A)',advance='no') 'N/A'
           else
              write(fout,'(I3)',advance='no') loop-1
           end if
        endif
        write(fout,'(3X)',advance='no')
        write(fout,'(I4)',advance='no') mode
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') abs(trace) !!<<<<<<
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') epsout
        write(fout,'(3X)',advance='no')
        write(fout,'(ES25.16)',advance='no') theta
        write(fout,*)
     end if


     !     do i=1,M0
     !print *,i,res(i),res(M0+i)
     !     end do

     !stop

     if (.not.testconv) then
        zAq(1,1)=trace ! dummy.....
         e=loop
        if (fpm(5)==1) e=e-1
        if (e==fpm(4)) then
           info=2 ! FEAST did not converge (#loop reaches maximum)
           testconv=.true. ! return final eigenvector anyway
        endif
     endif

     !! handle a special case (prevent convergence if no eigenvalue found)
     !! for mode=0 FEAST will exit after few iterations (see above)
     if (mode==0) testconv=.false. 

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
        allocate(zwork_loc(fpm(23)))
        call zfeast_grationalx(Zne,Wne,fpm(8),lambda,fpm(23),zwork_loc)
        do i=1,fpm(23)
           res(i)=abs(zwork_loc(i)) !! save them all          
           if (color(i)==1)  res(i)=200d0+res(i) 
        enddo
        deallocate(zwork_loc)
        deallocate(color)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Remark: q=x and work=Ax-eBx is known already, using FEAST-MPI q and work are already distributed   (lef,right vector components are known)     
        fpm(21)=1   ! prepare reentry- reloop (with contour integration)
        !fpm(21)=4 ! reloop (without contour integration) -in this case work=q (actually does not need "work")
        loop=loop+1
        if (fpm(10)==1) fpm(33)=0 ! initialize factorization id
        ijob=-2 ! do nothing
        return  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else    ! final eigenvectors/eigenvalues
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

#ifdef MPI       
        if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        if (fpm(15)==0) then
           if (fpm(23)/nb_procs2>=1) call MPI_ALLREDUCE(MPI_IN_PLACE,q(1,M0+1),N*fpm(23),MPI_DOUBLE_COMPLEX,MPI_SUM,L2_COMM_WORLD,code)
        endif
#endif

        !! reorder (shift) lambda, eigenvector and residual
        call ZCOPY( N*fpm(23),Q(1,1), 1, work(1,1), 1 )
        if(fpm(15)==0) call ZCOPY( N*fpm(23),Q(1,M0+1), 1, work(1,fpm(23)+1), 1 )
        allocate(zwork_loc(fpm(23)))
        if (fpm(15)==0) then
           allocate(work_loc2(2*fpm(23))) 
        else
           allocate(WORK_LOC2(fpm(23))) 
        end if
        call ZCOPY(fpm(23),lambda , 1, zwork_loc, 1 )
        call DCOPY(fpm(23),res , 1, work_loc2, 1 )
        if(fpm(15)==0) call DCOPY(fpm(23), res(M0+1), 1, work_loc2(fpm(23)+1), 1 )

        M0=fpm(23)  ! update value of M0 (new subspace)

        Q(1:N,1:fpm(23))=ZZERO
        if(fpm(15)==0) Q(1:N,M0+1:M0+fpm(23))=ZZERO
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
                 if(fpm(15)==0) then
                    q(1:N,M0+e)=work(1:N,M0+i)
                    res(M0+e)=work_loc2(M0+i)
                 endif
              endif
           else
              j=j+1
              q(1:N,mode+j)=work(1:N,i)
              lambda(mode+j)=zwork_loc(i)
              res(mode+j)=work_loc2(i)
              if(fpm(15)==0) then
                 q(1:N,M0+mode+j)=work(1:N,M0+i)
                 res(M0+mode+j)=work_loc2(M0+i)
              endif
           end if
           if (work_loc2(i)==-DONE) then ! spurious at the end
              k=k+1
              q(1:N,fpm(23)-k+1)=work(1:N,i)
              lambda(fpm(23)-k+1)=zwork_loc(i)
              res(fpm(23)-k+1)=work_loc2(i)
              if(fpm(15)==0) then
                 q(1:N,M0+fpm(23)-k+1)=work(1:N,M0+i)
                 res(M0+fpm(23)-k+1)=work_loc2(M0+i)
              end if
           endif
        end do

        deallocate(zwork_loc)
        deallocate(work_loc2)
        deallocate(color)

!!!!!!!!!!!!!!!!!!!!!!!!!!
        fpm(21)=100 ! The End

        if ((info==0).and.(fpm(37)==1)) info=6

     end if ! test convergence
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (fpm(21)==100) then !! THE END (ijob=0) 


!!!! L2 barrier to guarantee synchronization in printing
!#ifdef MPI
! call MPI_BARRIER(L2_COMM_WORLD,code)
!#endif 
     
     ijob=0 !! exit FEAST

     if (com) then !! Print  Information

        if (info>=200) then
           write(fout,'(A)',advance='no') 'PROBLEM with input parameters'
           write(fout,*) 
        end if

        if (info==-3) then
           write(fout,'(A)',advance='no') 'ERROR with reduced system'
           write(fout,*) 
        end if

!!$          if (info==-2) then
!!$             write(fout,'(A)',advance='no') 'ERROR from Inner Linear System Solver in FEAST driver'  
!!$             write(fout,*) 
!!$          end if
!!$
!!$          if (info==-1) then
!!$             write(fout,'(A)',advance='no') 'ERROR with Internal memory allocation'  
!!$             write(fout,*) 
!!$          end if

        if (info==1) then
           write(fout,'(A)',advance='no') '==>WARNING: No eigenvalue has been found in the proposed search interval'
           write(fout,*)
        endif

        if (info==3) then
           write(fout,'(A)',advance='no') '==>WARNING: Size subspace M0 too small'  
           write(fout,*)
        end if

        if (info==4) then
           write(fout,'(A)',advance='no') '==>WARNING: Only the subspace has been returned'  
           write(fout,*)
        end if

        if (info==5) then
           write(fout,'(A)',advance='no') '==>WARNING: Only stochastic estimation of #eigenvalues returned'  
           write(fout,*)
        end if

        if (info==6) then
           if (fpm(15)==0) then
              write(fout,'(A)',advance='no') '==>WARNING: FEAST converges but left/right subspaces are not bi-orthonormal'
           else
              write(fout,'(A)',advance='no') '==>WARNING: FEAST converges but right subspace is not orthonormal'
           end if
           write(fout,*) 
        end if

        if (info==2) then
           write(fout,'(A)',advance='no') '==>WARNING: FEAST did not converge "yet" (#loop reaches maximum allowed)'  
           write(fout,*)
        end if


        if (info==0) then
            write(fout,*)
            write(fout,'(A)',advance='no') '==>FEAST has successfully converged'
           if (fpm(6)==0) then  
              write(fout,'(A)',advance='no') ' with Trace tolerance <1E-'
           else
               write(fout,'(A)',advance='no') ' with Residual tolerance <1E-'
            end if
             if (fpm(3)<10) then
               write(fout,'(I1)') fpm(3)
            else
               write(fout,'(I2)') fpm(3)
            end if
             write(fout,'(A)',advance='no') '   # FEAST outside it. '
             write(fout,'(I8)') loop
             if (mod(fpm(30),10000)/1000==2) then !! ifeast
                 write(fout,'(A)',advance='no') '   # Inner BiCGstab it.'
                write(fout,'(I8)',advance='no') fpm(60)
                if (fpm(30)/100000==2)    write(fout,'(A)',advance='no') ' for rank2=0'  !pfeast
                 write(fout,*)
              end if
                write(fout,'(A)',advance='no') '   # Eigenvalue found  '
           write(fout,'(I8)') mode
        else
           write(fout,'(A)',advance='no') '==>INFO code = '
           write(fout,'(I4)',advance='no') info
           write(fout,*)
        end if
        

        
        if (info>=0) then             
           call system_clock(tend,tim) !! total time
           write(fout,'(A)') '----------------------------------------------------------------------------------------------'
           write(fout,*)
           write(fout,'(A)') '.-------------------.'
           write(fout,'(A)') '| FEAST-RCI timing  |'
           write(fout,'(A)') '--------------------.------------------.'
           write(fout,'(A)',advance='no') '| Fact. cases(10,20)|'
           write(fout,'(F12.4)',advance='no') (t10*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Solve cases(11,12)|'
           write(fout,'(F12.4)',advance='no') (t11*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| T(E)*x   cases(30)|'
           write(fout,'(F12.4)',advance='no') (t30*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| T(E)^H*x cases(31)|'
           write(fout,'(F12.4)',advance='no') (t31*1.0d0/tim)
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Misc. time        |'
           write(fout,'(F12.4)',advance='no') (tend-tini-t10-t11-t30-t31)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)',advance='no') '| Total time (s)    |'
           write(fout,'(F12.4)',advance='no') (tend-tini)*1.0d0/tim
           write(fout,'(A)') '      |'
           write(fout,'(A)') '--------------------------------------- '
           write(fout,*)
        endif
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,'(A)',advance='no') '*********** FEAST- END*************************'
        write(fout,*) 
        write(fout,'(A)',advance='no') '***********************************************'  
        write(fout,*) 
        write(fout,*)
        !! close file
        if (fpm(1)<0) close(fout)
     endif
#ifdef MPI
     call MPI_BARRIER(L2_COMM_WORLD,code)
     !if (fpm(49)>0)  call MPI_BARRIER(L3_COMM_WORLD,code)
#endif

  end if


end subroutine zfeast_grcipevx



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_grcipev(ijob,dmax,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve polynomial eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      P(lambda)qr=0 and P(lambda)^Hql=0                            
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  ! 
  !  Expert comments: This is also a Kernel routine (contains several cases depending of fpm(15))           
  !
  !                   fpm(15)=0 ! 2 sided algo (default)
  !                   fpm(15)=1 ! 1 sided algo- right eigenvector
  !                   fpm(15)=2 ! 1 sided algo with right/left conjugate- default for complex symmetric
  !
  !  A can be NON-HERMITIAN (could be real or complex)  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,2*M0):   Workspace
  !                            expert comment: size of (N,M0) if fpm(15)>0  
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !                            expert comment: size of (N,M0) if fpm(15)>0 (only right vectors)
  !  
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(2*M0) : Relative Residual of the solution (1-norm) 
  !                                                          contains right (1:mode) and left (M0+1:M0+mode)
  !                          
  !                            expert comment: size of (M0) if fpm(15)>0 (only rigth vectors)
  !  info       (output)       INTEGER: Error handling (0: successful exit) -- see FEAST documentation
  !=========================================================================
  !  Eric Polizzi 2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0,dmax
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0*dmax,*):: zAq,ZBq
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
!!!    Call expert routine zfeast_grcipevx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 


  if (ijob==-1) then
     fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     !!fpm(15)=2 ! ID for 1 sided feast for Complex Symmetric Matrix (default with name)
     if (fpm(30)==-111) then
        fpm(30)=141316 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  call zfeast_grcipevx(ijob,dmax,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne)


end subroutine zfeast_grcipev











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine zfeast_srcipevx(ijob,dmax,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) -  Includes option for custom integration nodes/weight
  !  Solve polynomial eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      P(lambda)qr=0 and P(lambda)^Hql=0                            
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  !
  !  A are COMPLEX SYMMETRIC 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !  Eric Polizzi  2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0,dmax
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0*dmax,*):: zAq,zBq
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
!!!    Call expert routine zfeast_grcipevx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ijob==-1) then
     !!fpm(15)=2 ! ID for 1-sided feast for Complex Symmetric Matrix (default with name)
     if (fpm(30)==-111) then
        fpm(30)=141117 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if

  call zfeast_grcipevx(ijob,dmax,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,Q,mode,res,info,Zne,Wne)

end subroutine zfeast_srcipevx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine zfeast_srcipev(ijob,dmax,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose 
  !  =======
  !  FEAST RCI (Reverse Communication Interfaces) - 
  !  Solve generalized eigenvalue problems and obtain eigenvalues lambda, and right qr 
  !      Aqr=lambda Bqr (Remark ql=conjg(qr))  within a search interval                        
  !
  !  Solve polynomial eigenvalue problems and obtain eigenvalues lambda, right qr and left ql eigenvectors
  !      P(lambda)qr=0 and P(lambda)^Hql=0                            
  !  wih P(E)=E^(dmax)*A_dmax+...+E*A_1+A_0
  ! 
  !  A are COMPLEX SYMMETRIC 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  ijob       (input/output) INTEGER :: ID of the RCI
  !                            INPUT on first entry: ijob=-1 
  !                            OUTPUT Return values (0,10,20,21,30,31,40,41)-- see FEAST documentation
  ! dmax        (input)        INTEGER: The degree of the polynomial
  !  N          (input)        INTEGER: Size system
  !  Ze         (output)       COMPLEX DOUBLE PRECISION        : Integration node
  !  work       (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace
  !  workc      (input/output) COMPLEX DOUBLE PRECISION (N,M0):   Workspace 
  !  zAq,zBq    (input/output) COMPLEX DOUBLE PRECISION (M0,M0) : Workspace for Reduced Eigenvalue System
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
  !  Eric Polizzi - 2019
  !=====================================================================

  implicit none
  integer :: ijob,N,M0,dmax
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  complex(kind=(kind(1.0d0))),dimension(N,*):: work,workc
  complex(kind=(kind(1.0d0))),dimension(M0*dmax,*):: zAq,zBq
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
!!!    Call expert routine zfeast_grcipevx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex(kind=(kind(1.0d0))), dimension(1:fpm(8)) :: Zne,Wne 


  if (ijob==-1) then
     fpm(29)=1 ! flag to inform expert routine that Zne,Wne generated using default contour
     !!fpm(15)=2 ! ID for 1 sided feast for Complex Symmetric Matrix (default with name)
     if (fpm(30)==-111) then
        fpm(30)=141116 ! code name
        call feastdefault(fpm,info)
     end if
     if (info/=0) then
        ijob=0 ! leave rci
        return
     end if
  end if

  call zfeast_gcontour(Emid,r,fpm(8),fpm(16),fpm(18),fpm(19),Zne,Wne)
  call zfeast_grcipevx(ijob,dmax,N,Ze,work,workc,zAq,zBq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info,Zne,Wne)

end subroutine zfeast_srcipev


