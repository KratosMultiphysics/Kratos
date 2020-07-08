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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Utility Routines for FEAST (Documented)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! List of routines:
! --------------------
!
!111111111111!!!!!! FEASTINIT Routines 
!
! feastinit                      ! initialize FEAST parameters fpm
! feastdefault                   ! choose default parameters for fpm (and error handling for user fpm parameters)
! feastinit_driver (undocumented)! feastinit with additional parameters for predefined drivers (if nedeed by user)


!22222222222!!!!!! DEFAULT CONTOUR FEAST Routines
!
! zfeast_contour             ! Return nodes/weights from ellipsoid contour symmetric with real axis- half contour only  
! cfeast_contour             ! Return nodes/weights from ellipsoid contour symmetric with real axis- half contour only  
! zfeast_gcontour            ! Return nodes/weights from ellipsoid contour located in complex plane (general)- return full contour
! cfeast_gcontour            ! Return nodes/weights from ellipsoid contour located in complex plane (general)- return full contour


!33333333333!!!!! CUSTOM CONTOUR FEAST Routines (assist the user in generating weight/modes for custom contour in compelx plane)
!
! zfeast_customcontour
! cfeast_customcontour


!44444444444!!!!! RATIONAL FEAST Routines
!
! dfeast_rational     ! Return the values of rational/selection function for default contour symmetric with real axis (eigenvalue real)
! sfeast_rational
!     
! dfeast_rationalx    ! Return the values of rational/selection function for custom contour symmetric with real axis (eigenvalue real)
! sfeast_rationalx
!
! zfeast_grational    ! Return the values of rational/selection function for default contour located in complex plane (eigenvalue complex)
! cfeast_grational
!
! zfeast_grationalx   ! Return the values of rational/selection function for custom contour (eigenvalue complex) 
! cfeast_grationalx
!


!555555555!!!!! Routines used by PFEAST (MPI needed)

! pfeastinit                :: initialize FEAST parameters fpm, and L2/L3 communicators
! pfeast_distribution_type  :: determine if a matrix is row distributed or not



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!111111111111!!!!!! FEASTINIT Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine feastinit(fpm)
  !  Purpose
  !  =======
  !
  !  Initialize all the input FEAST parameters
  !
  !  Arguments
  !  =========
  !
  !  fpm  (output) INTEGER(*): FEAST parameters (size array at least 64)
  !=====================================================================
  ! Eric Polizzi 2009-2019
  ! ====================================================================
  implicit none
  !-------------------------------------
  integer,dimension(*) :: fpm
  integer :: i
  do i=1,64
     fpm(i)=-111
  enddo
end subroutine feastinit




subroutine feastdefault(fpm,info)
  !  Purpose
  !  =======
  !
  !  Define the default values for the input FEAST parameters.
  !  Modified some fpm Values that are equal to -111 (the ones that are not modified explicitly by users after feastinit call)
  !  These default values are dependent on the id of FEAST scheme stored in fpm(30)
  !  The routine is also checking of the user fpm values are valid or not, and returns the value of info accordingly 
  !
  !  Arguments
  !  =========
  !
  !  fpm  (input/output) INTEGER(*): FEAST parameters (size array at least 64)
  !=====================================================================
  ! Eric Polizzi 2009-2019
  ! ====================================================================
  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  !  include "f90_noruntime_interface.fi"
  integer,dimension(*) :: fpm
  integer :: i,info
  logical :: test
  character(len=3) :: ctemp
  integer(8),save :: fout
  logical,save :: com
  integer :: rank2,rank3,code,NEW_COMM_WORLD
  integer,parameter :: max=5
  integer, dimension(max):: tnbe=(/24,32,40,48,56/)
  integer,dimension(6) :: dig
  integer :: rem

!!!! extract digit from fpm(30) code routine name
  rem = fpm(30)
  DO i = 1, 6
     dig(6-i+1) = rem - (rem/10)*10  ! Take advantage of integer division
     rem = rem/10
  END DO
!!!!


  info=0 ! initialization

!!!!! fpm(1)
!!!!! (0/1) comments off/on and <0 write on file feast<|fpm(1)|>.log
  if (fpm(1)==-111) then
     fpm(1)=0 
  elseif (fpm(1)>1) then
     info=101
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! fpm(14)
!!!! FEAST execution options
  ! 0- normal FEAST execution
  ! 1- return only subspace Q size M0 after 1 contour          
  ! 2- return stochastic estimates of the eigenvalue count
!!!! must be defined here for setting up fpm(2) and fpm(8) 
  if (fpm(14)==-111) then
     fpm(14)=0 ! (0,1,2)
  elseif   ((fpm(14)<0).or.(fpm(14)>2)) then
     info=114
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! fpm(16)
!!!! Integration type  (0: Gauss, 1: Trapezoidal, 2: Zolotarev)
!!!! must be defined here for setting up fpm(2) and fpm(8)
  if (fpm(16)==-111) then
     fpm(16)=0 ! Gauss default for FEAST symmetric
     if (dig(3)==2) fpm(16)=1 ! Trapezoid default for ifeast
     !! choose Trapezoid for non-Hermitian problem (old default value flag fpm(17))
     if (dig(4)==3) fpm(16)=1 ! Trapezoid for non-symmetric eigenvalue
     if ((dig(4)==1).and.(dig(2)==4)) fpm(16)=1 ! Trapezoid for complex symmetric eigenvalue
  elseif   ((fpm(16)<0).or.(fpm(16)>2)) then
     info=116
  end if
  if (fpm(16)==2) then
     if ((dig(4)==3).or.((dig(4)==1).and.(dig(2)==4))) info=116 ! No zolotarev for non-Hermitian problems
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  
!!!! fpm(2)
!!!! # of contour integration nodes (half-contour) - Hermitian 
  if (fpm(2)==-111) then
     fpm(2)=8 ! half-contour- Symmetric/Hermitian problem.     
     if (dig(3)==2) fpm(2)=4 ! ifeast - half nodes
     if (fpm(14)==2) fpm(2)=3 ! stochastic estimate
  elseif (fpm(2)<1) then
     info=102
  elseif ((fpm(16)==0).or.(fpm(16)==2)) then ! Gauss or Zolotarev restrictions
     test=.false.
     if (fpm(2)>20) then
        test=.true.
        do i=1,max
           if (fpm(2)==tnbe(i)) test=.false. 
        enddo
     end if
     if (test) info=102  
  end if

!!!! fpm(3)
!!!! tol double precision  
  if (fpm(3)==-111) then
     fpm(3)=12
  elseif  ((fpm(3)<0).or.(fpm(3)>16)) then
     info=103
  end if

!!!! fpm(4)
!!!! maxloop 
  if (fpm(4)==-111) then
     fpm(4)=20
     if (dig(3)==2) fpm(4)=50 ! ifeast
  elseif (fpm(4)<0) then
     info=104
  end if

!!!! fpm(5)
!!!! (0,1) Initial Guess Subspace/EigenValues no/yes   
  if (fpm(5)==-111) then
     fpm(5)=0
  elseif ((fpm(5)/=0).and.(fpm(5)/=1)) then
     info=105
  end if

!!!! fpm(6)
!!!! Convergence criteria on trace (0) / eigenvectors relative residual (1)
  if (fpm(6)==-111) then
     fpm(6)=1
  elseif  ((fpm(6)/=0).and.(fpm(6)/=1)) then
     info=106
  endif

!!!! fpm(7)
!!!! tol single precision <<<<<<<<<<<<< Deprecated in V4.0
  if (fpm(7)==-111) then
     fpm(7)=5
  elseif ((fpm(7)<0).or.(fpm(7)>7)) then
     info=107
  endif

!!!! fpm(8)
!!!! # of contour points (full-contour) - Non-Hermitian 
  if (fpm(8)==-111) then
     fpm(8)=16
     if (dig(3)==2) fpm(8)=8 ! ifeast - half nodes
     if (fpm(14)==2) fpm(8)=6 ! stochastic estimate
  elseif (fpm(8)<2) then
     info=108
  elseif (fpm(16)==0) then ! Gauss  for Non-Hermitian (#contour points - full-contour==2 half-Gauss contour)
     test=.false.
     if (fpm(8)>40) then
        test=.true.
        do i=1,max
           if (fpm(8)==2*tnbe(i)) test=.false. 
        enddo
     end if
     if (test) info=108
  end if

!!!! fpm(9)
!!!! default value !L2 communicator 
  !---------------------------------------
#ifdef MPI
  if (fpm(9)==-111) fpm(9)=MPI_COMM_WORLD 
#endif
!----------------------------------------

!!!! fpm(10)
!!!! FEAST used with stored factorization (0: no, 1: yes-store all the factorizations)  
  if (fpm(10)==-111) then
     fpm(10)=1 
     if (dig(5)==1) fpm(10)=0 ! direct call to RCI routines
     !if (dig(1)==2) fpm(10)=0 ! pfeast default
  elseif  ((fpm(10)/=0).and.(fpm(10)/=1)) then
     info=110
  end if

  !fpm(11)
  !fpm(12)

!!!! fpm(13)
!!!! customize rci interface (expert users)
  !0: default (options 10,11,20,21,30,31,40,41)
  !1: reduced eigenvalue (option 60)
  !2: inner product (option 50)
  !3: inner product+reduced eigenvalue (options 50,60)
  if (fpm(13)==-111) then
     fpm(13)=0
  elseif ((fpm(13)<0).or.(fpm(13)>3)) then
     info=113
  end if

!!!! fpm(14)
!!!! FEAST execution options
  ! 0- normal FEAST execution
  ! 1- return only subspace Q size M0 after 1 contour          
  ! 2- return stochastic estimates of the eigenvalue count
  
!!!!>>>>>>>>>>>> must be defined above for setting up fpm(2) and fpm(8)
 
!!!! fpm(15)
!!!! Contour schemes for FEAST for non-Hermitian. Choose between:
  !0- two-sided contour (default for complex /real non-symmetric)
  !1- 1-sided standard
  !2- 1-sided complex symmetric (default for complex symmetric)        
  if (fpm(15)==-111) then
     fpm(15)=0 ! 2 contour default non-symmetric or non-linear 
   !  if (dig(4)==3) fpm(15)=0 ! 2 contours (non-symmetric, linear or not)
   !  if ((dig(4)==2).and.(dig(6)>=6)) fpm(15)=0 ! 2 contours (polynomial Hermitian)
     if (dig(4)==1) fpm(15)=2 ! 1 contour left=right*
     
  elseif   ((fpm(15)<0).or.(fpm(15)>2)) then
     info=115
  end if
  if (fpm(14)==2) fpm(15)=1 ! Compute only right contour if estimate is on

!!!! fpm(16)
!!!! Integration type  (0: Gauss, 1: Trapezoidal, 2: Zolotarev)

!!!!>>>>>>>>>>>> must be defined above for setting up fpm(2) and fpm(8)
  
  
!!!! fpm(17) ! deprecated in v4.0
!!!! Integration type for non-symmetric (0: Gauss, 1: Trapezoidal)       
  !if (fpm(17)==-111) then
  !   fpm(17)=1
  !elseif   ((fpm(17)<0).or.(fpm(17)>1)) then
  !   info=117
  !end if

!!!! fpm(18)
!!!! ellipsoid contour - for Symmetric - fpm(18)/100 is ratio a/b (b is [Emin-Emax]) - Rq: circle is fpm(18)=100        
  if (fpm(18)==-111) then
     fpm(18)=100
     if ((dig(3)==1).and.(dig(6)<=5)) then ! feast and linear eig.
        if (dig(4)==2) fpm(18)=30 ! hermitian
        if ((dig(4)==1).and.(dig(2)/=3).and.(dig(2)/=4)) fpm(18)=30 ! real symmetric
     endif
  !if ((dig(3)==2).and.(dig(6)<=5).and.(fpm(14)==2)) then ! ifeast and linear eig. and stochastic
  !      if (dig(4)==2) fpm(18)=30 ! hermitian
  !      if ((dig(4)==1).and.(dig(2)/=3).and.(dig(2)/=4)) fpm(18)=30 ! real symmetric
  !   endif     
  elseif  (fpm(18)<0) then
     info=118
  end if

!!!! fpm(19)
!!!! Rotation angle in degree for ellispoid contour - for Non-symmetric [-180,180]        
  if (fpm(19)==-111) then
     fpm(19)=0
  elseif  ((fpm(19)<-180).or.(fpm(19)>180)) then
     info=119
  end if


!!!!!!!!! from 20 to 29
  !fpm(20)  ! current node/shift number (local with L2)
  !fpm(21)  ! internal RCI probe
  !fpm(22)  ! number of linear systems to solve
  !fpm(23)  ! effective number of rhs to solve (<=M0)  
  !fpm(24)  ! local origin of the first column (rhs submatrix) using L2 (case 30/31-40/41)
  !fpm(25)  ! local size of rhs block using L2 (case 30/31-40/41)
  !fpm(26)  ! contains calculated Ntotal from distributed L3 systems 
  !fpm(27)  
  !fpm(28)
  if (fpm(29)==-111) fpm(29)=0  ! 0/1 (internal) custom contour type (1:generated using default contour, 0 otherwise-default)  

!!!!!!!!! from 30 to 39

  !fpm(30) !(internal) name of feast routine
  fpm(31)=40 !(internal) FEAST version fpm(31)/10
  if (fpm(32)==-111) fpm(32)=10 ! Stochastic estimate- number of steps (trials) 
  !fpm(33) !id of the factorization (if fact. saving used)-works with fpm(10)
  !fpm(34) !same idea that fpm(24) in the case of left eigevectors stored in [M0+1:2*M0]
  !fpm(35) !same idea that fpm(25) in the case of left eigevectors stored in [M0+1:2*M0]
  if (fpm(36)==-111) fpm(36)=1 !(0:No,1:Yes) bi-orthogonalization flag (Yes/No) (non-symmetric feast) --undocumented --  
  if (fpm(37)==-111) fpm(37)=0 ! (internal)- used for bi-orthogonality test 
  if (fpm(38)==-111) fpm(38)=1 ! (0:No,1:Yes) spurious detection (Yes,No) (non-symmetric feast) --undocumented --
  if (fpm(39)==-111) fpm(39)=0 !fpm(39) (0:default,1)- solve the standard/1 vs generalized/0 reduced system (non-symmetric feast) -- undocumented --

!!!!!!!!!!!  
!!!!!!!!!!! from 40 to 49       !!!!!!!!!!!! DRIVERS OPTIONS
!!!!!!!!!!!  
  if (fpm(40)==-111) then
     fpm(40)=0
     ! Search interval option for Real-Symmetric/ Hermitian drivers
     !0  - undefined- normal FEAST use with user defined Emin,Emax
     !1  - search M0/2 largest eigenvalues- FEAST will find Emin,Emax
     !-1 - search M0/2 smallest eigenvalues- FEAST will find Emin,Emax
  endif

  if (fpm(41)==-111) then
     fpm(41)=1 !scale matrix 0:No 1:Yes (default 1) --> useful for single precision version, mixed precision, or cluster pardiso bug (see system1)--- works only for sparse drivers 
  end if
     
  if (fpm(42)==-111) fpm(42)=1 ! mixed precision/double precision
  !0- (off) full double precision
  !1- (on) mixed precision (single precision solver)- (new) default
  if (fpm(43)==-111) fpm(43)=0 !switch feast_sparse to ifeast_sparse interfaces (use of ifeast using old interfaces)
  !0- feast (default)
  !1- ifeast
  if (fpm(44)==-111) fpm(44)=0 ! Type of iterative solvers for IFEAST
  !0- bicgstab  no preconditioner (default)
  if (fpm(45)==-111) fpm(45)=1 !Accuracy of iterative solver 10^(-fpm(45))
  !2- default
  if (fpm(46)==-111) fpm(46)=40 !Maximum of iteration for inner solver (100 default)
  if (fpm(47)==-111) fpm(47)=0 !Enforced load balancing (# inner it by contour points)
  !0- Off-  Each linear systems must converge to fpm(45) accuracy with fpm(46) max iterations
  !1- On- Fastest linear system to converge enforces max iterations for all (default)
  if (fpm(48)==-111) fpm(48)=0 !Vector convergence criteria
  !1- convergence sets on "max residual"...last rhs to cv
  !0- convergence sets on "max residual" but within search interval
  !-1- convergence sets on "min residual"...first rhs to cv

  if (fpm(49)==-111) then
     fpm(49)=0 ! L3 communicator (0:default)
  end if
  ! special case- calling pfeasinit follows by a feast routine 
#ifdef MPI 
  if ((fpm(49)/=0).and.(dig(5)/=1)) then ! it means that L3 has been initialized (using pfeastinit for example) 
                                         ! AND we are not using the RCI interfaces
     if ((dig(1)==1).and.(fpm(43)==0)) then! if no pfeast routines are used, there should be no L3, and everthing goes to L2
           fpm(49)=0
           fpm(9)=fpm(59) ! this should be the communicator that regroups fpm(49) and fpm(9) (NEW_COMM_WORLD) used in pfeastinit
        endif
     end if
#endif
     
!!!!!!!!!!!!!! from 50 to 59
!fpm(50)=0,1 ! type of partitioning (drivers only)---< does not seem to be used anymore
  !fpm(51) ! local origin of the first matrix/rhs row using L3 (used only for initial guess- for reproducibility)
  !fpm(52)
  !fpm(53) ! local size of matrix/rhs row block using L3 (used only for initial guess - for reproducibility) 
  !fpm(54)
  !fpm(55)
  !fpm(56)
  !fpm(57) ! non-linear feast polynomial (degree order process for RCI) 
  !fpm(58)
  !fpm(59) ! contains global L1 communicator (NEW_COMM_WORLD) defined by user in  pfeastinit

  !!!! Possible output fpm   
  if (fpm(60)==-111) fpm(60)=0! This is an "output" fpm that counts the number of BicGstab iterations 
  !fpm(61)
  !fpm(62)
  !fpm(63) 
  if (fpm(64)==-111) fpm(64)=0 ! Additional feast parameters for driver interfaces (i.e size fpm>64) (0,1)  --undocumented--



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Error Handling !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! comments?
     com=.false. !default
     if (fpm(1)/=0) com=.true.
     rank2=0
     rank3=0
#ifdef MPI
     call MPI_COMM_RANK(fpm(9),rank2,code)
     if (fpm(49)/=0) call MPI_COMM_RANK(fpm(49),rank3,code)
#endif
     if (com.and.((rank2/=0).or.(rank3/=0))) com=.false. ! comment only in rank 0 if any

     if (com) then !! print information
        !file or screen output?
        if (fpm(1)<0) then
           write(ctemp,'(I3)') abs(fpm(1))
           fout=abs(fpm(1))+200!fpm(60) ! file number default
           open(fout,file="feast"//trim(adjustl(ctemp))//".log",action='write',position='append')
           write(fout,'(A)',advance='no') '############## NEW FEAST RUN ###################'
           write(fout,*)  
        elseif (fpm(1)>0) then
           fout=6 !screen
        endif
     end if

     if ((com).and.(info/=0)) then
        write(fout,'(A)',advance='no') 'PROBLEM with FEAST input parameter fpm('
        write(fout,'(I2)',advance='no') info-100
        write(fout,'(A)',advance='no') ')'
        write(fout,*)
     end if

if ((com).and.(fpm(1)<0)) close(fout) ! close file
     
  

end subroutine feastdefault




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine feastinit_driver(fpm,N)
  !  Purpose
  !  =======
  !
  !  Define the default values for the input FEAST parameters from 1-64
  !  and 65-N are optional user defined inputs that can be used for FEAST predefined interfaces,
  !  here fpm(65:N) is initialized at -111.
  !
  !  Arguments
  !  =========
  !
  !  fpm  (output) INTEGER(*): FEAST parameters (size array at least 64)
  !  N    (input)  INTEGER: size of array fpm (>=64)
  !=====================================================================
  ! Eric Polizzi 2009-2019
  ! ====================================================================
  implicit none
  integer,dimension(*) :: fpm
  integer :: N
  integer :: i

  call feastinit(fpm)

  if (N>64) then
     fpm(64)=1
     do i=65,N
        fpm(i)=-111 ! default values for the drivers inputs
     enddo
  end if
end subroutine feastinit_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!22222222222!!!!!! DEFAULT CONTOUR FEAST Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_contour(Emin,Emax,fpm2,fpm16,fpm18,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN FEAST INTEGRATION NODES/WEIGHTS
  !                         FROM an INPUT ELLIPSOID CONTOUR SYMMETRIC with REAL AXIS- return half-contour only
  !  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emin,Emax  (input)   REAL DOUBLE PRECISION    :search interval (horizontal ellipsoid axe b)
  !  fpm2       (input)   INTEGER                  :# nodes in half-contour 
  !  fpm16      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal, 2: Zolotarev)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !                                                 Rq: For Zolotarev fpm16=2,  ellipse is always a circle (fpm18 obsolete)
  !  Zne        (output)  COMPLEX DOUBLE PRECISION (fpm2): Complex coordinates of the  integration nodes for FEAST                                                      
  !  Wne        (output)  COMPLEX DOUBLE PRECISION (fpm2): Complex weights of the integration nodes for FEAST
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  !include "f90_noruntime_interface.fi"
  double precision :: Emin,Emax
  integer :: fpm2,fpm16,fpm18
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  double precision :: theta,xe,we,r,Emid     
  complex(kind=(kind(1.0d0))) :: zxe,zwe,jac 
  integer :: e
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  double precision, Parameter :: pi=3.1415926535897932d0
  double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0

  r=(Emax-Emin)/2.0d0
  Emid=Emin+r

!!!!!!!!!!!!!! Create Wne, Zne (half-contour)
  do e=1,fpm2 
     if(fpm16==0)then
        call dset_gauss_legendre(fpm2,e,xe,we) !! Gauss-points 
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+r*ONEC*cos(theta)+r*(DZERO,DONE)*(fpm18*1d-2)*sin(theta)
        jac=(r*(DZERO,DONE)*sin(theta)+ONEC*r*(fpm18*1d-2)*cos(theta))
        Wne(e)= (DONE/4.0d0)*we*jac
     elseif(fpm16==2)then
        call dset_zolotarev(fpm2,e,zxe,zwe)
        Zne(e)=zxe*r+Emid*ONEC
        Wne(e)=zwe*r
     elseif(fpm16==1)then
        theta=pi-(pi/fpm2)/2.0d0-(pi/fpm2)*(e-1)
        Zne(e)=Emid*ONEC+r*ONEC*cos(theta)+r*(DZERO,DONE)*(fpm18*1d-2)*sin(theta)
        jac=(r*(DZERO,DONE)*sin(theta)+ONEC*r*(fpm18*1d-2)*cos(theta))
        Wne(e)= (DONE/(2.0d0*fpm2))*jac
     endif
  enddo

end subroutine zfeast_contour



subroutine cfeast_contour(Emin,Emax,fpm2,fpm16,fpm18,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN FEAST INTEGRATION NODES/WEIGHTS
  !                         FROM an INPUT ELLIPSOID CONTOUR SYMMETRIC with REAL AXIS- return half-contour only
  !  
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emin,Emax  (input)   REAL SINGLE PRECISION    :search interval (horizontal ellipsoid axe b)
  !  fpm2       (input)   INTEGER                  :# nodes in half-contour 
  !  fpm16      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal, 2: Zolotarev)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !                                                 Rq: For Zolotarev fpm16=2,  ellipse is always a circle (fpm18 obsolete)
  !  Zne        (output)  COMPLEX SINGLE PRECISION (fpm2): Complex coordinates of the  integration nodes for FEAST                                                      
  !  Wne        (output)  COMPLEX SINGLE PRECISION (fpm2): Complex weights of the integration nodes for FEAST
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  !include "f90_noruntime_interface.fi"
  real :: Emin,Emax
  integer :: fpm2,fpm16,fpm18
  complex,dimension(*) :: Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  real :: theta,xe,we,r,Emid

  complex :: zxe,zwe,jac 
  integer :: e
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
  real, Parameter :: pi=3.1415926535897932e0
  real, parameter :: ba=-pi/2.0e0, ab=pi/2.0e0

  r=(Emax-Emin)/2.0e0
  Emid=Emin+r

!!!!!!!!!!!!!! Create Wne, Zne (half-contour)
  do e=1,fpm2 
     if(fpm16==0)then
        call sset_gauss_legendre(fpm2,e,xe,we) !! Gauss-points 
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+r*ONEC*cos(theta)+r*(SZERO,SONE)*(fpm18*1e-2)*sin(theta)
        jac=(r*(SZERO,SONE)*sin(theta)+ONEC*r*(fpm18*1e-2)*cos(theta))
        Wne(e)= (SONE/4.0e0)*we*jac
     elseif(fpm16==2)then
        call sset_zolotarev(fpm2,e,zxe,zwe)
        Zne(e)=zxe*r+Emid*ONEC
        Wne(e)=zwe*r
     elseif(fpm16==1)then
        theta=pi-(pi/fpm2)/2.0e0-(pi/fpm2)*(e-1)
        Zne(e)=Emid*ONEC+r*ONEC*cos(theta)+r*(SZERO,SONE)*(fpm18*1e-2)*sin(theta)
        jac=(r*(SZERO,SONE)*sin(theta)+ONEC*r*(fpm18*1e-2)*cos(theta))
        Wne(e)= (SONE/(2.0e0*fpm2))*jac
     endif
  enddo

end subroutine cfeast_contour


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_gcontour(Emid,r,fpm8,fpm17,fpm18,fpm19,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN FEAST INTEGRATION NODES/WEIGHTS
  !                         FROM an INPUT ELLIPSOID CONTOUR LOCATED in COMPLEX PLANE- return full contour
  !  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emid,r     (input)   COMPLEX/REAL DOUBLE PRECISION :search interval (middle,radius)- horizontal ellipsoid axe b=dble([Emid-r+Emid+r])
  !  fpm8       (input)   INTEGER                  :# nodes in full contour 
  !  fpm17      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !  fpm19      (input)   INTEGER                  : Rotation angle for the ellipse in degree [-180:180]
  !  Zne        (output)  COMPLEX DOUBLE PRECISION (fpm8): Complex coordinates of the  integration nodes for FEAST                                                      
  !  Wne        (output)  COMPLEX DOUBLE PRECISION (fpm8): Complex weights of the integration nodes for FEAST
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  !include "f90_noruntime_interface.fi"
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0,fpm8,fpm17,fpm18,fpm19
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  double precision :: theta,xe,we
  complex(kind=(kind(1.0d0))) :: Ze, jac, aux
  complex(kind=(kind(1.0d0))) :: zxe,zwe
  integer :: e,i,j,k,infoloc
  complex(kind=(kind(1.0d0))) :: res
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  double precision, Parameter :: pi=3.1415926535897932d0
  !  double precision, parameter :: ba=-pi, ab=pi
  double precision, parameter :: ba=-pi/2.0d0, ab=pi/2.0d0
  complex(kind=(kind(1.0d0))) :: nr

  !! ellipse axis rotation  
  theta=(fpm19/180d0)*pi
  nr=r*(ONEC*cos(theta)+(DZERO,DONE)*sin(theta))


  if(fpm17==0)then


     do e=1,fpm8/2
        call dset_gauss_legendre(fpm8/2,e,xe,we) !! Gauss-points - upper half
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+nr*ONEC*cos(theta)+nr*(DZERO,DONE)*(fpm18*1d-2)*sin(theta)
        jac=(nr*(DZERO,DONE)*sin(theta)+ONEC*nr*(fpm18*1d-2)*cos(theta))
        Wne(e) = (ONEC/(4.0d0))*we*jac
     enddo

     do e=fpm8/2+1,fpm8
        call dset_gauss_legendre(fpm8-fpm8/2,e-fpm8/2,xe,we) !! Gauss-points - lower half
        theta=-ba*xe-ab
        Zne(e)=Emid*ONEC+nr*ONEC*cos(theta)+nr*(DZERO,DONE)*(fpm18*1d-2)*sin(theta)
        jac=nr*(DZERO,DONE)*sin(theta)+ONEC*nr*(fpm18*1d-2)*cos(theta)
        Wne(e) = (ONEC/(4.0d0))*we*jac
     enddo



  elseif(fpm17==1)then

     do e=1,fpm8
        theta=pi-(2*pi/fpm8)/2.0d0-(2*pi/fpm8)*(e-1)
        Zne(e)=Emid*ONEC+nr*ONEC*cos(theta)+nr*(DZERO,DONE)*(fpm18*1d-2)*sin(theta)
        jac=(nr*(DZERO,DONE)*sin(theta)+ONEC*nr*(fpm18*1d-2)*cos(theta))
        Wne(e)=(ONEC/fpm8)*jac
     enddo

  endif





end subroutine zfeast_gcontour




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cfeast_gcontour(Emid,r,fpm8,fpm17,fpm18,fpm19,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN FEAST INTEGRATION NODES/WEIGHTS
  !                         FROM an INPUT ELLIPSOID CONTOUR LOCATED in COMPLEX PLANE- return full contour
  !  
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emid,r     (input)   COMPLEX/REAL SINGLE PRECISION :search interval (middle,radius)- horizontal ellipsoid axe b=dble([Emid-r+Emid+r])
  !  fpm8       (input)   INTEGER                  :# nodes in full contour 
  !  fpm17      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !  fpm19      (input)   INTEGER                  : Rotation angle for the ellipse in degree [-180:180]
  !  Zne        (output)  COMPLEX SINGLE PRECISION (fpm8): Complex coordinates of the  integration nodes for FEAST                                                      
  !  Wne        (output)  COMPLEX SINGLE PRECISION (fpm8): Complex weights of the integration nodes for FEAST
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  !include "f90_noruntime_interface.fi"
  complex :: Emid
  real :: r
  integer :: M0,fpm8,fpm17,fpm18,fpm19
  complex,dimension(*) :: Zne,Wne
  real :: theta,xe,we
  complex :: Ze, jac, aux
  complex :: zxe,zwe
  integer :: fpm22,e,i,j,k,infoloc
  complex :: res
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
  real, Parameter :: pi=3.1415926535897932e0
  real, parameter :: ba=-pi/2.0e0, ab=pi/2.0e0
  complex :: nr

  !! ellipse axis rotation  
  theta=(fpm19/180e0)*pi
  nr=r*(ONEC*cos(theta)+(SZERO,SONE)*sin(theta))



  if(fpm17==0)then

     do e=1,fpm8/2
        call sset_gauss_legendre(fpm8/2,e,xe,we) !! Gauss-points - upper half
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+nr*ONEC*cos(theta)+nr*(SZERO,SONE)*(fpm18*1e-2)*sin(theta)
        jac=(nr*(SZERO,SONE)*sin(theta)+ONEC*nr*(fpm18*1e-2)*cos(theta))
        Wne(e) = (ONEC/(4.0e0))*we*jac
     enddo

     do e=fpm8/2+1,fpm8
        call sset_gauss_legendre(fpm8-fpm8/2,e-fpm8/2,xe,we)!! Gauss-points - lower half
        theta=-ba*xe-ab
        Zne(e)=Emid*ONEC+nr*ONEC*cos(theta)+nr*(SZERO,SONE)*(fpm18*1e-2)*sin(theta)
        jac=(nr*(SZERO,SONE)*sin(theta)+ONEC*nr*(fpm18*1e-2)*cos(theta))
        Wne(e) = (ONEC/(4.0e0))*we*jac
     enddo


  elseif(fpm17==1)then

     do e=1,fpm8
        theta=pi-(2*pi/fpm8)/2.0e0-(2*pi/fpm8)*(e-1)
        Zne(e)=Emid*ONEC+nr*ONEC*cos(theta)+nr*(SZERO,SONE)*(fpm18*1e-2)*sin(theta)
        jac=(nr*(SZERO,SONE)*sin(theta)+ONEC*nr*(fpm18*1e-2)*cos(theta))
        Wne(e)=(ONEC/fpm8)*jac
     enddo

  endif

end subroutine cfeast_gcontour


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!33333333333!!!!! CUSTOM CONTOUR FEAST Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_customcontour(fpm8,N,Nedge,Tedge,Zedge,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: ASSIST the USER in GENERATING FEAST INTEGRATION NODES/WEIGHTS
  !                         FROM an INPUT CUSTOM CONTOUR 
  !  
  !  COMPLEX DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  fpm8  (input)        INTEGER    : # of integration nodes- should be equal to sum(Nedge(1:N)) 
  !  N     (input)        INTEGER    : # of geometrical pieces comprising the entire contour
  !  Nedge (input)        INTEGER(N) :# local number of interval by geometrical piece
  !  Tedge (input)        INTEGER(N) : type (nature) of each geometrical piece
  !                                     0:line segment, 100: half-circle, other: ellipse ratio=Tedge(k)/100
  !                                      
  !  Zedge (input)        COMPLEX DOUBLE PRECISION(N): First complex coordinate of each
  !                                                     geometrical piece  
  !  Zne   (output)       COMPLEX DOUBLE PRECISION(fpm8): Complex coordinates of the  
  !                                                         integration nodes for FEAST
  !                                                           
  !  Wne   (output)       COMPLEX DOUBLE PRECISION(fpm8): Complex weights of the
  !                                                         integration nodes for FEAST
  !=====================================================================
  ! James Kestyn 2013-2015
  ! ====================================================================   
  implicit none
  !include "f90_noruntime_interface.fi"
  integer :: fpm8
  integer :: N
  integer,dimension(*) :: Nedge, Tedge
  complex(kind=(kind(1.0d0))),dimension(*) :: Zedge
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  complex(kind=(kind(1.0d0))),dimension(1:fpm8+N) :: tZne,tWne
  complex(kind=(kind(1.0d0))) :: a,b,Emid
  double precision :: r,r2,theta,angle
  integer :: Nz
  integer :: i,j,cnt
  double precision, Parameter :: pi=3.1415926535897932d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Find total # point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Nz = fpm8+N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Find Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tWne(1:Nz)=(0.0d0,0.0d0)
  cnt = 0
  do i=1,N
     a = Zedge(i)
     if(i < N)then
        b = Zedge(i+1)
     else
        b= Zedge(1)
     endif

     select case(Tedge(i))

     case(0) !line segment
        do j=1,Nedge(i)+1
           tZne(cnt+j) = a + (j-1)*(b-a)/(Nedge(i)) 
           tWne(cnt+j) = (0.0d0,1.0d0)*(b-a)/( 2.0d0 * pi * (Nedge(i)) )
           if( j==1 .or. j==Nedge(i)+1 ) tWne(cnt+j) = 0.5d0*tWne(cnt+j)
        enddo

     case default ! Half Ellipse
        angle = atan2(aimag(b-a),dble(b-a))
        Emid = (a+b)/2.0d0
        r = abs(a-b)/2.0d0
        r2= Tedge(i)*r/100.0d0
        do j=1,Nedge(i)+1
           theta = pi-1.0d0*(j-1)*pi/(Nedge(i))
           tZne(cnt+j) = Emid +  (1.0d0,0.0d0)*(r*cos(theta)*cos(angle)-r2*sin(theta)*sin(angle)) + (0.0d0,1.0d0)*(r*cos(theta)*sin(angle)+r2*sin(theta)*cos(angle))
           tWne(cnt+j) = -r*sin(theta)*cos(angle) - r2*cos(theta)*sin(angle) + (0.0d0,1.0d0)*(-r*sin(theta)*sin(angle) + r2*cos(theta)*cos(angle))
           tWne(cnt+j) = tWne(cnt+j)*(0.0d0,-1.0d0)/(2*(Nedge(i))) 
           if( j==1 .or. j==Nedge(i)+1 ) tWne(cnt+j) = 0.5d0*tWne(cnt+j)
        enddo

     end select

     cnt = cnt+Nedge(i)+1
  enddo
!!! reorder array
  Wne(1:fpm8) = (0.0d0,0.0d0)
  Wne(1:Nedge(1)+1) = Wne(1:Nedge(1)+1)+tWne(1:Nedge(1)+1)
  Zne(1:fpm8) = (0.0d0,0.0d0)
  Zne(1:Nedge(1)) = tZne(1:Nedge(1))
  cnt = Nedge(1)
  do i=2,N-1
     Zne(cnt+1:cnt+1+Nedge(i)) = tZne(cnt+i:cnt+i+Nedge(i))
     Wne(cnt+1:cnt+1+Nedge(i)) = Wne(cnt+1:cnt+1+Nedge(i))+tWne(cnt+i:cnt+i+Nedge(i))
     cnt = cnt + Nedge(i)
  enddo
  Zne(cnt+1:cnt+1+Nedge(N)-1) = tZne(cnt+i:cnt+i+Nedge(N)-1)
  Wne(cnt+1:cnt+1+Nedge(N)-1) = Wne(cnt+1:cnt+1+Nedge(N)-1)+tWne(cnt+N:cnt+N+Nedge(N)-1)
  Wne(1) = Wne(1) + tWne(Nz)

end subroutine zfeast_customcontour


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cfeast_customcontour(fpm8,N,Nedge,Tedge,Zedge,Zne,Wne)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: ASSIST the USER in GENERATING FEAST INTEGRATION NODES/WEIGHTS
  !                         FROM an INPUT CUSTOM CONTOUR 
  !  
  !  COMPLEX DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  fpm8  (input)        INTEGER    : # of integration nodes- should be equal to sum(Nedge(1:N)) 
  !  N     (input)        INTEGER    : # of geometrical pieces comprising the entire contour
  !  Nedge (input)        INTEGER(N) :# local number of interval by geometrical piece
  !  Tedge (input)        INTEGER(N) : type (nature) of each geometrical piece
  !                                     0:line segment, 100: half-circle, other: ellipse ratio=Tedge(k)/100
  !  Zedge (input)        COMPLEX SINGLE PRECISION(N): First complex coordinate of each
  !                                                     geometrical piece  
  !  Zne   (output)       COMPLEX SINGLE PRECISION(fpm8): Complex coordinates of the  
  !                                                         integration nodes for FEAST
  !                                                           
  !  Wne   (output)       COMPLEX SINGLE PRECISION(fpm8): Complex weights of the
  !                                                         integration nodes for FEAST
  !=====================================================================
  ! James Kestyn 2013-2015
  ! ====================================================================   
  implicit none
  !include "f90_noruntime_interface.fi"
  integer :: fpm8
  integer :: N
  integer,dimension(*) :: Nedge, Tedge
  complex,dimension(*) :: Zedge
  complex,dimension(*) :: Zne,Wne
  complex,dimension(1:fpm8+N) :: tZne,tWne
  complex :: a,b,Emid
  real :: r,r2,theta,angle
  integer :: Nz
  integer :: i,j,cnt
  real, Parameter :: pi=3.14159265e0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Find total # point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Nz = fpm8+N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Find Zne, Wne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tWne(1:Nz)=(0.0e0,0.0e0)
  cnt = 0
  do i=1,N
     a = Zedge(i)
     if(i < N)then
        b = Zedge(i+1)
     else
        b= Zedge(1)
     endif

     select case(Tedge(i))

     case(0) !line segment
        do j=1,Nedge(i)+1
           tZne(cnt+j) = a + (j-1)*(b-a)/(Nedge(i)) 
           tWne(cnt+j) = (0.0e0,1.0e0)*(b-a)/( 2.0e0 * pi * (Nedge(i)) )
           if( j==1 .or. j==Nedge(i)+1 ) tWne(cnt+j) = 0.5e0*tWne(cnt+j)
        enddo

     case default ! Half Ellipse
        angle = atan2(aimag(b-a),real(b-a))
        Emid = (a+b)/2.0e0
        r = abs(a-b)/2.0e0
        r2= Tedge(i)*r/100.0e0
        do j=1,Nedge(i)+1
           theta = pi-1.0*(j-1)*pi/(Nedge(i))
           tZne(cnt+j) = Emid +  (1.0e0,0.0e0)*(r*cos(theta)*cos(angle)-r2*sin(theta)*sin(angle)) + (0.0e0,1.0e0)*(r*cos(theta)*sin(angle)+r2*sin(theta)*cos(angle))
           tWne(cnt+j) = -r*sin(theta)*cos(angle) - r2*cos(theta)*sin(angle) + (0.0e0,1.0e0)*(-r*sin(theta)*sin(angle) + r2*cos(theta)*cos(angle))
           tWne(cnt+j) = tWne(cnt+j)*(0.0e0,-1.0e0)/(2*(Nedge(i))) 
           if( j==1 .or. j==Nedge(i)+1 ) tWne(cnt+j) = 0.5*tWne(cnt+j)
        enddo

     end select

     cnt = cnt+Nedge(i)+1
  enddo
!!! reorder array
  Wne(1:fpm8) = (0.0e0,0.0e0)
  Wne(1:Nedge(1)+1) = Wne(1:Nedge(1)+1)+tWne(1:Nedge(1)+1)
  Zne(1:fpm8) = (0.0e0,0.0e0)
  Zne(1:Nedge(1)) = tZne(1:Nedge(1))
  cnt = Nedge(1)
  do i=2,N-1
     Zne(cnt+1:cnt+1+Nedge(i)) = tZne(cnt+i:cnt+i+Nedge(i))
     Wne(cnt+1:cnt+1+Nedge(i)) = Wne(cnt+1:cnt+1+Nedge(i))+tWne(cnt+i:cnt+i+Nedge(i))
     cnt = cnt + Nedge(i)
  enddo
  Zne(cnt+1:cnt+1+Nedge(N)-1) = tZne(cnt+i:cnt+i+Nedge(N)-1)
  Wne(cnt+1:cnt+1+Nedge(N)-1) = Wne(cnt+1:cnt+1+Nedge(N)-1)+tWne(cnt+N:cnt+N+Nedge(N)-1)
  Wne(1) = Wne(1) + tWne(Nz)

end subroutine cfeast_customcontour



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!44444444444!!!!! RATIONAL FEAST Routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfeast_rational(Emin,Emax,fpm2,fpm16,fpm18,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of REAL eigenvalues
  !                         USING an INPUT ELLIPSOID CONTOUR SYMMETRIC with REAL AXIS
  !  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emin,Emax  (input)   REAL DOUBLE PRECISION    :search interval (horizontal ellipsoid axe b)
  !  fpm2       (input)   INTEGER                  :# nodes in half-contour 
  !  fpm16      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal, 2: Zolotarev)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !                                                 Rq: For Zolotarev fpm16=2,  ellipse is always a circle (fpm18 obsolete)
  !  eig        (input)   REAL DOUBLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  REAL DOUBLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  double precision :: Emin,Emax
  integer :: M0,fpm2,fpm16,fpm18
  double precision,dimension(*) :: eig,f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  complex(kind=(kind(1.0d0))),dimension(1:fpm2) :: Zne, Wne
  complex(kind=(kind(1.0d0))) :: zxe,zwe
  !double precision :: r

  call zfeast_contour(Emin,Emax,fpm2,fpm16,fpm18,Zne,Wne)
  call dfeast_rationalx(Zne,Wne,fpm2,Eig,M0,f)

  if (fpm16==2) then ! zolotarev initialization
     !r=(Emax-Emin)/2.0d0
     call dset_zolotarev(fpm2,0,zxe,zwe)
     f(1:M0)=f(1:M0)+zwe!*r
  endif

end subroutine dfeast_rational




subroutine sfeast_rational(Emin,Emax,fpm2,fpm16,fpm18,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of REAL eigenvalues
  !                         USING an INPUT ELLIPSOID CONTOUR SYMMETRIC with REAL AXIS
  !  
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emin,Emax  (input)   REAL SINGLE PRECISION    :search interval (horizontal ellipsoid axe b)
  !  fpm2       (input)   INTEGER                  :# nodes in half-contour 
  !  fpm16      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal, 2: Zolotarev)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !                                                 Rq: For Zolotarev fpm16=2,  ellipse is always a circle (fpm18 obsolete)
  !  eig        (input)   REAL SINGLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  REAL SINGLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  real :: Emin,Emax
  integer :: M0,fpm2,fpm16,fpm18
  real,dimension(*) :: eig,f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  complex,dimension(1:fpm2) :: Zne, Wne
  complex :: zxe,zwe
  !real :: r

  call cfeast_contour(Emin,Emax,fpm2,fpm16,fpm18,Zne,Wne)
  call sfeast_rationalx(Zne,Wne,fpm2,Eig,M0,f)

  if (fpm16==2) then ! zolotarev initialization
     !r=(Emax-Emin)/2.0e0
     call sset_zolotarev(fpm2,0,zxe,zwe)
     f(1:M0)=f(1:M0)+zwe!*r
  endif

end subroutine sfeast_rational





subroutine dfeast_rationalx(Zne,Wne,fpm2,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of REAL eigenvalues
  !                         USING an a CUSTOM CONTOUR Symmetric with Real axis
  !  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Zne        (input)   COMPLEX DOUBLE PRECISION (fpm2): Complex coordinates of the  integration nodes for FEAST                             
  !  Wne        (input)   COMPLEX DOUBLE PRECISION (fpm2): Complex weights of the integration nodes for FEAST
  !  fpm2       (input)   INTEGER                  :# nodes in half-contour 
  !  eig        (input)   REAL DOUBLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  REAL DOUBLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  integer :: M0,fpm2
  double precision,dimension(*) :: eig,f
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!
  integer :: i,e
  double precision :: res

  do i=1,M0
     f(i)=0.0d0
     do e=1,fpm2 ! loop half-contour
        f(i)=f(i)+2.0d0*dble(Wne(e)/(Zne(e)-eig(i)))   ! factor 2 for the two contours       
     enddo
  enddo

end subroutine dfeast_rationalx




subroutine sfeast_rationalx(Zne,Wne,fpm2,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of REAL eigenvalues
  !                         USING an a CUSTOM CONTOUR Symmetric with Real axis
  !  
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Zne        (input)   COMPLEX SINGLE PRECISION (fpm2): Complex coordinates of the  integration nodes for FEAST                             
  !  Wne        (input)   COMPLEX SINGLE PRECISION (fpm2): Complex weights of the integration nodes for FEAST
  !  fpm2       (input)   INTEGER                  :# nodes in half-contour 
  !  eig        (input)   REAL SINGLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  REAL SINGLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  integer :: M0,fpm2
  real,dimension(*) :: eig,f
  complex,dimension(*) :: Zne,Wne
!!!!!!!!!!!!!!!!
  integer :: i,e
  real :: res

  do i=1,M0
     f(i)=0.0e0
     do e=1,fpm2 ! loop half-contour
        f(i)=f(i)+2.0e0*dble(Wne(e)/(Zne(e)-eig(i)))   ! factor 2 for the two contours       
     enddo
  enddo

end subroutine sfeast_rationalx





subroutine zfeast_grational(Emid,r,fpm8,fpm17,fpm18,fpm19,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of COMPLEX eigenvalues
  !                         USING an INPUT ELLIPSOID CONTOUR located in COMPLEX PLANE
  !  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emid,r     (input)   COMPLEX/REAL DOUBLE PRECISION :search interval (middle,radius)- horizontal ellipsoid axe b=dble([Emid-r+Emid+r])
  !  fpm8       (input)   INTEGER                  :# nodes in full contour 
  !  fpm17      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !  fpm19      (input)   INTEGER                  : Rotation angle for the ellipse in degree [-180:180]
  !  eig        (input)   COMPLEX DOUBLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  COMPLEX DOUBLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer :: M0,fpm8,fpm17,fpm18,fpm19
  complex(kind=(kind(1.0d0))),dimension(*) :: f,eig
!!!!
  complex(kind=(kind(1.0d0))),dimension(1:fpm8) :: Zne,Wne

  call zfeast_gcontour(Emid,r,fpm8,fpm17,fpm18,fpm19,Zne,Wne)
  call zfeast_grationalx(Zne,Wne,fpm8,Eig,M0,f)

end subroutine zfeast_grational



subroutine cfeast_grational(Emid,r,fpm8,fpm17,fpm18,fpm19,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of COMPLEX eigenvalues
  !                         USING an INPUT ELLIPSOID CONTOUR located in COMPLEX PLANE
  !  
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Emid,r     (input)   COMPLEX/REAL SINGLE PRECISION :search interval (middle,radius)- horizontal ellipsoid axe b=dble([Emid-r+Emid+r])
  !  fpm8       (input)   INTEGER                  :# nodes in full contour 
  !  fpm17      (input)   INTEGER                  :Type of integration  (0: Gauss, 1: Trapezoidal)
  !  fpm18      (input)   INTEGER                  :fpm18*0.01 is ratio a/b (a is vertical axis) (e.g. circle is 100)
  !  fpm19      (input)   INTEGER                  : Rotation angle for the ellipse in degree [-180:180]
  !  eig        (input)   COMPLEX SINGLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  COMPLEX SINGLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  complex :: Emid
  real :: r
  integer :: M0,fpm8,fpm17,fpm18,fpm19
  complex,dimension(*) :: f,eig
!!!!
  complex,dimension(1:fpm8) :: Zne,Wne

  call cfeast_gcontour(Emid,r,fpm8,fpm17,fpm18,fpm19,Zne,Wne)
  call cfeast_grationalx(Zne,Wne,fpm8,Eig,M0,f)

end subroutine cfeast_grational




subroutine zfeast_grationalx(Zne,Wne,fpm8,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of COMPLEX eigenvalues
  !                         USING an a CUSTOM CONTOUR located in Complex plane
  !  
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Zne        (input)   COMPLEX DOUBLE PRECISION (fpm8): Complex coordinates of the  integration nodes for FEAST                             
  !  Wne        (input)   COMPLEX DOUBLE PRECISION (fpm8): Complex weights of the integration nodes for FEAST
  !  fpm8       (input)   INTEGER                  :# nodes in full-contour 
  !  eig        (input)   COMPLEX DOUBLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  COMPLEX DOUBLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  integer :: M0,fpm8
  complex(kind=(kind(1.0d0))),dimension(*) :: Zne,Wne
  complex(kind=(kind(1.0d0))),dimension(*) :: eig,f
!!!!!!!!!!!!!!!!
  integer :: i,e
  do i=1,M0
     f(i)=(0.0d0,0.0d0)
     do e=1,fpm8 ! loop-full contour
        f(i)=f(i)+Wne(e)/(Zne(e)-eig(i))   !
     enddo
  enddo
end subroutine zfeast_grationalx


subroutine cfeast_grationalx(Zne,Wne,fpm8,Eig,M0,f)
  !  Purpose 
  !  =======
  !  FEAST UTILITY ROUTINE: RETURN VALUES of RATIONAL/SELECTION FUNCTION of a set of COMPLEX eigenvalues
  !                         USING an a CUSTOM CONTOUR located in Complex plane
  !  
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========   
  !  Zne        (input)   COMPLEX SINGLE PRECISION (fpm8): Complex coordinates of the  integration nodes for FEAST                             
  !  Wne        (input)   COMPLEX SINGLE PRECISION (fpm8): Complex weights of the integration nodes for FEAST
  !  fpm8       (input)   INTEGER                  :# nodes in full-contour 
  !  eig        (input)   COMPLEX SINGLE PRECISION (M0) : Contains list of eigenvalues 
  !  M0         (input)   INTEGER                       :# of eigenvalues
  !  f          (output)  COMPLEX SINGLE PRECISION (M0): Values of the rational functions                                                      
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ==================================================================== 
  implicit none
  integer :: M0,fpm8
  complex,dimension(*) :: Zne,Wne
  complex,dimension(*) :: eig,f
!!!!!!!!!!!!!!!!
  integer :: i,e
  do i=1,M0
     f(i)=(0.0e0,0.0e0)
     do e=1,fpm8 ! loop-full contour
        f(i)=f(i)+Wne(e)/(Zne(e)-eig(i))   !
     enddo
  enddo
end subroutine cfeast_grationalx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!555555555!!!!! Routines used by PFEAST (MPI needed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPI




subroutine pfeastinit(fpm,NEW_COMM_WORLD,nb_procs3)
  implicit none
  !  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Define Communication Worlds for PFEAST
!!!!!
!!!!! Inputs:
!!!!  ------
!!!!! NEW_COMM_WORLD (communicator)
!!!!!  * Custom Communicator for a given search interval obtained after L1 user parallelism
!!!!!  * MPI_COMM_WORLD otherwise
!!!!!  Rq: includes nb_procs for a given interval
!!!!!
!!!!! nb_procs3 (input/ouput) (integer): # of processor to consider by the user at L3 level
!!!!!                                    if nb_procs<nb_procs3 then ==> nb_procs3=nb_procs
!!!!!                                    if nb_procs>nb_procs3 then (remark: nb_procs=nb_procs2*nb_procs3)
!!!!!                                          if nb_procs3 is not a divider of nb_procs ==> nb_procs3=nb_procs  
!!!!!                                                
!!!!!  
!!!!! Outputs:
!!!!! -------
!!!!! fpm- integer(64) : * default feast parameters as returned by feastinit routine   
!!!!!                    * in addition:
!!!!                       fpm(49) == SOLVER_COMM_WORLD  -- Process associated with matrix partition (L3)
!!!!!                      fpm(9) == RCI_COMM_WORLD     -- Process associated with a given contour point (L2)
!!!!!
!!!!! Rq: optimally we must have  nb_procs=nb_procs2*nb_procs3
!!!!!      Performance increase along with memory requirement if 
!!!!!      *nb_procs2 from 1 to fpm(2) for (symmetric) or fpm(8) (for non-symmetric)
!!!!!  
  !=====================================================================
  ! James Kestyn 2016-2018
  ! Eric Polizzi 2019
  ! ====================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer,dimension(*) :: fpm
  integer :: NEW_COMM_WORLD
  integer :: nb_procs, nb_procs2, nb_procs3, n3
  integer :: rank,code,i
  integer :: SOLVER_COMM_WORLD, RCI_COMM_WORLD
  integer :: solver_rank, rci_rank, color, key
  integer :: NEW_NEW_COMM_WORLD,color2,new_rank

  n3=nb_procs3 !Rq: cannot change nb_procs3 in argument if needed (think of a number constant)
  
  call feastinit(fpm)
  rank=0
  call MPI_COMM_RANK(NEW_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(NEW_COMM_WORLD,nb_procs,code)

  !! global L1 communicator
  fpm(59)=NEW_COMM_WORLD

  !! recast value of nb_procs3 if needed
  if (nb_procs<n3) n3=nb_procs

  !! find nb_procs2 (optimally nbprocs=nb_procs2*nbprocs3)
  nb_procs2=nb_procs/n3 

  !! recast value of nb_procs3 if needed
  if (nb_procs2*n3/=nb_procs) n3=nb_procs
  !! Rq: we use this constraint since it would difficult 
  !!     to ignore some extra mpi process (if nbprocs/=nb_procs2*nbprocs3)
  !!     from an user prespective the unused mpi process would not contain any solutions 
  !!     which would be problematic.
  !!     we decide then on having nb_procs (or L1) to be fixed


!!! Each color is associated with a single contour point
!!! Each mpi process within a color is for the solver 
  color=0
  do i=1,nb_procs2
     if(rank >= (i-1)*n3) then
        color = i
     endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Split by contour point - local solver rank
  ! key=0
  key=rank
  call MPI_COMM_SPLIT(NEW_COMM_WORLD,color,key,SOLVER_COMM_WORLD,code)
  call MPI_COMM_RANK(SOLVER_COMM_WORLD,solver_rank,code) ! local rank
  fpm(49) = SOLVER_COMM_WORLD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Split by solver_rank - local rci rank
  !key=1
  key=rank
  call MPI_COMM_SPLIT(NEW_COMM_WORLD,solver_rank,key,RCI_COMM_WORLD,code)
  call MPI_COMM_RANK(RCI_COMM_WORLD,rci_rank,code) ! local rank
  fpm(9) = RCI_COMM_WORLD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine pfeastinit



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pfeast_distribution_type(N,isa,jsa,L3_COMM_WORLD,distype)
  implicit none
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine is used to determined if the matrix isa,jsa
  !               is stored locally in parallel among L3 mpi processes (row distribution)
  !               or stored globally in all L3 processors
  !
  !  isa,jsa (input) INTEGER(:) - CSR format of a given matrix A 
  !            
  !  L3_COMM_WORLD (input) INTEGER - MPI communicator
  !
  ! distype (output) INTEGER- distribution type
  !                           0: A is stored globally on all MPI
  !                           1: A is locall row distributed among MPI
  !=====================================================================
  ! Eric Polizzi 2019
  ! ====================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: n,L3_COMM_WORLD,distype
  integer,dimension(*) :: isa,jsa
!!!!!!!!!!!!!!!!
  integer :: Ncol,i,code

  !! find the maximum column number for jsa
  Ncol=1
  do i=1,isa(N+1)-1 ! N local or global
     if (jsa(i)>Ncol) Ncol=jsa(i)
  enddo
  if (L3_COMM_WORLD/=0)  call MPI_ALLREDUCE(MPI_IN_PLACE,Ncol,1,MPI_INTEGER,MPI_MAX,L3_COMM_WORLD,code) !Ncol is global
  distype=0 ! global distribution (default) 
  if (Ncol/=N) then !means that N is not global
     distype=1 !local distribution
  end if

end subroutine pfeast_distribution_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



#endif
