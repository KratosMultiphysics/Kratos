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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! NUMERICAL FORTRAN LIBRARY
!!!!!! E. Polizzi's research lab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Remark: compilation
!!!         possible directives: MKL, MPI
!!!         dependencies: dzlsprim.f90, sclsprim.f90



!!11111111111111111111!!!!!!! Quadrature

!!! Integration nodes and weights
! sset_gauss_legendre
! dset_gauss_legendre
! sset_zolotarev
! dset_zolotarev


!!22222222222222222222!!!!!!!! Mat-vec multiplication

!!! List of wrapper mat-vec routines (double/single precision)
!!!   call MKL with compiler directives, libsprim otherwise

!wdcsrmm, wzcsrmm    :: mat-vec multiplication =>  b=alpha*A*x+beta*b
!wzhcsrmm            :: mat-vec multiplication =>  b=alpha*A*x+beta*b (A is complex Hermitian if UPLO=L,U)
!wscsrmm, wccsrmm    :: mat-vec multiplication =>  b=alpha*A*x+beta*b
!wchcsrmm            :: mat-vec multiplication =>  b=alpha*A*x+beta*b (A is complex Hermitian if UPLO=L,U)


!!33333333333333333333!!!!!! Iterative Solvers and other stuff

!dbicgstab
!dbicgstab_rci

!zbicgstab
!zhbicgstab
!zbicgstab_rci

!cbicgstab
!cbicgstab_rci

!dqrgs  -  qr Modified GS-fact
!zqrgs  -  qr Modified GS-fact

!darnoldi
!zarnoldi
!zharnoldi
!dgarnoldi
!zhgarnoldi



!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI from here

!44444444444444444444 !!!!! MPI matrix distribution and mat-vec
!!!! MPI mat-vec routines using block csr storage  (mpi is needed)

!dcsr_distribute_row :: distribute the global matrix by local rows
!zcsr_distribute_row :: distribute the global matrix by local rows (complex)

!dget_local_block_csr:: redistribute a row distributed csr matrix into row/column csr local blocks among MPI processors
!zget_local_block_csr:: redistribute a row distributed csr matrix into row/column csr local blocks among MPI processors (complex double)
!cget_local_block_csr:: redistribute a row distributed csr matrix into row/column csr local blocks among MPI processors (complex single)

!dlbcsr_transpose :: transpose a matrix in local block csr format
!zlbcsr_transpose :: complex version

!csaddlbcsr ::     C=alpha*A+beta*B in local block csr format
!zdaddlbcsr ::     C=alpha*A+beta*B in local block csr format
!caddlbcsr ::     C=alpha*A+beta*B in local block csr format
!zaddlbcsr ::     C=alpha*A+beta*B in local block csr format
!zinc_addlbcsr ::  Increment matrix B by adding matrix A times alpha
!cinc_addlbcsr ::  Increment matrix B by adding matrix A times alpha

!dlbcsr2s      :: copy/convert real double precision lbcsr matrix into single precision
!zlbcsr2c      :: copy/convert complex double precision lbcsr matrix into single precision


!dlbcsr_uplo_to_csr:: convert a matrix in local block csr format from lower/upper to full
!zhlbcsr_uplo_to_csr :: (complex Hermitian version)
!zlbcsr_uplo_to_csr :: (complex Symmetric version)


!dlbcsr_distribute_row ::  redistribute a  row/column csr local blocks into a row distributed csr matrix among MPI processors
!zlbcsr_distribute_row :: complex version


!dlbcsrmm            :: MPI mat-vec multiplication =>  b=alpha*A*x+beta*b - real double precision version 
!zlbcsrmm            :: complex version
!zhlbcsrmm           :: complex Hermitian version
!clbcsrmm            :: complex single precision version 


!pdqrgs  -  qr Modified GS-fact
!pzqrgs  -  qr Modified GS-fact

!pdbicgstab
!pdbicgstab_rci
!pzbicgstab
!pzhbicgstab
!pzbicgstab_rci
!pcbicgstab
!pcbicgstab_rci


!pdarnoldi
!pzarnoldi
!pzharnoldi
!pdgarnoldi
!pzhgarnoldi



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!11111111111111111111!!!!!!! Quadrature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine sset_gauss_legendre(nbe,e,xe,we)
  !  Purpose 
  !  =======
  !  GAUSS-LEGENDRE Quadrature ROUTINE- Return weight and coordinate of the e^th node from nbe total # nodes
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  nbe        (input)        INTEGER: total # nodes for Gauss quadrature
  !  e          (input)        INTEGER: e^th node
  !  xe         (output)       REAL SINGLE PRECISION:  Gauss coordinate for node e 
  !  we         (output)       REAL SINGLE PRECISION:  Gauss weight for node e
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================
  implicit none
  integer :: nbe,e
  real :: xe,we
!!!!!!!!!!
  double precision ::dxe,dwe


  call  dset_gauss_legendre(nbe,e,dxe,dwe)

  xe=real(dxe)
  we=real(dwe)

end subroutine sset_gauss_legendre



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dset_gauss_legendre(nbe,e,xe,we)
  !  Purpose 
  !  =======
  !  GAUSS-LEGENDRE Quadrature ROUTINE- Return weight and coordinate of the e^th node from nbe total # nodes
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  nbe        (input)        INTEGER: total # nodes for Gauss quadrature
  !  e          (input)        INTEGER: e^th node
  !  xe         (output)       REAL DOUBLE PRECISION:  Gauss coordinate for node e 
  !  we         (output)       REAL DOUBLE PRECISION:  Gauss weight for node e
  !=====================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================
  implicit none
  integer :: nbe,e
  double precision :: xe,we
!!!!!!!!!!


  select case(nbe)


  case(1)
     we=2.0d0
     xe=0.0d0


  case(2)
     if (e==1) then
        we=1.0000000000000000d0
        xe=-0.5773502691896257d0
     elseif(e==2) then
        we=1.0000000000000000d0
        xe=0.5773502691896257d0
     endif


  case(3)
     select case(e)
     case(1)
        we=8.0d0/9.0d0
        xe=0.0d0
     case(2)
        we=5.0d0/9.0d0
        xe=sqrt(3.0d0/5.0d0)
     case(3)
        we=5.0d0/9.0d0 
        xe=-sqrt(3.0d0/5.0d0)
     end select


  case(4)
     select case(e)

     case(1,2)
        we=(18.0d0+sqrt(30.0d0))/36.0d0
        !6.52145154862546142644d-01
        xe=sqrt((3.0d0-2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
        !3.39981043584856264792d-01
     case(3,4)
        we=(18.0d0-sqrt(30.0d0))/36.0d0
        ! 3.47854845137453857383d-01
        xe=sqrt((3.0d0+2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
        !8.61136311594052575248d-01

     end select

     if (ibclr(e,0)==e) xe=-xe



  case(5)

     select case(e)
     case(1)
        we=128.0d0/225.0d0
        xe=0.0d0
     case(2,3)
        !      we=9.06179845938663992811d-01
        we=(322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
        !     xe= 5.38469310105683091018d-01 
        xe=(1.0d0/3.0d0)*sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))
        if (e==3) xe=-xe
     case(4,5)
        !     we=2.36926885056189087515d-01
        we=(322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
        !        xe=9.06179845938663992811d-01
        xe=(1.0d0/3.0d0)*sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))
        if (e==5) xe=-xe
     end select





  case(6)
     select case(e)
     case(1,2)
        xe=6.61209386466264513688d-01
        we= 3.60761573048138607569d-01
     case(3,4)
        xe=2.38619186083196908630d-01
        we=4.67913934572691047389d-01
     case(5,6)
        xe=9.32469514203152027832d-01
        we=1.71324492379170345043d-01
     end select

     if (ibclr(e,0)==e) xe=-xe



  case(7)
     if(e==1) then
        we=0.4179591836734694d0
        xe=0.0000000000000000d0
     elseif(e==2) then
        we=0.3818300505051189d0
        xe=0.4058451513773972d0
     elseif(e==3) then
        we=0.3818300505051189d0
        xe=-0.4058451513773972d0
     elseif(e==4) then
        we=0.2797053914892766d0
        xe=-0.7415311855993945d0
     elseif(e==5) then
        we=0.2797053914892766d0
        xe=0.7415311855993945d0
     elseif(e==6) then
        we=0.1294849661688697d0
        xe=-0.9491079123427585d0
     elseif(e==7) then
        we=0.1294849661688697d0
        xe=0.9491079123427585d0
     endif




  case(8)
     select case(e)
     case(1,2)
        xe=1.83434642495649804936d-01
        we=3.62683783378361982976d-01
     case(3,4)
        xe=5.25532409916328985830d-01
        we=3.13706645877887287338d-01
     case(5,6)
        xe=7.96666477413626739567d-01
        we=2.22381034453374470546d-01
     case(7,8)
        xe=9.60289856497536231661d-01
        we=1.01228536290376259154d-01
     end select
     if (ibclr(e,0)==e) xe=-xe




  case(9)
     if(e==1) then
        we=0.3302393550012598d0	
        xe=0.0000000000000000d0
     elseif(e==2) then
        we=0.1806481606948574d0	
        xe=-0.8360311073266358d0
     elseif(e==3) then
        we=0.1806481606948574d0	
        xe=0.8360311073266358d0
     elseif(e==4) then
        we=0.0812743883615744d0	
        xe=-0.9681602395076261d0
     elseif(e==5) then
        we=0.0812743883615744d0	
        xe=0.9681602395076261d0
     elseif(e==6) then
        we=0.3123470770400029d0 
        xe=-0.3242534234038089d0
     elseif(e==7) then
        we=0.3123470770400029d0	
        xe=0.3242534234038089d0
     elseif(e==8) then
        we=0.2606106964029354d0	
        xe=-0.6133714327005904d0
     elseif(e==9) then
        we=0.2606106964029354d0	
        xe=0.6133714327005904d0
     end if






  case(10)

     select case(e)
     case(1,2)
        we=2.95524224714752870187d-01
        xe=1.48874338981631210881d-01
     case(3,4)
        we=2.69266719309996355105d-01
        xe=4.33395394129247190794d-01
     case(5,6)
        we=2.19086362515982044000d-01
        xe=6.79409568299024406207d-01
     case(7,8)
        we=1.49451349150580593150d-01
        xe=8.65063366688984510759d-01
     case(9,10)
        we=6.66713443086881375920d-02
        xe=9.73906528517171720066d-01
     end select

     if (ibclr(e,0)==e) xe=-xe



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  case(11)
     if (e==1)  then
        we=  0.272925086777901d0
        xe=  0.000000000000000d0
     elseif (e==2)  then
        we=  0.262804544510247d0
        xe= -0.269543155952345d0
     elseif (e==3)  then
        we=  0.262804544510247d0
        xe=  0.269543155952345d0
     elseif (e==4)  then
        we=  0.233193764591991d0
        xe= -0.519096129206812d0
     elseif (e==5)  then
        we=  0.233193764591991d0
        xe=  0.519096129206812d0
     elseif (e==6)  then
        we=  0.186290210927734d0
        xe= -0.730152005574049d0
     elseif (e==7)  then
        we=  0.186290210927734d0
        xe=  0.730152005574049d0
     elseif (e==8)  then
        we=  0.125580369464905d0
        xe= -0.887062599768095d0
     elseif (e==9)  then
        we=  0.125580369464905d0
        xe=  0.887062599768095d0
     elseif (e==10)  then
        we=  5.566856711617370d-2
        xe= -0.978228658146057d0
     elseif (e==11)  then
        we=  5.566856711617370d-2
        xe=  0.978228658146057d0
     endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case(12)

     select case(e)
     case(1,2)
        we=2.49147045813402785006d-01
        xe=1.25233408511468915478d-01
     case(3,4)
        we=2.33492536538354808758d-01
        xe=3.67831498998180193757d-01
     case(5,6)
        we=2.03167426723065921743d-01
        xe=5.87317954286617447312d-01
     case(7,8)
        we=1.60078328543346226338d-01
        xe=7.69902674194304687059d-01
     case(9,10)
        we=1.06939325995318430960d-01
        xe=9.04117256370474856682d-01
     case(11,12)
        we=4.71753363865118271952d-02
        xe=9.81560634246719250712d-01
     end select

     if (ibclr(e,0)==e) xe=-xe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(13)
     if (e==1)  then
        we=  0.232551553230874d0
        xe=  0.000000000000000d0
     elseif (e==2)  then
        we=  0.226283180262897d0
        xe= -0.230458315955135d0
     elseif (e==3)  then
        we=  0.226283180262897d0
        xe=  0.230458315955135d0
     elseif (e==4)  then
        we=  0.207816047536889d0
        xe= -0.448492751036447d0
     elseif (e==5)  then
        we=  0.207816047536889d0
        xe=  0.448492751036447d0
     elseif (e==6)  then
        we=  0.178145980761946d0
        xe= -0.642349339440340d0
     elseif (e==7)  then
        we=  0.178145980761946d0
        xe=  0.642349339440340d0
     elseif (e==8)  then
        we=  0.138873510219787d0
        xe= -0.801578090733310d0
     elseif (e==9)  then
        we=  0.138873510219787d0
        xe=  0.801578090733310d0
     elseif (e==10)  then
        we=  9.212149983772849d-2
        xe= -0.917598399222978d0
     elseif (e==11)  then
        we=  9.212149983772849d-2
        xe=  0.917598399222978d0
     elseif (e==12)  then
        we=  4.048400476531590d-2
        xe= -0.984183054718588d0
     elseif (e==13)  then
        we=  4.048400476531590d-2
        xe=  0.984183054718588d0
     endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(14)
     if (e==1)  then
        we=  0.215263853463158d0
        xe= -0.108054948707344d0
     elseif (e==2)  then
        we=  0.215263853463158d0
        xe=  0.108054948707344d0
     elseif (e==3)  then
        we=  0.205198463721296d0
        xe= -0.319112368927890d0
     elseif (e==4)  then
        we=  0.205198463721296d0
        xe=  0.319112368927890d0
     elseif (e==5)  then
        we=  0.185538397477938d0
        xe= -0.515248636358154d0
     elseif (e==6)  then
        we=  0.185538397477938d0
        xe=  0.515248636358154d0
     elseif (e==7)  then
        we=  0.157203167158193d0
        xe= -0.687292904811685d0
     elseif (e==8)  then
        we=  0.157203167158193d0
        xe=  0.687292904811685d0
     elseif (e==9)  then
        we=  0.121518570687903d0
        xe= -0.827201315069765d0
     elseif (e==10)  then
        we=  0.121518570687903d0
        xe=  0.827201315069765d0
     elseif (e==11)  then
        we=  8.015808715976019d-2
        xe= -0.928434883663574d0
     elseif (e==12)  then
        we=  8.015808715976019d-2
        xe=  0.928434883663574d0
     elseif (e==13)  then
        we=  3.511946033175190d-2
        xe= -0.986283808696812d0
     elseif (e==14)  then
        we=  3.511946033175190d-2
        xe=  0.986283808696812d0
     endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(15)
     if (e==1)  then
        we=  0.202578241925561d0     
        xe=  0.000000000000000d0
     elseif (e==2)  then
        we=  0.198431485327112d0     
        xe= -0.201194093997434d0     
     elseif (e==3)  then
        we=  0.198431485327112d0     
        xe=  0.201194093997434d0     
     elseif (e==4)  then
        we=  0.186161000015562d0     
        xe= -0.394151347077563d0     
     elseif (e==5)  then
        we=  0.186161000015562d0     
        xe=  0.394151347077563d0     
     elseif (e==6)  then
        we=  0.166269205816994d0     
        xe= -0.570972172608539d0     
     elseif (e==7)  then
        we=  0.166269205816994d0     
        xe=  0.570972172608539d0     
     elseif (e==8)  then
        we=  0.139570677926154d0     
        xe= -0.724417731360170d0     
     elseif (e==9)  then
        we=  0.139570677926154d0     
        xe=  0.724417731360170d0     
     elseif (e==10)  then
        we=  0.107159220467172d0     
        xe= -0.848206583410427d0     
     elseif (e==11)  then
        we=  0.107159220467172d0     
        xe=  0.848206583410427d0     
     elseif (e==12)  then
        we=  7.036604748810810d-2
        xe= -0.937273392400706d0     
     elseif (e==13)  then
        we=  7.036604748810810d-2
        xe=  0.937273392400706d0     
     elseif (e==14)  then
        we=  3.075324199611730d-2
        xe= -0.987992518020485d0     
     elseif (e==15)  then
        we=  3.075324199611730d-2
        xe=  0.987992518020485d0
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(16)


     select case(e)
     case(1,2)
        we=1.89450610455068496287d-01
        xe=9.50125098376374401877d-02
     case(3,4)
        we=1.82603415044923588872d-01
        xe=2.81603550779258913231d-01
     case(5,6)
        we=1.69156519395002538183d-01
        xe=4.58016777657227386350d-01
     case(7,8)
        we=1.49595988816576732080d-01
        xe=6.17876244402643748452d-01
     case(9,10)
        we=1.24628971255533872056d-01
        xe=7.55404408355003033891d-01
     case(11,12)
        we=9.51585116824927848073d-02
        xe=8.65631202387831743866d-01
     case(13,14)
        we=6.22535239386478928628d-02
        xe=9.44575023073232576090d-01
     case(15,16)
        we=2.71524594117540948514d-02
        xe=9.89400934991649932601d-01
     end select

     if (ibclr(e,0)==e) xe=-xe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case(17)
     if (e==1)  then
        we=  0.179446470356207d0     
        xe=  0.000000000000000d0
     elseif (e==2)  then
        we=  0.176562705366993d0     
        xe= -0.178484181495848d0     
     elseif (e==3)  then
        we=  0.176562705366993d0     
        xe=  0.178484181495848d0     
     elseif (e==4)  then
        we=  0.168004102156450d0     
        xe= -0.351231763453876d0     
     elseif (e==5)  then
        we=  0.168004102156450d0     
        xe=  0.351231763453876d0     
     elseif (e==6)  then
        we=  0.154045761076810d0     
        xe= -0.512690537086477d0     
     elseif (e==7)  then
        we=  0.154045761076810d0     
        xe=  0.512690537086477d0     
     elseif (e==8)  then
        we=  0.135136368468526d0     
        xe= -0.657671159216691d0     
     elseif (e==9)  then
        we=  0.135136368468526d0     
        xe=  0.657671159216691d0     
     elseif (e==10)  then
        we=  0.111883847193404d0     
        xe= -0.781514003896801d0     
     elseif (e==11)  then
        we=  0.111883847193404d0     
        xe=  0.781514003896801d0     
     elseif (e==12)  then
        we=  8.503614831717921d-2
        xe= -0.880239153726986d0     
     elseif (e==13)  then
        we=  8.503614831717921d-2
        xe=  0.880239153726986d0     
     elseif (e==14)  then
        we=  5.545952937398720d-2
        xe= -0.950675521768768d0     
     elseif (e==15)  then
        we=  5.545952937398720d-2
        xe=  0.950675521768768d0     
     elseif (e==16)  then
        we=  2.414830286854790d-2
        xe= -0.990575475314417d0     
     elseif (e==17)  then
        we=  2.414830286854790d-2
        xe=  0.990575475314417d0  
     endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(18)
     if (e==1)  then
        we=  0.169142382963144d0     
        xe= -8.477501304173531d-2
     elseif (e==2)  then
        we=  0.169142382963144d0     
        xe=  8.477501304173531d-2
     elseif (e==3)  then
        we=  0.164276483745833d0    
        xe= -0.251886225691505d0    
     elseif (e==4)  then
        we=  0.164276483745833d0   
        xe=  0.251886225691505d0   
     elseif (e==5)  then
        we=  0.154684675126265d0   
        xe= -0.411751161462843d0   
     elseif (e==6)  then
        we=  0.154684675126265d0   
        xe=  0.411751161462843d0   
     elseif (e==7)  then
        we=  0.140642914670651d0   
        xe= -0.559770831073948d0   
     elseif (e==8)  then
        we=  0.140642914670651d0   
        xe=  0.559770831073948d0   
     elseif (e==9)  then
        we=  0.122555206711479d0   
        xe= -0.691687043060353d0   
     elseif (e==10)  then
        we=  0.122555206711479d0   
        xe=  0.691687043060353d0   
     elseif (e==11)  then
        we=  0.100942044106287d0   
        xe= -0.803704958972523d0   
     elseif (e==12)  then
        we=  0.100942044106287d0    
        xe=  0.803704958972523d0     
     elseif (e==13)  then
        we=  7.642573025488909d-2
        xe= -0.892602466497556d0     
     elseif (e==14)  then
        we=  7.642573025488909d-2
        xe=  0.892602466497556d0     
     elseif (e==15)  then
        we=  4.971454889496980d-2
        xe= -0.955823949571398d0     
     elseif (e==16)  then
        we=  4.971454889496980d-2
        xe=  0.955823949571398d0     
     elseif (e==17)  then
        we=  2.161601352648330d-2
        xe= -0.991565168420931d0     
     elseif (e==18)  then
        we=  2.161601352648330d-2
        xe=  0.991565168420931d0
     endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case(19)

     if (e==1)  then
        we=  0.161054449848784d0     
        xe=  0.000000000000000d0
     elseif (e==2)  then
        we=  0.158968843393954d0     
        xe= -0.160358645640225d0   
     elseif (e==3)  then
        we=  0.158968843393954d0   
        xe=  0.160358645640225d0   
     elseif (e==4)  then
        we=  0.152766042065860d0   
        xe= -0.316564099963630d0   
     elseif (e==5)  then
        we=  0.152766042065860d0   
        xe=  0.316564099963630d0   
     elseif (e==6)  then
        we=  0.142606702173607d0   
        xe= -0.464570741375961d0   
     elseif (e==7)  then
        we=  0.142606702173607d0   
        xe=  0.464570741375961d0   
     elseif (e==8)  then
        we=  0.128753962539336d0   
        xe= -0.600545304661681d0   
     elseif (e==9)  then
        we=  0.128753962539336d0   
        xe=  0.600545304661681d0   
     elseif (e==10)  then
        we=  0.111566645547334d0   
        xe= -0.720966177335229d0    
     elseif (e==11)  then
        we=  0.111566645547334d0     
        xe=  0.720966177335229d0     
     elseif (e==12)  then
        we=  9.149002162245000d-2
        xe= -0.822714656537143d0     
     elseif (e==13)  then
        we=  9.149002162245000d-2
        xe=  0.822714656537143d0     
     elseif (e==14)  then
        we=  6.904454273764120d-2
        xe= -0.903155903614818d0     
     elseif (e==15)  then
        we=  6.904454273764120d-2
        xe=  0.903155903614818d0     
     elseif (e==16)  then
        we=  4.481422676569960d-2
        xe= -0.960208152134830d0     
     elseif (e==17)  then
        we=  4.481422676569960d-2
        xe=  0.960208152134830d0     
     elseif (e==18)  then
        we=  1.946178822972650d-2
        xe= -0.992406843843584d0     
     elseif (e==19)  then
        we=  1.946178822972650d-2
        xe=  0.992406843843584d0 
     endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







  case(20)


     select case(e)
     case(1,2)
        we=1.52753387130725850699d-01
        xe=7.65265211334973337513d-02
     case(3,4)
        we=1.49172986472603746785d-01
        xe=2.27785851141645078076d-01
     case(5,6)
        we=1.42096109318382051326d-01
        xe=3.73706088715419560662d-01
     case(7,8)
        we=1.31688638449176626902d-01
        xe=5.10867001950827097985d-01
     case(9,10)
        we=1.18194531961518417310d-01
        xe=6.36053680726515025467d-01  
     case(11,12)
        we=1.01930119817240435039d-01
        xe=7.46331906460150792634d-01
     case(13,14)
        we=8.32767415767047487264d-02
        xe=8.39116971822218823420d-01 
     case(15,16)
        we=6.26720483341090635663d-02
        xe=9.12234428251325905857d-01
     case(17,18)
        we=4.06014298003869413320d-02 
        xe=9.63971927277913791287d-01
     case(19,20)
        we=1.76140071391521183115d-02
        xe=9.93128599185094924776d-01

     end select

     if (ibclr(e,0)==e) xe=-xe




  case(24)


     select case(e)
     case(1,2)
        we=1.27938195346752156976d-01
        xe=6.40568928626056260827d-02
     case(3,4)
        we=1.25837456346828296117d-01
        xe=1.91118867473616309153d-01
     case(5,6)
        we=1.21670472927803391202d-01
        xe=3.15042679696163374398d-01
     case(7,8)
        we=1.15505668053725601353d-01
        xe=4.33793507626045138478d-01
     case(9,10)
        we=1.07444270115965634785d-01
        xe=5.45421471388839535649d-01
     case(11,12)
        we=9.76186521041138882720d-02
        xe=6.48093651936975569268d-01
     case(13,14)
        we=8.61901615319532759152d-02
        xe=7.40124191578554364260d-01
     case(15,16)
        we=7.33464814110803057346d-02
        xe=8.20001985973902921981d-01
     case(17,18)
        we=5.92985849154367807461d-02
        xe=8.86415527004401034190d-01
     case(19,20)
        we=4.42774388174198061695d-02
        xe=9.38274552002732758539d-01
     case(21,22)
        we=2.85313886289336631809d-02
        xe=9.74728555971309498199d-01
     case(23,24)
        we=1.23412297999871995469d-02
        xe=9.95187219997021360195d-01
     end select

     if (ibclr(e,0)==e) xe=-xe


  case(32)


     select case(e)
     case(1,2)
        we=9.65400885147278005666d-02
        xe=4.83076656877383162364d-02
     case(3,4)
        we=9.56387200792748594185d-02
        xe=1.44471961582796493484d-01
     case(5,6)
        we=9.38443990808045656367d-02
        xe=2.39287362252137074544d-01
     case(7,8)
        we=9.11738786957638847129d-02
        xe=3.31868602282127649782d-01
     case(9,10)
        we=8.76520930044038111450d-02
        xe=4.21351276130635345353d-01
     case(11,12)
        we=8.33119242269467552223d-02
        xe=5.06899908932229390044d-01
     case(13,14)
        we=7.81938957870703064685d-02
        xe=5.87715757240762329066d-01
     case(15,16)
        we=7.23457941088485062287d-02
        xe=6.63044266930215200960d-01
     case(17,18)
        we=6.58222227763618468406d-02
        xe=7.32182118740289680412d-01
     case(19,20)
        we=5.86840934785355471448d-02
        xe=7.94483795967942406965d-01
     case(21,22)
        we=5.09980592623761761959d-02
        xe=8.49367613732569970160d-01
     case(23,24)
        we=4.28358980222266806557d-02
        xe=8.96321155766052123971d-01
     case(25,26)
        we=3.42738629130214331033d-02
        xe=9.34906075937739689159d-01
     case(27,28)
        we=2.53920653092620594561d-02
        xe=9.64762255587506430761d-01
     case(29,30)
        we=1.62743947309056706058d-02
        xe=9.85611511545268335400d-01
     case(31,32)
        we=7.01861000947009660028d-03
        xe=9.97263861849481563534d-01
     end select

     if (ibclr(e,0)==e) xe=-xe





  case(40)

     select case(e)
     case(1,2)
        we=7.75059479784248112668d-02
        xe=3.87724175060508219329d-02
     case(3,4)
        we=7.70398181642479655914d-02
        xe=1.16084070675255208481d-01
     case(5,6)
        we=7.61103619006262423723d-02
        xe=1.92697580701371099719d-01
     case(7,8)
        we=7.47231690579682641980d-02
        xe=2.68152185007253681152d-01
     case(9,10)
        we=7.28865823958040590609d-02
        xe=3.41994090825758473008d-01
     case(11,12)
        we=7.06116473912867796979d-02
        xe=4.13779204371605001525d-01
     case(13,14)
        we=6.79120458152339038265d-02
        xe=4.83075801686178712903d-01
     case(15,16)
        we=6.48040134566010380719d-02
        xe=5.49467125095128202056d-01
     case(17,18)
        we=6.13062424929289391679d-02    
        xe=6.12553889667980237972d-01
     case(19,20)
        we=5.74397690993915513665d-02
        xe=6.71956684614179548364d-01
     case(21,22)
        we=5.32278469839368243566d-02
        xe=7.27318255189927103277d-01
     case(23,24)
        we=4.86958076350722320604d-02
        xe=7.78305651426519387712d-01
     case(25,26)
        we=4.38709081856732719923d-02
        xe=8.24612230833311663197d-01
     case(27,28)
        we=3.87821679744720176413d-02
        xe=8.65959503212259503824d-01
     case(29,30)
        we=3.34601952825478473933d-02
        xe=9.02098806968874296732d-01
     case(31,32)
        we=2.79370069800234010984d-02
        xe=9.32812808278676533383d-01
     case(33,34)
        we=2.22458491941669572615d-02
        xe=9.57916819213791655824d-01
     case(35,36)
        we=1.64210583819078887131d-02
        xe=9.77259949983774262679d-01
     case(37,38)
        xe=9.90726238699457006464d-01
        we=1.04982845311528136146d-02
     case(39,40)
        we=4.52127709853319125846d-03
        xe=9.98237709710559200369d-01
     end select

     if (ibclr(e,0)==e) xe=-xe



  case(48)


     select case(e)
     case(1,2)
        we=6.47376968126839225006d-02
        xe=3.23801709628693620343d-02
     case(3,4)
        we=6.44661644359500822082d-02
        xe=9.70046992094626989322d-02
     case(5,6)
        we=6.39242385846481866207d-02
        xe=1.61222356068891718055d-01
     case(7,8)
        we=6.31141922862540256548d-02
        xe=2.24763790394689061224d-01
     case(9,10)
        we=6.20394231598926639029d-02
        xe=2.87362487355455576728d-01
     case(11,12)
        we=6.07044391658938800517d-02
        xe=3.48755886292160738148d-01
     case(13,14)
        we=5.91148396983956357477d-02
        xe=4.08686481990716729925d-01
     case(15,16)
        we=5.72772921004032157044d-02
        xe=4.66902904750958404535d-01
     case(17,18)
        we=5.51995036999841628676d-02
        xe=5.23160974722233033658d-01
     case(19,20)
        we=5.28901894851936670964d-02
        xe=5.77224726083972703838d-01
     case(21,22)
        we=5.03590355538544749590d-02
        xe=6.28867396776513624013d-01
     case(23,24)
        we=4.76166584924904748267d-02
        xe=6.77872379632663905208d-01
     case(25,26)
        we=4.46745608566942804201d-02
        xe=7.24034130923814654658d-01
     case(27,28)
        we=4.15450829434647492133d-02
        xe=7.67159032515740339276d-01
     case(29,30)
        we=3.82413510658307063158d-02
        xe=8.07066204029442627087d-01
     case(31,32)
        we=3.47772225647704388909d-02
        xe=8.43588261624393530704d-01
     case(33,34)
        we=3.11672278327980889025d-02
        xe=8.76572020274247885885d-01
     case(35,36)
        we=2.74265097083569482001d-02
        xe=9.05879136715569672805d-01
     case(37,38)
        we=2.35707608393243791410d-02
        xe=9.31386690706554333107d-01
     case(39,40)
        we=1.96161604573555278142d-02
        xe=9.52987703160430860724d-01
     case(41,42)
        we=1.55793157229438487279d-02
        xe=9.70591592546247250472d-01
     case(43,44)
        we=1.14772345792345394895d-02
        xe=9.84124583722826857765d-01
     case(45,46)
        we=7.32755390127626210220d-03
        xe=9.93530172266350757526d-01
     case(47,48)
        we=3.15334605230583863260d-03
        xe=9.98771007252426118580d-01

     end select

     if (ibclr(e,0)==e) xe=-xe






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  case(56)
     if (e==1)  then
        we=  5.557974630651440d-2
        xe= -2.779703528727540d-2
     elseif (e==2)  then
        we=  5.557974630651440d-2
        xe=  2.779703528727540d-2
     elseif (e==3)  then
        we=  5.540795250324510d-2
        xe= -8.330518682243540d-2
     elseif (e==4)  then
        we=  5.540795250324510d-2
        xe=  8.330518682243540d-2
     elseif (e==5)  then
        we=  5.506489590176240d-2
        xe= -0.138555846810376d0   
     elseif (e==6)  then
        we=  5.506489590176240d-2
        xe=  0.138555846810376d0   
     elseif (e==7)  then
        we=  5.455163687088940d-2
        xe= -0.193378238635275d0   
     elseif (e==8)  then
        we=  5.455163687088940d-2
        xe=  0.193378238635275d0   
     elseif (e==9)  then
        we=  5.386976186571450d-2
        xe= -0.247602909434337d0   
     elseif (e==10)  then
        we=  5.386976186571450d-2
        xe=  0.247602909434337d0   
     elseif (e==11)  then
        we=  5.302137852401080d-2
        xe= -0.301062253867221d0   
     elseif (e==12)  then
        we=  5.302137852401080d-2
        xe=  0.301062253867221d0   
     elseif (e==13)  then
        we=  5.200910915174140d-2
        xe= -0.353591032174955d0   
     elseif (e==14)  then
        we=  5.200910915174140d-2
        xe=  0.353591032174955d0   
     elseif (e==15)  then
        we=  5.083608261779850d-2
        xe= -0.405026880927091d0   
     elseif (e==16)  then
        we=  5.083608261779850d-2
        xe=  0.405026880927091d0   
     elseif (e==17)  then
        we=  4.950592468304760d-2
        xe= -0.455210814878460d0   
     elseif (e==18)  then
        we=  4.950592468304760d-2
        xe=  0.455210814878460d0   
     elseif (e==19)  then
        we=  4.802274679360030d-2
        xe= -0.503987718384382d0   
     elseif (e==20)  then
        we=  4.802274679360030d-2
        xe=  0.503987718384382d0   
     elseif (e==21)  then
        we=  4.639113337300190d-2
        xe= -0.551206824855535d0   
     elseif (e==22)  then
        we=  4.639113337300190d-2
        xe=  0.551206824855535d0   
     elseif (e==23)  then
        we=  4.461612765269230d-2
        xe= -0.596722182770663d0   
     elseif (e==24)  then
        we=  4.461612765269230d-2
        xe=  0.596722182770663d0   
     elseif (e==25)  then
        we=  4.270321608466710d-2
        xe= -0.640393106807007d0   
     elseif (e==26)  then
        we=  4.270321608466710d-2
        xe=  0.640393106807007d0   
     elseif (e==27)  then
        we=  4.065831138474450d-2
        xe= -0.682084612694470d0   
     elseif (e==28)  then
        we=  4.065831138474450d-2
        xe=  0.682084612694470d0   
     elseif (e==29)  then
        we=  3.848773425924770d-2
        xe= -0.721667834450188d0   
     elseif (e==30)  then
        we=  3.848773425924770d-2
        xe=  0.721667834450188d0   
     elseif (e==31)  then
        we=  3.619819387231520d-2
        xe= -0.759020422705129d0   
     elseif (e==32)  then
        we=  3.619819387231520d-2
        xe=  0.759020422705129d0   
     elseif (e==33)  then
        we=  3.379676711561180d-2
        xe= -0.794026922893866d0   
     elseif (e==34)  then
        we=  3.379676711561180d-2
        xe=  0.794026922893866d0   
     elseif (e==35)  then
        we=  3.129087674731040d-2
        xe= -0.826579132142882d0   
     elseif (e==36)  then
        we=  3.129087674731040d-2
        xe=  0.826579132142882d0   
     elseif (e==37)  then
        we=  2.868826847382270d-2
        xe= -0.856576433762749d0   
     elseif (e==38)  then
        we=  2.868826847382270d-2
        xe=  0.856576433762749d0   
     elseif (e==39)  then
        we=  2.599698705839200d-2
        xe= -0.883926108327828d0   
     elseif (e==40)  then
        we=  2.599698705839200d-2
        xe=  0.883926108327828d0   
     elseif (e==41)  then
        we=  2.322535156256530d-2
        xe= -0.908543620420656d0   
     elseif (e==42)  then
        we=  2.322535156256530d-2
        xe=  0.908543620420656d0   
     elseif (e==43)  then
        we=  2.038192988240260d-2
        xe= -0.930352880247496d0   
     elseif (e==44)  then
        we=  2.038192988240260d-2
        xe=  0.930352880247496d0   
     elseif (e==45)  then
        we=  1.747551291140090d-2
        xe= -0.949286479561963d0   
     elseif (e==46)  then
        we=  1.747551291140090d-2
        xe=  0.949286479561963d0   
     elseif (e==47)  then
        we=  1.451508927802150d-2
        xe= -0.965285901905490d0   
     elseif (e==48)  then
        we=  1.451508927802150d-2
        xe=  0.965285901905490d0   
     elseif (e==49)  then
        we=  1.150982434038340d-2
        xe= -0.978301709140256d0   
     elseif (e==50)  then
        we=  1.150982434038340d-2
        xe=  0.978301709140256d0   
     elseif (e==51)  then
        we=  8.469063163307901d-3
        xe= -0.988293715540162d0   
     elseif (e==52)  then
        we=  8.469063163307901d-3
        xe=  0.988293715540162d0     
     elseif (e==53)  then
        we=  5.402522246015300d-3
        xe= -0.995231226081070d0     
     elseif (e==54)  then
        we=  5.402522246015300d-3
        xe=  0.995231226081070d0     
     elseif (e==55)  then
        we=  2.323855375773200d-3
        xe= -0.999094343801466d0   
     elseif (e==56)  then
        we=  2.323855375773200d-3
        xe=  0.999094343801466d0     
     end if

  end select!!!!!!!!!!!!!!!!!!!!!!! nbe


end subroutine dset_gauss_legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine sset_zolotarev(nbe,e,xe,we)
  !  Purpose 
  !  =======
  !  Zolotarev Quadrature ROUTINE- Return weight and coordinate of the e^th node from nbe total # nodes
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  nbe        (input)        INTEGER: total # nodes for Zolotarev quadrature
  !  e          (input)        INTEGER: e^th node
  !  xe         (output)       COMPLEX SINGLE PRECISION:  Zolotarev coordinate for node e 
  !  we         (output)       COMPLEX SINGLE PRECISION:  Zolotarev weight for node e
  !=====================================================================
  ! Stefan Guettel, Eric Polizzi 2013-2015
  !=====================================================================
  implicit none
  integer :: nbe,e
  complex :: xe,we
!!!!!!!!!!

  complex(kind=(kind(1.0d0))) :: zxe,zwe

  call dset_zolotarev(nbe,e,zxe,zwe)

  xe=cmplx(zxe)
  we=cmplx(zwe)

end subroutine sset_zolotarev








subroutine dset_zolotarev(nbe,e,xe,we)
  !  Purpose 
  !  =======
  !  Zolotarev Quadrature ROUTINE- Return weight and coordinate of the e^th node from nbe total # nodes
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  nbe        (input)        INTEGER: total # nodes for Zolotarev quadrature
  !  e          (input)        INTEGER: e^th node
  !  xe         (output)       COMPLEX DOUBLE PRECISION:  Zolotarev coordinate for node e 
  !  we         (output)       COMPLEX DOUBLE PRECISION:  Zolotarev weight for node e
  !=====================================================================
  ! Stefan Guettel, Eric Polizzi 2013-2015
  !=====================================================================
  implicit none
  integer :: nbe,e
  complex(kind=(kind(1.0d0))) :: xe,we
!!!!!!!!!!



  select case(nbe)


  case(1)
     select case(e)
        !------------------------   n =  1, S = 1.0d06, rate = 9.92d-1  ------------------------
     case(0)
        we = (-4.9800399400799011d-1)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (9.9800399400799011d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     end select


  case(2)
     select case(e)
        !------------------------   n =  2, S = 1.0d06, rate = 7.18d-1  ------------------------
     case(0)
        we = (4.1805096443248230d-1)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-4.0933604666346268d-2)*(1.0d0,0.0d0) + (1.8306055366585177d-3)*(0.0d0,1.0d0) 
        xe = (-9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(2)
        we = (4.0933604666346268d-2)*(1.0d0,0.0d0) + (1.8306055366585177d-3)*(0.0d0,1.0d0) 
        xe = (9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     end select


  case(3)
     select case(e)
        !------------------------   n =  3, S = 1.0d06, rate = 3.58d-1  ------------------------
     case(0)
        we = (-2.6356075833756432d-1)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-9.4403480868925395d-3)*(1.0d0,0.0d0) + (1.1819628351771121d-4)*(0.0d0,1.0d0) 
        xe = (-9.9992162986865651d-1)*(1.0d0,0.0d0) + (1.2519349855705481d-2)*(0.0d0,1.0d0) 
     case(2)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (7.4467858236516826d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(3)
        we = (9.4403480868925395d-3)*(1.0d0,0.0d0) + (1.1819628351771121d-4)*(0.0d0,1.0d0) 
        xe = (9.9992162986865651d-1)*(1.0d0,0.0d0) + (1.2519349855705481d-2)*(0.0d0,1.0d0) 
     end select


  case(4)
     select case(e)
        !------------------------   n =  4, S = 1.0d06, rate = 1.71d-1  ------------------------
     case(0)
        we = (1.4575533545956809d-1)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-4.0449611415998383d-3)*(1.0d0,0.0d0) + (2.6445705300228226d-5)*(0.0d0,1.0d0) 
        xe = (-9.9997862836826079d-1)*(1.0d0,0.0d0) + (6.5377983092055223d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.6550649051159086d-1)*(1.0d0,0.0d0) + (5.0629517778791600d-2)*(0.0d0,1.0d0) 
        xe = (-9.5625772508631735d-1)*(1.0d0,0.0d0) + (2.9252549156054936d-1)*(0.0d0,1.0d0) 
     case(3)
        we = (1.6550649051159086d-1)*(1.0d0,0.0d0) + (5.0629517778791600d-2)*(0.0d0,1.0d0) 
        xe = (9.5625772508631735d-1)*(1.0d0,0.0d0) + (2.9252549156054936d-1)*(0.0d0,1.0d0) 
     case(4)
        we = (4.0449611415998383d-3)*(1.0d0,0.0d0) + (2.6445705300228226d-5)*(0.0d0,1.0d0) 
        xe = (9.9997862836826079d-1)*(1.0d0,0.0d0) + (6.5377983092055223d-3)*(0.0d0,1.0d0) 
     end select


  case(5)
     select case(e)
        !------------------------   n =  5, S = 1.0d06, rate = 8.39d-2  ------------------------
     case(0)
        we = (-7.7374834133583259d-2)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-2.3046399755639315d-3)*(1.0d0,0.0d0) + (1.0035321417280235d-5)*(0.0d0,1.0d0) 
        xe = (-9.9999051974059949d-1)*(1.0d0,0.0d0) + (4.3543574641737573d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-4.5688379331906018d-2)*(1.0d0,0.0d0) + (4.3789522408392290d-3)*(0.0d0,1.0d0) 
        xe = (-9.9543837740544383d-1)*(1.0d0,0.0d0) + (9.5406691528514928d-2)*(0.0d0,1.0d0) 
     case(3)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (4.8097001541667800d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(4)
        we = (4.5688379331906018d-2)*(1.0d0,0.0d0) + (4.3789522408392290d-3)*(0.0d0,1.0d0) 
        xe = (9.9543837740544383d-1)*(1.0d0,0.0d0) + (9.5406691528514928d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (2.3046399755639315d-3)*(1.0d0,0.0d0) + (1.0035321417280235d-5)*(0.0d0,1.0d0) 
        xe = (9.9999051974059949d-1)*(1.0d0,0.0d0) + (4.3543574641737573d-3)*(0.0d0,1.0d0) 
     end select


  case(6)
     select case(e)
        !------------------------   n =  6, S = 1.0d06, rate = 4.23d-2  ------------------------
     case(0)
        we = (4.0601929796073799d-2)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-1.5423558061466516d-3)*(1.0d0,0.0d0) + (5.0401387941688375d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999466072397813d-1)*(1.0d0,0.0d0) + (3.2678010245045879d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.7985701250813335d-2)*(1.0d0,0.0d0) + (8.0434461022422783d-4)*(0.0d0,1.0d0) 
        xe = (-9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.7924652624997453d-1)*(1.0d0,0.0d0) + (1.0970398051397839d-1)*(0.0d0,1.0d0) 
        xe = (-8.5293349191434553d-1)*(1.0d0,0.0d0) + (5.2201959577280344d-1)*(0.0d0,1.0d0) 
     case(4)
        we = (1.7924652624997453d-1)*(1.0d0,0.0d0) + (1.0970398051397839d-1)*(0.0d0,1.0d0) 
        xe = (8.5293349191434553d-1)*(1.0d0,0.0d0) + (5.2201959577280344d-1)*(0.0d0,1.0d0) 
     case(5)
        we = (1.7985701250813335d-2)*(1.0d0,0.0d0) + (8.0434461022422783d-4)*(0.0d0,1.0d0) 
        xe = (9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (1.5423558061466516d-3)*(1.0d0,0.0d0) + (5.0401387941688375d-6)*(0.0d0,1.0d0) 
        xe = (9.9999466072397813d-1)*(1.0d0,0.0d0) + (3.2678010245045879d-3)*(0.0d0,1.0d0) 
     end select


  case(7)
     select case(e)
        !------------------------   n =  7, S = 1.0d06, rate = 2.17d-2  ------------------------

     case(0)
        we = (-2.1237706090261987d-2)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-1.1399193084610077d-3)*(1.0d0,0.0d0) + (2.9915174055686861d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999655648000429d-1)*(1.0d0,0.0d0) + (2.6243147931738105d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-8.9861107855352777d-3)*(1.0d0,0.0d0) + (2.3319989914114602d-4)*(0.0d0,1.0d0) 
        xe = (-9.9966343892137466d-1)*(1.0d0,0.0d0) + (2.5942414765997054d-2)*(0.0d0,1.0d0) 
     case(3)
        we = (-7.5755925208015523d-2)*(1.0d0,0.0d0) + (1.7497105287481524d-2)*(0.0d0,1.0d0) 
        xe = (-9.7434899659958563d-1)*(1.0d0,0.0d0) + (2.2504229119296795d-1)*(0.0d0,1.0d0) 
     case(4)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (3.4547899051275072d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(5)
        we = (7.5755925208015523d-2)*(1.0d0,0.0d0) + (1.7497105287481524d-2)*(0.0d0,1.0d0) 
        xe = (9.7434899659958563d-1)*(1.0d0,0.0d0) + (2.2504229119296795d-1)*(0.0d0,1.0d0) 
     case(6)
        we = (8.9861107855352777d-3)*(1.0d0,0.0d0) + (2.3319989914114602d-4)*(0.0d0,1.0d0) 
        xe = (9.9966343892137466d-1)*(1.0d0,0.0d0) + (2.5942414765997054d-2)*(0.0d0,1.0d0) 
     case(7)
        we = (1.1399193084610077d-3)*(1.0d0,0.0d0) + (2.9915174055686861d-6)*(0.0d0,1.0d0) 
        xe = (9.9999655648000429d-1)*(1.0d0,0.0d0) + (2.6243147931738105d-3)*(0.0d0,1.0d0) 
     end select


  case(8)
     select case(e)
        !------------------------   n =  8, S = 1.0d06, rate = 1.12d-2  ------------------------
     case(0)
        we = (1.1099137041258145d-2)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-8.9892014626439772d-4)*(1.0d0,0.0d0) + (1.9770010320296091d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-5.2457912271928649d-3)*(1.0d0,0.0d0) + (9.0422169329207065d-5)*(0.0d0,1.0d0) 
        xe = (-9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(3)
        we = (-3.4625385252140740d-2)*(1.0d0,0.0d0) + (4.0175404307143140d-3)*(0.0d0,1.0d0) 
        xe = (-9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(4)
        we = (-1.5051737271560608d-1)*(1.0d0,0.0d0) + (1.3687697801045523d-1)*(0.0d0,1.0d0) 
        xe = (-7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(5)
        we = (1.5051737271560608d-1)*(1.0d0,0.0d0) + (1.3687697801045523d-1)*(0.0d0,1.0d0) 
        xe = (7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(6)
        we = (3.4625385252140740d-2)*(1.0d0,0.0d0) + (4.0175404307143140d-3)*(0.0d0,1.0d0) 
        xe = (9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (5.2457912271928649d-3)*(1.0d0,0.0d0) + (9.0422169329207065d-5)*(0.0d0,1.0d0) 
        xe = (9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(8)
        we = (8.9892014626439772d-4)*(1.0d0,0.0d0) + (1.9770010320296091d-6)*(0.0d0,1.0d0) 
        xe = (9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     end select


  case(9)
     select case(e)
        !------------------------   n =  9, S = 1.0d06, rate = 5.83d-3  ------------------------
     case(0)
        we = (-5.7991882972451281d-3)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-7.4104052419255275d-4)*(1.0d0,0.0d0) + (1.4058775564853451d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999820038372844d-1)*(1.0d0,0.0d0) + (1.8971634891048340d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-3.4078286502537015d-3)*(1.0d0,0.0d0) + (4.2667142950419761d-5)*(0.0d0,1.0d0) 
        xe = (-9.9992162986865651d-1)*(1.0d0,0.0d0) + (1.2519349855705481d-2)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.8278566233392055d-2)*(1.0d0,0.0d0) + (1.2481296809871830d-3)*(0.0d0,1.0d0) 
        xe = (-9.9767678352825928d-1)*(1.0d0,0.0d0) + (6.8125146669250583d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-8.9685989000547653d-2)*(1.0d0,0.0d0) + (3.4297391012017266d-2)*(0.0d0,1.0d0) 
        xe = (-9.3403206804331040d-1)*(1.0d0,0.0d0) + (3.5718915978335170d-1)*(0.0d0,1.0d0) 
     case(5)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (2.6881816060764269d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(6)
        we = (8.9685989000547653d-2)*(1.0d0,0.0d0) + (3.4297391012017266d-2)*(0.0d0,1.0d0) 
        xe = (9.3403206804331040d-1)*(1.0d0,0.0d0) + (3.5718915978335170d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (1.8278566233392055d-2)*(1.0d0,0.0d0) + (1.2481296809871830d-3)*(0.0d0,1.0d0) 
        xe = (9.9767678352825928d-1)*(1.0d0,0.0d0) + (6.8125146669250583d-2)*(0.0d0,1.0d0) 
     case(8)
        we = (3.4078286502537015d-3)*(1.0d0,0.0d0) + (4.2667142950419761d-5)*(0.0d0,1.0d0) 
        xe = (9.9992162986865651d-1)*(1.0d0,0.0d0) + (1.2519349855705481d-2)*(0.0d0,1.0d0) 
     case(9)
        we = (7.4104052419255275d-4)*(1.0d0,0.0d0) + (1.4058775564853451d-6)*(0.0d0,1.0d0) 
        xe = (9.9999820038372844d-1)*(1.0d0,0.0d0) + (1.8971634891048340d-3)*(0.0d0,1.0d0) 
     end select


  case(10)
     select case(e)
        !------------------------   n = 10, S = 1.0d06, rate = 3.04d-3  ------------------------
     case(0)
        we = (3.0298206493477586d-3)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-6.3052218358594105d-4)*(1.0d0,0.0d0) + (1.0535055580202612d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999860413950714d-1)*(1.0d0,0.0d0) + (1.6708438099384962d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-2.3906025222613010d-3)*(1.0d0,0.0d0) + (2.3134572353828454d-5)*(0.0d0,1.0d0) 
        xe = (-9.9995317824297625d-1)*(1.0d0,0.0d0) + (9.6768446184941678d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.0809177648659182d-2)*(1.0d0,0.0d0) + (4.8340087836508665d-4)*(0.0d0,1.0d0) 
        xe = (-9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-4.7956018043409793d-2)*(1.0d0,0.0d0) + (9.9110250490149158d-3)*(0.0d0,1.0d0) 
        xe = (-9.7930459779150147d-1)*(1.0d0,0.0d0) + (2.0239195820097600d-1)*(0.0d0,1.0d0) 
     case(5)
        we = (-1.1904157484648165d-1)*(1.0d0,0.0d0) + (1.4249242081357993d-1)*(0.0d0,1.0d0) 
        xe = (-6.4113075594177837d-1)*(1.0d0,0.0d0) + (7.6743166066140622d-1)*(0.0d0,1.0d0) 
     case(6)
        we = (1.1904157484648165d-1)*(1.0d0,0.0d0) + (1.4249242081357993d-1)*(0.0d0,1.0d0) 
        xe = (6.4113075594177837d-1)*(1.0d0,0.0d0) + (7.6743166066140622d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (4.7956018043409793d-2)*(1.0d0,0.0d0) + (9.9110250490149158d-3)*(0.0d0,1.0d0) 
        xe = (9.7930459779150147d-1)*(1.0d0,0.0d0) + (2.0239195820097600d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (1.0809177648659182d-2)*(1.0d0,0.0d0) + (4.8340087836508665d-4)*(0.0d0,1.0d0) 
        xe = (9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(9)
        we = (2.3906025222613010d-3)*(1.0d0,0.0d0) + (2.3134572353828454d-5)*(0.0d0,1.0d0) 
        xe = (9.9995317824297625d-1)*(1.0d0,0.0d0) + (9.6768446184941678d-3)*(0.0d0,1.0d0) 
     case(10)
        we = (6.3052218358594105d-4)*(1.0d0,0.0d0) + (1.0535055580202612d-6)*(0.0d0,1.0d0) 
        xe = (9.9999860413950714d-1)*(1.0d0,0.0d0) + (1.6708438099384962d-3)*(0.0d0,1.0d0) 
     end select


  case(11)
     select case(e)
        !------------------------   n = 11, S = 1.0d06, rate = 1.59d-3  ------------------------
     case(0)
        we = (-1.5829197980627985d-3)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-5.4916092156993054d-4)*(1.0d0,0.0d0) + (8.2078370994785510d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999888306607199d-1)*(1.0d0,0.0d0) + (1.4946125278547568d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.7758360089787139d-3)*(1.0d0,0.0d0) + (1.3891886287560548d-5)*(0.0d0,1.0d0) 
        xe = (-9.9996940384898791d-1)*(1.0d0,0.0d0) + (7.8224910290555765d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-6.9646144107369046d-3)*(1.0d0,0.0d0) + (2.2031308203370500d-4)*(0.0d0,1.0d0) 
        xe = (-9.9950004532269476d-1)*(1.0d0,0.0d0) + (3.1617390783098893d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-2.7405037741489228d-2)*(1.0d0,0.0d0) + (3.4688676881267449d-3)*(0.0d0,1.0d0) 
        xe = (-9.9208403387564725d-1)*(1.0d0,0.0d0) + (1.2557575295025525d-1)*(0.0d0,1.0d0) 
     case(5)
        we = (-9.1577068203742937d-2)*(1.0d0,0.0d0) + (4.9080286015786748d-2)*(0.0d0,1.0d0) 
        xe = (-8.8139525357848125d-1)*(1.0d0,0.0d0) + (4.7237951582316168d-1)*(0.0d0,1.0d0) 
     case(6)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (2.1994897763355334d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(7)
        we = (9.1577068203742937d-2)*(1.0d0,0.0d0) + (4.9080286015786748d-2)*(0.0d0,1.0d0) 
        xe = (8.8139525357848125d-1)*(1.0d0,0.0d0) + (4.7237951582316168d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (2.7405037741489228d-2)*(1.0d0,0.0d0) + (3.4688676881267449d-3)*(0.0d0,1.0d0) 
        xe = (9.9208403387564725d-1)*(1.0d0,0.0d0) + (1.2557575295025525d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (6.9646144107369046d-3)*(1.0d0,0.0d0) + (2.2031308203370500d-4)*(0.0d0,1.0d0) 
        xe = (9.9950004532269476d-1)*(1.0d0,0.0d0) + (3.1617390783098893d-2)*(0.0d0,1.0d0) 
     case(10)
        we = (1.7758360089787139d-3)*(1.0d0,0.0d0) + (1.3891886287560548d-5)*(0.0d0,1.0d0) 
        xe = (9.9996940384898791d-1)*(1.0d0,0.0d0) + (7.8224910290555765d-3)*(0.0d0,1.0d0) 
     case(11)
        we = (5.4916092156993054d-4)*(1.0d0,0.0d0) + (8.2078370994785510d-7)*(0.0d0,1.0d0) 
        xe = (9.9999888306607199d-1)*(1.0d0,0.0d0) + (1.4946125278547568d-3)*(0.0d0,1.0d0) 
     end select


  case(12)
     select case(e)
        !------------------------   n = 12, S = 1.0d06, rate = 8.28d-4  ------------------------
     case(0)
        we = (8.2698721091084559d-4)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-4.8687164580549528d-4)*(1.0d0,0.0d0) + (6.5885540028861806d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999908436863361d-1)*(1.0d0,0.0d0) + (1.3532412550538205d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.3784196581295379d-3)*(1.0d0,0.0d0) + (9.0120223119170229d-6)*(0.0d0,1.0d0) 
        xe = (-9.9997862836826079d-1)*(1.0d0,0.0d0) + (6.5377983092055223d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-4.7923676064111223d-3)*(1.0d0,0.0d0) + (1.1357135190625232d-4)*(0.0d0,1.0d0) 
        xe = (-9.9971931159972960d-1)*(1.0d0,0.0d0) + (2.3691728821737014d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-1.6899924452230135d-2)*(1.0d0,0.0d0) + (1.4262499855059792d-3)*(0.0d0,1.0d0) 
        xe = (-9.9645774817202859d-1)*(1.0d0,0.0d0) + (8.4094923199501848d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-5.6400393497717298d-2)*(1.0d0,0.0d0) + (1.7253249201870668d-2)*(0.0d0,1.0d0) 
        xe = (-9.5625772508631735d-1)*(1.0d0,0.0d0) + (2.9252549156054936d-1)*(0.0d0,1.0d0) 
     case(6)
        we = (-9.3578819652718012d-2)*(1.0d0,0.0d0) + (1.3830296710345763d-1)*(0.0d0,1.0d0) 
        xe = (-5.6039535450830469d-1)*(1.0d0,0.0d0) + (8.2822523907781964d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (9.3578819652718012d-2)*(1.0d0,0.0d0) + (1.3830296710345763d-1)*(0.0d0,1.0d0) 
        xe = (5.6039535450830469d-1)*(1.0d0,0.0d0) + (8.2822523907781964d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (5.6400393497717298d-2)*(1.0d0,0.0d0) + (1.7253249201870668d-2)*(0.0d0,1.0d0) 
        xe = (9.5625772508631735d-1)*(1.0d0,0.0d0) + (2.9252549156054936d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (1.6899924452230135d-2)*(1.0d0,0.0d0) + (1.4262499855059792d-3)*(0.0d0,1.0d0) 
        xe = (9.9645774817202859d-1)*(1.0d0,0.0d0) + (8.4094923199501848d-2)*(0.0d0,1.0d0) 
     case(10)
        we = (4.7923676064111223d-3)*(1.0d0,0.0d0) + (1.1357135190625232d-4)*(0.0d0,1.0d0) 
        xe = (9.9971931159972960d-1)*(1.0d0,0.0d0) + (2.3691728821737014d-2)*(0.0d0,1.0d0) 
     case(11)
        we = (1.3784196581295379d-3)*(1.0d0,0.0d0) + (9.0120223119170229d-6)*(0.0d0,1.0d0) 
        xe = (9.9997862836826079d-1)*(1.0d0,0.0d0) + (6.5377983092055223d-3)*(0.0d0,1.0d0) 
     case(12)
        we = (4.8687164580549528d-4)*(1.0d0,0.0d0) + (6.5885540028861806d-7)*(0.0d0,1.0d0) 
        xe = (9.9999908436863361d-1)*(1.0d0,0.0d0) + (1.3532412550538205d-3)*(0.0d0,1.0d0) 
     end select



  case(13)
     select case(e)
        !------------------------   n = 13, S = 1.0d06, rate = 4.32d-4  ------------------------
     case(0)
        we = (-4.3205406996127405d-4)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)       
     case(1)
        we = (-4.3767780889597465d-4)*(1.0d0,0.0d0) + (5.4147060148906031d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999923473788310d-1)*(1.0d0,0.0d0) + (1.2371433417836721d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.1075297345821931d-3)*(1.0d0,0.0d0) + (6.2076540794421900d-6)*(0.0d0,1.0d0) 
        xe = (-9.9998429261104993d-1)*(1.0d0,0.0d0) + (5.6048667404373113d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-3.4716788620964018d-3)*(1.0d0,0.0d0) + (6.4409116887187564d-5)*(0.0d0,1.0d0) 
        xe = (-9.9982794254032570d-1)*(1.0d0,0.0d0) + (1.8549536802280383d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-1.1122499655942791d-2)*(1.0d0,0.0d0) + (6.6672796922738465d-4)*(0.0d0,1.0d0) 
        xe = (-9.9820818130463496d-1)*(1.0d0,0.0d0) + (5.9836667491538521d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-3.4914922416108364d-2)*(1.0d0,0.0d0) + (6.7980844778368485d-3)*(0.0d0,1.0d0) 
        xe = (-9.8156757516467719d-1)*(1.0d0,0.0d0) + (1.9111539808538686d-1)*(0.0d0,1.0d0) 
     case(6)
        we = (-8.6881729971417979d-2)*(1.0d0,0.0d0) + (5.9724667568325514d-2)*(0.0d0,1.0d0) 
        xe = (-8.2407080158656099d-1)*(1.0d0,0.0d0) + (5.6648681712153093d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.8611110496879130d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(8)
        we = (8.6881729971417979d-2)*(1.0d0,0.0d0) + (5.9724667568325514d-2)*(0.0d0,1.0d0) 
        xe = (8.2407080158656099d-1)*(1.0d0,0.0d0) + (5.6648681712153093d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (3.4914922416108364d-2)*(1.0d0,0.0d0) + (6.7980844778368485d-3)*(0.0d0,1.0d0) 
        xe = (9.8156757516467719d-1)*(1.0d0,0.0d0) + (1.9111539808538686d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (1.1122499655942791d-2)*(1.0d0,0.0d0) + (6.6672796922738465d-4)*(0.0d0,1.0d0) 
        xe = (9.9820818130463496d-1)*(1.0d0,0.0d0) + (5.9836667491538521d-2)*(0.0d0,1.0d0) 
     case(11)
        we = (3.4716788620964018d-3)*(1.0d0,0.0d0) + (6.4409116887187564d-5)*(0.0d0,1.0d0) 
        xe = (9.9982794254032570d-1)*(1.0d0,0.0d0) + (1.8549536802280383d-2)*(0.0d0,1.0d0) 
     case(12)
        we = (1.1075297345821931d-3)*(1.0d0,0.0d0) + (6.2076540794421900d-6)*(0.0d0,1.0d0) 
        xe = (9.9998429261104993d-1)*(1.0d0,0.0d0) + (5.6048667404373113d-3)*(0.0d0,1.0d0) 
     case(13)
        we = (4.3767780889597465d-4)*(1.0d0,0.0d0) + (5.4147060148906031d-7)*(0.0d0,1.0d0) 
        xe = (9.9999923473788310d-1)*(1.0d0,0.0d0) + (1.2371433417836721d-3)*(0.0d0,1.0d0) 
     end select



  case(14)
     select case(e)
        !------------------------   n = 14, S = 1.0d06, rate = 2.26d-4  ------------------------
     case(0)
        we = (2.2572374690005281d-4)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0) 
     case(1)
        we = (-3.9783815941610621d-4)*(1.0d0,0.0d0) + (4.5352749686171894d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999935022357100d-1)*(1.0d0,0.0d0) + (1.1399791383407472d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-9.1483445273240691d-4)*(1.0d0,0.0d0) + (4.4839989288394089d-6)*(0.0d0,1.0d0) 
        xe = (-9.9998798819824153d-1)*(1.0d0,0.0d0) + (4.9013731987690712d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-2.6203167326490271d-3)*(1.0d0,0.0d0) + (3.9392030922682924d-5)*(0.0d0,1.0d0) 
        xe = (-9.9988701896134347d-1)*(1.0d0,0.0d0) + (1.5031610445923209d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-7.7209116615673213d-3)*(1.0d0,0.0d0) + (3.4528949382600389d-4)*(0.0d0,1.0d0) 
        xe = (-9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-2.2593959355968814d-2)*(1.0d0,0.0d0) + (3.0058501258880678d-3)*(0.0d0,1.0d0) 
        xe = (-9.9126623989398643d-1)*(1.0d0,0.0d0) + (1.3187585695053436d-1)*(0.0d0,1.0d0) 
     case(6)
        we = (-6.0452583497664004d-2)*(1.0d0,0.0d0) + (2.4667289940949678d-2)*(0.0d0,1.0d0) 
        xe = (-9.2588640299049352d-1)*(1.0d0,0.0d0) + (3.7780202323085216d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (-7.4349151462284435d-2)*(1.0d0,0.0d0) + (1.3043927963654731d-1)*(0.0d0,1.0d0) 
        xe = (-4.9519682077916088d-1)*(1.0d0,0.0d0) + (8.6878081740460389d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (7.4349151462284435d-2)*(1.0d0,0.0d0) + (1.3043927963654731d-1)*(0.0d0,1.0d0) 
        xe = (4.9519682077916088d-1)*(1.0d0,0.0d0) + (8.6878081740460389d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (6.0452583497664004d-2)*(1.0d0,0.0d0) + (2.4667289940949678d-2)*(0.0d0,1.0d0) 
        xe = (9.2588640299049352d-1)*(1.0d0,0.0d0) + (3.7780202323085216d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (2.2593959355968814d-2)*(1.0d0,0.0d0) + (3.0058501258880678d-3)*(0.0d0,1.0d0) 
        xe = (9.9126623989398643d-1)*(1.0d0,0.0d0) + (1.3187585695053436d-1)*(0.0d0,1.0d0) 
     case(11)
        we = (7.7209116615673213d-3)*(1.0d0,0.0d0) + (3.4528949382600389d-4)*(0.0d0,1.0d0) 
        xe = (9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(12)
        we = (2.6203167326490271d-3)*(1.0d0,0.0d0) + (3.9392030922682924d-5)*(0.0d0,1.0d0) 
        xe = (9.9988701896134347d-1)*(1.0d0,0.0d0) + (1.5031610445923209d-2)*(0.0d0,1.0d0) 
     case(13)
        we = (9.1483445273240691d-4)*(1.0d0,0.0d0) + (4.4839989288394089d-6)*(0.0d0,1.0d0) 
        xe = (9.9998798819824153d-1)*(1.0d0,0.0d0) + (4.9013731987690712d-3)*(0.0d0,1.0d0) 
     case(14)
        we = (3.9783815941610621d-4)*(1.0d0,0.0d0) + (4.5352749686171894d-7)*(0.0d0,1.0d0) 
        xe = (9.9999935022357100d-1)*(1.0d0,0.0d0) + (1.1399791383407472d-3)*(0.0d0,1.0d0) 
     end select


  case(15)
     select case(e)
        !------------------------   n = 15, S = 1.0d06, rate = 1.18d-4  ------------------------
     case(0)
        we = (-1.1792784382658184d-4)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)       
     case(1)
        we = (-3.6490219268281893d-4)*(1.0d0,0.0d0) + (3.8584221393928272d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999944096863402d-1)*(1.0d0,0.0d0) + (1.0573847074050477d-3)*(0.0d0,1.0d0) 
     case(2)
        we = (-7.7287550656539642d-4)*(1.0d0,0.0d0) + (3.3654081358322242d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999051974059949d-1)*(1.0d0,0.0d0) + (4.3543574641737573d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-2.0447659312819748d-3)*(1.0d0,0.0d0) + (2.5601146431954891d-5)*(0.0d0,1.0d0) 
        xe = (-9.9992162986865651d-1)*(1.0d0,0.0d0) + (1.2519349855705481d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-5.5985191613083579d-3)*(1.0d0,0.0d0) + (1.9423338474605119d-4)*(0.0d0,1.0d0) 
        xe = (-9.9939871608767983d-1)*(1.0d0,0.0d0) + (3.4672846469494513d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-1.5321885281304451d-2)*(1.0d0,0.0d0) + (1.4685091672664223d-3)*(0.0d0,1.0d0) 
        xe = (-9.9543837740544383d-1)*(1.0d0,0.0d0) + (9.5406691528514928d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (-4.0349347191103273d-2)*(1.0d0,0.0d0) + (1.0819051287289454d-2)*(0.0d0,1.0d0) 
        xe = (-9.6588107005927215d-1)*(1.0d0,0.0d0) + (2.5898601989519698d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (-7.9372797116653598d-2)*(1.0d0,0.0d0) + (6.6361266089171284d-2)*(0.0d0,1.0d0) 
        xe = (-7.6718747568358492d-1)*(1.0d0,0.0d0) + (6.4142293157810371d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.6129631883911560d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(9)
        we = (7.9372797116653598d-2)*(1.0d0,0.0d0) + (6.6361266089171284d-2)*(0.0d0,1.0d0) 
        xe = (7.6718747568358492d-1)*(1.0d0,0.0d0) + (6.4142293157810371d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (4.0349347191103273d-2)*(1.0d0,0.0d0) + (1.0819051287289454d-2)*(0.0d0,1.0d0) 
        xe = (9.6588107005927215d-1)*(1.0d0,0.0d0) + (2.5898601989519698d-1)*(0.0d0,1.0d0) 
     case(11)
        we = (1.5321885281304451d-2)*(1.0d0,0.0d0) + (1.4685091672664223d-3)*(0.0d0,1.0d0) 
        xe = (9.9543837740544383d-1)*(1.0d0,0.0d0) + (9.5406691528514928d-2)*(0.0d0,1.0d0) 
     case(12)
        we = (5.5985191613083579d-3)*(1.0d0,0.0d0) + (1.9423338474605119d-4)*(0.0d0,1.0d0) 
        xe = (9.9939871608767983d-1)*(1.0d0,0.0d0) + (3.4672846469494513d-2)*(0.0d0,1.0d0) 
     case(13)
        we = (2.0447659312819748d-3)*(1.0d0,0.0d0) + (2.5601146431954891d-5)*(0.0d0,1.0d0) 
        xe = (9.9992162986865651d-1)*(1.0d0,0.0d0) + (1.2519349855705481d-2)*(0.0d0,1.0d0) 
     case(14)
        we = (7.7287550656539642d-4)*(1.0d0,0.0d0) + (3.3654081358322242d-6)*(0.0d0,1.0d0) 
        xe = (9.9999051974059949d-1)*(1.0d0,0.0d0) + (4.3543574641737573d-3)*(0.0d0,1.0d0) 
     case(15)
        we = (3.6490219268281893d-4)*(1.0d0,0.0d0) + (3.8584221393928272d-7)*(0.0d0,1.0d0) 
        xe = (9.9999944096863402d-1)*(1.0d0,0.0d0) + (1.0573847074050477d-3)*(0.0d0,1.0d0) 
     end select


  case(16)
     select case(e)
        !------------------------   n = 16, S = 1.0d06, rate = 6.16d-5  ------------------------
     case(0)
        we = (6.1610602189565711d-5)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)       
     case(1)
        we = (-3.3720271377405092d-4)*(1.0d0,0.0d0) + (3.3256791542695477d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999951365073736d-1)*(1.0d0,0.0d0) + (9.8625467737379399d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-6.6519220239713778d-4)*(1.0d0,0.0d0) + (2.6062953723407513d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999232430038310d-1)*(1.0d0,0.0d0) + (3.9180786512652823d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.6400960309060023d-3)*(1.0d0,0.0d0) + (1.7485654672736178d-5)*(0.0d0,1.0d0) 
        xe = (-9.9994317254293952d-1)*(1.0d0,0.0d0) + (1.0660754418012596d-2)*(0.0d0,1.0d0) 
     case(4)
        we = (-4.2082918803256564d-3)*(1.0d0,0.0d0) + (1.1690372030085671d-4)*(0.0d0,1.0d0) 
        xe = (-9.9961437632196248d-1)*(1.0d0,0.0d0) + (2.7768663101670359d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-1.0834902847825161d-2)*(1.0d0,0.0d0) + (7.8006735739269233d-4)*(0.0d0,1.0d0) 
        xe = (-9.9741833534204793d-1)*(1.0d0,0.0d0) + (7.1809918002307196d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (-2.7404178905123672d-2)*(1.0d0,0.0d0) + (5.1408358392403394d-3)*(0.0d0,1.0d0) 
        xe = (-9.8285560240079251d-1)*(1.0d0,0.0d0) + (1.8437696393360870d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (-6.1233481833219840d-2)*(1.0d0,0.0d0) + (3.1256621441574164d-2)*(0.0d0,1.0d0) 
        xe = (-8.9067323792710495d-1)*(1.0d0,0.0d0) + (4.5464401815095595d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (-5.9982325962616087d-2)*(1.0d0,0.0d0) + (1.2163640524928188d-1)*(0.0d0,1.0d0) 
        xe = (-4.4227617016573273d-1)*(1.0d0,0.0d0) + (8.9687891563105204d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (5.9982325962616087d-2)*(1.0d0,0.0d0) + (1.2163640524928188d-1)*(0.0d0,1.0d0) 
        xe = (4.4227617016573273d-1)*(1.0d0,0.0d0) + (8.9687891563105204d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (6.1233481833219840d-2)*(1.0d0,0.0d0) + (3.1256621441574164d-2)*(0.0d0,1.0d0) 
        xe = (8.9067323792710495d-1)*(1.0d0,0.0d0) + (4.5464401815095595d-1)*(0.0d0,1.0d0) 
     case(11)
        we = (2.7404178905123672d-2)*(1.0d0,0.0d0) + (5.1408358392403394d-3)*(0.0d0,1.0d0) 
        xe = (9.8285560240079251d-1)*(1.0d0,0.0d0) + (1.8437696393360870d-1)*(0.0d0,1.0d0) 
     case(12)
        we = (1.0834902847825161d-2)*(1.0d0,0.0d0) + (7.8006735739269233d-4)*(0.0d0,1.0d0) 
        xe = (9.9741833534204793d-1)*(1.0d0,0.0d0) + (7.1809918002307196d-2)*(0.0d0,1.0d0) 
     case(13)
        we = (4.2082918803256564d-3)*(1.0d0,0.0d0) + (1.1690372030085671d-4)*(0.0d0,1.0d0) 
        xe = (9.9961437632196248d-1)*(1.0d0,0.0d0) + (2.7768663101670359d-2)*(0.0d0,1.0d0) 
     case(14)
        we = (1.6400960309060023d-3)*(1.0d0,0.0d0) + (1.7485654672736178d-5)*(0.0d0,1.0d0) 
        xe = (9.9994317254293952d-1)*(1.0d0,0.0d0) + (1.0660754418012596d-2)*(0.0d0,1.0d0) 
     case(15)
        we = (6.6519220239713778d-4)*(1.0d0,0.0d0) + (2.6062953723407513d-6)*(0.0d0,1.0d0) 
        xe = (9.9999232430038310d-1)*(1.0d0,0.0d0) + (3.9180786512652823d-3)*(0.0d0,1.0d0) 
     case(16)
        we = (3.3720271377405092d-4)*(1.0d0,0.0d0) + (3.3256791542695477d-7)*(0.0d0,1.0d0) 
        xe = (9.9999951365073736d-1)*(1.0d0,0.0d0) + (9.8625467737379399d-4)*(0.0d0,1.0d0) 
     end select


  case(17)
     select case(e)
        !------------------------   n = 17, S = 1.0d06, rate = 3.22d-5  ------------------------
     case(0)
        we = (-3.2188041018121893d-5)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)       
     case(1)
        we = (-3.1356830506300385d-4)*(1.0d0,0.0d0) + (2.8983656153429189d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999957281917840d-1)*(1.0d0,0.0d0) + (9.2431675345324153d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-5.8146393205053479d-4)*(1.0d0,0.0d0) + (2.0715601964910134d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999365377560656d-1)*(1.0d0,0.0d0) + (3.5626406656038721d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.3460694593148770d-3)*(1.0d0,0.0d0) + (1.2444383027989803d-5)*(0.0d0,1.0d0) 
        xe = (-9.9995726792672457d-1)*(1.0d0,0.0d0) + (9.2445833070427582d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-3.2596721764997800d-3)*(1.0d0,0.0d0) + (7.4412679397007801d-5)*(0.0d0,1.0d0) 
        xe = (-9.9973953682146444d-1)*(1.0d0,0.0d0) + (2.2822324947379177d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-7.9431247622769856d-3)*(1.0d0,0.0d0) + (4.4442256920342060d-4)*(0.0d0,1.0d0) 
        xe = (-9.9843843074930583d-1)*(1.0d0,0.0d0) + (5.5863225854435435d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (-1.9193908884236641d-2)*(1.0d0,0.0d0) + (2.6371573117971421d-3)*(0.0d0,1.0d0) 
        xe = (-9.9069280131322413d-1)*(1.0d0,0.0d0) + (1.3611676394242028d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (-4.3784116406463298d-2)*(1.0d0,0.0d0) + (1.5064213717605706d-2)*(0.0d0,1.0d0) 
        xe = (-9.4559743883263248d-1)*(1.0d0,0.0d0) + (3.2533902881942373d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (-7.1148642777824761d-2)*(1.0d0,0.0d0) + (6.9872895653441694d-2)*(0.0d0,1.0d0) 
        xe = (-7.1347440722420030d-1)*(1.0d0,0.0d0) + (7.0068129005709578d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.4232028316042250d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(10)
        we = (7.1148642777824761d-2)*(1.0d0,0.0d0) + (6.9872895653441694d-2)*(0.0d0,1.0d0) 
        xe = (7.1347440722420030d-1)*(1.0d0,0.0d0) + (7.0068129005709578d-1)*(0.0d0,1.0d0) 
     case(11)
        we = (4.3784116406463298d-2)*(1.0d0,0.0d0) + (1.5064213717605706d-2)*(0.0d0,1.0d0) 
        xe = (9.4559743883263248d-1)*(1.0d0,0.0d0) + (3.2533902881942373d-1)*(0.0d0,1.0d0) 
     case(12)
        we = (1.9193908884236641d-2)*(1.0d0,0.0d0) + (2.6371573117971421d-3)*(0.0d0,1.0d0) 
        xe = (9.9069280131322413d-1)*(1.0d0,0.0d0) + (1.3611676394242028d-1)*(0.0d0,1.0d0) 
     case(13)
        we = (7.9431247622769856d-3)*(1.0d0,0.0d0) + (4.4442256920342060d-4)*(0.0d0,1.0d0) 
        xe = (9.9843843074930583d-1)*(1.0d0,0.0d0) + (5.5863225854435435d-2)*(0.0d0,1.0d0) 
     case(14)
        we = (3.2596721764997800d-3)*(1.0d0,0.0d0) + (7.4412679397007801d-5)*(0.0d0,1.0d0) 
        xe = (9.9973953682146444d-1)*(1.0d0,0.0d0) + (2.2822324947379177d-2)*(0.0d0,1.0d0) 
     case(15)
        we = (1.3460694593148770d-3)*(1.0d0,0.0d0) + (1.2444383027989803d-5)*(0.0d0,1.0d0) 
        xe = (9.9995726792672457d-1)*(1.0d0,0.0d0) + (9.2445833070427582d-3)*(0.0d0,1.0d0) 
     case(16)
        we = (5.8146393205053479d-4)*(1.0d0,0.0d0) + (2.0715601964910134d-6)*(0.0d0,1.0d0) 
        xe = (9.9999365377560656d-1)*(1.0d0,0.0d0) + (3.5626406656038721d-3)*(0.0d0,1.0d0) 
     case(17)
        we = (3.1356830506300385d-4)*(1.0d0,0.0d0) + (2.8983656153429189d-7)*(0.0d0,1.0d0) 
        xe = (9.9999957281917840d-1)*(1.0d0,0.0d0) + (9.2431675345324153d-4)*(0.0d0,1.0d0) 
     end select




  case(18)
     select case(e)
        !------------------------   n = 18, S = 1.0d06, rate = 1.68d-5  ------------------------
     case(0)
        we = (1.6816423564713912d-5)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)       
     case(1)
        we = (-2.9315283498838234d-4)*(1.0d0,0.0d0) + (2.5500445397162234d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999962166450085d-1)*(1.0d0,0.0d0) + (8.6986829755370020d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-5.1496929311266781d-4)*(1.0d0,0.0d0) + (1.6828261686953906d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999466072397813d-1)*(1.0d0,0.0d0) + (3.2678010245045879d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.1264034015766962d-3)*(1.0d0,0.0d0) + (9.1672826634493694d-6)*(0.0d0,1.0d0) 
        xe = (-9.9996688370298181d-1)*(1.0d0,0.0d0) + (8.1382736097585021d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-2.5895909782401231d-3)*(1.0d0,0.0d0) + (4.9639523202098925d-5)*(0.0d0,1.0d0) 
        xe = (-9.9981632789817620d-1)*(1.0d0,0.0d0) + (1.9165345501885981d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-6.0051538188240234d-3)*(1.0d0,0.0d0) + (2.6855850880543866d-4)*(0.0d0,1.0d0) 
        xe = (-9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (-1.3875670060744669d-2)*(1.0d0,0.0d0) + (1.4477284130137483d-3)*(0.0d0,1.0d0) 
        xe = (-9.9460106547849614d-1)*(1.0d0,0.0d0) + (1.0377244600104668d-1)*(0.0d0,1.0d0) 
     case(7)
        we = (-3.1151203067143453d-2)*(1.0d0,0.0d0) + (7.6552207011467777d-3)*(0.0d0,1.0d0) 
        xe = (-9.7110718365323745d-1)*(1.0d0,0.0d0) + (2.3864374673784650d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (-5.9847706053290449d-2)*(1.0d0,0.0d0) + (3.6628501070755985d-2)*(0.0d0,1.0d0) 
        xe = (-8.5293349191434553d-1)*(1.0d0,0.0d0) + (5.2201959577280344d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (-4.9161393686735257d-2)*(1.0d0,0.0d0) + (1.1303165332996631d-1)*(0.0d0,1.0d0) 
        xe = (-3.9884344260344645d-1)*(1.0d0,0.0d0) + (9.1701903376769189d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (4.9161393686735257d-2)*(1.0d0,0.0d0) + (1.1303165332996631d-1)*(0.0d0,1.0d0) 
        xe = (3.9884344260344645d-1)*(1.0d0,0.0d0) + (9.1701903376769189d-1)*(0.0d0,1.0d0) 
     case(11)
        we = (5.9847706053290449d-2)*(1.0d0,0.0d0) + (3.6628501070755985d-2)*(0.0d0,1.0d0) 
        xe = (8.5293349191434553d-1)*(1.0d0,0.0d0) + (5.2201959577280344d-1)*(0.0d0,1.0d0) 
     case(12)
        we = (3.1151203067143453d-2)*(1.0d0,0.0d0) + (7.6552207011467777d-3)*(0.0d0,1.0d0) 
        xe = (9.7110718365323745d-1)*(1.0d0,0.0d0) + (2.3864374673784650d-1)*(0.0d0,1.0d0) 
     case(13)
        we = (1.3875670060744669d-2)*(1.0d0,0.0d0) + (1.4477284130137483d-3)*(0.0d0,1.0d0) 
        xe = (9.9460106547849614d-1)*(1.0d0,0.0d0) + (1.0377244600104668d-1)*(0.0d0,1.0d0) 
     case(14)
        we = (6.0051538188240234d-3)*(1.0d0,0.0d0) + (2.6855850880543866d-4)*(0.0d0,1.0d0) 
        xe = (9.9900149850137365d-1)*(1.0d0,0.0d0) + (4.4676682867128663d-2)*(0.0d0,1.0d0) 
     case(15)
        we = (2.5895909782401231d-3)*(1.0d0,0.0d0) + (4.9639523202098925d-5)*(0.0d0,1.0d0) 
        xe = (9.9981632789817620d-1)*(1.0d0,0.0d0) + (1.9165345501885981d-2)*(0.0d0,1.0d0) 
     case(16)
        we = (1.1264034015766962d-3)*(1.0d0,0.0d0) + (9.1672826634493694d-6)*(0.0d0,1.0d0) 
        xe = (9.9996688370298181d-1)*(1.0d0,0.0d0) + (8.1382736097585021d-3)*(0.0d0,1.0d0) 
     case(17)
        we = (5.1496929311266781d-4)*(1.0d0,0.0d0) + (1.6828261686953906d-6)*(0.0d0,1.0d0) 
        xe = (9.9999466072397813d-1)*(1.0d0,0.0d0) + (3.2678010245045879d-3)*(0.0d0,1.0d0) 
     case(18)
        we = (2.9315283498838234d-4)*(1.0d0,0.0d0) + (2.5500445397162234d-7)*(0.0d0,1.0d0) 
        xe = (9.9999962166450085d-1)*(1.0d0,0.0d0) + (8.6986829755370020d-4)*(0.0d0,1.0d0) 
     end select




  case(19)
     select case(e)
        !------------------------   n = 19, S = 1.0e+06, rate = 8.79d-6  ------------------------
     case(0)
        we = (-8.7856263509822341d-6)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)  
     case(1)
        we = (-2.7533033809645487d-4)*(1.0d0,0.0d0) + (2.2621346757974878d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999966248073835d-1)*(1.0d0,0.0d0) + (8.2160721114975059d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-4.6118603432212222d-4)*(1.0d0,0.0d0) + (1.3925261592991639d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999544150519182d-1)*(1.0d0,0.0d0) + (3.0194318731327206d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-9.5833058171406892d-4)*(1.0d0,0.0d0) + (6.9532666952665401d-6)*(0.0d0,1.0d0) 
        xe = (-9.9997367914835145d-1)*(1.0d0,0.0d0) + (7.2554124975633188d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-2.1020598481729306d-3)*(1.0d0,0.0d0) + (3.4453656580025595d-5)*(0.0d0,1.0d0) 
        xe = (-9.9986570402423880d-1)*(1.0d0,0.0d0) + (1.6388224922586798d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-4.6611742142131079d-3)*(1.0d0,0.0d0) + (1.7059346848525883d-4)*(0.0d0,1.0d0) 
        xe = (-9.9933093531757289d-1)*(1.0d0,0.0d0) + (3.6574331399287265d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (-1.0327323497992511d-2)*(1.0d0,0.0d0) + (8.4288136091312184d-4)*(0.0d0,1.0d0) 
        xe = (-9.9668591057932443d-1)*(1.0d0,0.0d0) + (8.1346147128569815d-2)*(0.0d0,1.0d0) 
     case(7)
        we = (-2.2535809698160335d-2)*(1.0d0,0.0d0) + (4.1214125318933632d-3)*(0.0d0,1.0d0) 
        xe = (-9.8368506773852704d-1)*(1.0d0,0.0d0) + (1.7989910368940024d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (-4.5528622674970827d-2)*(1.0d0,0.0d0) + (1.9161285626193562d-2)*(0.0d0,1.0d0) 
        xe = (-9.2169806995615144d-1)*(1.0d0,0.0d0) + (3.8790806622072871d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (-6.3225124338941188d-2)*(1.0d0,0.0d0) + (7.1180657114887469d-2)*(0.0d0,1.0d0) 
        xe = (-6.6409053886766345d-1)*(1.0d0,0.0d0) + (7.4765216256388667d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.2733920084458664d-1)*(0.0d0,1.0d0) 
        xe = (0.0000000000000000d0)*(1.0d0,0.0d0) + (1.0000000000000000d0)*(0.0d0,1.0d0) 
     case(11)
        we = (6.3225124338941188d-2)*(1.0d0,0.0d0) + (7.1180657114887469d-2)*(0.0d0,1.0d0) 
        xe = (6.6409053886766345d-1)*(1.0d0,0.0d0) + (7.4765216256388667d-1)*(0.0d0,1.0d0) 
     case(12)
        we = (4.5528622674970827d-2)*(1.0d0,0.0d0) + (1.9161285626193562d-2)*(0.0d0,1.0d0) 
        xe = (9.2169806995615144d-1)*(1.0d0,0.0d0) + (3.8790806622072871d-1)*(0.0d0,1.0d0) 
     case(13)
        we = (2.2535809698160335d-2)*(1.0d0,0.0d0) + (4.1214125318933632d-3)*(0.0d0,1.0d0) 
        xe = (9.8368506773852704d-1)*(1.0d0,0.0d0) + (1.7989910368940024d-1)*(0.0d0,1.0d0) 
     case(14)
        we = (1.0327323497992511d-2)*(1.0d0,0.0d0) + (8.4288136091312184d-4)*(0.0d0,1.0d0) 
        xe = (9.9668591057932443d-1)*(1.0d0,0.0d0) + (8.1346147128569815d-2)*(0.0d0,1.0d0) 
     case(15)
        we = (4.6611742142131079d-3)*(1.0d0,0.0d0) + (1.7059346848525883d-4)*(0.0d0,1.0d0) 
        xe = (9.9933093531757289d-1)*(1.0d0,0.0d0) + (3.6574331399287265d-2)*(0.0d0,1.0d0) 
     case(16)
        we = (2.1020598481729306d-3)*(1.0d0,0.0d0) + (3.4453656580025595d-5)*(0.0d0,1.0d0) 
        xe = (9.9986570402423880d-1)*(1.0d0,0.0d0) + (1.6388224922586798d-2)*(0.0d0,1.0d0) 
     case(17)
        we = (9.5833058171406892d-4)*(1.0d0,0.0d0) + (6.9532666952665401d-6)*(0.0d0,1.0d0) 
        xe = (9.9997367914835145d-1)*(1.0d0,0.0d0) + (7.2554124975633188d-3)*(0.0d0,1.0d0) 
     case(18)
        we = (4.6118603432212222d-4)*(1.0d0,0.0d0) + (1.3925261592991639d-6)*(0.0d0,1.0d0) 
        xe = (9.9999544150519182d-1)*(1.0d0,0.0d0) + (3.0194318731327206d-3)*(0.0d0,1.0d0) 
     case(19)
        we = (2.7533033809645487d-4)*(1.0d0,0.0d0) + (2.2621346757974878d-7)*(0.0d0,1.0d0) 
        xe = (9.9999966248073835d-1)*(1.0d0,0.0d0) + (8.2160721114975059d-4)*(0.0d0,1.0d0) 

     end select



  case(20)
     select case(e)
        !------------------------   n = 20, S = 1.0e+06, rate = 4.59d-6  ------------------------
     case(0)
        we = (4.5899908547308854d-6)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)  
     case(1)
        we = (-2.5962796229929012d-4)*(1.0d0,0.0d0) + (2.0212559742093714d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999969695333402d-1)*(1.0d0,0.0d0) + (7.7851990354738988d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-4.1698448757787214d-4)*(1.0d0,0.0d0) + (1.1706511996116833d-6)*(0.0d0,1.0d0) 
        xe = (-9.9999605921566304d-1)*(1.0d0,0.0d0) + (2.8074103982325071d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-8.2705236048838100d-4)*(1.0d0,0.0d0) + (5.4072170850776506d-6)*(0.0d0,1.0d0) 
        xe = (-9.9997862836826079d-1)*(1.0d0,0.0d0) + (6.5377983092055223d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-1.7381938336616146d-3)*(1.0d0,0.0d0) + (2.4737279114863458d-5)*(0.0d0,1.0d0) 
        xe = (-9.9989874616333840d-1)*(1.0d0,0.0d0) + (1.4230158853077478d-2)*(0.0d0,1.0d0) 
     case(5)
        we = (-3.7009135179021358d-3)*(1.0d0,0.0d0) + (1.1308608306421491d-4)*(0.0d0,1.0d0) 
        xe = (-9.9953348408535481d-1)*(1.0d0,0.0d0) + (3.0542007009883753d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (-7.8878951883124531d-3)*(1.0d0,0.0d0) + (5.1628614143891432d-4)*(0.0d0,1.0d0) 
        xe = (-9.9786481263071369d-1)*(1.0d0,0.0d0) + (6.5313212395890802d-2)*(0.0d0,1.0d0) 
     case(7)
        we = (-1.6672937540863027d-2)*(1.0d0,0.0d0) + (2.3430845492115548d-3)*(0.0d0,1.0d0) 
        xe = (-9.9026924750858347d-1)*(1.0d0,0.0d0) + (1.3916471333920749d-1)*(0.0d0,1.0d0) 
     case(8)
        we = (-3.3840259241556050d-2)*(1.0d0,0.0d0) + (1.0351956600694706d-2)*(0.0d0,1.0d0) 
        xe = (-9.5625772508631735d-1)*(1.0d0,0.0d0) + (2.9252549156054936d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (-5.7167528994478772d-2)*(1.0d0,0.0d0) + (4.0725963384848053d-2)*(0.0d0,1.0d0) 
        xe = (-8.1446048593714426d-1)*(1.0d0,0.0d0) + (5.8021902489235122d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (-4.0893404149440093d-2)*(1.0d0,0.0d0) + (1.0505409220020073d-1)*(0.0d0,1.0d0) 
        xe = (-3.6274701643022350d-1)*(1.0d0,0.0d0) + (9.3188765528413942d-1)*(0.0d0,1.0d0) 
     case(11)
        we = (4.0893404149440093d-2)*(1.0d0,0.0d0) + (1.0505409220020073d-1)*(0.0d0,1.0d0) 
        xe = (3.6274701643022350d-1)*(1.0d0,0.0d0) + (9.3188765528413942d-1)*(0.0d0,1.0d0) 
     case(12)
        we = (5.7167528994478772d-2)*(1.0d0,0.0d0) + (4.0725963384848053d-2)*(0.0d0,1.0d0) 
        xe = (8.1446048593714426d-1)*(1.0d0,0.0d0) + (5.8021902489235122d-1)*(0.0d0,1.0d0) 
     case(13)
        we = (3.3840259241556050d-2)*(1.0d0,0.0d0) + (1.0351956600694706d-2)*(0.0d0,1.0d0) 
        xe = (9.5625772508631735d-1)*(1.0d0,0.0d0) + (2.9252549156054936d-1)*(0.0d0,1.0d0) 
     case(14)
        we = (1.6672937540863027d-2)*(1.0d0,0.0d0) + (2.3430845492115548d-3)*(0.0d0,1.0d0) 
        xe = (9.9026924750858347d-1)*(1.0d0,0.0d0) + (1.3916471333920749d-1)*(0.0d0,1.0d0) 
     case(15)
        we = (7.8878951883124531d-3)*(1.0d0,0.0d0) + (5.1628614143891432d-4)*(0.0d0,1.0d0) 
        xe = (9.9786481263071369d-1)*(1.0d0,0.0d0) + (6.5313212395890802d-2)*(0.0d0,1.0d0) 
     case(16)
        we = (3.7009135179021358d-3)*(1.0d0,0.0d0) + (1.1308608306421491d-4)*(0.0d0,1.0d0) 
        xe = (9.9953348408535481d-1)*(1.0d0,0.0d0) + (3.0542007009883753d-2)*(0.0d0,1.0d0) 
     case(17)
        we = (1.7381938336616146d-3)*(1.0d0,0.0d0) + (2.4737279114863458d-5)*(0.0d0,1.0d0) 
        xe = (9.9989874616333840d-1)*(1.0d0,0.0d0) + (1.4230158853077478d-2)*(0.0d0,1.0d0) 
     case(18)
        we = (8.2705236048838100d-4)*(1.0d0,0.0d0) + (5.4072170850776506d-6)*(0.0d0,1.0d0) 
        xe = (9.9997862836826079d-1)*(1.0d0,0.0d0) + (6.5377983092055223d-3)*(0.0d0,1.0d0) 
     case(19)
        we = (4.1698448757787214d-4)*(1.0d0,0.0d0) + (1.1706511996116833d-6)*(0.0d0,1.0d0) 
        xe = (9.9999605921566304d-1)*(1.0d0,0.0d0) + (2.8074103982325071d-3)*(0.0d0,1.0d0) 
     case(20)
        we = (2.5962796229929012d-4)*(1.0d0,0.0d0) + (2.0212559742093714d-7)*(0.0d0,1.0d0) 
        xe = (9.9999969695333402d-1)*(1.0d0,0.0d0) + (7.7851990354738988d-4)*(0.0d0,1.0d0) 

     end select



  case(24)
     select case(e)
        !------------------------   n = 24, S = 1.0e+06, rate = 3.42d-7  ------------------------
     case(0)
        we = (3.4195439124751204d-7)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)  
     case(1)
        we = (-2.1181634653955156d-4)*(1.0d0,0.0d0) + (1.3642042022693250d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999979259955241d-1)*(1.0d0,0.0d0) + (6.4405034909712736d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-2.9967697189998575d-4)*(1.0d0,0.0d0) + (6.5908154932765347d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-5.1183808502514039d-4)*(1.0d0,0.0d0) + (2.3886747768399971d-6)*(0.0d0,1.0d0) 
        xe = (-9.9998911040415039d-1)*(1.0d0,0.0d0) + (4.6668054508287246d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-9.3627144787552391d-4)*(1.0d0,0.0d0) + (8.4923098454392950d-6)*(0.0d0,1.0d0) 
        xe = (-9.9995886691232916d-1)*(1.0d0,0.0d0) + (9.0699770347401325d-3)*(0.0d0,1.0d0) 
     case(5)
        we = (-1.7488125466066549d-3)*(1.0d0,0.0d0) + (3.0144437200358861d-5)*(0.0d0,1.0d0) 
        xe = (-9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(6)
        we = (-3.2850345253919943d-3)*(1.0d0,0.0d0) + (1.0695910797948011d-4)*(0.0d0,1.0d0) 
        xe = (-9.9947036029338032d-1)*(1.0d0,0.0d0) + (3.2542263212942395d-2)*(0.0d0,1.0d0) 
     case(7)
        we = (-6.1723524218736843d-3)*(1.0d0,0.0d0) + (3.7914247632396840d-4)*(0.0d0,1.0d0) 
        xe = (-9.9811874967746239d-1)*(1.0d0,0.0d0) + (6.1310370593392005d-2)*(0.0d0,1.0d0) 
     case(8)
        we = (-1.1543217321753032d-2)*(1.0d0,0.0d0) + (1.3393451640454283d-3)*(0.0d0,1.0d0) 
        xe = (-9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(9)
        we = (-2.1198732713608778d-2)*(1.0d0,0.0d0) + (4.6742873692502919d-3)*(0.0d0,1.0d0) 
        xe = (-9.7654221667714647d-1)*(1.0d0,0.0d0) + (2.1532602965569464d-1)*(0.0d0,1.0d0) 
     case(10)
        we = (-3.6499594635495881d-2)*(1.0d0,0.0d0) + (1.5642032244570850d-2)*(0.0d0,1.0d0) 
        xe = (-9.1915096715797229d-1)*(1.0d0,0.0d0) + (3.9390544496435209d-1)*(0.0d0,1.0d0) 
     case(11)
        we = (-5.0178640072982915d-2)*(1.0d0,0.0d0) + (4.5631281558717415d-2)*(0.0d0,1.0d0) 
        xe = (-7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(12)
        we = (-2.9413295974037672d-2)*(1.0d0,0.0d0) + (9.1338388378113380d-2)*(0.0d0,1.0d0) 
        xe = (-3.0652417782216373d-1)*(1.0d0,0.0d0) + (9.5186287269251479d-1)*(0.0d0,1.0d0) 
     case(13)
        we = (2.9413295974037672d-2)*(1.0d0,0.0d0) + (9.1338388378113380d-2)*(0.0d0,1.0d0) 
        xe = (3.0652417782216373d-1)*(1.0d0,0.0d0) + (9.5186287269251479d-1)*(0.0d0,1.0d0) 
     case(14)
        we = (5.0178640072982915d-2)*(1.0d0,0.0d0) + (4.5631281558717415d-2)*(0.0d0,1.0d0) 
        xe = (7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(15)
        we = (3.6499594635495881d-2)*(1.0d0,0.0d0) + (1.5642032244570850d-2)*(0.0d0,1.0d0) 
        xe = (9.1915096715797229d-1)*(1.0d0,0.0d0) + (3.9390544496435209d-1)*(0.0d0,1.0d0) 
     case(16)
        we = (2.1198732713608778d-2)*(1.0d0,0.0d0) + (4.6742873692502919d-3)*(0.0d0,1.0d0) 
        xe = (9.7654221667714647d-1)*(1.0d0,0.0d0) + (2.1532602965569464d-1)*(0.0d0,1.0d0) 
     case(17)
        we = (1.1543217321753032d-2)*(1.0d0,0.0d0) + (1.3393451640454283d-3)*(0.0d0,1.0d0) 
        xe = (9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(18)
        we = (6.1723524218736843d-3)*(1.0d0,0.0d0) + (3.7914247632396840d-4)*(0.0d0,1.0d0) 
        xe = (9.9811874967746239d-1)*(1.0d0,0.0d0) + (6.1310370593392005d-2)*(0.0d0,1.0d0) 
     case(19)
        we = (3.2850345253919943d-3)*(1.0d0,0.0d0) + (1.0695910797948011d-4)*(0.0d0,1.0d0) 
        xe = (9.9947036029338032d-1)*(1.0d0,0.0d0) + (3.2542263212942395d-2)*(0.0d0,1.0d0) 
     case(20)
        we = (1.7488125466066549d-3)*(1.0d0,0.0d0) + (3.0144437200358861d-5)*(0.0d0,1.0d0) 
        xe = (9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(21)
        we = (9.3627144787552391d-4)*(1.0d0,0.0d0) + (8.4923098454392950d-6)*(0.0d0,1.0d0) 
        xe = (9.9995886691232916d-1)*(1.0d0,0.0d0) + (9.0699770347401325d-3)*(0.0d0,1.0d0) 
     case(22)
        we = (5.1183808502514039d-4)*(1.0d0,0.0d0) + (2.3886747768399971d-6)*(0.0d0,1.0d0) 
        xe = (9.9998911040415039d-1)*(1.0d0,0.0d0) + (4.6668054508287246d-3)*(0.0d0,1.0d0) 
     case(23)
        we = (2.9967697189998575d-4)*(1.0d0,0.0d0) + (6.5908154932765347d-7)*(0.0d0,1.0d0) 
        xe = (9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     case(24)
        we = (2.1181634653955156d-4)*(1.0d0,0.0d0) + (1.3642042022693250d-7)*(0.0d0,1.0d0) 
        xe = (9.9999979259955241d-1)*(1.0d0,0.0d0) + (6.4405034909712736d-4)*(0.0d0,1.0d0) 

     end select




  case(32)
     select case(e)
        !------------------------   n = 32, S = 1.0e+06, rate = 1.90d-9  ------------------------

     case(0)
        we = (1.8979331439794578d-9)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)  
     case(1)
        we = (-1.5550096942091167d-4)*(1.0d0,0.0d0) + (7.4568333949205092d-8)*(0.0d0,1.0d0) 
        xe = (-9.9999988502256720d-1)*(1.0d0,0.0d0) + (4.7953608040659319d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.9125886486417976d-4)*(1.0d0,0.0d0) + (2.9623720666278892d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999880048610534d-1)*(1.0d0,0.0d0) + (1.5488790625862625d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-2.7099661183330458d-4)*(1.0d0,0.0d0) + (8.0605312562603178d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999557649285342d-1)*(1.0d0,0.0d0) + (2.9743898072795362d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-4.1304701073769696d-4)*(1.0d0,0.0d0) + (2.0998990288487514d-6)*(0.0d0,1.0d0) 
        xe = (-9.9998707711605239d-1)*(1.0d0,0.0d0) + (5.0838568916089094d-3)*(0.0d0,1.0d0) 
     case(5)
        we = (-6.5006154636645475d-4)*(1.0d0,0.0d0) + (5.4362011418772597d-6)*(0.0d0,1.0d0) 
        xe = (-9.9996503534011050d-1)*(1.0d0,0.0d0) + (8.3623021502164611d-3)*(0.0d0,1.0d0) 
     case(6)
        we = (-1.0364859003141988d-3)*(1.0d0,0.0d0) + (1.4059533908081803d-5)*(0.0d0,1.0d0) 
        xe = (-9.9990801328136669d-1)*(1.0d0,0.0d0) + (1.3563368892360707d-2)*(0.0d0,1.0d0) 
     case(7)
        we = (-1.6609433872528264d-3)*(1.0d0,0.0d0) + (3.6353506595416335d-5)*(0.0d0,1.0d0) 
        xe = (-9.9976055985963364d-1)*(1.0d0,0.0d0) + (2.1882023424535879d-2)*(0.0d0,1.0d0) 
     case(8)
        we = (-2.6660466676895595d-3)*(1.0d0,0.0d0) + (9.3974723940552110d-5)*(0.0d0,1.0d0) 
        xe = (-9.9937934228484449d-1)*(1.0d0,0.0d0) + (3.5226839402814472d-2)*(0.0d0,1.0d0) 
     case(9)
        we = (-4.2788366843057764d-3)*(1.0d0,0.0d0) + (2.4277947085922267d-4)*(0.0d0,1.0d0) 
        xe = (-9.9839418530232105d-1)*(1.0d0,0.0d0) + (5.6648484132539159d-2)*(0.0d0,1.0d0) 
     case(10)
        we = (-6.8532505272229727d-3)*(1.0d0,0.0d0) + (6.2623118524491391d-4)*(0.0d0,1.0d0) 
        xe = (-9.9585106320729500d-1)*(1.0d0,0.0d0) + (9.0998131348397709d-2)*(0.0d0,1.0d0) 
     case(11)
        we = (-1.0911580563159749d-2)*(1.0d0,0.0d0) + (1.6088283108259130d-3)*(0.0d0,1.0d0) 
        xe = (-9.8930445916784426d-1)*(1.0d0,0.0d0) + (1.4586530454710303d-1)*(0.0d0,1.0d0) 
     case(12)
        we = (-1.7104615708036029d-2)*(1.0d0,0.0d0) + (4.0907396646800102d-3)*(0.0d0,1.0d0) 
        xe = (-9.7257239010253027d-1)*(1.0d0,0.0d0) + (2.3260039984972461d-1)*(0.0d0,1.0d0) 
     case(13)
        we = (-2.5757868072832625d-2)*(1.0d0,0.0d0) + (1.0133058423378140d-2)*(0.0d0,1.0d0) 
        xe = (-9.3058053398093854d-1)*(1.0d0,0.0d0) + (3.6608724339390908d-1)*(0.0d0,1.0d0) 
     case(14)
        we = (-3.5010080480968125d-2)*(1.0d0,0.0d0) + (2.3540889306474175d-2)*(0.0d0,1.0d0) 
        xe = (-8.2984664904974115d-1)*(1.0d0,0.0d0) + (5.5799152239161820d-1)*(0.0d0,1.0d0) 
     case(15)
        we = (-3.6601603226060431d-2)*(1.0d0,0.0d0) + (4.7261950405979622d-2)*(0.0d0,1.0d0) 
        xe = (-6.1229562589776565d-1)*(1.0d0,0.0d0) + (7.9062890568614008d-1)*(0.0d0,1.0d0) 
     case(16)
        we = (-1.7142800732740431d-2)*(1.0d0,0.0d0) + (7.1497358882311499d-2)*(0.0d0,1.0d0) 
        xe = (-2.3315991108837755d-1)*(1.0d0,0.0d0) + (9.7243840723269459d-1)*(0.0d0,1.0d0) 
     case(17)
        we = (1.7142800732740431d-2)*(1.0d0,0.0d0) + (7.1497358882311499d-2)*(0.0d0,1.0d0) 
        xe = (2.3315991108837755d-1)*(1.0d0,0.0d0) + (9.7243840723269459d-1)*(0.0d0,1.0d0) 
     case(18)
        we = (3.6601603226060431d-2)*(1.0d0,0.0d0) + (4.7261950405979622d-2)*(0.0d0,1.0d0) 
        xe = (6.1229562589776565d-1)*(1.0d0,0.0d0) + (7.9062890568614008d-1)*(0.0d0,1.0d0) 
     case(19)
        we = (3.5010080480968125d-2)*(1.0d0,0.0d0) + (2.3540889306474175d-2)*(0.0d0,1.0d0) 
        xe = (8.2984664904974115d-1)*(1.0d0,0.0d0) + (5.5799152239161820d-1)*(0.0d0,1.0d0) 
     case(20)
        we = (2.5757868072832625d-2)*(1.0d0,0.0d0) + (1.0133058423378140d-2)*(0.0d0,1.0d0) 
        xe = (9.3058053398093854d-1)*(1.0d0,0.0d0) + (3.6608724339390908d-1)*(0.0d0,1.0d0) 
     case(21)
        we = (1.7104615708036029d-2)*(1.0d0,0.0d0) + (4.0907396646800102d-3)*(0.0d0,1.0d0) 
        xe = (9.7257239010253027d-1)*(1.0d0,0.0d0) + (2.3260039984972461d-1)*(0.0d0,1.0d0) 
     case(22)
        we = (1.0911580563159749d-2)*(1.0d0,0.0d0) + (1.6088283108259130d-3)*(0.0d0,1.0d0) 
        xe = (9.8930445916784426d-1)*(1.0d0,0.0d0) + (1.4586530454710303d-1)*(0.0d0,1.0d0) 
     case(23)
        we = (6.8532505272229727d-3)*(1.0d0,0.0d0) + (6.2623118524491391d-4)*(0.0d0,1.0d0) 
        xe = (9.9585106320729500d-1)*(1.0d0,0.0d0) + (9.0998131348397709d-2)*(0.0d0,1.0d0) 
     case(24)
        we = (4.2788366843057764d-3)*(1.0d0,0.0d0) + (2.4277947085922267d-4)*(0.0d0,1.0d0) 
        xe = (9.9839418530232105d-1)*(1.0d0,0.0d0) + (5.6648484132539159d-2)*(0.0d0,1.0d0) 
     case(25)
        we = (2.6660466676895595d-3)*(1.0d0,0.0d0) + (9.3974723940552110d-5)*(0.0d0,1.0d0) 
        xe = (9.9937934228484449d-1)*(1.0d0,0.0d0) + (3.5226839402814472d-2)*(0.0d0,1.0d0) 
     case(26)
        we = (1.6609433872528264d-3)*(1.0d0,0.0d0) + (3.6353506595416335d-5)*(0.0d0,1.0d0) 
        xe = (9.9976055985963364d-1)*(1.0d0,0.0d0) + (2.1882023424535879d-2)*(0.0d0,1.0d0) 
     case(27)
        we = (1.0364859003141988d-3)*(1.0d0,0.0d0) + (1.4059533908081803d-5)*(0.0d0,1.0d0) 
        xe = (9.9990801328136669d-1)*(1.0d0,0.0d0) + (1.3563368892360707d-2)*(0.0d0,1.0d0) 
     case(28)
        we = (6.5006154636645475d-4)*(1.0d0,0.0d0) + (5.4362011418772597d-6)*(0.0d0,1.0d0) 
        xe = (9.9996503534011050d-1)*(1.0d0,0.0d0) + (8.3623021502164611d-3)*(0.0d0,1.0d0) 
     case(29)
        we = (4.1304701073769696d-4)*(1.0d0,0.0d0) + (2.0998990288487514d-6)*(0.0d0,1.0d0) 
        xe = (9.9998707711605239d-1)*(1.0d0,0.0d0) + (5.0838568916089094d-3)*(0.0d0,1.0d0) 
     case(30)
        we = (2.7099661183330458d-4)*(1.0d0,0.0d0) + (8.0605312562603178d-7)*(0.0d0,1.0d0) 
        xe = (9.9999557649285342d-1)*(1.0d0,0.0d0) + (2.9743898072795362d-3)*(0.0d0,1.0d0) 
     case(31)
        we = (1.9125886486417976d-4)*(1.0d0,0.0d0) + (2.9623720666278892d-7)*(0.0d0,1.0d0) 
        xe = (9.9999880048610534d-1)*(1.0d0,0.0d0) + (1.5488790625862625d-3)*(0.0d0,1.0d0) 
     case(32)
        we = (1.5550096942091167d-4)*(1.0d0,0.0d0) + (7.4568333949205092d-8)*(0.0d0,1.0d0) 
        xe = (9.9999988502256720d-1)*(1.0d0,0.0d0) + (4.7953608040659319d-4)*(0.0d0,1.0d0) 

     end select



  case(40)
     select case(e)
        !------------------------   n = 40, S = 1.0e+06, rate = 1.05d-11  ------------------------
     case(0)
        we = (1.0534184635702104d-11)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)  
     case(1)
        we = (-1.2316275157188810d-4)*(1.0d0,0.0d0) + (4.7089605228083172d-8)*(0.0d0,1.0d0) 
        xe = (-9.9999992690943984d-1)*(1.0d0,0.0d0) + (3.8233638973867931d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.4116668350821769d-4)*(1.0d0,0.0d0) + (1.6980940530390962d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999927651656551d-1)*(1.0d0,0.0d0) + (1.2028991418483869d-3)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.7980618314001243d-4)*(1.0d0,0.0d0) + (3.9544892959663819d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-2.4472892912057994d-4)*(1.0d0,0.0d0) + (8.6076449501894848d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999381465333870d-1)*(1.0d0,0.0d0) + (3.5171942033385439d-3)*(0.0d0,1.0d0) 
     case(5)
        we = (-3.4542323313400009d-4)*(1.0d0,0.0d0) + (1.8477702426206758d-6)*(0.0d0,1.0d0) 
        xe = (-9.9998569283677707d-1)*(1.0d0,0.0d0) + (5.3492169287681366d-3)*(0.0d0,1.0d0) 
     case(6)
        we = (-4.9660191456992288d-4)*(1.0d0,0.0d0) + (3.9546355655326282d-6)*(0.0d0,1.0d0) 
        xe = (-9.9996829370430018d-1)*(1.0d0,0.0d0) + (7.9631392120442512d-3)*(0.0d0,1.0d0) 
     case(7)
        we = (-7.2034321732435960d-4)*(1.0d0,0.0d0) + (8.4580968155409520d-6)*(0.0d0,1.0d0) 
        xe = (-9.9993107266903825d-1)*(1.0d0,0.0d0) + (1.1740950172220682d-2)*(0.0d0,1.0d0) 
     case(8)
        we = (-1.0492875279641152d-3)*(1.0d0,0.0d0) + (1.8086662320217425d-5)*(0.0d0,1.0d0) 
        xe = (-9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(9)
        we = (-1.5313128096202157d-3)*(1.0d0,0.0d0) + (3.8671524949083791d-5)*(0.0d0,1.0d0) 
        xe = (-9.9968127428812481d-1)*(1.0d0,0.0d0) + (2.5245788513552012d-2)*(0.0d0,1.0d0) 
     case(10)
        we = (-2.2362336254057352d-3)*(1.0d0,0.0d0) + (8.2667929440118144d-5)*(0.0d0,1.0d0) 
        xe = (-9.9931740213171416d-1)*(1.0d0,0.0d0) + (3.6942249481074656d-2)*(0.0d0,1.0d0) 
     case(11)
        we = (-3.2650710840920827d-3)*(1.0d0,0.0d0) + (1.7664539309457733d-4)*(0.0d0,1.0d0) 
        xe = (-9.9853971624832527d-1)*(1.0d0,0.0d0) + (5.4022542282959551d-2)*(0.0d0,1.0d0) 
     case(12)
        we = (-4.7619077018430604d-3)*(1.0d0,0.0d0) + (3.7712311473102285d-4)*(0.0d0,1.0d0) 
        xe = (-9.9687868716768624d-1)*(1.0d0,0.0d0) + (7.8948610316017617d-2)*(0.0d0,1.0d0) 
     case(13)
        we = (-6.9259303930526263d-3)*(1.0d0,0.0d0) + (8.0360709842735071d-4)*(0.0d0,1.0d0) 
        xe = (-9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(14)
        we = (-1.0012930025861545d-2)*(1.0d0,0.0d0) + (1.7055229265479789d-3)*(0.0d0,1.0d0) 
        xe = (-9.8580171086824264d-1)*(1.0d0,0.0d0) + (1.6791362913487898d-1)*(0.0d0,1.0d0) 
     case(15)
        we = (-1.4290148409187728d-2)*(1.0d0,0.0d0) + (3.5889567295828485d-3)*(0.0d0,1.0d0) 
        xe = (-9.6987971939271034d-1)*(1.0d0,0.0d0) + (2.4358433839374305d-1)*(0.0d0,1.0d0) 
     case(16)
        we = (-1.9841360316044478d-2)*(1.0d0,0.0d0) + (7.4183032266655303d-3)*(0.0d0,1.0d0) 
        xe = (-9.3667338807504108d-1)*(1.0d0,0.0d0) + (3.5020417483522881d-1)*(0.0d0,1.0d0) 
     case(17)
        we = (-2.5990647961833875d-2)*(1.0d0,0.0d0) + (1.4779087936364021d-2)*(0.0d0,1.0d0) 
        xe = (-8.6928879621566224d-1)*(1.0d0,0.0d0) + (4.9430455063040313d-1)*(0.0d0,1.0d0) 
     case(18)
        we = (-3.0107184043793259d-2)*(1.0d0,0.0d0) + (2.7378768935233638d-2)*(0.0d0,1.0d0) 
        xe = (-7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(19)
        we = (-2.6715352778149914d-2)*(1.0d0,0.0d0) + (4.4418199008234377d-2)*(0.0d0,1.0d0) 
        xe = (-5.1540949947995840d-1)*(1.0d0,0.0d0) + (8.5694401675128040d-1)*(0.0d0,1.0d0) 
     case(20)
        we = (-1.1155359122546567d-2)*(1.0d0,0.0d0) + (5.8353572671816387d-2)*(0.0d0,1.0d0) 
        xe = (-1.8776815977090994d-1)*(1.0d0,0.0d0) + (9.8221337711122936d-1)*(0.0d0,1.0d0) 
     case(21)
        we = (1.1155359122546567d-2)*(1.0d0,0.0d0) + (5.8353572671816387d-2)*(0.0d0,1.0d0) 
        xe = (1.8776815977090994d-1)*(1.0d0,0.0d0) + (9.8221337711122936d-1)*(0.0d0,1.0d0) 
     case(22)
        we = (2.6715352778149914d-2)*(1.0d0,0.0d0) + (4.4418199008234377d-2)*(0.0d0,1.0d0) 
        xe = (5.1540949947995840d-1)*(1.0d0,0.0d0) + (8.5694401675128040d-1)*(0.0d0,1.0d0) 
     case(23)
        we = (3.0107184043793259d-2)*(1.0d0,0.0d0) + (2.7378768935233638d-2)*(0.0d0,1.0d0) 
        xe = (7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(24)
        we = (2.5990647961833875d-2)*(1.0d0,0.0d0) + (1.4779087936364021d-2)*(0.0d0,1.0d0) 
        xe = (8.6928879621566224d-1)*(1.0d0,0.0d0) + (4.9430455063040313d-1)*(0.0d0,1.0d0) 
     case(25)
        we = (1.9841360316044478d-2)*(1.0d0,0.0d0) + (7.4183032266655303d-3)*(0.0d0,1.0d0) 
        xe = (9.3667338807504108d-1)*(1.0d0,0.0d0) + (3.5020417483522881d-1)*(0.0d0,1.0d0) 
     case(26)
        we = (1.4290148409187728d-2)*(1.0d0,0.0d0) + (3.5889567295828485d-3)*(0.0d0,1.0d0) 
        xe = (9.6987971939271034d-1)*(1.0d0,0.0d0) + (2.4358433839374305d-1)*(0.0d0,1.0d0) 
     case(27)
        we = (1.0012930025861545d-2)*(1.0d0,0.0d0) + (1.7055229265479789d-3)*(0.0d0,1.0d0) 
        xe = (9.8580171086824264d-1)*(1.0d0,0.0d0) + (1.6791362913487898d-1)*(0.0d0,1.0d0) 
     case(28)
        we = (6.9259303930526263d-3)*(1.0d0,0.0d0) + (8.0360709842735071d-4)*(0.0d0,1.0d0) 
        xe = (9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(29)
        we = (4.7619077018430604d-3)*(1.0d0,0.0d0) + (3.7712311473102285d-4)*(0.0d0,1.0d0) 
        xe = (9.9687868716768624d-1)*(1.0d0,0.0d0) + (7.8948610316017617d-2)*(0.0d0,1.0d0) 
     case(30)
        we = (3.2650710840920827d-3)*(1.0d0,0.0d0) + (1.7664539309457733d-4)*(0.0d0,1.0d0) 
        xe = (9.9853971624832527d-1)*(1.0d0,0.0d0) + (5.4022542282959551d-2)*(0.0d0,1.0d0) 
     case(31)
        we = (2.2362336254057352d-3)*(1.0d0,0.0d0) + (8.2667929440118144d-5)*(0.0d0,1.0d0) 
        xe = (9.9931740213171416d-1)*(1.0d0,0.0d0) + (3.6942249481074656d-2)*(0.0d0,1.0d0) 
     case(32)
        we = (1.5313128096202157d-3)*(1.0d0,0.0d0) + (3.8671524949083791d-5)*(0.0d0,1.0d0) 
        xe = (9.9968127428812481d-1)*(1.0d0,0.0d0) + (2.5245788513552012d-2)*(0.0d0,1.0d0) 
     case(33)
        we = (1.0492875279641152d-3)*(1.0d0,0.0d0) + (1.8086662320217425d-5)*(0.0d0,1.0d0) 
        xe = (9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(34)
        we = (7.2034321732435960d-4)*(1.0d0,0.0d0) + (8.4580968155409520d-6)*(0.0d0,1.0d0) 
        xe = (9.9993107266903825d-1)*(1.0d0,0.0d0) + (1.1740950172220682d-2)*(0.0d0,1.0d0) 
     case(35)
        we = (4.9660191456992288d-4)*(1.0d0,0.0d0) + (3.9546355655326282d-6)*(0.0d0,1.0d0) 
        xe = (9.9996829370430018d-1)*(1.0d0,0.0d0) + (7.9631392120442512d-3)*(0.0d0,1.0d0) 
     case(36)
        we = (3.4542323313400009d-4)*(1.0d0,0.0d0) + (1.8477702426206758d-6)*(0.0d0,1.0d0) 
        xe = (9.9998569283677707d-1)*(1.0d0,0.0d0) + (5.3492169287681366d-3)*(0.0d0,1.0d0) 
     case(37)
        we = (2.4472892912057994d-4)*(1.0d0,0.0d0) + (8.6076449501894848d-7)*(0.0d0,1.0d0) 
        xe = (9.9999381465333870d-1)*(1.0d0,0.0d0) + (3.5171942033385439d-3)*(0.0d0,1.0d0) 
     case(38)
        we = (1.7980618314001243d-4)*(1.0d0,0.0d0) + (3.9544892959663819d-7)*(0.0d0,1.0d0) 
        xe = (9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     case(39)
        we = (1.4116668350821769d-4)*(1.0d0,0.0d0) + (1.6980940530390962d-7)*(0.0d0,1.0d0) 
        xe = (9.9999927651656551d-1)*(1.0d0,0.0d0) + (1.2028991418483869d-3)*(0.0d0,1.0d0) 
     case(40)
        we = (1.2316275157188810d-4)*(1.0d0,0.0d0) + (4.7089605228083172d-8)*(0.0d0,1.0d0) 
        xe = (9.9999992690943984d-1)*(1.0d0,0.0d0) + (3.8233638973867931d-4)*(0.0d0,1.0d0) 

     end select



  case(48)
     select case(e)
        !------------------------   n = 48, S = 1.0e+06, rate = 5.86d-14  ------------------------
     case(0)
        we = (5.8564264548977008d-14)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)  
     case(1)
        we = (-1.0207666216851158d-4)*(1.0d0,0.0d0) + (3.2463391730449853d-8)*(0.0d0,1.0d0) 
        xe = (-9.9999994942861647d-1)*(1.0d0,0.0d0) + (3.1802950252369703d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-1.1240090501800907d-4)*(1.0d0,0.0d0) + (1.1085597222977934d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999951365073736d-1)*(1.0d0,0.0d0) + (9.8625467737379399d-4)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.3409354176612142d-4)*(1.0d0,0.0d0) + (2.3523148062411967d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999846133468917d-1)*(1.0d0,0.0d0) + (1.7542315281065331d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-1.6934841332539216d-4)*(1.0d0,0.0d0) + (4.5718038200357800d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999635598178738d-1)*(1.0d0,0.0d0) + (2.6996338911749935d-3)*(0.0d0,1.0d0) 
     case(5)
        we = (-2.2173073497403943d-4)*(1.0d0,0.0d0) + (8.6876512741129979d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999232430038310d-1)*(1.0d0,0.0d0) + (3.9180786512652823d-3)*(0.0d0,1.0d0) 
     case(6)
        we = (-2.9653714734123929d-4)*(1.0d0,0.0d0) + (1.6407037867488601d-6)*(0.0d0,1.0d0) 
        xe = (-9.9998469398327960d-1)*(1.0d0,0.0d0) + (5.5327930710187233d-3)*(0.0d0,1.0d0) 
     case(7)
        we = (-4.0133002060141742d-4)*(1.0d0,0.0d0) + (3.0931700215852137d-6)*(0.0d0,1.0d0) 
        xe = (-9.9997030010278387d-1)*(1.0d0,0.0d0) + (7.7070689855653837d-3)*(0.0d0,1.0d0) 
     case(8)
        we = (-5.4669867904386235d-4)*(1.0d0,0.0d0) + (5.8285515797031256d-6)*(0.0d0,1.0d0) 
        xe = (-9.9994317254293952d-1)*(1.0d0,0.0d0) + (1.0660754418012596d-2)*(0.0d0,1.0d0) 
     case(9)
        we = (-7.4732098511848374d-4)*(1.0d0,0.0d0) + (1.0981144405960075d-5)*(0.0d0,1.0d0) 
        xe = (-9.9989206046724421d-1)*(1.0d0,0.0d0) + (1.4692427116338712d-2)*(0.0d0,1.0d0) 
     case(10)
        we = (-1.0234234700723219d-3)*(1.0d0,0.0d0) + (2.0687031904371676d-5)*(0.0d0,1.0d0) 
        xe = (-9.9979576857513042d-1)*(1.0d0,0.0d0) + (2.0209431938186778d-2)*(0.0d0,1.0d0) 
     case(11)
        we = (-1.4027639654332564d-3)*(1.0d0,0.0d0) + (3.8967906914869193d-5)*(0.0d0,1.0d0) 
        xe = (-9.9961437632196248d-1)*(1.0d0,0.0d0) + (2.7768663101670359d-2)*(0.0d0,1.0d0) 
     case(12)
        we = (-1.9232760208135970d-3)*(1.0d0,0.0d0) + (7.3391372929623125d-5)*(0.0d0,1.0d0) 
        xe = (-9.9927271804559314d-1)*(1.0d0,0.0d0) + (3.8131810470687581d-2)*(0.0d0,1.0d0) 
     case(13)
        we = (-2.6364847596615638d-3)*(1.0d0,0.0d0) + (1.3818198636717554d-4)*(0.0d0,1.0d0) 
        xe = (-9.9862934321886543d-1)*(1.0d0,0.0d0) + (5.2339610834027481d-2)*(0.0d0,1.0d0) 
     case(14)
        we = (-3.6116342963176664d-3)*(1.0d0,0.0d0) + (2.6002245345124118d-4)*(0.0d0,1.0d0) 
        xe = (-9.9741833534204793d-1)*(1.0d0,0.0d0) + (7.1809918002307196d-2)*(0.0d0,1.0d0) 
     case(15)
        we = (-4.9399091552901342d-3)*(1.0d0,0.0d0) + (4.8877131695888886d-4)*(0.0d0,1.0d0) 
        xe = (-9.9514075291324122d-1)*(1.0d0,0.0d0) + (9.8462591329232213d-2)*(0.0d0,1.0d0) 
     case(16)
        we = (-6.7365640264612303d-3)*(1.0d0,0.0d0) + (9.1691183433771899d-4)*(0.0d0,1.0d0) 
        xe = (-9.9086381972995408d-1)*(1.0d0,0.0d0) + (1.3486619572808120d-1)*(0.0d0,1.0d0) 
     case(17)
        we = (-9.1347263363820881d-3)*(1.0d0,0.0d0) + (1.7136119529180878d-3)*(0.0d0,1.0d0) 
        xe = (-9.8285560240079251d-1)*(1.0d0,0.0d0) + (1.8437696393360870d-1)*(0.0d0,1.0d0) 
     case(18)
        we = (-1.2254995789377469d-2)*(1.0d0,0.0d0) + (3.1800985784853704d-3)*(0.0d0,1.0d0) 
        xe = (-9.6794166373051727d-1)*(1.0d0,0.0d0) + (2.5117510946468846d-1)*(0.0d0,1.0d0) 
     case(19)
        we = (-1.6114837644635151d-2)*(1.0d0,0.0d0) + (5.8251138870511737d-3)*(0.0d0,1.0d0) 
        xe = (-9.4044447327612435d-1)*(1.0d0,0.0d0) + (3.3994733810458494d-1)*(0.0d0,1.0d0) 
     case(20)
        we = (-2.0411160688551314d-2)*(1.0d0,0.0d0) + (1.0418873853406705d-2)*(0.0d0,1.0d0) 
        xe = (-8.9067323792710495d-1)*(1.0d0,0.0d0) + (4.5464401815095595d-1)*(0.0d0,1.0d0) 
     case(21)
        we = (-2.4109735989028633d-2)*(1.0d0,0.0d0) + (1.7861402690121105d-2)*(0.0d0,1.0d0) 
        xe = (-8.0351990684851848d-1)*(1.0d0,0.0d0) + (5.9527788409964322d-1)*(0.0d0,1.0d0) 
     case(22)
        we = (-2.4987503287012616d-2)*(1.0d0,0.0d0) + (2.8488189217690068d-2)*(0.0d0,1.0d0) 
        xe = (-6.5940595627400689d-1)*(1.0d0,0.0d0) + (7.5178706082930324d-1)*(0.0d0,1.0d0) 
     case(23)
        we = (-1.9994108730100320d-2)*(1.0d0,0.0d0) + (4.0545468570332463d-2)*(0.0d0,1.0d0) 
        xe = (-4.4227617016573273d-1)*(1.0d0,0.0d0) + (8.9687891563105204d-1)*(0.0d0,1.0d0) 
     case(24)
        we = (-7.8174922097638110d-3)*(1.0d0,0.0d0) + (4.9162005002073170d-2)*(0.0d0,1.0d0) 
        xe = (-1.5704185299553061d-1)*(1.0d0,0.0d0) + (9.8759194833075170d-1)*(0.0d0,1.0d0) 
     case(25)
        we = (7.8174922097638110d-3)*(1.0d0,0.0d0) + (4.9162005002073170d-2)*(0.0d0,1.0d0) 
        xe = (1.5704185299553061d-1)*(1.0d0,0.0d0) + (9.8759194833075170d-1)*(0.0d0,1.0d0) 
     case(26)
        we = (1.9994108730100320d-2)*(1.0d0,0.0d0) + (4.0545468570332463d-2)*(0.0d0,1.0d0) 
        xe = (4.4227617016573273d-1)*(1.0d0,0.0d0) + (8.9687891563105204d-1)*(0.0d0,1.0d0) 
     case(27)
        we = (2.4987503287012616d-2)*(1.0d0,0.0d0) + (2.8488189217690068d-2)*(0.0d0,1.0d0) 
        xe = (6.5940595627400689d-1)*(1.0d0,0.0d0) + (7.5178706082930324d-1)*(0.0d0,1.0d0) 
     case(28)
        we = (2.4109735989028633d-2)*(1.0d0,0.0d0) + (1.7861402690121105d-2)*(0.0d0,1.0d0) 
        xe = (8.0351990684851848d-1)*(1.0d0,0.0d0) + (5.9527788409964322d-1)*(0.0d0,1.0d0) 
     case(29)
        we = (2.0411160688551314d-2)*(1.0d0,0.0d0) + (1.0418873853406705d-2)*(0.0d0,1.0d0) 
        xe = (8.9067323792710495d-1)*(1.0d0,0.0d0) + (4.5464401815095595d-1)*(0.0d0,1.0d0) 
     case(30)
        we = (1.6114837644635151d-2)*(1.0d0,0.0d0) + (5.8251138870511737d-3)*(0.0d0,1.0d0) 
        xe = (9.4044447327612435d-1)*(1.0d0,0.0d0) + (3.3994733810458494d-1)*(0.0d0,1.0d0) 
     case(31)
        we = (1.2254995789377469d-2)*(1.0d0,0.0d0) + (3.1800985784853704d-3)*(0.0d0,1.0d0) 
        xe = (9.6794166373051727d-1)*(1.0d0,0.0d0) + (2.5117510946468846d-1)*(0.0d0,1.0d0) 
     case(32)
        we = (9.1347263363820881d-3)*(1.0d0,0.0d0) + (1.7136119529180878d-3)*(0.0d0,1.0d0) 
        xe = (9.8285560240079251d-1)*(1.0d0,0.0d0) + (1.8437696393360870d-1)*(0.0d0,1.0d0) 
     case(33)
        we = (6.7365640264612303d-3)*(1.0d0,0.0d0) + (9.1691183433771899d-4)*(0.0d0,1.0d0) 
        xe = (9.9086381972995408d-1)*(1.0d0,0.0d0) + (1.3486619572808120d-1)*(0.0d0,1.0d0) 
     case(34)
        we = (4.9399091552901342d-3)*(1.0d0,0.0d0) + (4.8877131695888886d-4)*(0.0d0,1.0d0) 
        xe = (9.9514075291324122d-1)*(1.0d0,0.0d0) + (9.8462591329232213d-2)*(0.0d0,1.0d0) 
     case(35)
        we = (3.6116342963176664d-3)*(1.0d0,0.0d0) + (2.6002245345124118d-4)*(0.0d0,1.0d0) 
        xe = (9.9741833534204793d-1)*(1.0d0,0.0d0) + (7.1809918002307196d-2)*(0.0d0,1.0d0) 
     case(36)
        we = (2.6364847596615638d-3)*(1.0d0,0.0d0) + (1.3818198636717554d-4)*(0.0d0,1.0d0) 
        xe = (9.9862934321886543d-1)*(1.0d0,0.0d0) + (5.2339610834027481d-2)*(0.0d0,1.0d0) 
     case(37)
        we = (1.9232760208135970d-3)*(1.0d0,0.0d0) + (7.3391372929623125d-5)*(0.0d0,1.0d0) 
        xe = (9.9927271804559314d-1)*(1.0d0,0.0d0) + (3.8131810470687581d-2)*(0.0d0,1.0d0) 
     case(38)
        we = (1.4027639654332564d-3)*(1.0d0,0.0d0) + (3.8967906914869193d-5)*(0.0d0,1.0d0) 
        xe = (9.9961437632196248d-1)*(1.0d0,0.0d0) + (2.7768663101670359d-2)*(0.0d0,1.0d0) 
     case(39)
        we = (1.0234234700723219d-3)*(1.0d0,0.0d0) + (2.0687031904371676d-5)*(0.0d0,1.0d0) 
        xe = (9.9979576857513042d-1)*(1.0d0,0.0d0) + (2.0209431938186778d-2)*(0.0d0,1.0d0) 
     case(40)
        we = (7.4732098511848374d-4)*(1.0d0,0.0d0) + (1.0981144405960075d-5)*(0.0d0,1.0d0) 
        xe = (9.9989206046724421d-1)*(1.0d0,0.0d0) + (1.4692427116338712d-2)*(0.0d0,1.0d0) 
     case(41)
        we = (5.4669867904386235d-4)*(1.0d0,0.0d0) + (5.8285515797031256d-6)*(0.0d0,1.0d0) 
        xe = (9.9994317254293952d-1)*(1.0d0,0.0d0) + (1.0660754418012596d-2)*(0.0d0,1.0d0) 
     case(42)
        we = (4.0133002060141742d-4)*(1.0d0,0.0d0) + (3.0931700215852137d-6)*(0.0d0,1.0d0) 
        xe = (9.9997030010278387d-1)*(1.0d0,0.0d0) + (7.7070689855653837d-3)*(0.0d0,1.0d0) 
     case(43)
        we = (2.9653714734123929d-4)*(1.0d0,0.0d0) + (1.6407037867488601d-6)*(0.0d0,1.0d0) 
        xe = (9.9998469398327960d-1)*(1.0d0,0.0d0) + (5.5327930710187233d-3)*(0.0d0,1.0d0) 
     case(44)
        we = (2.2173073497403943d-4)*(1.0d0,0.0d0) + (8.6876512741129979d-7)*(0.0d0,1.0d0) 
        xe = (9.9999232430038310d-1)*(1.0d0,0.0d0) + (3.9180786512652823d-3)*(0.0d0,1.0d0) 
     case(45)
        we = (1.6934841332539216d-4)*(1.0d0,0.0d0) + (4.5718038200357800d-7)*(0.0d0,1.0d0) 
        xe = (9.9999635598178738d-1)*(1.0d0,0.0d0) + (2.6996338911749935d-3)*(0.0d0,1.0d0) 
     case(46)
        we = (1.3409354176612142d-4)*(1.0d0,0.0d0) + (2.3523148062411967d-7)*(0.0d0,1.0d0) 
        xe = (9.9999846133468917d-1)*(1.0d0,0.0d0) + (1.7542315281065331d-3)*(0.0d0,1.0d0) 
     case(47)
        we = (1.1240090501800907d-4)*(1.0d0,0.0d0) + (1.1085597222977934d-7)*(0.0d0,1.0d0) 
        xe = (9.9999951365073736d-1)*(1.0d0,0.0d0) + (9.8625467737379399d-4)*(0.0d0,1.0d0) 
     case(48)
        we = (1.0207666216851158d-4)*(1.0d0,0.0d0) + (3.2463391730449853d-8)*(0.0d0,1.0d0) 
        xe = (9.9999994942861647d-1)*(1.0d0,0.0d0) + (3.1802950252369703d-4)*(0.0d0,1.0d0) 

     end select



  case(56)
     select case(e)
        !------------------------   n = 56, S = 1.0e+06, rate = 3.33d-16  ------------------------
     case(0)
        we = (3.3306690738754696d-16)*(1.0d0,0.0d0) + (0.0000000000000000d0)*(0.0d0,1.0d0)  
     case(1)
        we = (-8.7205817070324285d-5)*(1.0d0,0.0d0) + (2.3745715284134577d-8)*(0.0d0,1.0d0) 
        xe = (-9.9999996292769566d-1)*(1.0d0,0.0d0) + (2.7229507390174712d-4)*(0.0d0,1.0d0) 
     case(2)
        we = (-9.3671614351072461d-5)*(1.0d0,0.0d0) + (7.8410128732040906d-8)*(0.0d0,1.0d0) 
        xe = (-9.9999964965324450d-1)*(1.0d0,0.0d0) + (8.3707430265302118d-4)*(0.0d0,1.0d0) 
     case(3)
        we = (-1.0708258208781082d-4)*(1.0d0,0.0d0) + (1.5676024389248878d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999892847206506d-1)*(1.0d0,0.0d0) + (1.4639175939086017d-3)*(0.0d0,1.0d0) 
     case(4)
        we = (-1.2843298795715179d-4)*(1.0d0,0.0d0) + (2.8246352114045597d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     case(5)
        we = (-1.5930565536197383d-4)*(1.0d0,0.0d0) + (4.9349124722739312d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999520196512837d-1)*(1.0d0,0.0d0) + (3.0977486537930887d-3)*(0.0d0,1.0d0) 
     case(6)
        we = (-2.0198918524432522d-4)*(1.0d0,0.0d0) + (8.5358791379195321d-7)*(0.0d0,1.0d0) 
        xe = (-9.9999107096619733d-1)*(1.0d0,0.0d0) + (4.2258712566424312d-3)*(0.0d0,1.0d0) 
     case(7)
        we = (-2.5964733455728753d-4)*(1.0d0,0.0d0) + (1.4715244324426047d-6)*(0.0d0,1.0d0) 
        xe = (-9.9998394069450114d-1)*(1.0d0,0.0d0) + (5.6673056293498336d-3)*(0.0d0,1.0d0) 
     case(8)
        we = (-3.3655289344432129d-4)*(1.0d0,0.0d0) + (2.5339479248362398d-6)*(0.0d0,1.0d0) 
        xe = (-9.9997165736303506d-1)*(1.0d0,0.0d0) + (7.5289089929910773d-3)*(0.0d0,1.0d0) 
     case(9)
        we = (-4.3840291264360994d-4)*(1.0d0,0.0d0) + (4.3617417606661609d-6)*(0.0d0,1.0d0) 
        xe = (-9.9995051075810371d-1)*(1.0d0,0.0d0) + (9.9486699918864014d-3)*(0.0d0,1.0d0) 
     case(10)
        we = (-5.7273754900334433d-4)*(1.0d0,0.0d0) + (7.5068928472465222d-6)*(0.0d0,1.0d0) 
        xe = (-9.9991411384269779d-1)*(1.0d0,0.0d0) + (1.3105912336512029d-2)*(0.0d0,1.0d0) 
     case(11)
        we = (-7.4949109140293982d-4)*(1.0d0,0.0d0) + (1.2919044514441024d-5)*(0.0d0,1.0d0) 
        xe = (-9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(12)
        we = (-9.8171039842577979d-4)*(1.0d0,0.0d0) + (2.2231801069936816d-5)*(0.0d0,1.0d0) 
        xe = (-9.9974367821792698d-1)*(1.0d0,0.0d0) + (2.2640182492421911d-2)*(0.0d0,1.0d0) 
     case(13)
        we = (-1.2864813335720538d-3)*(1.0d0,0.0d0) + (3.8254542936105366d-5)*(0.0d0,1.0d0) 
        xe = (-9.9955818427461929d-1)*(1.0d0,0.0d0) + (2.9722655494190512d-2)*(0.0d0,1.0d0) 
     case(14)
        we = (-1.6861033154932578d-3)*(1.0d0,0.0d0) + (6.5816167525760470d-5)*(0.0d0,1.0d0) 
        xe = (-9.9923902416215127d-1)*(1.0d0,0.0d0) + (3.9004776521238349d-2)*(0.0d0,1.0d0) 
     case(15)
        we = (-2.2095347979132953d-3)*(1.0d0,0.0d0) + (1.1320927252413921d-4)*(0.0d0,1.0d0) 
        xe = (-9.9868997900096679d-1)*(1.0d0,0.0d0) + (5.1169579273710977d-2)*(0.0d0,1.0d0) 
     case(16)
        we = (-2.8940726679749628d-3)*(1.0d0,0.0d0) + (1.9465226644598720d-4)*(0.0d0,1.0d0) 
        xe = (-9.9774576269728499d-1)*(1.0d0,0.0d0) + (6.7107324634595802d-2)*(0.0d0,1.0d0) 
     case(17)
        we = (-3.7870681227707027d-3)*(1.0d0,0.0d0) + (3.3445787196232090d-4)*(0.0d0,1.0d0) 
        xe = (-9.9612282741186675d-1)*(1.0d0,0.0d0) + (8.7973363633478865d-2)*(0.0d0,1.0d0) 
     case(18)
        we = (-4.9470931378947352d-3)*(1.0d0,0.0d0) + (5.7400507030525071d-4)*(0.0d0,1.0d0) 
        xe = (-9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(19)
        we = (-6.4430868535296079d-3)*(1.0d0,0.0d0) + (9.8314819657128148d-4)*(0.0d0,1.0d0) 
        xe = (-9.8855763469680935d-1)*(1.0d0,0.0d0) + (1.5084363719643493d-1)*(0.0d0,1.0d0) 
     case(20)
        we = (-8.3481034040599009d-3)*(1.0d0,0.0d0) + (1.6781503300194228d-3)*(0.0d0,1.0d0) 
        xe = (-9.8038756649065151d-1)*(1.0d0,0.0d0) + (1.9707922130589609d-1)*(0.0d0,1.0d0) 
     case(21)
        we = (-1.0720459860124735d-2)*(1.0d0,0.0d0) + (2.8477287489747280d-3)*(0.0d0,1.0d0) 
        xe = (-9.6648274241621845d-1)*(1.0d0,0.0d0) + (2.5673158865169970d-1)*(0.0d0,1.0d0) 
     case(22)
        we = (-1.3558373424360122d-2)*(1.0d0,0.0d0) + (4.7846871211399715d-3)*(0.0d0,1.0d0) 
        xe = (-9.4300389145487351d-1)*(1.0d0,0.0d0) + (3.3278170127121637d-1)*(0.0d0,1.0d0) 
     case(23)
        we = (-1.6705709101320253d-2)*(1.0d0,0.0d0) + (7.9063190596150665d-3)*(0.0d0,1.0d0) 
        xe = (-9.0388246359439139d-1)*(1.0d0,0.0d0) + (4.2778089252154988d-1)*(0.0d0,1.0d0) 
     case(24)
        we = (-1.9687576811665372d-2)*(1.0d0,0.0d0) + (1.2710819929445210d-2)*(0.0d0,1.0d0) 
        xe = (-8.4011858758926650d-1)*(1.0d0,0.0d0) + (5.4240276436151758d-1)*(0.0d0,1.0d0) 
     case(25)
        we = (-2.1505131459852337d-2)*(1.0d0,0.0d0) + (1.9556263525166893d-2)*(0.0d0,1.0d0) 
        xe = (-7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(26)
        we = (-2.0591627499276469d-2)*(1.0d0,0.0d0) + (2.8132163898937537d-2)*(0.0d0,1.0d0) 
        xe = (-5.9064294946629659d-1)*(1.0d0,0.0d0) + (8.0693302463448213d-1)*(0.0d0,1.0d0) 
     case(27)
        we = (-1.5387226218949800d-2)*(1.0d0,0.0d0) + (3.6764252330077289d-2)*(0.0d0,1.0d0) 
        xe = (-3.8608554905294024d-1)*(1.0d0,0.0d0) + (9.2246297964335111d-1)*(0.0d0,1.0d0) 
     case(28)
        we = (-5.7751262738652926d-3)*(1.0d0,0.0d0) + (4.2418103286963021d-2)*(0.0d0,1.0d0) 
        xe = (-1.3490312400076723d-1)*(1.0d0,0.0d0) + (9.9085879273226096d-1)*(0.0d0,1.0d0) 
     case(29)
        we = (5.7751262738652926d-3)*(1.0d0,0.0d0) + (4.2418103286963021d-2)*(0.0d0,1.0d0) 
        xe = (1.3490312400076723d-1)*(1.0d0,0.0d0) + (9.9085879273226096d-1)*(0.0d0,1.0d0) 
     case(30)
        we = (1.5387226218949800d-2)*(1.0d0,0.0d0) + (3.6764252330077289d-2)*(0.0d0,1.0d0) 
        xe = (3.8608554905294024d-1)*(1.0d0,0.0d0) + (9.2246297964335111d-1)*(0.0d0,1.0d0) 
     case(31)
        we = (2.0591627499276469d-2)*(1.0d0,0.0d0) + (2.8132163898937537d-2)*(0.0d0,1.0d0) 
        xe = (5.9064294946629659d-1)*(1.0d0,0.0d0) + (8.0693302463448213d-1)*(0.0d0,1.0d0) 
     case(32)
        we = (2.1505131459852337d-2)*(1.0d0,0.0d0) + (1.9556263525166893d-2)*(0.0d0,1.0d0) 
        xe = (7.3983485714849262d-1)*(1.0d0,0.0d0) + (6.7278851368618764d-1)*(0.0d0,1.0d0) 
     case(33)
        we = (1.9687576811665372d-2)*(1.0d0,0.0d0) + (1.2710819929445210d-2)*(0.0d0,1.0d0) 
        xe = (8.4011858758926650d-1)*(1.0d0,0.0d0) + (5.4240276436151758d-1)*(0.0d0,1.0d0) 
     case(34)
        we = (1.6705709101320253d-2)*(1.0d0,0.0d0) + (7.9063190596150665d-3)*(0.0d0,1.0d0) 
        xe = (9.0388246359439139d-1)*(1.0d0,0.0d0) + (4.2778089252154988d-1)*(0.0d0,1.0d0) 
     case(35)
        we = (1.3558373424360122d-2)*(1.0d0,0.0d0) + (4.7846871211399715d-3)*(0.0d0,1.0d0) 
        xe = (9.4300389145487351d-1)*(1.0d0,0.0d0) + (3.3278170127121637d-1)*(0.0d0,1.0d0) 
     case(36)
        we = (1.0720459860124735d-2)*(1.0d0,0.0d0) + (2.8477287489747280d-3)*(0.0d0,1.0d0) 
        xe = (9.6648274241621845d-1)*(1.0d0,0.0d0) + (2.5673158865169970d-1)*(0.0d0,1.0d0) 
     case(37)
        we = (8.3481034040599009d-3)*(1.0d0,0.0d0) + (1.6781503300194228d-3)*(0.0d0,1.0d0) 
        xe = (9.8038756649065151d-1)*(1.0d0,0.0d0) + (1.9707922130589609d-1)*(0.0d0,1.0d0) 
     case(38)
        we = (6.4430868535296079d-3)*(1.0d0,0.0d0) + (9.8314819657128148d-4)*(0.0d0,1.0d0) 
        xe = (9.8855763469680935d-1)*(1.0d0,0.0d0) + (1.5084363719643493d-1)*(0.0d0,1.0d0) 
     case(39)
        we = (4.9470931378947352d-3)*(1.0d0,0.0d0) + (5.7400507030525071d-4)*(0.0d0,1.0d0) 
        xe = (9.9333587640998278d-1)*(1.0d0,0.0d0) + (1.1525552757595411d-1)*(0.0d0,1.0d0) 
     case(40)
        we = (3.7870681227707027d-3)*(1.0d0,0.0d0) + (3.3445787196232090d-4)*(0.0d0,1.0d0) 
        xe = (9.9612282741186675d-1)*(1.0d0,0.0d0) + (8.7973363633478865d-2)*(0.0d0,1.0d0) 
     case(41)
        we = (2.8940726679749628d-3)*(1.0d0,0.0d0) + (1.9465226644598720d-4)*(0.0d0,1.0d0) 
        xe = (9.9774576269728499d-1)*(1.0d0,0.0d0) + (6.7107324634595802d-2)*(0.0d0,1.0d0) 
     case(42)
        we = (2.2095347979132953d-3)*(1.0d0,0.0d0) + (1.1320927252413921d-4)*(0.0d0,1.0d0) 
        xe = (9.9868997900096679d-1)*(1.0d0,0.0d0) + (5.1169579273710977d-2)*(0.0d0,1.0d0) 
     case(43)
        we = (1.6861033154932578d-3)*(1.0d0,0.0d0) + (6.5816167525760470d-5)*(0.0d0,1.0d0) 
        xe = (9.9923902416215127d-1)*(1.0d0,0.0d0) + (3.9004776521238349d-2)*(0.0d0,1.0d0) 
     case(44)
        we = (1.2864813335720538d-3)*(1.0d0,0.0d0) + (3.8254542936105366d-5)*(0.0d0,1.0d0) 
        xe = (9.9955818427461929d-1)*(1.0d0,0.0d0) + (2.9722655494190512d-2)*(0.0d0,1.0d0) 
     case(45)
        we = (9.8171039842577979d-4)*(1.0d0,0.0d0) + (2.2231801069936816d-5)*(0.0d0,1.0d0) 
        xe = (9.9974367821792698d-1)*(1.0d0,0.0d0) + (2.2640182492421911d-2)*(0.0d0,1.0d0) 
     case(46)
        we = (7.4949109140293982d-4)*(1.0d0,0.0d0) + (1.2919044514441024d-5)*(0.0d0,1.0d0) 
        xe = (9.9985147448075562d-1)*(1.0d0,0.0d0) + (1.7234528675274002d-2)*(0.0d0,1.0d0) 
     case(47)
        we = (5.7273754900334433d-4)*(1.0d0,0.0d0) + (7.5068928472465222d-6)*(0.0d0,1.0d0) 
        xe = (9.9991411384269779d-1)*(1.0d0,0.0d0) + (1.3105912336512029d-2)*(0.0d0,1.0d0) 
     case(48)
        we = (4.3840291264360994d-4)*(1.0d0,0.0d0) + (4.3617417606661609d-6)*(0.0d0,1.0d0) 
        xe = (9.9995051075810371d-1)*(1.0d0,0.0d0) + (9.9486699918864014d-3)*(0.0d0,1.0d0) 
     case(49)
        we = (3.3655289344432129d-4)*(1.0d0,0.0d0) + (2.5339479248362398d-6)*(0.0d0,1.0d0) 
        xe = (9.9997165736303506d-1)*(1.0d0,0.0d0) + (7.5289089929910773d-3)*(0.0d0,1.0d0) 
     case(50)
        we = (2.5964733455728753d-4)*(1.0d0,0.0d0) + (1.4715244324426047d-6)*(0.0d0,1.0d0) 
        xe = (9.9998394069450114d-1)*(1.0d0,0.0d0) + (5.6673056293498336d-3)*(0.0d0,1.0d0) 
     case(51)
        we = (2.0198918524432522d-4)*(1.0d0,0.0d0) + (8.5358791379195321d-7)*(0.0d0,1.0d0) 
        xe = (9.9999107096619733d-1)*(1.0d0,0.0d0) + (4.2258712566424312d-3)*(0.0d0,1.0d0) 
     case(52)
        we = (1.5930565536197383d-4)*(1.0d0,0.0d0) + (4.9349124722739312d-7)*(0.0d0,1.0d0) 
        xe = (9.9999520196512837d-1)*(1.0d0,0.0d0) + (3.0977486537930887d-3)*(0.0d0,1.0d0) 
     case(53)
        we = (1.2843298795715179d-4)*(1.0d0,0.0d0) + (2.8246352114045597d-7)*(0.0d0,1.0d0) 
        xe = (9.9999758153396057d-1)*(1.0d0,0.0d0) + (2.1993013049440135d-3)*(0.0d0,1.0d0) 
     case(54)
        we = (1.0708258208781082d-4)*(1.0d0,0.0d0) + (1.5676024389248878d-7)*(0.0d0,1.0d0) 
        xe = (9.9999892847206506d-1)*(1.0d0,0.0d0) + (1.4639175939086017d-3)*(0.0d0,1.0d0) 
     case(55)
        we = (9.3671614351072461d-5)*(1.0d0,0.0d0) + (7.8410128732040906d-8)*(0.0d0,1.0d0) 
        xe = (9.9999964965324450d-1)*(1.0d0,0.0d0) + (8.3707430265302118d-4)*(0.0d0,1.0d0) 
     case(56)
        we = (8.7205817070324285d-5)*(1.0d0,0.0d0) + (2.3745715284134577d-8)*(0.0d0,1.0d0) 
        xe = (9.9999996292769566d-1)*(1.0d0,0.0d0) + (2.7229507390174712d-4)*(0.0d0,1.0d0) 

     end select


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end select



end subroutine dset_zolotarev






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!2222222222222222 MAT-VEC MULTIPLICATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!! DOUBLE PRECISION VERSION


SUBROUTINE wdcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !
  !      REAL DOUBLE PRECISION version           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) DOUBLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  double precision :: alpha,beta
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*):: x
  double precision,dimension(*) ::b
!!!!!!!!
  character, dimension(6) :: matdescra

#ifdef MKL
  if(UPLO=='F') then
     matdescra(1)='G'
  else
     matdescra(1)='S' !M=N requirement
  end if
  matdescra(2)=UPLO
  matdescra(3)='N'
  matdescra(4)='F'
  if (TRANS=='T') then ! transpose
     call mkl_dcsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, N, beta, b, M)
  else
     call mkl_dcsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, M, beta, b, N)
  endif
#else
  !! default mat-vec in libsprim (if mkl not present)
  call dcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
#endif

end SUBROUTINE wdcsrmm






SUBROUTINE wzcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !      
  !      COMPLEX DOUBLE PRECISION version
  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) DOUBLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex (kind=kind(1.0d0)) :: alpha,beta
  complex (kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex (kind=kind(1.0d0)),dimension(*):: x
  complex (kind=kind(1.0d0)),dimension(*) ::b
!!!!!!!!
  character, dimension(6) :: matdescra

  !! mat-vec using intel mkl
#ifdef MKL
  if(UPLO=='F') then
     matdescra(1)='G'
  else
     matdescra(1)='S' !M=N requirement
  end if
  matdescra(2)=UPLO
  matdescra(3)='N'
  matdescra(4)='F'
  if (TRANS=='T'.or.TRANS=='C') then ! transpose
     call mkl_zcsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, N, beta, b, M)
  else
     call mkl_zcsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, M, beta, b, N)
  endif
#else
  !! default mat-vec in libsprim (if mkl not present)
  call zcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
#endif

end SUBROUTINE wzcsrmm



SUBROUTINE wzhcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !
  !      Note: A is complex Hermitian if UPLO='L or UPLO='U'     
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- hermitian CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- hermitian CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) COMPLEX DOUBLE PRECISION
  ! 
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  X     (input) COMPLEX DOUBLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  B     (input/output) COMPLEX DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex(kind=kind(1.0d0)) :: alpha,beta
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) :: x
  complex(kind=kind(1.0d0)),dimension(*) :: b
!!!!!!!!

  character, dimension(6) :: matdescra

  !! mat-vec using intel mkl
#ifdef MKL
  if(UPLO=='F') then
     matdescra(1)='G'
  else
     matdescra(1)='H' ! Hermitian Matrix  !M=N requirement
  end if
  matdescra(2)=UPLO
  matdescra(3)='N'
  matdescra(4)='F'
  if (TRANS=='T'.or.TRANS=='C') then ! transpose
     call mkl_zcsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, N, beta, b, M)
  else
     call mkl_zcsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, M, beta, b, N)
  endif
#else
  !! default mat-vec in libsprim (if mkl not present)
  call zhcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
#endif

end SUBROUTINE wzhcsrmm


!!!!!!!! SINGLE PRECISION VERSION

SUBROUTINE wscsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !
  !      REAL SINGLE PRECISION version             
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) SINGLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a SINGLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) SINGLE PRECISION
  !
  !  B     (input/output) SINGLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  real :: alpha,beta
  real,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  real,dimension(*):: x
  real,dimension(*) ::b

!!!!!!!!
  character, dimension(6) :: matdescra

  !!default mat-vec using intel mkl
#ifdef MKL
  if(UPLO=='F') then
     matdescra(1)='G'
  else
     matdescra(1)='S' !M=N requirement
  end if
  matdescra(2)=UPLO
  matdescra(3)='N'
  matdescra(4)='F'
  if (TRANS=='T') then ! transpose
     call mkl_scsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, N, beta, b, M)
  else
     call mkl_scsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, M, beta, b, N)
  endif
#else
  !! default mat-vec in libsprim (if mkl not present)
  call scsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
#endif


end SUBROUTINE wscsrmm



SUBROUTINE wccsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !      
  !      COMPLEX SINGLE PRECISION version
  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) SINGLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a SINGLE PRECISION)
  !
  !  X     (input) SINGLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) SINGLE PRECISION
  !
  !  B     (input/output) SINGLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex  :: alpha,beta
  complex ,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex ,dimension(*):: x
  complex ,dimension(*) ::b
!!!!!!!!

  character, dimension(6) :: matdescra

  !! mat-vec using intel mkl
#ifdef MKL
  if(UPLO=='F') then
     matdescra(1)='G'
  else
     matdescra(1)='S' !M=N requirement
  end if
  matdescra(2)=UPLO
  matdescra(3)='N'
  matdescra(4)='F'
  if (TRANS=='T'.or.TRANS=='C') then ! transpose
     call mkl_ccsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, N, beta, b, M)
  else
     call mkl_ccsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, M, beta, b, N)
  endif
#else
  !! default mat-vec in libsprim (if mkl not present)
  call ccsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
#endif


end SUBROUTINE wccsrmm





SUBROUTINE wchcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !
  !      Note: A is complex Hermitian if UPLO='L or UPLO='U'     
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric/hermitian CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric/hermitian CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) COMPLEX SINGLE PRECISION
  ! 
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX SINGLE PRECISION)
  !
  !  X     (input) COMPLEX SINGLE PRECISION
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) COMPLEX SINGLE PRECISION
  !
  !  B     (input/output) COMPLEX SINGLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex :: alpha,beta
  complex,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex,dimension(*) :: x
  complex,dimension(*) :: b
!!!!!!!!

  character, dimension(6) :: matdescra


  !! mat-vec using intel mkl
#ifdef MKL
  if(UPLO=='F') then
     matdescra(1)='G'
  else
     matdescra(1)='H' ! Hermitian Matrix  !M=N requirement
  end if
  matdescra(2)=UPLO
  matdescra(3)='N'
  matdescra(4)='F'
  if (TRANS=='T'.or.TRANS=='C') then ! transpose
     call mkl_ccsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, N, beta, b, M)
  else
     call mkl_ccsrmm(TRANS, N, rhs, M, alpha, matdescra, a, ja, ia, ia(2), X, M, beta, b, N)
  endif
#else
  !! default mat-vec in libsprim (if mkl not present)
  call chcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
#endif


end SUBROUTINE wchcsrmm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!3333333333333 ITERATIVE SOLVER 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!! DOUBLE PRECISION VERSION



subroutine  dbicgstab(opt,uplo,trans,sa,isa,jsa,N,M0,fj,xj,res,nbitmax,epso,comd,info) 
  implicit none
  integer :: opt
  character(len=1) :: uplo,trans
  double precision, dimension(*) :: sa
  integer, dimension(*) :: isa
  integer, dimension(*) :: jsa

  integer :: N,M0
   double precision, dimension(N,*) :: fj
   double precision, dimension(N,*):: xj
  double precision, dimension(M0) :: res

  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   double precision,dimension(:,:),allocatable :: work,work2!,temp1,temp2
  double precision,dimension(:),allocatable :: zdummy,sam
  integer :: ijob,j1,j2,k,i
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)

     allocate(zdummy(1:N))
     allocate(sam(1:isa(N+1)-1))

     do i=1,N
        do k=isa(i),isa(i+1)-1
           if (jsa(k)==i) zdummy(i)=sa(k)
        enddo
     enddo

     do i=1,N
        do k=isa(i),isa(i+1)-1
           ! scale matrix         
           sam(k)=sa(k)/(sqrt(zdummy(i))*sqrt(zdummy(jsa(k))))
        enddo
        ! scale rhs
        fj(i,1:M0)=fj(i,1:M0)/sqrt(zdummy(i))  !!!<< to optimize
     enddo

  end if
!!!!!!!!!!!!!!!!!!!!!


  allocate(work(N,7*M0))
  allocate(work2(M0,4))   



  ijob=-1
  do while (ijob/=0)
     call dbicgstab_rci(ijob,N,M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,info) 

     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input


     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)

!        print *,'workj1',sum(abs(work(1:N,j1)))
        
        if (opt==1) then ! jacobi preconditioner (2 sided)
           call wdcsrmm(UPLO,TRANS,N,N,M0,1.0d0,sam,isa,jsa,work(1,j1),0.0d0,work(1,j2))
        else
           call wdcsrmm(UPLO,TRANS,N,N,M0,1.0d0,sa,isa,jsa,work(1,j1),0.0d0,work(1,j2))
        end if
! print *,'workj2',sum(abs(work(1:N,j2)))!,work(1000,j2)
     end select
  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)
     do i=1,N
        ! scale solution
        xj(i,1:M0)=xj(i,1:M0)/sqrt(zdummy(i))
     enddo
     deallocate(zdummy)
     deallocate(sam)     
  end if
!!!!!!!!!!!!!!!!!!!!!

  deallocate(work,work2)


end subroutine dbicgstab



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  dbicgstab_rci(ijob,N,M0,fj,xj,work1,work2,j1,j2,res,nbitmax,epso,comd,info) 
  implicit none

  integer :: j1,j2 
  integer :: ijob
  integer :: N,M0
   double precision, dimension(N,*) :: fj,xj,work1
   double precision, dimension(M0,*):: work2 
  double precision, dimension(M0) :: res


  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   double precision, parameter :: ONE=1.0d0,ZERO=0.0d0

  double precision ::ares
  integer :: rank,i
  integer :: rj,rb,pp,pb,v,sb,t  !!! column index for the work1 array 
  integer:: rho_1,alpha,omega,ind !!! column index for the work2 array


  double precision,dimension(M0) :: resf
   double precision,dimension(M0) :: rho_2,aux0,aux1,beta
  integer,save :: nbit
  integer,save :: probe

  double precision :: ddot
  double precision :: DNRM2

  logical :: loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: b1,b2


  b1=1
  b2=M0

!!!!!!!!!! column index for the work1 array
  rj=1
  rb=M0+1
  pp=2*M0+1
  pb=3*M0+1
  v=4*M0+1
  sb=5*M0+1
  t=6*M0+1
!!!!!!!!!! column index for the work2 array
  rho_1=1
  alpha=2
  omega=3
  ind=4 ! if /= (0,0) do not include residual calculations in the max(norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rank=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !    rho_1(M0)   ! global  >> work2(:,1)
  !    rho_2(M0))  ! local
  !    beta(M0)    ! local
  !    alpha(M0)   ! global >> work2(:,2)
  !    omega(M0)   ! global >> work2(:,3)
  !    aux0(M0)    ! local
  !    aux1(M0)    ! local   
  !    resf(M0)    ! local (use res to store it)

!!!!!!!!!!!!

  info=0
  loop=.true.   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ijob==-1) then
!!!!!!!!!!!!!!! initialize index for residual calculation
     if (M0==1) then
        work2(1,ind)=ONE
     else
        do i=1,M0
           work2(i,ind)=res(i)*ONE
        enddo
!!! special case (user enter 0 everywhere ==> 1 everywhere)
        if (sum(work2(:,ind))==ZERO) work2(:,ind)=ONE

     end if

!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax      
!!!! let us temporaly save xj into work1(1,rj)
     call DCOPY(N*M0,xj(1,1),1,work1(1,rj),1)


!!!! Compute A*xj  ----- A*work1(1,rj)=>work1(1,rb)
     j1=rj
     j2=rb
     ijob=2
     probe=0
     return
     ! work1(1:N,rb:rb+M0-1)=ZEROC ! if initial guess 0
  end if


  if (probe==0) then

!!!! let us temporaly save fj into work1(1,rj)
     call DCOPY(N*M0,fj(1,1),1,work1(1,rj),1)

!!!!!!! compute residual r=f-A*x
     do i=1,M0   
        call DAXPY(N,-ONE,work1(1,i+rb-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!! CHECK THE  NORM of rj
     do i=1,M0
        res(i)=DNRM2(N,work1(1,i+rj-1),1)
        resf(i)=DNRM2(N,fj(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1))) !sqrt(sum(abs(work1(1:N,i+rj-1))**2))
        !resf(i)=sum(abs(fj(1:N,i+rj-1))) !sqrt(sum(abs(fj(1:N,i))**2)) ! cte all along
     enddo
  

     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))

     res=resf ! save resf value
     if (comd) print *,'rel. residual before iteration',ares,res(1),resf(1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! BIPCG-STAB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! choose rb
     call DCOPY(N*M0,work1(1,rj),1,work1(1,rb),1)

     nbit=0

     probe=1

  endif


  if (probe==1) then
     nbit=nbit+1
     if (nbit>1) rho_2=work2(1:M0,rho_1) 

!!!!!!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0
        work2(i,rho_1)=ddot(N,work1(1,rb+i-1),1,work1(1,rj+i-1),1)
        !work2(i,rho_1)=sum(work1(1:N,rb+i-1)*work1(1:N,rj+i-1))
     end do
    
     
!!!!!!!!!!!! TEST
     do i=1,M0
        if (work2(i,rho_1)==ZERO) then
           info= -201
           if ((comd).and.(rank==0)) print*,'ATTENTION ----- BICG-STAB OUT FAILED, rho_1=0 !!'
           nbitmax=nbit 
           ijob=0
           return              
        end if
     end do

!!!!!!!!!!!! CONDITION
     IF (nbit==1) THEN
        call DCOPY(N*M0,work1(1,rj),1,work1(1,pp),1)
     ELSE
        beta=(work2(1:M0,rho_1)/rho_2)*(work2(1:M0,alpha)/work2(1:M0,omega))
        do i=1,M0   !!!!!!<<< to optimize
           call DAXPY(N,-work2(i,omega),work1(1,v+i-1),1,work1(1,pp+i-1),1)
           call DSCAL(N,beta(i),work1(1,pp+i-1),1)
           call DAXPY(N,ONE,work1(1,i+rj-1),1,work1(1,i+pp-1),1)
           !pp(:,i)=r(:,i)+beta(i)*(pp(:,i)-omega(i)*v(:,i))
        end do
     END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve  M*pb=pp(k) ---> pb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call DCOPY(N*M0,work1(1,pp),1,work1(1,pb),1)
     j1=pb
     j2=pp
     ijob=1
     probe=2
     return
  endif

  if (probe==2) then  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by pb, results in v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=pb
     j2=v
     ijob=2
     probe=3
     return
  end if

  if (probe==3) then
!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0!! different RHS
        aux0(i)=ddot(N,work1(1,i+rb-1),1,work1(1,i+v-1),1)
        !aux0(i)=sum(work1(1:N,i+rb-1)*work1(1:N,i+v-1))
     end do

!!!!!!
     work2(1:M0,alpha)=work2(1:M0,rho_1)/aux0
     do i=1,M0   
        !ss(:,i)=r(:,i)-alpha(i)*v(:,i) !!!!!!!!! use rj to store ss !!!!!!!!!!!
        call DAXPY(N,-work2(i,alpha),work1(1,i+v-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!!! CHECK THE  NORM of ss
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DNRM2(N,work1(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1))) !sqrt(sum(abs(work1(1:N,i+rj-1))**2))
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))



     IF (ares<epso) then
        do i=1,M0  
           !      xj(:,i)=xj(:,i)+alpha(i)*pb(:,i) 
           call DAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        end do
        if ((comd).and.(rank==0)) print *,(nbit-1)*1.0d0+0.5d0,ares
        info=0
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares             
        return
     end IF

     res=resf ! save resf value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve M*sb=ss ---> sb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     call DCOPY(N*M0,work1(1,rj),1,work1(1,sb),1)
     j1=sb
     j2=rj
     ijob=1
     probe=4
     return
  end if


  if (probe==4) then          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by sb, results in t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=sb
     j2=t
     ijob=2
     probe=5
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! SCALAR PRODUCTS on different RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (probe==5) then

     do i=1,M0!! different RHS
        aux0(i)=ddot(N,work1(1,i+t-1),1,work1(1,i+rj-1),1)
        aux1(i)=ddot(N,work1(1,i+t-1),1,work1(1,i+t-1),1)
        !aux0(i)=sum(work1(1:N,i+t-1)*work1(1:N,i+rj-1))
        !aux1(i)=sum(work1(1:N,i+t-1)*work1(1:N,i+t-1))

     end do


!!!!!!
     work2(1:M0,omega)=aux0/aux1
     do i=1,M0
        if (work2(i,omega)==ZERO) then
           info= -202 
           if ((comd).and.(rank==0)) print *,'ATTENTION ----- BICG-STAB OUT FAILED, omega=0'
           nbitmax=nbit
           ijob=0
           return 
           !fail=.true.
        end if
     end do

     do i=1,M0
        ! xj(:,i)=xj(:,i)+alpha(i)*pb(:,i)+omega(i)*sb(:,i)
        call DAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        call DAXPY(N,work2(i,omega),work1(1,i+sb-1),1,xj(1,i),1) 
     end do

     do i=1,M0
        !   r(:,i)=ss(:,i)-omega(i)*t(:,i)
        !call ZCOPY(N,ss(1,i),1,rj(1,i),1)   !!! ss is rj
        call DAXPY(N,-work2(i,omega),work1(1,i+t-1),1,work1(1,i+rj-1),1) 
     end do


!!!!!!!!! CHECK THE  NORM of rj
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DNRM2(N,work1(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1)))!sqrt(sum(abs(work1(1:N,i+rj-1))**2))
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))
     if ((comd).and.(rank==0)) print *,nbit,ares


!!!! test for the other loop
     if (.not.((nbit<nbitmax).and.(ares.gt.epso))) loop=.false. 

     if (loop) then
        ijob=-2
        probe=1
        res=resf ! save resf 
     else
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares      
     end if
     return

  end if


end subroutine dbicgstab_rci


subroutine  zbicgstab(opt,uplo,trans,sa,isa,jsa,N,M0,fj,xj,res,nbitmax,epso,comd,info) 
  implicit none
  integer :: opt
  character(len=1) :: uplo,trans
  complex(kind=kind(1.0d0)), dimension(*) :: sa
  integer, dimension(*) :: isa
  integer, dimension(*) :: jsa

  integer :: N,M0
  complex(kind=kind(1.0d0)), dimension(N,*) :: fj
  complex(kind=kind(1.0d0)), dimension(N,*):: xj
  double precision, dimension(M0) :: res

  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: work,work2!,temp1,temp2
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: zdummy,sam
  integer :: ijob,j1,j2,k,i!,nnz,nbit0
  !character, dimension(6) :: matdescra
  !double precision,dimension(:),allocatable :: nres
  !complex(kind=kind(1.0d0)),dimension(:),allocatable :: dsa,alu
  !integer,dimension(:),allocatable :: jlu,ju,iw
  !integer :: ierr
  !double precision :: eps00
  !integer :: t1,t2,tim


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)

     allocate(zdummy(1:N))
     allocate(sam(1:isa(N+1)-1))

     do i=1,N
        do k=isa(i),isa(i+1)-1
           if (jsa(k)==i) zdummy(i)=sa(k)
        enddo
     enddo

     do i=1,N
        do k=isa(i),isa(i+1)-1
           ! scale matrix         
           sam(k)=sa(k)/(sqrt(zdummy(i))*sqrt(zdummy(jsa(k))))
        enddo
        ! scale rhs
        fj(i,1:M0)=fj(i,1:M0)/sqrt(zdummy(i))  !!!<< to optimize
     enddo

  end if
!!!!!!!!!!!!!!!!!!!!!


  allocate(work(N,7*M0))
  allocate(work2(M0,4))   



  ijob=-1
  do while (ijob/=0)
     call zbicgstab_rci(ijob,N,M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,info) 

     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input



        !print *,'call preconditioner for iteration'
        !! incomplete LU
!!$          call ilu0(n, sa, jsa, isa, alu, jlu, ju, iw, ierr)
!!$          if (ierr/=0) print *,'zero pivot!!!!'
!!$          do i=1,M0
!!$             !print *,i
!!$             call lusol(n, work(1,j2+i-1), work(1,j1+i-1), alu, jlu, ju)
!!$             
!!$          end do
!!$          
!!$                  print *,sum(work(:,j1+i-1)),sum(work(1,j1:j1+M0-1))
!!$          
        !          eps00=epso
        !    nbit0=3
        !   nres(1:M0)=work2(1:M0,4)
        !             call zbicgstabcopy(uplo,trans,dsa,isa,jsa,N,M0,work(1,j2),work(1,j1),nres,nbit0,eps00,comd,info)

        !            print *,'done' 
        !  call zSSOR(uplo,dsa,isa,jsa,n,m0,work(1,j1),work(1,j2),nres,nbit0)


     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)

!        print *,'workj1',sum(abs(work(1:N,j1)))
        
        if (opt==1) then ! jacobi preconditioner (2 sided)
           call wzcsrmm(UPLO,TRANS,N,N,M0,(1.0d0,0.0d0),sam,isa,jsa,work(1,j1),(0.0d0,0.0d0),work(1,j2))
        else
           call wzcsrmm(UPLO,TRANS,N,N,M0,(1.0d0,0.0d0),sa,isa,jsa,work(1,j1),(0.0d0,0.0d0),work(1,j2))
        end if
! print *,'workj2',sum(abs(work(1:N,j2)))!,work(1000,j2)
     end select
  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)
     do i=1,N
        ! scale solution
        xj(i,1:M0)=xj(i,1:M0)/sqrt(zdummy(i))
     enddo
     deallocate(zdummy)
     deallocate(sam)     
  end if
!!!!!!!!!!!!!!!!!!!!!

  deallocate(work,work2)


end subroutine zbicgstab



subroutine  zhbicgstab(opt,uplo,trans,sa,isa,jsa,N,M0,fj,xj,res,nbitmax,epso,comd,info) 
  implicit none
  integer :: opt
  character(len=1) :: uplo,trans
  complex(kind=kind(1.0d0)), dimension(*) :: sa
  integer, dimension(*) :: isa
  integer, dimension(*) :: jsa

  integer :: N,M0
  complex(kind=kind(1.0d0)), dimension(N,*) :: fj
  complex(kind=kind(1.0d0)), dimension(N,*):: xj
  double precision, dimension(M0) :: res

  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: work,work2!,temp1,temp2
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: zdummy,sam
  integer :: ijob,j1,j2,k,i!,nnz,nbit0
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)

     allocate(zdummy(1:N))
     allocate(sam(1:isa(N+1)-1))

     do i=1,N
        do k=isa(i),isa(i+1)-1
           if (jsa(k)==i) zdummy(i)=sa(k)
        enddo
     enddo

     do i=1,N
        do k=isa(i),isa(i+1)-1
           ! scale matrix         
           sam(k)=sa(k)/(sqrt(zdummy(i))*sqrt(zdummy(jsa(k))))
        enddo
        ! scale rhs
        fj(i,1:M0)=fj(i,1:M0)/sqrt(zdummy(i))  !!!<< to optimize
     enddo

  end if
!!!!!!!!!!!!!!!!!!!!!


  allocate(work(N,7*M0))
  allocate(work2(M0,4))   



  ijob=-1
  do while (ijob/=0)
     call zbicgstab_rci(ijob,N,M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,info) 

     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input



     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)

!        print *,'workj1',sum(abs(work(1:N,j1)))
        
        if (opt==1) then ! jacobi preconditioner (2 sided)
           call wzhcsrmm(UPLO,TRANS,N,N,M0,(1.0d0,0.0d0),sam,isa,jsa,work(1,j1),(0.0d0,0.0d0),work(1,j2))
        else
           call wzhcsrmm(UPLO,TRANS,N,N,M0,(1.0d0,0.0d0),sa,isa,jsa,work(1,j1),(0.0d0,0.0d0),work(1,j2))
        end if
     end select
  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)
     do i=1,N
        ! scale solution
        xj(i,1:M0)=xj(i,1:M0)/sqrt(zdummy(i))
     enddo
     deallocate(zdummy)
     deallocate(sam)     
  end if
!!!!!!!!!!!!!!!!!!!!!

  deallocate(work,work2)


end subroutine zhbicgstab


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  zbicgstab_rci(ijob,N,M0,fj,xj,work1,work2,j1,j2,res,nbitmax,epso,comd,info) 
  implicit none

  integer :: j1,j2 
  integer :: ijob
  integer :: N,M0
  complex(kind=kind(1.0d0)), dimension(N,*) :: fj
  complex(kind=kind(1.0d0)), dimension(N,*):: xj
  complex(kind=kind(1.0d0)), dimension(N,*):: work1
  complex(kind=kind(1.0d0)), dimension(M0,*):: work2 

  double precision, dimension(M0) :: res


  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=kind(1.0d0)), parameter :: ONEC=(1.0d0,0.0d0),ZEROC=(0.0d0,0.0d0)

  double precision ::ares
  integer :: rank,i
  integer :: rj,rb,pp,pb,v,sb,t  !!! column index for the work1 array 
  integer:: rho_1,alpha,omega,ind !!! column index for the work2 array


  double precision,dimension(M0) :: resf
  complex(kind=kind(1.0d0)),dimension(M0) :: rho_2,aux0,aux1,beta
  integer,save :: nbit
  integer,save :: probe

  !complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
  double precision :: DZNRM2

  logical :: loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: b1,b2


  b1=1
  b2=M0

!!!!!!!!!! column index for the work1 array
  rj=1
  rb=M0+1
  pp=2*M0+1
  pb=3*M0+1
  v=4*M0+1
  sb=5*M0+1
  t=6*M0+1
!!!!!!!!!! column index for the work2 array
  rho_1=1
  alpha=2
  omega=3
  ind=4 ! if /= (0,0) do not include residual calculations in the max(norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rank=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !    rho_1(M0)   ! global  >> work2(:,1)
  !    rho_2(M0))  ! local
  !    beta(M0)    ! local
  !    alpha(M0)   ! global >> work2(:,2)
  !    omega(M0)   ! global >> work2(:,3)
  !    aux0(M0)    ! local
  !    aux1(M0)    ! local   
  !    resf(M0)    ! local (use res to store it)

!!!!!!!!!!!!

  info=0
  loop=.true.   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ijob==-1) then
!!!!!!!!!!!!!!! initialize index for residual calculation
     if (M0==1) then
        work2(1,ind)=ONEC
     else
        do i=1,M0
           work2(i,ind)=res(i)*ONEC
        enddo
!!! special case (user enter 0 everywhere ==> 1 everywhere)
        if (sum(work2(:,ind))==ZEROC) work2(:,ind)=ONEC

     end if

!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax      
!!!! let us temporaly save xj into work1(1,rj)
     call ZCOPY(N*M0,xj(1,1),1,work1(1,rj),1)


!!!! Compute A*xj  ----- A*work1(1,rj)=>work1(1,rb)
     j1=rj
     j2=rb
     ijob=2
     probe=0
     return
     ! work1(1:N,rb:rb+M0-1)=ZEROC ! if initial guess 0
  end if


  if (probe==0) then

!!!! let us temporaly save fj into work1(1,rj)
     call ZCOPY(N*M0,fj(1,1),1,work1(1,rj),1)

!!!!!!! compute residual r=f-A*x
     do i=1,M0   
        call ZAXPY(N,-ONEC,work1(1,i+rb-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!! CHECK THE  NORM of rj
     do i=1,M0
        res(i)=DZNRM2(N,work1(1,i+rj-1),1)
        resf(i)=DZNRM2(N,fj(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1))) !sqrt(sum(abs(work1(1:N,i+rj-1))**2))
        !resf(i)=sum(abs(fj(1:N,i+rj-1))) !sqrt(sum(abs(fj(1:N,i))**2)) ! cte all along
     enddo
  

     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))

     res=resf ! save resf value
     if (comd) print *,'rel. residual before iteration',ares,res(1),resf(1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! BIPCG-STAB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! choose rb
     call ZCOPY(N*M0,work1(1,rj),1,work1(1,rb),1)

     nbit=0

     probe=1

  endif


  if (probe==1) then
     nbit=nbit+1
     if (nbit>1) rho_2=work2(1:M0,rho_1) 

!!!!!!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0
        !work2(i,rho_1)=zdotc(N,work1(1,rb+i-1),1,work1(1,rj+i-1),1)
         work2(i,rho_1)=dot_product(work1(1:N,rb+i-1),work1(1:N,rj+i-1))
        !work2(i,rho_1)=sum(work1(1:N,rb+i-1)*work1(1:N,rj+i-1))
     end do
    
     
!!!!!!!!!!!! TEST
     do i=1,M0
        if (work2(i,rho_1)==ZEROC) then
           info= -201
           if ((comd).and.(rank==0)) print*,'ATTENTION ----- BICG-STAB OUT FAILED, rho_1=0 !!'
           nbitmax=nbit 
           ijob=0
           return              
        end if
     end do

!!!!!!!!!!!! CONDITION
     IF (nbit==1) THEN
        call ZCOPY(N*M0,work1(1,rj),1,work1(1,pp),1)
     ELSE
        beta=(work2(1:M0,rho_1)/rho_2)*(work2(1:M0,alpha)/work2(1:M0,omega))
        do i=1,M0   !!!!!!<<< to optimize
           call ZAXPY(N,-work2(i,omega),work1(1,v+i-1),1,work1(1,pp+i-1),1)
           call ZSCAL(N,beta(i),work1(1,pp+i-1),1)
           call ZAXPY(N,ONEC,work1(1,i+rj-1),1,work1(1,i+pp-1),1)
           !pp(:,i)=r(:,i)+beta(i)*(pp(:,i)-omega(i)*v(:,i))
        end do
     END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve  M*pb=pp(k) ---> pb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call ZCOPY(N*M0,work1(1,pp),1,work1(1,pb),1)
     j1=pb
     j2=pp
     ijob=1
     probe=2
     return
  endif

  if (probe==2) then  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by pb, results in v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=pb
     j2=v
     ijob=2
     probe=3
     return
  end if

  if (probe==3) then
!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0!! different RHS
        !aux0(i)=zdotc(N,work1(1,i+rb-1),1,work1(1,i+v-1),1)
         aux0(i)=dot_product(work1(1:N,i+rb-1),work1(1:N,i+v-1))
        !aux0(i)=sum(work1(1:N,i+rb-1)*work1(1:N,i+v-1))
     end do

!!!!!!
     work2(1:M0,alpha)=work2(1:M0,rho_1)/aux0
     do i=1,M0   
        !ss(:,i)=r(:,i)-alpha(i)*v(:,i) !!!!!!!!! use rj to store ss !!!!!!!!!!!
        call ZAXPY(N,-work2(i,alpha),work1(1,i+v-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!!! CHECK THE  NORM of ss
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DZNRM2(N,work1(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1))) !sqrt(sum(abs(work1(1:N,i+rj-1))**2))
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))



     IF (ares<epso) then
        do i=1,M0  
           !      xj(:,i)=xj(:,i)+alpha(i)*pb(:,i) 
           call ZAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        end do
        if ((comd).and.(rank==0)) print *,(nbit-1)*1.0d0+0.5d0,ares
        info=0
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares             
        return
     end IF

     res=resf ! save resf value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve M*sb=ss ---> sb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     call ZCOPY(N*M0,work1(1,rj),1,work1(1,sb),1)
     j1=sb
     j2=rj
     ijob=1
     probe=4
     return
  end if


  if (probe==4) then          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by sb, results in t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=sb
     j2=t
     ijob=2
     probe=5
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! SCALAR PRODUCTS on different RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (probe==5) then

     do i=1,M0!! different RHS
        !aux0(i)=zdotc(N,work1(1,i+t-1),1,work1(1,i+rj-1),1)
        !aux1(i)=zdotc(N,work1(1,i+t-1),1,work1(1,i+t-1),1)
         aux0(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+rj-1))
        aux1(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+t-1))
        
        !aux0(i)=sum(work1(1:N,i+t-1)*work1(1:N,i+rj-1))
        !aux1(i)=sum(work1(1:N,i+t-1)*work1(1:N,i+t-1))

     end do


!!!!!!
     work2(1:M0,omega)=aux0/aux1
     do i=1,M0
        if (work2(i,omega)==ZEROC) then
           info= -202 
           if ((comd).and.(rank==0)) print *,'ATTENTION ----- BICG-STAB OUT FAILED, omega=0'
           nbitmax=nbit
           ijob=0
           return 
           !fail=.true.
        end if
     end do

     do i=1,M0
        ! xj(:,i)=xj(:,i)+alpha(i)*pb(:,i)+omega(i)*sb(:,i)
        call ZAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        call ZAXPY(N,work2(i,omega),work1(1,i+sb-1),1,xj(1,i),1) 
     end do

     do i=1,M0
        !   r(:,i)=ss(:,i)-omega(i)*t(:,i)
        !call ZCOPY(N,ss(1,i),1,rj(1,i),1)   !!! ss is rj
        call ZAXPY(N,-work2(i,omega),work1(1,i+t-1),1,work1(1,i+rj-1),1) 
     end do


!!!!!!!!! CHECK THE  NORM of rj
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DZNRM2(N,work1(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1)))!sqrt(sum(abs(work1(1:N,i+rj-1))**2))
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))
     if ((comd).and.(rank==0)) print *,nbit,ares


!!!! test for the other loop
     if (.not.((nbit<nbitmax).and.(ares.gt.epso))) loop=.false. 

     if (loop) then
        ijob=-2
        probe=1
        res=resf ! save resf 
     else
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares      
     end if
     return

  end if


end subroutine zbicgstab_rci




!!!!!!!! SINGLE PRECISION VERSION



subroutine  cbicgstab(opt,uplo,trans,sa,isa,jsa,N,M0,fj,xj,res,nbitmax,epso,comd,info) 
  implicit none
  integer :: opt
  character(len=1) :: uplo,trans
  complex, dimension(*) :: sa
  integer, dimension(*) :: isa
  integer, dimension(*) :: jsa

  integer :: N,M0
  complex, dimension(N,*) :: fj
  complex, dimension(N,*):: xj
  double precision, dimension(M0) :: res

  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  complex,dimension(:,:),allocatable :: work,work2
  complex,dimension(:),allocatable :: cdummy,sam
  integer :: ijob,j1,j2,k,i!,nnz,nbit0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)

     allocate(cdummy(1:N))
     allocate(sam(1:isa(N+1)-1))

     do i=1,N
        do k=isa(i),isa(i+1)-1
           if (jsa(k)==i) cdummy(i)=sa(k)
        enddo
     enddo

     do i=1,N
        do k=isa(i),isa(i+1)-1
           ! scale matrix         
           sam(k)=sa(k)/(sqrt(cdummy(i))*sqrt(cdummy(jsa(k))))
        enddo
        ! scale rhs
        fj(i,1:M0)=fj(i,1:M0)/sqrt(cdummy(i))
     enddo

  end if
!!!!!!!!!!!!!!!!!!!!!



  allocate(work(N,7*M0))
  allocate(work2(M0,4))

  !     allocate(nres(M0))



  ijob=-1
  do while (ijob/=0)
     call cbicgstab_rci(ijob,N,M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,info) 

     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input


     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)
!print *,'workj1',sum(abs(work(1:N,j1)))
        if (opt==1) then ! jacobi preconditioner (2 sided)           
           call wccsrmm(UPLO,TRANS,N,N,M0,(1.0e0,0.0e0),sam,isa,jsa,work(1,j1),(0.0e0,0.0e0),work(1,j2))   
        else
           call wccsrmm(UPLO,TRANS,N,N,M0,(1.0e0,0.0e0),sa,isa,jsa,work(1,j1),(0.0e0,0.0e0),work(1,j2))
        endif
!print *,'workj2',sum(abs(work(1:N,j2)))!

     end select

  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (opt==1) then ! jacobi preconditioner (2 sided)
     do i=1,N
        ! scale solution
        xj(i,1:M0)=xj(i,1:M0)/sqrt(cdummy(i))
     enddo
     deallocate(cdummy)
     deallocate(sam)
  end if
!!!!!!!!!!!!!!!!!!!!!



  deallocate(work,work2)

end subroutine cbicgstab



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  cbicgstab_rci(ijob,N,M0,fj,xj,work1,work2,j1,j2,res,nbitmax,epso,comd,info) 
  implicit none

  integer :: j1,j2 
  integer :: ijob
  integer :: N,M0
  complex, dimension(N,*) :: fj
  complex, dimension(N,*):: xj
  complex, dimension(N,*):: work1
  complex, dimension(M0,*):: work2 

  double precision, dimension(M0) :: res


  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex, parameter :: ONEC=(1.0e0,0.0e0),ZEROC=(0.0e0,0.0e0)

  double precision ::ares
  integer :: rank,i
  integer :: rj,rb,pp,pb,v,sb,t  !!! column index for the work1 array 
  integer:: rho_1,alpha,omega,ind !!! column index for the work2 array


  double precision,dimension(M0) :: resf
  complex,dimension(M0) :: rho_2,aux0,aux1,beta
  integer,save :: nbit
  integer,save :: probe


  !complex :: cdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
  real :: SCNRM2

  logical :: loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: b1,b2
  integer :: ii,tim,t1,t2


  b1=1
  b2=M0

!!!!!!!!!! column index for the work1 array
  rj=1
  rb=M0+1
  pp=2*M0+1
  pb=3*M0+1
  v=4*M0+1
  sb=5*M0+1
  t=6*M0+1
!!!!!!!!!! column index for the work2 array
  rho_1=1
  alpha=2
  omega=3
  ind=4 ! if /= (0,0) do not include residual calculations in the max(norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rank=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !    rho_1(M0)   ! global  >> work2(:,1)
  !    rho_2(M0))  ! local
  !    beta(M0)    ! local
  !    alpha(M0)   ! global >> work2(:,2)
  !    omega(M0)   ! global >> work2(:,3)
  !    aux0(M0)    ! local
  !    aux1(M0)    ! local   
  !    resf(M0)    ! local (use res to store it)

!!!!!!!!!!!!

  info=0
  loop=.true.   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ijob==-1) then
!!!!!!!!!!!!!!! initialize index for residual calculation
     if (M0==1) then
        work2(1,ind)=ONEC
     else
        do i=1,M0
           work2(i,ind)=res(i)*ONEC
        enddo
!!! special case (user enter 0 everywhere ==> 1 everywhere)
        if (sum(work2(:,ind))==ZEROC) work2(:,ind)=ONEC

     end if

!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax      
!!!! let us temporaly save xj into work1(1,rj)
     call CCOPY(N*M0,xj(1,1),1,work1(1,rj),1)

!!!! Compute A*xj  ----- A*work1(1,rj)=>work1(1,rb)
     j1=rj
     j2=rb
     ijob=2
     probe=0
     return
  end if


  if (probe==0) then

!!!! let us temporaly save fj into work1(1,rj)
     call CCOPY(N*M0,fj(1,1),1,work1(1,rj),1)

!!!!!!! compute residual r=f-A*x
     do i=1,M0   
        call CAXPY(N,-ONEC,work1(1,i+rb-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!! CHECK THE  NORM of rj
     do i=1,M0
        res(i)=SCNRM2(N,work1(1,i+rj-1),1)
        resf(i)=SCNRM2(N,fj(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1))) !sqrt(sum(abs(work1(1:N,i+rj-1))**2))
        !resf(i)=sum(abs(fj(1:N,i+rj-1))) !sqrt(sum(abs(fj(1:N,i))**2)) ! cte all along
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))

     res=resf ! save resf value
     if (comd) print *,'rel. residual before iteration',ares,res(1),resf(1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! BIPCG-STAB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! choose rb
     call CCOPY(N*M0,work1(1,rj),1,work1(1,rb),1)

     nbit=0

     probe=1

  endif


  if (probe==1) then
     nbit=nbit+1
     if (nbit>1) rho_2=work2(1:M0,rho_1) 


!!!!!!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0
        !work2(i,rho_1)=cdotc(N,work1(1,rb+i-1),1,work1(1,rj+i-1),1)
         work2(i,rho_1)=dot_product(work1(1:N,rb+i-1),work1(1:N,rj+i-1))
        !work2(i,rho_1)=sum(work1(1:N,rb+i-1)*work1(1:N,rj+i-1))
     end do

!!!!!!!!!!!! TEST
     do i=1,M0
        if (work2(i,rho_1)==ZEROC) then
           info= -201
           if ((comd).and.(rank==0)) print*,'ATTENTION ----- BICG-STAB OUT FAILED, rho_1=0 !!'
           ijob=0
           nbitmax=nbit 
           return              
        end if
     end do

!!!!!!!!!!!! CONDITION
     IF (nbit==1) THEN
        call CCOPY(N*M0,work1(1,rj),1,work1(1,pp),1)
     ELSE
        beta=(work2(1:M0,rho_1)/rho_2)*(work2(1:M0,alpha)/work2(1:M0,omega))
        do i=1,M0   !!!!!!<<< to optimize
           call CAXPY(N,-work2(i,omega),work1(1,v+i-1),1,work1(1,pp+i-1),1)
           call CSCAL(N,beta(i),work1(1,pp+i-1),1)
           call CAXPY(N,ONEC,work1(1,i+rj-1),1,work1(1,i+pp-1),1)
           !pp(:,i)=r(:,i)+beta(i)*(pp(:,i)-omega(i)*v(:,i))
        end do
     END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve M*pb=pp(k) ---> pb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call CCOPY(N*M0,work1(1,pp),1,work1(1,pb),1)
     j1=pb
     j2=pp
     ijob=1
     probe=2
     return
  endif

  if (probe==2) then  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by pb, results in v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=pb
     j2=v
     ijob=2
     probe=3
     return
  end if

  if (probe==3) then
!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0!! different RHS
        !aux0(i)=cdotc(N,work1(1,i+rb-1),1,work1(1,i+v-1),1)
         aux0(i)=dot_product(work1(1:N,i+rb-1),work1(1:N,i+v-1))
        !aux0(i)=sum(work1(1:N,i+rb-1)*work1(1:N,i+v-1))
     end do

!!!!!!
     work2(1:M0,alpha)=work2(1:M0,rho_1)/aux0
     do i=1,M0   
        !ss(:,i)=r(:,i)-alpha(i)*v(:,i) !!!!!!!!! use rj to store ss !!!!!!!!!!!
        call CAXPY(N,-work2(i,alpha),work1(1,i+v-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!!! CHECK THE  NORM of ss
     resf=res ! retrieve resf value
     do i=1,M0
        !res(i)=sum(abs((work1(1:N,i+rj-1)))) !sqrt(sum(abs(work1(1:N,i+rj-1))**2))
        res(i)=SCNRM2(N,work1(1,i+rj-1),1)
     enddo
     !enddo
     ! call system_clock(t2,tim)
     !print *,'time check norm ',(t2-t1)*1.0d0/tim

     !stop


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))



     IF (ares<epso) then
        do i=1,M0  
           !      xj(:,i)=xj(:,i)+alpha(i)*pb(:,i) 
           call CAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        end do
        if ((comd).and.(rank==0)) print *,(nbit-1)*1.0d0+0.5d0,ares
        info=0
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares             
        return
     end IF

     res=resf ! save resf value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve  M*sb=ss ---> sb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     call CCOPY(N*M0,work1(1,rj),1,work1(1,sb),1)
     j1=sb
     j2=rj
     ijob=1
     probe=4
     return
  end if


  if (probe==4) then          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by sb, results in t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=sb
     j2=t
     ijob=2
     probe=5
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! SCALAR PRODUCTS on different RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (probe==5) then
     !       print *,'b',work1(1:2,t),work1(1:2,rj)
     ! for system1
     ! b ( -1.22982194E-30,  2.56469992E-31) ( -8.80229881E-30,  1.61060142E-30) (  4.45166961E-13, -1.49033990E-14) (  2.93360995E-12, -5.34270097E-14)
     ! ==> in need of scaling or jacobi prec 
     !The "real*4" statement specifies the variable names to be single precision 4-byte re!al numbers which has 7 digits of accuracy and a magnitude range of 10 from -38 to +3!8. The "real" statement is the same as "real*4" statement in nearly all 32-bit compu!ters.

     do i=1,M0!! different RHS
        !aux0(i)=cdotc(N,work1(1,i+t-1),1,work1(1,i+rj-1),1)
        !aux1(i)=cdotc(N,work1(1,i+t-1),1,work1(1,i+t-1),1)
          aux0(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+rj-1))
        aux1(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+t-1))
        !aux0(i)=sum(work1(1:N,i+t-1)*work1(1:N,i+rj-1))
        !aux1(i)=sum(work1(1:N,i+t-1)*work1(1:N,i+t-1))          
     end do


!!!!!!
     work2(1:M0,omega)=aux0/aux1

     do i=1,M0
        if (work2(i,omega)==ZEROC) then
           info= -202 
           if ((comd).and.(rank==0)) print *,'ATTENTION ----- BICG-STAB OUT FAILED, omega=0'
           nbitmax=nbit
           ijob=0
           return 
           !fail=.true.
        end if
     end do

     do i=1,M0
        ! xj(:,i)=xj(:,i)+alpha(i)*pb(:,i)+omega(i)*sb(:,i)
        call CAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        call CAXPY(N,work2(i,omega),work1(1,i+sb-1),1,xj(1,i),1) 
     end do

     do i=1,M0
        !   r(:,i)=ss(:,i)-omega(i)*t(:,i)
        !call ZCOPY(N,ss(1,i),1,rj(1,i),1)   !!! ss is rj
        call CAXPY(N,-work2(i,omega),work1(1,i+t-1),1,work1(1,i+rj-1),1) 
     end do
     !print *,'a',work1(1:2,rj)

!!!!!!!!! CHECK THE  NORM of rj
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=SCNRM2(N,work1(1,i+rj-1),1)
        !res(i)=sum(abs(work1(1:N,i+rj-1)))!sqrt(sum(abs(work1(1:N,i+rj-1))**2))
     enddo

     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))
     if ((comd).and.(rank==0)) print *,nbit,ares


!!!! test for the other loop
     if (.not.((nbit<nbitmax).and.(ares.gt.epso))) loop=.false. 

     if (loop) then
        ijob=-2
        probe=1
        res=resf ! save resf 
     else
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares      
     end if
     return

  end if


end subroutine cbicgstab_rci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine  dqrgs(N,M,q,r)
  !! Purpose: Compute QR factorization using Gram-Schmidt algorithm
  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) REAL DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) REAL DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
  integer :: N,M
  double precision,dimension(N,*) :: q
  double precision,dimension(M,*) :: R
  !---
  double precision :: DNRM2
  double precision :: ddot
  integer :: k,j
   
do k=1,M
   R(k,k)=DNRM2(N,q(1,k),1)
   call DSCAL(N,1.0d0/R(k,k), Q(1,k), 1)
   do j=k+1,M
      R(k,j)=ddot(N,Q(1,k),1,Q(1,j),1)
      call DAXPY(N,-R(k,j),Q(1,k),1,Q(1,j),1)     
enddo
enddo



end subroutine dqrgs



subroutine  zqrgs(N,M,q,r)
  !! Purpose: Compute QR factorization using Gram-Schmidt algorithm
  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
  integer :: N,M
  complex(kind=kind(1.0d0)),dimension(N,*) :: q
  complex(kind=kind(1.0d0)),dimension(M,*) :: R
  !---
  double precision :: DZNRM2
  !complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
  integer :: k,j
   
do k=1,M
   R(k,k)=(1.0d0,0.0d0)*DZNRM2(N,q(1,k),1)
   call ZSCAL(N,(1.0d0,0.0d0)/R(k,k), Q(1,k), 1)
!R(k,k)=sqrt(sum(abs(q(1:N,k))**2))
!Q(1:N,k)=Q(1:N,k)/R(k,k)
   do j=k+1,M
      !R(k,j)=zdotc(N,Q(1,k),1,Q(1,j),1)
      R(k,j)=dot_product(Q(1:N,k),Q(1:N,j))
      call ZAXPY(N,-R(k,j),Q(1,k),1,Q(1,j),1)     
     ! R(k,j)=sum(conjg(Q(1:N,k))*Q(1:N,j))
!Q(1:N,j)=Q(1:N,j)-R(k,j)*Q(1:N,k)
enddo
enddo



end subroutine zqrgs









subroutine  darnoldi(UPLO,N,sa,isa,jsa,alpha,X,res,epse,itmax)
  !! Purpose: Lanczos algorithm 


  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa,alpha,res
  integer,dimension(*) :: isa,jsa
  double precision,dimension(N,*) :: X
  integer :: itmax
  double precision :: epse
!!!
  integer,dimension(4) :: iseed
   double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   double precision,dimension(:),allocatable :: work,u,beta,de,rese
    double precision,dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye
    integer :: i,j,lwork,info_lap,it
    double precision :: DNRM2,ddot,temp
     logical :: conv
  integer :: mine,maxe



allocate(de(2))
  allocate(rese(2))
  allocate(Xe(N,2))
  allocate(Ye(N,2))
 
    
   allocate(u(1:N))
   allocate(V(n,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors
 
!!! initialization
   !v1=random
     iseed=(/45,8,89,33/)
     call DLARNV(1,iseed , n, v(1,1) ) ! random 0-1:  uniform
    v(:,1)=v(:,1)/DNRM2(N,v(1,1),1) ! normalization


!!!!!!!!!!!!!! Arnoldi using Krylov subspace
! do i=2,itmax   
! call wdcsrmm(UPLO,'N',N,N,1,DONE,sa,isa,jsa,v(1,i-1),DZERO,v(1,i))
!end do
!call  dqrgs(N,itmax,V,R)
! call wdcsrmm(UPLO,'N',N,N,itmax,DONE,sa,isa,jsa,V,DZERO,X(1,1))
! H=matmul(transpose(V),X(1:N,1:itmax)) 


!!!!!!!!!!!!!!!!!!! Arnoldi direct
    H=DZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       call wdcsrmm(UPLO,'N',N,N,1,DONE,sa,isa,jsa,v(1,it),DZERO,u)
       do j=1,it
          H(j,it)=ddot(N,V(1,j),1,u,1)
          call DAXPY(N,-H(j,it),V(1,j),1,u,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DNRM2(N,u,1)
          call DCOPY(N,u,1,V(1,it+1),1)
          call DSCAL(N,DONE/H(it+1,it), v(1,it+1), 1)
       end if

  if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call DHSEQR('E', 'I', it, 1, it, H_temp, itmax, alpha, beta, Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (alpha(j)>alpha(maxe)) maxe=j
           if (alpha(j)<alpha(mine)) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call DGEMM('N', 'N', N, 1, it, DONE, V(1,1), N, Xred(1,mine), itmax, DZERO, Xe(1,1), N)
        call DGEMM('N', 'N', N, 1, it, DONE, V(1,1), N, Xred(1,maxe), itmax, DZERO, Xe(1,2), N)
        dE(1)=alpha(mine)
        dE(2)=alpha(maxe)

        call DLACPY('F',N,2,Xe,N,Ye,N)
        do j=1,2
           call DSCAL(N,dE(j), Ye(1,j), 1)
           rese(j)=DNRM2(N,Ye(1,j),1) ! placeholder denominateur
        enddo
        call wdcsrmm(UPLO,'N',N,N,2,DONE,sa,isa,jsa,Xe,-DONE,Ye)
 do j=1,2
           rese(j)=DNRM2(N,Ye(1,j),1)/rese(j)
           !   if (j==2)
          ! print *,it,j,de(j),rese(j)
enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do
       
     
itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)<alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo

!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
        call DGEMM('N', 'N', N, itmax, itmax, DONE, V(1,1), N, Xred, itmax, DZERO, X, N)
 
 

  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
      
   deallocate(u) 
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine darnoldi




subroutine  zarnoldi(UPLO,N,sa,isa,jsa,alpha,X,res,epse,itmax)
!! check amplitude of eigenvalues....

  !! Purpose: Lanczos algorithm 


  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
  character(len=1) :: UPLO
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) :: sa,alpha
  double precision,dimension(*) :: res
  integer,dimension(*) :: isa,jsa
  complex(kind=kind(1.0d0)),dimension(N,*) :: X
  integer :: itmax
  double precision :: epse
!!!
  integer,dimension(4) :: iseed
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   complex(kind=kind(1.0d0)),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
   complex(kind=kind(1.0d0)),dimension(:),allocatable :: work,u,beta,de
     double precision,dimension(:),allocatable :: rese
     complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye
    integer :: i,j,lwork,info_lap,it
    double precision :: DZNRM2,temp
    !complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
     logical :: conv
  integer :: mine,maxe



allocate(de(2))
  allocate(rese(2))
  allocate(Xe(N,2))
  allocate(Ye(N,2))
 
    
   allocate(u(1:N))
   allocate(V(n,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors
 
!!! initialization
   !v1=random
     iseed=(/45,8,89,33/)
     call ZLARNV(1,iseed , n, v(1,1) ) ! random 0-1:  uniform
    v(:,1)=v(:,1)/DZNRM2(N,v(1,1),1) ! normalization


!!!!!!!!!!!!!! Arnoldi using Krylov subspace
! do i=2,itmax   
! call wdcsrmm(UPLO,'N',N,N,1,DONE,sa,isa,jsa,v(1,i-1),DZERO,v(1,i))
!end do
!call  dqrgs(N,itmax,V,R)
! call wdcsrmm(UPLO,'N',N,N,itmax,DONE,sa,isa,jsa,V,DZERO,X(1,1))
! H=matmul(transpose(V),X(1:N,1:itmax)) 


!!!!!!!!!!!!!!!!!!! Arnoldi direct
    H=ZZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       call wzcsrmm(UPLO,'N',N,N,1,ZONE,sa,isa,jsa,v(1,it),ZZERO,u)
       do j=1,it
          !H(j,it)=zdotc(N,V(1,j),1,u,1)
          H(j,it)=dot_product(V(1:N,j),u(1:N))
          call ZAXPY(N,-H(j,it),V(1,j),1,u,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DZNRM2(N,u,1)
          call ZCOPY(N,u,1,V(1,it+1),1)
          call ZSCAL(N,ZONE/H(it+1,it), v(1,it+1), 1)
       end if

  if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call ZHSEQR('E', 'I', it, 1, it, H_temp, itmax, alpha, Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (abs(alpha(j)) > abs(alpha(maxe))) maxe=j
           if (abs(alpha(j)) < abs(alpha(mine))) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call DGEMM('N', 'N', N, 1, it, DONE, V(1,1), N, Xred(1,mine), itmax, DZERO, Xe(1,1), N)
        call DGEMM('N', 'N', N, 1, it, DONE, V(1,1), N, Xred(1,maxe), itmax, DZERO, Xe(1,2), N)
        dE(1)=alpha(mine)
        dE(2)=alpha(maxe)

        call ZLACPY('F',N,2,Xe,N,Ye,N)
        do j=1,2
           call ZSCAL(N,dE(j), Ye(1,j), 1)
           rese(j)=DZNRM2(N,Ye(1,j),1) ! placeholder denominateur
        enddo
        call wzcsrmm(UPLO,'N',N,N,2,ZONE,sa,isa,jsa,Xe,-ZONE,Ye)
 do j=1,2
           rese(j)=DZNRM2(N,Ye(1,j),1)/rese(j)
           !   if (j==2)
          ! print *,it,j,de(j),rese(j)
enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do
       
     
itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
do i=1,itmax
   do j=i+1,itmax
      if (abs(alpha(j))< abs(alpha(i))) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo

!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
        call ZGEMM('N', 'N', N, itmax, itmax, ZONE, V(1,1), N, Xred, itmax, ZZERO, X, N)
 
 

  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
      
   deallocate(u) 
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine zarnoldi






subroutine  zharnoldi(UPLO,N,sa,isa,jsa,alpha,X,res,epse,itmax)

  !! check and sort real part of eigenvalues....real eigenvalue here
  
  !! Purpose: Lanczos algorithm 


  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !===================================================================h=  
  implicit none
  character(len=1) :: UPLO
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) :: sa
  double precision,dimension(*) :: res,alpha
  integer,dimension(*) :: isa,jsa
  complex(kind=kind(1.0d0)),dimension(N,*) :: X
  integer :: itmax
  double precision :: epse
!!!
  integer,dimension(4) :: iseed
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   complex(kind=kind(1.0d0)),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
   complex(kind=kind(1.0d0)),dimension(:),allocatable :: work,u,beta,zalpha,dE
    double precision,dimension(:),allocatable :: rese
     complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye
    integer :: i,j,lwork,info_lap,it
    double precision :: DZNRM2,temp
    !complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
     logical :: conv
  integer :: mine,maxe



allocate(dE(2))
  allocate(rese(2))
  allocate(Xe(N,2))
  allocate(Ye(N,2))
 allocate(zalpha(itmax))
    
   allocate(u(1:N))
   allocate(V(n,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors
 
!!! initialization
   !v1=random
     iseed=(/45,8,89,33/)
     call ZLARNV(1,iseed , n, v(1,1) ) ! random 0-1:  uniform
    v(:,1)=v(:,1)/DZNRM2(N,v(1,1),1) ! normalization


!!!!!!!!!!!!!! Arnoldi using Krylov subspace
! do i=2,itmax   
! call wdcsrmm(UPLO,'N',N,N,1,DONE,sa,isa,jsa,v(1,i-1),DZERO,v(1,i))
!end do
!call  dqrgs(N,itmax,V,R)
! call wdcsrmm(UPLO,'N',N,N,itmax,DONE,sa,isa,jsa,V,DZERO,X(1,1))
! H=matmul(transpose(V),X(1:N,1:itmax)) 


!!!!!!!!!!!!!!!!!!! Arnoldi direct
    H=ZZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       call wzhcsrmm(UPLO,'N',N,N,1,ZONE,sa,isa,jsa,v(1,it),ZZERO,u)
       do j=1,it
          !H(j,it)=zdotc(N,V(1,j),1,u,1)
           H(j,it)=dot_product(V(1:N,j),u(1:N))
          call ZAXPY(N,-H(j,it),V(1,j),1,u,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DZNRM2(N,u,1)
          call ZCOPY(N,u,1,V(1,it+1),1)
          call ZSCAL(N,ZONE/H(it+1,it), v(1,it+1), 1)
       end if

  if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call ZHSEQR('E', 'I', it, 1, it, H_temp, itmax, zalpha, Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (dble(zalpha(j))> dble(zalpha(maxe))) maxe=j
           if (dble(zalpha(j))< dble(zalpha(mine))) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call ZGEMM('N', 'N', N, 1, it, ZONE, V(1,1), N, Xred(1,mine), itmax, ZZERO, Xe(1,1), N)
        call ZGEMM('N', 'N', N, 1, it, ZONE, V(1,1), N, Xred(1,maxe), itmax, ZZERO, Xe(1,2), N)
        dE(1)=zalpha(mine)
        dE(2)=zalpha(maxe)

        call ZLACPY('F',N,2,Xe,N,Ye,N)
        do j=1,2
           call ZSCAL(N,dE(j), Ye(1,j), 1)
           rese(j)=DZNRM2(N,Ye(1,j),1) ! placeholder denominateur
        enddo
        call wzhcsrmm(UPLO,'N',N,N,2,ZONE,sa,isa,jsa,Xe,-ZONE,Ye)
 do j=1,2
           rese(j)=DZNRM2(N,Ye(1,j),1)/rese(j)
           !   if (j==2)
!           print *,it,j,de(j),rese(j)
enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do
       
     
itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
alpha(1:itmax)=dble(zalpha(1:itmax))
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)< alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo

!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
        call ZGEMM('N', 'N', N, itmax, itmax, ZONE, V(1,1), N, Xred, itmax, ZZERO, X, N)
 

  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
  deallocate(zalpha)
      
   deallocate(u) 
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine zharnoldi





subroutine  dgarnoldi(UPLO,N,sa,isa,jsa,sb,isb,jsb,alpha,X,res,epse,itmax,epsb,itmaxb)
  !! Purpose: Lanczos algorithm 

  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
  character(len=1) :: UPLO
  integer :: N
  double precision,dimension(*) :: sa,alpha,res,sb
  integer,dimension(*) :: isa,jsa,isb,jsb
  double precision,dimension(N,*) :: X
  double precision :: epse,epsb
  integer :: itmax,itmaxb
!!!
  integer,dimension(4) :: iseed
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  double precision,dimension(:),allocatable :: work,u,w,beta,rese,de
  double precision,dimension(:,:),allocatable :: V,H,H_temp,Xred,aux,Xe,Ye
  integer :: i,j,it,lwork,info_lap,infoloc,linloops
  double precision :: DNRM2,ddot,temp,eps
  logical :: comb,conv
  integer :: mine,maxe


  allocate(de(2))
  allocate(rese(2))
  allocate(Xe(N,2))
  allocate(Ye(N,2))
  allocate(aux(N,2))

  allocate(u(1:N))
  allocate(w(1:N))
  allocate(V(n,itmax)) ! subspace
  allocate(H(itmax,itmax)) ! hessenberg matrix
  allocate(H_temp(itmax,itmax)) ! hessenberg matrix
  lwork=6*itmax!itmax!3*itmax-1
  allocate(work(lwork))
  allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
  allocate(beta(itmax)) ! imaginary part of vectors

!!! initialization
  !v1=random
  iseed=(/45,8,89,33/)
  call DLARNV(1,iseed , n, v(1,1) ) ! random 0-1:  uniform
  v(:,1)=v(:,1)/DNRM2(N,v(1,1),1) ! normalization

  comb=.false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Arnoldi direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  H=DZERO
  conv=.false.
  it=0
  do while ((.not.conv).and.(it<=itmax-1))
     it=it+1
     ! mat-vec
     call wdcsrmm(UPLO,'N',N,N,1,DONE,sa,isa,jsa,v(1,it),DZERO,w) 
     ! solve
     linloops=itmaxb
     eps=epsb
     u(1:n)=DZERO !initial guess
     call dbicgstab(0,UPLO,'N',sb,isb,jsb,N,1,w,u,res,linloops,eps,comb,infoloc)
     ! print *,linloops,eps
     do j=1,it
        H(j,it)=ddot(N,V(1,j),1,u,1)
        call DAXPY(N,-H(j,it),V(1,j),1,u,1)
     end do
    
     if (it<itmax) then
        H(it+1,it)=DNRM2(N,u,1)
        call DCOPY(N,u,1,V(1,it+1),1)
        call DSCAL(N,DONE/H(it+1,it), v(1,it+1), 1)
     end if

  if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call DHSEQR('E', 'I', it, 1, it, H_temp, itmax, alpha, beta, Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (alpha(j)>alpha(maxe)) maxe=j
           if (alpha(j)<alpha(mine)) mine=j
        enddo

!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
        call DGEMM('N', 'N', N, 1, it, DONE, V(1,1), N, Xred(1,mine), itmax, DZERO, Xe(1,1), N)
        call DGEMM('N', 'N', N, 1, it, DONE, V(1,1), N, Xred(1,maxe), itmax, DZERO, Xe(1,2), N)
        dE(1)=alpha(mine)
        dE(2)=alpha(maxe)

        call DLACPY('F',N,2,Xe,N,Ye,N)
        do j=1,2
           call DSCAL(N,dE(j), Ye(1,j), 1)
           rese(j)=DNRM2(N,Ye(1,j),1) ! placeholder denominateur
        enddo
        call wdcsrmm(UPLO,'N',N,N,2,DONE,sa,isa,jsa,Xe,DZERO,aux)

        linloops=itmaxb
        eps=epsb
        Xe=DZERO !initial guess
        call dbicgstab(0,UPLO,'N',sb,isb,jsb,N,2,aux,Xe,res,linloops,eps,comb,infoloc)
        do j=1,2
           call DAXPY(N,-DONE,Xe(1,j),1,Ye(1,j),1)
           rese(j)=DNRM2(N,Ye(1,j),1)/rese(j)
           !if (j==2)
          ! print *,it,j,de(j),rese(j)
        enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if

  end do

 ! itmax=50
 !  return
  
  itmax=it
   
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)<alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo

!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
     call DGEMM('N', 'N', N, itmax, itmax, DONE, V(1,1), N, Xred, itmax, DZERO, X, N)

  

  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
  deallocate(aux)


  deallocate(u)
  deallocate(w) 
  deallocate(V) 
  deallocate(H)
  deallocate(H_temp)
  deallocate(work)
  deallocate(Xred)
  deallocate(beta)

end subroutine dgarnoldi



subroutine  zhgarnoldi(UPLO,N,sa,isa,jsa,sb,isb,jsb,alpha,X,res,epse,itmax,epsb,itmaxb)
  !! Return Eigenvalue real

  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
  character(len=1) :: UPLO
  integer :: N
   complex(kind=kind(1.0d0)),dimension(*) :: sa,sb
  double precision,dimension(*) :: alpha,res
  integer,dimension(*) :: isa,jsa,isb,jsb
   complex(kind=kind(1.0d0)),dimension(N,*) :: X
  double precision :: epse,epsb
  integer :: itmax,itmaxb
!!!
  integer,dimension(4) :: iseed
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
    complex(kind=kind(1.0d0)),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
    complex(kind=kind(1.0d0)),dimension(:),allocatable :: work,u,w,beta,dE,zalpha
   double precision,dimension(:),allocatable ::  rese
   complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: V,H,H_temp,Xred,aux,Xe,Ye
  integer :: i,j,it,lwork,info_lap,infoloc,linloops
  double precision :: DZNRM2,temp,eps
  !complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
  logical :: comb,conv
  integer :: mine,maxe


  allocate(de(2))
  allocate(rese(2))
  allocate(Xe(N,2))
  allocate(Ye(N,2))
  allocate(aux(N,2))
  allocate(zalpha(itmax))
  
  allocate(u(1:N))
  allocate(w(1:N))
  allocate(V(n,itmax)) ! subspace
  allocate(H(itmax,itmax)) ! hessenberg matrix
  allocate(H_temp(itmax,itmax)) ! hessenberg matrix
  lwork=6*itmax!itmax!3*itmax-1
  allocate(work(lwork))
  allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
  allocate(beta(itmax)) ! imaginary part of vectors

!!! initialization
  !v1=random
  iseed=(/45,8,89,33/)
  call ZLARNV(1,iseed , n, v(1,1) ) ! random 0-1:  uniform
  v(:,1)=v(:,1)/DZNRM2(N,v(1,1),1) ! normalization

  comb=.false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Arnoldi direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  H=ZZERO
  conv=.false.
  it=0
  do while ((.not.conv).and.(it<=itmax-1))
     it=it+1
     ! mat-vec
     call wzhcsrmm(UPLO,'N',N,N,1,ZONE,sa,isa,jsa,v(1,it),ZZERO,w)  
     ! solve
     linloops=itmaxb
     eps=epsb
     u(1:n)=ZZERO !initial guess
      call zhbicgstab(0,UPLO,'N',sb,isb,jsb,N,1,w,u,res,linloops,eps,comb,infoloc)
     do j=1,it
        !H(j,it)=zdotc(N,V(1,j),1,u,1)
        H(j,it)=dot_product(V(1:N,j),u(1:N))
        call ZAXPY(N,-H(j,it),V(1,j),1,u,1)
     end do
    
     if (it<itmax) then
        H(it+1,it)=DZNRM2(N,u,1)
        call ZCOPY(N,u,1,V(1,it+1),1)
        call ZSCAL(N,ZONE/H(it+1,it), v(1,it+1), 1)
     end if

  if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call ZHSEQR('E', 'I', it, 1, it, H_temp, itmax, zalpha, Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (dble(zalpha(j))> dble(zalpha(maxe))) maxe=j
           if (dble(zalpha(j))< dble(zalpha(mine))) mine=j
        enddo

!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
        call ZGEMM('N', 'N', N, 1, it, ZONE, V(1,1), N, Xred(1,mine), itmax, ZZERO, Xe(1,1), N)
        call ZGEMM('N', 'N', N, 1, it, ZONE, V(1,1), N, Xred(1,maxe), itmax, ZZERO, Xe(1,2), N)
        dE(1)=zalpha(mine)
        dE(2)=zalpha(maxe)
        call ZLACPY('F',N,2,Xe,N,Ye,N)
        do j=1,2
           call ZSCAL(N,dE(j), Ye(1,j), 1)
           rese(j)=DZNRM2(N,Ye(1,j),1) ! placeholder denominateur
        enddo
        call wzhcsrmm(UPLO,'N',N,N,2,ZONE,sa,isa,jsa,Xe,ZZERO,aux)
        linloops=itmaxb
        eps=epsb
        Xe=ZZERO !initial guess
        call zhbicgstab(0,UPLO,'N',sb,isb,jsb,N,2,aux,Xe,res,linloops,eps,comb,infoloc)
        do j=1,2
           call ZAXPY(N,-ZONE,Xe(1,j),1,Ye(1,j),1)
           rese(j)=DZNRM2(N,Ye(1,j),1)/rese(j)
           !if (j==2)
!           print *,it,j,dE(j),rese(j)
        enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do
  
  itmax=it

 
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
alpha(1:itmax)=dble(zalpha(1:itmax))
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)< alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo

!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
     call ZGEMM('N', 'N', N, itmax, itmax, ZONE, V(1,1), N, Xred, itmax, ZZERO, X, N)

  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
  deallocate(aux)
  deallocate(zalpha)


  deallocate(u)
  deallocate(w) 
  deallocate(V) 
  deallocate(H)
  deallocate(H_temp)
  deallocate(work)
  deallocate(Xred)
  deallocate(beta)

end subroutine zhgarnoldi



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!444444444444444444444 MPI matrix distribution and mat-vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPI


subroutine  dcsr_distribute_row(N,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,MY_COMM_WORLD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine distribute the global matrix by rows  among all MPI processors- distribution performed using values of startj and endj 
  ! N           (input) size of the global matrix A
  ! sa,isa,jsa  (input) CSR format for the matrix A (DOUBLE PRECISION)
  !            commun to all processors
  ! matAj  (output) TYPE(DCSR)-
  !        local by rows
  ! startj,endj (output) INTEGER- local values with the start and end of each row block
  ! Nlocal (output) INTEGER- local size of row block (endj-startj-1)
  ! Nsize (output) INTEGER(:) - size nbprocs
  !                global - store the size of each row blocks
  ! MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================
  implicit none
  include 'mpif.h'  
  integer :: startj,endj,MY_COMM_WORLD,N,Nlocal
  double precision,dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),pointer :: sa
     integer,dimension(:),pointer :: isa,jsa     
  end type dcsr
  type(dcsr) :: matAj
  integer :: code,rank,nbprocs,Remain,Nj
  integer,dimension(*) :: Nsize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)

!!!! find startj and endj
  Nj = N / nbprocs
  Remain = N - Nj*nbprocs
  startj = rank*Nj + 1

  !! put remain at the end to simplify
  if ((Remain > 0).and.(rank==nbprocs-1)) Nj = Nj+Remain

  endj = Nj + startj - 1

!!!! global size list
  Nsize(1:nbprocs)=0
  Nsize(rank+1)=endj-startj+1
  call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 


!!!!!!! set up the local matrices

  matAj%n=endj-startj+1
  matAj%nnz= isa(endj+1) - isa(startj)  

!!! distribution      
  allocate(matAj%sa(1:matAj%nnz))
  allocate(matAj%jsa(1:matAj%nnz))
  allocate(matAj%isa(1:matAj%n+1))

  matAj%isa(1:matAj%n+1)=isa(startj:endj+1)-isa(startj)+1
  matAj%jsa(1:matAj%nnz)=jsa(isa(startj):isa(endj+1)-1)
  matAj%sa(1:matAj%nnz)=sa(isa(startj):isa(endj+1)-1)

  Nlocal=endj-startj+1 ! NLocal


end subroutine dcsr_distribute_row



subroutine  zcsr_distribute_row(N,sa,isa,jsa,matAj,startj,endj,Nlocal,Nsize,MY_COMM_WORLD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine distribute the global matrix by rows  among all MPI processors- distribution performed using values of startj and endj 
  ! N           (input) size of the global matrix A
  ! sa,isa,jsa  (input) CSR format for the matrix A (COMPLEX DOUBLE PRECISION)
  !            commun to all processors
  ! matAj  (output) TYPE(ZCSR)-
  !        local by rows
  ! startj,endj (output) INTEGER- local values with the start and end of each row block
  ! Nlocal (output) INTEGER- local size of row block (endj-startj-1)
  ! Nsize (output) INTEGER(:) - size nbprocs
  !                global - store the size of each row blocks
  ! MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================
  implicit none
  include 'mpif.h'  
  integer :: startj,endj,MY_COMM_WORLD,N,Nlocal
  complex(kind=kind(1.0d0)),dimension(*) :: sa
  integer,dimension(*) :: isa,jsa
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),pointer :: sa
     integer,dimension(:),pointer :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj
  integer :: code,rank,nbprocs,Remain,Nj
  integer,dimension(*) :: Nsize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)

!!!! find startj and endj
  Nj = N / nbprocs
  Remain = N - Nj*nbprocs
  startj = rank*Nj + 1

  !! put remain at the end to simplify
  if ((Remain > 0).and.(rank==nbprocs-1)) Nj = Nj+Remain

  endj = Nj + startj - 1

!!!! global size list
  Nsize(1:nbprocs)=0
  Nsize(rank+1)=endj-startj+1
  call MPI_ALLREDUCE(MPI_IN_PLACE,Nsize,nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 


!!!!!!! set up the local matrices

  matAj%n=endj-startj+1
  matAj%nnz= isa(endj+1) - isa(startj)  

!!! distribution      
  allocate(matAj%sa(1:matAj%nnz))
  allocate(matAj%jsa(1:matAj%nnz))
  allocate(matAj%isa(1:matAj%n+1))

  matAj%isa(1:matAj%n+1)=isa(startj:endj+1)-isa(startj)+1
  matAj%jsa(1:matAj%nnz)=jsa(isa(startj):isa(endj+1)-1)
  matAj%sa(1:matAj%nnz)=sa(isa(startj):isa(endj+1)-1)

  Nlocal=endj-startj+1 ! NLocal


end subroutine zcsr_distribute_row





subroutine dget_local_block_csr(Nsize,map,matAj,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  !-------------------------------------
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine distribute the matrix by row/column csr block among MPI processors- as a condition the matrix must first be distributed by rows. The output is used by (d)lbcsrmm to perform the mpi mat-vec.
  !   This is the REAL DOUBLE PRECISION VERSION
  !
  !   Example using 3 MPI:

  !    |_mpi1_|              A11 A12 A13
  !  A=|_mpi2_| ==>          A21 A22 A23
  !    |_mpi3_|              A31 A32 A33
  !
  !  The block Aij is stored in both process i and j (cross configuration). A given processor stores up to 2*N-1 blocks.
  !
  !  Nsize (input) INTEGER(:) - size nbprocs
  !        global - store the size of each row blocks
  !  map (output) INTEGER(:,:) - size nbprocs*nbprocs
  !        global - store the nnz of each submatrix blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAj  (input) TYPE(DCSR)-
  !        local by row- see explanations above
  !
  ! matAjb (output) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
  !        local- see explanations above
  !
  ! MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  ! nbprocs        (input) INTEGER
  !                  #MPI
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================  
  integer :: MY_COMM_WORLD
  integer,dimension(*) :: Nsize  
  integer,dimension(nbprocs,*) :: map
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: matAj
  type(dcsr),dimension(nbprocs,*) :: matAjb
!!!
  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer :: p,rank,i,j,k,idummy,sj,ej,nbprocs,code,tag,request
!!!
  tag=100

!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  !call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)
!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Extract row blocks (use for normal mat-vec)!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  p=rank+1
!!! initialization
  do i=1,nbprocs
     matAjb(p,i)%nnz=0
     matAjb(p,i)%n=matAj%n
  enddo




!!! get nnz
  do i=1,matAj%n
     do k=matAj%isa(i),matAj%isa(i+1)-1
        idummy=Nsize(1)
        j=1
        do while (matAj%jsa(k)>idummy)
           j=j+1
           idummy=idummy+Nsize(j)
        end do
        matAjb(p,j)%nnz=matAjb(p,j)%nnz+1
     end do
  enddo

!!! get the submatrices along the row
  sj=1
  ej=Nsize(1)
  do i=1,nbprocs
     if (matAjb(p,i)%nnz/=0) then
        allocate(matAjb(p,i)%sa(matAjb(p,i)%nnz))
        allocate(matAjb(p,i)%jsa(matAjb(p,i)%nnz))
        allocate(matAjb(p,i)%isa(matAjb(p,i)%n+1))
        call dsubmat(matAj%n,1,matAj%n,sj,ej,matAj%sa,matAj%jsa,matAj%isa,matAjb(p,i)%sa,matAjb(p,i)%jsa,matAjb(p,i)%isa)
     end if
     sj=sj+Nsize(i)
     if (i/=nbprocs) ej=ej+Nsize(i+1)
  enddo

!!! generate the global nnz matrix
  map(1:nbprocs,1:nbprocs)=0!.false.
  do i=1,nbprocs
     map(p,i)=matAjb(p,i)%nnz!.true.
  enddo
  !call MPI_ALLREDUCE(MPI_IN_PLACE,map,nbprocs*nbprocs,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,code) 
  call MPI_ALLREDUCE(MPI_IN_PLACE,map,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 


  !print *,'map',map(1:nbprocs*nbprocs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Add column blocks (use for transpose mat-vec)!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! memory allocations
  do i=1,nbprocs ! columns
     if ((i/=p).and.(map(i,p)/=0)) then
        allocate(matAjb(i,p)%sa(map(i,p)))
        allocate(matAjb(i,p)%jsa(map(i,p)))
        allocate(matAjb(i,p)%isa(Nsize(i)+1))
     end if
  end do

  do i=1,nbprocs
     if ((i/=p).and.(map(p,i)/=0)) then ! send
        !       print *,'send',p,i,map((p-1)*nbprocs+i)
        call MPI_ISEND(matAjb(p,i)%sa,map(p,i),MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD,request,code)
        call MPI_ISEND(matAjb(p,i)%jsa,map(p,i),MPI_INTEGER,i-1,tag+1, MY_COMM_WORLD,request,code)
        call MPI_ISEND(matAjb(p,i)%isa,Nsize(p)+1,MPI_INTEGER,i-1,tag+2, MY_COMM_WORLD,request,code)
     endif
  end do

  do i=1,nbprocs
     if ((i/=p).and.(map(i,p)/=0)) then !! receive
        !         print *,'send',p,i,map((i-1)*nbprocs+p)
        call MPI_RECV(matAjb(i,p)%sa,map(i,p),MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,status,code)
        call MPI_RECV(matAjb(i,p)%jsa,map(i,p),MPI_INTEGER,i-1,tag+1, MY_COMM_WORLD ,status,code)
        call MPI_RECV(matAjb(i,p)%isa,Nsize(i)+1,MPI_INTEGER,i-1,tag+2, MY_COMM_WORLD ,status,code)       
     end if
  end do
end subroutine dget_local_block_csr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine zget_local_block_csr(Nsize,map,matAj,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  !-------------------------------------
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine distribute the matrix by row/column csr block among MPI processors- as a condition the matrix must first be distributed by rows. The output is used by (d)lbcsrmm to perform the mpi mat-vec.
  ! This is the COMPLEX DOUBLE PRECISION VERSION
  !
  !   Example using 3 MPI:

  !    |_mpi1_|              A11 A12 A13
  !  A=|_mpi2_| ==>          A21 A22 A23
  !    |_mpi3_|              A31 A32 A33
  !
  !  The block Aij is stored in both process i and j (cross configuration). A given processor stores up to 2*N-1 blocks.
  !
  !  Nsize (input) INTEGER(:) - size nbprocs
  !        global - store the size of each row blocks
  !  map (output) INTEGER(:,:) - size nbprocs*nbprocs
  !        global - store the nnz of each submatrix blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAj  (input) TYPE(ZCSR)-
  !        local by row- see explanations above
  !
  ! matAjb (output) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- see explanations above
  !
  ! MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  ! nbprocs        (input) INTEGER
  !                  #MPI
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================  
  integer :: MY_COMM_WORLD
  integer,dimension(*) :: Nsize  
  integer,dimension(nbprocs,*) :: map
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj
  type(zcsr),dimension(nbprocs,*) :: matAjb
!!!
  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer :: p,rank,i,j,k,idummy,sj,ej,nbprocs,code,tag,request
!!!
  tag=100

!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  !call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)
!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Extract row blocks (use for normal mat-vec)!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  p=rank+1
!!! initialization
  do i=1,nbprocs
     matAjb(p,i)%nnz=0
     matAjb(p,i)%n=matAj%n
  enddo




!!! get nnz
  do i=1,matAj%n
     do k=matAj%isa(i),matAj%isa(i+1)-1
        idummy=Nsize(1)
        j=1
        do while (matAj%jsa(k)>idummy)
           j=j+1
           idummy=idummy+Nsize(j)
        end do
        matAjb(p,j)%nnz=matAjb(p,j)%nnz+1
     end do
  enddo

!!! get the submatrices along the row
  sj=1
  ej=Nsize(1)
  do i=1,nbprocs
     if (matAjb(p,i)%nnz/=0) then
        allocate(matAjb(p,i)%sa(matAjb(p,i)%nnz))
        allocate(matAjb(p,i)%jsa(matAjb(p,i)%nnz))
        allocate(matAjb(p,i)%isa(matAjb(p,i)%n+1))
        call zsubmat(matAj%n,1,matAj%n,sj,ej,matAj%sa,matAj%jsa,matAj%isa,matAjb(p,i)%sa,matAjb(p,i)%jsa,matAjb(p,i)%isa)
     end if
     sj=sj+Nsize(i)
     if (i/=nbprocs) ej=ej+Nsize(i+1)
  enddo

!!! generate the global nnz matrix
  map(1:nbprocs,1:nbprocs)=0!.false.
  do i=1,nbprocs
     map(p,i)=matAjb(p,i)%nnz!.true.
  enddo
  !call MPI_ALLREDUCE(MPI_IN_PLACE,map,nbprocs*nbprocs,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,code) 
  call MPI_ALLREDUCE(MPI_IN_PLACE,map,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 


  !print *,'map',map(1:nbprocs*nbprocs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Add column blocks (use for transpose mat-vec)!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! memory allocations
  do i=1,nbprocs ! columns
     if ((i/=p).and.(map(i,p)/=0)) then
        allocate(matAjb(i,p)%sa(map(i,p)))
        allocate(matAjb(i,p)%jsa(map(i,p)))
        allocate(matAjb(i,p)%isa(Nsize(i)+1))
     end if
  end do

  do i=1,nbprocs
     if ((i/=p).and.(map(p,i)/=0)) then ! send
        !       print *,'send',p,i,map((p-1)*nbprocs+i)
        call MPI_ISEND(matAjb(p,i)%sa,map(p,i),MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD,request,code)
        call MPI_ISEND(matAjb(p,i)%jsa,map(p,i),MPI_INTEGER,i-1,tag+1, MY_COMM_WORLD,request,code)
        call MPI_ISEND(matAjb(p,i)%isa,Nsize(p)+1,MPI_INTEGER,i-1,tag+2, MY_COMM_WORLD,request,code)
     endif
  end do

  do i=1,nbprocs
     if ((i/=p).and.(map(i,p)/=0)) then !! receive
        !         print *,'send',p,i,map((i-1)*nbprocs+p)
        call MPI_RECV(matAjb(i,p)%sa,map(i,p),MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
        call MPI_RECV(matAjb(i,p)%jsa,map(i,p),MPI_INTEGER,i-1,tag+1, MY_COMM_WORLD ,status,code)
        call MPI_RECV(matAjb(i,p)%isa,Nsize(i)+1,MPI_INTEGER,i-1,tag+2, MY_COMM_WORLD ,status,code)       
     end if
  end do
end subroutine zget_local_block_csr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine cget_local_block_csr(Nsize,map,matAj,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  !-------------------------------------
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine distribute the matrix by row/column csr block among MPI processors- as a condition the matrix must first be distributed by rows. The output is used by (d)lbcsrmm to perform the mpi mat-vec.
  ! This is the COMPLEX SINGLE PRECISION VERSION
  !
  !   Example using 3 MPI:

  !    |_mpi1_|              A11 A12 A13
  !  A=|_mpi2_| ==>          A21 A22 A23
  !    |_mpi3_|              A31 A32 A33
  !
  !  The block Aij is stored in both process i and j (cross configuration). A given processor stores up to 2*N-1 blocks.
  !
  !  Nsize (input) INTEGER(:) - size nbprocs
  !        global - store the size of each row blocks
  !  map (output) INTEGER(:,:) - size nbprocs*nbprocs
  !        global - store the nnz of each submatrix blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAj  (input) TYPE(CCSR)-
  !        local by row- see explanations above
  !
  ! matAjb (output) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- see explanations above
  !
  ! MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  ! nbprocs        (input) INTEGER
  !                  #MPI
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================  
  integer :: MY_COMM_WORLD
  integer,dimension(*) :: Nsize  
  integer,dimension(nbprocs,*) :: map
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr) :: matAj
  type(ccsr),dimension(nbprocs,*) :: matAjb
!!!
  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer :: p,rank,i,j,k,idummy,sj,ej,nbprocs,code,tag,request
!!!
  tag=100

!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  !call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)
!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Extract row blocks (use for normal mat-vec)!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  p=rank+1
!!! initialization
  do i=1,nbprocs
     matAjb(p,i)%nnz=0
     matAjb(p,i)%n=matAj%n
  enddo


!!! get nnz
  do i=1,matAj%n
     do k=matAj%isa(i),matAj%isa(i+1)-1
        idummy=Nsize(1)
        j=1
        do while (matAj%jsa(k)>idummy)
           j=j+1
           idummy=idummy+Nsize(j)
        end do
        matAjb(p,j)%nnz=matAjb(p,j)%nnz+1
     end do
  enddo

!!! get the submatrices along the row
  sj=1
  ej=Nsize(1)
  do i=1,nbprocs
     if (matAjb(p,i)%nnz/=0) then
        allocate(matAjb(p,i)%sa(matAjb(p,i)%nnz))
        allocate(matAjb(p,i)%jsa(matAjb(p,i)%nnz))
        allocate(matAjb(p,i)%isa(matAjb(p,i)%n+1))
        call csubmat(matAj%n,1,matAj%n,sj,ej,matAj%sa,matAj%jsa,matAj%isa,matAjb(p,i)%sa,matAjb(p,i)%jsa,matAjb(p,i)%isa)
     end if
     sj=sj+Nsize(i)
     if (i/=nbprocs) ej=ej+Nsize(i+1)
  enddo

!!! generate the global nnz matrix
  map(1:nbprocs,1:nbprocs)=0!.false.
  do i=1,nbprocs
     map(p,i)=matAjb(p,i)%nnz!.true.
  enddo
  !call MPI_ALLREDUCE(MPI_IN_PLACE,map,nbprocs*nbprocs,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,code) 
  call MPI_ALLREDUCE(MPI_IN_PLACE,map,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 

  !print *,'map',map(1:nbprocs*nbprocs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Add column blocks (use for transpose mat-vec)!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! memory allocations
  do i=1,nbprocs ! columns
     if ((i/=p).and.(map(i,p)/=0)) then
        allocate(matAjb(i,p)%sa(map(i,p)))
        allocate(matAjb(i,p)%jsa(map(i,p)))
        allocate(matAjb(i,p)%isa(Nsize(i)+1))
     end if
  end do

  do i=1,nbprocs
     if ((i/=p).and.(map(p,i)/=0)) then ! send
        !       print *,'send',p,i,map((p-1)*nbprocs+i)
        call MPI_ISEND(matAjb(p,i)%sa,map(p,i),MPI_COMPLEX,i-1,tag, MY_COMM_WORLD,request,code)
        call MPI_ISEND(matAjb(p,i)%jsa,map(p,i),MPI_INTEGER,i-1,tag+1, MY_COMM_WORLD,request,code)
        call MPI_ISEND(matAjb(p,i)%isa,Nsize(p)+1,MPI_INTEGER,i-1,tag+2, MY_COMM_WORLD,request,code)
     endif
  end do

  do i=1,nbprocs
     if ((i/=p).and.(map(i,p)/=0)) then !! receive
        !         print *,'send',p,i,map((i-1)*nbprocs+p)
        call MPI_RECV(matAjb(i,p)%sa,map(i,p),MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
        call MPI_RECV(matAjb(i,p)%jsa,map(i,p),MPI_INTEGER,i-1,tag+1, MY_COMM_WORLD ,status,code)
        call MPI_RECV(matAjb(i,p)%isa,Nsize(i)+1,MPI_INTEGER,i-1,tag+2, MY_COMM_WORLD ,status,code)       
     end if
  end do
end subroutine cget_local_block_csr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









subroutine dlbcsr_transpose(Nsize,mapA,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the in-place transpose of a matrix distributed using the local block csr format.
  !  
  ! This is the real double precision version  
  !
  !   Example using 3 MPI:

  !                  A11 A12 A13        A11=A11^T A12=A21^T A13=A31^T
  !                  A21 A22 A23    ==> A21=A12^T A22=A22^T A23=A32^T
  !                  A31 A32 A33        A31=A13^T A32=A23^T A33=A33^T
  !
  ! Nsize     (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input/output) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: transpose matrix
  !
  ! NEW_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: temp
  type(dcsr),dimension(nbprocs,*) :: matAjb
  integer :: i,j,n,m,nnz1,nnz2,rank,p,code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)

  p=rank+1

  do j=1,nbprocs

     if (mapA(j,p)>0) then !! non-zero block to transpose
        n=Nsize(j)
        m=Nsize(p)

        nnz1=mapA(j,p)       
        nnz2=mapA(p,j)
        !          if (p==2) print *,'heree',j,n,m,nnz1,nnz2
        if (nnz2>0) then !! temporary copy of transpose of matAjb(p,j)
           !     m=matAjb(p,j)%n
           allocate(temp%sa(nnz1))
           allocate(temp%jsa(nnz1))
           allocate(temp%isa(n+1))
           call dcsr_transpose(m,n,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,temp%sa,temp%isa,temp%jsa)

           if (nnz1==nnz2) then ! same memory space 
              ! transpose (if not diagonal block)       
              if (j/=p)  call dcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)
              matAjb(p,j)%nnz=nnz1

           else ! different memory space
              ! re-allocate
              deallocate(matAjb(p,j)%sa)
              deallocate(matAjb(p,j)%jsa)
              allocate(matAjb(p,j)%sa(1:nnz1))
              allocate(matAjb(p,j)%jsa(1:nnz1))
              matAjb(p,j)%nnz=nnz1
              ! transpose
              call dcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)
              ! re-allocate
              deallocate(matAjb(j,p)%sa)
              deallocate(matAjb(j,p)%jsa)
              allocate(matAjb(j,p)%sa(1:nnz2))
              allocate(matAjb(j,p)%jsa(1:nnz2))
              matAjb(j,p)%nnz=nnz2
           end if

           ! paste
           matAjb(j,p)%sa(1:nnz2)=temp%sa(1:nnz2)
           matAjb(j,p)%jsa(1:nnz2)=temp%jsa(1:nnz2)
           matAjb(j,p)%isa(1:n+1)=temp%isa(1:n+1)
           matAjb(j,p)%nnz=nnz2

           ! deallocate temp
           deallocate(temp%sa)
           deallocate(temp%jsa)
           deallocate(temp%isa)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        else !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
           matAjb(p,j)%n=m
           matAjb(p,j)%nnz=nnz1
           allocate(matAjb(p,j)%sa(nnz1))
           allocate(matAjb(p,j)%jsa(nnz1))
           allocate(matAjb(p,j)%isa(m+1))
           matAjb(p,j)%nnz=nnz1
           ! transpose
           call dcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)
           ! empty the old space
           matAjb(j,p)%n=0
           matAjb(j,p)%nnz=0
           deallocate(matAjb(j,p)%sa)
           deallocate(matAjb(j,p)%jsa)
           deallocate(matAjb(j,p)%isa)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else  !! (j,p) is a zero block so transpose (p,j) into (j,p) if not zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

        nnz2=mapA(p,j)
        !      if (p==2) print *,'heree',j,nnz2
        if (nnz2>0) then !! perform transpose in empty space
           n=Nsize(p)
           m=Nsize(j)
           allocate(matAjb(j,p)%sa(nnz2))
           allocate(matAjb(j,p)%jsa(nnz2))
           allocate(matAjb(j,p)%isa(m+1))
           matAjb(j,p)%nnz=nnz2
           call dcsr_transpose(n,m,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa)
           ! empty the old space
           matAjb(p,j)%n=0
           matAjb(p,j)%nnz=0
           deallocate(matAjb(p,j)%sa)
           deallocate(matAjb(p,j)%jsa)
           deallocate(matAjb(p,j)%isa)
        end if

     end if
  end do

!!! generate the global nnz matrix
  mapA(1:nbprocs,1:nbprocs)=0
  do i=1,nbprocs
     mapA(i,p)=matAjb(i,p)%nnz
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,mapA,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code)      

end subroutine dlbcsr_transpose



subroutine zlbcsr_transpose(Nsize,mapA,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the in-place transpose of a matrix distributed using the local block csr format.
  !  
  ! This is the complex double precision version  
  !
  !   Example using 3 MPI:

  !                  A11 A12 A13        A11=A11^T A12=A21^T A13=A31^T
  !                  A21 A22 A23    ==> A21=A12^T A22=A22^T A23=A32^T
  !                  A31 A32 A33        A31=A13^T A32=A23^T A33=A33^T
  !
  ! Nsize     (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input/output) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: transpose matrix
  !
  ! NEW_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: temp
  type(zcsr),dimension(nbprocs,*) :: matAjb
  integer :: i,j,n,m,nnz1,nnz2,rank,p,code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)

  p=rank+1

  do j=1,nbprocs

     if (mapA(j,p)>0) then !! non-zero block to transpose
        n=Nsize(j)
        m=Nsize(p)

        nnz1=mapA(j,p)       
        nnz2=mapA(p,j)
        !          if (p==2) print *,'heree',j,n,m,nnz1,nnz2
        if (nnz2>0) then !! temporary copy of transpose of matAjb(p,j)
           !     m=matAjb(p,j)%n
           allocate(temp%sa(nnz1))
           allocate(temp%jsa(nnz1))
           allocate(temp%isa(n+1))
           call zcsr_transpose(m,n,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,temp%sa,temp%isa,temp%jsa)

           if (nnz1==nnz2) then ! same memory space 
              ! transpose (if not diagonal block)       
              if (j/=p)  call zcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)
              matAjb(p,j)%nnz=nnz1

           else ! different memory space
              ! re-allocate
              deallocate(matAjb(p,j)%sa)
              deallocate(matAjb(p,j)%jsa)
              allocate(matAjb(p,j)%sa(1:nnz1))
              allocate(matAjb(p,j)%jsa(1:nnz1))
              matAjb(p,j)%nnz=nnz1
              ! transpose
              call zcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)
              ! re-allocate
              deallocate(matAjb(j,p)%sa)
              deallocate(matAjb(j,p)%jsa)
              allocate(matAjb(j,p)%sa(1:nnz2))
              allocate(matAjb(j,p)%jsa(1:nnz2))
              matAjb(j,p)%nnz=nnz2
           end if

           ! paste
           matAjb(j,p)%sa(1:nnz2)=temp%sa(1:nnz2)
           matAjb(j,p)%jsa(1:nnz2)=temp%jsa(1:nnz2)
           matAjb(j,p)%isa(1:n+1)=temp%isa(1:n+1)
           matAjb(j,p)%nnz=nnz2

           ! deallocate temp
           deallocate(temp%sa)
           deallocate(temp%jsa)
           deallocate(temp%isa)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        else !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
           matAjb(p,j)%n=m
           matAjb(p,j)%nnz=nnz1
           allocate(matAjb(p,j)%sa(nnz1))
           allocate(matAjb(p,j)%jsa(nnz1))
           allocate(matAjb(p,j)%isa(m+1))
           matAjb(p,j)%nnz=nnz1
           ! transpose
           call zcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)
           ! empty the old space
           matAjb(j,p)%n=0
           matAjb(j,p)%nnz=0
           deallocate(matAjb(j,p)%sa)
           deallocate(matAjb(j,p)%jsa)
           deallocate(matAjb(j,p)%isa)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else  !! (j,p) is a zero block so transpose (p,j) into (j,p) if not zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

        nnz2=mapA(p,j)
        !      if (p==2) print *,'heree',j,nnz2
        if (nnz2>0) then !! perform transpose in empty space
           n=Nsize(p)
           m=Nsize(j)
           allocate(matAjb(j,p)%sa(nnz2))
           allocate(matAjb(j,p)%jsa(nnz2))
           allocate(matAjb(j,p)%isa(m+1))
           matAjb(j,p)%nnz=nnz2
           call zcsr_transpose(n,m,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa)
           ! empty the old space
           matAjb(p,j)%n=0
           matAjb(p,j)%nnz=0
           deallocate(matAjb(p,j)%sa)
           deallocate(matAjb(p,j)%jsa)
           deallocate(matAjb(p,j)%isa)
        end if

     end if
  end do

!!! generate the global nnz matrix
  mapA(1:nbprocs,1:nbprocs)=0
  do i=1,nbprocs
     mapA(i,p)=matAjb(i,p)%nnz
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,mapA,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code)      

end subroutine zlbcsr_transpose



subroutine csaddlbcsr(Nsize,alpha,mapA,matAjb,beta,mapB,matBjb,mapC,matCjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !
  !  Addition of two  sparse local block CSR matrices A and B and return a new local block CSR matrix C 
  !  C=alpha*A+beta*B
  !  Remark: sparity pattern of A and B must already be included in C
  !
  !  C is complex single precision, A and B are real single precision, (alpha,beta are COMPLEX)
  !  
   !  Arguments
  !  =========
  !  Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha,beta  (input) COMPLEX SINGLE PRECISION
  !  mapA, mapB, mapC (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb,matBjb (input) Matrix of TYPE(SCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matCjb (input/output) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix C<=alpha*A+beta*B
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
 implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA,mapB,mapC
 type scsr
     integer ::n,m,nnz
     real,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type scsr
  type(scsr),dimension(nbprocs,*) :: matAjb,matBjb
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr),dimension(nbprocs,*) :: matCjb
  complex :: alpha,beta
  !--------------------
  integer :: rank,p,code

  integer :: j,nnz
  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
 

  !column scan
  do j=1,nbprocs
      if (mapC(p,j)/=0) then ! means that either A or B (or both) /=0
        if ((mapA(p,j)/=0).and.(mapB(p,j)/=0)) then
           call csaddcsr(Nsize(p),Nsize(j),3,alpha,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,beta,matBjb(p,j)%sa,matBjb(p,j)%isa,matBjb(p,j)%jsa,matCjb(p,j)%sa,matCjb(p,j)%isa,matCjb(p,j)%jsa)
        elseif (mapA(p,j)/=0) then ! so B block is equal to zero
           nnz=matAjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=alpha*matAjb(p,j)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=beta*matBjb(p,j)%sa(1:nnz)
        end if
      end if
   end do
   
   !row scan minus diagonal
  do j=1,nbprocs
      if ((mapC(j,p)/=0).and.(p/=j)) then ! means that either A or B (or both) /=0
        if ((mapA(j,p)/=0).and.(mapB(j,p)/=0)) then
           call csaddcsr(Nsize(j),Nsize(p),3,alpha,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,beta,matBjb(j,p)%sa,matBjb(j,p)%isa,matBjb(j,p)%jsa,matCjb(j,p)%sa,matCjb(j,p)%isa,matCjb(j,p)%jsa)
        elseif (mapA(j,p)/=0) then ! so B block is equal to zero
           nnz=matAjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=alpha*matAjb(j,p)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=beta*matBjb(j,p)%sa(1:nnz)
        end if
      end if
  end do
  
end subroutine csaddlbcsr



subroutine zdaddlbcsr(Nsize,alpha,mapA,matAjb,beta,mapB,matBjb,mapC,matCjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !
  !  Addition of two  sparse local block CSR matrices A and B and return a new local block CSR matrix C 
  !  C=alpha*A+beta*B
  !  Remark: sparity pattern of A and B must already be included in C
  !
  !  C is complex double precision, A and B are real double precision, (alpha,beta are COMPLEX)
  !  
   !  Arguments
  !  =========
  !  Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha,beta  (input) COMPLEX DOUBLE PRECISION
  !  mapA, mapB, mapC (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb,matBjb (input) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matCjb (input/output) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix C<=alpha*A+beta*B
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
 implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA,mapB,mapC
 type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr),dimension(nbprocs,*) :: matAjb,matBjb
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matCjb
  complex(kind=kind(1.0d0)) :: alpha,beta
  !--------------------
  integer :: rank,p,code,opt

  integer :: j,nnz
  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
  opt=3

  !column scan
  do j=1,nbprocs
      if (mapC(p,j)/=0) then ! means that either A or B (or both) /=0
        if ((mapA(p,j)/=0).and.(mapB(p,j)/=0)) then
          call zdaddcsr(Nsize(p),Nsize(j),opt,alpha,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,beta,matBjb(p,j)%sa,matBjb(p,j)%isa,matBjb(p,j)%jsa,matCjb(p,j)%sa,matCjb(p,j)%isa,matCjb(p,j)%jsa)
        elseif (mapA(p,j)/=0) then ! so B block is equal to zero
           nnz=matAjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=alpha*matAjb(p,j)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=beta*matBjb(p,j)%sa(1:nnz)
        end if
      end if
   end do
   
   !row scan minus diagonal
  do j=1,nbprocs
      if ((mapC(j,p)/=0).and.(p/=j)) then ! means that either A or B (or both) /=0
        if ((mapA(j,p)/=0).and.(mapB(j,p)/=0)) then
           call zdaddcsr(Nsize(j),Nsize(p),opt,alpha,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,beta,matBjb(j,p)%sa,matBjb(j,p)%isa,matBjb(j,p)%jsa,matCjb(j,p)%sa,matCjb(j,p)%isa,matCjb(j,p)%jsa)
        elseif (mapA(j,p)/=0) then ! so B block is equal to zero
           nnz=matAjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=alpha*matAjb(j,p)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=beta*matBjb(j,p)%sa(1:nnz)
        end if
      end if
  end do
  


end subroutine zdaddlbcsr



subroutine caddlbcsr(Nsize,alpha,mapA,matAjb,beta,mapB,matBjb,mapC,matCjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !
  !  Addition of two  sparse local block CSR matrices A and B and return a new local block CSR matrix C 
  !  C=alpha*A+beta*B
  !  Remark: sparity pattern of A and B must already be included in C
  !
  !  A, B, C are complex single precision, (alpha,beta are also COMPLEX)
  !  
   !  Arguments
  !  =========
  !  Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha,beta  (input) COMPLEX SINGLE PRECISION
  !  mapA, mapB, mapC (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb,matBjb (input) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matCjb (input/output) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix C<=alpha*A+beta*B
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
 implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA,mapB,mapC
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr),dimension(nbprocs,*) :: matAjb,matBjb,matCjb
  complex :: alpha,beta
  !--------------------
  integer :: rank,p,code

  integer :: j,nnz
  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
 

  !column scan
  do j=1,nbprocs
      if (mapC(p,j)/=0) then ! means that either A or B (or both) /=0
        if ((mapA(p,j)/=0).and.(mapB(p,j)/=0)) then
           call caddcsr(Nsize(p),Nsize(j),3,alpha,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,beta,matBjb(p,j)%sa,matBjb(p,j)%isa,matBjb(p,j)%jsa,matCjb(p,j)%sa,matCjb(p,j)%isa,matCjb(p,j)%jsa)
        elseif (mapA(p,j)/=0) then ! so B block is equal to zero
           nnz=matAjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=alpha*matAjb(p,j)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=beta*matBjb(p,j)%sa(1:nnz)
        end if
      end if
   end do
   
   !row scan minus diagonal
  do j=1,nbprocs
      if ((mapC(j,p)/=0).and.(p/=j)) then ! means that either A or B (or both) /=0
        if ((mapA(j,p)/=0).and.(mapB(j,p)/=0)) then
           call caddcsr(Nsize(j),Nsize(p),3,alpha,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,beta,matBjb(j,p)%sa,matBjb(j,p)%isa,matBjb(j,p)%jsa,matCjb(j,p)%sa,matCjb(j,p)%isa,matCjb(j,p)%jsa)
        elseif (mapA(j,p)/=0) then ! so B block is equal to zero
           nnz=matAjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=alpha*matAjb(j,p)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=beta*matBjb(j,p)%sa(1:nnz)
        end if
      end if
  end do
  
end subroutine caddlbcsr



subroutine zaddlbcsr(Nsize,alpha,mapA,matAjb,beta,mapB,matBjb,mapC,matCjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !
  !  Addition of two  sparse local block CSR matrices A and B and return a new local block CSR matrix C 
  !  C=alpha*A+beta*B
  !  Remark: sparity pattern of A and B must already be included in C
  !
  !  A,B,C are complex double precision (alpha,beta are also COMPLEX)
  !  
   !  Arguments
  !  =========
  !  Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha,beta  (input) COMPLEX DOUBLE PRECISION
  !  mapA, mapB, mapC (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb,matBjb (input) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matCjb (input/output) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix C<=alpha*A+beta*B
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
 implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA,mapB,mapC
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matAjb,matBjb,matCjb
  complex(kind=kind(1.0d0)) :: alpha,beta
  !--------------------
  integer :: rank,p,code,opt

  integer :: j,nnz
  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
  opt=3

  !column scan
  do j=1,nbprocs
      if (mapC(p,j)/=0) then ! means that either A or B (or both) /=0
        if ((mapA(p,j)/=0).and.(mapB(p,j)/=0)) then
          call zaddcsr(Nsize(p),Nsize(j),opt,alpha,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,beta,matBjb(p,j)%sa,matBjb(p,j)%isa,matBjb(p,j)%jsa,matCjb(p,j)%sa,matCjb(p,j)%isa,matCjb(p,j)%jsa)
        elseif (mapA(p,j)/=0) then ! so B block is equal to zero
           nnz=matAjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=alpha*matAjb(p,j)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(p,j)%isa(Nsize(p)+1)-1
           matCjb(p,j)%sa(1:nnz)=beta*matBjb(p,j)%sa(1:nnz)
        end if
      end if
   end do
   
   !row scan minus diagonal
  do j=1,nbprocs
      if ((mapC(j,p)/=0).and.(p/=j)) then ! means that either A or B (or both) /=0
        if ((mapA(j,p)/=0).and.(mapB(j,p)/=0)) then
           call zaddcsr(Nsize(j),Nsize(p),opt,alpha,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,beta,matBjb(j,p)%sa,matBjb(j,p)%isa,matBjb(j,p)%jsa,matCjb(j,p)%sa,matCjb(j,p)%isa,matCjb(j,p)%jsa)
        elseif (mapA(j,p)/=0) then ! so B block is equal to zero
           nnz=matAjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=alpha*matAjb(j,p)%sa(1:nnz)
        elseif (mapB(p,j)/=0) then ! so A block is equal to zero
           nnz=matBjb(j,p)%isa(Nsize(j)+1)-1
           matCjb(j,p)%sa(1:nnz)=beta*matBjb(j,p)%sa(1:nnz)
        end if
      end if
  end do
  


end subroutine zaddlbcsr



subroutine zinc_addlbcsr(Nsize,alpha,mapA,matAjb,mapB,matBjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !  
  !  Increment matrix B by adding matrix A times alpha
  !  B=B+alpha*A 
  !
  !  A,B are complex double precision
  !  We consider local block CSR matrices and that the lbcsr format of B "includes" the lbcsr format of A (matrices should be compatible, i.e. B contains *at least* all the same non-zeros elements than A)
  ! 
  !  Arguments
  !  =========
  ! Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha  (input) COMPLEX DOUBLE PRECISION
  !  mapA, mapB (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matBjb (input/output) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix B<=alpha*A+B
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA,mapB
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matAjb,matBjb
  complex(kind=kind(1.0d0)) :: alpha
  !--------------------
  integer :: rank,p,code
  integer :: j

  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1

  !column scan
  do j=1,nbprocs
  if (mapA(p,j)/=0) then
  call zinc_addcsr(Nsize(p),Nsize(j),alpha,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,matBjb(p,j)%sa,matBjb(p,j)%isa,matBjb(p,j)%jsa)
  endif
  enddo   
 !row scan minus diagonal
  do j=1,nbprocs
  if ((mapA(j,p)/=0).and.(j/=p)) then
  call zinc_addcsr(Nsize(j),Nsize(p),alpha,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matBjb(j,p)%sa,matBjb(j,p)%isa,matBjb(j,p)%jsa)
  endif
  enddo   


end subroutine zinc_addlbcsr



subroutine cinc_addlbcsr(Nsize,alpha,mapA,matAjb,mapB,matBjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !  
  !  Increment matrix B by adding matrix A times alpha
  !  B=B+alpha*A 
  !
  !  A,B are complex single precision
  !  We consider local block CSR matrices and that the lbcsr format of B "includes" the lbcsr format of A (matrices should be compatible, i.e. B contains *at least* all the same non-zeros elements than A)
  ! 
  !  Arguments
  !  =========
  ! Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha  (input) COMPLEX SINGLE PRECISION
  !  mapA, mapB (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matBjb (input/output) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix B<=alpha*A+B
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA,mapB
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr),dimension(nbprocs,*) :: matAjb,matBjb
  complex :: alpha
  !--------------------
  integer :: rank,p,code
  integer :: j

  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1

  !column scan
  do j=1,nbprocs
  if (mapA(p,j)/=0) then
  call cinc_addcsr(Nsize(p),Nsize(j),alpha,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,matBjb(p,j)%sa,matBjb(p,j)%isa,matBjb(p,j)%jsa)
  endif
  enddo   
 !row scan minus diagonal
  do j=1,nbprocs
  if ((mapA(j,p)/=0).and.(j/=p)) then
  call cinc_addcsr(Nsize(j),Nsize(p),alpha,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matBjb(j,p)%sa,matBjb(j,p)%isa,matBjb(j,p)%jsa)
  endif
  enddo   


end subroutine cinc_addlbcsr




subroutine dlbcsr2s(Nsize,mapA,matAjb,matBjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !  
  !  Copy/convert a local block csr format matrix from
  !  real double precision to single precision
  !
  !  Arguments
  !  =========
  ! Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha  (input) COMPLEX DOUBLE PRECISION
  !  mapA   (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matBjb (output) Matrix of TYPE(SCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix B<=real(A)
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr),dimension(nbprocs,*) :: matAjb
 type scsr
     integer ::n,m,nnz
     real,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type scsr
  type(scsr),dimension(nbprocs,*) :: matBjb
  !--------------------
  integer :: rank,p,code
  integer :: j,nnz,infoloc

  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1

  !column scan
  do j=1,nbprocs
     if (mapA(p,j)/=0) then
        nnz=matAjb(p,j)%isa(Nsize(p)+1)-1
        allocate(matBjb(p,j)%sa(nnz))
        allocate(matBjb(p,j)%jsa(nnz))
        allocate(matBjb(p,j)%isa(Nsize(p)+1))
        call DLAG2S(nnz,1,matAjb(p,j)%sa,nnz,matBjb(p,j)%sa,nnz,infoloc)
        matBjb(p,j)%jsa(1:nnz)= matAjb(p,j)%jsa(1:nnz)
        matBjb(p,j)%isa(1:Nsize(p)+1)= matAjb(p,j)%isa(1:Nsize(p)+1)
        matBjb(p,j)%nnz=nnz
        matBjb(p,j)%n=Nsize(p)
  endif
  enddo   
 !row scan minus diagonal
  do j=1,nbprocs
     if ((mapA(j,p)/=0).and.(j/=p)) then
        nnz=matAjb(j,p)%isa(Nsize(j)+1)-1  
        allocate(matBjb(j,p)%sa(nnz))
        allocate(matBjb(j,p)%jsa(nnz))
        allocate(matBjb(j,p)%isa(Nsize(j)+1))
        call DLAG2S(nnz,1,matAjb(j,p)%sa,nnz,matBjb(j,p)%sa,nnz,infoloc)
        matBjb(j,p)%jsa(1:nnz)= matAjb(j,p)%jsa(1:nnz)
        matBjb(j,p)%isa(1:Nsize(j)+1)= matAjb(j,p)%isa(1:Nsize(j)+1)
        matBjb(j,p)%nnz=nnz
        matBjb(j,p)%n=Nsize(j)
  endif
  enddo   


end subroutine dlbcsr2s







subroutine zlbcsr2c(Nsize,mapA,matAjb,matBjb,MY_COMM_WORLD,nbprocs)
  !  Purpose
  !  =======
  !  
  !  Copy/convert a local block csr format matrix from
  !  complex double precision to single precision
  !
  !  Arguments
  !  =========
  ! Nsize   (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  alpha  (input) COMPLEX DOUBLE PRECISION
  !  mapA   (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !              |23  12  0 |
  !        mapA= |12  18  0 |
  !              |0    0  16|
  !
  !        it means (for matrix A) that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  ! matBjb (output) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !        On exit: modified matrix B<=cmplx(A)
  !
  ! MY_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  implicit none
  include 'mpif.h'
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matAjb
 type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr),dimension(nbprocs,*) :: matBjb
  !--------------------
  integer :: rank,p,code
  integer :: j,nnz,infoloc

  
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1

  
  !column scan
  do j=1,nbprocs
     if (mapA(p,j)/=0) then
        nnz=matAjb(p,j)%isa(Nsize(p)+1)-1
        allocate(matBjb(p,j)%sa(nnz))
        allocate(matBjb(p,j)%jsa(nnz))
        allocate(matBjb(p,j)%isa(Nsize(p)+1))
        call ZLAG2C(nnz,1,matAjb(p,j)%sa,nnz,matBjb(p,j)%sa,nnz,infoloc)
        matBjb(p,j)%jsa(1:nnz)= matAjb(p,j)%jsa(1:nnz)
        matBjb(p,j)%isa(1:Nsize(p)+1)= matAjb(p,j)%isa(1:Nsize(p)+1)
        matBjb(p,j)%nnz=nnz
        matBjb(p,j)%n=Nsize(p)
  endif
  enddo   
  !row scan minus diagonal
  !print *,p,mapA(1:nbprocs,1:nbprocs)
  do j=1,nbprocs
     if ((mapA(j,p)/=0).and.(j/=p)) then
        nnz=matAjb(j,p)%isa(Nsize(j)+1)-1  
        allocate(matBjb(j,p)%sa(nnz))
        allocate(matBjb(j,p)%jsa(nnz))
        allocate(matBjb(j,p)%isa(Nsize(j)+1))
        call ZLAG2C(nnz,1,matAjb(j,p)%sa,nnz,matBjb(j,p)%sa,nnz,infoloc)
        matBjb(j,p)%jsa(1:nnz)= matAjb(j,p)%jsa(1:nnz)
        matBjb(j,p)%isa(1:Nsize(j)+1)= matAjb(j,p)%isa(1:Nsize(j)+1)
        matBjb(j,p)%nnz=nnz
        matBjb(j,p)%n=Nsize(j)
  endif
  enddo   


end subroutine zlbcsr2c







subroutine dlbcsr_uplo_to_csr(Nsize,mapA,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose: This subroutine convert a lower or upper matrix distributed using the lower/upper local block csr format into
  !           a full local block csr format.
  ! This is the real double precision version
  !
  !   Example using 3 MPI:

  !                  A11  0  0        A11=A11^T  A12=A21^T  A13=A31^T
  !                  A21 A22 0     ==>A21        A22=A22^T  A23=A32^T
  !                  A31 A32 A33      A31          A32      A33=A33^T
  !
  ! Nsize     (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input/output) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
  !        local- block lower/upper csr format
  !        On exit: block full csr format
  !
  !
  ! NEW_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: temp
  type(dcsr),dimension(nbprocs,*) :: matAjb
  integer :: i,j,n,m,nnz1,nnz2,rank,p,code,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)

  p=rank+1

  do j=1,nbprocs !! scan along the rows
     if (mapA(j,p)>0) then !! non-zero block to transpose
        n=Nsize(j)
        m=Nsize(p) 
        nnz1=mapA(j,p)       
        !     nnz2=mapA(p,j)
        !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (p/=j) then
           matAjb(p,j)%n=m
           matAjb(p,j)%nnz=nnz1
           allocate(matAjb(p,j)%sa(nnz1))
           allocate(matAjb(p,j)%jsa(nnz1))
           allocate(matAjb(p,j)%isa(m+1))
           matAjb(p,j)%nnz=nnz1
           call dcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)      
        else ! p==j 
           allocate(temp%sa(2*nnz1)) ! a bit overestimated
           allocate(temp%jsa(2*nnz1))! a bit overestimated
           allocate(temp%isa(m+1))
           call dcsr_uplo_to_csr(n,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,temp%sa,temp%isa,temp%jsa)

           ! copy result
           deallocate(matAjb(j,p)%sa)
           deallocate(matAjb(j,p)%jsa)
           nnz1=temp%isa(n+1)-1
           allocate(matAjb(j,p)%sa(nnz1))
           allocate(matAjb(j,p)%jsa(nnz1))
           matAjb(p,j)%sa(1:nnz1)= temp%sa(1:nnz1)
           matAjb(p,j)%jsa(1:nnz1)=temp%jsa(1:nnz1)
           matAjb(p,j)%isa(1:n+1)=temp%isa(1:n+1)
           matAjb(p,j)%nnz=nnz1
           deallocate(temp%sa)
           deallocate(temp%jsa)
           deallocate(temp%isa)
        end if
     end if
  end do

  do j=1,nbprocs   !! scan along the columns minus diagonal
     if ((mapA(p,j)>0).and.(j/=p)) then !! non-zero block to transpose
        n=Nsize(p)
        m=Nsize(j) 
        nnz1=mapA(p,j)       
        !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           matAjb(j,p)%n=m
           matAjb(j,p)%nnz=nnz1
           allocate(matAjb(j,p)%sa(nnz1))
           allocate(matAjb(j,p)%jsa(nnz1))
           allocate(matAjb(j,p)%isa(m+1))
           matAjb(j,p)%nnz=nnz1
           call dcsr_transpose(n,m,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa)      
        end if
     end do

!!! generate the global nnz matrix
  mapA(1:nbprocs,1:nbprocs)=0
  do i=1,nbprocs
     mapA(p,i)=matAjb(p,i)%nnz
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,mapA,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 


end subroutine dlbcsr_uplo_to_csr



subroutine zhlbcsr_uplo_to_csr(Nsize,mapA,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose: This subroutine convert a lower or upper matrix distributed using the lower/upper local block csr format into
  !           a full local block csr format.
  ! This is the complex double precision version-
  ! where input matrices are complex Hermitian
  !
  !   Example using 3 MPI:

  !                  A11  0  0        A11=A11^H A12=A21^H  A13=A31^H
  !                  A21 A22 0    ==>    A21    A22=A22^H  A23=A32^H
  !                  A31 A32 A33         A31       A23     A33=A33^H
  !
  ! Nsize     (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input/output) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block lower/upper csr format
  !        On exit: block full csr format
  !
  ! NEW_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: temp
  type(zcsr),dimension(nbprocs,*) :: matAjb
  integer :: i,j,n,m,nnz1,nnz2,rank,p,code,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)

  p=rank+1

  do j=1,nbprocs   !! scan along the rows
     if (mapA(j,p)>0) then !! non-zero block to transpose
        n=Nsize(j)
        m=Nsize(p) 
        nnz1=mapA(j,p)       
        !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (p/=j) then
           matAjb(p,j)%n=m
           matAjb(p,j)%nnz=nnz1
           allocate(matAjb(p,j)%sa(nnz1))
           allocate(matAjb(p,j)%jsa(nnz1))
           allocate(matAjb(p,j)%isa(m+1))
           matAjb(p,j)%nnz=nnz1
           call zcsr_htranspose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)      
        else ! p==j 
           allocate(temp%sa(2*nnz1)) ! a bit overestimated
           allocate(temp%jsa(2*nnz1))! a bit overestimated
           allocate(temp%isa(m+1))
           call zhcsr_uplo_to_csr(n,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,temp%sa,temp%isa,temp%jsa)

           ! copy result
           deallocate(matAjb(j,p)%sa)
           deallocate(matAjb(j,p)%jsa)
           nnz1=temp%isa(n+1)-1
           allocate(matAjb(j,p)%sa(nnz1))
           allocate(matAjb(j,p)%jsa(nnz1))
           matAjb(p,j)%sa(1:nnz1)= temp%sa(1:nnz1)
           matAjb(p,j)%jsa(1:nnz1)=temp%jsa(1:nnz1)
           matAjb(p,j)%isa(1:n+1)=temp%isa(1:n+1)
           matAjb(p,j)%nnz=nnz1
           deallocate(temp%sa)
           deallocate(temp%jsa)
           deallocate(temp%isa)
        end if
     end if
  end do

 do j=1,nbprocs   !! scan along the columns minus diagonal
     if ((mapA(p,j)>0).and.(j/=p)) then !! non-zero block to transpose
        n=Nsize(p)
        m=Nsize(j) 
        nnz1=mapA(p,j)       
        !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           matAjb(j,p)%n=n
           matAjb(j,p)%nnz=nnz1
           allocate(matAjb(j,p)%sa(nnz1))
           allocate(matAjb(j,p)%jsa(nnz1))
           allocate(matAjb(j,p)%isa(m+1))
           matAjb(j,p)%nnz=nnz1
           call zcsr_htranspose(n,m,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa) 
     end if
  end do

!!! generate the global nnz matrix
  mapA(1:nbprocs,1:nbprocs)=0
  do i=1,nbprocs
     mapA(p,i)=matAjb(p,i)%nnz
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,mapA,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 


end subroutine zhlbcsr_uplo_to_csr




subroutine zlbcsr_uplo_to_csr(Nsize,mapA,matAjb,MY_COMM_WORLD,nbprocs)
  implicit none
  include 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose: This subroutine convert a lower or upper matrix distributed using the lower/upper local block csr format into
  !           a full local block csr format.
  ! This is the complex double precision version-
  ! where input matrices are complex symmetric
  !
  !   Example using 3 MPI:
  !                  A11  0  0        A11=A11^T A12=A21^T  A13=A31^T
  !                  A21 A22 0    ==>    A21    A22=A22^T  A23=A32^T
  !                  A31 A32 A33         A31       A23     A33=A33^T
  !
  ! Nsize     (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input/output) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block lower/upper csr format
  !        On exit: block full csr format
  !
  ! NEW_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: temp
  type(zcsr),dimension(nbprocs,*) :: matAjb
  integer :: i,j,n,m,nnz1,nnz2,rank,p,code,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)

  p=rank+1

  do j=1,nbprocs  !! scan along the rows
     if (mapA(j,p)>0) then !! non-zero block to transpose
        n=Nsize(j)
        m=Nsize(p) 
        nnz1=mapA(j,p)       
        !     nnz2=mapA(p,j)
        !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (p/=j) then
           matAjb(p,j)%n=m
           matAjb(p,j)%nnz=nnz1
           allocate(matAjb(p,j)%sa(nnz1))
           allocate(matAjb(p,j)%jsa(nnz1))
           allocate(matAjb(p,j)%isa(m+1))
           matAjb(p,j)%nnz=nnz1
           call zcsr_transpose(n,m,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa)      
        else ! p==j 
           allocate(temp%sa(2*nnz1)) ! a bit overestimated
           allocate(temp%jsa(2*nnz1))! a bit overestimated
           allocate(temp%isa(m+1))
           call zcsr_uplo_to_csr(n,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa,temp%sa,temp%isa,temp%jsa)

           ! copy result
           deallocate(matAjb(j,p)%sa)
           deallocate(matAjb(j,p)%jsa)
           nnz1=temp%isa(n+1)-1
           allocate(matAjb(j,p)%sa(nnz1))
           allocate(matAjb(j,p)%jsa(nnz1))
           matAjb(p,j)%sa(1:nnz1)= temp%sa(1:nnz1)
           matAjb(p,j)%jsa(1:nnz1)=temp%jsa(1:nnz1)
           matAjb(p,j)%isa(1:n+1)=temp%isa(1:n+1)
           matAjb(p,j)%nnz=nnz1
           deallocate(temp%sa)
           deallocate(temp%jsa)
           deallocate(temp%isa)
        end if
     end if
  end do


 do j=1,nbprocs   !! scan along the columns minus diagonal
     if ((mapA(p,j)>0).and.(j/=p)) then !! non-zero block to transpose
        n=Nsize(p)
        m=Nsize(j) 
        nnz1=mapA(p,j)       
        !! create the new space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           matAjb(j,p)%n=m
           matAjb(j,p)%nnz=nnz1
           allocate(matAjb(j,p)%sa(nnz1))
           allocate(matAjb(j,p)%jsa(nnz1))
           allocate(matAjb(j,p)%isa(m+1))
           matAjb(j,p)%nnz=nnz1
           call zcsr_transpose(n,m,matAjb(p,j)%sa,matAjb(p,j)%isa,matAjb(p,j)%jsa,matAjb(j,p)%sa,matAjb(j,p)%isa,matAjb(j,p)%jsa) 
        end if
     end do
  
!!! generate the global nnz matrix
  mapA(1:nbprocs,1:nbprocs)=0
  do i=1,nbprocs
     mapA(p,i)=matAjb(p,i)%nnz
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,mapA,nbprocs*nbprocs,MPI_INTEGER,MPI_SUM,MY_COMM_WORLD,code) 


end subroutine zlbcsr_uplo_to_csr




subroutine  dlbcsr_distribute_row(Nsize,mapA,matAjb,matAj,MY_COMM_WORLD,nbprocs)
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine distribute a matrix by rows- as a condition the matrix must first be distributed using the local block csr format.
  ! This is the real double presision version
  !
  !   Example using 3 MPI:

  !    |_mpi1_|              A11 A12 A13
  !  A=|_mpi2_|       <==    A21 A22 A23    
  !    |_mpi3_|              A31 A32 A33
  !
  ! Nsize     (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !
  ! matAj  (output) TYPE(DCSR)-
  !        local by row- see explanations above
  !
  ! NEW_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr) :: matAj
  type(dcsr),dimension(nbprocs,*) :: matAjb
  integer :: i,j,rank,p,code,k,tok,row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)

  p=rank+1

  matAj%n=Nsize(p)
  matAj%nnz=sum(mapA(p,1:nbprocs))
  allocate(matAj%isa(Nsize(p)+1))
  allocate(matAj%jsa(matAj%nnz))
  allocate(matAj%sa(matAj%nnz))

  matAj%isa(1)=1
  tok=0
  do i=1,matAj%n ! loop over all rows
     row=0 !
     do j=1,nbprocs ! loop over blocks
        if (mapA(p,j)>0) then ! condition to include block in row distribution
           do k=matAjb(p,j)%isa(i),matAjb(p,j)%isa(i+1)-1 ! loop over element in row
              tok=tok+1
              row=row+1!
              matAj%sa(tok)=matAjb(p,j)%sa(k)
              matAj%jsa(tok)=matAjb(p,j)%jsa(k)
              if (j>1) matAj%jsa(tok)=matAj%jsa(tok)+sum(Nsize(1:j-1))
           end do
        end if
     enddo
     matAj%isa(i+1)=matAj%isa(i)+row
  enddo

end subroutine dlbcsr_distribute_row


subroutine  zlbcsr_distribute_row(Nsize,mapA,matAjb,matAj,MY_COMM_WORLD,nbprocs)
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine distribute a matrix by rows- as a condition the matrix must first be distributed using the local block csr format.
  ! This is the complex double presision version
  !
  !   Example using 3 MPI:

  !    |_mpi1_|              A11 A12 A13
  !  A=|_mpi2_|       <==    A21 A22 A23    
  !    |_mpi3_|              A31 A32 A33
  !
  ! Nsize     (input) INTEGER(:) - size nbprocs
  !            global - store the size of each row blocks
  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
  !            global - contains the nnz of each submatrix csr blocks
  !        Example:
  !        Each blocks contains a csr matrix, or could be equal to zero;
  !        the number of nnz is provided using the mapping map such that if
  !             |23  12  0 |
  !        map= |12  18  0 |
  !             |0    0  16|
  !
  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  ! matAjb (input) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- block csr format
  !
  ! matAj  (output) TYPE(ZCSR)-
  !        local by row- see explanations above
  !
  ! NEW_COMM_WORLD  (input) INTEGER - COMMUNICATOR
  !        
  ! nbprocs (input) INTEGER - # of MPI processors
  !
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================
  integer,dimension(*) :: Nsize 
  integer :: MY_COMM_WORLD,nbprocs
  integer,dimension(nbprocs,*) :: mapA
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr) :: matAj
  type(zcsr),dimension(nbprocs,*) :: matAjb
  integer :: i,j,rank,p,code,k,tok,row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)

  p=rank+1

  matAj%n=Nsize(p)
  matAj%nnz=sum(mapA(p,1:nbprocs))
  allocate(matAj%isa(Nsize(p)+1))
  allocate(matAj%jsa(matAj%nnz))
  allocate(matAj%sa(matAj%nnz))

  matAj%isa(1)=1
  tok=0
  do i=1,matAj%n ! loop over all rows
     row=0
     do j=1,nbprocs ! loop over blocks
        if (mapA(p,j)>0) then ! condition to include block in row distribution
           do k=matAjb(p,j)%isa(i),matAjb(p,j)%isa(i+1)-1 ! loop over element in row
              tok=tok+1
              row=row+1 !increment
              matAj%sa(tok)=matAjb(p,j)%sa(k)
              matAj%jsa(tok)=matAjb(p,j)%jsa(k)
              if (j>1) matAj%jsa(tok)=matAj%jsa(tok)+sum(Nsize(1:j-1))
           end do
        end if
     enddo
     matAj%isa(i+1)=matAj%isa(i)+row
  enddo

end subroutine zlbcsr_distribute_row





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dlbcsrmm(UPLO,TRANS,Nsize,map,M0,alpha,matAjb,Xj,beta,Yj,MY_COMM_WORLD,nbprocs)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      y=alpha*o(A)*x+beta*y   (square matrix)
  !
  !      This is the real double precision version
  !
  !      This is a MPI mat-vec version using MY_COMM_WORLD communicator with nbprocs for #MPI 
  ! 
  !      The square matrix A is stored locally in multi-block CSR format matAjb (data structure)
  !      and Xj, Yj  are the local row distributed solution and rhs, respectively.
  !
  !      Example using 3 MPI:
  !      A11 A12 A13
  !      A21 A22 A23
  !      A31 A32 A33
  !
  !      The block Aij is stored in both process i and j (cross configuration- useful for performing transpose mat-vec). A given processor stores up to 2*N-1 blocks.
  !     
  !      Each blocks contains a csr matrix, or could be equal to zero
  !      the number of nnz is provided using the mapping map such that if
  !           |23  12  0 |
  !      map= |12  18  0 |
  !           |0    0  16|
  !
  !      it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in block format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided (N must be equal to M - symmetric matrix)
  !       if UPLO='U' only the upper part of the matrix A is provided (N must be equal to M - symmetrix matrix) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !
  !  Nsize (input) INTEGER(:) - size nbprocs
  !        global - store the size of each row blocks
  !  map (input) INTEGER(:,:) - size nbprocs*nbprocs
  !        global - store the nnz of each submatrix blocks
  !  M0  (input) INTEGER
  !        The number of column of matrices X and Y
  !  alpha (input) DOUBLE PRECISION
  ! matAjb (input) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
  !        local- see explanations above
  !
  !  Xj   (input) DOUBLE PRECISION
  !        local (row distributed)- matrix of size (Nsize(j)xrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  Yj    (input/output) DOUBLE PRECISION
  !        local (row distributed) matrix of size (Nsize(j)xrhs).
  !        On exit contains the solution.
  !
  !  MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !  nbprocs        (input) INTEGER
  !                  #MPI
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-------------------------------------
  character(len=1) :: UPLO,TRANS
  integer :: MY_COMM_WORLD,M0,nbprocs
  integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
  double precision,dimension(*) :: Xj,Yj
  double precision :: alpha,beta
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr),dimension(nbprocs,*) :: matAjb
  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer ::i,p,rank,code,tag,request
  double precision, dimension(:),allocatable :: dauxj
  double precision, parameter :: DONE=1.0d0

  tag=100     

!!!!!!!
  !call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!
  !call wallocate_1d(dauxj,maxval(Nsize(1:nb_procs))*M0,infoloc)
  allocate(dauxj(1:maxval(Nsize(1:nbprocs))*M0))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=="F") then !!! Full format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (TRANS=="N") then !! normal mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','N',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wdcsrmm('F','N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,status,code)
              call wdcsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,DONE,Yj)
           end if
        end do


     elseif (TRANS=="T") then !! transpose mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','T',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wdcsrmm('F','T',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,status,code)
              call wdcsrmm('F','T',Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,DONE,Yj)
           end if
        end do

     end if


  else ! ====> symmetric matrix-- L/U format (using full format triangular) !!!!! A=L+L^T-D or A=U+U^T-D
     ! normal or transpose option ok (the same result)

     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do

     !call wdcsrmm(UPLO,'N',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
     call wdcsrmm(UPLO,'N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,status,code)
           call wdcsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,DONE,Yj)
        end if
     end do

!!!!

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do
!!!diag blocks done already

     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_PRECISION,i-1,tag, MY_COMM_WORLD ,status,code)
           call wdcsrmm('F','T',Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,DONE,Yj)
        end if
     end do


  end if


  deallocate(dauxj)


end subroutine dlbcsrmm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zlbcsrmm(UPLO,TRANS,Nsize,map,M0,alpha,matAjb,Xj,beta,Yj,MY_COMM_WORLD,nbprocs)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      y=alpha*o(A)*x+beta*y   (square matrix)
  !
  !      This is the complex double precision version
  !
  !      This is a MPI mat-vec version using MY_COMM_WORLD communicator with nbprocs for #MPI 
  ! 
  !      The square matrix A is stored locally in multi-block CSR format matAjb (data structure)
  !      and Xj, Yj  are the local row distributed solution and rhs, respectively.
  !
  !      Example using 3 MPI:
  !      A11 A12 A13
  !      A21 A22 A23
  !      A31 A32 A33
  !
  !      The block Aij is stored in both process i and j (cross configuration- useful for performing transpose mat-vec). A given processor stores up to 2*N-1 blocks.
  !     
  !      Each blocks contains a csr matrix, or could be equal to zero
  !      the number of nnz is provided using the mapping map such that if
  !           |23  12  0 |
  !      map= |12  18  0 |
  !           |0    0  16|
  !
  !      it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in block format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided (N must be equal to M - symmetric matrix)
  !       if UPLO='U' only the upper part of the matrix A is provided (N must be equal to M - symmetrix matrix) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode*
  !        * not fully tested for UPLO=L,U in particular
  !
  !  Nsize (input) INTEGER(:) - size nbprocs
  !        global - store the size of each row blocks
  !  map (input) INTEGER(:,:) - size nbprocs*nbprocs
  !        global - store the nnz of each submatrix blocks
  !  M0  (input) INTEGER
  !        The number of column of matrices X and Y
  !  alpha (input) COMPLEX DOUBLE PRECISION
  ! matAjb (input) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- see explanations above
  !
  !  Xj   (input) COMPLEX DOUBLE PRECISION
  !        local (row distributed)- matrix of size (Nsize(j)xrhs)
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  Yj    (input/output) COMPLEX DOUBLE PRECISION
  !        local (row distributed) matrix of size (Nsize(j)xrhs).
  !        On exit contains the solution.
  !
  !  MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !  nbprocs        (input) INTEGER
  !                  #MPI
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-------------------------------------
  character(len=1) :: UPLO,TRANS
  integer :: MY_COMM_WORLD,M0,nbprocs
  integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
  complex(kind=kind(1.0d0)),dimension(*) :: Xj,Yj
  complex(kind=kind(1.0d0)) :: alpha,beta
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matAjb
  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer ::i,p,rank,code,tag,request
  complex(kind=kind(1.0d0)), dimension(:),allocatable :: dauxj
  character(len=1) :: TRANS2
  complex(kind=kind(1.0d0)), parameter :: ZONE=(1.0d0,0.0d0)

  tag=100     

!!!!!!!
  !call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!
  !call wallocate_1d(dauxj,maxval(Nsize(1:nb_procs))*M0,infoloc)
  allocate(dauxj(1:maxval(Nsize(1:nbprocs))*M0))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=="F") then !!! Full format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (TRANS=="N") then !! normal mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','N',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wzcsrmm('F','N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
              call wzcsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,ZONE,Yj)
           end if
        end do


     elseif ((TRANS=="T").or.(TRANS=="C")) then !! transpose mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','T',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wzcsrmm('F',TRANS,Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
              call wzcsrmm('F',TRANS,Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,ZONE,Yj)
           end if
        end do

     end if


  else ! ====> symmetric matrix-- L/U format (using full format triangular) !!!!! A=L+L^T-D or A=U+U^T-D
     ! normal or transpose option ok (the same result)
     TRANS2='T' ! regular transpose
     if (TRANS=='C') TRANS2='C'

     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do


     ! diag blocks
        call wzcsrmm(UPLO,'N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
!     call wzcsrmm(UPLO,TRANS2,Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
           call wzcsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,ZONE,Yj)
        end if
     end do

!!!!

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do
!!!diag blocks done already
     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
           call wzcsrmm('F',TRANS2,Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,ZONE,Yj)
        end if
     end do


  end if


  deallocate(dauxj)


end subroutine zlbcsrmm




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zhlbcsrmm(UPLO,TRANS,Nsize,map,M0,alpha,matAjb,Xj,beta,Yj,MY_COMM_WORLD,nbprocs)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !

  !      y=alpha*o(A)*x+beta*y   (square matrix)
  !
  !      This is the complex double precision version
  !
  !      This is a MPI mat-vec version using MY_COMM_WORLD communicator with nbprocs for #MPI 
  ! 
  !      The square matrix A is stored locally in multi-block CSR format matAjb (data structure)
  !      and Xj, Yj  are the local row distributed solution and rhs, respectively.
  !
  !      Example using 3 MPI:
  !      A11 A12 A13
  !      A21 A22 A23
  !      A31 A32 A33
  !
  !      The block Aij is stored in both process i and j (cross configuration- useful for performing transpose mat-vec). A given processor stores up to 2*N-1 blocks.
  !     
  !      Each blocks contains a csr matrix, or could be equal to zero
  !      the number of nnz is provided using the mapping map such that if
  !           |23  12  0 |
  !      map= |12  18  0 |
  !           |0    0  16|
  !
  !      it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in block format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided (N must be equal to M - symmetric matrix)
  !       if UPLO='U' only the upper part of the matrix A is provided (N must be equal to M - symmetrix matrix) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode*
  !        * not fully tested for UPLO=L,U in particular
  !
  !  Nsize (input) INTEGER(:) - size nbprocs
  !        global - store the size of each row blocks
  !  map (input) INTEGER(:,:) - size nbprocs*nbprocs
  !        global - store the nnz of each submatrix blocks
  !  M0  (input) INTEGER
  !        The number of column of matrices X and Y
  !  alpha (input) COMPLEX DOUBLE PRECISION
  ! matAjb (input) Matrix of TYPE(ZCSR)- size nbprocs*nbprocs
  !        local- see explanations above
  !
  !  Xj   (input) COMPLEX DOUBLE PRECISION
  !        local (row distributed)- matrix of size (Nsize(j)xrhs)
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  Yj    (input/output) COMPLEX DOUBLE PRECISION
  !        local (row distributed) matrix of size (Nsize(j)xrhs).
  !        On exit contains the solution.
  !
  !  MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !  nbprocs        (input) INTEGER
  !                  #MPI
  !====================================================================
  ! Eric Polizzi 2018-2019
  !====================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-------------------------------------
  character(len=1) :: UPLO,TRANS
  integer :: MY_COMM_WORLD,M0,nbprocs
  integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
  complex(kind=kind(1.0d0)),dimension(*) :: Xj,Yj
  complex(kind=kind(1.0d0)) :: alpha,beta
  type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matAjb
  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer ::i,p,rank,code,tag,request
  complex(kind=kind(1.0d0)), dimension(:),allocatable :: dauxj
  character(len=1) :: TRANS2
  complex(kind=kind(1.0d0)), parameter :: ZONE=(1.0d0,0.0d0)

  tag=100     

!!!!!!!
  !call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!
  !call wallocate_1d(dauxj,maxval(Nsize(1:nb_procs))*M0,infoloc)
 
  allocate(dauxj(1:maxval(Nsize(1:nbprocs))*M0))
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=="F") then !!! Full format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (TRANS=="N") then !! normal mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','N',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wzcsrmm('F','N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
              call wzcsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,ZONE,Yj)
           end if
        end do


     elseif ((TRANS=="T").or.(TRANS=="C")) then !! transpose mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','T',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wzcsrmm('F',TRANS,Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
              call wzcsrmm('F',TRANS,Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,ZONE,Yj)
           end if
        end do

     end if


  else ! ====> hermitian matrix-- L/U format (using full format triangular) !!!!! A=L+L^T-D or A=U+U^T-D
     ! normal or transpose option ok (the same result)
     TRANS2='C' ! transpose conjugate
     if (TRANS=='C') TRANS2='T'

     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do


     ! diag blocks
        call wzhcsrmm(UPLO,'N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
!     call wzcsrmm(UPLO,TRANS2,Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
           call wzcsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,ZONE,Yj)
        end if
     end do

!!!!

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do
!!!diag blocks done already
     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_DOUBLE_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
           call wzcsrmm('F',TRANS2,Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,ZONE,Yj)
        end if
     end do


  end if


  deallocate(dauxj)


end subroutine zhlbcsrmm




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine clbcsrmm(UPLO,TRANS,Nsize,map,M0,alpha,matAjb,Xj,beta,Yj,MY_COMM_WORLD,nbprocs)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      y=alpha*o(A)*x+beta*y   (square matrix)
  !
  !      This is the complex single precision version
  !
  !      This is a MPI mat-vec version using MY_COMM_WORLD communicator with nbprocs for #MPI 
  ! 
  !      The square matrix A is stored locally in multi-block CSR format matAjb (data structure)
  !      and Xj, Yj  are the local row distributed solution and rhs, respectively.
  !
  !      Example using 3 MPI:
  !      A11 A12 A13
  !      A21 A22 A23
  !      A31 A32 A33
  !
  !      The block Aij is stored in both process i and j (cross configuration- useful for performing transpose mat-vec). A given processor stores up to 2*N-1 blocks.
  !     
  !      Each blocks contains a csr matrix, or could be equal to zero
  !      the number of nnz is provided using the mapping map such that if
  !           |23  12  0 |
  !      map= |12  18  0 |
  !           |0    0  16|
  !
  !      it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in block format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided (N must be equal to M - symmetric matrix)
  !       if UPLO='U' only the upper part of the matrix A is provided (N must be equal to M - symmetrix matrix) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode*
  !        * not fully tested for UPLO=L,U in particular
  !
  !  Nsize (input) INTEGER(:) - size nbprocs
  !        global - store the size of each row blocks
  !  map (input) INTEGER(:,:) - size nbprocs*nbprocs
  !        global - store the nnz of each submatrix blocks
  !  M0  (input) INTEGER
  !        The number of column of matrices X and Y
  !  alpha (input) COMPLEX SINGLE PRECISION
  ! matAjb (input) Matrix of TYPE(CCSR)- size nbprocs*nbprocs
  !        local- see explanations above
  !
  !  Xj   (input) COMPLEX SINGLE PRECISION
  !        local (row distributed)- matrix of size (Nsize(j)xrhs)
  !
  !  beta  (input) COMPLEX SINGLE PRECISION
  !
  !  Yj    (input/output) COMPLEX SINGLE PRECISION
  !        local (row distributed) matrix of size (Nsize(j)xrhs).
  !        On exit contains the solution.
  !
  !  MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !  nbprocs        (input) INTEGER
  !                  #MPI
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-------------------------------------
  character(len=1) :: UPLO,TRANS
  integer :: MY_COMM_WORLD,M0,nbprocs
  integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
  complex,dimension(*) :: Xj,Yj
  complex :: alpha,beta
  type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr),dimension(nbprocs,*) :: matAjb
  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer ::i,p,rank,code,tag,request
  complex, dimension(:),allocatable :: dauxj
  character(len=1) :: TRANS2
  complex, parameter :: CONE=(1.0d0,0.0d0)

  tag=100     

!!!!!!!
  !call MPI_COMM_SIZE(MY_COMM_WORLD,nbprocs,code)
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!
  !call wallocate_1d(dauxj,maxval(Nsize(1:nb_procs))*M0,infoloc)
  allocate(dauxj(1:maxval(Nsize(1:nbprocs))*M0))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=="F") then !!! Full format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (TRANS=="N") then !! normal mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','N',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wccsrmm('F','N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
              call wccsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,CONE,Yj)
           end if
        end do


     elseif ((TRANS=="T").or.(TRANS=="C")) then !! transpose mat-vec

        do i=1,nbprocs
           if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
              call MPI_ISEND(Xj,Nsize(p)*M0,MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
           endif
        end do
        !call wdcsrmm('F','T',matAjb(p,p)%n,matAjb(p,p)%n,M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
        call wccsrmm('F',TRANS,Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

        do i=1,nbprocs
           if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
              call MPI_RECV(dauxj,M0*Nsize(i), MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
              call wccsrmm('F',TRANS,Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,CONE,Yj)
           end if
        end do

     end if


  else ! ====> symmetric matrix-- L/U format (using full format triangular) !!!!! A=L+L^T-D or A=U+U^T-D
     ! normal or transpose option ok (the same result)
     TRANS2='T' ! regular transpose
     if (TRANS=='C') TRANS2='C'

     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do


     ! diag blocks
        call wccsrmm(UPLO,'N',Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)
!     call wzcsrmm(UPLO,TRANS2,Nsize(p),Nsize(p),M0,alpha,matAjb(p,p)%sa,matAjb(p,p)%isa,matAjb(p,p)%jsa,Xj,beta,Yj)

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
           call wccsrmm('F','N',Nsize(p),Nsize(i),M0,alpha,matAjb(p,i)%sa,matAjb(p,i)%isa,matAjb(p,i)%jsa,dauxj,CONE,Yj)
        end if
     end do

!!!!

     do i=1,nbprocs
        if ((p/=i).and.(map(p,i)/=0)) then !! send along the non-zero row
           call MPI_ISEND(Xj,Nsize(p)*M0,MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,request,code)
        endif
     end do
!!!diag blocks done already
     do i=1,nbprocs
        if ((p/=i).and.(map(i,p)/=0)) then !! receive along the column
           call MPI_RECV(dauxj,M0*Nsize(i), MPI_COMPLEX,i-1,tag, MY_COMM_WORLD ,status,code)
           call wccsrmm('F',TRANS2,Nsize(i),Nsize(p),M0,alpha,matAjb(i,p)%sa,matAjb(i,p)%isa,matAjb(i,p)%jsa,dauxj,CONE,Yj)
        end if
     end do


  end if


  deallocate(dauxj)


end subroutine clbcsrmm







subroutine  pdqrgs(N,M,q,r,MY_COMM_WORLD)
  !! Purpose: Compute QR factorization using Gram-Schmidt algorithm
  !           using modified GS in-place 
  ! N    (input) INTEGER  local #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) REAL DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) REAL DOUBLE PRECISION (M,M)- R matrix factors
  !
  !  MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
   !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  integer :: N,M,MY_COMM_WORLD
  double precision,dimension(N,*) :: q
  double precision,dimension(M,*) :: R
  !---
  double precision :: DNRM2
  double precision:: ddot
  integer :: k,j,rank,p,code

  !!!!
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
  !!!! 
   
do k=1,M
   R(k,k)=DNRM2(N,q(1,k),1)**2
   call MPI_ALLREDUCE(MPI_IN_PLACE,R(k,k),1,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
   R(k,k)=sqrt(R(k,k))
   call DSCAL(N,1.0d0/R(k,k), Q(1,k), 1)
   do j=k+1,M
      R(k,j)=ddot(N,Q(1,k),1,Q(1,j),1)
       call MPI_ALLREDUCE(MPI_IN_PLACE,R(k,j),1,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
      call DAXPY(N,-R(k,j),Q(1,k),1,Q(1,j),1)     
enddo
enddo



end subroutine pdqrgs







subroutine  pzqrgs(N,M,q,r,MY_COMM_WORLD)
  !! Purpose: Compute QR factorization using Gram-Schmidt algorithm
  !           using modified GS in-place 
  ! N    (input) INTEGER  local #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors
  !
  !  MY_COMM_WORLD  (input) INTEGER
  !                  MPI communicator
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
   !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  integer :: N,M,MY_COMM_WORLD
  complex(kind=kind(1.0d0)),dimension(N,*) :: q
  complex(kind=kind(1.0d0)),dimension(M,*) :: R
  !---
  double precision :: DZNRM2
  !complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
  integer :: k,j,rank,p,code

  !!!!
 call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
  !!!! 
   
do k=1,M
   R(k,k)=(1.0d0,0.0d0)*DZNRM2(N,q(1,k),1)**2
   call MPI_ALLREDUCE(MPI_IN_PLACE,R(k,k),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
    R(k,k)=sqrt(R(k,k))
   call ZSCAL(N,(1.0d0,0.0d0)/R(k,k), Q(1,k), 1)
   do j=k+1,M
      !R(k,j)=zdotc(N,Q(1,k),1,Q(1,j),1)
       R(k,j)=dot_product(Q(1:N,k),Q(1:N,j))
       call MPI_ALLREDUCE(MPI_IN_PLACE,R(k,j),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
      call ZAXPY(N,-R(k,j),Q(1,k),1,Q(1,j),1)     
enddo
enddo



end subroutine pzqrgs







subroutine  pdbicgstab(opt,uplo,trans,Nsize,map,M0,matAjb,fj,xj,res,nbitmax,epso,comd,MY_COMM_WORLD,nbprocs,info)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  integer :: opt, MY_COMM_WORLD, M0, nbprocs
  character(len=1) :: uplo,trans
 integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
   type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr),dimension(nbprocs,*) :: matAjb
  double precision, dimension(*) :: fj
  double precision, dimension(*):: xj
  double precision, dimension(M0) :: res
  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
integer :: code,rank,p
  double precision,dimension(:,:),allocatable :: work,work2
  integer :: ijob,j1,j2,k,i
 
!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!

  allocate(work(Nsize(p),7*M0))
  allocate(work2(M0,4))   

  ijob=-1
  do while (ijob/=0)
     call pdbicgstab_rci(ijob,Nsize(p),M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,MY_COMM_WORLD,info) 
     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input

     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)
       
        call dlbcsrmm(UPLO,TRANS,Nsize,map,M0,1.0d0,matAjb,work(1,j1),0.0d0,work(1,j2),MY_COMM_WORLD,nbprocs)

     end select
  end do




  deallocate(work,work2)


end subroutine pdbicgstab



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  pdbicgstab_rci(ijob,N,M0,fj,xj,work1,work2,j1,j2,res,nbitmax,epso,comd,MY_COMM_WORLD,info) 
  implicit none
    !-------------------------------------
  include 'mpif.h'
  !-------------------------------------

  integer :: j1,j2,MY_COMM_WORLD 
  integer :: ijob
  integer :: N,M0
  double precision, dimension(N,*) :: fj,xj,work1
  double precision, dimension(M0,*):: work2 
  double precision, dimension(M0) :: res


  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision, parameter :: ONE=1.0d0,ZERO=0.0d0

  double precision ::ares
  integer :: rank,i
  integer :: rj,rb,pp,pb,v,sb,t  !!! column index for the work1 array 
  integer:: rho_1,alpha,omega,ind !!! column index for the work2 array

  double precision,dimension(M0) :: resf
  double precision,dimension(M0) :: rho_2,aux0,aux1,beta
  integer,save :: nbit
  integer,save :: probe

  double precision :: ddot
  double precision :: DNRM2

  logical :: loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: code
  
  integer :: b1,b2


  b1=1
  b2=M0

!!!!!!!!!! column index for the work1 array
  rj=1
  rb=M0+1
  pp=2*M0+1
  pb=3*M0+1
  v=4*M0+1
  sb=5*M0+1
  t=6*M0+1
!!!!!!!!!! column index for the work2 array
  rho_1=1
  alpha=2
  omega=3
  ind=4 ! if /= (0,0) do not include residual calculations in the max(norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rank=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !    rho_1(M0)   ! global  >> work2(:,1)
  !    rho_2(M0))  ! local
  !    beta(M0)    ! local
  !    alpha(M0)   ! global >> work2(:,2)
  !    omega(M0)   ! global >> work2(:,3)
  !    aux0(M0)    ! local
  !    aux1(M0)    ! local   
  !    resf(M0)    ! local (use res to store it)

!!!!!!!!!!!!

  info=0
  loop=.true.   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ijob==-1) then
!!!!!!!!!!!!!!! initialize index for residual calculation
     if (M0==1) then
        work2(1,ind)=ONE
     else
        do i=1,M0
           work2(i,ind)=res(i)*ONE
        enddo
!!! special case (user enter 0 everywhere ==> 1 everywhere)
        if (sum(work2(:,ind))==ZERO) work2(:,ind)=ONE

     end if

!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax      
!!!! let us temporaly save xj into work1(1,rj)
     call DCOPY(N*M0,xj(1,1),1,work1(1,rj),1)


!!!! Compute A*xj  ----- A*work1(1,rj)=>work1(1,rb)
     j1=rj
     j2=rb
     ijob=2
     probe=0
     !work1(1:N,rb:rb+M0-1)=ZEROC ! if initial guess 0
     return
  end if


  if (probe==0) then

!!!! let us temporaly save fj into work1(1,rj)
     call DCOPY(N*M0,fj(1,1),1,work1(1,rj),1)

!!!!!!! compute residual r=f-A*x
     do i=1,M0   
        call DAXPY(N,-ONE,work1(1,i+rb-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!! CHECK THE  NORM of rj
     do i=1,M0
        res(i)=DNRM2(N,work1(1,i+rj-1),1)**2
        resf(i)=DNRM2(N,fj(1,i+rj-1),1)**2
     enddo
     call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     call MPI_ALLREDUCE(MPI_IN_PLACE,resf,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
        resf(i)=sqrt(resf(i)) ! norm L2
     enddo
    
     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))

     res=resf ! save resf value
     if (comd) print *,'rel. residual before iteration',ares,res(1),resf(1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! BIPCG-STAB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! choose rb
     call DCOPY(N*M0,work1(1,rj),1,work1(1,rb),1)

     nbit=0

     probe=1

  endif


  if (probe==1) then
     nbit=nbit+1
     if (nbit>1) rho_2=work2(1:M0,rho_1) 

!!!!!!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0
        work2(i,rho_1)=ddot(N,work1(1,rb+i-1),1,work1(1,rj+i-1),1)
     end do
 call MPI_ALLREDUCE(MPI_IN_PLACE,work2(1,rho_1),M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)

      
!!!!!!!!!!!! TEST
     do i=1,M0
        if (work2(i,rho_1)==ZERO) then
           info= -201
           if ((comd).and.(rank==0)) print*,'ATTENTION ----- BICG-STAB OUT FAILED, rho_1=0 !!'
           nbitmax=nbit 
           ijob=0
           return              
        end if
     end do

!!!!!!!!!!!! CONDITION
     IF (nbit==1) THEN
        call DCOPY(N*M0,work1(1,rj),1,work1(1,pp),1)
     ELSE
        beta=(work2(1:M0,rho_1)/rho_2)*(work2(1:M0,alpha)/work2(1:M0,omega))
        do i=1,M0   !!!!!!<<< to optimize
           call DAXPY(N,-work2(i,omega),work1(1,v+i-1),1,work1(1,pp+i-1),1)
           call DSCAL(N,beta(i),work1(1,pp+i-1),1)
           call DAXPY(N,ONE,work1(1,i+rj-1),1,work1(1,i+pp-1),1)
           !pp(:,i)=r(:,i)+beta(i)*(pp(:,i)-omega(i)*v(:,i))
        end do
     END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve  M*pb=pp(k) ---> pb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call DCOPY(N*M0,work1(1,pp),1,work1(1,pb),1)
     j1=pb
     j2=pp
     ijob=1
     probe=2
     return
  endif

  if (probe==2) then  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by pb, results in v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=pb
     j2=v
     ijob=2
     probe=3
     return
  end if

  if (probe==3) then
!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0!! different RHS
        aux0(i)=ddot(N,work1(1,i+rb-1),1,work1(1,i+v-1),1)
        !aux0(i)=sum(work1(1:N,i+rb-1)*work1(1:N,i+v-1))
     end do
 call MPI_ALLREDUCE(MPI_IN_PLACE,aux0,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)

!!!!!!
     work2(1:M0,alpha)=work2(1:M0,rho_1)/aux0(1:M0)
     do i=1,M0   
        !ss(:,i)=r(:,i)-alpha(i)*v(:,i) !!!!!!!!! use rj to store ss !!!!!!!!!!!
        call DAXPY(N,-work2(i,alpha),work1(1,i+v-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!!! CHECK THE  NORM of ss
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DNRM2(N,work1(1,i+rj-1),1)**2
     enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))



     IF (ares<epso) then
        do i=1,M0  
           !      xj(:,i)=xj(:,i)+alpha(i)*pb(:,i) 
           call DAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        end do
        if ((comd).and.(rank==0)) print *,(nbit-1)*1.0d0+0.5d0,ares
        info=0
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares             
        return
     end IF

     res=resf ! save resf value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve M*sb=ss ---> sb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     call DCOPY(N*M0,work1(1,rj),1,work1(1,sb),1)
     j1=sb
     j2=rj
     ijob=1
     probe=4
     return
  end if


  if (probe==4) then          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by sb, results in t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=sb
     j2=t
     ijob=2
     probe=5
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! SCALAR PRODUCTS on different RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (probe==5) then

     do i=1,M0!! different RHS
        aux0(i)=ddot(N,work1(1,i+t-1),1,work1(1,i+rj-1),1)
        aux1(i)=ddot(N,work1(1,i+t-1),1,work1(1,i+t-1),1)
     end do
     call MPI_ALLREDUCE(MPI_IN_PLACE,aux0,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     call MPI_ALLREDUCE(MPI_IN_PLACE,aux1,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)

!!!!!!
     work2(1:M0,omega)=aux0/aux1
     do i=1,M0
        if (work2(i,omega)==ZERO) then
           info= -202 
           if ((comd).and.(rank==0)) print *,'ATTENTION ----- BICG-STAB OUT FAILED, omega=0'
           nbitmax=nbit
           ijob=0
           return 
           !fail=.true.
        end if
     end do

     do i=1,M0
        ! xj(:,i)=xj(:,i)+alpha(i)*pb(:,i)+omega(i)*sb(:,i)
        call DAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        call DAXPY(N,work2(i,omega),work1(1,i+sb-1),1,xj(1,i),1) 
     end do

     do i=1,M0
        !   r(:,i)=ss(:,i)-omega(i)*t(:,i)
        !call ZCOPY(N,ss(1,i),1,rj(1,i),1)   !!! ss is rj
        call DAXPY(N,-work2(i,omega),work1(1,i+t-1),1,work1(1,i+rj-1),1) 
     end do


!!!!!!!!! CHECK THE  NORM of rj
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DNRM2(N,work1(1,i+rj-1),1)**2
     enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
     enddo



     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZERO) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))
     if ((comd).and.(rank==0)) print *,nbit,ares


!!!! test for the other loop
     if (.not.((nbit<nbitmax).and.(ares.gt.epso))) loop=.false. 

     if (loop) then
        ijob=-2
        probe=1
        res=resf ! save resf 
     else
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares      
     end if
     return

  end if


end subroutine pdbicgstab_rci










subroutine  pzbicgstab(opt,uplo,trans,Nsize,map,M0,matAjb,fj,xj,res,nbitmax,epso,comd,MY_COMM_WORLD,nbprocs,info)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  integer :: opt, MY_COMM_WORLD, M0, nbprocs
  character(len=1) :: uplo,trans
 integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
   type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matAjb
  complex(kind=kind(1.0d0)), dimension(*) :: fj
  complex(kind=kind(1.0d0)), dimension(*):: xj
  double precision, dimension(M0) :: res
  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
integer :: code,rank,p
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: work,work2
  integer :: ijob,j1,j2,k,i
 
!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!


  
  allocate(work(Nsize(p),7*M0))
  allocate(work2(M0,4))   

  ijob=-1
  do while (ijob/=0)
     call pzbicgstab_rci(ijob,Nsize(p),M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,MY_COMM_WORLD,info) 
!print *,'bicg',ijob
     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input

     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)
        !print *,'w1',p,sum(abs(work(1:Nsize(p),j1)))
        !print *,'w2',p,sum(abs(work(1:Nsize(p),j2)))
        !work(1:Nsize(p),j2+1:j2+M0-1)=(0.0d0,0.0d0)
!        if (p==1) print *,'map',map(1:2,1:2)
!        print *,'Nsize',p,Nsize(p)
        ! print *,'workj1',sum(abs(work(1:Nsize(p),j1)))
     !   if (p==1) matAjb(2,1)%sa(matAjb(2,1)%isa(Nsize(2)+1)-1)=(0.0d0,0.0d0)
      !  if (p==2) matAjb(2,1)%sa(matAjb(2,1)%isa(Nsize(2)+1)-1)=(0.0d0,0.0d0)
        !if (p==1)
        
!         matAjb(1,2)%sa(1:matAjb(1,2)%isa(Nsize(1)+1)-1)=(1.0d0,0.0d0)
  !       matAjb(2,1)%sa(1:matAjb(2,1)%isa(Nsize(2)+1)-1)=(1.0d0,0.0d0)

     !   if (p==1) print *,'a21',  sum(abs(matAjb(2,1)%sa(1:matAjb(2,1)%isa(Nsize(2)+1)-1)))
      !  if (p==2) print *,'b21',  sum(abs(matAjb(2,1)%sa(1:matAjb(2,1)%isa(Nsize(2)+1)-1)))
        

        
! if (p==1) matAjb(1,1)%sa(1:matAjb(1,1)%isa(Nsize(1)+1)-1)=(0.0d0,0.0d0)
! if (p==2) matAjb(2,2)%sa(1:matAjb(2,2)%isa(Nsize(2)+1)-1)=(0.0d0,0.0d0)


        !if (p==2) matAjb(1,2)%sa(matAjb(1,2)%isa(Nsize(1)+1)-1)=(0.0d0,0.0d0)
        
!        do i=1,p
!          print *,p,i,sum(abs(matAjb(p,i)%sa(1:matAjb(p,i)%isa(Nsize(p)+1)-1)))
!enddo
       !  print *,'workj1',sum(abs(work(1:Nsize(p),j1))),work(1000,j1)
   
        call zlbcsrmm(UPLO,TRANS,Nsize,map,M0,(1.0d0,0.0d0),matAjb,work(1,j1),(0.0d0,0.0d0),work(1,j2),MY_COMM_WORLD,nbprocs)
        !print *,'done'
!         print *,'workj2',sum(abs(work(1:Nsize(p),j2))),work(1000,j2)

     end select
  end do




  deallocate(work,work2)


end subroutine pzbicgstab







subroutine  pzhbicgstab(opt,uplo,trans,Nsize,map,M0,matAjb,fj,xj,res,nbitmax,epso,comd,MY_COMM_WORLD,nbprocs,info)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  integer :: opt, MY_COMM_WORLD, M0, nbprocs
  character(len=1) :: uplo,trans
 integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
   type zcsr
     integer ::n,m,nnz
     complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr),dimension(nbprocs,*) :: matAjb
  complex(kind=kind(1.0d0)), dimension(*) :: fj
  complex(kind=kind(1.0d0)), dimension(*):: xj
  double precision, dimension(M0) :: res
  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
integer :: code,rank,p
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: work,work2
  integer :: ijob,j1,j2,k,i
 
!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!


  
  allocate(work(Nsize(p),7*M0))
  allocate(work2(M0,4))   

  ijob=-1
  do while (ijob/=0)
     call pzbicgstab_rci(ijob,Nsize(p),M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,MY_COMM_WORLD,info) 
!print *,'bicg',ijob
     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input

     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)
    
        call zhlbcsrmm(UPLO,TRANS,Nsize,map,M0,(1.0d0,0.0d0),matAjb,work(1,j1),(0.0d0,0.0d0),work(1,j2),MY_COMM_WORLD,nbprocs)
        !print *,'done'
!         print *,'workj2',sum(abs(work(1:Nsize(p),j2))),work(1000,j2)

     end select
  end do




  deallocate(work,work2)


end subroutine pzhbicgstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  pzbicgstab_rci(ijob,N,M0,fj,xj,work1,work2,j1,j2,res,nbitmax,epso,comd,MY_COMM_WORLD,info) 
  implicit none
    !-------------------------------------
  include 'mpif.h'
  !-------------------------------------

  integer :: j1,j2,MY_COMM_WORLD 
  integer :: ijob
  integer :: N,M0
  complex(kind=kind(1.0d0)), dimension(N,*) :: fj
  complex(kind=kind(1.0d0)), dimension(N,*):: xj
  complex(kind=kind(1.0d0)), dimension(N,*):: work1
  complex(kind=kind(1.0d0)), dimension(M0,*):: work2 

  double precision, dimension(M0) :: res


  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex(kind=kind(1.0d0)), parameter :: ONEC=(1.0d0,0.0d0),ZEROC=(0.0d0,0.0d0)

  double precision ::ares
  integer :: rank,i
  integer :: rj,rb,pp,pb,v,sb,t  !!! column index for the work1 array 
  integer:: rho_1,alpha,omega,ind !!! column index for the work2 array

  double precision,dimension(M0) :: resf
  complex(kind=kind(1.0d0)),dimension(M0) :: rho_2,aux0,aux1,beta
  integer,save :: nbit
  integer,save :: probe

  !complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
  double precision :: DZNRM2

  logical :: loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: code
  
  integer :: b1,b2


  b1=1
  b2=M0

!!!!!!!!!! column index for the work1 array
  rj=1
  rb=M0+1
  pp=2*M0+1
  pb=3*M0+1
  v=4*M0+1
  sb=5*M0+1
  t=6*M0+1
!!!!!!!!!! column index for the work2 array
  rho_1=1
  alpha=2
  omega=3
  ind=4 ! if /= (0,0) do not include residual calculations in the max(norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rank=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !    rho_1(M0)   ! global  >> work2(:,1)
  !    rho_2(M0))  ! local
  !    beta(M0)    ! local
  !    alpha(M0)   ! global >> work2(:,2)
  !    omega(M0)   ! global >> work2(:,3)
  !    aux0(M0)    ! local
  !    aux1(M0)    ! local   
  !    resf(M0)    ! local (use res to store it)

!!!!!!!!!!!!

  info=0
  loop=.true.   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ijob==-1) then
!!!!!!!!!!!!!!! initialize index for residual calculation
     if (M0==1) then
        work2(1,ind)=ONEC
     else
        do i=1,M0
           work2(i,ind)=res(i)*ONEC
        enddo
!!! special case (user enter 0 everywhere ==> 1 everywhere)
        if (sum(work2(:,ind))==ZEROC) work2(:,ind)=ONEC

     end if

!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax      
!!!! let us temporaly save xj into work1(1,rj)
     call ZCOPY(N*M0,xj(1,1),1,work1(1,rj),1)


!!!! Compute A*xj  ----- A*work1(1,rj)=>work1(1,rb)
     j1=rj
     j2=rb
     ijob=2
     probe=0
     !work1(1:N,rb:rb+M0-1)=ZEROC ! if initial guess 0
     return
  end if


  if (probe==0) then

!!!! let us temporaly save fj into work1(1,rj)
     call ZCOPY(N*M0,fj(1,1),1,work1(1,rj),1)

!!!!!!! compute residual r=f-A*x
     do i=1,M0   
        call ZAXPY(N,-ONEC,work1(1,i+rb-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!! CHECK THE  NORM of rj
     do i=1,M0
        res(i)=DZNRM2(N,work1(1,i+rj-1),1)**2
        resf(i)=DZNRM2(N,fj(1,i+rj-1),1)**2
     enddo
     call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     call MPI_ALLREDUCE(MPI_IN_PLACE,resf,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
        resf(i)=sqrt(resf(i)) ! norm L2
     enddo
    
     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))

     res=resf ! save resf value
     if (comd) print *,'rel. residual before iteration',ares,res(1),resf(1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! BIPCG-STAB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! choose rb
     call ZCOPY(N*M0,work1(1,rj),1,work1(1,rb),1)

     nbit=0

     probe=1

  endif


  if (probe==1) then
     nbit=nbit+1
     if (nbit>1) rho_2=work2(1:M0,rho_1) 

!!!!!!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0
        !work2(i,rho_1)=zdotc(N,work1(1,rb+i-1),1,work1(1,rj+i-1),1)
        work2(i,rho_1)=dot_product(work1(1:N,rb+i-1),work1(1:N,rj+i-1))
     end do
 call MPI_ALLREDUCE(MPI_IN_PLACE,work2(1,rho_1),M0,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)

! print *,'work2,rho_1',work2(1:3,rho_1)
      
!!!!!!!!!!!! TEST
     do i=1,M0
        if (work2(i,rho_1)==ZEROC) then
           info= -201
           if ((comd).and.(rank==0)) print*,'ATTENTION ----- BICG-STAB OUT FAILED, rho_1=0 !!'
           nbitmax=nbit 
           ijob=0
           return              
        end if
     end do

!!!!!!!!!!!! CONDITION
     IF (nbit==1) THEN
        call ZCOPY(N*M0,work1(1,rj),1,work1(1,pp),1)
     ELSE
        beta=(work2(1:M0,rho_1)/rho_2)*(work2(1:M0,alpha)/work2(1:M0,omega))
        do i=1,M0   !!!!!!<<< to optimize
           call ZAXPY(N,-work2(i,omega),work1(1,v+i-1),1,work1(1,pp+i-1),1)
           call ZSCAL(N,beta(i),work1(1,pp+i-1),1)
           call ZAXPY(N,ONEC,work1(1,i+rj-1),1,work1(1,i+pp-1),1)
           !pp(:,i)=r(:,i)+beta(i)*(pp(:,i)-omega(i)*v(:,i))
        end do
     END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve  M*pb=pp(k) ---> pb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call ZCOPY(N*M0,work1(1,pp),1,work1(1,pb),1)
     j1=pb
     j2=pp
     ijob=1
     probe=2
     return
  endif

  if (probe==2) then  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by pb, results in v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=pb
     j2=v
     ijob=2
     probe=3
     return
  end if

  if (probe==3) then
!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0!! different RHS
        !aux0(i)=zdotc(N,work1(1,i+rb-1),1,work1(1,i+v-1),1)
        aux0(i)=dot_product(work1(1:N,i+rb-1),work1(1:N,i+v-1))
        !aux0(i)=sum(work1(1:N,i+rb-1)*work1(1:N,i+v-1))
     end do
 call MPI_ALLREDUCE(MPI_IN_PLACE,aux0,M0,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)

!  print *,'work1pb',sum(abs(work1(1:N,pb)))
!      print *,'work1v',sum(abs(work1(1:N,v)))
!print *,'aux0',aux0(1:3)
   
 
!!!!!!
     work2(1:M0,alpha)=work2(1:M0,rho_1)/aux0(1:M0)
     do i=1,M0   
        !ss(:,i)=r(:,i)-alpha(i)*v(:,i) !!!!!!!!! use rj to store ss !!!!!!!!!!!
        call ZAXPY(N,-work2(i,alpha),work1(1,i+v-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!!! CHECK THE  NORM of ss
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DZNRM2(N,work1(1,i+rj-1),1)**2
     enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))



     IF (ares<epso) then
        do i=1,M0  
           !      xj(:,i)=xj(:,i)+alpha(i)*pb(:,i) 
           call ZAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        end do
        if ((comd).and.(rank==0)) print *,(nbit-1)*1.0d0+0.5d0,ares
        info=0
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares             
        return
     end IF

     res=resf ! save resf value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve M*sb=ss ---> sb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     call ZCOPY(N*M0,work1(1,rj),1,work1(1,sb),1)
     j1=sb
     j2=rj
     ijob=1
     probe=4
     return
  end if


  if (probe==4) then          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by sb, results in t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=sb
     j2=t
     ijob=2
     probe=5
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! SCALAR PRODUCTS on different RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (probe==5) then

     do i=1,M0!! different RHS
        !aux0(i)=zdotc(N,work1(1,i+t-1),1,work1(1,i+rj-1),1)
        !aux1(i)=zdotc(N,work1(1,i+t-1),1,work1(1,i+t-1),1)
         aux0(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+rj-1))
        aux1(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+t-1))
     end do
     call MPI_ALLREDUCE(MPI_IN_PLACE,aux0,M0,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
     call MPI_ALLREDUCE(MPI_IN_PLACE,aux1,M0,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)

!!!!!!
     work2(1:M0,omega)=aux0/aux1
     do i=1,M0
        if (work2(i,omega)==ZEROC) then
           info= -202 
           if ((comd).and.(rank==0)) print *,'ATTENTION ----- BICG-STAB OUT FAILED, omega=0'
           nbitmax=nbit
           ijob=0
           return 
           !fail=.true.
        end if
     end do

     do i=1,M0
        ! xj(:,i)=xj(:,i)+alpha(i)*pb(:,i)+omega(i)*sb(:,i)
        call ZAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        call ZAXPY(N,work2(i,omega),work1(1,i+sb-1),1,xj(1,i),1) 
     end do

     do i=1,M0
        !   r(:,i)=ss(:,i)-omega(i)*t(:,i)
        !call ZCOPY(N,ss(1,i),1,rj(1,i),1)   !!! ss is rj
        call ZAXPY(N,-work2(i,omega),work1(1,i+t-1),1,work1(1,i+rj-1),1) 
     end do


!!!!!!!!! CHECK THE  NORM of rj
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=DZNRM2(N,work1(1,i+rj-1),1)**2
     enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
     enddo



     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))
     if ((comd).and.(rank==0)) print *,nbit,ares


!!!! test for the other loop
     if (.not.((nbit<nbitmax).and.(ares.gt.epso))) loop=.false. 

     if (loop) then
        ijob=-2
        probe=1
        res=resf ! save resf 
     else
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares      
     end if
     return

  end if


end subroutine pzbicgstab_rci






subroutine  pcbicgstab(opt,uplo,trans,Nsize,map,M0,matAjb,fj,xj,res,nbitmax,epso,comd,MY_COMM_WORLD,nbprocs,info)
  implicit none
  !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  integer :: opt, MY_COMM_WORLD, M0, nbprocs
  character(len=1) :: uplo,trans
 integer,dimension(*) :: Nsize
  integer,dimension(nbprocs,*) :: map
   type ccsr
     integer ::n,m,nnz
     complex,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type ccsr
  type(ccsr),dimension(nbprocs,*) :: matAjb
  complex, dimension(*) :: fj
  complex, dimension(*):: xj
  double precision, dimension(M0) :: res
  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
integer :: code,rank,p
complex,dimension(:,:),allocatable :: work,work2
integer :: ijob,j1,j2,k,i
 
!!!!!!!
  call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
  p=rank+1
!!!!
  
  allocate(work(Nsize(p),7*M0))
  allocate(work2(M0,4))   

  ijob=-1
  do while (ijob/=0)
     call pcbicgstab_rci(ijob,Nsize(p),M0,fj,xj,work,work2,j1,j2,res,nbitmax,epso,comd,MY_COMM_WORLD,info) 

     select case(ijob)
        !         case(1) !!solve M0 rhs with preconditioner if any M*work(:,j1:j1+M0-1)=work(:,j2:j2+M0-1) result in work(:,j1:)
        !! j2 can be used or not,  since  work(:,j2:j2+M0-1)=work(:,j1:j1+M0-1) as input

     case(2) !! mat-vec M0 rhs      work(:,j2:)<=A*work(:,j1:)
       
        call clbcsrmm(UPLO,TRANS,Nsize,map,M0,(1.0e0,0.0e0),matAjb,work(1,j1),(0.0e0,0.0e0),work(1,j2),MY_COMM_WORLD,nbprocs)
       
     end select
  end do



  deallocate(work,work2)


end subroutine pcbicgstab



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  pcbicgstab_rci(ijob,N,M0,fj,xj,work1,work2,j1,j2,res,nbitmax,epso,comd,MY_COMM_WORLD,info) 
  implicit none
    !-------------------------------------
  include 'mpif.h'
  !-------------------------------------

  integer :: j1,j2,MY_COMM_WORLD 
  integer :: ijob
  integer :: N,M0
  complex, dimension(N,*) :: fj
  complex, dimension(N,*):: xj
  complex, dimension(N,*):: work1
  complex, dimension(M0,*):: work2 

  double precision, dimension(M0) :: res


  integer :: info
  integer :: nbitmax
  double precision ::epso
  logical :: comd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  complex, parameter :: ONEC=(1.0d0,0.0d0),ZEROC=(0.0d0,0.0d0)

  double precision ::ares
  integer :: rank,i
  integer :: rj,rb,pp,pb,v,sb,t  !!! column index for the work1 array 
  integer:: rho_1,alpha,omega,ind !!! column index for the work2 array

  double precision,dimension(M0) :: resf
  complex,dimension(M0) :: rho_2,aux0,aux1,beta
  integer,save :: nbit
  integer,save :: probe

  !complex :: cdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
  real :: SCNRM2

  logical :: loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: code
  
  integer :: b1,b2


  b1=1
  b2=M0

!!!!!!!!!! column index for the work1 array
  rj=1
  rb=M0+1
  pp=2*M0+1
  pb=3*M0+1
  v=4*M0+1
  sb=5*M0+1
  t=6*M0+1
!!!!!!!!!! column index for the work2 array
  rho_1=1
  alpha=2
  omega=3
  ind=4 ! if /= (0,0) do not include residual calculations in the max(norm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rank=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !    rho_1(M0)   ! global  >> work2(:,1)
  !    rho_2(M0))  ! local
  !    beta(M0)    ! local
  !    alpha(M0)   ! global >> work2(:,2)
  !    omega(M0)   ! global >> work2(:,3)
  !    aux0(M0)    ! local
  !    aux1(M0)    ! local   
  !    resf(M0)    ! local (use res to store it)

!!!!!!!!!!!!

  info=0
  loop=.true.   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Initial rj,xj,res and resf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ijob==-1) then
!!!!!!!!!!!!!!! initialize index for residual calculation
     if (M0==1) then
        work2(1,ind)=ONEC
     else
        do i=1,M0
           work2(i,ind)=res(i)*ONEC
        enddo
!!! special case (user enter 0 everywhere ==> 1 everywhere)
        if (sum(work2(:,ind))==ZEROC) work2(:,ind)=ONEC

     end if

!!!!!!!!!!!!!!! MATVEC for the initial residual r=f-Ax      
!!!! let us temporaly save xj into work1(1,rj)
     call CCOPY(N*M0,xj(1,1),1,work1(1,rj),1)


!!!! Compute A*xj  ----- A*work1(1,rj)=>work1(1,rb)
     j1=rj
     j2=rb
     ijob=2
     probe=0
     !work1(1:N,rb:rb+M0-1)=ZEROC ! if initial guess 0
     return
  end if


  if (probe==0) then

!!!! let us temporaly save fj into work1(1,rj)
     call CCOPY(N*M0,fj(1,1),1,work1(1,rj),1)

!!!!!!! compute residual r=f-A*x
     do i=1,M0   
        call CAXPY(N,-ONEC,work1(1,i+rb-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!! CHECK THE  NORM of rj
     do i=1,M0
        res(i)=SCNRM2(N,work1(1,i+rj-1),1)**2
        resf(i)=SCNRM2(N,fj(1,i+rj-1),1)**2
     enddo
     call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     call MPI_ALLREDUCE(MPI_IN_PLACE,resf,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
        resf(i)=sqrt(resf(i)) ! norm L2
     enddo
    
     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))

     res=resf ! save resf value
     if (comd) print *,'rel. residual before iteration',ares,res(1),resf(1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! BIPCG-STAB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!! choose rb
     call CCOPY(N*M0,work1(1,rj),1,work1(1,rb),1)

     nbit=0

     probe=1

  endif


  if (probe==1) then
     nbit=nbit+1
     if (nbit>1) rho_2=work2(1:M0,rho_1) 

!!!!!!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0
        !work2(i,rho_1)=cdotc(N,work1(1,rb+i-1),1,work1(1,rj+i-1),1)
        work2(i,rho_1)=dot_product(work1(1:N,rb+i-1),work1(1:N,rj+i-1))
     end do
 call MPI_ALLREDUCE(MPI_IN_PLACE,work2(1,rho_1),M0,MPI_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)

! print *,'work2,rho_1',work2(1:3,rho_1)
      
!!!!!!!!!!!! TEST
     do i=1,M0
        if (work2(i,rho_1)==ZEROC) then
           info= -201
           if ((comd).and.(rank==0)) print*,'ATTENTION ----- BICG-STAB OUT FAILED, rho_1=0 !!'
           nbitmax=nbit 
           ijob=0
           return              
        end if
     end do

!!!!!!!!!!!! CONDITION
     IF (nbit==1) THEN
        call CCOPY(N*M0,work1(1,rj),1,work1(1,pp),1)
     ELSE
        beta=(work2(1:M0,rho_1)/rho_2)*(work2(1:M0,alpha)/work2(1:M0,omega))
        do i=1,M0   !!!!!!<<< to optimize
           call CAXPY(N,-work2(i,omega),work1(1,v+i-1),1,work1(1,pp+i-1),1)
           call CSCAL(N,beta(i),work1(1,pp+i-1),1)
           call CAXPY(N,ONEC,work1(1,i+rj-1),1,work1(1,i+pp-1),1)
           !pp(:,i)=r(:,i)+beta(i)*(pp(:,i)-omega(i)*v(:,i))
        end do
     END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve  M*pb=pp(k) ---> pb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call CCOPY(N*M0,work1(1,pp),1,work1(1,pb),1)
     j1=pb
     j2=pp
     ijob=1
     probe=2
     return
  endif

  if (probe==2) then  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by pb, results in v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=pb
     j2=v
     ijob=2
     probe=3
     return
  end if

  if (probe==3) then
!!!!!! SCALAR PRODUCT on different RHS
     do i=1,M0!! different RHS
        !aux0(i)=cdotc(N,work1(1,i+rb-1),1,work1(1,i+v-1),1)
        aux0(i)=dot_product(work1(1:N,i+rb-1),work1(1:N,i+v-1))
        !aux0(i)=sum(work1(1:N,i+rb-1)*work1(1:N,i+v-1))
     end do
 call MPI_ALLREDUCE(MPI_IN_PLACE,aux0,M0,MPI_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)

!  print *,'work1pb',sum(abs(work1(1:N,pb)))
!      print *,'work1v',sum(abs(work1(1:N,v)))
!print *,'aux0',aux0(1:3)
   
 
!!!!!!
     work2(1:M0,alpha)=work2(1:M0,rho_1)/aux0(1:M0)
     do i=1,M0   
        !ss(:,i)=r(:,i)-alpha(i)*v(:,i) !!!!!!!!! use rj to store ss !!!!!!!!!!!
        call CAXPY(N,-work2(i,alpha),work1(1,i+v-1),1,work1(1,i+rj-1),1) 
     end do

!!!!!!!!!! CHECK THE  NORM of ss
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=SCNRM2(N,work1(1,i+rj-1),1)**2
     enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
     enddo


     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))



     IF (ares<epso) then
        do i=1,M0  
           !      xj(:,i)=xj(:,i)+alpha(i)*pb(:,i) 
           call CAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        end do
        if ((comd).and.(rank==0)) print *,(nbit-1)*1.0d0+0.5d0,ares
        info=0
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares             
        return
     end IF

     res=resf ! save resf value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! SOLVE  solve M*sb=ss ---> sb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
     call CCOPY(N*M0,work1(1,rj),1,work1(1,sb),1)
     j1=sb
     j2=rj
     ijob=1
     probe=4
     return
  end if


  if (probe==4) then          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAT-VEC A by sb, results in t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     j1=sb
     j2=t
     ijob=2
     probe=5
     return
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! SCALAR PRODUCTS on different RHS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (probe==5) then

     do i=1,M0!! different RHS
       ! aux0(i)=cdotc(N,work1(1,i+t-1),1,work1(1,i+rj-1),1)
       ! aux1(i)=cdotc(N,work1(1,i+t-1),1,work1(1,i+t-1),1)
         aux0(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+rj-1))
        aux1(i)=dot_product(work1(1:N,i+t-1),work1(1:N,i+t-1))
     end do
     call MPI_ALLREDUCE(MPI_IN_PLACE,aux0,M0,MPI_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
     call MPI_ALLREDUCE(MPI_IN_PLACE,aux1,M0,MPI_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)

!!!!!!
     work2(1:M0,omega)=aux0/aux1
     do i=1,M0
        if (work2(i,omega)==ZEROC) then
           info= -202 
           if ((comd).and.(rank==0)) print *,'ATTENTION ----- BICG-STAB OUT FAILED, omega=0'
           nbitmax=nbit
           ijob=0
           return 
           !fail=.true.
        end if
     end do

     do i=1,M0
        ! xj(:,i)=xj(:,i)+alpha(i)*pb(:,i)+omega(i)*sb(:,i)
        call CAXPY(N,work2(i,alpha),work1(1,i+pb-1),1,xj(1,i),1) 
        call CAXPY(N,work2(i,omega),work1(1,i+sb-1),1,xj(1,i),1) 
     end do

     do i=1,M0
        !   r(:,i)=ss(:,i)-omega(i)*t(:,i)
        !call ZCOPY(N,ss(1,i),1,rj(1,i),1)   !!! ss is rj
        call CAXPY(N,-work2(i,omega),work1(1,i+t-1),1,work1(1,i+rj-1),1) 
     end do


!!!!!!!!! CHECK THE  NORM of rj
     resf=res ! retrieve resf value
     do i=1,M0
        res(i)=SCNRM2(N,work1(1,i+rj-1),1)**2
     enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,res,M0,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
     do i=1,M0
        res(i)=sqrt(res(i)) ! norm L2
     enddo



     if (dble(sum(work2(:,ind)))<0.0d0) then !look for the min
        ares=1.0d0
        do i=1,M0
           if ((res(i)/resf(i)<ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     else  if (dble(sum(work2(:,ind)))>0.0d0) then!look for the max
        ares=-1.0d0
        do i=1,M0
           if ((res(i)/resf(i)>ares).and.work2(i,ind)/=ZEROC) ares=res(i)/resf(i)
        enddo
     end if
     !ares=maxval(res(b1:b2)/resf(b1:b2))
     if ((comd).and.(rank==0)) print *,nbit,ares


!!!! test for the other loop
     if (.not.((nbit<nbitmax).and.(ares.gt.epso))) loop=.false. 

     if (loop) then
        ijob=-2
        probe=1
        res=resf ! save resf 
     else
        ijob=0
        nbitmax=nbit
        res=res/resf
        epso=ares      
     end if
     return

  end if


end subroutine pcbicgstab_rci








subroutine  pdarnoldi(UPLO,Nsize,mapA,matAjb,alpha,X,res,epse,itmax,MY_COMM_WORLD,nb_procs)
  !! Purpose: Lanczos algorithm 


  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
   !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: MY_COMM_WORLD,nb_procs
  double precision,dimension(*) :: alpha,res
   integer,dimension(nb_procs,*) :: mapA
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr), dimension(nb_procs,*) :: matAjb  
  double precision,dimension(*) :: X
  integer :: itmax
  integer,dimension(*) :: Nsize
  double precision :: epse
!!!
  integer,dimension(4) :: iseed
   double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   double precision,dimension(:),allocatable :: work,uloc,beta,de,rese,uglo
    double precision,dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye
    integer :: i,j,lwork,info_lap,it
    integer:: p,rank,code,distype,Ntotal,Nlocal,startj,endj
    double precision :: DNRM2,ddot,temp
     logical :: conv
  integer :: mine,maxe


!!!!  
call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
p=rank+1
!!!! 
 
Nlocal=Nsize(p)
Ntotal=sum(Nsize(1:nb_procs))
 !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

allocate(de(2))
  allocate(rese(2))
  allocate(Xe(Nlocal,2))
  allocate(Ye(Nlocal,2))
 
  allocate(uloc(Nlocal))
  allocate(uglo(Ntotal)) 
   allocate(V(Nlocal,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !v1=random - use Notal for consistency in getting emin-emax for differnt L3 processors

 iseed=(/45,8,89,33/)
 call DLARNV(1,iseed , Ntotal, uglo ) ! random 0-1:  uniform
 uglo(:)=uglo(:)/DNRM2(Ntotal,uglo,1) ! normalization
 v(1:Nlocal,1)=uglo(startj:startj+Nlocal-1) ! distribute
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Arnoldi direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    H=DZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       call dlbcsrmm(UPLO,'N',Nsize,mapA,1,DONE,matAjb,v(1,it),DZERO,uloc,MY_COMM_WORLD,nb_procs)
       do j=1,it
          H(j,it)=ddot(Nlocal,V(1,j),1,uloc,1)
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(j,it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
          call DAXPY(Nlocal,-H(j,it),V(1,j),1,uloc,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DNRM2(Nlocal,uloc,1)**2
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(it+1,it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
          H(it+1,it)=sqrt(H(it+1,it))
          call DCOPY(Nlocal,uloc,1,V(1,it+1),1)
          call DSCAL(Nlocal,DONE/H(it+1,it), v(1,it+1), 1)
       end if

         if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call DHSEQR('E', 'I', it, 1, it, H_temp, itmax, alpha, beta, Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (alpha(j)>alpha(maxe)) maxe=j
           if (alpha(j)<alpha(mine)) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call DGEMM('N', 'N', Nlocal, 1, it, DONE, V(1,1), Nlocal, Xred(1,mine), itmax, DZERO, Xe(1,1), Nlocal)
        call DGEMM('N', 'N', Nlocal, 1, it, DONE, V(1,1), Nlocal, Xred(1,maxe), itmax, DZERO, Xe(1,2), Nlocal)
        dE(1)=alpha(mine)
        dE(2)=alpha(maxe)

        call DLACPY('F',Nlocal,2,Xe,Nlocal,Ye,Nlocal)
        do j=1,2
           call DSCAL(Nlocal,dE(j), Ye(1,j), 1)
           rese(j)=DNRM2(Nlocal,Ye(1,j),1)**2 ! placeholder denominateur
        enddo
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  call dlbcsrmm(UPLO,'N',Nsize,mapA,2,DONE,matAjb,Xe,-DONE,Ye,MY_COMM_WORLD,nb_procs)
 do j=1,2
    rese(j)=(DNRM2(Nlocal,Ye(1,j),1)**2)/rese(j)
 enddo
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  do j=1,2
     rese(j)=sqrt(rese(j)) ! norm L2
     !   if (j==2)
     !print *,it,j,de(j),rese(j)
     enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do

itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)<alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo


!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
  call DGEMM('N', 'N', Nlocal, itmax, itmax, DONE, V(1,1), Nlocal, Xred, itmax, DZERO, X, Nlocal)
 


  deallocate(uglo)
  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
      
   deallocate(uloc) 
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine pdarnoldi







subroutine  pzarnoldi(UPLO,Nsize,mapA,matAjb,alpha,X,res,epse,itmax,MY_COMM_WORLD,nb_procs)
!! check and sort amplitude of eigenvalues....

  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
   !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: MY_COMM_WORLD,nb_procs
  complex(kind=kind(1.0d0)),dimension(*) :: alpha
  double precision,dimension(*) :: res
   integer,dimension(nb_procs,*) :: mapA
  type zcsr
     integer ::n,m,nnz
      complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr), dimension(nb_procs,*) :: matAjb  
   complex(kind=kind(1.0d0)),dimension(*) :: X
  integer :: itmax
  integer,dimension(*) :: Nsize
  double precision :: epse
!!!
  integer,dimension(4) :: iseed
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   complex(kind=kind(1.0d0)),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
   complex(kind=kind(1.0d0)),dimension(:),allocatable :: work,uloc,beta,de,uglo
    double precision,dimension(:),allocatable :: rese
     complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye
    integer :: i,j,lwork,info_lap,it
    integer:: p,rank,code,distype,Ntotal,Nlocal,startj,endj
    double precision :: DZNRM2,temp
   !   complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
     logical :: conv
  integer :: mine,maxe


!!!!  
call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
p=rank+1
!!!! 
 
Nlocal=Nsize(p)
Ntotal=sum(Nsize(1:nb_procs))
 !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

allocate(de(2))
  allocate(rese(2))
  allocate(Xe(Nlocal,2))
  allocate(Ye(Nlocal,2))
 
  allocate(uloc(Nlocal))
  allocate(uglo(Ntotal)) 
   allocate(V(Nlocal,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !v1=random - use Notal for consistency in getting emin-emax for differnt L3 processors

 iseed=(/45,8,89,33/)
 call ZLARNV(1,iseed , Ntotal, uglo ) ! random 0-1:  uniform
 uglo(:)=uglo(:)/DZNRM2(Ntotal,uglo,1) ! normalization
 v(1:Nlocal,1)=uglo(startj:startj+Nlocal-1) ! distribute
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Arnoldi direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    H=ZZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       call zlbcsrmm(UPLO,'N',Nsize,mapA,1,ZONE,matAjb,v(1,it),ZZERO,uloc,MY_COMM_WORLD,nb_procs)
       do j=1,it
          !H(j,it)=zdotc(Nlocal,V(1,j),1,uloc,1)
           H(j,it)=dot_product(V(1:Nlocal,j),uloc(1:Nlocal))
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(j,it),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
          call ZAXPY(Nlocal,-H(j,it),V(1,j),1,uloc,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DZNRM2(Nlocal,uloc,1)**2
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(it+1,it),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
          H(it+1,it)=sqrt(H(it+1,it))
          call ZCOPY(Nlocal,uloc,1,V(1,it+1),1)
          call ZSCAL(Nlocal,ZONE/H(it+1,it), v(1,it+1), 1)
       end if

         if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call ZHSEQR('E', 'I', it, 1, it, H_temp, itmax, alpha,  Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (abs(alpha(j))> abs(alpha(maxe))) maxe=j
           if (abs(alpha(j))< abs(alpha(mine))) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call ZGEMM('N', 'N', Nlocal, 1, it, ZONE, V(1,1), Nlocal, Xred(1,mine), itmax, ZZERO, Xe(1,1), Nlocal)
        call ZGEMM('N', 'N', Nlocal, 1, it, ZONE, V(1,1), Nlocal, Xred(1,maxe), itmax, ZZERO, Xe(1,2), Nlocal)
        dE(1)=alpha(mine)
        dE(2)=alpha(maxe)

        call ZLACPY('F',Nlocal,2,Xe,Nlocal,Ye,Nlocal)
        do j=1,2
           call ZSCAL(Nlocal,dE(j), Ye(1,j), 1)
           rese(j)=DZNRM2(Nlocal,Ye(1,j),1)**2 ! placeholder denominateur
        enddo
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  call zlbcsrmm(UPLO,'N',Nsize,mapA,2,ZONE,matAjb,Xe,-ZONE,Ye,MY_COMM_WORLD,nb_procs)
 do j=1,2
    rese(j)=(DZNRM2(Nlocal,Ye(1,j),1)**2)/rese(j)
 enddo
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  do j=1,2
     rese(j)=sqrt(rese(j)) ! norm L2
     !   if (j==2)
     !print *,it,j,de(j),rese(j)
     enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do

itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
do i=1,itmax
   do j=i+1,itmax
      if (abs(alpha(j))< abs(alpha(i))) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo


!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
  call ZGEMM('N', 'N', Nlocal, itmax, itmax, ZONE, V(1,1), Nlocal, Xred, itmax, ZZERO, X, Nlocal)
 


  deallocate(uglo)
  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
      
   deallocate(uloc) 
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine pzarnoldi







subroutine  pzharnoldi(UPLO,Nsize,mapA,matAjb,alpha,X,res,epse,itmax,MY_COMM_WORLD,nb_procs)
!! check and sort real part of eigenvalues....real eigenvalue here

  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
   !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: MY_COMM_WORLD,nb_procs
  double precision,dimension(*) :: res,alpha
   integer,dimension(nb_procs,*) :: mapA
  type zcsr
     integer ::n,m,nnz
      complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr), dimension(nb_procs,*) :: matAjb  
   complex(kind=kind(1.0d0)),dimension(*) :: X
  integer :: itmax
  integer,dimension(*) :: Nsize
  double precision :: epse
!!!
  integer,dimension(4) :: iseed
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   complex(kind=kind(1.0d0)),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
   complex(kind=kind(1.0d0)),dimension(:),allocatable :: work,uloc,beta,dE,uglo,zalpha
    double precision,dimension(:),allocatable :: rese
     complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye
    integer :: i,j,lwork,info_lap,it
    integer:: p,rank,code,distype,Ntotal,Nlocal,startj,endj
    double precision :: DZNRM2,temp
!      complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
     logical :: conv
  integer :: mine,maxe


!!!!  
call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
p=rank+1
!!!! 
 
Nlocal=Nsize(p)
Ntotal=sum(Nsize(1:nb_procs))
 !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

allocate(dE(2))
  allocate(rese(2))
  allocate(Xe(Nlocal,2))
  allocate(Ye(Nlocal,2))
  allocate(zalpha(itmax))
  
  allocate(uloc(Nlocal))
  allocate(uglo(Ntotal)) 
   allocate(V(Nlocal,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !v1=random - use Notal for consistency in getting emin-emax for differnt L3 processors

 iseed=(/45,8,89,33/)
 call ZLARNV(1,iseed , Ntotal, uglo ) ! random 0-1:  uniform
 uglo(:)=uglo(:)/DZNRM2(Ntotal,uglo,1) ! normalization
 v(1:Nlocal,1)=uglo(startj:startj+Nlocal-1) ! distribute
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Arnoldi direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    H=ZZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       call zhlbcsrmm(UPLO,'N',Nsize,mapA,1,ZONE,matAjb,v(1,it),ZZERO,uloc,MY_COMM_WORLD,nb_procs)
       do j=1,it
         ! H(j,it)=zdotc(Nlocal,V(1,j),1,uloc,1)
           H(j,it)=dot_product(V(1:Nlocal,j),uloc(1:Nlocal))
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(j,it),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
          call ZAXPY(Nlocal,-H(j,it),V(1,j),1,uloc,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DZNRM2(Nlocal,uloc,1)**2
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(it+1,it),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
          H(it+1,it)=sqrt(H(it+1,it))
          call ZCOPY(Nlocal,uloc,1,V(1,it+1),1)
          call ZSCAL(Nlocal,ZONE/H(it+1,it), v(1,it+1), 1)
       end if

         if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call ZHSEQR('E', 'I', it, 1, it, H_temp, itmax, zalpha,  Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (dble(zalpha(j))> dble(zalpha(maxe))) maxe=j
           if (dble(zalpha(j))< dble(zalpha(mine))) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call ZGEMM('N', 'N', Nlocal, 1, it, ZONE, V(1,1), Nlocal, Xred(1,mine), itmax, ZZERO, Xe(1,1), Nlocal)
        call ZGEMM('N', 'N', Nlocal, 1, it, ZONE, V(1,1), Nlocal, Xred(1,maxe), itmax, ZZERO, Xe(1,2), Nlocal)
        dE(1)=zalpha(mine)
        dE(2)=zalpha(maxe)

        call ZLACPY('F',Nlocal,2,Xe,Nlocal,Ye,Nlocal)
        do j=1,2
           call ZSCAL(Nlocal,dE(j), Ye(1,j), 1)
           rese(j)=DZNRM2(Nlocal,Ye(1,j),1)**2 ! placeholder denominateur
        enddo
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  call zhlbcsrmm(UPLO,'N',Nsize,mapA,2,ZONE,matAjb,Xe,-ZONE,Ye,MY_COMM_WORLD,nb_procs)
 do j=1,2
    rese(j)=(DZNRM2(Nlocal,Ye(1,j),1)**2)/rese(j)
 enddo
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  do j=1,2
     rese(j)=sqrt(rese(j)) ! norm L2
     !   if (j==2)
     !print *,it,j,de(j),rese(j)
     enddo
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do

itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
alpha(1:itmax)=dble(zalpha(1:itmax))
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)< alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo


!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
  call ZGEMM('N', 'N', Nlocal, itmax, itmax, ZONE, V(1,1), Nlocal, Xred, itmax, ZZERO, X, Nlocal)
 


  deallocate(uglo)
  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
  deallocate(zalpha)

  
   deallocate(uloc) 
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine pzharnoldi






subroutine  pdgarnoldi(UPLO,Nsize,mapA,matAjb,mapB,matBjb,alpha,X,res,epse,itmax,epsb,itmaxb,MY_COMM_WORLD,nb_procs)
  !! Purpose: Lanczos algorithm 


  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
   !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: MY_COMM_WORLD,nb_procs,itmaxb
  double precision,dimension(*) :: alpha,res
   integer,dimension(nb_procs,*) :: mapA,mapB
  type dcsr
     integer ::n,m,nnz
     double precision,dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type dcsr
  type(dcsr), dimension(nb_procs,*) :: matAjb,matBjb  
  double precision,dimension(*) :: X
  integer :: itmax
  integer,dimension(*) :: Nsize
  double precision :: epse,epsb
!!!
  integer,dimension(4) :: iseed
   double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   double precision,dimension(:),allocatable :: work,uloc,beta,de,rese,uglo,w
    double precision,dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye,aux
    integer :: i,j,lwork,info_lap,it,linloops,infoloc
    integer:: p,rank,code,distype,Ntotal,Nlocal,startj,endj
    double precision :: DNRM2,ddot,temp,eps
     logical :: conv,comb
  integer :: mine,maxe


!!!!  
call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
p=rank+1
!!!! 
Nlocal=Nsize(p)
Ntotal=sum(Nsize(1:nb_procs))
 !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

allocate(de(2))
  allocate(rese(2))
  allocate(Xe(Nlocal,2))
  allocate(Ye(Nlocal,2))
    allocate(aux(Nlocal,2))
 
  allocate(uloc(Nlocal))
   allocate(w(Nlocal))
  allocate(uglo(Ntotal)) 
   allocate(V(Nlocal,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !v1=random - use Notal for consistency in getting emin-emax for differnt L3 processors

 iseed=(/45,8,89,33/)
 call DLARNV(1,iseed , Ntotal, uglo ) ! random 0-1:  uniform
 uglo(:)=uglo(:)/DNRM2(Ntotal,uglo,1) ! normalization
 v(1:Nlocal,1)=uglo(startj:startj+Nlocal-1) ! distribute

 comb=.false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Arnoldi direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    H=DZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       !mat-vec
       call dlbcsrmm(UPLO,'N',Nsize,mapA,1,DONE,matAjb,v(1,it),DZERO,w,MY_COMM_WORLD,nb_procs)
       !solve
     linloops=itmaxb
     eps=epsb
     uloc(1:Nlocal)=DZERO !initial guess
     !call dbicgstab(1,UPLO,'N',sb,isb,jsb,N,1,w,u,res,linloops,eps,comb,infoloc)
     call pdbicgstab(0,uplo,'N',Nsize,mapB,1,matBjb,w,uloc,res,linloops,eps,comb,MY_COMM_WORLD,nb_procs,infoloc)
       !
       do j=1,it
          H(j,it)=ddot(Nlocal,V(1,j),1,uloc,1)
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(j,it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
          call DAXPY(Nlocal,-H(j,it),V(1,j),1,uloc,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DNRM2(Nlocal,uloc,1)**2
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(it+1,it),1,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
          H(it+1,it)=sqrt(H(it+1,it))
          call DCOPY(Nlocal,uloc,1,V(1,it+1),1)
          call DSCAL(Nlocal,DONE/H(it+1,it), v(1,it+1), 1)
       end if

         if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call DHSEQR('E', 'I', it, 1, it, H_temp, itmax, alpha, beta, Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (alpha(j)>alpha(maxe)) maxe=j
           if (alpha(j)<alpha(mine)) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call DGEMM('N', 'N', Nlocal, 1, it, DONE, V(1,1), Nlocal, Xred(1,mine), itmax, DZERO, Xe(1,1), Nlocal)
        call DGEMM('N', 'N', Nlocal, 1, it, DONE, V(1,1), Nlocal, Xred(1,maxe), itmax, DZERO, Xe(1,2), Nlocal)
        dE(1)=alpha(mine)
        dE(2)=alpha(maxe)

        call DLACPY('F',Nlocal,2,Xe,Nlocal,Ye,Nlocal)
        do j=1,2
           call DSCAL(Nlocal,dE(j), Ye(1,j), 1)
           rese(j)=DNRM2(Nlocal,Ye(1,j),1)**2 ! placeholder denominateur
        enddo
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
 call dlbcsrmm(UPLO,'N',Nsize,mapA,2,DONE,matAjb,Xe,DZERO,aux,MY_COMM_WORLD,nb_procs)
 
 linloops=itmaxb
     eps=epsb
     Xe(1:Nlocal,1:2)=DZERO !initial guess
     call pdbicgstab(0,uplo,'N',Nsize,mapB,2,matBjb,aux,Xe,res,linloops,eps,comb,MY_COMM_WORLD,nb_procs,infoloc)
     
  do j=1,2
           call DAXPY(Nlocal,-DONE,Xe(1,j),1,Ye(1,j),1)
           rese(j)=(DNRM2(Nlocal,Ye(1,j),1)**2)/rese(j)
        enddo
        
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  do j=1,2
     rese(j)=sqrt(rese(j)) ! norm L2
     !   if (j==2)
     !print *,it,j,de(j),rese(j)
  enddo
  
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do
! itmax=50
!           return
  
itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)<alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo


!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
  call DGEMM('N', 'N', Nlocal, itmax, itmax, DONE, V(1,1), Nlocal, Xred, itmax, DZERO, X, Nlocal)
 


  deallocate(uglo)
  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
  deallocate(aux)
      
  deallocate(uloc)
    deallocate(w)
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine pdgarnoldi







subroutine  pzhgarnoldi(UPLO,Nsize,mapA,matAjb,mapB,matBjb,alpha,X,res,epse,itmax,epsb,itmaxb,MY_COMM_WORLD,nb_procs)
!! check and sort real part of eigenvalues....real eigenvalue here


  !           using modified GS in-place 
  ! N    (input) INTEGER  #rows of the matrix Q
  ! M    (input) INTEGER  #columns of the matrix Q
  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
  !                     on entry: matrix to orthonormalize
  !                     on exit: orthonormal matrix
  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
  !====================================================================
  ! Eric Polizzi 2019
  !====================================================================  
  implicit none
   !-------------------------------------
  include 'mpif.h'
  !-------------------------------------
  character(len=1) :: UPLO
  integer :: MY_COMM_WORLD,nb_procs,itmaxb
  double precision,dimension(*) :: alpha,res
   integer,dimension(nb_procs,*) :: mapA,mapB
  type zcsr
     integer ::n,m,nnz
      complex(kind=kind(1.0d0)),dimension(:),allocatable :: sa
     integer,dimension(:),allocatable :: isa,jsa     
  end type zcsr
  type(zcsr), dimension(nb_procs,*) :: matAjb,matBjb  
   complex(kind=kind(1.0d0)),dimension(*) :: X
  integer :: itmax
  integer,dimension(*) :: Nsize
  double precision :: epse,epsb
!!!
  integer,dimension(4) :: iseed
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
   complex(kind=kind(1.0d0)),parameter :: ZONE=(DONE,DZERO), ZZERO=(DZERO,DZERO)
   complex(kind=kind(1.0d0)),dimension(:),allocatable :: work,uloc,beta,dE,uglo,w,zalpha
   double precision,dimension(:),allocatable :: rese
   complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: V,H,H_temp,Xred,Xe,Ye,aux
    integer :: i,j,lwork,info_lap,it,linloops,infoloc
    integer:: p,rank,code,distype,Ntotal,Nlocal,startj,endj
    double precision :: DZNRM2,temp,eps
 !    complex(kind=kind(1.0d0)) :: zdotc ! MKL problem with pgf90 and other hybrid combination gfortran+icc
     logical :: conv,comb
  integer :: mine,maxe


!!!!  
call MPI_COMM_RANK(MY_COMM_WORLD,rank,code)
p=rank+1
!!!! 
Nlocal=Nsize(p)
Ntotal=sum(Nsize(1:nb_procs))
 !! get startj,endj
     startj=1  
     endj=Nsize(1)
     do i=1,rank
        startj=startj+Nsize(i)
        endj=endj+Nsize(i+1)
     enddo

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

allocate(de(2))
  allocate(rese(2))
  allocate(Xe(Nlocal,2))
  allocate(Ye(Nlocal,2))
  allocate(aux(Nlocal,2))
    allocate(zalpha(itmax))
 
  allocate(uloc(Nlocal))
   allocate(w(Nlocal))
  allocate(uglo(Ntotal)) 
   allocate(V(Nlocal,itmax)) ! subspace
   allocate(H(itmax,itmax)) ! hessenberg matrix
    allocate(H_temp(itmax,itmax)) ! hessenberg matrix
   lwork=6*itmax!itmax!3*itmax-1
   allocate(work(lwork))
   allocate(Xred(itmax,itmax)) ! Ritz vector- reduced system
   allocate(beta(itmax)) ! imaginary part of vectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!! initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !v1=random - use Notal for consistency in getting emin-emax for differnt L3 processors

 iseed=(/45,8,89,33/)
 call ZLARNV(1,iseed , Ntotal, uglo ) ! random 0-1:  uniform
 uglo(:)=uglo(:)/DZNRM2(Ntotal,uglo,1) ! normalization
 v(1:Nlocal,1)=uglo(startj:startj+Nlocal-1) ! distribute

 comb=.false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Arnoldi direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    H=ZZERO
    conv=.false.
    it=0
    do while ((.not.conv).and.(it<=itmax-1))
       it=it+1
       !mat-vec
         call zhlbcsrmm(UPLO,'N',Nsize,mapA,1,ZONE,matAjb,v(1,it),ZZERO,w,MY_COMM_WORLD,nb_procs)
       !solve
     linloops=itmaxb
     eps=epsb
     uloc(1:Nlocal)=ZZERO !initial guess
 call pzhbicgstab(0,uplo,'N',Nsize,mapB,1,matBjb,w,uloc,res,linloops,eps,comb,MY_COMM_WORLD,nb_procs,infoloc)!!<<
     
       do j=1,it
          !H(j,it)=zdotc(Nlocal,V(1,j),1,uloc,1)
           H(j,it)=dot_product(V(1:Nlocal,j),uloc(1:Nlocal))
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(j,it),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
          call ZAXPY(Nlocal,-H(j,it),V(1,j),1,uloc,1)
       end do
       if (it<itmax) then
          H(it+1,it)=DZNRM2(Nlocal,uloc,1)**2
          call MPI_ALLREDUCE(MPI_IN_PLACE,H(it+1,it),1,MPI_DOUBLE_COMPLEX,MPI_SUM,MY_COMM_WORLD,code)
          H(it+1,it)=sqrt(H(it+1,it))
          call ZCOPY(Nlocal,uloc,1,V(1,it+1),1)
          call ZSCAL(Nlocal,ZONE/H(it+1,it), v(1,it+1), 1)
       end if

         if (((it>=6).and.((it/3)*3==it)).or.(it==itmax)) then ! every 3
!!! solved reduced system to check cv of the extremal eigenvalue (either at 1 or at it) (usually largest amplitude)
        H_temp(1:it,1:it)=H(1:it,1:it)
        call ZHSEQR('E', 'I', it, 1, it, H_temp, itmax, zalpha,  Xred, itmax, WORK,LWORK, INFO_lap)
        !! find lowest/largest eigenvalue
        mine=1
        maxe=1       
        do j=2,it
           if (dble(zalpha(j))> dble(zalpha(maxe))) maxe=j
           if (dble(zalpha(j))< dble(zalpha(mine))) mine=j
        enddo
!!! compute eigenpairs 1 and it and their residuals  R=||AX-X*Lambda||/||X*lambda||
  call ZGEMM('N', 'N', Nlocal, 1, it, ZONE, V(1,1), Nlocal, Xred(1,mine), itmax, ZZERO, Xe(1,1), Nlocal)
        call ZGEMM('N', 'N', Nlocal, 1, it, ZONE, V(1,1), Nlocal, Xred(1,maxe), itmax, ZZERO, Xe(1,2), Nlocal)
        dE(1)=zalpha(mine)
        dE(2)=zalpha(maxe)
        call ZLACPY('F',Nlocal,2,Xe,Nlocal,Ye,Nlocal)
        do j=1,2
           call ZSCAL(Nlocal,dE(j), Ye(1,j), 1)
           rese(j)=DZNRM2(Nlocal,Ye(1,j),1)**2 ! placeholder denominateur
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
 call zhlbcsrmm(UPLO,'N',Nsize,mapA,2,ZONE,matAjb,Xe,ZZERO,aux,MY_COMM_WORLD,nb_procs)

 linloops=itmaxb
     eps=epsb
     Xe(1:Nlocal,1:2)=ZZERO !initial guess
 call pzhbicgstab(0,uplo,'N',Nsize,mapB,2,matBjb,aux,Xe,res,linloops,eps,comb,MY_COMM_WORLD,nb_procs,infoloc)!!<<
      
  do j=1,2
           call ZAXPY(Nlocal,-ZONE,Xe(1,j),1,Ye(1,j),1)
           rese(j)=(DZNRM2(Nlocal,Ye(1,j),1)**2)/rese(j)
        enddo
        
 call MPI_ALLREDUCE(MPI_IN_PLACE,rese,2,MPI_DOUBLE_PRECISION,MPI_SUM,MY_COMM_WORLD,code)
  do j=1,2
     rese(j)=sqrt(rese(j)) ! norm L2
     !   if (j==2)
     !print *,it,j,dE(j),rese(j)
  enddo
  
        if ((rese(1)<epse).or.(rese(2)<epse)) conv=.true.       
     end if
  end do
! itmax=50
!           return
  
itmax=it
    
!! sorting by increasing order (selection)-not optimal (create a function mapping instead)
alpha(1:itmax)=dble(zalpha(1:itmax))
do i=1,itmax
   do j=i+1,itmax
      if (alpha(j)<alpha(i)) then
         temp=alpha(i)
         alpha(i)=alpha(j)
         alpha(j)=temp
         beta(1:itmax)=Xred(1:itmax,i)
         Xred(1:itmax,i)=Xred(1:itmax,j)
         Xred(1:itmax,j)=beta(1:itmax)
      end if
   enddo
   enddo


!!!! output the extremal residual
  res(1)=rese(1)
  res(it)=rese(2)
  
!!!!!!!!!!! Form eigenvectors
  call ZGEMM('N', 'N', Nlocal, itmax, itmax, ZONE, V(1,1), Nlocal, Xred, itmax, ZZERO, X, Nlocal)
 


  deallocate(uglo)
  deallocate(de)
  deallocate(rese)
  deallocate(Xe)
  deallocate(Ye)
  deallocate(aux)
      
  deallocate(uloc)
    deallocate(w)
   deallocate(V) 
   deallocate(H)
    deallocate(H_temp)
   deallocate(work)
   deallocate(Xred)
   deallocate(beta)

end subroutine pzhgarnoldi











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$subroutine dget_distribute_row_from_block_csr(p,nb_procs3,Nsize,mapA,matAjb,matAj)
!!$   implicit none
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$   !      Purpose: This subroutine distribute the matrix by rows- as a condition the matrix must first be distributed using the local block csr format.
!!$  !
!!$  !   Example using 3 MPI:
!!$  
!!$  !    |_mpi1_|              A11 A12 A13
!!$  !  A=|_mpi2_|       <==    A21 A22 A23    
!!$  !    |_mpi3_|              A31 A32 A33
!!$  !
!!$  ! p         (input) INTEGER - rank+1 (row in question)- local
!!$  ! nb_procs3 (input) INTEGER - # of MPI processors
!!$  ! Nsize     (input) INTEGER(:) - size nbprocs
!!$  !            global - store the size of each row blocks
!!$  !  mapA     (input) INTEGER(:,:) - size nbprocs*nbprocs
!!$  !            global - contains the nnz of each submatrix csr blocks
!!$  !        Example:
!!$  !        Each blocks contains a csr matrix, or could be equal to zero;
!!$  !        the number of nnz is provided using the mapping map such that if
!!$  !             |23  12  0 |
!!$  !        map= |12  18  0 |
!!$  !             |0    0  16|
!!$  !
!!$  !        it means that A31=A32=A13=A23=0 and A11 has 23 non-zero elements, etc. 
!!$  !
!!$  ! matAjb (input) Matrix of TYPE(DCSR)- size nbprocs*nbprocs
!!$  !        local- block csr format
!!$  !
!!$  ! matAj  (output) TYPE(DCSR)-
!!$  !        local by row- see explanations above
!!$  !
!!$  !====================================================================
!!$  ! Eric Polizzi 2018-2019
!!$  !====================================================================
!!$  integer :: p,nb_procs3
!!$  integer,dimension(*) :: Nsize  
!!$  integer,dimension(nb_procs3,*) :: mapA
!!$  type dcsr
!!$     integer ::n,m,nnz
!!$     double precision,dimension(:),allocatable :: sa
!!$     integer,dimension(:),allocatable :: isa,jsa     
!!$  end type dcsr
!!$  type(dcsr) :: matAj,temp
!!$  type(dcsr),dimension(nb_procs3,*) :: matAjb
!!$  integer :: i,j,k1,k2,shift
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  
!!$!! create new distribution by columns     
!!$allocate(temp%isa(1:sum(Nsize(p:nb_procs3))+1))
!!$allocate(temp%jsa(1:sum(mapA(p:nb_procs3,p))))
!!$allocate(temp%sa(1:sum(mapA(p:nb_procs3,p))))
!!$
!!$! first column block
!!$do i=1,Nsize(p)
!!$   if (mapA(p,p)>0) then
!!$    k1=matAjb(p,p)%isa(i)
!!$    k2=matAjb(p,p)%isa(i+1)-1
!!$    if (k2>=k1) then
!!$    temp%jsa(k1:k2)=matAjb(p,p)%jsa(k1:k2)
!!$    temp%sa(k1:k2)=matAjb(p,p)%sa(k1:k2)
!!$    endif
!!$    temp%isa(i)=matAjb(p,p)%isa(i)
!!$ else
!!$    temp%isa(i)=1
!!$ end if
!!$ enddo  
!!$! all others column blocks
!!$do j=p+1,nb_procs3
!!$   do i=1,Nsize(j)
!!$      if (mapA(j,p)>0) then
!!$    k1=matAjb(j,p)%isa(i)
!!$    k2=matAjb(j,p)%isa(i+1)-1
!!$    shift=sum(mapA(p:j-1,p))
!!$    if (k2>=k1) then
!!$    temp%jsa(shift+k1:shift+k2)=matAjb(j,p)%jsa(k1:k2)
!!$    temp%sa(shift+k1:shift+k2)=matAjb(j,p)%sa(k1:k2)
!!$    endif
!!$    temp%isa(sum(Nsize(p:j-1))+i)=matAjb(j,p)%isa(i)+shift
!!$ else
!!$    temp%isa(sum(Nsize(p:j-1))+i)=sum(mapA(p:j-1,p))+1
!!$ end if
!!$   end do
!!$end do
!!$ temp%isa(sum(Nsize(p:nb_procs3))+1)=sum(mapA(p:nb_procs3,p))+1
!!$
!!$ 
!!$! transpose the columns into rows
!!$call dcsr_transpose(sum(Nsize(p:nb_procs3)),Nsize(p),temp%sa,temp%isa,temp%jsa,matAj%sa,matAj%isa,matAj%jsa)
!!$! shift the j indexes
!!$if (p>1) then
!!$   matAj%jsa(1:sum(mapA(p:nb_procs3,p)))=matAj%jsa(1:sum(mapA(p:nb_procs3,p)))+sum(Nsize(1:p-1))
!!$end if
!!$
!!$deallocate(temp%isa)
!!$deallocate(temp%jsa)
!!$deallocate(temp%sa)
!!$
!!$end subroutine dget_distribute_row_from_block_csr
!!$



#endif




!!$
!!$
!!$subroutine dget_Chebyshev(opt,Ntotal,M0,sj,ej,x)
!!$implicit none
!!$integer :: opt,Ntotal,M0,sj,ej
!!$double precision, dimension(ej-sj+1,*) :: x
!!$!!!!!
!!$integer ::i,j,Nlocal,k
!!$
!!$Nlocal=ej-sj+1
!!$
!!$if (opt==1) then ! first kind
!!$   x(1:Nlocal,1)=1.0d0
!!$   do i=sj,ej
!!$      x(i-sj+1,2)=-1.0d0+(i-1)*2.0d0/(Ntotal-1)
!!$   end do
!!$   do j=3,M0 ! recurence
!!$      x(1:Nlocal,j)=2.0d0*x(1:Nlocal,2)*x(1:Nlocal,j-1)-x(1:Nlocal,j-2)
!!$   enddo
!!$elseif (opt==2) then ! second kind
!!$   x(1:Nlocal,1)=1.0d0
!!$   do i=sj,ej
!!$      x(i-sj+1,2)=2.0d0*(-1.0d0+(i-1)*2.0d0/(Ntotal-1))
!!$   end do
!!$   do j=3,M0 ! recurence
!!$     ! print *,j,M0
!!$      x(1:Nlocal,j)=x(1:Nlocal,2)*x(1:Nlocal,j-1)-x(1:Nlocal,j-2)
!!$   enddo
!!$end if
!!$   
!!$  
!!$
!!$end subroutine dget_Chebyshev

!!$
!!$subroutine  dlanczos(UPLO,N,sa,isa,jsa,alpha,X,res,itmax)
!!$  !! Purpose: Lanczos algorithm 
!!$
!!$
!!$  !           using modified GS in-place 
!!$  ! N    (input) INTEGER  #rows of the matrix Q
!!$  ! M    (input) INTEGER  #columns of the matrix Q
!!$  ! Q    (input/output) COMPLEX DOUBLE PRECISON (N,M)
!!$  !                     on entry: matrix to orthonormalize
!!$  !                     on exit: orthonormal matrix
!!$  ! R    (output) COMPLEX DOUBLE PRECISION (M,M)- R matrix factors   
!!$  !====================================================================
!!$  ! Eric Polizzi 2019
!!$  !====================================================================  
!!$  implicit none
!!$  character(len=1) :: UPLO
!!$  integer :: N
!!$  double precision,dimension(*) :: sa,alpha,res
!!$  integer,dimension(*) :: isa,jsa
!!$  double precision,dimension(N,*) :: X
!!$  integer :: itmax
!!$!  double precision :: eps
!!$!!!
!!$  integer,dimension(4) :: iseed
!!$   double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
!!$   double precision,dimension(:),allocatable :: r,beta,alpha_tmp,beta_tmp,work
!!$    double precision,dimension(:,:),allocatable :: v,xred
!!$    integer :: i,it
!!$    integer :: info_lap
!!$    double precision :: DNRM2
!!$
!!$   allocate(r(1:N)) 
!!$   allocate(v(n,0:itmax+1)) ! subspace
!!$   allocate(beta(itmax+1))
!!$ !  allocate(alpha_tmp(itmax))
!!$ !  allocate(beta_tmp(itmax-1))
!!$   allocate(xred(itmax,itmax))  !reduced system solution 
!!$   allocate(work(2*itmax-2))
!!$
!!$!!! initialization
!!$   !v0=0
!!$   v(:,0)=0.0d0
!!$   !beta0=0
!!$   beta(1)=0.0d0
!!$   !v1=random
!!$     iseed=(/45,8,89,33/)
!!$     call DLARNV(1,iseed , n, v(1,1) ) ! random 0-1:  uniform
!!$    v(:,1)=v(:,1)/sqrt(dot_product(v(:,1),v(:,1))) ! normalization
!!$!!!!!!!!!!!!!! Lanczos procedure
!!$     it=0
!!$     do while (it<itmax)
!!$        it=it+1 ! # of iterations
!!$        !r=Avj-beta_j*v_{j-1}
!!$        r=-beta(it)*v(:,it-1)
!!$        call wdcsrmm(UPLO,'N',N,N,1,DONE,sa,isa,jsa,v(1,it),DONE,r)
!!$       ! alpha=r^T.v_j
!!$        alpha(it)=dot_product(r,v(:,it))
!!$!!!!!!!!!!!! Solve reduced system (to check convergence)
!!$!        if (it>1) then
!!$!        alpha_tmp(1:it)=alpha(1:it)
!!$!        beta_tmp(1:it-1)=beta(2:it)
!!$!        call  DSTEV('V', it, alpha_tmp, beta_tmp, xred, itmax, WORK, INFO_lap )
!!$!        print *,it,alpha_tmp(1),alpha_tmp(it)
!!$!     end if
!!$!!!!!!!!!!!! prepare re-entry
!!$        ! r=r-alphaj*vj
!!$         r=r-alpha(it)*v(:,it)
!!$         !beta_{j+1}=||r||         
!!$         beta(it+1)=sqrt(dot_product(r,r))
!!$         !v_{j+1}=r/beta_{j+1}
!!$         v(:,it+1)=r/beta(it+1)
!!$     end do
!!$
!!$
!!$!!!!!!!!!!!! Solve reduced system 
!!$        call  DSTEV('V', itmax, alpha, beta(2), xred, itmax, WORK, INFO_lap )
!!$!        print *,'Emin-Emax',alpha(1),alpha(itmax)
!!$     
!!$!!!!!!!!!!! Form eigenvectors
!!$        call DGEMM('N', 'N', N, itmax, itmax, DONE, V(1,1), N, xred, itmax, DZERO, X, N)
!!$
!!$!!!!!!!!!!! Compute residual R=AX-X*Lambda
!!$        call DLACPY('F',N,itmax,X,N,V(1,1),N)
!!$        do i=1,itmax
!!$           call DSCAL(N,alpha(i), V(1,i), 1)
!!$           res(i)=DNRM2(N,V(1,i),1) ! placeholder denominateur
!!$         enddo  
!!$         call wdcsrmm(UPLO,'N',N,N,itmax,DONE,sa,isa,jsa,X,-DONE,V(1,1))
!!$         do i=1,itmax
!!$            res(i)=DNRM2(N,V(1,i),1)/res(i)
!!$!            print *,i,alpha(i),res(i)
!!$         enddo   
!!$
!!$  deallocate(r) 
!!$   deallocate(v) 
!!$   deallocate(beta)
!!$   deallocate(xred)   
!!$   deallocate(work)
!!$   
!!$end subroutine dlanczos
