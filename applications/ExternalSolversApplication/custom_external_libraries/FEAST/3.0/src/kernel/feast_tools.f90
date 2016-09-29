!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Utility Routines for FEAST (Documented)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!! List of routines:
! --------------------
!
!!!!!!! FEASTINIT Routines 
!
! feastinit                      ! default parameters for FEAST input integer array fpm
! feastinit_driver               ! feastinit with additional parameters for predefined drivers (if nedeed by user)


!!!!!!! DEFAULT CONTOUR FEAST Routines
!
! zfeast_contour             ! Return nodes/weights from ellipsoid contour symmetric with real axis- half contour only  
! cfeast_contour             ! Return nodes/weights from ellipsoid contour symmetric with real axis- half contour only  
! zfeast_gcontour            ! Return nodes/weights from ellipsoid contour located in complex plane (general)- return full contour
! cfeast_gcontour            ! Return nodes/weights from ellipsoid contour located in complex plane (general)- return full contour


!!!!!! CUSTOM CONTOUR FEAST Routines (assist the user in generating weight/modes for custom contour in compelx plane)
!
! zfeast_customcontour
! cfeast_customcontour


!!!!!! RATIONAL FEAST Routines
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine feastinit(fpm)
  !  Purpose
  !  =======
  !
  !  Define the default values for the input FEAST parameters.
  !
  !  Arguments
  !  =========
  !
  !  fpm  (output) INTEGER(*): FEAST parameters (size array at least 64)
  !=====================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================
  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  integer,dimension(*) :: fpm
  integer :: i
  do i=1,64
     fpm(i)=0
  enddo

  fpm(1)=0 ! (0,1) comments on/off
  fpm(2)=8 ! # of contour integration nodes (half-contour) - Hermitian
  fpm(3)=12! tol double precision
  fpm(4)=20! maxloop
  fpm(5)=0 ! (0,1) Initial Subspace yes/no 
  fpm(6)=1 ! Convergence criteria on trace (0) / eigenvectors relative residual (1)
  fpm(7)=5 ! tol single precision

  fpm(8)=16 ! # of contour points (full-contour) - Non-Hermitian    

  !---------------------------------------
#ifdef MPI
  fpm(9)=MPI_COMM_WORLD ! default value
#endif
  !----------------------------------------

  fpm(10)=0 ! used for the FEAST Predefined driver interfaces (0: default, 1: store all the factorizations)
  ! Remark: work with mpi as well - store all factors associated to a given mpi process  

  fpm(12)=0 ! (0,1) customize eigenvalue solver for rci interface (No,Yes) --undocumented--

  fpm(13)=1 ! (0,1) spurious detection (No,Yes) --undocumented--

  fpm(14)=0 ! (0,1,2) 
  ! 1- return only subspace Q size M0 after 1 contour          
  ! 2- return stochastic estimates of the eigenvalue count

  fpm(15)=1 ! bi-orthogonalization flag for non-symmetric (No,Yes)  

  fpm(16)=0 ! Integration type for symmetric (0: Gauss, 1: Trapezoidal, 2: Zolotarev)
  fpm(17)=1 ! Integration type for non-symmetric (0: Gauss, 1: Trapezoidal)

  fpm(18)=100  ! ellipsoid contour - for Symmetric - fpm(18)/100 is ratio a/b (b is [Emin-Emax]) - Rq: circle is fpm(18)=100
  fpm(19)=0 ! Rotation angle in degree for ellispoid contour - for Non-symmetric [-180,180]



  fpm(64)=0 ! Additional feast parameters for driver interfaces (i.e size fpm>64) (0,1)  --undocumented--

end subroutine feastinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine feastinit_driver(fpm,N)
  !  Purpose
  !  =======
  !
  !  Define the default values for the input FEAST parameters from 1-64
  !  and 65-N are optional user defined inputs that can be used for FEAST predefined interfaces,
  !  here fpm(65:N) is initialized at -9999.
  !
  !  Arguments
  !  =========
  !
  !  fpm  (output) INTEGER(*): FEAST parameters (size array at least 64)
  !  N    (input)  INTEGER: size of array fpm (>=64)
  !=====================================================================
  ! Eric Polizzi 2009-2013
  ! ====================================================================
  implicit none
  !-------------------------------------
#ifdef MPI
  include 'mpif.h'
#endif
  !-------------------------------------
  integer,dimension(*) :: fpm
  integer :: N
  integer :: i

  call feastinit(fpm)

  if (N>64) then
     fpm(64)=1
     do i=65,N
        fpm(i)=-9999 ! default values for the drivers inputs
     enddo
  end if
end subroutine feastinit_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



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
  include "f90_noruntime_interface.fi"
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
        call dset_feast_gauss_legendre(fpm2,e,xe,we) !! Gauss-points 
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+r*ONEC*wdcos(theta)+r*(DZERO,DONE)*(fpm18*1d-2)*wdsin(theta)
        jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*(fpm18*1d-2)*wdcos(theta))
        Wne(e)= (DONE/4.0d0)*we*jac
     elseif(fpm16==2)then
        call dset_feast_zolotarev(fpm2,e,zxe,zwe)
        Zne(e)=zxe*r+Emid*ONEC
        Wne(e)=zwe*r
     elseif(fpm16==1)then
        theta=pi-(pi/fpm2)/2.0d0-(pi/fpm2)*(e-1)
        Zne(e)=Emid*ONEC+r*ONEC*wdcos(theta)+r*(DZERO,DONE)*(fpm18*1d-2)*wdsin(theta)
        jac=(r*(DZERO,DONE)*wdsin(theta)+ONEC*r*(fpm18*1d-2)*wdcos(theta))
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
  include "f90_noruntime_interface.fi"
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
        call sset_feast_gauss_legendre(fpm2,e,xe,we) !! Gauss-points 
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+r*ONEC*wscos(theta)+r*(SZERO,SONE)*(fpm18*1e-2)*wssin(theta)
        jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*(fpm18*1e-2)*wscos(theta))
        Wne(e)= (SONE/4.0e0)*we*jac
     elseif(fpm16==2)then
        call sset_feast_zolotarev(fpm2,e,zxe,zwe)
        Zne(e)=zxe*r+Emid*ONEC
        Wne(e)=zwe*r
     elseif(fpm16==1)then
        theta=pi-(pi/fpm2)/2.0e0-(pi/fpm2)*(e-1)
        Zne(e)=Emid*ONEC+r*ONEC*wscos(theta)+r*(SZERO,SONE)*(fpm18*1e-2)*wssin(theta)
        jac=(r*(SZERO,SONE)*wssin(theta)+ONEC*r*(fpm18*1e-2)*wscos(theta))
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
  include "f90_noruntime_interface.fi"
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
  nr=r*(ONEC*wdcos(theta)+(DZERO,DONE)*wdsin(theta))


  if(fpm17==0)then


     do e=1,fpm8/2
        call dset_feast_gauss_legendre(fpm8/2,e,xe,we) !! Gauss-points - upper half
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+nr*ONEC*wdcos(theta)+nr*(DZERO,DONE)*(fpm18*1d-2)*wdsin(theta)
        jac=(nr*(DZERO,DONE)*wdsin(theta)+ONEC*nr*(fpm18*1d-2)*wdcos(theta))
        Wne(e) = (ONEC/(4.0d0))*we*jac
     enddo

     do e=fpm8/2+1,fpm8
        call dset_feast_gauss_legendre(fpm8-fpm8/2,e-fpm8/2,xe,we) !! Gauss-points - lower half
        theta=-ba*xe-ab
        Zne(e)=Emid*ONEC+nr*ONEC*wdcos(theta)+nr*(DZERO,DONE)*(fpm18*1d-2)*wdsin(theta)
        jac=nr*(DZERO,DONE)*wdsin(theta)+ONEC*nr*(fpm18*1d-2)*wdcos(theta)
        Wne(e) = (ONEC/(4.0d0))*we*jac
     enddo



  elseif(fpm17==1)then

     do e=1,fpm8
        theta=pi-(2*pi/fpm8)/2.0d0-(2*pi/fpm8)*(e-1)
        Zne(e)=Emid*ONEC+nr*ONEC*wdcos(theta)+nr*(DZERO,DONE)*(fpm18*1d-2)*wdsin(theta)
        jac=(nr*(DZERO,DONE)*wdsin(theta)+ONEC*nr*(fpm18*1d-2)*wdcos(theta))
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
  include "f90_noruntime_interface.fi"
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
  nr=r*(ONEC*wscos(theta)+(SZERO,SONE)*wssin(theta))



  if(fpm17==0)then

     do e=1,fpm8/2
        call sset_feast_gauss_legendre(fpm8/2,e,xe,we) !! Gauss-points - upper half
        theta=ba*xe+ab
        Zne(e)=Emid*ONEC+nr*ONEC*wscos(theta)+nr*(SZERO,SONE)*(fpm18*1e-2)*wssin(theta)
        jac=(nr*(SZERO,SONE)*wssin(theta)+ONEC*nr*(fpm18*1e-2)*wscos(theta))
        Wne(e) = (ONEC/(4.0e0))*we*jac
     enddo

     do e=fpm8/2+1,fpm8
        call sset_feast_gauss_legendre(fpm8-fpm8/2,e-fpm8/2,xe,we)!! Gauss-points - lower half
        theta=-ba*xe-ab
        Zne(e)=Emid*ONEC+nr*ONEC*wscos(theta)+nr*(SZERO,SONE)*(fpm18*1e-2)*wssin(theta)
        jac=(nr*(SZERO,SONE)*wssin(theta)+ONEC*nr*(fpm18*1e-2)*wscos(theta))
        Wne(e) = (ONEC/(4.0e0))*we*jac
     enddo


  elseif(fpm17==1)then

     do e=1,fpm8
        theta=pi-(2*pi/fpm8)/2.0e0-(2*pi/fpm8)*(e-1)
        Zne(e)=Emid*ONEC+nr*ONEC*wscos(theta)+nr*(SZERO,SONE)*(fpm18*1e-2)*wssin(theta)
        jac=(nr*(SZERO,SONE)*wssin(theta)+ONEC*nr*(fpm18*1e-2)*wscos(theta))
        Wne(e)=(ONEC/fpm8)*jac
     enddo

  endif

end subroutine cfeast_gcontour



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
  include "f90_noruntime_interface.fi"
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
        angle = wdatan2(aimag(b-a),dble(b-a))
        Emid = (a+b)/2.0d0
        r = abs(a-b)/2.0d0
        r2= Tedge(i)*r/100.0d0
        do j=1,Nedge(i)+1
           theta = pi-1.0d0*(j-1)*pi/(Nedge(i))
           tZne(cnt+j) = Emid +  (1.0d0,0.0d0)*(r*wdcos(theta)*wdcos(angle)-r2*wdsin(theta)*wdsin(angle)) + (0.0d0,1.0d0)*(r*wdcos(theta)*wdsin(angle)+r2*wdsin(theta)*wdcos(angle))
           tWne(cnt+j) = -r*wdsin(theta)*wdcos(angle) - r2*wdcos(theta)*wdsin(angle) + (0.0d0,1.0d0)*(-r*wdsin(theta)*wdsin(angle) + r2*wdcos(theta)*wdcos(angle))
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
  include "f90_noruntime_interface.fi"
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
        angle = wsatan2(aimag(b-a),real(b-a))
        Emid = (a+b)/2.0e0
        r = abs(a-b)/2.0e0
        r2= Tedge(i)*r/100.0e0
        do j=1,Nedge(i)+1
           theta = pi-1.0*(j-1)*pi/(Nedge(i))
           tZne(cnt+j) = Emid +  (1.0e0,0.0e0)*(r*wscos(theta)*wscos(angle)-r2*wssin(theta)*wssin(angle)) + (0.0e0,1.0e0)*(r*wscos(theta)*wssin(angle)+r2*wssin(theta)*wscos(angle))
           tWne(cnt+j) = -r*wssin(theta)*wscos(angle) - r2*wscos(theta)*wssin(angle) + (0.0e0,1.0e0)*(-r*wssin(theta)*wssin(angle) + r2*wscos(theta)*wscos(angle))
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
     call dset_feast_zolotarev(fpm2,0,zxe,zwe)
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
     call sset_feast_zolotarev(fpm2,0,zxe,zwe)
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







