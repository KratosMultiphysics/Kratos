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
!!!!!! Auxiliary Routines for FEAST (undocumented)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!! List of routines:
! --------------------
!!!!!!!! Find Name of routines
! feast_name              ! from fpm(30) code integer, find name of routines


!!!!!!!! Check Inputs Routines
! dcheck_feast_srci_input             ! Check validity of some RCI FEAST inputs - Symmetric/Hermitian drivers- double precision 
! scheck_feast_srci_input             ! Check validity of some RCI FEAST inputs - Symmetric/Hermitian drivers- single precision 
! dcheck_feast_grci_input             ! Check validity of some RCI FEAST inputs - Non-Hermitian drivers- double precision 
! scheck_feast_grci_input             ! Check validity of some RCI FEAST inputs - Non-Hermitian drivers- single precision 

!!!!!!!  Check Location eigenvalue routines
! zfeast_inside_contourx_old   ! unused-
! zfeast_bary_coef             ! unused-
! zfeast_inside_contourx       ! tentative for alternative
! zfeast_inside_contour

! cfeast_inside_contourx_old ! unused-
! cfeast_bary_coef           ! unused- 
! cfeast_inside_contourx     ! tentative for alternative
! cfeast_inside_contour


!dfeast_info                ! print feast data for symmetric or Hermitian
!zfeast_info                ! print feast data for generalized

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine feast_name(code,name)
  !  Purpose 
  !  =======
  !  Find name of FEAST routine from the integer code
  !
  !  Arguments
  !  =========
  !
  !  code   (input) INTEGER : FEAST fpm(30)
  !                           up to 6 digits e.g. 123444 (from digit 1 to digit 6)
  !                         digit 1: 1- regular feast; 2- pfeast
  !                         digit 2: 1-s; 2-d; 3-c; 4-z
  !                         digit 3: 1- direct feast; 2- ifeast
  !                         digit 4: 1-S; 2-H; 3-G
  !                         digit 5: 1- RCI; 2- Y; 3- B; 4- CSR; 5- E
  !                         digit 6: 0-"empty"; 1-X; 2-EV; 3-EVX; 4-GV; 5-GVX; 6- PEP; 7- PEPX  
  !  name   (output) Character(len=*) : FEAST name
  !
  !=====================================================================
  ! Eric Polizzi 2018-2019
  ! ====================================================================
  implicit none
  integer :: code 
  character(len=25) :: name

  integer,dimension(6) :: dig
  integer :: i,rem

  rem = code
  DO i = 1, 6
     dig(6-i+1) = rem - (rem/10)*10  ! Take advantage of integer division
     rem = rem/10
  END DO

  name=""
!!!
  select case(dig(1))
  case(2)
     name=trim(name)//"P"
  end select
!!!
  select case(dig(2))
  case(1)
     name=trim(name)//"S"
  case(2)
     name=trim(name)//"D"
  case(3)
     name=trim(name)//"C"
  case(4)
     name=trim(name)//"Z"      
  end select
!!!

  select case(dig(3))
  case(1)
     name=trim(name)//"FEAST_"
  case(2)
     name=trim(name)//"IFEAST_"   
  end select
!!!
  select case(dig(4))
  case(1)
     name=trim(name)//"S"
  case(2)
     name=trim(name)//"H"
  case(3)
     name=trim(name)//"G"    
  end select
!!!
  select case(dig(5))
  case(1)
     name=trim(name)//"RCI"
  case(2)
     name=trim(name)//"Y"
  case(3)
     name=trim(name)//"B"
  case(4)
     name=trim(name)//"CSR"
  case(5)
     name=trim(name)//"E"           
  end select
!!!
  select case(dig(6))
  case(1)
     name=trim(name)//"X"
  case(2)
     name=trim(name)//"EV"
  case(3)
     name=trim(name)//"EVX"
  case(4)
     name=trim(name)//"GV"
  case(5)
     name=trim(name)//"GVX"
  case(6)
     name=trim(name)//"PEV"
  case(7)
     name=trim(name)//"PEVX"   
  end select

end subroutine feast_name




subroutine dcheck_feast_srci_input(Emin,Emax,M0,N,info)
  !  Purpose 
  !  =======
  !  Error handling for the FEAST RCI double precision interfaces input parameters.
  !  Check the values of Emin, Emax, M0, N, and return 
  !  info code error /=0 if incorrect values are found
  !
  !  Arguments
  !  =========
  !
  !  Emin,Emax   (input) REAL DOUBLE PRECISION: search interval
  !  M0          (input) INTEGER: Size subspace
  !  N           (input) INTEGER: Size system
  !  info (input/output) INTEGER
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================
  implicit none
  double precision :: Emin,Emax
  integer :: N,M0,info

  if (Emin>=Emax) info=200 ! problem with Emin, Emax
  if ((M0<=1).or.(M0>N)) info=201 ! problem with M0 
  if (N<=0) info=202 ! problem with N

end subroutine dcheck_feast_srci_input





subroutine scheck_feast_srci_input(Emin,Emax,M0,N,info)
  !  Purpose 
  !  =======
  !  Error handling for the FEAST RCI single precision interfaces input parameters.
  !  Check the values of Emin, Emax, M0, N, and return 
  !  info code error /=0 if incorrect values are found
  !
  !  Arguments
  !  =========
  !
  !  Emin,Emax   (input) REAL SINGLE PRECISION: search interval
  !  M0          (input) INTEGER: Size subspace
  !  N           (input) INTEGER: Size system
  !  info (input/output) INTEGER
  !=====================================================================
  ! Eric Polizzi 2009
  !=====================================================================
  implicit none
  real :: Emin,Emax
  integer :: N,M0,info

  if (Emin>=Emax) info=200 ! problem with Emin, Emax
  if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
  if (N<=0) info=202 ! problem with N

end subroutine scheck_feast_srci_input




subroutine dcheck_feast_grci_input(r,M0,N,info)
  !  Purpose 
  !  =======
  !  Error handling for the FEAST RCI double precision interfaces input parameters.
  !  Check the values of Emid, r, M0, N, and return 
  !  info code error /=0 if incorrect values are found
  !
  !  Arguments
  !  =========
  !
  !  r           (input) REAL DOUBLE PRECISION: search interval (radius)
  !  M0          (input) INTEGER: Size subspace
  !  N           (input) INTEGER: Size system
  !  info (input/output) INTEGER
  !=====================================================================
  ! Eric Polizzi 2015
  ! ====================================================================
  implicit none
  double precision :: r
  integer :: N,M0,info

  if (r<=0.0d0) info=200 ! problem with r
  if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
  if (N<=0) info=202 ! problem with N

end subroutine dcheck_feast_grci_input



subroutine scheck_feast_grci_input(r,M0,N,info)
  !  Purpose 
  !  =======
  !  Error handling for the FEAST RCI real interfaces input parameters.
  !  Check the values of Emid, r, M0,N, and return 
  !  info code error /=0 if incorrect values are found
  !
  !  Arguments
  !  =========
  !
  !  r      (input) REAL DOUBLE PRECISION: search interval (radius)
  !  M0          (input) INTEGER: Size subspace
  !  N           (input) INTEGER: Size system
  !  info (input/output) INTEGER
  !=====================================================================
  ! Eric Polizzi 2015
  ! ====================================================================
  implicit none
  real :: r
  integer :: N,M0,info

  if (r<=0.0d0) info=200 ! problem with r
  if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
  if (N<=0) info=202 ! problem with N

end subroutine scheck_feast_grci_input




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_inside_contourx_old(Zne,fpm2,eig,M0,col) 
  !  Purpose 
  !  =======
  !  Check if a set of eigenvalues are located inside a polytope custom contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  DOUBLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  Zne    (input)  COMPLEX DOUBLE PRECISION (fpm2)   : list of coordinate nodes for the polytope custom contour 
  !  fpm2   (input)  INTEGER                           : # of coordinates nodes
  !  eig    (input)  COMPLEX DOUBLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn 2015
  ! ====================================================================
  implicit none
  integer :: fpm2,M0
  integer,dimension(M0) :: col
  complex(kind=(kind(1.0d0))), dimension(M0) :: eig
  complex(kind=(kind(1.0d0))), dimension(fpm2) :: Zne
!!!!!!!!!!!!!!

  double precision :: a0,a1,a2,alphaj,summ
  double precision,dimension(1:fpm2) :: wj
  integer :: i,j


  col(1:M0)=1

  Do j=1,M0

     summ=0.0d0

     if(eig(j) /= eig(j)) then !NaN
        col(j) = 0
     endif

     if (col(j)==1) then
        do i=1,fpm2
           call zfeast_bary_coef(Zne(mod(i+fpm2-2,fpm2)+1),Zne(i),Zne(mod(i,fpm2)+1),a0)
           call zfeast_bary_coef(Zne(mod(i+fpm2-2,fpm2)+1),Zne(i),eig(j),a1)
           call zfeast_bary_coef(Zne(i),Zne(mod(i,fpm2)+1),eig(j),a2)
           if ((a1*a2)/=0.0d0) then
              wj(i) = a0/(a1*a2)
              summ=summ+wj(i)
           else
              col(j)=0
           end if
        enddo
     end if

     if (col(j)==1) then
        do i=1,fpm2
           alphaj = wj(i)/summ
           if(alphaj < -1.0d-12) then
              col(j) = 0
           endif
        enddo
     end if

  end Do

end subroutine zfeast_inside_contourx_old



subroutine zfeast_bary_coef(v1,v2,v3,coef)
  !  Purpose 
  !  =======
  !  Auxiliary routine for zfeasst_inside_contourx
  !  Compute barycentric coefficient
  !
  !  DOUBLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  v1,v2,v3  (input)  COMPLEX DOUBLE PRECISION 
  !  coef      (output) REAL DOUBLE PRECISION
  !=====================================================================
  ! James Kestyn 2015
  ! ====================================================================
  implicit none
  double precision :: coef
  complex(kind=(kind(1.0d0))) :: v1,v2,v3
  double precision :: adbc1,adbc2,adbc3
  adbc1 = 1.0d0 * ( dble(v2)*aimag(v3) - dble(v3)*aimag(v2) )
  adbc2 = 1.0d0 * ( dble(v1)*aimag(v3) - dble(v3)*aimag(v1) )
  adbc3 = 1.0d0 * ( dble(v1)*aimag(v2) - dble(v2)*aimag(v1) )
  coef = adbc1 - adbc2 + adbc3

end subroutine zfeast_bary_coef





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zfeast_inside_contourx(Zne,fpm2,eig,M0,col) 
  !  Purpose 
  !  =======
  !  Check if a set of eigenvalues are located inside a polytope custom contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  DOUBLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  Zne    (input)  COMPLEX DOUBLE PRECISION (fpm2)   : list of coordinate nodes for the polytope custom contour 
  !  fpm2   (input)  INTEGER                           : # of coordinates nodes
  !  eig    (input)  COMPLEX DOUBLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn 2015
  ! ====================================================================
  implicit none
  integer :: fpm2,M0
  integer,dimension(M0) :: col
  complex(kind=(kind(1.0d0))), dimension(M0) :: eig
  complex(kind=(kind(1.0d0))), dimension(fpm2) :: Zne
!!!!!!!!!!!!!!
  complex(kind=(kind(1.0d0))) :: z1i,z1j
  double precision :: xp,yp,x1,y1,x2,y2,x3,y3
  double precision :: lambda1,lambda2,lambda3
  integer :: i,j,k


  col(1:M0)=0
  x1 = dble(Zne(1))
  y1 = aimag(Zne(1))

  do i=2,fpm2
     x2 = dble(Zne(i))
     y2 = aimag(Zne(i))
     z1i = (Zne(i) - Zne(1))/abs(Zne(i) - Zne(1))
     do j=i+1,fpm2                  
        x3 = dble(Zne(j))
        y3 = aimag(Zne(j))
        z1j = (Zne(j) - Zne(1))/abs(Zne(j) - Zne(1))

        if( abs(1.0d0 - abs(dble(z1i)*dble(z1j) + aimag(z1i)*aimag(z1j))) > 1.0d-8) then
           Do k=1,M0
              xp = dble(eig(k))
              yp = aimag(eig(k))
              lambda1 = ((y2-y3)*(xp-x3)+(x3-x2)*(yp-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
              lambda2 = ((y3-y1)*(xp-x3)+(x1-x3)*(yp-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
              lambda3 = 1.0d0-lambda1-lambda2
              if(lambda1>-1.0d-14 .and. lambda2>-1.0d-14 .and. lambda3>-1.0d-14)then
                 col(k)=1
              endif
           enddo
        end if
     enddo
  enddo

end subroutine zfeast_inside_contourx




subroutine zfeast_inside_contour(Emid,r,fpm18,fpm19,eig,M0,col)
  !  Purpose
  !  =======
  !  Check if a set of eigenvalues are located inside an ellipsoid contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  DOUBLE PRECISION version
  !
  !  Arguments
  !  =========
  !
  !  Emid   (input)  COMPLEX DOUBLE PRECISION          : Center of the ellipse
  !  r      (input)  REAL DOUBLE PRECISION             : Radius of horizontal axis b
  !  fpm18  (input)  INTEGER                           : fpm(18)*0.01 is ratio a/b (a is vertical axis) -- Rq: circle is 100
  !  fpm19 (input)  INTEGER     Rotation angle for the ellipse in degree [-180:180]
  !  eig    (input)  COMPLEX DOUBLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ====================================================================
  implicit none
  !include "f90_noruntime_interface.fi"
  integer :: M0,fpm18,fpm19
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer,dimension(M0) :: col
  complex(kind=(kind(1.0d0))), dimension(M0) :: eig
!!!!!!!!!!!!!!
  double precision,parameter :: pi=3.1415926535897932d0
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  complex(kind=(kind(1.0d0))) :: nE
  double precision :: l, theta
  integer :: j

  !! ellipse axis rotation (rotate E)
  theta=-(fpm19/180d0)*pi

  col(1:M0)=1

  do j=1,M0
     if(eig(j) /= eig(j)) then !NaN
        col(j) = 0
     endif

     if (col(j)==1) then
        nE = (eig(j)-Emid)*(ONEC*cos(theta)+(DZERO,DONE)*sin(theta))+Emid
        l=abs(nE-Emid)
        if ((r*r*(fpm18*1d-2))<l*(sqrt((r*(fpm18*1d-2)*dble(nE-Emid)/l)**2+(r*aimag(nE-Emid)/l)**2))) col(j)=0 ! polar coordinates
     endif
  enddo

end subroutine zfeast_inside_contour



subroutine zfeast_inside_contour_old(Emid,r,fpm19,eig,M0,col) 
  !  Purpose 
  !  =======
  !  Check if a set of eigenvalues are located inside an ellipsoid contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  DOUBLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  Emid   (input)  COMPLEX DOUBLE PRECISION          : Center of the ellipse
  !  r      (input)  REAL DOUBLE PRECISION             : Radius of horizontal axis b 
  !  fpm19  (input)  INTEGER                           : fpm(19)*0.01 is ratio a/b (a is vertical axis) -- Rq: circle is 100
  !  eig    (input)  COMPLEX DOUBLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ====================================================================
  implicit none
  integer :: M0,fpm19
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r
  integer,dimension(M0) :: col
  complex(kind=(kind(1.0d0))), dimension(M0) :: eig
!!!!!!!!!!!!!!
  !double precision :: dZr ! distance from the center to the ellipse surface at a given angle (defined by eig)
  integer :: j
  double precision :: l

  col(1:M0)=1

  do j=1,M0
     if(eig(j) /= eig(j)) then !NaN
        col(j) = 0
     endif

     if (col(j)==1) then
        l=abs(eig(j)-Emid)
        if ((r*r*(fpm19*1d-2))<l*(sqrt((r*(fpm19*1d-2)*dble(eig(j)-Emid)/l)**2+(r*aimag(eig(j)-Emid)/l)**2))) col(j)=0 ! polar coordinates
        !dZr=((r*r*(fpm19*1d-2))/(sqrt((r*(fpm19*1d-2)*dble(eig(j)-Emid)/l)**2+(r*aimag(eig(j)-Emid)/l)**2))) ! polar coordinate
        !if (l>abs(dzr)) col(j)=0 
     end if
  enddo

end subroutine zfeast_inside_contour_old




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cfeast_inside_contourx_old(Zne,fpm2,eig,M0,col) 
  !  Purpose 
  !  =======
  !  Check if a set of eigenvalues are located inside a polytope custom contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  SINGLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  Zne    (input)  COMPLEX SINGLE PRECISION (fpm2)   : list of coordinate nodes for the polytope custom contour 
  !  fpm2   (input)  INTEGER                           : # of coordinates nodes
  !  eig    (input)  COMPLEX SINGLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn 2015
  ! ====================================================================
  implicit none
  integer :: fpm2,M0
  integer,dimension(M0) :: col
  complex, dimension(M0) :: eig
  complex, dimension(fpm2) :: Zne
!!!!!!!!!!!!!!

  real:: a0,a1,a2,alphaj,summ
  real,dimension(1:fpm2) :: wj
  integer :: i,j


  col(1:M0)=1

  Do j=1,M0

     summ=0.0e0

     if(eig(j) /= eig(j)) then !NaN
        col(j) = 0
     endif

     if (col(j)==1) then
        do i=1,fpm2
           call cfeast_bary_coef(Zne(mod(i+fpm2-2,fpm2)+1),Zne(i),Zne(mod(i,fpm2)+1),a0)
           call cfeast_bary_coef(Zne(mod(i+fpm2-2,fpm2)+1),Zne(i),eig(j),a1)
           call cfeast_bary_coef(Zne(i),Zne(mod(i,fpm2)+1),eig(j),a2)
           if ((a1*a2)/=0.0e0) then
              wj(i) = a0/(a1*a2)
              summ=summ+wj(i)
           else
              col(j)=0
           end if
        enddo
     end if

     if (col(j)==1) then
        do i=1,fpm2
           alphaj = wj(i)/summ
           if(alphaj < -1.0e-5) then
              col(j) = 0
           endif
        enddo
     end if

  end Do

end subroutine cfeast_inside_contourx_old




subroutine cfeast_bary_coef(v1,v2,v3,coef)
  !  Purpose 
  !  =======
  !  Auxiliary routine for zfeasst_inside_contourx
  !  Compute barycentric coefficient
  !
  !  SINGLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  v1,v2,v3  (input)  COMPLEX SINGLE PRECISION 
  !  coef      (output) REAL SINGLE PRECISION
  !=====================================================================
  ! James Kestyn 2015
  ! ====================================================================
  implicit none
  real :: coef
  complex :: v1,v2,v3
  real :: adbc1,adbc2,adbc3
  adbc1 = 1.0e0 * ( real(v2)*aimag(v3) - real(v3)*aimag(v2) )
  adbc2 = 1.0e0 * ( real(v1)*aimag(v3) - real(v3)*aimag(v1) )
  adbc3 = 1.0e0 * ( real(v1)*aimag(v2) - real(v2)*aimag(v1) )
  coef = adbc1 - adbc2 + adbc3

end subroutine cfeast_bary_coef





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cfeast_inside_contourx(Zne,fpm2,eig,M0,col) 
  !  Purpose 
  !  =======
  !  Check if a set of eigenvalues are located inside a polytope custom contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  SINGLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  Zne    (input)  COMPLEX SINGLE PRECISION (fpm2)   : list of coordinate nodes for the polytope custom contour 
  !  fpm2   (input)  INTEGER                           : # of coordinates nodes
  !  eig    (input)  COMPLEX SINGLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn 2015
  ! ====================================================================
  implicit none
  integer :: fpm2,M0
  integer,dimension(M0) :: col
  complex, dimension(M0) :: eig
  complex, dimension(fpm2) :: Zne
!!!!!!!!!!!!!!
  complex :: z1i,z1j
  real :: xp,yp,x1,y1,x2,y2,x3,y3
  real :: lambda1,lambda2,lambda3
  integer :: i,j,k


  col(1:M0)=0
  x1 = real(Zne(1))
  y1 = aimag(Zne(1))

  do i=2,fpm2
     x2 = real(Zne(i))
     y2 = aimag(Zne(i))
     z1i = (Zne(i) - Zne(1))/abs(Zne(i) - Zne(1))
     do j=i+1,fpm2
        x3 = real(Zne(j))
        y3 = aimag(Zne(j))
        z1j = (Zne(j) - Zne(1))/abs(Zne(j) - Zne(1))

        if( abs(1.0e0 - abs(real(z1i)*real(z1j) + aimag(z1i)*aimag(z1j))) > 1.0e-4) then
           Do k=1,M0
              xp = real(eig(k))
              yp = aimag(eig(k))

              lambda1 = ((y2-y3)*(xp-x3)+(x3-x2)*(yp-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
              lambda2 = ((y3-y1)*(xp-x3)+(x1-x3)*(yp-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
              lambda3 = 1.0e0-lambda1-lambda2
              if(lambda1>-1.0e-5 .and. lambda2>-1.0e-5 .and. lambda3>-1.0e-5)then
                 col(k)=1
              endif
           enddo
        end if
     enddo
  enddo


end subroutine cfeast_inside_contourx



subroutine cfeast_inside_contour(Emid,r,fpm18,fpm19,eig,M0,col) 
  !  Purpose 
  !  =======
  !  Check if a set of eigenvalues are located inside an ellipsoid contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  SINGLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  Emid   (input)  COMPLEX SINGLE PRECISION          : Center of the ellipse
  !  r      (input)  REAL SINGLE PRECISION              : Radius of horizontal axis b 
  !  fpm18  (input)  INTEGER                           : fpm(19)*0.01 is ratio a/b (a is vertical axis) -- Rq: circle is 100
  !  fpm19 (input)  INTEGER     Rotation angle for the ellipse in degree [-180:180]
  !  eig    (input)  COMPLEX SINGLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ====================================================================
  implicit none
  !include "f90_noruntime_interface.fi"
  integer :: M0,fpm18,fpm19
  complex :: Emid
  real :: r
  integer,dimension(M0) :: col
  complex, dimension(M0) :: eig
!!!!!!!!!!!!!!
  real, parameter :: pi=3.1415926535897932e0
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
  complex :: nE
  real :: l, theta
  integer :: j

  !! ellipse axis rotation (rotate E)
  theta=-(fpm19/180e0)*pi


  col(1:M0)=1
  do j=1,M0
     if(eig(j) /= eig(j)) then !NaN
        col(j) = 0
     endif

     if (col(j)==1) then
        nE = (eig(j)-Emid)*(ONEC*cos(theta)+(SZERO,SONE)*sin(theta))+Emid
        l=abs(nE-Emid)
        if ((r*r*(fpm18*1e-2))<l*(sqrt((r*(fpm18*1e-2)*real(nE-Emid)/l)**2+(r*aimag(nE-Emid)/l)**2))) col(j)=0 ! polar coordinates
     end if
  enddo

end subroutine cfeast_inside_contour
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine cfeast_inside_contour_old(Emid,r,fpm19,eig,M0,col) 
  !  Purpose 
  !  =======
  !  Check if a set of eigenvalues are located inside an ellipsoid contour
  !  Color of the eigenvalue is set to 1 if it is inside, 0 if not
  !
  !  SINGLE PRECISION version
  !
  !  Arguments
  !  =========
  !  
  !  Emid   (input)  COMPLEX SINGLE PRECISION          : Center of the ellipse
  !  r      (input)  REAL SINGLE PRECISION              : Radius of horizontal axis b 
  !  fpm19  (input)  INTEGER                           : fpm(19)*0.01 is ratio a/b (a is vertical axis) -- Rq: circle is 100
  !  eig    (input)  COMPLEX SINGLE PRECISION (M0)     : list of eigenvalues to consider
  !  M0     (input)  INTEGER                           : # of eigenvalues
  !  color  (output) INTEGER (M0)                      : color array - return 1 if eigenvalue inside, 0 otherwise
  !=====================================================================
  ! James Kestyn, Eric Polizzi 2015
  ! ====================================================================
  implicit none
  integer :: M0,fpm19
  complex :: Emid
  real :: r
  integer,dimension(M0) :: col
  complex, dimension(M0) :: eig
!!!!!!!!!!!!!!
  !  double precision :: dZr ! distance from the center to the ellipse surface at a given angle (defined by eig)
  integer :: j
  double precision :: l

  col(1:M0)=1
  do j=1,M0
     if(eig(j) /= eig(j)) then !NaN
        col(j) = 0
     endif

     if (col(j)==1) then
        l=abs(eig(j)-Emid)
        if ((r*r*(fpm19*1e-2))<l*(sqrt((r*(fpm19*1e-2)*real(eig(j)-Emid)/l)**2+(r*aimag(eig(j)-Emid)/l)**2))) col(j)=0 ! polar coordinates
        !dZr=r*r*(fpm19*1e-2)/(sqrt((r*(fpm19*1e-2)*real(eig(j)-Emid)/l)**2+(r*aimag(eig(j)-Emid)/l)**2)) ! polar coordinates
        !if (l>abs(dzr)) col(j)=0
     end if
  enddo

end subroutine cfeast_inside_contour_old
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine dfeast_info(fout,fpm,Emin,Emax,N,M0)
  !!! Eric Polizzi 2019
  implicit none
  integer(8) :: fout
  integer,dimension(*) :: fpm
  double precision :: Emin,Emax
  integer :: N,M0
!!!!!!!!!!!!!

  write(fout,*)
  write(fout,'(A)') '.-------------------.'
  write(fout,'(A)') '| FEAST data        |'
  write(fout,'(A)') '--------------------.-------------------------.'
  write(fout,'(A)',advance='no') '| Emin              |'
  write(fout,'(ES24.16)',advance='no') Emin
  write(fout,'(A)') ' |'
  write(fout,'(A)',advance='no') '| Emax              |'
  write(fout,'(ES24.16)',advance='no') Emax
  write(fout,'(A)') ' |'
  write(fout,'(A)',advance='no') '| #Contour nodes    |'
  write(fout,'(I3)',advance='no') fpm(2)
  write(fout,'(A)') '   (half-contour)     |'
  write(fout,'(A)',advance='no') '| Quadrature rule   |'
  if (fpm(29)==1) then
     if (fpm(16)==0) then
        write(fout,'(A)') '  Gauss                  |'
     elseif(fpm(16)==1) then
        write(fout,'(A)') '  Trapezoidal            |'
     elseif(fpm(16)==2) then
        write(fout,'(A)') '  Zolotarev              |'
     end if
  else
     write(fout,'(A)') '  User''s custom         |'
  endif

  if (fpm(16)/=2) then
     write(fout,'(A)',advance='no') '| Ellipse ratio y/x |'
     write(fout,'(F6.2)',advance='no') fpm(18)/100.0d0
     write(fout,'(A)') '                   |'
  endif

  write(fout,'(A)',advance='no') '| System solver     |'
  if (mod(fpm(30),10000)/1000==2) then !! ifeast
      if (fpm(30)/100000==1) then
           write(fout,'(A)') '  BiCGstab               |'
        else ! pfeast
           write(fout,'(A)') '  BiCGstab MPI           |'
        endif
     !write(fout,'(A)') '  BiCGstab               |'
     write(fout,'(A)',advance='no') '|                   |'
     write(fout,'(A)',advance='no') '  eps=1E-'
     if (fpm(45)<10) then
        write(fout,'(I1)',advance='no') fpm(45)
     else
        write(fout,'(I2)',advance='no') fpm(45)
     end if
     write(fout,'(A)',advance='no') '; maxit='
     write(fout,'(I4)',advance='no') fpm(46)
     write(fout,'(A)') '   |'
  else 
     select case (mod(fpm(30),100)/10) !! solver interface
     case(1) !! direct call to RCI
        write(fout,'(A)') '  User''s custom          |'
     case(2,5)
        write(fout,'(A)') '  LAPACK dense           |'
     case(3)
        write(fout,'(A)') '  SPIKE banded           |'
     case(4)
        if ((fpm(30)/100000==1).or.(fpm(22)==1))  then
           write(fout,'(A)') '  MKL-Pardiso            |'
        else ! pfeast
           write(fout,'(A)') '  MKL-Cluster-Pardiso    |'
        endif
     end select
  end if

  if (mod(fpm(30),100)/10/=1) then

     write(fout,'(A)',advance='no') '|                   |'           
     if (fpm(42)==1) then
        write(fout,'(A)')   '  Single precision       |'
     else
        write(fout,'(A)')   '  Double precision       |'
     endif

if (mod(fpm(30),100)/10==4) then     
     write(fout,'(A)',advance='no') '|                   |'          
     if (fpm(41)==1) then
        write(fout,'(A)')   '  Matrix scaled          |'
     else
        write(fout,'(A)')   '  Matrix not scaled      |'
     endif
  end if
     
  end if

  write(fout,'(A)',advance='no') '| FEAST uses MKL?   |'             
#ifdef MKL             
  write(fout,'(A)') '  Yes                    |'
#else
  write(fout,'(A)') '  No                     |'
#endif

  write(fout,'(A)',advance='no') '| Fact. stored?     |'  

  if (fpm(10)==1) then
     write(fout,'(A)')    '  Yes                    |'
  else
     write(fout,'(A)') '  No                     |'     
  endif
  write(fout,'(A)',advance='no') '| Initial Guess     |'  

  if (fpm(5)==1) then
     write(fout,'(A)')    '  User''s custom          |'
  else
     write(fout,'(A)') '  Random                 |'     
  endif

  write(fout,'(A)',advance='no') '| Size system       |'
  write(fout,'(I7)',advance='no') N
  write(fout,'(A)')    '                  |'
  write(fout,'(A)',advance='no') '| Size subspace     |'
  write(fout,'(I7)',advance='no') M0
  write(fout,'(A)')    '                  |'
  write(fout,'(A)') '-----------------------------------------------'




  write(fout,*)




end subroutine dfeast_info




subroutine zfeast_info(fout,fpm,Emid,r,N,M0)
  !!! Eric Polizzi 2019
  implicit none
  integer(8) :: fout
  integer,dimension(*) :: fpm
  complex(kind=kind(1.0d0)) :: Emid
  double precision :: r
  integer :: N,M0
!!!!!!!!!!!!!

  write(fout,*)
  write(fout,'(A)') '.-------------------.'
  write(fout,'(A)') '| FEAST data        |'
  write(fout,'(A)') '--------------------.-------------------------.'
  write(fout,'(A)',advance='no') '| Emid              |'
  write(fout,'(ES24.16)',advance='no') dble(Emid)
  write(fout,'(A)') '+|'
  write(fout,'(A)',advance='no') '|                   |'
  write(fout,'(ES24.16)',advance='no') aimag(Emid)
  write(fout,'(A)') 'i|'
  write(fout,'(A)',advance='no') '| Radius            |'
  write(fout,'(ES24.16)',advance='no') r
  write(fout,'(A)') ' |'
  write(fout,'(A)',advance='no') '| #Contour nodes    | '
  write(fout,'(I3)',advance='no') fpm(8)
  write(fout,'(A)') '   (full-contour)    |'
  write(fout,'(A)',advance='no') '| Quadrature rule   |'
  if (fpm(29)==1) then
     if (fpm(16)==0) then
        write(fout,'(A)') '  Gauss                  |'
     elseif(fpm(16)==1) then
        write(fout,'(A)') '  Trapezoidal            |'
     end if
  else
     write(fout,'(A)') '  User''s custom         |'
  endif
  write(fout,'(A)',advance='no') '| Ellipse ratio y/x |'
  write(fout,'(F6.2)',advance='no') fpm(18)/100.0d0
  write(fout,'(A)') '                   |'
  write(fout,'(A)',advance='no') '| Ellipse angle     |'
  write(fout,'(I4)',advance='no') fpm(19)
  write(fout,'(A)') '                     |'
  write(fout,'(A)',advance='no') '| System solver     |'
  if (mod(fpm(30),10000)/1000==2) then !! ifeast
      if (fpm(30)/100000==1) then
           write(fout,'(A)') '  BiCGstab               |'
        else ! pfeast
           write(fout,'(A)') '  BiCGstab MPI           |'
        endif
     !write(fout,'(A)') '  BiCGstab               |'
     write(fout,'(A)',advance='no') '|                   |'
     write(fout,'(A)',advance='no') '  eps=1E-'
     if (fpm(45)<10) then
        write(fout,'(I1)',advance='no') fpm(45)
     else
        write(fout,'(I2)',advance='no') fpm(45)
     end if
     write(fout,'(A)',advance='no') '; maxit='
     write(fout,'(I4)',advance='no') fpm(46)
     write(fout,'(A)') '   |'
  else 
     select case (mod(fpm(30),100)/10) !! solver interface
     case(1) !! direct call to RCI
        write(fout,'(A)') '  User''s custom          |'
     case(2,5)
        write(fout,'(A)') '  LAPACK dense           |'
     case(3)
        write(fout,'(A)') '  SPIKE banded           |'
     case(4)
        if ((fpm(30)/100000==1).or.(fpm(22)==1)) then
           write(fout,'(A)') '  MKL-Pardiso            |'
        else ! pfeast
           write(fout,'(A)') '  MKL-Cluster-Pardiso MPI|'
        endif
     end select
  end if

  if (mod(fpm(30),100)/10/=1) then

     write(fout,'(A)',advance='no') '|                   |'           
     if (fpm(42)==1) then
        write(fout,'(A)')   '  Single precision       |'
     else
        write(fout,'(A)')   '  Double precision       |'
     endif
if (mod(fpm(30),100)/10==4) then
     write(fout,'(A)',advance='no') '|                   |'          
     if (fpm(41)==1) then
        write(fout,'(A)')   '  Matrix scaled          |'
     else
        write(fout,'(A)')   '  Matrix not scaled      |'
     endif
  end if

  end if

  write(fout,'(A)',advance='no') '| FEAST uses MKL?   |'             
#ifdef MKL             
  write(fout,'(A)') '  Yes                    |'
#else
  write(fout,'(A)') '  No                     |'
#endif

  write(fout,'(A)',advance='no') '| Fact. stored?     |'  

  if (fpm(10)==1) then
     write(fout,'(A)')    '  Yes                    |'
  else
     write(fout,'(A)') '  No                     |'     
  endif

  write(fout,'(A)',advance='no') '| Initial Guess     |'  

  if (fpm(5)==1) then
     write(fout,'(A)')    '  User''s custom          |'
  else
     write(fout,'(A)') '  Random                 |'     
  endif

 !if (mod(fpm(30),1000)/100==3) then ! S, H or G (G here)                 
     write(fout,'(A)',advance='no') '| Eigenvectors      |'  
     if (fpm(15)==0) then
        write(fout,'(A)')    '  Right and Left         |'
     else
        write(fout,'(A)') '  Right only             |'     
     endif
 !end if
  write(fout,'(A)',advance='no') '| Size system       |'
  write(fout,'(I7)',advance='no') N
  write(fout,'(A)')    '                  |'
  write(fout,'(A)',advance='no') '| Size subspace     |'
  write(fout,'(I7)',advance='no') M0
  write(fout,'(A)')    '                  |'
  write(fout,'(A)') '-----------------------------------------------'




  write(fout,*)




end subroutine zfeast_info




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  
