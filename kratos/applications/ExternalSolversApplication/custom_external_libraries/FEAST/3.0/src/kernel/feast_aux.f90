!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Auxiliary Routines for FEAST (undocumented)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!! List of routines:
! --------------------

!!!!!!!! Check Inputs Routines
! check_feast_fpm_input               ! Check validity of FEAST input parameters fpm[1-19]  
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





!!!!!!!! Integration nodes
! sset_feast_gauss_legendre
! dset_feast_gauss_legendre
! sset_feast_zolotarev
! dset_feast_zolotarev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine check_feast_fpm_input(fp,info)
  !  Purpose 
  !  =======
  !  Error handling for input FEAST parameters.
  !  Check the values for the input FEAST parameters, and return 
  !  info code error /=0 if incorrect values are found
  !
  !  Arguments
  !  =========
  !
  !  fp   (input) INTEGER(*) : FEAST parameters
  !  info (input/output) INTEGER
  !=====================================================================
  ! Eric Polizzi 2009-2015
  ! ====================================================================
  implicit none
  integer,dimension(*) :: fp
  integer :: info
  integer:: i
  logical :: test
  integer,parameter :: max=5
  integer, dimension(max):: tnbe=(/24,32,40,48,56/)
  

  if ((fp(1)/=0).and.(fp(1)/=1)) info=101

  test=.false.
  !if ((((fp(16)==2).and.(fp(17)==0)).or.(fp(16)==0).and.(fp(17)==0))) then ! only  Gauss and Zolotarev
  if ((fp(16)==0).or.(fp(16)==2)) then ! Gauss or Zolotarev for Hermitian (#contour points - half-contour)
     if (fp(2)>20) then
        test=.true.
        do i=1,max
           if (fp(2)==tnbe(i)) test=.false. 
        enddo
     end if
  end if

  if (test) info=102
  if (fp(2)<1) info=102

  if ((fp(3)<0).or.(fp(3)>16)) info=103
  if (fp(4)<0) info=104
  if ((fp(5)/=0).and.(fp(5)/=1)) info=105
  if ((abs(fp(6))>1)) info=106
  if ((fp(7)<0).or.(fp(7)>7)) info=107

  test=.false.
  if (fp(17)==0) then ! Gauss  for Non-Hermitian (#contour points - full-contour==2 half-Gauss contour)
     if (fp(8)>40) then
        test=.true.
        do i=1,max
           if (fp(8)==2*tnbe(i)) test=.false. 
        enddo
     end if
  end if

  if (test) info=108
  if (fp(8)<2) info=108

  if ((fp(12)/=0).and.(fp(12)/=1)) info=112
  if ((fp(14)<0).or.(fp(14)>2)) info=114
  if ((fp(16)<-1).or.(fp(16)>2)) info=116
  if ((fp(17)<-1).or.(fp(17)>1)) info=117

  if (fp(18)<0) info=118
  if (fp(19)<0) info=119
end subroutine check_feast_fpm_input





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
  if ((M0<=0).or.(M0>N)) info=201 ! problem with M0 
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
  include "f90_noruntime_interface.fi"
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
        nE = (eig(j)-Emid)*(ONEC*wdcos(theta)+(DZERO,DONE)*wdsin(theta))+Emid
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
  include "f90_noruntime_interface.fi"
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
        nE = (eig(j)-Emid)*(ONEC*wscos(theta)+(SZERO,SONE)*wssin(theta))+Emid
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








subroutine sset_feast_gauss_legendre(nbe,e,xe,we)
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


  call  dset_feast_gauss_legendre(nbe,e,dxe,dwe)

  xe=real(dxe)
  we=real(dwe)

end subroutine sset_feast_gauss_legendre



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dset_feast_gauss_legendre(nbe,e,xe,we)
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


end subroutine dset_feast_gauss_legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine sset_feast_zolotarev(nbe,e,xe,we)
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

  call dset_feast_zolotarev(nbe,e,zxe,zwe)

  xe=cmplx(zxe)
  we=cmplx(zwe)

end subroutine sset_feast_zolotarev








subroutine dset_feast_zolotarev(nbe,e,xe,we)
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



end subroutine dset_feast_zolotarev


