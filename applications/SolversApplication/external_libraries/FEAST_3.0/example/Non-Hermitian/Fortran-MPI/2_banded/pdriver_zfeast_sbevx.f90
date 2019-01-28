!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Expert Driver Example - Banded Storage
!!!!!!! solving Ax=ex with A complex-symmetric (non-Hermitian)
!!!!!!! James Kestyn, Eric Polizzi 2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program pdriver

  implicit none
 
  include 'mpif.h'

!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: feastparam 
  integer :: loop
  character :: UPLO='F'
!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  integer :: n,nnz,kl,ku,LDA
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: A

!!!!!!!!!!!!!!!!! Custom Contour
  integer :: ccN
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: Zne, Wne, Zedge
  integer, dimension(:), allocatable :: Nedge, Tedge 

!!!!!!!!!!!!!!!!! FEAST
  integer :: M0,M,info
  complex(kind=kind(1.0d0)) :: Emid
  double precision :: r,epsout
  double precision,dimension(:),allocatable :: resr ! residual
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: XR ! eigenvectors
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: E ! eigenvalues

!!!!!!!!!!!!!!!!! Others
  double precision :: rea,img
  integer :: t1,t2,tim
  integer :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI!!!!!!!!!!!!!!!!!!!!
integer :: code,rank,nb_procs
  call MPI_INIT(code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! read input file in coordinate format!!!!!!!
!!!!!!!!!!!!!!!!transform  in banded format directly !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='../../system4'

  !! find kl (ku=kl)
  kl=0
  open(10,file=trim(name),status='old')
  read(10,*) n,nnz
  do k=1,nnz
     read(10,*) i,j
     if ((i-j)>kl) kl=(i-j)
  end do
  close(10)
  
  !! form the banded matrix A 
  allocate(A(2*kl+1,N))
  LDA=2*kl+1
  A=0.0d0
  open(10,file=trim(name),status='old')
  read(10,*) N,nnz
  do k=1,nnz
     read(10,*) i,j,rea,img
     A(kl+1+(i-j),j) = rea*(1.0d0,0.0d0) + img*(0.0d0,1.0d0)
  end do
  close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(rank==0) print *,'banded matrix -system4- size',n
  if(rank==0) print *,'bandwidth',kl


  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in banded format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call feastinit(feastparam)
  feastparam(1)=1
  feastparam(6)=1
  M0=50 !! M0>=M

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create Custom Contour

  ccN = 3     !! number of pieces that make up contour
  allocate(Zedge(1:ccN))
  allocate(Nedge(1:ccN))
  allocate(Tedge(1:ccN))

  !!! Example contour - triangle  
  Zedge = (/(0.10d0,0.410d0),(4.2d0,0.41d0),(4.2d0,-8.3d0)/)
  Tedge(:) = (/0,0,0/)
  Nedge(:) = (/6,6,18/)

  !! Note: user must specify total # of contour points and edit feastparam(8)
  feastparam(8) = sum(Nedge(1:ccN))
  allocate(Zne(1:feastparam(8))) !! Contains the complex valued contour points 
  allocate(Wne(1:feastparam(8))) !! Contains the complex valued integrations weights

  !! Fill Zne/Wne
  call zfeast_customcontour(feastparam(8),ccN,Nedge,Tedge,Zedge,Zne,Wne)

if(rank==0)then
  print *,'---- Printing Countour ----'
  do i=1,feastparam(8)
     write(*,*)i,dble(Zne(i)),aimag(Zne(i))
  enddo
endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(E(1:M0))     ! Eigenvalue
  allocate(XR(1:n,1:M0)) ! Eigenvectors
  allocate(resr(1:M0))   ! Residual (if needed)

!!!!!!!!!!!!!  FEAST
 
  if(rank==0) print *, 'Enter FEAST'
  call zfeast_sbevx(UPLO,N,kl,A,LDA,feastparam,epsout,loop,Emid,r,M0,E,XR,M,resr,info,Zne,Wne)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(t2,tim)
  if(rank==0) print *,'FEAST OUTPUT INFO',info
  if (info==0 .and. rank==0) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'
     print *,'SIMULATION TIME',(t2-t1)*1.0d0/tim
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),resr(i)
     enddo
  endif

  call MPI_FINALIZE(code)

end program pdriver



