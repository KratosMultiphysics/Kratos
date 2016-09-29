!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Expert Driver - Dense Storage
!!!!!!! solving Ax=ex with A complex-symmetric (non-Hermitian)
!!!!!!! James Kestyn, Eric Polizzi 2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program pdriver 

  implicit none

include 'mpif.h'

!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: feastparam 
  integer :: loop
  character(len=1) :: UPLO='F' ! 'L' or 'U' also fine

!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  integer :: n,nnz
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: A

!!!!!!!!!!!!!!!!! FEAST
  integer :: M0,M,info
  double precision :: r,epsout
  complex(kind=kind(1.0d0)) :: Emid 
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: XR ! eigenvectors
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: E ! eigenvalues
  double precision,dimension(:),allocatable :: resr ! residual

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,j,k
  double precision :: rea,img
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI!!!!!!!!!!!!!!!!!!!!
integer :: code,rank,nb_procs
  call MPI_INIT(code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Read Coordinate format and convert to dense format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='../../system4'


  open(10,file=trim(name),status='old')
  read(10,*) n,nnz
  allocate(A(1:n,1:n))
  A=0.0d0
  do k=1,nnz
     read(10,*) i,j,rea,img
     A(i,j)=rea*(1.0d0,0.0d0)+img*(0.0d0,1.0d0)
  enddo
  close(10)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(rank==0) print *,'dense matrix -system4- size',n

  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call feastinit(feastparam)
  feastparam(1)=1
  M0=50 !! M0>=M

  Emid = (4.0d0,0.00d0)
  r = 3.0e0

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(E(1:M0))     ! Eigenvalue
  allocate(XR(1:n,1:M0)) ! Right Eigenvectors ( XL = CONJG(XR) )
  allocate(resr(1:M0))   ! Residual (if needed)

!!!!!!!!!!!!  FEAST
  print *, 'Enter FEAST',rank
  call zfeast_syev(UPLO,N,A,N,feastparam,epsout,loop,Emid,r,M0,E,XR,M,resr,info)

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
        print *,i,dble(E(i)),aimag(E(i)),resr(i)
     enddo
  endif

  call MPI_FINALIZE(code)


end program pdriver 



