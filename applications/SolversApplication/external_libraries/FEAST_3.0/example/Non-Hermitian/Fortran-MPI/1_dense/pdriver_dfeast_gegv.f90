!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Driver - Dense Storage 
!!!!!!! solving Ax=eBx with A real non-symmetric (non-Hermitian)
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
  double precision,dimension(:,:),allocatable :: A,B,Ab,Bb

!!!!!!!!!!!!!!!!! FEAST
  complex(kind=kind(1.0d0)) :: Emid
  double precision :: r,epsout
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: X    ! eigenvectors
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: E ! eigenvalues
  double precision,dimension(:),allocatable :: res  ! residual

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,j,k
  integer :: M0,M,info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI!!!!!!!!!!!!!!!!!!!!
  integer :: code,rank,nb_procs
  call MPI_INIT(code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Read Coordinate format and convert to dense format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='../../system3'
  open(10,file=trim(name),status='old')
  read(10,*) n,nnz
  allocate(A(1:n,1:n))
  allocate(B(1:n,1:n))
  A=0.0d0
  B=0.0d0
  do k=1,nnz
     read(10,*) i,j,A(i,j),B(i,j)
  enddo
  close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(rank==0) print *,'dense matrix -system3- size',n

  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Emid=(0.59d0,0.0d0)
  r=0.4100d0 
  M0=30!400!25 !! M0>=M

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(E(1:M0))      ! Eigenvalue
  allocate(X(1:n,1:2*M0)) ! Eigenvectors
  allocate(res(1:2*M0))   ! Residual (if needed)

!!!!!!!!!!!!  FEAST
  call feastinit(feastparam)
  feastparam(1)=1
  print *, 'Enter Feast',rank
  call dfeast_gegv(N,A,N,B,N,feastparam,epsout,loop,Emid,r,M0,E,X,M,res,info)

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
     print *,'# Search interval [Emid,r]',Emid,r
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),max(res(i),res(M0+i))
     enddo
  endif

  call MPI_FINALIZE(code)

end program pdriver



