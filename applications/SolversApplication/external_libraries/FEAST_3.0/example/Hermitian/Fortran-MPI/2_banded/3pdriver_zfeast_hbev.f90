!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Driver banded example !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! solving Ax=ex with A complex Hermitian --- A banded matrix!!!!!!!
!!!!!!! Using two search intervals
!!!!!!! by Eric Polizzi- 2009-2012!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program pdriver 

  implicit none

include 'mpif.h'

!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: feastparam 
  double precision :: epsout
  integer :: loop
  character(len=1) :: UPLO='F';  ! 'L' or 'U' also fine

!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  integer :: n,nnz,kl,ku
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: A

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,j,k
  integer :: M0,M,info
  double precision :: Emin,Emax,rea,img
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: X ! eigenvectors
  double precision,dimension(:),allocatable :: E,res ! eigenvalue+residual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI!!!!!!!!!!!!!!!!!!!!
integer :: code,rank,lrank,nb_procs,lnb_procs,color,key,NEW_COMM_WORLD
  call MPI_INIT(code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! read input file in coordinate format!!!!!!!
!!!!!!!!!!!!!!!!transform  in banded format directly !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='../../system2'

  !! find kl,ku (kl==ku here)
  kl=0
  open(10,file=trim(name),status='old')
  read(10,*) n,nnz
  do k=1,nnz
     read(10,*) i,j
     if (abs(i-j)>kl) kl=abs(i-j)
  end do
  close(10)
  ku=kl

  !! form the banded matrix
  allocate(A(1+kl+ku,N))
  A=(0.0d0,0.0d0)
  open(10,file=trim(name),status='old')
  read(10,*) N,nnz
  do k=1,nnz
     read(10,*) i,j,rea,img
  A(ku+1+(i-j),j)=rea*(1.0d0,0.0d0)+img*(0.0d0,1.0d0)
  end do
  close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (rank==0) then
  print *,'banded matrix -system2- size',n
  print *,'bandwidth',kl+ku+1
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in banded format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Definition of the two intervals 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
if (rank<=nb_procs/2-1) then
color=1 ! first interval
else
color=2 ! second interval
endif 

!!!!!!!!!!!!!!!!! create new_mpi_comm_world
key=0
call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,NEW_COMM_WORLD,code)
call MPI_COMM_RANK(NEW_COMM_WORLD,lrank,code) ! local rank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! search interval [Emin,Emax] including M eigenpairs
if (color==1) then !! 1st interval
 Emin=-0.35d0
 Emax= 0.0d0 
 M0=40
elseif(color==2) then !! 2nd interval
 Emin= 0.0d0
 Emax= 0.23d0 
 M0=40
endif

!!!!!!!!!!!!!! RUN INTERVALS in PARALLEL

  call system_clock(t1,tim)

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(e(1:M0))     ! Eigenvalue
  allocate(X(1:n,1:M0)) ! Eigenvectors
  allocate(res(1:M0))   ! Residual 

!!!!!!!!!!!!!  FEAST
  call feastinit(feastparam)

  feastparam(9)=NEW_COMM_WORLD
call zfeast_hbev(UPLO,N,kl,A,kl+ku+1,feastparam,epsout,loop,Emin,Emax,M0,E,X,M,res,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(t2,tim)

if (lrank==0) then
print *,'interval #',color
  print *,'FEAST OUTPUT INFO',info
  if (info==0) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'
     call MPI_COMM_SIZE(NEW_COMM_WORLD,lnb_procs,code)
     print *,'# of processors',lnb_procs
     print *,'SIMULATION TIME',(t2-t1)*1.0d0/tim
     print *,'# Search interval [Emin,Emax]',Emin,Emax
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),res(i)
     enddo
  endif
endif


call MPI_FINALIZE(code)

end program pdriver



