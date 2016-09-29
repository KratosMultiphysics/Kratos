!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Driver Example - CSR Storage
!!!!!!! solving Ax=eBx with A and B real non-symmetric (non-Hermitian) 
!!!!!!! James Kestyn, Eric Polizzi 2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program driver

  implicit none
!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: feastparam 
  integer :: loop
!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  integer :: n,nnz
  double precision,dimension(:),allocatable :: sa,sb
  integer,dimension(:),allocatable :: isa,jsa

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,k
  integer :: M0,M,info
  complex(kind=(kind(1.0d0))) :: Emid
  double precision :: r,epsout
  complex(kind=(kind(1.0d0))),dimension(:),allocatable :: E ! eigenvectors
  complex(kind=(kind(1.0d0))),dimension(:,:),allocatable :: X ! eigenvectors
  double precision,dimension(:),allocatable :: res ! eigenvalue+residual


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! read input file in csr format!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='../../system3'

!!!!!!!!!! form CSR arrays isa,jsa,sa 
  open(10,file=trim(name),status='old')
  read(10,*) n,nnz
  allocate(isa(1:N+1))
  isa=0
  isa(1)=1
  allocate(jsa(1:nnz))
  allocate(sa(1:nnz))
  allocate(sb(1:nnz))
  do k=1,nnz
     read(10,*) i,jsa(k),sa(k),sb(k)
     isa(i+1)=isa(i+1)+1
  end do
  close(10)
  do i=2,n+1
     isa(i)=isa(i)+isa(i-1)
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print *,'sparse matrix -system3- size',n
  print *,'nnz',nnz


  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in sparse format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! search interval [Emin,Emax] including M eigenpairs
  Emid=(0.59d0,0.0d0)
  r=0.4100d0 
  M0=40 !! M0>=M

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(e(1:M0))     ! Eigenvalue
  allocate(X(1:n,1:2*M0)) ! Eigenvectors
  allocate(res(1:2*M0))   ! Residual (if needed)

!!!!!!!!!!!!!  FEAST
  call feastinit(feastparam)
  feastparam(1)=1
  call dfeast_gcsrgv(N,sa,isa,jsa,sb,isa,jsa,feastparam,epsout,loop,Emid,r,M0,E,X,M,res,info)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(t2,tim)
  print *,'FEAST OUTPUT INFO',info
  if (info==0) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'
     print *,'SIMULATION TIME',(t2-t1)*1.0d0/tim
     print *,'# Search interval [Emin,Emax]',Emid,r
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),max(res(i),res(M0+i))
     enddo
  endif

end program driver



