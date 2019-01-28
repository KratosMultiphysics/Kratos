!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST Driver sparse example !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! solving Ax=eBx with A real and B spd --- A and B sparse matrix!!!
!!!!!!! by Eric Polizzi- 2009-2011!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program driver

  implicit none
!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: feastparam 
  double precision :: epsout
  integer :: loop
  character(len=1) :: UPLO='F'
!!!!!!!!!!!!!!!!! Matrix declaration variable
  character(len=100) :: name
  integer :: n,nnz
  double precision,dimension(:),allocatable :: sa,sb
  integer,dimension(:),allocatable :: isa,jsa

!!!!!!!!!!!!!!!!! Others
  integer :: t1,t2,tim
  integer :: i,k
  integer :: M0,M,info
  double precision :: Emin,Emax
  double precision,dimension(:,:),allocatable :: X ! eigenvectors
  double precision,dimension(:),allocatable :: E,res ! eigenvalue+residual


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! read input file in csr format!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  name='../../system1'

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
  print *,'sparse matrix -system1- size',n
  print *,'nnz',nnz


  call system_clock(t1,tim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST in sparse format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! search interval [Emin,Emax] including M eigenpairs
  Emin= 0.18d0
  Emax= 1.00d0 
  M0=25 !! M0>=M

!!!!!!!!!!!!! ALLOCATE VARIABLE 
  allocate(e(1:M0))     ! Eigenvalue
  allocate(X(1:n,1:M0)) ! Eigenvectors
  allocate(res(1:M0))   ! Residual 

!!!!!!!!!!!!!  FEAST
  call feastinit(feastparam)
  feastparam(1)=1
  call dfeast_scsrgv(UPLO,N,sa,isa,jsa,sb,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,E,X,M,res,info)

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

end program driver



