program helloworld
  implicit none
  !! eigenvalue system 
  integer,parameter :: N=2,LDA=2
  character(len=1),parameter :: UPLO='F' ! 'L' or 'U' also fine
  double precision,dimension(N*N) :: A=(/2.0d0,-1.0d0,-1.0d0,2.0d0/)
  double precision :: Emin=-5.0d0, Emax=5.0d0
  integer :: M0=2 ! (Initial) subspace dimension
  !! input parameters for FEAST
  integer,dimension(64) :: feastparam
  !! output variables for FEAST
  double precision,dimension(:),allocatable :: E, res
  double precision,dimension(:,:),allocatable :: X
  double precision :: epsout
  integer :: i,loop,info,M

!!! Allocate memory for eigenvalues.eigenvectors/residual
  allocate(E(M0))
  allocate(res(M0))
  allocate(X(N,M0))

!!!!!!!!!! FEAST
  call feastinit(feastparam)
  feastparam(1)=1 !! change from default value
  call dfeast_syev(UPLO,N,A,LDA,feastparam,epsout,loop,Emin,Emax,M0,E,X,M,res,info)

!!!!!!!!!! REPORT
  print *,'FEAST OUTPUT INFO',info
  if (info==0) then
     print *,'*************************************************'
     print *,'************** REPORT ***************************'
     print *,'*************************************************'
     print *,'# Search interval [Emin,Emax]',Emin,Emax
     print *,'# mode found/subspace',M,M0
     print *,'# iterations',loop
     print *,'TRACE',sum(E(1:M))
     print *,'Relative error on the Trace',epsout
     print *,'Eigenvalues/Residuals'
     do i=1,M
        print *,i,E(i),res(i)
     enddo
     print *,'Eigenvectors'
     do i=1,M
        print *,i,"(",X(1,i),X(2,i),")"
     enddo
  endif

end program helloworld
