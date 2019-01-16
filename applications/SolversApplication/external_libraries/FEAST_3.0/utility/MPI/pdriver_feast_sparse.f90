!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! FEAST general Driver sparse !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! solving Ax=eBx or Ax=eX      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! symmetric, hermitian or general matrices !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! single or double precision !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! where A (and B if any), provided by user in coordinate format !!!!                         
!!!!!!! by Eric Polizzi- 2009-2015  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program pdriver_feast_sparse

  implicit none
  include "mpif.h"

!!!!!!!!!!!!!!!!! Feast declaration variable
  integer,dimension(64) :: fpm
  double precision :: depsout
  real :: sepsout
  integer :: loop

  character(len=100) :: name
  character(len=1) :: UPLO,PRE,SHG,EG 
  integer :: n,nnza,nnzb
  real,dimension(:),allocatable :: ssa,ssb,sca,scb
  double precision,dimension(:),allocatable :: dsa,dsb,dca,dcb
  complex,dimension(:),allocatable :: csa,csb,cca,ccb
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: zsa,zsb,zca,zcb
  integer,dimension(:),allocatable :: isa,jsa,isb,jsb,ica,jca,icb,jcb


  double precision :: t1,t2
  integer :: i,j,k,pc
  integer :: M0,M,info
  character(len=1) :: cc


  double precision :: dEmin,dEmax,dr
  complex(kind=kind(1.0d0)):: zEmid
  complex:: cEmid
  real :: sEmin,sEmax,sr
  double precision:: drea,dimg
  real:: srea,simg

  double precision,dimension(:,:),allocatable :: dX
  real,dimension(:,:),allocatable :: sX
  complex,dimension(:,:),allocatable :: cX
  complex(kind=kind(1.0d0)),dimension(:,:),allocatable :: zX

  double precision,dimension(:),allocatable :: dres
  real,dimension(:),allocatable :: sres

  double precision,dimension(:),allocatable :: dE
  real,dimension(:),allocatable :: sE
  complex,dimension(:),allocatable :: cE
  complex(kind=kind(1.0d0)),dimension(:),allocatable :: zE
  integer :: code,rank,nb_procs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPI!!!!!!!!!!!!!!!!!!!!
  call MPI_INIT(code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! read main input file !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1,name)



  call feastinit(fpm)


!!!!!!!!!!!! DRIVER_FEAST_SPARSE input file  
  open(10,file=trim(name)//'.in',status='old')
  read(10,*) SHG ! type of eigenvalue problem "General, Hermitian, Symmetric" 
  read(10,*) EG ! type of eigenvalue probl  g== sparse generalized, e== sparse standard
  read(10,*) PRE  ! "PRE"==(s,d,c,z) resp. (single real,double real,complex,double complex) 
  read(10,*) UPLO ! UPLO==(F,L,U) reps. (Full csr, Lower csr, Upper csr) 

  if (SHG=='s') then

     if (PRE=='d') then
        read(10,*) dEmin
        read(10,*) dEmax
     elseif (PRE=='s') then
        read(10,*) sEmin
        read(10,*) sEmax
     elseif (PRE=='z') then
        read(10,*) drea,dimg
        zEmid=drea*(1.0d0,0.0d0)+dimg*(0.0d0,1.0d0)
        read(10,*) dr
     elseif (PRE=='c') then
        read(10,*) srea,simg
        cEmid=srea*(1.0e0,0.0e0)+simg*(0.0e0,1.0e0)
        read(10,*) sr
     end if


  elseif (SHG=='h') then

     if (PRE=='z') then
        read(10,*) dEmin
        read(10,*) dEmax
     else
        read(10,*) sEmin
        read(10,*) sEmax
     end if

  elseif (SHG=='g') then

     if ((PRE=='d').or.(PRE=='z')) then
        read(10,*) drea,dimg
        zEmid=drea*(1.0d0,0.0d0)+dimg*(0.0d0,1.0d0)
        read(10,*) dr
     else
        read(10,*) srea,simg
        cEmid=srea*(1.0e0,0.0e0)+simg*(0.0e0,1.0e0)
        read(10,*) sr
     end if
  end if

  read(10,*) M0   ! size subspace

  read(10,*) pc ! Some changes from default for fpm
  do i=1,pc
     read(10,*) j,fpm(j)
  enddo

  close(10)


!!!!!!!!!!!read matrix A

  open(10,file=trim(name)//'_A.mtx',status='old')
  k=0
  cc='%'
  do while(cc=='%')
     k=k+1 
     read(10,'(A1)') cc
  end do
  close(10)

  open(10,file=trim(name)//'_A.mtx',status='old')
  do i=1,k-1
     read(10,'(A1)') cc
  enddo
  read(10,*) n,n,nnza
  allocate(ica(nnza))
  allocate(jca(nnza))
  if (PRE=='s') then
     allocate(sca(nnza))
     do i=1,nnza
        read(10,*) ica(i),jca(i),sca(i)
     end do
  elseif (PRE=='d') then
     allocate(dca(nnza))
     do i=1,nnza
        read(10,*) ica(i),jca(i),dca(i)
     end do
  elseif (PRE=='c') then
     allocate(cca(nnza))
     do i=1,nnza
        read(10,*) ica(i),jca(i),srea,simg 
        cca(i)=srea*(1.0e0,0.0e0)+simg*(0.0e0,1.0e0)
     end do
  elseif (PRE=='z') then
     allocate(zca(nnza))
     do i=1,nnza
        read(10,*) ica(i),jca(i),drea,dimg 
        zca(i)=drea*(1.0d0,0.0d0)+dimg*(0.0d0,1.0d0)
     end do
  end if
  close(10)

  !! create csr format
  allocate(isa(1:n+1))
  allocate(jsa(1:nnza))


  if (PRE=='s') then
     allocate(ssa(1:nnza))
     call scoo2csr(n,nnza,ica,jca,sca,isa,jsa,ssa)
  elseif (PRE=='d') then
     allocate(dsa(1:nnza))
     call dcoo2csr(n,nnza,ica,jca,dca,isa,jsa,dsa)
  elseif (PRE=='c') then
     allocate(csa(1:nnza))
     call ccoo2csr(n,nnza,ica,jca,cca,isa,jsa,csa)
  elseif (PRE=='z') then
     allocate(zsa(1:nnza))
     call zcoo2csr(n,nnza,ica,jca,zca,isa,jsa,zsa)
  end if



!!!!!!!!!!!read matrix B if any
  if (EG=='g') then

     open(10,file=trim(name)//'_B.mtx',status='old')
     k=0
     cc='%'
     do while(cc=='%')
        k=k+1 
        read(10,'(A1)') cc
     end do
     close(10)

     open(10,file=trim(name)//'_B.mtx',status='old')
     do i=1,k-1
        read(10,'(A1)') cc
     enddo
     read(10,*) n,n,nnzb
     allocate(icb(nnzb))
     allocate(jcb(nnzb))
     if (PRE=='s') then
        allocate(scb(nnzb))
        do i=1,nnzb
           read(10,*) icb(i),jcb(i),scb(i)
        end do
     elseif (PRE=='d') then
        allocate(dcb(nnzb))
        do i=1,nnzb
           read(10,*) icb(i),jcb(i),dcb(i)
        end do
     elseif (PRE=='c') then
        allocate(ccb(nnzb))
        do i=1,nnzb
           read(10,*) icb(i),jcb(i),srea,simg 
           ccb(i)=srea*(1.0e0,0.0e0)+simg*(0.0e0,1.0e0)
        end do
     elseif (PRE=='z') then
        allocate(zcb(nnzb))
        do i=1,nnzb
           read(10,*) icb(i),jcb(i),drea,dimg 
           zcb(i)=drea*(1.0d0,0.0d0)+dimg*(0.0d0,1.0d0)
        end do
     end if
     close(10)

     !! create csr format
     allocate(isb(1:n+1))
     allocate(jsb(1:nnzb))


     if (PRE=='s') then
        allocate(ssb(1:nnzb))
        call scoo2csr(n,nnzb,icb,jcb,scb,isb,jsb,ssb)
     elseif (PRE=='d') then
        allocate(dsb(1:nnzb))
        call dcoo2csr(n,nnzb,icb,jcb,dcb,isb,jsb,dsb)
     elseif (PRE=='c') then
        allocate(csb(1:nnzb))
        call ccoo2csr(n,nnzb,icb,jcb,ccb,isb,jsb,csb)
     elseif (PRE=='z') then
        allocate(zsb(1:nnzb))
        call zcoo2csr(n,nnzb,icb,jcb,zcb,isb,jsb,zsb)
     end if


  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank==0) then
     print *,'matrix name ',trim(name)
     print *,'matrix -coordinate format- size',n
     print *,'sparse matrix A- nnz',nnza
     if (EG=='g') print *,'sparse matrix B- nnz',nnzb
     print *,''
  endif

  info=-1
  t1=MPI_WTIME()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! FEAST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!  FEAST SYMMETRIC

  if (SHG=='s') then 

     if ((PRE=='d').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','dfeast_scsrgv'
        allocate(dE(1:M0))     ! Eigenvalue
        allocate(dres(1:M0))   ! Residual 
        allocate(dX(1:n,1:M0))
        call dfeast_scsrgv(UPLO,N,dsa,isa,jsa,dsb,isb,jsb,fpm,depsout,loop,dEmin,dEmax,M0,dE,dX,M,dres,info)

     elseif  ((PRE=='d').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','dfeast_scsrev'
        allocate(dE(1:M0))     ! Eigenvalue
        allocate(dres(1:M0))   ! Residual 
        allocate(dX(1:n,1:M0))
        call dfeast_scsrev(UPLO,N,dsa,isa,jsa,fpm,depsout,loop,dEmin,dEmax,M0,dE,dX,M,dres,info)

     elseif ((PRE=='s').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','sfeast_scsrgv'
        allocate(sE(1:M0))     ! Eigenvalue
        allocate(sres(1:M0))   ! Residual 
        allocate(sX(1:n,1:M0))
        call sfeast_scsrgv(UPLO,N,ssa,isa,jsa,ssb,isb,jsb,fpm,sepsout,loop,sEmin,sEmax,M0,sE,sX,M,sres,info)

     elseif ((PRE=='s').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','sfeast_scsrev'
        allocate(sE(1:M0))     ! Eigenvalue
        allocate(sres(1:M0))   ! Residual 
        allocate(sX(1:n,1:M0))
        call sfeast_scsrev(UPLO,N,ssa,isa,jsa,fpm,sepsout,loop,sEmin,sEmax,M0,sE,sX,M,sres,info)

     elseif ((PRE=='z').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','zfeast_scsrgv'
        allocate(zE(1:M0))     ! Eigenvalue
        allocate(dres(1:M0))   ! Residual 
        allocate(zX(1:n,1:M0))
        call zfeast_scsrgv(UPLO,N,zsa,isa,jsa,zsb,isb,jsb,fpm,depsout,loop,zEmid,dr,M0,zE,zX,M,dres,info)

     elseif ((PRE=='z').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','zfeast_scsrev'
        allocate(zE(1:M0))     ! Eigenvalue
        allocate(dres(1:M0))   ! Residual 
        allocate(zX(1:n,1:M0))
        call zfeast_scsrev(UPLO,N,zsa,isa,jsa,fpm,depsout,loop,zEmid,dr,M0,zE,zX,M,dres,info)

     elseif ((PRE=='c').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','cfeast_scsrgv'
        allocate(cE(1:M0))     ! Eigenvalue
        allocate(sres(1:M0))   ! Residual 
        allocate(cX(1:n,1:M0))
        call cfeast_scsrgv(UPLO,N,csa,isa,jsa,csb,isb,jsb,fpm,sepsout,loop,cEmid,sr,M0,cE,cX,M,sres,info)

     elseif ((PRE=='c').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','cfeast_scsrev'
        allocate(cE(1:M0))     ! Eigenvalue
        allocate(sres(1:M0))   ! Residual 
        allocate(cX(1:n,1:M0))
        call cfeast_scsrev(UPLO,N,csa,isa,jsa,fpm,sepsout,loop,cEmid,sr,M0,cE,cX,M,sres,info)

     end if



!!!!!!!!!!!!!  FEAST HERMITIAN

  elseif (SHG=='h') then 


     if ((PRE=='z').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','zfeast_hcsrgv'
        allocate(dE(1:M0))     ! Eigenvalue
        allocate(dres(1:M0))   ! Residual 
        allocate(zX(1:n,1:M0))
        call zfeast_hcsrgv(UPLO,N,zsa,isa,jsa,zsb,isb,jsb,fpm,depsout,loop,dEmin,dEmax,M0,dE,zX,M,dres,info)

     elseif ((PRE=='z').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','zfeast_hcsrev'
        allocate(dE(1:M0))     ! Eigenvalue
        allocate(dres(1:M0))   ! Residual 
        allocate(zX(1:n,1:M0))
        call zfeast_hcsrev(UPLO,N,zsa,isa,jsa,fpm,depsout,loop,dEmin,dEmax,M0,dE,zX,M,dres,info)

     elseif ((PRE=='c').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','cfeast_hcsrgv'
        allocate(sE(1:M0))     ! Eigenvalue
        allocate(sres(1:M0))   ! Residual 
        allocate(cX(1:n,1:M0))
        call cfeast_hcsrgv(UPLO,N,csa,isa,jsa,csb,isb,jsb,fpm,sepsout,loop,sEmin,sEmax,M0,sE,cX,M,sres,info)

     elseif ((PRE=='c').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','cfeast_hcsrev'
        allocate(sE(1:M0))     ! Eigenvalue
        allocate(sres(1:M0))   ! Residual 
        allocate(cX(1:n,1:M0))
        call cfeast_hcsrev(UPLO,N,csa,isa,jsa,fpm,sepsout,loop,sEmin,sEmax,M0,sE,cX,M,sres,info)

     end if



!!!!!!!!!!!!!  FEAST GENERAL

  elseif (SHG=='g') then 


     if ((PRE=='d').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','dfeast_gcsrgv'
        allocate(zE(1:M0))     ! Eigenvalue
        allocate(dres(1:2*M0))   ! Residual 
        allocate(zX(1:n,1:2*M0))
        call dfeast_gcsrgv(N,dsa,isa,jsa,dsb,isb,jsb,fpm,depsout,loop,zEmid,dr,M0,zE,zX,M,dres,info)


     elseif  ((PRE=='d').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','dfeast_gcsrev'
        allocate(zE(1:M0))     ! Eigenvalue
        allocate(dres(1:2*M0))   ! Residual 
        allocate(zX(1:n,1:2*M0))
        call dfeast_gcsrev(N,dsa,isa,jsa,fpm,depsout,loop,zEmid,dr,M0,zE,zX,M,dres,info)

     elseif ((PRE=='s').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','sfeast_gcsrgv'
        allocate(cE(1:M0))     ! Eigenvalue
        allocate(sres(1:2*M0))   ! Residual 
        allocate(cX(1:n,1:2*M0))
        call sfeast_gcsrgv(N,ssa,isa,jsa,ssb,isb,jsb,fpm,sepsout,loop,cEmid,sr,M0,cE,cX,M,sres,info)

     elseif ((PRE=='s').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','sfeast_gcsrev'
        allocate(cE(1:M0))     ! Eigenvalue
        allocate(sres(1:2*M0))   ! Residual 
        allocate(cX(1:n,1:2*M0))
        call sfeast_gcsrev(N,ssa,isa,jsa,fpm,sepsout,loop,cEmid,sr,M0,cE,cX,M,sres,info)

     elseif ((PRE=='z').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','zfeast_gcsrgv'
        allocate(zE(1:M0))     ! Eigenvalue
        allocate(dres(1:2*M0))   ! Residual 
        allocate(zX(1:n,1:2*M0))
        call zfeast_gcsrgv(N,zsa,isa,jsa,zsb,isb,jsb,fpm,depsout,loop,zEmid,dr,M0,zE,zX,M,dres,info)

     elseif ((PRE=='z').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','zfeast_gcsrev'
        allocate(zE(1:M0))     ! Eigenvalue
        allocate(dres(1:2*M0))   ! Residual 
        allocate(zX(1:n,1:2*M0))
        call zfeast_gcsrev(N,zsa,isa,jsa,fpm,depsout,loop,zEmid,dr,M0,zE,zX,M,dres,info)

     elseif ((PRE=='c').and.(EG=='g')) then 
        if (rank==0) print *,'Routine  ','cfeast_gcsrgv'
        allocate(cE(1:M0))     ! Eigenvalue
        allocate(sres(1:2*M0))   ! Residual 
        allocate(cX(1:n,1:2*M0))
        call cfeast_gcsrgv(N,csa,isa,jsa,csb,isb,jsb,fpm,sepsout,loop,cEmid,sr,M0,cE,cX,M,sres,info)

     elseif ((PRE=='c').and.(EG=='e')) then 
        if (rank==0) print *,'Routine  ','cfeast_gcsrev'
        allocate(cE(1:M0))     ! Eigenvalue
        allocate(sres(1:2*M0))   ! Residual 
        allocate(cX(1:n,1:2*M0))
        call cfeast_gcsrev(N,csa,isa,jsa,fpm,sepsout,loop,cEmid,sr,M0,cE,cX,M,sres,info)

     end if

  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! POST-PROCESSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t2=MPI_WTIME()

  if (rank==0) then
     print *,'FEAST OUTPUT INFO',info

     IF ((info==0).or.(info==5).or.(info==6)) then
        print *,'*************************************************'
        print *,'************** REPORT ***************************'
        print *,'*************************************************'

        if (fpm(14)==0) then 
           print *,'Eigenvalues/Residuals'

           if (SHG/='g') then
              print *,'inside interval'
              do i=1,M0
                 if ((SHG=='s').and.(PRE=='z')) then
                    print *,i,zE(i),dres(i)
                 elseif ((SHG=='s').and.(PRE=='c')) then
                    print *,i,cE(i),sres(i)
                 elseif ((SHG=='s').and.(PRE=='d')) then
                    print *,i,dE(i),dres(i)
                 elseif ((SHG=='s').and.(PRE=='s')) then
                    print *,i,sE(i),sres(i)
                 elseif ((SHG=='h').and.(PRE=='z')) then
                    print *,i,dE(i),dres(i)
                 elseif ((SHG=='h').and.(PRE=='c')) then
                    print *,i,sE(i),sres(i)
                 end if
                 if (i==M) then 
print *,''
print *,'outside interval'
endif
              enddo
           else
              print *,'inside interval'
              do i=1,M0
                 if ((PRE=='s').or.(PRE=='c')) then
                    print *,i,cE(i),sres(i),sres(i+M0)
                 else
                    print *,i,zE(i),dres(i),dres(i+M0)
                 end if
                 if (i==M) then
print *,''
print *,'outside interval'
endif
              enddo

           end if

        elseif (fpm(14)==2) then
           print *,'Running average for estimation (1 -> M0)'
           do i=1,M0
              if ((PRE=='s').or.(PRE=='c')) then
                 print *,i,sres(i)
              else
                 print *,i,dres(i)
              end if
           enddo
        endif

        print *,'------------------------------------'
        print *,'SIMULATION TIME',(t2-t1)*1.0d0


        if ((SHG=='g').or.((SHG=='s').and.((PRE=='c').or.(PRE=='z')))) then
           print *,'# Search interval [Emid,r]'
           if ((PRE=='s').or.(PRE=='c')) then
              print *,cEmid,sr
           else
              print *,zEmid,dr
           end if
        else
           print *,'# Search interval [Emin,Emax]'
           if ((PRE=='s').or.(PRE=='c')) then
              print *,sEmin,sEmax
           else
              print *,dEmin,dEmax
           end if
        end if

        if (fpm(14)==2) then
           print *,'Subspace size M0        ',M0
           print *,'# estimated eigenvalue M',M

        else

           print *,'Subspace size M0    ',M0
           print *,'# eigenvalue found M',M
           print *,'# FEAST iterations  ',loop
           if ((PRE=='s').or.(PRE=='c')) then
              print *,'Relative error on the Trace  ',sepsout
              if (SHG/='g') then
                 print *,'Maximum eigenvector residual ',maxval(sres(1:M))
              else
                 print *,'Maximum eigenvector residuals',maxval(sres(1:M)),maxval(sres(M0+1:M0+M))
              end if
           else
              print *,'Relative error on the Trace  ',depsout
              if (SHG/='g') then
                 print *,'Maximum eigenvector residual ',maxval(dres(1:M))
              else
                 print *,'Maximum eigenvector residuals',maxval(dres(1:M)),maxval(dres(M0+1:M0+M))
              end if
           end if


        end if

     end IF
  end if


  call MPI_FINALIZE(code)


end program pdriver_feast_sparse









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine dcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  double precision,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  double precision :: adum


  isa(1:N+1) = 0  
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do



  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do


  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine dcoo2csr




subroutine zcoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  complex(kind=kind(1.0d0)),dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  complex(kind=kind(1.0d0)):: adum

  isa(1:N+1) = 0  
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do

  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do

  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine zcoo2csr



subroutine scoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  real,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  real :: adum


  isa(1:N+1) = 0  
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do

  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do

  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine scoo2csr




subroutine ccoo2csr(n,nnz,ic,jc,c,isa,jsa,sa)
  implicit none
  integer :: n,nnz
  integer,dimension(*) :: ic,jc,isa,jsa
  complex,dimension(*) :: c,sa
!!!
  integer :: k,k1,i,j,idum
  integer,dimension(n) :: iloc
  complex:: adum

  isa(1:N+1) = 0  
  !find how many elements in row
  do  k=1, nnz
     isa(ic(k)) = isa(ic(k))+1
  end do

  !build isa
  k = 1
  do  i=1,N+1
     k1 = isa(i)
     isa(i) = k
     k = k+k1
  end do

  iloc(1:n)=isa(1:n)
  !Build jsa, sa - increment local row counter
  do  k=1, nnz
     sa(iloc(ic(k))) =  c(k)
     jsa(iloc(ic(k))) = jc(k)
     iloc(ic(k)) = iloc(ic(k))+1
  end do
  ! Reorder by increasing column
  do i=1,n
     do k=isa(i),isa(i+1)-1
        do k1=k,isa(i+1)-1
           if (jsa(k1)<jsa(k)) then
              idum=jsa(k)
              jsa(k)=jsa(k1)
              jsa(k1)=idum
              adum=sa(k)
              sa(k)=sa(k1)
              sa(k1)=adum
           endif
        enddo
     enddo
  enddo


end subroutine ccoo2csr


