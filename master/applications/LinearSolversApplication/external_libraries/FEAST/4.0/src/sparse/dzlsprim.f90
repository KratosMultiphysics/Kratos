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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! SPARSE PRIMITIVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

!!!!! List of routines 


!daddcsr, zaddcsr           :: add two square csr matrices and return a new one (C=alpha*A+beta*B)
!zdaddcsr                   :: same as zaddcsr but A, B are real;  alpha, beta, C are complex
!zinc_addcsr                :: Increment matrix B by adding matrix A times alpha (sparsity pattern of A included in B)
  

!dcsr2csr_up,  zcsr2csr_up  :: extract upper triangular part of a full csr matrix
!dcsr2csr_low, zcsr2csr_low :: extract lower triangular part of a full csr matrix

!dcsrmm, zcsrmm             :: mat-vec multiplication =>  b=alpha*A*x+beta*b
!zhcsrmm                    :: mat-vec multiplication =>  b=alpha*A*x+beta*b (A is complex Hermitian if UPLO=L,U)
!dzcsrmm                    :: same as zcsrmm but A is real;  alpha, beta, b are complex

!dcsr_transpose, zcsr_transpose :: transpose a csr matrix
!zcsr_htranspose                :: transpose conjugate a csr matrix

!dcsr_uplo_to_csr, zcsr_uplo_to_csr  :: convert a upper or lower csr matrix into a full csr one
!zhcsr_uplo_to_csr                   :: convert a upper or lower (Hermitian) csr matrix into a full csr one

!dcsr_convert_upper                 :: convert a upper,lower or full csr matrix into an upper csr matrix
!zcsr_convert_upper                 :: convert a upper,lower or full csr matrix into an upper csr matrix (complex symmetric)
!zhcsr_convert_full                 :: convert a upper,lower or full csr Hermitian matrix into a full csr matrix


!dcoo2csr :: transform coordinate format into csr format
!zcoo2csr :: transform coordinate format into csr format


!dsubmat                            :: extract csr submatrix from csr matrix
!zsubmat                            :: extract csr submatrix from csr matrix



subroutine daddcsr(N,M,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (NxM any M) sparse CSR matrices A and B and return a new (NxM) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  A,B,C are real double precision
  !
  !  Arguments
  !  =========
  !
  !  N,M      (input) INTEGER
  !          The number of rows/columns of the matrices A,B,C.  N,M >= 0.
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs addition, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine computes only values of the array isc of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value ic(n+1) - 1 is the actual 
  !         number of the elements in the arrays sc and jsc.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array jsc 
  !
  !         if opt=3, the routine has been called previously with the parameter opt=2. On exit it will
  !         return the array sc 
  !         
  !         
  !  alpha  (input) DOUBLE PRECISION
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa DOUBLE PRECISION)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  sb,isb,jsb (input) CSR format for the matrix B (sb DOUBLE PRECISION)
  !
  !  sc,isc,jsc (input/output) CSR format for the matrix C (sc DOUBLE PRECISION).
  !             see definition of 'opt' to determine input or output.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  double precision,dimension(*) :: sa,sb,sc
  integer,dimension(*) :: isa,isb,isc,jsa,jsb,jsc
  double precision :: alpha,beta
  integer :: opt,N,M


  integer :: i,ka,kb,kc,ja,jb
  logical :: next,row,col,val

  row=.false.
  col=.false.
  val=.false.

  if (opt==0) then
     row=.true.
     col=.true.
     val=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     col=.true.
  elseif (opt==3) then
     val=.true.
  endif



  kc=1
  if (row) isc(1)=kc
  do i=1,N
     ka=isa(i)
     kb=isb(i)

     next=.false.
     do while (.not.next)

        if (ka<=isa(i+1)-1) then
           ja=jsa(ka)
        else
           ja=m+1
        endif

        if (kb<=isb(i+1)-1) then
           jb=jsb(kb)
        else
           jb=m+1
        endif

        if ((ja==jb).and.(ja/=M+1)) then
           if (val) sc(kc)=alpha*sa(ka)+beta*sb(kb)
           if (col) jsc(kc)=ja
           ka=ka+1
           kb=kb+1
           kc=kc+1
        elseif (ja<jb) then 
           if (val) sc(kc)=alpha*sa(ka)
           if (col) jsc(kc)=ja
           ka=ka+1
           kc=kc+1
        elseif (jb<ja) then
           if (val) sc(kc)=beta*sb(kb)
           if (col) jsc(kc)=jb
           kb=kb+1
           kc=kc+1
        end if

        if ((ka>isa(i+1)-1).and.(kb>isb(i+1)-1)) next=.true.

     end do


     if (row) isc(i+1)=kc

  end do

end subroutine daddcsr










subroutine zaddcsr(N,M,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (NxM any M) sparse CSR matrices A and B and return a new (NxM) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  A,B,C are complex double precision
  !
  !  Arguments
  !  =========
  !
  !  N,M      (input) INTEGER
  !         The number of rows/columns of the matrices A,B,C.  N,M >= 0.
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs addition, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine computes only values of the array isc of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value ic(n+1) - 1 is the actual 
  !         number of the elements in the arrays sc and jsc.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array jsc 
  !
  !         if opt=3, the routine has been called previously with the parameter opt=2. On exit it will
  !         return the array sc 
  !         
  !         
  !  alpha  (input) COMPLEX DOUBLE PRECISION
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa COMPLEX DOUBLE PRECISION)
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  sb,isb,jsb (input) CSR format for the matrix B (sb COMPLEX DOUBLE PRECISION)
  !
  !  sc,isc,jsc (input/output) CSR format for the matrix C (sc COMPLEX DOUBLE PRECISION).
  !             see definition of 'opt' to determine input or output.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  complex(kind=kind(1.0d0)),dimension(*) :: sa,sb,sc
  integer,dimension(*) :: isa,isb,isc,jsa,jsb,jsc
  complex(kind=kind(1.0d0)) :: alpha,beta
  integer :: opt,N,M


  integer :: i,ka,kb,kc,ja,jb
  logical :: next,row,col,val

  row=.false.
  col=.false.
  val=.false.

  if (opt==0) then
     row=.true.
     col=.true.
     val=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     col=.true.
  elseif (opt==3) then
     val=.true.
  endif



  kc=1
  if (row) isc(1)=kc

  
  do i=1,N
     ka=isa(i)
     kb=isb(i)
     next=.false.
     do while (.not.next)

        if (ka<=isa(i+1)-1) then
           ja=jsa(ka)
        else
           ja=m+1
        endif

        if (kb<=isb(i+1)-1) then
           jb=jsb(kb)
        else
           jb=m+1
        endif
        if ((ja==jb).and.(ja/=M+1)) then
           if (val) sc(kc)=alpha*sa(ka)+beta*sb(kb)
           if (col) jsc(kc)=ja
           ka=ka+1
           kb=kb+1
           kc=kc+1
        elseif (ja<jb) then
           if (val) sc(kc)=alpha*sa(ka)
           if (col) jsc(kc)=ja
           ka=ka+1
           kc=kc+1
        elseif (jb<ja) then
           if (val) sc(kc)=beta*sb(kb)
           if (col) jsc(kc)=jb
           kb=kb+1
           kc=kc+1
        end if

        if ((ka>isa(i+1)-1).and.(kb>isb(i+1)-1)) next=.true.

     end do


     if (row) isc(i+1)=kc

  end do

end subroutine zaddcsr




subroutine zdaddcsr(N,M,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (NxM any M) sparse CSR matrices A and B and return a new (NxM) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  C is complex double precision, A and B are real double precision, (alpha,beta are COMPLEX)
  !  
  !  Arguments
  !  =========
  !
  !  N,M    (input) INTEGER
  !         The number of rows/columns of the matrices A,B,C.  N,M >= 0.
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs addition, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine computes only values of the array isc of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value ic(n+1) - 1 is the actual 
  !         number of the elements in the arrays sc and jsc.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array jsc 
  !
  !         if opt=3, the routine has been called previously with the parameter opt=2. On exit it will
  !         return the array sc 
  !         
  !         
  !  alpha  (input) COMPLEX DOUBLE PRECISION
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa DOUBLE PRECISION)
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  sb,isb,jsb (input) CSR format for the matrix B (sb DOUBLE PRECISION)
  !
  !  sc,isc,jsc (input/output) CSR format for the matrix C (sc COMPLEX DOUBLE PRECISION).
  !             see definition of 'opt' to determine input or output.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  double precision,dimension(*) :: sa,sb
  complex(kind=kind(1.0d0)),dimension(*) :: sc
  integer,dimension(*) :: isa,isb,isc,jsa,jsb,jsc
  complex(kind=kind(1.0d0)) :: alpha,beta
  integer :: opt,N,M
  integer :: i,ka,kb,kc,ja,jb
  logical :: next,row,col,val

  row=.false.
  col=.false.
  val=.false.

  if (opt==0) then
     row=.true.
     col=.true.
     val=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     col=.true.
  elseif (opt==3) then
     val=.true.
  endif



  kc=1
  if (row) isc(1)=kc

  do i=1,N
     ka=isa(i)
     kb=isb(i)
     next=.false.
     do while (.not.next)

        if (ka<=isa(i+1)-1) then
           ja=jsa(ka)
        else
           ja=M+1
        endif

        if (kb<=isb(i+1)-1) then
           jb=jsb(kb)
        else
           jb=M+1
        endif
        if ((ja==jb).and.(ja/=M+1)) then !!<<
           if (val) sc(kc)=alpha*sa(ka)+beta*sb(kb)
           if (col) jsc(kc)=ja
           ka=ka+1
           kb=kb+1
           kc=kc+1
        elseif (ja<jb) then
           if (val) sc(kc)=alpha*sa(ka)
           if (col) jsc(kc)=ja
           ka=ka+1
           kc=kc+1
        elseif (jb<ja) then
           if (val) sc(kc)=beta*sb(kb)
           if (col) jsc(kc)=jb
           kb=kb+1
           kc=kc+1
        end if

        if ((ka>isa(i+1)-1).and.(kb>isb(i+1)-1)) next=.true.

     end do


     if (row) isc(i+1)=kc

  end do

end subroutine zdaddcsr







subroutine zinc_addcsr(N,M,alpha,sa,isa,jsa,sb,isb,jsb)
  !  Purpose
  !  =======
  !  
  !  Increment matrix B by adding matrix A times alpha
  !  B=B+alpha*A 
  !
  !  A,B are complex double precision
  !  We consider (NxM) CSR matrices and that the csr format of B "includes" the csr format of A (matrices should be compatible, i.e. B contains *at least* all the same non-zeros elements than A)
  ! 
  !  Arguments
  !  =========
  !
  !  N,M      (input) INTEGER
  !         The number of rows/columns of the matrices A,B,C.  N,M >= 0.
  !         
  !  alpha  (input) COMPLEX DOUBLE PRECISION
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa COMPLEX DOUBLE PRECISION)
  !
  !
  !  sb,isb,jsb (input/output) CSR format for the matrix B (sb COMPLEX DOUBLE PRECISION)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  complex(kind=kind(1.0d0)),dimension(*) :: sa,sb
  integer,dimension(*) :: isa,isb,jsa,jsb
  complex(kind=kind(1.0d0)) :: alpha
  integer :: N,M

  integer :: i,ka,kb


  do i=1,N
     do ka=isa(i),isa(i+1)-1 ! element in row A
        do kb=isb(i),isb(i+1)-1 ! find where it belongs in row B 
           if (jsa(ka)==jsb(kb)) then
              sb(kb)=sb(kb)+alpha*sa(ka)
              cycle
           end if
        end do
    enddo
 end do

end subroutine zinc_addcsr





subroutine dcsr2csr_up(opt,N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine extract the upper triangular part (including diagonal) of a real CSR matrix A 
  !      into a new CSR matrix sA
  !      the size of sa, and sja for csr upper format may be overestimated 
  !
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs extraction, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine returns only values of the array sia of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value sia(n+1) - 1 is the actual 
  !         number of the elements in the arrays sa and sja.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array sa and sja
  !
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,j,k
  logical :: row,colval
!!!!
  row=.false.
  colval=.false.
  if (opt==0) then
     row=.true.
     colval=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     colval=.true.
  end if



  k=0
  if (row) sia(1)=1
  do i=1,n 
     if (row) sia(i+1)=sia(i)
     do j=ia(i),ia(i+1)-1
        if (ja(j)>=i) then
           k=k+1
           if (row) sia(i+1)=sia(i+1)+1
           if (colval) then
              sja(k)=ja(j)
              sa(k)=a(j)
           end if
        endif
     end do
  end do
end subroutine dcsr2csr_up


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine zcsr2csr_up(opt,N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine extract the upper triangular part (including diagonal) of a real CSR matrix A 
  !      into a new CSR matrix sA
  !      the size of sa, and sja for csr upper format may be overestimated 
  !
  !
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs extraction, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine returns only values of the array sia of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value sia(n+1) - 1 is the actual 
  !         number of the elements in the arrays sa and sja.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array sa and sja
  !
  !
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,j,k
  logical :: row,colval
!!!!
  row=.false.
  colval=.false.
  if (opt==0) then
     row=.true.
     colval=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     colval=.true.
  end if


  k=0
  if (row) sia(1)=1
  do i=1,n 
     if (row) sia(i+1)=sia(i)
     do j=ia(i),ia(i+1)-1
        if (ja(j)>=i) then
           k=k+1
           if (row) sia(i+1)=sia(i+1)+1
           if (colval) then
              sja(k)=ja(j)
              sa(k)=a(j)
           end if
        endif
     end do
  end do
end subroutine zcsr2csr_up




subroutine dcsr2csr_low(opt,N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine extract the lower triangular part (including diagonal) of a real CSR matrix A 
  !      into a new CSR matrix sA
  !      the size of sa, and sja for csr lower format may be overestimated 
  !
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs extraction, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine returns only values of the array sia of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value sia(n+1) - 1 is the actual 
  !         number of the elements in the arrays sa and sja.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array sa and sja
  !
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,j,k
  logical :: row,colval
!!!!
  row=.false.
  colval=.false.
  if (opt==0) then
     row=.true.
     colval=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     colval=.true.
  end if



  k=0
  if (row) sia(1)=1
  do i=1,n 
     if (row) sia(i+1)=sia(i)
     do j=ia(i),ia(i+1)-1
        if (ja(j)<=i) then
           k=k+1
           if (row) sia(i+1)=sia(i+1)+1
           if (colval) then
              sja(k)=ja(j)
              sa(k)=a(j)
           end if
        endif
     end do
  end do
end subroutine dcsr2csr_low
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine zcsr2csr_low(opt,N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine extract the lower triangular part (including diagonal) of a COMPLEX CSR matrix A 
  !      into a new CSR matrix sA
  !      the size of sa, and sja for csr lower format may be overestimated 
  !
  !
  !  opt    (input) INTEGER
  !         If opt=0, the routine performs extraction, the memory for the output arrays ic, jc, c  
  !         must be allocated beforehand.
  !
  !         if opt=1 the routine returns only values of the array sia of length n + 1, the memory 
  !         for this array must be allocated beforehand. On exit the value sia(n+1) - 1 is the actual 
  !         number of the elements in the arrays sa and sja.
  !  
  !         if opt=2, the routine has been called previously with the parameter opt=1. On exit it will
  !         return the array sa and sja
  !
  !
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,j,k
  logical :: row,colval
!!!!
  row=.false.
  colval=.false.
  if (opt==0) then
     row=.true.
     colval=.true.
  elseif (opt==1) then
     row=.true.
  elseif (opt==2) then
     colval=.true.
  end if


  k=0
  if (row) sia(1)=1
  do i=1,n 
     if (row) sia(i+1)=sia(i)
     do j=ia(i),ia(i+1)-1
        if (ja(j)<=i) then
           k=k+1
           if (row) sia(i+1)=sia(i+1)+1
           if (colval) then
              sja(k)=ja(j)
              sa(k)=a(j)
           end if
        endif
     end do
  end do
end subroutine zcsr2csr_low


SUBROUTINE dcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !
  !      REAL DOUBLE PRECISION version  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) DOUBLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T'
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T'. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  double precision :: alpha,beta
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*):: x
  double precision,dimension(*) ::b
!!!!!!!!        

  integer ::i,j,k,s,n1
  double precision,parameter :: DZERO=0.0d0

  if (TRANS=='N') then
     n1=N
  elseif (TRANS=='T') then
     n1=M
  endif

!!!!!!!! Initialization

if (beta/=DZERO) then
call DSCAL(n1*rhs,beta,b,1)
else
b(1:n1*rhs)=DZERO
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
              enddo
           end do
        end do

     elseif (TRANS=='T') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*a(k)*x(i+(s-1)*N)
              end do
           end do
        end do

     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif ((UPLO=='L').or.(UPLO=='U')) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
              b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
              if (j/=i) b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*a(k)*x(i+(s-1)*N)
           end do
           end do
        end do

  end if

end SUBROUTINE dcsrmm





SUBROUTINE zcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !      
  !      COMPLEX DOUBLE PRECISION version
  !           
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) DOUBLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T','C'.
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T','C'. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex (kind=kind(1.0d0)) :: alpha,beta
  complex (kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex (kind=kind(1.0d0)),dimension(*):: x
  complex (kind=kind(1.0d0)),dimension(*) ::b
!!!!!!!!

  integer ::i,j,k,s,n1
  complex (kind=kind(1.0d0)),parameter :: ZEROC=(0.0d0,0.0d0)
 
 if (TRANS=='N') then
     n1=N
  else
     n1=M
  endif

!!!!!!!! Initialization

if (beta/=ZEROC) then
call ZSCAL(n1*rhs,beta,b,1)
else
b(1:n1*rhs)=ZEROC
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then      
!call mkl_set_num_threads(1)
!        !$OMP PARALLEL DO SHARED(a,ia,ja,alpha,x,b,n,rhs) PRIVATE(i,k,j)
        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
              end do
              ! call ZAXPY(rhs,alpha*a(k),x(j,1),N,b(i,1),N) !!!<<< here
              !b(i,r)=b(i,r)+alpha*a(k)*x(j,r)
              !b(i,r)=b(i,r)+a(k)*x(j,r)
           end do
        end do
!        !$OMP END PARALLEL DO
!call mkl_set_num_threads(4)
     elseif (TRANS=='T') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*a(k)*x(i+(s-1)*N)
              end do
           end do
        end do

     elseif (TRANS=='C') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
                 enddo
           end do
        end do


     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif ((UPLO=='L').or.(UPLO=='U')) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((TRANS=='N').or.(TRANS=='T')) then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
              b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
              if (j/=i) b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*a(k)*x(i+(s-1)*N)
              enddo
           end do
        end do


 elseif (TRANS=='C') then

         do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
              b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
              if (j/=i) b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*conjg(a(k))*x(j+(s-1)*M)
              enddo
           end do
        end do

     end if

  end if
  
  
END SUBROUTINE zcsrmm





SUBROUTINE zhcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !
  !      Note: A is complex Hermitian if UPLO='L or UPLO='U'     
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric/hermitian CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric/hermitian CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) COMPLEX DOUBLE PRECISION
  ! 
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  X     (input) COMPLEX DOUBLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T','C'.
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  B     (input/output) COMPLEX DOUBLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T','C'. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex(kind=kind(1.0d0)) :: alpha,beta
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) :: x
  complex(kind=kind(1.0d0)),dimension(*) :: b
!!!!!!!!        

  integer ::i,j,k,s,n1
  double precision,parameter :: DZERO=0.0d0
  complex(kind=kind(1.0d0)),parameter :: ZEROC=(DZERO,DZERO)
 
 if (TRANS=='N') then
     n1=N
  else
     n1=M
  endif


!!!!!!!! Initialization

if (beta/=ZEROC) then
call ZSCAL(n1*rhs,beta,b,1)
else
b(1:n1*rhs)=ZEROC
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
                 enddo
           end do
        end do

     elseif (TRANS=='T') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*a(k)*x(i+(s-1)*N)
                 enddo
           end do
        end do

     elseif (TRANS=='C') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
                 enddo
           end do
        end do

     end if

  
  elseif ((UPLO=='L').or.(UPLO=='U')) then



if ((TRANS=='N').or.(TRANS=='C')) then


     do i=1,N
        do k=ia(i),ia(i+1)-1
           j=ja(k)
           do s=1,rhs
           b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
           if (j/=i) b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
        enddo
        end do
     end do

elseif (TRANS=='T') then

         do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
              b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
              if (j/=i) b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
              enddo
           end do
        end do

     end if



  end if

end SUBROUTINE zhcsrmm



SUBROUTINE dzcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !      
  !      DOUBLE PRECISION version
  !           
  !      Note: same as zcsrmm but A is real
  !
  !      Arguments: 
  !
  !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided- symmetric CSR (N must be equal to M)
  !       if UPLO='U' only the upper part of the matrix A is provided- symmetric CSR (N must be equal to M) 
  !
  !  TRANS (input)  CHARACTER*1
  !       if TRANS='N' o(A)=A   -- normal mode
  !       if TRANS='T' o(A)=A^T -- transpose mode 
  !       if TRANS='C' o(A)=A^H -- transpose conjugate mode
  ! 
  !  N    (input) INTEGER
  !        The number of row of the matrix A and row of matrix B.  N >= 0.
  !  M    (input) INTEGER 
  !        The number of column of matrix A and row of matrix X. M>=0; M=N (square matrix) for UPLO=L,U
  !  rhs  (input) INTEGER
  !        The number of column of matrices B and C
  ! 
  !  alpha (input) COMPLEX DOUBLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  X     (input) COMPLEX DOUBLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T','C'.
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  B     (input/output) COMPLEX DOUBLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T','C'. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex (kind=kind(1.0d0)) :: alpha,beta
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex (kind=kind(1.0d0)),dimension(*):: x
  complex (kind=kind(1.0d0)),dimension(*) ::b
!!!!!!!!

  integer ::i,j,k,s,n1
  complex (kind=kind(1.0d0)),parameter :: ZEROC=(0.0d0,0.0d0)
 
 
 if (TRANS=='N') then
     n1=N
  else
     n1=M
  endif


!!!!!!!! Initialization

if (beta/=ZEROC) then
call ZSCAL(n1*rhs,beta,b,1)
else
b(1:n1*rhs)=ZEROC
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
                 enddo
           end do
        end do

     elseif ((TRANS=='T').or.(TRANS=='C')) then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*a(k)*x(i+(s-1)*N)
                 enddo
           end do
        end do

     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif ((UPLO=='L').or.(UPLO=='U')) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
              b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
              if (j/=i) b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*a(k)*x(i+(s-1)*N)
           enddo
           end do
        end do

  end if
      
  
END SUBROUTINE dzcsrmm







subroutine dcsr_transpose(N,M,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine transpose a double real CSR matrix A  
  !      into a new CSR matrix sA
  !
  !  N,M    (input) INTEGER
  !          The number of rows/columns of the matrices A.  N,M >= 0.
  !          The number of columns/rows of the matrices sA. N,M >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N,M
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
  integer :: i,j,k,istart
  integer,dimension(:),allocatable :: tab
!!!!

  allocate(tab(M))
  !! #of elements by row
  tab=0
  istart=N
  do i=1,N
     do j=ia(i),ia(i+1)-1
        tab(ja(j))=tab(ja(j))+1
        if (ja(j)<istart) istart=ja(j)
     end do
  end do

  !! find sia
  do i=1,istart-1
     sia(i)=1 !!<<0 before
  end do
  sia(istart)=1
  do i=istart+1,M+1
     sia(i)=sia(i-1)+tab(i-1)
  end do

  tab=sia(1:M)
  !! find jsa,sa
  do i=1,N
     do k=ia(i),ia(i+1)-1
        sja(tab(ja(k)))=i
        sa(tab(ja(k)))=a(k)
        tab(ja(k))=tab(ja(k))+1
     end do
  end do

  deallocate(tab)

end subroutine dcsr_transpose



subroutine zcsr_transpose(N,M,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the transpose  of 
  !      a double complex CSR matrix A into a new CSR matrix sA
  !
  !
  !  N,M    (input) INTEGER
  !          The number of rows/columns of the matrices A.  N,M >= 0.
  !          The number of columns/rows of the matrices sA. N,M >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N,M
   complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
   complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
  integer :: i,j,k,istart
   integer,dimension(:),allocatable :: tab
!!!!

  allocate(tab(M))
!!!!

  !! #of elements by row
  tab=0
  istart=N
  do i=1,N
     do j=ia(i),ia(i+1)-1
        tab(ja(j))=tab(ja(j))+1
        if (ja(j)<istart) istart=ja(j)
     end do
  end do

  !! find sia
  do i=1,istart-1
     sia(i)=1 !!<<0
  end do
  sia(istart)=1
  do i=istart+1,M+1
     sia(i)=sia(i-1)+tab(i-1)
  end do

  tab=sia(1:M)
  !! find jsa,sa
  do i=1,N
     do k=ia(i),ia(i+1)-1
        sja(tab(ja(k)))=i
        sa(tab(ja(k)))=a(k)
        tab(ja(k))=tab(ja(k))+1
     end do
  end do

 deallocate(tab) 

end subroutine zcsr_transpose




subroutine zcsr_htranspose(N,M,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the transpose conjugate of 
  !      a double complex CSR matrix A into a new CSR matrix sA
  !
  !
  !  N,M    (input) INTEGER
  !          The number of rows/columns of the matrices A.  N,M >= 0.
  !          The number of columns/rows of the matrices sA. N,M >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N,M
   complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
   complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
  integer :: i,j,k,istart
   integer,dimension(:),allocatable :: tab
!!!!

  allocate(tab(M))
!!!!

  !! #of elements by row
  tab=0
  istart=N
  do i=1,N
     do j=ia(i),ia(i+1)-1
        tab(ja(j))=tab(ja(j))+1
        if (ja(j)<istart) istart=ja(j)
     end do
  end do

  !! find sia
  do i=1,istart-1
     sia(i)=1 !!<<0
  end do
  sia(istart)=1
  do i=istart+1,M+1
     sia(i)=sia(i-1)+tab(i-1)
  end do

  tab=sia(1:M)
  !! find jsa,sa
  do i=1,N
     do k=ia(i),ia(i+1)-1
        sja(tab(ja(k)))=i
        sa(tab(ja(k)))=conjg(a(k))
        tab(ja(k))=tab(ja(k))+1
     end do
  end do

deallocate(tab)
  
end subroutine zcsr_htranspose



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine dcsr_uplo_to_csr(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a double real upper or lower CSR matrix A  
  !      into a new full CSR matrix sA
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  !include 'f90_noruntime_interface.fi'
  integer :: N
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  double precision,dimension(:),allocatable ::sb
  integer,dimension(:),allocatable :: sib
  integer,dimension(:),allocatable :: sjb
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
!!!!
  nnza=ia(n+1)-1

  allocate(sb(nnza))
  allocate(sib(N+1))
  allocate(sjb(nnza))


!!! transpose
  call dcsr_transpose(N,N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=DZERO
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call daddcsr(N,N,opt,DONE,a,ia,ja,DONE,sb,sib,sjb,sa,sia,sja)


  deallocate(sb)
  deallocate(sib)
  deallocate(sjb)



end subroutine dcsr_uplo_to_csr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zcsr_uplo_to_csr(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a complex symetric upper or lower CSR matrix A  
  !      into a new full CSR matrix sA 
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  !include 'f90_noruntime_interface.fi'
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  complex(kind=kind(1.0d0)),dimension(:),allocatable ::sb
  integer,dimension(:),allocatable :: sib
  integer,dimension(:),allocatable :: sjb
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
!!!!
  nnza=ia(n+1)-1

  allocate(sb(nnza))
  allocate(sib(N+1))
  allocate(sjb(nnza))


!!! transpose
  call zcsr_transpose(N,N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=ZEROC
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call zaddcsr(N,N,opt,ONEC,a,ia,ja,ONEC,sb,sib,sjb,sa,sia,sja)


  deallocate(sb)
  deallocate(sib)
  deallocate(sjb)



end subroutine zcsr_uplo_to_csr




subroutine zhcsr_uplo_to_csr(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a complex Hermitian upper or lower CSR matrix A  
  !      into a new full CSR matrix sA 
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  !include 'f90_noruntime_interface.fi'
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  complex(kind=kind(1.0d0)),dimension(:),allocatable ::sb
  integer,dimension(:),allocatable :: sib
  integer,dimension(:),allocatable :: sjb
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
!!!!
  nnza=ia(n+1)-1

  allocate(sb(nnza)) 
  allocate(sib(N+1)) 
  allocate(sjb(nnza))

!!! transpose
  call zcsr_htranspose(N,N,a,ia,ja,sb,sib,sjb)
 
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=ZEROC
     end do
  enddo

  
  !! add the two matrices to form the full matrix
  opt=0
  call zaddcsr(N,N,opt,ONEC,a,ia,ja,ONEC,sb,sib,sjb,sa,sia,sja)

  deallocate(sb)
  deallocate(sib)
  deallocate(sjb)



end subroutine zhcsr_uplo_to_csr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dcsr_convert_upper(n,UPLO,sa,isa,jsa,ssa,sisa,sjsa)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Purpose: This subroutine convert a double real upper, loweror full CSR matrix A  
  !      into a new CSR upper matrix sA
    !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
 !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided
  !       if UPLO='U' only the upper part of the matrix A is provided
  !
    !
  !  sa,isa,jsa (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  ssa,sisa,sjsa (output) CSR format for the matrix sA (sa  DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: n
  character(len=1) :: UPLO
  double precision,dimension(*) :: sa,ssa
  integer,dimension(*) :: isa,jsa,sisa,sjsa
  integer :: opt,nnz

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
      
       opt=0
       call dcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
      
!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 
       nnz=isa(n+1)-1
       call DCOPY(nnz,sa,1,ssa,1)
              sisa(1:n+1)= isa(1:n+1)
              sjsa(1:nnz)= jsa(1:nnz)

    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr'
       call dcsr_transpose(N,N,sa,isa,jsa,ssa,sisa,sjsa)
    endif
  

end subroutine dcsr_convert_upper



  subroutine zcsr_convert_upper(n,UPLO,sa,isa,jsa,ssa,sisa,sjsa)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      Purpose: This subroutine convert a double complex upper, loweror full CSR matrix A  
  !      into a new CSR upper matrix sA (symmetric matrix)
    !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
 !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided
  !       if UPLO='U' only the upper part of the matrix A is provided
  !
    !
  !  sa,isa,jsa (input) CSR format for the matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  ssa,sisa,sjsa (output) CSR format for the matrix sA (sa  COMPLEX DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: n
  character(len=1) :: UPLO
  complex(kind=kind(1.0d0)),dimension(*) :: sa,ssa
  integer,dimension(*) :: isa,jsa,sisa,sjsa
  integer :: opt,nnz

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr
      
       opt=0
       call zcsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)
      
!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 
       nnz=isa(n+1)-1
       call ZCOPY(nnz,sa,1,ssa,1)
              sisa(1:n+1)= isa(1:n+1)
              sjsa(1:nnz)= jsa(1:nnz)

    elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr'
       call zcsr_transpose(N,N,sa,isa,jsa,ssa,sisa,sjsa)
    endif
  

end subroutine zcsr_convert_upper




  subroutine zhcsr_convert_full(n,UPLO,sa,isa,jsa,ssa,sisa,sjsa)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a double complex upper, lower or full Hermitian CSR matrix A  
  !      into a new full CSR  matrix sA
    !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
 !  UPLO (input) CHARACTER*1
  !       if UPLO='F' all the matrix A is provided in CSR format ('F'ull format)
  !       if UPLO='L' only the lower part of the matrix A is provided
  !       if UPLO='U' only the upper part of the matrix A is provided
  !
    !
  !  sa,isa,jsa (input) CSR format for the Hermitian matrix A (a COMPLEX DOUBLE PRECISION)
  !
  !  ssa,sisa,sjsa (output) CSR format for the Hermitian matrix sA (sa COMPLEX DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: n
  character(len=1) :: UPLO
  complex(kind=kind(1.0d0)),dimension(*) :: sa,ssa
  integer,dimension(*) :: isa,jsa,sisa,sjsa
  integer :: opt,nnz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((UPLO=='F').or.(UPLO=='f')) then !!! direct copy

    nnz=isa(n+1)-1
       call ZCOPY(nnz,sa,1,ssa,1)
              sisa(1:n+1)= isa(1:n+1)
              sjsa(1:nnz)= jsa(1:nnz)
      
!!!!!!!!!!!!!!!!!!!!!!!!!
           else!! upper-csr or lower-csr to full csr

   call zhcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
end if
  

end subroutine zhcsr_convert_full








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









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dsubmat (n,i1,i2,j1,j2,a,ja,ia,ao,jao,iao)
    !*******************************************************************************
    ! SUBMAT extracts the submatrix A(i1:i2,j1:j2) 
    !*******************************************************************************
    !a,ja,ia = original matrix in CSR format
    !ao,jao,iao = extracted matrix in CSR format
    implicit none
    integer :: i1,i2,j1,j2,n
    double precision,dimension(*):: a
    integer,dimension(*) :: ia,ja
    double precision,dimension(*) :: ao
    integer,dimension(*) :: iao,jao

    integer :: klen,i,ii,k1,k2,k,j,nr,nc

    nr = i2-i1+1
    nc = j2-j1+1
    klen = 0
    !!
    do  i = 1,nr
       ii = i1+i-1
       k1 = ia(ii)
       k2 = ia(ii+1)-1
       iao(i) = klen+1
       do  k=k1,k2
          j = ja(k)
          if ((j .ge. j1).and. (j .le. j2)) then
             klen = klen+1
             ao(klen) = a(k)
             jao(klen) = j - j1+1
          endif
end do
    end do
    iao(nr+1) = klen+1
  end subroutine dsubmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zsubmat (n,i1,i2,j1,j2,a,ja,ia,ao,jao,iao)
    !*******************************************************************************
    ! SUBMAT extracts the submatrix A(i1:i2,j1:j2) 
    !*******************************************************************************
    !a,ja,ia = original matrix in CSR format
    !ao,jao,iao = extracted matrix in CSR format
    implicit none
    integer :: i1,i2,j1,j2,n
     complex(kind=kind(1.0d0)),dimension(*):: a
    integer,dimension(*) :: ia,ja
     complex(kind=kind(1.0d0)),dimension(*) :: ao
    integer,dimension(*) :: iao,jao

    integer :: klen,i,ii,k1,k2,k,j,nr,nc

    nr = i2-i1+1
    nc = j2-j1+1
    klen = 0
    !!
    do  i = 1,nr
       ii = i1+i-1
       k1 = ia(ii)
       k2 = ia(ii+1)-1
       iao(i) = klen+1
       do  k=k1,k2
          j = ja(k)
          if ((j .ge. j1).and. (j .le. j2)) then
             klen = klen+1
             ao(klen) = a(k)
             jao(klen) = j - j1+1
          endif
end do
    end do
    iao(nr+1) = klen+1
  end subroutine zsubmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









