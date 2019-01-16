  !=========================================================================================
  !Copyright (c) 2009-2015, The Regents of the University of Massachusetts, Amherst.
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

!dcsr2csr_up,  zcsr2csr_up  :: extract upper triangular part of a full csr matrix
!dcsr2csr_low, zcsr2csr_low :: extract lower triangular part of a full csr matrix

!dcsrmm, zcsrmm             :: mat-vec multiplication =>  b=alpha*A*x+beta*b
!zhcsrmm                    :: mat-vec multiplication =>  b=alpha*A*x+beta*b (A is complex Hermitian if UPLO=L,U)
!dzcsrmm                    :: same as zcsrmm but A is real;  alpha, beta, b are complex

!dcsr_transpose, zcsr_transpose :: transpose a csr matrix
!zcsr_htranspose                :: transpose conjugate a csr matrix

!dcsr_uplo_to_csr, zcsr_uplo_to_csr  :: convert a upper or lower csr matrix into a full csr one
!zhcsr_uplo_to_csr                   :: convert a upper or lower (Hermitian) csr matrix into a full csr one



subroutine daddcsr(N,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (square) sparse CSR matrices A and B and return a new (square) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  A,B,C are real double precision
  !
  !  Arguments
  !  =========
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,B,C.  N >= 0.
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
  integer :: opt,N


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
           ja=n+1
        endif

        if (kb<=isb(i+1)-1) then
           jb=jsb(kb)
        else
           jb=n+1
        endif

        if (ja==jb) then
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










subroutine zaddcsr(N,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (square) sparse CSR matrices A and B and return a new (square) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  A,B,C are complex double precision
  !
  !  Arguments
  !  =========
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,B,C.  N >= 0.
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
  integer :: opt,N


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
           ja=n+1
        endif

        if (kb<=isb(i+1)-1) then
           jb=jsb(kb)
        else
           jb=n+1
        endif

        if (ja==jb) then
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




subroutine zdaddcsr(N,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (square) sparse CSR matrices A and B and return a new (square) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  C is complex double precision, A and B are real double precision, (alpha,beta are COMPLEX)
  !  
  !  Arguments
  !  =========
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,B,C.  N >= 0.
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
  integer :: opt,N
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
           ja=n+1
        endif

        if (kb<=isb(i+1)-1) then
           jb=jsb(kb)
        else
           jb=n+1
        endif

        if (ja==jb) then
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
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  double precision :: alpha,beta
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(M,*):: x
  double precision,dimension(N,*) ::b
!!!!!!!!        

  integer ::i,j,k
  double precision,parameter :: DZERO=0.0d0
  

!!!!!!!! Initialization

if (beta/=DZERO) then
call DSCAL(N*rhs,b,beta,1)
else
b(1:N,1:rhs)=DZERO
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
           end do
        end do

     elseif (TRANS=='T') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*a(k)*x(i,1:rhs)
           end do
        end do

     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif ((UPLO=='L').or.(UPLO=='U')) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
              if (j/=i) b(j,1:rhs)=b(j,1:rhs)+alpha*a(k)*x(i,1:rhs)               
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
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex (kind=kind(1.0d0)) :: alpha,beta
  complex (kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex (kind=kind(1.0d0)),dimension(M,*):: x
  complex (kind=kind(1.0d0)),dimension(N,*) ::b
!!!!!!!!

  integer ::i,j,k
  complex (kind=kind(1.0d0)),parameter :: ZEROC=(0.0d0,0.0d0)
 


!!!!!!!! Initialization

if (beta/=ZEROC) then
call ZSCAL(N*rhs,b,beta,1)
else
b(1:N,1:rhs)=ZEROC
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
           end do
        end do

     elseif (TRANS=='T') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*a(k)*x(i,1:rhs)
           end do
        end do

     elseif (TRANS=='C') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*conjg(a(k))*x(i,1:rhs)
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
              b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
              if (j/=i) b(j,1:rhs)=b(j,1:rhs)+alpha*a(k)*x(i,1:rhs)               
           end do
        end do


 elseif (TRANS=='C') then

         do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*conjg(a(k))*x(i,1:rhs)
              if (j/=i) b(i,1:rhs)=b(i,1:rhs)+alpha*conjg(a(k))*x(j,1:rhs)               
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
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) COMPLEX DOUBLE PRECISION
  !
  !  B     (input/output) COMPLEX DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex(kind=kind(1.0d0)) :: alpha,beta
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(M,*) :: x
  complex(kind=kind(1.0d0)),dimension(N,*) :: b
!!!!!!!!        

  integer ::i,j,k
  double precision,parameter :: DZERO=0.0d0
  complex(kind=kind(1.0d0)),parameter :: ZEROC=(DZERO,DZERO)
 

!!!!!!!! Initialization

if (beta/=ZEROC) then
call ZSCAL(N*rhs,b,beta,1)
else
b(1:N,1:rhs)=ZEROC
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
           end do
        end do

     elseif (TRANS=='T') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*a(k)*x(i,1:rhs)
           end do
        end do

     elseif (TRANS=='C') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*conjg(a(k))*x(i,1:rhs)
           end do
        end do

     end if

  
  elseif ((UPLO=='L').or.(UPLO=='U')) then



if ((TRANS=='N').or.(TRANS=='C')) then


     do i=1,N
        do k=ia(i),ia(i+1)-1
           j=ja(k)
           b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
   if (j/=i) b(j,1:rhs)=b(j,1:rhs)+alpha*conjg(a(k))*x(i,1:rhs)       
        end do
     end do

elseif (TRANS=='T') then

         do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*conjg(a(k))*x(i,1:rhs)
              if (j/=i) b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)               
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
  !      Note: A is real
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
  !        matrix of size (Mxrhs)
  !
  !  beta  (input) DOUBLE PRECISION
  !
  !  B     (input/output) DOUBLE PRECISION
  !        matrix of size (Nxrhs). On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex (kind=kind(1.0d0)) :: alpha,beta
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex (kind=kind(1.0d0)),dimension(M,*):: x
  complex (kind=kind(1.0d0)),dimension(N,*) ::b
!!!!!!!!

  integer ::i,j,k
  complex (kind=kind(1.0d0)),parameter :: ZEROC=(0.0d0,0.0d0)
 


!!!!!!!! Initialization

if (beta/=ZEROC) then
call ZSCAL(N*rhs,b,beta,1)
else
b(1:N,1:rhs)=ZEROC
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (UPLO=='F') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (TRANS=='N') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
           end do
        end do

     elseif ((TRANS=='T').or.(TRANS=='C')) then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(j,1:rhs)=b(j,1:rhs)+alpha*a(k)*x(i,1:rhs)
           end do
        end do

     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif ((UPLO=='L').or.(UPLO=='U')) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              b(i,1:rhs)=b(i,1:rhs)+alpha*a(k)*x(j,1:rhs)
              if (j/=i) b(j,1:rhs)=b(j,1:rhs)+alpha*a(k)*x(i,1:rhs)               
           end do
        end do

  end if
      
  
END SUBROUTINE dzcsrmm







subroutine dcsr_transpose(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine transpose a double real CSR matrix A  
  !      into a new CSR matrix sA
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
  integer :: i,j,k,istart
  integer,dimension(1:N) :: tab
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
     sia(i)=0
  end do
  sia(istart)=1
  do i=istart+1,N+1
     sia(i)=sia(i-1)+tab(i-1)
  end do

  tab=sia(1:N)
  !! find jsa,sa
  do i=1,N
     do k=ia(i),ia(i+1)-1
        sja(tab(ja(k)))=i
        sa(tab(ja(k)))=a(k)
        tab(ja(k))=tab(ja(k))+1
     end do
  end do


end subroutine dcsr_transpose


subroutine zcsr_transpose(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the transpose conjugate of 
  !      a double complex CSR matrix A into a new CSR matrix sA
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,j,k,istart
  integer,dimension(1:N) :: tab
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
     sia(i)=0
  end do
  sia(istart)=1
  do i=istart+1,N+1
     sia(i)=sia(i-1)+tab(i-1)
  end do

  tab=sia(1:N)
  !! find jsa,sa
  do i=1,N
     do k=ia(i),ia(i+1)-1
        sja(tab(ja(k)))=i
  !      if (ja(k)/=i) then
           sa(tab(ja(k)))=a(k)
  !      else
  !         sa(tab(ja(k)))=a(k)
  !      end if
        tab(ja(k))=tab(ja(k))+1
     end do
  end do


end subroutine zcsr_transpose



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine zcsr_htranspose(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the transpose conjugate of 
  !      a double complex CSR matrix A into a new CSR matrix sA
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a DOUBLE COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa DOUBLE COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,j,k,istart
  integer,dimension(1:N) :: tab
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
     sia(i)=0
  end do
  sia(istart)=1
  do i=istart+1,N+1
     sia(i)=sia(i-1)+tab(i-1)
  end do

  tab=sia(1:N)
  !! find jsa,sa
  do i=1,N
     do k=ia(i),ia(i+1)-1
        sja(tab(ja(k)))=i
  !      if (ja(k)/=i) then
           sa(tab(ja(k)))=conjg(a(k))
  !      else
  !         sa(tab(ja(k)))=a(k)
  !      end if
        tab(ja(k))=tab(ja(k))+1
     end do
  end do


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
  include 'f90_noruntime_interface.fi'
  integer :: N
  double precision,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  double precision,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  double precision,dimension(:),pointer ::sb
  integer,dimension(:),pointer :: sib
  integer,dimension(:),pointer :: sjb
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
!!!!
  nnza=ia(n+1)-1

  call wallocate_1d(sb, nnza, info)
  call wallocate_1i(sib, N+1, info)
  call wallocate_1i(sjb, nnza, info)


!!! transpose
  call dcsr_transpose(N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=DZERO
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call daddcsr(N,opt,DONE,a,ia,ja,DONE,sb,sib,sjb,sa,sia,sja)


  call wdeallocate_1d(sb)
  call wdeallocate_1i(sib)
  call wdeallocate_1i(sjb)



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
  include 'f90_noruntime_interface.fi'
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  complex(kind=kind(1.0d0)),dimension(:),pointer ::sb
  integer,dimension(:),pointer :: sib
  integer,dimension(:),pointer :: sjb
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
!!!!
  nnza=ia(n+1)-1

  call wallocate_1z(sb, nnza, info)
  call wallocate_1i(sib, N+1, info)
  call wallocate_1i(sjb, nnza, info)


!!! transpose
  call zcsr_transpose(N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=ZEROC
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call zaddcsr(N,opt,ONEC,a,ia,ja,ONEC,sb,sib,sjb,sa,sia,sja)


  call wdeallocate_1z(sb)
  call wdeallocate_1i(sib)
  call wdeallocate_1i(sjb)



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
  include 'f90_noruntime_interface.fi'
  integer :: N
  complex(kind=kind(1.0d0)),dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex(kind=kind(1.0d0)),dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  complex(kind=kind(1.0d0)),dimension(:),pointer ::sb
  integer,dimension(:),pointer :: sib
  integer,dimension(:),pointer :: sjb
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
!!!!
  nnza=ia(n+1)-1

  call wallocate_1z(sb, nnza, info)
  call wallocate_1i(sib, N+1, info)
  call wallocate_1i(sjb, nnza, info)


!!! transpose
  call zcsr_htranspose(N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=ZEROC
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call zaddcsr(N,opt,ONEC,a,ia,ja,ONEC,sb,sib,sjb,sa,sia,sja)


  call wdeallocate_1z(sb)
  call wdeallocate_1i(sib)
  call wdeallocate_1i(sjb)



end subroutine zhcsr_uplo_to_csr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










