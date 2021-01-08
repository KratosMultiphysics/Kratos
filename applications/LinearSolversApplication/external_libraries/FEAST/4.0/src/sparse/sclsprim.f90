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

!!!!!!!! SINGLE PRECISION VERSION

!!!!! List of routines 


!saddcsr, caddcsr           :: add two square csr matrices and return a new one (C=alpha*A+beta*B)
!csaddcsr                   :: same as zaddcsr but A, B are real;  alpha, beta, C are complex
!cinc_addcsr                :: Increment matrix B by adding matrix A times alpha (sparsity pattern of A included in B)


!scsr2csr_up,  ccsr2csr_up  :: extract upper triangular part of a full csr matrix
!scsr2csr_low, ccsr2csr_low :: extract lower triangular part of a full csr matrix

!scsrmm, ccsrmm             :: mat-vec multiplication =>  b=alpha*A*x+beta*b
!chcsrmm                    :: mat-vec multiplication =>  b=alpha*A*x+beta*b (A is complex Hermitian if UPLO=L,U)
!sccsrmm                    :: same as zcsrmm but A is real;  alpha, beta, b are complex

!scsr_transpose, ccsr_transpose :: transpose a csr matrix
!ccsr_htranspose                :: transpose conjugate a csr matrix

!scsr_uplo_to_csr, ccsr_uplo_to_csr  :: convert a upper or lower csr matrix into a full csr one
!chcsr_uplo_to_csr                   :: convert a upper or lower (Hermitian) csr matrix into a full csr one


!scsr_convert_upper                 :: convert a upper,lower or full csr matrix into an upper csr matrix
!chcsr_convert_full                 :: convert a upper,lower or full csr Hermitian matrix into a full csr matrix


!scoo2csr :: transform coordinate format into csr format
!ccoo2csr :: transform coordinate format into csr format



!ssubmat                            :: extract csr submatrix from csr matrix
!csubmat                            :: extract csr submatrix from csr matrix



subroutine saddcsr(N,M,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (NxM any M) sparse CSR matrices A and B and return a new (NxM) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  A,B,C are real
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
  !  alpha  (input) REAL
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa REAL)
  !
  !  beta  (input) REAL
  !
  !  sb,isb,jsb (input) CSR format for the matrix B (sb REAL)
  !
  !  sc,isc,jsc (input/output) CSR format for the matrix C (sc REAL).
  !             see definition of 'opt' to determine input or output.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  real,dimension(*) :: sa,sb,sc
  integer,dimension(*) :: isa,isb,isc,jsa,jsb,jsc
  real :: alpha,beta
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
        elseif (ja<jb)  then
           if (val) sc(kc)=alpha*sa(ka)
           if (col) jsc(kc)=ja
           ka=ka+1
           kc=kc+1
        elseif (jb<ja)  then
           if (val) sc(kc)=beta*sb(kb)
           if (col) jsc(kc)=jb
           kb=kb+1
           kc=kc+1
        end if

        if ((ka>isa(i+1)-1).and.(kb>isb(i+1)-1)) next=.true.

     end do


     if (row) isc(i+1)=kc

  end do

end subroutine saddcsr










subroutine caddcsr(N,M,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (NxM any M) sparse CSR matrices A and B and return a new (NxM) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  A,B,C are complex 
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
  !  alpha  (input) COMPLEX 
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa COMPLEX)
  !
  !  beta  (input) COMPLEX
  !
  !  sb,isb,jsb (input) CSR format for the matrix B (sb COMPLEX)
  !
  !  sc,isc,jsc (input/output) CSR format for the matrix C (sc COMPLEX).
  !             see definition of 'opt' to determine input or output.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  complex ,dimension(*) :: sa,sb,sc
  integer,dimension(*) :: isa,isb,isc,jsa,jsb,jsc
  complex  :: alpha,beta
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

end subroutine caddcsr




subroutine csaddcsr(N,M,opt,alpha,sa,isa,jsa,beta,sb,isb,jsb,sc,isc,jsc)
  !  Purpose
  !  =======
  !
  !  Addition of two (NxM any M) sparse CSR matrices A and B and return a new (NxM) CSR matrix C 
  !  C=alpha*A+beta*B 
  !
  !  A is complex, B and C are real, (alpha,beta are COMPLEX)
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
  !  alpha  (input) COMPLEX
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa REAL)
  !
  !  beta  (input) COMPLEX
  !
  !  sb,isb,jsb (input) CSR format for the matrix B (sb REAL)
  !
  !  sc,isc,jsc (input/output) CSR format for the matrix C (sc COMPLEX).
  !             see definition of 'opt' to determine input or output.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  real,dimension(*) :: sa,sb
  complex ,dimension(*) :: sc
  integer,dimension(*) :: isa,isb,isc,jsa,jsb,jsc
  complex  :: alpha,beta
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
        elseif (ja<jb)  then
           if (val) sc(kc)=alpha*sa(ka)
           if (col) jsc(kc)=ja
           ka=ka+1
           kc=kc+1
        elseif (jb<ja)  then
           if (val) sc(kc)=beta*sb(kb)
           if (col) jsc(kc)=jb
           kb=kb+1
           kc=kc+1
        end if

        if ((ka>isa(i+1)-1).and.(kb>isb(i+1)-1)) next=.true.

     end do


     if (row) isc(i+1)=kc

  end do

end subroutine csaddcsr










subroutine cinc_addcsr(N,M,alpha,sa,isa,jsa,sb,isb,jsb)
  !  Purpose
  !  =======
  !  
  !  Increment matrix B by adding matrix A times alpha
  !  B=B+alpha*A 
  !
  !  A,B are complex single precision
  !  We consider (NxM) CSR matrices and that the csr format of B "includes" the csr format of A (matrices should be compatible, i.e. B contains *at least* all the same non-zeros elements than A)
  ! 
  !  Arguments
  !  =========
  !
  !  N,M      (input) INTEGER
  !         The number of rows/columns of the matrices A,B,C.  N,M >= 0.
  !         
  !  alpha  (input) COMPLEX SINGLE PRECISION
  !
  !  sa,isa,jsa (input) CSR format for the matrix A (sa COMPLEX SINGLE PRECISION)
  !
  !
  !  sb,isb,jsb (input/output) CSR format for the matrix B (sb COMPLEX SINGLE PRECISION)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  complex,dimension(*) :: sa,sb
  integer,dimension(*) :: isa,isb,jsa,jsb
  complex :: alpha
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

end subroutine cinc_addcsr





subroutine scsr2csr_up(opt,N,a,ia,ja,sa,sia,sja)
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
  !  a,ia,ja (input) CSR format for the matrix A (a REAL)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa REAL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  real,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  real,dimension(*) ::sa
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
end subroutine scsr2csr_up


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine ccsr2csr_up(opt,N,a,ia,ja,sa,sia,sja)
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
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  complex ,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex ,dimension(*) ::sa
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
end subroutine ccsr2csr_up




subroutine scsr2csr_low(opt,N,a,ia,ja,sa,sia,sja)
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
  !  a,ia,ja (input) CSR format for the matrix A (a REAL)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa REAL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  real,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  real,dimension(*) ::sa
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
end subroutine scsr2csr_low
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine ccsr2csr_low(opt,N,a,ia,ja,sa,sia,sja)
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
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: opt,N
  complex ,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex ,dimension(*) ::sa
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
end subroutine ccsr2csr_low


SUBROUTINE scsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !
  !      REAL SINGLE PRECISION version             
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
  !  alpha (input) SINGLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a SINGLE PRECISION)
  !
  !  X     (input) DOUBLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T'
  !
  !  beta  (input) SINGLE PRECISION
  !
  !  B     (input/output) SINGLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T'. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  real :: alpha,beta
  real,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  real,dimension(*):: x
  real,dimension(*) ::b
!!!!!!!!        

  integer ::i,j,k,s,n1
  real,parameter :: SZERO=0.0e0



  if (TRANS=='N') then
     n1=N
  elseif (TRANS=='T') then
     n1=M
  endif

!!!!!!!! Initialization

  if (beta/=SZERO) then
     call SSCAL(n1*rhs,beta,b,1)
  else
     b(1:n1*rhs)=SZERO
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

end SUBROUTINE scsrmm





SUBROUTINE ccsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !      
  !      COMPLEX SINGLE PRECISION version
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
  !  alpha (input) SINGLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a SINGLE PRECISION)
  !
  !  X     (input) SINGLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T','C'.
  !
  !  beta  (input) SINGLE PRECISION
  !
  !  B     (input/output) SINGLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T','C'. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex  :: alpha,beta
  complex ,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex ,dimension(*):: x
  complex ,dimension(*) ::b
!!!!!!!!

  integer ::i,j,k,s,n1
  complex,parameter :: ZEROC=(0.0e0,0.0e0)



  if (TRANS=='N') then
     n1=N
  else
     n1=M
  endif

!!!!!!!! Initialization

  if (beta/=ZEROC) then
     call CSCAL(n1*rhs,beta,b,1)
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
              end do
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
              end do
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
              end do
           end do
        end do


     elseif (TRANS=='C') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
                 if (j/=i) b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*conjg(a(k))*x(j+(s-1)*M)
              end do
           end do
        end do

     end if

  end if


END SUBROUTINE ccsrmm





SUBROUTINE chcsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
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
  !  alpha (input) COMPLEX SINGLE PRECISION
  ! 
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX SINGLE PRECISION)
  !
  !  X     (input) COMPLEX SINGLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T.'
  !
  !  beta  (input) COMPLEX SINGLE PRECISION
  !
  !  B     (input/output) COMPLEX SINGLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T','C'. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex :: alpha,beta
  complex,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex,dimension(*) :: x
  complex,dimension(*) :: b
!!!!!!!!        

  integer ::i,j,k,s,n1
  real,parameter :: SZERO=0.0e0
  complex,parameter :: ZEROC=(SZERO,SZERO)

  if (TRANS=='N') then
     n1=N
  else
     n1=M
  endif


!!!!!!!! Initialization

  if (beta/=ZEROC) then
     call CSCAL(n1*rhs,beta,b,1)
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
              end do
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

     elseif (TRANS=='C') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
              end do
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
              end do
           end do
        end do

     elseif (TRANS=='T') then

        do i=1,N
           do k=ia(i),ia(i+1)-1
              j=ja(k)
              do s=1,rhs
                 b(j+(s-1)*M)=b(j+(s-1)*M)+alpha*conjg(a(k))*x(i+(s-1)*N)
                 if (j/=i) b(i+(s-1)*N)=b(i+(s-1)*N)+alpha*a(k)*x(j+(s-1)*M)
              end do
           end do
        end do

     end if



  end if

end SUBROUTINE chcsrmm



SUBROUTINE sccsrmm(UPLO,TRANS,N,M,rhs,alpha,a,ia,ja,X,beta,B)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine peforms
  !
  !      b=alpha*o(A)*x+beta*b      
  !  
  !      where the NxM matrix A is stored in CSR format (a,ia,ja), 
  !      and X is Mxrhs and B is Nxrhs 
  !      
  !      SINGLE PRECISION version
  !           
  !      Note: same as ccsrmm but A is real
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
  !  alpha (input) COMPLEX SINGLE PRECISION
  !
  !  a,ia,ja (input) CSR format for the matrix A (a SINGLE PRECISION)
  !
  !  X     (input) COMPLEX SINGLE PRECISION
  !        matrix of size (Mxrhs) or (Nxrhs) if TRANS='T','C'.
  !
  !  beta  (input) COMPLEX SINGLE PRECISION
  !
  !  B     (input/output) COMPLEX SINGLE PRECISION
  !        matrix of size (Nxrhs) or (Mxrhs) if TRANS='T','C'.. On exit contains the solution.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  implicit none
  character(len=1) :: UPLO,TRANS
  integer :: N,M,rhs
  complex :: alpha,beta
  real,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex,dimension(*):: x
  complex,dimension(*) ::b
!!!!!!!!

  integer ::i,j,k,s,n1
  complex,parameter :: ZEROC=(0.0e0,0.0e0)


  if (TRANS=='N') then
     n1=N
  else
     n1=M
  endif



!!!!!!!! Initialization

  if (beta/=ZEROC) then
     call CSCAL(n1*rhs,beta,b,1)
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
              end do
           end do
        end do

     elseif ((TRANS=='T').or.(TRANS=='C')) then

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


END SUBROUTINE sccsrmm





subroutine scsr_transpose(N,M,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine transpose a double real CSR matrix A  
  !      into a new CSR matrix sA
  !
  !  N,M    (input) INTEGER
  !          The number of rows/columns of the matrices A.  N,M >= 0.
  !          The number of columns/rows of the matrices sA. N,M >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a SINGLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa SINGLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N,M
  real,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  real,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
  integer :: i,j,k,istart
  integer,dimension(:),allocatable :: tab
!!!!

  allocate(tab(M))
!!!!
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

end subroutine scsr_transpose



subroutine ccsr_transpose(N,M,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the transpose  of 
  !      a double complex CSR matrix A into a new CSR matrix sA
  !
  !
  !  N,M    (input) INTEGER
  !          The number of rows/columns of the matrices A.  N,M >= 0.
  !          The number of columns/rows of the matrices sA. N,M >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX SINGLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX SINGLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N,M
  complex,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
  integer :: i,j,k,istart
  integer,dimension(:),allocatable :: tab
!!!!

  allocate(tab(M))
!!!!
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

end subroutine ccsr_transpose




subroutine ccsr_htranspose(N,M,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine performs the transpose conjugate of 
  !      a double complex CSR matrix A into a new CSR matrix sA
  !
  !
  !  N,M    (input) INTEGER
  !          The number of rows/columns of the matrices A.  N,M >= 0.
  !          The number of columns/rows of the matrices sA. N,M >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX SINGLE PRECISION)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX SINGLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: N,M
  complex,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
  integer :: i,j,k,istart
  integer,dimension(:),allocatable :: tab
!!!!

  allocate(tab(M))
!!!!
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
     sia(i)=1 !!<< integer,dimension(:),allocatable :: tab
!!!!

     allocate(tab(M))
!!!!
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

end subroutine ccsr_htranspose




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine scsr_uplo_to_csr(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a real upper or lower CSR matrix A  
  !      into a new full CSR matrix sA
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a REAL)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa REAL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  !  include 'f90_noruntime_interface.fi'
  integer :: N
  real,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  real,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  real,dimension(:),allocatable ::sb
  integer,dimension(:),allocatable :: sib
  integer,dimension(:),allocatable :: sjb
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
!!!!
  nnza=ia(n+1)-1

  allocate(sb(nnza)) 
  allocate(sib(N+1)) 
  allocate(sjb(nnza)) 


!!! transpose
  call scsr_transpose(N,N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=SZERO
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call saddcsr(N,N,opt,SONE,a,ia,ja,SONE,sb,sib,sjb,sa,sia,sja)


  deallocate(sb)
  deallocate(sib)
  deallocate(sjb)



end subroutine scsr_uplo_to_csr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ccsr_uplo_to_csr(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a complex symetric upper or lower CSR matrix A  
  !      into a new full CSR matrix sA 
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  !  include 'f90_noruntime_interface.fi'
  integer :: N
  complex,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  complex,dimension(:),allocatable ::sb
  integer,dimension(:),allocatable :: sib
  integer,dimension(:),allocatable :: sjb
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
!!!!
  nnza=ia(n+1)-1


  allocate(sb(nnza)) 
  allocate(sib(N+1)) 
  allocate(sjb(nnza)) 


!!! transpose
  call ccsr_transpose(N,N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=ZEROC
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call caddcsr(N,N,opt,ONEC,a,ia,ja,ONEC,sb,sib,sjb,sa,sia,sja)


  deallocate(sb)
  deallocate(sib)
  deallocate(sjb)



end subroutine ccsr_uplo_to_csr






subroutine chcsr_uplo_to_csr(N,a,ia,ja,sa,sia,sja)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a complex Hermitian upper or lower CSR matrix A  
  !      into a new full CSR matrix sA 
  !
  !  N      (input) INTEGER
  !          The number of rows/columns of the matrices A,sA.  N >= 0.
  !
  !  a,ia,ja (input) CSR format for the matrix A (a COMPLEX)
  !
  !  sa,isa,jsa (output) CSR format for the matrix sA (sa COMPLEX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  !  include 'f90_noruntime_interface.fi'
  integer :: N
  complex ,dimension(*) ::a
  integer,dimension(*) :: ia
  integer,dimension(*) :: ja
  complex ,dimension(*) ::sa
  integer,dimension(*) :: sia
  integer,dimension(*) :: sja
!!!!
  integer :: i,k,nnza,info,opt
  complex ,dimension(:),allocatable ::sb
  integer,dimension(:),allocatable :: sib
  integer,dimension(:),allocatable :: sjb
  real,parameter :: SONE=1.0e0,SZERO=0.0e0
  complex ,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
!!!!
  nnza=ia(n+1)-1


  allocate(sb(nnza)) 
  allocate(sib(N+1)) 
  allocate(sjb(nnza)) 


!!! transpose
  call ccsr_htranspose(N,N,a,ia,ja,sb,sib,sjb)
!!! nullify diagonal
  do i=1,N
     do k=sib(i),sib(i+1)-1
        if (sjb(k)==i) sb(k)=ZEROC
     end do
  enddo
  !! add the two matrices to form the full matrix
  opt=0
  call caddcsr(N,N,opt,ONEC,a,ia,ja,ONEC,sb,sib,sjb,sa,sia,sja)


  deallocate(sb)
  deallocate(sib)
  deallocate(sjb)



end subroutine chcsr_uplo_to_csr




subroutine scsr_convert_upper(n,UPLO,sa,isa,jsa,ssa,sisa,sjsa)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a double real upper, lower or full CSR matrix A  
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
  !  sa,isa,jsa (input) CSR format for the matrix A (a SINGLE PRECISION)
  !
  !  ssa,sisa,sjsa (output) CSR format for the matrix sA (sa SINGLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: n
  character(len=1) :: UPLO
  real,dimension(*) :: sa,ssa
  integer,dimension(*) :: isa,jsa,sisa,sjsa
  integer :: opt,nnz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  if ((UPLO=='F').or.(UPLO=='f')) then !!! full csr to upper-csr

     opt=0
     call scsr2csr_up(opt,N,sa,isa,jsa,ssa,sisa,sjsa)

!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif ((UPLO=='U').or.(UPLO=='u')) then !! upper-csr already 
     nnz=isa(n+1)-1
     call SCOPY(nnz,sa,1,ssa,1)
     sisa(1:n+1)= isa(1:n+1)
     sjsa(1:nnz)= jsa(1:nnz)

  elseif ((UPLO=='L').or.(UPLO=='l')) then !!! lower-csr to upper-csr'
     call scsr_transpose(N,N,sa,isa,jsa,ssa,sisa,sjsa)
  endif


end subroutine scsr_convert_upper




subroutine chcsr_convert_full(n,UPLO,sa,isa,jsa,ssa,sisa,sjsa)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Purpose: This subroutine convert a single complex upper, lower or full Hermitian CSR matrix A  
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
  !  sa,isa,jsa (input) CSR format for the Hermitian matrix A (a COMPLEX SINGLE PRECISION)
  !
  !  ssa,sisa,sjsa (output) CSR format for the Hermitian matrix sA (sa COMPLEX SINGLE PRECISION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  integer :: n
  character(len=1) :: UPLO
  complex,dimension(*) :: sa,ssa
  integer,dimension(*) :: isa,jsa,sisa,sjsa
  integer :: opt,nnz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((UPLO=='F').or.(UPLO=='f')) then !!! direct copy

     nnz=isa(n+1)-1
     call CCOPY(nnz,sa,1,ssa,1)
     sisa(1:n+1)= isa(1:n+1)
     sjsa(1:nnz)= jsa(1:nnz)

!!!!!!!!!!!!!!!!!!!!!!!!!
  else!! upper-csr or lower-csr to full csr

     call chcsr_uplo_to_csr(N,sa,isa,jsa,ssa,sisa,sjsa)
  end if


end subroutine chcsr_convert_full







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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ssubmat (n,i1,i2,j1,j2,a,ja,ia,ao,jao,iao)
  !*******************************************************************************
  ! SUBMAT extracts the submatrix A(i1:i2,j1:j2) 
  !*******************************************************************************
  !a,ja,ia = original matrix in CSR format
  !ao,jao,iao = extracted matrix in CSR format
  implicit none
  integer :: i1,i2,j1,j2,n
  real,dimension(*):: a
  integer,dimension(*) :: ia,ja
  real,dimension(*) :: ao
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
end subroutine ssubmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine csubmat (n,i1,i2,j1,j2,a,ja,ia,ao,jao,iao)
  !*******************************************************************************
  ! SUBMAT extracts the submatrix A(i1:i2,j1:j2) 
  !*******************************************************************************
  !a,ja,ia = original matrix in CSR format
  !ao,jao,iao = extracted matrix in CSR format
  implicit none
  integer :: i1,i2,j1,j2,n
  complex,dimension(*):: a
  integer,dimension(*) :: ia,ja
  complex,dimension(*) :: ao
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
end subroutine csubmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









