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
!!!!!!!!!!! BANDED PRIMITIVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! DOUBLE PRECISION VERSION

!!!!! List of routines 


!DGBALU, DGBALU2  :: Approximate LU (no pivot + boosting)
!ZGBALU, ZGBALU2  :: Approximate LU (no pivot + boosting)

!DBGAUL, DGBAUL2  :: Approximate UL (no pivot + boosting)
!ZBGAUL, ZGBAUL2  :: Approximate UL (no pivot + boosting)

!DTBSM            :: Solve triangular systems AX=F or A'X=F
!ZTBSM            :: Solve triangular systems AX=F or A'X=F or conjg(A)'X=F

!DSBMM            :: Perform matrix-matrix multiplication;  C=alpha*A*B+beta*C;  matrix is structural symmetric or symmetric
!ZHBMM            :: Perform matrix-matrix multiplication;  C=alpha*A*B+beta*C;  matrix is structural symmetric or Hermitian
!ZSBMM            :: Perform matrix-matrix multiplication;  C=alpha*A*B+beta*C;  matrix is structural symmetric or symmetric
!ZGBMM            :: Perform matrix-matrix multiplication;  C=alpha*op(A)*op(B)+beta*C;  matrix is non-Hermitian
!DZGBMM           :: Same than ZGBMM but A matrix is real --> works only with DZGEMM (MKL-BLAS)- **instructions must be uncommented**

!DTBSMPY          :: modified version of DTBSM
!ZTBSMPY          :: modified version of ZTBSM




subroutine DGBALU(N,kl,ku,A,LDA,nzero,norm,info)
  !  Purpose
  !  =======
  !
  !  it computes an ("approximate") LU factorization (blocked-BLAS3) of a real   n-by-n band matrix A
  !  without using partial pivoting with row interchanges.+ diagonal boosting if any
  !
  !  This is the blocked version of the algorithm for square matrices, calling Level 3 BLAS.
  !
  !  Arguments
  !  =========
  !
  !  N      (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !  KU     (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  A      (input/output) DOUBLE PRECISION   array, dimension (LDA,N)
  !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  !
  !          On exit, contains the LU factorization: U is stored as an
  !          upper triangular band matrix with KU superdiagonals in
  !          rows 1 to KU+1. L is stored as an lower triangular band matrix
  !          with KL subdiagonals in rows KU+2 to KL+KU+1. The unit diagonal of L
  !          is not stored.
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
  !           a boost will take place if pivot <nzero*norm
  !
  ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
  !          Value of the boosting will be (+/-) eps*norm
  !          norm is a dummy argument (unused) if nzero=0.0d0
  ! INFO    (output) INTEGER
  !          = 0: successful exit
  !          -6 < INFO < 0 : if INFO = -i, the i-th argument had an illegal value
  !          =-10: internal memory allocation problem
  !          > 0: if INFO = +i, The factorization has been completed and 
  !               +i "nzero" pivots have been found (so +i boost took place) 
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================

  implicit none
  include 'f90_noruntime_interface.fi'
  integer :: N,kl,ku,LDA
  double precision,dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info

!!!!!!!!!!!!!!!!!!
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  integer :: m,i,nbm,jmin,klr,kur,j,k,il,ju,shift,info2
  double precision,dimension(:,:),pointer :: taux1,taux2
  INTEGER:: ILAENV
  integer :: info_alloc
!!!!!!!!!!!!!!!!!!!
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( N.LT.0 ) THEN
     INFO = -1
  ELSE IF( KL.LT.0 ) THEN
     INFO = -2
  ELSE IF( KU.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGBALU', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN

  !     Determine the block size for this environment
  m = ILAENV( 1, 'DGBTRF', ' ', N, N, KL, KU )
  m = MIN( m, min(kl,ku) )

  if (m==1) then
     !! use unblocked algorithm

     CALL DGBALU2(n,n,KL,KU,A,LDA,nzero,norm,info2)
     if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'DGBALU', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

  else

     !! use block algorithm
     nbm=n/m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Start recursion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! we consider the (sub)matrix
!!! A B
!!! C D
     ! we want (M1 is also the end of L1---- LU is then done on a rectangular matrix)
!!! L1       U1 K1 
!!! M1 L2       U2    
!!!
!!!!!!!! M1 (resp. K1) may be composed by rectangular part and triangular parts.

     call wallocate_2d(taux1, m, m,info_alloc) 
     call wallocate_2d(taux2, m, m,info_alloc) 
     if (info_alloc.NE.0) then
        info=-10
        RETURN
     end if


     taux1=DZERO
     taux2=DZERO

     if (n-1<=max(kl,ku)) jmin=-m+1

     DO  k=1,nbm-1

        if (n-(k)*m-1>max(kl,ku)) then

           jmin=(k-1)*m+1

           !For M1 (rectangular remaining row)
           klr=min(kl,n-jmin)-m
           if (klr<0) klr=0
           !For K1 (rectangular remaining column)
           kur=min(ku,n-jmin)-m
           if (kur<0) kur=0
           !! total row for submatrix -- L1 (including M1)
           il=min(kl,n-(jmin+m-1))+m

!!!!!!! L1U1 factorization of the rectangular (il*m) block A

           CALL  DGBALU2(il,m,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
           if (info2<0) THEN
              info=info2+1
              if (INFO.NE.0) CALL XERBLA( 'DGBALU', -INFO )
              RETURN
           ELSE
              info=info+info2
           END IF

!!!!!
!!!!!!!!  L1K1=B so find K1
!!!!!!

           !! B may be composed by a rectangular part (m*kur) and triangular part (m*m)  
           ! rectangular part (update A directly)
           if (kur>0) CALL DTRSM('L','L','N','U',m,kur,DONE,A(ku+1,jmin),kl+ku,A(ku+1-m,jmin+m),kl+ku) 
           ! triangular part (update A indirectly)   
           do i=1,m
              do j=1,i
                 taux1(i,j)=A(ku+i-m-kur+1-j,jmin+m+kur-1+j) 
              end do
              do j=i+1,m
                 taux1(i,j)=DZERO
              end do
           enddo

           CALL DTRSM('L','L','N','U',m,m,DONE,A(ku+1,jmin),kl+ku,taux1,m) 
           do j=1,m
              do i=j,m
                 A(ku+1+i-j-m-kur,jmin+m+kur+(j-1))=taux1(i,j) 
              end do
           end do

!!!!!!
!!!!!!!!  UPDATE D to close the recursion D<=D-M1K1
!!!!!!
           do i=1,m
              do j=i,m
                 taux2(i,j)=A(ku+m+klr+i+1-j,jmin-1+j)
              enddo
           enddo


!!!! we can decompose here into 4 multiplication

           ! M1 (rectangular klr*m) * K1 (rectangular m*kur)
           if ((klr>0).and.(kur>0)) call DGEMM('N','N',klr,kur,m,-DONE,A(ku+1+m,jmin),kl+ku,A(ku+1-m,jmin+m),&
	&kl+ku,DONE,A(ku+1,jmin+m),kl+ku)
           ! M1 (rectangular klr*m) * K1 (lower triangular m*m)
           if (klr>0) then
              call DGEMM('N','N',klr,m,m,-DONE,A(ku+1+m,jmin),kl+ku,taux1,m,DONE,A(ku+1-kur,jmin+m+kur),kl+ku)
           end if
           ! M1 (upper triangular m*m) * K1 (rectangular m*kur)
           if (kur>0) then
              call DGEMM('N','N',m,kur,m,-DONE,taux2,m,A(ku+1-m,jmin+m),kl+ku,DONE,A(ku+1+klr,jmin+m),kl+ku)
           end if
           ! M1 (upper triangular m*m) * K1 (lower triangular m*m)
           if (klr>0) then
              if (kur>0) then 
                 shift=-kur+klr
              else
                 shift=klr
              endif
           else
              if (kur>0) then
                 shift=-kur
              else
                 shift=0
              endif
           endif

           call DGEMM('N','N',m,m,m,-DONE,taux2,m,taux1,m,DONE,A(ku+1+shift,jmin+m+kur),kl+ku)
        end if
     end DO

!!!!!!!!!! LU on the Last block

     jmin=jmin+m
     il=n-jmin+1
     ju=il

!!!!!!! L1U1 factorization of the rectangular (il*ju) block A
     CALL  DGBALU2(il,ju,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
     if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'DGBALU', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

     call wdeallocate_2d(taux1)  
     call wdeallocate_2d(taux2)  

  end if

end subroutine DGBALU






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DGBALU2(M,N,kl,ku,A,LDA,nzero,norm,info)
  !
  !  Purpose
  !  =======
  !
  !  it computes an ("approximate") LU factorization (unblocked-BLAS2) of a real  m-by-n band matrix A
  !  without using partial pivoting with row interchanges.+ diagonal boosting if any
  !
  !  This is the unblocked version of the algorithm, calling Level 2 BLAS, inspired from LAPACK DGBTF2.
  !
  !  Arguments
  !  =========
  !
  !  M      (input) INTEGER
  !          The number of rows of the matrix A.  M >= 0.
  !  N      (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !  KU     (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  A      (input/output) DOUBLE PRECISION   array, dimension (LDA,N)
  !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
  !           a boost will take place if pivot <nzero*norm
  !
  ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
  !          Value of the boosting will be (+/-) eps*norm
  !          norm is a dummy argument (unused) if nzero=0.0
  ! INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -i, the i-th argument had an illegal value
  !          > 0: if INFO = +i, The factorization has been completed and 
  !               +i "nzero" pivots have been found (so +i boost took place) 
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================
  implicit none
  integer :: M,N,kl,ku,LDA
  double precision,dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info

!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  double precision :: pivot
  integer :: j,JU,KM
!!!!!!!!!!!!!!!!

  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSEIF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGBALU2', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN
  IF (nzero.EQ.DZERO) then
     pivot=DZERO
  else
     pivot=nzero*norm
  end IF
  JU=1
  DO  J = 1,MIN(M,N)
     if (abs(A(KU+1,J))<=pivot) then
        if (nzero.EQ.DZERO) then 
           info=-7
           CALL XERBLA( 'DGBALU2', -INFO )
           RETURN
        end IF
        A(KU+1, J)=A(KU+1,J)+sign(nzero,A(KU+1,J))*norm
        info=info+1
     end if
     KM=min(KL,M-J)
     JU=max(JU,min(J+KU,N))
     IF( KM.GT.0 ) THEN
        CALL DSCAL(KM, DONE/A(KU+1, J), A(KU+2, J), 1)
        IF( JU.GT.J )   CALL DGER(KM, JU-J, -DONE, A(KU+2, J), 1,A(KU, J+1), KL+KU, A(KU+1, J+1), KL+KU)
     end IF
  END DO

end subroutine DGBALU2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZGBALU(N,kl,ku,A,LDA,nzero,norm,info)

  !  Purpose
  !  =======
  !
  !  it computes an ("approximate") LU factorization (blocked-BLAS3) of a complex  n-by-n band matrix A
  !  without using partial pivoting with row interchanges.+ diagonal boosting if any
  !
  !  This is the blocked version of the algorithm for square matrices, calling Level 3 BLAS.
  !
  !  Arguments
  !  =========
  !
  !  N      (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !  KU     (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  A      (input/output) COMPLEX DOUBLE PRECISION   array, dimension (LDA,N)
  !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  !
  !          On exit, contains the LU factorization: U is stored as an
  !          upper triangular band matrix with KU superdiagonals in
  !          rows 1 to KU+1. L is stored as an lower triangular band matrix
  !          with KL subdiagonals in rows KU+2 to KL+KU+1. The unit diagonal of L
  !          is not stored.
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
  !           a boost will take place if pivot <nzero*norm
  !
  ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
  !          Value of the boosting will be (+/-) eps*norm
  !          norm is a dummy argument (unused) if nzero=0.0d0
  ! INFO    (output) INTEGER
  !          = 0: successful exit
  !          -6 < INFO < 0 : if INFO = -i, the i-th argument had an illegal value
  !          =-10: internal memory allocation problem
  !          > 0: if INFO = +i, The factorization has been completed and 
  !               +i "nzero" pivots have been found (so +i boost took place) 
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================

  implicit none
  include 'f90_noruntime_interface.fi'
  integer :: N,kl,ku,LDA
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info

!!!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  integer :: m,i,nbm,jmin,klr,kur,j,k,il,ju,shift,info2
  COMPLEX(kind=kind(1.0d0)),dimension(:,:),pointer :: taux1,taux2
  INTEGER:: ILAENV
  integer :: info_alloc
!!!!!!!!!!!!!!!!!!!
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( N.LT.0 ) THEN
     INFO = -1
  ELSE IF( KL.LT.0 ) THEN
     INFO = -2
  ELSE IF( KU.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZGBALU', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN

  !     Determine the block size for this environment
  m = ILAENV( 1, 'ZGBTRF', ' ', N, N, KL, KU )
  m = MIN( m, min(kl,ku) )

  if (m==1) then
     !! use unblocked algorithm

     CALL ZGBALU2(n,n,KL,KU,A,LDA,nzero,norm,info2)
     if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'ZGBALU', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

  else

     !! use block algorithm

     nbm=n/m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Start recursion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! we consider the (sub)matrix
!!! A B
!!! C D
     ! we want (M1 is also the end of L1---- LU is then done on a rectangular matrix)
!!! L1       U1 K1 
!!! M1 L2       U2    
!!!
!!!!!!!! M1 (resp. K1) may be composed by rectangular part and triangular parts.

     call wallocate_2z(taux1, m, m,info_alloc) 
     call wallocate_2z(taux2, m, m,info_alloc) 
     if (info_alloc.NE.0) then
        info=-10
        RETURN
     end if


     taux1=ZEROC
     taux2=ZEROC

     if (n-1<=max(kl,ku)) jmin=-m+1

     DO  k=1,nbm-1

        if (n-(k)*m-1>max(kl,ku)) then

           jmin=(k-1)*m+1

           !For M1 (rectangular remaining row)
           klr=min(kl,n-jmin)-m
           if (klr<0) klr=0
           !For K1 (rectangular remaining column)
           kur=min(ku,n-jmin)-m
           if (kur<0) kur=0
           !! total row for submatrix -- L1 (including M1)
           il=min(kl,n-(jmin+m-1))+m

!!!!!!! L1U1 factorization of the rectangular (il*m) block A

           CALL  ZGBALU2(il,m,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
           if (info2<0) THEN
              info=info2+1
              if (INFO.NE.0) CALL XERBLA( 'ZGBALU', -INFO )
              RETURN
           ELSE
              info=info+info2
           END IF

!!!!!
!!!!!!!!  L1K1=B so find K1
!!!!!!

           !! B may be composed by a rectangular part (m*kur) and triangular part (m*m)  
           ! rectangular part (update A directly)
           if (kur>0) CALL ZTRSM('L','L','N','U',m,kur,ONEC,A(ku+1,jmin),kl+ku,A(ku+1-m,jmin+m),kl+ku) 
           ! triangular part (update A indirectly)   
           do i=1,m
              do j=1,i
                 taux1(i,j)=A(ku+i-m-kur+1-j,jmin+m+kur-1+j) 
              end do
              do j=i+1,m
                 taux1(i,j)=ZEROC
              end do
           enddo

           CALL ZTRSM('L','L','N','U',m,m,ONEC,A(ku+1,jmin),kl+ku,taux1,m) 
           do j=1,m
              do i=j,m
                 A(ku+1+i-j-m-kur,jmin+m+kur+(j-1))=taux1(i,j) 
              end do
           end do

!!!!!!
!!!!!!!!  UPDATE D to close the recursion D<=D-M1K1
!!!!!!
           do i=1,m
              do j=i,m
                 taux2(i,j)=A(ku+m+klr+i+1-j,jmin-1+j)
              enddo
           enddo


!!!! we can decompose here into 4 multiplication

           ! M1 (rectangular klr*m) * K1 (rectangular m*kur)
           if ((klr>0).and.(kur>0)) call ZGEMM('N','N',klr,kur,m,-ONEC,A(ku+1+m,jmin),kl+ku,A(ku+1-m,jmin+m),kl+ku,&
	&ONEC,A(ku+1,jmin+m),kl+ku)
           ! M1 (rectangular klr*m) * K1 (lower triangular m*m)
           if (klr>0) then
              call ZGEMM('N','N',klr,m,m,-ONEC,A(ku+1+m,jmin),kl+ku,taux1,m,ONEC,A(ku+1-kur,jmin+m+kur),kl+ku)
           end if
           ! M1 (upper triangular m*m) * K1 (rectangular m*kur)
           if (kur>0) then
              call ZGEMM('N','N',m,kur,m,-ONEC,taux2,m,A(ku+1-m,jmin+m),kl+ku,ONEC,A(ku+1+klr,jmin+m),kl+ku)
           end if
           ! M1 (upper triangular m*m) * K1 (lower triangular m*m)
           if (klr>0) then
              if (kur>0) then 
                 shift=-kur+klr
              else
                 shift=klr
              endif
           else
              if (kur>0) then
                 shift=-kur
              else
                 shift=0
              endif
           endif

           call ZGEMM('N','N',m,m,m,-ONEC,taux2,m,taux1,m,ONEC,A(ku+1+shift,jmin+m+kur),kl+ku)
        end if
     end DO

!!!!!!!!!! LU on the Last block

     jmin=jmin+m
     il=n-jmin+1
     ju=il

!!!!!!! L1U1 factorization of the rectangular (il*ju) block A
     CALL  ZGBALU2(il,ju,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
     if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'ZGBALU', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

     call wdeallocate_2z(taux1)  
     call wdeallocate_2z(taux2)  

  end if

end subroutine ZGBALU






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZGBALU2(M,N,kl,ku,A,LDA,nzero,norm,info)
  !
  !  Purpose
  !  =======
  !
  !  it computes an ("approximate") LU factorization (unblocked-BLAS2) of a complex  m-by-n band matrix A
  !  without using partial pivoting with row interchanges.+ diagonal boosting if any
  !
  !  This is the unblocked version of the algorithm, calling Level 2 BLAS, inspired from LAPACK ZGBTF2.
  !
  !  Arguments
  !  =========
  !
  !  M      (input) INTEGER
  !          The number of rows of the matrix A.  M >= 0.
  !  N      (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !  KU     (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  A      (input/output) COMPLEX DOUBLE PRECISION   array, dimension (LDA,N)
  !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
  !           a boost will take place if pivot <nzero*norm
  !
  ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
  !          Value of the boosting will be (+/-) eps*norm
  !          norm is a dummy argument (unused) if nzero=0.0
  ! INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -i, the i-th argument had an illegal value
  !          > 0: if INFO = +i, The factorization has been completed and 
  !               +i "nzero" pivots have been found (so +i boost took place) 
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================
  implicit none
  integer :: M,N,kl,ku,LDA
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info

!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)), Parameter :: ONEC=(DONE,DZERO)
  double precision :: pivot
  integer :: j,JU,KM
!!!!!!!!!!!!!!!!

  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSEIF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZGBALU2', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN
  IF (nzero.EQ.DZERO) then
     pivot=DZERO
  else
     pivot=nzero*norm
  end IF
  JU=1
  DO  J = 1,MIN(M,N)
     if (abs(A(KU+1,J))<=pivot) then
        if (nzero.EQ.DZERO) then 
           info=-7
           CALL XERBLA( 'ZGBALU2', -INFO )
           RETURN
        end IF
        A(KU+1, J)=A(KU+1,J)+ONEC*sign(nzero,abs(A(KU+1,J)))*norm
        info=info+1
     end if
     KM=min(KL,M-J)
     JU=max(JU,min(J+KU,N))
     IF( KM.GT.0 ) THEN
        CALL ZSCAL(KM, ONEC/A(KU+1, J), A(KU+2, J), 1)
        IF( JU.GT.J )   CALL ZGERU(KM, JU-J, -ONEC, A(KU+2, J), 1,A(KU, J+1), KL+KU, A(KU+1, J+1), KL+KU)
     end IF
  END DO

end subroutine ZGBALU2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      subroutine DGBAUL(N,kl,ku,A,LDA,nzero,norm,info)
        !
        !  Purpose
        !  =======
        !
        !  it computes an ("approximate")  UL factorization (blocked-BLAS3) of a real n-by-n band matrix A
        !  without using partial pivoting with row interchanges.+ diagonal boosting if any
        !
        !  This is the blocked version of the algorithm for square matrices, calling Level 3 BLAS.
        
        !  Arguments
        !  =========
        !
        !  N      (input) INTEGER
        !          The number of columns of the matrix A.  N >= 0.
        !  KL      (input) INTEGER
        !          The number of subdiagonals within the band of A.  KL >= 0.
        !  KU     (input) INTEGER
        !          The number of superdiagonals within the band of A.  KU >= 0.
        !  A      (input/output) DOUBLE PRECISION   array, dimension (LDA,N)
        !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
        !          The j-th column of A is stored in the j-th column of the
        !          array A as follows:
        !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
        !
        !          On exit, contains the UL factorization: U is stored as an
        !          upper triangular band matrix with KU superdiagonals in
        !          rows 1 to KU. L is stored as an lower triangular band matrix
        !          with KL subdiagonals in rows KU+1 to KL+KU+1. The unit diagonal of U
        !          is not stored.
        ! LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= KL+KU+1.
        ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
        !           a boost will take place if pivot <nzero*norm
        !
        ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
        !          Value of the boosting will be (+/-) eps*norm
        !          norm is a dummy argument (unused) if nzero=0.0d0
        ! INFO    (output) INTEGER
        !          = 0: successful exit
        !          -6 < INFO < 0 : if INFO = -i, the i-th argument had an illegal value
        !          =-10: internal memory allocation problem
        !          > 0: if INFO = +i, The factorization has been completed and 
        !               +i "nzero" pivots have been found (so +i boost took place) 
        !=====================================================================
        ! Eric Polizzi 2009
        ! ====================================================================
                
  implicit none
  include 'f90_noruntime_interface.fi'
  integer :: N,kl,ku,LDA
  double precision,dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info


!!!!!!!
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  integer :: m,i,nbm,jmin,klr,kur,j,k,il,ju,shift,info2
  double precision,dimension(:,:),pointer :: taux1,taux2  
  INTEGER::  ILAENV
  integer :: info_alloc
!!!!!!!

!!!!!!!!!!!!!!!!!!!
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( N.LT.0 ) THEN
     INFO = -1
  ELSE IF( KL.LT.0 ) THEN
     INFO = -2
  ELSE IF( KU.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGBAUL', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN

        !     Determine the block size for this environment
        m = ILAENV( 1, 'DGBTRF', ' ', N, N, KL, KU )
        m = MIN( m, min(kl,ku) )
     
        if (m==1) then
           !! use unblocked algorithm

           CALL DGBAUL2(n,n,KL,KU,A,LDA,nzero,norm,info2)
 if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'DGBAUL', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

        else
           !! use block algorithm

           nbm=n/m


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Start recursion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! we consider the (sub)matrix
!!! A B
!!! C D
           ! we want (K1 is also the end of U1---- U1L1 is then done on a rectangular matrix)
!!!  U2 K1     L2      
!!!     U1     M1 L1      
!!!
!!!!!!!! M1 (resp. K1) may be composed by rectangular part and triangular parts.


 call wallocate_2d(taux1, m, m,info_alloc) 
 call wallocate_2d(taux2, m, m,info_alloc)  
 if (info_alloc.NE.0) then
        info=-10
        RETURN
     end if

           taux1=DZERO
           taux2=DZERO

           if (n-m+1<=max(kl,ku)) jmin=1

           DO  k=1,nbm-1


              if (n-k*m+1>max(kl,ku)) then

                 jmin=n-k*m+1
                 !For M1 (rectangular remaining row)
                 klr=min(kl,jmin)-m
                 if (klr<0) klr=0
                 !For K1 (rectangular remaining column)
                 kur=min(ku,jmin)-m
                 if (kur<0) kur=0
                 !! total row for submatrix -- U1 (including K1)
                 ju=m+m+kur


!!!!!!! U1L1 factorization of the rectangular (ju*m) block (D,B)

                 CALL DGBAUL2(ju,m,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
 if (info2<0) THEN
              info=info2+1
              if (INFO.NE.0) CALL XERBLA( 'DGBAUL', -INFO )
              RETURN
           ELSE
              info=info+info2
           END IF
                
!!!!!
!!!!!!!!  U1M1=C so find M1
!!!!!!
                 !! C may be composed by a rectangular part (m*klr) and upper triangular part (m*m)  
                 ! rectangular part (update A directly)
                 if (klr>0) CALL DTRSM('L','U','N','U',m,klr,DONE,A(ku+1,jmin),kl+ku,A(ku+1+klr,jmin-klr),kl+ku) 
                 ! triangular part (update A indirectly)   
                 do i=1,m
                    do j=1,i-1
                       taux1(i,j)=DZERO 
                    end do
                    do j=i,m
                       taux1(i,j)=A(ku+i+m+klr+1-j,jmin-m-klr-1+j)           
                    end do
                 end do

                 CALL DTRSM('L','U','N','U',m,m,DONE,A(ku+1,jmin),kl+ku,taux1,m) 
                 do j=1,m
                    do i=1,j
                       A(ku+1+i-j+m+klr,jmin-m-klr+(j-1))=taux1(i,j) 
                    end do
                 end do

!!!!!!
!!!!!!!!  UPDATE D to close the recursion D<=D-K1M1
!!!!!!
                 do i=1,m
                    do j=1,i
                       taux2(i,j)=A(ku+1-m-kur+i-j,jmin-1+j) 
                    end do
                 end do

!!!! we can decompose here into 4 multiplication

                 ! K1 (rectangular kur*m) * M1 (rectangular m*klr)
                 if ((klr>0).and.(kur>0)) call DGEMM('N','N',kur,klr,m,-DONE,A(ku+1-kur,jmin),kl+ku,A(ku+1+klr,jmin-klr),&
	&kl+ku,DONE,A(ku+1+klr-kur,jmin-klr),kl+ku)
                 ! K1 (rectangular kur*m) * M1 (upper triangular m*m)
                 if (kur>0) then
                    if (klr>0) then 
                       shift=m+klr-kur
                    else 
                       shift=m-kur
                    endif
                    call DGEMM('N','N',kur,m,m,-DONE,A(ku+1-kur,jmin),kl+ku,taux1,m,DONE,A(ku+1+shift,jmin-m-klr),kl+ku)
                 end if
                 ! K1 (lower triangular m*m) * M1 (rectangular m*klr)
                 if (klr>0) then
                    if (kur>0) then 
                       shift=klr-kur-m
                    else
                       shift=klr-m
                    endif
                    call DGEMM('N','N',m,klr,m,-DONE,taux2,m,A(ku+1+klr,jmin-klr),kl+ku,DONE,A(ku+1+shift,jmin-klr),kl+ku)
                 end if
                 ! M1 (upper triangular m*m) * K1 (lower triangular m*m)
                 if (klr>0) then
                    if (kur>0) then 
                       shift=klr-kur
                    else
                       shift=klr
                    endif
                 else
                    if (kur>0) then
                       shift=-kur
                    else
                       shift=0
                    endif
                 endif
                 call DGEMM('N','N',m,m,m,-DONE,taux2,m,taux1,m,DONE,A(ku+1+shift,jmin-m-klr),kl+ku)
              endif
           end DO


!!!!!!!!!! LU on the Last block
           il=jmin-1
           ju=il

           CALL DGBAUL2(ju,il,Kl,KU,A(1,1),LDA,nzero,norm,info2)
 if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'DGBAUL', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

     call wdeallocate_2d(taux1)  
     call wdeallocate_2d(taux2)  

        end if

      end subroutine DGBAUL










  SUBROUTINE DGBAUL2(M,N,kl,ku,A,LDA,nzero,norm,info)
!
!  Purpose
!  =======
!
!  it computes an "approximate" UL factorization (unblocked-BLAS2) of a real m-by-n band matrix A
!  without using partial pivoting with row interchanges.+ diagonal boosting if any
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of A (dense format)
!  N       (input) INTEGER
!          The number of columns of A (dense format) 
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!  A      (input/output) DOUBLE PRECISION   array, dimension (LDA,N)
!          On entry, the matrix A in band storage, in rows 1 to
!          KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array A as follows:
!          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
! LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= KL+KU+1.
! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
!           a boost will take place if pivot <nzero*norm
!
! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
!          Value of the boosting will be (+/-) eps*norm
!          norm is a dummy argument (unused) if nzero=0.0
! INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, The factorization has been completed and 
!               +i "nzero" pivots have been found (so +i boost took place) 
!=====================================================================
! Eric Polizzi 2009
! ====================================================================
implicit none
        INTEGER :: KL,KU,M,N,LDA
        DOUBLE PRECISION,dimension(LDA,*) ::   A
        double precision :: norm,nzero
        integer :: info
!!!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  double precision :: pivot
  integer :: J,JU,KM

    
  !     Test the input parameters.
  !
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSEIF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGBAUL2', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN
  IF (nzero.EQ.DZERO) then
     pivot=DZERO
  else
     pivot=nzero*norm
  end IF

        JU=1
        DO  J = min(M,N),1,-1
           if (abs(A(KU+1,J))<=pivot) then
           if (nzero.EQ.DZERO) then 
           info=-7
           CALL XERBLA( 'DGBAUL2', -INFO )
           RETURN
        end IF
              A(KU+1, J)=A(KU+1,J)+sign(nzero,A(KU+1,J))*norm            
              info=info+1
           end if
         KM=min(KU,J-1+M-N)
           IF( KM.GT.0 ) THEN
              CALL DSCAL(KM, DONE/A(KU+1, J), A(KU+1-KM, J), 1)
              JU= MIN(KL,J-1)
              CALL DGER(KM, JU, -DONE, A(KU+1-KM, J), 1,A(KU+JU+1, J-JU), KL+KU, A(JU+KU+1-KM, J-JU), KL+KU)
           end IF
        END DO

      end subroutine DGBAUL2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      subroutine ZGBAUL(N,kl,ku,A,LDA,nzero,norm,info)
        !
        !  Purpose
        !  =======
        !
        !  it computes an ("approximate")  UL factorization (blocked-BLAS3) of a complex n-by-n band matrix A
        !  without using partial pivoting with row interchanges.+ diagonal boosting if any
        !
        !  This is the blocked version of the algorithm for square matrices, calling Level 3 BLAS.
        
        !  Arguments
        !  =========
        !
        !  N      (input) INTEGER
        !          The number of columns of the matrix A.  N >= 0.
        !  KL      (input) INTEGER
        !          The number of subdiagonals within the band of A.  KL >= 0.
        !  KU     (input) INTEGER
        !          The number of superdiagonals within the band of A.  KU >= 0.
        !  A      (input/output) COMPLEX DOUBLE PRECISION   array, dimension (LDA,N)
        !          On entry, the matrix A in band storage, in rows 1 to KL+KU+1;       
        !          The j-th column of A is stored in the j-th column of the
        !          array A as follows:
        !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
        !
        !          On exit, contains the UL factorization: U is stored as an
        !          upper triangular band matrix with KU superdiagonals in
        !          rows 1 to KU. L is stored as an lower triangular band matrix
        !          with KL subdiagonals in rows KU+1 to KL+KU+1. The unit diagonal of U
        !          is not stored.
        ! LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= KL+KU+1.
        ! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
        !           a boost will take place if pivot <nzero*norm
        !
        ! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
        !          Value of the boosting will be (+/-) eps*norm
        !          norm is a dummy argument (unused) if nzero=0.0d0
        ! INFO    (output) INTEGER
        !          = 0: successful exit
        !          -6 < INFO < 0 : if INFO = -i, the i-th argument had an illegal value
        !          =-10: internal memory allocation problem
        !          > 0: if INFO = +i, The factorization has been completed and 
        !               +i "nzero" pivots have been found (so +i boost took place) 
        !=====================================================================
        ! Eric Polizzi 2009
        ! ====================================================================
                
  implicit none
  include 'f90_noruntime_interface.fi'
  integer :: N,kl,ku,LDA
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*) :: A
  double precision:: norm,nzero
  integer :: info


!!!!!!!
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
  integer :: m,i,nbm,jmin,klr,kur,j,k,il,ju,shift,info2
  COMPLEX(kind=kind(1.0d0)),dimension(:,:),pointer :: taux1,taux2  
  INTEGER::  ILAENV
  integer :: info_alloc
!!!!!!!

!!!!!!!!!!!!!!!!!!!
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( N.LT.0 ) THEN
     INFO = -1
  ELSE IF( KL.LT.0 ) THEN
     INFO = -2
  ELSE IF( KU.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZGBAUL', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN

        !     Determine the block size for this environment
        m = ILAENV( 1, 'ZGBTRF', ' ', N, N, KL, KU )
        m = MIN( m, min(kl,ku) )
     
        if (m==1) then
           !! use unblocked algorithm

           CALL ZGBAUL2(n,n,KL,KU,A,LDA,nzero,norm,info2)
 if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'ZGBAUL', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

        else
           !! use block algorithm

           nbm=n/m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Start recursion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! we consider the (sub)matrix
!!! A B
!!! C D
           ! we want (K1 is also the end of U1---- U1L1 is then done on a rectangular matrix)
!!!  U2 K1     L2      
!!!     U1     M1 L1      
!!!
!!!!!!!! M1 (resp. K1) may be composed by rectangular part and triangular parts.


 call wallocate_2z(taux1, m, m,info_alloc) 
 call wallocate_2z(taux2, m, m,info_alloc)  
 if (info_alloc.NE.0) then
        info=-10
        RETURN
     end if

           taux1=ZEROC
           taux2=ZEROC

           if (n-m+1<=max(kl,ku)) jmin=1

           DO  k=1,nbm-1


              if (n-k*m+1>max(kl,ku)) then

                 jmin=n-k*m+1
                 !For M1 (rectangular remaining row)
                 klr=min(kl,jmin)-m
                 if (klr<0) klr=0
                 !For K1 (rectangular remaining column)
                 kur=min(ku,jmin)-m
                 if (kur<0) kur=0
                 !! total row for submatrix -- U1 (including K1)
                 ju=m+m+kur


!!!!!!! U1L1 factorization of the rectangular (ju*m) block (D,B)

                 CALL ZGBAUL2(ju,m,KL,KU,A(1,jmin),LDA,nzero,norm,info2)
 if (info2<0) THEN
              info=info2+1
              if (INFO.NE.0) CALL XERBLA( 'ZGBAUL', -INFO )
              RETURN
           ELSE
              info=info+info2
           END IF
                
!!!!!
!!!!!!!!  U1M1=C so find M1
!!!!!!
                 !! C may be composed by a rectangular part (m*klr) and upper triangular part (m*m)  
                 ! rectangular part (update A directly)
                 if (klr>0) CALL ZTRSM('L','U','N','U',m,klr,ONEC,A(ku+1,jmin),kl+ku,A(ku+1+klr,jmin-klr),kl+ku) 
                 ! triangular part (update A indirectly)   
                 do i=1,m
                    do j=1,i-1
                       taux1(i,j)=ZEROC
                    end do
                    do j=i,m
                       taux1(i,j)=A(ku+i+m+klr+1-j,jmin-m-klr-1+j)           
                    end do
                 end do

                 CALL ZTRSM('L','U','N','U',m,m,ONEC,A(ku+1,jmin),kl+ku,taux1,m) 
                 do j=1,m
                    do i=1,j
                       A(ku+1+i-j+m+klr,jmin-m-klr+(j-1))=taux1(i,j) 
                    end do
                 end do

!!!!!!
!!!!!!!!  UPDATE D to close the recursion D<=D-K1M1
!!!!!!
                 do i=1,m
                    do j=1,i
                       taux2(i,j)=A(ku+1-m-kur+i-j,jmin-1+j) 
                    end do
                 end do

!!!! we can decompose here into 4 multiplication

                 ! K1 (rectangular kur*m) * M1 (rectangular m*klr)
                 if ((klr>0).and.(kur>0)) call ZGEMM('N','N',kur,klr,m,-ONEC,A(ku+1-kur,jmin),kl+ku,A(ku+1+klr,jmin-klr),&
	&kl+ku,ONEC,A(ku+1+klr-kur,jmin-klr),kl+ku)
                 ! K1 (rectangular kur*m) * M1 (upper triangular m*m)
                 if (kur>0) then
                    if (klr>0) then 
                       shift=m+klr-kur
                    else 
                       shift=m-kur
                    endif
                    call ZGEMM('N','N',kur,m,m,-ONEC,A(ku+1-kur,jmin),kl+ku,taux1,m,ONEC,A(ku+1+shift,jmin-m-klr),kl+ku)
                 end if
                 ! K1 (lower triangular m*m) * M1 (rectangular m*klr)
                 if (klr>0) then
                    if (kur>0) then 
                       shift=klr-kur-m
                    else
                       shift=klr-m
                    endif
                    call ZGEMM('N','N',m,klr,m,-ONEC,taux2,m,A(ku+1+klr,jmin-klr),kl+ku,ONEC,A(ku+1+shift,jmin-klr),kl+ku)
                 end if
                 ! M1 (upper triangular m*m) * K1 (lower triangular m*m)
                 if (klr>0) then
                    if (kur>0) then 
                       shift=klr-kur
                    else
                       shift=klr
                    endif
                 else
                    if (kur>0) then
                       shift=-kur
                    else
                       shift=0
                    endif
                 endif
                 call ZGEMM('N','N',m,m,m,-ONEC,taux2,m,taux1,m,ONEC,A(ku+1+shift,jmin-m-klr),kl+ku)
              endif
           end DO


!!!!!!!!!! LU on the Last block
           il=jmin-1
           ju=il

           CALL ZGBAUL2(ju,il,Kl,KU,A(1,1),LDA,nzero,norm,info2)
 if (info2<0) THEN
        info=info2+1
        if (INFO.NE.0) CALL XERBLA( 'ZGBAUL', -INFO )
        RETURN
     ELSE
        info=info+info2
     END IF

     call wdeallocate_2z(taux1)  
     call wdeallocate_2z(taux2)  

        end if

      end subroutine ZGBAUL






  SUBROUTINE ZGBAUL2(M,N,kl,ku,A,LDA,nzero,norm,info)
!
!  Purpose
!  =======
!
!  it computes an "approximate" UL factorization (unblocked-BLAS2) of a complex m-by-n band matrix A
!  without using partial pivoting with row interchanges.+ diagonal boosting if any
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of A (dense format)
!  N       (input) INTEGER
!          The number of columns of A (dense format) 
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!  A      (input/output) COMPLEX DOUBLE PRECISION   array, dimension (LDA,N)
!          On entry, the matrix A in band storage, in rows 1 to
!          KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array A as follows:
!          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
! LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= KL+KU+1.
! nzero    (input) DOUBLE PRECISION  value of the new zero for the pivot
!           a boost will take place if pivot <nzero*norm
!
! norm    (input) DOUBLE PRECISION  norm of the matrix (for example norm 1) 
!          Value of the boosting will be (+/-) eps*norm
!          norm is a dummy argument (unused) if nzero=0.0
! INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, The factorization has been completed and 
!               +i "nzero" pivots have been found (so +i boost took place) 
!=====================================================================
! Eric Polizzi 2009
! ====================================================================
implicit none
        INTEGER :: KL,KU,M,N,LDA
         COMPLEX(kind=kind(1.0d0)),dimension(LDA,*) ::   A
        double precision :: norm,nzero
        integer :: info
!!!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)), Parameter :: ONEC=(DONE,DZERO)
  double precision :: pivot
  integer :: J,JU,KM

    
  !     Test the input parameters.
  !
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSEIF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDA.LT.KL+KU+1 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZGBAUL2', -INFO )
     RETURN
  END IF

  !     Quick return if possible
  IF( N.EQ.0) RETURN
  IF (nzero.EQ.DZERO) then
     pivot=DZERO
  else
     pivot=nzero*norm
  end IF

        JU=1
        DO  J = min(M,N),1,-1
           if (abs(A(KU+1,J))<=pivot) then
           if (nzero.EQ.DZERO) then 
           info=-7
           CALL XERBLA( 'ZGBAUL2', -INFO )
           RETURN
        end IF
              A(KU+1, J)=A(KU+1,J)+ONEC*sign(nzero,abs(A(KU+1,J)))*norm            
              info=info+1
           end if
         KM=min(KU,J-1+M-N)
           IF( KM.GT.0 ) THEN
              CALL ZSCAL(KM, ONEC/A(KU+1, J), A(KU+1-KM, J), 1)
              JU= MIN(KL,J-1)
              CALL ZGERU(KM, JU, -ONEC, A(KU+1-KM, J), 1,A(KU+JU+1, J-JU), KL+KU, A(JU+KU+1-KM, J-JU), KL+KU)
           end IF
        END DO

      end subroutine ZGBAUL2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DTBSM(UPLO,TRANS,DIAG,N,rhs,kd,A,LDA,B,LDB)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Purpose:solves one of the systems of equations in real double precision
  !
  !     A*x = b,   or   A'*x = b,  
  !
  !  where b and x are n by rhs element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular band matrix, with ( kd + 1)
  !  diagonals.
  !
  !  It is a  Blocked version (BLAS3) of the routine DTBSV which allow blocking with 
  !  multiple right-hand sides      
  !
  ! Arguments: 
  !       UPLO (input) character: if 'L' or 'l' lower triangular , if 'U' or 'u' upper triangular
  !       TRANS (input) character:specifies the equations to be solved as follows:
  !              TRANS = 'N' or 'n'   A*x = b.
  !              TRANS = 'T' or 't'   A'*x = b.
  !       DIAG (input) character: if 'N' or 'n' non-unit diagonal, if 'U' or 'u' unit diagonal
  !       N (input) integer: order of matrix A
  !       rhs (input) integer :number of right-hand-sides for vector B        
  !       kd (input) integer: lower bandwidth of the matrix A if UPLO='L', or
  !                           upper bandwidth of the matrix A if UPLO='U'
  !       A (input) double precision array, banded input matrix 
  !           Note: same format than the one used in BLAS DTBSV
  !
  !           Before entry with UPLO = 'U' or 'u', the leading ( kd + 1 )
  !           by n part of the array A must contain the upper triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row
  !           ( kd + 1 ) of the array, the first super-diagonal starting at
  !           position 2 in row kd, and so on. The top left kd by kd triangle
  !           of the array A is not referenced.
  !
  !           Before entry with UPLO = 'L' or 'l', the leading ( kd + 1 )
  !           by n part of the array A must contain the lower triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row 1 of
  !           the array, the first sub-diagonal starting at position 1 in
  !           row 2, and so on. The bottom right kd by kd triangle of the
  !           array A is not referenced.
  !           
  !       LDA (input) integer : The leading dimension of the array A.  LDA >= kd+1.
  !       B (input/output) double precision 2D array: right hand side contains the solution on exit 
  !       LDB (input) integer : The leading dimension of the array B.  LDB >= N.
  !
  !======================================================================
  ! Eric Polizzi 2009
  !======================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  integer :: N,LDA,LDB,rhs
  DOUBLE PRECISION,dimension(LDA,*),intent(in) ::   A
  character(len=1)  :: UPLO,DIAG,TRANS
  integer :: kd
  DOUBLE PRECISION,dimension(LDB,*) ::   B

!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,parameter :: DONE=1.0d0
  character(len=1) :: TRANSB,SIDE,UPLO2,DIAG2
  integer :: k,ijinit,i
  integer :: m,m1,max
  integer :: i1,i2,i3,i4,i5,i6,j1,j2,j3,j4,s1,s2
  double precision, pointer, dimension(:,:) :: aux
  integer :: infoloc
  integer :: ldaux
!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_2d(aux,kd,rhs,infoloc)


  !! Size of the Recursive diagonal block and initialization of UPLO2
  !! for triangular off-diagonal blocks
  m=min(kd,n)


!!!!!!!!!!!!!! Case n<kd
if (n<kd) then ! treat as dense
  if ((UPLO=='L').or.(UPLO=='l')) then
call DTRSM('L',UPLO,TRANS,DIAG,N,rhs,DONE,A(1,1),LDA-1,B,LDB)
  elseif ((UPLO=='U').or.(UPLO=='u')) then
call DTRSM('L',UPLO,TRANS,DIAG,N,rhs,DONE,A(kd+1,1),LDA-1,B,LDB)
endif
return
endif
!!!!!!!!!!!!


  if ((UPLO=='L').or.(UPLO=='l')) then
     UPLO2='U'
  elseif ((UPLO=='U').or.(UPLO=='u')) then
     UPLO2='L'
  end if

  !!initialization of DIAG2,SIDE
  DIAG2='N'
  SIDE='L'


  !! how many submatrices (m*m) in A
  max=n/kd
  !! size of the first diagonal block before recursion
  m1=n-max*kd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Initialisation before the loop
!!!! only once for the first block decomposition with size m1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (m1/=0) then !! if there is a first block

     ldaux=kd

     TRANSB = 'N'

     if ((TRANS=='N').or.(TRANS=='n')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=1
           i1=1
           j3=j1
           i3=i1+m1
           j4=j1
           i4=kd+1

           i5=j1+m1 !rhs
           i6=i5+(kd-m1)

        elseif((UPLO=='U').or.(UPLO=='u')) then
           j1=n-m1+1                
           i1=kd+1
           j3=j1
           i3=m1+1
           j4=j1
           i4=1

           i5=j1-(kd-m1)
           i6=i5-m1
        end if

     elseif ((TRANS=='T').or.(TRANS=='t')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=n-m1+1
           i1=1
           j3=n-kd+1
           i3=kd+1-m1  
           j4=j1-kd
           i4=kd+1

           i5=j1-(kd-m1)
           i6=i5-m1
        elseif ((UPLO=='U').or.(UPLO=='u')) then
           j1=1
           i1=kd+1
           j3=j1+m1
           i3=i1-m1
           j4=kd+1
           i4=1

           i5=j1+m1
           i6=i5+(kd-m1)
        end if
     end if

!!!!!! Solve using blas3 the initial block
     call DTRSM ( SIDE, UPLO, TRANS, DIAG, m1, rhs, DONE, A(i1,j1), LDA-1,B(j1,1), LDB )
!!!!!! update the right hand side using blas3 in 2 steps
     ! step 1
     call DGEMM ( TRANS, TRANSB, kd-m1, rhs, m1, -DONE, A(i3,j3), LDA-1, B(j1,1), LDB,DONE, B(i5,1), LDB ) 
     ! step 2
     do i=1,rhs
        call DCOPY(m1,B(j1,i),1,aux(1,i),1)
     end do
     call DTRMM ( SIDE, UPLO2, TRANS, DIAG2, m1, rhs, DONE, A(i4,j4), LDA-1,aux,ldaux ) 
     do i=1,rhs
        call DAXPY(m1,-DONE,aux(1,i),1,B(i6,i),1)
     end do

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! starting point origin !! new origin and other initialization before loop
  if ((TRANS=='N').or.(TRANS=='n')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=m1+1
        i1=1
        i2=kd+1
        s1=1
        s2=0
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=n-m1-kd+1
        i1=kd+1
        i2=1
        s1=-1
        s2=0
     end if
  elseif ((TRANS=='T').or.(TRANS=='t')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=n-m1-kd+1
        i1=1
        i2=kd+1
        s1=-1
        s2=-1
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=m1+1
        i1=kd+1
        i2=1
        s1=1
        s2=1
     end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! LOOP for all the submatrices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO k=1,max
     j1=ijinit+s1*(k-1)*m
     j2=j1+s2*kd
     i3=j1+s1*kd
     !! Solve using blas3 the initial bloc
     call DTRSM ( SIDE, UPLO, TRANS, DIAG, kd, rhs, DONE, A(i1,j1),LDA-1,B(j1,1), LDB ) 
     !! update the right hand side using blas3
     if (k/=max) then
        do i=1,rhs
           call DCOPY(kd,B(j1,i),1,aux(1,i),1)
        end do
        call DTRMM ( SIDE, UPLO2, TRANS, DIAG2, kd, rhs, DONE, A(i2,j2), LDA-1,aux,kd ) 
        do i=1,rhs
           call DAXPY(kd,-DONE,aux(1,i),1,B(i3,i),1)
        end do
     end if
  END DO


  call wdeallocate_2d(aux)

end subroutine DTBSM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine ZTBSM(UPLO,TRANS,DIAG,N,rhs,kd,A,LDA,B,LDB)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Purpose:solves one of the systems of equations in complex double precision
  !
  !     A*x = b,   or   A'*x = b,  or   conjg( A' )*x = b,
  !
  !  where b and x are n by rhs element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular band matrix, with ( kd + 1)
  !  diagonals.
  !
  !  It is a  Blocked version (BLAS3) of the routine DTBSV which allow blocking with 
  !  multiple right-hand sides      
  !
  ! Arguments: 
  !       UPLO (input) character: if 'L' or 'l' lower triangular , if 'U' or 'u' upper triangular
  !       TRANS (input) character:specifies the equations to be solved as follows:
  !              TRANS = 'N' or 'n'   A*x = b.
  !              TRANS = 'T' or 't'   A'*x = b.
  !              TRANS = 'C' or 'c'   conjg( A' )*x = b.
  !       DIAG (input) character: if 'N' or 'n' non-unit diagonal, if 'U' or 'u' unit diagonal
  !       N (input) integer: order of matrix A
  !       rhs (input) integer :number of right-hand-sides for vector B        
  !       kd (input) integer: lower bandwidth of the matrix A if UPLO='L', or
  !                           upper bandwidth of the matrix A if UPLO='U'
  !       A (input) complex double precision array, banded input matrix 
  !           Note: same format than the one used in BLAS DTBSV
  !
  !           Before entry with UPLO = 'U' or 'u', the leading ( kd + 1 )
  !           by n part of the array A must contain the upper triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row
  !           ( kd + 1 ) of the array, the first super-diagonal starting at
  !           position 2 in row kd, and so on. The top left kd by kd triangle
  !           of the array A is not referenced.
  !
  !           Before entry with UPLO = 'L' or 'l', the leading ( kd + 1 )
  !           by n part of the array A must contain the lower triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row 1 of
  !           the array, the first sub-diagonal starting at position 1 in
  !           row 2, and so on. The bottom right kd by kd triangle of the
  !           array A is not referenced.
  !           
  !       LDA (input) integer : The leading dimension of the array A.  LDA >= kd+1.
  !       B (input/output) complex double precision 2D array: right hand side contains the solution on exit 
  !       LDB (input) integer : The leading dimension of the array B.  LDB >= N.
  !
  !======================================================================
  ! Eric Polizzi 2009
  !======================================================================
  implicit none
    include 'f90_noruntime_interface.fi'
  integer :: N,LDA,LDB,rhs
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*),intent(in) ::   A
  character(len=1)  :: UPLO,DIAG,TRANS
  integer :: kd
  COMPLEX(kind=kind(1.0d0)),dimension(LDB,*) ::   B
  integer :: infoloc
!!!!!!!!!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO)
  character(len=1) :: TRANSB,SIDE,UPLO2,DIAG2
  integer :: k,ijinit,i
  integer :: m,m1,max
  integer :: i1,i2,i3,i4,i5,i6,j1,j2,j3,j4,s1,s2
 ! COMPLEX(kind=kind(1.0d0)),dimension(kd,rhs) :: aux
  COMPLEX(kind=kind(1.0d0)),pointer,dimension(:,:) :: aux
  integer :: ldaux
!!!!!!!!!!!!!!!!!!!!!!!!

  call wallocate_2z(aux,kd,rhs,infoloc)
  !! Size of the Recursive diagonal block and initialization of UPLO2
  !! for triangular off-diagonal blocks
  m=min(kd,n)



!!!!!!!!!!!!!! Case n<kd
if (n<kd) then ! treat as dense
  if ((UPLO=='L').or.(UPLO=='l')) then
call ZTRSM('L',UPLO,TRANS,DIAG,N,rhs,ONEC,A(1,1),LDA-1,B,LDB)
  elseif ((UPLO=='U').or.(UPLO=='u')) then
call ZTRSM('L',UPLO,TRANS,DIAG,N,rhs,ONEC,A(kd+1,1),LDA-1,B,LDB)
endif
return
endif
!!!!!!!!!!!!!!!!!!


  if ((UPLO=='L').or.(UPLO=='l')) then
     UPLO2='U'
  elseif ((UPLO=='U').or.(UPLO=='u')) then
     UPLO2='L'
  end if

  !!initialization of DIAG2,SIDE
  DIAG2='N'
  SIDE='L'

  !! how many submatrices (m*m) in A
  max=n/kd
  !! size of the first diagonal block before recursion
  m1=n-max*kd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Initialisation before the loop
!!!! only once for the first block decomposition with size m1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (m1/=0) then !! if there is a first block

     ldaux=kd

     TRANSB = 'N'

     if ((TRANS=='N').or.(TRANS=='n')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=1
           i1=1
           j3=j1
           i3=i1+m1
           j4=j1
           i4=kd+1

           i5=j1+m1 !rhs
           i6=i5+(kd-m1)

        elseif((UPLO=='U').or.(UPLO=='u')) then
           j1=n-m1+1                
           i1=kd+1
           j3=j1
           i3=m1+1
           j4=j1
           i4=1

           i5=j1-(kd-m1)
           i6=i5-m1
        end if

     elseif ((TRANS=='T').or.(TRANS=='t').or.(TRANS=='C').or.(TRANS=='c')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=n-m1+1
           i1=1
           j3=n-kd+1
           i3=kd+1-m1  
           j4=j1-kd
           i4=kd+1

           i5=j1-(kd-m1)
           i6=i5-m1
        elseif ((UPLO=='U').or.(UPLO=='u')) then
           j1=1
           i1=kd+1
           j3=j1+m1
           i3=i1-m1
           j4=kd+1
           i4=1

           i5=j1+m1
           i6=i5+(kd-m1)
        end if
     end if

!!!!!! Solve using blas3 the initial block
     call ZTRSM ( SIDE, UPLO, TRANS, DIAG, m1, rhs, ONEC, A(i1,j1), LDA-1,B(j1,1), LDB )
!!!!!! update the right hand side using blas3 in 2 steps
     ! step 1
     call ZGEMM ( TRANS, TRANSB, kd-m1, rhs, m1, -ONEC, A(i3,j3), LDA-1, B(j1,1), LDB,ONEC, B(i5,1), LDB ) 
     ! step 2
     do i=1,rhs
        call ZCOPY(m1,B(j1,i),1,aux(1,i),1)
     end do
     ! call ZLACPY( 'F', m1, rhs, B(j1,1), N, aux(1,1), ldaux )
     call ZTRMM ( SIDE, UPLO2, TRANS, DIAG2, m1, rhs, ONEC, A(i4,j4), LDA-1,aux,ldaux ) 
     do i=1,rhs
        call ZAXPY(m1,-ONEC,aux(1,i),1,B(i6,i),1)
     end do

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! starting point origin !! new origin and other initialization before loop
  if ((TRANS=='N').or.(TRANS=='n')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=m1+1
        i1=1
        i2=kd+1
        s1=1
        s2=0
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=n-m1-kd+1
        i1=kd+1
        i2=1
        s1=-1
        s2=0
     end if
  elseif ((TRANS=='T').or.(TRANS=='t').or.(TRANS=='C').or.(TRANS=='c')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=n-m1-kd+1
        i1=1
        i2=kd+1
        s1=-1
        s2=-1
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=m1+1
        i1=kd+1
        i2=1
        s1=1
        s2=1
     end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! LOOP for all the submatrices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  DO k=1,max
     j1=ijinit+s1*(k-1)*m
     j2=j1+s2*kd
     i3=j1+s1*kd
     !! Solve using blas3 the initial bloc
     call ZTRSM ( SIDE, UPLO, TRANS, DIAG, kd, rhs, ONEC, A(i1,j1),LDA-1,B(j1,1), LDB ) 
     !! update the right hand side using blas3
     if (k/=max) then
        do i=1,rhs
           call ZCOPY(kd,B(j1,i),1,aux(1,i),1)
        end do
        call ZTRMM ( SIDE, UPLO2, TRANS, DIAG2, kd, rhs, -ONEC, A(i2,j2), LDA-1,aux,kd ) 
        do i=1,rhs
           call ZAXPY(kd,ONEC,aux(1,i),1,B(i3,i),1)
        end do
     end if
  END DO
  call wdeallocate_2z(aux)

end subroutine ZTBSM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine DSBMM(UPLO,n,rhs,kd,alpha,A,LDA,B,LDB,beta,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose
  !  =======
  !  Perform the multiplication:
  !                              C=alpha*A*B+beta*C
  !  where alpha and beta are double real scalars, x and y are double real vectors and A is an
  !  n-by-n (structural symmetric) double real band matrix, with kd sub-diagonals and kd super-diagonals,
  !  B and C are n-by-rhs matrices.
  !
  !  It is a  Blocked version (BLAS3) of the routine DSBMV which allow blocking with 
  !  multiple right-hand sides  (with  generalization to structural symetric case).    
  !
  !
  !  Arguments
  !  =========
  !
  ! UPLO (input) character: specifies whether the upper or lower
  !                         triangular part of the band matrix A is being supplied as follows:
  !
  !             UPLO = 'U' or 'u'   The upper triangular part of A is
  !                                 being supplied. A is supposed symmetric.
  !
  !             UPLO = 'L' or 'l'   The lower triangular part of A is
  !                                 being supplied. A is supposed symmetric.
  !
  !             UPLO = 'F' or 'f'   All element of A are supplied. A is supposed structural symmetric.
  !
  !
  !  N      (input) INTEGER
  !          The number of row and columns of the matrix A.  
  !          The number of rows of matrices B and C.   N >= 0.
  !  rhs    (input) INTEGER
  !          The number of columns of matrices B and C
  !  KD      (input) INTEGER
  !          The number of subdiagonals kl=kd and superdiagonals ku=kd within the band of A.  KD >= 0.
  ! alpha   (input) DOUBLE PRECISION
  !  A      (input) DOUBLE PRECISION  
  !          array, dimension (LDA,N); the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! B       (input) DOUBLE PRECISION - 2D array
  ! LDB     (input) INTEGER      
  !          The leading dimension of the array B.  LDB >= N.
  ! beta    (input) DOUBLE PRECISION
  ! C       (input/output) DOUBLE PRECISION - 2D array
  !         On exit, contains the solutions 
  ! LDC     (input) INTEGER      
  !          The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! Eric Polizzi 2009
  !======================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  CHARACTER(len=1) :: UPLO
  INTEGER :: N, RHS, KD,LDA,LDB,LDC
  DOUBLE PRECISION :: alpha,beta
  DOUBLE PRECISION,dimension(LDA,*) ::  A
  DOUBLE PRECISION,dimension(LDB,*) ::  B
  DOUBLE PRECISION,dimension(LDC,*) ::  C
!!!!!!!!!!!!!!
  double precision,dimension(:,:),pointer :: helper
  integer :: info_alloc
  logical :: low,up,full
  integer :: bl, i , k , tmp, imin,jmin,iminb,iminc,imin2,jmin2
  character(len=1) :: UPLO2,DIAG,TRANSA,TRANSB, SIDE
  double precision,parameter :: DONE=1.0d0, DZERO=0.0d0


  low=.false.
  up=.false.
  full=.true.
  if ((UPLO=='L').or.(UPLO=='l')) then
     low=.true.
     full=.false.
  end if
  if ((UPLO=='U').or.(UPLO=='u')) then
     up=.true.
     full=.false.
  end if



  bl = n/kd
  tmp = n - bl*kd



  call wallocate_2d(helper,kd,rhs,info_alloc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! diagonal blocks first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TRANSA='N'
  TRANSB='N'

  SIDE= 'L'


  if (low) then
     imin=1
  else !up or full
     imin=kd+1
  end if

  if (tmp/=0) then
!!!! last block of size tmp
     jmin=bl*kd+1
     if (full) then
        call DGEMM(TRANSA,TRANSB,tmp,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(jmin,1), LDB, beta, C(jmin,1), LDC)
     else ! low or up
        call DSYMM (SIDE, UPLO, tmp, rhs, ALPHA, A(imin,jmin), LDA-1, B(jmin,1), LDB, BETA,C(jmin,1), LDC)
     end if
  end if


  if (bl==0) return

!!!! all other blocks of size kd
  if (bl>0) then
     do k=1,bl
        jmin=(k-1)*kd+1
        if (full) then
           call DGEMM(TRANSA,TRANSB,kd,rhs,kd, alpha, A(imin,jmin), LDA-1, B(jmin,1), LDB, beta, C(jmin,1), LDC)
        else ! low or up
           call DSYMM (SIDE, UPLO, kd, rhs, ALPHA, A(imin,jmin), LDA-1, B(jmin,1), LDB, BETA,C(jmin,1), LDC)
        end if
     end do
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Lower Subdiagonal blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SIDE='L'
  DIAG='N'

  if (full) then
     imin=2*kd+1
     TRANSA='N'
     UPLO2='U'
  elseif (low) then
     imin=kd+1
     TRANSA='N'
     UPLO2='U'
  elseif (up) then
     imin=1
     TRANSA='T'
     UPLO2='L'
  end if

  do k=1,bl-1
     if (full) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
     elseif (up) then
        jmin=k*kd+1
        iminb=(k-1)*kd+1
        iminc=iminb+kd
     end if
     do i=1,rhs
        call DCOPY(kd,B(iminb,i),1,helper(1,i),1)
     end do
     call DTRMM(SIDE,UPLO2,TRANSA,DIAG,kd,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)
     do i=1,rhs
        call DAXPY(kd,DONE,helper(1,i),1,C(iminc,i),1)
     end do
  end do

!!!! final block (triangular+rectangular)
  if (tmp/=0) then
     k=bl
     if (full) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (up) then
        jmin=k*kd+1
        iminb=(k-1)*kd+1
        iminc=iminb+kd
        imin2=imin+tmp
        jmin2=jmin
     end if

     ! triangular
     do i=1,rhs
        call DCOPY(tmp,B(iminb,i),1,helper(1,i),1)
     end do

     call DTRMM(SIDE,UPLO2,TRANSA,DIAG,tmp,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call DAXPY(tmp,DONE,helper(1,i),1,C(iminc,i),1)
     end do
     ! rectangular

     call DGEMM(TRANSA,TRANSB,tmp,rhs,kd-tmp, alpha, A(imin2,jmin2), LDA-1, B(iminb+tmp,1), LDB, DONE, C(iminc,1), LDC)


  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Upper Superdiagonal blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SIDE='L'
  DIAG='N'

  if (full) then
     imin=1
     TRANSA='N'
     UPLO2='L'
  elseif (low) then
     imin=kd+1
     TRANSA='T'
     UPLO2='U'
  elseif (up) then
     imin=1
     TRANSA='N'
     UPLO2='L'
  end if


  do k=1,bl-1
     if (full) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin+kd
        iminc=iminb-kd
     elseif (up) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
     end if
     do i=1,rhs
        call DCOPY(kd,B(iminb,i),1,helper(1,i),1)
     end do

     call DTRMM(SIDE,UPLO2,TRANSA,DIAG,kd,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call DAXPY(kd,DONE,helper(1,i),1,C(iminc,i),1)
     end do

  end do



!!!! final block (triangular+rectangular)
  if (tmp/=0) then
     k=bl
     if (full) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
        imin2=imin+tmp
        jmin2=jmin
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin+kd
        iminc=iminb-kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (up) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
        imin2=imin+tmp
        jmin2=jmin
     end if


     ! triangular
     do i=1,rhs
        call DCOPY(tmp,B(iminb,i),1,helper(1,i),1)
     end do

     call DTRMM(SIDE,UPLO2,TRANSA,DIAG,tmp,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call DAXPY(tmp,DONE,helper(1,i),1,C(iminc,i),1)
     end do

     ! rectangular
     call DGEMM(TRANSA,TRANSB,kd-tmp,rhs,tmp, alpha, A(imin2,jmin2), LDA-1, B(iminb,1), LDB, DONE, C(iminc+tmp,1), LDC)

  end if


  call wdeallocate_2d(helper)  

end subroutine DSBMM






subroutine ZHBMM(UPLO,n,rhs,kd,alpha,A,LDA,B,LDB,beta,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose
  !  =======
  !  Perform the multiplication:
  !                              C=alpha*A*B+beta*C
  !  where alpha and beta are double complex scalars, x and y are double complex vectors and A is an
  !  n-by-n (structural symmetric) double complex band matrix, with kd sub-diagonals and kd super-diagonals,
  !  B and C are n-by-rhs matrices.
  !
  !  It is a  Blocked version (BLAS3) of the routine ZHBMV which allow blocking with 
  !  multiple right-hand sides  (with  generalization to structural symetric case.. not necessarely Hermitian).    
  !
  !
  !  Arguments
  !  =========
  !
  ! UPLO (input) character: specifies whether the upper or lower
  !                         triangular part of the band matrix A is being supplied as follows:
  !
  !             UPLO = 'U' or 'u'   The upper triangular part of A is
  !                                 being supplied. A is Hermitian.
  !
  !             UPLO = 'L' or 'l'   The lower triangular part of A is
  !                                 being supplied. A is Hermitian.
  !
  !             UPLO = 'F' or 'f'   All element of A are supplied. A is supposed structural symmetric. 
  !
  !
  !  N      (input) INTEGER
  !          The number of row and columns of the matrix A.  
  !          The number of rows of matrices B and C.   N >= 0.
  !  rhs    (input) INTEGER
  !          The number of columns of matrices B and C
  !  KD      (input) INTEGER
  !          The number of subdiagonals kl=kd and superdiagonals ku=kd within the band of A.  KD >= 0.
  ! alpha   (input) COMPLEX DOUBLE PRECISION
  !  A      (input) COMPLEX DOUBLE PRECISION  
  !          array, dimension (LDA,N); the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! B       (input) COMPLEX DOUBLE PRECISION - 2D array
  ! LDB     (input) INTEGER      
  !          The leading dimension of the array B.  LDB >= N.
  ! beta    (input) COMPLEX DOUBLE PRECISION
  ! C       (input/output) COMPLEX DOUBLE PRECISION - 2D array
  !         On exit, contains the solutions 
  ! LDC     (input) INTEGER      
  !          The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! Eric Polizzi 2009
  !======================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  CHARACTER(len=1) :: UPLO
  INTEGER :: N, RHS, KD,LDA,LDB,LDC
  COMPLEX(kind=kind(1.0d0)):: alpha,beta
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*) ::  A
  COMPLEX(kind=kind(1.0d0)),dimension(LDB,*) ::  B
  COMPLEX(kind=kind(1.0d0)),dimension(LDC,*) ::  C
!!!!!!!!!!!!!!

  COMPLEX(kind=kind(1.0d0)),dimension(:,:),pointer :: helper

  integer :: info_alloc
  logical :: low,up,full

  integer :: bl, i , k , tmp, imin,jmin,iminb,iminc,imin2,jmin2
  character(len=1) :: UPLO2,DIAG,TRANSA,TRANSB, SIDE
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO)



  low=.false.
  up=.false.
  full=.true.
  if ((UPLO=='L').or.(UPLO=='l')) then
     low=.true.
     full=.false.
  end if
  if ((UPLO=='U').or.(UPLO=='u')) then
     up=.true.
     full=.false.
  end if



  bl = n/kd
  tmp = n - bl*kd



  call wallocate_2z(helper,kd,rhs,info_alloc)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! diagonal blocks first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TRANSA='N'
  TRANSB='N'

  SIDE= 'L'


  if (low) then
     imin=1
  else !up or full
     imin=kd+1
  end if

  if (tmp/=0) then
!!!! last block of size tmp
     jmin=bl*kd+1
     if (full) then
        call ZGEMM(TRANSA,TRANSB,tmp,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(jmin,1), LDB, beta, C(jmin,1), LDC)
     else ! low or up
        call ZHEMM (SIDE, UPLO, tmp, rhs, ALPHA, A(imin,jmin), LDA-1, B(jmin,1), LDB, BETA,C(jmin,1), LDC)
     end if
  end if



  if (bl==0) return

!!!! all other blocks of size kd
  if (bl>0) then
     do k=1,bl
        jmin=(k-1)*kd+1
        if (full) then
           call ZGEMM(TRANSA,TRANSB,kd,rhs,kd, alpha, A(imin,jmin), LDA-1, B(jmin,1), LDB, beta, C(jmin,1), LDC)
        else ! low or up
           call ZHEMM (SIDE, UPLO, kd, rhs, ALPHA, A(imin,jmin), LDA-1, B(jmin,1), LDB, BETA,C(jmin,1), LDC)
        end if
     end do
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Lower Subdiagonal blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SIDE='L'
  DIAG='N'

  if (full) then
     imin=2*kd+1
     TRANSA='N'
     UPLO2='U'
  elseif (low) then
     imin=kd+1
     TRANSA='N'
     UPLO2='U'
  elseif (up) then
     imin=1
     TRANSA='C'
     UPLO2='L'
  end if

  do k=1,bl-1
     if (full) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
     elseif (up) then
        jmin=k*kd+1
        iminb=(k-1)*kd+1
        iminc=iminb+kd
     end if
     do i=1,rhs
        call ZCOPY(kd,B(iminb,i),1,helper(1,i),1)
     end do
     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,kd,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)
     do i=1,rhs
        call ZAXPY(kd,ONEC,helper(1,i),1,C(iminc,i),1)
     end do
  end do

!!!! final block (triangular+rectangular)
  if (tmp/=0) then
     k=bl
     if (full) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (up) then
        jmin=k*kd+1
        iminb=(k-1)*kd+1
        iminc=iminb+kd
        imin2=imin+tmp
        jmin2=jmin
     end if

     ! triangular
     do i=1,rhs
        call ZCOPY(tmp,B(iminb,i),1,helper(1,i),1)
     end do

     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,tmp,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call ZAXPY(tmp,ONEC,helper(1,i),1,C(iminc,i),1)
     end do

     ! rectangular
     call ZGEMM(TRANSA,TRANSB,tmp,rhs,kd-tmp, alpha, A(imin2,jmin2), LDA-1, B(iminb+tmp,1), LDB, ONEC, C(iminc,1), LDC)



  end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Upper Superdiagonal blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SIDE='L'
  DIAG='N'

  if (full) then
     imin=1
     TRANSA='N'
     UPLO2='L'
  elseif (low) then
     imin=kd+1
     TRANSA='C'
     UPLO2='U'
  elseif (up) then
     imin=1
     TRANSA='N'
     UPLO2='L'
  end if


  do k=1,bl-1
     if (full) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin+kd
        iminc=iminb-kd
     elseif (up) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
     end if
     do i=1,rhs
        call ZCOPY(kd,B(iminb,i),1,helper(1,i),1)
     end do

     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,kd,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call ZAXPY(kd,ONEC,helper(1,i),1,C(iminc,i),1)
     end do

  end do




!!!! final block (triangular+rectangular)
  if (tmp/=0) then
     k=bl
     if (full) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
        imin2=imin+tmp
        jmin2=jmin
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin+kd
        iminc=iminb-kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (up) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
        imin2=imin+tmp
        jmin2=jmin
     end if



     ! triangular
     do i=1,rhs
        call ZCOPY(tmp,B(iminb,i),1,helper(1,i),1)
     end do

     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,tmp,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call ZAXPY(tmp,ONEC,helper(1,i),1,C(iminc,i),1)
     end do



     ! rectangular
     call ZGEMM(TRANSA,TRANSB,kd-tmp,rhs,tmp, alpha, A(imin2,jmin2), LDA-1, B(iminb,1), LDB, ONEC, C(iminc+tmp,1), LDC)


  end if


  call wdeallocate_2z(helper)  

end subroutine ZHBMM






subroutine ZSBMM(UPLO,n,rhs,kd,alpha,A,LDA,B,LDB,beta,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose
  !  =======
  !  Perform the multiplication:
  !                              C=alpha*A*B+beta*C
  !  where alpha and beta are double complex scalars, x and y are double complex vectors and A is an
  !  n-by-n (structural symmetric) double complex band matrix, with kd sub-diagonals and kd super-diagonals,
  !  B and C are n-by-rhs matrices.
  !
  !  It is a  Blocked version (BLAS3) of the routine ZSBMV which allow blocking with 
  !  multiple right-hand sides  (with  generalization to structural symetric case).    
  !
  !
  !  Arguments
  !  =========
  !
  ! UPLO (input) character: specifies whether the upper or lower
  !                         triangular part of the band matrix A is being supplied as follows:
  !
  !             UPLO = 'U' or 'u'   The upper triangular part of A is
  !                                 being supplied. A is supposed symmetric.
  !
  !             UPLO = 'L' or 'l'   The lower triangular part of A is
  !                                 being supplied. A is supposed symmetric.
  !
  !             UPLO = 'F' or 'f'   All element of A are supplied. A is supposed structural symmetric.
  !

  !
  !  N      (input) INTEGER
  !          The number of row and columns of the matrix A.  
  !          The number of rows of matrices B and C.   N >= 0.
  !  rhs    (input) INTEGER
  !          The number of columns of matrices B and C
  !  KD      (input) INTEGER
  !          The number of subdiagonals kl=kd and superdiagonals ku=kd within the band of A.  KD >= 0.
  ! alpha   (input) COMPLEX DOUBLE PRECISION
  !  A      (input) COMPLEX DOUBLE PRECISION  
  !          array, dimension (LDA,N); the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! B       (input) COMPLEX DOUBLE PRECISION - 2D array
  ! LDB     (input) INTEGER      
  !          The leading dimension of the array B.  LDB >= N.
  ! beta    (input) COMPLEX DOUBLE PRECISION
  ! C       (input/output) COMPLEX DOUBLE PRECISION - 2D array
  !         On exit, contains the solutions 
  ! LDC     (input) INTEGER      
  !          The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! Eric Polizzi 2009
  !======================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  CHARACTER(len=1) :: UPLO
  INTEGER :: N, RHS, KD,LDA,LDB,LDC
  COMPLEX(kind=kind(1.0d0)):: alpha,beta
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*) ::  A
  COMPLEX(kind=kind(1.0d0)),dimension(LDB,*) ::  B
  COMPLEX(kind=kind(1.0d0)),dimension(LDC,*) ::  C
!!!!!!!!!!!!!!
  COMPLEX(kind=kind(1.0d0)),dimension(:,:),pointer :: helper
  integer :: info_alloc
  logical :: low,up,full
  integer :: bl, i , k , tmp, imin,jmin,iminb,iminc,imin2,jmin2
  character(len=1) :: UPLO2,DIAG,TRANSA,TRANSB, SIDE
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO)


  low=.false.
  up=.false.
  full=.true.
  if ((UPLO=='L').or.(UPLO=='l')) then
     low=.true.
     full=.false.
  end if
  if ((UPLO=='U').or.(UPLO=='u')) then
     up=.true.
     full=.false.
  end if


  bl = n/kd
  tmp = n - bl*kd

  call wallocate_2z(helper,kd,rhs,info_alloc)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! diagonal blocks first
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TRANSA='N'
  TRANSB='N'

  SIDE= 'L'


  if (low) then
     imin=1
  else !up or full
     imin=kd+1
  end if

  if (tmp/=0) then
!!!! last block of size tmp
     jmin=bl*kd+1
     if (full) then
        call ZGEMM(TRANSA,TRANSB,tmp,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(jmin,1), LDB, beta, C(jmin,1), LDC)
     else ! low or up
        call ZSYMM (SIDE, UPLO, tmp, rhs, ALPHA, A(imin,jmin), LDA-1, B(jmin,1), LDB, BETA,C(jmin,1), LDC)
     end if
  end if


  if (bl==0) return

!!!! all other blocks of size kd
  if (bl>0) then
     do k=1,bl
        jmin=(k-1)*kd+1
        if (full) then
           call ZGEMM(TRANSA,TRANSB,kd,rhs,kd, alpha, A(imin,jmin), LDA-1, B(jmin,1), LDB, beta, C(jmin,1), LDC)
        else ! low or up
           call ZSYMM (SIDE, UPLO, kd, rhs, ALPHA, A(imin,jmin), LDA-1, B(jmin,1), LDB, BETA,C(jmin,1), LDC)
        end if
     end do
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Lower Subdiagonal blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SIDE='L'
  DIAG='N'

  if (full) then
     imin=2*kd+1
     TRANSA='N'
     UPLO2='U'
  elseif (low) then
     imin=kd+1
     TRANSA='N'
     UPLO2='U'
  elseif (up) then
     imin=1
     TRANSA='T'
     UPLO2='L'
  end if

  do k=1,bl-1
     if (full) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
     elseif (up) then
        jmin=k*kd+1
        iminb=(k-1)*kd+1
        iminc=iminb+kd
     end if
     do i=1,rhs
        call ZCOPY(kd,B(iminb,i),1,helper(1,i),1)
     end do
     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,kd,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)
     do i=1,rhs
        call ZAXPY(kd,ONEC,helper(1,i),1,C(iminc,i),1)
     end do
  end do

!!!! final block (triangular+rectangular)
  if (tmp/=0) then
     k=bl
     if (full) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin
        iminc=iminb+kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (up) then
        jmin=k*kd+1
        iminb=(k-1)*kd+1
        iminc=iminb+kd
        imin2=imin+tmp
        jmin2=jmin
     end if

     ! triangular
     do i=1,rhs
        call ZCOPY(tmp,B(iminb,i),1,helper(1,i),1)
     end do

     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,tmp,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call ZAXPY(tmp,ONEC,helper(1,i),1,C(iminc,i),1)
     end do

     ! rectangular
     call ZGEMM(TRANSA,TRANSB,tmp,rhs,kd-tmp, alpha, A(imin2,jmin2), LDA-1, B(iminb+tmp,1), LDB, ONEC, C(iminc,1), LDC)

  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Upper Superdiagonal blocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SIDE='L'
  DIAG='N'

  if (full) then
     imin=1
     TRANSA='N'
     UPLO2='L'
  elseif (low) then
     imin=kd+1
     TRANSA='T'
     UPLO2='U'
  elseif (up) then
     imin=1
     TRANSA='N'
     UPLO2='L'
  end if


  do k=1,bl-1
     if (full) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin+kd
        iminc=iminb-kd
     elseif (up) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
     end if
     do i=1,rhs
        call ZCOPY(kd,B(iminb,i),1,helper(1,i),1)
     end do

     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,kd,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call ZAXPY(kd,ONEC,helper(1,i),1,C(iminc,i),1)
     end do

  end do


!!!! final block (triangular+rectangular)
  if (tmp/=0) then
     k=bl
     if (full) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
        imin2=imin+tmp
        jmin2=jmin
     elseif (low) then
        jmin=(k-1)*kd+1
        iminb=jmin+kd
        iminc=iminb-kd
        imin2=imin-tmp
        jmin2=jmin+tmp
     elseif (up) then
        jmin=k*kd+1
        iminb=jmin
        iminc=iminb-kd
        imin2=imin+tmp
        jmin2=jmin
     end if

     ! triangular
     do i=1,rhs
        call ZCOPY(tmp,B(iminb,i),1,helper(1,i),1)
     end do

     call ZTRMM(SIDE,UPLO2,TRANSA,DIAG,tmp,rhs,alpha,A(imin,jmin),LDA-1,helper,kd)

     do i=1,rhs
        call ZAXPY(tmp,ONEC,helper(1,i),1,C(iminc,i),1)
     end do

     ! rectangular
     call ZGEMM(TRANSA,TRANSB,kd-tmp,rhs,tmp, alpha, A(imin2,jmin2), LDA-1, B(iminb,1), LDB, ONEC, C(iminc+tmp,1), LDC)

  end if


  call wdeallocate_2z(helper)  

end subroutine ZSBMM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ZGBMM(transa,transb,n,rhs,kl,ku,alpha,A,LDA,B,LDB,beta,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose
  !  =======
  !  Perform the multiplication:
  !                              C=alpha*op(A)*op(B)+beta*C
  !  where alpha and beta are double complex scalars, x and y are double complex vectors and A is an
  !  n-by-n (non-symmetric) double complex band matrix, with kl sub-diagonals and ku super-diagonals,
  !  B and C are n-by-rhs matrices.
  !
  !
  !  Arguments
  !  =========
  !  TRANSA (input) CHARACTER
  !          'N'    A is normal
  !          'T'    A is transposed
  !          'C'    A is complex conjugate tranposed
  !  TRANSB (input) CHARACTER
  !          'N'    B is normal
  !          'T'    B is transposed
  !          'C'    B is complex conjugate tranposed
  !  N      (input) INTEGER
  !          The number of row and columns of the matrix A.  
  !          The number of rows of matrices B and C.   N >= 0.
  !  rhs    (input) INTEGER
  !          The number of columns of matrices B and C
  !  KU      (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  KL      (input) INTEGER 
  !                 bandwidth = kl+ku+1
  !          The number of subdiagonals within the band of A.  KL >= 0.
  ! alpha   (input) COMPLEX DOUBLE PRECISION
  !  A      (input) COMPLEX DOUBLE PRECISION  
  !          array, dimension (LDA,N); the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! B       (input) COMPLEX DOUBLE PRECISION - 2D array
  ! LDB     (input) INTEGER      
  !          The leading dimension of the array B.  LDB >= N.
  ! beta    (input) COMPLEX DOUBLE PRECISION
  ! C       (input/output) COMPLEX DOUBLE PRECISION - 2D array
  !         On exit, contains the solutions 
  ! LDC     (input) INTEGER      
  !          The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! James Kestyn 2013
  !======================================================================
  implicit none
  CHARACTER :: TRANSA, TRANSB
  INTEGER :: N, RHS, KU,KL,LDA,LDB,LDC
  COMPLEX(kind=kind(1.0d0)):: alpha,beta
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*) ::  A
  COMPLEX(kind=kind(1.0d0)),dimension(LDB,*) ::  B
  COMPLEX(kind=kind(1.0d0)),dimension(LDC,*) ::  C
!!!!!!!!!!!!!!

  integer :: info_alloc
  logical :: low,up,full

  integer :: bl, i,j, k , kd, tmp, imin,jmin,iminb,iminc,imin2,jmin2
  character(len=1) :: UPLO2,DIAG,SIDE
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO)


  if( TRANSA == 'N' .or. TRANSA == 'n' ) then

     do j=1,kl !first kl rows of A
        tmp = (ku)+j
        imin= (ku)+j
        jmin= 1
        call ZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(1,1), LDB, beta, C(j,1), LDC)
     enddo

     do j=1,N-ku-kl !Middle rows of A (no offset)
        tmp = ku+kl+1
        imin = ku+kl+1
        call ZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,j), LDA-1, B(j,1), LDB, beta, C(kl+j,1), LDC)
     enddo

     do j=1,ku !last ku rows of A
        tmp = (ku+kl+1)-j
        imin= kl+ku+1
        jmin= (N-kl-ku)+j!(ku+2)-i
        call ZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(N-kl-ku+j,1), LDB, beta, C(N-ku+j,1), LDC)
     enddo


  else

     do j=1,ku !first ku rows of A^C
        tmp = (kl)+j
        imin= (ku+2)-j
        jmin= j
        call ZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(1,1), LDB, beta, C(jmin,1), LDC)
     enddo


     do j=1,N-ku-kl
        tmp = ku+kl+1
        imin= 1
        jmin= ku+j
        call ZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(j,1), LDB, beta, C(jmin,1), LDC)
     enddo

     do j=1,kl !last kl rows of A^C
        tmp = (ku+kl+1)-j
        imin= 1
        jmin= (N-kl)+j!(ku+2)-i
        call ZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(N-ku-kl+j,1), LDB, beta, C(jmin,1), LDC)
     enddo

  endif

end subroutine ZGBMM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DGBMM(transa,transb,n,rhs,kl,ku,alpha,A,LDA,B,LDB,beta,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose
  !  =======
  !  Perform the multiplication:
  !                              C=alpha*op(A)*op(B)+beta*C
  !  where alpha and beta are double scalars, x and y are double vectors and A is an
  !  n-by-n (non-symmetric) double band matrix, with kl sub-diagonals and ku super-diagonals,
  !  B and C are n-by-rhs matrices.
  !
  !
  !  Arguments
  !  =========
  !  TRANSA (input) CHARACTER
  !          'N'    A is normal
  !          'T'    A is transposed
  !  TRANSB (input) CHARACTER
  !          'N'    B is normal
  !          'T'    B is transposed
  !  N      (input) INTEGER
  !          The number of row and columns of the matrix A.  
  !          The number of rows of matrices B and C.   N >= 0.
  !  rhs    (input) INTEGER
  !          The number of columns of matrices B and C
  !  KU      (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  KL      (input) INTEGER 
  !                 bandwidth = kl+ku+1
  !          The number of subdiagonals within the band of A.  KL >= 0.
  ! alpha   (input) COMPLEX DOUBLE PRECISION
  !  A      (input) COMPLEX DOUBLE PRECISION  
  !          array, dimension (LDA,N); the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! B       (input) COMPLEX DOUBLE PRECISION - 2D array
  ! LDB     (input) INTEGER      
  !          The leading dimension of the array B.  LDB >= N.
  ! beta    (input) COMPLEX DOUBLE PRECISION
  ! C       (input/output) COMPLEX DOUBLE PRECISION - 2D array
  !         On exit, contains the solutions 
  ! LDC     (input) INTEGER      
  !          The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! James Kestyn 2013
  ! Braegan Spring 2015 -- converted from complex to double
  !======================================================================
  implicit none
  CHARACTER :: TRANSA, TRANSB
  INTEGER :: N, RHS, KU,KL,LDA,LDB,LDC
  Double Precision:: alpha,beta
  Double Precision,dimension(LDA,*) ::  A
  Double Precision,dimension(LDB,*) ::  B
  Double Precision,dimension(LDC,*) ::  C
!!!!!!!!!!!!!!

  integer :: info_alloc
  logical :: low,up,full

  integer :: bl, i,j, k , kd, tmp, imin,jmin,iminb,iminc,imin2,jmin2
  character(len=1) :: UPLO2,DIAG,SIDE
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0


  if( TRANSA == 'N' .or. TRANSA == 'n' ) then

     do j=1,kl !first kl rows of A
        tmp = (ku)+j
        imin= (ku)+j
        jmin= 1
        call DGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(1,1), LDB, beta, C(j,1), LDC)
     enddo

     do j=1,N-ku-kl !Middle rows of A (no offset)
        tmp = ku+kl+1
        imin = ku+kl+1
        call DGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,j), LDA-1, B(j,1), LDB, beta, C(kl+j,1), LDC)
     enddo

     do j=1,ku !last ku rows of A
        tmp = (ku+kl+1)-j
        imin= kl+ku+1
        jmin= (N-kl-ku)+j!(ku+2)-i
        call DGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(N-kl-ku+j,1), LDB, beta, C(N-ku+j,1), LDC)
     enddo


  else

     do j=1,ku !first ku rows of A^C
        tmp = (kl)+j
        imin= (ku+2)-j
        jmin= j
        call DGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(1,1), LDB, beta, C(jmin,1), LDC)
     enddo


     do j=1,N-ku-kl
        tmp = ku+kl+1
        imin= 1
        jmin= ku+j
        call DGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(j,1), LDB, beta, C(jmin,1), LDC)
     enddo

     do j=1,kl !last kl rows of A^C
        tmp = (ku+kl+1)-j
        imin= 1
        jmin= (N-kl)+j!(ku+2)-i
        call DGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(N-ku-kl+j,1), LDB, beta, C(jmin,1), LDC)
     enddo

  endif

end subroutine DGBMM




subroutine DZGBMM(TRANSA,TRANSB,n,rhs,kl,ku,alpha,A,LDA,B,LDB,beta,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Purpose
  !  =======
  !  Perform the multiplication:
  !                              C=alpha*op(A)*op(B)+beta*C
  !  where alpha and beta are double complex scalars, x and y are double complex vectors and A is an
  !  n-by-n (non-symmetric) double precision band matrix, with kl sub-diagonals and ku super-diagonals,
  !  B and C are complex double precision n-by-rhs matrices.
  !
  !  It is a  Blocked version (BLAS3) of the routine ZHBMV which allow blocking with 
  !  multiple right-hand sides  (with  generalization to structural symetric case.. not necessarely Hermitian).    
  !
  !
  !  Arguments
  !  =========
  !  TRANSA (input) CHARACTER
  !          'N'    A is normal
  !          'T'    A is transposed
  !          'C'    A is complex conjugate tranposed
  !  TRANSB (input) CHARACTER
  !          'N'    B is normal
  !          'T'    B is transposed
  !          'C'    B is complex conjugate tranposed
  !  N      (input) INTEGER
  !          The number of row and columns of the matrix A.  
  !          The number of rows of matrices B and C.   N >= 0.
  !  rhs    (input) INTEGER
  !          The number of columns of matrices B and C
  !  KU      (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !  KL      (input) INTEGER 
  !                 bandwidth = kl+ku+1
  !          The number of subdiagonals within the band of A.  KL >= 0.
  ! alpha   (input) COMPLEX DOUBLE PRECISION
  !  A      (input) DOUBLE PRECISION  
  !          array, dimension (LDA,N); the matrix A in band storage, in rows 1 to KL+KU+1;       
  !          The j-th column of A is stored in the j-th column of the
  !          array A as follows:
  !          A(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  ! LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= KL+KU+1.
  ! B       (input) COMPLEX DOUBLE PRECISION - 2D array
  ! LDB     (input) INTEGER      
  !          The leading dimension of the array B.  LDB >= N.
  ! beta    (input) COMPLEX DOUBLE PRECISION
  ! C       (input/output) COMPLEX DOUBLE PRECISION - 2D array
  !         On exit, contains the solutions 
  ! LDC     (input) INTEGER      
  !          The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! James Kestyn 2015
  !======================================================================
  implicit none
  CHARACTER :: TRANSA, TRANSB
  INTEGER :: N, RHS, KU,KL,LDA,LDB,LDC
  COMPLEX(kind=kind(1.0d0)):: alpha,beta
  DOUBLE PRECISION,dimension(LDA,*) ::  A
  COMPLEX(kind=kind(1.0d0)),dimension(LDB,*) ::  B
  COMPLEX(kind=kind(1.0d0)),dimension(LDC,*) ::  C
!!!!!!!!!!!!!!

  integer :: info_alloc
  logical :: low,up,full

  integer :: bl, i,j, k , kd, tmp, imin,jmin,iminb,iminc,imin2,jmin2
  character(len=1) :: UPLO2,DIAG,SIDE
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO)

  if( TRANSA == 'N' .or. TRANSA == 'n') then

     do j=1,kl !first kl rows of A
        tmp = (ku)+j
        imin= (ku)+j
        jmin= 1
        !call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(1,1), LDB, beta, C(j,1), LDC)
     enddo

     do j=1,N-ku-kl !Middle rows of A (no offset)
        tmp = ku+kl+1
        imin = ku+kl+1
        !call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,j), LDA-1, B(j,1), LDB, beta, C(kl+j,1), LDC)
     enddo

     do j=1,ku !last ku rows of A
        tmp = (ku+kl+1)-j
        imin= kl+ku+1
        jmin= (N-kl-ku)+j!(ku+2)-i
        !call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA-1, B(N-kl-ku+j,1), LDB, beta, C(N-ku+j,1), LDC)
     enddo

  else

     do j=1,ku !first ku rows of A^C
        tmp = (kl)+j
        imin= (ku+2)-j
        jmin= j
        !call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(1,1), LDB, beta, C(jmin,1), LDC)
     enddo


     do j=1,N-ku-kl
        tmp = ku+kl+1
        imin= 1
        jmin= ku+j
        !call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(j,1), LDB, beta, C(jmin,1), LDC)
     enddo

     do j=1,kl !last kl rows of A^C
        tmp = (ku+kl+1)-j
        imin= 1
        jmin= (N-kl)+j!(ku+2)-i
        !call DZGEMM(TRANSA,TRANSB,1,rhs,tmp, alpha, A(imin,jmin), LDA, B(N-ku-kl+j,1), LDB, beta, C(jmin,1), LDC)
     enddo

  endif

end subroutine DZGBMM





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine DTBSMPY(UPLO,TRANS,DIAG,N,rhs,kd,A,LDA,B,LDB,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Minimally modified version of DTBSM to include sum along with blocked vectors
  !Purpose:solves one of the systems of equations in real double precision
  !        Additionally, adds a vector to the solved vector. 
  !     A*x = b,   or   A'*x = b,  
  !     and x=x+c ; result in c
  !  where b, x, and c are n by rhs element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular band matrix, with ( kd + 1)
  !  diagonals.
  !
  !  It is a  Blocked version (BLAS3) of the routine DTBSV which allow blocking with 
  !  multiple right-hand sides      
  !
  ! Arguments: 
  !       UPLO (input) character: if 'L' or 'l' lower triangular , if 'U' or 'u' upper triangular
  !              Note : I find this argument confusing in the transpose 
  !              case -- if we have A=LU, A^T=U^T L^T, U^T is lower triangular, so
  !              the U is marked 'L.'
  !              
  !       TRANS (input) character:specifies the equations to be solved as follows:
  !              TRANS = 'N' or 'n'   A*x = b.
  !              TRANS = 'T' or 't'   A'*x = b.
  !       DIAG (input) character: if 'N' or 'n' non-unit diagonal, if 'U' or 'u' unit diagonal
  !       N (input) integer: order of matrix A
  !       rhs (input) integer :number of right-hand-sides for vector B        
  !       kd (input) integer: lower bandwidth of the matrix A if UPLO='L', or
  !                           upper bandwidth of the matrix A if UPLO='U'
  !       A (input) double precision array, banded input matrix 
  !           Note: same format than the one used in BLAS DTBSV
  !
  !           Before entry with UPLO = 'U' or 'u', the leading ( kd + 1 )
  !           by n part of the array A must contain the upper triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row
  !           ( kd + 1 ) of the array, the first super-diagonal starting at
  !           position 2 in row kd, and so on. The top left kd by kd triangle
  !           of the array A is not referenced.
  !
  !           Before entry with UPLO = 'L' or 'l', the leading ( kd + 1 )
  !           by n part of the array A must contain the lower triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row 1 of
  !           the array, the first sub-diagonal starting at position 1 in
  !           row 2, and so on. The bottom right kd by kd triangle of the
  !           array A is not referenced.
  !           
  !       LDA (input) integer : The leading dimension of the array A.  LDA >= kd+1.
  !       B (input/output) double precision 2D array: right hand side contains the solution on exit 
  !       LDB (input) integer : The leading dimension of the array B.  LDB >= N.
  !       C (input) : double precision 2D array:
  !       LDC (input) integer : The leading dimension of the array C.  LDB >= N.
  !======================================================================
  ! Eric Polizzi 2009
  ! Modified by Braegan Spring 2014
  !======================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  integer :: N,LDA,LDB,LDC,rhs
  DOUBLE PRECISION,dimension(LDA,*),intent(in) ::   A
  character(len=1)  :: UPLO,DIAG,TRANS
  integer :: kd
  DOUBLE PRECISION,dimension(LDB,*) ::   B
  DOUBLE PRECISION,dimension(LDC,*) ::   C

!!!!!!!!!!!!!!!!!!!!!!!!
  double precision,parameter :: DONE=1.0d0
  character(len=1) :: TRANSB,SIDE,UPLO2,DIAG2
  integer :: k,ijinit,i
  integer :: m,m1,max
  integer :: i1,i2,i3,i4,i5,i6,j1,j2,j3,j4,s1,s2
  !double precision,dimension(kd,rhs) :: aux
  double precision,pointer,dimension(:,:) :: aux
  integer :: ldaux
  integer :: infoloc
!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_2d(aux,kd,rhs,infoloc)


  !! Size of the Recursive diagonal block and initialization of UPLO2
  !! for triangular off-diagonal blocks
  m=min(kd,n)


!!!!!!!!!!!!!! Case n<kd
if (n<kd) then ! treat as dense
  if ((UPLO=='L').or.(UPLO=='l')) then
call DTRSM('L',UPLO,TRANS,DIAG,N,rhs,DONE,A(1,1),LDA-1,B,LDB)
  elseif ((UPLO=='U').or.(UPLO=='u')) then
call DTRSM('L',UPLO,TRANS,DIAG,N,rhs,DONE,A(kd+1,1),LDA-1,B,LDB)
endif
return
endif
!!!!!!!!!!!!


  if ((UPLO=='L').or.(UPLO=='l')) then
     UPLO2='U'
  elseif ((UPLO=='U').or.(UPLO=='u')) then
     UPLO2='L'
  end if

  !!initialization of DIAG2,SIDE
  DIAG2='N'
  SIDE='L'


  !! how many submatrices (m*m) in A
  max=n/kd
  !! size of the first diagonal block before recursion
  m1=n-max*kd

!print *,n,m1,max
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Initialisation before the loop
!!!! only once for the first block decomposition with size m1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (m1/=0) then !! if there is a first block

     ldaux=kd

     TRANSB = 'N'

     if ((TRANS=='N').or.(TRANS=='n')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=1
           i1=1
           j3=j1
           i3=i1+m1
           j4=j1
           i4=kd+1

           i5=j1+m1 !rhs
           i6=i5+(kd-m1)

        elseif((UPLO=='U').or.(UPLO=='u')) then
           j1=n-m1+1                
           i1=kd+1
           j3=j1
           i3=m1+1
           j4=j1
           i4=1

           i5=j1-(kd-m1)
           i6=i5-m1
        end if

     elseif ((TRANS=='T').or.(TRANS=='t')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=n-m1+1
           i1=1
           j3=n-kd+1
           i3=kd+1-m1  
           j4=j1-kd
           i4=kd+1

           i5=j1-(kd-m1)
           i6=i5-m1
        elseif ((UPLO=='U').or.(UPLO=='u')) then
           j1=1
           i1=kd+1
           j3=j1+m1
           i3=i1-m1
           j4=kd+1
           i4=1

           i5=j1+m1
           i6=i5+(kd-m1)
        end if
     end if

!!!!!! Solve using blas3 the initial block
     call DTRSM ( SIDE, UPLO, TRANS, DIAG, m1, rhs, DONE, A(i1,j1), LDA-1,B(j1,1), LDB )
!!!!!!!!!!! update c<=c+b
do i=1,rhs
call DAXPY(m1,DONE,B(j1,i),1,C(j1,i),1) 
enddo
!!!!!! update the right hand side using blas3 in 2 steps
     ! step 1
     call DGEMM ( TRANS, TRANSB, kd-m1, rhs, m1, -DONE, A(i3,j3), LDA-1, B(j1,1), LDB,DONE, B(i5,1), LDB ) 
     ! step 2
     do i=1,rhs
        call DCOPY(m1,B(j1,i),1,aux(1,i),1)
     end do
     call DTRMM ( SIDE, UPLO2, TRANS, DIAG2, m1, rhs, DONE, A(i4,j4), LDA-1,aux,ldaux ) 
     do i=1,rhs
        call DAXPY(m1,-DONE,aux(1,i),1,B(i6,i),1)
     end do

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! starting point origin !! new origin and other initialization before loop
  if ((TRANS=='N').or.(TRANS=='n')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=m1+1
        i1=1
        i2=kd+1
        s1=1
        s2=0
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=n-m1-kd+1
        i1=kd+1
        i2=1
        s1=-1
        s2=0
     end if
  elseif ((TRANS=='T').or.(TRANS=='t')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=n-m1-kd+1
        i1=1
        i2=kd+1
        s1=-1
        s2=-1
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=m1+1
        i1=kd+1
        i2=1
        s1=1
        s2=1
     end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! LOOP for all the submatrices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO k=1,max
     j1=ijinit+s1*(k-1)*m
     j2=j1+s2*kd
     i3=j1+s1*kd
     !! Solve using blas3 the initial bloc
     call DTRSM ( SIDE, UPLO, TRANS, DIAG, kd, rhs, DONE, A(i1,j1),LDA-1,B(j1,1), LDB ) 
!! update c <= c+b
do i=1,rhs
     call DAXPY(kd,DONE,B(j1,i),1,C(j1,i),1)
enddo
     !! update the right hand side using blas3
     if (k/=max) then
        do i=1,rhs
           call DCOPY(kd,B(j1,i),1,aux(1,i),1)
        end do
        call DTRMM ( SIDE, UPLO2, TRANS, DIAG2, kd, rhs, DONE, A(i2,j2), LDA-1,aux,kd ) 
        do i=1,rhs
          call DAXPY(kd,-DONE,aux(1,i),1,B(i3,i),1)
        end do
     end if
  END DO
!  do i=1,rhs
!    call DAXPY(m1,DONE,B(n-m1+1,i),1,C(n-m1+1,i),1)
!  enddo
!  do i=1,rhs
!    call DAXPY(kd,DONE,B(1,i),1,C(1,i),1)
!  enddo

  call wdeallocate_2d(aux)
end subroutine DTBSMPY



subroutine ZTBSMPY(UPLO,TRANS,DIAG,N,rhs,kd,A,LDA,B,LDB,C,LDC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Minimally modified version of ZTBSM to include sum along with blocked vectors
  !Purpose:solves one of the systems of equations in complex double precision
  !
  !     A*x = b,   or   A'*x = b,  or   conjg( A' )*x = b,
  !     and x=x+c ; result in c
  !  where b, x and c are n by rhs element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular band matrix, with ( kd + 1)
  !  diagonals.
  !
  !  It is a  Blocked version (BLAS3) of the routine DTBSV which allow blocking with 
  !  multiple right-hand sides      
  !
  ! Arguments: 
  !       UPLO (input) character: if 'L' or 'l' lower triangular , if 'U' or 'u' upper triangular
  !       TRANS (input) character:specifies the equations to be solved as follows:
  !              TRANS = 'N' or 'n'   A*x = b.
  !              TRANS = 'T' or 't'   A'*x = b.
  !              TRANS = 'C' or 'c'   conjg( A' )*x = b.
  !       DIAG (input) character: if 'N' or 'n' non-unit diagonal, if 'U' or 'u' unit diagonal
  !       N (input) integer: order of matrix A
  !       rhs (input) integer :number of right-hand-sides for vector B        
  !       kd (input) integer: lower bandwidth of the matrix A if UPLO='L', or
  !                           upper bandwidth of the matrix A if UPLO='U'
  !       A (input) complex double precision array, banded input matrix 
  !           Note: same format than the one used in BLAS DTBSV
  !
  !           Before entry with UPLO = 'U' or 'u', the leading ( kd + 1 )
  !           by n part of the array A must contain the upper triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row
  !           ( kd + 1 ) of the array, the first super-diagonal starting at
  !           position 2 in row kd, and so on. The top left kd by kd triangle
  !           of the array A is not referenced.
  !
  !           Before entry with UPLO = 'L' or 'l', the leading ( kd + 1 )
  !           by n part of the array A must contain the lower triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row 1 of
  !           the array, the first sub-diagonal starting at position 1 in
  !           row 2, and so on. The bottom right kd by kd triangle of the
  !           array A is not referenced.
  !           
  !       LDA (input) integer : The leading dimension of the array A.  LDA >= kd+1.
  !       B (input) complex double precision 2D array.
  !       LDB (input) integer : The leading dimension of the array B.  LDB >= N.
  !       C (input/output) complex double precision 2D array; contains the solution on exit
  !       LDC (input) integer : The leading dimension of the array C.  LDC >= N.
  !
  !======================================================================
  ! Eric Polizzi 2009
  ! Modified by Braegan Spring 2014
  !======================================================================
  implicit none
    include 'f90_noruntime_interface.fi'
  integer :: N,LDA,LDB,LDC,rhs
  COMPLEX(kind=kind(1.0d0)),dimension(LDA,*),intent(in) ::   A
  character(len=1)  :: UPLO,DIAG,TRANS
  integer :: kd
  COMPLEX(kind=kind(1.0d0)),dimension(LDB,*) ::   B
  COMPLEX(kind=kind(1.0d0)),dimension(LDC,*) ::   C

!!!!!!!!!!!!!!!!!!!!!!!!
  Double Precision, Parameter :: DONE=1.0d0, DZERO=0.0d0
  COMPLEX(kind=kind(1.0d0)),parameter :: ONEC=(DONE,DZERO)
  character(len=1) :: TRANSB,SIDE,UPLO2,DIAG2
  integer :: k,ijinit,i
  integer :: m,m1,max
  integer :: i1,i2,i3,i4,i5,i6,j1,j2,j3,j4,s1,s2
  !COMPLEX(kind=kind(1.0d0)),dimension(kd,rhs) :: aux
  COMPLEX(kind=kind(1.0d0)),pointer,dimension(:,:) :: aux
  integer :: ldaux
  integer :: infoloc
!!!!!!!!!!!!!!!!!!!!!!!!
  call wallocate_2z(aux,kd,rhs,infoloc) 

  !! Size of the Recursive diagonal block and initialization of UPLO2
  !! for triangular off-diagonal blocks
  m=min(kd,n)



!!!!!!!!!!!!!! Case n<kd
if (n<kd) then ! treat as dense
  if ((UPLO=='L').or.(UPLO=='l')) then
call ZTRSM('L',UPLO,TRANS,DIAG,N,rhs,ONEC,A(1,1),LDA-1,B,LDB)
  elseif ((UPLO=='U').or.(UPLO=='u')) then
call ZTRSM('L',UPLO,TRANS,DIAG,N,rhs,ONEC,A(kd+1,1),LDA-1,B,LDB)
endif
return
endif
!!!!!!!!!!!!!!!!!!


  if ((UPLO=='L').or.(UPLO=='l')) then
     UPLO2='U'
  elseif ((UPLO=='U').or.(UPLO=='u')) then
     UPLO2='L'
  end if

  !!initialization of DIAG2,SIDE
  DIAG2='N'
  SIDE='L'

  !! how many submatrices (m*m) in A
  max=n/kd
  !! size of the first diagonal block before recursion
  m1=n-max*kd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Initialisation before the loop
!!!! only once for the first block decomposition with size m1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (m1/=0) then !! if there is a first block

     ldaux=kd

     TRANSB = 'N'

     if ((TRANS=='N').or.(TRANS=='n')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=1
           i1=1
           j3=j1
           i3=i1+m1
           j4=j1
           i4=kd+1

           i5=j1+m1 !rhs
           i6=i5+(kd-m1)

        elseif((UPLO=='U').or.(UPLO=='u')) then
           j1=n-m1+1                
           i1=kd+1
           j3=j1
           i3=m1+1
           j4=j1
           i4=1

           i5=j1-(kd-m1)
           i6=i5-m1
        end if

     elseif ((TRANS=='T').or.(TRANS=='t').or.(TRANS=='C').or.(TRANS=='c')) then
        if ((UPLO=='L').or.(UPLO=='l')) then
           j1=n-m1+1
           i1=1
           j3=n-kd+1
           i3=kd+1-m1  
           j4=j1-kd
           i4=kd+1

           i5=j1-(kd-m1)
           i6=i5-m1
        elseif ((UPLO=='U').or.(UPLO=='u')) then
           j1=1
           i1=kd+1
           j3=j1+m1
           i3=i1-m1
           j4=kd+1
           i4=1

           i5=j1+m1
           i6=i5+(kd-m1)
        end if
     end if

!!!!!! Solve using blas3 the initial block
     call ZTRSM ( SIDE, UPLO, TRANS, DIAG, m1, rhs, ONEC, A(i1,j1), LDA-1,B(j1,1), LDB )
!!!!!!!!!!! update c<=c+b
do i=1,rhs
call ZAXPY(m1,ONEC,B(j1,i),1,C(j1,i),1) 
enddo
!!!!!! update the right hand side using blas3 in 2 steps
     ! step 1
     call ZGEMM ( TRANS, TRANSB, kd-m1, rhs, m1, -ONEC, A(i3,j3), LDA-1, B(j1,1), LDB,ONEC, B(i5,1), LDB ) 
     ! step 2
     do i=1,rhs
        call ZCOPY(m1,B(j1,i),1,aux(1,i),1)
     end do
     ! call ZLACPY( 'F', m1, rhs, B(j1,1), N, aux(1,1), ldaux )
     call ZTRMM ( SIDE, UPLO2, TRANS, DIAG2, m1, rhs, ONEC, A(i4,j4), LDA-1,aux,ldaux ) 
     do i=1,rhs
        call ZAXPY(m1,-ONEC,aux(1,i),1,B(i6,i),1)
     end do

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! starting point origin !! new origin and other initialization before loop
  if ((TRANS=='N').or.(TRANS=='n')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=m1+1
        i1=1
        i2=kd+1
        s1=1
        s2=0
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=n-m1-kd+1
        i1=kd+1
        i2=1
        s1=-1
        s2=0
     end if
  elseif ((TRANS=='T').or.(TRANS=='t').or.(TRANS=='C').or.(TRANS=='c')) then
     if ((UPLO=='L').or.(UPLO=='l')) then
        ijinit=n-m1-kd+1
        i1=1
        i2=kd+1
        s1=-1
        s2=-1
     elseif ((UPLO=='U').or.(UPLO=='u')) then
        ijinit=m1+1
        i1=kd+1
        i2=1
        s1=1
        s2=1
     end if
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! LOOP for all the submatrices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  DO k=1,max
     j1=ijinit+s1*(k-1)*m
     j2=j1+s2*kd
     i3=j1+s1*kd
     !! Solve using blas3 the initial bloc
     call ZTRSM ( SIDE, UPLO, TRANS, DIAG, kd, rhs, ONEC, A(i1,j1),LDA-1,B(j1,1), LDB )
!! update c <= c+b
do i=1,rhs
     call ZAXPY(kd,ONEC,B(j1,i),1,C(j1,i),1)
enddo 
     !! update the right hand side using blas3
     if (k/=max) then
        do i=1,rhs
           call ZCOPY(kd,B(j1,i),1,aux(1,i),1)
        end do
        call ZTRMM ( SIDE, UPLO2, TRANS, DIAG2, kd, rhs, -ONEC, A(i2,j2), LDA-1,aux,kd ) 
        do i=1,rhs
           call ZAXPY(kd,ONEC,aux(1,i),1,B(i3,i),1)
        end do
     end if
  END DO
  call wdeallocate_2z(aux) 

end subroutine ZTBSMPY



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


