c-----------------------------------------------------------------------
c     some routines extracted/ modified from SPARSKIT2 + one from blas
c-----------------------------------------------------------------------
      subroutine readmtc (nmax,nzmax,job,fname,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c     this  subroutine reads  a boeing/harwell matrix,  given the
c     corresponding file. handles right hand sides in full format
c     only (no sparse right hand sides). Also the matrix must  be
c     in assembled forms.
c     It differs from readmt, in that  the name of the file needs
c     to be passed, and then the file is opened and closed within
c     this routine.
c     Author: Youcef Saad - Date: Oct 31, 1989
c     updated Jul 20, 1998 by Irene Moulitsas
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax   =  max column dimension  allowed for matrix. The array ia should
c           be of length at least ncol+1 (see below) if job.gt.0
c nzmax  = max number of nonzeros elements allowed. the arrays a,
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c
c job    = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz, title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c                     rhs may contain initial guesses and exact
c                     solutions appended to the actual right hand sides.
c                     this will be indicated by the output parameter
c                     guesol [see below].
c
c fname = name of the file where to read the matrix from.
c
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c on return:
c----------
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs
c          is provided then it is rest to job=2 or job=1 depending on
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a      = the a matrix in the a, ia, ja (column) storage format
c ja     = column number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c          if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are appended to the right hand-sides in the array rhs.
c          if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are  appended to the right hand-sides
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol   = number of columns in matrix
c nnz    = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed.
c
c title  = character*72 = title of matrix test ( character a*72).
c key    = character*8  = key of matrix
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages
c         * ierr  =  0 means that  the matrix has been read normally.
c         * ierr  =  1 means that  the array matrix could not be read
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could not be read
c         because nnz .gt. nzmax
c         * ierr  =  3 means that  the array matrix could not be read
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...)
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of
c         rhs as specified by the input value of nrhs is not
c         sufficient to store them. The rest of the matrix may have
c         been read normally.
c
c Notes:
c-------
c 1) This routine can be interfaced with the C language, since only
c    the name of the  file needs to be passed and no iounti number.
c
c 2) Refer to the  documentation on  the Harwell-Boeing formats for
c    details on the format assumed by readmt.
c    We summarize the format here for convenience.
c
c    a) all  lines in  inout are  assumed to be 80  character long.
c    b) the file  consists of a header followed by the block of the
c       column start  pointers followed  by  the  block  of the row
c       indices,  followed  by  the  block  of the  real values and
c       finally the  numerical  values of  the right-hand-side if a
c       right hand side is supplied.
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first  line  contains  the  title  (72  characters  long)
c         followed  by  the  8-character  identifier (name  of  the
c         matrix, called key) [ A72,A8 ]
c       * second line  contains the number of lines for each of the
c         following data blocks (4 of them) and the total number of
c         lines excluding the header.  [5i4]
c       * the   third  line  contains  a  three   character  string
c         identifying the type of  matrices as they  are referenced
c         in  the Harwell  Boeing documentation [e.g., rua, rsa,..]
c         and the number of rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth  line contains the variable fortran format for
c         the following data blocks. [2A16,2A20]
c       * The fifth  line is  only present if  right-hand-sides are
c         supplied. It  consists  of  three  one  character-strings
c         containing the  storage  format for the  right-hand-sides
c         ('F'= full,'M'=sparse=same as matrix), an  initial  guess
c         indicator  ('G' for yes),  an  exact  solution  indicator
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices.  [A3,11X,2I14]
c     d) The three  following blocks follow the header as described
c        above.
c     e) In case the right hand-side are in sparse formats then the
c        fourth  block  uses the  same  storage  format as  for the
c        matrix to  describe  the NRHS right  hand  sides provided,
c        with a column being replaced by a right hand side.
c-----------------------------------------------------------------------
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     &     valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     &     nnz, neltvl, nrhs, nmax, nzmax, nrwindx
      integer ia (nmax+1), ja (nzmax)
      real*8 a(nzmax), rhs(*)
      character fname*100
c-----------------------------------------------------------------------
      ierr = 0
      lenrhs = nrhs
c
      iounit=15
      open(iounit,file = fname)
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd,
     &     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt,
     &     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs, nrwindx
 11   format (a3,11x,i14,i14)
c
c     anything else to read ?
c
      if (job .le. 0) goto 12
c     ---- check whether matrix is readable ------
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) goto 12
c     ---- read pointer and row numbers ----------
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)
c     --- reading values of matrix if required....
      if (job .le. 1)  goto 12
c     --- and if available -----------------------
      if (valcrd .le. 0) then
         job = 1
         goto 12
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)
c     --- reading rhs if required ----------------
      if (job .le. 2)  goto 12
c     --- and if available -----------------------
      if ( rhscrd .le. 0) then
         job = 2
         goto 12
      endif
c
c     --- read right-hand-side.--------------------
c
      if (rhstyp(1:1) .eq. 'M') then
         ierr = 4
         goto 12
      endif
c
      guesol = rhstyp(2:3)
c
      nvec = 1
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') nvec=nvec+1
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') nvec=nvec+1
c
      len = nrhs*nrow
c
      if (len*nvec .gt. lenrhs) then
        ierr = 5
        goto 12
      endif
c
c     read right-hand-sides
c
      next = 1
      iend = len
      read(iounit,rhsfmt) (rhs(i), i = next, iend)
c
c     read initial guesses if available
c
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
c     read exact solutions if available
c
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
 12   close(iounit)
      return
c---------end-of-readmt_c-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------

      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place.
c-----------------------------------------------------------------------
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n     = dimension of A.
c job   = integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
c         for any other normal usage, enter ipos=1.
c a     = real array of length nnz (nnz=number of nonzero elements in input
c         matrix) containing the nonzero elements.
c ja    = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia    = integer of size n+1. ia(k) contains the position in a, ja of
c         the beginning of the k-th row.
c
c on return:
c ----------
c output arguments:
c ao    = real array of size nzz containing the "a" part of the transpose
c jao   = integer array of size nnz containing the column indices.
c iao   = integer array of size n+1 containing the "ia" index array of
c         the transpose.
c
c-----------------------------------------------------------------------
      call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
c-----------------------------------------------------------------------
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place.
c-----------------------------------------------------------------------
c Rectangular version.  n is number of rows of CSR matrix,
c                       n2 (input) is number of columns of CSC matrix.
c-----------------------------------------------------------------------
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n     = number of rows of CSR matrix.
c n2    = number of columns of CSC matrix.
c job   = integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
c         for any other normal usage, enter ipos=1.
c a     = real array of length nnz (nnz=number of nonzero elements in input
c         matrix) containing the nonzero elements.
c ja    = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia    = integer of size n+1. ia(k) contains the position in a, ja of
c         the beginning of the k-th row.
c
c on return:
c ----------
c output arguments:
c ao    = real array of size nzz containing the "a" part of the transpose
c jao   = integer array of size nnz containing the column indices.
c iao   = integer array of size n+1 containing the "ia" index array of
c         the transpose.
c
c-----------------------------------------------------------------------
c----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue
 3    continue
c---------- compute pointers from lengths ------------------------------
      iao(1) = ipos
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c--------------- now do the actual copying -----------------------------
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1
            j = ja(k)
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c-------------------------- reshift iao and leave ----------------------
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
c--------------- end of csrcsc2 ----------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine gauss (n,a,ierr)
c-----------------------------------------------------------------------
      implicit none
      integer n, ierr
      real*8 a(n,n)
c
c     does the Gaussian factorization a := LU
c
      integer i, j, k
      real*8 piv
c-----------------------------------------------------------------------
      ierr = 0
      do k=1, n
         if (a(k,k) .eq. 0.0) then
            ierr = 1
            return
         endif
c
         a(k,k) = 1.0/a(k,k)
         do i=k+1, n
            piv = a(i,k) * a(k,k)
            do j=k+1, n
               a(i,j) = a(i,j) - piv*a(k,j)
            enddo
            a(i,k) = piv
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine bxinv (m, n, a, b, c)
      implicit none
      integer m, n
      real*8 a(n,n), b(m,n), c(m,n)
c
c     does the operation c :=  - b * inv(a)
c     where a has already been factored by Gauss.
c
      integer i, j, k
      real*8 sum
c
c     c = b (LU)**(-1) = b U\inv x L \inv
c     get c := b U \inv
c
      do i=1, m
c
c     U solve -- solve for row number k :   c(k,*) U = b(k,*)
c
         c(i,1) = - b(i,1) * a(1,1)
         do j=2, n
            sum = - b(i,j)
            do k=1, j-1
               sum = sum -  c(i,k) * a(k,j)
            enddo
            c(i,j) = sum*a(j,j)
         enddo
      enddo
c
      do i=1, m
c
c     L- solve -- solve for row number i :   c(i,*) U = b(i,*)
c
         do j=n-1, 1, -1
            sum = c(i,j)
            do k=j+1, n
               sum = sum -  c(i,k) * a(k,j)
            enddo
            c(i,j) = sum
         enddo
      enddo
c    
      end
c-----------------------------------------------------------------------

        subroutine qsplit(a,ind,n,ncut)
        real*8 a(n)
        integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end
c
      subroutine rnrms   (nrow, nrm, a, ia, diag) 
      real*8 a(*), diag(nrow), scal 
      integer ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each row of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c
c     compute the norm if each element.
c     
         scal = 0.0d0
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         if (nrm .eq. 0) then
            do 2 k=k1, k2
               scal = max(scal,abs(a(k) ) ) 
 2          continue
         elseif (nrm .eq. 1) then
            do 3 k=k1, k2
               scal = scal + abs(a(k) ) 
 3          continue
         else
            do 4 k=k1, k2
               scal = scal+a(k)**2
 4          continue
         endif 
         if (nrm .eq. 2) scal = sqrt(scal) 
         diag(ii) = scal
 1    continue
      return
c-----------------------------------------------------------------------
c-------------end-of-rnrms----------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine cnrms   (nrow, nrm, a, ja, ia, diag) 
      real*8 a(*), diag(nrow) 
      integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each column of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 10 k=1, nrow 
         diag(k) = 0.0d0
 10   continue
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            j = ja(k) 
c     update the norm of each column
            if (nrm .eq. 0) then
               diag(j) = max(diag(j),abs(a(k) ) ) 
            elseif (nrm .eq. 1) then
               diag(j) = diag(j) + abs(a(k) ) 
            else
               diag(j) = diag(j)+a(k)**2
            endif 
 2       continue
 1    continue
      if (nrm .ne. 2) return
      do 3 k=1, nrow
         diag(k) = sqrt(diag(k))
 3    continue
      return
c-----------------------------------------------------------------------
c------------end-of-cnrms-----------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine roscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
      real*8 a(*), b(*), diag(nrow) 
      integer nrow,job,nrm,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the rows of A such that their norms are one on return
c 3 choices of norms: 1-norm, 2-norm, max-norm.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the rows have been scaled, i.e., on return 
c        we have B = Diag*A.
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Row number i is a zero row.
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call rnrms (nrow,nrm,a,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0d0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call diamua(nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c-------end-of-roscal---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine coscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
c----------------------------------------------------------------------- 
      real*8 a(*),b(*),diag(nrow) 
      integer nrow,job,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the columns of A such that their norms are one on return
c result matrix written on b, or overwritten on A.
c 3 choices of norms: 1-norm, 2-norm, max-norm. in place.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the columns have been scaled, i.e., on return 
c        we have B = A * Diag
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Column number i is a zero row.
c Notes:
c-------
c 1)     The column dimension of A is not needed. 
c 2)     algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call cnrms (nrow,nrm,a,ja,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call amudia (nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c--------end-of-coscal-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = Diag * A  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c           in this case use job=0.
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     normalize each row 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii) 
         do 2 k=k1, k2
            b(k) = a(k)*scal
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c----------end-of-diamua------------------------------------------------
c-----------------------------------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine amudia (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = A * Diag  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     scale each element 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            b(k) = a(k)*diag(ja(k)) 
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c-----------------------------------------------------------------------
c-----------end-of-amudiag----------------------------------------------
      end 
c----------------------------------------------------------------------- 

      subroutine csrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
c-----------------------------------------------------------------------
      real*8 a(*),ao(*) 
      integer ir(*),jc(*),ja(*),ia(nrow+1) 
c----------------------------------------------------------------------- 
c  Compressed Sparse Row      to      Coordinate 
c----------------------------------------------------------------------- 
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry: 
c---------
c nrow	= dimension of the matrix.
c job   = integer serving as a job indicator. 
c         if job = 1 fill in only the array ir, ignore jc, and ao.
c         if job = 2 fill in ir, and jc but not ao 
c         if job = 3 fill in everything.
c         The reason why these options are provided is that on return 
c         ao and jc are the same as a, ja. So when job = 3, a and ja are
c         simply copied into ao, jc.  When job=2, only jc and ir are
c         returned. With job=1 only the array ir is returned. Moreover,
c         the algorithm is in place:
c	     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) 
c         will write the output matrix in coordinate format on a, ja,ia.
c
c a,
c ja,
c ia    = matrix in compressed sparse row format.
c nzmax = length of space available in ao, ir, jc.
c         the code will stop immediatly if the number of
c         nonzero elements found in input matrix exceeds nzmax.
c 
c on return:
c----------- 
c ao, ir, jc = matrix in coordinate format.
c
c nnz        = number of nonzero elements in matrix.
c ierr       = integer error indicator.
c         ierr .eq. 0 means normal retur
c         ierr .eq. 1 means that the the code stopped 
c         because there was no space in ao, ir, jc 
c         (according to the value of  nzmax).
c 
c NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with 
c         ao being the same array as as a, and jc the same array as ja. 
c         but ir CANNOT be the same as ia. 
c         2) note the order in the output arrays, 
c------------------------------------------------------------------------
      ierr = 0
      nnz = ia(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = 1
         return
      endif
c------------------------------------------------------------------------
      goto (3,2,1) job
 1    do 10 k=1,nnz
         ao(k) = a(k)
 10   continue
 2    do 11 k=1,nnz
         jc(k) = ja(k)
 11   continue
c
c     copy backward to allow for in-place processing. 
c
 3    do 13 i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do 12 k=k1,k2,-1
            ir(k) = i
 12      continue
 13   continue
      return
c------------- end-of-csrcoo ------------------------------------------- 
c----------------------------------------------------------------------- 
      end

c----------------------------------------------------------------------- 
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
c----------------------------------------------------------------------- 
      real*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
c-----------------------------------------------------------------------
c  Coordinate     to   Compressed Sparse Row 
c----------------------------------------------------------------------- 
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c--------- 
c nrow	= dimension of the matrix 
c nnz	= number of nonzero elements in matrix
c a,
c ir, 
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c 	  the elements, ir(k) = its row number and jc(k) = its column 
c	  number. The order of the elements is arbitrary. 
c
c on return:
c----------- 
c ir 	is destroyed
c
c ao, jao, iao = matrix in general sparse matrix format with ao 
c 	continung the real values, jao containing the column indices, 
c	and iao being the pointer to the beginning of the row, 
c	in arrays ao, jao.
c
c Notes:
c------ This routine is NOT in place.  See coicsr
c
c------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
c------------- end of coocsr ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
