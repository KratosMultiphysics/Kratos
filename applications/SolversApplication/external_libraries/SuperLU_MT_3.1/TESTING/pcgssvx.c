/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_cdefs.h"


void
pcgssvx(int_t nprocs, superlumt_options_t *superlumt_options, SuperMatrix *A, 
	int_t *perm_c, int_t *perm_r, equed_t *equed, float *R, float *C,
	SuperMatrix *L, SuperMatrix *U,
	SuperMatrix *B, SuperMatrix *X, float *recip_pivot_growth, 
	float *rcond, float *ferr, float *berr, 
	superlu_memusage_t *superlu_memusage, int_t *info)
{
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley, 
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Purpose
 * =======
 *
 * pcgssvx() solves the system of linear equations A*X=B or A'*X=B, using
 * the LU factorization from cgstrf(). Error bounds on the solution and
 * a condition estimate are also provided. It performs the following steps:
 *
 * 1. If A is stored column-wise (A->Stype = NC):
 *  
 *    1.1. If fact = EQUILIBRATE, scaling factors are computed to equilibrate
 *         the system:
 *           trans = NOTRANS: diag(R)*A*diag(C)*inv(diag(C))*X = diag(R)*B
 *           trans = TRANS:  (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
 *           trans = CONJ:   (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
 *         Whether or not the system will be equilibrated depends on the
 *         scaling of the matrix A, but if equilibration is used, A is
 *         overwritten by diag(R)*A*diag(C) and B by diag(R)*B 
 *         (if trans = NOTRANS) or diag(C)*B (if trans = TRANS or CONJ).
 *
 *    1.2. Permute columns of A, forming A*Pc, where Pc is a permutation matrix
 *         that usually preserves sparsity.
 *         For more details of this step, see csp_colorder.c.
 *
 *    1.3. If fact = DOFACT or EQUILIBRATE, the LU decomposition is used to 
 *         factor the matrix A (after equilibration if fact = EQUILIBRATE) as
 *         Pr*A*Pc = L*U, with Pr determined by partial pivoting.
 *
 *    1.4. Compute the reciprocal pivot growth factor.
 *
 *    1.5. If some U(i,i) = 0, so that U is exactly singular, then the routine
 *         returns with info = i. Otherwise, the factored form of A is used to
 *         estimate the condition number of the matrix A. If the reciprocal of
 *         the condition number is less than machine precision, 
 *         info = A->ncol+1 is returned as a warning, but the routine still
 *         goes on to solve for X and computes error bounds as described below.
 *
 *    1.6. The system of equations is solved for X using the factored form
 *         of A.
 *
 *    1.7. Iterative refinement is applied to improve the computed solution
 *         matrix and calculate error bounds and backward error estimates
 *         for it.
 *
 *    1.8. If equilibration was used, the matrix X is premultiplied by
 *         diag(C) (if trans = NOTRANS) or diag(R) (if trans = TRANS or CONJ)
 *         so that it solves the original system before equilibration.
 *
 * 2. If A is stored row-wise (A->Stype = NR), apply the above algorithm
 *    to the tranpose of A:
 *
 *    2.1. If fact = EQUILIBRATE, scaling factors are computed to equilibrate
 *         the system:
 *           trans = NOTRANS:diag(R)*A'*diag(C)*inv(diag(C))*X = diag(R)*B
 *           trans = TRANS: (diag(R)*A'*diag(C))**T *inv(diag(R))*X = diag(C)*B
 *           trans = CONJ:  (diag(R)*A'*diag(C))**H *inv(diag(R))*X = diag(C)*B
 *         Whether or not the system will be equilibrated depends on the
 *         scaling of the matrix A, but if equilibration is used, A' is
 *         overwritten by diag(R)*A'*diag(C) and B by diag(R)*B 
 *         (if trans = NOTRANS) or diag(C)*B (if trans = TRANS or CONJ).
 *
 *    2.2. Permute columns of transpose(A) (rows of A), 
 *         forming transpose(A)*Pc, where Pc is a permutation matrix that
 *         usually preserves sparsity.
 *         For more details of this step, see csp_colorder.c.
 *
 *    2.3. If fact = DOFACT or EQUILIBRATE, the LU decomposition is used to 
 *         factor the matrix A (after equilibration if fact = EQUILIBRATE) as
 *         Pr*transpose(A)*Pc = L*U, with the permutation Pr determined by
 *         partial pivoting.
 *
 *    2.4. Compute the reciprocal pivot growth factor.
 *
 *    2.5. If some U(i,i) = 0, so that U is exactly singular, then the routine
 *         returns with info = i. Otherwise, the factored form of transpose(A)
 *         is used to estimate the condition number of the matrix A.
 *         If the reciprocal of the condition number is less than machine
 *         precision, info = A->nrow+1 is returned as a warning, but the
 *         routine still goes on to solve for X and computes error bounds
 *         as described below.
 *
 *    2.6. The system of equations is solved for X using the factored form
 *         of transpose(A).
 *
 *    2.7. Iterative refinement is applied to improve the computed solution
 *         matrix and calculate error bounds and backward error estimates
 *         for it.
 *
 *    2.8. If equilibration was used, the matrix X is premultiplied by
 *         diag(C) (if trans = NOTRANS) or diag(R) (if trans = TRANS or CONJ)
 *         so that it solves the original system before equilibration.
 *
 * See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * nprocs (input) int_t
 *         Number of processes (or threads) to be spawned and used to perform
 *         the LU factorization by pcgstrf(). There is a single thread of
 *         control to call pcgstrf(), and all threads spawned by pcgstrf() 
 *         are terminated before returning from pcgstrf().
 *
 * superlumt_options (input) superlumt_options_t*
 *         The structure defines the input parameters and data structure
 *         to control how the LU factorization will be performed.
 *         The following fields should be defined for this structure:
 *
 *         o fact (fact_t)
 *           Specifies whether or not the factored form of the matrix
 *           A is supplied on entry, and if not, whether the matrix A should
 *           be equilibrated before it is factored.
 *           = DOFACT: The matrix A will be factored, and the factors will be
 *             stored in L and U.
 *           = EQUILIBRATE: The matrix A will be equilibrated if necessary,
 *             then factored into L and U.
 *           = FACTORED: On entry, L, U, perm_r and perm_c contain the 
 *             factored form of A. If equed is not NOEQUIL, the matrix A has
 *             been equilibrated with scaling factors R and C.
 *             A, L, U, perm_r are not modified.
 *
 *         o trans (trans_t)
 *           Specifies the form of the system of equations:
 *           = NOTRANS: A * X = B        (No transpose)
 *           = TRANS:   A**T * X = B     (Transpose)
 *           = CONJ:    A**H * X = B     (Transpose)
 *
 *         o refact (yes_no_t)
 *           Specifies whether this is first time or subsequent factorization.
 *           = NO:  this factorization is treated as the first one;
 *           = YES: it means that a factorization was performed prior to this
 *               one. Therefore, this factorization will re-use some
 *               existing data structures, such as L and U storage, column
 *               elimination tree, and the symbolic information of the
 *               Householder matrix.
 *
 *         o panel_size (int_t)
 *           A panel consists of at most panel_size consecutive columns.
 *
 *         o relax (int_t)
 *           To control degree of relaxing supernodes. If the number
 *           of nodes (columns) in a subtree of the elimination tree is less
 *           than relax, this subtree is considered as one supernode,
 *           regardless of the row structures of those columns.
 *
 *         o diag_pivot_thresh (float)
 *           Diagonal pivoting threshold. At step j of the Gaussian 
 *           elimination, if 
 *               abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
 *           use A_jj as pivot, else use A_ij with maximum magnitude. 
 *           0 <= diag_pivot_thresh <= 1. The default value is 1, 
 *           corresponding to partial pivoting.
 *
 *         o usepr (yes_no_t)
 *           Whether the pivoting will use perm_r specified by the user.
 *           = YES: use perm_r; perm_r is input, unchanged on exit.
 *           = NO:  perm_r is determined by partial pivoting, and is output.
 *
 *         o drop_tol (double) (NOT IMPLEMENTED)
 *	     Drop tolerance parameter. At step j of the Gaussian elimination,
 *           if abs(A_ij)/(max_i abs(A_ij)) < drop_tol, drop entry A_ij.
 *           0 <= drop_tol <= 1. The default value of drop_tol is 0,
 *           corresponding to not dropping any entry.
 *
 *         o work (void*) of size lwork
 *           User-supplied work space and space for the output data structures.
 *           Not referenced if lwork = 0;
 *
 *         o lwork (int_t)
 *           Specifies the length of work array.
 *           = 0:  allocate space internally by system malloc;
 *           > 0:  use user-supplied work array of length lwork in bytes,
 *                 returns error if space runs out.
 *           = -1: the routine guesses the amount of space needed without
 *                 performing the factorization, and returns it in
 *                 superlu_memusage->total_needed; no other side effects.
 *
 *    -------------------------------------------------------------------------
 *    ** The following 3 arrays must be allocated before calling this routine.
 *    -------------------------------------------------------------------------
 *         o etree  (int*)
 *           Elimination tree of A'*A, dimension A->ncol.
 *
 *         o colcnt_h (int*)
 *           Column colunts of the Householder matrix.
 *
 *         o part_super_h (int*)
 *           Partition of the supernodes in the Householder matrix.
 *	     part_super_h[k] = size of the supernode beginning at column k;
 * 	                     = 0, elsewhere.
 *
 * A       (input/output) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol), where
 *         A->nrow = A->ncol. Currently, the type of A can be:
 *         Stype = NC or NR, Dtype = _D, Mtype = GE. In the future,
 *         more general A will be handled.
 *
 *         On entry, If superlumt_options->fact = FACTORED and equed is not 
 *         NOEQUIL, then A must have been equilibrated by the scaling factors
 *         in R and/or C.  On exit, A is not modified 
 *         if superlumt_options->fact = FACTORED or DOFACT, or 
 *         if superlumt_options->fact = EQUILIBRATE and equed = NOEQUIL.
 *
 *         On exit, if superlumt_options->fact = EQUILIBRATE and equed is not
 *         NOEQUIL, A is scaled as follows:
 *         If A->Stype = NC:
 *           equed = ROW:  A := diag(R) * A
 *           equed = COL:  A := A * diag(C)
 *           equed = BOTH: A := diag(R) * A * diag(C).
 *         If A->Stype = NR:
 *           equed = ROW:  transpose(A) := diag(R) * transpose(A)
 *           equed = COL:  transpose(A) := transpose(A) * diag(C)
 *           equed = BOTH: transpose(A) := diag(R) * transpose(A) * diag(C).
 *
 * perm_c  (input/output) int_t*
 *	   If A->Stype = NC, Column permutation vector of size A->ncol,
 *         which defines the permutation matrix Pc; perm_c[i] = j means
 *         column i of A is in position j in A*Pc.
 *         On exit, perm_c may be overwritten by the product of the input
 *         perm_c and a permutation that postorders the elimination tree
 *         of Pc'*A'*A*Pc; perm_c is not changed if the elimination tree
 *         is already in postorder.
 *
 *         If A->Stype = NR, column permutation vector of size A->nrow,
 *         which describes permutation of columns of tranpose(A) 
 *         (rows of A) as described above.
 * 
 * perm_r  (input/output) int_t*
 *         If A->Stype = NC, row permutation vector of size A->nrow, 
 *         which defines the permutation matrix Pr, and is determined
 *         by partial pivoting.  perm_r[i] = j means row i of A is in 
 *         position j in Pr*A.
 *
 *         If A->Stype = NR, permutation vector of size A->ncol, which
 *         determines permutation of rows of transpose(A)
 *         (columns of A) as described above.
 *
 *         If superlumt_options->usepr = NO, perm_r is output argument;
 *         If superlumt_options->usepr = YES, the pivoting routine will try 
 *            to use the input perm_r, unless a certain threshold criterion
 *            is violated. In that case, perm_r is overwritten by a new
 *            permutation determined by partial pivoting or diagonal 
 *            threshold pivoting.
 * 
 * equed   (input/output) equed_t*
 *         Specifies the form of equilibration that was done.
 *         = NOEQUIL: No equilibration.
 *         = ROW:  Row equilibration, i.e., A was premultiplied by diag(R).
 *         = COL:  Column equilibration, i.e., A was postmultiplied by diag(C).
 *         = BOTH: Both row and column equilibration, i.e., A was replaced 
 *                 by diag(R)*A*diag(C).
 *         If superlumt_options->fact = FACTORED, equed is an input argument, 
 *         otherwise it is an output argument.
 *
 * R       (input/output) double*, dimension (A->nrow)
 *         The row scale factors for A or transpose(A).
 *         If equed = ROW or BOTH, A (if A->Stype = NC) or transpose(A)
 *            (if A->Stype = NR) is multiplied on the left by diag(R).
 *         If equed = NOEQUIL or COL, R is not accessed.
 *         If fact = FACTORED, R is an input argument; otherwise, R is output.
 *         If fact = FACTORED and equed = ROW or BOTH, each element of R must
 *            be positive.
 * 
 * C       (input/output) double*, dimension (A->ncol)
 *         The column scale factors for A or transpose(A).
 *         If equed = COL or BOTH, A (if A->Stype = NC) or trnspose(A)
 *            (if A->Stype = NR) is multiplied on the right by diag(C).
 *         If equed = NOEQUIL or ROW, C is not accessed.
 *         If fact = FACTORED, C is an input argument; otherwise, C is output.
 *         If fact = FACTORED and equed = COL or BOTH, each element of C must
 *            be positive.
 *         
 * L       (output) SuperMatrix*
 *	   The factor L from the factorization
 *             Pr*A*Pc=L*U              (if A->Stype = NC) or
 *             Pr*transpose(A)*Pc=L*U   (if A->Stype = NR).
 *         Uses compressed row subscripts storage for supernodes, i.e.,
 *         L has types: Stype = SCP, Dtype = _D, Mtype = TRLU.
 *
 * U       (output) SuperMatrix*
 *	   The factor U from the factorization
 *             Pr*A*Pc=L*U              (if A->Stype = NC) or
 *             Pr*transpose(A)*Pc=L*U   (if A->Stype = NR).
 *         Uses column-wise storage scheme, i.e., U has types:
 *         Stype = NCP, Dtype = _D, Mtype = TRU.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = DN, Dtype = _D, Mtype = GE.
 *         On entry, the right hand side matrix.
 *         On exit,
 *            if equed = NOEQUIL, B is not modified; otherwise
 *            if A->Stype = NC:
 *               if trans = NOTRANS and equed = ROW or BOTH, B is overwritten
 *                  by diag(R)*B;
 *               if trans = TRANS or CONJ and equed = COL of BOTH, B is
 *                  overwritten by diag(C)*B;
 *            if A->Stype = NR:
 *               if trans = NOTRANS and equed = COL or BOTH, B is overwritten
 *                  by diag(C)*B;
 *               if trans = TRANS or CONJ and equed = ROW of BOTH, B is
 *                  overwritten by diag(R)*B.
 *
 * X       (output) SuperMatrix*
 *         X has types: Stype = DN, Dtype = _D, Mtype = GE. 
 *         If info = 0 or info = A->ncol+1, X contains the solution matrix
 *         to the original system of equations. Note that A and B are modified
 *         on exit if equed is not NOEQUIL, and the solution to the 
 *         equilibrated system is inv(diag(C))*X if trans = NOTRANS and
 *         equed = COL or BOTH, or inv(diag(R))*X if trans = TRANS or CONJ
 *         and equed = ROW or BOTH.
 *
 * recip_pivot_growth (output) float*
 *         The reciprocal pivot growth factor computed as
 *             max_j ( max_i(abs(A_ij)) / max_i(abs(U_ij)) ).
 *         If recip_pivot_growth is much less than 1, the stability of the
 *         LU factorization could be poor.
 *
 * rcond   (output) float*
 *         The estimate of the reciprocal condition number of the matrix A
 *         after equilibration (if done). If rcond is less than the machine
 *         precision (in particular, if rcond = 0), the matrix is singular
 *         to working precision. This condition is indicated by a return
 *         code of info > 0.
 *
 * ferr    (output) float*, dimension (B->ncol)   
 *         The estimated forward error bound for each solution vector   
 *         X(j) (the j-th column of the solution matrix X).   
 *         If XTRUE is the true solution corresponding to X(j), FERR(j) 
 *         is an estimated upper bound for the magnitude of the largest 
 *         element in (X(j) - XTRUE) divided by the magnitude of the   
 *         largest element in X(j).  The estimate is as reliable as   
 *         the estimate for RCOND, and is almost always a slight   
 *         overestimate of the true error.
 *
 * berr    (output) float*, dimension (B->ncol)
 *         The componentwise relative backward error of each solution   
 *         vector X(j) (i.e., the smallest relative change in   
 *         any element of A or B that makes X(j) an exact solution).
 *
 * superlu_memusage (output) superlu_memusage_t*
 *         Record the memory usage statistics, consisting of following fields:
 *         - for_lu (float)
 *           The amount of space used in bytes for L\U data structures.
 *         - total_needed (float)
 *           The amount of space needed in bytes to perform factorization.
 *         - expansions (int_t)
 *           The number of memory expansions during the LU factorization.
 *
 * info    (output) int_t*
 *         = 0: successful exit   
 *         < 0: if info = -i, the i-th argument had an illegal value   
 *         > 0: if info = i, and i is   
 *              <= A->ncol: U(i,i) is exactly zero. The factorization has   
 *                    been completed, but the factor U is exactly   
 *                    singular, so the solution and error bounds   
 *                    could not be computed.   
 *              = A->ncol+1: U is nonsingular, but RCOND is less than machine
 *                    precision, meaning that the matrix is singular to
 *                    working precision. Nevertheless, the solution and
 *                    error bounds are computed because there are a number
 *                    of situations where the computed solution can be more
 *                    accurate than the value of RCOND would suggest.   
 *              > A->ncol+1: number of bytes allocated when memory allocation
 *                    failure occurred, plus A->ncol.
 *
 */

    NCformat  *Astore;
    DNformat  *Bstore, *Xstore;
    complex    *Bmat, *Xmat;
    int_t       ldb, ldx, nrhs;
    SuperMatrix *AA; /* A in NC format used by the factorization routine.*/
    SuperMatrix AC; /* Matrix postmultiplied by Pc */
    int_t       colequ, equil, dofact, notran, rowequ;
    char      norm[1];
    trans_t   trant;
    int_t     j, info1;
    int i;
    float amax, anorm, bignum, smlnum, colcnd, rowcnd, rcmax, rcmin;
    int_t       n, relax, panel_size;
    Gstat_t   Gstat;
    double    t0;      /* temporary time */
    double    *utime;
    flops_t   *ops, flopcnt;
   
    /* External functions */
    extern float clangs(char *, SuperMatrix *);
    extern double slamch_(char *);

    Astore = A->Store;
    Bstore = B->Store;
    Xstore = X->Store;
    Bmat   = Bstore->nzval;
    Xmat   = Xstore->nzval;
    n      = A->ncol;
    ldb    = Bstore->lda;
    ldx    = Xstore->lda;
    nrhs   = B->ncol;
    superlumt_options->perm_c = perm_c;
    superlumt_options->perm_r = perm_r;

    *info = 0;
    dofact = (superlumt_options->fact == DOFACT);
    equil = (superlumt_options->fact == EQUILIBRATE);
    notran = (superlumt_options->trans == NOTRANS);
    if (dofact || equil) {
	*equed = NOEQUIL;
	rowequ = FALSE;
	colequ = FALSE;
    } else {
	rowequ = (*equed == ROW) || (*equed == BOTH);
	colequ = (*equed == COL) || (*equed == BOTH);
	smlnum = slamch_("Safe minimum");
	bignum = 1. / smlnum;
    }

    /* ------------------------------------------------------------
       Test the input parameters.
       ------------------------------------------------------------*/
    if ( nprocs <= 0 ) *info = -1;
    else if ( (!dofact && !equil && (superlumt_options->fact != FACTORED))
	      || (!notran && (superlumt_options->trans != TRANS) && 
		 (superlumt_options->trans != CONJ))
	      || (superlumt_options->refact != YES && 
		  superlumt_options->refact != NO)
	      || (superlumt_options->usepr != YES &&
		  superlumt_options->usepr != NO)
	      || superlumt_options->lwork < -1 )
        *info = -2;
    else if ( A->nrow != A->ncol || A->nrow < 0 ||
	      (A->Stype != SLU_NC && A->Stype != SLU_NR) ||
	      A->Dtype != SLU_C || A->Mtype != SLU_GE )
	*info = -3;
    else if ((superlumt_options->fact == FACTORED) && 
	     !(rowequ || colequ || (*equed == NOEQUIL))) *info = -6;
    else {
	if (rowequ) {
	    rcmin = bignum;
	    rcmax = 0.;
	    for (j = 0; j < A->nrow; ++j) {
		rcmin = SUPERLU_MIN(rcmin, R[j]);
		rcmax = SUPERLU_MAX(rcmax, R[j]);
	    }
	    if (rcmin <= 0.) *info = -7;
	    else if ( A->nrow > 0)
		rowcnd = SUPERLU_MAX(rcmin,smlnum) / SUPERLU_MIN(rcmax,bignum);
	    else rowcnd = 1.;
	}
	if (colequ && *info == 0) {
	    rcmin = bignum;
	    rcmax = 0.;
	    for (j = 0; j < A->nrow; ++j) {
		rcmin = SUPERLU_MIN(rcmin, C[j]);
		rcmax = SUPERLU_MAX(rcmax, C[j]);
	    }
	    if (rcmin <= 0.) *info = -8;
	    else if (A->nrow > 0)
		colcnd = SUPERLU_MAX(rcmin,smlnum) / SUPERLU_MIN(rcmax,bignum);
	    else colcnd = 1.;
	}
	if (*info == 0) {
	    if ( B->ncol < 0 || Bstore->lda < SUPERLU_MAX(0, A->nrow) ||
		      B->Stype != SLU_DN || B->Dtype != SLU_C || 
		      B->Mtype != SLU_GE )
		*info = -11;
	    else if ( X->ncol < 0 || Xstore->lda < SUPERLU_MAX(0, A->nrow) ||
		      B->ncol != X->ncol || X->Stype != SLU_DN ||
		      X->Dtype != SLU_C || X->Mtype != SLU_GE )
		*info = -12;
	}
    }
    if (*info != 0) {
	i = -(*info);
	xerbla_("pcgssvx", &i);
	return;
    }
    
    
    /* ------------------------------------------------------------
       Allocate storage and initialize statistics variables. 
       ------------------------------------------------------------*/
    panel_size = superlumt_options->panel_size;
    relax = superlumt_options->relax;
    StatAlloc(n, nprocs, panel_size, relax, &Gstat);
    StatInit(n, nprocs, &Gstat);
    utime = Gstat.utime;
    ops = Gstat.ops;
    
    /* ------------------------------------------------------------
       Convert A to NC format when necessary.
       ------------------------------------------------------------*/
    if ( A->Stype == SLU_NR ) {
	NRformat *Astore = A->Store;
	AA = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	cCreate_CompCol_Matrix(AA, A->ncol, A->nrow, Astore->nnz, 
			       Astore->nzval, Astore->colind, Astore->rowptr,
			       SLU_NC, A->Dtype, A->Mtype);
	if ( notran ) { /* Reverse the transpose argument. */
	    trant = TRANS;
	    notran = 0;
	} else {
	    trant = NOTRANS;
	    notran = 1;
	}
    } else { /* A->Stype == NC */
	trant = superlumt_options->trans;
	AA = A;
    }

    /* ------------------------------------------------------------
       Diagonal scaling to equilibrate the matrix.
       ------------------------------------------------------------*/
    if ( equil ) {
	t0 = SuperLU_timer_();
	/* Compute row and column scalings to equilibrate the matrix A. */
	cgsequ(AA, R, C, &rowcnd, &colcnd, &amax, &info1);
	
	if ( info1 == 0 ) {
	    /* Equilibrate matrix A. */
	    claqgs(AA, R, C, rowcnd, colcnd, amax, equed);
	    rowequ = (*equed == ROW) || (*equed == BOTH);
	    colequ = (*equed == COL) || (*equed == BOTH);
	}
	utime[EQUIL] = SuperLU_timer_() - t0;
    }

    /* ------------------------------------------------------------
       Scale the right hand side.
       ------------------------------------------------------------*/
    if ( notran ) {
	if ( rowequ ) {
	    for (j = 0; j < nrhs; ++j)
		for (i = 0; i < A->nrow; ++i) {
                        cs_mult(&Bmat[i+j*ldb], &Bmat[i+j*ldb], R[i]);
		}
	}
    } else if ( colequ ) {
	for (j = 0; j < nrhs; ++j)
	    for (i = 0; i < A->nrow; ++i) {
                    cs_mult(&Bmat[i+j*ldb], &Bmat[i+j*ldb], C[i]);
	    }
    }

    
    /* ------------------------------------------------------------
       Perform the LU factorization.
       ------------------------------------------------------------*/
    if ( dofact || equil ) {
	
        /* Obtain column etree, the column count (colcnt_h) and supernode
	   partition (part_super_h) for the Householder matrix. */
	t0 = SuperLU_timer_();
	sp_colorder(AA, perm_c, superlumt_options, &AC);
	utime[ETREE] = SuperLU_timer_() - t0;

#if ( PRNTlevel >= 2 )    
	printf("Factor PA = LU ... relax %d\tw %d\tmaxsuper %d\trowblk %d\n", 
	       relax, panel_size, sp_ienv(3), sp_ienv(4));
	fflush(stdout);
#endif
	
	/* Compute the LU factorization of A*Pc. */
	t0 = SuperLU_timer_();
	pcgstrf(superlumt_options, &AC, perm_r, L, U, &Gstat, info);
	utime[FACT] = SuperLU_timer_() - t0;
	
	flopcnt = 0;
	for (i = 0; i < nprocs; ++i) flopcnt += Gstat.procstat[i].fcops;
	ops[FACT] = flopcnt;

	if ( superlumt_options->lwork == -1 ) {
	    superlu_memusage->total_needed = *info - A->ncol;
	    return;
	}
    }

    if ( *info > 0 ) {
	if ( *info <= A->ncol ) {
	    /* Compute the reciprocal pivot growth factor of the leading
	       rank-deficient *info columns of A. */
	    *recip_pivot_growth = cPivotGrowth(*info, AA, perm_c, L, U);
	}
    } else {

	/* ------------------------------------------------------------
	   Compute the reciprocal pivot growth factor *recip_pivot_growth.
	   ------------------------------------------------------------*/
	*recip_pivot_growth = cPivotGrowth(A->ncol, AA, perm_c, L, U);

	/* ------------------------------------------------------------
	   Estimate the reciprocal of the condition number of A.
	   ------------------------------------------------------------*/
	t0 = SuperLU_timer_();
	if ( notran ) {
	    *(unsigned char *)norm = '1';
	} else {
	    *(unsigned char *)norm = 'I';
	}
	anorm = clangs(norm, AA);
	cgscon(norm, L, U, anorm, rcond, info);
	utime[RCOND] = SuperLU_timer_() - t0;
    
	/* ------------------------------------------------------------
	   Compute the solution matrix X.
	   ------------------------------------------------------------*/
	for (j = 0; j < nrhs; j++)    /* Save a copy of the right hand sides */
	    for (i = 0; i < B->nrow; i++)
		Xmat[i + j*ldx] = Bmat[i + j*ldb];
    
	t0 = SuperLU_timer_();
	cgstrs(trant, L, U, perm_r, perm_c, X, &Gstat, info);
	utime[SOLVE] = SuperLU_timer_() - t0;
	ops[SOLVE] = ops[TRISOLVE];
    
	/* ------------------------------------------------------------
	   Use iterative refinement to improve the computed solution and
	   compute error bounds and backward error estimates for it.
	   ------------------------------------------------------------*/
	t0 = SuperLU_timer_();
	cgsrfs(trant, AA, L, U, perm_r, perm_c, *equed,
	       R, C, B, X, ferr, berr, &Gstat, info);
	utime[REFINE] = SuperLU_timer_() - t0;

	/* ------------------------------------------------------------
	   Transform the solution matrix X to a solution of the original
	   system.
	   ------------------------------------------------------------*/
	if ( notran ) {
	    if ( colequ ) {
		for (j = 0; j < nrhs; ++j)
		    for (i = 0; i < A->nrow; ++i) {
                        cs_mult(&Xmat[i+j*ldx], &Xmat[i+j*ldx], C[i]);
		    }
	    }
	} else if ( rowequ ) {
	    for (j = 0; j < nrhs; ++j)
		for (i = 0; i < A->nrow; ++i) {
                    cs_mult(&Xmat[i+j*ldx], &Xmat[i+j*ldx], R[i]);
		}
	}
	
	/* Set INFO = A->ncol+1 if the matrix is singular to 
	   working precision.*/
	if ( *rcond < slamch_("E") ) *info = A->ncol + 1;
	
    }

    superlu_cQuerySpace(nprocs, L, U, panel_size, superlu_memusage);

    /* ------------------------------------------------------------
       Deallocate storage after factorization.
       ------------------------------------------------------------*/
    if ( dofact || equil ) {
        Destroy_CompCol_Permuted(&AC);
    }
    if ( A->Stype == SLU_NR ) {
	Destroy_SuperMatrix_Store(AA);
	SUPERLU_FREE(AA);
    }

    /* ------------------------------------------------------------
       Print timings, then deallocate statistic variables.
       ------------------------------------------------------------*/
#ifdef PROFILE
    {
	SCPformat *Lstore = (SCPformat *) L->Store;
	ParallelProfile(n, Lstore->nsuper+1, Gstat.num_panels, nprocs, &Gstat);
    }
#endif
    /*PrintStat(&Gstat);*/
    StatFree(&Gstat);
}
