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
pcgssv(int_t nprocs, SuperMatrix *A, int_t *perm_c, int_t *perm_r, 
       SuperMatrix *L, SuperMatrix *U, SuperMatrix *B, int_t *info )
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
 * PCGSSV solves the system of linear equations A*X=B, using the parallel
 * LU factorization routine PCGSTRF. It performs the following steps:
 *
 *   1. If A is stored column-wise (A->Stype = NC):
 *
 *      1.1. Permute the columns of A, forming A*Pc, where Pc is a 
 *           permutation matrix. 
 *           For more details of this step, see sp_preorder.c.
 *
 *      1.2. Factor A as Pr*A*Pc=L*U with the permutation Pr determined
 *           by Gaussian elimination with partial pivoting.
 *           L is unit lower triangular with offdiagonal entries
 *           bounded by 1 in magnitude, and U is upper triangular.
 *
 *      1.3. Solve the system of equations A*X=B using the factored
 *           form of A.
 *
 *   2. If A is stored row-wise (A->Stype = NR), apply the above algorithm
 *      to the tranpose of A:
 *
 *      2.1. Permute columns of tranpose(A) (rows of A),
 *           forming transpose(A)*Pc, where Pc is a permutation matrix. 
 *           For more details of this step, see sp_preorder.c.
 *
 *      2.2. Factor A as Pr*transpose(A)*Pc=L*U with the permutation Pr
 *           determined by Gaussian elimination with partial pivoting.
 *           L is unit lower triangular with offdiagonal entries
 *           bounded by 1 in magnitude, and U is upper triangular.
 *
 *      2.3. Solve the system of equations A*X=B using the factored
 *           form of A.
 * 
 *   See supermatrix.h for the definition of "SuperMatrix" structure.
 *
 *
 * Arguments
 * =========
 *
 * nprocs (input) int_t
 *        Number of processes (or threads) to be spawned and used to perform
 *        the LU factorization by pcgstrf(). There is a single thread of
 *        control to call pcgstrf(), and all threads spawned by pcgstrf()
 *        are terminated before returning from pcgstrf().
 *
 * A      (input) SuperMatrix*
 *        Matrix A in A*X=B, of dimension (A->nrow, A->ncol), where
 *        A->nrow = A->ncol. Currently, the type of A can be:
 *        Stype = NC or NR; Dtype = _D; Mtype = GE. In the future,
 *        more general A will be handled.
 *
 * perm_c (input/output) int_t*
 *        If A->Stype=NC, column permutation vector of size A->ncol,
 *        which defines the permutation matrix Pc; perm_c[i] = j means 
 *        column i of A is in position j in A*Pc.
 *        On exit, perm_c may be overwritten by the product of the input
 *        perm_c and a permutation that postorders the elimination tree
 *        of Pc'*A'*A*Pc; perm_c is not changed if the elimination tree
 *        is already in postorder.
 *
 *        If A->Stype=NR, column permutation vector of size A->nrow
 *        which describes permutation of columns of tranpose(A) 
 *        (rows of A) as described above.
 * 
 * perm_r (output) int_t*,
 *        If A->Stype=NR, row permutation vector of size A->nrow, 
 *        which defines the permutation matrix Pr, and is determined 
 *        by partial pivoting.  perm_r[i] = j means row i of A is in 
 *        position j in Pr*A.
 *
 *        If A->Stype=NR, permutation vector of size A->ncol, which
 *        determines permutation of rows of transpose(A)
 *        (columns of A) as described above.
 *
 * L      (output) SuperMatrix*
 *        The factor L from the factorization 
 *            Pr*A*Pc=L*U              (if A->Stype=NC) or
 *            Pr*transpose(A)*Pc=L*U   (if A->Stype=NR).
 *        Uses compressed row subscripts storage for supernodes, i.e.,
 *        L has types: Stype = SCP, Dtype = _D, Mtype = TRLU.
 *
 * U      (output) SuperMatrix*
 *	  The factor U from the factorization
 *            Pr*A*Pc=L*U              (if A->Stype=NC) or
 *            Pr*transpose(A)*Pc=L*U   (if A->Stype=NR).
 *        Use column-wise storage scheme, i.e., U has types:
 *        Stype = NCP, Dtype = _D, Mtype = TRU.
 *
 * B      (input/output) SuperMatrix*
 *        B has types: Stype = DN, Dtype = _D, Mtype = GE.
 *        On entry, the right hand side matrix.
 *        On exit, the solution matrix if info = 0;
 *
 * info   (output) int_t*
 *	  = 0: successful exit
 *        > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 *   
 */
    trans_t  trans;
    NCformat *Astore;
    DNformat *Bstore;
    SuperMatrix *AA; /* A in NC format used by the factorization routine.*/
    SuperMatrix AC; /* Matrix postmultiplied by Pc */
    int_t  n, panel_size, relax;
    int i;
    fact_t   fact;
    yes_no_t refact, usepr;
    float diag_pivot_thresh, drop_tol;
    void *work;
    int_t lwork;
    superlumt_options_t superlumt_options;
    Gstat_t  Gstat;
    double   t; /* Temporary time */
    double   *utime;
    flops_t  *ops, flopcnt;

    /* ------------------------------------------------------------
       Test the input parameters.
       ------------------------------------------------------------*/
    Astore = A->Store;
    Bstore = B->Store;
    *info = 0;
    if ( nprocs <= 0 ) *info = -1;
    else if ( A->nrow != A->ncol || A->nrow < 0 || 
	      (A->Stype != SLU_NC && A->Stype != SLU_NR) ||
	      A->Dtype != SLU_C || A->Mtype != SLU_GE )
	*info = -2;
    else if ( B->ncol < 0 || Bstore->lda < SUPERLU_MAX(1, A->nrow) )*info = -7;
    if ( *info != 0 ) {
        i = -(*info);
	xerbla_("pcgssv", &i);
	return;
    }

#if 0
    /* Use the best sequential code. 
       if this part is commented out, we will use the parallel code 
       run on one processor. */
    if ( nprocs == 1 ) {
        return;
    }
#endif

    fact               = EQUILIBRATE;
    refact             = NO;
    trans              = NOTRANS;
    panel_size         = sp_ienv(1);
    relax              = sp_ienv(2);
    diag_pivot_thresh  = 1.0;
    usepr              = NO;
    drop_tol           = 0.0;
    work               = NULL;
    lwork              = 0;

    /* ------------------------------------------------------------
       Allocate storage and initialize statistics variables. 
       ------------------------------------------------------------*/
    n = A->ncol;
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
	trans = TRANS;
    } else if ( A->Stype == SLU_NC ) AA = A;

    /* ------------------------------------------------------------
       Initialize the option structure superlumt_options using the
       user-input parameters;
       Apply perm_c to the columns of original A to form AC.
       ------------------------------------------------------------*/
    pcgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
		 diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r,
		 work, lwork, AA, &AC, &superlumt_options, &Gstat);

    /* ------------------------------------------------------------
       Compute the LU factorization of A.
       The following routine will create nprocs threads.
       ------------------------------------------------------------*/
    pcgstrf(&superlumt_options, &AC, perm_r, L, U, &Gstat, info);

    flopcnt = 0;
    for (i = 0; i < nprocs; ++i) flopcnt += Gstat.procstat[i].fcops;
    ops[FACT] = flopcnt;

#if ( PRNTlevel==1 )
    printf("nprocs = %d, flops %e, Mflops %.2f\n",
	   nprocs, flopcnt, flopcnt/utime[FACT]*1e-6);
    printf("Parameters: w %d, relax %d, maxsuper %d, rowblk %d, colblk %d\n",
	   sp_ienv(1), sp_ienv(2), sp_ienv(3), sp_ienv(4), sp_ienv(5));
    fflush(stdout);
#endif

    /* ------------------------------------------------------------
       Solve the system A*X=B, overwriting B with X.
       ------------------------------------------------------------*/
    if ( *info == 0 ) {
        t = SuperLU_timer_();
	cgstrs (trans, L, U, perm_r, perm_c, B, &Gstat, info);
	utime[SOLVE] = SuperLU_timer_() - t;
	ops[SOLVE] = ops[TRISOLVE];
    }

    /* ------------------------------------------------------------
       Deallocate storage after factorization.
       ------------------------------------------------------------*/
    pxgstrf_finalize(&superlumt_options, &AC);
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
    PrintStat(&Gstat);
    StatFree(&Gstat);
}
