/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_zdefs.h"


void
pzgstrf_init(int_t nprocs, fact_t fact, trans_t trans, yes_no_t refact,
             int_t panel_size, int_t relax,
	     double diag_pivot_thresh, yes_no_t usepr, double drop_tol,
	     int_t *perm_c, int_t *perm_r, void *work, int_t lwork,
	     SuperMatrix *A, SuperMatrix *AC, 
	     superlumt_options_t *superlumt_options, Gstat_t *Gstat)
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
 * pzgstrf_init() initializes the option structure superlumt_options, using 
 * the user-input parameters. These options control how the factorization
 * will be performed by routine pzgstf().
 * 
 * In addition, it calls sp_colorder() to compute a postordered etree[], 
 * colcnt_h[] and super_part_h, and permute the columns of A using the 
 * permutation vector perm_c[]. See sp_colorder.c for details.
 *
 * Arguments
 * =========
 *
 * nprocs (input) int_t
 *        Number of processes used to perform LU factorization by pzgstrf().
 *
 * fact   (input) fact_t
 *        Specifies whether or not the factored form of the matrix is supplied.
 *
 * trans  (input) trans_t
 *        Specifies the form of the system of equations:
 *        = NOTRANS: A * X = B        (No transpose)
 *        = TRANS:   A**T * X = B     (Transpose)
 *        = CONJ:    A**H * X = B     (Transpose)
 *
 * refact (input) yes_no_t
 *        Specifies whether we want to use perm_r from a previous factor.
 *
 * panel_size (input) int_t
 *        A panel consists of at most panel_size consecutive columns.
 *
 * relax  (input) int_t
 *        To control degree of relaxing supernodes. If the number
 *        of nodes (columns) in a subtree of the elimination tree is less
 *        than relax, this subtree is considered as one supernode,
 *        regardless of the row structures of those columns.
 *
 * diag_pivot_thresh (input) double
 *        Diagonal pivoting threshold. At step j of the Gaussian elimination,
 *        if abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
 *        use A_jj as pivot. 0 <= diag_pivot_thresh <= 1. The default
 *        value is 1, corresponding to partial pivoting.
 *
 * drop_tol (input) double (NOT IMPLEMENTED)
 *	  Drop tolerance parameter. At step j of the Gaussian elimination,
 *        if abs(A_ij)/(max_i abs(A_ij)) < drop_tol, drop entry A_ij.
 *        0 <= drop_tol <= 1. The default value of drop_tol is 0, 
 *        corresponding to not dropping any entry.
 *
 * perm_c (input) int_t*, dimension A->ncol
 *	  Column permutation vector, which defines the 
 *        permutation matrix Pc; perm_c[i] = j means column i of A is 
 *        in position j in A*Pc.
 *        When search for diagonal, perm_c[*] is applied to the
 *        row subscripts of A, so that diagonal threshold pivoting
 *        can find the diagonal of A, instead of that of A*Pc.
 *
 * perm_r (input/output) int_t*, dimension A->nrow
 *        Row permutation vector which defines the permutation matrix Pr,
 *        perm_r[i] = j means row i of A is in position j in Pr*A.
 *        If usepr = NO, perm_r is output argument;
 *        If usepr = YES, the pivoting routine will try to use the input
 *           perm_r, unless a certain threshold criterion is violated.
 *           In that case, perm_r is overwritten by a new permutation
 *           determined by partial pivoting or diagonal threshold pivoting.
 *
 * work   (input) void* of size lwork
 *        User-supplied work space and space for the output data structures.
 *        Not referenced if lwork = 0;
 *
 * lwork  (input) int_t
 *        Specifies the length of work array.
 *        = 0:  allocate space internally by system malloc;
 *        > 0:  use user-supplied work array of length lwork in bytes,
 *              returns error if space runs out.
 *        = -1: the routine guesses the amount of space needed without
 *              performing the factorization, and returns it in
 *              superlu_memusage->total_needed; no other side effects.
 *
 * A      (input) SuperMatrix*
 *        Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *        of linear equations is A->nrow. Currently, the type of A can be:
 *        Stype = NC or NCP; Dtype = _D; Mtype = GE. In the future, more
 *        general A can be handled.
 *
 * AC     (output) SuperMatrix*
 *        The resulting matrix after applied the column permutation
 *        perm_c[] to matrix A. The type of AC can be:
 *        Stype = NCP; Dtype = _D; Mtype = GE.
 *
 * superlumt_options (output) superlumt_options_t*
 *        The structure defines the parameters to control how the sparse
 *        LU factorization is performed, and will be input to pzgstrf().
 *
 * Gstat  (output) Gstat_t*
 *        Record the time used in zsp_colorder phase.
 *
 */
    double t;

    superlumt_options->nprocs = nprocs;
    superlumt_options->refact = refact;
    superlumt_options->panel_size = panel_size;
    superlumt_options->relax = relax;
    superlumt_options->diag_pivot_thresh = diag_pivot_thresh;
    superlumt_options->usepr = usepr;
    superlumt_options->drop_tol = drop_tol;
    superlumt_options->SymmetricMode = NO;
    superlumt_options->PrintStat = NO;

    /* 
     * The following should be retained for repeated factorizations.
     */
    superlumt_options->perm_c = perm_c;
    superlumt_options->perm_r = perm_r;
    superlumt_options->work = work;
    superlumt_options->lwork = lwork;

    if ( refact == NO ) { /* First time factorization, need allocation. */
        if ( !(superlumt_options->etree = intMalloc(A->ncol)) )
	    SUPERLU_ABORT("Malloc fails for etree[].");
        if ( !(superlumt_options->colcnt_h = intMalloc(A->ncol)) )
	    SUPERLU_ABORT("Malloc fails for colcnt_h[].");
        if ( !(superlumt_options->part_super_h = intMalloc(A->ncol)) )
	    SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    }

    t = SuperLU_timer_();
    sp_colorder(A, perm_c, superlumt_options, AC);
    Gstat->utime[ETREE] = SuperLU_timer_() - t;

#if ( DEBUGlevel==1 )
    printf("** pzgstrf_init() called\n");
#endif
}

