/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_ddefs.h"

void
pdgstrf_thread_finalize(pdgstrf_threadarg_t *pdgstrf_threadarg, 
			pxgstrf_shared_t *pxgstrf_shared,
			SuperMatrix *A, int_t *perm_r,
			SuperMatrix *L, SuperMatrix *U
			)
{
/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 *
 * Purpose
 * =======
 * 
 * pdgstrf_thread_finalize() performs cleanups after the multithreaded 
 * factorization pdgstrf_thread(). It sets up the L and U data
 * structures, and deallocats the storage associated with the structures
 * pxgstrf_shared and pdgstrf_threadarg.
 *
 * Only Master thread calls this routine.
 *
 * Arguments
 * =========
 *
 * pdgstrf_threadarg (input) pdgstrf_threadarg_t*
 *          The structure contains the parameters to each thread.
 *
 * pxgstrf_shared (input) pxgstrf_shared_t*
 *          The structure contains the shared task queue, the 
 *          synchronization variables, and the L and U data structures.
 *
 * A        (input) SuperMatrix*
 *	    Original matrix A, permutated by columns, of dimension
 *          (A->nrow, A->ncol). The type of A can be:
 *          Stype = NCP; Dtype = _D; Mtype = GE.
 *
 * perm_r   (input) int_t*, dimension A->nrow
 *          Row permutation vector which defines the permutation matrix Pr,
 *          perm_r[i] = j means row i of A is in position j in Pr*A.
 *
 * L        (output) SuperMatrix*
 *          The factor L from the factorization Pr*A=L*U; use compressed row 
 *          subscripts storage for supernodes, i.e., L has type: 
 *          Stype = SCP, Dtype = _D, Mtype = TRLU.
 *
 * U        (output) SuperMatrix*
 *	    The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
 *          storage scheme, i.e., U has type:
 *          Stype = NCP, Dtype = _D, Mtype = TRU.
 *
 *
 */
    register int_t nprocs, n, i, iinfo;
    int_t       nnzL, nnzU;
    superlumt_options_t *superlumt_options;
    GlobalLU_t *Glu;
    extern ExpHeader *dexpanders;

    n = A->ncol;
    superlumt_options = pdgstrf_threadarg->superlumt_options;
    Glu = pxgstrf_shared->Glu;
    Glu->supno[n] = Glu->nsuper;

    countnz(n, pxgstrf_shared->xprune, &nnzL, &nnzU, Glu);
    fixupL(n, perm_r, Glu);
    
#ifdef COMPRESS_LUSUP
    compressSUP(n, pxgstrf_shared->Glu);
#endif

    if ( superlumt_options->refact == YES ) {
        /* L and U structures may have changed due to possibly different
	   pivoting, although the storage is available. */
        ((SCPformat *)L->Store)->nnz = nnzL;
	((SCPformat *)L->Store)->nsuper = Glu->supno[n];
	((NCPformat *)U->Store)->nnz = nnzU;
    } else {
	dCreate_SuperNode_Permuted(L, A->nrow, A->ncol, nnzL, Glu->lusup,
				   Glu->xlusup, Glu->xlusup_end, 
				   Glu->lsub, Glu->xlsub, Glu->xlsub_end,
				   Glu->supno, Glu->xsup, Glu->xsup_end,
				   SLU_SCP, SLU_D, SLU_TRLU);
	dCreate_CompCol_Permuted(U, A->nrow, A->ncol, nnzU, Glu->ucol,
				 Glu->usub, Glu->xusub, Glu->xusub_end,
				 SLU_NCP, SLU_D, SLU_TRU);
    }

    /* Combine the INFO returned from individual threads. */
    iinfo = 0;
    nprocs = superlumt_options->nprocs;
    for (i = 0; i < nprocs; ++i) {
        if ( pdgstrf_threadarg[i].info ) {
	    if (iinfo) iinfo=SUPERLU_MIN(iinfo, pdgstrf_threadarg[i].info);
	    else iinfo = pdgstrf_threadarg[i].info;
	}
    }
    *pxgstrf_shared->info = iinfo;

#if ( DEBUGlevel>=2 )
    printf("Last nsuper %d\n", Glu->nsuper);
    QueryQueue(&pxgstrf_shared->taskq);
    PrintGLGU(n, pxgstrf_shared->xprune, Glu);
    PrintInt10("perm_r", n, perm_r);
    PrintInt10("inv_perm_r", n, pxgstrf_shared->inv_perm_r);
#endif

    /* Deallocate the storage used by the parallel scheduling algorithm. */
    ParallelFinalize(pxgstrf_shared);
    SUPERLU_FREE(pdgstrf_threadarg);
    SUPERLU_FREE(pxgstrf_shared->inv_perm_r);
    SUPERLU_FREE(pxgstrf_shared->inv_perm_c);
    SUPERLU_FREE(pxgstrf_shared->xprune);
    SUPERLU_FREE(pxgstrf_shared->ispruned);
    SUPERLU_FREE(dexpanders);
    dexpanders = 0;

#if ( DEBUGlevel>=1 )
    printf("** pdgstrf_thread_finalize() called\n");
#endif
}
