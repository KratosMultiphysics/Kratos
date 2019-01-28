/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include <stdlib.h> /* for getenv and atoi */
#include "slu_mt_ddefs.h"

void
pdgstrf(superlumt_options_t *superlumt_options, SuperMatrix *A, int_t *perm_r,
	SuperMatrix *L, SuperMatrix *U,	Gstat_t *Gstat, int_t *info)
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
 * PDGSTRF computes an LU factorization of a general sparse nrow-by-ncol
 * matrix A using partial pivoting with row interchanges. The factorization
 * has the form
 *     Pr * A = L * U
 * where Pr is a row permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapezoidal if A->nrow > A->ncol), and U is
 * upper triangular (upper trapezoidal if A->nrow < A->ncol).
 *
 * Arguments
 * =========
 * 
 * superlumt_options (input) superlumt_options_t*
 *        The structure defines the parameters to control how the sparse
 *        LU factorization is performed. The following fields must be set 
 *        by the user:
 *
 *        o nprocs (int_t)
 *          Number of processes to be spawned and used for factorization.
 *
 *        o refact (yes_no_t)
 *          Specifies whether this is first time or subsequent factorization.
 *          = NO:  this factorization is treated as the first one;
 *          = YES: it means that a factorization was performed prior to this
 *                 one. Therefore, this factorization will re-use some
 *                 existing data structures, such as L and U storage, column
 *                 elimination tree, and the symbolic information of the
 *                 Householder matrix.
 *
 *        o panel_size (int_t)
 *          A panel consists of at most panel_size consecutive columns.
 *
 *        o relax (int_t)
 *          Degree of relaxing supernodes. If the number of nodes (columns)
 *          in a subtree of the elimination tree is less than relax, this 
 *          subtree is considered as one supernode, regardless of the row
 *          structures of those columns.
 *
 *        o diag_pivot_thresh (double)
 *	    Diagonal pivoting threshold. At step j of Gaussian elimination,
 *          if abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
 *          use A_jj as pivot. 0 <= diag_pivot_thresh <= 1. The default
 *          value is 1.0, corresponding to partial pivoting.
 *
 *        o usepr (yes_no_t)
 *          Whether the pivoting will use perm_r specified by the user.
 *          = YES: use perm_r; perm_r is input, unchanged on exit.
 *          = NO:  perm_r is determined by partial pivoting, and is output.
 *
 *        o drop_tol (double) (NOT IMPLEMENTED)
 *	    Drop tolerance parameter. At step j of the Gaussian elimination,
 *          if abs(A_ij)/(max_i abs(A_ij)) < drop_tol, drop entry A_ij.
 *          0 <= drop_tol <= 1. The default value of drop_tol is 0,
 *          corresponding to not dropping any entry.
 *
 *        o perm_c (int_t*)
 *	    Column permutation vector of size A->ncol, which defines the 
 *          permutation matrix Pc; perm_c[i] = j means column i of A is 
 *          in position j in A*Pc.
 *
 *        o perm_r (int_t*)
 *	    Column permutation vector of size A->nrow.
 *          If superlumt_options->usepr = NO, this is an output argument.
 *
 *        o work (void*) of size lwork
 *          User-supplied work space and space for the output data structures.
 *          Not referenced if lwork = 0;
 *
 *        o lwork (int_t)
 *          Specifies the length of work array.
 *            = 0:  allocate space internally by system malloc;
 *            > 0:  use user-supplied work array of length lwork in bytes,
 *                  returns error if space runs out.
 *            = -1: the routine guesses the amount of space needed without
 *                  performing the factorization, and returns it in
 *                  superlu_memusage->total_needed; no other side effects.
 *
 * A      (input) SuperMatrix*
 *	  Original matrix A, permuted by columns, of dimension
 *        (A->nrow, A->ncol). The type of A can be:
 *        Stype = NCP; Dtype = _D; Mtype = GE.
 *
 * perm_r (input/output) int_t*, dimension A->nrow
 *        Row permutation vector which defines the permutation matrix Pr,
 *        perm_r[i] = j means row i of A is in position j in Pr*A.
 *        If superlumt_options->usepr = NO, perm_r is output argument;
 *        If superlumt_options->usepr = YES, the pivoting routine will try 
 *           to use the input perm_r, unless a certain threshold criterion
 *           is violated. In that case, perm_r is overwritten by a new
 *           permutation determined by partial pivoting or diagonal 
 *           threshold pivoting.
 *
 * L      (output) SuperMatrix*
 *        The factor L from the factorization Pr*A=L*U; use compressed row 
 *        subscripts storage for supernodes, i.e., L has type: 
 *        Stype = SCP, Dtype = _D, Mtype = TRLU.
 *
 * U      (output) SuperMatrix*
 *	  The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
 *        storage scheme, i.e., U has types: Stype = NCP, Dtype = _D,
 *        Mtype = TRU.
 *
 * Gstat  (output) Gstat_t*
 *        Record all the statistics about the factorization; 
 *        See Gstat_t structure defined in slu_mt_util.h.
 *
 * info   (output) int_t*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                and division by zero will occur if it is used to solve a
 *                system of equations.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 *
 */
    pdgstrf_threadarg_t *pdgstrf_threadarg;
    pxgstrf_shared_t pxgstrf_shared;
    register int_t nprocs = superlumt_options->nprocs;
    register int_t i, iinfo;
    double    *utime = Gstat->utime;
    double    usrtime, wtime;
    double    usertimer_();
#if ( MACH==SUN )
    thread_t  *thread_id;
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_t *thread_id;
    void      *status;
#endif
    void      *pdgstrf_thread(void *);


    /* --------------------------------------------------------------
       Initializes the parallel data structures for pdgstrf_thread().
       --------------------------------------------------------------*/
    pdgstrf_threadarg = pdgstrf_thread_init(A, L, U, superlumt_options,
					    &pxgstrf_shared, Gstat, info);
    if ( *info ) return;

    /* Start timing factorization. */
    usrtime = usertimer_();
    wtime = SuperLU_timer_(); 

    /* ------------------------------------------------------------
       On a SUN multiprocessor system, use Solaris thread.
       ------------------------------------------------------------*/
#if ( MACH==SUN )
    
    /* Create nproc threads for concurrent factorization. */
    thread_id = (thread_t *) SUPERLU_MALLOC(nprocs * sizeof(thread_t));
    
    for (i = 1; i < nprocs; ++i) {
#if ( PRNTlevel==1 )
	printf(".. Create unbound threads: i " IFMT ", nprocs " IFMT "\n",
               i, nprocs);
#endif
	if ( (iinfo = thr_create(NULL, 0,
			       pdgstrf_thread, &(pdgstrf_threadarg[i]),
			       0, &thread_id[i])) )
	{
	    fprintf(stderr, "thr_create: %d\n", iinfo);
	    SUPERLU_ABORT("thr_creat()");
	}
    }
	 
    pdgstrf_thread( &(pdgstrf_threadarg[0]) );
    
    /* Wait for all threads to terminate. */
    for (i = 1; i < nprocs; i++) thr_join(thread_id[i], 0, 0);
    SUPERLU_FREE (thread_id);
/* _SOLARIS_2 */


    /* ------------------------------------------------------------
       On a DEC multiprocessor system, use pthread.
       ------------------------------------------------------------*/
#elif ( MACH==DEC )	/* Use DECthreads ... */
    
    /* Create nproc threads for concurrent factorization. */
    thread_id = (pthread_t *) SUPERLU_MALLOC(nprocs * sizeof(pthread_t));
    
    for (i = 0; i < nprocs; ++i) {
	if ( iinfo = pthread_create(&thread_id[i],
				    NULL,
				    pdgstrf_thread, 
				    &(pdgstrf_threadarg[i])) ) {
	    fprintf(stderr, "pthread_create: %d\n", iinfo);
	    SUPERLU_ABORT("pthread_create()");
	}
	/*	pthread_bind_to_cpu_np(thread_id[i], i);*/
    }
	 
    /* Wait for all threads to terminate. */
    for (i = 0; i < nprocs; i++)
	pthread_join(thread_id[i], &status);
    SUPERLU_FREE (thread_id);
/* _DEC */

    
    /* ------------------------------------------------------------
       On a SGI Power Challenge or Origin multiprocessor system,
       use parallel C.
       ------------------------------------------------------------*/
#elif ( MACH==SGI || MACH==ORIGIN )	/* Use parallel C ... */

    if ( getenv("MP_SET_NUMTHREADS") ) {
        i = atoi(getenv("MP_SET_NUMTHREADS"));
	if ( nprocs > i ) {
	    printf("nprocs=" IFMT "> environment allowed: MP_SET_NUMTHREADS=" IFMT "\n",
		   nprocs, i);
	    exit(-1);
	}
    }
#pragma parallel
#pragma shared  (pdgstrf_threadarg)
/*#pragma numthreads (max = nprocs)*/
#pragma numthreads (nprocs)
    {
      pdgstrf_thread( pdgstrf_threadarg );
    }
/* _SGI or _ORIGIN */

    /* ------------------------------------------------------------
       On a Cray PVP multiprocessor system, use microtasking.
       ------------------------------------------------------------*/
#elif ( MACH==CRAY_PVP )       /* Use C microtasking. */

    if ( getenv("NCPUS") ) {
        i = atoi(getenv("NCPUS"));
	if ( nprocs > i ) {
	    printf("nprocs=" IFMT "> environment allowed: NCPUS=" IFMT "\n",
		   nprocs, i);
	    exit(-1);
	}
    }
#pragma _CRI taskloop private (i,nprocs)  shared (pdgstrf_threadarg)
    /* Stand-alone task loop */
    for (i = 0; i < nprocs; ++i) {
	pdgstrf_thread( &(pdgstrf_threadarg[i]) );
    }
/* _CRAY_PVP */

    /* ------------------------------------------------------------
       Use POSIX threads.
       ------------------------------------------------------------*/
#elif ( MACH==PTHREAD )	/* Use pthread ... */
    
    /* Create nproc threads for concurrent factorization. */
    thread_id = (pthread_t *) SUPERLU_MALLOC(nprocs * sizeof(pthread_t));
    
    for (i = 0; i < nprocs; ++i) {
	if ( iinfo = pthread_create(&thread_id[i],
				    NULL,
				    pdgstrf_thread, 
				    &(pdgstrf_threadarg[i])) ) {
	    fprintf(stderr, "pthread_create: " IFMT "\n", iinfo);
	    SUPERLU_ABORT("pthread_create()");
	}
    }
	 
    /* Wait for all threads to terminate. */
    for (i = 0; i < nprocs; i++)
	pthread_join(thread_id[i], &status);
    SUPERLU_FREE (thread_id);
/* _PTHREAD */

    /* ------------------------------------------------------------
       Use openMP.
       ------------------------------------------------------------*/
#elif ( MACH==OPENMP ) /* Use openMP ... */

#pragma omp parallel for shared (pdgstrf_threadarg) private (i)
    /* Stand-alone task loop */
    for (i = 0; i < nprocs; ++i) {
        pdgstrf_thread( &(pdgstrf_threadarg[i]) );
    }

/* _OPENMP */

    /* ------------------------------------------------------------
       On all other systems, use single processor.
       ------------------------------------------------------------*/
#else

    printf("pdgstrf() is not parallelized on this machine.\n");
    printf("pdgstrf() will be run on single processor.\n");
    pdgstrf_thread( &(pdgstrf_threadarg[0]) );
    
#endif    

    wtime = SuperLU_timer_() - wtime;
    usrtime = usertimer_() - usrtime;
    utime[FACT] = wtime;

#if ( PRNTlevel==1 )
    printf(".. pdgstrf_thread() returns info " IFMT ", usrtime %.2f, wtime %.2f\n",
           *info, usrtime, wtime);
#endif

    /* check_mem_leak("after pdgstrf_thread()"); */
    
    /* ------------------------------------------------------------
       Clean up and free storage after multithreaded factorization.
       ------------------------------------------------------------*/
    pdgstrf_thread_finalize(pdgstrf_threadarg, &pxgstrf_shared, 
			    A, perm_r, L, U);

}
