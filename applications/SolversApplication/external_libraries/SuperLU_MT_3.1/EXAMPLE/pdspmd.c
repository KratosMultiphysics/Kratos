/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Last modified: 02/15/2013
 *                04/16/2015
 *
 * Purpose
 * =======
 *
 * This program illustrates how to integrate the SPMD mode
 * of the factor routine pdgstrf_thread() into a larger SPMD 
 * application code, and manage the threads yourself.
 * In this example, the threads creation happens only once.
 * 
 */
#include <stdlib.h> /* for getenv and atoi */
#include "slu_mt_ddefs.h"


/* Arguments passed to each dot product thread. */
typedef struct {
    int_t i;
    int_t len;
    int_t nprocs;
    double *global_dot;
    double *x;
} pddot_threadarg_t;

#if ( MACH==SUN )
mutex_t pddot_lock;
#elif ( MACH==DEC || MACH==PTHREAD )
pthread_mutex_t pddot_lock;
#elif ( MACH==SGI || MACH==ORIGIN || MACH==CRAY_PVP )
int_t pddot_lock = 1;
#endif

/* Arguments passed to each SPMD program. */
typedef struct {
    pdgstrf_threadarg_t *pdgstrf_threadarg; /* for sparse LU. */
    pddot_threadarg_t *pddot_threadarg; /* for dot product. */
} dspmd_arg_t;

void *dspmd(void *);
void *pddot_thread(void *);

main(int argc, char *argv[])
{
    SuperMatrix A, AC, L, U, B;
    NCformat    *Astore;
    SCPformat   *Lstore;
    NCPformat   *Ustore;
    superlumt_options_t superlumt_options;
    pxgstrf_shared_t pxgstrf_shared;
    pdgstrf_threadarg_t *pdgstrf_threadarg;
    pddot_threadarg_t *pddot_threadarg;
    dspmd_arg_t *dspmd_arg;
    int_t         nprocs;
    fact_t      fact;
    trans_t     trans;
    yes_no_t    refact, usepr;
    double      u, drop_tol;
    double      *a;
    int_t         *asub, *xa;
    int_t         *perm_c; /* column permutation vector */
    int_t         *perm_r; /* row permutations from partial pivoting */
    void        *work;
    int_t         info, lwork, nrhs, ldx; 
    int_t         m, n, nnz, permc_spec, panel_size, relax;
    int_t         i, firstfact;
    double      *rhsb, *xact;
    Gstat_t Gstat;
    int_t    vlength;
    double *xvector, xdot, temp;
    double zero = 0.0;
#if ( MACH==SUN )
    thread_t  *thread_id;
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_t *thread_id;
    void      *status;
#endif
    void parse_command_line();

    /* Default parameters to control factorization. */
    nprocs = 1;
    fact  = EQUILIBRATE;
    trans = NOTRANS;
    refact= NO;
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    u     = 1.0;
    usepr = NO;
    drop_tol = 0.0;
    work = NULL;
    lwork = 0;
    nrhs  = 1;

    /* Get the number of processes from command line. */
    parse_command_line(argc, argv, &nprocs);

    /* Read the input matrix stored in Harwell-Boeing format. */
    dreadhb(&m, &n, &nnz, &a, &asub, &xa);

    /* Set up the sparse matrix data structure for A. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    if (!(rhsb = doubleMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsb[].");
    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(trans, nrhs, xact, ldx, &A, &B);
    
    if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
    if ( !(superlumt_options.etree = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for etree[].");
    if ( !(superlumt_options.colcnt_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    if ( !(superlumt_options.part_super_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");


    /* ------------------------------------------------------------
       Get column permutation vector perm_c[], according to permc_spec:
       permc_spec = 0: natural ordering 
       permc_spec = 1: minimum degree ordering on structure of A'*A
       permc_spec = 2: minimum degree ordering on structure of A'+A
       permc_spec = 3: approximate minimum degree for unsymmetric matrices
       ------------------------------------------------------------*/ 	
    permc_spec = 1;
    get_perm_c(permc_spec, &A, perm_c);


    /* ------------------------------------------------------------
       Allocate storage and initialize statistics variables. 
       ------------------------------------------------------------*/
    StatAlloc(n, nprocs, panel_size, relax, &Gstat);
    StatInit(n, nprocs, &Gstat);

    /* ------------------------------------------------------------
       Initialize the option structure superlumt_options using the
       user-input parameters;
       Apply perm_c to the columns of original A to form AC.
       ------------------------------------------------------------*/
    pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
		 u, usepr, drop_tol, perm_c, perm_r,
		 work, lwork, &A, &AC, &superlumt_options, &Gstat);

    /* --------------------------------------------------------------
       Initializes the parallel data structures for pdgstrf_thread().
       --------------------------------------------------------------*/
    pdgstrf_threadarg = pdgstrf_thread_init(&AC, &L, &U, &superlumt_options,
					    &pxgstrf_shared, &Gstat, &info);
    
    /* ------------------------------------------------------------
       Initialization for inner product routine.
       ------------------------------------------------------------*/
    vlength = 1000;
    xdot = zero;
    xvector = doubleMalloc(vlength);
    for (i = 0; i < vlength; ++i){
        xvector[i] = i+1;
    }
    pddot_threadarg = (pddot_threadarg_t *) 
        SUPERLU_MALLOC(nprocs * sizeof(pddot_threadarg_t));
#if ( MACH==SUN )
    mutex_init(&pddot_lock, USYNC_THREAD, 0);
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_init(&pddot_lock, NULL);
#endif
    for (i = 0; i < nprocs; ++i) {
        pddot_threadarg[i].i = i;
        pddot_threadarg[i].len = vlength;
	pddot_threadarg[i].nprocs = nprocs;
        pddot_threadarg[i].global_dot = &xdot;
	pddot_threadarg[i].x = xvector;
    }

    /* ------------------------------------------------------------
       Create multiple threads
       ------------------------------------------------------------*/
       
#if ( MACH==SUN ) /* Use Solaris thread ... */
    
    /* Prepare arguments to each SPMD program. */
    dspmd_arg = (dspmd_arg_t *) 
        SUPERLU_MALLOC(nprocs * sizeof(dspmd_arg_t));
    for (i = 0; i < nprocs; ++i) {
	dspmd_arg[i].pdgstrf_threadarg = &pdgstrf_threadarg[i];
	dspmd_arg[i].pddot_threadarg = &pddot_threadarg[i];
    }

    thread_id = (thread_t *) SUPERLU_MALLOC(nprocs * sizeof(thread_t));
    
    /* Run the SPMD programs in parallel. */
    for (i = 0; i < nprocs; ++i)
	thr_create(NULL, 0, dspmd, &dspmd_arg[i], 0, &thread_id[i]);
	 
    /* Wait for all threads to terminate. */
    for (i = 0; i < nprocs; i++) pthread_join(thread_id[i], 0, 0);

    SUPERLU_FREE (dspmd_arg);
    SUPERLU_FREE (thread_id);

#elif ( MACH==DEC || MACH==PTHREAD ) /* Use pthread ... */
    
    /* Prepare arguments to each SPMD program. */
    dspmd_arg = (dspmd_arg_t *) 
        SUPERLU_MALLOC(nprocs * sizeof(dspmd_arg_t));
    for (i = 0; i < nprocs; ++i) {
	dspmd_arg[i].pdgstrf_threadarg = &pdgstrf_threadarg[i];
	dspmd_arg[i].pddot_threadarg = &pddot_threadarg[i];
    }

    thread_id = (pthread_t *) SUPERLU_MALLOC(nprocs * sizeof(pthread_t));
    
    /* Run the SPMD programs in parallel. */
    for (i = 0; i < nprocs; ++i)
	pthread_create(&thread_id[i], NULL,
		       dspmd, &dspmd_arg[i]);
	 
    /* Wait for all threads to terminate. */
    for (i = 0; i < nprocs; i++) pthread_join(thread_id[i], &status);

    SUPERLU_FREE (dspmd_arg);
    SUPERLU_FREE (thread_id);

/* _PTHREAD */

#elif ( MACH==SGI || MACH==ORIGIN ) /* Use parallel C ... */
    if ( getenv("MP_SET_NUMTHREADS") ) {
        i = atoi(getenv("MP_SET_NUMTHREADS"));
	if ( nprocs > i ) {
	    printf("nprocs=" IFMT " > environment allowed: MP_SET_NUMTHREADS=" IFMT "\n",
		   nprocs, i);
	    exit(-1);
	}
    }
#pragma parallel
#pragma shared (pdgstrf_threadarg, pddot_threadarg)
#pragma numthreads (nprocs)
    {
        pdgstrf_thread( pdgstrf_threadarg );
	pddot_thread( pddot_threadarg );
    }
/* _SGI or _ORIGIN */

#elif ( MACH==CRAY_PVP ) /* Use C microtasking ... */
    if ( getenv("NCPUS") ) {
        i = atoi(getenv("NCPUS"));
	if ( nprocs > i ) {
	    printf("nprocs=" IFMT " > environment allowed: NCPUS=" IFMT "\n",
		   nprocs, i);
	    exit(-1);
	}
    }
#pragma _CRI taskloop private(i,nprocs) \
                      shared(pdgstrf_threadarg, pddot_threadarg)
    /* stand-alone task loop */
    for (i = 0; i < nprocs; ++i) {
	pdgstrf_thread( &pdgstrf_threadarg[i] );
	pddot_thread( &pddot_threadarg[i] );
    }
/* _CRAY_PVP */

#endif /* _PTHREAD */

 
    /* ------------------------------------------------------------
       Clean up and free storage after multithreaded factorization.
       ------------------------------------------------------------*/
    pdgstrf_thread_finalize(pdgstrf_threadarg, &pxgstrf_shared, 
			    &A, perm_r, &L, &U);

     /* ------------------------------------------------------------
       Deallocate storage after factorization.
       ------------------------------------------------------------*/
    pxgstrf_finalize(&superlumt_options, &AC);


    /* ------------------------------------------------------------
       Solve the system A*X=B, overwriting B with X.
       ------------------------------------------------------------*/
    dgstrs(trans, &L, &U, perm_r, perm_c, &B, &Gstat, &info);
    
    printf("\n** Result of sparse LU **\n");
    dinf_norm_error(nrhs, &B, xact); /* Check inf. norm of the error */

    Lstore = (SCPformat *) L.Store;
    Ustore = (NCPformat *) U.Store;
    printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
    printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
    printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
    fflush(stdout);

/* DAN FIX THIS */

    printf("\n** Result of dot product **\n");
    temp = vlength*(vlength+1)*(2*vlength+1)/6;

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    if ( lwork == 0 ) {
        Destroy_SuperNode_SCP(&L);
        Destroy_CompCol_NCP(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }
    StatFree(&Gstat);

    SUPERLU_FREE (pddot_threadarg);
    SUPERLU_FREE (xvector);

#if ( MACH==SUN )
    mutex_destroy( &pddot_lock );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_destroy( &pddot_lock );
#endif
}


void
*dspmd(void *arg)
{
    pdgstrf_threadarg_t *pdgstrf_threadarg;
    pddot_threadarg_t *pddot_threadarg;
    
    pdgstrf_threadarg = ((dspmd_arg_t *)arg)->pdgstrf_threadarg;
    pddot_threadarg = ((dspmd_arg_t *)arg)->pddot_threadarg;

    /* Perform LU factorization. */
    pdgstrf_thread(pdgstrf_threadarg);

    /* Perform local inner product, and add the sum into the global one. */
    pddot_thread(pddot_threadarg);

    return 0;
}


void
*pddot_thread(void *arg)
{
#if ( MACH==SGI || MACH==ORIGIN )
#if ( MACH==SGI )
    int_t         i = mpc_my_threadnum();
#else
    int_t         i = mp_my_threadnum();
#endif
    pddot_threadarg_t *pddot_arg = &((pddot_threadarg_t *)arg)[i];
#else
    pddot_threadarg_t *pddot_arg = arg;
    int_t i = pddot_arg->i;
#endif
    int_t len = pddot_arg->len;
    int_t nprocs = pddot_arg->nprocs;
    double *global_dot = pddot_arg->global_dot;
    double *x = pddot_arg->x;
    int_t chunk, start, end;
    double temp;
    double zero = 0.0;

    chunk = len / nprocs;
    start = i * chunk;
    end = start + chunk;
    if ( i == nprocs-1 ) end = len;
    printf("(" IFMT ") nprocs " IFMT ",  chunk " IFMT ", start " IFMT ", end " IFMT "\n", 
	   i, nprocs, chunk, start, end);

    temp = zero;
    for (i = start; i < end; ++i){
        temp += x[i]*x[i];
    }
   
#if ( MACH==SUN )
    mutex_lock( &pddot_lock );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_lock( &pddot_lock );
#elif ( MACH==SGI || MACH==ORIGIN )
#pragma critical lock( pddot_lock );
#elif ( MACH==CRAY_PVP )
#pragma _CRI guard pddot_lock
#endif
    {
	*global_dot += temp;
    }
#if ( MACH==SUN )
    mutex_unlock( &pddot_lock );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_unlock( &pddot_lock );
#elif ( MACH==CRAY_PVP )
#pragma _CRI endguard pddot_lock
#endif

    return 0;
}


/*  
 * Parse command line to get nprocs, the number of processes.
 */
void
parse_command_line(int argc, char *argv[], int_t *nprocs)
{
    register int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "hp:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options: (default values are in parenthesis)\n");
	    printf("\t-p <int> - number of processes     ( " IFMT " )\n", *nprocs);
	    exit(1);
	    break;
	  case 'p': *nprocs = atoi(optarg); 
	            break;
  	}
    }
}

