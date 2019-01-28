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
 * of the factor routine pzgstrf_thread() into a larger SPMD 
 * application code, and manage the threads yourself.
 * In this example, the threads creation happens only once.
 * 
 */
#include <stdlib.h> /* for getenv and atoi */
#include "slu_mt_zdefs.h"


/* Arguments passed to each dot product thread. */
typedef struct {
    int_t i;
    int_t len;
    int_t nprocs;
    doublecomplex *global_dot;
    doublecomplex *x;
} pzdot_threadarg_t;

#if ( MACH==SUN )
mutex_t pzdot_lock;
#elif ( MACH==DEC || MACH==PTHREAD )
pthread_mutex_t pzdot_lock;
#elif ( MACH==SGI || MACH==ORIGIN || MACH==CRAY_PVP )
int_t pzdot_lock = 1;
#endif

/* Arguments passed to each SPMD program. */
typedef struct {
    pzgstrf_threadarg_t *pzgstrf_threadarg; /* for sparse LU. */
    pzdot_threadarg_t *pzdot_threadarg; /* for dot product. */
} zspmd_arg_t;

void *zspmd(void *);
void *pzdot_thread(void *);

main(int argc, char *argv[])
{
    SuperMatrix A, AC, L, U, B;
    NCformat    *Astore;
    SCPformat   *Lstore;
    NCPformat   *Ustore;
    superlumt_options_t superlumt_options;
    pxgstrf_shared_t pxgstrf_shared;
    pzgstrf_threadarg_t *pzgstrf_threadarg;
    pzdot_threadarg_t *pzdot_threadarg;
    zspmd_arg_t *zspmd_arg;
    int_t         nprocs;
    fact_t      fact;
    trans_t     trans;
    yes_no_t    refact, usepr;
    double      u, drop_tol;
    doublecomplex      *a;
    int_t         *asub, *xa;
    int_t         *perm_c; /* column permutation vector */
    int_t         *perm_r; /* row permutations from partial pivoting */
    void        *work;
    int_t         info, lwork, nrhs, ldx; 
    int_t         m, n, nnz, permc_spec, panel_size, relax;
    int_t         i, firstfact;
    doublecomplex      *rhsb, *xact;
    Gstat_t Gstat;
    int_t    vlength;
    doublecomplex *xvector, xdot, temp;
    doublecomplex zero = {0.0, 0.0};
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
    zreadhb(&m, &n, &nnz, &a, &asub, &xa);

    /* Set up the sparse matrix data structure for A. */
    zCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_Z, SLU_GE);

    if (!(rhsb = doublecomplexMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsb[].");
    zCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_Z, SLU_GE);
    xact = doublecomplexMalloc(n * nrhs);
    ldx = n;
    zGenXtrue(n, nrhs, xact, ldx);
    zFillRHS(trans, nrhs, xact, ldx, &A, &B);
    
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
    pzgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
		 u, usepr, drop_tol, perm_c, perm_r,
		 work, lwork, &A, &AC, &superlumt_options, &Gstat);

    /* --------------------------------------------------------------
       Initializes the parallel data structures for pzgstrf_thread().
       --------------------------------------------------------------*/
    pzgstrf_threadarg = pzgstrf_thread_init(&AC, &L, &U, &superlumt_options,
					    &pxgstrf_shared, &Gstat, &info);
    
    /* ------------------------------------------------------------
       Initialization for inner product routine.
       ------------------------------------------------------------*/
    vlength = 1000;
    xdot = zero;
    xvector = doublecomplexMalloc(vlength);
    for (i = 0; i < vlength; ++i){
        xvector[i].r = i+1;
        xvector[i].i = 0;
    }
    pzdot_threadarg = (pzdot_threadarg_t *) 
        SUPERLU_MALLOC(nprocs * sizeof(pzdot_threadarg_t));
#if ( MACH==SUN )
    mutex_init(&pzdot_lock, USYNC_THREAD, 0);
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_init(&pzdot_lock, NULL);
#endif
    for (i = 0; i < nprocs; ++i) {
        pzdot_threadarg[i].i = i;
        pzdot_threadarg[i].len = vlength;
	pzdot_threadarg[i].nprocs = nprocs;
        pzdot_threadarg[i].global_dot = &xdot;
	pzdot_threadarg[i].x = xvector;
    }

    /* ------------------------------------------------------------
       Create multiple threads
       ------------------------------------------------------------*/
       
#if ( MACH==SUN ) /* Use Solaris thread ... */
    
    /* Prepare arguments to each SPMD program. */
    zspmd_arg = (zspmd_arg_t *) 
        SUPERLU_MALLOC(nprocs * sizeof(zspmd_arg_t));
    for (i = 0; i < nprocs; ++i) {
	zspmd_arg[i].pzgstrf_threadarg = &pzgstrf_threadarg[i];
	zspmd_arg[i].pzdot_threadarg = &pzdot_threadarg[i];
    }

    thread_id = (thread_t *) SUPERLU_MALLOC(nprocs * sizeof(thread_t));
    
    /* Run the SPMD programs in parallel. */
    for (i = 0; i < nprocs; ++i)
	thr_create(NULL, 0, zspmd, &zspmd_arg[i], 0, &thread_id[i]);
	 
    /* Wait for all threads to terminate. */
    for (i = 0; i < nprocs; i++) pthread_join(thread_id[i], 0, 0);

    SUPERLU_FREE (zspmd_arg);
    SUPERLU_FREE (thread_id);

#elif ( MACH==DEC || MACH==PTHREAD ) /* Use pthread ... */
    
    /* Prepare arguments to each SPMD program. */
    zspmd_arg = (zspmd_arg_t *) 
        SUPERLU_MALLOC(nprocs * sizeof(zspmd_arg_t));
    for (i = 0; i < nprocs; ++i) {
	zspmd_arg[i].pzgstrf_threadarg = &pzgstrf_threadarg[i];
	zspmd_arg[i].pzdot_threadarg = &pzdot_threadarg[i];
    }

    thread_id = (pthread_t *) SUPERLU_MALLOC(nprocs * sizeof(pthread_t));
    
    /* Run the SPMD programs in parallel. */
    for (i = 0; i < nprocs; ++i)
	pthread_create(&thread_id[i], NULL,
		       zspmd, &zspmd_arg[i]);
	 
    /* Wait for all threads to terminate. */
    for (i = 0; i < nprocs; i++) pthread_join(thread_id[i], &status);

    SUPERLU_FREE (zspmd_arg);
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
#pragma shared (pzgstrf_threadarg, pzdot_threadarg)
#pragma numthreads (nprocs)
    {
        pzgstrf_thread( pzgstrf_threadarg );
	pzdot_thread( pzdot_threadarg );
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
                      shared(pzgstrf_threadarg, pzdot_threadarg)
    /* stand-alone task loop */
    for (i = 0; i < nprocs; ++i) {
	pzgstrf_thread( &pzgstrf_threadarg[i] );
	pzdot_thread( &pzdot_threadarg[i] );
    }
/* _CRAY_PVP */

#endif /* _PTHREAD */

 
    /* ------------------------------------------------------------
       Clean up and free storage after multithreaded factorization.
       ------------------------------------------------------------*/
    pzgstrf_thread_finalize(pzgstrf_threadarg, &pxgstrf_shared, 
			    &A, perm_r, &L, &U);

     /* ------------------------------------------------------------
       Deallocate storage after factorization.
       ------------------------------------------------------------*/
    pxgstrf_finalize(&superlumt_options, &AC);


    /* ------------------------------------------------------------
       Solve the system A*X=B, overwriting B with X.
       ------------------------------------------------------------*/
    zgstrs(trans, &L, &U, perm_r, perm_c, &B, &Gstat, &info);
    
    printf("\n** Result of sparse LU **\n");
    zinf_norm_error(nrhs, &B, xact); /* Check inf. norm of the error */

    Lstore = (SCPformat *) L.Store;
    Ustore = (NCPformat *) U.Store;
    printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
    printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
    printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
    fflush(stdout);

/* DAN FIX THIS */

    printf("\n** Result of dot product **\n");
    temp.r = vlength*(vlength+1)*(2*vlength+1)/6;
    temp.i = 0;
    printf("computed %e %e, correct %e\n", xdot.r, xdot.i, temp.r);

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

    SUPERLU_FREE (pzdot_threadarg);
    SUPERLU_FREE (xvector);

#if ( MACH==SUN )
    mutex_destroy( &pzdot_lock );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_destroy( &pzdot_lock );
#endif
}


void
*zspmd(void *arg)
{
    pzgstrf_threadarg_t *pzgstrf_threadarg;
    pzdot_threadarg_t *pzdot_threadarg;
    
    pzgstrf_threadarg = ((zspmd_arg_t *)arg)->pzgstrf_threadarg;
    pzdot_threadarg = ((zspmd_arg_t *)arg)->pzdot_threadarg;

    /* Perform LU factorization. */
    pzgstrf_thread(pzgstrf_threadarg);

    /* Perform local inner product, and add the sum into the global one. */
    pzdot_thread(pzdot_threadarg);

    return 0;
}


void
*pzdot_thread(void *arg)
{
#if ( MACH==SGI || MACH==ORIGIN )
#if ( MACH==SGI )
    int_t         i = mpc_my_threadnum();
#else
    int_t         i = mp_my_threadnum();
#endif
    pzdot_threadarg_t *pzdot_arg = &((pzdot_threadarg_t *)arg)[i];
#else
    pzdot_threadarg_t *pzdot_arg = arg;
    int_t i = pzdot_arg->i;
#endif
    int_t len = pzdot_arg->len;
    int_t nprocs = pzdot_arg->nprocs;
    doublecomplex *global_dot = pzdot_arg->global_dot;
    doublecomplex *x = pzdot_arg->x;
    int_t chunk, start, end;
    doublecomplex temp;
    doublecomplex zero = {0.0, 0.0};
    doublecomplex temp_comp;

    chunk = len / nprocs;
    start = i * chunk;
    end = start + chunk;
    if ( i == nprocs-1 ) end = len;
    printf("(" IFMT ") nprocs " IFMT ",  chunk " IFMT ", start " IFMT ", end " IFMT "\n", 
	   i, nprocs, chunk, start, end);

    temp = zero;
    for (i = start; i < end; ++i){
        zz_mult(&temp_comp, &x[i], &x[i]);
        z_add(&temp, &temp, &temp_comp);
    }
   
#if ( MACH==SUN )
    mutex_lock( &pzdot_lock );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_lock( &pzdot_lock );
#elif ( MACH==SGI || MACH==ORIGIN )
#pragma critical lock( pzdot_lock );
#elif ( MACH==CRAY_PVP )
#pragma _CRI guard pzdot_lock
#endif
    {
        z_add(global_dot, global_dot, &temp);
    }
#if ( MACH==SUN )
    mutex_unlock( &pzdot_lock );
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_unlock( &pzdot_lock );
#elif ( MACH==CRAY_PVP )
#pragma _CRI endguard pzdot_lock
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

