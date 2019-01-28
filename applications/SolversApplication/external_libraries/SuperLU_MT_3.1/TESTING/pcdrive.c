/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include <string.h>
#include "slu_mt_cdefs.h"

#define NTESTS    5    /* Number of test types */
#define NTYPES    11   /* Number of matrix types */
#define NTRAN     2    
#define THRESH    40.0

#if defined ( _LONGINT )
#define FMT1      "%10s:n=%lld, test(%d)=%12.5g\n"
#define	FMT2      "%10s:fact=%d, trans=%d, refact=%d, equed=%d, n=%lld, imat=%d, test(%d)=%12.5g\n"
#define FMT3      "%10s:info=%lld, izero=%lld, n=%lld, nrhs=%lld, imat=%d, nfail=%d\n"
#else
#define FMT1      "%10s:n=%d, test(%d)=%12.5g\n"
#define	FMT2      "%10s:fact=%d, trans=%d, refact=%d, equed=%d, n=%d, imat=%d, test(%d)=%12.5g\n"
#define FMT3      "%10s:info=%d, izero=%d, n=%d, nrhs=%d, imat=%d, nfail=%d\n"
#endif


main(int argc, char *argv[])
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
 * pcdrive() is the main test program for the DOUBLE linear 
 * equation driver routines pcgssv() and pcgssvx().
 * 
 * The program is invoked by a shell script file -- pctest.csh.
 * The output from the tests are written into a file -- pctest.out.
 *
 * =====================================================================
 */
    complex         *a, *a_save;
    int_t            *asub, *asub_save;
    int_t            *xa, *xa_save;
    SuperMatrix    A, B, X, L, U;
    SuperMatrix    ASAV, AC;
    int_t            *perm_r; /* row permutation from partial pivoting */
    int_t            *perm_c, *pc_save; /* column permutation */
    complex       zero = {0.0, 0.0};
    float         *R, *C;
    float         *ferr, *berr;
    float         *rwork;
    complex	   *wwork;
    float          diag_pivot_thresh, drop_tol;
    void           *work;
    int_t            info, lwork, nrhs, panel_size, relax;
    int_t            nprocs, m, n, nnz;
    int              n32, info32;
    complex        *xact;
    complex        *rhsb, *solx, *bsav;
    int_t            ldb, ldx;
    float          rpg, rcond;
    int            i, j, k1;
    float         rowcnd, colcnd, amax;
    int_t            maxsuper, rowblk, colblk;
    int_t            prefact, dofact, equil, iequed, norefact;
    int            nt, nrun, nfail, nerrs, imat, fimat, nimat;
    int_t            nfact, ifact, nrefact, irefact, nusepr, iusepr, itran;
    int_t            kl, ku, lda;
    int            kl32, ku32, mode, lda32;
    int_t            zerot, izero, ioff;
    float         anorm, cndnum;
    complex         *Afull;
    float         result[NTESTS];
    Gstat_t        Gstat;
    superlumt_options_t superlumt_options;
    superlu_memusage_t superlu_memusage;
    char matrix_type[8];
    char path[4], sym[1], dist[1];
    fact_t fact;
    trans_t trans;
    equed_t equed;
    yes_no_t refact, usepr;
    void parse_command_line();

    /* Fixed set of parameters */
    int      iseed[]  = {1994, 1995, 1996, 1997};
    equed_t  equeds[] = {NOEQUIL, ROW, COL, BOTH};
    fact_t   facts[]  = {FACTORED, DOFACT, EQUILIBRATE};
    yes_no_t refacts[]= {YES, NO};
    trans_t  transs[] = {NOTRANS, TRANS};
    yes_no_t useprs[] = {YES, NO};

    /* Some function prototypes */
    extern int_t pcgst01(int_t, int_t, SuperMatrix *, SuperMatrix *,
                         SuperMatrix *, int_t *, int_t *, float *);
    extern int_t pcgst02(trans_t, int_t, int_t, int_t, SuperMatrix *, complex *,
                         int_t, complex *, int_t, float *resid);
    extern int_t pcgst04(int_t, int_t, complex *, int_t,
                         complex *, int_t, float rcond, float *resid);
    extern int_t pcgst07(trans_t, int_t, int_t, SuperMatrix *, complex *, int_t,
                         complex *, int_t, complex *, int_t,
                         float *, float *, float *);
    extern int clatb4_(char *, int *, int *, int *, char *, int *, int *,
                       float *, int *, float *, char *);
    extern int clatms_(int *, int *, char *, int *, char *, float *d,
                       int *, float *, float *, int *, int *,
                       char *, complex *, int *, complex *, int *);
    extern int_t sp_cconvert(int_t, int_t, complex *, int_t, int_t, int_t,
                           complex *a, int_t *, int_t *, int_t *);


    /* Executable statements */

    strcpy(path, "CGE");
    nrun  = 0;
    nfail = 0;
    nerrs = 0;

    /* Defaults */
    nprocs     = 1;
    n          = 1;
    nrhs       = 1;
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    diag_pivot_thresh = 1.0;
    usepr = NO;
    drop_tol = 0.0;
    work     = NULL;
    lwork    = 0;
    strcpy(matrix_type, "LA");
    parse_command_line(argc, argv, matrix_type, &nprocs, &n,
		       &panel_size, &relax, &nrhs, &maxsuper,
		       &rowblk, &colblk, &lwork);
    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    fprintf(stderr, "expert: cannot allocate " IFMT " bytes\n", lwork);
	    exit (-1);
	}
	bzero(work, lwork);
    }

#if ( DEBUGlevel>=1 )
    printf("n = " IFMT ", nprocs = " IFMT ", w = " IFMT ", relax = " IFMT "\n",
	   n, nprocs, panel_size, relax);
#endif

    if ( strcmp(matrix_type, "LA") == 0 ) {
	/* Test LAPACK matrix suite. */
	m = n;
	lda = SUPERLU_MAX(n, 1);
	nnz = n * n;        /* upper bound */
	fimat = 1;
	nimat = NTYPES;
	Afull = complexCalloc(lda * n);
	callocateA(n, nnz, &a, &asub, &xa);
    } else {
	/* Read a sparse matrix */
	fimat = nimat = 0;
	creadhb(&m, &n, &nnz, &a, &asub, &xa);
    }

    ldb  = m;
    ldx  = n;

    callocateA(n, nnz, &a_save, &asub_save, &xa_save);
    rhsb = complexMalloc(m * nrhs);
    bsav = complexMalloc(m * nrhs);
    solx = complexMalloc(n * nrhs);
    xact = complexMalloc(n * nrhs);
    perm_r  = intMalloc(n);
    perm_c  = intMalloc(n);
    pc_save = intMalloc(n);
    R       = (float *) SUPERLU_MALLOC(m*sizeof(float));
    C       = (float *) SUPERLU_MALLOC(n*sizeof(float));
    ferr    = (float *) SUPERLU_MALLOC(nrhs*sizeof(float));
    berr    = (float *) SUPERLU_MALLOC(nrhs*sizeof(float));
    j = SUPERLU_MAX(m,n) * SUPERLU_MAX(4,nrhs);    
    rwork   = (float *) SUPERLU_MALLOC(j*sizeof(float));
    for (i = 0; i < j; ++i) rwork[i] = 0.;
    if ( !R ) SUPERLU_ABORT("SUPERLU_MALLOC fails for R");
    if ( !C ) SUPERLU_ABORT("SUPERLU_MALLOC fails for C");
    if ( !ferr ) SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr");
    if ( !berr ) SUPERLU_ABORT("SUPERLU_MALLOC fails for berr");
    if ( !rwork ) SUPERLU_ABORT("SUPERLU_MALLOC fails for rwork");
    wwork = complexCalloc( SUPERLU_MAX(m,n) * SUPERLU_MAX(4,nrhs) );

    cCreate_Dense_Matrix(&B, m, nrhs, rhsb, ldb, SLU_DN, SLU_C, SLU_GE);
    cCreate_Dense_Matrix(&X, n, nrhs, solx, ldx, SLU_DN, SLU_C, SLU_GE);

    /* Fill in options used by the factorization routine. */
    superlumt_options.panel_size = panel_size;
    superlumt_options.relax = relax;
    superlumt_options.diag_pivot_thresh = diag_pivot_thresh;
    superlumt_options.drop_tol = drop_tol;
    superlumt_options.perm_c = perm_c;
    superlumt_options.perm_r = perm_r;
    superlumt_options.work = work;
    superlumt_options.lwork = lwork;

    for (i = 0; i < n; ++i) perm_c[i] = pc_save[i] = i;

    /* All matrix types */
    for (imat = fimat; imat <= nimat; ++imat) {
	
	if ( imat ) {

	    /* Skip types 5, 6, or 7 if the matrix size is too small. */
	    zerot = (imat >= 5 && imat <= 7);
	    if ( zerot && n < imat-4 )
		continue;
	    
	    /* Set up parameters with CLATB4 and generate a test matrix
	       with CLATMS.  */
            n32 = n;  kl32 = kl; ku32 = ku;
	    clatb4_(path, &imat, &n32, &n32, sym, &kl32, &ku32, &anorm, &mode,
		    &cndnum, dist);
	    kl = kl32; ku = ku32;

            lda32 = lda;
	    clatms_(&n32, &n32, dist, iseed, sym, &rwork[0], &mode, &cndnum,
		    &anorm, &kl32, &ku32, "No packing", Afull, &lda32,
		    &wwork[0], &info32);
            info = info32;
            
	    if ( info ) {
		printf(FMT3, "CLATMS", info, izero, n, nrhs, imat, nfail);
		continue;
	    }

	    /* For types 5-7, zero one or more columns of the matrix
	       to test that INFO is returned correctly.   */
	    if ( zerot ) {
		if ( imat == 5 ) izero = 1;
		else if ( imat == 6 ) izero = n;
		else izero = n / 2 + 1;
		ioff = (izero - 1) * lda;
		if ( imat < 7 ) {
		    for (i = 0; i < n; ++i) Afull[ioff + i] = zero;
		} else {
		    for (j = 0; j < n - izero + 1; ++j)
			for (i = 0; i < n; ++i)
			    Afull[ioff + i + j*lda] = zero;
		}
	    } else {
		izero = 0;
	    }

	    /* Convert to sparse representation. */
	    sp_cconvert(n, n, Afull, lda, kl, ku, a, asub, xa, &nnz);

	} else {
	    izero = 0;
	    zerot = 0;
	}
	
	cCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_C, SLU_GE);
	/*Print_CompCol_NC(&A);*/

	/* Save a copy of matrix A in ASAV */
	cCreate_CompCol_Matrix(&ASAV, m, n, nnz, a_save, asub_save, xa_save,
			       SLU_NC, SLU_C, SLU_GE);
	cCopy_CompCol_Matrix(&A, &ASAV);
	
	/* Form exact solution. */
	cGenXtrue(n, nrhs, xact, ldx);
	
	for (iequed = 0; iequed < 4; ++iequed) {
	    equed = equeds[iequed];
	    if (iequed == 0) nfact = 3;
	    else nfact = 1;

	    for (ifact = 0; ifact < nfact; ++ifact) {
		fact = facts[ifact];
		if (ifact == 0) nrefact = nusepr = 1;
		else nrefact = nusepr = 2;

		for (irefact = 0; irefact < nrefact; ++irefact) {
		    refact = refacts[irefact];
		    for (iusepr = 0; iusepr < nusepr; ++iusepr) {
		        usepr = useprs[iusepr];
			norefact  = ( refact == NO );
			prefact   = ( fact == FACTORED || refact == YES 
				      || usepr == YES );
			dofact    = ( fact == DOFACT );
			equil     = ( fact == EQUILIBRATE );

			/* Restore the matrix A. */
			cCopy_CompCol_Matrix(&ASAV, &A);
			
#if ( DEBUGlevel>=1 )
    printf("imat=%2d, equed=%2d, fact=%2d, refact=%2d, usepr=%2d\n",
	   imat, equed, fact, refact, usepr);
    fflush(stdout);
#endif
			if ( zerot ) {
			    if ( prefact ) continue;
			} else if ( ! dofact ) {
			    if ( equil || iequed ) {
			      /* Compute row and column scale factors to
				 equilibrate matrix A.    */
			        cgsequ(&A, R, C, &rowcnd, &colcnd,&amax,&info);

				/* Force equilibration. */
				if ( !info && n > 0 ) {
				    if ( equed == ROW ) {
					rowcnd = 0.;
					colcnd = 1.;
				    } else if ( equed == COL ) {
					rowcnd = 1.;
					colcnd = 0.;
				    } else if ( equed == BOTH ) {
					rowcnd = 0.;
					colcnd = 0.;
				    }
				}
			
				/* Equilibrate the matrix. */
				claqgs(&A, R, C, rowcnd, colcnd, amax, &equed);
			    }
			}
		    
			if ( prefact ) {	/* First time factor */
			
			    StatAlloc(n, nprocs, panel_size, relax, &Gstat);
			    StatInit(n, nprocs, &Gstat);
			
			    /* Initialize options, reorder the matrix, 
			       obtain the column etree. 
                               First time always use partial pivoting. */
		            pcgstrf_init(nprocs, fact, NOTRANS, NO, 
                                         panel_size, relax, 1.0, NO, drop_tol,
					 perm_c, perm_r, work, lwork,
					 &A, &AC, &superlumt_options, &Gstat);

#if ( DEBUGlevel>=1 )
			    printf("Test PCGSTRF\n");
#endif
			    /* Factor the matrix AC. */
			    pcgstrf(&superlumt_options, &AC, perm_r, &L, &U, 
				    &Gstat, &info);

			    /* Restore parameters for subsequent factors. */
			    superlumt_options.refact = refact;
			    superlumt_options.usepr = usepr;
			    superlumt_options.diag_pivot_thresh=diag_pivot_thresh;

			    if ( info ) { 
			        printf("** First factor: info " IFMT ", equed %d\n",
				       info, equed);
				if ( lwork == -1 ) {
				  printf("** Estimated memory: " IFMT " bytes\n",
					 info - n);
				  exit(0);
				}
			    }
	
			    /* Deallocate storage. */
			    Destroy_CompCol_Permuted(&AC);
			    StatFree(&Gstat);
			
			} /* if .. first time factor */
		    
			for (itran = 0; itran < NTRAN; ++itran) {
			    trans = transs[itran];

			    /* Restore the matrix A. */
			    cCopy_CompCol_Matrix(&ASAV, &A);
			
			    /* Set the right hand side. */
			    cFillRHS(trans, nrhs, xact, ldx, &A, &B);
			    cCopy_Dense_Matrix(m, nrhs, rhsb, ldb, bsav, ldb);

			    /*----------------
			     * Test pcgssv
			     *----------------*/
			    if ( dofact && norefact && itran == 0) {
			        /* Not yet factored, and untransposed */
			        cCopy_Dense_Matrix(m, nrhs, rhsb, ldb,
						   solx, ldx);

				pcgssv(nprocs, &A, perm_c, perm_r, &L, &U, 
				       &X, &info);

#if ( DEBUGlevel>=1 )
			    printf("Test PCGSSV: info = " IFMT "\n", info);
#endif
				if ( info && info != izero ) {
                                    printf(FMT3, "pcgssv",
					   info, izero, n, nrhs, imat, nfail);
				} else {
                                    /* Reconstruct matrix from factors and
				       compute residual. */
                                    pcgst01(m, n, &A, &L, &U, perm_c, perm_r,
					    &result[0]);

				    nt = 1;
				    if ( izero == 0 ) {
				        /* Compute residual of the computed
					   solution. */
				        cCopy_Dense_Matrix(m, nrhs, rhsb, ldb,
							   wwork, ldb);
					pcgst02(trans, m, n, nrhs, &A, solx,
						ldx, wwork,ldb, &result[1]);
					nt = 2;
				    }
				
				    /* Print information about the tests that
				       did not pass the threshold.      */
				    for (i = 0; i < nt; ++i) {
				        if ( result[i] >= THRESH ) {
					    printf(FMT1, "pcgssv", n, i,
						   result[i]);
					    ++nfail;
					}
				    }
				    nrun += nt;
				} /* else .. info == 0 */

				/* Restore perm_c. */
				for (i = 0; i < n; ++i) perm_c[i] = pc_save[i];

				if ( lwork == 0 ) {
				  Destroy_SuperNode_SCP(&L);
				  Destroy_CompCol_NCP(&U);
				}
				
			    } /* if .. end of testing pcgssv */

			    /*----------------
			     * Test pcgssvx
			     *----------------*/
			    
			    /* Equilibrate the matrix if fact = FACTORED and
			       equed = ROW, COL or BOTH.   */
			    if ( iequed > 0 && n > 0 ) {
			        claqgs(&A, R, C, rowcnd, colcnd, amax, &equed);
			    }
			
			    /* Initialize more options that are changed from
			       iteration to iteration. */
			    superlumt_options.nprocs = nprocs;
			    superlumt_options.fact = fact;
			    superlumt_options.trans = trans;
			    superlumt_options.refact = refact;
			    superlumt_options.usepr = usepr;
			    if ( !prefact ) {
			       if ( !(superlumt_options.etree = intMalloc(n)) )
			       	  SUPERLU_ABORT("Malloc fails for etree[].");
			       if ( !(superlumt_options.colcnt_h = intMalloc(n)) )
				  SUPERLU_ABORT("Malloc fails for colcnt_h[].");
			       if ( !(superlumt_options.part_super_h = intMalloc(n)) )
				  SUPERLU_ABORT("Malloc fails for colcnt_h[].");
			    }

			    /* Solve the system and compute the condition 
			       number and error bounds using pcgssvx.   */
			    pcgssvx(nprocs, &superlumt_options, &A, perm_c, 
				    perm_r, &equed, R, C, &L, &U, &B, &X,
				    &rpg, &rcond, ferr, berr, 
				    &superlu_memusage, &info);

#if ( DEBUGlevel>=1 )
			    printf("Test PCGSSVX: info " IFMT "\n", info);
#endif
			    if ( info && info != izero ) {
			        printf(FMT3, "pcgssvx",
				       info, izero, n, nrhs, imat, nfail);
				if ( lwork == -1 ) {
				    printf("** Estimated memory: %.0f bytes\n",
					   superlu_memusage.total_needed);
				    exit(0);
				}
			    } else {
			        if ( !prefact ) {
				    /* Reconstruct matrix from factors and
				       compute residual. */
				    pcgst01(m, n, &A, &L, &U, perm_c, perm_r, 
					    &result[0]);
				    k1 = 0;
				} else {
				  k1 = 1;
				}

				if ( !info ) {
				    /* Compute residual of the computed 
				       solution.*/
				    cCopy_Dense_Matrix(m, nrhs, bsav, ldb,
						       wwork, ldb);
				    pcgst02(trans, m, n, nrhs, &ASAV, solx,
					    ldx, wwork, ldb, &result[1]);

				    /* Check solution from generated exact
				       solution. */
				    pcgst04(n, nrhs, solx, ldx, xact, ldx,
					    rcond, &result[2]);

				    /* Check the error bounds from iterative
				       refinement. */
				    pcgst07(trans, n, nrhs, &ASAV, bsav, ldb,
					    solx, ldx, xact, ldx, ferr, berr,
					    &result[3]);

				    /* Print information about the tests that
				       did not pass the threshold.    */
				    for (i = k1; i < NTESTS; ++i) {
				        if ( result[i] >= THRESH ) {
					  printf(FMT2, "pcgssvx",
						 fact, trans, refact, equed,
						 n, imat, i, result[i]);
					  ++nfail;
					}
				    }
				    nrun += NTESTS;
				} /* if .. info == 0 */
			    } /* else .. end of testing pcgssvx */

			    if ( !prefact ) {
				if ( lwork == 0 ) {
				    Destroy_SuperNode_SCP(&L);
				    Destroy_CompCol_NCP(&U);
				} else if ( lwork > 0 ) {
				    bzero(work, lwork);
				    /*for (i = 0; i < lwork; ++i) 
				        ((char*)work)[i] = 0;*/
	                            Destroy_SuperMatrix_Store(&L);
	                            Destroy_SuperMatrix_Store(&U);
				}
			    }
			} /* for itran ... */

			if ( prefact ) {
			    if ( lwork == 0 ) {
			        Destroy_SuperNode_SCP(&L);
				Destroy_CompCol_NCP(&U);
			    } else if ( lwork > 0 ) {
			        bzero(work,lwork);
			        /*for (i = 0; i < lwork; ++i)
				   ((char*)work)[i] = 0;*/
	                        Destroy_SuperMatrix_Store(&L);
	                        Destroy_SuperMatrix_Store(&U);
			    }
			}
			SUPERLU_FREE(superlumt_options.etree);
			SUPERLU_FREE(superlumt_options.colcnt_h);
			SUPERLU_FREE(superlumt_options.part_super_h);
		    } /* for iusepr ... */
		} /* for irefact ... */
	    } /* for ifact ... */
	} /* for iequed ... */
#if 0    
    if ( !info ) {
	PrintPerf(&L, &U, &superlu_memusage, rpg, rcond, ferr, berr, equed);
    }
#endif    
	/* Deallocate some storage to prepare the test of next matrix. */
	Destroy_SuperMatrix_Store(&A);
	Destroy_SuperMatrix_Store(&ASAV);

    } /* for imat ... */

    /* Print a summary of the results. */
    PrintSumm("CGE", nfail, nrun, nerrs);

    SUPERLU_FREE (a);
    SUPERLU_FREE (asub);
    SUPERLU_FREE (xa);
    SUPERLU_FREE (a_save);
    SUPERLU_FREE (asub_save);
    SUPERLU_FREE (xa_save);
    if ( strcmp(matrix_type, "LA") == 0 ) SUPERLU_FREE (Afull);
    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (bsav);
    SUPERLU_FREE (solx);    
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (pc_save);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    SUPERLU_FREE (rwork);
    SUPERLU_FREE (wwork);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if ( lwork > 0 ) {
	SUPERLU_FREE (work);
    }

    return 0;
}

/*  
 * Parse command line options to get relaxed snode size, panel size, etc.
 */
void
parse_command_line(int argc, char *argv[], char *matrix_type,
		   int *nprocs, int *n, int *w, int *relax, int *nrhs,
		   int *maxsuper, int *rowblk, int *colblk, int *lwork)
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "ht:n:p:w:r:s:m:b:c:l:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-t <char*> - matrix type (\"LA\" or \"SP\")\n");
	    printf("\t-p <int> - number of processes\n");
	    printf("\t-w <int> - panel size\n");
	    printf("\t-r <int> - granularity of relaxed supernodes\n");
	    printf("\t-s <int> - number of right-hand sides\n");
	    printf("\t-m <int> - maximum size of a supernode\n");
	    printf("\t-b <int> - row block size in 2D partition\n");
	    printf("\t-c <int> - column block size in 2D partition\n");
	    printf("\t-l <int> - length of work[] array\n");
	    exit(1);
	    break;
	  case 't': strcpy(matrix_type, optarg);
	            break;
	  case 'n': *n = atoi(optarg);
	            break;
	  case 'p': *nprocs = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'r': *relax = atoi(optarg); 
	            break;
	  case 's': *nrhs = atoi(optarg); 
	            break;
	  case 'm': *maxsuper = atoi(optarg); 
	            break;
	  case 'b': *rowblk = atoi(optarg); 
	            break;
	  case 'c': *colblk = atoi(optarg); 
	            break;
	  case 'l': *lwork = atoi(optarg); 
	            break;
  	}
    }
}
