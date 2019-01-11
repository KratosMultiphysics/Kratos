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
 * Purpose:
 * ========
 *
 * This example illustrates how to use PCGSSVX to solve systems repeatedly
 * with the same sparsity pattern of matrix A.
 * In this case, the column permutation vector perm_c is computed once.
 * The following data structures will be reused in the subsequent call to
 * PCGSSVX: perm_c, etree, colcnt_h, part_super_h.
 *
 */
#include "slu_mt_cdefs.h"


main(int argc, char *argv[])
{
    SuperMatrix A, A1, L, U;
    SuperMatrix B, B1, X;
    NCformat    *Astore;
    SCPformat   *Lstore;
    NCPformat   *Ustore;
    int_t         nprocs;
    fact_t      fact;
    trans_t     trans;
    yes_no_t    refact, usepr;
    equed_t     equed;
    complex      *a, *a1;
    int_t         *asub, *xa, *asub1, *xa1;
    int_t         *perm_c; /* column permutation vector */
    int_t         *perm_r; /* row permutations from partial pivoting */
    void        *work;
    superlumt_options_t superlumt_options;
    int_t         info, lwork, nrhs, ldx, panel_size, relax;
    int_t         m, n, nnz, permc_spec, i;
    complex      *rhsb, *rhsx, *xact;
    float      *R, *C;
    float      *ferr, *berr;
    float      u, drop_tol, rpg, rcond;
    superlu_memusage_t superlu_memusage;
    void parse_command_line();

    /* Default parameters to control factorization. */
    nprocs = 1;
    fact  = EQUILIBRATE;
    trans = NOTRANS;
    equed = NOEQUIL;
    refact= NO;
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    u     = 1.0;
    usepr = NO;
    drop_tol = 0.0;
    lwork = 0;
    nrhs  = 1;

    /* Command line options to modify default behavior. */
    parse_command_line(argc, argv, &nprocs, &lwork, &panel_size, &relax, 
		       &u, &fact, &trans, &refact, &equed);

    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
	if ( !work ) {
	    SUPERLU_ABORT("CLINSOLX: cannot allocate work[]");
	}
    }

#if ( PRNTlevel==1 )
    cpp_defs();
    printf("int_t %d bytes\n", sizeof(int_t));
#endif

#define HB
#if defined( DEN )
    m = n;
    nnz = n * n;
    cband(n, n, nnz, &a, &asub, &xa);
#elif defined( BAND )
    m = n;
    nnz = (2*b+1) * n;
    cband(n, b, nnz, &a, &asub, &xa);
#elif defined( BD )
    nb = 5;
    bs = 200;
    m = n = bs * nb;
    nnz = bs * bs * nb;
    cblockdiag(nb, bs, nnz, &a, &asub, &xa);
#elif defined( HB )
    creadhb(&m, &n, &nnz, &a, &asub, &xa);
#else    
    creadtriple(&m, &n, &nnz, &a, &asub, &xa);
#endif

    if ( !(a1 = complexMalloc(nnz)) ) SUPERLU_ABORT("Malloc fails for a1[].");
    if ( !(asub1 = intMalloc(nnz)) ) SUPERLU_ABORT("Malloc fails for asub1[].");
    if ( !(xa1 = intMalloc(n+1)) ) SUPERLU_ABORT("Malloc fails for xa1[].");
    for (i = 0; i < nnz; ++i) {
        a1[i] = a[i];
	asub1[i] = asub[i];
    }
    for (i = 0; i < n+1; ++i) xa1[i] = xa[i];

    cCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_C, SLU_GE);
    Astore = A.Store;
    printf("Dimension " IFMT "x" IFMT "; # nonzeros " IFMT "\n", A.nrow, A.ncol, Astore->nnz);
    
    if (!(rhsb = complexMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsb[].");
    if (!(rhsx = complexMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsx[].");
    cCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_C, SLU_GE);
    cCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_C, SLU_GE);
    xact = complexMalloc(n * nrhs);
    ldx = n;
    cGenXtrue(n, nrhs, xact, ldx);
    cFillRHS(trans, nrhs, xact, ldx, &A, &B);
    
    if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
    if (!(R = (float *) SUPERLU_MALLOC(A.nrow * sizeof(float)))) 
        SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (float *) SUPERLU_MALLOC(A.ncol * sizeof(float))) )
        SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (float *) SUPERLU_MALLOC(nrhs * sizeof(float))) )
        SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (float *) SUPERLU_MALLOC(nrhs * sizeof(float))) ) 
        SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
    permc_spec = 1;
    get_perm_c(permc_spec, &A, perm_c);

    superlumt_options.nprocs = nprocs;
    superlumt_options.fact = fact;
    superlumt_options.trans = trans;
    superlumt_options.refact = refact;
    superlumt_options.panel_size = panel_size;
    superlumt_options.relax = relax;
    superlumt_options.usepr = usepr;
    superlumt_options.drop_tol = drop_tol;
    superlumt_options.diag_pivot_thresh = u;
    superlumt_options.SymmetricMode = NO;
    superlumt_options.PrintStat = NO;
    superlumt_options.perm_c = perm_c;
    superlumt_options.perm_r = perm_r;
    superlumt_options.work = work;
    superlumt_options.lwork = lwork;
    if ( !(superlumt_options.etree = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for etree[].");
    if ( !(superlumt_options.colcnt_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    if ( !(superlumt_options.part_super_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    
    /* ------------------------------------------------------------
       WE SOLVE THE LINEAR SYSTEM FOR THE FIRST TIME: AX = B
       ------------------------------------------------------------*/
    pcgssvx(nprocs, &superlumt_options, &A, perm_c, perm_r,
	    &equed, R, C, &L, &U, &B, &X, &rpg, &rcond,
	    ferr, berr, &superlu_memusage, &info);

    if ( info == 0 || info == n+1 ) {

	printf("Recip. pivot growth = %e\n", rpg);
	printf("Recip. condition number = %e\n", rcond);
	printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
	for (i = 0; i < nrhs; ++i) {
	    printf(IFMT "%16e%16e\n", i+1, ferr[i], berr[i]);
	}
	       
        Lstore = (SCPformat *) L.Store;
        Ustore = (NCPformat *) U.Store;
	printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
	       superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
	       superlu_memusage.expansions);
	     
	fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: " IFMT " bytes\n", info - n);
    }

    printf("First system: pcgssvx(): info " IFMT "\n----\n", info);

    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);

    /* ------------------------------------------------------------
       NOW WE SOLVE ANOTHER LINEAR SYSTEM: A1*X = B1
       ONLY THE SPARSITY PATTERN OF A1 IS THE SAME AS THAT OF A.
       ------------------------------------------------------------*/
    superlumt_options.refact = YES;
    cCreate_CompCol_Matrix(&A1, m, n, nnz, a1, asub1, xa1, SLU_NC, SLU_C, SLU_GE);
    cCreate_Dense_Matrix(&B1, m, nrhs, rhsb, m, SLU_DN, SLU_C, SLU_GE);

    pcgssvx(nprocs, &superlumt_options, &A1, perm_c, perm_r,
	    &equed, R, C, &L, &U, &B1, &X, &rpg, &rcond,
	    ferr, berr, &superlu_memusage, &info);

    if ( info == 0 || info == n+1 ) {

	printf("Recip. pivot growth = %e\n", rpg);
	printf("Recip. condition number = %e\n", rcond);
	printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
	for (i = 0; i < nrhs; ++i) {
	    printf(IFMT "%16e%16e\n", i+1, ferr[i], berr[i]);
	}
	       
        Lstore = (SCPformat *) L.Store;
        Ustore = (NCPformat *) U.Store;
	printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
	       superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
	       superlu_memusage.expansions);
	     
	fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: " IFMT " bytes\n", info - n);
    }

    printf("Second system: pcgssvx(): info " IFMT "\n", info);

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (rhsx);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A1);
    Destroy_SuperMatrix_Store(&B1);
    Destroy_SuperMatrix_Store(&X);
    SUPERLU_FREE (superlumt_options.etree);
    SUPERLU_FREE (superlumt_options.colcnt_h);
    SUPERLU_FREE (superlumt_options.part_super_h);
    if ( lwork == 0 ) {
        Destroy_SuperNode_SCP(&L);
        Destroy_CompCol_NCP(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }
}

/*  
 * Parse command line options.
 */
void
parse_command_line(int argc, char *argv[], int_t *nprocs, int_t *lwork, 
		   int_t *w, int_t *relax, double *u, fact_t *fact, 
		   trans_t *trans, yes_no_t *refact, equed_t *equed)
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "hp:l:w:x:u:f:t:r:e:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-p <int> - number of processes\n");
	    printf("\t-l <int> - length of work[*] array\n");
	    printf("\t-w <int> - panel size\n");
	    printf("\t-x <int> - maximum size of relaxed supernodes\n");
	    printf("\t-u <int> - pivoting threshold\n");
	    printf("\t-f <FACTORED/DOFACT/EQUILIBRATE> - factor control\n");
	    printf("\t-t <NOTRANS/TRANS/CONJ> - transpose or not\n");
	    printf("\t-r <NO/YES> - refactor or not\n");
	    printf("\t-e <NOEQUIL/ROW/COL/BOTH> - equilibrate or not\n");
	    exit(1);
	    break;
	  case 'p': *nprocs = atoi(optarg);
	            break;
	  case 'l': *lwork = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'x': *relax = atoi(optarg); 
	            break;
	  case 'u': *u = atof(optarg); 
	            break;
	  case 'f': *fact = (fact_t) atoi(optarg);
	            break;
	  case 't': *trans = (trans_t) atoi(optarg);
	            break;
	  case 'r': *refact = (yes_no_t) atoi(optarg);
	            break;
	  case 'e': *equed = (equed_t) atoi(optarg);
	            break;
	  default: fprintf(stderr, "Invalid command line option.\n");
		   break;
  	}
    }
}
