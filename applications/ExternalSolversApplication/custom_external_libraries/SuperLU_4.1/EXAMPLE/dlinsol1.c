
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include "slu_ddefs.h"

main(int argc, char *argv[])
{
    SuperMatrix A;
    NCformat *Astore;
    double   *a;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    double   *xact, *rhs;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
     */
    set_default_options(&options);

    /* Now we modify the default options to use the symmetric mode. */
    options.SymmetricMode = YES;
    options.ColPerm = MMD_AT_PLUS_A;
    options.DiagPivotThresh = 0.001;

#if 1
    /* Read matrix A from a file in Harwell-Boeing format.*/
    if (argc < 2)
    {
	printf("Usage:\n%s [OPTION] < [INPUT] > [OUTPUT]\nOPTION:\n"
		"-h -hb:\n\t[INPUT] is a Harwell-Boeing format matrix.\n"
		"-r -rb:\n\t[INPUT] is a Rutherford-Boeing format matrix.\n"
		"-t -triplet:\n\t[INPUT] is a triplet format matrix.\n",
		argv[0]);
	return 0;
    }
    else
    {
	switch (argv[1][1])
	{
	    case 'H':
	    case 'h':
		printf("Input a Harwell-Boeing format matrix:\n");
		dreadhb(&m, &n, &nnz, &a, &asub, &xa);
		break;
	    case 'R':
	    case 'r':
		printf("Input a Rutherford-Boeing format matrix:\n");
		dreadrb(&m, &n, &nnz, &a, &asub, &xa);
		break;
	    case 'T':
	    case 't':
		printf("Input a triplet format matrix:\n");
		dreadtriple(&m, &n, &nnz, &a, &asub, &xa);
		break;
	    default:
		printf("Unrecognized format.\n");
		return 0;
	}
    }
#else
    /* Read the matrix in Harwell-Boeing format. */
    dreadhb(&m, &n, &nnz, &a, &asub, &xa);
#endif

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);
    ldx = n;
    dGenXtrue(n, nrhs, xact, ldx);
    dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    if ( info == 0 ) {

	/* This is how you could access the solution matrix. */
        double *sol = (double*) ((DNformat*) B.Store)->nzval; 

	 /* Compute the infinity norm of the error. */
	dinf_norm_error(nrhs, &B, xact);

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    	printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);
	
	dQuerySpace(&L, &U, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	
    } else {
	printf("dgssv() error returns INFO= %d\n", info);
	if ( info <= n ) { /* factorization completes */
	    dQuerySpace(&L, &U, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	}
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}

