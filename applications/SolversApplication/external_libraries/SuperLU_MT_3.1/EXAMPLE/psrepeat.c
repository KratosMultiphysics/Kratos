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
 * Last modified: 02/14/2013
 *                04/16/2015
 *
 * Purpose 
 * =======
 *
 * This program illustrates how to perform the multiple 
 * factorizations for the matrix with same sparsity pattern yet
 * different numerical values. Here, computing a fill-reducing
 * ordering and performing column permutation are done only once.
 * In addition, the memory for the factors L and U is allocated
 * once, and re-used in the subsequent factorizations.
 * 
 */
#include "slu_mt_sdefs.h"


main(int argc, char *argv[])
{
    SuperMatrix A, AC, L, U, B;
    NCformat    *Astore;
    SCPformat   *Lstore;
    NCPformat   *Ustore;
    superlumt_options_t superlumt_options;
    pxgstrf_shared_t pxgstrf_shared;
    psgstrf_threadarg_t *psgstrf_threadarg;
    int_t         nprocs;
    fact_t      fact;
    trans_t     trans;
    yes_no_t    refact, usepr;
    float      u, drop_tol;
    float      *a;
    int_t         *asub, *xa;
    int_t         *perm_c; /* column permutation vector */
    int_t         *perm_r; /* row permutations from partial pivoting */
    void        *work;
    int_t         info, lwork, nrhs, ldx; 
    int_t         m, n, nnz, permc_spec, panel_size, relax;
    int_t         i, firstfact;
    float      *rhsb, *xact;
    Gstat_t Gstat;
    flops_t     flopcnt;
    void parse_command_line();

    /* Default parameters to control factorization. */
    nprocs = 1;
    fact  = EQUILIBRATE;
    trans = NOTRANS;
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
    sreadhb(&m, &n, &nnz, &a, &asub, &xa);

    /* Set up the sparse matrix data structure for A. */
    sCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_S, SLU_GE);

    if (!(rhsb = floatMalloc(m * nrhs))) SUPERLU_ABORT("Malloc fails for rhsb[].");
    sCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_S, SLU_GE);
    xact = floatMalloc(n * nrhs);
    ldx = n;
    sGenXtrue(n, nrhs, xact, ldx);
    sFillRHS(trans, nrhs, xact, ldx, &A, &B);
    
    if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");
    if ( !(superlumt_options.etree = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for etree[].");
    if ( !(superlumt_options.colcnt_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");
    if ( !(superlumt_options.part_super_h = intMalloc(n)) )
	SUPERLU_ABORT("Malloc fails for colcnt_h[].");

    /********************************
     * THE FIRST TIME FACTORIZATION *
     ********************************/

    /* ------------------------------------------------------------
       Allocate storage and initialize statistics variables. 
       ------------------------------------------------------------*/
    StatAlloc(n, nprocs, panel_size, relax, &Gstat);
    StatInit(n, nprocs, &Gstat);

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
       Initialize the option structure superlumt_options using the
       user-input parameters;
       Apply perm_c to the columns of original A to form AC.
       ------------------------------------------------------------*/
    refact= NO;
    psgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
		 u, usepr, drop_tol, perm_c, perm_r,
		 work, lwork, &A, &AC, &superlumt_options, &Gstat);

    /* ------------------------------------------------------------
       Compute the LU factorization of A.
       The following routine will create nprocs threads.
       ------------------------------------------------------------*/
    psgstrf(&superlumt_options, &AC, perm_r, &L, &U, &Gstat, &info);
    
    flopcnt = 0;
    for (i = 0; i < nprocs; ++i) flopcnt += Gstat.procstat[i].fcops;
    Gstat.ops[FACT] = flopcnt;

    /* ------------------------------------------------------------
       Solve the system A*X=B, overwriting B with X.
       ------------------------------------------------------------*/
    sgstrs(trans, &L, &U, perm_r, perm_c, &B, &Gstat, &info);
    
    printf("\n** Result of sparse LU **\n");
    sinf_norm_error(nrhs, &B, xact); /* Check inf. norm of the error */

    Destroy_CompCol_Permuted(&AC); /* Free extra arrays in AC. */


    /*********************************
     * THE SUBSEQUENT FACTORIZATIONS *
     *********************************/

    /* ------------------------------------------------------------
       Re-initialize statistics variables and options used by the
       factorization routine psgstrf().
       ------------------------------------------------------------*/
    StatInit(n, nprocs, &Gstat);
    refact= YES;
    psgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
		 u, usepr, drop_tol, perm_c, perm_r,
		 work, lwork, &A, &AC, &superlumt_options, &Gstat);

    /* ------------------------------------------------------------
       Compute the LU factorization of A.
       The following routine will create nprocs threads.
       ------------------------------------------------------------*/
    psgstrf(&superlumt_options, &AC, perm_r, &L, &U, &Gstat, &info);
    
    flopcnt = 0;
    for (i = 0; i < nprocs; ++i) flopcnt += Gstat.procstat[i].fcops;
    Gstat.ops[FACT] = flopcnt;

    /* ------------------------------------------------------------
       Re-generate right-hand side B, then solve A*X= B.
       ------------------------------------------------------------*/
    sFillRHS(trans, nrhs, xact, ldx, &A, &B);
    sgstrs(trans, &L, &U, perm_r, perm_c, &B, &Gstat, &info);

    
     /* ------------------------------------------------------------
       Deallocate storage after factorization.
       ------------------------------------------------------------*/
    pxgstrf_finalize(&superlumt_options, &AC);

    printf("\n** Result of sparse LU **\n");
    sinf_norm_error(nrhs, &B, xact); /* Check inf. norm of the error */

    Lstore = (SCPformat *) L.Store;
    Ustore = (NCPformat *) U.Store;
    printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
    printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
    printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
    fflush(stdout);

    SUPERLU_FREE (rhsb);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (superlumt_options.etree);
    SUPERLU_FREE (superlumt_options.colcnt_h);
    SUPERLU_FREE (superlumt_options.part_super_h);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    if ( lwork == 0 ) {
        Destroy_SuperNode_SCP(&L);
        Destroy_CompCol_NCP(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }
    StatFree(&Gstat);
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

