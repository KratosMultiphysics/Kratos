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
sp_colorder(SuperMatrix *A, int_t *perm_c, superlumt_options_t *options,
	    SuperMatrix *AC)
{
/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 *
 * Purpose
 * =======
 *
 * sp_colorder() permutes the columns of the original matrix A into AC. 
 * It performs the following steps:
 *
 *    1. Apply column permutation perm_c[] to A's column pointers to form AC;
 *
 *    2. If options->refact = NO, then
 *       (1) Compute column etree etree[] of AC'AC;
 *       (2) Post order etree[] to get a postordered elimination tree etree[],
 *           and a postorder permutation post[];
 *       (3) Apply post[] permutation to columns of AC;
 *       (4) Overwrite perm_c[] with the product perm_c * post.
 *       (5) Compute the column count (colcnt_h) and the supernode 
 *           partition (part_super_h) for the Householder matrix H.
 *
 * Arguments
 * =========
 *
 * A      (input) SuperMatrix*
 *        Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *        of the linear equations is A->nrow. Currently, the type of A can be:
 *        Stype = NC or NCP; Dtype = _D; Mtype = GE.
 *
 * perm_c (input/output) int_t*
 *	  Column permutation vector of size A->ncol, which defines the 
 *        permutation matrix Pc; perm_c[i] = j means column i of A is 
 *        in position j in A*Pc.
 *
 * options (input/output) superlumt_options_t*
 *        If options->refact = YES, then options is an
 *        input argument. The arrays etree[], colcnt_h[] and part_super_h[]
 *        are available from a previous factor and will be re-used.
 *        If options->refact = NO, then options is an output argument. 
 *
 * AC     (output) SuperMatrix*
 *        The resulting matrix after applied the column permutation
 *        perm_c[] to matrix A. The type of AC can be:
 *        Stype = NCP; Dtype = _D; Mtype = GE.
 *
 */

    NCformat  *Astore;
    NCPformat *ACstore;
    int_t i, n, nnz, nlnz;
    yes_no_t  refact = options->refact;
    int_t *iwork, *post, *iperm;
    int_t *invp;
    int_t *part_super_ata;

    extern void at_plus_a(const int_t, const int_t, int_t *, int_t *,
			  int_t *, int_t **, int_t **, int_t);

    n     = A->ncol;

    /* Apply column permutation perm_c to A's column pointers so to
       obtain NCP format in AC = A*Pc.  */
    AC->Stype       = SLU_NCP;
    AC->Dtype       = A->Dtype;
    AC->Mtype       = A->Mtype;
    AC->nrow        = A->nrow;
    AC->ncol        = A->ncol;
    Astore          = A->Store;
    ACstore = AC->Store = (void *) malloc( sizeof(NCPformat) );
    ACstore->nnz    = Astore->nnz;
    ACstore->nzval  = Astore->nzval;
    ACstore->rowind = Astore->rowind;
    ACstore->colbeg = intMalloc(n);
    ACstore->colend = intMalloc(n);
    nnz             = Astore->nnz;

#ifdef CHK_COLORDER
    print_int_vec("pre_order:", n, perm_c);
    dcheck_perm("Initial perm_c", n, perm_c);
#endif      

    for (i = 0; i < n; i++) {
	ACstore->colbeg[perm_c[i]] = Astore->colptr[i]; 
	ACstore->colend[perm_c[i]] = Astore->colptr[i+1];
    }
	
    if ( refact == NO ) {
	int_t *b_colptr, *b_rowind, bnz, j;
	
	iwork = intMalloc(n+1);
	part_super_ata = intMalloc(n);
    
	if ( options->SymmetricMode ) {
	    /* Compute the etree of C = Pc*(A'+A)*Pc'. */
	    int_t *c_colbeg, *c_colend;

	    /* Form B = A + A'. */
	    at_plus_a(n, Astore->nnz, Astore->colptr, Astore->rowind,
		      &bnz, &b_colptr, &b_rowind, 1);
	    
	    /* Form C = Pc*B*Pc'. */
	    c_colbeg = (int_t*) intMalloc(n);
	    c_colend = (int_t*) intMalloc(n);
	    if (!(c_colbeg) || !(c_colend) )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for c_colbeg/c_colend");
	    for (i = 0; i < n; i++) {
		c_colbeg[perm_c[i]] = b_colptr[i]; 
		c_colend[perm_c[i]] = b_colptr[i+1];
	    }
	    for (j = 0; j < n; ++j) {
		for (i = c_colbeg[j]; i < c_colend[j]; ++i) {
		    b_rowind[i] = perm_c[b_rowind[i]];
		}
		iwork[perm_c[j]] = j; /* inverse perm_c */
	    }

	    /* Compute etree of C. */
	    sp_symetree( c_colbeg, c_colend, b_rowind, n, 
			 options->etree );

	    /* Restore B to be A+A', without column permutation */
	    for (i = 0; i < bnz; ++i)
		b_rowind[i] = iwork[b_rowind[i]];

	    SUPERLU_FREE(c_colbeg);
	    SUPERLU_FREE(c_colend);
	    
	} else {
	    /* Compute the column elimination tree. */
	    sp_coletree( ACstore->colbeg, ACstore->colend, ACstore->rowind,
			 A->nrow, A->ncol, options->etree );
	}

#ifdef CHK_COLORDER	
	print_int_vec("etree:", n, otpions->etree);
#endif	

	/* Post order etree. */
	post = (int_t *) TreePostorder(n, options->etree);
	invp  = intMalloc(n);
	for (i = 0; i < n; ++i) invp[post[i]] = i;

#ifdef CHK_COLORDER
	print_int_vec("post:", n+1, post);
	dcheck_perm("post", n, post);	
#endif	

	/* Renumber etree in postorder. */
	for (i = 0; i < n; ++i) iwork[post[i]] = post[options->etree[i]];
	for (i = 0; i < n; ++i) options->etree[i] = iwork[i];

#ifdef CHK_COLORDER	
	print_int_vec("postorder etree:", n, options->etree);
#endif

	/* Postmultiply A*Pc by post[]. */
	for (i = 0; i < n; ++i) iwork[post[i]] = ACstore->colbeg[i];
	for (i = 0; i < n; ++i) ACstore->colbeg[i] = iwork[i];
	for (i = 0; i < n; ++i) iwork[post[i]] = ACstore->colend[i];
	for (i = 0; i < n; ++i) ACstore->colend[i] = iwork[i];

	for (i = 0; i < n; ++i)
	    iwork[i] = post[perm_c[i]];  /* product of perm_c and post */
	for (i = 0; i < n; ++i) perm_c[i] = iwork[i];
	for (i = 0; i < n; ++i) invp[perm_c[i]] = i; /* inverse of perm_c */

	iperm = post; /* alias to the same address */

#ifdef ZFD_PERM
	/* Permute the rows of AC to have zero-free diagonal. */
	printf("** Permute the rows to have zero-free diagonal....\n");
	for (i = 0; i < n; ++i)
	    iwork[i] = ACstore->colend[i] - ACstore->colbeg[i];
	zfdperm(n, nnz, ACstore->rowind, ACstore->colbeg, iwork, iperm);
#else
	for (i = 0; i < n; ++i) iperm[i] = i;
#endif	

	/* NOTE: iperm is returned as column permutation so that
	 * the diagonal is nonzero. Since a symmetric permutation
	 * preserves the diagonal, we can do the following:
	 *     P'(AP')P = P'A
	 * That is, we apply the inverse of iperm to rows of A
	 * to get zero-free diagonal. But since iperm is defined
	 * in MC21A inversely as our definition of permutation,
	 * so it is indeed an inverse for our purpose. We can
	 * apply it directly.
	 */

	if ( options->SymmetricMode ) {
	    /* Determine column count in the Cholesky factor of B = A+A' */
#if 0
	    cholnzcnt(n, Astore->colptr, Astore->rowind,
		      invp, perm_c, options->etree, 
		      options->colcnt_h, &nlnz, options->part_super_h);
#else
	    cholnzcnt(n, b_colptr, b_rowind, invp, perm_c,
		      options->etree, options->colcnt_h, &nlnz,
		      options->part_super_h);
#endif
#if ( PRNTlevel>=1 ) 
	    printf(".. bnz %d\n", bnz);
#endif

	    SUPERLU_FREE(b_colptr);
	    if ( bnz ) SUPERLU_FREE(b_rowind);

	} else {
	    /* Determine the row and column counts in the QR factor. */
	    qrnzcnt(n, nnz, Astore->colptr, Astore->rowind, iperm,
		    invp, perm_c, options->etree, options->colcnt_h, &nlnz,
		    part_super_ata, options->part_super_h);
	}

#if ( PRNTlevel>=2 )
	dCheckZeroDiagonal(n, ACstore->rowind, ACstore->colbeg,
			   ACstore->colend, perm_c);
	print_int_vec("colcnt", n, options->colcnt_h);
	dPrintSuperPart("Hpart", n, options->part_super_h);
	print_int_vec("iperm", n, iperm);
#endif	
	
#ifdef CHK_COLORDER
	print_int_vec("Pc*post:", n, perm_c);
	dcheck_perm("final perm_c", n, perm_c);	
#endif

	SUPERLU_FREE (post);
	SUPERLU_FREE (invp);
	SUPERLU_FREE (iwork);
	SUPERLU_FREE (part_super_ata);

    } /* if refact == NO */

}

int_t
dCheckZeroDiagonal(int_t n, int_t *rowind, int_t *colbeg,
		  int_t *colend, int_t *perm)
{
    register int_t i, j, nzd, nd = 0;

    for (j = 0; j < n; ++j) {
	nzd = 0;
	for (i = colbeg[j]; i < colend[j]; ++i) {
	    if ( perm[rowind[i]] == j ) {
		nzd = 1;
		++nd;
		break;
	    }
	}
	if ( nzd == 0 ) printf("Zero diagonal at column " IFMT "\n", j);
    }

    printf(".. dCheckZeroDiagonal() -- # diagonals " IFMT "\n", nd);

    return 0;
}

int_t
dPrintSuperPart(char *pname, int_t n, int_t *part_super)
{
    register int_t i;
    FILE *fopen(), *fp;
    char fname[20];
    strcpy(fname, pname);
    strcat(fname, ".dat");
    fp = fopen(fname, "w");
    for (i = 0; i < n; ++i)
	if ( part_super[i] ) fprintf(fp, IFMT, i);
    fprintf(fp, IFMT, n);
    fclose(fp);
    return 0;
}

int_t dcheck_perm(char *what, int_t n, int_t *perm)
{
    register int_t i;
    int_t          *marker;
    marker = (int_t *) intCalloc(n);

    for (i = 0; i < n; ++i) {
	if ( marker[perm[i]] == 1 || perm[i] >= n ) {
	    printf("%s: Not a valid PERM[" IFMT "] = " IFMT "\n", what, i, perm[i]);
	    SUPERLU_ABORT("Invalid perm.");
	} else {
	    marker[perm[i]] = 1;
	}
    }

    return 0;
}
