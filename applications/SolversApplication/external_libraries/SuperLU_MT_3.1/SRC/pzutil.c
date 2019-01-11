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
 */
#include <math.h>
#include "slu_mt_zdefs.h"

void
zCreate_CompCol_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz, doublecomplex *nzval,
		      int_t *rowind, int_t *colptr,
		      Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NCformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NCformat) );
    Astore = (NCformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colptr = colptr;
}

void
zCreate_CompRow_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz, doublecomplex *nzval,
		      int_t *colind, int_t *rowptr,
		      Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NRformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NRformat) );
    Astore = (NRformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->colind = colind;
    Astore->rowptr = rowptr;
}

void
zCreate_CompCol_Permuted(SuperMatrix *A, int_t m, int_t n, int_t nnz, doublecomplex *nzval,
			 int_t *rowind, int_t *colbeg, int_t *colend,
			 Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NCPformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NCPformat) );
    Astore = (NCPformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colbeg = colbeg;
    Astore->colend = colend;
}
/*
 * Convert a row compressed storage int_to a column compressed storage.
 */
void
zCompRow_to_CompCol(int_t m, int_t n, int_t nnz, 
		    doublecomplex *a, int_t *colind, int_t *rowptr,
		    doublecomplex **at, int_t **rowind, int_t **colptr)
{
    register int_t i, j, col, relpos;
    int_t *marker;

    /* Allocate storage for another copy of the matrix. */
    *at = (doublecomplex *) doublecomplexMalloc(nnz);
    *rowind = intMalloc(nnz);
    *colptr = intMalloc(n+1);
    marker = intCalloc(n);
    
    /* Get counts of each column of A, and set up column pointers */
    for (i = 0; i < m; ++i)
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
    (*colptr)[0] = 0;
    for (j = 0; j < n; ++j) {
	(*colptr)[j+1] = (*colptr)[j] + marker[j];
	marker[j] = (*colptr)[j];
    }

    /* Transfer the matrix into the compressed column storage. */
    for (i = 0; i < m; ++i) {
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
	    col = colind[j];
	    relpos = marker[col];
	    (*rowind)[relpos] = i;
	    (*at)[relpos] = a[j];
	    ++marker[col];
	}
    }

    SUPERLU_FREE(marker);
}


/* Copy matrix A into matrix B. */
void
zCopy_CompCol_Matrix(SuperMatrix *A, SuperMatrix *B)
{
    NCformat *Astore, *Bstore;
    int_t      ncol, nnz, i;

    B->Stype = A->Stype;
    B->Dtype = A->Dtype;
    B->Mtype = A->Mtype;
    B->nrow  = A->nrow;;
    B->ncol  = ncol = A->ncol;
    Astore   = (NCformat *) A->Store;
    Bstore   = (NCformat *) B->Store;
    Bstore->nnz = nnz = Astore->nnz;
    for (i = 0; i < nnz; ++i)
	((doublecomplex *)Bstore->nzval)[i] = ((doublecomplex *)Astore->nzval)[i];
    for (i = 0; i < nnz; ++i) Bstore->rowind[i] = Astore->rowind[i];
    for (i = 0; i <= ncol; ++i) Bstore->colptr[i] = Astore->colptr[i];
}


int_t zPrint_CompCol_Matrix(SuperMatrix *A)
{
    NCformat     *Astore;
    register int_t i;
    doublecomplex       *dp;
    
    printf("\nCompCol matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (NCformat *) A->Store;
    dp = (doublecomplex *) Astore->nzval;
    printf("nrow " IFMT ", ncol " IFMT ", nnz " IFMT "\n", A->nrow,A->ncol,Astore->nnz);
    printf("\nnzval: ");
    for (i = 0; i < Astore->nnz; ++i) printf("%f %f ",
          dp[i].r, dp[i].i);
    printf("\nrowind: ");
    for (i = 0; i < Astore->nnz; ++i) printf(IFMT, Astore->rowind[i]);
    printf("\ncolptr: ");
    for (i = 0; i <= A->ncol; ++i) printf(IFMT, Astore->colptr[i]);
    printf("\nend CompCol matrix.\n");

    return 0;
}

int_t zPrint_CCS_to_triplets(SuperMatrix *A)
{
    NCformat     *Astore;
    register int_t i, j;
    doublecomplex       *dp;
    
    Astore = (NCformat *) A->Store;
    dp = (doublecomplex *) Astore->nzval;
    printf(IFMT IFMT IFMT, A->nrow,A->ncol,Astore->nnz);
    for (j = 0; j < A->ncol; ++j) {
	for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
	    printf(IFMT IFMT "%20.16e %20.16e\n", Astore->rowind[i], j,
	             dp[i].r, dp[i].i);
	}
    }
    return 0;
}

int_t zPrint_Dense_Matrix(SuperMatrix *A)
{
    DNformat     *Astore;
    register int_t i;
    double       *dp;
    
    printf("\nDense matrix: ");
    printf("Stype %d , Dtype %d , Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (DNformat *) A->Store;
    dp = (double *) Astore->nzval;
    printf("nrow " IFMT ", ncol " IFMT ", lda " IFMT "\n",
           A->nrow,A->ncol,Astore->lda);
    printf("\nnzval: ");
    for (i = 0; i < 2*A->nrow; ++i) printf("%f  ", dp[i]);
    printf("\nend Dense matrix.\n");

    return 0;
}

void
zCreate_Dense_Matrix(SuperMatrix *X, int_t m, int_t n, doublecomplex *x, int_t ldx,
		    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    DNformat    *Xstore;
    
    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    X->Store = (void *) SUPERLU_MALLOC( sizeof(DNformat) );
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (doublecomplex *) x;
}

void
zCopy_Dense_Matrix(int_t M, int_t N, doublecomplex *X, int_t ldx, doublecomplex *Y, int_t ldy)
{
/*
 *
 *  Purpose
 *  =======
 *
 *  Copies a two-dimensional matrix X to another matrix Y.
 */
    int_t    i, j;
    
    for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i)
            Y[i + j*ldy] = X[i + j*ldx];
}

void
zCreate_SuperNode_Matrix(SuperMatrix *L, int_t m, int_t n, int_t nnz, doublecomplex *nzval,
			int_t *nzval_colptr, int_t *rowind, int_t *rowind_colptr,
			int_t *col_to_sup, int_t *sup_to_col,
			Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    SCformat *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    L->Store = (void *) SUPERLU_MALLOC( sizeof(SCformat) );
    Lstore = L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colptr = nzval_colptr;
    Lstore->rowind = rowind;
    Lstore->rowind_colptr = rowind_colptr;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_col = sup_to_col;

}

void
zCreate_SuperNode_Permuted(SuperMatrix *L, int_t m, int_t n, int_t nnz,
			   doublecomplex *nzval, 
			   int_t *nzval_colbeg, int_t *nzval_colend,
			   int_t *rowind, int_t *rowind_colbeg, int_t *rowind_colend,
			   int_t *col_to_sup, 
			   int_t *sup_to_colbeg, int_t *sup_to_colend,
			   Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    SCPformat *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    L->Store = (void *) SUPERLU_MALLOC( sizeof(SCPformat) );
    Lstore = L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colbeg = nzval_colbeg;
    Lstore->nzval_colend = nzval_colend;
    Lstore->rowind = rowind;
    Lstore->rowind_colbeg = rowind_colbeg;
    Lstore->rowind_colend = rowind_colend;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_colbeg = sup_to_colbeg;
    Lstore->sup_to_colend = sup_to_colend;

}


/*
 * Diagnostic print of column "jcol" in the U/L factor.
 */
void
zprint_lu_col(int_t pnum, char *msg, int_t pcol, int_t jcol, int_t w, int_t pivrow,
	      int_t *xprune, GlobalLU_t *Glu)
{
    int_t     i, k, fsupc;
    int_t     *xsup, *supno;
    int_t     *xlsub, *xlsub_end, *lsub;
    doublecomplex  *lusup;
    int_t     *xlusup, *xlusup_end;

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    lusup   = Glu->lusup;
    xlusup  = Glu->xlusup;
    xlusup_end = Glu->xlusup_end;
    
    printf("("IFMT") %s fstcol " IFMT ",col " IFMT ",w " IFMT ": pivrow " IFMT ", supno " IFMT ", xprune " IFMT "\n", 
	   pnum, msg, pcol, jcol, w, pivrow, supno[jcol], xprune[jcol]);
    
    printf("(" IFMT ")\tU-col: xusub " IFMT " - " IFMT "\n",
	   pnum, Glu->xusub[jcol], Glu->xusub_end[jcol]);
    for (i = Glu->xusub[jcol]; i < Glu->xusub_end[jcol]; i++)
	printf("(" IFMT ")\t" IFMT "\t%8e\t%8e\n", pnum, Glu->usub[i],
	 Glu->ucol[i].r, Glu->ucol[i].i);
    fsupc = xsup[supno[jcol]];
    k = xlusup[jcol];
    printf("(" IFMT ")\tL-col in s-node: xlsub " IFMT " - " IFMT ", xlusup " IFMT "-" IFMT "\n",
	   pnum, xlsub[fsupc],xlsub_end[fsupc],xlusup[jcol],xlusup_end[jcol]);
    for (i = xlsub[fsupc]; i < xlsub_end[fsupc]; ++i)
	printf("(" IFMT ")\t" IFMT "\t%.8e\t%.8e\n", pnum, lsub[i], 
	lusup[k++].r, lusup[k++].i);

    fflush(stdout);
}

/*Dan fix above printf 's*/

/*
 * Check whether vec[*] == 0. For the two vectors dense[*] and tempv[*],
 * this invariant should be mantained before and after calling some
 * numeric update routines, such as "panel_bmod" and "column_bmod". 
 */
void
zcheck_zero_vec(int_t pnum, char *msg, int_t n, doublecomplex *vec)
{
    register int_t i, nonzero;

    nonzero = FALSE;
    for (i = 0; i < n; ++i) {
        if ((vec[i].r != 0.0) || (vec[i].i != 0.0))
        {
            printf("(" IFMT ") vec[" IFMT "] = %.10e\t%.10e; should be zero!\n",
                   pnum, i, vec[i].r, vec[i].i);
            nonzero = TRUE;
        }
    }
    if ( nonzero ) {
	printf("(" IFMT ") %s\n", pnum, msg);
	SUPERLU_ABORT("Not a zero vector.");
    } else {
        printf(".. Normal exit zcheck_zero_vec() ..\n");
    }
}


void
zGenXtrue(int_t n, int_t nrhs, doublecomplex *x, int_t ldx)
{
    int_t  i, j;
    for (j = 0; j < nrhs; ++j) {
	for (i = 0; i < n; ++i) {
            x[i + j*ldx].r = 1.0;
            x[i + j*ldx].i = 0.0;
        }
    }
}

/*
 * Let rhs[i] = sum of i-th row of A, so the solution vector is all 1's
 */
void
zFillRHS(trans_t trans, int_t nrhs, doublecomplex *x, int_t ldx, SuperMatrix *A, SuperMatrix *B)
{
    NCformat *Astore;
    doublecomplex   *Aval;
    DNformat *Bstore;
    doublecomplex   *rhs;
    doublecomplex one = {1.0, 0.0};
    doublecomplex zero = {0.0, 0.0};
    int_t      ldc;
    char     trans_c[1];

    Astore = A->Store;
    Aval   = (doublecomplex *) Astore->nzval;
    Bstore = B->Store;
    rhs    = Bstore->nzval;
    ldc    = Bstore->lda;
    
    if ( trans == NOTRANS ) *trans_c = 'N';
    else *trans_c = 'T';
    
    sp_zgemm(trans_c, A->nrow, nrhs, A->ncol, one, A,
	     x, ldx, zero, rhs, ldc);
}

/* 
 * Fills a double precision array with a given value.
 */
void 
zfill(doublecomplex *a, int_t alen, doublecomplex dval)
{
    register int_t i;
    for (i = 0; i < alen; i++) a[i] = dval;
}



/* 
 * Check the inf-norm of the error vector 
 */
void zinf_norm_error(int_t nrhs, SuperMatrix *X, doublecomplex *xtrue)
{
    DNformat *Xstore;
    double err, xnorm;
    doublecomplex *Xmat, *soln_work;
    doublecomplex temp;
    int_t i, j;

    Xstore = X->Store;
    Xmat = Xstore->nzval;

    for (j = 0; j < nrhs; j++) {
      soln_work = &Xmat[j*Xstore->lda];
      err = xnorm = 0.0;
      for (i = 0; i < X->nrow; i++) {
        z_sub(&temp, &soln_work[i], &xtrue[i]);
        err = SUPERLU_MAX(err, z_abs(&temp));
        xnorm = SUPERLU_MAX(xnorm, z_abs(&soln_work[i]));
      }
      err = err / xnorm;
      printf("||X - Xtrue||/||X|| = %e\n", err);
    }
}



/* Print performance of the code. */
void
zPrintPerf(SuperMatrix *L, SuperMatrix *U, superlu_memusage_t *superlu_memusage, 
	double rpg, double rcond, double *ferr,
	double *berr, char *equed, Gstat_t *Gstat)
{
    SCPformat *Lstore;
    NCPformat *Ustore;
    double   *utime;
    flops_t  *ops;
    
    utime = Gstat->utime;
    ops   = Gstat->ops;
    
    if ( utime[FACT] != 0. )
	printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT],
	       ops[FACT]*1e-6/utime[FACT]);
    printf("Identify relaxed snodes	= %8.2f\n", utime[RELAX]);
    if ( utime[SOLVE] != 0. )
	printf("Solve flops = %.0f, Mflops = %8.2f\n", ops[SOLVE],
	       ops[SOLVE]*1e-6/utime[SOLVE]);
    
    Lstore = (SCPformat *) L->Store;
    Ustore = (NCPformat *) U->Store;
    printf("\t#NZ in factor L = " IFMT "\n", Lstore->nnz);
    printf("\t#NZ in factor U = " IFMT "\n", Ustore->nnz);
    printf("\t#NZ in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - L->ncol);
	
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
	   superlu_memusage->for_lu/1e6, superlu_memusage->total_needed/1e6,
	   superlu_memusage->expansions);
	
    printf("\tFactor\tMflops\tSolve\tMflops\tEtree\tEquil\tRcond\tRefine\n");
    printf("PERF:%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
	   utime[FACT], ops[FACT]*1e-6/utime[FACT],
	   utime[SOLVE], ops[SOLVE]*1e-6/utime[SOLVE],
	   utime[ETREE], utime[EQUIL], utime[RCOND], utime[REFINE]);
    
    printf("\tRpg\t\tRcond\t\tFerr\t\tBerr\t\tEquil?\n");
    printf("NUM:\t%e\t%e\t%e\t%e\t%s\n",
	   rpg, rcond, ferr[0], berr[0], equed);
    
#if 0

    printf("\tTRSV (total%%)\tGEMV (total%%)\tfloat_time%%\tmax_n\tmax_m\tmin_n\tmin_m\tavg_n\tavg_m\n");
    printf("BLAS:\t%.0f  %.2f\t%.0f  %.2f\t%.2f\t\t" IFMT "\t" IFMT "\t" IFMT "\t" IFMT "\t%.0f\t%.0f\n",
	   ops[TRSV], ops[TRSV]/ops[FACT], ops[GEMV], ops[GEMV]/ops[FACT],
	   utime[FLOAT]/utime[FACT],
	   max_blas_n, max_gemv_m, min_blas_n, min_gemv_m,
	   (float)sum_blas_n/num_blas, (float)sum_gemv_m/num_blas);
    printf("\tRCOND\tREFINE\tFERR\n");
    printf("SOLVES:\t" IFMT "\t" IFMT "\t" IFMT "\n", no_solves[RCOND],
	   no_solves[REFINE], no_solves[FERR]);
    
    flops_dist_for_matlab();

#endif
    
}


int_t print_doublecomplex_vec(char *what, int_t n, int_t *ind, doublecomplex *vec)
{
    int_t i;
    printf("%s: n " IFMT "\n", what, n);
    for (i = 0; i < n; ++i) printf(IFMT "\t%f%f\n", ind[i], vec[i].r, vec[i].i);
    return 0;
}

