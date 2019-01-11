/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include "slu_mt_cdefs.h"

void
cgstrs(trans_t trans, SuperMatrix *L, SuperMatrix *U, 
       int_t *perm_r, int_t *perm_c, SuperMatrix *B, Gstat_t *Gstat, int_t *info)
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
 * cgstrs() solves a system of linear equations A*X=B or A'*X=B
 * with A sparse and B dense, using the LU factorization computed by
 * pcgstrf().
 *
 * Arguments
 * =========
 *
 * trans   (input) Specifies the form of the system of equations:
 *          = NOTRANS: A * X = B  (No transpose)
 *          = TRANS:   A'* X = B  (Transpose)
 *
 * L       (input) SuperMatrix*
 *         The factor L from the factorization Pr*A*Pc=L*U as computed by
 *         pcgstrf(). Use compressed row subscripts storage for supernodes,
 *         i.e., L has types: Stype = SCP, Dtype = _D, Mtype = TRLU.
 *
 * U       (input) SuperMatrix*
 *         The factor U from the factorization Pr*A*Pc=L*U as computed by
 *         pcgstrf(). Use column-wise storage scheme, i.e., U has types:
 *         Stype = NCP, Dtype = _D, Mtype = TRU.
 *
 * perm_r  (input) int_t*
 *         Row permutation vector of size L->nrow, which defines the
 *         permutation matrix Pr; perm_r[i] = j means row i of A is in
 *         position j in Pr*A.
 *
 * perm_c  (int_t*) dimension A->ncol
 *	   Column permutation vector, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = DN, Dtype = _D, Mtype = GE.
 *         On entry, the right hand side matrix.
 *         On exit, the solution matrix if info = 0;
 *
 * Gstat   (output) Gstat_t*
 *          Record all the statistics about the triangular solves; 
 *          See Gstat_t structure defined in slu_mt_util.h.
 *
 * info    (output) Diagnostics
 * 	   = 0: successful exit
 *	   < 0: if info = -i, the i-th argument had an illegal value
 *
 */
#if ( MACH==CRAY_PVP )
    _fcd ftcs1, ftcs2, ftcs3, ftcs4;
#endif

#ifdef USE_VENDOR_BLAS
    int      incx = 1, incy = 1;
    complex   alpha = {1.0, 0.0}, beta = {1.0, 0.0};
#endif

    complex   temp_comp;
    register int_t j, k, jcol, iptr, luptr, ksupno, istart, irow, bptr;
    register int_t fsupc, nsuper;
    int        i, n, nsupc, nsupr, nrow, nrhs, ldb;
    int_t      *supno;
    DNformat *Bstore;
    SCPformat *Lstore;
    NCPformat *Ustore;
    complex   *Lval, *Uval, *Bmat;
    complex   *work, *work_col, *rhs_work, *soln;
    flops_t  solve_ops;
    void cprint_soln();

    /* Test input parameters ... */
    *info = 0;
    Bstore = B->Store;
    ldb = Bstore->lda;
    nrhs = B->ncol;
    if ( trans != NOTRANS && trans != TRANS && trans != CONJ ) *info = -1;
    else if ( L->nrow != L->ncol || L->nrow < 0 ) *info = -3;
    else if ( U->nrow != U->ncol || U->nrow < 0 ) *info = -4;
    else if ( ldb < SUPERLU_MAX(0, L->nrow) ) *info = -6;
    if ( *info ) {
        i = -(*info);
	xerbla_("cgstrs", &i);
	return;
    }

    n = L->nrow;
    work = complexCalloc(n * nrhs);
    if ( !work ) SUPERLU_ABORT("Malloc fails for local work[].");
    soln = complexMalloc(n);
    if ( !soln ) SUPERLU_ABORT("Malloc fails for local soln[].");

    Bmat = Bstore->nzval;
    Lstore = L->Store;
    Lval = Lstore->nzval;
    Ustore = U->Store;
    Uval = Ustore->nzval;
    supno = Lstore->col_to_sup;
    nsuper = Lstore->nsuper;
    solve_ops = 0;
    
    if ( trans == NOTRANS ) {
	/* Permute right hand sides to form Pr*B */
	for (i = 0, bptr = 0; i < nrhs; i++, bptr += ldb) {
	    rhs_work = &Bmat[bptr];
	    for (k = 0; k < n; k++) soln[perm_r[k]] = rhs_work[k];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}
	
	/* Forward solve PLy=Pb. */
/*>>	for (k = 0; k < n; k += nsupc) {
	    ksupno = supno[k];
*/
	for (ksupno = 0; ksupno <= nsuper; ++ksupno) {
	    fsupc = L_FST_SUPC(ksupno);
	    istart = L_SUB_START(fsupc);
	    nsupr = L_SUB_END(fsupc) - istart;
	    nsupc = L_LAST_SUPC(ksupno) - fsupc;
	    nrow = nsupr - nsupc;

	    solve_ops += nsupc * (nsupc - 1) * nrhs;
	    solve_ops += 2 * nrow * nsupc * nrhs;
	    
	    if ( nsupc == 1 ) {
		for (j = 0, bptr = 0; j < nrhs; j++, bptr += ldb) {
		    rhs_work = &Bmat[bptr];
	    	    luptr = L_NZ_START(fsupc);
		    for (iptr=istart+1; iptr < L_SUB_END(fsupc); iptr++){
			irow = L_SUB(iptr);
			++luptr;
                        cc_mult(&temp_comp, &rhs_work[fsupc], &Lval[luptr]);
                        c_sub(&rhs_work[irow], &rhs_work[irow], &temp_comp);
		    }
		}
	    } else {
	    	luptr = L_NZ_START(fsupc);
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
		ftcs1 = _cptofcd("L", strlen("L"));
		ftcs2 = _cptofcd("N", strlen("N"));
		ftcs3 = _cptofcd("U", strlen("U"));
 		CTRSM(ftcs1, ftcs1, ftcs2, ftcs3, &nsupc, &nrhs, &alpha,
		      &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
		
		CGEMM(ftcs2, ftcs2,  &nrow, &nrhs, &nsupc, &alpha, 
		      &Lval[luptr+nsupc], &nsupr, &Bmat[fsupc], &ldb, 
		      &beta, &work[0], &n );
#else
 		ctrsm_("L", "L", "N", "U", &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
		
		cgemm_( "N", "N", &nrow, &nrhs, &nsupc, &alpha, 
			&Lval[luptr+nsupc], &nsupr, &Bmat[fsupc], &ldb, 
			&beta, &work[0], &n );
#endif
		for (j = 0, bptr = 0; j < nrhs; j++, bptr += ldb) {
		    rhs_work = &Bmat[bptr];
		    work_col = &work[j*n];
		    iptr = istart + nsupc;
		    for (i = 0; i < nrow; i++) {
			irow = L_SUB(iptr);
                        c_sub(&rhs_work[irow], &rhs_work[irow], &work_col[i]);
                        work_col[i].r = 0.0;
                        work_col[i].i = 0.0;
			iptr++;
		    }
		}
#else		
		for (j = 0, bptr = 0; j < nrhs; j++, bptr += ldb) {
		    rhs_work = &Bmat[bptr];
		    clsolve ((int_t)nsupr, (int_t)nsupc, &Lval[luptr], &rhs_work[fsupc]);
		    cmatvec ((int_t)nsupr, (int_t)nrow, (int_t)nsupc, &Lval[luptr+nsupc],
			     &rhs_work[fsupc], &work[0] );

		    iptr = istart + nsupc;
		    for (i = 0; i < nrow; i++) {
			irow = L_SUB(iptr);
                        c_sub(&rhs_work[irow], &rhs_work[irow], &work[i]);
                        work[i].r = 0.0;
                        work[i].i = 0.0;
			iptr++;
		    }
		}
#endif		    
	    } /* if-else: nsupc == 1 ... */
	} /* for L-solve */

#if ( DEBUGlevel>=2 )
  	printf("After L-solve: y=\n");
	cprint_soln(n, nrhs, Bmat);
#endif

	/*
	 * Back solve Ux=y.
	 */
/*>>	for (k = n-1; k >= 0; k -= nsupc) {
	    ksupno = supno[k];
*/
	for (ksupno = nsuper; ksupno >= 0; --ksupno) {
	    fsupc = L_FST_SUPC(ksupno);
	    istart = L_SUB_START(fsupc);
	    nsupr = L_SUB_END(fsupc) - istart;
	    nsupc = L_LAST_SUPC(ksupno) - fsupc;
	    luptr = L_NZ_START(fsupc);

	    solve_ops += nsupc * (nsupc + 1) * nrhs;

	    /* dense triangular matrix */
	    if ( nsupc == 1 ) {
		rhs_work = &Bmat[0];
		for (j = 0; j < nrhs; j++) {
                    c_div(&rhs_work[fsupc], &rhs_work[fsupc], &Lval[luptr]);
		    rhs_work += ldb;
		}
	    } else {
#ifdef USE_VENDOR_BLAS
#if ( MACH==CRAY_PVP )
		ftcs1 = _cptofcd("L", strlen("L"));
		ftcs2 = _cptofcd("U", strlen("U"));
		ftcs3 = _cptofcd("N", strlen("N"));
		CTRSM(ftcs1, ftcs2, ftcs3, ftcs3, &nsupc, &nrhs, &alpha,
		      &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#else
		ctrsm_("L", "U", "N", "N", &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#endif
#else		
		for (j = 0, bptr = fsupc; j < nrhs; j++, bptr += ldb) {
		    cusolve (nsupr, nsupc, &Lval[luptr], &Bmat[bptr]);
		}
#endif		
	    }

	    /* matrix-vector update */
	    for (j = 0, bptr = 0; j < nrhs; ++j, bptr += ldb) {
		rhs_work = &Bmat[bptr];
		for (jcol = fsupc; jcol < fsupc + nsupc; jcol++) {
                    solve_ops += 8*(U_NZ_END(jcol) - U_NZ_START(jcol));
		    for (i = U_NZ_START(jcol); i < U_NZ_END(jcol); i++ ){
			irow = U_SUB(i);
                        cc_mult(&temp_comp, &rhs_work[jcol], &Uval[i]);
                        c_sub(&rhs_work[irow], &rhs_work[irow], &temp_comp);
		    }
		}
	    }
	    
	} /* for U-solve */

#if ( DEBUGlevel>=2 )
  	printf("After U-solve: x=\n");
	cprint_soln(n, nrhs, Bmat);
#endif

	/* Compute the final solution X <= Pc*X. */
	for (i = 0, bptr = 0; i < nrhs; i++, bptr += ldb) {
	    rhs_work = &Bmat[bptr];
	    for (k = 0; k < n; k++) soln[k] = rhs_work[perm_c[k]];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}
	
    } else { /* Solve A'*X=B */
	/* Permute right hand sides to form Pc'*B. */
	for (i = 0, bptr = 0; i < nrhs; i++, bptr += ldb) {
	    rhs_work = &Bmat[bptr];
	    for (k = 0; k < n; k++) soln[perm_c[k]] = rhs_work[k];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}
	
        if (trans == TRANS) {
            for (k = 0; k < nrhs; ++k) {
                /* Multiply by inv(U'). */
                sp_ctrsv("U", "T", "N", L, U, &Bmat[k*ldb], info);

                /* Multiply by inv(L'). */
                sp_ctrsv("L", "T", "U", L, U, &Bmat[k*ldb], info);
            }
         } else { /* trans == CONJ */
            for (k = 0; k < nrhs; ++k) {
                /* Multiply by conj(inv(U')). */
                sp_ctrsv("U", "C", "N", L, U, &Bmat[k*ldb], info);

                /* Multiply by conj(inv(L')). */
                sp_ctrsv("L", "C", "U", L, U, &Bmat[k*ldb], info);
            }
         }
	/* Compute the final solution X <= Pr'*X (=inv(Pr)*X) */
	for (i = 0, bptr = 0; i < nrhs; i++, bptr += ldb) {
	    rhs_work = &Bmat[bptr];
	    for (k = 0; k < n; k++) soln[k] = rhs_work[perm_r[k]];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}

    } /* if-else trans */

    Gstat->ops[TRISOLVE] = solve_ops;
    SUPERLU_FREE(work);
    SUPERLU_FREE(soln);
}

/*
 * Diagnostic print of the solution vector
 */
void
cprint_soln(int_t n, int_t nrhs, complex *soln)
{
    int_t i;

    for (i = 0; i < n; i++)
	printf("\t" IFMT ": %.10f\t%.10f\n", i, soln[i].r, soln[i].i);
}
