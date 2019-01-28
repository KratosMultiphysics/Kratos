/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

#include <math.h>
#include "slu_mt_ddefs.h"


double
dPivotGrowth(int_t ncols, SuperMatrix *A, int_t *perm_c, 
             SuperMatrix *L, SuperMatrix *U)
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
 * Compute the reciprocal pivot growth factor of the leading ncols columns
 * of the matrix, using the formula:
 *     min_j ( max_i(abs(A_ij)) / max_i(abs(U_ij)) )
 *
 * Arguments
 * =========
 *
 * ncols    (input) int_t
 *          The number of columns of matrices A, L and U.
 *
 * A        (input) SuperMatrix*
 *          Original matrix A, permuted by columns, of dimension
 *          (A->nrow, A->ncol). The type of A can be:
 *          Stype = NC; Dtype = _D; Mtype = GE.
 *
 * L        (output) SuperMatrix*
 *          The factor L from the factorization Pr*A=L*U; use compressed row
 *          subscripts storage for supernodes, i.e., L has type:
 *          Stype = SC; Dtype = _D; Mtype = TRLU.
 *
 * U        (output) SuperMatrix*
 *          The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
 *          storage scheme, i.e., U has types: Stype = NC;
 *          Dtype = _D; Mtype = TRU.
 *
 */
    NCformat *Astore;
    SCPformat *Lstore;
    NCPformat *Ustore;
    double  *Aval, *Lval, *Uval;
    int_t      fsupc, nsupr, luptr, nz_in_U;
    int_t      i, j, k, oldcol;
    int_t      *inv_perm_c;
    double   rpg, maxaj, maxuj;
    extern   double dlamch_(char *);
    double   smlnum;
    double   *luval;
   
    /* Get machine constants. */
    smlnum = dlamch_("S");
    rpg = 1. / smlnum;

    Astore = A->Store;
    Lstore = L->Store;
    Ustore = U->Store;
    Aval = Astore->nzval;
    Lval = Lstore->nzval;
    Uval = Ustore->nzval;
    
    inv_perm_c = intMalloc(A->ncol);
    for (j = 0; j < A->ncol; ++j) inv_perm_c[perm_c[j]] = j;

    for (k = 0; k <= Lstore->nsuper; ++k) {
	fsupc = L_FST_SUPC(k);
	nsupr = L_SUB_END(fsupc) - L_SUB_START(fsupc);
	luptr = L_NZ_START(fsupc);
	luval = &Lval[luptr];
	nz_in_U = 1;
	
	for (j = fsupc; j < L_LAST_SUPC(k) && j < ncols; ++j) {
	    maxaj = 0.;
            oldcol = inv_perm_c[j];
	    for (i = Astore->colptr[oldcol]; i < Astore->colptr[oldcol+1]; i++)
		maxaj = SUPERLU_MAX( maxaj, fabs(Aval[i]) );
	
	    maxuj = 0.;
	    for (i = Ustore->colbeg[j]; i < Ustore->colend[j]; i++)
		maxuj = SUPERLU_MAX( maxuj, fabs(Uval[i]) );
	    
	    /* Supernode */
	    for (i = 0; i < nz_in_U; ++i)
		maxuj = SUPERLU_MAX( maxuj, fabs(luval[i]) );

	    ++nz_in_U;
	    luval += nsupr;

	    if ( maxuj == 0. )
		rpg = SUPERLU_MIN( rpg, 1.);
	    else
		rpg = SUPERLU_MIN( rpg, maxaj / maxuj );
	}
	
	if ( j >= ncols ) break;
    }

    SUPERLU_FREE(inv_perm_c);
    return (rpg);
}
