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
pxgstrf_finalize(superlumt_options_t *superlumt_options, SuperMatrix *AC)
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
 * pxgstrf_finalize() deallocates storage after factorization pxgstrf().
 *
 * Arguments
 * =========
 *
 * superlumt_options (input) superlumt_options_t*
 *        The structure contains the parameters to facilitate sparse
 *        LU factorization.
 *
 * AC (input) SuperMatrix*
 *        The original matrix with columns permuted.
 */
    SUPERLU_FREE(superlumt_options->etree);
    SUPERLU_FREE(superlumt_options->colcnt_h);
    SUPERLU_FREE(superlumt_options->part_super_h);
    Destroy_CompCol_Permuted(AC);
#if ( DEBUGlevel>=1 )
    printf("** pxgstrf_finalize() called\n");
#endif
}
