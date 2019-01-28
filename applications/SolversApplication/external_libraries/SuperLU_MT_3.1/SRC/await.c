/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * This routine should NOT be optimized.
 */
#include "slu_mt_ddefs.h"

int_t await(volatile int_t *status)
{
    register int_t i, j, k, randnum;

    /* randnum = ( random() & 0xff ); */
    randnum = 0;
    while ( *status ) ;
#if 0
    {
	/* Length better be adaptive to the number of processors */
	k = randnum;
	for (i = 0; i < randnum; ++i) {
	    j += k;
	    k = -k;
	}
    }
#endif
    return 0;
}
