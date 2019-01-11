/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/* 
 * -- SuperLU MT routine (version 2.2) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Last modified: August 18, 2014
 *
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */
#include "slu_mt_machines.h"

#if ( MACH==SUN )
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#else

#if ( MACH==OPENMP )

#include <omp.h>
double SuperLU_timer_()
{
    return omp_get_wtime();
}

#else

double SuperLU_timer_()
{
    double dclock();
    return (dclock());
}

#endif

#endif
