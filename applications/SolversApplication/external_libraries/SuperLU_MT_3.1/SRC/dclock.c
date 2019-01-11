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
 */

#ifdef WIN32

#include <windows.h>
double dclock()
{
    LARGE_INTEGER time, frequency;
    QueryPerformanceCounter(&time);
    QueryPerformanceFrequency(&frequency) ;
    return ((double)time.QuadPart /(double)frequency.QuadPart);
}

/* double dclock()
{
    return (double)clock() / (double)CLOCKS_PER_SEC;
    } */

double usertimer_()
{
    return dclock();
}


#else  /* default to unix */

#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double usertimer_()
{
    struct tms use;
    double tmp;
    int clocks_per_sec = sysconf(_SC_CLK_TCK);

    times ( &use );
    tmp = use.tms_utime;
    tmp += use.tms_stime;
    return (double)(tmp) / clocks_per_sec;
}

double extract(tv)
struct timeval *tv;
{
  double tmp;

  tmp = tv->tv_sec;
  tmp += tv->tv_usec/1.0e6;

  return(tmp);
}

double dclock()
{
    struct timeval tp;
    struct timezone tzp;

    /* wall-clock time */
    gettimeofday(&tp,&tzp);

    return(extract(&tp));
}

#endif
