/* 
 * -- SuperLU MT routine (alpha version) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */


#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double SuperLU_timer_() {
    return ( (double)gethrtime() / 1e9 );
}

#else

double SuperLU_timer_()
{
    double dclock();
    return (dclock());
}

#endif


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
  tmp += tv->tv_usec/1000000.0;
 
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

