#include <stdio.h> 
/* 
 * Purpose
 * ======= 
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */


#ifdef _GIVE_UP_THIS_ONE_
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double sys_timer()
{
    return ( (double)gethrtime() / 1e9 );
}

#else

#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>

#ifndef CLK_TCK
#define CLK_TCK 100
#endif

#ifndef CLOCKS_PER_SEC 
#define CLOCKS_PER_SEC 1000000
#endif 

double sys_timer_CLOCK() {
  clock_t tmp;
  tmp = clock();
  return (double) tmp/(CLOCKS_PER_SEC);
}

double sys_timer() {
    struct tms use;
    clock_t tmp;
    times(&use);
    tmp = use.tms_utime + use.tms_stime;
    return (double)(tmp) / CLK_TCK;
}

#endif
