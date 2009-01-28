/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * These macros define which machine will be used.
 */

#ifndef __SUPERLU_MACHINES /* allow multiple inclusions */
#define __SUPERLU_MACHINES

#define SGI	        0
#define ORIGIN	        1
#define DEC	        2
#define CRAY_PVP	3
#define SUN             4
#define PTHREAD         5
#define OPENMP		6

#ifdef _SGI
#define MACH SGI 
#endif

#ifdef _ORIGIN
#define MACH ORIGIN 
#endif

#ifdef _DEC
#define MACH DEC 
#endif

#ifdef _CRAY
#define MACH CRAY_PVP 
#endif

#ifdef _SOLARIS
#define MACH SUN 
#endif

#ifdef _PTHREAD
#define MACH PTHREAD
#endif

#ifdef __OPENMP
#define MACH OPENMP
#endif

#endif /* __SUPERLU_MACHINES */
