/**
 *
 * @file z_integer.c
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2018-07-16
 *
 * @precisions normal z -> c d s
 *
 **/
/* Copyright 2004,2007-2009 ENSEIRB, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : integer.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the generic integer **/
/**                type.                                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 sep 1998     **/
/**                                 to     22 sep 1998     **/
/**                # Version 0.1  : from : 07 jan 2002     **/
/**                                 to     17 jan 2003     **/
/**                # Version 1.0  : from : 23 aug 2005     **/
/**                                 to   : 19 dec 2006     **/
/**                # Version 2.0  : from : 26 feb 2008     **/
/**                                 to   : 26 feb 2008     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to   : 21 jan 2009     **/
/**                                                        **/
/************************************************************/
#include <ctype.h>
#include <limits.h>
#include <time.h>
#include "common.h"


/*
   Function: qsortIntFloatAsc

   Sort 2 arrays simultaneously, the first array is an
   array of pastix_int_t and used as key for sorting.
   The second array is an array of pastix_complex64_t.

   Parameters:
     pbase       - Array of pointers to the first element of each array to sort.
     total_elems - Number of element in each array.

   Returns:
     Nothing

*/

static size_t intsortsize[2] = { sizeof(pastix_int_t), sizeof(pastix_complex64_t) };
#define INTSORTNAME            z_qsortIntFloatAsc
#define INTSORTSIZE(x)         (intsortsize[x])
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {					\
    pastix_int_t     t;								\
    long    disp_p   = (((pastix_int_t*)p)-((pastix_int_t*)base_ptr));			\
    long    disp_q   = (((pastix_int_t*)q)-((pastix_int_t*)base_ptr));			\
    pastix_complex64_t * floatptr = *(pbase+1);					\
    pastix_complex64_t   f;								\
    /* swap integers */							\
    t = *((pastix_int_t *) (p));							\
    *((pastix_int_t *) (p)) = *((pastix_int_t *) (q));					\
    *((pastix_int_t *) (q)) = t;							\
    /* swap corresponding values */					\
    f = floatptr[disp_p];						\
    floatptr[disp_p] = floatptr[disp_q];				\
    floatptr[disp_q] = f;						\
  } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "common_sort2.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB

