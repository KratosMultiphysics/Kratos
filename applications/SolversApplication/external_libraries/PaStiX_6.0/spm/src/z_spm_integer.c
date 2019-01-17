/**
 *
 * @file z_spm_integer.c
 *
 * SParse Matrix package integer sorting routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @fn      void z_spmIntFltSortAsc(void ** const pbase, const spm_int_t n)
 * @ingroup spm_dev_integer
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * spm_int_t and used as key for sorting.  The second array is an array of
 * spm_complex64_t.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Couple of pointers to an array of integers and to an array of
 *          spm_complex64_t to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
static size_t intsortsize_if[2] = { sizeof(spm_int_t), sizeof(spm_complex64_t) };
#define INTSORTNAME            z_spmIntFltSortAsc
#define INTSORTSIZE(x)         (intsortsize_if[x])
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {                                     \
        spm_int_t     t;                                                \
        long    disp_p   = (((spm_int_t*)p)-((spm_int_t*)base_ptr));    \
        long    disp_q   = (((spm_int_t*)q)-((spm_int_t*)base_ptr));    \
        spm_complex64_t * floatptr = *(pbase+1);                        \
        spm_complex64_t   f;                                            \
        /* swap integers */                                             \
        t = *((spm_int_t *) (p));                                       \
        *((spm_int_t *) (p)) = *((spm_int_t *) (q));                    \
        *((spm_int_t *) (q)) = t;                                       \
        /* swap corresponding values */                                 \
        f = floatptr[disp_p];                                           \
        floatptr[disp_p] = floatptr[disp_q];                            \
        floatptr[disp_q] = f;                                           \
    } while (0)
#define INTSORTCMP(p,q)             (*((spm_int_t *) (p)) < *((spm_int_t *) (q)))
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @fn      void z_spmIntIntFltSortAsc(void ** const pbase, const spm_int_t n)
 * @ingroup spm_dev_integer
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * spm_int_t and used as key for sorting.  The second array is an array of
 * spm_complex64_t.
 *
 *******************************************************************************
 *
 * @param[inout] pbase
 *          Couple of pointers to an array of integers and to an array of
 *          spm_complex64_t to sort.
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 *******************************************************************************
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
static size_t intsortsize_iif[3] = { sizeof(spm_int_t), sizeof(spm_int_t), sizeof(spm_complex64_t) };
#define INTSORTNAME            z_spmIntIntFltSortAsc
#define INTSORTSIZE(x)         (intsortsize_iif[x])
#define INTSORTNTAB            3
#define INTSORTSWAP(p,q)       do {                                     \
        spm_int_t     t;                                                \
        long    disp_p   = (((spm_int_t*)p)-((spm_int_t*)base_ptr));    \
        long    disp_q   = (((spm_int_t*)q)-((spm_int_t*)base_ptr));    \
        spm_int_t       * intptr = *(pbase+1);                          \
        spm_complex64_t * fltptr = *(pbase+2);                          \
        spm_complex64_t   f;                                            \
        /* swap integers */                                             \
        t = *((spm_int_t *) (p));                                       \
        *((spm_int_t *) (p)) = *((spm_int_t *) (q));                    \
        *((spm_int_t *) (q)) = t;                                       \
        /* swap on second integer array */                              \
        t = intptr[disp_p];                                             \
        intptr[disp_p] = intptr[disp_q];                                \
        intptr[disp_q] = t;                                             \
        /* swap corresponding values */                                 \
        f = fltptr[disp_p];                                             \
        fltptr[disp_p] = fltptr[disp_q];                                \
        fltptr[disp_q] = f;                                             \
    } while (0)

static inline int
intsortcmp_iif( void ** const pbase, spm_int_t *p, spm_int_t *q ) {
    spm_int_t *int1ptr = pbase[0];
    spm_int_t *int2ptr = pbase[1];
    return ( *p < *q ) || (( *p == *q ) && ( int2ptr[ p - int1ptr ] < int2ptr[ q - int1ptr ] ));
}
#define INTSORTCMP(p,q) intsortcmp_iif( pbase, (spm_int_t*)p, (spm_int_t*)q )
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
