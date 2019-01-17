/**
 *
 * @file spm_datatypes.h
 *
 * @copyright 2013-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Definitions of the datatypes used in SPM
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @date 2017-01-17
 *
 */
#ifndef _spm_datatypes_h_
#define _spm_datatypes_h_

#include <inttypes.h>

/** ****************************************************************************
 * Integers
 */
#if defined(SPM_INT64)

typedef int64_t  spm_int_t;
typedef uint64_t spm_uint_t;
#define SPM_MPI_INT MPI_INTEGER8
#define SPM_INT_MAX INT64_MAX

#elif defined(SPM_INT32)

typedef int32_t  spm_int_t;
typedef uint32_t spm_uint_t;
#define SPM_MPI_INT MPI_INTEGER4
#define SPM_INT_MAX INT32_MAX

#elif defined(SPM_LONG)

typedef long          spm_int_t;
typedef unsigned long spm_uint_t;
#define SPM_MPI_INT MPI_LONG
#define SPM_INT_MAX LONG_MAX

#else

typedef int          spm_int_t;
typedef unsigned int spm_uint_t;
#define SPM_MPI_INT MPI_INT
#define SPM_INT_MAX INT_MAX

#endif

static inline spm_int_t spm_imin( spm_int_t a, spm_int_t b ) {
    return ( a < b ) ? a : b;
}

static inline spm_int_t spm_imax( spm_int_t a, spm_int_t b ) {
    return ( a > b ) ? a : b;
}

static inline spm_int_t spm_iceil( spm_int_t a, spm_int_t b ) {
    return ( a + b - 1 ) / b;
}

/** ****************************************************************************
 * Double that are not converted through precision generator functions
 **/
typedef double spm_fixdbl_t;

/** ****************************************************************************
 * Complex numbers (Extracted from PaRSEC project)
 **/
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
/* Windows and non-Intel compiler */
#include <complex>
typedef std::complex<float>  spm_complex32_t;
typedef std::complex<double> spm_complex64_t;
#else
typedef float  _Complex      spm_complex32_t;
typedef double _Complex      spm_complex64_t;
#endif

#if !defined(__cplusplus) && defined(HAVE_COMPLEX_H)
#include <complex.h>
#else

#ifdef __cplusplus
extern "C" {
#endif

/* These declarations will not clash with what C++ provides because
 * the names in C++ are name-mangled. */

extern double cabs     (spm_complex64_t z);
extern double creal    (spm_complex64_t z);
extern double cimag    (spm_complex64_t z);

extern float  cabsf    (spm_complex32_t z);
extern float  crealf   (spm_complex32_t z);
extern float  cimagf   (spm_complex32_t z);

extern spm_complex64_t conj  (spm_complex64_t z);
extern spm_complex64_t csqrt (spm_complex64_t z);

extern spm_complex32_t conjf (spm_complex32_t z);
extern spm_complex32_t csqrtf(spm_complex32_t z);

#ifdef __cplusplus
}
#endif

#endif /* HAVE_COMPLEX_H */


static inline size_t
spm_size_of(spm_coeftype_t type)
{
    switch(type) {
    case SpmFloat:     return   sizeof(float);
    case SpmDouble:    return   sizeof(double);
    case SpmComplex32: return 2*sizeof(float);
    case SpmComplex64: return 2*sizeof(double);
    default:
        fprintf(stderr, "spm_size_of: invalid type parameter\n");
        assert(0);
        return sizeof(double);
    }
}

#endif /* _spm_datatypes_h_ */
