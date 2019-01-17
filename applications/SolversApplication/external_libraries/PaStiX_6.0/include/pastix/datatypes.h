/**
 *
 * @file datatypes.h
 *
 * @copyright 2013-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * Definitions of the datatypes used in PaStiX
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 */
#ifndef _pastix_datatypes_h_
#define _pastix_datatypes_h_

#include <inttypes.h>

/** ****************************************************************************
 * Integers
 **/
#if defined(PASTIX_INT64)

typedef int64_t  pastix_int_t;
typedef uint64_t pastix_uint_t;
#define PASTIX_MPI_INT MPI_INTEGER8
#define PASTIX_INT_MAX INT64_MAX

#elif defined(PASTIX_INT32)

typedef int32_t  pastix_int_t;
typedef uint32_t pastix_uint_t;
#define PASTIX_MPI_INT MPI_INTEGER4
#define PASTIX_INT_MAX INT32_MAX

#elif defined(PASTIX_LONG)

typedef long          pastix_int_t;
typedef unsigned long pastix_uint_t;
#define PASTIX_MPI_INT MPI_LONG
#define PASTIX_INT_MAX LONG_MAX

#else

typedef int          pastix_int_t;
typedef unsigned int pastix_uint_t;
#define PASTIX_MPI_INT MPI_INT
#define PASTIX_INT_MAX INT_MAX

#endif


/** ****************************************************************************
 * Double that are not converted through precision generator functions
 **/
typedef double pastix_fixdbl_t;

/** ****************************************************************************
 * Complex numbers (Extracted from PaRSEC project)
 **/
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
/* Windows and non-Intel compiler */
#include <complex>
typedef std::complex<float>  pastix_complex32_t;
typedef std::complex<double> pastix_complex64_t;
#else
typedef float  _Complex      pastix_complex32_t;
typedef double _Complex      pastix_complex64_t;
#endif

#if !defined(__cplusplus) && defined(HAVE_COMPLEX_H)
#include <complex.h>
#else

#ifdef __cplusplus
extern "C" {
#endif

/* These declarations will not clash with what C++ provides because
 * the names in C++ are name-mangled. */

extern double cabs     (pastix_complex64_t z);
extern double creal    (pastix_complex64_t z);
extern double cimag    (pastix_complex64_t z);

extern float  cabsf    (pastix_complex32_t z);
extern float  crealf   (pastix_complex32_t z);
extern float  cimagf   (pastix_complex32_t z);

extern pastix_complex64_t conj  (pastix_complex64_t z);
extern pastix_complex64_t csqrt (pastix_complex64_t z);

extern pastix_complex32_t conjf (pastix_complex32_t z);
extern pastix_complex32_t csqrtf(pastix_complex32_t z);

#ifdef __cplusplus
}
#endif

#endif /* HAVE_COMPLEX_H */


static inline size_t
pastix_size_of(pastix_coeftype_t type)
{
    switch(type) {
    case PastixFloat:     return   sizeof(float);
    case PastixDouble:    return   sizeof(double);
    case PastixComplex32: return 2*sizeof(float);
    case PastixComplex64: return 2*sizeof(double);
    default:
        fprintf(stderr, "pastix_size_of: invalid type parameter\n");
        assert(0);
        return sizeof(double);
    }
}

/** ****************************************************************************
 * Pastix data structures
 **/

/* Sparse matrix */
struct spmatrix_s;
typedef struct spmatrix_s spmatrix_t;

/* To make it compatible with version with spm inside pastix */
typedef spmatrix_t pastix_spm_t;

/* Main structure of the pastix solver associated to a given problem */
struct pastix_data_s;
typedef struct pastix_data_s pastix_data_t;

/* Graph structure (No values) */
struct pastix_graph_s;
typedef struct pastix_graph_s pastix_graph_t;

/* Ordering structure */
struct pastix_order_s;
typedef struct pastix_order_s pastix_order_t;

/* Solver matrix structure to store L(U)*/
struct solver_matrix_s;
typedef struct solver_matrix_s SolverMatrix;

#endif /* _pastix_datatypes_h_ */
