/**
 * @file z_nan_check.h
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Manage nancheck for lowrank kernels.
 * This header describes all the LAPACKE functions used for low-rank kernels,
 * as well as some macros to manage nancheck.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Gregoire Pichon
 * @date 2018-07-16
 * @precisions normal z -> s d c
 *
 */
#ifndef _z_nan_check_h_
#define _z_nan_check_h_

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// #endif /* DOXYGEN_SHOULD_SKIP_THIS */

#if defined(PASTIX_DEBUG_LR_NANCHECK)
#define LAPACKE_zlacpy_work LAPACKE_zlacpy
#define LAPACKE_zlaset_work LAPACKE_zlaset

#define LAPACKE_zunmlq_work( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_, _w_, _ldw_ ) \
    LAPACKE_zunmlq( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_ )
#define LAPACKE_zunmqr_work( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_, _w_, _ldw_ ) \
    LAPACKE_zunmqr( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_ )

#define LAPACKE_zgeqrf_work( _layout_, _m_, _n_, _a_, _lda_, _tau_, _w_, _ldw_ ) \
    LAPACKE_zgeqrf( _layout_, _m_, _n_, _a_, _lda_, _tau_ )
#define LAPACKE_zgelqf_work( _layout_, _m_, _n_, _a_, _lda_, _tau_, _w_, _ldw_ ) \
    LAPACKE_zgelqf( _layout_, _m_, _n_, _a_, _lda_, _tau_ )

#if defined(PRECISION_z) || defined(PRECISION_c)
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, (double*)(_w_) )
#else
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, (double*)(_w_) )
#endif

#else

#if defined(PRECISION_z) || defined(PRECISION_c)
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ )
#else
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_ )
#endif

#endif /* defined(PASTIX_DEBUG_LR_NANCHECK) */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif /* _z_nan_check_h_ */
