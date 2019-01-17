/**
 * @file pastix_zlrcores.h
 *
 * PaStiX kernel header.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 */
#ifndef _pastix_zlrcores_h_
#define _pastix_zlrcores_h_

#include "pastix_lowrank.h"

/**
 *
 * @addtogroup kernel_lr
 * @{
 *    This module contains all the low-rank kernels working on pastix_lr_t
 *    matrix representations.
 *
 *    @name PastixComplex64 low-rank kernels
 *    @{
 */
void core_zlralloc( pastix_int_t M, pastix_int_t N, pastix_int_t rkmax, pastix_lrblock_t *A );
void core_zlrfree ( pastix_lrblock_t *A );
int  core_zlrsze  ( int copy, pastix_int_t M, pastix_int_t N, pastix_lrblock_t *A, pastix_int_t newrk, pastix_int_t newrkmax, pastix_int_t rklimit );
int  core_zlr2ge  ( pastix_trans_t trans, pastix_int_t M, pastix_int_t N, const pastix_lrblock_t *Alr, pastix_complex64_t *A, pastix_int_t lda );

void core_zlrcpy  ( const pastix_lr_t *lowrank,
                    pastix_trans_t transA, pastix_complex64_t alpha,
                    pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                    pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                    pastix_int_t offx, pastix_int_t offy );

void core_zlrconcatenate_u( pastix_complex64_t alpha,
                            pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                            pastix_int_t M2,                        pastix_lrblock_t *B,
                            pastix_int_t offx,
                            pastix_complex64_t *u1u2 );
void core_zlrconcatenate_v( pastix_trans_t transA1, pastix_complex64_t alpha,
                            pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                             pastix_int_t N2,       pastix_lrblock_t *B,
                            pastix_int_t offy,
                            pastix_complex64_t *v1v2 );
/**
 *     @}
 * @}
 *
 * @addtogroup kernel_lr_lrmm
 * @{
 *    This are the kernels to compute the low-rank updates
 *
 *    @name PastixComplex64 LRMM low-rank kernels
 *    @{
 */

/**
 * @brief Structure to store all the parameters of the core_zlrmm family functions
 */
typedef struct core_zlrmm_s {
    const pastix_lr_t      *lowrank;     /**< The lowrank structure                                                 */
    pastix_trans_t          transA;      /**< Specify op(A) and is equal to PastixNoTrans, PastixTrans, or PastixConjTrans */
    pastix_trans_t          transB;      /**< Specify op(B) and is equal to PastixNoTrans, PastixTrans, or PastixConjTrans */
    pastix_int_t            M;           /**< Number of rows     of the A matrix                                    */
    pastix_int_t            N;           /**< Number of columns  of the B matrix                                    */
    pastix_int_t            K;           /**< Number of columns  of the A matrix (= number of rows of the B matrix) */
    pastix_int_t            Cm;          /**< Number of rows     of the C matrix that receives the AB contribution  */
    pastix_int_t            Cn;          /**< Number of columns  of the C matrix that receives the AB contribution  */
    pastix_int_t            offx;        /**< Horizontal offsets of the AB product in the C matrix                  */
    pastix_int_t            offy;        /**< Vertical   offsets of the AB product in the C matrix                  */
    pastix_complex64_t      alpha;       /**< The alpha factor                                                      */
    const pastix_lrblock_t *A;           /**< The A matrix described in a low-rank structure                        */
    const pastix_lrblock_t *B;           /**< The B matrix described in a low-rank structure                        */
    pastix_complex64_t      beta;        /**< The beta factor                                                       */
    pastix_lrblock_t       *C;           /**< The C matrix described in a low-rank structure                        */
    pastix_complex64_t     *work;        /**< The pointer to an available workspace                                 */
    pastix_int_t            lwork;       /**< The size of the given workspace                                       */
    pastix_int_t            lwused;      /**< The size of the workspace that is already used                        */
    pastix_atomic_lock_t   *lock;        /**< The lock to protect the concurrent accesses on the C matrix           */
} core_zlrmm_t;

/**
 * @brief Initialize all the parameters of the core_zlrmm family functions to ease the access
 */
#define PASTE_CORE_ZLRMM_PARAMS(_a_)                   \
    const pastix_lr_t      *lowrank = (_a_)->lowrank;  \
    pastix_trans_t          transA  = (_a_)->transA;   \
    pastix_trans_t          transB  = (_a_)->transB;   \
    pastix_int_t            M       = (_a_)->M;        \
    pastix_int_t            N       = (_a_)->N;        \
    pastix_int_t            K       = (_a_)->K;        \
    pastix_int_t            Cm      = (_a_)->Cm;       \
    pastix_int_t            Cn      = (_a_)->Cn;       \
    pastix_int_t            offx    = (_a_)->offx;     \
    pastix_int_t            offy    = (_a_)->offy;     \
    pastix_complex64_t      alpha   = (_a_)->alpha;    \
    const pastix_lrblock_t *A       = (_a_)->A;        \
    const pastix_lrblock_t *B       = (_a_)->B;        \
    pastix_complex64_t      beta    = (_a_)->beta;     \
    pastix_lrblock_t       *C       = (_a_)->C;        \
    pastix_complex64_t     *work    = (_a_)->work;     \
    pastix_int_t            lwork   = (_a_)->lwork;    \
    pastix_atomic_lock_t   *lock    = (_a_)->lock;

/**
 * @brief Void all the parameters of the core_zlrmm family functions to silent warnings
 */
#define PASTE_CORE_ZLRMM_VOID                   \
    (void)lowrank;                              \
    (void)transA;                               \
    (void)transB;                               \
    (void)M;                                    \
    (void)N;                                    \
    (void)K;                                    \
    (void)Cm;                                   \
    (void)Cn;                                   \
    (void)offx;                                 \
    (void)offy;                                 \
    (void)alpha;                                \
    (void)A;                                    \
    (void)B;                                    \
    (void)beta;                                 \
    (void)C;                                    \
    (void)work;                                 \
    (void)lwork;                                \
    (void)lock

/**
 * @brief Function to get a workspace pointer if space is available in the one provided
 * @param[inout] params  The parameters structure for core_zlrmm family functions
 * @param[in]    newsize The required workspace size in number of elements
 * @return The pointer to the workspace if enough space available, NULL otherwise.
 */
static inline pastix_complex64_t *
core_zlrmm_getws( core_zlrmm_t *params,
                  ssize_t newsize )
{
    pastix_complex64_t *work = NULL;
    if ( (params->lwused + newsize) <= params->lwork )
    {
        work = params->work + params->lwused;
        params->lwused += newsize;
    }
    /* else */
    /* { */
    /*     if ( (params->work == NULL) || (params->lwused == 0) ) */
    /*     { */
    /*         params->work = realloc( params->work, newsize * sizeof(pastix_complex64_t) ); */
    /*         params->lwork  = newsize; */
    /*         params->lwused = newsize; */
    /*         work = params->work; */
    /*     } */
    /* } */
    return work;
}

/**
 *      @name update_fr Functions to perform the update on a full-rank matrix
 *      @{
 */
pastix_fixdbl_t core_zfrfr2fr( core_zlrmm_t *params );
pastix_fixdbl_t core_zfrlr2fr( core_zlrmm_t *params );
pastix_fixdbl_t core_zlrfr2fr( core_zlrmm_t *params );
pastix_fixdbl_t core_zlrlr2fr( core_zlrmm_t *params );

/**
 *      @}
 *      @name update_lr Functions to prepare the AB product for an update on a low-rank matrix
 *      @{
 */
pastix_fixdbl_t core_zfrfr2lr( core_zlrmm_t     *params,
                               pastix_lrblock_t *AB,
                               int              *infomask,
                               pastix_int_t      Kmax );
pastix_fixdbl_t core_zfrlr2lr( core_zlrmm_t     *params,
                               pastix_lrblock_t *AB,
                               int              *infomask,
                               pastix_int_t      Brkmin );
pastix_fixdbl_t core_zlrfr2lr( core_zlrmm_t     *params,
                               pastix_lrblock_t *AB,
                               int              *infomask,
                               pastix_int_t      Arkmin );
pastix_fixdbl_t core_zlrlr2lr( core_zlrmm_t     *params,
                               pastix_lrblock_t *AB,
                               int              *infomask );

/**
 *      @}
 *      @name add_lr Functions to add the AB contribution in a low-rank format to any C matrix
 *      @{
 */
pastix_fixdbl_t core_zlr2fr  ( core_zlrmm_t           *params,
                               const pastix_lrblock_t *AB,
                               pastix_trans_t          transV );
pastix_fixdbl_t core_zlr2lr  ( core_zlrmm_t           *params,
                               const pastix_lrblock_t *AB,
                               pastix_trans_t          transV );
pastix_fixdbl_t core_zlr2null( core_zlrmm_t           *params,
                               const pastix_lrblock_t *AB,
                               pastix_trans_t          transV,
                               int                     infomask );

/**
 *       @}
 */

pastix_fixdbl_t core_zlrmm( core_zlrmm_t *params );

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_lr_svd
 * @{
 *    This is the SVD implementation of the low-rank kernels based on the LAPACK
 *    GESVD function.
 *
 *    @name PastixComplex64 SVD low-rank kernels
 *    @{
 */

pastix_fixdbl_t core_zge2lr_svd( pastix_fixdbl_t tol, pastix_int_t rklimit, pastix_int_t m, pastix_int_t n,
                                 const void *Avoid, pastix_int_t lda, pastix_lrblock_t *Alr );
pastix_fixdbl_t core_zrradd_svd( const pastix_lr_t *lowrank, pastix_trans_t transA1, const void *alphaptr,
                                 pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                 pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                                 pastix_int_t offx, pastix_int_t offy);

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_lr_rrqr
 * @{
 *    This is the rank-revealing QR implementation of the low-rank kernels based
 *    on the modified LAPACK GEQP3 function.
 *
 *    @name PastixComplex64 RRQR low-rank kernels
 *    @{
 */

pastix_fixdbl_t core_zge2lr_rrqr( pastix_fixdbl_t tol, pastix_int_t rklimit, pastix_int_t m, pastix_int_t n,
                                  const void *Avoid, pastix_int_t lda, pastix_lrblock_t *Alr );
pastix_fixdbl_t core_zrradd_rrqr( const pastix_lr_t *lowrank, pastix_trans_t transA1, const void *alphaptr,
                                  pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                  pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                                  pastix_int_t offx, pastix_int_t offy);

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_lr_debug
 * @{
 *    This is the debug routines for the low rank kernels.
 *
 *    @name PastixComplex64 low-rank debug functions
 *    @{
 */

void core_zlrdbg_printsvd( pastix_int_t              M,
                           pastix_int_t              N,
                           const pastix_complex64_t *A,
                           pastix_int_t              lda );

int  core_zlrdbg_check_orthogonality( pastix_int_t              M,
                                      pastix_int_t              N,
                                      const pastix_complex64_t *A,
                                      pastix_int_t              lda );

int  core_zlrdbg_check_orthogonality_AB( pastix_int_t M, pastix_int_t NA, pastix_int_t NB,
                                         const pastix_complex64_t *A, pastix_int_t lda,
                                         const pastix_complex64_t *B, pastix_int_t ldb );

/**
 *    @}
 * @}
 *
 */
#endif /* _pastix_zlrcores_h_ */
