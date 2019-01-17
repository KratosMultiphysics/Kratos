/**
 *
 * @file gpu_zgemmsp.c
 *
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX GPU kernel routines
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2018-07-16
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "blend/solver.h"
#include "kernels_trace.h"
#include "pastix_zcuda.h"
#include "pastix_cuda.h"
#include <cublas.h>

static char transstr[3] = { 'N', 'T', 'C' };
static char sidestr[2] = { 'L', 'R' };
static char uplostr[3] = { 'U', 'L', 'A' };
static char diagstr[2] = { 'N', 'U' };

void
gpu_zgemmsp_fermi( const SolverMatrix *solvmatr,
                   pastix_uplo_t uplo, pastix_trans_t trans,
                   int *blocktab,
                   const SolverCblk      *cblk,
                   const SolverBlok      *blok,
                         SolverCblk      *fcblk,
                   const cuDoubleComplex *A,
                   const cuDoubleComplex *B,
                         cuDoubleComplex *C,
                         cudaStream_t stream )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex mzone = make_cuDoubleComplex(-1.0, 0.0);
    cuDoubleComplex zone  = make_cuDoubleComplex( 1.0, 0.0);
#else
    double mzone = -1.0;
    double zone  =  1.0;
#endif

    pastix_int_t lda, ldb, ldc;
    pastix_int_t M, N, K, s;
    pastix_int_t blocknbr, fblocknbr;

    s = (uplo == PastixUpper) ? 1 : 0;

    M = cblk->stride - blok[s].coefind;
    N = blok_rownbr( blok );
    K = cblk_colnbr( cblk );

    lda = cblk->stride;
    ldb = lda;
    ldc = fcblk->stride;

    blocknbr  = cblk[1].fblokptr - blok;
    blocknbr -= s;
    fblocknbr = fcblk[1].fblokptr - fcblk[0].fblokptr;

    C = C + ldc * ( blok->frownum - fcblk->fcolnum );

    pastix_fermi_zgemmsp( 'N', transstr[trans - PastixNoTrans], M, N, K,
                          mzone, A + blok[s].coefind, lda,
                                 B + blok[0].coefind, ldb,
                          zone,  C, ldc,
                          blocknbr,  blocktab + 2 * ((blok+s)        - solvmatr->bloktab),
                          fblocknbr, blocktab + 2 * (fcblk->fblokptr - solvmatr->bloktab),
                          stream );
}

/**
 *******************************************************************************
 *
 * @brief Compute the updates associated to one off-diagonal block on a GPU.
 *
 *******************************************************************************
 *
 * @param[in] sideA
 *          Specify if A and C belong to the lower part, or to the upper part.
 *          If sideA == PastixLCoef, the contribution of:
 *          (block .. (cblk[1].fblokptr-1)) -by- block is computed and added to
 *          C, otherwise the contribution:
 *          (block+1 .. (cblk[1].fblokptr-1)) -by- block is computed and added
 *          to C.
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] sideB
 *          Specify if B belongs to the L part, or to the U part. this is used
 *          internally in the kernel to select the correct data pointer.
 *          If PastixLCoef, B belongs to the L part, otherwise B belogns to the
 *          U part.
 *
 * @param[in] trans
 *          Specify the transposition used for the B matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok
 *          The block from which we compute the contributions.
 *
 * @param[inout] fcblk
 *          The pointer to the data structure that describes the panel on which
 *          we compute the contributions. The C pointer must be one of the
 *          coeftab from this fcblk. Next column blok must be accessible through
 *          fcblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] stream
 *          The CUDA stream that will execute the kernel.
 *
 *******************************************************************************/
void
gpucblk_zgemmsp(       pastix_coefside_t  sideA,
                       pastix_coefside_t  sideB,
                       pastix_trans_t     trans,
                 const SolverCblk        *cblk,
                 const SolverBlok        *blok,
                       SolverCblk        *fcblk,
                 const cuDoubleComplex   *A,
                 const cuDoubleComplex   *B,
                       cuDoubleComplex   *C,
                 const pastix_lr_t       *lowrank,
                       cudaStream_t       stream )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex mzone = make_cuDoubleComplex(-1.0, 0.0);
    cuDoubleComplex zone  = make_cuDoubleComplex( 1.0, 0.0);
#else
    double mzone = -1.0;
    double zone  =  1.0;
#endif
    gemm_params_t params;
    const SolverBlok *iterblok;
    const SolverBlok *fblok;
    const SolverBlok *lblok;

    pastix_int_t stride, stridef, indblok;
    pastix_int_t N, K, max_m = 0;
    int i, shift, count, ldb;

    pastix_fixdbl_t time = kernel_trace_start( PastixKernelGEMMCblk2d2d );

    assert( !(cblk->cblktype  & CBLK_COMPRESSED) );
    assert( !(fcblk->cblktype & CBLK_COMPRESSED) );

    shift = (sideA == PastixUCoef) ? 1 : 0;

    stride  = cblk->stride;
    stridef = fcblk->stride;
    K = cblk_colnbr( cblk );

    /* First blok */
    indblok = blok->coefind;

    N = blok_rownbr( blok );

    /* Move B to the right pointer */
    B = B + indblok;
    ldb = (cblk->cblktype & CBLK_LAYOUT_2D) ? N : stride;

    /* Get the first block of the distant panel */
    fblok = fcblk->fblokptr;

    /* Get the last block to stop the iteration */
    lblok = cblk[1].fblokptr;
    count = (lblok - blok) - shift;

    for (iterblok=blok+shift, i=0; iterblok<lblok; iterblok++, i++) {
        /* Find facing blok */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        stridef = (fcblk->cblktype  & CBLK_LAYOUT_2D) ? blok_rownbr( fblok ) : stridef;
        params.p[i].M    = blok_rownbr( iterblok );
        params.p[i].Aptr = A + iterblok->coefind;
        params.p[i].lda  = (cblk->cblktype  & CBLK_LAYOUT_2D) ? params.p[i].M : stride;
        params.p[i].Cptr = C +
            fblok->coefind + iterblok->frownum - fblok->frownum +
            (blok->frownum - fcblk->fcolnum) * stridef;
        params.p[i].ldc  = stridef;

        max_m = pastix_imax( max_m, params.p[i].M);

        if (i+1 == MAX_BATCH_COUNT) {
            pastix_zgemm_vbatched_nt(
                trans, N, K,
                /* alpha  */  mzone,
                /* B      */  B, ldb,
                /* beta   */  zone,
                max_m, MAX_BATCH_COUNT,
                stream, params );

            /* Restart the loop */
            i = -1;
            count -= MAX_BATCH_COUNT;
            max_m = 0;
        }
    }

    if (count > 0) {
        pastix_zgemm_vbatched_nt(
            trans, N, K,
            /* alpha  */  mzone,
            /* B      */  B, ldb,
            /* beta   */  zone,
            max_m, count,
            stream, params );
    }

#if defined(PASTIX_GENERATE_MODEL)
    cudaStreamSynchronize( stream );
#endif
    {
        pastix_int_t k = cblk_colnbr( cblk );
        pastix_int_t n = blok_rownbr( blok );
        pastix_int_t m = cblk->stride;

        m -= (cblk->cblktype & CBLK_LAYOUT_2D) ? blok->coefind / k : blok->coefind;
        m -= (sideA == PastixUCoef) ? blok_rownbr( blok ) : 0;

        kernel_trace_stop( PastixKernelGEMMCblk2d2d, m, n, k, FLOPS_ZGEMM( m, n, k ), time );
    }
    (void)sideB; (void)lowrank; (void)time;
}

/**
 *******************************************************************************
 *
 * @brief Compute the updates associated to one off-diagonal block on a GPU.
 *
 *******************************************************************************
 *
 *    C_l = C_l - A_l * op(B_s), with B_s = B_l, or B_u
 *  or
 *    C_u = C_u - A_u * op(B_s), with B_s = B_l, or B_u
 *
 *******************************************************************************
 *
 * @param[in] sideA
 *          Specify if A and C belong to the L part, or to the U part of the
 *          matrix. This is used internally in the kernels to select the correct
 *          data pointers.  If PastixLCoef, A and C belong to the L part,
 *          otherwise A and C belong to the U part.
 *
 * @param[in] sideB
 *          Specify if B belongs to the lower or upper part of the matrix. This
 *          is used internally in the kernels to select the correct data
 *          pointers.  If PastixLCoef, B belongs to the L part, otherwise B
 *          belongs to the U part.
 *
 * @param[in] transB
 *          Specify wheter B should be used as PastixNoTrans, PastixTrans, or
 *          PastixConjTrans in the computations.
 *
 * @param[in] cblk
 *          The cblk structure to which block A and B belong to. The A and B
 *          pointers must be one of the [lu]coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[inout] fcblk
 *          The pointer to the data structure that describes the panel on which
 *          we compute the contributions. The C pointer must be one of the
 *          [lu]coeftab from this fcblk.
 *          Next column blok must be accessible through fcblk[1].
 *
 * @param[in] blok_mk
 *          Specify the index of the A block in the cblk column. This index is
 *          0-based for the diagonal block.
 *
 * @param[in] blok_nk
 *          Specify the index of the B block in the cblk column. This index is
 *          0-based for the diagonal block.
 *
 * @param[in] blok_mn
 *          Specify the index of the C block in the fcblk column. This index is
 *          0-based for the diagonal block.
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] lowrank
 *          The structure with the low-rank parameters.
 *
 * @param[in] stream
 *          The CUDA stream that will execute the kernel.
 *
 *******************************************************************************/
void
gpublok_zgemmsp(       pastix_coefside_t  sideA,
                       pastix_coefside_t  sideB,
                       pastix_trans_t     trans,
                 const SolverCblk        *cblk,
                       SolverCblk        *fcblk,
                       pastix_int_t       blok_mk,
                       pastix_int_t       blok_nk,
                       pastix_int_t       blok_mn,
                 const cuDoubleComplex   *A,
                 const cuDoubleComplex   *B,
                       cuDoubleComplex   *C,
                 const pastix_lr_t       *lowrank,
                       cudaStream_t       stream )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex mzone = make_cuDoubleComplex(-1.0, 0.0);
    cuDoubleComplex zone  = make_cuDoubleComplex( 1.0, 0.0);
#else
    double mzone = -1.0;
    double zone  =  1.0;
#endif
    const SolverBlok *blokA, *blokB, *blokC;
    const SolverBlok *bA, *bB, *bC;
    const SolverBlok *fblokK, *lblokK;
    const SolverBlok *fblokN, *lblokN;

    const cuDoubleComplex *Aptr, *Bptr;
    cuDoubleComplex *Cptr;
    pastix_int_t M, N, K, lda, ldb, ldc, cblk_n, cblk_m;
    pastix_int_t full_m;
    size_t offsetA, offsetB, offsetC;

    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t time = kernel_trace_start( PastixKernelGEMMBlok2d2d );

    /* Both cblk and fcblk must be stored in 2D */
    assert( cblk->cblktype  & CBLK_LAYOUT_2D );
    assert( fcblk->cblktype & CBLK_LAYOUT_2D );

    /**
     * Blocs on column K
     */
    fblokK = cblk[0].fblokptr;
    lblokK = cblk[1].fblokptr;

    blokB = fblokK + blok_nk;
    offsetB = blokB->coefind;
    cblk_n = blokB->fcblknm;

    blokA = fblokK + blok_mk;
    offsetA = blokA->coefind;
    cblk_m = blokA->fcblknm;

    /**
     * Blocs on column N
     */
    fblokN = fcblk[0].fblokptr;
    lblokN = fcblk[1].fblokptr;

    blokC = fblokN + blok_mn;
    offsetC = blokC->coefind;
    assert( blokC->fcblknm == cblk_m );

    K = cblk_colnbr( cblk );
    full_m = 0;

    cublasSetKernelStream( stream );
    bC = blokC;
    for (bA = blokA; (bA < lblokK) && (bA->fcblknm == cblk_m); bA++) {
        M = blok_rownbr(bA);
        Aptr = A + bA->coefind - offsetA;
        lda = M;
        full_m += M;

        /* Find facing C blok */
        while (!is_block_inside_fblock( bA, bC )) {
            bC++;
            assert( bC < lblokN );
        }

        Cptr = C + bC->coefind - offsetC;
        ldc = blok_rownbr(bC);

        for (bB = blokB; (bB < lblokK) && (bB->fcblknm == cblk_n); bB++) {
            N = blok_rownbr( bB );
            Bptr = B + bB->coefind - offsetB;
            ldb = N;

            cublasZgemm( 'N', transstr[trans - PastixNoTrans],
                         M, N, K,
                         mzone, Aptr, lda,
                                Bptr, ldb,
                          zone, Cptr + (bA->frownum - bC->frownum)
                                     + (bB->frownum - fcblk->fcolnum) * ldc, ldc );

            flops += FLOPS_ZGEMM( M, N, K );
        }
    }

#if defined(PASTIX_GENERATE_MODEL)
    cudaStreamSynchronize( stream );
#endif
    kernel_trace_stop( PastixKernelGEMMBlok2d2d,
                       full_m, full_m, K, flops, time );

    (void)lblokN; (void)sideA; (void)sideB; (void)lowrank; (void)time;
}

/**
 *******************************************************************************
 *
 * @brief Compute the solve update of a block in a panel.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          - PastixLCoef, use the lower part of the off-diagonal blocks.
 *          - PastixUCoef, use the upper part of the off-diagonal blocks
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************/
void
gpublok_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                 pastix_trans_t trans, pastix_diag_t diag,
                       SolverCblk      *cblk,
                       pastix_int_t     blok_m,
                 const cuDoubleComplex *A,
                       cuDoubleComplex *C,
                 const pastix_lr_t     *lowrank,
                       cudaStream_t     stream )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zone  = make_cuDoubleComplex( 1.0, 0.0);
#else
    double zone  =  1.0;
#endif
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc, offset, cblk_m, full_m;
    cuDoubleComplex *Cptr;
    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t time = kernel_trace_start( PastixKernelTRSMBlok2d );

    assert( !(cblk->cblktype & CBLK_COMPRESSED));

    /* if ( cblk->cblktype & CBLK_COMPRESSED ) { */
    /*     core_ztrsmsp_lrsub( coef, side, uplo, trans, diag, */
    /*                         cblk, blok_m, lowrank ); */
    /* } */
    /* else { */

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    blok   = fblok + blok_m;
    offset = blok->coefind;
    cblk_m = blok->fcblknm;
    full_m = 0;

    cublasSetKernelStream( stream );
    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {

        Cptr = C + blok->coefind - offset;
        M   = blok_rownbr(blok);
        ldc = M;

        cublasZtrsm( sidestr[side - PastixLeft],
                     uplostr[uplo - PastixUpper],
                     transstr[trans - PastixNoTrans],
                     diagstr[diag - PastixNonUnit],
                     M, N, zone,
                     A, lda,
                     Cptr, ldc );
        full_m += M;
        flops += FLOPS_ZTRSM( side, M, N );
    }

    /* } */

#if defined(PASTIX_GENERATE_MODEL)
    cudaStreamSynchronize( stream );
#endif
    kernel_trace_stop( PastixKernelTRSMBlok2d,
                       full_m, N, 0, flops, time );

    (void)lowrank; (void)coef;
}
