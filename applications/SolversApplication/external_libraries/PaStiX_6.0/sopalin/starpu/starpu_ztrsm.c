/**
 *
 * @file starpu_ztrsm.c
 *
 * PaStiX ztrsm StarPU wrapper.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 * @precisions normal z -> s d c
 *
 * @addtogroup starpu_trsm_solve
 * @{
 *
 **/

#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"

/**
 *******************************************************************************
 *
 * @brief Apply a forward solve related to one cblk to all the right hand side.
 *        (StarPU version)
 *
 ******************************************************************************* *
 *
 * @param[in] mode
 *          Specify whether the schur complement and interface are applied to
 *          the right-hand-side. It has to be either PastixSolvModeLocal,
 *          PastixSolvModeInterface or PastixSolvModeSchur.
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right
 *          in the equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX, and
 *          descriptor of b (providing nrhs, b and ldb).
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] prio
 *          The priority of the task in th DAG.
 *
 *******************************************************************************/
void
starpu_cblk_ztrsmsp_forward( pastix_solv_mode_t  mode,
                             pastix_side_t       side,
                             pastix_uplo_t       uplo,
                             pastix_trans_t      trans,
                             pastix_diag_t       diag,
                             sopalin_data_t     *sopalin_data,
                             const SolverCblk   *cblk,
                             pastix_int_t        prio )
{
    pastix_coefside_t  cs;
    SolverMatrix      *datacode = sopalin_data->solvmtx;
    SolverCblk        *fcbk;
    SolverBlok        *blok;
    pastix_trans_t     tA;

    if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode != PastixSolvModeSchur) ) {
        return;
    }

    if ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixTrans;
        cs = PastixUCoef;

        /* Right is not handled yet */
        assert( 0 );
    }
    else if ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;

        /* Right is not handled yet */
        assert( 0 );
    }
    else if ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixNoTrans;
        cs = PastixUCoef;

        /* We do not handle conjtrans in complex as we store U^t */
        assert( trans != PastixConjTrans );
    }
    else if ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;
    }
    else {
        /* This correspond to case treated in backward TRSM */
        assert(0);
        return;
    }

    /* In sequential */
    assert( cblk->fcolnum == cblk->lcolidx );

    /* Solve the diagonal block */
    starpu_stask_blok_ztrsm(
        cs, side, PastixLower, tA, diag, cblk, sopalin_data, prio );

    /* Apply the update */
    for (blok = cblk[0].fblokptr+1; blok < cblk[1].fblokptr; blok++ ) {
        fcbk  = datacode->cblktab + blok->fcblknm;

        if ( (fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeLocal) ) {
            return;
        }

        starpu_stask_blok_zgemm( cs, PastixLeft, tA,
                                 cblk, blok, fcbk, sopalin_data, prio );
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply a backward solve related to one cblk to all the right hand side.
 *        (StarPU version)
 *
 *******************************************************************************
 *
 * @param[in] mode
 *          Specify whether the schur complement and interface are applied to
 *          the right-hand-side. It has to be either PastixSolvModeLocal,
 *          PastixSolvModeInterface or PastixSolvModeSchur.
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right
 *          in the equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX, and
 *          descriptor of b (providing nrhs, b and ldb).
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] prio
 *          The priority of the task in th DAG.
 *
 *******************************************************************************/
void
starpu_cblk_ztrsmsp_backward( pastix_solv_mode_t  mode,
                              pastix_side_t       side,
                              pastix_uplo_t       uplo,
                              pastix_trans_t      trans,
                              pastix_diag_t       diag,
                              sopalin_data_t     *sopalin_data,
                              const SolverCblk   *cblk,
                              pastix_int_t        prio )
{
    pastix_coefside_t  cs;
    SolverMatrix      *datacode = sopalin_data->solvmtx;
    SolverCblk        *fcbk;
    SolverBlok        *blok;
    pastix_int_t       j;
    pastix_trans_t     tA;

    /*
     *  Left / Upper / NoTrans (Backward)
     */
    if ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixTrans;
        cs = PastixUCoef;
    }
    else if ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;
    }
    else if ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixNoTrans;
        cs = PastixUCoef;

        /* Right is not handled yet */
        assert( 0 );

        /* We do not handle conjtrans in complex as we store U^t */
        assert( trans != PastixConjTrans );
    }
    else if ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;

        /* Right is not handled yet */
        assert( 0 );
    }
    else {
        /* This correspond to case treated in forward TRSM */
        assert(0);
        return;
    }

    if ( !(cblk->cblktype & CBLK_IN_SCHUR) || (mode == PastixSolvModeSchur) ) {
        /* Solve the diagonal block */
        starpu_stask_blok_ztrsm(
            cs, side, PastixLower, tA, diag, cblk, sopalin_data, prio );
    }

    /* Apply the update */
    for (j   = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
        blok = datacode->bloktab + datacode->browtab[j];
        fcbk = datacode->cblktab + blok->lcblknm;

        if ( (fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeInterface) ) {
            continue;
        }

        starpu_stask_blok_zgemm( cs, PastixRight, tA,
                                 cblk, blok, fcbk, sopalin_data, prio );
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply a TRSM on a problem with 1 dimension (StarPU version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The data that provide the mode.
 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX., and
 *          descriptor of b (providing nrhs, b and ldb).
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right
 *          in the equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 *******************************************************************************/
void
starpu_ztrsm_sp1dplus( pastix_data_t      *pastix_data,
                       sopalin_data_t     *sopalin_data,
                       int                 side,
                       int                 uplo,
                       int                 trans,
                       int                 diag )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk *cblk;
    pastix_int_t i, cblknbr;
    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

    /* Backward like */
    if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
    {
        cblknbr = (mode == PastixSolvModeLocal) ? datacode->cblkschur : datacode->cblknbr;

        cblk    = datacode->cblktab + cblknbr - 1;
        for (i=0; i<cblknbr; i++, cblk--){
            starpu_cblk_ztrsmsp_backward( mode, side, uplo, trans, diag,
                                          sopalin_data, cblk, cblknbr - i );
        }
    }
    /* Forward like */
    else
        /*
         * ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         * ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) )
         */
    {
        cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;

        cblk    = datacode->cblktab;
        for (i=0; i<cblknbr; i++, cblk++){
            starpu_cblk_ztrsmsp_forward( mode, side, uplo, trans, diag,
                                         sopalin_data, cblk, cblknbr - i );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply the TRSM solve (StarPU version).
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide informations about starpu and the schur solving mode.
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right
 *          in the equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX, and
 *          descriptor of b (providing nrhs, b and ldb).
 *
 * @param[in] nrhs
 *          The number of right hand side.
 *
 * @param[inout] b
 *          The pointer to vectors of the right hand side.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************/
void
starpu_ztrsm( pastix_data_t      *pastix_data,
              int                 side,
              int                 uplo,
              int                 trans,
              int                 diag,
              sopalin_data_t     *sopalin_data,
              int                 nrhs,
              pastix_complex64_t *b,
              int                 ldb )
{
    starpu_sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->starpu_desc;
    starpu_dense_matrix_desc_t  *ddesc;

    /*
     * Start StarPU if not already started
     */
    if (pastix_data->starpu == NULL) {
        int argc = 0;
        pastix_starpu_init( pastix_data, &argc, NULL, NULL );
    }

    if ( sdesc == NULL ) {
        /* Create the sparse matrix descriptor */
        starpu_sparse_matrix_init( sopalin_data->solvmtx,
                                   sizeof( pastix_complex64_t ), PastixHermitian,
                                   1, 0 );
        sdesc = sopalin_data->solvmtx->starpu_desc;
    }

    /* Create the dense matrix descriptor */
    starpu_dense_matrix_init( sopalin_data->solvmtx,
                              nrhs, (char*)b, ldb,
                              sizeof(pastix_complex64_t), 1, 0 );
    ddesc = sopalin_data->solvmtx->starpu_desc_rhs;

    starpu_resume();
    starpu_ztrsm_sp1dplus( pastix_data, sopalin_data, side, uplo, trans, diag );

    starpu_sparse_matrix_getoncpu( sdesc );
    starpu_dense_matrix_getoncpu( ddesc );
    starpu_task_wait_for_all();
#if defined(PASTIX_WITH_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#endif
    starpu_pause();

    return;
}

/**
 *@}
 */
