/**
 *
 * @file codelet_solve_zgemm.c
 *
 * StarPU codelet for gemm function
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "codelets.h"

/**
 * Block version
 */
static struct starpu_perfmodel starpu_solve_blok_zgemm_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "solve_blok_zgemm",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_solve_blok_zgemm_cpu(void *descr[], void *cl_arg)
{
    pastix_coefside_t   coef;
    pastix_side_t       side;
    pastix_trans_t      trans;
    SolverCblk         *cblk;
    SolverBlok         *blok;
    SolverCblk         *fcbk;
    pastix_complex64_t *B, *C;
    pastix_int_t        nrhs, ldb, ldc;

    B    = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    ldb  = (pastix_int_t)        STARPU_MATRIX_GET_LD (descr[1]);
    nrhs = (pastix_int_t)        STARPU_MATRIX_GET_NY (descr[1]);
    C    = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    ldc  = (pastix_int_t)        STARPU_MATRIX_GET_LD (descr[2]);

    starpu_codelet_unpack_args( cl_arg, &coef, &side, &trans,
                                &cblk, &blok, &fcbk );

    solve_blok_zgemm( coef, side, trans, nrhs,
                      cblk, blok, fcbk, B, ldb, C, ldc );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_blok_zgemm, 3 )

/**
 *******************************************************************************
 *
 * @brief Submit a task to perform a gemm.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          Specify whether the computation are made with the L part, or the U
 *          part of A. It has to be either PastixLCoef, or PastixUCoef.
 *
 * @param[in] side
 *          Specify the side parameter of the TRSM.
 *
 * @param[in] trans
 *          Specify the transposition used for the matrix A in the
 *          computation. It has to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure that corresponds to the A and B matrix.
 *
 * @param[in] blok
 *          The blok structure that corresponds to the A matrix, and that
 *          belongs either to cblk or fcbk depending on the side parameter.
 *
 * @param[inout] fcbk
 *          The cblk structure that corresponds to the C matrix.

 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX, and
 *          descriptor of b (providing nrhs, b and ldb).

 * @param[in] prio
 *          The priority of the task in th DAG.
 *
 *******************************************************************************/
void
starpu_stask_blok_zgemm( pastix_coefside_t coef,
                         pastix_side_t     side,
                         pastix_trans_t    trans,
                         const SolverCblk *cblk,
                         const SolverBlok *blok,
                         SolverCblk       *fcbk,
                         sopalin_data_t   *sopalin_data,
                         pastix_int_t      prio )
{
    SolverMatrix          *solvmtx = sopalin_data->solvmtx;
    pastix_int_t           cblknum = cblk - solvmtx->cblktab;
    pastix_int_t           fcbknum = fcbk - solvmtx->cblktab;
    struct starpu_codelet *codelet = &cl_solve_blok_zgemm_cpu;
    starpu_data_handle_t   handle;

    /* if ( cblk->cblktype & CBLK_TASKS_2D ) { */
    /*     handle = blok->handler[coef]; */
    /* } */
    /* else { */
        handle = cblk->handler[coef];
    /* } */

    starpu_insert_task(
        pastix_codelet(codelet),
        STARPU_VALUE, &coef,         sizeof(pastix_coefside_t),
        STARPU_VALUE, &side,         sizeof(pastix_side_t),
        STARPU_VALUE, &trans,        sizeof(pastix_trans_t),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_VALUE, &blok,         sizeof(SolverBlok*),
        STARPU_VALUE, &fcbk,         sizeof(SolverCblk*),
        STARPU_R,      handle,
        STARPU_R,      solvmtx->starpu_desc_rhs->handletab[cblknum],
        STARPU_RW,     solvmtx->starpu_desc_rhs->handletab[fcbknum],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "solve_blok_zgemm",
#endif
        0);
    (void)prio;
}

/**
 * @}
 */
