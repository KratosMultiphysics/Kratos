/**
 *
 * @file codelet_solve_ztrsm.c
 *
 * StarPU codelet for TRSM function
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
 * @addtogroup starpu_trsm_solve
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
static struct starpu_perfmodel starpu_solve_blok_ztrsm_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "solve_blok_ztrsm",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_solve_blok_ztrsm_cpu(void *descr[], void *cl_arg)
{
    pastix_coefside_t   coef;
    pastix_side_t       side;
    pastix_uplo_t       uplo;
    pastix_trans_t      trans;
    pastix_diag_t       diag;
    SolverCblk         *cblk;
    pastix_complex64_t *B;
    pastix_int_t        nrhs, ldb;

    B    = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    ldb  = (pastix_int_t)        STARPU_MATRIX_GET_LD (descr[1]);
    nrhs = (pastix_int_t)        STARPU_MATRIX_GET_NY (descr[1]);

    starpu_codelet_unpack_args( cl_arg, &coef, &side, &uplo, &trans, &diag, &cblk );

    solve_blok_ztrsm( coef, side, uplo, trans, diag,
                      cblk, nrhs, B, ldb );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_blok_ztrsm, 2 )

/**
 *******************************************************************************
 *
 * @brief Submit a task to do a trsm related to a diagonal
 * block of the matrix A.
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
 * @param[in] uplo
 *          Specify the uplo parameter of the TRSM.
 *
 * @param[in] trans
 *          Specify the transposition used for the matrix A in the
 *          computation. It has to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure that corresponds to the A and B matrix.
 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX, and
 *          descriptor of b (providing nrhs, b and ldb).

 * @param[in] prio
 *          The priority of the task in th DAG.
 *
 *******************************************************************************/
void
starpu_stask_blok_ztrsm( pastix_coefside_t coef,
                         pastix_side_t     side,
                         pastix_uplo_t     uplo,
                         pastix_trans_t    trans,
                         pastix_diag_t     diag,
                         const SolverCblk *cblk,
                         sopalin_data_t   *sopalin_data,
                         pastix_int_t      prio )
{
    SolverMatrix          *solvmtx = sopalin_data->solvmtx;
    pastix_int_t           cblknum = cblk - solvmtx->cblktab;
    struct starpu_codelet *codelet = &cl_solve_blok_ztrsm_cpu;
    starpu_data_handle_t   handle;

    /* if ( cblk->cblktype & CBLK_TASKS_2D ) { */
    /*     handle = cblk->fblokptr->handler[coef]; */
    /* } */
    /* else { */
        handle = cblk->handler[coef];
    /* } */

    starpu_insert_task(
        pastix_codelet(codelet),
        STARPU_VALUE, &coef,         sizeof(pastix_coefside_t),
        STARPU_VALUE, &side,         sizeof(pastix_side_t),
        STARPU_VALUE, &uplo,         sizeof(pastix_uplo_t),
        STARPU_VALUE, &trans,        sizeof(pastix_trans_t),
        STARPU_VALUE, &diag,         sizeof(pastix_diag_t),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_R,      handle,
        STARPU_RW,     solvmtx->starpu_desc_rhs->handletab[cblknum],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "solve_blok_ztrsm",
#endif
        0);
    (void)prio;
}

/**
 * @}
 */
