/**
 *
 * @file codelet_solve_zdiagsp.c
 *
 * StarPU codelets for diag functions.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-13
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup starpu_diag_solve
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
 * Cblk version
 */
static struct starpu_perfmodel starpu_solve_cblk_zdiag_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "solve_cblk_zdiag",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_solve_cblk_zdiag_cpu(void *descr[], void *cl_arg)
{
    SolverCblk         *cblk;
    int                 nrhs;
    pastix_complex64_t *b;
    int                 ldb;

    b    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(descr[1]);
    ldb  = (int)STARPU_MATRIX_GET_LD(descr[1]);
    nrhs = (int)STARPU_MATRIX_GET_NY(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &cblk );

    solve_cblk_zdiag( cblk, nrhs,
                      b, ldb, NULL );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_cblk_zdiag, 2 )
/**
 *******************************************************************************
 *
 * @brief Submit a task to perform a diagonal solve related to one cblk
 * to all the right hand side.
 *
 *******************************************************************************
 *
 * @param[inout] sopalin_data
 *          The structure providing descriptor about b and the A matrix.
 *
 * @param[in] cblk
 *          The cblk structure to which diagonal block belongs to.
 *
 * @param[in] prio
 *          The priority of the task in the DAG.
 *
 *******************************************************************************/
void
starpu_stask_cblk_zdiag( sopalin_data_t *sopalin_data,
                        SolverCblk     *cblk,
                        int             prio )
{
    SolverMatrix *solvmtx = sopalin_data->solvmtx;
    pastix_int_t  cblknum = cblk - solvmtx->cblktab;
    starpu_data_handle_t   handle;

    /* if ( cblk->cblktype & CBLK_TASKS_2D ) { */
    /*     handle = cblk->fblokptr->handler[0]; */
    /* } */
    /* else { */
        handle = cblk->handler[0];
    /* } */

    starpu_insert_task(
        pastix_codelet(&cl_solve_cblk_zdiag_cpu),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_R,      handle,
        STARPU_RW,     solvmtx->starpu_desc_rhs->handletab[cblknum],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "solve_cblk_zdiag",
#endif
        0);
    (void)prio;
}
/**
 * @}
 */
