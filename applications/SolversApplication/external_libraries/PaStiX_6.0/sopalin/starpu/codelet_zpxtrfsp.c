/**
 *
 * @file codelet_zpxtrfsp.c
 *
 * StarPU codelets for complex LL^t functions
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @precisions normal z -> z c
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
#include "pastix_starpu_model.h"

/**
 * Cblk version
 */
static struct starpu_perfmodel starpu_cblk_zpxtrfsp1d_panel_model =
{
    .type = STARPU_PER_ARCH,
    .symbol = "cblk_zpxtrfsp1d_panel",
    .arch_cost_function = cblk_pxtrf_cost,
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_cblk_zpxtrfsp1d_panel_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
    pastix_complex64_t *L;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);

    starpu_codelet_unpack_args( cl_arg, &sopalin_data, &cblk );

    solvmtx = sopalin_data->solvmtx;
    nbpivot = cpucblk_zpxtrfsp1d_panel( solvmtx, cblk, L );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( cblk_zpxtrfsp1d_panel, 1 )

void
starpu_task_cblk_zpxtrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk,
                                   int             prio )
{
    starpu_insert_task(
        pastix_codelet(&cl_cblk_zpxtrfsp1d_panel_cpu),
        STARPU_VALUE, &sopalin_data, sizeof(sopalin_data_t*),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_RW,     cblk->handler[0],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "cblk_zpxtrfsp1d_panel",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * Blok version
 */
static struct starpu_perfmodel starpu_blok_zpxtrfsp_model =
{
    .type = STARPU_PER_ARCH,
    .symbol = "blok_zpxtrfsp",
    .arch_cost_function = blok_pxtrf_cost,
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_blok_zpxtrfsp_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
    pastix_complex64_t *L;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);

    starpu_codelet_unpack_args( cl_arg, &sopalin_data, &cblk );

    assert(cblk->cblktype & CBLK_TASKS_2D);

    solvmtx = sopalin_data->solvmtx;
    nbpivot = cpucblk_zpxtrfsp1d_pxtrf( solvmtx, cblk, L );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zpxtrfsp, 1 )

void
starpu_task_blok_zpxtrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk,
                         int             prio )
{
    starpu_insert_task(
        pastix_codelet(&cl_blok_zpxtrfsp_cpu),
        STARPU_VALUE, &sopalin_data, sizeof(sopalin_data_t*),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_RW,     cblk->fblokptr->handler[0],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "blok_zpxtrfsp",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * @}
 */
