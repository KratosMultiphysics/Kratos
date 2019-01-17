/**
 *
 * @file codelet_zscalo.c
 *
 * StarPU codelets for blas-like functions
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
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
static struct starpu_perfmodel starpu_blok_zscalo_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "blok_zscalo",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_blok_zscalo_cpu(void *descr[], void *cl_arg)
{
    pastix_trans_t  trans;
    SolverCblk     *cblk;
    pastix_int_t    blok_m;

    const pastix_complex64_t *A;
    const pastix_complex64_t *D;
    pastix_complex64_t *B;

    A = (const pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);
    D = (const pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[1]);
    B = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[2]);

    starpu_codelet_unpack_args( cl_arg, &trans, &cblk, &blok_m );

    assert( cblk->cblktype & CBLK_TASKS_2D );

    cpublok_zscalo( trans, cblk, blok_m, A, D, B );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zscalo, 3 )

void
starpu_task_blok_zscalo( pastix_trans_t    trans,
                         SolverCblk       *cblk,
                         SolverBlok       *blok,
                         sopalin_data_t   *sopalin_data,
                         int               prio )
{
    starpu_data_handle_t *handler = (starpu_data_handle_t*)(blok->handler);
    SolverBlok  *blokA = blok;
    pastix_int_t blok_m = blok - cblk->fblokptr;
    pastix_int_t M = blok_rownbr( blokA );

    while( (blokA < cblk[1].fblokptr) &&
           (blokA[0].fcblknm == blokA[1].fcblknm) &&
           (blokA[0].lcblknm == blokA[1].lcblknm) )
    {
        blokA++;
        M += blok_rownbr( blokA );
    }

    starpu_vector_data_register( handler + 1, -1, (uintptr_t)NULL, M * cblk_colnbr( cblk ),
                                 sopalin_data->solvmtx->starpu_desc->typesze );

    starpu_insert_task(
        pastix_codelet(&cl_blok_zscalo_cpu),
        STARPU_VALUE, &trans,  sizeof(pastix_trans_t),
        STARPU_VALUE, &cblk,   sizeof(SolverCblk*),
        STARPU_VALUE, &blok_m, sizeof(pastix_int_t),
        STARPU_R,      blok->handler[0],
        STARPU_R,      cblk->fblokptr->handler[0],
        STARPU_W,      blok->handler[1],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "blok_zscalo",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * @}
 */
