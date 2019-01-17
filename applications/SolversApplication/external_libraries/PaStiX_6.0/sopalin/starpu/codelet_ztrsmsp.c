/**
 *
 * @file codelet_ztrsmsp.c
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
#if defined(PASTIX_WITH_CUDA)
#include "pastix_zcuda.h"
#endif
#include "pastix_starpu.h"
#include "codelets.h"
#include "pastix_starpu_model.h"

/**
 * Block version
 */
static struct starpu_perfmodel starpu_blok_ztrsmsp_model =
{
    .type = STARPU_PER_ARCH,
    .symbol = "blok_ztrsmsp",
    .arch_cost_function = blok_trsmsp_cost,
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_blok_ztrsmsp_cpu(void *descr[], void *cl_arg)
{
    pastix_coefside_t coef;
    pastix_side_t     side;
    pastix_uplo_t     uplo;
    pastix_trans_t    trans;
    pastix_diag_t     diag;
    SolverCblk       *cblk;
    pastix_int_t      blok_m;
    sopalin_data_t   *sopalin_data;

    const pastix_complex64_t *A;
    pastix_complex64_t *C;

    A = (const pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);
    C = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &coef, &side, &uplo, &trans, &diag,
                                &cblk, &blok_m, &sopalin_data );

    assert( cblk->cblktype & CBLK_TASKS_2D );

    cpublok_ztrsmsp( coef, side, uplo, trans, diag,
                     cblk, blok_m, A, C,
                     &(sopalin_data->solvmtx->lowrank) );
}

#if defined(PASTIX_WITH_CUDA)
static void fct_blok_ztrsmsp_gpu(void *descr[], void *cl_arg)
{
    pastix_coefside_t coef;
    pastix_side_t     side;
    pastix_uplo_t     uplo;
    pastix_trans_t    trans;
    pastix_diag_t     diag;
    SolverCblk       *cblk;
    pastix_int_t      blok_m;
    sopalin_data_t   *sopalin_data;

    const cuDoubleComplex *A;
    cuDoubleComplex *C;

    A = (const cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[0]);
    C = (cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &coef, &side, &uplo, &trans, &diag,
                                &cblk, &blok_m, &sopalin_data );

    assert( cblk->cblktype & CBLK_TASKS_2D );

    gpublok_ztrsmsp( coef, side, uplo, trans, diag,
                     cblk, blok_m, A, C,
                     &(sopalin_data->solvmtx->lowrank),
                     starpu_cuda_get_local_stream() );
}
#endif /* defined(PASTIX_WITH_CUDA) */
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_ztrsmsp, 2 )
CODELETS_GPU( blok_ztrsmsp, 2, STARPU_CUDA_ASYNC )

void
starpu_task_blok_ztrsmsp( pastix_coefside_t coef,
                          pastix_side_t     side,
                          pastix_uplo_t     uplo,
                          pastix_trans_t    trans,
                          pastix_diag_t     diag,
                          SolverCblk       *cblk,
                          SolverBlok       *blok,
                          sopalin_data_t   *sopalin_data,
                          int               prio )
{
    pastix_int_t blok_m = blok - cblk->fblokptr;
    struct starpu_codelet *codelet = &cl_blok_ztrsmsp_cpu;

#if defined(PASTIX_WITH_CUDA)
    if ( !(cblk->cblktype & CBLK_COMPRESSED) ) {
        codelet = &cl_blok_ztrsmsp_gpu;
    }
#endif

    starpu_insert_task(
        pastix_codelet(codelet),
        STARPU_VALUE, &coef,         sizeof(pastix_coefside_t),
        STARPU_VALUE, &side,         sizeof(pastix_side_t),
        STARPU_VALUE, &uplo,         sizeof(pastix_uplo_t),
        STARPU_VALUE, &trans,        sizeof(pastix_trans_t),
        STARPU_VALUE, &diag,         sizeof(pastix_diag_t),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_VALUE, &blok_m,       sizeof(pastix_int_t),
        STARPU_VALUE, &sopalin_data, sizeof(sopalin_data_t*),
        STARPU_R,      cblk->fblokptr->handler[coef],
        STARPU_RW,     blok->handler[coef],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "blok_ztrsmsp",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * @}
 */
