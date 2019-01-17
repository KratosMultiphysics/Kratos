/**
 *
 * @file pastix_starpu_model.c
 *
 * Model function for StarPU codelets.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2018-07-16
 *
 * @{
 *
 **/

#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "codelets.h"
#if !defined(PASTIX_WITH_STARPU)
#error "This file should not be compiled if Starpu is not enabled"
#endif
#include "pastix_starpu_model.h"

double
cblk_gemmsp_cost( struct starpu_task           *task,
                  struct starpu_perfmodel_arch *arch,
                  unsigned                      nimpl )
{
    pastix_coefside_t sideA;
    pastix_coefside_t sideB;
    pastix_trans_t    trans;
    SolverCblk       *cblk;
    SolverBlok       *blok;
    SolverCblk       *fcblk;
    sopalin_data_t   *sopalin_data;
    pastix_int_t shift, M, N, K;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sideA, &sideB, &trans,
                                &cblk, &blok, &fcblk, &sopalin_data );

    shift = (sideA == PastixUCoef) ? 1 : 0;
    K = cblk_colnbr( cblk );
    N = blok_rownbr( blok );
    M = cblk->stride - (cblk->cblktype & CBLK_LAYOUT_2D ? (blok + shift)->coefind / K : (blok + shift)->coefind);

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelGEMMCblk2d2d][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelGEMMCblk2d2d][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    /* Get cost in us */
    cost = modelsGetCost3Param( coefs, M, N, K ) * 1e6;

    (void)nimpl;
    return cost;
}

double
blok_gemmsp_cost( struct starpu_task           *task,
                  struct starpu_perfmodel_arch *arch,
                  unsigned                      nimpl )
{
    pastix_coefside_t sideA;
    pastix_coefside_t sideB;
    pastix_trans_t    trans;
    const SolverCblk *cblk;
    SolverCblk       *fcblk;
    pastix_int_t      blok_mk, blok_nk, blok_mn;
    sopalin_data_t   *sopalin_data;
    pastix_int_t      M, N, K;
    SolverBlok *blokA, *blokB, *lblk;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sideA, &sideB, &trans, &cblk, &fcblk,
                                &blok_mk, &blok_nk, &blok_mn, &sopalin_data );

    M = 0;
    N = 0;
    K = cblk_colnbr( cblk );

    blokA = cblk->fblokptr + blok_mk;
    blokB = cblk->fblokptr + blok_nk;
    lblk  = cblk[1].fblokptr;

    M = blok_rownbr(blokA);
    while( (blokA < lblk) &&
           (blokA[0].fcblknm == blokA[1].fcblknm) &&
           (blokA[0].lcblknm == blokA[1].lcblknm) )
    {
        blokA++;
        M += blok_rownbr(blokA);
    }

    N = blok_rownbr(blokB);
    while( (blokB < lblk) &&
           (blokB[0].fcblknm == blokB[1].fcblknm) &&
           (blokB[0].lcblknm == blokB[1].lcblknm) )
    {
        blokB++;
        N += blok_rownbr(blokB);
    }

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelGEMMBlok2d2d][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelGEMMBlok2d2d][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    /* Get cost in us */
    cost = modelsGetCost3Param( coefs, M, N, K ) * 1e6;

    (void)nimpl;
    return cost;
}

double
blok_trsmsp_cost( struct starpu_task           *task,
                  struct starpu_perfmodel_arch *arch,
                  unsigned                      nimpl)
{
    pastix_coefside_t coef;
    pastix_side_t     side;
    pastix_uplo_t     uplo;
    pastix_trans_t    trans;
    pastix_diag_t     diag;
    SolverCblk       *cblk;
    pastix_int_t      blok_m;
    sopalin_data_t   *sopalin_data;
    SolverBlok       *blok;
    SolverBlok       *lblk;
    pastix_int_t      M, N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &coef, &side, &uplo, &trans, &diag,
                                &cblk, &blok_m, &sopalin_data );

    N = cblk_colnbr( cblk );
    blok = cblk->fblokptr + blok_m;
    lblk = cblk[1].fblokptr;

    M = blok_rownbr(blok);
    while( (blok < lblk) &&
           (blok[0].fcblknm == blok[1].fcblknm) &&
           (blok[0].lcblknm == blok[1].lcblknm) )
    {
        blok++;
        M += blok_rownbr(blok);
    }

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelTRSMBlok2d][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelTRSMBlok2d][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    /* Get cost in us */
    cost = modelsGetCost2Param( coefs, M, N ) * 1e6;

    (void)nimpl;
    return cost;
}

double
cblk_getrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    M, N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );
    M = cblk->stride - N;

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelGETRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += 2. * modelsGetCost2Param( coefs, M, N );
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelGETRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += 2. * modelsGetCost2Param( coefs, M, N );
    }
    else {
        coefs = NULL;
        assert(0);
    }

    (void)nimpl;

    /* return cost in us */
    return cost * 1.e6;
}

double
blok_getrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelGETRF][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelGETRF][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    (void)nimpl;

    /* Compute cost in us */
    cost = modelsGetCost1Param( coefs, N ) * 1e6;

    return cost;
}

double
cblk_potrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    M, N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );
    M = cblk->stride - N;

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelPOTRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelPOTRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else {
        assert(0);
    }

    (void)nimpl;

    /* return cost in us */
    return cost * 1.e6;
}

double
blok_potrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelPOTRF][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelPOTRF][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    (void)nimpl;

    /* Compute cost in us */
    cost = modelsGetCost1Param( coefs, N ) * 1e6;

    return cost;
}

double
cblk_pxtrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    M, N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );
    M = cblk->stride - N;

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelPXTRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelPXTRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else {
        assert(0);
    }

    (void)nimpl;

    /* return cost in us */
    return cost * 1.e6;
}

double
blok_pxtrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelPXTRF][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelPXTRF][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    (void)nimpl;

    /* Compute cost in us */
    cost = modelsGetCost1Param( coefs, N ) * 1e6;

    return cost;
}

double
cblk_hetrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    M, N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );
    M = cblk->stride - N;

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelHETRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelHETRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else {
        assert(0);
    }

    (void)nimpl;

    /* return cost in us */
    return cost * 1.e6;
}

double
blok_hetrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelHETRF][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelHETRF][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    (void)nimpl;

    /* Compute cost in us */
    cost = modelsGetCost1Param( coefs, N ) * 1e6;

    return cost;
}

double
cblk_sytrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    M, N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );
    M = cblk->stride - N;

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelSYTRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelSYTRF][0]);
        cost = modelsGetCost1Param( coefs, N );

        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelTRSMCblk2d][0]);
        cost += modelsGetCost2Param( coefs, M, N );
    }
    else {
        assert(0);
    }

    (void)nimpl;

    /* return cost in us */
    return cost * 1.e6;
}

double
blok_sytrf_cost( struct starpu_task           *task,
                 struct starpu_perfmodel_arch *arch,
                 unsigned                      nimpl )
{
    sopalin_data_t *sopalin_data;
    SolverCblk     *cblk;
    pastix_int_t    N;
    double *coefs, cost = 0.;

    starpu_codelet_unpack_args( task->cl_arg, &sopalin_data, &cblk );

    N = cblk_colnbr( cblk );

    if (arch->devices->type == STARPU_CPU_WORKER) {
        coefs = &((*(sopalin_data->cpu_coefs))[PastixKernelSYTRF][0]);
    }
    else if(arch->devices->type == STARPU_CUDA_WORKER) {
        coefs = &((*(sopalin_data->gpu_coefs))[PastixKernelSYTRF][0]);
    }
    else {
        coefs = NULL;
        assert(0);
    }

    (void)nimpl;

    /* Compute cost in us */
    cost = modelsGetCost1Param( coefs, N ) * 1e6;

    return cost;
}
