/**
 *
 * @file sequential_zhetrf.c
 *
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @precisions normal z -> c
 *
 **/
#include "common.h"
#include "isched.h"
#include "solver.h"
#include "sopalin_data.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "parsec/pastix_zparsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

void
sequential_zhetrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  N, i, lwork1, lwork2;
    (void)sopalin_data;

    lwork1 = datacode->offdmax;
    lwork2 = pastix_imax( datacode->gemmmax, datacode->blokmax );
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
        lwork2 = pastix_imax( lwork2, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work1, lwork1, pastix_complex64_t );
    MALLOC_INTERN( work2, lwork2, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        N = cblk_colnbr( cblk );

        /* Compute */
        cpucblk_zhetrfsp1d( datacode, cblk,
                            /*
                             * Workspace size has been computed without the
                             * diagonal block, thus in order to work with generic
                             * TRSM and GEMM kernels, we must shift the DLh workspace
                             * by the diagonal block size
                             */
                            work1 - (N*N), work2, lwork2 );
    }

    memFree_null( work1 );
    memFree_null( work2 );
}

void
thread_pzhetrf( isched_thread_t *ctx, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  N, i, ii, lwork1, lwork2;
    pastix_int_t  tasknbr, *tasktab;
    int rank = ctx->rank;

    lwork1 = datacode->offdmax;
    lwork2 = pastix_imax( datacode->gemmmax, datacode->blokmax );
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
        lwork2 = pastix_imax( lwork2, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work1, lwork1, pastix_complex64_t );
    MALLOC_INTERN( work2, lwork2, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            continue;

        N = cblk_colnbr( cblk );

        /* Wait */
        do { } while( cblk->ctrbcnt );

        /* Compute */
        cpucblk_zhetrfsp1d( datacode, cblk,
                            /*
                             * Workspace size has been computed without the
                             * diagonal block, thus in order to work with generic
                             * TRSM and GEMM kernels, we must shift the DLh workspace
                             * by the diagonal block size
                             */
                            work1 - (N*N), work2, lwork2 );
    }

    memFree_null( work1 );
    memFree_null( work2 );
}

void
thread_zhetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzhetrf, sopalin_data );
}

static void (*zhetrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zhetrf,
    thread_zhetrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zhetrf,
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_zhetrf
#else
    NULL
#endif
};

void
sopalin_zhetrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zhetrf)(pastix_data_t *, sopalin_data_t *) = zhetrf_table[ sched ];

    if (zhetrf == NULL) {
        zhetrf = thread_zhetrf;
    }
    zhetrf( pastix_data, sopalin_data );

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( pastix_data, sopalin_data->solvmtx, "hetrf.txt" );
#endif
}
