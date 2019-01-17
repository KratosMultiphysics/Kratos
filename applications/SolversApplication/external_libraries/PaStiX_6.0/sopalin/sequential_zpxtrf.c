/**
 *
 * @file sequential_zpxtrf.c
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
 * @precisions normal z -> z c
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
sequential_zpxtrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    pastix_complex64_t *work;
    pastix_int_t  i, lwork;
    (void)sopalin_data;

    lwork = datacode->gemmmax;
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
        lwork = pastix_imax( lwork, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work, lwork, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        /* Compute */
        cpucblk_zpxtrfsp1d( datacode, cblk,
                            work, lwork );
    }

    memFree_null( work );
}

void
thread_pzpxtrf( isched_thread_t *ctx, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work;
    pastix_int_t  i, ii, lwork;
    pastix_int_t  tasknbr, *tasktab;
    int rank = ctx->rank;

    lwork = datacode->gemmmax;
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
        lwork = pastix_imax( lwork, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work, lwork, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            continue;

        /* Wait */
        do { } while( cblk->ctrbcnt );

        /* Compute */
        cpucblk_zpxtrfsp1d( datacode, cblk,
                            work, lwork );
    }

    memFree_null( work );
}

void
thread_zpxtrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzpxtrf, sopalin_data );
}

static void (*zpxtrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zpxtrf,
    thread_zpxtrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zpxtrf,
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_zpxtrf
#else
    NULL
#endif
};

void
sopalin_zpxtrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zpxtrf)(pastix_data_t *, sopalin_data_t *) = zpxtrf_table[ sched ];

    if (zpxtrf == NULL) {
        zpxtrf = thread_zpxtrf;
    }
    zpxtrf( pastix_data, sopalin_data );

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( pastix_data, sopalin_data->solvmtx, "pxtrf.txt" );
#endif
}
