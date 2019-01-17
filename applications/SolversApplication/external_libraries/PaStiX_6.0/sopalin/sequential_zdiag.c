/**
 *
 * @file sequential_zdiag.c
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
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

void
sequential_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data,
                  int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t  i, cblknbr;
    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

    cblk = datacode->cblktab;
    cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
    for (i=0; i<cblknbr; i++, cblk++) {
        solve_cblk_zdiag( cblk, nrhs,
                          b + cblk->lcolidx, ldb, NULL );
    }
}

struct args_zdiag_t
{
    pastix_data_t  *pastix_data;
    sopalin_data_t *sopalin_data;
    int nrhs;
    pastix_complex64_t *b;
    int ldb;
};

void
thread_pzdiag( isched_thread_t *ctx, void *args )
{
    struct args_zdiag_t *arg = (struct args_zdiag_t*)args;
    pastix_data_t      *pastix_data  = arg->pastix_data;
    sopalin_data_t     *sopalin_data = arg->sopalin_data;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    pastix_complex64_t *b = arg->b;
    int nrhs  = arg->nrhs;
    int ldb   = arg->ldb;
    SolverCblk *cblk;
    Task       *t;
    pastix_int_t i, ii, cblknbr;
    pastix_int_t tasknbr, *tasktab;
    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
    int rank = ctx->rank;

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];
    cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;

        if ( t->cblknum >= cblknbr ) {
            continue;
        }
        cblk = datacode->cblktab + t->cblknum;
        solve_cblk_zdiag( cblk, nrhs,
                          b + cblk->lcolidx, ldb, NULL );
    }
}

void
thread_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data,
              int nrhs, pastix_complex64_t *b, int ldb )
{
    struct args_zdiag_t args_zdiag = {pastix_data, sopalin_data, nrhs, b, ldb};
    isched_parallel_call( pastix_data->isched, thread_pzdiag, &args_zdiag );
}

static void (*zdiag_table[4])(pastix_data_t *, sopalin_data_t *,
                              int, pastix_complex64_t *, int) = {
    sequential_zdiag,
    thread_zdiag,
#if defined(PASTIX_WITH_PARSEC)
    NULL, /* parsec_zdiag not yet implemented */
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_zdiag
#else
    NULL
#endif
};

void
sopalin_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data,
               int nrhs, pastix_complex64_t *b, int ldb )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zdiag)(pastix_data_t *, sopalin_data_t *, int, pastix_complex64_t *, int) = zdiag_table[ sched ];

    if (zdiag == NULL) {
        zdiag = thread_zdiag;
    }
    zdiag( pastix_data, sopalin_data, nrhs, b, ldb );
}
