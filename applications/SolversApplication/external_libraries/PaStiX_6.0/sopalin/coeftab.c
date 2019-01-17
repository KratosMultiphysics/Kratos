/**
 *
 * @file coeftab.c
 *
 * PaStiX coefficient array initialization and free routines.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "bcsc.h"
#include "isched.h"
#include "solver.h"
#include "coeftab.h"
#include "pastix_zcores.h"
#include "pastix_ccores.h"
#include "pastix_dcores.h"
#include "pastix_scores.h"
#include "pastix_zlrcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "sopalin/parsec/pastix_parsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "sopalin/starpu/pastix_starpu.h"
#endif

coeftab_fct_memory_t coeftabMemory[4] =
{
    coeftab_smemory, coeftab_dmemory, coeftab_cmemory, coeftab_zmemory
};

/**
 * @brief Internal structure specific to the parallel call of pcoeftabInit()
 */
struct coeftabinit_s {
    const SolverMatrix  *datacode; /**< The sovler matrix                         */
    const pastix_bcsc_t *bcsc;     /**< The internal block CSC                    */
    char               **dirtemp;  /**< The pointer to the output directory       */
    pastix_coefside_t    side;     /**< The side of the matrix beeing initialized */
};

/**
 *******************************************************************************
 *
 * @brief Internal routine called by each static thread to Initialize the solver
 * matrix structure.
 *
 * This routine is the routine called by each thread in the static scheduler and
 * launched by the coeftabinit().
 *
 *******************************************************************************
 *
 * @param[inout] ctx
 *          The internal scheduler context
 *
 * @param[in] args
 *          The data structure specific to the function cpucblk_zinit()
 *
 *******************************************************************************/
void
pcoeftabInit( isched_thread_t *ctx,
              void            *args )
{
    struct coeftabinit_s *ciargs   = (struct coeftabinit_s*)args;
    const SolverMatrix   *datacode = ciargs->datacode;
    const pastix_bcsc_t  *bcsc     = ciargs->bcsc;
    char                **dirtemp  = ciargs->dirtemp;
    pastix_coefside_t     side     = ciargs->side;
    pastix_int_t i, itercblk;
    pastix_int_t task;
    int rank = ctx->rank;

    void (*initfunc)( pastix_coefside_t, const SolverMatrix*,
                      const pastix_bcsc_t*, pastix_int_t, char **) = NULL;

    switch( bcsc->flttype ) {
    case PastixComplex32:
        initfunc = cpucblk_cinit;
        break;
    case PastixComplex64:
        initfunc = cpucblk_zinit;
        break;
    case PastixFloat:
        initfunc = cpucblk_sinit;
        break;
    case PastixDouble:
    case PastixPattern:
    default:
        initfunc = cpucblk_dinit;
    }

    for (i=0; i < datacode->ttsknbr[rank]; i++)
    {
        task = datacode->ttsktab[rank][i];
        itercblk = datacode->tasktab[task].cblknum;

        /* Init as full rank */
        initfunc( side, datacode, bcsc, itercblk, dirtemp );
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the solver matrix structure
 *
 * This routine is a parallel routine to initialize the solver matrix structure
 * through the internal static scheduler
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that hold the solver matrix to initialize.
 *
 * @param[in] side
 *          Describe the side(s) of the matrix that must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 *******************************************************************************/
void
coeftabInit( pastix_data_t    *pastix_data,
             pastix_coefside_t side )
{
    struct coeftabinit_s args;

    args.datacode   = pastix_data->solvmatr;
    args.bcsc       = pastix_data->bcsc;
    args.side       = side;
    args.dirtemp    = &(pastix_data->dirtemp);

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    /* Make sure dirtemp is initialized before calling it with multiple threads */
    pastix_gendirtemp( args.dirtemp );
#endif

    isched_parallel_call( pastix_data->isched, pcoeftabInit, &args );
}

/**
 *******************************************************************************
 *
 * @brief Free the solver matrix structure
 *
 * This routine free all data structure refereing to the solver matrix L, even
 * the runtime descriptors if present.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure of the problem.
 *
 *******************************************************************************/
void
coeftabExit( SolverMatrix *solvmtx )
{
    pastix_int_t i;

#if defined(PASTIX_WITH_PARSEC)
    {
        if ( solvmtx->parsec_desc != NULL ) {
            parsec_sparse_matrix_destroy( solvmtx->parsec_desc );
            free( solvmtx->parsec_desc );
        }
        solvmtx->parsec_desc = NULL;
    }
#endif
#if defined(PASTIX_WITH_STARPU)
    {
        if ( solvmtx->starpu_desc != NULL ) {
            starpu_sparse_matrix_destroy( solvmtx->starpu_desc );
            free( solvmtx->starpu_desc );
        }
        solvmtx->starpu_desc = NULL;
    }
#endif

    /* Free arrays of solvmtx */
    if(solvmtx->cblktab)
    {
        for (i = 0; i < solvmtx->cblknbr; i++)
        {
            if (solvmtx->cblktab[i].lcoeftab) {
                memFree_null(solvmtx->cblktab[i].lcoeftab);
            }

            if (solvmtx->cblktab[i].ucoeftab) {
                memFree_null(solvmtx->cblktab[i].ucoeftab);
            }

            if (solvmtx->cblktab[i].cblktype & CBLK_COMPRESSED) {
                SolverBlok *blok  = solvmtx->cblktab[i].fblokptr;
                SolverBlok *lblok = solvmtx->cblktab[i+1].fblokptr;

                if ( blok->LRblock != NULL ) {
                    for (; blok<lblok; blok++) {
                        core_zlrfree(blok->LRblock);
                        if (solvmtx->factotype == PastixFactLU)
                            core_zlrfree(blok->LRblock+1);
                    }
                    free(solvmtx->cblktab[i].fblokptr->LRblock);
                }
            }
        }
    }
}

/**
 * @brief Internal structure specific to the parallel call of pcoeftabComp()
 */
struct coeftabcomp_s {
    SolverMatrix        *solvmtx; /**< The sovler matrix               */
    pastix_coeftype_t    flttype; /**< The arithmetic type             */
    pastix_atomic_lock_t lock;    /**< Lock to protect the gain update */
    pastix_int_t         gain;    /**< The memory gain on output       */
};

/**
 *******************************************************************************
 *
 * @brief Internal routine called by each static thread to Initialize the solver
 * matrix structure.
 *
 * This routine is the routine called by each thread in the static scheduler and
 * launched by the coeftabCompress().
 *
 *******************************************************************************
 *
 * @param[inout] ctx
 *          The internal scheduler context
 *
 * @param[in] args
 *          The data structure specific to the function cpucblk_zcompress()
 *
 *******************************************************************************/
static void
pcoeftabComp( isched_thread_t *ctx,
              void            *args )
{
    struct coeftabcomp_s *ccargs   = (struct coeftabcomp_s*)args;
    SolverMatrix         *solvmtx  = ccargs->solvmtx;
    pastix_coeftype_t     flttype  = ccargs->flttype;
    pastix_atomic_lock_t *lock     = &(ccargs->lock);
    pastix_int_t         *fullgain = &(ccargs->gain);
    SolverCblk           *cblk;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    pastix_int_t i, itercblk;
    pastix_int_t task, gain = 0;
    int rank = ctx->rank;

    pastix_int_t (*compfunc)( pastix_coefside_t, SolverCblk*, pastix_lr_t ) = NULL;

    switch( flttype ) {
    case PastixComplex32:
        compfunc = cpucblk_ccompress;
        break;
    case PastixComplex64:
        compfunc = cpucblk_zcompress;
        break;
    case PastixFloat:
        compfunc = cpucblk_scompress;
        break;
    case PastixDouble:
    case PastixPattern:
    default:
        compfunc = cpucblk_dcompress;
    }

    for (i=0; i < solvmtx->ttsknbr[rank]; i++)
    {
        task     = solvmtx->ttsktab[rank][i];
        itercblk = solvmtx->tasktab[task].cblknum;
        cblk     = solvmtx->cblktab + itercblk;

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            gain += compfunc( side, cblk, solvmtx->lowrank );
        }
    }

    pastix_atomic_lock( lock );
    *fullgain += gain;
    pastix_atomic_unlock( lock );
}

/**
 *******************************************************************************
 *
 * @brief Compress the factorized matrix structure if not already done.
 *
 * This routine compress all column blocks that are marked for compression, and
 * return the amount of memory saved by the compression.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that holds the problem
 *
 *******************************************************************************
 *
 * @return The memory gain resulting from the compression to low-rank format in
 *         number of elements.
 *
 *******************************************************************************/
pastix_int_t
coeftabCompress( pastix_data_t *pastix_data )
{
    struct coeftabcomp_s args;
    pastix_lr_t  *lr;

    args.solvmtx   = pastix_data->solvmatr;
    args.flttype   = pastix_data->bcsc->flttype;
    args.lock      = PASTIX_ATOMIC_UNLOCKED;
    args.gain      = 0;

    /* Set the lowrank properties */
    lr = &(pastix_data->solvmatr->lowrank);
    lr->compress_method     = pastix_data->iparm[IPARM_COMPRESS_METHOD];
    lr->compress_min_width  = pastix_data->iparm[IPARM_COMPRESS_MIN_WIDTH];
    lr->compress_min_height = pastix_data->iparm[IPARM_COMPRESS_MIN_HEIGHT];
    lr->tolerance           = sqrt( pastix_data->dparm[DPARM_COMPRESS_TOLERANCE] );

    isched_parallel_call( pastix_data->isched, pcoeftabComp, (void*)(&args) );

    return args.gain;
}
