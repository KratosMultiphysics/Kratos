/**
 *
 * @file solver_backup.c
 *
 * PaStiX solver structure routines.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#include "solver.h"

/**
 * @ingroup blend_dev_solver
 * @brief Structure to store backup of counter modified during numerical factorization and solve steps.
 */
struct SolverBackup_s {
    pastix_int_t  arftmax;        /**< Maximum area of FanIn Target: doubled by LU decomposition                   */
    pastix_int_t  nbftmax;        /**< Maximum number of FanIn Target: doubled by LU decomposition                 */
    pastix_int_t *task_ctrbcnt;   /**< Number of contribution: the counter is decreased to 0 during factorization  */
    pastix_int_t *task_ftgtcnt;   /**< Number of FanIn contrib: the counter is decreased to 0 during factorization */
    pastix_int_t *fanin_ctrbnbr;  /**< Number of contribution to FanIn: decreased during facto and solve           */
    pastix_int_t *fanin_prionum;  /**< Replaced by the number of msg packed during factorization sends             */
    pastix_int_t *symbol_cblknum; /**< Replaced by the negative FanIn index during facto and solve                 */
    pastix_int_t  symbol_nodenbr; /**< ???                                                                         */
};

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Initialize the backup structure.
 *
 * This function saves the initial values of the counters before
 * modifications. See structure description for more information about what is
 * saved and why.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure holding information for factorization
 *          and solve steps.
 *
 *******************************************************************************
 *
 * @retval Pointer to the allocated backup structure with the copy of
 *         information that might be destroyed.
 *
 *******************************************************************************/
SolverBackup_t *
solverBackupInit( const SolverMatrix *solvmtx )
{
    SolverBackup_t *b;
    pastix_int_t i;

    MALLOC_INTERN( b, 1, SolverBackup_t );
    memset( b, 0, sizeof(SolverBackup_t) );

    b->arftmax = solvmtx->arftmax;
    b->nbftmax = solvmtx->nbftmax;

    if (solvmtx->tasknbr)
    {
        Task *task = solvmtx->tasktab;

        MALLOC_INTERN(b->task_ctrbcnt, solvmtx->tasknbr, pastix_int_t);
        MALLOC_INTERN(b->task_ftgtcnt, solvmtx->tasknbr, pastix_int_t);

        for (i=0; i<solvmtx->tasknbr; i++, task++)
        {
            b->task_ctrbcnt[i] = task->ctrbcnt;
            b->task_ftgtcnt[i] = task->ftgtcnt;
        }
    }

    if (solvmtx->ftgtnbr)
    {
        solver_ftgt_t *ftgt = solvmtx->ftgttab;

        MALLOC_INTERN(b->fanin_ctrbnbr, solvmtx->ftgtnbr, pastix_int_t);
        MALLOC_INTERN(b->fanin_prionum, solvmtx->ftgtnbr, pastix_int_t);

        for (i=0; i<solvmtx->ftgtnbr; i++, ftgt++)
        {
            b->fanin_ctrbnbr[i] = ftgt->infotab[FTGT_CTRBNBR];
            b->fanin_prionum[i] = ftgt->infotab[FTGT_PRIONUM];
        }
    }

    if (solvmtx->bloknbr) {
        SolverBlok *blok = solvmtx->bloktab;

        MALLOC_INTERN(b->symbol_cblknum, solvmtx->bloknbr, pastix_int_t);

        for (i=0; i<solvmtx->bloknbr; i++, blok++) {
            b->symbol_cblknum[i] = blok->fcblknm;
        }
    }

    b->symbol_nodenbr = solvmtx->nodenbr;

    {
        SolverCblk *cblk = solvmtx->cblktab;
        for (i=0; i<solvmtx->cblknbr; i++, cblk++)
        {
            cblk->ctrbcnt = cblk[1].brownum - cblk[0].brownum;
        }
    }
    return b;
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Restore initial values.
 *
 * Restore counter values to be able to call a second factorization or solve
 * step. The amount of information restored depends on the value of
 * solvmtx->restore. If it is equal to:
 *     - 0: Nothing is restored
 *     - 1: A solve step has been performed and partial information is restored.
 *     - 2: A factorization step has been performed and full information is restored.
 * The value of solvmtx->restore is noramally initialized to 0 during the
 * structure creation, and then set to the correct value by the routine
 * modifying the solvmtx structure.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure holding information for factorization
 *          and solve steps.
 *          On exit, the counters have been restored to their original value
 *          stored in the backup structure.
 *
 * @param[in] b
 *          The backup structure pointer returned by the call to solverBackupInit().
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the data has been restored successfuly.
 * @retval PASTIX_ERR_BADPARAMETER if one of the parameter is incorrect.
 *
 *******************************************************************************/
int
solverBackupRestore( SolverMatrix         *solvmtx,
                     const SolverBackup_t *b       )
{
    pastix_int_t i;

    if ( solvmtx == NULL || b == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( solvmtx->restore == 0) {
        return PASTIX_SUCCESS;
    }

    /* After factorization */
    if ( solvmtx->restore == 2 ) {
        solvmtx->arftmax = b->arftmax;
        solvmtx->nbftmax = b->nbftmax;

        if (solvmtx->tasknbr)
        {
            Task *task = solvmtx->tasktab;

            for (i=0; i<solvmtx->tasknbr; i++, task++)
            {
                task->ctrbcnt = b->task_ctrbcnt[i];
                task->ftgtcnt = b->task_ftgtcnt[i];
            }
        }
    }

    if (solvmtx->ftgtnbr)
    {
        solver_ftgt_t *ftgt = solvmtx->ftgttab;

        for (i=0; i<solvmtx->ftgtnbr; i++, ftgt++)
        {
            ftgt->infotab[FTGT_CTRBNBR] = b->fanin_ctrbnbr[i];
            ftgt->infotab[FTGT_CTRBCNT] = b->fanin_ctrbnbr[i];
            ftgt->infotab[FTGT_PRIONUM] = b->fanin_prionum[i];
        }
    }

    if (solvmtx->bloknbr) {
        SolverBlok *blok = solvmtx->bloktab;

        for (i=0; i<solvmtx->bloknbr; i++, blok++) {
            blok->fcblknm = b->symbol_cblknum[i];
            blok->gpuid = GPUID_UNDEFINED;
        }
    }
    solvmtx->nodenbr = b->symbol_nodenbr;

    if (solvmtx->cblknbr) {
        SolverCblk *cblk = solvmtx->cblktab;

        for (i=0; i<solvmtx->cblknbr; i++, cblk++) {
            cblk->gpuid = GPUID_UNDEFINED;
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Free the solver backup data structure.
 *
 * Clean the data structure holding the information backup and free the given
 * pointer because it has necessarily been allocated in solverBackupInit().
 *
 *******************************************************************************
 *
 * @param[inout] b
 *          The backup structure to destroy. On exit, b cannot be used anymore.
 *
 *******************************************************************************/
void
solverBackupExit( SolverBackup_t *b )
{
    if (b->task_ctrbcnt)
    {
        memFree_null(b->task_ctrbcnt);
        memFree_null(b->task_ftgtcnt);
    }

    if (b->fanin_ctrbnbr)
    {
        memFree_null(b->fanin_ctrbnbr);
        memFree_null(b->fanin_prionum);
    }

    if (b->symbol_cblknum) {
        memFree_null(b->symbol_cblknum);
    }
    memFree_null(b);
}
