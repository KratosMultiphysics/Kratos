/**
 *
 * @file simu_run.c
 *
 * PaStiX simulation functions.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "common.h"
#include "symbol.h"
#include "extendVector.h"
#include "queue.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "blendctrl.h"
#include "solver.h"
#include "simu.h"
#include "perf.h"

#if defined(PASTIX_BLEND_GENTRACE)
#include <GTG.h>
#endif

/**
 * @addtogroup blend_dev_simu
 * @{
 *
 */

/**
 *******************************************************************************
 *
 * @brief Compute the cost of a communication and its update.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blend control structure that describes the architecture and the
 *          cost of the communication between nodes.
 *
 * @param[in] ftgt
 *          The fan-in for which the cost is computed.
 *
 * @param[in] clustsrc
 *          The index of the source pastix process.
 *
 * @param[in] sync_comm_nbr
 *          The number of simultaneous communication.
 *
 * @param[out] send
 *          The time cost of the send operation.
 *
 * @param[out] add
 *          The time cost of the addition operation.
 *
 *******************************************************************************/
static inline void
simu_computeFtgtCosts( const BlendCtrl     *ctrl,
                       const solver_ftgt_t *ftgt,
                       pastix_int_t         clustsrc,
                       pastix_int_t         sync_comm_nbr,
                       double              *send,
                       double              *add )
{
    pastix_int_t M, N;
    pastix_int_t clustdst = ctrl->core2clust[ftgt->infotab[FTGT_PROCDST]];
    double startup, bandwidth, addcost;

    *send = 0.;
    *add  = 0.;

    if( clustsrc == clustdst ) {
        return;
    }

    assert( (clustsrc >= 0) && (clustdst >= 0) );

    N = (ftgt->infotab[FTGT_LCOLNUM] - ftgt->infotab[FTGT_FCOLNUM] + 1);
    M = (ftgt->infotab[FTGT_LROWNUM] - ftgt->infotab[FTGT_FROWNUM] + 1);

    assert( (N > 0) && (M > 0) );

    getCommunicationCosts( ctrl, clustsrc, clustdst, sync_comm_nbr, &startup, &bandwidth );

    *send = (startup + bandwidth * (M * N * sizeof(double) + FTGT_MAXINFO * sizeof(pastix_int_t)));
    addcost = PERF_GEAM( M, N );
    *add = addcost > 0. ? addcost : 0.0;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Compute the number of contributions to each block.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The symbol matrix structure describing the problem.
 *
 * @param[inout] simuctrl
 *          The main simulation structure. On exit, the ctrbcnt field of the
 *          blocks is updated with the number of contribution that they each
 *          should receive.
 *
 * @param[in] ricar
 *          True if ILU(k) factorization is applied to change the algorithm to
 *          compute the number of contribution.
 *
 *******************************************************************************/
static inline void
simu_computeBlockCtrbNbr(const symbol_matrix_t *symbptr,
                               SimuCtrl        *simuctrl,
                               pastix_int_t     ricar )
{
    pastix_int_t i, j, k;
    pastix_int_t facebloknum, firstbloknum;

    /*
     * Compute the number of contribution per block to each block.
     * Might be optimized if we computed the input graph before.
     */
    {
        symbol_cblk_t *curcblk;

        curcblk = symbptr->cblktab;
        for(i=0; i<symbptr->cblknbr; i++, curcblk++)
        {
            pastix_int_t fbloknum = curcblk[0].bloknum + 1;
            pastix_int_t lbloknum = curcblk[1].bloknum;

            /* 1D cblk computed */
            for(j=fbloknum; j<lbloknum; j++)
            {
                firstbloknum = 0;

                /* Add contribution due to E2 */
                for(k=j; k<lbloknum; k++)
                {
                    facebloknum = pastixSymbolGetFacingBloknum( symbptr, j, k, firstbloknum, ricar );
                    if(facebloknum >= 0) {
                        simuctrl->bloktab[facebloknum].ctrbcnt++;
                        firstbloknum = facebloknum;
                    }
                }
            }
        }
    }

    /* Set up the task ctrbcnt and cblkcnt */
    {
        SimuTask *task = simuctrl->tasktab;

        for(i=0;i<simuctrl->tasknbr;i++)
        {
            pastix_int_t fbloknum = symbptr->cblktab[task->cblknum  ].bloknum;
            pastix_int_t lbloknum = symbptr->cblktab[task->cblknum+1].bloknum;

            task->ctrbcnt = 0;
            for(j=fbloknum; j<lbloknum; j++) {
                task->ctrbcnt += simuctrl->bloktab[j].ctrbcnt;
            }

            simuctrl->cblktab[task->cblknum].ctrbcnt = task->ctrbcnt;
            task++;
        }
    }
}


/**
 *******************************************************************************
 *
 * @brief Print the number of contribution per cblk and block for debug.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blendctrl structure with the simulation parameters.
 *
 * @param[in] symbptr
 *          The symbol matrix structure describing the problem.
 *
 * @param[in] simuctrl
 *          The main simulation structure.
 *
 *******************************************************************************/
static inline void
simu_printBlockCtrbNbr( const BlendCtrl       *ctrl,
                        const symbol_matrix_t *symbptr,
                        const SimuCtrl        *simuctrl )
{
    FILE *fd1 = NULL;
    FILE *fd2 = NULL;
    pastix_int_t i, j;
    symbol_cblk_t *curcblk;

    fd1 = pastix_fopenw( ctrl->dirtemp, "contribblok.txt", "w" );
    if ( fd1 == NULL ) {
        return;
    }
    fd2 = pastix_fopenw( ctrl->dirtemp, "contribcblk.txt", "w" );
    if ( fd2 == NULL ) {
        return;
    }

    curcblk = symbptr->cblktab;
    for(i=0; i<symbptr->cblknbr; i++, curcblk++)
    {
        pastix_int_t fbloknum = curcblk[0].bloknum + 1;
        pastix_int_t lbloknum = curcblk[1].bloknum;

        /* 1D cblk computed */
        for(j=fbloknum; j<lbloknum; j++)
        {
            fprintf(fd1, "%ld %ld\n", (long)j, (long)simuctrl->bloktab[j].ctrbcnt);
        }
#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
        fprintf(fd2, "%ld %ld %ld\n", (long)i, (long)simuctrl->cblktab[i].ctrbcnt, (long)curcblk->split_cblk);
#else
        fprintf(fd2, "%ld %ld\n", (long)i, (long)simuctrl->cblktab[i].ctrbcnt);
#endif
    }

    fclose( fd1 );
    fclose( fd2 );
}

/**
 *******************************************************************************
 *
 * @brief Insert a task in the ready task queues of all its candidates.
 *
 * This function according to the ready date of a task put this task on the
 * ready queue of a processor.
 * When the ready date of a task is inferior to the proc timer then he
 * task is ordered according to its priorities in the elimination tree.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blend control structure to provide the candtab and the
 *          core2clust arrays.
 *
 * @param[inout] simuctrl
 *          The main simulation structure. On exit, the ready tasks queues of
 *          the candidates for tasknum are updated.
 *
 * @param[in] tasknum
 *          The index of the task to insert as a ready task.
 *
 *******************************************************************************/
static inline void
simu_putInAllReadyQueues( const BlendCtrl *ctrl,
                          SimuCtrl        *simuctrl,
                          pastix_int_t     tasknum )
{
    const SimuTask *task     = simuctrl->tasktab + tasknum;
    const Cand     *cblkcand = ctrl->candtab + task->cblknum;
    SimuProc  *sproc;
    double ready_date = 0.0;
    pastix_int_t procnum;
    pastix_int_t bloknum = task->bloknum;
    pastix_int_t treelevel = cblkcand->treelevel;

    assert( tasknum != -1 );

    /* Get the ready date of the task on the processor passed in parameter */
    if( cblkcand->fccandnum == cblkcand->lccandnum )
    {
        ready_date = timerVal( &(task->time) );
        sproc = &(simuctrl->proctab[cblkcand->fcandnum]);

        for(procnum =  cblkcand->fcandnum;
            procnum <= cblkcand->lcandnum; procnum++, sproc++)
        {
            if(ready_date > timerVal(TIMER(procnum))) {
                pqueuePush2( sproc->futuretask, tasknum, ready_date, treelevel);
            }
            else {
                pqueuePush2( sproc->readytask, tasknum, treelevel, bloknum);
            }
        }
    }
    else
    {
        sproc = &(simuctrl->proctab[cblkcand->fcandnum]);

        for(procnum =  cblkcand->fcandnum;
            procnum <= cblkcand->lcandnum; procnum++, sproc++)
        {
            ready_date = timerVal( simuctrl->ftgttimetab + CLUST2INDEX(bloknum, ctrl->core2clust[procnum]) );

            if(ready_date > timerVal(TIMER(procnum))) {
                pqueuePush2( sproc->futuretask, tasknum, ready_date, treelevel);
            }
            else {
                pqueuePush2( sproc->readytask, tasknum, treelevel, bloknum);
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Look for the best next couple (tasknum, corenum) that is ready to be
 * executed.
 *
 * This function is the main and more costly one. It looks for each worker,
 * which tasks is the first one available for execution, and from all those
 * couples, which one is the first one to finish.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blend control structure to provide the candtab and the
 *          core2clust arrays.
 *
 * @param[inout] simuctrl
 *          The main simulation structure. On exit, the structure is updated
 *          with the extraction of the next best task to run.
 *
 * @param[out] procnumptr
 *          The index of the candidate to run the task.
 *
 *******************************************************************************
 *
 * @return The next task selected for execution in the simulator. The worker
 *         selected is treturned in the procnumptr field.
 *
 *******************************************************************************/
static inline pastix_int_t
simu_getNextTaskNextProc( const BlendCtrl *ctrl,
                          SimuCtrl        *simuctrl,
                          pastix_int_t    *procnumptr )
{
    pastix_int_t p;
    pastix_int_t procnum = -1;
    pastix_int_t tasknum;
    double earlytimeready = PASTIX_INT_MAX;
    double earlyproctimer = PASTIX_INT_MAX;
    double timeready;
    pastix_int_t earlytask = -1;

    /* Find the earlier task in the processor heaps */
    for(p=0;p<ctrl->total_nbcores;p++)
    {
        tasknum = -1;
        /*
         * First we search the earlier task in the set of task whose ready date
         * is < proc timer
         */
        while(pqueueSize(simuctrl->proctab[p].readytask)>0)
        {
            tasknum = pqueueRead(simuctrl->proctab[p].readytask);
            if( simuctrl->bloktab[simuctrl->tasktab[tasknum].bloknum].ownerclust >= 0 )
            {
                /* This task have to be remove from the heap (already mapped) */
                pqueuePop(simuctrl->proctab[p].readytask);
                tasknum = -1;
            }
            else
                break;
        }
        /*
         * We found no task which ready date is < proc timer so we search one
         * that minimizes ready date - proc-timer
         */
        if(tasknum == -1)
        {
            while(pqueueSize(simuctrl->proctab[p].futuretask)>0)
            {
                tasknum = pqueueRead(simuctrl->proctab[p].futuretask);
                if( simuctrl->bloktab[simuctrl->tasktab[tasknum].bloknum].ownerclust >= 0 )
                {
                    /* This task have to be remove from the heap (already mapped) */
                    pqueuePop(simuctrl->proctab[p].futuretask);
                    tasknum = -1;
                }
                else {
                    break;
                }
            }
        }

        if(tasknum != -1)
        {
            timeready = MAX(timerVal(TIMER(p)),
                            timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(simuctrl->tasktab[tasknum].bloknum, ctrl->core2clust[p])])));

            /*
             * We prevent to distribute on the same processor set when all time
             * are equals
             */
            if((timeready == earlytimeready) && (timerVal(TIMER(p)) < earlyproctimer))
            {
                procnum = p;
                earlyproctimer = timerVal(TIMER(p));
                earlytask = tasknum;
                earlytimeready = timeready;
            }

            if(timeready < earlytimeready)
            {
                procnum  = p;
                earlytask = tasknum;
                earlytimeready = timeready;
            }
        }
    }
    if(procnum != -1)
    {
        if(pqueueSize(simuctrl->proctab[procnum].readytask)>0)
        {
            assert(earlytask == pqueuePop(simuctrl->proctab[procnum].readytask));
        }
        else
        {
            assert(earlytask == pqueuePop(simuctrl->proctab[procnum].futuretask));
        }
    }
    *procnumptr = procnum;
    return earlytask;
}

/**
 *******************************************************************************
 *
 * @brief Compute the instant t where the task will be received by a node.
 *
 * Compute the time the cblk would have RECEIVED and ADDED all its contributions
 * if it was mapped on a given cand CLUSTER.
 * @warning These times do not include add time for fan in target
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The blend control structure to provide the candtab and the
 *          core2clust arrays.
 *
 * @param[in] symbptr
 *          The symbol matrix structure describing the problem.
 *
 * @param[inout] simuctrl
 *          The main simulation structure. On exit, the ready tasks queues of
 *          the candidates for tasknum are updated.
 *
 * @param[in] tasknum
 *          The index of the task to insert as a ready task.
 *
 *******************************************************************************/
static inline void
simu_computeTaskReceiveTime( const BlendCtrl       *ctrl,
                             const symbol_matrix_t *symbptr,
                                   SimuCtrl        *simuctrl,
                                   pastix_int_t     tasknum )
{
    pastix_int_t i, j;
    double lftgttime = 0;
    double sftgttime = 0;
    pastix_int_t   lftgtnum  = -1;
    pastix_int_t   cblknum;
    pastix_int_t   bloknum;
    pastix_int_t   clustdst;

    bloknum = simuctrl->tasktab[tasknum].bloknum;
    cblknum = simuctrl->tasktab[tasknum].cblknum;

    /* If the task is local, all sons sending contributions are local => no treatment */
    if( ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum ) {
        return;
    }

    /*
     * Compute the cblk on proc timer that is time the cblk would have received
     * all its contributions if it was mapped on a given cand processor These
     * times INCLUDE add time for fan in target !!
     */

    /* Compute receive time (time at which a non-local processor should received the target) */
    /* find the latest ftgt receive time and the second latest*/
    for(i=simuctrl->bloktab[bloknum].ftgtnum; i<simuctrl->bloktab[bloknum+1].ftgtnum; i++)
    {
        /* Source of this ftgt */
        clustdst = INDEX2CLUST(i, bloknum);

        /* Task with several cand proc */
        /* The information about ftgt costs are in the ftgt of the diagonal block;
         this loop sums the cost of all the ftgt received by the blocks in this column block */
        if(simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR]>0) {
            for(j=bloknum;j<symbptr->cblktab[cblknum+1].bloknum;j++)
            {
                if(simuctrl->ftgttab[simuctrl->bloktab[j].ftgtnum + i-simuctrl->bloktab[bloknum].ftgtnum].ftgt.infotab[FTGT_CTRBNBR]>0)
                {
                    double send, add;

                    simu_computeFtgtCosts( ctrl, &(simuctrl->ftgttab[CLUST2INDEX(j, clustdst)].ftgt), clustdst,
                                           ctrl->candtab[cblknum].lccandnum - ctrl->candtab[cblknum].fccandnum + 1,
                                           &send, &add );

                    simuctrl->ftgttab[i].costadd  += add;
                    simuctrl->ftgttab[i].costsend += send;
                }
            }
        }

#if defined(PASTIX_DEBUG_BLEND)
        if(!(simuctrl->ftgttab[i].costsend >= 0.0)) {
            errorPrint("ftgt %ld costsend %f", (long)i, simuctrl->ftgttab[i].costsend);
        }
        if(!(simuctrl->ftgttab[i].costadd >= 0.0)) {
            errorPrint("ftgt %ld costadd %f", (long)i, simuctrl->ftgttab[i].costadd);
        }

        assert(simuctrl->ftgttab[i].costsend >= 0.0);
        assert(simuctrl->ftgttab[i].costadd >= 0.0);
#endif

        /* ftgttab[].timerecv is the time this ftgt will be receive */
        timerSet(&(simuctrl->ftgttab[i].timerecv), timerVal(&(simuctrl->ftgttimetab[i])) + simuctrl->ftgttab[i].costsend + simuctrl->ftgttab[i].costadd);

        /* If this ftgt the last reveived or the second last received ?? */
        if(timerVal(&(simuctrl->ftgttab[i].timerecv)) > lftgttime)
        {
            lftgttime = timerVal(&(simuctrl->ftgttab[i].timerecv));
            lftgtnum  = i;
        }
        else {
            if(timerVal(&(simuctrl->ftgttab[i].timerecv)) > sftgttime) {
                sftgttime = timerVal(&(simuctrl->ftgttab[i].timerecv));
            }
        }
    }


    /*
     * Put in ftgttimetab[] the date at which the cluster would have received
     * and add all the ftgt if the task was mapped on it.
     */
    for(i=simuctrl->bloktab[bloknum].ftgtnum; i<simuctrl->bloktab[bloknum+1].ftgtnum;i++)
    {
        if(i != lftgtnum) {
            timerSet(&(simuctrl->ftgttimetab[i]), lftgttime);
        }
        else {
            timerSetMax( &(simuctrl->ftgttimetab[i]), sftgttime );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Update the Fan In target structure
 *
 * Increment the contribution counter of the fan-in and integrate to the ftgt
 * area the new contribution.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The pointer to the symbolic matrix structure.
 *
 * @param[inout] simuctrl
 *          The pointer to the simulation structure. On exit, data regarding the
 *          computational unit pr are updated.
 *
 * @param[in] ftgtnum
 *          Index of the fanin target to update.
 *
 * @param[in] bloknum
 *          Index of the first off-diagonal block generating a contribution to
 *          the ftgtnum Fan In.
 *
 * @param[in] fbloknum
 *          Index of the off-diagonal block that is multiplied by blocknum to
 *          produce the update.
 *
 *******************************************************************************/
static inline void
simu_updateFtgt( const symbol_matrix_t *symbptr,
                       SimuCtrl        *simuctrl,
                       pastix_int_t     ftgtnum,
                       pastix_int_t     bloknum,
                       pastix_int_t     fbloknum )
{
    solver_ftgt_t *ftgt     = &(simuctrl->ftgttab[ftgtnum].ftgt);
    pastix_int_t  *infotab  = ftgt->infotab;
    symbol_blok_t *blokptr  = (symbptr->bloktab) + bloknum;
    symbol_blok_t *fblokptr = (symbptr->bloktab) + fbloknum;

    infotab[FTGT_CTRBNBR]++;

    /* Update ftgt dimensions to the maximum area covering all contributions */
    if( blokptr->frownum < infotab[FTGT_FCOLNUM] ) {
        infotab[FTGT_FCOLNUM] = blokptr->frownum;
    }

    if( blokptr->lrownum > infotab[FTGT_LCOLNUM] ) {
        infotab[FTGT_LCOLNUM] = blokptr->lrownum;
    }

    if( fblokptr->frownum < infotab[FTGT_FROWNUM] ) {
        infotab[FTGT_FROWNUM] = fblokptr->frownum;
    }

    if( fblokptr->lrownum > infotab[FTGT_LROWNUM] ) {
        infotab[FTGT_LROWNUM] = fblokptr->lrownum;
    }

    assert( (infotab[FTGT_LCOLNUM] - infotab[FTGT_FCOLNUM] + 1) > 0 );
    assert( (infotab[FTGT_LROWNUM] - infotab[FTGT_FROWNUM] + 1) > 0 );
}

/**
 *******************************************************************************
 *
 * @brief Simulate the task execution.
 *
 * Update the timers of the selected worker, as well as those of the current
 * cblk, and the targeted cblks by the update.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The pointer to the global blend control structure.
 *
 * @param[in] symbptr
 *          The pointer to the symbolic matrix structure.
 *
 * @param[inout] simuctrl
 *          The pointer to the simulation structure. On exit, data regarding the
 *          computational unit pr are updated.
 *
 * @param[in] tasknum
 *          The task index of the one, we want to simulate the execution.
 *
 *******************************************************************************
 *
 * @remark In this function, we use the standard [f|l]blocknum for first and
 * last bloknum, and facingcblk, facingblok for the facing block and column
 * block.
 *
 *******************************************************************************/
static inline void
simu_computeTask( const BlendCtrl       *ctrl,
                  const symbol_matrix_t *symbptr,
                        SimuCtrl        *simuctrl,
                        pastix_int_t     tasknum )
{
    pastix_int_t  i, j;
    pastix_int_t  cblknum;
    pastix_int_t  fbloknum;
    pastix_int_t  lbloknum;
    pastix_int_t  firstfacingblok;
    pastix_int_t  facingblok;
    pastix_int_t  facingcblk;
    pastix_int_t  local;
    pastix_int_t  ftgtnum;
    pastix_int_t  procnum;
    pastix_int_t  clustnum;
    SimuProc     *sproc;
    CostMatrix   *costmtx;

    cblknum  = simuctrl->tasktab[tasknum].cblknum;
    procnum  = simuctrl->ownetab[cblknum];
    clustnum = ctrl->core2clust[procnum];
    sproc    = &(simuctrl->proctab[procnum]);
    costmtx  = ctrl->costmtx;

    fbloknum = symbptr->cblktab[cblknum  ].bloknum;
    lbloknum = symbptr->cblktab[cblknum+1].bloknum;

    assert( (procnum >= ctrl->candtab[cblknum].fcandnum) &&
            (procnum <= ctrl->candtab[cblknum].lcandnum) );

    /* Add factorization time of the diagonal blok + cost of the TRSM operation on the cblk*/
    timerAdd(&(sproc->timer), costmtx->blokcost[fbloknum]);

    for(i=fbloknum+1; i<lbloknum; i++)
    {
        /* Add cost of the GEMM update related to this off-diagonal block */
        timerAdd(&(sproc->timer), costmtx->blokcost[i]);

        facingcblk = symbptr->bloktab[i].fcblknm;

        /*
         * If only one candidate cluster, we can consider the facingcblk as
         * local because it is an ancestor of the current cblk in the
         * elimination tree.
         */
        local = ( ctrl->candtab[facingcblk].fccandnum == ctrl->candtab[facingcblk].lccandnum ) ? 1 : 0;

        firstfacingblok = symbptr->cblktab[facingcblk].bloknum;

        for(j=i; j<lbloknum; j++)
        {
            /* TODO: symbolGetFacingBloknum is too expensive !! */
            facingblok = pastixSymbolGetFacingBloknum(symbptr, i, j, firstfacingblok, ctrl->ricar);

            /* If the couple (i, j) generates a contribution, applies it */
            if( facingblok >= 0 ) {
                pastix_int_t facingdiagblok;
                pastix_int_t facingtask;

                /* Decrease contributions on block and column block */
                simuctrl->cblktab[facingcblk].ctrbcnt--;
                simuctrl->bloktab[facingblok].ctrbcnt--;

                /* Checks */
                assert(simuctrl->cblktab[facingcblk].ctrbcnt >= 0);
                assert(simuctrl->bloktab[facingblok].ctrbcnt >= 0);

                /* Update to start next search from the last facing block */
                firstfacingblok = facingblok;

                facingdiagblok = symbptr->cblktab[facingcblk].bloknum;
                facingtask     = simuctrl->bloktab[facingdiagblok].tasknum;

                assert( facingcblk == simuctrl->tasktab[facingtask].cblknum );
                assert( facingtask < simuctrl->tasknbr );

                if(!local)
                {
                    ftgtnum = CLUST2INDEX(facingblok, clustnum);
                    simu_updateFtgt( symbptr, simuctrl, ftgtnum, i, j );

                    /* Update timer ready for receiver of the ftgt */
                    ftgtnum = CLUST2INDEX( facingdiagblok, clustnum );
                    timerSetMax( &(simuctrl->ftgttimetab[ftgtnum]),
                                 timerVal(&(sproc->timer)) );

                }
                else {

                    /* Update timer of the task (associated to the diagonal block) */
                    timerSetMax( &(simuctrl->tasktab[facingtask].time),
                                 timerVal(&(sproc->timer)) );
                }

                if( simuctrl->cblktab[facingcblk].ctrbcnt == 0 ) {
                    if (!local) {
                        simu_computeTaskReceiveTime(ctrl, symbptr, simuctrl, facingtask );
                    }

                    /* Put the task in the ready heap of its local candidat processor */
                    simu_putInAllReadyQueues( ctrl, simuctrl, facingtask );
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Push all task from future to ready
 *
 * This routine pushes all future task from the future task heap to the ready
 * one, if the time at which the task will be ready is already passed by the
 * computation unit.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The pointer to the global blend control structure.
 *
 * @param[inout] simuctrl
 *          The pointer to the simulation structure. On exit, data regarding the
 *          computational unit pr are updated.
 *
 * @param[in] procnum
 *          The computational unit index for which, the data need to be transfer
 *          from the future task heap to ready task heap if the computatuional
 *          unit timer is more advanced than the ready time of the tasks.
 *
 *******************************************************************************/
static inline void
simu_pushToReadyHeap( const BlendCtrl *ctrl,
                      SimuCtrl        *simuctrl,
                      pastix_int_t     procnum )
{
    SimuProc    *sproc;
    pastix_int_t tasknum;
    pastix_int_t cblknum;
    pastix_int_t clustnum;

    clustnum = ctrl->core2clust[procnum];
    sproc    = &(simuctrl->proctab[procnum]);

    /*
     * Move each task from future task heap to ready heap if the timer is
     * further in the future than the ready time
     */
    while( pqueueSize(sproc->futuretask) > 0 )
    {
        tasknum = pqueueRead(sproc->futuretask);

        if(! timerComp( &(sproc->timer),
                        &(simuctrl->ftgttimetab[CLUST2INDEX(simuctrl->tasktab[tasknum].bloknum, clustnum )])) )
        {
            tasknum = pqueuePop(sproc->futuretask);
            cblknum = simuctrl->tasktab[tasknum].cblknum;

            pqueuePush2(sproc->readytask, tasknum,
                        ctrl->candtab[cblknum].treelevel,
                        simuctrl->tasktab[tasknum].bloknum );
        }
        else {
            break;
        }
    }
}

/**
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Run the simulation to map the data on the nodes
 *
 * This routine simulates the numerical factorization to generate the static
 * scheduling and the final mapping of the column block onto the PaStiX
 * processes.
 *
 *******************************************************************************
 *
 * @param[inout] simuctrl
 *          The pointer to the simulation structure initialized by simuInit().
 *
 * @param[in] ctrl
 *          The pointer to the blend control structure which contains the
 *          required data, such as the worker distribution among the processes,
 *          the candidates array for each column block, and the cost of the
 *          computations.
 *
 * @param[in] symbptr
 *          The block symbol structure of the problem.
 *
 *******************************************************************************/
void
simuRun( SimuCtrl              *simuctrl,
         const BlendCtrl       *ctrl,
         const symbol_matrix_t *symbptr )
{

    pastix_int_t i, j, b;
    pastix_int_t cblknum, bloknum;
    /*pastix_int_t c;*/
    pastix_int_t pr;

#if defined(PASTIX_BLEND_GENTRACE)
    static volatile pastix_atomic_lock_t trace_lock = PASTIX_ATOMIC_UNLOCKED;
    char **procnames = NULL;

    if (ctrl->clustnum == 0)
    {
        pastix_atomic_lock( &trace_lock );
        char *tracename = NULL;
        int rc;

        pastix_gendirtemp( ctrl->dirtemp );
        if ( *(ctrl->dirtemp) == NULL ) {
            tracename = strdup( "blend" );
        }
        else {
            rc = asprintf( &tracename, "%s/blend", *(ctrl->dirtemp) );
        }

        setTraceType (PAJE);
        initTrace (tracename, 0, GTG_FLAG_NONE);
        free(tracename);

        addContType ("CT_Appli", "0",        "Application" );
        addContType ("CT_P",     "CT_Appli", "Process"     );
        addContType ("CT_T",     "CT_P",     "Thread"      );
        addStateType("ST_TS",    "CT_T",     "Thread State");

        addLinkType ("LT_TL", "Split Event Link", "CT_P", "CT_T", "CT_T");

        /* Create root container of the application */
        addContainer (0.00000, "Appli", "CT_Appli", "0", "PaStiX Blend Simulation", "");

        /* Add all possible states */
        addEntityValue ("Wait", "ST_TS", "Waiting",   GTG_LIGHTGREY);
        addEntityValue ("Comp", "ST_TS", "Computing", GTG_RED);

        /* Add each process and thread */
        procnames = (char**) malloc ( ctrl->total_nbthrds * sizeof(char*) );
        pr = 0;
        for (i=0; i<ctrl->clustnbr; i++) {
            char *clustname;
            char *clustalias;

            rc = asprintf( &clustname, "Process %02d", (int)i); assert(rc!=-1);
            rc = asprintf( &clustalias, "P%d", (int)i); assert(rc!=-1);
            addContainer (0.00000, clustalias, "CT_P", "Appli", clustname, "");

            for (j=0; j<ctrl->local_nbthrds; j++, pr++) {
                char *procname;
                char *procalias;

                rc = asprintf( &procname, "Thread %02d", (int)pr); assert(rc!=-1);
                rc = asprintf( &procalias, "T%d", (int)pr); assert(rc!=-1);
                addContainer (0.00000, procalias, "CT_T", clustname, procname, "");

                procnames[pr] = procalias;
                free(procname);
            }

            free(clustname); free(clustalias);
        }
    }
#endif /* defined(PASTIX_BLEND_GENTRACE) */

    /* Compute number of contributions per blocks, cblks, tasks */
    simu_computeBlockCtrbNbr( symbptr, simuctrl, ctrl->ricar );

    if ( ctrl->iparm[IPARM_VERBOSE] > 4 ) {
        simu_printBlockCtrbNbr( ctrl, symbptr, simuctrl );
    }

    /*
     * All ready tasks are put in the task heaps of their respective candidates
     */
    for(i=0;i<symbptr->cblknbr;i++)
    {
        pastix_int_t tasknum;
        if(simuctrl->cblktab[i].ctrbcnt == 0)
        {
            tasknum = simuctrl->bloktab[symbptr->cblktab[i].bloknum].tasknum;
            assert(ctrl->candtab[i].treelevel < 0);

            if( ctrl->costlevel ) {
                assert(ctrl->candtab[i].costlevel <= 0);
            }

            assert(simuctrl->tasktab[tasknum].cblknum == i);
            //assert(ctrl->candtab[i].cblktype == CBLK_1D);

            simu_putInAllReadyQueues( ctrl, simuctrl, tasknum );
        }
    }

    /*
     * Run simulation and map the task onto a single candidate
     */
    while(1)
    {
        SimuTask    *task;
        pastix_int_t clustnum;

        /* Get the next earlier task index and the processor on which it is mapped */
        i = simu_getNextTaskNextProc(ctrl, simuctrl, &pr);

        /* No more tasks */
        if( i == -1 ) {
            break;
        }

        task     = &(simuctrl->tasktab[i]);
        bloknum  = task->bloknum;
        cblknum  = task->cblknum;
        clustnum = ctrl->core2clust[pr];

        assert(cblknum < symbptr->cblknbr);
        assert(bloknum < symbptr->bloknbr);

        /* Make sure the cblk is not already atributed to someone and give it to the selected proc */
        assert( simuctrl->ownetab[cblknum] < 0 );
        simuctrl->ownetab[cblknum] = pr;
        for(j = symbptr->cblktab[cblknum].bloknum;
            j < symbptr->cblktab[cblknum+1].bloknum; j++)
        {
            simuctrl->bloktab[j].ownerclust = clustnum;
        }
        task->prionum = simuctrl->clustab[clustnum].prionum;
        simuctrl->clustab[clustnum].prionum++;

        /* Add task to the selected processor list */
        extendint_Add(simuctrl->proctab[pr].tasktab, i);

        /* Backup which cluster will get the data for the second run of proportionnal mapping */
        ctrl->candtab[cblknum].cluster = clustnum;

        /*
         * Compute the time at which each proc cand will have added its ftgt and
         * received block target if the task is mapped on
         */
        if( ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum ) {
            /*
             * All contributions come from the same node
             * Time do not depend on the reception of a ftgt
             */
            timerSetMax( TIMER(pr), timerVal(&(task->time)) );
        }
        else {
            /*
             * Contributions might come from different nodes
             * Time depends on the reception of a ftgt
             */
            timerSetMax( TIMER(pr),
                         timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(bloknum, clustnum)])) );
        }

#if defined(PASTIX_BLEND_GENTRACE)
        if (ctrl->clustnum == 0) {
            assert( (procnames != NULL) && (pr < ctrl->total_nbthrds) );
            assert( procnames[pr] != NULL );
            setState( timerVal( TIMER(pr) ), "ST_TS", procnames[pr], "Comp" );
        }
#endif

        /*
         * Fill some fanintarget info (task of type E2 does not have any ftgt)
         */
        if(simuctrl->bloktab[bloknum].ftgtnum < simuctrl->bloktab[bloknum+1].ftgtnum)
        {
            /* Task with several cand cluster */
            for(b=bloknum; b<symbptr->cblktab[cblknum+1].bloknum; b++)
            {
                for(j=simuctrl->bloktab[b].ftgtnum; j<simuctrl->bloktab[b+1].ftgtnum; j++)
                {
                    if( (simuctrl->ftgttab[j].ftgt.infotab[FTGT_CTRBNBR] >0)
                        && (j != CLUST2INDEX(b, clustnum)))
                    {
                        simuctrl->ftgttab[j].clustnum = INDEX2CLUST(j, b);
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_PRIONUM] = task->prionum;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_PROCDST] = pr;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_BLOKDST] = b;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_TASKDST] = simuctrl->bloktab[bloknum].tasknum;
#if defined(PASTIX_WITH_STARPU_DIST)
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_GCBKDST] = simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].cblknum;
#endif
                        extendint_Add(&(simuctrl->clustab[INDEX2CLUST(j,b)].ftgtsend[clustnum]), j);

                        simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].ftgtcnt++;

                        if (clustnum == ctrl->clustnum) {
                            simuctrl->ftgtcnt++;
                        }
                    }
                }
            }
            simuctrl->ftgtprio++;
        }
        else {
            assert(ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum);
        }

        /* Simulate the task computation */
        simu_computeTask( ctrl, symbptr, simuctrl, i );

#if defined(PASTIX_BLEND_GENTRACE)
        if (ctrl->clustnum == 0) {
            assert( (procnames != NULL) && (pr < ctrl->total_nbthrds) );
            assert( procnames[pr] != NULL );
            setState( timerVal( TIMER(pr) ), "ST_TS", procnames[pr], "Wait" );
        }
#endif
        simu_pushToReadyHeap(ctrl, simuctrl, pr);
    }

    /* Compute maximum time */
    {
        double maxtime = 0;
        for(pr=0; pr<ctrl->total_nbcores; pr++)
        {
            if(timerVal(TIMER(pr)) > maxtime) {
                maxtime = timerVal(TIMER(pr));
            }
        }
        set_dparm(ctrl->dparm, DPARM_PRED_FACT_TIME, maxtime);
    }

#if defined(PASTIX_BLEND_GENTRACE)
    if (ctrl->clustnum == 0) {
        assert( procnames != NULL );
        for(pr=0; pr<ctrl->total_nbthrds; pr++) {
            free(procnames[pr]);
        }
        free(procnames);

        endTrace();
        pastix_atomic_unlock( &trace_lock );
    }
#endif

#if defined(PASTIX_DEBUG_BLEND)
    for(i=0;i<simuctrl->cblknbr;i++) {
        /* Check valid for 1D distribution only */
        assert( simuctrl->ownetab[i] >= 0 );
    }
    for(i=0;i<symbptr->bloknbr;i++) {
        assert( simuctrl->bloktab[i].ownerclust >= 0 );
    }
#endif
}
