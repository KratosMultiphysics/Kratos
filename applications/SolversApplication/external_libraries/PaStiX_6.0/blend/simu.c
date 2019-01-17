/**
 *
 * @file simu.c
 *
 * PaStiX simulation basic functions.
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
 * @addtogroup blend_dev_simu
 * @{
 *
 **/
#include <stdio.h>
#include <strings.h>
#include <math.h>

#include "common.h"
#include "symbol.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "queue.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "solver.h"
#include "simu.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the simulation structures.
 *
 *******************************************************************************
 *
 * @param[inout] simuctrl
 *          The allocated pointer to main control structure for the simulation.
 *
 * @param[in] symbptr
 *          The pointer to the symbol matrix structure studied.
 *
 * @param[in] candtab
 *          The pointer to the candidate information associated to the symbol
 *          structure and that will guide the simulation.
 *
 * @param[in] clustnbr
 *          The number of PaStiX processes involved in the simulation.
 *
 * @param[in] procnbr
 *          The total number of workers involved in the simulation.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
pastix_int_t
simuInit( SimuCtrl              *simuctrl,
          const symbol_matrix_t *symbptr,
          const Cand            *candtab,
          pastix_int_t           clustnbr,
          pastix_int_t           procnbr )
{
    pastix_int_t i, j;
    pastix_int_t p;
    pastix_int_t ftgtcur;
    pastix_int_t candnbr;
    pastix_int_t step;
    pastix_int_t cblknbr = symbptr->cblknbr;
    pastix_int_t bloknbr = symbptr->bloknbr;

    simuctrl->cblknbr  = cblknbr;
    simuctrl->ftgtprio = 0;
    simuctrl->tasktab  = NULL;
    simuctrl->ftgtcnt  = 0;

    /* Processor initialisation */
    MALLOC_INTERN(simuctrl->proctab, procnbr, SimuProc);
    for(i=0;i<procnbr;i++)
    {
        timerSet(TIMER(i), 0.0); /* for paragraph numeric tolerance */
        MALLOC_INTERN(simuctrl->proctab[i].futuretask, 1, pastix_queue_t);
        MALLOC_INTERN(simuctrl->proctab[i].readytask,  1, pastix_queue_t);
        pqueueInit(simuctrl->proctab[i].futuretask, 100);
        pqueueInit(simuctrl->proctab[i].readytask,  100);

        MALLOC_INTERN(simuctrl->proctab[i].tasktab, 1, ExtendVectorINT);
        extendint_Init(simuctrl->proctab[i].tasktab, bloknbr/procnbr + 1);
    }

    /* Cluster initialization */
    MALLOC_INTERN(simuctrl->clustab, clustnbr, SimuCluster);
    step = procnbr / clustnbr;
    for(i=0;i<clustnbr;i++)
    {
        simuctrl->clustab[i].fprocnum = i*step;
        simuctrl->clustab[i].lprocnum = simuctrl->clustab[i].fprocnum + step - 1;
        MALLOC_INTERN(simuctrl->clustab[i].ftgtsend, clustnbr, ExtendVectorINT);
        simuctrl->clustab[i].prionum  = 0;
        for(p=0;p<clustnbr;p++) {
            extendint_Init(&(simuctrl->clustab[i].ftgtsend[p]), cblknbr/(2*clustnbr)+1);
        }
    }
    simuctrl->clustab[clustnbr-1].lprocnum = procnbr-1;

    MALLOC_INTERN(simuctrl->ownetab, cblknbr, pastix_int_t);

    /* Affect a negative value to cblk not mapped */
    for(i=0;i<cblknbr;i++) {
        simuctrl->ownetab[i] = -1;
    }

    MALLOC_INTERN(simuctrl->cblktab, cblknbr+1, SimuCblk);
    MALLOC_INTERN(simuctrl->bloktab, bloknbr+1, SimuBlok);
    ftgtcur = 0;

    for(i=0;i<cblknbr;i++)
    {
        candnbr = candtab[i].lccandnum - candtab[i].fccandnum + 1;
        simuctrl->cblktab[i].ctrbcnt = 0;

        for(j=symbptr->cblktab[i].bloknum;j<symbptr->cblktab[i+1].bloknum;j++)
        {
            simuctrl->bloktab[j].ftgtnum    = ftgtcur;
            simuctrl->bloktab[j].tasknum    = -1;
            simuctrl->bloktab[j].fccandnum  = candtab[i].fccandnum;
            simuctrl->bloktab[j].ctrbcnt    = 0;
            simuctrl->bloktab[j].ownerclust = -1;
            /*if(candnbr > 1)*/
            ftgtcur += candnbr;
        }
    }
    /* One extracblk for avoiding side effect */
    simuctrl->bloktab[bloknbr].ftgtnum = ftgtcur;
    simuctrl->ftgtnbr = ftgtcur;

    if(simuctrl->ftgtnbr > 0)
    {
        /* Allocate and Initialize the timer for the reception of each ftgt on a candidate cluster */
        MALLOC_INTERN(simuctrl->ftgttimetab, simuctrl->ftgtnbr, SimuTimer);
        for(i=0;i<simuctrl->ftgtnbr;i++) {
            timerSet(&(simuctrl->ftgttimetab[i]), 0.0);
        }

        MALLOC_INTERN(simuctrl->ftgttab, ftgtcur, SimuFtgt);
        for(i=0;i<simuctrl->ftgtnbr;i++)
        {
            simuctrl->ftgttab[i].clustnum = -1;
            timerSet(&(simuctrl->ftgttab[i].timerecv), 0.0);
            simuctrl->ftgttab[i].costsend = 0.0;
            simuctrl->ftgttab[i].costadd  = 0.0;
            bzero(simuctrl->ftgttab[i].ftgt.infotab,FTGT_MAXINFO*sizeof(pastix_int_t));
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_FCOLNUM] = PASTIX_INT_MAX;
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_FROWNUM] = PASTIX_INT_MAX;
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR] = 0;
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBCNT] = 0;
        }
    }
    else
    {
        simuctrl->ftgttab     = NULL;
        simuctrl->ftgttimetab = NULL;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Reallocate the simulation structures to compact them.
 *
 * All data that were globally allocated and replicated on every node is freed,
 * and only the local information are allocated. This is used only when a second
 * proportionnal mapping with local information is performed for the dynamic
 * scheduling.
 *
 *******************************************************************************
 *
 * @param[inout] simuctrl
 *          The simulation structure to reallocate.
 *
 * @param[in] procnbr
 *          The total number of workers.
 *
 * @param[in] local_nbthrds
 *          The number of local workers per process.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
pastix_int_t
simuRealloc(SimuCtrl     *simuctrl,
            pastix_int_t  procnbr,
            pastix_int_t  local_nbthrds)
{
    pastix_int_t i;

    /* Free processor structure */
    for(i=0;i<procnbr;i++)
    {
        pqueueExit    (simuctrl->proctab[i].readytask);
        memFree_null  (simuctrl->proctab[i].readytask);
        pqueueExit    (simuctrl->proctab[i].futuretask);
        memFree_null  (simuctrl->proctab[i].futuretask);
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null  (simuctrl->proctab[i].tasktab);
    }
    memFree_null(simuctrl->proctab);

    /* Initialization for local threads */
    MALLOC_INTERN(simuctrl->proctab, local_nbthrds, SimuProc);
    for(i=0;i<local_nbthrds;i++)
    {
        MALLOC_INTERN(simuctrl->proctab[i].tasktab, 1, ExtendVectorINT);
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Free the simulation structure.
 *
 *******************************************************************************
 *
 * @param[inout] simuctrl
 *          The simulation structure to free.
 *
 * @param[in] clustnbr
 *          The total number of PaStiX processes.
 *
 * @param[in] procnbr
 *          The total number of workers.
 *
 * @param[in] local_nbthrds
 *          The number of local workers per process.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
void
simuExit(SimuCtrl     *simuctrl,
         pastix_int_t  clustnbr,
         pastix_int_t  procnbr,
         pastix_int_t  local_nbthrds)
{
    pastix_int_t i,j;
    (void)local_nbthrds; (void)procnbr;

#ifndef PASTIX_DYNSCHED
    for(i=0;i<procnbr;i++)
    {
        pqueueExit    (simuctrl->proctab[i].readytask);
        memFree_null  (simuctrl->proctab[i].readytask);
        pqueueExit    (simuctrl->proctab[i].futuretask);
        memFree_null  (simuctrl->proctab[i].futuretask);
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null  (simuctrl->proctab[i].tasktab);
    }
#else
    for(i=0;i<local_nbthrds;i++)
    {
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null(simuctrl->proctab[i].tasktab);
    }
#endif

    for(i=0;i<clustnbr;i++)
    {
        for(j=0;j<clustnbr;j++) {
            extendint_Exit(&(simuctrl->clustab[i].ftgtsend[j]));
        }
        memFree_null(simuctrl->clustab[i].ftgtsend);
    }

    if(simuctrl->ftgttab != NULL)
    {
        memFree_null(simuctrl->ftgttab);
        memFree_null(simuctrl->ftgttimetab);
    }
    memFree_null(simuctrl->tasktab);
    memFree_null(simuctrl->proctab);
    memFree_null(simuctrl->clustab);
    memFree_null(simuctrl->ownetab);
    memFree_null(simuctrl->cblktab);
    memFree_null(simuctrl->bloktab);
    memFree_null(simuctrl);
}

/**
 *@}
 */
