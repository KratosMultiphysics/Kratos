/**
 *
 * @file simu_task.c
 *
 * PaStiX simulation task basic functions.
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
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "symbol.h"
#include "queue.h"
#include "extendVector.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "blendctrl.h"
#include "solver.h"
#include "simu.h"

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_simu
 *
 * @brief Initialize the tasktab array of the simulation structure.
 *
 *******************************************************************************
 *
 * @param[inout] simuctrl
 *          The main simulation structure. On exit, the tasktab array is built
 *          and initialized with default values.
 *
 * @param[in] symbptr
 *          The pointer to the symbol matrix structure studied.
 *
 *******************************************************************************/
void
simuTaskBuild( SimuCtrl              *simuctrl,
               const symbol_matrix_t *symbptr )
{
    pastix_int_t i, j;
    pastix_int_t tasknbr = 0;
    SimuTask *task = NULL;

    /* One task per cblk */
    simuctrl->tasknbr = symbptr->cblknbr;

    MALLOC_INTERN( simuctrl->tasktab, simuctrl->tasknbr, SimuTask );
    assert( simuctrl->tasktab );

    task = simuctrl->tasktab;
    tasknbr = 0;
    for(i=0;i<symbptr->cblknbr;i++)
    {
        task->prionum  = -1;
        task->cblknum  = i;
        task->bloknum  = symbptr->cblktab[i].bloknum;
        task->bloknum2 = -1;
        task->ctrbcnt  = 0;
        task->ftgtcnt  = 0;
        task->facebloknum = -1;
        task->cost     = -1;
        timerSet(&(task->time), 0.0);
        task->mesglen  = 0.0;
        task->tasknext = -1;
        for(j = symbptr->cblktab[i].bloknum;
            j < symbptr->cblktab[i+1].bloknum; j++ )
        {
            simuctrl->bloktab[j].tasknum = tasknbr;
        }
        tasknbr++;
        task++;
    }
}
