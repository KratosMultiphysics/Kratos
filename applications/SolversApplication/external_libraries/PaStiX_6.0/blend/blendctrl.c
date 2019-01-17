/**
 *
 * @file blendctrl.c
 *
 * PaStiX analyse control parameters function.
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"
#include "elimintree.h"
#include "cost.h"
#include "extendVector.h"
#include "cand.h"
#include "blendctrl.h"
#include "perf.h"

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_ctrl
 *
 * @brief Return the communication cost between two cores
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The Blend control data structure that holds the cluster architecture
 *
 * @param[in] clustsrc
 *          The source cluster of the communication
 *
 * @param[in] clustdst
 *          The remote cluster of the communication
 *
 * @param[in] sync_comm_nbr
 *          The number of simultaneous communication between clustsrc and clustdst
 *
 * @param[out] startup
 *          On exit, holds the startup/latency cost of the communication between
 *          the two nodes.
 *
 * @param[out] bandwidth
 *          On exit, holds the bandwidth of the communication between the two
 *          nodes.
 *
 *******************************************************************************/
void
getCommunicationCosts( const BlendCtrl *ctrl,
                       pastix_int_t     clustsrc,
                       pastix_int_t     clustdst,
                       pastix_int_t     sync_comm_nbr,
                       double          *startup,
                       double          *bandwidth)
{
    assert((clustsrc >= 0) && (clustsrc < ctrl->clustnbr));
    assert((clustdst >= 0) && (clustdst < ctrl->clustnbr));
    assert((sync_comm_nbr > 0) && (sync_comm_nbr <= ctrl->clustnbr));

    if(clustsrc == clustdst)
    {
        *startup   = 0.;
        *bandwidth = 0.;
        return;
    }

    /* Shared Memory */
    if( ctrl->clust2smp[clustsrc] == ctrl->clust2smp[clustdst] )
    {
        switch (sync_comm_nbr)
        {
        case 1:
        case 2:
            *startup   = SHARED_STARTUP_1;
            *bandwidth = SHARED_BANDWIDTH_1;
            return;
        case 3:
        case 4:
            *startup   = SHARED_STARTUP_2;
            *bandwidth = SHARED_BANDWIDTH_2;
            return;
        case 5:
        case 6:
        case 7:
        case 8:
            *startup   = SHARED_STARTUP_4;
            *bandwidth = SHARED_BANDWIDTH_4;
            return;
        default:
            *startup   = SHARED_STARTUP_8;
            *bandwidth = SHARED_BANDWIDTH_8;
            return;
        }
    }
    else
    {
        switch (sync_comm_nbr)
        {
        case 1:
        case 2:
            *startup   = CLUSTER_STARTUP_1;
            *bandwidth = CLUSTER_BANDWIDTH_1;
            return;
        case 3:
        case 4:
            *startup   = CLUSTER_STARTUP_2;
            *bandwidth = CLUSTER_BANDWIDTH_2;
            return;
        case 5:
        case 6:
        case 7:
        case 8:
            *startup   = CLUSTER_STARTUP_4;
            *bandwidth = CLUSTER_BANDWIDTH_4;
            return;
        default:
            *startup   = CLUSTER_STARTUP_8;
            *bandwidth = CLUSTER_BANDWIDTH_8;
            return;
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_ctrl
 *
 * @brief Initialize the Blend control structure.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure of the problem instance.  Integer
 *          parameters that are used to build the Blend control
 *          structure. IPARM_ABS, IPARM_TASKS2D_LEVEL, IPARM_TASKS2D_WIDTH,
 *          IPARM_COMPRESS_WHEN, IPARM_COMPRESS_MIN_WIDTH, IPARM_INCOMPLETE,
 *          IPARM_MAX_BLOCKSIZE, IPARM_MIN_BLOCKSIZE, IPARM_THREAD_NBR, and
 *          IPARM_VERBOSE are used in this function.
 *
 * @param[inout] ctrl
 *          The Blend control data structure to initialize.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
blendCtrlInit( pastix_data_t *pastix_data,
               BlendCtrl     *ctrl )
{
    pastix_int_t  procnum = pastix_data->inter_node_procnum;
    pastix_int_t  procnbr = pastix_data->inter_node_procnbr;
    pastix_int_t *iparm   = pastix_data->iparm;
    double       *dparm   = pastix_data->dparm;
    pastix_int_t  local_coresnbr = iparm[IPARM_THREAD_NBR];
    pastix_int_t  local_thrdsnbr = iparm[IPARM_THREAD_NBR];
    pastix_int_t  i;

    /* Check parameters */
    if( ctrl == NULL )
    {
        errorPrint("blendCtrlInit: Illegal ctrl parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( procnum < 0 )
    {
        errorPrint("blendCtrlInit: Illegal procnum parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( procnbr < 1 )
    {
        errorPrint("blendCtrlInit: Illegal procnbr parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( local_coresnbr < 1 )
    {
        errorPrint("blendCtrlInit: Illegal local_coresnbr parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( local_thrdsnbr < 1 )
    {
        errorPrint("blendCtrlInit: Illegal local_thrdsnbr parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( procnum >= procnbr )
    {
        errorPrint("blendCtrlInit: Incompatible values of procnum(%d) and procnbr (%d)\n",
                   (int) procnum, (int) procnbr);
        return PASTIX_ERR_BADPARAMETER;
    }
    if( ctrl == NULL )
    {
        errorPrint("blendCtrlInit: Illegal ctrl parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Initialize options */
    ctrl->count_ops = 1;
#if defined(PASTIX_DEBUG_BLEND)
    ctrl->debug     = 1;
#else
    ctrl->debug     = 0;
#endif
    ctrl->timer     = 1;
    ctrl->ricar     = iparm[IPARM_INCOMPLETE];
    ctrl->leader    = 0;

    /* Proportional Mapping options */
    ctrl->allcand     = iparm[IPARM_ALLCAND];
    ctrl->nocrossproc = 0;
    ctrl->costlevel   = 1;

    /* Spliting options */
    ctrl->blcolmin = iparm[IPARM_MIN_BLOCKSIZE];
    ctrl->blcolmax = iparm[IPARM_MAX_BLOCKSIZE];
    ctrl->abs      = iparm[IPARM_ABS];
    ctrl->updatecandtab = 0; /* TODO: Set to 1 to match former version of PaStiX, move back to 0 */
    if(ctrl->blcolmin > ctrl->blcolmax)
    {
        errorPrint("Parameter error : blocksize max < blocksize min (cf. iparm.txt).");
        assert(ctrl->blcolmin <= ctrl->blcolmax);
    }

    /* 2D options */
    ctrl->level_tasks2d = iparm[IPARM_TASKS2D_LEVEL];
    ctrl->width_tasks2d = iparm[IPARM_TASKS2D_WIDTH];

    /* Save iparm for other options */
    ctrl->iparm = iparm;
    ctrl->dparm = dparm;
    ctrl->dirtemp = &(pastix_data->dirtemp);

    /*
     * Initialize architecture description
     */
    /* Id and number of MPI processes */
    ctrl->clustnum = procnum;
    ctrl->clustnbr = procnbr;

    /* Local informations */
    ctrl->local_nbcores = local_coresnbr;
    ctrl->local_nbthrds = local_thrdsnbr;
    ctrl->local_nbctxts = ctrl->local_nbthrds;

    /* Total information (should require a MPI_Reduce if different informations on each node) */
    ctrl->total_nbcores = ctrl->local_nbcores * procnbr;
    ctrl->total_nbthrds = ctrl->local_nbthrds * procnbr;

    /* Create the array of associativity bewteen MPI process ids and SMP node ids */
    /* Rq: We could use a MPI reduction for irregular pattern                     */
    /* TODO: insert back the number of MPI processes per node                     */
    MALLOC_INTERN(ctrl->clust2smp, ctrl->clustnbr, pastix_int_t);
    for(i=0; i < ctrl->clustnbr; i++) {
        ctrl->clust2smp[i] = i;
    }

    /* Create the array of associativity bewteen core ids and MPI process ids */
    /* Rq: We could use a MPI reduction for irregular pattern             */
    MALLOC_INTERN(ctrl->core2clust, ctrl->total_nbcores, pastix_int_t);
    for(i=0; i < ctrl->total_nbcores; i++) {
        ctrl->core2clust[i] = i / ctrl->local_nbcores;
    }

    ctrl->etree   = NULL;
    ctrl->costmtx = NULL;
    ctrl->candtab = NULL;

#ifdef PASTIX_DYNSCHED
    MALLOC_INTERN(ctrl->btree, 1, BubbleTree);
#endif

    return PASTIX_SUCCESS;
}


/**
 *******************************************************************************
 *
 * @ingroup blend_dev_ctrl
 *
 * @brief Finalize the Blend control structure.
 *
 *******************************************************************************
 *
 * @param[inout] ctrl
 *          The Blend control data structure to free.
 *
 *******************************************************************************/
void
blendCtrlExit(BlendCtrl *ctrl)
{
    if(ctrl->clust2smp) {
        memFree_null(ctrl->clust2smp);
    }
    if(ctrl->core2clust) {
        memFree_null(ctrl->core2clust);
    }
    if(ctrl->candtab) {
        candExit( ctrl->candtab );
        ctrl->candtab = NULL;
    }
}
