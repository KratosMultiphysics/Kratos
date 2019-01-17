/**
 *
 * @file simu.h
 *
 * PaStiX simulation structure and basic functions.
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
#ifndef _simu_h_
#define _simu_h_

#include "simu_timer.h"

/**
 * @brief Process structure for the simulation.
 */
typedef struct simu_cluster_s {
    pastix_int_t     fprocnum;   /**< Global index of the first processor belonging to the cluster (Check is it is not redundant) */
    pastix_int_t     lprocnum;   /**< Global index of the last processor belonging to the cluster (Check is it is not redundant)  */
    ExtendVectorINT *ftgtsend;   /**< Arrays of ftgt sent by this proc (one vector per processor)                                 */
    pastix_int_t     prionum;    /**< Counter to order tasks on one cluster                                                       */
} SimuCluster;

/**
 * @brief Thread structure for the simulation.
 */
typedef struct simu_proc_s {
    SimuTimer        timer;      /**< Simulated clock of the processor                                  */
    pastix_queue_t  *readytask;  /**< Heap of tasks ready to be executed                                */
    pastix_queue_t  *futuretask; /**< Heap of tasks ready to be executed in a near future (after timer) */
    ExtendVectorINT *tasktab;    /**< Vector to store tasks affected to the candidate                   */
} SimuProc;

/**
 * @brief Fan-in structure for the simulation.
 */
typedef struct simu_ftgt_s {
    solver_ftgt_t ftgt;         /**< Fan-in informations                            */
    pastix_int_t  clustnum;     /**< Cluster sending the contribution               */
    SimuTimer     timerecv;     /**< Simulated clock of the reception time          */
    double        costsend;     /**< Cost to send the contribution                  */
    double        costadd;      /**< Cost to add the contribution to its final cblk */
} SimuFtgt;

/**
 * @brief Column block structure for the simulation.
 */
typedef struct simu_cblk_s {
    pastix_int_t ctrbcnt;       /**< Counter of remaining contributions for the cblk */
} SimuCblk;

/**
 * @brief Block structure for the simulation.
 */
typedef struct simu_blok_s {
    pastix_int_t tasknum;       /**< Task index opeating on this block (stored per block for 2D computations)   */
    pastix_int_t ftgtnum;       /**< Index of the first fanin destinated to this
                                     block in the ftgttab. This index is also used to find the first cblk timer
                                     (one per cand proc) in the timetab array                                   */
    pastix_int_t ctrbcnt;       /**< Counter of remaining contributions                                         */
    int          fccandnum;     /**< First candidate that is attributed to the cblk of the block                */
    int          ownerclust;    /**< Processor on which the block is distributed                                */
} SimuBlok;

/**
 * @brief Task structure for the simulation.
 */
typedef struct simu_task_s {
    pastix_int_t prionum;       /**< priority of the task                                      */
    pastix_int_t cblknum;       /**< Number of the cblknum the task deals with                 */
    pastix_int_t bloknum;       /**< number of the block that the task deals with              */
    pastix_int_t bloknum2;      /**< */
    pastix_int_t facebloknum;   /**< Number of the facing block for E2                         */
    SimuTimer    time;          /**< Time the task is ready if it doesn't need message         */
    pastix_int_t mesglen;       /**< Time to send the block target                             */
    double       cost;          /**< Cost of the task                                          */
    pastix_int_t ctrbcnt;       /**< nbr ftgt + le btgt (+ E1 pret si E2)                      */
    pastix_int_t ftgtcnt;       /**< nbr of contrib from fan-in target                         */
    pastix_int_t tasknext;      /**< chainage des E1 ou E2, si fin = -1 => liberer les btagptr */
} SimuTask;

/**
 * @brief Control structure for the simulation.
 */
typedef struct simuctrl_s {
    pastix_int_t  cblknbr;      /**< Number of cblk                                            */
    pastix_int_t  ftgtprio;     /**< Priority to assign to current ftgts                       */
    pastix_int_t  tasknbr;      /**< Number of tasks                                           */
    pastix_int_t  ftgtcnt;      /**< Number of received communication                          */
    SimuTask     *tasktab;      /**< SimuTask vector                                           */
    SimuProc     *proctab;      /**< Virtual processor tab                                     */
    SimuCluster  *clustab;      /**< Virtual cluster tab                                       */
    pastix_int_t *ownetab;      /**< Vector containing the distribution of the diagonal blok   */
    SimuCblk     *cblktab;      /**< SimuCblk vector                                           */
    SimuBlok     *bloktab;      /**< SimuBlok vector                                           */
    SimuFtgt     *ftgttab;      /**< Vector containing the fan in target                       */
    pastix_int_t  ftgtnbr;      /**< The number of fan-in contribution                         */
    SimuTimer    *ftgttimetab;  /**< Vector containing a timer for each cluster on each ftgt   */
} SimuCtrl;


pastix_int_t simuInit        ( SimuCtrl *, const symbol_matrix_t *, const Cand *, pastix_int_t, pastix_int_t );
pastix_int_t simuRealloc     ( SimuCtrl *, pastix_int_t, pastix_int_t );
void         simuExit        ( SimuCtrl *, pastix_int_t, pastix_int_t, pastix_int_t );
void         simuTaskBuild   ( SimuCtrl *, const symbol_matrix_t * );

#define CLUST2INDEX(n,c) ((c) + simuctrl->bloktab[n].ftgtnum - simuctrl->bloktab[n].fccandnum)
#define INDEX2CLUST(r,s) ((r) - simuctrl->bloktab[s].ftgtnum + simuctrl->bloktab[s].fccandnum)
#define TIMER(pr)        (&(simuctrl->proctab[pr].timer))

#endif /* _simu_h_ */

/**
 *@}
 */
