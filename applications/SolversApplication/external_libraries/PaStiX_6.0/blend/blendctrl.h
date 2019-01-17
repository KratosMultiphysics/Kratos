/**
 *
 * @file blendctrl.h
 *
 * PaStiX analyse control parameters structure.
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_ctrl
 * @{
 *    This module handles all option to parameterize the final analyze step that
 *    performs proportional mapping and the generation of the final solver
 *    matrix structure per MPI process.
 *
 **/
#ifndef _blendctrl_h_
#define _blendctrl_h_

/**
 * @brief The type and structure definitions.
*/
typedef struct blendctrl_s {
    pastix_int_t    count_ops ;      /**< Print costs in term of number of elementary operations            */
    pastix_int_t    debug ;          /**< Make additional checks after each step                            */
    pastix_int_t    timer;           /**< Print execution times                                             */
    pastix_int_t    ooc;             /**< Enable the out-of-core version of Pastix (Option unused for now)  */
    pastix_int_t    ricar;           /**< Enable the ILU(k) dedicated steps                                 */
    pastix_int_t    leader;          /**< Leader for sequential tasks                                       */

    /**
     * @name Proportional Mapping
     * @{
     */
    pastix_int_t    allcand;         /**< All processors are candidate for each cblk                        */
    pastix_int_t    nocrossproc;     /**< Forbid a processor to be candidate in two
                                          different branches shared with different partners                 */
    pastix_int_t    costlevel;       /**< Enable/disable computation and use of subtree cost                */

    /**
     * @}
     * @name Symbol split
     * @{
     */
    pastix_int_t    blcolmin ;       /**< Minimun number of columns for a good use of BLAS primitives       */
    pastix_int_t    blcolmax;        /**< Maximum number of columns for a good use of BLAS primitives       */
    pastix_int_t    abs;             /**< Adaptative block size:
                                            - 0, all block are cut to blcolmin
                                            - >0, try to make (ncand*abs) cblk                              */
    pastix_int_t    updatecandtab;   /**< Update the candtab array after splitting, otherwise each new cblk
                                          has the same properties as the original one                       */

    /**
     * @}
     * @name 2D
     * @{
     */
    pastix_int_t    level_tasks2d;   /**< Level to shift from 1D to 2D. Automaticaly computed if < 0, only 1D if 0 */
    pastix_int_t    width_tasks2d;   /**< Minimal width to consider a cblk 2D if autolevel (level_tasks2d < 0)     */

    /**
     * @}
     * @name Architecture
     * @{
     */
    pastix_int_t    clustnum;        /**< Id of current MPI process                                         */
    pastix_int_t    clustnbr;        /**< Number of MPI processes                                           */
    pastix_int_t    total_nbcores;   /**< Total number of physical cores used for the simulation            */
    pastix_int_t    total_nbthrds;   /**< Total number of threads used for the simulation                   */
    pastix_int_t    local_nbcores;   /**< Local number of physical cores used by the current MPI process    */
    pastix_int_t    local_nbthrds;   /**< Local number of threads used by the current MPI process           */
    pastix_int_t    local_nbctxts;   /**< Local number of contexts (used for dynamic scheduler and runtimes)*/
    pastix_int_t   *clust2smp;       /**< clust2smp[i] = SMP node on which i_th MPI
                                          process is running, if multiple MPI processes per node            */
    pastix_int_t   *core2clust;      /**< core2clust[i] = cluster owning the core i                         */

    /**
     * @}
     * @name Parameters arrays
     * @{
     */
    pastix_int_t   *iparm;           /**< In/Out Integer parameters                                         */
    double         *dparm;           /**< In/Out Float parameters                                           */
    char          **dirtemp;         /**< Temporary unique directory to store output files
                                          (must point to pastix_data->dirtemp                               */
    /**
     * @}
     * @name Other
     * @{
     */
    //BubbleTree        *btree;        /**< arbre de bulles                                                   */
    EliminTree        *etree;        /**< the elimination tree                                              */
    CostMatrix        *costmtx;      /**< the cost bounded to each cblk and blok                            */
    Cand              *candtab;      /**< processor candidate tab                                           */
    FILE              *tracefile;    /**< File holding the simulated trace                                  */
    /**
     * @}
     */
} BlendCtrl;

int  blendCtrlInit ( pastix_data_t *pastix_data,
                     BlendCtrl     *ctrl );

void blendCtrlExit (BlendCtrl *);

void getCommunicationCosts( const BlendCtrl *ctrl,
                            pastix_int_t clustsrc,
                            pastix_int_t clustdst,
                            pastix_int_t sync_comm_nbr,
                            double *startup,
                            double *bandwidth);

#endif /* _blendctrl_h_ */

/**
 * @}
 */
