/**
 *
 * @file symbol_reorder.c
 *
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Vincent Bridonneau
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "symbol.h"
#include "symbol_reorder.h"
#include "pastix/order.h"
#include "blend/queue.h"
#include "blend/extendVector.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Sequential version for reordering
 *
 *  This algorithm reorders cblks. (Sequential version).
 *
 *******************************************************************************
 *
 * @param[int, out] pastix_data
 *          The pastix_data providing the scheduler, the symbolic structure,
 *          and the ordering providing by Scotch that will be updated with
 *          the new rows permutation for each supernode. It will also gives
 *          a srop criteria the reordering of each supernode. The split_level
 *          field activates the split level heuristic, dividing distances
 *          computations into two stages: for upper and for lower
 *          contruibuting supernodes.
 *
 * @param[out] levels
 *          The pointer to the levels structure, giving the level of
 *          each supernode in the elimination tree. To be computed inside.
 *
 * @param[in] maxdepth
 *          The maximum depth in the elimination tree.
 *
 *******************************************************************************/
static inline void
sequential_reorder( pastix_data_t         *pastix_data,
                    pastix_int_t           maxdepth,
                    pastix_int_t          *levels )
{
    symbol_matrix_t *symbptr = pastix_data->symbmtx;
    symbol_cblk_t   *cblk;
    pastix_int_t     cblknbr;
    pastix_int_t    *iparm = pastix_data->iparm;
    pastix_order_t  *order = pastix_data->ordemesh;
    pastix_int_t     itercblk;
    pastix_int_t    *depthweight;

    cblk = symbptr->cblktab;
    cblknbr = symbptr->cblknbr;

    /*
     * Solves the Traveler Salesman Problem on each cblk to minimize the number
     * of off-diagonal blocks per row
     */
    MALLOC_INTERN( depthweight, maxdepth, pastix_int_t );

    for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++) {

        if ( cblk->fcolnum >= symbptr->schurfcol )
            continue;

        memset( depthweight, 0, maxdepth * sizeof(pastix_int_t) );

        symbol_reorder_cblk( symbptr, cblk, order,
                             levels,
                             depthweight, maxdepth,
                             iparm[IPARM_REORDERING_SPLIT],
                             iparm[IPARM_REORDERING_STOP] );
    }

    memFree_null( depthweight );
}

/**
 *
 * @brief The argument for reorder parallel functions data structure
 *
 * This structure provides parallel functions for reordering arguments
 * to computes the same thing sequential functions do but in parallel.
 *
 */
struct args_reorder_t
{
    pastix_data_t      *pastix_data;
    pastix_int_t        maxdepth;
    const pastix_int_t *levels;
    ExtendVectorINT    *tasktab;
};

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Parallel basic version for reordering
 *
 *  This algorithm reorders cblks. Tasks are cut such as every thread
 *  has the same number of task. (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread. This provides information
 *          about total number of thread used and the rank of the current thread.
 *
 * @param[int, out] args
 *          The argument for reordering functions.
 *
 *******************************************************************************/
static inline void
thread_preorder_basic_stategy( isched_thread_t *ctx, void *args )
{
    struct args_reorder_t *arg = (struct args_reorder_t*)args;
    pastix_data_t      *pastix_data = arg->pastix_data;
    symbol_matrix_t    *symbptr = pastix_data->symbmtx;
    symbol_cblk_t      *cblk;
    pastix_int_t       *iparm = pastix_data->iparm;
    pastix_order_t     *order = pastix_data->ordemesh;
    pastix_int_t        maxdepth = arg->maxdepth;
    pastix_int_t        ii, cblknbr, tasknbr;
    pastix_int_t       *depthweight;
    const pastix_int_t *levels = arg->levels;
    pastix_int_t        rank = (pastix_int_t)ctx->rank;
    pastix_int_t        size = (pastix_int_t)ctx->global_ctx->world_size;

    cblknbr = symbptr->cblknbr;
    tasknbr = cblknbr / size;

    if (rank < (cblknbr % size)) {
        tasknbr ++;
    }

    /*
     * Solves the Traveler Salesman Problem on each cblk to minimize the number
     * of off-diagonal blocks per row
     */
    MALLOC_INTERN( depthweight, maxdepth, pastix_int_t );

    cblk = symbptr->cblktab + rank;
    for (ii=0; ii<tasknbr; ii++, cblk += size) {

        memset( depthweight, 0, maxdepth * sizeof(pastix_int_t) );

        symbol_reorder_cblk( symbptr, cblk, order,
                             levels,
                             depthweight, maxdepth,
                             iparm[IPARM_REORDERING_SPLIT],
                             iparm[IPARM_REORDERING_STOP] );
    }

    memFree_null( depthweight );
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Parallel improved version for reordering
 *
 *  This algorithm reorders cblks after ordering
 *  them according to their size. Priority queue are used
 *  to sort cblks increasing their cost and process decreasing
 *  their load. (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread. This provides information
 *          about total number of thread used and the rank of the current thread.
 *
 * @param[int, out] args
 *          The argument for reordering functions.
 *
 *******************************************************************************/
static inline void
thread_preorder_zigzag_stategy( isched_thread_t *ctx, void *args )
{
    struct args_reorder_t *arg = (struct args_reorder_t*)args;
    pastix_data_t      *pastix_data = arg->pastix_data;
    symbol_matrix_t    *symbptr = pastix_data->symbmtx;
    symbol_cblk_t      *cblk;
    pastix_int_t       *iparm = pastix_data->iparm;
    pastix_order_t     *order = pastix_data->ordemesh;
    pastix_int_t        maxdepth = arg->maxdepth;
    pastix_int_t        ii;
    ExtendVectorINT    *tasktab;
    pastix_int_t        tasknbr;
    pastix_int_t       *depthweight;
    const pastix_int_t *levels = arg->levels;
    pastix_int_t        rank = (pastix_int_t)ctx->rank;

    /*
     * Solves the Traveler Salesman Problem on each cblk to minimize the number
     * of off-diagonal blocks per row
     */
    MALLOC_INTERN( depthweight, maxdepth, pastix_int_t );

    tasktab = arg->tasktab + rank;
    tasknbr = extendint_Size( tasktab );

    for ( ii=0; ii<tasknbr; ii++ ) {
        cblk = symbptr->cblktab + extendint_Read(tasktab, ii);

        memset( depthweight, 0, maxdepth * sizeof(pastix_int_t) );
        symbol_reorder_cblk( symbptr, cblk, order,
                             levels,
                             depthweight, maxdepth,
                             iparm[IPARM_REORDERING_SPLIT],
                             iparm[IPARM_REORDERING_STOP] );
    }

    memFree_null( depthweight );
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Function called by each thread.
 *
 *  This function select and call a strategy to reorder cblks.
 *  (Parallel version)
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread. This provides information
 *          about total number of thread used and the rank of the current thread.
 *
 * @param[int, out] args
 *          The argument for reordering functions.
 *
 *******************************************************************************/
static inline void
thread_preorder( isched_thread_t *ctx, void *args )
{
#ifdef BASIC_REORDERING_STRATEGY
    thread_preorder_basic_stategy ( ctx, args );
#else
    thread_preorder_zigzag_stategy( ctx, args );
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Computes the cost of a cblk
 *
 *  This function computes the value of a cblk using the
 *  formula : cost = a.(cblk->lcolcum - cblk->fcolnum + 1)^2
 *  where a = (cblk[1].bloknum - cblk[0].bloknum) / 2 + 1
 *  It represents the load (number of calcul) a process
 *  needs to reorder the cblk.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          The cblk we want to compute the cost.
 *
 *******************************************************************************/
static inline double
cost( symbol_cblk_t *cblk )
{
    double n = (double)(cblk->lcolnum - cblk->fcolnum + 1);
    return n*n * ((double)(cblk[1].bloknum - cblk[0].bloknum) / 2.0 + 1.0);
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Order cblks for each process
 *
 *  This function sorts cblks increasing their cost
 *  and process decreasingly their load using priority queue.
 *  After that every process has a load nearly equal, except
 *  the first one wich has a greater load than others.
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          The context of the current thread. This provides information
 *          about total number of thread used and the rank of the current thread.
 *
 * @param[int, out] args
 *          The argument for reordering functions.
 *
 *******************************************************************************/
static inline void
order_tasks( isched_t              *ctx,
             struct args_reorder_t *args )
{
    pastix_data_t   *pastix_data = args->pastix_data;
    symbol_matrix_t *symbmtx = pastix_data->symbmtx;
    pastix_queue_t   cblks;
    pastix_queue_t   procs;
    pastix_int_t     cblknbr = symbmtx->cblknbr;
    pastix_int_t     size = ctx->world_size;
    pastix_int_t     itercblk, iterproc;
    pastix_int_t     cblk_id, proc_id;
    double           cblk_cost, proc_cost;
    symbol_cblk_t   *cblk;

    pqueueInit( &cblks, cblknbr );
    pqueueInit( &procs, size );

    /*
     * Sort the cblks decreasing
     */
    cblk = symbmtx->cblktab;
    for ( itercblk=0; itercblk < cblknbr; itercblk++, cblk++ ) {

        if ( cblk->fcolnum >= symbmtx->schurfcol )
            continue;

        pqueuePush1( &cblks, itercblk, -cost(cblk) );
    }

    for ( iterproc=0; iterproc < size; ++iterproc) {
        pqueuePush1( &procs, iterproc, 0. );
    }

    while ( pqueueSize( &cblks ) > 0 ) {
        cblk_id = pqueuePop1( &cblks, &cblk_cost );
        proc_id = pqueuePop1( &procs, &proc_cost );
        proc_cost += -cblk_cost; // Negative because of reverse sort
        pqueuePush1 ( &procs, proc_id, proc_cost );

        extendint_Add( args->tasktab + proc_id, cblk_id );
    }

    pqueueExit( &cblks );
    pqueueExit( &procs );
}


/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Prepare arguments for parallel subroutines and order cblks.
 *
 *  This function creates array for each process that countains cblks
 *  the process has to reorder. These cblks are provided by order_tasks
 *  functions. Afterwards, cblks are reordered by reordering functions.
 *
 *******************************************************************************
 *
 * @param[int, out] pastix_data
 *          The pastix_data providing the scheduler, the symbolic structure,
 *          and the ordering providing by Scotch that will be updated with
 *          the new rows permutation for each supernode. It will also gives
 *          a srop criteria the reordering of each supernode. The split_level
 *          field activates the split level heuristic, dividing distances
 *          computations into two stages: for upper and for lower
 *          contruibuting supernodes.
 *
 * @param[out] levels
 *          The pointer to the levels structure, giving the level of
 *          each supernode in the elimination tree. To be computed inside.
 *
 * @param[in] maxdepth
 *          The maximum depth in the elimination tree.
 *
 *******************************************************************************/
static inline void
thread_reorder( pastix_data_t *pastix_data,
                pastix_int_t   maxdepth,
                pastix_int_t  *levels )
{
    struct args_reorder_t args_reorder = { pastix_data, maxdepth, levels, NULL };
    pastix_int_t          size = pastix_data->isched->world_size;
    pastix_int_t          iterproc;
    pastix_int_t          cblknbr = pastix_data->symbmtx->cblknbr;
    pastix_int_t          cblkavg = pastix_imax(1, (pastix_int_t)(cblknbr / size));

    MALLOC_INTERN( args_reorder.tasktab, size, ExtendVectorINT );

    for ( iterproc=0; iterproc < size; ++iterproc ) {
        extendint_Init( args_reorder.tasktab + iterproc, cblkavg );
    }
    order_tasks( pastix_data->isched, &args_reorder );

    isched_parallel_call( pastix_data->isched, thread_preorder, &args_reorder );

    for ( iterproc=0; iterproc < size; ++iterproc ) {
        extendint_Exit( args_reorder.tasktab + iterproc );
    }

    memFree_null( args_reorder.tasktab );
}

static void (*reorder_table[4])(pastix_data_t *, pastix_int_t , pastix_int_t *) = {
    sequential_reorder,
    thread_reorder,
    NULL,
    NULL
};

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Reorder all node
 *
 * This function computes the set each cblks, and then call a TSP heuristic
 * to minimize the Hamiltonian Path.
 *
 *******************************************************************************
 *
 * @param[int, out] pastix_data
 *          The pastix_data providing the scheduler, the symbolic structure,
 *          and the ordering providing by Scotch that will be updated with
 *          the new rows permutation for each supernode. It will also gives
 *          a srop criteria the reordering of each supernode. The split_level
 *          field activates the split level heuristic, dividing distances
 *          computations into two stages: for upper and for lower
 *          contruibuting supernodes.
 *
 * @param[out] levels
 *          The pointer to the levels structure, giving the level of
 *          each supernode in the elimination tree. To be computed inside.
 *
 * @param[in] maxdepth
 *          The maximum depth in the elimination tree.
 *
 *******************************************************************************/
void
symbol_reorder( pastix_data_t *pastix_data,
                pastix_int_t   maxdepth,
                pastix_int_t  *levels )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*reorder)(pastix_data_t *, pastix_int_t , pastix_int_t * ) = reorder_table[ sched ];

    if (reorder == NULL) {
        reorder = thread_reorder;
    }
    reorder( pastix_data, maxdepth, levels );
}
