/**
 *
 * @file propmap.c
 *
 * PaStiX analyse proportionnal mapping functions.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_propmap
 * @{
 *
 **/
#include "common.h"
#include "symbol.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "queue.h"

/**
 * @brief Minimal work ratio to accept the contribution from one additional candidate
 */
#define CROSS_TOLERANCE 0.1

/**
 * @brief Proportional mapping structure to forward the arguments throught the
 *        recursive calls.
 */
typedef struct propmap_s {
    const EliminTree *etree;       /**< Elimination tree to map                    */
    Cand             *candtab;     /**< Candidate array for each node of the etree */
    pastix_int_t      candnbr;     /**< Number of candidates available             */
    int               nocrossproc; /**< Enable/disable the distribution of one candidate on multiple branches */
} propmap_t;


/**
 *******************************************************************************
 *
 * @brief Set the given candidates to all the subtree w/o conditions.
 *
 *******************************************************************************
 *
 * @param[in] pmptr
 *          Pointer to the parameters of the proportional mapping algorithm.
 *
 * @param[in] rootnum
 *          Index of the root of the subtree to be given to [fcandnum,lcandnum]
 *
 * @param[in] fcandnum
 *          Rank of the first candidate attributed to this subtree.
 *
 * @param[in] lcandnum
 *          Rank of the last candidate attributed to this subtree.
 *
 * @param[in] cluster
 *          The cluster value for the subtree.
 *
 *******************************************************************************/
static inline void
propMappSubtreeOn1P( const propmap_t *pmptr,
                     pastix_int_t     rootnum,
                     pastix_int_t     fcandnum,
                     pastix_int_t     lcandnum,
                     pastix_int_t     cluster )
{
    pastix_int_t i;
    pastix_int_t sonsnbr;

    pmptr->candtab[rootnum].fcandnum = fcandnum;
    pmptr->candtab[rootnum].lcandnum = lcandnum;
    pmptr->candtab[rootnum].cluster  = cluster;

    /* Recursively apply the affectation to the sons */
    sonsnbr = pmptr->etree->nodetab[rootnum].sonsnbr;
    for(i=0;i<sonsnbr;i++) {
        propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                             fcandnum, lcandnum, cluster );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Apply the proportional mapping algorithm to a subtree.
 *
 *******************************************************************************
 *
 * @param[in] pmptr
 *          Pointer to the parameters of the proportional mapping algorithm.
 *
 * @param[in] rootnum
 *          Index of the root of the subtree to process.
 *
 * @param[in] fcandnum
 *          Rank of the first candidate attributed to this subtree.
 *
 * @param[in] lcandnum
 *          Rank of the last candidate attributed to this subtree.
 *
 * @param[in] cluster
 *          The cluster value for the subtree.
 *
 * @param[inout] cost_remain
 *          Array of size (lcandnum-fcandnum+1).
 *          Stores the remaining idle time of each candidate to balance the load
 *          among them. On exit, the cost is updated with the affected subtrees.
 *
 *******************************************************************************/
static inline void
propMappSubtree( const propmap_t *pmptr,
                 pastix_int_t     rootnum,
                 pastix_int_t     fcandnum,
                 pastix_int_t     lcandnum,
                 pastix_int_t     cluster,
                 double          *cost_remain )
{
    pastix_queue_t *queue_tree;
    pastix_int_t p;
    pastix_int_t candnbr, scandnbr;
    pastix_int_t fcand = 0;
    pastix_int_t lcand = 0;
    pastix_int_t i;
    pastix_int_t sonsnbr;
    double isocost;
    double aspt_cost;
    double cumul_cost;
    double *sub_cost_remain = NULL;
    double epsilon;

    candnbr = lcandnum - fcandnum + 1;

    /* If only one candidate, all the subtree belongs to it */
    if (candnbr == 1)
    {
        propMappSubtreeOn1P( pmptr, rootnum, fcandnum, lcandnum, cluster );
        return;
    }

    /* Set the cand group for this tree node */
    pmptr->candtab[rootnum].fcandnum = fcandnum;
    pmptr->candtab[rootnum].lcandnum = lcandnum;
    pmptr->candtab[rootnum].cluster  = cluster;

    /* This treenode is a leave, return */
    if(pmptr->etree->nodetab[rootnum].sonsnbr == 0)
    {
        return;
    }

    /* Work that each processor is intended to get from this treenode */
    isocost = pmptr->etree->nodetab[rootnum].total / candnbr;
    for(p=0;p<candnbr;p++)
        cost_remain[p] -= isocost;

    /* Get cost remaining in the descendance of the treenode */
    aspt_cost = pmptr->etree->nodetab[rootnum].subtree - pmptr->etree->nodetab[rootnum].total;

    /*
     * If the first and last candidate have already received more work that they
     * should, we remove them from the set
     */
    if(cost_remain[0] <= 0)
        fcand = 1;
    else
        fcand = 0;
    if (cost_remain[candnbr-1] <= 0)
        candnbr--;

    assert(fcand < candnbr);
    assert(candnbr >= 1);

    /* Make sure that the sum of cost_remain in used proc is at least equals to aspt_cost */
    cumul_cost = 0;
    for(i=fcand; i<candnbr; i++)
        cumul_cost += cost_remain[i];

    if (cumul_cost < aspt_cost) {
        double ratio = aspt_cost / cumul_cost;
        for(i=fcand; i<candnbr; i++)
            cost_remain[i] *= ratio;
    }

    /* Compute the minimun participation rate of a candidat processor */
    epsilon = (CROSS_TOLERANCE * cumul_cost) / ((double)candnbr);

    /*
     * Compute the cand group for each proc
     */
    sonsnbr = pmptr->etree->nodetab[rootnum].sonsnbr;

    /* Create the list of sons sorted by descending order of cost */
    MALLOC_INTERN(queue_tree, 1, pastix_queue_t);
    pqueueInit(queue_tree, sonsnbr);
    for(i=0; i<sonsnbr; i++)
    {
        double soncost;

        /* Cost in the current subtree to be mapped */
        cumul_cost = -pmptr->etree->nodetab[eTreeSonI(pmptr->etree, rootnum, i)].subtree;

        /* Cost of the root node in the subtree */
        soncost    = -pmptr->etree->nodetab[eTreeSonI(pmptr->etree, rootnum, i)].total;

        pqueuePush2(queue_tree, i, cumul_cost, (pastix_int_t)soncost);
    }

    /* Proportionnal mapping of the subtree on remaining candidates           */
    /* The first stage deals only with nodes that require multiple candidates */
    while (pqueueSize(queue_tree) > 0)
    {
        i = pqueuePop2( queue_tree, &cumul_cost, NULL );
        cumul_cost = -cumul_cost;

        /*
         * If the subtree cost is less than what the first candidate can
         * process, folowing nodes will be smaller and following candidates have
         * at least the same amount of ressources available so, we skip to the
         * second stage
         */
        if (cumul_cost <= cost_remain[fcand]) {
            cost_remain[fcand] -= cumul_cost;
            propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                                 fcandnum+fcand, fcandnum+fcand, cluster);
            break;
        }

        /*
         * If a candidate cannot participate to multiple subtree, we check with
         * epsilon to avoid overflow
         */
        if ( pmptr->nocrossproc &&
             (cumul_cost <= (cost_remain[fcand]+epsilon)) )
        {
            cost_remain[fcand] -= cumul_cost;
            propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                                 fcandnum+fcand, fcandnum+fcand, cluster );
            break;
        }

        /* If the first candidate doesn't have enough ressources, we skip it */
        if( (fcand < candnbr-1) && (cost_remain[fcand] <= epsilon) )
        {
            fcand++;
        }

        /*
         * Computes how many candidate will participate to this node.  We add
         * candidate as long as we have some and they all have more than espilon
         * extra work.
         */
        lcand = fcand;
        scandnbr = 1;
        cumul_cost -= cost_remain[fcand];
        while ((( pmptr->nocrossproc && (cumul_cost > ((double)scandnbr)*epsilon)) ||
                (!pmptr->nocrossproc && (cumul_cost >                    epsilon))) &&
               (lcand < candnbr - 1))
        {
            lcand++; scandnbr++;
            assert( cost_remain[lcand] > 0 );
            cumul_cost -= cost_remain[lcand];
        }

        /*
         * Prepare the sub_cost_remain array.
         * If some cost is remaining, we distribute equally the work to each
         * candidate, otherwise we give the cost remaining to the last candidate
         * and set to zero the remaining cost of the first candidates.
         */
        MALLOC_INTERN(sub_cost_remain, lcand-fcand+1, double);
        if (cumul_cost > 0)
        {
            isocost = cumul_cost / scandnbr;
            for(p=0; p<scandnbr; p++)
            {
                sub_cost_remain[p]   = cost_remain[fcand+p] + isocost;
                cost_remain[fcand+p] = - isocost;
            }
        }
        else
        {
            for(p=0; p<scandnbr-1; p++)
            {
                sub_cost_remain[p] = cost_remain[fcand+p];
                cost_remain[fcand+p] = 0.0;
            }
            sub_cost_remain[scandnbr-1] = cost_remain[lcand] + cumul_cost;
            cost_remain[fcand+scandnbr-1] = -cumul_cost;  /* cumul_cost <= 0 */
        }
        /* Go on to subtree */
        propMappSubtree( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                         fcandnum+fcand, fcandnum+lcand, cluster, sub_cost_remain);
        memFree_null(sub_cost_remain);

        if ( (lcand < candnbr - 1) &&
             ( pmptr->nocrossproc ||
               (cost_remain[lcand] < epsilon) ) )
        {
            fcand = lcand+1;
        }
        else
        {
            if ( pmptr->nocrossproc )
                break;
            fcand = lcand;
        }
    }

    if (pqueueSize(queue_tree) > 0)
    {
        pastix_queue_t *queue_proc;

        /*
         * Second stage:
         * Distribute the single candidate subtree. At each step of the algorithm we
         * associate together the candidate with maximum cost_remain and the largest
         * son.
         */

        /* Fill queue proc order by remain cost descending */
        MALLOC_INTERN(queue_proc, 1, pastix_queue_t);
        pqueueInit(queue_proc, candnbr);
        for (i=0; i<candnbr; i++)
            pqueuePush1(queue_proc, i, -cost_remain[i]);

        while (pqueueSize(queue_tree) > 0)
        {
            /* Get the largest node */
            i = pqueuePop2( queue_tree, &cumul_cost, NULL );
            cumul_cost = -cumul_cost;

            /* Get the candidate with the largest cost_remain */
            fcand = pqueuePop(queue_proc);

            /* Map them together */
            propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                                 fcandnum+fcand, fcandnum+fcand, cluster);

            /* Update cost_remain and re-insert into the sorted queue */
            cost_remain[fcand] -= cumul_cost;
            pqueuePush1(queue_proc, fcand, -cost_remain[fcand]);
        }

        pqueueExit(queue_proc);
        memFree(queue_proc);
    }

    pqueueExit(queue_tree);
    memFree(queue_tree);
    return;
}


/**
 * @}
 */

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Apply the proportional mapping algorithm.
 *
 * This function computes the proportionnal mapping of the elimination tree. The
 * result is a set of potential candidates to compute each node of the
 * elimination tree. The real candidate will be affected during the simulation
 * with simuRun(). It is then important to reduce as much as possible the number
 * of candidates per node, while keeping enough freedom for the scheduling to
 * allow a good load balance and few idle times in the final static decision.
 *
 *******************************************************************************
 *
 * @param[inout] candtab
 *          On entry, the candtab array must conatins the cost of each node of
 *          the elimination tree, and there depth in the tree as computed by
 *          candBuild().
 *          On exit, the fields fcandnum, and lcandnum are computed with the
 *          proportional mapping algorithm that tries to balance the load
 *          between the candidate and distribute the branches to everyone
 *          according to their cost.
 *
 * @param[in] etree
 *          The elimination tree to map on the ressources.
 *
 * @param[in] candnbr
 *          The total number of candidate to distribute over the elimination
 *          tree.
 *
 * @param[in] nocrossproc
 *          If nocrossproc is enabled, candidates can NOT be part of two
 *          subranches with different co-workers in each branch.
 *          If nocrossproc is disabled, candidate can be shared between two
 *          subranches if the amount of extra work exceeds 10%.
 *
 * @param[in] allcand
 *          No proportional mapping is performed and everyone is candidate to
 *          everything. This will have a large performance impact on the
 *          simulation.
 *
 *******************************************************************************/
void
propMappTree( Cand               *candtab,
              const EliminTree   *etree,
              pastix_int_t        candnbr,
              int nocrossproc, int allcand )
{
    propmap_t pmdata;
    pastix_int_t p;

    /* Prepare the stucture */
    pmdata.candtab     = candtab;
    pmdata.etree       = etree;
    pmdata.candnbr     = candnbr;
    pmdata.nocrossproc = nocrossproc;

    if (allcand) {
        propMappSubtreeOn1P( &pmdata, eTreeRoot(etree),
                             0, candnbr-1, 0 );
    }
    else {
        double *cost_remain = NULL;
        double isocost;

        /* Prepare the initial cost_remain array */
        MALLOC_INTERN(cost_remain, candnbr, double);
        isocost = etree->nodetab[ eTreeRoot(etree) ].subtree / candnbr;

        for(p=0; p<candnbr; p++) {
            cost_remain[p] = isocost;
        }

        propMappSubtree( &pmdata, eTreeRoot(etree),
                         0, candnbr-1,
                         0, cost_remain);

        memFree_null(cost_remain);
    }
}
