/**
 *
 * @file order_apply_level_order.c
 *
 * PaStiX order function that apply reverse level ordering to the elimination
 * tree.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include <string.h>
#include "common.h"
#include "elimintree.h"
#include "pastix/order.h"

/**
 *******************************************************************************
 *
 * @ingroup order_dev
 *
 * @brief This routine build the elimination tree associated to an ordering.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          The ordering structure.
 *
 *******************************************************************************
 *
 * @return the elimination tree structure associated to the ordering.
 *
 *******************************************************************************/
static inline EliminTree *
pastixOrderBuildEtree( const pastix_order_t *order )
{
    EliminTree *etree = NULL;
    eTreeNode_t *enode;
    pastix_int_t i, fathnum;

    etree = eTreeInit( order->cblknbr );

    /* Compute the fathers and the number of sons */
    enode = etree->nodetab;
    for(i=0; i<order->cblknbr; i++, enode++)
    {
        fathnum = order->treetab[i];
        enode->fathnum = fathnum;

        assert(fathnum < (order->cblknbr+1) );
        etree->nodetab[ fathnum ].sonsnbr ++;
    }

    /* Set the index of the first sons */
    eTreeSetSons( etree );

    return etree;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine reorder the elimination tree nodes per level.
 *
 *******************************************************************************
 *
 * @param[in] order
 *          The ordering structure to reorder.
 *
 * @param[in] level_tasks2d
 *          Define the ways 2D tasks are decided. If < 0, autolvel will be made
 *          based on all blocks above the minimal width_tasks2d criteria. If 0,
 *          1D tasks will be used, and if > 0, only the first level_tasks2d lvel
 *          of the elimination tree will be considered as 2D tasks.
 *
 * @param[in] width_tasks2d
 *          Define the minimal width for the supernodes that are considered as 2D
 *          blocks if level_tasks2d < 0. Unused otherwise.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if the ordering structure is incorrect.
 *
 *******************************************************************************/
int
pastixOrderApplyLevelOrder( pastix_order_t *order,
                            pastix_int_t    level_tasks2d,
                            pastix_int_t    width_tasks2d )
{
    pastix_order_t  oldorder;
    EliminTree     *etree;
    pastix_int_t    baseval;                  /* Node base value            */
    pastix_int_t    i, s, node, sonsnbr;
    pastix_int_t    nfcol, ofcol, size;

    /* Parameter checks */
    if ( order == NULL ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order pointer");
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( (order->permtab == NULL) && (order->vertnbr > 0) ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order->permtab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( order->rangtab == NULL ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order->rangtab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( (order->treetab == NULL) && (order->cblknbr > 0) ) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid order->treetab pointer");
        return PASTIX_ERR_BADPARAMETER;
    }

    if (order->cblknbr < 0) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid nunber of column blocks");
        return PASTIX_ERR_BADPARAMETER;
    }
    baseval = order->baseval;
    if (baseval < 0) {
        errorPrint ("pastixOrderApplyLevelOrder: invalid vertex node base number");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Quick return */
    if (order->cblknbr == 0) {
        return PASTIX_SUCCESS;
    }

    assert(baseval == order->rangtab[0]);

    memcpy( &oldorder, order, sizeof(pastix_order_t) );
    pastixOrderAlloc( order,
                      oldorder.vertnbr,
                      oldorder.cblknbr );

    /*
     * Build the elimination tree from top to bottom
     */
    etree = pastixOrderBuildEtree( &oldorder );

    /*
     * Build the sorted array per level
     * If autolevel is enabled for 2D, we need to sort the 2D cblks first and
     * then the 1D.
     */
    if ( level_tasks2d < 0 )
    {
        pastix_int_t  pos_1D, pos_2D;
        pastix_int_t  tot_nb_2D = 0;
        pastix_int_t  fathnum;
        pastix_int_t *sorted = order->permtab;
        int8_t *is_2D;

        MALLOC_INTERN(is_2D, order->cblknbr, int8_t);
        memset(is_2D, 0, order->cblknbr * sizeof(int8_t));

#if defined(PASTIX_BLEND_DEEPEST_DISTRIB)
        for(node=0; node<order->cblknbr; node++ ){
            size = oldorder.rangtab[ node+1 ] - oldorder.rangtab[ node ];
            if (is_2D[node] == 1) {
                continue;
            }

            fathnum = etree->nodetab[node].fathnum;

            /* Mark as 2D only if the father is a real supernode, otherwise we don't reorder */
            if ( size >= width_tasks2d )
            {
                /* Force all brothers to be considered as 2D */
                sonsnbr = etree->nodetab[fathnum].sonsnbr;
                for(i=0; i<sonsnbr; i++) {
                    s = eTreeSonI(etree, fathnum, i);
                    if (is_2D[s] == 0) {
                        is_2D[s] = 1;
                        tot_nb_2D++;
                    }
                }

                /* Force parent and thus uncles to be 2D too */
                while( (fathnum != -1) && (is_2D[fathnum] == 0) ) {
                    fathnum = etree->nodetab[fathnum].fathnum;

                    sonsnbr = etree->nodetab[fathnum].sonsnbr;
                    for(i=0; i<sonsnbr; i++) {
                        s = eTreeSonI(etree, fathnum, i);
                        if (is_2D[s] == 0) {
                            is_2D[s] = 1;
                            tot_nb_2D++;
                        }
                    }
                }
            }
        }

#else /* defined(PASTIX_BLEND_DEEPEST_DISTRIB) */

        /* First pass to choose which nodes are 2D from the top to bottom */
        for(node=order->cblknbr-1; node>-1; node--) {

            fathnum = etree->nodetab[node].fathnum;
            size = oldorder.rangtab[ node+1 ] - oldorder.rangtab[ node ];

            /* Mark as 2D only if the father is a real supernode, otherwise we don't reorder */
            if ( (size >= width_tasks2d) &&
                 ((fathnum == -1) ||
                  ((is_2D[fathnum] == 1) && (etree->nodetab[fathnum].sonsnbr == 2))) )
            {
                is_2D[node] = 1;
                tot_nb_2D++;
            }
        }
#endif

        /* Lets start by inserting the roots */
        pos_2D = 0;
        pos_1D = tot_nb_2D;

        sonsnbr = etree->nodetab[-1].sonsnbr;
        assert( etree->nodetab[-1].fsonnum == 0 );
        for (i=0; i<sonsnbr; i++) {
            node = etree->sonstab[i];
            if ( is_2D[node] ) {
                sorted[pos_2D] = node;
                pos_2D++;
            }
            else {
                sorted[pos_1D] = node;
                pos_1D++;
            }
        }

        /* Second pass to sort nodes: firstly by type (1D/2D) and then by levels */
        for(i=0; i<order->cblknbr; i++) {
            pastix_int_t current_1D = 0;
            pastix_int_t current_2D = 0;
            pastix_int_t sons1D, sons2D;

            node    = sorted[i];
            sonsnbr = etree->nodetab[node].sonsnbr;
            sons2D  = 0;

            /* Count the number of 2D sons */
            for(s=0; s<sonsnbr; s++) {
                pastix_int_t son = eTreeSonI(etree, node, s);
                if (is_2D[son]) {
                    sons2D++;
                }
            }
            sons1D = sonsnbr - sons2D;

            /*
             * We put the sons in reverse order to keep the original order
             * betwen the brothers. This matters for the Minimum Degree part of
             * the ordering algorithm.
             */
            for(s=0; s<sonsnbr; s++) {
                pastix_int_t son = eTreeSonI(etree, node, s);
                if (is_2D[son]) {
                    sorted[pos_2D + sons2D - current_2D - 1] = son;
                    current_2D++;
                }
                else{
                    sorted[pos_1D + sons1D - current_1D - 1] = son;
                    current_1D++;
                }
                etree->nodetab[ son ].fathnum = order->cblknbr - i - 1;
            }
            pos_2D += sons2D;
            pos_1D += sons1D;
        }
        memFree_null(is_2D);
    }
    else
    {
        pastix_int_t  pos;
        pastix_int_t *sorted = order->permtab;

        /* Start with the roots */
        pos = etree->nodetab[-1].sonsnbr;
        memcpy( sorted, etree->sonstab, pos * sizeof(pastix_int_t) );

        for(i=0; i<order->cblknbr; i++) {
            node = sorted[i];
            sonsnbr = etree->nodetab[node].sonsnbr;
            /*
             * We put the sons in reverse order to keep the original order
             * betwen the brothers. This matters for the Minimum Degree part of
             * the ordering algorithm.
             */
            for(s=0; s<sonsnbr; s++) {
                pastix_int_t son = eTreeSonI(etree, node, s);
                sorted[ pos + sonsnbr-1-s ] = son;
                etree->nodetab[ son ].fathnum = order->cblknbr - i - 1;
            }
            pos += sonsnbr;
        }
        assert(pos == order->cblknbr);
    }

    /* Let's rebuild peritab, treetab, and rangtab */
    order->rangtab[0] = 0;

    for(i=0; i<order->cblknbr; i++ ) {
        node = order->permtab[order->cblknbr - i - 1];
        size = oldorder.rangtab[ node+1 ] - oldorder.rangtab[ node ];

        ofcol = oldorder.rangtab[ node ];
        nfcol = order->rangtab[ i ];
        order->rangtab[ i+1 ] = nfcol + size;
        order->treetab[ i ] = etree->nodetab[ node ].fathnum;

        memcpy( order->peritab + nfcol, oldorder.peritab + ofcol,
                size * sizeof( pastix_int_t ) );
    }

    /* Update the permutation */
    for (i=0; i<order->vertnbr; i++) {
        order->permtab[ order->peritab[i] ] = i;
    }

    pastixOrderExit( &oldorder );

    eTreeExit( etree );

    return PASTIX_SUCCESS;
}
