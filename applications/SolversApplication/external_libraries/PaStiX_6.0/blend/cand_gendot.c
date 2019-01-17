/**
 *
 * @file cand_gendot.c
 *
 *  PaStiX analyse routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains functions to generate dot files of the elimination tree and the
 * compressed elimination tree.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_elim
 * @{
 *
 **/
#include "common.h"
#include "symbol.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"

/**
 *******************************************************************************
 *
 * @brief Compute the number of nodes in the compressed tree. The compression is
 * based on identical candidates for computations.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 * @param[in] candtab
 *          The candidate array associated to the elimination tree for
 *          additional information.
 *
 * @param[in] rootnum
 *          The root index of the subtree to compress.
 *
 *******************************************************************************
 *
 * @return The number of nodes of the subtree.
 *
 *******************************************************************************/
static inline pastix_int_t
compress_getNodeNbr( const EliminTree *etree,
                     const Cand       *candtab,
                     pastix_int_t      rootnum )
{
    pastix_int_t i, merge, fcand, lcand;
    pastix_int_t sonsnbr, nbnodes;

    fcand = candtab[rootnum].fcandnum;
    lcand = candtab[rootnum].lcandnum;

    sonsnbr = etree->nodetab[rootnum].sonsnbr;
    merge   = 1;
    nbnodes = 1;
    for( i=0; i<sonsnbr; i++ )
    {
        pastix_int_t son = eTreeSonI(etree, rootnum, i);
        nbnodes += compress_getNodeNbr( etree, candtab, son );

        if ( (fcand != candtab[son].fcandnum) ||
             (lcand != candtab[son].lcandnum) )
        {
            merge = 0;
        }
    }

    /* If all sons have the same candidate set, then they will be merged within the current node */
    if (merge) {
        nbnodes -= sonsnbr;
    }
    return nbnodes;
}

/**
 *******************************************************************************
 *
 * @brief Compress a subtree. The compression is based on identical candidates
 * for computations.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 * @param[in] candtab
 *          The candidate array associated to the elimination tree for
 *          additional information. (optionnal)
 *
 * @param[in] rootnum
 *          The root index of the subtree to compress.
 *
 * @param[inout] ctree
 *          The pointer to the compressed elimination tree.
 *
 * @param[inout] ccand
 *          The pointer to the compressed candidate array associated to ctree.
 *
 * @param[in] fathnum
 *          The index of the father in the compressed elimination tree.
 *
 * @param[inout] cnodeidx
 *          The index of the next available spot in the compressed elimation
 *          tree array to store the nodes.
 *
 * @param[inout] tmp
 *          A temporary array of the size of the number of nodes in etree for
 *          the initial root. It is used to store the elements visited in the
 *          tree.
 *
 *******************************************************************************/
static inline void
compress_setSonsNbr( const EliminTree *etree,
                     const Cand       *candtab,
                     pastix_int_t      rootnum,
                     EliminTree       *ctree,
                     Cand             *ccand,
                     pastix_int_t      fathnum,
                     pastix_int_t     *cnodeidx,
                     pastix_int_t     *tmp )
{
    eTreeNode_t *cnode;
    double total;
    pastix_int_t i, merge, fcand, lcand;
    pastix_int_t sonsnbr, gdsonsnbr;
    (*cnodeidx)++;
    assert( *cnodeidx < ctree->nodenbr );

    fcand = candtab[rootnum].fcandnum;
    lcand = candtab[rootnum].lcandnum;

    /* Copy current node */
    cnode = ctree->nodetab + (*cnodeidx);
    cnode->total   = etree->nodetab[rootnum].total;
    cnode->subtree = etree->nodetab[rootnum].subtree;
    cnode->cripath = etree->nodetab[rootnum].cripath;
    cnode->fathnum = fathnum;

    ccand[ *cnodeidx ].fcandnum = fcand;
    ccand[ *cnodeidx ].lcandnum = lcand;

    sonsnbr = etree->nodetab[rootnum].sonsnbr;
    if (sonsnbr == 0) {
        return;
    }

    memcpy( tmp, etree->sonstab + etree->nodetab[rootnum].fsonnum,
            sonsnbr * sizeof(pastix_int_t) );
    do {
        total     = 0.;
        gdsonsnbr = 0;
        merge     = 1;
        for( i=0; i<sonsnbr; i++ )
        {
            pastix_int_t son = tmp[i];

            /* Backup grandsons after the sons */
            memcpy( tmp + sonsnbr + gdsonsnbr,
                    etree->sonstab + etree->nodetab[ son ].fsonnum,
                    etree->nodetab[ son ].sonsnbr * sizeof(pastix_int_t) );
            gdsonsnbr += etree->nodetab[son].sonsnbr;
            total     += etree->nodetab[son].total;

            if ( (fcand != candtab[son].fcandnum) ||
                 (lcand != candtab[son].lcandnum) )
            {
                merge = 0;
            }
        }

        /* If all sons have the same candidate set, then they will be merged within the current node */
        if ( merge ) {
            /* Grandsons become sons */
            for (i=0; i<gdsonsnbr; i++) {
                tmp[i] = tmp[sonsnbr+i];
            }
            sonsnbr = gdsonsnbr;
            cnode->total += total;
        }
    }
    while( merge && (sonsnbr>0) );

    /* Recurse on sons */
    merge = *cnodeidx;
    for(i=0; i<sonsnbr; i++) {
        pastix_int_t son = (*cnodeidx) + 1;
        compress_setSonsNbr( etree, candtab, tmp[i],
                             ctree, ccand,
                             merge, cnodeidx, tmp+sonsnbr );
        tmp[i] = son;
    }
    cnode->sonsnbr = sonsnbr;

    /* Compress the single candidate nodes */
    for(i=0; i<sonsnbr; i++) {
        pastix_int_t j, soni, sonj;

        soni = tmp[i];
        if ( ccand[soni].fcandnum != ccand[soni].lcandnum )
            continue;

        for( j=i+1; j<sonsnbr; j++ ) {
            sonj = tmp[j];

            if ( (ccand[sonj].fcandnum != ccand[soni].fcandnum) ||
                 (ccand[sonj].lcandnum != ccand[soni].lcandnum) )
                continue;

            ctree->nodetab[soni].total   += ctree->nodetab[sonj].total;
            ctree->nodetab[soni].subtree += ctree->nodetab[sonj].subtree;
            assert( ctree->nodetab[sonj].sonsnbr == 0 );
            assert( ctree->nodetab[sonj].fathnum == ctree->nodetab[soni].fathnum );
            ctree->nodetab[sonj].fathnum = -2;

            sonsnbr--;
            tmp[j] = tmp[sonsnbr]; j--;
        }
    }
    cnode->sonsnbr = sonsnbr;

    return;
}

/**
 *******************************************************************************
 *
 * @brief Print the elimination tree in a dot file.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 * @param[in] candtab
 *          The candidate array associated to the elimination tree for
 *          additional information. (optionnal)
 *
 * @param[inout] stream
 *          The file to which write the elimination tree in the dot format.
 *
 *******************************************************************************/
void
candGenDot( const EliminTree *etree,
            const Cand       *candtab,
            FILE             *stream )
{
    pastix_int_t i;

    fprintf(stream,
            "digraph G {\n"
            "\tcolor=white\n"
            "\trankdir=BT;\n");

    for (i=0; i < etree->nodenbr; i++)
    {
        if ((etree->nodetab[i]).fathnum == -2)
            continue;

        if ( candtab == NULL ) {
            fprintf( stream, "\t\"%ld\" [label=\"#%ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\"]\n",
                     (long)i, (long)i,
                     etree->nodetab[i].subtree,
                     etree->nodetab[i].total,
                     etree->nodetab[i].cripath );
        }
        else {
            if ( candtab[i].lcandnum != candtab[i].fcandnum ) {
                fprintf( stream, "\t\"%ld\" [label=\"#%ld\\nCand: %ld - %ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\"]\n",
                         (long)i, (long)i,
                         (long)(candtab[i].fcandnum),
                         (long)(candtab[i].lcandnum),
                         etree->nodetab[i].subtree,
                         etree->nodetab[i].total,
                         etree->nodetab[i].cripath );
            }
            else {
                fprintf(stream, "\t\"%ld\" [label=\"#%ld\\nCand: %ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\" colorscheme=set312 style=filled fillcolor=%ld]\n",
                        (long)i, (long)i,
                        (long)(candtab[i].fcandnum),
                        etree->nodetab[i].subtree,
                        etree->nodetab[i].total,
                        etree->nodetab[i].cripath,
                        (long)((candtab[i].lcandnum % 12) + 1));
            }
        }
        if ((etree->nodetab[i]).fathnum == -1)
            continue;
        fprintf(stream, "\t\"%ld\"->\"%ld\"\n", (long)i, (long)((etree->nodetab[i]).fathnum));
    }
    fprintf(stream, "}\n");
}

/**
 *******************************************************************************
 *
 * @brief Print one level of the elimination subtree in a dot file.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 * @param[in] candtab
 *          The candidate array associated to the elimination tree for
 *          additional information. (optionnal)
 *
 * @param[inout] stream
 *          The file to which write the elimination tree in the dot format.
 *
 * @param[in] nblevel
 *          The number of level remaining to be printed.
 *
 * @param[in] rootnum
 *          The root of the subtree to print.
 *
 *******************************************************************************/
static inline void
candGenDotLevelSub( const EliminTree *etree,
                    const Cand       *candtab,
                    FILE             *stream,
                    pastix_int_t      nblevel,
                    pastix_int_t      rootnum )
{
    pastix_int_t i, son;

    assert( (etree->nodetab[rootnum]).fathnum != -2 );

    /* Print current node informations */
    if ( candtab == NULL ) {
        fprintf( stream, "\t\"%ld\" [label=\"#%ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\"]\n",
                 (long)rootnum, (long)rootnum,
                 etree->nodetab[rootnum].subtree,
                 etree->nodetab[rootnum].total,
                 etree->nodetab[rootnum].cripath );
    }
    else {
        if ( candtab[rootnum].lcandnum != candtab[rootnum].fcandnum ) {
            fprintf( stream, "\t\"%ld\" [label=\"#%ld\\nCand: %ld - %ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\"]\n",
                     (long)rootnum, (long)rootnum,
                     (long)(candtab[rootnum].fcandnum),
                     (long)(candtab[rootnum].lcandnum),
                     etree->nodetab[rootnum].subtree,
                     etree->nodetab[rootnum].total,
                     etree->nodetab[rootnum].cripath );
        }
        else {
            fprintf(stream, "\t\"%ld\" [label=\"#%ld\\nCand: %ld\\nSubtree cost: %e\\nNode cost: %e\\nNode CP: %e\" colorscheme=set312 style=filled fillcolor=%ld]\n",
                    (long)rootnum, (long)rootnum,
                    (long)(candtab[rootnum].fcandnum),
                    etree->nodetab[rootnum].subtree,
                    etree->nodetab[rootnum].total,
                    etree->nodetab[rootnum].cripath,
                    (long)((candtab[rootnum].lcandnum % 12) + 1));
        }
    }

    if ( nblevel > 0 ) {
        nblevel--;

        for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
        {
            son = eTreeSonI(etree, rootnum, i);
            candGenDotLevelSub( etree, candtab, stream, nblevel, son );

            fprintf(stream, "\t\"%ld\"->\"%ld\"\n", (long)son, (long)rootnum);
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Print the first levels of the elimination tree in a dot file.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 * @param[in] candtab
 *          The candidate array associated to the elimination tree for
 *          additional information. (optionnal)
 *
 * @param[inout] stream
 *          The file to which write the elimination tree in the dot format.
 *
 * @param[in] nblevel
 *          The number of levels of the elimination tree to print.
 *
 *******************************************************************************/
void
candGenDotLevel( const EliminTree *etree,
                 const Cand       *candtab,
                 FILE             *stream,
                 pastix_int_t      nblevel )
{
    fprintf(stream,
            "digraph G {\n"
            "\tcolor=white\n"
            "\trankdir=BT;\n");

    candGenDotLevelSub( etree, candtab, stream,
                        nblevel, eTreeRoot( etree ) );

    fprintf(stream, "}\n");
}

/**
 *******************************************************************************
 *
 * @brief Print the compressed elimination tree in a dot file, where all nodes
 * with the same candidates are merged together.
 *
 *******************************************************************************
 *
 * @param[in] etree
 *          The pointer to the elimination tree.
 *
 * @param[in] candtab
 *          The candidate array associated to the elimination tree for
 *          additional information.
 *
 * @param[inout] stream
 *          The file to which write the elimination tree in the dot format.
 *
 *******************************************************************************/
void
candGenCompressedDot( const EliminTree *etree,
                      const Cand       *candtab,
                      FILE             *stream )
{
    EliminTree  *ctree;
    Cand        *ccand;
    pastix_int_t cnodesnbr, cnodeidx, *tmp;

    cnodesnbr = compress_getNodeNbr( etree, candtab, eTreeRoot(etree)  );

    /* Let's create a second compressed elimination tree, and the associated candtab */
    ctree = eTreeInit( cnodesnbr );

    ccand = candInit( cnodesnbr );

    MALLOC_INTERN(tmp, etree->nodenbr, pastix_int_t);
    cnodeidx = -1;
    compress_setSonsNbr( etree, candtab, eTreeRoot( etree ),
                         ctree, ccand, -1, &cnodeidx, tmp );
    memFree_null(tmp);

    /* Write the dot file */
    candGenDot( ctree, ccand, stream );

    /* Free temporary ressources */
    candExit( ccand );
    eTreeExit( ctree );
}

/**
 * @}
 */
