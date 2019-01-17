/**
 *
 * @file cand.c
 *
 * PaStiX analyse functions to manipulate candidates on the elimination tree
 * structure.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_cand
 * @{
 *    This module contains all subroutines to initialize the candidates array
 *    for each supernode, as well as supernode properties that are defined by
 *    level such as 2D layouts and 2D tasks.
 *
 **/
#include "common.h"
#include "symbol.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "solver.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the candtab array with default values.
 *
 *******************************************************************************
 *
 * @param[in] cblknbr
 *          The size of the candtab array.
 *
 *******************************************************************************
 *
 * @return The array of size cblknbr+1 of Cand structure initialized.
 *
 *******************************************************************************/
Cand *
candInit( pastix_int_t cblknbr )
{
    Cand *candtab, *cand;
    pastix_int_t i;

    MALLOC_INTERN( candtab, cblknbr+1, Cand );
    cand = candtab;

    for(i=-1;i<cblknbr;i++, cand++)
    {
        cand->costlevel = 0.0;
        cand->treelevel = 0;
        cand->fcandnum  = -1;
        cand->lcandnum  = -1;
        cand->fccandnum = -1;
        cand->lccandnum = -1;
        cand->cluster   = -1;
        cand->cblktype  = CBLK_LAYOUT_2D | CBLK_TASKS_2D;
    }

    return candtab+1;
}

/**
 *******************************************************************************
 *
 * @brief Exit and free the candtab structure given
 *
 *******************************************************************************
 *
 * @param[inout] candtab
 *          The candtab array to free.
 *
 *******************************************************************************/
void
candExit( Cand *candtab )
{
    candtab--;
    memFree_null( candtab );
}

/**
 *******************************************************************************
 *
 * @brief Print the candidates array into the candtab.txt file
 *
 *******************************************************************************
 *
 * @param[in] candtab
 *          The array of size cblknbr to print in the file.
 *
 * @param[in] cblknbr
 *          The size of the candtab array.
 *
 * @param[inout] directory
 *          Directory where to store the file. If NULL, initialized to
 *          pastix-XXXXXX, with XXXXXX a random generated string.
 *
 *******************************************************************************/
void
candSave( const Cand    *candtab,
          pastix_int_t   cblknbr,
          char         **directory )
{
    pastix_int_t i;
    FILE *f = NULL;

    f = pastix_fopenw( directory, "candtab.txt", "w" );

    fprintf(f, "%ld\n", (long)cblknbr );
    for(i=-1;i<cblknbr;i++)
    {
        fprintf(f, "%lf %ld %ld %ld %ld %ld %ld %ld\n",
                (double)candtab[i].costlevel,
                (long)candtab[i].treelevel,
                (long)candtab[i].fcandnum,
                (long)candtab[i].lcandnum,
                (long)candtab[i].fccandnum,
                (long)candtab[i].lccandnum,
                (long)candtab[i].cluster,
                (long)candtab[i].cblktype );
    }
    fclose(f);
}

/**
 *******************************************************************************
 *
 * @brief Set the clusters candidates from the cores canditates
 *
 *******************************************************************************
 *
 * @param[inout] candtab
 *          On entry, the array of candidates with the first and last core
 *          candidates initialized.
 *          On exit, the array of candidate with te first and last cluster
 *          candidate information updated.
 *
 * @param[in] cblknbr
 *          The size of the candtab array.
 *
 * @param[in] core2clust
 *          An array that defines the cluster (MPI process) that owns each core
 *          candidate.
 *
 * @param[in] coresnbr
 *          The size of the core2clust array.
 *
 *******************************************************************************/
void
candSetClusterCand(       Cand         *candtab,
                          pastix_int_t  cblknbr,
                    const pastix_int_t *core2clust,
                          pastix_int_t  coresnbr )
{
    pastix_int_t i;
    (void)coresnbr;

    assert( candtab[-1].fcandnum == 0 );
    assert( candtab[-1].lcandnum == coresnbr-1 );
    candtab[-1].fccandnum = core2clust[ candtab[-1].fcandnum ];
    candtab[-1].lccandnum = core2clust[ candtab[-1].lcandnum ];

    for(i=0; i<cblknbr; i++) {
        assert( candtab[i].fcandnum >= 0 );
        assert( candtab[i].lcandnum >= 0 );
        assert( candtab[i].fcandnum < coresnbr );
        assert( candtab[i].lcandnum < coresnbr );
        candtab[i].fccandnum = core2clust[ candtab[i].fcandnum ];
        candtab[i].lccandnum = core2clust[ candtab[i].lcandnum ];
    }
}

/**
 *******************************************************************************
 *
 * @brief Check the correctness of the computed candidates
 *
 * Each node of the elimination tree must have a set of candidates included in
 * its father's set.
 *
 *******************************************************************************
 *
 * @param[in] candtab
 *          On entry, the array of candidates to check.
 *
 * @param[in] symbmtx
 *          The symbol matrix structure associated to the candidate array.
 *
 *******************************************************************************
 *
 * @retval 0 if bad candidat set appear.
 * @retval 1 if success.
 *
 *******************************************************************************/
int
candCheck( const Cand            *candtab,
           const symbol_matrix_t *symbmtx )
{
    pastix_int_t i, j;
    pastix_int_t facecblknum;

    for(i=0; i<symbmtx->cblknbr; i++)
    {
        for(j = symbmtx->cblktab[i].bloknum;
            j < symbmtx->cblktab[i+1].bloknum; j++)
        {
            facecblknum = symbmtx->bloktab[j].fcblknm;

            if( (candtab[i].fcandnum < candtab[facecblknum].fcandnum) ||
                (candtab[i].lcandnum > candtab[facecblknum].lcandnum) )
            {
                errorPrint("bad processor candidat sets : cblk %ld candidat =[%ld %ld] father %ld candidat = [%ld %ld].",
                           (long)i, (long)candtab[i].fcandnum, (long)candtab[i].lcandnum,
                           (long)facecblknum, (long)candtab[facecblknum].fcandnum,
                           (long)candtab[facecblknum].lcandnum);
                return 0;
            }
        }
    }
    return 1;
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to update the cost fields of the both the candtab
 *        array, and the elimination tree structure.
 *
 *******************************************************************************
 *
 * @param[in]    rootnum
 *               Root of the subtree.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where fields treelevel and
 *               costlevel are updated. Treelevel represent the depth of the
 *               node in the elimination tree, and costlevel the cost of the
 *               path from rootnum to the root of the tree.
 *
 * @param[inout] etree
 *               Pointer to the global elimination tree structure. The node
 *               fields total and subtree are updated. Total represents the
 *               total cost of of the current node only, and subtree the cost of
 *               the current node and all its descendents.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 * @param[in]    costmtx
 *               Pointer to the cost matrix associated to the symbol matrix and
 *               that holds the cost of each cblk and blok.
 *
 * @param[out]   cripath
 *               On exit, contains the length of the critical path of the
 *               subtree.
 *
 *******************************************************************************
 *
 * @return The cost of the subtree.
 *
 *******************************************************************************/
static inline double
candSubTreeBuild( pastix_int_t           rootnum,
                  Cand                  *candtab,
                  EliminTree            *etree,
                  const symbol_matrix_t *symbmtx,
                  const CostMatrix      *costmtx,
                  double                *cripath )
{
    double cost, mycp = 0.0;
    pastix_int_t i, son;

    /* Get cost of current node */
    if ( rootnum == -1 ) {
        cost = 0.;
    }
    else {
        cost = costmtx->cblkcost[rootnum];
    }
    etree->nodetab[ rootnum ].total   = cost;
    etree->nodetab[ rootnum ].subtree = cost;

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        double soncp = 0.0;

        son = eTreeSonI(etree, rootnum, i);
        candtab[ son ].treelevel = candtab[ rootnum ].treelevel - 1;
        candtab[ son ].costlevel = candtab[ rootnum ].costlevel - cost;

        etree->nodetab[ rootnum ].subtree +=
            candSubTreeBuild( son, candtab, etree, symbmtx, costmtx, &soncp );

        mycp = (mycp > soncp) ? mycp : soncp;
    }

    /* Update local critical path */
    if (rootnum >= 0)
    {
        pastix_int_t bloknum = symbmtx->cblktab[ rootnum ].bloknum;
        pastix_int_t fcblknm;

        /* Add Facto and solve */
        mycp += costmtx->blokcost[ bloknum ];

        /* Add first GEMM */
        bloknum++;
        if (bloknum <  symbmtx->cblktab[ rootnum+1 ].bloknum) {
            fcblknm = symbmtx->bloktab[ bloknum ].fcblknm;

            while( (bloknum <  symbmtx->cblktab[ rootnum+1 ].bloknum) &&
                   (fcblknm == symbmtx->bloktab[ bloknum   ].fcblknm) )
            {
                mycp += costmtx->blokcost[ bloknum ];
                bloknum++;
            }
        }
    }
    etree->nodetab[ rootnum ].cripath = mycp;
    *cripath = mycp;

    return etree->nodetab[ rootnum ].subtree;
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to compute the distribution of the nodes among the
 * different levels.
 *
 * This function defines which cblk are candidates to be stored with a 2D
 * layout, computed as 2D tasks, compressible, and/or be part of the Schur
 * complement. The criteria to remove the 2D task flags or the LR flag is the
 * width of the cblk. As soon as one cblk that do not match the criterai is
 * found, all its descendant loose the flag too.
 *
 *******************************************************************************
 *
 * @param[in]    rootnum
 *               Root of the subtree.
 *
 * @param[in]    cblktype
 *               List of flags that can be forwarded to rootnum and its
 *               descendents.
 *
 * @param[in]    ratiolimit2D
 *               Ratio that defines the minimal size to allow the flag 2D to be
 *               forwarded to the sons.
 *
 * @param[in]    ratiolimitLR
 *               Ratio that defines the minimal size to allow the flag LR to be
 *               forwarded to the sons.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where field cblktype is
 *               updated. cblktype defines the optimization/properties that are
 *               defined on each cblk and which are defined by level in the
 *               tree.
 *
 * @param[in]    etree
 *               Pointer to the global elimination tree structure that is used
 *               to travel through the cblk, and affect the properies with the
 *               correct filiation property.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 *******************************************************************************/
static inline void
candSubTreeDistribFirstWidth( pastix_int_t           rootnum,
                              pastix_int_t           cblktype,
                              pastix_int_t           ratiolimit2D,
                              pastix_int_t           ratiolimitLR,
                              Cand                  *candtab,
                              const EliminTree      *etree,
                              const symbol_matrix_t *symbmtx )
{
    pastix_int_t i, son;
    pastix_int_t width = symbmtx->cblktab[ rootnum ].lcolnum - symbmtx->cblktab[ rootnum ].fcolnum + 1;

    if ( (cblktype & CBLK_IN_SCHUR) &&
         (symbmtx->cblktab[ rootnum ].lcolnum < symbmtx->schurfcol) )
    {
        cblktype &= ~(CBLK_IN_SCHUR);
    }

    if( (cblktype & CBLK_TASKS_2D) && (width < ratiolimit2D) ) {
        cblktype = cblktype & (~CBLK_TASKS_2D);
    }

    if( (cblktype & CBLK_COMPRESSED) && (width < ratiolimitLR) ) {
        cblktype = cblktype & (~CBLK_COMPRESSED);
    }

    candtab[ rootnum ].cblktype = cblktype;

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candSubTreeDistribFirstWidth( son, candtab[ rootnum ].cblktype,
                                      ratiolimit2D, ratiolimitLR,
                                      candtab, etree, symbmtx );
    }
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to compute the distribution of the nodes among the
 * different levels.
 *
 * This function defines which cblk are candidates to be stored with a 2D
 * layout, computed as 2D tasks, compressible, and/or be part of the Schur
 * complement. The criteria to remove the 2D task flags and the LR flag is the
 * width of the cblk.
 * This function enables the flags to the deepest nodes that respect the
 * condition, as well as to all their ascendants in the elimination tree.
 *
 *******************************************************************************
 *
 * @param[in]    rootnum
 *               Root of the subtree.
 *
 * @param[in]    cblktype
 *               List of flags that can be forwarded to rootnum and its
 *               descendents.
 *
 * @param[in]    ratiolimit2D
 *               Ratio that defines the minimal size to allow the flag 2D for a
 *               cblk and its ascendants.
 *
 * @param[in]    ratiolimitLR
 *               Ratio that defines the minimal size to allow the flag LR for a
 *               cblk and its ascendants.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where field cblktype is
 *               updated. cblktype defines the optimization/properties that are
 *               defined on each cblk and which are defined by level in the
 *               tree.
 *
 * @param[in]    etree
 *               Pointer to the global elimination tree structure that is used
 *               to travel through the cblk, and affect the properies with the
 *               correct filiation property.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 *******************************************************************************
 *
 * @return the cblktype flag of the root of the subtree.
 *
 *******************************************************************************/
static inline pastix_int_t
candSubTreeDistribDeepestWidth( pastix_int_t           rootnum,
                                pastix_int_t           cblktype,
                                pastix_int_t           ratiolimit2D,
                                pastix_int_t           ratiolimitLR,
                                Cand                  *candtab,
                                const EliminTree      *etree,
                                const symbol_matrix_t *symbmtx )
{
    pastix_int_t i, son;
    pastix_int_t sonstype = 0;
    pastix_int_t width = symbmtx->cblktab[ rootnum ].lcolnum - symbmtx->cblktab[ rootnum ].fcolnum + 1;

    if ( (cblktype & CBLK_IN_SCHUR) &&
         (symbmtx->cblktab[ rootnum ].lcolnum < symbmtx->schurfcol) )
    {
        cblktype &= ~(CBLK_IN_SCHUR);
    }

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        sonstype |= candSubTreeDistribDeepestWidth( son, cblktype,
                                                    ratiolimit2D, ratiolimitLR,
                                                    candtab, etree, symbmtx );
    }

    if( (cblktype & CBLK_TASKS_2D) && (width < ratiolimit2D) ) {
        cblktype = cblktype & (~CBLK_TASKS_2D);
    }

    if( (cblktype & CBLK_COMPRESSED) && (width < ratiolimitLR) ) {
        cblktype = cblktype & (~CBLK_COMPRESSED);
    }

    candtab[ rootnum ].cblktype = cblktype | sonstype;
    return candtab[ rootnum ].cblktype;
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to compute the distribution of the nodes among the
 * different levels based on depth.
 *
 * This function defines which cblk are candidates to be stored with a 2D layout,
 * computed as 2D tasks, and/or be part of the Schur complement. The criteria to
 * remove the 2D task flags is the depth on the elimination tree.
 *
 *******************************************************************************
 *
 * @param[in]    rootnum
 *               Root of the subtree.
 *
 * @param[in]    cblktype
 *               List of flags that can be forwarded to rootnum and its
 *               descendents.
 *
 * @param[in]    level2D
 *               The number of levels of the tree that will be flagged as 2D
 *               tasks.
 *
 * @param[in]    ratiolimitLR
 *               Ratio that defines the minimal size to allow the flag LR for a
 *               cblk and its ascendants.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where field cblktype is
 *               updated. cblktype defines the optimization/properties that are
 *               defined on each cblk and which are defined by level in the
 *               tree.
 *
 * @param[in]    etree
 *               Pointer to the global elimination tree structure that is used
 *               to travel through the cblk, and affect the properies with the
 *               correct filiation property.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 *******************************************************************************
 *
 * @return the cblktype flag of the root of the subtree.
 *
 *******************************************************************************/
static inline pastix_int_t
candSubTreeDistribDeepestLevel( pastix_int_t           rootnum,
                                pastix_int_t           cblktype,
                                pastix_int_t           level2D,
                                pastix_int_t           ratiolimitLR,
                                Cand                  *candtab,
                                const EliminTree      *etree,
                                const symbol_matrix_t *symbmtx )
{
    pastix_int_t i, son;
    pastix_int_t sonstype = 0;
    pastix_int_t width = symbmtx->cblktab[ rootnum ].lcolnum - symbmtx->cblktab[ rootnum ].fcolnum + 1;

    if ( (cblktype & CBLK_IN_SCHUR) &&
         (symbmtx->cblktab[ rootnum ].lcolnum < symbmtx->schurfcol) )
    {
        cblktype &= ~(CBLK_IN_SCHUR);
    }

    if( (cblktype & CBLK_TASKS_2D) && (level2D <= 0) ) {
        cblktype = cblktype & (~CBLK_TASKS_2D);
    }

    level2D--;
    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        sonstype |= candSubTreeDistribDeepestLevel( son, cblktype,
                                                    level2D, ratiolimitLR,
                                                    candtab, etree, symbmtx );
    }

    if( (cblktype & CBLK_COMPRESSED) && (width < ratiolimitLR) ) {
        cblktype = cblktype & (~CBLK_COMPRESSED);
    }

    candtab[ rootnum ].cblktype = cblktype | sonstype;
    return candtab[ rootnum ].cblktype;
}

/**
 *******************************************************************************
 *
 * @brief Recursive function to compute the distribution of the nodes among the
 * different levels based on depth.
 *
 * This function defines which cblk are candidates to be stored with a 2D layout,
 * computed as 2D tasks, and/or be part of the Schur complement. The criteria to
 * remove the 2D task flags is the depth on the elimination tree.
 *
 *******************************************************************************
 *
 * @param[in]    rootnum
 *               Root of the subtree.
 *
 * @param[in]    cblktype
 *               List of flags that can be forwarded to rootnum and its
 *               descendents.
 *
 * @param[in]    level2D
 *               The number of levels of the tree that will be flagged as 2D
 *               tasks.
 *
 * @param[in]    ratiolimitLR
 *               Ratio that defines the minimal size to allow the flag LR for a
 *               cblk and its ascendants.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array where field cblktype is
 *               updated. cblktype defines the optimization/properties that are
 *               defined on each cblk and which are defined by level in the
 *               tree.
 *
 * @param[in]    etree
 *               Pointer to the global elimination tree structure that is used
 *               to travel through the cblk, and affect the properies with the
 *               correct filiation property.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 *******************************************************************************/
static inline void
candSubTreeDistribFirstLevel( pastix_int_t           rootnum,
                              pastix_int_t           cblktype,
                              pastix_int_t           level2D,
                              pastix_int_t           ratiolimitLR,
                              Cand                  *candtab,
                              const EliminTree      *etree,
                              const symbol_matrix_t *symbmtx )
{
    pastix_int_t i, son;
    pastix_int_t width = symbmtx->cblktab[ rootnum ].lcolnum - symbmtx->cblktab[ rootnum ].fcolnum + 1;

    if ( (cblktype & CBLK_IN_SCHUR) &&
         (symbmtx->cblktab[ rootnum ].lcolnum < symbmtx->schurfcol) )
    {
        cblktype &= ~(CBLK_IN_SCHUR);
    }

    if( (cblktype & CBLK_TASKS_2D) && (level2D <= 0) ) {
        cblktype = cblktype & (~CBLK_TASKS_2D);
    }

    if( (cblktype & CBLK_COMPRESSED) && (width < ratiolimitLR) ) {
        cblktype = cblktype & (~CBLK_COMPRESSED);
    }

    level2D--;
    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candSubTreeDistribFirstLevel( son, cblktype,
                                      level2D, ratiolimitLR,
                                      candtab, etree, symbmtx );
    }

    candtab[ rootnum ].cblktype = cblktype;
}

/**
 *******************************************************************************
 *
 * @brief Finish to build the candtab array for the proportionnal mapping.
 *
 * This function defines which cblk are candidates to be stored with a 2D
 * layout, computed as 2D tasks, and/or be part of the Schur complement. It also
 * set the cost of each cblk and the total cost of each subtree.
 *
 *******************************************************************************
 *
 * @param[in]    level_tasks2d
 *               Defines the level in the elimination tree to switch from 1D to
 *               2D tasks in the cblk flag.
 *               If < 0, automatic level is defined based on the width of the
 *               cblk with respect to the minimal width_tasks2d width.
 *               If 0, all cblks are tagged as 1D.
 *               If > 0, only the first level_tasks2d level of the tree are
 *               tagged as 2D.
 *
 * @param[in]    width_tasks2d
 *               If level_tasks2d < 0, defines the minimal width to keep the 2D
 *               flag on a goven cblk.
 *
 * @param[in]    lr_when
 *               Define if compression technics will be used or not. If not
 *               PastixCompressNever, then all cblkw with a larger width than
 *               lr_width are tagged as compressible.
 *
 * @param[in]    lr_width
 *               Define the minimal width to flag a cblk as compressible if
 *               lr_when != PastixCompressNever.
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array that needs to be initialized.
 *
 * @param[inout] etree
 *               Pointer to the elimination tree that needs to be construct on entry.
 *               On exit, the cost of each node, and the total cost of its
 *               associated subtree is initialized.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 * @param[in]    costmtx
 *               Pointer to the cost matrix associated to the symbol matrix and
 *               that holds the cost of each cblk and blok.
 *
 *******************************************************************************/
void
candBuild( pastix_int_t level_tasks2d, pastix_int_t width_tasks2d,
           pastix_compress_when_t lr_when, pastix_int_t lr_width,
           Cand                  *candtab,
           EliminTree            *etree,
           const symbol_matrix_t *symbmtx,
           const CostMatrix      *costmtx )
{
    double cp = 0.0;
    pastix_int_t i, son, root = eTreeRoot(etree);
    pastix_int_t cblktype = CBLK_LAYOUT_2D | CBLK_TASKS_2D | CBLK_IN_SCHUR | CBLK_COMPRESSED;

#if defined(PASTIX_CUDA_FERMI)
    cblktype = CBLK_IN_SCHUR | CBLK_COMPRESSED;
#endif

    /* Let's start with the root */
    candtab[ root ].costlevel = 0.;
    candtab[ root ].treelevel = 0;

    candSubTreeBuild( root, candtab, etree, symbmtx, costmtx, &cp );

    if ( lr_when == PastixCompressNever ) {
        lr_width = PASTIX_INT_MAX;
    }

    for(i=0; i<etree->nodetab[root].sonsnbr; i++)
    {
        son = eTreeSonI(etree, root, i);
#if defined(PASTIX_BLEND_DEEPEST_DISTRIB)
        /*
         * Find the deepest node that matches the criteria for a flag, and assign
         * the flag to all its ancestors to the root
         */
        if( level_tasks2d < 0 )
        {
            candSubTreeDistribDeepestWidth( son, cblktype,
                                            width_tasks2d, lr_width,
                                            candtab, etree, symbmtx );
        }
        else
        {
            candSubTreeDistribDeepestLevel( son, cblktype,
                                            level_tasks2d, lr_width,
                                            candtab, etree, symbmtx );
        }
#else
        /*
         * Propagate the flags to all the sons as long as the node matches the
         * criteria to keep them. Stops earlier than previous case with btterfly
         * like meshes.
         */
        if( level_tasks2d < 0 )
        {
            candSubTreeDistribFirstWidth( son, cblktype,
                                          width_tasks2d, lr_width,
                                          candtab, etree, symbmtx );
        }
        else
        {
            candSubTreeDistribFirstLevel( son, cblktype,
                                          level_tasks2d, lr_width,
                                          candtab, etree, symbmtx );
        }
#endif
    }
}

/**
 *******************************************************************************
 *
 * @brief Update the candtab array costs after the symbol split algorithm has
 * been applied.
 *
 * This function update the costs and critical path of each node after the
 * symbol matrix has been split to generate more parallelism.
 *
 *******************************************************************************
 *
 * @param[inout] candtab
 *               Pointer to the global candtab array that needs to be updated.
 *
 * @param[inout] etree
 *               Pointer to the elimination tree that needs to be construct on entry.
 *               On exit, the cost of each node, and the total cost of its
 *               associated subtree is updated.
 *
 * @param[in]    symbmtx
 *               Pointer to the symbol matrix we are working with.
 *
 * @param[in]    costmtx
 *               Pointer to the cost matrix associated to the symbol matrix and
 *               that holds the cost of each cblk and blok.
 *
 *******************************************************************************/
void
candUpdate( Cand                  *candtab,
            EliminTree            *etree,
            const symbol_matrix_t *symbmtx,
            const CostMatrix      *costmtx )
{
    double cp = 0.0;
    pastix_int_t root = eTreeRoot(etree);

    /* Let's start with the root */
    candtab[ root ].costlevel = 0.;
    candtab[ root ].treelevel = 0;

    candSubTreeBuild( root, candtab, etree, symbmtx, costmtx, &cp );
}

/**
 *@}
 */

