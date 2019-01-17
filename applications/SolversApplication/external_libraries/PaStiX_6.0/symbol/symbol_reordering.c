/**
 *
 * @file symbol_reordering.c
 *
 * PaStiX symbol structure reordering routines
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include "symbol.h"
#include "symbol_reorder.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Compute the level of supernode cblknum, with Scotch treetab.
 *
 *******************************************************************************
 *
 * @param[in] treetab
 *          The pointer to the elimination tree from Scotch.
 *
 * @param[in] levels
 *          The supernode array which contains levels. Used to check if
 *          the level was already computed.
 *
 * @param[in] cblknum
 *          The supernode for which the level is computed.
 *
 *******************************************************************************
 *
 * @return the level of cblknum.
 *
 *******************************************************************************/
static inline pastix_int_t
compute_cblklevel( const pastix_int_t *treetab,
                   const pastix_int_t *levels,
                   pastix_int_t        cblknum )
{
    /* If cblknum level has already been computed */
    if ( levels[cblknum] != 0 ) {
        return levels[cblknum];
    }
    else {
        pastix_int_t father = treetab[cblknum];

        if ( father == -1 ) {
            return 1;
        }
        else {
            return compute_cblklevel( treetab, levels, father ) + 1;
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Compute the distance between two rows of a same supernode.
 *
 *******************************************************************************
 *
 * @param[in] vectors
 *          The pointer to the sets of contributing supernodes for
 *          each row of the current supernode.
 *
 * @param[in] vectors_size
 *          The pointer to the sizes of each set of contributing
 *          supernode, to stop the computation when a row have been totally
 *          covered.
 *
 * @param[in] xi
 *          The index of the first row.
 *
 * @param[in] xj
 *          The index of the second row.
 *
 * @param[in] stop
 *          The stop criteria to disregard rows that are far away.
 *
 *******************************************************************************
 *
 * @return The distance between rows xi and xj.
 *
 *******************************************************************************/
static inline pastix_int_t
hamming_distance( pastix_int_t **vectors,
                  pastix_int_t  *vectors_size,
                  pastix_int_t   xi,
                  pastix_int_t   xj,
                  pastix_int_t   stop )
{
    /* For the fictive vertex */
    if ( xi == -1 ) {
        return vectors_size[xj];
    }
    if ( xj == -1 ) {
        return vectors_size[xi];
    }

    pastix_int_t sum   = 0;
    pastix_int_t *set1 = vectors[xi];
    pastix_int_t *set2 = vectors[xj];
    pastix_int_t *end1 = vectors[xi] + vectors_size[xi];
    pastix_int_t *end2 = vectors[xj] + vectors_size[xj];

    if ( vectors_size[xi] - vectors_size[xj] >= stop ) {
        return stop;
    }
    if ( vectors_size[xj] - vectors_size[xi] >= stop ) {
        return stop;
    }

    while( ( set1 < end1 ) && ( set2 < end2 ) ) {
        if( *set1 == *set2 ) {
            set1++;
            set2++;
        }
        else if( *set1 < *set2 ) {
            while ( ( set1 < end1 ) && ( *set1 < *set2 ) ) {
                sum ++;
                set1++;
            }
        }
        else if( *set1 > *set2 ) {
            while ( ( set2 < end2 ) && ( *set1 > *set2 ) ) {
                sum ++;
                set2++;
            }
        }
        else {
            errorPrint("reordering: fatal error occured");
        }

        /* The computation is stopped if sum overlapped a given limit (stop criteria) */
        if ( sum >= stop ) {
            return stop;
        }
    }

    sum += end1 - set1;
    sum += end2 - set2;

    if ( sum >= stop ) {
        return stop;
    }

    return sum;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Reorder rows of a supernode with the nearest insertion TSP heuristic.
 *
 * See reordering paper: http://epubs.siam.org/doi/10.1137/16M1062454.
 *
 *******************************************************************************
 *
 * @param[in] size
 *          Number of rows in the current supernode.
 *
 * @param[in, out] order
 *          The pointer to the ordering structure. At exit, this
 *          ordering is updated with the new ordering for the current supernode
 *          being reordered.
 *
 * @param[in] sn_id
 *          Identifier for the current supernode.
 *
 * @param[in] lw_vectors
 *          The pointer to the sets of lower contributing
 *          supernodes for each row of the current supernode. Those lower
 *          contributing supernodes correspond to supernodes with a level higher
 *          than split_level criteria.
 *
 * @param[in] lw_vectors_size
 *          The pointer to the sizes of each set of lower contributing
 *          supernode, to stop the computation when a row have been totally
 *          covered.
 *
 * @param[in] up_vectors
 *          The pointer to the sets of upper contributing
 *          supernodes for each row of the current supernode. Those upper
 *          contributing supernodes correspond to supernodes with a level smaller
 *          than split_level criteria.
 *
 * @param[in] up_vectors_size
 *          The pointer to the sizes of each set of upper contributing
 *          supernode, to stop the computation when a row have been totally
 *          covered.
 *
 * @param[in] stop_criteria
 *          The stop criteria to disregard rows that are far away.
 *
 *******************************************************************************/
static inline void
symbol_reorder_tsp( pastix_int_t size, pastix_order_t *order, pastix_int_t sn_id,
                    pastix_int_t **lw_vectors, pastix_int_t *lw_vectors_size,
                    pastix_int_t **up_vectors, pastix_int_t *up_vectors_size,
                    pastix_int_t stop_criteria )
{
    pastix_int_t  i, j, k, l, elected;
    pastix_int_t *tmpinvp;
    pastix_int_t *tmplen;
    pastix_int_t  distance;


    if ( size < 3 ) {
        return;
    }

    MALLOC_INTERN( tmpinvp, size+1, pastix_int_t );
    MALLOC_INTERN( tmplen,  size+1, pastix_int_t );
    memset( tmplen, 0, ( size + 1 ) * sizeof(pastix_int_t) );

    /* Insert a ghost element with no connexion to any supernodes */
    tmpinvp[0] = -1;
    tmpinvp[1] = 0;

    distance = hamming_distance( lw_vectors, lw_vectors_size, 0, -1, stop_criteria );

    tmplen[0] = distance;
    tmplen[1] = distance;

    for (i=1; i<size; i++) {
        pastix_int_t first_pos;
        pastix_int_t last_pos;

        pastix_int_t lw_before_pos;
        pastix_int_t lw_after_pos;

        pastix_int_t up_before_pos;
        pastix_int_t up_after_pos;

        pastix_int_t minl;
        pastix_int_t mpos;
        pastix_int_t min_cut;

        /* Start by adding the row in first position */
        lw_before_pos = hamming_distance( lw_vectors, lw_vectors_size, i,
                                          tmpinvp[0], stop_criteria );
        lw_after_pos  = hamming_distance( lw_vectors, lw_vectors_size, i,
                                          tmpinvp[1], stop_criteria );
        up_after_pos  = hamming_distance( up_vectors, up_vectors_size, i,
                                          tmpinvp[1], 1 );

        minl    = lw_before_pos + lw_after_pos - tmplen[0];
        mpos    = 1;
        min_cut = -1;

        for (j=1; j<i; j++) {
            up_before_pos = up_after_pos;
            up_after_pos  = hamming_distance( up_vectors, up_vectors_size, i,
                                              tmpinvp[j+1], 1 );

            if ( up_before_pos < 1 ||
                 up_after_pos  < 1 ) {

                /* If split was used previously, this first distance may not be already computed */
                if ( lw_after_pos == -1 )
                    lw_before_pos = hamming_distance( lw_vectors, lw_vectors_size, i,
                                                      tmpinvp[j], stop_criteria );
                else
                    lw_before_pos = lw_after_pos;


                lw_after_pos = hamming_distance( lw_vectors, lw_vectors_size, i,
                                                 tmpinvp[j+1], stop_criteria );

                l = lw_before_pos + lw_after_pos - tmplen[j];


                /* Minimize the cut between two lines, for the same TSP result */
                if ( l == minl ) {
                    if ( lw_before_pos < min_cut ) {
                        min_cut = lw_before_pos;
                        minl    = l;
                        mpos    = j + 1;
                    }
                    if ( lw_after_pos < min_cut ) {
                        min_cut = lw_after_pos;
                        minl    = l;
                        mpos    = j + 1;
                    }
                }

                /* Position that minimizes TSP */
                if ( l < minl ) {
                    minl    = l;
                    mpos   =  j + 1;
                    min_cut = lw_before_pos;
                    if ( lw_after_pos < min_cut ) {
                        min_cut = lw_after_pos;
                    }
                }

                if ( l < minl ) {
                    minl    = l;
                    mpos    = j + 1;
                    min_cut = lw_before_pos;
                    if ( lw_after_pos < min_cut ) {
                        min_cut = lw_after_pos;
                    }
                }


                /* Stop if two lines are equal (already done tmpinvp[j]) */
                if ( lw_after_pos == 0 ) {
                    min_cut = 0;
                    minl    = l;
                    mpos    = j + 1;
                    j       = i;
                }
            }
            else {
                lw_after_pos = -1;
            }

        }

        /* Test between last and first */
        first_pos = hamming_distance( lw_vectors, lw_vectors_size, i,
                                      tmpinvp[0], stop_criteria );
        last_pos  = hamming_distance( lw_vectors, lw_vectors_size, i,
                                      tmpinvp[i], stop_criteria );

        lw_before_pos = hamming_distance( lw_vectors, lw_vectors_size, i,
                                          tmpinvp[mpos-1], stop_criteria );
        lw_after_pos  = hamming_distance( lw_vectors, lw_vectors_size, i,
                                          tmpinvp[mpos  ], stop_criteria);

        l = first_pos + last_pos - tmplen[i];
        if ( l < minl ) {
            minl = l;
            mpos = i + 1;
        }

        if ( mpos > 0 ) {
            tmplen[mpos-1] = lw_before_pos;
        }

        if ( mpos < ( i + 1 ) ) {
            pastix_int_t tmpi, tmpl;
            k = i;
            l = lw_after_pos;

            /* Insert the line in the tmpinvp/tmplen arrays */
            for (j=mpos; j<i+2; j++) {
                tmpi = tmpinvp[j];
                tmpl = tmplen[j];

                tmpinvp[j] = k;
                tmplen[j]  = l;

                k = tmpi;
                l = tmpl;
            }
        }
        else {
            tmpinvp[i+1] = i;
            tmplen[i+1]  = first_pos;
        }
    }

    /* Look for the ghost element */
    elected = 0;
    for (i=0; i<size; i++) {
        if ( tmpinvp[i] == -1 ) {
            elected = i;
        }
    }

    /* Apply the local permutation to the global one */
    {
        pastix_int_t *sn_connected;
        pastix_int_t *peritab = order->peritab + order->rangtab[sn_id];

        MALLOC_INTERN( sn_connected, size, pastix_int_t );
        for (i=0; i<size; i++) {
            sn_connected[i] = peritab[ tmpinvp[(i + 1 + elected)%(size+1)] ];
        }
        memcpy( peritab, sn_connected, size * sizeof(pastix_int_t) );
        memFree_null( sn_connected );
    }

    memFree_null( tmpinvp );
    memFree_null( tmplen );
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_reordering
 *
 * @brief Reorder a supernode
 *
 * This function computes the set of contributing supernodes for each row, and
 * then call a TSP heuristic to minimize the Hamiltonian Path.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The pointer to the symbolic structure.
 *
 * @param[in] cblk
 *          The pointer to the current supernode being reordered.
 *
 * @param[in, out] order
 *          The ordering providing by Scotch. This ordering will be
 *          updated with the new rows permutation for the current supernode.
 *
 * @param[out] levels
 *          The pointer to the levels structure, giving the level of
 *          each supernode in the elimination tree. To be computed inside.
 *
 * @param[in, out] depthweight
 *          This array provides the number of supernodes
 *          corresponding to depth from 1 to depthmax.
 *
 * @param[in] depthmax
 *          The maximum depth in the elimination tree.
 *
 * @param[in] split_level
 *          Parameter to activate the split level heuristic,
 *          dividing distances computations into two stages: for upper and for
 *          lower contruibuting supernodes. If a resulting distance for upper
 *          supernodes is large enough, the computation is stopped, as long as
 *          it will be large taking into account lower contributing supernodes.
 *
 * @param[in] stop_criteria
 *          The stop criteria to disregard rows that are far away.
 *
 *******************************************************************************/
void
symbol_reorder_cblk( const symbol_matrix_t *symbptr,
                     const symbol_cblk_t   *cblk,
                     pastix_order_t     *order,
                     const pastix_int_t *levels,
                     pastix_int_t       *depthweight,
                     pastix_int_t        depthmax,
                     pastix_int_t        split_level,
                     pastix_int_t        stop_criteria )
{
    symbol_blok_t *blok;
    pastix_int_t **up_vectors, *up_vectors_size;
    pastix_int_t **lw_vectors, *lw_vectors_size;

    pastix_int_t size              = cblk->lcolnum - cblk->fcolnum + 1;
    pastix_int_t local_split_level = split_level;
    pastix_int_t i, iterblok;
    pastix_int_t *brow = symbptr->browtab;

    /**
     * Compute hamming vectors in two subsets:
     *   - The upper subset contains the cblk with level higher than the split_level
     *     in the elimination tree, (or depth lower than levels[cblk])
     *   - The lower subset contains the cblk with level lower than the split_level
     *     in the elimination tree, (or depth higher than levels[cblk])
     *
     * The delimitation between the lower and upper levels is made such that
     * the upper level represents 17% to 25% of the total number of cblk.
     */
    {
        pastix_int_t blokweight;
        pastix_int_t weight = 0;

        /* Compute the weigth of each level */
        for (iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++) {
            blok       = symbptr->bloktab + brow[iterblok];
            blokweight = blok->lrownum - blok->frownum + 1;

            depthweight[ levels[ blok->lcblknm ] - 1 ] += blokweight;
            weight += blokweight;
        }

        /**
         * Compute the split_level:
         *    We start with the given split_level parameter
         *    and we try to correct it to minimize the following iterative process
         */
        {
            /* Current for each line within the current cblk the number of contributions */
            pastix_int_t up_total = 0;
            pastix_int_t lw_total = 0;
            pastix_int_t sign     = 0;

          split:
            up_total = 0;
            lw_total = 0;

            for (i=0; i<local_split_level; i++) {
                up_total += depthweight[i];
            }
            for (; i<depthmax; i++) {
                lw_total += depthweight[i];
            }

            /* If there are too many upper bloks */
            if ( (lw_total < (5 * up_total)) &&
                 (lw_total > 10) && (up_total > 10) && (sign <= 0)) {
                local_split_level--;
                sign--;
                goto split;
            }

            /* If there are too many lower bloks */
            if ( (lw_total > (3 * up_total)) &&
                 (lw_total > 10) && (up_total > 10) && (sign >= 0) ) {
                local_split_level++;
                sign++;
                goto split;
            }
        }

        /* Compute the Hamming vector size for each row of the cblk */
        MALLOC_INTERN( up_vectors_size, size, pastix_int_t );
        memset( up_vectors_size, 0, size * sizeof(pastix_int_t) );
        MALLOC_INTERN( lw_vectors_size, size, pastix_int_t );
        memset( lw_vectors_size, 0, size * sizeof(pastix_int_t) );

        for (iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++) {
            blok = symbptr->bloktab + brow[iterblok];

            /* For upper levels in nested dissection */
            if (levels[blok->lcblknm] <= local_split_level) {
                for (i=blok->frownum; i<=blok->lrownum; i++) {
                    pastix_int_t index = i - cblk->fcolnum;
                    up_vectors_size[index]++;
                }
            }
            else {
                for (i=blok->frownum; i<=blok->lrownum; i++) {
                    pastix_int_t index = i - cblk->fcolnum;
                    lw_vectors_size[index]++;
                }
            }
        }

        /* Initiate Hamming vectors structure */
        MALLOC_INTERN(lw_vectors, size, pastix_int_t*);
        MALLOC_INTERN(up_vectors, size, pastix_int_t*);
        for (i=0; i<size; i++) {
            MALLOC_INTERN(lw_vectors[i], lw_vectors_size[i], pastix_int_t);
            MALLOC_INTERN(up_vectors[i], up_vectors_size[i], pastix_int_t);
            memset(lw_vectors[i], 0, lw_vectors_size[i] * sizeof(pastix_int_t));
            memset(up_vectors[i], 0, up_vectors_size[i] * sizeof(pastix_int_t));
        }
        memset( lw_vectors_size, 0, size * sizeof(pastix_int_t) );
        memset( up_vectors_size, 0, size * sizeof(pastix_int_t) );

        /* Fill-in vectors structure with contributing cblks */
        for (iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++)
        {
            blok = symbptr->bloktab + brow[iterblok];

            /* For upper levels in nested dissection */
            if (levels[blok->lcblknm] <= local_split_level) {
                for (i=blok->frownum; i<=blok->lrownum; i++) {
                    pastix_int_t index = i - cblk->fcolnum;
                    up_vectors[index][up_vectors_size[index]] = blok->lcblknm;
                    up_vectors_size[index]++;
                }
            }
            else{
                for (i=blok->frownum; i<=blok->lrownum; i++) {
                    pastix_int_t index = i - cblk->fcolnum;
                    lw_vectors[index][lw_vectors_size[index]] = blok->lcblknm;
                    lw_vectors_size[index]++;
                }
            }
        }
    }

    /* Apply the pseudo-TSP algorithm to the rows in the current supernode */
    symbol_reorder_tsp( size, order, cblk - symbptr->cblktab,
                        lw_vectors, lw_vectors_size,
                        up_vectors, up_vectors_size,
                        stop_criteria );

    for (i=0; i<size; i++) {
        memFree_null( lw_vectors[i] );
        memFree_null( up_vectors[i] );
    }
    memFree_null( lw_vectors );
    memFree_null( up_vectors );
    memFree_null( lw_vectors_size );
    memFree_null( up_vectors_size );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Compute the reordering on the complete matrix.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure that describes the solver instance.  On
 *          exit, the field symbmtx is updated with the new symbol matrix, and
 *          the field ordemesh is updated with the new ordering. The split_level
 *          field activates the split level heuristic,
 *          dividing distances computations into two stages: for upper and for
 *          lower contruibuting supernodes. If a resulting distance for upper
 *          supernodes is large enough, the computation is stopped, as long as
 *          it will be large taking into account lower contributing supernodes.
 *          The stop_criteria field disregards rows that are far away.
 *
 *******************************************************************************/
void
pastixSymbolReordering( pastix_data_t *pastix_data )
{
    symbol_matrix_t *symbptr = pastix_data->symbmtx;
    pastix_int_t cblknbr = symbptr->cblknbr;

    pastix_int_t i, maxdepth;
    pastix_int_t *levels;
    pastix_order_t *order = pastix_data->ordemesh;

    /* Create the levels array to compute the depth of each cblk and the maximum depth */
    {
        maxdepth = 0;
        levels   = calloc( cblknbr, sizeof(pastix_int_t) );

        for (i=0; i<cblknbr; i++) {
            levels[i] = compute_cblklevel( order->treetab, levels, i );
            maxdepth = pastix_imax( maxdepth, levels[i] );
        }
    }

    /**
     * Compute the reordering using either sequential or parallel method
     */
    symbol_reorder( pastix_data, maxdepth, levels );

    /* Update the permutation */
    for (i=0; i<symbptr->nodenbr; i++) {
        order->permtab[ order->peritab[i] ] = i;
    }
    memFree_null( levels );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * @brief Compute the number of operations required to compute the reordering on
 * the complete matrix.
 *
 * The number of operation is compuyted and then printed on the standard output.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The pointer to the symbolic structure.
 *
 *******************************************************************************/
void
pastixSymbolReorderingPrintComplexity( const symbol_matrix_t *symbptr )
{
    symbol_cblk_t  *cblk;
    pastix_int_t itercblk, iterblok;
    pastix_int_t cblknbr;
    pastix_int_t nbiops, schur_nbiops = 0;

    cblk    = symbptr->cblktab;
    cblknbr = symbptr->cblknbr;
    nbiops = 0;

    /*
     * nbcblk is the number of non zeroes intersection between indivudal rows
     * and block columns.
     */
    for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++) {
        pastix_int_t width;
        pastix_int_t nbcblk = 0;

        if (cblk->fcolnum >= symbptr->schurfcol )
            continue;

        for (iterblok=cblk[0].brownum; iterblok<cblk[1].brownum; iterblok++) {
            symbol_blok_t *blok = symbptr->bloktab + symbptr->browtab[iterblok];
            assert( blok->fcblknm == itercblk );

            nbcblk += blok->lrownum - blok->frownum + 1;
        }
        width = cblk->lcolnum - cblk->fcolnum + 1;
        nbiops += nbcblk * (width-1);

        if ( itercblk == (cblknbr-1) ) {
            schur_nbiops = nbcblk * (width-1);
        }
    }
    fprintf( stdout, OUT_REORDERING_OPS,
             (long)schur_nbiops, (double)schur_nbiops / (double)(nbiops) * 100.,
             (long)nbiops );
}
