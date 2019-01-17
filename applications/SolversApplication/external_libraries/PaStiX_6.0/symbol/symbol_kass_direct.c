/**
 *
 * @file symbol_kass_direct.c
 *
 * This file contains routines to create the graph of the direct
 * factorization of a graph A.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "symbol_kass.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_kass
 *
 * @brief Compute the non zero pattern of the direct factorization of a matrix
 * A, given the supernode partition associated.
 *
 *******************************************************************************
 *
 * @param[in] graphA
 *          The graph structure of the original matrix A to factorize.
 *
 * @param[in] cblknbr
 *          The number of supernodes.
 *
 * @param[in] rangtab
 *          Integer array of size cblknbr.
 *          The supernode partition for the entry graph.
 *
 * @param[out] treetab
 *          Integer array of size cblknbr.
 *          On exit, the array contains the tree associated to the supernode
 *          partition.
 *
 * @param[out] graphL
 *          The graph the structure of the non zero pattern of the factorized
 *          matrix.  On entry, a pointer to a graph structure. No need for
 *          initialization.  On exit, the structure contains the computed graph.
 *
 *******************************************************************************
 *
 * @retval >=0, the number of non zero entries in the generated graph.
 * @retval -i, if the i^th parameter is incorrect
 *
 *******************************************************************************/
pastix_int_t
kassFactDirect( const kass_csr_t   *graphA,
                      pastix_int_t  cblknbr,
                const pastix_int_t *rangtab,
                      pastix_int_t *treetab,
                      kass_csr_t   *graphL )
{
    pastix_int_t i, j, k, nnz;
    pastix_int_t nnznbr, father;
    pastix_int_t *tmpj      = NULL;
    pastix_int_t *tmp       = NULL;
    pastix_int_t *tmp2      = NULL;
    pastix_int_t *ja        = NULL;
    pastix_int_t *node2cblk = NULL;

    /* Check parameters */
    if( graphA == NULL ) {
        return -1;
    }
    if( (cblknbr < 0) || (cblknbr > graphA->n) ) {
        return -2;
    }
    if( rangtab == NULL ) {
        return -3;
    }
    if( graphL == NULL ) {
        return -5;
    }

    /* Quick return */
    if( graphA->n == 0 )
        return 0;

    MALLOC_INTERN(tmpj,      graphA->n, pastix_int_t);
    MALLOC_INTERN(tmp,       graphA->n, pastix_int_t);
    MALLOC_INTERN(node2cblk, graphA->n, pastix_int_t);

    for(k=0;k<cblknbr;k++)
        for(i=rangtab[k];i<rangtab[k+1];i++)
            node2cblk[i] = k;

    /* Compute the nnz structure of each supernode in A */
    kass_csrInit( cblknbr, graphL );
    for(k=0; k<cblknbr; k++)
    {
        /* Put the diagonal elements (In case A does not contains them) */
        j = 0;
        for(i=rangtab[k]; i<rangtab[k+1]; i++)
            tmpj[j++] = i;
        nnznbr = j;

        for(i=rangtab[k]; i<rangtab[k+1]; i++)
        {
            j = 0;

            /* We take only the elements greater than i */
            while( (j < graphA->nnz[i]) && (graphA->rows[i][j] <= i) )
                j++;

            /* Merge the actual list with the edges of the ith vertex */
            nnznbr = pastix_intset_union( nnznbr,           tmpj,
                                          graphA->nnz[i]-j, graphA->rows[i] + j,
                                          tmp);

            /* Swap tmpj and the merged set tmp */
            tmp2 = tmpj;
            tmpj = tmp;
            tmp  = tmp2;
        }

#if !defined(NDEBUG) && defined(PASTIX_DEBUG_SYMBOL)
        /* Check that the first elements are the diagonal ones */
        {
            pastix_int_t ind;
            ind = 0;
            assert(nnznbr >= (rangtab[k+1]-rangtab[k]));
            for(j=rangtab[k]; j<rangtab[k+1]; j++) {
                assert(tmpj[ind] == j);
                ind++;
            }
            assert(nnznbr > 0);
        }
#endif

        /* Update graphL */
        graphL->nnz[k] = nnznbr;
        MALLOC_INTERN(graphL->rows[k], nnznbr, pastix_int_t);
        memcpy(graphL->rows[k], tmpj, sizeof(pastix_int_t)*nnznbr);
    }

    /* Compute the symbolic factorization */
    for(k=0;k<cblknbr;k++)
    {
        i = 0;
        ja = graphL->rows[k];
        while( (i < graphL->nnz[k]) && (node2cblk[ja[i]] <= k) )
            i++;

        if( i < graphL->nnz[k] )
            father = node2cblk[ja[i]];
        else
            father = -1;
        treetab[k] = father;

        /* Merge son's nodes into father's list */
        if( (father != k) && (father > 0) )
        {
            nnznbr = pastix_intset_union( graphL->nnz[k] - i,  graphL->rows[k] + i,
                                          graphL->nnz[father], graphL->rows[father],
                                          tmpj );

            memFree( graphL->rows[father] );
            MALLOC_INTERN( graphL->rows[father], nnznbr, pastix_int_t);
            memcpy( graphL->rows[father], tmpj, sizeof(pastix_int_t)*nnznbr );
            graphL->nnz[father] = nnznbr;
        }
    }

#if defined(PASTIX_DEBUG_SYMBOL)
    /* Check that all terms of A are in the pattern */
    {
        pastix_int_t ind;
        for(k=0;k<cblknbr;k++)
        {
            /* Put the diagonal elements (A does not contains them) */
            for(i=rangtab[k];i<rangtab[k+1];i++)
            {
                j = 0;
                while( (j < graphA->nnz[i]) && (graphA->rows[i][j] < i) )
                    j++;

                for(ind = j; ind < graphA->nnz[i]; ind++)
                    assert(graphA->rows[i][ind] >= i);
                for(ind = j+1; ind < graphA->nnz[i]; ind++)
                    assert(graphA->rows[i][ind] > graphA->rows[i][ind-1]);

                ind = pastix_intset_union( graphL->nnz[k],   graphL->rows[k],
                                           graphA->nnz[i]-j, graphA->rows[i]+j,
                                           tmp );

                assert(ind <= graphL->nnz[k]);
            }
        }
    }
#endif

    memFree(node2cblk);
    memFree(tmpj);
    memFree(tmp);

    /*
     * Computes the nnz in graphL
     */
    nnz = 0;
    for(i=0; i<cblknbr; i++)
    {
        pastix_int_t ncol, nrow;
        ncol = rangtab[i+1]-rangtab[i];
        nrow = graphL->nnz[i];

        assert( nrow >= ncol );
        assert( nrow <= graphA->n );

        nnz += (ncol*(ncol+1))/2;
        nnz += (ncol*(nrow-ncol));
    }
    return nnz;
}
