/**
 *
 * @file graph_symmetrize.c
 *
 * PaStiX graph symmetrize routine
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "graph.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Symmetrize a given graph
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of vertex of the original graph.
 *
 * @param[in] ia
 *          Array of size n+1.
 *          Index of first edge for each vertex in ja array.
 *
 * @param[in] ja
 *          Array of size nnz = ia[n] - ia[0].
 *          Edges for each vertex.
 *
 * @param[in] loc2glob
 *          Array of size n.
 *          Global numbering of each local vertex.
 *
 * @param[inout] newgraph
 *          The initialized graph structure where the symmetrized graph will be
 *          stored. The allocated data must be freed with graphClean.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success,
 * @retval PASTIX_ERR_ALLOC if allocation went wrong,
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int graphSymmetrize(       pastix_int_t    n,
                     const pastix_int_t   *ia,
                     const pastix_int_t   *ja,
                     const pastix_int_t   *loc2glob,
                           pastix_graph_t *newgraph )
{
    pastix_int_t *nbrEltCol = NULL; /* nbrEltCol[i] = Number of elt to add in column i */
    pastix_int_t  itercol, iterrow, iterrow2; /* iterators */
    pastix_int_t *newia;
    pastix_int_t *newja;
    pastix_int_t  nnz;
    int baseval = ia[0];

    if (loc2glob == NULL) {
        /*
         * The graph is centralized and/or duplicated on each node.
         * Everyone work on the same set of data and compute the symmetric graph
         */

        MALLOC_INTERN(nbrEltCol, n, pastix_int_t);
        memset( nbrEltCol, 0, n*sizeof(pastix_int_t) );

        /*
         * Compute number of elements by col to add for correcting the CSC
         */
        for (itercol=0; itercol<n; itercol++)
        {
            pastix_int_t frow = ia[itercol]   - baseval;
            pastix_int_t lrow = ia[itercol+1] - baseval;
            for (iterrow=frow; iterrow<lrow; iterrow++)
            {
                pastix_int_t rowidx = ja[iterrow]  - baseval;
                if ( rowidx != itercol )
                {
                    /* It is not a diagonal element, so we have (i,j) and we look for (j,i) element */
                    /* i = itercol+1, j=ja[iterrow] */
                    pastix_int_t frow2 = ia[rowidx]   - baseval;
                    pastix_int_t lrow2 = ia[rowidx+1] - baseval;
                    int flag = 0;

                    for (iterrow2=frow2; iterrow2<lrow2; iterrow2++)
                    {
                        if (ja[iterrow2] == itercol+baseval)
                        {
                            /* We found (j,i) so let's stop this madness */
                            flag = 1;
                            break;
                        }
                    }

                    if (flag == 0)
                    {
                        /* We never found (j,i) so we increase nbrEltCol[j] */
                        (nbrEltCol[rowidx])++;
                    }
                }
            }
        }

        /* Let's compute the new ia in C numbering */
        MALLOC_INTERN(newia, n+1, pastix_int_t);

        newia[0] = ia[0];
        for (itercol=0;itercol<n;itercol++)
        {
            newia[itercol+1] = newia[itercol] + ia[itercol+1] - ia[itercol] + nbrEltCol[itercol];
        }

        assert( newia[n] >= ia[n] );
        nnz = newia[n] - baseval;

        /* Let's build the new ja */
        MALLOC_INTERN(newja, nnz, pastix_int_t);

        if ( newia[n] > ia[n])
        {
            for (itercol=0;itercol<nnz;itercol++)
            {
                newja[itercol] = -1;
            }

            for (itercol=0;itercol<n;itercol++)
            {
                pastix_int_t frow  = ia[itercol]   - baseval;
                pastix_int_t lrow  = ia[itercol+1] - baseval;
                pastix_int_t nbelt = ia[itercol+1] - ia[itercol];

                assert( newia[itercol] >= ia[itercol] );
                assert( newia[itercol] <= newia[itercol+1] );

                if ( newia[itercol] == newia[itercol+1] ) {
                    pastix_print_error( "graphSymmetrize: unknown %ld is fully disconnected (even from itself)\n",
                                        (long)itercol );
                    continue;
                }

                /* Let's copy what we already have at the end of the space reserved,
                 * we will add the new elements in front of it */
                memcpy( newja + (newia[itercol+1] - baseval) - nbelt,
                        ja    + frow,
                        nbelt * sizeof(pastix_int_t) );

                /* We add (j, i) edges missing */
                for (iterrow=frow; iterrow<lrow; iterrow++)
                {
                    pastix_int_t rowidx = ja[iterrow] - baseval;
                    if ( rowidx != itercol )
                    {
                        /* It is not a diagonal element, so we have (i,j) and we look for (j,i) element */
                        /* i = itercol+1, j=ja[iterrow] */
                        pastix_int_t frow2 = ia[rowidx]   - baseval;
                        pastix_int_t lrow2 = ia[rowidx+1] - baseval;
                        int flag = 0;

                        for (iterrow2=frow2; iterrow2<lrow2; iterrow2++)
                        {
                            if (ja[iterrow2] == itercol+baseval)
                            {
                                /* We found (j,i) so let's stop this madness */
                                flag = 1;
                                break;
                            }
                        }

                        if (flag == 0)
                        {
                            /* We never found (j,i) so we increase nbrEltCol[j] */
                            nbrEltCol[rowidx]--;
                            assert( nbrEltCol[rowidx] >= 0 );
                            newja[ newia[rowidx] - baseval + nbrEltCol[rowidx] ] = itercol + baseval;
                        }
                    }
                }
            }

            for (itercol=0;itercol<nnz;itercol++)
            {
                assert( newja[itercol] != -1);
            }

            /* Sort in place each subset */
            for (itercol=0;itercol<n;itercol++)
            {
                pastix_int_t frow = newia[itercol]   - baseval;
                pastix_int_t lrow = newia[itercol+1] - baseval;

                intSort1asc1( (newja+frow), (lrow-frow));
                assert( nbrEltCol[itercol] == 0);
            }
        }
        else
        {
            memcpy( newja, ja, nnz * sizeof(pastix_int_t) );
        }

        memFree_null( nbrEltCol );

        newgraph->gN = n;
        newgraph->n  = n;
        newgraph->colptr   = newia;
        newgraph->rows     = newja;
        newgraph->loc2glob = NULL;
    }
#if defined(PASTIX_WITH_MPI) && defined(PASTIX_DISTRIBUTED)
    else {
        /*
         * The distributed interface is used, everyone collaborate to get a
         * symmetric graph from its local subset of data.
         */

        /* Compute the global2local array if needed */
        if ( newgraph->glob2loc == NULL ) {
        }

    }
#endif /*defined(PASTIX_WITH_MPI) && defined(PASTIX_DISTRIBUTED)*/

    return EXIT_SUCCESS;
}
