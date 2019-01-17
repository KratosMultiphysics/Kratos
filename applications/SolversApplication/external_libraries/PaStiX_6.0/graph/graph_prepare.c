/**
 *
 * @file graph_prepare.c
 *
 * PaStiX graph construction routines
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
#include "spm.h"
#include "graph.h"
#if defined(PASTIX_DISTRIBUTED) && 0
#include "cscd_utils_intern.h"
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine removes the diagonal edges from a centralized graph.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          On entry, the pointer to the graph structure with possible diagonal
 *          edges (i,i).
 *          On exit, all entries of type (i,i) are removed from the graph.
 *
 *******************************************************************************/
void
graphNoDiag( pastix_graph_t *graph )
{
    pastix_int_t  i, j, indj;
    pastix_int_t  n   = graph->n;
    pastix_int_t *ia  = graph->colptr;
    pastix_int_t *ja  = graph->rows;
    pastix_int_t *ja2 = graph->rows;
    int baseval = ia[0];

    indj = ia[0];
    for(i=0; i<n; i++, ia++)
    {
        for (j = ia[0]; j < ia[1]; j++, ja++ )
        {
            /* If diagonal element, we skip it */
            if ( (ja[0]-baseval) == i ) {
                continue;
            }
            /* Otherwise we save it */
            else {
                *ja2 = *ja;
                ja2++;
            }
        }
        /* Save the new ia[i] */
        ia[0] = indj;

        /* Store the new ia[i+1] */
        indj = (ja2 - graph->rows) + baseval;
    }
    ia[0] = indj;

    graph->rows =
        (pastix_int_t *) memRealloc ( graph->rows,
                                      (ia[0]-baseval)*sizeof (pastix_int_t) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine sortes the subarray of edges of each vertex.
 *
 * WARNING: The sort is always performed, do not call this routine
 * when it is not required.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          On entry, the pointer to the graph structure.
 *          On exit, the same graph with subarrays of edges sorted by ascending
 *          order.
 *
 *******************************************************************************/
void
graphSort( pastix_graph_t *graph )
{
    pastix_int_t *ia = graph->colptr;
    pastix_int_t *ja = graph->rows;
    pastix_int_t  n = graph->n;
    pastix_int_t  itercol;
    int baseval = ia[0];

    /* Sort in place each subset */
    for (itercol=0;itercol<n;itercol++)
    {
        pastix_int_t frow = ia[itercol]   - baseval;
        pastix_int_t lrow = ia[itercol+1] - baseval;

        intSort1asc1( (ja+frow), (lrow-frow) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine initializes the graph.
 *
 * This routine will also symmetrize the graph, remove duplicates,
 * remove the diagonal elements, and keep only the lower part.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pointer to the solver instance. On exit, the fields n, cols,
 *          rows and loc2glob are initialized for future steps of the solver.
 *
 * @param[in] spm
 *          The initial user spm that needs to be transformed in a
 *          correct graph for future call in ordering and symbolic factorization
 *          routines.
 *
 * @param[out] graph
 *          On exit, the pointer to the allocated graph structure is returned.
 *          The graph can then be used with ordering and symbol factorization
 *          tools.
 *          The graph is symmetrized without diagonal elements and rows array is
 *          sorted.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success,
 * @retval !0 on failure.
 *
 *******************************************************************************/
int
graphPrepare(      pastix_data_t   *pastix_data,
             const spmatrix_t      *spm,
                   pastix_graph_t **graph )
{
    pastix_graph_t *tmpgraph  = NULL;
    pastix_int_t *iparm   = pastix_data->iparm;
    pastix_int_t  procnum = pastix_data->procnum;
    pastix_int_t  n       = spm->gN;
    pastix_int_t *colptr  = spm->colptr;
    pastix_int_t *rows    = spm->rowptr;
    pastix_int_t *loc2glob= spm->loc2glob;
    int io_strategy = iparm[IPARM_IO_STRATEGY];

    MALLOC_INTERN( tmpgraph, 1, pastix_graph_t );
    memset( tmpgraph, 0, sizeof(pastix_graph_t) );

    if (iparm[IPARM_VERBOSE] > PastixVerboseNo)
        pastix_print(procnum, 0, "%s", OUT_SUBSTEP_GRAPH);

    if (io_strategy & PastixIOLoadGraph)
    {
        graphLoad( pastix_data, tmpgraph );
    }
    else
    {
        /*
         * Centralized graph
         */
        if (loc2glob == NULL)
        {
            tmpgraph->gN = spm->gN;

            /*
             * TODO: change test for requirement from the user to correct his
             * mistakes
             */
            if ( (spm->mtxtype == SpmSymmetric) ||
                 (spm->mtxtype == SpmHermitian) )
            {
                if (iparm[IPARM_VERBOSE] > PastixVerboseNo)
                    pastix_print(procnum, 0, "%s", OUT_ORDER_SYMGRAPH);

                graphSymmetrize( n, colptr, rows, loc2glob, tmpgraph );
                assert( n == tmpgraph->n );
            }
            else
            {
                pastix_int_t nnz = colptr[n]-colptr[0];
                tmpgraph->n = n;
                MALLOC_INTERN(tmpgraph->colptr, (n+1), pastix_int_t);
                MALLOC_INTERN(tmpgraph->rows,   nnz,   pastix_int_t);
                memcpy(tmpgraph->colptr, colptr, (n+1)*sizeof(pastix_int_t));
                memcpy(tmpgraph->rows,   rows,     nnz*sizeof(pastix_int_t));

                if (iparm[IPARM_VERBOSE] > PastixVerboseNo) {
                    pastix_print(procnum, 0, "%s", OUT_ORDER_SORT);
                }
                graphSort( tmpgraph );
            }

            {
                if (iparm[IPARM_VERBOSE] > PastixVerboseNo) {
                    pastix_print(procnum, 0, "%s", OUT_ORDER_NODIAG);
                }

                graphNoDiag( tmpgraph );
            }
        }
#if defined(PASTIX_DISTRIBUTED) && 0
        /*
         * Distributed graph
         */
        else
        {
            MPI_Comm     pastix_comm = pastix_data->pastix_comm;
            pastix_int_t gN = 0;
            int copy_l2g = 1;

            /* Check that we use Fortran ordering */
            assert( colptr[0] == 1 );

            MPI_Allreduce(&n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
            if ( (spm->mtxtype == SpmSymmetric) ||
                 (spm->mtxtype == SpmHermitian) )
            {
                cscd_symgraph_int(n, colptr, rows, NULL,
                                  &(tmpgraph->n),
                                  &(tmpgraph->colptr),
                                  &(tmpgraph->rows), NULL,
                                  loc2glob,
                                  pastix_comm, 1);
                assert( n == tmpgraph->n );
            }
            else
            {
                pastix_int_t nnz = colptr[n]-colptr[0];
                tmpgraph->n = n;
                MALLOC_INTERN(tmpgraph->colptr, (n+1), pastix_int_t);
                MALLOC_INTERN(tmpgraph->rows,   nnz,   pastix_int_t);
                memcpy(tmpgraph->colptr, colptr, (n+1)*sizeof(pastix_int_t));
                memcpy(tmpgraph->rows,   rows,     nnz*sizeof(pastix_int_t));
            }
            MALLOC_INTERN(tmpgraph->loc2glob,   n,   pastix_int_t);
            memcpy(tmpgraph->loc2glob, loc2glob, n*sizeof(pastix_int_t));

            cscd_noDiag(tmpgraph->n,
                        tmpgraph->colptr,
                        tmpgraph->rows,
                        NULL,
                        loc2glob);

            /* Create contiguous partitions for ordering tools */
            {
                pastix_int_t i;
                int ok  = 0;
                int gok = 0;

                /* Check if matrix is allready partitionned in contiguous blocks */
                for (i = 0; i < n-1; i++)
                    if (loc2glob[i] != (loc2glob[i+1] - 1) )
                        ok = 1;

                MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_SUM, pastix_comm);

                /*
                 * If the partition is incorrect, we create a permutation to linearize the sets
                 */
                if ( !gok ) {
                    pastix_int_t  ldisp;
                    pastix_int_t *all_n;
                    pastix_int_t *displs;

                    /* Gather the locals n */
                    MALLOC_INTERN(all_n,  pastix_data->procnbr, pastix_int_t);
                    MALLOC_INTERN(displs, pastix_data->procnbr, pastix_int_t);

                    MPI_Allgather(&n,    1, PASTIX_MPI_INT,
                                  all_n, 1, PASTIX_MPI_INT,
                                  pastix_comm);

                    displs[0] = 0;
                    for (i = 1; i < pastix_data->procnbr; i++)
                        displs[i] = displs[i-1] + all_n[i-1];
                    ldisp = displs[ pastix_data->procnum ] + 1;

                    /* Collect the locals loc2glob */
                    MALLOC_INTERN(pastix_data->PTS_peritab, gN, pastix_int_t);
                    MPI_Allgatherv((void*)loc2glob, n, PASTIX_MPI_INT,
                                   pastix_data->PTS_peritab, all_n, displs, PASTIX_MPI_INT,
                                   pastix_comm);

                    memFree_null(displs);
                    memFree_null(all_n);

                    MALLOC_INTERN(pastix_data->PTS_permtab, gN, pastix_int_t);
                    for (i = 0; i < gN; i++)
                        pastix_data->PTS_permtab[pastix_data->PTS_peritab[i]-1] = i+1;

                    /* Apply the new permutation to the local graph */
                    for (i = 0; i < (tmpgraph->colptr)[n] - 1; i++)
                        tmpgraph->rows[i] = pastix_data->PTS_permtab[(tmpgraph->rows)[i]-1];

                    /* Initialize loc2glob */
                    copy_l2g = 0;
                    MALLOC_INTERN(tmpgraph->loc2glob, n, pastix_int_t);
                    for (i = 0; i < n; i++,ldisp++)
                        tmpgraph->loc2glob[i] = ldisp;
                }

                tmpgraph->gN = gN;
            }

            if (copy_l2g)
            {
                MALLOC_INTERN(tmpgraph->loc2glob, n, pastix_int_t);
                memcpy(tmpgraph->loc2glob, loc2glob, n*sizeof(pastix_int_t));
            }
        }
#else
        assert(loc2glob == NULL );
#endif
    }

    tmpgraph->dof  = spm->dof;
    tmpgraph->dofs = spm->dofs;
    *graph = tmpgraph;
    return PASTIX_SUCCESS;
}
