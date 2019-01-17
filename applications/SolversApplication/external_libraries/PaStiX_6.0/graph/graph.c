/**
 *
 * @file graph.c
 *
 * PaStiX graph structure routines
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
 * @addtogroup pastix_graph
 * @{
 *
 **/
#include "common.h"
#include "graph.h"

/**
 *******************************************************************************
 *
 * @brief Free the content of the graph structure.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The pointer graph structure to free.
 *
 *******************************************************************************/
void graphExit( pastix_graph_t *graph )
{
    /* Parameter checks */
    if ( graph == NULL ) {
        errorPrint("graphClean: graph pointer is NULL");
        return;
    }

    graph->gN = 0;
    graph->n  = 0;

    if ( (graph->colptr == NULL) ||
         (graph->rows   == NULL) )
    {
        errorPrint("graphClean: graph pointer is not correctly initialized");
        return;
    }

    memFree_null(graph->colptr);
    memFree_null(graph->rows);

    if (graph->loc2glob != NULL)
    {
        memFree_null( graph->loc2glob );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Rebase the graph to the given value.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The graph to rebase.
 *
 * @param[in] baseval
 *          The base value to use in the graph (0 or 1).
 *
 *******************************************************************************/
void graphBase( pastix_graph_t *graph,
                int             baseval )
{
    pastix_int_t baseadj;
    pastix_int_t i, n, nnz;

    /* Parameter checks */
    if ( graph == NULL ) {
        errorPrint("graphBase: graph pointer is NULL");
        return;
    }
    if ( (graph->colptr == NULL) ||
         (graph->rows   == NULL) )
    {
        errorPrint("graphBase: graph pointer is not correctly initialized");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        errorPrint("graphBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - graph->colptr[0];
    if (baseadj == 0)
	return;

    n   = graph->n;
    nnz = graph->colptr[n] - graph->colptr[0];

    for (i = 0; i <= n; i++) {
        graph->colptr[i]   += baseadj;
    }
    for (i = 0; i < nnz; i++) {
        graph->rows[i] += baseadj;
    }

    if (graph->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            graph->loc2glob[i] += baseadj;
        }
    }
    return;
}

/**
 * @}
 */
