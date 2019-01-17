/**
 *
 * @file graph_io.c
 *
 * PaStiX graph IO routines
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
#include "spm.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Load a graph from a file
 *
 * The file is named 'graphname' in the current directory.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pointer to the solver instance to get options as rank,
 *          communicators, ...
 *
 * @param[inout] graph
 *          The graph structure to store the loaded graph.
 *          The graph is read from the file named by the environment variable
 *          PASTIX_FILE_GRAPH, and if PASTIX_FILE_GRAPH is not defined, the
 *          default filename "graphname" in the current directory is used.
 *
 *******************************************************************************/
void
graphLoad( const pastix_data_t *pastix_data,
           pastix_graph_t      *graph )
{
    spmatrix_t spm;
    FILE *stream = NULL;
    char *filename = NULL;
    int env = 1;

    assert( pastix_data->procnbr == 1 );

    /* Parameter checks */
    if ( graph == NULL ) {
        return;
    }

    /*
     * Get the environment variable as second option
     */
    filename = pastix_getenv( "PASTIX_FILE_GRAPH" );
    env = 1;

    /*
     * Get the default name as third option
     */
    if ( filename == NULL ) {
        filename = "graphname";
        env = 0;
    }

    stream = pastix_fopen( filename );
    if ( stream ) {
        spmLoad( &spm, stream );
        fclose(stream);

        spmConvert( SpmCSC, &spm );
        graph->gN       = spm.gN;
        graph->n        = spm.n;
        graph->dof      = spm.dof;
        assert( spm.dof == 1 );
        graph->colptr   = spm.colptr;
        graph->rows     = spm.rowptr;
        graph->loc2glob = spm.loc2glob;
    }

    if (env) {
        pastix_cleanenv( filename );
    }
    (void)pastix_data;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Save a graph to file.
 *
 * The file is named 'graphgen' in the current directory.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pointer to the solver instance to get options as rank,
 *          communicators, ...
 *
 * @param[in] graph
 *          The graph structure to store the loaded graph.
 *          The graph is written to the file named by the environment variable
 *          PASTIX_FILE_GRAPH, and if PASTIX_FILE_GRAPH is not defined, the
 *          default filename "graphname" in the current directory is used.
 *
 *******************************************************************************/
void
graphSave( pastix_data_t        *pastix_data,
           const pastix_graph_t *graph )
{
    spmatrix_t spm;
    FILE *stream   = NULL;
    int   env      = 1;
    char *filename = NULL;

    /* Parameter checks */
    if ( graph == NULL ) {
        return;
    }

    /*
     * Get the environment variable as first option
     */
    filename = pastix_getenv( "PASTIX_FILE_GRAPH" );

    /*
     * Get the default name as second option
     */
    if ( filename == NULL ) {
        filename = "graphgen";
        env = 0;
    }

    assert( pastix_data->procnbr == 1 );
    spm.n   = graph->n;
    spm.nnz = graph->colptr[ graph->n ] - graph->colptr[ 0 ];
    spm.dof = graph->dof;
    assert( spm.dof == 1 );
    spm.colptr   = graph->colptr;
    spm.rowptr   = graph->rows;
    spm.loc2glob = graph->loc2glob;
    spm.dofs     = NULL;

    spmUpdateComputedFields( &spm );

    stream = pastix_fopenw( &(pastix_data->dirtemp), filename, "w" );
    if ( stream ) {
        spmSave( &spm, stream );
        fclose(stream);
    }

    if (env) {
        pastix_cleanenv( filename );
    }
}
