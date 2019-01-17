/**
 *
 * @file graph.h
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
 *
 * @addtogroup pastix_graph
 * @{
 *   @brief Functions to generate and manipulate the graph structure.
 *
 *   This module provides the set of function to prepare the graph structure
 *   associated to a given sparse matrix.
 *   It is possible to symmetrize a graph, to extract a subgraph and to apply a
 *   new permutation.
 *
 **/
#ifndef _graph_h_
#define _graph_h_

/**
 * @brief Graph structure.
 *
 * This structure describes the adjacency graph of a sparse matrix.
 */
struct pastix_graph_s {
    pastix_int_t  gN;       /**< Global number of vertices in compressed graph    */
    pastix_int_t  n;        /**< Number of local vertices in compressed graph     */
    pastix_int_t *colptr;   /**< List of indirections to rows for each vertex     */
    pastix_int_t *rows;     /**< List of edges for each vertex                    */
    pastix_int_t *loc2glob; /**< Corresponding numbering from local to global     */
    pastix_int_t *glob2loc; /**< Corresponding numbering from global to local     */
    pastix_int_t  dof;      /**< Degre of freedom to move to uncompressed graph   */
    pastix_int_t *dofs;     /**< Array of the first column of each element in the
                                 expanded matrix [+1,based]                       */
};

/**
 * @name Graph basic subroutines
 * @{
 */
int  graphPrepare(       pastix_data_t   *pastix_data,
                   const spmatrix_t      *spm,
                         pastix_graph_t  **graph );
void graphBase   (       pastix_graph_t  *graph, int baseval );
void graphExit   (       pastix_graph_t  *graph );

/**
 * @}
 * @name Graph IO subroutines
 * @{
 */
void graphLoad( const pastix_data_t  *pastix_data,
                pastix_graph_t       *graph );
void graphSave( pastix_data_t        *pastix_data,
                const pastix_graph_t *graph );

/**
 * @}
 * @name Graph manipulation subroutines
 * @{
 */
void graphSort      ( pastix_graph_t *graph );
void graphNoDiag    ( pastix_graph_t *graph );
int  graphSymmetrize(       pastix_int_t    n,
                      const pastix_int_t   *ia,
                      const pastix_int_t   *ja,
                      const pastix_int_t   *loc2glob,
                            pastix_graph_t *newgraph );

int  graphIsolate   (       pastix_int_t    n,
                      const pastix_int_t   *colptr,
                      const pastix_int_t   *rows,
                            pastix_int_t    isolate_n,
                            pastix_int_t   *isolate_list,
                            pastix_int_t   **new_colptr,
                            pastix_int_t   **new_rows,
                            pastix_int_t   **new_perm,
                            pastix_int_t   **new_invp );

int  graphApplyPerm ( const pastix_graph_t *graphA,
                      const pastix_int_t   *perm,
                            pastix_graph_t *graphPA );

/**
 * @}
 */

#endif /* _graph_h_ */

/**
 * @}
 */
