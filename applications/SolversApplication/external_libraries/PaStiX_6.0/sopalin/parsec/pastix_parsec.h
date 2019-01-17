/**
 *
 * @file pastix_parsec.h
 *
 * PaRSEC support for the numerical factorization and solve of PaStiX.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2018-07-16
 *
 * @addtogroup pastix_parsec
 * @{
 *   This module describes the functionnality provided by the runtime system
 *   PaRSEC for the numerical factorization and solve.
 *
 **/
#ifndef _pastix_parsec_h_
#define _pastix_parsec_h_

#include <parsec.h>
#include <parsec/data_distribution.h>

/**
 * @name PaRSEC sparse matrix descriptor
 * @{
 *
 *   @brief PaRSEC descriptor stucture for the sparse matrix.
 */
typedef struct parsec_sparse_matrix_desc_s {
    parsec_data_collection_t  super;      /**< Every PaRSEC descriptors must inherit from parsec_desc_t                        */
    int                       typesze;    /**< Arithmetic size                                                                 */
    int                       mtxtype;    /**< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.            */
    SolverMatrix             *solvmtx;    /**< Solver matrix structure that describes the problem and stores the original data */
    void                    **d_blocktab; /**< Pointer to GPU arrays that contains frownum,lrownum of each block for Fermi     */
} parsec_sparse_matrix_desc_t;

void parsec_sparse_matrix_init( SolverMatrix *solvmtx,
                                int typesize, int mtxtype,
                                int nodes, int myrank );
void parsec_sparse_matrix_destroy( parsec_sparse_matrix_desc_t *desc );

/**
 * @}
 * @name PaRSEC vector descriptor
 * @{
 *
 *   @brief PaRSEC descriptor for the vectors linked to a given sparse matrix.
 */
typedef struct parsec_sparse_vector_desc_t {
    parsec_data_collection_t  super;  /**< Every PaRSEC descriptors must inherit from parsec_desc_t                        */
    int                       typesze; /**< Arithmetic size                                                                 */
    SolverMatrix             *solvmtx; /**< Solver matrix structure that describes the problem and stores the original data */
} parsec_sparse_vector_desc_t;


void parsec_sparse_vector_init( parsec_sparse_vector_desc_t *desc,
                                int typesze, int nodes, int myrank );
void parsec_sparse_vector_destroy( parsec_sparse_vector_desc_t *desc );

/**
 * @}
 * @name PaRSEC control function
 * @{
 *
 */
void pastix_parsec_init( pastix_data_t *pastix,
                         int *argc, char **argv[],
                         const int *bindtab );
void pastix_parsec_finalize( pastix_data_t *pastix );

/**
 * @}
 */
#endif /* _pastix_parsec_h_ */

/**
 *@}
 */
