/**
 *
 * @file symbol_kass_csr.c
 *
 * PaStiX symbolic factorization routines
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
#include "graph.h"
#include "symbol_kass.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_kass
 *
 * @brief Initialize the data structure by doing the first allocations within
 * the structure and initializing the fields.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The size of the graph that needs to be initialized.
 *
 * @param[out] csr
 *          The graph to initialize.
 *
 *******************************************************************************/
void
kass_csrInit( pastix_int_t n, kass_csr_t *csr )
{
    csr->n = n;
    MALLOC_INTERN( csr->nnz,  n, pastix_int_t  );
    MALLOC_INTERN( csr->rows, n, pastix_int_t* );

    memset( csr->nnz,  0, n*sizeof(pastix_int_t ) );
    memset( csr->rows, 0, n*sizeof(pastix_int_t*) );
}


/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_kass
 *
 * @brief Free the data store in the structure.
 *
 *******************************************************************************
 *
 * @param[inout] csr
 *          The graph to clean.
 *
 *******************************************************************************/
void
kass_csrClean( kass_csr_t *csr )
{
    pastix_int_t i;
    for(i=0; i< csr->n; i++) {
        if ( csr->nnz[i] != 0 ) {
            memFree_null( csr->rows[i] );
        }
    }
    memFree_null( csr->rows );
    memFree_null( csr->nnz  );
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_kass
 *
 * @brief Computes the number of non zero entries in the graph.
 *
 * It is using the following formula: nnz = sum( i=0..n, nnz[n] )
 * The formula must be post computed to adapt to presence of diagonal elements
 * or not, and to the symmetry of the graph.
 *
 *******************************************************************************
 *
 * @param[in] csr
 *          The graph on which the number of non zero entries is computed.
 *
 *******************************************************************************
 *
 * @retval The number of non zero entries.
 *
 *******************************************************************************/
pastix_int_t
kass_csrGetNNZ( const kass_csr_t *csr )
{
    pastix_int_t i, nnz;
    nnz = 0;
    for(i=0; i< csr->n; i++) {
        nnz += csr->nnz[i];
    }
    return nnz;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_kass
 *
 * @brief Compact a compressed graph.
 *
 * All nodes with no non zero entries are removed from the graph, the allocated
 * space is not adjusted.
 *
 *******************************************************************************
 *
 * @param[inout] csr
 *          The graph to compact.
 *          On entry, graph which might contain nodes with no non zero entries.
 *          On exit, all those nodes are suppressed and the compressed graph is
 *          returned.
 *
 *******************************************************************************/
void
kass_csrCompact( kass_csr_t *csr )
{
    pastix_int_t n = csr->n;
    pastix_int_t i, j;

    /* Look for first deleted node */
    for(i=0,j=0; i<n; i++,j++) {
        if (csr->nnz[i] == 0)
            break;
    }

    /* Compact the nodes */
    for(; i<n; i++) {
        if ( csr->nnz[i] > 0 ) {
            assert( j < i );
            csr->nnz[j]  = csr->nnz[i];
            csr->rows[j] = csr->rows[i];
            csr->nnz[i]  = 0;
            csr->rows[i] = NULL;
            j++;
        }
    }

    csr->n = j;
}

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_kass
 *
 * @brief Generate the graph of P*A from the graph of A and the
 * permutation vector.
 *
 *******************************************************************************
 *
 * @param[in] graphA
 *          The original graph Aon which the permutation will be applied.
 *
 * @param[in] perm
 *          Integer array of size graphA->n. Contains the permutation to apply to A.
 *
 * @param[inout] graphPA
 *          On entry, the initialized graph with size graphA->n.
 *          On exit, contains the graph of P A.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success.
 * @retval PASTIX_ERR_ALLOC if allocation went wrong.
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int
kass_csrGenPA( const pastix_graph_t *graphA,
               const pastix_int_t   *perm,
                     kass_csr_t     *graphPA )
{
    pastix_int_t *rowsPA, *rowsA;
    pastix_int_t i, j, ip;
    pastix_int_t baseval;
    pastix_int_t n = graphPA->n = graphA->n;

    baseval = graphA->colptr[0];

    /* Compute the number of nnz per vertex */
    for(i=0; i<n; i++)
    {
        ip = perm[i];
        /* Add diagonal (could be removed fro direct) */
        graphPA->nnz[ip] = graphA->colptr[i+1] - graphA->colptr[i] + 1;
    }

    /* Create the row vector */
    for(i=0; i<n; i++)
    {
        ip = perm[i] - baseval;

        MALLOC_INTERN( graphPA->rows[ip],
                       graphPA->nnz[ip],
                       pastix_int_t );

        rowsPA = graphPA->rows[ip];
        rowsA  = graphA->rows + graphA->colptr[i] - baseval;

        /* Add diagonal */
        *rowsPA = ip;
        rowsPA++;

        for(j=1; j<graphPA->nnz[ip]; j++, rowsPA++) {
            *rowsPA = perm[ *rowsA ];
            rowsA++;
        }

        intSort1asc1( graphPA->rows[ip],
                      graphPA->nnz[ip]);
    }
    return PASTIX_SUCCESS;
}
