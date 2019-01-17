/**
 *
 * @file graph_isolate.c
 *
 * PaStiX graph isolate routine
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

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Isolate a subset of vertices from a given graph.
 *
 * Return a new graph cleaned from those vertices.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of columns of the original GRAPH matrix.
 *
 * @param[in] colptr
 *          Array of size n+1.
 *          Index of first element of each column in rows array.
 *
 * @param[in] rows
 *          Array of size nnz = colptr[n] - colptr[0].
 *          Rows of each non zero entries.
 *
 * @param[in] isolate_n
 *          The number of columns to isolate from the original graph.
 *
 * @param[inout] isolate_list
 *          Array of size isolate_n.
 *          List of columns to isolate. On exit, the list is sorted by ascending
 *          indexes. Must be based as the graph.
 *
 * @param[out] new_colptr
 *          Array of size n-isolate_n+1.
 *          Index of first element of each column in rows array for the new graph.
 *          If new_colptr == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure based as the input colptr.
 *
 * @param[out] new_rows
 *          Array of size new_nnz = (*new_colptr)[n] - (*new_colptr)[0].
 *          Rows of each non zero entries for the new graph.
 *          If new_rows == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure based as the input rows.
 *
 * @param[out] new_perm
 *          Array of size n-isolate_n.
 *          Contains permutation generated to isolate the columns at the end of
 *          the graph that is 0-based.
 *          If new_perm == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_invp
 *          Array of size n-isolate_n.
 *          Contains the inverse permutation generated to isolate the columns
 *          at the end of the graph that is 0-based.
 *          If new_invp == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success,
 * @retval PASTIX_ERR_ALLOC if allocation went wrong,
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int graphIsolate(       pastix_int_t   n,
                  const pastix_int_t  *colptr,
                  const pastix_int_t  *rows,
                        pastix_int_t   isolate_n,
                        pastix_int_t  *isolate_list,
                        pastix_int_t **new_colptr,
                        pastix_int_t **new_rows,
                        pastix_int_t **new_perm,
                        pastix_int_t **new_invp )
{
    pastix_int_t *tmpcolptr = NULL;
    pastix_int_t *tmprows   = NULL;
    pastix_int_t *tmpperm   = NULL;
    pastix_int_t *tmpinvp   = NULL;
    pastix_int_t  baseval = colptr[0];
    pastix_int_t  nnz = colptr[n] - baseval;
    pastix_int_t  new_n = n - isolate_n;
    pastix_int_t  new_nnz;
    pastix_int_t  i, j, ip, k;
    pastix_int_t  iter_isolate = 0;
    pastix_int_t  iter_non_isolate  = 0;

    if ( (isolate_n > n)  || (isolate_n < 0) ) {
        errorPrintW( "Number of columns to isolate greater than the columns in the GRAPH matrix\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Quick Return */
    if (isolate_n == 0) {
        if (new_colptr != NULL) *new_colptr = (pastix_int_t*)colptr;
        if (new_rows   != NULL) *new_rows   = (pastix_int_t*)rows;
        return PASTIX_SUCCESS;
    }

    if (isolate_n == n) {
        if (new_colptr != NULL) {
            MALLOC_INTERN(*new_colptr, n, pastix_int_t);
            memcpy( *new_colptr, colptr, n*sizeof(pastix_int_t) );
        }
        if (new_rows != NULL) {
            MALLOC_INTERN(*new_rows, nnz, pastix_int_t);
            memcpy( *new_rows, rows, nnz*sizeof(pastix_int_t) );
        }
        if (new_perm != NULL) {
            MALLOC_INTERN(*new_perm, n, pastix_int_t);
            for (i = 0; i < n; i++) {
                (*new_perm)[i] = i;
            }
        }
        if (new_invp != NULL) {
            MALLOC_INTERN(*new_invp, n, pastix_int_t);
            for (i = 0; i < n; i++) {
                (*new_invp)[i] = i;
            }
        }
        return PASTIX_SUCCESS;
    }

    /* Sort the lost of vertices */
    intSort1asc1(isolate_list, isolate_n);

    /* Init invp array */
    MALLOC_INTERN(tmpinvp, n, pastix_int_t);
    for (i = 0; i <n; i++) {
        if ((iter_isolate < isolate_n) &&
            (i == isolate_list[iter_isolate]-baseval))
        {
            tmpinvp[new_n+iter_isolate] = i;
            iter_isolate++;
        }
        else
        {
            tmpinvp[iter_non_isolate] = i;
            iter_non_isolate++;
        }
    }

    assert(iter_non_isolate == new_n    );
    assert(iter_isolate     == isolate_n);

    /* Init perm array */
    MALLOC_INTERN(tmpperm, n, pastix_int_t);
    for(i = 0; i < n; i++)
        tmpperm[tmpinvp[i]] = i;

#if defined(PASTIX_DEBUG_GRAPH)
    for(i = 0; i < n; i++)
    {
        assert(tmpperm[i] < n );
        assert(tmpperm[i] > -1);
    }
#endif

    /* Create the new_colptr array */
    MALLOC_INTERN(tmpcolptr, new_n + 1, pastix_int_t);
    memset(tmpcolptr, 0, (new_n + 1)*sizeof(pastix_int_t));

    tmpcolptr[0] = baseval;
    for (i=0; i<n; i++)
    {
        ip = tmpperm[i];
        if (ip < new_n)
        {
            for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j++)
            {
                /* Count edges in each column of the new graph */
                if (tmpperm[rows[j]-baseval] < new_n)
                {
                    tmpcolptr[ip+1]++;
                }
            }
        }
    }

    for (i = 0; i < new_n; i++)
        tmpcolptr[i+1] += tmpcolptr[i];

    new_nnz = tmpcolptr[new_n] - tmpcolptr[0];

    /* Create the new rows array */
    if ( new_nnz != 0 ) {
        MALLOC_INTERN(tmprows, new_nnz, pastix_int_t);
        for (i = 0; i <n; i++)
        {
            ip = tmpperm[i];
            if (ip < new_n)
            {
                k = tmpcolptr[ip]-baseval;
                for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j ++)
                {
                    /* Count edges in each column of the new graph */
                    if (tmpperm[rows[j]-baseval] < new_n)
                    {
                        tmprows[k] = tmpperm[rows[j]-baseval] + baseval;
                        k++;
                    }
                }
                assert( k == tmpcolptr[ip+1]-baseval );
            }
        }
    }

    if (new_colptr != NULL) {
        *new_colptr = tmpcolptr;
    } else {
        memFree_null( tmpcolptr );
    }
    if (new_rows != NULL) {
        *new_rows = tmprows;
    } else {
        memFree_null( tmprows );
    }
    if (new_perm != NULL) {
        *new_perm = tmpperm;
    } else {
        memFree_null( tmpperm );
    }
    if (new_invp != NULL) {
        *new_invp = tmpinvp;
    } else {
        memFree_null( tmpinvp );
    }

    return PASTIX_SUCCESS;
}

