/**
 *
 * @file order_compute_metis.c
 *
 * PaStiX order driver to perform ordering with Metis library
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
#include "pastix/order.h"
#include <metis.h>

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Compute the ordering of the graph given as parameter
 * with Metis library.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_ORDERING_DEFAULT, IPARM_METIS_CTYPE,
 *   IPARM_METIS_RTYPE, IPARM_METIS_NO2HOP, IPARM_METIS_NSEPS,
 *   IPARM_METIS_NITER, IPARM_METIS_UFACTOR, IPARM_METIS_COMPRESS,
 *   IPARM_METIS_CCORDER, IPARM_METIS_PFACTOR, IPARM_METIS_SEED,
 *   IPARM_METIS_DBGLVL.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field oerdemesh is initialize with the result of the
 *          ordering realized by Scotch.
 *
 * @param[inout] graph
 *          The graph prepared by graphPrepare function on which wwe want to
 *          perform the ordering. On exit, the graph might be rebased.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed,
 * @retval PASTIX_ERR_INTEGER_TYPE if Metis integer type is not the
 *         same size as PaStiX ones,
 * @retval PASTIX_ERR_INTERNAL if an error occurs internally to Metis.
 *
 *******************************************************************************/
int
pastixOrderComputeMetis( pastix_data_t  *pastix_data,
                         pastix_graph_t *graph )
{
    pastix_int_t   *iparm    = pastix_data->iparm;
    pastix_order_t *ordemesh = pastix_data->ordemesh;
    pastix_int_t    n;
    pastix_int_t    baseval = graph->colptr[0];
    idx_t opt[METIS_NOPTIONS];
    int   rc;

    if ( sizeof(pastix_int_t) != sizeof(idx_t)) {
        errorPrint("Inconsistent integer type between PaStiX and Metis\n");
        return PASTIX_ERR_INTEGER_TYPE;
    }

    /* Set of valid options for METIS_NodeND */
    METIS_SetDefaultOptions(opt);
    if (iparm[IPARM_ORDERING_DEFAULT] != 1) {
        opt[METIS_OPTION_CTYPE   ] = iparm[IPARM_METIS_CTYPE   ];
        opt[METIS_OPTION_RTYPE   ] = iparm[IPARM_METIS_RTYPE   ];
        opt[METIS_OPTION_NO2HOP  ] = iparm[IPARM_METIS_NO2HOP  ];
        opt[METIS_OPTION_NSEPS   ] = iparm[IPARM_METIS_NSEPS   ];
        opt[METIS_OPTION_NITER   ] = iparm[IPARM_METIS_NITER   ];
        opt[METIS_OPTION_UFACTOR ] = iparm[IPARM_METIS_UFACTOR ];
        opt[METIS_OPTION_COMPRESS] = iparm[IPARM_METIS_COMPRESS];
        opt[METIS_OPTION_CCORDER ] = iparm[IPARM_METIS_CCORDER ];
        opt[METIS_OPTION_PFACTOR ] = iparm[IPARM_METIS_PFACTOR ];
    }
    opt[METIS_OPTION_SEED     ] = iparm[IPARM_METIS_SEED];
    opt[METIS_OPTION_NUMBERING] = baseval;
    opt[METIS_OPTION_DBGLVL   ] = iparm[IPARM_METIS_DBGLVL];

    n = graph->n;
    rc = pastixOrderAlloc( ordemesh, graph->n, 0 );
    if (rc != PASTIX_SUCCESS )
    {
        errorPrint("pastixOrderComputeMetis: Error during odering initialization\n");
        return rc;
    }
    ordemesh->baseval = baseval;
    rc = METIS_NodeND( &n, graph->colptr, graph->rows, NULL,
                       opt, ordemesh->peritab, ordemesh->permtab);

    assert( n == graph->n );
    if (rc != METIS_OK )
    {
        errorPrint("pastixOrderComputeMetis: Invalid code returned by METIS_NodeND (%d)\n", rc);
        return PASTIX_ERR_INTERNAL;
    }

#if defined(PASTIX_DEBUG_ORDERING)
    {
        pastix_int_t i;
        for(i=0; i<n; i++) {
            assert( ordemesh->permtab[i] >= baseval );
            assert( ordemesh->permtab[i] <  (n+baseval) );
            assert( ordemesh->peritab[i] >= baseval );
            assert( ordemesh->peritab[i] <  (n+baseval) );
        }
    }
#endif
    return PASTIX_SUCCESS;
}
