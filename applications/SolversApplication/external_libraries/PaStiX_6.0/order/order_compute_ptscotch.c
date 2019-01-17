/**
 *
 * @file order_compute_ptscotch.c
 *
 * PaStiX order driver to perform ordering with PT-Scotch library.
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
#if defined(PASTIX_ORDERING_PTSCOTCH)
#include <ptscotch.h>
#endif
#include "order_scotch_strats.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Compute the ordering of the graph given as parameter
 * with PT-Scotch library.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_ORDERING_DEFAULT, IPARM_SCOTCH_SWITCH_LEVEL,
 *   IPARM_SCOTCH_CMIN, IPARM_SCOTCH_CMAX, IPARM_SCOTCH_FRAT
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
 * @retval PASTIX_ERR_INTEGER_TYPE if Scotch integer type is not the
 *         same size as PaStiX ones,
 * @retval PASTIX_ERR_INTERNAL if an error occurs internally to Scotch.
 *
 *******************************************************************************/
int
pastixOrderComputePTScotch( pastix_data_t  *pastix_data,
                            pastix_graph_t *graph )
{
    SCOTCH_Dordering ordedat;
    SCOTCH_Ordering  ordering;
    SCOTCH_Dgraph    dgraph;
    SCOTCH_Strat     stratdat;
#if defined(PERSONAL_PTSCOTCH_STRATEGY)
    char             strat[1024];
#endif
    MPI_Comm         pastix_comm;
    pastix_order_t  *ordemesh = pastix_data->ordemesh;
    Clock            timer;
    pastix_int_t    *colptr;
    pastix_int_t    *rows;
    pastix_int_t    *iparm = pastix_data->iparm;
    pastix_int_t     procnum;
    pastix_int_t     n, gN, nnz, baseval;

    procnum     = pastix_data->procnum;
    pastix_comm = pastix_data->pastix_comm;

    /* Check integer compatibility */
    if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num)) {
        errorPrint("pastixOrderComputePTScotch: Inconsistent integer type between Pastix and PT-Scotch\n");
        return PASTIX_ERR_INTEGER_TYPE;
    }

     /* Enable this define to fix the SCOTCH random generator */
#if defined(SCOTCH_FIX_SEED)
    _SCOTCHintRandInit();
#endif

    gN     = graph->gN;
    n      = graph->n;
    colptr = graph->colptr;
    rows   = graph->rows;
    baseval= colptr[0];
    nnz    = colptr[n] - baseval;

    /* Build distributed graph */
    SCOTCH_dgraphInit(&dgraph, pastix_comm);

    if ( SCOTCH_dgraphBuild (&dgraph,
                             baseval,      /* baseval */
                             n,            /* number of local vertices */
                             n,            /* Maximum number of local vertices     */
                             colptr,
                             NULL,
                             NULL,         /* Local vertex load array (if any)     */
                             NULL,         /* Local vertex label array (if any)    */
                             nnz,
                             nnz,
                             rows,         /* Local edge array                     */
                             NULL,         /* Ghost edge array (if any); not const */
                             NULL))
    {
        errorPrint("SCOTCH_dgraphBuild");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    if (SCOTCH_dgraphCheck(&dgraph)) {
        errorPrint("pastix: SCOTCH_dgraphCheck");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    if (SCOTCH_stratInit(&stratdat))
    {
        errorPrint("pastix : SCOTCH_stratInit");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    /*
     * Create Strategy string for Scotch
     */
    /* TODO : Add default strategies for PT-Scotch */
#if defined(PERSONAL_PTSCOTCH_STRATEGY)
    if (iparm[IPARM_ORDERING_DEFAULT] == 1)
    {
        if (iparm[IPARM_INCOMPLETE] == 0)
            sprintf(strat, PTSCOTCH_STRAT_DIRECT);
        else
            sprintf(strat, PTSCOTCH_STRAT_INCOMP);
    }
    else /* Personal strategy */
    {
        sprintf(strat,  PTSCOTCH_STRAT_PERSO,
                (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                (long)  iparm[IPARM_SCOTCH_CMIN],
                (long)  iparm[IPARM_SCOTCH_CMAX],
                ((float)iparm[IPARM_SCOTCH_FRAT])/100.,
                (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                (long)  iparm[IPARM_SCOTCH_CMIN],
                (long)  iparm[IPARM_SCOTCH_CMAX],
                ((float)iparm[IPARM_SCOTCH_FRAT])/100.);

        if (iparm[IPARM_VERBOSE] > PastixVerboseNo)
            pastix_print(procnum, 0, "PT-Scotch Strategy |%s|\n", strat);
    }

    if (SCOTCH_stratDgraphOrder(&stratdat, strat))
    {
        errorPrint("pastix : SCOTCH_stratDgraphOrder");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }
#else
    if (iparm[IPARM_VERBOSE] > PastixVerboseNo)
        pastix_print(procnum, 0, "PaStiX works only with PT-Scotch default strategy %s", "");
#endif

    clockStart(timer);
    if (0 != SCOTCH_dgraphOrderInit(&dgraph, &ordedat))
    {
        errorPrint("pastix : SCOTCH_dgraphOrderInit");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    if (0 != SCOTCH_dgraphOrderCompute(&dgraph, &ordedat, &stratdat))
    {
        errorPrint("pastix : SCOTCH_dgraphOrderCompute");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    SCOTCH_stratExit(&stratdat);

    /*       if (0 != SCOTCH_dgraphOrderPerm(dgraph, ordedat, perm)) */
    /*  { */
    /*    errorPrint("pastix : SCOTCH_dgraphOrderPerm"); */
    /*    EXIT(MOD_SOPALIN,INTERNAL_ERR); */
    /*  } */

    /* TODO: pastixOrderInit prototype has changed !!! */
    assert(0);
    pastixOrderInit(ordemesh, gN, gN);
    memset( ordemesh->rangtab, 0, (gN+1)*sizeof(pastix_int_t));

    SCOTCH_dgraphCorderInit (&dgraph,
                             &ordering,
                             (SCOTCH_Num *)(ordemesh->permtab),
                             (SCOTCH_Num *)(ordemesh->peritab),
                             &ordemesh->cblknbr,
                             ordemesh->rangtab,
                             NULL);

    if (procnum == 0) {
        SCOTCH_dgraphOrderGather (&dgraph, &ordedat, &ordering);
    }
    else {
        SCOTCH_dgraphOrderGather (&dgraph, &ordedat, NULL);
    }

    ordemesh->baseval = dgraph->baseval;
    MPI_Bcast(&ordemesh->cblknbr, 1,                     PASTIX_MPI_INT, 0, pastix_comm);
    MPI_Bcast( ordemesh->rangtab, (ordemesh->cblknbr+1), PASTIX_MPI_INT, 0, pastix_comm);
    MPI_Bcast( ordemesh->permtab, gN,                    PASTIX_MPI_INT, 0, pastix_comm);
    MPI_Bcast( ordemesh->peritab, gN,                    PASTIX_MPI_INT, 0, pastix_comm);

    SCOTCH_dgraphCorderExit( &dgraph, &ordering );
    SCOTCH_dgraphOrderExit( &dgraph, &ordedat );
    SCOTCH_dgraphExit( &dgraph );

    pastixOrderBase(ordemesh, 0);

    if (iparm[IPARM_VERBOSE] > PastixVerboseNot)
        pastix_print(procnum, 0, TIME_COMPUTE_ORDERING, clockVal(timer));

    clockStop(timer);
    if (iparm[IPARM_VERBOSE] > PastixVerboseYes)
        pastix_print(procnum, 0, "%s", OUT_ORDERINIT);

    return PASTIX_SUCCESS;
}
