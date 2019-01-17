/**
 *
 * @file symbol_kass.c
 *
 * PaStiX Kass symbolic factorization routines
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
#include "pastix/order.h"
#include "symbol.h"
#include "symbol_kass.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol

 * @brief Generate the symbol matrix whith an amalgamation algorithm.
 *
 * This version of the symbolic factorization almagamate small blocks together
 * in larger ones. The symbol matrix and the associated supernode partition are
 * returned.  See after for the different parameters.
 *
 *******************************************************************************
 *
 * @param[in] verbose
 *          Adjust the level of verbosity of the function
 *
 * @param[in] ilu
 *          - 1: incomplete factorization will be performed.
 *          - 0 : direct factorization will be performed.
 *
 * @param[in] levelk
 *          Unused if ilu == 0.
 *          - k >= 0: symbol matrix for ILU(k) factorization will be generated.
 *          - < 0: symbol matrix for direct factorization will be generated.
 *
 * @param[in] rat_cblk
 *          Must be >= 0. Fill ratio that limits the amalgamation process based
 *          on the graph structure.
 *
 * @param[in] rat_blas
 *          Must be >= rat_cblk. Fill ratio that limits the amalgamation process
 *          that merges blocks in order to reduce the BLAS computational time
 *          (see amalgamate() for further informations).
 *
 * @param[inout] symbmtx
 *          The symbol matrix structure to construct. On entry, the initialized
 *          structure (see pastixSymbolInit()). On exit, the symbol matrix generated
 *          after the amalgamation process.
 *
 * @param[inout] csc
 *          The original csc for which the symbol matrix needs to be generated.
 *          Rebase to C numbering on exit.
 *
 * @param[inout] orderptr
 *          The oder structure that contains the perm and invp array generated
 *          by the ordering step. The supernode partition might be initialized
 *          or not.
 *          On exit, it is rebased to c numbering and contains the updated
 *          perm/invp arrays as well as the supernode partition.
 *
 * @param[in] pastix_comm
 *          The PaStiX instance communicator.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on success.
 *          \retval PASTIX_ERR_ALLOC if allocation went wrong.
 *          \retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int
pastixSymbolKass(int verbose, int ilu, int levelk, int rat_cblk, int rat_blas,
                 symbol_matrix_t   *symbmtx,
                 pastix_graph_t *csc,
                 pastix_order_t *orderptr,
                 MPI_Comm        pastix_comm)
{
    kass_csr_t graphPA, graphL;
    pastix_int_t snodenbr;
    pastix_int_t *snodetab = NULL;
    pastix_int_t *streetab = NULL;
    pastix_int_t *ia       = NULL;
    pastix_int_t  i, j, n;
    pastix_int_t *perm;
    pastix_int_t *invp;
    pastix_int_t *invp2;
    pastix_int_t  newcblknbr;
    pastix_int_t *newrangtab = NULL;
    pastix_int_t *newtreetab = NULL;
    pastix_int_t  nnzA, nnzL;
    Clock timer;
    double nnzS;
    int procnum;

    MPI_Comm_rank(pastix_comm, &procnum);

    /* Check parameters correctness */
    if ( (ilu == 0) || (levelk < 0) ) {
        /* Forces levelk to -1 */
        levelk = -1;
    }
    if ( symbmtx == NULL ) {
        errorPrintW("symbolKass: wrong parameter symbmtx");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( csc == NULL ) {
        errorPrintW("symbolKass: wrong parameter csc");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( orderptr == NULL ) {
        errorPrintW("symbolKass: wrong parameter orderptr");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Convert Fortran to C numbering */
    graphBase( csc, 0 );
    pastixOrderBase( orderptr, 0 );

    n  = csc->n;
    ia = csc->colptr;
    perm     = orderptr->permtab;
    invp     = orderptr->peritab;
    snodenbr = orderptr->cblknbr;
    snodetab = orderptr->rangtab;
    streetab = orderptr->treetab;

    /* Create the graph of P A */
    kass_csrInit( n, &graphPA );
    kass_csrGenPA( csc, perm, &graphPA );

    if (verbose > PastixVerboseYes)
        pastix_print(procnum, 0,
                     "Level of fill = %ld\n"
                     "Amalgamation ratio: cblk = %d, blas = %d\n",
                     (long)levelk, rat_cblk, rat_blas);

    /*
     * Compute the graph of the factorized matrix L
     */
    /* Direct Factorization */
    if((ilu == 0) || (levelk == -1))
    {
        /*
         * (Re)compute the streetab
         */
        if (streetab == NULL) { MALLOC_INTERN( streetab, n, pastix_int_t ); }

        clockStart(timer);
        nnzL = kassFactDirect( &graphPA, snodenbr, snodetab, streetab, &graphL );
        clockStop(timer);
        if (verbose > PastixVerboseYes)
            pastix_print(procnum, 0,
                         "Time to compute scalar symbolic direct factorization  %.3g s\n",
                         clockVal(timer));
    }
    /* ILU(k) Factorization */
    else
    {
        clockStart(timer);
        nnzL = kassFactLevel( &graphPA, levelk, &graphL );
        clockStop(timer);
        if (verbose > PastixVerboseYes)
            pastix_print(procnum, 0,
                         "Time to compute scalar symbolic factorization of ILU(%ld) %.3g s\n",
                         (long)levelk, clockVal(timer));

        /*
         * Compute the treetab array for amalgamation step
         */
        {
            pastix_int_t *treetab;
            assert(streetab != NULL);

            MALLOC_INTERN(treetab, graphL.n, pastix_int_t);
            memcpy( treetab, streetab, graphL.n * sizeof(pastix_int_t) );
            for(j=0;j<snodenbr;j++)
            {
                for(i=snodetab[j]; i<snodetab[j+1]-1; i++)
                    streetab[i] = i+1;

                /* Generic version */
                if( (treetab[j] == -1) ||
                    (treetab[j] == j ) )
                {
                    streetab[i] = -1;
                }
                else
                {
                    streetab[i] = snodetab[treetab[j]];
                }
            }
            memFree(treetab);
        }

        snodenbr = -1;
        snodetab = NULL;
    }

    nnzA = ( kass_csrGetNNZ( &graphPA ) + n ) / 2;
    kass_csrClean( &graphPA );

    if (verbose > PastixVerboseYes)
        pastix_print( procnum, 0,
                      "Scalar nnza = %ld nnzlk = %ld, fillrate0 = %.3g \n",
                      (long)nnzA, (long)nnzL, (double)nnzL / (double)nnzA );

    /*
     * Amalgamate the blocks
     */
    clockStart(timer);
    MALLOC_INTERN(invp2, n, pastix_int_t);

    amalgamate( (double)rat_cblk / 100.,
                (double)rat_blas / 100.,
                &graphL, nnzL,
                snodenbr, snodetab, streetab,
                &newcblknbr, &newrangtab, &newtreetab,
                invp2, pastix_comm );

    if( orderptr->rangtab != NULL ) {
        memFree(orderptr->rangtab);
        orderptr->cblknbr = 0;
    }
    if( orderptr->treetab != NULL ) {
        memFree(orderptr->treetab);
    }
    orderptr->cblknbr = newcblknbr;
    orderptr->rangtab = newrangtab;
    orderptr->treetab = newtreetab;

    /*
     * invp2 is the invp vector generated for PA
     * we compute the invp and perm vectors of A
     */
    for(i=0;i<n;i++)
        invp2[i] = invp[invp2[i]];

    memcpy(invp, invp2, sizeof(pastix_int_t)*n);
    for(i=0;i<n;i++)
        perm[invp[i]] = i;
    memFree(invp2);

    /* Let's build the symbol matrix */
    kassBuildSymbol( &graphL, newcblknbr, newrangtab, symbmtx );
    kass_csrClean( &graphL );

    /* Patch the symbol matrix to have a real elimination tree */
    if (levelk != -1) {
        nnzS = pastixSymbolGetNNZ( symbmtx );

        if (verbose > PastixVerboseYes) {
            pastix_print(procnum, 0, "Number of blocks in the non patched symbol matrix = %ld \n",
                         (long)symbmtx->bloknbr);
            pastix_print(procnum, 0, "Number of non zeroes in the non patched symbol matrix = %g, fillrate1 %.3g \n",
                         nnzS+n, (nnzS+n)/(ia[n]/2.0 +n) );
        }
        if(pastixSymbolCheck(symbmtx) != 0) {
            errorPrint("SymbolCheck on symbol matrix before patch failed !!!");
            assert(0);
        }

        /* Patch the symbol matrix */
        kassPatchSymbol( symbmtx );
    }

    nnzS = pastixSymbolGetNNZ( symbmtx );

    clockStop(timer);
    if (verbose > PastixVerboseYes) {
        pastix_print(procnum, 0, "Time to compute the amalgamation of supernodes %.3g s\n", clockVal(timer));
        pastix_print(procnum, 0, "Number of cblk in the amalgamated symbol matrix = %ld \n", (long)newcblknbr);
        pastix_print(procnum, 0, "Number of block in final symbol matrix = %ld \n",
                     (long)symbmtx->bloknbr);
        pastix_print(procnum, 0, "Number of non zero in final symbol matrix = %g, fillrate2 %.3g \n",
                     nnzS+n, (nnzS+n)/(ia[n]/2.0 +n));
    }
    return PASTIX_SUCCESS;
}
