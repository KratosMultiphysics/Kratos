/**
 *
 * @file order_check.c
 *
 * PaStiX order function to check correctness of a ordering structure.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "pastix/order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine checks the correctness of the ordering structure.
 *
 *******************************************************************************
 *
 * @param[in] ordeptr
 *          The ordering structure to check.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if the ordering structure is incorrect.
 *
 *******************************************************************************/
int
pastixOrderCheck (const pastix_order_t * const  ordeptr)
{
    pastix_int_t         baseval;                  /* Node base value            */
    pastix_int_t         cblkmax;                  /* Maximum supernode value    */
    pastix_int_t         vnodmax;                  /* Maximum node value         */
    pastix_int_t         vnodnum;                  /* Number of current node     */
    pastix_int_t         rangnum;                  /* Current column block index */
    const pastix_int_t * peritax;                  /* Based access to peritab    */
    const pastix_int_t * permtax;                  /* Based access to permtab    */

    /* Parameter checks */
    if (ordeptr == NULL) {
        errorPrint ("pastixOrderCheck: invalid ordeptr pointer");
        return PASTIX_ERR_BADPARAMETER;
    }

    if (ordeptr->cblknbr < 0) {
        errorPrint ("pastixOrderCheck: invalid nunber of column blocks");
        return PASTIX_ERR_BADPARAMETER;
    }

    baseval = ordeptr->baseval;
    if (baseval < 0) {
        errorPrint ("pastixOrderCheck: invalid vertex node base number");
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( (ordeptr->permtab == NULL) && (ordeptr->vertnbr > 0) ) {
        errorPrint ("pastixOrderCheck: permtab array is missing");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( (ordeptr->peritab == NULL) && (ordeptr->vertnbr > 0) ) {
        errorPrint ("pastixOrderCheck: peritab array is missing");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (ordeptr->rangtab == NULL) {
        errorPrint ("pastixOrderCheck: rangtab array is missing");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( (ordeptr->treetab == NULL) && (ordeptr->cblknbr > 0) ) {
        errorPrint ("pastixOrderCheck: treetab array is missing");
        return PASTIX_ERR_BADPARAMETER;
    }

    assert(baseval == ordeptr->rangtab[0]);

    peritax = ordeptr->peritab - baseval; /* Set based accesses */
    vnodmax = ordeptr->rangtab[ordeptr->cblknbr] - 1;

    assert((ordeptr->rangtab[ordeptr->cblknbr] - baseval) == ordeptr->vertnbr);

    /*
     * Check the values in rangtab
     */
    for (rangnum = 0; rangnum < ordeptr->cblknbr; rangnum ++)
    {
        if ((ordeptr->rangtab[rangnum] <  baseval) ||
            (ordeptr->rangtab[rangnum] >  vnodmax) ||
            (ordeptr->rangtab[rangnum] >= ordeptr->rangtab[rangnum + 1]))
        {
            errorPrint ("pastixOrderCheck: invalid range array");
            return PASTIX_ERR_BADPARAMETER;
        }
    }

    permtax = ordeptr->permtab - baseval;

    /*
     * Check perm and invp, as well as the symmetry between the two
     */
    for (vnodnum = baseval;
         vnodnum <= vnodmax; vnodnum ++)
    {
        pastix_int_t vnodold;

        vnodold = peritax[vnodnum];
        if ((vnodold < baseval) ||
            (vnodold > vnodmax) ||
            (permtax[vnodold] != vnodnum))
        {
            errorPrint ("pastixOrderCheck: invalid permutation arrays");
            return PASTIX_ERR_BADPARAMETER;
        }
    }

    /*
     * Check the treetab
     */
    cblkmax = ordeptr->cblknbr - baseval;
    for (rangnum = 0; rangnum < ordeptr->cblknbr-1; rangnum ++)
    {
        if ((ordeptr->treetab[rangnum] > cblkmax    ) ||
            ((ordeptr->treetab[rangnum] != -1)       &&
             (ordeptr->treetab[rangnum]  < (baseval+rangnum)) ) )
        {
            errorPrint ("pastixOrderCheck: invalid range array in treetab");
            return PASTIX_ERR_BADPARAMETER;
        }
    }
    if ((rangnum > 0) && (ordeptr->treetab[rangnum] != (baseval-1)))
    {
        errorPrint ("pastixOrderCheck: invalid father for cblknbr-1 node");
        return PASTIX_ERR_BADPARAMETER;
    }

    return PASTIX_SUCCESS;
}
