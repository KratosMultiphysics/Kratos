/**
 *
 * @file cost.c
 *
 * PaStiX analyse functions for the cost matrix arrays.
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_cost
 * @{
 *
 **/
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "symbol.h"
#include "cost.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the cost matrix structure
 *
 *******************************************************************************
 *
 * @param[inout] costmtx
 *          The cost matrix structure to initialize.
 *
 *******************************************************************************/
void
costMatrixInit( CostMatrix *costmtx )
{
    costmtx->blokcost = NULL;
    costmtx->cblkcost = NULL;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Free the cost matrix structure
 *
 *******************************************************************************
 *
 * @param[inout] costmtx
 *          The cost matrix structure to free.
 *
 *******************************************************************************/
void
costMatrixExit( CostMatrix *costmtx )
{
    if(costmtx->blokcost != NULL) {
	memFree_null(costmtx->blokcost);
    }

    if(costmtx->cblkcost != NULL) {
	memFree_null(costmtx->cblkcost);
    }
}

/**
 *******************************************************************************
 *
 * @brief Build the cost matrix structure from the symbol matrix structure.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The symbol matrix structure.
 *
 * @param[in] flttype
 *          The floating point arithmetic that will be used to adapt the cost models.
 *
 * @param[in] factotype
 *          The factorization that will be appplied to adapt the cost models.
 *
 *******************************************************************************
 *
 * @return The cost matrix structure initialized with the cost of each block
 *         operation, and each cblk operation.
 *
 *******************************************************************************/
CostMatrix *
costMatrixBuild( const symbol_matrix_t *symbmtx,
                 pastix_coeftype_t      flttype,
                 pastix_factotype_t     factotype )
{
    CostMatrix *costmtx = NULL;

    MALLOC_INTERN(costmtx, 1, CostMatrix);
    costMatrixInit(costmtx);

    MALLOC_INTERN( costmtx->cblkcost, symbmtx->cblknbr, double );
    MALLOC_INTERN( costmtx->blokcost, symbmtx->bloknbr, double );

    pastixSymbolGetTimes( symbmtx, flttype, factotype,
                          costmtx->cblkcost, costmtx->blokcost );

    return costmtx;
}

/**
 * @}
 */
