/**
 *
 * @file extracblk.h
 *
 * PaStiX analyse headers for extra symbolic structure functions.
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_split
 * @{
 *    This module contains all subroutines to initialize the candidates array
 *    for each supernode, as well as supernode properties that are defined by
 *    level such as 2D layouts and 2D tasks.
 *
 **/
#ifndef _extracblk_h_
#define _extracblk_h_

/**
 * @brief Extra symbol cblk structure
 */
typedef struct extracblk_s {
    pastix_int_t   cblknbr; /**< Number of cblk allocated                          */
    pastix_int_t   addcblk; /**< Number of cblk created                            */
    pastix_int_t   addblok; /**< Number of blok created                            */
    pastix_int_t   addblof; /**< Number of blok created due to facing cblk splited */
    pastix_int_t  *sptcblk; /**< Index for splitted cblk in the cblktab            */
    pastix_int_t  *sptcbnb; /**< Number of splitted cblk for a cblk                */
    pastix_int_t   curcblk; /**< Cursor for cblktab                                */
    pastix_int_t   sizcblk; /**< Size of allocated cblktab                         */
    symbol_cblk_t *cblktab; /**< Array of column blocks [+1,based]                 */
} ExtraCblk_t;

void extraCblkInit ( pastix_int_t        cblknbr,
                     ExtraCblk_t        *extracblk );
void extraCblkExit ( ExtraCblk_t        *extracblk );
void extraCblkAdd  ( ExtraCblk_t        *extracblk,
                     pastix_int_t        fcolnum,
                     pastix_int_t        lcolnum   );
void extraCblkMerge( const ExtraCblk_t  *extracblk,
                     symbol_matrix_t    *newsymb,
                     Cand              **candtab   );

#endif /* _extracblk_h_ */

/**
 *@}
 */
