/**
 *
 * @file extendVector.h
 *
 * PaStiX analyse header for extend vector structure.
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_extint
 * @{
 *    This module contains all subroutines to manage the extend int vector
 *    arrays. This structure manage a vector of integer automatically extensible
 *    when elements are added. This is used in the simulation algorithm to store
 *    the available tasks per worker.
 *
 **/
#ifndef _extendvector_h_
#define _extendvector_h_

/**
 * @brief The extend integer array structure.
*/
typedef struct ExtendVectorINT_s {
    pastix_int_t  vecsize; /**< The size of the vector             */
    pastix_int_t  eltnbr;  /**< The number of elements stored      */
    pastix_int_t *inttab;  /**< The actual array with the elements */
} ExtendVectorINT;

pastix_int_t *extendint_Init  (       ExtendVectorINT *, pastix_int_t );
void          extendint_Exit  (       ExtendVectorINT * );
void          extendint_Add   (       ExtendVectorINT *, pastix_int_t );
pastix_int_t  extendint_Size  ( const ExtendVectorINT * );
pastix_int_t  extendint_Read  ( const ExtendVectorINT *, pastix_int_t );
void          extendint_Clear (       ExtendVectorINT * );
void          extendint_ToSize(       ExtendVectorINT *, pastix_int_t );
void          extendint_incr  (       ExtendVectorINT * );

#endif /* _extendvector_h_ */

/**
 *@}
 */
