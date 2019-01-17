/**
 *
 * @file extendVector.c
 *
 * PaStiX analyse functions for the extend vector structure.
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
 *
 **/
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "extendVector.h"

/**
 *******************************************************************************
 *
 * @brief Initialize the extendVector structure with the initial size given.
 *
 *******************************************************************************
 *
 * @param[inout] vec
 *          The allocated pointer to initialize.
 *
 * @param[in] size
 *          The initial size of the vector.
 *
 *******************************************************************************
 *
 * @return The pointer to the internal integer array.
 *
 *******************************************************************************/
pastix_int_t *
extendint_Init(ExtendVectorINT *vec, pastix_int_t size)
{
    vec->vecsize = size;
    vec->eltnbr  = 0;
    vec->inttab  = NULL;
    MALLOC_INTERN(vec->inttab, size, pastix_int_t);
    return vec->inttab;
}

/**
 *******************************************************************************
 *
 * @brief Free the extendVector structure.
 *
 *******************************************************************************
 *
 * @param[inout] vec
 *          The allocated pointer to free.
 *
 *******************************************************************************/
void
extendint_Exit(ExtendVectorINT *vec)
{
    if(vec->inttab != NULL) {
        memFree_null(vec->inttab);
    }
    /*memFree_null(vec);*/
}

/**
 *******************************************************************************
 *
 * @brief Add an element elt to the end of the vector.
 *
 *******************************************************************************
 *
 * @param[inout] vec
 *          The pointer to the extend vector. Might be extended if there is not
 *          enough space..
 *
 * @param[in] elt
 *          The value to store in the vector
 *
 *******************************************************************************/
void
extendint_Add(ExtendVectorINT *vec, pastix_int_t elt)
{
    vec->inttab[vec->eltnbr] = elt;
    extendint_incr(vec);
}

/**
 *******************************************************************************
 *
 * @brief Return the number of element stored in the vector.
 *
 *******************************************************************************
 *
 * @param[in] vec
 *          The extend vector structure.
 *
 *******************************************************************************
 *
 * @return The number of element stored.
 *
 *******************************************************************************/
pastix_int_t
extendint_Size(const ExtendVectorINT *vec)
{
  return vec->eltnbr;
}

/**
 *******************************************************************************
 *
 * @brief Return the element of index eltnum.
 *
 *******************************************************************************
 *
 * @param[in] vec
 *          The extend vector structure.
 *
 * @param[in] eltnum
 *          The index of the elemnt to return.
 *
 *******************************************************************************
 *
 * @return The value of the eltnum^th element.
 *
 *******************************************************************************/
pastix_int_t
extendint_Read(const ExtendVectorINT *vec, pastix_int_t eltnum)
{
    assert(eltnum <= vec->eltnbr);
    return vec->inttab[eltnum];
}

/**
 *******************************************************************************
 *
 * @brief Reallocate the vector to the given size.
 *
 *******************************************************************************
 *
 * @param[inout] vec
 *          The extend vector structure to realloc.
 *
 * @param[in] size
 *          The new size of the vector. If the size is smaller than the actual
 *          size of the vector, nothing is done. Otherwise the vectore is
 *          enlarged and the values are kept.
 *
 *******************************************************************************/
void
extendint_ToSize(ExtendVectorINT *vec, pastix_int_t size)
{
    extendint_Clear(vec);

    if(size <= vec->vecsize) {  /* there 's enough space */
        return;
    }

    if(vec->inttab != NULL) {
        memFree_null(vec->inttab);
    }

    MALLOC_INTERN(vec->inttab, size, pastix_int_t);
    vec->vecsize = size;
}

/**
 *******************************************************************************
 *
 * @brief Increment the number of element stored.
 *
 * This is the internal function that reallocate the vector.
 *
 *******************************************************************************
 *
 * @param[inout] vec
 *          The extend vector structure to realloc.
 *
 *******************************************************************************/
void
extendint_incr(ExtendVectorINT *vec)
{
    vec->eltnbr++;
    /** if the vector is not big enough, make it bigger !! **/
    if(!(vec->eltnbr < vec->vecsize))
    {
        pastix_int_t *tmp;
        tmp = vec->inttab;

        vec->vecsize = vec->vecsize + vec->vecsize/2 +1;

        MALLOC_INTERN(vec->inttab, vec->vecsize, pastix_int_t);
        memcpy(vec->inttab, tmp, sizeof(pastix_int_t)*vec->eltnbr);
        memFree_null(tmp);
    }
}

/**
 *******************************************************************************
 *
 * @brief Cleanup the vector.
 *
 * Set the number of elements to 0. Does NµOT free the sturcture.
 *
 *******************************************************************************
 *
 * @param[inout] vec
 *          The extend vector to clear.
 *
 *******************************************************************************/
void
extendint_Clear(ExtendVectorINT *vec)
{
    vec->eltnbr = 0;
}

/**
 *@}
 */
