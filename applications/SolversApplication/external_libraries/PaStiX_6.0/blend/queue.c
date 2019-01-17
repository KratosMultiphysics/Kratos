/**
 *
 * @file queue.c
 *
 * PaStiX queue structure.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @addtogroup blend_dev_queue
 * @{
 *
 **/
#include <stdio.h>
#include "common.h"
#include "queue.h"

/**
 *******************************************************************************
 *
 * @brief Compare two items of the queue based on the couple (key1, key2).
 *
 *******************************************************************************
 *
 * @param[in] item1
 *          The first item to compare.
 *
 * @param[in] item2
 *          The second item to compare.
 *
 *******************************************************************************
 *
 * @retval 1 if item1 < item2
 * @retval 0 otherwise
 *
 *******************************************************************************/
static inline int
pqueueItemComparison(const pastix_queue_item_t *item1,
                     const pastix_queue_item_t *item2)
{
    /* if elt1 < elt2 return 1  */
    /* if elt1 = elt2 return 0  */
    /* if elt1 > elt2 return 0 */

    if ( item1->key1 == item2->key1) {
        return item1->key2 < item2->key2;
    }
    else {
        return item1->key1 < item2->key1;
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the queue structure with an initial space to store the
 * elements.
 *
 *******************************************************************************
 *
 * @param[inout] q
 *          The qallocated pointer to the queue to initialize.
 *
 * @param[in] size
 *          The inital size for the queue.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if malloc failed.
 *
 *******************************************************************************/
int
pqueueInit(pastix_queue_t *q,
           pastix_int_t    size)
{
    q->size = size;
    q->used = 0;
    if (q->size != 0)
    {
        MALLOC_INTERN(q->elttab, size, pastix_queue_item_t);
    }
    else
    {
        q->elttab  = NULL;
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Free the structure associated to the queue.
 *
 *******************************************************************************
 *
 * @param[inout] q
 *          The pointer to the queue to free.
 *
 *******************************************************************************/
void
pqueueExit(pastix_queue_t *q)
{
    if(q->size != 0)
    {
        memFree_null(q->elttab);
    }
    q->size = 0;
}

/**
 *******************************************************************************
 *
 * @brief Return the size of the queue.
 *
 *******************************************************************************
 *
 * @param[in] q
 *          The pointer to the queue.
 *
 *******************************************************************************
 *
 * @return The size of the queue.
 *
 *******************************************************************************/
pastix_int_t
pqueueSize(const pastix_queue_t *q)
{
  return q->used;
}

/**
 *******************************************************************************
 *
 * @brief Reset the number of used element to 0.
 *
 *******************************************************************************
 *
 * @param[inout] q
 *          The pointer to the queue.
 *
 *******************************************************************************/
void
pqueueClear(pastix_queue_t *q)
{
  q->used = 0;
}

/**
 *******************************************************************************
 *
 * @brief Insert an element into the sorted queue.
 *
 *******************************************************************************
 *
 * @param[inout] q
 *          The pointer to the queue.
 *
 * @param[in] elt
 *          The element to insert in the queue.
 *
 * @param[in] key1
 *          The first key of the element.
 *
 * @param[in] key2
 *          The second key of the element.
 *
 *******************************************************************************/
void
pqueuePush2(pastix_queue_t *q,
            pastix_int_t    elt,
            double          key1,
            double          key2)
{
    pastix_int_t i, hi;

    /* Allocate more space if necessary */
    if(q->size == q->used)
    {
        pastix_queue_item_t *tmp;
        tmp = q->elttab;

        assert( (q->size == 0) || (tmp != NULL) );

        MALLOC_INTERN(q->elttab, q->size*2+1, pastix_queue_item_t);
        memcpy(q->elttab, tmp, q->size * sizeof(pastix_queue_item_t));

        q->size = q->size*2 +1;
        memFree_null(tmp);
    }

    q->elttab[q->used].key1   = key1;
    q->elttab[q->used].key2   = key2;
    q->elttab[q->used].eltptr = elt;
    q->used++;

    i = q->used - 1;
    hi= (i+1)/2-1;

    while( (i > 0) &&
           pqueueItemComparison(q->elttab + i,
                                q->elttab + hi) )
    {
        pastix_queue_item_t swap = q->elttab[i];

        q->elttab[i ] = q->elttab[hi];
        q->elttab[hi] = swap;

        i = hi; hi = (i+1)/2-1;
    }
}

/**
 *******************************************************************************
 *
 * @brief Read the first element of the queue.
 *
 *******************************************************************************
 *
 * @param[in] q
 *          The pointer to the queue.
 *
 *******************************************************************************
 *
 * @return The value of the first element sorted by (key1, key2).
 *
 *******************************************************************************/
pastix_int_t
pqueueRead(const pastix_queue_t *q)
{
    return q->elttab[0].eltptr;
}

/**
 *******************************************************************************
 *
 * @brief Remove the first element of the queue and return its keys if needed.
 *
 *******************************************************************************
 *
 * @param[inout] q
 *          The pointer to the queue. On exit, the queue without its head.
 *
 * @param[out] key1
 *          If key1 != NULL, stores the associated key1 to the first element on
 *          exit.
 *
 * @param[out] key2
 *          If key2 != NULL, stores the associated key2 to the first element on
 *          exit.
 *
 *******************************************************************************
 *
 * @return The value of the first element sorted by (key1, key2).
 *
 *******************************************************************************/
pastix_int_t
pqueuePop2(pastix_queue_t *q, double *key1, double*key2)
{
    pastix_int_t i, j;
    pastix_int_t return_elt;

    if (q->used == 0) {
        return -1;
    }

    return_elt = q->elttab[0].eltptr;
    if (key1 != NULL) { *key1 = q->elttab[0].key1; }
    if (key2 != NULL) { *key2 = q->elttab[0].key2; }

    q->elttab[0] = q->elttab[q->used-1];
    q->used--;

    i = 1;
    while(i <= (q->used/2))
    {
        if( (2*i == q->used)
            || pqueueItemComparison(q->elttab + 2*i-1,
                                    q->elttab + 2*i  ) )     /*(q->keytab[2*i-1] < q->keytab[2*i]))*/
        {
            j = 2*i;
        }
        else
        {
            j = 2*i+1;
        }
        if (!pqueueItemComparison(q->elttab + i-1,
                                  q->elttab + j-1))         /*(q->keytab[i-1] >= q->keytab[j-1])*/
        {
            pastix_queue_item_t swap;

            swap           = q->elttab[i-1];
            q->elttab[i-1] = q->elttab[j-1];
            q->elttab[j-1] = swap;

            i = j;
        }
        else {
            break;
        }
    }
    return return_elt;
}

/**
 *******************************************************************************
 *
 * @brief Print the queue.
 *
 * Print the queue elements for debug. They are not printed in a sorted order,
 * but in the order they are stored in the queue.
 *
 *******************************************************************************
 *
 * @param[in] q
 *          The pointer to the queue.
 *
 *******************************************************************************/
void
pqueuePrint(const pastix_queue_t *q)
{
    pastix_queue_item_t *item = q->elttab;
    pastix_int_t i;

    fprintf(stderr, "Queue :\n");
    for (i = 0; i < q->used; i++, item++) {
        fprintf(stderr, "(%ld %ld %ld) ",
                (long)(item->eltptr),
                (long)(item->key1),
                (long)(item->key2) );
        if (i%4 == 3) {
            fprintf(stderr, "\n");
        }
    }
    fprintf(stderr, "\n");
}

/**
 *@}
 */
