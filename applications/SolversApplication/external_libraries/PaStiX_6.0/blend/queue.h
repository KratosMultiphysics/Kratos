/**
 *
 * @file queue.h
 *
 * PaStiX queue structure header.
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
 *    This module describes the queue structure used in the analyze part of the solver.
 *    The sorting is based on a balanced tree that is partially updated at
 *    insertion and suppression.
 *
 **/
#ifndef _queue_h_
#define _queue_h_

/**
 * @brief Queue item structure.
 */
typedef struct pastix_queue_item_s {
    double       key1;   /**< Key 1 of the element   */
    double       key2;   /**< Key 2 of the element   */
    pastix_int_t eltptr; /**< Pointer to the element */
} pastix_queue_item_t;

/**
 * @brief Queue structure.
 */
typedef struct pastix_queue_s {
    pastix_int_t         size;   /**< Allocated memory size          */
    pastix_int_t         used;   /**< Number of element in the queue */
    pastix_queue_item_t *elttab; /**< Array of the element           */
} pastix_queue_t;

int          pqueueInit(        pastix_queue_t *, pastix_int_t );
void         pqueueExit(        pastix_queue_t * );
pastix_int_t pqueueSize(  const pastix_queue_t * );
void         pqueueClear(       pastix_queue_t * );
void         pqueuePush2(       pastix_queue_t *, pastix_int_t, double, double );
pastix_int_t pqueueRead ( const pastix_queue_t * );
pastix_int_t pqueuePop2 (       pastix_queue_t *, double *, double * );
void         pqueuePrint( const pastix_queue_t * );

/**
 * @brief Push an element with a single key.
 * @param[inout] q
 *               The queue structure.
 * @param[in]    elt
 *               The element to insert.
 * @param[in]    key1
 *               The first key of the element to insert (the second will be 0.).
 */
static inline void
pqueuePush1(pastix_queue_t *q, pastix_int_t elt, double key1) {
    pqueuePush2( q, elt, key1, 0. );
}

/**
 * @brief Pop the head of the queue whithout returning the keys.
 * @param[inout] q
 *               The queue structure.
 * @return The element at the head of the queue.
 */
static inline pastix_int_t
pqueuePop(pastix_queue_t *q){
    return pqueuePop2(q, NULL, NULL);
}

/**
 * @brief Pop the head of the queue and get the associated first key.
 * @param[inout] q
 *               The queue structure.
 * @param[out]   key1
 *               The first key of the element removed from the head of the queue.
 * @return The element at the head of the queue.
 */
static inline pastix_int_t
pqueuePop1(pastix_queue_t *q, double *key1){
    return pqueuePop2(q, key1, NULL);
}

#endif /* _queue_h_ */

/**
 * @}
 */
