/**
 *
 * @file isched.c
 *
 * @copyright 2008-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2010-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 *
 * PaStiX Internal Thread System routines
 * Contains basic functions to handle threads for internal schedulings.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 */
#include "common.h"
#include <pthread.h>
#include "isched_barrier.h"
#include "isched.h"

#if defined(HAVE_HWLOC)
#include "isched_hwloc.h"
#elif defined(HAVE_SCHED_SETAFFINITY)
#  include <linux/unistd.h>
#  include <sched.h>
#elif defined(MAC_OS_X)
#  include <mach/mach_init.h>
#  include <mach/thread_policy.h>
/**
 * Expose the hidden kernel interface.
 */
extern kern_return_t thread_policy_set( thread_t               thread,
                                        thread_policy_flavor_t flavor,
                                        thread_policy_t        policy_info,
                                        mach_msg_type_number_t count);
#endif  /* define(HAVE_HWLOC) */

#if ISCHED_IMPLEMENT_BARRIERS

int
isched_barrier_init(isched_barrier_t* barrier,
                    const void* attr,
                    unsigned int count)
{
    int rc;

    if( 0 != (rc = pthread_mutex_init(&(barrier->mutex), attr)) ) {
        return rc;
    }

    barrier->count      = count;
    barrier->curcount   = 0;
    barrier->generation = 0;
    if( 0 != (rc = pthread_cond_init(&(barrier->cond), NULL)) ) {
        pthread_mutex_destroy( &(barrier->mutex) );
        return rc;
    }
    return 0;
}

int
isched_barrier_wait(isched_barrier_t* barrier)
{
    int generation;

    pthread_mutex_lock( &(barrier->mutex) );
    if( (barrier->curcount + 1) == barrier->count) {
        barrier->generation++;
        barrier->curcount = 0;
        pthread_cond_broadcast( &(barrier->cond) );
        pthread_mutex_unlock( &(barrier->mutex) );
        return 1;
    }
    barrier->curcount++;
    generation = barrier->generation;
    for(;;) {
        pthread_cond_wait( &(barrier->cond), &(barrier->mutex) );
        if( generation != barrier->generation ) {
            break;
        }
    }
    pthread_mutex_unlock( &(barrier->mutex) );
    return 0;
}

int
isched_barrier_destroy(isched_barrier_t* barrier)
{
    pthread_mutex_destroy( &(barrier->mutex) );
    pthread_cond_destroy( &(barrier->cond) );
    barrier->count    = 0;
    barrier->curcount = 0;
    return 0;
}

#endif  /* ISCHED_IMPLEMENT_BARRIERS */

/**
 * Temporary stucture used at thread creation
 */
typedef struct __isched_init_s {
    isched_t        *global_ctx;
    int              rank;
    int              bindto;
} __isched_init_t;

void *isched_thread_init(void *ptr);
void *isched_thread_destroy(isched_thread_t *ptr);
void *isched_parallel_section(isched_thread_t *ctx);

/***************************************************************************//**
 *  Returns core id
 **/
void *
isched_thread_init(void *ptr)
{
    __isched_init_t *isched = (__isched_init_t*)(ptr);
    isched_thread_t *thread_ctx;

    MALLOC_INTERN( thread_ctx, 1, isched_thread_t );

    thread_ctx->global_ctx = isched->global_ctx;
    thread_ctx->rank       = isched->rank;

    /* Set thread affinity for the worker */
    if ( isched->bindto != -1 ) {
        isched_topo_bind_on_core_index( isched->bindto );
    }

    if ( thread_ctx->rank != 0 ) {
        return isched_parallel_section( thread_ctx );
    }
    else {
        return (void*)thread_ctx;
    }
}

void *
isched_thread_destroy(isched_thread_t *ctx)
{
    /* thread 0 is not bound */
    if ( ctx->rank != 0 ) {
        isched_topo_unbind();
    }
    memFree_null( ctx );

    return NULL;
}

/***************************************************************************//**
 *  Main thread control run by each working thread
 **/
void *
isched_parallel_section(isched_thread_t *ctx)
{
    isched_t *isched = ctx->global_ctx;
    int action;

    /* Wait for everyone but the master thread to enter this section */
    isched_barrier_wait( &(isched->barrier) );

    while(1) {
        pthread_mutex_lock( &(isched->statuslock) );
        while ((action = isched->status) == ISCHED_ACT_STAND_BY)
            pthread_cond_wait( &(isched->statuscond), &(isched->statuslock) );
        pthread_mutex_unlock( &(isched->statuslock) );
        isched_barrier_wait( &(isched->barrier) );

        switch (action) {
            case ISCHED_ACT_PARALLEL:
                isched->pfunc( ctx, isched->pargs );
                break;
            case ISCHED_ACT_FINALIZE:
                return isched_thread_destroy( ctx );
            default:
                fprintf(stderr, "isched_parallel_section: undefined action\n");
                return NULL;
        }
        isched_barrier_wait(&(isched->barrier) );
    }

    isched_thread_destroy( ctx );
    return NULL;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_isched
 *
 * @brief Initialize the stuctures of the internal thread scheduler.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = PASTIX_NUM_THREADS if it is set, the
 *          system number of core otherwise.
 *
 * @param[in] coresbind
 *          Array to specify where to bind each thread.
 *          Each thread i is binded to coresbind[hwloc(i)] if hwloc is
 *          provided, or to coresbind[i] otherwise.
 *          If coresbind = NULL, the identity function distribution is used.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS successful exit
 *
 ******************************************************************************/
isched_t *ischedInit(int cores, const int *coresbind)
{
    __isched_init_t *initdata;
    isched_t *isched;
    int core;

    /* Init the isched structure */
    MALLOC_INTERN(isched, 1, isched_t);
    if (isched == NULL) {
        fprintf(stderr, "ischedInit: isched allocation failed\n");
        return NULL;
    }

    pthread_mutex_init(&(isched->statuslock), NULL);
    pthread_cond_init( &(isched->statuscond), NULL);
    isched->status = ISCHED_ACT_STAND_BY;

    isched->pfunc = NULL;
    isched->pargs = NULL;

    /* Init number of cores and topology */
    isched_topo_init();

    /* Set number of cores */
    if ( cores < 1 ) {
        isched->world_size = pastix_getenv_get_value_int("PASTIX_NUM_THREADS", -1);
        if ( isched->world_size == -1 ) {
            isched->world_size = isched_topo_world_size();
            fprintf(stderr, "ischedInit: Could not find the number of cores: the thread number is set to %d\n", isched->world_size);
        }
    }
    else
        isched->world_size = cores;

    if (isched->world_size <= 0) {
        fprintf(stderr, "ischedInit: failed to get system size, set to 1\n");
        isched->world_size = 1;
    }

    /* Initialize barrier */
    isched_barrier_init( &(isched->barrier), NULL, isched->world_size );

    /* If the given coresbind is NULL, try to get one from the environment */
    if (coresbind == NULL) {
        //coresbind = pastix_getenv_get_array_int("PASTIX_THREADS_AFFINITY", NULL );
    }

    MALLOC_INTERN( initdata, isched->world_size, __isched_init_t );

    /* Initialize default thread attributes */
    if ( isched->world_size > 1 ) {
        /*  Launch threads */
        MALLOC_INTERN(isched->tids, isched->world_size, pthread_t);

        for (core = 1; core < isched->world_size; core++) {
            initdata[core].global_ctx = isched;
            initdata[core].rank       = core;
            initdata[core].bindto     = (coresbind == NULL) ? core : coresbind[core];

            pthread_create(
                &isched->tids[core],
                NULL,
                isched_thread_init,
                (void*)(initdata + core));
        }
    }
    else {
        isched->tids = NULL;
    }

    initdata[0].global_ctx = isched;
    initdata[0].rank       = 0;
    initdata[0].bindto     = (coresbind == NULL) ? -1 : coresbind[0];
    isched->master = (isched_thread_t*)isched_thread_init( initdata );

    /* Wait for the other threads to finish their initialization */
    isched_barrier_wait( &(isched->barrier) );

    /* Free temporary informations */
    memFree_null(initdata);

    return isched;
}

/**
 *****************************************************************************
 *
 * @ingroup pastix_isched
 *
 *  ischedFinalize - Finalize the structures associated to the internal threads
 *  scheduler.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS successful exit
 *
 ******************************************************************************/
int ischedFinalize(isched_t *isched)
{
    int core;
    int status;
    void *exitcodep;

    /* Make sure the threads exit the parallel function */
    pthread_mutex_lock(&isched->statuslock);
    isched->status = ISCHED_ACT_FINALIZE;
    pthread_mutex_unlock(&isched->statuslock);
    pthread_cond_broadcast(&isched->statuscond);
    isched_barrier_wait(&(isched->barrier));
    isched->status = ISCHED_ACT_STAND_BY;

    // Join threads
    for (core = 1; core < isched->world_size; core++) {
        status = pthread_join(isched->tids[core], &exitcodep);
        if (status != 0) {
            fprintf(stderr, "ischedFinalize: pthread_join() failed\n");
            return status;
        }
    }

    /* Unbind main thread */
    isched_thread_destroy( isched->master );

    isched_barrier_destroy(&(isched->barrier));

    /* Destroy topology */
    isched_topo_destroy();

    memFree_null(isched->tids);
    memFree_null(isched);

    return PASTIX_SUCCESS;
}

