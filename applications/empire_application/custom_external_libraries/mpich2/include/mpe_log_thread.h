/*
   (C) 2007 by Argonne National Laboratory.
       See COPYRIGHT in top-level directory.
*/

#define HAVE_PTHREAD_IN_MPI 1

#define MPE_ThreadID_t   int

typedef struct {
    MPE_ThreadID_t   thdID;
    int              is_log_on;
} MPE_ThreadStm_t;

#if defined( HAVE_PTHREAD_IN_MPI )
#include <pthread.h>
#include "mpe_callstack.h"

#define MPE_LOG_THREAD_PRINTSTACK() \
        do { \
            MPE_CallStack_t  cstk; \
            MPE_CallStack_init( &cstk ); \
            MPE_CallStack_fancyprint( &cstk, 2, \
                                      "\t", 1, MPE_CALLSTACK_UNLIMITED ); \
        } while (0)

/* MPE coarse-grained lock support mechanism */
extern pthread_mutex_t  MPE_Thread_mutex;
extern pthread_key_t    MPE_ThreadStm_key;
extern MPE_ThreadID_t   MPE_Thread_count;

#define MPE_LOG_THREADSTM_DECL \
    MPE_ThreadStm_t *thdstm = NULL;

#define MPE_LOG_THREAD_LOCK  \
    if ( pthread_mutex_lock( &MPE_Thread_mutex ) != 0 ) { \
        perror( "MPE_LOG_THREAD: pthread_mutex_lock() fails!\n" ); \
        MPE_LOG_THREAD_PRINTSTACK(); \
        pthread_exit( NULL ); \
    }

#define MPE_LOG_THREAD_UNLOCK \
    if ( pthread_mutex_unlock( &MPE_Thread_mutex ) != 0 ) { \
        perror( "MPE_LOG_THREAD: pthread_mutex_unlock() fails!\n" ); \
        MPE_LOG_THREAD_PRINTSTACK(); \
        pthread_exit( NULL ); \
    }

/* MPE_Thread_count is global variable and should be protected by a mutex */
#define MPE_LOG_THREADSTM_SET \
    MPE_LOG_THREAD_LOCK \
    thdstm  = (MPE_ThreadStm_t *) malloc( sizeof(MPE_ThreadStm_t) ); \
    thdstm->thdID     = MPE_Thread_count ; \
    thdstm->is_log_on = 1; \
    if ( pthread_setspecific( MPE_ThreadStm_key, thdstm ) != 0 ) { \
        perror( "MPE_LOG_THREAD: pthread_setspecific() fails!\n" ); \
        MPE_LOG_THREAD_PRINTSTACK(); \
        pthread_exit( NULL ); \
    } \
    MPE_Thread_count++; \
    MPE_LOG_THREAD_UNLOCK

#define MPE_LOG_THREADSTM_GET \
    thdstm = (MPE_ThreadStm_t *) pthread_getspecific( MPE_ThreadStm_key ); \
    if ( thdstm == NULL ) { \
        MPE_LOG_THREADSTM_SET \
    }

#define MPE_LOG_THREAD_INIT      MPE_Log_thread_init();
#define MPE_LOG_THREAD_FINALIZE  MPE_Log_thread_sync( MPE_Thread_count );

void MPE_Log_thread_init( void );

#define THREADID        thdstm->thdID
#define IS_MPELOG_ON    thdstm->is_log_on

#else

/* When there is no thread support, thdID has to be a constant, i.e. 0 */
#define MPE_LOG_THREADSTM_DECL

/* Define everything else to NULL */
#define MPE_LOG_THREAD_LOCK
#define MPE_LOG_THREAD_UNLOCK
#define MPE_LOG_THREADSTM_SET
#define MPE_LOG_THREADSTM_GET
#define MPE_LOG_THREAD_INIT
#define MPE_LOG_THREAD_FINALIZE

#define THREADID        0
#define IS_MPELOG_ON    is_mpelog_on

#endif
