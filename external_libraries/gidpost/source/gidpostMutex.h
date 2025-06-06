#pragma once

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#ifndef GP_MUTEX_ENV
#define GP_MUTEX_ENV extern
#endif

#ifdef WIN32
// defined in  gidpostHash.c
GP_MUTEX_ENV CRITICAL_SECTION  cs;

#define _INIT_MUTEX_    InitializeCriticalSection( &cs)
#define _DESTROY_MUTEX_ DeleteCriticalSection( &cs);
#define _LOCK_          EnterCriticalSection(&cs)
#define _UNLOCK_        LeaveCriticalSection(&cs)
#else
// defined in  gidpostHash.c
GP_MUTEX_ENV pthread_mutex_t lock;

#define _INIT_MUTEX_    pthread_mutex_init(&lock, NULL);
#define _DESTROY_MUTEX_ pthread_mutex_destroy(&lock);
#define _LOCK_          pthread_mutex_lock(&lock);
#define _UNLOCK_        pthread_mutex_unlock(&lock);
#endif
