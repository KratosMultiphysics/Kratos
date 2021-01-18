/* gidpost 2.0 */
/* -*- mode: c++ -*-
 *
 *  gidpostHash.c --
 *
 *    This file implement the interface to the pool of file handlers
 *    which is a mapping from unsigned integers (INT*4) to
 *    CPostFile*. The external API will work with unsigned integers as
 *    file handlers.
 */

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#include <string.h>

#include "gidpostHash.h"
#include "hashtab.h"
#include "recycle.h"

#include <assert.h>

static GiD_FILE lastKey = 0;
static htab *hashTable = NULL;
static reroot *fd_pool = NULL;

#ifdef WIN32
CRITICAL_SECTION  cs;

#define _INIT_MUTEX_    InitializeCriticalSection( &cs)
#define _DESTROY_MUTEX_ DeleteCriticalSection( &cs);
#define _LOCK_          EnterCriticalSection(&cs)
#define _UNLOCK_        LeaveCriticalSection(&cs)
#else
pthread_mutex_t lock;

#define _INIT_MUTEX_    pthread_mutex_init(&lock, NULL);
#define _DESTROY_MUTEX_ pthread_mutex_destroy(&lock);
#define _LOCK_          pthread_mutex_lock(&lock);
#define _UNLOCK_        pthread_mutex_unlock(&lock);
#endif

int GiD_HashInit()
{
  //printf("ENTER: GiD_HashInit\n");
  if (fd_pool)
    {
    //printf("LEAVE: GiD_HashInit Already inited\n"); 
    return 0;
    }
  _INIT_MUTEX_;

  fd_pool = remkroot(sizeof(GiD_FILE));
           
  if (!hashTable) {
    hashTable = gid_hcreate(16);
    assert(hashTable);
    //printf("LEAVE: GiD_HashInit %p\n", hashTable);  
    return 0;
  }
  //printf("LEAVE: GiD_HashInit\n");  
  return -1;

}

int GiD_HashDone()
{
  _LOCK_;
  
  htabdestroy(hashTable);
  refree(fd_pool);

  hashTable = NULL;
  fd_pool=NULL;
  
  _UNLOCK_;
  _DESTROY_MUTEX_;
  return 0;
}

GiD_FILE  GiD_HashAdd(void *data)
{
  GiD_FILE fd;
  ub1 *key;
  
  _LOCK_;
  
  /* must check overflow and possibly reuse removed keys */
  fd = ++lastKey;
  /* create a new key */
  key = (ub1*)renew(fd_pool);
  /* with fd as value */
  memcpy(key, &fd, sizeof(fd));
  /* insert data into hash table under key */
  hadd(hashTable, key, sizeof(GiD_FILE), data);

  _UNLOCK_;
  
  return fd;
}

void *GiD_HashFind  (GiD_FILE fd)
{
  void *data = NULL;
  
  if (!hashTable) {
    return NULL;
  }
  _LOCK_;

  if (hfind(hashTable, (ub1*)&fd, sizeof(GiD_FILE))) {
    data = hstuff(hashTable);
  }

  _UNLOCK_;
  return data;
}

void     GiD_HashRemove(GiD_FILE fd)
{
  if (!hashTable) {
    return;
  }

  _LOCK_;
  
  if (hfind(hashTable, (ub1*)&fd, sizeof(GiD_FILE))) {
    /* free the key from the pool */
    redel(fd_pool, hkey(hashTable));
    /* free the hash entry */
    hdel(hashTable);
  }

  _UNLOCK_;
}
