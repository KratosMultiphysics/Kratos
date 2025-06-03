/* gidpost */
/* -*- mode: c++ -*-
 *
 *  gidpostHash.c --
 *
 *    This file implement the interface to the pool of file handlers
 *    which is a mapping from unsigned integers (INT*4) to
 *    CPostFile*. The external API will work with unsigned integers as
 *    file handlers.
 */

#define GP_MUTEX_ENV 
#include "gidpostMutex.h"

#include <string.h>

#include "gidpostHash.h"

#include <assert.h>

#define OLD_HASH_TABLE

#ifdef OLD_HASH_TABLE

#include "hashtab.h"
#include "recycle.h"

static GiD_FILE lastKey = 0;
static htab *hashTable = NULL;
static reroot *fd_pool = NULL;

int GiD_HashInit()
{
  if ( fd_pool )
    return -1; // already initialized

  _INIT_MUTEX_;

  if(!fd_pool){
    fd_pool = remkroot(sizeof(GiD_FILE));
  }
           
  if (!hashTable) {
    hashTable = hcreate(16);
    assert(hashTable);
    return 0;
  }
  return -1;
}

int GiD_HashDone()
{
  _LOCK_;
  
  hdestroy(hashTable);
  refree(fd_pool);
  hashTable = NULL;
  fd_pool=NULL;
  
  _UNLOCK_;
  _DESTROY_MUTEX_;
  return 0;
}

GiD_FILE  GiD_HashAdd( CPostFile *file)
{
  void *data = ( void * )file;
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

CPostFile *GiD_HashFind  (GiD_FILE fd)
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
  return ( CPostFile *)data;
}

int     GiD_HashRemove(GiD_FILE fd)
{
  if (!hashTable) {
    return -1;
  }

  _LOCK_;
  
  if (hfind(hashTable, (ub1*)&fd, sizeof(GiD_FILE))) {
    /* free the key from the pool */
    redel(fd_pool, hkey(hashTable));
    /* free the hash entry */
    hdel(hashTable);
  }

  _UNLOCK_;

  return 0;
}

#else // OLD_HASH_TABLE

#include "gidpostInt.h"

static GiD_FILE G_last_file_id = 0;
static CPostFile **G_file_vector = NULL;
static size_t G_file_vector_size = 0; // number of entries in G_file_vector

// returns 0 if ok, -1 if not
static int _FileVectorResize( size_t new_size ) {
  CPostFile **tmp_vector = ( CPostFile ** )realloc( G_file_vector, new_size * sizeof( CPostFile *) );
  int fail = -1;
  if ( tmp_vector ) {
    long long int to_zero_size = new_size - G_file_vector_size;
    if ( to_zero_size > 0 ) {
      memset( &tmp_vector[ G_file_vector_size ], 0, ( size_t )to_zero_size * sizeof( CPostFile * ) );
    }
    G_file_vector = tmp_vector;
    G_file_vector_size = new_size;
  }
  return fail;
}

int GiD_HashInit( void ) {
  if ( G_file_vector != NULL ) {
    // already initialized
    return -1;
  }

  int fail = -1;
  _INIT_MUTEX_;
  _LOCK_;
  G_last_file_id = 0;
  size_t initial_size = 1024;
  fail = _FileVectorResize( initial_size );
  _UNLOCK_;

  return fail;
}

int GiD_HashDone( void ) {
  _LOCK_;
  G_last_file_id = 0;
  G_file_vector_size = 0;
  free( G_file_vector );
  _UNLOCK_;
  _DESTROY_MUTEX_;
  return 0;
}

static int _IsValidFileTable() {
  assert( G_file_vector != NULL );
  if ( G_file_vector == NULL ) {
    fprintf( stderr, "GiDPost: files open table not initialized.\n" );
    return 0;
  }
  assert( G_last_file_id < G_file_vector_size );
  if ( G_last_file_id == G_file_vector_size ) {
    fprintf( stderr, "GiDPost: maxium number of open files reached = %d.\n", ( int)G_file_vector_size );
    return 0;
  }
  return 1;
}

GiD_FILE GiD_HashAdd( CPostFile *file_pointer ) {
  int ok = _IsValidFileTable();
  if ( !ok )
    return 0;

  int new_file_id = 0;
  _LOCK_;
  // TODO: in theory i should look for NULL holes inside G_file_vector_size until G_last_file_id
  // but for simplicity i just add it at the end of the table.
  G_file_vector[ G_last_file_id ] = file_pointer;
  G_last_file_id++;
  new_file_id = G_last_file_id;
  _UNLOCK_;
  return new_file_id;
}

CPostFile *GiD_HashFind( GiD_FILE fd ) {
  int ok = _IsValidFileTable();
  if ( !ok )
    return NULL;
  assert( fd > 0 );
  if ( fd <= 0 ) {
    fprintf( stderr, "GiDPost: invalid file id = %d.\n", fd );
    return NULL;
  }
  return G_file_vector[ fd - 1 ];
}

int GiD_HashRemove( GiD_FILE fd ) {
  int ok = _IsValidFileTable();
  if ( !ok )
    return -1;
  assert( fd > 0 );
  if ( fd <= 0 ) {
    fprintf( stderr, "GiDPost: invalid file id = %d.\n", fd );
    return -1;
  }

  _LOCK_;
  G_file_vector[ fd - 1 ] = NULL;
  if ( fd == G_last_file_id ) {
    // if it was the last used id, then reuse it
    G_last_file_id--;
  }
  _UNLOCK_;

  return 0;
}

#endif // OLD_HASH_TABLE
