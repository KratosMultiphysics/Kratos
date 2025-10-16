/*
--------------------------------------------------------------------
By Bob Jenkins, 1996.  hash.h.  Public Domain.

This implements a hash table.
* Keys are unique.  Adding an item fails if the key is already there.
* Keys and items are pointed at, not copied.  If you change the value
  of the key after it is inserted then _gp_hfind will not be able to find it.
* The hash table maintains a position that can be set and queried.
* The table length doubles dynamically and never shrinks.  The insert
  that causes table doubling may take a long time.
* The table length splits when the table length equals the number of items
  Comparisons usually take 7 instructions.
  Computing a hash value takes 35+6n instructions for an n-byte key.

  _gp_hcreate  - create a hash table
  _gp_hdestroy - destroy a hash table
   _gp_hcount  - The number of items in the hash table
   _gp_hkey    - key at the current position
   _gp_hkeyl   - key length at the current position
   _gp_hstuff  - stuff at the current position
  _gp_hfind    - find an item in the table
   _gp_hadd    - insert an item into the table
   _gp_hdel    - delete an item from the table
  _gp_hstat    - print statistics about the table
   _gp_hfirst  - position at the first item in the table
   _gp_hnext   - move the position to the next item in the table
--------------------------------------------------------------------
*/

#ifndef _HASHTAB_H_
#define _HASHTAB_H_

#include "standard.h"

/* PRIVATE TYPES AND DEFINITIONS */

struct _gp_hitem
{
  ub1          *key;      /* key that is hashed */
  ub4           keyl;     /* length of key */
  void         *stuff;    /* stuff stored in this _gp_hitem */
  ub4           hval;     /* hash value */
  struct _gp_hitem *next;     /* next _gp_hitem in list */
};
typedef  struct _gp_hitem  G_hitem;


struct _gp_htab
{
  struct _gp_hitem **table;   /* hash table, array of size 2^logsize */
  word           logsize; /* log of size of table */
  size_t         mask;    /* (hashval & mask) is position in table */
  ub4            count;   /* how many items in this hash table so far? */
  ub4            apos;    /* position in the array */
  struct _gp_hitem  *ipos;    /* current item in the array */
  struct reroot *space;   /* space for the hitems */
  ub4            bcount;  /* # hitems useable in current block */
};
typedef  struct _gp_htab  G_htab;


/* PUBLIC FUNCTIONS */

/* _gp_hcreate - create a hash table
   ARGUMENTS:
     logsize - 1<<logsize will be the initial table length
   RETURNS:
     the new table
 */
G_htab *_gp_hcreate( word logsize);


/* _gp_hdestroy - destroy a hash table
   ARGUMENTS:
     t - the hash table to be destroyed.  Note that the items and keys
         will not be freed, the user created them and must destroy
         them himself.
   RETURNS:
     nothing
 */
void  _gp_hdestroy( G_htab *t );


/* _gp_hcount, _gp_hkey, _gp_hkeyl, _gp_hstuff
     ARGUMENTS:
     t - the hash table
   RETURNS:
     _gp_hcount - (ub4)    The number of items in the hash table
     _gp_hkey   - (ub1 *)  key for the current item
     _gp_hkeyl  - (ub4)    key length for the current item
     _gp_hstuff - (void *) stuff for the current item
   NOTE:
     The current position always has an item as long as there
       are items in the table, so hexist can be used to test if the
       table is empty.
     _gp_hkey, _gp_hkeyl, and _gp_hstuff will crash if _gp_hcount returns 0
 */
#define _gp_hcount(t) ((t)->count)
#define _gp_hkey(t)   ((t)->ipos->key)
#define _gp_hkeyl(t)  ((t)->ipos->keyl)
#define _gp_hstuff(t) ((t)->ipos->stuff)



/* _gp_hfind - move the current position to a given key
   ARGUMENTS:
     t    - the hash table
     key  - the key to look for
     keyl - length of the key
   RETURNS:
     TRUE if the item exists, FALSE if it does not.
     If the item exists, moves the current position to that item.
 */
word  _gp_hfind( G_htab *t, ub1 *key, ub4 keyl );


/* _gp_hadd - add a new item to the hash table
          change the position to point at the item with the key
   ARGUMENTS:
     t     - the hash table
     key   - the key to look for
     keyl  - length of the key
     stuff - other stuff to be stored in this item
   RETURNS:
     FALSE if the operation fails (because that key is already there).
 */
word  _gp_hadd( G_htab *t, ub1 *key, ub4 keyl, void *stuff );


/* _gp_hdel - delete the item at the current position
          change the position to the following item
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if there is no current item (meaning the table is empty)
  NOTE:
    This frees the item, but not the key or stuff stored in the item.
    If you want these then deal with them first.  For example:
      if (_gp_hfind(tab, key, keyl))
      {
        free(_gp_hkey(tab));
        free(_gp_hstuff(tab));
        _gp_hdel(tab);
      }
 */
word  _gp_hdel( G_htab *t );


/* _gp_hfirst - move position to the first item in the table
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if there is no current item (meaning the table is empty)
  NOTE:
 */
word _gp_hfirst( G_htab *t );


/* _gp_hnext - move position to the next item in the table
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if the position wraps around to the beginning of the table
  NOTE:
    To see every item in the table, do
      if (_gp_hfirst(t)) do
      {
        key   = _gp_hkey(t);
        stuff = _gp_hstuff(t);
      }
      while (_gp_hnext(t));
 */
/* word _gp_hnext(/o_ _gp_htab *t _o/); */
#define _gp_hnext(t) \
  ((!(t)->ipos) ? FALSE :  \
   ((t)->ipos=(t)->ipos->next) ? TRUE : _gp_hnbucket(t))

/* _gp_hnbucket - PRIVATE - move to first item in the next nonempty bucket
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if the position wraps around to the beginning of the table
  NOTE:
    This is private to hashtab; do not use it externally.
 */
word _gp_hnbucket( G_htab *t );


/* _gp_hstat - print statistics about the hash table
  ARGUMENTS:
    t    - the hash table
  NOTE:
    items <0>:  <#buckets with zero items> buckets
    items <1>:  <#buckets with 1 item> buckets
    ...
    buckets: #buckets  items: #items  existing: x
    ( x is the average length of the list when you look for an
      item that exists.  When the item does not exists, the average
      length is #items/#buckets. )

    If you put n items into n buckets, expect 1/(n!)e buckets to
    have n items.  That is, .3678 0, .3678 1, .1839 2, ...
    Also expect "existing" to be about 2.
 */
void _gp_hstat( G_htab *t );

#endif   /* HASHTAB */
