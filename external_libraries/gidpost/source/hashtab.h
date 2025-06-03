/*
--------------------------------------------------------------------
By Bob Jenkins, 1996.  hash.h.  Public Domain.

This implements a hash table.
* Keys are unique.  Adding an item fails if the key is already there.
* Keys and items are pointed at, not copied.  If you change the value
  of the key after it is inserted then hfind will not be able to find it.
* The hash table maintains a position that can be set and queried.
* The table length doubles dynamically and never shrinks.  The insert
  that causes table doubling may take a long time.
* The table length splits when the table length equals the number of items
  Comparisons usually take 7 instructions.
  Computing a hash value takes 35+6n instructions for an n-byte key.

  hcreate  - create a hash table
  hdestroy - destroy a hash table
   hcount  - The number of items in the hash table
   hkey    - key at the current position
   hkeyl   - key length at the current position
   hstuff  - stuff at the current position
  hfind    - find an item in the table
   hadd    - insert an item into the table
   hdel    - delete an item from the table
  hstat    - print statistics about the table
   hfirst  - position at the first item in the table
   hnext   - move the position to the next item in the table
--------------------------------------------------------------------
*/

#ifndef _HASHTAB_H_
#define _HASHTAB_H_

#include "standard.h"

/* PRIVATE TYPES AND DEFINITIONS */

struct hitem
{
  ub1          *key;      /* key that is hashed */
  ub4           keyl;     /* length of key */
  void         *stuff;    /* stuff stored in this hitem */
  ub4           hval;     /* hash value */
  struct hitem *next;     /* next hitem in list */
};
typedef  struct hitem  hitem;


struct htab
{
  struct hitem **table;   /* hash table, array of size 2^logsize */
  word           logsize; /* log of size of table */
  size_t         mask;    /* (hashval & mask) is position in table */
  ub4            count;   /* how many items in this hash table so far? */
  ub4            apos;    /* position in the array */
  struct hitem  *ipos;    /* current item in the array */
  struct reroot *space;   /* space for the hitems */
  ub4            bcount;  /* # hitems useable in current block */
};
typedef  struct htab  htab;





/* PUBLIC FUNCTIONS */

/* hcreate - create a hash table
   ARGUMENTS:
     logsize - 1<<logsize will be the initial table length
   RETURNS:
     the new table
 */
htab *hcreate( word logsize);


/* hdestroy - destroy a hash table
   ARGUMENTS:
     t - the hash table to be destroyed.  Note that the items and keys
         will not be freed, the user created them and must destroy
         them himself.
   RETURNS:
     nothing
 */
void  hdestroy( htab *t );


/* hcount, hkey, hkeyl, hstuff
     ARGUMENTS:
     t - the hash table
   RETURNS:
     hcount - (ub4)    The number of items in the hash table
     hkey   - (ub1 *)  key for the current item
     hkeyl  - (ub4)    key length for the current item
     hstuff - (void *) stuff for the current item
   NOTE:
     The current position always has an item as long as there
       are items in the table, so hexist can be used to test if the
       table is empty.
     hkey, hkeyl, and hstuff will crash if hcount returns 0
 */
#define hcount(t) ((t)->count)
#define hkey(t)   ((t)->ipos->key)
#define hkeyl(t)  ((t)->ipos->keyl)
#define hstuff(t) ((t)->ipos->stuff)



/* hfind - move the current position to a given key
   ARGUMENTS:
     t    - the hash table
     key  - the key to look for
     keyl - length of the key
   RETURNS:
     TRUE if the item exists, FALSE if it does not.
     If the item exists, moves the current position to that item.
 */
word  hfind( htab *t, ub1 *key, ub4 keyl );


/* hadd - add a new item to the hash table
          change the position to point at the item with the key
   ARGUMENTS:
     t     - the hash table
     key   - the key to look for
     keyl  - length of the key
     stuff - other stuff to be stored in this item
   RETURNS:
     FALSE if the operation fails (because that key is already there).
 */
word  hadd( htab *t, ub1 *key, ub4 keyl, void *stuff );


/* hdel - delete the item at the current position
          change the position to the following item
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if there is no current item (meaning the table is empty)
  NOTE:
    This frees the item, but not the key or stuff stored in the item.
    If you want these then deal with them first.  For example:
      if (hfind(tab, key, keyl))
      {
        free(hkey(tab));
        free(hstuff(tab));
        hdel(tab);
      }
 */
word  hdel( htab *t );


/* hfirst - move position to the first item in the table
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if there is no current item (meaning the table is empty)
  NOTE:
 */
word hfirst( htab *t );


/* hnext - move position to the next item in the table
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if the position wraps around to the beginning of the table
  NOTE:
    To see every item in the table, do
      if (hfirst(t)) do
      {
        key   = hkey(t);
        stuff = hstuff(t);
      }
      while (hnext(t));
 */
/* word hnext(/o_ htab *t _o/); */
#define hnext(t) \
  ((!(t)->ipos) ? FALSE :  \
   ((t)->ipos=(t)->ipos->next) ? TRUE : hnbucket(t))

/* hnbucket - PRIVATE - move to first item in the next nonempty bucket
  ARGUMENTS:
    t    - the hash table
  RETURNS:
    FALSE if the position wraps around to the beginning of the table
  NOTE:
    This is private to hashtab; do not use it externally.
 */
word hnbucket( htab *t );


/* hstat - print statistics about the hash table
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
void hstat( htab *t );

#endif   /* HASHTAB */
