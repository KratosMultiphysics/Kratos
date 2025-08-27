/*
--------------------------------------------------------------------
By Bob Jenkins, September 1996.  _gp_recycle.h
You may use this code in any way you wish, and it is free.  No warranty.

This manages memory for commonly-allocated structures.
It allocates RESTART to REMAX items at a time.
Timings have shown that, if malloc is used for every new structure,
  malloc will consume about 90% of the time in a program.  This
  module cuts down the number of mallocs by an order of magnitude.
This also decreases memory fragmentation, and freeing all structures
  only requires freeing the root.
--------------------------------------------------------------------
*/

#ifndef _RECYCLE_H_
#define _RECYCLE_H_

#include "standard.h"

#define RESTART    0
#define REMAX      32000

struct _gp_recycle
{
   struct _gp_recycle *next;
};
typedef  struct _gp_recycle  G_recycle;

struct reroot
{
   struct _gp_recycle *list;     /* list of malloced blocks */
   struct _gp_recycle *trash;    /* list of deleted items */
   size_t          size;     /* size of an item */
   size_t          logsize;  /* log_2 of number of items in a block */
   word            numleft;  /* number of bytes left in this block */
};
typedef  struct reroot  reroot;

/* make a new recycling root */
reroot  *_gp_remkroot( size_t mysize );

/* free a recycling root and all the items it has made */
void     _gp_refree( struct reroot *r );

/* get a new (cleared) item from the root */
#define _gp_renew(r) ((r)->numleft ? \
   (((char *)((r)->list+1))+((r)->numleft-=(word)(r)->size)) : _gp_renewx(r))

char    *_gp_renewx( struct reroot *r );

/* delete an item; let the root _gp_recycle it */
/* void     _gp_redel(/o_ struct reroot *r, struct _gp_recycle *item _o/); */
#define _gp_redel(root,item) { \
   ((G_recycle *)item)->next=(root)->trash; \
   (root)->trash=(G_recycle *)(item); \
}

/* malloc, but complain to stderr and exit program if no joy */
/* use plain free() to free memory allocated by _gp_remalloc() */
char    *_gp_remalloc( size_t len, char *purpose );

#endif  /* RECYCLE */
