/*
--------------------------------------------------------------------
By Bob Jenkins, September 1996.  _gp_recycle.c
You may use this code in any way you wish, and it is free.  No warranty.

This manages memory for commonly-allocated structures.
It allocates RESTART to REMAX items at a time.
Timings have shown that, if malloc is used for every new structure,
  malloc will consume about 90% of the time in a program.  This
  module cuts down the number of mallocs by an order of magnitude.
This also decreases memory fragmentation, and freeing structures
  only requires freeing the root.
--------------------------------------------------------------------
*/

#ifndef STANDARD
# include "standard.h"
#endif
#ifndef RECYCLE
# include "recycle.h"
#endif

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif // _WIN32

reroot *_gp_remkroot( size_t  size)
{
   reroot *r = (reroot *)_gp_remalloc(sizeof(reroot), "recycle.c, root");
   r->list = (G_recycle *)0;
   r->trash = (G_recycle *)0;
   r->size = align(size);
   r->logsize = RESTART;
   r->numleft = 0;
   return r;
}

void  _gp_refree( struct reroot *r)
{
   G_recycle *temp;
   if ((temp = r->list)) while (r->list)
   {
      temp = r->list->next;
      free((char *)r->list);
      r->list = temp;
   }
   free((char *)r);
   return;
}

/* to be called from the macro _gp_renew only */
char  *_gp_renewx( struct reroot *r)
{
   G_recycle *temp;
   if (r->trash)
   {  /* pull a node off the trash heap */
      temp = r->trash;
      r->trash = temp->next;
      /*BOB: (void)memset((void *)temp, 0, r->size);*/
      memset(temp, 0, r->size);
   }
   else
   {  /* allocate a new block of G_nodes */
      r->numleft =(word)(r->size*(size_t)((ub4)1<<r->logsize));
      if (r->numleft < REMAX) ++r->logsize;
      temp = (G_recycle *)_gp_remalloc(sizeof(G_recycle) + r->numleft, 
				 "recycle.c, data");
      temp->next = r->list;
      r->list = temp;
      r->numleft-=(word)(r->size);
      temp = (G_recycle *)((char *)(r->list+1)+r->numleft);
   }
   return (char *)temp;
}

char   *_gp_remalloc(
    size_t  len,
    char   *purpose
                 )
{
  char *x = (char *)malloc(len);
  if (!x)
  {
    fprintf(stderr, "GiDPost: malloc of %u failed for %s\n", 
	    ( unsigned int)len, purpose);
    exit(SUCCESS);
  }
  return x;
}

#ifndef _WIN32
#pragma GCC diagnostic pop
#endif // _WIN32
