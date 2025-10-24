/*
--------------------------------------------------------------------
By Bob Jenkins, 1996.  hashtab.c.  Public Domain.

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

#ifndef STANDARD
#include "standard.h"
#endif
#ifndef LOOKUPA
#include "lookupa.h"
#endif
#ifndef HASHTAB
#include "hashtab.h"
#endif
#ifndef RECYCLE
#include "recycle.h"
#endif

#ifdef HSANITY
/* sanity check -- make sure ipos, apos, and count make sense */
static void  hsanity( _gp_htab *t)
{
  ub4    i, end, counter;
  _gp_hitem *h;

  /* test that apos makes sense */
  end = (ub4)1<<(t->logsize);
  if (end < t->apos)
    printf("error:  end %ld  apos %ld\n", end, t->apos);

  /* test that ipos is in bucket apos */
  if (t->ipos)
  {
    for (h=t->table[t->apos];  h && h != t->ipos;  h = h->next)
      ;
    if (h != t->ipos)
      printf("error:ipos not in apos, apos is %ld\n", t->apos);
  }

  /* test that t->count is the number of elements in the table */
  counter=0;
  for (counter=0, i=0;  i<end;  ++i)
    for (h=t->table[i];  h;  h=h->next)
      ++counter;
  if (counter != t->count)
    printf("error: counter %ld  t->count %ld\n", counter, t->count);
}
#endif

/*
 * hgrow - Double the size of a hash table.
 * Allocate a new, 2x bigger array,
 * move everything from the old array to the new array,
 * then free the old array.
 */
static void hgrow( G_htab  *t)
    /* table */
{
  register ub4     newsize = (ub4)1<<(++t->logsize);
  register ub4     newmask = newsize-1;
  register ub4     i;
  register G_hitem **oldtab = t->table;
  register G_hitem **newtab = (G_hitem **)malloc(newsize*sizeof(G_hitem *));

  /* make sure newtab is cleared */
  for (i=0; i<newsize; ++i) newtab[i] = (G_hitem *)0;
  t->table = newtab;
  t->mask = newmask;

  /* Walk through old table putting entries in new table */
  for (i=newsize>>1; i--;)
  {
    register G_hitem *this, *that, **newplace;
    for (this = oldtab[i]; this;)
    {
      that = this;
      this = this->next;
      newplace = &newtab[(that->hval & newmask)];
      that->next = *newplace;
      *newplace = that;
    }
  }

  /* position the hash table on some existing item */
  _gp_hfirst(t);

  /* free the old array */
  free((char *)oldtab);

}

/* _gp_hcreate - create a hash table initially of size power(2,logsize) */
G_htab *_gp_hcreate( word  logsize)    /* log base 2 of the size of the hash table */
{
  ub4 i,len;
  G_htab *t = (G_htab *)malloc(sizeof(G_htab));

  len = ((ub4)1<<logsize);
  t->table = (G_hitem **)malloc(sizeof(G_hitem *)*(ub4)len);
  for (i=0; i<len; ++i) t->table[i] = (G_hitem *)0;
  t->logsize = logsize;
  t->mask = len-1;
  t->count = 0;
  t->apos = (ub4)0;
  t->ipos = (G_hitem *)0;
  t->space = _gp_remkroot(sizeof(G_hitem));
  t->bcount = 0;
  return t;
}

/* _gp_hdestroy - destroy the hash table and free all its memory */
void _gp_hdestroy( G_htab  *t)    /* the table */
{
  /* _gp_hitem *h; */
  _gp_refree(t->space);
  free((char *)t->table);
  free((char *)t);
}

/* _gp_hcount() is a macro, see hashtab.h */
/* _gp_hkey() is a macro, see hashtab.h */
/* _gp_hkeyl() is a macro, see hashtab.h */
/* _gp_hstuff() is a macro, see hashtab.h */

/* _gp_hfind - find an item with a given key in a hash table */
word   _gp_hfind( 
    G_htab  *t,     /* table */
    ub1   *key,   /* key to find */
    ub4    keyl  /* key length */
              )
{
  G_hitem *h;
  ub4    x = _gp_lookup(key,keyl,0);
  ub4    y;
  for (h = t->table[y=(x&t->mask)]; h; h = h->next)
  {
    if ((x == h->hval) && 
        (keyl == h->keyl) && 
        !memcmp(key, h->key, keyl))
    {
      t->apos = y;
      t->ipos = h;
      return TRUE;
    }
  }
  return FALSE;
}

/*
 * _gp_hadd - add an item to a hash table.
 * return FALSE if the key is already there, otherwise TRUE.
 */
word _gp_hadd( 
    G_htab  *t,      /* table */
    ub1   *key,    /* key to add to hash table */
    ub4    keyl,   /* key length */
    void  *stuff  /* stuff to associate with this key */
           )
{
  register G_hitem  *h,**hp;
  register ub4     y, x = _gp_lookup(key,keyl,0);

  /* make sure the key is not already there */
  for (h = t->table[(y=(x&t->mask))]; h; h = h->next)
  {
    if ((x == h->hval) && 
        (keyl == h->keyl) && 
        !memcmp(key, h->key, keyl))
    {
      t->apos = y;
      t->ipos = h;
      return FALSE;
    }
  }

  /* find space for a new item */
  h = (G_hitem *)_gp_renew(t->space);

  /* make the hash table bigger if it is getting full */
  if (++t->count > (ub4)1<<(t->logsize))
  {
    hgrow(t);
    y = (x&t->mask);
  }

  /* add the new key to the table */
  h->key   = key;
  h->keyl  = keyl;
  h->stuff = stuff;
  h->hval  = x;
  hp = &t->table[y];
  h->next = *hp;
  *hp = h;
  t->ipos = h;
  t->apos = y;

#ifdef HSANITY
  hsanity(t);
#endif  /* HSANITY */

  return TRUE;
}

/* _gp_hdel - delete the item at the current position */
word  _gp_hdel( G_htab *t)      /* the hash table */
{
  G_hitem  *h;    /* item being deleted */
  G_hitem **ip;   /* a counter */

  /* check for item not existing */
  if (!(h = t->ipos)) return FALSE;

  /* remove item from its list */
  for (ip = &t->table[t->apos]; *ip != h; ip = &(*ip)->next)
    ;
  *ip = (*ip)->next;
  --(t->count);

  /* adjust position to something that exists */
  if (!(t->ipos = h->next)) _gp_hnbucket(t);

  /* _gp_recycle the deleted _gp_hitem node */
  _gp_redel(t->space, h);

#ifdef HSANITY
  hsanity(t);
#endif  /* HSANITY */

  return TRUE;
}

/* _gp_hfirst - position on the first element in the table */
word _gp_hfirst( G_htab  *t)    /* the hash table */
{
  t->apos = (ub4)(t->mask);
  (void)_gp_hnbucket(t);
  return (t->ipos != (G_hitem *)0);
}

/* _gp_hnext() is a macro, see hashtab.h */

/*
 * _gp_hnbucket - Move position to the first item in the next bucket.
 * Return TRUE if we did not wrap around to the beginning of the table
 */
word _gp_hnbucket( G_htab *t)
{
  ub4  oldapos = t->apos;
  ub4  end = (ub4)1<<(t->logsize);
  ub4  i;

  /* see if the element can be found without wrapping around */
  for (i=oldapos+1; i<end; ++i)
  {
    if (t->table[i&t->mask])
    {
      t->apos = i;
      t->ipos = t->table[i];
      return TRUE;
    }
  }

  /* must have to wrap around to find the last element */
  for (i=0; i<=oldapos; ++i)
  {
    if (t->table[i])
    {
      t->apos = i;
      t->ipos = t->table[i];
      return FALSE;
    }
  }

  return FALSE;
}

void _gp_hstat( G_htab  *t)
{
  ub4     i,j;
  double  total = 0.0;
  G_hitem  *h;
  G_hitem  *walk, *walk2, *stat = (G_hitem *)0;

  /* in stat, keyl will store length of list, hval the number of buckets */
  for (i=0; i<=t->mask; ++i)
  {
    for (h=t->table[i], j=0; h; ++j, h=h->next)
      ;
    for (walk=stat; walk && (walk->keyl != j); walk=walk->next)
      ;
    if (walk)
    {
      ++(walk->hval);
    }
    else
    {
      walk = (G_hitem *)_gp_renew(t->space);
      walk->keyl = j;
      walk->hval = 1;
      if (!stat || stat->keyl > j) {walk->next=stat; stat=walk;}
      else
      {
        for (walk2=stat;
             walk2->next && (walk2->next->keyl<j);
             walk2=walk2->next)
          ;
        walk->next = walk2->next;
        walk2->next = walk;
      }
    }
  }

  /* figure out average list length for existing elements */
  for (walk=stat; walk; walk=walk->next)
  {
    total+=(double)walk->hval*(double)walk->keyl*(double)walk->keyl;
  }
  if (t->count) total /= (double)t->count;
  else          total  = (double)0;

  /* print statistics */
  printf("\n");
  for (walk=stat; walk; walk=walk->next)
  {
    printf("items %ld:  %ld buckets\n", walk->keyl, walk->hval);
  }
  printf("\nbuckets: %ld  items: %ld  existing: %g\n\n",
         ((ub4)1<<t->logsize), t->count, total);

  /* clean up */
  while (stat)
  {
    walk = stat->next;
    _gp_redel(t->space, stat);
    stat = walk;
  }
}

#if 0

// this example is for linux/macos
// as they use different versions of _gp_hcreate/_gp_hdestroy as the ones above implemented
// the above ones are similar to hcreate_r / hsearch_r / hdestroy_r
#include <search.h>

void pp() {
  static char *data[] = { "alpha",  "bravo", "charlie", "delta",  "echo",     "foxtrot", "golf",   "hotel",  "india",
                          "juliet", "kilo",  "lima",    "mike",   "november", "oscar",   "papa",   "quebec", "romeo",
                          "sierra", "tango", "uniform", "victor", "whisky",   "x-ray",   "yankee", "zulu" };

  ENTRY e, *ep;
  int i;

  _gp_hcreate( 30 );

  for ( i = 0; i < 24; i++ ) {
    e.key = data[ i ];
    /* data is just an integer, instead of a
       pointer to something */
    e.data = ( void * )i;
    ep = hsearch( e, ENTER );
    /* there should be no failures */
    if ( ep == NULL ) {
      fprintf( stderr, "entry failed\n" );
      exit( EXIT_FAILURE );
    }
  }

  for ( i = 22; i < 26; i++ ) {
    /* print two entries from the table, and
       show that two are not in the table */
    e.key = data[ i ];
    ep = hsearch( e, FIND );
    printf( "%9.9s -> %9.9s:%d\n", e.key, ep ? ep->key : "NULL", ep ? ( int )( ep->data ) : 0 );
  }
  _gp_hdestroy();
  exit( EXIT_SUCCESS );
}
#endif // if 0
