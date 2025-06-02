/*
--------------------------------------------------------------------
By Bob Jenkins, 1996.  hashtab.c.  Public Domain.

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
static void  hsanity(t)
htab *t;
{
  ub4    i, end, counter;
  hitem *h;

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
static void hgrow( htab  *t)
    /* table */
{
  register ub4     newsize = (ub4)1<<(++t->logsize);
  register ub4     newmask = newsize-1;
  register ub4     i;
  register hitem **oldtab = t->table;
  register hitem **newtab = (hitem **)malloc(newsize*sizeof(hitem *));

  /* make sure newtab is cleared */
  for (i=0; i<newsize; ++i) newtab[i] = (hitem *)0;
  t->table = newtab;
  t->mask = newmask;

  /* Walk through old table putting entries in new table */
  for (i=newsize>>1; i--;)
  {
    register hitem *this, *that, **newplace;
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
  hfirst(t);

  /* free the old array */
  free((char *)oldtab);

}

/* hcreate - create a hash table initially of size power(2,logsize) */
htab *hcreate(logsize)
word  logsize;    /* log base 2 of the size of the hash table */
{
  ub4 i,len;
  htab *t = (htab *)malloc(sizeof(htab));

  len = ((ub4)1<<logsize);
  t->table = (hitem **)malloc(sizeof(hitem *)*(ub4)len);
  for (i=0; i<len; ++i) t->table[i] = (hitem *)0;
  t->logsize = logsize;
  t->mask = len-1;
  t->count = 0;
  t->apos = (ub4)0;
  t->ipos = (hitem *)0;
  t->space = remkroot(sizeof(hitem));
  t->bcount = 0;
  return t;
}

/* hdestroy - destroy the hash table and free all its memory */
void hdestroy( t)
htab  *t;    /* the table */
{
  /* hitem *h; */
  refree(t->space);
  free((char *)t->table);
  free((char *)t);
}

/* hcount() is a macro, see hashtab.h */
/* hkey() is a macro, see hashtab.h */
/* hkeyl() is a macro, see hashtab.h */
/* hstuff() is a macro, see hashtab.h */

/* hfind - find an item with a given key in a hash table */
word   hfind( t, key, keyl )
htab  *t;     /* table */
ub1   *key;   /* key to find */
ub4    keyl;  /* key length */
{
  hitem *h;
  ub4    x = lookup(key,keyl,0);
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
 * hadd - add an item to a hash table.
 * return FALSE if the key is already there, otherwise TRUE.
 */
word hadd( t, key, keyl, stuff)
htab  *t;      /* table */
ub1   *key;    /* key to add to hash table */
ub4    keyl;   /* key length */
void  *stuff;  /* stuff to associate with this key */
{
  register hitem  *h,**hp;
  register ub4     y, x = lookup(key,keyl,0);

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
  h = (hitem *)renew(t->space);

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

/* hdel - delete the item at the current position */
word  hdel(t)
htab *t;      /* the hash table */
{
  hitem  *h;    /* item being deleted */
  hitem **ip;   /* a counter */

  /* check for item not existing */
  if (!(h = t->ipos)) return FALSE;

  /* remove item from its list */
  for (ip = &t->table[t->apos]; *ip != h; ip = &(*ip)->next)
    ;
  *ip = (*ip)->next;
  --(t->count);

  /* adjust position to something that exists */
  if (!(t->ipos = h->next)) hnbucket(t);

  /* recycle the deleted hitem node */
  redel(t->space, h);

#ifdef HSANITY
  hsanity(t);
#endif  /* HSANITY */

  return TRUE;
}

/* hfirst - position on the first element in the table */
word hfirst(t)
htab  *t;    /* the hash table */
{
  t->apos = (ub4)(t->mask);
  (void)hnbucket(t);
  return (t->ipos != (hitem *)0);
}

/* hnext() is a macro, see hashtab.h */

/*
 * hnbucket - Move position to the first item in the next bucket.
 * Return TRUE if we did not wrap around to the beginning of the table
 */
word hnbucket(t)
htab *t;
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

void hstat(t)
htab  *t;
{
  ub4     i,j;
  double  total = 0.0;
  hitem  *h;
  hitem  *walk, *walk2, *stat = (hitem *)0;

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
      walk = (hitem *)renew(t->space);
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
    redel(t->space, stat);
    stat = walk;
  }
}

#if 0

// this example is for linux/macos
// as they use different versions of hcreate/hdestroy as the ones above implemented
// the above ones are similar to hcreate_r / hsearch_r / hdestroy_r
#include <search.h>

void pp() {
  static char *data[] = { "alpha",  "bravo", "charlie", "delta",  "echo",     "foxtrot", "golf",   "hotel",  "india",
                          "juliet", "kilo",  "lima",    "mike",   "november", "oscar",   "papa",   "quebec", "romeo",
                          "sierra", "tango", "uniform", "victor", "whisky",   "x-ray",   "yankee", "zulu" };

  ENTRY e, *ep;
  int i;

  hcreate( 30 );

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
  hdestroy();
  exit( EXIT_SUCCESS );
}
#endif // if 0
