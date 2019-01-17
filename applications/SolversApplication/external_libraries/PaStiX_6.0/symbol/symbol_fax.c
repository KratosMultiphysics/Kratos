/**
 *
 * @file symbol_fax.c
 *
 * PaStiX fax symbol structure routines issued from Scotch esmumps library.
 * This is the generic block symbolic factorization routine.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @date 2018-07-16
 *
 *   Dates:
 *     Version 0.0 - from 22 jul 1998 to 29 sep 1998
 *     Version 0.1 - from 04 apr 1999 to 21 apr 1999
 *     Version 0.2 - from 08 may 2000 to 09 may 2000
 *     Version 1.0 - from 13 mar 2002 to 08 jun 2002
 *     Version 1.2 - from 23 aug 2002 to 23 aug 2002
 *     Version 2.0 - from 21 mar 2003 to 21 mar 2003
 *
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*
**  The defines and includes.
*/
/*
  Macros:

  SYMBOL_FAX_INCLUDED - Has to be defined if the code his included
                         in an external file function (see <faxi_graph.c>)
  SYMBOL_FAX          - Defined if SYMBOL_FAXI_INCLUDED is not...
                         But does not seem to be used...
*/
#define SYMBOL_FAX_ITERATOR_END }
#ifndef SYMBOL_FAX_INCLUDED                       /* If included from other file */
#define SYMBOL_FAX

#include "common.h"
#include "symbol.h"
#include "pastix/order.h"
#include "fax.h"
#include "symbol_fax.h"

/*
  Macro: SYMBOL_FAX_ITERATOR

  Loop for all adjacent edges, used in <symbolFaxi>.
  Must be defined in including file if SYMBOL_FAXI_INCLUDED is defined.

  Parameters:
    ngbdptr - Neighbour pointer.
    vertnum - Vertex index.
    vertend - Iterator.
*/
#define SYMBOL_FAX_ITERATOR(ngbdptr, vertnum, vertend)			\
  for (vertend  = ngbfrst ((ngbdptr), (vertnum));			\
       vertend >= baseval;						\
       vertend  = ngbnext (ngbdptr)) {
/*
  Macro: SYMBOL_FAX_VERTEX_DEGREE

  Computes the number of adjacent edges to a vertex.

  Parameters:
    ngbdptr - Neighbour pointer.
    vertnum - Vertex index.
*/
#define SYMBOL_FAX_VERTEX_DEGREE(ngbdptr, vertnum)	\
  (ngbdegr ((ngbdptr), (vertnum)))

/*
  Function: symbolFax

  Symbolic factorization routine.

  This routine computes the block symbolic
  factorization of the given matrix
  according to the given vertex ordering.

  Algorithm:

  The algorithm is implemented in a
  cache-friendly manner, by using a single
  dynamic array which grows along with the
  number of computed blocks. The array is
  decomposed in the following manner:

    - In a first phase, a hash table and a
      sort area are reserved at the end of
      the space of already computed blocks.
      The sort area is created far enough from
      the end of the array of already computed
      blocks such that if there are no contributing
      blocks all new blocks can be created without
      colliding with the sort area.

    - Then, in a second phase, if the current
      column block does have contributing column
      blocks, an area for simply-linked temporary
      blocks is reserved at least after the sort area,
      leaving enough space to create all of the
      corresponding potential new blocks
      just after all the blocks of the previous
      column block (right picture).


  >     |ccccccccccc| <- bloktab (bloktax)
  >     |ccccccccccc|
  >     |ccccccccccc|                                :ccccccccccc:
  >     |ccccccccccc| >- Computed blocks ----------< |ccccccccccc|
  >     |ccccccccccc|                                |ccccccccccc|
  >     |-----------|                                |:::::::::::|
  >     |hhhhhhhhhhh| <- hashtab = bloknum --------> |bcbcbcbcbcb|
  >     |hhhhhhhhhhh|                 |              |cbcbcbcbcbc|
  >     |hhhhhhhhhhh|                 |              |bcbcbcbcbcb|
  >     |hhhhhhhhhhh|                 |              |cbcbcbcbcbc|
  >     |-----------|                 |              |bcbcbcbcbcb|
  >     |           |                 |              |-----------|
  >     |-----------| <- sorttab...... ------------> |           |
  >     |sssssssssss|                                |           |
  >     |sssssssssss|                                |           |
  >     |-----------| <- ............................|           |
  >     |           |                     tloktab -> |-----------|
  >     |           |                                |ttttttttttt|
  >     |           |                                |ttttttttttt|
  >     :           :                                |-----------|
  >     :___________:                                :___________:
  >                   <- bloktab + blokmax


  Parameters:
    symbptr - Symbolic block matrix [based]
    vertnbr - Number of vertices
    edgenbr - Number of edges
    baseval - Base value
    ngbdptr - Neighbor bookkeeping area
    ngbfrst - First neighbor function
    ngbnext - Next neighbor function
    ngbdegr - Vertex degree function (upper bound)
    ordeptr - Matrix ordering

  Returns:
    0  - on success.
    !0 - on error.

*/

int
pastixSymbolFax( symbol_matrix_t * const symbptr,
                 const pastix_int_t      vertnbr,
                 const pastix_int_t      edgenbr,
                 const pastix_int_t      baseval,
                 void * const            ngbdptr,
                 pastix_int_t            ngbfrst (void * const, const pastix_int_t),
                 pastix_int_t            ngbnext (void * const),
                 pastix_int_t            ngbdegr (void * const, const pastix_int_t),
                 const pastix_order_t * const ordeptr )
#endif /* SYMBOL_FAX_INCLUDED */
{
    pastix_int_t                       vertnum;  /* Vertex number of current column                   */
    pastix_int_t                       vertend;  /* Current end vertex number                         */
    const pastix_int_t * restrict      permtax;  /* Based access to direct permutation array          */
    const pastix_int_t * restrict      peritax;  /* Based access to inverse permutation array         */
    const pastix_int_t * restrict      rangtax;  /* Based access to column block range array          */
    pastix_int_t * restrict            ctrbtax;  /* Based access to array of contribution chains      */
    symbol_cblk_t * restrict           cblktax;  /* Based access to column block array                */
    pastix_int_t                       cblknum;  /* Based number of current column block              */
    pastix_int_t                       cblkctr;  /* Based number of current contributing column block */
    symbol_blok_t * restrict           bloktax;  /* Based access to block array                       */
    pastix_int_t                       bloknum;  /* Based number of current first free block slot     */
    pastix_int_t                       blokmax;  /* Maximum number of blocks in array                 */
    SymbolFaxTlok * restrict           tloktab;  /* Beginning of array of temporary blocks            */
    pastix_int_t                       ctrbsum;  /* Number of contributing blocks for column block    */
    pastix_int_t * restrict            sorttab;  /* Beginning of sort area                            */
    pastix_int_t                       sortnbr;  /* Number of vertices in sort area and hash table    */
    pastix_int_t * restrict            hashtab;  /* Hash vertex table                                 */
    pastix_int_t                       hashmsk;  /* Mask for access to hash table                     */
    pastix_int_t                       colend;   /* Column number of vertex neighbor                  */

    permtax = ordeptr->permtab - baseval;           /* Compute array bases */
    peritax = ordeptr->peritab - baseval;
    rangtax = ordeptr->rangtab - baseval;

    /* Estimate size of initial block array */
    blokmax  = ordeptr->cblknbr * (2 + edgenbr / vertnbr) + 2;

    /* Allocate arrays for factoring   */
    {
        pastix_int_t  *ctrbtab = NULL; /* Array for contribution chaining */
        symbol_cblk_t *cblktab = NULL; /* Column block array              */
        symbol_blok_t *bloktab = NULL; /* Block array                     */

        MALLOC_INTERN(ctrbtab, ordeptr->cblknbr,     pastix_int_t);
        MALLOC_INTERN(cblktab, ordeptr->cblknbr + 1, symbol_cblk_t);
        MALLOC_INTERN(bloktab, blokmax,              symbol_blok_t);

        cblktax = cblktab - baseval;                  /* Set based accesses */
        bloktax = bloktab - baseval;
        ctrbtax = ctrbtab - baseval;

        memset (ctrbtab, ~0, ordeptr->cblknbr * sizeof (pastix_int_t)); /* Initialize column block contributions link array */
    }

    bloknum = baseval;
    for (cblknum = baseval; cblknum < baseval + ordeptr->cblknbr; cblknum ++) { /* For all column blocks */
        pastix_int_t                 colnum;                   /* Number of current column [based]                  */
        pastix_int_t                 colmax;                   /* Maximum column index for current column block     */
        pastix_int_t                 degrsum;
        pastix_int_t                 tlokmax;

        /* Compute offsets and check for array size */
        {
            pastix_int_t                 hashsiz;
            pastix_int_t                 hashmax;
            pastix_int_t                 ctrbtmp;
            pastix_int_t                 sortoft;                /* Offset of sort array                   */
            pastix_int_t                 tlokoft;                /* Offset of temporary block array        */
            pastix_int_t                 tlndoft;                /* Offset of end of temporary block array */

            colnum = rangtax[cblknum];
            colmax = rangtax[cblknum + 1];              /* Get maximum column value */

            cblktax[cblknum].fcolnum = colnum;          /* Set column block data */
            cblktax[cblknum].lcolnum = colmax - 1;
            cblktax[cblknum].bloknum = bloknum;
            cblktax[cblknum].brownum = -1;

            degrsum = 0;
            for ( ; colnum < colmax; colnum ++) { /* For all columns                                  */
                degrsum += SYMBOL_FAX_VERTEX_DEGREE (ngbdptr, peritax[colnum]); /* Add column degrees */
            }

            for (hashmax = 256; hashmax < degrsum; hashmax *= 2); /* Get upper bound on hash table size */
            hashsiz = hashmax << 2;                               /* Fill hash table at 1/4 of capacity */
            hashmsk = hashsiz - 1;

            for (ctrbsum = 0, ctrbtmp = ctrbtax[cblknum]; /* Follow chain of contributing column blocks */
                 ctrbtmp != ~0; ctrbtmp = ctrbtax[ctrbtmp])
            {
                ctrbsum += cblktax[ctrbtmp + 1].bloknum - cblktax[ctrbtmp].bloknum - 2; /* Sum contributing column blocks */
            }

            tlokmax = degrsum + ctrbsum;
            sortoft = tlokmax * sizeof (symbol_blok_t);
            if ((hashsiz * (pastix_int_t)sizeof(pastix_int_t)) > sortoft) {  /* Compute offset of sort area */
                sortoft = (hashsiz * sizeof (pastix_int_t));
            }
            tlokoft = sortoft +  degrsum    * sizeof (pastix_int_t); /* Compute offset of temporary block area */
            tlndoft = tlokoft + (tlokmax+1) * sizeof (SymbolFaxTlok); /* Compute end of area          */

            if (((char *) (bloktax + bloknum) + tlndoft) > /* If not enough room */
                ((char *) (bloktax + blokmax)))
            {
                symbol_blok_t *        bloktmp;              /* Temporary pointer for array resizing */

                do {
                    blokmax = blokmax + (blokmax >> 2) + 4; /* Increase block array size by 25% as long as it does not fit */
                }
                while (((char *) (bloktax + bloknum) + tlndoft) > ((char *) (bloktax + blokmax)));

                if ((bloktmp = (symbol_blok_t *) memRealloc (bloktax + baseval, (blokmax * sizeof (symbol_blok_t)))) == NULL) {
                    errorPrint ("symbolFax: out of memory (2)");
                    memFree    (bloktax + baseval);
                    memFree    (cblktax + baseval);
                    memFree    (ctrbtax + baseval);
                    return     (1);
                }
                bloktax = bloktmp - baseval;
            }

            hashtab = (pastix_int_t *)           (bloktax + bloknum);
            sorttab = (pastix_int_t *)  ((char *) hashtab + sortoft);
            tloktab = (SymbolFaxTlok *) ((char *) hashtab + tlokoft);

            memset (hashtab, ~0, hashsiz * sizeof (pastix_int_t)); /* Initialize hash table */
        }

        sortnbr = 0;                                  /* No vertices yet                 */
        for (colnum = rangtax[cblknum]; colnum < colmax; colnum ++) { /* For all columns */
            pastix_int_t                 hashnum;

            vertnum = peritax[colnum];                  /* Get associated vertex      */
            SYMBOL_FAX_ITERATOR (ngbdptr, vertnum, vertend) /* For all adjacent edges */
            {
                colend = permtax[vertend];                /* Get end column number      */

                if (colend < colmax) {                      /* If end vertex number in left columns */
                    continue;                               /* Skip to next neighbor                */
                }

                for (hashnum = (colend * SYMBOL_FAX_HASHPRIME) & hashmsk; ; /* Search end column in hash table */
                     hashnum = (hashnum + 1) & hashmsk)
                {
                    pastix_int_t *               hashptr;

                    hashptr = hashtab + hashnum;            /* Point to hash slot           */
                    if (*hashptr == colend) {               /* If end column in hash table  */
                        break;                              /* Skip to next end column      */
                    }
                    if (*hashptr == ~0) {                   /* If slot is empty             */
                        *hashptr = colend;                  /* Set column in hash table     */
                        sorttab[sortnbr ++] = colend;       /* Add end column to sort array */
                        break;
                    }
                }
            }
            SYMBOL_FAX_ITERATOR_END;                     /* End of loop on neighbors */
        }                                             /* End of loop on columns   */
        assert( sortnbr <= degrsum );
        intSort1asc1 (sorttab, sortnbr);              /* Sort neighbor array */

        cblkctr = cblknum;
        if (ctrbtax[cblknum] == ~0) {                 /* If column is not to be updated */
            pastix_int_t                 sortnum;

            bloktax[bloknum].frownum = cblktax[cblknum].fcolnum; /* Build diagonal block */
            bloktax[bloknum].lrownum = cblktax[cblknum].lcolnum;
            bloktax[bloknum].lcblknm = cblknum;
            bloktax[bloknum].fcblknm = cblknum;
            bloknum ++;

            for (sortnum = 0; sortnum < sortnbr; ) {    /* For all entries in sorted array */

                colend = sorttab[sortnum];
                if (colend >= rangtax[cblkctr + 1]) {     /* If column block number to be found */
                    pastix_int_t                 cblktmm;            /* Median value                       */
                    pastix_int_t                 cblktmx;            /* Maximum value                      */

                    for (cblkctr ++,                        /* Find new column block by dichotomy */
                             cblktmx = ordeptr->cblknbr + baseval;
                         cblktmx - cblkctr > 1; )
                    {
                        cblktmm = (cblktmx + cblkctr) >> 1;
                        if (rangtax[cblktmm] <= colend) {
                            cblkctr = cblktmm;
                        }
                        else {
                            cblktmx = cblktmm;
                        }
                    }
                }

                bloktax[bloknum].frownum = colend;        /* Set beginning of new block */
                while ((++ sortnum < sortnbr) &&          /* Scan extent of block       */
                       (sorttab[sortnum] - 1 == sorttab[sortnum - 1]) &&
                       (sorttab[sortnum] < rangtax[cblkctr + 1])) ;
                bloktax[bloknum].lrownum = sorttab[sortnum - 1]; /* Set end of block */
                bloktax[bloknum].lcblknm = cblknum;
                bloktax[bloknum].fcblknm = cblkctr;
                bloknum ++;                               /* One more block */
            }
        }
        else {                                   /* Column will be updated           */
            pastix_int_t sortnum;                /* Current index in sort array      */
            pastix_int_t tloknum;                /* Current index on temporary block */
            pastix_int_t tlokfre;                /* Index of first free block        */

            tloktab->frownum = cblktax[cblknum].fcolnum; /* Build diagonal chained block */
            tloktab->lrownum = cblktax[cblknum].lcolnum;
            tloktab->fcblknm = cblknum;
            tloktab->nextnum = 1;
            tloknum = 1;
            assert( tloknum < tlokmax );

            for (sortnum = 0; sortnum < sortnbr; ) {    /* For all entries in sorted array */

                colend = sorttab[sortnum];
                if (colend >= rangtax[cblkctr + 1]) {     /* If column block number to be found */
                    pastix_int_t                 cblktmm;            /* Median value                       */
                    pastix_int_t                 cblktmx;            /* Maximum value                      */

                    for (cblkctr ++,                        /* Find new column block by dichotomy */
                             cblktmx = ordeptr->cblknbr + baseval;
                         cblktmx - cblkctr > 1; )
                    {
                        cblktmm = (cblktmx + cblkctr) >> 1;
                        if (rangtax[cblktmm] <= colend) {
                            cblkctr = cblktmm;
                        }
                        else {
                            cblktmx = cblktmm;
                        }
                    }
                }
                tloktab[tloknum].frownum = colend;        /* Set beginning of new block */
                while ((++ sortnum < sortnbr) &&          /* Scan extent of block       */
                       (sorttab[sortnum] - 1 == sorttab[sortnum - 1]) &&
                       (sorttab[sortnum] < rangtax[cblkctr + 1])) ;
                tloktab[tloknum].lrownum = sorttab[sortnum - 1]; /* Set end of block */
                tloktab[tloknum].fcblknm = cblkctr;
                tloktab[tloknum].nextnum = tloknum + 1;   /* Chain block */
                tloknum = tloknum + 1;
                assert( tloknum < tlokmax );
            }
            tloktab[tloknum].frownum =                  /* Build trailing block */
                tloktab[tloknum].lrownum = vertnbr + baseval;
            tloktab[tloknum].fcblknm = ordeptr->cblknbr + baseval;
            tloktab[tloknum].nextnum = 0;               /* Set end of chain (never chain to diagonal block) */

            tlokfre = ++ tloknum;                       /* Build free chain for possible contributing blocks */
            for ( ; tloknum < tlokmax; tloknum = tloknum + 1)
            {
                tloktab[tloknum].nextnum = tloknum + 1;
            }
            tloktab[tloknum].nextnum = ~0;              /* Set end of free chain */

            for (cblkctr = ctrbtax[cblknum]; cblkctr != ~0; cblkctr = ctrbtax[cblkctr]) { /* Follow chain */
                pastix_int_t                 blokctr;              /* Current index of contributing column block     */
                pastix_int_t                 tloklst;              /* Index of previous temporary block              */

                tloklst = 0;                              /* Previous is diagonal block */
                tloknum = 0;                              /* Current is diagonal block  */

                for (blokctr = cblktax[cblkctr].bloknum + 2; /* For all blocks in contributing column block */
                     blokctr < cblktax[cblkctr + 1].bloknum; blokctr ++) {
                    while ((tloktab[tloknum].fcblknm < bloktax[blokctr].fcblknm) || /* Skip unmatched chained blocks */
                           (tloktab[tloknum].lrownum < bloktax[blokctr].frownum - 1)) {
                        tloklst = tloknum;
                        tloknum = tloktab[tloknum].nextnum;
                    }

                    if ((bloktax[blokctr].fcblknm < tloktab[tloknum].fcblknm) || /* If contributing block has no mate */
                        (bloktax[blokctr].lrownum < tloktab[tloknum].frownum - 1)) {
                        pastix_int_t                 tloktmp;

#ifdef FAX_DEBUG
                        if (tlokfre == ~0) {
                            errorPrint ("symbolFax: internal error (1)");
                            memFree    (bloktax + baseval);
                            memFree    (cblktax + baseval);
                            memFree    (ctrbtax + baseval);
                            return     (1);
                        }
#endif /* FAX_DEBUG */
                        tloktmp                  =
                            tloktab[tloklst].nextnum = tlokfre;   /* Chain new block                */
                        tloktab[tlokfre].frownum = bloktax[blokctr].frownum; /* Copy block data */
                        tloktab[tlokfre].lrownum = bloktax[blokctr].lrownum;
                        tloktab[tlokfre].fcblknm = bloktax[blokctr].fcblknm;
                        tlokfre                  = tloktab[tlokfre].nextnum;
                        tloktab[tloktmp].nextnum = tloknum;   /* Complete chaining                    */
                        tloknum                  = tloktab[tloklst].nextnum; /* Resume from new block */
                        continue;                             /* Process next block                   */
                    }

                    if ((bloktax[blokctr].lrownum >= tloktab[tloknum].frownum - 1) && /* Update chained block lower bound */
                        (bloktax[blokctr].frownum <  tloktab[tloknum].frownum))
                    {
                        tloktab[tloknum].frownum = bloktax[blokctr].frownum;
                    }

                    if ((bloktax[blokctr].frownum <= tloktab[tloknum].lrownum + 1) && /* Update chained block upper bound */
                        (bloktax[blokctr].lrownum >  tloktab[tloknum].lrownum))
                    {
                        pastix_int_t                 tloktmp;

                        tloktab[tloknum].lrownum = bloktax[blokctr].lrownum;

                        for (tloktmp = tloktab[tloknum].nextnum; /* Aggregate following chained blocks */
                             (tloktab[tloktmp].fcblknm == tloktab[tloknum].fcblknm) &&
                                 (tloktab[tloktmp].frownum <= tloktab[tloknum].lrownum + 1);
                             tloktmp = tloktab[tloknum].nextnum )
                        {
                            if (tloktab[tloktmp].lrownum > tloktab[tloknum].lrownum) { /* Merge aggregated block */
                                tloktab[tloknum].lrownum = tloktab[tloktmp].lrownum;
                            }
                            tloktab[tloknum].nextnum = tloktab[tloktmp].nextnum; /* Unlink aggregated block */
                            tloktab[tloktmp].nextnum = tlokfre;
                            tlokfre                  = tloktmp;
                        }
                    }
                }
            }

            for (tloknum = 0;                           /* For all chained blocks                    */
                 tloktab[tloknum].nextnum != 0;         /* Until trailer block is reached            */
                 tloknum = tloktab[tloknum].nextnum, bloknum ++) { /* Copy block data to block array */
                bloktax[bloknum].frownum = tloktab[tloknum].frownum;
                bloktax[bloknum].lrownum = tloktab[tloknum].lrownum;
                bloktax[bloknum].lcblknm = cblknum;
                bloktax[bloknum].fcblknm = tloktab[tloknum].fcblknm;
            }
        }
        if ((bloknum - cblktax[cblknum].bloknum) > 2) { /* If more than one extra-diagonal blocks exist                 */
            ctrbtax[cblknum] = ctrbtax[bloktax[cblktax[cblknum].bloknum + 1].fcblknm]; /* Link contributing column blocks */
            ctrbtax[bloktax[cblktax[cblknum].bloknum + 1].fcblknm] = cblknum;
        }
    }
    cblktax[cblknum].fcolnum =                      /* Set last column block data */
        cblktax[cblknum].lcolnum = vertnbr + baseval;
    cblktax[cblknum].bloknum = bloknum;
    cblktax[cblknum].brownum = -1;

    memFree (ctrbtax + baseval);                    /* Free contribution link array */

    symbptr->baseval = baseval;                     /* Fill in matrix fields */
    symbptr->cblknbr = ordeptr->cblknbr;
    symbptr->bloknbr = bloknum - baseval;
    symbptr->cblktab = cblktax + baseval;
    symbptr->bloktab = (symbol_blok_t *) memRealloc (bloktax + baseval, (bloknum - baseval) * sizeof (symbol_blok_t)); /* Set array to its exact size */
    symbptr->nodenbr = vertnbr;
    symbptr->browtab = NULL;

#ifdef FAX_DEBUG
    if (pastixSymbolCheck (symbptr) != 0) {
        errorPrint ("symbolFax: internal error (2)");
        pastixSymbolExit (symbptr);
        return     (1);
    }
#endif /* FAX_DEBUG */

    return PASTIX_SUCCESS;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
