#include <stdio.h>
#include <stdlib.h>
#include "globheads.h"
#include "protos.h"

typedef struct __KeyType
{
  int var;   /* row number */
  int key;   /* hash value */
} KeyType;

int KeyComp( const void *vfst, const void *vsnd )
{
  KeyType *fst = (KeyType *)vfst, *snd = (KeyType *)vsnd;
  if( fst->key == snd->key ) {
    if( fst->var < snd->var )
      return -1;
    return 1;
  }
  if( fst->key < snd->key )
    return -1;
  return 1;
}

int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm,
                 double eps, double *t_hash, double *t_angle )
{
/*----------------------------------------------------------------------------
 * Setup Blocks ( rows and columns might be permuted to get better results )
 *----------------------------------------------------------------------------
 * Na Li, Aug 2001
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * csmat   = a matrix stored in SpaFmt format
 * eps     = parameter for deciding when to do a union of two rows
 *           into the same group.  Two rows u and v are merged into a 
 *           block  when cos(<u,v>) == (u,v)/(|u|*|v|), is > eps. 
 *           eps should be <= 1. 
 *----------------------------------------------------------------------------
 * on return:
 * ==========
 * csmat   = matrix stored in SpaFmt format after permutation
 * pnBlock = dimension of the block matrix
 * pnB     = dimension of each block
 *
 *----------------------------------------------------------------------------
 * Combination of hash method and angle method:
 *----------------------------------------------------------------------------
 * Designed for the matrices with symmetric patterns
 * (1) Hash method
 *     a. Calculate hash values
 *     b. qsort rows according to their hash values
 *     c. Get compressed graph as the following format:
 * (2) Angle method
 *     a. Calculate A^T
 *     b. for i-th row, calculate dot product (row_i, row_j) using A*A^T
 *        algorithm where j = i+1, ..., n-1 and group[j] == -1
 *        if cos( <row_i, row_j> ) = (row_i,row_j)/|row_i||row_j| is > eps,
 *        we merge row_i and row_j by resetting
 *        group[j] = i and size[i] = size[i]+size[j]
 *--------------------------------------------------------------------------*/
  int n = csmat->n, nBlock = 0, i, j, k;
  csptr at = NULL;
  KeyType *group = NULL;
  CompressType *compress = NULL;
  int *perm = NULL, *nB = NULL;
  int nzcount0, nzcount, key0, key, *ja0, *ja, row0, row, newblock;
  int *iw = NULL, *jbuf = NULL;
  int cnt, pos, nnz_i, row_j, col, bkcnt;
  int nextBlockID, nextBlockPos, belongTo, grp;
  double eps_2 = eps * eps, t1, t2;

  t1 = sys_timer(); /* begin Hash method timer */
  group = (KeyType *)Malloc( n*sizeof(KeyType), "init_blocks" );
  compress = (CompressType *)Malloc( n*sizeof(CompressType), "init_blocks" );
  perm = (int *)Malloc( n * sizeof(int), "init_blocks" );
  iw = perm; /* iw and perm array can share memory here because they will
	      * never be used at the same time */
  for( i = 0; i < n; i++ ) {
    iw[i] = 0;
    compress[i].grp = -1;
  }
/*-------------------- compress matrix based on hash algorithm */
/*-------------------- get hash value of each row */
  for( i = 0; i < n; i++ ) {
    nzcount = csmat->nzcount[i];
    key = 0;
    ja = csmat->ja[i];
    for( j = 0; j < nzcount; j++ )
      key += ja[j]+1;
    group[i].key = key;
    group[i].var = i;
  }
/*-------------------- sort rows -- uses function KeyComp */
  qsort( group, n, sizeof(KeyType), KeyComp );

/*-------------------- compress matrix */
  for( i = 0; i < n; i++ ) {
    row0 = group[i].var;
    if( compress[row0].grp != -1 ) continue; /* already assigned */
    key0 = group[i].key;
    nzcount0 = csmat->nzcount[row0];
    ja0 = csmat->ja[row0];
/*-------------------- beginning of new block. set .grp and .count */
    compress[row0].grp = -1;
    compress[row0].count = 1;
/*-------------------- loop over all rows having same check-sum keys */
    for( j = i + 1; j < n; j++ ) {
      key = group[j].key;
      if( key != key0 ) break;
      row = group[j].var;
      if( compress[row].grp != -1 ) continue; /* already assigned */
      nzcount = csmat->nzcount[row];
      if( nzcount != nzcount0 ) continue;
      ja = csmat->ja[row];
      newblock = 0;
/*-------------------- compare patterns of the rows             */
      for( k = 0; k < nzcount; k++ ) iw[ ja0[k] ] = 1;
      for( k = 0; k < nzcount; k++ ) {
        if( iw[ ja[k] ] == 0 ) {
          newblock = 1;
          break;
        }
      }
      for( k = 0; k < nzcount; k++ ) iw[ ja0[k] ] = 0; /* reset iw */
/*-------------------- row belongs to group row0                    */ 
      if( !newblock ) {
        compress[row].grp = row0;
        compress[row0].count++;
      }
    }
  }
  t2 = sys_timer(); /* end Hash method timer */
  *t_hash = t2 - t1;

  t1 = sys_timer(); /* begin angle method timer */
  nB = (int *)Malloc( n * sizeof(int), "init_blocks" );
  jbuf = (int *)Malloc( n * sizeof(int), "init_blocks" );

/*-------------------- compress matrix based on angle algorithm */
/*-------------------- calculate compressed A^T                 */
  at = (csptr)Malloc( sizeof(SparMat), "init_blocks" );
  setupCS( at, n, 0 );
  if( CSparTran( csmat, at, compress ) != 0 ) 
    return -1;

/*----------------------------------------------------------------------------
 * only the row representing beginning of block satisfies:
 *    compress[row].grp = -1, so far.
 * how many such rows is up to the compression rate of hash compression
 * algorithm we did above
 *--------------------------------------------------------------------------*/

/*---------------------------------------------------------------
 * use group to backup original compressed matrix by Hash method.
 * It is very important because the array 'compress' will be changed
 * during Angle method. Or, we'll get incorrect inner product.
 *--------------------------------------------------------------*/
  for( i = 0; i < n; i++ ) {
    group[i].var = compress[i].grp;
    group[i].key = compress[i].count;
  }

  for( i = 0; i < n; i++ ) {
    if( compress[i].grp != -1 ) continue;
    nB[nBlock] = compress[i].count; /* !!! not 1 here */
    cnt = 0;
/*-------------------- calculate (u,v_j ), j = i+1,...,n, using product
 *-------------------- algorithm of A * A_T */
    nnz_i = csmat->nzcount[i];
    for( j = 0; j < nnz_i; j++ ) {
      row_j = csmat->ja[i][j];
      if( group[row_j].var != -1 ) /* i.e. original compress[row_j].grp */
        continue;
      bkcnt = group[row_j].key;    /* i.e. original compress[row_j].count */
      for( k = at->nzcount[row_j] - 1; k >= 0; k-- ) {
	col = at->ja[row_j][k];
	if( col <= i ) break;
        if( compress[col].grp != -1 ) continue; /* needed because compress
                                           array is dynamically updated */
	if( iw[col] == 0 ) { /* new nonzero of (u,v_j) */
	  jbuf[cnt] = col;
	  cnt++;
	}
        iw[col] += bkcnt; /* correct for matrix with symmetric pattern */
      }
    }
/*-------------------- set group for row i and reset iw */
    for( j = 0; j < cnt; j++ ) {
      pos = jbuf[j];
      if( iw[pos] * iw[pos] >= eps_2 * nnz_i * csmat->nzcount[pos] ) {
	compress[pos].grp = i;
	nB[nBlock] += compress[pos].count; /* !!! not 1 here */
      }
      iw[pos] = 0; /* reset iw */
    }
    nBlock++; /* begin new block, add block count by 1 */
  } /* end loop i */

/*-------------------- free group                                   */
  if( group ) {
    /* no need group array any more */
    free( group );
    group = NULL;
  }

  *pnBlock = nBlock;
  *pnB = (int *)Malloc( nBlock * sizeof(int), "init_blocks" );
  for( i = 0; i < nBlock; i++ ) {
    if( nB[i] > MAX_BLOCK_SIZE ) {
      fprintf( stderr, "Block of size = %d exceeds MAX_BLOCK_SIZE\n", nB[i] );
      return -1;
    }
    (*pnB)[i] = nB[i];
  }

/*-------------------- calculate permutation array -  Array nB will
 * be used to store next  available position in each  block */
  nextBlockID = 0;
  nextBlockPos = 0;
  for( i = 0; i < n; i++ ) {
    if( compress[i].grp == -1 ) {
      perm[i] = nextBlockPos;
      nextBlockPos += (*pnB)[nextBlockID++];
      nB[i] = 1;
    } else {
      belongTo = compress[i].grp;
      grp = compress[belongTo].grp;
      if( grp != -1 ) /* should find the final beginning of block */
	belongTo = grp;
      perm[i] = perm[belongTo] + nB[belongTo];
      nB[belongTo]++;
    }
  }
  t2 = sys_timer(); /* end angle method timer */
  *t_angle = t2 - t1;
  *pperm = perm;

  cleanCS( at );
  free( nB );
  free( jbuf );
  free( compress );

  return 0;
}
