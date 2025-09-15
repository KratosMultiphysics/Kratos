/*
--------------------------------------------------------------------
lookupa.h, by Bob Jenkins, December 1996.  Same as lookup2.c
Use this code however you wish.  Public Domain.  No warranty.
Source is http://burtleburtle.net/bob/c/lookupa.c
--------------------------------------------------------------------
*/

#ifndef _LOOKUPA_H_
#define _LOOKUPA_H_

/*
--------------------------------------------------------------------
lookup() -- hash a variable-length key into a 32-bit value
  k     : the key (the unaligned variable-length array of bytes)
  len   : the length of the key, counting by bytes
  level : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = lookup( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

ub4 lookup( /* register ub1 *k,
               register ub4  length,
               register ub4  level */);


/*
--------------------------------------------------------------------
checksum() -- hash a variable-length key into a 256-bit value
  k     : the key (the unaligned variable-length array of bytes)
  len   : the length of the key, counting by bytes
  state : an array of CHECKSTATE 4-byte values (256 bits)
The state is the checksum.  Every bit of the key affects every bit of
the state.  There are no funnels.  About 112+6.875len instructions.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0; i<8; ++i) state[i] = 0x9e3779b9;
  for (i=0, h=0; i<n; ++i) checksum( k[i], len[i], state);

(c) Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial, as long
as this whole comment accompanies it.

See http://burtleburtle.net/bob/hash/evahash.html
Use to detect changes between revisions of documents, assuming nobody
is trying to cause collisions.  Do NOT use for cryptography.
--------------------------------------------------------------------
*/
void  checksum( /* register ub1 *k,
                   register ub4  len,
                   register ub4 *state */ );

#endif
