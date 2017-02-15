/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/
/**
 * \file mmg2d/queue_2d.c
 * \brief Queue data structure
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "mmg2d.h"


pQueue MMG2_kiuini(MMG5_pMesh mesh,int nbel,double declic,int base) {
  pQueue   q;
  MMG5_pTria    pt;
  int      k;

  _MMG5_SAFE_CALLOC(q,1,Queue);
  _MMG5_SAFE_CALLOC(q->stack,(1+nbel),int);

  q->cur = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) || pt->qual < declic )     continue;
    else if ( base > 0 && pt->flag < base-1 )  continue;
    q->stack[q->cur] = k;
    q->cur++;
  }
  return(q);
}


void MMG2_kiufree(pQueue q) {
  _MMG5_SAFE_FREE(q->stack);
  _MMG5_SAFE_FREE(q);
}


int MMG2_kiudel(pQueue q,int iel) {
  int   k;
  if ( !q->stack[0] )
    return(0);

  else if ( q->cur != iel && !q->stack[iel] )
    return(0);

  else {
    /*recherche un peu longue!!!!*/
    for (k=q->cur-1; k>0; k--)
      if ( q->stack[k] == iel )  break;
    if(!k) return(0);
    assert(q->stack[k] == iel);
    if ( k == q->cur-1 ) { /*dernier elt de la pile*/
      q->cur--;
      q->stack[k]   = 0;
      return(1);
    }
    else { /*on met le dernier elt et on le met dans la case*/
      q->stack[k]   = q->stack[q->cur-1];
      q->stack[q->cur-1] = 0;
      q->cur--;
      return(1);
    }
  }

  return(0);
}


int MMG2_kiuput(pQueue q,int iel) {
  int    k;
  puts("on put");
  if ( !q->stack[0] )
    return(0);

  else if ( iel == q->cur || q->stack[iel] )
    return(0);

  else if ( iel > q->cur ) {
    q->stack[q->cur] = iel;
    q->stack[iel]    = 0;
    q->cur = iel;
    return(1);
  }

  else if ( iel < q->stack[0] ) {
    q->stack[iel] = q->stack[0];
    q->stack[0]   = iel;
    return(1);
  }

  else {
    for (k=iel-1; k>=0; k--)
      if ( q->stack[k] )  break;
    q->stack[iel] = q->stack[k];
    q->stack[k]   = iel;
    return(1);
  }

  return(0);
}


int MMG2_kiupop(pQueue q) {
  int  cur;

  /*cur = q->stack[0];
    q->stack[0]   = q->stack[q->cur-1];
    q->stack[q->cur] = 0;
    if(q->cur) q->cur--;*/

  /*on depile par le haut*/
  if(!q->cur) return(0);
  cur = q->stack[q->cur-1];
  q->stack[q->cur-1] = 0;
  q->cur--;

  return(cur);
}


