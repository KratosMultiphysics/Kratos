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
 * \file mmgs/gentools_s.c
 * \brief Generic algebraic and algorithmic tools.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"


/* Delete all triangle references in mesh */
int delref(MMG5_pMesh mesh) {
  MMG5_pTria    pt;
  int      k;

  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    pt->ref = 0;
  }

  return(1);
}

/* Start from triangle start, and pile up triangles by adjacency, till a GEO or REF curve is met ;
   pass all references of travelled faces to ref ; putreq = 1 if boundary edges met must
   be set to MG_REQ, 0 otherwise. */
int setref(MMG5_pMesh mesh,int start,int ref,int putreq) {
  MMG5_pTria      pt,pt1;
  int        *list,*adja,cur,base,k,iel,jel,ilist;
  char       j,voy;

  ilist = cur = 0;
  _MMG5_SAFE_CALLOC(list,mesh->nt+1,int);
  base = ++mesh->base;

  /* Pile up triangles from start, till a GEO boundary is met */
  pt = &mesh->tria[start];
  list[ilist] = start;
  ilist++;
  assert( ilist <= mesh->nt );
  pt->flag = base;

  do {
    iel = list[cur];
    pt = &mesh->tria[iel];
    adja = &mesh->adja[3*(iel-1)+1];

    for(j=0; j<3; j++) {
      if( MG_EDG(pt->tag[j]) ) {
        if( putreq ) {
          pt->tag[j] |= MG_REQ;
          jel = adja[j] / 3;
          voy = adja[j] % 3;
          if( !jel ) continue;
          pt1 = &mesh->tria[jel];
          pt1->tag[voy] |= MG_REQ;
        }
        continue;
      }
      jel = adja[j] / 3;
      assert(jel);
      pt1 = &mesh->tria[jel];
      if ( pt1->flag == base )  continue;

      list[ilist] = jel;
      ilist++;
      assert( ilist <= mesh->nt );
      pt1->flag = base;
    }
    cur++;
  }
  while( cur < ilist );

  /* Set all references of triangles of list to ref */
  for (k=0; k<ilist; k++) {
    iel = list[k];
    pt  = &mesh->tria[iel];
    pt->ref = ref;
  }

  return(1);
}

/** find the element number in packed numerotation */
int _MMGS_indElt(MMG5_pMesh mesh, int kel) {
  MMG5_pTria pt;
  int    ne, k;

  ne = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( MG_EOK(pt) ) {
      ne++;
      if ( k == kel )  return(ne);
    }
  }
  return(0);
}

/** find the point number in packed numerotation */
int _MMGS_indPt(MMG5_pMesh mesh, int kp) {
  MMG5_pPoint ppt;
  int         np, k;

  np = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      np++;
      if ( k == kp )  return(np);
    }
  }
  return(0);
}
