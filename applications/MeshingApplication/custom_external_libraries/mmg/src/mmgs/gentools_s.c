/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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

#include "libmmgs_private.h"


/* Delete all triangle references in mesh */
int delref(MMG5_pMesh mesh) {
  MMG5_pTria    pt;
  MMG5_int      k;

  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    pt->ref = 0;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param start index of the tetra from which we start
 * \param ref reference to set
 * \param putreq 1 if boundary edges must be set to required
 *
 * \return 1 if success, 0 if fail
 *
 * Start from triangle start, and pile up triangles by adjacency, till a GEO or
 * REF curve is met ; pass all references of travelled faces to ref ; putreq = 1
 * if boundary edges met must be set to MG_REQ, 0 otherwise.
 *
 */
int setref(MMG5_pMesh mesh,MMG5_int start,MMG5_int ref,int putreq) {
  MMG5_pTria pt,pt1;
  MMG5_int   base,*list,*adja,cur,k,iel,jel;
  int        ilist;
  int8_t     j,voy;

  ilist = cur = 0;
  MMG5_SAFE_CALLOC(list,mesh->nt+1,MMG5_int,return 0);
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
  MMG5_SAFE_FREE(list);
  return 1;
}

/** find the element number in packed numerotation */
MMG5_int MMGS_indElt(MMG5_pMesh mesh, MMG5_int kel) {
  MMG5_pTria pt;
  MMG5_int   ne, k;

  ne = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( MG_EOK(pt) ) {
      ne++;
      if ( k == kel )  return ne;
    }
  }
  return 0;
}

/** find the point number in packed numerotation */
MMG5_int MMGS_indPt(MMG5_pMesh mesh, MMG5_int kp) {
  MMG5_pPoint ppt;
  MMG5_int    np, k;

  np = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      np++;
      if ( k == kp )  return np;
    }
  }
  return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param nsd index of subdomain to keep.
 *
 * Keep only subdomain of index \a nsd and remove other subdomains.
 *
 */
void MMGS_keep_only1Subdomain ( MMG5_pMesh mesh,MMG5_int nsd ) {

  if ( !nsd ) {
    return;
  }

  if ( mesh->info.imprim > 4 || mesh->info.ddebug ) {
    fprintf(stdout,"\n  -- ONLY KEEP DOMAIN OF REF %"MMG5_PRId"\n",nsd );
  }

  MMG5_mark_verticesAsUnused ( mesh );

  MMG5_keep_subdomainElts ( mesh, nsd, MMGS_delElt );

  MMG5_mark_usedVertices ( mesh,MMGS_delPt );

  return;
}
