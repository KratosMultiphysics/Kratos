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
#include "mmg2d.h"


static unsigned char inxt[3]  = {1,2,0};
static unsigned char iprev[3] = {2,0,1};


/* find all triangles sharing P
   in:  ifirst    : triangle containing p
   iploc     : index of p in start
   out: list  : list of triangles */
int MMG2_boulep(MMG5_pMesh mesh, int ifirst, int iploc, int * list) {
  MMG5_pTria  pt;
  MMG5_pPoint ppt;
  int    ip,voy,ilist,iel,*adja,i,iadr;

  if ( ifirst < 1 ) return(0);
  pt = &mesh->tria[ifirst];
  if ( !M_EOK(pt) ) return(0);
  ip = pt->v[iploc];
  ppt = &mesh->point[ip];
  if ( !M_VOK(ppt) ) return(0);

  /* init list */
  ilist       = 1;
  list[ilist] = 3*ifirst + iploc;

  iadr = 3*(ifirst-1) + 1;
  adja = &mesh->adja[iadr];
  iel  = adja[inxt[iploc]]/3;
  voy  = adja[inxt[iploc]]%3;
  i    = inxt[voy];

  while ( iel && (iel != ifirst) && mesh->tria[iel].v[0]){
    if(ilist==MMG2D_LMAX-1) return(0);
    list[++ilist] = 3*iel + i;
    assert( ip==(&mesh->tria[iel])->v[i] );
    iadr = 3*(iel-1) + 1;
    adja = &mesh->adja[iadr];
    iel  = adja[inxt[i]]/3;
    voy = adja[inxt[i]]%3;
    i   = inxt[voy];
  }

  if ( iel!=ifirst ) {
    iadr = 3*(ifirst-1) + 1;
    adja = &mesh->adja[iadr];
    iel  = adja[iprev[iploc]]/3;
    voy  = adja[iprev[iploc]]%3;
    i    = iprev[voy];

    while ( iel && (iel != ifirst) && mesh->tria[iel].v[0]) {
      if(ilist==MMG2D_LMAX-1) return(0);
      list[++ilist] = 3*iel + i;
      assert( ip==(&mesh->tria[iel])->v[i] );
      iadr = 3*(iel-1) + 1;
      adja = &mesh->adja[iadr];
      iel  = adja[iprev[i]]/3;
      if (!iel) break;
      voy = adja[iprev[i]]%3;
      i   = iprev[voy];

    }
  }

  return(ilist);
}
