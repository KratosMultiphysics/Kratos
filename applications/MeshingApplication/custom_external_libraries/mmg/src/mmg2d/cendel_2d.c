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



int MMG2_cendel(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base) {
  MMG5_pTria      pt,pt1;
  pQueue     queue;
  double      crit;
  int       *adja,*list,adj,iadr,i,k,ns,np;

  /* queue on quality */
  queue = MMG2_kiuini(mesh,mesh->nt,declic,-1);
  assert(queue);
  list  = (int*)malloc(MMG2D_LMAX*sizeof(int));
  assert(list);
  ns = 0;
  np = 0;
  do {
    k = MMG2_kiupop(queue);
    if ( !k )  break;
    np++;
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;

    /* base internal edges */
    iadr  = 3*(k-1) + 1;
    adja  = &mesh->adja[iadr];
    for (i=0; i<3; i++) {
      adj = adja[i] / 3;
      if ( !adj || pt->ref != mesh->tria[adj].ref )  continue;
      //check required
      if((mesh->point[pt->v[MMG2_iare[i][0]]].tag & M_REQUIRED)
         && (mesh->point[pt->v[MMG2_iare[i][1]]].tag & M_REQUIRED)) {
        continue;
      }

      pt1  = &mesh->tria[adj];
      crit = 0.99 * M_MAX(pt->qual,pt1->qual);
      if ( MMG2_swapar(mesh,sol,k,i,crit,list) ) {
        ns++;
        break;
      }

    }
  }
  while ( k );

  if ( mesh->info.imprim < - 4 )
    fprintf(stdout,"     %7d PROPOSED  %7d SWAPPED\n",np,ns);

  MMG2_kiufree(queue);
  free(list);
  return(ns);
}
