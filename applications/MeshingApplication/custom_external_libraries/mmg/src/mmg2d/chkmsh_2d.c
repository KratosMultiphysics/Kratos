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

/**
 * \param mesh pointer toward the mesh structure.
 * \param severe level of performed check
 * \param base unused argument.
 * \return 0 if fail, 1 if success.
 *
 * Check the mesh validity
 *
 */
int _MMG5_mmg2dChkmsh(MMG5_pMesh mesh, int severe,int base) {
//  MMG5_pPoint ppt;
  MMG5_pTria  pt1,pt2;
  MMG5_pEdge         ped;
  int   *adja,*adja1,adj,adj1,k,i,iadr;
//  int   kk,l,nk,j,ip,lon,len;
  int     *list;
  unsigned char voy,voy1;

  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if ( !M_EOK(pt1) )  continue;
    iadr = (k-1)*3 + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<3; i++) {
      adj = adja[i] / 3;
      voy = adja[i] % 3;
      if ( !adj )  continue;

      if ( adj == k ) {
        fprintf(stdout,"  1. Wrong adjacency %d %d\n",k,adj);
        printf("vertices of %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
        printf("adj of %d: %d %d %d \n",
               k,adja[0]/3,adja[1]/3,adja[2]/3);
        return(0);
      }
      pt2 = &mesh->tria[adj];
      if ( !M_EOK(pt2) ) {
        fprintf(stdout,"  4. Invalid adjacent %d %d\n",adj,k);
        printf("vertices of %d: %d %d %d\n",
               k,pt1->v[0],pt1->v[1],pt1->v[2]);
        printf("vertices adj %d: %d %d %d \n",
               adj,pt2->v[0],pt2->v[1],pt2->v[2]);
        printf("adj of %d: %d %d %d\n",k,adja[0]/3,adja[1]/3,adja[2]/3);
        return(0);
      }
      iadr  = (adj-1)*3 + 1;
      adja1 = &mesh->adja[iadr];
      adj1  = adja1[voy] / 3;
      voy1  = adja1[voy] % 3;
      if ( adj1 != k || voy1 != i ) {
        fprintf(stdout,"  2. Wrong adjacency %d %d\n",k,adj1);
        printf("vertices of %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
        printf("adj(k) %d: %d %d %d \n",adj,pt2->v[0],pt2->v[1],pt2->v[2]);
        printf("adj(%d): %d %d %d\n",
               k,adja[0]/3,adja[1]/3,adja[2]/3);
        printf("adj(%d): %d %d %d %d\n",
               adj,adja1[0]/3,adja1[1]/3,adja1[2]/3,adja1[3]/3);
        return(0);
      }

      /*chk edge*/
      if(pt1->edg[i]) {
        ped = &mesh->edge[pt1->edg[i]];
        if(!(((ped->a==pt1->v[MMG2_iare[i][0]]) || (ped->a==pt1->v[MMG2_iare[i][1]]))
             || ((ped->b==pt1->v[MMG2_iare[i][0]]) || (ped->b==pt1->v[MMG2_iare[i][1]])))) {
          printf("  3. Wrong edge in triangle %d\n",k);
          printf("vertices of %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
          printf("edge %d : %d %d\n",i,ped->a,ped->b);
          return(0);
        }
      }

    }
  }

  if ( !severe )  return(1);

  _MMG5_SAFE_CALLOC(list,MMG2D_LMAX,int);

  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if ( !M_EOK(pt1) )  continue;
    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[iadr];

    /*for (i=0; i<3; i++) {
      adj = (adja[i]-1) / 3 + 1;
      voy = (adja[i]-1) % 3;
      if ( !adj )  continue;

      ip  = pt1->v[i];
      ppt = &mesh->point[ip];
      if ( !M_VOK(ppt) ) {
      fprintf(stdout,"  6. Unused vertex %d  %d\n",k,ip);
      printf("%d %d %d\n",pt1->v[0],pt1->v[1],pt1->v[2]);
      return(0);
      }
      lon = boulep(mesh,k,i,list);
      for (l=1; l<=lon; l++) {
      kk  = list[l] / 3;
      nk  = list[l] % 3;
      pt2 = &mesh->tria[kk];
      if ( pt2->v[nk] != ip ) {
      fprintf(stdout,"  5. Wrong ball %d, %d\n",ip,pt2->v[nk]);
      return(0);
      }
      }
      if ( lon < 1 )  continue;
      len = 0;
      for (kk=1; kk<=mesh->nt; kk++) {
      pt2 = &mesh->tria[kk];
      if ( !pt2->v[0] )  continue;
      for (j=0; j<3; j++)
      if ( pt2->v[j] == ip ) {
      len++;
      break;
      }
      }
      if ( len != lon ) {
      fprintf(stdout,"  7. Incorrect ball %d: %d %d\n",pt1->v[i],lon,len);
      return(0);
      }
      } */
  }
  _MMG5_SAFE_FREE(list);
  return(1);
}
