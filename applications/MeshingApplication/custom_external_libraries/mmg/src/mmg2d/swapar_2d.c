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
 * \file mmg2d/swapar_2d.c
 * \brief Functions for swapping process.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg2d.h"

int MMG2_swapar(MMG5_pMesh mesh,MMG5_pSol sol,int k,int i,double crit,int *list) {
  MMG5_pTria   pt,pt1;
  MMG5_Tria    tmp;
  double  cal1,cal2,air1,air2,airn1,airn2;
  int    *adj,*adj1,k1,k2,k3;
  int    i1,i2,vo1,vo2,vo3;
  int    num,num1,numa1;

  adj = &mesh->adja[3*(k-1)+1];
  pt  = &mesh->tria[k];
  pt1 = &mesh->tria[adj[i] / 3];

  if(pt->ref!=pt1->ref) {
    return(0);
  }

  air1 = MMG2_quickarea(mesh->point[pt->v[0]].c,mesh->point[pt->v[1]].c,mesh->point[pt->v[2]].c);
  air2 = MMG2_quickarea(mesh->point[pt1->v[0]].c,mesh->point[pt1->v[1]].c,mesh->point[pt1->v[2]].c);
  /* simul */
  i1 = MMG2_idir[i+1];
  i2 = MMG2_idir[i+2];

  tmp.v[0] = pt->v[i];
  tmp.v[1] = pt->v[i1];
  tmp.v[2] = pt1->v[adj[i] % 3];
  cal1 = MMG2_caltri_in(mesh,sol,&tmp);
  airn1 = MMG2_quickarea(mesh->point[tmp.v[0]].c,mesh->point[tmp.v[1]].c,mesh->point[tmp.v[2]].c);
  if ( cal1 > crit )  return(0);

  tmp.v[0] = pt->v[i];
  tmp.v[1] = pt1->v[adj[i] % 3];
  tmp.v[2] = pt->v[i2];
  cal2 = MMG2_caltri_in(mesh,sol,&tmp);
  airn2 = MMG2_quickarea(mesh->point[tmp.v[0]].c,mesh->point[tmp.v[1]].c,mesh->point[tmp.v[2]].c);
  if ( cal2 > crit )  return(0);

  if( airn1 < 0 || airn2 < 0 || fabs((air1+air2)-(airn1+airn2)) > EPSD) {
    if(mesh->info.ddebug) printf("  ## Warning: non convex configuration\n");
    return(0);
  }

  /* update structure */
  k1   = adj[i] / 3;
  assert(k1);
  vo1  = adj[i] % 3;
  k2   = adj[i1] / 3;
  vo2  = adj[i1] % 3;
  adj1 = &mesh->adja[3*(k1-1)+1];
  k3   = adj1[MMG2_idir[vo1+1]] / 3;
  vo3  = adj1[MMG2_idir[vo1+1]] % 3;

  pt->v[i2] = pt1->v[vo1];
  pt->qual  = cal1;
  list[1] = k;
  pt1->v[MMG2_idir[vo1+2]] = pt->v[i];
  pt1->qual = cal2;
  list[2] = k1;

  /*update edge*/
  num = pt->edg[i];
  assert(!num);
  num1 = pt->edg[i1];
  numa1 = pt1->edg[MMG2_idir[vo1+1]];

  if(k2) (&mesh->adja[3*(k2-1)+1])[vo2] = 3*k1+vo1;
  pt1->edg[vo1] = num1;
  if(k1) mesh->adja[3*(k1-1)+1+vo1] = 3*k2+vo2;
  vo1  = MMG2_idir[vo1+1];
  mesh->adja[3*(k-1)+1+i]    = 3*k3+vo3;
  pt->edg[i] = numa1;
  if(k3) mesh->adja[3*(k3-1)+1+vo3] = 3*k+i;
  mesh->adja[3*(k-1)+1+i1]   = 3*k1+vo1;
  pt->edg[i1] = 0;
  if(k1) mesh->adja[3*(k1-1)+1+vo1] = 3*k+i1;
  pt1->edg[vo1] = 0;


  return(1);
}
