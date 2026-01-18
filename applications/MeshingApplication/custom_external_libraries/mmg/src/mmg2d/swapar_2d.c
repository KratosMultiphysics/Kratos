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
 * \file mmg2d/swapar_2d.c
 * \brief Functions for swapping process.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg2d_private.h"

/* Version of edge swapping specific to the boundary enforcement stage in Delaunay meshing;
   the quality of the resulting criterion should be > crit.
   list returns both modified triangles */
int MMG2D_swapdelone(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int k,int8_t i,double crit,MMG5_int *list) {
  MMG5_pTria         pt,pt1,pt0;
  double             cal1,cal2,area1,area2,arean1,arean2;
  MMG5_int           *adja,*adja1,k1,k2,k3,vo2,vo3,num1,numa1;
  int8_t             i1,i2,j,j1,j2;

  adja = &mesh->adja[3*(k-1)+1];
  k1  = adja[i] / 3;
  if ( !k1 ) return 0;

  j    = adja[i] % 3;
  j1    = MMG5_inxt2[j];
  j2    = MMG5_iprv2[j];
  pt0  = &mesh->tria[0];
  pt   = &mesh->tria[k];
  pt1  = &mesh->tria[k1];

  if(pt->ref!=pt1->ref) {
    return 0;
  }

  area1 = MMG2D_quickarea(mesh->point[pt->v[0]].c,mesh->point[pt->v[1]].c,mesh->point[pt->v[2]].c);
  area2 = MMG2D_quickarea(mesh->point[pt1->v[0]].c,mesh->point[pt1->v[1]].c,mesh->point[pt1->v[2]].c);

  /* Simulate the alternate configuration */
  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];

  pt0->v[0] = pt->v[i];
  pt0->v[1] = pt->v[i1];
  pt0->v[2] = pt1->v[j];
  cal1 = MMG2D_caltri_iso(mesh,sol,pt0);
  arean1 = MMG2D_quickarea(mesh->point[pt0->v[0]].c,mesh->point[pt0->v[1]].c,mesh->point[pt0->v[2]].c);
  if ( cal1 > crit )  return 0;

  pt0->v[0] = pt->v[i];
  pt0->v[1] = pt1->v[j];
  pt0->v[2] = pt->v[i2];
  cal2 = MMG2D_caltri_iso(mesh,sol,pt0);
  arean2 = MMG2D_quickarea(mesh->point[pt0->v[0]].c,mesh->point[pt0->v[1]].c,mesh->point[pt0->v[2]].c);
  if ( cal2 > crit )  return 0;

  if ( arean1 < 0.0 || arean2 < 0.0 || fabs((area1+area2)-(arean1+arean2)) > MMG2D_EPSD ) {
    if(mesh->info.ddebug) printf("  ## Warning: non convex configuration\n");
    return 0;
  }

  /* Update vertices of both triangles */
  k2    = adja[i1] / 3;
  vo2   = adja[i1] % 3;
  adja1 = &mesh->adja[3*(k1-1)+1];
  k3    = adja1[j1] / 3;
  vo3   = adja1[j1] % 3;

  pt->v[i2] = pt1->v[j];
  pt->qual  = cal1;
  list[1] = k;
  pt1->v[j2] = pt->v[i];
  pt1->qual = cal2;
  list[2] = k1;

  /* Update edge references */
#ifndef NDEBUG
  MMG5_int num = pt->edg[i];
  assert ( !num );
#endif
  num1 = pt->edg[i1];
  numa1 = pt1->edg[j1];

  /* Update adjacencies */
  mesh->adja[3*(k1-1)+1+j] = 3*k2+vo2;
  pt1->edg[j] = num1;
  if ( k2 )
    mesh->adja[3*(k2-1)+1+vo2] = 3*k1+j;

  mesh->adja[3*(k-1)+1+i]    = 3*k3+vo3;
  pt->edg[i] = numa1;
  if ( k3 )
    mesh->adja[3*(k3-1)+1+vo3] = 3*k+i;

  mesh->adja[3*(k-1)+1+i1]   = 3*k1+j1;
  pt->edg[i1] = 0;
  mesh->adja[3*(k1-1)+1+j1] = 3*k+i1;
  pt1->edg[j1] = 0;

  return 1;
}

/* Check whether swap of edge i in triangle k is valid, and suitable for the mesh */
int MMG2D_chkswp(MMG5_pMesh mesh, MMG5_pSol met,MMG5_int k,int8_t i,int8_t typchk) {
  MMG5_pTria          pt,pt0,pt1;
  double              /*loni,lona,*/cal1,cal2,calnat,calchg;
  MMG5_int            *adja,ip,ip1,ip2,iq,kk;
  uint8_t             i1,i2,ii,ii1,ii2;

  pt0 = &mesh->tria[0];
  pt  = &mesh->tria[k];
  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];
  if ( MG_EDG(pt->tag[i]) || MG_SIN(pt->tag[i]) )  return 0;

  ip = pt->v[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  adja = &mesh->adja[3*(k-1)+1];
  if ( !adja[i] )  return 0;

  kk = adja[i] / 3;
  ii = adja[i] % 3;
  ii1 = MMG5_inxt2[ii];
  ii2 = MMG5_iprv2[ii];

  pt1 = &mesh->tria[kk];
  iq = pt1->v[ii];

  /* If mesh->info.fem : avoid creating a non BDY edge with BDY endpoints */
  if ( mesh->info.fem ) {
    if ( (mesh->point[ip].tag & MG_BDY) && (mesh->point[iq].tag & MG_BDY) ) return 0;
  }

  /* Check length in typchk = 2 mode ; prevent swap if the created edge is
   longer than the swapped one, and than the maximum authorized edge length */
  /* I believe this test is hindering the reach of good quality */
  /* if ( typchk == 2 && met->m ) { */
  /*   loni = MMG2D_lencurv(mesh,met,ip1,ip2); */
  /*   lona = MMG2D_lencurv(mesh,met,ip,iq); */
  /*   if ( loni > 1.0 )  loni = MG_MIN(1.0 / loni,MMG2D_LSHRT); */
  /*   if ( lona > 1.0 )  lona = 1.0 / lona; */
  /*   if ( lona < loni )  return 0; */
  /* } */

  /* Check qualities: see possible bug in mmgs + correct for metric (use in anisotropic?) */
  if ( typchk == 2 && met->m && met->size == 3 ) {
    /* initial quality */
    pt0->v[0]= ip;  pt0->v[1]= ip1;  pt0->v[2]= ip2;
    pt0->tag[0] = pt->tag[i];
    pt0->tag[1] = pt->tag[i1];
    pt0->tag[2] = pt->tag[i2];
    cal1 = MMG2D_caltri_ani(mesh,met,pt0);

    pt0->v[0]= ip1;  pt0->v[1]= iq;   pt0->v[2]= ip2;
    pt0->tag[0] = pt1->tag[ii2];
    pt0->tag[1] = pt1->tag[ii];
    pt0->tag[2] = pt1->tag[ii1];
    cal2 = MMG2D_caltri_ani(mesh,met,pt0);

    calnat = MG_MIN(cal1,cal2);
    assert(calnat > 0.);

    /* quality after swap */
    pt0->v[0]= ip;  pt0->v[1]= ip1;  pt0->v[2]= iq;
    pt0->tag[0] = pt1->tag[ii1];
    pt0->tag[1] = MG_NUL;
    pt0->tag[2] = pt->tag[i2];
    cal1 = MMG2D_caltri_ani(mesh,met,pt0);

    pt0->v[0]= ip;  pt0->v[1]= iq;   pt0->v[2]= ip2;
    pt0->tag[0] = pt1->tag[ii2];
    pt0->tag[1] = pt->tag[i1];
    pt0->tag[2] = MG_NUL;
    cal2 = MMG2D_caltri_ani(mesh,met,pt0);

    calchg = MG_MIN(cal1,cal2);
  }
  else {
    pt0->v[0]= ip;  pt0->v[1]= ip1;  pt0->v[2]= ip2;
    cal1 = MMG2D_caltri_iso(mesh,NULL,pt0);
    pt0->v[0]= ip1;  pt0->v[1]= iq;   pt0->v[2]= ip2;
    cal2 = MMG2D_caltri_iso(mesh,NULL,pt0);
    calnat = MG_MIN(cal1,cal2);
    pt0->v[0]= ip;  pt0->v[1]= ip1;  pt0->v[2]= iq;
    cal1 = MMG2D_caltri_iso(mesh,NULL,pt0);
    pt0->v[0]= ip;  pt0->v[1]= iq;   pt0->v[2]= ip2;
    cal2 = MMG2D_caltri_iso(mesh,NULL,pt0);
    calchg = MG_MIN(cal1,cal2);
  }

  return calchg > 1.01 * calnat;
}

/* Effective swap of edge i in triangle k */
int MMG2D_swapar(MMG5_pMesh mesh,MMG5_int k,int8_t i) {
  MMG5_pTria    pt,pt1;
  MMG5_int      *adja,adj,k11,k21;
  int8_t        i1,i2,j,jj,j2,v11,v21;

  pt   = &mesh->tria[k];
  if ( MG_EDG(pt->tag[i]) || MG_SIN(pt->tag[i]) )  return 0;

  adja = &mesh->adja[3*(k-1)+1];
  assert(adja[i]);

  adj = adja[i] / 3;
  j   = adja[i] % 3;
  pt1 = &mesh->tria[adj];

  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];

  /* update structure */
  k11 = adja[i1] / 3;
  v11 = adja[i1] % 3;
  adja = &mesh->adja[3*(adj-1)+1];
  jj  = MMG5_inxt2[j];
  j2  = MMG5_iprv2[j];
  k21 = adja[jj] / 3;
  v21 = adja[jj] % 3;

  pt->v[i2]  = pt1->v[j];
  pt1->v[j2] = pt->v[i];

  /* update info */
  pt->tag[i] = pt1->tag[jj];
  pt->edg[i] = pt1->edg[jj];
  pt->base   = mesh->base;
  pt1->tag[j] = pt->tag[i1];
  pt1->edg[j] = pt->edg[i1];
  pt->tag[i1] = 0;
  pt->edg[i1] = 0;
  pt1->tag[jj] = 0;
  pt1->edg[jj] = 0;
  pt1->base    = mesh->base;

  /* update adjacent */
  mesh->adja[3*(k-1)+1+i]     = 3*k21+v21;
  if ( k21 )
    mesh->adja[3*(k21-1)+1+v21] = 3*k+i;
  mesh->adja[3*(k-1)+1+i1]    = 3*adj+jj;
  mesh->adja[3*(adj-1)+1+jj]  = 3*k+i1;
  if ( k11 )
    mesh->adja[3*(k11-1)+1+v11] = 3*adj+j;
  mesh->adja[3*(adj-1)+1+j]   = 3*k11+v11;

  return 1;
}



