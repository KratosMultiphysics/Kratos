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
#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"

extern uint8_t ddb;

/**
 * \param mesh pointer to the mesh
 * \param met pointer to the metric
 * \param k triangle index
 * \param i local index of the edge that we want to test in the triangle \a k
 * \param list edge's shell (to fill)
 * \param typchk type eof check to perform.
 *
 * \return 1 if we must collapse, 0 otherwise
 *
 * Check whether the validity and the geometry of the mesh are
 * preserved when collapsing edge i (p1->p2)
 *
 */
int MMG2D_chkcol(MMG5_pMesh mesh, MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int *list,int8_t typchk) {
  MMG5_pTria     pt0,pt,pt1,pt2;
  MMG5_pPoint    ppt,p2;
  double         lon,len,calold,calnew,caltmp;
  MMG5_int       ip1,ip2,ipb,l,ll,lj,jel,kel,*adja;
  int            ilist;
  uint8_t        i1,i2,j,jj,j2,voy,open;

  pt0 = &mesh->tria[0];
  pt = &mesh->tria[k];

  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];

  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  open = 1;
  adja = &mesh->adja[3*(k-1)+1];

  calold = calnew = DBL_MAX;

  /* If typchk == 2, avoid recreating long edges */
  lon = 0.0;
  if ( typchk == 2 && met->m ) {
    lon = MMG2D_lencurv(mesh,met,ip1,ip2);
    lon = MG_MAX(2.-lon,1.6);
  }

  /* Avoid collapsing a boundary point over a regular one (leads to boundary degeneration) */
  if ( MG_EDG(mesh->point[ip1].tag) && !MG_EDG(mesh->point[ip2].tag) ) return 0;

  /* Avoid subdivising one ref. component consisting of only one layer of elements */
  if ( MG_EDG(mesh->point[ip1].tag) && (pt->tag[i1]||pt->tag[i2]) ) return 0;

  /* Avoid closing one ref. component consisting of only one element (maybe
   * redondant with previous test)*/
  if ( MG_EDG(pt->tag[i1]) && MG_EDG(pt->tag[i2]) ) return 0;

  jel = adja[i] / 3;
  if ( jel ) {
    open = 0;
    pt1 = &mesh->tria[jel];
    j = adja[i] % 3;
    jj = MMG5_inxt2[j];
    j2 = MMG5_iprv2[j];
    if ( MG_EDG(pt1->tag[jj]) && MG_EDG(pt1->tag[j2]) ) return 0;
  }

  /* collect all triangles around vertex i1
   simplification w.r.to the surface case: no need to use spacial ball boulechknm (impossible situation in 2d) */
  int8_t dummy;
  ilist = MMG5_boulet(mesh,k,i1,list,0,&dummy);
  if ( ilist <= 0 )
    return 0;

  /* Avoid collapsing a "void triangle" */
  if ( open ) {
    jel = list[ilist-1] / 3;
    pt1 = &mesh->tria[jel];
    j = list[ilist-1] % 3;
    j2 = MMG5_iprv2[j];
    ipb = pt1->v[j2];

    /* Travel the ball of ip2 until a boundary is met */
    jel = k;
    pt1 = &mesh->tria[jel];
    j = i1;

    while ( jel && !MG_EDG(pt1->tag[j]) ) {
      adja = &mesh->adja[3*(jel-1)+1];
      jel = adja[j] / 3;
      pt1 = &mesh->tria[jel];
      j = MMG5_inxt2[adja[j] % 3];
    }
    j = MMG5_iprv2[j];

    if ( ipb == pt1->v[j] ) return 0;
  }

  /* Avoid creating edges with two BDY endpoints, which are not themselves BDY */
  if ( mesh->info.fem ) {
    p2 = &mesh->point[ip2];
    if ( p2->tag & MG_BDY ) {
      /* Travel all edges in the ball but for the first and last */
      for (l=0; l<ilist-1; l++) {
        jel = list[l] / 3;
        j = list[l] % 3;
        jj = MMG5_inxt2[j];
        j2 = MMG5_iprv2[j];
        pt1 = &mesh->tria[jel];
        p2 = &mesh->point[pt1->v[j2]];
        if ( (p2->tag & MG_BDY) && !(pt1->tag[jj] & MG_BDY) ) return 0;
      }
    }
  }

  if ( ilist > 3 || ( ilist == 3 && open ) ) {      // ADD : Second test should work too
    /* Avoid a collapse that would close one ref. component that consists only of one element */
    /* Useless test I believe */
    /*if ( MG_EDG(pt->tag[i2]) ) {
      jel = list[1] / 3;
      if ( ! jel ) return 0;
      pt1 = &mesh->tria[jel];
      if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;
    }*/

    /* Travel the ball of i1 (but for the two elements 0 and ilist-1 (the last one in the case of
     a closed ball) which will disappear */
    for (l=1; l<ilist-1+open; l++) {
      jel = list[l] / 3;
      j = list[l] % 3;
      jj = MMG5_inxt2[j];
      j2 = MMG5_iprv2[j];
      pt1 = &mesh->tria[jel];

      memcpy(pt0,pt1,sizeof(MMG5_Tria));
      pt0->v[j] = ip2;

      /* Check length to avoid recreating long elements */
      if ( typchk == 2 && met->m && !MG_EDG(mesh->point[ip2].tag) ) {
        ip1 = pt1->v[j2];
        len = MMG2D_lencurv(mesh,met,ip1,ip2);
        if ( len > lon )  return 0;
      }

      /* Check that the newly created triangles will not have to be split again */
      if ( l == 1 ) {
        pt0->tag[j2] |= pt->tag[i1];
      }
      else if ( l == ilist-2 && !open ) {
        ll = list[ilist-1+open] / 3;

        if ( ll > mesh->nt )  return 0;
        lj = list[ilist-1+open] % 3;
        pt0->tag[jj] |= mesh->tria[ll].tag[lj];
      }

      if ( typchk == 1 )
        if ( MMG2D_chkedg(mesh,0) )  return 0;


      /* Check quality and volume inversion */
      if ( typchk == 2 && met->m && met->size == 3 )
        caltmp = MMG2D_ALPHAD*MMG2D_caltri_ani(mesh,met,pt1);
      else
        caltmp = MMG2D_ALPHAD*MMG2D_caltri_iso(mesh,NULL,pt1);

      calold = MG_MIN(calold,caltmp);

      if ( typchk == 2 && met->m && met->size == 3 )
        caltmp = MMG2D_ALPHAD*MMG2D_caltri_ani(mesh,met,pt0);
      else
        caltmp = MMG2D_ALPHAD*MMG2D_caltri_iso(mesh,NULL,pt0);

      if ( caltmp < MMG2D_NULKAL )  return 0;
      calnew = MG_MIN(calnew,caltmp);
      if ( calold < MMG2D_NULKAL && calnew <= calold )  return 0;
      else if ( calnew < MMG2D_NULKAL || calnew < 0.001*calold )  return 0;
    }
  }

  /* Specific test: no collapse if any interior edge is EDG */
  else if ( ilist == 3 ) {
    ppt = &mesh->point[pt->v[i1]];
    if ( MG_SIN(ppt->tag) )  return 0;
    else if ( MG_EDG(pt->tag[i2]) && !MG_EDG(pt->tag[i]) )  return 0;
    else if ( !MG_EDG(pt->tag[i2]) && MG_EDG(pt->tag[i]) )  return 0; // What is the use of that?
    else if ( MG_EDG(pt->tag[i2]) && MG_EDG(pt->tag[i]) && MG_EDG(pt->tag[i1]) )  return 0;

    /* Check geometric approximation */
    jel = list[1] / 3;
    j   = list[1] % 3;
    jj  = MMG5_inxt2[j];
    j2  = MMG5_iprv2[j];
    pt0 = &mesh->tria[0];
    pt1 = &mesh->tria[jel];
    memcpy(pt0,pt1,sizeof(MMG5_Tria));
    pt0->v[j] = ip2;

    jel = list[2] / 3;
    j   = list[2] % 3;
    pt1 = &mesh->tria[jel];
    pt0->tag[jj] |= pt1->tag[j];
    pt0->tag[j2] |= pt1->tag[MMG5_inxt2[j]];

    if ( typchk == 1 )
      if ( MMG2D_chkedg(mesh,0) )  return 0;

  }

  /* Particular case when there are two triangles in the ball of the collapsed point ip1 */
  else {
    assert ( ilist == 2 );  // Not necessarily! Check for that case too!
    if ( ilist !=2 ) return 0;
    if ( !open )  return 0;

    jel = list[1] / 3;
    j   = list[1] % 3;

    /* Topological test : avoid creating twice the same triangle */
    adja = &mesh->adja[3*(jel-1)+1];
    kel = adja[j] / 3;
    voy = adja[j] % 3;
    pt2 = &mesh->tria[kel];
    if ( pt2->v[voy] == ip2 ) return 0;

    jj  = MMG5_inxt2[j];
    pt1 = &mesh->tria[jel];
    if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;
    else if ( !(pt1->tag[jj] & MG_GEO) )  return 0;

    /* Check quality and geometric approximation: elements with two trias in the
     * ball should be removed */
    pt0 = &mesh->tria[0];
    memcpy(pt0,pt1,sizeof(MMG5_Tria));
    pt0->v[j] = ip2;

    calold = MMG2D_ALPHAD*MMG2D_caltri_iso(mesh,NULL,pt1);
    calnew = MMG2D_ALPHAD*MMG2D_caltri_iso(mesh,NULL,pt0);
    if ( calnew < MMG2D_NULKAL )  return 0;
    if ( calold < MMG2D_NULKAL && calnew <= calold )  return 0;

    if ( typchk == 1 )
      if ( MMG2D_chkedg(mesh,0) )  return 0;
  }

  return ilist;
}

/* Perform effective collapse of edge i in tria k, i1->i2 */
int MMG2D_colver(MMG5_pMesh mesh,int ilist,MMG5_int *list) {
  MMG5_pTria   pt,pt1,pt2;
  MMG5_int     iel,jel,ip1,ip2,k,kel,*adja;
  uint8_t      i,j,jj,i1,i2,open;

  iel = list[0] / 3;
  i1 =   list[0] % 3;
  i = MMG5_iprv2[i1];
  i2 = MMG5_inxt2[i1];
  pt = &mesh->tria[iel];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  /* Check for open ball */
  adja = &mesh->adja[3*(iel-1)+1];
  open = adja[i] == 0;

  /* Update of the vertex ip1 -> ip2 */
  mesh->point[ip2].tag |= mesh->point[ip1].tag;

  for (k=1; k<ilist-1+open; k++) {
    jel = list[k] / 3;
    jj  = list[k] % 3;
    pt1 = &mesh->tria[jel];
    pt1->v[jj] = ip2;
    pt1->base  = mesh->base;    // What is the use of updating the base field?
  }

  /* Update adjacent with the first element (which necessarily exists) */
  jel = list[1] / 3;
  jj  = list[1] % 3;
  j   = MMG5_iprv2[jj];
  pt1 = &mesh->tria[jel];

  pt1->tag[j] |= pt->tag[i1];
  pt1->edg[j] = MG_MAX(pt1->edg[j],pt->edg[i1]);

  if ( adja[i1] ) {
    kel = adja[i1] / 3;
    k   = adja[i1] % 3;
    mesh->adja[3*(kel-1)+1+k] = 3*jel+j;
    mesh->adja[3*(jel-1)+1+j] = 3*kel+k;
    pt2 = &mesh->tria[kel];
    pt2->tag[k] |= pt1->tag[j];
    pt2->edg[k] = MG_MAX(pt1->edg[j],pt2->edg[k]);
  }
  else
    mesh->adja[3*(jel-1)+1+j] = 0;

  /* When the boundary is not open, update the adjacency relation with the second element */
  if ( !open ) {
    iel = list[ilist-1] / 3;
    i1 = list[ilist-1] % 3;
    pt = &mesh->tria[iel];

    jel = list[ilist-2] / 3;
    jj  = list[ilist-2] % 3;
    j   = MMG5_inxt2[jj];
    pt1 = &mesh->tria[jel];
    pt1->tag[j] |= pt->tag[i1];
    pt1->edg[j] = MG_MAX(pt->edg[i1],pt1->edg[j]);
    adja = &mesh->adja[3*(iel-1)+1];

    if ( adja[i1] ) {
      kel = adja[i1] / 3;
      k   = adja[i1] % 3;
      mesh->adja[3*(kel-1)+1+k] = 3*jel + j;
      mesh->adja[3*(jel-1)+1+j] = 3*kel + k;
      pt2 = &mesh->tria[kel];
      pt2->tag[k] |= pt1->tag[j];
      pt2->edg[k] = MG_MAX(pt1->edg[j],pt2->edg[k]);
    }
    else
      mesh->adja[3*(jel-1)+1+j] = 0;
  }

  /* Delete removed point and elements */
  MMG2D_delPt(mesh,ip1);
  MMG2D_delElt(mesh,list[0] / 3);
  if ( !open )  MMG2D_delElt(mesh,list[ilist-1] / 3);

  return 1;
}

/* Perform effective collapse of edge i in tria k, i1->i2
   in the particular case where only three elements are in the ball of i */
int MMG2D_colver3(MMG5_pMesh mesh,MMG5_int *list) {
  MMG5_pTria  pt,pt1,pt2;
  MMG5_int    iel,jel,kel,mel,ip,*adja;
  uint8_t     i,i1,j,j1,j2,k,m;

  /* Update of the new point for triangle list[0] */
  iel = list[0] / 3;
  i = list[0] % 3;
  i1 = MMG5_inxt2[i];
  pt = &mesh->tria[iel];
  ip = pt->v[i];
  jel = list[1] / 3;
  j   = list[1] % 3;
  j1  = MMG5_inxt2[j];
  j2  = MMG5_iprv2[j];
  pt1 = &mesh->tria[jel];

  kel = list[2] / 3;
  k   = list[2] % 3;
  pt2 = &mesh->tria[kel];

  /* update info */
  pt1->v[j]     = pt->v[i1];
  pt1->tag[j1] |= pt2->tag[k];
  pt1->edg[j1]  = MG_MAX(pt1->edg[j1],pt2->edg[k]);
  pt1->tag[j2] |= pt->tag[i];
  pt1->edg[j2]  = MG_MAX(pt1->edg[j2],pt->edg[i]);
  pt1->base     = mesh->base;

  /* Update adjacency relations */
  adja = &mesh->adja[3*(jel-1)+1];
  adja[j1] = mesh->adja[3*(kel-1)+1+k];
  adja[j2] = mesh->adja[3*(iel-1)+1+i];

  mel  = adja[j2] / 3;
  if ( mel ) {
    m    = adja[j2] % 3;
    pt   = &mesh->tria[mel];
    pt->tag[m]  = pt1->tag[j2];
    pt->edg[m]  = pt1->edg[j2];
    mesh->adja[3*(mel-1)+1+m] = 3*jel + j2;
  }

  mel = adja[j1] / 3;
  if ( mel ) {
    m    = adja[j1] % 3;
    pt   = &mesh->tria[mel];
    pt->tag[m]  = pt1->tag[j1];
    pt->edg[m]  = pt1->edg[j1];
    mesh->adja[3*(mel-1)+1+m] = 3*jel + j1;
  }

  /* remove vertex + elements */
  MMG2D_delPt(mesh,ip);
  MMG2D_delElt(mesh,iel);
  MMG2D_delElt(mesh,kel);

  return 1;
}

/* Perform effective collapse of edge i in tria k, i1->i2
 in the particular case where only two elements are in the ball of i */
int MMG2D_colver2(MMG5_pMesh mesh,MMG5_int *list) {
  MMG5_pTria   pt,pt1;
  MMG5_int     *adja,iel,jel,kel,ip1,ip2;
  int8_t       i1,i2,jj,j2,k;

  /* update of new point for triangle list[0] */
  iel = list[0] / 3;
  i1  = list[0] % 3;
  i2  = MMG5_inxt2[i1];
  pt  = &mesh->tria[iel];
  ip1  = pt->v[i1];
  ip2  = pt->v[i2];

  jel = list[1] / 3;
  j2  = list[1] % 3;
  jj  = MMG5_iprv2[j2];
  pt1 = &mesh->tria[jel];

  /* update info */
  pt1->v[j2] = ip2;
  pt1->tag[jj] |= pt->tag[i1];
  pt1->edg[jj] = pt->edg[i1];
  pt1->base = mesh->base;

  /* update neighbours of new triangle */
  adja = &mesh->adja[3*(jel-1)+1];
  adja[jj] = mesh->adja[3*(iel-1)+1+i1];
  adja = &mesh->adja[3*(iel-1)+1];
  kel  = adja[i1] / 3;
  k    = adja[i1] % 3;
  if ( kel )
    mesh->adja[3*(kel-1)+1+k] = 3*jel + jj;

  /* remove vertex + element */
  MMG2D_delPt(mesh,ip1);
  MMG2D_delElt(mesh,iel);

  return 1;
}
