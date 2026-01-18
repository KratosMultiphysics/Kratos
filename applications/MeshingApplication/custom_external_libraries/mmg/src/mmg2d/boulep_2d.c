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

#include "libmmg2d.h"
#include "libmmg2d_private.h"


static uint8_t inxt[3]  = {1,2,0};
static uint8_t iprev[3] = {2,0,1};

/* Function boulep should be eventually removed (after all its occurences are removed) */
/* Find all triangles sharing P
   in:  ifirst    : triangle containing p
   iploc     : index of p in start
   out: list  : list of triangles */
int MMG2D_boulep(MMG5_pMesh mesh, MMG5_int ifirst, int iploc, MMG5_int * list) {
  MMG5_pTria  pt;
  MMG5_pPoint ppt;
  MMG5_int    ip,iel,*adja,iadr;
  int         voy,ilist,i;

  if ( ifirst < 1 ) return 0;
  pt = &mesh->tria[ifirst];
  if ( !MG_EOK(pt) ) return 0;
  ip = pt->v[iploc];
  ppt = &mesh->point[ip];
  if ( !MG_VOK(ppt) ) return 0;

  /* init list */
  ilist       = 1;
  list[ilist] = 3*ifirst + iploc;

  iadr = 3*(ifirst-1) + 1;
  adja = &mesh->adja[iadr];
  iel  = adja[inxt[iploc]]/3;
  voy  = adja[inxt[iploc]]%3;
  i    = inxt[voy];

  while ( iel && (iel != ifirst) && mesh->tria[iel].v[0]){
    if(ilist==MMG2D_LMAX-1) return 0;
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
      if(ilist==MMG2D_LMAX-1) return 0;
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

  return ilist;
}

/**
 * \param mesh pointer to the mesh
 * \param start index of the first tria
 * \param ip local index of the point on which we work in tria \a start
 * \param pleft index of the left point of ip along the geom edge
 * (local index of the point inside a triangle)
 * \param pright index of the right point of ip along the geom edge
 * (local index of the point inside a triangle)
 * \param nn computed normal at ip
 *
 * \return 0 if fail, 1 if success
 *
 * Travel the ball of point \a ip in triangle \a start, which is assumed to lie
 * either on the external or on an internal boundary of the mesh, and return the
 * normal vector convention: the normal vector is oriented from the half ball it
 * starts with towards its exterior return pright = 3*kk+ii, where kk = last
 * triangle in the first travel, and ii = local index of ip in kk pleft =
 * 3*kk+ii, where kk = last triangle in the second travel, and ii = local index
 * of ip in kk.
 *
 */
int MMG2D_boulen(MMG5_pMesh mesh, MMG5_int start,int8_t ip, MMG5_int *pleft, MMG5_int *pright, double *nn) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            ux,uy,dd,n1[2],n2[2];
  MMG5_int          *adja,k,kk,refs;
  int8_t            notedg;
  int8_t            i,ii,i1,i2;

  /* First travel of the ball of ip; initialization */
  kk = start;
  ii = MMG5_iprv2[ip];
  refs = mesh->tria[start].ref;

  do {
    k = kk;
    i = MMG5_iprv2[ii];
    adja = &mesh->adja[3*(k-1)+1];
    kk = adja[i] / 3;
    ii = adja[i] % 3;

    notedg = mesh->info.opnbdy ?
      (!mesh->tria[k].tag[i]) : (mesh->tria[kk].ref == refs);
  }
  while ( kk && (kk != start) && notedg );

  if ( kk == start ) {
    fprintf(stderr,"  ## Error: %s: Unable to find a boundary edge in"
            " the ball of point %" MMG5_PRId ".\n",__func__,MMG2D_indPt(mesh,mesh->tria[start].v[ip]));
    return 0;
  }

  /* Calculation of the first normal vector */
  pt = &mesh->tria[k];
  i1 = MMG5_iprv2[i];
  i2 = MMG5_inxt2[i];

  p1 = &mesh->point[pt->v[i1]];
  p2 = &mesh->point[pt->v[i2]];
  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  dd = ux*ux + uy*uy;

  if ( dd < MMG5_EPSD ) {
    fprintf(stderr,"\n  ## Error: %s: Null edge"
            " length (%e).\n",__func__,dd);
    return 0;
  }

  dd = 1.0 / sqrt(dd);
  n1[0] = -uy*dd;
  n1[1] = ux*dd;

  *pright = 3*k+i1;

  /* Second travel */
  kk = start;
  ii = MMG5_inxt2[ip];

  do {
    k = kk;
    i = MMG5_inxt2[ii];
    adja = &mesh->adja[3*(k-1)+1];
    kk = adja[i] / 3;
    ii = adja[i] % 3;

    notedg = mesh->info.opnbdy ?
      !mesh->tria[k].tag[i] : mesh->tria[kk].ref == refs;
  }
  while ( kk && (kk != start) && notedg );

  /* Calculation of the second normal vector */
  pt = &mesh->tria[k];
  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];

  p1 = &mesh->point[pt->v[i1]];
  p2 = &mesh->point[pt->v[i2]];
  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  dd = ux*ux + uy*uy;

  if ( dd < MMG5_EPSD ) {
    fprintf(stderr,"\n  ## Error: %s: Null edge length"
            " (%e).\n",__func__,dd);
    return 0;
  }

  dd = 1.0 / sqrt(dd);
  n2[0] = uy*dd;
  n2[1] = -ux*dd;

  *pleft = 3*k+i1;

  /* Averaging of normal vectors */
  nn[0] = n1[0] + n2[0];
  nn[1] = n1[1] + n2[1];
  dd = nn[0]*nn[0] + nn[1]*nn[1];
  if ( dd > MMG5_EPSD ){
    dd = 1.0 / sqrt(dd);
    nn[0] *= dd;
    nn[1] *= dd;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of triangle to start.
 * \param ip index of point for wich we compute the ball.
 * \return 1 if success, 0 if fail.
 *
 * Find the two endpoints of the boundary curves joining ip, fill \a ip1 and
 * \a ip2 with their indices and list neighbouring triangles
 *
 */
int MMG2D_bouleendp(MMG5_pMesh mesh,MMG5_int start,int8_t ip,MMG5_int *ip1,MMG5_int *ip2,MMG5_int *list) {
  MMG5_pTria    pt;
  MMG5_int      *adja,k;
  int8_t        i,i1,i2;
  static int8_t mmgWarn0=0;
  int           ilist;

  *ip1 = 0;
  *ip2 = 0;
  if ( start < 1 ) return 0;
  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) ) return 0;

  /* Travel elements in the forward sense */
  ilist= 0;
  k = start;
  i = ip;
  do {

    if ( ilist > MMG2D_LMAX-2 )  return -ilist;
    list[ilist] = k;
    ++ilist;

    pt = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    if ( MG_EDG(pt->tag[i1]) ) {
      if ( *ip1 == 0 ) *ip1 = pt->v[i2];
      else {
        if ( *ip2 != 0 && *ip1 != pt->v[i2] && *ip2 != pt->v[i2] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 non singular"
                    " point at the intersection of 3 edges.\n",
                    __func__);
          }
          return 0;
        }
        if ( *ip1 != pt->v[i2] ) *ip2 = pt->v[i2];
      }
    }

    if ( MG_EDG(pt->tag[i2]) ) {
      if ( *ip1 == 0 )
        *ip1 = pt->v[i1];
      else {
        if ( *ip2 != 0 && *ip1 != pt->v[i1] && *ip2 != pt->v[i1] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 non singular"
                    " point at the intersection of 3 edges.\n",
                    __func__);
          }
          return 0;
        }
        if ( *ip1 != pt->v[i1] ) *ip2 = pt->v[i1];
      }
    }

    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];

  }
  while ( k && k != start );
  if ( k > 0 ) return ilist;

  /* Travel the ball in the reverse sense when a boundary is hit, starting from the neighbor of k */
  k = start;
  i = ip;
  adja = &mesh->adja[3*(k-1)+1];
  i2 = MMG5_iprv2[i];
  k = adja[i2] / 3;
  i = adja[i2] % 3;
  i = MMG5_iprv2[i];

  if ( !k ) return ilist;

  do {

    if ( ilist > MMG2D_LMAX-2 )  return -ilist;
    list[ilist] = k ;
    ++ilist;

    pt = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    if ( MG_EDG(pt->tag[i1]) ) {
      if ( *ip1 == 0 )
        *ip1 = pt->v[i2];
      else {
        if ( *ip2 != 0 && *ip1 != pt->v[i2] && *ip2 != pt->v[i2] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 non singular"
                    " point at the intersection of 3 edges.\n",
                    __func__);
          }
          return 0;
        }
        if ( *ip1 != pt->v[i2] ) *ip2 = pt->v[i2];
      }
    }

    if ( MG_EDG(pt->tag[i2]) ) {
      if ( *ip1 == 0 ) *ip1 = pt->v[i1];
      else {
        if ( *ip2 != 0 && *ip1 != pt->v[i1] && *ip2 != pt->v[i1] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 non singular"
                    " point at the intersection of 3 edges.\n",
                    __func__);
          }
          return 0;
        }
        if ( *ip1 != pt->v[i1] ) *ip2 = pt->v[i1];
      }
    }

    k  = adja[i2] / 3;
    if ( k == 0 )  break;
    i  = adja[i2] % 3;
    i  = MMG5_iprv2[i];

  }
  while ( k );

  return ilist;
}

