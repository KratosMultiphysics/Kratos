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
 * \file common/boulep.c
 * \brief Functions for ball of points computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon.h"

extern MMG5_Info  info;

/**
 * \param mesh pointer toward the mesh structure.
 * \param adjt pointer toward the table of triangle adjacency.
 * \param start index of triangle where we start to work.
 * \param ip index of vertex where the normal is computed.
 * \param nn pointer toward the computed tangent.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute average normal of triangles sharing P without crossing ridge.
 *
 */
int _MMG5_boulen(MMG5_pMesh mesh,int *adjt,int start,int ip,double *nn) {
  MMG5_pTria    pt;
  double        n[3],dd;
  int           *adja,k;
  char          i,i1,i2;

  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return(0);
  nn[0] = nn[1] = nn[2] = 0.0;

  /* store neighbors */
  k  = start;
  i  = ip;
  i1 = _MMG5_inxt2[i];
  do {
    pt = &mesh->tria[k];
    _MMG5_nortri(mesh,pt,n);
    nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];

    if ( pt->tag[i1] & MG_GEO ) {
      k = 0;
      break;
    }
    adja = &adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i2 = adja[i1] % 3;
    i1 = _MMG5_iprv2[i2];
  }
  while ( k && k != start );

  if ( k == 0 ) {
    k  = start;
    i  = ip;
    i2 = _MMG5_iprv2[i];
    pt = &mesh->tria[k];
    do {
      if ( pt->tag[i2] & MG_GEO )  break;

      adja = &adjt[3*(k-1)+1];
      k  = adja[i2] / 3;
      if ( k == 0 )  break;
      i1 = adja[i2] % 3;
      i2 = _MMG5_inxt2[i1];
      pt = &mesh->tria[k];

      _MMG5_nortri(mesh,pt,n);

      nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];
    }
    while ( k && k != start );
  }

  /* normalize */
  dd = nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2];
  if ( dd > _MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    nn[0] *= dd;
    nn[1] *= dd;
    nn[2] *= dd;
    return(1);
  }

  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param adjt pointer toward the table of triangle adjacency.
 * \param start index of triangle where we start to work.
 * \param ip index of vertex where the tangent is computed.
 * \param tt pointer toward the computed tangent.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute the tangent to the curve at point \a ip.
 *
 */
int _MMG5_boulec(MMG5_pMesh mesh,int *adjt,int start,int ip,double *tt) {
  MMG5_pTria    pt;
  MMG5_pPoint   p0,p1,p2;
  double        dd;
  int           *adja,k;
  char          i,i1,i2;

  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) )       return(0);
  p0 = &mesh->point[pt->v[ip]];
  if ( !MG_EDG(p0->tag) )  return(0);

  /* check other triangle vertices */
  k  = start;
  i  = ip;
  i1 = _MMG5_inxt2[i];
  i2 = _MMG5_iprv2[i];
  p1 = p2 = 0;
  do {
    pt = &mesh->tria[k];
    if ( MG_EDG(pt->tag[i1]) ) {
      p1 = &mesh->point[pt->v[i2]];
      k  = 0;
      break;
    }
    adja = &adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i2 = adja[i1] % 3;
    i1 = _MMG5_iprv2[i2];
  }
  while ( k && k != start );

  /* check if open boundary hit */
  if ( k == 0 ) {
    k  = start;
    i  = ip;
    i1 = _MMG5_inxt2[i];
    i2 = _MMG5_iprv2[i];
    do {
      pt = &mesh->tria[k];
      if ( MG_EDG(pt->tag[i2]) ) {
        p2 = &mesh->point[pt->v[i1]];
        break;
      }
      adja = &adjt[3*(k-1)+1];
      k  = adja[i2] / 3;
      i1 = adja[i2] % 3;
      i2 = _MMG5_inxt2[i1];
    }
    while ( k );
  }

  if ( !p1 || !p2 )
    return(0);
  else if ( p1 == p2 )
    p2 = p0;

  /* tangent approx */
  tt[0] = p2->c[0] - p1->c[0];
  tt[1] = p2->c[1] - p1->c[1];
  tt[2] = p2->c[2] - p1->c[2];
  dd = tt[0]*tt[0] + tt[1]*tt[1] + tt[2]*tt[2];
  if ( dd > _MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    tt[0] *= dd;
    tt[1] *= dd;
    tt[2] *= dd;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param adjt pointer toward the table of triangle adjacency.
 * \param start index of triangle where we start to work.
 * \param ip index of vertex on which we work.
 * \param list pointer toward the computed list of edges incident to \a ip.
 * \param ng pointer toward the number of ridges.
 * \param nr pointer toward the number of reference edges.
 * \param lmax maxmum size for the ball of the point \a ip.
 * \return The number of edges incident to the vertex \a ip.
 *
 * Store edges and count the number of ridges and reference edges incident to
 * the vertex \a ip.
 *
 */
int _MMG5_bouler(MMG5_pMesh mesh,int *adjt,int start,int ip,
                 int *list,int *ng,int *nr,int lmax) {
  MMG5_pTria    pt;
  int           *adja,k,ns;
  char          i,i1,i2;

  pt  = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return(0);

  /* check other triangle vertices */
  k  = start;
  i  = ip;
  *ng = *nr = ns = 0;
  do {
    i1 = _MMG5_inxt2[i];
    if ( MG_EDG(pt->tag[i1])) {
      i2 = _MMG5_iprv2[i];
      if ( pt->tag[i1] & MG_GEO )
        *ng = *ng + 1;
      else if ( pt->tag[i1] & MG_REF )
        *nr = *nr + 1;
      ns++;
      list[ns] = pt->v[i2];
      if ( ns > lmax-2 )  return(-ns);
    }
    adja = &adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = _MMG5_inxt2[i];
    pt = &mesh->tria[k];
  }
  while ( k && k != start );

  /* reverse loop */
  if ( k != start ) {
    k = start;
    i = ip;
    do {
      pt = &mesh->tria[k];
      i2 = _MMG5_iprv2[i];
      if ( MG_EDG(pt->tag[i2]) ) {
        i1 = _MMG5_inxt2[i];
        if ( pt->tag[i2] & MG_GEO )
          *ng = *ng + 1;
        else if ( pt->tag[i2] & MG_REF )
          *nr = *nr + 1;
        ns++;
        list[ns] = pt->v[i1];
        if ( ns > lmax-2 )  return(-ns);
      }
      adja = &adjt[3*(k-1)+1];
      k = adja[i2] / 3;
      i = adja[i2] % 3;
      i = _MMG5_iprv2[i];
    }
    while ( k && k != start );
  }
  return(ns);
}
