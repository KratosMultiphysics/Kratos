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
 * \file mmgs/boulep_s.c
 * \brief Functions for ball of points computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmgs_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of tetra to start to compute the ball.
 * \param ip index of point in tetra \a start for which we want to compute
 * the ball.
 * \param list pointer to the computed ball of point.
 *
 * \return size of list if success, -size if overflow, 0 if cfg is non-manifold.
 *
 * Find all triangles sharing \a ip, \f$list[0] = start\f$ . Do not stop when
 * crossing ridge. Check whether resulting configuration is manifold.
 *
 */
int boulechknm(MMG5_pMesh mesh,MMG5_int start,int ip,MMG5_int *list) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  MMG5_int      *adja,k,iel,base;
  int           ilist;
  int8_t        i,i1,i2,ia,iq,voy;

  base = ++mesh->base;

  pt = &mesh->tria[start];
  ia = MMG5_iprv2[ip];
  iq = MMG5_inxt2[ip];
  if ( !MG_EOK(pt) )  return 0;
  ppt = &mesh->point[pt->v[ip]];
  if ( ppt->tag & MG_NOM )  return 0;
  ilist = 0;

  /* store neighbors */
  k = start;
  i = ip;
  do {
    if ( ilist > MMG5_TRIA_LMAX-2 )  return -ilist;
    list[ilist] = 3*k + i;
    ++ilist;

    pt = &mesh->tria[k];

    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];
    ppt = &mesh->point[pt->v[i1]];
    ppt->s = base;
    ppt = &mesh->point[pt->v[i2]];
    ppt->s = base;

    adja = &mesh->adja[3*(k-1)+1];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];
  }
  while ( k && k != start );

  /* check if boundary hit */
  if ( k <= 0 ) {
    k = start;
    i = ip;
    do {
      adja = &mesh->adja[3*(k-1)+1];
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      pt = &mesh->tria[k];
      ppt = &mesh->point[pt->v[i1]];
      ppt->s = base;
      ppt = &mesh->point[pt->v[i2]];
      ppt->s = base;

      k  = adja[i2] / 3;
      if ( k == 0 )  break;
      i  = adja[i2] % 3;
      i  = MMG5_iprv2[i];

      if ( ilist > MMG5_TRIA_LMAX-2 )  return -ilist;
      list[ilist] = 3*k + i;
      ilist++;
    }
    while ( k );
  }

  pt = &mesh->tria[start];
  i1 = MMG5_inxt2[ip];
  i2 = MMG5_iprv2[ip];
  ppt = &mesh->point[pt->v[i1]];
  ppt->s = 0;
  ppt = &mesh->point[pt->v[i2]];
  ppt->s = 0;

  adja = &mesh->adja[3*(start-1)+1];
  iel = adja[ia] / 3;
  voy = adja[ia] % 3;

  if( iel ) {
    pt = &mesh->tria[iel];
    ppt = &mesh->point[pt->v[voy]];
    ppt->s = 0;
  }

  /* check if a collapse may lead to a non-convex situation */
  k = start;
  i = iq;
  do {
    pt = &mesh->tria[k];

    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];
    ppt = &mesh->point[pt->v[i1]];
    if ( ppt->s == base ) return 0;
    ppt = &mesh->point[pt->v[i2]];
    if ( ppt->s == base ) return 0;

    adja = &mesh->adja[3*(k-1)+1];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];
  }
  while ( k && k != start );
  if( k > 0 ) return ilist;

  /* check if boundary hit */
  k = start;
  i = iq;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    pt = &mesh->tria[k];
    ppt = &mesh->point[pt->v[i1]];
    if ( ppt->s == base ) return 0;
    ppt = &mesh->point[pt->v[i2]];
    if ( ppt->s == base ) return 0;

    k  = adja[i2] / 3;
    if ( k == 0 )  break;
    i  = adja[i2] % 3;
    i  = MMG5_iprv2[i];
  }
  while ( k );

  return ilist;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of the starting triangle.
 * \param ip index of the looked ridge point.
 * \param il1 pointer to the first ball size.
 * \param l1 pointer to the first computed ball (associated to \a n1's
 * side).
 * \param il2 pointer to the second ball size.
 * \param l2 pointer to the second computed ball (associated to \a n2's
 * side).
 * \param global ip0 index of the first extremity of the ridge.
 * \param global ip1 index of the second extremity of the ridge.
 * \return 0 if fail, 1 otherwise.
 *
 * Computation of the two balls of a ridge point: the list \a l1 is associated
 * to normal \a n1's side. \a ip0 and \a ip1 are the indices of the 2 ending
 * point of the ridge. Both lists are returned enumerated in direct order.
 *
 */
int bouletrid(MMG5_pMesh mesh,MMG5_int start,MMG5_int ip,int *il1,MMG5_int *l1,int *il2,MMG5_int *l2,MMG5_int *ip0,MMG5_int *ip1) {
  MMG5_pTria   pt;
  MMG5_pPoint  ppt;
  MMG5_int     idp,k,kold,*adja,iel,*list1,*list2,aux;
  int          *ilist1,*ilist2;
  uint8_t      i,iold,i1,i2,ipn;
  double       *n1,*n2,nt[3],ps1,ps2;

  pt = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return 0;

  idp = pt->v[ip];
  ppt = &mesh->point[idp];
  assert( ppt->tag & MG_GEO );

  /* set pointers: first manifold is on side of triangle */
  if ( !MMG5_nortri(mesh,pt,nt) )  return 0;

  n1 = &(mesh->xpoint[ppt->xp].n1[0]);
  n2 = &(mesh->xpoint[ppt->xp].n2[0]);
  ps1 = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
  ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];

  if ( fabs(ps1) < fabs(ps2) ) {
    list1  = l2;
    list2  = l1;
    ilist1 = il2;
    ilist2 = il1;
  }
  else {
    list1  = l1;
    list2  = l2;
    ilist1 = il1;
    ilist2 = il2;
  }
  *ilist1 = 0;

  /* First ball, first side (via i1)*/
  k = start;
  i = ip;
  do {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];
    kold = k;
    iold = i;
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  // Remark: here the test k!=start is a security bound: theorically it is
  // useless but in case of bad edge tag, it ensure that the loop is not
  // infinite.
  while ( k && !(pt->tag[i1] & MG_GEO) && k != start );
  *ip0 = pt->v[i2];

  /* Store the needed elements to start back in the new area,
     and complete first ball, second side (via i2) */
  iel = k;
  ipn = i;
  k  = kold;
  i  = iold;
  do {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    if ( (*ilist1) > MMG5_TRIA_LMAX-2 )  return 0;
    list1[(*ilist1)] = 3*k+i;
    (*ilist1)++;
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];
    k = adja[i2] / 3;
    i = adja[i2] % 3;
    i = MMG5_iprv2[i];
  }
  while ( k && !(MG_GEO & pt->tag[i2]) );
  *ip1 = pt->v[i1];

  /* Invert the order in list1, so that it is enumerated in DIRECT order */
  for (k=0; k<(*ilist1) / 2;k++) {
    aux = list1[k];
    list1[k] = list1[*ilist1-1-k];
    list1[*ilist1-1-k] = aux;
  }
  /* At this point, either something has been stored in iel or an open boundary has been hit */
  *ilist2 = 0;
  if ( !iel )  return 1;

  /* Else, start back from the hit boundary, until another boundary is hit */
  k  = iel;
  i  = ipn;
  do {
    pt   = &mesh->tria[k];
    adja = &mesh->adja[3*(k-1)+1];
    if ( *ilist2 > MMG5_TRIA_LMAX-2 )  return 0;
    list2[*ilist2] = 3*k+i;
    (*ilist2)++;
    i1 = MMG5_inxt2[i];
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  while ( k && !(MG_GEO & pt->tag[i1]) );

  if ( !(MG_GEO & pt->tag[i1]) )
    return 0;
  else
    return 1;
}
