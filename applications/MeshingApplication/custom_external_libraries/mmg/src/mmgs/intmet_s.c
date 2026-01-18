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
 * \file mmgs/intmet_s.c
 * \brief Metric interpolations.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmgs_private.h"

extern int8_t ddb;




/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param s interpolation parameter.
 * \param mr computed metric.
 * \return call to MMG5_interpreg_ani (thus, 0 if fail, 1 otherwise).
 *
 * Metric interpolation on edge \a i in elt \a it at
 * parameter \f$ 0 <= s0 <= 1 \f$ from \a p1 result is stored in \a mr. edge
 * \f$ p_1-p_2 \f$ must not be a ridge.
 *
 * */
int intregmet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,double s,double mr[6]) {
  MMG5_pTria     pt;

  pt  = &mesh->tria[k];

  return MMG5_interpreg_ani(mesh,met,pt,i,s,mr);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k triangle in which we interpole the metrics.
 * \param i edge along which we interpole the metrics.
 * \param ip index of point in which we compute the interpolated metric.
 * \param s parameter at which we compute the interpolation.
 * \return 1 if success, 0 otherwise.
 *
 * Linear interpolation of sizemap along edge i of tria k.
 *
 */
int intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,double s) {
  MMG5_pTria pt;
  MMG5_int   ip1,ip2;
  int8_t     i1,i2;

  pt  = &mesh->tria[k];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  met->m[ip] = s * (met->m[ip1] + met->m[ip2]);
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param ip global index of the new point in which we want to compute the metric.
 * \param s interpolation parameter (between 0 and 1).
 * \return 0 if fail, 1 otherwise.
 *
 * Interpolation of anisotropic sizemap at parameter \a s along edge \a i of elt
 * \a k for special storage of ridges metrics (after defsiz call).
 *
 */
int intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,double s) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  MMG5_pxPoint  go;
  double        *m;
  int           i1, i2;
  MMG5_int      ip1, ip2;

  pt  = &mesh->tria[k];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  m  = &met->m[6*ip];
  if ( pt->tag[i] & MG_GEO ) {
    ppt = &mesh->point[ip];
    assert(ppt->xp);
    go = &mesh->xpoint[ppt->xp];
    return MMG5_intridmet(mesh,met,ip1,ip2,s,go->n1,m);
  }
  else {
    return intregmet(mesh,met,k,i,s,m);
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param ip global index of the new point in which we want to compute the metric.
 * \param s interpolation parameter (between 0 and 1).
 * \return 0 if fail, 1 otherwise.
 *
 * Interpolation of anisotropic sizemap at parameter \a s along edge \a i of elt
 * \a k for classic storage of ridges metrics (before defsiz call).
 *
 */
int MMGS_intmet33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,double s) {
  MMG5_pTria    pt;
  double        *mr,*m,*n;
  int           i1, i2;
  MMG5_int      ip1, ip2;

  pt  = &mesh->tria[k];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  m   = &met->m[6*ip1];
  n   = &met->m[6*ip2];
  mr  = &met->m[6*ip];

  return MMG5_mmgIntmet33_ani(m,n,mr,s);
}
