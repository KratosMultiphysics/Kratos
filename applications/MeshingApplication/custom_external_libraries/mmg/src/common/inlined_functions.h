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
 * \file common/inlined_functions.h
 * \brief inlined Functions
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon.h"

#ifndef _INLINED_FUNC_H
#define _INLINED_FUNC_H

/**
 * \param mesh pointer toward the mesh structure.
 * \param np0 index of edge's extremity.
 * \param np1 index of edge's extremity.
 * \param m0 metric at point np0.
 * \param m1 metric at point np1.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of surface edge \f$[np0;np1]\f$ according to the prescribed
 * aniso metrics \a m0 and \a m1.
 *
 */
static inline
double _MMG5_lenEdg(MMG5_pMesh mesh,int np0,int np1,
                    double *m0,double *m1,char isedg) {
  MMG5_pPoint   p0,p1;
  double        gammaprim0[3],gammaprim1[3],t[3],*n1,*n2,ux,uy,uz,ps1,ps2,l0,l1;

  p0 = &mesh->point[np0];
  p1 = &mesh->point[np1];

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  /* computation of the two tangent vectors to the underlying curve of [i0i1] */
  if ( MG_SIN(p0->tag) || (MG_NOM & p0->tag) ) {
    gammaprim0[0] = ux;
    gammaprim0[1] = uy;
    gammaprim0[2] = uz;
  }
  else if ( isedg ) {
    memcpy(t,p0->n,3*sizeof(double));
    ps1 = ux*t[0] + uy*t[1] + uz*t[2];
    gammaprim0[0] = ps1*t[0];
    gammaprim0[1] = ps1*t[1];
    gammaprim0[2] = ps1*t[2];
  }
  else {
    if ( MG_GEO & p0->tag ) {
      //assert(p0->xp);
      n1 = &mesh->xpoint[p0->xp].n1[0];
      n2 = &mesh->xpoint[p0->xp].n2[0];
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
      ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

      if ( fabs(ps2) < fabs(ps1) ) {
        n1  = &mesh->xpoint[p0->xp].n2[0];
        ps1 = ps2;
      }
    }
    else if ( MG_REF & p0->tag || MG_BDY & p0->tag ) {
      // ( MG_BDY  & p0->tag ) => mmg3d
      n1  = &mesh->xpoint[p0->xp].n1[0];
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    }
    else {
      // we come from mmgs because in mmg3d the boundary points are tagged
      // MG_BDY.
      n1  = &(p0->n[0]);
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    }
    gammaprim0[0] = ux - ps1*n1[0];
    gammaprim0[1] = uy - ps1*n1[1];
    gammaprim0[2] = uz - ps1*n1[2];
  }

  if ( MG_SIN(p1->tag) || (MG_NOM & p1->tag) ) {
    gammaprim1[0] = -ux;
    gammaprim1[1] = -uy;
    gammaprim1[2] = -uz;
  }
  else if ( isedg ) {
    memcpy(t,p1->n,3*sizeof(double));
    ps1 = -ux*t[0] - uy*t[1] - uz*t[2];
    gammaprim1[0] = ps1*t[0];
    gammaprim1[1] = ps1*t[1];
    gammaprim1[2] = ps1*t[2];
  }
  else {
    if ( MG_GEO & p1->tag ) {
      n1 = &mesh->xpoint[p1->xp].n1[0];
      n2 = &mesh->xpoint[p1->xp].n2[0];
      ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
      ps2 = -ux*n2[0] - uy*n2[1] - uz*n2[2];

      if ( fabs(ps2) < fabs(ps1) ) {
        n1  = &mesh->xpoint[p1->xp].n2[0];
        ps1 = ps2;
      }
    }
    else if ( MG_REF & p1->tag || MG_BDY & p1->tag ) {
      // ( MG_BDY  & p1->tag ) => mmg3d )
      n1  = &mesh->xpoint[p1->xp].n1[0];
      ps1 = - ux*n1[0] - uy*n1[1] - uz*n1[2];
    }
    else {
      // we come from mmgs because in mmg3d the boundary points are tagged
      // MG_BDY.
      n1  = &(p1->n[0]);
      ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
    }
    gammaprim1[0] = - ux - ps1*n1[0];
    gammaprim1[1] = - uy - ps1*n1[1];
    gammaprim1[2] = - uz - ps1*n1[2];
  }

  /* computation of the length of the two tangent vectors in their respective tangent plane */
  /* l_ab = int_a^b sqrt(m_ij d_t x_i(t) d_t x_j(t) ) : evaluated by a 2-point quadrature method. */
  l0 = m0[0]*gammaprim0[0]*gammaprim0[0] + m0[3]*gammaprim0[1]*gammaprim0[1] \
    + m0[5]*gammaprim0[2]*gammaprim0[2] \
    + 2.0*m0[1]*gammaprim0[0]*gammaprim0[1]  + 2.0*m0[2]*gammaprim0[0]*gammaprim0[2] \
    + 2.0*m0[4]*gammaprim0[1]*gammaprim0[2];

  l1 = m1[0]*gammaprim1[0]*gammaprim1[0] + m1[3]*gammaprim1[1]*gammaprim1[1] \
    + m1[5]*gammaprim1[2]*gammaprim1[2] \
    +2.0*m1[1]*gammaprim1[0]*gammaprim1[1]  + 2.0*m1[2]*gammaprim1[0]*gammaprim1[2] \
    + 2.0*m1[4]*gammaprim1[1]*gammaprim1[2];

  if(l0 < 0) {
    printf("%s:%d:Error: negative edge length (%e)\n",__FILE__,__LINE__,l0);
    exit(EXIT_FAILURE);
  }
  if(l1 < 0) {
    printf("%s:%d:Error: negative edge length (%e)\n",__FILE__,__LINE__,l1);
    exit(EXIT_FAILURE);
  }
  l0 = 0.5*(sqrt(l0) + sqrt(l1));

  return(l0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param np1 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of surface edge \f$[i0;i1]\f$ according to the prescribed
 * aniso metric (for special storage of metrics at ridges points). Here the
 * length is computed taking into account the curve nature of the surface edge.
 *
 */
static inline
double _MMG5_lenSurfEdg_ani(MMG5_pMesh mesh,MMG5_pSol met,int np0,int np1,char isedg) {
  MMG5_pPoint   p0,p1;
  double        *m0,*m1,met0[6],met1[6],ux,uy,uz;

  p0 = &mesh->point[np0];
  p1 = &mesh->point[np1];

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  /* Set metrics */
  if ( MG_SIN(p0->tag) || (MG_NOM & p0->tag)) {
    m0 = &met->m[6*np0];
  }
  else if ( MG_GEO & p0->tag ) {
    if ( !_MMG5_buildridmet(mesh,met,np0,ux,uy,uz,met0) )  {
      printf("%s:%d:Error: Unable to compute the metric along the ridge.\n "
             "Exit program.\n",__FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }
    m0 = met0;
  }
  else {
    m0 = &met->m[6*np0];
  }

  if ( MG_SIN(p1->tag) || (MG_NOM & p1->tag)) {
    m1 = &met->m[6*np1];
  }
  else if ( MG_GEO & p1->tag ) {
    if ( !_MMG5_buildridmet(mesh,met,np1,ux,uy,uz,met1) )  {
      printf("%s:%d:Error: Unable to compute the metric along the ridge.\n "
             "Exit program.\n",__FILE__,__LINE__);
      exit(EXIT_FAILURE);
    }
    m1 = met1;
  }
  else {
    m1 = &met->m[6*np1];
  }

  return(_MMG5_lenEdg(mesh,np0,np1,m0,m1,isedg));
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param np1 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of surface edge \f$[i0;i1]\f$ according to the prescribed
 * aniso metric (for classic storage of metrics at ridges points).
 *
 */
static inline
double _MMG5_lenSurfEdg33_ani(MMG5_pMesh mesh,MMG5_pSol met,
                              int np0,int np1,char isedg) {
  double        *m0,*m1;

  /* Set metrics */
  m0 = &met->m[6*np0];
  m1 = &met->m[6*np1];

  return(_MMG5_lenEdg(mesh,np0,np1,m0,m1,isedg));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param ip1 index of edge's extremity.
 * \param ip2 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise (dummy arg for
 * compatibility with \a lenedg_ani).
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of surface edge \f$[i0;i1]\f$ according to the prescribed iso
 * metric.
 *
 */
static
inline double _MMG5_lenSurfEdg_iso(MMG5_pMesh mesh,MMG5_pSol met,int ip1,int ip2, char isedg) {
  MMG5_pPoint   p1,p2;
  double   h1,h2,l,r,len;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  h1 = met->m[ip1];
  h2 = met->m[ip2];
  l = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]) \
    + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < _MMG5_EPS ? l / h1 : l / (h2-h1) * log(r+1.0);

  return(len);
}

#endif
