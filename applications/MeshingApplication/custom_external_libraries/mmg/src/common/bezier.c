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
 * \file common/bezier.c
 * \brief Functions for Bezier surface computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgcommon_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param i0 index of the first extremity of the edge.
 * \param i1 index of the second extremity of the edge.
 * \param b0 first computed bezier coefficient.
 * \param b1 second computer bezier coefficient.
 * \param isrid is \f$[p0;p1]\f$ a special edge?
 * \param v normal to the triangle from which we come.
 *
 * Computes the Bezier coefficients associated to the underlying curve to
 * \f$[p0;p1]\f$.
 *
 */
inline void MMG5_bezierEdge(MMG5_pMesh mesh,MMG5_int i0,MMG5_int i1,
                             double b0[3],double b1[3], int8_t isrid,double v[3])
{
  MMG5_pPoint    p0,p1;
  double         ux,uy,uz,*n1,*n2,*t,ps1,ps2;

  p0 = &mesh->point[i0];
  p1 = &mesh->point[i1];

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  if ( isrid ) {
    if ( MG_SIN(p0->tag) || (p0->tag & MG_NOM) ) {
      b0[0] = p0->c[0] + MMG5_ATHIRD*ux;
      b0[1] = p0->c[1] + MMG5_ATHIRD*uy;
      b0[2] = p0->c[2] + MMG5_ATHIRD*uz;
    }
    else {
      t = &p0->n[0];
      ps1 = t[0]*ux + t[1]*uy + t[2]*uz;
      b0[0] = p0->c[0] + MMG5_ATHIRD*ps1*t[0];
      b0[1] = p0->c[1] + MMG5_ATHIRD*ps1*t[1];
      b0[2] = p0->c[2] + MMG5_ATHIRD*ps1*t[2];
    }

    if (MG_SIN(p1->tag) || (p1->tag & MG_NOM) ) {
      b1[0] = p1->c[0] - MMG5_ATHIRD*ux;
      b1[1] = p1->c[1] - MMG5_ATHIRD*uy;
      b1[2] = p1->c[2] - MMG5_ATHIRD*uz;
    }
    else {
      t = &p1->n[0];
      ps1 = -(t[0]*ux + t[1]*uy + t[2]*uz);
      b1[0] = p1->c[0] + MMG5_ATHIRD*ps1*t[0];
      b1[1] = p1->c[1] + MMG5_ATHIRD*ps1*t[1];
      b1[2] = p1->c[2] + MMG5_ATHIRD*ps1*t[2];
    }
  }

  /* regular edge */
  else {
    if ( MG_SIN(p0->tag) || (p0->tag & MG_NOM) ) {
      b0[0] = p0->c[0] + MMG5_ATHIRD*ux;
      b0[1] = p0->c[1] + MMG5_ATHIRD*uy;
      b0[2] = p0->c[2] + MMG5_ATHIRD*uz;
    }
    else {
      if ( MG_GEO & p0->tag ) {
        n1 = &mesh->xpoint[p0->xp].n1[0];
        n2 = &mesh->xpoint[p0->xp].n2[0];
        ps1 = v[0]*n1[0] + v[1]*n1[1] + v[2]*n1[2];
        ps2 = v[0]*n2[0] + v[1]*n2[1] + v[2]*n2[2];
        if ( ps1 < ps2 ) {
          n1 = &mesh->xpoint[p0->xp].n2[0];
          ps1 = ps2;
        }
      }
      else if ( (MG_REF & p0->tag) || (MG_BDY & p0->tag) ) {
        // if MG_BDY, we are in mmg3d: the normal is stored in the xPoint
        n1 = &mesh->xpoint[p0->xp].n1[0];
        ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
      }
      else {
        n1 = &p0->n[0];
        ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
      }
      b0[0] = MMG5_ATHIRD*(2.0*p0->c[0] + p1->c[0] - ps1*n1[0]);
      b0[1] = MMG5_ATHIRD*(2.0*p0->c[1] + p1->c[1] - ps1*n1[1]);
      b0[2] = MMG5_ATHIRD*(2.0*p0->c[2] + p1->c[2] - ps1*n1[2]);
    }

    if ( MG_SIN(p1->tag) || (p1->tag & MG_NOM) ) {
      b1[0] = p1->c[0] - MMG5_ATHIRD*ux;
      b1[1] = p1->c[1] - MMG5_ATHIRD*uy;
      b1[2] = p1->c[2] - MMG5_ATHIRD*uz;
    }
    else {
      if ( MG_GEO & p1->tag ) {
        n1 = &mesh->xpoint[p1->xp].n1[0];
        n2 = &mesh->xpoint[p1->xp].n2[0];
        ps1 = -(v[0]*n1[0] + v[1]*n1[1] + v[2]*n1[2]);
        ps2 = -(v[0]*n2[0] + v[1]*n2[1] + v[2]*n2[2]);
        if ( fabs(ps2) < fabs(ps1) ) {
          n1 = &mesh->xpoint[p1->xp].n2[0];
          ps1 = ps2;
        }
      }
      else if ( (MG_REF & p1->tag ) || (MG_BDY & p1->tag) ) {
        // if MG_BDY, we are in mmg3d: the normal is stored in the xPoint
        n1 = &mesh->xpoint[p1->xp].n1[0];
        ps1 = -(ux*n1[0] + uy*n1[1] + uz*n1[2]);
      }
      else {
        n1 = &p1->n[0];
        ps1 = -(ux*n1[0] + uy*n1[1] + uz*n1[2]);
      }
      b1[0] = MMG5_ATHIRD*(2.0*p1->c[0] + p0->c[0] - ps1*n1[0]);
      b1[1] = MMG5_ATHIRD*(2.0*p1->c[1] + p0->c[1] - ps1*n1[1]);
      b1[2] = MMG5_ATHIRD*(2.0*p1->c[2] + p0->c[2] - ps1*n1[2]);
    }
  }
}
