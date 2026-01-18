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
 * \file mmg3d/cenrad_3d.c
 * \brief Compute radius and center of circumscribing circle to the element.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \date 2013
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \remark Delaunay mode only (\a MMG_PATTERN flag set to \a OFF).
 */

#include "libmmg3d_private.h"
/**
 * \param mesh pointer to the mesh structure.
 * \param ct coordinates of vertices of the element.
 * \param c center of circumscribing circle to the element.
 * \param rad squared radius of circumscribing circle to the element.
 * \return 0 if failed, 1 otherwise.
 *
 * Compute radius (squared) and center of circumscribing circle to the element.
 *
 */
int MMG5_cenrad_iso(MMG5_pMesh mesh,double *ct,double *c,double *rad) {
  double      dd,ux,uy,uz,n1[3],n2[3],n3[3],*c1,*c2,*c3,*c4,pl1,pl2,pl3;
  double      cc1,cc2,cc3;

  c1 = &ct[0];
  c2 = &ct[3];
  c3 = &ct[6];
  c4 = &ct[9];

  ux = c4[0] - c1[0];
  uy = c4[1] - c1[1];
  uz = c4[2] - c1[2];

  dd = ux*ux + uy*uy + uz*uz;
  if ( dd < MMG5_EPSD2 )  return 0;

  dd = 1.0 / sqrt(dd);

  n1[0] = ux*dd;
  n1[1] = uy*dd;
  n1[2] = uz*dd;

  /* plan: vecteur directeur passant par milieu(1,4) */
  pl1 = n1[0]*(c4[0]+c1[0]) \
    + n1[1]*(c4[1]+c1[1]) + n1[2]*(c4[2]+c1[2]);

  ux = c4[0] - c2[0];
  uy = c4[1] - c2[1];
  uz = c4[2] - c2[2];

  dd = ux*ux + uy*uy + uz*uz;
  if ( dd < MMG5_EPSD2 )  return 0;

  dd = 1.0 / sqrt(dd);

  n2[0] = ux*dd;
  n2[1] = uy*dd;
  n2[2] = uz*dd;
  pl2 = n2[0]*(c4[0]+c2[0]) \
    + n2[1]*(c4[1]+c2[1]) + n2[2]*(c4[2]+c2[2]);

  ux = c4[0] - c3[0];
  uy = c4[1] - c3[1];
  uz = c4[2] - c3[2];

  dd = ux*ux + uy*uy + uz*uz;
  if ( dd < MMG5_EPSD2 )  return 0;

  dd = 1.0 / sqrt(dd);

  n3[0] = ux*dd;
  n3[1] = uy*dd;
  n3[2] = uz*dd;
  pl3 = n3[0]*(c4[0]+c3[0]) \
    + n3[1]*(c4[1]+c3[1]) + n3[2]*(c4[2]+c3[2]);

  /* center = intersection of 3 planes */
  ux = n2[1]*n3[2] - n2[2]*n3[1];
  uy = n1[2]*n3[1] - n1[1]*n3[2];
  uz = n1[1]*n2[2] - n1[2]*n2[1];

  dd = n1[0]*ux + n2[0]*uy + n3[0]*uz;
  if(fabs((dd))<1e-12)  return 0;
  dd = 0.5 / dd;

  cc1 = ux*pl1 + uy*pl2 + uz*pl3;
  cc2 = pl1 * (n2[2]*n3[0] - n2[0]*n3[2]) \
    + pl2 * (n1[0]*n3[2] - n3[0]*n1[2]) \
    + pl3 * (n2[0]*n1[2] - n2[2]*n1[0]);
  cc3 = pl1 * (n2[0]*n3[1] - n2[1]*n3[0]) \
    + pl2 * (n3[0]*n1[1] - n3[1]*n1[0]) \
    + pl3 * (n1[0]*n2[1] - n2[0]*n1[1]);

  c[0] = dd * cc1;
  c[1] = dd * cc2;
  c[2] = dd * cc3;

  /* radius (squared) */
  *rad = (c[0] - c4[0]) * (c[0] - c4[0]) \
    + (c[1] - c4[1]) * (c[1] - c4[1]) \
    + (c[2] - c4[2]) * (c[2] - c4[2]);

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param ct coordinates of vertices of the element.
 * \param m metric at the point for which we compute the cavity.
 * \param c center of circumscribing circle to the element.
 * \param rad squared radius of circumscribing circle to the element.
 * \return 0 if failed, 1 otherwise.
 *
 * Compute radius (squared) and center of circumscribing circle to the element
 * for an anisotropic metric \a m.
 *
 */
int MMG5_cenrad_ani(MMG5_pMesh mesh,double *ct,double *m,double *c,double *rad) {
  double      d1,d2,d3,det,dd,ux,uy,uz,vx,vy,vz,wx,wy,wz;
  double      ax,ay,az,bx,by,bz,cx,cy,cz;


  dd =      m[0]*ct[0]*ct[0] + m[3]*ct[1]*ct[1] + m[5]*ct[2]*ct[2] \
    + 2.0*(m[1]*ct[0]*ct[1] + m[2]*ct[0]*ct[2] + m[4]*ct[1]*ct[2]);

  /* MMG_lengths */
  d1 =      m[0]*ct[3]*ct[3] + m[3]*ct[4]*ct[4] + m[5]*ct[5]*ct[5] \
    + 2.0*(m[1]*ct[3]*ct[4] + m[2]*ct[3]*ct[5] + m[4]*ct[4]*ct[5]) - dd;

  d2 =      m[0]*ct[6]*ct[6] + m[3]*ct[7]*ct[7] + m[5]*ct[8]*ct[8] \
    + 2.0*(m[1]*ct[6]*ct[7] + m[2]*ct[6]*ct[8] + m[4]*ct[7]*ct[8]) - dd;

  d3 =      m[0]*ct[9]*ct[9] + m[3]*ct[10]*ct[10] + m[5]*ct[11]*ct[11] \
    + 2.0*(m[1]*ct[9]*ct[10] + m[2]*ct[9]*ct[11] + m[4]*ct[10]*ct[11]) - dd;

  ux = ct[3] - ct[0];
  uy = ct[4] - ct[1];
  uz = ct[5] - ct[2];

  vx = ct[6] - ct[0];
  vy = ct[7] - ct[1];
  vz = ct[8] - ct[2];

  wx = ct[9] - ct[0];
  wy = ct[10] - ct[1];
  wz = ct[11] - ct[2];

  /* M.u */
  ax = m[0]*ux + m[1]*uy + m[2]*uz;
  ay = m[1]*ux + m[3]*uy + m[4]*uz;
  az = m[2]*ux + m[4]*uy + m[5]*uz;

  bx = m[0]*vx + m[1]*vy + m[2]*vz;
  by = m[1]*vx + m[3]*vy + m[4]*vz;
  bz = m[2]*vx + m[4]*vy + m[5]*vz;

  cx = m[0]*wx + m[1]*wy + m[2]*wz;
  cy = m[1]*wx + m[3]*wy + m[4]*wz;
  cz = m[2]*wx + m[4]*wy + m[5]*wz;

  /* center */
  c[0] = d1 *(by*cz - bz*cy) - d2 * (ay*cz - az*cy) + d3 * (ay*bz - az*by); 
  c[1] = d1 *(bz*cx - bx*cz) - d2 * (az*cx - ax*cz) + d3 * (az*bx - ax*bz);
  c[2] = d1 *(bx*cy - by*cx) - d2 * (ax*cy - ay*cx) + d3 * (ax*by - ay*bx);

  det = ax * (by*cz - bz*cy) - ay * (bx*cz - bz*cx) + az * (bx*cy - cx*by);
  det = 1.0 / (2.0*det);

  c[0] *= det;
  c[1] *= det;
  c[2] *= det;

  /* radius (squared) */
  ux = ct[0] - c[0];
  uy = ct[1] - c[1];
  uz = ct[2] - c[2];
  *rad =      m[0]*ux*ux + m[3]*uy*uy + m[5]*uz*uz \
    + 2.0*(m[1]*ux*uy + m[2]*ux*uz + m[4]*uy*uz);

  return 1;
}
