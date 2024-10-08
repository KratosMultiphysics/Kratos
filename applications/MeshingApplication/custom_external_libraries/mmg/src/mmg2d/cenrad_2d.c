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
 * \file mmg2d/cenrad_2d.c
 * \brief Compute radius and center of circumscribing circle to the element.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \date 2015
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg2d_private.h"
/**
 * \param mesh Pointer toward the mesh structure.
 * \param *ct coordinates of vertices of the element.
 * \param *c center of circumscribing circle to the element.
 * \param *rad radius of circumscribing circle to the element.
 * \return 0 if failed, 1 otherwise.
 *
 * Compute radius and center of circumscribing circle to the element.
 *
 */
int MMG2D_cenrad_iso(MMG5_pMesh mesh,double *ct,double *c,double *rad) {
  double      dd,ux,uy,n1[2],n2[2],*c1,*c2,*c3,pl1,pl2;
  double      cc1,cc2;

  c1 = &ct[0];
  c2 = &ct[2];
  c3 = &ct[4];

  ux = c3[0] - c1[0];
  uy = c3[1] - c1[1];

  dd = 1.0 / sqrt(ux*ux + uy*uy);
  n1[0] = ux*dd;
  n1[1] = uy*dd;

  /* droite passant par le milieu de c1c3 */
  pl1 = 0.5*(n1[0]*(c3[0]+c1[0])+ n1[1]*(c3[1]+c1[1])) ;

  ux = c3[0] - c2[0];
  uy = c3[1] - c2[1];

  dd = 1.0 / sqrt(ux*ux + uy*uy);
  n2[0] = ux*dd;
  n2[1] = uy*dd;
  pl2 = 0.5*(n2[0]*(c3[0]+c2[0])+ n2[1]*(c3[1]+c2[1]));

  /* center = intersection of 3 mediatrice */
  dd = n1[0]*n2[1] - n2[0]*n1[1] ;
  if(fabs((dd))<1e-12)  return 0;
  dd = 1./dd;

  cc1 = n2[1]*pl1 - n1[1]*pl2;
  cc2 = -n2[0]*pl1 + n1[0]*pl2;

  c[0] = dd * cc1;
  c[1] = dd * cc2;

  /* radius (squared) */
  *rad = (c[0] - c1[0]) * (c[0] - c1[0]) \
    + (c[1] - c1[1]) * (c[1] - c1[1]);

  /* printf("check rad %e -- %e %e\n",*rad, (c[0] - c2[0]) * (c[0] - c2[0]) \ */
  /*     + (c[1] - c2[1]) * (c[1] - c2[1]), (c[0] - c3[0]) * (c[0] - c3[0]) \ */
  /*     + (c[1] - c3[1]) * (c[1] - c3[1])); */

  return 1;
}
