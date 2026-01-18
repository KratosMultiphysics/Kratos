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
 * \file mmg2d/anisomovpt_2d.c
 * \brief Node relocation routines
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"

/* Relocate internal vertex whose ball is passed */
int MMG2D_movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int ilist,MMG5_int *list,int8_t improve) {
  MMG5_pTria         pt,pt0;
  MMG5_pPoint        ppt0,p0,p1,p2;
  double             calold,calnew,area,det,alpha,ps,ps1,ps2,step,sqdetm1,sqdetm2;
  double             gr[2],grp[2],*m0,*m1,*m2;
  MMG5_int           k,iel,ip0,ip1,ip2;
  int8_t             i,i1,i2;
  static int8_t      mmgWarn0=0;

  pt0 = &mesh->tria[0];
  ppt0 = &mesh->point[0];
  gr[0] = gr[1] = 0.0;
  calold = calnew = DBL_MAX;
  step = 0.1;

  /* Step 1: Calculation of the gradient of the variance function; store the quality of the previous
     configuration */
  p0 = p1 = p2 = NULL;
  m0 = NULL;
  assert ( ilist );
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    pt = &mesh->tria[iel];

    /* Quality of pt */
    calold = MG_MIN(MMG2D_caltri(mesh,met,pt),calold);

    i = list[k] % 3;
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    ip0 = pt->v[i];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];

    p0 = &mesh->point[ip0];
    p1 = &mesh->point[ip1];
    p2 = &mesh->point[ip2];

    area = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
    area = 0.5*fabs(area);

    m0 = &met->m[3*ip0];
    m1 = &met->m[3*ip1];
    m2 = &met->m[3*ip2];

    sqdetm1 = sqrt(m1[0]*m1[2]-m1[1]*m1[1]);
    sqdetm2 = sqrt(m2[0]*m2[2]-m2[1]*m2[1]);

    gr[0] += MMG5_ATHIRD*area*((p1->c[0]-p0->c[0])*sqdetm1 + (p2->c[0]-p0->c[0])*sqdetm2);
    gr[1] += MMG5_ATHIRD*area*((p1->c[1]-p0->c[1])*sqdetm1 + (p2->c[1]-p0->c[1])*sqdetm2);
  }

  /* Preconditionning of the gradient gr = M^{-1}gr */
  assert(m0 && (m0+1) && (m0+2));
  det = m0[0]*m0[2]-m0[1]*m0[1];

  if ( det < MMG5_EPSD ) return 0;
  det = 1.0 / det;

  grp[0] = det*(m0[2]*gr[0]-m0[1]*gr[1]);
  grp[1] = det*(-m0[1]*gr[0]+m0[0]*gr[1]);

  /* Step 2: Identification of the triangle such that gr is comprised in the associated angular sector */
  ps1 = ps2 = 0.;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    pt = &mesh->tria[iel];

    i = list[k] % 3;
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    ip0 = pt->v[i];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];

    p0 = &mesh->point[ip0];
    p1 = &mesh->point[ip1];
    p2 = &mesh->point[ip2];

    ps1 = (p1->c[0]-p0->c[0])*grp[1] - (p1->c[1]-p0->c[1])*grp[0];
    ps2 = grp[0]*(p2->c[1]-p0->c[1]) - grp[1]*(p2->c[0]-p0->c[0]);

    if ( ps1 >= 0.0 && ps2 >= 0.0 ) break;
  }

  if ( k == ilist ) {
    if ( !mmgWarn0 ) {
      mmgWarn0=1;
      fprintf(stderr,"\n  ## Error: %s: impossible to locate at least"
              " 1 gradient - abort.\n",__func__);
    }
    return 0;
  }

  /* coordinates of the proposed position for relocation p = p0 + alpha*set*grp, so that
     the new point is inside the triangle */
  det = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
  ps = ps1+ps2;
  if ( ps < MMG5_EPSD ) return 0;
  alpha = det / ps;

  ppt0->c[0] = p0->c[0] + alpha*step*grp[0];
  ppt0->c[1] = p0->c[1] + alpha*step*grp[1];

  /* Step 3: Vertex relocation and checks */
  i  = 0;
  pt = NULL;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i   = list[k] % 3;
    pt  = &mesh->tria[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i] = 0;

    calnew = MG_MIN(MMG2D_caltri(mesh,met,pt0),calnew);
  }

  if (calold < MMG2D_NULKAL && calnew <= calold) return 0;
  else if (calnew < MMG2D_NULKAL) return 0;
  else if ( improve && calnew < 1.02 * calold ) return 0;
  else if ( calnew < 0.3 * calold ) return 0;

  /* Update of the coordinates of the point */
  p0 = &mesh->point[pt->v[i]];
  p0->c[0] = ppt0->c[0];
  p0->c[1] = ppt0->c[1];

  return 1;
}
