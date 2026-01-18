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
 * \file mmgs/bezier_s.c
 * \brief Functions for Bezier surface computation.
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
 * \param pt pointer to the triangle structure.
 * \param pb pointer to the computed Bezier structure.
 * \param ori triangle orientation (unused but here for compatibility
 * with the MMG5_bezierCP interface).
 * \return 1.
 *
 * Compute Bezier control points on triangle \a pt (cf. \cite vlachos2001curved).
 *
 * \todo merge with the MMG5_mmg3dBezierCP function and remove the pointer
 * toward this functions.
 *
 */
int MMG5_mmgsBezierCP(MMG5_pMesh mesh,MMG5_Tria *pt,MMG5_pBezier pb,
                      int8_t ori) {
  MMG5_pPoint    p[3];
  double         *n1,*n2,nt[3],ps,ps2,dd,ux,uy,uz,ll;
  MMG5_int       ia,ib,ic;
  int8_t         i,i1,i2;

  ia   = pt->v[0];
  ib   = pt->v[1];
  ic   = pt->v[2];
  p[0] = &mesh->point[ia];
  p[1] = &mesh->point[ib];
  p[2] = &mesh->point[ic];

  memset(pb,0,sizeof(MMG5_Bezier));

  /* first 3 CP = vertices, normals */
  for (i=0; i<3; i++) {
    memcpy(&pb->b[i],p[i]->c,3*sizeof(double));
    pb->p[i] = p[i];

    if ( MS_SIN(p[i]->tag) ) {
      MMG5_nortri(mesh,pt,pb->n[i]);
    }
    else if ( MG_EDG(p[i]->tag) ) {
      MMG5_nortri(mesh,pt,nt);

      n1 = &mesh->xpoint[p[i]->xp].n1[0];
      n2 = &mesh->xpoint[p[i]->xp].n2[0];

      ps  = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
      ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];
      if ( fabs(ps) > fabs(ps2) )
        memcpy(&pb->n[i],n1,3*sizeof(double));
      else
        memcpy(&pb->n[i],n2,3*sizeof(double));
      memcpy(&pb->t[i],p[i]->n,3*sizeof(double));
    }
    else
      memcpy(&pb->n[i],p[i]->n,3*sizeof(double));
  }

  /* compute control points along edges */
  for (i=0; i<3; i++) {
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    ux = p[i2]->c[0] - p[i1]->c[0];
    uy = p[i2]->c[1] - p[i1]->c[1];
    uz = p[i2]->c[2] - p[i1]->c[2];

    ll = ux*ux + uy*uy + uz*uz;   // A PROTEGER !!!!

    /* choose normals */
    n1 = pb->n[i1];
    n2 = pb->n[i2];

    /* check for boundary curve */
    if ( MG_EDG(pt->tag[i]) ) {
      if ( MS_SIN(p[i1]->tag) ) {
        dd = 1.0 / 3.0;
        pb->b[2*i+3][0] = p[i1]->c[0] + dd * ux;
        pb->b[2*i+3][1] = p[i1]->c[1] + dd * uy;
        pb->b[2*i+3][2] = p[i1]->c[2] + dd * uz;
      }
      else {
        dd = (ux*pb->t[i1][0] + uy*pb->t[i1][1] + uz*pb->t[i1][2]) / 3.0;
        pb->b[2*i+3][0] = p[i1]->c[0] + dd * pb->t[i1][0];
        pb->b[2*i+3][1] = p[i1]->c[1] + dd * pb->t[i1][1];
        pb->b[2*i+3][2] = p[i1]->c[2] + dd * pb->t[i1][2];
      }
      if ( MS_SIN(p[i2]->tag) ) {
        dd = 1.0 / 3.0;
        pb->b[2*i+4][0] = p[i2]->c[0] - dd * ux;
        pb->b[2*i+4][1] = p[i2]->c[1] - dd * uy;
        pb->b[2*i+4][2] = p[i2]->c[2] - dd * uz;
      }
      else {
        dd = -(ux*pb->t[i2][0] + uy*pb->t[i2][1] + uz*pb->t[i2][2]) / 3.0;
        pb->b[2*i+4][0] = p[i2]->c[0] + dd * pb->t[i2][0];
        pb->b[2*i+4][1] = p[i2]->c[1] + dd * pb->t[i2][1];
        pb->b[2*i+4][2] = p[i2]->c[2] + dd * pb->t[i2][2];
      }

      /* tangent evaluation */
      ps = ux*(pb->t[i1][0]+pb->t[i2][0]) + uy*(pb->t[i1][1]+pb->t[i2][1]) + uz*(pb->t[i1][2]+pb->t[i2][2]);
      ps = 2.0 * ps / ll;
      pb->t[i+3][0] = pb->t[i1][0] + pb->t[i2][0] - ps*ux;
      pb->t[i+3][1] = pb->t[i1][1] + pb->t[i2][1] - ps*uy;
      pb->t[i+3][2] = pb->t[i1][2] + pb->t[i2][2] - ps*uz;
      dd = pb->t[i+3][0]*pb->t[i+3][0] + pb->t[i+3][1]*pb->t[i+3][1] + pb->t[i+3][2]*pb->t[i+3][2];
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        pb->t[i+3][0] *= dd;
        pb->t[i+3][1] *= dd;
        pb->t[i+3][2] *= dd;
      }
    }

    else { /* internal edge */
      ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
      pb->b[2*i+3][0] = (2.0*p[i1]->c[0] + p[i2]->c[0] - ps*n1[0]) / 3.0;
      pb->b[2*i+3][1] = (2.0*p[i1]->c[1] + p[i2]->c[1] - ps*n1[1]) / 3.0;
      pb->b[2*i+3][2] = (2.0*p[i1]->c[2] + p[i2]->c[2] - ps*n1[2]) / 3.0;

      ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
      pb->b[2*i+4][0] = (2.0*p[i2]->c[0] + p[i1]->c[0] - ps*n2[0]) / 3.0;
      pb->b[2*i+4][1] = (2.0*p[i2]->c[1] + p[i1]->c[1] - ps*n2[1]) / 3.0;
      pb->b[2*i+4][2] = (2.0*p[i2]->c[2] + p[i1]->c[2] - ps*n2[2]) / 3.0;
    }

    /* normal evaluation */
    ps = ux*(n1[0]+n2[0]) + uy*(n1[1]+n2[1]) + uz*(n1[2]+n2[2]);
    ps = 2.0*ps / ll;
    pb->n[i+3][0] = n1[0] + n2[0] - ps*ux;
    pb->n[i+3][1] = n1[1] + n2[1] - ps*uy;
    pb->n[i+3][2] = n1[2] + n2[2] - ps*uz;
    dd = pb->n[i+3][0]*pb->n[i+3][0] + pb->n[i+3][1]*pb->n[i+3][1] + pb->n[i+3][2]*pb->n[i+3][2];
    if ( dd > MMG5_EPSD2 ) {
      dd = 1.0 / sqrt(dd);
      pb->n[i+3][0] *= dd;
      pb->n[i+3][1] *= dd;
      pb->n[i+3][2] *= dd;
    }
  }

  /* Central Bezier coefficient */
  for (i=0; i<3; i++) {
    dd = 0.5 / 3.0;
    pb->b[9][0] -= dd * pb->b[i][0];
    pb->b[9][1] -= dd * pb->b[i][1];
    pb->b[9][2] -= dd * pb->b[i][2];
  }
  for (i=0; i<3; i++) {
    pb->b[9][0] += 0.25 * (pb->b[2*i+3][0] + pb->b[2*i+4][0]);
    pb->b[9][1] += 0.25 * (pb->b[2*i+3][1] + pb->b[2*i+4][1]);
    pb->b[9][2] += 0.25 * (pb->b[2*i+3][2] + pb->b[2*i+4][2]);
  }

  return 1;
}

/**
 * \param pb pointer to the Bezier structure.
 * \param uv coordinates of the point in the parametric space.
 * \param o computed coordinates of the point in the real space.
 * \param no computed normal.
 * \param to computed tangent.
 * \return 1.
 *
 * Compute \a o, \a no and \a to at \f$(u,v)\f$ in Bezier patch.
 *
 */
int MMGS_bezierInt(MMG5_pBezier pb,double uv[2],double o[3],double no[3],double to[3]) {
  double    dd,u,v,w,ps,ux,uy,uz;
  int8_t    i;

  memset(to,0,3*sizeof(double));
  u = uv[0];
  v = uv[1];
  w = 1 - u - v;

  /* coordinates + normals */
  for (i=0; i<3; i++) {
    o[i]  = pb->b[0][i]*w*w*w + pb->b[1][i]*u*u*u + pb->b[2][i]*v*v*v \
      + 3.0 * (pb->b[3][i]*u*u*v + pb->b[4][i]*u*v*v + pb->b[5][i]*w*v*v \
               + pb->b[6][i]*w*w*v + pb->b[7][i]*w*w*u + pb->b[8][i]*w*u*u)\
      + 6.0 * pb->b[9][i]*u*v*w;

    /* quadratic interpolation of normals */
    no[i] =        pb->n[0][i]*w*w + pb->n[1][i]*u*u + pb->n[2][i]*v*v \
      + 2.0 * (pb->n[3][i]*u*v + pb->n[4][i]*v*w + pb->n[5][i]*u*w);

    /* linear interpolation, not used here
       no[i] = pb->n[0][i]*w + pb->n[1][i]*u + pb->n[2][i]*v; */
  }

  /* tangent */
  if ( w < MMG5_EPSD2 ) {
    ux = pb->b[2][0] - pb->b[1][0];
    uy = pb->b[2][1] - pb->b[1][1];
    uz = pb->b[2][2] - pb->b[1][2];
    dd = ux*ux + uy*uy + uz*uz;
    if ( dd > MMG5_EPSD2 ) {
      dd = 1.0 / sqrt(dd);
      ux *= dd;
      uy *= dd;
      uz *= dd;
    }

    /* corners */
    if ( MG_SIN(pb->p[1]->tag) ) {
      pb->t[1][0] = ux;
      pb->t[1][1] = uy;
      pb->t[1][2] = uz;
    }
    if ( MG_SIN(pb->p[2]->tag) ) {
      pb->t[2][0] = ux;
      pb->t[2][1] = uy;
      pb->t[2][2] = uz;
    }

    ps = pb->t[1][0]* pb->t[2][0] + pb->t[1][1]* pb->t[2][1] + pb->t[1][2]* pb->t[2][2];
    if ( ps > 0.0 ) {
      to[0] = pb->t[1][0]*u + pb->t[2][0]*v;
      to[1] = pb->t[1][1]*u + pb->t[2][1]*v;
      to[2] = pb->t[1][2]*u + pb->t[2][2]*v;
    }
    else {
      to[0] = -pb->t[1][0]*u + pb->t[2][0]*v;
      to[1] = -pb->t[1][1]*u + pb->t[2][1]*v;
      to[2] = -pb->t[1][2]*u + pb->t[2][2]*v;
    }
  }

  if ( u < MMG5_EPSD2 ) {
    ux = pb->b[2][0] - pb->b[0][0];
    uy = pb->b[2][1] - pb->b[0][1];
    uz = pb->b[2][2] - pb->b[0][2];
    dd = ux*ux + uy*uy + uz*uz;
    if ( dd > MMG5_EPSD2 ) {
      dd = 1.0 / sqrt(dd);
      ux *= dd;
      uy *= dd;
      uz *= dd;
    }

    /* corners */
    if ( MG_SIN(pb->p[0]->tag) ) {
      pb->t[0][0] = ux;
      pb->t[0][1] = uy;
      pb->t[0][2] = uz;
    }
    if ( MG_SIN(pb->p[2]->tag) ) {
      pb->t[2][0] = ux;
      pb->t[2][1] = uy;
      pb->t[2][2] = uz;
    }

    ps = pb->t[0][0]* pb->t[2][0] + pb->t[0][1]* pb->t[2][1] + pb->t[0][2]* pb->t[2][2];
    if ( ps > 0.0 ) {
      to[0] = pb->t[0][0]*w + pb->t[2][0]*v;
      to[1] = pb->t[0][1]*w + pb->t[2][1]*v;
      to[2] = pb->t[0][2]*w + pb->t[2][2]*v;
    }
    else {
      to[0] = -pb->t[0][0]*w + pb->t[2][0]*v;
      to[1] = -pb->t[0][1]*w + pb->t[2][1]*v;
      to[2] = -pb->t[0][2]*w + pb->t[2][2]*v;
    }
  }

  if ( v < MMG5_EPSD2 ) {
    ux = pb->b[1][0] - pb->b[0][0];
    uy = pb->b[1][1] - pb->b[0][1];
    uz = pb->b[1][2] - pb->b[0][2];
    dd = ux*ux + uy*uy + uz*uz;
    if ( dd > MMG5_EPSD2 ) {
      dd = 1.0 / sqrt(dd);
      ux *= dd;
      uy *= dd;
      uz *= dd;
    }

    /* corners */
    if ( MG_SIN(pb->p[0]->tag) ) {
      pb->t[0][0] = ux;
      pb->t[0][1] = uy;
      pb->t[0][2] = uz;
    }
    if ( MG_SIN(pb->p[1]->tag) ) {
      pb->t[1][0] = ux;
      pb->t[1][1] = uy;
      pb->t[1][2] = uz;
    }

    ps = pb->t[0][0]* pb->t[1][0] + pb->t[0][1]* pb->t[1][1] + pb->t[0][2]* pb->t[1][2];
    if ( ps > 0.0 ) {
      to[0] = pb->t[0][0]*w + pb->t[1][0]*u;
      to[1] = pb->t[0][1]*w + pb->t[1][1]*u;
      to[2] = pb->t[0][2]*w + pb->t[1][2]*u;
    }
    else {
      to[0] = -pb->t[0][0]*w + pb->t[1][0]*u;
      to[1] = -pb->t[0][1]*w + pb->t[1][1]*u;
      to[2] = -pb->t[0][2]*w + pb->t[1][2]*u;
    }
  }

  dd = no[0]*no[0] + no[1]*no[1] + no[2]*no[2];
  if ( dd > MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    no[0] *= dd;
    no[1] *= dd;
    no[2] *= dd;
  }

  dd = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
  if ( dd > MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    to[0] *= dd;
    to[1] *= dd;
    to[2] *= dd;
  }

  return 1;
}
