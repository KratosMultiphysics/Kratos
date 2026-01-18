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
#include "libmmg2d_private.h"

// extern int8_t ddb;

/* Check if triangle k should be split based on geometric and rough edge length considerations */
int MMG2D_chkedg(MMG5_pMesh mesh, MMG5_int k) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  MMG5_pPar         par;
  double            hausd[3],hmax[3],ps,cosn,ux,uy,ll,li,t1[2],t2[2];
  int               l;
  int8_t            i,i1,i2,fill;

  pt = &mesh->tria[k];
  hausd[0] = hausd[1] = hausd[2] = mesh->info.hausd;
  hmax[0]  = hmax[1]  = hmax[2]  = mesh->info.hmax;

  /* Local parameters for pt and k */
  int8_t isloc = 0;

  if ( mesh->info.parTyp & MG_Tria ) {
    for ( l=0; l<mesh->info.npar; ++l ) {
      par = &mesh->info.par[l];

      if ( par->elt != MMG5_Triangle )  continue;
      if ( par->ref != pt->ref ) continue;

      hmax[0]  = hmax[1]  = hmax[2]  = par->hmax;
      hausd[0] = hausd[1] = hausd[2] = par->hausd;
      isloc = 1;
      break;
    }
  }

  fill = 0;
  if ( mesh->info.parTyp & MG_Edge ) {
    if ( isloc ) {
      for ( l=0; l<mesh->info.npar; ++l ) {
        par = &mesh->info.par[l];

        if ( par->elt != MMG5_Edg )  continue;

        for (i=0; i<3; i++) {
          if ( par->ref != pt->edg[i] ) continue;

          hmax[i]  = MG_MIN(hmax[i],par->hmax);
          hausd[i] = MG_MIN(hausd[i],par->hausd);

          MG_SET(fill,i);
        }

        /* Stop the loop over local parameters if all edges have been found */
        if ( fill == 7 ) {
          break;
        }

      }
    }
    else {
      for ( l=0; l<mesh->info.npar; ++l ) {
        par = &mesh->info.par[l];

        if ( par->elt != MMG5_Edg )  continue;

        for (i=0; i<3; i++) {
          if ( par->ref != pt->edg[i] ) continue;

          hmax[i]  = par->hmax;
          hausd[i] = par->hausd;

          MG_SET(fill,i);
        }

        /* Stop the loop over local parameters if all edges have been found */
        if ( fill == 7 ) {
          break;
        }
      }
    }
  }


  /* Analyze the three edges of k */
  for (i=0; i<3; i++) {
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt->v[i2]];

    /* Check the lengths of the edges */
    ux = p2->c[0] - p1->c[0];
    uy = p2->c[1] - p1->c[1];
    ll = ux*ux + uy*uy;

    /* Long edges should be split */
    if ( ll > hmax[i]*hmax[i] ) {
      MG_SET(pt->flag,i);
      continue;
    }
    /* Do not split very short edges */
    else if ( ll < MMG5_EPSD ) continue;

    /* Split non geometric edges connecting two parts of the border */
    else if ( mesh->info.fem &&
              ( (!MG_EDG(pt->tag[i])) && (p1->tag > MG_NOTAG) && (p2->tag > MG_NOTAG)) ) {
      MG_SET(pt->flag,i);
      continue;
    }

    /* Test the distance with respect to the continuous support in the case of a geometric edge */
    if ( !MG_EDG(pt->tag[i]) ) continue;

    /* Collect tangent vectors at both endpoints; remark t1 and t2 need not be oriented in the same fashion */
    if ( (MG_CRN & p1->tag) || (p1->tag & MG_NOM) ) {
      li = 1.0 / sqrt(ll);
      t1[0] = li*ux;
      t1[1] = li*uy;
    }
    else {
      t1[0] = -p1->n[1];
      t1[1] = p1->n[0];
    }

    if ( (MG_CRN & p2->tag) || (p2->tag & MG_NOM) ) {
      li = 1.0 / sqrt(ll);
      t2[0] = li*ux;
      t2[1] = li*uy;
    }
    else {
      t2[0] = -p2->n[1];
      t2[1] = p2->n[0];
    }

    /* Evaluate Hausdorff distance with respect to the geometric support */
    ps = t1[0]*ux + t1[1]*uy;
    ps *= ps;
    cosn = ps/ll ;
    cosn *= (1.0-cosn);
    cosn *= ll;
    if ( cosn > 9.0*hausd[i]*hausd[i] ) {   // Not so sure about that 9.0
      MG_SET(pt->flag,i);
      continue;
    }

    ps = -(t2[0]*ux + t2[1]*uy);
    ps *= ps;
    cosn = ps/ll ;
    cosn *= (1.0-cosn);
    cosn *= ll;
    if ( cosn > 9.0*hausd[i]*hausd[i] ) {
      MG_SET(pt->flag,i);
      continue;
    }
  }

  return pt->flag;
}


/* Calculate coordinates o[2] and interpolated normal vector no[2] of a new point
 situated at parametric distance s from i1 = inxt2[i] */
int MMG2D_bezierCurv(MMG5_pMesh mesh,MMG5_int k,int8_t i,double s,double *o,double *no) {
  MMG5_pTria         pt;
  MMG5_pPoint        p1,p2;
  double             b1[2],b2[2],t1[2],t2[2],n1[2],n2[2],bn[2],ux,uy,ll,li,ps;
  int8_t             i1,i2;

  pt = &mesh->tria[k];
  if ( !MG_EOK(pt) ) return 0;

  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];
  p1 = &mesh->point[pt->v[i1]];
  p2 = &mesh->point[pt->v[i2]];

  /* When the edge is not geometric, simply take the midpoint */
  if ( !MG_EDG(pt->tag[i]) ) {
    o[0] = (1-s)*p1->c[0]+s*p2->c[0];
    o[1] = (1-s)*p1->c[1]+s*p2->c[1];
    memset(no,0,2*sizeof(double));
    return 1;
  }

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  ll = ux*ux + uy*uy;

  if ( ll < MMG5_EPSD ) return 0;

  /* Recover normal and tangent vectors */
  if ( (MG_CRN & p1->tag) || (p1->tag & MG_NOM) ) {
    li = 1.0 / sqrt(ll);
    t1[0] = li*ux;
    t1[1] = li*uy;

    n1[0] = t1[1];
    n1[1] = - t1[0];
  }
  else {
    n1[0] = p1->n[0];
    n1[1] = p1->n[1];

    t1[0] = -p1->n[1];
    t1[1] = p1->n[0];
  }

  if ( (MG_CRN & p2->tag) || (p2->tag & MG_NOM) ) {
    li = 1.0 / sqrt(ll);
    t2[0] = li*ux;
    t2[1] = li*uy;

    n2[0] = t2[1];
    n2[1] = - t2[0];
  }
  else {
    n2[0] = p2->n[0];
    n2[1] = p2->n[1];

    t2[0] = -p2->n[1];
    t2[1] = p2->n[0];
  }

  /* When either p1 or p2 is singular, make orientation of both normal vectors consistent
   (otherwise, it is already the case) */
  if ( (MG_CRN & p1->tag) || (p1->tag & MG_NOM) ){
    ps = n1[0]*n2[0] + n1[1]*n2[1];
    if ( ps < 0.0 ) {
      n1[0] *= -1.0;
      n1[1] *= -1.0;
    }
  }
  else if ( (MG_CRN & p2->tag) || (p2->tag & MG_NOM) ) {
    ps = n1[0]*n2[0] + n1[1]*n2[1];
    if ( ps < 0.0 ) {
      n2[0] *= -1.0;
      n2[1] *= -1.0;
    }
  }

  /* Calculation of control points */
  ps = ux*t1[0]+uy*t1[1];
  b1[0] = p1->c[0] + MMG5_ATHIRD*ps*t1[0];
  b1[1] = p1->c[1] + MMG5_ATHIRD*ps*t1[1];

  ps = ux*t2[0]+uy*t2[1];
  b2[0] = p2->c[0] - MMG5_ATHIRD*ps*t2[0];
  b2[1] = p2->c[1] - MMG5_ATHIRD*ps*t2[1];

  ps = ux*(n1[0]+n2[0]) + uy*(n1[1]+n2[1]);
  ps = 2.0*ps/ll;

  bn[0] = n1[0]+n2[0] - ps*ux;
  bn[1] = n1[1]+n2[1] - ps*uy;
  ps = bn[0]*bn[0] + bn[1]*bn[1];
  if ( ps > MMG5_EPSD2 ) {
    ps = 1.0 / sqrt(ps);
    bn[0] *= ps;
    bn[1] *= ps;
  }

  /* Interpolation of the position and normal vector */
  o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p1->c[0] + 3.0*(1.0-s)*(1.0-s)*s*b1[0] + 3.0*(1.0-s)*s*s*b2[0] + s*s*s*p2->c[0];
  o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p1->c[1] + 3.0*(1.0-s)*(1.0-s)*s*b1[1] + 3.0*(1.0-s)*s*s*b2[1] + s*s*s*p2->c[1];

  no[0] = (1.0-s)*(1.0-s)*n1[0] + 2*(1.0-s)*s*bn[0] + s*s*n2[0];
  no[1] = (1.0-s)*(1.0-s)*n1[1] + 2*(1.0-s)*s*bn[1] + s*s*n2[1];

  return 1;
}
