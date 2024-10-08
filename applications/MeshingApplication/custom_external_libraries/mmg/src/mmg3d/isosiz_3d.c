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
 * \file mmg3d/isosiz_3d.c
 * \brief Fonctions for isotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "inlined_functions_private.h"
#include "mmgexterns_private.h"

extern int8_t ddb;

#define MAXLEN    1.0e9
#define A64TH     0.015625
#define A16TH     0.0625
#define A32TH     0.03125


/**
 * \brief Compute edge length from edge's coordinates.
 * \param *ca pointer toward the coordinates of the first edge's extremity.
 * \param *cb pointer toward the coordinates of the second edge's extremity.
 * \param *ma pointer toward the metric associated to the first edge's extremity.
 * \param *mb pointer toward the metric associated to the second edge's extremity.
 * \return edge length.
 *
 * Compute length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge extremities) according to the isotropic size
 * prescription.
 *
 */
inline double MMG5_lenedgCoor_iso(double *ca,double *cb,double *ma,double *mb) {
  double   h1,h2,l,r,len;

  h1 = *ma;
  h2 = *mb;
  l = (cb[0]-ca[0])*(cb[0]-ca[0]) + (cb[1]-ca[1])*(cb[1]-ca[1]) \
    + (cb[2]-ca[2])*(cb[2]-ca[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < MMG5_EPS ? l / h1 : l / (h2-h1) * log1p(r);

  return len;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param nump index of point in which the size must be computed.
 * \param lists pointer toward the surfacic ball of \a nump.
 * \param ilists size of surfacic ball of \a nump.
 * \param hmin minimal edge size.
 * \param hmax maximal edge size.
 * \param hausd hausdorff value.
 * \return the isotropic size at the point if success, FLT_MAX if fail.
 *
 * Define isotropic size at regular point nump, whose surfacic ball is provided
 * and update metric at 'regular' non-manifold points of the surfacic ball of
 * point.
 *
 */
static double
MMG5_defsizreg(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int nump,MMG5_int *lists,
                int ilists, double hmin,double hmax,double hausd) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  MMG5_Tria     tt;
  MMG5_Bezier   b;
  double        ux,uy,uz,det2d,h,isqhmin,isqhmax,ll,lmin,lmax,hnm,s;
  double        *n,*t,r[3][3],lispoi[3*MMG3D_LMAX+1],intm[3],b0[3],b1[3],c[3],tAA[6],tAb[3],d[3];
  double        kappa[2],vp[2][2];
  MMG5_int      k,na,nb,ntempa,ntempb,iel,ip0;
  int8_t        iface,i,j,i0;
  static int8_t mmgWarn0=0,mmgWarn1=0,mmgWarn2=0,mmgWarn3=0;

  p0 = &mesh->point[nump];

  if ( (!p0->xp) || MG_EDG_OR_NOM(p0->tag) || (p0->tag & MG_REQ))  {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Error: %s: at least 1 wrong point"
              " qualification : xp ? %" MMG5_PRId ".\n",__func__,p0->xp);
    }
    return FLT_MAX;
  }

  /* Check that we don't pass here for a corner (not directly tested previously
   * because corners doesn't have xpoints) */
  assert ( (!(MG_CRN & p0->tag)) && "We should not pass here with corner points" );

  isqhmin = 1.0 / (hmin*hmin);
  isqhmax = 1.0 / (hmax*hmax);

  n = &mesh->xpoint[p0->xp].n1[0];

  /* Step 1 : rotation matrix that sends normal n to the third coordinate vector
   * of R^3 */
  if ( !MMG5_rotmatrix(n,r) ) {
    if ( !mmgWarn1 ) {
      mmgWarn1 = 1;
      fprintf(stderr,"\n  ## Warning: %s: function MMG5_rotmatrix return 0.\n",
              __func__);
    }
    return FLT_MAX;
  }

  /* Step 2 : rotation of the oriented surfacic ball with r : lispoi[k] is the
     common edge between faces lists[k-1] and lists[k] */
  iel   = lists[0] / 4;
  iface = lists[0] % 4;
  pt    = &mesh->tetra[iel];
  lmin  = MAXLEN;
  lmax  = 0.0;

  na = nb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[MMG5_idir[iface][i]] != nump ) {
      if ( !na )
        na = pt->v[MMG5_idir[iface][i]];
      else
        nb = pt->v[MMG5_idir[iface][i]];
    }
  }

  for (k=1; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    ntempa = ntempb = 0;
    for (i=0; i<3; i++) {
      if ( pt->v[MMG5_idir[iface][i]] != nump ) {
        if ( !ntempa )
          ntempa = pt->v[MMG5_idir[iface][i]];
        else
          ntempb = pt->v[MMG5_idir[iface][i]];
      }
    }
    if ( ntempa == na )
      p1 = &mesh->point[na];
    else if ( ntempa == nb )
      p1 = &mesh->point[nb];
    else if ( ntempb == na )
      p1 = &mesh->point[na];
    else {
      assert(ntempb == nb);
      p1 = &mesh->point[nb];
    }
    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

    ll = lispoi[3*k+1]*lispoi[3*k+1] + lispoi[3*k+2]*lispoi[3*k+2] + lispoi[3*k+3]*lispoi[3*k+3];
    lmin = MG_MIN(lmin,ll);
    lmax = MG_MAX(lmax,ll);

    na = ntempa;
    nb = ntempb;
  }

  /* Finish with point 0 */
  iel   = lists[0] / 4;
  iface = lists[0] % 4;
  pt    = &mesh->tetra[iel];
  ntempa = ntempb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[MMG5_idir[iface][i]] != nump ) {
      if ( !ntempa )
        ntempa = pt->v[MMG5_idir[iface][i]];
      else
        ntempb = pt->v[MMG5_idir[iface][i]];
    }
  }
  if ( ntempa == na )
    p1 = &mesh->point[na];
  else if ( ntempa == nb )
    p1 = &mesh->point[nb];
  else if ( ntempb == na )
    p1 = &mesh->point[na];
  else {
    assert(ntempb == nb);
    p1 = &mesh->point[nb];
  }

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  lispoi[1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
  lispoi[2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
  lispoi[3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

  ll = lispoi[1]*lispoi[1] + lispoi[2]*lispoi[2] + lispoi[3]*lispoi[3];
  lmin = MG_MIN(lmin,ll);
  lmax = MG_MAX(lmax,ll);

  /* list goes modulo ilist */
  lispoi[3*ilists+1] = lispoi[1];
  lispoi[3*ilists+2] = lispoi[2];
  lispoi[3*ilists+3] = lispoi[3];

  /* At this point, lispoi contains the oriented surface ball of point p0, that has been rotated
     through r, with the convention that triangle l has edges lispoi[l]; lispoi[l+1] */
  if ( lmax/lmin > 4.0*hmax*hmax/(hmin*hmin) )  return hmax;

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 )  return hmax;
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 )    return hmax;

  /* Reconstitution of the curvature tensor at p0 in the tangent plane,
     with a quadric fitting approach */
  memset(intm,0x00,3*sizeof(double));
  memset(tAA,0x00,6*sizeof(double));
  memset(tAb,0x00,3*sizeof(double));

  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;

    assert( 0<=iface && iface<4 && "unexpected local face idx");
    MMG5_tet2tri(mesh,iel,iface,&tt);

    pxt   = &mesh->xtetra[mesh->tetra[iel].xt];
    if ( !MMG5_bezierCP(mesh,&tt,&b,MG_GET(pxt->ori,iface)) ) {
      if ( !mmgWarn2 ) {
        mmgWarn2 = 1;
        fprintf(stderr,"\n  ## Warning: %s: function MMG5_bezierCP return 0.\n",
                __func__);
      }
      return FLT_MAX;
    }

    for (i0=0; i0<3; i0++) {
      if ( tt.v[i0] == nump )  break;
    }
    assert(i0 < 3);

    for (j=0; j<10; j++) {
      c[0] = b.b[j][0] - p0->c[0];
      c[1] = b.b[j][1] - p0->c[1];
      c[2] = b.b[j][2] - p0->c[2];

      b.b[j][0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
      b.b[j][1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
      b.b[j][2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];
    }

    /* Mid-point along left edge and endpoint in the rotated frame */
    if ( i0 == 0 ) {
      memcpy(b0,&(b.b[7][0]),3*sizeof(double));
      memcpy(b1,&(b.b[8][0]),3*sizeof(double));
    }
    else if ( i0 == 1 ) {
      memcpy(b0,&(b.b[3][0]),3*sizeof(double));
      memcpy(b1,&(b.b[4][0]),3*sizeof(double));
    }
    else {
      memcpy(b0,&(b.b[5][0]),3*sizeof(double));
      memcpy(b1,&(b.b[6][0]),3*sizeof(double));
    }
    s = 0.5;

    /* At this point, the two control points are expressed in the rotated frame */
    c[0] = 3.0*s*(1.0-s)*(1.0-s)*b0[0] + 3.0*s*s*(1.0-s)*b1[0] + s*s*s*lispoi[3*k+1];
    c[1] = 3.0*s*(1.0-s)*(1.0-s)*b0[1] + 3.0*s*s*(1.0-s)*b1[1] + s*s*s*lispoi[3*k+2];
    c[2] = 3.0*s*(1.0-s)*(1.0-s)*b0[2] + 3.0*s*s*(1.0-s)*b1[2] + s*s*s*lispoi[3*k+3];

    /* Fill matric tAA and second member tAb*/
    tAA[0] += c[0]*c[0]*c[0]*c[0];
    tAA[1] += c[0]*c[0]*c[1]*c[1];
    tAA[2] += c[0]*c[0]*c[0]*c[1];
    tAA[3] += c[1]*c[1]*c[1]*c[1];
    tAA[4] += c[0]*c[1]*c[1]*c[1];
    tAA[5] += c[0]*c[0]*c[1]*c[1];

    tAb[0] += c[0]*c[0]*c[2];
    tAb[1] += c[1]*c[1]*c[2];
    tAb[2] += c[0]*c[1]*c[2];

    s = 1.0;
    /* At this point, the two control points are expressed in the rotated frame */
    c[0] = 3.0*s*(1.0-s)*(1.0-s)*b0[0] + 3.0*s*s*(1.0-s)*b1[0] + s*s*s*lispoi[3*k+1];
    c[1] = 3.0*s*(1.0-s)*(1.0-s)*b0[1] + 3.0*s*s*(1.0-s)*b1[1] + s*s*s*lispoi[3*k+2];
    c[2] = 3.0*s*(1.0-s)*(1.0-s)*b0[2] + 3.0*s*s*(1.0-s)*b1[2] + s*s*s*lispoi[3*k+3];

    /* Fill matric tAA and second member tAb*/
    tAA[0] += c[0]*c[0]*c[0]*c[0];
    tAA[1] += c[0]*c[0]*c[1]*c[1];
    tAA[2] += c[0]*c[0]*c[0]*c[1];
    tAA[3] += c[1]*c[1]*c[1]*c[1];
    tAA[4] += c[0]*c[1]*c[1]*c[1];
    tAA[5] += c[0]*c[0]*c[1]*c[1];

    tAb[0] += c[0]*c[0]*c[2];
    tAb[1] += c[1]*c[1]*c[2];
    tAb[2] += c[0]*c[1]*c[2];

    /* Mid-point along median edge and endpoint in the rotated frame */
    if ( i0 == 0 ) {
      c[0] = A64TH*(b.b[1][0] + b.b[2][0] + 3.0*(b.b[3][0] + b.b[4][0])) \
        + 3.0*A16TH*(b.b[6][0] + b.b[7][0] + b.b[9][0]) + A32TH*(b.b[5][0] + b.b[8][0]);
      c[1] = A64TH*(b.b[1][1] + b.b[2][1] + 3.0*(b.b[3][1] + b.b[4][1])) \
        + 3.0*A16TH*(b.b[6][1] + b.b[7][1] + b.b[9][1]) + A32TH*(b.b[5][1] + b.b[8][1]);
      c[2] = A64TH*(b.b[1][2] + b.b[2][2] + 3.0*(b.b[3][2] + b.b[4][2])) \
        + 3.0*A16TH*(b.b[6][2] + b.b[7][2] + b.b[9][2]) + A32TH*(b.b[5][2] + b.b[8][2]);

      d[0] = 0.125*b.b[1][0] + 0.375*(b.b[3][0] + b.b[4][0]) + 0.125*b.b[2][0];
      d[1] = 0.125*b.b[1][1] + 0.375*(b.b[3][1] + b.b[4][1]) + 0.125*b.b[2][1];
      d[2] = 0.125*b.b[1][2] + 0.375*(b.b[3][2] + b.b[4][2]) + 0.125*b.b[2][2];
    }
    else if (i0 == 1) {
      c[0] = A64TH*(b.b[0][0] + b.b[2][0] + 3.0*(b.b[5][0] + b.b[6][0])) \
        + 3.0*A16TH*(b.b[3][0] + b.b[8][0] + b.b[9][0]) + A32TH*(b.b[4][0] + b.b[7][0]);
      c[1] = A64TH*(b.b[0][1] + b.b[2][1] + 3.0*(b.b[5][1] + b.b[6][1])) \
        + 3.0*A16TH*(b.b[3][1] + b.b[8][1] + b.b[9][1]) + A32TH*(b.b[4][1] + b.b[7][1]);
      c[2] = A64TH*(b.b[0][2] + b.b[2][2] + 3.0*(b.b[5][2] + b.b[6][2])) \
        + 3.0*A16TH*(b.b[3][2] + b.b[8][2] + b.b[9][2]) + A32TH*(b.b[4][2] + b.b[7][2]);

      d[0] = 0.125*b.b[2][0] + 0.375*(b.b[5][0] + b.b[6][0]) + 0.125*b.b[0][0];
      d[1] = 0.125*b.b[2][1] + 0.375*(b.b[5][1] + b.b[6][1]) + 0.125*b.b[0][1];
      d[2] = 0.125*b.b[2][2] + 0.375*(b.b[5][2] + b.b[6][2]) + 0.125*b.b[0][2];
    }
    else {
      c[0] = A64TH*(b.b[0][0] + b.b[1][0] + 3.0*(b.b[7][0] + b.b[8][0])) \
        + 3.0*A16TH*(b.b[4][0] + b.b[5][0] + b.b[9][0]) + A32TH*(b.b[3][0] + b.b[6][0]);
      c[1] = A64TH*(b.b[0][1] + b.b[1][1] + 3.0*(b.b[7][1] + b.b[8][1])) \
        + 3.0*A16TH*(b.b[4][1] + b.b[5][1] + b.b[9][1]) + A32TH*(b.b[3][1] + b.b[6][1]);
      c[2] = A64TH*(b.b[0][2] + b.b[1][2] + 3.0*(b.b[7][2] + b.b[8][2])) \
        + 3.0*A16TH*(b.b[4][2] + b.b[5][2] + b.b[9][2]) + A32TH*(b.b[3][2] + b.b[6][2]);

      d[0] = 0.125*b.b[0][0] + 0.375*(b.b[7][0] + b.b[8][0]) + 0.125*b.b[1][0];
      d[1] = 0.125*b.b[0][1] + 0.375*(b.b[7][1] + b.b[8][1]) + 0.125*b.b[1][1];
      d[2] = 0.125*b.b[0][2] + 0.375*(b.b[7][2] + b.b[8][2]) + 0.125*b.b[1][2];
    }

    /* Fill matric tAA and second member tAb*/
    tAA[0] += c[0]*c[0]*c[0]*c[0];
    tAA[1] += c[0]*c[0]*c[1]*c[1];
    tAA[2] += c[0]*c[0]*c[0]*c[1];
    tAA[3] += c[1]*c[1]*c[1]*c[1];
    tAA[4] += c[0]*c[1]*c[1]*c[1];
    tAA[5] += c[0]*c[0]*c[1]*c[1];

    tAb[0] += c[0]*c[0]*c[2];
    tAb[1] += c[1]*c[1]*c[2];
    tAb[2] += c[0]*c[1]*c[2];

    tAA[0] += d[0]*d[0]*d[0]*d[0];
    tAA[1] += d[0]*d[0]*d[1]*d[1];
    tAA[2] += d[0]*d[0]*d[0]*d[1];
    tAA[3] += d[1]*d[1]*d[1]*d[1];
    tAA[4] += d[0]*d[1]*d[1]*d[1];
    tAA[5] += d[0]*d[0]*d[1]*d[1];

    tAb[0] += d[0]*d[0]*d[2];
    tAb[1] += d[1]*d[1]*d[2];
    tAb[2] += d[0]*d[1]*d[2];
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !MMG5_sys33sym(tAA,tAb,c) )  return hmax;

  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
     are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
  if( !MMG5_eigensym(intm,kappa,vp) ){
    if ( !mmgWarn3 ) {
      mmgWarn3 = 1;
      fprintf(stderr,"\n  # Warning: %s: function MMG5_eigensym return 0.\n",
              __func__);
    }
    return FLT_MAX;
  }

  /* h computation : h(x) = sqrt( 9*hausd / (2 * max(kappa1(x),kappa2(x)) ) */
  kappa[0] = 2.0/9.0 * fabs(kappa[0]) / hausd;
  kappa[0] = MG_MIN(kappa[0],isqhmin);
  kappa[0] = MG_MAX(kappa[0],isqhmax);

  kappa[1] = 2.0/9.0 * fabs(kappa[1]) / hausd;
  kappa[1] = MG_MIN(kappa[1],isqhmin);
  kappa[1] = MG_MAX(kappa[1],isqhmax);

  kappa[0] = 1.0 / sqrt(kappa[0]);
  kappa[1] = 1.0 / sqrt(kappa[1]);

  h = MG_MIN(kappa[0],kappa[1]);

  /* Travel surfacic ball one last time and update metric of non manifold
   * regular points of ball */
  for (k=0; k<ilists; k++) {
    iface = lists[k] % 4;

    for (j=0; j<3; j++) {
      i0  = MMG5_idir[iface][j];
      ip0 = pt->v[i0];
      p1  = &mesh->point[ip0];

      if ( MG_SIN(p1->tag) ) {
        /* Do not treat singular points */
        continue;
      }

      if ( !(p1->tag & MG_NOM) ) {
        /* Do not treat manifold points */
        continue;
      }

      assert(p1->xp);
      t = &p1->n[0];
      memcpy(c,t,3*sizeof(double));

      d[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
      d[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];

      hnm = intm[0]*d[0]*d[0] + 2.0*intm[1]*d[0]*d[1] + intm[2]*d[1]*d[1];
      hnm = 2.0/9.0 * fabs(hnm) / hausd;
      hnm = MG_MIN(hnm,isqhmin);
      hnm = MG_MAX(hnm,isqhmax);
      hnm = 1.0 / sqrt(hnm);
      met->m[ip0] = MG_MIN(met->m[ip0],hnm);
    }
  }
  return h;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param nump index of point in which the size must be computed.
 * \param lists pointer toward the surfacic ball of \a nump.
 * \param ilists size of surfacic ball of \a nump.
 * \param hmin minimal edge size.
 * \param hmax maximal edge size.
 * \return the isotropic size at the point.
 *
 * For -nosurf option : define isotropic size at regular point nump, whose
 * surfacic ball is provided.  The size is computed as the mean of the length of
 * the surface edges passing through \a nump.
 *
 */
double MMG5_meansizreg_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int nump,MMG5_int *lists,
                int ilists, double hmin,double hmax) {
  MMG5_pTetra       pt;
  MMG5_pPoint       p0,p1;
  double            len,ux,uy,uz;
  MMG5_int          k,iel,ip1;
  int8_t            i,iface;

  p0 = &mesh->point[nump];

  len = 0;
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;
    pt    = &mesh->tetra[iel];

    for (i=0; i<3; i++) {
      if ( pt->v[MMG5_idir[iface][i]] == nump ) {
        break;
      }
    }
    assert(i!=3);

    ip1 = pt->v[MMG5_idir[iface][MMG5_inxt2[i]]];
    p1 = &mesh->point[ip1];

    ux  = p1->c[0] - p0->c[0];
    uy  = p1->c[1] - p0->c[1];
    uz  = p1->c[2] - p0->c[2];
    len += sqrt(ux*ux + uy*uy + uz*uz);
  }
  len /=ilists;

  return MG_MIN(hmax,MG_MAX(hmin,len));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param hash edge hashtable.
 * \param pt tetra to process.
 * \param i index of the edge of the tetra \a pt that we process.
 *
 * \return 1 if success, 0 if fail.
 *
 * If the edge \a i of the tetra \a pt is seen for the first time, compute its
 * euclidean length, add this length to the metric of the edge extremities and
 * increment the count of times we have processed this extremities.
 *
 */
static inline
int MMG3D_sum_reqEdgeLengthsAtPoint(MMG5_pMesh mesh,MMG5_pSol met,MMG5_Hash *hash,
                                  MMG5_pTetra pt,int8_t i) {
  MMG5_int         ip0,ip1;

  ip0 = pt->v[MMG5_iare[i][0]];
  ip1 = pt->v[MMG5_iare[i][1]];

  /* Check if the edge is already treated */
  if ( MMG5_hashGet(hash,ip0,ip1) ) return 1;

  /* Mark the edge as treated */
  if ( !MMG5_hashEdge(mesh,hash,ip0,ip1,1) ) return 0;

  if ( !MMG5_sum_reqEdgeLengthsAtPoint(mesh,met,ip0,ip1) )
    return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param ismet 1 if user provided metric
 *
 * \return 0 if fail, 1 otherwise
 *
 * Compute the metric at points on trequired adges as the mean of the lengths of
 * the required eges to which belongs the point. The processeed points are
 * marked with flag 3.
 *
 */
int MMG3D_set_metricAtPointsOnReqEdges ( MMG5_pMesh mesh,MMG5_pSol met,int8_t ismet) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_Hash    hash;
  MMG5_int     k,j,ip0,ip1,iad0,iad1;
  int          i;

  /* Reset the input metric at required edges extremities */
  if ( ismet ) {
    for ( k=1; k<=mesh->ne; k++ ) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      if ( pt->tag & MG_REQ ) {
        for ( i=0; i<6; i++ ) {
          ip0 = pt->v[MMG5_iare[i][0]];
          ip1 = pt->v[MMG5_iare[i][1]];
          iad0 = met->size*ip0;
          iad1 = met->size*ip1;
          for ( j=0; j<met->size; ++j ) {
            met->m[iad0+j] = 0.;
            met->m[iad1+j] = 0.;
          }
        }
      }
      else {
        if ( !pt->xt ) continue;
        pxt = &mesh->xtetra[pt->xt];

        for ( i=0; i<6; i++ ) {
          if ( (pxt->tag[i] & MG_REQ) || (pxt->tag[i] & MG_NOSURF) ||
               (pxt->tag[i] & MG_PARBDY) ) {
            ip0 = pt->v[MMG5_iare[i][0]];
            ip1 = pt->v[MMG5_iare[i][1]];
            iad0 = met->size*ip0;
            iad1 = met->size*ip1;
            for ( j=0; j<met->size; ++j ) {
              met->m[iad0+j] = 0.;
              met->m[iad1+j] = 0.;
            }
          }
        }
      }
    }
  }

  /* Process the required edges and add the edge length to the metric of the
   * edge extremities */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return 0;

  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    if ( pt->tag & MG_REQ ) {
      for ( i=0; i<6; i++ ) {
        if ( !MMG3D_sum_reqEdgeLengthsAtPoint(mesh,met,&hash,pt,i) ) {
          MMG5_DEL_MEM(mesh,hash.item);
          return 0;
        }
      }
    }
    else {
      if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];

      for ( i=0; i<6; i++ ) {
        if ( (pxt->tag[i] & MG_REQ) || (pxt->tag[i] & MG_NOSURF) ||
             (pxt->tag[i] & MG_PARBDY) ) {
          if ( !MMG3D_sum_reqEdgeLengthsAtPoint(mesh,met,&hash,pt,i) ) {
            MMG5_DEL_MEM(mesh,hash.item);
            return 0;
          }
        }
      }
    }
  }
  MMG5_DEL_MEM(mesh,hash.item);

  /* Travel the points and compute the metric of the points belonging to
   * required edges as the mean of the required edges length */
  if ( !MMG5_compute_meanMetricAtMarkedPoints ( mesh,met ) ) {
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Define isotropic size map at all boundary vertices of the mesh, associated
 * with geometric approx, and prescribe hmax at the internal vertices Field h of
 * Point is used, to store the prescribed size (not inverse, squared,...)
 *
 */
int MMG3D_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra    pt,ptloc;
  MMG5_pPrism    pp;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    p0,p1;
  double         hp,v[3],b0[3],b1[3],b0p0[3],b1b0[3],p1b1[3],hausd,hmin,hmax;
  double         secder0[3],secder1[3],kappa,tau[3],gammasec[3],ntau2,intau,ps,lm;
  MMG5_int       lists[MMG3D_LMAX+2],k,ip0,ip1;
  int64_t        listv[MMG3D_LMAX+2];
  int            ilists,ilistv,l;
  int            kk,isloc;
  int8_t         ismet;
  int8_t         i,j,ia,ised,i0,i1;
  MMG5_pPar      par;

  if ( !MMG5_defsiz_startingMessage (mesh,met,__func__) ) {
    return 0;
  }

  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
    p0->s    = 0;
  }

  /** 1) Size at internal points */
  hmax = DBL_MAX;
  hmin = 0.;

  /* alloc structure */
  if ( !met->m ) {
    ismet = 0;

    /* Allocate and store the header informations for each solution */
    if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,MMG5_Scalar) ) {
      return 0;
    }
  }
  else {
    ismet = 1;
    assert ( met->np );
   }

  /** Step 1: Set metric at points belonging to a required edge: compute the
   * metric as the mean of the length of the required eges passing through the
   * point */
  if ( !mesh->info.nosizreq ) {
    if ( !MMG3D_set_metricAtPointsOnReqEdges ( mesh,met,ismet ) ) {
      return 0;
    }
  }

  /** Step 2: size at non required internal points */
  if ( !ismet ) {
    /* init constant size */
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<4; i++) {
        ip0 = pt->v[i];
        p0  = &mesh->point[ip0];

        if ( p0->flag ) continue;

        /** First step: search for local parameters */
        isloc   = 0;
        hmax = mesh->info.hmax;

        /* Local param at vertex */
        if ( mesh->info.parTyp & MG_Vert ) {
          for (l=0; l<mesh->info.npar; l++) {

            par = &mesh->info.par[l];
            if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) {
              hmax = par->hmax;
              isloc   = 1;
              break;
            }
          }
        }

        /* Local param at tetrahedra */
        if ( mesh->info.parTyp & MG_Tetra ) {
          ilistv = MMG5_boulevolp(mesh,k,i,listv);
          l = 0;
          do
          {
            if ( isloc )  break;

            par = &mesh->info.par[l];
            if ( par->elt != MMG5_Tetrahedron )  continue;

            for ( kk=0; kk<ilistv; ++kk ) {
              ptloc = &mesh->tetra[listv[kk]/4];
              if ( par->ref == ptloc->ref ) {
                hmax = par->hmax;
                isloc   = 1;
                break;
              }
            }
          } while ( ++l<mesh->info.npar );

          for ( ; l<mesh->info.npar; ++l ) {
            par = &mesh->info.par[l];
            if ( par->elt != MMG5_Tetrahedron ) continue;

            for ( kk=0; kk<ilistv; ++kk ) {
              ptloc = &mesh->tetra[listv[kk]/4];
              if ( par->ref == ptloc->ref ) {
                hmax = MG_MIN(hmax,par->hmax);
                break;
              }
            }
          }
        }
        /** Second step: set the metric */
        met->m[ip0] = hmax;
        p0->flag    = 1;
      }
    }

    /** Set size at points that cannot be reached from the tetra */
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) )  continue;

      for (i=0; i<6; i++) {
        ip0 = pp->v[i];
        p0  = &mesh->point[ip0];

        if ( p0->flag ) continue;

        met->m[ip0] = hmax;
        p0->flag    = 1;
      }
    }
  }
  else {

    /* size truncation */
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<4; i++) {
        ip0 = pt->v[i];
        p0  = &mesh->point[ip0];

        if ( p0->flag ) continue;

        /** First step: search for local parameters */
        isloc   = 0;
        hmin = mesh->info.hmin;
        hmax = mesh->info.hmax;

        /* Local param at vertex */
        if ( mesh->info.parTyp & MG_Vert ) {
          for (l=0; l<mesh->info.npar; l++) {
            par = &mesh->info.par[l];
            if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) {
              hmin = par->hmin;
              hmax = par->hmax;
              isloc   = 1;
              break;
            }
          }
        }

        /* Local param ar tetrahedra */
        if ( mesh->info.parTyp & MG_Tetra ) {
          ilistv = MMG5_boulevolp(mesh,k,i,listv);
          l = 0;
          do
          {
            if ( isloc )  break;

            par = &mesh->info.par[l];
            if ( par->elt != MMG5_Tetrahedron )  continue;

            for ( kk=0; kk<ilistv; ++kk ) {
              ptloc = &mesh->tetra[listv[kk]/4];
              if ( par->ref == ptloc->ref ) {
                hmin = par->hmin;
                hmax = par->hmax;
                isloc   = 1;
                break;
              }
            }
          } while ( ++l<mesh->info.npar );

          for ( ; l<mesh->info.npar; ++l ) {
            par = &mesh->info.par[l];
            if ( par->elt != MMG5_Tetrahedron ) continue;

            for ( kk=0; kk<ilistv; ++kk ) {
              ptloc = &mesh->tetra[listv[kk]/4];
              if ( par->ref == ptloc->ref ) {
                hmin = MG_MAX(hmin,par->hmin);
                hmax = MG_MIN(hmax,par->hmax);
                break;
              }
            }
          }
        }
        /** Second step: set the metric */
        met->m[ip0] = MG_MIN(hmax,MG_MAX(hmin,met->m[ip0]));
        p0->flag    = 1;
      }
    }

    /** Set size at points that cannot be reached from the tetra */
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) )  continue;

      for (i=0; i<6; i++) {
        ip0 = pp->v[i];
        p0  = &mesh->point[ip0];

        if ( p0->flag ) continue;

        met->m[ip0] = MG_MIN(hmax,MG_MAX(hmin,met->m[ip0]));
        p0->flag    = 1;
      }
    }
  }

  /** Step 3: size at regular surface points */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    // Warning: why are we skipped the tetra with negative refs ?
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    else if ( !pt->xt )  continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( !(pxt->ftag[i] & MG_BDY) ) continue;
      if ( !MG_GET(mesh->xtetra[mesh->tetra[k].xt].ori,i) ) continue;

      for (j=0; j<3; j++) {
        i0  = MMG5_idir[i][j];
        ip0 = pt->v[i0];
        p0  = &mesh->point[ip0];

        if ( p0->flag>1 ) continue;

        if ( MG_SIN_OR_NOM(p0->tag) || MG_EDG(p0->tag) )
          continue;

        /** First step: search for local parameters */
        if ( MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0) != 1 )
          continue;

        if ( !MMG3D_localParamReg(mesh,ip0,listv,ilistv,lists,ilists,
                                   &hausd,&hmin,&hmax) ) {
          hmin = mesh->info.hmin;
          hmax = mesh->info.hmax;
          hausd = mesh->info.hausd;
        }

        /** Second step: set the metric */
        /* Define size coming from the hausdorff approximation at regular
         * surface point and update metric at non-singular non-manifold points
         * of surfacic ball*/
        hp  = MMG5_defsizreg(mesh,met,ip0,lists,ilists,hmin,hmax,hausd);

        met->m[ip0] = MG_MIN(met->m[ip0],hp);
        p0->flag = 2;
      }
    }
  }

  /** 3) Travel all boundary faces to update size prescription for points on
   * ridges/edges */
  /* Warning: here we pass more than once per each point because we see it from
     all the edges to which it belongs */
  for (k=1; k<=mesh->ne; k++) {

    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    else if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++) {
      if ( !(pxt->ftag[i] & MG_BDY) )  continue;
      else if ( !MMG5_norface(mesh,k,i,v) )  continue;

      for (j=0; j<3; j++) {
        ia = MMG5_iarf[i][j];
        i0 = MMG5_iare[ia][0];
        i1 = MMG5_iare[ia][1];
        ip0 = pt->v[i0];
        ip1 = pt->v[i1];
        p0  = &mesh->point[ip0];
        p1  = &mesh->point[ip1];

        /* Skip this step if both points are on a required edge */
        if ( p0->flag == 3 && p1->flag == 3 ) continue;

        /* Skip regular edges */
        if ( !MG_EDG(p0->tag) && !MG_EDG(p1->tag) )  continue;

        /** First step: search for local parameters */
        if ( !MMG3D_localParamNm(mesh,k,i,ia,&hausd,&hmin,&hmax) ) {
          hausd = mesh->info.hausd;
          hmin  = mesh->info.hmin;
          hmax  = mesh->info.hmax;
        }

        /** Second step: set metric */
        ised = MG_EDG_OR_NOM(pxt->tag[ia]);

        MMG5_BezierEdge(mesh,ip0,ip1,b0,b1,ised,v);

        b0p0[0] = b0[0] - p0->c[0];
        b0p0[1] = b0[1] - p0->c[1];
        b0p0[2] = b0[2] - p0->c[2];

        b1b0[0] = b1[0] - b0[0];
        b1b0[1] = b1[1] - b0[1];
        b1b0[2] = b1[2] - b0[2];

        p1b1[0] = p1->c[0] - b1[0];
        p1b1[1] = p1->c[1] - b1[1];
        p1b1[2] = p1->c[2] - b1[2];

        secder0[0] = p0->c[0] + b1[0] - 2.0*b0[0];
        secder0[1] = p0->c[1] + b1[1] - 2.0*b0[1];
        secder0[2] = p0->c[2] + b1[2] - 2.0*b0[2];

        secder1[0] = p1->c[0] + b0[0] - 2.0*b1[0];
        secder1[1] = p1->c[1] + b0[1] - 2.0*b1[1];
        secder1[2] = p1->c[2] + b0[2] - 2.0*b1[2];

        kappa = 0.0;
        for (l=0; l<4; l++) {
          tau[0] = 3.0*(1.0-MMG5_ATHIRD*l)*(1.0-MMG5_ATHIRD*l)*b0p0[0] + 6.0*MMG5_ATHIRD*l*(1.0-MMG5_ATHIRD*l)*b1b0[0] \
            + 3.0*MMG5_ATHIRD*l*MMG5_ATHIRD*l*p1b1[0];
          tau[1] = 3.0*(1.0-MMG5_ATHIRD*l)*(1.0-MMG5_ATHIRD*l)*b0p0[1] + 6.0*MMG5_ATHIRD*l*(1.0-MMG5_ATHIRD*l)*b1b0[1] \
            + 3.0*MMG5_ATHIRD*l*MMG5_ATHIRD*l*p1b1[1];
          tau[2] = 3.0*(1.0-MMG5_ATHIRD*l)*(1.0-MMG5_ATHIRD*l)*b0p0[2] + 6.0*MMG5_ATHIRD*l*(1.0-MMG5_ATHIRD*l)*b1b0[2] \
            + 3.0*MMG5_ATHIRD*l*MMG5_ATHIRD*l*p1b1[2];

          gammasec[0] = 6.0*((1.0-MMG5_ATHIRD*l)*secder0[0] + MMG5_ATHIRD*l*secder1[0]);
          gammasec[1] = 6.0*((1.0-MMG5_ATHIRD*l)*secder0[1] + MMG5_ATHIRD*l*secder1[1]);
          gammasec[2] = 6.0*((1.0-MMG5_ATHIRD*l)*secder0[2] + MMG5_ATHIRD*l*secder1[2]);

          ntau2 = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
          if ( ntau2 < MMG5_EPSD )  continue;
          intau = 1.0/sqrt(ntau2);
          ntau2 = 1.0/ntau2;
          tau[0] *= intau;
          tau[1] *= intau;
          tau[2] *= intau;

          ps = gammasec[0]*tau[0] + gammasec[1]*tau[1] + gammasec[2]*tau[2];
          gammasec[0] = gammasec[0]*ntau2 - ps*ntau2*tau[0];
          gammasec[1] = gammasec[1]*ntau2 - ps*ntau2*tau[1];
          gammasec[2] = gammasec[2]*ntau2 - ps*ntau2*tau[2];
          kappa = MG_MAX(kappa,gammasec[0]*gammasec[0] + gammasec[1]*gammasec[1] + gammasec[2]*gammasec[2] );
        }
        kappa = sqrt(kappa);
        if ( kappa < MMG5_EPSD )
          lm = MAXLEN;
        else
          lm = sqrt(8.0*hausd / kappa);

        if ( MG_EDG(p0->tag) && !MG_SIN_OR_NOM(p0->tag) && p0->flag != 3 )
          met->m[ip0] = MG_MAX(hmin,MG_MIN(met->m[ip0],lm));
        if ( MG_EDG(p1->tag) && !MG_SIN_OR_NOM(p1->tag) && p1->flag != 3 )
          met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Set the s field of the points that belongs to a required edge to 4*ne+3, set it to
 * 0 otherwise.
 *
 */
void MMG3D_mark_pointsOnReqEdge_fromTetra (  MMG5_pMesh mesh ) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt;
  MMG5_int     k;
  int8_t       i;

  for ( k=1; k<=mesh->np; k++ ) {
    ppt = &mesh->point[k];
    ppt->s = 0;
  }

  /* Mark the points that belongs to a required edge */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( (!MG_EOK(pt)) || !pt->xt ) { continue; }

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<6; i++) {
      if ( pxt->tag[i] & MG_REQ ) {
        mesh->point[pt->v[MMG5_iare[i][0]]].s = 4*mesh->ne+3;
        mesh->point[pt->v[MMG5_iare[i][1]]].s = 4*mesh->ne+3;
      }
    }
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Enforce mesh gradation by truncating size map.
 *
 */
int MMG3D_gradsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra    pt;
  MMG5_pPoint    p0,p1;
  double         l,hn,ux,uy,uz;
  int            it,maxit;
  MMG5_int       ip0,ip1,k,nu,nup;
  int8_t         i,j,ia,i0,i1;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Grading mesh\n");

  MMG3D_mark_pointsOnReqEdge_fromTetra ( mesh );

  for (k=1; k<=mesh->np; k++) {
    mesh->point[k].flag = mesh->base;
  }

  it = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;

      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          ia  = MMG5_iarf[i][j];
          i0  = MMG5_iare[ia][0];
          i1  = MMG5_iare[ia][1];
          ip0 = pt->v[i0];
          ip1 = pt->v[i1];
          p0  = &mesh->point[ip0];
          p1  = &mesh->point[ip1];
          if ( p0->flag < mesh->base-1 && p1->flag < mesh->base-1 )  continue;

          /* Skip points belonging to a required edge */
          if ( p0->s || p1->s ) continue;

          ux = p1->c[0]-p0->c[0];
          uy = p1->c[1]-p0->c[1];
          uz = p1->c[2]-p0->c[2];

          l = ux*ux + uy*uy + uz*uz;
          l = sqrt(l);

          if ( met->m[ip0] < met->m[ip1] ) {
            if ( met->m[ip0] < MMG5_EPSD )  continue;
            hn = met->m[ip0] + mesh->info.hgrad*l;
            if ( met->m[ip1] > hn ) {
              met->m[ip1] = hn;
              p1->flag = mesh->base;
              nu++;
            }
          }
          else {
            if ( met->m[ip1] < MMG5_EPSD )  continue;
            hn = met->m[ip1] + mesh->info.hgrad*l;
            if ( met->m[ip0] > hn ) {
              met->m[ip0] = hn;
              p0->flag = mesh->base;
              nu++;
            }
          }
        }
      }
    }
    nup += nu;
  }
  while( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     gradation: %7" MMG5_PRId " updated, %d iter.\n",nup,it);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Enforce mesh gradation by truncating size map.
 *
 */
int MMG3D_gradsizreq_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra    pt;
  MMG5_pPoint    p0,p1;
  double         hgrad,ll,h0,h1,hn,ux,uy,uz;
  int            it,maxit;
  MMG5_int       ip0,ip1,ipslave,ipmaster,k,nu,nup;
  int8_t         i,j,ia,i0,i1;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  ** Grading required points.\n");
  }

  if ( mesh->info.hgrad < 0. ) {
    /** Mark the edges belonging to a required entity */
    MMG3D_mark_pointsOnReqEdge_fromTetra ( mesh );
  }

  /** Update the sizes and mark the treated points */
  hgrad = mesh->info.hgradreq;
  it = nup = 0;
  maxit = 100;

  do {
    nu = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) {
        continue;
      }

      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          ia  = MMG5_iarf[i][j];
          i0  = MMG5_iare[ia][0];
          i1  = MMG5_iare[ia][1];
          ip0 = pt->v[i0];
          ip1 = pt->v[i1];
          p0  = &mesh->point[ip0];
          p1  = &mesh->point[ip1];

          if ( MMG5_abs ( p0->s - p1->s ) < 2 ) {
            /* No size to propagate */
            continue;
          }
          else if ( p0->s > p1->s ) {
            ipmaster = ip0;
            ipslave  = ip1;
          }
          else {
            assert ( p1->s > p0->s );
            ipmaster = ip1;
            ipslave  = ip0;
          }

          ux = p1->c[0]-p0->c[0];
          uy = p1->c[1]-p0->c[1];
          uz = p1->c[2]-p0->c[2];

          ll = ux*ux + uy*uy + uz*uz;
          ll = sqrt(ll);

          h0 = met->m[ipmaster];
          h1 = met->m[ipslave];
          if ( h0 < h1 ) {
            if ( h0 < MMG5_EPSD ) {
              continue;
            }
            hn  = h0 + hgrad*ll;
            if ( h1 <= hn ) {
              continue;
            }
          }
          else {
            hn = h0 - hgrad*ll;
            if ( h1 >= hn ) {
              continue;
            }
          }
          met->m[ipslave]           = hn;
          mesh->point[ipslave].s    = mesh->point[ipmaster].s - 1;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 && nup ) {
    fprintf(stdout,"     gradation (required): %7" MMG5_PRId " updated, %d iter.\n",nup,it);
  }

  return nup;
}
