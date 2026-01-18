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
 * \file common/intmet.c
 * \brief Functions to compute metric interpolation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"
#include "mmgexterns_private.h"

/**
 * \param m input metric.
 * \param n input metric.
 * \param mr computed output metric.
 * \param s parameter coordinate for the interpolation of metrics \a m and \a n.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute the interpolated \f$(3 x 3)\f$ metric from metrics \a m and \a n, at
 * parameter \a s : \f$ mr = (1-s)*m +s*n \f$, both metrics being expressed in
 * the simultaneous reduction basis: linear interpolation of sizes.
 *
 */
int MMG5_mmgIntmet33_ani(double *m,double *n,double *mr,double s) {
  int     order;
  double  lambda[3],vp[3][3],mu[3],is[6],isnis[6],mt[9],P[9],dd;
  int8_t  i;
  static int8_t mmgWarn=0;

  /* Compute inverse of square root of matrix M : is =
   * P*diag(1/sqrt(lambda))*{^t}P */
  order = MMG5_eigenv3d(1,m,lambda,vp);
  if ( !order ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to diagonalize at least"
              " 1 metric.\n",__func__);
      mmgWarn = 1;
    }
    return 0;
  }

  for (i=0; i<3; i++) {
    if ( lambda[i] < MMG5_EPSD ) return 0;
    lambda[i] = sqrt(lambda[i]);
    lambda[i] = 1.0 / lambda[i];
  }

  is[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0]
    + lambda[2]*vp[2][0]*vp[2][0];
  is[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1]
    + lambda[2]*vp[2][0]*vp[2][1];
  is[2] = lambda[0]*vp[0][0]*vp[0][2] + lambda[1]*vp[1][0]*vp[1][2]
    + lambda[2]*vp[2][0]*vp[2][2];
  is[3] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1]
    + lambda[2]*vp[2][1]*vp[2][1];
  is[4] = lambda[0]*vp[0][1]*vp[0][2] + lambda[1]*vp[1][1]*vp[1][2]
    + lambda[2]*vp[2][1]*vp[2][2];
  is[5] = lambda[0]*vp[0][2]*vp[0][2] + lambda[1]*vp[1][2]*vp[1][2]
    + lambda[2]*vp[2][2]*vp[2][2];

  mt[0] = n[0]*is[0] + n[1]*is[1] + n[2]*is[2];
  mt[1] = n[0]*is[1] + n[1]*is[3] + n[2]*is[4];
  mt[2] = n[0]*is[2] + n[1]*is[4] + n[2]*is[5];
  mt[3] = n[1]*is[0] + n[3]*is[1] + n[4]*is[2];
  mt[4] = n[1]*is[1] + n[3]*is[3] + n[4]*is[4];
  mt[5] = n[1]*is[2] + n[3]*is[4] + n[4]*is[5];
  mt[6] = n[2]*is[0] + n[4]*is[1] + n[5]*is[2];
  mt[7] = n[2]*is[1] + n[4]*is[3] + n[5]*is[4];
  mt[8] = n[2]*is[2] + n[4]*is[4] + n[5]*is[5];

  isnis[0] = is[0]*mt[0] + is[1]*mt[3] + is[2]*mt[6];
  isnis[1] = is[0]*mt[1] + is[1]*mt[4] + is[2]*mt[7];
  isnis[2] = is[0]*mt[2] + is[1]*mt[5] + is[2]*mt[8];
  isnis[3] = is[1]*mt[1] + is[3]*mt[4] + is[4]*mt[7];
  isnis[4] = is[1]*mt[2] + is[3]*mt[5] + is[4]*mt[8];
  isnis[5] = is[2]*mt[2] + is[4]*mt[5] + is[5]*mt[8];

  order = MMG5_eigenv3d(1,isnis,lambda,vp);
  if ( !order ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to diagonalize at least"
              " 1 metric.\n",__func__);
      mmgWarn = 1;
    }
    return 0;
  }

  /* P = is * (vp) */
  P[0] = is[0]*vp[0][0] + is[1]*vp[0][1] + is[2]*vp[0][2];
  P[1] = is[0]*vp[1][0] + is[1]*vp[1][1] + is[2]*vp[1][2];
  P[2] = is[0]*vp[2][0] + is[1]*vp[2][1] + is[2]*vp[2][2];
  P[3] = is[1]*vp[0][0] + is[3]*vp[0][1] + is[4]*vp[0][2];
  P[4] = is[1]*vp[1][0] + is[3]*vp[1][1] + is[4]*vp[1][2];
  P[5] = is[1]*vp[2][0] + is[3]*vp[2][1] + is[4]*vp[2][2];
  P[6] = is[2]*vp[0][0] + is[4]*vp[0][1] + is[5]*vp[0][2];
  P[7] = is[2]*vp[1][0] + is[4]*vp[1][1] + is[5]*vp[1][2];
  P[8] = is[2]*vp[2][0] + is[4]*vp[2][1] + is[5]*vp[2][2];

  /* At this point, theory states that ^tPMP = I, {^t}PNP=\Lambda */
  /* Linear interpolation between sizes */
  for(i=0; i<3; i++) {
    if ( lambda[i] < 0.0 ) return 0;
    dd = s*sqrt(lambda[i]) + (1.0-s);
    dd = dd*dd;
    if ( dd < MMG5_EPSD )  return 0;
    mu[i] = lambda[i]/dd;
  }

  if ( !MMG5_invmatg(P,mt) )  return 0;

  /* Resulting matrix = ^tP^{-1} diag(mu) P^{-1} */
  mr[0] = mu[0]*mt[0]*mt[0] + mu[1]*mt[3]*mt[3] + mu[2]*mt[6]*mt[6];
  mr[1] = mu[0]*mt[0]*mt[1] + mu[1]*mt[3]*mt[4] + mu[2]*mt[6]*mt[7];
  mr[2] = mu[0]*mt[0]*mt[2] + mu[1]*mt[3]*mt[5] + mu[2]*mt[6]*mt[8];
  mr[3] = mu[0]*mt[1]*mt[1] + mu[1]*mt[4]*mt[4] + mu[2]*mt[7]*mt[7];
  mr[4] = mu[0]*mt[1]*mt[2] + mu[1]*mt[4]*mt[5] + mu[2]*mt[7]*mt[8];
  mr[5] = mu[0]*mt[2]*mt[2] + mu[1]*mt[5]*mt[5] + mu[2]*mt[8]*mt[8];

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param ip1 global index of ridge extremity.
 * \param ip2 global index of ridge extremity.
 * \param s interpolation parameter (between 0 and 1).
 * \param v normal at the point at which we want to compute the metric.
 * \param mr computed anisotropic size.
 * \return 1 if success, 0 otherwise.
 *
 * Anisotropic metric interpolation between two points \f$p_1\f$ and \f$p_2\f$
 * such that \f$edge_0 = (p_1p_2)\f$ is ridge. \a v is a direction vector, aimed
 * at pointing towards direction of n1 at interpolated point.
 *
 */
int MMG5_intridmet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int ip1, MMG5_int ip2,double s,
                    double v[3],double mr[6]) {
  MMG5_pxPoint   go1,go2;
  MMG5_pPoint    p1,p2;
  double         *m1,*m2,*n11,*n12,*n21,*n22,ps11,ps12,dd;
  double         hu1,hu2,hn1,hn2;

  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];
  m1  = &met->m[6*ip1];
  m2  = &met->m[6*ip2];

  /* Case when both endpoints are singular */
  if ( (MG_SIN(p1->tag) || (p1->tag & MG_NOM)) &&
       (MG_SIN(p2->tag) || (p2->tag & MG_NOM)) ) {
    /* m1 and m2 are isotropic metrics */
    /* Remark (Algiane): Perspective of improvement can be:
       - 1. to not force isotropy at singular point
       - 2. to not force the metric to be along tangent and normal dir at ridge point
    */

    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      if ( s < 0.5 ) {
        mr[0] = m1[0];
        mr[1] = m1[0];
        mr[2] = m1[0];
        mr[3] = m1[0];
        mr[4] = m1[0];
      }
      else {
        mr[0] = m2[0];
        mr[1] = m2[0];
        mr[2] = m2[0];
        mr[3] = m2[0];
        mr[4] = m2[0];
      }
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
      mr[1] = mr[0];
      mr[2] = mr[0];
      mr[3] = mr[0];
      mr[4] = mr[0];
    }
  }
  /* vertex p1 is singular, p2 is regular */
  else if (  (MG_SIN(p1->tag) || (p1->tag & MG_NOM)) &&
            !(MG_SIN(p2->tag) || (p2->tag & MG_NOM)) ) {
    /* m1 is an isotropic metric and m2 is a "ridge" metric that respect our
     * storage convention. */
    go2 = &mesh->xpoint[p2->xp];
    n21 = &go2->n1[0];
    n22 = &go2->n2[0];

    /* Interpolation of the eigenvalue associated to tangent vector */
    dd = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      mr[0] = s < 0.5 ? m1[0] : m2[0];
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
    }

    /* Interpolation of the two other eigenvalues for each configuration. */
    /* 1. For the surface ruled by n1. */
    dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hu1 = s < 0.5 ? m1[0] : m2[1];
    }
    else {
      hu1 = m1[0]*m2[1] / dd;
    }
    dd  = (1-s)*sqrt(m2[3]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hn1 = s < 0.5 ? m1[0] : m2[3];
    }
    else {
      hn1 = m1[0]*m2[3] / dd;
    }

    /* 2. For the surface ruled by n2. */
    dd = (1-s)*sqrt(m2[2]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hu2 = s < 0.5 ? m1[0] : m2[2];
    }
    else {
      hu2 = m1[0]*m2[2] / dd;
    }
    dd  = (1-s)*sqrt(m2[4]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hn2 = s < 0.5 ? m1[0] : m2[4];
    }
    else {
      hn2 = m1[0]*m2[4] / dd;
    }

    /* Decision of the ordering of hu1 and hu2 */
    ps11 = n21[0]*v[0] + n21[1]*v[1] + n21[2]*v[2];
    ps12 = n22[0]*v[0] + n22[1]*v[1] + n22[2]*v[2];
    if ( fabs(ps11) > fabs(ps12) ) {
      mr[1] = hu1;
      mr[2] = hu2;
      mr[3] = hn1;
      mr[4] = hn2;
    }
    else {
      mr[1] = hu2;
      mr[2] = hu1;
      mr[3] = hn2;
      mr[4] = hn1;
    }
  }
  /* vertex p2 is singular, p1 is regular */
  else if ( ( MG_SIN(p2->tag) || (p2->tag & MG_NOM)) &&
            !(MG_SIN(p1->tag) || (p1->tag & MG_NOM)) ) {
    /* m2 is an isotropic metric and m1 is a "ridge" metric that respect our
     * storage convention. */
    go1 = &mesh->xpoint[p1->xp];
    n11 = &go1->n1[0];
    n12 = &go1->n2[0];

    /* Interpolation of the eigenvalue associated to tangent vector */
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      mr[0] = s < 0.5 ? m1[0] : m2[0];
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
    }
    /* Interpolation of the two other eigenvalues for each configuration. */
    /* 1. For the surface ruled by n1. */
    dd = (1-s)*sqrt(m2[0]) + s*sqrt(m1[1]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hu1 = s < 0.5 ? m1[1] : m2[0];
    }
    else {
      hu1 = m1[1]*m2[0] / dd;
    }
    dd = (1-s)*sqrt(m2[0]) + s*sqrt(m1[3]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hn1 = s < 0.5 ? m1[3] : m2[0];
    }
    else {
      hn1 = m1[3]*m2[0] / dd;
    }

    /* 2. For the surface ruled by n2. */
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[2]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hu2 = s < 0.5 ? m1[2] : m2[0];
    }
    else {
      hu2 = m1[2]*m2[0] / dd;
    }
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[4]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      hn2 = s < 0.5 ? m1[4] : m2[0];
    }
    else {
      hn2 = m1[4]*m2[0] / dd;
    }

    /* Decision of the ordering of hu1 and hu2 */
    ps11 = n11[0]*v[0] + n11[1]*v[1] + n11[2]*v[2];
    ps12 = n12[0]*v[0] + n12[1]*v[1] + n12[2]*v[2];
    if ( fabs(ps11) > fabs(ps12) ) {
      mr[1] = hu1;
      mr[2] = hu2;
      mr[3] = hn1;
      mr[4] = hn2;
    }
    else {
      mr[1] = hu2;
      mr[2] = hu1;
      mr[3] = hn2;
      mr[4] = hn1;
    }
  }
  /* p1,p2 : nonsingular vertices */
  else {
    go1 = &mesh->xpoint[p1->xp];
    go2 = &mesh->xpoint[p2->xp];

    /* Interpolation of the eigenvalue associated to tangent vector */
    dd  = (1-s)*sqrt(m2[0]) + s*sqrt(m1[0]);
    dd *= dd;
    if ( dd < MMG5_EPSD ) {
      mr[0] = s < 0.5 ? m1[0] : m2[0];
    }
    else {
      mr[0] = m1[0]*m2[0] / dd;
    }

    /* Pairing of normal vectors at p1 and p2 */
    n11 = &go1->n1[0];
    n12 = &go1->n2[0];
    n21 = &go2->n1[0];
    n22 = &go2->n2[0];
    ps11 = n11[0]*n21[0] + n11[1]*n21[1] + n11[2]*n21[2];
    ps12 = n11[0]*n22[0] + n11[1]*n22[1] + n11[2]*n22[2];
    if ( fabs(ps11) > fabs(ps12) ) {   //n11 and n21 go together
      /* 1. For the surface ruled by n1. */
      dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[1]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hu1 = s < 0.5 ? m1[1] : m2[1];
      }
      else {
        hu1 = m1[1]*m2[1] / dd;
      }
      dd  = (1-s)*sqrt(m2[3]) + s*sqrt(m1[3]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hn1 = s < 0.5 ? m1[3] : m2[3];
      }
      else {
        hn1 = m1[3]*m2[3] / dd;
      }
      /* 2. For the surface ruled by n2. */
      dd = (1-s)*sqrt(m2[2]) + s*sqrt(m1[2]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hu2 = s < 0.5 ? m1[2] : m2[2];
      }
      else {
        hu2 = m1[2]*m2[2] / dd;
      }
      dd = (1-s)*sqrt(m2[4]) + s*sqrt(m1[4]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hn2 = s < 0.5 ? m1[4] : m2[4];
      }
      else {
        hn2 = m1[4]*m2[4] / dd;
      }

    }
    else {
      /* 1. */
      dd  = (1-s)*sqrt(m2[2]) + s*sqrt(m1[1]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hu1 = s < 0.5 ? m1[1] : m2[2];
      }
      else {
        hu1 = m1[1]*m2[2] / dd;
      }
      dd  = (1-s)*sqrt(m2[4]) + s*sqrt(m1[3]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hn1 = s < 0.5 ? m1[3] : m2[4];
      }
      else {
        hn1 = m1[3]*m2[4] / dd;
      }

      /* 2. */
      dd  = (1-s)*sqrt(m2[1]) + s*sqrt(m1[2]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hu2 = s < 0.5 ? m1[2] : m2[1];
      }
      else {
        hu2 = m1[2]*m2[1] / dd;
      }
      dd  = (1-s)*sqrt(m2[3]) + s*sqrt(m1[4]);
      dd *= dd;
      if ( dd < MMG5_EPSD ) {
        hn2 = s < 0.5 ? m1[4] : m2[3];
      }
      else {
        hn2 = m1[4]*m2[3] / dd;
      }
    }

    /* Now, hu1 is the eigenvalue associated to the direction at interpolated
       point, closest to n11 (hu2 -> n12) ; one may need a different
       orientation, and put eigenvalue of direction closest to v (= interpolated
       normal) first */
    ps11 = n11[0]*v[0] + n11[1]*v[1] + n11[2]*v[2];
    ps12 = n12[0]*v[0] + n12[1]*v[1] + n12[2]*v[2];
    if ( fabs(ps11) > fabs(ps12) ) {
      mr[1] = hu1;
      mr[2] = hu2;
      mr[3] = hn1;
      mr[4] = hn2;
    }
    else {
      mr[1] = hu2;
      mr[2] = hu1;
      mr[3] = hn2;
      mr[4] = hn1;
    }
  }
  mr[5] = 0.0;

  return 1;
}

/**
 * \param ma pointer on a metric
 * \param mb pointer on a metric
 * \param mp pointer on the computed interpolated metric
 * \param t interpolation parameter (comprise between 0 and 1)
 *
 *
 * Linear interpolation of isotropic sizemap along an edge
 *
 */
int MMG5_interp_iso(double *ma,double *mb,double *mp,double t) {

  *mp = (1.0-t)*(*ma) + t*(*mb);

  return 1;

}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param pt pointer to the triangle structure.
 * \param i edge of the triangle pt
 * \param s interpolated parameter (comprise between 0 and 1)
 * \param mr computed interpolated metric
 * \return 0 if fail, 1 otherwise.
 *
 * Metric interpolation between points p1 and p2, in tria \a pt at parameter 0
 * <= \a s <= 1 from p1 result is stored in \a mr. edge p1p2 must not be a ridge
 */
int MMG5_interpreg_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt,int8_t i,
                        double s,double mr[6]) {
  MMG5_pPoint    p1,p2;
  MMG5_Bezier    b;
  double         b1[3],b2[3],bn[3],c[3],nt[3],cold[3],nold[3],n[3];
  double         m1old[6],m2old[6],m1[6],m2[6],rbasis[3][3];
  double         *n1,*n2,step,u,r[3][3],dd,ddbn;
  int            nstep,l;
  MMG5_int       ip1,ip2;
  int8_t         i1,i2;
  static int     warn=0,warnnorm=0;

  /* Number of steps for parallel transport */
  nstep = 4;
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];

  if ( !MMG5_bezierCP(mesh,pt,&b,1) )  return 0;

  n1 = &b.n[i1][0];
  n2 = &b.n[i2][0];
  memcpy(bn,&b.n[i+3][0],3*sizeof(double));
  memcpy(b1,&b.b[2*i+3][0],3*sizeof(double));
  memcpy(b2,&b.b[2*i+4][0],3*sizeof(double));

  ddbn = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];

  /* Parallel transport of metric at p1 to point p(s) */
  step = s / nstep;
  cold[0] = p1->c[0];
  cold[1] = p1->c[1];
  cold[2] = p1->c[2];

  nold[0] = n1[0];
  nold[1] = n1[1];
  nold[2] = n1[2];

  if ( MG_SIN(p1->tag) || (p1->tag & MG_NOM) || (ddbn < MMG5_EPSD) ) {
    memcpy(m1,&met->m[6*ip1],6*sizeof(double));
  }
  else {
    if ( MG_GEO & p1->tag ) {
      MMG5_nortri(mesh,pt,nt);
      if ( !MMG5_buildridmetnor(mesh,met,pt->v[i1],nt,m1,rbasis) )  return 0;
    }
    else {
      memcpy(m1,&met->m[6*ip1],6*sizeof(double));
    }
    memcpy(m1old,m1,6*sizeof(double));

    /* Go from point (l-1)step, to point l step */
    for (l=1; l<=nstep; l++) {
      u    = l*step;
      c[0] = p1->c[0]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[0]\
        + 3.0*u*u*(1.0-u)*b2[0] + u*u*u*p2->c[0];
      c[1] = p1->c[1]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[1]\
        + 3.0*u*u*(1.0-u)*b2[1] + u*u*u*p2->c[1];
      c[2] = p1->c[2]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[2]\
        + 3.0*u*u*(1.0-u)*b2[2] + u*u*u*p2->c[2];

      n[0] = (1.0-u)*(1.0-u)*n1[0] + 2.0*u*(1.0-u)*bn[0] + u*u*n2[0];
      n[1] = (1.0-u)*(1.0-u)*n1[1] + 2.0*u*(1.0-u)*bn[1] + u*u*n2[1];
      n[2] = (1.0-u)*(1.0-u)*n1[2] + 2.0*u*(1.0-u)*bn[2] + u*u*n2[2];
      dd   = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd < MMG5_EPSD )  return 0;
      dd = 1.0 / sqrt(dd);
      n[0] *= dd;
      n[1] *= dd;
      n[2] *= dd;

      if ( !MMG5_paratmet(cold,nold,m1old,c,n,m1) )  return 0;

      memcpy(cold,c,3*sizeof(double));
      memcpy(nold,n,3*sizeof(double));
      memcpy(m1old,m1,6*sizeof(double));
    }
  }

  /* Parallel transport of metric at p2 to point p(s) */
  step = (1.0-s) / nstep;
  cold[0] = p2->c[0];
  cold[1] = p2->c[1];
  cold[2] = p2->c[2];

  nold[0] = n2[0];
  nold[1] = n2[1];
  nold[2] = n2[2];

  if ( MG_SIN(p2->tag) || (p2->tag & MG_NOM) || (ddbn < MMG5_EPSD) ) {
    memcpy(m2,&met->m[6*ip2],6*sizeof(double));

    /* In this pathological case, n is empty */
    if ( MG_SIN(p1->tag) || (p1->tag & MG_NOM) ) {
      memcpy(n,n2,3*sizeof(double));
      assert( MMG5_EPSD < (n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]) && "normal at p2 is 0" );
    }
    else if (ddbn < MMG5_EPSD) {
      /* Other case where n is empty: bezier normal is 0 */
      if ( !warnnorm ) {
        fprintf(stderr,"  ## Warning: %s: %d: unexpected case (null normal),"
                " impossible interpolation.\n",__func__,__LINE__);
        warnnorm = 1;
      }
      return 0;
    }
  }
  else {
    if ( p2->tag & MG_GEO ) {
      MMG5_nortri(mesh,pt,nt);
      if ( !MMG5_buildridmetnor(mesh,met,pt->v[i2],nt,m2,rbasis))  return 0;
    }
    else {
      memcpy(m2,&met->m[6*ip2],6*sizeof(double));
    }
    memcpy(m2old,m2,6*sizeof(double));

    /* Go from point (l-1)step, to point l step */
    for (l=1; l<=nstep; l++) {
      u    = 1.0 - l*step;
      c[0] = p1->c[0]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[0]\
        + 3.0*u*u*(1.0-u)*b2[0] + u*u*u*p2->c[0];
      c[1] = p1->c[1]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[1]\
        + 3.0*u*u*(1.0-u)*b2[1] + u*u*u*p2->c[1];
      c[2] = p1->c[2]*(1.0-u)*(1.0-u)*(1.0-u) + 3.0*(1.0-u)*(1.0-u)*u*b1[2]\
        + 3.0*u*u*(1.0-u)*b2[2] + u*u*u*p2->c[2];

      n[0] = (1.0-u)*(1.0-u)*n1[0] + 2.0*u*(1.0-u)*bn[0] + u*u*n2[0];
      n[1] = (1.0-u)*(1.0-u)*n1[1] + 2.0*u*(1.0-u)*bn[1] + u*u*n2[1];
      n[2] = (1.0-u)*(1.0-u)*n1[2] + 2.0*u*(1.0-u)*bn[2] + u*u*n2[2];
      dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
      if ( dd < MMG5_EPSD )  return 0;
      dd = 1.0 / sqrt(dd);
      n[0] *= dd;
      n[1] *= dd;
      n[2] *= dd;

      if ( !MMG5_paratmet(cold,nold,m2old,c,n,m2) )  return 0;

      memcpy(cold,c,3*sizeof(double));
      memcpy(nold,n,3*sizeof(double));
      memcpy(m2old,m2,6*sizeof(double));
    }
  }
  /* At this point, c is point p(s), n is the normal at p(s), m1 and m2 are the 3*3
     transported metric tensors from p1 and p2 to p(s) */

  /* Rotate both matrices to the tangent plane */
  if ( !MMG5_rotmatrix(n,r) )  return 0;
  MMG5_rmtr(r,m1,m1old);
  MMG5_rmtr(r,m2,m2old);

  /* Interpolate both metrics expressed in the same tangent plane. */
  if ( !MMG5_mmgIntmet33_ani(m1old,m2old,mr,s) ) {
    if ( !warn ) {
      ++warn;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 impossible metric"
              " interpolation.\n", __func__);

      if ( mesh->info.ddebug ) {
        fprintf(stderr," points: %"MMG5_PRId": %e %e %e (tag %s)\n",MMG5_indPt(mesh,ip1),
                mesh->point[ip1].c[0],mesh->point[ip1].c[1],mesh->point[ip1].c[2],
                MMG5_Get_tagName(mesh->point[ip1].tag));
        fprintf(stderr,"         %"MMG5_PRId": %e %e %e (tag %s)\n",MMG5_indPt(mesh,ip2),
                mesh->point[ip1].c[0],mesh->point[ip1].c[1],mesh->point[ip1].c[2],
                MMG5_Get_tagName(mesh->point[ip2].tag));

        fprintf(stderr,"\n BEFORE ROTATION:\n");
        fprintf(stderr,"\n metric %e %e %e %e %e %e\n",
                m1old[0],m1old[1],m1old[2],m1old[3],m1old[4],m1old[5]);
        fprintf(stderr,"     %e %e %e %e %e %e\n",
                m2old[0],m2old[1],m2old[2],m2old[3],m2old[4],m2old[5]);

        fprintf(stderr,"\n AFTER ROTATION (to %e %e %e):\n",n[0],n[1],n[2]);
        fprintf(stderr,"\n metric %e %e %e %e %e %e\n",
                m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]);
        fprintf(stderr,"     %e %e %e %e %e %e\n",
                m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]);
      }
    }
    return 0;
  }

  /* End rotating back tangent metric into canonical basis of R^3 : mr =
     {^t}R*mr*R m1old serve for nothing, let it be the temporary resulting
     matrix. */
  m1old[0] = mr[0]*r[0][0]*r[0][0] + mr[3]*r[1][0]*r[1][0] + mr[5]*r[2][0]*r[2][0]
    + 2.*( mr[1]*r[0][0]*r[1][0] + mr[2]*r[0][0]*r[2][0] + mr[4]*r[1][0]*r[2][0] );
  m1old[3] = mr[0]*r[0][1]*r[0][1] + mr[3]*r[1][1]*r[1][1] + mr[5]*r[2][1]*r[2][1]
    + 2.*( mr[1]*r[0][1]*r[1][1] + mr[2]*r[0][1]*r[2][1] + mr[4]*r[1][1]*r[2][1] );
  m1old[5] = mr[0]*r[0][2]*r[0][2] + mr[3]*r[1][2]*r[1][2] + mr[5]*r[2][2]*r[2][2]
    + 2.*( mr[1]*r[0][2]*r[1][2] + mr[2]*r[0][2]*r[2][2] + mr[4]*r[1][2]*r[2][2] );

  m1old[1] = mr[0]*r[0][0]*r[0][1] + mr[1]*r[0][0]*r[1][1] + mr[2]*r[0][0]*r[2][1]
    + mr[1]*r[1][0]*r[0][1] + mr[3]*r[1][0]*r[1][1] + mr[4]*r[1][0]*r[2][1]
    + mr[2]*r[2][0]*r[0][1] + mr[4]*r[2][0]*r[1][1] + mr[5]*r[2][0]*r[2][1];
  m1old[2] = mr[0]*r[0][0]*r[0][2] + mr[1]*r[0][0]*r[1][2] + mr[2]*r[0][0]*r[2][2]
    + mr[1]*r[1][0]*r[0][2] + mr[3]*r[1][0]*r[1][2] + mr[4]*r[1][0]*r[2][2]
    + mr[2]*r[2][0]*r[0][2] + mr[4]*r[2][0]*r[1][2] + mr[5]*r[2][0]*r[2][2];
  m1old[4] = mr[0]*r[0][1]*r[0][2] + mr[1]*r[0][1]*r[1][2] + mr[2]*r[0][1]*r[2][2]
    + mr[1]*r[1][1]*r[0][2] + mr[3]*r[1][1]*r[1][2] + mr[4]*r[1][1]*r[2][2]
    + mr[2]*r[2][1]*r[0][2] + mr[4]*r[2][1]*r[1][2] + mr[5]*r[2][1]*r[2][2];

  memcpy(mr,m1old,6*sizeof(double));

  return 1;

}
