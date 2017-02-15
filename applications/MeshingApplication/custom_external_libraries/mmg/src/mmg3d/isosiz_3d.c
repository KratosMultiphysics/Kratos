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
 * \file mmg3d/isosiz_3d.c
 * \brief Fonctions for isotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg3d.h"
#include "inlined_functions.h"

extern char   ddb;

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
inline double _MMG5_lenedgCoor_iso(double *ca,double *cb,double *ma,double *mb) {
  double   h1,h2,l,r,len;

  h1 = *ma;
  h2 = *mb;
  l = (cb[0]-ca[0])*(cb[0]-ca[0]) + (cb[1]-ca[1])*(cb[1]-ca[1]) \
    + (cb[2]-ca[2])*(cb[2]-ca[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < _MMG5_EPS ? l / h1 : l / (h2-h1) * log(r+1.0);

  return(len);
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
 * \return the isotropic size at the point.
 *
 * Define isotropic size at regular point nump, whose surfacic ball is provided.
 *
 */
static double
_MMG5_defsizreg(MMG5_pMesh mesh,MMG5_pSol met,int nump,int *lists,
                int ilists, double hmin,double hmax,double hausd) {
  MMG5_pTetra       pt;
  MMG5_pxTetra      pxt;
  MMG5_pPoint       p0,p1;
  MMG5_Tria         tt;
  _MMG5_Bezier      b;
  double       ux,uy,uz,det2d,h,isqhmin,isqhmax,ll,lmin,lmax,hnm,s;
  double       *n,*t,r[3][3],lispoi[3*MMG3D_LMAX+1],intm[3],b0[3],b1[3],c[3],tAA[6],tAb[3],d[3];
  double       kappa[2],vp[2][2];
  int          k,na,nb,ntempa,ntempb,iel,ip0;
  char         iface,i,j,i0;

  p0 = &mesh->point[nump];

  if ( !p0->xp || MG_EDG(p0->tag) || (p0->tag & MG_NOM) || (p0->tag & MG_REQ))  {
    fprintf(stderr,"    ## Func. _MMG5_defsizreg : wrong point qualification : xp ? %d\n",p0->xp);
    return(0);
  }
  isqhmin = 1.0 / (hmin*hmin);
  isqhmax = 1.0 / (hmax*hmax);

  n = &mesh->xpoint[p0->xp].n1[0];

  /* Step 1 : rotation matrix that sends normal n to the third coordinate vector of R^3 */
  if ( !_MMG5_rotmatrix(n,r) ) {
    fprintf(stderr,"%s:%d: Error: function _MMG5_rotmatrix return 0\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  /* Step 2 : rotation of the oriented surfacic ball with r : lispoi[k] is the common edge
     between faces lists[k-1] and lists[k] */
  iel   = lists[0] / 4;
  iface = lists[0] % 4;
  pt    = &mesh->tetra[iel];
  lmin  = MAXLEN;
  lmax  = 0.0;

  na = nb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != nump ) {
      if ( !na )
        na = pt->v[_MMG5_idir[iface][i]];
      else
        nb = pt->v[_MMG5_idir[iface][i]];
    }
  }

  for (k=1; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    ntempa = ntempb = 0;
    for (i=0; i<3; i++) {
      if ( pt->v[_MMG5_idir[iface][i]] != nump ) {
        if ( !ntempa )
          ntempa = pt->v[_MMG5_idir[iface][i]];
        else
          ntempb = pt->v[_MMG5_idir[iface][i]];
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
    if ( pt->v[_MMG5_idir[iface][i]] != nump ) {
      if ( !ntempa )
        ntempa = pt->v[_MMG5_idir[iface][i]];
      else
        ntempb = pt->v[_MMG5_idir[iface][i]];
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
  if ( lmax/lmin > 4.0*hmax*hmax/(hmin*hmin) )  return(hmax);

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 )  return(hmax);
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 )    return(hmax);

  /* Reconstitution of the curvature tensor at p0 in the tangent plane,
     with a quadric fitting approach */
  memset(intm,0x00,3*sizeof(double));
  memset(tAA,0x00,6*sizeof(double));
  memset(tAb,0x00,3*sizeof(double));

  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;

    _MMG5_tet2tri(mesh,iel,iface,&tt);

    pxt   = &mesh->xtetra[mesh->tetra[iel].xt];
    if ( !_MMG5_bezierCP(mesh,&tt,&b,MG_GET(pxt->ori,iface)) ) {
      fprintf(stderr,"%s:%d: Error: function _MMG5_bezierCP return 0\n",
              __FILE__,__LINE__);
      exit(EXIT_FAILURE);
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
  if ( !_MMG5_sys33sym(tAA,tAb,c) )  return(hmax);

  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
     are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
  if( !_MMG5_eigensym(intm,kappa,vp) ){
    fprintf(stderr,"%s:%d: Error: function _MMG5_eigensym return 0\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
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

  /* Travel surfacic ball one last time and update non manifold point metric */
  for (k=0; k<ilists; k++) {
    iel = lists[k] / 4;
    iface = lists[k] % 4;

    for (j=0; j<3; j++) {
      i0  = _MMG5_idir[iface][j];
      ip0 = pt->v[i0];
      p1  = &mesh->point[ip0];
      if( !(p1->tag & MG_NOM) || MG_SIN(p1->tag) ) continue;
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
  return(h);
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
double _MMG5_meansizreg_iso(MMG5_pMesh mesh,MMG5_pSol met,int nump,int *lists,
                int ilists, double hmin,double hmax) {
  MMG5_pTetra       pt;
  MMG5_pPoint       p0,p1;
  double            len,ux,uy,uz;
  int               k,iel,ip1;
  char              i,iface;

  p0 = &mesh->point[nump];

  len = 0;
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    iface = lists[k] % 4;
    pt    = &mesh->tetra[iel];

    for (i=0; i<3; i++) {
      if ( pt->v[_MMG5_idir[iface][i]] == nump ) {
        break;
      }
    }
    assert(i!=3);

    ip1 = pt->v[_MMG5_idir[iface][_MMG5_inxt2[i]]];
    p1 = &mesh->point[ip1];

    ux  = p1->c[0] - p0->c[0];
    uy  = p1->c[1] - p0->c[1];
    uz  = p1->c[2] - p0->c[2];
    len += sqrt(ux*ux + uy*uy + uz*uz);
  }
  len /=ilists;

  return(MG_MIN(hmax,MG_MAX(hmin,len)));
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
int _MMG3D_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra    pt,ptloc;
  MMG5_pPrism    pp;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    p0,p1;
  double         hp,v[3],b0[3],b1[3],b0p0[3],b1b0[3],p1b1[3],hausd,hmin,hmax;
  double         secder0[3],secder1[3],kappa,tau[3],gammasec[3],ntau2,intau,ps,lm;
  int            lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv,k,ip0,ip1,l;
  int            kk,isloc,ismet;
  char           i,j,ia,ised,i0,i1;
  MMG5_pPar      par;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining isotropic map\n");

  if ( mesh->info.hmax < 0.0 ) {
    //  mesh->info.hmax = 0.5 * mesh->info.delta;
    fprintf(stderr,"%s:%d:Error: negative hmax value.\n",__FILE__,__LINE__);
    return(0);
  }

  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
  }

  /** 1) Size at internal points */
  /* alloc structure */
  if ( !met->m ) {
    ismet      = 0;
    met->np    = mesh->np;
    met->npmax = mesh->npmax;
    met->size  = 1;
    met->dim   = 3;
    _MMG5_ADD_MEM(mesh,(met->npmax+1)*sizeof(double),"solution",return(0));
    _MMG5_SAFE_MALLOC(met->m,(mesh->npmax+1),double);

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
          ilistv = _MMG5_boulevolp(mesh,k,i,listv);
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
    ismet = 1;
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
          ilistv = _MMG5_boulevolp(mesh,k,i,listv);
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

  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
  }

  /** 2) size at regular surface points */
  if ( mesh->info.nosurf && ismet ) return(1);

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
        i0  = _MMG5_idir[i][j];
        ip0 = pt->v[i0];
        p0  = &mesh->point[ip0];

        if ( p0->flag ) continue;

        if ( !mesh->info.nosurf ) {
          if ( MG_SIN(p0->tag) || MG_EDG(p0->tag) || (p0->tag & MG_NOM) )
            continue;
        }
        else
          if ( p0->tag & MG_NOM ) continue;

        /** First step: search for local parameters */
        if ( _MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0) != 1 )
          continue;

        if ( !_MMG3D_localParamReg(mesh,ip0,listv,ilistv,lists,ilists,
                                   &hausd,&hmin,&hmax) ) {
          hmin = mesh->info.hmin;
          hmax = mesh->info.hmax;
          hausd = mesh->info.hausd;
        }

        /** Second step: set the metric */
        if ( !mesh->info.nosurf ) {
          /* Define size coming from the hausdorff approximation at regular
           * surface point */
          hp  = _MMG5_defsizreg(mesh,met,ip0,lists,ilists,hmin,hmax,hausd);
        }
        else {
          /* Define size at regular surface point for the -nosurf option (ie a
           * manifold point): If a metric is provided, we preserve it. Without
           * initial metric, the size is computed as the mean of the length of
           * edges passing through the point */
            hp = _MMG5_meansizreg_iso(mesh,met,ip0,lists,ilists,hmin,hmax);
        }

        met->m[ip0] = MG_MIN(met->m[ip0],hp);
        p0->flag = 1;
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
      else if ( !_MMG5_norface(mesh,k,i,v) )  continue;

      for (j=0; j<3; j++) {
        ia = _MMG5_iarf[i][j];
        i0 = _MMG5_iare[ia][0];
        i1 = _MMG5_iare[ia][1];
        ip0 = pt->v[i0];
        ip1 = pt->v[i1];
        p0  = &mesh->point[ip0];
        p1  = &mesh->point[ip1];

        if ( !MG_EDG(p0->tag) && !MG_EDG(p1->tag) )  continue;

        /** First step: search for local parameters */
        if ( !_MMG3D_localParamNm(mesh,k,i,ia,&hausd,&hmin,&hmax) ) {
          hausd = mesh->info.hausd;
          hmin  = mesh->info.hmin;
          hmax  = mesh->info.hmax;
        }

        /** Second step: set metric */
        if ( !mesh->info.nosurf ) {
          ised = MG_EDG(pxt->tag[ia]) || ( pxt->tag[ia] & MG_NOM );

          _MMG5_BezierEdge(mesh,ip0,ip1,b0,b1,ised,v);

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
            tau[0] = 3.0*(1.0-_MMG5_ATHIRD*l)*(1.0-_MMG5_ATHIRD*l)*b0p0[0] + 6.0*_MMG5_ATHIRD*l*(1.0-_MMG5_ATHIRD*l)*b1b0[0] \
              + 3.0*_MMG5_ATHIRD*l*_MMG5_ATHIRD*l*p1b1[0];
            tau[1] = 3.0*(1.0-_MMG5_ATHIRD*l)*(1.0-_MMG5_ATHIRD*l)*b0p0[1] + 6.0*_MMG5_ATHIRD*l*(1.0-_MMG5_ATHIRD*l)*b1b0[1] \
              + 3.0*_MMG5_ATHIRD*l*_MMG5_ATHIRD*l*p1b1[1];
            tau[2] = 3.0*(1.0-_MMG5_ATHIRD*l)*(1.0-_MMG5_ATHIRD*l)*b0p0[2] + 6.0*_MMG5_ATHIRD*l*(1.0-_MMG5_ATHIRD*l)*b1b0[2] \
              + 3.0*_MMG5_ATHIRD*l*_MMG5_ATHIRD*l*p1b1[2];

            gammasec[0] = 6.0*((1.0-_MMG5_ATHIRD*l)*secder0[0] + _MMG5_ATHIRD*l*secder1[0]);
            gammasec[1] = 6.0*((1.0-_MMG5_ATHIRD*l)*secder0[1] + _MMG5_ATHIRD*l*secder1[1]);
            gammasec[2] = 6.0*((1.0-_MMG5_ATHIRD*l)*secder0[2] + _MMG5_ATHIRD*l*secder1[2]);

            ntau2 = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
            if ( ntau2 < _MMG5_EPSD )  continue;
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
          if ( kappa < _MMG5_EPSD )
            lm = MAXLEN;
          else
            lm = sqrt(8.0*hausd / kappa);

          if ( MG_EDG(p0->tag) && !(p0->tag & MG_NOM) && !MG_SIN(p0->tag) )
            met->m[ip0] = MG_MAX(hmin,MG_MIN(met->m[ip0],lm));
          if ( MG_EDG(p1->tag) && !(p1->tag & MG_NOM) && !MG_SIN(p1->tag) )
            met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
        }
        else {
          if ( !p0->flag ) {
            /* -nosurf option: very rough eval of the metric over non-manifold
             * points: take the non-manifold edge length */
            lm  = (p0->c[0]-p1->c[0])*(p0->c[0]-p1->c[0]);
            lm += (p0->c[1]-p1->c[1])*(p0->c[1]-p1->c[1]);
            lm += (p0->c[2]-p1->c[2])*(p0->c[2]-p1->c[2]);

            lm = sqrt(lm);

            lm = MG_MIN(hmax,MG_MAX(hmin,lm));
            met->m[ip0] = MG_MIN(met->m[ip0],lm);
          }
        }
      }
    }
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Enforce mesh gradation by truncating size map.
 *
 */
int _MMG5_gradsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra    pt;
  MMG5_pPoint    p0,p1;
  double    l,hn;
  int       ip0,ip1,it,maxit,nu,nup,k;
  char      i,j,ia,i0,i1;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Grading mesh\n");

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

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
          ia  = _MMG5_iarf[i][j];
          i0  = _MMG5_iare[ia][0];
          i1  = _MMG5_iare[ia][1];
          ip0 = pt->v[i0];
          ip1 = pt->v[i1];
          p0  = &mesh->point[ip0];
          p1  = &mesh->point[ip1];
          if ( p0->flag < mesh->base-1 && p1->flag < mesh->base-1 )  continue;

          l = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])\
            + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
          l = sqrt(l);

          if ( met->m[ip0] < met->m[ip1] ) {
            if ( met->m[ip0] < _MMG5_EPSD )  continue;
            hn = met->m[ip0] + mesh->info.hgrad*l;
            if ( met->m[ip1] > hn ) {
              met->m[ip1] = hn;
              p1->flag = mesh->base;
              nu++;
            }
          }
          else {
            if ( met->m[ip1] < _MMG5_EPSD )  continue;
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
    fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
  return(1);
}
