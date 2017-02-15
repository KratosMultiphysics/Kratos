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
 * \file common/mettools.c
 * \brief Metric tools for the mmg applications.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param t tangent at the ridge point.
 * \param n normal at the ridge point.
 * \param dtan metric size along the tangent direction.
 * \param dv metric size along the \f$t^{}n\f$ direction.
 * \param dn metric size along the normal direction.
 * \param m computed metric at the ridge point.
 * \return 1
 *
 * Build metric tensor at a fictitious ridge point, whose normal and tangent are
 * provided.
 *
 */
inline int
_MMG5_buildridmetfic(MMG5_pMesh mesh,double t[3],double n[3],double dtan,
                     double dv,double dn,double m[6]) {
  double u[3],r[3][3];

  u[0] = n[1]*t[2] - n[2]*t[1];
  u[1] = n[2]*t[0] - n[0]*t[2];
  u[2] = n[0]*t[1] - n[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(dtan,dv,dn)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n[2];

  m[0] = dtan*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1] + dn*r[0][2]*r[0][2];
  m[1] = dtan*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1] + dn*r[0][2]*r[1][2];
  m[2] = dtan*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1] + dn*r[0][2]*r[2][2];
  m[3] = dtan*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1] + dn*r[1][2]*r[1][2];
  m[4] = dtan*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1] + dn*r[1][2]*r[2][2];
  m[5] = dtan*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1] + dn*r[2][2]*r[2][2];

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m pointer toward the first metric to intersect.
 * \param n pointer toward the second metric to intersect.
 * \param mr pointer toward the computed intersected metric.
 * \return 1.
 *
 * Compute the intersected (2 x 2) metric between metrics \a m and \a n,
 * PRESERVING the directions of \a m. Result is stored in \a mr.
 *
 */
int _MMG5_intmetsavedir(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  int    i;
  double lambda[2],vp[2][2],siz,isqhmin;

  isqhmin = 1.0 / (mesh->info.hmin * mesh->info.hmin);
  _MMG5_eigensym(m,lambda,vp);

  for (i=0; i<2; i++) {
    siz = n[0]*vp[i][0]*vp[i][0] + 2.0*n[1]*vp[i][0]*vp[i][1]
      + n[2]*vp[i][1]*vp[i][1];
    lambda[i] = MG_MAX(lambda[i],siz);
    lambda[i] = MG_MIN(lambda[i],isqhmin);
  }
  mr[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
  mr[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
  mr[2] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param ux distance \f$[p0;p1]\f$ along x axis.
 * \param uy distance \f$[p0;p1]\f$ along y axis.
 * \param uz distance \f$[p0;p1]\f$ along z axis.
 * \param mr computed metric tensor.
 *
 * \return 1 if success
 *
 * Build metric tensor at ridge point p0, when computations with respect to p1
 * are to be held.
 *
 */
int _MMG5_buildridmet(MMG5_pMesh mesh,MMG5_pSol met,int np0,
                      double ux,double uy,double uz,double mr[6]) {
  MMG5_pPoint  p0;
  MMG5_pxPoint go;
  double       ps1,ps2,*n1,*n2,*t,*m,dv,dn,u[3],r[3][3];

  p0 = &mesh->point[np0];
  if ( !(MG_GEO & p0->tag) )  return(0);
  m = &met->m[6*np0];
  t = &p0->n[0];
  go = &mesh->xpoint[p0->xp];

  /* Decide between the two possible configurations */
  n1 = &go->n1[0];
  n2 = &go->n2[0];

  ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
  ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

  if ( fabs(ps2)<fabs(ps1) ) {
    n1 = &go->n2[0];
    dv = m[2];
    dn = m[4];
  }
  else{
    dv = m[1];
    dn = m[3];
  }

  u[0] = n1[1]*t[2] - n1[2]*t[1];
  u[1] = n1[2]*t[0] - n1[0]*t[2];
  u[2] = n1[0]*t[1] - n1[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,dn)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

  mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1] + dn*r[0][2]*r[0][2];
  mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1] + dn*r[0][2]*r[1][2];
  mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1] + dn*r[0][2]*r[2][2];
  mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1] + dn*r[1][2]*r[1][2];
  mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1] + dn*r[1][2]*r[2][2];
  mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1] + dn*r[2][2]*r[2][2];
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param nt normal direction at the ridge point.
 * \param mr computed metric tensor.
 *
 * Build metric tensor at ridge point \a p0, when the 'good' normal direction is
 * given by \a nt.
 *
 */
int _MMG5_buildridmetnor(MMG5_pMesh mesh,MMG5_pSol met,int np0,double nt[3],double mr[6]) {
  MMG5_pPoint  p0;
  MMG5_pxPoint go;
  double       ps1,ps2,*n1,*n2,*t,*m,dv,dn,u[3],r[3][3];

  p0 = &mesh->point[np0];
  if ( !(MG_GEO & p0->tag) )  return(0);
  m = &met->m[6*np0];
  t = &p0->n[0];
  go = &mesh->xpoint[p0->xp];

  /* Decide between the two possible configurations */
  n1 = &go->n1[0];
  n2 = &go->n2[0];

  ps1 = nt[0]*n1[0] + nt[1]*n1[1] + nt[2]*n1[2];
  ps2 = nt[0]*n2[0] + nt[1]*n2[1] + nt[2]*n2[2];

  if ( fabs(ps2) > fabs(ps1) ) {
    n1 = &go->n2[0];
    dv = m[2];
    dn = m[4];
  }
  else{
    dv = m[1];
    dn = m[3];
  }

  u[0] = n1[1]*t[2] - n1[2]*t[1];
  u[1] = n1[2]*t[0] - n1[0]*t[2];
  u[2] = n1[0]*t[1] - n1[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,0)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

  mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1] + dn*r[0][2]*r[0][2];
  mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1] + dn*r[0][2]*r[1][2];
  mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1] + dn*r[0][2]*r[2][2];
  mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1] + dn*r[1][2]*r[1][2];
  mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1] + dn*r[1][2]*r[2][2];
  mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1] + dn*r[2][2]*r[2][2];

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m pointer toward a \f$(2x2)\f$ metric.
 * \param n pointer toward a \f$(2x2)\f$ metric.
 * \param mr computed \f$(2x2)\f$ metric.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute the intersected (2 x 2) metric from metrics \a m and \a n : take
 * simultaneous reduction, and proceed to truncation in sizes.
 *
 */
int _MMG5_intersecmet22(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  double  det,imn[4],dd,sqDelta,trimn,lambda[2],vp0[2],vp1[2],dm[2],dn[2],vnorm,d0,d1,ip[4];
  double  isqhmin,isqhmax;

  isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < _MMG5_EPS*_MMG5_EPS ) {
    fprintf(stderr,"  ## Function intersecmet : null metric det : %E \n",det);
    return(0);
  }
  det = 1.0 / det;

  imn[0] = det * ( m[2]*n[0] - m[1]*n[1]);
  imn[1] = det * ( m[2]*n[1] - m[1]*n[2]);
  imn[2] = det * (-m[1]*n[0] + m[0]*n[1]);
  imn[3] = det * (-m[1]*n[1] + m[0]*n[2]);
  dd = imn[0] - imn[3];
  sqDelta = sqrt(fabs(dd*dd + 4.0*imn[1]*imn[2]));
  trimn = imn[0] + imn[3];

  lambda[0] = 0.5 * (trimn - sqDelta);
  if ( lambda[0] < 0.0 ) {
    fprintf(stderr," ## Eigenvalues : %f \n",lambda[0]);
    return(0);
  }

  /* First case : matrices m and n are homothetic : n = lambda0*m */
  if ( sqDelta < _MMG5_EPS ) {
    /* Diagonalize m and truncate eigenvalues : trimn, det, etc... are reused */
    if (fabs(m[1]) < _MMG5_EPS) {
      dm[0]   = m[0];
      dm[1]   = m[2];
      vp0[0] = 1;
      vp0[1] = 0;
      vp1[0] = 0;
      vp1[1] = 1;
    } else {
      dd    = m[0] - m[2];
      trimn = m[0] + m[2];
      det   = m[0]*m[2] - m[1]*m[1];

      sqDelta = sqrt(fabs(dd*dd +4*0*m[1]*m[1]));
      dm[0]   = 0.5 * (trimn + sqDelta);
      dm[1]   = 0.5 * (trimn - sqDelta);

      vp0[0] = m[1];
      vp0[1] = (dm[0]-m[0]);
      vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
      if ( vnorm < _MMG5_EPS ) {
        vp0[0] = (dm[0] - m[2]);
        vp0[1] = m[1];
        vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);

        if ( vnorm < _MMG5_EPS ) return(0);
      }

      vnorm   = 1.0 / vnorm;
      vp0[0] *= vnorm;
      vp0[1] *= vnorm;

      vp1[0] = m[1];
      vp1[1] = (dm[1]-m[0]);
      vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);

      if ( vnorm < _MMG5_EPS ) {
        vp1[0] = (dm[1] - m[2]);
        vp1[1] = m[1];
        vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);

        if ( vnorm < _MMG5_EPS ) return(0);
      }

      vnorm   = 1.0 / vnorm;
      vp1[0] *= vnorm;
      vp1[1] *= vnorm;
    }
    /* Eigenvalues of the resulting matrix*/
    dn[0] = MG_MAX(dm[0],lambda[0]*dm[0]);
    dn[0] = MG_MIN(isqhmin,MG_MAX(isqhmax,dn[0]));
    dn[1] = MG_MAX(dm[1],lambda[0]*dm[1]);
    dn[1] = MG_MIN(isqhmin,MG_MAX(isqhmax,dn[1]));

    /* Intersected metric = P diag(d0,d1){^t}P, P = (vp0, vp1) stored in columns */
    mr[0] = dn[0]*vp0[0]*vp0[0] + dn[1]*vp1[0]*vp1[0];
    mr[1] = dn[0]*vp0[0]*vp0[1] + dn[1]*vp1[0]*vp1[1];
    mr[2] = dn[0]*vp0[1]*vp0[1] + dn[1]*vp1[1]*vp1[1];

    return(1);
  }

  /* Second case : both eigenvalues of imn are distinct ; theory says qf associated to m and n
     are diagonalizable in basis (vp0, vp1) - the coreduction basis */
  else {
    lambda[1] = 0.5 * (trimn + sqDelta);
    assert(lambda[1] >= 0.0);

    vp0[0] = imn[1];
    vp0[1] = (lambda[0] - imn[0]);
    vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);

    if ( vnorm < _MMG5_EPS ) {
      vp0[0] = (lambda[0] - imn[3]);
      vp0[1] = imn[2];
      vnorm  = sqrt(vp0[0]*vp0[0] + vp0[1]*vp0[1]);
    }

    vnorm   = 1.0 / vnorm;
    vp0[0] *= vnorm;
    vp0[1] *= vnorm;

    vp1[0] = imn[1];
    vp1[1] = (lambda[1] - imn[0]);
    vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);

    if ( vnorm < _MMG5_EPS ) {
      vp1[0] = (lambda[1] - imn[3]);
      vp1[1] = imn[2];
      vnorm  = sqrt(vp1[0]*vp1[0] + vp1[1]*vp1[1]);
    }

    vnorm   = 1.0 / vnorm;
    vp1[0] *= vnorm;
    vp1[1] *= vnorm;

    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp0[0]*vp0[0] + 2.0*m[1]*vp0[0]*vp0[1] + m[2]*vp0[1]*vp0[1];
    dm[1] = m[0]*vp1[0]*vp1[0] + 2.0*m[1]*vp1[0]*vp1[1] + m[2]*vp1[1]*vp1[1];
    dn[0] = n[0]*vp0[0]*vp0[0] + 2.0*n[1]*vp0[0]*vp0[1] + n[2]*vp0[1]*vp0[1];
    dn[1] = n[0]*vp1[0]*vp1[0] + 2.0*n[1]*vp1[0]*vp1[1] + n[2]*vp1[1]*vp1[1];

    /* Diagonal values of the intersected metric */
    d0 = MG_MAX(dm[0],dn[0]);
    d0 = MG_MIN(isqhmin,MG_MAX(d0,isqhmax));

    d1 = MG_MAX(dm[1],dn[1]);
    d1 = MG_MIN(isqhmin,MG_MAX(d1,isqhmax));

    /* Intersected metric = tP^-1 diag(d0,d1)P^-1, P = (vp0, vp1) stored in columns */
    det = vp0[0]*vp1[1] - vp0[1]*vp1[0];
    if ( fabs(det) < _MMG5_EPS )  return(0);
    det = 1.0 / det;

    ip[0] =  vp1[1]*det;
    ip[1] = -vp1[0]*det;
    ip[2] = -vp0[1]*det;
    ip[3] =  vp0[0]*det;

    mr[0] = d0*ip[0]*ip[0] + d1*ip[2]*ip[2];
    mr[1] = d0*ip[0]*ip[1] + d1*ip[2]*ip[3];
    mr[2] = d0*ip[1]*ip[1] + d1*ip[3]*ip[3];
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param np global index of vertex in which we intersect the metrics.
 * \param me physical metric at point \a np.
 * \param n normal or tangent at point np.
 * \return 0 if fail, 1 otherwise.
 *
 * Intersect the surface metric held in np (supported in tangent plane of \a np)
 * with 3*3 physical metric in \a me. For ridge points, this function fill the
 * \f$ p_0->m[3]\f$ and \f$ p_0->m[4]\f$ fields that contains respectively the
 * specific sizes in the \f$n_1\f$ and \f$n_2\f$ directions.
 *
 */
int _MMG5_mmgIntextmet(MMG5_pMesh mesh,MMG5_pSol met,int np,double me[6],
                       double n[3]) {
  MMG5_pPoint         p0;
  MMG5_pxPoint        go;
  double              hu,isqhmin,isqhmax,dd,alpha1,alpha2,alpha3,u[3];
  double              lambda[3],vp[3][3];
  double              *m,*n1,*n2,*t,r[3][3],mrot[6],mr[3],mtan[3],metan[3];
  char                i;

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  p0 = &mesh->point[np];
  m  = &met->m[6*np];

  /* Case of a singular point : take smallest size prescribed by met, or me in
   * every direction */
  if ( MG_SIN(p0->tag) || (p0->tag & MG_NOM) ) {
    /* Characteristic polynomial of me */
    _MMG5_eigenv(1,me,lambda,vp);

    hu = m[0];
    for(i=0; i<3; i++) {
        hu = MG_MAX(hu,lambda[i]);
    }
    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[0] = hu;
    m[3] = hu;
    m[5] = hu;
  }
  /* Case of a ridge point : take sizes in 3 directions t,n1,u */
  else if ( p0->tag & MG_GEO ) {
    /* Size prescribed by metric me in direction t */
    t = n;
    hu = me[0]*t[0]*t[0] + me[3]*t[1]*t[1] + me[5]*t[2]*t[2] \
      + 2.0*me[1]*t[0]*t[1] + 2.0*me[2]*t[0]*t[2] + 2.0*me[4]*t[1]*t[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[0] = MG_MAX(m[0],hu);

    /* Size prescribed by metric me in direction u1 = n1 ^ t */
    assert ( p0->xp );
    go = &mesh->xpoint[p0->xp];
    n1 = &go->n1[0];
    n2 = &go->n2[0];

    u[0] = n1[1]*t[2] - n1[2]*t[1];
    u[1] = n1[2]*t[0] - n1[0]*t[2];
    u[2] = n1[0]*t[1] - n1[1]*t[0];
    dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if ( dd > _MMG5_EPSD ) {
      dd = 1.0 / sqrt(dd);

      u[0] *= dd;
      u[1] *= dd;
      u[2] *= dd;

      hu = me[0]*u[0]*u[0] + me[3]*u[1]*u[1] + me[5]*u[2]*u[2]          \
        + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

      hu = MG_MIN(isqhmin,hu);
      hu = MG_MAX(isqhmax,hu);
      m[1] = MG_MAX(m[1],hu);
    }
    /* Size prescribed by metric me in direction u2 = n2 ^ t */
    u[0] = n2[1]*t[2] - n2[2]*t[1];
    u[1] = n2[2]*t[0] - n2[0]*t[2];
    u[2] = n2[0]*t[1] - n2[1]*t[0];
    dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if ( dd > _MMG5_EPSD ) {
      dd = 1.0 / sqrt(dd);

      u[0] *= dd;
      u[1] *= dd;
      u[2] *= dd;

      hu =     me[0]*u[0]*u[0] +     me[3]*u[1]*u[1] +     me[5]*u[2]*u[2] \
        + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

      hu = MG_MIN(isqhmin,hu);
      hu = MG_MAX(isqhmax,hu);
      m[2] = MG_MAX(m[2],hu);
    }

    /* Size prescribed by metric me in direction n1 */
    hu = me[0]*n1[0]*n1[0] + me[3]*n1[1]*n1[1] + me[5]*n1[2]*n1[2]
      + 2.0*me[1]*n1[0]*n1[1] + 2.0*me[2]*n1[0]*n1[2] + 2.0*me[4]*n1[1]*n1[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[3] = hu;

    /* Size prescribed by metric me in direction n2 */
    hu = me[0]*n2[0]*n2[0] + me[3]*n2[1]*n2[1] + me[5]*n2[2]*n2[2]
      + 2.0*me[1]*n2[0]*n2[1] + 2.0*me[2]*n2[0]*n2[2] + 2.0*me[4]*n2[1]*n2[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[4] = hu;
  }
  /* Case of a ref, or regular point : intersect metrics in tangent plane */
  else {
    _MMG5_rotmatrix(n,r);

    /* Expression of both metrics in tangent plane */
    _MMG5_rmtr(r,m,mrot);
    mtan[0] = mrot[0];
    mtan[1] = mrot[1];
    mtan[2] = mrot[3];


    _MMG5_rmtr(r,me,mrot);
    metan[0] = mrot[0];
    metan[1] = mrot[1];
    metan[2] = mrot[3];

    /* Intersection of metrics in the tangent plane */
    if ( !_MMG5_intersecmet22(mesh,mtan,metan,mr) ) {
      fprintf(stderr,"WARNING IMPOSSIBLE INTERSECTION : SURFACIC METRIC SKIPPED \n");
      m[0] = me[0];
      m[1] = me[1];
      m[2] = me[2];
      m[3] = me[3];
      m[4] = me[4];
      m[5] = me[5];

      return(0);
    }

    /* Back to the canonical basis of \mathbb{R}^3 : me = ^tR*mr*R : mtan and
     * metan are reused */
    mtan[0]  = mr[0]*r[0][0] + mr[1]*r[1][0];
    mtan[1]  = mr[0]*r[0][1] + mr[1]*r[1][1];
    mtan[2]  = mr[0]*r[0][2] + mr[1]*r[1][2];
    metan[0] = mr[1]*r[0][0] + mr[2]*r[1][0];
    metan[1] = mr[1]*r[0][1] + mr[2]*r[1][1];
    metan[2] = mr[1]*r[0][2] + mr[2]*r[1][2];
 
    alpha1 = r[2][0]*mrot[5];
    alpha2 = r[2][1]*mrot[5];
    alpha3 = r[2][2]*mrot[5];

    m[0] = r[0][0] * mtan[0] + r[1][0] * metan[0] + r[2][0]*alpha1;
    m[1] = r[0][0] * mtan[1] + r[1][0] * metan[1] + r[2][0]*alpha2;
    m[2] = r[0][0] * mtan[2] + r[1][0] * metan[2] + r[2][0]*alpha3;
    m[3] = r[0][1] * mtan[1] + r[1][1] * metan[1] + r[2][1]*alpha2;
    m[4] = r[0][1] * mtan[2] + r[1][1] * metan[2] + r[2][1]*alpha3;
    m[5] = r[0][2] * mtan[2] + r[1][2] * metan[2] + r[2][2]*alpha3;

    /* Truncate the metric in the third direction (because me was not
     * truncated) */
    _MMG5_eigenv(1,m,lambda,vp);

    for (i=0; i<3; i++) {
      if(lambda[i]<=0) {
        fprintf(stderr,"%s:%d:Error: wrong metric at point %d -- eigenvalues :"
               " %e %e %e\n",__FILE__,__LINE__,
               np,lambda[0],lambda[1],lambda[2]);
        fprintf(stderr,"  ## Surfacic metric skipped. \n");
        m[0] = me[0];
        m[1] = me[1];
        m[2] = me[2];
        m[3] = me[3];
        m[4] = me[4];
        m[5] = me[5];
        return(0);
      }
      lambda[i]=MG_MIN(isqhmin,lambda[i]);
      lambda[i]=MG_MAX(isqhmax,lambda[i]);
    }

    m[0] = vp[0][0]*vp[0][0]*lambda[0] + vp[1][0]*vp[1][0]*lambda[1]
      + vp[2][0]*vp[2][0]*lambda[2];
    m[1] = vp[0][0]*vp[0][1]*lambda[0] + vp[1][0]*vp[1][1]*lambda[1]
      + vp[2][0]*vp[2][1]*lambda[2];
    m[2] = vp[0][0]*vp[0][2]*lambda[0] + vp[1][0]*vp[1][2]*lambda[1]
      + vp[2][0]*vp[2][2]*lambda[2];
    m[3] = vp[0][1]*vp[0][1]*lambda[0] + vp[1][1]*vp[1][1]*lambda[1]
      + vp[2][1]*vp[2][1]*lambda[2];
    m[4] = vp[0][1]*vp[0][2]*lambda[0] + vp[1][1]*vp[1][2]*lambda[1]
      + vp[2][1]*vp[2][2]*lambda[2];
    m[5] = vp[0][2]*vp[0][2]*lambda[0] + vp[1][2]*vp[1][2]*lambda[1]
      + vp[2][2]*vp[2][2]*lambda[2];
  }

  return(1);
}

/**
 * \param c0 table of the coordinates of the starting point.
 * \param n0 normal at the starting point.
 * \param m metric to be transported.
 * \param c1 table of the coordinates of the ending point.
 * \param n1 normal at the ending point.
 * \param mt computed metric.
 * \return 0 if fail, 1 otherwise.
 *
 * Parallel transport of a metric tensor field, attached to point \a c0, with
 * normal \a n0, to point \a c1, with normal \a n1.
 *
 */
int _MMG5_paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]) {
  double  r[3][3],mrot[6],mtan[3],lambda[2],vp[2][2],u[3],ps,ll;

  /* Take the induced metric tensor in the tangent plane by change of basis : R * M * {^t}R*/
  if ( !_MMG5_rotmatrix(n0,r) )  return(0);
  _MMG5_rmtr(r,m,mrot);
  mtan[0] = mrot[0];
  mtan[1] = mrot[1];
  mtan[2] = mrot[3];

  /* Take eigenvectors of metric tensor in tangent plane */
  _MMG5_eigensym(mtan,lambda,vp);

  /* Eigenvector in canonical basis = {t}R*vp[0] */
  u[0] = r[0][0]*vp[0][0] + r[1][0]*vp[0][1];
  u[1] = r[0][1]*vp[0][0] + r[1][1]*vp[0][1];
  u[2] = r[0][2]*vp[0][0] + r[1][2]*vp[0][1];

  /* Projection in the tangent plane of c1 */
  ps = u[0]*n1[0] + u[1]*n1[1] + u[2]*n1[2];
  u[0] -= ps*n1[0];
  u[1] -= ps*n1[1];
  u[2] -= ps*n1[2];
  ll = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
  if ( ll < _MMG5_EPSD )  return(0);
  ll = 1.0 / sqrt(ll);
  u[0] *= ll;
  u[1] *= ll;
  u[2] *= ll;

  /* And the transported metric is diag(lambda[0], lambda[1], mrot[5]) in basis
   * (u,n1^u,n1) */
  r[0][0] = u[0];
  r[1][0] = u[1];
  r[2][0] = u[2];

  r[0][1] = n1[1]*u[2] - n1[2]*u[1];
  r[1][1] = n1[2]*u[0] - n1[0]*u[2];
  r[2][1] = n1[0]*u[1] - n1[1]*u[0];

  ll = r[0][1]*r[0][1] + r[1][1]*r[1][1] + r[2][1]*r[2][1];
  if ( ll < _MMG5_EPSD )  return(0);
  ll = 1.0 / sqrt(ll);
  r[0][1] *= ll;
  r[1][1] *= ll;
  r[2][1] *= ll;

  r[0][2] = n1[0];
  r[1][2] = n1[1];
  r[2][2] = n1[2];

  /*mt = R * diag(lambda[0], lambda[1], mrot[5])*{^t}R */
  mt[0] = lambda[0]*r[0][0]*r[0][0] + lambda[1]*r[0][1]*r[0][1]
    + mrot[5]*r[0][2]*r[0][2];

  mt[1] = lambda[0]*r[0][0]*r[1][0]
    + lambda[1]*r[0][1]*r[1][1] + mrot[5]*r[0][2]*r[1][2];

  mt[2] = lambda[0]*r[0][0]*r[2][0]
    + lambda[1]*r[0][1]*r[2][1] + mrot[5]*r[0][2]*r[2][2];

  mt[3] = lambda[0]*r[1][0]*r[1][0] + lambda[1]*r[1][1]*r[1][1]
    + mrot[5]*r[1][2]*r[1][2];

  mt[4] = lambda[0]*r[2][0]*r[1][0]
    + lambda[1]*r[2][1]*r[1][1] + mrot[5]*r[2][2]*r[1][2];

  mt[5] = lambda[0]*r[2][0]*r[2][0] + lambda[1]*r[2][1]*r[2][1]
    + mrot[5]*r[2][2]*r[2][2];

  return(1);
}
