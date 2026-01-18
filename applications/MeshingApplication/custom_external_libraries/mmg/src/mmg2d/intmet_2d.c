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
 * \file mmg2d/intmet_2d.c
 * \brief Interpolation of metrics
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "libmmg2d_private.h"

/* Interpolation of isotropic metric met along edge i of triangle k, according to parameter s;
   ip = index of the new point */
int MMG2D_intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,double s) {
  MMG5_pTria  pt;
  MMG5_int    ip1,ip2;
  int8_t      i1,i2;
  
  pt  = &mesh->tria[k];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  met->m[ip] = (1.0-s)*met->m[ip1] + s*met->m[ip2];
  
  return 1;
}

/* Interpolation of metrics m and n according to parameter s; result is stored in mr. A simultaneous reduction of both matrices is performed and the sizes are interpolated. */
int MMG5_interpmet22(MMG5_pMesh mesh,double *m,double *n,double s,double *mr) {
  double        det,imn[4],dd,den,sqDelta,trimn,lambda[2],vp[2][2],dm[2];
  double        dn[2],vnorm,d0,d1,ip[4];
  static int8_t mmgWarn0=0,mmgWarn1=0;

  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < MMG5_EPS*MMG5_EPS ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Error: %s: null metric det : %E \n",
              __func__,det);
    }
    return 0;
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
    if ( !mmgWarn1 ) {
      mmgWarn1 = 1;
      fprintf(stderr,"\n  ## Error: %s: at least 1 negative eigenvalue: %f \n",
              __func__,lambda[0]);
    }
    return 0;
  }

  /* First case : matrices m and n are homothetic: n = lambda0*m */
  if ( sqDelta < MMG5_EPS ) {
    /* Diagonalize m and truncate eigenvalues : trimn, det, etc... are reused */
    if ( fabs(m[1]) <= MMG5_EPS || fabs(n[1]) <= MMG5_EPS ) {
      dm[0]   = m[0];
      dm[1]   = m[2];
      vp[0][0] = 1;
      vp[0][1] = 0;
      vp[1][0] = 0;
      vp[1][1] = 1;
    }
    else {
      MMG5_eigensym( m, dm, vp );

      vp[0][0] = vp[0][0];
      vp[0][1] = vp[0][1];

      vp[1][0] = vp[1][0];
      vp[1][1] = vp[1][1];
    }
    
    /* Eigenvalues of metric n */
    dn[0] = lambda[0]*dm[0];
    dn[1] = lambda[0]*dm[1];
    
    /* Diagonal values of the intersected metric */
    dd = dn[0]*dm[0];
    den = (1.0-s)*(1.0-s)*dn[0] + s*s*dm[0] + 2.0*s*(1.0-s)*sqrt(dd);
    
    /* If den is too small (should not happen) simply interpolate diagonal values; else interpolate sizes */
    if ( den < MMG5_EPS ) d0 = (1.0-s)*dm[0] + s*dn[0];
    else d0 = dd / den;
    
    dd = dn[1]*dm[1];
    den = (1.0-s)*(1.0-s)*dn[1] + s*s*dm[1] + 2.0*s*(1.0-s)*sqrt(dd);

    if ( den < MMG5_EPS ) d1 = (1.0-s)*dm[1] + s*dn[1];
    else d1 = dd / den;

    /* Intersected metric = P diag(d0,d1){^t}P, P = (vp[0], vp[1]) stored in columns */
    mr[0] = d0*vp[0][0]*vp[0][0] + d1*vp[1][0]*vp[1][0];
    mr[1] = d0*vp[0][0]*vp[0][1] + d1*vp[1][0]*vp[1][1];
    mr[2] = d0*vp[0][1]*vp[0][1] + d1*vp[1][1]*vp[1][1];
    
    return 1;
  }
  
  /* Second case: both eigenvalues of imn are distinct ; theory says qf associated to m and n
   are diagonalizable in basis (vp[0], vp[1]) - the coreduction basis */
  else {
    lambda[1] = 0.5 * (trimn + sqDelta);
    assert(lambda[1] >= 0.0);
    
    vp[0][0] = imn[1];
    vp[0][1] = (lambda[0] - imn[0]);
    vnorm  = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);
    
    if ( vnorm < MMG5_EPS ) {
      vp[0][0] = (lambda[0] - imn[3]);
      vp[0][1] = imn[2];
      vnorm  = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);
    }
    
    vnorm   = 1.0 / vnorm;
    vp[0][0] *= vnorm;
    vp[0][1] *= vnorm;
    
    vp[1][0] = imn[1];
    vp[1][1] = (lambda[1] - imn[0]);
    vnorm  = sqrt(vp[1][0]*vp[1][0] + vp[1][1]*vp[1][1]);
    
    if ( vnorm < MMG5_EPS ) {
      vp[1][0] = (lambda[1] - imn[3]);
      vp[1][1] = imn[2];
      vnorm  = sqrt(vp[1][0]*vp[1][0] + vp[1][1]*vp[1][1]);
    }
    
    vnorm   = 1.0 / vnorm;
    vp[1][0] *= vnorm;
    vp[1][1] *= vnorm;

    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp[0][0]*vp[0][0] + 2.0*m[1]*vp[0][0]*vp[0][1] + m[2]*vp[0][1]*vp[0][1];
    dm[1] = m[0]*vp[1][0]*vp[1][0] + 2.0*m[1]*vp[1][0]*vp[1][1] + m[2]*vp[1][1]*vp[1][1];
    dn[0] = n[0]*vp[0][0]*vp[0][0] + 2.0*n[1]*vp[0][0]*vp[0][1] + n[2]*vp[0][1]*vp[0][1];
    dn[1] = n[0]*vp[1][0]*vp[1][0] + 2.0*n[1]*vp[1][0]*vp[1][1] + n[2]*vp[1][1]*vp[1][1];
    
    /* Diagonal values of the intersected metric */
    dd = dn[0]*dm[0];
    den = (1.0-s)*(1.0-s)*dn[0] + s*s*dm[0] + 2.0*s*(1.0-s)*sqrt(dd);
    
    if ( den < MMG5_EPS ) d0 = (1.0-s)*dm[0] + s*dn[0];
    else d0 = dd / den;
    
    dd = dn[1]*dm[1];
    den = (1.0-s)*(1.0-s)*dn[1] + s*s*dm[1] + 2.0*s*(1.0-s)*sqrt(dd);
    
    if ( den < MMG5_EPS ) d1 = (1.0-s)*dm[1] + s*dn[1];
    else d1 = dd / den;
    
    /* Intersected metric = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored in columns */
    det = vp[0][0]*vp[1][1] - vp[0][1]*vp[1][0];
    if ( fabs(det) < MMG5_EPS )  return 0;
    det = 1.0 / det;
    
    ip[0] =  vp[1][1]*det;
    ip[1] = -vp[1][0]*det;
    ip[2] = -vp[0][1]*det;
    ip[3] =  vp[0][0]*det;
    
    mr[0] = d0*ip[0]*ip[0] + d1*ip[2]*ip[2];
    mr[1] = d0*ip[0]*ip[1] + d1*ip[2]*ip[3];
    mr[2] = d0*ip[1]*ip[1] + d1*ip[3]*ip[3];
  }
  
  return 1;
}

/* Interpolation of anisotropic metric met along edge i of triangle k, according to parameter s;
 ip = index of the new point */
int MMG2D_intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,double s) {
  MMG5_pTria     pt;
  double         *m1,*m2,*mr;
#ifndef NDEBUG
  double         det;
#endif
  MMG5_int       ip1,ip2;
  int8_t         i1,i2;
  static int8_t  mmgWarn=0;

  pt = &mesh->tria[k];
  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  m1 = &met->m[3*ip1];
  m2 = &met->m[3*ip2];
  mr = &met->m[3*ip];

  if ( !MMG5_interpmet22(mesh,m1,m2,s,mr) ) {
    if ( !mmgWarn ) {
      mmgWarn=1;
      fprintf(stderr,"  ## Error: %s: at least 1 naive interpolation.\n",
              __func__);
    }
    mr[0] = (1.0-s)*m1[0] + s*m2[0];
    mr[1] = (1.0-s)*m1[1] + s*m2[1];
    mr[2] = (1.0-s)*m1[2] + s*m2[2];
  }

#ifndef NDEBUG
  // Check the result
  det = mr[0]*mr[2] - mr[1]*mr[1];
  if ( fabs(det) < MMG5_EPS*MMG5_EPS ) {
    fprintf(stderr,"\n  ## Error: %s: interpolation results in null metric det:"
            " %e \n",__func__,det);

    printf("\nInterpolated metrics:\n");

    double lambda[2],vp[2][2];
    printf("Metric %" MMG5_PRId " (tag %d): %e %e %e\n",ip1,mesh->point[ip1].tag,m1[0],m1[1],m1[2]);

    MMG5_eigen2(m1,lambda,vp);
    printf ("eigenval: %e %e\n",lambda[0],lambda[1] );
    printf ("size:     %e %e\n",1./sqrt(lambda[0]),1./sqrt(lambda[1]) );
    printf ("eigenvec: %e %e\n",vp[0][0],vp[0][1] );
    printf ("eigenvec: %e %e\n\n",vp[1][0],vp[1][1] );

    printf("Metric %" MMG5_PRId " (tag %d): %e %e %e\n",ip2,mesh->point[ip2].tag,m2[0],m2[1],m2[2]);
    MMG5_eigen2(m2,lambda,vp);
    printf ("eigenval: %e %e\n",lambda[0],lambda[1] );
    printf ("size:     %e %e\n",1./sqrt(lambda[0]),1./sqrt(lambda[1]) );
    printf ("eigenvec: %e %e\n",vp[0][0],vp[0][1] );
    printf ("eigenvec: %e %e\n\n",vp[1][0],vp[1][1] );

    assert(0);
  }
#endif

  return 1;
}
