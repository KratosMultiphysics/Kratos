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
#include "mmg2dexterns_private.h"

/* Compute isotropic edge length */
double long_iso(double *ca,double *cb,double *ma,double *mb) {
  double   ha,hb,ux,uy,dd,rap,len;

  ha = *ma;
  hb = *mb;
  ux = cb[0] - ca[0];
  uy = cb[1] - ca[1];
  dd = sqrt(ux*ux + uy*uy);

  rap = (hb - ha) / ha;
  if ( fabs(rap) < MMG2D_EPSD )
    len = dd / ha;
  else
    len = dd * (1.0/ha + 1.0/hb + 8.0 / (ha+hb)) / 6.0;

  return len;
}


/* compute aniso edge length */
double long_ani(double *ca,double *cb,double *ma,double *mb) {
  double   ux,uy,dd1,dd2,len;
  ux = cb[0] - ca[0];
  uy = cb[1] - ca[1];

  dd1 = ma[0]*ux*ux + ma[2]*uy*uy + 2.0*ma[1]*ux*uy;
  if ( dd1 <= 0.0 )  dd1 = 0.0;
  dd2 = mb[0]*ux*ux + mb[2]*uy*uy + 2.0*mb[1]*ux*uy;
  if ( dd2 <= 0.0 )  dd2 = 0.0;

  len = (sqrt(dd1)+sqrt(dd2)+4.0*sqrt(0.5*(dd1+dd2))) / 6.0;

  return len;
}

/** Calculate length of a curve in the considered isotropic metric */
double MMG2D_lencurv_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int ip1,MMG5_int ip2) {
  MMG5_pPoint     p1,p2;
  double          h1,h2,len,l,r;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  h1 = met->m[ip1];
  h2 = met->m[ip2];

  l = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = ( fabs(r) < MMG5_EPS ) ? ( l/h1 ) : ( l / (h2-h1) * log1p(r) );

  return len;
}

/* Calculate length of a curve in the considered anisotropic metric by using a two-point quadrature formula */
double MMG2D_lencurv_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int ip1,MMG5_int ip2) {
  MMG5_pPoint      p1,p2;
  double           len,*m1,*m2,ux,uy,l1,l2;
  static int8_t    mmgWarn0=0,mmgWarn1=0;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  m1 = &met->m[3*ip1];
  m2 = &met->m[3*ip2];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];

  l1 = m1[0]*ux*ux + 2.0*m1[1]*ux*uy + m1[2]*uy*uy;
  l2 = m2[0]*ux*ux + 2.0*m2[1]*ux*uy + m2[2]*uy*uy;

  if ( l1 < 0.0 ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Error: %s: at least 1 negative edge length"
              " (l1: %e).\n",__func__,l1);
    }
    return 0.;
  }
  if ( l2 < 0.0 ) {
    if ( !mmgWarn1 ) {
      mmgWarn1 = 1;
      fprintf(stderr,"\n  ## Error: %s: at least 1 negative edge length"
              " (l2: %e)\n",__func__,l2);
    }
    return 0.;
  }

  l1 = sqrt(l1);
  l2 = sqrt(l2);

  len = 0.5*(l1+l2);

  return len;
}

/* print histo of edge lengths */
int MMG2D_prilen(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria    pt;
  double        lavg,len,lmin,lmax;
  int           ia,ipa,ipb;
  MMG5_int      iamin,ibmin,iamax,ibmax,hl[9],nullEdge,navg;
  static double bd[9] = {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};
  MMG5_int      k,l;

  navg  = 0;
  lavg  = 0.0;
  lmin  = 1.e20;
  lmax  = 0.0;
  iamin = 0;
  ibmin = 0;
  iamax = 0;
  ibmax = 0;
  nullEdge = 0;

  for (k=0; k<9; k++)  hl[k] = 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (ia=0; ia<3; ia++) {
      l = (&mesh->adja[3*(k-1)+1])[ia];
      if ( l < 3*k )  continue;

      ipa = MMG2D_iare[ia][0];
      ipb = MMG2D_iare[ia][1];

      if ( sol->m )
        len = MMG2D_lencurv(mesh,sol,pt->v[ipa],pt->v[ipb]);
      else
        len = MMG2D_lencurv_iso(mesh,sol,pt->v[ipa],pt->v[ipb]);

      navg++;
      lavg += len;

      /* find largest, smallest edge */
      if (len < lmin) {
        lmin  = len;
        iamin = pt->v[ipa];
        ibmin = pt->v[ipb];
      }
      if (len > lmax) {
        lmax  = len;
        iamax = pt->v[ipa];
        ibmax = pt->v[ipb];
      }

      /* update histogram */
      if (len < bd[3]) {
        if (len > bd[2])       hl[2]++;
        else if (len > bd[1])  hl[1]++;
        else                   hl[0]++;
      }
      else if (len < bd[5]) {
        if (len > bd[4])       hl[4]++;
        else if (len > bd[3])  hl[3]++;
      }
      else if (len < bd[6])    hl[5]++;
      else if (len < bd[7])    hl[6]++;
      else if (len < bd[8])    hl[7]++;
      else                     hl[8]++;
    }
  }
  MMG5_displayLengthHisto(mesh, navg, &lavg, iamin, ibmin, lmin,
                          iamax, ibmax, lmax,nullEdge, &bd[0], &hl[0],0);


  return 1;
}
