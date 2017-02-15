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
#include "mmg2d.h"


/* compute iso edge length */
double long_iso(double *ca,double *cb,double *ma,double *mb) {
  double   ha,hb,ux,uy,dd,rap,len;

  ha = *ma;
  hb = *mb;
  ux = cb[0] - ca[0];
  uy = cb[1] - ca[1];
  dd = sqrt(ux*ux + uy*uy);

  rap = (hb - ha) / ha;
  if ( fabs(rap) < EPSD )
    len = dd / ha;
  else
    len = dd * (1.0/ha + 1.0/hb + 8.0 / (ha+hb)) / 6.0;

  return(len);
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

  return(len);
}
/* print histo of edge lengths */
int MMG2_prilen(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria       pt;
  double      lavg,len,ecart,som,lmin,lmax,*ca,*cb,*ma,*mb;
  int         k,l,navg,ia,ipa,ipb,iamin,ibmin,iamax,ibmax,hl[9];
  int       iadr;
  static double bd[9] = {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};
//{0.0, 0.2, 0.5, 0.7071, 0.9, 1.111, 1.4142, 2.0, 5.0 };
  navg  = 0;
  lavg  = 0.0;
  lmin  = 1.e20;
  lmax  = 0.0;
  som   = 0.0;
  iamin = 0;
  ibmin = 0;
  iamax = 0;
  ibmax = 0;

  for (k=0; k<9; k++)  hl[k] = 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;

    for (ia=0; ia<3; ia++) {
      l = (&mesh->adja[3*(k-1)+1])[ia];
      if ( l < 3*k )  continue;

      ipa = MMG2_iare[ia][0];
      ipb = MMG2_iare[ia][1];
      ca  = &mesh->point[pt->v[ipa]].c[0];
      cb  = &mesh->point[pt->v[ipb]].c[0];

      iadr = pt->v[ipa]*sol->size;
      ma   = &sol->m[iadr];
      iadr = pt->v[ipb]*sol->size;
      mb   = &sol->m[iadr];

      len = MMG2_length(ca,cb,ma,mb);
      navg++;
      ecart = len;
      lavg += len;

      /* update efficiency index */
      if ( ecart > 1.0 )  ecart = 1.0 / ecart;

      som  += (ecart - 1.0);

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
  _MMG5_displayHisto(mesh, navg, &lavg, iamin, ibmin, lmin,
                     iamax, ibmax, lmax, &bd[0], &hl[0]);


  return(1);
}
