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
 * \file mmg2d/quality_2d.c
 * \brief Functions to compute the quality.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"

/**
 * \param mesh pointer to the mesh
 * \param pt pointer to the tria
 *
 * \return the oriented area of the triangle.
 *
 * Compute oriented area of tria pt
 *
 */
double MMG2D_quickcal(MMG5_pMesh mesh, MMG5_pTria pt) {
  MMG5_pPoint        p0,p1,p2;
  double             cal;

  p0 = &mesh->point[pt->v[0]];
  p1 = &mesh->point[pt->v[1]];
  p2 = &mesh->point[pt->v[2]];

  cal = MMG2D_quickarea(p0->c,p1->c,p2->c);
  return cal;
}

/**
 * \param a coordinates of first vertex of tria
 * \param b coordinates of second vertex of tria
 * \param c coordinates of third vertex of tria
 *
 * \return non-normalized quality if success, 0 if triangle is null or inverted.
 *
 * Compute quality of a triangle from the datum of its 3 vertices.
 *
 */
double MMG2D_caltri_iso_3pt(double *a,double *b,double *c) {
  double        abx,aby,acx,acy,bcx,bcy,area,h1,h2,h3,hm;

  abx = b[0] - a[0];
  aby = b[1] - a[1];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];

  /* orientation */
  area = abx*acy - aby*acx;
  if ( area <= 0.0 ) return 0.0;

  /* edge lengths */
  h1 = abx*abx + aby*aby;
  h2 = acx*acx + acy*acy;
  h3 = bcx*bcx + bcy*bcy;

  hm = h1 + h2 + h3;

  if ( hm > 0. ) {
    return  area / hm;
  }
  else {
    return 0.0;
  }
}

/**
 * \param pointer to the mesh
 * \param pointer to the metric (for compatibility with aniso interface)
 * \param pt pointer to the tria
 *
 * \return non-normalized quality if success, 0 if triangle is null or inverted.
 *
 * Compute quality of the triangle pt when the supplied metric is isotropic;
 * return 0 in the case that the triangle has inverted orientation.
 */
double MMG2D_caltri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt) {
  double    *a,*b,*c;

  a  = mesh->point[pt->v[0]].c;
  b  = mesh->point[pt->v[1]].c;
  c  = mesh->point[pt->v[2]].c;

  return MMG2D_caltri_iso_3pt(a,b,c);
}



/* Compute quality of the triangle pt when the supplied metric is anisotropic;
 return 0 in the case that the triangle has inverted orientation */
double MMG2D_caltri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt) {
  double     abx,aby,acx,acy,bcx,bcy;
  double     *a,*b,*c,*ma,*mb,*mc;
  double     area,aream,hm,m[6],h1,h2,h3;
  MMG5_int   ipa,ipb,ipc;
  int        i;

  ipa = pt->v[0];
  ipb = pt->v[1];
  ipc = pt->v[2];

  a  = mesh->point[ipa].c;
  b  = mesh->point[ipb].c;
  c  = mesh->point[ipc].c;

  ma = &met->m[3*ipa];
  mb = &met->m[3*ipb];
  mc = &met->m[3*ipc];

  /* Isotropic area of pt */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  area = abx*acy - aby*acx;
  if ( area <= 0.0 ) return 0.0;

  for (i=0; i<3; i++)  m[i] = (ma[i]+mb[i]+mc[i]) / 3.0;

  /*  Anisotropic edge lengths */
  h1 = m[0]*abx*abx + m[2]*aby*aby + 2.0*m[1]*abx*aby;
  h1 = h1 > 0.0 ? sqrt(h1) : 0.0;
  h2 = m[0]*acx*acx + m[2]*acy*acy + 2.0*m[1]*acx*acy;
  h2 = h2 > 0.0 ? sqrt(h2) : 0.0;
  h3 = m[0]*bcx*bcx + m[2]*bcy*bcy + 2.0*m[1]*bcx*bcy;
  h3 = h3 > 0.0 ? sqrt(h3) : 0.0;

  hm = h1*h1 + h2*h2 + h3*h3;

  /* Anisotropic volume of pt */
  aream = sqrt(m[0]*m[2]-m[1]*m[1])*area;

  /* Quality measure = (Vol_M(T) / (l(ab)^2+l(ac)^2+l(bc)^2)) */
  if ( hm > 0. ) {
    return  aream/hm;
  }
  else {
    return (0.0);
  }
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 *
 * \return 0 if the worst element has a nul quality, 1 otherwise.
 *
 * Print histogram of mesh qualities.
 *
 */
int MMG2D_outqua(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  double        rap,rapmin,rapmax,rapavg,med,good;
  int           i,ir,imax,his[5];
  MMG5_int      k,iel,ok,nex;
  static int8_t mmgWarn0;

  if ( !mesh->nt ) {
    return 1;
  }

  /* Compute triangle quality*/
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if( !MG_EOK(pt) )   continue;

    if ( !met->m ) {
      pt->qual = MMG2D_caltri_iso(mesh,met,pt);
    }
    else
      pt->qual = MMG2D_caltri(mesh,met,pt);
  }
  if ( mesh->info.imprim <= 0 ) return 1;

  rapmin  = 2.0;
  rapmax  = 0.0;
  rapavg  = med = good = 0.0;
  iel     = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;
    if ( (!mmgWarn0) && (MMG2D_quickcal(mesh,pt) < 0.0) ) {
      mmgWarn0 = 1;
      fprintf(stderr,"  ## Warning: %s: at least 1 negative area\n",__func__);
    }

    if ( !met->m ) {
      rap = MMG2D_ALPHAD * MMG2D_caltri_iso(mesh,met,pt);
    }
    else
      rap = MMG2D_ALPHAD * MMG2D_caltri(mesh,met,pt);

    if ( rap < rapmin ) {
      rapmin = rap;
      iel    = ok;
    }
    if ( rap > 0.5 )  med++;
    if ( rap > 0.12 ) good++;
    if ( rap < MMG2D_BADKAL )  mesh->info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

#ifndef DEBUG
  fprintf(stdout,"\n  -- MESH QUALITY   %" MMG5_PRId "\n",mesh->nt - nex);
  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%" MMG5_PRId ")\n",
          rapmax,rapavg / (mesh->nt-nex),rapmin,iel);
#else
  fprintf(stdout,"     BEST   %e  AVRG.   %e  WRST.   %e (%" MMG5_PRId ")\n => %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
          rapmax,rapavg / (mesh->nt-nex),rapmin,iel,
          MMG5_indPt(mesh,mesh->tria[iel].v[0]),MMG5_indPt(mesh,mesh->tria[iel].v[1]),
          MMG5_indPt(mesh,mesh->tria[iel].v[2]));
#endif

  /* print histo */
  fprintf(stdout,"     HISTOGRAMM:");
  fprintf(stdout,"  %6.2f %% > 0.12\n",100.0*(good/(float)(mesh->nt-nex)));
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"                  %6.2f %% >  0.5\n",100.0*( med/(float)(mesh->nt-nex)));
    imax = MG_MIN(4,(int)(5.*rapmax));
    for (i=imax; i>=(int)(5*rapmin); i--) {
      fprintf(stdout,"     %5.1f < Q < %5.1f   %7d   %6.2f %%\n",
              i/5.,i/5.+0.2,his[i],100.*(his[i]/(float)(mesh->nt-nex)));
    }
  }

  return  MMG5_minQualCheck(iel,rapmin,1.);
}
