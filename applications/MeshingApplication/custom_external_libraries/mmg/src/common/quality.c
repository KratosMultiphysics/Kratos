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
 * \file common/quality.c
 * \brief Functions to compute elements quality and edge lengths.
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
 * \param met pointer toward the meric structure.
 * \param pt pointer toward the triangle structure.
 * \return The computed quality.
 *
 * Compute the quality of the surface triangle \a ptt with respect to
 * an anisotropic metric and a classic storage of the ridges metrics.
 *
 */
double _MMG5_caltri33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt) {
  double   anisurf,dd,abx,aby,abz,acx,acy,acz,bcx,bcy,bcz;
  double  *a,*b,*c,*ma,*mb,*mc,m[6],l0,l1,l2,rap;
  int      ia,ib,ic;
  char     i;

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];

  ma = &met->m[6*ia];
  mb = &met->m[6*ib];
  mc = &met->m[6*ic];

  /* 2*area */
  anisurf  = _MMG5_surftri33_ani(mesh,pt,ma,mb,mc);
  if ( anisurf <= _MMG5_EPSD ) return(0.0);

  dd  = 1.0 / 3.0;
  for (i=0; i<6; i++)
    m[i] = dd * (ma[i] + mb[i] + mc[i]);

  a = &mesh->point[ia].c[0];
  b = &mesh->point[ib].c[0];
  c = &mesh->point[ic].c[0];

  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  /* length */
  l0 = m[0]*abx*abx + m[3]*aby*aby + m[5]*abz*abz
    + 2.0*(m[1]*abx*aby + m[2]*abx*abz + m[4]*aby*abz);

  l1 = m[0]*acx*acx + m[3]*acy*acy + m[5]*acz*acz
      + 2.0*(m[1]*acx*acy + m[2]*acx*acz + m[4]*acy*acz);

  l2 = m[0]*bcx*bcx + m[3]*bcy*bcy + m[5]*bcz*bcz
      + 2.0*(m[1]*bcx*bcy + m[2]*bcx*bcz + m[4]*bcy*bcz);

  rap = l0 + l1 + l2;

  /* quality = 2*area/length */
  if ( rap > _MMG5_EPSD ) {
    return( anisurf / rap);
  }
  else
    return(0.0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The computed quality.
 *
 * Compute the quality of the surface triangle \a ptt with respect to
 * an anisotropic metric.
 *
 * \warning The quality is computed as if the triangle is a "straight" triangle.
 *
 */
inline double _MMG5_caltri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  MMG5_pPoint   p[3];
  double        rap,anisurf,l0,l1,l2,m[6],mm[6];
  double        abx,aby,abz,acx,acy,acz,bcy,bcx,bcz;
  int           np[3],i,j;
  char          i1,i2;

  for (i=0; i<3; i++) {
    np[i] = ptt->v[i];
    p[i]  = &mesh->point[np[i]];
  }

  /* Set metric tensors at vertices of tria iel */
  for ( j=0; j<6; ++j) {
    mm[j] = 0;
  }

  for(i=0; i<3; i++) {

    if ( MG_SIN(p[i]->tag) || (MG_NOM & p[i]->tag) ) {
      memcpy(&m[0],&met->m[6*np[i]],6*sizeof(double));
    }
    else if ( p[i]->tag & MG_GEO ) {
      i1 = _MMG5_inxt2[i];
      i2 = _MMG5_iprv2[i];
      abx = 0.5*(p[i1]->c[0]+p[i2]->c[0]) - p[i]->c[0];
      aby = 0.5*(p[i1]->c[1]+p[i2]->c[1]) - p[i]->c[1];
      abz = 0.5*(p[i1]->c[2]+p[i2]->c[2]) - p[i]->c[2];
      if ( !_MMG5_buildridmet(mesh,met,np[i],abx,aby,abz,&m[0]) )  return(0.0);
    }
    else {
      memcpy(&m[0],&met->m[6*np[i]],6*sizeof(double));
    }

    for ( j=0; j<6; ++j) {
      mm[j] += _MMG5_ATHIRD*m[j];
    }
  }

  anisurf = _MMG5_surftri33_ani(mesh,ptt,mm,mm,mm);

  /* length */
  abx = p[1]->c[0] - p[0]->c[0];
  aby = p[1]->c[1] - p[0]->c[1];
  abz = p[1]->c[2] - p[0]->c[2];
  acx = p[2]->c[0] - p[0]->c[0];
  acy = p[2]->c[1] - p[0]->c[1];
  acz = p[2]->c[2] - p[0]->c[2];
  bcx = p[2]->c[0] - p[1]->c[0];
  bcy = p[2]->c[1] - p[1]->c[1];
  bcz = p[2]->c[2] - p[1]->c[2];


  l0 = mm[0]*abx*abx + mm[3]*aby*aby + mm[5]*abz*abz
    + 2.0*(mm[1]*abx*aby + mm[2]*abx*abz + mm[4]*aby*abz);

  l1 = mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz
      + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);

  l2 = mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz
      + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);

  rap = l0 + l1 + l2;

  if ( rap < _MMG5_EPSD ) return(0.0);

  /* quality = 2*area/length */
  return (anisurf / rap);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The computed quality.
 *
 * Compute the quality of the surface triangle \a ptt with respect to
 * an isotropic metric.
 *
 */
inline double _MMG5_caltri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  double   *a,*b,*c,cal,abx,aby,abz,acx,acy,acz,bcx,bcy,bcz,rap;

  a = &mesh->point[ptt->v[0]].c[0];
  b = &mesh->point[ptt->v[1]].c[0];
  c = &mesh->point[ptt->v[2]].c[0];

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  cal  = (aby*acz - abz*acy) * (aby*acz - abz*acy);
  cal += (abz*acx - abx*acz) * (abz*acx - abx*acz);
  cal += (abx*acy - aby*acx) * (abx*acy - aby*acx);

  if ( cal < _MMG5_EPSD2 )  return(0.0);

  /* qual = 2.*surf / length */
  rap  = abx*abx + aby*aby + abz*abz;
  rap += acx*acx + acy*acy + acz*acz;
  rap += bcx*bcx + bcy*bcy + bcz*bcz;

  if ( rap < _MMG5_EPSD2 )  return(0.0);

  return(sqrt(cal) / rap);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param ned edges number.
 * \param avlen pointer toward the average edges lengths.
 * \param amin index of first extremity of the smallest edge.
 * \param bmin index of second extremity of the smallest edge.
 * \param lmin smallest edge length.
 * \param amax index of first extremity of the largest edge.
 * \param bmax index of second extremity of the largest edge.
 * \param lmax largest edge length.
 * \param bd pointer toward the table of the quality span.
 * \param hl pointer toward the table that store the number of edges for each
 * span of quality
 *
 * Display histogram of edge length.
 *
 */
void _MMG5_displayHisto(MMG5_pMesh mesh, int ned, double *avlen,
                        int amin, int bmin, double lmin,
                        int amax, int bmax, double lmax, double *bd, int *hl )
{
  double dned;
  int    k;

  dned     = (double)ned;
  (*avlen) = (*avlen) / dned;

  fprintf(stdout,"\n  -- RESULTING EDGE LENGTHS  %d\n",ned);
  fprintf(stdout,"     AVERAGE LENGTH         %12.4f\n",(*avlen));
  fprintf(stdout,"     SMALLEST EDGE LENGTH   %12.4f   %6d %6d\n",
          lmin,amin,bmin);
  fprintf(stdout,"     LARGEST  EDGE LENGTH   %12.4f   %6d %6d \n",
          lmax,amax,bmax);
  if ( abs(mesh->info.imprim) < 3 ) return;

  if ( hl[2]+hl[3]+hl[4] )
    fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n",
            bd[2],bd[5],hl[2]+hl[3]+hl[4],100.*(hl[2]+hl[3]+hl[4])/(double)ned);
/*  if ( hl[3]+hl[4]+hl[5] ) */
/*     fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n", */
/*             bd[3],bd[6],hl[3]+hl[4]+hl[5],100.*(hl[3]+hl[4]+hl[5])/(double)ned); */

  if ( abs(mesh->info.imprim) < 4 ) return;

  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"\n     HISTOGRAMM:\n");
    if ( hl[0] )
      fprintf(stdout,"     0.00 < L < 0.30  %8d   %5.2f %%  \n",
              hl[0],100.*(hl[0]/(float)ned));
    if ( lmax > 0.2 ) {
      for (k=2; k<9; k++) {
        if ( hl[k-1] > 0 )
          fprintf(stdout,"   %6.2f < L <%5.2f  %8d   %5.2f %%  \n",
                  bd[k-1],bd[k],hl[k-1],100.*(hl[k-1]/(float)ned));
      }
      if ( hl[8] )
        fprintf(stdout,"     5.   < L         %8d   %5.2f %%  \n",
                hl[8],100.*(hl[8]/(float)ned));
    }
  }
}
