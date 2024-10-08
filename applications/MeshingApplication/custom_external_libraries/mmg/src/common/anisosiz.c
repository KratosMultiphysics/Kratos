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
 * \file common/anisosiz.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"
#include "mmgexterns_private.h"
#include "inlined_functions_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param m pointer toward the metric at triangle vertices.
 * \param ptt pointer toward the triangle structure.
 * \return The double of the triangle area.
 *
 * Compute the double of the area of the surface triangle \a ptt with respect to
 * the anisotropic metric \a m.
 *
 */
static inline
double MMG5_surf(MMG5_pMesh mesh,double m[3][6],MMG5_pTria ptt) {
  MMG5_Bezier    b;
  double         surf,dens,J[3][2],mJ[3][2],tJmJ[2][2];
  int8_t         i,nullDens;
  static int8_t  mmgErr=0;

  surf = 0.0;

  if ( !MMG5_bezierCP(mesh,ptt,&b,1) ) return 0.0;

  /* Compute density integrand of volume at the 3 vertices of T */
  nullDens = 0;
  for (i=0; i<3; i++) {
    if ( i == 0 ) {
      J[0][0] = 3.0*( b.b[7][0] - b.b[0][0] ) ; J[0][1] = 3.0*( b.b[6][0] - b.b[0][0] );
      J[1][0] = 3.0*( b.b[7][1] - b.b[0][1] ) ; J[1][1] = 3.0*( b.b[6][1] - b.b[0][1] );
      J[2][0] = 3.0*( b.b[7][2] - b.b[0][2] ) ; J[2][1] = 3.0*( b.b[6][2] - b.b[0][2] );
    }
    else if ( i == 1 ) {
      J[0][0] = 3.0*( b.b[1][0] - b.b[8][0] ) ; J[0][1] = 3.0*( b.b[3][0] - b.b[8][0] );
      J[1][0] = 3.0*( b.b[1][1] - b.b[8][1] ) ; J[1][1] = 3.0*( b.b[3][1] - b.b[8][1] );
      J[2][0] = 3.0*( b.b[1][2] - b.b[8][2] ) ; J[2][1] = 3.0*( b.b[3][2] - b.b[8][2] );
    }
    else {
      J[0][0] = 3.0*( b.b[4][0] - b.b[5][0] ) ; J[0][1] = 3.0*( b.b[2][0] - b.b[5][0] );
      J[1][0] = 3.0*( b.b[4][1] - b.b[5][1] ) ; J[1][1] = 3.0*( b.b[2][1] - b.b[5][1] );
      J[2][0] = 3.0*( b.b[4][2] - b.b[5][2] ) ; J[2][1] = 3.0*( b.b[2][2] - b.b[5][2] );
    }

    mJ[0][0] = m[i][0]*J[0][0] + m[i][1]*J[1][0] + m[i][2]*J[2][0];
    mJ[1][0] = m[i][1]*J[0][0] + m[i][3]*J[1][0] + m[i][4]*J[2][0];
    mJ[2][0] = m[i][2]*J[0][0] + m[i][4]*J[1][0] + m[i][5]*J[2][0];

    mJ[0][1] = m[i][0]*J[0][1] + m[i][1]*J[1][1] + m[i][2]*J[2][1];
    mJ[1][1] = m[i][1]*J[0][1] + m[i][3]*J[1][1] + m[i][4]*J[2][1];
    mJ[2][1] = m[i][2]*J[0][1] + m[i][4]*J[1][1] + m[i][5]*J[2][1];

    /* dens = sqrt(tJacsigma * M * Jacsigma )*/
    tJmJ[0][0] = J[0][0]*mJ[0][0] + J[1][0]*mJ[1][0] + J[2][0]*mJ[2][0];
    tJmJ[0][1] = J[0][0]*mJ[0][1] + J[1][0]*mJ[1][1] + J[2][0]*mJ[2][1];
    tJmJ[1][0] = J[0][1]*mJ[0][0] + J[1][1]*mJ[1][0] + J[2][1]*mJ[2][0];
    tJmJ[1][1] = J[0][1]*mJ[0][1] + J[1][1]*mJ[1][1] + J[2][1]*mJ[2][1];

    dens = tJmJ[0][0]*tJmJ[1][1] - tJmJ[1][0]*tJmJ[0][1];
    if ( dens <= MMG5_EPSD2 ) {
#ifndef NDEBUG
      if ( !mmgErr ) {
        fprintf(stderr,"\n  ## Warning: %s: at least 1 negative or null density.\n",
                __func__);
        mmgErr = 1;
      }
#endif
      ++nullDens;
      continue;
    }
    surf += sqrt(fabs(dens));
  }

  if ( nullDens==3 ) return 0;

  surf *= MMG5_ATHIRD;
  return surf;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The double of the triangle area.
 *
 * Compute the double of the area of the surface triangle \a ptt with respect to
 * the anisotropic metric \a met (for special storage of ridges metrics).
 *
 */
double MMG5_surftri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  MMG5_pPoint    p[3];
  MMG5_int       np[3];
  double         ux,uy,uz,m[3][6],rbasis[3][3];
  int8_t         i1,i2;
  int            i;

  for (i=0; i<3; i++) {
    np[i] = ptt->v[i];
    p[i]  = &mesh->point[np[i]];
  }

  /* Set metric tensors at vertices of tria iel */
  for(i=0; i<3; i++) {

    if ( MG_SIN(p[i]->tag) || (MG_NOM & p[i]->tag) ) {
      memcpy(&m[i][0],&met->m[6*np[i]],6*sizeof(double));
    }
    else if ( p[i]->tag & MG_GEO ) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      ux = 0.5*(p[i1]->c[0]+p[i2]->c[0]) - p[i]->c[0];
      uy = 0.5*(p[i1]->c[1]+p[i2]->c[1]) - p[i]->c[1];
      uz = 0.5*(p[i1]->c[2]+p[i2]->c[2]) - p[i]->c[2];
      /* Note that rbasis is unused in this function */
      if ( !MMG5_buildridmet(mesh,met,np[i],ux,uy,uz,&m[i][0],rbasis) )  return 0.0;
    }
    else {
      memcpy(&m[i][0],&met->m[6*np[i]],6*sizeof(double));
    }
  }
  return MMG5_surf(mesh,m,ptt);

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param ptt pointer toward the triangle structure.
 * \param ma metric at triangle vertex.
 * \param mb metric at triangle vertex.
 * \param mc metric at triangle vertex.
 * \return The double of the triangle area.
 *
 * Compute the double of the area of the surface triangle \a ptt with respect to
 * the anisotropic metric \a met (for classic storage of ridges metrics).
 *
 */
double MMG5_surftri33_ani(MMG5_pMesh mesh,MMG5_pTria ptt,
                           double ma[6], double mb[6], double mc[6]) {
  double         mm[6];
  double         *a,*b,*c,abx,aby,abz,acx,acy,acz,dens[3],surf;
  int            i;
  MMG5_int       ia,ib,ic;

  ia = ptt->v[0];
  ib = ptt->v[1];
  ic = ptt->v[2];

  a = &mesh->point[ia].c[0];
  b = &mesh->point[ib].c[0];
  c = &mesh->point[ic].c[0];

  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  /* Compute the mean of the metrics over the triangle */
  for (i=0; i<6; i++)
    mm[i] = MMG5_ATHIRD * (ma[i] + mb[i]+ mc[i]);

  /* Compute sqrt(det(t^JmJ))  (= int_T sqrt(t^JmJ) for a non-curve element) */
  dens[0] = (abx*abx*mm[0]+abx*aby*mm[1]+abx*abz*mm[2])
    + (aby*abx*mm[1]+aby*aby*mm[3]+aby*abz*mm[4])
    + (abz*abx*mm[2]+abz*aby*mm[4]+abz*abz*mm[5]);

  dens[1] = (abx*acx*mm[0]+abx*acy*mm[1]+abx*acz*mm[2])
    + (aby*acx*mm[1]+aby*acy*mm[3]+aby*acz*mm[4])
    + (abz*acx*mm[2]+abz*acy*mm[4]+abz*acz*mm[5]);

  dens[2] = (acx*acx*mm[0]+acx*acy*mm[1]+acx*acz*mm[2])
    + (acy*acx*mm[1]+acy*acy*mm[3]+acy*acz*mm[4])
    + (acz*acx*mm[2]+acz*acy*mm[4]+acz*acz*mm[5]);

  surf = dens[0]*dens[2]-dens[1]*dens[1];

  if ( surf < MMG5_EPSD ) return 0.0;

  surf = sqrt(surf);

  return surf;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ismet 1 if user provided metric.
 *
 * Search for points with unintialized metric and define anisotropic size at
 * this points.
 *
 */
void MMG5_defUninitSize(MMG5_pMesh mesh,MMG5_pSol met,int8_t ismet )
{
  MMG5_pPoint   ppt;
  double        *m,*n,r[3][3],isqhmax;
  MMG5_int      k;

  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) || ppt->flag > 0 )  continue;

    m = &met->m[6*k];
    if(ismet) {
      if ( !(MG_SIN(ppt->tag) || (ppt->tag & MG_NOM)) && (ppt->tag & MG_GEO) ) {
        m[0] = m[1] = m[2] = m[3] = m[4] = isqhmax;
        m[5] = 0;
      }

      ppt->flag = 1;
      continue;
    }

    memset(m,0,6*sizeof(double));
    if (  (MG_SIN(ppt->tag) || (MG_NOM & ppt->tag)) ) {
      m[0] = m[3] = m[5] = isqhmax;
    }
    else if ( ppt->tag & MG_GEO ) {
      /* We store the size in the tangent dir in m[0], in the n1 dir in m[1] and
       * in the n2 dir in m[2]. */
      m[0] = m[1] = m[2] = m[3] = m[4] = isqhmax;
    }
    else {
      n = ppt->tag & MG_REF ? &mesh->xpoint[ppt->xp].n1[0] : ppt->n;
      MMG5_rotmatrix(n,r);
      m[0] = isqhmax*(r[0][0]*r[0][0]+r[1][0]*r[1][0]+r[2][0]*r[2][0]);
      m[1] = isqhmax*(r[0][0]*r[0][1]+r[1][0]*r[1][1]+r[2][0]*r[2][1]);
      m[2] = isqhmax*(r[0][0]*r[0][2]+r[1][0]*r[1][2]+r[2][0]*r[2][2]);
      m[3] = isqhmax*(r[0][1]*r[0][1]+r[1][1]*r[1][1]+r[2][1]*r[2][1]);
      m[4] = isqhmax*(r[0][1]*r[0][2]+r[1][1]*r[1][2]+r[2][1]*r[2][2]);
      m[5] = isqhmax*(r[0][2]*r[0][2]+r[1][2]*r[1][2]+r[2][2]*r[2][2]);
    }
    ppt->flag = 2;
  }
}

/**
 * \param k index of the tetrahedra from which we come.
 * \param p0 pointer toward the point on which we want to def the metric.
 * \param i0 pointer toward the local index of the point in tria.
 * \param b control polygon of triangle.
 * \param r rotation matrix.
 * \param c physical coordinates of the curve edge mid-point.
 * \param lispoi list of incident vertices to p0
 * \param tAA matrix to fill
 * \param tAb second member
 *
 * Fill matrice \sum tAA and second member \sum tAb with \f$ A=( X_{P_i}^2
 * Y_{P_i}^2 X_{P_i}Y_{P_i}) \f$ and \f$ b= Z_{P_i}\f$ with P_i the
 * physical points at edge [i0;i1] extremities and middle.  Compute the physical
 * coor \a c of the curve edge's mid-point for a regular or reference point.
 *
 */
void MMG5_fillDefmetregSys( MMG5_int k, MMG5_pPoint p0, int i0, MMG5_Bezier b,
                             double r[3][3], double c[3],
                             double *lispoi,
                             double tAA[6], double tAb[3] )
{
  double b0[3],b1[3],d[3];
  int    j;

  for(j=0; j<10; j++){
    c[0] = b.b[j][0] - p0->c[0];
    c[1] = b.b[j][1] - p0->c[1];
    c[2] = b.b[j][2] - p0->c[2];

    b.b[j][0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
    b.b[j][1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
    b.b[j][2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];
  }

/* Mid-point along edge [i0;i1] and endpoint in the rotated frame */
  if(i0 == 0){
    memcpy(b0,&(b.b[7][0]),3*sizeof(double));
    memcpy(b1,&(b.b[8][0]),3*sizeof(double));
  }
  else if(i0 == 1){
    memcpy(b0,&(b.b[3][0]),3*sizeof(double));
    memcpy(b1,&(b.b[4][0]),3*sizeof(double));
  }
  else{
    memcpy(b0,&(b.b[5][0]),3*sizeof(double));
    memcpy(b1,&(b.b[6][0]),3*sizeof(double));
  }

/* At this point, the two control points b0 and b1 are expressed in the
 * rotated frame. We compute the physical coor of the curve edge's
 * mid-point. */
  c[0] = 3.0/8.0*b0[0] + 3.0/8.0*b1[0] + 1.0/8.0*lispoi[3*k+1];
  c[1] = 3.0/8.0*b0[1] + 3.0/8.0*b1[1] + 1.0/8.0*lispoi[3*k+2];
  c[2] = 3.0/8.0*b0[2] + 3.0/8.0*b1[2] + 1.0/8.0*lispoi[3*k+3];

/* Fill matrice \sum tAA and second member \sum tAb with \f$A=( X_{P_i}^2
 * Y_{P_i}^2 X_{P_i}Y_{P_i})\f$ and \f$b= Z_{P_i}\f$ with P_i the
 * physical points at edge [i0;i1] extremities and middle. */
  tAA[0] += c[0]*c[0]*c[0]*c[0];
  tAA[1] += c[0]*c[0]*c[1]*c[1];
  tAA[2] += c[0]*c[0]*c[0]*c[1];
  tAA[3] += c[1]*c[1]*c[1]*c[1];
  tAA[4] += c[0]*c[1]*c[1]*c[1];
  tAA[5] += c[0]*c[0]*c[1]*c[1];

  tAb[0] += c[0]*c[0]*c[2];
  tAb[1] += c[1]*c[1]*c[2];
  tAb[2] += c[0]*c[1]*c[2];

  tAA[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1];
  tAA[1] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];
  tAA[2] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2];
  tAA[3] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
  tAA[4] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+2];
  tAA[5] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+2];

  tAb[0] += lispoi[3*k+1]*lispoi[3*k+1]*lispoi[3*k+3];
  tAb[1] += lispoi[3*k+2]*lispoi[3*k+2]*lispoi[3*k+3];
  tAb[2] += lispoi[3*k+1]*lispoi[3*k+2]*lispoi[3*k+3];

/* Mid-point along median edge (coor of mid-point in parametric space : (1/4
 * 1/4)) and endpoint in the rotated frame (coor of end-point in parametric
 * space (1/2 1/2)). */
  if(i0 == 0){
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
  else if(i0 == 1){
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
  else{
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

/* Fill matrice tAA and second member tAb*/
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

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param r pointer toward the rotation matrix.
 * \param c physical coordinates of the curve edge mid-point.
 * \param tAA matrix of the system to solve.
 * \param tAb second member.
 * \param m pointer toward the metric.
 * \param isqhmax maximum size for edge.
 * \param isqhmin minimum size for edge.
 * \param hausd hausdorff value at point.
 * \return 1 if success, 0 if fail.
 *
 * Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) for a regular
 * point.
 *
 */
int MMG5_solveDefmetregSys( MMG5_pMesh mesh, double r[3][3], double c[3],
                             double tAA[6], double tAb[3], double *m,
                             double isqhmin, double isqhmax, double hausd)
{
  double intm[3], kappa[2], vp[2][2], b0[3], b1[3], b2[3];
  static int mmgWarn0=0;

  memset(intm,0x0,3*sizeof(double));

  /* case planar surface : tAb = 0 => no curvature */
  /* isotropic metric with hmax size*/

  if((tAb[0]*tAb[0] + tAb[1]*tAb[1] + tAb[2]*tAb[2]) < MMG5_EPSD) {
    m[0] = isqhmax;
    m[1] = 0;
    m[2] = 0;
    m[3] = isqhmax;
    m[4] = 0;
    m[5] = isqhmax;
    return 1;
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !MMG5_sys33sym(tAA,tAb,c) ) {
    if ( !mmgWarn0 ) {
       fprintf(stderr,"\n  ## Warning: %s: unable to solve the system on at"
               " least 1 point.\n",__func__);
      mmgWarn0 = 1;
    }
    return 0;
  }
  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
     are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
  MMG5_eigensym(intm,kappa,vp);

  /* Truncation of eigenvalues */
  kappa[0] = 2.0/9.0 * fabs(kappa[0])/hausd;
  kappa[0] = MG_MIN(kappa[0],isqhmin);
  kappa[0] = MG_MAX(kappa[0],isqhmax);

  kappa[1] = 2.0/9.0 * fabs(kappa[1])/hausd;
  kappa[1] = MG_MIN(kappa[1],isqhmin);
  kappa[1] = MG_MAX(kappa[1],isqhmax);

  /* Send back the metric to the canonical basis of tangent plane :
     diag(lambda) = {^t}vp * M * vp, M = vp * diag(lambda) * {^t}vp */
  intm[0] = kappa[0]*vp[0][0]*vp[0][0] + kappa[1]*vp[1][0]*vp[1][0];
  intm[1] = kappa[0]*vp[0][0]*vp[0][1] + kappa[1]*vp[1][0]*vp[1][1];
  intm[2] = kappa[0]*vp[0][1]*vp[0][1] + kappa[1]*vp[1][1]*vp[1][1];

  /* At this point, intm (with 0 replaced by isqhmax in the z direction) is the
     desired metric, except it is expressed in the rotated bc, that is intm = R
     * metric in bc * ^t R, so metric in bc = ^tR*intm*R */
  /* b0, b1 and b2 are the lines of matrix intm*R  */
  // intm = intm[0]  intm[1]    0
  //        intm[1]  intm[2]    0
  //           0       0     isqhmax

  b0[0] = intm[0]*r[0][0] + intm[1]*r[1][0];
  b0[1] = intm[0]*r[0][1] + intm[1]*r[1][1];
  b0[2] = intm[0]*r[0][2] + intm[1]*r[1][2];
  b1[0] = intm[1]*r[0][0] + intm[2]*r[1][0];
  b1[1] = intm[1]*r[0][1] + intm[2]*r[1][1];
  b1[2] = intm[1]*r[0][2] + intm[2]*r[1][2];
  b2[0] = isqhmax*r[2][0];
  b2[1] = isqhmax*r[2][1];
  b2[2] = isqhmax*r[2][2];

  m[0] = r[0][0] * b0[0] + r[1][0] * b1[0] + r[2][0] * b2[0];
  m[1] = r[0][0] * b0[1] + r[1][0] * b1[1] + r[2][0] * b2[1];
  m[2] = r[0][0] * b0[2] + r[1][0] * b1[2] + r[2][0] * b2[2];

  m[3] = r[0][1] * b0[1] + r[1][1] * b1[1] + r[2][1] * b2[1];
  m[4] = r[0][1] * b0[2] + r[1][1] * b1[2] + r[2][1] * b2[2];

  m[5] = r[0][2] * b0[2] + r[1][2] * b1[2] + r[2][2] * b2[2];

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param p0 pointer toward the point on which we want to define the metric.
 * \param ipref table containing the indices of the edge extremities.
 * \param r pointer toward the rotation matrix.
 * \param c physical coordinates of the curve edge mid-point.
 * \param tAA matrix of the system to solve.
 * \param tAb second member.
 * \param m pointer toward the metric.
 * \param isqhmax maximum size for edge.
 * \param isqhmin minimum size for edge.
 * \param hausd hausdorff value at point (unused).
 * \return 1 if success, 0 if fail.
 *
 * Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) for a ref
 * point.
 *
 */
int MMG5_solveDefmetrefSys( MMG5_pMesh mesh, MMG5_pPoint p0, MMG5_int ipref[2],
                             double r[3][3], double c[3],
                             double tAA[6], double tAb[3], double *m,
                             double isqhmin, double isqhmax, double hausd)
{
  MMG5_pPoint   p1;
  double        intm[3], kappa[2], vp[2][2], b0[3], b1[3], b2[3], kappacur;
  double        gammasec[3],tau[2], ux, uy, uz, ps1, l, ll, *t, *t1;
  int           i;
  static int8_t mmgWarn=0;

  (void)hausd;

  memset(intm,0x0,3*sizeof(double));

  /* case planar surface : tAb = 0 => no curvature */
  /* isotropic metric with hmax size*/
  if((tAb[0]*tAb[0] + tAb[1]*tAb[1] + tAb[2]*tAb[2]) < MMG5_EPSD) {
    m[0] = isqhmax;
    m[1] = 0;
    m[2] = 0;
    m[3] = isqhmax;
    m[4] = 0;
    m[5] = isqhmax;
    return 1;
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !MMG5_sys33sym(tAA,tAb,c) ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to solve the system on at"
              " least 1 point.\n", __func__);
      mmgWarn = 1;
    }
    return 0;
  }
  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach :
     vp[0] and vp[1] are the two pr. directions of curvature, and the two
     curvatures can be inferred from lambdas*/
  MMG5_eigensym(intm,kappa,vp);

  /* Truncation of eigenvalues */
  kappa[0] = 2.0/9.0 * fabs(kappa[0])/mesh->info.hausd;
  kappa[0] = MG_MIN(kappa[0],isqhmin);
  kappa[0] = MG_MAX(kappa[0],isqhmax);

  kappa[1] = 2.0/9.0 * fabs(kappa[1])/mesh->info.hausd;
  kappa[1] = MG_MIN(kappa[1],isqhmin);
  kappa[1] = MG_MAX(kappa[1],isqhmax);

  /* Send back the metric to the canonical basis of tangent plane :
     diag(lambda) = {^t}vp * M * vp, M = vp * diag(lambda) * {^t}vp */
  intm[0] = kappa[0]*vp[0][0]*vp[0][0] + kappa[1]*vp[1][0]*vp[1][0];
  intm[1] = kappa[0]*vp[0][0]*vp[0][1] + kappa[1]*vp[1][0]*vp[1][1];
  intm[2] = kappa[0]*vp[0][1]*vp[0][1] + kappa[1]*vp[1][1]*vp[1][1];

  /* Now express metric with respect to the approx of the underlying ref
   * curve */
  t = &p0->n[0];
  kappacur = 0.0;

  for (i=0; i<2; i++) {
    p1 = &mesh->point[ipref[i]];
    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    ps1 =  ux*t[0] + uy*t[1] + uz*t[2];
    c[0] = MMG5_ATHIRD*ps1*t[0];
    c[1] = MMG5_ATHIRD*ps1*t[1];
    c[2] = MMG5_ATHIRD*ps1*t[2];

    b0[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
    b0[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
    b0[2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];

    if ( (MG_CRN & p1->tag) || (MG_NOM & p1->tag) ) {
      c[0] = p1->c[0] - MMG5_ATHIRD*ux;
      c[1] = p1->c[1] - MMG5_ATHIRD*uy;
      c[2] = p1->c[2] - MMG5_ATHIRD*uz;
    }
    else {
      assert(MG_REF & p1->tag);
      t1 = &(p1->n[0]);
      ps1 =  -(ux*t1[0] + uy*t1[1] + uz*t1[2]);
      c[0] = p1->c[0] + MMG5_ATHIRD*ps1*t1[0];
      c[1] = p1->c[1] + MMG5_ATHIRD*ps1*t1[1];
      c[2] = p1->c[2] + MMG5_ATHIRD*ps1*t1[2];
    }
    c[0] -= p0->c[0];
    c[1] -= p0->c[1];
    c[2] -= p0->c[2];

    b1[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
    b1[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
    b1[2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];

    /* Everything is expressed in the rotated frame */
    tau[0] = 3.0*b0[0];
    tau[1] = 3.0*b0[1];
    ll = tau[0]*tau[0] + tau[1]*tau[1];
    if ( ll < MMG5_EPSD ) {
      kappacur = isqhmax;
      continue;
    }
    l = 1.0 / sqrt(ll);
    tau[0] *= l;
    tau[1] *= l;

    gammasec[0] = -12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = -12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = -12.0*b0[2] + 6.0*b1[2];

    ps1 = tau[0]*gammasec[0] + tau[1]*gammasec[1];
    c[0] = gammasec[0] - ps1*tau[0];
    c[1] = gammasec[1] - ps1*tau[1];
    c[2] = gammasec[2];

    // p.s. with normal at p0
    kappacur = MG_MAX(kappacur,MG_MAX(0.0,1.0/ll*fabs(c[2])));
  }

  /* Rotation of tangent vector : tau is reused */
  c[0] =  r[0][0]*t[0] + r[0][1]*t[1] + r[0][2]*t[2];
  c[1] =  r[1][0]*t[0] + r[1][1]*t[1] + r[1][2]*t[2];
  c[2] =  r[2][0]*t[0] + r[2][1]*t[1] + r[2][2]*t[2];
  memcpy(tau,&c[0],2*sizeof(double));

  /* Truncation of curvature */
  kappacur = 1.0/8.0 * kappacur/mesh->info.hausd;
  kappacur = MG_MIN(kappacur,isqhmin);
  kappacur = MG_MAX(kappacur,isqhmax);

  /* The associated matrix in basis (rt, orth rt) */
  c[0] = kappacur*tau[0]*tau[0] + isqhmax*tau[1]*tau[1];
  c[1] = (kappacur - isqhmax)*tau[0]*tau[1];
  c[2] = kappacur*tau[1]*tau[1] + isqhmax*tau[0]*tau[0];

  /* Reuse b0 for commodity */
  MMG5_intmetsavedir(mesh,c,intm,b0);
  memcpy(intm,b0,3*sizeof(double));

  /* At this point, intm (with 0 replaced by isqhmax in the z direction) is the
     desired metric, except it is expressed in the rotated bc, that is intm = R
     * metric in bc * ^t R, so metric in bc = ^tR*intm*R */
  // intm = intm[0]  intm[1]    0
  //        intm[1]  intm[2]    0
  //           0       0     isqhmax

  /* b0 and b1 serve now for nothing : let them be the lines of matrix intm*R */
  b0[0] = intm[0]*r[0][0] + intm[1]*r[1][0];
  b0[1] = intm[0]*r[0][1] + intm[1]*r[1][1];
  b0[2] = intm[0]*r[0][2] + intm[1]*r[1][2];
  b1[0] = intm[1]*r[0][0] + intm[2]*r[1][0];
  b1[1] = intm[1]*r[0][1] + intm[2]*r[1][1];
  b1[2] = intm[1]*r[0][2] + intm[2]*r[1][2];
  b2[0] = isqhmax*r[2][0];
  b2[1] = isqhmax*r[2][1];
  b2[2] = isqhmax*r[2][2];

  m[0] = r[0][0] * b0[0] + r[1][0] * b1[0] + r[2][0] * b2[0];
  m[1] = r[0][0] * b0[1] + r[1][0] * b1[1] + r[2][0] * b2[1];
  m[2] = r[0][0] * b0[2] + r[1][0] * b1[2] + r[2][0] * b2[2];

  m[3] = r[0][1] * b0[1] + r[1][1] * b1[1] + r[2][1] * b2[1];
  m[4] = r[0][1] * b0[2] + r[1][1] * b1[2] + r[2][1] * b2[2];

  m[5] = r[0][2] * b0[2] + r[1][2] * b1[2] + r[2][2] * b2[2];

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param p0 pointer toward the point at which we define the metric.
 * \param idp global index of the point at which we define the metric.
 * \param iprid pointer toward the two extremities of the ridge.
 * \param isqhmin minimum edge size.
 * \param isqhmax maximum edge size.
 * \return the computed ridge size in the tangent direction.
 *
 * Compute the specific size that we want to apply to a ridge in the direction
 * of the tangent of the ridge.  This wanted size is computed as using the
 * majoration of the haudorff distance between a triangle and its ideal curve
 * approximation by the Hessian of the signed distance function to the ideal
 * surface.
 *
 *
 \f[
 d^H(\partial \Omega,S_T) \leq \displaystyle\frac{1}{2}\left(
 \frac{d-1}{d} \right)^2 \max\limits_{T\in S_T} \max\limits_{x \in T
 } \max\limits_{y,z \in T} \langle \left| H(d_\omega)(x) \right|
 yz,yz\rangle.
 \f]
 *
 * where \f$ d^H(\partial \Omega,S_T) \f$ is the distance between the
 * triangle \f$S_T\f$ and the ideal boundary (reconstructed using cubic Bezier
 * patches) \f$ \partial \Omega \f$, \f$ d\f$ is the mesh dimension and \f$
 * H(d_\Omega) \f$ is the hessian matrix of the signed distance function to \f$
 * \Omega \f$.
 *
 * For all \f$ x \in \partial \Omega \f$, \f$H(d_\omega)(x)\f$ is the second
 * fundamental form whose eigenvalues are the principal curvatures (\f$\kappa_1
 * \f$ and \f$ \kappa_2 \f$) of \f$ \partial \Omega \f$ at \f$ x \f$ so the
 * previous formula can be rewritten as:
 \f[
 d^H(\partial \Omega,S_T) \leq \displaystyle\frac{1}{2}\left(
 \frac{d-1}{d} \right)^2 \max( \left|\kappa_1\right|,\left|\kappa_2\right|).
 \f]
 *
 * As we want to respect the imposed the \a hausd threshold and we are
 * interessed to get the wanted size along the tangent direction at the ridge
 * point only, finally, m is computed as
 *
 \f[
 m = \kappa \frac{1}{2\textrm{hausd}}
 \left(\frac{d-1}{d} \right)^2.
 \f]
 *
 * See Theorem 1 of \cite dapogny2014three.
 *
 **/
double MMG5_ridSizeInTangentDir(MMG5_pMesh mesh, MMG5_pPoint p0, MMG5_int idp,
                                 MMG5_int* iprid, double isqhmin,double isqhmax)
{
  int    i;
  double n0[3],tau[3],gammasec[3],c[3],ps,ll,l,m;
  double b0[3],b1[3],kappacur;

  m = isqhmax;
  for (i=0; i<2; i++) {
    // Remark: bezierEdge don't use n0 in case of a ridge so it's ok to call it
    // with an undefined n0.
    MMG5_bezierEdge(mesh,idp,iprid[i],b0,b1,1,n0);

    /* tau is the bezier curve derivative in p0 (parametric coor t=0) */
    tau[0] = 3.0*(b0[0] - p0->c[0]);
    tau[1] = 3.0*(b0[1] - p0->c[1]);
    tau[2] = 3.0*(b0[2] - p0->c[2]);
    ll = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
    if ( ll < MMG5_EPSD )  continue;
    l = 1.0 / sqrt(ll);
    tau[0] *= l;
    tau[1] *= l;
    tau[2] *= l;

    /* Acceleration of the bezier curve in p0 (parameteric coor t=0) in
     * canonical basis */
    gammasec[0] = 6.0*p0->c[0] -12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = 6.0*p0->c[1] -12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = 6.0*p0->c[2] -12.0*b0[2] + 6.0*b1[2];

    /* Scalar component of acceleration along tangent field (ie along velocity) */
    ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];

    /* Normal component of acceleration along velocity (ps*tau[] being the vector
     * component along velocity) */
    c[0] = gammasec[0] - ps*tau[0];
    c[1] = gammasec[1] - ps*tau[1];
    c[2] = gammasec[2] - ps*tau[2];

    /* Evaluation of curvature (derivative of unit tangent vector, i.e. change
     * in direction): be T the unit tangent vector, s. a. T = v/||v||, kappa =
     * ||T'||/||v|| = sqrt(c[0]) /||v||^2 */
    kappacur = MG_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));

    /* Computation of size from the mesh dimension and the curvature */
    kappacur = 1.0/8.0*kappacur/mesh->info.hausd;
    kappacur = MG_MIN(kappacur,isqhmin);
    kappacur = MG_MAX(kappacur,isqhmax);
    m = MG_MAX(m,kappacur);
  }
  return m;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param i0 local index in the face of the point on which we want to compute
 * the metric
 * \param bcu pointer toward the barycentric coordinates of vector \a u in the
 * looked face.
 * \param b bezier control polygon for the looked face.
 * \param isqhmin minimum edge size.
 * \param isqhmax maximum edge size.
 * \return the computed ridge size in first or second normal direction
 * (depending of i0).
 *
 * Compute the specific size of a ridge in the direction of the normal of the
 * looked face.  This wanted size is computed as using the
 * majoration of the haudorff distance between a triangle and its ideal curve
 * approximation by the Hessian of the signed distance function to the ideal
 * surface. See documentation of \ref MMG5_ridSizeInTangentDir function.
 *
 **/
double MMG5_ridSizeInNormalDir(MMG5_pMesh mesh,int i0,double* bcu,
                                MMG5_Bezier *b,double isqhmin,double isqhmax)
{
  double lambda[2],Jacb[3][2],Hb[3][3],tau[3],ll,l,gammasec[3],c[3];
  double ps,kappacur;

  if ( i0 == 0 ) { // w = 1, u,v = 0
    lambda[0] = bcu[1];
    lambda[1] = bcu[2];

    /* Jacb[i] = Jacobian matrix of i-th component of b at point p0.
       Jacb[i] = (\partial_u S_p0[i], \partial_v S_p0[i] ) */
    Jacb[0][0] = 3.0*(b->b[7][0]-b->b[0][0]);
    Jacb[1][0] = 3.0*(b->b[7][1]-b->b[0][1]);
    Jacb[2][0] = 3.0*(b->b[7][2]-b->b[0][2]);

    Jacb[0][1] = 3.0*(b->b[6][0]-b->b[0][0]);
    Jacb[1][1] = 3.0*(b->b[6][1]-b->b[0][1]);
    Jacb[2][1] = 3.0*(b->b[6][2]-b->b[0][2]);

    /* Hb[i] = hessian matrix of i-th component of b at point p0.
       Hb[i] = (\partial_u^2 S_p0[i], \partial_v \partial_u S_p0[i],\partial_v^2 S_p0[i]) */
    Hb[0][0] = 6.0*(b->b[0][0] - 2.0*b->b[7][0] + b->b[8][0]);
    Hb[1][0] = 6.0*(b->b[0][1] - 2.0*b->b[7][1] + b->b[8][1]);
    Hb[2][0] = 6.0*(b->b[0][2] - 2.0*b->b[7][2] + b->b[8][2]);

    Hb[0][1] = 6.0*(b->b[0][0] - b->b[7][0] - b->b[6][0] + b->b[9][0]);
    Hb[1][1] = 6.0*(b->b[0][1] - b->b[7][1] - b->b[6][1] + b->b[9][1]);
    Hb[2][1] = 6.0*(b->b[0][2] - b->b[7][2] - b->b[6][2] + b->b[9][2]);

    Hb[0][2] = 6.0*(b->b[0][0] - 2.0*b->b[6][0] + b->b[5][0]);
    Hb[1][2] = 6.0*(b->b[0][1] - 2.0*b->b[6][1] + b->b[5][1]);
    Hb[2][2] = 6.0*(b->b[0][2] - 2.0*b->b[6][2] + b->b[5][2]);
  }
  else if ( i0 == 1 ) {  //w = v = 0, u = 1;
    lambda[0] = bcu[0];
    lambda[1] = bcu[1];

    Jacb[0][0] = 3.0*(b->b[1][0]-b->b[8][0]);
    Jacb[1][0] = 3.0*(b->b[1][1]-b->b[8][1]);
    Jacb[2][0] = 3.0*(b->b[1][2]-b->b[8][2]);

    Jacb[0][1] = 3.0*(b->b[3][0]-b->b[8][0]);
    Jacb[1][1] = 3.0*(b->b[3][1]-b->b[8][1]);
    Jacb[2][1] = 3.0*(b->b[3][2]-b->b[8][2]);

    Hb[0][0] = 6.0*(b->b[1][0] - 2.0*b->b[8][0] + b->b[7][0]);
    Hb[1][0] = 6.0*(b->b[1][1] - 2.0*b->b[8][1] + b->b[7][1]);
    Hb[2][0] = 6.0*(b->b[1][2] - 2.0*b->b[8][2] + b->b[7][2]);

    Hb[0][1] = 6.0*(b->b[7][0] - b->b[8][0] - b->b[9][0] + b->b[3][0]);
    Hb[1][1] = 6.0*(b->b[7][1] - b->b[8][1] - b->b[9][1] + b->b[3][1]);
    Hb[2][1] = 6.0*(b->b[7][2] - b->b[8][2] - b->b[9][2] + b->b[3][2]);

    Hb[0][2] = 6.0*(b->b[4][0] - 2.0*b->b[9][0] + b->b[7][0]);
    Hb[1][2] = 6.0*(b->b[4][1] - 2.0*b->b[9][1] + b->b[7][1]);
    Hb[2][2] = 6.0*(b->b[4][2] - 2.0*b->b[9][2] + b->b[7][2]);
  }
  else {   //w =u = 0, v =1
    lambda[0] = bcu[2];
    lambda[1] = bcu[0];

    Jacb[0][0] = 3.0*(b->b[4][0]-b->b[5][0]);
    Jacb[1][0] = 3.0*(b->b[4][1]-b->b[5][1]);
    Jacb[2][0] = 3.0*(b->b[4][2]-b->b[5][2]);

    Jacb[0][1] = 3.0*(b->b[2][0]-b->b[5][0]);
    Jacb[1][1] = 3.0*(b->b[2][1]-b->b[5][1]);
    Jacb[2][1] = 3.0*(b->b[2][2]-b->b[5][2]);

    Hb[0][0] = 6.0*(b->b[3][0] - 2.0*b->b[9][0] + b->b[6][0]);
    Hb[1][0] = 6.0*(b->b[3][1] - 2.0*b->b[9][1] + b->b[6][1]);
    Hb[2][0] = 6.0*(b->b[3][2] - 2.0*b->b[9][2] + b->b[6][2]);

    Hb[0][1] = 6.0*(b->b[4][0] - b->b[5][0] - b->b[9][0] + b->b[6][0]);
    Hb[1][1] = 6.0*(b->b[4][1] - b->b[5][1] - b->b[9][1] + b->b[6][1]);
    Hb[2][1] = 6.0*(b->b[4][2] - b->b[5][2] - b->b[9][2] + b->b[6][2]);

    Hb[0][2] = 6.0*(b->b[2][0] - 2.0*b->b[5][0] + b->b[6][0]);
    Hb[1][2] = 6.0*(b->b[2][1] - 2.0*b->b[5][1] + b->b[6][1]);
    Hb[2][2] = 6.0*(b->b[2][2] - 2.0*b->b[5][2] + b->b[6][2]);
  }

  /* tau = jacb *(lambda0,lambda1)*/
  tau[0] = Jacb[0][0]*lambda[0] + Jacb[0][1]*lambda[1];
  tau[1] = Jacb[1][0]*lambda[0] + Jacb[1][1]*lambda[1];
  tau[2] = Jacb[2][0]*lambda[0] + Jacb[2][1]*lambda[1];
  ll = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
  if ( ll < MMG5_EPSD )  return 0;

  l = 1.0 / sqrt(ll);
  tau[0] *= l;
  tau[1] *= l;
  tau[2] *= l;

  gammasec[0] = Hb[0][0]*lambda[0]*lambda[0] + 2.0*Hb[0][1]*lambda[0]*lambda[1] + Hb[0][2]*lambda[1]*lambda[1];
  gammasec[1] = Hb[1][0]*lambda[0]*lambda[0] + 2.0*Hb[1][1]*lambda[0]*lambda[1] + Hb[1][2]*lambda[1]*lambda[1];
  gammasec[2] = Hb[2][0]*lambda[0]*lambda[0] + 2.0*Hb[2][1]*lambda[0]*lambda[1] + Hb[2][2]*lambda[1]*lambda[1];

  ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];

  /* Normal component of acceleration along velocity (ps*tau[] being the vector
   * component along velocity) */
  c[0] = gammasec[0] - ps*tau[0];
  c[1] = gammasec[1] - ps*tau[1];
  c[2] = gammasec[2] - ps*tau[2];

  /* Evaluation of curvature (derivative of unit tangent vector, i.e. change
   * in direction) */
  kappacur = MG_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));

  /* Computation of size from the mesh dimension and the curvature */
  kappacur = 1.0/8.0 * kappacur/mesh->info.hausd;
  kappacur = MG_MIN(kappacur,isqhmin);
  kappacur = MG_MAX(kappacur,isqhmax);

  return kappacur;
}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param pt pointer toward a triangle.
 * \param np1 global index of the first extremity of the edge.
 * \param np2 global index of the second extremity of the edge.
 * \return -1 if no gradation is needed, else index of graded point.
 *
 * Enforces gradation of metric in one extremity of edge \$f[ np1; np2]\$f in
 * tria \a pt with respect to the other, along the direction of the associated
 * support curve first, then along the normal direction.
 *
 * \warning The gradation along the direction normal to the surface is made in
 * an "isotropic way".
 *
 * \remark ALGIANE: a mettre à plat : dans le cas d'une métrique très aniso avec
 * la grande taille quasiment dans la direction de l'arête on se retrouve à
 * modifier la grande taille uniquement (car proche de l'arête) sauf que cette
 * modification n'a quasi pas d'influence sur le calcul de la longueur d'arête.
 */
MMG5_int MMG5_grad2metSurf(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pTria pt, MMG5_int np1,
                           MMG5_int np2)
{
  MMG5_pPoint  p1,p2;
  double       *mm1,*mm2,*nn1,*nn2,ps1,ps2,ux,uy,uz,m1[6],m2[6],n1[3],n2[3],nt[3];
  double       r1[3][3],r2[3][3],t1[2],t2[2],c[3],mtan1[3],mtan2[3],mr1[6],mr2[6];
  double       mtmp[3][3],val,rbasis[3][3];
  double       /*,l1,l2*/l,dd;
  double       lambda[2],vp[2][2],alpha,beta,mu[3];
  int          kmin,idx;
  int8_t       ichg;

  p1 = &mesh->point[np1];
  p2 = &mesh->point[np2];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  uz = p2->c[2] - p1->c[2];

  mm1 = &met->m[6*np1];
  mm2 = &met->m[6*np2];

  if( !MMG5_nortri(mesh,pt,nt) )
    return -1;

  /* Recover normal and metric associated to p1 */
  if( MG_SIN(p1->tag) || (MG_NOM & p1->tag)){
    memcpy(n1,nt,3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }
  else if( MG_GEO & p1->tag ){
    nn1 = &mesh->xpoint[p1->xp].n1[0];
    nn2 = &mesh->xpoint[p1->xp].n2[0];
    ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
    ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

    if( fabs(ps1) < fabs(ps2))
      memcpy(n1,nn2,3*sizeof(double));
    else
      memcpy(n1,nn1,3*sizeof(double));

    /* Note that rbasis is unused in this function */
    if( !MMG5_buildridmet(mesh,met,np1,ux,uy,uz,m1,rbasis) )
      return -1;
  }
  else if( ( MG_REF & p1->tag ) || (MG_BDY & p1->tag) ){
    // if MG_BDY, we are in mmg3d: the normal is stored in the xPoint
    memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }
  else{
    // mmgs only
    memcpy(n1,p1->n,3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }

  /* Recover normal and metric associated to p2 */
  if ( MG_SIN(p2->tag) || (MG_NOM & p2->tag)) {
    memcpy(n2,nt,3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }
  else if ( MG_GEO & p2->tag ) {
    nn1 = &mesh->xpoint[p2->xp].n1[0];
    nn2 = &mesh->xpoint[p2->xp].n2[0];
    ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
    ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

    if( fabs(ps1) < fabs(ps2))
      memcpy(n2,nn2,3*sizeof(double));
    else
      memcpy(n2,nn1,3*sizeof(double));

    /* Note that rbasis is unused in this function */
    if( !MMG5_buildridmet(mesh,met,np2,ux,uy,uz,m2,rbasis) )
      return -1;
  }
  else if( (MG_REF & p2->tag) || (MG_BDY & p2->tag) ){
    // if MG_BDY, we are in mmg3d: the normal is stored in the xPoint
    memcpy(n2,&(mesh->xpoint[p2->xp].n1[0]),3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }
  else{
    //mmgs only
    memcpy(n2,p2->n,3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }

  /* Rotation matrices mapping n1/n2 to e_3 */
  MMG5_rotmatrix(n1,r1);
  MMG5_rotmatrix(n2,r2);

  /* Geodesic length of support curve to edge i */
  l = sqrt(ux*ux+uy*uy+uz*uz);

  /* Characteristic sizes in direction of support curve */
  MMG5_rmtr(r1,m1,mr1);

  mtan1[0] = mr1[0];
  mtan1[1] = mr1[1];
  mtan1[2] = mr1[3];

  c[0] = r1[0][0]*ux + r1[0][1]*uy + r1[0][2]*uz;
  c[1] = r1[1][0]*ux + r1[1][1]*uy + r1[1][2]*uz;

  memcpy(t1,c,2*sizeof(double));
  /* Here we work in the tangent plane (thus in 2d) */
  dd = t1[0]*t1[0] + t1[1]*t1[1];
  if(dd < MMG5_EPSD2)
    return -1;

  dd = 1.0/sqrt(dd);
  t1[0] *= dd;
  t1[1] *= dd;

  /* edge length in metric mtan1: sqrt(t^(t1) * mtan1 * t1). */
  ps1 =  mtan1[0]*t1[0]*t1[0] + 2.0*mtan1[1]*t1[0]*t1[1] + mtan1[2]*t1[1]*t1[1];
  assert ( ps1  > 0. );
  ps1 = sqrt(ps1);

  MMG5_rmtr(r2,m2,mr2);

  mtan2[0] = mr2[0];
  mtan2[1] = mr2[1];
  mtan2[2] = mr2[3];

  c[0] = - ( r2[0][0]*ux + r2[0][1]*uy + r2[0][2]*uz );
  c[1] = - ( r2[1][0]*ux + r2[1][1]*uy + r2[1][2]*uz );
  memcpy(t2,c,2*sizeof(double));

  dd = t2[0]*t2[0] + t2[1]*t2[1];
  if(dd < MMG5_EPSD2)
    return -1;

  dd = 1.0/sqrt(dd);
  t2[0] *= dd;
  t2[1] *= dd;

  /* edge length: sqrt(t^(t2) * mtan2 * t2) */
  ps2 = mtan2[0]*t2[0]*t2[0] + 2.0*mtan2[1]*t2[0]*t2[1] + mtan2[2]*t2[1]*t2[1];
  assert ( ps2  > 0. );
  ps2 = sqrt(ps2);

  /* Metric in p1 has to be changed ( M1 > M2 ) */
  if( ps2 > ps1 ){
    /* compute alpha = h2 + hgrad*l */
    alpha = ps2 /(1.0+mesh->info.hgrad*l*ps2);
    if( ps1 >= alpha - MMG5_EPS ) {
      /* No need to graduate: l_{M2+hgrad*l} < l_{M1} => M2+hgrad*l > M1 */
      return -1;
    }
    MMG5_eigensym(mtan1,lambda,vp);
    /* Project the vector t1 along the main directions of the metric */
    /* Remark: along the third direction mr1 is already diagonal,
     * thus vp[2][.] =( 0 0 1) and vp[.][2] = 0. */
    c[0] = t1[0]*vp[0][0] + t1[1]*vp[0][1] ;
    c[1] = t1[0]*vp[1][0] + t1[1]*vp[1][1] ;

    /* Find index of the maximum value of c: this allow to detect which of the
     * main directions of the metric is closest to our edge direction. We want
     * that our new metric respect the gradation related to the size associated
     * to this main direction (the ichg direction). */
    ichg = 0;
    val  = fabs(c[ichg]);
    for (idx = 1; idx<2; ++idx) {
      if ( fabs(c[idx]) > val ) {
        val = fabs(c[idx]);
        ichg = idx;
      }
    }
    assert(c[ichg]*c[ichg] > MMG5_EPS );
   /* Compute beta coef such as lambda_1 = beta*lambda_1 => h1 = h2 + hgrad*l
    * with h1 = 1/ps1 and h2 = 1/ps2 (see p317 of Charles Dapogny Thesis). */
    beta = (alpha*alpha - ps1*ps1)/(c[ichg]*c[ichg]);

    /* Metric update */
    if( MG_SIN(p1->tag) || (MG_NOM & p1->tag)){
      /* lambda_new = 0.5 lambda_1 + 0.5 beta lambda_1: here we choose to not
       * respect the gradation in order to restric the influence of the singular
       * points. */
      mm1[0] += 0.5*beta;
      mm1[3] += 0.5*beta;
      mm1[5] += 0.5*beta;
    }
    else if( p1->tag & MG_GEO ){
      /* lambda[ichg] is the metric eigenvalue associated to the main metric
       * direction closest to our edge direction. Find were is stored this
       * eigenvalue in our special storage of ridge metric (mm-lambda = 0) and
       * update it. */
      c[0] = fabs(mm1[0]-lambda[ichg]);
      c[1] = fabs(mm1[1]-lambda[ichg]);
      c[2] = fabs(mm1[2]-lambda[ichg]);

      /* Find index of the minimum value of c */
      kmin = 0;
      val = fabs(c[kmin]);
      for (idx = 1; idx<3; ++idx) {
        if ( fabs(c[idx]) < val ) {
          val = fabs(c[idx]);
          kmin = idx;
        }
      }
      mm1[kmin] += beta;
    }
    else{
      /* Update the metric eigenvalue associated to the main metric direction
       * which is closest to our edge direction (because this is the one that is
       * the more influent on our edge length). */
      mu[0] = lambda[0];
      mu[1] = lambda[1];
      mu[2] = mr1[5];

      mu[ichg] += beta;

      mtan1[0] = mu[0]*vp[0][0]*vp[0][0] + mu[1]*vp[1][0]*vp[1][0];
      mtan1[1] = mu[0]*vp[0][0]*vp[0][1] + mu[1]*vp[1][0]*vp[1][1];
      mtan1[2] = mu[0]*vp[0][1]*vp[0][1] + mu[1]*vp[1][1]*vp[1][1];

      /* Return in initial basis */
      /* Because of the rotation, we know that:
       * mr.[2] = mr.[4]= 0 */
      mtmp[0][0] = mtan1[0]*r1[0][0] + mtan1[1]*r1[1][0];
      mtmp[0][1] = mtan1[0]*r1[0][1] + mtan1[1]*r1[1][1];
      mtmp[0][2] = mtan1[0]*r1[0][2] + mtan1[1]*r1[1][2];

      mtmp[1][0] = mtan1[1]*r1[0][0] + mtan1[2]*r1[1][0];
      mtmp[1][1] = mtan1[1]*r1[0][1] + mtan1[2]*r1[1][1];
      mtmp[1][2] = mtan1[1]*r1[0][2] + mtan1[2]*r1[1][2];

      mtmp[2][0] = mr1[5]*r1[2][0];
      mtmp[2][1] = mr1[5]*r1[2][1];
      mtmp[2][2] = mr1[5]*r1[2][2];

      m1[0] = r1[0][0]*mtmp[0][0] + r1[1][0]*mtmp[1][0] + r1[2][0]*mtmp[2][0];
      m1[1] = r1[0][0]*mtmp[0][1] + r1[1][0]*mtmp[1][1] + r1[2][0]*mtmp[2][1];
      m1[2] = r1[0][0]*mtmp[0][2] + r1[1][0]*mtmp[1][2] + r1[2][0]*mtmp[2][2];

      m1[3] = r1[0][1]*mtmp[0][1] + r1[1][1]*mtmp[1][1] + r1[2][1]*mtmp[2][1];
      m1[4] = r1[0][1]*mtmp[0][2] + r1[1][1]*mtmp[1][2] + r1[2][1]*mtmp[2][2];

      m1[5] = r1[0][2]*mtmp[0][2] + r1[1][2]*mtmp[1][2] + r1[2][2]*mtmp[2][2];

      memcpy(mm1,m1,6*sizeof(double));
    }
    return np1;
  }
  /* Metric in p2 has to be changed ( M2 > M1 ) */
  else{
    alpha = ps1 /(1.0+mesh->info.hgrad*l*ps1);
    if( ps2 >= alpha - MMG5_EPS) {
      /* No need to graduate: l_{M1+hgrad*l} < l_{M2} => M1+hgrad*l > M2 */
      return -1;
    }
    MMG5_eigensym(mtan2,lambda,vp);

    c[0] = t2[0]*vp[0][0] + t2[1]*vp[0][1] ;
    c[1] = t2[0]*vp[1][0] + t2[1]*vp[1][1] ;

    /* Detect which of the main directions of the metric is closest to our edge
     * direction. */
    ichg = 0;
    val  = fabs(c[ichg]);
    for (idx = 1; idx<2; ++idx) {
      if ( fabs(c[idx]) > val ) {
        val = fabs(c[idx]);
        ichg = idx;
      }
    }
    assert(c[ichg]*c[ichg] > MMG5_EPS );

   /* Compute beta coef such as lambda_2new = beta*lambda_2 => h2 = h1 + hgrad*l
    * (see p317 of Charles Dapogny Thesis). */
    beta = (alpha*alpha - ps2*ps2)/(c[ichg]*c[ichg]);

    /* Metric update: update the metric eigenvalue associated to the main metric
       * direction which is closest to our edge direction (because this is the
       * one that is the more influent on our edge length). */
    if( MG_SIN(p2->tag) || (MG_NOM & p2->tag)){
      /* lambda_new = 0.5 lambda_1 + 0.5 beta lambda_1: here we choose to not
       * respect the gradation in order to restric the influence of the singular
       * points. */
      mm2[0] += 0.5*beta;
      mm2[3] += 0.5*beta;
      mm2[5] += 0.5*beta;
    }
    else if( p2->tag & MG_GEO ){
      c[0] = fabs(mm2[0]-lambda[ichg]);
      c[1] = fabs(mm2[1]-lambda[ichg]);
      c[2] = fabs(mm2[2]-lambda[ichg]);

      kmin = 0;
      val = fabs(c[kmin]);
      for (idx = 1; idx<3; ++idx) {
        if ( fabs(c[idx]) < val ) {
          val = fabs(c[idx]);
          kmin = idx;
        }
      }
      mm2[kmin] += beta;
    }
    else{
      mu[0] = lambda[0];
      mu[1] = lambda[1];
      mu[2] = mr2[5];

      mu[ichg] += beta;

      mtan2[0] = mu[0]*vp[0][0]*vp[0][0] + mu[1]*vp[1][0]*vp[1][0];
      mtan2[1] = mu[0]*vp[0][0]*vp[0][1] + mu[1]*vp[1][0]*vp[1][1];
      mtan2[2] = mu[0]*vp[0][1]*vp[0][1] + mu[1]*vp[1][1]*vp[1][1];

      /* Return in initial basis */
      mtmp[0][0] = mtan2[0]*r2[0][0] + mtan2[1]*r2[1][0];
      mtmp[0][1] = mtan2[0]*r2[0][1] + mtan2[1]*r2[1][1];
      mtmp[0][2] = mtan2[0]*r2[0][2] + mtan2[1]*r2[1][2];

      mtmp[1][0] = mtan2[1]*r2[0][0] + mtan2[2]*r2[1][0];
      mtmp[1][1] = mtan2[1]*r2[0][1] + mtan2[2]*r2[1][1];
      mtmp[1][2] = mtan2[1]*r2[0][2] + mtan2[2]*r2[1][2];

      mtmp[2][0] =  mr2[5]*r2[2][0];
      mtmp[2][1] =  mr2[5]*r2[2][1];
      mtmp[2][2] =  mr2[5]*r2[2][2];

      m2[0] = r2[0][0]*mtmp[0][0] + r2[1][0]*mtmp[1][0] + r2[2][0]*mtmp[2][0];
      m2[1] = r2[0][0]*mtmp[0][1] + r2[1][0]*mtmp[1][1] + r2[2][0]*mtmp[2][1];
      m2[2] = r2[0][0]*mtmp[0][2] + r2[1][0]*mtmp[1][2] + r2[2][0]*mtmp[2][2];

      m2[3] = r2[0][1]*mtmp[0][1] + r2[1][1]*mtmp[1][1] + r2[2][1]*mtmp[2][1];
      m2[4] = r2[0][1]*mtmp[0][2] + r2[1][1]*mtmp[1][2] + r2[2][1]*mtmp[2][2];

      m2[5] = r2[0][2]*mtmp[0][2] + r2[1][2]*mtmp[1][2] + r2[2][2]*mtmp[2][2];

      memcpy(mm2,m2,6*sizeof(double));
    }
    return np2;
  }
}

/**
 * \param mesh pointer toward the mesh
 * \param m first matrix
 * \param n second matrix
 * \param dm eigenvalues of m in the coreduction basis (to fill)
 * \param dn eigenvalues of n in the coreduction basis (to fill)
 * \param vp coreduction basis (to fill)
 *
 * \return 0 if fail 1 otherwise.
 *
 * Perform simultaneous reduction of matrices \a m and \a n.
 *
 */
int MMG5_simred2d(MMG5_pMesh mesh,double *m,double *n,double dm[2],
                  double dn[2],double vp[2][2] ) {

  double         det,lambda[2],imn[4];
  int            order;
  static int8_t  mmgWarn0=0;

  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < MMG5_EPS*MMG5_EPS ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 null metric det : %E \n",
              __func__,det);
    }
    return 0;
  }
  det = 1.0 / det;

  imn[0] = det * ( m[2]*n[0] - m[1]*n[1]);
  imn[1] = det * ( m[2]*n[1] - m[1]*n[2]);
  imn[2] = det * (-m[1]*n[0] + m[0]*n[1]);
  imn[3] = det * (-m[1]*n[1] + m[0]*n[2]);

  /* Find eigenvalues of imn */
  order = MMG5_eigenv2d(0,imn,lambda,vp);

  if ( !order ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 failing"
              " simultaneous reduction.\n",__func__);
    }
    return 0;
  }

  /* First case : matrices m and n are homothetic: n = lambda0*m */
  if ( order == 2 ) {

    /* Subcase where m is diaonal */
    if ( fabs(m[1]) < MMG5_EPS ) {
      dm[0]   = m[0];
      dm[1]   = m[2];
      vp[0][0] = 1;
      vp[0][1] = 0;
      vp[1][0] = 0;
      vp[1][1] = 1;
    }
    /* Subcase where m is not diagonal; dd,trimn,... are reused */
    else
      MMG5_eigensym(m,dm,vp);

    /* Eigenvalues of metric n */
    dn[0] = lambda[0]*dm[0];
    dn[1] = lambda[0]*dm[1];

  }
  /* Second case: both eigenvalues of imn are distinct ; theory says qf associated to m and n
   are diagonalizable in basis (vp[0], vp[1]) - the coreduction basis */
  else if( order == 1 ) {

    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp[0][0]*vp[0][0] + 2.0*m[1]*vp[0][0]*vp[0][1] + m[2]*vp[0][1]*vp[0][1];
    dm[1] = m[0]*vp[1][0]*vp[1][0] + 2.0*m[1]*vp[1][0]*vp[1][1] + m[2]*vp[1][1]*vp[1][1];
    dn[0] = n[0]*vp[0][0]*vp[0][0] + 2.0*n[1]*vp[0][0]*vp[0][1] + n[2]*vp[0][1]*vp[0][1];
    dn[1] = n[0]*vp[1][0]*vp[1][0] + 2.0*n[1]*vp[1][0]*vp[1][1] + n[2]*vp[1][1]*vp[1][1];
  }

  assert ( dm[0] >= MMG5_EPSD2 &&  dm[1] >= MMG5_EPSD2 && "positive eigenvalue" );
  assert ( dn[0] >= MMG5_EPSD2 &&  dn[1] >= MMG5_EPSD2 && "positive eigenvalue" );

  if ( dm[0] < MMG5_EPSOK || dn[0] < MMG5_EPSOK ) { return 0; }
  if ( dm[1] < MMG5_EPSOK || dn[1] < MMG5_EPSOK ) { return 0; }

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param m first matrix
 * \param n second matrix
 * \param dm eigenvalues of m in the coreduction basis (to fill)
 * \param dn eigenvalues of n in the coreduction basis (to fill)
 * \param vp coreduction basis (to fill)
 *
 * \return 0 if fail 1 otherwise.
 *
 * Perform simultaneous reduction of matrices \a m and \a n.
 *
 */
int MMG5_simred3d(MMG5_pMesh mesh,double *m,double *n,double dm[3],
                  double dn[3],double vp[3][3] ) {

  double        lambda[3],im[6],imn[9];
  int           order;
  static int8_t mmgWarn0=0;

  /* Compute imn = M^{-1}N */
  if ( !MMG5_invmat ( m,im ) ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: unable to invert the matrix.\n",__func__);
    }
    return 0;
  }

  MMG5_mn(im,n,imn);

  /* Find eigenvalues of imn */
  order = MMG5_eigenv3d(0,imn,lambda,vp);

  if ( !order ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 failing"
              " simultaneous reduction.\n",__func__);
    }
    return 0;
  }

  if ( order == 3 ) {
    /* First case : matrices m and n are homothetic: n = lambda0*m */
    if ( (fabs(m[1]) < MMG5_EPS && fabs(m[2]) < MMG5_EPS
          && fabs(m[4]) < MMG5_EPS) ) {
      /* Subcase where m is diaonal */
        dm[0]   = m[0];
        dm[1]   = m[3];
        dm[2]   = m[5];
        vp[0][0] = 1;
        vp[0][1] = 0;
        vp[0][2] = 0;
        vp[1][0] = 0;
        vp[1][1] = 1;
        vp[1][2] = 0;
        vp[2][0] = 0;
        vp[2][1] = 0;
        vp[2][2] = 1;
    }
    else {
      /* Subcase where m is not diagonal; dd,trimn,... are reused */
      MMG5_eigenv3d(1,m,dm,vp);
    }
    /* Eigenvalues of metric n */
    dn[0] = lambda[0]*dm[0];
    dn[1] = lambda[0]*dm[1];
    dn[2] = lambda[0]*dm[2];
  }
  else if( order == 2 ) {
    /* Second case: two eigenvalues of imn are coincident (first two entries of
     * the lambda array) and one is distinct (last entry).
     * Simultaneous reduction gives a block diagonalization. The 2x2 blocks are
     * homothetic and can be diagonalized through the eigenvectors of one of the
     * two blocks. */
    double mred[6],nred[6];
    /* project matrices on the coreduction basis: they should have the
     * block-diagonal form [ m0, m1, 0, m3, 0, m5 ] */
    MMG5_rmtr(vp,m,mred);
    MMG5_rmtr(vp,n,nred);
    /* compute projections on the last eigenvector (that with multiplicity 1) */
    dm[2] = mred[5];
    dn[2] = nred[5];
    /* re-arrange matrices so that the first three entries describe the
     * 2x2 blocks to be diagonalized (the two blocks are homothetic) */
    mred[2] = mred[3];
    nred[2] = nred[3];
    /* diagonalization of the first 2x2 block */
    if( fabs(mred[1]) < MMG5_EPS ) {
     /* first case: the blocks are diagonal, basis vp is unchanged */
      dm[0] = mred[0];
      dm[1] = mred[2];
    } else {
      /* second case: the blocks are not diagonal */
      double wp[2][2],up[2][3];
      int8_t i,j,k;
      MMG5_eigensym(mred,dm,wp);
      /* change the basis vp (vp[2] is unchanged) */
      for( j = 0; j < 2; j++ ) {
        for( i = 0; i < 3; i++ ) {
          up[j][i] = 0.;
          for( k = 0; k < 2; k++ ) {
            up[j][i] += vp[k][i]*wp[j][k];
          }
        }
      }
      for( j = 0; j < 2; j++ ) {
        for( i = 0; i < 3; i++ ) {
          vp[j][i] = up[j][i];
        }
      }
    }
    /* homothetic diagonalization of the second 2x2 block */
    dn[0] = lambda[0]*dm[0];
    dn[1] = lambda[0]*dm[1];
  }
  else {
    /* Third case: eigenvalues of imn are distinct ; theory says qf associated
       to m and n are diagonalizable in basis (vp[0], vp[1], vp[2]) - the
       coreduction basis */
    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp[0][0]*vp[0][0] + 2.0*m[1]*vp[0][0]*vp[0][1] + 2.0*m[2]*vp[0][0]*vp[0][2]
          + m[3]*vp[0][1]*vp[0][1] + 2.0*m[4]*vp[0][1]*vp[0][2]     + m[5]*vp[0][2]*vp[0][2];
    dm[1] = m[0]*vp[1][0]*vp[1][0] + 2.0*m[1]*vp[1][0]*vp[1][1] + 2.0*m[2]*vp[1][0]*vp[1][2]
          + m[3]*vp[1][1]*vp[1][1] + 2.0*m[4]*vp[1][1]*vp[1][2]     + m[5]*vp[1][2]*vp[1][2];
    dm[2] = m[0]*vp[2][0]*vp[2][0] + 2.0*m[1]*vp[2][0]*vp[2][1] + 2.0*m[2]*vp[2][0]*vp[2][2]
          + m[3]*vp[2][1]*vp[2][1] + 2.0*m[4]*vp[2][1]*vp[2][2]     + m[5]*vp[2][2]*vp[2][2];

    dn[0] = n[0]*vp[0][0]*vp[0][0] + 2.0*n[1]*vp[0][0]*vp[0][1] + 2.0*n[2]*vp[0][0]*vp[0][2]
          + n[3]*vp[0][1]*vp[0][1] + 2.0*n[4]*vp[0][1]*vp[0][2]     + n[5]*vp[0][2]*vp[0][2];
    dn[1] = n[0]*vp[1][0]*vp[1][0] + 2.0*n[1]*vp[1][0]*vp[1][1] + 2.0*n[2]*vp[1][0]*vp[1][2]
          + n[3]*vp[1][1]*vp[1][1] + 2.0*n[4]*vp[1][1]*vp[1][2]     + n[5]*vp[1][2]*vp[1][2];
    dn[2] = n[0]*vp[2][0]*vp[2][0] + 2.0*n[1]*vp[2][0]*vp[2][1] + 2.0*n[2]*vp[2][0]*vp[2][2]
          + n[3]*vp[2][1]*vp[2][1] + 2.0*n[4]*vp[2][1]*vp[2][2]     + n[5]*vp[2][2]*vp[2][2];
  }

  assert ( dm[0] >= MMG5_EPSD2 && dm[1] >= MMG5_EPSD2 && dm[2] >= MMG5_EPSD2 && "positive eigenvalue" );
  assert ( dn[0] >= MMG5_EPSD2 && dn[1] >= MMG5_EPSD2 && dn[2] >= MMG5_EPSD2 && "positive eigenvalue" );

  if ( dm[0] < MMG5_EPSOK || dn[0] < MMG5_EPSOK ) { return 0; }
  if ( dm[1] < MMG5_EPSOK || dn[1] < MMG5_EPSOK ) { return 0; }
  if ( dm[2] < MMG5_EPSOK || dn[2] < MMG5_EPSOK ) { return 0; }

  return 1;
}

/**
 * \param dim square matrix size
 * \param dm diagonal values array
 * \param dn diagonal values array
 * \param vp basis vectors array
 * \param swap swap array
 * \param perm permutation array
 *
 * Sort and permute diagonal values (and basis vectors) in increasing order
 * with respect to the first matrix.
 *
 */
void MMG5_sort_simred( int8_t dim,double *dm,double *dn,double *vp,
                       double *swap,int8_t *perm ) {
  MMG5_nsort(dim,dm,perm);
  MMG5_nperm(dim,0,1,dm,swap,perm);
  MMG5_nperm(dim,0,1,dn,swap,perm);
  for( int8_t i = 0; i < dim; i++ )
    MMG5_nperm(dim,i,dim,vp,swap,perm);
}

/**
 * \param mesh pointer toward the mesh structure
 * \param mex first symmetric test matrix
 * \param nex second symmetric test matrix
 * \param dm diagonalization of the first matrix on the reduction basis
 * \param dn diagonalization of the second matrix on the reduction basis
 * \param vp simultaneous reduction basis (stored by columns)
 * \return 1 if success, 0 if fail
 *
 * For a couple of 2x2 symmetric matrices, Test:
 * - the computation of the simultaneous reduction values of the matrices;
 * - the computation of the simultaneous reduction basis vectors..
 *
 */
int MMG5_test_simred2d(MMG5_pMesh mesh,double *mex,double *nex,double *dmex,double *dnex,double vpex[][2]) {
  double dmnum[2],dnnum[2],vpnum[2][2],ivpnum[2][2],mnum[3],nnum[3]; /* Numerical quantities */
  double swap[2],maxerr,err;
  int8_t perm[2]; /* permutation array */

  /** Compute simultaneous reduction */
  if( !MMG5_simred2d(mesh,mex,nex,dmnum,dnnum,vpnum ) )
    return 0;

  /* Naively sort eigenpairs in increasing order */
  MMG5_sort_simred(2,dmnum,dnnum,(double *)vpnum,swap,perm);

  /* Check diagonal values error in norm inf */
  maxerr = MMG5_test_mat_error(2,(double *)dmex,(double *)dmnum);
  if( maxerr > 100*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error first matrix coreduction values: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }
  maxerr = MMG5_test_mat_error(2,(double *)dnex,(double *)dnnum);
  if( maxerr > 1000*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error second matrix coreduction values: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /* Check eigenvectors error through scalar product */
  maxerr = 0.;
  for( int8_t i = 0; i < 2; i++ ) {
    err = 0.;
    for( int8_t j = 0; j < 2; j++ )
      err += vpex[i][j] * vpnum[i][j];
    err = 1.-fabs(err);
    maxerr = MG_MAX(maxerr,err);
  }
  if( maxerr > MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix coreduction vectors: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /** Recompose matrices from diagonal values */
  if( !MMG5_invmat22(vpnum,ivpnum) )
    return 0;
  MMG5_simredmat(2,mnum,dmnum,(double *)ivpnum);
  MMG5_simredmat(2,nnum,dnnum,(double *)ivpnum);

  /* Check matrices in norm inf */
  maxerr = MMG5_test_mat_error(3,mex,mnum);
  if( maxerr > 1.e2*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error first matrix coreduction recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }
  maxerr = MMG5_test_mat_error(3,nex,nnum);
  if( maxerr > 1.e4*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error second matrix coreduction recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param mex first symmetric test matrix
 * \param nex second symmetric test matrix
 * \param dm diagonalization of the first matrix on the reduction basis
 * \param dn diagonalization of the second matrix on the reduction basis
 * \param vp simultaneous reduction basis (stored by columns)
 * \return 1 if success, 0 if fail
 *
 * For a couple of 3x3 symmetric matrices, Test:
 * - the computation of the simultaneous reduction values of the matrices;
 * - the computation of the simultaneous reduction basis vectors..
 *
 */
int MMG5_test_simred3d(MMG5_pMesh mesh,double *mex,double *nex,double *dmex,double *dnex,double vpex[][3]) {
  double dmnum[3],dnnum[3],vpnum[3][3],ivpnum[3][3],mnum[6],nnum[6]; /* Numerical quantities */
  double swap[3],maxerr,err;
  int8_t perm[3]; /* permutation array */

  /** Compute simultaneous reduction */
  if( !MMG5_simred3d(mesh,mex,nex,dmnum,dnnum,vpnum ) )
    return 0;

  /* Naively sort eigenpairs in increasing order */
  MMG5_sort_simred(3,dmnum,dnnum,(double *)vpnum,swap,perm);

  /* Check diagonal values error in norm inf */
  maxerr = MMG5_test_mat_error(3,(double *)dmex,(double *)dmnum);
  if( maxerr > 100*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error first matrix coreduction values: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }
  maxerr = MMG5_test_mat_error(3,(double *)dnex,(double *)dnnum);
  if( maxerr > 1000*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error second matrix coreduction values: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /* Check eigenvectors error through scalar product */
  maxerr = 0.;
  for( int8_t i = 0; i < 3; i++ ) {
    err = 0.;
    for( int8_t j = 0; j < 3; j++ )
      err += vpex[i][j] * vpnum[i][j];
    err = 1.-fabs(err);
    maxerr = MG_MAX(maxerr,err);
  }
  if( maxerr > MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix coreduction vectors: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /** Recompose matrices from diagonal values */
  if( !MMG5_invmat33(vpnum,ivpnum) )
    return 0;
  MMG5_simredmat(3,mnum,dmnum,(double *)ivpnum);
  MMG5_simredmat(3,nnum,dnnum,(double *)ivpnum);

  /* Check matrices in norm inf */
  maxerr = MMG5_test_mat_error(6,mex,mnum);
  if( maxerr > 1.e2*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error first matrix coreduction recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }
  maxerr = MMG5_test_mat_error(6,nex,nnum);
  if( maxerr > 1.e3*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error second matrix coreduction recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param m first matrix
 * \param n second matrix
 * \param dm eigenvalues of m in the coreduction basis
 * \param dn eigenvalues of n in the coreduction basis
 * \param vp coreduction basis
 * \param ier flag of the updated sizes: (ier & 1) if we dm has been modified, (ier & 2) if dn has been modified.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored in
 * columns in 2D.
 *
 */
inline
int MMG5_updatemet2d_ani(double *m,double *n,double dm[2],double dn[2],
                         double vp[2][2],int8_t ier ) {
  double det,ip[4];

  det = vp[0][0]*vp[1][1] - vp[0][1]*vp[1][0];
  if ( fabs(det) < MMG5_EPS )  return 0;
  det = 1.0 / det;

  ip[0] =  vp[1][1]*det;
  ip[1] = -vp[1][0]*det;
  ip[2] = -vp[0][1]*det;
  ip[3] =  vp[0][0]*det;

  if ( ier & 1 ) {
    m[0] = dm[0]*ip[0]*ip[0] + dm[1]*ip[2]*ip[2];
    m[1] = dm[0]*ip[0]*ip[1] + dm[1]*ip[2]*ip[3];
    m[2] = dm[0]*ip[1]*ip[1] + dm[1]*ip[3]*ip[3];
  }
  if ( ier & 2 ) {
    n[0] = dn[0]*ip[0]*ip[0] + dn[1]*ip[2]*ip[2];
    n[1] = dn[0]*ip[0]*ip[1] + dn[1]*ip[2]*ip[3];
    n[2] = dn[0]*ip[1]*ip[1] + dn[1]*ip[3]*ip[3];
  }
  return 1;
}

/**
 * \param m first matrix
 * \param n second matrix
 * \param dm eigenvalues of m in the coreduction basis
 * \param dn eigenvalues of n in the coreduction basis
 * \param vp coreduction basis
 * \param ier flag of the updated sizes: (ier & 1) if we dm has been modified, (ier & 2) if dn has been modified.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored in
 * columns in 3D.
 *
 */
inline
int MMG5_updatemet3d_ani(double *m,double *n,double dm[3],double dn[3],
                         double vp[3][3],int8_t ier ) {
  double ivp[3][3];

  if( !MMG5_invmat33(vp,ivp) )
    return 0;

  if ( ier & 1 ) {
    MMG5_simredmat(3,m,dm,(double *)ivp);
  }
  if ( ier & 2 ) {
    MMG5_simredmat(3,n,dn,(double *)ivp);
  }
  return 1;
}

/**
 *
 * Test Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored
 * in columns in 2D.
 *
 */
int MMG5_test_updatemet2d_ani() {
  double mex[3] = { 508., -504,  502.}; /* Test matrix 1 */
  double nex[3] = {4020.,-2020.,1020.}; /* Test matrix 2 */
  double dm[2] = {  1., 100. }; /* Exact cobasis projection 1 */
  double dn[2] = {500.,   4. }; /* Exact cobasis projection 2 */
  double vp[2][2] = {{ 1./sqrt(2.),1./sqrt(2.)},
                     {1./sqrt(5.),2./sqrt(5.)}}; /* Exact cobasis vectors */
  double mnum[3],nnum[3],maxerr; /* Numerical quantities */

  /** Recompose matrices from exact simultaneous diagonalization */
  if( !MMG5_updatemet2d_ani(mnum,nnum,dm,dn,vp,3) )
    return 0;

  /* Check values error in norm inf */
  maxerr = MMG5_test_mat_error(3,(double *)mex,(double *)mnum);
  if( maxerr > 1.e4*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error first matrix recomposition from simultaneous reduction: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }
  maxerr = MMG5_test_mat_error(3,(double *)nex,(double *)nnum);
  if( maxerr > 1.e4*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error second matrix recomposition from simultaneous reduction: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 *
 * Test Update of the metrics = tP^-1 diag(d0,d1)P^-1, P = (vp[0], vp[1]) stored
 * in columns in 3D.
 *
 */
int MMG5_test_updatemet3d_ani() {
  double mex[6] = {111./2.,-109./2.,  89./2.,111./2.,-91./2.,111./2.}; /* Test matrix 1 */
  double nex[6] = {409./2.,-393./2.,-407./2.,409./2.,391./2.,409./2.}; /* Test matrix 2 */
  double dm[3] = {1., 10.,100.}; /* Exact cobasis projection 1 */
  double dn[3] = {8.,400.,  1.}; /* Exact cobasis projection 2 */
  double vp[3][3] = {{1./sqrt(2.),1./sqrt(2.),0.},
                     {0.,         1./sqrt(2.),1./sqrt(2.)},
                     {1./sqrt(2.),         0.,1./sqrt(2.)}}; /* Exact cobasis vectors */
  double mnum[6],nnum[6],maxerr; /* Numerical quantities */

  /** Recompose matrices from exact simultaneous diagonalization */
  if( !MMG5_updatemet3d_ani(mnum,nnum,dm,dn,vp,3) )
    return 0;

  /* Check values error in norm inf */
  maxerr = MMG5_test_mat_error(6,(double *)mex,(double *)mnum);
  if( maxerr > 100*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error first matrix recomposition from simultaneous reduction: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }
  maxerr = MMG5_test_mat_error(6,(double *)nex,(double *)nnum);
  if( maxerr > 100*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error second matrix recomposition from simultaneous reduction: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param dm eigenvalues of the first matrix (not modified)
 * \param dn eigenvalues of the second matrix (modified)
 * \param difsiz maximal size gap authorized by the gradation.
 * \param dir direction in which the sizes are graded.
 * \param ier 2 if dn has been updated, 0 otherwise.
 *
 *  Gradation of size dn = 1/sqrt(eigenv of the tensor) for required points in
 *  the \a idir direction.
 *
 */
void MMG5_gradEigenvreq(double *dm,double *dn,double difsiz,int8_t dir,int8_t *ier) {
  double hm,hn;

  hm = 1.0 / sqrt(dm[dir]);
  hn = 1.0 / sqrt(dn[dir]);

  if ( hn > hm + difsiz + MMG5_EPSOK ) {
    /* Decrease the size in \a ipslave */
    hn = hm+difsiz;
    dn[dir] = 1.0 / (hn*hn);
    (*ier) = 2;
  }
  else if ( hn + MMG5_EPSOK < hm - difsiz ) {
    /* Increase the size in \a ipslave */
    hn = hm-difsiz;
    dn[dir] = 1.0 / (hn*hn);
    (*ier) = 2;
  }
}

/**
 * \param n  matrix to update
 * \param dn eigenvalues of n in the coreduction basis
 * \param vp coreduction basis
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update of the metric n = tP^-1 diag(dn0,dn1)P^-1, P = (vp[0], vp[1]) stored in
 * columns
 *
 */
int MMG5_updatemetreq_ani(double *n,double dn[2],double vp[2][2]) {
  double det,ip[4];

  det = vp[0][0]*vp[1][1] - vp[0][1]*vp[1][0];
  if ( fabs(det) < MMG5_EPS )  return 0;
  det = 1.0 / det;

  ip[0] =  vp[1][1]*det;
  ip[1] = -vp[1][0]*det;
  ip[2] = -vp[0][1]*det;
  ip[3] =  vp[0][0]*det;

  n[0] = dn[0]*ip[0]*ip[0] + dn[1]*ip[2]*ip[2];
  n[1] = dn[0]*ip[0]*ip[1] + dn[1]*ip[2]*ip[3];
  n[2] = dn[0]*ip[1]*ip[1] + dn[1]*ip[3]*ip[3];

  return 1;
}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param pt pointer toward the processed triangle.
 * \param npmaster edge extremity that cannot be modified
 * \param npslave edge extremity to modify to respect the gradation.
 *
 * \return 0 if no gradation is needed, 1 otherwise.
 *
 * Enforces gradation of metric of the extremity ±a npslave of edge \$f[
 * npmaster; npslave]\$f in tria \a pt with respect to the other, along the
 * direction of the associated support curve first, then along the normal
 * direction.
 *
 * \warning The gradation along the direction normal to the surface is made in
 * an "isotropic way".
 *
 */
int MMG5_grad2metSurfreq(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pTria pt, MMG5_int npmaster,
                         MMG5_int npslave)
{

  MMG5_pPoint  p1,p2;
  double       *mm1,*mm2,*nn1,*nn2,ps1,ps2,ux,uy,uz,m1[6],m2[6],n1[3],n2[3],nt[3];
  double       r1[3][3],r2[3][3],mtan1[3],mtan2[3],mr1[6],mr2[6];
  double       mtmp[3][3],rbasis1[3][3],rbasis2[3][3];
  double       l,difsiz,rmet3D[6];
  double       lambda[2],vp[2][2],beta,mu[3];
  int          cfg_m2;
  int8_t       ier;

  p1 = &mesh->point[npmaster];
  p2 = &mesh->point[npslave];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  uz = p2->c[2] - p1->c[2];

  mm1 = &met->m[6*npmaster];
  mm2 = &met->m[6*npslave];

  cfg_m2 = 0;
  ier = 0;

  if( !MMG5_nortri(mesh,pt,nt) )
    return 0;

  /* Recover normal and metric associated to p1 */
  if( MG_SIN(p1->tag) || (MG_NOM & p1->tag)){
    memcpy(n1,nt,3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }
  else if( MG_GEO & p1->tag ){
    nn1 = &mesh->xpoint[p1->xp].n1[0];
    nn2 = &mesh->xpoint[p1->xp].n2[0];
    ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
    ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

    if( fabs(ps1) < fabs(ps2))
      memcpy(n1,nn2,3*sizeof(double));
    else
      memcpy(n1,nn1,3*sizeof(double));

    if( !MMG5_buildridmetnor(mesh,met,npmaster,nt,m1,rbasis1) ) { return 0; }
  }
  else if( ( MG_REF & p1->tag ) || (MG_BDY & p1->tag) ){
    // if MG_BDY, we are in mmg3d: the normal is stored in the xPoint
    memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }
  else{
    // mmgs only
    memcpy(n1,p1->n,3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }

  /* Recover normal and metric associated to p2 */
  if ( MG_SIN(p2->tag) || (MG_NOM & p2->tag)) {
    memcpy(n2,nt,3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }
  else if ( MG_GEO & p2->tag ) {
    nn1 = &mesh->xpoint[p2->xp].n1[0];
    nn2 = &mesh->xpoint[p2->xp].n2[0];
    ps1 = nt[0]*nn1[0] + nt[1]*nn1[1] + nt[2]*nn1[2];
    ps2 = nt[0]*nn2[0] + nt[1]*nn2[1] + nt[2]*nn2[2];

    if( fabs(ps1) < fabs(ps2))
      memcpy(n2,nn2,3*sizeof(double));
    else
      memcpy(n2,nn1,3*sizeof(double));

    cfg_m2 = MMG5_buildridmetnor(mesh,met,npslave,nt,m2,rbasis2);
    if( !cfg_m2 ) { return 0; }
  }
  else if( (MG_REF & p2->tag) || (MG_BDY & p2->tag) ){
    // if MG_BDY, we are in mmg3d: the normal is stored in the xPoint
    memcpy(n2,&(mesh->xpoint[p2->xp].n1[0]),3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }
  else{
   //mmgs Only
    memcpy(n2,p2->n,3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }

  /* Rotation matrices mapping n1/n2 to e_3 */
  MMG5_rotmatrix(n1,r1);
  MMG5_rotmatrix(n2,r2);

  /* Geodesic length of support curve to edge i */
  l = sqrt(ux*ux+uy*uy+uz*uz);

  /* Characteristic sizes in direction of support curve */
  MMG5_rmtr(r1,m1,mr1);

  mtan1[0] = mr1[0];
  mtan1[1] = mr1[1];
  mtan1[2] = mr1[3];

  MMG5_rmtr(r2,m2,mr2);

  mtan2[0] = mr2[0];
  mtan2[1] = mr2[1];
  mtan2[2] = mr2[3];

  difsiz = mesh->info.hgradreq*l;

  /* Simultaneous reduction of mtan1 and mtan2 */
  if ( !MMG5_simred2d(mesh,mtan1,mtan2,lambda,mu,vp) ) {
    return 0;
  }

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the first direction */
  MMG5_gradEigenvreq(lambda,mu,difsiz,0,&ier);

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the second direction */
  MMG5_gradEigenvreq(lambda,mu,difsiz,1,&ier);

  if ( !ier ) {
    return 0;
  }

  /* Metric update using the simultaneous reduction technique */
  if( MG_SIN(p2->tag) || (MG_NOM & p2->tag)){
    /* We choose to not respect the gradation in order to restrict the influence
     * of the singular points. Thus:
     * lambda_new = = 0.5 lambda_1 + 0.5 lambda_new = lambda_1 + 0.5 beta.
     * with beta the smallest variation of the eigenvalues (lambda_new-lambda_1). */

    /* This point can have an anisotropic metric if a user-provided metric is
     * found. So, compute the eigenvalues. */
    double ll[3],rr[3][3],llmin;
    int i;
    if( !MMG5_eigenv3d(1,mm2,ll, rr) ) {
      return 0;
    }
    llmin = DBL_MAX;
    for( i = 0; i < 3; i++ )
      if( ll[i] < llmin )
        llmin = ll[i];


    beta = mu[0] - llmin;

    if ( fabs(beta) < fabs(llmin-mu[1]) ) {
      beta = mu[1] - llmin;
    }
    ll[0] += 0.5*beta;
    ll[1] += 0.5*beta;
    ll[2] += 0.5*beta;
    assert ( ll[0]>0. && ll[1]>0. && ll[2]>0. );
    mm2[0] = ll[0]*rr[0][0]*rr[0][0] + ll[1]*rr[1][0]*rr[1][0] + ll[2]*rr[2][0]*rr[2][0];
    mm2[1] = ll[0]*rr[0][0]*rr[0][1] + ll[1]*rr[1][0]*rr[1][1] + ll[2]*rr[2][0]*rr[2][1];
    mm2[2] = ll[0]*rr[0][0]*rr[0][2] + ll[1]*rr[1][0]*rr[1][2] + ll[2]*rr[2][0]*rr[2][2];
    mm2[3] = ll[0]*rr[0][1]*rr[0][1] + ll[1]*rr[1][1]*rr[1][1] + ll[2]*rr[2][1]*rr[2][1];
    mm2[4] = ll[0]*rr[0][1]*rr[0][2] + ll[1]*rr[1][1]*rr[1][2] + ll[2]*rr[2][1]*rr[2][2];
    mm2[5] = ll[0]*rr[0][2]*rr[0][2] + ll[1]*rr[1][2]*rr[1][2] + ll[2]*rr[2][2]*rr[2][2];
  }
  else if ( p2->tag & MG_GEO ){
    if ( !MMG5_updatemetreq_ani(mtan2,mu,vp) ) { return 0; }

    /* Here mtan2 contains the gradated metric in the coreduction basis: compute
     * the sizes in the directions (t,u=t^n,n): Computation in 3D but it is
     * maybe more efficient to work in the tangent plane (but we need to compute
     * the basis of the ridge metric in the tangent plane) */
    rmet3D[0] = mtan2[0];
    rmet3D[1] = mtan2[1];
    rmet3D[2] = 0;
    rmet3D[3] = mtan2[2];
    rmet3D[4] = 0;
    rmet3D[5] = mr2[5];

    mu[0] = rmet3D[0]*rbasis2[0][0]*rbasis2[0][0] + 2. * rmet3D[1]*rbasis2[1][0]*rbasis2[0][0]
      + 2. * rmet3D[2]*rbasis2[2][0]*rbasis2[0][0]
      + rmet3D[3]*rbasis2[1][0]*rbasis2[1][0] + 2. * rmet3D[4]*rbasis2[2][0]*rbasis2[1][0]
      + rmet3D[5]*rbasis2[2][0]*rbasis2[2][0];

    /* h = 1/sqrt(t_e M e) */
    assert ( mu[0] > MMG5_EPSD2 );

    mu[1] = rmet3D[0]*rbasis2[0][1]*rbasis2[0][1] + 2. * rmet3D[1]*rbasis2[1][1]*rbasis2[0][1]
      + 2. * rmet3D[2]*rbasis2[2][1]*rbasis2[0][1]
      + rmet3D[3]*rbasis2[1][1]*rbasis2[1][1] + 2. * rmet3D[4]*rbasis2[2][1]*rbasis2[1][1]
      + rmet3D[5]*rbasis2[2][1]*rbasis2[2][1];

    /* h = 1/sqrt(t_e M e) */
    assert ( mu[1] > MMG5_EPSD2 );

    /* Update the ridge metric */
    mm2[0] =  mu[0];

    assert ( cfg_m2 );
    mm2[cfg_m2] = mu[1];

  }
  else{
    /* Update of the metrics */
    mu[2] = mr2[5];

    if ( !MMG5_updatemetreq_ani(mtan2,mu,vp) ) { return 0; }

    /* Return in initial basis */
    mtmp[0][0] = mtan2[0]*r2[0][0] + mtan2[1]*r2[1][0];
    mtmp[0][1] = mtan2[0]*r2[0][1] + mtan2[1]*r2[1][1];
    mtmp[0][2] = mtan2[0]*r2[0][2] + mtan2[1]*r2[1][2];

    mtmp[1][0] = mtan2[1]*r2[0][0] + mtan2[2]*r2[1][0];
    mtmp[1][1] = mtan2[1]*r2[0][1] + mtan2[2]*r2[1][1];
    mtmp[1][2] = mtan2[1]*r2[0][2] + mtan2[2]*r2[1][2];

    mtmp[2][0] =  mr2[5]*r2[2][0];
    mtmp[2][1] =  mr2[5]*r2[2][1];
    mtmp[2][2] =  mr2[5]*r2[2][2];

    m2[0] = r2[0][0]*mtmp[0][0] + r2[1][0]*mtmp[1][0] + r2[2][0]*mtmp[2][0];
    m2[1] = r2[0][0]*mtmp[0][1] + r2[1][0]*mtmp[1][1] + r2[2][0]*mtmp[2][1];
    m2[2] = r2[0][0]*mtmp[0][2] + r2[1][0]*mtmp[1][2] + r2[2][0]*mtmp[2][2];

    m2[3] = r2[0][1]*mtmp[0][1] + r2[1][1]*mtmp[1][1] + r2[2][1]*mtmp[2][1];
    m2[4] = r2[0][1]*mtmp[0][2] + r2[1][1]*mtmp[1][2] + r2[2][1]*mtmp[2][2];

    m2[5] = r2[0][2]*mtmp[0][2] + r2[1][2]*mtmp[1][2] + r2[2][2]*mtmp[2][2];

#ifndef NDEBUG
    /* Check the validity of the output metric */
    ier = MMG5_eigenv3d(1,m2,mu, r2);

    assert ( ier );
    assert ( mu[0] > 0.);
    assert ( mu[1] > 0.);
    assert ( mu[2] > 0.);
#endif

    memcpy(mm2,m2,6*sizeof(double));
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Compute the mean metric at mesh points with a non-nul \a s field. At the
 * beginning, for a given point \a ip, \f$ met->m[met->size * ip] \f$ contains
 * the sum of n metrics and the \a s field of \a ip contains the number of
 * metrics summed in the point. Set the flag of the processed points to 3.
 *
 */
int MMG5_compute_meanMetricAtMarkedPoints_ani ( MMG5_pMesh mesh,MMG5_pSol met ) {
  MMG5_pPoint p0;
  double      lm;
  MMG5_int    k,iadr;
  int         mmgWarn = 0;

  for ( k=1; k<=mesh->np; k++ ) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) )  continue;

    if ( !p0->s ) continue;

    iadr = met->size*k;
    lm   = p0->s/met->m[iadr];
    met->m[iadr] = lm*lm;

    if ( mesh->dim==2 ) {
      met->m[iadr+2] = met->m[iadr];
    }
    else if ( !MG_RID(p0->tag) ) {
      /* Classic metric */
      met->m[iadr+3] = met->m[iadr+5] = met->m[iadr];
    }
    else {
      /* Ridge metric */
      met->m[iadr+2] = met->m[iadr+1] = met->m[iadr];
      met->m[iadr+4] = met->m[iadr+3] = met->m[iadr];
    }

    p0->flag = 3;

    /* Warn the user that edge size is erased */
    if ( !mmgWarn ) {
      mmgWarn = 1;
      if ( mesh->info.ddebug || (mesh->info.imprim > 4) ) {
        printf("\n  -- SIZEMAP CORRECTION : overwritten of sizes at required vertices\n");
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param it number of performed iteration (to fill)
 *
 * \return nup, the number of points updated.
 *
 *
 * Standard gradation procedure.
 *
 */
MMG5_int MMG5_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met,int *it) {
  MMG5_pTria   pt;
  MMG5_pPoint  p1,p2;
  int          maxit;
  MMG5_int     ier,k,np1,np2,nup,nu;
  int8_t       i;

  /** Mark the edges belonging to a required entity */
  MMG5_mark_pointsOnReqEdge_fromTria ( mesh );

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  (*it) = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        np1 = pt->v[MMG5_inxt2[i]];
        np2 = pt->v[MMG5_iprv2[i]];
        p1 = &mesh->point[np1];
        p2 = &mesh->point[np2];

        if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;
        /* Skip points belonging to a required edge */
        if ( p1->s || p2->s ) continue;

        ier = MMG5_grad2met_ani(mesh,met,pt,np1,np2);
        if ( ier == np1 ) {
          p1->flag = mesh->base;
          nu++;
        }
        else if ( ier == np2 ) {
          p2->flag = mesh->base;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while( ++(*it) < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     gradation: %7"MMG5_PRId" updated, %d iter.\n",nup,(*it));

  return nup;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 1
 *
 *
 * Enforces mesh gradation (on required entities) by truncating metric field.
 *
 */
int MMG5_gradsizreq_ani(MMG5_pMesh mesh,MMG5_pSol met) {

  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  int               it,maxit,ier;
  MMG5_int          k,np1,np2,npslave,npmaster,nup,nu;
  int8_t            i;


  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  ** Grading required points.\n");
  }

  if ( mesh->info.hgrad < 0. ) {
    /** Mark the edges belonging to a required entity (already done if the
     * classic gradation is enabled) */
    MMG5_mark_pointsOnReqEdge_fromTria ( mesh );
  }

  it = nup = 0;
  maxit = 100;

  do {
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        np1 = pt->v[MMG5_inxt2[i]];
        np2 = pt->v[MMG5_iprv2[i]];
        p1 = &mesh->point[np1];
        p2 = &mesh->point[np2];

        if ( MMG5_abs ( p1->s - p2->s ) < 2 ) {
          /* No size to propagate */
          continue;
        }
        else if ( p1->s > p2->s ) {
          npmaster = np1;
          npslave  = np2;
        }
        else {
          assert ( p2->s > p1->s );
          npmaster = np2;
          npslave  = np1;
        }

        /* Impose the gradation to npslave from npmaster: coming from mmgs,
         * MMG5_grad2metreq_ani is a pointer toward MMG5_grad2metSurfreq,
         * comming from mmg2d, it is a pointer toward MMG2D_grad2metreq_ani */
        ier = MMG5_grad2metreq_ani(mesh,met,pt,npmaster,npslave);

        if ( ier ) {
          mesh->point[npslave].s = mesh->point[npmaster].s - 1;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 && nup ) {
    fprintf(stdout,"     gradation (required): %7"MMG5_PRId" updated, %d iter.\n",nup,it);
  }

  return 1;
}
