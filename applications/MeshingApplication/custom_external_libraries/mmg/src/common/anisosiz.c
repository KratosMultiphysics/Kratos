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
 * \file common/anisosiz.c
 * \brief Fonctions for anisotropic size map computation.
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
 * \param m pointer toward the metric at triangle vertices.
 * \param ptt pointer toward the triangle structure.
 * \return The double of the triangle area.
 *
 * Compute the double of the area of the surface triangle \a ptt with respect to
 * the anisotropic metric \a m.
 *
 */
static inline
double _MMG5_surf(MMG5_pMesh mesh,double m[3][6],MMG5_pTria ptt) {
  _MMG5_Bezier   b;
  double         surf,dens,J[3][2],mJ[3][2],tJmJ[2][2];
  char           i;

  surf = 0.0;

  if ( !_MMG5_bezierCP(mesh,ptt,&b,1) ) return(0.0);


  /* Compute density integrand of volume at the 3 vertices of T */
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
    /* if ( dens < 0.0 ) { */
    /*   fprintf(stdout,"  ## Density should be positive : %E for elt %d %d %d \n",dens,ptt->v[0],ptt->v[1],ptt->v[2]); */
    /*   return(0.0); */
    /* } */
    surf += sqrt(fabs(dens));
  }

  surf *= _MMG5_ATHIRD;
  return(surf);
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
double _MMG5_surftri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  MMG5_pPoint    p[3];
  int            np[3];
  double         ux,uy,uz,m[3][6];
  char           i1,i2;
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
      i1 = _MMG5_inxt2[i];
      i2 = _MMG5_iprv2[i];
      ux = 0.5*(p[i1]->c[0]+p[i2]->c[0]) - p[i]->c[0];
      uy = 0.5*(p[i1]->c[1]+p[i2]->c[1]) - p[i]->c[1];
      uz = 0.5*(p[i1]->c[2]+p[i2]->c[2]) - p[i]->c[2];
      if ( !_MMG5_buildridmet(mesh,met,np[i],ux,uy,uz,&m[i][0]) )  return(0.0);
    }
    else {
      memcpy(&m[i][0],&met->m[6*np[i]],6*sizeof(double));
    }
  }
  return(_MMG5_surf(mesh,m,ptt));

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
double _MMG5_surftri33_ani(MMG5_pMesh mesh,MMG5_pTria ptt,
                           double ma[6], double mb[6], double mc[6]) {
  double         mm[6];
  double         *a,*b,*c,abx,aby,abz,acx,acy,acz,dens[3],surf;
  int            i,ia,ib,ic;

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
    mm[i] = _MMG5_ATHIRD * (ma[i] + mb[i]+ mc[i]);

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

  if ( surf < _MMG5_EPSD ) return(0.0);

  surf = sqrt(surf);

  return(surf);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ismet has the user provided a metric?
 *
 * Search for points with unintialized metric and define anisotropic size at
 * this points.
 *
 */
void _MMG5_defUninitSize(MMG5_pMesh mesh,MMG5_pSol met,char ismet)
{
  MMG5_pPoint   ppt;
  double        *m,*n,r[3][3],isqhmax;
  int           k;

  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) || ppt->flag == 1 )  continue;

    m = &met->m[6*k];
    if(ismet) {
      if ( !MG_SIN(ppt->tag) && (ppt->tag & MG_GEO) ) {
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
      _MMG5_rotmatrix(n,r);
      m[0] = isqhmax*(r[0][0]*r[0][0]+r[1][0]*r[1][0]+r[2][0]*r[2][0]);
      m[1] = isqhmax*(r[0][0]*r[0][1]+r[1][0]*r[1][1]+r[2][0]*r[2][1]);
      m[2] = isqhmax*(r[0][0]*r[0][2]+r[1][0]*r[1][2]+r[2][0]*r[2][2]);
      m[3] = isqhmax*(r[0][1]*r[0][1]+r[1][1]*r[1][1]+r[2][1]*r[2][1]);
      m[4] = isqhmax*(r[0][1]*r[0][2]+r[1][1]*r[1][2]+r[2][1]*r[2][2]);
      m[5] = isqhmax*(r[0][2]*r[0][2]+r[1][2]*r[1][2]+r[2][2]*r[2][2]);
    }
    ppt->flag = 1;
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
void _MMG5_fillDefmetregSys( int k, MMG5_pPoint p0, int i0, _MMG5_Bezier b,
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
int _MMG5_solveDefmetregSys( MMG5_pMesh mesh, double r[3][3], double c[3],
                             double tAA[6], double tAb[3], double *m,
                             double isqhmin, double isqhmax, double hausd)
{
  double intm[3], kappa[2], vp[2][2], b0[3], b1[3], b2[3];

  memset(intm,0x0,3*sizeof(double));

  /* case planar surface : tAb = 0 => no curvature */
  /* isotropic metric with hmax size*/

  if((tAb[0]*tAb[0] + tAb[1]*tAb[1] + tAb[2]*tAb[2]) < _MMG5_EPSD) {
    m[0] = isqhmax;
    m[1] = 0;
    m[2] = 0;
    m[3] = isqhmax;
    m[4] = 0;
    m[5] = isqhmax;
    return(1);
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !_MMG5_sys33sym(tAA,tAb,c) ) {
    printf("%s:%d: Warning: unable to solve the system.\n",__FILE__,__LINE__);
    return(0);
  }
  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach : vp[0] and vp[1]
     are the two pr. directions of curvature, and the two curvatures can be inferred from lambdas*/
  _MMG5_eigensym(intm,kappa,vp);

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

  /* Security check : normal in the kernel */
  /*if((fabs(p0->m[0]*n[0] + p0->m[1]*n[1] + p0->m[2]*n[2] ) > 1.0e-4)){
    printf("VALEUR ETRANGE... %f \n",fabs(p0->m[0]*n[0] + p0->m[1]*n[1] + p0->m[2]*n[2] ));
    }
    if((fabs(p0->m[1]*n[0] + p0->m[3]*n[1] + p0->m[4]*n[2] ) > 1.0e-4)){
    printf("VALEUR ETRANGE... %f \n",fabs(p0->m[1]*n[0] + p0->m[3]*n[1] + p0->m[4]*n[2] ));
    }

    if((fabs(p0->m[2]*n[0] + p0->m[4]*n[1] + p0->m[5]*n[2] ) > 1.0e-4)){
    printf("VALEUR ETRANGE... %f \n",fabs(p0->m[2]*n[0] + p0->m[4]*n[1] + p0->m[5]*n[2]));
    } */

  /*if(ddb) {
    printf("La matrice %f %f %f\n",p0->m[0],p0->m[1],p0->m[2]);
    printf("            %f %f %f\n",p0->m[1],p0->m[3],p0->m[4]);
    printf("            %f %f %f\n",p0->m[2],p0->m[4],p0->m[5]);

    }*/
  return(1);
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
 * \param hausd hausdorff value at point.
 * \return 1 if success, 0 if fail.
 *
 * Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) for a ref
 * point.
 *
 */
int _MMG5_solveDefmetrefSys( MMG5_pMesh mesh, MMG5_pPoint p0, int ipref[2],
                             double r[3][3], double c[3],
                             double tAA[6], double tAb[3], double *m,
                             double isqhmin, double isqhmax, double hausd)
{
  MMG5_pPoint  p1;
  double       intm[3], kappa[2], vp[2][2], b0[3], b1[3], b2[3], kappacur;
  double       gammasec[3],tau[2], ux, uy, uz, ps1, l, ll, *t, *t1;
  int          i;

  memset(intm,0x0,3*sizeof(double));

  /* case planar surface : tAb = 0 => no curvature */
  /* isotropic metric with hmax size*/
  if((tAb[0]*tAb[0] + tAb[1]*tAb[1] + tAb[2]*tAb[2]) < _MMG5_EPSD) {
    m[0] = isqhmax;
    m[1] = 0;
    m[2] = 0;
    m[3] = isqhmax;
    m[4] = 0;
    m[5] = isqhmax;
    return(1);
  }

  /* solve now (a b c) = tAA^{-1} * tAb */
  if ( !_MMG5_sys33sym(tAA,tAb,c) ) {
    printf("%s:%d: Warning: unable to solve the system.\n",__FILE__,__LINE__);
    return(0);
  }
  intm[0] = 2.0*c[0];
  intm[1] = c[2];
  intm[2] = 2.0*c[1];

  /* At this point, intm stands for the integral matrix of Taubin's approach :
     vp[0] and vp[1] are the two pr. directions of curvature, and the two
     curvatures can be inferred from lambdas*/
  _MMG5_eigensym(intm,kappa,vp);

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
    c[0] = _MMG5_ATHIRD*ps1*t[0];
    c[1] = _MMG5_ATHIRD*ps1*t[1];
    c[2] = _MMG5_ATHIRD*ps1*t[2];

    b0[0] =  r[0][0]*c[0] + r[0][1]*c[1] + r[0][2]*c[2];
    b0[1] =  r[1][0]*c[0] + r[1][1]*c[1] + r[1][2]*c[2];
    b0[2] =  r[2][0]*c[0] + r[2][1]*c[1] + r[2][2]*c[2];

    if ( (MG_CRN & p1->tag) || (MG_NOM & p1->tag) ) {
      c[0] = p1->c[0] - _MMG5_ATHIRD*ux;
      c[1] = p1->c[1] - _MMG5_ATHIRD*uy;
      c[2] = p1->c[2] - _MMG5_ATHIRD*uz;
    }
    else {
      assert(MG_REF & p1->tag);
      t1 = &(p1->n[0]);
      ps1 =  -(ux*t1[0] + uy*t1[1] + uz*t1[2]);
      c[0] = p1->c[0] + _MMG5_ATHIRD*ps1*t1[0];
      c[1] = p1->c[1] + _MMG5_ATHIRD*ps1*t1[1];
      c[2] = p1->c[2] + _MMG5_ATHIRD*ps1*t1[2];
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
    if ( ll < _MMG5_EPSD ) {
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
  _MMG5_intmetsavedir(mesh,c,intm,b0);
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

  return(1);
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
 * Compute the specific size of a ridge in the direction of the tangent of the
 * ridge.
 *
 **/
double _MMG5_ridSizeInTangentDir(MMG5_pMesh mesh, MMG5_pPoint p0, int idp,
                                 int* iprid, double isqhmin,double isqhmax)
{
  int    i;
  double n0[3],tau[3],gammasec[3],c[3],ps,ll,l,m;
  double b0[3],b1[3],kappacur;

  m = isqhmax;
  for (i=0; i<2; i++) {
    kappacur = 0.0;
    // Remark: bezierEdge don't use n0 in case of a ridge so it's ok to call it
    // with an undefined n0.
    _MMG5_bezierEdge(mesh,idp,iprid[i],b0,b1,1,n0);

    /* tau is the bezier curve derivative in p0 (parametric coor t=0) */
    tau[0] = 3.0*(b0[0] - p0->c[0]);
    tau[1] = 3.0*(b0[1] - p0->c[1]);
    tau[2] = 3.0*(b0[2] - p0->c[2]);
    ll = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];
    if ( ll < _MMG5_EPSD )  continue;
    l = 1.0 / sqrt(ll);
    tau[0] *= l;
    tau[1] *= l;
    tau[2] *= l;

    gammasec[0] = 6.0*p0->c[0] -12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = 6.0*p0->c[1] -12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = 6.0*p0->c[2] -12.0*b0[2] + 6.0*b1[2];

    ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];
    c[0] = gammasec[0] - ps*tau[0];
    c[1] = gammasec[1] - ps*tau[1];
    c[2] = gammasec[2] - ps*tau[2];

    kappacur = MG_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));
    kappacur = 1.0/8.0*kappacur/mesh->info.hausd;
    kappacur = MG_MIN(kappacur,isqhmin);
    kappacur = MG_MAX(kappacur,isqhmax);
    m = MG_MAX(m,kappacur);
  }
  return(m);
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
 * looked face.
 *
 **/
double _MMG5_ridSizeInNormalDir(MMG5_pMesh mesh,int i0,double* bcu,
                                _MMG5_Bezier *b,double isqhmin,double isqhmax)
{
  double lambda[2],Jacb[3][2],Hb[3][3],tau[3],ll,l,gammasec[3],c[3];
  double ps,kappacur;

  if ( i0 == 0 ) { // w = 1, u,v = 0
    lambda[0] = bcu[1];
    lambda[1] = bcu[2];

    Jacb[0][0] = 3.0*(b->b[7][0]-b->b[0][0]);
    Jacb[1][0] = 3.0*(b->b[7][1]-b->b[0][1]);
    Jacb[2][0] = 3.0*(b->b[7][2]-b->b[0][2]);

    Jacb[0][1] = 3.0*(b->b[6][0]-b->b[0][0]);
    Jacb[1][1] = 3.0*(b->b[6][1]-b->b[0][1]);
    Jacb[2][1] = 3.0*(b->b[6][2]-b->b[0][2]);

    /* Hb[i] = hessian matrix of i-th component of b at point p0 */
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
  if ( ll < _MMG5_EPSD )  return(0);

  l = 1.0 / sqrt(ll);
  tau[0] *= l;
  tau[1] *= l;
  tau[2] *= l;

  gammasec[0] = Hb[0][0]*lambda[0]*lambda[0] + 2.0*Hb[0][1]*lambda[0]*lambda[1] + Hb[0][2]*lambda[1]*lambda[1];
  gammasec[1] = Hb[1][0]*lambda[0]*lambda[0] + 2.0*Hb[1][1]*lambda[0]*lambda[1] + Hb[1][2]*lambda[1]*lambda[1];
  gammasec[2] = Hb[2][0]*lambda[0]*lambda[0] + 2.0*Hb[2][1]*lambda[0]*lambda[1] + Hb[2][2]*lambda[1]*lambda[1];

  ps = tau[0]*gammasec[0] + tau[1]*gammasec[1] + tau[2]*gammasec[2];
  c[0] = gammasec[0] - ps*tau[0];
  c[1] = gammasec[1] - ps*tau[1];
  c[2] = gammasec[2] - ps*tau[2];

  kappacur = MG_MAX(0.0,1.0/ll*sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]));
  kappacur = 1.0/8.0 * kappacur/mesh->info.hausd;
  kappacur = MG_MIN(kappacur,isqhmin);
  kappacur = MG_MAX(kappacur,isqhmax);

  return(kappacur);
}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param pt pointer toward a triangle.
 * \param i edge index in triangle \a pt.
 * \return -1 if no gradation is needed, else index of graded point.
 *
 * Enforces gradation of metric in one extremity of edge \a i in tria \a pt
 * with respect to the other, along the direction of the associated support
 * curve first, then along the normal direction.
 *
 * \warning The gradation along the direction normal to the surface is made in
 * an "isotropic way".
 *
 */
int _MMG5_grad2metSurf(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pTria pt, int i)
{
  MMG5_pPoint   p1,p2;
  double   *mm1,*mm2,*nn1,*nn2,ps1,ps2,ux,uy,uz,m1[6],m2[6],n1[3],n2[3],nt[3];
  double   r1[3][3],r2[3][3],t1[2],t2[2],c[3],mtan1[3],mtan2[3],mr1[6],mr2[6];
  double   mtmp[3][3],val;
  double   /*,l1,l2*/l,dd;
  double   lambda[2],vp[2][2],alpha,beta,mu[3];
  int      np1,np2,kmin,idx;
  char     i1,i2,ichg;

  i1 = _MMG5_inxt2[i];
  i2 = _MMG5_iprv2[i];
  np1 = pt->v[i1];
  np2 = pt->v[i2];

  p1 = &mesh->point[np1];
  p2 = &mesh->point[np2];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  uz = p2->c[2] - p1->c[2];

  mm1 = &met->m[6*np1];
  mm2 = &met->m[6*np2];

  if( !_MMG5_nortri(mesh,pt,nt) )
    return(-1);

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

    if( !_MMG5_buildridmet(mesh,met,np1,ux,uy,uz,m1) )
      return(-1);
  }
  else if( ( MG_REF & p1->tag ) ){
    memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
    memcpy(m1,mm1,6*sizeof(double));
  }
  else{
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

    if( !_MMG5_buildridmet(mesh,met,np2,ux,uy,uz,m2) )
      return(-1);
  }
  else if( (MG_REF & p2->tag) ){
    memcpy(n2,&(mesh->xpoint[p2->xp].n1[0]),3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }
  else{
    memcpy(n2,p2->n,3*sizeof(double));
    memcpy(m2,mm2,6*sizeof(double));
  }

  /* Rotation matrices mapping n1/n2 to e_3 */
  _MMG5_rotmatrix(n1,r1);
  _MMG5_rotmatrix(n2,r2);

  /* Geodesic length of support curve to edge i */

  // Commentated because the next line overwrite it... to check!
  /* ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2]; */
  /* t1[0] = ux - ps1*n1[0]; */
  /* t1[1] = uy - ps1*n1[1]; */
  /* t1[2] = uz - ps1*n1[2]; */

  /* // warning : t2 seems to be wrong calculated */
  /* ps2 = - (ux*n2[0] + uy*n2[1] + uz*n2[2]); */
  /* t2[0] = -ux - ps1*n2[0]; */
  /* t2[1] = -uy - ps1*n2[1]; */
  /* t2[2] = -uz - ps1*n2[2]; */

  /* l1 = m1[0]*t1[0]*t1[0] + m1[3]*t1[1]*t1[1] + m1[5]*t1[2]*t1[2] \ */
  /*   + 2.0 * ( m1[1]*t1[0]*t1[1] + m1[2]*t1[0]*t1[2] + m1[4]*t1[1]*t1[2] ) ; */
  /* l2 = m2[0]*t2[0]*t2[0] + m2[3]*t2[1]*t2[1] + m2[5]*t2[2]*t2[2] \ */
  /*   + 2.0 * ( m2[1]*t2[0]*t2[1] + m2[2]*t2[0]*t2[2] + m2[4]*t2[1]*t2[2] ) ; */
  /* l = 0.5* ( sqrt(l1) + sqrt(l2) ) ; */

  l = sqrt(ux*ux+uy*uy+uz*uz);

  /* Characteristic sizes in direction of support curve */
  _MMG5_rmtr(r1,m1,mr1);

  mtan1[0] = mr1[0];
  mtan1[1] = mr1[1];
  mtan1[2] = mr1[3];

  c[0] = r1[0][0]*ux + r1[0][1]*uy + r1[0][2]*uz;
  c[1] = r1[1][0]*ux + r1[1][1]*uy + r1[1][2]*uz;

  memcpy(t1,c,2*sizeof(double));
  // Here we work in the tangent plane (thus in 2d)
  dd = t1[0]*t1[0] + t1[1]*t1[1];
  if(dd < _MMG5_EPSD2)
    return(-1);

  dd = 1.0/sqrt(dd);
  t1[0] *= dd;
  t1[1] *= dd;

  // edge length in metric mtan1: sqrt(t^(t1) * mtan1 * t1).
  ps1 =  mtan1[0]*t1[0]*t1[0] + 2.0*mtan1[1]*t1[0]*t1[1] + mtan1[2]*t1[1]*t1[1];
  ps1 = sqrt(ps1);

  _MMG5_rmtr(r2,m2,mr2);

  mtan2[0] = mr2[0];
  mtan2[1] = mr2[1];
  mtan2[2] = mr2[3];

  c[0] = - ( r2[0][0]*ux + r2[0][1]*uy + r2[0][2]*uz );
  c[1] = - ( r2[1][0]*ux + r2[1][1]*uy + r2[1][2]*uz );
  memcpy(t2,c,2*sizeof(double));

  dd = t2[0]*t2[0] + t2[1]*t2[1];
  if(dd < _MMG5_EPSD2)
    return(-1);

  dd = 1.0/sqrt(dd);
  t2[0] *= dd;
  t2[1] *= dd;

  // edge length: sqrt(t^(t2) * mtan2 * t2)
  ps2 = mtan2[0]*t2[0]*t2[0] + 2.0*mtan2[1]*t2[0]*t2[1] + mtan2[2]*t2[1]*t2[1];
  ps2 = sqrt(ps2);

  /* Metric in p1 has to be changed */
  if( ps2 > ps1 ){
    /* compute alpha = h2 + hgrad*l */
    alpha = ps2 /(1.0+mesh->info.hgrad*l*ps2);
    if( ps1 >= alpha -_MMG5_EPS )
      return(-1);

    _MMG5_eigensym(mtan1,lambda,vp);
    /* Project the vector t1 along the main directions of the metric */
    // Remark: along the third direction mr1 is already diagonal,
    // thus vp[2][.] =( 0 0 1) and vp[.][2] = 0.
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
    assert(c[ichg]*c[ichg] > _MMG5_EPS );
   /* Compute beta coef such as lambda_1 = beta*lambda_1 => h1 = h2 + hgrad*l
    * (see p317 of Charles Dapogny Thesis). */
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

      // Find index of the minimum value of c
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
      // Because of the rotation, we know that:
      // mr.[2] = mr.[4]= 0
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
    return(i1);
  }
  /* Metric in p2 has to be changed */
  else{
    alpha = ps1 /(1.0+mesh->info.hgrad*l*ps1);
    if( ps2 >= alpha - _MMG5_EPS)
      return(-1);

    _MMG5_eigensym(mtan2,lambda,vp);

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
    assert(c[ichg]*c[ichg] > _MMG5_EPS );

   /* Compute beta coef such as lambda_1 = beta*lambda_1 => h1 = h2 + hgrad*l
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
    return(i2);
  }
}
