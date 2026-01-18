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
 * \brief inlined Functions
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg3d_private.h"
#include "inlined_functions_private.h"
#include "mmg3dexterns_private.h"

#ifndef _INLINED_FUNCT_3D_H
#define _INLINED_FUNCT_3D_H


/**
 * \brief Compute edge length from edge's coordinates.
 * \param ca pointer to the coordinates of the first edge's extremity.
 * \param cb pointer to the coordinates of the second edge's extremity.
 * \param sa pointer to the metric associated to the first edge's extremity.
 * \param sb pointer to the metric associated to the second edge's extremity.
 * \return edge length.
 *
 * Compute length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge extremities) according to the anisotropic size
 * prescription.
 *
 */
static
inline double MMG5_lenedgCoor_ani(double *ca,double *cb,double *sa,double *sb) {
  double   ux,uy,uz,dd1,dd2,len;

  ux = cb[0] - ca[0];
  uy = cb[1] - ca[1];
  uz = cb[2] - ca[2];

  dd1 =      sa[0]*ux*ux + sa[3]*uy*uy + sa[5]*uz*uz \
    + 2.0*(sa[1]*ux*uy + sa[2]*ux*uz + sa[4]*uy*uz);
  if ( dd1 <= 0.0 )  dd1 = 0.0;

  dd2 =      sb[0]*ux*ux + sb[3]*uy*uy + sb[5]*uz*uz \
    + 2.0*(sb[1]*ux*uy + sb[2]*ux*uz + sb[4]*uy*uz);
  if ( dd2 <= 0.0 )  dd2 = 0.0;

  /*longueur approchee*/
  /*precision a 3.5 10e-3 pres*/
  if(fabs(dd1-dd2) < 0.05 ) {
    len = sqrt(0.5*(dd1+dd2));
    return len;
  }
  len = (sqrt(dd1)+sqrt(dd2)+4.0*sqrt(0.5*(dd1+dd2))) / 6.0;

  return len;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param ia index of edge in tetra \a pt .
 * \param pt pointer to the tetra from which we come.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of edge \f$[i0;i1]\f$ according to the prescribed aniso
 * metric (for classic storage of metrics at ridges points).
 *
 * \warning in this function we may erroneously approximate the length of a
 * curve boundary edge by the length of the straight edge if the "MG_BDY" tag is
 * missing along the edge.
 */
static
inline double MMG5_lenedg33_ani(MMG5_pMesh mesh ,MMG5_pSol met, int ia,
                                 MMG5_pTetra pt)
{
  MMG5_int    ip1,ip2;
  int8_t      isedg;

  ip1 = pt->v[MMG5_iare[ia][0]];
  ip2 = pt->v[MMG5_iare[ia][1]];

  if ( pt->xt && (mesh->xtetra[pt->xt].tag[ia] & MG_BDY)) {
    isedg = ( mesh->xtetra[pt->xt].tag[ia] & MG_GEO);
    // Computation of the length of a curve edge with 33 aniso metric
    return MMG5_lenSurfEdg33_ani(mesh, met, ip1, ip2, isedg);
  } else {
    // Computation for an internal edge with 33 aniso metric
    return MMG5_lenedgCoor_ani(mesh->point[ip1].c,mesh->point[ip2].c,
                                &met->m[6*ip1],&met->m[6*ip2]);
  }
  return 0.0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param ia index of edge in tetra \a pt .
 * \param pt pointer to the tetra from which we come.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of edge \f$[i0;i1]\f$ according to the prescribed aniso
 * metric (for classic storage of metrics at ridges points).
 *
 */
static
inline double MMG5_lenedgspl33_ani(MMG5_pMesh mesh ,MMG5_pSol met, int ia,
                                    MMG5_pTetra pt)
{
  MMG5_pPoint pp1,pp2;
  double      *m1,*m2;
  MMG5_int    ip1,ip2;

  ip1 = pt->v[MMG5_iare[ia][0]];
  ip2 = pt->v[MMG5_iare[ia][1]];

  pp1 = &mesh->point[ip1];
  pp2 = &mesh->point[ip2];

  m1 = &met->m[6*ip1];
  m2 = &met->m[6*ip2];

  return MMG5_lenedgCoor_ani(pp1->c,pp2->c,m1,m2);
}



/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param ia index of edge in tetra \a pt .
 * \param pt pointer to the tetra from which we come.
 * \return length of edge according to the prescribed metric, 0 if fail.
 *
 * Compute length of a straight edge \f$[i0;i1]\f$ according to the prescribed
 * aniso metric (for special storage of metrics at ridges points).
 *
 */
static
inline double MMG5_lenedgspl_ani(MMG5_pMesh mesh ,MMG5_pSol met, int ia,
                                  MMG5_pTetra pt)
{
  MMG5_pPoint pp1,pp2;
  double      m1[6],m2[6];
  MMG5_int    ip1,ip2;
  int         i;

  ip1 = pt->v[MMG5_iare[ia][0]];
  ip2 = pt->v[MMG5_iare[ia][1]];

  pp1 = &mesh->point[ip1];
  pp2 = &mesh->point[ip2];

  if ( MG_RID(pp1->tag) ) {
    if ( !MMG5_moymet(mesh,met,pt,m1) ) return 0;
  } else {
    for ( i=0; i<6; ++i )
      m1[i] = met->m[6*ip1+i];
  }

  if ( MG_RID(pp2->tag) ) {
    if ( !MMG5_moymet(mesh,met,pt,m2) ) return 0;
  } else {
    for ( i=0; i<6; ++i )
      m2[i] = met->m[6*ip2+i];
  }

  return MMG5_lenedgCoor_ani(pp1->c,pp2->c,m1,m2);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param ia index of edge in tetra \a pt .
 * \param pt pointer to the tetra from which we come.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of edge \f$[i0;i1]\f$ according to the prescribed aniso
 * metric (for special storage of metrics at ridges points).
 *
 * \warning in this function we may erroneously approximate the length of a
 * curve boundary edge by the length of the straight edge if the "MG_BDY" tag is
 * missing along the edge.
 *
 */
static
inline double MMG5_lenedg_ani(MMG5_pMesh mesh ,MMG5_pSol met, int ia,
                               MMG5_pTetra pt)
{
  MMG5_int    ip1,ip2;
  int8_t      isedg;

  ip1 = pt->v[MMG5_iare[ia][0]];
  ip2 = pt->v[MMG5_iare[ia][1]];

  if ( pt->xt && (mesh->xtetra[pt->xt].tag[ia] & MG_BDY)) {
    isedg = ( mesh->xtetra[pt->xt].tag[ia] & MG_GEO);
    // Computation of the length of a curve edge with ridge metric
    return MMG5_lenSurfEdg_ani(mesh, met, ip1, ip2, isedg);
  } else {
    // Computation for an internal edge with ridge metric
    return MMG5_lenedgspl_ani(mesh ,met, ia, pt);
  }
  return 0.0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param ia index of edge in tetra \a pt .
 * \param pt pointer to the tetra from which we come.
 * \return length of edge according to the prescribed metric.
 *
 * Compute length of edge \f$[i0;i1]\f$ according to the prescribed iso
 * metric (length identic for internal edges than for surface edges).
 *
 */
static
inline double MMG5_lenedg_iso(MMG5_pMesh mesh,MMG5_pSol met,int ia,
                               MMG5_pTetra pt) {
  MMG5_int ip1,ip2;

  ip1 = pt->v[MMG5_iare[ia][0]];
  ip2 = pt->v[MMG5_iare[ia][1]];
  return MMG5_lenSurfEdg_iso(mesh,met,ip1,ip2,0);
}

static
inline double MMG5_lenedgspl_iso(MMG5_pMesh mesh ,MMG5_pSol met, int ia,
                                  MMG5_pTetra pt) {
  MMG5_int ip1,ip2;

  ip1 = pt->v[MMG5_iare[ia][0]];
  ip2 = pt->v[MMG5_iare[ia][1]];

  return MMG5_lenSurfEdg_iso(mesh,met,ip1,ip2,0);

}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the meric structure.
 * \param iel index of element.
 * \return The oriented quality of element \a iel or 0.0 if \a iel is inverted.
 *
 * Compute tetra oriented quality of iel.
 *
 */
static
inline double MMG5_orcal(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int iel) {
  MMG5_pTetra     pt;

  pt = &mesh->tetra[iel];

  return MMG5_caltet(mesh,met,pt);
}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the meric structure.
 * \param pt pointer to a tetrahedra.
 * \return The isotropic quality of the tet in LES measure, 0 if fail
 *
 * Compute the quality of the tet pt with respect to the LES quality measure.
 * \f$Q=\frac{V}{V_{ref}}\f$ with \f$V_{ref}=8\frac{\sqrt{3}}{27}R
 * \sqrt{R}\f$ and R the radius of the circumscribe circle.
 *
 */
static
inline double MMG3D_caltetLES_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt) {
  double    ct[12],cs[3],rad,Vref,V,cal;
  int j,l;

  for (j=0,l=0; j<4; j++,l+=3) {
    memcpy(&ct[l],mesh->point[pt->v[j]].c,3*sizeof(double));
  }

  if(!MMG5_cenrad_iso(mesh,ct,cs,&rad)) {
    return 0.0;
  }

  assert(rad>0.);

  /* Vref volume */
  Vref = 8.*sqrt(3)/27.*rad*sqrt(rad);

  V = MMG5_orvol(mesh->point,pt->v) * MMG3D_DET2VOL;

  if ( V<0. ) {
    return 0.0;
  }

  cal =  V/Vref; //1-Qles in order to have the best quality equal to 1

  /* Prevent calities > 1 due to the machin epsilon */
  cal = MG_MIN (1., cal);

  // the normalization by ALPHAD
  //in order to be coherent with the other quality measure
  return cal/MMG3D_ALPHAD;
}

/**
 * \param a pointer to the coor of the first tetra vertex.
 * \param b pointer to the coor of the second tetra vertex.
 * \param c pointer to the coor of the third tetra vertex.
 * \param d pointer to the coor of the fourth tetra vertex.
 * \return The isotropic quality of the tet.
 *
 * Compute the quality of a tetra given by 4 points a,b,c,d with respect to the
 * isotropic metric \a met.
 *
 */
static
inline double MMG5_caltet_iso_4pt(double *a, double *b, double *c, double *d) {
  double       abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz;
  double       cdx,cdy,cdz;
  double       vol,v1,v2,v3,rap;

  /* volume: (ac^ad).ab/6. Note that here we compute 6*vol. */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  rap = abx*abx + aby*aby + abz*abz;

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  rap += acx*acx + acy*acy + acz*acz;

  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  rap += adx*adx + ady*ady + adz*adz;

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol < MMG5_EPSD2 )  return 0.0;

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];
  rap += bcx*bcx + bcy*bcy + bcz*bcz;

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];
  rap += bdx*bdx + bdy*bdy + bdz*bdz;

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];
  rap += cdx*cdx + cdy*cdy + cdz*cdz;
  if ( rap < MMG5_EPSD2 )  return 0.0;

  /* quality = 6*vol / len^3/2. Q = 1/(12 sqrt(3)) for the regular tetra of length 1. */
  rap = rap * sqrt(rap);
  return vol / rap;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param pt pointer to a tetrahedra.
 * \return The isotropic quality of the tet.
 *
 * Compute the quality of the tet pt with respect to the isotropic metric \a
 * met.
 *
 */
static
inline double MMG5_caltet_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra  pt) {
  double       *a,*b,*c,*d;
  MMG5_int     ia, ib, ic, id;

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];

  a = mesh->point[ia].c;
  b = mesh->point[ib].c;
  c = mesh->point[ic].c;
  d = mesh->point[id].c;

  return MMG5_caltet_iso_4pt(a,b,c,d);

}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the meric structure.
 * \param pt pointer to a tetrahedra.
 * \return The anisotropic quality of the tet or 0.0 if fail.
 *
 * Compute the quality of the tet pt with respect to the anisotropic metric \a
 * met. \f$ Q = V_met(K) / (sum(len(edge_K)^2)^(3/2) \f$.
 *
 * \todo test with the square of this measure
 */
static
inline double MMG5_caltet_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt) {
  double       cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double       h1,h2,h3,h4,h5,h6,det,vol,rap,v1,v2,v3,num;
  double       bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double       *a,*b,*c,*d;
  double       mm[6];
  MMG5_int     ip[4];

  ip[0] = pt->v[0];
  ip[1] = pt->v[1];
  ip[2] = pt->v[2];
  ip[3] = pt->v[3];

  /* average metric */
  if ( !MMG5_moymet(mesh,met,pt,&mm[0]) )
    return (0.0);

  a = mesh->point[ip[0]].c;
  b = mesh->point[ip[1]].c;
  c = mesh->point[ip[2]].c;
  d = mesh->point[ip[3]].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol <= 0. )  return 0.0;

  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
      - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
      + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
  if ( det < MMG5_EPSD2 )   {
    return 0.0;
  }
  det = sqrt(det) * vol;

  /* edge lengths */
  h1 = mm[0]*abx*abx + mm[3]*aby*aby + mm[5]*abz*abz
    + 2.0*(mm[1]*abx*aby + mm[2]*abx*abz + mm[4]*aby*abz);
  h2 =  mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz
    + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);
  h3 = mm[0]*adx*adx + mm[3]*ady*ady + mm[5]*adz*adz
    + 2.0*(mm[1]*adx*ady + mm[2]*adx*adz + mm[4]*ady*adz);
  h4 =  mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz
    + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);
  h5 =  mm[0]*bdx*bdx + mm[3]*bdy*bdy + mm[5]*bdz*bdz
    + 2.0*(mm[1]*bdx*bdy + mm[2]*bdx*bdz + mm[4]*bdy*bdz);
  h6 =  mm[0]*cdx*cdx + mm[3]*cdy*cdy + mm[5]*cdz*cdz
    + 2.0*(mm[1]*cdx*cdy + mm[2]*cdx*cdz + mm[4]*cdy*cdz);

  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;

  num = sqrt(rap) * rap;

  cal = det / num;

  return cal;
}

#endif
