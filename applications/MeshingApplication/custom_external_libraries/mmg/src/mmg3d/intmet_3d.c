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
 * \file mmg3d/intmet_3d.c
 * \brief Metric interpolations.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmg3d_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param ip global index of the new point in which we want to compute the metric.
 * \param s interpolation parameter (between 0 and 1).
 * \return 0 if fail, 1 otherwise.
 *
 * Interpolation of anisotropic sizemap at parameter \a s along edge \a i of elt
 * \a k for a special storage of ridges metric (after defsiz call).
 *
 */
int MMG5_intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,
                      double s) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  MMG5_pxPoint  pxp;
  double        *m;
  MMG5_int      ip1,ip2;

  pt = &mesh->tetra[k];
  m  = &met->m[6*ip];
  ip1 = pt->v[MMG5_iare[i][0]];
  ip2 = pt->v[MMG5_iare[i][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( pxt->tag[i] & MG_GEO && !(pxt->tag[i] & MG_NOM)  ) {
      ppt = &mesh->point[ip];
      assert(ppt->xp);
      pxp = &mesh->xpoint[ppt->xp];
      return MMG5_intridmet(mesh,met,ip1,ip2,s,pxp->n1,m);
    }
    else if ( pxt->tag[i] & MG_BDY ) {
     return MMG5_intregmet(mesh,met,k,i,s,m);
    }
    else {
      /* The edge is an internal edge. */
      return MMG5_intvolmet(mesh,met,k,i,s,m);
    }
  }
  else {
    /* The edge is an internal edge. */
    return MMG5_intvolmet(mesh,met,k,i,s,m);
  }
  return 0;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param ip global index of the new point in which we want to compute the metric.
 * \param s interpolation parameter (between 0 and 1).
 * \return 0 if fail, 1 otherwise.
 *
 * Interpolation of anisotropic sizemap at parameter \a s along edge \a i of elt
 * \a k for a classic storage of ridges metrics (before defsiz call).
 *
 */
int MMG3D_intmet33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,
                      double s) {
  MMG5_pTetra   pt;
  double        *m,*n,*mr;
  MMG5_int      ip1,ip2;

  pt = &mesh->tetra[k];
  ip1 = pt->v[MMG5_iare[i][0]];
  ip2 = pt->v[MMG5_iare[i][1]];

  m   = &met->m[6*ip1];
  n   = &met->m[6*ip2];
  mr  = &met->m[6*ip];

  return MMG5_mmgIntmet33_ani(m,n,mr,s);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param ip global index of the new point in which we want to compute the metric.
 * \param s interpolation parameter (between 0 and 1).
 * \return 0 if fail, 1 otherwise.
 *
 * Interpolation of anisotropic sizemap at parameter \a s along edge \a i of elt
 * \a k.
 *
 */
int MMG5_intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int ip,
                      double s) {
  MMG5_pTetra   pt;
  MMG5_int      ip1, ip2;
  double        *m1,*m2,*mm;

  pt = &mesh->tetra[k];
  ip1 = pt->v[MMG5_iare[i][0]];
  ip2 = pt->v[MMG5_iare[i][1]];

  m1 = &met->m[met->size*ip1];
  m2 = &met->m[met->size*ip2];
  mm = &met->m[met->size*ip];

  return MMG5_interp_iso(m1,m2,mm,s);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param s interpolation parameter.
 * \param mr computed metric.
 * \return  0 if fail, 1 otherwise.
 *
 * Metric interpolation on edge \a i in elt \a it at
 * parameter \f$ 0 <= s0 <= 1 \f$ from \a p1 result is stored in \a mr. edge
 * \f$ p_1-p_2 \f$ must not be a ridge.
 *
 * */
int MMG5_intregmet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,double s,
                    double mr[6]) {
  MMG5_pTetra     pt;
  MMG5_pxTetra    pxt;
  MMG5_Tria       ptt;
  int             ifa0, ifa1, iloc, ier;

  pt   = &mesh->tetra[k];
  pxt  = &mesh->xtetra[pt->xt];
  ifa0 = MMG5_ifar[i][0];
  ifa1 = MMG5_ifar[i][1];

  ier  = -1;

  if ( pxt->ftag[ifa0] & MG_BDY ) {
    MMG5_tet2tri( mesh,k,ifa0,&ptt);
    iloc = MMG5_iarfinv[ifa0][i];
    assert(iloc >= 0);
    ier = MMG5_interpreg_ani(mesh,met,&ptt,iloc,s,mr);
  }
  else if ( pxt->ftag[ifa1] & MG_BDY ) {
    MMG5_tet2tri( mesh,k,ifa1,&ptt);
    iloc = MMG5_iarfinv[ifa1][i];
    assert(iloc >= 0);
    ier = MMG5_interpreg_ani(mesh,met,&ptt,iloc,s,mr);
  }

  /* if ier = -1, then i is a boundary edge but the tet has no bdy face. Don't
   * do anything, the edge will be split via a boundary tetra. Otherwise, if
   * ier=0, interpreg_ani has failed, if ier=1, interpreg_ani succeed. */
  if ( mesh->info.ddebug && !ier ) {
    fprintf(stderr," %s: %d: interpreg_ani error.\n",__func__,__LINE__);
    fprintf(stderr," Elt %"MMG5_PRId": %"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId"\n",
            MMG3D_indElt(mesh,k),MMG3D_indPt(mesh,pt->v[0]),
            MMG3D_indPt(mesh,pt->v[1]),MMG3D_indPt(mesh,pt->v[2]),
            MMG3D_indPt(mesh,pt->v[3]));
  }

  return ier;
}


/**
 * \param ma pointer on a metric
 * \param mb pointer on a metric
 * \param mp pointer on the computed interpolated metric
 * \param t interpolation parameter (comprise between 0 and 1)
 * \return 1 if success, 0 if fail.
 *
 * Linear interpolation of anisotropic sizemap along an internal edge.
 *
 */
static inline int
MMG5_intregvolmet(double *ma,double *mb,double *mp,double t) {
  double        dma[6],dmb[6],mai[6],mbi[6],mi[6];
  int           i;
  static int8_t mmgWarn=0;

  for (i=0; i<6; i++) {
    dma[i] = ma[i];
    dmb[i] = mb[i];
  }
  if ( !MMG5_invmat(dma,mai) || !MMG5_invmat(dmb,mbi) ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 invalid metric.\n",__func__);
    }
    return 0;
  }
  for (i=0; i<6; i++)
    mi[i] = (1.0-t)*mai[i] + t*mbi[i];

  if ( !MMG5_invmat(mi,mai) ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 invalid metric.\n",__func__);
    }
    return 0;
  }

  for (i=0; i<6; i++)  mp[i] = mai[i];
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param i local index of edge in \a k.
 * \param s interpolation parameter.
 * \param mr computed metric.
 * \return  0 if fail, 1 otherwise.
 *
 * Metric interpolation on edge \a i in elt \a it at
 * parameter \f$ 0 <= s0 <= 1 \f$ from \a p1 result is stored in \a mr. edge
 * \f$ p_1-p_2 \f$ is an internal edge.
 *
 */
int MMG5_intvolmet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,double s,
                    double mr[6]) {
  MMG5_pTetra     pt;
  MMG5_pPoint     pp1,pp2;
  double          m1[6],m2[6];
  MMG5_int        ip1,ip2;
  int             l,ier;

  pt  = &mesh->tetra[k];

  ip1 = pt->v[MMG5_iare[i][0]];
  ip2 = pt->v[MMG5_iare[i][1]];

  pp1 = &mesh->point[ip1];
  pp2 = &mesh->point[ip2];

  // build metric at ma and mb points
  if ( MG_RID(pp1->tag) ) {
    if (!MMG5_moymet(mesh,met,pt,m1)) return 0;
  } else {
    for ( l=0; l<6; ++l )
      m1[l] = met->m[6*ip1+l];
  }
  if ( MG_RID(pp2->tag) ) {
    if (!MMG5_moymet(mesh,met,pt,m2)) return 0;
  } else {
    for ( l=0; l<6; ++l )
      m2[l] = met->m[6*ip2+l];
  }

  ier = MMG5_intregvolmet(m1,m2,mr,s);
  if ( mesh->info.ddebug && ( (!ier) || (fabs(mr[5]) < 1e-6) ) ) {
    fprintf(stderr,"  ## Error: %s:\n",__func__);
    fprintf(stderr,"            pp1 : %d %d \n",
            MG_SIN(pp1->tag) || (MG_NOM & pp1->tag),pp1->tag & MG_GEO);
    fprintf(stderr,"            m1 %e %e %e %e %e %e\n",
            m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]);
    fprintf(stderr,"            pp2 : %d %d \n",
            MG_SIN(pp2->tag) || (MG_NOM & pp2->tag),pp2->tag & MG_GEO);
    fprintf(stderr,"            m2 %e %e %e %e %e %e\n",
            m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]);
    fprintf(stderr,"            mr %e %e %e %e %e %e\n",
            mr[0],mr[1],mr[2],mr[3],mr[4],mr[5]);
    return 0;
  }

  return 1;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the tetra.
 * \param ip index of the point on which we compute the metric.
 * \param cb barycentric coordinates of \a ip in \a k.
 * \return 1.
 *
 * Linear interpolation of isotropic sizemap in a tetra given the barycentric
 * coordinates of the new point in \a k.
 *
 */
int MMG5_interp4bar_iso(MMG5_pMesh mesh, MMG5_pSol met, MMG5_int k, MMG5_int ip,
                         double cb[4]) {
  MMG5_pTetra pt;

  pt = &mesh->tetra[k];

  met->m[ip] = cb[0]*met->m[pt->v[0]]+cb[1]*met->m[pt->v[1]] +
    cb[2]*met->m[pt->v[2]]+cb[3]*met->m[pt->v[3]];

  return 1;

}

/**
 * \param met pointer toward the metric structure.
 * \param ip index of the point on which we compute the metric.
 * \param cb barycentric coordinates of \a ip in the tetra.
 * \param dm0 metric of the first vertex of the tet.
 * \param dm1 metric of the second vertex of the tet.
 * \param dm2 metric of the third vertex of the tet.
 * \param dm3 metric of the fourth vertex of the tet.
 * \return 1 if success, 0 if fail.
 *
 * Linear interpolation of anisotropic sizemap in a tetra given the barycentric
 * coordinates of the new point in a tetra.
 *
 */
static inline
int MMG5_interp4barintern(MMG5_pSol met,MMG5_int ip,double cb[4],double dm0[6],
                           double dm1[6],double dm2[6],double dm3[6]) {
  double        m0i[6],m1i[6],m2i[6],m3i[6],mi[6];
  int           i;
  static int8_t mmgWarn=0;

 if ( !MMG5_invmat(dm0,m0i) || !MMG5_invmat(dm1,m1i) ||
       !MMG5_invmat(dm2,m2i) || !MMG5_invmat(dm3,m3i) ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 invalid metric.\n",__func__);
    }
    return 0;
  }
  for (i=0; i<6; i++)
    mi[i] = cb[0]*m0i[i] + cb[1]*m1i[i] + cb[2]*m2i[i] + cb[3]*m3i[i];

  if ( !MMG5_invmat(mi,m0i) ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 invalid metric.\n",__func__);
    }
    return 0;
  }

  for (i=0; i<6; i++)  met->m[met->size*ip+i] = m0i[i];

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the tetra.
 * \param ip index of the point on which we compute the metric.
 * \param cb barycentric coordinates of \a ip in \a k.
 * \return 1 if success, 0 if fail.
 *
 * Linear interpolation of anisotropic sizemap in a tetra given the barycentric
 * coordinates of the new point in \a k.
 *
 */
int MMG5_interp4bar_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG5_int k, MMG5_int ip,
                         double cb[4]) {
  MMG5_pTetra   pt;
  MMG5_pPoint   pp1,pp2,pp3,pp4;
  double        dm0[6],dm1[6],dm2[6],dm3[6];
  int           i;

  pt  = &mesh->tetra[k];
  pp1 = &mesh->point[pt->v[0]];
  if ( MG_SIN_OR_NOM(pp1->tag) ) {
    for (i=0; i<6; i++) {
      dm0[i] = met->m[met->size*pt->v[0]+i];
    }
  } else if(pp1->tag & MG_GEO) {
    if (!MMG5_moymet(mesh,met,pt,&dm0[0])) return 0;
  } else{
    for (i=0; i<6; i++) {
      dm0[i] = met->m[met->size*pt->v[0]+i];
    }
  }
  pp2 = &mesh->point[pt->v[1]];
  if ( MG_SIN_OR_NOM(pp2->tag) ) {
    for (i=0; i<6; i++) {
      dm1[i] = met->m[met->size*pt->v[1]+i];
    }
  } else if(pp2->tag & MG_GEO) {
    if (!MMG5_moymet(mesh,met,pt,&dm1[0])) return 0;
  } else{
    for (i=0; i<6; i++) {
      dm1[i] = met->m[met->size*pt->v[1]+i];
    }
  }
  pp3 = &mesh->point[pt->v[2]];
  if ( MG_SIN_OR_NOM(pp3->tag) ) {
    for (i=0; i<6; i++) {
      dm2[i] = met->m[met->size*pt->v[2]+i];
    }
  } else if(pp3->tag & MG_GEO) {
    if (!MMG5_moymet(mesh,met,pt,&dm2[0])) return 0;
  } else{
    for (i=0; i<6; i++) {
      dm2[i] = met->m[met->size*pt->v[2]+i];
    }
  }
  pp4 = &mesh->point[pt->v[3]];
  if ( MG_SIN_OR_NOM(pp4->tag) ) {
    for (i=0; i<6; i++) {
      dm3[i] = met->m[met->size*pt->v[3]+i];
    }
  } else if(pp4->tag & MG_GEO) {
    if (!MMG5_moymet(mesh,met,pt,&dm3[0])) return 0;
  } else{
    for (i=0; i<6; i++) {
      dm3[i] = met->m[met->size*pt->v[3]+i];
    }
  }

  return MMG5_interp4barintern(met,ip,cb,dm0,dm1,dm2,dm3);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the tetra.
 * \param ip index of the point on which we compute the metric.
 * \param cb barycentric coordinates of \a ip in \a k.
 * \return 1 if success, 0 if fail.
 *
 * Linear interpolation of anisotropic sizemap in a tetra given the barycentric
 * coordinates of the new point in \a k.
 *
 */
int MMG5_interp4bar33_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG5_int k, MMG5_int ip,
                           double cb[4]) {
  MMG5_pTetra   pt;
  double        dm0[6],dm1[6],dm2[6],dm3[6];
  int           i;

  pt  = &mesh->tetra[k];
  for (i=0; i<6; i++) {
    dm0[i] = met->m[met->size*pt->v[0]+i];
  }

  for (i=0; i<6; i++) {
    dm1[i] = met->m[met->size*pt->v[1]+i];
  }

  for (i=0; i<6; i++) {
    dm2[i] = met->m[met->size*pt->v[2]+i];
  }

  for (i=0; i<6; i++) {
    dm3[i] = met->m[met->size*pt->v[3]+i];
  }
  return MMG5_interp4barintern(met,ip,cb,dm0,dm1,dm2,dm3);
}
