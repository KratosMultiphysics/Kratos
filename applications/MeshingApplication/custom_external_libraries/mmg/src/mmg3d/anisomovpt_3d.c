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
 * \file mmg3d/anisomovpt_3d.c
 * \brief Functions to move a point in the mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmg3d.h"
#include "inlined_functions_3d_private.h"
#include "mmg3dexterns_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param list pointer to the volumic ball of the point.
 * \param ilist size of the volumic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.
 *
 * \return 0 if we can't move the point, 1 if we can.
 *
 * Move internal point whose volumic is passed.
 *
 * \remark the metric is not interpolated at the new position.
 * \remark we don't check if we break the hausdorff criterion.
 *
 */
int MMG5_movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met, MMG3D_pPROctree PROctree, int64_t *list,int ilist,
                       int improve) {


  MMG5_pTetra          pt,pt0;
  MMG5_pPoint          p0,p1,p2,p3,ppt0;
  double               vol,totvol,m[6];
  double               calold,calnew,callist[MMG3D_LMAX+2],det;
  MMG5_int             iel;
  int                  k,i0;

  assert ( ilist > 0 );
  if ( ilist <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s:"
            " volumic ball has null or negative size (%d)\n",
            __func__,ilist);
    return 0;
  }

  pt0    = &mesh->tetra[0];
  ppt0   = &mesh->point[0];
  memset(ppt0,0,sizeof(MMG5_Point));

  if ( met->m ) {
    iel = list[0] / 4;
    i0  = list[0] % 4;
    memcpy(&met->m[0],&met->m[met->size*mesh->tetra[iel].v[i0]],met->size*sizeof(double));
  }

  /* Coordinates of optimal point */
  calold = DBL_MAX;
  totvol = 0.0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    pt = &mesh->tetra[iel];
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];
    p3 = &mesh->point[pt->v[3]];
    vol= MMG5_det4pt(p0->c,p1->c,p2->c,p3->c);

    if ( !MMG5_moymet(mesh,met,pt,m) ) {
      // MMG5_moymet must succeed because we have at least 1 point of the tet
      // that is internal.
      return 0;
    }

    det = m[0] * ( m[3]*m[5] - m[4]*m[4]) - m[1] * ( m[1]*m[5] - m[2]*m[4])
      + m[2] * ( m[1]*m[4] - m[2]*m[3]);
    if ( det < MMG5_EPSD2 ) {
      return 0;
    }

    vol *= sqrt(det);

    totvol += vol;
    /* barycenter */
    ppt0->c[0] += 0.25 * vol*(p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]);
    ppt0->c[1] += 0.25 * vol*(p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]);
    ppt0->c[2] += 0.25 * vol*(p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]);
    calold = MG_MIN(calold, pt->qual);
  }
  if (totvol < MMG5_EPSD2) {
    return 0;
  }

  totvol = 1.0 / totvol;
  ppt0->c[0] *= totvol;
  ppt0->c[1] *= totvol;
  ppt0->c[2] *= totvol;

  /* Check new position validity */
  calnew = DBL_MAX;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    i0  = list[k] % 4;
    pt  = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    callist[k] = MMG5_orcal(mesh,met,0);
    if (callist[k] < MMG5_NULKAL) {
      return 0;
    }
    calnew = MG_MIN(calnew,callist[k]);
  }
  if (calold < MMG5_EPSOK && calnew <= calold) {
    return 0;
  }
  else if (calnew < MMG5_EPSOK) {
    return 0;
  }
  else if ( improve && calnew < 1.02* calold ) {
    return 0;
  }
  else if ( calnew < 0.3 * calold ) {
    return 0;
  }

  /* update position */
  if ( PROctree )
    MMG3D_movePROctree(mesh, PROctree, pt->v[i0], ppt0->c, p0->c);

  p0 = &mesh->point[pt->v[i0]];
  p0->c[0] = ppt0->c[0];
  p0->c[1] = ppt0->c[1];
  p0->c[2] = ppt0->c[2];
  for (k=0; k<ilist; k++) {
    (&mesh->tetra[list[k]/4])->qual=callist[k];
    (&mesh->tetra[list[k]/4])->mark=mesh->mark;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param listv pointer to the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer to the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.

 * \return 0 if we can't move the point, 1 if we can, -1 if we fail.
 *
 * \remark we don't check if we break the hausdorff criterion.
 * \remark the metric is not interpolated at the new position.
 *
 * Move boundary regular point, whose volumic and surfacic balls are passed.
 *
 */
int MMG5_movbdyregpt_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree,
                         int64_t *listv,int ilistv,MMG5_int *lists,int ilists,
                         int improveSurf, int improveVol) {
  MMG5_pTetra       pt,pt0;
  MMG5_pxTetra      pxt;
  MMG5_pPoint       p0;
  MMG5_Tria         tt;
  MMG5_pxPoint      pxp;
  MMG5_Bezier       b;
  double            n[3],r[3][3],lispoi[3*MMG3D_LMAX+1],det2d;
  double            detloc,gv[2],step,lambda[3];
  double            o[3],no[3],*m0,ncur[3],nprev[3],nneighi[3];
  double            calold,calnew,caltmp,callist[MMG3D_LMAX+2];
  MMG5_int          k,kel,iel,ip0,nxp;
  int               ier,l;
  uint8_t           i0,iface,i;
  static int        warn = 0;
  static int8_t     mmgErr0=0;

  step = 0.1;
  if ( ilists < 2 )      return 0;

  k  = listv[0] / 4;
  i0 = listv[0] % 4;
  pt = &mesh->tetra[k];
  ip0 = pt->v[i0];
  p0 = &mesh->point[ip0];
  m0 = &met->m[6*ip0];
  assert( p0->xp && !MG_EDG(p0->tag) );

  n[0] = mesh->xpoint[p0->xp].n1[0];
  n[1] = mesh->xpoint[p0->xp].n1[1];
  n[2] = mesh->xpoint[p0->xp].n1[2];

  /** Step 1 : rotation matrix that sends normal n to the third coordinate vector of R^3 */
  if ( !MMG5_rotmatrix(n,r) ) {
    return 0;
  }

  /** Step 2 : rotation of the oriented surfacic ball with r : lispoi[k] is the common edge
      between faces lists[k-1] and lists[k] */
  if ( !MMG3D_rotate_surfacicBall(mesh,lists,ilists,ip0,r,lispoi) ) {
    return 0;
  }

  /** Step 3 :  Compute gradient towards optimal position = centre of mass of the
      ball, projected to tangent plane */
  gv[0] = 0.0;
  gv[1] = 0.0;

  for (k=0; k<ilists; k++) {
    iel    = lists[k] / 4;
    iface  = lists[k] % 4;
    pt     = &mesh->tetra[iel];
    pxt    = &mesh->xtetra[pt->xt];

    assert( 0<=iface && iface<4 && "unexpected local face idx");
    MMG5_tet2tri(mesh,iel,iface,&tt);

    if(!MMG5_bezierCP(mesh,&tt,&b,MG_GET(pxt->ori,iface))){
      if( !mmgErr0 ) {
        mmgErr0 = 1;
        fprintf(stderr,"\n  ## Error: %s: function MMG5_bezierCP return 0.\n",
                __func__);
      }
      return -1;
    }

    /* Compute integral of sqrt(T^J(xi)  M(P(xi)) J(xi)) * P(xi) over the triangle */
    if ( !MMG5_elementWeight(mesh,met,&tt,p0,&b,r,gv) ) {
      if ( !warn ) {
        ++warn;
        fprintf(stderr,"\n  ## Warning: %s:"
                " unable to compute optimal position for at least"
                " 1 point.\n",__func__);
      }
      return 0;
    }
  }

  /* At this point : gv = - gradient of V = direction to follow */
  /** Step 4 : locate new point in the ball, and compute its barycentric coordinates */
  det2d = lispoi[1]*gv[1] - lispoi[2]*gv[0];
  kel = 0;
  if ( det2d >= 0.0 ) {
    for (k=0; k<ilists; k++) {
      detloc = gv[0]*lispoi[3*(k+1)+2] - gv[1]*lispoi[3*(k+1)+1];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == ilists ) {
      return 0;
    }
  }
  else {
    for (k=ilists-1; k>=0; k--) {
      detloc = lispoi[3*k+1]*gv[1] - lispoi[3*k+2]*gv[0];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == -1 ) {
      return 0;
    }
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  det2d = -gv[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) + \
    gv[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(det2d) < MMG5_EPSD2 ) {
    return 0;
  }

  det2d = 1.0 / det2d;
  step *= det2d;

  det2d = lispoi[3*(kel)+1]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]) - \
    lispoi[3*(kel)+2 ]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]);
  step *= det2d;
  step  = fabs(step);
  gv[0] *= step;
  gv[1] *= step;

  /* Computation of the barycentric coordinates of the new point in the corresponding triangle. */
  det2d = lispoi[3*kel+1]*lispoi[3*(kel+1)+2] - lispoi[3*kel+2]*lispoi[3*(kel+1)+1];
  if ( det2d < MMG5_EPSD2 ) {
    return 0;
  }
  det2d = 1.0 / det2d;
  lambda[1] = lispoi[3*(kel+1)+2]*gv[0] - lispoi[3*(kel+1)+1]*gv[1];
  lambda[2] = -lispoi[3*(kel)+2]*gv[0] + lispoi[3*(kel)+1]*gv[1];
  lambda[1]*= (det2d);
  lambda[2]*= (det2d);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /** Step 5 : come back to original problem, compute patch in triangle iel and
   * check that geometric approx has not been degraded too much */
  // Remark: if we call following function with a pointer for n, we have to set
  // the pointer again after the function call as it may invalidate it if it
  // reallocates the xpoint array
  nxp = MMG3D_movbdyregpt_geom(mesh,lists,kel,ip0,n,lambda,o,no);
  if ( nxp < 0 ) {
    return -1;
  }
  else if ( !nxp ) {
    return 0;
  }
  pxp = &mesh->xpoint[nxp];

  /* parallel transport of metric at p0 to new point. */
  if ( !MMG5_paratmet(p0->c,n,m0,o,no,&met->m[0]) ) {
    return 0;
  }

  /* For each surfacic triangle build a virtual displaced triangle for check
   * purposes :
   *      - check the new triangle qualities;
   *      - check normal deviation with the adjacent through the edge facing ip0
   *        and the previous one */
  k           = lists[ilists-1] / 4;
  iface       = lists[ilists-1] % 4;

  assert( 0<=iface && iface<4 && "unexpected local face idx");
  MMG5_tet2tri(mesh,k,iface,&tt);
  for( i=0 ; i<3 ; i++ )
    if ( tt.v[i] == ip0 )      break;
  assert ( i<3 );
  if ( i>=3 ) return 0;
  tt.v[i] = 0;

  if ( !MMG5_nortri(mesh, &tt, nprev) ) return 0;

  calold = calnew = DBL_MAX;
  for (l=0; l<ilists; l++) {
    k     = lists[l] / 4;
    iface = lists[l] % 4;

    assert( 0<=iface && iface<4 && "unexpected local face idx");
    MMG5_tet2tri(mesh,k,iface,&tt);
    calold = MG_MIN(calold,MMG5_caltri(mesh,met,&tt));

    for( i=0 ; i<3 ; i++ )
      if ( tt.v[i] == ip0 )      break;

    assert ( i<3 );
    if ( i>=3 ) return 0;
    tt.v[i] = 0;

    caltmp = MMG5_caltri(mesh,met,&tt);
    if ( caltmp < MMG5_EPSD2 ) {
      /* We don't check the input triangle qualities, thus we may have a very
       * bad triangle in our mesh */
      return 0;
    }
    calnew = MG_MIN(calnew,caltmp);

    if ( !MMG5_nortri(mesh, &tt, ncur) ) return 0;

    if ( !MG_GEO_OR_NOM(tt.tag[i]) ) {
      /* Check normal deviation between k and the triangle facing ip0 */
      /* Note that we don't check the ridge creation along non-manifold edges
       * because nm edges are skipped while we add ridge tags during analysis
       * step. */
      ier = MMG3D_normalAdjaTri(mesh,k,iface, i,nneighi);
      if ( ier <= 0 ) {
        return 0;
      }
      ier =  MMG5_devangle( ncur, nneighi, mesh->info.dhd );
      if ( ier <= 0 ) {
        return 0;
      }
    }

    i = MMG5_iprv2[i];

    if ( !MG_GEO_OR_NOM(tt.tag[i]) ) {
      /* Check normal deviation between k and the previous triangle */
      ier =  MMG5_devangle( ncur, nprev, mesh->info.dhd );
      if ( ier<=0 ) {
        return 0;
      }
    }
    memcpy(nprev, ncur, 3*sizeof(double));

  }
  if ( calold < MMG5_EPSOK && calnew <= calold ) {
    return 0;
  }  else if (calnew < MMG5_EPSOK) {
    return 0;
  }
  else if (improveSurf && calnew < 1.02*calold) {
    return 0;
  }
  else if ( calnew < 0.3*calold ) {
    return 0;
  }
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    k    = listv[l] / 4;
    i0 = listv[l] % 4;
    pt = &mesh->tetra[k];
    pt0 = &mesh->tetra[0];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    calold = MG_MIN(calold, pt->qual);
    callist[l]=MMG5_orcal(mesh,met,0);
    if ( callist[l] < MMG5_NULKAL )  {
      return 0;
    }
    calnew = MG_MIN(calnew,callist[l]);
  }

  if ( calold < MMG5_EPSOK && calnew <= calold ) {
    return 0;
  }
  else if (calnew < MMG5_EPSOK) {
    return 0;
  }
  else if (improveVol && calnew < calold) {
    return 0;
  }
  else if ( calnew < 0.3*calold ) {
    return 0;
  }

  /* When all tests have been carried out, update coordinates, normals and metrics*/
  if ( PROctree )
    MMG3D_movePROctree(mesh, PROctree, ip0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  n[0] = no[0];
  n[1] = no[1];
  n[2] = no[2];

  memcpy(m0,&met->m[0],6*sizeof(double));

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual= callist[l];
    (&mesh->tetra[listv[l]/4])->mark= mesh->mark;
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param listv pointer to the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer to the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.
 * \param edgTag MG_REF, MG_NOM or MG_GEO depending on type of edge on which we move.
 *
 * \return 0 if fail, 1 if success.
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 * Move boundary reference, ridge or non-manifold point, whose volumic and
 * surfacic balls are passed.
 *
 * \todo Factorization of this function with the iso one.
 */
static inline
int MMG3D_movbdycurvept_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree,int64_t *listv,
                            int ilistv,MMG5_int *lists, int ilists,int improve,const uint16_t edgTag){
  MMG5_pTetra           pt;
  MMG5_pPoint           p0;
  MMG5_Tria             tt;
  double                ll1old,ll2old,l1new,l2new;
  double                o[3],no[3],no2[3],to[3], ncur[3],nprev[3],nneighi[3];
  double                calold,calnew,caltmp;
  MMG5_int              l,iel,ip0,ip1,ip2,ip;
  uint8_t               i,iface,isrid;

  ip1 = ip2 = 0;
  pt        = &mesh->tetra[listv[0]/4];
  ip0       = pt->v[listv[0]%4];
  p0        = &mesh->point[ip0];

  /** Step 0: Compute if the edge is a simple ridge to know if we have to
   * compute a second normal at point */
  isrid     = ((MG_GEO & edgTag) && !(MG_NOM & edgTag));

  assert ( edgTag & p0->tag );

  /** Step 1: Travel surfacic ball and recover the two ending points of curve
     (that will be stored in \a ip1 and \a ip2) */
  int ier = MMG3D_curveEndingPts(mesh,lists,ilists,edgTag,ip0,&ip1,&ip2);
  if ( !ier ) {
    return 0;
  }

  /** Step 2: At this point, we get the point extremities of the curve passing
     through ip0 : ip1, ip2, along with support tets it1,it2, the surface faces
     iface1,iface2, and the associated edges ie1,ie2. Computation of
     displacement along curve and checks */

  /** a. Changes needed for choice of time step : see manuscript notes */
  ll1old = MMG5_lenSurfEdg(mesh,met,ip0,ip1,isrid);
  ll2old = MMG5_lenSurfEdg(mesh,met,ip0,ip2,isrid);

  if ( (!ll1old) || (!ll2old) ) {
    return 0;
  }

  /** b. Check sense of displacement, compute support of the associated edge and
   * features of the new position */
  ip = MMG3D_movbdycurvept_newPosForSimu( mesh,p0,ip0,ip1,ip2,ll1old,ll2old,
                                          isrid,MMG3D_MOVSTEP,o,no,no2,to,edgTag );
  if ( !ip ) {
    return 0;
  }

  /** Metric interpolation */
  if ( (MG_GEO & edgTag) && !(MG_NOM & edgTag) ) {
    /* Interpolation of metric between ip0 and ip2 along ridge */
    if ( !MMG5_intridmet(mesh,met,ip0,ip,MMG3D_MOVSTEP,no,&met->m[0]) ) {
      return 0;
    }
  }
  else {
    /* Interpolation of metric between ip0 and ip2 along non manifold or ref edge */
    if ( !MMG5_paratmet(p0->c,mesh->xpoint[p0->xp].n1,&met->m[6*ip0],o,no,&met->m[0]) ) {
      return 0;
    }
  }

  /** Check whether proposed move is admissible under consideration of distances */
  l1new = MMG5_lenSurfEdg(mesh,met,0,ip1,isrid);
  l2new = MMG5_lenSurfEdg(mesh,met,0,ip2,isrid);

  if ( (!l1new) || (!l2new) ) {
    return 0;
  }

  if ( fabs(l2new -l1new) >= fabs(ll2old -ll1old) ) {
    return 0;
  }

  /** For each surfacic triangle build a virtual displaced triangle for check
   * purposes :
   *      - check the new triangle qualities;
   *      - check normal deviation with the adjacent through the edge facing ip0
   *        and the previous one */
  iel         = lists[ilists-1] / 4;
  iface       = lists[ilists-1] % 4;

  assert( 0<=iface && iface<4 && "unexpected local face idx");
  MMG5_tet2tri(mesh,iel,iface,&tt);
  for( i=0 ; i<3 ; i++ ) {
    if ( tt.v[i] == ip0 ) {
      break;
    }
  }

  assert ( i<3 );
  if ( i>=3 ) {
    return 0;
  }
  tt.v[i] = 0;

  if ( !MMG5_nortri(mesh, &tt, nprev) ) {
    return 0;
  }

  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilists ; l++ ){
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;

    assert( 0<=iface && iface<4 && "unexpected local face idx");
    MMG5_tet2tri(mesh,iel,iface,&tt);
    caltmp = MMG5_caltri(mesh,met,&tt);
    calold = MG_MIN(calold,caltmp);

    for( i=0 ; i<3 ; i++ ) {
      if ( tt.v[i] == ip0 ) {
        break;
      }
    }
    assert(i<3);
    if ( i==3 ) {
      return 0;
    }

    tt.v[i] = 0;

    caltmp = MMG5_caltri(mesh,met,&tt);
    if ( caltmp < MMG5_EPSD2 ) {
      /* We don't check the input triangle qualities, thus we may have a very
       * bad triangle in our mesh */
      return 0;
    }
    calnew = MG_MIN(calnew,caltmp);

    if ( !MMG5_nortri(mesh, &tt, ncur) ) {
      return 0;
    }

    if ( !MG_GEO_OR_NOM(tt.tag[i]) ) {
      /* Check normal deviation between iel and the triangle facing ip0 */
      /* Note that we don't check the ridge creation along non-manifold edges
       * because nm edges are skipped while we add ridge tags during analysis
       * step. */
      int16_t ier2 = MMG3D_normalAdjaTri(mesh,iel,iface, i,nneighi);
      if ( ier2 <=0 ) {
        return 0;
      }
      ier2 =  MMG5_devangle( ncur, nneighi, mesh->info.dhd );
      if ( ier2 <= 0 ) {
        return 0;
      }
    }

    i = MMG5_iprv2[i];

    if ( !MG_GEO_OR_NOM(tt.tag[i]) ) {
      /* Check normal deviation between k and the previous triangle */
      int16_t ier2 =  MMG5_devangle( ncur, nprev, mesh->info.dhd );
      if ( ier2<=0 ) {
        return 0;
      }
    }
    memcpy(nprev, ncur, 3*sizeof(double));
  }
  if ( calold < MMG5_EPSOK && calnew <= calold ) {
    return 0;
  }
  else if ( calnew < calold ) {
    return 0;
  }
  memset(&mesh->xpoint[mesh->point[0].xp],0,sizeof(MMG5_xPoint));

  /** d. Check whether all volumes remain positive with new position of the
   * point and update coor, normals, tangents and qualities if move is
   * accepted. */
  ier =  MMG3D_movbdycurvept_chckAndUpdate(mesh,met,PROctree,listv,ilistv,
                                           improve,p0,ip0,isrid,o,no,no2,to);

  if ( ier ) {
    memcpy(&met->m[6*ip0],&met->m[0],6*sizeof(double));
  }

  return ier;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param listv pointer to the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer to the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.

 * \return 0 if fail, 1 if success.
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 * Move boundary reference point, whose volumic and surfacic balls are passed.
 *
 */
int MMG5_movbdyrefpt_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree,int64_t *listv,
                          int ilistv, MMG5_int *lists, int ilists,int improve){

  return MMG3D_movbdycurvept_ani(mesh,met,PROctree,listv,ilistv,lists,ilists,improve,MG_REF);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param listv pointer to the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer to the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.
 * \return 0 if fail, 1 if success.
 *
 * Move boundary non manifold point, whose volumic and (exterior)
 * surfacic balls are passed
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 */
int MMG5_movbdynompt_ani(MMG5_pMesh mesh,MMG5_pSol met, MMG3D_pPROctree PROctree,
                         int64_t *listv,int ilistv, MMG5_int *lists, int ilists,
                         int improve){

  return MMG3D_movbdycurvept_ani(mesh,met,PROctree,listv,ilistv,lists,ilists,improve,MG_NOM);
}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param listv pointer to the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer to the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * \return 0 if fail, 1 if success.
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 * Move boundary ridge point, whose volumic and surfacic balls are passed.
 *
 */
int MMG5_movbdyridpt_ani(MMG5_pMesh mesh, MMG5_pSol met, MMG3D_pPROctree PROctree,
                         int64_t *listv,int ilistv,MMG5_int *lists,int ilists,
                         int improve) {

  return MMG3D_movbdycurvept_ani(mesh,met,PROctree,listv,ilistv,lists,ilists,improve,MG_GEO);
}
