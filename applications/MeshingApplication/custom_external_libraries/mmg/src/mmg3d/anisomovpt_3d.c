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

#include "inlined_functions_3d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param list pointer toward the volumic ball of the point.
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
int _MMG5_movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met, _MMG3D_pOctree octree, int *list,int ilist,
                       int improve) {


  MMG5_pTetra          pt,pt0;
  MMG5_pPoint          p0,p1,p2,p3,ppt0;
  double               vol,totvol,m[6];
  double               calold,calnew,*callist,det;
  int                  k,iel,i0;
  // Dynamic alloc for windows comptibility
  _MMG5_SAFE_MALLOC(callist, ilist, double);

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
    vol= _MMG5_det4pt(p0->c,p1->c,p2->c,p3->c);

    if ( !_MMG5_moymet(mesh,met,pt,m) ) {
      // _MMG5_moymet must succeed because we have at least 1 point of the tet
      // that is internal.
      _MMG5_SAFE_FREE(callist);
      return(0);
    }

    det = m[0] * ( m[3]*m[5] - m[4]*m[4]) - m[1] * ( m[1]*m[5] - m[2]*m[4])
      + m[2] * ( m[1]*m[4] - m[2]*m[3]);
    if ( det < _MMG5_EPSOK ) {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }

    vol *= sqrt(det);

    totvol += vol;
    /* barycenter */
    ppt0->c[0] += 0.25 * vol*(p0->c[0] + p1->c[0] + p2->c[0] + p3->c[0]);
    ppt0->c[1] += 0.25 * vol*(p0->c[1] + p1->c[1] + p2->c[1] + p3->c[1]);
    ppt0->c[2] += 0.25 * vol*(p0->c[2] + p1->c[2] + p2->c[2] + p3->c[2]);
    calold = MG_MIN(calold, pt->qual);
  }
  if (totvol < _MMG5_EPSD2) {
    _MMG5_SAFE_FREE(callist);
    return(0);
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
    callist[k] = _MMG5_orcal(mesh,met,0);
    if (callist[k] < _MMG5_EPSD2) {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
    calnew = MG_MIN(calnew,callist[k]);
  }
  if (calold < _MMG5_NULKAL && calnew <= calold) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if (calnew < _MMG5_NULKAL) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if ( improve && calnew < 1.02* calold ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if ( calnew < 0.3 * calold ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* update position */
  if ( octree )
    _MMG3D_moveOctree(mesh, octree, pt->v[i0], ppt0->c, p0->c);

  p0 = &mesh->point[pt->v[i0]];
  p0->c[0] = ppt0->c[0];
  p0->c[1] = ppt0->c[1];
  p0->c[2] = ppt0->c[2];
  for (k=0; k<ilist; k++) {
    (&mesh->tetra[list[k]/4])->qual=callist[k];
  }

  _MMG5_SAFE_FREE(callist);
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 1.02 of the old minimum element quality.

 * \return 0 if we can't move the point, 1 if we can.
 *
 * \remark we don't check if we break the hausdorff criterion.
 * \remark the metric is not interpolated at the new position.
 *
 * Move boundary regular point, whose volumic and surfacic balls are passed.
 *
 */
int _MMG5_movbdyregpt_ani(MMG5_pMesh mesh, MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv,int *lists,int ilists,
                          int improve) {
  MMG5_pTetra       pt,pt0;
  MMG5_pxTetra      pxt;
  MMG5_pPoint       p0,p1,p2,ppt0;
  MMG5_Tria         tt;
  MMG5_pxPoint      pxp;
  _MMG5_Bezier      pb;
  double            *n,r[3][3],lispoi[3*MMG3D_LMAX+1],ux,uy,uz,det2d;
  double            detloc,gv[2],step,lambda[3];
  double            uv[2],o[3],no[3],to[3],*m0;
  double            calold,calnew,caltmp,*callist;
  int               k,kel,iel,l,n0,na,nb,ntempa,ntempb,ntempc,nxp;
  unsigned char     i0,iface,i;
  static int        warn = 0;

  // Dynamic alloc for windows comptibility
  _MMG5_SAFE_MALLOC(callist, ilistv, double);

  step = 0.1;
  if ( ilists < 2 )      return(0);

  k  = listv[0] / 4;
  i0 = listv[0] % 4;
  pt = &mesh->tetra[k];
  n0 = pt->v[i0];
  p0 = &mesh->point[n0];
  m0 = &met->m[6*n0];
  assert( p0->xp && !MG_EDG(p0->tag) );

  n = &(mesh->xpoint[p0->xp].n1[0]);

  /** Step 1 : rotation matrix that sends normal n to the third coordinate vector of R^3 */
  if ( !_MMG5_rotmatrix(n,r) ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /** Step 2 : rotation of the oriented surfacic ball with r : lispoi[k] is the common edge
      between faces lists[k-1] and lists[k] */
  k                     = lists[0] / 4;
  iface = lists[0] % 4;
  pt            = &mesh->tetra[k];
  na = nb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != n0 ) {
      if ( !na )
        na = pt->v[_MMG5_idir[iface][i]];
      else
        nb = pt->v[_MMG5_idir[iface][i]];
    }
  }

  for (l=1; l<ilists; l++) {
    k                   = lists[l] / 4;
    iface = lists[l] % 4;
    pt          = &mesh->tetra[k];
    ntempa = ntempb = 0;
    for (i=0; i<3; i++) {
      if ( pt->v[_MMG5_idir[iface][i]] != n0 ) {
        if ( !ntempa )
          ntempa = pt->v[_MMG5_idir[iface][i]];
        else
          ntempb = pt->v[_MMG5_idir[iface][i]];
      }
    }
    if ( ntempa == na )
      p1 = &mesh->point[na];
    else if ( ntempa == nb )
      p1 = &mesh->point[nb];
    else if ( ntempb == na )
      p1 = &mesh->point[na];
    else {
      assert(ntempb == nb);
      p1 = &mesh->point[nb];
    }
    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*l+1] =      r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*l+2] =      r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*l+3] =      r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

    na = ntempa;
    nb = ntempb;
  }

  /* Finish with point 0 */;
  k      = lists[0] / 4;
  iface  = lists[0] % 4;
  pt     = &mesh->tetra[k];
  ntempa = ntempb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != n0 ) {
      if ( !ntempa )
        ntempa = pt->v[_MMG5_idir[iface][i]];
      else
        ntempb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  if ( ntempa == na )
    p1 = &mesh->point[na];
  else if ( ntempa == nb )
    p1 = &mesh->point[nb];
  else if ( ntempb == na )
    p1 = &mesh->point[na];
  else {
    assert(ntempb == nb);
    p1 = &mesh->point[nb];
  }

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  lispoi[1] =    r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
  lispoi[2] =    r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
  lispoi[3] =    r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

  /* list goes modulo ilist */
  lispoi[3*ilists+1] = lispoi[1];
  lispoi[3*ilists+2] = lispoi[2];
  lispoi[3*ilists+3] = lispoi[3];

  /* At this point, lispoi contains the oriented surface ball of point p0, that has been rotated
     through r, with the convention that triangle l has edges lispoi[l]; lispoi[l+1] */

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 ) {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
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

    _MMG5_tet2tri(mesh,iel,iface,&tt);

    if(!_MMG5_bezierCP(mesh,&tt,&pb,MG_GET(pxt->ori,iface))){
      _MMG5_SAFE_FREE(callist);
      return(0);
    }

    /* Compute integral of sqrt(T^J(xi)  M(P(xi)) J(xi)) * P(xi) over the triangle */
    if ( !_MMG5_elementWeight(mesh,met,&tt,p0,&pb,r,gv) ) {
      if ( !warn ) {
        ++warn;
        printf("  ## Warning: unable to compute optimal position for at least"
               " 1 point.\n" );
      }
      _MMG5_SAFE_FREE(callist);
      return(0);
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
      _MMG5_SAFE_FREE(callist);
      return(0);
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
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  det2d = -gv[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) + \
    gv[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(det2d) < _MMG5_EPSD ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
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
  if ( det2d < _MMG5_EPSD ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  det2d = 1.0 / det2d;
  lambda[1] = lispoi[3*(kel+1)+2]*gv[0] - lispoi[3*(kel+1)+1]*gv[1];
  lambda[2] = -lispoi[3*(kel)+2]*gv[0] + lispoi[3*(kel)+1]*gv[1];
  lambda[1]*= (det2d);
  lambda[2]*= (det2d);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /** Step 5 : come back to original problem, and compute patch in triangle iel */
  iel    = lists[kel] / 4;
  iface  = lists[kel] % 4;
  pt     = &mesh->tetra[iel];
  pxt    = &mesh->xtetra[pt->xt];

  _MMG5_tet2tri(mesh,iel,iface,&tt);

  if(!_MMG5_bezierCP(mesh,&tt,&pb,MG_GET(pxt->ori,iface))){
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* Now, for Bezier interpolation, one should identify which of i,i1,i2 is 0,1,2
     recall uv[0] = barycentric coord associated to pt->v[1], uv[1] associated to pt->v[2] and
     1 - uv[0] - uv[1] is associated to pt->v[0]. For this, use the fact that kel, kel + 1 is
     positively oriented with respect to n */
  na = nb = 0;
  for( i=0 ; i<4 ; i++ ){
    if ( (pt->v[i] != n0) && (pt->v[i] != pt->v[iface]) ) {
      if ( !na )
        na = pt->v[i];
      else
        nb = pt->v[i];
    }
  }
  p1 = &mesh->point[na];
  p2 = &mesh->point[nb];
  detloc = _MMG5_det3pt1vec(p0->c,p1->c,p2->c,n);

  /* ntempa = point to which is associated 1 -uv[0] - uv[1], ntempb = uv[0], ntempc = uv[1] */
  ntempa = pt->v[_MMG5_idir[iface][0]];
  ntempb = pt->v[_MMG5_idir[iface][1]];
  ntempc = pt->v[_MMG5_idir[iface][2]];

  /* na = lispoi[kel] -> lambda[1], nb = lispoi[kel+1] -> lambda[2] */
  if ( detloc > 0.0 ) {
    if ( ntempb == na )
      uv[0] = lambda[1];
    else if (ntempb == nb )
      uv[0] = lambda[2];
    else {
      assert(ntempb == n0);
      uv[0] = lambda[0];
    }
    if ( ntempc == na )
      uv[1] = lambda[1];
    else if (ntempc == nb )
      uv[1] = lambda[2];
    else {
      assert(ntempc == n0);
      uv[1] = lambda[0];
    }
  }
  /* nb = lispoi[kel] -> lambda[1], na = lispoi[kel+1] -> lambda[2] */
  else {
    if ( ntempb == na )
      uv[0] = lambda[2];
    else if (ntempb == nb )
      uv[0] = lambda[1];
    else {
      assert(ntempb == n0);
      uv[0] = lambda[0];
    }
    if ( ntempc == na )
      uv[1] = lambda[2];
    else if (ntempc == nb )
      uv[1] = lambda[1];
    else {
      assert(ntempc == n0);
      uv[1] = lambda[0];
    }
  }
  if(!_MMG3D_bezierInt(&pb,uv,o,no,to)){
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* Test : make sure that geometric approximation has not been degraded too much */
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];

  ppt0->tag      = p0->tag;
  ppt0->ref      = p0->ref;


  nxp = mesh->xp + 1;
  if ( nxp > mesh->xpmax ) {
    _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                       "larger xpoint table",
                       _MMG5_SAFE_FREE(callist);
                       return(0));
    n = &(mesh->xpoint[p0->xp].n1[0]);
  }
  ppt0->xp = nxp;
  pxp = &mesh->xpoint[nxp];
  memcpy(pxp,&(mesh->xpoint[p0->xp]),sizeof(MMG5_xPoint));
  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  // parallel transport of metric at p0 to new point.
  if ( !_MMG5_paratmet(p0->c,n,m0,o,no,&met->m[0]) ) {
    _MMG5_SAFE_FREE(callist);
    return 0;
  }

  /* For each surfacic triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilists; l++) {
    k                   = lists[l] / 4;
    iface = lists[l] % 4;
    pt          = &mesh->tetra[k];
    _MMG5_tet2tri(mesh,k,iface,&tt);
    calold = MG_MIN(calold,_MMG5_caltri(mesh,met,&tt));
    for( i=0 ; i<3 ; i++ )
      if ( tt.v[i] == n0 )      break;
    assert(i<3);
    tt.v[i] = 0;
    caltmp = _MMG5_caltri(mesh,met,&tt);
    if ( caltmp < _MMG5_EPSD ) {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if (calnew < _MMG5_NULKAL) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if (improve && calnew < 1.02*calold) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if ( calnew < 0.3*calold ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
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
    callist[l]=_MMG5_orcal(mesh,met,0);
    if ( callist[l] < _MMG5_EPSD )  {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
    calnew = MG_MIN(calnew,callist[l]);
  }

  if ( calold < _MMG5_NULKAL && calnew <= calold ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if (calnew < _MMG5_NULKAL) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if (improve && calnew < calold) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
  else if ( calnew < 0.3*calold ) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* When all tests have been carried out, update coordinates, normals and metrics*/
  if ( octree )
    _MMG3D_moveOctree(mesh, octree, n0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  n[0] = no[0];
  n[1] = no[1];
  n[2] = no[2];

  memcpy(m0,&met->m[0],6*sizeof(double));

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual= callist[l];
  }
  _MMG5_SAFE_FREE(callist);
  return(1);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
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
int _MMG5_movbdyrefpt_ani(MMG5_pMesh mesh, MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv, int *lists, int ilists,
                          int improve){
  MMG5_pTetra           pt,pt0;
  MMG5_pPoint           p0,ppt0;
  MMG5_Tria             tt;
  MMG5_pxPoint          pxp;
  double                step,ll1old,ll2old,l1new,l2new;
  double                o[3],no[3],to[3];
  double                calold,calnew,caltmp,*callist;
  int                   l,iel,ip0,ipa,ipb,iptmpa,iptmpb,it1,it2,ip1,ip2,ip,nxp;
  int16_t               tag;
  unsigned char         i,i0,ie,iface,iface1,iface2,iea,ieb,ie1,ie2;

  step = 0.1;
  ip1 = ip2 = 0;
  pt    = &mesh->tetra[listv[0]/4];
  ip0 = pt->v[listv[0]%4];
  p0    = &mesh->point[ip0];

  assert ( MG_REF & p0->tag );

  /* Travel surfacic ball and recover the two ending points of ref curve :
     two senses must be used */
  iel = lists[0]/4;
  iface = lists[0]%4;
  pt = &mesh->tetra[iel];
  ipa = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[_MMG5_idir[iface][i]];
      else
        ipb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=1; l<ilists; l++) {
    iel = lists[l]/4;
    iface = lists[l]%4;
    pt = &mesh->tetra[iel];
    iea = ieb = 0;
    for (i=0; i<3; i++) {
      ie = _MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[_MMG5_iare[ie][0]] == ip0) || (pt->v[_MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[_MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[_MMG5_iare[iea][0]];
    else {
      assert(pt->v[_MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[_MMG5_iare[iea][1]];
    }
    if ( pt->v[_MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[_MMG5_iare[ieb][0]];
    else {
      assert(pt->v[_MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[_MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( MG_REF & tag ) {
        it1 = iel;
        ip1 = iptmpa;
        ie1 = iea;
        iface1 = iface;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[ieb];
      else  tag = 0;
      if ( MG_REF & tag ) {
        it1 = iel;
        ip1 = iptmpb;
        ie1 = ieb;
        iface1 = iface;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }

  /* Now travel surfacic list in the reverse sense so as to get the second ridge */
  iel = lists[0]/4;
  iface = lists[0]%4;
  pt = &mesh->tetra[iel];
  ipa = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[_MMG5_idir[iface][i]];
      else
        ipb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=ilists-1; l>0; l--) {
    iel         = lists[l] / 4;
    iface = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    iea         = ieb = 0;
    for (i=0; i<3; i++) {
      ie = _MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[_MMG5_iare[ie][0]] == ip0) || (pt->v[_MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[_MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[_MMG5_iare[iea][0]];
    else {
      assert(pt->v[_MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[_MMG5_iare[iea][1]];
    }
    if ( pt->v[_MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[_MMG5_iare[ieb][0]];
    else {
      assert(pt->v[_MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[_MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( MG_REF & tag ) {
        it2 = iel;
        ip2 = iptmpa;
        ie2 = iea;
        iface2 = iface;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      assert(pt->xt);
      tag = mesh->xtetra[pt->xt].tag[ieb];
      if ( MG_REF & tag ) {
        it2 = iel;
        ip2 = iptmpb;
        ie2 = ieb;
        iface2 = iface;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }
  if ( !(ip1 && ip2 && (ip1 != ip2)) )  return(0);

  /* At this point, we get the point extremities of the ref limit curve passing through ip0 :
     ip1, ip2, along with support tets it1,it2, the surface faces iface1,iface2, and the
     associated edges ie1,ie2.*/

  /* Changes needed for choice of time step : see manuscript notes */
  ll1old = _MMG5_lenSurfEdg(mesh,met,ip0,ip1,0);
  ll2old = _MMG5_lenSurfEdg(mesh,met,ip0,ip2,0);
  if ( ll1old < ll2old ) { //move towards p2
    iel = it2;
    ie  = ie2;
    iface = iface2;
    ip = ip2;
  }
  else {
    iel = it1;
    ie  = ie1;
    iface = iface1;
    ip = ip1;
  }

  /* Compute support of the associated edge, and features of the new position */
  if ( !(_MMG5_BezierRef(mesh,ip0,ip,step,o,no,to)) )  return(0);

  /* Test : make sure that geometric approximation has not been degraded too much */
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];
  ppt0->tag  = p0->tag;
  ppt0->ref  = p0->ref;


  nxp = mesh->xp + 1;
  if ( nxp > mesh->xpmax ) {
    _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                       "larger xpoint table",
                       return(0));
  }
  ppt0->xp = nxp;
  pxp = &mesh->xpoint[nxp];
  memcpy(pxp,&(mesh->xpoint[p0->xp]),sizeof(MMG5_xPoint));

  ppt0->n[0] = to[0];
  ppt0->n[1] = to[1];
  ppt0->n[2] = to[2];

  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  /* Interpolation of metric between ip0 and ip2 */
  if ( !_MMG5_paratmet(p0->c,mesh->xpoint[p0->xp].n1,&met->m[6*ip0],o,no,&met->m[0]) )
    return(0);

  /* Check whether proposed move is admissible under consideration of distances */
  l1new = _MMG5_lenSurfEdg(mesh,met,0,ip1,0);
  l2new = _MMG5_lenSurfEdg(mesh,met,0,ip2,0);
  if ( fabs(l2new -l1new) >= fabs(ll2old -ll1old) )
    return(0);

  /* For each surface triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilists ; l++ ){
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    _MMG5_tet2tri(mesh,iel,iface,&tt);
    calold = MG_MIN(calold,_MMG5_caltri(mesh,met,&tt));
    for( i=0 ; i<3 ; i++ )
      if ( tt.v[i] == ip0 )      break;
    assert(i<3);
    tt.v[i] = 0;
    caltmp = _MMG5_caltri(mesh,met,&tt);
    if ( caltmp < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )    return(0);
  else if ( calnew < calold )    return(0);
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */
  // Dynamic allocations for windows compatibility
  _MMG5_SAFE_MALLOC(callist, ilistv, double);

  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilistv ; l++ ){
    iel = listv[l] / 4;
    i0  = listv[l] % 4;
    pt  = &mesh->tetra[iel];
    pt0 = &mesh->tetra[0];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    calold = MG_MIN(calold, pt->qual);
    callist[l] = _MMG5_orcal(mesh,met,0);
    if (callist[l] < _MMG5_EPSD) {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
    calnew = MG_MIN(calnew,callist[l]);
  }
  if ((calold < _MMG5_NULKAL && calnew <= calold) ||
      (calnew < _MMG5_NULKAL) || (calnew <= 0.3*calold)) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  } else if (improve && calnew < calold) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* Update coordinates, normals, for new point */
  if ( octree )
    _MMG3D_moveOctree(mesh, octree, ip0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  pxp = &mesh->xpoint[p0->xp];
  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  p0->n[0] = to[0];
  p0->n[1] = to[1];
  p0->n[2] = to[2];

  memcpy(&met->m[6*ip0],&met->m[0],6*sizeof(double));

  for( l=0 ; l<ilistv ; l++ ){
    (&mesh->tetra[listv[l]/4])->qual = callist[l];
  }
  _MMG5_SAFE_FREE(callist);
  return(1);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
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
int _MMG5_movbdynompt_ani(MMG5_pMesh mesh,MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv, int *lists, int ilists,
                          int improve){
  MMG5_pTetra       pt,pt0;
  MMG5_pPoint       p0,ppt0;
  MMG5_pxPoint      pxp;
  MMG5_Tria         tt;
  double            step,ll1old,ll2old,l1new,l2new;
  double            calold,calnew,caltmp,*callist;
  double            o[3],no[3],to[3];
  int               ip0,ip1,ip2,ip,iel,ipa,ipb,l,iptmpa,iptmpb,it1,it2,nxp;
  int16_t           tag;
  char              iface,i,i0,iea,ieb,ie,ie1,ie2,iface1,iface2;

  step = 0.1;
  ip1 = ip2 = 0;
  pt = &mesh->tetra[listv[0]/4];
  ip0 = pt->v[listv[0]%4];
  p0 = &mesh->point[ip0];

  assert ( p0->tag & MG_NOM );

  /* Travel surfacic ball and recover the two ending points of non manifold curve :
     two senses must be used */
  iel = lists[0]/4;
  iface = lists[0]%4;
  pt = &mesh->tetra[iel];
  ipa = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[_MMG5_idir[iface][i]];
      else
        ipb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=1; l<ilists; l++) {
    iel = lists[l]/4;
    iface = lists[l]%4;
    pt = &mesh->tetra[iel];
    iea = ieb = 0;
    for (i=0; i<3; i++) {
      ie = _MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[_MMG5_iare[ie][0]] == ip0) || (pt->v[_MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[_MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[_MMG5_iare[iea][0]];
    else {
      assert(pt->v[_MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[_MMG5_iare[iea][1]];
    }
    if ( pt->v[_MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[_MMG5_iare[ieb][0]];
    else {
      assert(pt->v[_MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[_MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( MG_NOM & tag ) {
        it1 = iel;
        ip1 = iptmpa;
        ie1 = iea;
        iface1 = iface;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[ieb];
      else  tag = 0;
      if ( MG_NOM & tag ) {
        it1 = iel;
        ip1 = iptmpb;
        ie1 = ieb;
        iface1 = iface;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }

  /* Now travel surfacic list in the reverse sense so as to get the second non manifold point */
  iel = lists[0]/4;
  iface = lists[0]%4;
  pt = &mesh->tetra[iel];
  ipa = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[_MMG5_idir[iface][i]];
      else
        ipb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=ilists-1; l>0; l--) {
    iel         = lists[l] / 4;
    iface = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    iea         = ieb = 0;
    for (i=0; i<3; i++) {
      ie = _MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[_MMG5_iare[ie][0]] == ip0) || (pt->v[_MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[_MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[_MMG5_iare[iea][0]];
    else {
      assert(pt->v[_MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[_MMG5_iare[iea][1]];
    }
    if ( pt->v[_MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[_MMG5_iare[ieb][0]];
    else {
      assert(pt->v[_MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[_MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( MG_NOM & tag ) {
        it2 = iel;
        ip2 = iptmpa;
        ie2 = iea;
        iface2 = iface;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      assert(pt->xt);
      tag = mesh->xtetra[pt->xt].tag[ieb];
      if ( MG_NOM & tag ) {
        it2 = iel;
        ip2 = iptmpb;
        ie2 = ieb;
        iface2 = iface;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }
  if ( !(ip1 && ip2 && (ip1 != ip2)) )  return(0);

  /* At this point, we get the point extremities of the non manifold curve passing through ip0 :
     ip1, ip2, along with support tets it1,it2, the surface faces iface1,iface2, and the
     associated edges ie1,ie2.*/
  ll1old = _MMG5_lenSurfEdg(mesh,met,ip0,ip1,0);
  ll2old = _MMG5_lenSurfEdg(mesh,met,ip0,ip2,0);

  if ( ll1old < ll2old ) { //move towards p2
    iel = it2;
    ie  = ie2;
    iface = iface2;
    ip = ip2;
  }
  else {
    iel = it1;
    ie  = ie1;
    iface = iface1;
    ip = ip1;
  }

  /* Compute support of the associated edge, and features of the new position */
  if ( !(_MMG5_BezierNom(mesh,ip0,ip,step,o,no,to)) )  return(0);

  /* Test : make sure that geometric approximation has not been degraded too much */
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];
  ppt0->tag  = p0->tag;
  ppt0->ref  = p0->ref;

  nxp = mesh->xp + 1;
  if ( nxp > mesh->xpmax ) {
    _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                       "larger xpoint table",
                       return(0));
  }
  ppt0->xp = nxp;
  pxp = &mesh->xpoint[nxp];
  memcpy(pxp,&(mesh->xpoint[p0->xp]),sizeof(MMG5_xPoint));

  ppt0->n[0] = to[0];
  ppt0->n[1] = to[1];
  ppt0->n[2] = to[2];

  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  /* Interpolation of metric between ip0 and ip2 */
  if ( !_MMG5_paratmet(p0->c,mesh->xpoint[p0->xp].n1,&met->m[6*ip0],o,no,&met->m[0]) )
    return(0);

  /* Check whether proposed move is admissible under consideration of distances */
  l1new = _MMG5_lenSurfEdg(mesh,met,0,ip1,0);
  l2new = _MMG5_lenSurfEdg(mesh,met,0,ip2,0);
  if ( fabs(l2new -l1new) >= fabs(ll2old -ll1old) )
    return(0);

  /* For each surface triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilists ; l++ ){
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    _MMG5_tet2tri(mesh,iel,iface,&tt);
    caltmp = _MMG5_caltri(mesh,met,&tt);
    calold = MG_MIN(calold,caltmp);
    for( i=0 ; i<3 ; i++ )
      if ( tt.v[i] == ip0 )      break;
    assert(i<3);

    tt.v[i] = 0;
    caltmp = _MMG5_caltri(mesh,met,&tt);
    if ( caltmp < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )    return(0);
  else if ( calnew < calold )    return(0);
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */
  // Dynamic allocations for windows compatibility
  _MMG5_SAFE_MALLOC(callist, ilistv, double);

  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilistv ; l++ ){
    iel = listv[l] / 4;
    i0  = listv[l] % 4;
    pt  = &mesh->tetra[iel];
    pt0 = &mesh->tetra[0];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    calold = MG_MIN(calold, pt->qual);
    callist[l]= _MMG5_orcal(mesh,met,0);
    if (callist[l] < _MMG5_EPSD) {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
    calnew = MG_MIN(calnew,callist[l]);
  }
  if ((calold < _MMG5_NULKAL && calnew <= calold) ||
      (calnew < _MMG5_NULKAL) || (calnew <= 0.3*calold)) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  } else if (improve && calnew < calold) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* Update coordinates, normals, for new point */
  if ( octree )
    _MMG3D_moveOctree(mesh, octree, ip0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  pxp = &mesh->xpoint[p0->xp];
  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

  p0->n[0] = to[0];
  p0->n[1] = to[1];
  p0->n[2] = to[2];

  memcpy(&met->m[6*ip0],&met->m[0],6*sizeof(double));

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual = callist[l];
  }
  _MMG5_SAFE_FREE(callist);
  return(1);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param listv pointer toward the volumic ball of the point.
 * \param ilistv size of the volumic ball.
 * \param lists pointer toward the surfacic ball of the point.
 * \param ilists size of the surfacic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * \return 0 if fail, 1 if success.
 *
 * \remark we don't check if we break the hausdorff criterion.
 *
 * Move boundary ridge point, whose volumic and surfacic balls are passed.
 *
 */
int _MMG5_movbdyridpt_ani(MMG5_pMesh mesh, MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv,int *lists,int ilists,
                          int improve) {
  MMG5_pTetra          pt,pt0;
  MMG5_pPoint          p0,ppt0;
  MMG5_Tria            tt;
  MMG5_pxPoint         pxp;
  double               step,l1old,l2old,l1new,l2new;
  double               o[3],no1[3],no2[3],to[3];
  double               calold,calnew,caltmp,*callist;
  int                  l,iel,ip0,ipa,ipb,iptmpa,iptmpb,it1,it2,ip1,ip2,ip,nxp;
  int16_t              tag;
  unsigned char        i,i0,ie,iface,iface1,iface2,iea,ieb,ie1,ie2;

  step = 0.1;
  ip1 = ip2 = 0;
  pt    = &mesh->tetra[listv[0]/4];
  ip0 = pt->v[listv[0]%4];
  p0    = &mesh->point[ip0];

  assert ( MG_GEO & p0->tag );

  /* Travel surfacic ball an recover the two ending points of ridge : two senses
     must be used POSSIBLE OPTIMIZATION HERE : One travel only is needed */
  iel           = lists[0] / 4;
  iface         = lists[0] % 4;
  pt            = &mesh->tetra[iel];
  ipa           = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[_MMG5_idir[iface][i]];
      else
        ipb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=1; l<ilists; l++) {
    iel         = lists[l] / 4;
    iface = lists[l] % 4;
    pt  = &mesh->tetra[iel];
    iea = ieb = 0;
    for (i=0; i<3; i++) {
      ie = _MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[_MMG5_iare[ie][0]] == ip0) || (pt->v[_MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[_MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[_MMG5_iare[iea][0]];
    else {
      assert(pt->v[_MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[_MMG5_iare[iea][1]];
    }
    if ( pt->v[_MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[_MMG5_iare[ieb][0]];
    else {
      assert(pt->v[_MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[_MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( MG_GEO & tag ) {
        it1 = iel;
        ip1 = iptmpa;
        ie1 = iea;
        iface1 = iface;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[ieb];
      else  tag = 0;
      if ( MG_GEO & tag ) {
        it1 = iel;
        ip1 = iptmpb;
        ie1 = ieb;
        iface1 = iface;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }

  /* Now travel surfacic list in the reverse sense so as to get the second ridge */
  iel           = lists[0] / 4;
  iface = lists[0] % 4;
  pt    = &mesh->tetra[iel];
  ipa = ipb = 0;
  for (i=0; i<3; i++) {
    if ( pt->v[_MMG5_idir[iface][i]] != ip0 ) {
      if ( !ipa )
        ipa = pt->v[_MMG5_idir[iface][i]];
      else
        ipb = pt->v[_MMG5_idir[iface][i]];
    }
  }
  assert(ipa && ipb);

  for (l=ilists-1; l>0; l--) {
    iel         = lists[l]/4;
    iface = lists[l]%4;
    pt  = &mesh->tetra[iel];
    iea = ieb = 0;
    for (i=0; i<3; i++) {
      ie = _MMG5_iarf[iface][i]; //edge i on face iface
      if ( (pt->v[_MMG5_iare[ie][0]] == ip0) || (pt->v[_MMG5_iare[ie][1]] == ip0) ) {
        if ( !iea )
          iea = ie;
        else
          ieb = ie;
      }
    }
    if ( pt->v[_MMG5_iare[iea][0]] != ip0 )
      iptmpa = pt->v[_MMG5_iare[iea][0]];
    else {
      assert(pt->v[_MMG5_iare[iea][1]] != ip0);
      iptmpa = pt->v[_MMG5_iare[iea][1]];
    }
    if ( pt->v[_MMG5_iare[ieb][0]] != ip0 )
      iptmpb = pt->v[_MMG5_iare[ieb][0]];
    else {
      assert(pt->v[_MMG5_iare[ieb][1]] != ip0);
      iptmpb = pt->v[_MMG5_iare[ieb][1]];
    }
    if ( (iptmpa == ipa) || (iptmpa == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[iea];
      else  tag = 0;
      if ( MG_GEO & tag ) {
        it2 = iel;
        ip2 = iptmpa;
        ie2 = iea;
        iface2 = iface;
        break;
      }
    }
    if ( (iptmpb == ipa) || (iptmpb == ipb) ) {
      if ( pt->xt )  tag = mesh->xtetra[pt->xt].tag[ieb];
      else  tag = 0;
      if ( MG_GEO & tag ) {
        it2 = iel;
        ip2 = iptmpb;
        ie2 = ieb;
        iface2 = iface;
        break;
      }
    }
    ipa = iptmpa;
    ipb = iptmpb;
  }
  if ( !(ip1 && ip2 && (ip1 != ip2)) ) return(0);

  /* At this point, we get the point extremities of the ridge curve passing through ip0 :
     ip1, ip2, along with support tets it1,it2, the surface faces iface1,iface2, and the
     associated edges ie1,ie2.*/

  /* Changes needed for choice of time step : see manuscript notes */
  l1old = _MMG5_lenSurfEdg(mesh,met,ip0,ip1,1);
  l2old = _MMG5_lenSurfEdg(mesh,met,ip0,ip2,1);
  l1old = l1old*l1old;
  l2old = l2old*l2old;

  if ( l1old < l2old ) { //move towards p2
    iel = it2;
    ie  = ie2;
    iface = iface2;
    ip = ip2;
  }
  else {
    iel = it1;
    ie  = ie1;
    iface = iface1;
    ip = ip1;
  }

  /* Compute support of the associated edge, and features of the new position */
  if ( !(_MMG5_BezierRidge(mesh,ip0,ip,step,o,no1,no2,to)) )  return(0);

  /* Test : make sure that geometric approximation has not been degraded too much */
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];
  ppt0->tag      = p0->tag;
  ppt0->ref      = p0->ref;

  nxp = mesh->xp+1;
  if ( nxp > mesh->xpmax ) {
    _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                       "larger xpoint table",
                       return(0));
  }
  ppt0->xp = nxp;
  pxp = &mesh->xpoint[nxp];
  memcpy(pxp,&(mesh->xpoint[p0->xp]),sizeof(MMG5_xPoint));

  ppt0->n[0] = to[0];
  ppt0->n[1] = to[1];
  ppt0->n[2] = to[2];

  pxp->n1[0] = no1[0];
  pxp->n1[1] = no1[1];
  pxp->n1[2] = no1[2];

  pxp->n2[0] = no2[0];
  pxp->n2[1] = no2[1];
  pxp->n2[2] = no2[2];

  /* Interpolation of metric between ip0 and ip2 */
  if ( !_MMG5_intridmet(mesh,met,ip0,ip,step,no1,&met->m[0]) ) return 0;

  /* Check whether proposed move is admissible under consideration of distances */
  l1new = _MMG5_lenSurfEdg(mesh,met,0,ip1,1);
  l2new = _MMG5_lenSurfEdg(mesh,met,0,ip2,1);
  if ( fabs(l2new -l1new) >= fabs(l2old -l1old) )
    return(0);

  /* For each surfacic triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilists; l++) {
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    _MMG5_tet2tri(mesh,iel,iface,&tt);
    calold = MG_MIN(calold,_MMG5_caltri(mesh,met,&tt));
    for (i=0; i<3; i++) {
      if ( tt.v[i] == ip0 )      break;
    }
    assert(i<3);
    tt.v[i] = 0;
    caltmp = _MMG5_caltri(mesh,met,&tt);
    if ( caltmp < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )    return(0);
  else if ( calnew <= calold )  return(0);
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */
  // Dynamic allocations for windows compatibility
  _MMG5_SAFE_MALLOC(callist, ilistv, double);

  calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    iel = listv[l] / 4;
    i0  = listv[l] % 4;
    pt  = &mesh->tetra[iel];
    pt0 = &mesh->tetra[0];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[i0] = 0;
    calold = MG_MIN(calold, pt->qual);
    callist[l]=_MMG5_orcal(mesh,met,0);
    if (callist[l] < _MMG5_EPSD) {
      _MMG5_SAFE_FREE(callist);
      return(0);
    }
    calnew = MG_MIN(calnew,callist[l]);
  }
  if ((calold < _MMG5_NULKAL && calnew <= calold) ||
      (calnew < _MMG5_NULKAL) || (calnew <= 0.3*calold)) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  } else if (improve && calnew < calold) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* Update coordinates, normals, for new point */
  if ( octree )
    _MMG3D_moveOctree(mesh, octree, ip0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  pxp = &mesh->xpoint[p0->xp];
  pxp->n1[0] = no1[0];
  pxp->n1[1] = no1[1];
  pxp->n1[2] = no1[2];

  pxp->n2[0] = no2[0];
  pxp->n2[1] = no2[1];
  pxp->n2[2] = no2[2];

  p0->n[0] = to[0];
  p0->n[1] = to[1];
  p0->n[2] = to[2];

  memcpy(&met->m[6*ip0],&met->m[0],6*sizeof(double));

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual = callist[l];
  }
  _MMG5_SAFE_FREE(callist);
  return(1);
}
