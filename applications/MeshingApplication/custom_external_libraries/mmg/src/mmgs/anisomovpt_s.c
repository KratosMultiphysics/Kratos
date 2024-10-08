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
 * \file mmgs/anisomovpt_s.c
 * \brief Functions to move a point in the mesh (with anisotropic metric).
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmgs_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param list ball of point.
 * \param ilist size of the point ball.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute movement of an internal point whose ball is passed.
 *
 */
int movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *list,int ilist) {
  MMG5_pTria     pt,pt0;
  MMG5_pPoint    p0,ppt0;
  MMG5_Bezier    pb;
  double         r[3][3],lispoi[3*MMG5_TRIA_LMAX+1],*m0;//,m[6],mo[6];
  double         gv[2],area,detloc,step,lambda[3],o[3],no[3],to[3],uv[2];
  double         calold,calnew,caltmp;
  MMG5_int       k,iel,kel,nump,nbeg,nend;
  int8_t         i0,i1,i2,ier;
  static int     warn=0;
  static int8_t  mmgErr0=0,mmgErr1=0;

  step = 0.1;

  /* Make sure ball of point is closed */
  iel = list[0] / 3;
  i0  = list[0] % 3;
  i1  = MMG5_inxt2[i0];

  pt   = &mesh->tria[iel];
  nump = pt->v[i0];
  nbeg = pt->v[i1];
  p0   = &mesh->point[nump];
  m0   = &met->m[6*nump];
  assert( !p0->tag );

  iel = list[ilist-1] / 3;
  i0  = list[ilist-1] % 3;
  i2  = MMG5_iprv2[i0];

  pt   = &mesh->tria[iel];
  nend = pt->v[i2];
  if ( nbeg != nend )  return 0;

  /* Rotation of the ball of p0 */
  if ( !MMGS_surfballRotation(mesh,p0,list,ilist,r,lispoi,p0->n)  ) {
    return 0;
  }

  /** Step 2 : Compute gradient towards optimal position = centre of mass of the
     ball, projected to tangent plane */
  gv[0] = 0.0;
  gv[1] = 0.0;

  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    pt = &mesh->tria[iel];
    if ( !MMG5_bezierCP(mesh,pt,&pb,1) )  return 0;

    /* Compute integral of sqrt(T^J(xi)  M(P(xi)) J(xi)) * P(xi) over the triangle */
    if ( !MMG5_elementWeight(mesh,met,pt,p0,&pb,r,gv) ) {
      if ( !warn ) {
        ++warn;
        fprintf(stderr,"\n  ## Warning: %s: unable to compute optimal position for at least"
                " 1 point.\n",__func__ );
      }
      return 0;
    }
  }

  /* At this point : gv = - gradient of V = direction to follow */
  /** Step 3 : locate new point in the ball, and compute its barycentric coordinates */
  area = lispoi[1]*gv[1] - lispoi[2]*gv[0];
  kel = 0;
  if ( area >= 0.0 ) {
    for (k=0; k<ilist; k++) {
      detloc = gv[0]*lispoi[3*(k+1)+2] - gv[1]*lispoi[3*(k+1)+1]; /*orientation with respect to p2 */
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == ilist )  return 0;
  }
  else {
    for (k=ilist-1; k>=0; k--) {
      detloc = lispoi[3*k+1]*gv[1] - lispoi[3*k+2]*gv[0]; //orientation with respect to p2
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == -1 )  return 0;
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  area = - gv[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) \
    + gv[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(area) < MMG5_EPSD2 )  return 0;
  area = 1.0 / area;
  step *= area;

  area = lispoi[3*(kel)+1]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]) \
    - lispoi[3*(kel)+2 ]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]);
  step *= area;
  step = fabs(step);
  gv[0] *= step;
  gv[1] *= step;

  /* Computation of the barycentric coordinates of the new point in the corresponding triangle. */
  area = lispoi[3*kel+1]*lispoi[3*(kel+1)+2] - lispoi[3*kel+2]*lispoi[3*(kel+1)+1];
  if ( area < MMG5_EPSD2 )  return 0;
  area = 1.0 / area;
  lambda[1] = lispoi[3*(kel+1)+2]*gv[0] - lispoi[3*(kel+1)+1]*gv[1];
  lambda[2] = -lispoi[3*(kel)+2]*gv[0] + lispoi[3*(kel)+1]*gv[1];
  lambda[1]*= (area);
  lambda[2]*= (area);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /** Step 4 : come back to original problem, and compute patch in triangle iel */
  iel = list[kel] / 3;
  i0  = list[kel] % 3;
  pt  = &mesh->tria[iel];

  ier = MMG5_bezierCP(mesh,pt,&pb,1);
  if ( !ier ) {
    if( !mmgErr0 ) {
      mmgErr0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: function MMG5_bezierCP return 0.\n",
              __func__);
    }
    return 0;
  }

  /* Now, for Bezier interpolation, one should identify which of i,i1,i2 is 0,1,2
     recall uv[0] = barycentric coord associated to pt->v[1], uv[1] associated to pt->v[2] */
  if ( i0 == 0 ) {
    uv[0] = lambda[1];
    uv[1] = lambda[2];
  }
  else if ( i0 == 1 ) {
    uv[0] = lambda[0];
    uv[1] = lambda[1];
  }
  else{
    uv[0] = lambda[2];
    uv[1] = lambda[0];
  }

  ier = MMGS_bezierInt(&pb,uv,o,no,to);
  if ( !ier ) {
    if( !mmgErr1 ) {
      mmgErr1 = 1;
      fprintf(stderr,"  ## Warning: %s: function MMGS_bezierInt return 0.\n",
              __func__);
    }
    return 0;
  }

  /* Second test : check whether geometric approximation has not been too much degraded */
  ppt0 = &mesh->point[0];
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];

  ppt0->n[0] = no[0];
  ppt0->n[1] = no[1];
  ppt0->n[2] = no[2];
  ppt0->tag  = 0;

  // parallel transport of metric at p0 to new point
  MMG5_paratmet(p0->c,p0->n,m0,o,no,&met->m[0]);

  calold = calnew = DBL_MAX;
  for (k= 0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    pt  = &mesh->tria[iel];
    pt0 = &mesh->tria[0];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i0] = 0;

    caltmp = caleltsig_ani(mesh,met,iel);
    calold = MG_MIN(calold,caltmp);
    caltmp = caleltsig_ani(mesh,met,0);
    if ( caltmp < MMG5_EPSD2 )  {
      /* We don't check the input triangle qualities, thus we may have a very
       * bad triangle in our mesh */
      return 0;
    }
    calnew = MG_MIN(calnew,caltmp);

    if ( calold < MMG5_EPSOK && calnew <= calold )  return 0;
    else if (calnew < MMG5_EPSOK)  return 0;
    else if ( calnew < 0.3*calold )      return 0;
    /* if ( chkedg(mesh,0) )  return 0; */
  }

  /* Finally, update coordinates and normals of point, if new position is accepted :*/
  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  p0->n[0] = no[0];
  p0->n[1] = no[1];
  p0->n[2] = no[2];

  memcpy(m0,&met->m[0],6*sizeof(double));

  return 1;
}

/* Compute movement of a ref, or ridge point whose ball is passed */
int movridpt_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *list,int ilist) {
  MMG5_pTria    pt,pt0;
  MMG5_pPoint   p0,p1,p2,ppt0;
  MMG5_pxPoint  go;
  double        *m0,*m00,step,l1old,l2old,ll1old,ll2old;
  double        lam0,lam1,lam2,o[3],nn1[3],nn2[3],to[3],mo[6];
  double        l1new,l2new,calold,calnew;
  MMG5_int      it,it1,it2,ip,ip0,ip1,ip2,k,iel;
  int8_t        voy1,voy2,isrid,isrid1,isrid2,i0,i1,i2;
  static int8_t mmgWarn0 = 0;

  step  = 0.2;
  isrid1 = isrid2 = 0;
  it1   = it2 = 0;
  ip1   = ip2 = 0;
  voy1  = voy2 = 0;

  /* First step : make sure 2 ridge edges pass through point, and recoved them,
     along with neighbouring triangles */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = MMG5_inxt2[i0];
    i2  = MMG5_inxt2[i1];
    pt  = &mesh->tria[iel];

    if ( MG_EDG(pt->tag[i1]) ) {
      if ( !it1 ) {
        it1  = iel;
        ip1  = pt->v[i2]; // edge(i1) = (p0p2)
        voy1 = i1;
        if ( pt->tag[i1] & MG_GEO )  isrid1 = 1;
      }
      else if ( it1 && !it2 ) {
        if ( ip1 != pt->v[i2] ) {
          it2  = iel;
          ip2  = pt->v[i2]; // edge(i1) = (p0p2)
          voy2 = i1;
          if ( pt->tag[i1] & MG_GEO )  isrid2 = 1;
        }
      }
      else if ( it1 && it2 && (pt->v[i2] != ip1) && (pt->v[i2] != ip2) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the"
                  " intersection of 3 ridge edges\n",__func__);
        }
        return 0;
      }
    }

    if ( MG_EDG(pt->tag[i2]) ) {
      if ( !it1 ) {
        it1  = iel;
        ip1  = pt->v[i1]; // edge(i2) = (p0p1)
        voy1 = i2;
        if ( pt->tag[i2] & MG_GEO )  isrid1 = 1;
      }
      else if ( it1 && !it2 ) {
        if ( ip1 != pt->v[i1] ) {
          it2  = iel;
          ip2  = pt->v[i1]; // edge(i1) = (p0p2)
          voy2 = i2;
          if ( pt->tag[i2] & MG_GEO )  isrid2 = 1;
        }
      }
      else if ( it1 && it2 && (pt->v[i1] != ip1) && (pt->v[i1] != ip2) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the"
                  " intersection of 3 ridge edges\n",__func__);
        }
        return 0;
      }
    }
  }

  /* Second step : compute lengths of edges (p0p1),(p0p2) */
  iel = list[0] / 3;
  i0  = list[0] % 3;
  pt  = &mesh->tria[iel];
  ip0 = pt->v[i0];
  p0  = &mesh->point[ip0];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];
  m0  = &met->m[6*ip0];

  l1old = MMG5_lenSurfEdg(mesh,met,ip0,ip1,1);
  l2old = MMG5_lenSurfEdg(mesh,met,ip0,ip2,1);

  if ( (!l1old) || (!l2old) ) return 0;

  if ( l1old < l2old ) { //move towards p2
    ip    = ip2;
    isrid = isrid2;
    it    = it2;
  }
  else {
    ip    = ip1;
    isrid = isrid1;
    it    = it1;
  }

  /* Third step : infer arc length of displacement, parameterized over edges */
  if ( !MMGS_paramDisp ( mesh,it,isrid,ip0,ip,step,o) ) {
    return 0;
  }

  /* check whether proposed move is admissible uner consideration of distances */
  /* Normal and tangent vector updates */
  lam0 = (1.0-step)*(1.0-step);
  lam1 = 2.0*step*(1.0-step);
  lam2 = step*step;

  /* Move is made towards p2 */
  ll1old = l1old*l1old;
  ll2old = l2old*l2old;

  if ( l2old > l1old ) {
    if ( !MMGS_moveTowardPoint(mesh,p0,p2,ll2old,lam0,lam1,lam2,nn1,nn2,to) ) {
      return 0;
    }

    /* Interpolation of metric between ip0 and ip2 */
    if ( isrid ) {
      if(!MMG5_intridmet(mesh,met,mesh->tria[it2].v[MMG5_inxt2[voy2]],
                          mesh->tria[it2].v[MMG5_iprv2[voy2]],(1.0-step),nn1,mo))
        return 0;
    }
    else {
      if ( !MMG5_paratmet(p0->c,p0->n,m0,o,nn1,mo) )  return 0;
    }
  }

  /* Move along p1 */
  else {
    if ( !MMGS_moveTowardPoint(mesh,p0,p1,ll1old,lam0,lam1,lam2,nn1,nn2,to) ) {
      return 0;
    }

    /* Interpolation of metric between ip0 and ip1 */
    if ( isrid ) {
      MMG5_intridmet(mesh,met,mesh->tria[it1].v[MMG5_inxt2[voy1]],
                      mesh->tria[it1].v[MMG5_iprv2[voy1]],(1.0-step),nn1,mo);
    }
    else {
      if ( !MMG5_paratmet(p0->c,p0->n,m0,o,nn1,mo) )  return 0;
    }
  }

  /* Check proposed motion */
  ppt0 = &mesh->point[0];
  ppt0->xp = 0;
  ppt0->tag = p0->tag;

  go = &mesh->xpoint[0];
  go->n1[0] = nn1[0];
  go->n1[1] = nn1[1];
  go->n1[2] = nn1[2];

  if ( isrid ) {
    go->n2[0] = nn2[0];
    go->n2[1] = nn2[1];
    go->n2[2] = nn2[2];
  }

  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];

  ppt0->n[0] = to[0];
  ppt0->n[1] = to[1];
  ppt0->n[2] = to[2];

  m00 = &met->m[0];
  memcpy(m00,mo,6*sizeof(double));

  /* Check whether proposed move is admissible under consideration of distances */
  l1new = MMG5_lenSurfEdg(mesh,met,0,ip1,1);
  l2new = MMG5_lenSurfEdg(mesh,met,0,ip2,1);

  if ( (!l1new) || (!l2new) ) return 0;

  if ( fabs(l2new -l1new) >= fabs(l2old -l1old) ) {
    ppt0->tag = 0;
    return 0;
  }

  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    pt  = &mesh->tria[iel];
    pt0 = &mesh->tria[0];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i0] = 0;

    calold = caleltsig_ani(mesh,met,iel);
    calnew = caleltsig_ani(mesh,met,0);
    if ( (calnew < 0.001) && (calnew<calold) ) {
      ppt0->tag = 0;
      return 0;
    }
    else if ( calnew < 0.3*calold ) {
      ppt0->tag = 0;
      return 0;
    }
    /* if ( chkedg(mesh,0) )  return 0;*/
  }

  /* Coordinates, normals, tangents update */
  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  mesh->xpoint[p0->xp].n1[0] = nn1[0];
  mesh->xpoint[p0->xp].n1[1] = nn1[1];
  mesh->xpoint[p0->xp].n1[2] = nn1[2];

  if ( isrid ) {
    mesh->xpoint[p0->xp].n2[0] = nn2[0];
    mesh->xpoint[p0->xp].n2[1] = nn2[1];
    mesh->xpoint[p0->xp].n2[2] = nn2[2];
  }

  p0->n[0] = to[0];
  p0->n[1] = to[1];
  p0->n[2] = to[2];
  memcpy(m0,mo,6*sizeof(double));

  return 1;
}
