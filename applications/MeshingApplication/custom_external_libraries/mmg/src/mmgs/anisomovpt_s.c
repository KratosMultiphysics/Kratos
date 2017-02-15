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

#include "mmgs.h"

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
int movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist) {
  MMG5_pTria     pt,pt0;
  MMG5_pPoint    p0,p1,ppt0;
  _MMG5_Bezier   pb;
  double         r[3][3],ux,uy,uz,*n,area,lispoi[3*_MMG5_LMAX+1],*m0;//,m[6],mo[6];
  double         gv[2],detloc,step,lambda[3],o[3],no[3],to[3],uv[2];
  double         calold,calnew,caltmp;
  int            k,iel,kel,nump,nbeg,nend;
  char           i0,i1,i2,ier;
  static int     warn=0;
  step = 0.1;

  /* Make sure ball of point is closed */
  iel = list[0] / 3;
  i0  = list[0] % 3;
  i1  = _MMG5_inxt2[i0];

  pt   = &mesh->tria[iel];
  nump = pt->v[i0];
  nbeg = pt->v[i1];
  p0   = &mesh->point[nump];
  m0   = &met->m[6*nump];
  assert( !p0->tag );

  iel = list[ilist-1] / 3;
  i0  = list[ilist-1] % 3;
  i2  = _MMG5_iprv2[i0];

  pt   = &mesh->tria[iel];
  nend = pt->v[i2];
  if ( nbeg != nend )  return(0);

  /** Step 1 : Rotation matrix that sends normal at p0 to e_z */
  n = &(p0->n[0]);
  if ( !_MMG5_rotmatrix(n,r) )  return(0);

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt  = &mesh->tria[iel];
    p1  = &mesh->point[pt->v[i1]];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] = r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] = r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] = r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
  }

  /* list goes modulo ilist */
  lispoi[3*ilist+1] = lispoi[1];
  lispoi[3*ilist+2] = lispoi[2];
  lispoi[3*ilist+3] = lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilist-1; k++) {
    area = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( area < 0.0 )  return(0);
  }
  area = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
  if ( area < 0.0 )  return(0);

  /** Step 2 : Compute gradient towards optimal position = centre of mass of the
     ball, projected to tangent plane */
  gv[0] = 0.0;
  gv[1] = 0.0;

  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    pt = &mesh->tria[iel];
    if ( !_MMG5_bezierCP(mesh,pt,&pb,1) )  return(0);

    /* Compute integral of sqrt(T^J(xi)  M(P(xi)) J(xi)) * P(xi) over the triangle */
    if ( !_MMG5_elementWeight(mesh,met,pt,p0,&pb,r,gv) ) {
      if ( !warn ) {
        ++warn;
        fprintf(stderr,"  ## Warning: unable to compute optimal position for at least"
               " 1 point.\n" );
      }
      return(0);
    }
  }

  /* At this point : gv = - gradient of V = direction to follow */
  /* Step 3 : locate new point in the ball, and compute its barycentric coordinates */
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
    if ( k == ilist )  return(0);
  }
  else {
    for (k=ilist-1; k>=0; k--) {
      detloc = lispoi[3*k+1]*gv[1] - lispoi[3*k+2]*gv[0]; //orientation with respect to p2
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == -1 )  return(0);
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  iel = list[kel] / 3;

  area = - gv[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) \
    + gv[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(area) < _MMG5_EPSD )  return(0);
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
  if ( area < _MMG5_EPSD )  return(0);
  area = 1.0 / area;
  lambda[1] = lispoi[3*(kel+1)+2]*gv[0] - lispoi[3*(kel+1)+1]*gv[1];
  lambda[2] = -lispoi[3*(kel)+2]*gv[0] + lispoi[3*(kel)+1]*gv[1];
  lambda[1]*= (area);
  lambda[2]*= (area);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /* Step 4 : come back to original problem, and compute patch in triangle iel */
  iel = list[kel] / 3;
  i0  = list[kel] % 3;
  i1  = _MMG5_inxt2[i0];
  i2  = _MMG5_inxt2[i1];
  pt  = &mesh->tria[iel];

  ier = _MMG5_bezierCP(mesh,pt,&pb,1);
  assert(ier);

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

  ier = _MMGS_bezierInt(&pb,uv,o,no,to);
  assert(ier);

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
  _MMG5_paratmet(p0->c,p0->n,m0,o,no,&met->m[0]);

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
    if ( caltmp < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,caltmp);

    if ( calold < NULKAL && calnew <= calold )  return(0);
    else if ( calnew < 0.3*calold )      return(0);
    /* if ( chkedg(mesh,0) )  return(0); */
  }

  /* Finally, update coordinates and normals of point, if new position is accepted :*/
  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  p0->n[0] = no[0];
  p0->n[1] = no[1];
  p0->n[2] = no[2];

  memcpy(m0,&met->m[0],6*sizeof(double));

  return(1);
}

/* Compute movement of a ref, or ridge point whose ball is passed */
int movridpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist) {
  MMG5_pTria    pt,pt0;
  MMG5_pPoint   p0,p1,p2,ppt0;
  MMG5_pxPoint    go;
  _MMG5_Bezier   b;
  double  *m0,*m00,step,l1old,l2old,ll1old,ll2old,uv[2],o[3],nn1[3],nn2[3],to[3],mo[6];
  double   lam0,lam1,lam2,*no1,*no2,*np1,*np2;
  double   psn11,psn12,ps2,l1new,l2new,dd1,dd2,ddt,calold,calnew;
  int      it1,it2,ip0,ip1,ip2,k,iel,ier;
  char     voy1,voy2,isrid,isrid1,isrid2,i0,i1,i2;
//#warning this step is different than the one used on iso or for int pts in aniso
  step  = 0.2;
  isrid = isrid1 = isrid2 = 0;
  it1   = it2 = 0;
  ip1   = ip2 = 0;
  voy1  = voy2 = 0;

  /* First step : make sure 2 ridge edges pass through point, and recoved them,
     along with neighbouring triangles */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_inxt2[i1];
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
        fprintf(stderr,"   *** function movridptaniso : 3 ridge edges landing on point %d\n",pt->v[i0]);
        return(0);
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
        fprintf(stderr,"   *** function movridptaniso : 3 ridge edges landing on point %d\n",pt->v[i0]);
        return(0);
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

  l1old = _MMG5_lenSurfEdg(mesh,met,ip0,ip1,1);
  l2old = _MMG5_lenSurfEdg(mesh,met,ip0,ip2,1);
  ll1old = l1old*l1old;
  ll2old = l2old*l2old;

  /* Third step : infer arc length of displacement, parameterized over edges */
  /* Move is made towards p2 */
  if ( l2old > l1old ) {
    isrid = isrid2;
    pt = &mesh->tria[it2];

    ier = _MMG5_bezierCP(mesh,pt,&b,1);
    assert(ier);

    /* fill table uv */
    if ( pt->v[0] == ip0 ) {
      if ( pt->v[1] == ip2 ) {
        uv[0] = step;
        uv[1] = 0.0;
      }
      else if ( pt->v[2] == ip2 ) {
        uv[0] = 0.0;
        uv[1] = step;
      }
    }
    else if ( pt->v[0] == ip2 ) {
      if ( pt->v[1] == ip0 ) {
        uv[0] = 1.0 - step;
        uv[1] = 0.0;
      }
      else if ( pt->v[2] == ip0 ) {
        uv[0] = 0.0;
        uv[1] = 1.0-step;
      }
    }
    else {
      if ( pt->v[1] == ip0 ) {
        uv[0] = 1.0 - step;
        uv[1] = step;
      }
      else if ( pt->v[2] == ip0 ) {
        uv[0] = step;
        uv[1] = 1.0-step;
      }
    }
    ier = _MMGS_bezierInt(&b,uv,o,nn1,to);
    assert(ier);
  }
  /* move towards p1 */
  else {
    isrid = isrid1;
    pt = &mesh->tria[it1];

    ier = _MMG5_bezierCP(mesh,pt,&b,1);
    assert(ier);

    /* fill table uv */
    if ( pt->v[0] == ip0 ) {
      if ( pt->v[1] == ip1 ) {
        uv[0] = step;
        uv[1] = 0.0;
      }
      else if ( pt->v[2] == ip1 ) {
        uv[0] = 0.0;
        uv[1] = step;
      }
    }
    else if ( pt->v[0] == ip1 ) {
      if ( pt->v[1] == ip0 ) {
        uv[0] = 1.0 - step;
        uv[1] = 0.0;
      }
      else if ( pt->v[2] == ip0 ) {
        uv[0] = 0.0;
        uv[1] = 1.0-step;
      }
    }
    else {
      if ( pt->v[1] == ip0 ) {
        uv[0] = 1.0-step;
        uv[1] = step;
      }
      else if ( pt->v[2] == ip0 ) {
        uv[0] = step;
        uv[1] = 1.0-step;
      }
    }
    ier = _MMGS_bezierInt(&b,uv,o,nn1,to);
    assert(ier);
  }

  /* check whether proposed move is admissible uner consideration of distances */
  /* Normal and tangent vector updates */
  lam0 = (1.0-step)*(1.0-step);
  lam1 = 2.0*step*(1.0-step);
  lam2 = step*step;

  /* Move is made towards p2 */
  if ( l2old > l1old ) {
    no1 = &mesh->xpoint[p0->xp].n1[0];
    no2 = &mesh->xpoint[p0->xp].n2[0];

    if ( MS_SIN(p2->tag) ) {
      np1 = &mesh->xpoint[p0->xp].n1[0];
      np2 = &mesh->xpoint[p0->xp].n2[0];
    }
    else {
      np1 = &mesh->xpoint[p2->xp].n1[0];
      np2 = &mesh->xpoint[p2->xp].n2[0];
    }
    psn11 = no1[0]*np1[0] + no1[1]*np1[1] + no1[2]*np1[2];
    psn12 = no1[0]*np2[0] + no1[1]*np2[1] + no1[2]*np2[2];

    /* no1 goes with np1, no2 with np2 */
    if ( fabs(1.0-fabs(psn11)) < fabs(1.0-fabs(psn12)) ){
      nn1[0] = no1[0]+np1[0];
      nn1[1] = no1[1]+np1[1];
      nn1[2] = no1[2]+np1[2];

      nn2[0] = no2[0]+np2[0];
      nn2[1] = no2[1]+np2[1];
      nn2[2] = no2[2]+np2[2];

      ps2 = (p2->c[0]-p0->c[0])*nn1[0]+(p2->c[1]-p0->c[1])*nn1[1]+(p2->c[2]-p0->c[2])*nn1[2];
      if ( ll2old < _MMG5_EPSD )  return(0);
      ps2 *= (2.0 / ll2old);
      nn1[0] -= ps2*(p2->c[0]-p0->c[0]);
      nn1[1] -= ps2*(p2->c[1]-p0->c[1]);
      nn1[2] -= ps2*(p2->c[2]-p0->c[2]);

      ps2 = (p2->c[0]-p0->c[0])*nn2[0]+(p2->c[1]-p0->c[1])*nn2[1]+(p2->c[2]-p0->c[2])*nn2[2];
      ps2 *= (2.0/ll2old);
      nn2[0] -=  ps2*(p2->c[0]-p0->c[0]);
      nn2[1] -=  ps2*(p2->c[1]-p0->c[1]);
      nn2[2] -=  ps2*(p2->c[2]-p0->c[2]);

      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];
      if ( (dd1 < _MMG5_EPSD2) || (dd2<_MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      /* At this point, the control points associated to the interpolation formula for normals
         have been computed .*/
      nn1[0] = lam0*no1[0] + lam1*nn1[0] + lam2*np1[0];
      nn1[1] = lam0*no1[1] + lam1*nn1[1] + lam2*np1[1];
      nn1[2] = lam0*no1[2] + lam1*nn1[2] + lam2*np1[2];

      nn2[0] = lam0*no2[0] + lam1*nn2[0] + lam2*np2[0];
      nn2[1] = lam0*no2[1] + lam1*nn2[1] + lam2*np2[1];
      nn2[2] = lam0*no2[2] + lam1*nn2[2] + lam2*np2[2];

      to[0] = nn1[1]*nn2[2]-nn1[2]*nn2[1];
      to[1] = nn1[2]*nn2[0]-nn1[0]*nn2[2];
      to[2] = nn1[0]*nn2[1]-nn1[1]*nn2[0];

      ddt = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];

      if ( (dd1 < _MMG5_EPSD2) || (dd2<_MMG5_EPSD2) || (ddt < _MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      ddt = 1.0 / sqrt(ddt);
      to[0] = ddt*to[0];
      to[1] = ddt*to[1];
      to[2] = ddt*to[2];
    }

    /* no1 goes with np2 and no2 with np1 */
    else {
      nn1[0] = no1[0]+np2[0];
      nn1[1] = no1[1]+np2[1];
      nn1[2] = no1[2]+np2[2];

      nn2[0] = no2[0]+np1[0];
      nn2[1] = no2[1]+np1[1];
      nn2[2] = no2[2]+np1[2];

      ps2 = (p2->c[0]-p0->c[0])*nn1[0]+(p2->c[1]-p0->c[1])*nn1[1]+(p2->c[2]-p0->c[2])*nn1[2];
      if ( ll2old < _MMG5_EPSD )  return(0);
      ps2 *= (2.0 / ll2old);
      nn1[0] -= ps2*(p2->c[0]-p0->c[0]);
      nn1[1] -= ps2*(p2->c[1]-p0->c[1]);
      nn1[2] -= ps2*(p2->c[2]-p0->c[2]);

      ps2 = (p2->c[0]-p0->c[0])*nn2[0]+(p2->c[1]-p0->c[1])*nn2[1]+(p2->c[2]-p0->c[2])*nn2[2];
      ps2 *= (2.0 / ll2old);
      nn2[0] -=  ps2*(p2->c[0]-p0->c[0]);
      nn2[1] -=  ps2*(p2->c[1]-p0->c[1]);
      nn2[2] -=  ps2*(p2->c[2]-p0->c[2]);

      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];

      if (( dd1 < _MMG5_EPSD2 ) || (dd2<_MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      /* At this point, the control points associated to the interpolation formula for normals
         have been computed .*/
      nn1[0] = lam0*no1[0] + lam1*nn1[0] + lam2*np2[0];
      nn1[1] = lam0*no1[1] + lam1*nn1[1] + lam2*np2[1];
      nn1[2] = lam0*no1[2] + lam1*nn1[2] + lam2*np2[2];

      nn2[0] = lam0*no2[0] + lam1*nn2[0] + lam2*np1[0];
      nn2[1] = lam0*no2[1] + lam1*nn2[1] + lam2*np1[1];
      nn2[2] = lam0*no2[2] + lam1*nn2[2] + lam2*np1[2];

      to[0] = nn1[1]*nn2[2]-nn1[2]*nn2[1];
      to[1] = nn1[2]*nn2[0]-nn1[0]*nn2[2];
      to[2] = nn1[0]*nn2[1]-nn1[1]*nn2[0];

      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];
      ddt = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];

      if ( (dd1 < _MMG5_EPSD2) || (dd2<_MMG5_EPSD2) || (ddt < _MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      ddt = 1.0 / sqrt(ddt);
      to[0] = ddt*to[0];
      to[1] = ddt*to[1];
      to[2] = ddt*to[2];
    }

    /* Interpolation of metric between ip0 and ip2 */
    if ( isrid ) {
      if(!_MMG5_intridmet(mesh,met,mesh->tria[it2].v[_MMG5_inxt2[voy2]],
                          mesh->tria[it2].v[_MMG5_iprv2[voy2]],(1.0-step),nn1,mo))
        return 0;
    }
    else {
      if ( !_MMG5_paratmet(p0->c,p0->n,m0,o,nn1,mo) )  return(0);
    }
  }

  /* Move along p1 */
  else {
    no1 = &mesh->xpoint[p0->xp].n1[0];
    no2 = &mesh->xpoint[p0->xp].n2[0];
    if ( MS_SIN(p1->tag) ) {
      np1 = &mesh->xpoint[p0->xp].n1[0];
      np2 = &mesh->xpoint[p0->xp].n2[0];
    }
    else {
      np1 = &mesh->xpoint[p1->xp].n1[0];
      np2 = &mesh->xpoint[p1->xp].n2[0];
    }
    psn11 = no1[0]*np1[0] + no1[1]*np1[1] + no1[2]*np1[2];
    psn12 = no1[0]*np2[0] + no1[1]*np2[1] + no1[2]*np2[2];

    /* no1 goes with np1, no2 with np2 */
    if ( fabs(1.0-fabs(psn11)) < fabs(1.0-fabs(psn12)) ) {
      nn1[0] = no1[0]+np1[0];
      nn1[1] = no1[1]+np1[1];
      nn1[2] = no1[2]+np1[2];

      nn2[0] = no2[0]+np2[0];
      nn2[1] = no2[1]+np2[1];
      nn2[2] = no2[2]+np2[2];

      ps2 = (p1->c[0]-p0->c[0])*nn1[0]+(p1->c[1]-p0->c[1])*nn1[1]+(p1->c[2]-p0->c[2])*nn1[2];
      if ( ll1old < _MMG5_EPSD )  return(0);
      ps2 *= (2.0 / ll1old);
      nn1[0] -= ps2*(p1->c[0]-p0->c[0]);
      nn1[1] -= ps2*(p1->c[1]-p0->c[1]);
      nn1[2] -= ps2*(p1->c[2]-p0->c[2]);

      ps2 = (p1->c[0]-p0->c[0])*nn2[0]+(p1->c[1]-p0->c[1])*nn2[1]+(p1->c[2]-p0->c[2])*nn2[2];
      ps2 *= (2.0 / ll1old);
      nn2[0] -=  ps2*(p1->c[0]-p0->c[0]);
      nn2[1] -=  ps2*(p1->c[1]-p0->c[1]);
      nn2[2] -=  ps2*(p1->c[2]-p0->c[2]);

      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];

      if ( (dd1 < _MMG5_EPSD2 ) || (dd2<_MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      /* Interpolation formula from the control point */
      nn1[0] = lam0*no1[0] + lam1*nn1[0] + lam2*np1[0];
      nn1[1] = lam0*no1[1] + lam1*nn1[1] + lam2*np1[1];
      nn1[2] = lam0*no1[2] + lam1*nn1[2] + lam2*np1[2];

      nn2[0] = lam0*no2[0] + lam1*nn2[0] + lam2*np2[0];
      nn2[1] = lam0*no2[1] + lam1*nn2[1] + lam2*np2[1];
      nn2[2] = lam0*no2[2] + lam1*nn2[2] + lam2*np2[2];

      to[0] = nn1[1]*nn2[2]-nn1[2]*nn2[1];
      to[1] = nn1[2]*nn2[0]-nn1[0]*nn2[2];
      to[2] = nn1[0]*nn2[1]-nn1[1]*nn2[0];

      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];
      ddt = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];

      if ( (dd1 < _MMG5_EPSD2) || (dd2<_MMG5_EPSD2) || (ddt < _MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      ddt = 1.0 / sqrt(ddt);
      to[0] = ddt*to[0];
      to[1] = ddt*to[1];
      to[2] = ddt*to[2];
    }
    /* no1 goes with np2 and no2 with np1 */
    else {
      nn1[0] = no1[0]+np2[0];
      nn1[1] = no1[1]+np2[1];
      nn1[2] = no1[2]+np2[2];

      nn2[0] = no2[0]+np1[0];
      nn2[1] = no2[1]+np1[1];
      nn2[2] = no2[2]+np1[2];

      ps2 = (p1->c[0]-p0->c[0])*nn1[0]+(p1->c[1]-p0->c[1])*nn1[1]+(p1->c[2]-p0->c[2])*nn1[2];
      if ( ll1old < _MMG5_EPSD )  return(0);
      ps2 *= (2.0 / ll1old);
      nn1[0] -= ps2*(p1->c[0]-p0->c[0]);
      nn1[1] -= ps2*(p1->c[1]-p0->c[1]);
      nn1[2] -= ps2*(p1->c[2]-p0->c[2]);

      ps2 = (p1->c[0]-p0->c[0])*nn2[0]+(p1->c[1]-p0->c[1])*nn2[1]+(p1->c[2]-p0->c[2])*nn2[2];
      ps2 *= (2.0 / ll1old);
      nn2[0] -=  ps2*(p1->c[0]-p0->c[0]);
      nn2[1] -=  ps2*(p1->c[1]-p0->c[1]);
      nn2[2] -=  ps2*(p1->c[2]-p0->c[2]);

      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];

      if ( (dd1 < _MMG5_EPSD2) || (dd2<_MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      /* Interpolation formula from the control point */
      nn1[0] = lam0*no1[0] + lam1*nn1[0] + lam2*np2[0];
      nn1[1] = lam0*no1[1] + lam1*nn1[1] + lam2*np2[1];
      nn1[2] = lam0*no1[2] + lam1*nn1[2] + lam2*np2[2];

      nn2[0] = lam0*no2[0] + lam1*nn2[0] + lam2*np1[0];
      nn2[1] = lam0*no2[1] + lam1*nn2[1] + lam2*np1[1];
      nn2[2] = lam0*no2[2] + lam1*nn2[2] + lam2*np1[2];

      to[0] = nn1[1]*nn2[2]-nn1[2]*nn2[1];
      to[1] = nn1[2]*nn2[0]-nn1[0]*nn2[2];
      to[2] = nn1[0]*nn2[1]-nn1[1]*nn2[0];

      dd1 = nn1[0]*nn1[0] + nn1[1]*nn1[1] + nn1[2]*nn1[2];
      dd2 = nn2[0]*nn2[0] + nn2[1]*nn2[1] + nn2[2]*nn2[2];
      ddt = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];

      if ( (dd1 < _MMG5_EPSD2) || (dd2<_MMG5_EPSD2) || (ddt < _MMG5_EPSD2) )  return(0);
      dd1 = 1.0 / sqrt(dd1);
      nn1[0] = dd1*nn1[0];
      nn1[1] = dd1*nn1[1];
      nn1[2] = dd1*nn1[2];

      dd2 = 1.0 / sqrt(dd2);
      nn2[0] = dd2*nn2[0];
      nn2[1] = dd2*nn2[1];
      nn2[2] = dd2*nn2[2];

      ddt = 1.0 / sqrt(ddt);
      to[0] = ddt*to[0];
      to[1] = ddt*to[1];
      to[2] = ddt*to[2];
    }

    /* Interpolation of metric between ip0 and ip1 */
    if ( isrid ) {
      _MMG5_intridmet(mesh,met,mesh->tria[it1].v[_MMG5_inxt2[voy1]],
                      mesh->tria[it1].v[_MMG5_iprv2[voy1]],(1.0-step),nn1,mo);
    }
    else {
      if ( !_MMG5_paratmet(p0->c,p0->n,m0,o,nn1,mo) )  return(0);
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
  l1new = _MMG5_lenSurfEdg(mesh,met,0,ip1,1);
  l2new = _MMG5_lenSurfEdg(mesh,met,0,ip2,1);
  if ( fabs(l2new -l1new) >= fabs(l2old -l1old) ) {
    ppt0->tag = 0;
    return(0);
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
      return(0);
    }
    else if ( calnew < 0.3*calold ) {
      ppt0->tag = 0;
      return(0);
    }
    /* if ( chkedg(mesh,0) )  return(0);*/
  }

  /* Coordinates, normals, tangents update */
  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  no1[0] = nn1[0];
  no1[1] = nn1[1];
  no1[2] = nn1[2];

  if ( isrid ) {
    no2[0] = nn2[0];
    no2[1] = nn2[1];
    no2[2] = nn2[2];
  }

  p0->n[0] = to[0];
  p0->n[1] = to[1];
  p0->n[2] = to[2];
  memcpy(m0,mo,6*sizeof(double));

  return(1);
}
