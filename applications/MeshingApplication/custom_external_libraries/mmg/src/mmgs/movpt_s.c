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
 * \file mmgs/movpt_s.c
 * \brief Functions to move a point in the mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"
#include <math.h>


/* compute movement of an internal point whose ball is passed */
int movintpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist) {
  MMG5_pPoint   p0,p1,ppt0;
  MMG5_pTria    pt,pt0;
  _MMG5_Bezier   b;
  double   aa,bb,ab,ll,l,mlon,devmean,GV[3],gv[2],cosalpha,sinalpha,r[3][3],*n,lispoi[3*_MMG5_LMAX+1];
  double   ux,uy,uz,det2d,detloc,step,lambda[3],uv[2],o[3],no[3],to[3],Vold,Vnew,calold,calnew,caltmp;
  int      ier,iel,ipp,k,kel,npt,ibeg,iend;
  char     i0,i1,i2;

  step = 0.1;
  Vold = 0.0;
  Vnew = 0.0;

  k  = list[0] / 3;
  i0 = list[0] % 3;
  i1 = _MMG5_inxt2[i0];
  pt = &mesh->tria[k];
  ibeg = pt->v[i1];
  ipp  = pt->v[i0];
  p0   = &mesh->point[ipp]; /* point to move */

  k  = list[ilist-1] / 3;
  i0 = list[ilist-1] % 3;
  i1 = _MMG5_inxt2[i0];
  i2 = _MMG5_inxt2[i1];
  pt = &mesh->tria[k];
  iend = pt->v[i2];

  /* check for open ball */
  if ( iend != ibeg )  return(0);

  npt = ilist; // number of POINTS in the ball = number of triangles. Each point is counted as the
  // i1 of its associated triangle

  /* Step 1 : computation of the gradient of variance as a function from R^3 to R
     compute mlon := 1/n \sum_{i=1,...,n}{||p-p_i||^2} */
  mlon = 0.0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt  = &mesh->tria[iel];
    p1  = &mesh->point[pt->v[i1]];
    mlon += (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1]) \
      + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
  }
  mlon = mlon / npt;

  /* compute GV := 4/n * sum_{i=1,...,n}{(||p-p_i||^2 -m)(p-pi)} and
     Vold = sum_{i=1,...,n}{(||p-p_i||^2 -m)^2} */
  GV[0] = GV[1] = GV[2] = 0.0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt  = &mesh->tria[iel];
    p1  = &mesh->point[pt->v[i1]];

    devmean = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1]) \
      + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]) - mlon;
    GV[0] += devmean*(p1->c[0]-p0->c[0]);
    GV[1] += devmean*(p1->c[1]-p0->c[1]);
    GV[2] += devmean*(p1->c[2]-p0->c[2]);
    Vold  += devmean*devmean;
  }

  GV[0] *= (4.0 / npt);
  GV[1] *= (4.0 / npt);
  GV[2] *= (4.0 / npt);
  Vold  *= (1.0 / npt);

  /* Step 2 : computation of the rotation matrix T_p0 S -> [z = 0] */
  n  = &p0->n[0]; //once again, that depends on what kind of point is considered.
  aa = n[0]*n[0];
  bb = n[1]*n[1];
  ab = n[0]*n[1];
  ll = aa+bb;
  cosalpha = n[2];
  sinalpha = sqrt(1.0- MG_MIN(1.0,cosalpha*cosalpha));

  /* No rotation needed in this case */
  if ( ll < _MMG5_EPS ) {
    if ( n[2] > 0.0 ) {
      r[0][0] = 1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = 1.0;
    }
    else {
      r[0][0] = -1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = -1.0;
    }
  }
  else {
    l = sqrt(ll);

    r[0][0] = (aa*cosalpha + bb)/ll;
    r[0][1] = ab*(cosalpha-1)/ll;
    r[0][2] = -n[0]*sinalpha/l;
    r[1][0] = r[0][1];
    r[1][1] = (bb*cosalpha + aa)/ll;
    r[1][2] = -n[1]*sinalpha/l;
    r[2][0] = n[0]*sinalpha/l;
    r[2][1] = n[1]*sinalpha/l;
    r[2][2] = cosalpha;
  }

  /* apply rotation to vector GV */
  gv[0] =  r[0][0]*GV[0] + r[0][1]*GV[1] + r[0][2]*GV[2];
  gv[1] =  r[1][0]*GV[0] + r[1][1]*GV[1] + r[1][2]*GV[2];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt = &mesh->tria[iel];
    p1 = &mesh->point[pt->v[i1]];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
  }

  /* list goes modulo ilist */
  lispoi[3*ilist+1] =  lispoi[1];
  lispoi[3*ilist+2] =  lispoi[2];
  lispoi[3*ilist+3] =  lispoi[3];

  gv[0] = 0.0;
  gv[1] = 0.0;
  for (k=0; k<ilist; k++) {
    gv[0] += lispoi[3*k+1];
    gv[1] += lispoi[3*k+2];
  }
  gv[0] *= (1.0 / npt);
  gv[1] *= (1.0 / npt);

  /* Check all projections over tangent plane. */
  for (k=0; k<ilist-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 )  return(0);
  }
  det2d = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 )  return(0);

  /* Step 3 : locate new point in the ball, and compute its barycentric coordinates */
  det2d = lispoi[1]*gv[1] - lispoi[2]*gv[0];
  kel = 0;
  if ( det2d >= 0.0 ) {
    for (k=0; k<ilist; k++) {
      detloc = gv[0]*lispoi[3*(k+1)+2] - gv[1]*lispoi[3*(k+1)+1]; //orientation with respect to p2
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == ilist )  return(0);
  }
  else {
    for (k=ilist-1; k>=0; k--) {
      detloc = gv[1]*lispoi[3*k+1] - gv[0]*lispoi[3*k+2]; //orientation with respect to p2
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == -1 )  return(0);
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  iel = list[kel] / 3;
  det2d = -gv[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) \
    +  gv[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(det2d) < _MMG5_EPSD )  return(0);

  det2d = 1/det2d;
  step *= det2d;
  det2d = lispoi[3*(kel)+1]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]) - \
    lispoi[3*(kel)+2 ]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]);
  step *= det2d;
  step  = fabs(step);
  gv[0] *= step;
  gv[1] *= step;

  /* Computation of the barycentric coordinates of the new point in the corresponding triangle. */
  det2d = lispoi[3*kel+1]*lispoi[3*(kel+1)+2] - lispoi[3*kel+2]*lispoi[3*(kel+1)+1];
  if ( det2d < _MMG5_EPSD )  return(0);
  det2d = 1.0 / det2d;
  lambda[1] = lispoi[3*(kel+1)+2]*gv[0] - lispoi[3*(kel+1)+1]*gv[1];
  lambda[2] = -lispoi[3*(kel)+2]*gv[0] + lispoi[3*(kel)+1]*gv[1];
  lambda[1]*= (det2d);
  lambda[2]*= (det2d);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /* Step 4 : come back to original problem, and compute patch in triangle iel */
  iel  = list[kel]/3;
  i0 = list[kel]%3;
  i1 = _MMG5_inxt2[i0];
  i2 = _MMG5_inxt2[i1];
  pt = &mesh->tria[iel];

  ier = _MMG5_bezierCP(mesh,pt,&b,1);
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
  else {
    uv[0] = lambda[2];
    uv[1] = lambda[0];
  }
  ier = _MMGS_bezierInt(&b,uv,o,no,to);
  assert(ier);

  /* First test : check whether variance has been decreased */
  mlon = 0.0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1 = _MMG5_inxt2[i0];
    pt = &mesh->tria[iel];
    p1 = &mesh->point[pt->v[i1]];
    mlon += (p1->c[0]-o[0])*(p1->c[0]-o[0]) + (p1->c[1]-o[1])*(p1->c[1]-o[1]) \
      + (p1->c[2]-o[2])*(p1->c[2]-o[2]);
  }
  mlon = mlon / npt;

  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt  = &mesh->tria[iel];
    p1  = &mesh->point[pt->v[i1]];
    devmean = (p1->c[0]-o[0])*(p1->c[0]-o[0]) + (p1->c[1]-o[1])*(p1->c[1]-o[1]) \
      + (p1->c[2]-o[2])*(p1->c[2]-o[2]) - mlon;
    Vnew   += devmean*devmean;
  }
  Vnew  *= (1.0 / npt);
  /* if ( Vold < Vnew )  return(0); */

  /* Second test : check whether geometric approximation has not been too much degraded */
  ppt0 = &mesh->point[0];
  ppt0->tag = p0->tag;
  ppt0->c[0] = o[0];
  ppt0->c[1] = o[1];
  ppt0->c[2] = o[2];

  ppt0->n[0] = no[0];
  ppt0->n[1] = no[1];
  ppt0->n[2] = no[2];

  calold = calnew = DBL_MAX;
  for (k= 0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    pt  = &mesh->tria[iel];
    pt0 = &mesh->tria[0];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i0] = 0;
    caltmp = caleltsig_iso(mesh,NULL,iel);
    calold = MG_MIN(calold,caltmp);
    caltmp = caleltsig_iso(mesh,NULL,0);
    if ( caltmp < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,caltmp);
    /*if ( (calnew < BADKAL) && (calnew<=calold) )  return(0);
      if ( chkedg(mesh,0) )  return(0); */
  }
  if ( calold < NULKAL && calnew <= calold )      return(0);
  else if ( calnew < 0.3*calold )  return(0);

  /* Finally, update coordinates and normals of point, if new position is accepted : */
  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  p0->n[0] = no[0];
  p0->n[1] = no[1];
  p0->n[2] = no[2];

  return(1);
}

/* compute movement of a ridge point whose ball (consisting of triangles) is passed */
int movridpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist) {
  MMG5_pTria   pt,pt0;
  MMG5_pxPoint   go;
  _MMG5_Bezier  b;
  MMG5_pPoint  p0,p1,p2,ppt0;
  double  step,dd1,dd2,ddt,ps2,ll1old,ll1new,ll2old,ll2new,uv[2],o[3],*no1,*no2,*np1,*np2;
  double  nn1[3],nn2[3],psn11,psn12,to[3],calold,calnew,lam0,lam1,lam2;
  int     ier,k,iel,ip0,ip1,ip2,it1,it2;
  char    i0,i1,i2,isrid1,isrid2,isrid;

  step   = 0.1;
  isrid  = 0  ;
  isrid1 = 0  ;  isrid2 = 0;
  it1    = 0  ;  it2    = 0;
  ip1    = 0  ;  ip2    = 0;

  /* First step : make sure 2 ridge edges pass through point, and recoved them, with neighbouring triangles */
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
        if ( pt->tag[i1] & MG_GEO )  isrid1 = 1;
      }
      else if ( it1 && !it2 ) {
        if ( ip1 != pt->v[i2] ) {
          it2  = iel;
          ip2  = pt->v[i2]; // edge(i1) = (p0p2)
          if ( pt->tag[i1] & MG_GEO )  isrid2 = 1;
        }
      }
      else if ( it1 && it2 && (pt->v[i2] != ip1) && (pt->v[i2] != ip2) ) {
        //fprintf(stderr,"WRONG CONFIGURATION : 3 ridge edges landing on point %d\n",pt->v[i0]);
        return(0);
      }
    }

    if ( MG_EDG(pt->tag[i2]) ) {
      if ( !it1 ) {
        it1  = iel;
        ip1  = pt->v[i1]; // edge(i2) = (p0p1)
        if ( pt->tag[i2] & MG_GEO )  isrid1 = 1;
      }
      else if ( it1 && !it2 ) {
        if ( ip1 != pt->v[i1] ) {
          it2  = iel;
          ip2  = pt->v[i1]; // edge(i1) = (p0p2)
          if ( pt->tag[i2] & MG_GEO )  isrid2 = 1;
        }
      }
      else if ( it1 && it2 && (pt->v[i1] != ip1) && (pt->v[i1] != ip2) ) {
        //fprintf(stderr,"WRONG CONFIGURATION : 3 ridge edges landing on point %d\n",pt->v[i0]);
        return(0);
      }
    }
  }

  /* Second step : compute variance of lengths of edges (p0p1),(p0p2) and its gradient */
  iel = list[0] / 3;
  i0  = list[0] % 3;
  pt  = &mesh->tria[iel];
  ip0 = pt->v[i0];
  p0  = &mesh->point[ip0];
  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];

  ll1old = (p1->c[0] -p0->c[0])* (p1->c[0] -p0->c[0])  + (p1->c[1] -p0->c[1])* (p1->c[1] -p0->c[1]) \
    + (p1->c[2] -p0->c[2])* (p1->c[2] -p0->c[2]);
  ll2old = (p2->c[0] -p0->c[0])* (p2->c[0] -p0->c[0])  + (p2->c[1] -p0->c[1])* (p2->c[1] -p0->c[1]) \
    + (p2->c[2] -p0->c[2])* (p2->c[2] -p0->c[2]);

  /* Third step : infer arc length of displacement, parameterized over edges */
  if ( ll2old > ll1old ) {  /* move is made towards p2 */
    isrid = isrid2;
    pt = &mesh->tria[it2];

    ier = _MMG5_bezierCP(mesh,pt,&b,1);
    assert(ier);

    /* fill table uv with barycentric coor in the tria */
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
        uv[1] = 1.0 - step;
      }
    }
    else {
      if ( pt->v[1] == ip0 ) {
        uv[0] = 1.0 - step;
        uv[1] = step;
      }
      else if ( pt->v[2] == ip0 ) {
        uv[0] = step;
        uv[1] = 1.0 - step;
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
        uv[1] = 1.0 - step;
      }
    }
    else {
      if ( pt->v[1] == ip0 ) {
        uv[0] = 1.0 - step;
        uv[1] = step;
      }
      else if ( pt->v[2] == ip0 ) {
        uv[0] = step;
        uv[1] = 1.0 - step;
      }
    }
    ier = _MMGS_bezierInt(&b,uv,o,nn1,to);
    assert(ier);
  }

  /* check whether proposed move is admissible */
  ll1new = (p1->c[0] -o[0])* (p1->c[0] -o[0])  + (p1->c[1] -o[1])* (p1->c[1] -o[1]) \
    + (p1->c[2] -o[2])* (p1->c[2] -o[2]);
  ll2new = (p2->c[0] -o[0])* (p2->c[0] -o[0])  + (p2->c[1] -o[1])* (p2->c[1] -o[1]) \
    + (p2->c[2] -o[2])* (p2->c[2] -o[2]);
  if( fabs(ll2new -ll1new) >= fabs(ll2old -ll1old) )  return(0);

  /* normal and tangent updates */
  // Bezier basis function of order 2
  lam0 = (1.0-step)*(1.0-step);
  lam1 = 2.0*step*(1.0-step);
  lam2 = step*step;

  if ( ll2old > ll1old ) {  //move is made towards p2
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

    if ( fabs(1.0-fabs(psn11)) < fabs(1.0-fabs(psn12)) ){    //no1 goes with np1, no2 with np2
      nn1[0] = no1[0]+np1[0];
      nn1[1] = no1[1]+np1[1];
      nn1[2] = no1[2]+np1[2];

      nn2[0] = no2[0]+np2[0];
      nn2[1] = no2[1]+np2[1];
      nn2[2] = no2[2]+np2[2];

      // p0p2 . nn1
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

      /* At this point, we got the control points associated to the interpolaton formula for normals.*/
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
    else {     /* no1 goes with np2 and no2 with np1 */
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

      /* At this point, we got the control points associated to the interpolation formula for normals.*/
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
  }
  else {   /* move along p1 */
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

    if ( fabs(1.0-fabs(psn11)) < fabs(1.0-fabs(psn12)) ) { /* no1 goes with np1, no2 with np2 */
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
    else {   /* no1 goes with np2 and no2 with np1 */
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
  }

  /* check proposed motion */
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

  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    pt  = &mesh->tria[iel];
    pt0 = &mesh->tria[0];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i0] = 0;
    calold = caleltsig_iso(mesh,NULL,iel);
    calnew = caleltsig_iso(mesh,NULL,0);
    if ( (calnew < 0.001) && (calnew<calold) )  return(0);
    //if ( chkedg(mesh,0) )  return(0);
  }

  /* coordinates, normals, tangents update */
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

  return(1);
}

