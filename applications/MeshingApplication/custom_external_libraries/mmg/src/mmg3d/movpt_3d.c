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
 * \file mmg3d/movpt_3d.c
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
 * than 0.9 of the old minimum element quality.
 *
 * \return 0 if we can't move the point, 1 if we can.
 *
 * Move internal point whose volumic is passed.
 *
 * \remark the metric is not interpolated at the new position.
 * \remark we don't check if we break the hausdorff criterion.
 *
 */
int _MMG5_movintpt_iso(MMG5_pMesh mesh,MMG5_pSol met, _MMG3D_pOctree octree,
                       int *list,int ilist,int improve) {
  MMG5_pTetra               pt,pt0;
  MMG5_pPoint               p0,p1,p2,p3,ppt0;
  double               vol,totvol;
  double               calold,calnew,*callist;
  int                  k,iel,i0;

  // Dynamic alloc for windows comptibility
  _MMG5_SAFE_MALLOC(callist, ilist, double);

  pt0    = &mesh->tetra[0];
  ppt0   = &mesh->point[0];
  memset(ppt0,0,sizeof(MMG5_Point));

  iel = list[0] / 4;
  i0  = list[0] % 4;

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
  else if ( improve && calnew < 1.02 * calold ) {
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
    (&mesh->tetra[list[k]/4])->mark=mesh->mark;
  }


  _MMG5_SAFE_FREE(callist);
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param list pointer toward the volumic ball of the point.
 * \param ilist size of the volumic ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 0.9 of the old minimum element quality.
 *
 * \return 0 if we can't move the point, 1 if we can.
 *
 * Move internal point whose volumic ball is passed (for LES optimization).  The
 * optimal point position is computed as the barycenter of the optimal point
 * position for each tetra. The optimal point position for a tetra is the point
 * located over the normal of the face at the face barycenter and at the
 * distance 1 of the face.
 *
 * \remark the metric is not interpolated at the new position.
 * \remark we don't check if we break the hausdorff criterion.
 *
 */
int _MMG5_movintptLES_iso(MMG5_pMesh mesh,MMG5_pSol met, _MMG3D_pOctree octree,
                          int *list,int ilist,int improve) {
  MMG5_pTetra               pt,pt0;
  MMG5_pPoint               p0,p1,p2,p3,ppt0;
  double               vol,totvol;
  double               calold,calnew,*callist;
  double               x21,y21,z21,x31,y31,z31,nx,ny,nz,bary[3],dd,len;
  double               u10[3],u20[3],u30[3],oldc[3],coe;
  int                  k,iel,ifac,iter,maxtou;

  // Dynamic alloc for windows comptibility
  _MMG5_SAFE_MALLOC(callist, ilist, double);

  pt0    = &mesh->tetra[0];
  ppt0   = &mesh->point[0];
  memset(ppt0,0,sizeof(MMG5_Point));

  /* Coordinates of optimal point */
  calold = DBL_MAX;
  totvol = 0.0;
  for (k=0; k<ilist; k++) {
    iel  = list[k] / 4;
    ifac = list[k] % 4;
    pt = &mesh->tetra[iel];

    p0 = &mesh->point[pt->v[ifac]];
    memcpy(oldc,p0->c,3*sizeof(double));

    p1 = &mesh->point[pt->v[_MMG5_idir[ifac][0]]];
    p2 = &mesh->point[pt->v[_MMG5_idir[ifac][1]]];
    p3 = &mesh->point[pt->v[_MMG5_idir[ifac][2]]];

    vol= _MMG5_det4pt(p0->c,p1->c,p2->c,p3->c);
    totvol += vol;

    /* local optimal point : length 1 over the normal of the face */
    x21 = p2->c[0]-p1->c[0];
    y21 = p2->c[1]-p1->c[1];
    z21 = p2->c[2]-p1->c[2];

    x31 = p3->c[0]-p1->c[0];
    y31 = p3->c[1]-p1->c[1];
    z31 = p3->c[2]-p1->c[2];

    nx = (y31*z21 - z31*y21);
    ny = (z31*x21 - x31*z21);
    nz = (x31*y21 - y31*x21);

    dd = sqrt(nx*nx+ny*ny+nz*nz);
    dd = 1./dd;
    nx *= dd;
    ny *= dd;
    nz *= dd;

    u10[0] = p1->c[0]-p0->c[0];
    u10[1] = p1->c[1]-p0->c[1];
    u10[2] = p1->c[2]-p0->c[2];

    u20[0] = p2->c[0]-p0->c[0];
    u20[1] = p2->c[1]-p0->c[1];
    u20[2] = p2->c[2]-p0->c[2];

    u30[0] = p3->c[0]-p0->c[0];
    u30[1] = p3->c[1]-p0->c[1];
    u30[2] = p3->c[2]-p0->c[2];

// mmg3d4
    len =  sqrt(u10[0]*u10[0]+u10[1]*u10[1]+u10[2]*u10[2])/met->m[pt->v[_MMG5_idir[ifac][0]]]
      + sqrt(u20[0]*u20[0]+u20[1]*u20[1]+u20[2]*u20[2])/met->m[pt->v[_MMG5_idir[ifac][1]]]
      + sqrt(u30[0]*u30[0]+u30[1]*u30[1]+u30[2]*u30[2])/met->m[pt->v[_MMG5_idir[ifac][2]]];

    len /= 3.;

    len = 1./len;

    /* face barycenter */
    bary[0] = (p1->c[0]+p2->c[0]+p3->c[0])/3.;
    bary[1] = (p1->c[1]+p2->c[1]+p3->c[1])/3.;
    bary[2] = (p1->c[2]+p2->c[2]+p3->c[2])/3.;


// vecteur bary_face->bary_tet projeté sur la normale
    /* len = (0.25*(p0->c[0]+p1->c[0]+p2->c[0]+p3->c[0]) - bary[0])*nx */
    /*   +  (0.25*(p0->c[1]+p1->c[1]+p2->c[1]+p3->c[1]) - bary[1])*ny */
    /*   +  (0.25*(p0->c[2]+p1->c[2]+p2->c[2]+p3->c[2]) - bary[2])*nz; */

    /* len = (p0->c[0]-bary[0])*nx +  (p0->c[1]-bary[1])*ny + (p0->c[2]-bary[2])*nz; */


    /* Optimal point: barycenter of local optimal points */
    vol = 1;
    ppt0->c[0] += vol* ( bary[0] + nx * len );
    ppt0->c[1] += vol* ( bary[1] + ny * len );
    ppt0->c[2] += vol* ( bary[2] + nz * len );
    calold = MG_MIN(calold, pt->qual);
  }
  if (totvol < _MMG5_EPSD2) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* totvol = 1.0 / totvol; */
  /* ppt0->c[0] *= totvol; */
  /* ppt0->c[1] *= totvol; */
  /* ppt0->c[2] *= totvol; */
  ppt0->c[0] *= 1./(double) ilist;
  ppt0->c[1] *= 1./(double) ilist;
  ppt0->c[2] *= 1./(double) ilist;

   coe = 0.9;
   maxtou = 10;
   iter = 0;
   do {
     p0->c[0] = (1. - coe) *oldc[0] + coe * ppt0->c[0] ;
     p0->c[1] = (1. - coe) *oldc[1] + coe * ppt0->c[1];
     p0->c[2] = (1. - coe) *oldc[2] + coe * ppt0->c[2];

     /* Check new position validity */
     calnew = DBL_MAX;
     for (k=0; k<ilist; k++) {
       iel = list[k] / 4;
       pt  = &mesh->tetra[iel];
       memcpy(pt0,pt,sizeof(MMG5_Tetra));
       callist[k] = _MMG5_caltet(mesh,met,pt0);//_MMG5_orcal(mesh,met,0);
       if (calold < _MMG5_NULKAL && callist[k] <= calold) {
         break;
       } else if ((callist[k] < _MMG5_EPSD2) || (callist[k] < _MMG5_NULKAL)) {
         break;
       } else if( callist[k] < 0.1) {
         if (callist[k] < 1.01*calold) break;
       } else {
         if (callist[k] < 1.02*calold) break;
       }
       calnew = MG_MIN(calnew,callist[k]);
     }
     if(k<ilist) {
       coe *= 0.5;
       continue;
     }
     break;
   } while(++iter <= maxtou );

   if ( iter > maxtou ) {
     memcpy(p0->c,oldc,3*sizeof(double));
     _MMG5_SAFE_FREE(callist);
     return(0);
   }

   /* update position */
  if ( octree )
    _MMG3D_moveOctree(mesh, octree, mesh->tetra[list[0]/4].v[list[0]%4], ppt0->c, p0->c);

   for (k=0; k<ilist; k++) {
     (&mesh->tetra[list[k]/4])->qual=callist[k];
     (&mesh->tetra[list[k]/4])->mark=mesh->mark;
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
 * Move boundary regular point, whose volumic and surfacic balls are passed.
 *
 * \remark the metric is not interpolated at the new position.
 */
int _MMG5_movbdyregpt_iso(MMG5_pMesh mesh, MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv,int *lists,int ilists,
                          int improve) {
  MMG5_pTetra       pt,pt0;
  MMG5_pxTetra      pxt;
  MMG5_pPoint       p0,p1,p2,ppt0;
  MMG5_Tria         tt;
  MMG5_pxPoint      pxp;
  _MMG5_Bezier      b;
  double            *n,r[3][3],lispoi[3*MMG3D_LMAX+1],ux,uy,uz,det2d;
  double            detloc,oppt[2],step,lambda[3];
  double            ll,m[2],uv[2],o[3],no[3],to[3];
  double            calold,calnew,caltmp,*callist;
  int               k,kel,iel,l,n0,na,nb,ntempa,ntempb,ntempc,nut,nxp;
  unsigned char     i0,iface,i;

  step = 0.1;
  nut    = 0;
  oppt[0] = 0.0;
  oppt[1] = 0.0;
  if ( ilists < 2 )      return(0);

  k      = listv[0] / 4;
  i0 = listv[0] % 4;
  pt = &mesh->tetra[k];
  n0 = pt->v[i0];
  p0 = &mesh->point[n0];
  assert( p0->xp && !MG_EDG(p0->tag) );

  n = &(mesh->xpoint[p0->xp].n1[0]);

  /** Step 1 : rotation matrix that sends normal n to the third coordinate vector of R^3 */
  if ( !_MMG5_rotmatrix(n,r) ) return(0);

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
    if ( det2d < 0.0 )  return(0);
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 )    return(0);

  /** Step 3 : Compute optimal position to make current triangle equilateral, and average of
      these positions*/
  for (k=0; k<ilists; k++) {
    m[0] = 0.5*(lispoi[3*(k+1)+1] + lispoi[3*k+1]);
    m[1] = 0.5*(lispoi[3*(k+1)+2] + lispoi[3*k+2]);
    ux = lispoi[3*(k+1)+1] - lispoi[3*k+1];
    uy = lispoi[3*(k+1)+2] - lispoi[3*k+2];
    ll = ux*ux + uy*uy;
    if ( ll < _MMG5_EPSD )    continue;
    nut++;
    oppt[0] += (m[0]-_MMG5_SQR32*uy);
    oppt[1] += (m[1]+_MMG5_SQR32*ux);
  }
  oppt[0] *= (1.0 / nut);
  oppt[1] *= (1.0 / nut);

  /** Step 4 : locate new point in the ball, and compute its barycentric coordinates */
  det2d = lispoi[1]*oppt[1] - lispoi[2]*oppt[0];
  kel = 0;
  if ( det2d >= 0.0 ) {
    for (k=0; k<ilists; k++) {
      detloc = oppt[0]*lispoi[3*(k+1)+2] - oppt[1]*lispoi[3*(k+1)+1];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == ilists ) return(0);
  }
  else {
    for (k=ilists-1; k>=0; k--) {
      detloc = lispoi[3*k+1]*oppt[1] - lispoi[3*k+2]*oppt[0];
      if ( detloc >= 0.0 ) {
        kel = k;
        break;
      }
    }
    if ( k == -1 ) return(0);
  }

  /* Sizing of time step : make sure point does not go out corresponding triangle. */
  det2d = -oppt[1]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]) + \
    oppt[0]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]);
  if ( fabs(det2d) < _MMG5_EPSD ) return(0);

  det2d = 1.0 / det2d;
  step *= det2d;

  det2d = lispoi[3*(kel)+1]*(lispoi[3*(kel+1)+2] - lispoi[3*(kel)+2]) - \
    lispoi[3*(kel)+2 ]*(lispoi[3*(kel+1)+1] - lispoi[3*(kel)+1]);
  step *= det2d;
  step  = fabs(step);
  oppt[0] *= step;
  oppt[1] *= step;

  /* Computation of the barycentric coordinates of the new point in the corresponding triangle. */
  det2d = lispoi[3*kel+1]*lispoi[3*(kel+1)+2] - lispoi[3*kel+2]*lispoi[3*(kel+1)+1];
  if ( det2d < _MMG5_EPSD )    return(0);
  det2d = 1.0 / det2d;
  lambda[1] = lispoi[3*(kel+1)+2]*oppt[0] - lispoi[3*(kel+1)+1]*oppt[1];
  lambda[2] = -lispoi[3*(kel)+2]*oppt[0] + lispoi[3*(kel)+1]*oppt[1];
  lambda[1]*= (det2d);
  lambda[2]*= (det2d);
  lambda[0] = 1.0 - lambda[1] - lambda[2];

  /** Step 5 : come back to original problem, and compute patch in triangle iel */
  iel    = lists[kel] / 4;
  iface  = lists[kel] % 4;
  pt     = &mesh->tetra[iel];
  pxt    = &mesh->xtetra[pt->xt];

  _MMG5_tet2tri(mesh,iel,iface,&tt);

  if(!_MMG5_bezierCP(mesh,&tt,&b,MG_GET(pxt->ori,iface))){
    fprintf(stderr,"%s:%d: Error: function _MMG5_bezierCP return 0\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
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
  if(!_MMG3D_bezierInt(&b,uv,o,no,to)){
    fprintf(stderr,"%s:%d: Error: function _MMG3D_bezierInt return 0\n",
            __FILE__,__LINE__);
    exit(EXIT_FAILURE);
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
                       return(0));
    n = &(mesh->xpoint[p0->xp].n1[0]);
  }
  ppt0->xp = nxp;
  pxp = &mesh->xpoint[nxp];
  memcpy(pxp,&(mesh->xpoint[p0->xp]),sizeof(MMG5_xPoint));
  pxp->n1[0] = no[0];
  pxp->n1[1] = no[1];
  pxp->n1[2] = no[2];

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
    if ( caltmp < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )    return(0);
  else if (calnew < _MMG5_NULKAL) return(0);
  else if (improve && calnew < 1.02*calold) return(0);
  else if ( calnew < 0.3*calold )        return(0);
  memset(pxp,0,sizeof(MMG5_xPoint));

  /* Test : check whether all volumes remain positive with new position of the point */

  // Dynamic allocations for windows compatibility
  _MMG5_SAFE_MALLOC(callist, ilistv, double);

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
  if (callist[l] < _MMG5_EPSD) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }
    calnew = MG_MIN(calnew,callist[l]);
  }
  if (calold < _MMG5_NULKAL && calnew <= calold) {
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
  else if (calnew < 0.3*calold) {
    _MMG5_SAFE_FREE(callist);
    return(0);
  }

  /* When all tests have been carried out, update coordinates and normals */
  if ( octree )
    _MMG3D_moveOctree(mesh, octree, n0, o, p0->c);

  p0->c[0] = o[0];
  p0->c[1] = o[1];
  p0->c[2] = o[2];

  n[0] = no[0];
  n[1] = no[1];
  n[2] = no[2];

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual= callist[l];
    (&mesh->tetra[listv[l]/4])->mark=mesh->mark;
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
 * Move boundary reference point, whose volumic and surfacic balls are passed.
 *
 * \remark the metric is not interpolated at the new position.
 */
int _MMG5_movbdyrefpt_iso(MMG5_pMesh mesh, MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv, int *lists, int ilists,
                          int improve){
  MMG5_pTetra           pt,pt0;
  MMG5_pxTetra          pxt;
  MMG5_pPoint           p0,p1,p2,ppt0;
  MMG5_Tria             tt;
  MMG5_pxPoint          pxp;
  MMG5_pPar             par;
  double                step,ll1old,ll2old,o[3],no[3],to[3];
  double                calold,calnew,caltmp,*callist,hmax,hausd;
  int                   l,iel,ip0,ipa,ipb,iptmpa,iptmpb,it1,it2,ip1,ip2,ip,nxp;
  int                   isloc,j;
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
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  ll1old = (p1->c[0] -p0->c[0])* (p1->c[0] -p0->c[0]) \
    + (p1->c[1] -p0->c[1])* (p1->c[1] -p0->c[1])      \
    + (p1->c[2] -p0->c[2])* (p1->c[2] -p0->c[2]);
  ll2old = (p2->c[0] -p0->c[0])* (p2->c[0] -p0->c[0]) \
    + (p2->c[1] -p0->c[1])* (p2->c[1] -p0->c[1])      \
    + (p2->c[2] -p0->c[2])* (p2->c[2] -p0->c[2]);

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

  /* For each surface triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilists ; l++ ){
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    pxt         = &mesh->xtetra[pt->xt];
    _MMG5_tet2tri(mesh,iel,iface,&tt);
    calold = MG_MIN(calold,_MMG5_caltri(mesh,met,&tt));
    for( i=0 ; i<3 ; i++ )
      if ( tt.v[i] == ip0 )      break;
    assert(i<3);
    tt.v[i] = 0;
    caltmp = _MMG5_caltri(mesh,met,&tt);
    if ( caltmp < _MMG5_EPSD )        return(0);
    calnew = MG_MIN(calnew,caltmp);

    /* Local parameters for tt and iel */
    hmax  = mesh->info.hmax;
    hausd = mesh->info.hausd;

    isloc = 0;
    if ( mesh->info.parTyp & MG_Tetra ) {
      for ( j=0; j<mesh->info.npar; ++j ) {
        par = &mesh->info.par[j];

        if ( par->elt != MMG5_Tetrahedron )  continue;
        if ( par->ref != pt->ref ) continue;

        hmax = par->hmax;
        hausd = par->hausd;
        isloc = 1;
        break;
      }
    }
    if ( mesh->info.parTyp & MG_Tria ) {
      if ( isloc ) {
        for ( j=0; j<mesh->info.npar; ++j ) {
          par = &mesh->info.par[j];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax = MG_MIN(hmax,par->hmax);
          hausd = MG_MIN(hausd,par->hausd);
          break;
        }
      }
      else {
        for ( j=0; j<mesh->info.npar; ++j ) {
          par = &mesh->info.par[j];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax  = par->hmax;
          hausd = par->hausd;
          isloc = 1;
          break;
        }
      }
    }

    if ( _MMG5_chkedg(mesh,&tt,MG_GET(pxt->ori,iface),hmax,hausd,isloc) ) {
      memset(pxp,0,sizeof(MMG5_xPoint));
      return(0);
    }
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

  for( l=0 ; l<ilistv ; l++ ){
    (&mesh->tetra[listv[l]/4])->qual = callist[l];
    (&mesh->tetra[listv[l]/4])->mark = mesh->mark;
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
 * \remark the metric is not interpolated at the new position.
 */
int _MMG5_movbdynompt_iso(MMG5_pMesh mesh,MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv, int *lists, int ilists,
                          int improve){
  MMG5_pTetra       pt,pt0;
  MMG5_pxTetra      pxt;
  MMG5_pPoint       p0,p1,p2,ppt0;
  MMG5_pxPoint      pxp;
  MMG5_Tria         tt;
  MMG5_pPar         par;
  double            step,ll1old,ll2old,calold,calnew,caltmp,*callist;
  double            o[3],no[3],to[3],hmax,hausd;
  int               ip0,ip1,ip2,ip,iel,ipa,ipb,l,iptmpa,iptmpb,it1,it2,nxp;
  int               j,isloc;
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

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  ll1old = (p1->c[0] -p0->c[0])* (p1->c[0] -p0->c[0]) \
    + (p1->c[1] -p0->c[1])* (p1->c[1] -p0->c[1])      \
    + (p1->c[2] -p0->c[2])* (p1->c[2] -p0->c[2]);
  ll2old = (p2->c[0] -p0->c[0])* (p2->c[0] -p0->c[0]) \
    + (p2->c[1] -p0->c[1])* (p2->c[1] -p0->c[1])      \
    + (p2->c[2] -p0->c[2])* (p2->c[2] -p0->c[2]);

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

  /* For each surface triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for( l=0 ; l<ilists ; l++ ){
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    pxt         = &mesh->xtetra[pt->xt];
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


    /* Local parameters for tt and iel */
    hmax  = mesh->info.hmax;
    hausd = mesh->info.hausd;

    isloc = 0;
    if ( mesh->info.parTyp & MG_Tetra ) {
      for ( j=0; j<mesh->info.npar; ++j ) {
        par = &mesh->info.par[j];

        if ( par->elt != MMG5_Tetrahedron )  continue;
        if ( par->ref != pt->ref ) continue;

        hmax = par->hmax;
        hausd = par->hausd;
        isloc = 1;
        break;
      }
    }
    if ( mesh->info.parTyp & MG_Tria ) {
      if ( isloc ) {
        for ( j=0; j<mesh->info.npar; ++j ) {
          par = &mesh->info.par[j];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax = MG_MIN(hmax,par->hmax);
          hausd = MG_MIN(hausd,par->hausd);
          break;
        }
      }
      else {
        for ( j=0; j<mesh->info.npar; ++j ) {
          par = &mesh->info.par[j];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax  = par->hmax;
          hausd = par->hausd;
          isloc = 1;
          break;
        }
      }
    }

    if ( _MMG5_chkedg(mesh,&tt,MG_GET(pxt->ori,iface),hmax,hausd,isloc) ) {
      memset(pxp,0,sizeof(MMG5_xPoint));
      return(0);
    }
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

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual = callist[l];
    (&mesh->tetra[listv[l]/4])->mark = mesh->mark;
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
 * Move boundary ridge point, whose volumic and surfacic balls are passed.
 *
 */
int _MMG5_movbdyridpt_iso(MMG5_pMesh mesh, MMG5_pSol met, _MMG3D_pOctree octree, int *listv,
                          int ilistv,int *lists,int ilists,
                          int improve) {
  MMG5_pTetra          pt,pt0;
  MMG5_pxTetra         pxt;
  MMG5_pPoint          p0,p1,p2,ppt0;
  MMG5_Tria            tt;
  MMG5_pxPoint         pxp;
  MMG5_pPar            par;
  double               step,ll1old,ll2old,o[3],no1[3],no2[3],to[3];
  double               calold,calnew,caltmp,*callist,hmax,hausd;
  int                  l,iel,ip0,ipa,ipb,iptmpa,iptmpb,it1,it2,ip1,ip2,ip,nxp;
  int                  j,isloc;
  int16_t              tag;
  unsigned char        i,i0,ie,iface,iface1,iface2,iea,ieb,ie1,ie2;

  step = 0.1;
  ip1 = ip2 = 0;
  pt    = &mesh->tetra[listv[0]/4];
  ip0 = pt->v[listv[0]%4];
  p0    = &mesh->point[ip0];

  assert ( MG_GEO & p0->tag );

  /* Travel surfacic ball an recover the two ending points of ridge : two senses must be used
     POSSIBLE OPTIMIZATION HERE : One travel only is needed */
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
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  ll1old = (p1->c[0] -p0->c[0])* (p1->c[0] -p0->c[0]) \
    + (p1->c[1] -p0->c[1])* (p1->c[1] -p0->c[1])      \
    + (p1->c[2] -p0->c[2])* (p1->c[2] -p0->c[2]);
  ll2old = (p2->c[0] -p0->c[0])* (p2->c[0] -p0->c[0]) \
    + (p2->c[1] -p0->c[1])* (p2->c[1] -p0->c[1])      \
    + (p2->c[2] -p0->c[2])* (p2->c[2] -p0->c[2]);

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

  /* For each surfacic triangle, build a virtual displaced triangle for check purposes */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilists; l++) {
    iel         = lists[l] / 4;
    iface       = lists[l] % 4;
    pt          = &mesh->tetra[iel];
    pxt         = &mesh->xtetra[pt->xt];
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

    /* Local parameters for tt and iel */
    hmax  = mesh->info.hmax;
    hausd = mesh->info.hausd;

    isloc = 0;
    if ( mesh->info.parTyp & MG_Tetra ) {
      for ( j=0; j<mesh->info.npar; ++j ) {
        par = &mesh->info.par[j];

        if ( par->elt != MMG5_Tetrahedron )  continue;
        if ( par->ref != pt->ref ) continue;

        hmax = par->hmax;
        hausd = par->hausd;
        isloc = 1;
        break;
      }
    }
    if ( mesh->info.parTyp & MG_Tria ) {
      if ( isloc ) {
        for ( j=0; j<mesh->info.npar; ++j ) {
          par = &mesh->info.par[j];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax = MG_MIN(hmax,par->hmax);
          hausd = MG_MIN(hausd,par->hausd);
          break;
        }
      }
      else {
        for ( j=0; j<mesh->info.npar; ++j ) {
          par = &mesh->info.par[j];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax  = par->hmax;
          hausd = par->hausd;
          isloc = 1;
          break;
        }
      }
    }

    if ( _MMG5_chkedg(mesh,&tt,MG_GET(pxt->ori,iface),hmax,hausd,isloc) ) {
      memset(pxp,0,sizeof(MMG5_xPoint));
      return(0);
    }
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

  for(l=0; l<ilistv; l++){
    (&mesh->tetra[listv[l]/4])->qual = callist[l];
    (&mesh->tetra[listv[l]/4])->mark = mesh->mark;
  }
  _MMG5_SAFE_FREE(callist);
  return(1);
}


int _MMG3D_movv_ani(MMG5_pMesh mesh,MMG5_pSol sol,int k,int ib) {
  MMG5_pTetra   pt,pt1;
  MMG5_pPoint   ppa,ppb,p1,p2,p3;
  int           j,iadr,ipb,iter,maxiter,l,lon,iel,i1,i2,i3,list[MMG3D_LMAX+2];
  double        *mp,coe,qualtet[MMG3D_LMAX+2];
  double        ax,ay,az,bx,by,bz,nx,ny,nz,dd,len,qual,oldc[3];
  assert(k);
  assert(ib<4);
  pt = &mesh->tetra[k];
  ppa  = &mesh->point[pt->v[ib]];
  if(ppa->tag & MG_BDY) return(0);
  iadr = pt->v[ib]*sol->size + 0;
  mp   = &sol->m[iadr];

  /*compute normal*/
  i1 = pt->v[_MMG5_idir[ib][0]];
  i2 = pt->v[_MMG5_idir[ib][1]];
  i3 = pt->v[_MMG5_idir[ib][2]];
  p1 = &mesh->point[i1];
  p2 = &mesh->point[i2];
  p3 = &mesh->point[i3];

  ax = p3->c[0] - p1->c[0];
  ay = p3->c[1] - p1->c[1];
  az = p3->c[2] - p1->c[2];

  bx = p2->c[0] - p1->c[0];
  by = p2->c[1] - p1->c[1];
  bz = p2->c[2] - p1->c[2];

  nx = (ay*bz - az*by);
  ny = (az*bx - ax*bz);
  nz = (ax*by - ay*bx);

  dd = sqrt(nx*nx+ny*ny+nz*nz);
  dd = 1./dd;
  nx *=dd;
  ny *=dd;
  nz *=dd;
  len = 0;
  for (j=0; j<3; j++) {
    ipb = pt->v[ _MMG5_idir[ib][j] ];
    ppb = &mesh->point[ipb];

    ax  = ppb->c[0] - ppa->c[0];
    ay  = ppb->c[1] - ppa->c[1];
    az  = ppb->c[2] - ppa->c[2];

    dd =       mp[0]*ax*ax + mp[3]*ay*ay + mp[5]*az*az \
      + 2.0*(mp[1]*ax*ay + mp[2]*ax*az + mp[4]*ay*az);
    assert(dd>0);
    len += sqrt(dd);
  }

  dd  = 1.0 / (double)3.;
  len = 1.0 / len;
  len *= dd;
  memcpy(oldc,ppa->c,3*sizeof(double));

  lon = _MMG5_boulevolp(mesh,k,ib,&list[0]);
  if(mesh->info.imprim < 0 ) if(lon < 4 && lon) printf("lon petit : %d\n",lon);
  if(!lon) return(0);

  coe     = 1.;
  iter    = 0;
  maxiter = 20;
  do {
    ppa->c[0] = oldc[0] + coe * nx * len;
    ppa->c[1] = oldc[1] + coe * ny * len;
    ppa->c[2] = oldc[2] + coe * nz * len;

    for (l=0; l<lon; l++) {
      iel = list[l] / 4 ;
      pt1 = &mesh->tetra[iel];

      qual = _MMG5_caltet(mesh,sol,pt1);
      /*warning if we increase the coefficient (ex 1.4), the mesh quality becomes poor very quickly*/
      if ( qual*1.01 <= pt1->qual) break;
      qualtet[l] = qual;

    }
    if ( l >= lon )  break;
    coe *= 0.5;
  }
  while ( ++iter <= maxiter );
  if ( iter > maxiter) {
    memcpy(ppa->c,oldc,3*sizeof(double));
    return(0);
  }

  for (l=0; l<lon; l++) {
    iel = list[l] / 4;
    pt1 = &mesh->tetra[iel];
    pt1->qual = qualtet[l];
    pt1->mark = mesh->mark;
    //    if ( pt1->qual < declic )
    //  MMG_kiudel(queue,iel);
  }
  return(1);

}


int _MMG3D_movv_iso(MMG5_pMesh mesh,MMG5_pSol sol,int k,int ib) {
  MMG5_pTetra pt,pt1;
  MMG5_pPoint ppa,ppb,p1,p2,p3;
  int    j,iadr,ipb,iter,maxiter,l,lon,iel,i1,i2,i3,list[MMG3D_LMAX+2];;
  double  hp,coe,crit,qualtet[MMG3D_LMAX+2];;
  double ax,ay,az,bx,by,bz,nx,ny,nz,dd,len,qual,oldc[3];

  assert(k);
  assert(ib<4);
  pt = &mesh->tetra[k];

  ppa  = &mesh->point[pt->v[ib]];
  if(ppa->tag & MG_BDY) return(0);

  iadr = (pt->v[ib])*sol->size;
  hp   = sol->m[iadr];

  /*compute normal*/
  i1 = pt->v[_MMG5_idir[ib][0]];
  i2 = pt->v[_MMG5_idir[ib][1]];
  i3 = pt->v[_MMG5_idir[ib][2]];
  p1 = &mesh->point[i1];
  p2 = &mesh->point[i2];
  p3 = &mesh->point[i3];

  ax = p3->c[0] - p1->c[0];
  ay = p3->c[1] - p1->c[1];
  az = p3->c[2] - p1->c[2];

  bx = p2->c[0] - p1->c[0];
  by = p2->c[1] - p1->c[1];
  bz = p2->c[2] - p1->c[2];

  nx = (ay*bz - az*by);
  ny = (az*bx - ax*bz);
  nz = (ax*by - ay*bx);

  dd = sqrt(nx*nx+ny*ny+nz*nz);
  dd = 1./dd;
  nx *=dd;
  ny *=dd;
  nz *=dd;
  len = 0;
  for (j=0; j<3; j++) {
    ipb = pt->v[ _MMG5_idir[ib][j] ];
    ppb = &mesh->point[ipb];

    ax  = ppb->c[0] - ppa->c[0];
    ay  = ppb->c[1] - ppa->c[1];
    az  = ppb->c[2] - ppa->c[2];

    dd    =   sqrt(ax*ax +ay*ay +az*az);
    len  +=   dd/hp;
  }

  dd  = 1.0 / (double)3.;
  len *= dd;
  if(len > 0.) len = 1.0 / len;
  else printf("MMG_movevertex len %e\n",len);

  memcpy(oldc,ppa->c,3*sizeof(double));

  lon = _MMG5_boulevolp(mesh,k,ib,&list[0]);
  if(mesh->info.imprim < 0) if(lon < 4 && lon) printf("lon petit : %d\n",lon);
  if(!lon) return(0);

  /*qual crit*/
  crit = pt->qual;
  for (l=1; l<lon; l++) {
    iel = list[l] / 4;
    pt1 = &mesh->tetra[iel];
    if ( crit > pt1->qual )  crit = pt1->qual;
  }
  /* if ( (crit > 100/ALPHAD) ) {
     crit *= 1.1;
     } else
  */ crit *= 1.01;

  coe     = 1.;
  iter    = 0;
  maxiter = 20;
  do {

    ppa->c[0] = oldc[0] + coe * nx * len;
    ppa->c[1] = oldc[1] + coe * ny * len;
    ppa->c[2] = oldc[2] + coe * nz * len;
    for (l=0; l<lon; l++) {
      iel = list[l] / 4;
      pt1 = &mesh->tetra[iel];
      qual = _MMG5_caltet(mesh,sol,pt1);
      if ( qual < crit )  break;
      qualtet[l] = qual;

    }
    if ( l >= lon )  break;
    coe *= 0.5;
  }
  while ( ++iter <= maxiter );
  if ( iter > maxiter) {
    memcpy(ppa->c,oldc,3*sizeof(double));
    return(0);
  }

  for (l=0; l<lon; l++) {
    iel = list[l] / 4;
    pt1 = &mesh->tetra[iel];
    pt1->qual = qualtet[l];
    pt1->mark = mesh->mark;
    //    if ( pt1->qual < declic )
    //  MMG_kiudel(queue,iel);
  }
  return(1);

}
