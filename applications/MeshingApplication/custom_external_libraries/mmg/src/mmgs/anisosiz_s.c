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
 * \file mmgs/anisosiz_s.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"
#include "inlined_functions.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a SINGULARITY of the geometry, associated to
 * the geometric approx of the surface. metric \f$=\alpha*Id\f$, \f$\alpha =\f$
 * size.
 *
 */
static int _MMG5_defmetsin(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria         pt;
  MMG5_pPoint        p0;
  MMG5_pPar          par;
  double             *m,n[3],isqhmin,isqhmax,b0[3],b1[3],ps1,tau[3];
  double             ntau2,gammasec[3];
  double             c[3],kappa,maxkappa,alpha,hausd,hausd_v;
  int                ilist,list[_MMG5_LMAX+2],k,i,iel,idp,isloc,init_s;
  unsigned char      i0,i1,i2;

  pt  = &mesh->tria[it];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex: useless for now because new points are created
   * without reference (inside the domain) */
  hausd_v = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc     = 0;
  /* for (i=0; i<mesh->info.npar; i++) { */
  /*   par = &mesh->info.par[i]; */
  /*   if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) { */
  /*     hausd_v = par->hausd; */
  /*     isqhmin = par->hmin; */
  /*     isqhmax = par->hmax; */
  /*     isloc   = 1; */
  /*   } */
  /* } */

  ilist = boulet(mesh,it,ip,list);
  if ( !ilist ) {
    fprintf(stderr,"%s:%d:Error: unable to compute the ball af the point %d.\n",
           __FILE__,__LINE__, idp);
    return(0);
  }

  maxkappa = 0.0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];
    pt  = &mesh->tria[iel];

    /* Computation of the two control points associated to edge p0p1 with
     * p1=mesh->point[pt->v[i1]]: p0 is singular */
    _MMG5_nortri(mesh,pt,n);
    if ( MG_EDG(pt->tag[i2]) )
      _MMG5_bezierEdge(mesh,idp,pt->v[i1],b0,b1,1,n);
    else
      _MMG5_bezierEdge(mesh,idp,pt->v[i1],b0,b1,0,n);

    /* tangent vector */
    tau[0] = 3.0*(b0[0] - p0->c[0]);
    tau[1] = 3.0*(b0[1] - p0->c[1]);
    tau[2] = 3.0*(b0[2] - p0->c[2]);
    ntau2  = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];

    /* 2nd order derivative */
    gammasec[0] = 6.0*p0->c[0] - 12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = 6.0*p0->c[1] - 12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = 6.0*p0->c[2] - 12.0*b0[2] + 6.0*b1[2];
    if ( ntau2 < _MMG5_EPSD )  continue;
    ntau2 = 1.0 / ntau2;

    /* derivative via the normal parametrization */
    ps1  = gammasec[0]*tau[0] + gammasec[1]*tau[1] + gammasec[2]*tau[2];
    c[0] = gammasec[0] - ps1*tau[0]*ntau2;
    c[1] = gammasec[1] - ps1*tau[1]*ntau2;
    c[2] = gammasec[2] - ps1*tau[2]*ntau2;

    kappa = ntau2 * sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

    /* local parameters for triangle */
    hausd   = hausd_v;
    init_s  = 0;
    for (i=0; i<mesh->info.npar; i++) {
      par = &mesh->info.par[i];
      if ( (par->elt == MMG5_Triangle) && (pt->ref == par->ref ) ) {
        if ( !isloc ) {
          hausd   = par->hausd;
          if ( !init_s ) {
            isqhmin = par->hmin;
            isqhmax = par->hmax;
            init_s  = 1;
          }
          else {
            isqhmin = MG_MAX(par->hmin,isqhmin);
            isqhmax = MG_MIN(par->hmax,isqhmax);
          }
        }
        else {
          hausd   = MG_MIN(par->hausd,hausd);
          isqhmin = MG_MAX(par->hmin,isqhmin);
          isqhmax = MG_MIN(par->hmax,isqhmax);
        }
      }
    }
    maxkappa = MG_MAX(kappa/hausd,maxkappa);
  }

  isqhmin  = 1.0 / (isqhmin*isqhmin);
  isqhmax  = 1.0 / (isqhmax*isqhmax);

  alpha = 1.0 / 8.0 * maxkappa;
  alpha = MG_MIN(alpha,isqhmin);
  alpha = MG_MAX(alpha,isqhmax);

  m = &met->m[6*idp];
  memset(m,0,6*sizeof(double));
  m[0] = m[3] = m[5] = alpha;

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
 * \return 1 if success, 0 otherwise.
 *
 * Compute metric tensor associated to a ridge point : convention is a bit weird
 * here :
 * \a p->m[0] is the specific size in direction \a t,
 * \a p->m[1] is the specific size in direction \f$ u_1 = n_1^{}t\f$,
 * \a p->m[2] is the specific size in direction \f$ u_2 = n_2^{}t\f$,
 * \a p->m[3] is the specific size in direction \f$ n_1\f$
 * (computed by the \a _MMG5_intextmet function),
 * \a p->m[4] is the specific size in direction \f$ n_2\f$,
 * (computed by the \a _MMG5_intextmet function),
 * and at each time, metric tensor has to be recomputed, depending on the side.
 *
 */
static int _MMG5_defmetrid(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria     pt;
  MMG5_pPoint    p0,p1,p2;
  _MMG5_Bezier   b;
  MMG5_pPar      par;
  int            k,iel,idp,ilist1,ilist2,ilist,*list,list1[_MMG5_LMAX+2];
  int            list2[_MMG5_LMAX+2],iprid[2],ier,isloc;
  double         *m,isqhmin,isqhmax,*n1,*n2,*n,*t,trot[2],u[2];
  double         r[3][3],lispoi[3*_MMG5_LMAX+1],ux,uy,uz,det,bcu[3];
  double         detg,detd;
  unsigned char  i,i0,i1,i2;

  pt  = &mesh->tria[it];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters */
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc   = 0;
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    if ( /*( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) )
           || */( (par->elt == MMG5_Triangle) && (pt->ref == par->ref ) ) ) {
      if ( !isloc ) {
        isqhmin = par->hmin;
        isqhmax = par->hmax;
        isloc = 1;
      }
      else {
        isqhmin = MG_MAX(isqhmin,par->hmin);
        isqhmax = MG_MIN(isqhmax,par->hmax);
      }
    }
  }

  isqhmin = 1.0 / (isqhmin*isqhmin);
  isqhmax = 1.0 / (isqhmax*isqhmax);

  n1 = &mesh->xpoint[p0->xp].n1[0];
  n2 = &mesh->xpoint[p0->xp].n2[0];
  t  = p0->n;

  m = &met->m[6*idp];
  memset(m,0,6*sizeof(double));
  m[0] = isqhmax;
  m[1] = isqhmax;
  m[2] = isqhmax;
  m[3] = isqhmax;
  m[4] = isqhmax;

  ier = bouletrid(mesh,it,ip,&ilist1,list1,&ilist2,list2,&iprid[0],&iprid[1]);
  if ( !ier ) {
    fprintf(stderr,"%s:%d:Error: unable to compute the two balls at the ridge"
           " point %d.\n",__FILE__,__LINE__, idp);
    return(0);
  }

  /* Specific size in direction of t */
  m[0] = MG_MAX(m[0],_MMG5_ridSizeInTangentDir(mesh,p0,idp,iprid,isqhmin,isqhmax));

  /* Characteristic sizes in directions u1 and u2 */
  for (i=0; i<2; i++) {
    if ( i==0 ) {
      n = n1;
      ilist = ilist1;
      list  = &list1[0];
    }
    else {
      n = n2;
      ilist = ilist2;
      list  = &(list2[0]);
    }
    _MMG5_rotmatrix(n,r);

    /* Apply rotation to the half-ball under consideration */
    i1 = 0;
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

    /* last point : the half-ball is open : ilist tria, and ilist +1 points ;
       lists are enumerated in direct order */
    i2 = _MMG5_inxt2[i1];
    p2 = &mesh->point[pt->v[i2]];

    ux = p2->c[0] - p0->c[0];
    uy = p2->c[1] - p0->c[1];
    uz = p2->c[2] - p0->c[2];

    lispoi[3*ilist+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*ilist+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*ilist+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

    /* At this point, lispoi contains all the points of the half-ball of p0, rotated
       so that t_{p_0}S = [z = 0] */

    /* Rotated tangent vector (trot[2] = 0), and third direction */
    trot[0] = r[0][0]*t[0] + r[0][1]*t[1] + r[0][2]*t[2];
    trot[1] = r[1][0]*t[0] + r[1][1]*t[1] + r[1][2]*t[2];

    u[0] = -trot[1];
    u[1] =  trot[0];

    /* Travel half-ball at p0 and stop at first triangle containing u */
    for (k=0; k<ilist; k++) {
      detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
      detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
      if ( detg > 0.0 && detd > 0.0 )  break;
    }

    /* If triangle not found, try with -u */
    if ( k == ilist ) {
      u[0] *= -1.0;
      u[1] *= -1.0;

      for (k=0; k<ilist; k++) {
        detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
        detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
        if ( detg > 0.0 && detd > 0.0 )  break;
      }
    }
    if ( k == ilist )  continue;

    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];
    pt = &mesh->tria[iel];
    if ( !_MMG5_bezierCP(mesh,pt,&b,1) )  continue;

    /* Barycentric coordinates of vector u in tria iel */
    detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
    detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
    det = detg + detd;
    if ( det < _MMG5_EPSD )  continue;

    det = 1.0 / det;
    bcu[0] = 0.0;
    bcu[1] = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
    bcu[1] *= det;
    assert(bcu[1] <= 1.0);
    bcu[2] = 1.0 - bcu[1];

    /* Computation of tangent vector and second derivative of curve t \mapsto
       b(tbcu) (not in rotated frame) */
    m[i+1] = MG_MAX(m[i+1],
                    _MMG5_ridSizeInNormalDir(mesh,i0,bcu,&b,isqhmin,isqhmax));
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REF vertex of the mesh, associated to the
 * geometric approx of the surface.
 *
 */
static int _MMG5_defmetref(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria         pt;
  MMG5_pPoint        p0,p1;
  _MMG5_Bezier       b;
  MMG5_pPar          par;
  int                i,ilist,list[_MMG5_LMAX+2],k,iel,ipref[2],idp,isloc;
  double             *m,isqhmin,isqhmax,*n,r[3][3],lispoi[3*_MMG5_LMAX+1];
  double             ux,uy,uz,det2d,intm[3],c[3];
  double             tAA[6],tAb[3],hausd;
  unsigned char      i0,i1,i2;

  ipref[0] = ipref[1] = 0;
  pt  = &mesh->tria[it];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex: useless for now because new points are created
   * without reference (inside the domain) */
  hausd   = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc = 0;
  /* for (i=0; i<mesh->info.npar; i++) { */
  /*   par = &mesh->info.par[i]; */
  /*   if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) { */
  /*     hausd   = par->hausd; */
  /*     isqhmin = par->hmin; */
  /*     isqhmax = par->hmax; */
  /*     isloc = 1; */
  /*   } */
  /* } */

  ilist = boulet(mesh,it,ip,list);
  if ( !ilist ) {
    fprintf(stderr,"%s:%d:Error: unable to compute the ball af the point %d.\n",
           __FILE__,__LINE__, idp);
    return(0);
  }

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  n  = &mesh->xpoint[p0->xp].n1[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  _MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];
    pt = &mesh->tria[iel];
    p1 = &mesh->point[pt->v[i1]];

    /* Store the two ending points of ref curves */
    if ( MG_REF & pt->tag[i1] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i2];
      }
      else if ( !ipref[1] && (pt->v[i2] != ipref[0]) ) {
        ipref[1] = pt->v[i2];
      }
      else if ( (pt->v[i2] != ipref[0]) && (pt->v[i2] != ipref[1]) ) {
        fprintf(stderr,"%s:%d:Error: three adjacent ref at a non singular point\n",
               __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }

    if ( MG_REF & pt->tag[i2] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i1];
      }
      else if ( !ipref[1] && (pt->v[i1] != ipref[0]) ) {
        ipref[1] = pt->v[i1];
      }
      else if ( (pt->v[i1] != ipref[0]) && (pt->v[i1] != ipref[1]) ) {
        fprintf(stderr,"%s:%d:Error: three adjacent ref at a non singular point\n",
               __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }

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

  /* Check all projections over tangent plane. */
  for (k=0; k<ilist-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d < 0.0 ) {
      //fprintf(stderr,"PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
  if ( det2d < 0.0 ) {
    //fprintf(stderr,"PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
    return(0);
  }
  assert(ipref[0] && ipref[1]);

  /* At this point, lispoi contains all the points of the ball of p0, rotated
     so that t_{p_0}S = [z = 0], ipref1 and ipref2 are the indices of other ref points. */

  /* Second step : reconstitution of the curvature tensor at p0 in the tangent plane,
     with a quadric fitting approach */
  memset(intm,0x00,3*sizeof(double));
  memset(tAA,0x00,6*sizeof(double));
  memset(tAb,0x00,3*sizeof(double));

  hausd   = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc = 0;
  for (k=0; k<ilist; k++) {
    /* Approximation of the curvature in the normal section associated to tau : by assumption,
       p1 is either regular, either on a ridge (or a singularity), but p0p1 is not ridge*/
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt = &mesh->tria[iel];
    _MMG5_bezierCP(mesh,pt,&b,1);


    /* 1. Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
     * points at edge [i0;i1] extremities and middle.
     * 2. Compute the physical coor \a c of the curve edge's
     * mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i0,b,r,c,lispoi,tAA,tAb);

    /* local parameters */
    for (i=0; i<mesh->info.npar; i++) {
      par = &mesh->info.par[i];
      if ( (par->elt == MMG5_Triangle) && (pt->ref == par->ref ) ) {
        if ( !isloc ) {
          hausd = par->hausd;
          isqhmin = par->hmin;
          isqhmax = par->hmax;
          isloc = 1;
        }
        else {
          // take the wanted value between the two local parameters asked
          // by the user.
          hausd = MG_MIN(hausd,par->hausd);
          isqhmin = MG_MAX(isqhmin,par->hmin);
          isqhmax = MG_MIN(isqhmax,par->hmax);
        }
      }
    }
  }

  isqhmin = 1.0 / (isqhmin*isqhmin);
  isqhmax = 1.0 / (isqhmax*isqhmax);

  return( _MMG5_solveDefmetrefSys(mesh,p0,ipref,r,c,tAA,tAb,m,
                                  isqhmin,isqhmax,hausd) );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param it index of the triangle in which we work.
 * \param ip index of the point on which we want to compute the metric in \a it.
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REGULAR vertex of the mesh, associated to
 * the geometric approx of the surface.
 *
 */
static int _MMG5_defmetreg(MMG5_pMesh mesh,MMG5_pSol met,int it,int ip) {
  MMG5_pTria          pt;
  MMG5_pPoint         p0,p1;
  _MMG5_Bezier        b;
  MMG5_pPar           par;
  int                 ilist,list[_MMG5_LMAX+2],k,iel,idp,isloc,i;
  double              *n,*m,r[3][3],ux,uy,uz,lispoi[3*_MMG5_LMAX+1];
  double              det2d,c[3],isqhmin,isqhmax;
  double              tAA[6],tAb[3],hausd;
  unsigned char       i0,i1;

  pt  = &mesh->tria[it];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex: useless for now because new points are created
   * without reference (inside the domain) */
  hausd   = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc     = 0;
  /* for (i=0; i<mesh->info.npar; i++) { */
  /*   par = &mesh->info.par[i]; */
  /*   if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) { */
  /*     hausd   = par->hausd; */
  /*     isqhmin = par->hmin; */
  /*     isqhmax = par->hmax; */
  /*     isloc   = 1; */
  /*   } */
  /* } */

  ilist = boulet(mesh,it,ip,list);
  if ( !ilist ) {
    fprintf(stderr,"%s:%d:Error: unable to compute the ball af the point %d.\n",
           __FILE__,__LINE__, idp);
    return(0);
  }

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  n  = &p0->n[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  if ( !_MMG5_rotmatrix(n,r) ) return(0);
  m = &met->m[6*idp];

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

  /* Check all projections over tangent plane. */
  for (k=0; k<ilist-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    if ( det2d <= 0.0 ) {
      //fprintf(stderr,"PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilist-1)+1]*lispoi[3*0+2] - lispoi[3*(ilist-1)+2]*lispoi[3*0+1];
  if ( det2d <= 0.0 ) {
    //fprintf(stderr,"PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
    return(0);
  }

  /* At this point, lispoi contains all the points of the ball of p0, rotated
     so that t_{p_0}S = [z = 0] */

  /* Second step : reconstitution of the curvature tensor at p0 in the tangent
     plane, with a quadric fitting approach */
  memset(tAA,0x00,6*sizeof(double));
  memset(tAb,0x00,3*sizeof(double));

  hausd   = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc = 0;
  for (k=0; k<ilist; k++) {
    /* Approximation of the curvature in the normal section associated to tau :
       by assumption, p1 is either regular, either on a ridge (or a
       singularity), but p0p1 is not ridge*/
    iel = list[k] / 3;
    i0  = list[k] % 3;
    i1  = _MMG5_inxt2[i0];
    pt = &mesh->tria[iel];
    _MMG5_bezierCP(mesh,pt,&b,1);

    /* 1. Fill matrice tAA and second member tAb with A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i}) and b=\sum Z_{P_i} with P_i the physical
     * points at edge [i0;i1] extremities and middle.
     * 2. Compute the physical coor \a c of the curve edge's
     * mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i0,b,r,c,lispoi,tAA,tAb);

    /* local parameters */
    for (i=0; i<mesh->info.npar; i++) {
      par = &mesh->info.par[i];
      if ( (par->elt == MMG5_Triangle) && (pt->ref == par->ref ) ) {
        if ( !isloc ) {
          hausd = par->hausd;
          isqhmin = par->hmin;
          isqhmax = par->hmax;
          isloc = 1;
        }
        else {
          // take the minimum value between the two local parameters asked by
          // the user.
          hausd = MG_MIN(hausd,par->hausd);
          isqhmin = MG_MAX(isqhmin,par->hmin);
          isqhmax = MG_MIN(isqhmax,par->hmax);
        }
      }
    }
  }

  isqhmin = 1.0 / (isqhmin*isqhmin);
  isqhmax = 1.0 / (isqhmax*isqhmax);

  /* 2. Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) */
  return(_MMG5_solveDefmetregSys( mesh,r, c, tAA, tAb, m, isqhmin, isqhmax,
                                  hausd));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param np global index of vertex in which we intersect the metrics.
 * \param me physical metric at point \a np.
 * \return 0 if fail, 1 otherwise.
 *
 * Intersect the surface metric held in np (supported in tangent plane of \a np)
 * with 3*3 physical metric in \a me. For ridge points, this function fill the
 * \f$ p_0->m[3]\f$ and \f$ p_0->m[4]\f$ fields that contains respectively the
 * specific sizes in the \f$n_1\f$ and \f$n_2\f$ directions.
 *
 */
static inline
int _MMGS_intextmet(MMG5_pMesh mesh,MMG5_pSol met,int np,double me[6]) {
  MMG5_pPoint         p0;
  double              *n;
  double              dummy_n[3];

   p0 = &mesh->point[np];

   dummy_n[0] = dummy_n[1] = dummy_n[2] = 0.;

   if ( MG_SIN(p0->tag) || (p0->tag & MG_NOM) ) {
     n = &dummy_n[0];
   }
   else {
     n = &p0->n[0];
   }

  return(_MMG5_mmgIntextmet(mesh,met,np,me,n));

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric stucture.
 * \return 0 if fail, 1 otherwise.
 *
 * Define size at points by intersecting the surfacic metric and the
 * physical metric.
 *
 */
int _MMGS_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  double        mm[6];
  int           k;
  char          i,ismet;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining anisotropic map\n");

  if ( mesh->info.hmax < 0.0 ) {
    //  mesh->info.hmax = 0.5 * mesh->info.delta;
    fprintf(stdout,"%s:%d:Error: negative hmax value.\n",__FILE__,__LINE__);
    return(0);
  }

  if ( met->m )
    ismet = 1;
  else {
    ismet = 0;

     _MMG5_calelt     = _MMG5_caltri_ani;
     _MMG5_lenSurfEdg = _MMG5_lenSurfEdg_ani;

     met->np    = mesh->np;
     met->npmax = mesh->npmax;
     met->size  = 6;
     met->dim   = 3;
     _MMG5_ADD_MEM(mesh,(6*(met->npmax+1))*sizeof(double),"solution",return(0));
     _MMG5_SAFE_CALLOC(met->m,6*(mesh->npmax+1),double);
  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->flag = 0;
  }

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->flag || !MG_VOK(ppt) )  continue;
      if ( ismet )  memcpy(mm,&met->m[6*(pt->v[i])],6*sizeof(double));

      if ( MS_SIN(ppt->tag) ) {
        if ( !_MMG5_defmetsin(mesh,met,k,i) )  continue;
      }
      else if ( ppt->tag & MG_GEO ) {
        if ( !_MMG5_defmetrid(mesh,met,k,i))  continue;
      }
      else if ( ppt->tag & MG_REF ) {
        if ( !_MMG5_defmetref(mesh,met,k,i) )  continue;
      }
      else if ( ppt->tag )  continue;
      else {
        if ( !_MMG5_defmetreg(mesh,met,k,i) )  continue;
      }
      if ( ismet ) {
        if ( !_MMGS_intextmet(mesh,met,pt->v[i],mm) ) {
          fprintf(stdout,"%s:%d:Error: unable to intersect metrics"
                  " at point %d.\n",__FILE__,__LINE__, pt->v[i]);
          return(0);
        }
      }
      ppt->flag = 1;
    }
  }

  /* search for unintialized metric */
  _MMG5_defUninitSize(mesh,met,ismet);

  return(1);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 1
 *
 *
 * Enforces mesh gradation by truncating metric field.
 *
 */
int gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria   pt;
  MMG5_pPoint  p1,p2;
  double  *m,mv;
  int     k,it,nup,nu,maxit;
  char    i,ier,i1,i2;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Anisotropic mesh gradation\n");

  mesh->base = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  /* First step : make ridges iso in each apairing direction */
  for (k=1; k<= mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !MG_VOK(p1) ) continue;
    if ( MS_SIN(p1->tag) ) continue;
    if ( !(p1->tag & MG_GEO) ) continue;

    m = &met->m[6*k];
    mv = MG_MAX(m[1],m[2]);
    m[1] = mv;
    m[2] = mv;
    mv = MG_MAX(m[3],m[4]);
    m[3] = mv;
    m[4] = mv;
  }

  /* Second step : standard gradation procedure */
  it = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_iprv2[i];
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];

        if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;
        ier = _MMG5_grad2metSurf(mesh,met,pt,i);
        if ( ier == i1 ) {
          p1->flag = mesh->base;
          nu++;
        }
        else if ( ier == i2 ) {
          p2->flag = mesh->base;
          nu++;
        }
      }
    }
    nup += nu;
  }
  while( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 )  fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
  return(1);
}
