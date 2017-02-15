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
 * \file mmg3d/anisosiz_3d.c
 * \brief Fonctions for anisotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "inlined_functions_3d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param pt pointer toward a tetra.
 * \param m1 computed metric.
 * \return the number of vertices used for the mean computation, 0 if fail.
 *
 * Compute mean metric over the internal tetra \a pt. Do not take into account
 * the metric values at ridges points (because we don't know how to build it).
 *
 */
inline int _MMG5_moymet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt,double *m1) {
  MMG5_pPoint  ppt;
  double  mm[6],*mp;
  double  dd;
  int     i,k,n;

  n = 0;
  for (k=0; k<6; ++k) mm[k] = 0.;
  for(i=0 ; i<4 ; i++) {
    ppt = &mesh->point[pt->v[i]];
    if(!(MG_SIN(ppt->tag) || MG_NOM & ppt->tag) && (ppt->tag & MG_GEO)) continue;
    n++;
    mp = &met->m[6*pt->v[i]];
    for (k=0; k<6; ++k) {
      mm[k] += mp[k];
    }
  }

  if(!n) {
    // fprintf(stderr,"  ## Warning: 4 ridges points... Unable to compute metric.\n");
    return(0);
  }
  dd = 1./n;
  for (k=0; k<6; ++k) m1[k] = mm[k]*dd;
  return(n);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param kel index of the tetra in which we work.
 * \param iface face of the tetra on which we work.
 * \param ip index of the point on which we want to compute the metric
 * (in tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a SINGULARITY of the geometry, associated to
 * the geometric approx of the surface. metric \f$=\alpha*Id\f$, \f$\alpha =\f$
 * size.
 *
 */
static int _MMG5_defmetsin(MMG5_pMesh mesh,MMG5_pSol met,int kel, int iface, int ip) {
  MMG5_pTetra        pt;
  MMG5_pxTetra       pxt;
  MMG5_pPoint        p0;
  MMG5_pPar          par;
  double             *m,n[3],isqhmin,isqhmax,b0[3],b1[3],ps1,tau[3];
  double             ntau2,gammasec[3];
  double             c[3],kappa,maxkappa,alpha, hausd,hausd_v;
  int                lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilist,ilists,ilistv;
  int                k,iel,idp,ifac,isloc,init_s;
  unsigned char      i,i0,i1,i2;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex: useless for now because new points are created
   * without reference (inside the domain) */
  hausd_v = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc   = 0;

  /* if ( mesh->info.parTyp & MG_Vert ) { */
  /*   for (i=0; i<mesh->info.npar; i++) { */
  /*     par = &mesh->info.par[i]; */
  /*     if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) { */
  /*       hausd_v = par->hausd; */
  /*       isqhmin = par->hmin; */
  /*       isqhmax = par->hmax; */
  /*       isloc   = 1; */
  /*       break; */
  /*     } */
  /*   } */
  /*  } */

  if ( mesh->adja[4*(kel-1)+iface+1] ) return(0);
  ilist = _MMG5_boulesurfvolp(mesh,kel,ip,iface,
                              listv,&ilistv,lists,&ilists,(p0->tag & MG_NOM));

  if ( ilist!=1 ) {
    fprintf(stderr,"Error; unable to compute the ball af the point %d.\n", idp);
    fprintf(stderr,"Exit program.\n");
    exit(EXIT_FAILURE);
  }


  /* travel across the ball of ip to find the minimal local params imposed on
   * tetras */
  if ( mesh->info.parTyp & MG_Tetra ) {
    i = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[i];
      if ( par->elt != MMG5_Tetrahedron )  continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd_v = par->hausd;
          isqhmin = par->hmin;
          isqhmax = par->hmax;
          isloc   = 1;
          break;
        }
      }
    } while ( ++i<mesh->info.npar );

    for ( ; i<mesh->info.npar; ++i ) {
      par = &mesh->info.par[i];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd_v = MG_MIN(hausd_v,par->hausd);
          isqhmin = MG_MAX(isqhmin,par->hmin);
          isqhmax = MG_MIN(isqhmax,par->hmax);
          break;
        }
      }
    }
  }


  maxkappa = 0.0;
  init_s = 0;

  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );
    pxt = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i0  = _MMG5_idir[ifac][i];
    i1  = _MMG5_idir[ifac][_MMG5_inxt2[i]];
    i   = _MMG5_iprv2[i];
    i2  = _MMG5_idir[ifac][i];

    /* Computation of the two control points associated to edge p0p1 with
     * p1=mesh->point[pt->v[i1]]: p0 is singular */
    _MMG5_norpts(mesh,pt->v[i0],pt->v[i1],pt->v[i2],n);

    _MMG5_bezierEdge(mesh,idp,pt->v[i1],b0,b1,
                     MG_EDG(pxt->tag[_MMG5_iarf[ifac][i]]),n);

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
    hausd  = hausd_v;

    if ( mesh->info.parTyp & MG_Tria ) {
      if ( !isloc ) {
        for ( i=0; i<mesh->info.npar; ++i ) {
          par = &mesh->info.par[i];

          if ( (par->elt != MMG5_Triangle) || (pxt->ref[ifac] != par->ref ) )
            continue;

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
          break;
        }
      }
      else {
        for ( i=0 ; i<mesh->info.npar; ++i) {
          par = &mesh->info.par[i];

          if ( (par->elt != MMG5_Triangle) || (pxt->ref[ifac] != par->ref ) )
            continue;

          hausd   = MG_MIN(par->hausd,hausd);
          isqhmin = MG_MAX(par->hmin,isqhmin);
          isqhmax = MG_MIN(par->hmax,isqhmax);
          break;
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
 * \param kel index of the tetra in which we work.
 * \param iface face of the tetra on which we work.
 * \param ip index of the point on which we want to compute the metric
 * (in tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Compute metric tensor associated to a ridge point : convention is a bit weird
 * here :
 * \a p->m[0] is the specific size in direction \a t,
 * \a p->m[1] is the specific size in direction \f$ u_1 = n_1^t\f$
 * \a p->m[2] is the specific size in direction \f$ u_2 = n_2^t\f$,
 * and at each time, metric tensor has to be recomputed, depending on the side.
 *
 */
static int _MMG5_defmetrid(MMG5_pMesh mesh,MMG5_pSol met,int kel,
                           int iface, int ip)
{
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_Tria      ptt;
  MMG5_pPoint    p0,p1,p2;
  MMG5_pPar      par;
  _MMG5_Bezier   b;
  int            k,iel,idp,ilist1,ilist2,ilist,*list;
  int            list1[MMG3D_LMAX+2],list2[MMG3D_LMAX+2],iprid[2],ier;
  double         *m,isqhmin,isqhmax,*n1,*n2,*n,*t;
  double         trot[2],u[2],ux,uy,uz,det,bcu[3];
  double         r[3][3],lispoi[3*MMG3D_LMAX+1];
  double         detg,detd;
  int            i,i0,i1,i2,ifac,isloc;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];
  pxt = &mesh->xtetra[pt->xt];

  /* local parameters */
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;

  if ( mesh->info.parTyp ) {
    isloc   = 0;

    i = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[i];
      if ( /*( (par->elt != MMG5_Vertex) || (p0->ref != par->ref ) ) && */
        ( (par->elt != MMG5_Triangle)    || (pxt->ref[iface] != par->ref ) ) &&
        ( (par->elt != MMG5_Tetrahedron) || (pt->ref != par->ref )         ) )
        continue;

      isqhmin = par->hmin;
      isqhmax = par->hmax;
      isloc = 1;
    } while ( ++i<mesh->info.npar );

    for ( ; i<mesh->info.npar; ++i) {
      par = &mesh->info.par[i];

      if ( /*( (par->elt != MMG5_Vertex) || (p0->ref != par->ref ) ) && */
        ( (par->elt != MMG5_Triangle)    || (pxt->ref[iface] != par->ref ) ) &&
        ( (par->elt != MMG5_Tetrahedron) || (pt->ref != par->ref )         ) )
        continue;

      isqhmin = MG_MAX(isqhmin,par->hmin);
      isqhmax = MG_MIN(isqhmax,par->hmax);
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

  // Call bouletrid that construct the surfacic ball
  ier = _MMG5_bouletrid(mesh,kel,iface,ip,&ilist1,list1,&ilist2,list2,
                        &iprid[0],&iprid[1] );
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
      iel  = list[k] / 4;
      ifac = list[k] % 4;
      pt = &mesh->tetra[iel];
      for ( i0=0; i0!=3; ++i0 ) {
        if ( pt->v[_MMG5_idir[ifac][i0]]==idp ) break;
      }
      assert(i0!=3);
      i1   = _MMG5_inxt2[i0];
      p1 = &mesh->point[pt->v[_MMG5_idir[ifac][i1]]];

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
    p2 = &mesh->point[pt->v[_MMG5_idir[ifac][i2]]];

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

    iel  = list[k] / 4;
    ifac = list[k] % 4;
    pt = &mesh->tetra[iel];
    for ( i0=0; i0!=3; ++i0 ) {
      if ( pt->v[_MMG5_idir[ifac][i0]]==idp ) break;
    }
    assert(i0!=3);
    i1  = _MMG5_inxt2[i0];
    i2  = _MMG5_iprv2[i0];

    _MMG5_tet2tri(mesh,iel,ifac,&ptt);
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];
    if ( !_MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,i)) )  continue;

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
 * \param kel index of the triangle in which we work.
 * \param iface face of the tetra on which we work.
 * \param ip index of the point on which we want to compute the metric
 * (in tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REF vertex of the mesh, associated to the
 * geometric approx of the surface.
 *
 */
static int _MMG5_defmetref(MMG5_pMesh mesh,MMG5_pSol met,int kel, int iface, int ip) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_Tria     ptt;
  MMG5_pPoint   p0,p1;
  MMG5_pxPoint  px0;
  _MMG5_Bezier  b;
  MMG5_pPar     par;
  int           lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv,ilist;
  int           k,iel,ipref[2],idp,ifac,isloc;
  double        *m,isqhmin,isqhmax,*n,r[3][3],lispoi[3*MMG3D_LMAX+1];
  double        ux,uy,uz,det2d,c[3];
  double        tAA[6],tAb[3], hausd;
  unsigned char i1,i2,itri1,itri2,i;

  ipref[0] = ipref[1] = 0;
  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex: useless for now because new points are created
   * without reference (inside the domain) */
  hausd   = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc = 0;

  /* if ( mesh->info.parTyp & MG_Vert ) { */
  /* for (i=0; i<mesh->info.npar; i++) { */
  /*   par = &mesh->info.par[i]; */
  /*   if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) { */
  /*     hausd   = par->hausd; */
  /*     isqhmin = par->hmin; */
  /*     isqhmax = par->hmax; */
  /*     isloc = 1; */
  /*   } */
  /* } */
  /* } */

  ilist = _MMG5_boulesurfvolp(mesh,kel,ip,iface,listv,&ilistv,lists,&ilists,0);

  if ( ilist!=1 ) {
    fprintf(stderr,"%s:%d:Error: unable to compute the ball af the point %d.\n",
           __FILE__,__LINE__, idp);
    fprintf(stderr,"Exit program.\n");
    exit(EXIT_FAILURE);
  }

  /* travel across the ball of ip to find the minimal local params imposed on
   * tetras */
  if ( mesh->info.parTyp & MG_Tetra ) {
    i = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[i];
      if ( par->elt != MMG5_Tetrahedron )  continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd   = par->hausd;
          isqhmin = par->hmin;
          isqhmax = par->hmax;
          isloc   = 1;
          break;
        }
      }
    } while ( ++i<mesh->info.npar );

    for ( ; i<mesh->info.npar; ++i ) {
      par = &mesh->info.par[i];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd   = MG_MIN(hausd  ,par->hausd);
          isqhmin = MG_MAX(isqhmin,par->hmin);
          isqhmax = MG_MIN(isqhmax,par->hmax);
          break;
        }
      }
    }
  }

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  assert( p0->xp && !MG_SIN(p0->tag) && MG_EDG(p0->tag) && !(MG_NOM & p0->tag) );

  px0 = &mesh->xpoint[p0->xp];

  n  = &px0->n1[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  _MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );
    pxt   = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0    = _MMG5_idir[ifac][i];
    itri1 = _MMG5_inxt2[i];
    i1    = _MMG5_idir[ifac][itri1];
    itri2 = _MMG5_iprv2[i];
    i2    = _MMG5_idir[ifac][itri2];
    p1    = &mesh->point[pt->v[i1]];

    /* Store the two ending points of ref curves */
    if ( MG_REF & pxt->tag[_MMG5_iarf[ifac][itri1]] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i2];
      }
      else if ( !ipref[1] && (pt->v[i2] != ipref[0]) ) {
        ipref[1] = pt->v[i2];
      }
      else if ( (pt->v[i2] != ipref[0]) && (pt->v[i2] != ipref[1]) ) {
        fprintf(stderr,"%s:%d:Error: three adjacent ref at a non singular point.\n",
               __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }

    if ( MG_REF & pxt->tag[_MMG5_iarf[ifac][itri2]] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i1];
      }
      else if ( !ipref[1] && (pt->v[i1] != ipref[0]) ) {
        ipref[1] = pt->v[i1];
      }
      else if ( (pt->v[i1] != ipref[0]) && (pt->v[i1] != ipref[1]) ) {
        fprintf(stderr,"%s:%d:Error: three adjacent ref at a non singular point.\n",
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
  lispoi[3*ilists+1] =  lispoi[1];
  lispoi[3*ilists+2] =  lispoi[2];
  lispoi[3*ilists+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    assert(det2d);
    if ( det2d <= 0.0 ) {
      //fprintf(stderr,"PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  assert(det2d);
  if ( det2d <= 0.0 ) {
    //fprintf(stderr,"PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
    return(0);
  }
  assert(ipref[0] && ipref[1]);

  /* At this point, lispoi contains all the points of the ball of p0, rotated so
     that t_{p_0}S = [z = 0], ipref1 and ipref2 are the indices of other ref
     points. */

  /* Second step : reconstitution of the curvature tensor at p0 in the tangent
     plane, with a quadric fitting approach */
  memset(tAA,0x00,6*sizeof(double));
  memset(tAb,0x00,3*sizeof(double));

  for (k=0; k<ilists; k++) {
    /* Approximation of the curvature in the normal section associated to tau :
       by assumption, p1 is either regular, either on a ridge (or a
       singularity), but p0p1 is not ridge*/
    iel  = lists[k] / 4;
    ifac = lists[k] % 4;
    pt  = &mesh->tetra[iel];
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0  = _MMG5_idir[ifac][i];
    // i1  = _MMG5_idir[ifac][_MMG5_inxt2[i]];

    _MMG5_tet2tri(mesh,iel,ifac,&ptt);

    _MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,ifac));

    /* 1. Fill matrice tAA and second member tAb with \f$A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i})\f$ and \f$b=\sum Z_{P_i}\f$ with \f$P_i\f$
     * the physical points at edge [i0;i1] extremities and middle.
     *
     * 2. Compute the physical coor \a c of the curve edge's mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i,b,r,c,lispoi,tAA,tAb);

    /* local parameters */
    if ( mesh->info.parTyp & MG_Tria ) {
      if ( isloc ) {
        for ( i=0; i<mesh->info.npar; ++i ) {
          par = &mesh->info.par[i];

          if ( (par->elt != MMG5_Triangle) || (pxt->ref[ifac] != par->ref ) )
            continue;

          hausd   = MG_MIN(par->hausd,hausd);
          isqhmin = MG_MAX(par->hmin,isqhmin);
          isqhmax = MG_MIN(par->hmax,isqhmax);
          break;
        }
      }
      else {
        for ( i=0; i<mesh->info.npar; ++i ) {
          par = &mesh->info.par[i];

          if ( (par->elt != MMG5_Triangle) || (pxt->ref[ifac] != par->ref ) )
            continue;

          hausd   = par->hausd;
          isqhmin = par->hmin;
          isqhmax = par->hmax;
          isloc   = 1;
          break;
        }
      }
    }
  }

  isqhmin = 1.0 / (isqhmin*isqhmin);
  isqhmax = 1.0 / (isqhmax*isqhmax);

  /* Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) */
  return(_MMG5_solveDefmetrefSys( mesh, p0, ipref, r, c, tAA, tAb, m,
                                  isqhmin, isqhmax, hausd));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param kel index of the triangle in which we work.
 * \param iface working face.
 * \param ip index of the point on which we want to compute the metric
 * in (tetra \a kel).
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a REGULAR vertex of the mesh, associated to
 * the geometric approx of the surface.
 *
 */
static int _MMG5_defmetreg(MMG5_pMesh mesh,MMG5_pSol met,int kel,int iface, int ip) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_Tria      ptt;
  MMG5_pPoint    p0,p1;
  MMG5_pxPoint   px0;
  _MMG5_Bezier   b;
  MMG5_pPar      par;
  int            lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv,ilist;
  int            k,iel,idp,ifac,isloc;
  double         *n,*m,r[3][3],ux,uy,uz,lispoi[3*MMG3D_LMAX+1];
  double         det2d,c[3],isqhmin,isqhmax;
  double         tAA[6],tAb[3],hausd;
  unsigned char  i1,i;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex */
  hausd   = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc     = 0;
  /* if ( mesh->info.parTyp & MG_Vert ) { */
  /* for (i=0; i<mesh->info.npar; i++) { */
  /*   par = &mesh->info.par[i]; */
  /*   if ( (par->elt == MMG5_Vertex) && (p0->ref == par->ref ) ) { */
  /*     hausd   = par->hausd; */
  /*     isqhmin = par->hmin; */
  /*     isqhmax = par->hmax; */
  /*     isloc   = 1; */
  /*   } */
  /* } */
  /*  } */

  ilist = _MMG5_boulesurfvolp(mesh,kel,ip,iface,listv,&ilistv,lists,&ilists,0);

  if ( ilist!=1 ) {
    fprintf(stderr,"%s:%d:Error: unable to compute the ball af the point %d.\n",
           __FILE__,__LINE__, idp);
    fprintf(stderr,"Exit program.\n");
    exit(EXIT_FAILURE);
  }

  /* travel across the ball of ip to find the minimal local params imposed on
   * tetras */
  if ( mesh->info.parTyp & MG_Tetra ) {
    i = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[i];
      if ( par->elt != MMG5_Tetrahedron )  continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd   = par->hausd;
          isqhmin = par->hmin;
          isqhmax = par->hmax;
          isloc   = 1;
          break;
        }
      }
    } while ( ++i<mesh->info.npar );

    for ( ; i<mesh->info.npar; ++i ) {
      par = &mesh->info.par[i];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd   = MG_MIN(hausd  ,par->hausd);
          isqhmin = MG_MAX(isqhmin,par->hmin);
          isqhmax = MG_MIN(isqhmax,par->hmax);
          break;
        }
      }
    }
  }

  /* Computation of the rotation matrix T_p0 S -> [z = 0] */
  assert( !(p0->tag & MG_NOSURF) );
  assert( p0->xp && !MG_SIN(p0->tag) && !MG_EDG(p0->tag)
              && !(MG_NOM & p0->tag) );
  px0 = &mesh->xpoint[p0->xp];

  n  = &px0->n1[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > _MMG5_EPSD2 );

  _MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0 = _MMG5_idir[ifac][i];
    i1 = _MMG5_idir[ifac][_MMG5_inxt2[i]];
    p1 = &mesh->point[pt->v[i1]];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
  }

  /* list goes modulo ilist */
  lispoi[3*ilists+1] =  lispoi[1];
  lispoi[3*ilists+2] =  lispoi[2];
  lispoi[3*ilists+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    assert(det2d);
    if ( det2d <= 0.0 ) {
      //fprintf(stderr,"PROBLEM : BAD PROJECTION OVER TANGENT PLANE %f \n", det2d);
      return(0);
    }
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  assert(det2d);
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

  for (k=0; k<ilists; k++) {
    /* Approximation of the curvature in the normal section associated to tau :
       by assumption, p1 is either regular, either on a ridge (or a
       singularity), but p0p1 is not ridge*/
    iel  = lists[k] / 4;
    ifac = lists[k] % 4;
    pt  = &mesh->tetra[iel];
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[_MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0  = _MMG5_idir[ifac][i];
    // i1  = _MMG5_idir[ifac][_MMG5_inxt2[i]];

    _MMG5_tet2tri(mesh,iel,ifac,&ptt);

    _MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,ifac));

    /* 1. Fill matrice tAA and second member tAb with \f$A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i})\f$ and \f$b=\sum Z_{P_i}\f$ with P_i the
     * physical points at edge [i0;i1] extremities and middle.
     *
     * 2. Compute the physical coor \a c of the curve edge's mid-point.
     */
    _MMG5_fillDefmetregSys(k,p0,i,b,r,c,lispoi,tAA,tAb);

    /* local parameters */
    if ( mesh->info.parTyp & MG_Tria ) {
      if ( isloc ) {
        for ( i=0; i<mesh->info.npar; ++i ) {
          par = &mesh->info.par[i];

          if ( (par->elt != MMG5_Triangle) || (pxt->ref[ifac] != par->ref ) )
            continue;

          hausd   = MG_MIN(par->hausd,hausd);
          isqhmin = MG_MAX(par->hmin,isqhmin);
          isqhmax = MG_MIN(par->hmax,isqhmax);
          break;
        }
      }
      else {
        for ( i=0; i<mesh->info.npar; ++i ) {
          par = &mesh->info.par[i];

          if ( (par->elt != MMG5_Triangle) || (pxt->ref[ifac] != par->ref ) )
            continue;

          hausd   = par->hausd;
          isqhmin = par->hmin;
          isqhmax = par->hmax;
          isloc   = 1;
          break;
        }
      }
    }
  }

  isqhmin = 1.0 / (isqhmin*isqhmin);
  isqhmax = 1.0 / (isqhmax*isqhmax);

  /* Solve tAA * tmp_m = tAb and fill m with tmp_m (after rotation) */
  return(_MMG5_solveDefmetregSys( mesh,r, c, tAA, tAb, m, isqhmin, isqhmax,
                                  hausd));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a non-boundary vertex of the mesh.
 * Allocate the metric if needed.
 * Truncate the metric at the hmin/hmax values.
 *
 */
static inline
int _MMG5_defmetvol(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt,ptloc;
  MMG5_pPoint   ppt;
  MMG5_pPar     par;
  double        v[3][3],lambda[3],isqhmax,isqhmin,*m;
  int           list[MMG3D_LMAX+2],ilist,k,l,i,j,isloc,ip;

  isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->flag = 0;
  }

  if ( !met->m ) {
    met->np    = mesh->np;
    met->npmax = mesh->npmax;
    met->size  = 6;
    met->dim   = 3;
    _MMG5_ADD_MEM(mesh,(6*(met->npmax+1))*sizeof(double),"solution",return(0));
    _MMG5_SAFE_MALLOC(met->m,(6*(mesh->npmax+1)),double);


    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;

      for ( l=0; l<4; ++l ) {

        ip = pt->v[l];
        ppt = &mesh->point[ip];
        if ( ppt->flag || ppt->tag & MG_BDY || !MG_VOK(ppt) ) continue;

        /** 1. no metric is provided: Set isotropic hmax size at the vertex */

        /** First step: search for local parameters */
        if ( mesh->info.parTyp ) {
          isqhmax = mesh->info.hmax;
          isloc   = 0;

          /* Local parameters at vertex */
          /* if ( mesh->info.parTyp & MG_Vert ) { */
             /* for (i=0; i<mesh->info.npar; i++) { */
          /*   par = &mesh->info.par[i]; */
          /*   if ( (par->elt == MMG5_Vertex) && (ppt->ref == par->ref ) ) { */
          /*     isqhmax = par->hmax; */
          /*     isloc   = 1; */
          /*     break;       */
          /*   } */
          /* } */
          /* } */

          /* Local parameters at tetra */
          if ( mesh->info.parTyp & MG_Tetra ) {
            ilist = _MMG5_boulevolp(mesh,k,l,list);

            i = 0;
            do
            {
              if ( isloc )  break;

              par = &mesh->info.par[i];
              if ( par->elt != MMG5_Tetrahedron )  continue;

              for ( j=0; j<ilist; ++j ) {
                ptloc = &mesh->tetra[list[j]/4];
                if ( par->ref == ptloc->ref ) {
                  isqhmax = par->hmax;
                  isloc   = 1;
                  break;
                }
              }
            } while ( ++i<mesh->info.npar );

            for ( ; i<mesh->info.npar; ++i ) {
              par = &mesh->info.par[i];
              if ( par->elt != MMG5_Tetrahedron ) continue;

              for ( j=0; j<ilist; ++j ) {
                ptloc = &mesh->tetra[list[j]/4];
                if ( par->ref == ptloc->ref ) {
                  isqhmax = MG_MIN(isqhmax,par->hmax);
                  break;
                }
              }
            }
          }
          isqhmax = 1./(isqhmax*isqhmax);
        }

        /** Second step: set metric */
        m = &met->m[met->size*ip];
        m[0] = isqhmax;
        m[1] = 0.;
        m[2] = 0.;
        m[3] = isqhmax;
        m[4] = 0.;
        m[5] = isqhmax;
        ppt->flag = 1;
      }
    }
    return(1);
  }

  /** 2. A metric is provided: truncate it by hmax/hmin */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for ( l=0; l<4; ++l ) {

      ip = pt->v[l];
      ppt = &mesh->point[ip];
      if ( ppt->flag || ppt->tag & MG_BDY || !MG_VOK(ppt) ) continue;

      /** First step: search for local parameters */
      if ( mesh->info.parTyp ) {
        isqhmin = mesh->info.hmin;
        isqhmax = mesh->info.hmax;
        isloc   = 0;

        /* Local parameters at vertex */
        /* if ( mesh->info.parTyp & MG_Vert ) { */
           /* for (i=0; i<mesh->info.npar; i++) { */
        /*   par = &mesh->info.par[i]; */
        /*   if ( (par->elt == MMG5_Vertex) && (ppt->ref == par->ref ) ) { */
        /*     isqhmin = par->hmin; */
        /*     isqhmax = par->hmax; */
        /*     isloc   = 1; */
        /*     break;       */
        /*   } */
        /* } */
        /* } */

        /* Local parameters at tetra */
        if ( mesh->info.parTyp & MG_Tetra ) {
          ilist = _MMG5_boulevolp(mesh,k,l,list);

          i = 0;
          do
          {
            if ( isloc )  break;

            par = &mesh->info.par[i];
            if ( par->elt != MMG5_Tetrahedron )  continue;

            for ( j=0; j<ilist; ++j ) {
              ptloc = &mesh->tetra[list[j]/4];
              if ( par->ref == ptloc->ref ) {
                isqhmin = par->hmin;
                isqhmax = par->hmax;
                isloc   = 1;
                break;
              }
            }
          } while ( ++i<mesh->info.npar );

          for ( ; i<mesh->info.npar; ++i ) {
            par = &mesh->info.par[i];
            if ( par->elt != MMG5_Tetrahedron ) continue;

            for ( j=0; j<ilist; ++j ) {
              ptloc = &mesh->tetra[list[j]/4];
              if ( par->ref == ptloc->ref ) {
                isqhmin = MG_MAX(isqhmin,par->hmin);
                isqhmax = MG_MIN(isqhmax,par->hmax);
                break;
              }
            }
          }
        }
        isqhmin = 1./(isqhmin*isqhmin);
        isqhmax = 1./(isqhmax*isqhmax);
      }


      /** Second step: set metric */
      m = &met->m[met->size*ip];
      _MMG5_eigenv(1,m,lambda,v);

      for (i=0; i<3; i++) {
        if(lambda[i]<=0) {
          fprintf(stderr,"%s:%d:Error: wrong metric at point %d -- eigenvalues :"
                  " %e %e %e\n",__FILE__,__LINE__,
                  ip,lambda[0],lambda[1],lambda[2]);
          return(0);
        }
        lambda[i]=MG_MIN(isqhmin,lambda[i]);
        lambda[i]=MG_MAX(isqhmax,lambda[i]);
      }

      m[0] = v[0][0]*v[0][0]*lambda[0] + v[1][0]*v[1][0]*lambda[1]
        + v[2][0]*v[2][0]*lambda[2];
      m[1] = v[0][0]*v[0][1]*lambda[0] + v[1][0]*v[1][1]*lambda[1]
        + v[2][0]*v[2][1]*lambda[2];
      m[2] = v[0][0]*v[0][2]*lambda[0] + v[1][0]*v[1][2]*lambda[1]
        + v[2][0]*v[2][2]*lambda[2];
      m[3] = v[0][1]*v[0][1]*lambda[0] + v[1][1]*v[1][1]*lambda[1]
        + v[2][1]*v[2][1]*lambda[2];
      m[4] = v[0][1]*v[0][2]*lambda[0] + v[1][1]*v[1][2]*lambda[1]
        + v[2][1]*v[2][2]*lambda[2];
      m[5] = v[0][2]*v[0][2]*lambda[0] + v[1][2]*v[1][2]*lambda[1]
        + v[2][2]*v[2][2]*lambda[2];
    }
  }

  return(1);
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
int _MMG3D_intextmet(MMG5_pMesh mesh,MMG5_pSol met,int np,double me[6]) {
  MMG5_pPoint         p0;
  MMG5_pxPoint        go;
  double              *n;
  double              dummy_n[3];

   p0 = &mesh->point[np];

   dummy_n[0] = dummy_n[1] = dummy_n[2] = 0.;

   if ( MG_SIN(p0->tag) || (p0->tag & MG_NOM) ) {
     n = &dummy_n[0];
   }
   else if ( p0->tag & MG_GEO ) {
     /* Take the tangent at point */
     n = &p0->n[0];
   }
   else {
     /* Take the normal at point */
     assert(p0->xp);
     go = &mesh->xpoint[p0->xp];
     n = &go->n1[0];
   }

  return(_MMG5_mmgIntextmet(mesh,met,np,me,n));

}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param iel tetra in which we work
 * \param iploc local index of point where metric must be defined
 * \param iface face to which belong \a iploc
 * \return 0 if fail, 1 otherwise.
 *
 * Define metric at point \a iploc of tetra \a iel for \a iploc a surface point
 * marked by option nosurf.
 *
 * If a metric is provided, preserve it except on
 * ridges points where a metric with hmax size and directions normal to the
 * ridges is intersected with the user metric.
 *
 * Without initial metric, an isotropic metric is computed from the mean of the
 * lengths of the edges passing through the point (except on ridges points).
 *
 */
static inline
int _MMG3D_nosurfsiz_ani(MMG5_pMesh mesh,MMG5_pSol met,int iel, int iploc,
                         int iface,int ismet) {

  MMG5_pTetra   pt;
  MMG5_pPoint   ppt,p0,p1;
  double        *m,isqhmin,isqhmax,ux,uy,uz,lm,lambda[3],v[3][3];
  int           lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv;
  int           i,iadr,i0,ip0,ip1,i1,ia,j;

  pt    = &mesh->tetra[iel];
  ip0   = pt->v[iploc];
  iadr  = ip0*6;
  ppt   = &mesh->point[ip0];

  if ( ismet ) {

    if ( ppt->flag ) return(1);

    /** A metric is provided: metric truncature */
    /* Step 1: recover the local hmin and hmax values */
    if ( !(ppt->tag & MG_NOM) ) {
      if ( _MMG5_boulesurfvolp(mesh,iel,iploc,iface,listv,&ilistv,lists,&ilists,0)!=1 ||
           !_MMG3D_localParamReg(mesh,ip0,listv,ilistv,lists,ilists,NULL,&isqhmin,&isqhmax) )
      {
        isqhmin = mesh->info.hmin;
        isqhmax = mesh->info.hmax;
      }
    }
    else {
      i0 = _MMG5_idirinv[iface][iploc];
      j = _MMG5_inxt2[i0];
      ia = _MMG5_iarf[iface][j];

      if ( !_MMG3D_localParamNm(mesh,iel,iface,ia,NULL,&isqhmin,&isqhmax) ) {
        isqhmin  = mesh->info.hmin;
        isqhmax  = mesh->info.hmax;
      }
    }
    isqhmin = 1./(isqhmin*isqhmin);
    isqhmax = 1./(isqhmax*isqhmax);

    /* Step 2: size truncature */
    m = &met->m[iadr];
    _MMG5_eigenv(1,m,lambda,v);

    for (i=0; i<3; i++) {
      if(lambda[i]<=0) {
        fprintf(stderr,"%s:%d:Error: wrong metric at point %d -- eigenvalues :"
                " %e %e %e\n",__FILE__,__LINE__,
                ip0,lambda[0],lambda[1],lambda[2]);
        return(0);
      }
      lambda[i]=MG_MIN(isqhmin,lambda[i]);
      lambda[i]=MG_MAX(isqhmax,lambda[i]);
    }

    m[0] = v[0][0]*v[0][0]*lambda[0] + v[1][0]*v[1][0]*lambda[1]
      + v[2][0]*v[2][0]*lambda[2];
    m[1] = v[0][0]*v[0][1]*lambda[0] + v[1][0]*v[1][1]*lambda[1]
      + v[2][0]*v[2][1]*lambda[2];
    m[2] = v[0][0]*v[0][2]*lambda[0] + v[1][0]*v[1][2]*lambda[1]
      + v[2][0]*v[2][2]*lambda[2];
    m[3] = v[0][1]*v[0][1]*lambda[0] + v[1][1]*v[1][1]*lambda[1]
      + v[2][1]*v[2][1]*lambda[2];
    m[4] = v[0][1]*v[0][2]*lambda[0] + v[1][1]*v[1][2]*lambda[1]
      + v[2][1]*v[2][2]*lambda[2];
    m[5] = v[0][2]*v[0][2]*lambda[0] + v[1][2]*v[1][2]*lambda[1]
      + v[2][2]*v[2][2]*lambda[2];
  }
  else {
    /** No metric is provided: Define size at regular surface point for the
     * -nosurf option (ie a manifold point): the size is computed as the mean of
     * the length of edges passing through the point */
    if ( !(ppt->tag & MG_NOM) ) {

      if ( ppt->flag ) return(1);

      /* First step: search for local parameters */
      // To improve: compute the volume and surface ball af the non-manifold
      // point to apply local parameters
      if ( _MMG5_boulesurfvolp(mesh,iel,iploc,iface,listv,&ilistv,lists,&ilists,0) !=1 )
        return(0);

      if ( !_MMG3D_localParamReg(mesh,ip0,listv,ilistv,lists,ilists,
                                 NULL,&isqhmin,&isqhmax) ) {
        isqhmin = mesh->info.hmin;
        isqhmax = mesh->info.hmax;
      }

      /* Step 2: edge length computation */
      lm = _MMG5_meansizreg_iso(mesh,met,ip0,lists,ilists,isqhmin,isqhmax);
      met->m[iadr] = 1./(lm*lm);
      met->m[iadr+3] = met->m[iadr+5] = met->m[iadr];
      met->m[iadr+1] = met->m[iadr+2] = met->m[iadr+4] = 0;
    }
    else {
      i0 = _MMG5_idirinv[iface][iploc];
      j = _MMG5_inxt2[i0];
      ia = _MMG5_iarf[iface][j];
      i0 = _MMG5_iare[ia][0];
      i1 = _MMG5_iare[ia][1];
      ip0 = pt->v[i0];
      ip1 = pt->v[i1];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( !MG_EDG(mesh->xtetra[pt->xt].tag[ia]) ) return(0);

      /** First step: search for local parameters */
      if ( !_MMG3D_localParamNm(mesh,iel,iface,ia,NULL,&isqhmin,&isqhmax) ) {
        isqhmin  = mesh->info.hmin;
        isqhmax  = mesh->info.hmax;
      }
      isqhmin =  1.0 / (isqhmin*isqhmin);
      isqhmax =  1.0 / (isqhmax*isqhmax);

      /** Second step: Very rough eval of the metric at non-manifold point: take
       * the non-manifold edge length */
      ux = (p0->c[0]-p1->c[0]);
      uy = (p0->c[1]-p1->c[1]);
      uz = (p0->c[2]-p1->c[2]);
      lm = ux*ux + uy*uy + uz*uz;

      lm = MG_MAX(isqhmax,MG_MIN(isqhmin,1./lm));
      met->m[iadr] = MG_MAX(met->m[iadr],lm);

      met->m[iadr+3] = met->m[iadr+5] = met->m[iadr];
      met->m[iadr+1] = met->m[iadr+2] = met->m[iadr+4] = 0;
    }
  }

  return(1);

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
int _MMG3D_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  double        mm[6];
  int           k,l,iploc;
  char          i,ismet;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining anisotropic map\n");

  if ( mesh->info.hmax < 0.0 ) {
    //  mesh->info.hmax = 0.5 * mesh->info.delta;
    fprintf(stderr,"%s:%d:Error: negative hmax value.\n",__FILE__,__LINE__);
    return(0);
  }

  if ( met->m )
    ismet = 1;
  else {
    ismet = 0;

    _MMG5_caltet         = _MMG5_caltet_ani;
    _MMG5_caltri         = _MMG5_caltri_ani;
    _MMG5_lenedg         = _MMG5_lenedg_ani;
    MMG3D_lenedgCoor     = _MMG5_lenedgCoor_ani;
    _MMG5_lenSurfEdg     = _MMG5_lenSurfEdg_ani;
  }

  if ( !_MMG5_defmetvol(mesh,met) )  return(0);

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->flag = 0;
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    // Warning: why are we skipped the tetra with negative refs ?
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    else if ( !pt->xt )  continue;

    pxt = &mesh->xtetra[pt->xt];
    for (l=0; l<4; l++) {
      if ( !(pxt->ftag[l] & MG_BDY) ) continue;
      // In multidomain case, acces the face through a tetra for which it is
      // well oriented.
      if ( !(MG_GET(pxt->ori,l)) ) continue;

      for (i=0; i<3; i++) {
        iploc = _MMG5_idir[l][i];
        ppt   = &mesh->point[pt->v[iploc]];

        if ( !MG_VOK(ppt) )  continue;

        if ( mesh->info.nosurf ) {
          if ( !_MMG3D_nosurfsiz_ani(mesh,met,k,iploc,l,ismet) ) continue;
        }
        else {
          if ( ppt->flag ) continue;

          if ( ismet )  memcpy(mm,&met->m[6*(pt->v[iploc])],6*sizeof(double));

          if ( (MG_SIN(ppt->tag) || (ppt->tag & MG_NOM) ) ) {
            if ( !_MMG5_defmetsin(mesh,met,k,l,iploc) )  continue;
          }
          else if ( ppt->tag & MG_GEO ) {
            if ( !_MMG5_defmetrid(mesh,met,k,l,iploc))  continue;
          }
          else if ( ppt->tag & MG_REF ) {
            if ( !_MMG5_defmetref(mesh,met,k,l,iploc) )  continue;
          } else {
            if ( !_MMG5_defmetreg(mesh,met,k,l,iploc) )  continue;
          }
          if ( ismet ) {
            if ( !_MMG3D_intextmet(mesh,met,pt->v[iploc],mm) ) {
              fprintf(stderr,"%s:%d:Error: unable to intersect metrics"
                      " at point %d.\n",__FILE__,__LINE__, pt->v[iploc]);
              return(0);
            }
          }
        }
        ppt->flag = 1;
      }
    }
  }

  /* search for unintialized metric */
  _MMG5_defUninitSize(mesh,met,ismet);

  return(1);
}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param pt pointer toward a tetra.
 * \param ia edge index in tetra \a pt.
 * \return -1 if no gradation is needed, else index of graded point.
 *
 * Enforces gradation of metric in one extremity of edge \a ia in tetra \a pt
 * with respect to the other.
 *
 */
static inline
int _MMG5_grad2metVol(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt,int ia) {
  MMG5_pPoint    p1,p2;
  double         *mm1,*mm2,m1[6],m2[6],ps1,ps2,ux,uy,uz;
  double         c[5],l,val,t[3];
  double         lambda[3],vp[3][3],alpha,beta,mu[3];
  int            ip1,ip2,kmin,i;
  char           i1,i2,ichg;

  i1  = _MMG5_iare[ia][0];
  i2  = _MMG5_iare[ia][1];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  p1  = &mesh->point[ip1];
  p2  = &mesh->point[ip2];

  mm1  = &met->m[6*ip1];
  mm2  = &met->m[6*ip2];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  uz = p2->c[2] - p1->c[2];

  if ( (!( MG_SIN(p1->tag) || (p1->tag & MG_NOM) )) &&  p1->tag & MG_GEO ) {
    /* Recover normal and metric associated to p1 */
    if( !_MMG5_buildridmet(mesh,met,ip1,ux,uy,uz,m1) )
      return(-1);
  }
  else
    memcpy(m1,mm1,6*sizeof(double));

  if ( (!( MG_SIN(p2->tag) || (p2->tag & MG_NOM) )) && p2->tag & MG_GEO ) {
    /* Recover normal and metric associated to p2 */
    if( !_MMG5_buildridmet(mesh,met,ip2,ux,uy,uz,m2) )
      return(-1);
  }
  else
    memcpy(m2,mm2,6*sizeof(double));

  l = sqrt(ux*ux+uy*uy+uz*uz);

  t[0] = ux/l;
  t[1] = uy/l;
  t[2] = uz/l;

  // edge length in metric m1: sqrt(t^u * m1 * u).
  ps1 =  m1[0]*t[0]*t[0] + 2.0*m1[1]*t[0]*t[1] + m1[3]*t[1]*t[1]
    + 2.0*m1[2]*t[0]*t[2] + 2.0*m1[4]*t[1]*t[2] + m1[5]*t[2]*t[2];
  ps1 = sqrt(ps1);

  // edge length in metric m2: sqrt(t^u * m2 * u).
  ps2 =  m2[0]*t[0]*t[0] + 2.0*m2[1]*t[0]*t[1] + m2[3]*t[1]*t[1]
    + 2.0*m2[2]*t[0]*t[2] + 2.0*m2[4]*t[1]*t[2] + m2[5]*t[2]*t[2];
  ps2 = sqrt(ps2);

  /* Metric in p1 has to be changed */
  if ( ps2 > ps1 ){
    /* compute alpha = h2 + hgrad*l */
    alpha = ps2 /(1.0+mesh->info.hgrad*l*ps2);
    if( ps1 >= alpha -_MMG5_EPS )
      return(-1);

    _MMG5_eigenv(1,m1,lambda,vp);
    /* Project the vector t1 along the main directions of the metric */
    c[0] = t[0]*vp[0][0] + t[1]*vp[0][1] + t[2]*vp[0][2];
    c[1] = t[0]*vp[1][0] + t[1]*vp[1][1] + t[2]*vp[1][2];
    c[2] = t[0]*vp[2][0] + t[1]*vp[2][1] + t[2]*vp[2][2];

    /* Find index of the maximum value of c: this allow to detect which of the
     * main directions of the metric is closest to our edge direction. We want
     * that our new metric respect the gradation related to the size associated
     * to this main direction (the ichg direction). */
    ichg = 0;
    val  = fabs(c[ichg]);
    for (i = 1; i<3; ++i) {
      if ( fabs(c[i]) > val ) {
        val = fabs(c[i]);
        ichg = i;
      }
    }
    assert(c[ichg]*c[ichg] > _MMG5_EPS );
   /* Compute beta coef such as lambda_1 = beta*lambda_1 => h1 = h2 + hgrad*l
    * (see p317 of Charles Dapogny Thesis). */
    beta = (alpha*alpha - ps1*ps1)/(c[ichg]*c[ichg]);

    /* Metric update */
    if( MG_SIN(p1->tag) || (p1->tag & MG_NOM) ){
      /* lambda_new = 0.5 lambda_1 + 0.5 beta lambda_1: here we choose to not
       * respect the gradation in order to restric the influence of the singular
       * points. */
      mm1[0] += 0.5*beta;
      mm1[3] += 0.5*beta;
      mm1[5] += 0.5*beta;
    }
    else if( p1->tag & MG_GEO ) {
      /* lambda[ichg] is the metric eigenvalue associated to the main metric
       * direction closest to our edge direction. Find were is stored this
       * eigenvalue in our special storage of ridge metric (mm-lambda = 0) and
       * update it. */
      c[0] = fabs(mm1[0]-lambda[ichg]);
      c[1] = fabs(mm1[1]-lambda[ichg]);
      c[2] = fabs(mm1[2]-lambda[ichg]);
      c[3] = fabs(mm1[3]-lambda[ichg]);
      c[4] = fabs(mm1[4]-lambda[ichg]);

      // Find index af the minimum value of c
      kmin = 0;
      val = fabs(c[kmin]);
      for (i = 1; i<5; ++i) {
        if ( fabs(c[i]) < val ) {
          val = fabs(c[i]);
          kmin = i;
        }
      }
      mm1[kmin] += beta;
    }
    else {
      /* Update the metric eigenvalue associated to the main metric direction
       * which is closest to our edge direction (because this is the one that is
       * the more influent on our edge length). */
      mu[0] = lambda[0];
      mu[1] = lambda[1];
      mu[2] = lambda[2];

      mu[ichg] += beta;

      m1[0] = mu[0]*vp[0][0]*vp[0][0] + mu[1]*vp[1][0]*vp[1][0] + vp[2][0]*vp[2][0]*mu[2];
      m1[1] = mu[0]*vp[0][0]*vp[0][1] + mu[1]*vp[1][0]*vp[1][1] + vp[2][1]*vp[2][0]*mu[2];
      m1[2] = mu[0]*vp[0][0]*vp[0][2] + mu[1]*vp[1][0]*vp[1][2] + vp[2][0]*vp[2][2]*mu[2];
      m1[3] = mu[0]*vp[0][1]*vp[0][1] + mu[1]*vp[1][1]*vp[1][1] + vp[2][1]*vp[2][1]*mu[2];
      m1[4] = mu[0]*vp[0][1]*vp[0][2] + mu[1]*vp[1][1]*vp[1][2] + vp[2][1]*vp[2][2]*mu[2];
      m1[5] = mu[0]*vp[0][2]*vp[0][2] + mu[1]*vp[1][2]*vp[1][2] + vp[2][2]*vp[2][2]*mu[2];

      memcpy(mm1,m1,6*sizeof(double));
    }
    return(i1);
  }
  /* Metric in p2 has to be changed */
  else{
    alpha = ps1 /(1.0+mesh->info.hgrad*l*ps1);
    if( ps2 >= alpha - _MMG5_EPS)
      return(-1);

    _MMG5_eigenv(1,m2,lambda,vp);

    c[0] = t[0]*vp[0][0] + t[1]*vp[0][1] + t[2]*vp[0][2];
    c[1] = t[0]*vp[1][0] + t[1]*vp[1][1] + t[2]*vp[1][2];
    c[2] = t[0]*vp[2][0] + t[1]*vp[2][1] + t[2]*vp[2][2];

    /* Detect which of the main directions of the metric is closest to our edge
     * direction. */
    ichg = 0;
    val  = fabs(c[ichg]);
    for (i = 1; i<3; ++i) {
      if ( fabs(c[i]) > val ) {
        val = fabs(c[i]);
        ichg = i;
      }
    }
    assert(c[ichg]*c[ichg] > _MMG5_EPS );
    /* Compute beta coef such as lambda_1 = beta*lambda_1 => h1 = h2 + hgrad*l
     * (see p317 of Charles Dapogny Thesis). */
    beta = (alpha*alpha - ps2*ps2)/(c[ichg]*c[ichg]);

    /* Metric update: update the metric eigenvalue associated to the main metric
     * direction which is closest to our edge direction (because this is the
     * one that is the more influent on our edge length). */
    if( MG_SIN(p2->tag) || (p2->tag & MG_NOM) ){
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
      c[3] = fabs(mm2[3]-lambda[ichg]);
      c[4] = fabs(mm2[4]-lambda[ichg]);

      kmin = 0;
      val = fabs(c[kmin]);
      for (i = 1; i<5; ++i) {
        if ( c[i] < fabs(val) ) {
          val = fabs(c[i]);
          kmin = i;
        }
      }
      mm2[kmin] += beta;
    }
    else{
      mu[0] = lambda[0];
      mu[1] = lambda[1];
      mu[2] = lambda[2];

      mu[ichg] += beta;

      m2[0] = mu[0]*vp[0][0]*vp[0][0] + mu[1]*vp[1][0]*vp[1][0] + vp[2][0]*vp[2][0]*mu[2];
      m2[1] = mu[0]*vp[0][0]*vp[0][1] + mu[1]*vp[1][0]*vp[1][1] + vp[2][1]*vp[2][0]*mu[2];
      m2[2] = mu[0]*vp[0][0]*vp[0][2] + mu[1]*vp[1][0]*vp[1][2] + vp[2][0]*vp[2][2]*mu[2];
      m2[3] = mu[0]*vp[0][1]*vp[0][1] + mu[1]*vp[1][1]*vp[1][1] + vp[2][1]*vp[2][1]*mu[2];
      m2[4] = mu[0]*vp[0][1]*vp[0][2] + mu[1]*vp[1][1]*vp[1][2] + vp[2][1]*vp[2][2]*mu[2];
      m2[5] = mu[0]*vp[0][2]*vp[0][2] + mu[1]*vp[1][2]*vp[1][2] + vp[2][2]*vp[2][2]*mu[2];
      memcpy(mm2,m2,6*sizeof(double));
    }
    return(i2);
  }
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
int _MMG5_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_Tria     ptt;
  MMG5_pPoint   p0,p1;
  double        *m,mv;
  int           k,it,itv,nup,nu,nupv,maxit;
  int           i,j,ip0,ip1;
  char          ier,i0,i1;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Anisotropic mesh gradation\n");

  mesh->base = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  /* First step : make ridges iso in each apairing direction */
  for (k=1; k<= mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !MG_VOK(p1) ) continue;
    if ( MG_SIN(p1->tag) || (p1->tag & MG_NOM) ) continue;
    if ( !(p1->tag & MG_GEO) ) continue;

    m = &met->m[6*k];
    mv = MG_MAX(m[1],m[2]);
    m[1] = mv;
    m[2] = mv;
    mv = MG_MAX(m[3],m[4]);
    m[3] = mv;
    m[4] = mv;
  }

  it = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

      if ( pxt ) {
        for (i=0; i<4; i++) {
          if ( pxt->ftag[i] & MG_BDY) {
            /* Gradation along a surface edge */
            /* virtual triangle */
            _MMG5_tet2tri(mesh,k,i,&ptt);
            for (j=0; j<3; j++) {
              i0  = _MMG5_inxt2[j];
              i1  = _MMG5_iprv2[j];
              ip0 = ptt.v[i0];
              ip1 = ptt.v[i1];
              p0  = &mesh->point[ip0];
              p1  = &mesh->point[ip1];
              if ( (p0->flag < mesh->base-1) && (p1->flag < mesh->base-1) )
                continue;
              /* gradation along the tangent plane */
              ier = _MMG5_grad2metSurf(mesh,met,&ptt,j);
              if ( ier == i0 ) {
                p0->flag = mesh->base;
                nu++;
              }
              else if ( ier == i1 ) {
                p1->flag = mesh->base;
                nu++;
              }
            }
          }

          else continue;
        }
      }
    }
    nup += nu;
  }
  while( ++it < maxit && nu > 0 );

  mesh->base = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  nupv = itv = 0;
  maxit = 500;

  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      for (i=0; i<4; i++) {
        /* Gradation along a volume edge */
        i0  = _MMG5_iare[i][0];
        i1  = _MMG5_iare[i][1];
        ip0 = pt->v[i0];
        ip1 = pt->v[i1];
        p0  = &mesh->point[ip0];
        p1  = &mesh->point[ip1];
        if ( p0->flag < mesh->base-1 && p1->flag < mesh->base-1 )  continue;

        ier = _MMG5_grad2metVol(mesh,met,pt,i);
        if ( ier == i0 ) {
          p0->flag = mesh->base;
          nu++;
        }
        else if ( ier == i1 ) {
          p1->flag = mesh->base;
          nu++;
        }
      }
    }
    nupv += nu;
  }
  while( ++itv < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 3 ) {
    if ( abs(mesh->info.imprim) < 5 && !mesh->info.ddebug ) {
      fprintf(stdout,"    gradation: %7d updated, %d iter\n",nup+nupv,it+itv);
    }
    else {
      fprintf(stdout,"    surface gradation: %7d updated, %d iter\n"
              "    volume gradation:  %7d updated, %d iter\n",nup,it,nupv,itv);
    }
  }
  return(1);
}
