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

#include "libmmg3d.h"
#include "inlined_functions_3d_private.h"
#include "mmg3dexterns_private.h"
#include "mmgexterns_private.h"

/**
 * \param dm matrix eigenvalues (1x3 array).
 * \param vp eigenvectors matrix (3x3 array, eigenvectors stored by lines).
 * \return 1 if success, 0 if fail.
 *
 * Print eigendecomposition.
 */
int MMG3D_printEigenv(double dm[3],double vp[3][3]) {

  printf("--- Eigenvalues:\n");
  printf("%e %e %e\n",dm[0],dm[1],dm[2]);
  printf("---Eigenvectors (visualization by columns):\n");
  printf("%e %e %e\n",vp[0][0],vp[1][0],vp[2][0]);
  printf("%e %e %e\n",vp[0][1],vp[1][1],vp[2][1]);
  printf("%e %e %e\n",vp[0][2],vp[1][2],vp[2][2]);

  return 1;
}

/**
 * \param symmat flag for symmetric(1) or non-symmetric(0) matrix..
 * \param m matrix (1x6 or 1x9 array).
 * \return 1 if success, 0 if fail.
 *
 * Print matrix entries.
 */
int MMG3D_printMat(int8_t symmat,double *m) {

  if( symmat ) {
    printf("%e %e %e\n",m[0],m[1],m[2]);
    printf("%e %e %e\n",m[1],m[3],m[4]);
    printf("%e %e %e\n",m[2],m[4],m[5]);
  } else {
    printf("%e %e %e\n",m[0],m[1],m[2]);
    printf("%e %e %e\n",m[3],m[4],m[5]);
    printf("%e %e %e\n",m[6],m[7],m[8]);
  }

  return 1;
}

/**
 * \param symmat flag for symmetric(1) or non-symmetric(0) matrix..
 * \param m first matrix (1x6 or 1x9 array).
 * \param mr second matrix (1x6 or 1x9 array).
 * \return 1 if success, 0 if fail.
 *
 * Print relative error between two matrices, for each matrix entry.
 */
int MMG3D_printErrorMat(int8_t symmat,double *m,double *mr) {
  double dm[9],dd;
  int i,dim;

  if( symmat )
    dim = 6;
  else
    dim = 9;

  dd = 0.0;
  for( i = 0; i < dim; i++ )
    if( fabs(m[i]) > dd )
      dd = fabs(m[i]);
  dd = 1.0 / dd;

  for( i = 0; i < dim; i++ )
      dm[i] = (m[i]-mr[i])*dd;

  if( !MMG3D_printMat(symmat,dm) ) return 0;

  return 1;
}

int MMG3D_chk4ridVertices(MMG5_pMesh mesh, MMG5_pTetra pt) {
  MMG5_pPoint  ppt;
  int          i;
  int          n;

  n = 0;
  for(i=0 ; i<4 ; i++) {
    ppt = &mesh->point[pt->v[i]];
    if ( MG_RID(ppt->tag) ) continue;
    n++;
  }

  if(!n) {
    //fprintf(stderr,"\n  ## Warning: 4 ridges points... Unable to compute metric.\n");
    return 0;
  }

  return 1;
}

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
inline int MMG5_moymet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt,double *m1) {
  MMG5_pPoint   ppt;
  double        mm[6],*mp;
  double        dd;
  int           i,k,n;
  int8_t        ddebug = 0;
  static int8_t mmgWarn=0;

  n = 0;
  for (k=0; k<6; ++k) mm[k] = 0.;
  for(i=0 ; i<4 ; i++) {
    ppt = &mesh->point[pt->v[i]];
    if ( MG_RID(ppt->tag) ) continue;
    n++;
    mp = &met->m[6*pt->v[i]];
    for (k=0; k<6; ++k) {
      mm[k] += mp[k];
    }
  }

  if(!n) {
    if ( ddebug && !mmgWarn ) {
      mmgWarn=1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 tetra with 4 ridges vertices"
              "... Unable to compute metric.\n",__func__);
    }
    return 0;
  }
  dd = 1./n;
  for (k=0; k<6; ++k) m1[k] = mm[k]*dd;
  return n;
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
 * Define metric map at a SINGULARITY (corner, required or non-manifold points)
 * of the geometry, associated to the geometric approx of the surface. metric
 * \f$=\alpha*Id\f$, \f$\alpha =\f$ size.
 *
 */
static int MMG5_defmetsin(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int kel, int iface, int ip) {
  MMG5_pTetra        pt;
  MMG5_pxTetra       pxt;
  MMG5_pPoint        p0;
  MMG5_pPar          par;
  double             *m,n[3],isqhmin,isqhmax,b0[3],b1[3],ps1,tau[3];
  double             ntau2,gammasec[3];
  double             c[3],kappa,maxkappa,alpha, hausd,hausd_v;
  MMG5_int           lists[MMG3D_LMAX+2];
  int64_t            listv[MMG3D_LMAX+2];
  int                k,ilist,ifac,isloc,init_s,ilists,ilistv;
  MMG5_int           idp,iel;
  uint8_t            i,i0,i1,i2;
  static int8_t      mmgWarn = 0;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex: useless for now because new points are created
   * without reference (inside the domain) */
  hausd_v = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc   = 0;

  if ( mesh->adja[4*(kel-1)+iface+1] ) return 0;
  ilist = MMG5_boulesurfvolp(mesh,kel,ip,iface,
                              listv,&ilistv,lists,&ilists,(p0->tag & MG_NOM));

  if ( ilist!=1 ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: at least 1 metric not computed:"
              " unable to compute the ball of point\n",__func__);
      mmgWarn = 1;
    }
    return 0;
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
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i0  = MMG5_idir[ifac][i];
    i1  = MMG5_idir[ifac][MMG5_inxt2[i]];
    i   = MMG5_iprv2[i];
    i2  = MMG5_idir[ifac][i];

    /* Computation of the two control points associated to edge p0p1 with
     * p1=mesh->point[pt->v[i1]]: p0 is singular */
    MMG5_norpts(mesh,pt->v[i0],pt->v[i1],pt->v[i2],n);

    MMG5_BezierEdge(mesh,idp,pt->v[i1],b0,b1,MG_EDG_OR_NOM(pxt->tag[MMG5_iarf[ifac][i]]),n);

    /* tangent vector */
    tau[0] = 3.0*(b0[0] - p0->c[0]);
    tau[1] = 3.0*(b0[1] - p0->c[1]);
    tau[2] = 3.0*(b0[2] - p0->c[2]);
    ntau2  = tau[0]*tau[0] + tau[1]*tau[1] + tau[2]*tau[2];

    /* 2nd order derivative */
    gammasec[0] = 6.0*p0->c[0] - 12.0*b0[0] + 6.0*b1[0];
    gammasec[1] = 6.0*p0->c[1] - 12.0*b0[1] + 6.0*b1[1];
    gammasec[2] = 6.0*p0->c[2] - 12.0*b0[2] + 6.0*b1[2];
    if ( ntau2 < MMG5_EPSD )  continue;
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

  return 1;
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
static int MMG5_defmetrid(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int kel,
                           int iface, MMG5_int ip)
{
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_Tria      ptt;
  MMG5_pPoint    p0,p1,p2;
  MMG5_pPar      par;
  MMG5_Bezier    b;
  MMG5_int       iel,idp,*list;
  int            k,ilist1,ilist2,ilist;
  MMG5_int       list1[MMG3D_LMAX+2],list2[MMG3D_LMAX+2],iprid[2],ier;
  double         *m,isqhmin,isqhmax,*n1,*n2,*n,*t;
  double         trot[2],u[2],ux,uy,uz,det,bcu[3];
  double         r[3][3],lispoi[3*MMG3D_LMAX+1];
  double         detg,detd;
  int            i,i0,i1,i2,ifac,isloc;
  static int8_t  mmgWarn = 0;

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
      if ( ( (par->elt != MMG5_Triangle)    || (pxt->ref[iface] != par->ref ) ) &&
           ( (par->elt != MMG5_Tetrahedron) || (pt->ref != par->ref )         ) )
        continue;

      isqhmin = par->hmin;
      isqhmax = par->hmax;
      isloc = 1;
    } while ( ++i<mesh->info.npar );

    for ( ; i<mesh->info.npar; ++i) {
      par = &mesh->info.par[i];

      if ( ( (par->elt != MMG5_Triangle)    || (pxt->ref[iface] != par->ref ) ) &&
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
  ier = MMG5_bouletrid(mesh,kel,iface,ip,&ilist1,list1,&ilist2,list2,
                        &iprid[0],&iprid[1] );
  if ( !ier ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: at least 1 metric not computed:"
              " unable to compute the ball of point\n",__func__);
      mmgWarn = 1;
    }
  }

  /* Specific size in direction of t */
  m[0] = MG_MAX(m[0],MMG5_ridSizeInTangentDir(mesh,p0,idp,iprid,isqhmin,isqhmax));

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
    MMG5_rotmatrix(n,r);

    /* Apply rotation to the half-ball under consideration */
    i1   = 0;
    ifac = -1; // Remove uninitialized warning
    for (k=0; k<ilist; k++) {
      iel  = list[k] / 4;
      ifac = list[k] % 4;
      pt = &mesh->tetra[iel];
      for ( i0=0; i0!=3; ++i0 ) {
        if ( pt->v[MMG5_idir[ifac][i0]]==idp ) break;
      }
      assert(i0!=3);
      i1   = MMG5_inxt2[i0];
      p1 = &mesh->point[pt->v[MMG5_idir[ifac][i1]]];

      ux = p1->c[0] - p0->c[0];
      uy = p1->c[1] - p0->c[1];
      uz = p1->c[2] - p0->c[2];

      lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
      lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
      lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
    }

    /* last point : the half-ball is open : ilist tria, and ilist +1 points ;
       lists are enumerated in direct order */
    i2 = MMG5_inxt2[i1];
    p2 = &mesh->point[pt->v[MMG5_idir[ifac][i2]]];

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
      if ( pt->v[MMG5_idir[ifac][i0]]==idp ) break;
    }
    assert(i0!=3);

    assert( 0<=ifac && ifac<4 && "unexpected local face idx");
    MMG5_tet2tri(mesh,iel,ifac,&ptt);
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];
    if ( !MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,ifac)) )  continue;

    /* Barycentric coordinates of vector u in tria iel */
    detg = lispoi[3*k+1]*u[1] - lispoi[3*k+2]*u[0];
    detd = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
    det = detg + detd;
    if ( det < MMG5_EPSD )  continue;

    det = 1.0 / det;
    bcu[0] = 0.0;
    bcu[1] = u[0]*lispoi[3*(k+1)+2] - u[1]*lispoi[3*(k+1)+1];
    bcu[1] *= det;
    assert(bcu[1] <= 1.0);
    bcu[2] = 1.0 - bcu[1];

    /* Computation of tangent vector and second derivative of curve t \mapsto
       b(tbcu) (not in rotated frame) */
    m[i+1] = MG_MAX(m[i+1],
                    MMG5_ridSizeInNormalDir(mesh,i0,bcu,&b,isqhmin,isqhmax));
  }

  return 1;
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
static int MMG5_defmetref(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int kel, int iface, int ip) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_Tria     ptt;
  MMG5_pPoint   p0,p1;
  MMG5_pxPoint  px0;
  MMG5_Bezier   b;
  MMG5_pPar     par;
  MMG5_int      lists[MMG3D_LMAX+2];
  int64_t       listv[MMG3D_LMAX+2];
  int           k,ilists,ilistv,ilist;
  MMG5_int      iel,ipref[2],idp;
  int           ifac,isloc;
  double        *m,isqhmin,isqhmax,*n,r[3][3],lispoi[3*MMG3D_LMAX+1];
  double        ux,uy,uz,det2d,c[3];
  double        tAA[6],tAb[3], hausd;
  uint8_t       i1,i2,itri1,itri2,i;
  static int8_t mmgWarn0=0,mmgWarn1=0;

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

  ilist = MMG5_boulesurfvolp(mesh,kel,ip,iface,listv,&ilistv,lists,&ilists,0);

  if ( ilist!=1 ) {
    if ( !mmgWarn0 ) {
      fprintf(stderr,"\n  ## Warning: %s: at least 1 metric not computed:"
              " unable to compute the ball of point\n",__func__);
      mmgWarn0 = 1;
    }
    return 0;
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
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > MMG5_EPSD2 );

  MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt    = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );
    pxt   = &mesh->xtetra[pt->xt];

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0    = MMG5_idir[ifac][i];
    itri1 = MMG5_inxt2[i];
    i1    = MMG5_idir[ifac][itri1];
    itri2 = MMG5_iprv2[i];
    i2    = MMG5_idir[ifac][itri2];
    p1    = &mesh->point[pt->v[i1]];

    /* Store the two ending points of ref curves */
    if ( MG_REF & pxt->tag[MMG5_iarf[ifac][itri1]] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i2];
      }
      else if ( !ipref[1] && (pt->v[i2] != ipref[0]) ) {
        ipref[1] = pt->v[i2];
      }
      else if ( (pt->v[i2] != ipref[0]) && (pt->v[i2] != ipref[1]) ) {
        if ( !mmgWarn1 ) {
          fprintf(stderr,"\n  ## Warning: %s: at least 1 metric not computed:"
                  " non singular point at intersection of 3 ref edges.\n",__func__);
          mmgWarn1 = 1;
        }
        return 0;
      }
    }

    if ( MG_REF & pxt->tag[MMG5_iarf[ifac][itri2]] ) {
      if ( !ipref[0] ) {
        ipref[0] = pt->v[i1];
      }
      else if ( !ipref[1] && (pt->v[i1] != ipref[0]) ) {
        ipref[1] = pt->v[i1];
      }
      else if ( (pt->v[i1] != ipref[0]) && (pt->v[i1] != ipref[1]) ) {
        if ( !mmgWarn1 ) {
          fprintf(stderr,"\n  ## Warning: %s: at least 1 metric not computed:"
                  " non singular point at intersection of 3 ref edges.\n",__func__);
          mmgWarn1 = 1;
        }
        return 0;
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
  assert ( ilists >= 1 );
  lispoi[3*ilists+1] =  lispoi[1];
  lispoi[3*ilists+2] =  lispoi[2];
  lispoi[3*ilists+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    assert(det2d);
    if ( det2d <= 0.0 ) {
      return 0;
    }
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  assert(det2d);
  if ( det2d <= 0.0 ) {
    return 0;
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
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0  = MMG5_idir[ifac][i];
    // i1  = MMG5_idir[ifac][MMG5_inxt2[i]];
    assert( 0<=ifac && ifac<4 && "unexpected local face idx");
    MMG5_tet2tri(mesh,iel,ifac,&ptt);

    MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,ifac));

    /* 1. Fill matrice tAA and second member tAb with \f$A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i})\f$ and \f$b=\sum Z_{P_i}\f$ with \f$P_i\f$
     * the physical points at edge [i0;i1] extremities and middle.
     *
     * 2. Compute the physical coor \a c of the curve edge's mid-point.
     */
    MMG5_fillDefmetregSys(k,p0,i,b,r,c,lispoi,tAA,tAb);

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
  return MMG5_solveDefmetrefSys( mesh, p0, ipref, r, c, tAA, tAb, m,
                                  isqhmin, isqhmax, hausd);
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
static int MMG5_defmetreg(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int kel,int iface, int ip) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_Tria      ptt;
  MMG5_pPoint    p0,p1;
  MMG5_pxPoint   px0;
  MMG5_Bezier    b;
  MMG5_pPar      par;
  MMG5_int       lists[MMG3D_LMAX+2];
  int64_t        listv[MMG3D_LMAX+2];
  int            k,ilist,ilists,ilistv;
  int            ifac,isloc;
  MMG5_int       iel,idp;
  double         *n,*m,r[3][3],ux,uy,uz,lispoi[3*MMG3D_LMAX+1];
  double         det2d,c[3],isqhmin,isqhmax;
  double         tAA[6],tAb[3],hausd;
  uint8_t        i1,i;
  static int8_t  mmgWarn = 0;

  pt  = &mesh->tetra[kel];
  idp = pt->v[ip];
  p0  = &mesh->point[idp];

  /* local parameters at vertex */
  hausd   = mesh->info.hausd;
  isqhmin = mesh->info.hmin;
  isqhmax = mesh->info.hmax;
  isloc     = 0;

  ilist = MMG5_boulesurfvolp(mesh,kel,ip,iface,listv,&ilistv,lists,&ilists,0);

  if ( ilist!=1 ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: at least 1 metric not computed:"
              " unable to compute the ball of point.\n",
              __func__);
      mmgWarn = 1;
    }
    return 0;
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
  assert( p0->xp && !MG_SIN(p0->tag) && !MG_EDG(p0->tag) && !(MG_NOM & p0->tag) );
  px0 = &mesh->xpoint[p0->xp];

  n  = &px0->n1[0];
  assert ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] > MMG5_EPSD2 );

  MMG5_rotmatrix(n,r);
  m = &met->m[6*idp];

  /* Apply rotation \circ translation to the whole ball */
  for (k=0; k<ilists; k++) {
    iel   = lists[k] / 4;
    ifac  = lists[k] % 4;
    pt = &mesh->tetra[iel];
    assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac] & MG_BDY) );

    for ( i = 0; i < 3; i++ ) {
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0 = MMG5_idir[ifac][i];
    i1 = MMG5_idir[ifac][MMG5_inxt2[i]];
    p1 = &mesh->point[pt->v[i1]];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    lispoi[3*k+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*k+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*k+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;
  }

  /* list goes modulo ilist */
  assert ( ilists >= 1 );
  lispoi[3*ilists+1] =  lispoi[1];
  lispoi[3*ilists+2] =  lispoi[2];
  lispoi[3*ilists+3] =  lispoi[3];

  /* Check all projections over tangent plane. */
  for (k=0; k<ilists-1; k++) {
    det2d = lispoi[3*k+1]*lispoi[3*(k+1)+2] - lispoi[3*k+2]*lispoi[3*(k+1)+1];
    assert(det2d);
    if ( det2d <= 0.0 ) {
      return 0;
    }
  }
  det2d = lispoi[3*(ilists-1)+1]*lispoi[3*0+2] - lispoi[3*(ilists-1)+2]*lispoi[3*0+1];
  assert(det2d);
  if ( det2d <= 0.0 ) {
    return 0;
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
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    // i0  = MMG5_idir[ifac][i];
    // i1  = MMG5_idir[ifac][MMG5_inxt2[i]];
    assert( 0<=ifac && iface<4 && "unexpected local face idx");
    MMG5_tet2tri(mesh,iel,ifac,&ptt);

    MMG5_bezierCP(mesh,&ptt,&b,MG_GET(pxt->ori,ifac));

    /* 1. Fill matrice tAA and second member tAb with \f$A=(\sum X_{P_i}^2 \sum
     * Y_{P_i}^2 \sum X_{P_i}Y_{P_i})\f$ and \f$b=\sum Z_{P_i}\f$ with P_i the
     * physical points at edge [i0;i1] extremities and middle.
     *
     * 2. Compute the physical coor \a c of the curve edge's mid-point.
     */
    MMG5_fillDefmetregSys(k,p0,i,b,r,c,lispoi,tAA,tAb);

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
  return MMG5_solveDefmetregSys( mesh,r, c, tAA, tAb, m, isqhmin, isqhmax,
                                  hausd);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ismet 1 if user provided metric
 * \return 1 if success, 0 otherwise.
 *
 * Define metric map at a non-boundary vertex of the mesh.
 * Allocate the metric if needed.
 * Truncate the metric at the hmin/hmax values.
 *
 */
static inline
int MMG5_defmetvol(MMG5_pMesh mesh,MMG5_pSol met,int8_t ismet) {
  MMG5_pTetra   pt,ptloc;
  MMG5_pPoint   ppt;
  MMG5_pPar     par;
  double        isqhmax,isqhmin,*m;
  MMG5_int      k,ip;
  int64_t       list[MMG3D_LMAX+2];
  int           l,i,j,isloc,ilist;

  isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

  if ( !ismet ) {

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

          /* Local parameters at tetra */
          if ( mesh->info.parTyp & MG_Tetra ) {
            ilist = MMG5_boulevolp(mesh,k,l,list);

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
    return 1;
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

        /* Local parameters at tetra */
        if ( mesh->info.parTyp & MG_Tetra ) {
          ilist = MMG5_boulevolp(mesh,k,l,list);

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
      if ( !MMG5_truncate_met3d(met,ip,isqhmin,isqhmax) ) {
        return 0;
      }
    }
  }

  return 1;
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
int MMG3D_intextmet(MMG5_pMesh mesh,MMG5_pSol met,int np,double me[6]) {
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

  return MMG5_mmgIntextmet(mesh,met,np,me,n);

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric stucture.
 * \return 0 if fail, 1 otherwise.
 *
 * Define size at points by intersecting the surfacic metric and the
 * physical metric.
 *
 *
 * 1. On singular (CRN, REQ, NOM) points, the metric on P is made isotropic.
 * 2. On non-singular ridge points, the metric is forced to be aligned with the
 *    ridge directions and surface normals.
 * 3. On regular boundary points, the metric can be anisotropic on the tangent
 *    plane, but it is forced to be aligned to the normal direction.
 */
int MMG3D_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  double        mm[6];
  MMG5_int      k,l;
  int           iploc;
  int8_t        ismet;
  int8_t        i;
  static int8_t mmgErr = 0;

  if ( !MMG5_defsiz_startingMessage (mesh,met,__func__) ) {
    return 0;
  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->flag = 0;
    ppt->s    = 0;
  }

  if ( !met->m ) {
    ismet = 0;

    /* Allocate and store the header informations for each solution */
    if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,3) ) {
      return 0;
    }

    MMG5_caltet         = MMG5_caltet_ani;
    MMG5_caltri         = MMG5_caltri_ani;
    MMG5_lenedg         = MMG5_lenedg_ani;
    MMG3D_lenedgCoor    = MMG5_lenedgCoor_ani;
    MMG5_lenSurfEdg     = MMG5_lenSurfEdg_ani;
  }
  else {
    ismet = 1;
  }

  /** Step 1: Set metric at points belonging to a required edge: compute the
   * metric as the mean of the length of the required eges passing through the
   * point */
  if ( !mesh->info.nosizreq ) {
    if ( !MMG3D_set_metricAtPointsOnReqEdges ( mesh,met,ismet ) ) {
      return 0;
    }
  }

  /* Step 2: metric definition at internal points */
  if ( !MMG5_defmetvol(mesh,met,ismet) )  return 0;

  /* Step 3: metric definition at boundary points */
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
        iploc = MMG5_idir[l][i];
        ppt   = &mesh->point[pt->v[iploc]];

        if ( !MG_VOK(ppt) )  continue;

        if ( ppt->flag > 1 ) continue;

        if ( ismet )  memcpy(mm,&met->m[6*(pt->v[iploc])],6*sizeof(double));

        if ( MG_SIN_OR_NOM(ppt->tag) ) {
          if ( !MMG5_defmetsin(mesh,met,k,l,iploc) )  continue;
        }
        else if ( ppt->tag & MG_GEO ) {
          if ( !MMG5_defmetrid(mesh,met,k,l,iploc))  continue;
        }
        else if ( ppt->tag & MG_REF ) {
          if ( !MMG5_defmetref(mesh,met,k,l,iploc) )  continue;
        } else {
          if ( !MMG5_defmetreg(mesh,met,k,l,iploc) )  continue;
        }
        if ( ismet ) {
          if ( !MMG3D_intextmet(mesh,met,pt->v[iploc],mm) ) {
            if ( !mmgErr ) {
              fprintf(stderr,"\n  ## Error: %s: unable to intersect metrics"
                      " at point %" MMG5_PRId ".\n",__func__,
                      MMG3D_indPt(mesh,pt->v[iploc]));
              mmgErr = 1;
            }
            return 0;
          }
        }
        ppt->flag = 2;
      }
    }
  }
  /* Now the metric storage at ridges is the "mmg" one. */
  mesh->info.metRidTyp = 1;

  /* search for unintialized metric */
  MMG5_defUninitSize(mesh,met,ismet);

  return 1;
}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param ip global point index.
 * \param ux edge vector x-component.
 * \param uy edge vector y-component.
 * \param uz edge vector z-component.
 * \param m point metric (copy, to be filled).
 * \param ridgedir pointer to the index of the normal direction used for ridge
 * metric (on ridge points only).
 *
 * \return 0 on failure, 1 if success.
 *
 * Get metric tensor from metric structure (pass from ridge storage to classical
 * storage on non-singular ridge points, copy metric on all other points).
 * See: \cite borouchaki1998mesh. The Hv-correction
 * is used (gradation with respect to H-variation measure).
 */
static inline
int MMG5_grad2metVol_getmet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int ip,double ux,double uy,double uz,double *m,int8_t *ridgedir) {
  MMG5_pPoint  ppt;
  MMG5_pxPoint pxp;
  double *mm,*nn1,*nn2,rbasis[3][3],ps1,ps2;

  ppt = &mesh->point[ip];
  mm  = &met->m[6*ip];


  if( MG_SIN_OR_NOM(ppt->tag) ){

    /* no normal, no tangent plane */
    memcpy(m,mm,6*sizeof(double));

  }
  else if( ppt->tag & MG_GEO ) {

    /* Recover normal and metric */
    pxp = &mesh->xpoint[ppt->xp];
    nn1 = pxp->n1;
    nn2 = pxp->n2;
    ps1 = ux*nn1[0] + uy*nn1[1] + uz*nn1[2];
    ps2 = ux*nn2[0] + uy*nn2[1] + uz*nn2[2];
    if ( fabs(ps2)<fabs(ps1) ) {
      *ridgedir = 1;
    }
    else{
      *ridgedir = 0;
    }
    /* Note that rbasis is not used in this function */
    if( !MMG5_buildridmet(mesh,met,ip,ux,uy,uz,m,rbasis) )
      return 0;

  }
  else if( ppt->tag & MG_BDY ) {
    memcpy(m,mm,6*sizeof(double));
  }
  else {

    /* internal point */
    memcpy(m,mm,6*sizeof(double));

  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh.
 * \param ppt pointer toward the P point structure.
 * \param l Euclidean length of the edge PQ.
 * \param m metric tensor on point P.
 * \param mext metric tensor extended on Q.
 *
 * Extend metric tensor from point P to point Q according to an h-variation law
 * applied on each eigenvalue.
 *
 */
static inline
void MMG5_grad2metVol_extmet(MMG5_pMesh mesh,MMG5_pPoint ppt,double l,double *m,double *mext) {
  double lambda[3],vp[3][3];

  MMG5_eigenv3d(1,m,lambda,vp);
  for( int8_t i = 0; i < 3; i++ ) {
    double hext = 1./sqrt(lambda[i]) + mesh->info.hgrad*l + MMG5_EPSOK;
    lambda[i] = 1./(hext*hext);
  }
  MMG5_eigenvmatsym3d(mesh,mext,lambda,vp);

}

/**
 * \param mesh pointer toward the mesh.
 * \param ppt pointer toward the P point structure.
 * \param m metric tensor on point P.
 * \param mext extended metric tensor from Q to P.
 * \param ridgedir normal direction for metric reconstruction on P (on ridge only).
 * \param iloc local index of point P on the edge.
 * \param ier pointer to the local indices (on the edge) of updated points, with bitwise encoding.
 *
 * Use metric intersection to gradate the anisotropic metric on point P, given
 * the metric extended from point Q on the edge PQ.
 * 1. On singular (CRN, REQ, NOM) points, the metric on P is isotropic (as in MMG5_defmetsin)
 *    and should remain isotropic.
 * 2. On non-singular ridge points, the metric on P should remain aligned with
 *    the ridge directions (thus it is not possible to apply metric intersection,
 *    sizes are truncated instead).
 * 3. On regular boundary points, the metric can be anisotropic on the tangent
 *    plane, but it should remain aligned to the normal direction (thus,
 *    intersection is used only in the tangent plane and sizes are truncated in
 *    the normal direction).
 * 4. On interior volume points, 3D metric intersection can be used.
 *
 */
static inline
void MMG3D_gradSimred(MMG5_pMesh mesh,MMG5_pPoint ppt,double m[6],double mext[6],int8_t ridgedir,int8_t iloc,int *ier) {
  double tol = MMG5_EPS;

  if ( MG_SIN_OR_NOM(ppt->tag) ) {
    double dm[3],dmext[3],vp[3][3],beta;

    /* Simultaneous reduction basis */
    if( !MMG5_simred3d(mesh,m,mext,dm,dmext,vp) ) {
      *ier = -1;
      return;
    }
    beta = 1.0;
    /* Compute maximum size variation */
    for( int8_t i = 0; i < 3; i++ ) {
      beta = MG_MAX(beta,dmext[i]/dm[i]);
    }
    if( beta > 1.0 + tol ) {
      (*ier) = (*ier) | iloc;
    }
    if( (*ier) & iloc ) {
      /* You have to homothetically scale the metric in order to make it
       * consistent with the above check and for the gradation loop to converge.
       * Change the gradation law in the metric extension in
       * MMG5_grad2metVol_extmet if you want to restrict the influence of
       * singular points. */
      m[0] *= beta;
      m[3] *= beta;
      m[5] *= beta;
    }

  }
  else if ( ppt->tag & MG_GEO ) {
    MMG5_pxPoint pxp;
    double mr[6],mrext[6],dm[3],dmext[3];
    double u[3],r[3][3],*t,*n;

    /* Compute ridge orthonormal basis (t, n x t, n) */
    t = ppt->n;
    pxp = &mesh->xpoint[ppt->xp];
    if( ridgedir )
      n = pxp->n2;
    else
      n = pxp->n1;
    MMG5_crossprod3d(n,t,u);
    /* Store basis in rotation matrix */
    for( int8_t i = 0; i < 3; i++ ) {
      r[0][i] = t[i];
      r[1][i] = u[i];
      r[2][i] = n[i];
    }
    /* Rotate matrices to the ridge reference system */
    MMG5_rmtr(r,m,mr);
    dm[0] = mr[0];
    dm[1] = mr[3];
    dm[2] = mr[5];
    MMG5_rmtr(r,mext,mrext);
    dmext[0] = mrext[0];
    dmext[1] = mrext[3];
    dmext[2] = mrext[5];
    /* The ridge metrics is diagonal and must remain diagonal, while the
     * extended metric from the other point  is not necessarily diagonal in the
     * ridge basis. Do not perform intersection, but simply evaluate lengths on
     * the ridge basis and truncate sizes. */
    for( int8_t i = 0; i < 3; i++ ) {
      if( dmext[i] > dm[i]*(1.0 + tol) ) {
        dm[i] = dmext[i];
        (*ier) = (*ier) | iloc;
      }
    }
    /* Update metric */
    if( (*ier) & iloc ) {
      /* Re-build matrix */
      MMG5_eigenvmatsym3d(mesh,m,dm,r);
    }

  }
  else if( ppt->tag & MG_BDY ) {
    MMG5_pxPoint pxp;
    double mr[6],mrext[6],mtan[3],mtanext[3],r[3][3];
    double dm[3],dmext[3],vp[2][2];

    /* Rotation matrices mapping n to e_3 */
    pxp = &mesh->xpoint[ppt->xp];
    MMG5_rotmatrix(pxp->n1,r);

    /* Rotate metrics to the tangent plane */
    MMG5_rmtr(r,m,mr);
    mtan[0] = mr[0];
    mtan[1] = mr[1];
    mtan[2] = mr[3];
    MMG5_rmtr(r,mext,mrext);
    mtanext[0] = mrext[0];
    mtanext[1] = mrext[1];
    mtanext[2] = mrext[3];
    /* The current point metric is aligned with the surface normal direction,
     * while the extended metric from the other point need not.
     * Perform intersection only in the tangent plane, while simply truncate
     * sizes in the normal direction.
     *
     * Simultaneous reduction basis */
    if( !MMG5_simred2d(mesh,mtan,mtanext,dm,dmext,vp) ) {
      *ier = -1;
      return;
    }
    /* Intersection in the tangent plane */
    for( int8_t i = 0; i < 2; i++ ) {
      if( dmext[i] > dm[i]*(1.0 + tol) ) {
        dm[i] = dmext[i];
        (*ier) = (*ier) | iloc;
      }
    }
    /* The current point metric is aligned with the surface normal direction,
     * while the extended metric from the other point need not.
     * Simply evaluate lengths in the normal direction and truncate sizes. */
    dm[2]    = mr[5];
    dmext[2] = mrext[5];
    if( dmext[2] > dm[2]*(1.0 + tol) ) {
      dm[2] = dmext[2];
      (*ier) = (*ier) | iloc;
    }
    /* Update metric */
    if( (*ier) & iloc ) {
      /* Simultaneous reduction basis is non-orthogonal, so invert it for the
       * inverse transformation for the tangent-plane metric. */
      double ivp[2][2];
      if( !MMG5_invmat22(vp,ivp) ) {
        *ier = -1;
        return;
      }
      MMG5_simredmat(2,mtan,dm,(double *)ivp);
      /* Re-assemble 3D metric: use the intersected metrics in the tangent
       * plane, and the truncated size in the normal direction. */
      mr[0] = mtan[0];
      mr[1] = mtan[1];
      mr[3] = mtan[2];
      mr[2] = mr[4] = 0.0;
      mr[5] = dm[2];
      /* Transpose rotation matrix and rotate back into the point metric*/
      MMG5_transpose3d(r);
      MMG5_rmtr(r,mr,m);
    }

  }
  else { /* internal point */
    double dm[3],dmext[3],vp[3][3];

    /* Simultaneous reduction basis */
    if( !MMG5_simred3d(mesh,m,mext,dm,dmext,vp) ) {
      *ier = -1;
      return;
    }

    /* Gradation of sizes in the simultaneous reduction basis */
    for( int8_t i = 0; i< 3; i++ ) {
      if( dmext[i] > dm[i]*(1.0 + tol) ) {
        dm[i] = dmext[i];
        (*ier) = (*ier) | iloc;
      }
    }
    /* Update metric */
    if( (*ier) & iloc ) {
      /* Simultaneous reduction basis is non-orthogonal, so invert it for the
       * inverse transformation */
      double ivp[3][3];
      if( !MMG5_invmat33(vp,ivp) ) {
        *ier = -1;
        return;
      }
      MMG5_simredmat(3,m,dm,(double *)ivp);
    }

  }

}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param ip global index of the point.
 * \param m metric tensor on the point (copy).
 * \param ridgedir normal direction for metric reconstruction (on ridge only).
 *
 * Set new metric tensor to the metric structure (pass from classical storage to
 * ridge storage on non-singular ridge points, copy metric on all other
 * points).
 *
 */
static inline
void MMG5_grad2metVol_setmet(MMG5_pMesh mesh,MMG5_pSol met,int ip,double *m,int8_t ridgedir) {
  MMG5_pPoint ppt;
  double *mm;

  ppt = &mesh->point[ip];
  mm = &met->m[6*ip];

  if ( MG_RID(ppt->tag) ) {
    double mr[6];
    double u[3],r[3][3],*t,*n;

    /* Compute ridge orthonormal basis (t, n x t, n) */
    t = ppt->n;
    MMG5_pxPoint pxp = &mesh->xpoint[ppt->xp];
    if( ridgedir )
      n = pxp->n2;
    else
      n = pxp->n1;
    MMG5_crossprod3d(n,t,u);
    /* Store basis in rotation matrix */
    for( int8_t i = 0; i < 3; i++ ) {
      r[0][i] = t[i];
      r[1][i] = u[i];
      r[2][i] = n[i];
    }
    /* Rotate matrices to the ridge reference system */
    MMG5_rmtr(r,m,mr);
    /* Store matrix in ridge format */
    mm[0]          = mr[0];
    mm[1+ridgedir] = mr[3];
    mm[3+ridgedir] = mr[5];

  } else {

    memcpy(mm,m,6*sizeof(double));

  }

  return;
}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param np1 global index of the first edge extremity.
 * \param np2 global index of the second edge extremity.
 *
 * \return -1 on failure, else local indices of graded point on the edge
 * (bitwise encoded, so 0 for no update, 1 if the first point iis updated, 2 if
 * the second one is updated, 3 if both). Use the third bit to switch on warning
 * information at output.
 *
 * Enforces gradation of metric in both extremities of an edge with respect to
 * one other.
 *
 */
static inline
int MMG5_grad2metVol(MMG5_pMesh mesh,MMG5_pSol met,int np1,int np2) {
  MMG5_pPoint    p1,p2;
  double         m1[6],m2[6],mext1[6],mext2[6];
  double         ux,uy,uz,l;
  int8_t         ridgedir1,ridgedir2;
  int            ier = 0;

  p1  = &mesh->point[np1];
  p2  = &mesh->point[np2];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  uz = p2->c[2] - p1->c[2];
  l = sqrt(ux*ux+uy*uy+uz*uz);


  /** Recover normal and metric associated to p1 and p2 (metric can be in ridge
   * storage) */
  if( !MMG5_grad2metVol_getmet(mesh,met,np1,ux,uy,uz,m1,&ridgedir1) ) {
    return -1;
  }
  if( !MMG5_grad2metVol_getmet(mesh,met,np2,ux,uy,uz,m2,&ridgedir2) ) {
    return -1;
  }
  /* (metric now follows standard 3D storage) */


  /** Extend p2 metric and gradate p1 */
  if( p2->flag >= mesh->base-1 ) {
    /* Extend p2 metrics */
    MMG5_grad2metVol_extmet(mesh,p2,l,m2,mext2);

    /* Gradate p1 metrics */
    MMG3D_gradSimred(mesh,p1,m1,mext2,ridgedir1,1,&ier);
    if( ier == -1 )
      return ier;
#ifndef NDEBUG
    double mtmp[6];
    int iertmp = 0;
    memcpy(mtmp,m1,6*sizeof(double));
    MMG3D_gradSimred(mesh,p1,mtmp,mext2,ridgedir1,1,&iertmp);
    if( iertmp & 1 )
      ier |= 4;
#endif
  }


  /** Extend p1 metric and gradate p2 (p1 has already been updated) */
  if( p1->flag >= mesh->base-1 ) {
    /* Expand p1 metrics */
    MMG5_grad2metVol_extmet(mesh,p1,l,m1,mext1);

    /* Gradate p2 metrics */
    MMG3D_gradSimred(mesh,p2,m2,mext1,ridgedir2,2,&ier);
    if( ier == -1 )
      return ier;
#ifndef NDEBUG
    double mtmp[6];
    int iertmp = 0;
    memcpy(mtmp,m2,6*sizeof(double));
    MMG3D_gradSimred(mesh,p2,mtmp,mext1,ridgedir2,2,&iertmp);
    if( iertmp & 2 )
      ier |= 4;
#endif
  }


  /* Set metrics to the met structure, back to ridge storage */
  if( ier & 1 )
    MMG5_grad2metVol_setmet(mesh,met,np1,m1,ridgedir1);
  if( ier & 2 )
    MMG5_grad2metVol_setmet(mesh,met,np2,m2,ridgedir2);

  return ier;
}


/**
 * \param n  matrix to update
 * \param dn eigenvalues of n in the coreduction basis
 * \param vp coreduction basis
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update of the metric n = tP^-1 diag(dn0,dn1,dn2)P^-1, P = (vp[0],vp[1],vp[2])
 * stored in columns
 *
 */
int MMG3D_updatemetreq_ani(double *n,double dn[3],double vp[3][3]) {
  double ip[3][3];

  /* P^-1 */
  if ( !MMG5_invmat33(vp,ip) ) {
    return 0;
  }

  /* tp^-1 * dn * P^-1 */
  n[0] = dn[0]*ip[0][0]*ip[0][0] + dn[1]*ip[1][0]*ip[1][0] + dn[2]*ip[2][0]*ip[2][0];
  n[1] = dn[0]*ip[0][0]*ip[0][1] + dn[1]*ip[1][0]*ip[1][1] + dn[2]*ip[2][0]*ip[2][1];
  n[2] = dn[0]*ip[0][0]*ip[0][2] + dn[1]*ip[1][0]*ip[1][2] + dn[2]*ip[2][0]*ip[2][2];

  n[3] = dn[0]*ip[0][1]*ip[0][1] + dn[1]*ip[1][1]*ip[1][1] + dn[2]*ip[2][1]*ip[2][1];
  n[4] = dn[0]*ip[0][1]*ip[0][2] + dn[1]*ip[1][1]*ip[1][2] + dn[2]*ip[2][1]*ip[2][2];

  n[5] = dn[0]*ip[0][2]*ip[0][2] + dn[1]*ip[1][2]*ip[1][2] + dn[2]*ip[2][2]*ip[2][2];

  return 1;
}

/**
 * \param mesh pointer toward the mesh.
 * \param met pointer toward the metric structure.
 * \param pt pointer toward a tetra.
 * \param npmaster edge extremity that cannot be modified
 * \param npslave edge extremity to modify to respect the gradation.
 *
 * \return 1 if the graddation succeed, 0 otherwise
 *
 * Enforces gradation of metric from required entity toward other using the
 * simultaneous reduction technique (note that as the gradation is propagated,
 * we can be on an edge without a required extremity).
 *
 */
static inline
int MMG5_grad2metVolreq(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt,MMG5_int npmaster,
                        MMG5_int npslave) {
  MMG5_pPoint    p1,p2;
  double         *mm1,*mm2,m1[6],m2[6],ux,uy,uz;
  double         l,difsiz,rbasis1[3][3],rbasis2[3][3];
  double         lambda[3],vp[3][3],beta,mu[3];
  int            cfg_m2;
  int8_t         ier;

  p1  = &mesh->point[npmaster];
  p2  = &mesh->point[npslave];

  ux = p2->c[0] - p1->c[0];
  uy = p2->c[1] - p1->c[1];
  uz = p2->c[2] - p1->c[2];

  mm1  = &met->m[6*npmaster];
  mm2  = &met->m[6*npslave];

  cfg_m2 = 0;
  ier = 0;

  if ( MG_RID(p1->tag) ) {
    if ( MG_RID(p2->tag) ) {
      // The volume gradation from ridge point toward another ridge point is
      // bugged...
      return 0;
    }

    /* Recover normal and metric associated to p1 */
    if( !MMG5_buildridmet(mesh,met,npmaster,ux,uy,uz,m1,rbasis1) ) { return 0; }
  }
  else {
    memcpy(m1,mm1,6*sizeof(double));
  }

  if ( MG_RID(p2->tag) ) {
    /* Recover normal and metric associated to p2 */
    cfg_m2 = MMG5_buildridmet(mesh,met,npslave,ux,uy,uz,m2,rbasis2);
    if( !cfg_m2 ) { return 0; }
  }
  else {
    memcpy(m2,mm2,6*sizeof(double));
  }

  l = sqrt(ux*ux+uy*uy+uz*uz);

  difsiz = mesh->info.hgradreq*l;

  /* Simultaneous reduction of mtan1 and mtan2 */
  if ( !MMG5_simred3d(mesh,m1,m2,lambda,mu,vp) ) {
    return 0;
  }

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the first direction */
  MMG5_gradEigenvreq(lambda,mu,difsiz,0,&ier);

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the second direction */
  MMG5_gradEigenvreq(lambda,mu,difsiz,1,&ier);

  /* Gradation of sizes = 1/sqrt(eigenv of the tensors) in the third direction */
  MMG5_gradEigenvreq(lambda,mu,difsiz,2,&ier);

  if ( !ier ) {
    return 0;
  }

  /* Metric update using the simultaneous reduction technique */
  if( MG_SIN_OR_NOM(p2->tag) ) {
    /* We choose to not respect the gradation in order to restrict the influence
     * of the singular points. Thus:
     * lambda_new = = 0.5 lambda_1 + 0.5 lambda_new = lambda_1 + 0.5 beta.
     * with beta the smallest variation of the eigenvalues (lambda_new-lambda_1). */

    /* This point can have an anisotropic metric if a user-provided metric is
     * found. So, compute the eigenvalues. */
    double ll[3],rr[3][3],llmin;
    int i;
    if( !MMG5_eigenv3d(1,mm2,ll, rr) )
      return 0;
    llmin = DBL_MAX;
    for( i = 0; i < 3; i++ )
      if( ll[i] < llmin )
        llmin = ll[i];


    beta = mu[0] - llmin;

    if ( fabs(beta) < fabs(llmin-mu[1]) ) {
      beta = mu[1] - llmin;
    }
    ll[0] += 0.5*beta;
    ll[1] += 0.5*beta;
    ll[2] += 0.5*beta;
    assert ( ll[0]>0. && ll[1]>0. && ll[2]>0. );
    mm2[0] = ll[0]*rr[0][0]*rr[0][0] + ll[1]*rr[1][0]*rr[1][0] + ll[2]*rr[2][0]*rr[2][0];
    mm2[1] = ll[0]*rr[0][0]*rr[0][1] + ll[1]*rr[1][0]*rr[1][1] + ll[2]*rr[2][0]*rr[2][1];
    mm2[2] = ll[0]*rr[0][0]*rr[0][2] + ll[1]*rr[1][0]*rr[1][2] + ll[2]*rr[2][0]*rr[2][2];
    mm2[3] = ll[0]*rr[0][1]*rr[0][1] + ll[1]*rr[1][1]*rr[1][1] + ll[2]*rr[2][1]*rr[2][1];
    mm2[4] = ll[0]*rr[0][1]*rr[0][2] + ll[1]*rr[1][1]*rr[1][2] + ll[2]*rr[2][1]*rr[2][2];
    mm2[5] = ll[0]*rr[0][2]*rr[0][2] + ll[1]*rr[1][2]*rr[1][2] + ll[2]*rr[2][2]*rr[2][2];
  }
  else if( p2->tag & MG_GEO ){

    if ( !MMG3D_updatemetreq_ani(m2,mu,vp) ) { return 0; }

    /* Here mtan2 contains the gradated metric in the coreduction basis: compute
     * the sizes in the directions (t,u=t^n,n) */
    mu[0] = m2[0]*rbasis2[0][0]*rbasis2[0][0] + 2. * m2[1]*rbasis2[1][0]*rbasis2[0][0]
      + 2. * m2[2]*rbasis2[2][0]*rbasis2[0][0]
      + m2[3]*rbasis2[1][0]*rbasis2[1][0] + 2. * m2[4]*rbasis2[2][0]*rbasis2[1][0]
      + m2[5]*rbasis2[2][0]*rbasis2[2][0];

    /* h = 1/sqrt(t_e M e) */
    assert ( mu[0] > MMG5_EPSD2 );

    mu[1] = m2[0]*rbasis2[0][1]*rbasis2[0][1] + 2. * m2[1]*rbasis2[1][1]*rbasis2[0][1]
      + 2. * m2[2]*rbasis2[2][1]*rbasis2[0][1]
      + m2[3]*rbasis2[1][1]*rbasis2[1][1] + 2. * m2[4]*rbasis2[2][1]*rbasis2[1][1]
      + m2[5]*rbasis2[2][1]*rbasis2[2][1];

    /* h = 1/sqrt(t_e M e) */
    assert ( mu[1] > MMG5_EPSD2 );

    mu[2] = m2[0]*rbasis2[0][2]*rbasis2[0][2] + 2. * m2[1]*rbasis2[1][2]*rbasis2[0][2]
      + 2. * m2[2]*rbasis2[2][2]*rbasis2[0][2]
      + m2[3]*rbasis2[1][2]*rbasis2[1][2] + 2. * m2[4]*rbasis2[2][2]*rbasis2[1][2]
      + m2[5]*rbasis2[2][2]*rbasis2[2][2];

    /* h = 1/sqrt(t_e M e) */
    assert ( mu[2] > MMG5_EPSD2 );

    /* Update the ridge metric */
    mm2[0] =  mu[0];

    assert ( cfg_m2 );

    mm2[cfg_m2] = mu[1];
    mm2[cfg_m2+2] = mu[2];
  }
  else{

    if ( !MMG3D_updatemetreq_ani(m2,mu,vp) ) { return 0; }
    memcpy(mm2,m2,6*sizeof(double));
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 1
 *
 *
 * Enforces mesh gradation by truncating metric field.
 *
 * 1. On singular (CRN, REQ, NOM) points, the metric on P is isotropic
 *    (as in MMG5_defmetsin) and should remain isotropic.
 * 2. On non-singular ridge points, the metric on P should remain aligned with
 *    the ridge directions (thus it is not possible to apply metric intersection,
 *    sizes are truncated instead).
 * 3. On regular boundary points, the metric can be anisotropic on the tangent
 *    plane, but it should remain aligned to the normal direction (thus,
 *    intersection is used only in the tangent plane and sizes are truncated in
 *    the normal direction).
 * 4. On interior volume points, 3D metric intersection can be used.
 *
 */
int MMG3D_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_Hash     edgeTable;
  MMG5_hedge    *pht;
  MMG5_pTetra   pt;
  MMG5_pPoint   p0,p1;
  double        *m,mv;
  int           i,itv,maxit,ier;
  MMG5_int      k,np0,np1,nu,nupv;
  static int    mmgWarn = 0;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Anisotropic mesh gradation\n");

  /* First step : make ridges iso in each apairing direction */
  for (k=1; k<= mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !MG_VOK(p1) ) continue;
    if ( MG_SIN_OR_NOM(p1->tag) ) continue;
    if ( !(p1->tag & MG_GEO) ) continue;

    m = &met->m[6*k];
    mv = MG_MAX(m[1],m[2]);
    m[1] = mv;
    m[2] = mv;
    mv = MG_MAX(m[3],m[4]);
    m[3] = mv;
    m[4] = mv;
  }

  /** Mark the edges belonging to a required entity */
  MMG3D_mark_pointsOnReqEdge_fromTetra ( mesh );


  /* alloc hashtable */
  if ( !MMG5_hashNew(mesh,&edgeTable,mesh->nemax,3*mesh->nemax) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate hash table.\n",__func__);
    return 0;
  }

  /* build edge table */
  for(k=1 ; k<=mesh->ne ; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) {
      continue;
    }
    for(i=0 ; i<6 ; i++) {
      np0 = pt->v[MMG5_iare[i][0]];
      np1 = pt->v[MMG5_iare[i][1]];

      ier = MMG5_hashEdge(mesh,&edgeTable,np0,np1,k);
      if ( !ier ) {
        if ( !mmgWarn ) {
          mmgWarn = 1;
          fprintf(stderr,"\n  ## Warning: %s: unable to hash at least one edge"
                  " (tria %" MMG5_PRId ", edge %d).\n",__func__,MMG3D_indElt(mesh,k),i);
        }
      }
    }
  }

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  nupv = itv = 0;
  maxit = 500;
  mmgWarn = 0;

  /* analyze mesh edges via hash table */
  do {
    ++mesh->base;
    nu = 0;
    for (k=0; k<edgeTable.siz; k++) {
      pht = &edgeTable.item[k];
      /* analyze linked list */
      while ( pht ) {
        if ( !pht->a )  break;
        np0  = pht->a;
        np1  = pht->b;
        p0 = &mesh->point[np0];
        p1 = &mesh->point[np1];

        /* Skip edge if both nodes have been updated more than 1 iteration ago */
        if ( (p0->flag < mesh->base-1) && (p1->flag < mesh->base-1) ) {
          pht = pht->nxt ? &edgeTable.item[pht->nxt] : 0;
          continue;
        }

        /* Skip points belonging to a required edge */
        if ( p0->s || p1->s ) {
          pht = pht->nxt ? &edgeTable.item[pht->nxt] : 0;
          continue;
        }

        ier = MMG5_grad2metVol(mesh,met,np0,np1);
        if( ier == -1 ) {
          break;
        } else {
          if ( ier & 1 ) {
            p0->flag = mesh->base;
            nu++;
          }
          if ( ier & 2 ) {
            p1->flag = mesh->base;
            nu++;
          }
          if ( !mmgWarn && (ier & 4) ) {
            mmgWarn = itv;
          }
        }

        /* next edge */
        pht = pht->nxt ? &edgeTable.item[pht->nxt] : 0;
      }
    }
    nupv += nu;
  } while ( ++itv < maxit && nu > 0 );
  MMG5_SAFE_FREE(edgeTable.item);

  if ( abs(mesh->info.imprim) > 3 ) {
    if( mmgWarn ) {
      fprintf(stderr,"\n      ## Warning: %s: Non-idempotent metric"
                     " intersections since iteration %d.\n",__func__,mmgWarn);
    }

    fprintf(stdout,"    gradation: %7" MMG5_PRId " updated, %d iter\n",nupv,itv);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Enforce mesh gradation by truncating size map.
 *
 */
int MMG3D_gradsizreq_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_Tria     ptt;
  MMG5_pPoint   p0,p1;
  int           it,itv,maxit;
  int           i,j,ier;
  MMG5_int      nup,nu,nupv,k,np0,np1,npmaster,npslave;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  ** Grading required points.\n");
  }

  if ( mesh->info.hgrad < 0. ) {
    /** Mark the edges belonging to a required entity (already done if the
     * classic gradation is enabled) */
    MMG3D_mark_pointsOnReqEdge_fromTetra ( mesh );
  }

  it = nup = 0;
  maxit = 100;
  do {
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
            MMG5_tet2tri(mesh,k,i,&ptt);
            for (j=0; j<3; j++) {
              np0 = ptt.v[MMG5_inxt2[j]];
              np1 = ptt.v[MMG5_iprv2[j]];
              p0  = &mesh->point[np0];
              p1  = &mesh->point[np1];

              if ( MMG5_abs ( p0->s - p1->s ) < 2 ) {
                /* No size to propagate */
                continue;
              }
              else if ( p0->s > p1->s ) {
                npmaster = np0;
                npslave  = np1;
              }
              else {
                assert ( p1->s > p0->s );
                npmaster = np1;
                npslave  = np0;
              }

              /* Impose the gradation to npslave from npmaster */
              /* gradation along the tangent plane */
              ier = MMG5_grad2metSurfreq(mesh,met,&ptt,npmaster,npslave);
              if ( ier ) {
                mesh->point[npslave].s = mesh->point[npmaster].s - 1;
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

  nupv = itv = 0;
  maxit = 500;

  do {
    nu = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      for (i=0; i<4; i++) {
        /* Gradation along a volume edge */
        np0  = pt->v[MMG5_iare[i][0]];
        np1  = pt->v[MMG5_iare[i][1]];
        p0  = &mesh->point[np0];
        p1  = &mesh->point[np1];

        if ( MMG5_abs ( p0->s - p1->s ) < 2 ) {
          /* No size to propagate */
          continue;
        }
        else if ( p0->s > p1->s ) {
          npmaster = np0;
          npslave  = np1;
        }
        else {
          assert ( p1->s > p0->s );
          npmaster = np1;
          npslave  = np0;
        }

        /* Impose the gradation to npslave from npmaster */
        ier =  MMG5_grad2metVolreq(mesh,met,pt,npmaster,npslave);
        if ( ier ) {
          mesh->point[npslave].s = mesh->point[npmaster].s - 1;

          nu++;
        }
      }
    }
    nupv += nu;
  }
  while( ++itv < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 3 ) {
    if ( abs(mesh->info.imprim) < 5 && !mesh->info.ddebug ) {
      fprintf(stdout,"    gradation: %7" MMG5_PRId " updated, %d iter\n",nup+nupv,it+itv);
    }
    else {
      fprintf(stdout,"    surface gradation: %7" MMG5_PRId " updated, %d iter\n"
              "    volume gradation:  %7" MMG5_PRId " updated, %d iter\n",nup,it,nupv,itv);
    }
  }

  return 1;
}
