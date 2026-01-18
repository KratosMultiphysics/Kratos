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
 * \file mmg3d/swap_3d.c
 * \brief Functions for swapping process over boundary.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg3d.h"
#include "inlined_functions_3d_private.h"
#include "mmg3dexterns_private.h"

extern int8_t ddb;

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param list pointer to the shell of the edge.
 * \param ilist pointer to the size of the shell of the edge.
 * \param it1 first element of the open shell.
 * \param it2 last element of the open shell.
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 * \return -1 if fail, 0 if we can not swap the edge, 1 otherwise.
 *
 * Check whether edge whose shell is provided should be swapped for
 * geometric approximation purposes (the 2 surface triangles are also
 * provided).
 *
 */
int MMG5_chkswpbdy(MMG5_pMesh mesh, MMG5_pSol met, int64_t *list,int ilist,
                    MMG5_int it1,MMG5_int it2,int8_t typchk) {
  MMG5_pTetra   pt,pt0;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1,ppt0;
  MMG5_Tria     tt1,tt2;
  MMG5_pPar     par;
  double        b0[3],b1[3],n[3],v[3],c[3],ux,uy,uz,ps,disnat,dischg;
  double        cal1,cal2,calnat,calchg,calold,calnew,caltmp,hausd;
  MMG5_int      iel,iel1,iel2,np,nq,na1,na2,k,nminus,nplus,info,refm;
  int           isloc,l;
  int8_t        ifa1,ifa2,ia,ip,iq,ia1,ia2,j,isshell,ier;

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];
  pt0 = &mesh->tetra[0];
  ppt0= &mesh->point[0];
  memset(ppt0,0,sizeof(MMG5_Point));

  np = pt->v[MMG5_iare[ia][0]];
  nq = pt->v[MMG5_iare[ia][1]];

  // Algiane 05/04/24: I think that the assumption that was previously made that
  // we can arrive from a tetrahedra without a boundary face (i.e. without an
  // xtetra) never happens
  assert ( pt->xt && "Boundary edges have to be swapped from a boundary face" );

  /* No swap of geometric edge */
  pxt = &mesh->xtetra[pt->xt];
  if ( (pxt->edg[ia]>0) || MG_EDG_OR_NOM(pxt->tag[ia]) || (pxt->tag[ia] & MG_REQ) ) {
    return 0;
  }

  /* No swap when either internal or external component has only 1 element (as
   * we can't swap geometric edges here we know that the edge shares at most 2
   * domains).*/
  nminus = nplus = 0;
  refm   = pt->ref;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    pt = &mesh->tetra[iel];
    if ( pt->ref == refm ) {
      nminus++;
    }
    else {
      nplus++;
    }
  }
  if ( nplus == 1 || nminus == 1 )  return 0;

  iel1 = it1 / 4;
  ifa1 = it1 % 4;

  assert(it2);
  iel2 = it2 / 4;
  ifa2 = it2 % 4;
  assert( 0<=ifa1 && ifa1<4 && "unexpected local face idx");
  assert( 0<=ifa2 && ifa2<4 && "unexpected local face idx");
  MMG5_tet2tri(mesh,iel1,ifa1,&tt1);
  MMG5_tet2tri(mesh,iel2,ifa2,&tt2);

  for (ia1=0; ia1<3; ia1++) {
    if ( (tt1.v[ia1] != np) && (tt1.v[ia1] != nq) )  break;
  }
  assert( ia1 < 3 );
  if ( ia1==3 ) return 0;

  assert( (tt1.v[MMG5_inxt2[ia1]] == np && tt1.v[MMG5_iprv2[ia1]] == nq) ||
          (tt1.v[MMG5_inxt2[ia1]] == nq && tt1.v[MMG5_iprv2[ia1]] == np) );
  na1 = tt1.v[ia1];

  for (ia2=0; ia2<3; ia2++) {
    if ( (tt2.v[ia2] != np) && (tt2.v[ia2] != nq) )  break;
  }

  assert ( ia2 < 3 );
  if ( ia2 ==3 ) return 0;

  assert ( (tt2.v[MMG5_inxt2[ia2]] == np && tt2.v[MMG5_iprv2[ia2]] == nq) ||
           (tt2.v[MMG5_inxt2[ia2]] == nq && tt2.v[MMG5_iprv2[ia2]] == np) );
  na2 = tt2.v[ia2];

  /* Check non convexity (temporarily use b0,b1)*/
  MMG5_norpts(mesh,tt1.v[ia1],tt1.v[MMG5_inxt2[ia1]],tt2.v[ia2],b0);
  MMG5_norpts(mesh,tt2.v[ia2],tt2.v[MMG5_inxt2[ia2]],tt1.v[ia1],b1);
  ps = b0[0]*b1[0] + b0[1]*b1[1] + b0[2]*b1[2];

  /* Here we put ANGEDG because in nr mode the test over dhd may create inverted
   * tetra */
  if ( ps < MMG5_ANGEDG ) {
    return 0;
  }

  /* Check normal deviation with neighbours */
  if ( !MG_GEO_OR_NOM( tt1.tag[MMG5_iprv2[ia1]] ) ) {
    ier = MMG3D_normalAdjaTri(mesh,iel1,ifa1,MMG5_iprv2[ia1],n);
    if ( ier < 0 ) return -1;
    else if ( !ier ) return 0;
    ps = b0[0]*n[0] + b0[1]*n[1] + b0[2]*n[2];

    if ( ps < mesh->info.dhd )  return 0;
  }

  if ( !MG_GEO_OR_NOM( tt2.tag[MMG5_inxt2[ia2]]) ) {
    ier = MMG3D_normalAdjaTri(mesh,iel2,ifa2,MMG5_inxt2[ia2],n);
    if ( ier<0 ) return -1;
    else if ( !ier ) return 0;
    ps = b0[0]*n[0] + b0[1]*n[1] + b0[2]*n[2];

    if ( ps < mesh->info.dhd )  return 0;
  }

  if ( !MG_GEO_OR_NOM( tt1.tag[MMG5_inxt2[ia1]] ) ) {
    ier = MMG3D_normalAdjaTri(mesh,iel1,ifa1,MMG5_inxt2[ia1],n);
    if ( ier<0 ) return -1;
    else if ( !ier ) return 0;
    ps = b1[0]*n[0] + b1[1]*n[1] + b1[2]*n[2];

    if ( ps < mesh->info.dhd )  return 0;
  }

  if ( !MG_GEO_OR_NOM(tt2.tag[MMG5_iprv2[ia2]]) ) {
    ier = MMG3D_normalAdjaTri(mesh,iel2,ifa2,MMG5_iprv2[ia2],n);
    if ( ier<0 ) return -1;
    else if ( !ier ) return 0;
    ps = b1[0]*n[0] + b1[1]*n[1] + b1[2]*n[2];

    if ( ps < mesh->info.dhd )  return 0;
  }

  /* Compare contributions to Hausdorff distance in both configurations */
  MMG5_norface(mesh,iel1,ifa1,v);

  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];

  /* local parameters */
  hausd = mesh->info.hausd;
  isloc = 0;

  /* Local params at triangles containing the edge */
  if ( mesh->info.parTyp & MG_Tria ) {
    if ( tt1.ref == tt2.ref ) {
      for ( l=0; l<mesh->info.npar; ++l ) {
        par = &mesh->info.par[l];
        if ( par->elt != MMG5_Triangle ) continue;

        hausd   = par->hausd;
        isloc   = 1;
        break;
      }
    }
    else {
      l = 0;
      info = -1000;
      do {
        if ( isloc ) break;
        par = &mesh->info.par[l];
        if ( par->elt != MMG5_Triangle ) continue;

        if ( tt1.ref!=par->ref && tt2.ref !=par->ref )  continue;

        hausd   = par->hausd;
        isloc   = 1;
        info = par->ref;
      } while ( ++l < mesh->info.npar );

      for ( ; l<mesh->info.npar; ++l ) {
        par = &mesh->info.par[l];
        if ( par->elt != MMG5_Triangle || par->ref==info ) continue;

        if ( tt1.ref!=par->ref && tt2.ref !=par->ref )  continue;

        hausd = MG_MIN(hausd,par->hausd);
        break;
      }
    }
  }

  /* Local params at tetra of the edge shell */
  if ( mesh->info.parTyp & MG_Tetra ) {
    l = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilist; ++k ) {
        pt = &mesh->tetra[list[k]/6];
        if ( par->ref != pt->ref ) continue;

        hausd   = par->hausd;
        isloc   = 1;
      }
    } while ( ++l<mesh->info.npar );

    for ( ; l<mesh->info.npar; ++l ) {
      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilist; ++k ) {
        pt = &mesh->tetra[list[k]/6];
        if ( par->ref != pt->ref ) continue;

        hausd = MG_MIN(hausd,par->hausd);
        break;
      }
    }
  }

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  MMG5_BezierEdge(mesh,np,nq,b0,b1,0,v);
  c[0] = b0[0] - (p0->c[0] + MMG5_ATHIRD*ux);
  c[1] = b0[1] - (p0->c[1] + MMG5_ATHIRD*uy);
  c[2] = b0[2] - (p0->c[2] + MMG5_ATHIRD*uz);

  disnat = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

  c[0] = b1[0] - (p1->c[0] - MMG5_ATHIRD*ux);
  c[1] = b1[1] - (p1->c[1] - MMG5_ATHIRD*uy);
  c[2] = b1[2] - (p1->c[2] - MMG5_ATHIRD*uz);

  disnat = MG_MAX(disnat, c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

  /* local param at vertices */
  // hausd = min (hausd_ref, hausd_np,hausd_nq)
  disnat = MG_MAX(disnat,hausd * hausd);

  p0 = &mesh->point[na1];
  p1 = &mesh->point[na2];
  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  /* local param at vertices */
  // hausd = min (hausd_ref, hausd_na1,hausd_na2)
  MMG5_BezierEdge(mesh,na1,na2,b0,b1,0,v);
  c[0] = b0[0] - (p0->c[0] + MMG5_ATHIRD*ux);
  c[1] = b0[1] - (p0->c[1] + MMG5_ATHIRD*uy);
  c[2] = b0[2] - (p0->c[2] + MMG5_ATHIRD*uz);

  dischg = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];

  c[0] = b1[0] - (p1->c[0] - MMG5_ATHIRD*ux);
  c[1] = b1[1] - (p1->c[1] - MMG5_ATHIRD*uy);
  c[2] = b1[2] - (p1->c[2] - MMG5_ATHIRD*uz);

  dischg = MG_MAX(dischg,c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  dischg = MG_MAX(dischg,hausd * hausd);

  if ( dischg > disnat )   return 0;

  if ( typchk==1 && met->size > 1 && met->m ) {
    cal1 = MMG5_caltri33_ani(mesh,met,&tt1);
    cal2 = MMG5_caltri33_ani(mesh,met,&tt2);
  }
  else {
    cal1 = MMG5_caltri(mesh,met,&tt1);
    cal2 = MMG5_caltri(mesh,met,&tt2);
  }

  calnat = MG_MIN(cal1,cal2);
  for (j=0; j<3; j++) {
    if ( tt1.v[j] == nq )  tt1.v[j] = na2;
    if ( tt2.v[j] == np )  tt2.v[j] = na1;
  }

  if ( typchk==1 && met->size > 1 && met->m ) {
    cal1 = MMG5_caltri33_ani(mesh,met,&tt1);
    cal2 = MMG5_caltri33_ani(mesh,met,&tt2);
  }
  else {
    cal1 = MMG5_caltri(mesh,met,&tt1);
    cal2 = MMG5_caltri(mesh,met,&tt2);
  }

  calchg = MG_MIN(cal1,cal2);
  if ( calchg < 1.01 * calnat )  return 0;

  /* Check mechanical validity of forthcoming operations */
  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];
  ppt0->c[0] = 0.5*(p0->c[0] + p1->c[0]);
  ppt0->c[1] = 0.5*(p0->c[1] + p1->c[1]);
  ppt0->c[2] = 0.5*(p0->c[2] + p1->c[2]);

#ifndef NDEBUG
  /* Security check: ensure that the edge is boundary */
  uint16_t  tag = 0;
  MMG5_int ref = 0;
  if ( !MMG3D_get_shellEdgeTag(mesh,list[0]/6,list[0]%6,&tag,&ref) ) {
    fprintf(stderr,"\n  ## Warning: %s: 0. unable to get edge info"
            " (tetra %d).\n",__func__,MMG3D_indElt(mesh,list[0]/6));
    return 0;
  }
  assert ( (tag & MG_BDY)  && "Edge should be boundary but is not");
#endif

  if ( met->m ) {
    pt  = &mesh->tetra[list[0]/6];
    assert ( pt->xt && "Boundary edge interpolated from non-boundary face");

    /* Mark edge as boundary to ensure suitable detection of bdy edge during
     * interpolation */
    mesh->xtetra[pt->xt].tag[list[0]%6] |= MG_BDY;

    if ( typchk == 1 && (met->size>1) ) {
      if ( MMG3D_intmet33_ani(mesh,met,list[0]/6,list[0]%6,0,0.5) <= 0 )
        return 0;
    }
    else {
      if ( MMG5_intmet(mesh,met,list[0]/6,list[0]%6,0,0.5) <= 0 )
        return 0;
    }
  }

  /* Check validity of insertion of midpoint on edge (pq), then collapse of m on a1 */
  calold = calnew = DBL_MAX;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    pt  = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    calold = MG_MIN(calold, pt->qual);
    assert ( isfinite(calold) );

    ia1 = ia2 = ip = iq = -1;
    for (j=0; j< 4; j++) {
      if (pt->v[j] == np)  ip = j;
      else if (pt->v[j] == nq) iq = j;
      else if ( ia1 < 0 ) ia1 = j;
      else ia2 = j;
    }
    assert((ip >= 0) && (iq >= 0) && (ia1 >= 0) && (ia2 >= 0));
    isshell = (pt->v[ia1] == na1 || pt->v[ia2] == na1);

    /* 2 elts resulting from split and collapse */
    pt0->v[ip] = 0;

    if ( typchk==1 && met->size > 1 && met->m )
      caltmp = MMG5_caltet33_ani(mesh,met,pt0);
    else
      caltmp = MMG5_orcal(mesh,met,0);


    if ( caltmp < MMG5_NULKAL )  return 0;

    if ( !isshell ) {
      /* Test that we don't recreate an existing elt */
      MMG5_int adj = mesh->adja[4*(iel-1)+1+ip];
      if ( adj ) {
        int8_t voy  = adj%4;
        adj /= 4;

        if ( mesh->tetra[adj].v[voy] == na1 ) {
          return 0;
        }
      }

      /* Test future quality */
      pt0->v[ip] = na1;

      if ( typchk==1 && met->size > 1 && met->m )
        caltmp = MMG5_caltet33_ani(mesh,met,pt0);
      else
        caltmp = MMG5_orcal(mesh,met,0);

      calnew = MG_MIN(calnew,caltmp);
    }
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[iq] = 0;

    if ( typchk==1 && met->size > 1 && met->m )
      caltmp = MMG5_caltet33_ani(mesh,met,pt0);
    else
      caltmp = MMG5_orcal(mesh,met,0);

    if ( caltmp < MMG5_NULKAL )  return 0;

    if ( !isshell ) {
      /* Test that we don't recreate an existing elt */
      MMG5_int adj = mesh->adja[4*(iel-1)+1+iq];
      if ( adj ) {
        int8_t voy  = adj%4;
        adj /= 4;

        if ( mesh->tetra[adj].v[voy] == na1 ) {
          return 0;
        }
      }

      /* Test future quality */
      pt0->v[iq] = na1;

      if ( typchk==1 && met->size > 1 && met->m )
        caltmp = MMG5_caltet33_ani(mesh,met,pt0);
      else
        caltmp = MMG5_orcal(mesh,met,0);

      if ( caltmp < MMG5_NULKAL )  return 0;

      calnew = MG_MIN(calnew,caltmp);
    }
  }
  if ( calold < MMG5_EPSOK && calnew <= calold ) return 0;

  else if ( calnew < 0.3 * calold )  return 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param met pointer to the solution structure
 * \param list pointer to the shell of the edge
 * \param ret dobble of the number of tetrahedra in the shell
 * \param it1 boundary face carrying the beforehand tested terminal
 * point for collapse
 * \param PROctree pointer to the PROctree structure in Delaunay mode,
 * NULL pointer in pattern mode.
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 * \return -1 if lack of memory, 0 if fail to swap, 1 otherwise
 *
 * Swap boundary edge whose shell is provided.
 *
 */
int MMG5_swpbdy(MMG5_pMesh mesh,MMG5_pSol met,int64_t *list,int ret,MMG5_int it1,
                 MMG3D_pPROctree PROctree, int8_t typchk) {
  MMG5_pTetra   pt,pt1;
  MMG5_pPoint   p0,p1;
  int           ilist;
  MMG5_int      iel,np,nq,nm,src,iel1;
  double        c[3];
  int8_t        ia,iface1,j,ipa,im;
  int           ier;
#ifndef NDEBUG
  MMG5_int      na;
#endif

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  np = pt->v[MMG5_iare[ia][0]];
  nq = pt->v[MMG5_iare[ia][1]];
#ifndef NDEBUG
  na = 0;
#endif

  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];

  /* search for na = the point on quadrangle surfacic configuration on which collapse
     validity has been checked in MMG5_chkswpbdy */
  iel1 = it1 / 4;
  iface1 = it1 % 4;
  pt1 = &mesh->tetra[iel1];

  for (j=0; j<3;j++) {
    ipa = MMG5_idir[iface1][j];
    if ( (pt1->v[ipa] != np)&&(pt1->v[ipa] != nq) ) {
#ifndef NDEBUG
      na = pt1->v[ipa];
#endif
      break;
    }
  }
  assert(na);

  /* Create midpoint m on edge (pq), then split edge */
  c[0] = 0.5*( p0->c[0] + p1->c[0]);
  c[1] = 0.5*( p0->c[1] + p1->c[1]);
  c[2] = 0.5*( p0->c[2] + p1->c[2]);
#ifdef USE_POINTMAP
  src = mesh->point[np].src;
#else
  src = 1;
#endif
  nm = MMG3D_newPt(mesh,c,MG_BDY,src);
  if ( !nm ) {
    MMG3D_POINT_REALLOC(mesh,met,nm,mesh->gap,
                         fprintf(stderr,"\n  ## Error: %s: unable to allocate a"
                                 " new point\n",__func__);
                         MMG5_INCREASE_MEM_MESSAGE();
                         return -1
                         ,c,MG_BDY,src);
  }
  assert ( met );
  if ( met->m ) {
    if ( typchk == 1 && (met->size>1) ) {
      if ( MMG3D_intmet33_ani(mesh,met,iel,ia,nm,0.5)<=0 )  return 0;
    }
    else {
      if ( MMG5_intmet(mesh,met,iel,ia,nm,0.5)<=0 )  return 0;
    }
  }

  ier = MMG5_split1b(mesh,met,list,ret,nm,0,typchk-1,0);
  /* pointer adress may change if we need to realloc memory during split */
  pt1 = &mesh->tetra[iel1];

  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to swap boundary edge.\n",
      __func__);
    return -1;
  }
  else if ( !ier )  {
    MMG3D_delPt(mesh,nm);
    return 0;
  }

  /* Collapse m on na after taking (new) ball of m */
  memset(list,0,(MMG3D_LMAX+2)*sizeof(MMG5_int));
  for (j=0; j<3; j++) {
    im = MMG5_idir[iface1][j];
    if ( pt1->v[im] == nm )  break;
  }
  if ( pt1->v[im] != nm ){
    MMG3D_delPt(mesh,nm);
    fprintf(stderr,"\n  # Warning: %s: pt1->v[im] != nm.\n",__func__);
    return 0;
  }
  ilist = MMG5_boulevolp(mesh,iel1,im,list);

  assert(list[0]/4 == iel1);
  assert(pt1->v[ipa] == na);

  ier = MMG5_colver(mesh,met,list,ilist,ipa,typchk);
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to swap boundary edge.\n",
      __func__);
    return -1;
  }
  else if ( ier ) {
    MMG3D_delPt(mesh,ier);
    ier = 1;
  }

  /* Check for non convex situation */
  assert ( ier && "Unable to collapse the point created during the boundary swap");

  return ier;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param k index of the tetrahedron with multiple boundary faces (to be swapped).
 * \param metRidTyp metric storage (classic or special)
 * \param ifac face of the tetra \a k that give the best results for the swap23
 * \param conf0 detected configuration for the swap23 of the tetra \a k
 * \param adj neighbour of the tetra k through the face \a ifac (4*k1+ifac1)
 * \param conf1 detected configuration for the swap23 of the tetra \a adj/4
 * \return -1 if lack of memory, 0 if fail to swap, 1 otherwise.
 *
 * Search an adjacent to the tetra \a k and perform swap 2->3 (the common face
 * of the 2 tetra is destroyed and replaced by a common edge used by the three
 * new elts).
 *
 * \remark used in anatet4 to remove the tetra with multiple boundary faces.
 *
 */
int MMG3D_swap23(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t metRidTyp,
                 int ifac,int conf0,MMG5_int adj,int conf1) {
  MMG5_pTetra   pt0,pt1,ptnew;
  MMG5_xTetra   xt[3];
  MMG5_pxTetra  pxt0,pxt1;
  MMG5_int      xt1,k1,*adja,iel,np;
  MMG5_int      adj0_2,adj0_3,adj1_1,adj1_2,adj1_3;
  int8_t        i,isxt[3];
  uint8_t       tau0[4],tau1[4];
  const uint8_t *taued0,*taued1;

  pt0     = &mesh->tetra[k];

  assert ( pt0->xt );
  assert ( ifac>=0 && adj>0 );

  /** Neighbouring element with which we will try to swap */
  k1   = adj/4;

  assert(k1);

  /* Search in which configurations are the tetrahedra (default is case 0-0)
   *
   *           3                    2------------- 0
   *         ,/|`\                  |`\          /|
   *       ,/  |  `\                |  `\       /.|
   *     ,/    '.   `\              '.   `\    / |
   *   ,/       |     `\             |     `\ / .|
   * ,/         |       `\           |       /\.|
   * 0-----------'.--------2          '     /  3
   * `\.         |      ,/            |    / ,/
   *    `\.      |    ,/              |   /,/
   *       `\.   '. ,/                '. ,/
   *          `\. |/                   |/
   *             `1                    1
   */

  /* k may be in configuration 0, 3, 6 or 9. Default is case 0 */
  switch(conf0) {
  case 3:
    tau0[0] = 1; tau0[1] = 0; tau0[2] = 3; tau0[3] = 2;
    taued0 = &MMG5_permedge[3][0];
    break;
  case 6:
    tau0[0] = 2; tau0[1] = 0; tau0[2] = 1; tau0[3] = 3;
    taued0 = &MMG5_permedge[6][0];
    break;
  case 9:
    tau0[0] = 3; tau0[1] = 0; tau0[2] = 2; tau0[3] = 1;
    taued0 = &MMG5_permedge[9][0];
    break;
  default:
    assert ( !conf0 );

    tau0[0] = 0; tau0[1] = 1; tau0[2] = 2; tau0[3] = 3;
    taued0 = &MMG5_permedge[0][0];
    break;
  }

  /* k1 may be in configuration adj%4, adj%4+1, adj%4+2. Default case is case 0 */
  pt1 = &mesh->tetra[k1];

  assert(pt0->ref == pt1->ref);

  switch(conf1) {
  case 1:
    tau1[0] = 0; tau1[1] = 2; tau1[2] = 3; tau1[3] = 1;
    taued1 = &MMG5_permedge[1][0];
    break;
  case 2:
    tau1[0] = 0; tau1[1] = 3; tau1[2] = 1; tau1[3] = 2;
    taued1 = &MMG5_permedge[2][0];
    break;
  case 3:
    tau1[0] = 1; tau1[1] = 0; tau1[2] = 3; tau1[3] = 2;
    taued1 = &MMG5_permedge[3][0];
    break;
  case 4:
    tau1[0] = 1; tau1[1] = 3; tau1[2] = 2; tau1[3] = 0;
    taued1 = &MMG5_permedge[5][0];
    break;
  case 5:
    tau1[0] = 1; tau1[1] = 2; tau1[2] = 0; tau1[3] = 3;
    taued1 = &MMG5_permedge[4][0];
    break;
  case 6:
    tau1[0] = 2; tau1[1] = 0; tau1[2] = 1; tau1[3] = 3;
    taued1 = &MMG5_permedge[6][0];
    break;
  case 7:
    tau1[0] = 2; tau1[1] = 1; tau1[2] = 3; tau1[3] = 0;
    taued1 = &MMG5_permedge[7][0];
    break;
  case 8:
    tau1[0] = 2; tau1[1] = 3; tau1[2] = 0; tau1[3] = 1;
    taued1 = &MMG5_permedge[8][0];
    break;
  case 9:
    tau1[0] = 3; tau1[1] = 0; tau1[2] = 2; tau1[3] = 1;
    taued1 = &MMG5_permedge[9][0];
    break;
  case 10:
    tau1[0] = 3; tau1[1] = 2; tau1[2] = 1; tau1[3] = 0;
    taued1 = &MMG5_permedge[11][0];
    break;
  case 11:
    tau1[0] = 3; tau1[1] = 1; tau1[2] = 0; tau1[3] = 2;
    taued1 = &MMG5_permedge[10][0];
    break;
  default:
    assert(!conf1);
    tau1[0] = 0; tau1[1] = 1; tau1[2] = 2; tau1[3] = 3;
    taued1 = &MMG5_permedge[0][0];
    break;
  }

  /** Swap */
  /* Store useful information from pt1 before overwrite by memcpy*/
  xt1 = pt1->xt;

  np    = pt1->v[tau1[0]];

  MMG5_int ref[6] = {0};
  uint16_t  tag[6] = {0};
  for (i=0;i<6;i++) {
    if ( !MMG3D_get_shellEdgeTag(mesh,k1,taued1[i],&tag[i],&ref[i]) ) {
      fprintf(stderr,"\n  ## Error: %s: %d. unable to get edge info.\n",__func__,i);
      return 0;
    }
  }

  memcpy(pt1,pt0,sizeof(MMG5_Tetra));

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return -1);
    pt0 = &mesh->tetra[k];
    pt1 = &mesh->tetra[k1];
  }
  ptnew = &mesh->tetra[iel];
  memcpy(ptnew,pt0,sizeof(MMG5_Tetra));

  /* First tetra: k */
  pt0->v[tau0[1]] = np;

  /* Second tetra: k1 */
  pt1->v[tau0[2]] = np;

  /* Third tetra: iel */
  ptnew->v[tau0[3]] = np;

  /* xtetra and adjacency update */
  pxt0 = &mesh->xtetra[pt0->xt];
  memcpy(&xt[0],pxt0,sizeof(MMG5_xTetra));
  memcpy(&xt[1],pxt0,sizeof(MMG5_xTetra));
  memcpy(&xt[2],pxt0,sizeof(MMG5_xTetra));

  /* Store the old adja */
  adja = &mesh->adja[4*(k-1) +1];
  adj0_2 = adja[tau0[2]];
  adj0_3 = adja[tau0[3]];

  adja = &mesh->adja[4*(k1-1) +1];
  adj1_1 = adja[tau1[1]];
  adj1_2 = adja[tau1[2]];
  adj1_3 = adja[tau1[3]];

  /* New adja for the new tets */
  adja = &mesh->adja[4*(k-1) +1];
  adja[tau0[0]] = adj1_1;
  adja[tau0[2]] = 4*k1  + tau0[1] ;
  adja[tau0[3]] = 4*iel + tau0[1] ;
  if ( adj1_1 )
    mesh->adja[4*(adj1_1/4-1) + 1 + adj1_1%4] = 4*k + tau0[0];

  adja = &mesh->adja[4*(k1-1) +1];
  adja[tau0[0]] = adj1_3;
  adja[tau0[1]] = 4*k   + tau0[2] ;
  adja[tau0[2]] = adj0_2;
  adja[tau0[3]] = 4*iel + tau0[2] ;
  if ( adj1_3 )
    mesh->adja[4*(adj1_3/4-1) + 1 + adj1_3%4] = 4*k1 + tau0[0];
  if ( adj0_2 )
    mesh->adja[4*(adj0_2/4-1) + 1 + adj0_2%4] = 4*k1 + tau0[2];

  adja = &mesh->adja[4*(iel-1) +1];
  adja[tau0[0]] = adj1_2;
  adja[tau0[1]] = 4*k   + tau0[3] ;
  adja[tau0[2]] = 4*k1  + tau0[3] ;
  adja[tau0[3]] = adj0_3;
  if ( adj1_2 )
    mesh->adja[4*(adj1_2/4-1) + 1 + adj1_2%4] = 4*iel + tau0[0];
  if ( adj0_3 )
    mesh->adja[4*(adj0_3/4-1) + 1 + adj0_3%4] = 4*iel + tau0[3];

  pxt1 = NULL;
  if ( !pt1->xt ) {
    /* Assignation of the xt fields to the appropriate tets */
    /* xt[0] */
    xt[0].tag[taued0[0]] = 0;
    xt[0].tag[taued0[3]] = 0;
    xt[0].tag[taued0[4]] = 0;

    xt[0].edg[taued0[0]] = 0;
    xt[0].edg[taued0[3]] = 0;
    xt[0].edg[taued0[4]] = 0;

    xt[0].ref[ tau0[0]] = 0;
    xt[0].ref[ tau0[2]] = 0;
    xt[0].ref[ tau0[3]] = 0;
    xt[0].ftag[tau0[0]] = 0;
    xt[0].ftag[tau0[2]] = 0;
    xt[0].ftag[tau0[3]] = 0;

    MG_SET(xt[0].ori, tau0[0]);
    MG_SET(xt[0].ori, tau0[2]);
    MG_SET(xt[0].ori, tau0[3]);

    /* xt[1] */
    xt[1].tag[taued0[1]] = 0;
    xt[1].tag[taued0[3]] = 0;
    xt[1].tag[taued0[5]] = 0;

    xt[1].edg[taued0[1]] = 0;
    xt[1].edg[taued0[3]] = 0;
    xt[1].edg[taued0[5]] = 0;

    xt[1].ref[ tau0[0]] = 0;
    xt[1].ref[ tau0[1]] = 0;
    xt[1].ref[ tau0[3]] = 0;
    xt[1].ftag[tau0[0]] = 0;
    xt[1].ftag[tau0[1]] = 0;
    xt[1].ftag[tau0[3]] = 0;

    MG_SET(xt[1].ori, tau0[0]);
    MG_SET(xt[1].ori, tau0[1]);
    MG_SET(xt[1].ori, tau0[3]);

    /* xt[2] */
    xt[1].tag[taued0[2]] = 0;
    xt[1].tag[taued0[4]] = 0;
    xt[1].tag[taued0[5]] = 0;

    xt[1].edg[taued0[2]] = 0;
    xt[1].edg[taued0[4]] = 0;
    xt[1].edg[taued0[5]] = 0;

    xt[1].ref[ tau0[0]] = 0;
    xt[1].ref[ tau0[1]] = 0;
    xt[1].ref[ tau0[2]] = 0;
    xt[1].ftag[tau0[0]] = 0;
    xt[1].ftag[tau0[1]] = 0;
    xt[1].ftag[tau0[2]] = 0;

    MG_SET(xt[1].ori, tau0[0]);
    MG_SET(xt[1].ori, tau0[1]);
    MG_SET(xt[1].ori, tau0[2]);

  }
  else {
    pxt1 = &mesh->xtetra[xt1];

    /* Assignation of the xt fields to the appropriate tets */
    /* Warning: after collapses, some boundary edges not connected to boundary
     * faces may have a 0 tag inside a xtetra (see \ref MMG5_colver when a
     * xtetra is assigned to one of the neighbours of the tetra of the edge
     * shell). In consequence, we cannot simply use the stored tags. */

    /* xt[0] */
    xt[0].tag[taued0[0]] = 0;

    xt[0].tag[taued0[3]] = tag[2];
    xt[0].tag[taued0[4]] = tag[1];
    /* As the edge tag of tetra 0 may be erroneous if the edge doesn't belong to
     * a boundary face */
    xt[0].tag[taued0[5]] = tag[5];


    xt[0].edg[taued0[0]] = 0;
    xt[0].edg[taued0[3]] = ref[2];
    xt[0].edg[taued0[4]] = ref[1];
    xt[0].edg[taued0[5]] = ref[5];

    xt[0].ref[ tau0[0]] = pxt1->ref[tau1[1]];
    xt[0].ref[ tau0[2]] = 0;
    xt[0].ref[ tau0[3]] = 0;
    xt[0].ftag[tau0[0]] = pxt1->ftag[tau1[1]];
    xt[0].ftag[tau0[2]] = 0;
    xt[0].ftag[tau0[3]] = 0;

    if ( MG_GET(pxt1->ori,tau1[1]) ) MG_SET(xt[0].ori, tau0[0]);
    MG_SET(xt[0].ori, tau0[2]);
    MG_SET(xt[0].ori, tau0[3]);

    /* xt[1] */
    xt[1].tag[taued0[1]] = 0;

    xt[1].tag[taued0[3]] = tag[0];
    xt[1].tag[taued0[5]] = tag[1];

    /* As the edge tag of tetra 0 may be erroneous if the edge doesn't belong to
     * a boundary face */
    xt[1].tag[taued0[4]] = tag[3];

    xt[1].edg[taued0[1]] = 0;
    xt[1].edg[taued0[3]] = ref[0];
    xt[1].edg[taued0[5]] = ref[1];
    xt[1].edg[taued0[4]] = ref[3];

    xt[1].ref[ tau0[0]] = pxt1->ref[tau1[3]];
    xt[1].ref[ tau0[1]] = 0;
    xt[1].ref[ tau0[3]] = 0;
    xt[1].ftag[tau0[0]] = pxt1->ftag[tau1[3]];
    xt[1].ftag[tau0[1]] = 0;
    xt[1].ftag[tau0[3]] = 0;

    if ( MG_GET(pxt1->ori,tau1[3]) ) MG_SET(xt[1].ori, tau0[0]);
    MG_SET(xt[1].ori, tau0[1]);
    MG_SET(xt[1].ori, tau0[3]);

    /* xt[2] */
    xt[2].tag[taued0[2]] = 0;

    xt[2].tag[taued0[4]] = tag[0];
    xt[2].tag[taued0[5]] = tag[2];

    /* As the edge tag of tetra 0 may be erroneous if the edge doesn't belong to
     * a boundary face */
    xt[2].tag[taued0[3]] = tag[4];

    xt[2].edg[taued0[2]] = 0;
    xt[2].edg[taued0[4]] = ref[0];
    xt[2].edg[taued0[5]] = ref[2];
    xt[2].edg[taued0[3]] = ref[4];

    xt[2].ref[ tau0[0]] = pxt1->ref[tau1[2]];
    xt[2].ref[ tau0[1]] = 0;
    xt[2].ref[ tau0[2]] = 0;
    xt[2].ftag[tau0[0]] = pxt1->ftag[tau1[2]];
    xt[2].ftag[tau0[1]] = 0;
    xt[2].ftag[tau0[2]] = 0;

    if ( MG_GET(pxt1->ori,tau1[2]) ) MG_SET(xt[2].ori, tau0[0]);
    MG_SET(xt[2].ori, tau0[1]);
    MG_SET(xt[2].ori, tau0[2]);
  }

    /* Assignation of the xt fields to the appropriate tets */
  isxt[0] = isxt[1] = isxt[2] = 0;
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
  }

  assert ( (isxt[0] && isxt[1]) || (isxt[1] && isxt[2]) || (isxt[2] && isxt[0]) );

  if ( isxt[0] ) {
    memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));

    if ( xt1 ) {
      if ( isxt[1] ) {
        pt1->xt = xt1;
        memcpy(pxt1,&xt[1],sizeof(MMG5_xTetra));
        if ( isxt[2] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return -1);
          }
          ptnew->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[2],sizeof(MMG5_xTetra));
        }
        else ptnew->xt = 0;
      }
      else {
        pt1->xt   = 0;
        ptnew->xt = xt1;
        memcpy(pxt1,&xt[2],sizeof(MMG5_xTetra));
      }
    }
    else {
      if ( isxt[1] ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return -1);
        }
        pt1->xt = mesh->xt;
        pxt0 = &mesh->xtetra[mesh->xt];
        memcpy(pxt0,&xt[1],sizeof(MMG5_xTetra));
      }
      else pt1->xt = 0;

      if ( isxt[2] ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return -1);
        }
        ptnew->xt = mesh->xt;
        pxt0 = &mesh->xtetra[mesh->xt];
        memcpy(pxt0,&xt[2],sizeof(MMG5_xTetra));
      }
      else ptnew->xt = 0;
    }
  }
  else {
    pt0->xt = 0;
    memcpy(pxt0 ,&xt[2],sizeof(MMG5_xTetra));

    if ( xt1 ) {
      pt1->xt = xt1;
      memcpy(pxt1,&xt[1],sizeof(MMG5_xTetra));
    }
    else {
      mesh->xt++;
      if ( mesh->xt > mesh->xtmax ) {
        /* realloc of xtetras table */
        MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                           "larger xtetra table",
                           mesh->xt--;
                           fprintf(stderr,"  Exit program.\n");
                           return -1);
      }
      pt1->xt = mesh->xt;
      pxt0 = &mesh->xtetra[mesh->xt];
      memcpy(pxt0,&xt[1],sizeof(MMG5_xTetra));
    }
  }

  /** Quality Update */
  if ( (!metRidTyp) && met->m && met->size>1 )
    pt0->qual   = MMG5_caltet33_ani(mesh,met,pt0);
  else
    pt0->qual   = MMG5_orcal(mesh,met,k);

  if ( (!metRidTyp) && met->m && met->size>1 )
    pt1->qual   = MMG5_caltet33_ani(mesh,met,pt1);
  else
    pt1->qual   = MMG5_orcal(mesh,met,k1);

  if ( (!metRidTyp) && met->m && met->size>1 )
    ptnew->qual   = MMG5_caltet33_ani(mesh,met,ptnew);
  else
    ptnew->qual   = MMG5_orcal(mesh,met,iel);

  pt0->mark   = mesh->mark;
  pt1->mark   = mesh->mark;
  ptnew->mark = mesh->mark;

  return 1;
}
