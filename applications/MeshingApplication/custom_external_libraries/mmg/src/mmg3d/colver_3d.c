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
 * \file mmg3d/colver_3d.c
 * \brief Functions for vertices collapsing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "inlined_functions_3d_private.h"
#include "mmg3dexterns_private.h"

extern int8_t  ddb;

/** Check whether collapse ip -> iq could be performed, ip internal ;
 *  'mechanical' tests (positive jacobian) are not performed here */
int MMG5_chkcol_int(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t iface,
                    int8_t iedg,int64_t *list,int ilist,int8_t typchk) {
  MMG5_pTetra   pt,pt0;
  MMG5_pPoint   p0;
  double        calold,calnew,caltmp,ll,lon;
  MMG5_int      j,iel,nq,nr;
  int8_t        i,jj,ip,iq;

  iq  = MMG5_idir[iface][MMG5_iprv2[iedg]];
  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  nq  = pt->v[iq];

  lon = 1.6;
  if ( typchk == 2 && met->m ) {
    lon = MMG5_lenedg(mesh,met,MMG5_iarf[iface][iedg],pt);

    if ( !lon ) return 0;
    /*on cherche a se rapprocher de 1*/
    //lon = MG_MAX(0.8/lon,1.6);// test ok but less good than the next one
    lon = MG_MAX(2.-lon,1.6);
  }

  calold = calnew = DBL_MAX;
  for (j=0; j<ilist; j++) {
    iel = list[j] / 4;
    ip  = list[j] % 4;
    pt  = &mesh->tetra[iel];

    /* exclude elements from shell */
    for (jj=0; jj<4; jj++)  if ( pt->v[jj] == nq )  break;
    if ( jj < 4 )  continue;
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    /* Update edges tag for pt0 (needed by lenedg). */

    /* prevent from recreating internal edge between boundaries */
    if ( mesh->info.fem==typchk ) {
      p0 = &mesh->point[nq];
      if ( (p0->tag & MG_BDY) && !(p0->tag & MG_PARBDY) ) {
        i = ip;
        for (jj=0; jj<3; jj++) {
          i = MMG5_inxt3[i];
          p0 = &mesh->point[pt->v[i]];
          if ( (p0->tag & MG_BDY) && !(p0->tag & MG_PARBDY) ){
            return 0;
          }
        }

        /* Prevent from creating a tetra with 4 bdy vertices */
        // Algiane (2022) this test is useless I think (because we forbid the
        // creation of internal edges between boundary points)
#ifndef NDEBUG
        i  = ip;
        nr = 0;
        for (jj=0; jj<3; jj++) {
          i = MMG5_inxt3[i];
          p0 = &mesh->point[pt->v[i]];
          if ( (p0->tag & MG_BDY) && !(p0->tag & MG_PARBDY) ) ++nr;
        }
        if ( nr==3 ) {
          assert ( 0 && "Uncomment this test, it is not useless");
          return 0;
        }
#endif
      }
    }
    else {
      /* In aniso : prevent from creating a tetra with 4 ridges vertices or
       * internal edges between two ridges (unable to split it after because of
       * the ridge metric) */
      if ( met->size==6 ) {
        p0 = &mesh->point[nq];

        if ( (p0->tag & MG_GEO) && !MG_SIN_OR_NOM(p0->tag) ) {
          i = ip;
          for (jj=0; jj<3; jj++) {
            i = MMG5_inxt3[i];
            p0 = &mesh->point[pt->v[i]];
            if ( p0->tag & MG_GEO && !MG_SIN_OR_NOM(p0->tag) ) {
              return 0;
            }
          }

          // Algiane (2022) this test is useless because we forbid the creation of
          // internal edges between boundary points but we can keep it in case we
          // comment the previous test
#ifndef NDEBUG
          i  = ip;
          nr = 0;
          for (jj=0; jj<3; jj++) {
            i = MMG5_inxt3[i];
            p0 = &mesh->point[pt->v[i]];
            if ( p0->tag & MG_GEO && !MG_SIN_OR_NOM(p0->tag) ) {
              ++nr;
            }
          }
          if ( nr==3 ) {
            assert ( 0 && "Uncomment this test, it is not useless");
            return 0;
          }
#endif
        }
      }
    }

    pt0->v[ip] = nq;

    calold = MG_MIN(calold,pt->qual);
    if ( typchk==1 && met->m && met->size > 1 )
      caltmp = MMG5_caltet33_ani(mesh,met,pt0);
    else
      caltmp = MMG5_orcal(mesh,met,0);

    if ( caltmp < MMG5_NULKAL )  return 0;
    calnew = MG_MIN(calnew,caltmp);
    /* check length */
    if ( typchk == 2 && met->m ) {
      for (jj=0; jj<6; jj++) {
        /* Rough evaluation of edge length (doesn't take into account if some of
         * the modified edges of pt0 are boundaries): for a more precise
         * computation, we need to update the edge tags of pt0.  */
        ll = MMG5_lenedgspl(mesh,met,jj,pt0);
        if ( (!ll) || (ll > lon) )//LOPTL too small, we need to put greater than 1.41
          return 0;
      }
    }
  }
  if ( calold < MMG5_EPSOK && calnew <= calold )  return 0;
  else if ( calnew < MMG5_EPSOK || calnew < 0.3*calold )  return 0;

  return ilist;
}

/**
 * \param mesh pointer toward the mesh
 * \param start tetra from which we start to travel
 * \param end tetra at which we stop the travel
 * \param na edge vertex
 * \param nb edge vertex
 * \param piv global index of the pivot to set the sense of travel
 * \param iel pointer toward the last element of the shell
 * \param iopp pointer toward the ending boundary face of the shell
 *
 * \return -1 if fail, \a piv otherwise.
 *
 * Unfold the shell of the edge \a na \a nb from tetra \a start in the direction
 * given by the pivot \a piv).
 *
 */
static inline
MMG5_int MMG3D_unfold_shell(MMG5_pMesh  mesh,MMG5_int start,MMG5_int end, MMG5_int na, MMG5_int nb,MMG5_int piv,
                       MMG5_int *iel,int8_t *iopp) {
  MMG5_pTetra  pt;
  MMG5_int     adj,*adja;
  int8_t       i,ipiv,isface;

  adj = start;
  do {
    *iel  = adj;
    pt   = &mesh->tetra[*iel];
    adja = &mesh->adja[4*(*iel-1)+1];

    /* Identification of edge number in tetra iel */
    if ( !MMG3D_findEdge(mesh,pt,*iel,na,nb,1,NULL,&i) ) return -1;

    /* set sense of travel */
    if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
      adj  = adja[ MMG5_ifar[i][0] ] / 4;
      ipiv = MMG5_ifar[i][1];
      *iopp = MMG5_ifar[i][0];
      piv  = pt->v[ipiv];
    }
    else {
      adj  = adja[ MMG5_ifar[i][1] ] / 4;
      ipiv = MMG5_ifar[i][0];
      *iopp = MMG5_ifar[i][1];
      piv  = pt->v[ipiv];
    }

    isface = 0;
    if ( pt->xt ) {
      isface = (MG_BDY & mesh->xtetra[pt->xt].ftag[*iopp]);
    }
  }
  while ( adj && ( adj != end ) && !isface );

  return piv;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param k index of the starting tetra.
 * \param iface local index of the starting face in the tetra \a k.
 * \param ideg local index of the starting edge in the face \a iface.
 * \param lists surfacic ball of p.
 * \param ilists number of elements in the surfacic ball of p.
 *
 * \return 0 if the check of the topology fails, 1 if success, -1 if the
 * function fails.
 *
 * Topological check on the surface ball of np and nq in collapsing np->nq ;
 *  iface = boundary face on which lie edge iedg - in local face num.  (pq, or
 *  ia in local tet notation). See the Mmg Google
 *  Drive/Documentation/mmg3d/topchkcol_bdy3D.pdf for a picture of the
 *  configuration.
 *
 */
static int
MMG5_topchkcol_bdy(MMG5_pMesh mesh,MMG5_int k,int iface,int8_t iedg,MMG5_int *lists,
                    int ilists) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  double        n0[3],n1[3],devnew;
  MMG5_int      nump,numq,iel,jel,jel1,nap,nbp,naq,nbq,nro;
  int8_t        ip,iq,iopp,i,j,j1,jface,jface1;

  pt = &mesh->tetra[k];
  ip = MMG5_idir[iface][MMG5_inxt2[iedg]];
  iq = MMG5_idir[iface][MMG5_iprv2[iedg]];
  nump = pt->v[ip];
  numq = pt->v[iq];

  /* Surface ball has been enumerated as f1,...,f2 - f1,f2 = both triangles of
   * surface shell, point nap, facing the first vanishing face in surface ball
   * of p */
  nro = pt->v[MMG5_idir[iface][iedg]];

  jel = lists[1] / 4;
  jface = lists[1] % 4;

  pt = &mesh->tetra[jel];
  for (j=0; j<3; j++) {
    i = MMG5_idir[jface][j];
    if ( pt->v[i] != nump && pt->v[i] != nro ) break;
  }
  assert(j<3);

  nap = pt->v[i];

  /* Unfold shell of (numq,nro), starting from (jel,jface), with pivot nump */
  naq = MMG3D_unfold_shell(mesh,k,k,numq,nro,nump,&iel,&iopp);

  if ( naq < 0 ) {
    return -1;
  }
  else if ( nap == naq ) {
    return 0;
  }

  assert ( mesh->tetra[k].xt && "initial tetra is not boundary");
  pxt = &mesh->xtetra[mesh->tetra[k].xt];

  if ( !( MG_GEO_OR_NOM(pxt->tag[MMG5_iarf[iface][MMG5_inxt2[iedg]]]) ||
          MG_GEO_OR_NOM(pxt->tag[MMG5_iarf[iface][MMG5_iprv2[iedg]]])   ) ) {

    /* Check the normal deviation between the boundary faces sharing the edge
     * numq (or nump)-nro */
    if ( !MMG5_norpts (mesh,numq,nro,nap,n0) )  return 0;
    if ( !MMG5_norface(mesh,iel,iopp    ,n1) )  return 0;

    devnew = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
    if ( devnew < mesh->info.dhd ) return 0;
  }

  /*  Point nbp, facing the second vanishing face in surface ball of p */
  jel1   = lists[ilists-1] / 4;
  jface1 = lists[ilists-1] % 4;
  pt     = &mesh->tetra[jel1];
  for (j1=0; j1<3; j1++) {
    i = MMG5_idir[jface1][j1];
    if ( pt->v[i] != nump && pt->v[i] != numq )  break;
  }
  assert(j1<3);

  nro   = pt->v[i];
  jel   = lists[ilists-2] / 4;
  jface = lists[ilists-2] % 4;
  pt    = &mesh->tetra[jel];
  for (j=0; j<3; j++) {
    i = MMG5_idir[jface][j];
    if ( pt->v[i] != nump && pt->v[i] != nro )  break;
  }
  assert(j<3);

  nbp = pt->v[i];

  /* Unfold shell of (numq,nro), starting from (jel,jface), with pivot nump */
  nbq = MMG3D_unfold_shell(mesh,lists[ilists-1]/4,k,numq,nro,nump,&iel,&iopp);

  if ( nbq < 0 ) {
    return -1;
  }
  else if ( nbp == nbq ) {
    return 0;
  }

  pxt      = &mesh->xtetra[mesh->tetra[jel1].xt];

  if ( !( MG_GEO_OR_NOM(pxt->tag[MMG5_iarf[jface1][MMG5_iprv2[j1]]]) ||
          MG_GEO_OR_NOM(pxt->tag[MMG5_iarf[jface1][MMG5_inxt2[j1]]])   ) ) {

    /* Check the normal deviation between the boundary faces sharing the edge
     * numq (or nump)-nro */
    if ( !MMG5_norpts (mesh,nro,numq,nbp,n0) )  return 0;
    if ( !MMG5_norface(mesh,iel,iopp    ,n1) )  return 0;

    devnew = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
    if ( devnew < mesh->info.dhd )  return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param start tetra from which we start to travel
 * \param na edge vertex
 * \param nb edge vertex
 * \param tag new edge tag
 * \param ref new edge ref
 * \param piv global index of the pivot to set the sense of travel
 * \param adj index of adjacent tetra for the travel
 * \param filled 1 if an xtetra has been found (so tag and ref are filled)
 *
 * \return -1 if fail, \a start if shell has been completely travelled without
 * founding an xtetra, adj if an xtetra has been found, 0 otherwise
 *
 * Get tag and ref of the edge \a na \a nb from tetra \a start by traveling its
 * shell in one direction (given by the pivot \a piv).  Stop when meeting the
 * first xtetra with a non 0 tag (it is sufficient if tags and refs are
 * consistent through the edge shell except for edges with null tags);
 *
 */
static inline
int MMG3D_get_shellEdgeTag_oneDir(MMG5_pMesh  mesh,MMG5_int start, MMG5_int na, MMG5_int nb,
                                  int16_t *tag,MMG5_int *ref, MMG5_int piv,MMG5_int adj,
                                  int8_t *filled) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     *adja;
  int8_t       i;

  *filled = 0;
  while ( adj && (adj != start) ) {
    pt     = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,1,NULL,&i) ) return -1;

    /* update edge ref and tag */
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      *ref  = pxt->edg[i];
      if ( pxt->tag[i] & MG_BDY ) {
        *tag |= pxt->tag[i];
        *filled = 1;
        return adj;
      }
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ MMG5_ifar[i][1] ] /4;
      piv = pt->v[ MMG5_ifar[i][0] ];
    }
  }

  return adj;
}

/**
 * \param mesh pointer toward the mesh
 * \param start tetra from which we start to travel
 * \param ia local index of edge that must be updated
 * \param tag new edge tag
 * \param ref new edge ref
 * \return 1 if success, 0 if fail.
 *
 * Get tag and ref of the edge \ia of tetra \a start by traveling its shell.
 * Stop when meeting the first xtetra (it is sufficient if tags and refs are
 * consistent through the edge shell);
 *
 */
static inline
int MMG3D_get_shellEdgeTag(MMG5_pMesh  mesh,MMG5_int start, int8_t ia,int16_t *tag,MMG5_int *ref) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     piv,na,nb,adj,*adja;
  int8_t       filled;

  pt   = &mesh->tetra[start];

  assert( start >= 1 &&  MG_EOK(pt) );

  pxt  = NULL;
  na   = pt->v[MMG5_iare[ia][0]];
  nb   = pt->v[MMG5_iare[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( pxt->tag[ia] & MG_BDY ) {
      *tag |= pxt->tag[ia];
      *ref = pxt->edg[ia];
      return 1;
    }
  }

  /* Travel in one direction */
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][0]] / 4;
  piv = pt->v[MMG5_ifar[ia][1]];

  adj = MMG3D_get_shellEdgeTag_oneDir(mesh,start,na,nb,tag,ref,piv,adj,&filled);

  /* If an xtetra has been found or if all shell has been travelled, stop, else,
   * travel it the other sense */
  if ( adj > 0 ) {
    assert ( (adj == start) || filled );
    return 1;
  }
  else if ( adj < 0 ) return 0;

  assert(!adj);

  pt = &mesh->tetra[start];
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][1]] / 4;
  piv = pt->v[MMG5_ifar[ia][0]];

  adj = MMG3D_get_shellEdgeTag_oneDir(mesh,start,na,nb,tag,ref,piv,adj,&filled);

  if ( adj < 0 ) return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element in which we collapse.
 * \param iface face through wich we perform the collapse
 * \param iedg edge to collapse (in local face num)
 * \param listv pointer toward the list of the tetra in the ball of \a p0.
 * \param ilistv number of tetra in the ball of \a p0.
 * \param lists pointer toward the surfacic ball of \a p0.
 * \param ilists number of tetra in the surfacic ball of \a p0.
 * \param refmin reference of one of the two subdomains in presence
 * \param refplus reference of the other subdomain in presence
 * \param typchk  typchk type of checking permformed for edge length
 * (hmax or MMG3D_LLONG criterion).
 * \param isnm 1 if edge is non-manifold
 * \param isnmint 1 if ip is an internal non manifold point;
 *
 * \return ilistv if success, 0 if the point cannot be collapsed, -1 if fail.
 *
 * Check whether collapse ip -> iq could be performed, ip boundary point;
 *  'mechanical' tests (positive jacobian) are not performed here ;
 *  iface = boundary face on which lie edge iedg - in local face num.
 *  (pq, or ia in local tet notation).
 * If isnm is 1, the collapse occurs along an external MG_NOM edge.
 * If isnmint is 1, ip is an internal non manifold point and dont have normals.
 *  In this last case, \a lists, \a ilists \a refmin, \a refplus and \a isnm
 *  variables aren't used (we neither have a surfacic ball nor "positive" and
 *  "negative" volumes)
 *
 * \remark we don't check edge lengths.
 */
int MMG5_chkcol_bdy(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t iface,
                    int8_t iedg,int64_t *listv,int ilistv,MMG5_int *lists,int ilists,
                    MMG5_int refmin,MMG5_int refplus, int8_t typchk,int isnm,int8_t isnmint) {
  MMG5_pTetra  pt,pt0,pt1;
  MMG5_pxTetra pxt;
  MMG5_Tria    tt;
  MMG5_pPar    par;
  double       calold,calnew,caltmp,nadja[3],nprvold[3],nprvnew[3],ncurold[3],ncurnew[3];
  double       ps,devold,devnew,hmax,hausd;
  MMG5_int     nump,numq,ndepmin,ndepplus,l,kk,iel;
  int          nr,nbbdy,isloc,iedgeOpp,ipp;
  int16_t      tag;
  int8_t       iopp,iopp2,ia,ip,i,iq,i0,i1,ier,isminp,isplp;
#ifndef NDEBUG
  MMG5_pPoint  p0;
#endif

  pt   = &mesh->tetra[k];
  pxt  = 0;
  pt0  = &mesh->tetra[0];
  ip   = MMG5_idir[iface][MMG5_inxt2[iedg]];
  nump = pt->v[ip];
  numq = pt->v[MMG5_idir[iface][MMG5_iprv2[iedg]]];

#ifndef NDEBUG
  p0   = &mesh->point[nump];
  assert(p0->tag & MG_BDY);
  assert(p0->xp);

  /* Remove maybe-uninitialized value warning */
  nprvold[0] = nprvold[1] = nprvold[2] = 0.;
  ncurold[0] = ncurold[1] = ncurold[2]  = 0.;
  nprvnew[0] = nprvnew[1] = nprvnew[2] = 0.;
  ncurnew[0] = ncurnew[1] = ncurnew[2]  = 0.;
#endif

  ndepmin = ndepplus = 0;
  isminp  = isplp = 0;

  if ( !isnm ) {
    refmin  = MG_MINUS;
    refplus = MG_PLUS;
  }

  if ( !isnmint ) {
    /* prevent collapse in case surface ball has 3 triangles */
    if ( ilists <= 2 )  return 0;  // ATTENTION, Normalement, avec 2 c est bon !

    /* Surfacic ball is enumerated with first tet having (pq) as edge n° MMG5_iprv2[ip] on face iopp */
    MMG5_startedgsurfball(mesh,nump,numq,lists,ilists);
  }

  /* check tetra quality */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    iel = listv[l] / 4;
    ipp = listv[l] % 4;
    assert(iel && ipp >=0 && ipp < 4 && "unexpected tetra or local vertex idx");

    pt  = &mesh->tetra[iel];

    if ( !isnmint ) {
      if ( pt->ref == refmin ) isminp = 1;
      else if ( pt->ref == refplus ) isplp = 1;
    }

    /* Topological test for tetras of the shell */
    for (iq=0; iq<4; iq++)
      if ( pt->v[iq] == numq )  break;

    if ( iq < 4 ) {
      nbbdy = 0;
      if ( pt->xt )  {
        pxt = &mesh->xtetra[pt->xt];
      }

      for (i=0; i<4; i++) {
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) )  nbbdy++;
      }

      if ( nbbdy == 4 ) {
        /* Only one element in the domain: we don't want to delete it */
        return 0;
      }
      else if ( nbbdy == 3 ) {
        /* Topological problem triggered when one of the two faces of collapsed
         edge is the only internal one : closing a part of the domain */

        /* Identification of edge number in tetra iel */
        if ( !MMG3D_findEdge(mesh,pt,iel,numq,nump,1,NULL,&ia) ) return -1;

        i0 = MMG5_ifar[ia][0];
        i1 = MMG5_ifar[ia][1];
        if ( pt->xt && (!(pxt->ftag[i0] & MG_BDY) || !(pxt->ftag[i1] & MG_BDY)) )
          return 0;
      }

      /* Now check that the 2 faces identified by collapse are not boundary */
      if ( pt->xt && (pxt->ftag[ipp] & MG_BDY) && (pxt->ftag[iq] & MG_BDY) )
        return 0;

      if ( pt->xt )  {
        for (i=0; i<4; i++) {
          if ( i==ipp || i==iq ) {
            continue;
          }

          /*  Avoid surface crimping: check that the collapse doesn't merge 3
           *  bdy edge along a non bdy face: we have to check the edge of each
           *  shell because some MG_BDY tags may be missings due to the creation
           *  of an xtetra during a previous collapse */
          if ( !(pxt->ftag[i] & MG_BDY) ) {
            int16_t  tag0,tag1,tag2;
            MMG5_int ref0,ref1,ref2;

            tag0 = tag1 = tag2 = 0;
            ref0 = ref1 = ref2 = 0;

            if ( !MMG3D_get_shellEdgeTag(mesh,iel,MMG5_iarf[i][0],&tag0,&ref0) ) {
              fprintf(stderr,"\n  ## Error: %s: 0. unable to get edge info.\n",__func__);
              return 0;
            }
            if ( !MMG3D_get_shellEdgeTag(mesh,iel,MMG5_iarf[i][1],&tag1,&ref1) ) {
              fprintf(stderr,"\n  ## Error: %s: 1. unable to get edge info.\n",__func__);
              return 0;
            }
            if ( !MMG3D_get_shellEdgeTag(mesh,iel,MMG5_iarf[i][2],&tag2,&ref2) ) {
              fprintf(stderr,"\n  ## Error: %s: 2. unable to get edge info.\n",__func__);
              return 0;
            }
            if ( (tag0 & MG_BDY) && (tag1 & MG_BDY) && (tag2 & MG_BDY ) ) {
              return 0;
            }
          }
        }
      }

      continue;
    }

    if ( !isnmint ) {
      if ( isnm || mesh->info.iso ) {
        /* Volume test for tetras outside the shell */
        if ( (!ndepmin) && (pt->ref == refmin) ) {
          ndepmin = iel;
        }
        else if ( (!ndepplus) && (pt->ref == refplus) ) {
          ndepplus = iel;
        }
      }
    }

    /* Prevent from creating a tetra with 4 ridges metrics in aniso mode */
    if ( met && met->m && met->size == 6 ) {
      if ( (mesh->point[numq].tag & MG_GEO) && !MG_SIN_OR_NOM(mesh->point[numq].tag) ) {
        i  = ipp;
        nr = 0;
        for (iq=0; iq<3; iq++) {
          i = MMG5_inxt3[i];
          if ( (mesh->point[pt->v[i]].tag & MG_GEO) &&
               !MG_SIN_OR_NOM(mesh->point[pt->v[i]].tag) ) {
            ++nr;
          }
        }
        if ( nr==3 ) return 0;
      }
    }

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ipp] = numq;

    calold = MG_MIN(calold, pt->qual);
    if ( typchk==1 && met->m && met->size > 1 )
      caltmp = MMG5_caltet33_ani(mesh,met,pt0);
    else
      caltmp = MMG5_orcal(mesh,met,0);

    if ( caltmp < MMG5_NULKAL )  return 0;
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < MMG5_EPSOK && calnew <= calold )  return 0;
  else if ( calnew < MMG5_EPSOK || calnew < 0.3*calold )  return 0;

  if ( isnmint ) {
    return ilistv;
  }

  /* analyze surfacic ball of p */
  for (l=1; l<ilists-1; l++) {
    iel  = lists[l] / 4;
    iopp = lists[l] % 4;
    assert(iel && iopp >=0 && iopp < 4 && "unexpected tetra or local vertex idx");

    pt   = &mesh->tetra[iel];
    pxt  = &mesh->xtetra[pt->xt];
    assert(pt->xt);

    /* retrieve vertex in tetra */
    for (ip=0; ip<4; ip++)
      if ( pt->v[ip] == nump )  break;

    assert(ip<4);
    if ( ip==4 ) return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip] = numq;

    if ( !MMG5_norface(mesh,iel,iopp,ncurold) )  return 0;
    if ( !MMG5_norface(mesh,0,iopp,ncurnew) )    return 0;

    /* check normal flipping */
    ps = ncurold[0]*ncurnew[0] + ncurold[1]*ncurnew[1] + ncurold[2]*ncurnew[2];
    if ( ps < 0.0 )  return 0;

    /* check normal deviation between l and its neighbour through the edge ia */
    ia       = MMG5_idirinv[iopp][ip]; /* index of p in tria iopp */

    iedgeOpp =  MMG5_iarf[iopp][ia];

    if ( ! ( MG_GEO_OR_NOM(mesh->xtetra[pt->xt].tag[iedgeOpp])) ) {

      ier = MMG3D_normalAdjaTri(mesh,iel,iopp,ia,nadja);
      if ( ier < 0 )  return -1;
      else if (!ier ) return 0;

      devnew = nadja[0]*ncurnew[0] + nadja[1]*ncurnew[1] + nadja[2]*ncurnew[2];

      if ( devnew < mesh->info.dhd ) return 0;
    }

    if ( l == 1 ) {
      /* check normal deviation between the new first tri of the surfacic ball
       * and the surfacic triangle facing ip in the old first tri of the
       * surfacic ball (that vanishes due to the collapse)
       */
      kk    = lists[0] / 4;
      iopp2 = lists[0] % 4;
      pt1   = &mesh->tetra[kk];
      assert(pt1->xt);

      /* retrieve vertex ip in tria iopp2 */
      for (ipp=0; ipp<3; ipp++)
        if ( pt1->v[MMG5_idir[iopp2][ipp]] == nump )  break;
      assert(ipp<3);

      iedgeOpp =  MMG5_iarf[iopp2][ipp];

      if ( ! ( MG_GEO_OR_NOM(mesh->xtetra[pt1->xt].tag[iedgeOpp] ) ) ) {
        ier = MMG3D_normalAdjaTri(mesh,kk,iopp2,ipp,nadja);

        if ( ier < 0 )  return -1;
        else if ( !ier )  return 0;

        devnew = nadja[0]*ncurnew[0] + nadja[1]*ncurnew[1] + nadja[2]*ncurnew[2];
        if ( devnew < mesh->info.dhd ) {
          return 0;
        }
      }
    }
    else {
      /* check normal deviation between l and l-1 */
      ia = MMG5_iprv2[ia];         /* edge between l-1 and l, in local num of tria */
      ia = MMG5_iarf[iopp][ia];    /* edge between l-1 and l in local num of tetra */

      if ( !MG_GEO_OR_NOM(mesh->xtetra[pt->xt].tag[ia]) ) {
        devold = nprvold[0]*ncurold[0] + nprvold[1]*ncurold[1] + nprvold[2]*ncurold[2];
        devnew = nprvnew[0]*ncurnew[0] + nprvnew[1]*ncurnew[1] + nprvnew[2]*ncurnew[2];

        if ( devold < MMG5_ANGEDG ) {
          if ( devnew < devold )  {
            return 0;
          }
        }
        else if ( devnew < MMG5_ANGEDG )  {
          return 0;
        }
      }
    }

    /* check Hausdorff distance to geometric support */
    MMG5_tet2tri(mesh,iel,iopp,&tt);
    if ( l == 1 ) {
      for (i=0; i<3; i++) {
        if ( tt.v[i] == nump )  break;
      }

      assert(i<3);
      if ( i==3 ) return 0;

      /* Index of the third point of the first collapsed triangle */
      i  = MMG5_inxt2[i];
      ia = MMG5_inxt2[i];
      tag = pxt->tag[MMG5_iarf[iopp][ia]];
      tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
    }
    else if ( l == ilists-2 ) {
      for (i=0; i<3; i++) {
        if ( tt.v[i] == nump )  break;
      }
      assert(i<3);
      /* Index of the third point of the second collapsed triangle */
      i  = MMG5_iprv2[i];
      ia = MMG5_iprv2[i];
      tag = pxt->tag[MMG5_iarf[iopp][ia]];
      tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
    }

    for (i=0; i<3; i++) {
      if ( tt.v[i] == nump )  break;
    }
    assert(i<3);
    if ( i==3 ) return 0;

    tt.v[i] = numq;

    /* Local parameters for tt and iel */
    hmax  = mesh->info.hmax;
    hausd = mesh->info.hausd;
    isloc = 0;
    if ( mesh->info.parTyp & MG_Tetra ) {
      for ( kk=0; kk<mesh->info.npar; ++kk ) {
        par = &mesh->info.par[kk];

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
        for ( kk=0; kk<mesh->info.npar; ++kk ) {
          par = &mesh->info.par[kk];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax = MG_MIN(hmax,par->hmax);
          hausd = MG_MIN(hausd,par->hausd);
          break;
        }
      }
      else {
        for ( kk=0; kk<mesh->info.npar; ++kk ) {
          par = &mesh->info.par[kk];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != tt.ref ) continue;

          hmax  = par->hmax;
          hausd = par->hausd;
          isloc = 1;
          break;
        }
      }
    }

    /* If some edges must be splitted, refuse the collapse */
    if ( MMG5_chkedg(mesh,&tt,MG_GET(pxt->ori,iopp),hmax,hausd,isloc) > 0 )  return 0;

    memcpy(nprvold,ncurold,3*sizeof(double));
    memcpy(nprvnew,ncurnew,3*sizeof(double));
  }

  /* check normal deviation between the new last tri of the surfacic ball
   * and the surfacic triangle facing ip in the old last tri of the
   * surfacic ball (that vanishes due to the collapse)
   */
  kk    = lists[ilists-1] / 4;
  iopp2 = lists[ilists-1] % 4;
  pt1   = &mesh->tetra[kk];
  assert(pt1->xt);

  /* retrieve vertex ip in tria iopp2 */
  for (ipp=0; ipp<3; ipp++)
    if ( pt1->v[MMG5_idir[iopp2][ipp]] == nump )  break;
  assert(ipp<3);

  iedgeOpp =  MMG5_iarf[iopp2][ipp];

  if ( ! ( MG_GEO_OR_NOM(mesh->xtetra[pt1->xt].tag[iedgeOpp]) ) ) {

    ier = MMG3D_normalAdjaTri(mesh,kk,iopp2,ipp,nadja);
    if ( ier < 0 )  return -1;
    else if ( !ier ) return 0;

    devnew = nadja[0]*ncurnew[0] + nadja[1]*ncurnew[1] + nadja[2]*ncurnew[2];

    if ( devnew < mesh->info.dhd ) {
      return 0;
    }
  }

  if ( isnm || mesh->info.iso ) {
    ier = MMG3D_chkmanicoll(mesh,k,iface,iedg,ndepmin,ndepplus,refmin,refplus,isminp,isplp);
    if ( !ier )  return 0;
  }
  else {
    /* Topological check for surface ball */
    ier = MMG5_topchkcol_bdy(mesh,k,iface,iedg,lists,ilists);
    if ( ier<0 ) return -1;
    else if ( !ier )  return 0;
  }

  return ilistv;
}

/**
 * \param pt tetra of the shell of the edge to collapse
 * \param pxt xtetra associated to \a pt
 * \param np global index of point to collapse
 * \param nq global index of point on which we collapse
 * \param ip local index of \a np in \a pt
 * \param pt1 tetra neighbour to \a pt through \a nq
 * \param pxt1 xtetra associated to \a pt1
 * \param voyp point facing tetra \a pt in \a pt1
 *
 * Update tag and ref of the edges of \a pxt1 that belongs to the face sharing
 * \a pt1 and \a pt (face \a ip in \a pt).
 *
 */
static inline
void MMG3D_update_edgeTag(MMG5_pTetra pt,MMG5_pxTetra pxt,MMG5_int np, MMG5_int nq,
                          uint8_t ip, MMG5_pTetra pt1,MMG5_pxTetra pxt1,
                          uint8_t voyp) {

  int      i,j;
  MMG5_int p0,p1;
  uint8_t  ia,iav;
  int16_t  tag,tag1;

  /* update tags for edges */
  for ( j=0; j<3; j++ ) {
    ia = MMG5_iarf[ip][j];
    p0 = pt->v[MMG5_iare[ia][0]];
    p1 = pt->v[MMG5_iare[ia][1]];
    if ( pxt->tag[ia] ) {
      for ( i=0; i<3; i++ ) {
        iav=MMG5_iarf[voyp][i];
        if ( p0==nq ) {
          if ( ((pt1->v[MMG5_iare[iav][0]]==np) && (pt1->v[MMG5_iare[iav][1]]==p1)) ||
               ((pt1->v[MMG5_iare[iav][0]]==p1) && (pt1->v[MMG5_iare[iav][1]]==np)) )
            break;
        }
        else if ( p1==nq ) {
          if ( ((pt1->v[MMG5_iare[iav][0]]==np ) && (pt1->v[MMG5_iare[iav][1]]==p0)) ||
               ((pt1->v[MMG5_iare[iav][0]]==p0) && (pt1->v[MMG5_iare[iav][1]]==np )) )
            break;
        }
        else {
          if ( ((pt1->v[MMG5_iare[iav][0]]==p0) && (pt1->v[MMG5_iare[iav][1]]==p1)) ||
               ((pt1->v[MMG5_iare[iav][0]]==p1) && (pt1->v[MMG5_iare[iav][1]]==p0)) )
            break;
        }
      }
      assert(i!=3);
      tag1 = pxt1->tag[iav];
      tag  =  pxt->tag[ia];
      pxt1->tag[iav] |= pxt->tag[ia];
      if( ((tag1 & MG_REQ) && !(tag1 & MG_NOSURF)) ||
          (( tag & MG_REQ) && !( tag & MG_NOSURF)) )
        pxt1->tag[iav] &= ~MG_NOSURF;
      pxt1->edg[iav] = MG_MAX ( pxt1->edg[iav],pxt->edg[ia] );
    }
  }
}

/**
 * \param mesh pointer toward the mesh
 * \param start tetra from which we start to travel
 * \param na edge vertex
 * \param nb edge vertex
 * \param tag new edge tag
 * \param ref new edge ref
 * \param piv global index of the pivot to set the sense of travel
 * \param adj index of adjacent tetra for the travel
 *
 * \return -1 if fail, \a start if shell has been completely travelled, 0 otherwise
 *
 * Update tag and ref of the edge \a na \a nb from tetra \a start by traveling
 * its shell in one direction (given by the pivot \a piv).
 *
 */
static inline
MMG5_int MMG3D_update_shellEdgeTag_oneDir(MMG5_pMesh  mesh,MMG5_int start, MMG5_int na, MMG5_int nb,
                                     int16_t tag,MMG5_int ref, MMG5_int piv,MMG5_int adj) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     *adja;
  int16_t      xtag;
  int8_t       i;

  assert ( tag & MG_BDY && "Unexpected non boundary tag");

  while ( adj && (adj != start) ) {
    pt     = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,1,NULL,&i) ) return -1;

    /* update edge ref and tag */
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];

      /* if tag and edge are already consistent, no need to update the shell (do
       * not consider edges with tag 0 as they com from the creation of a xtetra
       * during a previous collapse and are not updated) */
      xtag = pxt->tag[i] | MG_BDY;

      if ( pxt->tag[i] & MG_BDY ) {
        if ( xtag == tag &&  pxt->edg[i] == ref ) {
          return start;
        }
      }

      pxt->edg[i] = ref;
      pxt->tag[i] |= tag;
      if( ((xtag & MG_REQ) && !(xtag & MG_NOSURF)) ||
          (( tag & MG_REQ) && !( tag & MG_NOSURF)))
        pxt->tag[i] &= ~MG_NOSURF;
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ MMG5_ifar[i][1] ] /4;
      piv = pt->v[ MMG5_ifar[i][0] ];
    }
  }

  return adj;
}

/**
 * \param mesh pointer toward the mesh
 * \param start tetra from which we start to travel
 * \param ia local index of edge that must be updated
 * \param tag new edge tag
 * \param ref new edge ref
 * \return 1 if success, 0 if fail.
 *
 * Update tag and ref of the edge \ia of tetra \a start by traveling its shell.
 *
 */
static inline
int MMG3D_update_shellEdgeTag(MMG5_pMesh  mesh,MMG5_int start, int8_t ia,int16_t tag,MMG5_int ref) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     piv,na,nb,adj,*adja;
  int16_t      xtag;

  pt   = &mesh->tetra[start];

  assert( start >= 1 &&  MG_EOK(pt) && "invalid tetra" );
  assert ( tag & MG_BDY && "Unexpected non boundary tag");

  pxt  = NULL;
  na   = pt->v[MMG5_iare[ia][0]];
  nb   = pt->v[MMG5_iare[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];

    /* if tag and edge are already consistent, no need to update the shell (do
     * not consider edges with tag 0 as they com from the creation of a xtetra
     * during a previous collapse and are not updated) */
    xtag = pxt->tag[ia] | MG_BDY;

    if ( pxt->tag[ia] & MG_BDY ) {
      if ( xtag == tag &&  pxt->edg[ia] == ref ) {
        return 1;
      }
    }

    pxt->tag[ia] |= tag;
    if( ((xtag & MG_REQ) && !(xtag & MG_NOSURF)) ||
        (( tag & MG_REQ) && !( tag & MG_NOSURF)))
      pxt->tag[ia] &= ~MG_NOSURF;
    pxt->edg[ia]  = ref;
  }

  /* Travel in one direction */
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][0]] / 4;
  piv = pt->v[MMG5_ifar[ia][1]];

  adj =   MMG3D_update_shellEdgeTag_oneDir(mesh,start,na,nb,tag,ref,piv,adj);

  /* If all shell has been travelled, stop, else, travel it the other sense */
  if ( adj == start )  return 1;
  else if ( adj < 0 ) return 0;

  assert(!adj);

  pt = &mesh->tetra[start];
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][1]] / 4;
  piv = pt->v[MMG5_ifar[ia][0]];

  adj = MMG3D_update_shellEdgeTag_oneDir(mesh,start,na,nb,tag,ref,piv,adj);

  if ( adj < 0 ) return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param list pointer toward the ball of the point
 * \param ilist number of elements in the ball of the point
 * \param indq local index of the point on which we collapse
 * \param typchk type of check performed depending on the remeshing step
 *
 * \return np the index of the collpased point if success, 0 if we cannot
 * collapse, -1 if we fail.
 *
 * Collapse vertex p = list[0]%4 of tetra list[0]/4 over vertex indq of tetra
 * list[0]/4.  Only physical tests (positive jacobian) are done
 * (i.e. approximation of the surface, etc... must be performed outside).
 *
 */
MMG5_int MMG5_colver(MMG5_pMesh mesh,MMG5_pSol met,int64_t *list,int ilist,int8_t indq,int8_t typchk) {
  MMG5_pTetra          pt,pt1;
  MMG5_pxTetra         pxt,pxt1;
  MMG5_xTetra          xt,xts;
  MMG5_int             iel,jel,pel,qel,k,np,nq,*adja;
  uint8_t              ip,iq,j,voy,voyp,voyq;

  /* coledge[i] contains the local indices of edges that will be merged by the
   * collapse corresponding with the configuration i. The edge coledge[i][0] is
   * merged with the edge coledge[i][1] a,d the edge coledge[i][2] is merged
   * with edge coledge[i][3].
   * Config 0: merge of vertices 0 and 1
   * config 1: merge of vertices 0 and 2
   * config 2: merge of vertices 0 and 3
   * config 3: merge of vertices 1 and 2
   * config 4: merge of vertices 1 and 3
   * config 5: merge of vertices 2 and 3
   */
  const int8_t MMG5_coledge[6][4] = {
  {1,3,2,4}, {0,3,2,5}, {0,4,1,5},{0,1,4,5}, {0,2,3,5}, {4,3,2,1} };

  iel = list[0] / 4;
  ip  = list[0] % 4;
  pt  = &mesh->tetra[iel];
  np  = pt->v[ip];
  nq  = pt->v[indq];

  /* Mark elements of the shell of edge (pq) */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];

    for (j=0; j<3; j++) {
      iq = MMG5_idir[ip][j];
      if ( pt->v[iq] == nq ) {

        list[k] *= -1;
        break;
      }
    }
  }

  /* Avoid recreating existing elt */
  for (k=0; k<ilist; k++) {

    if ( list[k] < 0 )  continue;

    iel = list[k] / 4;
    ip  = list[k] % 4;

    adja = &mesh->adja[4*(iel-1)+1];
    jel  = adja[ip];
    if ( !jel )  continue;

    jel /= 4;
    voy  = adja[ip] % 4;
    pt = &mesh->tetra[jel];
    if (pt->v[voy] == nq) {
      return 0;
    }
  }

  /** Merge tags and references of edges that will merge due to the collapse
   * (the shell of each edge is travelled so each xtetra of the shell is
   * updated). Note that it can't be done in the previous loop because the mesh
   * would be corrupted if we stop the collapse. It can't neither be done in the
   * next loop because we start to delete the elements of the shell. */
  for (k=0; k<ilist; k++) {

    if ( list[k] > 0 )  continue;

    iel = -list[k] / 4;
    ip  = -list[k] % 4;
    pt  = &mesh->tetra[iel];

    for (j=0; j<3; j++) {
      iq = MMG5_idir[ip][j];
      if ( pt->v[iq] == nq ) {
        break;
      }
    }

    /* The configuration is computed by setting the ip and iq bits to 1 */
    int8_t cfg = 0;
    MG_SET(cfg,ip);
    MG_SET(cfg,iq);

    const int8_t *coled;
    switch(cfg) {
    case 3:
      /* collapse of vertices 0 and 1 */
      coled = MMG5_coledge[0];
      break;
    case 5:
      /* collapse of vertices 0 and 2 */
      coled = MMG5_coledge[1];
      break;
    case 9:
      /* collapse of vertices 0 and 3 */
      coled = MMG5_coledge[2];
      break;
    case 6:
      /* collapse of vertices 1 and 2 */
      coled = MMG5_coledge[3];
      break;
    case 10:
      /* collapse of vertices 1 and 3 */
      coled = MMG5_coledge[4];
      break;
    case 12:
      /* collapse of vertices 2 and 3 */
      coled = MMG5_coledge[5];
      break;
    default:
      assert ( 0 && "Unexpected collapse configuration");
    }

    int l = 0;
    for ( l=0; l<3; l+=2 ) {
      /* when l=0 we update 2 edges that will be merged together, when l=2 we
       * update the two others. In practice we don't care which edge comes
       * from ip and which one from iq, it only matters that iped and iqed
       * will be merged at the end of the collapse. */
      int iped = coled[l+0];
      int iqed = coled[l+1];

      int16_t  tagip = 0;
      MMG5_int refip = 0;
      int16_t  tagiq = 0;
      MMG5_int refiq = 0;

      if ( !MMG3D_get_shellEdgeTag(mesh,iel,iped,&tagip,&refip) ) {
        fprintf(stderr,"\n  ## Error: %s: 1. unable to get edge info.\n",__func__);
        return 0;
      }
      if ( !MMG3D_get_shellEdgeTag(mesh,iel,iqed,&tagiq,&refiq) ) {
        fprintf(stderr,"\n  ## Error: %s: 2. unable to get edge info.\n",__func__);
        return 0;
      }

      if ( (tagip != tagiq) || (refip != refiq) ) {
        /* Update edge infos */
        tagip |= tagiq;
        refip = (refip > refiq )? refip : refiq;

        /* If the xtetra info are consistent when entering in colver taged and
         * refed contains suitable values and we can travel the shell of
         * iped and iqed to update all the needed info on edges */
        if ( !MMG3D_update_shellEdgeTag(mesh,iel,iped,tagip,refip) ) {
          fprintf(stderr,"\n  ## Error: %s: 1. unable to update edge info.\n",__func__);
          return 0;
        }
        if ( !MMG3D_update_shellEdgeTag(mesh,iel,iqed,tagip,refip) ) {
          fprintf(stderr,"\n  ## Error: %s: 1. unable to update edge info.\n",__func__);
          return 0;
        }
      }
#ifndef NDEBUG
      else {
        assert ( tagip == tagiq && refip == refiq );
        if ( mesh->info.ddebug && (tagip & MG_BDY) ) {
          MMG3D_chk_shellEdgeTag(mesh,iel,iped,tagip,refip);
          MMG3D_chk_shellEdgeTag(mesh,iel,iqed,tagip,refip);
        }
      }
#endif
    }
  }

  /* deal with the shell of edge (pq) and the implied updates */
  for (k=0; k<ilist; k++) {
    if ( list[k] > 0 )  continue;
    iel = (-list[k]) / 4;
    ip  = (-list[k]) % 4;
    pt  = &mesh->tetra[iel];

    iq  = MMG5_inxt3[ip];
    for (j=0; j<3; j++) {
      if ( pt->v[iq] == nq )  break;
      iq = MMG5_inxt3[iq];
    }
    assert(j<3);


    adja = &mesh->adja[4*(iel-1)+1];

    /* pel = neighbour of iel that belongs to ball of p \setminus shell, same for qel */
    pel  = adja[iq] / 4;
    voyp = adja[iq] % 4;
    qel  = adja[ip] / 4;
    voyq = adja[ip] % 4;

    /* Update adjacency relations */
    if ( pel ) {
      adja = &mesh->adja[4*(pel-1)+1];
      adja[voyp] = 4*qel+voyq;
    }
    if ( qel ) {
      adja = &mesh->adja[4*(qel-1)+1];
      adja[voyq] = 4*pel+voyp;
    }

    /* Update references for faces (one in pel) ;
       possibly, creation of a new field pxt for pel must be carried out */
    if ( pel ) {
      pt1 = &mesh->tetra[pel];
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        memcpy(&xts,pxt,sizeof(MMG5_xTetra));
        if ( pt1->xt > 0 ) {
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyp] = MG_MAX(pxt1->ref[voyp],pxt->ref[ip]);
          pxt1->ftag[voyp] = pxt1->ftag[voyp] | pxt->ftag[ip];

          /* the face voyp of pxt1 is projected into the face ip of pxt so it
           * takes its orientation */
          if ( pxt->ftag[ip] & MG_BDY ) {
            if ( MG_GET(pxt->ori,ip) )
              MG_SET( pxt1->ori,voyp );
            else
              MG_CLR( pxt1->ori,voyp );
          }
          /* correct the orientation if the new adjacent has a different
           * reference */
          if( qel ) {
            if( pt1->ref > mesh->tetra[qel].ref )
              MG_SET( pxt1->ori,voyp );
            else if( pt1->ref < mesh->tetra[qel].ref )
              MG_CLR( pxt1->ori,voyp );
          }

#ifndef NDEBUG
          else {
            /* Check that a non parallel external boundary face has always a good
             * orientation */
            if ( (pxt1->ftag[voyp] & MG_BDY) && !(pxt1->ftag[voyp] & MG_PARBDY) ) {
              assert ( MG_GET( pxt1->ori,voyp ) );
            }
          }
#endif

          /* update tags for edges */
          MMG3D_update_edgeTag(pt,pxt,np,nq,ip,pt1,pxt1,voyp);
        }
        else {
          /* shell tet has a xtetra and pel exists but pel don't have a xtetra */
          pxt1 = &xt;
          memset(pxt1,0,sizeof(MMG5_xTetra));
          pxt1->ref[voyp] = pxt->ref[ip];
          pxt1->ftag[voyp] = pxt->ftag[ip];
          pxt1->ori = 15;
          if ( !MG_GET(pxt->ori,ip) )  MG_CLR(pxt1->ori,voyp);

          /* update tags for edges */
          MMG3D_update_edgeTag(pt,pxt,np,nq,ip,pt1,pxt1,voyp);

          /* Recover the already used place by pxt: now pel has a xtetra but the
           * edges coming from voyp (not directly involved in the collapse) will
           * have the tag 0 that may be non consistent with other tags in their
           * respective shells. */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(MMG5_xTetra));
        }
      }
      else {
        /* Shell tet don't have a xtetra: only the values of pel corresponding
         * to pt become 0 */
        if ( pt1->xt > 0 ) {
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyp]  = 0;
          pxt1->ftag[voyp] = 0;
          MG_SET(pxt1->ori,voyp);
        }
      }

      if ( qel ) {
        pt1 = &mesh->tetra[qel];
        if ( pt->xt ) {
          pxt = &xts;
          if ( pt1->xt > 0 ) {
            pxt1 = &mesh->xtetra[pt1->xt];
            pxt1->ref[voyq]  = MG_MAX(pxt1->ref[voyq],pxt->ref[iq]);
            pxt1->ftag[voyq] = (pxt1->ftag[voyq] | pxt->ftag[iq]);

            /* the face voyq of pxt1 is projected into the face iq of pxt so it
             * takes its orientation */
            if ( pxt->ftag[iq] & MG_BDY ) {
              if ( MG_GET(pxt->ori,iq) )
                MG_SET( pxt1->ori,voyq );
              else
                MG_CLR( pxt1->ori,voyq );
            }
            /* correct the orientation if the new adjacent has a different
             * reference */
            if( pel ) {
              if( pt1->ref > mesh->tetra[pel].ref )
                MG_SET( pxt1->ori,voyq );
              else if( pt1->ref < mesh->tetra[pel].ref )
                MG_CLR( pxt1->ori,voyq );
            }

#ifndef NDEBUG
            else {
              /* Check that a non parallel external boundary face has always a good
               * orientation */
              if ( (pxt1->ftag[voyq] & MG_BDY) && !(pxt1->ftag[voyq] & MG_PARBDY) ) {
                assert ( MG_GET( pxt1->ori,voyq ) );
              }
            }
#endif

            /* update tags for edges */
            MMG3D_update_edgeTag(pt,pxt,nq,np,iq,pt1,pxt1,voyq);
          }
          else {
            /* pel exists, shell tet has a xtetra but qel doesn't have boundary
             * tetra: create it */
            pxt1 = &xt;
            memset(pxt1,0,sizeof(MMG5_xTetra));
            pxt1->ref[voyq] = pxt->ref[iq];
            pxt1->ftag[voyq] = pxt->ftag[iq];
            pxt1->ori = 15;
            if ( !MG_GET(pxt->ori,iq) )  MG_CLR(pxt1->ori,voyq);

            /* update tags for edges */
            MMG3D_update_edgeTag(pt,pxt,nq,np,iq,pt1,pxt1,voyq);

            /* Create new field xt */
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;return -1;);
            }
            /* Now qel has a xtetra but the edges coming from voyq (not directly
             * involved in the collapse) will have the tag 0 that may be non
             * consistent with other tags in their respective shells. */
            pt1->xt = mesh->xt;
            pxt = &mesh->xtetra[pt1->xt];
            memcpy(pxt,pxt1,sizeof(MMG5_xTetra));
          }
        }
        else {
          /* pel exist but shell tet doesn't have a boundary tetra: only the
           * values of qel corresponding to pt become 0 */
          if ( pt1->xt > 0 ) {
            pxt1 = &mesh->xtetra[pt1->xt];
            pxt1->ref[voyq]  = 0;
            pxt1->ftag[voyq] = 0;
            MG_SET(pxt1->ori,voyq);
          }
        }
      }
    }
    else {
      /* pel==0: No adjacent through face iq */
      assert(pt->xt);
      pxt = &mesh->xtetra[pt->xt];
      if ( qel ) {
        pt1 = &mesh->tetra[qel];
        if ( pt1->xt > 0 ) {
          pxt1 = &mesh->xtetra[pt1->xt];
          pxt1->ref[voyq]  = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq];

          MG_SET(pxt1->ori,voyq);

          /* update tags for edges */
          MMG3D_update_edgeTag(pt,pxt,nq,np,iq,pt1,pxt1,voyq);
        }
        else {
          pxt1 = &xt;
          memset(pxt1,0,sizeof(MMG5_xTetra));
          pxt1->ref[voyq]  = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq];
          pxt1->ori = 15;

          /* update tags for edges */
          MMG3D_update_edgeTag(pt,pxt,nq,np,iq,pt1,pxt1,voyq);

          /* Recover the already used place by pxt */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(MMG5_xTetra));
        }
      }
    }
    if ( !MMG3D_delElt(mesh,iel) ) {
      return -1;
    }
  }

  /* Update vertices coordinates for elements that do not belong to the shell of (pq) */
  for (k=0; k<ilist;  k++) {
    if ( list[k] < 0 )  continue;
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];
    pt->v[ip] = nq;
    if ( typchk==1 && met->m && met->size > 1 )
      pt->qual=MMG5_caltet33_ani(mesh,met,pt);
    else
      pt->qual=MMG5_orcal(mesh,met,iel);
    pt->mark=mesh->mark;
  }

  return np;
}
