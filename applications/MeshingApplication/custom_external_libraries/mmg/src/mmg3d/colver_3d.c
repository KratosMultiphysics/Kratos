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

#include "inlined_functions_3d.h"

extern char  ddb;

/** Check whether collapse ip -> iq could be performed, ip internal ;
 *  'mechanical' tests (positive jacobian) are not performed here */
int _MMG5_chkcol_int(MMG5_pMesh mesh,MMG5_pSol met,int k,char iface,
                     char iedg,int *list,int ilist,char typchk) {
  MMG5_pTetra   pt,pt0;
  MMG5_pPoint   p0;
  double   calold,calnew,caltmp,lon,ll;
  int      j,iel,nq,nr;
  char     i,jj,ip,iq;

  ip  = _MMG5_idir[iface][_MMG5_inxt2[iedg]];
  iq  = _MMG5_idir[iface][_MMG5_iprv2[iedg]];
  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  nq  = pt->v[iq];

  lon = 1.e20;
  if ( typchk == 2 && met->m ) {
    lon = _MMG5_lenedg(mesh,met,_MMG5_iarf[iface][iedg],pt);

    if ( !lon ) return(0);

    lon = MG_MIN(lon,_MMG5_LSHRT);
    lon = MG_MAX(1.0/lon,_MMG5_LLONG);
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
    if ( mesh->info.fem ) {
      p0 = &mesh->point[nq];
      if ( p0->tag & MG_BDY ) {
        i = ip;
        for (jj=0; jj<3; jj++) {
          i = _MMG5_inxt3[i];
          p0 = &mesh->point[pt->v[i]];
          if ( p0->tag & MG_BDY )  return(0);
        }
      }
    }

    /* Prevent from creating a tetra with 4 ridges vertices */
    p0 = &mesh->point[nq];
    if ( p0->tag & MG_GEO ) {
      i  = ip;
      nr = 0;
      for (jj=0; jj<3; jj++) {
        i = _MMG5_inxt3[i];
        p0 = &mesh->point[pt->v[i]];
        if ( p0->tag & MG_GEO ) ++nr;
      }
      if ( nr==3 ) return(0);
    }

    pt0->v[ip] = nq;

    calold = MG_MIN(calold,pt->qual);
    if ( typchk==1 && met->m && met->size > 1 )
      caltmp = _MMG5_caltet33_ani(mesh,met,pt0);
    else
      caltmp = _MMG5_orcal(mesh,met,0);

    if ( caltmp < _MMG5_EPSD )  return(0);
    calnew = MG_MIN(calnew,caltmp);
    /* check length */
    if ( typchk == 2 && met->m ) {
      for (jj=0; jj<6; jj++) {
        /* Rough evaluation of edge length (doesn't take into account if some of
         * the modified edges of pt0 are boundaries): for a more precise
         * computation, we need to update the edge tags of pt0.  */
        ll = _MMG5_lenedgspl(mesh,met,jj,pt0);
        if ( (!ll) || (ll > lon) )
          return(0);
      }
    }
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )  return(0);
  else if ( calnew < _MMG5_NULKAL || calnew < 0.3*calold )  return(0);

  return(ilist);
}

/** Topological check on the surface ball of np and nq in collapsing np->nq ;
 *  iface = boundary face on which lie edge iedg - in local face num.
 *  (pq, or ia in local tet notation) */
static int
_MMG5_topchkcol_bdy(MMG5_pMesh mesh,int k,int iface,char iedg,int *lists,int ilists) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  int      nump,numq,piv0,piv,iel,jel,nap,nbp,naq,nbq,nro,adj,*adja;
  char     ip,iq,ipiv,iopp,i,j,jface,ipa,ipb,isface;

  pt = &mesh->tetra[k];
  ip = _MMG5_idir[iface][_MMG5_inxt2[iedg]];
  iq = _MMG5_idir[iface][_MMG5_iprv2[iedg]];
  nump = pt->v[ip];
  numq = pt->v[iq];

  /* Pivot in enumeration of the surface ball of np */
  ipiv = _MMG5_idir[iface][_MMG5_inxt2[_MMG5_idirinv[iface][ip]]];
  piv0 = pt->v[ipiv];

  /* Surface ball has been enumerated as f1,...,f2 - f1,f2 = both triangles of surface shell */
  if ( piv0 == numq ) {
    /*  Point nap, facing the first vanishing face in surface ball of p */
    nro = pt->v[_MMG5_idir[iface][_MMG5_iprv2[_MMG5_idirinv[iface][ip]]]];

    jel = lists[1] / 4;
    jface = lists[1] % 4;

    pt = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = _MMG5_idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro ) break;
    }
    assert(j<3);

    nap = pt->v[i];

    /* Unfold shell of (nq,nro), starting from (k,iface), with pivot np */
    adj = k;
    piv = nump;
    do {
      iel = adj;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = _MMG5_iare[i][0];
        ipb = _MMG5_iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro)  && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        adj  = adja[ _MMG5_ifar[i][0] ] / 4;
        ipiv = _MMG5_ifar[i][1];
        iopp = _MMG5_ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ _MMG5_ifar[i][1] ] / 4;
        ipiv = _MMG5_ifar[i][0];
        iopp = _MMG5_ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    naq = piv;
    if ( nap == naq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nap=%d (%d in saveMesh)\n",
        k,_MMG3D_indElt(mesh,k),nap,_MMG3D_indPt(mesh,nap));*/
      return(0);
    }

    /*  Point nbp, facing the second vanishing face in surface ball of p */
    jel   = lists[ilists-1] / 4;
    jface = lists[ilists-1] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = _MMG5_idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != numq )  break;
    }
    assert(j<3);

    nro   = pt->v[i];
    jel   = lists[ilists-2] / 4;
    jface = lists[ilists-2] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = _MMG5_idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro )  break;
    }
    assert(j<3);

    nbp = pt->v[i];

    /* Unfold shell of (nq,nro), starting from (jel,jface), with pivot np */
    adj = lists[ilists-1] / 4;
    piv=  nump;
    do {
      iel  = adj;
      pt   = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = _MMG5_iare[i][0];
        ipb = _MMG5_iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro) && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        adj  = adja[ _MMG5_ifar[i][0] ] / 4;
        ipiv = _MMG5_ifar[i][1];
        iopp = _MMG5_ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ _MMG5_ifar[i][1] ] / 4;
        ipiv = _MMG5_ifar[i][0];
        iopp = _MMG5_ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt    = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    nbq = piv;
    if ( nbp == nbq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nbp=%d (%d in saveMesh)\n",
        k,_MMG3D_indElt(mesh,k),nbp,_MMG3D_indPt(mesh,nbp));*/
      return(0);
    }
  }
  /* Surface ball has been enumerated as f1,f2,... */
  else {
    /*  Point nap, facing the fist vanishing face in surface ball of p */
    nro   = piv0;
    jel   = lists[ilists-1] / 4;
    jface = lists[ilists-1] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = _MMG5_idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro )  break;
    }
    assert(j<3);

    nap = pt->v[i];

    /* Unfold shell of (nq,nro), starting from (k,iface), with pivot np */
    adj = k;
    piv = nump;
    do {
      iel = adj;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = _MMG5_iare[i][0];
        ipb = _MMG5_iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro) && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        adj  = adja[ _MMG5_ifar[i][0] ] / 4;
        ipiv = _MMG5_ifar[i][1];
        iopp = _MMG5_ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ _MMG5_ifar[i][1] ] / 4;
        ipiv = _MMG5_ifar[i][0];
        iopp = _MMG5_ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt    = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    naq = piv;
    if ( nap == naq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nap=%d (%d in saveMesh)\n",
        k,_MMG3D_indElt(mesh,k),nap,_MMG3D_indPt(mesh,nap));*/
      return(0);
    }

    /*  Point nbp, facing the second vanishing face in surface ball of p */
    jel   = lists[1] / 4;
    jface = lists[1] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = _MMG5_idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != numq )  break;
    }
    assert(j<3);

    nro   = pt->v[i];
    jel   = lists[2] / 4;
    jface = lists[2] % 4;
    pt    = &mesh->tetra[jel];
    for (j=0; j<3; j++) {
      i = _MMG5_idir[jface][j];
      if ( pt->v[i] != nump && pt->v[i] != nro )  break;
    }
    assert(j<3);

    nbp = pt->v[i];

    /* Unfold shell of (nq,nro), starting from lists[1], with pivot np */
    adj = lists[1] / 4;
    piv =  nump;
    do {
      iel  = adj;
      pt   = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      /* Identification of edge number in tetra iel */
      for (i=0; i<6; i++) {
        ipa = _MMG5_iare[i][0];
        ipb = _MMG5_iare[i][1];
        if ( ((pt->v[ipa] == numq) && (pt->v[ipb] == nro)) || ((pt->v[ipa] == nro) && (pt->v[ipb] == numq))  ) break;
      }
      assert(i<6);

      /* set sense of travel */
      if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
        adj  = adja[ _MMG5_ifar[i][0] ] / 4;
        ipiv = _MMG5_ifar[i][1];
        iopp = _MMG5_ifar[i][0];
        piv  = pt->v[ipiv];
      }
      else {
        adj  = adja[ _MMG5_ifar[i][1] ] / 4;
        ipiv = _MMG5_ifar[i][0];
        iopp = _MMG5_ifar[i][1];
        piv  = pt->v[ipiv];
      }

      isface = 0;
      if ( pt->xt ) {
        pxt    = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && ( adj != k ) && !isface );

    nbq = piv;
    if ( nbp == nbq ) {
      /*printf("%s: %d: On devrait rarement passer ici:",__FILE__,__LINE__);
        printf(" k=%d (%d in saveMesh), nap=%d (%d in saveMesh)\n",
        k,_MMG3D_indElt(mesh,k),nap,_MMG3D_indPt(mesh,nap));*/
      return(0);
    }
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element in which we collapse.
 * \param iface face through wich we perform the collapse
 * \param iedg edge to collapse
 * \param listv pointer toward the list of the tetra in the ball of \a p0.
 * \param typchk  typchk type of checking permformed for edge length
 * (hmax or _MMG5_LLONG criterion).
 *
 * Check whether collapse ip -> iq could be performed, ip boundary point ;
 *  'mechanical' tests (positive jacobian) are not performed here ;
 *  iface = boundary face on which lie edge iedg - in local face num.
 *  (pq, or ia in local tet notation).
 *
 */
int _MMG5_chkcol_bdy(MMG5_pMesh mesh,MMG5_pSol met,int k,char iface,
                     char iedg,int *listv,int ilistv,int *lists,int ilists,
                     char typchk) {
  MMG5_pTetra        pt,pt0;
  MMG5_pxTetra       pxt;
  MMG5_pPoint        p0;
  MMG5_Tria          tt;
  MMG5_pPar          par;
  double        calold,calnew,caltmp,nprvold[3],nprvnew[3],ncurold[3],ncurnew[3];
  double        ps,devold,devnew,hmax,hausd;
  int           ipp,nump,numq,l,iel,kk;
  int           nr,nbbdy,ndepmin,ndepplus,isloc;
  int16_t       tag;
  char          iopp,ia,ip,i,iq,i0,i1,ier,isminp,isplp;

  pt   = &mesh->tetra[k];
  pxt  = 0;
  pt0  = &mesh->tetra[0];
  ia   = _MMG5_iarf[iface][iedg];
  ip   = _MMG5_idir[iface][_MMG5_inxt2[iedg]];
  nump = pt->v[ip];
  numq = pt->v[_MMG5_idir[iface][_MMG5_iprv2[iedg]]];
  p0   = &mesh->point[nump];
  assert(p0->tag & MG_BDY);
  assert(p0->xp);

  ndepmin = ndepplus = 0;
  isminp  = isplp = 0;

  /* prevent collapse in case surface ball has 3 triangles */
  if ( ilists <= 2 )  return(0);  // ATTENTION, Normalement, avec 2 c est bon !

  /* Surfacic ball is enumerated with first tet having (pq) as edge n° _MMG5_iprv2[ip] on face iopp */
  _MMG5_startedgsurfball(mesh,nump,numq,lists,ilists);

  /* check tetra quality */
  calold = calnew = DBL_MAX;
  for (l=0; l<ilistv; l++) {
    iel = listv[l] / 4;
    ipp = listv[l] % 4;
    pt  = &mesh->tetra[iel];

    if ( pt->ref == MG_MINUS ) isminp = 1;
    else if ( pt->ref == MG_PLUS ) isplp = 1;

    /* Topological test for tetras of the shell */
    for (iq=0; iq<4; iq++)
      if ( pt->v[iq] == numq )  break;

    if ( iq < 4 ) {
      nbbdy = 0;
      if ( pt->xt )  pxt = &mesh->xtetra[pt->xt];
      for (i=0; i<4; i++) {
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) )  nbbdy++;
      }

      /* Topological problem triggered when one of the two faces of collapsed edge is the only
         internal one : closing a part of the domain */
      if ( nbbdy == 4 )
        return(0);
      else if ( nbbdy == 3 ) {
        for (ia=0; ia<6; ia++) {
          i0 = _MMG5_iare[ia][0];
          i1 = _MMG5_iare[ia][1];
          if ( ((pt->v[i0] == nump) && (pt->v[i1] == numq)) ||
               ((pt->v[i0] == numq) && (pt->v[i1] == nump)) )
            break;
        }
        assert(ia < 6);
        i0 = _MMG5_ifar[ia][0];
        i1 = _MMG5_ifar[ia][1];
        if ( pt->xt && (!(pxt->ftag[i0] & MG_BDY) || !(pxt->ftag[i1] & MG_BDY)) )
          return(0);
      }

      /* Now check that the 2 faces identified by collapse are not boundary */
      if ( pt->xt && (pxt->ftag[ipp] & MG_BDY) && (pxt->ftag[iq] & MG_BDY) )
        return(0);

      continue;
    }

    /* Volume test for tetras outside the shell */
    if ( mesh->info.iso ) {
      if ( !ndepmin && pt->ref == MG_MINUS )
        ndepmin = iel;
      else if ( !ndepplus && pt->ref == MG_PLUS )
        ndepplus = iel;
    }

    /* Prevent from creating a tetra with 4 ridges vertices */
    if ( mesh->point[numq].tag & MG_GEO ) {
      i  = ipp;
      nr = 0;
      for (iq=0; iq<3; iq++) {
        i = _MMG5_inxt3[i];
        if ( mesh->point[pt->v[i]].tag & MG_GEO ) ++nr;
      }
      if ( nr==3 ) return(0);
    }

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ipp] = numq;

    calold = MG_MIN(calold, pt->qual);
    if ( typchk==1 && met->m && met->size > 1 )
      caltmp = _MMG5_caltet33_ani(mesh,met,pt0);
    else
      caltmp = _MMG5_orcal(mesh,met,0);

    if ( caltmp < _MMG5_EPSD )  return(0);
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calold < _MMG5_NULKAL && calnew <= calold )  return(0);
  else if ( calnew < _MMG5_NULKAL || calnew < 0.3*calold )  return(0);

  /* analyze surfacic ball of p */
  for (l=1; l<ilists-1; l++) {
    iel  = lists[l] / 4;
    iopp = lists[l] % 4;
    pt   = &mesh->tetra[iel];
    pxt = &mesh->xtetra[pt->xt];
    assert(pt->xt);

    /* retrieve vertex in tetra */
    for (ip=0; ip<4; ip++)
      if ( pt->v[ip] == nump )  break;
    assert(ip<4);

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip] = numq;

    if ( !_MMG5_norface(mesh,iel,iopp,ncurold) )  return(0);
    if ( !_MMG5_norface(mesh,0,iopp,ncurnew) )    return(0);

    /* check normal flipping */
    ps = ncurold[0]*ncurnew[0] + ncurold[1]*ncurnew[1] + ncurold[2]*ncurnew[2];
    if ( ps < 0.0 )  return(0);

    /* check normal deviation */
    if ( l > 1 ) {
      ia = _MMG5_idirinv[iopp][ip]; /* index of p in tria iopp */
      ia = _MMG5_iprv2[ia];         /* edge between l-1 and l, in local num of tria */
      ia = _MMG5_iarf[iopp][ia];    /* edge between l-1 and l in local num of tetra */

      if ( (!pt->xt) || (!(mesh->xtetra[pt->xt].tag[ia] & MG_GEO)) ) {

        devold = nprvold[0]*ncurold[0] + nprvold[1]*ncurold[1] + nprvold[2]*ncurold[2];
        devnew = nprvnew[0]*ncurnew[0] + nprvnew[1]*ncurnew[1] + nprvnew[2]*ncurnew[2];
        if ( devold < _MMG5_ANGEDG ) {
          if ( devnew < devold )  return(0);
        }
        else if ( devnew < _MMG5_ANGEDG )  return(0);
      }
    }

    /* check Hausdorff distance to geometric support */
    _MMG5_tet2tri(mesh,iel,iopp,&tt);
    if ( l == 1 ) {
      for (i=0; i<3; i++) {
        if ( tt.v[i] == nump )  break;
      }
      assert(i<3);
      /* Index of the third point of the first collapsed triangle */
      i  = _MMG5_inxt2[i];
      ia = _MMG5_inxt2[i];
      tag = pxt->tag[_MMG5_iarf[iopp][ia]];
      tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
    }
    else if ( l == ilists-2 ) {
      for (i=0; i<3; i++) {
        if ( tt.v[i] == nump )  break;
      }
      assert(i<3);
      /* Index of the third point of the first collapsed triangle */
      i  = _MMG5_iprv2[i];
      ia = _MMG5_iprv2[i];
      tag = pxt->tag[_MMG5_iarf[iopp][ia]];
      tt.tag[ia] = MG_MAX(tt.tag[ia],tag);
    }

    for (i=0; i<3; i++) {
      if ( tt.v[i] == nump )  break;
    }
    assert(i<3);
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

    if ( _MMG5_chkedg(mesh,&tt,MG_GET(pxt->ori,iopp),hmax,hausd,isloc) )  return(0);

    memcpy(nprvold,ncurold,3*sizeof(double));
    memcpy(nprvnew,ncurnew,3*sizeof(double));
  }

  /* Ensure collapse does not lead to a non manifold configuration (case of implicit surface)*/
  if ( mesh->info.iso ) {
    ier = _MMG5_chkmanicoll(mesh,k,iface,iedg,ndepmin,ndepplus,isminp,isplp);
    if ( !ier )  return(0);
  }
  /* Topological check for surface ball */
  else {
    ier = _MMG5_topchkcol_bdy(mesh,k,iface,iedg,lists,ilists);
    if ( !ier )  return(0);
  }

  return(ilistv);
}

/** Collapse vertex p = list[0]%4 of tetra list[0]/4 over vertex indq of tetra list[0]/4.
 *  Only physical tests (positive jacobian) are done (i.e. approximation of the surface,
 *  etc... must be performed outside). */
int _MMG5_colver(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist,char indq,char typchk) {
  MMG5_pTetra          pt,pt1;
  MMG5_pxTetra         pxt,pxt1;
  MMG5_xTetra          xt,xts;
  int             i,iel,jel,pel,qel,k,np,nq,*adja,p0,p1;
  unsigned char   ip,iq,j,voy,voyp,voyq,ia,iav;
  unsigned char   (*ind)[2];
  int             *p0_c,*p1_c;
  char            indar[4][4][2] = {
    /* indar[ip][iq][0/1]: indices of edges which have iq for extremity but not ip*/
    { {-1,-1}, { 3, 4}, { 3, 5}, { 4, 5} },
    { { 1, 2}, {-1,-1}, { 1, 5}, { 2, 5} },
    { { 0, 2}, { 0, 4}, {-1,-1}, { 2, 4} },
    { { 0, 1}, { 0, 3}, { 1, 3}, {-1,-1} } };

  // Dynamic allocations for windows compatibility
  if (!(ind = malloc(ilist * sizeof(unsigned char[2])))) {
	  perror("  ## Memory problem: malloc");
	  exit(EXIT_FAILURE);
  }
  _MMG5_SAFE_CALLOC(p0_c, ilist, int);
  _MMG5_SAFE_CALLOC(p1_c, ilist, int);

  iel = list[0] / 4;
  ip  = list[0] % 4;
  pt  = &mesh->tetra[iel];
  np  = pt->v[ip];
  nq  = pt->v[indq];

  /* Mark elements of the shell of edge (pq) */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 4;
    i   = list[k] % 4;
    pt  = &mesh->tetra[iel];

    for (j=0; j<3; j++) {
      i = _MMG5_inxt3[i];
      if ( pt->v[i] == nq ) {
        /* list edges that we need to update */
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          ip  = list[k]%4;
          ind[k][0] = indar[ip][i][0];
          if ( pxt->tag[ind[k][0]] || pxt->edg[ind[k][0]] ) {
            if ( _MMG5_iare[ind[k][0]][0]==i )  p0_c[k] = pt->v[_MMG5_iare[ind[k][0]][1]];
            else  p0_c[k] = pt->v[_MMG5_iare[ind[k][0]][0]];
          }
          ind[k][1] = indar[ip][i][1];
          if ( pxt->tag[ind[k][1]] || pxt->edg[ind[k][1]] ) {
            if ( _MMG5_iare[ind[k][1]][0]==i )  p1_c[k] = pt->v[_MMG5_iare[ind[k][1]][1]];
            else  p1_c[k] = pt->v[_MMG5_iare[ind[k][1]][0]];
          }
        }
        list[k] *= -1;
        break;
      }
    }
  }

  /* avoid recreating existing elt */
  for (k=0; k<ilist; k++) {
    if ( list[k] < 0 )  continue;
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];

    /* update edges of elements that do not belong to the shell of pq */
    if ( !pt->xt ) {
      continue;
    }
    pxt = &mesh->xtetra[pt->xt];
    for ( i=0; i<ilist; i++ ) {
      if ( (list[i]>0) || (!(mesh->tetra[-list[i]/4].xt)) )  continue;
      pt1  = &mesh->tetra[-list[i]/4];
      pxt1 = &mesh->xtetra[pt1->xt];
      if ( p0_c[i] ) {
        for ( j=0; j<3; j++) {
          ia = _MMG5_idir[ip][j];
          if ( pt->v[ia]==p0_c[i] ) {
            pxt->tag[_MMG5_arpt[ip][j]] |= pxt1->tag[ind[i][0]];
            if ( !pxt->edg[_MMG5_arpt[ip][j]] )
              pxt->edg[_MMG5_arpt[ip][j]] = pxt1->edg[ind[i][0]];
            else if ( pxt1->edg[ind[i][0]] )
              pxt->edg[_MMG5_arpt[ip][j]] =
                MG_MAX(pxt->edg[_MMG5_arpt[ip][j]],pxt1->edg[ind[i][0]]);
            break;
          }
        }
      }
      if ( p1_c[i] ) {
        for ( j=0; j<3; j++) {
          ia = _MMG5_idir[ip][j];
          if ( pt->v[ia]==p1_c[i] ) {
            pxt->tag[_MMG5_arpt[ip][j]] |= pxt1->tag[ind[i][1]];
            if ( !pxt->edg[_MMG5_arpt[ip][j]] )
              pxt->edg[_MMG5_arpt[ip][j]] = pxt1->edg[ind[i][1]];
            else if ( pxt1->edg[ind[i][1]] )
              pxt->edg[_MMG5_arpt[ip][j]] =
                MG_MAX(pxt->edg[_MMG5_arpt[ip][j]],pxt1->edg[ind[i][1]]);
            break;
          }
        }
      }
    }
    adja = &mesh->adja[4*(iel-1)+1];
    jel  = adja[ip];
    if ( !jel )  continue;

    jel /= 4;
    voy  = adja[ip] % 4;
    pt = &mesh->tetra[jel];
    if (pt->v[voy] == nq) {
      _MMG5_SAFE_FREE(ind); _MMG5_SAFE_FREE(p0_c); _MMG5_SAFE_FREE(p1_c);
      return(0);
    }
  }

  /* deal with the shell of edge (pq) and the implied updates */
  for (k=0; k<ilist; k++) {
    if ( list[k] > 0 )  continue;
    iel = (-list[k]) / 4;
    ip  = (-list[k]) % 4;
    pt  = &mesh->tetra[iel];

    iq  = ip;
    for (j=0; j<3; j++) {
      iq = _MMG5_inxt3[iq];
      if ( pt->v[iq] == nq )  break;
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

          if ( qel && (pt1->ref < mesh->tetra[qel].ref) )  MG_CLR( pxt1->ori,voyp );
          else   MG_SET(pxt1->ori,voyp);


          /* update tags for edges */
          for ( j=0; j<3; j++ ) {
            ia = _MMG5_iarf[ip][j];
            p0 = pt->v[_MMG5_iare[ia][0]];
            p1 = pt->v[_MMG5_iare[ia][1]];

            for ( i=0; i<3; i++ ) {
              iav=_MMG5_iarf[voyp][i];
              if ( p0==nq ) {
                if ( ((pt1->v[_MMG5_iare[iav][0]]==np) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                     ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==np)) )
                  break;
              }
              else if ( p1==nq ) {
                if ( ((pt1->v[_MMG5_iare[iav][0]]==np) && (pt1->v[_MMG5_iare[iav][1]]==p0)) ||
                     ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==np)) )
                  break;
              }
              else {
                if ( ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                     ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==p0)) )
                  break;
              }
            }
            assert(i!=3);
            pxt1->tag[iav] = pxt1->tag[iav] | pxt->tag[ia];
            pxt1->edg[iav] = MG_MAX(pxt1->edg[iav],pxt->edg[ia]);
          }
        }
        else {
          pxt1 = &xt;
          memset(pxt1,0,sizeof(MMG5_xTetra));
          pxt1->ref[voyp] = pxt->ref[ip];
          pxt1->ftag[voyp] = pxt->ftag[ip];
          pxt1->ori = 15;
          if ( !MG_GET(pxt->ori,ip) )  MG_CLR(pxt1->ori,voyp);

          /* update tags for edges */
          for ( j=0; j<3; j++ ) {
            ia = _MMG5_iarf[ip][j];
            p0 = pt->v[_MMG5_iare[ia][0]];
            p1 = pt->v[_MMG5_iare[ia][1]];
            if ( pxt->tag[ia] ) {
              for ( i=0; i<3; i++ ) {
                iav=_MMG5_iarf[voyp][i];
                if ( p0==nq ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==np) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==np)) )
                    break;
                }
                else if ( p1==nq ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==np ) && (pt1->v[_MMG5_iare[iav][1]]==p0)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==np )) )
                    break;
                }
                else {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt->tag[ia];
              pxt1->edg[iav] = pxt->edg[ia];
            }
          }
          /* Recover the already used place by pxt */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(MMG5_xTetra));
        }
      }
      else {
        /* Only the values corresponding to pt become 0 */
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

            if ( pel && (pt1->ref < mesh->tetra[pel].ref) )  MG_CLR( pxt1->ori,voyq );
            else   MG_SET(pxt1->ori,voyq);

            /* update tags for edges */
            for ( j=0; j<3; j++ ) {
              ia = _MMG5_iarf[iq][j];
              p0 = pt->v[_MMG5_iare[ia][0]];
              p1 = pt->v[_MMG5_iare[ia][1]];
              for ( i=0; i<3; i++ ) {
                iav=_MMG5_iarf[voyq][i];
                if ( p0==np ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==nq) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==nq)) )
                    break;
                }
                else if ( p1==np ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==nq ) && (pt1->v[_MMG5_iare[iav][1]]==p0)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==nq )) )
                    break;
                }
                else {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt1->tag[iav] | pxt->tag[ia];
              pxt1->edg[iav] = MG_MAX(pxt1->edg[iav],pxt->edg[ia]);
            }
          }
          else {
            pxt1 = &xt;
            memset(pxt1,0,sizeof(MMG5_xTetra));
            pxt1->ref[voyq] = pxt->ref[iq];
            pxt1->ftag[voyq] = pxt->ftag[iq];
            pxt1->ori = 15;
            if ( !MG_GET(pxt->ori,iq) )  MG_CLR(pxt1->ori,voyq);
            /* update tags for edges */
            for ( j=0; j<3; j++ ) {
              ia = _MMG5_iarf[iq][j];
              p0 = pt->v[_MMG5_iare[ia][0]];
              p1 = pt->v[_MMG5_iare[ia][1]];
              if ( pxt->tag[ia] ) {
                for ( i=0; i<3; i++ ) {
                  iav=_MMG5_iarf[voyq][i];
                  if ( p0==np ) {
                    if ( ((pt1->v[_MMG5_iare[iav][0]]==nq) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                         ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==nq)) )
                      break;
                  }
                  else if ( p1==np ) {
                    if ( ((pt1->v[_MMG5_iare[iav][0]]==nq ) && (pt1->v[_MMG5_iare[iav][1]]==p0)) ||
                         ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==nq )) )
                      break;
                  }
                  else {
                    if ( ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                         ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==p0)) )
                      break;
                  }
                }
                assert(i!=3);
                pxt1->tag[iav] = pxt->tag[ia];
                pxt1->edg[iav] = pxt->edg[ia];
              }
            }
            /* Create new field xt */
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
			                     _MMG5_SAFE_FREE(ind); _MMG5_SAFE_FREE(p0_c); _MMG5_SAFE_FREE(p1_c);
                                 return(-1));
            }
            pt1->xt = mesh->xt;
            pxt = &mesh->xtetra[pt1->xt];
            memcpy(pxt,pxt1,sizeof(MMG5_xTetra));
          }
        }
        else {
          /* Only the values corresponding to pt become 0 */
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
          for ( j=0; j<3; j++ ) {
            ia = _MMG5_iarf[iq][j];
            p0 = pt->v[_MMG5_iare[ia][0]];
            p1 = pt->v[_MMG5_iare[ia][1]];
            if ( pxt->tag[ia] ) {
              for ( i=0; i<3; i++ ) {
                iav=_MMG5_iarf[voyq][i];
                if ( p0==np ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==nq) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==nq)) )
                    break;
                }
                else if ( p1==np ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==nq ) && (pt1->v[_MMG5_iare[iav][1]]==p0)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==nq )) )
                    break;
                }
                else {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt->tag[ia];
              pxt1->edg[iav] = pxt->edg[ia];
            }
          }
        }
        else {
          pxt1 = &xt;
          memset(pxt1,0,sizeof(MMG5_xTetra));
          pxt1->ref[voyq]  = pxt->ref[iq];
          pxt1->ftag[voyq] = pxt->ftag[iq];
          pxt1->ori = 15;

          /* update tags for edges */
          for ( j=0; j<3; j++ ) {
            ia = _MMG5_iarf[iq][j];
            p0 = pt->v[_MMG5_iare[ia][0]];
            p1 = pt->v[_MMG5_iare[ia][1]];
            if ( pxt->tag[ia] ) {
              for ( i=0; i<3; i++ ) {
                iav = _MMG5_iarf[voyq][i];
                if ( p0==np ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==nq) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==nq)) )
                    break;
                }
                else if ( p1==np ) {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==nq ) && (pt1->v[_MMG5_iare[iav][1]]==p0)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==nq )) )
                    break;
                }
                else {
                  if ( ((pt1->v[_MMG5_iare[iav][0]]==p0) && (pt1->v[_MMG5_iare[iav][1]]==p1)) ||
                       ((pt1->v[_MMG5_iare[iav][0]]==p1) && (pt1->v[_MMG5_iare[iav][1]]==p0)) )
                    break;
                }
              }
              assert(i!=3);
              pxt1->tag[iav] = pxt->tag[ia];
              pxt1->edg[iav] = pxt->edg[ia];
            }
          }
          /* Recover the already used place by pxt */
          pt1->xt = pt->xt;
          memcpy(pxt,pxt1,sizeof(MMG5_xTetra));
        }
      }
    }
    _MMG3D_delElt(mesh,iel);
  }

  /* Update vertices coordinates for elements that do not belong to the shell of (pq) */
  for (k=0; k<ilist;  k++) {
    if ( list[k] < 0 )  continue;
    iel = list[k] / 4;
    ip  = list[k] % 4;
    pt  = &mesh->tetra[iel];
    pt->v[ip] = nq;
    if ( typchk==1 && met->m && met->size > 1 )
      pt->qual=_MMG5_caltet33_ani(mesh,met,pt);
    else
      pt->qual=_MMG5_orcal(mesh,met,iel);
    pt->mark=mesh->mark;
  }

  _MMG5_SAFE_FREE(ind); _MMG5_SAFE_FREE(p0_c); _MMG5_SAFE_FREE(p1_c);
  return(np);
}
