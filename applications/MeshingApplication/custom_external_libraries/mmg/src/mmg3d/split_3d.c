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
 * \file mmg3d/split_3d.c
 * \brief Functions to create new points.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "inlined_functions_3d.h"

extern char  ddb;

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \return 0 if split leads to invalid situation, else 1.
 *
 * Simulate the splitting of 1 edge of element
 *
 */
int MMG3D_split1_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4];
  const unsigned char *taued;

  /* tau = sigma^-1 = permutation that sends the reference config (edge 01 split) to the current */
  pt = &mesh->tetra[k];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK )  return 0;

  pt0 = &mesh->tetra[0];

  /* default is case 1 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(pt->flag) {
  case 2:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &MMG5_permedge[6][0];
    break;
  case 4:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  case 8:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;
  case 16:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  case 32:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &MMG5_permedge[11][0];
    break;
  }

  /* Test volume of the two created tets */
  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[1]] = vx[taued[0]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[0]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split 1 edge of tetra \a k.
 *
 */
int MMG5_split1(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra         pt,pt1;
  MMG5_xTetra         xt,xt1;
  MMG5_pxTetra        pxt0;
  int                 iel;
  char                i,isxt,isxt1;
  unsigned char       tau[4];
  const unsigned char *taued;

  /* create a new tetra */
  pt  = &mesh->tetra[k];
  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt = &mesh->tetra[k];
  }

  pt1 = &mesh->tetra[iel];
  pt1 = memcpy(pt1,pt,sizeof(MMG5_Tetra));
  pxt0 = NULL;
  if ( pt->xt ) {
    pxt0 = &mesh->xtetra[pt->xt];
    memcpy(&xt,pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt1,pxt0,sizeof(MMG5_xTetra));
  }
  else {
    memset(&xt,0,sizeof(MMG5_xTetra));
    memset(&xt1,0,sizeof(MMG5_xTetra));
  }

  /* default is case 1 */
  tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(pt->flag) {
  case 2:
    tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
    taued = &MMG5_permedge[6][0];
    break;
  case 4:
    tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  case 8:
    tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;
  case 16:
    tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  case 32:
    tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
    taued = &MMG5_permedge[11][0];
    break;
  }

  /* Generic formulation of split of 1 edge */
  pt->v[tau[1]] = pt1->v[tau[0]] = vx[taued[0]];

  if ( pt->xt ) {
    /* Reset edge tag */
    xt.tag [taued[3]] = 0;  xt.tag [taued[4]] = 0;
    xt1.tag[taued[1]] = 0;  xt1.tag[taued[2]] = 0;
    xt.edg [taued[3]] = 0;  xt.edg [taued[4]] = 0;
    xt1.edg[taued[1]] = 0;  xt1.edg[taued[2]] = 0;
    xt.ref [  tau[0]] = 0;  xt.ftag [ tau[0]] = 0;  MG_SET( xt.ori, tau[0]);
    xt1.ref[  tau[1]] = 0;  xt1.ftag[ tau[1]] = 0;  MG_SET(xt1.ori, tau[1]);
  }

  pt->flag = pt1->flag = 0;
  isxt  = 0 ;
  isxt1 = 0;
  for (i=0; i<4; i++) {
    if ( xt.ref[i]  || xt.ftag[i] )  isxt = 1;
    if ( xt1.ref[i] || xt1.ftag[i] ) isxt1 = 1;
    if ( isxt && isxt1 )  goto nextstep1;
  }

nextstep1:
  if ( pt->xt ) {
    if ( isxt && !isxt1 ) {
      pt1->xt = 0;
      memcpy(pxt0,&xt,sizeof(MMG5_xTetra));
    }
    else if ( !isxt && isxt1 ) {
      pt1->xt = pt->xt;
      pt->xt  = 0;
      pxt0 = &mesh->xtetra[pt1->xt];
      memcpy(pxt0,&xt1,sizeof(MMG5_xTetra));
    }
    else if ( isxt && isxt1 ) {
      mesh->xt++;
      if ( mesh->xt > mesh->xtmax ) {
        /* realloc of xtetras table */
        MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                           "larger xtetra table",
                           mesh->xt--;
                           fprintf(stderr,"  Exit program.\n");
                           return 0);
        pxt0 = &mesh->xtetra[pt->xt];
      }
      pt1->xt = mesh->xt;

      assert ( pxt0 );
      memcpy(pxt0,&xt,sizeof(MMG5_xTetra));
      pxt0 = &mesh->xtetra[mesh->xt];
      memcpy(pxt0,&xt1,sizeof(MMG5_xTetra));
    }
    else {
      pt->xt = 0;
      pt1->xt = 0;
    }
  }
  /* Quality update */
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    pt->qual=MMG5_caltet33_ani(mesh,met,pt);
    pt1->qual=MMG5_caltet33_ani(mesh,met,pt1);
  }
  else
  {
    pt->qual=MMG5_orcal(mesh,met,k);
    pt1->qual=MMG5_orcal(mesh,met,iel);
  }
  return 1;
}
/**
 * \param mesh  pointer toward the mesh structure
 * \param start index of the tetra that we want to split
 * \param iface local index of the boundary face that we want to split
 * \param ia    local index of the boundary edge that we want to split
 * \param idx   local index of the new tetra that we want to study after the
 *              splitting of the tetra \a start (idx=0 or 1)
 * \param ip    new point index
 * \param n0    normal of the new boundary face in the tetra idx.
 *
 * \return 1 if success (no new sharp angle), 0 if we create a sharp angle,
 *         -1 if fail.
 *
 * Check that the split of the edge \a ia of the tetra \a start does not create
 * a ridge along the \f$ idx^{th} \f$ edge opposite to \a ip in the boundary
 * triangle \a iface. Store the normal of the \f$ idx^{th} \f$ boundary triangle
 * in \a n0.
 *
 */
static inline
int MMG3D_normalDeviation(MMG5_pMesh mesh , int  start, char   iface, char ia,
                           int        idx  , int  ip   , double n0[3])
{
  MMG5_Tria tt0;
  double    n1[3];
  int       iedge,iploc,ier;

  /** Store the first boundary triangle (the one that is created in the boundary
   * face that we split) */
  MMG5_tet2tri(mesh,start,iface,&tt0);

  iedge = MMG5_iarfinv[iface][ia];

  switch (idx)
  {
  case 0:
    iploc = MMG5_iprv2[iedge];
    break;
  case 1:
    iploc = MMG5_inxt2[iedge];
    break;
  }

  tt0.v[iploc] = ip;

  /** Compute the normal of the first triangle */
  if ( !MMG5_nortri(mesh, &tt0, n0) ) return -1;

  if ( tt0.tag[iploc] & MG_GEO || tt0.tag[iploc] & MG_NOM ) return 1;

  /** Compute the normal of the second triangle (triangle adjacent to the first
   * through the edge iploc) */
  ier = MMG3D_normalAdjaTri(mesh,start,iface,iploc,n1);
  if ( ier<0 ) return -1;
  else if ( !ier ) return 0;

  ier =  MMG5_devangle( n0, n1, mesh->info.dhd );

  return  ier;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric.
 * \param list pointer toward the edge shell.
 * \param ret size of the edge shell.
 * \param ip new point index.
 *
 * \return 1 if all checks are ok
 * \return 0 if fail due to a very bad quality elt
 * \return 2 if fail due to a ridge angle creation
 *
 * Simulate at the same time creation and bulging of one point, with new
 * position o and tag \a tag, to be inserted at an edge, whose shell is passed.
 *
 */
int MMG3D_simbulgept(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ret,int ip) {
  MMG5_pTetra    pt,pt0;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    ppt0;
  double         calold,calnew,caltmp;
  double         n0[6],n1[6];
  int            j,k,iel,ilist,idx,iface,ier,sum1,sum2,mins1,mins2,maxs1,maxs2;
  int            is0,is1,is2;
  char           ie,ia,ib,complete,wrongOri;

  ilist = ret / 2;
  pt0  = &mesh->tetra[0];
  ppt0 = &mesh->point[0];
  memcpy(ppt0->c  ,&mesh->point[ip].c  , 3*sizeof(double));
  ppt0->tag = mesh->point[ip].tag;

  memcpy(&met->m[0],&met->m[met->size*ip], met->size*sizeof(double));

  calold = calnew = DBL_MAX;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    ie  = list[k] % 6;
    ia = MMG5_iare[ie][0];
    ib = MMG5_iare[ie][1];

    pt = &mesh->tetra[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ia] = 0;
    calold = MG_MIN(calold,pt->qual);
    caltmp = MMG5_orcal(mesh,met,0);
    if ( caltmp < MMG5_EPSOK )  return 0;
    calnew = MG_MIN(calnew,caltmp);

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ib] = 0;
    caltmp = MMG5_orcal(mesh,met,0);
    if ( caltmp < MMG5_EPSOK )  return 0;
    calnew = MG_MIN(calnew,caltmp);
  }
  if ( calnew <= MMG5_EPSOK ) {
    return 0;
  }

  /* if ( calnew < 0.3*calold )  return 0;*/

  /** Check the deviation for new triangles */
  /* analyze surfacic ball of p */
  wrongOri = complete = idx = 0;
  maxs1 = mins1 = sum1 = 0;
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    ie  = list[k] % 6;

    pt   = &mesh->tetra[iel];
    if(!pt->xt) continue;

    pxt  = &mesh->xtetra[pt->xt];

    for ( j=0; j<2; ++j ) {
      iface = MMG5_ifar[ie][j];
      if ( !(pxt->ftag[iface] & MG_BDY) ) continue;

      is0 = pt->v[MMG5_idir[iface][0]];
      is1 = pt->v[MMG5_idir[iface][1]];
      is2 = pt->v[MMG5_idir[iface][2]];

      /* Normal deviation between the two new triangles and their neighbors */
      ier = MMG3D_normalDeviation(mesh,iel,iface,ie,0,ip,&n0[idx]);
      if ( ier < 0 ) return -1;
      else if ( ier == 0 ) return 2;

      ier = MMG3D_normalDeviation(mesh,iel,iface,ie,1,ip,&n1[idx]);
      if ( ier < 0 ) return -1;
      else if ( ier == 0 ) return 2;

      /* Test sharp angle creation along the new edge */
      if ( !MMG5_devangle(&n0[idx],&n1[idx],mesh->info.dhd) ) {
        return 2;
      }

      if ( !idx ) {
        sum1 = is0 + is1 + is2;
        mins1 = MG_MIN(is0, MG_MIN(is1,is2));
        maxs1 = MG_MAX(is0, MG_MAX(is1,is2));
        idx = 3;
      }
      else {
        /* don't check if it is a ridge edge or if we have already cross 2
         * boundaries */
        if ( complete || pxt->tag[ie] & MG_GEO || pxt->tag[ie] & MG_NOM )
          continue;

        /* We are manifold thus we have exactly two faces in our shell: check
         * that we don't see twice a boundary face (multidomain case) */
        sum2 = is0 + is1 + is2;
        mins2 = MG_MIN(is0, MG_MIN(is1,is2));
        maxs2 = MG_MAX(is0, MG_MAX(is1,is2));

        if ( (sum2 == sum1 && mins2 == mins1 && maxs2 == maxs1) ) {
          /* Multidomain: we see the tria for the second time (and from another
           * side), this means that the next seen tria will be seen from the
           * wrong side. */
          wrongOri = 1;
          continue;
        }
        else if ( wrongOri ) {
          /* We skeep this tria because it is seen from a wrong side. The next
           * will be ok */
          wrongOri = 0;
          continue;
        }
        complete = 1;

        /* Test sharp angle creation along the splitted edge */
        if ( !MMG5_devangle(&n0[0],&n1[idx],mesh->info.dhd) ) {
          return 2;
        }
        if ( !MMG5_devangle(&n1[0],&n0[idx],mesh->info.dhd) ) {
          return 2;
        }
      }
    }
  }

  return 1;
}

/**
 * \param mesh  pointer toward the mesh structure
 * \param start index of the working tetra
 * \param iface local index of the boundary face of the tetra \a start
 * \param ia    local index on face \a iface of the edge through which we seek
 *              the adjacent triangle of the triangle \a iface of \a start.
 * \param n     normal of the new boundary face in the tetra idx.
 *
 * \return 1 if success, 0 if we want to refuse the collapse, -1 if fail.
 *
 * Compute the normal of the adjacent triangle of the triangle \a iface of the
 * tetra \a start through the edge \a ia (in local numbering of the face).
 *
 */
int MMG3D_normalAdjaTri(MMG5_pMesh mesh , int start, char iface, int ia,
                         double n[3]                                     )
{
  MMG5_Tria tt;
  int       iedgeOpp,list[MMG3D_LMAX+2],it1,it2,it;

  iedgeOpp = MMG5_iarf[iface][ia];

  /** Store the adjacent boundary triangle (triangle adjacent to \a iface
   * through the edge ia */
  if ( MMG5_coquilface( mesh,start,iface,iedgeOpp,list,&it1,&it2,0) <= 0 )
    return -1;

  if ( it1/4 != start || it1%4 != iface ) {
    //assert ( it2/4==start && it2%4==iface );
    if ( it2/4!=start || it2%4!=iface ) return 0;

    it = it1;
  }
  else {
    it = it2;
  }
  MMG5_tet2tri(mesh,it/4,it%4,&tt);

  /** Compute the normal of the second triangle */
  if ( !MMG5_nortri(mesh, &tt, n) ) return 0;

  return 1;
}



/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param list pointer toward the shell of edge.
 * \param ret size of the shell of edge.
 * \param ip idex of new point.
 * \param cas flag to watch the length of the new edges.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage,
 * 1 for special storage.
 * \param chkRidTet if 1, avoid the creation of a tet with 4 ridge vertices
 * \return -1 if we fail, 0 if we don't split the edge, 1 if success.
 *
 * Split edge \f$list[0]\%6\f$, whose shell list is passed, introducing point \a
 * ip Beware : shell has to be enumerated in ONLY ONE TRAVEL (always same
 * sense).
 *
 */
int MMG5_split1b(MMG5_pMesh mesh, MMG5_pSol met,int *list, int ret, int ip,
                  int cas,char metRidTyp,char chkRidTet){
  MMG5_pTetra          pt,pt1,pt0;
  MMG5_xTetra          xt,xt1;
  MMG5_pxTetra         pxt0;
  double               lmin,lmax,len;
  int                  ilist,k,open,iel,jel,*newtet,nump,*adja,j;
  int                 *adjan,nei2,nei3,mel;
  char                 ie,tau[4],isxt,isxt1,i,voy;
  const unsigned char *taued;

  ilist = ret / 2;
  open  = ret % 2;


  if ( cas && met->m ) {
    lmin = 0.6;
    lmax = 1.3;
    for (j=0; j<ilist; j++) {
      for (i=0; i<6; i++) {
        pt   = &mesh->tetra[list[j]/6];
        if ( (!metRidTyp) && met->m && met->size>1 )
          len = MMG5_lenedg33_ani(mesh,met,i,pt);
        else
          len  = MMG5_lenedg(mesh,met,i,pt);
        if ( len < lmin) {
          lmin = len;
        }
        else if ( len > lmax) {
          lmax = len;
        }
      }
    }
    assert( lmin!=0 );

    /** 2 different checks :
        1) are we creating a too small edge  (BUG_Split1b_SpereIso_0.125h_met)
        2) in aniso and from the last wave of anatet(typchk=1): avoid the
        creation of a tetra with 4 ridge vertices.
    **/
    for (j=0; j<ilist; j++) {
      iel = list[j] / 6;
      pt  = &mesh->tetra[iel];
      ie  = list[j] % 6;
      pt0 = &mesh->tetra[0];
      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      /* tau = sigma^-1 = permutation that sends the ref config (edge 01 split) to current */
      tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
      taued = &MMG5_permedge[0][0];
      switch(ie){
      case 1:
        tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
        taued = &MMG5_permedge[6][0];
        break;
      case 2:
        tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
        taued = &MMG5_permedge[2][0];
        break;
      case 3:
        tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
        taued = &MMG5_permedge[4][0];
        break;
      case 4:
        tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
        taued = &MMG5_permedge[10][0];
        break;
      case 5:
        tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
        taued = &MMG5_permedge[11][0];
        break;
      }

      pt0->v[MMG5_isar[ie][1]] = ip;
      if ( chkRidTet ) {
        if ( !MMG3D_chk4ridVertices(mesh,pt0) ) {
          break;
        }
      }
      if ( (!metRidTyp) && met->m && met->size>1 )
        len = MMG5_lenedgspl33_ani(mesh,met,taued[5],pt0);
      else
        len = MMG5_lenedgspl(mesh,met,taued[5],pt0);
      if ( len < lmin )  break;
      memcpy(pt0,pt,sizeof(MMG5_Tetra));

      pt0->v[MMG5_isar[ie][0]] = ip;
      if ( chkRidTet ) {
        if ( !MMG3D_chk4ridVertices(mesh,pt0) ) {
          break;
        }
      }

      if ( (!metRidTyp) && met->m && met->size>1 )
        len = MMG5_lenedgspl33_ani(mesh,met,taued[5],pt0);
      else
        len = MMG5_lenedgspl(mesh,met,taued[5],pt0);
      if ( len < lmin )  break;
    }
    if ( j < ilist )  return 0;
  }

  MMG5_SAFE_CALLOC(newtet,ilist,int,return -1);

  iel = list[0] / 6;
  ie  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  nump = pt->v[MMG5_iare[ie][0]];

  /* Fill list newtet[k] = +_created tetra for list[k]/6 : + if kept tetra (= one associated to
     pt->v[tau[0]]) is associated with nump, - if with numq */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    ie  = list[k] % 6;
    pt  = &mesh->tetra[iel];
    /* identity : case 0 */
    tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
    switch(ie) {
    case 1:
      tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
      break;
    case 2:
      tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
      break;
    case 3:
      tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
      break;
    case 4:
      tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
      break;
    case 5:
      tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
      break;
    }
    jel = MMG3D_newElt(mesh);
    if ( !jel ) {
      MMG3D_TETRA_REALLOC(mesh,jel,mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                  " a new element.\n",__func__);
                          MMG5_INCREASE_MEM_MESSAGE();
                          k--;
                          for ( ; k>=0 ; --k ) {
                            if ( !MMG3D_delElt(mesh,abs(newtet[k])) ) return -1;
                          }
                          return -1);
      pt  = &mesh->tetra[iel];
    }
    pt1 = &mesh->tetra[jel];
    memcpy(pt1,pt,sizeof(MMG5_Tetra));

    if ( pt->v[tau[0]] == nump )
      newtet[k] = jel;
    else
      newtet[k] = -jel;
  }

  /* Special case : only one element in the shell */
  if ( ilist == 1 ) {
    assert(open);
    iel = list[0] / 6;
    ie  = list[0] % 6;
    pt = &mesh->tetra[iel];
    jel = abs(newtet[0]);
    pt1 = &mesh->tetra[jel];

    pxt0 = 0;
    if ( pt->xt ) {
      pxt0 = &mesh->xtetra[pt->xt];
      memcpy(&xt,pxt0,sizeof(MMG5_xTetra));
      memcpy(&xt1,pxt0,sizeof(MMG5_xTetra));
    }
    else {
      memset(&xt,0, sizeof(MMG5_xTetra));
      memset(&xt1,0, sizeof(MMG5_xTetra));
    }

    /* tau = sigma^-1 = permutation that sends the reference config (edge 01 split) to current */
    tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
    taued = &MMG5_permedge[0][0];
    switch(ie){
    case 1:
      tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
      taued = &MMG5_permedge[6][0];
      break;
    case 2:
      tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
      taued = &MMG5_permedge[2][0];
      break;
    case 3:
      tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
      taued = &MMG5_permedge[4][0];
      break;
    case 4:
      tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
      taued = &MMG5_permedge[10][0];
      break;
    case 5:
      tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
      taued = &MMG5_permedge[11][0];
      break;
    }

    /* Generic formulation of split of 1 edge */
    pt->v[tau[1]] = pt1->v[tau[0]] = ip;
    if ( pt->xt ) {
      /* Reset edge tag */
      xt.tag [taued[3]] = 0;  xt.tag [taued[4]] = 0;
      xt1.tag[taued[1]] = 0;  xt1.tag[taued[2]] = 0;
      xt.edg [taued[3]] = 0;  xt.edg [taued[4]] = 0;
      xt1.edg[taued[1]] = 0;  xt1.edg[taued[2]] = 0;
      xt.ref [  tau[0]] = 0;  xt.ftag [ tau[0]] = 0;  MG_SET( xt.ori, tau[0]);
      xt1.ref[  tau[1]] = 0;  xt1.ftag[ tau[1]] = 0;  MG_SET(xt1.ori, tau[1]);
    }

    pt->flag = pt1->flag = 0;

    isxt = 0 ;
    isxt1 = 0;

    for (i=0; i<4; i++) {
      if ( xt.ref[i]  || xt.ftag[i] )  isxt = 1;
      if ( xt1.ref[i] || xt1.ftag[i])  isxt1 = 1;
    }

    if ( pt->xt ) {
      if ( (isxt) && (!isxt1) ) {
        pt1->xt = 0;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(MMG5_xTetra));
      }
      else if ( (!isxt) && (isxt1) ) {
        pt1->xt = pt->xt;
        pt->xt = 0;
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(MMG5_xTetra));
      }
      else if ( isxt && isxt1 ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             return -1);
        }
        pt1->xt = mesh->xt;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(MMG5_xTetra));
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(MMG5_xTetra));
      }
      else {
        pt->xt = 0;
        pt1->xt = 0;
      }
    }

    /* Update of adjacency relations */
    adja = &mesh->adja[4*(iel-1)+1];
    adjan = &mesh->adja[4*(jel-1)+1];

    adja[tau[2]] = 0;
    adja[tau[3]] = 0;
    adjan[tau[2]] = 0;
    adjan[tau[3]] = 0;

    mel = adja[tau[0]] / 4;
    voy = adja[tau[0]] % 4;
    adja[tau[0]] = 4*jel + tau[1];
    adjan[tau[0]] = 4*mel + voy;
    adjan[tau[1]] = 4*iel + tau[0];

    if ( mel ) {
      adjan = &mesh->adja[4*(mel -1) +1];
      adjan[voy] = 4*jel + tau[0];
    }
    /* Quality update */
    if ( (!metRidTyp) && met->m && met->size>1 ) {
      pt->qual=MMG5_caltet33_ani(mesh,met,pt);
      pt1->qual=MMG5_caltet33_ani(mesh,met,pt1);
    }
    else {
      pt->qual=MMG5_orcal(mesh,met,iel);
      pt1->qual=MMG5_orcal(mesh,met,jel);
    }
    pt->mark  = mesh->mark;
    pt1->mark = mesh->mark;


    MMG5_SAFE_FREE(newtet);
    return 1;
  }

  /* General case : update each element of the shell */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 6;
    ie  = list[k] % 6;
    pt = &mesh->tetra[iel];
    jel = abs(newtet[k]);
    pt1 = &mesh->tetra[jel];

    pxt0 = 0;
    if ( pt->xt ) {
      pxt0 = &mesh->xtetra[pt->xt];
      memcpy(&xt,pxt0,sizeof(MMG5_xTetra));
      memcpy(&xt1,pxt0,sizeof(MMG5_xTetra));
    }
    else {
      memset(&xt,0, sizeof(MMG5_xTetra));
      memset(&xt1,0, sizeof(MMG5_xTetra));
    }

    /* tau = sigma^-1 = permutation that sends the reference config (edge 01 split) to current */
    tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
    taued = &MMG5_permedge[0][0];
    switch(ie){
    case 1:
      tau[0] = 2; tau[1] = 0; tau[2] = 1; tau[3] = 3;
      taued = &MMG5_permedge[6][0];
      break;
    case 2:
      tau[0] = 0; tau[1] = 3; tau[2] = 1; tau[3] = 2;
      taued = &MMG5_permedge[2][0];
      break;
    case 3:
      tau[0] = 1; tau[1] = 2; tau[2] = 0; tau[3] = 3;
      taued = &MMG5_permedge[4][0];
      break;
    case 4:
      tau[0] = 3; tau[1] = 1; tau[2] = 0; tau[3] = 2;
      taued = &MMG5_permedge[10][0];
      break;
    case 5:
      tau[0] = 3; tau[1] = 2; tau[2] = 1; tau[3] = 0;
      taued = &MMG5_permedge[11][0];
      break;
    }

    /* Generic formulation of split of 1 edge */
    pt->v[tau[1]] = pt1->v[tau[0]] = ip;
    if ( pt->xt ) {
      /* Reset edge tag */
      xt.tag [taued[3]] = 0;  xt.tag [taued[4]] = 0;
      xt1.tag[taued[1]] = 0;  xt1.tag[taued[2]] = 0;
      xt.edg [taued[3]] = 0;  xt.edg [taued[4]] = 0;
      xt1.edg[taued[1]] = 0;  xt1.edg[taued[2]] = 0;
      xt.ref [  tau[0]] = 0;  xt.ftag [ tau[0]] = 0;  MG_SET( xt.ori, tau[0]);
      xt1.ref[  tau[1]] = 0;  xt1.ftag[ tau[1]] = 0;  MG_SET(xt1.ori, tau[1]);
    }

    pt->flag = pt1->flag = 0;

    isxt = 0 ;
    isxt1 = 0;

    for (i=0; i<4; i++) {
      if ( xt.ref[i]  || xt.ftag[i]  )  isxt  = 1;
      if ( xt1.ref[i] || xt1.ftag[i] )  isxt1 = 1;
    }

    if ( pt->xt ) {
      if ( (isxt)&&(!isxt1) ) {
        pt1->xt = 0;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(MMG5_xTetra));
      }
      else if ((!isxt)&&(isxt1) ) {
        pt1->xt = pt->xt;
        pt->xt = 0;
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(MMG5_xTetra));
      }
      else if (isxt && isxt1 ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             return -1);
        }
        pt1->xt = mesh->xt;
        pxt0 = &mesh->xtetra[pt->xt];
        memcpy(pxt0,&xt,sizeof(MMG5_xTetra));
        pxt0 = &mesh->xtetra[pt1->xt];
        memcpy(pxt0,&xt1,sizeof(MMG5_xTetra));
      }
      else {
        pt->xt = 0;
        pt1->xt = 0;
      }
    }

    /* Update of adjacency relations */
    adja = &mesh->adja[4*(iel-1)+1];
    adjan = &mesh->adja[4*(jel-1)+1];

    nei2 = adja[tau[2]];
    nei3 = adja[tau[3]];

    /* Adjacency relations through both splitted faces */
    if ( k == 0 ) {
      if ( (list[1] / 6) == (nei2 / 4) ) {
        if ( MG_SMSGN(newtet[0],newtet[1]) ) {  //new elt of list[0] goes with new elt of list[1]
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*abs(newtet[1])+(nei2 %4);
        }
        else {
          adja[tau[2]] = 4*abs(newtet[1])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }

        if ( open ) {
          adja[tau[3]] = 0;
          adjan[tau[3]] = 0;
        }

        else {
          assert((list[ilist-1] / 6) == (nei3 / 4));
          if ( MG_SMSGN(newtet[0],newtet[ilist-1]) ) {
            adja[tau[3]] = nei3;
            adjan[tau[3]] = 4*abs(newtet[ilist-1])+(nei3 %4);
          }
          else {
            adja[tau[3]] = 4*abs(newtet[ilist-1])+(nei3 %4);
            adjan[tau[3]] = nei3;
          }
        }
      }

      else {
        assert((list[1] / 6) == (nei3 / 4));
        if ( MG_SMSGN(newtet[0],newtet[1]) ) {
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*abs(newtet[1])+(nei3 %4);
        }
        else {
          adja[tau[3]] = 4*abs(newtet[1])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }

        if ( open ) {
          adja[tau[2]] = 0;
          adjan[tau[2]] = 0;
        }

        else {
          assert((list[ilist-1]) / 6 == (nei2 / 4));
          if ( MG_SMSGN(newtet[0],newtet[ilist-1]) ) {
            adja[tau[2]] = nei2;
            adjan[tau[2]] = 4*abs(newtet[ilist-1])+(nei2 %4);
          }
          else {
            adja[tau[2]] = 4*abs(newtet[ilist-1])+(nei2 %4);
            adjan[tau[2]] = nei2;
          }
        }
      }
    }

    else if ( k==ilist-1 ) {
      if ( (list[ilist-2] / 6) == (nei2 / 4) ) {
        if ( MG_SMSGN(newtet[ilist-1],newtet[ilist-2]) ) {
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*abs(newtet[ilist-2])+(nei2 %4);
        }
        else {
          adja[tau[2]] = 4*abs(newtet[ilist-2])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }

        if ( open ) {
          adja[tau[3]] = 0;
          adjan[tau[3]] = 0;
        }

        else {
          assert((list[0]) / 6 == (nei3 / 4));
          if ( MG_SMSGN(newtet[ilist-1],newtet[0]) ) {
            adja[tau[3]] = nei3;
            adjan[tau[3]] = 4*abs(newtet[0])+(nei3 %4);
          }
          else {
            adja[tau[3]] = 4*abs(newtet[0])+(nei3 %4);
            adjan[tau[3]] = nei3;
          }
        }
      }

      else {
        assert((list[ilist-2] / 6) == (nei3 / 4));
        if ( MG_SMSGN(newtet[ilist-1],newtet[ilist-2]) ) {
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*abs(newtet[ilist-2])+(nei3 %4);
        }
        else {
          adja[tau[3]] = 4*abs(newtet[ilist-2])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }

        if ( open ) {
          adja[tau[2]] = 0;
          adjan[tau[2]] = 0;
        }

        else {
          assert((list[0]) / 6 == (nei2 / 4));
          if ( MG_SMSGN(newtet[ilist-1],newtet[0]) ) {
            adja[tau[2]] = nei2;
            adjan[tau[2]] = 4*abs(newtet[0])+(nei2 %4);
          }
          else {
            adja[tau[2]] = 4*abs(newtet[0])+(nei2 %4);
            adjan[tau[2]] = nei2;
          }
        }
      }
    }

    else {
      if ( (list[k-1] / 6) == (nei2 / 4) ) {
        if ( MG_SMSGN(newtet[k],newtet[k-1]) ) {
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*abs(newtet[k-1])+(nei2 %4);
        }
        else {
          adja[tau[2]] = 4*abs(newtet[k-1])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }

        assert((list[k+1]) / 6 == (nei3 / 4));
        if ( MG_SMSGN(newtet[k],newtet[k+1]) ) {
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*abs(newtet[k+1])+(nei3 %4);
        }
        else {
          adja[tau[3]] = 4*abs(newtet[k+1])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }
      }

      else {
        assert((list[k-1] / 6) == (nei3 / 4));
        if ( MG_SMSGN(newtet[k],newtet[k-1]) ) {
          adja[tau[3]] = nei3;
          adjan[tau[3]] = 4*abs(newtet[k-1])+(nei3 %4);
        }
        else {
          adja[tau[3]] = 4*abs(newtet[k-1])+(nei3 %4);
          adjan[tau[3]] = nei3;
        }

        assert((list[k+1]) / 6 == (nei2 / 4));
        if ( MG_SMSGN(newtet[k],newtet[k+1]) ) {
          adja[tau[2]] = nei2;
          adjan[tau[2]] = 4*abs(newtet[k+1])+(nei2 %4);
        }
        else {
          adja[tau[2]] = 4*abs(newtet[k+1])+(nei2 %4);
          adjan[tau[2]] = nei2;
        }
      }
    }

    /* Internal adjacency relations update */
    mel = adja[tau[0]] / 4;
    voy = adja[tau[0]] % 4;
    adja[tau[0]] = 4*jel + tau[1];
    adjan[tau[0]] = 4*mel + voy;
    adjan[tau[1]] = 4*iel + tau[0];

    if ( mel ) {
      adjan = &mesh->adja[4*(mel -1) +1];
      adjan[voy] = 4*jel + tau[0];
    }
    /* Quality update */
    if ( (!metRidTyp) && met->m && met->size>1 ) {
      pt->qual=MMG5_caltet33_ani(mesh,met,pt);
      pt1->qual=MMG5_caltet33_ani(mesh,met,pt1);
    }
    else {
      pt->qual=MMG5_orcal(mesh,met,iel);
      pt1->qual=MMG5_orcal(mesh,met,jel);
    }
    pt->mark  = mesh->mark;
    pt1->mark = mesh->mark;
  }

  MMG5_SAFE_FREE(newtet);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 * Simulate split of two edges that belong to a common face
 *
 */
int MMG3D_split2sf_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]){
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4],imin;
  const unsigned char *taued;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* identity is case 48 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(pt->flag){
  case 24 :
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    taued = &MMG5_permedge[1][0];
    break;
  case 40 :
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  case 6 :
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &MMG5_permedge[5][0];
    break;
  case 34 :
    tau[0] = 1 ; tau[1] = 0 ; tau[2] = 3 ; tau[3] = 2;
    taued = &MMG5_permedge[3][0];
    break;
  case 36 :
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;
  case 20 :
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &MMG5_permedge[6][0];
    break;
  case 5 :
    tau[0] = 2 ; tau[1] = 1 ; tau[2] = 3 ; tau[3] = 0;
    taued = &MMG5_permedge[7][0];
    break;
  case 17 :
    tau[0] = 2 ; tau[1] = 3 ; tau[2] = 0 ; tau[3] = 1;
    taued = &MMG5_permedge[8][0];
    break;
  case 9 :
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &MMG5_permedge[9][0];
    break;
  case 3 :
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &MMG5_permedge[11][0];
    break;
  case 10 :
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  }

  /* Test orientation of the three tets to be created */
  imin = (pt->v[tau[1]] < pt->v[tau[2]]) ? tau[1] : tau[2] ;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[1]] = vx[taued[4]];
  pt0->v[tau[2]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  if ( imin == tau[1] ) {
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[2]] = vx[taued[5]];
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[3]] = vx[taued[5]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else {
    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[1]] = vx[taued[4]];
    pt0->v[tau[3]] = vx[taued[5]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split of two edges that belong to a common face : 1 tetra becomes 3
 *
 */
int MMG5_split2sf(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp){
  MMG5_pTetra         pt[3];
  MMG5_xTetra         xt[3];
  MMG5_pxTetra        pxt0;
  int                 iel,i;
  int                 newtet[3];
  char                flg,imin,firstxt,isxt[3];
  unsigned char       tau[4];
  const unsigned char *taued;

  pt[0] = &mesh->tetra[k];
  flg   = pt[0]->flag;
  pt[0]->flag = 0;
  newtet[0]=k;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
  }
  pt[1] = &mesh->tetra[iel];
  memcpy(pt[1],pt[0],sizeof(MMG5_Tetra));
  newtet[1]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
  }
  pt[2] = &mesh->tetra[iel];
  memcpy(pt[2],pt[0],sizeof(MMG5_Tetra));
  newtet[2]=iel;

  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[1],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[2],pxt0,sizeof(MMG5_xTetra));
  }
  else {
    pxt0 = 0;
    memset(&xt[0],0,sizeof(MMG5_xTetra));
    memset(&xt[1],0,sizeof(MMG5_xTetra));
    memset(&xt[2],0,sizeof(MMG5_xTetra));
  }
  /* identity is case 48 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(flg){
  case 24 :
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    taued = &MMG5_permedge[1][0];
    break;
  case 40 :
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  case 6 :
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &MMG5_permedge[5][0];
    break;
  case 34 :
    tau[0] = 1 ; tau[1] = 0 ; tau[2] = 3 ; tau[3] = 2;
    taued = &MMG5_permedge[3][0];
    break;
  case 36 :
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;
  case 20 :
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &MMG5_permedge[6][0];
    break;
  case 5 :
    tau[0] = 2 ; tau[1] = 1 ; tau[2] = 3 ; tau[3] = 0;
    taued = &MMG5_permedge[7][0];
    break;
  case 17 :
    tau[0] = 2 ; tau[1] = 3 ; tau[2] = 0 ; tau[3] = 1;
    taued = &MMG5_permedge[8][0];
    break;
  case 9 :
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &MMG5_permedge[9][0];
    break;
  case 3 :
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    taued = &MMG5_permedge[11][0];
    break;
  case 10 :
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  }

  /* Generic formulation for the split of 2 edges belonging to a common face */
  imin = (pt[0]->v[tau[1]] < pt[0]->v[tau[2]]) ? tau[1] : tau[2] ;
  pt[0]->v[tau[1]]  = vx[taued[4]] ;  pt[0]->v[tau[2]] = vx[taued[5]];
  xt[0].tag[taued[0]] = 0;  xt[0].tag[taued[1]] = 0;
  xt[0].tag[taued[3]] = 0;  xt[0].edg[taued[0]] = 0;
  xt[0].edg[taued[1]] = 0;  xt[0].edg[taued[3]] = 0;
  xt[0].ref[  tau[3]] = 0;  xt[0].ftag[ tau[3]] = 0;  MG_SET(xt[0].ori, tau[3]);

  if ( imin == tau[1] ) {
    pt[1]->v[tau[2]] = vx[taued[5]];  pt[1]->v[tau[3]] = vx[taued[4]];
    pt[2]->v[tau[3]] = vx[taued[5]];

    xt[1].tag[taued[1]] = 0;  xt[1].tag[taued[2]] = 0;
    xt[1].tag[taued[3]] = 0;  xt[1].tag[taued[5]] = 0;
    xt[1].edg[taued[1]] = 0;  xt[1].edg[taued[2]] = 0;
    xt[1].edg[taued[3]] = 0;  xt[1].edg[taued[5]] = 0;
    xt[1].ref [ tau[1]] = 0;  xt[1].ref [ tau[3]] = 0;
    xt[1].ftag[ tau[1]] = 0;  xt[1].ftag[ tau[3]] = 0;
    MG_SET(xt[1].ori, tau[1]);  MG_SET(xt[1].ori, tau[3]);

    xt[2].tag[taued[2]] = 0;  xt[2].tag[taued[4]] = 0;
    xt[2].edg[taued[2]] = 0;  xt[2].edg[taued[4]] = 0;
    xt[2].ref[  tau[2]] = 0;  xt[2].ftag[ tau[2]] = 0;  MG_SET(xt[2].ori, tau[2]);
  }
  else {
    pt[1]->v[tau[3]] = vx[taued[4]];
    pt[2]->v[tau[1]] = vx[taued[4]];  pt[2]->v[tau[3]] = vx[taued[5]];

    xt[1].tag[taued[2]] = 0;  xt[1].tag[taued[5]] = 0;
    xt[1].edg[taued[2]] = 0;  xt[1].edg[taued[5]] = 0;
    xt[1].ref[  tau[1]] = 0;  xt[1].ftag[ tau[1]] = 0;  MG_SET(xt[1].ori, tau[1]);

    xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[2]] = 0;
    xt[2].tag[taued[3]] = 0;  xt[2].tag[taued[4]] = 0;
    xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[2]] = 0;
    xt[2].edg[taued[3]] = 0;  xt[2].edg[taued[4]] = 0;
    xt[2].ref [ tau[2]] = 0;  xt[2].ref [ tau[3]] = 0;
    xt[2].ftag[ tau[2]] = 0;  xt[2].ftag[ tau[3]] = 0;
    MG_SET(xt[2].ori, tau[2]);  MG_SET(xt[2].ori, tau[3]);
  }

  /* Assignation of the xt fields to the appropriate tets */
  isxt[0] = isxt[1] = isxt[2] = 0;
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      pt[1]->xt = pt[2]->xt = 0;
      for (i=1; i<3; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&(xt[i]),sizeof(MMG5_xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = 0;
      for (i=1; i<3; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[pt[i]->xt];
            memcpy(pxt0,&(xt[i]),sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    pt[0]->qual=MMG5_caltet33_ani(mesh,met,pt[0]);
    pt[1]->qual=MMG5_caltet33_ani(mesh,met,pt[1]);
    pt[2]->qual=MMG5_caltet33_ani(mesh,met,pt[2]);
  }
  else
  {
    pt[0]->qual=MMG5_orcal(mesh,met,newtet[0]);
    pt[1]->qual=MMG5_orcal(mesh,met,newtet[1]);
    pt[2]->qual=MMG5_orcal(mesh,met,newtet[2]);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 *  Simulate split of two opposite edges.
 *
 */
int MMG3D_split2_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]){
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4];
  const unsigned char *taued;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* identity is case 33 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(pt->flag){
  case 18:
    tau[0] = 3;  tau[1] = 1;  tau[2] = 0;  tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  case 12:
    tau[0] = 0;  tau[1] = 3;  tau[2] = 1;  tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  }

  /* Test orientation of the 4 tets to be created */
  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[1]] = vx[taued[0]];  pt0->v[tau[2]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[1]] = vx[taued[0]];  pt0->v[tau[3]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[0]];  pt0->v[tau[2]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[0]];  pt0->v[tau[3]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split of two OPPOSITE edges
 *
 */
int MMG5_split2(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra         pt[4];
  MMG5_xTetra         xt[4];
  MMG5_pxTetra        pxt0;
  int                 i,iel;
  int                 newtet[4];
  char                flg,firstxt,isxt[4];
  unsigned char       tau[4];
  const unsigned char *taued;

  pt[0] = &mesh->tetra[k];
  flg   = pt[0]->flag;
  pt[0]->flag = 0;
  newtet[0]=k;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
  }
  pt[1] = &mesh->tetra[iel];
  memcpy(pt[1],pt[0],sizeof(MMG5_Tetra));
  newtet[1]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
  }
  pt[2] = &mesh->tetra[iel];
  memcpy(pt[2],pt[0],sizeof(MMG5_Tetra));
  newtet[2]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
    pt[2] = &mesh->tetra[newtet[2]];
  }
  pt[3] = &mesh->tetra[iel];
  memcpy(pt[3],pt[0],sizeof(MMG5_Tetra));
  newtet[3]=iel;

  pxt0 = 0;
  if ( pt[0]->xt) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[1],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[2],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[3],pxt0,sizeof(MMG5_xTetra));
  }
  else {
    memset(&xt[0],0,sizeof(MMG5_xTetra));
    memset(&xt[1],0,sizeof(MMG5_xTetra));
    memset(&xt[2],0,sizeof(MMG5_xTetra));
    memset(&xt[3],0,sizeof(MMG5_xTetra));
  }
  /* identity : case 33 */
  tau[0] = 0;  tau[1] = 1;  tau[2] = 2;  tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(flg){
  case 18:
    tau[0] = 3;  tau[1] = 1;  tau[2] = 0;  tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  case 12:
    tau[0] = 0;  tau[1] = 3;  tau[2] = 1;  tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  }

  /* Generic formulation for the split of 2 opposite edges */
  pt[0]->v[tau[1]] = vx[taued[0]];  pt[0]->v[tau[2]] = vx[taued[5]];
  pt[1]->v[tau[1]] = vx[taued[0]];  pt[1]->v[tau[3]] = vx[taued[5]];
  pt[2]->v[tau[0]] = vx[taued[0]];  pt[2]->v[tau[2]] = vx[taued[5]];
  pt[3]->v[tau[0]] = vx[taued[0]];  pt[3]->v[tau[3]] = vx[taued[5]];

  xt[0].tag[taued[1]] = 0;  xt[0].tag[taued[3]] = 0;
  xt[0].tag[taued[4]] = 0;  xt[0].edg[taued[1]] = 0;
  xt[0].edg[taued[3]] = 0;  xt[0].edg[taued[4]] = 0;
  xt[0].ref [ tau[0]] = 0;  xt[0].ref [ tau[3]] = 0;
  xt[0].ftag[ tau[0]] = 0;  xt[0].ftag[ tau[3]] = 0;
  MG_SET(xt[0].ori, tau[0]);  MG_SET(xt[0].ori, tau[3]);

  xt[1].tag[taued[2]] = 0;  xt[1].tag[taued[3]] = 0;
  xt[1].tag[taued[4]] = 0;  xt[1].edg[taued[2]] = 0;
  xt[1].edg[taued[3]] = 0;  xt[1].edg[taued[4]] = 0;
  xt[1].ref [ tau[0]] = 0;  xt[1].ref [ tau[2]] = 0;
  xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[2]] = 0;
  MG_SET(xt[1].ori, tau[0]);  MG_SET(xt[1].ori, tau[2]);

  xt[2].tag[taued[1]] = 0;  xt[2].tag[taued[2]] = 0;
  xt[2].tag[taued[3]] = 0;  xt[2].edg[taued[1]] = 0;
  xt[2].edg[taued[2]] = 0;  xt[2].edg[taued[3]] = 0;
  xt[2].ref [ tau[1]] = 0;  xt[2].ref [ tau[3]] = 0;
  xt[2].ftag[ tau[1]] = 0;  xt[2].ftag[ tau[3]] = 0;
  MG_SET(xt[2].ori, tau[1]);  MG_SET(xt[2].ori, tau[3]);

  xt[3].tag[taued[1]] = 0;  xt[3].tag[taued[2]] = 0;
  xt[3].tag[taued[4]] = 0;  xt[3].edg[taued[1]] = 0;
  xt[3].edg[taued[2]] = 0;  xt[3].edg[taued[4]] = 0;
  xt[3].ref [ tau[1]] = 0;  xt[3].ref [ tau[2]] = 0;
  xt[3].ftag[ tau[1]] = 0;  xt[3].ftag[ tau[2]] = 0;
  MG_SET(xt[3].ori, tau[1]);  MG_SET(xt[3].ori, tau[2]);

  /* Assignation of the xt fields to the appropriate tets */
  memset(isxt,0,4*sizeof(char));
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] )  isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] )  isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] )  isxt[2] = 1;
    if ( xt[3].ref[i] || xt[3].ftag[i] )  isxt[3] = 1;
  }

  if ( pt[0]->xt) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
        }
        else {
          pt[i]->xt = 0;
        }
      }
    }
    else {
      firstxt = 1;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[pt[i]->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
        else {
          pt[i]->xt = 0;
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    pt[0]->qual=MMG5_caltet33_ani(mesh,met,pt[0]);
    pt[1]->qual=MMG5_caltet33_ani(mesh,met,pt[1]);
    pt[2]->qual=MMG5_caltet33_ani(mesh,met,pt[2]);
    pt[3]->qual=MMG5_caltet33_ani(mesh,met,pt[3]);
  }
  else {
    pt[0]->qual=MMG5_orcal(mesh,met,newtet[0]);
    pt[1]->qual=MMG5_orcal(mesh,met,newtet[1]);
    pt[2]->qual=MMG5_orcal(mesh,met,newtet[2]);
    pt[3]->qual=MMG5_orcal(mesh,met,newtet[3]);
  }
  return 1;
}

/** Simulate split of 1 face (3 edges) */
int MMG3D_split3_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4];
  const unsigned char *taued;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* identity is case 11 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(pt->flag) {
  case 21:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  case 38:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &MMG5_permedge[9][0];
    break;
  case 56:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &MMG5_permedge[5][0];
    break;
  }

  /* Check orientation of the 4 newly created tets */
  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[1]] = vx[taued[0]];
  pt0->v[tau[2]] = vx[taued[1]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[0]];
  pt0->v[tau[2]] = vx[taued[3]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[1]];
  pt0->v[tau[1]] = vx[taued[3]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[0]];
  pt0->v[tau[1]] = vx[taued[3]];
  pt0->v[tau[2]] = vx[taued[1]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * 1 face (3 edges) subdivided
 *
 */
int MMG5_split3(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra         pt[4];
  MMG5_xTetra         xt[4];
  MMG5_pxTetra        pxt0;
  int                 iel,i;
  int                 newtet[4];
  char                flg,firstxt,isxt[4];
  unsigned char       tau[4];
  const unsigned char *taued;

  pt[0] = &mesh->tetra[k];
  flg   = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

  /* create 3 new tetras */
  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
  }
  pt[1] = &mesh->tetra[iel];
  pt[1] = memcpy(pt[1],pt[0],sizeof(MMG5_Tetra));
  newtet[1]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate a"
                                " new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
  }
  pt[2] = &mesh->tetra[iel];
  pt[2] = memcpy(pt[2],pt[0],sizeof(MMG5_Tetra));
  newtet[2]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
    pt[2] = &mesh->tetra[newtet[2]];
  }
  pt[3] = &mesh->tetra[iel];
  pt[3] = memcpy(pt[3],pt[0],sizeof(MMG5_Tetra));
  newtet[3]=iel;

  pxt0 = 0;
  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[1],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[2],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[3],pxt0, sizeof(MMG5_xTetra));
  }
  else {
    memset(&xt[0],0, sizeof(MMG5_xTetra));
    memset(&xt[1],0, sizeof(MMG5_xTetra));
    memset(&xt[2],0, sizeof(MMG5_xTetra));
    memset(&xt[3],0, sizeof(MMG5_xTetra));
  }

  /* update vertices, case 11 is default */
  tau[0] = 0; tau[1] = 1; tau[2] = 2; tau[3] = 3;
  taued = &MMG5_permedge[0][0];
  switch(flg) {
  case 21:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    taued = &MMG5_permedge[2][0];
    break;
  case 38:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    taued = &MMG5_permedge[9][0];
    break;
  case 56:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &MMG5_permedge[5][0];
    break;
  }

  /* Generic formulation of split of 3 edges */
  pt[0]->v[tau[1]] = vx[taued[0]];  pt[0]->v[tau[2]] = vx[taued[1]];
  pt[1]->v[tau[0]] = vx[taued[0]];  pt[1]->v[tau[2]] = vx[taued[3]];
  pt[2]->v[tau[0]] = vx[taued[1]];  pt[2]->v[tau[1]] = vx[taued[3]];
  pt[3]->v[tau[0]] = vx[taued[0]];  pt[3]->v[tau[1]] = vx[taued[3]];  pt[3]->v[tau[2]] = vx[taued[1]];

  xt[0].tag[taued[3]] = 0;  xt[0].tag[taued[4]] = 0;
  xt[0].tag[taued[5]] = 0;  xt[0].edg[taued[3]] = 0;
  xt[0].edg[taued[4]] = 0;  xt[0].edg[taued[5]] = 0;
  xt[0].ref[  tau[0]] = 0;  xt[0].ftag[ tau[0]] = 0;  MG_SET(xt[0].ori, tau[0]);

  xt[1].tag[taued[1]] = 0;  xt[1].tag[taued[2]] = 0;
  xt[1].tag[taued[5]] = 0;  xt[1].edg[taued[1]] = 0;
  xt[1].edg[taued[2]] = 0;  xt[1].edg[taued[5]] = 0;
  xt[1].ref[  tau[1]] = 0;  xt[1].ftag[ tau[1]] = 0;  MG_SET(xt[1].ori, tau[1]);

  xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[2]] = 0;
  xt[2].tag[taued[4]] = 0;  xt[2].edg[taued[0]] = 0;
  xt[2].edg[taued[2]] = 0;  xt[2].edg[taued[4]] = 0;
  xt[2].ref[  tau[2]] = 0;  xt[2].ftag[ tau[2]] = 0;  MG_SET(xt[2].ori, tau[2]);

  xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[1]] = 0;
  xt[3].tag[taued[2]] = 0;  xt[3].tag[taued[3]] = 0;
  xt[3].tag[taued[4]] = 0;  xt[3].tag[taued[5]] = 0;
  xt[3].edg[taued[0]] = 0;  xt[3].edg[taued[1]] = 0;
  xt[3].edg[taued[2]] = 0;  xt[3].edg[taued[3]] = 0;
  xt[3].edg[taued[4]] = 0;  xt[3].edg[taued[5]] = 0;
  xt[3].ref [ tau[0]] = 0;  xt[3].ref [ tau[1]] = 0;  xt[3].ref [tau[2]] = 0;
  xt[3].ftag[ tau[0]] = 0;  xt[3].ftag[ tau[1]] = 0;  xt[3].ftag[tau[2]] = 0;
  MG_SET(xt[3].ori, tau[0]);  MG_SET(xt[3].ori, tau[1]);  MG_SET(xt[3].ori, tau[2]);

  /* Assignation of the xt fields to the appropriate tets */
  memset(isxt,0,4*sizeof(char));
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    if ( xt[3].ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&(xt[i]),sizeof(MMG5_xTetra));
          }
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    pt[0]->qual=MMG5_caltet33_ani(mesh,met,pt[0]);
    pt[1]->qual=MMG5_caltet33_ani(mesh,met,pt[1]);
    pt[2]->qual=MMG5_caltet33_ani(mesh,met,pt[2]);
    pt[3]->qual=MMG5_caltet33_ani(mesh,met,pt[3]);
  }
  else {
    pt[0]->qual=MMG5_orcal(mesh,met,newtet[0]);
    pt[1]->qual=MMG5_orcal(mesh,met,newtet[1]);
    pt[2]->qual=MMG5_orcal(mesh,met,newtet[2]);
    pt[3]->qual=MMG5_orcal(mesh,met,newtet[3]);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 *  Simulate split of 3 edges in cone configuration.
 *
 */
int MMG3D_split3cone_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4],ia,ib;
  const unsigned char *taued;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* identity is case 7 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];

  switch(pt->flag) {
  case 25:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;

  case 42:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &MMG5_permedge[6][0];
    break;

  case 52:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  }

  /* Generic formulation of split of 3 edges in cone configuration (edges 0,1,2 splitted) */
  /* Fill ia,ib,ic so that pt->v[ia] < pt->v[ib] < pt->v[ic] */
  if ( pt->v[tau[1]] < pt->v[tau[2]] ) {
    ia = tau[1];
    ib = tau[2];
  }
  else {
    ia = tau[2];
    ib = tau[1];
  }

  if ( pt->v[tau[3]] < pt->v[ia] ) {
    ib = ia;
    ia = tau[3];
  }
  else {
    if ( pt->v[tau[3]] < pt->v[ib] ) {
      ib = tau[3];
    }
    else {
    }
  }

  /* Check orientation of the 4 newly created tets */
  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[1]] = vx[taued[0]];
  pt0->v[tau[2]] = vx[taued[1]];
  pt0->v[tau[3]] = vx[taued[2]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  if ( ia == tau[3] ) {
    pt0->v[tau[0]] = vx[taued[2]];
    pt0->v[tau[1]] = vx[taued[0]];
    pt0->v[tau[2]] = vx[taued[1]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    if ( ib == tau[1] ) {
      pt0->v[tau[0]] = vx[taued[0]];
      pt0->v[tau[2]] = vx[taued[1]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;

      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[tau[0]] = vx[taued[1]] ;
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;
    }
    else {
      assert(ib == tau[2]);
      pt0->v[tau[0]] = vx[taued[1]];
      pt0->v[tau[1]] = vx[taued[0]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;

      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[tau[0]] = vx[taued[0]] ;
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;
    }
  }
  else if (ia == tau[2] ) {
    pt0->v[tau[0]] = vx[taued[1]];
    pt0->v[tau[1]] = vx[taued[0]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    if ( ib == tau[3] ) {
      pt0->v[tau[0]] = vx[taued[2]];
      pt0->v[tau[1]] = vx[taued[0]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;

      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[tau[0]] = vx[taued[0]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;
    }
    else {
      assert(ib == tau[1]);
      pt0->v[tau[0]] = vx[taued[0]];
      pt0->v[tau[3]] = vx[taued[2]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;

      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[tau[0]] = vx[taued[2]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;
    }
  }
  else {
    assert(ia == tau[1]);

    pt0->v[tau[0]] = vx[taued[0]];
    pt0->v[tau[2]] = vx[taued[1]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    if ( ib == tau[2] ) {
      pt0->v[tau[0]] = vx[taued[1]];
      pt0->v[tau[3]] = vx[taued[2]] ;
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;

      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[tau[0]] = vx[taued[2]] ;
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;

    }
    else {
      assert(ib == tau[3]);

      pt0->v[tau[0]] = vx[taued[2]];
      pt0->v[tau[2]] = vx[taued[1]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;

      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[tau[0]] = vx[taued[1]];
      vnew = MMG5_orvol(mesh->point,pt0->v);
      if ( vnew < MMG5_EPSOK )  return 0;
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split 3 edge in cone configuration
 *
 */
int MMG5_split3cone(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra         pt[4];
  MMG5_xTetra         xt[4];
  MMG5_pxTetra        pxt0;
  int                 iel,i;
  int                 newtet[4];
  char                flg,firstxt,isxt[4],ia,ib;
  unsigned char       tau[4];
  const unsigned char *taued;

  pt[0]  = &mesh->tetra[k];
  flg = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

  /* create 3 new tetras */
  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
  }
  pt[1] = &mesh->tetra[iel];
  memcpy(pt[1],pt[0],sizeof(MMG5_Tetra));
  newtet[1]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
  }
  pt[2] = &mesh->tetra[iel];
  memcpy(pt[2],pt[0],sizeof(MMG5_Tetra));
  newtet[2]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate a"
                                " new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
    pt[2] = &mesh->tetra[newtet[2]];
  }
  pt[3] = &mesh->tetra[iel];
  memcpy(pt[3],pt[0],sizeof(MMG5_Tetra));
  newtet[3]=iel;

  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[1],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[2],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[3],pxt0, sizeof(MMG5_xTetra));
  }
  else {
    pxt0 = 0;
    memset(&xt[0],0, sizeof(MMG5_xTetra));
    memset(&xt[1],0, sizeof(MMG5_xTetra));
    memset(&xt[2],0, sizeof(MMG5_xTetra));
    memset(&xt[3],0, sizeof(MMG5_xTetra));
  }

  /* Set permutation of vertices : reference configuration is 7 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];

  switch(flg) {
  case 25:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;

  case 42:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    taued = &MMG5_permedge[6][0];
    break;

  case 52:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    taued = &MMG5_permedge[10][0];
    break;
  }

  /* Generic formulation of split of 3 edges in cone configuration (edges 0,1,2 splitted) */
  /* Fill ia,ib,ic so that pt->v[ia] < pt->v[ib] < pt->v[ic] */
  if ( (pt[0])->v[tau[1]] < (pt[0])->v[tau[2]] ) {
    ia = tau[1];
    ib = tau[2];
  }
  else {
    ia = tau[2];
    ib = tau[1];
  }

  if ( (pt[0])->v[tau[3]] < (pt[0])->v[ia] ) {
    ib = ia;
    ia = tau[3];
  }
  else {
    if ( (pt[0])->v[tau[3]] < (pt[0])->v[ib] ) {
      ib = tau[3];
    }
    else {
    }
  }

  pt[0]->v[tau[1]] = vx[taued[0]] ; pt[0]->v[tau[2]] = vx[taued[1]] ; pt[0]->v[tau[3]] = vx[taued[2]];
  xt[0].tag[taued[3]] = 0;  xt[0].tag[taued[4]] = 0;
  xt[0].tag[taued[5]] = 0;  xt[0].edg[taued[3]] = 0;
  xt[0].edg[taued[4]] = 0;  xt[0].edg[taued[5]] = 0;
  xt[0].ref [ tau[0]] = 0;
  xt[0].ftag[ tau[0]] = 0;
  MG_SET(xt[0].ori, tau[0]);

  if ( ia == tau[3] ) {
    pt[1]->v[tau[0]] = vx[taued[2]] ; pt[1]->v[tau[1]] = vx[taued[0]] ; pt[1]->v[tau[2]] = vx[taued[1]];
    xt[1].tag[taued[0]] = 0;  xt[1].tag[taued[1]] = 0;
    xt[1].tag[taued[3]] = 0;  xt[1].tag[taued[4]] = 0;
    xt[1].tag[taued[5]] = 0;  xt[1].edg[taued[0]] = 0;
    xt[1].edg[taued[1]] = 0;  xt[1].edg[taued[3]] = 0;
    xt[1].edg[taued[4]] = 0;  xt[1].edg[taued[5]] = 0;
    xt[1].ref [ tau[0]] = 0;  xt[1].ref [ tau[3]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[3]] = 0;
    MG_SET(xt[1].ori, tau[0]);  MG_SET(xt[1].ori, tau[3]);

    if ( ib == tau[1] ) {
      pt[2]->v[tau[0]] = vx[taued[0]] ; pt[2]->v[tau[2]] = vx[taued[1]] ;
      xt[2].tag[taued[1]] = 0;  xt[2].tag[taued[2]] = 0;
      xt[2].tag[taued[3]] = 0;  xt[2].tag[taued[5]] = 0;
      xt[2].edg[taued[1]] = 0;  xt[2].edg[taued[2]] = 0;
      xt[2].edg[taued[3]] = 0;  xt[2].edg[taued[5]] = 0;
      xt[2].ref [ tau[0]] = 0;  xt[2].ref [ tau[1]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[1]] = 0;
      MG_SET(xt[2].ori, tau[0]);  MG_SET(xt[2].ori, tau[1]);

      pt[3]->v[tau[0]] = vx[taued[1]] ;
      xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[2]] = 0;
      xt[3].edg[taued[0]] = 0;  xt[3].edg[taued[2]] = 0;
      xt[3].ref [ tau[2]] = 0;
      xt[3].ftag[ tau[2]] = 0;
      MG_SET(xt[3].ori, tau[2]);
    }
    else {
      assert(ib == tau[2]);

      pt[2]->v[tau[0]] = vx[taued[1]] ; pt[2]->v[tau[1]] = vx[taued[0]] ;
      xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[2]] = 0;
      xt[2].tag[taued[3]] = 0;  xt[2].tag[taued[4]] = 0;
      xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[2]] = 0;
      xt[2].edg[taued[3]] = 0;  xt[2].edg[taued[4]] = 0;
      xt[2].ref [ tau[0]] = 0;  xt[2].ref [ tau[2]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[2]] = 0;
      MG_SET(xt[2].ori, tau[0]);  MG_SET(xt[2].ori, tau[2]);

      pt[3]->v[tau[0]] = vx[taued[0]] ;
      xt[3].tag[taued[1]] = 0;  xt[3].tag[taued[2]] = 0;
      xt[3].edg[taued[1]] = 0;  xt[3].edg[taued[2]] = 0;
      xt[3].ref [ tau[1]] = 0;
      xt[3].ftag[ tau[1]] = 0;
      MG_SET(xt[3].ori, tau[1]);
    }
  }

  else if (ia == tau[2] ) {
    pt[1]->v[tau[0]] = vx[taued[1]] ; pt[1]->v[tau[1]] = vx[taued[0]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].tag[taued[0]] = 0;  xt[1].tag[taued[2]] = 0;
    xt[1].tag[taued[3]] = 0;  xt[1].tag[taued[4]] = 0;
    xt[1].tag[taued[5]] = 0;  xt[1].edg[taued[0]] = 0;
    xt[1].edg[taued[2]] = 0;  xt[1].edg[taued[3]] = 0;
    xt[1].edg[taued[4]] = 0;  xt[1].edg[taued[5]] = 0;
    xt[1].ref [ tau[0]] = 0;  xt[1].ref [ tau[2]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[2]] = 0;
    MG_SET(xt[1].ori, tau[0]);  MG_SET(xt[1].ori, tau[2]);

    if ( ib == tau[3] ) {
      pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[1]] = vx[taued[0]] ;
      xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[1]] = 0;
      xt[2].tag[taued[3]] = 0;  xt[2].tag[taued[4]] = 0;
      xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[1]] = 0;
      xt[2].edg[taued[3]] = 0;  xt[2].edg[taued[4]] = 0;
      xt[2].ref [ tau[0]] = 0;  xt[2].ref [ tau[3]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[3]] = 0;
      MG_SET(xt[2].ori, tau[0]);  MG_SET(xt[2].ori, tau[3]);

      pt[3]->v[tau[0]] = vx[taued[0]] ;
      xt[3].tag[taued[1]] = 0;  xt[3].tag[taued[2]] = 0;
      xt[3].edg[taued[1]] = 0;  xt[3].edg[taued[2]] = 0;
      xt[3].ref [ tau[1]] = 0;
      xt[3].ftag[ tau[1]] = 0;
      MG_SET(xt[3].ori, tau[1]);
    }
    else {
      assert(ib == tau[1]);

      pt[2]->v[tau[0]] = vx[taued[0]] ; pt[2]->v[tau[3]] = vx[taued[2]] ;
      xt[2].tag[taued[1]] = 0;  xt[2].tag[taued[2]] = 0;
      xt[2].tag[taued[4]] = 0;  xt[2].tag[taued[5]] = 0;
      xt[2].edg[taued[1]] = 0;  xt[2].edg[taued[2]] = 0;
      xt[2].edg[taued[4]] = 0;  xt[2].edg[taued[5]] = 0;
      xt[2].ref [ tau[0]] = 0;  xt[2].ref [ tau[1]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[1]] = 0;
      MG_SET(xt[2].ori, tau[0]);  MG_SET(xt[2].ori, tau[1]);

      pt[3]->v[tau[0]] = vx[taued[2]] ;
      xt[3].tag[taued[0]] = 0;    xt[3].tag[taued[1]] = 0;
      xt[3].edg[taued[0]] = 0;    xt[3].edg[taued[1]] = 0;
      xt[3].ref [ tau[3]] = 0;
      xt[3].ftag[ tau[3]] = 0;
      MG_SET(xt[3].ori, tau[3]);
    }
  }
  else {
    assert(ia == tau[1]);

    pt[1]->v[tau[0]] = vx[taued[0]] ; pt[1]->v[tau[2]] = vx[taued[1]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].tag[taued[1]] = 0;  xt[1].tag[taued[2]] = 0;
    xt[1].tag[taued[3]] = 0;  xt[1].tag[taued[4]] = 0;
    xt[1].tag[taued[5]] = 0;  xt[1].edg[taued[1]] = 0;
    xt[1].edg[taued[2]] = 0;  xt[1].edg[taued[3]] = 0;
    xt[1].edg[taued[4]] = 0;  xt[1].edg[taued[5]] = 0;
    xt[1].ref [ tau[0]] = 0;  xt[1].ref [ tau[1]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[1]] = 0;
    MG_SET(xt[1].ori, tau[0]);  MG_SET(xt[1].ori, tau[1]);

    if ( ib == tau[2] ) {
      pt[2]->v[tau[0]] = vx[taued[1]] ; pt[2]->v[tau[3]] = vx[taued[2]] ;
      xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[2]] = 0;
      xt[2].tag[taued[4]] = 0;  xt[2].tag[taued[5]] = 0;
      xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[2]] = 0;
      xt[2].edg[taued[4]] = 0;  xt[2].edg[taued[5]] = 0;
      xt[2].ref [ tau[0]] = 0;  xt[2].ref [ tau[2]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[2]] = 0;
      MG_SET(xt[2].ori, tau[0]);  MG_SET(xt[2].ori, tau[2]);

      pt[3]->v[tau[0]] = vx[taued[2]] ;
      xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[1]] = 0;
      xt[3].edg[taued[0]] = 0;  xt[3].edg[taued[1]] = 0;
      xt[3].ref [ tau[3]] = 0;
      xt[3].ftag[ tau[3]] = 0;
      MG_SET(xt[3].ori, tau[3]);
    }
    else {
      assert(ib == tau[3]);

      pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[2]] = vx[taued[1]] ;
      xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[1]] = 0;
      xt[2].tag[taued[3]] = 0;  xt[2].tag[taued[5]] = 0;
      xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[1]] = 0;
      xt[2].edg[taued[3]] = 0;  xt[2].edg[taued[5]] = 0;
      xt[2].ref [ tau[0]] = 0;  xt[2].ref [ tau[3]] = 0;
      xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[3]] = 0;
      MG_SET(xt[2].ori, tau[0]);  MG_SET(xt[2].ori, tau[3]);

      pt[3]->v[tau[0]] = vx[taued[1]] ;
      xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[2]] = 0;
      xt[3].edg[taued[0]] = 0;  xt[3].edg[taued[2]] = 0;
      xt[3].ref [ tau[2]] = 0;
      xt[3].ftag[ tau[2]] = 0;
      MG_SET(xt[3].ori, tau[2]);
    }
  }

  /* Assignation of the xt fields to the appropriate tets */
  isxt[0] = isxt[1] = isxt[2] = isxt[3] = 0;

  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    if ( xt[3].ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
  }

  if ( (pt[0])->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;
      for ( i=1; i<4; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
      }
      (pt[0])->xt = 0;
    }
  }
  /* Quality update */
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    pt[0]->qual=MMG5_caltet33_ani(mesh,met,pt[0]);
    pt[1]->qual=MMG5_caltet33_ani(mesh,met,pt[1]);
    pt[2]->qual=MMG5_caltet33_ani(mesh,met,pt[2]);
    pt[3]->qual=MMG5_caltet33_ani(mesh,met,pt[3]);
  }
  else {
    pt[0]->qual=MMG5_orcal(mesh,met,newtet[0]);
    pt[1]->qual=MMG5_orcal(mesh,met,newtet[1]);
    pt[2]->qual=MMG5_orcal(mesh,met,newtet[2]);
    pt[3]->qual=MMG5_orcal(mesh,met,newtet[3]);
  }
  return 1;
}

/**
 * \param pt initial tetra
 * \param vx index of points to insert along edges
 * \param tau vertices permutation
 * \param taued edges permutation
 * \param sym vertices symmetry
 * \param symed edges symmetry
 * \param ip0 vertex 0 for reference config
 * \param ip1 vertex 1 for reference config
 * \param ip2 vertex 2 for reference config
 * \param ip3 vertex 3 for reference config
 * \param ie0 edge 0 for reference config
 * \param ie1 edge 1 for reference config
 * \param ie2 edge 2 for reference config
 * \param ie3 edge 3 for reference config
 * \param ie4 edge 4 for reference config
 * \param ie5 edge 5 for reference config
 * \param imin03 minimal index of vertices ip0 and ip3
 * \param imin12 minimal index of vertices ip1 and ip2
 *
 * Set permutation /symmetry of vertices for 3 opposite edges config:
 * generic case : 35
 *
 */
static inline
void MMG3D_configSplit3op(MMG5_pTetra pt,int vx[6],unsigned char tau[4],
                           const unsigned char **taued,
                           unsigned char sym[4],  unsigned char symed[6],
                           unsigned char *ip0,unsigned char *ip1,
                           unsigned char *ip2,unsigned char *ip3,
                           unsigned char *ie0,unsigned char *ie1,
                           unsigned char *ie2,unsigned char *ie3,
                           unsigned char *ie4,unsigned char *ie5,
                           unsigned char *imin03,unsigned char *imin12) {

  /* Set permutation /symmetry of vertices : generic case : 35 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  (*taued) = &MMG5_permedge[0][0];

  sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
  symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
  symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;

  switch(pt->flag) {
  case 19:
    tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
    (*taued) = &MMG5_permedge[0][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 13:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    (*taued) = &MMG5_permedge[2][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 37:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    (*taued) = &MMG5_permedge[2][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 22:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    (*taued) = &MMG5_permedge[10][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 28:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    (*taued) = &MMG5_permedge[10][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 26:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    (*taued) = &MMG5_permedge[6][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 14:
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    (*taued) = &MMG5_permedge[1][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 49:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    (*taued) = &MMG5_permedge[11][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 50:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    (*taued) = &MMG5_permedge[11][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;

  case 44:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    (*taued) = &MMG5_permedge[9][0];

    sym[0] = 0;  sym[1] = 1 ; sym[2] = 2 ; sym[3] = 3;
    symed[0] = 0 ; symed[1] = 1 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 4 ; symed[5] = 5;
    break;

  case 41:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    (*taued) = &MMG5_permedge[4][0];

    sym[0] = 0;  sym[1] = 2 ; sym[2] = 1 ; sym[3] = 3;
    symed[0] = 1 ; symed[1] = 0 ; symed[2] = 2;
    symed[3] = 3 ; symed[4] = 5 ; symed[5] = 4;
    break;
  }

  (*ip0) = tau[sym[0]];
  (*ip1) = tau[sym[1]];
  (*ip2) = tau[sym[2]];
  (*ip3) = tau[sym[3]];

  (*ie0) = (*taued)[symed[0]];
  (*ie1) = (*taued)[symed[1]];
  (*ie2) = (*taued)[symed[2]];
  (*ie3) = (*taued)[symed[3]];
  (*ie4) = (*taued)[symed[4]];
  (*ie5) = (*taued)[symed[5]];

  /* Test : to be removed eventually */
  assert(vx[(*ie0)] > 0);
  assert(vx[(*ie1)] > 0);
  assert(vx[(*ie5)] > 0);
  assert(vx[(*ie2)] <= 0);
  assert(vx[(*ie3)] <= 0);
  assert(vx[(*ie4)] <= 0);

  (*imin03) = (pt->v[(*ip0)] < pt->v[(*ip3)]) ? (*ip0) : (*ip3);
  (*imin12) = (pt->v[(*ip1)] < pt->v[(*ip2)]) ? (*ip1) : (*ip2);

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 *  Simulate split of 3 edges in opposite configuration.
 *
 */
int MMG3D_split3op_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4],sym[4],symed[6],ip0,ip1,ip2,ip3,ie0,ie1,ie2,ie3;
  unsigned char       ie4,ie5,imin03,imin12;
  const unsigned char *taued=NULL;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* Set permutation /symmetry of vertices : generic case : 35 */
  MMG3D_configSplit3op(pt,vx,tau,&taued,sym,symed,&ip0,&ip1,&ip2,&ip3,
                        &ie0,&ie1,&ie2,&ie3,&ie4,&ie5,&imin03,&imin12);

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  if ( (imin12 == ip2) && (imin03 == ip0) ) {
    pt0->v[ip0] = vx[ie1] ;  pt0->v[ip1] = vx[ie0] ; pt0->v[ip3] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie0] ; pt0->v[ip3] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie0] ; pt0->v[ip2] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip1] = vx[ie0] ; pt0->v[ip2] = vx[ie1] ; pt0->v[ip3] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip1] = vx[ie0] ; pt0->v[ip2] = vx[ie5];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }

  else if ( (imin12 == ip1) && (imin03 == ip0) ) {
    pt0->v[ip0] = vx[ie1] ; pt0->v[ip3] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie0] ; pt0->v[ip2] = vx[ie1] ; pt0->v[ip3] = vx[ie5];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie0] ; pt0->v[ip2] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip1] = vx[ie0] ; pt0->v[ip2] = vx[ie5];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip1] = vx[ie0] ; pt0->v[ip2] = vx[ie1]; pt0->v[ip3] = vx[ie5];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else if ( (imin12 == ip2) && (imin03 == ip3) ) {
    pt0->v[ip1] = vx[ie0] ; pt0->v[ip2] = vx[ie1] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie1] ; pt0->v[ip1] = vx[ie0] ; pt0->v[ip2] = vx[ie5];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie0] ; pt0->v[ip2] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie1] ; pt0->v[ip1] = vx[ie0]; pt0->v[ip3] = vx[ie5];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie0] ; pt0->v[ip3] = vx[ie5];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else {
    assert((imin12 == ip1) && (imin03 == ip3)) ;

    pt0->v[ip1] = vx[ie0] ; pt0->v[ip2] = vx[ie1] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie1] ; pt0->v[ip3] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie0] ; pt0->v[ip2] = vx[ie1] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[ip0] = vx[ie1] ; pt0->v[ip2] = vx[ie5] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split 3 opposite edges in a tetra
 *
 */
int MMG5_split3op(MMG5_pMesh mesh, MMG5_pSol met, int k, int vx[6],char metRidTyp){
  MMG5_pTetra          pt[5];
  MMG5_xTetra          xt[5];
  MMG5_pxTetra         pxt0;
  int                  iel;
  int                  newtet[5];
  unsigned char        imin12,imin03,tau[4],sym[4],symed[6],ip0,ip1,ip2,ip3,ie0,ie1;
  unsigned char        ie2,ie3,ie4,ie5,isxt[5],firstxt,i;
  const unsigned char *taued=NULL;

  pt[0]  = &mesh->tetra[k];
  newtet[0]=k;

  // To avoid warning about potentially uninitialized value for newtet[4]
  newtet[4] = 0;

  /* Set permutation /symmetry of vertices : generic case : 35 */
  MMG3D_configSplit3op(pt[0],vx,tau,&taued,sym,symed,&ip0,&ip1,&ip2,&ip3,
                        &ie0,&ie1,&ie2,&ie3,&ie4,&ie5,&imin03,&imin12);
  pt[0]->flag  = 0;

  /* Create new elements according to the current configuration */
  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
  }

  pt[1] = &mesh->tetra[iel];
  pt[1] = memcpy(pt[1],pt[0],sizeof(MMG5_Tetra));
  newtet[1]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate a"
                                " new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
  }
  pt[2] = &mesh->tetra[iel];
  pt[2] = memcpy(pt[2],pt[0],sizeof(MMG5_Tetra));
  newtet[2]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
    pt[2] = &mesh->tetra[newtet[2]];
  }
  pt[3] = &mesh->tetra[iel];
  pt[3] = memcpy(pt[3],pt[0],sizeof(MMG5_Tetra));
  newtet[3]=iel;

  if ( (pt[0])->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    memcpy(&xt[0],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[1],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[2],pxt0, sizeof(MMG5_xTetra));
    memcpy(&xt[3],pxt0, sizeof(MMG5_xTetra));
  }
  else {
    pxt0 = 0;
    memset(&xt[0],0, sizeof(MMG5_xTetra));
    memset(&xt[1],0, sizeof(MMG5_xTetra));
    memset(&xt[2],0, sizeof(MMG5_xTetra));
    memset(&xt[3],0, sizeof(MMG5_xTetra));
  }

  if ( !((imin12 == ip1) && (imin03 == ip3)) ) {
    iel = MMG3D_newElt(mesh);
    if ( !iel ) {
      MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                  " a new element.\n",__func__);
                          MMG5_INCREASE_MEM_MESSAGE();
                          fprintf(stderr,"  Exit program.\n");
                          return 0);
      pt[0] = &mesh->tetra[newtet[0]];
      pt[1] = &mesh->tetra[newtet[1]];
      pt[2] = &mesh->tetra[newtet[2]];
      pt[3] = &mesh->tetra[newtet[3]];
    }
    pt[4] = &mesh->tetra[iel];
    pt[4] = memcpy(pt[4],pt[0],sizeof(MMG5_Tetra));
    newtet[4]=iel;

    if ( pt[0]->xt ) {
      pxt0 = &mesh->xtetra[(pt[0])->xt];
      memcpy(&xt[4],pxt0, sizeof(MMG5_xTetra));
    }

    else {
      pxt0 = 0;
      memset(&xt[4],0, sizeof(MMG5_xTetra));
    }
  }

  /* Generic formulation of split of 3 edges in op configuration (edges 0,1,5 splitted) */
  if ( (imin12 == ip2) && (imin03 == ip0) ) {
    pt[0]->v[ip0] = vx[ie1] ;  pt[0]->v[ip1] = vx[ie0] ; pt[0]->v[ip3] = vx[ie5] ;
    xt[0].tag[ie0] = 0;  xt[0].tag[ie2] = 0;
    xt[0].tag[ie3] = 0;  xt[0].tag[ie4] = 0;
    xt[0].edg[ie0] = 0;  xt[0].edg[ie2] = 0;
    xt[0].edg[ie3] = 0;  xt[0].edg[ie4] = 0;
    xt[0].ref [ip0] = 0 ; xt[0].ref [ip2] = 0 ;
    xt[0].ftag[ip0] = 0 ; xt[0].ftag[ip2] = 0 ;
    MG_SET(xt[0].ori, ip0); MG_SET(xt[0].ori, ip2);

    pt[1]->v[ip0] = vx[ie0] ; pt[1]->v[ip3] = vx[ie5] ;
    xt[1].tag[ie1] = 0;  xt[1].tag[ie2] = 0;
    xt[1].tag[ie4] = 0;  xt[1].edg[ie1] = 0;
    xt[1].edg[ie2] = 0;  xt[1].edg[ie4] = 0;
    xt[1].ref [ip1] = 0 ; xt[1] .ref[ip2] = 0 ;
    xt[1].ftag[ip1] = 0 ; xt[1].ftag[ip2] = 0 ;
    MG_SET(xt[1].ori, ip1); MG_SET(xt[1].ori, ip2);

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie5] ;
    xt[2].tag[ie1] = 0;  xt[2].tag[ie2] = 0;
    xt[2].tag[ie3] = 0;  xt[2].edg[ie2] = 0;
    xt[2].edg[ie3] = 0;
    xt[2].ref [ip1] = 0 ; xt[2].ref [ip3] = 0 ;
    xt[2].ftag[ip1] = 0 ; xt[2].ftag[ip3] = 0 ;
    MG_SET(xt[2].ori, ip1); MG_SET(xt[2].ori, ip3);

    pt[3]->v[ip1] = vx[ie0] ; pt[3]->v[ip2] = vx[ie1] ; pt[3]->v[ip3] = vx[ie5] ;
    xt[3].tag[ie2] = 0;  xt[3].tag[ie3] = 0;
    xt[3].tag[ie4] = 0;  xt[3].tag[ie5] = 0;
    xt[3].edg[ie2] = 0;  xt[3].edg[ie3] = 0;
    xt[3].edg[ie4] = 0;  xt[3].edg[ie5] = 0;
    xt[3].ref [ip0] = 0 ; xt[3].ref [ip2] = 0 ;
    xt[3].ftag[ip0] = 0 ; xt[3].ftag[ip2] = 0 ;
    MG_SET(xt[3].ori, ip0); MG_SET(xt[3].ori, ip2);

    pt[4]->v[ip1] = vx[ie0] ; pt[4]->v[ip2] = vx[ie5];
    xt[4].tag[ie1] = 0;  xt[4].tag[ie3] = 0;
    xt[4].tag[ie4] = 0;  xt[4].edg[ie1] = 0;
    xt[4].edg[ie3] = 0;  xt[4].edg[ie4] = 0;
    xt[4].ref [ip0] = 0 ; xt[4].ref [ip3] = 0 ;
    xt[4].ftag[ip0] = 0 ; xt[4].ftag[ip3] = 0 ;
    MG_SET(xt[4].ori, ip0); MG_SET(xt[4].ori, ip3);
  }

  else if ( (imin12 == ip1) && (imin03 == ip0) ) {
    pt[0]->v[ip0] = vx[ie1] ; pt[0]->v[ip3] = vx[ie5] ;
    xt[0].tag[ie0] = 0;  xt[0].tag[ie2] = 0;
    xt[0].tag[ie4] = 0;  xt[0].edg[ie0] = 0;
    xt[0].edg[ie2] = 0;  xt[0].edg[ie4] = 0;
    xt[0].ref[ip2]  = 0 ;
    xt[0].ftag[ip2] = 0 ;
    MG_SET(xt[0].ori, ip2);

    pt[1]->v[ip0] = vx[ie0] ; pt[1]->v[ip2] = vx[ie1] ; pt[1]->v[ip3] = vx[ie5];
    xt[1].tag[ie1] = 0;  xt[1].tag[ie2] = 0;
    xt[1].tag[ie3] = 0;  xt[1].tag[ie4] = 0;
    xt[1].tag[ie5] = 0;  xt[1].edg[ie1] = 0;
    xt[1].edg[ie2] = 0;  xt[1].edg[ie3] = 0;
    xt[1].edg[ie4] = 0;  xt[1].edg[ie5] = 0;
    xt[1].ref [ip0] = 0 ; xt[1].ref [ip1] = 0 ; xt[1].ref [ip2] = 0 ;
    xt[1].ftag[ip0] = 0 ; xt[1].ftag[ip1] = 0 ; xt[1].ftag[ip2] = 0 ;
    MG_SET(xt[1].ori, ip0); MG_SET(xt[1].ori, ip1); MG_SET(xt[1].ori, ip2);

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie5] ;
    xt[2].tag[ie1] = 0;  xt[2].tag[ie2] = 0;
    xt[2].tag[ie3] = 0;  xt[2].edg[ie1] = 0;
    xt[2].edg[ie2] = 0;  xt[2].edg[ie3] = 0;
    xt[2].ref [ip1] = 0 ; xt[2].ref [ip3] = 0 ;
    xt[2].ftag[ip1] = 0 ; xt[2].ftag[ip3] = 0 ;
    MG_SET(xt[2].ori, ip1); MG_SET(xt[2].ori, ip3);

    pt[3]->v[ip1] = vx[ie0] ; pt[3]->v[ip2] = vx[ie5];
    xt[3].tag[ie1] = 0;  xt[3].tag[ie3] = 0;
    xt[3].tag[ie4] = 0;  xt[3].edg[ie1] = 0;
    xt[3].edg[ie3] = 0;  xt[3].edg[ie4] = 0;
    xt[3].ref [ip0] = 0 ; xt[3].ref [ip3] = 0 ;
    xt[3].ftag[ip0] = 0 ; xt[3].ftag[ip3] = 0 ;
    MG_SET(xt[3].ori, ip0); MG_SET(xt[3].ori, ip3);

    pt[4]->v[ip1] = vx[ie0] ; pt[4]->v[ip2] = vx[ie1]; pt[4]->v[ip3] = vx[ie5];
    xt[4].tag[ie2] = 0;  xt[4].tag[ie3] = 0;
    xt[4].tag[ie4] = 0;  xt[4].tag[ie5] = 0;
    xt[4].edg[ie2] = 0;  xt[4].edg[ie3] = 0;
    xt[4].edg[ie4] = 0;  xt[4].edg[ie5] = 0;
    xt[4].ref [ip0] = 0 ; xt[4].ref [ip2] = 0 ;
    xt[4].ftag[ip0] = 0 ; xt[4].ftag[ip2] = 0 ;
    MG_SET(xt[4].ori, ip0); MG_SET(xt[4].ori, ip2);
  }

  else if ( (imin12 == ip2) && (imin03 == ip3) ) {
    pt[0]->v[ip1] = vx[ie0] ; pt[0]->v[ip2] = vx[ie1] ;
    xt[0].tag[ie3] = 0;  xt[0].tag[ie4] = 0;
    xt[0].tag[ie5] = 0;  xt[0].edg[ie3] = 0;
    xt[0].edg[ie4] = 0;  xt[0].edg[ie5] = 0;
    xt[0].ref[ip0]  = 0 ;
    xt[0].ftag[ip0] = 0 ;
    MG_SET(xt[0].ori, ip0);

    pt[1]->v[ip0] = vx[ie1] ; pt[1]->v[ip1] = vx[ie0] ; pt[1]->v[ip2] = vx[ie5];
    xt[1].tag[ie0] = 0;  xt[1].tag[ie1] = 0;
    xt[1].tag[ie2] = 0;  xt[1].tag[ie3] = 0;
    xt[1].tag[ie4] = 0;  xt[1].edg[ie0] = 0;
    xt[1].edg[ie1] = 0;  xt[1].edg[ie2] = 0;
    xt[1].edg[ie3] = 0;  xt[1].edg[ie4] = 0;
    xt[1].ref [ip0] = 0 ; xt[1].ref [ip2] = 0 ; xt[1].ref [ip3] = 0 ;
    xt[1].ftag[ip0] = 0 ; xt[1].ftag[ip2] = 0 ; xt[1].ftag[ip3] = 0 ;
    MG_SET(xt[1].ori, ip1); MG_SET(xt[1].ori, ip2); MG_SET(xt[1].ori, ip3);

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie5] ;
    xt[2].tag[ie1] = 0;  xt[2].tag[ie2] = 0;
    xt[2].tag[ie3] = 0;  xt[2].edg[ie1] = 0;
    xt[2].edg[ie2] = 0;  xt[2].edg[ie3] = 0;
    xt[2].ref [ip1] = 0 ; xt[2].ref [ip3] = 0 ;
    xt[2].ftag[ip1] = 0 ; xt[2].ftag[ip3] = 0 ;
    MG_SET(xt[2].ori, ip1); MG_SET(xt[2].ori, ip3);

    pt[3]->v[ip0] = vx[ie1] ; pt[3]->v[ip1] = vx[ie0]; pt[3]->v[ip3] = vx[ie5];
    xt[3].tag[ie0] = 0;  xt[3].tag[ie2] = 0;
    xt[3].tag[ie3] = 0;  xt[3].tag[ie4] = 0;
    xt[3].edg[ie0] = 0;  xt[3].edg[ie2] = 0;
    xt[3].edg[ie3] = 0;  xt[3].edg[ie4] = 0;
    xt[3].ref [ip0] = 0 ; xt[3].ref [ip2] = 0 ;
    xt[3].ftag[ip0] = 0 ; xt[3].ftag[ip2] = 0 ;
    MG_SET(xt[3].ori, ip0); MG_SET(xt[3].ori, ip2);

    pt[4]->v[ip0] = vx[ie0] ; pt[4]->v[ip3] = vx[ie5];
    xt[4].tag[ie1] = 0;  xt[4].tag[ie2] = 0;
    xt[4].tag[ie4] = 0;  xt[4].edg[ie1] = 0;
    xt[4].edg[ie2] = 0;  xt[4].edg[ie4] = 0;
    xt[4].ref [ip1] = 0 ; xt[4].ref [ip2] = 0 ;
    xt[4].ftag[ip1] = 0 ; xt[4].ftag[ip2] = 0 ;
    MG_SET(xt[4].ori, ip1); MG_SET(xt[4].ori, ip2);
  }
  else {
    assert((imin12 == ip1) && (imin03 == ip3)) ;

    pt[0]->v[ip1] = vx[ie0] ; pt[0]->v[ip2] = vx[ie1] ;
    xt[0].tag[ie3] = 0;  xt[0].tag[ie4] = 0;
    xt[0].tag[ie5] = 0;  xt[0].edg[ie3] = 0;
    xt[0].edg[ie4] = 0;  xt[0].edg[ie5] = 0;
    xt[0].ref [ip0] = 0 ;
    xt[0].ftag[ip0] = 0 ;
    MG_SET(xt[0].ori, ip0);

    pt[1]->v[ip0] = vx[ie1] ; pt[1]->v[ip3] = vx[ie5] ;
    xt[1].tag[ie0] = 0;  xt[1].tag[ie2] = 0;
    xt[1].tag[ie4] = 0;  xt[1].edg[ie0] = 0;
    xt[1].edg[ie2] = 0;  xt[1].edg[ie4] = 0;
    xt[1].ref [ip2] = 0 ;
    xt[1].ftag[ip2] = 0 ;
    MG_SET(xt[1].ori, ip2);

    pt[2]->v[ip0] = vx[ie0] ; pt[2]->v[ip2] = vx[ie1] ;
    xt[2].tag[ie1] = 0;  xt[2].tag[ie2] = 0;
    xt[2].tag[ie3] = 0;  xt[2].tag[ie5] = 0;
    xt[2].edg[ie1] = 0;  xt[2].edg[ie2] = 0;
    xt[2].edg[ie3] = 0;  xt[2].edg[ie5] = 0;
    xt[2].ref [ip0] = 0 ; xt[2].ref [ip1] = 0 ;
    xt[2].ftag[ip0] = 0 ; xt[2].ftag[ip1] = 0 ;
    MG_SET(xt[2].ori, ip0); MG_SET(xt[2].ori, ip1);

    pt[3]->v[ip0] = vx[ie1] ; pt[3]->v[ip2] = vx[ie5] ;
    xt[3].tag[ie0] = 0;  xt[3].tag[ie1] = 0;
    xt[3].tag[ie2] = 0;  xt[3].tag[ie3] = 0;
    xt[3].edg[ie0] = 0;  xt[3].edg[ie1] = 0;
    xt[3].edg[ie2] = 0;  xt[3].edg[ie3] = 0;
    xt[3].ref [ip2] = 0 ; xt[3].ref [ip3] = 0 ;
    xt[3].ftag[ip2] = 0 ; xt[3].ftag[ip3] = 0 ;
    MG_SET(xt[3].ori, ip2); MG_SET(xt[3].ori, ip3);
  }

  /* Assignation of the xt fields to the appropriate tets */
  if ( (imin12 == ip1) && (imin03 == ip3) ) {
    isxt[0] = isxt[1] = isxt[2] = isxt[3] = 0;

    for (i=0; i<4; i++) {
      if ( (xt[0]).ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
      if ( (xt[1]).ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
      if ( (xt[2]).ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
      if ( (xt[3]).ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
    }

    if ( pt[0]->xt ) {
      if ( isxt[0] ) {
        memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
        pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;

        for (i=1; i<4; i++) {
          if ( isxt[i] ) {
            mesh->xt++;
            if ( mesh->xt >= mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
      }
      else {
        firstxt = 1;
        pt[1]->xt = pt[2]->xt = pt[3]->xt = 0;

        for (i=1; i<4; i++) {
          if ( isxt[i] ) {
            if ( firstxt ) {
              firstxt = 0;
              pt[i]->xt = pt[0]->xt;
              pxt0 = &mesh->xtetra[(pt[i])->xt];
              memcpy(pxt0,&(xt[i]),sizeof(MMG5_xTetra));
            }
            else {
              mesh->xt++;
              if ( mesh->xt > mesh->xtmax ) {
                /* realloc of xtetras table */
                MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                   "larger xtetra table",
                                   mesh->xt--;
                                   fprintf(stderr,"  Exit program.\n");
                                   return 0);
              }
              pt[i]->xt = mesh->xt;
              pxt0 = &mesh->xtetra[mesh->xt];
              memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
            }
          }
        }
        pt[0]->xt = 0;
      }
    }

  }
  else {
    isxt[0] = isxt[1] = isxt[2] = isxt[3] = isxt[4] = 0;

    for (i=0; i<4; i++) {
      if ( (xt[0]).ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
      if ( (xt[1]).ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
      if ( (xt[2]).ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
      if ( (xt[3]).ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
      if ( (xt[4]).ref[i] || xt[4].ftag[i] ) isxt[4] = 1;
    }

    if ( pt[0]->xt ) {
      if ( isxt[0] ) {
        memcpy(pxt0,&(xt[0]),sizeof(MMG5_xTetra));
        pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = 0;

        for(i=1; i<5; i++) {
          if ( isxt[i] ) {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
      }
      else {
        firstxt = 1;
        pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = 0;

        for (i=1; i<5; i++) {
          if ( isxt[i] ) {
            if ( firstxt ) {
              firstxt = 0;
              pt[i]->xt = pt[0]->xt;
              pxt0 = &mesh->xtetra[pt[i]->xt];
              memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
            }
            else {
              mesh->xt++;
              if ( mesh->xt > mesh->xtmax ) {
                /* realloc of xtetras table */
                MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                   "larger xtetra table",
                                   mesh->xt--;
                                   fprintf(stderr,"  Exit program.\n");
                                   return 0);
              }
              pt[i]->xt = mesh->xt;
              pxt0 = &mesh->xtetra[mesh->xt];
              memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
            }
          }
        }
        pt[0]->xt = 0;
      }
    }
  }
  /* Quality update */
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    pt[0]->qual=MMG5_caltet33_ani(mesh,met,pt[0]);
    pt[1]->qual=MMG5_caltet33_ani(mesh,met,pt[1]);
    pt[2]->qual=MMG5_caltet33_ani(mesh,met,pt[2]);
    pt[3]->qual=MMG5_caltet33_ani(mesh,met,pt[3]);
    if ( !((imin12 == ip1) && (imin03 == ip3)) ) {
      pt[4]->qual=MMG5_caltet33_ani(mesh,met,pt[4]);
    }
  }
  else {
    pt[0]->qual=MMG5_orcal(mesh,met,newtet[0]);
    pt[1]->qual=MMG5_orcal(mesh,met,newtet[1]);
    pt[2]->qual=MMG5_orcal(mesh,met,newtet[2]);
    pt[3]->qual=MMG5_orcal(mesh,met,newtet[3]);
    if ( !((imin12 == ip1) && (imin03 == ip3)) ) {
      pt[4]->qual=MMG5_orcal(mesh,met,newtet[4]);
    }
  }
  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k tetra index.
 * \param metRidTyp metric storage (classic or special)
 * \return 0 if fail, index of created point otherwise (\a ib)
 *
 * Split a tetra in 4 tetras by introducing its barycenter. FOR NOW : flags,
 * that tell which edge should be split, are not updated (erased) : UPDATE
 * NEEDED ?
 *
 */
int MMG5_split4bar(MMG5_pMesh mesh, MMG5_pSol met, int k,char metRidTyp) {
  MMG5_pTetra   pt[4];
  MMG5_pPoint   ppt;
  MMG5_xTetra   xt[4];
  MMG5_pxTetra  pxt0;
  double        o[3],cb[4];
  int           i,ib,iel,iadr,*adja,adj1,adj2,adj3;
  int           newtet[4];
  unsigned char isxt[4],firstxt;

  pt[0] = &mesh->tetra[k];
  pt[0]->flag = 0;
  newtet[0]=k;

  o[0] = o[1] = o[2] = 0.0;
  for (i=0; i<4; i++) {
    ib    = pt[0]->v[i];
    ppt   = &mesh->point[ib];
    o[0] += ppt->c[0];
    o[1] += ppt->c[1];
    o[2] += ppt->c[2];
  }
  o[0] *= 0.25;
  o[1] *= 0.25;
  o[2] *= 0.25;

  cb[0] = 0.25; cb[1] = 0.25;  cb[2] = 0.25;  cb[3] = 0.25;
  ib = MMG3D_newPt(mesh,o,0);
  if ( !ib ) {
    MMG3D_POINT_REALLOC(mesh,met,ib,mesh->gap,
                         fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                 " a new point\n",__func__);
                         MMG5_INCREASE_MEM_MESSAGE();
                         return 0
                         ,o,0);
  }
  if ( met->m ) {
    if ( !metRidTyp && met->size > 1 )
      MMG5_interp4bar33_ani(mesh,met,k,ib,cb);
    else
      MMG5_interp4bar(mesh,met,k,ib,cb);
  }

  /* create 3 new tetras */
  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        MMG3D_delPt(mesh,ib);
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
  }
  pt[1] = &mesh->tetra[iel];
  pt[1] = memcpy(pt[1],pt[0],sizeof(MMG5_Tetra));
  newtet[1]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        MMG3D_delPt(mesh,ib);
                        if ( !MMG3D_delElt(mesh,newtet[1]) ) return 0;
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
  }
  pt[2] = &mesh->tetra[iel];
  pt[2] = memcpy(pt[2],pt[0],sizeof(MMG5_Tetra));
  newtet[2]=iel;

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        MMG3D_delPt(mesh,ib);
                        if ( !MMG3D_delElt(mesh,newtet[1]) ) return 0;
                        if ( !MMG3D_delElt(mesh,newtet[2]) ) return 0;
                        return 0);
    pt[0] = &mesh->tetra[newtet[0]];
    pt[1] = &mesh->tetra[newtet[1]];
    pt[2] = &mesh->tetra[newtet[2]];
  }
  pt[3] = &mesh->tetra[iel];
  pt[3] = memcpy(pt[3],pt[0],sizeof(MMG5_Tetra));
  newtet[3]=iel;

  memset(&xt[0],0, sizeof(MMG5_xTetra));
  memset(&xt[1],0, sizeof(MMG5_xTetra));
  memset(&xt[2],0, sizeof(MMG5_xTetra));
  memset(&xt[3],0, sizeof(MMG5_xTetra));
  pxt0 = 0;
  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[pt[0]->xt];
    memcpy(&xt[0],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[1],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[2],pxt0,sizeof(MMG5_xTetra));
    memcpy(&xt[3],pxt0,sizeof(MMG5_xTetra));
  }

  /* Update adjacency */
  if ( mesh->adja ) {
    iadr  = 4*(newtet[0]-1)+1;
    adja  = &mesh->adja[iadr];

    /* Store the old adjacents */
    adj1 = mesh->adja[iadr+1];
    adj2 = mesh->adja[iadr+2];
    adj3 = mesh->adja[iadr+3];

    /* Update the new ones */
    adja[1] = 4*newtet[1];
    adja[2] = 4*newtet[2];
    adja[3] = 4*newtet[3];

    iadr    = 4*(newtet[1]-1)+1;
    adja    = &mesh->adja[iadr];
    adja[0] = 4*newtet[0] + 1;
    adja[1] = adj1;
    adja[2] = 4*newtet[2] + 1;
    adja[3] = 4*newtet[3] + 1;
    if ( adj1 )
      mesh->adja[4*(adj1/4-1) + 1+adj1%4] = 4*newtet[1]+1;

    iadr    = 4*(newtet[2]-1)+1;
    adja    = &mesh->adja[iadr];
    adja[0] = 4*newtet[0] + 2;
    adja[1] = 4*newtet[1] + 2;
    adja[2] = adj2;
    adja[3] = 4*newtet[3] + 2;
    if ( adj2 )
      mesh->adja[4*(adj2/4-1) + 1+adj2%4] = 4*newtet[2]+2;

    iadr    = 4*(newtet[3]-1)+1;
    adja    = &mesh->adja[iadr];
    adja[0] = 4*newtet[0] + 3;
    adja[1] = 4*newtet[1] + 3;
    adja[2] = 4*newtet[2] + 3;
    adja[3] = adj3;
    if ( adj3 )
      mesh->adja[4*(adj3/4-1) + 1+adj3%4] = 4*newtet[3]+3;
  }

  /* Update vertices and xt fields */
  pt[0]->v[0] = pt[1]->v[1] = pt[2]->v[2] = pt[3]->v[3] = ib;

  xt[0].tag[0]  = 0;  xt[0].edg[0]  = 0;
  xt[0].tag[1]  = 0;  xt[0].edg[1]  = 0;
  xt[0].tag[2]  = 0;  xt[0].edg[2]  = 0;
  xt[0].ref [1] = 0;  xt[0].ref [2] = 0;  xt[0].ref [3] = 0;
  xt[0].ftag[1] = 0;  xt[0].ftag[2] = 0;  xt[0].ftag[3] = 0;
  MG_SET(xt[0].ori, 1);  MG_SET(xt[0].ori, 2);  MG_SET(xt[0].ori, 3);

  xt[1].tag[0]  = 0;  xt[1].edg[0]  = 0;
  xt[1].tag[3]  = 0;  xt[1].edg[3]  = 0;
  xt[1].tag[4]  = 0;  xt[1].edg[4]  = 0;
  xt[1].ref [0] = 0;  xt[1].ref [2] = 0;  xt[1].ref [3] = 0;
  xt[1].ftag[0] = 0;  xt[1].ftag[2] = 0;  xt[1].ftag[3] = 0;
  MG_SET(xt[1].ori, 0);  MG_SET(xt[1].ori, 2);  MG_SET(xt[1].ori, 3);

  xt[2].tag[1]  = 0;  xt[2].edg[1]  = 0;
  xt[2].tag[3]  = 0;  xt[2].edg[3]  = 0;
  xt[2].tag[5]  = 0;  xt[2].edg[5]  = 0;
  xt[2].ref [0] = 0;  xt[2].ref [1] = 0;  xt[2].ref [3] = 0;
  xt[2].ftag[0] = 0;  xt[2].ftag[1] = 0;  xt[2].ftag[3] = 0;
  MG_SET(xt[2].ori, 0);  MG_SET(xt[2].ori, 1);  MG_SET(xt[2].ori, 3);

  xt[3].tag[2]  = 0;  xt[3].edg[2]  = 0;
  xt[3].tag[4]  = 0;  xt[3].edg[4]  = 0;
  xt[3].tag[5]  = 0;  xt[3].edg[5]  = 0;
  xt[3].ref [0] = 0;  xt[3].ref [1] = 0;  xt[3].ref [2] = 0;
  xt[3].ftag[0] = 0;  xt[3].ftag[1] = 0;  xt[3].ftag[2] = 0;
  MG_SET(xt[3].ori, 0);  MG_SET(xt[3].ori, 1);  MG_SET(xt[3].ori, 2);

  /* Assignation of the xt fields to the appropriate tets */
  memset(isxt,0,4*sizeof(char));
  for (i=0; i<4; i++) {
    if ( xt[0].ref[i] || xt[0].ftag[i] ) isxt[0] = 1;
    if ( xt[1].ref[i] || xt[1].ftag[i] ) isxt[1] = 1;
    if ( xt[2].ref[i] || xt[2].ftag[i] ) isxt[2] = 1;
    if ( xt[3].ref[i] || xt[3].ftag[i] ) isxt[3] = 1;
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
        }
        else {
          pt[i]->xt = 0;
        }
      }
    }
    else {
      firstxt = 1;
      for (i=1; i<4; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
        else {
          pt[i]->xt = 0;
        }
      }
      pt[0]->xt = 0;
    }
  }
  /* Quality update */
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    pt[0]->qual=MMG5_caltet33_ani(mesh,met,pt[0]);
    pt[1]->qual=MMG5_caltet33_ani(mesh,met,pt[1]);
    pt[2]->qual=MMG5_caltet33_ani(mesh,met,pt[2]);
    pt[3]->qual=MMG5_caltet33_ani(mesh,met,pt[3]);
  }
  else {
    pt[0]->qual=MMG5_orcal(mesh,met,newtet[0]);
    pt[1]->qual=MMG5_orcal(mesh,met,newtet[1]);
    pt[2]->qual=MMG5_orcal(mesh,met,newtet[2]);
    pt[3]->qual=MMG5_orcal(mesh,met,newtet[3]);
  }

  return ib;
}

/**
 * \param pt initial tetra
 * \param vx index of points to insert along edges
 * \param tau vertices permutation
 * \param taued edges permutation
 * \param imin23 minimal index of vertices ip0 and ip3
 * \param imin12 minimal index of vertices ip1 and ip2
 *
 * Set permutation of vertices for the split of 4 edges when 3 lie on the same
 * face. Reference configuration 23
 *
 */
static inline
void MMG3D_configSplit4sf(MMG5_pTetra pt,int vx[6],unsigned char tau[4],
                           const unsigned char **taued,
                           unsigned char *imin23,unsigned char *imin12) {

  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  (*taued) = &MMG5_permedge[0][0];
  switch(pt->flag){
  case 29:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    (*taued) = &MMG5_permedge[5][0];
    break;

  case 53:
    tau[0] = 3 ; tau[1] = 0 ; tau[2] = 2 ; tau[3] = 1;
    (*taued) = &MMG5_permedge[9][0];
    break;

  case 60:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    (*taued) = &MMG5_permedge[10][0];
    break;

  case 57:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    (*taued) = &MMG5_permedge[4][0];
    break;

  case 58:
    tau[0] = 2 ; tau[1] = 3 ; tau[2] = 0 ; tau[3] = 1;
    (*taued) = &MMG5_permedge[8][0];
    break;

  case 27:
    tau[0] = 1 ; tau[1] = 0 ; tau[2] = 3 ; tau[3] = 2;
    (*taued) = &MMG5_permedge[3][0];
    break;

  case 15:
    tau[0] = 0 ; tau[1] = 2 ; tau[2] = 3 ; tau[3] = 1;
    (*taued) = &MMG5_permedge[1][0];
    break;

  case 43:
    tau[0] = 2 ; tau[1] = 1 ; tau[2] = 3 ; tau[3] = 0;
    (*taued) = &MMG5_permedge[7][0];
    break;

  case 39:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    (*taued) = &MMG5_permedge[2][0];
    break;

  case 54:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    (*taued) = &MMG5_permedge[11][0];
    break;

  case 46:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    (*taued) = &MMG5_permedge[6][0];
    break;
  }

  (*imin23) = (pt->v[tau[2]] < pt->v[tau[3]]) ? tau[2] : tau[3];
  (*imin12) = (pt->v[tau[1]] < pt->v[tau[2]]) ? tau[1] : tau[2];
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 *  Simulate split of 4 edges in a configuration when 3 lie on the same face.
 *
 */
int MMG3D_split4sf_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4];
  unsigned char       imin23,imin12;
  const unsigned char *taued = NULL;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* Set permutation of vertices : reference configuration : 23 */
  MMG3D_configSplit4sf(pt,vx,tau,&taued,&imin23,&imin12);

  /* Generic formulation of split of 4 edges (with 3 on same face) */
  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[1]] = vx[taued[0]];
  pt0->v[tau[2]] = vx[taued[1]];
  pt0->v[tau[3]] = vx[taued[2]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[2]];
  pt0->v[tau[1]] = vx[taued[0]];
  pt0->v[tau[2]] = vx[taued[1]];
  pt0->v[tau[3]] = vx[taued[4]] ;
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  if ( imin12 == tau[1] ) {
    pt0->v[tau[0]] = vx[taued[0]];
    pt0->v[tau[2]] = vx[taued[1]];
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[1]];
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else {
    pt0->v[tau[0]] = vx[taued[1]];
    pt0->v[tau[1]] = vx[taued[0]];
    pt0->v[tau[3]] = vx[taued[4]] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[0]];
    pt0->v[tau[3]] = vx[taued[4]] ;
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  if ( imin23 == tau[2] ) {
    pt0->v[tau[0]] = vx[taued[1]];
    pt0->v[tau[1]] = vx[taued[4]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[2]];
    pt0->v[tau[1]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else {
    pt0->v[tau[0]] = vx[taued[2]];
    pt0->v[tau[1]] = vx[taued[4]];
    pt0->v[tau[2]] = vx[taued[1]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[1]];
    pt0->v[tau[1]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split 4 edges in a configuration when 3 lie on the same face
 *
 */
int MMG5_split4sf(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra         pt[6];
  MMG5_xTetra         xt[6];
  MMG5_pxTetra        pxt0;
  int                 iel;
  int                 newtet[6];
  char                firstxt,isxt[6],j,i;
  unsigned char       tau[4],imin23,imin12;
  const unsigned char *taued = NULL;

  pt[0]  = &mesh->tetra[k];
  newtet[0]=k;

  /* Set permutation of vertices : reference configuration : 23 */
  MMG3D_configSplit4sf(pt[0],vx,tau,&taued,&imin23,&imin12);
  pt[0]->flag  = 0;

  /* create 5 new tetras */
  for (j=1; j<6; j++) {
    iel = MMG3D_newElt(mesh);
    if ( !iel ) {
      MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                  " a new element.\n",__func__);
                          MMG5_INCREASE_MEM_MESSAGE();
                          fprintf(stderr,"  Exit program.\n");
                          return 0);
      for ( i=0; i<j; i++)
        pt[i] = &mesh->tetra[newtet[i]];
    }
    pt[j] = &mesh->tetra[iel];
    pt[j] = memcpy(pt[j],pt[0],sizeof(MMG5_Tetra));
    newtet[j]=iel;
  }

  if ( (pt[0])->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    for (j=0; j<6; j++) {
      memcpy(&xt[j],pxt0, sizeof(MMG5_xTetra));
    }
  }
  else {
    pxt0 = 0;
    for (j=0; j<6; j++) {
      memset(&xt[j],0, sizeof(MMG5_xTetra));
    }
  }

  /* Generic formulation of split of 4 edges (with 3 on same face) */
  pt[0]->v[tau[1]] = vx[taued[0]] ;   pt[0]->v[tau[2]] = vx[taued[1]] ;   pt[0]->v[tau[3]] = vx[taued[2]];
  xt[0].tag[taued[3]] = 0;  xt[0].tag[taued[4]] = 0;
  xt[0].tag[taued[5]] = 0;  xt[0].edg[taued[3]] = 0;
  xt[0].edg[taued[4]] = 0;  xt[0].edg[taued[5]] = 0;
  xt[0].ref [ tau[0]] = 0 ;
  xt[0].ftag[ tau[0]] = 0 ;
  MG_SET(xt[0].ori, tau[0]);

  pt[1]->v[tau[0]] = vx[taued[2]] ; pt[1]->v[tau[1]] = vx[taued[0]] ;
  pt[1]->v[tau[2]] = vx[taued[1]] ; pt[1]->v[tau[3]] = vx[taued[4]] ;
  xt[1].tag[taued[0]] = 0;  xt[1].tag[taued[1]] = 0;
  xt[1].tag[taued[2]] = 0;  xt[1].tag[taued[3]] = 0;
  xt[1].tag[taued[4]] = 0;  xt[1].tag[taued[5]] = 0;
  xt[1].edg[taued[0]] = 0;  xt[1].edg[taued[1]] = 0;
  xt[1].edg[taued[2]] = 0;  xt[1].edg[taued[3]] = 0;
  xt[1].edg[taued[4]] = 0;  xt[1].edg[taued[5]] = 0;
  xt[1].ref [ tau[0]] = 0 ; xt[1].ref [ tau[1]] = 0 ; xt[1].ref [tau[3]] = 0 ;
  xt[1].ftag[ tau[0]] = 0 ; xt[1].ftag[ tau[1]] = 0 ; xt[1].ftag[tau[3]] = 0 ;
  MG_SET(xt[1].ori, tau[0]); MG_SET(xt[1].ori, tau[1]); MG_SET(xt[1].ori, tau[3]);

  if ( imin12 == tau[1] ) {
    pt[2]->v[tau[0]] = vx[taued[0]] ; pt[2]->v[tau[2]] = vx[taued[1]] ; pt[2]->v[tau[3]] = vx[taued[4]] ;
    xt[2].tag[taued[1]] = 0;  xt[2].tag[taued[2]] = 0;
    xt[2].tag[taued[3]] = 0;  xt[2].tag[taued[5]] = 0;
    xt[2].edg[taued[1]] = 0;  xt[2].edg[taued[2]] = 0;
    xt[2].edg[taued[3]] = 0;  xt[2].edg[taued[5]] = 0;
    xt[2].ref [ tau[0]] = 0 ; xt[2].ref [ tau[1]] = 0 ;
    xt[2].ftag[ tau[0]] = 0 ; xt[2].ftag[ tau[1]] = 0 ;
    MG_SET(xt[2].ori, tau[0]); MG_SET(xt[2].ori, tau[1]);

    pt[3]->v[tau[0]] = vx[taued[1]] ; pt[3]->v[tau[3]] = vx[taued[4]] ;
    xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[2]] = 0;
    xt[3].tag[taued[5]] = 0;  xt[3].edg[taued[0]] = 0;
    xt[3].edg[taued[2]] = 0;  xt[3].edg[taued[5]] = 0;
    xt[3].ref [ tau[1]] = 0 ; xt[3].ref [ tau[2]] = 0 ;
    xt[3].ftag[ tau[1]] = 0 ; xt[3].ftag[ tau[2]] = 0 ;
    MG_SET(xt[3].ori, tau[1]); MG_SET(xt[3].ori, tau[2]);
  }
  else {
    pt[2]->v[tau[0]] = vx[taued[1]] ; pt[2]->v[tau[1]] = vx[taued[0]] ; pt[2]->v[tau[3]] = vx[taued[4]] ;
    xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[2]] = 0;
    xt[2].tag[taued[3]] = 0;  xt[2].tag[taued[4]] = 0;
    xt[2].tag[taued[5]] = 0;  xt[2].edg[taued[0]] = 0;
    xt[2].edg[taued[2]] = 0;  xt[2].edg[taued[3]] = 0;
    xt[2].edg[taued[4]] = 0;  xt[2].edg[taued[5]] = 0;
    xt[2].ref [ tau[0]] = 0 ; xt[2].ref [ tau[1]] = 0 ; xt[2].ref [tau[2]] = 0 ;
    xt[2].ftag[ tau[0]] = 0 ; xt[2].ftag[ tau[1]] = 0 ; xt[2].ftag[tau[2]] = 0 ;
    MG_SET(xt[2].ori, tau[0]); MG_SET(xt[2].ori, tau[1]); MG_SET(xt[2].ori, tau[2]);

    pt[3]->v[tau[0]] = vx[taued[0]] ; pt[3]->v[tau[3]] = vx[taued[4]] ;
    xt[3].tag[taued[1]] = 0;  xt[3].tag[taued[2]] = 0;
    xt[3].tag[taued[5]] = 0;  xt[3].edg[taued[1]] = 0;
    xt[3].edg[taued[2]] = 0;  xt[3].edg[taued[5]] = 0;
    xt[3].ref [ tau[1]] = 0 ;
    xt[3].ftag[ tau[1]] = 0 ;
    MG_SET(xt[3].ori, tau[1]);
  }

  if ( imin23 == tau[2] ) {
    pt[4]->v[tau[0]] = vx[taued[1]] ; pt[4]->v[tau[1]] = vx[taued[4]] ; pt[4]->v[tau[3]] = vx[taued[2]] ;
    xt[4].tag[taued[0]] = 0;  xt[4].tag[taued[2]] = 0;
    xt[4].tag[taued[3]] = 0;  xt[4].tag[taued[4]] = 0;
    xt[4].tag[taued[5]] = 0;
    xt[4].edg[taued[0]] = 0;  xt[4].edg[taued[2]] = 0;
    xt[4].edg[taued[3]] = 0;  xt[4].edg[taued[4]] = 0;
    xt[4].edg[taued[5]] = 0;
    xt[4].ref [ tau[0]] = 0;  xt[4].ref [ tau[2]] = 0 ;
    xt[4].ref [ tau[3]] = 0 ;
    xt[4].ftag[ tau[0]] = 0;  xt[4].ftag[ tau[2]] = 0 ;
    xt[4].ftag[ tau[3]] = 0 ;
    MG_SET(xt[4].ori, tau[0]); MG_SET(xt[4].ori, tau[2]); MG_SET(xt[4].ori, tau[3]);

    pt[5]->v[tau[0]] = vx[taued[2]] ; pt[5]->v[tau[1]] = vx[taued[4]] ;
    xt[5].tag[taued[0]] = 0;  xt[5].tag[taued[1]] = 0;
    xt[5].tag[taued[3]] = 0;  xt[5].edg[taued[0]] = 0;
    xt[5].edg[taued[1]] = 0;  xt[5].edg[taued[3]] = 0;
    xt[5].ref [ tau[3]] = 0 ;
    xt[5].ftag[ tau[3]] = 0 ;
    MG_SET(xt[5].ori, tau[3]);
  }
  else {
    pt[4]->v[tau[0]] = vx[taued[2]] ; pt[4]->v[tau[1]] = vx[taued[4]] ; pt[4]->v[tau[2]] = vx[taued[1]] ;
    xt[4].tag[taued[0]] = 0;  xt[4].tag[taued[1]] = 0;
    xt[4].tag[taued[3]] = 0;  xt[4].tag[taued[5]] = 0;
    xt[4].edg[taued[0]] = 0;  xt[4].edg[taued[1]] = 0;
    xt[4].edg[taued[3]] = 0;  xt[4].edg[taued[5]] = 0;
    xt[4].ref [ tau[0]] = 0;  xt[4].ref [ tau[3]] = 0 ;
    xt[4].ftag[ tau[0]] = 0;  xt[4].ftag[ tau[3]] = 0 ;
    MG_SET(xt[4].ori, tau[0]); MG_SET(xt[4].ori, tau[3]);

    pt[5]->v[tau[0]] = vx[taued[1]] ; pt[5]->v[tau[1]] = vx[taued[4]] ;
    xt[5].tag[taued[0]] = 0;  xt[5].tag[taued[2]] = 0;
    xt[5].tag[taued[3]] = 0;  xt[5].edg[taued[0]] = 0;
    xt[5].edg[taued[2]] = 0;  xt[5].edg[taued[3]] = 0;
    xt[5].ref [ tau[2]] = 0;  xt[5].ref [ tau[3]] = 0 ;
    xt[5].ftag[ tau[2]] = 0;  xt[5].ftag[ tau[3]] = 0 ;
    MG_SET(xt[5].ori, tau[2]); MG_SET(xt[5].ori, tau[3]);
  }

  /* Assignation of the xt fields to the appropriate tets */
  for (j=0; j<6;j++) {
    isxt[j] = 0;
  }

  for (i=0; i<4; i++) {
    for (j=0; j<6; j++) {
      if ( (xt[j]).ref[i] || xt[j].ftag[i] ) isxt[j] = 1;
    }
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for (i=1; i<6; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for (i=1; i<6; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
      }
      pt[0]->xt = 0;
    }
  }

  if ( (!metRidTyp) && met->m && met->size>1 ) {
    for (i=0; i<6; i++) {
      pt[i]->qual=MMG5_caltet33_ani(mesh,met,pt[i]);
    }
  }
  else {
    for (i=0; i<6; i++) {
      pt[i]->qual=MMG5_orcal(mesh,met,newtet[i]);
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 *  Simulate split of 4 edges in opposite configuration.
 *
 */
int MMG3D_split4op_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4];
  unsigned char       imin01,imin23;
  const unsigned char *taued;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* Set permutation of vertices : reference configuration 30 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];

  switch(pt->flag){
  case 45:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &MMG5_permedge[5][0];
    break;

  case 51:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;
  }

  imin01 = (pt->v[tau[0]] < pt->v[tau[1]]) ? tau[0] : tau[1];
  imin23 = (pt->v[tau[2]] < pt->v[tau[3]]) ? tau[2] : tau[3];

  /* Generic formulation for split of 4 edges, with no 3 edges lying on the same
   * face */
  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  if ( imin01 == tau[0] ) {
    pt0->v[tau[2]] = vx[taued[3]]; pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[1]] = vx[taued[4]]; pt0->v[tau[2]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[1]] = vx[taued[3]]; pt0->v[tau[2]] = vx[taued[1]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else {
    pt0->v[tau[2]] = vx[taued[1]]; pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[1]]; pt0->v[tau[2]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[2]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  if ( imin23 == tau[2] ) {
    pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[1]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[1]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[1]]; pt0->v[tau[1]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else {
    pt0->v[tau[0]] = vx[taued[1]]; pt0->v[tau[1]] = vx[taued[3]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[1]] = vx[taued[3]];
    pt0->v[tau[2]] = vx[taued[1]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[1]] = vx[taued[4]];
    pt0->v[tau[2]] = vx[taued[3]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split 4 edges in a configuration when no 3 edges lie on the same face
 *
 */
int MMG5_split4op(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra         pt[6];
  MMG5_xTetra         xt[6];
  MMG5_pxTetra        pxt0;
  int                 iel;
  int                 newtet[6];
  char                flg,firstxt,isxt[6],i,j,imin01,imin23;
  unsigned char       tau[4];
  const unsigned char *taued;

  pt[0]  = &mesh->tetra[k];
  flg = pt[0]->flag;
  pt[0]->flag  = 0;
  newtet[0]=k;

  /* Set permutation of vertices : reference configuration 30 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];

  switch(flg){
  case 45:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &MMG5_permedge[5][0];
    break;

  case 51:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;
  }

  imin01 = ((pt[0])->v[tau[0]] < (pt[0])->v[tau[1]]) ? tau[0] : tau[1];
  imin23 = ((pt[0])->v[tau[2]] < (pt[0])->v[tau[3]]) ? tau[2] : tau[3];

  /* create 5 new tetras */
  for (j=1; j<6; j++) {
    iel = MMG3D_newElt(mesh);
    if ( !iel ) {
      MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                  " a new element.\n",__func__);
                          MMG5_INCREASE_MEM_MESSAGE();
                          fprintf(stderr,"  Exit program.\n");
                          return 0);
      for ( i=0; i<j; i++)
        pt[i] = &mesh->tetra[newtet[i]];
    }
    pt[j] = &mesh->tetra[iel];
    pt[j] = memcpy(pt[j],pt[0],sizeof(MMG5_Tetra));
    newtet[j]=iel;
  }

  if ( (pt[0])->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    for (j=0; j<6; j++) {
      memcpy(&xt[j],pxt0, sizeof(MMG5_xTetra));
    }
  }
  else {
    pxt0 = 0;
    for (j=0; j<6; j++) {
      memset(&xt[j],0, sizeof(MMG5_xTetra));
    }
  }

  /* Generic formulation for split of 4 edges, with no 3 edges lying on the same face */
  if ( imin01 == tau[0] ) {
    pt[0]->v[tau[2]] = vx[taued[3]] ; pt[0]->v[tau[3]] = vx[taued[4]];
    xt[0].tag[taued[1]] = 0;  xt[0].tag[taued[5]] = 0;
    xt[0].tag[taued[2]] = 0;  xt[0].edg[taued[1]] = 0;
    xt[0].edg[taued[5]] = 0;  xt[0].edg[taued[2]] = 0;
    xt[0].ref [ tau[1]] = 0;
    xt[0].ftag[ tau[1]] = 0;
    MG_SET(xt[0].ori, tau[1]);

    pt[1]->v[tau[1]] = vx[taued[4]] ; pt[1]->v[tau[2]] = vx[taued[3]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].tag[taued[0]] = 0;  xt[1].tag[taued[1]] = 0;
    xt[1].tag[taued[3]] = 0;  xt[1].tag[taued[4]] = 0;
    xt[1].tag[taued[5]] = 0;  xt[1].edg[taued[0]] = 0;
    xt[1].edg[taued[1]] = 0;  xt[1].edg[taued[3]] = 0;
    xt[1].edg[taued[4]] = 0;  xt[1].edg[taued[5]] = 0;
    xt[1].ref [ tau[0]] = 0;  xt[1].ref [ tau[1]] = 0;  xt[1].ref [tau[3]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[1]] = 0;  xt[1].ftag[tau[3]] = 0;
    MG_SET(xt[1].ori, tau[0]);  MG_SET(xt[1].ori, tau[1]);  MG_SET(xt[1].ori, tau[3]);

    pt[2]->v[tau[1]] = vx[taued[3]] ; pt[2]->v[tau[2]] = vx[taued[1]] ; pt[2]->v[tau[3]] = vx[taued[2]];
    xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[3]] = 0;
    xt[2].tag[taued[4]] = 0;  xt[2].tag[taued[5]] = 0;
    xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[3]] = 0;
    xt[2].edg[taued[4]] = 0;  xt[2].edg[taued[5]] = 0;
    xt[2].ref [ tau[0]] = 0;  xt[2].ref [ tau[2]] = 0;
    xt[2].ftag[ tau[0]] = 0;  xt[2].ftag[ tau[2]] = 0;
    MG_SET(xt[2].ori, tau[0]);  MG_SET(xt[2].ori, tau[2]);
  }
  else {
    pt[0]->v[tau[2]] = vx[taued[1]] ; pt[0]->v[tau[3]] = vx[taued[2]];
    xt[0].tag[taued[3]] = 0;  xt[0].tag[taued[4]] = 0;
    xt[0].tag[taued[5]] = 0;  xt[0].edg[taued[3]] = 0;
    xt[0].edg[taued[4]] = 0;  xt[0].edg[taued[5]] = 0;
    xt[0].ref [ tau[0]] = 0;
    xt[0].ftag[ tau[0]] = 0;
    MG_SET(xt[0].ori, tau[0]);

    pt[1]->v[tau[0]] = vx[taued[1]] ; pt[1]->v[tau[2]] = vx[taued[3]] ; pt[1]->v[tau[3]] = vx[taued[2]];
    xt[1].tag[taued[0]] = 0;  xt[1].tag[taued[1]] = 0;
    xt[1].tag[taued[2]] = 0;  xt[1].tag[taued[4]] = 0;
    xt[1].tag[taued[5]] = 0;  xt[1].edg[taued[0]] = 0;
    xt[1].edg[taued[1]] = 0;  xt[1].edg[taued[2]] = 0;
    xt[1].edg[taued[4]] = 0;  xt[1].edg[taued[5]] = 0;
    xt[1].ref [ tau[0]] = 0;  xt[1].ref [ tau[1]] = 0;  xt[1].ref [tau[2]] = 0;
    xt[1].ftag[ tau[0]] = 0;  xt[1].ftag[ tau[1]] = 0;  xt[1].ftag[tau[2]] = 0;
    MG_SET(xt[1].ori, tau[0]);  MG_SET(xt[1].ori, tau[1]);  MG_SET(xt[1].ori, tau[2]);

    pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[2]] = vx[taued[3]] ; pt[2]->v[tau[3]] = vx[taued[4]];
    xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[1]] = 0;
    xt[2].tag[taued[2]] = 0;  xt[2].tag[taued[5]] = 0;
    xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[1]] = 0;
    xt[2].edg[taued[2]] = 0;  xt[2].edg[taued[5]] = 0;
    xt[2].ref [ tau[1]] = 0;  xt[2].ref [ tau[3]] = 0;
    xt[2].ftag[ tau[1]] = 0;  xt[2].ftag[ tau[3]] = 0;
    MG_SET(xt[2].ori, tau[1]);  MG_SET(xt[2].ori, tau[3]);
  }

  if ( imin23 == tau[2] ) {
    pt[3]->v[tau[0]] = vx[taued[2]] ; pt[3]->v[tau[1]] = vx[taued[4]];
    xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[1]] = 0;
    xt[3].tag[taued[3]] = 0;  xt[3].edg[taued[0]] = 0;
    xt[3].edg[taued[1]] = 0;  xt[3].edg[taued[3]] = 0;
    xt[3].ref [ tau[3]] = 0;
    xt[3].ftag[ tau[3]] = 0;
    MG_SET(xt[3].ori, tau[3]);

    pt[4]->v[tau[0]] = vx[taued[2]] ; pt[4]->v[tau[1]] = vx[taued[3]] ; pt[4]->v[tau[3]] = vx[taued[4]];
    xt[4].tag[taued[0]] = 0;  xt[4].tag[taued[1]] = 0;
    xt[4].tag[taued[2]] = 0;  xt[4].tag[taued[4]] = 0;
    xt[4].tag[taued[5]] = 0;  xt[4].edg[taued[0]] = 0;
    xt[4].edg[taued[1]] = 0;  xt[4].edg[taued[2]] = 0;
    xt[4].edg[taued[4]] = 0;  xt[4].edg[taued[5]] = 0;
    xt[4].ref [ tau[1]] = 0;  xt[4].ref [ tau[2]] = 0;  xt[4].ref [tau[3]] = 0;
    xt[4].ftag[ tau[1]] = 0;  xt[4].ftag[ tau[2]] = 0;  xt[4].ftag[tau[3]] = 0;
    MG_SET(xt[4].ori, tau[1]);  MG_SET(xt[4].ori, tau[2]);  MG_SET(xt[4].ori, tau[3]);

    pt[5]->v[tau[0]] = vx[taued[1]] ; pt[5]->v[tau[1]] = vx[taued[3]] ; pt[5]->v[tau[3]] = vx[taued[2]];
    xt[5].tag[taued[0]] = 0;  xt[5].tag[taued[2]] = 0;
    xt[5].tag[taued[4]] = 0;  xt[5].tag[taued[5]] = 0;
    xt[5].edg[taued[0]] = 0;  xt[5].edg[taued[2]] = 0;
    xt[5].edg[taued[4]] = 0;  xt[5].edg[taued[5]] = 0;
    xt[5].ref [ tau[0]] = 0;  xt[5].ref [ tau[2]] = 0;
    xt[5].ftag[ tau[0]] = 0;  xt[5].ftag[ tau[2]] = 0;
    MG_SET(xt[5].ori, tau[0]);  MG_SET(xt[5].ori, tau[2]);
  }
  else {
    pt[3]->v[tau[0]] = vx[taued[1]] ; pt[3]->v[tau[1]] = vx[taued[3]];
    xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[2]] = 0;
    xt[3].tag[taued[4]] = 0;  xt[3].edg[taued[0]] = 0;
    xt[3].edg[taued[2]] = 0;  xt[3].edg[taued[4]] = 0;
    xt[3].ref [ tau[2]] = 0;
    xt[3].ftag[ tau[2]] = 0;
    MG_SET(xt[3].ori, tau[2]);

    pt[4]->v[tau[0]] = vx[taued[2]] ; pt[4]->v[tau[1]] = vx[taued[3]] ; pt[4]->v[tau[2]] = vx[taued[1]];
    xt[4].tag[taued[0]] = 0;  xt[4].tag[taued[1]] = 0;
    xt[4].tag[taued[3]] = 0;  xt[4].tag[taued[4]] = 0;
    xt[4].tag[taued[5]] = 0;  xt[4].edg[taued[0]] = 0;
    xt[4].edg[taued[1]] = 0;  xt[4].edg[taued[3]] = 0;
    xt[4].edg[taued[4]] = 0;  xt[4].edg[taued[5]] = 0;
    xt[4].ref [ tau[0]] = 0;  xt[4].ref [ tau[2]] = 0;  xt[4].ref [tau[3]] = 0;
    xt[4].ftag[ tau[0]] = 0;  xt[4].ftag[ tau[2]] = 0;  xt[4].ftag[tau[3]] = 0;
    MG_SET(xt[4].ori, tau[0]);  MG_SET(xt[4].ori, tau[2]);  MG_SET(xt[4].ori, tau[3]);

    pt[5]->v[tau[0]] = vx[taued[2]] ; pt[5]->v[tau[1]] = vx[taued[4]] ; pt[5]->v[tau[2]] = vx[taued[3]];
    xt[5].tag[taued[0]] = 0;  xt[5].tag[taued[1]] = 0;
    xt[5].tag[taued[3]] = 0;  xt[5].tag[taued[5]] = 0;
    xt[5].edg[taued[0]] = 0;  xt[5].edg[taued[1]] = 0;
    xt[5].edg[taued[3]] = 0;  xt[5].edg[taued[5]] = 0;
    xt[5].ref [ tau[1]] = 0;  xt[5].ref [ tau[3]] = 0;
    xt[5].ftag[ tau[1]] = 0;  xt[5].ftag[ tau[3]] = 0;
    MG_SET(xt[5].ori, tau[1]); MG_SET(xt[5].ori, tau[3]);
  }

  /* Assignation of the xt fields to the appropriate tets */
  for (j=0; j<6; j++) {
    isxt[j] = 0;
  }

  for (i=0; i<4; i++) {
    for(j=0;j<6;j++ ) {
      if ( (xt[j]).ref[i] || xt[j].ftag[i] ) isxt[j] = 1;
    }
  }

  // In this case, at least one of the 4 created tets must have a special field
  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for (i=1; i<6; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = 0;

      for (i=1; i<6; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[ pt[i]->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&(xt[i]),sizeof(MMG5_xTetra));
          }
        }
      }
      pt[0]->xt = 0;

    }
  }
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    for (i=0; i<6; i++) {
      pt[i]->qual=MMG5_caltet33_ani(mesh,met,pt[i]);
    }
  }
  else {
    for (i=0; i<6; i++) {
      pt[i]->qual=MMG5_orcal(mesh,met,newtet[i]);
    }
  }
  return 1;
}

/**
 * \param pt initial tetra
 * \param vx index of points to insert along edges
 * \param tau vertices permutation
 * \param taued edges permutation
 * \param imin minimal index of vertices \a tau[0] and \a tau[1]
 *
 * Set permutation of vertices for the split of 5 edges.
 * Reference configuration is 62.
 *
 */
static inline
void MMG3D_configSplit5(MMG5_pTetra pt,int vx[6],unsigned char tau[4],
                         const unsigned char **taued,unsigned char *imin) {

  /* set permutation of vertices and edges ; reference configuration : 62 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  (*taued)  = &MMG5_permedge[0][0];

  switch(pt->flag) {
  case 61:
    tau[0] = 2 ; tau[1] = 0 ; tau[2] = 1 ; tau[3] = 3;
    (*taued)  = &MMG5_permedge[6][0];
    break;

  case 59:
    tau[0] = 0 ; tau[1] = 3 ; tau[2] = 1 ; tau[3] = 2;
    (*taued)  = &MMG5_permedge[2][0];
    break;

  case 55:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    (*taued)  = &MMG5_permedge[4][0];
    break;

  case 47:
    tau[0] = 3 ; tau[1] = 1 ; tau[2] = 0 ; tau[3] = 2;
    (*taued)  = &MMG5_permedge[10][0];
    break;

  case 31:
    tau[0] = 3 ; tau[1] = 2 ; tau[2] = 1 ; tau[3] = 0;
    (*taued)  = &MMG5_permedge[11][0];
    break;
  }

  /* Generic formulation of split of 5 edges */
  (*imin) = (pt->v[tau[0]] < pt->v[tau[1]]) ? tau[0] : tau[1];
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 *  Simulate split of 5 edges.
 *
 */
int MMG3D_split5_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;
  unsigned char       tau[4];
  unsigned char       imin;
  const unsigned char *taued=NULL;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* Set permutation of vertices : reference configuration : 62 */
  MMG3D_configSplit5(pt,vx,tau,&taued,&imin);

  /* Generic formulation of split of 5 edges */
  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[1]] = vx[taued[4]];
  pt0->v[tau[2]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[1]]; pt0->v[tau[1]] = vx[taued[3]];
  pt0->v[tau[3]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[1]] = vx[taued[4]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[1]] = vx[taued[3]];
  pt0->v[tau[2]] = vx[taued[1]]; pt0->v[tau[3]] = vx[taued[5]];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  memcpy(pt0,pt,sizeof(MMG5_Tetra));
  if ( imin == tau[0] ) {
    pt0->v[tau[2]] = vx[taued[3]]; pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[1]] = vx[taued[4]]; pt0->v[tau[2]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[1]] = vx[taued[3]]; pt0->v[tau[2]] = vx[taued[1]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }
  else {
    pt0->v[tau[2]] = vx[taued[1]]; pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[2]]; pt0->v[tau[2]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[4]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;

    memcpy(pt0,pt,sizeof(MMG5_Tetra));
    pt0->v[tau[0]] = vx[taued[1]]; pt0->v[tau[2]] = vx[taued[3]];
    pt0->v[tau[3]] = vx[taued[2]];
    vnew = MMG5_orvol(mesh->point,pt0->v);
    if ( vnew < MMG5_EPSOK )  return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * Split 5 edges
 *
 */
int MMG5_split5(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra         pt[7];
  MMG5_xTetra         xt[7];
  MMG5_pxTetra        pxt0;
  int                 iel,i,j;
  int                 newtet[7];
  char                firstxt,isxt[7];
  unsigned char       tau[4],imin;
  const unsigned char *taued=NULL;

  pt[0]  = &mesh->tetra[k];
  newtet[0]=k;

  /* set permutation of vertices and edges ; reference configuration : 62 */
  MMG3D_configSplit5(pt[0],vx,tau,&taued,&imin);
  pt[0]->flag  = 0;

  /* create 6 new tetras */
  for (i=1; i<7; i++) {
    iel = MMG3D_newElt(mesh);
    if ( !iel ) {
      MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                  " a new element.\n",__func__);
                          MMG5_INCREASE_MEM_MESSAGE();
                          fprintf(stderr,"  Exit program.\n");
                          return 0);
      for ( j=0; j<i; j++)
        pt[j] = &mesh->tetra[newtet[j]];
    }
    pt[i] = &mesh->tetra[iel];
    pt[i] = memcpy(pt[i],pt[0],sizeof(MMG5_Tetra));
    newtet[i]=iel;
  }

  if ( pt[0]->xt ) {
    pxt0 = &mesh->xtetra[(pt[0])->xt];
    for (i=0; i<7; i++) {
      memcpy(&xt[i],pxt0, sizeof(MMG5_xTetra));
    }
  }
  else {
    pxt0 = 0;
    for (i=0; i<7; i++) {
      memset(&xt[i],0, sizeof(MMG5_xTetra));
    }
  }

  /* Generic formulation of split of 5 edges */
  pt[0]->v[tau[0]] = vx[taued[2]] ;   pt[0]->v[tau[1]] = vx[taued[4]] ;   pt[0]->v[tau[2]] = vx[taued[5]];
  xt[0].tag[taued[0]] = 0;  xt[0].tag[taued[1]] = 0;
  xt[0].tag[taued[3]] = 0;  xt[0].edg[taued[0]] = 0;
  xt[0].edg[taued[1]] = 0;  xt[0].edg[taued[3]] = 0;
  xt[0].ref [ tau[3]] = 0;
  xt[0].ftag[ tau[3]] = 0;
  MG_SET(xt[0].ori, tau[3]);

  pt[1]->v[tau[0]] = vx[taued[1]] ; pt[1]->v[tau[1]] = vx[taued[3]] ; pt[1]->v[tau[3]] = vx[taued[5]];
  xt[1].tag[taued[0]] = 0;  xt[1].tag[taued[2]] = 0;
  xt[1].tag[taued[4]] = 0;  xt[1].edg[taued[0]] = 0;
  xt[1].edg[taued[2]] = 0;  xt[1].edg[taued[4]] = 0;
  xt[1].ref [ tau[2]] = 0;
  xt[1].ftag[ tau[2]] = 0;
  MG_SET(xt[1].ori, tau[2]);

  pt[2]->v[tau[0]] = vx[taued[2]] ; pt[2]->v[tau[1]] = vx[taued[4]];
  pt[2]->v[tau[2]] = vx[taued[3]] ; pt[2]->v[tau[3]] = vx[taued[5]];
  xt[2].tag[taued[0]] = 0;  xt[2].tag[taued[1]] = 0;
  xt[2].tag[taued[2]] = 0;  xt[2].tag[taued[3]] = 0;
  xt[2].tag[taued[4]] = 0;  xt[2].tag[taued[5]] = 0;
  xt[2].edg[taued[0]] = 0;  xt[2].edg[taued[1]] = 0;
  xt[2].edg[taued[2]] = 0;  xt[2].edg[taued[3]] = 0;
  xt[2].edg[taued[4]] = 0;  xt[2].edg[taued[5]] = 0;
  xt[2].ref [tau[1]] = 0 ;  xt[2].ref [ tau[2]] = 0;  xt[2].ref [tau[3]] = 0 ;
  xt[2].ftag[tau[1]] = 0 ;  xt[2].ftag[ tau[2]] = 0;  xt[2].ftag[tau[3]] = 0 ;
  MG_SET(xt[2].ori, tau[1]);  MG_SET(xt[2].ori, tau[2]);  MG_SET(xt[2].ori, tau[3]);

  pt[3]->v[tau[0]] = vx[taued[2]] ; pt[3]->v[tau[1]] = vx[taued[3]];
  pt[3]->v[tau[2]] = vx[taued[1]] ; pt[3]->v[tau[3]] = vx[taued[5]];
  xt[3].tag[taued[0]] = 0;  xt[3].tag[taued[1]] = 0;
  xt[3].tag[taued[2]] = 0;  xt[3].tag[taued[3]] = 0;
  xt[3].tag[taued[4]] = 0;  xt[3].tag[taued[5]] = 0;
  xt[3].edg[taued[0]] = 0;  xt[3].edg[taued[1]] = 0;
  xt[3].edg[taued[2]] = 0;  xt[3].edg[taued[3]] = 0;
  xt[3].edg[taued[4]] = 0;  xt[3].edg[taued[5]] = 0;
  xt[3].ref [ tau[0]] = 0;  xt[3].ref [ tau[2]] = 0;  xt[3].ref [tau[3]] = 0 ;
  xt[3].ftag[ tau[0]] = 0;  xt[3].ftag[ tau[2]] = 0;  xt[3].ftag[tau[3]] = 0 ;
  MG_SET(xt[3].ori, tau[0]);  MG_SET(xt[3].ori, tau[2]);  MG_SET(xt[3].ori, tau[3]);

  if ( imin == tau[0] ) {
    pt[4]->v[tau[2]] = vx[taued[3]] ; pt[4]->v[tau[3]] = vx[taued[4]];
    xt[4].tag[taued[1]] = 0;  xt[4].tag[taued[2]] = 0;
    xt[4].tag[taued[5]] = 0;  xt[4].edg[taued[1]] = 0;
    xt[4].edg[taued[2]] = 0;  xt[4].edg[taued[5]] = 0;
    xt[4].ref [ tau[1]] = 0;
    xt[4].ftag[ tau[1]] = 0;
    MG_SET(xt[4].ori, tau[1]);

    pt[5]->v[tau[1]] = vx[taued[4]] ; pt[5]->v[tau[2]] = vx[taued[3]]; pt[5]->v[tau[3]] = vx[taued[2]];
    xt[5].tag[taued[0]] = 0;
    xt[5].tag[taued[1]] = 0;  xt[5].tag[taued[3]] = 0;
    xt[5].tag[taued[4]] = 0;  xt[5].tag[taued[5]] = 0;
    xt[5].edg[taued[0]] = 0;
    xt[5].edg[taued[1]] = 0;  xt[5].edg[taued[3]] = 0;
    xt[5].edg[taued[4]] = 0;  xt[5].edg[taued[5]] = 0;
    xt[5].ref [ tau[0]] = 0;  xt[5].ref [ tau[1]] = 0; xt[5].ref [tau[3]] = 0 ;
    xt[5].ftag[ tau[0]] = 0;  xt[5].ftag[ tau[1]] = 0; xt[5].ftag[tau[3]] = 0 ;
    MG_SET(xt[5].ori, tau[0]); MG_SET(xt[5].ori, tau[1]); MG_SET(xt[5].ori, tau[3]);

    pt[6]->v[tau[1]] = vx[taued[3]] ; pt[6]->v[tau[2]] = vx[taued[1]]; pt[6]->v[tau[3]] = vx[taued[2]];
    xt[6].tag[taued[0]] = 0;
    xt[6].tag[taued[3]] = 0;  xt[6].tag[taued[4]] = 0;
    xt[6].tag[taued[5]] = 0;  xt[6].edg[taued[0]] = 0;
    xt[6].edg[taued[3]] = 0;
    xt[6].edg[taued[4]] = 0;  xt[6].edg[taued[5]] = 0;
    xt[6].ref [ tau[0]] = 0;  xt[6].ref [ tau[2]] = 0;
    xt[6].ftag[ tau[0]] = 0;  xt[6].ftag[ tau[2]] = 0;
    MG_SET(xt[6].ori, tau[0]); MG_SET(xt[6].ori, tau[2]);

  }
  else {
    pt[4]->v[tau[2]] = vx[taued[1]] ; pt[4]->v[tau[3]] = vx[taued[2]];
    xt[4].tag[taued[3]] = 0;  xt[4].tag[taued[4]] = 0;
    xt[4].tag[taued[5]] = 0;  xt[4].edg[taued[3]] = 0;
    xt[4].edg[taued[4]] = 0;  xt[4].edg[taued[5]] = 0;
    xt[4].ref [ tau[0]] = 0;
    xt[4].ftag[ tau[0]] = 0;
    MG_SET(xt[4].ori, tau[0]);

    pt[5]->v[tau[0]] = vx[taued[2]] ; pt[5]->v[tau[2]] = vx[taued[3]]; pt[5]->v[tau[3]] = vx[taued[4]];
    xt[5].tag[taued[0]] = 0;  xt[5].tag[taued[1]] = 0;
    xt[5].tag[taued[2]] = 0;  xt[5].tag[taued[5]] = 0;
    xt[5].edg[taued[0]] = 0;  xt[5].edg[taued[1]] = 0;
    xt[5].edg[taued[2]] = 0;  xt[5].edg[taued[5]] = 0;
    xt[5].ref [ tau[1]] = 0; xt[5].ref [ tau[3]] = 0;
    xt[5].ftag[ tau[1]] = 0; xt[5].ftag[ tau[3]] = 0;
    MG_SET(xt[5].ori, tau[1]); MG_SET(xt[5].ori, tau[3]);

    pt[6]->v[tau[0]] = vx[taued[1]] ; pt[6]->v[tau[2]] = vx[taued[3]]; pt[6]->v[tau[3]] = vx[taued[2]];
    xt[6].tag[taued[0]] = 0;  xt[6].tag[taued[1]] = 0;
    xt[6].tag[taued[2]] = 0;  xt[6].tag[taued[4]] = 0;
    xt[6].tag[taued[5]] = 0;  xt[6].edg[taued[0]] = 0;
    xt[6].edg[taued[1]] = 0;  xt[6].edg[taued[2]] = 0;
    xt[6].edg[taued[4]] = 0;  xt[6].edg[taued[5]] = 0;
    xt[6].ref [ tau[0]] = 0;  xt[6].ref [ tau[1]] = 0; xt[6].ref [tau[2]] = 0 ;
    xt[6].ftag[ tau[0]] = 0;  xt[6].ftag[ tau[1]] = 0; xt[6].ftag[tau[2]] = 0 ;
    MG_SET(xt[6].ori, tau[0]); MG_SET(xt[6].ori, tau[1]); MG_SET(xt[6].ori, tau[2]);
  }

  /* Assignation of the xt fields to the appropriate tets */
  for (j=0; j<7; j++) {
    isxt[j] = 0;
  }

  for (i=0; i<4; i++) {
    for (j=0; j<7; j++) {
      if ( (xt[j]).ref[i] || xt[j].ftag[i] ) isxt[j] = 1;
    }
  }

  if ( pt[0]->xt ) {
    if ( isxt[0] ) {
      memcpy(pxt0,&xt[0],sizeof(MMG5_xTetra));
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = pt[6]->xt = 0;

      for (i=1; i<7; i++) {
        if ( isxt[i] ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return 0);
          }
          pt[i]->xt = mesh->xt;
          pxt0 = &mesh->xtetra[mesh->xt];
          memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
        }
      }
    }
    else {
      firstxt = 1;
      pt[1]->xt = pt[2]->xt = pt[3]->xt = pt[4]->xt = pt[5]->xt = pt[6]->xt = 0;

      for (i=1; i<7; i++) {
        if ( isxt[i] ) {
          if ( firstxt ) {
            firstxt = 0;
            pt[i]->xt = pt[0]->xt;
            pxt0 = &mesh->xtetra[(pt[i])->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
          else {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              /* realloc of xtetras table */
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");
                                 return 0);
            }
            pt[i]->xt = mesh->xt;
            pxt0 = &mesh->xtetra[mesh->xt];
            memcpy(pxt0,&xt[i],sizeof(MMG5_xTetra));
          }
        }
      }
      pt[0]->xt = 0;

    }
  }
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    for (i=0; i<7; i++) {
      pt[i]->qual=MMG5_caltet33_ani(mesh,met,pt[i]);
    }
  }
  else {
    for (i=0; i<7; i++) {
      pt[i]->qual=MMG5_orcal(mesh,met,newtet[i]);
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 *
 * \return 0 if the split fail, 1 otherwise
 *
 *  Simulate split of 6 edges.
 *
 */
int MMG3D_split6_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]) {
  MMG5_pTetra         pt,pt0;
  double              vold,vnew;

  pt  = &mesh->tetra[k];
  pt0 = &mesh->tetra[0];
  vold = MMG5_orvol(mesh->point,pt->v);

  if ( vold < MMG5_EPSOK ) return 0;

  /* Modify first tetra */
  pt0->v[1] = vx[0]; pt0->v[2] = vx[1]; pt0->v[3] = vx[2];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Modify second tetra */
  pt0->v[0] = vx[0]; pt0->v[2] = vx[3]; pt0->v[3] = vx[4];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Modify 3rd tetra */
  pt0->v[0] = vx[1]; pt0->v[1] = vx[3]; pt0->v[3] = vx[5];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Modify 4th tetra */
  pt0->v[0] = vx[2]; pt0->v[1] = vx[4]; pt0->v[2] = vx[5];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Modify 5th tetra */
  pt0->v[0] = vx[0]; pt0->v[1] = vx[3]; pt0->v[2] = vx[1];
  pt0->v[3] = vx[2];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Modify 6th tetra */
  pt0->v[0] = vx[2]; pt0->v[1] = vx[0]; pt0->v[2] = vx[3];
  pt0->v[3] = vx[4];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Modify 7th tetra */
  pt0->v[0] = vx[2]; pt0->v[1] = vx[3]; pt0->v[2] = vx[1];
  pt0->v[3] = vx[5];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Modify last tetra */
  pt0->v[0] = vx[2]; pt0->v[1] = vx[3]; pt0->v[2] = vx[5];
  pt0->v[3] = vx[4];
  vnew = MMG5_orvol(mesh->point,pt0->v);
  if ( vnew < MMG5_EPSOK )  return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param metRidTyp metric storage (classic or special)
 *
 * \return 0 if fail, 1 otherwise
 *
 * split all faces (6 edges)
 *
 */
int MMG5_split6(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp) {
  MMG5_pTetra    pt[8];
  MMG5_xTetra    xt0,xt;
  MMG5_pxTetra   pxt;
  int       i,j,iel,nxt0;
  int       newtet[8];
  char      isxt0,isxt;

  pt[0]  = &mesh->tetra[k];
  pt[0]->flag  = 0;
  newtet[0]=k;

  nxt0 = pt[0]->xt;
  pxt = &mesh->xtetra[nxt0];
  memcpy(&xt0,pxt,sizeof(MMG5_xTetra));

  /* create 7 new tetras */
  for (i=1; i<8; i++) {
    iel = MMG3D_newElt(mesh);
    if ( !iel ) {
      MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                  " a new element.\n",__func__);
                          MMG5_INCREASE_MEM_MESSAGE();
                          fprintf(stderr,"  Exit program.\n");
                          return 0);
      for ( j=0; j<i; j++ )
        pt[j] = &mesh->tetra[newtet[j]];
    }
    pt[i] = &mesh->tetra[iel];
    pt[i] = memcpy(pt[i],pt[0],sizeof(MMG5_Tetra));
    newtet[i]=iel;
  }

  /* Modify first tetra */
  pt[0]->v[1] = vx[0] ; pt[0]->v[2] = vx[1]; pt[0]->v[3] = vx[2];
  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[3] = 0;  xt.tag[4] = 0;
    xt.tag[5] = 0;  xt.edg[3] = 0;
    xt.edg[4] = 0;  xt.edg[5] = 0;
    xt.ref[0] = 0;  xt.ftag[0] = 0; MG_SET(xt.ori, 0);
    isxt0 = 0;
    for(i=0;i<4;i++ ) {
      if ( (xt.ref[i]) || xt.ftag[i] ) isxt0 = 1;
    }

    if ( isxt0 ) {
      memcpy(pxt,&xt,sizeof(MMG5_xTetra));
    }
    else {
      pt[0]->xt = 0;
    }
  }

  /* Modify second tetra */
  pt[1]->v[0] = vx[0] ; pt[1]->v[2] = vx[3]; pt[1]->v[3] = vx[4];

  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[1] = 0;  xt.tag[2] = 0;
    xt.tag[5] = 0;  xt.edg[1] = 0;
    xt.edg[2] = 0;  xt.edg[5] = 0;
    xt.ref[1] = 0;  xt.ftag[1] = 0; MG_SET(xt.ori, 1);

    isxt = 0;

    for (i=0; i<4; i++) {
      if ( (xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[1]->xt = 0;
    if ( isxt ) {
      if ( !isxt0 ) {
        isxt0 = 1;
        pt[1]->xt = nxt0;
        pxt = &mesh->xtetra[pt[1]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return 0);
        }
        pt[1]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[1]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
    }
  }

  /* Modify 3rd tetra */
  pt[2]->v[0] = vx[1] ; pt[2]->v[1] = vx[3]; pt[2]->v[3] = vx[5];

  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[0] = 0;  xt.tag[2] = 0;
    xt.tag[4] = 0;  xt.edg[0] = 0;
    xt.edg[2] = 0;  xt.edg[4] = 0;
    xt.ref[2] = 0;  xt.ftag[2] = 0;  MG_SET(xt.ori, 2);
    isxt = 0;

    for (i=0; i<4;i++) {
      if ( (xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[2]->xt = 0;
    if ( isxt ) {
      if ( !isxt0 ) {
        isxt0 = 1;
        pt[2]->xt = nxt0;
        pxt = &mesh->xtetra[pt[2]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return 0);
        }
        pt[2]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[2]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
    }
  }

  /* Modify 4th tetra */
  pt[3]->v[0] = vx[2] ; pt[3]->v[1] = vx[4]; pt[3]->v[2] = vx[5];

  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[0] = 0;  xt.tag[1] = 0;
    xt.tag[3] = 0;  xt.edg[0] = 0;
    xt.edg[1] = 0;  xt.edg[3] = 0;
    xt.ref[3] = 0;  xt.ftag[3] = 0;  MG_SET(xt.ori, 3);

    isxt = 0;

    for (i=0; i<4; i++) {
      if ( (xt.ref[i]) || xt.ftag[i]) isxt = 1;
    }

    pt[3]->xt = 0;
    if ( isxt ) {
      if ( !isxt0 ) {
        isxt0 = 1;
        pt[3]->xt = nxt0;
        pxt = &mesh->xtetra[pt[3]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return 0);
        }
        pt[3]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[3]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
    }
  }

  /* Modify 5th tetra */
  pt[4]->v[0] = vx[0] ; pt[4]->v[1] = vx[3]; pt[4]->v[2] = vx[1] ; pt[4]->v[3] = vx[2];

  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[0] = 0;  xt.tag[1] = 0;
    xt.tag[2] = 0;  xt.tag[3] = 0;
    xt.edg[0] = 0;  xt.edg[1] = 0;
    xt.edg[2] = 0;  xt.edg[3] = 0;
    xt.tag[4] = 0;  xt.edg[4] = 0;
    xt.tag[5] = 0;  xt.edg[5] = 0;
    xt.ref [0] = 0 ; xt.ref [1] = 0 ; xt.ref [2] = 0;
    xt.ftag[0] = 0 ; xt.ftag[1] = 0 ; xt.ftag[2] = 0;
    MG_SET(xt.ori, 0); MG_SET(xt.ori, 1); MG_SET(xt.ori, 2);

    isxt = 0;

    if ( (xt.ref[3]) || xt.ftag[3]) isxt = 1;

    pt[4]->xt = 0;
    if ( isxt ) {
      if ( !isxt0 ) {
        isxt0 = 1;
        pt[4]->xt = nxt0;
        pxt = &mesh->xtetra[(pt[4])->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return 0);
        }
        pt[4]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[4]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
    }
  }

  /* Modify 6th tetra */
  pt[5]->v[0] = vx[2] ; pt[5]->v[1] = vx[0]; pt[5]->v[2] = vx[3] ; pt[5]->v[3] = vx[4];

  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[0] = 0;  xt.tag[1] = 0;
    xt.tag[2] = 0;  xt.tag[3] = 0;
    xt.tag[4] = 0;  xt.tag[5] = 0;
    xt.edg[0] = 0;  xt.edg[1] = 0;
    xt.edg[2] = 0;  xt.edg[3] = 0;
    xt.edg[4] = 0;  xt.edg[5] = 0;
    xt.ref [0] = 0 ; xt.ref [1] = 0 ; xt.ref [3] = 0;
    xt.ftag[0] = 0 ; xt.ftag[1] = 0 ; xt.ftag[3] = 0;
    MG_SET(xt.ori, 0); MG_SET(xt.ori, 1); MG_SET(xt.ori, 3);

    isxt = 0;

    if ( (xt.ref[2]) || xt.ftag[2]) isxt = 1;

    pt[5]->xt = 0;
    if ( isxt ) {
      if ( !isxt0 ) {
        isxt0 = 1;
        pt[5]->xt = nxt0;
        pxt = &mesh->xtetra[pt[5]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return 0);
        }
        pt[5]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[5]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
    }
  }

  /* Modify 7th tetra */
  pt[6]->v[0] = vx[2] ; pt[6]->v[1] = vx[3]; pt[6]->v[2] = vx[1] ; pt[6]->v[3] = vx[5];

  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[0] = 0;  xt.edg[0] = 0;
    xt.tag[1] = 0;  xt.tag[2] = 0;
    xt.tag[3] = 0;  xt.tag[4] = 0;
    xt.edg[1] = 0;  xt.edg[2] = 0;
    xt.edg[3] = 0;  xt.edg[4] = 0;
    xt.tag[5] = 0;  xt.edg[5] = 0;
    xt.ref [0] = 0 ; xt.ref [2] = 0 ; xt.ref [3] = 0;
    xt.ftag[0] = 0 ; xt.ftag[2] = 0 ; xt.ftag[3] = 0;
    MG_SET(xt.ori, 0); MG_SET(xt.ori, 2); MG_SET(xt.ori, 3);

    isxt = 0;

    if ( (xt.ref[1]) || xt.ftag[1]) isxt = 1;

    pt[6]->xt = 0;
    if ( isxt ) {
      if ( !isxt0 ) {
        isxt0 = 1;
        pt[6]->xt = nxt0;
        pxt = &mesh->xtetra[pt[6]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return 0);
        }
        pt[6]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[6]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
    }
  }

  /* Modify last tetra */
  pt[7]->v[0] = vx[2] ; pt[7]->v[1] = vx[3]; pt[7]->v[2] = vx[5] ; pt[7]->v[3] = vx[4];

  if ( nxt0 ) {
    memcpy(&xt,&xt0,sizeof(MMG5_xTetra));
    xt.tag[0] = 0;  xt.tag[1] = 0;
    xt.tag[2] = 0;  xt.tag[3] = 0;
    xt.tag[4] = 0;  xt.tag[5] = 0;
    xt.edg[0] = 0;  xt.edg[1] = 0;
    xt.edg[2] = 0;  xt.edg[3] = 0;
    xt.edg[4] = 0;  xt.edg[5] = 0;
    xt.ref [1] = 0 ; xt.ref [2] = 0 ; xt.ref [3] = 0;
    xt.ftag[1] = 0 ; xt.ftag[2] = 0 ; xt.ftag[3] = 0;
    MG_SET(xt.ori, 1); MG_SET(xt.ori, 2); MG_SET(xt.ori, 3);

    isxt = 0;

    if ( (xt.ref[0]) || xt.ftag[0]) isxt = 1;

    pt[7]->xt = 0;
    if ( isxt ) {
      if ( !isxt0 ) {
        pt[7]->xt = nxt0;
        pxt = &mesh->xtetra[pt[7]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
      else {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          /* realloc of xtetras table */
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh->xt--;
                             fprintf(stderr,"  Exit program.\n");
                             return 0);
        }
        pt[7]->xt = mesh->xt;
        pxt = &mesh->xtetra[pt[7]->xt];
        memcpy(pxt,&xt,sizeof(MMG5_xTetra));
      }
    }
  }
  if ( (!metRidTyp) && met->m && met->size>1 ) {
    for (i=0; i<8; i++) {
      pt[i]->qual=MMG5_caltet33_ani(mesh,met,pt[i]);
    }
  }
  else {
    for (i=0; i<8; i++) {
      pt[i]->qual=MMG5_orcal(mesh,met,newtet[i]);
    }
  }
  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ip index of new point.
 * \param list pointer toward the shell of edge.
 * \param ret size of the shell of edge.
 * \param crit quality threshold.
 * \return 0 if fail, 1 otherwise.
 *
 * Check quality before split.
 *
 */
static inline
int MMG3D_chksplit(MMG5_pMesh mesh, MMG5_pSol met,int ip,
                    int* list,int ret,double crit) {
  MMG5_pTetra   pt0,pt1;
  double        cal,critloc;
  int           l,jel,na,ipb,lon;

  lon = ret/2;
  critloc = 1.;
  for (l=0; l<lon; l++) {
    jel = list[l] / 6;
    pt1 = &mesh->tetra[jel];
    if(pt1->qual < critloc) critloc = pt1->qual;
  }
  critloc *= crit;

  pt0  = &mesh->tetra[0];
  for (l=0; l<lon; l++) {
    jel = list[l] / 6;
    na  = list[l] % 6;
    pt1 = &mesh->tetra[jel];

    memcpy(pt0->v,pt1->v,4*sizeof(int));
    ipb = MMG5_iare[na][0];
    pt0->v[ipb] = ip;
    cal = MMG5_caltet(mesh,met,pt0);
    if ( cal < critloc ) {
      MMG3D_delPt(mesh,ip);
      return 0;
    }

    memcpy(pt0->v,pt1->v,4*sizeof(int));
    ipb = MMG5_iare[na][1];
    pt0->v[ipb] = ip;
    cal = MMG5_caltet(mesh,met,pt0);
    if ( cal < critloc ) {
      MMG3D_delPt(mesh,ip);
      return 0;
    }
  }
  return 1;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param iel tetra index
 * \param iar edge index of iel
 * \param crit quality threshold.
 * \return -1 if lack of memory, 0 if we don't split the edge, ip if success.
 *
 * Split edge iar of iel and verify that every new tet have a better quality than crit
 *
 */
int MMG5_splitedg(MMG5_pMesh mesh, MMG5_pSol met,int iel, int iar, double crit){
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_pPoint  p0,p1;
  double       o[3];
  int          list[MMG3D_LMAX+2],i0,i1,ip,warn,lon,ier;
  int16_t      tag;

  warn = 0;
  pt = &mesh->tetra[iel];
  lon = MMG5_coquil(mesh,iel,iar,list);
  if ( (!lon || lon<0) )
    return 0;
  if(lon%2) return 0;

  i0 = pt->v[MMG5_iare[iar][0]];
  i1 = pt->v[MMG5_iare[iar][1]];
  p0  = &mesh->point[i0];
  p1  = &mesh->point[i1];

  tag = MG_NOTAG;
  if ( pt->xt ){
    pxt  = &mesh->xtetra[pt->xt];
    if ( (pxt->ftag[MMG5_ifar[iar][0]] & MG_BDY) ||
         (pxt->ftag[MMG5_ifar[iar][1]] & MG_BDY) ) {
      tag  = pxt->tag[iar];
      tag |= MG_BDY;
    }
  }

  if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) {
    return 0;
  }

  o[0] = 0.5*(p0->c[0] + p1->c[0]);
  o[1] = 0.5*(p0->c[1] + p1->c[1]);
  o[2] = 0.5*(p0->c[2] + p1->c[2]);

  ip = MMG3D_newPt(mesh,o,tag);

  if ( !ip )  {
    assert ( mesh );
    /* reallocation of point table */
    MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                         warn=1;
                         break
                         ,o,tag);
  }

  if ( warn ) {
    fprintf(stderr,"\n  ## Warning: %s:",__func__);
    fprintf(stderr," unable to allocate a new point in last call"
            " of MMG5_adpspl.\n");
    MMG5_INCREASE_MEM_MESSAGE();
  }

  ier = MMG5_intmet(mesh,met,iel,iar,ip,0.5);
  if ( !ier ) {
    MMG3D_delPt(mesh,ip);
    return 0;
  }
  else if (ier < 0 ) {
    MMG3D_delPt(mesh,ip);
    return 0;
  }

  ier = MMG3D_simbulgept(mesh,met,list,lon,ip);
  assert ( (!mesh->info.ddebug) || (mesh->info.ddebug && ier != -1) );
  if ( ier <= 0 || ier == 2 ) return 0;

  ier = MMG3D_chksplit(mesh,met,ip,&list[0],lon,crit);
  if(!ier) return 0;
  ier = MMG5_split1b(mesh,met,list,lon,ip,0,1,0);
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to split.\n",__func__);
    return -1;
  }
  else if ( !ier ) {
    MMG3D_delPt(mesh,ip);
    return 0;
  }

  return ip;
}
