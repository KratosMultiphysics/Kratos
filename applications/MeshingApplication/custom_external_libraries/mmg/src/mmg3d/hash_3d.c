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
 * \file mmg3d/hash_3d.c
 * \brief Functions for hash tables management and tetrahedra packing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"

#define MMG5_KC    13

extern int8_t  ddb;

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 1 if success, 0 if fail
 *
 * tetra packing.
 *
 */
int MMG5_paktet(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  MMG5_int      k;

  k = 1;
  do {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) {
      pt1 = &mesh->tetra[mesh->ne];
      assert( pt && pt1 && MG_EOK(pt1) );
      memcpy(pt,pt1,sizeof(MMG5_Tetra));
      if ( !MMG3D_delElt(mesh,mesh->ne) )  return 0;
    }
  }
  while ( ++k < mesh->ne );

  /* Recreate nil chain */
  assert(mesh->ne<=mesh->nemax);

  if ( mesh->ne == mesh->nemax )
    mesh->nenil = 0;
  else {
    mesh->nenil = mesh->ne + 1;

    for(k=mesh->nenil; k<=mesh->nemax-1; k++){
      mesh->tetra[k].v[3] = k+1;
    }

    mesh->tetra[mesh->nemax].v[3] = 0;
  }
  return 1;
}

/** return index of triangle ia ib ic */
MMG5_int MMG5_hashGetFace(MMG5_Hash *hash,MMG5_int ia,MMG5_int ib,MMG5_int ic) {
  MMG5_hedge  *ph;
  MMG5_int    key;
  MMG5_int    mins,maxs,sum;

  if ( !hash->item )  return 0;

  mins = MG_MIN(ia,MG_MIN(ib,ic));
  maxs = MG_MAX(ia,MG_MAX(ib,ic));

  /* compute key */
  sum = ia + ib + ic;
  key = (MMG5_KA*(int64_t)mins + MMG5_KB*(int64_t)maxs) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == mins && ph->b == maxs && ph->s == sum )
      return ph->k;
    else {
      while ( ph->nxt ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return ph->k;
      }
    }
  }

  return 0;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param pack we pack the mesh at function begining if \f$pack=1\f$.
 * \return 0 if failed, 1 otherwise.
 *
 * Create table of adjacency. Set pack variable to 0 for a compact
 * mesh and to 1 for a mesh that need to be packed.
 *
 */
int MMG3D_hashTetra(MMG5_pMesh mesh, int pack) {
  MMG5_pTetra    pt,pt1;
  MMG5_int       key;
  MMG5_int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,sum,sum1,iadr;
  MMG5_int       *hcode,*link,hsize,inival;
  uint8_t        i,ii,i1,i2,i3;

  /* default */
  if ( mesh->adja ) {
    return 1;
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING STRUCTURE\n");

  /* packing : if not hash does not work */
  if ( pack )  {
    if ( ! MMG5_paktet(mesh) ) return 0;
  }

  /* memory alloc */
  MMG5_ADD_MEM(mesh,(4*mesh->nemax+5)*sizeof(MMG5_int),"adjacency table",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->adja,4*mesh->nemax+5,MMG5_int,return 0);
  MMG5_SAFE_CALLOC(hcode,mesh->ne+5,MMG5_int,return 0);

  link  = mesh->adja;
  hsize = mesh->ne;

  /* init */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");

  if ( sizeof(MMG5_int) == 8 ) {
    inival = LONG_MAX;
  }
  else {
    inival = INT_MAX;
  }

  iadr   = 0;
  for (k=0; k<=mesh->ne; k++)
    hcode[k] = -inival;

  /* hash tetras */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<4; i++) {
      i1 = MMG5_idir[i][0];
      i2 = MMG5_idir[i][1];
      i3 = MMG5_idir[i][2];
      mins = MG_MIN(pt->v[i1],MG_MIN(pt->v[i2],pt->v[i3]));
      maxs = MG_MAX(pt->v[i1],MG_MAX(pt->v[i2],pt->v[i3]));

      /* compute key and insert */
      sum = pt->v[i1] + pt->v[i2] + pt->v[i3];
      key = (MMG5_KA*(int64_t)mins+MMG5_KB*(int64_t)maxs+MMG5_KC*(int64_t)sum)%hsize+1;
      iadr++;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 2: adjacencies\n");
  for (l=iadr; l>0; l--) {
    if ( link[l] >= 0 )  continue;

    /* current element */
    k = (l-1) / 4 + 1;
    i = (l-1) % 4;
    i1 = MMG5_idir[i][0];
    i2 = MMG5_idir[i][1];
    i3 = MMG5_idir[i][2];
    pt = &mesh->tetra[k];
    mins = MG_MIN(pt->v[i1],MG_MIN(pt->v[i2],pt->v[i3]));
    maxs = MG_MAX(pt->v[i1],MG_MAX(pt->v[i2],pt->v[i3]));
    sum  = pt->v[i1] + pt->v[i2] + pt->v[i3];

    /* accross link */
    ll      = -link[l];
    pp      = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 4 + 1;
      ii = (ll-1) % 4;
      i1 = MMG5_idir[ii][0];
      i2 = MMG5_idir[ii][1];
      i3 = MMG5_idir[ii][2];
      pt1  = &mesh->tetra[kk];
      sum1 = pt1->v[i1] + pt1->v[i2] + pt1->v[i3];
      if ( sum1 == sum ) {
        mins1 = MG_MIN(pt1->v[i1],MG_MIN(pt1->v[i2],pt1->v[i3]));
        maxs1 = MG_MAX(pt1->v[i1],MG_MAX(pt1->v[i2],pt1->v[i3]));

        /* adjacent found */
        if ( mins1 == mins && maxs1 == maxs ) {
          if ( pp != 0 )  link[pp] = link[ll];
          link[l]  = 4*kk + ii;
          link[ll] = 4*k + i;
          break;
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  MMG5_SAFE_FREE(hcode);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Create partial table of adjacency for prisms (prism \f$ <-> \f$ prism).
 *
 * \remark Adjacencies between prisms and tetra are not filled here.
 *
 * \warning check the hashtable efficiency
 */
int MMG3D_hashPrism(MMG5_pMesh mesh) {
  MMG5_pPrism    pp,pp1;
  MMG5_int       key;
  MMG5_int       k,kk,l,ll,jj;
  MMG5_int       max12,min12,max34,min34,mins,mins1,mins_b, mins_b1,maxs,maxs1;
  MMG5_int       iadr;
  MMG5_int       *hcode,*link,hsize,inival;
  uint8_t        i,ii,i1,i2,i3,i4;

  if ( !mesh->nprism ) return 1;

  /* default */
  if ( mesh->adjapr ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: no re-build of adjacencies of prisms. "
              "mesh->adjapr must be freed to enforce analysis.\n",__func__);
    }
    return 1;
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING PRISMS ADJACENCY\n");

  /* memory alloc */
  MMG5_ADD_MEM(mesh,(5*mesh->nprism+6)*sizeof(MMG5_int),"prism adjacency table",
                printf("  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->adjapr,5*mesh->nprism+6,MMG5_int,return 0);
  MMG5_SAFE_CALLOC(hcode,mesh->nprism+6,MMG5_int,return 0);

  link  = mesh->adjapr;
  hsize = mesh->nprism;

  /* init */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");

  if ( sizeof(MMG5_int) == 8 ) {
    inival = LONG_MAX;
  }
  else {
    inival = INT_MAX;
  }

  iadr   = 0;
  for (k=0; k<=mesh->nprism; k++)
    hcode[k] = -inival;

  /* hash prism */
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    assert ( MG_EOK(pp) );
    for (i=0; i<2; i++) {
      /* Triangular face */
      i1 = MMG5_idir_pr[i][0];
      i2 = MMG5_idir_pr[i][1];
      i3 = MMG5_idir_pr[i][2];

      min12 = MG_MIN(pp->v[i1],pp->v[i2]);
      /* mins = minimum index of triangle vertices */
      mins = MG_MIN(min12,pp->v[i3]);

      max12  = MG_MAX(pp->v[i1],pp->v[i2]);
      /* maxs = maximum index of triangle vertices */
      maxs   = MG_MAX(max12,pp->v[i3]);

      /* mins_b = second minimum index of triangle vertices */
      mins_b = pp->v[i1] + pp->v[i2] + pp->v[i3] -mins -maxs;

      /* compute key and insert */
      key = (MMG5_KA*(int64_t)mins+MMG5_KB*(int64_t)mins_b+MMG5_KC*(int64_t)maxs)%hsize+1;
      iadr++;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
    for ( ; i<5; ++i) {
      /* Quadrilateral face */
      i1 = MMG5_idir_pr[i][0];
      i2 = MMG5_idir_pr[i][1];
      i3 = MMG5_idir_pr[i][2];
      i4 = MMG5_idir_pr[i][3];

      min12 = MG_MIN(pp->v[i1],pp->v[i2]);
      min34 = MG_MIN(pp->v[i3],pp->v[i4]);
      /* mins = minimum index of quadrilateral vertices */
      mins = MG_MIN(min12,min34);

      max12  = MG_MAX(pp->v[i1],pp->v[i2]);
      max34  = MG_MAX(pp->v[i3],pp->v[i4]);
      /* maxs = maximum index of quadrilateral vertices */
      maxs   = MG_MAX(max12,max34);

      /* mins_b = second minimum index of quadrilateral vertices */
      mins_b = MG_MIN( MG_MIN(max12,max34),MG_MAX(min12,min34));

      /* compute key and insert */
      key = (MMG5_KA*(int64_t)mins+MMG5_KB*(int64_t)mins_b+MMG5_KC*(int64_t)maxs)%hsize+1;
      iadr++;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 2: adjacencies\n");
  for (l=iadr; l>0; l--) {
    if ( link[l] >= 0 )  continue;

    /* current element */
    k = (l-1) / 5 + 1;
    i = (l-1) % 5;

    switch (i)
    {
    case 0:
    case 1:
      i1 = MMG5_idir_pr[i][0];
      i2 = MMG5_idir_pr[i][1];
      i3 = MMG5_idir_pr[i][2];
      pp = &mesh->prism[k];

      min12 = MG_MIN(pp->v[i1],pp->v[i2]);
      mins  = MG_MIN(min12,pp->v[i3]);

      max12 = MG_MAX(pp->v[i1],pp->v[i2]);
      maxs  = MG_MAX(max12,pp->v[i3]);

      mins_b = pp->v[i1] + pp->v[i2] + pp->v[i3] - mins - maxs;

      break;

    default:
      i1 = MMG5_idir_pr[i][0];
      i2 = MMG5_idir_pr[i][1];
      i3 = MMG5_idir_pr[i][2];
      i4 = MMG5_idir_pr[i][3];
      pp = &mesh->prism[k];

      min12 = MG_MIN(pp->v[i1],pp->v[i2]);
      min34 = MG_MIN(pp->v[i3],pp->v[i4]);
      mins  = MG_MIN(min12,min34);

      max12 = MG_MAX(pp->v[i1],pp->v[i2]);
      max34 = MG_MAX(pp->v[i3],pp->v[i4]);
      maxs  = MG_MAX(max12,max34);

      mins_b = MG_MIN( MG_MIN(max12,max34),MG_MAX(min12,min34));
    }

    /* accross link */
    ll      = -link[l];
    jj      = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 5 + 1;
      ii = (ll-1) % 5;

      switch (ii)
      {
      case 0:
      case 1:
        i1 = MMG5_idir_pr[ii][0];
        i2 = MMG5_idir_pr[ii][1];
        i3 = MMG5_idir_pr[ii][2];

        pp1  = &mesh->prism[kk];

        min12 = MG_MIN(pp1->v[i1],pp1->v[i2]);
        mins1  = MG_MIN(min12,pp1->v[i3]);

        max12 = MG_MAX(pp1->v[i1],pp1->v[i2]);
        maxs1  = MG_MAX(max12,pp1->v[i3]);

        mins_b1 =  pp1->v[i1] + pp1->v[i2] + pp1->v[i3] - mins1 - maxs1;

        break;

      default:
        i1 = MMG5_idir_pr[ii][0];
        i2 = MMG5_idir_pr[ii][1];
        i3 = MMG5_idir_pr[ii][2];
        i4 = MMG5_idir_pr[ii][3];
        pp1  = &mesh->prism[kk];

        min12  = MG_MIN(pp1->v[i1],pp1->v[i2]);
        min34  = MG_MIN(pp1->v[i3],pp1->v[i4]);
        mins1  = MG_MIN(min12,min34);

        max12  = MG_MAX(pp1->v[i1],pp1->v[i2]);
        max34  = MG_MAX(pp1->v[i3],pp1->v[i4]);
        maxs1  = MG_MAX(max12,max34);

        mins_b1 = MG_MIN( MG_MIN(max12,max34),MG_MAX(min12,min34));
      }

      /* adjacent found */
      if ( mins1 == mins && maxs1 == maxs && mins_b1 == mins_b ) {
        if ( jj != 0 )  link[jj] = link[ll];
        link[l]  = 5*kk + ii;
        link[ll] = 5*k + i;
        break;
      }

      jj = ll;
      ll = -link[ll];
    }
  }
  MMG5_SAFE_FREE(hcode);
  return 1;
}

/**
 * \param mesh pointer towar the mesh structure.
 * \param hash edges hash table.
 * \return 1 if success, 0 if failed.
 *
 * Set non-manifold tag at extremities of a non-manifold edge.  Check if a
 * non-manifold edge is at the interface of several distinct domains (thus we
 * can't travel from a domain through another one by adjacency) and in this
 * case, mark the edge and its extremities as required. Free the edge hash table
 * \a hash if success.
 *
 * \warning if fail, the edge hash table \a hash is not freed.
 *
 */
static inline
int MMG5_setEdgeNmTag(MMG5_pMesh mesh, MMG5_Hash *hash) {
  MMG5_pTetra         pt;
  MMG5_pxTetra        pxt;
  MMG5_pTria          ptt;
  MMG5_hedge          *ph;
  MMG5_int            adj,pradj,piv;
  int64_t             list[MMG3D_LMAX+2];
  MMG5_int            key;
  MMG5_int            k,l,i1,i2,na,nb,ia,it1,it2, nr;
  MMG5_int            start;
  int                 ilist,nbdy,ipa,ipb;
  int8_t              iface,hasadja,i;
  static int8_t       mmgWarn0=0,mmgWarn1=0;

  nr = 0;

  /* First: seek edges at the interface of two distinct domains and mark it as
   * required */
  for (k=1; k<=mesh->nt; k++) {
    ptt  = &mesh->tria[k];

    if ( !MG_EOK(ptt) ) continue;

    for (l=0; l<3; l++) {

      /* Skip parallel edges */
      if ( (ptt->tag[l] & MG_PARBDY) || (ptt->tag[l] & MG_BDY) ) continue;

      if ( ptt->tag[l] & MG_NOM ) {
        i1 = MMG5_inxt2[l];
        i2 = MMG5_iprv2[l];

        /* compute key */
        na  = MG_MIN(ptt->v[i1],ptt->v[i2]);
        nb  = MG_MAX(ptt->v[i1],ptt->v[i2]);
        key = (MMG5_KA*na + MMG5_KB*nb) % hash->siz;
        ph  = &hash->item[key];

        assert(ph->a);
        while ( ph->a ) {
          if ( ph->a == na && ph->b == nb ) break;
          assert(ph->nxt);
          ph = &hash->item[ph->nxt];
        }
        /* Set edge tag and point tags to MG_REQ if the non-manifold edge shared
         * separated domains */
        if ( ph->s > 3 ) {
          start = ptt->cc/4;
          assert(start);
          pt = &mesh->tetra[start];


          for (ia=0; ia<6; ++ia) {
            ipa = MMG5_iare[ia][0];
            ipb = MMG5_iare[ia][1];
            if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
                 (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
          }
          assert(ia<6);


          /* Travel throug the shell of the edge until reaching a tetra without adjacent
           * or until reaching the starting tetra */
          iface = ptt->cc%4;
          MMG3D_coquilFaceFirstLoop(mesh,start,na,nb,iface,ia,list,&ilist,&it1,&it2,
                                    &piv,&adj,&hasadja,&nbdy,1);

          /* At this point, the first travel, in one direction, of the shell is
             complete. Now, analyze why the travel ended. */
          if ( adj == start ) {
            if ( !it2 ) {
              if ( !mmgWarn0 ) {
                mmgWarn0 = 1;
                fprintf(stderr,"\n  ## Warning: %s: at least 1 wrong boundary tag:"
                       " Only 0 or 1 boundary triangles founded in the shell of the edge\n",
                       __func__);
              }
            }
            if ( nbdy < 2 )
              MMG5_coquilFaceErrorMessage(mesh, it1/4, it2/4);
          }
          else {
            /* A boundary has been detected : slightly different configuration */
            if ( hasadja ) {

              /* Start back everything from this tetra adj */
              MMG3D_coquilFaceSecondLoopInit(mesh,piv,&iface,&i,list,&ilist,&it1,
                                              &pradj,&adj);

              nbdy = 1;
              while ( adj ) {
                pradj = adj;

                if ( MMG5_openCoquilTravel(mesh,na,nb,&adj,&piv,&iface,&i)<0 ) {
                  return 0;
                }

                /* overflow */
                if ( ++ilist > MMG3D_LMAX-2 ) {
                  if ( !mmgWarn1 ) {
                    mmgWarn1 = 1;
                    fprintf(stderr,"\n  ## Warning: %s: problem in surface remesh"
                            " process. At least 1 shell of edge (%" MMG5_PRId "-%" MMG5_PRId ") contains"
                            " too many elts.\n",__func__,MMG3D_indPt(mesh,na),
                            MMG3D_indPt(mesh,nb));
                    fprintf(stderr,"\n  ##          Try to modify the hausdorff"
                            " number, or/and the maximum mesh.\n");
                  }
                  return 0;
                }

                pt = &mesh->tetra[pradj];
                if ( pt->xt ) {
                  pxt = &mesh->xtetra[pt->xt];
                  if ( pxt->ftag[iface] & MG_BDY ) ++nbdy;
                }
              }

              assert(!adj);
              it2 = 4*pradj + iface;

              if ( (!it1 || !it2) || (it1 == it2) ) {
                MMG5_coquilFaceErrorMessage(mesh, it1/4, it2/4);
                return 0;
              }
            }
          }

          /* If ph->s do not match the number of encountred boundaries we have
             separated domains. */
          if ( nbdy != ph->s ) {
            if ( !(ptt->tag[l] & MG_REQ) ) {
              ptt->tag[l] |= MG_REQ;
              ptt->tag[l] &= ~MG_NOSURF;
              ++nr;
            }
            mesh->point[ptt->v[MMG5_inxt2[l]]].tag |= MG_REQ;
            mesh->point[ptt->v[MMG5_iprv2[l]]].tag |= MG_REQ;
            mesh->point[ptt->v[MMG5_inxt2[l]]].tag &= ~MG_NOSURF;
            mesh->point[ptt->v[MMG5_iprv2[l]]].tag &= ~MG_NOSURF;
          }

          /* Work done for this edge: reset ph->s/ */
          ph->s = 0;
        }
      }
    }
  }
  if ( mesh->info.ddebug || abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"     %" MMG5_PRId " required edges added\n",nr);

  /* Free the edge hash table */
  MMG5_DEL_MEM(mesh,hash->item);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 1 if success, 0 if fail.
 * Seek the non-required non-manifold points and try to analyse whether they are
 * corner or required.
 *
 * \remark We don't know how to travel through the shell of a non-manifold point
 * by triangle adjacency. Thus the work done here can't be performed in the \ref
 * MMG5_singul function.
 */
static inline
int MMG5_setVertexNmTag(MMG5_pMesh mesh) {
  MMG5_pTetra         ptet;
  MMG5_pPoint         ppt;
  MMG5_Hash           hash;
  int                 i,ier;
  MMG5_int            k,base,np,nc,nm,nre, ng, nrp;

  /* Second: seek the non-required non-manifold points and try to analyse
   * whether they are corner or required. */

  /* Hash table used by boulernm to store the special edges passing through
   * a given point */
  np = 0;
  for (k=1; k<=mesh->np; ++k) {
    ppt = &mesh->point[k];
    if ( (!(ppt->tag & MG_NOM)) || (ppt->tag & MG_REQ) ) continue;
    ++np;
  }

  if ( ! MMG5_hashNew(mesh,&hash,np,(MMG5_int)(3.71*np)) ) return 0;

  nc = nre = 0;
  base = ++mesh->base;
  for (k=1; k<=mesh->ne; ++k) {
    ptet = &mesh->tetra[k];
    if ( !MG_EOK(ptet) ) continue;

    for ( i=0; i<4; ++i ) {
      ppt = &mesh->point[ptet->v[i]];

      /* Skip parallel points */
      if ( ppt->tag & MG_PARBDY ) continue;

      if ( (!MG_VOK(ppt)) || (ppt->flag==base)  ) continue;
      ppt->flag = base;

      if ( (!(ppt->tag & MG_NOM)) || (ppt->tag & MG_REQ) ) continue;

      ier = MMG5_boulernm(mesh,&hash, k, i, &ng, &nrp, &nm);
      if ( ier < 0 ) return 0;
      else if ( !ier ) continue;

      if ( (ng+nrp+nm) > 2 ) {
        /* More than 2 feature edges are passing through the point: point is
         * marked as corner */
        ppt->tag |= MG_CRN + MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
        nc++;
      }
      else if ( (ng == 2) || (nrp == 2) || (nm == 2) ) {
        /* Exactly 2 edges of same type are passing through the point: do
         * nothing */
        continue;
      }
      else if ( (ng+nrp+nm) == 2 ) {
        /* 2 edges of different type are passing through the point: point is
         * marked as required */
        ppt->tag |= MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
      }
      else if ( ng == 1 && !nrp ){
        ppt->tag |= MG_CRN + MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
        nc++;
      }
      else if ( (ng+nrp+nm) == 1 ){
        /* Only 1 feature edge is passing through the point: point is
         * marked as corner */
        assert ( (ng == 1) || (nrp==1) || (nm==1) );
        ppt->tag |= MG_CRN + MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
        nc++;
      }
      else {
        assert ( 0 && "unexpected case");
      }
    }
  }

  /* Free the edge hash table */
  MMG5_DEL_MEM(mesh,hash.item);

  if ( mesh->info.ddebug || abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"     %" MMG5_PRId " corner and %" MMG5_PRId " required vertices added\n",nc,nre);

  return 1;
}

/**
 * \param mesh pointer towar the mesh structure.
 * \param hash edges hash table.
 * \return 1 if success, 0 if failed.
 *
 * Set tags to non-manifold edges and vertices. Not done before because we need
 * the \ref MMG5_xTetra table.
 *
 * \warning if fail, the edge hash table \a hash is not freed.
 *
 */
int MMG5_setNmTag(MMG5_pMesh mesh, MMG5_Hash *hash) {

  /* First: seek edges at the interface of two distinct domains and mark it as
   * required */
  if ( !MMG5_setEdgeNmTag(mesh,hash) ) return 0;

  /* Second: seek the non-required non-manifold points and try to analyse
   * whether they are corner or required. */
  if ( !MMG5_setVertexNmTag(mesh) ) return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash Edges hash table.
 * \return 1 if success, 0 if failed.
 *
 * Create surface adjacency table. Allocate the edge hash table \a hash but
 * don't free it.
 *
 */
int MMG3D_hashTria(MMG5_pMesh mesh, MMG5_Hash *hash) {

  MMG5_DEL_MEM(mesh,mesh->adjt);

  MMG5_ADD_MEM(mesh,(3*mesh->nt+4)*sizeof(MMG5_int),"surfacic adjacency table",return 0);
  MMG5_SAFE_CALLOC(mesh->adjt,3*mesh->nt+4,MMG5_int,return 0);

  return  MMG5_mmgHashTria(mesh, mesh->adjt, hash, mesh->info.iso) ;
}


/** remove edge from hash table */
int MMG5_hashPop(MMG5_Hash *hash,MMG5_int a,MMG5_int b) {
  MMG5_hedge  *ph,*php;
  MMG5_int    key;
  MMG5_int    ia,ib,iph,iphp;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->a ) return 0;
  else if ( ph->a == ia && ph->b == ib ) {
    if ( !ph->nxt ) {
      memset(ph,0,sizeof(MMG5_hedge));
      return 1;
    }
    else {
      iph = ph->nxt;
      php = ph;
      ph  = &hash->item[ph->nxt];
      memcpy(php,ph,sizeof(MMG5_hedge));
      memset(ph,0,sizeof(MMG5_hedge));
      ph->nxt   = hash->nxt;
      hash->nxt = iph;
      return 1;
    }
  }
  while ( ph->nxt ) {
    php = ph;
    ph  = &hash->item[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      if ( !ph->nxt ) {
        memset(ph,0,sizeof(MMG5_hedge));
        ph->nxt   = hash->nxt;
        hash->nxt = php->nxt;
        php->nxt  = 0;
      }
      else {
        iph  = ph->nxt;
        iphp = php->nxt;
        php->nxt = iph;
        memset(ph,0,sizeof(MMG5_hedge));
        ph->nxt   = hash->nxt;
        hash->nxt = iphp;
      }
      return 1;
    }
  }
  return 0;
}


/**
 * \param hash pointer toward the hash table in which edges are stored
 * \param a first edge extremity
 * \param b second edge extremity
 * \param ref reference to assign to the edge
 * \param tag tag to assign
 *
 * \return 0 if the edge is not in the hash table, 1 otherwise
 *
 * set tag to edge on geometry
 *
 */
int MMG5_hTag(MMG5_HGeom *hash,MMG5_int a,MMG5_int b,MMG5_int ref,int16_t tag) {
  MMG5_hgeom  *ph;
  MMG5_int    key;
  MMG5_int    ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )
    return 0;
  else if ( ph->a == ia && ph->b == ib ) {
    ph->tag |= tag;
    if ( ref ) {
      ph->ref  = ref;
    }
    return 1;
  }
  while ( ph->nxt ) {
    ph = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      ph->tag |= tag;
      if ( ref ) {
        ph->ref  = ref;
      }
      return 1;
    }
  }
  return 0;
}

/** remove edge from hash table */
int MMG5_hPop(MMG5_HGeom *hash,MMG5_int a,MMG5_int b,MMG5_int *ref,int16_t *tag) {
  MMG5_hgeom  *ph,*php;
  MMG5_int    key;
  MMG5_int    ia,ib,iph,iphp;

  *ref = 0;
  *tag = 0;

  assert ( hash->siz );

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )  return 0;
  else if ( ph->a == ia && ph->b == ib ) {
    *ref = ph->ref;
    *tag = ph->tag;
    if ( !ph->nxt ) {
      memset(ph,0,sizeof(MMG5_hgeom));
    }
    else {
      iph = ph->nxt;
      php = ph;
      ph  = &hash->geom[ph->nxt];
      memcpy(php,ph,sizeof(MMG5_hgeom));
      memset(ph,0,sizeof(MMG5_hgeom));
      ph->nxt   = hash->nxt;
      hash->nxt = iph;
    }
    return 1;
  }
  while ( ph->nxt ) {
    php = ph;
    ph  = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      *ref = ph->ref;
      *tag = ph->tag;
      if ( !ph->nxt ) {
        memset(ph,0,sizeof(MMG5_hgeom));
        ph->nxt   = hash->nxt;
        hash->nxt = php->nxt;
        php->nxt  = 0;
      }
      else {
        iph  = ph->nxt;
        iphp = php->nxt;
        php->nxt = iph;
        memset(ph,0,sizeof(MMG5_hgeom));
        ph->nxt   = hash->nxt;
        hash->nxt = iphp;
      }
      return 1;
    }
  }
  return 0;
}

/** get ref and tag to edge on geometry */
int MMG5_hGet(MMG5_HGeom *hash,MMG5_int a,MMG5_int b,MMG5_int *ref,int16_t *tag) {
  MMG5_hgeom  *ph;
  MMG5_int    key;
  MMG5_int    ia,ib;

  *tag = 0;
  *ref = 0;

  assert ( hash->siz );

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( !ph->a )  return 0;
  else if ( ph->a == ia && ph->b == ib ) {
    *ref = ph->ref;
    *tag = ph->tag;
    return 1;
  }
  while ( ph->nxt ) {
    ph = &hash->geom[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      *ref = ph->ref;
      *tag = ph->tag;
      return 1;
    }
  }
  return 0;
}

/** store edge on geometry */
int MMG5_hEdge(MMG5_pMesh mesh,MMG5_HGeom *hash,MMG5_int a,MMG5_int b,MMG5_int ref,int16_t tag) {
  MMG5_hgeom  *ph;
  MMG5_int    key;
  MMG5_int    ia,ib,j;

  assert ( hash->siz );

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->geom[key];

  if ( ph->a == ia && ph->b == ib )
    return 1;
  else if ( ph->a ) {
    while ( ph->nxt ) {
      ph = &hash->geom[ph->nxt];
      if ( ph->a == ia && ph->b == ib )  return 1;
    }
    ph->nxt = hash->nxt;
    ph      = &hash->geom[hash->nxt];
    ph->a   = ia;   ph->b   = ib;
    ph->ref = ref;  ph->tag = tag;
    hash->nxt = ph->nxt;
    ph->nxt = 0;
    if ( hash->nxt >= hash->max ) {
      if ( mesh->info.ddebug )
        fprintf(stderr,"\n  ## Memory alloc problem (edge): %" MMG5_PRId "\n",hash->max);
      MMG5_TAB_RECALLOC(mesh,hash->geom,hash->max,MMG5_GAP,MMG5_hgeom,
                         "larger htab table",
                         fprintf(stderr,"  Exit program.\n");return 0;);
      for (j=hash->nxt; j<hash->max; j++)  hash->geom[j].nxt = j+1;
    }
    return 1;
  }
  /* insert new edge */
  ph->a   = ia;   ph->b   = ib;
  ph->ref = ref;  ph->tag = tag;
  ph->nxt = 0;
  return 1;
}

/** to store edge on geometry */
int MMG5_hNew(MMG5_pMesh mesh,MMG5_HGeom *hash,MMG5_int hsiz,MMG5_int hmax) {
  MMG5_int   k;

  /* adjust hash table params */
  hash->siz  = hsiz + 1;
  hash->max  = hmax + 2;
  hash->nxt  = hash->siz;

  MMG5_ADD_MEM(mesh,(hash->max+1)*sizeof(MMG5_hgeom),"Edge hash table",return 0);
  MMG5_SAFE_CALLOC(hash->geom,(hash->max+1),MMG5_hgeom,return 0);

  if ( !hash->geom ) {
    perror("  ## Memory problem: calloc");
    return 0;
  }
  for (k=hash->siz; k<hash->max; k++)
    hash->geom[k].nxt = k+1;

  return 1;
}

/**
 * \param mesh pointer toward he mesh structure.
 * \return 0 if failed, 1 otherwise
 *
 * Build hashtable for initial mesh edges.
 *
 */
int MMG5_hGeom(MMG5_pMesh mesh) {
  MMG5_pTria   pt;
  MMG5_pEdge   pa;
  MMG5_Hash    hash;
  MMG5_int     edg,*adja,k,kk;
  int          ier;
  int16_t      tag;
  int8_t       i,i1,i2;

  /* if edges exist in mesh, hash special edges from existing field */
  if ( mesh->na ) {
    if ( !mesh->htab.geom ) {
      mesh->namax = MG_MAX((MMG5_int)(1.5*mesh->na),MMG3D_NAMAX);
      if ( !MMG5_hNew(mesh,&mesh->htab,mesh->na,3*mesh->namax) )
        return 0;
    }
    else {
      if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
        fprintf(stderr,"\n  ## Warning: %s: no re-hash of edges of mesh. ",
                __func__);
        fprintf(stderr,"mesh->htab.geom must be freed to enforce analysis.\n");
      }
      MMG5_DEL_MEM(mesh,mesh->edge);
      mesh->na   = 0;
      return 1;
    }

    /* store initial edges */
    for (k=1; k<=mesh->na; k++) {
      pa = &mesh->edge[k];
      if ( !MMG5_hEdge(mesh,&mesh->htab,pa->a,pa->b,pa->ref,pa->tag) )
        return 0;
    }

    /* now check triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        if( (pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY) ) continue;
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        /* transfer non manifold tag to edges */
        if ( pt->tag[i] & MG_NOM ) {
          ier = MMG5_hTag(&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
          if ( !ier ) {
            /* The edge is marked as non manifold but doesn't exist in the mesh */
            ier = MMG5_hEdge(mesh,&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]);
            if ( !ier )
              return 0;
          }
        }
        MMG5_hGet(&mesh->htab,pt->v[i1],pt->v[i2],&edg,&tag);
        pt->edg[i]  = edg;

        /* If we use the nosurf option and the edge is required, we don't want
         * to detect it as an edge whose tag has been modified for the option */
        if ( mesh->info.nosurf && (tag & MG_REQ) )
          pt->tag[i] &= ~MG_NOSURF;

        /* Store the edge tag inside the triangle */
        pt->tag[i] |= tag;

        MMG5_hTag(&mesh->htab,pt->v[i1],pt->v[i2],edg,pt->tag[i]);
      }
    }
    MMG5_DEL_MEM(mesh,mesh->edge);
    mesh->na   = 0;
  }
  /* else, infer special edges from information carried by triangles */
  else {
    if ( !mesh->adjt ) {
      memset(&hash,0x0,sizeof(MMG5_Hash));
      ier = MMG3D_hashTria(mesh,&hash);
      MMG5_DEL_MEM(mesh,hash.item);
      if ( !ier ) return 0;
    }

    for (k=1; k<=mesh->nt; k++) {
      pt   = &mesh->tria[k];
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        if( (pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY) ) continue;
        kk  = adja[i] / 3;
        if ( !kk || pt->tag[i] & MG_NOM )
          mesh->na++;
         else if ( (k < kk) && ( pt->edg[i] || pt->tag[i] ) )  mesh->na++;
      }
    }

    if ( mesh->htab.geom )
      MMG5_DEL_MEM(mesh,mesh->htab.geom);

    mesh->namax = MG_MAX((MMG5_int)(1.5*mesh->na),MMG3D_NAMAX);
    if ( !MMG5_hNew(mesh,&mesh->htab,mesh->na,3*mesh->namax) )
      return 0;

    mesh->na = 0;

    /* build hash for edges */
    for (k=1; k<=mesh->nt; k++) {
      pt   = &mesh->tria[k];
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        if( (pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY) ) continue;
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        kk  = adja[i] / 3;
        if ( (!kk) || pt->tag[i] & MG_NOM ) {
          if ( pt->tag[i] & MG_NOM ) {
            if ( mesh->info.iso )
              pt->edg[i] = ( pt->edg[i] != 0 ) ?  -MMG5_abs(pt->edg[i]) : mesh->info.isoref;
          }
          if ( !MMG5_hEdge(mesh,&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]) )
            return 0;
        }
        else if ( k < kk && ( pt->edg[i] || pt->tag[i] ) ) {
          if ( !MMG5_hEdge(mesh,&mesh->htab,pt->v[i1],pt->v[i2],pt->edg[i],pt->tag[i]))
            return 0;
        }
      }
    }
    /* now check triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        if( (pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY) ) continue;
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        MMG5_hGet(&mesh->htab,pt->v[i1],pt->v[i2],&edg,&tag);
        pt->edg[i]  = edg;
        pt->tag[i] |= tag;
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param ntmesh number of boundary tria found in the mesh.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Fill tria array with missing triangles:
 * - if called from the first time, xtetra is not allocated so missing triangles
 *    may be founded at domains interface or external boundary.
 * - Otherwise, the tria array is reconstructed from the xtetra
 *   infos.
 *
 * \remark mesh->xtetra is not allocated when \ref MMG5_bdryTria is called by
 * \ref MMG3D_analys in mesh adp mode but it is allocated when called by
 * packMesh or by analys in ls discretization mode.
 *
 */
static inline
int MMG5_bdryTria(MMG5_pMesh mesh, MMG5_int ntmesh) {
  MMG5_pTetra    pt,pt1;
  MMG5_pPrism    pp;
  MMG5_pTria     ptt;
  MMG5_pPoint    ppt;
  MMG5_pxTetra   pxt;
  MMG5_pxPrism   pxpr;
  MMG5_Hash      hash;
  MMG5_int       ref,*adja,adj,k,ia,ib,ic,kt,ntinit;
  int            tofree=0;
  int8_t         i,j;

  hash.item = NULL;

  ntinit = mesh->nt;

  if  ( mesh->nprism && (ntmesh!=ntinit) ) {
    /* If a triangle at the interface between a prism and a tetra is not
     * provided, the hashtable is used to recover from the prism a boundary tria
     * created by tetra */
    if ( ! MMG5_hashNew(mesh,&hash,0.51*ntmesh,1.51*ntmesh) ) return 0;
    tofree=1;
  }
  else if ( mesh->nt ) {
    /* Hash given bdry triangles */
    if ( ! MMG5_hashNew(mesh,&hash,(MMG5_int)(0.51*mesh->nt),(MMG5_int)(1.51*mesh->nt)) ) return 0;
    tofree=1;
  }

  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ptt->v[i]];
      if ( !mesh->info.iso ) ppt->tag |= MG_BDY;
    }
  }

  /* Add boundary triangles stored on tetra */
  if ( ntmesh != ntinit ) {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      adja = &mesh->adja[4*(k-1)+1];
      pxt = 0;
      if ( pt->xt )  pxt = &mesh->xtetra[pt->xt];
      for (i=0; i<4; i++) {
        adj = adja[i] / 4;
        pt1 = &mesh->tetra[adj];

        if ( (!mesh->info.opnbdy) || (!mesh->xtetra) ) {
          /* Classic mode (no open bdy) or no info stored in xtetra (first
           * call) */
          if ( adj && ( pt->ref <= pt1->ref) )  continue;
        } else {
          /* info stored in xtetra and opnbdy mode */
          if ( adj && ( (pt->ref<pt1->ref) || (!pt->xt) ||
                        (!(pxt->ftag[i] & MG_BDY)) || (!MG_GET(pxt->ori,i) ) ) )
            continue;
        }

        ia = pt->v[MMG5_idir[i][0]];
        ib = pt->v[MMG5_idir[i][1]];
        ic = pt->v[MMG5_idir[i][2]];

        kt = MMG5_hashGetFace(&hash,ia,ib,ic);
        if ( kt ) {
          /* Face is already stored */
          continue;
        }
        else if ( mesh->nprism ) {
          /* Update the list of boundary trias to be able to recover tria at the
           * interface between tet and prisms */
          if ( !MMG5_hashFace(mesh,&hash,ia,ib,ic,mesh->nt+1) ) {
            MMG5_DEL_MEM(mesh,hash.item);
            return 0;
          }
        }

        /* face does not exists: add it in tria array */
        mesh->nt++;
        ptt = &mesh->tria[mesh->nt];
        ptt->v[0] = pt->v[MMG5_idir[i][0]];
        mesh->point[ptt->v[0]].tag |= MG_BDY;
        ptt->v[1] = pt->v[MMG5_idir[i][1]];
        mesh->point[ptt->v[1]].tag |= MG_BDY;
        ptt->v[2] = pt->v[MMG5_idir[i][2]];
        mesh->point[ptt->v[2]].tag |= MG_BDY;

        /* the cc field is used to be able to recover the tetra (and its face)
         * from which comes a boundary triangle (when called by packmesh =>
         * mesh->nt=0 at the beginning of the function) */
        ptt->cc = 4*k + i;

        if ( pxt ) {
          /* useful only when saving mesh or in ls mode */
          for( j = 0; j < 3; j++ ) {
            if ( pxt->tag[MMG5_iarf[i][j]] ) {
              ptt->tag[j] = pxt->tag[MMG5_iarf[i][j]];
              /* Remove redundant boundary tag */
              ptt->tag[j] &= ~MG_BDY;
            }
            if ( pxt->edg[MMG5_iarf[i][j]] )
              ptt->edg[j] = pxt->edg[MMG5_iarf[i][j]];
          }
        }
        if ( adj ) {
          if ( mesh->info.iso ) {
            /* Triangle at the interface between two tets is set to the user-defined ref if any, or else to mesh->info.isoref ref */
            if ( pxt && pxt->ftag[i] & MG_BDY )
              ptt->ref = pxt->ref[i];
            else if( MMG5_isLevelSet(mesh,pt->ref,pt1->ref) )
              ptt->ref = mesh->info.isoref;
            else
              ptt->ref = MG_MIN(pt->ref,pt1->ref);
          }
          /* useful only when saving mesh or in ls mode */
          else {
            /* Triangle at the interface between two tet is set its init ref or
             * to the min ref of the adja tetra */
            ptt->ref  = pxt ? pxt->ref[i] : MG_MIN(pt->ref,pt1->ref);
          }
        }
        else {
          /* useful only when saving mesh */
          ptt->ref  = pxt ? pxt->ref[i] : pt->ref;
        }
      }
    }
  }

  /* Add boundary triangles stored on prisms */
  if ( mesh->nprism ) {
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) )  continue;

      adja = &mesh->adjapr[5*(k-1)+1];
      pxpr = 0;
      if ( pp->xpr )  pxpr = &mesh->xprism[pp->xpr];

      for (i=0; i<2; i++) {
        adj = adja[i]/5;


        if ( adj < 0 ) {
          if ( !mesh->nt ) continue;

          /* Tria at the interface of a prism and a tetra: mark it as required */
          ia = pp->v[0+i*3];
          ib = pp->v[1+i*3];
          ic = pp->v[2+i*3];

          kt = MMG5_hashGetFace(&hash,ia,ib,ic);
          if ( !kt ) continue;

          ptt = &mesh->tria[kt];

          if ( !(ptt->tag[0] & MG_REQ) ) {
            ptt->tag[0] |= MG_REQ;
            ptt->tag[0] |= MG_NOSURF;
          }
          if ( !(ptt->tag[1] & MG_REQ) ) {
            ptt->tag[1] |= MG_REQ;
            ptt->tag[1] |= MG_NOSURF;
          }
          if ( !(ptt->tag[2] & MG_REQ) ) {
            ptt->tag[2] |= MG_REQ;
            ptt->tag[2] |= MG_NOSURF;
          }

          continue;
        }

        ref = mesh->prism[adj].ref;
        if ( adj && ( pp->ref <= ref) )  continue;

        ia = pp->v[0+i*3];
        ib = pp->v[1+i*3];
        ic = pp->v[2+i*3];

        kt = MMG5_hashGetFace(&hash,ia,ib,ic);
        if ( kt ) {
          continue;
        }

        mesh->nt++;

        ptt = &mesh->tria[mesh->nt];
        ptt->v[0] = ia;
        ptt->v[1] = ib;
        ptt->v[2] = ic;
        mesh->point[ptt->v[0]].tag |= MG_BDY;
        mesh->point[ptt->v[1]].tag |= MG_BDY;
        mesh->point[ptt->v[2]].tag |= MG_BDY;

        /* the cc field is used to be able to recover the prism (and its face)
         * from which comes a boundary triangle (when called by packmesh =>
         * mesh->nt=0 at the beginning of the function) */
        ptt->cc = 5*k + i;
        if ( pxpr ) {
          /* useful only when saving mesh */
          if ( pxpr->tag[MMG5_iarf_pr[i][0]] )  ptt->tag[0] = pxpr->tag[MMG5_iarf_pr[i][0]];
          if ( pxpr->tag[MMG5_iarf_pr[i][1]] )  ptt->tag[1] = pxpr->tag[MMG5_iarf_pr[i][1]];
          if ( pxpr->tag[MMG5_iarf_pr[i][2]] )  ptt->tag[2] = pxpr->tag[MMG5_iarf_pr[i][2]];
        }
        if ( adj ) {
          if ( mesh->info.iso ) ptt->ref = mesh->info.isoref;
          /* useful only when saving mesh */
          else ptt->ref  = pxpr ? pxpr->ref[i] : 0;
        }
        else {
          /* useful only when saving mesh */
          ptt->ref  = pxpr ? pxpr->ref[i] : 0;
        }
      }
    }
  }

  assert(mesh->nt==ntmesh);

  if ( ntmesh != ntinit ) {
    /* set point tag */
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        ppt = &mesh->point[ptt->v[i]];
        ppt->tag |= MG_BDY;
      }
    }
  }

  if ( tofree ) MMG5_DEL_MEM(mesh,hash.item);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 otherwise.
 *
 * - Remove double triangles from tria array.
 *
 * - Remove triangles that do not belong to a boundary (non opnbdy mode) from
 *   tria array.
 *
 * - Check the matching between actual and given number of faces in the mesh:
 * Count the number of faces in mesh and compare this number to the number of
 *   given triangles.
 *
 * - If the founded number exceed the given one, add the missing
 *   boundary triangles (call to MMG5_bdryTria). Do nothing otherwise.
 *
 * - Fill the adjacency relationship between prisms and tetra (fill adjapr with
 *   a negative value to mark this special faces).
 *
 * - Set to required the triangles at interface betwen prisms and tet.
 *
 */
int MMG5_chkBdryTria(MMG5_pMesh mesh) {
  MMG5_pTetra    pt,pt1;
  MMG5_pPrism    pp,pp1;
  MMG5_pTria     ptt,pttnew;
  MMG5_int       *adja,adj,k,kk,i,j,ntmesh;
  MMG5_int       ia,ib,ic, nbl,nt,ntpres;
  int            iface;
  MMG5_Hash      hashElt, hashTri;

  /** Step 1: scan the mesh and count the boundaries */
  ntmesh = ntpres = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i];

      if ( !adj ) {
        ++ntmesh;
        continue;
      }
      adj /= 4;
      pt1 = &mesh->tetra[adj];
      if ( pt->ref > pt1->ref )
        ++ntmesh;
    }
  }

  if ( mesh->info.opnbdy && mesh->xtetra ) {
    /* We want to preserve internal triangle and we came from bdryBuild: we need
     * to count the preserved boudaries */
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || !pt->xt )  continue;
      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i];

        if ( !adj ) continue;

        adj /= 4;
        pt1 = &mesh->tetra[adj];

        if ( pt->ref != pt1->ref ) continue;

        if ( (mesh->xtetra[pt->xt].ftag[i] & MG_BDY) &&
             (MG_GET(mesh->xtetra[pt->xt].ori,i) ) ) ++ntpres;
      }
    }
  }

  if ( mesh->nprism ) {
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) )  continue;

      adja = &mesh->adjapr[5*(k-1)+1];
      for (i=0; i<2; i++) {
        adj = adja[i];

        if ( !adj ) {
          ++ntmesh;
          continue;
        }
        else if ( adj<0 ) {
          continue;
        }

        adj /= 5;
        pp1 = &mesh->prism[adj];
        if ( pp->ref > pp1->ref) {
          ++ntmesh;
        }
      }
    }

    /* Detect the triangles at the interface of the prisms and tetra (they have been
     * counted twice) */
    if ( ! MMG5_hashNew(mesh,&hashTri,0.51*ntmesh,1.51*ntmesh) ) return 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i];
        if ( adj ) continue;

        ia = pt->v[MMG5_idir[i][0]];
        ib = pt->v[MMG5_idir[i][1]];
        ic = pt->v[MMG5_idir[i][2]];
        if ( !MMG5_hashFace(mesh,&hashTri,ia,ib,ic,5*k+i) ) {
          MMG5_DEL_MEM(mesh,hashTri.item);
          return 0;
        }
      }
    }

    /* Fill the adjacency relationship between prisms and tetra (fill adjapr with
     * a negative value to mark this special faces) */
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) )  continue;

      adja = &mesh->adjapr[5*(k-1)+1];
      for (i=0; i<2; i++) {
        adj = adja[i];
        if ( adj ) continue;

        ia = pp->v[MMG5_idir_pr[i][0]];
        ib = pp->v[MMG5_idir_pr[i][1]];
        ic = pp->v[MMG5_idir_pr[i][2]];

        j = MMG5_hashGetFace(&hashTri,ia,ib,ic);
        if ( !j ) continue;

        --ntmesh;
        adja[i] = -j;
      }
    }
    MMG5_DEL_MEM(mesh,hashTri.item);
  }

  /** Step 2: detect the extra boundaries (that will be ignored) provided by the
   * user */
  if ( mesh->nt ) {
    if ( ! MMG5_hashNew(mesh,&hashElt,0.51*ntmesh,1.51*ntmesh) ) return 0;
    // Hash the boundaries found in the mesh
    if ( mesh->info.opnbdy) {
      /* We want to keep the internal triangles: we must hash all the tetra faces */
      for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !MG_EOK(pt) )  continue;

        for (i=0; i<4; i++) {
          ia = pt->v[MMG5_idir[i][0]];
          ib = pt->v[MMG5_idir[i][1]];
          ic = pt->v[MMG5_idir[i][2]];
          if ( !MMG5_hashFace(mesh,&hashElt,ia,ib,ic,4*k+i) ) return 0;
        }
      }
    } else {
      for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !MG_EOK(pt) )  continue;
        adja = &mesh->adja[4*(k-1)+1];
        for (i=0; i<4; i++) {
          adj = adja[i];
          if ( !adj ) {
            ia = pt->v[MMG5_idir[i][0]];
            ib = pt->v[MMG5_idir[i][1]];
            ic = pt->v[MMG5_idir[i][2]];
            if ( !MMG5_hashFace(mesh,&hashElt,ia,ib,ic,4*k+i) ) return 0;
          }
          adj /= 4;

          pt1 = &mesh->tetra[adj];
          if ( pt->ref > pt1->ref ) {
            ia = pt->v[MMG5_idir[i][0]];
            ib = pt->v[MMG5_idir[i][1]];
            ic = pt->v[MMG5_idir[i][2]];
            if ( !MMG5_hashFace(mesh,&hashElt,ia,ib,ic,4*k+i) ) return 0;
          }
        }
      }
    }
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) )  continue;
      adja = &mesh->adjapr[5*(k-1)+1];
      for (i=0; i<2; i++) {
        adj = adja[i];
        if ( !adj ) {
          ia = pp->v[MMG5_idir_pr[i][0]];
          ib = pp->v[MMG5_idir_pr[i][1]];
          ic = pp->v[MMG5_idir_pr[i][2]];
          if ( !MMG5_hashFace(mesh,&hashElt,ia,ib,ic,5*k+i) ) return 0;
        }
        else if ( adj<0 ) continue;

        adj /= 5;

        pp1 = &mesh->prism[MMG5_abs(adj)];
        if ( pp->ref > pp1->ref ) {
          ia = pp->v[MMG5_idir_pr[i][0]];
          ib = pp->v[MMG5_idir_pr[i][1]];
          ic = pp->v[MMG5_idir_pr[i][2]];
          if ( !MMG5_hashFace(mesh,&hashElt,ia,ib,ic,5*k+i) ) return 0;
        }
      }
    }


    // Travel through the tria, delete those that are not in the hash tab or
    // that are stored more that once.
    nt=0; nbl=1;

    if ( ! MMG5_hashNew(mesh,&hashTri,0.51*mesh->nt,1.51*mesh->nt) ) return 0;

    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];

      ia = ptt->v[0];
      ib = ptt->v[1];
      ic = ptt->v[2];

      i = MMG5_hashGetFace(&hashElt,ia,ib,ic);
      j = MMG5_hashFace(mesh,&hashTri,ia,ib,ic,k);

      ptt->cc = i;

      if ( !j ) {
        MMG5_DEL_MEM(mesh,hashElt.item);
        MMG5_DEL_MEM(mesh,hashTri.item);
        return 0;
      }
      else if ( j > 0 ) {
        /* the face already exists in the tria table */
        continue;
      }

      if ( !i ) {
        /* the triangle is not a boundary tri or a tri at the interface of two
         * subdomains with different references and the user don't ask to keep
         * it. */
        continue;
      }

      if ( mesh->info.opnbdy ) {
        kk    = i/4;
        iface = i%4;
        adj   = mesh->adja[4*(kk-1)+1+iface];
        /* Check if we have found a triangle at the interface of 2 doms of same
         * ref */
        if ( adj && mesh->tetra[kk].ref == mesh->tetra[adj/4].ref ) {
          ++ntpres;
        }
      }

      ++nt;
      if ( k!=nbl ) {
        pttnew = &mesh->tria[nbl];
        memcpy(pttnew,ptt,sizeof(MMG5_Tria));
      }
      ++nbl;
    }
    nbl = mesh->nt-nt;
    if ( nbl ) {
      fprintf(stderr,"\n  ## Warning: %s: %" MMG5_PRId " extra boundaries provided."
              " Ignored\n",__func__,nbl);
      MMG5_ADD_MEM(mesh,(-nbl)*sizeof(MMG5_Tria),"triangles",return 0);
      MMG5_SAFE_REALLOC(mesh->tria,mesh->nt+1,nt+1,MMG5_Tria,"triangles",return 0);
      mesh->nt = nt;
    }
    MMG5_DEL_MEM(mesh,hashElt.item);
    MMG5_DEL_MEM(mesh,hashTri.item);
  }
  ntmesh +=ntpres;

  /** Step 3: add the missing boundary triangles or, if the mesh contains
   * prisms, set to required the triangles at interface betwen prisms and tet */
  if ( ntpres && (mesh->info.imprim > 5 || mesh->info.ddebug) )
    printf("     %" MMG5_PRId " triangles between 2 tetrahdra with same"
           " references\n",ntpres);

  if ( mesh->nt==ntmesh && !mesh->nprism ) {
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      for (i=0; i<3; i++) {
        if ( !mesh->info.iso ) mesh->point[ptt->v[i]].tag |= MG_BDY;
      }
    }
    return 1;
  }

  nbl = 0;
  if ( !mesh->nt ) {
    MMG5_ADD_MEM(mesh,(ntmesh+1)*sizeof(MMG5_Tria),"triangles",return 0);
    MMG5_SAFE_CALLOC(mesh->tria,ntmesh+1,MMG5_Tria,return 0);
  }
  else {
    assert((!mesh->nprism && ntmesh>mesh->nt)||(mesh->nprism && ntmesh>=mesh->nt));
    if ( ntmesh > mesh->nt ) {
      MMG5_ADD_MEM(mesh,(ntmesh-mesh->nt)*sizeof(MMG5_Tria),"triangles",return 0);
      MMG5_SAFE_RECALLOC(mesh->tria,mesh->nt+1,ntmesh+1,MMG5_Tria,"triangles",return 0);
      nbl = ntmesh-mesh->nt;
    }
  }
  if ( nbl && (mesh->info.imprim > 5 || mesh->info.ddebug) )
    fprintf(stderr,"\n  ## Warning: %s: %" MMG5_PRId " extra boundaries founded\n",
            __func__,nbl);

  /* Fill missing bdy triangles */
  return MMG5_bdryTria(mesh,ntmesh);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 if success.
 *
 * Set the triangles references to the tetrahedra faces and edges.
 *
 */
int MMG5_bdrySet(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  MMG5_pPrism   pp;
  MMG5_pTria    ptt;
  MMG5_pxTetra  pxt;
  MMG5_pxPrism  pxp;
  MMG5_Hash     hash;
  MMG5_int      ref,*adja,adj,k,ia,ib,ic,kt,initedg[3];
  int           j;
  int16_t       tag,inittag[3];
  int8_t        i,i1,i2;

  if ( !mesh->nt )  return 1;

  if ( mesh->xtetra ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Error: %s: mesh->xtetra must be freed.\n",__func__);
    }
    return 0;
  }
  if ( mesh->xprism ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Error: %s: mesh->xprism must be freed.\n",__func__);
    }
    return 0;
  }

  if ( ! MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return 0;
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) ) return 0;
  }

  mesh->xt     = 0;
  mesh->xtmax  = mesh->ntmax;
  assert(mesh->xtmax);

  MMG5_ADD_MEM(mesh,(mesh->xtmax+1)*sizeof(MMG5_xTetra),"boundary tetrahedra",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->xtetra,mesh->xtmax+1,MMG5_xTetra,return 0);

  /* assign references to tetras faces */
  if ( !mesh->info.opnbdy ) {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i] / 4;
        pt1 = &mesh->tetra[adj];
        if ( !adj || ( pt->ref != pt1->ref) ) {
          ia = pt->v[MMG5_idir[i][0]];
          ib = pt->v[MMG5_idir[i][1]];
          ic = pt->v[MMG5_idir[i][2]];
          kt = MMG5_hashGetFace(&hash,ia,ib,ic);
          assert(kt);
          if ( !pt->xt ) {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");return 0;);
            }
            pt->xt = mesh->xt;
          }
          ptt = &mesh->tria[kt];
          pxt = &mesh->xtetra[pt->xt];
          pxt->ref[i]   = ptt->ref;
          pxt->ftag[i] |= MG_BDY;
          pxt->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);
        }
      }
    }
  }
  else {
    /* Internal triangles preservations */
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<4; i++) {
        ia = pt->v[MMG5_idir[i][0]];
        ib = pt->v[MMG5_idir[i][1]];
        ic = pt->v[MMG5_idir[i][2]];
        kt = MMG5_hashGetFace(&hash,ia,ib,ic);

        if ( !kt ) continue;

        if ( !pt->xt ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");return 0;);
          }
          pt->xt = mesh->xt;
        }
        ptt = &mesh->tria[kt];
        pxt = &mesh->xtetra[mesh->xt];
        pxt->ref[i]   = ptt->ref;
        pxt->ftag[i] |= MG_BDY;
        pxt->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);
      }
    }
  }

  if ( !mesh->info.opnbdy ) {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      if ( !pt->xt )  continue;
      pxt = &mesh->xtetra[pt->xt];
      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i] / 4;
        pt1 = &mesh->tetra[adj];
        /* Set edge tag */
        if ( pxt->ftag[i] ) {
          if ( adj && (pt->ref == pt1->ref ) ) {
            continue;
          }
          else {
            ia = pt->v[MMG5_idir[i][0]];
            ib = pt->v[MMG5_idir[i][1]];
            ic = pt->v[MMG5_idir[i][2]];
            kt = MMG5_hashGetFace(&hash,ia,ib,ic);
            ptt = &mesh->tria[kt];

            /* Set flag to know if tetra has the same orientation than the
             * triangle */
            if ( ptt->v[0] == ia && ptt->v[1] == ib && ptt->v[2] == ic ) {
              MG_SET(pxt->ori,i);
              for (j=0; j<3; j++) {
                tag = pxt->ftag[i] | ptt->tag[j];
                if ( tag ) {
                  if ( !MMG5_settag(mesh,k,MMG5_iarf[i][j],tag,ptt->edg[j]) )
                    return 0;
                }
              }
            }
            else
              MG_CLR(pxt->ori,i);
          }
        }
      }
    }
  }
  else {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      if ( !pt->xt )  continue;
      pxt = &mesh->xtetra[pt->xt];
      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i] / 4;
        pt1 = &mesh->tetra[adj];

        ia = pt->v[MMG5_idir[i][0]];
        ib = pt->v[MMG5_idir[i][1]];
        ic = pt->v[MMG5_idir[i][2]];
        kt = MMG5_hashGetFace(&hash,ia,ib,ic);

        if ( !kt ) continue;

        ptt = &mesh->tria[kt];

        /* Set flag to know if tetra has the same orientation than the triangle
         * + force the triangle numbering to match the tetra face numbering */
        if ( adj ) {
          for ( j=0; j<3; ++j ) {
            i1 = MMG5_inxt2[j];
            i2 = MMG5_inxt2[i1];
            if (  ptt->v[j]==ia && ptt->v[i1]==ib && ptt->v[i2]==ic )
              break;
          }
          if ( j<3 ) {
            MG_SET(pxt->ori,i);
            if ( j!=0 ) {
              /* Triangle vertices+tag/edg reordering */
              ptt->v[0] = ia;
              ptt->v[1] = ib;
              ptt->v[2] = ic;

              inittag[0]  = ptt->tag[0];
              inittag[1]  = ptt->tag[1];
              inittag[2]  = ptt->tag[2];
              ptt->tag[0] = inittag[j];
              ptt->tag[1] = inittag[i1];
              ptt->tag[2] = inittag[i2];

              initedg[0]  = ptt->edg[0];
              initedg[1]  = ptt->edg[1];
              initedg[2]  = ptt->edg[2];
              ptt->edg[0] = initedg[j];
              ptt->edg[1] = initedg[i1];
              ptt->edg[2] = initedg[i2];
            }
          }
          else {
            MG_CLR(pxt->ori,i);
          }
        }
        else  MG_SET(pxt->ori,i);

        /* Set edge tag */
        if ( pxt->ftag[i] ) {
          if ( adj && ( (pt->ref < pt1->ref) || !MG_GET(pxt->ori,i) ) ) {
            continue;
          }
          else {
            for (j=0; j<3; j++) {
              tag = pxt->ftag[i] | ptt->tag[j];
              if ( tag ) {
                if ( !MMG5_settag(mesh,k,MMG5_iarf[i][j],tag,ptt->edg[j]) )
                  return 0;
              }
            }
          }
        }
      }
    }
  }

  if ( !mesh->nprism ) {
    MMG5_DEL_MEM(mesh,hash.item);
    return 1;
  }

  mesh->xpr     = 0;
  MMG5_ADD_MEM(mesh,(mesh->nprism+1)*sizeof(MMG5_xPrism),"boundary prisms",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->xprism,mesh->nprism+1,MMG5_xPrism,return 0);

  /* assign references to prism faces */
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) )  continue;
    adja = &mesh->adjapr[5*(k-1)+1];
    for (i=0; i<2; i++) {
      adj = adja[i] / 5;
      if ( adj<0 ) {
        ref = mesh->tetra[MMG5_abs(adj)].ref;
      } else {
        ref = mesh->prism[adj].ref;
      }
      if ( adj && (pp->ref == ref) ) continue;

      ia = pp->v[MMG5_idir_pr[i][0]];
      ib = pp->v[MMG5_idir_pr[i][1]];
      ic = pp->v[MMG5_idir_pr[i][2]];
      kt = MMG5_hashGetFace(&hash,ia,ib,ic);
      assert(kt);
      if ( !pp->xpr ) {
        mesh->xpr++;
        pp->xpr = mesh->xpr;
      }
      ptt = &mesh->tria[kt];
      pxp = &mesh->xprism[mesh->xpr];
      pxp->ref[i]   = ptt->ref;
      pxp->ftag[i] |= MG_BDY;
      pxp->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);

      for (j=0; j<3; j++) {
        pxp->tag[MMG5_iarf[i][j]] |= pxp->ftag[i] | ptt->tag[j];
        pxp->edg[MMG5_iarf[i][j]] = ptt->edg[j];
      }
    }
  }
  MMG5_ADD_MEM(mesh,(mesh->xpr-mesh->nprism)*sizeof(MMG5_xPrism),"boundary prisms",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_RECALLOC(mesh->xprism,mesh->nprism+1,mesh->xpr+1,MMG5_xPrism,
                      "boundary prisms",return 0);

  MMG5_DEL_MEM(mesh,hash.item);
  return 1;
}

/**
 * \param mesh
 * \return 1 if successful, 0 if not.
 *
 * Update tag and refs of tetra edges. If tetra is required, set the
 * faces/edges to required.
 *
 * \remark While remeshing:
 *   - A tetra with a boundary face has a xtetra:
 *     - a boundary edge has consistent tags as soon as it has a non 0 tag;
#warning check this case
 *     - some boundary edges may have a 0 tag
 *
 * - a tetra with a boundary edge but no bdy faces may or may not have a xtetra
 *     and boundary edges may have or may not have tags so we can't guess
 *     nothing with such xtetra.
 */
int MMG5_bdryUpdate(MMG5_pMesh mesh) {
  MMG5_pTetra   pt;
  MMG5_pTria    ptt;
  MMG5_pxTetra  pxt;
  MMG5_Hash     hash;
  MMG5_int      ia,ib,ic,k,kt;
  int           j;
  int16_t       tag;
  int8_t        i;

  if ( !mesh->nt )  return 1;
  if ( !MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) )  return 0;
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( pt->tag & MG_REQ ) {
      mesh->point[mesh->tetra[k].v[0]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[1]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[2]].tag |= MG_REQ;
      mesh->point[mesh->tetra[k].v[3]].tag |= MG_REQ;
      if ( !MMG5_settag(mesh,k,0,MG_REQ,0) ) return 0;
      if ( !MMG5_settag(mesh,k,1,MG_REQ,0) ) return 0;
      if ( !MMG5_settag(mesh,k,2,MG_REQ,0) ) return 0;
      if ( !MMG5_settag(mesh,k,3,MG_REQ,0) ) return 0;
      if ( !MMG5_settag(mesh,k,4,MG_REQ,0) ) return 0;
      if ( !MMG5_settag(mesh,k,5,MG_REQ,0) ) return 0;
      mesh->point[mesh->tetra[k].v[0]].tag &= ~MG_NOSURF;
      mesh->point[mesh->tetra[k].v[1]].tag &= ~MG_NOSURF;
      mesh->point[mesh->tetra[k].v[2]].tag &= ~MG_NOSURF;
      mesh->point[mesh->tetra[k].v[3]].tag &= ~MG_NOSURF;
      if ( !MMG5_deltag(mesh,k,0,MG_NOSURF) ) return 0;
      if ( !MMG5_deltag(mesh,k,1,MG_NOSURF) ) return 0;
      if ( !MMG5_deltag(mesh,k,2,MG_NOSURF) ) return 0;
      if ( !MMG5_deltag(mesh,k,3,MG_NOSURF) ) return 0;
      if ( !MMG5_deltag(mesh,k,4,MG_NOSURF) ) return 0;
      if ( !MMG5_deltag(mesh,k,5,MG_NOSURF) ) return 0;
    }

    if ( !pt->xt )  continue;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++) {
      /* Set edge tag */
      if ( ! MG_GET(pxt->ori,i) ) continue;
      if ( pxt->ftag[i] & MG_BDY ) {
        ia = pt->v[MMG5_idir[i][0]];
        ib = pt->v[MMG5_idir[i][1]];
        ic = pt->v[MMG5_idir[i][2]];
        kt = MMG5_hashGetFace(&hash,ia,ib,ic);
        assert(kt);
        ptt = &mesh->tria[kt];
        if ( pt->tag & MG_REQ ) {
          pxt->ftag[i] |= MG_REQ;
          ptt->tag[0]   = MG_REQ;
          ptt->tag[1]   = MG_REQ;
          ptt->tag[2]   = MG_REQ;
          pxt->ftag[i] &= ~MG_NOSURF;
          ptt->tag[0]  &= ~MG_NOSURF;
          ptt->tag[1]  &= ~MG_NOSURF;
          ptt->tag[2]  &= ~MG_NOSURF;
        }
        for ( j=0; j<3; j++ ) {
          tag = ptt->tag[j];
          if ( tag || ptt->edg[j] ) {
            if ( !MMG5_settag(mesh,k,MMG5_iarf[i][j],tag,ptt->edg[j]) )
              return 0;
          }
        }
      }
    }
  }
  MMG5_DEL_MEM(mesh,hash.item);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Make orientation of triangles compatible with tetra faces for external tria
 * and with domain of max ref for interface tria.
 *
 */
int MMG5_bdryPerm(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  MMG5_pTria    ptt;
  MMG5_Hash     hash;
  MMG5_int      *adja,adj,k,kt,ia,ib,ic,nf;
  int8_t        i;

  if ( !mesh->nt ) return 1;

  /* store triangles temporarily */
  if ( !MMG5_hashNew(mesh,&hash,MG_MAX(0.51*mesh->nt,100),MG_MAX(1.51*mesh->nt,300)) )
    return 0;

  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
  }

  /* check orientation */
  nf = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      adj = adja[i] / 4;
      pt1 = &mesh->tetra[adj];

      if ( adj && (pt->ref <= pt1->ref || pt->ref == MG_PLUS) )
        continue;

      ia = pt->v[MMG5_idir[i][0]];
      ib = pt->v[MMG5_idir[i][1]];
      ic = pt->v[MMG5_idir[i][2]];
      kt = MMG5_hashGetFace(&hash,ia,ib,ic);
      if ( kt ) {
        /* check orientation */
        ptt = &mesh->tria[kt];
        if ( ptt->v[0] == ia && ptt->v[1] == ib && ptt->v[2] == ic )
          continue;
        else {
          ptt->v[0] = ia;
          ptt->v[1] = ib;
          ptt->v[2] = ic;
          nf++;
        }
      }
    }
  }
  if ( mesh->info.ddebug && nf > 0 )
    fprintf(stdout,"  ## %" MMG5_PRId " faces reoriented\n",nf);

  MMG5_DEL_MEM(mesh,hash.item);

  return 1;
}
