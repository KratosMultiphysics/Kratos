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
#include "libmmg2d_private.h"

#define KTA     7
#define KTB    11

/**
 * \param mesh pointer toward the mesh
 * \return 1 if success, 0 if fail
 *
 * Create adjacency relations between the triangles dein the mesh
 *
 */
int MMG2D_hashTria(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  MMG5_int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1;
  MMG5_int       *hcode,*link,inival,hsize,iadr;
  uint8_t        i,ii,i1,i2;
  unsigned int   key;

  if ( mesh->adja )  return 1;
  if ( !mesh->nt )  return 0;

  /* memory alloc */
  MMG5_SAFE_CALLOC(hcode,mesh->nt+1,MMG5_int,return 0);

  /* memory alloc */
  MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(MMG5_int),"adjacency table",
                printf("  Exit program.\n");
                return 0;);
  MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,MMG5_int,return 0);

  link  = mesh->adja;
  hsize = mesh->nt;

  /* init */
  if ( sizeof(MMG5_int) == 8 ) {
    inival = LONG_MAX;
  }
  else {
    inival = INT_MAX;
  }

  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !pt->v[0] )  continue;
    for (i=0; i<3; i++) {
      i1 = MMG2D_idir[i+1];
      i2 = MMG2D_idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }

      /* compute key */
      key = (KTA*(int64_t)mins + KTB*(int64_t)maxs)%hsize+1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = MMG2D_idir[i+1];
    i2 = MMG2D_idir[i+2];
    pt = &mesh->tria[k];

    mins = M_MIN(pt->v[i1],pt->v[i2]);
    maxs = M_MAX(pt->v[i1],pt->v[i2]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;

    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = MMG2D_idir[ii+1];
      i2 = MMG2D_idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }

      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
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
 * Create full table of adjacency for quadrangles (quad \f$ <-> \f$ quad
 * adjacencies and quad \f$ -> \f$ tri adjacencies): 1) if the edge \f$ i1 \f$
 * of quad \f$ k1 \f$ is adja to quad \f$ k2 \f$ through edge \f$ i2 \f$, \f$
 * adja[4*(k1-1)+1+i1] = 4*k2+i2 \f$.  2) if the edge \f$ i1 \f$ of quad \f$ k1
 * \f$ is adja to tria \f$ k2 \f$ through edge \f$ i2 \f$, \f$
 * adja[4*(k1-1)+1+i1] = -(3*k2+i2) \f$.
 *
 */
int MMG2D_hashQuad(MMG5_pMesh mesh) {
  MMG5_pQuad     pq,pq1;
  MMG5_pTria     pt;
  MMG5_Hash      hash;
  MMG5_int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1,iadr;
  MMG5_int       *hcode,*link,hsize,inival;
  uint8_t        i,ii,i1,i2;
  unsigned int   key;

  /** Step 1: Fill adjacendies between quadrangles */
  if ( !mesh->nquad ) {
    return 1;
  }

  /* default */
  if ( mesh->adjq ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: no re-build of adjacencies of quadrangles. "
              "mesh->adjq must be freed to enforce analysis.\n",__func__);
    }
    return 1;
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING QUAD ADJACENCY\n");

  /* memory alloc */
  MMG5_ADD_MEM(mesh,(4*mesh->nquad+5)*sizeof(MMG5_int),"quad adjacency table",
               fprintf(stderr,"  Exit program.\n");
               return 0);
  MMG5_SAFE_CALLOC(mesh->adjq,4*mesh->nquad+5,MMG5_int,return 0);
  MMG5_SAFE_CALLOC(hcode,mesh->nquad+5,MMG5_int,return 0);

  link  = mesh->adjq;
  hsize = mesh->nquad;

  /* init */
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");

  if ( sizeof(MMG5_int) == 8 ) {
    inival = LONG_MAX;
  }
  else {
    inival = INT_MAX;
  }

  iadr   = 0;
  for (k=0; k<=mesh->nquad; k++)
    hcode[k] = -inival;

  /* hash quads */
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;
    for (i=0; i<4; i++) {
      i1 = MMG2D_idir_q[i][0];
      i2 = MMG2D_idir_q[i][1];
      if ( pq->v[i1] < pq->v[i2] ) {
        mins = pq->v[i1];
        maxs = pq->v[i2];
      }
      else {
        mins = pq->v[i2];
        maxs = pq->v[i1];
      }

      /* compute key */
      key = (KTA*(int64_t)mins + KTB*(int64_t)maxs)%hsize+1;

      /* insert */
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
    i1 = MMG2D_idir_q[i][0];
    i2 = MMG2D_idir_q[i][1];

    pq = &mesh->quadra[k];
    if ( pq->v[i1] < pq->v[i2] ) {
      mins = pq->v[i1];
      maxs = pq->v[i2];
    }
    else {
      mins = pq->v[i2];
      maxs = pq->v[i1];
    }

    /* accross link */
    ll      = -link[l];
    pp      = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 4 + 1;
      ii = (ll-1) % 4;
      i1 = MMG2D_idir_q[ii][0];
      i2 = MMG2D_idir_q[ii][1];
      pq1  = &mesh->quadra[kk];

      if ( pq1->v[i1] < pq1->v[i2] ) {
        mins1 = pq1->v[i1];
        maxs1 = pq1->v[i2];
      }
      else {
        mins1 = pq1->v[i2];
        maxs1 = pq1->v[i1];
      }

      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l]  = 4*kk + ii;
        link[ll] = 4*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  MMG5_SAFE_FREE(hcode);

  /** Step 2: Fill adjacencies between quadrangles and triangles */
  /* Temporarily allocate a hash structure for storing edges */
  if ( !MMG5_hashNew( mesh,&hash,0.51*mesh->nt, 1.51*mesh->nt) ) return 0;

  /* Hash edge belonging to only one triangle (hcode==0) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    hcode = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( hcode[i] ) {
        continue;
      }
      /* The edge belongs to one tria only */
      MMG5_hashEdge(mesh,&hash,pt->v[MMG5_iprv2[i]],pt->v[MMG5_inxt2[i]],3*k+i);
    }
  }

  /* Search edges that belong to one quad only and update adjacency array if the
   * edge is founded in the hash table */
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;

    hcode = &mesh->adjq[4*(k-1)+1];
    for (i=0; i<4; i++) {
      assert ( hcode[i] >= 0 );

      if ( hcode[i] ) {
        continue;
      }

      /* The edge belongs to one quad only */
      i1 = MMG2D_idir_q[i][0];
      i2 = MMG2D_idir_q[i][1];

      kk = MMG5_hashGet(&hash,pq->v[i1],pq->v[i2]);
      if ( kk ) {
        /* The edge is at the interface between a quad and a triangle */
        hcode[i] = -kk;
      }
    }
  }
  MMG5_DEL_MEM(mesh,hash.item);

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 0 if fail, 1 otherwise
 *
 * Transfer some input edge data to the corresponding triangles fields
 *
 */
int MMG2D_assignEdge(MMG5_pMesh mesh) {
  MMG5_Hash       hash;
  MMG5_pTria      pt;
  MMG5_pQuad      pq;
  MMG5_pEdge      pa;
  MMG5_int        ia;
  MMG5_int        k;
  int8_t          ier;
  uint8_t         i,i1,i2;

  /**
      The cleaning of required tags inside triangles has been initially added by
      commit da4b099c. It probably followed the report of a bug arising when
      several library functions are successively called without cleaning the
      mesh structure but on december 2022 I am not able to reproduce this bug.
      Due to tags cleaning, input required edges are lost in ls discretization
      mode (see issue #171).
      The metRidtyp field (previsouly not used in 2D) is
      now used to mark if \a MMG2D_assignEdge function is called for the first
      time inside the library and if we have to clean triangle tags (in order to
      fix issue #171 without breaking again the initial fix). */
  if ( !mesh->info.metRidTyp ) {
    /* Try to clean triangle structure (in case where mmg2dlib is called after
     * mmg2dmesh) */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      /* If all edges are required, the triangle is maybe required by the user */
      if ( (pt->tag[0] & MG_REQ) && (pt->tag[1] & MG_REQ) && (pt->tag[2] & MG_REQ) ) {
        continue;
      }

      /* Otherwise there is no reason to have required tags on edges */
      for ( ia = 0; ia < 3; ++ia ) {
        pt->tag[ia] &= ~MG_REQ;
      }
    }
    mesh->info.metRidTyp = 1;
  }

  if ( !mesh->na ) return 1;

  /* Temporarily allocate a hash structure for storing edges */
  ier = MMG5_hashNew ( mesh,&hash, mesh->na,3*mesh->na );
  if ( !ier ) {
    printf("  ## Error: %s: Unable to allocate edge hash table\n.",__func__);
    return 0;
  }

  /* hash mesh edges */
  for (k=1; k<=mesh->na; k++) {
    ier = MMG5_hashEdge(mesh,&hash,mesh->edge[k].a,mesh->edge[k].b,k);
    if ( !ier ) {
      fprintf(stderr,"\n  ## Error: %s: unable to hash edge %" MMG5_PRId " %" MMG5_PRId ".\n",__func__,
              MMG2D_indPt(mesh,mesh->edge[k].a),MMG2D_indPt(mesh,mesh->edge[k].b));
      return 0;
    }
  }

  /* set references to triangles */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      ia = MMG5_hashGet(&hash,pt->v[i],pt->v[i1]);
      if ( ia ) {
        i2 = MMG5_inxt2[i1];
        pa = &mesh->edge[ia];
        pt->edg[i2] = pa->ref;
        pt->tag[i2] |= pa->tag;

      }
    }
  }

  /* set references to quadrangles */
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;

    for (i=0; i<4; i++) {
      i1 = MMG2D_idir_q[i][0];
      i2 = MMG2D_idir_q[i][1];
      ia = MMG5_hashGet(&hash,pq->v[i1],pq->v[i2]);
      if ( ia ) {
        pa = &mesh->edge[ia];
        pq->edg[i]  = pa->ref;
        pq->tag[i] |= pa->tag;
      }
    }
  }

  /* Delete the hash for edges */
  MMG5_DEL_MEM(mesh,hash.item);
  MMG5_DEL_MEM(mesh,mesh->edge);
  mesh->na = 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 1 if success, 0 if fail
 *
 * Create the edges in the mesh from the information stored in the triangles, or
 * by identifying the different components of the mesh.
 *
 * \remark Possible extension needed to take into account constrained edges
 * \remark Call in debug mode only
 *
 */
int MMG2D_bdryEdge(MMG5_pMesh mesh) {
  MMG5_pTria      pt,pt1;
  MMG5_pEdge      pa;
  MMG5_pPoint     p0;
  MMG5_int        k,*adja,natmp,iel;
  int8_t          i,i1,i2;

  natmp = 0;
  mesh->na = 0;

  /* First step: Count number of boundary edges */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      iel = adja[i] / 3;
      pt1 = &mesh->tria[iel];

      if ( iel && pt->ref <= pt1->ref ) continue;
      natmp++;
    }
  }

  /* Second step: Create edge mesh and store the corresponding edges */
  MMG5_ADD_MEM(mesh,(natmp+1)*sizeof(MMG5_Edge),"edges",return 0);
  MMG5_SAFE_CALLOC(mesh->edge,natmp+1,MMG5_Edge,return 0);

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      iel = adja[i] / 3;
      pt1 = &mesh->tria[iel];

      if ( iel && pt->ref <= pt1->ref ) continue;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];

      mesh->na++;
      pa = &mesh->edge[mesh->na];
      pa->a = pt->v[i1];
      pa->b = pt->v[i2];

      pa->tag = pt->tag[i];
      pa->ref = pt->edg[i];

    }
  }

  /* Set point tags */
  for (k=1; k<=mesh->na; k++) {
    pa = &mesh->edge[k];
    p0 = &mesh->point[pa->a];
    p0->tag |= MG_BDY;

    p0 = &mesh->point[pa->b];
    p0->tag |= MG_BDY;
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a solution structure.
 * \param met pointer toward a solution structure.
 *
 * \return 0 if memory problem (uncomplete mesh), 1 otherwise.
 *
 * Pack the mesh and metric and create explicitly all the mesh structures
 * (edges).
 *
 */
int MMG2D_pack(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  MMG5_pTria         pt,ptnew,pt1;
  MMG5_pQuad         pq,pq1;
  MMG5_pEdge         ped;
  MMG5_pPoint        ppt,pptnew;
  MMG5_int           np,ned,nt,k,iel,nbl,isol,isolnew,memWarn,nc;
  MMG5_int           iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;
  int8_t             i,i1,i2;

  /* Keep only one domain if asked */
  MMG2D_keep_only1Subdomain ( mesh, mesh->info.nsd );

  /** Recreate adjacencies if need be */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Warning: %s: hashing problem. Exit program.\n",
            __func__);
    return 0;
  }

  /** Pack vertex indices */
  np = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++np;

    if ( ppt->tag & MG_CRN )  nc++;

    if ( ppt->tag & MG_NOSURF ) {
      ppt->tag &= ~MG_NOSURF;
      ppt->tag &= ~MG_REQ;
    }
  }

  /** Count the number of edges in the mesh */
  memWarn = 0;
  ned = 0;

  if ( mesh->edge ) {
    fprintf(stderr,"\n  ## Warning: %s: unexpected edge table..."
            " Ignored data.\n",__func__);
    MMG5_DEL_MEM(mesh,mesh->edge);
    mesh->na = 0;
  }

  if ( mesh->nquad && mesh->quadra ) {
    if ( !MMG2D_hashQuad(mesh) ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to build quad adjacencies."
              " Quad edges will be ignored.\n",__func__);
    }
  }

  mesh->na = 0;
  /** Count edges stored in triangles */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      iel = adja[i] / 3;

      if ( pt->tag[i] & MG_NOSURF ) {
        pt->tag[i] &= ~MG_REQ;
        pt->tag[i] &= ~MG_NOSURF;
      }

      pt1 = &mesh->tria[iel];
      if ( (!iel) || (pt->ref > pt1->ref) ) {
        ++mesh->na;
      }
      else if ( (pt->ref==pt1->ref) && MG_SIN(pt->tag[i]) ) {
        ++mesh->na;
      }
      else if ( mesh->info.opnbdy ) {
        if ( (pt->tag[i] & MG_REF) || pt->tag[i] & MG_BDY ) {
          assert ( pt->tag[i] & (MG_REF+MG_BDY) );
          ++mesh->na;
        }
      }
    }
  }
  /** Count edges stored in quadrangles */
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) ) continue;
    adja = &mesh->adjq[4*(k-1)+1];

    for (i=0; i<4; i++) {
      iel = adja[i] / 4;

      if ( iel < 0) {
        /* Edge at the interface between a quad and a tria: treated from the tria */
        continue;
      }

      pq1 = &mesh->quadra[iel];
      if ( (!iel) || (pq->ref > pq1->ref) ) {
        ++mesh->na;
      }
      else if ( (pq->ref==pq1->ref) && MG_SIN(pq->tag[i]) ) {
        ++mesh->na;
      }
      else if ( mesh->info.opnbdy ) {
        if ( (pq->tag[i] & MG_REF) || pq->tag[i] & MG_BDY ) {
          assert ( pq->tag[i] & (MG_REF+MG_BDY) );
          ++mesh->na;
        }
      }
    }
  }

  /** Pack edges */
  mesh->namax = mesh->na;
  if ( mesh->na ) {

    MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"final edges", memWarn=1);

    if ( memWarn ) {
      if ( mesh->info.ddebug )
        printf("  -- Attempt to allocate a smallest edge table...\n");
      mesh->namax = mesh->na;
      memWarn = 0;
      MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"final edges",
                    fprintf(stderr,"\n  ## Warning: %s: uncomplete mesh.\n",
                            __func__);
                    memWarn=1);
    }

    if ( memWarn )
      mesh->na = 0;
    else {
      /* We have enough memory to allocate the edge table */
      MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge, return 0);

      nt = 0;
      /* Edges stored in triangles */
      for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MG_EOK(pt) ) continue;
        ++nt;
        adja = &mesh->adja[3*(k-1)+1];

        for (i=0; i<3; i++) {
          i1 = MMG5_inxt2[i];
          i2 = MMG5_iprv2[i];
          iel = adja[i] / 3;
          pt1 = &mesh->tria[iel];
          if ( !iel || (pt->ref > pt1->ref) ||
               ((pt->ref==pt1->ref) && MG_SIN(pt->tag[i])) ||
               (mesh->info.opnbdy && ((pt->tag[i] & MG_REF) || (pt->tag[i] & MG_BDY)))) {
            ++ned;
            ped = &mesh->edge[ned];
            ped->a = pt->v[i1];
            ped->b = pt->v[i2];
            /* the base field is used to be able to recover the tria (and its face)
             * from which comes a boundary edge */
            ped->base = 3*nt+i;
            ped->ref = pt->edg[i];
            ped->tag = pt->tag[i];
          }
        }
      }

      /* Edges stored in quadrangles */
      nt = 0;
      for (k=1; k<=mesh->nquad; k++) {
        pq = &mesh->quadra[k];
        if ( !MG_EOK(pq) ) continue;
        ++nt;
        adja = &mesh->adjq[4*(k-1)+1];

        for (i=0; i<4; i++) {
          i1 = MMG2D_idir_q[i][0];
          i2 = MMG2D_idir_q[i][1];

          iel = adja[i] / 4;

          if ( iel < 0) {
            /* Edge at the interface between a quad and a tria: treated from the tria */
            continue;
          }

          pq1 = &mesh->quadra[iel];
          if ( !iel || (pq->ref > pq1->ref) ||
               ((pq->ref==pq1->ref) && MG_SIN(pq->tag[i])) ||
               (mesh->info.opnbdy && ((pq->tag[i] & MG_REF) || (pq->tag[i] & MG_BDY)))) {
            ++ned;
            ped = &mesh->edge[ned];
            ped->a = pq->v[i1];
            ped->b = pq->v[i2];
            /* the base field is used to be able to recover the quad (and its face)
             * from which comes a boundary edge */
            ped->base = 4*nt+i;
            ped->ref = pq->edg[i];
            ped->tag = pq->tag[i];
          }
        }
      }
    }
  }

  for (k=1; k<=mesh->na; k++) {
    ped  = &mesh->edge[k];
    if ( !ped->a ) continue;
    ped->a = mesh->point[ped->a].tmp;
    ped->b = mesh->point[ped->b].tmp;
  }

  /** Pack triangles */
  nt  = 0;
  nbl = 1;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    nt++;

    if ( k != nbl ) {
      ptnew = &mesh->tria[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tria));

      /* Update the adjacency */
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadrnew = 3*(nbl-1) + 1;
      adjanew = &mesh->adja[iadrnew];

      for(i=0; i<3; i++) {
        adjanew[i] = adja[i];
        if ( !adja[i] ) continue;
        iadrv = 3*(adja[i]/3-1)+1;
        adjav = &mesh->adja[iadrv];
        voy = i;
        adjav[adja[i]%3] = 3*nbl + voy;
        adja[i] = 0;
      }
      memset(pt,0,sizeof(MMG5_Tria));
    }
    nbl++;
  }
  mesh->nt = nt;

  /** Pack quadrangles */
  if ( mesh->quadra ) {
    k = 1;
    do {
      pq = &mesh->quadra[k];
      if ( !MG_EOK(pq) ) {
        pq1 = &mesh->quadra[mesh->nquad];
        assert( pq && pq1 && MG_EOK(pq1) );
        memcpy(pq,pq1,sizeof(MMG5_Quad));
        --mesh->nquad;
      }
    }
    while ( ++k < mesh->nquad );
  }

  if ( mesh->quadra ) {
    for (k=1; k<=mesh->nquad; k++) {
      pq = &mesh->quadra[k];
      if ( !MG_EOK(pq) )  continue;
      pq->v[0] = mesh->point[pq->v[0]].tmp;
      pq->v[1] = mesh->point[pq->v[1]].tmp;
      pq->v[2] = mesh->point[pq->v[2]].tmp;
      pq->v[3] = mesh->point[pq->v[3]].tmp;
    }
  }

  /** Pack solutions (metric map, displacement, ...) */
  if ( sol && sol->m ) {
    nbl = 1;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      isol    = k * sol->size;
      isolnew = nbl * sol->size;

      for (i=0; i<sol->size; i++)
        sol->m[isolnew + i] = sol->m[isol + i];
      ++nbl;
    }
  }

  if ( met && met->m ) {
    nbl = 1;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      isol    = k * met->size;
      isolnew = nbl * met->size;

      for (i=0; i<met->size; i++)
        met->m[isolnew + i] = met->m[isol + i];
      ++nbl;
    }
  }

  /** Pack vertices */
  np  = 0;
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;

    if ( k != nbl ) {
      pptnew = &mesh->point[nbl];
      memcpy(pptnew,ppt,sizeof(MMG5_Point));
      ppt->tag   = 0;
      assert ( ppt->tmp == nbl );
    }
    np++;
    if ( k != nbl ) {
      ppt = &mesh->point[k];
      memset(ppt,0,sizeof(MMG5_Point));
      ppt->tag    = 0;
    }
    nbl++;
  }
  mesh->np = np;
  if ( sol && sol->m ) sol->np  = np;

  /** Reset ppt->tmp field */
  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  if(mesh->np < mesh->npmax - 3) {
    mesh->npnil = mesh->np + 1;
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;
  }
  else {
    mesh->npnil = 0;
  }

  /** Reset garbage collector */
  if ( mesh->nt < mesh->ntmax - 3 ) {
    mesh->nenil = mesh->nt + 1;
    for (k=mesh->nenil; k<mesh->ntmax-1; k++)
      mesh->tria[k].v[2] = k+1;
  }
  else {
    mesh->nenil = 0;
  }

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"     NUMBER OF VERTICES       %8" MMG5_PRId "   CORNERS %8" MMG5_PRId "\n",mesh->np,nc);
    fprintf(stdout,"     NUMBER OF TRIANGLES      %8" MMG5_PRId "\n",mesh->nt);
    if ( mesh->nquad ) {
      fprintf(stdout,"     NUMBER OF QUADRILATERALS %8" MMG5_PRId "\n",mesh->nquad);
    }
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES          %8" MMG5_PRId "\n",mesh->na);
  }

  if ( memWarn ) return 0;

  return 1;
}
