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
 * \file common/hash.c
 * \brief Functions for hash tables management and tetrahedra packing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param adjt pointer to the adjacency table of the surfacic mesh.
 * \param hash pointer to the edge hash table.
 * \param chkISO flag to say if we check ISO references
 * (so if we come from mmg3d).
 * \return 1 if success, 0 otherwise.
 *
 * Create surface adjacency
 *
 * \remark the ph->s field computation is useless in mmgs.
 *
 *
 * \remark: as all triangles are mesh boundaries, we do not need to mark their
 * adges as MG_BDY so the MG_BDY tag may be used inside geometrical triangles
 * (external non-parallel, or internal parallel) to tag edges on the
 * intersection with purely parallel (non-geometrical) triangles.
 * The MG_PARBDYBDY tag is also added, as it does not have a supporting triangle
 * to inherit this tag from.
 *
 */
int MMG5_mmgHashTria(MMG5_pMesh mesh, MMG5_int *adjt, MMG5_Hash *hash, int chkISO) {
  MMG5_pTria     pt,pt1,pt2;
  MMG5_hedge     *ph;
  MMG5_int       kk;
  MMG5_int       *adja,hmax,k,ia,ib,jel,lel,dup,nmf;
  int8_t         i,i1,i2,j,l;
  MMG5_int       key;

  /* adjust hash table params */
  hmax =(MMG5_int)(3.71*mesh->np);
  hash->siz  = mesh->np;
  hash->max  = hmax + 1;
  hash->nxt  = hash->siz;
  MMG5_ADD_MEM(mesh,(hash->max+1)*sizeof(MMG5_hedge),"hash table",return 0);
  MMG5_SAFE_CALLOC(hash->item,hash->max+1,MMG5_hedge,return 0);

  for (k=hash->siz; k<hash->max; k++)
    hash->item[k].nxt = k+1;

  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");

  /* hash triangles */
  ++mesh->base;
  dup = nmf = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    pt->flag = 0;
    // base field of triangles has to be setted because it is used in setadj (mmgs
    // and mmg3d) to detect moebius strip and to flip tria orientation
    pt->base = mesh->base;

    adja = &adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {

      if ( !MG_EOK(pt) )  continue;

      /* Skip parallel edges */
      if( (pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY) ) continue;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      /* compute key */
      ia  = MG_MIN(pt->v[i1],pt->v[i2]);
      ib  = MG_MAX(pt->v[i1],pt->v[i2]);
      key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
      ph  = &hash->item[key];

      /* store edge */
      if ( ph->a == 0 ) {
        ph->a = ia;
        ph->b = ib;
        ph->k = 3*k + i;
        ph->nxt = 0;
        ++ph->s;
        continue;
      }
      /* update info about adjacent */
#ifndef NDEBUG
      int8_t ok = 0;
#endif
      while ( ph->a ) {
        if ( ph->a == ia && ph->b == ib ) {
          jel = ph->k / 3;
          j   = ph->k % 3;
          pt1 = &mesh->tria[jel];
          /* discard duplicate face */
          if ( pt1->v[j] == pt->v[i] ) {
            pt->v[0] = 0;
            dup++;
          }
          /* update adjacent */
          else if ( !adjt[3*(jel-1)+1+j] ) {
            adja[i] = 3*jel + j;
            adjt[3*(jel-1)+1+j] = 3*k + i;
            ++ph->s;
          }
          /* non-manifold case */
          else if ( adja[i] != 3*jel+j ) {
            lel = adjt[3*(jel-1)+1+j]/3;
            l   = adjt[3*(jel-1)+1+j]%3;
            pt2 = &mesh->tria[lel];

            if ( chkISO && ( (pt->ref == mesh->info.isoref) || (pt->ref < 0)) ) {
              adjt[3*(lel-1)+1+l] = 0;
              adja[i] = 3*jel+j;
              adjt[3*(jel-1)+1+j] = 3*k + i;
            }
            pt->tag[i] |= MG_NOM;
            pt1->tag[j] |= MG_NOM;
            pt2->tag[l] |= MG_NOM;
            nmf++;
            ++ph->s;
          }
#ifndef NDEBUG
          ok = 1;
#endif
          break;
        }
        else if ( !ph->nxt ) {
          ph->nxt = hash->nxt;
          ph = &hash->item[ph->nxt];
          assert(ph);

          if ( hash->nxt >= hash->max-1 ) {
            if ( mesh->info.ddebug ) {
              fprintf(stderr,"\n  ## Warning: %s: memory alloc problem (edge):"
                      " %" MMG5_PRId "\n",__func__,hash->max);
            }
            MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,MMG5_GAP,MMG5_hedge,
                               "MMG5_edge",
                               MMG5_DEL_MEM(mesh,hash->item);
                               return 0);

            ph = &hash->item[hash->nxt];

            for (kk=ph->nxt; kk<hash->max; kk++)
              hash->item[kk].nxt = kk+1;
          }

          hash->nxt = ph->nxt;
          ph->a = ia;
          ph->b = ib;
          ph->k = 3*k + i;
          ph->nxt = 0;
          ++ph->s;
#ifndef NDEBUG
          ok = 1;
#endif
          break;
        }
        else
          ph = &hash->item[ph->nxt];
      }
      assert(ok);
    }
  }

  /* Now loop on "only" parallel edges in order to add a MG_BDY tag on those
   * that are found in the hash table and their adjacents (if manifold; in the
   * non-manifold case the MG_NOM tag suffices).
   *
   * Rationale:
   *   - put MG_BDY tags on edges that are contact edges between a "true"
   *     geometrical boundary (!MG_PARBDY || (MG_PARBDY && MG_PARBDYBDY) and a
   *     purele parallel (non-geometrical) one (MG_PARBDY && !MG_PARBDYBDY);
   *   - add also a MG_PARBDYBDY tag on those edges (as it cannot be inherited
   *     from any triangle there) and their extremities.
   */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      if( !(pt->tag[i] & MG_PARBDY) || (pt->tag[i] & MG_PARBDYBDY) ) continue;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      /* compute key */
      ia  = MG_MIN(pt->v[i1],pt->v[i2]);
      ib  = MG_MAX(pt->v[i1],pt->v[i2]);
      key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
      ph  = &hash->item[key];

      if ( ph->a == 0 )  continue;
      while( ph->a ) {
        if ( ph->a == ia && ph->b == ib ) {
          jel = ph->k / 3;
          j   = ph->k % 3;
          pt1 = &mesh->tria[jel];
          pt1->tag[j] |= MG_BDY + MG_PARBDYBDY;
          mesh->point[ia].tag |= MG_PARBDYBDY;
          mesh->point[ib].tag |= MG_PARBDYBDY;
          /* update adjacent */
          lel = adjt[3*(jel-1)+1+j]/3;
          l   = adjt[3*(jel-1)+1+j]%3;
          if( lel ) {
            pt2 = &mesh->tria[lel];
            pt2->tag[l] |= MG_BDY + MG_PARBDYBDY;
            mesh->point[ia].tag |= MG_PARBDYBDY;
            mesh->point[ib].tag |= MG_PARBDYBDY;
          }
          break;
        } else if ( !ph->nxt ) {
          break;
        } else {
          ph = &hash->item[ph->nxt];
        }
      }
    }
  }

  /* set tag */
  for (k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_NOM ) {
        mesh->point[pt->v[MMG5_inxt2[i]]].tag |= MG_NOM;
        mesh->point[pt->v[MMG5_iprv2[i]]].tag |= MG_NOM;
      }
    }
  }

  if ( abs(mesh->info.imprim) > 5 && dup+nmf > 0 ) {
    fprintf(stdout,"  ## ");  fflush(stdout);
    if ( nmf > 0 )  fprintf(stdout,"[non-manifold model]  ");
    if ( dup > 0 )  fprintf(stdout," %" MMG5_PRId " duplicate removed",dup);
    fprintf(stdout,"\n");
  }
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- completed.\n");

  return 1;
}

/**
 * \param mesh pointer to the mesh.
 * \param hash pointer to the hash table to fill.
 * \param ia first vertex of face to hash.
 * \param ib second vertex of face to hash.
 * \param ic third vertex of face to hash.
 * \param k index of face to hash.
 *
 * \return 0 if fail, -1 if the face is newly hashed, index of the first face
 * hashed if another face with same vertices exist.
 *
 *
 **/
MMG5_int MMG5_hashFace(MMG5_pMesh mesh,MMG5_Hash *hash,MMG5_int ia,MMG5_int ib,MMG5_int ic,MMG5_int k) {
  MMG5_hedge     *ph;
  MMG5_int       key;
  MMG5_int       mins,maxs,sum,j;

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
      while ( ph->nxt && ph->nxt < hash->max ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return ph->k;
      }
    }
    ph->nxt = hash->nxt;
    ph      = &hash->item[hash->nxt];
    ph->a   = mins;
    ph->b   = maxs;
    ph->s   = sum;
    ph->k   = k;
    hash->nxt = ph->nxt;
    ph->nxt = 0;

    if ( hash->nxt >= hash->max ) {
      MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,MMG5_GAP,MMG5_hedge,"face",return 0;);
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    return -1;
  }

  /* insert new face */
  ph->a = mins;
  ph->b = maxs;
  ph->s = sum;
  ph->k = k;
  ph->nxt = 0;

  return -1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param hash pointer to the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \param k index of point along the edge.
 *
 * \return 2 if a new edge has been added, 1 if edge was already listed, 0 if fail.
 *
 * Add edge \f$[a;b]\f$ to the hash table.
 *
 */
int MMG5_hashEdge(MMG5_pMesh mesh,MMG5_Hash *hash, MMG5_int a,MMG5_int b,MMG5_int k) {
  MMG5_hedge  *ph;
  MMG5_int     key;
  MMG5_int    ia,ib,j;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a == ia && ph->b == ib )
    return 1;
  else if ( ph->a ) {
    while ( ph->nxt && ph->nxt < hash->max ) {
      ph = &hash->item[ph->nxt];
      if ( ph->a == ia && ph->b == ib )  return 1;
    }
    ph->nxt   = hash->nxt;
    ph        = &hash->item[hash->nxt];

    if ( hash->nxt >= hash->max-1 ) {
      if ( mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: memory alloc problem (edge):"
                " %" MMG5_PRId "\n",__func__,hash->max);

      MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,MMG5_GAP,MMG5_hedge,
                         "MMG5_edge",return 0);
      /* ph pointer may be false after realloc */
      ph        = &hash->item[hash->nxt];

      for (j=ph->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    hash->nxt = ph->nxt;
  }

  /* insert new edge */
  ph->a = ia;
  ph->b = ib;
  ph->k = k;
  ph->nxt = 0;

  return 2;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param hash pointer to the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \param k new index of point along the edge.
 * \return 1 if success, 0 if fail (edge is not found).
 *
 * Update the index of the point stored along the edge \f$[a;b]\f$
 * \note In ParMmg, hash_pmmg.c: PMMG_hashUpdate_all updates k and s at the same time;
 * \note PMMG_hashUpdate_all might be moved here if needed one day in mmg.
 *
 */
int MMG5_hashUpdate(MMG5_Hash *hash, MMG5_int a,MMG5_int b,MMG5_int k) {
  MMG5_hedge  *ph;
  MMG5_int     key;
  MMG5_int    ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->item[key];

  while ( ph->a ) {
    if ( ph->a == ia && ph->b == ib ) {
      ph->k = k;
      return 1;
    }

    if ( !ph->nxt ) return 0;

    ph = &hash->item[ph->nxt];

  }

  return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param hash pointer to the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \param tag edge tag
 * \return the edge tag if success, 0 if fail.
 *
 * Add edge \f$[a;b]\f$ to the hash table if it doesn't exist and store the edge
 * tag. If the edge exist, add the new tag to the already stored tags.
 *
 */
int MMG5_hashEdgeTag(MMG5_pMesh mesh,MMG5_Hash *hash, MMG5_int a,MMG5_int b,uint16_t tag) {
  MMG5_hedge  *ph;
  MMG5_int     key;
  MMG5_int    ia,ib,j;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == ia && ph->b == ib ) {
      ph->k |= tag;
      return ph->k;
    }
    else {
      while ( ph->nxt && ph->nxt < hash->max ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == ia && ph->b == ib )  {
          ph->k |= tag;
          return ph->k;
        }
      }
    }
    ph->nxt   = hash->nxt;
    ph        = &hash->item[hash->nxt];
    ph->a     = ia;
    ph->b     = ib;
    ph->k     = tag;
    hash->nxt = ph->nxt;
    ph->nxt   = 0;

    if ( hash->nxt >= hash->max ) {
      MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,mesh->gap,MMG5_hedge,
                        "edge hash table",return 0);
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    return tag;
  }

  /* insert new edge */
  ph->a     = ia;
  ph->b     = ib;
  ph->k     = tag;
  ph->nxt   = 0;

  return tag;
}

/**
 * \param hash pointer to the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \return the index of point stored along \f$[a;b]\f$.
 *
 * Find the index of point stored along \f$[a;b]\f$.
 * \note In ParMmg, hash_pmmg.c: PMMG_hashGet_all gets k and s at the same time;
 * \note PMMG_hashGet_all might be moved here if needed one day in mmg.
 *
 */
MMG5_int MMG5_hashGet(MMG5_Hash *hash,MMG5_int a,MMG5_int b) {
  MMG5_hedge  *ph;
  MMG5_int    key;
  MMG5_int    ia,ib;

  if ( !hash->item ) return 0;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->a )  return 0;
  if ( ph->a == ia && ph->b == ib )  return ph->k;
  while ( ph->nxt ) {
    ph = &hash->item[ph->nxt];
    if ( ph->a == ia && ph->b == ib )  return ph->k;
  }
  return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param hash pointer to the hash table of edges.
 * \param hsiz initial size of hash table.
 * \param hmax maximal size of hash table.
 * \return 1 if success, 0 if fail.
 *
 * Hash edges or faces.
 *
 */
int MMG5_hashNew(MMG5_pMesh mesh,MMG5_Hash *hash,MMG5_int hsiz,MMG5_int hmax) {
  MMG5_int   k;

  /* adjust hash table params */
  hash->siz  = hsiz+1;
  hash->max  = hmax + 2;
  hash->nxt  = hash->siz;

  MMG5_ADD_MEM(mesh,(hash->max+1)*sizeof(MMG5_hedge),"hash table",
                return 0);
  MMG5_SAFE_CALLOC(hash->item,(hash->max+1),MMG5_hedge,return 0);

  for (k=hash->siz; k<hash->max; k++)
    hash->item[k].nxt = k+1;

  return 1;
}
