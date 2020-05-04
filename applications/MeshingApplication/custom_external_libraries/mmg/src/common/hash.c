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

#include "mmgcommon.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param adjt pointer toward the adjacency table of the surfacic mesh.
 * \param hash pointer toward the edge hash table.
 * \param chkISO flag to say if we check ISO references
 * (so if we come from mmg3d).
 * \return 1 if success, 0 otherwise.
 *
 * Create surface adjacency
 *
 * \remark the ph->s field computation is useless in mmgs.
 *
 */
int MMG5_mmgHashTria(MMG5_pMesh mesh, int *adjt, MMG5_Hash *hash, int chkISO) {
  MMG5_pTria     pt,pt1;
  MMG5_hedge    *ph;
  int            *adja,k,kk,jel,lel,hmax,dup,nmf,ia,ib;
  char           i,i1,i2,j,l;
  unsigned int   key;

  /* adjust hash table params */
  hmax = 3.71*mesh->np;
  hash->siz  = mesh->np;
  hash->max  = hmax + 1;
  hash->nxt  = hash->siz;
  MMG5_ADD_MEM(mesh,(hash->max+1)*sizeof(MMG5_hedge),"hash table",return 0);
  MMG5_SAFE_CALLOC(hash->item,hash->max+1,MMG5_hedge,return 0);

  for (k=hash->siz; k<hash->max; k++)
    hash->item[k].nxt = k+1;

  if ( mesh->info.ddebug )  fprintf(stdout,"  h- stage 1: init\n");

  /* hash triangles */
  mesh->base = 1;
  dup = nmf = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    pt->flag = 0;
    pt->base = mesh->base;
    adja = &adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      /* compute key */
      ia  = MG_MIN(pt->v[i1],pt->v[i2]);
      ib  = MG_MAX(pt->v[i1],pt->v[i2]);
      key = (MMG5_KA*ia + MMG5_KB*ib) % hash->siz;
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
      char ok = 0;
#endif
      while ( ph->a ) {
        if ( ph->a == ia && ph->b == ib ) {
          jel = ph->k / 3;
          j   = ph->k % 3;
          pt1 = &mesh->tria[jel];
          /* discard duplicate face */
          if ( pt1->v[j] == pt->v[i] ) {
            pt1->v[0] = 0;
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
            if ( chkISO && ( (pt->ref == MG_ISO) || (pt->ref < 0)) ) {
              lel = adjt[3*(jel-1)+1+j]/3;
              l   = adjt[3*(jel-1)+1+j]%3;
              adjt[3*(lel-1)+1+l] = 0;
              adja[i] = 3*jel+j;
              adjt[3*(jel-1)+1+j] = 3*k + i;
              (mesh->tria[jel]).tag[j] |= MG_GEO + MG_NOM;
              (mesh->tria[lel]).tag[l] |= MG_GEO + MG_NOM;
            }
            else {
              lel = adjt[3*(jel-1)+1+j]/3;
              l   = adjt[3*(jel-1)+1+j]%3;
              (mesh->tria[lel]).tag[l] |= MG_GEO + MG_NOM;
              pt1->tag[j] |= MG_GEO + MG_NOM;
            }
            pt->tag[i] |= MG_GEO + MG_NOM;
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
                      " %d\n",__func__,hash->max);
            }
            MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,0.2,MMG5_hedge,
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
    if ( dup > 0 )  fprintf(stdout," %d duplicate removed",dup);
    fprintf(stdout,"\n");
  }
  if ( mesh->info.ddebug )  fprintf(stdout,"  h- completed.\n");
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \param k index of point along the edge.
 * \return 1 if success, 0 if fail.
 *
 * Add edge \f$[a;b]\f$ to the hash table.
 *
 */
int MMG5_hashEdge(MMG5_pMesh mesh,MMG5_Hash *hash, int a,int b,int k) {
  MMG5_hedge  *ph;
  int          key,ia,ib,j;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*ia + MMG5_KB*ib) % hash->siz;
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
                " %d\n",__func__,hash->max);

      MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,0.2,MMG5_hedge,
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

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \param k new index of point along the edge.
 * \return 1 if success, 0 if fail (edge is not found).
 *
 * Update the index of the point stored along the edge \f$[a;b]\f$
 *
 */
int MMG5_hashUpdate(MMG5_Hash *hash, int a,int b,int k) {
  MMG5_hedge  *ph;
  int          key,ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*ia + MMG5_KB*ib) % hash->siz;
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
 * \param hash pointer toward the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \return the index of point stored along \f$[a;b]\f$.
 *
 * Find the index of point stored along  \f$[a;b]\f$.
 *
 */
int MMG5_hashGet(MMG5_Hash *hash,int a,int b) {
  MMG5_hedge  *ph;
  int          key,ia,ib;

  if ( !hash->item ) return 0;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*ia + MMG5_KB*ib) % hash->siz;
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
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \param hsiz initial size of hash table.
 * \param hmax maximal size of hash table.
 * \return 1 if success, 0 if fail.
 *
 * Hash edges or faces.
 *
 */
int MMG5_hashNew(MMG5_pMesh mesh,MMG5_Hash *hash,int hsiz,int hmax) {
  int   k;

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
