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
 * \file mmgs/hash_s.c
 * \brief Functions for hash tables management and triangle packing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmgs_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 if fail.
 *
 * tria packing
 *
 */
static int paktri(MMG5_pMesh mesh) {
  MMG5_pTria   pt,pt1;
  MMG5_int     k;

  k = 1;
  do {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) {
      pt1 = &mesh->tria[mesh->nt];
      assert ( pt && pt1 );
      memcpy(pt,pt1,sizeof(MMG5_Tria));
      if ( !MMGS_delElt(mesh,mesh->nt) )  return 0;
    }
  }
  while ( ++k < mesh->nt );

  /* Recreate nil chain */
  mesh->nenil = mesh->nt + 1;

  for(k=mesh->nenil; k<=mesh->ntmax-1; k++){
    mesh->tria[k].v[2] = k+1;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 if fail.
 *
 * Create adjacency table.
 *
 */
int MMGS_hashTria(MMG5_pMesh mesh) {
  MMG5_Hash          hash;

  if ( mesh->adja )  return 1;
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING STRUCTURE\n");

  /* tassage */
  if ( !paktri(mesh) )  return 0;

  MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(MMG5_int),"adjacency table",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,MMG5_int,return 0);

  if ( !MMG5_mmgHashTria(mesh, mesh->adja, &hash, 0) ) return 0;

  MMG5_DEL_MEM(mesh,hash.item);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0.
 *
 * Copy the properties (ref and tag) of the declared edges to the triangles,
 * where they are assigned to the individual corners of the triangle. First a
 * hash is created for rapid lookup of the edges. Then in a loop over all edges
 * of all triangles, the hash is probed for each edge, and if it exists its
 * properties are copied. Thus, declared edges that do not occur in any triangle
 * will be silently ignored.
 *
 * \remark this function handle all the provided edges.
 *
 */
int MMGS_assignEdge(MMG5_pMesh mesh) {
  MMG5_Hash  hash;
  MMG5_pTria  pt;
  MMG5_pEdge  pa;
  MMG5_int    k,ia;
  int8_t      i,i1,i2;

  if ( !mesh->na ) return 1;

  /* adjust hash table params */
  hash.siz  = mesh->na;
  hash.max  = 3*mesh->na+1;
  MMG5_ADD_MEM(mesh,(hash.max+1)*sizeof(MMG5_hedge),"hash table",return 0);
  MMG5_SAFE_CALLOC(hash.item,hash.max+1,MMG5_hedge,return 0);

  hash.nxt  = mesh->na;
  for (k=mesh->na; k<hash.max; k++)
    hash.item[k].nxt = k+1;

  /* hash mesh edges */
  for (k=1; k<=mesh->na; k++)
    MMG5_hashEdge(mesh,&hash,mesh->edge[k].a,mesh->edge[k].b,k);

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

  /* reset edge structure */
  MMG5_DEL_MEM(mesh,hash.item);
  MMG5_DEL_MEM(mesh,mesh->edge);
  mesh->na = 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0.
 *
 * Copy the edge tags stored in triangles in the other triangles sharing the
 * edge.
 *
 */
int MMGS_bdryUpdate(MMG5_pMesh mesh) {
  MMG5_Hash   hash;
  MMG5_pTria  pt;
  MMG5_int    k,nad;
  int         tag;
  int8_t      i,i1,i2;

  /* adjust hash table params */
  /* Euler formula : na ~ 3np */
  if ( !MMG5_hashNew(mesh,&hash,3*mesh->np,9*mesh->np) ) {
    printf("  # Error: %s: Not enough memory to allocate edge hash table",__func__);
  }

  /* hash tagged edges */
  nad = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for ( i=0; i<3; ++i ) {
      if ( mesh->tria[k].tag[i] ) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        if ( !MMG5_hashEdgeTag(mesh,&hash,pt->v[i1],pt->v[i2],pt->tag[i]) ) {
          printf("  # Error: %s: Lack of memory.",__func__);
          return 0;
        }
        ++nad;
      }
    }
  }

  /* update tags */
  if ( nad ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];

        tag = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( tag ) {
          pt->tag[i] |= tag;
        }
      }
    }
  }

  /* free hash structure */
  MMG5_DEL_MEM(mesh,hash.item);

  return 1;
}
