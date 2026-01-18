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
 * \file mmgs/librnbg_s.c
 * \brief Functions for scotch renumerotation.
 * \author Cedric Lachat (Inria/UBordeaux)
 * \version 5
 * \date 2013
 * \copyright GNU Lesser General Public License.
 */

#include "libmmgs_private.h"
#include "libmmgs.h"

#ifdef USE_SCOTCH

#include "librnbg_private.h"

/**
 * \param trias pointer to a table containing the tetra structures.
 * \param *perm pointer to the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first tria to swap.
 * \param ind2 index of the second tria to swap.
 *
 * Swap two tetras in the table of tetrahedras.
 *
 */
static inline
void MMG5_swapTri(MMG5_pTria trias, MMG5_int* perm, MMG5_int ind1, MMG5_int ind2) {
  MMG5_Tria pttmp;
  MMG5_int  tmp;

  /* 2-- swap the triangles */
  memcpy(&pttmp       ,&trias[ind2],sizeof(MMG5_Tria));
  memcpy(&trias[ind2],&trias[ind1],sizeof(MMG5_Tria));
  memcpy(&trias[ind1],&pttmp       ,sizeof(MMG5_Tria));

  /* 3-- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param boxVertNbr number of vertices by box.
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure
 * \param fields pointer to an array of solution fields
 * \param permNodGlob array to store the global permutation of nodes (if provided)
 *
 * \return 0 if the renumbering fail and we can't rebuild tetrahedra hashtable,
 * 1 if the renumbering fail but we can rebuild tetrahedra hashtable or
 * if the renumbering success.
 *
 * Modifies the node indicies to prevent from cache missing.
 *
 */
int MMG5_mmgsRenumbering(int boxVertNbr, MMG5_pMesh mesh, MMG5_pSol sol,
                         MMG5_pSol fields,MMG5_int* permNodGlob) {
  MMG5_pPoint  ppt;
  MMG5_pTria   ptri;
  SCOTCH_Num   edgeNbr;
  SCOTCH_Num   *vertTab, *edgeTab, *permVrtTab;
  SCOTCH_Graph graf ;
  MMG5_int     vertNbr, nodeGlbIdx, triaIdx, ballTriIdx;
  MMG5_int     j, k, edgeSiz;
  MMG5_int     *vertOldTab, *permNodTab, ntreal, npreal;
  MMG5_int     *adja,iadr;
  int          i;


  /* Computing the number of vertices and a contiguous tabular of vertices */
  vertNbr = 0;

  MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_int),"vertOldTab",return 1);
  MMG5_SAFE_CALLOC(vertOldTab,mesh->nt+1,MMG5_int,return 1);

  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {

    /* Testing if the tria exists */
    if (!mesh->tria[triaIdx].v[0]) continue;
    vertOldTab[triaIdx] = ++vertNbr;
  }

  if ( vertNbr/2 < MMG5_BOXSIZE ) {
    /* not enough tetra to renum */
    MMG5_DEL_MEM(mesh,vertOldTab);
    return 1;
  }
  /* Allocating memory to compute adjacency lists */
  MMG5_ADD_MEM(mesh,(vertNbr+2)*sizeof(SCOTCH_Num),"vertTab",
                MMG5_DEL_MEM(mesh,vertOldTab);
                return 1);
  MMG5_SAFE_CALLOC(vertTab,vertNbr+2,SCOTCH_Num,return 1);

  if (!memset(vertTab, ~0, sizeof(SCOTCH_Num)*(vertNbr + 2))) {
    perror("  ## Memory problem: memset");
    MMG5_DEL_MEM(mesh,vertOldTab);
    MMG5_DEL_MEM(mesh,vertTab);
    return 1;
  }

  edgeNbr = 1;
  /* Euler-Poincare formulae edgeSiz = 20*mesh->np~4*mesh->nt;
     (2*(12*mesh->np (triangles)-2*mesh->np (boundary triangles))) */
  edgeSiz = vertNbr*4;

  MMG5_ADD_MEM(mesh,edgeSiz*sizeof(SCOTCH_Num),"edgeTab",
                MMG5_DEL_MEM(mesh,vertOldTab);
                MMG5_DEL_MEM(mesh,vertTab);
                return 1);
  MMG5_SAFE_CALLOC(edgeTab,edgeSiz,SCOTCH_Num,return 1);


  /* Computing the adjacency list for each vertex */
  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {

    /* Testing if the tria exists */
    if (!mesh->tria[triaIdx].v[0]) continue;

    iadr = 3*(triaIdx-1) + 1;
    adja = &mesh->adja[iadr];
    for (i=0; i<3; i++) {
      ballTriIdx = adja[i] / 3;

      if (!ballTriIdx) continue;

      /* Testing if one neighbour of triaIdx has already been added */
      if (vertTab[vertOldTab[triaIdx]] < 0)
        vertTab[vertOldTab[triaIdx]] = edgeNbr;

      /* Testing if edgeTab memory is enough */
      if (edgeNbr >= edgeSiz) {
        int oldsize = edgeSiz;
        MMG5_ADD_MEM(mesh,0.2*sizeof(SCOTCH_Num),"edgeTab",
                      MMG5_DEL_MEM(mesh,vertOldTab);
                      MMG5_DEL_MEM(mesh,vertTab);
                      return 1);
        edgeSiz *= 1.2;
        MMG5_SAFE_REALLOC(edgeTab,oldsize,edgeSiz,SCOTCH_Num,"scotch table",return 1);
      }

      edgeTab[edgeNbr++] = vertOldTab[ballTriIdx];
    }
  }
  vertTab[vertNbr+1] = edgeNbr;
  edgeNbr--;
  /*check if some tria are alone*/
  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {

    /* Testing if the tria exists */
    if (!mesh->tria[triaIdx].v[0]) continue;
    if (vertTab[vertOldTab[triaIdx]] < 0) {
      if(vertOldTab[triaIdx] == vertNbr) {
        fprintf(stderr,"  ## Warning: %s: graph error, no renumbering.\n",
                __func__);
        MMG5_DEL_MEM(mesh,edgeTab);
        MMG5_DEL_MEM(mesh,vertTab);
        return 1;
      }
      if(vertTab[vertOldTab[triaIdx] + 1] > 0)
        vertTab[vertOldTab[triaIdx]] = vertTab[vertOldTab[triaIdx] + 1];
      else {
        if(vertOldTab[triaIdx]+1 == vertNbr) {
          fprintf(stderr,"  ## Warning: %s: graph error, no renumbering.\n",
                  __func__);
          MMG5_DEL_MEM(mesh,edgeTab);
          MMG5_DEL_MEM(mesh,vertTab);
          return 1;
        }
        i = 1;
        do  {
          i++;
        } while((vertTab[vertOldTab[triaIdx] + i] < 0) && ((vertOldTab[triaIdx] + i) < vertNbr));
        if(vertOldTab[triaIdx] + i == vertNbr) {
          fprintf(stderr,"  ## Warning: %s: graph error, no renumbering.\n",
                  __func__);
          MMG5_DEL_MEM(mesh,edgeTab);
          MMG5_DEL_MEM(mesh,vertTab);
          return 1;
        }
        vertTab[vertOldTab[triaIdx]] = vertTab[vertOldTab[triaIdx] + i];
      }
    }
  }

  /* free adjacents to gain memory space */
  MMG5_DEL_MEM(mesh,mesh->adja);

  /* Building the graph by calling Scotch functions */
  SCOTCH_graphInit(&graf) ;
  CHECK_SCOTCH(SCOTCH_graphBuild(&graf, (SCOTCH_Num) 1, vertNbr, vertTab+1,
                                 NULL, NULL, NULL, edgeNbr, edgeTab+1, NULL),
               "scotch_graphbuild", 0) ;

#ifndef NDEBUG
  /* don't check in release mode */
  if ( mesh->info.imprim > 6 || mesh->info.ddebug )
    fprintf(stdout,"checking graph...\n");

  CHECK_SCOTCH(SCOTCH_graphCheck(&graf), "scotch_graphcheck", 0);
#endif

  MMG5_ADD_MEM(mesh,(vertNbr+1)*sizeof(SCOTCH_Num),"permVrtTab",
               MMG5_DEL_MEM(mesh,vertOldTab);
               MMG5_DEL_MEM(mesh,vertTab);
               MMG5_DEL_MEM(mesh,edgeTab);
               if( !MMGS_hashTria(mesh) ) return 0;
               return 1);
  MMG5_SAFE_CALLOC(permVrtTab,vertNbr+1,SCOTCH_Num,return 1);

  CHECK_SCOTCH(MMG5_kPartBoxCompute(&graf, vertNbr, boxVertNbr, permVrtTab, mesh),
               "boxCompute", 0);

  SCOTCH_graphExit(&graf) ;

  MMG5_DEL_MEM(mesh,edgeTab);
  MMG5_DEL_MEM(mesh,vertTab);

  /* Computing the new point list and modifying the adja strcuture */
  ntreal = 0;
  npreal = 0;

  /* Create the final permutation table for trias (stored in vertOldTab) */
  for (j=1; j<= mesh->nt; j++) {
    if ( !mesh->tria[triaIdx].v[0] )  continue;

    vertOldTab[triaIdx] = permVrtTab[vertOldTab[triaIdx]];
  }
  MMG5_DEL_MEM(mesh,permVrtTab);

  /* Permute triangles */
  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {
    while ( vertOldTab[triaIdx] != triaIdx  && vertOldTab[triaIdx] )
      MMG5_swapTri(mesh->tria,vertOldTab,triaIdx,vertOldTab[triaIdx]);
  }
  MMG5_DEL_MEM(mesh,vertOldTab);

  MMG5_ADD_MEM(mesh,(mesh->np+1)*sizeof(MMG5_int),"permNodTab",
                if( !MMGS_hashTria(mesh) ) return 0;
                return 1);
  MMG5_SAFE_CALLOC(permNodTab,mesh->np+1,MMG5_int,return 1);

  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {
    ptri = &mesh->tria[triaIdx];

    /* Testing if the tria exists */
    if (!ptri->v[0]) continue;

    ntreal++;

    for(j = 0 ; j <= 2 ; j++) {

      nodeGlbIdx = ptri->v[j];

      if ( permNodTab[nodeGlbIdx] ) continue;

      ppt = &mesh->point[nodeGlbIdx];

      if ( !(ppt->tag & MG_NUL) )
        /* Building the new point list */
        permNodTab[nodeGlbIdx] = ++npreal;
    }
  }

  /* Append unseen required points for orphan points preservation */
  for ( k=1; k<=mesh->np; ++k) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) {
      continue;
    }
    if ( permNodTab[k] ) continue;


    if ( ppt->tag & MG_REQ ) {
      /* Add orphan required point to permnodtab */
      permNodTab[k] = ++npreal;
    }
  }

  /* Create the final permutation table for trias (stored in vertOldTab) and *
     modify the numbering of the nodes of each tria */
  for( triaIdx = 1; triaIdx < ntreal + 1; triaIdx++) {
    for(j = 0 ; j < 3 ; j++) {
      mesh->tria[triaIdx].v[j] = permNodTab[mesh->tria[triaIdx].v[j]];
    }
  }

  /* If needed, store update the global permutation for point array */
  if ( permNodGlob ) {
    for ( k=1; k<=mesh->np; ++k ) {
      permNodGlob[k] = permNodTab[permNodGlob[k]];
    }
  }

  /* Permute nodes and sol */
  for (j=1; j<= mesh->np; j++) {
    while ( permNodTab[j] != j && permNodTab[j] )
      MMG5_swapNod(mesh,mesh->point,sol->m,fields,permNodTab,j,permNodTab[j],sol->size);
  }
  MMG5_DEL_MEM(mesh,permNodTab);

  mesh->nt = ntreal;
  mesh->np = npreal;

  if ( mesh->np == mesh->npmax )
    mesh->npnil = 0;
  else
    mesh->npnil = mesh->np + 1;

  if ( mesh->nt == mesh->ntmax )
    mesh->nenil = 0;
  else
    mesh->nenil = mesh->nt + 1;

  if ( mesh->npnil ) {
    for (k=mesh->npnil; k<mesh->npmax-1; k++) {
      mesh->point[k].tmp  = k+1;
    }
    mesh->point[mesh->npmax-1].tmp = 0;
    mesh->point[mesh->npmax  ].tmp = 0;
  }

  if ( mesh->nenil ) {
    for (k=mesh->nenil; k<mesh->ntmax-1; k++) {
      mesh->tria[k].v[2] = k+1;
    }
    mesh->tria[mesh->ntmax-1].v[2] = 0;
    mesh->tria[mesh->ntmax  ].v[2] = 0;
  }

  if( !MMGS_hashTria(mesh) ) return 0;

  return 1;
}
#endif
