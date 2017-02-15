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
 * \file mmgs/librnbg_s.c
 * \brief Functions for scotch renumerotation.
 * \author Cedric Lachat (Inria/UBordeaux)
 * \version 5
 * \date 2013
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"

#ifdef USE_SCOTCH

#include "librnbg.h"

/**
 * \param trias pointer toward a table containing the tetra structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first tria to swap.
 * \param ind2 index of the second tria to swap.
 *
 * Swap two tetras in the table of tetrahedras.
 *
 */
static inline
void _MMG5_swapTri(MMG5_pTria trias, int* perm, int ind1, int ind2) {
  MMG5_Tria pttmp;
  int   tmp;

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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward he solution structure
 * \return 0 if the renumbering fail and we can't rebuild tetrahedra hashtable,
 * 1 if the renumbering fail but we can rebuild tetrahedra hashtable or
 * if the renumbering success.
 *
 * Modifies the node indicies to prevent from cache missing.
 *
 */
int _MMG5_mmgsRenumbering(int boxVertNbr, MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pPoint ppt;
  MMG5_pTria ptri;
  SCOTCH_Num edgeNbr;
  SCOTCH_Num *vertTab, *edgeTab, *permVrtTab;
  SCOTCH_Graph graf ;
  int    vertNbr, nodeGlbIdx, triaIdx, ballTriIdx;
  int    i, j, k;
  int    edgeSiz;
  int    *vertOldTab, *permNodTab, ntreal, npreal;
  int    *adja,iadr;


  /* Computing the number of vertices and a contiguous tabular of vertices */
  vertNbr = 0;

  _MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(int),"vertOldTab",return(1));
  _MMG5_SAFE_CALLOC(vertOldTab,mesh->nt+1,int);

  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {

    /* Testing if the tria exists */
    if (!mesh->tria[triaIdx].v[0]) continue;
    vertOldTab[triaIdx] = ++vertNbr;
  }

  if ( vertNbr/2 < _MMG5_BOXSIZE ) {
    /* not enough tetra to renum */
    _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->nt+1)*sizeof(int));
    return(1);
  }
  /* Allocating memory to compute adjacency lists */
  _MMG5_ADD_MEM(mesh,(vertNbr+2)*sizeof(SCOTCH_Num),"vertTab",
                _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                return(1));
  _MMG5_SAFE_CALLOC(vertTab,vertNbr+2,SCOTCH_Num);

  if (!memset(vertTab, ~0, sizeof(SCOTCH_Num)*(vertNbr + 2))) {
    perror("  ## Memory problem: memset");
    _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->nt+1)*sizeof(int));
    _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
    return 1;
  }

  edgeNbr = 1;
  /* Euler-Poincare formulae edgeSiz = 20*mesh->np~4*mesh->nt;
     (2*(12*mesh->np (triangles)-2*mesh->np (boundary triangles))) */
  edgeSiz = vertNbr*4;

  _MMG5_ADD_MEM(mesh,edgeSiz*sizeof(SCOTCH_Num),"edgeTab",
                _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->nt+1)*sizeof(int));
                _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                return(1));
  _MMG5_SAFE_CALLOC(edgeTab,edgeSiz,SCOTCH_Num);


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
        _MMG5_ADD_MEM(mesh,0.2*sizeof(SCOTCH_Num),"edgeTab",
                      _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                      _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                      return(1));
        edgeSiz *= 1.2;
        _MMG5_SAFE_REALLOC(edgeTab,edgeSiz,SCOTCH_Num,"scotch table");
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
        fprintf(stdout,"WARNING graph problem, no renum\n");
        _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
        _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
        return(1);
      }
      if(vertTab[vertOldTab[triaIdx] + 1] > 0)
        vertTab[vertOldTab[triaIdx]] = vertTab[vertOldTab[triaIdx] + 1];
      else {
        if(vertOldTab[triaIdx]+1 == vertNbr) {
          fprintf(stdout,"WARNING graph problem, no renum\n");
          _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
          _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
          return(1);
        }
        i = 1;
        do  {
          i++;
        } while((vertTab[vertOldTab[triaIdx] + i] < 0) && ((vertOldTab[triaIdx] + i) < vertNbr));
        if(vertOldTab[triaIdx] + i == vertNbr) {
          fprintf(stdout,"WARNING graph problem, no renum\n");
          _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
          _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
          return(1);
        }
        vertTab[vertOldTab[triaIdx]] = vertTab[vertOldTab[triaIdx] + i];
      }
    }
  }

  /* free adjacents to gain memory space */
  _MMG5_DEL_MEM(mesh,mesh->adja,(3*mesh->ntmax+5)*sizeof(int));

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

  _MMG5_ADD_MEM(mesh,(vertNbr+1)*sizeof(SCOTCH_Num),"permVrtTab",
                _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
                if( !_MMGS_hashTria(mesh) ) return(0);
                return(1));
  _MMG5_SAFE_CALLOC(permVrtTab,vertNbr+1,SCOTCH_Num);

  CHECK_SCOTCH(_MMG5_kPartBoxCompute(graf, vertNbr, boxVertNbr, permVrtTab, mesh),
               "boxCompute", 0);

  SCOTCH_graphExit(&graf) ;

  _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
  _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));

  /* Computing the new point list and modifying the adja strcuture */
  _MMG5_ADD_MEM(mesh,(mesh->np+1)*sizeof(int),"permNodTab",
                _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                _MMG5_DEL_MEM(mesh,permVrtTab,(vertNbr+1)*sizeof(SCOTCH_Num));
                _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
                if( !_MMGS_hashTria(mesh) ) return(0);
                return(1));
  _MMG5_SAFE_CALLOC(permNodTab,mesh->np+1,int);

  ntreal = 0;
  npreal = 0;

  /* Create the final permutation table for trias (stored in vertOldTab) */
  for (j=1; j<= mesh->nt; j++) {
    if ( !mesh->tria[triaIdx].v[0] )  continue;

    vertOldTab[triaIdx] = permVrtTab[vertOldTab[triaIdx]];
  }
  _MMG5_DEL_MEM(mesh,permVrtTab,(vertNbr+1)*sizeof(SCOTCH_Num));

  /* Permute triangles */
  for(triaIdx = 1 ; triaIdx < mesh->nt + 1 ; triaIdx++) {
    while ( vertOldTab[triaIdx] != triaIdx  && vertOldTab[triaIdx] )
      _MMG5_swapTri(mesh->tria,vertOldTab,triaIdx,vertOldTab[triaIdx]);
  }

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
  _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));

  /* Create the final permutation table for trias (stored in vertOldTab) and *
     modify the numbering of the nodes of each tria */
  for( triaIdx = 1; triaIdx < ntreal + 1; triaIdx++) {
    for(j = 0 ; j < 3 ; j++) {
      mesh->tria[triaIdx].v[j] = permNodTab[mesh->tria[triaIdx].v[j]];
    }
  }

  /* Permute nodes and sol */
  for (j=1; j<= mesh->np; j++) {
    while ( permNodTab[j] != j && permNodTab[j] )
      _MMG5_swapNod(mesh->point,sol->m,permNodTab,j,permNodTab[j],sol->size);
  }
  _MMG5_DEL_MEM(mesh,permNodTab,(mesh->np+1)*sizeof(int));

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

  if ( mesh->npnil )
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;

  if ( mesh->nenil )
    for (k=mesh->nenil; k<mesh->ntmax-1; k++)
      mesh->tria[k].v[2] = k+1;

  if( !_MMGS_hashTria(mesh) ) return(0);

  return 1;
}
#endif
