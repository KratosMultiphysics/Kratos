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
 * \file mmg3d/librnbg_3d.c
 * \brief Functions for scotch renumerotation.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \author Cedric Lachat (Inria/UBordeaux)
 * \version 5
 * \date 2013
 * \copyright GNU Lesser General Public License.
 */

#include "mmg3d.h"

#ifdef USE_SCOTCH

#include "librnbg.h"


/**
 * \param graf the input graph structure.
 * \param vertNbr the number of vertices.
 * \param boxVertNbr the number of vertices of each box.
 * \param permVrtTab the new numbering.
 * \param mesh pointer toward the mesh structure.
 * \return 0 if ok, 1 otherwise.
 *
 * Internal function that computes a new numbering of graph vertices using a
 * bipartitioning.
 *
 * \warning Not used.
 **/
static inline
int _MMG5_biPartBoxCompute(SCOTCH_Graph graf, int vertNbr, int boxVertNbr, SCOTCH_Num *permVrtTab,MMG5_pMesh mesh) {
  int boxNbr, vertIdx, boxIdx;
  SCOTCH_Num tmp, tmp2, *partTab, *partNumTab, *partPrmTab;
  SCOTCH_Strat strat ;

  /* Computing the number of boxes */
  boxNbr = vertNbr / boxVertNbr;
  if (boxNbr * boxVertNbr != vertNbr) {
    boxNbr = boxNbr + 1;
  }


  /* Initializing SCOTCH functions */
  CHECK_SCOTCH(SCOTCH_stratInit(&strat), "scotch_stratInit", 0) ;
#ifdef SCOTCH_6
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, "r{job=t,map=t,poli=S,sep=m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}|m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}}"), "scotch_stratGraphMap", 0) ;
#else
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, "r{job=t,map=t,poli=S,sep=m{type=h,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}|m{,vert=80,low=h{pass=10}f{bal=0.005,move=0},asc=b{bnd=f{bal=0.05,move=0},org=f{bal=0.05,move=0}}}}"), "scotch_stratGraphMap", 0) ;
#endif

  _MMG5_ADD_MEM(mesh,vertNbr*sizeof(SCOTCH_Num),"partTab",return(1));
  _MMG5_SAFE_CALLOC(partTab,vertNbr,SCOTCH_Num);

  /* Partionning the graph */
  CHECK_SCOTCH(SCOTCH_graphPart(&graf, boxNbr, &strat, partTab), "scotch_graphPart", 0);

  _MMG5_ADD_MEM(mesh,boxNbr*sizeof(SCOTCH_Num),"partNumTab",return(1));
  _MMG5_SAFE_CALLOC(partNumTab,boxNbr,SCOTCH_Num);

  if (!memset(partNumTab, 0, boxNbr*sizeof(SCOTCH_Num))) {
    perror("  ## Memory problem: memset");
    return 1;
  }

  /* Computing the number of elements of each box */
  for( vertIdx = 0 ; vertIdx< vertNbr ;vertIdx++)
    partNumTab[partTab[vertIdx]] += 1;


  /* partition permutation tabular */
  _MMG5_ADD_MEM(mesh,(vertNbr+1)*sizeof(SCOTCH_Num),"partPrmTab",return(1));
  _MMG5_SAFE_CALLOC(partPrmTab,vertNbr+1,SCOTCH_Num);

  /* Copying the previous tabular in order to have the index of the first
   * element of each box
   * */
  tmp = partNumTab[0];
  partNumTab[0] = 0;
  for(boxIdx = 1; boxIdx < boxNbr ; boxIdx++) {
    tmp2 = partNumTab[boxIdx];
    partNumTab[boxIdx] = partNumTab[boxIdx-1] + tmp;
    tmp = tmp2;
  }

  /* partPrmTab is built such as each vertex belongs to his box */
  for( vertIdx = 0;vertIdx< vertNbr;vertIdx++)
    partPrmTab[partNumTab[partTab[vertIdx]]++] = vertIdx;


  /* Infering the new numbering */
  for (vertIdx = 0; vertIdx < vertNbr ; vertIdx++)
    permVrtTab[partPrmTab[vertIdx] + 1] = vertIdx + 1;

  _MMG5_DEL_MEM(mesh,partTab,vertNbr*sizeof(SCOTCH_Num));
  _MMG5_DEL_MEM(mesh,partNumTab,boxNbr*sizeof(SCOTCH_Num));
  _MMG5_DEL_MEM(mesh,partPrmTab,(vertNbr+1)*sizeof(SCOTCH_Num));

  SCOTCH_stratExit(&strat) ;
  return 0;
}


/**
 * \param tetras pointer toward a table containing the tetra structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first tetra to swap.
 * \param ind2 index of the second tetra to swap.
 *
 * Swap two tetras in the table of tetrahedras.
 *
 */
static inline
void _MMG5_swapTet(MMG5_pTetra tetras/*, int* adja*/, int* perm, int ind1, int ind2) {
  MMG5_Tetra pttmp;
  int        tmp;

  /* Commentated part: swap for adja table if we don't free it in renumbering *
   * function (faster but need of 4*mesh->nemax*sizeof(int) extra bytes ) */
  /* int   adjatmp,kadj,ifadj,j; */

  /* /\* 1-- swap the adja table *\/ */
  /* /\* First: replace ind2 by ind1 in adjacents tetras of ind2 *\/ */
  /* for ( j=1; j<5; j++ ) { */
  /*   if ( adja[4*(ind2-1)+j]/4 ) { */
  /*     kadj  = adja[4*(ind2-1)+j]/4; */
  /*     ifadj = adja[4*(ind2-1)+j]%4; */
  /*     adja[4*(kadj-1)+1+ifadj] = 4*ind1 + adja[4*(kadj-1)+1+ifadj]%4; */
  /*   } */
  /* } */

  /* /\* Second: replace ind1 by ind2 in adjacents tetras of ind1*\/ */
  /* for ( j=1; j<5; j++ ) { */
  /*   if ( adja[4*(ind1-1)+j]/4 ) { */
  /*     kadj  = adja[4*(ind1-1)+j]/4; */
  /*     ifadj = adja[4*(ind1-1)+j]%4; */
  /*     if ( kadj == ind1 ) */
  /*       adja[4*(ind2-1)+1+ifadj] = 4*ind2 + adja[4*(ind2-1)+1+ifadj]%4; */
  /*     else */
  /*       adja[4*(kadj-1)+1+ifadj] = 4*ind2 + adja[4*(kadj-1)+1+ifadj]%4; */
  /*   } */
  /* } */

  /* /\* Third: swap adjacents for ind1 and ind2 *\/ */
  /* for ( j=1; j<5; j++ ) { */
  /*   adjatmp = adja[4*(ind2-1)+j]; */
  /*   adja[4*(ind2-1)+j] = adja[4*(ind1-1)+j]; */
  /*   adja[4*(ind1-1)+j] = adjatmp; */
  /* } */

  /* 2-- swap the tetrahedras */
  memcpy(&pttmp       ,&tetras[ind2],sizeof(MMG5_Tetra));
  memcpy(&tetras[ind2],&tetras[ind1],sizeof(MMG5_Tetra));
  memcpy(&tetras[ind1],&pttmp       ,sizeof(MMG5_Tetra));

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
int _MMG5_mmg3dRenumbering(int boxVertNbr, MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pPoint ppt;
  MMG5_pTetra ptet;
  MMG5_pPrism pp;
  SCOTCH_Num  edgeNbr;
  SCOTCH_Num  *vertTab, *edgeTab, *permVrtTab;
  SCOTCH_Graph graf ;
  int    vertNbr, nodeGlbIdx, tetraIdx, ballTetIdx;
  int    i, j, k;
  int    edgeSiz;
  int    *vertOldTab, *permNodTab, nereal, npreal;
  int    *adja,iadr;


  /* Computing the number of vertices and a contiguous tabular of vertices */
  vertNbr = 0;

  _MMG5_ADD_MEM(mesh,(mesh->ne+1)*sizeof(int),"vertOldTab",return(1));
  _MMG5_SAFE_CALLOC(vertOldTab,mesh->ne+1,int);

  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {

    /* Testing if the tetra exists */
    if (!mesh->tetra[tetraIdx].v[0]) continue;
    vertOldTab[tetraIdx] = ++vertNbr;
  }

  if ( vertNbr/2 < _MMG5_BOXSIZE ) {
    /* not enough tetra to renum */
    _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
    return(1);
  }
  /* Allocating memory to compute adjacency lists */
  _MMG5_ADD_MEM(mesh,(vertNbr+2)*sizeof(SCOTCH_Num),"vertTab",
                _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                return(1));
  _MMG5_SAFE_CALLOC(vertTab,vertNbr+2,SCOTCH_Num);

  if (!memset(vertTab, ~0, sizeof(SCOTCH_Num)*(vertNbr + 2))) {
    perror("  ## Memory problem: memset");
    _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
    _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
    return 1;
  }

  edgeNbr = 1;
  /* Euler-Poincare formulae edgeSiz = 20*mesh->np~4*mesh->ne;
     (2*(12*mesh->np (triangles)-2*mesh->np (boundary triangles))) */
  edgeSiz = vertNbr*4;

  _MMG5_ADD_MEM(mesh,edgeSiz*sizeof(SCOTCH_Num),"edgeTab",
                _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                return(1));
  _MMG5_SAFE_CALLOC(edgeTab,edgeSiz,SCOTCH_Num);


  /* Computing the adjacency list for each vertex */
  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {

    /* Testing if the tetra exists */
    if (!mesh->tetra[tetraIdx].v[0]) continue;

    iadr = 4*(tetraIdx-1) + 1;
    adja = &mesh->adja[iadr];
    for (i=0; i<4; i++) {
      ballTetIdx = adja[i];

      if (!ballTetIdx) continue;

      ballTetIdx /= 4;

      /* Testing if one neighbour of tetraIdx has already been added */
      if (vertTab[vertOldTab[tetraIdx]] < 0)
        vertTab[vertOldTab[tetraIdx]] = edgeNbr;

      /* Testing if edgeTab memory is enough */
      if (edgeNbr >= edgeSiz) {
        _MMG5_ADD_MEM(mesh,0.2*sizeof(SCOTCH_Num),"edgeTab",
                      _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                      _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                      return(1));
        edgeSiz *= 1.2;
        _MMG5_SAFE_REALLOC(edgeTab,edgeSiz,SCOTCH_Num,"scotch table");
      }

      edgeTab[edgeNbr++] = vertOldTab[ballTetIdx];
    }
  }
  vertTab[vertNbr+1] = edgeNbr;
  edgeNbr--;
  /* Check if some tetras are alone */
  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {

    /* Testing if the tetra exists */
    if (!mesh->tetra[tetraIdx].v[0]) continue;
    if (vertTab[vertOldTab[tetraIdx]] < 0) {
      if(vertOldTab[tetraIdx] == vertNbr) {
        fprintf(stdout,"WARNING graph problem, no renum\n");
        _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
        _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
        return(1);
      }
      if(vertTab[vertOldTab[tetraIdx] + 1] > 0)
        vertTab[vertOldTab[tetraIdx]] = vertTab[vertOldTab[tetraIdx] + 1];
      else {
        if(vertOldTab[tetraIdx]+1 == vertNbr) {
          fprintf(stdout,"WARNING graph problem, no renum\n");
          _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
          _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
          return(1);
        }
        i = 1;
        do  {
          i++;
        } while((vertTab[vertOldTab[tetraIdx] + i] < 0) && ((vertOldTab[tetraIdx] + i) < vertNbr));
        if(vertOldTab[tetraIdx] + i == vertNbr) {
          fprintf(stdout,"WARNING graph problem, no renum\n");
          _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
          _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
          return(1);
        }
        vertTab[vertOldTab[tetraIdx]] = vertTab[vertOldTab[tetraIdx] + i];
      }
    }
  }

  /* free adjacents to gain memory space */
  _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  /* Building the graph by calling Scotch functions */
  SCOTCH_graphInit(&graf) ;
  CHECK_SCOTCH(SCOTCH_graphBuild(&graf, (SCOTCH_Num) 1, vertNbr, vertTab+1,
                                 NULL, NULL, NULL, edgeNbr, edgeTab+1, NULL),
               "scotch_graphbuild", 0) ;
#ifndef NDEBUG
  /* don't check in release mode */
  if ( mesh->info.imprim > 6 || mesh->info.ddebug )
    fprintf(stdout,"** Checking scotch graph.\n");
  CHECK_SCOTCH(SCOTCH_graphCheck(&graf), "scotch_graphcheck", 0);
#endif

  _MMG5_ADD_MEM(mesh,(vertNbr+1)*sizeof(SCOTCH_Num),"permVrtTab",
                _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));
                _MMG5_DEL_MEM(mesh,vertTab,(vertNbr+2)*sizeof(SCOTCH_Num));
                _MMG5_DEL_MEM(mesh,edgeTab,edgeSiz*sizeof(SCOTCH_Num));
                if( !MMG3D_hashTetra(mesh,1) ) return(0);
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
                if( !MMG3D_hashTetra(mesh,1) ) return(0);
                return(1));
  _MMG5_SAFE_CALLOC(permNodTab,mesh->np+1,int);

  nereal = 0;
  npreal = 0;
  /* Create the final permutation table for tetras (stored in vertOldTab) */
  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {
    if ( !mesh->tetra[tetraIdx].v[0] )  continue;
    vertOldTab[tetraIdx] = permVrtTab[vertOldTab[tetraIdx]];
  }
  _MMG5_DEL_MEM(mesh,permVrtTab,(vertNbr+1)*sizeof(SCOTCH_Num));

  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {
    while ( vertOldTab[tetraIdx] != tetraIdx && vertOldTab[tetraIdx] )
      _MMG5_swapTet(mesh->tetra/*,mesh->adja*/,vertOldTab,tetraIdx,vertOldTab[tetraIdx]);
  }

  for(tetraIdx = 1 ; tetraIdx < mesh->ne + 1 ; tetraIdx++) {
    ptet = &mesh->tetra[tetraIdx];

    /* Testing if the tetra exists */
    if (!ptet->v[0]) continue;

    nereal++;

    for(j = 0 ; j <= 3 ; j++) {

      nodeGlbIdx = ptet->v[j];

      if ( permNodTab[nodeGlbIdx] ) continue;

      ppt = &mesh->point[nodeGlbIdx];

      if ( !(ppt->tag & MG_NUL) )
        /* Building the new point list */
        permNodTab[nodeGlbIdx] = ++npreal;
    }
  }
  for( k = 1; k <= mesh->nprism; ++k) {
    pp = &mesh->prism[k];

    if ( !MG_EOK(pp) ) continue;

    for(j = 0; j<6 ; j++) {
      nodeGlbIdx = pp->v[j];

      if ( permNodTab[nodeGlbIdx] ) continue;

      ppt = &mesh->point[nodeGlbIdx];

      if ( !(ppt->tag & MG_NUL) )
        /* Building the new point list */
        permNodTab[nodeGlbIdx] = ++npreal;
    }
  }

  _MMG5_DEL_MEM(mesh,vertOldTab,(mesh->ne+1)*sizeof(int));

  /* Modify the numbering of the nodes of each tetra */
  for( tetraIdx = 1; tetraIdx < nereal + 1; tetraIdx++) {
    for(j = 0 ; j <= 3 ; j++) {
      mesh->tetra[tetraIdx].v[j] = permNodTab[mesh->tetra[tetraIdx].v[j]];
    }
  }

  /* Modify the numbering of the nodes of each prism */
  for( k = 1; k <= mesh->nprism; ++k) {
    for(j = 0; j<6 ; j++) {
      mesh->prism[k].v[j] = permNodTab[mesh->prism[k].v[j]];
    }
  }

  /* Modify the numbering of the nodes of each quad */
  for( k = 1; k <= mesh->nquad; ++k) {
    for(j = 0; j<4 ; j++) {
      mesh->quad[k].v[j] = permNodTab[mesh->quad[k].v[j]];
    }
  }

  /* Permute nodes and sol */
  for (j=1; j<= mesh->np; j++) {
    while ( permNodTab[j] != j && permNodTab[j] )
      _MMG5_swapNod(mesh->point,sol->m,permNodTab,j,permNodTab[j],sol->size);
  }
  _MMG5_DEL_MEM(mesh,permNodTab,(mesh->np+1)*sizeof(int));


  mesh->ne = nereal;
  mesh->np = npreal;

  if ( mesh->np == mesh->npmax )
    mesh->npnil = 0;
  else
    mesh->npnil = mesh->np + 1;

  if ( mesh->ne == mesh->nemax )
    mesh->nenil = 0;
  else
    mesh->nenil = mesh->ne + 1;

  if ( mesh->npnil )
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;

  if ( mesh->nenil )
    for (k=mesh->nenil; k<mesh->nemax-1; k++)
      mesh->tetra[k].v[3] = k+1;

  if( !MMG3D_hashTetra(mesh,0) ) return(0);

  return 1;
}
#endif

