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
 * \file common/librnbg.c
 * \brief Functions for scotch renumerotation.
 * \author Cedric Lachat  (Inria/UBordeaux)
 * \author Algiane Froehly  (Inria/UBordeaux)
 * \version 5
 * \date 2013
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon.h"

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
 * k-partitioning and assuming that baseval of the graph is 1.
 *
 **/
int _MMG5_kPartBoxCompute(SCOTCH_Graph graf, int vertNbr, int boxVertNbr,
                          SCOTCH_Num *permVrtTab,MMG5_pMesh mesh) {
  int boxNbr, vertIdx;
#ifdef SCOTCH_5
  SCOTCH_Num logMaxVal, SupMaxVal, InfMaxVal, maxVal;
#endif
  char s[200];
  SCOTCH_Num *sortPartTb;
  SCOTCH_Strat strat ;
  SCOTCH_Arch arch;

  /* Computing the number of boxes */
  boxNbr = vertNbr / boxVertNbr;
  if (boxNbr * boxVertNbr != vertNbr) {
    boxNbr = boxNbr + 1;
  }


  /* Initializing SCOTCH functions */
  CHECK_SCOTCH(SCOTCH_stratInit(&strat), "scotch_stratInit", 0) ;
#ifdef SCOTCH_6
  CHECK_SCOTCH(SCOTCH_archCmplt(&arch, boxNbr), "scotch_archCmplt", 0) ;
#else
  CHECK_SCOTCH(SCOTCH_archVcmplt(&arch), "scotch_archVcmplt", 0) ;
#endif
  sprintf(s, "m{vert=%d,low=r{job=t,map=t,poli=S,sep=m{vert=80,low=h{pass=10}f{bal=0.0005,move=80},asc=f{bal=0.005,move=80}}}}", vertNbr / boxVertNbr);
  CHECK_SCOTCH(SCOTCH_stratGraphMap(&strat, s), "scotch_stratGraphMap", 0) ;

  _MMG5_ADD_MEM(mesh,2*vertNbr*sizeof(SCOTCH_Num),"sortPartTb",return(1));
  _MMG5_SAFE_CALLOC(sortPartTb,2*vertNbr,SCOTCH_Num);

  /* Partionning the graph */
  CHECK_SCOTCH(SCOTCH_graphMap(&graf, &arch, &strat, sortPartTb), "scotch_graphMap", 0);


#ifdef SCOTCH_6
  // Looking for the max value in sortPartTb and computing sortPartTb as
  // followed :
  //  - sortPartTb[2i] is the box value
  //  - sortPartTb[2i+1] is the vertex number
  maxVal = sortPartTb[0];
#endif
  for (vertIdx = vertNbr - 1 ; vertIdx >= 0 ; vertIdx--) {
    sortPartTb[2*vertIdx] = sortPartTb[vertIdx];
    sortPartTb[2*vertIdx+1] = vertIdx + 1;
#ifdef SCOTCH_5
    if (sortPartTb[vertIdx] > maxVal)
      maxVal = sortPartTb[vertIdx];
#endif
  }

#ifdef SCOTCH_5
  // Determining the log of MaxVal
  logMaxVal = 0;
  while ( maxVal > 0) {
    logMaxVal++;
    maxVal >>= 1;
  }

  // Infering the interval in which box values will be
  InfMaxVal = logMaxVal << logMaxVal;
  SupMaxVal = (logMaxVal << (logMaxVal + 1)) - 1;

  // Increasing box values until they are in the previous interval
  for (vertIdx = 0 ; vertIdx < vertNbr ; vertIdx++) {
    while (!(sortPartTb[2*vertIdx] >= InfMaxVal && sortPartTb[2*vertIdx] <= SupMaxVal)) {
      sortPartTb[2*vertIdx] <<= 1;
    }
  }
#endif

  // Sorting the tabular, which contains box values and vertex numbers
  _SCOTCHintSort2asc1(sortPartTb, vertNbr);


  /* Infering the new numbering */
  for (vertIdx = 0; vertIdx < vertNbr ; vertIdx++) {
    permVrtTab[sortPartTb[2*vertIdx + 1]] = vertIdx + 1;
  }

  SCOTCH_stratExit(&strat) ;
  SCOTCH_archExit(&arch) ;

  _MMG5_DEL_MEM(mesh,sortPartTb,2*vertNbr*sizeof(SCOTCH_Num));

  return 0;
}

/**
 * \param points pointer toward a table containing the point structures.
 * \param sols pointer toward a table containing the solution structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first tetra to swap.
 * \param ind2 index of the second tetra to swap.
 * \param solsize size of the solution.
 *
 * Swap two nodes in the table of vertices.
 *
 */
void _MMG5_swapNod(MMG5_pPoint points, double* sols, int* perm,
                   int ind1, int ind2, int solsiz) {
  MMG5_Point ptttmp;
  MMG5_Sol   soltmp;
  int        tmp,addr2,addr1;

  /* swap the points */
  memcpy(&ptttmp      ,&points[ind2],sizeof(MMG5_Point));
  memcpy(&points[ind2],&points[ind1],sizeof(MMG5_Point));
  memcpy(&points[ind1],&ptttmp      ,sizeof(MMG5_Point));

  /* swap the sols */
  if ( sols ) {
    addr1 = ind1*solsiz;
    addr2 = ind2*solsiz;
    memcpy(&soltmp     ,&sols[addr2],solsiz*sizeof(double));
    memcpy(&sols[addr2],&sols[addr1],solsiz*sizeof(double));
    memcpy(&sols[addr1],&soltmp     ,solsiz*sizeof(double));
  }

  /* swap the permutaion table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}
#endif

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 * \return 0 if \a _MMG5_renumbering fail (non conformal mesh), 1 otherwise
 * (renumerotation success of renumerotation fail but the mesh is still
 *  conformal).
 *
 * Call scotch renumbering.
 *
 **/
int _MMG5_scotchCall(MMG5_pMesh mesh, MMG5_pSol met)
{
#ifdef USE_SCOTCH
  /*check enough vertex to renum*/
  if ( mesh->info.renum && (mesh->np/2. > _MMG5_BOXSIZE) && mesh->np>100000 ) {
    /* renumbering begin */
    if ( mesh->info.imprim > 5 )
      fprintf(stdout,"  -- RENUMBERING. \n");

    if ( !_MMG5_renumbering(_MMG5_BOXSIZE,mesh, met) ) {
      fprintf(stdout,"  ## Unable to renumbering mesh. \n");
      fprintf(stdout,"  ## Try to run without renumbering option (-rn 0)\n");
      return(0);
    }

    if ( mesh->info.imprim > 5) {
      fprintf(stdout,"  -- PHASE RENUMBERING COMPLETED. \n");
    }

    if ( mesh->info.ddebug )  _MMG5_chkmsh(mesh,1,0);
    /* renumbering end */
  }
  return(1);
#else
  return(1);
#endif
}
