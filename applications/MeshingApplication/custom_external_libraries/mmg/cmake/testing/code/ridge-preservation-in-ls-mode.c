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
 * Example of use of the mmg3d library (basic use of mesh adaptation)
 *
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/** Include the mmg3d library hader file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "libmmg3d_private.h"
#include "mmg/mmg3d/libmmg3d.h"


int main(int argc,char *argv[]) {
  MMG5_pMesh      mesh1,mesh2;
  MMG5_pSol       sol1,sol2;
  char            *file1, *file2;

  fprintf(stdout,"  -- CHECK FOR PRESERVATION OF INPUT RIDGES IN LS MODE \n");

  if ( argc != 3 ) {
    printf(" Usage: %s filein fileout \n",argv[0]);
    printf("        filout has to come from the call of Mmg over the"
           " filein mesh with -noinsert -noswap -rn 0 arguments.\n");
    return(1);
  }

  /* Name and path of the mesh file */
  file1 = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( file1 == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(file1,argv[1]);

  file2 = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( file2 == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(file2,argv[2]);

  /** Read meshes */
  mesh1 = NULL;
  mesh2 = NULL;
  sol1  = NULL;
  sol2  = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh1,
                  MMG5_ARG_ppMet,&sol1,
                  MMG5_ARG_end);

  if ( MMG3D_loadMesh(mesh1,file1) != 1 )  {
    fprintf(stderr,"Error: %s: %d: File not found %s\n.",__func__,__LINE__,file1);
    exit(EXIT_FAILURE);
  }


  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh2,
                  MMG5_ARG_ppMet,&sol2,
                  MMG5_ARG_end);

  if ( MMG3D_loadMesh(mesh2,file2) != 1 )  {
    fprintf(stderr,"Error: %s: %d: File not found %s\n.",__func__,__LINE__,file2);
    exit(EXIT_FAILURE);
  }

  /** We want to check that all ridges of input mesh are still ridges in the
   * output of Mmg, that is, the second mesh of the current app) */

  /* Step 1: Hash edges of mesh2 and store their tags */
  MMG5_int k;
  MMG5_HGeom hash;
  MMG5_hNew(mesh2,&hash,mesh2->na,3*mesh2->na);

  for ( k=1; k<=mesh2->na; ++k ) {
    MMG5_pEdge ped = &mesh2->edge[k];
    if ( !MMG5_hEdge(mesh2,&hash,ped->a,ped->b,ped->ref,ped->tag) ) {
      fprintf(stderr,"Error: %s: %d: Unable to hash edge %" MMG5_PRId
              ": %" MMG5_PRId " %" MMG5_PRId ".\n",
              __func__,__LINE__,k,ped->a,ped->b);
      exit(EXIT_FAILURE);
    }
  }

  /* Step 2: Travel ridges of mesh1 and check that, if it is found in mesh2, it
   * is a ridge too (as level-set splitting may splits boudary edges, we can't
   * ensure that all ridges of input mesh will be present in the output mesh
   * (they may have been splitted) */
  MMG5_int ier = 0;
  for ( k=1; k<=mesh1->na; ++k ) {
    MMG5_pEdge ped = &mesh1->edge[k];
    if ( ped->tag & MG_GEO ) {
      MMG5_int ref;
      uint16_t tag;
      if ( !MMG5_hGet(&hash,ped->a,ped->b,&ref,&tag) ) {
        continue;
      }
      if ( ! (tag & MG_GEO) ) {
        /* ridge of mesh1 exists in mesh2 but is not ridge anymore */
        fprintf(stderr,"Error: %s: %d: Ridge %" MMG5_PRId
                " (%" MMG5_PRId " %" MMG5_PRId ") of first mesh is not"
                " ridge in second mesh.\n",__func__,__LINE__,k,ped->a,ped->b);
        ++ier;
      }
    }
  }

  if ( ier ) {
    fprintf(stderr,"Error: %s: %d: At least %" MMG5_PRId " missing ridges.\n",
            __func__,__LINE__,ier);
    exit(EXIT_FAILURE);
  }

  MMG5_SAFE_FREE(hash.geom);

  /** 3) Free the MMG3D5 structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh1,
                 MMG5_ARG_ppMet,&sol1,
                 MMG5_ARG_end);
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh2,
                 MMG5_ARG_ppMet,&sol2,
                 MMG5_ARG_end);

  free(file1);
  file1 = NULL;

  free(file2);
  file2 = NULL;

  return 0;
}
