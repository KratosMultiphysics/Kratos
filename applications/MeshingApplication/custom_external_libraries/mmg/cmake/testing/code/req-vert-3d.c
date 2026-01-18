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
 * Test preservation of required vertex: vertex 1 is required and shoule not move.
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
  MMG5_pMesh      mesh;
  MMG5_pSol       sol;
  char            *file;
  double          c[3];
  int             k,ier;

  fprintf(stdout,"  -- CHECK PRESERVATION OF REQUIRED VERTICES \n");

  if ( argc != 2 ) {
    printf(" Usage: %s filein \n",argv[0]);
    return(1);
  }

  /* Name and path of the mesh file */
  file = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( file == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(file,argv[1]);

  /** Read mesh */
  mesh = NULL;
  sol  = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh,
                  MMG5_ARG_ppMet,&sol,
                  MMG5_ARG_end);

  if ( MMG3D_loadMesh(mesh,file) != 1 )  {
    fprintf(stderr,"Error: %s: %d: File not found %s\n.",__func__,__LINE__,file);
    exit(EXIT_FAILURE);
  }

  /** Check that vertex number 1 is required and store its coordinates */
  if ( !(mesh->point[1].tag & MG_REQ) ) {
    fprintf(stderr,"Error: %s: %d: This test expects that vertex of index 1 is required\n.",__func__,__LINE__);
    exit(EXIT_FAILURE);
  }
  c[0] = mesh->point[1].c[0];
  c[1] = mesh->point[1].c[1];
  c[2] = mesh->point[1].c[2];

  /** Enable vertex regularisation */
  if ( MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_xreg,1) != 1 )
    exit(EXIT_FAILURE);

  /** remesh function */
  ier = MMG3D_mmg3dlib(mesh,sol);


  /** Check that coordinates of vertex 1 (that is required) have not changed. */
  double dd[3];
  dd[0] =  c[0] -  mesh->point[1].c[0];
  dd[1] =  c[1] -  mesh->point[1].c[1];
  dd[2] =  c[2] -  mesh->point[1].c[2];


  int j;
  for (j=0; j<3; ++j) {
    printf("%.15lf %.15lf\n",c[j],mesh->point[1].c[j]);
    if ( fabs(c[j]- mesh->point[1].c[j]) > 1e-5 ) {
      fprintf(stderr,"Error: %s: %d:"
              " Modification of coordinates of vertex 1 (required):"
              "   input  coor: %15lf %15lf %15lf\n"
              "   output coor: %15lf %15lf %15lf\n",__func__,__LINE__,
                  c[0],c[1],c[2], mesh->point[1].c[0], mesh->point[1].c[1], mesh->point[1].c[2]);
      exit(EXIT_FAILURE);
    }
  }

  fprintf(stdout,"MMG3D: REQUIRED VERTEX SUCCESFULLY PRESERVED.\n");

  /** 3) Free the MMG3D5 structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh,
                 MMG5_ARG_ppMet,&sol,
                 MMG5_ARG_end);

  free(file);
  file = NULL;

  return 0;
}
