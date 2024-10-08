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
 * Example of use of the mmg3d library (migrate from the mmg3d4 to the mmg3d
 * library).
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
#include "mmg/mmg3d/libmmg3d.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  int             ier,k;
  char            *outname;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");
  if ( argc != 2 ) {
    printf(" Usage: %s  fileout \n",argv[0]);
    return(1);
  }

  /* Name and path of the mesh files */
  outname = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( outname == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(outname,argv[1]);

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh:
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */

  mmgMesh = NULL;
  mmgSol  = NULL;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG3D_Set* functions */

  /** Manually set of the mesh */
  /** a) give the size of the mesh: 12 vertices, 12 tetra,0 prisms, 20
   * triangles, 0 quads, 0 edges */
  if ( MMG3D_Set_meshSize(mmgMesh,12,12,0,20,0,0) != 1 )  exit(EXIT_FAILURE);

  /** b) give the vertices: for each vertex, give the coordinates, the reference
      and the position in mesh of the vertex */
  mmgMesh->point[1].c[0]  = 0.;  mmgMesh->point[1].c[1]  = 0.; mmgMesh->point[1].c[2]  = 0.; mmgMesh->point[1].ref  = 0;
  /* or with the api function :
     if ( MMG3D_Set_vertex(mmgMesh,0  ,0  ,0  ,0,  1) != 1 )  exit(EXIT_FAILURE); */
  mmgMesh->point[2].c[0]  = 0.5; mmgMesh->point[2].c[1]  = 0;  mmgMesh->point[2].c[2]  = 0;  mmgMesh->point[2].ref  = 0;
  mmgMesh->point[3].c[0]  = 0.5; mmgMesh->point[3].c[1]  = 0;  mmgMesh->point[3].c[2]  = 1;  mmgMesh->point[3].ref  = 0;
  mmgMesh->point[4].c[0]  = 0;   mmgMesh->point[4].c[1]  = 0;  mmgMesh->point[4].c[2]  = 1;  mmgMesh->point[4].ref  = 0;
  mmgMesh->point[5].c[0]  = 0;   mmgMesh->point[5].c[1]  = 1;  mmgMesh->point[5].c[2]  = 0;  mmgMesh->point[5].ref  = 0;
  mmgMesh->point[6].c[0]  = 0.5; mmgMesh->point[6].c[1]  = 1;  mmgMesh->point[6].c[2]  = 0;  mmgMesh->point[6].ref  = 0;
  mmgMesh->point[7].c[0]  = 0.5; mmgMesh->point[7].c[1]  = 1;  mmgMesh->point[7].c[2]  = 1;  mmgMesh->point[7].ref  = 0;
  mmgMesh->point[8].c[0]  = 0;   mmgMesh->point[8].c[1]  = 1;  mmgMesh->point[8].c[2]  = 1;  mmgMesh->point[8].ref  = 0;
  mmgMesh->point[9].c[0]  = 1;   mmgMesh->point[9].c[1]  = 0;  mmgMesh->point[9].c[2]  = 0;  mmgMesh->point[9].ref  = 0;
  mmgMesh->point[10].c[0] = 1;   mmgMesh->point[10].c[1] = 1;  mmgMesh->point[10].c[2] = 0;  mmgMesh->point[10].ref = 0;
  mmgMesh->point[11].c[0] = 1;   mmgMesh->point[11].c[1] = 0;  mmgMesh->point[11].c[2] = 1;  mmgMesh->point[11].ref = 0;
  mmgMesh->point[12].c[0] = 1;   mmgMesh->point[12].c[1] = 1;  mmgMesh->point[12].c[2] = 1;  mmgMesh->point[12].ref = 0;

  /*tetra*/
  mmgMesh->tetra[1].v[0]  = 1;  mmgMesh->tetra[1].v[1]  = 2;  mmgMesh->tetra[1].v[2]  = 4;  mmgMesh->tetra[1].v[3]  = 8;  mmgMesh->tetra[1].ref  = 1;
  /* or with the api function :
     if ( MMG3D_Set_tetrahedra(mmgMesh,1 ,2 ,4 ,8, 1) != 1 )  exit(EXIT_FAILURE); */
  mmgMesh->tetra[2].v[0]  = 8;  mmgMesh->tetra[2].v[1]  = 3;  mmgMesh->tetra[2].v[2]  = 2;  mmgMesh->tetra[2].v[3]  = 7;  mmgMesh->tetra[2].ref  = 1;
  mmgMesh->tetra[3].v[0]  = 2;  mmgMesh->tetra[3].v[1]  = 5;  mmgMesh->tetra[3].v[2]  = 6;  mmgMesh->tetra[3].v[3]  = 8;  mmgMesh->tetra[3].ref  = 1;
  mmgMesh->tetra[4].v[0]  = 8;  mmgMesh->tetra[4].v[1]  = 5;  mmgMesh->tetra[4].v[2]  = 1;  mmgMesh->tetra[4].v[3]  = 2;  mmgMesh->tetra[4].ref  = 1;
  mmgMesh->tetra[5].v[0]  = 2;  mmgMesh->tetra[5].v[1]  = 7;  mmgMesh->tetra[5].v[2]  = 8;  mmgMesh->tetra[5].v[3]  = 6;  mmgMesh->tetra[5].ref  = 1;
  mmgMesh->tetra[6].v[0]  = 2;  mmgMesh->tetra[6].v[1]  = 4;  mmgMesh->tetra[6].v[2]  = 3;  mmgMesh->tetra[6].v[3]  = 8;  mmgMesh->tetra[6].ref  = 1;
  mmgMesh->tetra[7].v[0]  = 2;  mmgMesh->tetra[7].v[1]  = 9;  mmgMesh->tetra[7].v[2]  = 3;  mmgMesh->tetra[7].v[3]  = 7;  mmgMesh->tetra[7].ref  = 2;
  mmgMesh->tetra[8].v[0]  = 7;  mmgMesh->tetra[8].v[1]  = 11; mmgMesh->tetra[8].v[2]  = 9;  mmgMesh->tetra[8].v[3]  = 12; mmgMesh->tetra[8].ref  = 2;
  mmgMesh->tetra[9].v[0]  = 9;  mmgMesh->tetra[9].v[1]  = 6;  mmgMesh->tetra[9].v[2]  = 10; mmgMesh->tetra[9].v[3]  = 7;  mmgMesh->tetra[9].ref  = 2;
  mmgMesh->tetra[10].v[0] = 7;  mmgMesh->tetra[10].v[1] = 6;  mmgMesh->tetra[10].v[2] = 2;  mmgMesh->tetra[10].v[3] = 9;  mmgMesh->tetra[10].ref = 2;
  mmgMesh->tetra[11].v[0] = 9;  mmgMesh->tetra[11].v[1] = 12; mmgMesh->tetra[11].v[2] = 7;  mmgMesh->tetra[11].v[3] = 10; mmgMesh->tetra[11].ref = 2;
  mmgMesh->tetra[12].v[0] = 9;  mmgMesh->tetra[12].v[1] = 3;  mmgMesh->tetra[12].v[2] = 11; mmgMesh->tetra[12].v[3] = 7;  mmgMesh->tetra[12].ref = 2;

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG3D_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG3D_Set* functions */

  /** Manually set of the sol */
  /** a) give info for the sol structure: sol applied on vertex entities,
      number of vertices=12, the sol is scalar*/
  if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,12,MMG5_Scalar) != 1 )
    exit(EXIT_FAILURE);

  /** b) give solutions values and positions */
  for(k=1 ; k<=12 ; k++) {
    mmgSol->m[k] = 0.5;
    /* or with the api function :
       if ( MMG3D_Set_scalarSol(mmgSol,0.5,k) != 1 ) exit(EXIT_FAILURE); */
  }
  /** 4) If you don't use the API functions, you MUST call
      the MMG3D_Set_handGivenMesh() function. Don't call it if you use
      the API functions */
  MMG3D_Set_handGivenMesh(mmgMesh);

  /** 5) (not mandatory): check if the number of given entities match with mesh size */
  if ( MMG3D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

  /** ------------------------------ STEP  II -------------------------- */
  /** remesh function */
  // WARNING: the MMG3D_mmg3dlib function returns 1 if success, 0 if fail.
  // The MMG3D4 library was working opposite.
  ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");


  /** ------------------------------ STEP III -------------------------- */
  /** get results */
  /** Two solutions: just use the MMG3D_saveMesh/MMG3D_saveSol functions
      that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
      using the MMG3D_getMesh/MMG3D_getSol functions */

  /** 1) Automatically save the mesh */
  if ( MMG3D_saveMesh(mmgMesh,outname) != 1 ) exit(EXIT_FAILURE);

  /** 2) Automatically save the solution */
  if ( MMG3D_saveSol(mmgMesh,mmgSol,outname) !=1 ) exit(EXIT_FAILURE);

  /** 3) Free the MMG3D5 structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);
  free(outname);
  outname = NULL;

  return(ier);
}
