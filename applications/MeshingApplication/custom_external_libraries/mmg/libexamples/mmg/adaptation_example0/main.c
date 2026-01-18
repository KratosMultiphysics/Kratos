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
 * Example of use of the mmg library that gathers the mmg2d, mmgs and mmg3d
 * libraries (basic use of mesh adaptation).
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

/** Include the mmg library hader file */
// if the "include/mmg" dir is in your include path
//#include "libmmg.h"
// if your include path do not contain the "mmg" subdirectories
#include "mmg/libmmg.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  int             ier;
  char            *filename,*filename_o2d,*filename_os,*filename_o3d,*ptr;

  fprintf(stdout,"  -- TEST MMGLIB \n");

  if ( argc != 4 ) {
    printf(" Usage: %s 2d_filein 3d_filein fileout\n",argv[0]);
    return(1);
  }

  /** ================== 2d remeshing using the mmg2d library ========== */

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

 /* Name and path of the mesh file */

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
  MMG2D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh, MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG2D_Set* functions */

  /** with MMG2D_loadMesh function */
  if ( MMG2D_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);
  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG2D_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG2D_Set* functions */
  if ( MMG2D_loadSol(mmgMesh,mmgSol,filename) != 1 )  exit(EXIT_FAILURE);

  /** ------------------------------ STEP  II -------------------------- */
  /** remesh function */
  ier = MMG2D_mmg2dlib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG2DLIB\n");

  /** ------------------------------ STEP III -------------------------- */
  /** get results */
  /** Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
      that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
      using the MMG2D_getMesh/MMG2D_getSol functions */

  /** 1) Automatically save the mesh */
  /*save result*/
  filename_o2d = (char *) calloc(strlen(argv[3]) + 4, sizeof(char));
  if ( filename_o2d == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename_o2d,argv[3]);
  ptr = strstr(filename_o2d,".mesh");
  if ( !ptr )  ptr = strstr(filename_o2d,".msh");
  if ( ptr ) *ptr = '\0';
  strcat(filename_o2d,".2d");

  if ( MMG2D_saveMesh(mmgMesh,filename_o2d) != 1 )  exit(EXIT_FAILURE);
  /*save metric*/
  if ( MMG2D_saveSol(mmgMesh,mmgSol,filename_o2d) != 1 )  exit(EXIT_FAILURE);

  /** 2) Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh, MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free(filename);
  filename = NULL;
  /** ================ surface remeshing using the mmgs library ======== */

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[2]);

  filename_os = (char *) calloc(strlen(argv[3]) + 3, sizeof(char));
  if ( filename_os == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename_os,argv[3]);
  ptr = strstr(filename_os,".mesh");
  if ( !ptr )  ptr = strstr(filename_os,".msh");
  if ( ptr ) *ptr = '\0';
  strcat(filename_os,".s");

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */

  mmgMesh = NULL;
  mmgSol  = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh, MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMGS_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMGS_Set* functions */

  /** with MMGS_loadMesh function */
  /** function calling */
  if ( MMGS_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);
  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMGS_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMGS_Set* functions */

  /** With MMGS_loadSol function */
  if ( MMGS_loadSol(mmgMesh,mmgSol,filename) != 1 )  exit(EXIT_FAILURE);

  /** ------------------------------ STEP  II -------------------------- */
  /** remesh function */
  ier = MMGS_mmgslib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMGSLIB\n");

  /** ------------------------------ STEP III -------------------------- */
  /** get results */
  /** Two solutions: just use the MMGS_saveMesh/MMGS_saveSol functions
      that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
      using the MMGS_getMesh/MMGS_getSol functions */

  /** 1) Automatically save the mesh */
  if ( MMGS_saveMesh(mmgMesh,filename_os) != 1 )  exit(EXIT_FAILURE);

  /** 2) Automatically save the solution */
  if ( MMGS_saveSol(mmgMesh,mmgSol,filename_os) != 1 )  exit(EXIT_FAILURE);

  /** 3) Free the MMGS structures */
  MMGS_Free_all(MMG5_ARG_start,
                MMG5_ARG_ppMesh,&mmgMesh, MMG5_ARG_ppMet,&mmgSol,
                MMG5_ARG_end);

  /** ================== 3d remeshing using the mmg3d library ========== */
  filename_o3d = (char *) calloc(strlen(argv[3]) + 4, sizeof(char));
  if ( filename_o3d == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename_o3d,argv[3]);
  ptr = strstr(filename_o3d,".mesh");
  if ( !ptr )  ptr = strstr(filename_o3d,".msh");
  if ( ptr ) *ptr = '\0';
  strcat(filename_o3d,".3d");

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
                  MMG5_ARG_ppMesh,&mmgMesh, MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG3D_Set* functions */

  /** with MMG3D_loadMesh function */
  if ( MMG3D_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG3D_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG3D_Set* functions */

  /** With MMG3D_loadSol function */
  if ( MMG3D_loadSol(mmgMesh,mmgSol,filename) != 1 )  exit(EXIT_FAILURE);

  /** ------------------------------ STEP  II -------------------------- */
  /** remesh function */
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
  if ( MMG3D_saveMesh(mmgMesh,filename_o3d) != 1 )  exit(EXIT_FAILURE);

  /** 2) Automatically save the solution */
  if ( MMG3D_saveSol(mmgMesh,mmgSol,filename_o3d) != 1 )  exit(EXIT_FAILURE);

  /** 3) Free the MMG3D structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free(filename);
  free(filename_o2d);
  free(filename_os);
  free(filename_o3d);

  return(0);
}
