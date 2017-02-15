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
 * Example of use of the mmg3d library (basic use of lagrangian motion option)
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
  MMG5_pSol       mmgSol,mmgDisp;
  int             k,ier;
  char            *pwd,*inname,*outname;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  /* Name and path of the mesh files */
  pwd = getenv("PWD");
  inname = (char *) calloc(strlen(pwd) + 40, sizeof(char));
  if ( inname == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  outname = (char *) calloc(strlen(pwd) + 49, sizeof(char));
  if ( outname == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  sprintf(inname, "%s%s%s", pwd, "/../libexamples/mmg3d/example4/", "tinyBoxt");

  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh: mesh=&mmgMesh, sol=&mmgSol */
  mmgMesh = NULL;
  mmgSol  = NULL;
  mmgDisp = NULL; //Useless here: just needed forthe lagrangian motion option
  MMG5_Init_mesh(&mmgMesh,&mmgSol,&mmgDisp);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG5_loadMesh function that will read a .mesh(b)
     file formatted or manually set your mesh using the MMG5_Set* functions */

  /** with MMG5_loadMesh function */
  /** a) (not mandatory): give the mesh name
     (by default, the "mesh.mesh" file is oppened)*/
  if ( !MMG5_Set_inputMeshName(mmgMesh,inname) )
    exit(EXIT_FAILURE);
  /** b) function calling */
  if ( !MMG5_loadMesh(mmgMesh) )  exit(EXIT_FAILURE);

  /** 3) Build displacement in MMG5 format */
  /** Two solutions: just use the MMG5_loadMet function that will read a .sol(b)
      file formatted or manually set your sol using the MMG5_Set* functions */

  /**------------------- Lagrangian motion option ----------------------------*/
  /* Ask for lagrangian motion (mode 1) */
  if ( !MMG5_Set_iparameter(mmgMesh,mmgDisp,MMG5_IPARAM_lag, 1) )
    exit(EXIT_FAILURE);

  /** With MMG5_loadMet function */
  /** a) (not mandatory): give the sol name
     (by default, the "mesh.sol" file is oppened)*/
  if ( !MMG5_Set_inputSolName(mmgMesh,mmgDisp,inname) )
    exit(EXIT_FAILURE);

  /** b) function calling */
  if ( !MMG5_loadMet(mmgMesh,mmgDisp) )
    exit(EXIT_FAILURE);

  /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( !MMG5_Chk_meshData(mmgMesh,mmgDisp) ) exit(EXIT_FAILURE);

  /** 5) (not mandatory): set your global parameters using the
      MMG5_Set_iparameter and MMG5_Set_dparameter function
      (resp. for integer parameters and double param)*/


  /**------------------- Lagrangian motion computation ---------------------*/

  /* debug mode ON (default value = OFF) */
  if ( !MMG5_Set_iparameter(mmgMesh,mmgDisp,MMG5_IPARAM_debug, 1) )
    exit(EXIT_FAILURE);

  /** remesh function */
  ier = MMG5_mmg3dmov(mmgMesh,mmgSol,mmgDisp);
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /* (Not mandatory) Automatically save the mesh */
  sprintf(outname, "%s%s%s", pwd, "/../libexamples/mmg3d/example4/", "tinyBoxt.o.mesh");
  if ( !MMG5_Set_outputMeshName(mmgMesh,outname) )
    exit(EXIT_FAILURE);

  MMG5_saveMesh(mmgMesh);

  /* (Not mandatory) Automatically save the solution */
  if ( !MMG5_Set_outputSolName(mmgMesh,mmgSol,outname) )
    exit(EXIT_FAILURE);

  MMG5_saveMet(mmgMesh,mmgSol);

  /* 9) free the MMG3D5 structures */
  MMG5_Free_all(mmgMesh,mmgSol,mmgDisp);

  free(inname);
  inname = NULL;
  free(outname);
  outname = NULL;

  return(ier);
}
