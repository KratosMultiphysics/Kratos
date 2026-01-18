/** =============================================================================
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
**/

/**
 * Example of use of the mmg3d library (advanced use of mesh adaptation)
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
  int             k,ier;
  char            *inname, *outname1, *outname2;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  if ( argc != 4 ) {
    printf(" Usage: %s filein fileout1 filout2 \n",argv[0]);
    return(1);
  }

  /* Name and path of the mesh files */
  inname = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( inname == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(inname,argv[1]);

  outname1 = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( outname1 == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(outname1,argv[2]);

  outname2 = (char *) calloc(strlen(argv[3]) + 1, sizeof(char));
  if ( outname2 == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(outname2,argv[3]);

  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh:
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer to your MMG5_pMesh (that stores your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer to your MMG5_pSol (that stores your metric) */

  mmgMesh = NULL;
  mmgSol  = NULL;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
     file formatted or manually set your mesh using the MMG3D_Set* functions */

  /** with MMG3D_loadMesh function */
  if ( MMG3D_loadMesh(mmgMesh,inname) != 1 )  exit(EXIT_FAILURE);

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG3D_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG3D_Set* functions */

  /** With MMG3D_loadSol function */
  if ( MMG3D_loadSol(mmgMesh,mmgSol,inname) != 1 )
    exit(EXIT_FAILURE);

  /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( MMG3D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

  /** 5) (not mandatory): set your global parameters using the
      MMG3D_Set_iparameter and MMG3D_Set_dparameter function
      (resp. for integer parameters and double param)*/


  /**------------------- First wave of refinment---------------------*/

  /* debug mode ON (default value = OFF) */
  if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_debug, 1) != 1 )
    exit(EXIT_FAILURE);

  /* maximal memory size (default value = 50/100*ram) */
  if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_mem, 600) != 1 )
    exit(EXIT_FAILURE);

  /* Maximal mesh size (default FLT_MAX)*/
  if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmax,40) != 1 )
    exit(EXIT_FAILURE);

  /* Minimal mesh size (default 0)*/
  if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hmin,0.001) != 1 )
    exit(EXIT_FAILURE);

  /* Global hausdorff value (default value = 0.01) applied on the whole boundary */
  if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hausd, 0.1) != 1 )
    exit(EXIT_FAILURE);

  /* Gradation control*/
  if ( MMG3D_Set_dparameter(mmgMesh,mmgSol,MMG3D_DPARAM_hgrad, 2) != 1 )
    exit(EXIT_FAILURE);

  /** remesh function */
  ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /* (Not mandatory) Automatically save the mesh */
  MMG3D_saveMesh(mmgMesh,outname1);

  if ( MMG3D_saveSol(mmgMesh,mmgSol,outname1) != 1 )
    exit(EXIT_FAILURE);


  /**------------------- Second wave of refinment---------------------*/
  /* We add different local hausdorff numbers on boundary componants (this
     local values are used instead of the global hausdorff number) */

  /* verbosity (default value = 4)*/
  if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, 4) != 1 )
    exit(EXIT_FAILURE);

  if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_mem, 1000) != 1 )
    exit(EXIT_FAILURE);
  if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_debug, 0) != 1 )
    exit(EXIT_FAILURE);


  /** 6) (not mandatory): set your local parameters */
  /* use 2 local hausdorff numbers on ref 36 (hausd = 0.01) and 38 (hausd = 1) */
  if ( MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_numberOfLocalParam,2) != 1 )
    exit(EXIT_FAILURE);

  /** for each local parameter: give the type and reference of the element on which
      you will apply a particular hausdorff number and the hausdorff number
      to apply. The global hausdorff number is applied on all boundary triangles
      without local hausdorff number */

  /* Be careful if you change the hausdorff number (or gradation value)
     between 2 run: the information of the previous hausdorff number
     (resp. gradation) is contained in the metric computed during
     the previous run.
     Then, you can not grow up the hausdorff value (resp. gradation) without
     resetting this metric (but you can decrease this value). */

  if ( MMG3D_Set_localParameter(mmgMesh,mmgSol,MMG5_Triangle,36,0.00001,40.,0.01) != 1 )
    exit(EXIT_FAILURE);
  if ( MMG3D_Set_localParameter(mmgMesh,mmgSol,MMG5_Triangle,38,0.00001,40.,1) != 1 )
    exit(EXIT_FAILURE);

  /** remesh function */
  ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");

  /**------------------- Last wave of refinment---------------------*/
  /** remesh function */

  /* New metric to see the effect of the local hausdorff number on triangles
     of ref 38: constant and of size 10 */
  for ( k=1; k<=mmgSol->np; k++ ) {
    if ( MMG3D_Set_scalarSol(mmgSol,10,k) != 1 ) exit(EXIT_FAILURE);
  }

  ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");


  /* 7) Automatically save the mesh */
  if ( MMG3D_saveMesh(mmgMesh,outname2) != 1 )
    exit(EXIT_FAILURE);

  /* 8) Automatically save the solution */
  if ( MMG3D_saveSol(mmgMesh,mmgSol,outname2) != 1 )
    exit(EXIT_FAILURE);


  /* 9) free the MMG3D5 structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free(inname);
  inname = NULL;

  free(outname1);
  outname1 = NULL;

  free(outname2);
  outname2 = NULL;

  return(ier);
}
