/*Authors Cécile Dobrzynski

  Example for using mmg2dlib

*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/** Include the mmg2d library hader file */
// if the header file is in the "include" directory
// #include "libmmg2d.h"
// if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2d.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  char            *filename, *outname;

  int             ier;

  fprintf(stdout,"  -- TEST MMG2DMESH \n");


  if ( argc != 3 ) {
    printf(" Usage: %s filein fileout \n",argv[0]);
    return(1);
  }

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  outname = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( outname == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(outname,argv[2]);

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
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG2D_Set* functions */

  /** read the mesh in a mesh file */
  MMG2D_loadMesh(mmgMesh,filename);

  /** Set parameters : for example set the maximal edge size to 0.1 */
  MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmax,0.1);

  /** Higher verbosity level */
  MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose,5);


  /** Generate the mesh */
  ier = MMG2D_mmg2dmesh(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG2DMESH: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG2DMESH\n");

  /*save result*/
  if ( MMG2D_saveMesh(mmgMesh,outname) != 1 )
    exit(EXIT_FAILURE);

  /*save metric*/
  if ( MMG2D_saveSol(mmgMesh,mmgSol,outname) != 1 )
    exit(EXIT_FAILURE);

  /** 3) Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free(filename);
  filename = NULL;


  free(outname);
  outname = NULL;

  return(0);
}
