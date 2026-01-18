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
  MMG5_int        k,np;
  int             ier;

  fprintf(stdout,"  -- TEST MMG2DLIB \n");

  if ( argc != 3 ) {
    printf(" Usage: %s filein  fileout\n",argv[0]);
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

  /** with MMG2D_loadMesh function */
  if ( MMG2D_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);

  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG2D_loadMet function that will read a .sol(b)
      file formatted or manually set your sol using the MMG2D_Set* functions */

  /** Manually set of the sol */
  /** a) Get np the number of vertex */
  if ( MMG2D_Get_meshSize(mmgMesh,&np,NULL,NULL,NULL) != 1 )
    exit(EXIT_FAILURE);

  /** b) give info for the sol structure: sol applied on vertex entities,
      number of vertices=np, the sol is scalar*/
  if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,np,MMG5_Scalar) != 1 )
    exit(EXIT_FAILURE);

  /** c) give solutions values and positions */
  for(k=1 ; k<=np ; k++) {
    if ( MMG2D_Set_scalarSol(mmgSol,0.01,k) != 1 ) exit(EXIT_FAILURE);
  }

  /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( MMG2D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

  ier = MMG2D_mmg2dlib(mmgMesh,mmgSol);

  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG2DLIB\n");

  /*save result*/
  if ( MMG2D_saveMesh(mmgMesh,outname) != 1 )  exit(EXIT_FAILURE);

  /*save metric*/
  if ( MMG2D_saveSol(mmgMesh,mmgSol,outname) != 1 )  exit(EXIT_FAILURE);

  /** 5) Free the MMG3D5 structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  free(filename);
  filename = NULL;

  free(outname);
  outname = NULL;

  return(0);
}
