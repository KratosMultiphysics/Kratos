/**
 * Example of input output for the mmgs library for multiple solutions at mesh
 * vertices
 *
 * \author Algiane Froehly (InriaSoft)
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

/** Include the mmgs library hader file */
// if the header file is in the "include" directory
// #include "libmmgs.h"
// if the header file is in "include/mmg/mmgs"
#include "mmg/mmgs/libmmgs.h"

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol,mmgMet,tmpSol;
  int             i,j,opt;

  /* To manually recover the mesh */
  MMG5_int        np;
  int             nsol,typSol[MMG5_NSOLS_MAX];
  double          *sols;

  /* Filenames */
  char            *filename, *fileout;

  fprintf(stdout,"  -- TEST MMGSLIB \n");

  if ( argc != 4 ) {
    printf(" Usage: %s filein fileout io_option\n",argv[0]);
    printf("     io_option = 0 to Get/Set the solution field by field\n");
    printf("     io_option = 1 to Get/Set the solution field by field"
           " and vertex by vertex\n");
    return(1);
  }

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  fileout = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[2]);

  opt = atoi(argv[3]);


  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh:
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */

  mmgMesh = NULL;
  mmgMet  = NULL;
  mmgSol  = NULL;
  tmpSol  = NULL;
  MMGS_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgMet,
                  MMG5_ARG_end);

  /** 2) Build initial mesh and solutions in MMG5 format */
  /** Two solutions: just use the MMGS_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMGS_Set* functions */

  /** Automatic loading of the mesh and multiple solutions */
  if ( MMGS_loadMesh(mmgMesh,filename) != 1 )  exit(EXIT_FAILURE);

  if ( MMGS_loadAllSols(mmgMesh,&mmgSol,filename) != 1 )
    exit(EXIT_FAILURE);

  /** ------------------------------ STEP II --------------------------- */

  /** 3) Transfer the solutions in a new solutions array */
  /** a) Get the solutions sizes */
  if ( MMGS_Get_solsAtVerticesSize(mmgMesh,&mmgSol,&nsol,&np,typSol) != 1 )
    exit(EXIT_FAILURE);

  /** b) Manually set the size of the new solution: give info for the sol
      structure: number of solutions, type of entities on which applied the
      solutions, number of vertices, type of the solution */
  if ( MMGS_Set_solsAtVerticesSize(mmgMesh,&tmpSol,nsol,np,typSol) != 1 )
    exit(EXIT_FAILURE);

  /** c) Get each solution and set it in the new structure */

  /** b) give solutions values and positions */
  /* Get the entire field of a given solution */
  for ( i=1; i<=nsol; ++i ) {
    if ( !opt ) {
      /* Get the ith solution array */
      if ( typSol[i-1] == MMG5_Scalar )
        sols = (double*) calloc(np, sizeof(double));
      else if ( typSol[i-1] == MMG5_Vector )
        sols = (double*) calloc(np*3, sizeof(double));
      else if ( typSol[i-1] == MMG5_Tensor ) {
        sols = (double*) calloc(np*6, sizeof(double));
      }
      else {
        printf("Unexpected value for solution type\n");
        exit(EXIT_FAILURE);
      }

      if ( MMGS_Get_ithSols_inSolsAtVertices(mmgSol,i,sols) !=1 ) exit(EXIT_FAILURE);

      /* Set the ith solution in the new structure */
      if ( MMGS_Set_ithSols_inSolsAtVertices(tmpSol,i,sols) !=1 ) exit(EXIT_FAILURE);
    }
    else {
      /* Get the ith solution array vertex by vertex */
      if ( typSol[i-1] == MMG5_Scalar )
        sols = (double*) calloc(1, sizeof(double));
      else if ( typSol[i-1] == MMG5_Vector )
        sols = (double*) calloc(3, sizeof(double));
      else if ( typSol[i-1] == MMG5_Tensor ) {
        sols = (double*) calloc(6, sizeof(double));
      }
      else {
        printf("Unexpected value for solution type\n");
        exit(EXIT_FAILURE);
      }

      for ( j=1; j<=np; ++j ) {
        if ( MMGS_Get_ithSol_inSolsAtVertices(mmgSol,i,sols,j) !=1 ) exit(EXIT_FAILURE);

        /* Set the ith solution vertex by vertex in the new structure */
        if ( MMGS_Set_ithSol_inSolsAtVertices(tmpSol,i,sols,j) !=1 ) exit(EXIT_FAILURE);
      }

    }
    free(sols); sols = NULL;
  }


  /** ------------------------------ STEP III -------------------------- */
  /** Save the new data */
  /** Use the MMGS_saveMesh/MMGS_saveAllSols functions */
  /* save the mesh */
  if ( MMGS_saveMesh(mmgMesh,fileout) != 1 )
    exit(EXIT_FAILURE);

  /*s ave the solutions array */
  if ( MMGS_saveAllSols(mmgMesh,&tmpSol,fileout) != 1 )
    exit(EXIT_FAILURE);

  /** 3) Free the MMGS structures */
  MMGS_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppSols,&tmpSol,
                 MMG5_ARG_ppSols,&mmgSol,
                 MMG5_ARG_end);

  free(filename);
  filename = NULL;

  free(fileout);
  fileout = NULL;

  return(0);
}
