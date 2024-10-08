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
 * Example of use of the mmg3d library to load an save files whose format is
 * detected using the extension), get the mesh from Mmg to store it into user
 * data structures and get the adjacency relationship between tetra.
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
  int             ier,silent=0;
  /* To save final mesh in a file */
  FILE*           inm;
  char            *fileout_generic,*fileout_medit,*filein;


  fprintf(stdout,"  -- TEST LOAD AND GET MESH DATA \n");

  if ( argc != 4 ) {
    printf(" Usage: %s filein fileout silent_mode\n",argv[0]);
    return(EXIT_FAILURE);
  }

  /* Name and path of the input mesh file */
  filein = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filein == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(filein,argv[1]);

  fileout_generic = (char *) calloc(strlen(argv[2])+1, sizeof(char));
  if ( fileout_generic == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(fileout_generic,argv[2]);

  fileout_medit = (char *) calloc(strlen(argv[2]) + 6, sizeof(char));
  if ( fileout_medit == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(fileout_medit,argv[2]);
  strcat(fileout_medit,".mesh");

  silent = atoi(argv[3]);

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh  structure of Mmg */
  /* args of InitMesh:
   * MMG5_ARG_start: we start to give the args of a variadic func
   * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
   * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
   * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
   * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */
  mmgMesh = NULL;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,
                  MMG5_ARG_end);


  /** 2) Mesh loading depending on detected extension */
  ier = MMG3D_loadGenericMesh(mmgMesh,NULL,filein);

  if ( ier<1 ) {
    if ( ier==0 ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",filein);
      fprintf(stderr,"  ** UNABLE TO OPEN INPUT FILE.\n");
      return EXIT_FAILURE;
    }
    else {
      assert ( ier == -1 );
      fprintf(stderr,"  ** UNABLE TO READ INPUT FILE.\n");
      return EXIT_FAILURE;
    }
  }

  /** 3) Get the mesh from Mmg data structure into local data structure */

  /** a) get the size of the mesh: np = #vertices, ne = #tetra, npr = #prisms,
   * nt = #triangles, nq = #quadrangls, na = #edges */
  MMG5_int  np, ne, nt, na;
  if ( MMG3D_Get_meshSize(mmgMesh,&np,&ne,NULL,&nt,NULL,&na) !=1 ) {
    fprintf(stderr,"  ERROR: Unable to get mesh size.\n");
    return EXIT_FAILURE;
  }

  /** b) get point coordinates and references: i^th point coordinates are stored
   * in points[3*(i-1)]@3. Reference of this point is in ref[i-1] */
  double   *points=NULL; // point coordinates
  MMG5_int *ref=NULL; // point references (==gmsh tags, ==colors)

  points = (double *)malloc(3*np*sizeof(double));
  assert(points);
  ref    = (MMG5_int *)malloc(np*sizeof(MMG5_int));
  assert(ref);

  ier = MMG3D_Get_vertices(mmgMesh,points, ref,NULL,NULL);
  if ( !ier ) {
    fprintf(stderr,"  ERROR: Unable to get mesh vertices.\n");
    return EXIT_FAILURE;
  }

  /** c) get tetra and tetra references:  i^th tetra vertices are stored
   * in tetra[4*(i-1)]@4. Reference of this tetra is in tetref[i-1]*/
  MMG5_int    *tetref=NULL, *tetra=NULL;

  tetra = (MMG5_int *)malloc(4*ne*sizeof(MMG5_int));
  assert(tetra);
  tetref = (MMG5_int *)malloc(ne*sizeof(MMG5_int));
  assert(tetref);

  ier = MMG3D_Get_tetrahedra(mmgMesh,tetra, tetref,NULL);
  if ( !ier ) {
    fprintf(stderr,"  ERROR: Unable to get mesh tetra.\n");
    return EXIT_FAILURE;
  }

  /** ... etc... */

  /** 4) Mesh saving at Medit format to check it */
  MMG5_int k;

  if( !(inm = fopen(fileout_medit,"w")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT MESH FILE.\n");
    exit(EXIT_FAILURE);
  }
  /* Header */
  fprintf(inm,"MeshVersionFormatted 2\n");
  fprintf(inm,"\nDimension 3\n");

  /* Mesh vertices */
  fprintf(inm,"\nVertices\n%"MMG5_PRId"\n",np);
  for(k=1; k<=np; k++) {
    MMG5_int address = 3*(k-1);
    fprintf(inm,"%.15lg %.15lg %.15lg %"MMG5_PRId" \n",
            points[address],points[address+1],points[address+2],ref[k-1]);
  }

  /* Mesh tetra */
  fprintf(inm,"\nTetrahedra\n%"MMG5_PRId"\n",ne);
  for(k=1; k<=ne; k++) {
    MMG5_int address = 4*(k-1);
    fprintf(inm,"%"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId" %"
            MMG5_PRId" \n",tetra[address],tetra[address+1],
            tetra[address+2],tetra[address+3],tetref[k-1]);
  }

  fprintf(inm,"\nEnd\n");
  fclose(inm);

  /* Memory release */
  free(points); free(ref);
  free(tetra);  free(tetref);
  free(fileout_medit);
  fileout_medit = NULL;

  /** 5) get and print tetra adjacencies in terminal */
  for ( k=1; k<=ne; ++k ) {
    MMG5_int adja[4];
    ier = MMG3D_Get_adjaTet(mmgMesh,k,adja);
    if ( !ier ) {
      fprintf(stderr,"  ERROR: Unable to get adjacents of tetra %"MMG5_PRId".\n",k);
      return EXIT_FAILURE;
    }

    if ( !silent ) {
      fprintf(stdout,"Tetra %"MMG5_PRId" is adjacent to tetras %"MMG5_PRId" %"
              MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId"\n",
              k,adja[0],adja[1],adja[2],adja[3]);
    }
  }

  /** 6) Save output file depending on the detected extension */
  ier = MMG3D_saveGenericMesh(mmgMesh,NULL,fileout_generic);

  if ( ier<1 ) {
    if ( ier==0 ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",fileout_generic);
      fprintf(stderr,"  ** UNABLE TO OPEN INPUT FILE.\n");
      return EXIT_FAILURE;
    }
  }
  free(fileout_generic);
  fileout_generic = NULL;

  /** 7) Free the MMG3D5 structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,
                 MMG5_ARG_end);

  return(EXIT_SUCCESS);
}
