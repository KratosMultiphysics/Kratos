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
#include "mmg/mmg3d/libmmg3d.h"

#define MAX0(a,b)     (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX0(a,b)) > (MAX0(c,d))) ? (MAX0(a,b)) : (MAX0(c,d)))

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  int             ier,k;
  /* To save final mesh in a file */
  FILE*           inm;
  /* To manually recover the mesh */
  int             np, ne, nt, na, nc, nr, nreq, typEntity, typSol;
  int             ref, Tetra[4], Tria[3], Edge[2], *corner, *required, *ridge;
  double          Point[3],Sol;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of mesh and sol structures */
  /* args of InitMesh: mesh=&mmgMesh, sol=&mmgSol */
  mmgMesh = NULL;
  mmgSol  = NULL;
  MMG5_Init_mesh(&mmgMesh,&mmgSol,NULL);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the MMG5_loadMesh function that will read a .mesh(b)
      file formatted or manually set your mesh using the MMG5_Set* functions */

  /** Manually set of the mesh */
  /** a) give the size of the mesh: 12 vertices, 12 tetra, 20 triangles, 0 edges */
  if ( !MMG5_Set_meshSize(mmgMesh,12,12,20,0) )  exit(EXIT_FAILURE);

  /** b) give the vertices: for each vertex, give the coordinates, the reference
      and the position in mesh of the vertex */
  if ( !MMG5_Set_vertex(mmgMesh,0  ,0  ,0  ,0,  1) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,0  ,0  ,0,  2) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,0  ,1  ,0,  3) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0  ,0  ,1  ,0,  4) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0  ,1  ,0  ,0,  5) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,1  ,0  ,0,  6) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0.5,1  ,1  ,0,  7) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,0  ,1  ,1  ,0,  8) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,0  ,0  ,0,  9) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,1  ,0  ,0, 10) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,0  ,1  ,0, 11) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_vertex(mmgMesh,1  ,1  ,1  ,0, 12) )  exit(EXIT_FAILURE);

  /** c) give the tetrahedras: for each tetrahedra,
      give the vertices index, the reference and the position of the tetra */
  if ( !MMG5_Set_tetrahedron(mmgMesh,  1,  4,  2,  8,1, 1) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  8,  3,  2,  7,1, 2) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  5,  2,  6,  8,1, 3) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  5,  8,  1,  2,1, 4) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  7,  2,  8,  6,1, 5) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  2,  4,  3,  8,1, 6) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  9,  2,  3,  7,2, 7) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  7, 11,  9, 12,2, 8) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  6,  9, 10,  7,2, 9) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  6,  7,  2,  9,2,10) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh, 12,  9,  7, 10,2,11) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_tetrahedron(mmgMesh,  9,  3, 11,  7,2,12) )  exit(EXIT_FAILURE);

  /** d) give the triangles (not mandatory): for each triangle,
      give the vertices index, the reference and the position of the triangle */
  if ( !MMG5_Set_triangle(mmgMesh,  1,  4,  8, 3, 1) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  1,  2,  4, 3, 2) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  8,  3,  7, 3, 3) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  8,  6, 3, 4) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  6,  2, 3, 5) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  2,  1, 3, 6) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  5,  1,  8, 3, 7) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  7,  6,  8, 3, 8) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  4,  3,  8, 3, 9) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  2,  3,  4, 3,10) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  9,  3,  2, 4,11) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh, 11,  9, 12, 4,12) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  7, 11, 12, 4,13) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  6,  7, 10, 4,14) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  6, 10,  9, 4,15) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  6,  9,  2, 4,16) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh, 12, 10,  7, 4,17) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh, 12,  9, 10, 4,18) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  3, 11,  7, 4,19) )  exit(EXIT_FAILURE);
  if ( !MMG5_Set_triangle(mmgMesh,  9, 11,  3, 4,20) )  exit(EXIT_FAILURE);


  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG5_loadMet function that will read a .sol(b)
      file formatted or manually set your sol using the MMG5_Set* functions */

  /** Manually set of the sol */
  /** a) give info for the sol structure: sol applied on vertex entities,
      number of vertices=12, the sol is scalar*/
  if ( !MMG5_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,12,MMG5_Scalar) )
    exit(EXIT_FAILURE);

  /** b) give solutions values and positions */
  for(k=1 ; k<=12 ; k++) {
    if ( !MMG5_Set_scalarSol(mmgSol,0.5,k) ) exit(EXIT_FAILURE);
  }

  /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( !MMG5_Chk_meshData(mmgMesh,mmgSol) ) exit(EXIT_FAILURE);

  /** ------------------------------ STEP  II -------------------------- */
  /** remesh function */
  ier = MMG5_mmg3dlib(mmgMesh,mmgSol);
  if ( ier == MMG5_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH\n");
    return(ier);
  } else if ( ier == MMG5_LOWFAILURE )
    fprintf(stdout,"BAD ENDING OF MMG3DLIB\n");


  /** ------------------------------ STEP III -------------------------- */
  /** get results */
  /** Two solutions: just use the MMG5_saveMesh/MMG5_saveMet functions
      that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
      using the MMG5_getMesh/MMG5_getSol functions */

  /** 1) Manually get the mesh (in this example we show how to save the mesh
      in the mesh.o.mesh file) */
  if( !(inm = fopen("mesh.o.mesh","w")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN mesh.o.mesh FILE.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(inm,"MeshVersionFormatted 2\n");
  fprintf(inm,"\nDimension 3\n");

  /** a) get the size of the mesh: vertices, tetra, triangles, edges */
  if ( !MMG5_Get_meshSize(mmgMesh,&np,&ne,&nt,&na) )  exit(EXIT_FAILURE);

  /* Table to know if a vertex is corner */
  corner = (int*)calloc(np+1,sizeof(int));
  if ( !corner ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  /* Table to know if a vertex/tetra/tria/edge is required */
  required = (int*)calloc(MAX4(np,ne,nt,na)+1 ,sizeof(int));
  if ( !required ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  /* Table to know if a coponant is corner and/or required */
  ridge = (int*)calloc(na+1 ,sizeof(int));
  if ( !ridge ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }

  nreq = 0; nc = 0;
  fprintf(inm,"\nVertices\n%d\n",np);
  for(k=1; k<=np; k++) {
    /** b) Vertex recovering */
    if ( !MMG5_Get_vertex(mmgMesh,&(Point[0]),&(Point[1]),&(Point[2]),
			  &ref,&(corner[k]),&(required[k])) )  exit(EXIT_FAILURE);
    fprintf(inm,"%.15lg %.15lg %.15lg %d \n",Point[0],Point[1],Point[2],ref);
    if ( corner[k] )  nc++;
    if ( required[k] )  nreq++;
  }
  fprintf(inm,"\nCorners\n%d\n",nc);
  for(k=1; k<=np; k++) {
    if ( corner[k] )  fprintf(inm,"%d \n",k);
  }
  fprintf(inm,"\nRequiredVertices\n%d\n",nreq);
  for(k=1; k<=np; k++) {
    if ( required[k] )  fprintf(inm,"%d \n",k);
  }
  free(corner);
  corner = NULL;

  nreq = 0;
  fprintf(inm,"\nTriangles\n%d\n",nt);
  for(k=1; k<=nt; k++) {
    /** d) Triangles recovering */
    if ( !MMG5_Get_triangle(mmgMesh,&(Tria[0]),&(Tria[1]),&(Tria[2]),
			    &ref,&(required[k])) )  exit(EXIT_FAILURE);
    fprintf(inm,"%d %d %d %d \n",Tria[0],Tria[1],Tria[2],ref);
    if ( required[k] )  nreq++;
  }
  fprintf(inm,"\nRequiredTriangles\n%d\n",nreq);
  for(k=1; k<=nt; k++) {
    if ( required[k] )  fprintf(inm,"%d \n",k);
  }

  nreq = 0;nr = 0;
  fprintf(inm,"\nEdges\n%d\n",na);
  for(k=1; k<=na; k++) {
    /** e) Edges recovering */
    if ( !MMG5_Get_edge(mmgMesh,&(Edge[0]),&(Edge[1]),&ref,
			&(ridge[k]),&(required[k])) )  exit(EXIT_FAILURE);
    fprintf(inm,"%d %d %d \n",Edge[0],Edge[1],ref);
    if ( ridge[k] )  nr++;
    if ( required[k] )  nreq++;
  }
  fprintf(inm,"\nRequiredEdges\n%d\n",nreq);
  for(k=1; k<=na; k++) {
    if ( required[k] )  fprintf(inm,"%d \n",k);
  }
  fprintf(inm,"\nRidges\n%d\n",nr);
  for(k=1; k<=na; k++) {
    if ( ridge[k] )  fprintf(inm,"%d \n",k);
  }

  nreq = 0;
  fprintf(inm,"\nTetrahedra\n%d\n",ne);
  for(k=1; k<=ne; k++) {
    /** c) Tetra recovering */
    if ( !MMG5_Get_tetrahedron(mmgMesh,&(Tetra[0]),&(Tetra[1]),&(Tetra[2]),&(Tetra[3]),
			       &ref,&(required[k])) )  exit(EXIT_FAILURE);
    fprintf(inm,"%d %d %d %d %d \n",Tetra[0],Tetra[1],Tetra[2],Tetra[3],ref);
    if ( required[k] )  nreq++;
  }
  fprintf(inm,"\nRequiredTetrahedra\n%d\n",nreq);
  for(k=1; k<=ne; k++) {
    if ( required[k] )  fprintf(inm,"%d \n",k);
  }

  fprintf(inm,"\nEnd\n");
  fclose(inm);

  free(required);
  required = NULL;
  free(ridge);
  ridge    = NULL;

  /** 2) Manually get the solution (in this example we show how to save the
      solution in the mesh.o.sol file) */
  if( !(inm = fopen("mesh.o.sol","w")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN mesh.o.sol FILE.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(inm,"MeshVersionFormatted 2\n");
  fprintf(inm,"\nDimension 3\n");

  /** a) get the size of the sol: type of entity (SolAtVertices,...),
      number of sol, type of solution (scalar, tensor...) */
  if ( !MMG5_Get_solSize(mmgMesh,mmgSol,&typEntity,&np,&typSol) )
    exit(EXIT_FAILURE);

  if ( ( typEntity != MMG5_Vertex )  || ( typSol != MMG5_Scalar ) )
    exit(EXIT_FAILURE);

  fprintf(inm,"\nSolAtVertices\n%d\n",np);
  fprintf(inm,"1 1 \n\n");
  for(k=1; k<=np; k++) {
    /** b) Vertex recovering */
    if ( !MMG5_Get_scalarSol(mmgSol,&Sol) )  exit(EXIT_FAILURE);
    fprintf(inm,"%.15lg \n",Sol);
  }
  fprintf(inm,"\nEnd\n");
  fclose(inm);

  /** 3) Free the MMG3D5 structures */
  MMG5_Free_all(mmgMesh,mmgSol,NULL);

  return(ier);
}
