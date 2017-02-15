/**
 * Example of use of the mmg2d library (basic use of mesh adaptation)
 *
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
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

/** Include the mmg2d library hader file */
// if the header file is in the "include" directory
// #include "libmmg2d.h"
// if the header file is in "include/mmg/mmg2d"
#include "mmg/mmg2d/libmmg2d.h"

#define MAX0(a,b)     (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX0(a,b)) > (MAX0(c,d))) ? (MAX0(a,b)) : (MAX0(c,d)))

int main(int argc,char *argv[]) {
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  int             ier,k;
  /* To save final mesh in a file */
  FILE*           inm;
  /* To manually recover the mesh */
  int             np, nt, na, nc, nr, nreq, typEntity, typSol;
  int             ref, Tria[3], Edge[2], *corner, *required, *ridge;
  double          Point[3],Sol;


  fprintf(stdout,"  -- TEST MMG2DLIB \n");

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

  /** Manually set of the mesh */
  /** a) give the size of the mesh: 4 vertices, 2 triangles, 4 edges */
  /* allocation */
  if ( MMG2D_Set_meshSize(mmgMesh,4,2,4) != 1 )  exit(EXIT_FAILURE);

/** b) give the vertices: for each vertex, give the coordinates, the reference
      and the position in mesh of the vertex */
  if ( MMG2D_Set_vertex(mmgMesh,0  ,0  ,0  ,  1) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_vertex(mmgMesh,1  ,0  ,0  ,  2) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_vertex(mmgMesh,1  ,1  ,0  ,  3) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_vertex(mmgMesh,0  ,1  ,0  ,  4) != 1 )  exit(EXIT_FAILURE);

 /** c) give the triangles: for each triangle,
      give the vertices index, the reference and the position of the triangle */
  if ( MMG2D_Set_triangle(mmgMesh,  1,  2,  4, 1, 1) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_triangle(mmgMesh,  2,  3,  4, 1, 2) != 1 )  exit(EXIT_FAILURE);


  /** d) give the edges (not mandatory): for each edge,
      give the vertices index, the reference and the position of the edge */
  if ( MMG2D_Set_edge(mmgMesh,  1,  2, 1, 1) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_edge(mmgMesh,  2,  3, 2, 2) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_edge(mmgMesh,  3,  4, 3, 3) != 1 )  exit(EXIT_FAILURE);
  if ( MMG2D_Set_edge(mmgMesh,  4,  1, 4, 4) != 1 )  exit(EXIT_FAILURE);


  /** 3) Build sol in MMG5 format */
  /** Two solutions: just use the MMG2D_loadSol function that will read a .sol(b)
      file formatted or manually set your sol using the MMG2D_Set* functions */

  /** Manually set of the sol */
  /** a) give info for the sol structure: sol applied on vertex entities,
      number of vertices=4, the sol is scalar*/
  if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,4,MMG5_Scalar) != 1 )
    exit(EXIT_FAILURE);

  /** b) give solutions values and positions */
  for(k=1 ; k<=4 ; k++) {
    if ( MMG2D_Set_scalarSol(mmgSol,0.01,k) != 1 ) exit(EXIT_FAILURE);
  }

 /** 4) (not mandatory): check if the number of given entities match with mesh size */
  if ( MMG2D_Chk_meshData(mmgMesh,mmgSol) != 1 ) exit(EXIT_FAILURE);

  /*save init mesh*/
  if ( MMG2D_saveMesh(mmgMesh,"init.mesh") != 1 )
    exit(EXIT_FAILURE);
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"init.sol") != 1 )
    exit(EXIT_FAILURE);

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

  /** 1) Manually get the mesh (in this example we show how to save the mesh
      in the result.mesh file) */

  if( !(inm = fopen("result.mesh","w")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN result.mesh FILE.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(inm,"MeshVersionFormatted 2\n");
  fprintf(inm,"\nDimension 2\n");
  /** a) get the size of the mesh: vertices, tetra, triangles, edges */
  if ( MMG2D_Get_meshSize(mmgMesh,&np,&nt,&na) !=1 )  exit(EXIT_FAILURE);

  /* Table to know if a vertex is corner */
  corner = (int*)calloc(np+1,sizeof(int));
  if ( !corner ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  /* Table to know if a vertex/tetra/tria/edge is required */
  required = (int*)calloc(MAX4(np,0,nt,na)+1 ,sizeof(int));
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
    if ( MMG2D_Get_vertex(mmgMesh,&(Point[0]),&(Point[1]),
                          &ref,&(corner[k]),&(required[k])) != 1 )
      exit(EXIT_FAILURE);
    fprintf(inm,"%.15lg %.15lg %d \n",Point[0],Point[1],ref);
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
    if ( MMG2D_Get_triangle(mmgMesh,&(Tria[0]),&(Tria[1]),&(Tria[2]),
                            &ref,&(required[k])) != 1 )
      exit(EXIT_FAILURE);
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
    if ( MMG2D_Get_edge(mmgMesh,&(Edge[0]),&(Edge[1]),&ref,
                        &(ridge[k]),&(required[k])) != 1 )  exit(EXIT_FAILURE);
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
  if ( MMG2D_Get_solSize(mmgMesh,mmgSol,&typEntity,&np,&typSol) != 1 )
    exit(EXIT_FAILURE);

  if ( ( typEntity != MMG5_Vertex )  || ( typSol != MMG5_Scalar ) )
    exit(EXIT_FAILURE);

  fprintf(inm,"\nSolAtVertices\n%d\n",np);
  fprintf(inm,"1 1 \n\n");
  for(k=1; k<=np; k++) {
    /** b) Vertex recovering */
    if ( MMG2D_Get_scalarSol(mmgSol,&Sol) != 1 )  exit(EXIT_FAILURE);
    fprintf(inm,"%.15lg \n",Sol);
  }
  fprintf(inm,"\nEnd\n");
  fclose(inm);


  /*save result*/
  if ( MMG2D_saveMesh(mmgMesh,"result.mesh") != 1 )
    exit(EXIT_FAILURE);

  /*save metric*/
  if ( MMG2D_saveSol(mmgMesh,mmgSol,"result") != 1 )
    exit(EXIT_FAILURE);

  /** 3) Free the MMG2D structures */
  MMG2D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return(0);
}
