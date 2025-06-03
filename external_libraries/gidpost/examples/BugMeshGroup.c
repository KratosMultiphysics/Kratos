#include <stdlib.h>
#include "gidpost.h"


//#define ASCII_FORMAT
//#define USE_MESH_GROUP


typedef struct {
  int id;
  float x, y, z;
} SCoord;

typedef struct {
  int id;
  int n1, n2, n3;
} SElem;

SCoord G_nodes[9];
SElem G_elems[8];

void GenMesh()
{
  int i, j, idx;
  SCoord * node;
  SElem * elm1, * elm2;

  idx = 1;
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++, idx++ ) {
      node = G_nodes + idx -1;
      node->id = idx;
      node->x = (float)(i);
      node->y = (float)(j);
      node->z = 0;
    }
  }
  idx = 0;
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++, idx+=2 ) {
      elm1 = G_elems + idx;
      elm2 = elm1 + 1;
      elm1->id = idx+1;
      elm1->n1 = i*3 + j + 1;
      elm1->n3 = i*3 + j + 1 + 3 ;
      elm1->n2 = elm1->n3 + 1;
      elm2->id = idx+2;
      elm2->n1 = elm1->n1;
      elm2->n2 = elm1->n1 + 1;
      elm2->n3 = elm1->n2;
    }
  }
}

float Random()
{
  return rand()/(float)(RAND_MAX);
}

int main()
{
  int i;
  SCoord * node;
  SElem * elm;
  int elmi[4];
  /* int elmi2[10];*/
  GiD_FILE fdm, fdr;

  GenMesh();

  GiD_PostInit();

#if defined(ASCII_FORMAT)
  fdm = GiD_fOpenPostMeshFile( "test_fd.post.msh", GiD_PostAscii);
  fdr = GiD_fOpenPostResultFile( "test_fd.post.res", GiD_PostAscii);
#else
  fdm = fdr = GiD_fOpenPostResultFile( "test_fd.post.bin", GiD_PostBinary);
#endif







  /* write mesh info */

  // TODO START MESHGROUP
#ifdef USE_MESH_GROUP
  GiD_fBeginMeshGroup(fdm, "Mesh 1");
#endif


  GiD_fBeginMeshColor(fdm, "TestMsh", GiD_2D, GiD_Triangle, 3, 0, 0.99, 0);
  /* coordinates */
  GiD_fBeginCoordinates(fdm);
  for ( i = 0; i < 9; i++ ) {
    node = G_nodes + i;
    GiD_fWriteCoordinates(fdm, node->id, node->x, node->y, node->z );
  }
  GiD_fEndCoordinates(fdm);
  /* elements */
  GiD_fBeginElements(fdm);
  for ( i = 0; i < 8; i++ ) {
    elm = G_elems + i;
    elmi[0] = elm->n1;
    elmi[1] = elm->n2;
    elmi[2] = elm->n3;
    elmi[3] = 2;
    GiD_fWriteElementMat(fdm, elm->id, elmi);
  }
  GiD_fEndElements(fdm);
  GiD_fEndMesh(fdm);


  // TODO END MESHGROUP
#ifdef USE_MESH_GROUP
  GiD_fEndMeshGroup(fdm);
#endif






  /* now results info */


  // TODO Start ONMESHGROUP
#ifdef USE_MESH_GROUP
  GiD_fBeginOnMeshGroup(fdr, "Mesh 1");
#endif

  GiD_fBeginResult(fdr, "EscalarNodos", "Analysis", 1.0, GiD_Scalar, GiD_OnNodes, NULL, NULL, 0, NULL);
  for ( i = 0; i < 9; i++ ) {
    GiD_fWriteScalar(fdr, G_nodes[i].id, Random());
  }
  GiD_fEndResult(fdr);

  // TODO END MESHGROUP
#ifdef USE_MESH_GROUP
  GiD_fEndOnMeshGroup(fdr);
#endif




#ifdef ASCII_FORMAT
  GiD_fClosePostMeshFile(fdm);
#endif
  GiD_fClosePostResultFile(fdr);

  GiD_PostDone();

  return 0;
}
