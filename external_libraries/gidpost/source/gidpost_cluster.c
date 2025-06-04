/* gidpost 2.11 */
/*
 *  gidpost_cluster.c--
 *
 *    This file implement a C interface for generating postprocessing
 *    results in the 'New postprocess format' of GiD. See declaration
 *    in gidpost.h. This is just a part to give service to DEM Clusters.
 *
 */

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "gidpostInt.h"
#include "gidpost.h"
#include "gidpostHash.h"
#include "gidpostFILES.h"
  
#ifdef HDF5
  #include "gidpostHDF5.h"
#endif


#ifdef WIN32
#define snprintf _snprintf
#endif

/* defined in gidpostFILES.c */
extern CPostFile *G_MeshFile;
extern CPostFile *G_ResultFile;
extern CPostFile *G_outputMesh;
// extern GiD_PostMode G_PostMode;

/*
 *  Write a cluster element member at the current Elements Block.
 *  A cluster element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *  
 */

int _GiD_WriteCluster(CPostFile *File, int id, int nid)
{
  /* state checking */
  assert(_GiDfiles_CheckState(POST_MESH_ELEM, File));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id, 0);
  CPostFile_WriteInteger(File, nid, 1);
  if (CPostFile_IsBinary(File)) {
    CPostFile_WriteInteger(File, 1, 1);    
  }
  return 0;
}

int GiD_WriteCluster(int id, int nid)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteCluster_HDF5(id,nid);
  }
#endif

  return _GiD_WriteCluster(G_outputMesh, id, nid);
}

int GiD_fWriteCluster(GiD_FILE fd, int id, int nid)
{

  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}

  return _GiD_WriteCluster(File, id, nid);
}

/*
 *  Write a cluster element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  A cluster element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     mat: material identification.
 *  
 */

int _GiD_WriteClusterMat(CPostFile * File, int id, int nid, int mat)
{
  /* state checking */
  assert(_GiDfiles_CheckState(POST_MESH_ELEM, File));    
  /* keep on the same state */
  CPostFile_WriteInteger(File, id,  0);
  CPostFile_WriteInteger(File, nid, 1);
  CPostFile_WriteInteger(File, mat, 2);
  return 0;
}

int GiD_WriteClusterMat(int id, int nid, int mat)
{
#ifdef HDF5
  if(PostMode==GiD_PostHDF5){
    return GiD_WriteClusterMat_HDF5(id,nid,mat);
  }
#endif

  return _GiD_WriteClusterMat(G_outputMesh, id, nid, mat);
}

int GiD_fWriteClusterMat(GiD_FILE fd, int id, int nid, int mat)
{
  
  CPostFile *File = NULL;
  if ( ( File = FD2FILE( fd)) == NULL) { return -8;}

  return _GiD_WriteClusterMat(File, id, nid, mat);
}
