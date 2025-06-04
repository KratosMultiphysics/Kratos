#pragma once __GIDPOST_CLUSTER_FUNCTIONS__

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

GIDPOST_API
int GiD_WriteCluster(int id, int nid);

GIDPOST_API
int GiD_fWriteCluster(GiD_FILE fd, int id, int nid);

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

GIDPOST_API
int GiD_WriteClusterMat(int id, int nid, int mat);

GIDPOST_API
int GiD_fWriteClusterMat(GiD_FILE fd, int id, int nid, int mat);
