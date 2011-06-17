#ifndef GID_MESH_LIBRARY_H
#define GID_MESH_LIBRARY_H

#include <stdlib.h>

typedef void *GiDInput_Handle;
typedef void *GiDOutput_Handle;

//TYPES OF ELEMENTS
//  GiD_NoElement=0
//  GiD_Point=1
//  GiD_Linear=2
//  GiD_Triangle=3
//  GiD_Quadrilateral=4
//  GiD_Tetrahedra=5
//  GiD_Hexahedra=6
//  GiD_Prism=7
//  GiD_Pyramid=8
//  GiD_Sphere=9
//  GiD_Circle=10


//FUNCTIONS OF THE LIBRARY

//this function initializies an input Handle with a first type of elements and its nodes.
GiDInput_Handle GetGidInputHandleWithBoundaryMesh(int ndime,int nlen,double* coords,int flen,int* faces,int type_faces_mesh);
GiDInput_Handle GetGidInputHandleWithBoundaryMesh(int ndime,int nlen,double* coords,int flen,int* faces,int number_of_types,int* flen_by_types,int* type_faces_mesh);
GiDInput_Handle GetGidInputHandleWithVolumeMesh(int ndime,int nlen,double* coords,int elen,int* elems,int type_elem_mesh);
GiDInput_Handle GetGidInputHandleWithVolumeMesh(int ndime,int nlen,double* coords,int elen,int* elems,int number_of_types,int* elen_by_types,int* type_elem_mesh);
void DeleteInputHandle(GiDInput_Handle hdl_gin);
//Add nodes to in the entry structure because user wants that the mesher use
// in the generation.
void AddNodesInMesh(GiDInput_Handle hdl_gin,int nlen,double* coord);
//Add attributes to the nodes of the entry
void AddAttributeMeshNodes(GiDInput_Handle hdl_gin, double* attributes);
//Add markers to the nodes of the entry
void AddMarkersMeshNodes(GiDInput_Handle hdl_gin,double* markers);

void* GetGidOutputHandle();

void DeleteOutputHandle(GiDOutput_Handle hdl_gout);
double* GetCoordNodes(GiDOutput_Handle hdl_gout);
int GetNumberOfNodes(GiDOutput_Handle hdl_gout);
double* GetAttributesOfNodes(GiDOutput_Handle hdl_gout);
int GetNumberOfAttributesByNode(GiDOutput_Handle hdl_gout);
int* GetMarkersOfNodes(GiDOutput_Handle hdl_gout);
int GetNumberOfMarkersByNode(GiDOutput_Handle hdl_gout);
int* GetConnecElements(GiDOutput_Handle hdl_gout);
int GetNumberOfElements(GiDOutput_Handle hdl_gout);
double* GetCoordAddedNodes(GiDOutput_Handle hdl_gout);
int GetNumberOfAddedNodes(GiDOutput_Handle hdl_gout);
int* GetBaseNodesForInterpolation(GiDOutput_Handle hdl_gout);
double* GetWightsNodesForInterpolation(GiDOutput_Handle hdl_gout);
int GetNumberOfNodesForInterpolation(GiDOutput_Handle hdl_gout);
int* GetMarkersOfAddedNodes(GiDOutput_Handle hdl_gout);
int GeNumberOftMarkersByAddedNode(GiDOutput_Handle hdl_gout);
int* GetNeighborsByFaces(GiDOutput_Handle hdl_gout);

//functions that returns several results: mesh nodes,mesh elems,attributes nodes,markers nodes, added nodes in the final mesh
//....


//************************************************************************
//mesh generation functions
//************************************************************************

//Mesh generation using Delaunay method
//int GiDML_DelaunayConstrainedTetrahedraMesher(GiDInput_Handle hdl_gin, GiDOutput_Handle hdl_gout,double size_user,void* interface_gid=0);
//int DelaunayConstrainedTetrahedraMesher_CheckConsistency(GiDInput_Handle hdl_gin);


//Mesh generation using IsoStuffing method
void GiDML_IsoSurfaceTetrahedraMesher(GiDInput_Handle hdl_gin, GiDOutput_Handle hdl_gout, double size_user,void* interface_gid=0);
int IsoSurfaceTetrahedraMesher_CheckConsistency(GiDInput_Handle hdl_gin);
const char *IsoSurfaceTetrahedraMesher_CheckState(int state);


//************************************************************************
//mesh editing functions
//************************************************************************

//Improves local connectivities of bad quality tetras
//int GiDML_ImproveConectivitiesOfTetrahedralMesh(GiDInput_Handle hdl_gin, GiDOutput_Handle hdl_gout);


//************************************************************************
//mesh analysis functions
//************************************************************************

//Gets nodes distance (in a tetrahedra mesh) to a defined triangles skin
void GiDML_GetNodesDistanceRadiusInTetrahedraMesh(GiDInput_Handle hdl_gin, GiDOutput_Handle hdl_gout, double Radius,int cuadratic_distances);
void GiDML_GetNodesDistanceInTetrahedraMesh(GiDInput_Handle hdl_gin, GiDOutput_Handle hdl_gout,int cuadratic_distances);


//Detects if the input nodes are inside or outside a certain skin of triangles.
void GiDML_DetectInsideAndOutsideNodesInTetrahedraMesh(GiDInput_Handle hdl_gin, GiDOutput_Handle hdl_gout);

#endif //GID_MESH_LIBRARY_H
