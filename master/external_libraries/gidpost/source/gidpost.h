/* gidpost 2.1 */
/*
 *  gidpost.h --
 *
 *    This file declare a C interface for generating postprocessing
 *    results in the 'New postprocess format' of GiD. The interface is
 *    composed by a set of functions according to the specification of
 *    the GiD postprocess file format from version 7.2. Althought the
 *    specification only describe the ascii format, with this
 *    interface is possible to write in ascii, ascii zipped and binary
 *    format. Binary files are always zipped. All functions returns an
 *    integer indicating fail with a non zero (!0) value or success
 *    with a zero value (0).
 *
 */

#ifndef __GIDPOST__
#define __GIDPOST__

/* if this is undefined, it is not necessary to link with library HDF5 and format GiD_PostHDF5 gets deactivated
   this macro should be defined in the project configuration
*/
/*
#define HDF5
*/

#ifdef HAVE_GIDPOST_CONFIG_H
#include "gidpost_config.h"
/* build GIDPOST_VERSION stringificing GP_VERSION_MAJOR and GP_VERSION_MINOR */
#define __GP_V_STR( str)     __GP_V_STR_( str)
#define __GP_V_STR_( str)     #str
#define GIDPOST_VERSION     "" __GP_V_STR( GP_VERSION_MAJOR) "." __GP_V_STR( GP_VERSION_MINOR) ""
#else
#define GIDPOST_VERSION     "2.1"
#endif

#define GP_CONST const

/*
#if defined(_MSC_VER) && defined(GIDPOST_SHARED)
# ifndef WIN32
#  define WIN32
# endif
# if defined(GIDPOST_EXPORTS)
#  define GIDPOST_API __declspec(dllexport)
# else
#  define GIDPOST_API __declspec(dllimport)
# endif
#else
# define GIDPOST_API
#endif
*/

#if defined (WIN32) && defined (GIDPOST_SHARED)
  #if defined (GIDPOST_EXPORTS)
    #define  GIDPOST_API __declspec(dllexport)
  #else
    #define  GIDPOST_API __declspec(dllimport)
  #endif /* defined (GIDPOST_EXPORTS) */
#else /* defined (WIN32) && defined (GIDPOST_SHARED) */
 #define GIDPOST_API
#endif


#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
  
/* ################################################################################
 * #    command to review contents of the file in GiD_PostHDF5 format:
 * #           h5dump --string FILE.flavia.res > FILE.txt
 * #
 * ################################################################################
*/

/* opening mode */

typedef enum { 
  GiD_PostAscii=0, 
  GiD_PostAsciiZipped, 
  GiD_PostBinary,
  GiD_PostHDF5
} GiD_PostMode;

#define GiD_PostBinary2 GiD_PostHDF5 /* back compatibility */

/* domain dimension */

typedef enum { GiD_2D = 2, GiD_3D = 3 } GiD_Dimension;

/* post element types */

typedef enum {
  GiD_NoElement = 0,
  GiD_Point,
  GiD_Linear,
  GiD_Triangle,
  GiD_Quadrilateral,
  GiD_Tetrahedra,
  GiD_Hexahedra,
  GiD_Prism,
  GiD_Pyramid,
  GiD_Sphere,
  GiD_Circle,
  GiD_Cluster
} GiD_ElementType;

typedef enum {
  GiD_Scalar = 0,
  GiD_Vector,
  GiD_Matrix,
  GiD_PlainDeformationMatrix,
  GiD_MainMatrix,
  GiD_LocalAxes,
  GiD_ComplexScalar,
  GiD_ComplexVector
} GiD_ResultType;

typedef enum { GiD_OnNodes=0, GiD_OnGaussPoints } GiD_ResultLocation;

typedef unsigned int GiD_FILE;

/*
  GiD_PostInit -- Initialization of gidpost library, must be called
  from the main thread of the program.
 */
GIDPOST_API int GiD_PostInit();
GIDPOST_API int GiD_PostDone();
  
/* ---------------------------------------------------------------------------
 *
 *  Post Mesh Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post mesh file.
 */

GIDPOST_API
int GiD_OpenPostMeshFile(GP_CONST char * FileName,
		          GiD_PostMode Mode );

GIDPOST_API
GiD_FILE GiD_fOpenPostMeshFile(GP_CONST char * FileName,
		               GiD_PostMode Mode );

/*
 *  Close the current post mesh file
 */

GIDPOST_API
int GiD_ClosePostMeshFile();

GIDPOST_API
int GiD_fClosePostMeshFile(GiD_FILE fd);

/*
 *  Begin a new mesh. After that you should write the nodes and the
 *  elements. When writing the elements it is assumed a connectivity
 *  of size given in the paramenter NNode.
 */

GIDPOST_API
int GiD_BeginMesh(GP_CONST char * MeshName,
		  GiD_Dimension Dim, GiD_ElementType EType,
		  int NNode);
  
GIDPOST_API
int GiD_fBeginMesh(GiD_FILE fd, GP_CONST char * MeshName,
		   GiD_Dimension Dim, GiD_ElementType EType,
		   int NNode );
/*
 *  Begin a new mesh. After that you should write the nodes and the
 *  elements. When writing the elements it is assumed a connectivity
 *  of size given in the paramenter NNode. With this function you can
 *  specify a color for the mesh by its RGB components where each
 *  component take values on the interval [0,1].
 */

GIDPOST_API
int GiD_BeginMeshColor(GP_CONST char * MeshName,
		       GiD_Dimension Dim, GiD_ElementType EType,
		       int NNode,
		       double Red, double Green, double Blue);

GIDPOST_API
int GiD_fBeginMeshColor(GiD_FILE fd, GP_CONST char * MeshName,
		        GiD_Dimension Dim, GiD_ElementType EType,
		        int NNode,
		        double Red, double Green, double Blue);
/*
 *  End current mesh. 
 */

GIDPOST_API
int GiD_EndMesh();

GIDPOST_API
int GiD_fEndMesh(GiD_FILE fd);

/* 
 *    Declares the units of the mesh
 */
GIDPOST_API
int GiD_MeshUnit(GP_CONST char * UnitName);

GIDPOST_API
int GiD_fMeshUnit(GiD_FILE fd,GP_CONST char * UnitName);

  /* 
     Associates a local axes result to the current mesh
   */
GIDPOST_API int GiD_MeshLocalAxes(GP_CONST char * Result, GP_CONST char * Analysis,
		              double step);

/*
 *  Start a coordinate block in the current mesh. All nodes
 *  coordinates must be written between a to this function and
 *  GiD_EndCoordinates().
 */

GIDPOST_API
int GiD_BeginCoordinates();

GIDPOST_API
int GiD_fBeginCoordinates(GiD_FILE fd);

/*
 *  Close the current coordinate block
 */

GIDPOST_API
int GiD_EndCoordinates();

GIDPOST_API
int GiD_fEndCoordinates(GiD_FILE fd);

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes withing the group.
 */

GIDPOST_API
int GiD_BeginMeshGroup( GP_CONST char* Name );

GIDPOST_API
int GiD_fBeginMeshGroup( GiD_FILE fd, GP_CONST char* Name );

/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */

GIDPOST_API
int GiD_EndMeshGroup();

GIDPOST_API
int GiD_fEndMeshGroup(GiD_FILE fd);

/*
 *  Write a coordinate member at the current Coordinates Block 
 */

GIDPOST_API
int GiD_WriteCoordinates(int id, double x, double y, double z);

GIDPOST_API
int GiD_fWriteCoordinates(GiD_FILE fd, int id, double x, double y, double z);

GIDPOST_API
int GiD_WriteCoordinates2D(int id, double x, double y);

GIDPOST_API
int GiD_fWriteCoordinates2D(GiD_FILE fd, int id, double x, double y);

/*
 *  Start a elements block in the current mesh
 */

GIDPOST_API
int GiD_BeginElements();

GIDPOST_API
int GiD_fBeginElements(GiD_FILE fd);

/*
 *  Close the current elements block
 */

GIDPOST_API
int GiD_EndElements();

GIDPOST_API
int GiD_fEndElements(GiD_FILE fd);

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to GiD_BeginMesh.
 *  
 */

GIDPOST_API
int GiD_WriteElement(int id, int nid[]);

GIDPOST_API
int GiD_fWriteElement(GiD_FILE fd, int id, int nid[]);

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */

GIDPOST_API
int GiD_WriteElementMat(int id, int nid[]);

GIDPOST_API
int GiD_fWriteElementMat(GiD_FILE fd, int id, int nid[]);

/*
 *  Write an sphere element member at the current Elements Block.
 *  An sphere element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the sphere element.
 *  
 */

GIDPOST_API
int GiD_WriteSphere(int id, int nid, double r);

GIDPOST_API
int GiD_fWriteSphere(GiD_FILE fd, int id, int nid, double r);

/*
 *  Write an sphere element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  An sphere element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the sphere element.
 *
 *     mat: material identification.
 *  
 */

GIDPOST_API
int GiD_WriteSphereMat(int id, int nid, double r, int mat);

GIDPOST_API
int GiD_fWriteSphereMat(GiD_FILE fd, int id, int nid, double r, int mat);

/*
 *  Write a circle element member at the current Elements Block.
 *  A circle element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the circle element.      
 *
 *     nx, ny, nz : normal to the plane containing the circle.
 *  
 */

GIDPOST_API
int GiD_WriteCircle(int id, int nid, double r,
		    double nx, double ny, double nz);

GIDPOST_API
int GiD_fWriteCircle(GiD_FILE fd, int id, int nid, double r,
		     double nx, double ny, double nz);

/*
 *  Write a circle element member at the current Elements
 *  Block. Providing also a material identification.
 *  
 *  A circle element is defined by:
 *
 *     id: element id
 *
 *     nid: node center given by the node id specified previously in
 *          the coordinate block.
 *
 *     r : radius of the circle element.      
 *
 *     nx, ny, nz : normal to the plane containing the circle.
 *  
 */

GIDPOST_API
int GiD_WriteCircleMat(int id, int nid, double r,
		       double nx, double ny, double nz, int mat);

GIDPOST_API
int GiD_fWriteCircleMat(GiD_FILE fd, int id, int nid, double r,
		        double nx, double ny, double nz, int mat);

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

/* ---------------------------------------------------------------------------
 *
 *  Post Result Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post result file. All subsequent call to functions
 *  write the information into this file. If there is no mesh file
 *  opened then the output of the mesh is writen into this file also.
 *
 */

GIDPOST_API
int GiD_OpenPostResultFile(GP_CONST char * FileName, GiD_PostMode Mode);

GIDPOST_API
GiD_FILE GiD_fOpenPostResultFile(GP_CONST char * FileName, GiD_PostMode Mode);

/*
 *  Close the current post result file
 */

GIDPOST_API
int GiD_ClosePostResultFile();

GIDPOST_API
int GiD_fClosePostResultFile(GiD_FILE fd);

/*
 *  Begin Gauss Points definition. The gauss point definition should
 *  have a name wich may be referenced in futher results blocks. The
 *  gauss points could be internal (InternalCoord=1) or given
 *  (InternalCoord=0). If the gauss points are given then the list of
 *  its natural coordinates should be written using the function
 *  GiD_WriteGaussPoint2D or GiD_WriteGaussPoint3D depending on the
 *  dimension of the element type.
 */

GIDPOST_API
int GiD_BeginGaussPoint(GP_CONST char * name, GiD_ElementType EType,
		        GP_CONST char * MeshName,
		        int GP_number, int NodesIncluded, int InternalCoord);

GIDPOST_API
int GiD_fBeginGaussPoint(GiD_FILE fd, GP_CONST char * name,
		         GiD_ElementType EType,
		         GP_CONST char * MeshName,
		         int GP_number, int NodesIncluded, int InternalCoord);

/*
 *  End current Gauss Points definition
 */

GIDPOST_API
int GiD_EndGaussPoint();

GIDPOST_API
int GiD_fEndGaussPoint(GiD_FILE fd);

/*
 *  Write internal gauss point coordinate.
 */

GIDPOST_API
int GiD_WriteGaussPoint2D(double x, double y);

GIDPOST_API
int GiD_fWriteGaussPoint2D(GiD_FILE fd, double x, double y);

GIDPOST_API
int GiD_WriteGaussPoint3D(double x, double y, double z);

GIDPOST_API
int GiD_fWriteGaussPoint3D(GiD_FILE fd, double x, double y, double z);

/*
 *  Begin a Range Table definition. With a range table you can group
 *  the result values into intervals and label each interval with a
 *  name. Inside GiD this can be visualized with a contour range. 
 */

GIDPOST_API
int GiD_BeginRangeTable(GP_CONST char * name);

GIDPOST_API
int GiD_fBeginRangeTable(GiD_FILE fd, GP_CONST char * name);

/*
 *  End a Range Table definition
 */

GIDPOST_API
int GiD_EndRangeTable();

GIDPOST_API
int GiD_fEndRangeTable(GiD_FILE fd);

/*
 *  Write Range functions --
 *
 *   WriteMinRange : write a range with an implicit minimum value, the
 *   minimum absolute in the result set.
 *
 *   WriteRange : write an explicit range.
 *
 *   WritemaxRange: write a range with an implicit maximum value, the
 *   maximum absolute in the result set.
 */

GIDPOST_API
int GiD_WriteMinRange(double max, GP_CONST char * name);

GIDPOST_API
int GiD_fWriteMinRange(GiD_FILE fd, double max, GP_CONST char * name);

GIDPOST_API
int GiD_WriteRange(double min, double max, GP_CONST char * name);

GIDPOST_API
int GiD_fWriteRange(GiD_FILE fd, double min, double max, GP_CONST char * name);

GIDPOST_API
int GiD_WriteMaxRange(double min, GP_CONST char * name);

GIDPOST_API
int GiD_fWriteMaxRange(GiD_FILE fd, double min, GP_CONST char * name);

/*
 *  Begin Result Block. This function open a result block. All the result
 *  information is provided the function call.
 */

GIDPOST_API
int GiD_BeginResult(GP_CONST char * Result, GP_CONST char * Analysis,
		    double step,
		    GiD_ResultType Type, GiD_ResultLocation Where,
		    GP_CONST char * GaussPointsName,
		    GP_CONST char * RangeTable, 
		    int compc, GP_CONST char * compv[]); 
  
GIDPOST_API
int GiD_fBeginResult(GiD_FILE fd, GP_CONST char * Result,
		     GP_CONST char * Analysis,
		     double step,
		     GiD_ResultType Type, GiD_ResultLocation Where,
		     GP_CONST char * GaussPointsName,
		     GP_CONST char * RangeTable, 
		     int compc, GP_CONST char * compv[]); 
  
/*
 *  Begin Result Block. This function open a result block. Only the
 *  result, analisys, location and location name are provided in the
 *  function call. The other result attributes as range table or
 *  components names are provided in a separated function calls.
 */

GIDPOST_API
int GiD_BeginResultHeader(GP_CONST char * Result, GP_CONST char * Analysis, 
		          double step,
		          GiD_ResultType Type, GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName);

GIDPOST_API
int GiD_fBeginResultHeader(GiD_FILE fd, GP_CONST char * Result,
		           GP_CONST char * Analysis, 
		           double step,
		           GiD_ResultType Type, GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName);

/*
 *  Define the range table associated to the current result, either a
 *  single result block or the current result defined in a result
 *  group.
 */

GIDPOST_API
int GiD_ResultRange(GP_CONST char * RangeTable);

GIDPOST_API
int GiD_fResultRange(GiD_FILE fd, GP_CONST char * RangeTable);

/*
 *  Define the components names associated to the current result,
 *  either a single result block or the current result defined in a
 *  result group.
 */

GIDPOST_API
int GiD_ResultComponents(int compc, GP_CONST char * compv[]);

GIDPOST_API
int GiD_fResultComponents(GiD_FILE fd, int compc, GP_CONST char * compv[]);

/*
 *  Define the unit string associated to the current result, either a
 *  single result block or the current result defined in a result
 *  group.
 */

GIDPOST_API
int GiD_ResultUnit(GP_CONST char * UnitName);

GIDPOST_API
int GiD_fResultUnit(GiD_FILE fd, GP_CONST char * UnitName);

/*
 *  Include property "Name" directly associated to current result
 *  with value "Value"
 */

GIDPOST_API
int GiD_ResultUserDefined(GP_CONST char * Name,GP_CONST char * Value);

GIDPOST_API
int GiD_fResultUserDefined(GiD_FILE fd, GP_CONST char * Name,GP_CONST char * Value);


/*
 *  Begin a result group. All grouped in the same analysis and step. Also the
 *  result location and location name is provided.
 */

GIDPOST_API
int GiD_BeginResultGroup(GP_CONST char * Analysis, double step,
		         GiD_ResultLocation Where,
		         GP_CONST char * GaussPointsName);

GIDPOST_API
int GiD_fBeginResultGroup(GiD_FILE fd, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName);

/*
 *  Define a result member of a result group given the name and result
 *  type. The second prototype enable us to specify the dimension of the
 *  result types. Most of the types do not allow more than one
 *  dimension. Bellow if the set of valid dimensions:
 *
 *  -  Scalar : 1 (GiD_WriteScalar)
 *  -  Vector : 2 (GiD_Write2DVector), 3(GiD_WriteVector) or 4 (GiD_WriteVectorModule)
 *  -  Matrix : 3 (GiD_Write2DMatrix) or 6 (GiD_Write3DMatrix)
 *  -  PlainDeformationMatrix : 4 (GiD_WritePlainDefMatrix)
 *  -  MainMatrix : 12 (GiD_WriteMainMatrix)
 *  -  LocalAxes : 3 (GiD_WriteLocalAxes)
 *  -  ComplexScalar : 2 ( GiD_WriteComplexScalar)
 *  -  ComplexVector : 4 ( GiD_Write2DComplexVector), 6 ( GiD_WriteComplexVector)
 */

GIDPOST_API
int GiD_ResultDescription(GP_CONST char * Result, GiD_ResultType Type);

GIDPOST_API
int GiD_fResultDescription(GiD_FILE fd, GP_CONST char * Result,
		           GiD_ResultType Type);

GIDPOST_API
int GiD_ResultDescriptionDim(GP_CONST char * Result, GiD_ResultType Type,
		              size_t dim);

GIDPOST_API
int GiD_fResultDescriptionDim(GiD_FILE fd, GP_CONST char * Result,
		              GiD_ResultType Type,
		              size_t dim);


  /*
    Associates the given local axes result to the current result.
    If result is a vector or a matrix, it means that the result is expressed in
    these local axes. Vector vx,vy,vz is not used
    If result is a scalar, local axes are a help for drawing the LineDiagrams
    the direction of the line diagram is given by vx,vy,vz expressed in Local Axes
    Analysis can be NULL and means the same analysis
   */
GIDPOST_API int GiD_ResultLocalAxes(GP_CONST char * Result,GP_CONST char * Analysis,
		                    double step,double vx,double vy,double vz);

  /*
    declares a vectorial result to be a deformation vector
   */
GIDPOST_API int GiD_ResultIsDeformationVector(int boolean);


/*
 *  This function is not needed anymore. Mark the starting point for
 *  writing the values of the current result either single or group.
 */

GIDPOST_API
int GiD_ResultValues();

GIDPOST_API
int GiD_fResultValues(GiD_FILE fd);

/*
 *  Close a previously opened result either single or group.
 */

GIDPOST_API
int GiD_EndResult();

GIDPOST_API
int GiD_fEndResult(GiD_FILE fd);

/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */

GIDPOST_API
int GiD_BeginOnMeshGroup(char * Name);

GIDPOST_API
int GiD_fBeginOnMeshGroup(GiD_FILE fd, char * Name);

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */

GIDPOST_API
int GiD_EndOnMeshGroup();

GIDPOST_API
int GiD_fEndOnMeshGroup(GiD_FILE fd);

/*
 *  Flushes all pending output into the postprocess file. This
 *  function should be called only when strictly necessary when
 *  writing in GiD_PostAsciiZipped * or GiD_PostBinary modes because
 *  it can degrade compression.
 */

GIDPOST_API
int GiD_FlushPostFile();

GIDPOST_API
int GiD_fFlushPostFile(GiD_FILE fd);

/*
 *  Write result functions
 */

GIDPOST_API
int GiD_WriteScalar(int id, double v);

GIDPOST_API
int GiD_fWriteScalar(GiD_FILE fd, int id, double v);

GIDPOST_API
int GiD_Write2DVector(int id, double x, double y);

GIDPOST_API
int GiD_fWrite2DVector(GiD_FILE fd, int id, double x, double y);

GIDPOST_API
int GiD_WriteVector(int id, double x, double y, double z);

GIDPOST_API
int GiD_fWriteVector(GiD_FILE fd, int id, double x, double y, double z);

GIDPOST_API
int GiD_WriteVectorModule(int id, double x, double y, double z, double mod);

GIDPOST_API
int GiD_fWriteVectorModule(GiD_FILE fd, int id,
		           double x, double y, double z, double mod);

GIDPOST_API
int GiD_Write2DMatrix(int id, double Sxx, double Syy, double Sxy);

GIDPOST_API
int GiD_fWrite2DMatrix(GiD_FILE fd, int id, double Sxx, double Syy, double Sxy);

GIDPOST_API
int GiD_Write3DMatrix(int id, double Sxx, double Syy, double Szz,
		      double Sxy, double Syz, double Sxz);

GIDPOST_API
int GiD_fWrite3DMatrix(GiD_FILE fd, int id, double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz);

GIDPOST_API
int GiD_WritePlainDefMatrix(int id, double Sxx, double Syy, double Sxy,
		            double Szz);

GIDPOST_API
int GiD_fWritePlainDefMatrix(GiD_FILE fd, int id,
		             double Sxx, double Syy, double Sxy, double Szz);

GIDPOST_API
int GiD_WriteMainMatrix(int id,
		        double Si, double Sii, double Siii,
		        double Vix, double Viy, double Viz,
		        double Viix, double Viiy, double Viiz,
		        double Viiix, double Viiiy, double Viiiz);
  
GIDPOST_API
int GiD_fWriteMainMatrix(GiD_FILE fd, int id,
		         double Si, double Sii, double Siii,
		         double Vix, double Viy, double Viz,
		         double Viix, double Viiy, double Viiz,
		         double Viiix, double Viiiy, double Viiiz);
  
GIDPOST_API
int GiD_WriteLocalAxes(int id, double euler_1, double euler_2, double euler_3);

GIDPOST_API
int GiD_fWriteLocalAxes(GiD_FILE fd, int id,
		        double euler_1, double euler_2, double euler_3);

GIDPOST_API
int GiD_WriteComplexScalar( int id, double complex_real, double complex_imag);

GIDPOST_API
int GiD_fWriteComplexScalar( GiD_FILE fd, int id, double complex_real, double complex_imag);

GIDPOST_API
int GiD_Write2DComplexVector( int id,
                              double x_real, double x_imag,
                              double y_real, double y_imag);

GIDPOST_API
int GiD_fWrite2DComplexVector( GiD_FILE fd, int id,
                               double x_real, double x_imag,
                               double y_real, double y_imag);

GIDPOST_API
int GiD_WriteComplexVector( int id,
                            double x_real, double x_imag,
                            double y_real, double y_imag,
                            double z_real, double z_imag);

GIDPOST_API
int GiD_fWriteComplexVector( GiD_FILE fd, int id,
                             double x_real, double x_imag,
                             double y_real, double y_imag,
                             double z_real, double z_imag);

#ifdef __cplusplus
}
#endif

#endif
