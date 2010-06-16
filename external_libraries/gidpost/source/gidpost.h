/* gidpost 1.7 */
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

#define USE_CONST
#if defined(USE_CONST)
#define GP_CONST const
#else
#define GP_CONST /* empty */
#endif

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

#ifdef __cplusplus
extern "C" {
#endif

/* opening mode */

typedef enum { 
  GiD_PostAscii, 
  GiD_PostAsciiZipped, 
  GiD_PostBinary } GiD_PostMode;

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
  GiD_Prism
} GiD_ElementType;

typedef enum {
  GiD_Scalar = 0,
  GiD_Vector,
  GiD_Matrix,
  GiD_PlainDeformationMatrix,
  GiD_MainMatrix,
  GiD_LocalAxes
} GiD_ResultType;

typedef enum { GiD_OnNodes, GiD_OnGaussPoints } GiD_ResultLocation;

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

/*
 *  Close the current post mesh file
 */

GIDPOST_API
int GiD_ClosePostMeshFile();

/*
 *  Begin a new mesh. After that you should write the nodes and the
 *  elements. When writing the elements it is assumed a connectivity
 *  of size given in the paramenter NNode.
 */

GIDPOST_API
int GiD_BeginMesh( GP_CONST char * MeshName,
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
int GiD_BeginMeshColor( GP_CONST char * MeshName,
			GiD_Dimension Dim, GiD_ElementType EType,
			int NNode,
			double Red, double Green, double Blue );

/*
 *  End current mesh. 
 */

GIDPOST_API
int GiD_EndMesh();

/*
 *  Start a coordinate block in the current mesh. All nodes
 *  coordinates must be written between a to this function and
 *  GiD_EndCoordinates().
 */

GIDPOST_API
int GiD_BeginCoordinates();

/*
 *  Close the current coordinate block
 */

GIDPOST_API
int GiD_EndCoordinates();

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes withing the group.
 */

GIDPOST_API
int GiD_BeginMeshGroup(GP_CONST char* Name);

/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */

GIDPOST_API
int GiD_EndMeshGroup();

/*
 *  Write a coordinate member at the current Coordinates Block 
 */

GIDPOST_API
int GiD_WriteCoordinates(int id, double x, double y, double z);

GIDPOST_API
int GiD_WriteCoordinates2D(int id, double x, double y);

/*
 *  Start a elements block in the current mesh
 */

GIDPOST_API
int GiD_BeginElements();

/*
 *  Close the current elements block
 */

GIDPOST_API
int GiD_EndElements();

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to GiD_BeginMesh.
 *  
 */

GIDPOST_API
int GiD_WriteElement( int id, int nid[] );

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */

GIDPOST_API
int GiD_WriteElementMat( int id, int nid[] );

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
int GiD_OpenPostResultFile( GP_CONST char * FileName, GiD_PostMode Mode );

/*
 *  Close the current post result file
 */

GIDPOST_API
int GiD_ClosePostResultFile();

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
int GiD_BeginGaussPoint( GP_CONST char * name, GiD_ElementType EType,
			 GP_CONST char * MeshName,
			 int GP_number, int NodesIncluded, int InternalCoord );

/*
 *  End current Gauss Points definition
 */

GIDPOST_API
int GiD_EndGaussPoint();

/*
 *  Write internal gauss point coordinate.
 */

GIDPOST_API
int GiD_WriteGaussPoint2D( double x, double y );

GIDPOST_API
int GiD_WriteGaussPoint3D( double x, double y, double z );

/*
 *  Begin a Range Table definition. With a range table you can group
 *  the result values into intervals and label each interval with a
 *  name. Inside GiD this can be visualized with a contour range. 
 */

GIDPOST_API
int GiD_BeginRangeTable( GP_CONST char * name );

/*
 *  End a Range Table definition
 */

GIDPOST_API
int GiD_EndRangeTable();

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
int GiD_WriteMinRange( double max, GP_CONST char * name );

GIDPOST_API
int GiD_WriteRange( double min, double max, GP_CONST char * name );

GIDPOST_API
int GiD_WriteMaxRange(double min, GP_CONST char * name);

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

/*
 *  Define the range table associated to the current result, either a
 *  single result block or the current result defined in a result
 *  group.
 */

GIDPOST_API
int GiD_ResultRange( GP_CONST char * RangeTable);

/*
 *  Define the components names associated to the current result,
 *  either a single result block or the current result defined in a
 *  result group.
 */

GIDPOST_API
int GiD_ResultComponents( int compc, GP_CONST char * compv[]);

/*
 *  Define the unit string associated to the current result, either a
 *  single result block or the current result defined in a result
 *  group.
 */

/* This function is not supported yet inside GiD */
/* int GiD_ResultUnit(char * UnitName); */

/*
 *  Begin a result group. All grouped in the same analysis and step. Also the
 *  result location and location name is provided.
 */

GIDPOST_API
int GiD_BeginResultGroup(GP_CONST char * Analysis, double step,
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
 */

GIDPOST_API
int GiD_ResultDescription(GP_CONST char * Result, GiD_ResultType Type);

GIDPOST_API
int GiD_ResultDescriptionDim(GP_CONST char * Result, GiD_ResultType Type,
			      size_t dim);

/*
 *  This function is not needed anymore. Mark the starting point for
 *  writing the values of the current result either single or group.
 */

GIDPOST_API
int GiD_ResultValues();

/*
 *  Close a previously opened result either single or group.
 */

GIDPOST_API
int GiD_EndResult();

/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */

GIDPOST_API
int GiD_BeginOnMeshGroup(char * Name);

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */

GIDPOST_API
int GiD_EndOnMeshGroup();

/*
 *  Flushes all pending output into the postprocess file. This function should
 *  be called only when strictly necessary when writing in GiD_PostAsciiZipped *
 *  or GiD_PostBinary modes because it can degrade compression.
 */

GIDPOST_API
int GiD_FlushPostFile();

/*
 *  Write result functions
 */

GIDPOST_API
int GiD_WriteScalar(int id, double v);

GIDPOST_API
int GiD_Write2DVector(int id, double x, double y);

GIDPOST_API
int GiD_WriteVector(int id, double x, double y, double z);

GIDPOST_API
int GiD_WriteVectorModule(int id, double x, double y, double z, double mod);

GIDPOST_API
int GiD_Write2DMatrix(int id, double Sxx, double Syy, double Sxy);

GIDPOST_API
int GiD_Write3DMatrix(int id, double Sxx, double Syy, double Szz,
                      double Sxy, double Syz, double Sxz);

GIDPOST_API
int GiD_WritePlainDefMatrix(int id, double Sxx, double Syy, double Sxy,
			    double Szz);

GIDPOST_API
int GiD_WriteMainMatrix(int id,
                        double Si, double Sii, double Siii,
                        double Vix, double Viy, double Viz,
                        double Viix, double Viiy, double Viiz,
                        double Viiix, double Viiiy, double Viiiz);
  
GIDPOST_API
int GiD_WriteLocalAxes(int id, double euler_1, double euler_2, double euler_3);

#ifdef __cplusplus
}
#endif

#endif
