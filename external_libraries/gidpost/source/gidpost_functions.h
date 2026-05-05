#ifndef __GIDPOST_FUNCTIONS__
#define __GIDPOST_FUNCTIONS__

#include "gidpost_types.h"

#include "gidpost_cluster_functions.h"

/*
  GiD_PostInit -- Initialization of gidpost library, must be called
  from the main thread of the program.
 */
GIDPOST_API int GiD_PostInit( void );
GIDPOST_API int GiD_PostDone( void );
GIDPOST_API GP_CONST char *GiD_PostGetVersion( void );
GIDPOST_API int GiD_PostIsThreadSafe( GiD_PostMode Mode ); // returns -1 on error
  
GIDPOST_API int GiD_PostSetFormatReal(GP_CONST char * format_real);
GIDPOST_API GP_CONST char *GiD_PostGetFormatReal( void );
GIDPOST_API GP_CONST char *GiD_PostGetFormatStep( void );

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
GiD_FILE GiD_fOpenPostMeshFile(GP_CONST char * FileName, GiD_PostMode Mode );
GIDPOST_API
GiD_FILE GiD_fOpenPostMeshFile_utf8(GP_CONST char * FileName, GiD_PostMode Mode );

/*
 *  Close the current post mesh file
 */

GIDPOST_API
int GiD_fClosePostMeshFile(GiD_FILE fd);

/*
 *  Begin a new mesh. After that you should write the nodes and the
 *  elements. When writing the elements it is assumed a connectivity
 *  of size given in the parameter NNode.
 */
GIDPOST_API
int GiD_fBeginMesh(GiD_FILE fd, GP_CONST char * MeshName,
		   GiD_Dimension Dim, GiD_ElementType EType,
		   int NNode );

/*
 *  Begin a new mesh. After that you should write the nodes and the
 *  elements. When writing the elements it is assumed a connectivity
 *  of size given in the parameter NNode. With this function you can
 *  specify a color for the mesh by its RGB components where each
 *  component take values on the interval [0,1].
 */
GIDPOST_API
int GiD_fBeginMeshColor(GiD_FILE fd, GP_CONST char * MeshName,
		        GiD_Dimension Dim, GiD_ElementType EType,
		        int NNode,
		        double Red, double Green, double Blue);
/*
 *  End current mesh. 
 */
GIDPOST_API
int GiD_fEndMesh(GiD_FILE fd);

/* 
 *    Declares the units of the mesh
 */
GIDPOST_API
int GiD_fMeshUnit(GiD_FILE fd,GP_CONST char * UnitName);

  /* 
     Associates a local axes result to the current mesh
   */
// Only used for Compassis files
GIDPOST_API int GiD_fMeshLocalAxes(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis,
		              double step);

/*
 *  Start a coordinate block in the current mesh. All nodes
 *  coordinates must be written between a to this function and
 *  GiD_EndCoordinates( void ).
 */
GIDPOST_API
int GiD_fBeginCoordinates(GiD_FILE fd);

/*
 *  Close the current coordinate block
 */
GIDPOST_API
int GiD_fEndCoordinates(GiD_FILE fd);

/*
 * This function open a group of mesh. This makes possible specifying
 * multiples meshes within the group.
 */
GIDPOST_API
int GiD_fBeginMeshGroup( GiD_FILE fd, GP_CONST char* Name );

/*
 * This function close the previously opened group of mesh. See
 * GiD_BeginMeshGroup.
 */
GIDPOST_API
int GiD_fEndMeshGroup(GiD_FILE fd);

/*
 *  Write a coordinate member at the current Coordinates Block 
 */
GIDPOST_API
int GiD_fWriteCoordinates(GiD_FILE fd, int id, double x, double y, double z);

GIDPOST_API
int GiD_fWriteCoordinates2D(GiD_FILE fd, int id, double x, double y);

// GiD_fWriteCoordinatesBlock includes BeginCoordinates() and EndCoordinates()
GIDPOST_API
int GiD_fWriteCoordinatesBlock( GiD_FILE fd, int num_points, GP_CONST double *xyz_array );

GIDPOST_API
int GiD_fWriteCoordinatesIdBlock( GiD_FILE fd, int num_points, GP_CONST int *list_ids, GP_CONST double *xyz_array );

/*
 *  Start a elements block in the current mesh
 */
GIDPOST_API
int GiD_fBeginElements(GiD_FILE fd);

/*
 *  Close the current elements block
 */
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
int GiD_fWriteElement(GiD_FILE fd, int id, int nid[]);

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *  
 */
GIDPOST_API
int GiD_fWriteElementMat(GiD_FILE fd, int id, int nid[]);

// GiD_fWriteElementsBlock includes BeginElements() and EndElements()
GIDPOST_API
int GiD_fWriteElementsBlock( GiD_FILE fd, int num_elements, GP_CONST int *connectivities );
GIDPOST_API
int GiD_fWriteElementsIdBlock( GiD_FILE fd, int num_elements, GP_CONST int *list_ids, GP_CONST int *connectivities );
GIDPOST_API
int GiD_fWriteElementsMatBlock( GiD_FILE fd, int num_elements, GP_CONST int *connectivities,
                                GP_CONST int *lst_material_id );
GIDPOST_API
int GiD_fWriteElementsIdMatBlock( GiD_FILE fd, int num_elements, GP_CONST int *list_ids, GP_CONST int *connectivities,
                                  GP_CONST int *lst_material_id );

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
int GiD_fWriteCircleMat(GiD_FILE fd, int id, int nid, double r,
		        double nx, double ny, double nz, int mat);

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
GiD_FILE GiD_fOpenPostResultFile(GP_CONST char * FileName, GiD_PostMode Mode);
GIDPOST_API
GiD_FILE GiD_fOpenPostResultFile_utf8(GP_CONST char * FileName, GiD_PostMode Mode);

/*
 *  Close the current post result file
 */
GIDPOST_API
int GiD_fClosePostResultFile(GiD_FILE fd);

/*
 *  Begin Gauss Points definition. The gauss point definition should
 *  have a name which may be referenced in further results blocks. The
 *  gauss points could be internal (InternalCoord=1) or given
 *  (InternalCoord=0). If the gauss points are given then the list of
 *  its natural coordinates should be written using the function
 *  GiD_WriteGaussPoint2D or GiD_WriteGaussPoint3D depending on the
 *  dimension of the element type.
 */
GIDPOST_API
int GiD_fBeginGaussPoint(GiD_FILE fd, GP_CONST char * name,
		         GiD_ElementType EType,
		         GP_CONST char * MeshName,
		         int GP_number, int NodesIncluded, int InternalCoord);

/*
 *  End current Gauss Points definition
 */
GIDPOST_API
int GiD_fEndGaussPoint(GiD_FILE fd);

/*
 *  Write internal gauss point coordinate.
 */
GIDPOST_API
int GiD_fWriteGaussPoint2D(GiD_FILE fd, double x, double y);
GIDPOST_API
int GiD_fWriteGaussPoint3D(GiD_FILE fd, double x, double y, double z);

/*
 *  Begin a Range Table definition. With a range table you can group
 *  the result values into intervals and label each interval with a
 *  name. Inside GiD this can be visualized with a contour range. 
 */
GIDPOST_API
int GiD_fBeginRangeTable(GiD_FILE fd, GP_CONST char * name);

/*
 *  End a Range Table definition
 */
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
int GiD_fWriteMinRange(GiD_FILE fd, double max, GP_CONST char * name);
GIDPOST_API
int GiD_fWriteRange(GiD_FILE fd, double min, double max, GP_CONST char * name);
GIDPOST_API
int GiD_fWriteMaxRange(GiD_FILE fd, double min, GP_CONST char * name);


// DEPRECATED, use GiD_fBeginResultHeader
// GIDPOST_API_DEPRECATED
GIDPOST_API
int // GCC_GIDPOST_API_DEPRECATED 
GiD_fBeginResult(GiD_FILE fd, GP_CONST char * Result,
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
int GiD_fResultRange(GiD_FILE fd, GP_CONST char * RangeTable);

/*
 *  Define the components names associated to the current result,
 *  either a single result block or the current result defined in a
 *  result group.
 */
GIDPOST_API
int GiD_fResultComponents(GiD_FILE fd, int compc, GP_CONST char * compv[]);

/*
 *  Define the unit string associated to the current result, either a
 *  single result block or the current result defined in a result
 *  group.
 */
GIDPOST_API
int GiD_fResultUnit(GiD_FILE fd, GP_CONST char * UnitName);

/*
 *  Include property "Name" directly associated to current result
 *  with value "Value"
 */
GIDPOST_API
int GiD_fResultUserDefined(GiD_FILE fd, GP_CONST char * Name,GP_CONST char * Value);


/*
 *  Begin a result group. All grouped in the same analysis and step. Also the
 *  result location and location name is provided.
 */
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
 *  -  ComplexMatrix : 6 ( GiD_Write2DComplexMatrix), 12 ( GiD_WriteComplexMatrix)
 */
GIDPOST_API
int GiD_fResultDescription(GiD_FILE fd, GP_CONST char * Result,
		           GiD_ResultType Type);

GIDPOST_API
int GiD_fResultDescriptionDim(GiD_FILE fd, GP_CONST char * Result,
		              GiD_ResultType Type,
		              int dim);

  /*
    Associates the given local axes result to the current result.
    If result is a vector or a matrix, it means that the result is expressed in
    these local axes. Vector vx,vy,vz is not used
    If result is a scalar, local axes are a help for drawing the LineDiagrams
    the direction of the line diagram is given by vx,vy,vz expressed in Local Axes
    Analysis can be NULL and means the same analysis
   */
// Only used for Compassis files
GIDPOST_API int GiD_fResultLocalAxes(GiD_FILE fd,GP_CONST char * Result,GP_CONST char * Analysis,
		                    double step,double vx,double vy,double vz);

  /*
    declares a vectorial result to be a deformation vector
   */
// Only used for Compassis files
GIDPOST_API int GiD_fResultIsDeformationVector(GiD_FILE fd,int boolean);

/*
 *  This function is not needed anymore. Mark the starting point for
 *  writing the values of the current result either single or group.
 */
GIDPOST_API
int GiD_fResultValues(GiD_FILE fd);

/*
 *  Close a previously opened result either single or group.
 */
GIDPOST_API
int GiD_fEndResult(GiD_FILE fd);


// GiD_fWriteResultBlock includes BeginResult()/BeginValues() and EndValues()/EndResult()
GIDPOST_API
int GiD_fWriteResultBlock( GiD_FILE fd, GP_CONST char *result_name, GP_CONST char *analysis_name, double step_value,
                           GiD_ResultType result_type, GiD_ResultLocation result_location,
                           GP_CONST char *gauss_point_name, GP_CONST char *range_table_name, int num_component_names,
                           GP_CONST char **list_component_names, 
			   GP_CONST char *unit_name,
			   int num_result_values, GP_CONST int *list_result_ids,
			   // list_component_values === 
			   //   num_component_values == 1 --> scalar_array
			   //   num_component_values == 2 --> VxVyVz_array
			   //   num_component_values == 6 --> SxxSyySxySzzSxzSyz_array
			   //   ...
                           int num_component_values, GP_CONST double *list_component_values );


/*
 * Results which are referred to a mesh group (see GiD_BeginMeshGroup)
 * should be written between a call to this function and
 * GiD_EndOnMeshGroup.
 */
GIDPOST_API
int GiD_fBeginOnMeshGroup(GiD_FILE fd, char * Name);

/*
 * This function close a previously opened block of result over a mesh
 * group.
 */
GIDPOST_API
int GiD_fEndOnMeshGroup(GiD_FILE fd);

/*
 *  Flushes all pending output into the postprocess file. This
 *  function should be called only when strictly necessary when
 *  writing in GiD_PostAsciiZipped * or GiD_PostBinary modes because
 *  it can degrade compression.
 */
GIDPOST_API
int GiD_fFlushPostFile(GiD_FILE fd);

/*
 *  Write result functions
 */
GIDPOST_API
int GiD_fWriteScalar(GiD_FILE fd, int id, double v);
GIDPOST_API
int GiD_fWrite2DVector(GiD_FILE fd, int id, double x, double y);
GIDPOST_API
int GiD_fWriteVector(GiD_FILE fd, int id, double x, double y, double z);
GIDPOST_API
int GiD_fWriteVectorModule(GiD_FILE fd, int id,
		           double x, double y, double z, double mod);
GIDPOST_API
int GiD_fWrite2DMatrix(GiD_FILE fd, int id, double Sxx, double Syy, double Sxy);
GIDPOST_API
int GiD_fWrite3DMatrix(GiD_FILE fd, int id, double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz);
GIDPOST_API
int GiD_fWritePlainDefMatrix(GiD_FILE fd, int id,
		             double Sxx, double Syy, double Sxy, double Szz);
GIDPOST_API
int GiD_fWriteMainMatrix(GiD_FILE fd, int id,
		         double Si, double Sii, double Siii,
		         double Vix, double Viy, double Viz,
		         double Viix, double Viiy, double Viiz,
		         double Viiix, double Viiiy, double Viiiz);
GIDPOST_API
int GiD_fWriteLocalAxes(GiD_FILE fd, int id,
		        double euler_1, double euler_2, double euler_3);
GIDPOST_API
int GiD_fWriteComplexScalar( GiD_FILE fd, int id, double complex_real, double complex_imag);
GIDPOST_API
int GiD_fWrite2DComplexVector( GiD_FILE fd, int id,
                               double x_real, double x_imag,
                               double y_real, double y_imag);
GIDPOST_API
int GiD_fWriteComplexVector( GiD_FILE fd, int id,
                             double x_real, double x_imag,
                             double y_real, double y_imag,
                             double z_real, double z_imag);
GIDPOST_API
int GiD_fWrite2DComplexMatrix(GiD_FILE fd, int id,
                              double Sxx_real, double Syy_real, double Sxy_real,
                              double Sxx_imag, double Syy_imag, double Sxy_imag);
GIDPOST_API
int GiD_fWrite3DComplexMatrix(GiD_FILE fd, int id,
                              double Sxx_real, double Syy_real, double Szz_real,
                              double Sxy_real, double Syz_real, double Sxz_real,
                              double Sxx_imag, double Syy_imag, double Szz_imag,
                              double Sxy_imag, double Syz_imag, double Sxz_imag);
GIDPOST_API
int GiD_fWriteNurbsSurface( GiD_FILE fd, int id, int n, double* v );
GIDPOST_API
int GiD_fWriteNurbsSurfaceVector( GiD_FILE fd, int id, int n, int num_comp, double* v );

/* User defined properties defined inside Mesh or Result blocks
   ENABLE_HDF5: stored as properties/attributes (Name, value) of the current Mesh/N or Result/N folder
   ASCII / raw binary: stored as comments
     # Name: value
   Define the macro COMPASSIS_USER_ATTRIBUTES_FORMAT
   to have it like compassis wants:
     # ResultUserDefined \"%s\" \"%s\"      or
     # ResultUserDefined \"%s\" %s
*/
/*
 *  Include property "Name" directly associated to current mesh
 *  with value "Value"
 */
GIDPOST_API
int GiD_fWriteMeshUserAttribute( GiD_FILE fd, GP_CONST char *Name, GP_CONST char *Value );
GIDPOST_API
int GiD_fWriteResultUserAttribute( GiD_FILE fd, GP_CONST char *Name, GP_CONST char *Value );

#endif // __GIDPOST_FUNCTIONS__
