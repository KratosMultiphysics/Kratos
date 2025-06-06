#pragma once

#ifdef ENABLE_HDF5

#define MAX_CONCURRENT_DATASETS 20

typedef struct _CurrentHdf5WriteData CurrentHdf5WriteData;

const char *GiD_GetHDF5Version( void );
int GiD_IsThreadSafe_HDF5( void );

CurrentHdf5WriteData *new_CurrentHdf5WriteData( void );
void delete_CurrentHdf5WriteData( CurrentHdf5WriteData *obj );

/*################################################################################
*#    Mesh file
*################################################################################*/

int GiD_OpenPostMeshFile_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * FileName);
int GiD_ClosePostMeshFile_HDF5( CurrentHdf5WriteData *obj );

/*################################################################################
*#    Write Mesh
*################################################################################*/

int GiD_BeginMesh_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,
		       int NNode);
int GiD_BeginMeshColor_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,
		            int NNode,double Red, double Green, double Blue);
int GiD_EndMesh_HDF5( CurrentHdf5WriteData *obj );
int GiD_MeshUnit_HDF5( CurrentHdf5WriteData *obj,  GP_CONST char * UnitName );
int GiD_MeshLocalAxes_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis,double step);
int GiD_BeginCoordinates_HDF5( CurrentHdf5WriteData *obj );
int GiD_EndCoordinates_HDF5( CurrentHdf5WriteData *obj );
int GiD_WriteCoordinates_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y, double z);
int GiD_WriteCoordinates2D_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y);
// GiD_WriteCoordinatesBlock_HDF5 includes BeginCoordinates and EndCoordinates
int GiD_WriteCoordinatesBlock_HDF5( CurrentHdf5WriteData *obj, int num_points, const double *xyz_array);
int GiD_WriteCoordinatesIdBlock_HDF5( CurrentHdf5WriteData *obj, int num_points, const int *list_ids, const double *xyz_array );
int GiD_BeginElements_HDF5( CurrentHdf5WriteData *obj );
int GiD_EndElements_HDF5( CurrentHdf5WriteData *obj );
int GiD_WriteElement_HDF5( CurrentHdf5WriteData *obj, int id, int nid[]);
int GiD_WriteElementMat_HDF5( CurrentHdf5WriteData *obj, int id, int nid[]);
// GiD_WriteElementsBlock_HDF5 includes BeginElements and EndElements
int GiD_WriteElementsBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *connectivities );
int GiD_WriteElementsIdBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *list_ids, const int *connectivities );
int GiD_WriteElementsMatBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *connectivities, const int *lst_material_id );
int GiD_WriteElementsIdMatBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *list_ids, const int *connectivities, const int *lst_material_id );
int GiD_WriteSphere_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r);
int GiD_WriteSphereMat_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r, int mat);
int GiD_WriteCircle_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r,double nx, double ny, double nz);
int GiD_WriteCircleMat_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r,double nx, double ny, double nz, int mat);

/*################################################################################
*#    Results file
*################################################################################*/

int GiD_OpenPostResultFile_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * FileName);
int GiD_ClosePostResultFile_HDF5( CurrentHdf5WriteData *obj );

/*################################################################################
*#    Write Gauss points
*################################################################################*/

int GiD_BeginGaussPoint_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * name, GiD_ElementType EType,GP_CONST char * MeshName,
  int GP_number, int NodesIncluded, int InternalCoord);
int GiD_EndGaussPoint_HDF5( CurrentHdf5WriteData *obj );
int GiD_WriteGaussPoint2D_HDF5( CurrentHdf5WriteData *obj, double x, double y);
int GiD_WriteGaussPoint3D_HDF5( CurrentHdf5WriteData *obj, double x, double y, double z);

/*################################################################################
*#    Write range table
*################################################################################*/

int GiD_BeginRangeTable_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * name);
int GiD_EndRangeTable_HDF5( CurrentHdf5WriteData *obj );
int GiD_WriteMinRange_HDF5( CurrentHdf5WriteData *obj, double max, GP_CONST char * name);
int GiD_WriteRange_HDF5( CurrentHdf5WriteData *obj, double min, double max, GP_CONST char * name);
int GiD_WriteMaxRange_HDF5( CurrentHdf5WriteData *obj, double min, GP_CONST char * name);

/*################################################################################
*#    Write results
*################################################################################*/

int GiD_BeginResult_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis,double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName,
  GP_CONST char * RangeTable, 
  int compc, GP_CONST char * compv[]);
int GiD_BeginResultHeader_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis, double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName);
int GiD_BeginResultGroup_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Analysis, double step,GiD_ResultLocation Where,
			      GP_CONST char * GaussPointsName);
int GiD_ResultDescription_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GiD_ResultType Type);

/*
    Associates the given local axes result to the current result.
    If result is a vector or a matrix, it means that the result is expressed in
    these local axes. Vector vx,vy,vz is not used
    If result is a scalar, local axes are a help for drawing the LineDiagrams
    the direction of the line diagram is given by vx,vy,vz expressed in Local Axes
    Analysis can be NULL and means the same analysis
*/
int GiD_ResultLocalAxes_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis,
			     double step,double vx,double vy,double vz);
int GiD_ResultIsDeformationVector_HDF5( CurrentHdf5WriteData *obj, int boolean);
int GiD_ResultRange_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * RangeTable);
int GiD_ResultComponents_HDF5( CurrentHdf5WriteData *obj, int compc, GP_CONST char * compv[]);
int GiD_ResultUnit_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * UnitName);
int GiD_ResultUserDefined_HDF5( CurrentHdf5WriteData *obj, GP_CONST char* Name,GP_CONST char* Value);
int GiD_ResultValues_HDF5( CurrentHdf5WriteData *obj );
int GiD_EndResult_HDF5( CurrentHdf5WriteData *obj );
int GiD_FlushPostFile_HDF5( CurrentHdf5WriteData *obj );


int GiD_WriteScalar_HDF5( CurrentHdf5WriteData *obj, int id, double v);
int GiD_Write2DVector_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y);
int GiD_WriteVector_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y, double z);
int GiD_WriteVectorModule_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y, double z, double mod);
int GiD_Write2DMatrix_HDF5( CurrentHdf5WriteData *obj, int id, double Sxx, double Syy, double Sxy);
int GiD_Write3DMatrix_HDF5( CurrentHdf5WriteData *obj, int id, double Sxx, double Syy, double Szz,
  double Sxy, double Syz, double Sxz);
int GiD_WritePlainDefMatrix_HDF5( CurrentHdf5WriteData *obj, int id, double Sxx, double Syy, double Sxy,
  double Szz);
int GiD_WriteMainMatrix_HDF5( CurrentHdf5WriteData *obj, int id,
  double Si, double Sii, double Siii,
  double Vix, double Viy, double Viz,
  double Viix, double Viiy, double Viiz,
  double Viiix, double Viiiy, double Viiiz);
int GiD_WriteLocalAxes_HDF5( CurrentHdf5WriteData *obj, int id, double euler_1, double euler_2, double euler_3);
/* complex results */
int GiD_WriteComplexScalar_HDF5( CurrentHdf5WriteData *obj, int id, double complex_real, double complex_imag);
int GiD_Write2DComplexVector_HDF5( CurrentHdf5WriteData *obj, int id, 
				  double x_real, double x_imag,
				  double y_real, double y_imag);
int GiD_WriteComplexVector_HDF5( CurrentHdf5WriteData *obj, int id, 
				double x_real, double x_imag,
				double y_real, double y_imag,
				double z_real, double z_imag);
int GiD_Write2DComplexMatrix_HDF5( CurrentHdf5WriteData *obj, int id,
                                  double Sxx_real, double Syy_real, double Sxy_real,
                                  double Sxx_imag, double Syy_imag, double Sxy_imag);
int GiD_WriteComplexMatrix_HDF5( CurrentHdf5WriteData *obj, int id,
                                double Sxx_real, double Syy_real, double Szz_real,
                                double Sxy_real, double Syz_real, double Sxz_real,
                                double Sxx_imag, double Syy_imag, double Szz_imag,
                                double Sxy_imag, double Syz_imag, double Sxz_imag);
/* OnNurbsSurface results */
int GiD_WriteNurbsSurface_HDF5( CurrentHdf5WriteData *obj, int id, int n, double* v);//scalar
int GiD_WriteNurbsSurfaceVector_HDF5( CurrentHdf5WriteData *obj,  int id, int n, int num_comp, double* v );//vector3

/* MeshGroups: */
int GiD_BeginMeshGroup_HDF5( CurrentHdf5WriteData *obj,  const char *Name);
int GiD_EndMeshGroup_HDF5( CurrentHdf5WriteData *obj );
int GiD_BeginOnMeshGroup_HDF5( CurrentHdf5WriteData *obj,  const char *Name);
int GiD_EndOnMeshGroup_HDF5( CurrentHdf5WriteData *obj );

/* User defined properties defined inside Mesh or Result blocks
   HDF5: stored as properties/attributes of the Mesh/Result folder (Name, value)
   ASCII / raw binary: stored as comments
     # Name: value
   Define the macro COMPASSIS_USER_ATTRIBUTES_FORMAT
   to have it like compassis wants:
     # ResultUserDefined \"%s\" \"%s\"      or
     # ResultUserDefined \"%s\" %s
*/
int GiD_WriteMeshUserAttribute_HDF5( CurrentHdf5WriteData *obj,  GP_CONST char *Name, GP_CONST char *Value);
int GiD_WriteResultUserAttribute_HDF5( CurrentHdf5WriteData *obj,  GP_CONST char *Name, GP_CONST char *Value);

// GiD_fWriteResultBlock includes BeginResult()/BeginValues() and EndValues()/EndResult()
int GiD_WriteResultBlock_HDF5( CurrentHdf5WriteData *obj, GP_CONST char *result_name, GP_CONST char *analysis_name,
                               double step_value, GiD_ResultType result_type, GiD_ResultLocation result_location,
                               GP_CONST char *gauss_point_name, GP_CONST char *range_table_name,
                               int num_component_names, GP_CONST char *list_component_names[], GP_CONST char *unit_name, 
                               int num_result_values,
                               GP_CONST int *list_result_ids, int num_component_values,
                               GP_CONST double *list_component_values );

#endif // ENABLE_HDF5
