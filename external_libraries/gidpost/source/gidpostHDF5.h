
#ifdef HDF5

#define MAX_CONCURRENT_DATASETS 20


/*################################################################################
*#    Mesh file
*################################################################################*/

int GiD_OpenPostMeshFile_HDF5(GP_CONST char * FileName);
int GiD_ClosePostMeshFile_HDF5();

/*################################################################################
*#    Write Mesh
*################################################################################*/

int GiD_BeginMesh_HDF5(GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,
		       int NNode);
int GiD_BeginMeshColor_HDF5(GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,
		            int NNode,double Red, double Green, double Blue);
int GiD_EndMesh_HDF5();
int GiD_MeshUnit_HDF5();
int GiD_MeshLocalAxes_HDF5(GP_CONST char * Result, GP_CONST char * Analysis,double step);
int GiD_BeginCoordinates_HDF5();
int GiD_EndCoordinates_HDF5();
int GiD_WriteCoordinates_HDF5(int id, double x, double y, double z);
int GiD_WriteCoordinates2D_HDF5(int id, double x, double y);
int GiD_BeginElements_HDF5();
int GiD_EndElements_HDF5();
int GiD_WriteElement_HDF5(int id, int nid[]);
int GiD_WriteElementMat_HDF5(int id, int nid[]);
int GiD_WriteSphere_HDF5(int id, int nid, double r);
int GiD_WriteSphereMat_HDF5(int id, int nid, double r, int mat);
int GiD_WriteCircle_HDF5(int id, int nid, double r,double nx, double ny, double nz);
int GiD_WriteCircleMat_HDF5(int id, int nid, double r,double nx, double ny, double nz, int mat);
int GiD_WriteCluster_HDF5(int id, int nid);
int GiD_WriteClusterMat_HDF5(int id, int nid, int mat);

/*################################################################################
*#    Results file
*################################################################################*/

int GiD_OpenPostResultFile_HDF5(GP_CONST char * FileName);
int GiD_ClosePostResultFile_HDF5();

/*################################################################################
*#    Write Gauss points
*################################################################################*/

int GiD_BeginGaussPoint_HDF5(GP_CONST char * name, GiD_ElementType EType,GP_CONST char * MeshName,
  int GP_number, int NodesIncluded, int InternalCoord);
int GiD_EndGaussPoint_HDF5();
int GiD_WriteGaussPoint2D_HDF5(double x, double y);
int GiD_WriteGaussPoint3D_HDF5(double x, double y, double z);

/*################################################################################
*#    Write range table
*################################################################################*/

int GiD_BeginRangeTable_HDF5(GP_CONST char * name);
int GiD_EndRangeTable_HDF5();
int GiD_WriteMinRange_HDF5(double max, GP_CONST char * name);
int GiD_WriteRange_HDF5(double min, double max, GP_CONST char * name);
int GiD_WriteMaxRange_HDF5(double min, GP_CONST char * name);

/*################################################################################
*#    Write results
*################################################################################*/

int GiD_BeginResult_HDF5(GP_CONST char * Result, GP_CONST char * Analysis,double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName,
  GP_CONST char * RangeTable, 
  int compc, GP_CONST char * compv[]);
int GiD_BeginResultHeader_HDF5(GP_CONST char * Result, GP_CONST char * Analysis, double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName);
int GiD_BeginResultGroup_HDF5(GP_CONST char * Analysis, double step,GiD_ResultLocation Where,
			      GP_CONST char * GaussPointsName);
int GiD_ResultDescription_HDF5(GP_CONST char * Result, GiD_ResultType Type);

/*
    Associates the given local axes result to the current result.
    If result is a vector or a matrix, it means that the result is expressed in
    these local axes. Vector vx,vy,vz is not used
    If result is a scalar, local axes are a help for drawing the LineDiagrams
    the direction of the line diagram is given by vx,vy,vz expressed in Local Axes
    Analysis can be NULL and means the same analysis
*/
int GiD_ResultLocalAxes_HDF5(GP_CONST char * Result, GP_CONST char * Analysis,
			     double step,double vx,double vy,double vz);
int GiD_ResultIsDeformationVector_HDF5(int boolean);
int GiD_ResultRange_HDF5(GP_CONST char * RangeTable);
int GiD_ResultComponents_HDF5(int compc, GP_CONST char * compv[]);
int GiD_ResultUnit_HDF5(GP_CONST char * UnitName);
int GiD_ResultUserDefined_HDF5(GP_CONST char* Name,GP_CONST char* Value);
int GiD_ResultValues_HDF5();
int GiD_EndResult_HDF5();
int GiD_FlushPostFile_HDF5();


int GiD_WriteScalar_HDF5(int id, double v);
int GiD_Write2DVector_HDF5(int id, double x, double y);
int GiD_WriteVector_HDF5(int id, double x, double y, double z);
int GiD_WriteVectorModule_HDF5(int id, double x, double y, double z, double mod);
int GiD_Write2DMatrix_HDF5(int id, double Sxx, double Syy, double Sxy);
int GiD_Write3DMatrix_HDF5(int id, double Sxx, double Syy, double Szz,
  double Sxy, double Syz, double Sxz);
int GiD_WritePlainDefMatrix_HDF5(int id, double Sxx, double Syy, double Sxy,
  double Szz);
int GiD_WriteMainMatrix_HDF5(int id,
  double Si, double Sii, double Siii,
  double Vix, double Viy, double Viz,
  double Viix, double Viiy, double Viiz,
  double Viiix, double Viiiy, double Viiiz);
int GiD_WriteLocalAxes_HDF5(int id, double euler_1, double euler_2, double euler_3);
/* complex results */
int GiD_WriteComplexScalar_HDF5(int id, double complex_real, double complex_imag);
int GiD_Write2DComplexVector_HDF5(int id, 
				  double x_real, double x_imag,
				  double y_real, double y_imag);
int GiD_WriteComplexVector_HDF5(int id, 
				double x_real, double x_imag,
				double y_real, double y_imag,
				double z_real, double z_imag);

#endif



