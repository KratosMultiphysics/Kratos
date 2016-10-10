
/*###############################################################################
 *    command to review contents of the file:
 *           h5dump --string FILE.flavia.res > FILE.txt
 *
 *###############################################################################*/

#include "gidpost.h"

#ifdef HDF5

#include <stdio.h>
#include <string.h>

#include "hdf5c.h"

static int current_mesh_num;
static int current_mesh_dataset_id;
static int current_mesh_nnode;
static GiD_ElementType current_mesh_etype;
static int current_gauss_points_num;
static int current_gauss_points_internal_coord;
static int current_range_table_num;
static int current_range_table_idx_num;

struct Myresult
{
  int num;
  int dataset_id;
  char name[2048];
};

struct MyresultGroup
{
  char Analysis[2048];
  double step;
  GiD_ResultLocation Where;
  char GaussPointsName[2048];
};

static int num_results_total,curr_result_group,num_results_group;
static struct Myresult myresults[MAX_CONCURRENT_DATASETS]; 
static struct MyresultGroup curr_my_result_group;

int GiD_OpenPostMeshFile_HDF5(GP_CONST char * FileName)
{
  int ret;
  current_mesh_num=0;
  ret=hdf5c_init(FileName); 
  if(ret>=0){
    hdf5c_set_attribute("/","GiD Post Results File","1.1");
    return 0;
  }
  return ret;
}

int GiD_ClosePostMeshFile_HDF5()
{
  int ret;
  ret=hdf5c_end();
  if(ret<0) return ret;
  return 0;
}

int GiD_BeginMesh_HDF5(GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,int NNode)
{
  char meshN[1024],buf[1024];
  char* enames[]={"NoElement","Point","Linear","Triangle","Quadrilateral","Tetrahedra","Hexahedra","Prism","Pyramid",
      "Sphere","Circle","Point"};
  
  hdf5c_create_group("Meshes");
  current_mesh_num++;
  sprintf(meshN,"Meshes/%d",current_mesh_num);
  hdf5c_create_group(meshN);
  if(MeshName) hdf5c_set_attribute(meshN,"Name",MeshName);

  switch(Dim){
    case GiD_2D: strcpy(buf,"2"); break;
    case GiD_3D: strcpy(buf,"3"); break;
  }
  hdf5c_set_attribute(meshN,"Dimension",buf);
  
  hdf5c_set_attribute(meshN,"ElemType",enames[EType]);
  sprintf(buf,"%d",NNode);
  hdf5c_set_attribute(meshN,"Nnode",buf);
  current_mesh_nnode=NNode;
  current_mesh_etype=EType;
  return 0;
}

int GiD_BeginMeshColor_HDF5(GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,
  int NNode,double Red, double Green, double Blue)
{
  int ret;
  char meshN[1024],buf[1024];
  ret=GiD_BeginMesh_HDF5(MeshName,Dim,EType,NNode);
  sprintf(meshN,"Meshes/%d",current_mesh_num);
  sprintf(buf,"%g %g %g",Red,Green,Blue);
  hdf5c_set_attribute(meshN,"Color",buf);
  return ret;
}

int GiD_EndMesh_HDF5()
{
  return 0;
}

int GiD_MeshUnit_HDF5(char * UnitName)
{
  char meshN[1024];
  if(UnitName){
    sprintf(meshN,"Meshes/%d",current_mesh_num);
     hdf5c_set_attribute(meshN,"UnitName",UnitName); 
  }
  return 0;
}


int GiD_MeshLocalAxes_HDF5(GP_CONST char * Result, GP_CONST char * Analysis,double step)
{
  char meshN[1024],buf[1024];
  sprintf(meshN,"Meshes/%d",current_mesh_num);
  if(Analysis){
    hdf5c_set_attribute(meshN,"LocalAxes Analysis",Analysis);
  }
  hdf5c_set_attribute(meshN,"LocalAxes Result",Result);
  sprintf(buf,"%g",step);
  hdf5c_set_attribute(meshN,"LocalAxes step",buf);
  return 0;
}

int GiD_BeginCoordinates_HDF5()
{
  char setN[1024];
  sprintf(setN,"Meshes/%d/Coordinates",current_mesh_num);
  current_mesh_dataset_id=hdf5c_start_dataset(setN,1,3);
  return current_mesh_dataset_id;
}

int GiD_EndCoordinates_HDF5()
{
  return hdf5c_end_dataset(current_mesh_dataset_id);
}

int GiD_WriteCoordinates_HDF5(int id, double x, double y, double z)
{
  int intvalues[1];
  double doublevalues[3];
  intvalues[0]=id;
  doublevalues[0]=x; doublevalues[1]=y; doublevalues[2]=z;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteCoordinates2D_HDF5(int id, double x, double y)
{
  return GiD_WriteCoordinates(id,x,y,0.0);
}

int GiD_BeginElements_HDF5()
  {
  char setN[1024];
  int num_int,num_real=0;
  
  switch(current_mesh_etype){
    case GiD_Sphere: num_int=3; num_real=1; break;
    case GiD_Circle: num_int=3; num_real=4; break;
    case GiD_Cluster: num_int=3; break;
    default: num_int=current_mesh_nnode+2; break;
  }
  sprintf(setN,"Meshes/%d/Elements",current_mesh_num);
  current_mesh_dataset_id=hdf5c_start_dataset(setN,num_int,num_real);
  return current_mesh_dataset_id;
}

int GiD_EndElements_HDF5()
{
  return hdf5c_end_dataset(current_mesh_dataset_id);
}

int GiD_WriteElement_HDF5(int id, int nid[])
{
  int intvalues[30];
  double doublevalues[1];
  intvalues[0]=id;
  memcpy(&intvalues[1],&nid[0],current_mesh_nnode*sizeof(int));
  intvalues[current_mesh_nnode+1]=0;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteElementMat_HDF5(int id, int nid[])
{
  int intvalues[30];
  double doublevalues[1];
  intvalues[0]=id;
  memcpy(&intvalues[1],&nid[0],(current_mesh_nnode+1)*sizeof(int));
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}
 
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

int GiD_WriteSphere_HDF5(int id, int nid, double r)
{
  int intvalues[3];
  double doublevalues[1];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=0;
  doublevalues[0]=r;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteSphereMat_HDF5(int id, int nid, double r, int mat)
{
  int intvalues[3];
  double doublevalues[1];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=mat;
  doublevalues[0]=r;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteCircle_HDF5(int id, int nid, double r,double nx, double ny, double nz)
{
  int intvalues[3];
  double doublevalues[4];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=0;
  doublevalues[0]=r;
  doublevalues[1]=nx;
  doublevalues[2]=ny;
  doublevalues[3]=nz;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteCircleMat_HDF5(int id, int nid, double r,double nx, double ny, double nz, int mat)
{
  int intvalues[3];
  double doublevalues[4];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=mat;
  doublevalues[0]=r;
  doublevalues[1]=nx;
  doublevalues[2]=ny;
  doublevalues[3]=nz;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

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

int GiD_WriteCluster_HDF5(int id, int nid)
{
  int intvalues[3];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=0;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues);
  return 0;
}

int GiD_WriteClusterMat_HDF5(int id, int nid, int mat)
{
  int intvalues[3];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=mat;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues);
  return 0;
}

int GiD_OpenPostResultFile_HDF5(GP_CONST char * FileName)
{
  int ret;
  current_mesh_num=0;
  current_gauss_points_num=0;
  current_range_table_num=0;
  num_results_total=0;
  curr_result_group=0;
  num_results_group=0;
  ret=hdf5c_init(FileName); 
  if(ret>=0){
    hdf5c_set_attribute("/","GiD Post Results File","1.1");
    return 0;
  }
  return ret;
}

int GiD_ClosePostResultFile_HDF5()
{
  int ret;
  ret=hdf5c_end();
  if(ret<0) return ret;
  return 0;
}

int GiD_BeginGaussPoint_HDF5(GP_CONST char * name, GiD_ElementType EType,GP_CONST char * MeshName,
  int GP_number, int NodesIncluded, int InternalCoord)
{
  char gpN[1024],buf[1024];
  char* enames[]={"NoElement","Point","Linear","Triangle","Quadrilateral","Tetrahedra","Hexahedra","Prism","Pyramid",
      "Sphere","Circle","Point"};

  hdf5c_create_group("GaussPoints");
  current_gauss_points_num++;
  sprintf(gpN,"GaussPoints/%d",current_gauss_points_num);
  hdf5c_create_group(gpN);
  hdf5c_set_attribute(gpN,"Name",name);
  hdf5c_set_attribute(gpN,"ElemType",enames[EType]);
  if(MeshName) hdf5c_set_attribute(gpN,"MeshName",MeshName);
  sprintf(buf,"%d",GP_number);
  hdf5c_set_attribute(gpN,"GP_number",buf);
  sprintf(buf,"%d",NodesIncluded);
  hdf5c_set_attribute(gpN,"NodesIncluded",buf);
  sprintf(buf,"%d",InternalCoord);
  hdf5c_set_attribute(gpN,"InternalCoord",buf);

  if(InternalCoord==0){
    sprintf(gpN,"GaussPoints/%d",current_gauss_points_num);
    current_gauss_points_internal_coord=0;
    current_mesh_dataset_id=hdf5c_start_dataset(gpN,0,3);
    return current_mesh_dataset_id;
  }
  else {
    current_gauss_points_internal_coord=1;
    return 0;
  }
}

int GiD_EndGaussPoint_HDF5()
{
  if(current_gauss_points_internal_coord==0)
    return hdf5c_end_dataset(current_mesh_dataset_id);
  else
    return 0;
}

int GiD_WriteGaussPoint2D_HDF5(double x, double y)
{
  int intvalues[1];
  double doublevalues[3];
  doublevalues[0]=x; doublevalues[1]=y; doublevalues[2]=0.0;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteGaussPoint3D_HDF5(double x, double y, double z)
{
  int intvalues[1];
  double doublevalues[3];
  doublevalues[0]=x; doublevalues[1]=y; doublevalues[2]=z;
  hdf5c_addto_dataset(current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_BeginRangeTable_HDF5(GP_CONST char * name)
{
  char rtN[1024];  
  hdf5c_create_group("ResultRangesTable");
  current_range_table_num++;
  current_range_table_idx_num=0;
  sprintf(rtN,"ResultRangesTable/%d",current_range_table_num);
  hdf5c_create_group(rtN);
  hdf5c_set_attribute(rtN,"Name",name);
  return 0;
}

int GiD_EndRangeTable_HDF5()
{
  return 0;
}
 
int GiD_WriteMinRange_HDF5(double max, GP_CONST char * name)
{
  char rtN[1024],buf[1024];
  current_range_table_idx_num++;
  sprintf(rtN,"ResultRangesTable/%d/%d",current_range_table_num,current_range_table_idx_num);
  hdf5c_create_group(rtN);
  hdf5c_set_attribute(rtN,"Name",name);
  sprintf(buf,"%g",max);
  hdf5c_set_attribute(rtN,"Max",buf);
  return 0;
}

int GiD_WriteRange_HDF5(double min, double max, GP_CONST char * name)
{
  char rtN[1024],buf[1024];
  current_range_table_idx_num++;
  sprintf(rtN,"ResultRangesTable/%d/%d",current_range_table_num,current_range_table_idx_num);
  hdf5c_create_group(rtN);
  hdf5c_set_attribute(rtN,"Name",name);
  sprintf(buf,"%g",min);
  hdf5c_set_attribute(rtN,"Min",buf);
  sprintf(buf,"%g",max);
  hdf5c_set_attribute(rtN,"Max",buf);
  return 0;
}

int GiD_WriteMaxRange_HDF5(double min, GP_CONST char * name)
{
  char rtN[1024],buf[1024];
  current_range_table_idx_num++;
  sprintf(rtN,"ResultRangesTable/%d/%d",current_range_table_num,current_range_table_idx_num);
  hdf5c_create_group(rtN);
  hdf5c_set_attribute(rtN,"Name",name);
  sprintf(buf,"%g",min);
  hdf5c_set_attribute(rtN,"Min",buf);
  return 0;
}

typedef enum {
  RH_start,
  RH_add
} ResultHeaderMode;

int _GiD_BeginResultHeader_HDF5_init(GP_CONST char * Result, GP_CONST char * Analysis, double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName,ResultHeaderMode RHmode)
{  
  char resN[2048],buf[2048];
  char* rtnames[]={"Scalar","Vector","Matrix","PlainDeformationMatrix","MainMatrix","LocalAxes", 
		   "ComplexScalar", "ComplexVector"};
  
  if(RHmode==RH_start){
    num_results_group=1;
  } else {
    num_results_group++;
  }
  curr_result_group=0;

  myresults[num_results_group-1].num=++num_results_total;
  myresults[num_results_group-1].dataset_id=-1;
  sprintf(resN,"Results/%d",num_results_total);
  strcpy(myresults[num_results_group-1].name,resN);

  hdf5c_create_group("Results");
  hdf5c_create_group(resN);

  hdf5c_set_attribute(resN,"Name",Result);
  hdf5c_set_attribute(resN,"Analysis",Analysis);

  sprintf(buf,"%g",step);
  hdf5c_set_attribute(resN,"Step",buf);
  
  hdf5c_set_attribute(resN,"ResultType",rtnames[Type]);
  switch(Where){
    case GiD_OnNodes: hdf5c_set_attribute(resN,"ResultLocation","OnNodes"); break;
    case GiD_OnGaussPoints: hdf5c_set_attribute(resN,"ResultLocation","OnGaussPoints"); break;
  }
  if(GaussPointsName) hdf5c_set_attribute(resN,"GaussPointsName",GaussPointsName); 
  return 0;
}

int GiD_BeginResult_HDF5(GP_CONST char * Result, GP_CONST char * Analysis,double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName,
  GP_CONST char * RangeTable, 
  int compc, GP_CONST char * compv[])
{
  int fail,i;
  char* resN,buf[2048];
  
  fail=_GiD_BeginResultHeader_HDF5_init(Result,Analysis,step,Type,Where,GaussPointsName,RH_start);
  if(fail<0) return 1;
  resN=myresults[num_results_group-1].name;

  if(RangeTable) hdf5c_set_attribute(resN,"RangeTable",RangeTable); 
  if(compc){
    sprintf(buf,"%d",compc);
    hdf5c_set_attribute(resN,"NumComponents",buf); 
  }
  for(i=0;i<compc;i++){
    sprintf(buf,"Component %d",i+1);
    hdf5c_set_attribute(resN,buf,compv[i]); 
  }
  return 0;
}


int GiD_BeginResultHeader_HDF5(GP_CONST char * Result, GP_CONST char * Analysis, double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName)
{
  int fail;
  fail=_GiD_BeginResultHeader_HDF5_init(Result,Analysis,step,Type,Where,GaussPointsName,RH_start);
  if(fail<0) return 1;
  return 0;
}

int GiD_BeginResultGroup_HDF5(GP_CONST char * Analysis, double step,GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName)
{
  strcpy(curr_my_result_group.Analysis,Analysis);
  curr_my_result_group.step=step;
  curr_my_result_group.Where=Where;
  strcpy(curr_my_result_group.GaussPointsName,GaussPointsName);
  
  num_results_group=0;
  return 0;
}

int GiD_ResultDescription_HDF5(GP_CONST char * Result, GiD_ResultType Type)
{
  int fail;
  fail=_GiD_BeginResultHeader_HDF5_init(Result,curr_my_result_group.Analysis,curr_my_result_group.step,
    Type,curr_my_result_group.Where,curr_my_result_group.GaussPointsName,RH_add);
  if(fail<0) return 1;
  return 0;
}

/*
    Associates the given local axes result to the current result.
    If result is a vector or a matrix, it means that the result is expressed in
    these local axes. Vector vx,vy,vz is not used
    If result is a scalar, local axes are a help for drawing the LineDiagrams
    the direction of the line diagram is given by vx,vy,vz expressed in Local Axes
    Analysis can be NULL and means the same analysis
*/
int GiD_ResultLocalAxes_HDF5(GP_CONST char * Result, GP_CONST char * Analysis,
		                  double step,double vx,double vy,double vz)
{
  char bufres[2048],*resN;
  resN=myresults[num_results_group-1].name;
  if(Analysis){
    hdf5c_set_attribute(resN,"LocalAxes Analysis",Analysis);
  }
  hdf5c_set_attribute(resN,"LocalAxes Result",Result);
  sprintf(bufres,"%g",step);
  hdf5c_set_attribute(resN,"LocalAxes step",bufres);
  sprintf(bufres,"%g,%g,%g",vx,vy,vz);
  hdf5c_set_attribute(resN,"LocalAxes vector",bufres);
  return 0;
}

int GiD_ResultIsDeformationVector_HDF5(int boolean)
{
  char* resN=myresults[num_results_group-1].name;

  if(boolean) hdf5c_set_attribute(resN,"IsDeformationVector","1");
  else hdf5c_set_attribute(resN,"IsDeformationVector","0");
  return 0;
}

int GiD_ResultRange_HDF5(GP_CONST char * RangeTable)
{
  char* resN=myresults[num_results_group-1].name;
  if(RangeTable) hdf5c_set_attribute(resN,"RangeTable",RangeTable); 
  return 0;
}

int GiD_ResultComponents_HDF5(int compc, GP_CONST char * compv[])
{
  int i;
  char buf[2048];
  char* resN=myresults[num_results_group-1].name;
  if(compc){
    sprintf(buf,"%d",compc);
    hdf5c_set_attribute(resN,"NumComponents",buf); 
  }
  for(i=0;i<compc;i++){
    sprintf(buf,"Component %d",i+1);
    hdf5c_set_attribute(resN,buf,compv[i]); 
  }
  return 0;
}

int GiD_ResultUnit_HDF5(GP_CONST char * UnitName)
{
  char* resN=myresults[num_results_group-1].name;
  if(UnitName) hdf5c_set_attribute(resN,"UnitName",UnitName); 
  return 0;
}

int GiD_ResultUserDefined_HDF5(GP_CONST char* Name,GP_CONST char* Value)
{
  char* resN=myresults[num_results_group-1].name;
  hdf5c_set_attribute(resN,Name,Value);
  return 0;
}

int GiD_ResultValues_HDF5()
{
  return 0;
}

int GiD_EndResult_HDF5()
{
  int i,fail;
  for(i=0;i<num_results_group;i++){
    fail=hdf5c_end_dataset(myresults[i].dataset_id);
    if(fail<0) return fail;
  }
  num_results_group=0;
  return 0;
}

int GiD_FlushPostFile_HDF5()
{
  int ret;
  ret=hdf5c_flush();
  if(ret<0) return ret;
  return 0;
}

int GiD_WriteScalar_HDF5(int id, double v)
{
  int intvalues[1];
  double doublevalues[1];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,1);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=v;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_Write2DVector_HDF5(int id, double x, double y)
{
  int intvalues[1];
  double doublevalues[2];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,2);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x;  doublevalues[1]=y;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_WriteVector_HDF5(int id, double x, double y, double z)
{
  int intvalues[1];
  double doublevalues[3];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,3);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x;  doublevalues[1]=y; doublevalues[2]=z;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_WriteVectorModule_HDF5(int id, double x, double y, double z, double mod)
{
  int intvalues[1];
  double doublevalues[4];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,4);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x;  doublevalues[1]=y; doublevalues[2]=z; doublevalues[3]=mod;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_Write2DMatrix_HDF5(int id, double Sxx, double Syy, double Sxy)
{
  int intvalues[1];
  double doublevalues[3];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,3);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Sxx;  doublevalues[1]=Syy; doublevalues[2]=Sxy;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_Write3DMatrix_HDF5(int id, double Sxx, double Syy, double Szz,
  double Sxy, double Syz, double Sxz)
{
  int intvalues[1];
  double doublevalues[6];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,6);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Sxx;  doublevalues[1]=Syy; doublevalues[2]=Szz;
  doublevalues[3]=Sxy;  doublevalues[4]=Syz; doublevalues[5]=Sxz;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}
 
int GiD_WritePlainDefMatrix_HDF5(int id, double Sxx, double Syy, double Sxy,
  double Szz)
{
  int intvalues[1];
  double doublevalues[4];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,4);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Sxx;  doublevalues[1]=Syy; doublevalues[2]=Sxy;
  doublevalues[3]=Szz;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_WriteMainMatrix_HDF5(int id,
  double Si, double Sii, double Siii,
  double Vix, double Viy, double Viz,
  double Viix, double Viiy, double Viiz,
  double Viiix, double Viiiy, double Viiiz)
{
  int intvalues[1];
  double doublevalues[12];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,12);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Si;  doublevalues[1]=Sii; doublevalues[2]=Siii;
  doublevalues[3]=Vix;  doublevalues[4]=Viy; doublevalues[5]=Viz;
  doublevalues[6]=Viix;  doublevalues[7]=Viiy; doublevalues[8]=Viiz;
  doublevalues[9]=Viiix;  doublevalues[10]=Viiiy; doublevalues[11]=Viiiz;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_WriteLocalAxes_HDF5(int id, double euler_1, double euler_2, double euler_3)
{
  int intvalues[1];
  double doublevalues[3];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,3);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=euler_1;  doublevalues[1]=euler_2; doublevalues[2]=euler_3;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_WriteComplexScalar_HDF5( int id, double complex_real, double complex_imag) {
  int intvalues[1];
  double doublevalues[2];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,2);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=complex_real;
  doublevalues[1]=complex_imag;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_Write2DComplexVector_HDF5( int id,
				   double x_real, double x_imag,
				   double y_real, double y_imag) {
  int intvalues[1];
  double doublevalues[4];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,4);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x_real;
  doublevalues[1]=x_imag;
  doublevalues[2]=y_real;
  doublevalues[3]=y_imag;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

int GiD_WriteComplexVector_HDF5( int id,
				 double x_real, double x_imag,
				 double y_real, double y_imag,
				 double z_real, double z_imag) {
  int intvalues[1];
  double doublevalues[6];
  char* resN=myresults[curr_result_group].name;
  if(myresults[curr_result_group].dataset_id==-1){
    myresults[curr_result_group].dataset_id=hdf5c_start_dataset(resN,1,6);
    if(myresults[curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x_real;
  doublevalues[1]=x_imag;
  doublevalues[2]=y_real;
  doublevalues[3]=y_imag;
  doublevalues[4]=z_real;
  doublevalues[5]=z_imag;
  hdf5c_addto_dataset(myresults[curr_result_group].dataset_id,intvalues,doublevalues);
  curr_result_group++;
  if(curr_result_group==num_results_group) curr_result_group=0;
  return 0;
}

#endif
