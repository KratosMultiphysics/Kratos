
/*###############################################################################
 *    command to review contents of the file:
 *           h5dump --string FILE.flavia.res > FILE.txt
 *
 *###############################################################################*/

#include <assert.h>

#include "gidpost.h"
#include "gidpostHDF5.h"
#include "gidpostMutex.h"

// #ifdef ENABLE_HDF5

#include <stdio.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#define strdup( S)    _strdup( S)
#endif
#include "hdf5c.h"


/* private structures */

struct Myresult {
  int num;
  int dataset_id;
  char name[ 2048 ];
};

struct MyresultGroup {
  char Analysis[ 2048 ];
  double step;
  GiD_ResultLocation Where;
  char GaussPointsName[ 2048 ];
};

typedef struct _MeshGroupData {
  int current_mesh_group_id;
  int using_mesh_group;
  int current_on_mesh_group_id;
  int using_on_mesh_group;
  int num_mesh_groups;
  char **lst_mesh_groups_names;
  char current_mesh_group_name[ 2048 ];
  char current_on_mesh_group_name[ 2048 ];
} MeshGroupData;

static MeshGroupData *new_MeshGroupData( void ) {
  MeshGroupData *obj = ( MeshGroupData * )malloc( sizeof( MeshGroupData ) );
  if ( obj ) {
    memset( obj, 0, sizeof( MeshGroupData ) );
  }
  return obj;
}

static void delete_MeshGroupData( MeshGroupData *obj ) {
  free( obj );
}

static int resetMeshGroupData( MeshGroupData *obj ) {
  if ( obj == NULL )
    return -1;
  obj->current_mesh_group_id = 0;
  obj->using_mesh_group = 0;
  obj->current_on_mesh_group_id = 0;
  obj->using_on_mesh_group = 0;
  int i = 0;
  if ( obj->lst_mesh_groups_names ) {
    for ( i = 0; i < obj->num_mesh_groups; i++ ) {
      free( obj->lst_mesh_groups_names[ i ] );
    }
  }
  free( obj->lst_mesh_groups_names );
  obj->lst_mesh_groups_names = NULL;
  obj->num_mesh_groups = 0;
  return 0;
}

/* public structure */

#define LOCAL_AXES_MAX_LEN    500

typedef struct _CurrentHdf5WriteData {
  int current_mesh_num;
  int current_mesh_dataset_id;
  int current_mesh_nnode;
  GiD_ElementType current_mesh_etype;
  int current_gauss_points_num;
  int current_gauss_points_internal_coord;
  int current_range_table_num;
  int current_range_table_idx_num;
  int num_results_total;
  int curr_result_group;
  int num_results_group;
  struct Myresult myresults[ MAX_CONCURRENT_DATASETS ];
  struct MyresultGroup curr_my_result_group;
  char local_axes_format[ LOCAL_AXES_MAX_LEN ];
  int local_axes_format_created;
  CPostHdf5File *post_h5_file;
  MeshGroupData *mesh_group_data;
} CurrentHdf5WriteData;

CurrentHdf5WriteData *new_CurrentHdf5WriteData( void ) {
  static int G_num_HDF5_files = 0;

  /* if we are opening more than one HDF5 file, check for thread safety */
  _LOCK_;
  if ( G_num_HDF5_files == 1 ) {
    if ( GiD_IsThreadSafe_HDF5() <= 0 ) {
      fprintf( stderr, "GiDPost: HDF5 library is not thread safe. Using %s\n", GiD_GetHDF5Version() );
      return NULL;
    }
  }
  G_num_HDF5_files++;
  _UNLOCK_;

  CurrentHdf5WriteData *obj = ( CurrentHdf5WriteData * )malloc( sizeof( CurrentHdf5WriteData ) );
  if ( obj ) {
    memset( obj, 0, sizeof( CurrentHdf5WriteData ) );
    obj->local_axes_format_created = 0;
  }
  return obj;
}

void delete_CurrentHdf5WriteData( CurrentHdf5WriteData *obj ) {
  delete_MeshGroupData( obj->mesh_group_data );
  free( obj );
}

/* private functions */

static const char *getCurrentMeshGroupPath( MeshGroupData *obj, char *dst_current_mesh_group_path, size_t max_len ) {
  const char *ret = NULL;
  if ( obj->using_mesh_group && ( obj->current_mesh_group_id >= 0 ) ) {
    snprintf( dst_current_mesh_group_path, max_len, "MeshGroup%d", obj->current_mesh_group_id + 1 );
    ret = dst_current_mesh_group_path;
  }
  return ret;
}

static const char *getCurrentMeshPath( CurrentHdf5WriteData *obj, char *dst_mesh_path, size_t max_len ) {
  MeshGroupData *mgd = obj->mesh_group_data;
  char current_mesh_group_path[ 900 ];
  const char *mesh_group = getCurrentMeshGroupPath( mgd, current_mesh_group_path, 900 );
  if ( mesh_group ) {
    hdf5c_f_create_group( obj->post_h5_file, mesh_group );
    hdf5c_f_set_attribute( obj->post_h5_file, mesh_group, "Name", mgd->current_mesh_group_name );
    snprintf( dst_mesh_path, max_len, "%s/Meshes", mesh_group );
  } else {
    strcpy( dst_mesh_path, "Meshes" );
  }
  return dst_mesh_path;
}

static const char *getCurrentOnMeshGroupPath( MeshGroupData *obj, char *dst_current_on_mesh_group_path, size_t max_len ) {
  const char *ret = NULL;
  if ( obj->using_on_mesh_group && ( obj->current_on_mesh_group_id >= 0 ) ) {
    snprintf( dst_current_on_mesh_group_path, max_len, "MeshGroup%d", obj->current_on_mesh_group_id + 1 );
    ret = dst_current_on_mesh_group_path;
  }
  return ret;
}

static const char *getCurrentResultPath( CurrentHdf5WriteData *obj, char *result_path, size_t max_len ) {
  MeshGroupData *mgd = obj->mesh_group_data;
  char dst_mesh_group_path[ 900 ];
  const char *result_group = getCurrentOnMeshGroupPath( mgd, dst_mesh_group_path, 900 );
  if ( result_group ) {
    hdf5c_f_create_group( obj->post_h5_file, result_group );
    hdf5c_f_set_attribute( obj->post_h5_file, result_group, "Name", mgd->current_mesh_group_name );
    snprintf( result_path, max_len, "%s/Results", result_group );
  } else {
    strcpy( result_path, "Results" );
  }
  return result_path;
}

static int initGlobalHdf5File( CurrentHdf5WriteData *h5_wd) {
  if ( h5_wd == NULL )
    return -1;
  if ( h5_wd->post_h5_file != NULL )
    return -1;
  h5_wd->post_h5_file = new_CPostHdf5File();
  return h5_wd->post_h5_file ? 0 : -1;
}

static int endGlobalHdf5File( CurrentHdf5WriteData *h5_wd ) {
  if ( h5_wd == NULL )
    return -1;
  if ( h5_wd->post_h5_file == NULL )
    return -1;
  delete_CPostHdf5File( h5_wd->post_h5_file );
  h5_wd->post_h5_file = NULL;
  return 0;
}


/* public functions */
const char *GiD_GetHDF5Version( void ) {
  return hdf5c_f_get_version();
}

int GiD_IsThreadSafe_HDF5( void ) {
  return hdf5c_f_is_thread_safe();
}

int GiD_OpenPostMeshFile_HDF5( CurrentHdf5WriteData *obj,  GP_CONST char *FileName ) {
  if ( !obj ) {
    return -1;
  }
  int ret;
  ret = initGlobalHdf5File( obj );
  if ( ret < 0 )
    return ret;
  obj->current_mesh_num = 0;
  obj->mesh_group_data = new_MeshGroupData();
  if ( !resetMeshGroupData( obj->mesh_group_data ) )
    return -1;
  ret = hdf5c_f_init( obj->post_h5_file, FileName );
  if ( ret >= 0 ) {
    hdf5c_f_set_attribute( obj->post_h5_file, "/", "GiD Post Results File", "1.1" );
    return 0;
  }
  return ret;
}

int GiD_ClosePostMeshFile_HDF5( CurrentHdf5WriteData *obj ) {
  if ( !obj ) {
    return -1;
  }
  int ret;
  ret=hdf5c_f_end( obj->post_h5_file);
  if(ret<0) 
    return ret;
  ret = endGlobalHdf5File( obj );
  if ( ret < 0 )
    return ret;
  return 0;
}

int GiD_BeginMesh_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,int NNode)
{
  char meshN[1024],buf[1024];
  char* enames[]={"NoElement","Point","Linear","Triangle","Quadrilateral","Tetrahedra","Hexahedra","Prism","Pyramid",
      "Sphere","Circle"};

  char dst_mesh_path[ 1024 ];
  const char *mesh_path = getCurrentMeshPath( obj, dst_mesh_path, 1024 );
  hdf5c_f_create_group( obj->post_h5_file, mesh_path);
  obj->current_mesh_num++;
  snprintf(meshN, 1024, "%s/%d", mesh_path, obj->current_mesh_num);
  hdf5c_f_create_group( obj->post_h5_file, meshN);
  if(MeshName) 
    hdf5c_f_set_attribute( obj->post_h5_file, meshN,"Name",MeshName);

  switch(Dim){
    case GiD_2D: strcpy(buf,"2"); break;
    case GiD_3D: strcpy(buf,"3"); break;
  }
  hdf5c_f_set_attribute( obj->post_h5_file, meshN,"Dimension",buf);
  
  hdf5c_f_set_attribute( obj->post_h5_file, meshN,"ElemType",enames[EType]);
  snprintf(buf, 1024, "%d",NNode);
  hdf5c_f_set_attribute( obj->post_h5_file, meshN,"Nnode",buf);
  obj->current_mesh_nnode=NNode;
  obj->current_mesh_etype=EType;
  return 0;
}

int GiD_BeginMeshColor_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * MeshName,GiD_Dimension Dim, GiD_ElementType EType,
  int NNode,double Red, double Green, double Blue)
{
  int ret;
  char meshN[1024],buf[1024];
  ret=GiD_BeginMesh_HDF5( obj, MeshName,Dim,EType,NNode);
  char dst_mesh_path[ 1024 ];
  const char *mesh_path = getCurrentMeshPath( obj, dst_mesh_path, 1024 );
  snprintf(meshN,1024, "%s/%d", mesh_path, obj->current_mesh_num);
  snprintf(buf,1024, "%f %f %f",Red,Green,Blue);
  hdf5c_f_set_attribute( obj->post_h5_file, meshN,"Color",buf);
  return ret;
}

int GiD_EndMesh_HDF5( CurrentHdf5WriteData *obj )
{
  return 0;
}

int GiD_MeshUnit_HDF5( CurrentHdf5WriteData *obj,  GP_CONST char * UnitName)
{
  if ( UnitName ) {
    char meshN[ 1024 ];
    char dst_mesh_path[ 1024 ];
    snprintf( meshN, 1024, "%s/%d", getCurrentMeshPath( obj, dst_mesh_path, 1024 ), obj->current_mesh_num );
    hdf5c_f_set_attribute( obj->post_h5_file, meshN, "UnitName", UnitName );
  }
  return 0;
}


int GiD_MeshLocalAxes_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis,double step)
{
  char meshN[1024],buf[1024];
  char dst_mesh_path[ 1024 ];
  snprintf(meshN, 1024, "%s/%d", getCurrentMeshPath( obj, dst_mesh_path, 1024 ), obj->current_mesh_num);
  if(Analysis){
    hdf5c_f_set_attribute( obj->post_h5_file, meshN,"LocalAxes Analysis",Analysis);
  }
  hdf5c_f_set_attribute( obj->post_h5_file, meshN,"LocalAxes Result",Result);
  snprintf(buf, 1024, GiD_PostGetFormatStep(),step);
  hdf5c_f_set_attribute( obj->post_h5_file, meshN,"LocalAxes step",buf);
  return 0;
}

int GiD_BeginCoordinates_HDF5( CurrentHdf5WriteData *obj )
{
  char setN[1024];
  char dst_mesh_path[ 1024 ];
  snprintf(setN, 1024, "%s/%d/Coordinates", getCurrentMeshPath( obj, dst_mesh_path, 1024 ), obj->current_mesh_num);
  obj->current_mesh_dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, setN,1,3);
  return obj->current_mesh_dataset_id;
}

int GiD_EndCoordinates_HDF5( CurrentHdf5WriteData *obj )
{
  return hdf5c_f_end_dataset( obj->post_h5_file, obj->current_mesh_dataset_id);
}

int GiD_WriteCoordinates_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y, double z)
{
  int intvalues[1];
  double doublevalues[3];
  intvalues[0]=id;
  doublevalues[0]=x; doublevalues[1]=y; doublevalues[2]=z;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteCoordinates2D_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y)
{
  return GiD_WriteCoordinates_HDF5( obj, id, x, y, 0.0);
}


// WriteCoordinatesBlock includes BeginCoordinates() and EndCoordinates()
int GiD_WriteCoordinatesBlock_HDF5( CurrentHdf5WriteData *obj, int num_points,
                                    const double *xyz_array ) {
  int ret_val = 0;
  if ( num_points > 0 ) {
    int *list_ids = ( int * )malloc( ( size_t)num_points * sizeof( int ) );
    if ( list_ids ) {
      for ( int idx_point = 0; idx_point < num_points; idx_point++ ) {
        list_ids[ idx_point ] = idx_point + 1;
      }
      ret_val = GiD_WriteCoordinatesIdBlock_HDF5( obj, num_points, list_ids, xyz_array );
      free( list_ids );
    }
  }
  return ret_val;
}

// WriteCoordinatesBlock includes BeginCoordinates() and EndCoordinates()
int GiD_WriteCoordinatesIdBlock_HDF5( CurrentHdf5WriteData *obj, int num_points, 
                                      const int *list_ids, const double *xyz_array ) {
  // in GiD Post hdf5 ( = compassis hdf5) data is stored column-wise in different datasets
  // for instance: 1 data set with 1 column for id's
  // 1 data set with 1 column for x-coordinate
  // 2 data set with 1 column for x-coordinate
  // 3 data set with 1 column for x-coordinate
  double *list_x = ( double * )malloc( ( size_t)num_points * sizeof( double ) );
  double *list_y = ( double * )malloc( ( size_t)num_points * sizeof( double ) );
  double *list_z = ( double * )malloc( ( size_t)num_points * sizeof( double ) );
  if ( list_x && list_y && list_z ) {
    for ( int idx_point = 0; idx_point < num_points; idx_point++ ) {
      list_x[ idx_point ] = xyz_array[ idx_point * 3 + 0 ];
      list_y[ idx_point ] = xyz_array[ idx_point * 3 + 1 ];
      list_z[ idx_point ] = xyz_array[ idx_point * 3 + 2 ];
    }

    char dataset_name[ 1024 ];
    char dst_mesh_path[ 1024 ];
    snprintf( dataset_name, 1024, "%s/%d/Coordinates", getCurrentMeshPath( obj, dst_mesh_path, 1024 ), obj->current_mesh_num );

    t_hdf5c_dataset_type lst_types[] = { HDF5C_DATASET_INTEGER, HDF5C_DATASET_DOUBLE, HDF5C_DATASET_DOUBLE, HDF5C_DATASET_DOUBLE };
    const void *lst_data[] = { list_ids, list_x, list_y, list_z};
    const int num_datasets = 4; // ids, x, y, z
    hdf5c_write_dataset_list( obj->post_h5_file, dataset_name, num_datasets, num_points, lst_types, lst_data );
    free( list_z ); list_z = NULL;
    free( list_y ); list_y = NULL;
    free( list_x ); list_z = NULL;
  }
  return 0;
}


int GiD_BeginElements_HDF5( CurrentHdf5WriteData *obj )
  {
  char setN[1024];
  int num_int,num_real=0;
  
  switch(obj->current_mesh_etype){
    case GiD_Sphere: num_int=3; num_real=1; break;
    case GiD_Circle: num_int=3; num_real=4; break;
    case GiD_Cluster: num_int=3; num_real=0; break;
    default: num_int=obj->current_mesh_nnode+2; break;
  }
  char dst_mesh_path[ 1024 ];
  snprintf(setN, 1024, "%s/%d/Elements", getCurrentMeshPath( obj, dst_mesh_path, 1024 ), obj->current_mesh_num);
  obj->current_mesh_dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, setN,num_int,num_real);
  return obj->current_mesh_dataset_id;
}

int GiD_EndElements_HDF5( CurrentHdf5WriteData *obj )
{
  return hdf5c_f_end_dataset( obj->post_h5_file, obj->current_mesh_dataset_id);
}

int GiD_WriteElement_HDF5( CurrentHdf5WriteData *obj, int id, int nid[] ) {
  int intvalues[ 30 ];
  double doublevalues[ 1 ] = { -1.0 };
  intvalues[ 0 ] = id;
  memcpy( &intvalues[ 1 ], &nid[ 0 ], ( size_t)( obj->current_mesh_nnode) * sizeof( int ) );
  intvalues[ obj->current_mesh_nnode + 1 ] = 0;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id, intvalues,
                          doublevalues );
  return 0;
}

int GiD_WriteElementMat_HDF5( CurrentHdf5WriteData *obj, int id, int nid[] ) {
  int intvalues[ 30 ];
  double doublevalues[ 1 ] = { -1.0 };
  intvalues[ 0 ] = id;
  memcpy( &intvalues[ 1 ], &nid[ 0 ], ( size_t)( obj->current_mesh_nnode + 1 ) * sizeof( int ) );
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id, intvalues,
                          doublevalues );
  return 0;
}

// GiD_WriteElementsBlock_HDF5 includes BeginElements and EndElements
int GiD_WriteElementsBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *connectivities ) {
  int ret_val = 0;
  if ( num_elements > 0 ) {
    int *list_ids = ( int * )malloc( ( size_t)num_elements * sizeof( int ) );
    int *list_mat = ( int * )malloc( ( size_t)num_elements * sizeof( int ) );
    if ( list_ids && list_mat ) {
      for ( int idx_elem = 0; idx_elem < num_elements; idx_elem++ ) {
        list_ids[ idx_elem ] = idx_elem + 1;
        list_mat[ idx_elem ] = 0;
      }
      ret_val = GiD_WriteElementsIdMatBlock_HDF5( obj, num_elements, list_ids, connectivities, list_mat );
      free( list_mat );
      free( list_ids );
    }
  }
  return ret_val;
}
int GiD_WriteElementsIdBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *list_ids, const int *connectivities ) {
  int ret_val = 0;
  if ( num_elements > 0 ) {
    int *list_mat = ( int * )malloc( ( size_t)num_elements * sizeof( int ) );
    if ( list_mat ) {
      for ( int idx_elem = 0; idx_elem < num_elements; idx_elem++ ) {
        list_mat[ idx_elem ] = 0;
      }
      ret_val = GiD_WriteElementsIdMatBlock_HDF5( obj, num_elements, list_ids, connectivities, list_mat );
      free( list_mat );
    }
  }
  return ret_val;
  return 0;
}
int GiD_WriteElementsMatBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *connectivities, const int *lst_material_id ) {
  int ret_val = 0;
  if ( num_elements > 0 ) {
    int *list_ids = ( int * )malloc( ( size_t)num_elements * sizeof( int ) );
    if ( list_ids ) {
      for ( int idx_elem = 0; idx_elem < num_elements; idx_elem++ ) {
        list_ids[ idx_elem ] = idx_elem + 1;
      }
      ret_val = GiD_WriteElementsIdMatBlock_HDF5( obj, num_elements, list_ids, connectivities, lst_material_id );
      free( list_ids );
    }
  }
  return ret_val;
}
int GiD_WriteElementsIdMatBlock_HDF5( CurrentHdf5WriteData *obj, int num_elements, const int *list_ids, const int *connectivities, const int *lst_material_id ) {
  // in GiD Post hdf5 ( = compassis hdf5) data is stored column-wise in different datasets
  // for instance: 1 data set with 1 column for id's
  // 1 data set with 1 column for x-coordinate
  // 2 data set with 1 column for x-coordinate
  // 3 data set with 1 column for x-coordinate
  const int num_nodes = obj->current_mesh_nnode;
  const int num_components = num_nodes + 1 + 1; // num_nodes per element + id + mat
  t_hdf5c_dataset_type *list_data_type = ( t_hdf5c_dataset_type * )malloc( ( size_t)num_components * sizeof( t_hdf5c_dataset_type ) );
  int **list_data_values = ( int ** )malloc( ( size_t)num_components * sizeof( int * ) );
  for ( int idx = 0; idx < num_components; idx++ ) {
    list_data_type[ idx ] = HDF5C_DATASET_INTEGER;
  }
  list_data_values[ 0 ] = ( int *)list_ids;
  int fail = 0;
  for ( int idx = 0; idx < num_nodes; idx++ ) {
    list_data_values[ idx + 1 ] = ( int * )malloc( ( size_t)num_elements * sizeof( int ) );
    if ( list_data_values[ idx + 1 ] == NULL ) {
      fail = -1;
      return fail;
    }
  }
  list_data_values[ num_nodes + 1 ] = ( int * )lst_material_id;
  for ( int idx_elem = 0; idx_elem < num_elements; idx_elem++ ) {
    for ( int idx = 0; idx < num_nodes; idx++ ) {
      list_data_values[ idx + 1 ][ idx_elem ] = connectivities[ idx_elem * num_nodes + idx ];
    }
  }

  char dataset_name[ 1024 ];
  char dst_mesh_path[ 1024 ];
  snprintf( dataset_name, 1024, "%s/%d/Elements", getCurrentMeshPath( obj, dst_mesh_path, 1024 ), obj->current_mesh_num );
  hdf5c_write_dataset_list( obj->post_h5_file, dataset_name, num_components, num_elements, list_data_type, ( const void **)list_data_values );

  for ( int idx = 0; idx < num_nodes; idx++ ) {
    free( list_data_values[ idx + 1 ] );
    list_data_values[ idx + 1 ] = NULL;
  }
  free( list_data_values );
  list_data_values = NULL;

  free( list_data_type);
  list_data_type = NULL;

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

int GiD_WriteSphere_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r)
{
  int intvalues[3];
  double doublevalues[1];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=0;
  doublevalues[0]=r;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteSphereMat_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r, int mat)
{
  int intvalues[3];
  double doublevalues[1];
  intvalues[0]=id;
  intvalues[1]=nid;
  intvalues[2]=mat;
  doublevalues[0]=r;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteCircle_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r,double nx, double ny, double nz)
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
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteCircleMat_HDF5( CurrentHdf5WriteData *obj, int id, int nid, double r,double nx, double ny, double nz, int mat)
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
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id,intvalues,doublevalues);
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

int GiD_OpenPostResultFile_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * FileName) {
  if ( !obj ) {
    return -1;
  }
#ifdef _WIN32
  wchar_t wname[MAX_PATH];
  char externalname[MAX_PATH];
#endif  
  int ret;
  ret = initGlobalHdf5File( obj );
  if ( ret < 0 )
    return -1;
  obj->current_mesh_num=0;
  obj->current_gauss_points_num=0;
  obj->current_range_table_num=0;
  obj->num_results_total=0;
  obj->curr_result_group=0;
  obj->num_results_group=0;
  obj->mesh_group_data = new_MeshGroupData();
  int fail = resetMeshGroupData( obj->mesh_group_data );
  if ( fail )
    return -1;
#ifdef _WIN32  
  ret=hdf5c_f_init( obj->post_h5_file, FileName);
  if(ret<0){
    /* try again converting from utf-8 to external*/
    MultiByteToWideChar(CP_UTF8,0,FileName,-1,wname,MAX_PATH);    
    WideCharToMultiByte(CP_ACP,0,wname,-1,externalname,MAX_PATH,0,0);  
    ret=hdf5c_f_init( obj->post_h5_file, externalname);
  }
#else  
  ret=hdf5c_f_init( obj->post_h5_file, FileName);
#endif

  
  if(ret>=0){
    hdf5c_f_set_attribute( obj->post_h5_file, "/","GiD Post Results File","1.1");
    hdf5c_f_set_attribute( obj->post_h5_file, "/","WriteStatus","Writing");
    return 0;
  }
  return ret;
}

int GiD_ClosePostResultFile_HDF5( CurrentHdf5WriteData *obj )
{
  if ( !obj ) {
    return -1;
  }
  int ret;
  hdf5c_f_set_attribute( obj->post_h5_file, "/","WriteStatus","Finished");
  ret=hdf5c_f_end( obj->post_h5_file);
  if(ret<0) 
    return ret;
  ret = endGlobalHdf5File( obj );
  if ( ret < 0 )
    return ret;
  return 0;
}

int GiD_BeginGaussPoint_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * name, GiD_ElementType EType,GP_CONST char * MeshName,
  int GP_number, int NodesIncluded, int InternalCoord)
{
  char gpN[1024],buf[1024];
  char* enames[]={"NoElement","Point","Linear","Triangle","Quadrilateral","Tetrahedra","Hexahedra","Prism","Pyramid",
      "Sphere","Circle"};

  hdf5c_f_create_group( obj->post_h5_file, "GaussPoints");
  obj->current_gauss_points_num++;
  snprintf(gpN, 1024, "GaussPoints/%d",obj->current_gauss_points_num);
  hdf5c_f_create_group( obj->post_h5_file, gpN);
  hdf5c_f_set_attribute( obj->post_h5_file, gpN,"Name",name);
  hdf5c_f_set_attribute( obj->post_h5_file, gpN,"ElemType",enames[EType]);
  if(MeshName) hdf5c_f_set_attribute( obj->post_h5_file, gpN,"MeshName",MeshName);
  snprintf(buf, 1024, "%d",GP_number);
  hdf5c_f_set_attribute( obj->post_h5_file, gpN,"GP_number",buf);
  snprintf(buf, 1024, "%d",NodesIncluded);
  hdf5c_f_set_attribute( obj->post_h5_file, gpN,"NodesIncluded",buf);
  snprintf(buf, 1024, "%d",InternalCoord);
  hdf5c_f_set_attribute( obj->post_h5_file, gpN,"InternalCoord",buf);

  if(InternalCoord==0){
    snprintf(gpN, 1024, "GaussPoints/%d",obj->current_gauss_points_num);
    obj->current_gauss_points_internal_coord=0;
    obj->current_mesh_dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, gpN,0,3);
    return obj->current_mesh_dataset_id;
  }
  else {
    obj->current_gauss_points_internal_coord=1;
    return 0;
  }
}

int GiD_EndGaussPoint_HDF5( CurrentHdf5WriteData *obj )
{
  if(obj->current_gauss_points_internal_coord==0)
    return hdf5c_f_end_dataset( obj->post_h5_file, obj->current_mesh_dataset_id);
  else
    return 0;
}

int GiD_WriteGaussPoint2D_HDF5( CurrentHdf5WriteData *obj, double x, double y)
{
  int intvalues[1] = { 1};
  double doublevalues[3];
  doublevalues[0]=x; doublevalues[1]=y; doublevalues[2]=0.0;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_WriteGaussPoint3D_HDF5( CurrentHdf5WriteData *obj, double x, double y, double z)
{
  int intvalues[1] = { 1};
  double doublevalues[3];
  doublevalues[0]=x; doublevalues[1]=y; doublevalues[2]=z;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->current_mesh_dataset_id,intvalues,doublevalues);
  return 0;
}

int GiD_BeginRangeTable_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * name)
{
  char rtN[1024];  
  hdf5c_f_create_group( obj->post_h5_file, "ResultRangesTable");
  obj->current_range_table_num++;
  obj->current_range_table_idx_num=0;
  snprintf(rtN, 1024, "ResultRangesTable/%d",obj->current_range_table_num);
  hdf5c_f_create_group( obj->post_h5_file, rtN);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Name",name);
  return 0;
}

int GiD_EndRangeTable_HDF5( CurrentHdf5WriteData *obj )
{
  return 0;
}
 
int GiD_WriteMinRange_HDF5( CurrentHdf5WriteData *obj, double max, GP_CONST char * name)
{
  char rtN[1024],buf[1024];
  obj->current_range_table_idx_num++;
  snprintf(rtN, 1024, "ResultRangesTable/%d/%d",obj->current_range_table_num,obj->current_range_table_idx_num);
  hdf5c_f_create_group( obj->post_h5_file, rtN);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Name",name);
  snprintf(buf, 1024, GiD_PostGetFormatReal(),max);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Max",buf);
  return 0;
}

int GiD_WriteRange_HDF5( CurrentHdf5WriteData *obj, double min, double max, GP_CONST char * name)
{
  char rtN[1024],buf[1024];
  obj->current_range_table_idx_num++;
  snprintf(rtN, 1024, "ResultRangesTable/%d/%d",obj->current_range_table_num,obj->current_range_table_idx_num);
  hdf5c_f_create_group( obj->post_h5_file, rtN);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Name",name);
  snprintf(buf, 1024, GiD_PostGetFormatReal(),min);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Min",buf);
  snprintf(buf, 1024, GiD_PostGetFormatReal(),max);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Max",buf);
  return 0;
}

int GiD_WriteMaxRange_HDF5( CurrentHdf5WriteData *obj, double min, GP_CONST char * name)
{
  char rtN[1024],buf[1024];
  obj->current_range_table_idx_num++;
  snprintf(rtN, 1024, "ResultRangesTable/%d/%d",obj->current_range_table_num,obj->current_range_table_idx_num);
  hdf5c_f_create_group( obj->post_h5_file, rtN);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Name",name);
  snprintf(buf, 1024, GiD_PostGetFormatReal(),min);
  hdf5c_f_set_attribute( obj->post_h5_file, rtN,"Min",buf);
  return 0;
}

typedef enum {
  RH_start,
  RH_add
} ResultHeaderMode;

int _GiD_BeginResultHeader_HDF5_init( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis, double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName,ResultHeaderMode RHmode)
{  
  char resN[2048],buf[2048];
  char* rtnames[]={"Scalar","Vector","Matrix","PlainDeformationMatrix","MainMatrix","LocalAxes", 
		   "ComplexScalar", "ComplexVector"};
  
  if(RHmode==RH_start){
    obj->num_results_group=1;
  } else {
    obj->num_results_group++;
  }
  obj->curr_result_group=0;

  obj->myresults[obj->num_results_group-1].num=++obj->num_results_total;
  obj->myresults[obj->num_results_group-1].dataset_id=-1;
  
  char dst_results_path[ 1024 ];
  const char *results_path = getCurrentResultPath( obj, dst_results_path, 1024 );
  snprintf(resN, 2048, "%s/%d", results_path, obj->num_results_total);
  strcpy(obj->myresults[obj->num_results_group-1].name,resN);

  hdf5c_f_create_group( obj->post_h5_file, results_path);
  hdf5c_f_create_group( obj->post_h5_file, resN);

  hdf5c_f_set_attribute( obj->post_h5_file, resN,"Name",Result);
  hdf5c_f_set_attribute( obj->post_h5_file, resN,"Analysis",Analysis);

  snprintf(buf, 2048, GiD_PostGetFormatStep(),step);
  hdf5c_f_set_attribute( obj->post_h5_file, resN,"Step",buf);
  
  hdf5c_f_set_attribute( obj->post_h5_file, resN,"ResultType",rtnames[Type]);
  switch(Where){
    case GiD_OnNodes: 
      hdf5c_f_set_attribute( obj->post_h5_file, resN,"ResultLocation","OnNodes"); 
      break;
    case GiD_OnGaussPoints: 
      hdf5c_f_set_attribute( obj->post_h5_file, resN,"ResultLocation","OnGaussPoints"); 
      break;
    case GiD_OnNurbsLine: 
      hdf5c_f_set_attribute( obj->post_h5_file, resN,"ResultLocation","OnNurbsLine"); 
      break;
    case GiD_OnNurbsSurface: 
      hdf5c_f_set_attribute( obj->post_h5_file, resN,"ResultLocation","OnNurbsSurface"); 
      break;
    case GiD_OnNurbsVolume: 
      hdf5c_f_set_attribute( obj->post_h5_file, resN,"ResultLocation","OnNurbsVolume"); 
      break;
  }
  if(GaussPointsName && GaussPointsName[ 0]) {
    hdf5c_f_set_attribute( obj->post_h5_file, resN,"GaussPointsName",GaussPointsName);
  }
  return 0;
}

int GiD_BeginResult_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis,double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName,
  GP_CONST char * RangeTable, 
  int compc, GP_CONST char * compv[])
{
  int fail,i;
  char* resN,buf[2048];
  
  fail=_GiD_BeginResultHeader_HDF5_init (obj, Result,Analysis,step,Type,Where,GaussPointsName,RH_start);
  if(fail<0) return 1;
  resN=obj->myresults[obj->num_results_group-1].name;

  if ( RangeTable && *RangeTable) {
    hdf5c_f_set_attribute( obj->post_h5_file, resN,"RangeTable",RangeTable);
  }
  if(compc){
    snprintf(buf, 2048, "%d",compc);
    hdf5c_f_set_attribute( obj->post_h5_file, resN,"NumComponents",buf);
  }
  for(i=0;i<compc;i++){
    snprintf(buf, 2048, "Component %d",i+1);
    hdf5c_f_set_attribute( obj->post_h5_file, resN,buf,compv[i]);
  }
  return 0;
}


int GiD_BeginResultHeader_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis, double step,
  GiD_ResultType Type, GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName)
{
  int fail;
  fail=_GiD_BeginResultHeader_HDF5_init( obj, Result,Analysis,step,Type,Where,GaussPointsName,RH_start);
  if(fail<0) return 1;
  return 0;
}

int GiD_BeginResultGroup_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Analysis, double step,GiD_ResultLocation Where,
  GP_CONST char * GaussPointsName)
{
  strcpy(obj->curr_my_result_group.Analysis,Analysis);
  obj->curr_my_result_group.step=step;
  obj->curr_my_result_group.Where=Where;
  if ( GaussPointsName) {
    strcpy( obj->curr_my_result_group.GaussPointsName,GaussPointsName);
  } else {
    obj->curr_my_result_group.GaussPointsName[ 0] = '\0';
  }
  
  obj->num_results_group=0;
  return 0;
}

int GiD_ResultDescription_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GiD_ResultType Type)
{
  int fail;
  fail=_GiD_BeginResultHeader_HDF5_init( obj, Result,obj->curr_my_result_group.Analysis,obj->curr_my_result_group.step,
    Type,obj->curr_my_result_group.Where,obj->curr_my_result_group.GaussPointsName,RH_add);
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
int GiD_ResultLocalAxes_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * Result, GP_CONST char * Analysis,
		                  double step,double vx,double vy,double vz)
{
  char bufres[2048];
  char* resN;
  if( obj->local_axes_format_created == 0){
    snprintf( obj->local_axes_format, LOCAL_AXES_MAX_LEN,"%s,%s,%s",GiD_PostGetFormatReal(),GiD_PostGetFormatReal(),GiD_PostGetFormatReal());
    obj->local_axes_format_created = 1;
  }
  resN=obj->myresults[obj->num_results_group-1].name;
  if(Analysis){
    hdf5c_f_set_attribute( obj->post_h5_file, resN,"LocalAxes Analysis",Analysis);
  }
  hdf5c_f_set_attribute( obj->post_h5_file, resN,"LocalAxes Result",Result);
  snprintf(bufres, 2048, GiD_PostGetFormatStep(),step);
  hdf5c_f_set_attribute( obj->post_h5_file, resN,"LocalAxes step",bufres);
  snprintf(bufres, 2048, obj->local_axes_format,vx,vy,vz);
  hdf5c_f_set_attribute( obj->post_h5_file, resN,"LocalAxes vector",bufres);
  return 0;
}

int GiD_ResultIsDeformationVector_HDF5( CurrentHdf5WriteData *obj, int boolean)
{
  char* resN=obj->myresults[obj->num_results_group-1].name;

  if(boolean) hdf5c_f_set_attribute( obj->post_h5_file, resN,"IsDeformationVector","1");
  else hdf5c_f_set_attribute( obj->post_h5_file, resN,"IsDeformationVector","0");
  return 0;
}

int GiD_ResultRange_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * RangeTable)
{
  char* resN=obj->myresults[obj->num_results_group-1].name;
  if(RangeTable) hdf5c_f_set_attribute( obj->post_h5_file, resN,"RangeTable",RangeTable); 
  return 0;
}

int GiD_ResultComponents_HDF5( CurrentHdf5WriteData *obj, int compc, GP_CONST char * compv[])
{
  int i;
  char buf[2048];
  char* resN=obj->myresults[obj->num_results_group-1].name;
  if(compc){
    snprintf(buf, 2048, "%d",compc);
    hdf5c_f_set_attribute( obj->post_h5_file, resN,"NumComponents",buf); 
  }
  for(i=0;i<compc;i++){
    snprintf(buf, 2048, "Component %d",i+1);
    hdf5c_f_set_attribute( obj->post_h5_file, resN,buf,compv[i]); 
  }
  return 0;
}

int GiD_ResultUnit_HDF5( CurrentHdf5WriteData *obj, GP_CONST char * UnitName)
{
  char* resN=obj->myresults[obj->num_results_group-1].name;
  if(UnitName) hdf5c_f_set_attribute( obj->post_h5_file, resN,"UnitName",UnitName); 
  return 0;
}

int GiD_ResultUserDefined_HDF5( CurrentHdf5WriteData *obj, GP_CONST char* Name,GP_CONST char* Value)
{
  char* resN=obj->myresults[obj->num_results_group-1].name;
  hdf5c_f_set_attribute( obj->post_h5_file, resN,Name,Value);
  return 0;
}

int GiD_ResultValues_HDF5( CurrentHdf5WriteData *obj )
{
  return 0;
}

int GiD_EndResult_HDF5( CurrentHdf5WriteData *obj )
{
  int i,fail;
  for(i=0;i<obj->num_results_group;i++){
    fail=hdf5c_f_end_dataset( obj->post_h5_file, obj->myresults[i].dataset_id);
    if(fail<0) return fail;
  }
  obj->num_results_group=0;
  return 0;
}

int GiD_FlushPostFile_HDF5( CurrentHdf5WriteData *obj )
{
  int ret;
  ret=hdf5c_f_flush( obj->post_h5_file);
  if(ret<0) return ret;
  return 0;
}

int GiD_WriteScalar_HDF5( CurrentHdf5WriteData *obj, int id, double v)
{
  int intvalues[1];
  double doublevalues[1];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,1);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=v;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_Write2DVector_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y)
{
  int intvalues[1];
  double doublevalues[2];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,2);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x;  doublevalues[1]=y;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteVector_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y, double z)
{
  int intvalues[1];
  double doublevalues[3];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,3);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x;  doublevalues[1]=y; doublevalues[2]=z;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteVectorModule_HDF5( CurrentHdf5WriteData *obj, int id, double x, double y, double z, double mod)
{
  int intvalues[1];
  double doublevalues[4];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,4);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x;  doublevalues[1]=y; doublevalues[2]=z; doublevalues[3]=mod;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_Write2DMatrix_HDF5( CurrentHdf5WriteData *obj, int id, double Sxx, double Syy, double Sxy)
{
  int intvalues[1];
  double doublevalues[3];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,3);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Sxx;  doublevalues[1]=Syy; doublevalues[2]=Sxy;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_Write3DMatrix_HDF5( CurrentHdf5WriteData *obj, int id, double Sxx, double Syy, double Szz,
  double Sxy, double Syz, double Sxz)
{
  int intvalues[1];
  double doublevalues[6];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,6);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Sxx;  doublevalues[1]=Syy; doublevalues[2]=Szz;
  doublevalues[3]=Sxy;  doublevalues[4]=Syz; doublevalues[5]=Sxz;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}
 
int GiD_WritePlainDefMatrix_HDF5( CurrentHdf5WriteData *obj, int id, double Sxx, double Syy, double Sxy,
  double Szz)
{
  int intvalues[1];
  double doublevalues[4];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,4);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Sxx;  doublevalues[1]=Syy; doublevalues[2]=Sxy;
  doublevalues[3]=Szz;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteMainMatrix_HDF5( CurrentHdf5WriteData *obj, int id,
  double Si, double Sii, double Siii,
  double Vix, double Viy, double Viz,
  double Viix, double Viiy, double Viiz,
  double Viiix, double Viiiy, double Viiiz)
{
  int intvalues[1];
  double doublevalues[12];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,12);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Si;  doublevalues[1]=Sii; doublevalues[2]=Siii;
  doublevalues[3]=Vix;  doublevalues[4]=Viy; doublevalues[5]=Viz;
  doublevalues[6]=Viix;  doublevalues[7]=Viiy; doublevalues[8]=Viiz;
  doublevalues[9]=Viiix;  doublevalues[10]=Viiiy; doublevalues[11]=Viiiz;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteLocalAxes_HDF5( CurrentHdf5WriteData *obj, int id, double euler_1, double euler_2, double euler_3)
{
  int intvalues[1];
  double doublevalues[3];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,3);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=euler_1;  doublevalues[1]=euler_2; doublevalues[2]=euler_3;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteComplexScalar_HDF5( CurrentHdf5WriteData *obj,  int id, double complex_real, double complex_imag) {
  int intvalues[1];
  double doublevalues[2];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,2);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=complex_real;
  doublevalues[1]=complex_imag;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_Write2DComplexVector_HDF5( CurrentHdf5WriteData *obj,  int id,
				   double x_real, double x_imag,
				   double y_real, double y_imag) {
  int intvalues[1];
  double doublevalues[4];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,4);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x_real;
  doublevalues[1]=x_imag;
  doublevalues[2]=y_real;
  doublevalues[3]=y_imag;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteComplexVector_HDF5( CurrentHdf5WriteData *obj,  int id,
				 double x_real, double x_imag,
				 double y_real, double y_imag,
				 double z_real, double z_imag) {
  int intvalues[1];
  double doublevalues[6];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,6);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=x_real;
  doublevalues[1]=x_imag;
  doublevalues[2]=y_real;
  doublevalues[3]=y_imag;
  doublevalues[4]=z_real;
  doublevalues[5]=z_imag;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_Write2DComplexMatrix_HDF5( CurrentHdf5WriteData *obj, int id,
                                  double Sxx_real, double Syy_real, double Sxy_real,
                                  double Sxx_imag, double Syy_imag, double Sxy_imag) {
  int intvalues[1];
  double doublevalues[6];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,6);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[0]=Sxx_real;  doublevalues[1]=Syy_real; doublevalues[2]=Sxy_real;
  doublevalues[3]=Sxx_imag;  doublevalues[4]=Syy_imag; doublevalues[5]=Sxy_imag;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteComplexMatrix_HDF5( CurrentHdf5WriteData *obj, int id,
                                double Sxx_real, double Syy_real, double Szz_real,
                                double Sxy_real, double Syz_real, double Sxz_real,
                                double Sxx_imag, double Syy_imag, double Szz_imag,
                                double Sxy_imag, double Syz_imag, double Sxz_imag) {
  int intvalues[1];
  double doublevalues[12];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,12);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;
  doublevalues[ 0]=Sxx_real;  doublevalues[ 1]=Syy_real; doublevalues[ 2]=Szz_real;
  doublevalues[ 3]=Sxy_real;  doublevalues[ 4]=Syz_real; doublevalues[ 5]=Sxz_real;
  doublevalues[ 6]=Sxx_imag;  doublevalues[ 7]=Syy_imag; doublevalues[ 8]=Szz_imag;
  doublevalues[ 9]=Sxy_imag;  doublevalues[10]=Syz_imag; doublevalues[11]=Sxz_imag;
  hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,doublevalues);
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

int GiD_WriteNurbsSurface_HDF5( CurrentHdf5WriteData *obj, int id, int num_control_points, double* v) {
  int i_control_point;
  int intvalues[1];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,1);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;/* repeat num_control_points times the surface id to be analogous to OnGaussPoints values */
  for(i_control_point=0;i_control_point<num_control_points;i_control_point++){
    hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,&v[i_control_point]);
  }
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

/* v has components interleaved, e.g for 3 components {v1x v1y v1z ... vnum_control_pointsx vnum_control_pointsy vnum_control_pointsz} */
int GiD_WriteNurbsSurfaceVector_HDF5( CurrentHdf5WriteData *obj,  int id, int num_control_points, int num_comp, double* v ) {
  int i_control_point;
  int intvalues[1];
  char* resN=obj->myresults[obj->curr_result_group].name;
  if(obj->myresults[obj->curr_result_group].dataset_id==-1){
    obj->myresults[obj->curr_result_group].dataset_id=hdf5c_f_start_dataset( obj->post_h5_file, resN,1,num_comp);
    if(obj->myresults[obj->curr_result_group].dataset_id==-1) return -1;
  }
  intvalues[0]=id;/* repeat num_control_points times the surface id to be analogous to OnGaussPoints values */
  for(i_control_point=0;i_control_point<num_control_points;i_control_point++){
    hdf5c_f_add_to_dataset( obj->post_h5_file, obj->myresults[obj->curr_result_group].dataset_id,intvalues,&v[i_control_point*num_comp]);
  }
  obj->curr_result_group++;
  if(obj->curr_result_group==obj->num_results_group) obj->curr_result_group=0;
  return 0;
}

/* MeshGroups: */
int AddMeshGroupName( MeshGroupData *obj, const char *name) {
  if ( !name) return -1;
  int ret_idx = obj->num_mesh_groups;
  /* may be we should check for repetitions ... */
  obj->lst_mesh_groups_names = ( char **)realloc( obj->lst_mesh_groups_names, ( size_t)( obj->num_mesh_groups + 1) * sizeof( char *));
  if ( obj->lst_mesh_groups_names == NULL )
    return -1;
  obj->lst_mesh_groups_names[ ret_idx] = strdup( name);
  obj->num_mesh_groups++;
  return ret_idx ;
}

int GetMeshGroupNameIdx( MeshGroupData *obj, const char *name) {
  if ( !name) return -1;
  int ret_idx = -1;
  int i = 0;
  for ( i = 0; i < obj->num_mesh_groups; i++) {
    if ( !strcmp( obj->lst_mesh_groups_names[ i], name)) {
      ret_idx = i;
      break;
    }
  }
  return ret_idx;
}


int GiD_BeginMeshGroup_HDF5( CurrentHdf5WriteData *obj,  const char *Name) {
  MeshGroupData *mgd = obj->mesh_group_data;
  mgd->current_mesh_group_id = AddMeshGroupName( mgd, Name);
  strcpy( mgd->current_mesh_group_name, Name);
  mgd->using_mesh_group = 1;
  return 0;
}

int GiD_EndMeshGroup_HDF5( CurrentHdf5WriteData *obj ) {
  MeshGroupData *mgd = obj->mesh_group_data;
  mgd->using_mesh_group = 0;
  return 0;
}

int GiD_BeginOnMeshGroup_HDF5( CurrentHdf5WriteData *obj,  const char *Name) {
  MeshGroupData *mgd = obj->mesh_group_data;
  mgd->current_on_mesh_group_id = GetMeshGroupNameIdx( mgd, Name);
  if ( mgd->current_on_mesh_group_id == -1)
    return -1;
  strcpy( mgd->current_mesh_group_name, Name);
  mgd->using_on_mesh_group = 1;
  return 0;
}

int GiD_EndOnMeshGroup_HDF5( CurrentHdf5WriteData *obj ) {
  MeshGroupData *mgd = obj->mesh_group_data;
  mgd->using_on_mesh_group = 0;
  return 0;
}

/* User defined properties defined inside Mesh or Result blocks
   HDF5: stored as properties/attributes of the Mesh/Result folder (Name, value)
   ASCII / raw binary: stored as comments
     # Name: value
   Define the macro COMPASSIS_USER_ATTRIBUTES_FORMAT
   to have it like compassis wants:
     # ResultUserDefined \"%s\" \"%s\"      or
     # ResultUserDefined \"%s\" %s
*/
int GiD_WriteMeshUserAttribute_HDF5( CurrentHdf5WriteData *obj, GP_CONST char *Name, GP_CONST char *Value ) {
  if ( obj->current_mesh_num > 0 ) {
    char meshN[ 1024 ];
    char dest_path[ 1024 ];
    snprintf( meshN, 1024, "%s/%d", getCurrentMeshPath( obj, dest_path, 1024 ), obj->current_mesh_num );
    hdf5c_f_set_attribute( obj->post_h5_file, meshN, Name, Value );
    return 0;
  }
  return -1;
}

int GiD_WriteResultUserAttribute_HDF5( CurrentHdf5WriteData *obj,  GP_CONST char *Name, GP_CONST char *Value ) {
  if ( obj->num_results_group > 0 ) {
    const char *resN = obj->myresults[ obj->num_results_group - 1 ].name;
    hdf5c_f_set_attribute( obj->post_h5_file, resN, Name, Value );
    return 0;
  }
  return -1;
}


// GiD_fWriteResultBlock includes BeginResult()/BeginValues() and EndValues()/EndResult()
int GiD_WriteResultBlock_HDF5( CurrentHdf5WriteData *obj, GP_CONST char *result_name, GP_CONST char *analysis_name,
                               double step_value, GiD_ResultType result_type, GiD_ResultLocation result_location,
                               GP_CONST char *gauss_point_name, GP_CONST char *range_table_name,
                               int num_component_names, GP_CONST char *list_component_names[], 
                               GP_CONST char *unit_name,
                               int num_result_values,
                               GP_CONST int *list_result_ids, int num_component_values,
                               GP_CONST double *list_component_values ) {
  // at least there should be result values
  if ( list_component_values == NULL )
    return -1;

  int _fail = 0;

  _fail = GiD_BeginResult_HDF5( obj, result_name, analysis_name, step_value, result_type, result_location,
                                 gauss_point_name, range_table_name, num_component_names, list_component_names );

  if ( !_fail ) {
    if ( unit_name && *unit_name ) {
      _fail = GiD_ResultUnit_HDF5( obj, unit_name );
    }
  }

  if ( !_fail ) {

    // in GiD Post hdf5 ( = compassis hdf5) data is stored column-wise in different datasets
    // for instance: 1 data set with 1 column for id's
    // 1 data set with 1 column for x-coordinate
    // 2 data set with 1 column for x-coordinate
    // 3 data set with 1 column for x-coordinate
    const int *list_ids = list_result_ids;
    if ( list_ids == NULL ) {
      int *my_list_ids = ( int * )malloc( ( size_t)num_result_values * sizeof( int ) );
      if ( my_list_ids ) {
        for ( int id = 0; id < num_result_values; id++ ) {
          my_list_ids[ id ] = id + 1;
        }
        list_ids = my_list_ids;
      }
    }
    double **list_values = ( double ** )malloc( ( size_t)num_component_values * sizeof( double *));
    if ( list_values == NULL )
      return -1;
    for ( int idx_component = 0; idx_component < num_component_values; idx_component++ ) {
      list_values[ idx_component ] = ( double * )malloc( ( size_t)num_result_values * sizeof( double ) );
      if ( list_values[ idx_component ] == NULL )
        return -1;
    }

    if ( list_values ) {
      int idx_array_result_values = 0;
      for ( int id = 0; id < num_result_values; id++ ) {
        for ( int idx_component = 0; idx_component < num_component_values; idx_component++ ) {
          list_values[ idx_component ][ id] = list_component_values[ idx_array_result_values ];
          idx_array_result_values++;
        }
      }

      char dataset_name[ 2048 ];

      
      char dst_results_path[ 1024 ];
      const char *results_path = getCurrentResultPath( obj, dst_results_path, 1024 );
      snprintf( dataset_name, 2048, "%s/%d", results_path, obj->num_results_total );
      

      const int num_datasets = num_component_values + 1; // id + Scalar | id + Vx + Vy + Vz  | id + Sxx + Syy + ...
      t_hdf5c_dataset_type *lst_types = ( t_hdf5c_dataset_type *)malloc( ( size_t) num_datasets * sizeof( t_hdf5c_dataset_type));
      const void **lst_data = ( const void ** )malloc( ( size_t )num_datasets * sizeof( const void * ) );
      if ( lst_types && lst_data ) {
        lst_types[ 0 ] = HDF5C_DATASET_INTEGER;
        lst_data[ 0 ] = list_ids;
        for ( int idx_dataset = 1; idx_dataset < num_datasets; idx_dataset++ ) {
          lst_types[ idx_dataset ] = HDF5C_DATASET_DOUBLE;
          lst_data[ idx_dataset ] = ( const void *)( list_values[ idx_dataset - 1 ]);
        }
      } else {
        return -1;
      }
      hdf5c_write_dataset_list( obj->post_h5_file, dataset_name, num_datasets, num_result_values, lst_types, lst_data );
      free( lst_data );
      lst_data = NULL;
      free( lst_types );
      lst_types = NULL;
      if ( list_ids != list_result_ids ) {
        // it was allocated by us
        int *my_list_ids = ( int *)list_ids;
        free( my_list_ids );
      }
      list_ids = NULL;
      for ( int idx_component = 0; idx_component < num_component_values; idx_component++ ) {
        free( list_values[ idx_component ] );
        list_values[ idx_component ] = NULL;
      }
      free( list_values );
      list_values = NULL;
    } else {
      return -1;
    }
    return 0;
  }

  //  no need for GiD_EndResult_HDF5( obj ); as it does nothing (i.e. dataset has been already written above

  return _fail;
}

// #endif // ENABLE_HDF5
