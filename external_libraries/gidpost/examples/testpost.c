#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "gidpost.h"

#ifdef WIN32
#define strcasecmp  _stricmp
#define strdup  _strdup
#include <windows.h>   // needed to include correctly synchapi.h
#include <synchapi.h>  // for sleep()
#define sleep( sec)   Sleep( sec * 1000)
#else // WIN32
#include <strings.h>
#include <unistd.h>  // for sleep()
#endif // WIN32

typedef struct {
  int id;
  float x, y, z;
} SCoord;

typedef struct {
  int id;
  int n1, n2, n3;
} SElem;

SCoord G_nodes[9 * 2]; // 1st 9 nodes and 2nd 9 nodes are the same
SElem G_elems[8 * 2]; // 1st 8 elements and 2nd 8 elements are the same

#define NUM_NODES_2   9
#define NUM_ELEMS_2   8
#define NUM_NODES    18
#define NUM_ELEMS    16

void GenMesh()
{
  int i, j, idx;

  idx = 1;
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++, idx++ ) {
      SCoord *node = G_nodes + idx -1;
      node->id = idx;
      node->x = (float)(i);
      node->y = (float)(j);
      node->z = 0.0f;
    }
  }
  for ( i = 0; i < NUM_NODES_2; i++) {
    j = i + NUM_NODES_2;
    G_nodes[ j] = G_nodes[ i];
    G_nodes[ j].id = j + 1;
    G_nodes[ j].z += 3.0f;
  }
  idx = 0;
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++, idx+=2 ) {
      SElem *elem1 = G_elems + idx;
      SElem *elem2 = elem1 + 1;
      elem1->id = idx+1;
      elem1->n1 = i*3 + j + 1;
      elem1->n3 = i*3 + j + 1 + 3 ;
      elem1->n2 = elem1->n3 + 1;
      elem2->id = idx+2;
      elem2->n1 = elem1->n1;
      elem2->n2 = elem1->n1 + 1;
      elem2->n3 = elem1->n2;
    }
  }
  for ( i = 0; i < NUM_ELEMS_2; i++) {
    j = i + NUM_ELEMS_2;
    G_elems[ j] = G_elems[ i];
    G_elems[ j].id += NUM_ELEMS_2;
    G_elems[ j].n1 += NUM_NODES_2;
    G_elems[ j].n2 += NUM_NODES_2;
    G_elems[ j].n3 += NUM_NODES_2;
  }
}

float Random()
{
  return rand()/(float)(RAND_MAX);
}

void write_results( GiD_FILE fdr, int use_mesh_groups, const char *analysis, double step_value, int last_elem_id, const char *range_table_name ) {

  if ( use_mesh_groups )
    GiD_fBeginOnMeshGroup( fdr, "MyGroup" );

  GiD_fBeginResult( fdr, "debug nodes", analysis, step_value, GiD_Scalar, GiD_OnNodes, NULL, NULL, 0, NULL );
  GiD_fWriteResultUserAttribute( fdr, "myResultAttr", "my Result value" );
  int i = 0;
  int j = 0;
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWriteScalar( fdr, G_nodes[ i ].id, ( double )( G_nodes[ i ].id ) );
  }
  GiD_fEndResult( fdr);

  GiD_fBeginResult( fdr, "debug elements", analysis, step_value, GiD_Scalar, GiD_OnGaussPoints,
                   "element_1gp", NULL, 0, NULL );
  // writing results for 1st and 2nd meshes
  for ( i = 0; i < NUM_ELEMS; i++ ) {
    for ( j = 0; j < 1; j++ )
      GiD_fWriteScalar( fdr, G_elems[ i ].id, ( double )G_elems[ i ].id );
  }
  // writing results for spheres and circles
  for ( i = NUM_ELEMS; i < last_elem_id; i++) {
    GiD_fWriteScalar( fdr, i + 1, ( double )( i + 1 ) );
  }
  GiD_fEndResult( fdr);

  GiD_fBeginResult( fdr, "almost node_scalar (r)", analysis, step_value, GiD_Scalar, GiD_OnNodes,
                   NULL, range_table_name, 0, NULL );
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWriteScalar( fdr, G_nodes[ i ].id, Random() );
  }
  GiD_fEndResult( fdr);

  /* scalar result over G_nodes */

  // result group over nodes
  GiD_fBeginResultGroup( fdr, analysis, step_value, GiD_OnNodes, NULL );
  GiD_fResultDescription( fdr, "nodal//ScalarG", GiD_Scalar );
  GiD_fResultUnit( fdr, "m" );
  GiD_fResultRange( fdr, range_table_name );
  GiD_fResultDescriptionDim( fdr, "nodal//VectorG", GiD_Vector, 4 );
  GiD_fResultUnit( fdr, "m/s" );
  GiD_fResultDescription( fdr, "nodal//MatrixG", GiD_Matrix );
  GiD_fResultDescription( fdr, "nodal//Local AxesG", GiD_LocalAxes );
  GiD_fResultUnit( fdr, "N·m" );
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWriteScalar( fdr, G_nodes[ i ].id, Random() );
    GiD_fWriteVectorModule( fdr, G_nodes[ i ].id, Random(), Random(), Random(), -1 );
    GiD_fWrite3DMatrix( fdr, G_nodes[ i ].id, Random(), Random(), Random(),
                       Random(), Random(), Random() );
    GiD_fWriteLocalAxes( fdr, G_nodes[ i ].id, Random(), Random(), Random() );
  }
  GiD_fEndResult( fdr);

  // result group over gauss points
  GiD_fBeginResultGroup( fdr, analysis, step_value, GiD_OnGaussPoints, "element_3gp" );
  GiD_fResultDescription( fdr, "gaussian//ScalarG", GiD_Scalar );
  GiD_fResultUnit( fdr, "m" );
  GiD_fResultRange( fdr, range_table_name );
  GiD_fResultDescriptionDim( fdr, "gaussian//VectorG", GiD_Vector, 4 );
  GiD_fResultUnit( fdr, "m/s" );
  GiD_fResultDescription( fdr, "gaussian//MatrixG", GiD_Matrix );
  GiD_fResultDescription( fdr, "gaussian//Local AxesG", GiD_LocalAxes );
  GiD_fResultUnit( fdr, "N·m" );
  for ( i = 0; i < NUM_ELEMS; i++ ) {
    for ( j = 0; j < 3; j++ ) { // 3 gauss points per element
      GiD_fWriteScalar( fdr, G_elems[ i ].id, Random() );
      GiD_fWriteVectorModule( fdr, G_elems[ i ].id, Random(), Random(), Random(), -1 );
      GiD_fWrite3DMatrix( fdr, G_elems[ i ].id, Random(), Random(), Random(),
                         Random(), Random(), Random() );
      GiD_fWriteLocalAxes( fdr, G_elems[ i ].id, Random(), Random(), Random() );
    }
  }
  GiD_fEndResult( fdr);

  /* scalar results over gauss points: Triangle + 3 GP */

  GiD_fBeginResult( fdr, "Over_GP Triangle3", analysis, step_value, GiD_Scalar, GiD_OnGaussPoints,
                   "element_3gp", range_table_name, 0, NULL );
  for ( i = 0; i < NUM_ELEMS; i++ ) {
    for ( j = 0; j < 3; j++ )
      GiD_fWriteScalar( fdr, G_elems[ i ].id, Random() );
  }
  GiD_fEndResult( fdr);
  GiD_fBeginResult( fdr, "Nodal_Scalar with spaces ( r2)",
                   analysis, step_value, GiD_Scalar, GiD_OnNodes,
                   NULL, range_table_name, 0, NULL );
  for ( i = 0; i < NUM_NODES; i++ ) {
    float rr = Random();
    GiD_fWriteScalar( fdr, G_nodes[ i ].id, rr * rr );
  }
  GiD_fEndResult( fdr);

  /* vector result over G_nodes */

  GiD_fBeginResultHeader( fdr, "VELOCITIES (m/s)", analysis, step_value, GiD_Vector, GiD_OnNodes, NULL );
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWriteVector( fdr, G_nodes[ i ].id, Random(), Random(), Random() );
  }
  GiD_fEndResult( fdr);

  /* matrix result */

  GiD_fBeginResult( fdr, "Matrix", analysis, step_value, GiD_Matrix, GiD_OnNodes,
                   NULL, NULL, 0, NULL );
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWrite3DMatrix( fdr, G_nodes[ i ].id, Random(), Random(), Random(),
                       Random(), Random(), Random() );
  }
  GiD_fEndResult( fdr);

  /* local axes result */

  GiD_fBeginResult( fdr, "Local Axes", analysis, step_value, GiD_LocalAxes, GiD_OnNodes,
                   NULL, NULL, 0, NULL );
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWriteLocalAxes( fdr, G_nodes[ i ].id, Random(), Random(), Random() );
  }
  GiD_fEndResult( fdr);

  GiD_fBeginResult( fdr, "Complex Scalar", analysis, step_value, GiD_ComplexScalar, GiD_OnNodes,
                   NULL, NULL, 0, NULL );
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWriteComplexScalar( fdr, G_nodes[ i ].id, Random(), Random() );
  }
  GiD_fEndResult( fdr);

  GiD_fBeginResult( fdr, "Complex Vector (last result)", analysis, step_value, GiD_ComplexVector, GiD_OnNodes,
                   NULL, NULL, 0, NULL );
  for ( i = 0; i < NUM_NODES; i++ ) {
    GiD_fWriteComplexVector( fdr, G_nodes[ i ].id,
                            Random(), Random(), Random(),
                            Random(), Random(), Random() );
  }
  GiD_fEndResult( fdr);

  if ( use_mesh_groups )
    GiD_fEndOnMeshGroup( fdr);

  GiD_fFlushPostFile( fdr);
}

// returns a malloc'ed vector with created dummy results which is to be free'd outside
static double *create_dummy_result_block( int num_components, int num_values, int offset_value ) {
  int num_total_values = num_components * num_values;
  double *values = ( double * )malloc( ( size_t)( num_total_values ) * sizeof( double ) );
  if ( values == NULL )
    return NULL;
  // define predictable results
  int offset_values = 0;
  for ( int idx_result = 0; idx_result < num_values; idx_result++ ) {
    for ( int idx_component = 0; idx_component < num_components; idx_component++ ) {
      double my_value = 100 * idx_result + idx_component + 1;
      values[ offset_values ] = my_value + offset_value;
      offset_values++;
    }
  }
  return values;
}

#define C_NEW( type, size)   ( type *)malloc( ( size_t)( size) * sizeof( type))
#define C_DELETE( ptr)       { free( ptr); ptr = NULL;}

void write_results_blocks( GiD_FILE fdr, int use_mesh_groups, const char *analysis, double step_value, int last_elem_id, const char *range_table_name ) {
    int list_node_ids[ NUM_NODES];
    int *list_elem_ids = C_NEW( int, last_elem_id);
    const int num_gp_per_element = 3;
    int *list_elems_ids_3gp = C_NEW( int, last_elem_id * num_gp_per_element);
    double list_node_scalar[ NUM_NODES];
    double *list_elem_scalar = C_NEW( double, last_elem_id);
    for ( int in = 0; in < NUM_NODES; in++) {
        list_node_ids[ in] = G_nodes[ in].id;
        list_node_scalar[ in] = ( double)G_nodes[ in].id;
    }
    int idx_elem_3gp = 0;
    for ( int ie = 0; ie < NUM_ELEMS; ie++) {
        list_elem_ids[ ie] = G_elems[ ie].id;
        list_elem_scalar[ ie] = ( double)G_elems[ ie].id;
        for ( int igp = 0; igp < num_gp_per_element; igp++) {
            list_elems_ids_3gp[ idx_elem_3gp] = G_elems[ ie].id;
            idx_elem_3gp++;
        }
    }
    for ( int ie = NUM_ELEMS; ie < last_elem_id; ie++) {
        list_elem_ids[ ie] = ie + 1;
        list_elem_scalar[ ie] = ( double)( ie + 1);
    }


    char *unit_name = "m"; // just an example
    // GiD_fWriteResultBlock includes BeginResult()/BeginValues() and EndValues()/EndResult()
    GiD_fWriteResultBlock( fdr, "debug nodes", analysis, step_value, GiD_Scalar, GiD_OnNodes, NULL, 
                           NULL, 0, NULL, unit_name,
                           NUM_NODES, list_node_ids, 1, list_node_scalar );
    GiD_fWriteResultBlock( fdr, "debug elements", analysis, step_value, GiD_Scalar, GiD_OnGaussPoints, "element_1gp",  
                           NULL, 0, NULL, unit_name,
                           last_elem_id, list_elem_ids, 1, list_elem_scalar );

    int num_components = 1;
    int offset_value = ( int)( 0.5f + step_value * 1000.0);
    double *nodal_scalar = create_dummy_result_block( num_components, NUM_NODES, offset_value );
    GiD_fWriteResultBlock( fdr, "almost node_scalar (r)", analysis, step_value, GiD_Scalar, GiD_OnNodes, NULL,
                           range_table_name, 0, NULL, NULL,
                           NUM_NODES, list_node_ids, num_components, nodal_scalar );
    free( nodal_scalar ); nodal_scalar = NULL;

    // result group over nodes
    // printf( "no support for Result Groups when writing ResultBlocks.\n" );
    // result group over gauss points
    // printf( "no support for Result Groups when writing ResultBlocks.\n" );
    num_components = 1;
    double *elem_3gp_scalar = create_dummy_result_block( num_components, NUM_ELEMS * num_gp_per_element, offset_value );
    GiD_fWriteResultBlock( fdr, "Over_GP Triangle3", analysis, step_value, GiD_Scalar, GiD_OnGaussPoints, "element_3gp",
                           range_table_name, 0, NULL, NULL,
                           NUM_ELEMS * num_gp_per_element, list_elems_ids_3gp, num_components, elem_3gp_scalar );
    free( elem_3gp_scalar ); elem_3gp_scalar = NULL;

    num_components = 3; // VECTOR
    unit_name = "m/s";
    double *nodal_vector = create_dummy_result_block( num_components, NUM_NODES, offset_value );
    GiD_fWriteResultBlock( fdr, "VELOCITIES", analysis, step_value, GiD_Vector, GiD_OnNodes, NULL,
                           NULL, 0, NULL, unit_name,
                           NUM_NODES, list_node_ids, num_components, nodal_vector );
    free( nodal_vector ); nodal_vector = NULL;

    num_components = 6; // MATRIX
    double *nodal_matrix = create_dummy_result_block( num_components, NUM_NODES, offset_value );
    GiD_fWriteResultBlock( fdr, "Matrix", analysis, step_value, GiD_Matrix, GiD_OnNodes, NULL,
                           NULL, 0, NULL, NULL,
                           NUM_NODES, list_node_ids, num_components, nodal_matrix );
    free( nodal_matrix ); nodal_matrix = NULL;

    num_components = 3; // Local Axis
    nodal_vector = create_dummy_result_block( num_components, NUM_NODES, offset_value );
    GiD_fWriteResultBlock( fdr, "Local Axis", analysis, step_value, GiD_LocalAxes, GiD_OnNodes, NULL,
                           NULL, 0, NULL, NULL,
                           NUM_NODES, list_node_ids, num_components, nodal_vector );
    free( nodal_vector ); nodal_vector = NULL;

    num_components = 2; // Complex scalar
    nodal_scalar = create_dummy_result_block( num_components, NUM_NODES, offset_value );
    GiD_fWriteResultBlock( fdr, "Complex Scalar", analysis, step_value, GiD_ComplexScalar, GiD_OnNodes, NULL,
                           NULL, 0, NULL, NULL,
                           NUM_NODES, list_node_ids, num_components, nodal_scalar );
    free( nodal_scalar ); nodal_scalar = NULL;

    num_components = 6; // Complex Vector
    nodal_vector = create_dummy_result_block( num_components, NUM_NODES, offset_value );
    GiD_fWriteResultBlock( fdr, "Complex Vector (last result)", analysis, step_value, GiD_ComplexVector, GiD_OnNodes, NULL,
                           NULL, 0, NULL, unit_name,
                           NUM_NODES, list_node_ids, num_components, nodal_vector );
    free( nodal_vector ); nodal_vector = NULL;

    C_DELETE( list_elem_ids );
    C_DELETE( list_elems_ids_3gp );
    C_DELETE( list_elem_scalar);

    GiD_fFlushPostFile( fdr);
}

void print_help_and_exit( const char *full_prog_name) {
  // get only the name of the executable
  const char *prog_name = &full_prog_name[ strlen( full_prog_name ) - 1 ];
  for ( ; prog_name > full_prog_name; prog_name-- ) {
    if ( ( *prog_name == '/' ) || ( *prog_name == '\\' ) ) {
      break;
    }
  }
  prog_name++;
  printf( "Usage: %s [ -h] [ -abort] [ -array] [ -f ascii|bin|hdf5] [ -g ] [ -n number_of_steps_to_write] filename_prefix\n", prog_name);
  printf( "    will create filename_prefix.post.{ msh,res | bin | h5} depending on the choosen format.\n");
  printf( "    Use: %s -- filename_prefix          if filename_prefix begins with a '-'\n", prog_name);
  printf( "            -abort                does not close the file properly, after last result it aborts" );
  printf( "                                      ( to simulate calculation program abort )\n " );
  printf( "            -array                Write Results arrays using GiD_fWriteResultsBlock()\n" );
  printf( "            -g                    Write files within MeshGroup / OnMeshGroup\n" );
  printf( "            -n N                  Write N steps, default 1, should be N >= 0\n" );
  exit( 0);
}

int main( int argc, char *argv[])
{
  int i = 0;
  int last_node_id, last_elem_id;
  const char *analysis = "Analysis_of whatever";
  const char *range_table_name = "simple range table";
#define BUF_SIZE 1024
  char buf[ BUF_SIZE ];
  char test_filename[ BUF_SIZE ];
  char *prefix = NULL;
  char *format = NULL;
  int skip_options = 0;
  int early_terminate = 0;
  int use_mesh_groups = 0;
  int write_arrays = 0;
  int number_of_steps_to_write = 1;

  for ( int ia = 1; ia < argc; ia++) {
    if ( !skip_options && ( argv[ ia][ 0] == '-')) {
      char opt = ( char)tolower( argv[ ia][ 1]);
      if ( ( opt == 'h') || !strcasecmp( argv[ ia], "--help") ) {
        print_help_and_exit( argv[ 0 ] );
      } else if ( opt == 'f' ) {
        ia++;
        if ( ia < argc ) {
          format = strdup( argv[ ia ] );
        } else {
          printf( "Missing arguments.\n" );
          print_help_and_exit( argv[ 0 ] );
        }
      } else if ( opt == 'n') {
        ia++;
        int error = 1;
        if ( ia < argc ) {
          int ns = 0;
          if ( sscanf( argv[ ia ], "%d", &ns ) == 1 ) {
            if ( ns >= 0 ) {
              number_of_steps_to_write = ns;
              error = 0;
            }
          }
        }
        if ( error) {
          printf( "Missing arguments.\n" );
          print_help_and_exit( argv[ 0 ] );
        }
      } else if ( ( opt == 'a') && !strcasecmp( argv[ ia], "-abort")  ) {
        early_terminate = 1;
      } else if ( ( opt == 'a') && !strcasecmp( argv[ ia], "-array")  ) {
        write_arrays = 1;
      } else if ( opt == 'g' ) {
        use_mesh_groups = 1;
      } else if ( opt == '-' ) {
        skip_options = 1;
      } else {
        printf( "Unknown option '%s'.\n", argv[ ia]);
        print_help_and_exit( argv[ 0]);
      }
    } else {
      prefix = strdup( argv[ ia]);
      break;
    }
  }

  printf( "version = %s\n", GiD_PostGetVersion());

  GiD_PostInit();
  GenMesh();

  // default values
  if ( !prefix) prefix = strdup( "test");
  if ( !format) format = strdup( "ascii");

  GiD_FILE fdm = 0, fdr = 0;
  int fail = 0;
  if ( !strcasecmp( format, "ascii")) {
    GiD_PostSetFormatReal("%.8g");
    snprintf( buf, BUF_SIZE, "%s.post.msh", prefix);
    fdm = GiD_fOpenPostMeshFile( buf, GiD_PostAscii );
    fail = ( fdm == 0 ) ? 1 : 0;
    snprintf( buf, BUF_SIZE, "%s.post.res", prefix);
    // fdr = GiD_fOpenPostResultFile( buf, GiD_PostAsciiZipped );
    fdr = GiD_fOpenPostResultFile( buf, GiD_PostAscii );
    fail += ( fdr == 0 ) ? 1 : 0;
    printf( "Creating ASCII gid post files '%s.post.msh' and '%s.post.res'.\n",
            prefix, prefix);
  } else if ( !strcasecmp( format, "bin")) {
    snprintf( buf, BUF_SIZE, "%s.post.bin", prefix);
    fdm = fdr = GiD_fOpenPostResultFile( buf, GiD_PostBinary );
    fail = ( fdm == 0 ) ? 1 : 0;
    printf( "Creating BINary gid post file '%s'.\n", buf);
  } else if ( !strcasecmp( format, "hdf5")) {
    snprintf( buf, BUF_SIZE, "%s.post.h5", prefix);
    printf( "Creating HDF5 gid post file '%s'.\n", buf );
    fdm = fdr = GiD_fOpenPostResultFile( buf, GiD_PostHDF5 );
    fail = ( fdm == 0 ) ? 1 : 0;
  } else {
    printf( "Unkown format '%s'.\n", format);
    print_help_and_exit( argv[ 0]);
  }
  strcpy( test_filename, buf);

  if ( fail) {
    printf( "Error opening file with format '%s', error code = %d.\n", format, fail);
    exit( 1);
  }

  if ( use_mesh_groups )
    GiD_fBeginMeshGroup( fdm, "MyGroup");
  
  /* write mesh info */
  GiD_fBeginMeshColor( fdm, "TestMesh", GiD_2D, GiD_Triangle, 3, 0, 0.99, 0);
  GiD_fWriteMeshUserAttribute( fdm, "myMeshAttr", "myMeshValue" );
  /* coordinates */
  GiD_fBeginCoordinates( fdm);
  last_node_id = 0;
  for ( i = 0; i < NUM_NODES_2; i++ ) {
    const SCoord *node = G_nodes + i;
    last_node_id = last_node_id < node->id ? node->id : last_node_id;
    GiD_fWriteCoordinates( fdm, node->id, node->x, node->y, node->z ); 
  }
  GiD_fEndCoordinates( fdm);
  /* elements */
  GiD_fBeginElements( fdm);
  last_elem_id = 0;
  for ( i = 0; i < NUM_ELEMS_2; i++ ) {
    const SElem *elem = G_elems + i;
    int elemi[4];
    elemi[0] = elem->n1;
    elemi[1] = elem->n2;
    elemi[2] = elem->n3;
    elemi[3] = 1;
    last_elem_id = last_elem_id < elem->id ? elem->id : last_elem_id;
    GiD_fWriteElementMat( fdm, elem->id, elemi);
  }
  GiD_fEndElements( fdm);
  GiD_fEndMesh( fdm);

  /* write mesh info */
  GiD_fBeginMeshColor( fdm, "TestMesh 2", GiD_2D, GiD_Triangle, 3, 0, 0, 0.99);
  /* coordinates */
  GiD_fBeginCoordinates( fdm);
  for ( i = NUM_NODES_2; i < NUM_NODES; i++ ) {
    const SCoord *node = G_nodes + i;
    GiD_fWriteCoordinates( fdm, node->id, node->x, node->y, node->z); 
  }
  GiD_fEndCoordinates( fdm);
  /* elements */
  GiD_fBeginElements( fdm);
  // repeating connectivities
  for ( i = NUM_ELEMS_2; i < NUM_ELEMS; i++ ) {
    const SElem *elem = G_elems + i;
    int elemi[4];
    elemi[0] = elem->n1;
    elemi[1] = elem->n2;
    elemi[2] = elem->n3;
    elemi[3] = 2;
    last_elem_id = last_elem_id < elem->id ? elem->id : last_elem_id;
    GiD_fWriteElementMat( fdm, elem->id, elemi);
  }
  GiD_fEndElements( fdm);
  GiD_fEndMesh( fdm);

  /* write mesh info */
  // defining a mesh pf spheres over 1st mesh
  GiD_fBeginMeshColor( fdm, "Spheres", GiD_3D, GiD_Sphere, 1, 1.0, 0.7, 0.3 );
  /* coordinates */
  GiD_fBeginCoordinates( fdm);
  GiD_fEndCoordinates( fdm);
  /* elements */
  GiD_fBeginElements( fdm);
  for ( i = 0; i < NUM_NODES_2; i++ ) {
    double radius = Random() * 0.3 + 0.3;
    last_elem_id++;
    GiD_fWriteSphere( fdm, last_elem_id, i + 1, radius );
  }
  GiD_fEndElements( fdm);
  GiD_fEndMesh( fdm);

  // defining a mesh of circles over the 2nd mesh
  GiD_fBeginMeshColor( fdm, "Circles", GiD_3D, GiD_Circle, 1, 0.8, 1.0, 0.2 );
  /* coordinates */
  GiD_fBeginCoordinates( fdm);
  GiD_fEndCoordinates( fdm);
  /* elements */
  GiD_fBeginElements( fdm);
  for ( i = 0; i < NUM_NODES_2; i++ ) {
    double radius = Random() * 0.3 + 0.3;
    last_elem_id++;
    GiD_fWriteCircle( fdm, last_elem_id, i + NUM_NODES_2 + 1, radius, 0.0, 0.0, 1.0 );
  }
  GiD_fEndElements( fdm);
  GiD_fEndMesh( fdm);
  if ( use_mesh_groups )
    GiD_fEndMeshGroup( fdm);

  // last_node_id += last_node_id;
  // last_elem_id += last_elem_id;

  /* now results info */
  /* write the gauss points */
  GiD_fBeginGaussPoint( fdr, "element_3gp", GiD_Triangle, NULL, 3, 0, 1);  
  GiD_fEndGaussPoint( fdr);

  GiD_fBeginGaussPoint( fdr, "element_1gp", GiD_Triangle, NULL, 1, 0, 1);  
  GiD_fEndGaussPoint( fdr);

  /* results are between 0.0 and 1.0 */
  GiD_fBeginRangeTable( fdr, range_table_name);
  GiD_fWriteMinRange( fdr, 0.1, "abs min");
  GiD_fWriteRange( fdr, 0.1, 0.4, "min");
  GiD_fWriteRange( fdr, 0.4, 0.6, "middle");
  GiD_fWriteRange( fdr, 0.6, 0.9, "max");
  GiD_fWriteMaxRange( fdr, 0.9, "abs max");
  GiD_fEndRangeTable( fdr);

  for ( int idx_step = 0; idx_step < number_of_steps_to_write; idx_step++ ) {
    const double step_value = idx_step + 1.0;
    if ( write_arrays) {
      printf( "Writing results blocks of step %g\n", step_value );
      write_results_blocks( fdr, use_mesh_groups, analysis, step_value, last_elem_id, range_table_name );
    } else {
      printf( "Writing results of step %g\n", step_value );
      write_results( fdr, use_mesh_groups, analysis, step_value, last_elem_id, range_table_name );
    }
    sleep( 1 ); // sleep 1 second between steps to write

    if ( early_terminate ) {
      printf( "aborting...\n" );
      abort();
    }
  }

  if ( !strcasecmp( format, "ascii")) {
    GiD_fClosePostMeshFile( fdm);
    fdm = 0;
  }
  GiD_fClosePostResultFile( fdr);
  fdr = 0;
  GiD_PostDone();
  printf( "... done.\n");

  printf( "%s written with format '%s'", test_filename, format );
  if ( write_arrays) {
    printf( " using *Block() API" );
  }
  printf( "\n" );
  
  free( prefix); prefix = NULL;
  free( format); format = NULL;
  
  return 0;
}
