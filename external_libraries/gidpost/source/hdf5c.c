
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "hdf5c.h"
#include "hdf5.h"
#include "H5public.h"
#include "H5LTpublic.h"

// #define MAX_CONCURRENT_DATASETS 20

typedef enum _OpenType {
  Read_OT,
  Write_OT,
  Create_OT
} OpenType;

struct hdf5c_buffer {
  int used;
  hid_t file_id;
  char datasetname[ 1024 ];
  int memsize, size, num_int, num_real;
  int *intarray[ 20 ];
  double *doublearray[ 20 ];
};

#ifdef NDEBUG
#define print_error( str)
#else
static void _print_stderr( const char *str ) {
  fprintf( stderr, "gidpost hdf5 error: %s\n", str );
}
#define print_error( str)      _print_stderr( str);
#define HDF5_DEBUG
#endif

// typedef struct _CPostHdf5File CPostHdf5File; is done in hdf5.h
struct _CPostHdf5File {
   hid_t file_id;
   hid_t group_id;
   char group_path[ 1024 ];
   struct hdf5c_buffer file_buffer[ MAX_CONCURRENT_DATASETS ];
   char error_buffer[ 8192 ];
};


const char *hdf5c_f_get_version() {
  static char *version_str = NULL;

  if ( !version_str ) {
    char buf[ 1024 ];
    unsigned major_number, minor_number, patch_number;
    hbool_t is_thread_safe = 0;

    // retrieve the library version
    herr_t h5_error = H5get_libversion( &major_number, &minor_number, &patch_number );
    if ( h5_error >= 0 ) {
      // is this a thread-safe library build?
      h5_error = H5is_library_threadsafe( &is_thread_safe );
    }
    if ( h5_error >= 0 ) {
      snprintf( buf, 1024, "HDF5 %d.%d.%d  ( Thread-safety %s )", major_number, minor_number, patch_number, ( is_thread_safe > 0 ) ? "enabled" : "disabled" );
      version_str = strdup( buf );
    }
  }
  return version_str;
}

// -1 on error
int hdf5c_f_is_thread_safe( void ) {
  int ret_value = -1;
  hbool_t is_thread_safe = 0;
  // is this a thread-safe library build?
  herr_t h5_error = H5is_library_threadsafe( &is_thread_safe );
  if ( h5_error >= 0 ) { // ok
    ret_value = is_thread_safe ? 1 : 0;
  } else {
    ret_value = -1;
  }
  return ret_value;
}

// returns -1 if error
int checkHDF5threadSafe( void ) {
    unsigned major_number, minor_number, patch_number;
    static hbool_t is_thread_safe = 0;
    static int checked = 0;

    if ( !checked ) {
      // retrieve the library version
      herr_t h5_error = H5get_libversion( &major_number, &minor_number, &patch_number );
      if ( h5_error < 0 ) {
	    return -1;
      }
      // is this a thread-safe library build?
      h5_error = H5is_library_threadsafe( &is_thread_safe );
      if ( h5_error < 0 ) {
	    return -1;
      }

      printf( "Welcome to HDF5 %d.%d.%d  ( Thread-safety %s )\n", major_number, minor_number, patch_number, ( is_thread_safe > 0 ) ? "enabled" : "disabled" );
      checked = 1;
    }
    return ( is_thread_safe > 0 ) ? 0 : -1;
}

CPostHdf5File *new_CPostHdf5File( void ) {
  // thrread safety checking is done when opening second hdf5 file in gidpostHDF5.c new_CurrentHdf5WriteData()
  // int fail = checkHDF5threadSafe();
  // if ( fail ) {
  //   return NULL;
  // }
   CPostHdf5File *obj = ( CPostHdf5File * )malloc( sizeof( CPostHdf5File ));
   if ( obj ) {
     memset( obj, 0, sizeof( CPostHdf5File ) );
   }
   // obj->file_id = 0;
   // obj->file_id = 0;
   // obj->group_path[ 0] = '\0';
   // for ( int i = 0; i < MAX_CONCURRENT_DATASETS; i++ ) {
   //   memset( &( obj->file_buffer[ i ] ), 0, sizeof( struct hdf5c_buffer ) );
   // }
   // obj->error_buffer[ 0] = '\0';
   return obj;
}


void delete_CPostHdf5File( CPostHdf5File *obj ) {
  free( obj );
}

#define H5_GID_COMPRESSION_LEVEL   1
/*
  // with the H5Cpp interface:
  // Using level 1 because:
  // level 6 too costly, is 4 times slower than level 1 and only gets 6,6 % reduction
  // level 2 is 12 % slower and only 1,9 % size reduction
  #define H5_GID_COMPRESSION_LEVEL   6
*/

int hdf5c_f_open( CPostHdf5File *obj, const char *filename_in, OpenType otype ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i, is_hdf5;
  hid_t file_id;
  unsigned int flags = H5F_ACC_RDWR;

  H5Eset_auto( H5E_DEFAULT, NULL, NULL );

  is_hdf5 = H5Fis_hdf5( filename_in );
  if ( is_hdf5 <= 0 && otype == Read_OT ) {
    sprintf( obj->error_buffer, "File '%s' does not exists or is not in hdf5 format", filename_in );
    print_error( obj->error_buffer );
    return -1;
  }
  if ( is_hdf5 && otype != Create_OT ) {
    if ( otype == Read_OT ) flags = H5F_ACC_RDONLY;
    file_id = H5Fopen( filename_in, flags, H5P_DEFAULT );
    if ( file_id < 0 ) {
      sprintf( obj->error_buffer, "Could not open file '%s'", filename_in );
      print_error( obj->error_buffer );
      return -1;
    }
  } else { // Create_OT
#if H5_VERSION_GE( 1, 10, 0)
    // https://docs.hdfgroup.org/archive/support/HDF5/Tutor/swmr.html
    hid_t fapl = H5Pcreate( H5P_FILE_ACCESS );
    herr_t status = H5Pset_libver_bounds( fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST );
    if ( status < 0 ) {
      sprintf( obj->error_buffer, "could not create file '%s' (set_libver_bounds failed) ", filename_in );
      print_error( obj->error_buffer );
      return -1;
    }
    // it is one of opening with flag H5F_ACC_SWMR_WRITE or using H5Fstart_swmr_write() once after H5Dcreate().
    file_id = H5Fcreate( filename_in, H5F_ACC_TRUNC | H5F_ACC_SWMR_WRITE, H5P_DEFAULT, fapl );
    // file_id = H5Fcreate( filename_in, H5F_ACC_TRUNC , H5P_DEFAULT, fapl );
    // H5Pclose( fapl);
#else // H5_VERSION_GE( 1, 10, 0)
    file_id = H5Fcreate( filename_in, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
#endif // H5_VERSION_GE( 1, 10, 0)
    if ( file_id < 0 ) {
      sprintf( obj->error_buffer, "could not create file '%s'", filename_in );
      print_error( obj->error_buffer );
      return -1;
    }
  }
  obj->file_id = file_id;
  for ( i = 0; i < MAX_CONCURRENT_DATASETS; i++ ) {
    obj->file_buffer[ i ].used = 0;
    obj->file_buffer[ i ].file_id = file_id;
    memset( obj->file_buffer[ i ].intarray, 0, 20 * sizeof( int * ) );
    memset( obj->file_buffer[ i ].doublearray, 0, 20 * sizeof( double * ) );
  }
  strcpy( obj->error_buffer, "" );
  return 0;
}

int hdf5c_f_init( CPostHdf5File *obj, const char *filename ) {
  /* this is equal to hdf5c_open(filename_in,Create_OT) for compatibility */
  return hdf5c_f_open( obj, filename, Create_OT);
}

int hdf5c_f_end( CPostHdf5File *obj ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i, j;
  herr_t ht_ret = 0;
  if ( obj->group_id > 0 ) {
    H5Gclose( obj->group_id );
    obj->group_id = 0;
  }
  if ( obj->file_buffer[ 0 ].file_id > 0 ) {
    ht_ret = H5Fclose( obj->file_buffer[ 0 ].file_id );
    if ( ht_ret < 0) {
      sprintf( obj->error_buffer, "Error closing hdf5 file");
      print_error( obj->error_buffer);
    }
  }
  for ( i = 0; i < MAX_CONCURRENT_DATASETS; i++ ) {
    for ( j = 0; j < 20; j++ ) {
      if ( obj->file_buffer[ i ].intarray[ j ] ) {
        free( obj->file_buffer[ i ].intarray[ j ] );
        obj->file_buffer[ i ].intarray[ j ] = NULL;
      }
      if ( obj->file_buffer[ i ].doublearray[ j ] ) {
        free( obj->file_buffer[ i ].doublearray[ j ] );
        obj->file_buffer[ i ].doublearray[ j ] = NULL;
      }
    }
  }

  return ht_ret;
}

int hdf5c_f_flush( CPostHdf5File *obj ) {
  if ( obj == NULL ) {
    return -1;
  }
  herr_t status = H5Fflush( obj->file_id, H5F_SCOPE_GLOBAL );
  if ( status < 0 ) {
    print_error( "Can not flush hdf5 file" );
  }
  return status;
}

int hdf5c_f_flush_dataset( CPostHdf5File *obj, hid_t dataset_id, const char *name, const char *str_type ) {
  // not using flush data set, it seems that ?corrupts' hdf5's tables and program crashes at the end
  return 0;
  if ( obj == NULL ) {
    return -1;
  }
  herr_t status = H5Fflush(  obj->file_buffer[ dataset_id ].file_id, H5F_SCOPE_GLOBAL );
  if ( status < 0 ) {
    char buf[ 1024];
    sprintf( buf, "Can not flush hdf5 dataset/array name '%s' of type '%s'", name, str_type );
    print_error( buf );
  }
  return status;
}

char *hdf5c_f_last_error( CPostHdf5File *obj ) {
  return obj ? obj->error_buffer : NULL;
}

/*################################################################################
 *    write functions
 *################################################################################*/

int hdf5c_f_create_group( CPostHdf5File *obj, const char *name ) {
  if ( obj == NULL ) {
    return -1;
  }

  H5G_info_t g_info;
  herr_t herr = H5Gget_info_by_name( obj->file_id, name, &g_info, H5P_DEFAULT);
  if ( herr >=0) // it exists
    return 0;

  hid_t ht_group = 0;
  // ht_group=H5Gopen(obj->file_id,name,H5P_DEFAULT);
  // if(ht_group>=0){
  //   H5Gclose(ht_group);
  //   return 0;
  // }
  ht_group=H5Gcreate(obj->file_id,name,0,H5P_DEFAULT,H5P_DEFAULT);
  if(ht_group<0){
    sprintf(obj->error_buffer,"name '%s'. error creating group",name);
    print_error( obj->error_buffer);
    return -1;
  } else {
    H5Gclose (ht_group);
  }
  return 1;
}

int hdf5c_f_set_attribute( CPostHdf5File *obj, const char *name, const char *attname, const char *value ) {
  if ( obj == NULL ) {
    return -1;
  }

  herr_t ht_err = 0;

  /*
    ht_err=H5LTset_attribute_char(obj->file_buffer.file_id,name,attname,value,strlen(value));
  */
  if ( value) {
    ht_err = H5LTset_attribute_string( obj->file_id, name, attname, value );
  } else {
    ht_err = -1;
  }
  if ( ht_err < 0 ) {
    sprintf( obj->error_buffer, "cannot set attribute '%s' to '%s'", name, value ? value : "(NULL)" );
    print_error( obj->error_buffer);
    return -1;
  }
  return 0;
}

int hdf5c_f_start_dataset( CPostHdf5File *obj, const char *name, int num_int, int num_real ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i,dataset_id;
  
  for(dataset_id=0;dataset_id<MAX_CONCURRENT_DATASETS;dataset_id++)
    if(!obj->file_buffer[dataset_id].used) break;
  if(dataset_id==MAX_CONCURRENT_DATASETS){
    sprintf(obj->error_buffer,"Too many concurrent datasets n=%d",MAX_CONCURRENT_DATASETS);
    print_error( obj->error_buffer);
    return -1;
  }
  
  if(strlen(name)>1023){
    sprintf(obj->error_buffer,"name '%s'. too long",name);
    print_error( obj->error_buffer);
    return -1;
  }
  strcpy(obj->file_buffer[dataset_id].datasetname,name);
  // hid_t ht_group = 0;
  // ht_group=H5Gopen(obj->file_buffer[dataset_id].file_id,name,H5P_DEFAULT);
  H5G_info_t g_info;
  herr_t herr = H5Gget_info_by_name( obj->file_buffer[dataset_id].file_id, name, &g_info, H5P_DEFAULT);
  if ( herr < 0) {
    hid_t ht_group=H5Gcreate(obj->file_buffer[dataset_id].file_id,name,0,H5P_DEFAULT,H5P_DEFAULT);
    if(ht_group<0){
      sprintf(obj->error_buffer,"name '%s'. error creating group",name);
      print_error( obj->error_buffer);
      return -1;
    } else {
      H5Gclose (ht_group);
    }
  }
  
  obj->file_buffer[dataset_id].used=1;
  if(obj->file_buffer[dataset_id].memsize==0){
    obj->file_buffer[dataset_id].memsize=10000;
  }
  for(i=obj->file_buffer[dataset_id].num_int;i<num_int;i++){
    if(!obj->file_buffer[dataset_id].intarray[i])
      obj->file_buffer[dataset_id].intarray[i]=(int*) malloc(( size_t)obj->file_buffer[dataset_id].memsize*sizeof(int));
  }
  obj->file_buffer[dataset_id].num_int=num_int;
  for(i=obj->file_buffer[dataset_id].num_real;i<num_real;i++){
    if(!obj->file_buffer[dataset_id].doublearray[i])
      obj->file_buffer[dataset_id].doublearray[i]=(double*) malloc(( size_t)obj->file_buffer[dataset_id].memsize*sizeof(double));
  }
  obj->file_buffer[dataset_id].num_real=num_real;
  obj->file_buffer[dataset_id].size=0;
  return dataset_id;
}

int hdf5c_f_add_to_dataset( CPostHdf5File *obj, int dataset_id, const int intvalues[],
                            const double doublevalues[] ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i;
  if ( obj->file_buffer[ dataset_id ].size == obj->file_buffer[ dataset_id ].memsize ) {
    obj->file_buffer[ dataset_id ].memsize *= 2;
    for ( i = 0; i < 20; i++ ) {
      if ( obj->file_buffer[ dataset_id ].intarray[ i ] ) {
        obj->file_buffer[ dataset_id ].intarray[ i ] =
            realloc( obj->file_buffer[ dataset_id ].intarray[ i ], ( size_t )obj->file_buffer[ dataset_id ].memsize * sizeof( int ) );
      }
      if ( obj->file_buffer[ dataset_id ].doublearray[ i ] ) {
        obj->file_buffer[ dataset_id ].doublearray[ i ] =
            realloc( obj->file_buffer[ dataset_id ].doublearray[ i ], ( size_t )obj->file_buffer[ dataset_id ].memsize * sizeof( double ) );
      }
    }
  }
  for ( i = 0; i < obj->file_buffer[ dataset_id ].num_int; i++ ) {
    obj->file_buffer[ dataset_id ].intarray[ i ][ obj->file_buffer[ dataset_id ].size ] = intvalues[ i ];
  }
  for ( i = 0; i < obj->file_buffer[ dataset_id ].num_real; i++ ) {
    obj->file_buffer[ dataset_id ].doublearray[ i ][ obj->file_buffer[ dataset_id ].size ] = doublevalues[ i ];
  }
  obj->file_buffer[ dataset_id ].size++;
  return 0;
}

static int hdf5c_f_create_int_dataset_from_buffer( CPostHdf5File *obj, int dataset_id, char *name,
                                                   int compress, int length, int *intarray ) {
  if ( obj == NULL ) {
    return -1;
  }
  if ( length <= 0) {
    char *buf = ( char *)malloc( ( strlen( name) + 1024) * sizeof( char));
    sprintf( buf, "create_int_dataset: Empty array, data set not created for '%s'", name);
    print_error( buf);
    free( buf);
    return 0; // -1;
  }

  herr_t ht_err = 0;
  hid_t plist, sid, dataset;
  hsize_t dims[ 1 ];
  dims[ 0 ] = ( hsize_t )length;

  if ( ( compress > 0) && ( length > 0) ) {
    hsize_t chunk_size[ 1 ];
    chunk_size[ 0 ] = ( hsize_t )length;
    plist = H5Pcreate( H5P_DATASET_CREATE );
    ht_err = H5Pset_chunk( plist, 1, chunk_size );
    if ( ht_err < 0 ) {
      sprintf( obj->error_buffer, "cannot H5Pset_chunk size to %d for array name '%s'",
               ( int )chunk_size[ 0 ], name );
      print_error( obj->error_buffer );
    }
    H5Pset_deflate( plist, ( unsigned int )compress );
    sid = H5Screate_simple( 1, dims, NULL );
    dataset = H5Dcreate( obj->file_buffer[ dataset_id ].file_id, name, H5T_NATIVE_INT, sid,
                         H5P_DEFAULT, plist, H5P_DEFAULT );
    // ht_err = H5Fstart_swmr_write( obj->file_buffer[ dataset_id ].file_id );
    ht_err = H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intarray );
    if ( ht_err < 0 ) {
      sprintf( obj->error_buffer, "error writing int dataset/array name '%s'", name );
      print_error( obj->error_buffer );
    }
    // hid_t sid_write = H5Screate_simple( 1, dims, NULL );
    // ht_err = H5Dwrite( dataset, H5T_NATIVE_INT, sid_write, sid, H5P_DEFAULT, intarray );
    // H5Oenable_mdc_flushes( dataset );
    // H5Sclose( sid_write );
    H5Sclose( sid );
    H5Pclose( plist );
    H5Dclose( dataset );
    // hdf5c_f_flush_dataset( obj, dataset_id, name, "H5T_NATIVE_INT");
  } else {
    // no compression or empty dataset
    ht_err =
        H5LTmake_dataset_int( obj->file_buffer[ dataset_id ].file_id, name, 1, dims, intarray );
  }
  if ( ht_err < 0 ) {
    sprintf( obj->error_buffer, "cannot set array name '%s'", name );
    print_error( obj->error_buffer );
    return -1;
  }
  return 0;
}

static int hdf5c_f_create_double_dataset_from_buffer( CPostHdf5File *obj, int dataset_id,
                                                      char *name, int compress, int length,
                                                      double *doublearray ) {
  if ( obj == NULL ) {
    return -1;
  }
  if ( length <= 0) {
    char *buf = ( char *)malloc( ( strlen( name) + 1024) * sizeof( char));
    sprintf( buf, "create_double_dataset: Empty array, data set not created for '%s'", name);
    print_error( buf);
    free( buf);
    return 0; // -1;
  }

  herr_t ht_err = 0;
  hid_t plist, sid, dataset;
  hsize_t dims[ 1 ];
  dims[ 0 ] = ( hsize_t )length;

  if ( ( compress > 0) && ( length > 0) ) {
    hsize_t chunk_size[ 1 ];
    chunk_size[ 0 ] = ( hsize_t )length;
    plist = H5Pcreate( H5P_DATASET_CREATE );
    ht_err = H5Pset_chunk( plist, 1, chunk_size );
    if ( ht_err < 0 ) {
      sprintf( obj->error_buffer, "cannot H5Pset_chunk size to %d for array name '%s'",
               ( int )chunk_size[ 0 ], name );
      print_error( obj->error_buffer );
    }
    H5Pset_deflate( plist, ( unsigned int )compress );
    sid = H5Screate_simple( 1, dims, NULL );
    dataset = H5Dcreate( obj->file_buffer[ dataset_id ].file_id, name, H5T_NATIVE_DOUBLE, sid,
                         H5P_DEFAULT, plist, H5P_DEFAULT );
    // ht_err = H5Fstart_swmr_write( obj->file_buffer[ dataset_id ].file_id );
    ht_err = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, doublearray );
      if ( ht_err < 0 ) {
        sprintf( obj->error_buffer, "error writing double dataset/array name '%s'", name );
        print_error( obj->error_buffer );
      }
    // hid_t sid_write = H5Screate_simple( 1, dims, NULL );
    // ht_err = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, sid_write, sid, H5P_DEFAULT, doublearray );
    // H5Oenable_mdc_flushes( dataset );
    // H5Sclose( sid_write );
    H5Sclose( sid );
    H5Pclose( plist );
    H5Dclose( dataset );
    // hdf5c_f_flush_dataset( obj, dataset_id, name, "H5T_NATIVE_DOUBLE");
  } else {
    // no compression or empty dataset
    ht_err = H5LTmake_dataset_double( obj->file_buffer[ dataset_id ].file_id, name, 1, dims,
                                      doublearray);
  }
  if ( ht_err < 0 ) {
    sprintf( obj->error_buffer, "cannot set array name '%s'", name );
    print_error( obj->error_buffer );
    return -1;
  }
  return 0;
}

int hdf5c_f_end_dataset( CPostHdf5File *obj, int dataset_id ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i, fail, inum = 1;
  int compress = H5_GID_COMPRESSION_LEVEL;
  char buf[ 2048 ];
  for ( i = 0; i < obj->file_buffer[ dataset_id ].num_int; i++ ) {
    sprintf( buf, "%s/%d", obj->file_buffer[ dataset_id ].datasetname, inum );
    fail = hdf5c_f_create_int_dataset_from_buffer( obj, dataset_id, buf, compress,
                                                   obj->file_buffer[ dataset_id ].size,
                                                   obj->file_buffer[ dataset_id ].intarray[ i ] );
    if ( fail < 0 )
      return fail;
    inum++;
  }
  for ( i = 0; i < obj->file_buffer[ dataset_id ].num_real; i++ ) {
    sprintf( buf, "%s/%d", obj->file_buffer[ dataset_id ].datasetname, inum );
    fail = hdf5c_f_create_double_dataset_from_buffer(
        obj, dataset_id, buf, compress, obj->file_buffer[ dataset_id ].size,
        obj->file_buffer[ dataset_id ].doublearray[ i ] );
    if ( fail < 0 )
      return fail;
    inum++;
  }
  obj->file_buffer[ dataset_id ].used = 0;
  return 0;
}

const char *getDatasetTypeString( const t_hdf5c_dataset_type ds_type) {
    const char *ret_str = "Unknown type";
    switch ( ds_type ) {
    case HDF5C_DATASET_INTEGER:
      ret_str = "HDF5C_DATASET_INTEGER";
      break;
    case HDF5C_DATASET_DOUBLE:
      ret_str = "HDF5C_DATASET_DOUBLE";
      break;
    case HDF5C_DATASET_NONE:
      ret_str = "HDF5C_DATASET_NONE";
      break;
    }
    return ret_str;
}

static int _hdf5c_f_create_dataset( CPostHdf5File *obj, const char *dataset_name,
                                    int dataset_subindex, int compress, int length,
                                    t_hdf5c_dataset_type ds_type, const void *data_array ) {
  if ( obj == NULL ) {
    return -1;
  }
  if ( length <= 0) {
    char buf[ 1024];
    sprintf( buf, "create_dataset: Empty array, data set not created for '%s'", dataset_name);
    print_error( buf);
    return 0; // -1;
  }

  char name_dataset_subtable[ 2048 ];
  snprintf( name_dataset_subtable, 2048, "%s/%d", dataset_name, dataset_subindex );

  hsize_t dims[ 1 ];
  herr_t ht_err = 0;
  int fail = 0;
  
  dims[ 0 ] = ( hsize_t )length;
  if (( compress > 0) && ( length > 0)) {
    hid_t dataset_id = 0;
    hid_t type_id = H5T_NATIVE_INT;
    switch ( ds_type ) {
    case HDF5C_DATASET_INTEGER:
      type_id = H5T_NATIVE_INT;
      break;
    case HDF5C_DATASET_DOUBLE:
      type_id = H5T_NATIVE_DOUBLE;
      break;
    case HDF5C_DATASET_NONE:
      fail = -1;
      ht_err = -1;
      break;
    }
    if ( !fail ) {
      hsize_t chunk_size[ 1 ];
      chunk_size[ 0 ] = ( hsize_t )length;
      hid_t plist = H5Pcreate( H5P_DATASET_CREATE );
      ht_err = H5Pset_chunk( plist, 1, chunk_size );
      if ( ht_err < 0 ) {
        sprintf( obj->error_buffer, "cannot H5Pset_chunk size to %d for array name '%s'",
                 ( int )chunk_size[ 0 ], name_dataset_subtable );
        print_error( obj->error_buffer );
      }
      H5Pset_deflate( plist, ( unsigned int )compress );
      hid_t sid = H5Screate_simple( 1, dims, NULL );
      dataset_id = H5Dcreate( obj->file_id, name_dataset_subtable, type_id, sid, H5P_DEFAULT,
                              plist, H5P_DEFAULT );
      // ht_err = H5Fstart_swmr_write( obj->file_buffer[ dataset_id ].file_id );
      ht_err = H5Dwrite( dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_array );
      if ( ht_err < 0 ) {
        sprintf( obj->error_buffer, "error writing '%s' dataset/array name '%s'",
                 getDatasetTypeString( ds_type), name_dataset_subtable );
        print_error( obj->error_buffer );
      }
      // hid_t sid_write = H5Screate_simple( 1, dims, NULL );
      // ht_err = H5Dwrite( dataset_id, type_id, sid_write, sid, H5P_DEFAULT, data_array );
      // H5Oenable_mdc_flushes( dataset_id );
      // H5Sclose( sid_write );
      H5Sclose( sid );
      H5Pclose( plist );
      H5Dclose( dataset_id );
      // hdf5c_f_flush_dataset( obj, dataset_id, name_dataset_subtable, getDatasetTypeString( type_id));
    }
  } else {
    // no compression or empty dataset
    switch ( ds_type ) {
    case HDF5C_DATASET_INTEGER:
      ht_err = H5LTmake_dataset_int( obj->file_id, name_dataset_subtable, 1, dims,
                                     ( const int * )data_array );
      break;
    case HDF5C_DATASET_DOUBLE:
      ht_err = H5LTmake_dataset_double( obj->file_id, name_dataset_subtable, 1, dims,
                                        ( const double * )data_array );
      break;
    case HDF5C_DATASET_NONE:
      fail = -1;
      ht_err = -1;
      break;
    }
  }
  if ( ht_err < 0 ) {
    sprintf( obj->error_buffer, "cannot set array name '%s' with type '%s'", name_dataset_subtable, getDatasetTypeString( ds_type) );
    print_error( obj->error_buffer );
    fail = -1;
  } else {
    fail = 0;
  }
  return fail;
}

// hdf5c_write_dataset_list = hdf5c_f_start_dataset() + hdf5c_f_add_to_dataset() +
// hdf5c_f_end_dataset()
int hdf5c_write_dataset_list( CPostHdf5File *obj, const char *name, int num_datasets, int num_values, 
			      t_hdf5c_dataset_type *data_types, const void **data_values ) {
  if ( obj == NULL ) {
    return -1;
  }

  if ( strlen( name ) > 1023 ) {
    sprintf( obj->error_buffer, "name '%s'. too long", name );
    print_error( obj->error_buffer);
    return -1;
  }

  // check if name exists, if not, create it.
  // hid_t ht_group = 0;
  // ht_group = H5Gopen( obj->file_id, name, H5P_DEFAULT );
  H5G_info_t g_info;
  herr_t herr = H5Gget_info_by_name( obj->file_id, name, &g_info, H5P_DEFAULT);
  if ( herr < 0) {
    hid_t ht_group = H5Gcreate( obj->file_id, name, 0, H5P_DEFAULT, H5P_DEFAULT );
    if ( ht_group < 0 ) {
      sprintf( obj->error_buffer, "name '%s'. error creating group", name );
      print_error( obj->error_buffer);
      return -1;
    } else {
      H5Gclose( ht_group );
    }
  }

  // write datasets: ids, x, y, z
  int compress = H5_GID_COMPRESSION_LEVEL;
  int dataset_subindex = 1; 
  int fail = 0;
  for ( int idx_data = 0; ( idx_data < num_datasets) && ( data_types[ idx_data] != HDF5C_DATASET_NONE) && data_values[ idx_data]; idx_data++ ) {
    t_hdf5c_dataset_type current_type = data_types[ idx_data ];
    const void *current_data = data_values[ idx_data ];
    fail = _hdf5c_f_create_dataset( obj, name, dataset_subindex, compress, num_values, current_type, current_data );
    if ( fail ) {
      break;
    }
    dataset_subindex++;
  }

  return fail;
}

/*################################################################################
 *    Read functions
 *################################################################################*/

int hdf5c_f_num_children( CPostHdf5File *obj, const char *group_path ) {
  if ( obj == NULL ) {
    return -1;
  }

  herr_t ht_err = 0;
  hsize_t num_obj;

  if( obj->group_id > 0) {
    H5Gclose( obj->group_id);
    obj->group_id = 0;
  }
  obj->group_id=H5Gopen(obj->file_id,group_path,H5P_DEFAULT);
  if(obj->group_id<0){
    sprintf(obj->error_buffer,"Could not find group '%s'",group_path);
    print_error( obj->error_buffer);
    return -1;
  }
  strcpy(obj->group_path,group_path);
  
  ht_err=H5Gget_num_objs(obj->group_id,&num_obj);

  if(ht_err<0){
    sprintf(obj->error_buffer,"group '%s' is not a correct group",group_path);
    print_error( obj->error_buffer);
    return -1;
  }
  return (int)num_obj;
}

int hdf5c_f_get_children_name( CPostHdf5File *obj, const char *group_path, int index, char *name, int max_len ) {
  if ( obj == NULL ) {
    return -1;
  }

  ssize_t len,len0;
  if(obj->group_id<=0 || strcmp(group_path,obj->group_path)!=0){
    if ( obj->group_id > 0) {
        H5Gclose(obj->group_id);
        obj->group_id = 0;
    }
    obj->group_id=H5Gopen(obj->file_id,group_path,H5P_DEFAULT);
    if(obj->group_id<0){
      sprintf(obj->error_buffer,"Could not find group '%s'",group_path);
      print_error( obj->error_buffer);
      return -1;
    }
    strcpy(obj->group_path,group_path);
  }
  if((int)strlen(group_path)>=max_len-10) return -1;
  strcpy(name,group_path);
  if(name[strlen(name)-1]!='/'){
    strcat(name,"/");
  }
  len0= ( ssize_t)strlen(name);
  len=H5Gget_objname_by_idx(obj->group_id, ( hsize_t)index,&name[len0],( size_t)max_len-( size_t)len0);
  return ( int)( len+len0);
}

int hdf5c_f_get_attribute( CPostHdf5File *obj, const char *obj_name, char *att_name, char *value, int max_len ) {
  if ( obj == NULL ) {
    return -1;
  }

  int rank;
  herr_t ht_err = 0;
  hsize_t dims[10];
  size_t  type_size;
  H5T_class_t class_id;
  
  ht_err=H5LTget_attribute_info(obj->file_id,obj_name,att_name, dims,&class_id,&type_size);
  if(ht_err<0){
    sprintf(obj->error_buffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
    print_error( obj->error_buffer);
    return -1;
  }  
  if(class_id==H5T_INTEGER) {
    if(type_size!=1){
      sprintf( obj->error_buffer,
               "name '%s' cannot be retrieved. Bad type %d. Has to be one array of char",
	           obj_name,class_id);
      print_error( obj->error_buffer);
      return -1;
    }
    ht_err=H5LTget_attribute_ndims(obj->file_id,obj_name,att_name,&rank);
    /*
      rank is 1 for H5T_INTEGER
    */
    if(ht_err<0 || rank!=1){
      sprintf(obj->error_buffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
      print_error( obj->error_buffer);
      return -1;
    }
    if(dims[0]>= ( hsize_t)max_len) return -1;
    ht_err=H5LTget_attribute_char(obj->file_id,obj_name,att_name,value);
  } else if (class_id==H5T_STRING){
    if(type_size<1){
      sprintf( obj->error_buffer,
               "name '%s' cannot be retrieved. Bad type %d. Has to be one array of char",
	           obj_name,class_id);
      print_error( obj->error_buffer);
      return -1;
    }
    ht_err=H5LTget_attribute_ndims(obj->file_id,obj_name,att_name,&rank);
    /*
      rank is 0 for H5T_STRING
    */
    if(ht_err<0 || rank!=0){
      sprintf(obj->error_buffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
      print_error( obj->error_buffer);
      return -1;
    }
    if((int)type_size>=max_len) return -1;
    ht_err=H5LTget_attribute_string(obj->file_id,obj_name,att_name,value);
  } else {
    sprintf( obj->error_buffer,
             "name '%s' cannot be retrieved. Bad type %d Has to be one array of char or a string",
             obj_name,class_id);
    print_error( obj->error_buffer);
    return -1;
  }
  return 0;
}

int hdf5c_f_open_dataset_list( CPostHdf5File *obj, const char *obj_name, int *num_int, int *num_real ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i,rank,idx;
  hsize_t num;
  char name[1024];
  
  herr_t ht_err = 0;
  hid_t ht_dataset,ht_type;
  hsize_t dims[1];
  H5T_class_t class_id;
  size_t type_size;

  num=( hsize_t)hdf5c_f_num_children( obj, obj_name);
  if(num<=0) return -1;
  
  obj->file_buffer[0].num_int=obj->file_buffer[0].num_real=0;
  
  for(i=0;i< ( int)num;i++){
    hdf5c_f_get_children_name( obj, obj_name,i,name,1024);
    ht_err=H5LTget_dataset_ndims(obj->file_id,name,&rank);
    if(ht_err<0 || rank!=1){
      sprintf(obj->error_buffer,"dataset name '%s' cannot be retrieved",name);
      print_error( obj->error_buffer);
      return -1;
    }
    ht_err=H5LTget_dataset_info (obj->file_id,name,dims,&class_id,NULL);
    if(ht_err<0){
      sprintf(obj->error_buffer,"dataset name '%s' cannot be retrieved",name);
      print_error( obj->error_buffer);
      return -1;
    }
    ht_dataset=H5Dopen(obj->file_id,name,H5P_DEFAULT);
    ht_type=H5Dget_type(ht_dataset);
    type_size=H5Tget_size(ht_type);
    H5Dclose(ht_dataset);
    H5Tclose(ht_type);
    
    if(class_id==H5T_INTEGER){
      idx=obj->file_buffer[0].num_int;
      obj->file_buffer[0].intarray[idx]=malloc(((size_t)dims[0])*sizeof(int));
      ht_err=H5LTread_dataset_int (obj->file_id,name,obj->file_buffer[0].intarray[idx]);
      if(ht_err<0){
	    sprintf(obj->error_buffer,"array name '%s' cannot be retrieved",name);
        print_error( obj->error_buffer);
	    return -1;
      }
      obj->file_buffer[0].num_int++;
    }
    else if(class_id==H5T_FLOAT && type_size==4){
      sprintf(obj->error_buffer,"4 byte floats not supported");
      print_error( obj->error_buffer);
      return -1;
    }
    else if(class_id==H5T_FLOAT && type_size==8){
      idx=obj->file_buffer[0].num_real;
      obj->file_buffer[0].doublearray[idx]=malloc(((size_t)dims[0])*sizeof(double));
      ht_err=H5LTread_dataset_double(obj->file_id,name,obj->file_buffer[0].doublearray[idx]);
      if(ht_err<0){
	    sprintf(obj->error_buffer,"array name '%s' cannot be retrieved",name);
        print_error( obj->error_buffer);
	    return -1;
      }
      obj->file_buffer[0].num_real++;
    }
    else {
      sprintf( obj->error_buffer,
               "array name '%s' cannot be retrieved. unknown type=%d,%llu",
	           name,class_id,( long long unsigned int)type_size);
      print_error( obj->error_buffer);
      return -1;
    }
    obj->file_buffer[0].size=(int)dims[0];
  }
  obj->file_buffer[0].used=1;
  *num_int=obj->file_buffer[0].num_int;
  *num_real=obj->file_buffer[0].num_real;
  return obj->file_buffer[0].size;
}

int hdf5c_f_give_dataset_list_values( CPostHdf5File *obj, int index, int int_values[], double double_values[] ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i;
  if(index<0 || index>=obj->file_buffer[0].size){
    return -1;
  }
  for(i=0;i<obj->file_buffer[0].num_int;i++){
    int_values[i]=obj->file_buffer[0].intarray[i][index];
  }
  for(i=0;i<obj->file_buffer[0].num_real;i++){
    double_values[i]=obj->file_buffer[0].doublearray[i][index];
  }
  return 0;
}

int hdf5c_f_close_dataset_list( CPostHdf5File *obj ) {
  if ( obj == NULL ) {
    return -1;
  }

  int i;
  for ( i = 0; i < obj->file_buffer[ 0 ].num_int; i++ ) {
    free( obj->file_buffer[ 0 ].intarray[ i ] );
    obj->file_buffer[ 0 ].intarray[ i ] = NULL;
  }
  for ( i = 0; i < obj->file_buffer[ 0 ].num_real; i++ ) {
    free( obj->file_buffer[ 0 ].doublearray[ i ] );
    obj->file_buffer[ 0 ].doublearray[ i ] = NULL;
  }
  obj->file_buffer[ 0 ].num_int = obj->file_buffer[ 0 ].num_real = 0;
  obj->file_buffer[ 0 ].used = 0;
  return 0;
}
