#ifndef __HDF5C__
#define __HDF5C__

#ifdef __cplusplus
extern "C" {
#endif

  /* used in hdf5c.c and gidpostHDF5.c */
#define MAX_CONCURRENT_DATASETS 20

  // only used in hdf5c_f_open() which is never called from outside hdf5c.c
  // typedef enum _OpenType {
  //   Read_OT,
  //   Write_OT,
  //   Create_OT
  // } OpenType;

  typedef enum _t_hdf5c_dataset_type {
    HDF5C_DATASET_NONE = 0,
    HDF5C_DATASET_INTEGER,
    HDF5C_DATASET_DOUBLE,
  } t_hdf5c_dataset_type;
  const char *getDatasetTypeString( const t_hdf5c_dataset_type type_id);

  typedef struct _CPostHdf5File CPostHdf5File;

  CPostHdf5File *new_CPostHdf5File( void );
  void delete_CPostHdf5File( CPostHdf5File *obj );

  const char *hdf5c_f_get_version( void );
  int hdf5c_f_is_thread_safe( void ); // returns -1 on error

  /*################################################################################
   *    file functions
   ################################################################################*/

  /* returns -1 or NULL on error */
  /* multiple file version */
  // hdf5c_f_open() is never called from outside hdf5c.c
  // int hdf5c_f_open( CPostHdf5File *obj, const char *zFile, OpenType otype );
  int hdf5c_f_init( CPostHdf5File *obj, const char *zFile );  /* this is equal to hdf5c_open(zFile,Create_OT) for compatibility */
  int hdf5c_f_end( CPostHdf5File *obj );
  int hdf5c_f_flush( CPostHdf5File *obj );
  char *hdf5c_f_last_error( CPostHdf5File *obj );

  /*################################################################################
   *    write functions
   ################################################################################*/

  /* returns -1 or NULL on error */
  /* multiple file version */
  int hdf5c_f_create_group( CPostHdf5File *obj, const char *name );
  int hdf5c_f_set_attribute( CPostHdf5File *obj, const char *name, const char *attname, const char *value );
  int hdf5c_f_start_dataset( CPostHdf5File *obj, const char *name, int num_int, int num_real );
  int hdf5c_f_add_to_dataset( CPostHdf5File *obj, int dataset_id, const int intvalues[], const double doublevalues[] );
  int hdf5c_f_end_dataset( CPostHdf5File *obj, int dataset_id );

  // hdf5c_write_dataset_list = hdf5c_f_start_dataset() + hdf5c_f_add_to_dataset() + hdf5c_f_end_dataset()
  int hdf5c_write_dataset_list( CPostHdf5File *obj, const char *name, int num_datasets, int num_values, 
                                t_hdf5c_dataset_type *data_types, const void **data_values );

  /*################################################################################
   *    Read functions
   ################################################################################*/

  /* returns -1 or NULL on error */
  /* multiple file version */
  int hdf5c_f_num_children( CPostHdf5File *obj, const char *group_path );
  int hdf5c_f_get_children_name( CPostHdf5File *obj, const char *group_path, int index, char *name, int max_len );
  int hdf5c_f_get_attribute( CPostHdf5File *obj, const char *obj_name, char *att_name, char *value, int max_len );
  int hdf5c_f_open_dataset_list( CPostHdf5File *obj, const char *obj_name, int *num_int, int *num_real );
  int hdf5c_f_give_dataset_list_values( CPostHdf5File *obj, int index, int int_values[], double double_values[] );
  int hdf5c_f_close_dataset_list( CPostHdf5File *obj);

#ifdef __cplusplus
}
#endif

#endif /* HDF5C */
