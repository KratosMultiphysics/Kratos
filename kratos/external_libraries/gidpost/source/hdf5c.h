  

#ifndef __HDF5C__
#define __HDF5C__

#ifdef __cplusplus
extern "C" {
#endif

#include "hdf5.h"
  
#define MAX_CONCURRENT_DATASETS 20
  
typedef enum _OpenType {
  Read_OT,
  Write_OT,
  Create_OT
} OpenType;

struct hdf5c_buffer
{
  int used;
  hid_t file_id;
  char datasetname[1024];
  int memsize,size,num_int,num_real;
  int* intarray[20];
  double* doublearray[20];
};
  
int hdf5c_open(const char* zFile,OpenType otype);
int hdf5c_init(const char* zFile);  /* this is equal to hdf5c_open(zFile,Create_OT) fro compatibility */
int hdf5c_end();
int hdf5c_flush();
char* hdf5c_last_error();

/*################################################################################
 *    write functions
 ################################################################################*/

int hdf5c_create_group(const char* name);
int hdf5c_set_attribute(const char* name,const char* attname,const char* value);

/* returns dataset_id; -1 on error */
int hdf5c_start_dataset(const char* name,int num_int,int num_real);
int hdf5c_addto_dataset(int dataset_id,int intvalues[],double doublevalues[]);
int hdf5c_end_dataset(int dataset_id);

/*################################################################################
 *    Read functions
 ################################################################################*/

int hdf5c_num_children(const char* group_path);
int hdf5c_get_children_name(const char* group_path,int index,char* name,int max_len);
int hdf5c_get_attribute(const char* obj_name,char* att_name,char* value,int max_len);
int hdf5c_open_dataset_list(const char* obj_name,int* num_int,int* num_real);
int hdf5c_give_dataset_list_values(int index,int int_values[],double double_values[]);
int hdf5c_close_dataset_list();
  
  
/*################################################################################
 *    Funtions that do not use global variables
 ################################################################################*/

hid_t hdf5cM_open(const char* zFile,OpenType otype);
int hdf5cM_end(hid_t file_id);
int hdf5cM_flush(hid_t file_id);
int hdf5cM_get_attribute(hid_t file_id,const char* obj_name,char* att_name,char* value,
  int max_len);
int hdf5cM_open_dataset_list(hid_t file_id,struct hdf5c_buffer* hdf5c_buffer,
  const char* obj_name,int* num_int,int* num_real);
int hdf5cM_close_dataset_list(struct hdf5c_buffer* hdf5c_buffer);
int hdf5cM_give_dataset_list_values(struct hdf5c_buffer* hdf5c_buffer,int index,
    int int_values[],double double_values[]);
int hdf5cM_num_children(hid_t file_id,const char* group_path);
int hdf5cM_get_children_name(hid_t file_id,const char* group_path,
    int index,char* name,int max_len);
  
  /* NOTE: the functions below are extracted from HDF5LT and are to be used internally only */

extern herr_t H5LTfind_dataset( hid_t loc_id, const char *dset_name );
extern herr_t  H5LTmake_dataset_char( hid_t loc_id,
		              const char *dset_name,
		              int rank,
		              const hsize_t *dims,
		              const char *buffer );

extern herr_t H5LTmake_dataset_int( hid_t loc_id,
		            const char *dset_name,
		            int rank,
		            const hsize_t *dims,
		            const int *data );
extern herr_t H5LTmake_dataset_double( hid_t loc_id,
		               const char *dset_name,
		               int rank,
		               const hsize_t *dims,
		               const double *data );
extern herr_t H5LTmake_dataset_float( hid_t loc_id,
		              const char *dset_name,
		              int rank,
		              const hsize_t *dims,
		              const float *data );
extern herr_t H5LTget_dataset_ndims( hid_t loc_id,
		             const char *dset_name,
		             int *rank );
extern herr_t H5LTget_dataset_info( hid_t loc_id,
		            const char *dset_name,
		            hsize_t *dims,
		            H5T_class_t *type_class,
		            size_t *type_size );
extern herr_t  H5LTread_dataset_char( hid_t loc_id,
		              const char *dset_name,
		              char *buffer );
extern herr_t H5LTread_dataset_int( hid_t loc_id,
		            const char *dset_name,
		            int *data );
extern herr_t H5LTread_dataset_float( hid_t loc_id,
		              const char *dset_name,
		              float *data );
extern herr_t H5LTread_dataset_double( hid_t loc_id,
  const char *dset_name,
  double *data );
extern herr_t H5LTset_attribute_string( hid_t loc_id,
		                const char *obj_name,
		                const char *attr_name,
		                const char *attr_data );
extern herr_t H5LTset_attribute_char( hid_t loc_id,
		              const char *obj_name,
		              const char *attr_name,
		              const char *data,
		              size_t size );
extern herr_t H5LTget_attribute_info( hid_t loc_id,
		              const char *obj_name,
		              const char *attr_name,
		              hsize_t *dims,
		              H5T_class_t *type_class,
		              size_t *type_size );
extern herr_t H5LTget_attribute_ndims( hid_t loc_id,
		               const char *obj_name,
		               const char *attr_name,
		               int *rank );
extern herr_t H5LTget_attribute_string( hid_t loc_id,
		                const char *obj_name,
		                const char *attr_name,
		                char *data );
extern herr_t H5LTget_attribute_char( hid_t loc_id,
		              const char *obj_name,
		              const char *attr_name,
		              char *data );
extern herr_t H5LTmake_dataset_int( hid_t loc_id,
		            const char *dset_name,
		            int rank,
		            const hsize_t *dims,
		            const int *data );
extern herr_t H5LTmake_dataset_double( hid_t loc_id,
		               const char *dset_name,
		               int rank,
		               const hsize_t *dims,
		               const double *data );
extern herr_t H5LTset_attribute_string( hid_t loc_id,
		                const char *obj_name,
		                const char *attr_name,
		                const char *attr_data );
extern herr_t H5LTmake_dataset_int( hid_t loc_id,
		            const char *dset_name,
		            int rank,
		            const hsize_t *dims,
		            const int *data );
extern herr_t H5LTmake_dataset_double( hid_t loc_id,
		               const char *dset_name,
		               int rank,
		               const hsize_t *dims,
		               const double *data );
extern herr_t H5LTfind_attribute( hid_t loc_id, const char* attr_name );
extern herr_t H5LT_find_attribute( hid_t loc_id, const char* attr_name );

#ifdef __cplusplus
}
#endif

#endif /* HDF5C */
