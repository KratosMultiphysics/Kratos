
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "hdf5c.h"

static hid_t global_file_id=0;
static hid_t global_group_id=0;
static char global_group_path[1024];
static struct hdf5c_buffer global_hdf5c_buffer[MAX_CONCURRENT_DATASETS];
static char errorBuffer[1024];

int hdf5c_open(const char* zFile,OpenType otype)
{
  int i,is_hdf5;
  hid_t file_id;
  unsigned int flags=H5F_ACC_RDWR;
  
  H5Eset_auto(H5E_DEFAULT,NULL,NULL);
  
  is_hdf5=H5Fis_hdf5(zFile);
  if(is_hdf5<=0 && otype==Read_OT){
    sprintf(errorBuffer,"File '%s' does not exists or is not in hdf5 format",zFile);
    return -1;
  }
  if(is_hdf5 && otype!=Create_OT){
    if(otype==Read_OT) flags=H5F_ACC_RDONLY;
    file_id=H5Fopen(zFile,flags,H5P_DEFAULT);
    if(file_id<0){
      sprintf(errorBuffer,"Could not open file '%s'",zFile);
      return -1;
    }
  } else {
    file_id=H5Fcreate(zFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(file_id<0){
      sprintf(errorBuffer,"could not create file '%s'",zFile);
      return -1;
    }
  }
  global_file_id=file_id;
  for(i=0;i<MAX_CONCURRENT_DATASETS;i++){
    global_hdf5c_buffer[i].used=0;
    global_hdf5c_buffer[i].file_id=file_id;
    memset(global_hdf5c_buffer[i].intarray,0,20*sizeof(int*));
    memset(global_hdf5c_buffer[i].doublearray,0,20*sizeof(double*));
  }
  strcpy(errorBuffer,"");
  return 0;
}

hid_t hdf5cM_open(const char* zFile,OpenType otype)
{
  int is_hdf5;
  hid_t file_id;
  unsigned int flags=H5F_ACC_RDWR;
  
  H5Eset_auto(H5E_DEFAULT,NULL,NULL);
  
  is_hdf5=H5Fis_hdf5(zFile);
  if(is_hdf5<=0 && otype==Read_OT){
    sprintf(errorBuffer,"File '%s' does not exists or is not in hdf5 format",zFile);
    return -1;
  }
  if(is_hdf5 && otype!=Create_OT){
    if(otype==Read_OT) flags=H5F_ACC_RDONLY;
    file_id=H5Fopen(zFile,flags,H5P_DEFAULT);
    if(file_id<0){
      sprintf(errorBuffer,"Could not open file '%s'",zFile);
      return -1;
    }
  } else {
    file_id=H5Fcreate(zFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(file_id<0){
      sprintf(errorBuffer,"could not create file '%s'",zFile);
      return -1;
    }
  }
  strcpy(errorBuffer,"");
  return file_id;
}

int hdf5c_init(const char* zFile)
{
  return hdf5c_open(zFile,Create_OT);
}

int hdf5c_end()
{
  int i,j;
  herr_t ht_ret;
  if(global_group_id>=0){
    H5Gclose(global_group_id);
    global_group_id=0;
  }
  if(global_hdf5c_buffer[0].file_id>0)
    ht_ret=H5Fclose (global_hdf5c_buffer[0].file_id);
  for(i=0;i<MAX_CONCURRENT_DATASETS;i++){
    for(j=0;j<20;j++){
      if(global_hdf5c_buffer[i].intarray[j]) free(global_hdf5c_buffer[i].intarray[j]);
      if(global_hdf5c_buffer[i].doublearray[j]) free(global_hdf5c_buffer[i].doublearray[j]);
    }
  }
  return ht_ret;
}

int hdf5cM_end(hid_t file_id)
{
  herr_t ht_ret;
  
  ht_ret=H5Fclose(file_id);
  return ht_ret;
}

int hdf5c_flush()
{
  return H5Fflush(global_file_id,H5F_SCOPE_GLOBAL);
}

int hdf5cM_flush(hid_t file_id)
{
  return H5Fflush(file_id,H5F_SCOPE_GLOBAL);
}

char* hdf5c_last_error()
{
  return errorBuffer;
}

/*################################################################################
 *    write functions
 *################################################################################*/

int hdf5c_create_group(const char* name)
{
  hid_t ht_group;

  ht_group=H5Gopen(global_file_id,name,H5P_DEFAULT);
  if(ht_group>=0){
    H5Gclose(ht_group);
    return 0;
  }
  ht_group=H5Gcreate(global_file_id,name,0,H5P_DEFAULT,H5P_DEFAULT);
  if(ht_group<0){
    sprintf(errorBuffer,"name '%s'. error creating group",name);
    return -1;
  }
  H5Gclose (ht_group);
  return 1;
}

int hdf5c_start_dataset(const char* name,int num_int,int num_real)
{
  int i,dataset_id;
  hid_t ht_group;
  
  for(dataset_id=0;dataset_id<MAX_CONCURRENT_DATASETS;dataset_id++)
    if(!global_hdf5c_buffer[dataset_id].used) break;
  if(dataset_id==MAX_CONCURRENT_DATASETS){
    sprintf(errorBuffer,"Too many concurrent datasets n=%d",MAX_CONCURRENT_DATASETS);
    return -1;
  }
  
  if(strlen(name)>1023){
    sprintf(errorBuffer,"name '%s'. too long",name);
    return -1;
  }
  strcpy(global_hdf5c_buffer[dataset_id].datasetname,name);
  ht_group=H5Gopen(global_hdf5c_buffer[dataset_id].file_id,name,H5P_DEFAULT);
  if(ht_group<0){
    ht_group=H5Gcreate(global_hdf5c_buffer[dataset_id].file_id,name,0,H5P_DEFAULT,H5P_DEFAULT);
    if(ht_group<0){
      sprintf(errorBuffer,"name '%s'. error creating group",name);
      return -1;
    }
  }
  H5Fclose (ht_group);
  
  global_hdf5c_buffer[dataset_id].used=1;
  if(global_hdf5c_buffer[dataset_id].memsize==0){
    global_hdf5c_buffer[dataset_id].memsize=10000;
  }
  for(i=global_hdf5c_buffer[dataset_id].num_int;i<num_int;i++){
    if(!global_hdf5c_buffer[dataset_id].intarray[i])
      global_hdf5c_buffer[dataset_id].intarray[i]=(int*) malloc(( size_t)global_hdf5c_buffer[dataset_id].memsize*sizeof(int));
  }
  global_hdf5c_buffer[dataset_id].num_int=num_int;
  for(i=global_hdf5c_buffer[dataset_id].num_real;i<num_real;i++){
    if(!global_hdf5c_buffer[dataset_id].doublearray[i])
      global_hdf5c_buffer[dataset_id].doublearray[i]=(double*) malloc(( size_t)global_hdf5c_buffer[dataset_id].memsize*sizeof(double));
  }
  global_hdf5c_buffer[dataset_id].num_real=num_real;
  global_hdf5c_buffer[dataset_id].size=0;
  return dataset_id;
}

int hdf5c_addto_dataset(int dataset_id,int intvalues[],double doublevalues[])
{
  int i;
  if(global_hdf5c_buffer[dataset_id].size==global_hdf5c_buffer[dataset_id].memsize){
    global_hdf5c_buffer[dataset_id].memsize*=2;
    for(i=0;i<20;i++){
      if(global_hdf5c_buffer[dataset_id].intarray[i])
	global_hdf5c_buffer[dataset_id].intarray[i]=realloc(global_hdf5c_buffer[dataset_id].intarray[i],
							    ( size_t)global_hdf5c_buffer[dataset_id].memsize*sizeof(int));
      if(global_hdf5c_buffer[dataset_id].doublearray[i])
	global_hdf5c_buffer[dataset_id].doublearray[i]=realloc(global_hdf5c_buffer[dataset_id].doublearray[i],
							       ( size_t)global_hdf5c_buffer[dataset_id].memsize*sizeof(double));
    }
  }
  for(i=0;i<global_hdf5c_buffer[dataset_id].num_int;i++)
    global_hdf5c_buffer[dataset_id].intarray[i][global_hdf5c_buffer[dataset_id].size]=intvalues[i];
  for(i=0;i<global_hdf5c_buffer[dataset_id].num_real;i++)
    global_hdf5c_buffer[dataset_id].doublearray[i][global_hdf5c_buffer[dataset_id].size]=doublevalues[i];
  global_hdf5c_buffer[dataset_id].size++;
  return 0;
}

int hdf5c_create_int_dataset(int dataset_id,char* name,int compress,int length,int* intarray)
{
  hsize_t dims[1],chunk_size[1] ;
  hid_t ht_err,plist,sid,dataset;

  dims[0]= ( hsize_t)length;
  if(compress>0){
    chunk_size[0]= ( hsize_t)length;
    plist=H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, 1, chunk_size);
    H5Pset_deflate(plist, ( unsigned int)compress);
    sid = H5Screate_simple(1, dims, NULL);
    dataset = H5Dcreate(global_hdf5c_buffer[dataset_id].file_id,name,H5T_NATIVE_INT,sid,
		        H5P_DEFAULT,plist,H5P_DEFAULT);
    ht_err=H5Dwrite(dataset,H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,intarray);
    H5Dclose(dataset);
    H5Sclose(sid);
    H5Pclose(plist);
  } else {
    ht_err=H5LTmake_dataset_int(global_hdf5c_buffer[dataset_id].file_id,name,1,dims,intarray);
  }
  if(ht_err<0){
    sprintf(errorBuffer,"cannot set array name '%s'",name);
    return -1;
  }
  return 0;
}

int hdf5c_create_double_dataset(int dataset_id,char* name,int compress,int length,double* doublearray)
{
  hsize_t dims[1],chunk_size[1] ;
  hid_t ht_err,plist,sid,dataset;

  dims[0]=( hsize_t)length;
  if(compress>0){
    chunk_size[0]=( hsize_t)length;
    plist=H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, 1, chunk_size);
    H5Pset_deflate(plist, ( unsigned int)compress);
    sid = H5Screate_simple(1, dims, NULL);
    dataset = H5Dcreate(global_hdf5c_buffer[dataset_id].file_id,name,H5T_NATIVE_DOUBLE,sid,
		        H5P_DEFAULT,plist,H5P_DEFAULT);
    ht_err=H5Dwrite(dataset,H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,doublearray);
    H5Dclose(dataset);
    H5Sclose(sid);
    H5Pclose(plist);
  } else {
    ht_err=H5LTmake_dataset_double(global_hdf5c_buffer[dataset_id].file_id,name,1,dims,doublearray);
  }
  if(ht_err<0){
    sprintf(errorBuffer,"cannot set array name '%s'",name);
    return -1;
  }
  return 0;
}
  
int hdf5c_end_dataset(int dataset_id)
{
  int i,fail,inum=1;
  int compress=6;
  char buf[2048];
  for(i=0;i<global_hdf5c_buffer[dataset_id].num_int;i++){
    sprintf(buf,"%s/%d",global_hdf5c_buffer[dataset_id].datasetname,inum);
    fail=hdf5c_create_int_dataset(dataset_id,buf,compress,global_hdf5c_buffer[dataset_id].size,
      global_hdf5c_buffer[dataset_id].intarray[i]);
    if(fail<0) return fail;
    inum++;
  }
  for(i=0;i<global_hdf5c_buffer[dataset_id].num_real;i++){
    sprintf(buf,"%s/%d",global_hdf5c_buffer[dataset_id].datasetname,inum);
    fail=hdf5c_create_double_dataset(dataset_id,buf,compress,global_hdf5c_buffer[dataset_id].size,
      global_hdf5c_buffer[dataset_id].doublearray[i]);
    if(fail<0) return fail;
    inum++;
  }
  global_hdf5c_buffer[dataset_id].used=0;
  return 0;
}
 
int hdf5c_set_attribute(const char* name,const char* attname,const char* value)
{
  hid_t ht_err;
  
  /*
    ht_err=H5LTset_attribute_char(global_hdf5c_buffer.file_id,name,attname,value,strlen(value));
  */
  ht_err=H5LTset_attribute_string(global_file_id,name,attname,value);
  if(ht_err<0){
    sprintf(errorBuffer,"cannot set attribute to '%s'",name);
    return -1;
  }
  return 0;
}

/*################################################################################
 *    Read functions
 *################################################################################*/

int hdf5c_num_children(const char* group_path)
{
  hid_t ht_err;
  hsize_t num_obj;

  if(global_group_id>=0) H5Gclose(global_group_id);
  global_group_id=H5Gopen(global_file_id,group_path,H5P_DEFAULT);
  if(global_group_id<0){
    sprintf(errorBuffer,"Could not find group '%s'",group_path);
    return -1;
  }
  strcpy(global_group_path,group_path);
  
  ht_err=H5Gget_num_objs(global_group_id,&num_obj);

  if(ht_err<0){
    sprintf(errorBuffer,"group '%s' is not a correct group",group_path);
    return -1;
  }
  return (int)num_obj;
}

int hdf5cM_num_children(hid_t file_id,const char* group_path)
{
  hid_t ht_err;
  hsize_t num_obj;
  hid_t group_id;
  
  group_id=H5Gopen(file_id,group_path,H5P_DEFAULT);
  if(group_id<0){
    sprintf(errorBuffer,"Could not find group '%s'",group_path);
    return -1;
  }

  ht_err=H5Gget_num_objs(group_id,&num_obj);

  if(ht_err<0){
    sprintf(errorBuffer,"group '%s' is not a correct group",group_path);
    return -1;
  }
  H5Gclose(group_id);
  return (int)num_obj;
}

int hdf5c_get_children_name(const char* group_path,int index,char* name,int max_len)
{
  ssize_t len,len0;
  if(global_group_id<=0 || strcmp(group_path,global_group_path)!=0){
    if(global_group_id>=0) H5Gclose(global_group_id);
    global_group_id=H5Gopen(global_file_id,group_path,H5P_DEFAULT);
    if(global_group_id<0){
      sprintf(errorBuffer,"Could not find group '%s'",group_path);
      return -1;
    }
    strcpy(global_group_path,group_path);
  }
  if((int)strlen(group_path)>=max_len-10) return -1;
  strcpy(name,group_path);
  if(name[strlen(name)-1]!='/'){
    strcat(name,"/");
  }
  len0= ( ssize_t)strlen(name);
  len=H5Gget_objname_by_idx(global_group_id,index,&name[len0],max_len-len0);
  return ( int)( len+len0);
}

int hdf5cM_get_children_name(hid_t file_id,const char* group_path,int index,char* name,int max_len)
{  
  hid_t group_id;
  ssize_t len,len0;
  group_id=H5Gopen(file_id,group_path,H5P_DEFAULT);

  if(group_id<=0 ){
    if(group_id>=0) H5Gclose(group_id);
    group_id=H5Gopen(file_id,group_path,H5P_DEFAULT);
    if(group_id<0){
      sprintf(errorBuffer,"Could not find group '%s'",group_path);
      return -1;
    }
  }  
  if((int)strlen(group_path)>=max_len-10){
    H5Gclose(group_id);
    return -1;
  }
  strcpy(name,group_path);
  if(name[strlen(name)-1]!='/'){
    strcat(name,"/");
  }
  len0= ( ssize_t)strlen(name);
  len=H5Gget_objname_by_idx(group_id,index,&name[len0],max_len-len0);
  H5Gclose(group_id);
  return ( int)( len+len0);
}

int hdf5c_get_attribute(const char* obj_name,char* att_name,char* value,int max_len)
{
  int rank;
  hid_t ht_err;
  hsize_t dims[10];
  size_t  type_size;
  H5T_class_t class_id;
  
  ht_err=H5LTget_attribute_info(global_file_id,obj_name,att_name, dims,&class_id,&type_size);
  if(ht_err<0){
    sprintf(errorBuffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
    return -1;
  }  
  if(class_id==H5T_INTEGER) {
    if(type_size!=1){
      sprintf(errorBuffer,"name '%s' cannot be retrieved. Bad type %d. Has to be one array of char",
	obj_name,class_id);
      return -1;
    }
    ht_err=H5LTget_attribute_ndims(global_file_id,obj_name,att_name,&rank);
    /*
      rank is 1 for H5T_INTEGER
    */
    if(ht_err<0 || rank!=1){
      sprintf(errorBuffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
      return -1;
    }
    if(dims[0]>= ( hsize_t)max_len) return -1;
    ht_err=H5LTget_attribute_char(global_file_id,obj_name,att_name,value);
  } else if (class_id==H5T_STRING){
    if(type_size<1){
      sprintf(errorBuffer,"name '%s' cannot be retrieved. Bad type %d. Has to be one array of char",
	obj_name,class_id);
      return -1;
    }
    ht_err=H5LTget_attribute_ndims(global_file_id,obj_name,att_name,&rank);
    /*
      rank is 0 for H5T_STRING
    */
    if(ht_err<0 || rank!=0){
      sprintf(errorBuffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
      return -1;
    }
    if((int)type_size>=max_len) return -1;
    ht_err=H5LTget_attribute_string(global_file_id,obj_name,att_name,value);
  } else {
    sprintf(errorBuffer,"name '%s' cannot be retrieved. Bad type %d Has to be one array of char or a string",
      obj_name,class_id);
    return -1;
  }
  return 0;
}

int hdf5cM_get_attribute(hid_t file_id,const char* obj_name,char* att_name,char* value,int max_len)
{
  int rank;
  hid_t ht_err;
  hsize_t dims[10];
  size_t  type_size;
  H5T_class_t class_id;
  
  ht_err=H5LTget_attribute_info(file_id,obj_name,att_name, dims,&class_id,&type_size);
  if(ht_err<0){
    sprintf(errorBuffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
    return -1;
  }  
  if(class_id==H5T_INTEGER) {
    if(type_size!=1){
      sprintf(errorBuffer,"name '%s' cannot be retrieved. Bad type %d. Has to be one array of char",
	obj_name,class_id);
      return -1;
    }
    ht_err=H5LTget_attribute_ndims(file_id,obj_name,att_name,&rank);
    /*
      rank is 1 for H5T_INTEGER
    */
    if(ht_err<0 || rank!=1){
      sprintf(errorBuffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
      return -1;
    }
    if(dims[0]>=( hsize_t)max_len) return -1;
    ht_err=H5LTget_attribute_char(file_id,obj_name,att_name,value);
  } else if (class_id==H5T_STRING){
    if(type_size<1){
      sprintf(errorBuffer,"name '%s' cannot be retrieved. Bad type %d. Has to be one array of char",
	obj_name,class_id);
      return -1;
    }
    ht_err=H5LTget_attribute_ndims(file_id,obj_name,att_name,&rank);
    /*
      rank is 0 for H5T_STRING
    */
    if(ht_err<0 || rank!=0){
      sprintf(errorBuffer,"name '%s'. Attribute cannot be retrieved '%s'",obj_name,att_name);
      return -1;
    }
    if((int)type_size>=max_len) return -1;
    ht_err=H5LTget_attribute_string(file_id,obj_name,att_name,value);
  } else {
    sprintf(errorBuffer,"name '%s' cannot be retrieved. Bad type %d Has to be one array of char or a string",
      obj_name,class_id);
    return -1;
  }
  return 0;
}

int hdf5c_open_dataset_list(const char* obj_name,int* num_int,int* num_real)
{
  int i,rank,idx;
  hsize_t num;
  char name[1024];
  
  hid_t ht_err,ht_dataset,ht_type;
  hsize_t dims[1];
  H5T_class_t class_id;
  size_t type_size;

  num=( hsize_t)hdf5c_num_children(obj_name);
  if(num<=0) return -1;
  
  global_hdf5c_buffer[0].num_int=global_hdf5c_buffer[0].num_real=0;
  
  for(i=0;i< ( int)num;i++){
    hdf5c_get_children_name(obj_name,i,name,1024);
    ht_err=H5LTget_dataset_ndims(global_file_id,name,&rank);
    if(ht_err<0 || rank!=1){
      sprintf(errorBuffer,"dataset name '%s' cannot be retrieved",name);
      return -1;
    }
    ht_err=H5LTget_dataset_info (global_file_id,name,dims,&class_id,NULL);
    if(ht_err<0){
      sprintf(errorBuffer,"dataset name '%s' cannot be retrieved",name);
      return -1;
    }
    ht_dataset=H5Dopen(global_file_id,name,H5P_DEFAULT);
    ht_type=H5Dget_type(ht_dataset);
    type_size=H5Tget_size(ht_type);
    H5Dclose(ht_dataset);
    H5Tclose(ht_type);
    
    if(class_id==H5T_INTEGER){
      idx=global_hdf5c_buffer[0].num_int;
      global_hdf5c_buffer[0].intarray[idx]=malloc(((size_t)dims[0])*sizeof(int));
      ht_err=H5LTread_dataset_int (global_file_id,name,global_hdf5c_buffer[0].intarray[idx]);
      if(ht_err<0){
	sprintf(errorBuffer,"array name '%s' cannot be retrieved",name);
	return -1;
      }
      global_hdf5c_buffer[0].num_int++;
    }
    else if(class_id==H5T_FLOAT && type_size==4){
      sprintf(errorBuffer,"4 byte floats not supported");
      return -1;
    }
    else if(class_id==H5T_FLOAT && type_size==8){
      idx=global_hdf5c_buffer[0].num_real;
      global_hdf5c_buffer[0].doublearray[idx]=malloc(((size_t)dims[0])*sizeof(double));
      ht_err=H5LTread_dataset_double(global_file_id,name,global_hdf5c_buffer[0].doublearray[idx]);
      if(ht_err<0){
	sprintf(errorBuffer,"array name '%s' cannot be retrieved",name);
	return -1;
      }
      global_hdf5c_buffer[0].num_real++;
    }
    else {
      sprintf(errorBuffer,"array name '%s' cannot be retrieved. unknown type=%d,%lud",
	name,class_id,type_size);
      return -1;
    }
    global_hdf5c_buffer[0].size=(int)dims[0];
  }
  global_hdf5c_buffer[0].used=1;
  *num_int=global_hdf5c_buffer[0].num_int;
  *num_real=global_hdf5c_buffer[0].num_real;
  return global_hdf5c_buffer[0].size;
}

int hdf5cM_open_dataset_list(hid_t file_id,struct hdf5c_buffer* hdf5c_buffer,const char* obj_name,
  int* num_int,int* num_real)
{
  int i,rank,idx;
  hsize_t num;
  char name[1024];
  
  hid_t ht_err,ht_dataset,ht_type;
  hsize_t dims[1];
  H5T_class_t class_id;
  size_t type_size;

  num=( hsize_t)hdf5cM_num_children(file_id,obj_name);
  if(num<=0) return -1;  
  
  hdf5c_buffer->used=0;
  hdf5c_buffer->file_id=file_id;
  memset(hdf5c_buffer->intarray,0,20*sizeof(int*));
  memset(hdf5c_buffer->doublearray,0,20*sizeof(double*));
  hdf5c_buffer->num_int=hdf5c_buffer->num_real=0;
  
  for(i=0;i< ( int)num;i++){
    hdf5cM_get_children_name(file_id,obj_name,i,name,1024);
    ht_err=H5LTget_dataset_ndims(file_id,name,&rank);
    if(ht_err<0 || rank!=1){
      sprintf(errorBuffer,"dataset name '%s' cannot be retrieved",name);
      return -1;
    }
    ht_err=H5LTget_dataset_info (file_id,name,dims,&class_id,NULL);
    if(ht_err<0){
      sprintf(errorBuffer,"dataset name '%s' cannot be retrieved",name);
      return -1;
    }
    ht_dataset=H5Dopen(file_id,name,H5P_DEFAULT);
    ht_type=H5Dget_type(ht_dataset);
    type_size=H5Tget_size(ht_type);
    H5Dclose(ht_dataset);
    H5Tclose(ht_type);
    
    if(class_id==H5T_INTEGER){
      idx=hdf5c_buffer->num_int;
      hdf5c_buffer->intarray[idx]=malloc(((size_t)dims[0])*sizeof(int));
      ht_err=H5LTread_dataset_int (file_id,name,hdf5c_buffer->intarray[idx]);
      if(ht_err<0){
	sprintf(errorBuffer,"array name '%s' cannot be retrieved",name);
	return -1;
      }
      hdf5c_buffer->num_int++;
    }
    else if(class_id==H5T_FLOAT && type_size==4){
      sprintf(errorBuffer,"4 byte floats not supported");
      return -1;
    }
    else if(class_id==H5T_FLOAT && type_size==8){
      idx=hdf5c_buffer->num_real;
      hdf5c_buffer->doublearray[idx]=malloc(((size_t)dims[0])*sizeof(double));
      ht_err=H5LTread_dataset_double(file_id,name,hdf5c_buffer->doublearray[idx]);
      if(ht_err<0){
	sprintf(errorBuffer,"array name '%s' cannot be retrieved",name);
	return -1;
      }
      hdf5c_buffer->num_real++;
    }
    else {
      sprintf(errorBuffer,"array name '%s' cannot be retrieved. unknown type=%d,%lud",
	name,class_id,type_size);
      return -1;
    }
    hdf5c_buffer->size=(int)dims[0];
  }
  hdf5c_buffer->used=1;
  *num_int=hdf5c_buffer->num_int;
  *num_real=hdf5c_buffer->num_real;
  return hdf5c_buffer->size;
}

int hdf5c_give_dataset_list_values(int index,int int_values[],double double_values[])
{
  int i;
  if(index<0 || index>=global_hdf5c_buffer[0].size){
    return -1;
  }
  for(i=0;i<global_hdf5c_buffer[0].num_int;i++){
    int_values[i]=global_hdf5c_buffer[0].intarray[i][index];
  }
  for(i=0;i<global_hdf5c_buffer[0].num_real;i++){
    double_values[i]=global_hdf5c_buffer[0].doublearray[i][index];
  }
  return 0;
}

int hdf5cM_give_dataset_list_values(struct hdf5c_buffer* hdf5c_buffer,int index,int int_values[],double double_values[])
{
  int i;
  if(index<0 || index>=hdf5c_buffer->size){
    return -1;
  }
  for(i=0;i<hdf5c_buffer->num_int;i++){
    int_values[i]=hdf5c_buffer->intarray[i][index];
  }
  for(i=0;i<hdf5c_buffer->num_real;i++){
    double_values[i]=hdf5c_buffer->doublearray[i][index];
  }
  return 0;
}

int hdf5c_close_dataset_list()
{
  int i;
  for(i=0;i<global_hdf5c_buffer[0].num_int;i++){
    free(global_hdf5c_buffer[0].intarray[i]);
    global_hdf5c_buffer[0].intarray[i]=NULL;
  }
  for(i=0;i<global_hdf5c_buffer[0].num_real;i++){
    free(global_hdf5c_buffer[0].doublearray[i]);
    global_hdf5c_buffer[0].doublearray[i]=NULL;
  }
  global_hdf5c_buffer[0].num_int=global_hdf5c_buffer[0].num_real=0;
  global_hdf5c_buffer[0].used=0;
  return 0;
}

int hdf5cM_close_dataset_list(struct hdf5c_buffer* hdf5c_buffer)
{
  int i,j;
  for(i=0;i<hdf5c_buffer->num_int;i++){
    free(hdf5c_buffer->intarray[i]);
    hdf5c_buffer->intarray[i]=NULL;
  }
  for(i=0;i<hdf5c_buffer->num_real;i++){
    free(hdf5c_buffer->doublearray[i]);
    hdf5c_buffer->doublearray[i]=NULL;
  }
  hdf5c_buffer->num_int=hdf5c_buffer->num_real=0;
  hdf5c_buffer->used=0;
  
  for(j=0;j<20;j++){
    if(hdf5c_buffer->intarray[j]) free(hdf5c_buffer->intarray[j]);
    if(hdf5c_buffer->doublearray[j]) free(hdf5c_buffer->doublearray[j]);
  }
  return 0;
}

/*###############################################################################
 *    All functions until the end are literally copied from library: HDF5_HL
 *###############################################################################*/

/*-------------------------------------------------------------------------
* Function: H5LTset_attribute_string
*
* Purpose: Creates and writes a string attribute named attr_name and attaches
*          it to the object specified by the name obj_name.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: July 23, 2001
*
* Comments: If the attribute already exists, it is overwritten
*
* Modifications:
*
*-------------------------------------------------------------------------
*/

herr_t H5LTset_attribute_string( hid_t loc_id,
		                const char *obj_name,
		                const char *attr_name,
		                const char *attr_data )
{
    hid_t      attr_type;
    hid_t      attr_space_id;
    hid_t      attr_id;
    hid_t      obj_id;
    int        has_attr;
    size_t     attr_size;

    /* Open the object */
    if ((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
	return -1;

    /* Create the attribute */
    if ( (attr_type = H5Tcopy( H5T_C_S1 )) < 0 )
	goto out;

    attr_size = strlen( attr_data ) + 1; /* extra null term */

    if ( H5Tset_size( attr_type, (size_t)attr_size) < 0 )
	goto out;

    if ( H5Tset_strpad( attr_type, H5T_STR_NULLTERM ) < 0 )
	goto out;

    if ( (attr_space_id = H5Screate( H5S_SCALAR )) < 0 )
	goto out;

    /* Verify if the attribute already exists */
    has_attr = H5LT_find_attribute(obj_id, attr_name);

    /* The attribute already exists, delete it */
    if(has_attr == 1)
	if(H5Adelete(obj_id, attr_name) < 0)
	    goto out;

    /* Create and write the attribute */

    if((attr_id = H5Acreate2(obj_id, attr_name, attr_type, attr_space_id, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	goto out;

    if(H5Awrite(attr_id, attr_type, attr_data) < 0)
	goto out;

    if(H5Aclose(attr_id) < 0)
	goto out;

    if(H5Sclose(attr_space_id) < 0)
	goto out;

    if(H5Tclose(attr_type) < 0)
	goto out;

    /* Close the object */
    if(H5Oclose(obj_id) < 0)
	return -1;

    return 0;

out:

    H5Oclose(obj_id);
    return -1;
}





/*-------------------------------------------------------------------------
* Function: H5LT_set_attribute_numerical
*
* Purpose: Private function used by H5LTset_attribute_int and H5LTset_attribute_float
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: July 25, 2001
*
* Comments:
*
*-------------------------------------------------------------------------
*/


herr_t H5LT_set_attribute_numerical( hid_t loc_id,
		                    const char *obj_name,
		                    const char *attr_name,
		                    size_t size,
		                    hid_t tid,
		                    const void *data )
{

    hid_t      obj_id, sid, attr_id;
    hsize_t    dim_size=size;
    int        has_attr;

    /* Open the object */
    if ((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
	return -1;

    /* Create the data space for the attribute. */
    if ( (sid = H5Screate_simple( 1, &dim_size, NULL )) < 0 )
	goto out;

    /* Verify if the attribute already exists */
    has_attr = H5LT_find_attribute(obj_id, attr_name);

    /* The attribute already exists, delete it */
    if(has_attr == 1)
	if(H5Adelete(obj_id, attr_name) < 0)
	    goto out;

    /* Create the attribute. */
    if((attr_id = H5Acreate2(obj_id, attr_name, tid, sid, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	goto out;

    /* Write the attribute data. */
    if(H5Awrite(attr_id, tid, data) < 0)
	goto out;

    /* Close the attribute. */
    if(H5Aclose(attr_id) < 0)
	goto out;

    /* Close the dataspace. */
    if(H5Sclose(sid) < 0)
	goto out;

    /* Close the object */
    if(H5Oclose(obj_id) < 0)
	return -1;

    return 0;

out:
    H5Oclose(obj_id);
    return -1;
}

/*-------------------------------------------------------------------------
* Function: H5LT_make_dataset
*
* Purpose: Creates and writes a dataset of a type tid
*
* Return: Success: 0, Failure: -1
*
* Programmer: Quincey Koziol, koziol@hdfgroup.org
*
* Date: October 10, 2007
*
*-------------------------------------------------------------------------
*/

static herr_t
H5LT_make_dataset_numerical( hid_t loc_id,
		            const char *dset_name,
		            int rank,
		            const hsize_t *dims,
		            hid_t tid,
		            const void *data )
{
    hid_t   did = -1, sid = -1;

    /* Create the data space for the dataset. */
    if((sid = H5Screate_simple(rank, dims, NULL)) < 0)
	return -1;

    /* Create the dataset. */
    if((did = H5Dcreate2(loc_id, dset_name, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)
	goto out;

    /* Write the dataset only if there is data to write */
    if(data)
	if(H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
	    goto out;

    /* End access to the dataset and release resources used by it. */
    if(H5Dclose(did) < 0)
	return -1;

    /* Terminate access to the data space. */
    if(H5Sclose(sid) < 0)
	return -1;

    return 0;

out:
    H5E_BEGIN_TRY {
	H5Dclose(did);
	H5Sclose(sid);
    } H5E_END_TRY;
    return -1;
}

/*-------------------------------------------------------------------------
* Function: H5LTget_dataset_ndims
*
* Purpose: Gets the dimensionality of a dataset.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 4, 2001
*
*-------------------------------------------------------------------------
*/

herr_t H5LTget_dataset_ndims( hid_t loc_id,
		             const char *dset_name,
		             int *rank )
{
    hid_t       did = -1;
    hid_t       sid = -1;

    /* Open the dataset. */
    if((did = H5Dopen2(loc_id, dset_name, H5P_DEFAULT)) < 0)
	return -1;

    /* Get the dataspace handle */
    if((sid = H5Dget_space(did)) < 0)
	goto out;

    /* Get rank */
    if((*rank = H5Sget_simple_extent_ndims(sid)) < 0)
	goto out;

    /* Terminate access to the dataspace */
    if(H5Sclose(sid) < 0)
	goto out;

    /* End access to the dataset */
    if(H5Dclose(did))
	return -1;

    return 0;

out:
    H5E_BEGIN_TRY {
	H5Dclose(did);
	H5Sclose(sid);
    } H5E_END_TRY;
    return -1;
}

/*-------------------------------------------------------------------------
* Function: H5LTget_dataset_info
*
* Purpose: Gets information about a dataset.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 4, 2001
*  Modified: February 28, 2006: checked for NULL parameters
*
*-------------------------------------------------------------------------
*/

herr_t H5LTget_dataset_info( hid_t loc_id,
		            const char *dset_name,
		            hsize_t *dims,
		            H5T_class_t *type_class,
		            size_t *type_size )
{
    hid_t       did = -1;
    hid_t       tid = -1;
    hid_t       sid = -1;

    /* open the dataset. */
    if((did = H5Dopen2(loc_id, dset_name, H5P_DEFAULT)) < 0)
	return -1;

    /* get an identifier for the datatype. */
    tid = H5Dget_type(did);

    /* get the class. */
    if(type_class != NULL)
	*type_class = H5Tget_class(tid);

    /* get the size. */
    if(type_size!=NULL)
	*type_size = H5Tget_size(tid);

    if(dims != NULL) {
	/* get the dataspace handle */
	if((sid = H5Dget_space(did)) < 0)
	    goto out;

	/* get dimensions */
	if(H5Sget_simple_extent_dims(sid, dims, NULL) < 0)
	    goto out;

	/* terminate access to the dataspace */
	if(H5Sclose(sid) < 0)
	    goto out;
    } /* end if */

    /* release the datatype. */
    if(H5Tclose(tid))
	return -1;

    /* end access to the dataset */
    if(H5Dclose(did))
	return -1;

    return 0;

out:
    H5E_BEGIN_TRY {
	H5Tclose(tid);
	H5Sclose(sid);
	H5Dclose(did);
    } H5E_END_TRY;
    return -1;

}


/*-------------------------------------------------------------------------
* Function: find_attr
*
* Purpose: operator function used by H5LT_find_attribute
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: June 21, 2001
*
* Comments:
*
* Modifications:
*
*-------------------------------------------------------------------------
*/
static herr_t
find_attr(hid_t loc_id, const char *name, const H5A_info_t *ainfo,
	  void *op_data)
{
    int ret = H5_ITER_CONT;

    /* Shut compiler up */
    loc_id = loc_id; ainfo = ainfo;

    /* Define a positive value for return value if the attribute was found. This will
    * cause the iterator to immediately return that positive value,
    * indicating short-circuit success
    */
    if(strcmp(name, (char *)op_data) == 0)
	ret = H5_ITER_STOP;

    return ret;
}

/*-------------------------------------------------------------------------
* Function: H5LTfind_attribute
*
* Purpose: Inquires if an attribute named attr_name exists attached to
*          the object loc_id.
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: May 17, 2006
*
* Comments:
*  Calls the private version of the function
*
*-------------------------------------------------------------------------
*/

herr_t H5LTfind_attribute( hid_t loc_id, const char* attr_name )
{
    return H5LT_find_attribute(loc_id,attr_name);
}



/*-------------------------------------------------------------------------
* Function: H5LT_find_attribute
*
* Purpose: Inquires if an attribute named attr_name exists attached to the object loc_id.
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: June 21, 2001
*
* Comments:
*  The function uses H5Aiterate2 with the operator function find_attr
*
* Return:
*  Success: The return value of the first operator that
*              returns non-zero, or zero if all members were
*              processed with no operator returning non-zero.
*
*  Failure: Negative if something goes wrong within the
*              library, or the negative value returned by one
*              of the operators.
*
*-------------------------------------------------------------------------
*/

herr_t
H5LT_find_attribute( hid_t loc_id, const char* attr_name )
{
    return H5Aiterate2(loc_id, H5_INDEX_NAME, H5_ITER_INC, NULL, find_attr, (void *)attr_name);
}

/*-------------------------------------------------------------------------
* Function: H5LTget_attribute_ndims
*
* Purpose: Gets the dimensionality of an attribute.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 4, 2001
*
*-------------------------------------------------------------------------
*/

herr_t H5LTget_attribute_ndims( hid_t loc_id,
		               const char *obj_name,
		               const char *attr_name,
		               int *rank )
{
    hid_t       attr_id;
    hid_t       sid;
    hid_t       obj_id;

    /* Open the object */
    if((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
	return -1;

    /* Open the attribute. */
    if((attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT)) < 0)
    {
	H5Oclose(obj_id);
	return -1;
    }

    /* Get the dataspace handle */
    if((sid = H5Aget_space(attr_id)) < 0)
	goto out;

    /* Get rank */
    if((*rank = H5Sget_simple_extent_ndims(sid)) < 0)
	goto out;

    /* Terminate access to the attribute */
    if ( H5Sclose( sid ) < 0 )
	goto out;

    /* End access to the attribute */
    if ( H5Aclose( attr_id ) )
	goto out;;

    /* Close the object */
    if(H5Oclose(obj_id) < 0 )
	return -1;

    return 0;

out:
    H5Aclose( attr_id );
    H5Oclose(obj_id);
    return -1;

}

/*-------------------------------------------------------------------------
* Function: H5LTget_attribute_info
*
* Purpose: Gets information about an attribute.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 4, 2001
*
*-------------------------------------------------------------------------
*/

herr_t H5LTget_attribute_info( hid_t loc_id,
		              const char *obj_name,
		              const char *attr_name,
		              hsize_t *dims,
		              H5T_class_t *type_class,
		              size_t *type_size )
{
    hid_t       attr_id;
    hid_t       tid;
    hid_t       sid;
    hid_t       obj_id;

    /* Open the object */
    if((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
	return -1;

    /* Open the attribute. */
    if((attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT)) < 0)
    {
	H5Oclose(obj_id);
	return -1;
    }

    /* Get an identifier for the datatype. */
    tid = H5Aget_type(attr_id);

    /* Get the class. */
    *type_class = H5Tget_class(tid);

    /* Get the size. */
    *type_size = H5Tget_size( tid );

    /* Get the dataspace handle */
    if ( (sid = H5Aget_space( attr_id )) < 0 )
	goto out;

    /* Get dimensions */
    if ( H5Sget_simple_extent_dims( sid, dims, NULL) < 0 )
	goto out;

    /* Terminate access to the dataspace */
    if ( H5Sclose( sid ) < 0 )
	goto out;

    /* Release the datatype. */
    if ( H5Tclose( tid ) )
	goto out;

    /* End access to the attribute */
    if ( H5Aclose( attr_id ) )
	goto out;

    /* Close the object */
    if(H5Oclose(obj_id) < 0 )
	return -1;

    return 0;

out:
    H5Tclose(tid);
    H5Aclose(attr_id);
    H5Oclose(obj_id);
    return -1;

}

  /*-------------------------------------------------------------------------
* Function: H5LTset_attribute_char
*
* Purpose: Create and write an attribute.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: November 7, 2001
*
* Comments:
*
*-------------------------------------------------------------------------
*/

herr_t H5LTset_attribute_char( hid_t loc_id,
		              const char *obj_name,
		              const char *attr_name,
		              const char *data,
		              size_t size )
{

    if ( H5LT_set_attribute_numerical( loc_id, obj_name, attr_name, size,
	H5T_NATIVE_CHAR, data ) < 0 )
	return -1;

    return 0;
}

/*-------------------------------------------------------------------------
* Function: H5LT_get_attribute_mem
*
* Purpose: Reads an attribute named attr_name with the memory type mem_type_id
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 19, 2002
*
* Comments: Private function
*
* Modifications:
*
*-------------------------------------------------------------------------
*/


static herr_t H5LT_get_attribute_mem(hid_t loc_id,
		                     const char *obj_name,
		                     const char *attr_name,
		                     hid_t mem_type_id,
		                     void *data)
{
    /* identifiers */
    hid_t obj_id = -1;
    hid_t attr_id = -1;

    /* Open the object */
    if((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
	goto out;

    if((attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT)) < 0)
	goto out;

    if(H5Aread(attr_id, mem_type_id, data) < 0)
	goto out;

    if(H5Aclose(attr_id) < 0)
	goto out;
    attr_id = -1;

    /* Close the object */
    if(H5Oclose(obj_id) < 0)
	goto out;
    obj_id = -1;

    return 0;

out:
    if(attr_id > 0)
	H5Aclose(attr_id);
    return -1;
}

/*-------------------------------------------------------------------------
* Function: H5LT_get_attribute_disk
*
* Purpose: Reads an attribute named attr_name with the datatype stored on disk
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 19, 2002
*
* Comments:
*
* Modifications:
*
*-------------------------------------------------------------------------
*/

herr_t H5LT_get_attribute_disk( hid_t loc_id,
		               const char *attr_name,
		               void *attr_out )
{
    /* identifiers */
    hid_t attr_id;
    hid_t attr_type;

    if(( attr_id = H5Aopen(loc_id, attr_name, H5P_DEFAULT)) < 0)
	return -1;

    if((attr_type = H5Aget_type(attr_id)) < 0)
	goto out;

    if(H5Aread(attr_id, attr_type, attr_out) < 0)
	goto out;

    if(H5Tclose(attr_type) < 0)
	goto out;

    if ( H5Aclose( attr_id ) < 0 )
	return -1;;

    return 0;

out:
    H5Tclose( attr_type );
    H5Aclose( attr_id );
    return -1;
}


/*-------------------------------------------------------------------------
* Function: H5LTget_attribute_string
*
* Purpose: Reads an attribute named attr_name
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 19, 2002
*
* Comments:
*
* Modifications:
*
*-------------------------------------------------------------------------
*/


herr_t H5LTget_attribute_string( hid_t loc_id,
		                const char *obj_name,
		                const char *attr_name,
		                char *data )
{
    /* identifiers */
    hid_t      obj_id;

    /* Open the object */
    if ((obj_id = H5Oopen( loc_id, obj_name, H5P_DEFAULT)) < 0)
	return -1;

    /* Get the attribute */
    if ( H5LT_get_attribute_disk( obj_id, attr_name, data ) < 0 )
	return -1;

    /* Close the object */
    if(H5Oclose(obj_id) < 0)
	return -1;

    return 0;

}


/*-------------------------------------------------------------------------
* Function: H5LTget_attribute_char
*
* Purpose: Reads an attribute named attr_name
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 19, 2002
*
* Comments:
*
* Modifications:
*
*-------------------------------------------------------------------------
*/
herr_t H5LTget_attribute_char( hid_t loc_id,
		              const char *obj_name,
		              const char *attr_name,
		              char *data )
{
    /* Get the attribute */
    if(H5LT_get_attribute_mem(loc_id, obj_name, attr_name, H5T_NATIVE_CHAR, data) < 0)
	return -1;

    return 0;
}


/*-------------------------------------------------------------------------
* Function: H5LTmake_dataset_char
*
* Purpose: Creates and writes a dataset of H5T_NATIVE_CHAR type
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 14, 2001
*
* Comments:
*
* Modifications:
*
*
*-------------------------------------------------------------------------
*/

herr_t H5LTmake_dataset_char( hid_t loc_id,
		             const char *dset_name,
		             int rank,
		             const hsize_t *dims,
		             const char *data )
{
    return(H5LT_make_dataset_numerical(loc_id, dset_name, rank, dims, H5T_NATIVE_CHAR, data));
}

/*-------------------------------------------------------------------------
* Function: H5LTmake_dataset_int
*
* Purpose: Creates and writes a dataset of H5T_NATIVE_INT type
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 14, 2001
*
* Comments:
*
* Modifications:
*
*
*-------------------------------------------------------------------------
*/



herr_t H5LTmake_dataset_int( hid_t loc_id,
		            const char *dset_name,
		            int rank,
		            const hsize_t *dims,
		            const int *data )
{
    return(H5LT_make_dataset_numerical(loc_id, dset_name, rank, dims, H5T_NATIVE_INT, data));
}

/*-------------------------------------------------------------------------
* Function: H5LTmake_dataset_double
*
* Purpose: Creates and writes a dataset of H5T_NATIVE_DOUBLE type
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 14, 2001
*
* Comments:
*
* Modifications:
*
*
*-------------------------------------------------------------------------
*/


herr_t H5LTmake_dataset_double( hid_t loc_id,
		               const char *dset_name,
		               int rank,
		               const hsize_t *dims,
		               const double *data )
{
    return(H5LT_make_dataset_numerical(loc_id, dset_name, rank, dims, H5T_NATIVE_DOUBLE, data));
}

/*-------------------------------------------------------------------------
* Function: H5LTmake_dataset_float
*
* Purpose: Creates and writes a dataset of H5T_NATIVE_FLOAT type
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: September 14, 2001
*
* Comments:
*
* Modifications:
*
*
*-------------------------------------------------------------------------
*/


herr_t H5LTmake_dataset_float( hid_t loc_id,
		              const char *dset_name,
		              int rank,
		              const hsize_t *dims,
		              const float *data )
{
    return(H5LT_make_dataset_numerical(loc_id, dset_name, rank, dims, H5T_NATIVE_FLOAT, data));
}

/*-------------------------------------------------------------------------
* Function: H5LT_read_dataset
*
* Purpose: Reads a dataset from disk.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Quincey Koziol, koziol@hdfgroup.org
*
* Date: October 8, 2007
*
*-------------------------------------------------------------------------
*/

static herr_t
H5LT_read_dataset_numerical(hid_t loc_id, const char *dset_name, hid_t tid, void *data)
{
    hid_t   did;

    /* Open the dataset. */
    if((did = H5Dopen2(loc_id, dset_name, H5P_DEFAULT)) < 0)
	return -1;

    /* Read */
    if(H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
	goto out;

    /* End access to the dataset and release resources used by it. */
    if(H5Dclose(did))
	return -1;

    return 0;

out:
    H5Dclose(did);
    return -1;
}

/*-------------------------------------------------------------------------
* Function: H5LTread_dataset_char
*
* Purpose: Reads a dataset from disk.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: November 5, 2001
*
*-------------------------------------------------------------------------
*/

herr_t H5LTread_dataset_char( hid_t loc_id,
		             const char *dset_name,
		             char *data )
{
    return(H5LT_read_dataset_numerical(loc_id, dset_name, H5T_NATIVE_CHAR, data));
}

/*-------------------------------------------------------------------------
* Function: H5LTread_dataset_int
*
* Purpose: Reads a dataset from disk.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: November 5, 2001
*
*-------------------------------------------------------------------------
*/

herr_t H5LTread_dataset_int( hid_t loc_id,
		            const char *dset_name,
		            int *data )
{
    return(H5LT_read_dataset_numerical(loc_id, dset_name, H5T_NATIVE_INT, data));
}

/*-------------------------------------------------------------------------
* Function: H5LTread_dataset_float
*
* Purpose: Reads a dataset from disk.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: November 5, 2001
*
*-------------------------------------------------------------------------
*/

herr_t H5LTread_dataset_float( hid_t loc_id,
		              const char *dset_name,
		              float *data )
{
    return(H5LT_read_dataset_numerical(loc_id, dset_name, H5T_NATIVE_FLOAT, data));
}

/*-------------------------------------------------------------------------
* Function: H5LTread_dataset_double
*
* Purpose: Reads a dataset from disk.
*
* Return: Success: 0, Failure: -1
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: November 5, 2001
*
*-------------------------------------------------------------------------
*/

herr_t H5LTread_dataset_double( hid_t loc_id,
		               const char *dset_name,
		               double *data )
{
    return(H5LT_read_dataset_numerical(loc_id, dset_name, H5T_NATIVE_DOUBLE, data));
}

/*-------------------------------------------------------------------------
* Function: find_dataset
*
* Purpose: operator function used by H5LTfind_dataset
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: June 21, 2001
*
* Comments:
*
* Modifications:
*
*-------------------------------------------------------------------------
*/

static herr_t
find_dataset(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *op_data)
{
    /* Define a default zero value for return. This will cause the iterator to continue if
    * the dataset is not found yet.
    */
    int ret = 0;

    /* Shut the compiler up */
    loc_id = loc_id;
    linfo = linfo;

    /* Define a positive value for return value if the dataset was found. This will
    * cause the iterator to immediately return that positive value,
    * indicating short-circuit success
    */
    if(strcmp(name, (char *)op_data) == 0)
	ret = 1;

    return ret;
}


/*-------------------------------------------------------------------------
* Function: H5LTfind_dataset
*
* Purpose:  Inquires if a dataset named dset_name exists attached
*           to the object loc_id.
*
* Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
*
* Date: July 15, 2001
*
* Return:
*     Success: The return value of the first operator that
*              returns non-zero, or zero if all members were
*              processed with no operator returning non-zero.
*
*      Failure:    Negative if something goes wrong within the
*              library, or the negative value returned by one
*              of the operators.
*
*-------------------------------------------------------------------------
*/

herr_t
H5LTfind_dataset( hid_t loc_id, const char *dset_name )
{
    return H5Literate(loc_id, H5_INDEX_NAME, H5_ITER_INC, 0, find_dataset, (void *)dset_name);
}
  
  
  
  
