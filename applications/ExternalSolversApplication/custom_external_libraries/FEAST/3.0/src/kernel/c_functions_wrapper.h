/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FILE created  by Huiying Xu 2007
!!      modified by Eric Polizzi 2009 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#define _GNU_SOURCE
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>


#define F90_WRITE_NOTATION      -2 
#define F90_LOGICAL             -1 
#define F90_CHAR                0 
#define F90_STRING              1 
#define F90_INTEGER             2 
#define F90_REAL                3 
#define F90_DOUBLE_PRECISION    4 
#define F90_INTEGER_STAR_8      5 
#define F90_DOUBLE_COMPLEX      8
#define F90_COMPLEX             9 


#define ALLOC_SUCCESS	0
#define ALLOC_FAIL	-1
 
 typedef struct {  
   long int nelts;  
   long int stride;  
   long int lower_bound;  
 } dimension_info;    
   
 typedef struct {  
   void * base;    
   long int size;     
   void * offset;  
   long int bitsets;   
   long int ndim;  
   long int reserved;  
   dimension_info dim_info[1];  
 } array_dim1;  
 
 typedef struct {  
   void * base;  
   long int size;    
   void * offset;  
   long int bitsets;  
   long int ndim;  
   long int reserved;  
   dimension_info dim_info[2];  
 } array_dim2;  
   
   
typedef struct {
  void * base;
  long int size;
  void * offset;
  long int bitsets;
  long int ndim;
  long int reserved;
  dimension_info dim_info[3];
} array_dim3;
 
 
typedef struct {
  void * base;
  long int size;
  void * offset;
  long int bitsets;
  long int ndim;
  long int reserved;
  dimension_info dim_info[];
} array_all_dim;

//void ess_(int *i);

double wdatan2_(double *x,double *y);
double wdsin_(double *x);
double wdcos_(double *x);
float wsatan2_(float *x, float *y);
float wssin_(float *x);
float wscos_(float *x);

//void wstop_();
//void char2int_(char *s, int *value);

void wdeallocate_all_type(array_all_dim *p);
void wdeallocate_1i_(array_dim1 *p);
void wdeallocate_2i_(array_dim2 *p);
void wdeallocate_3i_(array_dim3 *p);
void wdeallocate_1d_(array_dim1 *p);
void wdeallocate_2d_(array_dim2 *p);
void wdeallocate_3d_(array_dim3 *p);
void wdeallocate_1z_(array_dim1 *p);
void wdeallocate_2z_(array_dim2 *p);
void wdeallocate_3z_(array_dim3 *p);
void wdeallocate_1s_(array_dim1 *p);
void wdeallocate_2s_(array_dim2 *p);
void wdeallocate_3s_(array_dim3 *p);
void wdeallocate_1c_(array_dim1 *p);
void wdeallocate_2c_(array_dim2 *p);
void wdeallocate_3c_(array_dim3 *p);


void wallocate_all_type(array_all_dim *p, int * firstNo, int *secondNo, int *thirdNo, int type, int *alloc_info);
void wallocate_1i_(array_dim1 *p, int * element_No, int *alloc_info);
void wallocate_2i_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info);
void wallocate_3i_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info);
void wallocate_1d_(array_dim1 *p, int * element_No, int *alloc_info);
void wallocate_2d_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info);
void wallocate_3d_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info);
void wallocate_1z_(array_dim1 *p, int * element_No, int *alloc_info);
void wallocate_2z_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info);
void wallocate_3z_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info);
void wallocate_1s_(array_dim1 *p, int * element_No, int *alloc_info);
void wallocate_2s_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info);
void wallocate_3s_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info);
void wallocate_1c_(array_dim1 *p, int * element_No, int *alloc_info);
void wallocate_2c_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info);
void wallocate_3c_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info);


//void wwrite_(unsigned long int *file_id, void *buffer, int *type);

void wwrite_n_(unsigned long int *file_id);
void wwrite_t_(unsigned long int *file_id);
void wwrite_i_(unsigned long int *file_id, void *buffer);
void wwrite_s_(unsigned long int *file_id, void *buffer);
void wwrite_f_(unsigned long int *file_id, void *buffer);
void wwrite_d_(unsigned long int *file_id, void *buffer);






