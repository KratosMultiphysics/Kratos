/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FILE created  by Huiying Xu 2007
!!      modified by Eric Polizzi 2009 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/


#include "c_functions_wrapper.h"


double wdatan2_(double *x,double *y){
return atan2(*x,*y);
}

double wdsin_(double *x){
return sin(*x);
}

double wdcos_(double *x){
return cos(*x);
}

float wsatan2_(float *x,float *y){
return atan2(*x,*y);
}

float wssin_(float *x){
return sin(*x);
}

float wscos_(float *x){
return cos(*x);
}





void wdeallocate_all_type(array_all_dim *p) {

/*
 * This function deallocates the space of an array.
 * It works for all dimensions of arrays.
 */

        long int i, j;
/*
        if (p->base != (void *) 0) free(p->base);
*/
        free(p->base);
        p->base = (void *) 0;

	return;
   
        j = p->ndim;
  
        p->size = 0;
        p->offset = (void *) 0;
        p->bitsets = 0;
        p->ndim = 0;
        p->reserved = 0;
  
        for (i = 0; i < j; i++)
                memset(p->dim_info+i, 0, sizeof(dimension_info));
}


void wdeallocate_1i_(array_dim1 *p) {

/*
 * This function deallocates the space of one-dimensional integer array.
 */

	array_all_dim *p1 = (array_all_dim *) p;

	wdeallocate_all_type(p1);

}
 
 
void wdeallocate_2i_(array_dim2 *p) {

/*
 * This function deallocates the space of two-dimensional integer array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);

}


void wdeallocate_3i_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional integer array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
  
        wdeallocate_all_type(p1);
 
}


void wdeallocate_1d_(array_dim1 *p) {

/*
 * This function deallocates the space of one-dimensional double array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2d_(array_dim2 *p) {
   
/*
 * This function deallocates the space of two-dimensional double array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);
}

void wdeallocate_3d_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional double array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
            
        wdeallocate_all_type(p1);
}



void wdeallocate_1z_(array_dim1 *p) {

/*
 * This function deallocates the space of one-dimensional double complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2z_(array_dim2 *p) {
   
/*
 * This function deallocates the space of two-dimensional double complex array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);
}

void wdeallocate_3z_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional double complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
            
        wdeallocate_all_type(p1);
}



void wdeallocate_1s_(array_dim1 *p) {

/*
 * This function deallocates the space of one-dimensional real array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2s_(array_dim2 *p) {
   
/*
 * This function deallocates the space of two-dimensional real array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);
}

void wdeallocate_3s_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional real array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
            
        wdeallocate_all_type(p1);
}



void wdeallocate_1c_(array_dim1 *p) {

/*
 * This function deallocates the space of one-dimensional complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
 
	wdeallocate_all_type(p1);
	
}

void wdeallocate_2c_(array_dim2 *p) {
   
/*
 * This function deallocates the space of two-dimensional complex array.
 */
        
	array_all_dim *p1 = (array_all_dim *) p;

        wdeallocate_all_type(p1);
}

void wdeallocate_3c_(array_dim3 *p) {

/*
 * This function deallocates the space of three-dimensional complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
            
        wdeallocate_all_type(p1);
}





void wallocate_all_type(array_all_dim *p, int * firstNo, int *secondNo, int *thirdNo, int type, int *alloc_info) {

/*
 * This function allocates the space for an array.
 * It works for all dimensions of arrays.
 */

        dimension_info dim;
        void *base = NULL;
        /* int element_size; */
	size_t element_size;

	switch(type) {
	    case F90_INTEGER:
		element_size = sizeof(int);
		break;
	    case F90_REAL:
		element_size = sizeof(float);
		break;
            case F90_COMPLEX:
		element_size = 2*sizeof(float);
		break;
	    case F90_DOUBLE_COMPLEX:
		element_size = 2*sizeof(double);
		break;
	    default:   /*case F90_DOUBLE_PRECISION: */
		element_size = sizeof(double);
	}


	if ((*thirdNo) != 0) {
		base = (void *) malloc(element_size * (*firstNo) * (*secondNo) * (*thirdNo));
		if (base == NULL) {
			*alloc_info = ALLOC_FAIL;
			return;
		}
		else
			*alloc_info = ALLOC_SUCCESS;
		p->ndim = 3;
	}
	else if ((*secondNo) != 0) {
		base = (void *) malloc(element_size * (*firstNo) * (*secondNo));
                if (base == NULL) {
                        *alloc_info = ALLOC_FAIL;
			return;
		}
                else
                        *alloc_info = ALLOC_SUCCESS;
                p->ndim = 2;
        }
	else {
                base = (void *) malloc((*firstNo) * element_size);
                if (base == NULL) {
                        *alloc_info = ALLOC_FAIL;
			return;
		}
                else
                        *alloc_info = ALLOC_SUCCESS;
                p->ndim = 1;
	}


        p->base = base;
        p->size = (long int) element_size;
        p->offset = 0;
        p->bitsets = 5;
        p->reserved = 0;

        dim.nelts = (long int)(*firstNo);
        dim.stride = (long int) element_size;
        dim.lower_bound = 1;
        memcpy(p->dim_info, &dim, sizeof(dimension_info));

	if ((*secondNo) != 0) {
		dim.nelts = (long int)(*secondNo);
        	dim.stride = (long int)((*firstNo) * element_size);
        	dim.lower_bound = 1;
		memcpy(p->dim_info+1, &dim, sizeof(dimension_info));
	}

        if ((*thirdNo) != 0) {
                dim.nelts = (long int) (*thirdNo);
                dim.stride = (long int) ((*firstNo) * (*secondNo) * element_size);
                dim.lower_bound = 1;
                memcpy(p->dim_info+2, &dim, sizeof(dimension_info));
        }
}




void wallocate_1i_(array_dim1 *p, int * element_No, int *alloc_info) {

/*
 * This function allocates the space for one-dimensional integer array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *)p;

	wallocate_all_type(p1, element_No, &i, &j, F90_INTEGER, alloc_info);

}


 
void wallocate_2i_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {
	
/*
 * This function allocates the space for two-dimensional integer array.
 */

	int j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, rowNo, columNo, &j, F90_INTEGER, alloc_info);

}


void wallocate_3i_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) {

/*
 * This function allocates the space for three-dimensional integer array.
 */
      
	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_INTEGER, alloc_info);
}




void wallocate_1d_(array_dim1 *p, int * element_No, int *alloc_info) {

/*
 * This function allocates the space for one-dimensional double array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_DOUBLE_PRECISION, alloc_info);

}



void wallocate_2d_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional double array.
 */

	int j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_DOUBLE_PRECISION, alloc_info);

}


void wallocate_3d_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) {

/*
 * This function allocates the space for three-dimensional double array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_DOUBLE_PRECISION, alloc_info);
 
}


void wallocate_1z_(array_dim1 *p, int * element_No, int *alloc_info) {

/*
 * This function allocates the space for one-dimensional double array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_DOUBLE_COMPLEX, alloc_info);

}



void wallocate_2z_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional double array.
 */

	int j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_DOUBLE_COMPLEX, alloc_info);

}


void wallocate_3z_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) {

/*
 * This function allocates the space for three-dimensional double array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_DOUBLE_COMPLEX, alloc_info);
 
}




void wallocate_1s_(array_dim1 *p, int * element_No, int *alloc_info) {

/*
 * This function allocates the space for one-dimensional real array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_REAL, alloc_info);

}



void wallocate_2s_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional real array.
 */

	int j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_REAL, alloc_info);

}


void wallocate_3s_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) {

/*
 * This function allocates the space for three-dimensional real array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_REAL, alloc_info);
 
}


void wallocate_1c_(array_dim1 *p, int * element_No, int *alloc_info) {

/*
 * This function allocates the space for one-dimensional complex array.
 */

	int i = 0, j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
 
	wallocate_all_type(p1, element_No, &i, &j, F90_COMPLEX, alloc_info);

}



void wallocate_2c_(array_dim2 *p, int * rowNo, int *columNo, int *alloc_info) {

/*
 * This function allocates the space for two-dimensional complex array.
 */

	int j = 0;
	array_all_dim *p1 = (array_all_dim *) p;
	
	wallocate_all_type(p1, rowNo, columNo, &j, F90_COMPLEX, alloc_info);

}


void wallocate_3c_(array_dim3 *p, int * firstNo, int *secondNo, int *thirdNo, int *alloc_info) {

/*
 * This function allocates the space for three-dimensional complex array.
 */

	array_all_dim *p1 = (array_all_dim *) p;
        wallocate_all_type(p1, firstNo, secondNo, thirdNo, F90_COMPLEX, alloc_info);
 
}








void wwrite_n_(unsigned long int *file_id) {
	int flag = 0;
        if (*file_id == (unsigned long int) 6) flag = -1;

		 
	if (flag == 0) {fflush((FILE*) *file_id); fprintf((FILE *) *file_id, "\n");}
			else {fflush(stdout); printf("\n");}
}
 



void wwrite_t_(unsigned long int *file_id){

	int flag = 0;

        if (*file_id == (unsigned long int) 6) flag = -1;

	
	if (flag == 0) {fflush((FILE *) *file_id);fprintf((FILE *) *file_id, "\t");}
			else {fflush(stdout); printf("\t");}
}




void wwrite_i_(unsigned long int *file_id, void *buffer){

	int flag = 0;
        if (*file_id == (unsigned long int) 6) flag = -1;
	if (flag == 0) {fflush((FILE *) *file_id);fprintf((FILE *) *file_id, "%d", * (int *)buffer);}
		else {fflush(stdout); printf("%d", * (int *)buffer);}
}
	       

void wwrite_s_(unsigned long int *file_id, void *buffer){

	int flag = 0;
if (*file_id == (unsigned long int) 6) flag = -1;
 if (flag == 0) {fflush((FILE *) *file_id); fprintf((FILE *) *file_id, "%s", (char *) buffer);}
		else {fflush(stdout); printf("%s", (char *) buffer);}
	
}



void wwrite_f_(unsigned long int *file_id, void *buffer){

	int flag = 0;
if (*file_id == (unsigned long int) 6) flag = -1;
 if (flag == 0) {fflush((FILE *) *file_id); fprintf((FILE *) *file_id, "%f", * (float *)buffer);}
		else {fflush(stdout); printf("%.7e", * (float *)buffer);}
}

void wwrite_d_(unsigned long int *file_id, void *buffer){

	int flag = 0;
if (*file_id == (unsigned long int) 6) flag = -1;
 if (flag == 0) {fflush((FILE *) *file_id); fprintf((FILE *) *file_id, "%e", * (double *)buffer);}
		else {fflush(stdout); printf("%.15e", * (double *)buffer);}
}



