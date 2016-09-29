/*
!=========================================================================================
!Copyright (c) 2009-2011, The Regents of the University of Massachusetts, Amherst.
!Developed by E. Polizzi
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification, 
!are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list of conditions 
!   and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
!   and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
!    products derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
!ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!==========================================================================================
*/

extern double wdsin_(double *x);
extern double wdcos_(double *x);
extern float wssin_(float *x);
extern float wscos_(float *x);

extern void wallocate_3i_(int *A, int *firstNo, int *secondNo, int *thirdNo, int *alloc_info); 
/*extern void wallocate_2i_(int *A, int *firstNo, int *secondNo,  int *alloc_info); 
extern void wallocate_1i_(int *A, int *firstNo, int *alloc_info); 
extern void wallocate_3d_(double *A, int *firstNo, int *secondNo, int *thirdNo, int *alloc_info); 
extern void wallocate_2d_(double *A, int *firstNo, int *secondNo,  int *alloc_info); 
extern void wallocate_1d_(double *A, int *firstNo, int *alloc_info); 
extern void wallocate_3z_(double *A, int *firstNo, int *secondNo, int *thirdNo, int *alloc_info); 
extern void wallocate_2z_(double *A, int *firstNo, int *secondNo,  int *alloc_info); 
extern void wallocate_1z_(double *A, int *firstNo, int *alloc_info); 
extern void wallocate_3s_(float *A, int *firstNo, int *secondNo, int *thirdNo, int *alloc_info); 
extern void wallocate_2s_(float *A, int *firstNo, int *secondNo,  int *alloc_info); 
extern void wallocate_1s_(float *A, int *firstNo, int *alloc_info); 
extern void wallocate_3c_(float *A, int *firstNo, int *secondNo, int *thirdNo, int *alloc_info); 
extern void wallocate_2c_(float *A, int *firstNo, int *secondNo,  int *alloc_info); 
extern void wallocate_1c_(float *A, int *firstNo, int *alloc_info); 
extern void wdeallocate_3i_(int *A);    
extern void wdeallocate_2i_(int *A);    
extern void wdeallocate_1i_(int *A);    
extern void wdeallocate_3d_(double *A);    
extern void wdeallocate_2d_(double *A);    
extern void wdeallocate_1d_(double *A);    
extern void wdeallocate_3z_(double *A);    
extern void wdeallocate_2z_(double *A);    
extern void wdeallocate_1z_(double *A);   
extern void wdeallocate_3s_(float *A);    
extern void wdeallocate_2s_(float *A);    
extern void wdeallocate_1s_(float *A);    
extern void wdeallocate_3c_(float *A);    
extern void wdeallocate_2c_(float *A);    
extern void wdeallocate_1c_(float *A);    
extern void wwrite_n_(long int *file_id);       
extern void wwrite_t_(long int *file_id);       
extern void wwrite_s_(long int *file_id, void *buffer);       
extern void wwrite_i_(long int *file_id, int *buffer);       
extern void wwrite_f_(long int *file_id, float *buffer);       
extern void wwrite_d_(long int *file_id, double *buffer);   */    



double wdsin(double *x){
  double wdsin_(double *x);}

double wdcos(double *x){
  double wdcos_(double *x);}

float wssin(float *x){
  float wssin_(float *x);}

float wscos(float *x){
  float wscos_(float *x);}

void wallocate_3i(int *A, int *firstNo, int *secondNo, int *thirdNo, int *alloc_info){
wallocate_3i_(A, firstNo, secondNo, thirdNo, alloc_info);} 

