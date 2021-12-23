/*
!=========================================================================================
!Copyright (c) 2009-2019, The Regents of the University of Massachusetts, Amherst.
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// External declarations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern void feastinit_(int *fpm);
extern void feastinit_driver_(int *fpm,int *N);
extern void cfeast_customcontour_(int *fpm2,int *N,int *Nedge,int *Tedge,float *Zedge,float *Zne,float *Wne);
extern void zfeast_customcontour_(int *fpm2,int *N,int *Nedge,int *Tedge,double *Zedge,double *Zne,double *Wne);
extern void zfeast_contour_(double *Emin, double *Emax, int *fpm2, int *fpm16, int *fpm18, double *Zne, double *Wne);
extern void cfeast_contour_(float *Emin, float *Emax, int *fpm2, int *fpm16, int *fpm18, float *Zne, float *Wne);
extern void zfeast_gcontour_(double *Emid, double *r, int *fpm2, int *fpm17, int *fpm19, double *Zne, double *Wne);
extern void dfeast_rational_(double *Emin, double *Emax, int *fpm2, int *fpm16, int *fpm18,double *Eig, int *M0,double *f);
extern void sfeast_rational_(float *Emin,float *Emax,int *fpm2,int *fpm16,int *fpm18,float *Eig,int *M0,float *f);
extern void sfeast_rationalx_(float *Zne,float *Wne,int *fpm2,float *Eig,int *M0,float *f);
extern void zfeast_grational_(double *Emid,double *r,int *fpm2,int *fpm17,int *fpm19,double *Eig,int *M0,double *f);
extern void cfeast_grational_(float *Emid,float*r,int *fpm2,int *fpm17,int *fpm19,float *Eig,int *M0,float*f);
extern void zfeast_grationalx_(double *Zne,double *Wne,int *fpm2,double *Eig,int *M0,double *f);
extern void cfeast_grationalx_(float *Zne,float *Wne,int *fpm2,float*Eig,int *M0,float *f);
extern void dfeast_rationalx_(double *Zne,double *Wne,int *fpm2,double *Eig,int *M0,double *f);
extern void cfeast_gcontour_(float *Emid, float *r, int *fpm2, int *fpm17, int *fpm19, float *Zne, float *Wne);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FEAST interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void feastinit(int *fpm) {
     feastinit_(fpm);
}
void FEASTINIT(int *fpm) {
     feastinit_(fpm);
}
void feast_init_driver(int *fpm,int *N) {
     feastinit_driver_(fpm,N);
}
void cfeast_customcontour(int *fpm2,int *N,int *Nedge,int *Tedge,float *Zedge,float *Zne,float *Wne){
     cfeast_customcontour_(fpm2,N,Nedge,Tedge,Zedge,Zne,Wne);
}
void zfeast_customcontour(int *fpm2,int *N,int *Nedge,int *Tedge,double *Zedge,double *Zne,double *Wne){
     zfeast_customcontour_(fpm2,N,Nedge,Tedge,Zedge,Zne,Wne);
}
void zfeast_contour(double *Emin, double *Emax, int *fpm2, int *fpm16, int *fpm18, double *Zne, double *Wne){
     zfeast_contour_(Emin,Emax,fpm2,fpm16,fpm18,Zne,Wne);
}
void cfeast_contour(float *Emin,float *Emax, int *fpm2, int *fpm16, int *fpm18, float *Zne, float *Wne){
     cfeast_contour_(Emin,Emax,fpm2,fpm16,fpm18,Zne,Wne);
}
void zfeast_gcontour(double *Emid, double *r, int *fpm2, int *fpm17, int *fpm19, double *Zne, double *Wne){
     zfeast_gcontour_(Emid,r,fpm2,fpm17,fpm19,Zne,Wne);
}
void cfeast_gcontour(float *Emid, float *r, int *fpm2, int *fpm17, int *fpm19, float *Zne, float *Wne){
     cfeast_gcontour_(Emid,r,fpm2,fpm17,fpm19,Zne,Wne);
}
void dfeast_rational(double *Emin, double *Emax, int *fpm2, int *fpm16, int *fpm18,double *Eig, int *M0,double *f){
     dfeast_rational_(Emin,Emax,fpm2,fpm16,fpm18,Eig,M0,f);
}
void sfeast_rational(float *Emin,float *Emax,int *fpm2,int *fpm16,int *fpm18,float *Eig,int *M0,float *f){
     sfeast_rational_(Emin,Emax,fpm2,fpm16,fpm18,Eig,M0,f);
}
void dfeast_rationalx(double *Zne,double *Wne,int *fpm2,double *Eig,int *M0,double *f){
            dfeast_rationalx_(Zne,Wne,fpm2,Eig,M0,f);
}
void sfeast_rationalx(float *Zne,float *Wne,int *fpm2,float *Eig,int *M0,float *f){
     sfeast_rationalx_(Zne,Wne,fpm2,Eig,M0,f);
}
void zfeast_grational(double *Emid,double *r,int *fpm2,int *fpm17,int *fpm19,double *Eig,int *M0,double *f){
     zfeast_grational_(Emid,r,fpm2,fpm17,fpm19,Eig,M0,f);
}
void cfeast_grational(float *Emid,float*r,int *fpm2,int *fpm17,int *fpm19,float *Eig,int *M0,float*f){
     cfeast_grational_(Emid,r,fpm2,fpm17,fpm19,Eig,M0,f);
}
void zfeast_grationalx(double *Zne,double *Wne,int *fpm2,double *Eig,int *M0,double *f){
     zfeast_grationalx_(Zne,Wne,fpm2,Eig,M0,f);
}
void cfeast_grationalx(float *Zne,float *Wne,int *fpm2,float*Eig,int *M0,float *f){
     cfeast_grationalx_(Zne,Wne,fpm2,Eig,M0,f);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

















