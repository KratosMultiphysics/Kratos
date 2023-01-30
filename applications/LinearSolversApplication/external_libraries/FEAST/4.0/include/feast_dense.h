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


/////////////////////////////////////////////////////////////////////////////////////////////////////
// External declarations
/////////////////////////////////////////////////////////////////////////////////////////////////////


//FEAST Hermitian and Non-Hermitian and Non-linear (polynomial)

// Double precision
extern void dfeast_sygvx_(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void dfeast_sygv_(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_syevx_(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void dfeast_syev_(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hegvx_(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void zfeast_hegv_(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_heevx_(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void zfeast_heev_(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gegvx_(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gegv_(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_geevx_(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_geev_(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_sygvx_(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_sygv_(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_syevx_(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_syev_(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gegvx_(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gegv_(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_geevx_(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_geev_(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

extern void dfeast_sypevx_(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void dfeast_sypev_(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hepevx_(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void zfeast_hepev_(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gepevx_(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gepev_(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_sypevx_(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_sypev_(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gepevx_(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gepev_(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);




/////////////////////////////////////////////////////////////////////////////////////////////////////
// FEAST interfaces
/////////////////////////////////////////////////////////////////////////////////////////////////////

// Double precision
void dfeast_sygvx(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  dfeast_sygvx_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_sygv(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sygv_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void dfeast_syevx(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  dfeast_syevx_(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_syev(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_syev_(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hegvx(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  zfeast_hegvx_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hegv(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hegv_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_heevx(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  zfeast_heevx_(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_heev(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_heev_(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_gegvx(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_gegvx_(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gegv(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_gegv_(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_geevx(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_geevx_(N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_geev(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_geev_(N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_sygvx(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_sygvx_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_sygv(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_sygv_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_syevx(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_syevx_(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_syev(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_syev_(UPLO,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gegvx(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_gegvx_(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gegv(int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_gegv_(N,A,LDA,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_geevx(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_geevx_(N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_geev(int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_geev_(N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void DFEAST_SYGV(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sygv_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SYEV(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_syev_(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HEGV(char *UPLO,int *N,double *A,int *LDA,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hegv_(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HEEV(char *UPLO,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_heev_(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}


void dfeast_sypevx(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  dfeast_sypevx_(UPLO,d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_sypev(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sypev_(UPLO,d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void zfeast_hepevx(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  zfeast_hepevx_(UPLO,d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hepev(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hepev_(UPLO,d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_gepevx(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_gepevx_(d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gepev(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_gepev_(d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_sypevx(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_sypevx_(UPLO,d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_sypev(char *UPLO,int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_sypev_(UPLO,d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gepevx(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_gepevx_(d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gepev(int *d,int *N,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_gepev_(d,N,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
