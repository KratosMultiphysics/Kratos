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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// External declarations
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Double precision
extern void dfeast_sbgvx_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void dfeast_sbgv_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_sbevx_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void dfeast_sbev_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hbgvx_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void zfeast_hbgv_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hbevx_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void zfeast_hbev_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gbgvx_(int *N,int *kla,int *kua,double *A,int *LDA,int *klb,int *kub,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gbgv_(int *N,int *kla,int *kua,double *A,int *LDA,int *klb, int *kub, double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gbevx_(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gbev_(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_sbgvx_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_sbgv_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_sbevx_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_sbev_(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gbgvx_(int *N,int *kla,int *kua,double *A,int *LDA,int *klb,int *kub,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gbgv_(int *N,int *kla,int *kua,double *A,int *LDA,int *klb, int *kub, double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gbevx_(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gbev_(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// FEAST interfaces
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Double precision
void dfeast_sbgvx(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_sbgvx_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_sbgv(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sbgv_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void dfeast_sbevx(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_sbevx_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_sbev(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sbev_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hbgvx(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_hbgvx_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hbgv(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hbgv_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hbevx(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_hbevx_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hbev(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hbev_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SBGV(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sbgv_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SBEV(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_sbev_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HBGV(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hbgv_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HBEV(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hbev_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_gbgvx(int *N,int *kla,int *kua,double *A,int *LDA,int *klb,int *kub,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_gbgvx_(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gbgv(int *N,int *kla,int *kua,double *A,int *LDA,int *klb,int *kub,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_gbgv_(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_gbevx(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_gbevx_(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gbev(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_gbev_(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_sbgvx(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_sbgvx_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_sbgv(char *UPLO,int *N,int *kla,double *A,int *LDA,int *klb,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_sbgv_(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_sbevx(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_sbevx_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_sbev(char *UPLO,int *N,int *kla,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_sbev_(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gbgvx(int *N,int *kla,int *kua,double *A,int *LDA,int *klb,int *kub,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_gbgvx_(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gbgv(int *N,int *kla,int *kua,double *A,int *LDA,int *klb,int *kub,double *B,int *LDB,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_gbgv_(N,kla,kua,A,LDA,klb,kub,B,LDB,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gbevx(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_gbevx_(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gbev(int *N,int *kla,int *kua,double *A,int *LDA,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_gbev_(N,kla,kua,A,LDA,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
