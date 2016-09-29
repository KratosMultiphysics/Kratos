/*
!=========================================================================================
!Copyright (c) 2009-2015, The Regents of the University of Massachusetts, Amherst.
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

////////////////////////////////////////////////////////////////////////////////////////////////////
// External declarations
////////////////////////////////////////////////////////////////////////////////////////////////////

// Double precision
extern void dfeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hcsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_hcsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hcsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_hcsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

// Single precision
extern void sfeast_scsrgvx_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void sfeast_scsrgv_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void sfeast_scsrevx_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void sfeast_scsrev_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_hcsrgvx_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void cfeast_hcsrgv_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_hcsrevx_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void cfeast_hcsrev_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_gcsrgvx_(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void cfeast_gcsrgv_(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_gcsrevx_(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void cfeast_gcsrev_(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_scsrgvx_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void cfeast_scsrgv_(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void cfeast_scsrevx_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void cfeast_scsrev_(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void sfeast_gcsrgvx_(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void sfeast_gcsrgv_(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info);
extern void sfeast_gcsrevx_(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne);
extern void sfeast_gcsrev_(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info);




////////////////////////////////////////////////////////////////////////////////////////////////////
// FEAST interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////

// Double precision

void dfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void dfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_scsrevx_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hcsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_hcsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hcsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hcsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_hcsrevx_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hcsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_gcsrevx_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_gcsrev_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_scsrevx_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_scsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_gcsrevx_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *feastparam,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_gcsrev_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}


// Single precision

void sfeast_scsrgvx(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
  sfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void sfeast_scsrgv(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void sfeast_scsrevx(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
  sfeast_scsrevx_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void sfeast_scsrev(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_scsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void cfeast_hcsrgvx(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
  cfeast_hcsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void cfeast_hcsrgv(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void cfeast_hcsrevx(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
  cfeast_hcsrevx_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void cfeast_hcsrev(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_hcsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void SFEAST_SCSRGV(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void SFEAST_SCSREV(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  sfeast_scsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void CFEAST_HCSRGV(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void CFEAST_HCSREV(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emin,float *Emax,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
  cfeast_hcsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void cfeast_gcsrgvx(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
     cfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void cfeast_gcsrgv(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
     cfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void cfeast_gcsrevx(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
     cfeast_gcsrevx_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void cfeast_gcsrev(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
     cfeast_gcsrev_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void cfeast_scsrgvx(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
     cfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void cfeast_scsrgv(char *UPLO,int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
     cfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void cfeast_scsrevx(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
     cfeast_scsrevx_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void cfeast_scsrev(char *UPLO,int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
     cfeast_scsrev_(UPLO,N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void sfeast_gcsrgvx(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
     sfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void sfeast_gcsrgv(int *N,float *sa,int *isa,int *jsa,float *sb,int *isb,int *jsb,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
     sfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void sfeast_gcsrevx(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info,float *Zne,float *Wne){
     sfeast_gcsrevx_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void sfeast_gcsrev(int *N,float *sa,int *isa,int *jsa,int *feastparam,float *epsout,int *loop,float *Emid,float *r,int *M0,float *lambda,float *q,int *mode,float *res,int *info){
     sfeast_gcsrev_(N,sa,isa,jsa,feastparam,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}


