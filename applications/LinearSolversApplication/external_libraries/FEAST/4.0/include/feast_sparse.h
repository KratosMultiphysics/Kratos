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

////////////////////////////////////////////////////////////////////////////////////////////////////
// External declarations
////////////////////////////////////////////////////////////////////////////////////////////////////

// FEAST Double precision
extern void dfeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hcsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_hcsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hcsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_hcsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

extern void difeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void difeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void difeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void difeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_hcsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_hcsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_hcsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_hcsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void difeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void difeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void difeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void difeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);


// non-linear polynomial
extern void dfeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hcsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_hcsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zfeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void dfeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void dfeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

extern void difeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void difeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_hcsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_hcsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zifeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void zifeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void difeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void difeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);





////////////////////////////////////////////////////////////////////////////////////////////////////
// FEAST interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////

// FEAST Double precision

void dfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void dfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hcsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_hcsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hcsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hcsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_hcsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hcsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DFEAST_SCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zfeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void dfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}


// IFEAST Double precision

void difeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  difeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void difeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  difeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void difeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  difeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void difeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  difeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zifeast_hcsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zifeast_hcsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_hcsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zifeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zifeast_hcsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zifeast_hcsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_hcsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zifeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DIFEAST_SCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  difeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void DIFEAST_SCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  difeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZIFEAST_HCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zifeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZIFEAST_HCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zifeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zifeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zifeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zifeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zifeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zifeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zifeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zifeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zifeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zifeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zifeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     zifeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zifeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void difeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     difeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void difeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     difeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void difeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     difeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void difeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     difeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}



// Non-linear Polynomial

// FEAST Double precision

void dfeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void zfeast_hcsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_hcsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hcsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_hcsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void zfeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void zfeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zfeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void dfeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  dfeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void dfeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  dfeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

// IFEAST Double precision

void difeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  difeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void difeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  difeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void zifeast_hcsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zifeast_hcsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_hcsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zifeast_hcsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void zifeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zifeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zifeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void zifeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  zifeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zifeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zifeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void difeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  difeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void difeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  difeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
