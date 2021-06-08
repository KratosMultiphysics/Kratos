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


#include "feast_sparse.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// External declarations
////////////////////////////////////////////////////////////////////////////////////////////////////

// PFEAST Double precision
extern void pdfeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdfeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdfeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdfeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_hcsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_hcsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_hcsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_hcsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdfeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdfeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdfeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdfeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

// PIFEAST Double precision
extern void pdifeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdifeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdifeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdifeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_hcsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_hcsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_hcsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_hcsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_scsrgvx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_scsrgv_(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_scsrevx_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_scsrev_(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdifeast_gcsrgvx_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdifeast_gcsrgv_(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdifeast_gcsrevx_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdifeast_gcsrev_(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);


// non-linear polynomial
extern void pdfeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdfeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_hcsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_hcsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzfeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzfeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdfeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdfeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

extern void pdifeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdifeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_hcsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_hcsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pzifeast_scsrpevx_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pzifeast_scsrpev_(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void pdifeast_gcsrpevx_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne);
extern void pdifeast_gcsrpev_(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

////////////////////////////////////////////////////////////////////////////////////////////////////
// FEAST interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////

// PFEAST Double precision

void pdfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pdfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdfeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdfeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pzfeast_hcsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzfeast_hcsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_hcsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pzfeast_hcsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzfeast_hcsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_hcsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzfeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PDFEAST_SCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PDFEAST_SCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdfeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PZFEAST_HCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzfeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PZFEAST_HCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzfeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pzfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pzfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzfeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzfeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pzfeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzfeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzfeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pzfeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzfeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzfeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pdfeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pdfeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdfeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pdfeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pdfeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pdfeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdfeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pdfeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}



// PIFEAST Double precision

void pdifeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdifeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdifeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdifeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pdifeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdifeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdifeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdifeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pzifeast_hcsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzifeast_hcsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_hcsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzifeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pzifeast_hcsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzifeast_hcsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_hcsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzifeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PDIFEAST_SCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdifeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PDIFEAST_SCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdifeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PZIFEAST_HCSRGV(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzifeast_hcsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void PZIFEAST_HCSREV(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzifeast_hcsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void pzifeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzifeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzifeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pzifeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzifeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzifeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pzifeast_scsrgvx(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzifeast_scsrgvx_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_scsrgv(char *UPLO,int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzifeast_scsrgv_(UPLO,N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pzifeast_scsrevx(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pzifeast_scsrevx_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_scsrev(char *UPLO,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pzifeast_scsrev_(UPLO,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pdifeast_gcsrgvx(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pdifeast_gcsrgvx_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdifeast_gcsrgv(int *N,double *sa,int *isa,int *jsa,double *sb,int *isb,int *jsb,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pdifeast_gcsrgv_(N,sa,isa,jsa,sb,isb,jsb,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void pdifeast_gcsrevx(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     pdifeast_gcsrevx_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdifeast_gcsrev(int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     pdifeast_gcsrev_(N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}



// Non-linear Polynomial

// PFEAST Double precision

void pdfeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdfeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdfeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdfeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pzfeast_hcsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzfeast_hcsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_hcsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzfeast_hcsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pzfeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzfeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzfeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pzfeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzfeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzfeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzfeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pdfeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdfeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdfeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdfeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

// IFEAST Double precision

void pdifeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdifeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdifeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdifeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pzifeast_hcsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzifeast_hcsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_hcsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzifeast_hcsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pzifeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzifeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzifeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pzifeast_scsrpevx(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pzifeast_scsrpevx_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pzifeast_scsrpev(char *UPLO,int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pzifeast_scsrpev_(UPLO,d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}

void pdifeast_gcsrpevx(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
  pdifeast_gcsrpevx_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void pdifeast_gcsrpev(int *d,int *N,double *sa,int *isa,int *jsa,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  pdifeast_gcsrpev_(d,N,sa,isa,jsa,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
