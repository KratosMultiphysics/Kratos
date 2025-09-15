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


#include "feast_tools.h"


//////////////////////////////////////////////////////////////////////////
// External function delarations
//////////////////////////////////////////////////////////////////////////

// Double precision
extern void dfeast_srci_(char *ijob,int *N,double *Ze,double *work,double   *workc,double   *Aq,double   *Sq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_hrci_(char *ijob,int *N,double   *Ze,double   *work,double   *workc,double   *zAq,double   *zSq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double   *q,int *mode,double *res,int *info);
extern void dfeast_srcix_(char *ijob,int *N,double *Ze,double *work,double   *workc,double   *Aq,double   *Sq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void zfeast_hrcix_(char *ijob,int *N,double   *Ze,double   *work,double   *workc,double   *zAq,double   *zSq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double   *q,int *mode,double *res,int *info,double *Zne, double *Wne);
extern void zfeast_srci_(char *ijob,int *N,double *Ze,double *work ,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,int *mode,double *resr,int *info);
extern void dfeast_grci_(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res, int *info);
extern void zfeast_grci_(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);
extern void zfeast_srcix_(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,int *mode,double *resr,int *info,double *Zne, double *Wne);
extern void dfeast_grcix_(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info, double *Zne, double *Wne);
extern void zfeast_grcix_(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info, double *Zne, double *Wne);

extern void zfeast_grcipevx_(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info, double *Zne, double *Wne);
extern void zfeast_grcipev_(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);

extern void zfeast_srcipevx_(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info, double *Zne, double *Wne);
extern void zfeast_srcipev_(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info);







//////////////////////////////////////////////////////////////////////////
// FEAST interfaces
//////////////////////////////////////////////////////////////////////////

// Double precision

void dfeast_srci(char *ijob,int *N,double   *Ze,double *work,double   *workc,double   *Aq,double   *Sq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_srci_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void dfeast_srcix(char *ijob,int *N,double   *Ze,double *work,double   *workc,double   *Aq,double   *Sq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
     dfeast_srcix_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_hrci(char *ijob,int *N,double   *Ze,double   *work,double   *workc,double   *zAq,double   *zSq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double   *q,int *mode,double *res,int *info){
    zfeast_hrci_(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void zfeast_hrcix(char *ijob,int *N,double   *Ze,double   *work,double   *workc,double   *zAq,double   *zSq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double   *q,int *mode,double *res,int *info,double *Zne, double *Wne){
    zfeast_hrcix_(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_srci(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,int *mode,double *resr,int *info){
     zfeast_srci_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info);
}
void dfeast_grci(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_grci_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_grci(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     zfeast_grci_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_srcix(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *qr,int *mode,double *resr,int *info,double *Zne,double *Wne){
     zfeast_srcix_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,qr,mode,resr,info,Zne,Wne);
}
void dfeast_grcix(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne,double *Wne){
     dfeast_grcix_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_grcix(char *ijob,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
     zfeast_grcix_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void DFEAST_SRCI(char *ijob,int *N,double   *Ze,double *work,double   *workc,double   *Aq,double   *Sq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
     dfeast_srci_(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}
void ZFEAST_HRCI(char *ijob,int *N,double   *Ze,double   *work,double   *workc,double   *zAq,double   *zSq,int *fpm,double *epsout,int *loop,double *Emin,double *Emax,int *M0,double *lambda,double   *q,int *mode,double *res,int *info){
    zfeast_hrci_(ijob,N,Ze,work,workc,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,lambda,q,mode,res,info);
}


void zfeast_grcipevx(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  zfeast_grcipevx_(ijob,d,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_grcipev(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_grcipev_(ijob,d,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}
void zfeast_srcipevx(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info,double *Zne, double *Wne){
  zfeast_grcipevx_(ijob,d,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info,Zne,Wne);
}
void zfeast_srcipev(char *ijob,int *d,int *N,double *Ze,double *work,double *workc,double *Aq,double *Sq,int *fpm,double *epsout,int *loop,double *Emid,double *r,int *M0,double *lambda,double *q,int *mode,double *res,int *info){
  zfeast_grcipev_(ijob,d,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emid,r,M0,lambda,q,mode,res,info);
}




//////////////////////////////////////////////////////////////////////////


