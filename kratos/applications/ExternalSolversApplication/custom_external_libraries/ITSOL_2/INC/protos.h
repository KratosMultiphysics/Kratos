#ifndef __ITSOL_INCLUDED_PROTOS_H__
#define __ITSOL_INCLUDED_PROTOS_H__

#include <stdio.h>
#include "globheads.h"

#if defined(_SGI) || defined(_LINUX)
#define dnrm2   dnrm2_
#define ddot    ddot_
#define daxpy   daxpy_
#define qsplit  qsplit_
#define dscal   dscal_
#define dgemv   dgemv_
#define dgemm   dgemm_
#define dgetrf  dgetrf_
#define dgetri  dgetri_
#define dgesvd  dgesvd_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#define qsplit  qsplit_
#elif defined(_IBM)
#include <essl.h>
#define dnrm2   dnrm2
#define ddot    ddot
#define daxpy   daxpy
#define qsplit  qsplit
#define dscal   dscal
#define dgemv   dgemv
#define dgemm   dgemm
#define dgetrf  dgetrf
#define dgetri  dgetri
#define dgesvd  dgesvd
#define readmtc readmtc
#define csrcsc  csrcsc 
#define roscal  roscal
#define coscal  coscal
#define qsplit  qsplit
#else
#define dnrm2   dnrm2_
#define ddot    ddot_
#define daxpy   daxpy_
#define qsplit  qsplit_
#define dscal   dscal_
#define dgemv   dgemv_
#define dgemm   dgemm_
#define dgetrf  dgetrf_
#define dgetri  dgetri_
#define dgesvd  dgesvd_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#define qsplit  qsplit_
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

#define MAX_LINE        256
#define MAX_HBNAME      64
#define MAX_MAT			100

/* FORTRAN routines */
extern void readmtc(int*,  int*,  int*,  char*,  double*,  int*,
	     int*,  double*, int*,  char*,  int*,  int*,  int*,
	     char*,  char*, char*,  int*) ;
extern void csrcsc (int*, int*, int*, double*, int*, int*, double*,
		    int*, int*) ; 
extern void qsplit(double *a, int *ind, int *n, int *ncut);	
extern void dgesvd(char*, char*, int*, int*, double*, int*, double*,
		   double *, int*, double*, int*, double*, int*,
		   int*); 
extern void csrcoo ( int *, int *, int *, double *, int *, int *, int *,
		double *, int *, int *, int *);    

#if defined(_IBM)
#define DDOT(n,x,incx,y,incy)        ddot((n), (x), (incx), (y), (incy)) 
#define DCOPY(n,x,incx,y,incy)       dcopy((n), (x), (incx), (y), \
					   (incy)) 
#define DSCAL(n,alpha,x,incx)        dscal((n), (alpha), (x), (incx)) 
#define DAXPY(n,alpha,x,incx,y,incy) daxpy((n), (alpha), (x), (incx), \
					   (y), (incy)) 
#define DNRM2(n,x,incx)              dnrm2((n), (x), (incx))

#define IDMIN(n,sx,incx)             idmin((n), (sx), (incx))
#define DGEMV(transa,m,n,alpha,a,lda,x,incx,beta,y,incy)		\
  dgemv((transa), (m), (n),						\
	(alpha), (a), (lda), (x), (incx),				\
	(beta), (y), (incy))

#define DGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)		\
  dgemm((transa),(transb),						\
	(l),(n),(m),(alpha),(a),					\
	(lda),(b),(ldb),(beta),(c),(ldc))
#define DGETRF(m, n, a, lda, ipvt, info)  \
  dgetrf((m), (n), (a), (lda), (ipvt), (info))
#define DGETRI(n, a, lda, ipvt, work, lwork, info)		\
  dgetri((n), (a), (lda), (ipvt), (work), (lwork), (info))
#else
#define DDOT(n,x,incx,y,incy)        ddot(&(n),(x),&(incx),(y),&(incy))
#define DCOPY(n,x,incx,y,incy)       dcopy(&(n),(x),&(incx),(y),&(incy))
#define DSCAL(n,alpha,x,incx)        dscal(&(n),&(alpha),(x), &(incx))
#define DAXPY(n,alpha,x,incx,y,incy) daxpy(&(n), &(alpha), (x), \
					 &(incx), y, &(incy))
#define DNRM2(n, x, incx)            dnrm2(&(n), (x), &(incx))
#define IDMIN(n, sx, incx)           idmin((&(n), (sx), &(incx))
#define DGEMV(transa, m, n, alpha, a, lda, x, incx, beta, y, incy)  \
  dgemv((transa), &(m), &(n), &(alpha), (a), &(lda), (x), &(incx), \
	 &(beta), (y), &(incy))
#define DGEMM(transa,transb,l,n,m,alpha,a,lda,b,ldb,beta,c,ldc)	\
  dgemm((transa), (transb), &(l), &(n), &(m), &(alpha), (a),	\
	 &(lda), b, &(ldb), &(beta), (c), &(ldc)) 
#define DGETRF(m, n, a, lda, ipvt, info)		\
  dgetrf(&(m), &(n), (a), &(lda), (ipvt), (info))
#define DGETRI(n, a, lda, ipvt, work, lwork, info)			\
  dgetri(&(n), (a), &(lda), (ipvt), (work), &(lwork), (info))

extern double ddot(int *n, double *x, int *incx, double *y, int
		    *incy);  
extern void   dcopy(int *n, double *x, int *incx, double *y, int
		    *incy); 
extern void   dscal(int *n, double *alpha, double *x, int *incx);
extern void   daxpy(int *n, double *alpha, double *x, int *incx,
		    double *y, int *incy);
extern double dnrm2(int *n, double *x, int *incx);
extern void   idmin(int *n, double *sx, int *incx);
extern void   dgemv(char *transa, int *m, int *n, double *alpha,
		    double *a, int *lda, double *x, int *incx, double
		    *beta, double *y, int *incy);
extern void   dgemm(char *transa, char *transb, int *l, int *m, int
		    *n, double *alpha, double *a, int *lda, double
		    *b, int *ldb, double *beta, double *c, int *ldc);       
extern void   dgetrf(int *m, int *n, double *a, int *lda, int *ipvt,
		     int *info); 
extern void   dgetri(int *n, double *a, int *lda, int *ipvt, double
		     *work, int *lwork, int *info);
#endif 

/* sets.c */
extern int nnz_arms (arms PreSt,  FILE *ft);
extern void errexit(char *f_str, ...);
extern void *Malloc(int nbytes, char *msg); 
extern int setupCS(csptr amat, int len, int job); 
extern int cleanCS(csptr amat);
extern int nnz_cs (csptr A) ;
extern int cscpy(csptr amat, csptr bmat);
extern int setupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E);
extern int setupVBMat(vbsptr vbmat, int n, int *nB);
extern int setupILUT(ilutptr amat, int len);
extern int cleanVBMat(vbsptr vbmat); 
extern int nnzVBMat(vbsptr vbmat) ;
extern int memVBMat(vbsptr vbmat); 
extern int setupVBILU(vbiluptr lu, int n, int *bsz);
extern int cleanVBILU(vbiluptr lu); 
extern int cleanILU( iluptr lu );
extern int cleanILUT(ilutptr amat, int indic);
extern int cleanP4(p4ptr amat);
extern int mallocVBRow(vbiluptr lu, int nrow); 
extern int mallocRow( iluptr lu, int nrow );
extern void zrmC(int m, int n, BData data); 
extern void copyBData(int m, int n, BData dst, BData src, int isig); 
extern int CSRcs(int n, double *a, int *ja, int *ia, csptr mat, int 
		 rsa); 
extern int csrvbsrC(int job, int nBlk, int *nB, csptr csmat, vbsptr
		    vbmat);  
extern int col2vbcol( int col, vbsptr vbmat );
extern int nnz_vbilu(vbiluptr lu); 
extern int nnz_lev4(p4ptr levmat, int *lev, FILE *ft);
extern int setupILU( iluptr lu, int n );
extern int CS2lum( int n, csptr Amat, iluptr mat, int typ);
extern int COOcs(int n, int nnz,  double *a, int *ja, int *ia, csptr bmat);
void coocsr_(int*, int*, double*, int*, int*, double*, int*, int*);

/* MatOps.c */
extern int diag_scal(vbsptr vbmat);
extern int diagvec(vbsptr vbmat, BData x, BData y);
extern void matvec(csptr mata, double *x, double *y); 
extern void matvecCSR(SMatptr mat, double *x, double *y);
extern void matvecz(csptr mata, double *x, double *y, double *z);
extern void vbmatvec(vbsptr vbmat, double *x, double *y);
extern void luinv(int n, double *a, double *x, double *y); 
extern int vblusolC(double *y, double *x, vbiluptr lu); 
extern int lusolC( double *y, double *x, iluptr lu ); 
extern int rpermC(csptr mat, int *perm); 
extern int cpermC(csptr mat, int *perm) ; 
extern int dpermC(csptr mat, int *perm) ; 
extern int CSparTran(csptr amat, csptr bmat, CompressType *compress);
extern double vbnorm2(int sz, double *a);
extern void Lsol(csptr mata, double *b, double *x);
extern void Usol(csptr mata, double *b, double *x);
extern int ascend (p4ptr levmat, double *x, double *wk);
extern int descend(p4ptr levmat, double *x, double *wk);
extern int armsol2(double *x,  arms Prec);
extern int condestArms(arms armspre, double *y, FILE *fp );
extern int VBcondestC( vbiluptr, FILE *fp ); 
extern int CondestLUM(iluptr lu, double *y, double *x, FILE *fp );
extern void matvecVBR(SMatptr mat, double *x, double *y);
extern void matvecLDU(SMatptr mat, double *x, double *y);
extern int  preconILU(double *x, double *y, SPreptr mat);
extern int  preconVBR(double *x, double *y, SPreptr mat);
extern int  preconLDU(double *x, double *y, SPreptr mat);
extern int  preconARMS(double *x, double *y, SPreptr mat);
extern p4ptr Lvsol2(double *x, int nlev, p4ptr levmat, ilutptr ilusch) ;
extern int   Uvsol2(double *x, int nlev, int n, p4ptr levmat, ilutptr
		    ilusch); 
extern void  SchLsol(ilutptr ilusch, double *y) ;
extern void  SchUsol(ilutptr ilusch, double *y) ;
#ifndef _IBM
extern void dgemv(char*, int *, int*, double*, double *, int*,
		  double*, int*, double*, double*, int*);  
extern void dgemm(char*, char*, int*, int*, int*, double*, double*,
		  int*, double*, int*, double*, double*, int*) ;  
extern void dgetrf(int*, int*, double*, int*, int*, int*); 
extern void dgetri(int*, double*, int*, int*, double*,  int*, int* );
extern double dnrm2( int *, double *, int * );
extern void dscal(int*, double*, double*, int*); 
#endif 
extern int invGauss(int nn, double *A); 
extern int invSVD(int nn, double *A) ;

/* setblks.c */ 
extern int KeyComp(const void *vfst, const void *vsnd);
extern int init_blocks(csptr, int *, int **, int **, double, double *,
		       double *);  

/* upper directory */
extern int vbilukC( int lofM, vbsptr vbmat, vbiluptr lu, FILE *fp ); 
extern int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol,
		    BData *w, FILE *fp ); 
extern int ilutc(iluptr mt, iluptr lu, int lfil, double tol, int drop,
		 FILE *fp ); 
extern int ilukC( int lofM, csptr csmat, iluptr lu, FILE *fp );
extern int ilut( csptr csmat, iluptr lu, int lfil, double tol,
		 FILE *fp );
extern   int fgmr(SMatptr Amat, SPreptr lu, double *rhs, double *sol, double tol,
	   int im, int *itmax, FILE *fits ); 
extern int arms2(csptr Amat, int *ipar, double *droptol, int *lfil, 
	  double tolind, arms PreMat, FILE *ft) ;
extern int condestLU( iluptr, FILE *);
extern int nnz_ilu( iluptr lu ); 
extern void roscal (int* nrow, int* job, int* nrm, double *a, int *ja,
		    int *ia, double *diag, double *b, int *jb, int
		    *ib, int *ierr) ;  
extern void coscal (int* nrow, int* job, int* nrm, double *a, int *ja,
		    int *ia, double *diag, double *b, int *jb, int
		    *ib, int *ierr) ;  
extern int outputLU( iluptr lu, char *filename );
extern int lumsolC(double *y, double *x, iluptr lu );
extern void lumatvec(iluptr mat, double *x, double *y);
extern int CSClum( int n, double *a, int *ja, int *ia, iluptr mat, int rsa ); 
extern int CSClumC(csptr amat, iluptr mat, int rsa ); 
extern void setup_arms (arms Levmat);
extern int cleanARMS(arms ArmsPre);
extern int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F,
		    csptr E, csptr C); 
/* misc.c */
extern int SparTran(csptr amat, csptr bmat, int job, int flag); 
extern int coscalC(csptr mata, double *diag, int nrm);
extern void dscale(int n, double *dd, double *x, double * y);
extern void hilosort(csptr mat, int abval, int hilo);
extern void printmat(FILE *ft, csptr A, int i0, int i1);
extern void qqsort(int *ja, double *ma, int left, int right);
extern void qsort2C(int *ja, double *ma, int left, int right, int
		    abval); 
extern void qsort3i(int *wa, int *cor1, int *cor2, int left, int
		    right); 
extern void qsortC(int *ja, double *ma, int left, int right, int
		   abval); 
extern void qsortR2I(double *wa, int *cor1, int *cor2, int left, int
		     right); 
extern int qsplitC(double *a, int *ind, int n, int ncut);
extern int roscalC(csptr mata, double *diag, int nrm);
extern void swapj(int v[], int i, int j);
extern void swapm(double v[], int i, int j);

/* piluNEW.c */
extern int pilu(p4ptr amat, csptr B, csptr C, double *droptol, int
		*lfil, csptr schur);

/* ilutpC.c */
extern int ilutD(csptr amat, double *droptol, int *lfil, ilutptr
		 ilusch);
extern int ilutpC(csptr amat, double *droptol, int *lfil, double
		  permtol, int mband, ilutptr ilusch);

/* PQ.c */
extern int PQperm(csptr mat, int bsize, int *Pord, int *Qord, int
		  *nnod, double tol);
extern int add2com(int *nback, int nod, int *iord, int *riord);
extern int add2is(int *last, int nod, int *iord, int *riord);
extern int indsetC(csptr mat, int bsize, int *iord, int *nnod, double
		   tol); 
extern int preSel(csptr mat, int *icor, int *jcor, int job, double
		  tol, int *count);
/* indsetC.c */
extern int weightsC(csptr mat, double *w);

/* setblks.c */
extern int KeyComp( const void *vfst, const void *vsnd );
extern int init_blocks( csptr csmat, int *pnBlock, int **pnB, int
			**pperm, double eps, double *t_hash, double
			*t_angle );

/* systimer.c */
extern double sys_timer(void);

/*auxill.c */
extern void randvec (double *v, int n);		   	       

#endif 
