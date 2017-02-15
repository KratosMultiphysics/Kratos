/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

#ifndef _MMGCOMMON_H
#define _MMGCOMMON_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <complex.h>

#define POSIX
#define GNU

#if (defined(__APPLE__) && defined(__MACH__))
#include <sys/sysctl.h>
#elif defined(__unix__) || defined(__unix) || defined(unix)
#include <unistd.h>
#elif defined(_WIN16) || defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) || defined(__TOS_WIN__) || defined(__WINDOWS__)
#define _WIN32_WINNT 0x0500
#include <windows.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

#include "eigenv.h"
#include "libmmgcommon.h"

#define MG_VER   "5.2.2"
#define MG_REL   "Feb 8, 2017"
#define MG_CPY   "Copyright (c) IMB-LJLL, 2004-"
#define MG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

/** Check if \a a and \a b have the same sign */
#define MG_SMSGN(a,b)  (((double)(a)*(double)(b) > (0.0)) ? (1) : (0))

/** size of box for renumbering with scotch. */
#define _MMG5_BOXSIZE 500

/** Maximal memory used if available memory compitation fail. */
#define _MMG5_MEMMAX  800

/* Domain refs in iso mode */
#define MG_PLUS    2
#define MG_MINUS   3

/* numerical accuracy */
#define _MMG5_ANGEDG    0.707106781186548   /*0.573576436351046 */
#define _MMG5_ANGLIM   -0.999999
#define _MMG5_ATHIRD    0.333333333333333

#define _MMG5_EPSD      1.e-30
#define _MMG5_EPSD2     1.0e-200
#define _MMG5_EPS       1.e-06
#define _MMG5_EPSOK     1.e-18
#define _MMG5_EPS3      1.e-42

#define _MMG5_SQR32     0.866025403784439

#ifndef M_PI
#define M_PI            3.14159265358979323846   /**< pi   */
#define M_PI_2          1.57079632679489661923   /**< pi/2 */
#endif

#define A64TH     0.015625
#define A16TH     0.0625
#define A32TH     0.03125

/* Macros */
#define MG_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MG_MIN(a,b) (((a) < (b)) ? (a) : (b))

/* tags */
#define  MG_NOTAG     (0)
#define  MG_REF       (1 << 0)        /**< 1  edge reference  */
#define  MG_GEO       (1 << 1)        /**< 2  geometric ridge */
#define  MG_REQ       (1 << 2)        /**< 4  required entity */
#define  MG_NOM       (1 << 3)        /**< 8  non manifold    */
#define  MG_BDY       (1 << 4)        /**< 16  boundary entity */
#define  MG_CRN       (1 << 5)        /**< 32  corner         */
#define  MG_NOSURF    (1 << 6)        /**< 64  freezed boundary */
// 1 << 13 reserveed for MG_INTERFACE tag of parMmg
#define  MG_NUL       (1 << 14)       /**< 16384 vertex removed */

/* binary tags for local parameters */
#define  MG_Vert   (1 << 0 )  /**< 1 local parameter applied over vertex */
#define  MG_Tria   (1 << 1 )  /**< 2 local parameter applied over triangle */
#define  MG_Tetra  (1 << 2 )  /**< 4 local parameter applied over tetrahedron */

#define MG_VOK(ppt)      (ppt && ((ppt)->tag < MG_NUL)) /**< Vertex OK */
#define MG_EOK(pt)       (pt && ((pt)->v[0] > 0))       /**< Element OK */

#define MG_EDG(tag) ((tag & MG_GEO) || (tag & MG_REF)) /**< Edge or Ridge */
#define MG_SIN(tag) ((tag & MG_CRN) || (tag & MG_REQ)) /**< Corner or Required */

#define MG_SET(flag,bit) ((flag) |= (1 << (bit)))  /**< bit number bit is set to 1 */
#define MG_CLR(flag,bit) ((flag) &= ~(1 << (bit))) /**< bit number bit is set to 0 */
#define MG_GET(flag,bit) ((flag) & (1 << (bit)))   /**< return bit number bit value */

#define _MMG5_KA 7 /*!< Key for hash tables. */
#define _MMG5_KB 11  /*!< Key for hash tables. */


/** Reset the customized signals and set the internal counters of points, edges,
 * tria and tetra to the suitable value (needed by users to recover their mesh
 * using the API) */
#define _LIBMMG5_RETURN(mesh,met,val)do          \
{                                                \
  signal(SIGABRT,SIG_DFL);                       \
  signal(SIGFPE,SIG_DFL);                        \
  signal(SIGILL,SIG_DFL);                        \
  signal(SIGSEGV,SIG_DFL);                       \
  signal(SIGTERM,SIG_DFL);                       \
  signal(SIGINT,SIG_DFL);                        \
  mesh->npi = mesh->np;                          \
  mesh->nti = mesh->nt;                          \
  mesh->nai = mesh->na;                          \
  mesh->nei = mesh->ne;                          \
  met->npi  = met->np;                           \
  return(val);                                   \
}while(0)

/* Macros for memory management */
/** Check if used memory overflow maximal authorized memory.
    Execute the command law if lack of memory. */
#define _MMG5_CHK_MEM(mesh,size,string,law) do                          \
  {                                                                     \
    if ( ((mesh)->memCur) > ((mesh)->memMax) ||                         \
         ((mesh)->memCur < 0 )) {                                       \
      fprintf(stderr,"  ## Error:");                                    \
      fprintf(stderr," unable to allocate %s.\n",string);               \
      fprintf(stderr,"  ## Check the mesh size or ");                   \
      fprintf(stderr,"increase maximal authorized memory with the -m option.\n"); \
      (mesh)->memCur -= (long long)(size);                              \
      law;                                                              \
    }                                                                   \
  }while(0)

/** Free pointer ptr of mesh structure and compute the new used memory.
    size is the size of the pointer */
#define _MMG5_DEL_MEM(mesh,ptr,size) do         \
  {                                             \
    (mesh)->memCur -= (long long)(size);        \
    free(ptr);                                  \
    ptr = NULL;                                 \
  }while(0)

/** Increment memory counter memCur and check if we don't overflow
    the maximum authorizied memory memMax. */
#define _MMG5_ADD_MEM(mesh,size,message,law) do \
  {                                             \
    (mesh)->memCur += (long long)(size);        \
    _MMG5_CHK_MEM(mesh,size,message,law);       \
  }while(0)

/** Safe deallocation */
#define _MMG5_SAFE_FREE(ptr) do                 \
  {                                             \
    free(ptr);                                  \
    ptr = NULL;                                 \
  }while(0)

/** Safe allocation with calloc */
#define _MMG5_SAFE_CALLOC(ptr,size,type) do     \
  {                                             \
    ptr = (type *)calloc((size),sizeof(type));  \
    if ( !ptr ) {                               \
      perror("  ## Memory problem: calloc");    \
      exit(EXIT_FAILURE);                       \
    }                                           \
  }while(0)

/** Safe allocation with malloc */
#define _MMG5_SAFE_MALLOC(ptr,size,type) do     \
  {                                             \
    ptr = (type *)malloc((size)*sizeof(type));  \
    if ( !ptr ) {                               \
      perror("  ## Memory problem: malloc");    \
      exit(EXIT_FAILURE);                       \
    }                                           \
  }while(0)

/** Safe reallocation */
#define _MMG5_SAFE_REALLOC(ptr,size,type,message) do        \
  {                                                         \
    type* tmp;                                              \
    tmp = (type *)realloc((ptr),(size)*sizeof(type));       \
    if ( !tmp ) {                                           \
      _MMG5_SAFE_FREE(ptr);                                 \
      perror(" ## Memory problem: realloc");                \
      exit(EXIT_FAILURE);                                   \
    }                                                       \
                                                            \
    if ( abs(mesh->info.imprim) > 6 || mesh->info.ddebug )  \
      fprintf(stdout,                                       \
              "  ## Warning: %s:%d: %s reallocation.\n",    \
              __FILE__,__LINE__,message);                   \
                                                            \
                                                            \
    (ptr) = tmp;                                            \
  }while(0)

/** safe reallocation with memset at 0 for the new values of tab */
#define _MMG5_SAFE_RECALLOC(ptr,prevSize,newSize,type,message) do \
  {                                                               \
    type* tmp;                                                    \
    int k;                                                        \
                                                                  \
    tmp = (type *)realloc((ptr),(newSize)*sizeof(type));          \
    if ( !tmp ) {                                                 \
      _MMG5_SAFE_FREE(ptr);                                       \
      perror(" ## Memory problem: realloc");                      \
      exit(EXIT_FAILURE);                                         \
    }                                                             \
                                                                  \
    if ( abs(mesh->info.imprim) > 6 || mesh->info.ddebug )        \
      fprintf(stdout,                                             \
              "  ## Warning: %s:%d: %s reallocation.\n",          \
              __FILE__,__LINE__,message);                         \
                                                                  \
    (ptr) = tmp;                                                  \
    for ( k=prevSize; k<newSize; k++) {                           \
      memset(&ptr[k],0,sizeof(type));                             \
    }                                                             \
  }while(0)

/** Reallocation of ptr of type type at size (initSize+wantedGap*initSize)
    if possible or at maximum available size if not. Execute the command law
    if reallocation failed. Memset to 0 for the new values of table. */
#define _MMG5_TAB_RECALLOC(mesh,ptr,initSize,wantedGap,type,message,law) do \
  {                                                                     \
    int gap;                                                            \
                                                                        \
    if ( (mesh->memMax-mesh->memCur) <                                  \
         (long long) (wantedGap*initSize*sizeof(type)) ) {              \
      gap = (int)((mesh->memMax-mesh->memCur)/sizeof(type));            \
      if(gap<1) {                                                       \
        fprintf(stderr,"  ## Error:");                                  \
        fprintf(stderr," unable to allocate %s.\n",message);            \
        fprintf(stderr,"  ## Check the mesh size or ");                 \
        fprintf(stderr,"increase maximal authorized memory with the -m option.\n"); \
        law;                                                            \
      }                                                                 \
    }                                                                   \
    else                                                                \
      gap = (int)(wantedGap*initSize);                                  \
                                                                        \
    _MMG5_ADD_MEM(mesh,gap*sizeof(type),message,law);                   \
    _MMG5_SAFE_RECALLOC((ptr),initSize+1,initSize+gap+1,type,message);  \
    initSize = initSize+gap;                                            \
  }while(0);

/** Error message when lack of memory */
#define _MMG5_INCREASE_MEM_MESSAGE() do                     \
  {                                                         \
    printf("  ## Check the mesh size or increase maximal"); \
    printf(" authorized memory with the -m option.\n");     \
  } while(0)


/** Inlined functions for libraries and executables */
#ifdef USE_SCOTCH
/** Warn user that we overflow asked memory during scotch call */
static inline
void _MMG5_warnScotch(MMG5_pMesh mesh) {
  if ( mesh->info.imprim > 4 || mesh->info.ddebug ) {
    if ( mesh->info.mem >= 0 ) {
      fprintf(stdout,"  ## Warning: we will overflow the memory asked with \"-m\"");
      fprintf(stdout," option during Scotch call.\n" );
    }
  }
}
#endif
/**
 * \param sigid signal number.
 *
 * Signal handling: specify error messages depending from catched signal.
 *
 */
static inline
void _MMG5_excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
  case SIGABRT:
    fprintf(stdout,"  *** potential lack of memory.\n");  break;
  case SIGFPE:
    fprintf(stdout,"  *** Floating-point exception\n"); break;
  case SIGILL:
    fprintf(stdout,"  *** Illegal instruction\n"); break;
  case SIGSEGV:
    fprintf(stdout,"  *** Segmentation fault\n");  break;
  case SIGTERM:
  case SIGINT:
    fprintf(stdout,"  *** Program killed\n");  break;
  }
  exit(EXIT_FAILURE);
}

/* Macro for fortran function generation */
/**
 * \def FORTRAN_NAME(nu,nl,pl,pc)
 * \brief Adds function definitions.
 * \param nu function name in upper case.
 * \param nl function name in lower case.
 * \param pl type of arguments.
 * \param pc name of arguments.
 * \note Macro coming from Scotch library.
 *
 * Adds function definitions with upcase, underscore and double
 * underscore to match any fortran compiler.
 *
 */
#define FORTRAN_NAME(nu,nl,pl,pc)               \
  void nu pl;                                   \
  void nl pl                                    \
  { nu pc; }                                    \
  void nl##_ pl                                 \
  { nu pc; }                                    \
  void nl##__ pl                                \
  { nu pc; }                                    \
  void nu pl

/**
 * \def FORTRAN_VARIADIC(nu,nl,pl,body)
 * \brief Adds function definitions.
 * \param nu function name in upper case.
 * \param nl function name in lower case.
 * \param pl type of arguments.
 * \param body body of the function.
 *
 * Adds function definitions with upcase, underscore and double
 * underscore to match any fortran compiler.
 *
 */
#define FORTRAN_VARIADIC(nu,nl,pl,body)           \
  void nu pl                                      \
  { body }                                        \
  void nl pl                                      \
  { body }                                        \
  void nl##_ pl                                   \
  { body }                                        \
  void nl##__ pl                                  \
  { body }                                        \


/* Global variables */
  static const unsigned char _MMG5_inxt2[6] = {1,2,0,1,2}; /*!< next vertex of triangle: {1,2,0} */
static const unsigned char _MMG5_iprv2[3] = {2,0,1}; /*!< previous vertex of triangle: {2,0,1} */

/* Private structures */
/**
 * \struct _MMG5_Bezier
 *
 * Store the Bezier definition of a surface triangle.
 *
 */
typedef struct {
  double       b[10][3];/*!< Bezier basis functions */
  double       n[6][3]; /*!< Normals at points */
  double       t[6][3]; /*!< Tangents at points */
  MMG5_pPoint  p[3];    /*!< Triangle vertices */
} _MMG5_Bezier;
typedef _MMG5_Bezier * _MMG5_pBezier;

/**
 * \struct _MMG5_hedge
 * \brief Used to hash edges (memory economy compared to \ref MMG5_hgeom).
 */
typedef struct {
  int   a,b,nxt;
  int   k; /*!< k = point along edge a b or triangle index */
  int   s;
} _MMG5_hedge;

/**
 * \struct _MMG5_Hash
 * \brief Identic as \ref MMG5_HGeom but use \ref _MMG5_hedge to store edges
 * instead of \ref MMG5_hgeom (memory economy).
 */
typedef struct {
  int     siz,max,nxt;
  _MMG5_hedge  *item;
} _MMG5_Hash;


/**
 * \struct _MMG5_iNode
 * \brief Cell for linked list of integer value.
 */
typedef struct _MMG5_iNode_s {
  int val;
  struct _MMG5_iNode_s *nxt;
} _MMG5_iNode;

/**
 * \struct _MMG5_dNode
 * \brief Cell for linked list of double value.
 */
typedef struct _MMG5_dNode_s {
  int    k;
  double val;
  struct _MMG5_dNode_s *nxt;
} _MMG5_dNode;


/* Functions declarations */
extern double _MMG5_det3pt1vec(double c0[3],double c1[3],double c2[3],double v[3]);
extern double _MMG5_det4pt(double c0[3],double c1[3],double c2[3],double c3[3]);
extern double _MMG5_orvol(MMG5_pPoint point,int *v);
extern int _MMG5_Add_inode( MMG5_pMesh mesh, _MMG5_iNode **liLi, int val );
extern int _MMG5_Alloc_inode( MMG5_pMesh mesh, _MMG5_iNode **node );
extern int _MMG5_Add_dnode( MMG5_pMesh mesh, _MMG5_dNode **liLi, int, double);
extern int _MMG5_Alloc_dnode( MMG5_pMesh mesh, _MMG5_dNode **node );
extern void   _MMG5_bezierEdge(MMG5_pMesh, int, int, double*, double*, char,double*);
int    _MMG5_buildridmet(MMG5_pMesh,MMG5_pSol,int,double,double,double,double*);
extern int    _MMG5_buildridmetfic(MMG5_pMesh,double*,double*,double,double,double,double*);
int    _MMG5_buildridmetnor(MMG5_pMesh, MMG5_pSol, int,double*, double*);
int    _MMG5_paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]);
extern int    _MMG5_rmtr(double r[3][3],double m[6], double mr[6]);
int    _MMG5_boundingBox(MMG5_pMesh mesh);
int    _MMG5_boulec(MMG5_pMesh, int*, int, int i,double *tt);
int    _MMG5_boulen(MMG5_pMesh, int*, int, int i,double *nn);
int    _MMG5_bouler(MMG5_pMesh, int*, int, int i,int *,int *,int *, int);
double _MMG5_caltri33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt);
extern double _MMG5_caltri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
extern double _MMG5_caltri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
void   _MMG5_defUninitSize(MMG5_pMesh mesh,MMG5_pSol met, char ismet);
void   _MMG5_displayHisto(MMG5_pMesh,int, double*, int, int, double, int, int,
                          double, double*, int*);
int    _MMG5_elementWeight(MMG5_pMesh,MMG5_pSol,MMG5_pTria,MMG5_pPoint,
                           _MMG5_Bezier*,double r[3][3],double gv[2]);
void   _MMG5_fillDefmetregSys( int, MMG5_pPoint, int, _MMG5_Bezier,double r[3][3],
                               double *, double *, double *, double *);
extern void _MMG5_Free_ilinkedList( MMG5_pMesh mesh, _MMG5_iNode *liLi );
extern void _MMG5_Free_dlinkedList( MMG5_pMesh mesh, _MMG5_dNode *liLi );
int    _MMG5_grad2metSurf(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pTria pt, int i);
int    _MMG5_hashEdge(MMG5_pMesh mesh,_MMG5_Hash *hash,int a,int b,int k);
int    _MMG5_hashGet(_MMG5_Hash *hash,int a,int b);
int    _MMG5_hashNew(MMG5_pMesh mesh, _MMG5_Hash *hash,int hsiz,int hmax);
int    _MMG5_intmetsavedir(MMG5_pMesh mesh, double *m,double *n,double *mr);
int    _MMG5_intridmet(MMG5_pMesh,MMG5_pSol,int,int,double,double*,double*);
int    _MMG5_mmgIntmet33_ani(double*,double*,double*,double);
int    _MMG5_mmgIntextmet(MMG5_pMesh,MMG5_pSol,int,double *,double *);
long long _MMG5_memSize(void);
void   _MMG5_mmgDefaultValues(MMG5_pMesh mesh);
int    _MMG5_mmgHashTria(MMG5_pMesh mesh, int *adja, _MMG5_Hash*, int chkISO);
void   _MMG5_mmgInit_parameters(MMG5_pMesh mesh);
void   _MMG5_mmgUsage(char *prog);
extern int    _MMG5_nonUnitNorPts(MMG5_pMesh,int,int,int,double*);
extern double _MMG5_nonorsurf(MMG5_pMesh mesh,MMG5_pTria pt);
extern int    _MMG5_norpts(MMG5_pMesh,int,int,int,double *);
extern int    _MMG5_nortri(MMG5_pMesh mesh,MMG5_pTria pt,double *n);
void   _MMG5_printTria(MMG5_pMesh mesh,char* fileName);
extern int    _MMG5_rotmatrix(double n[3],double r[3][3]);
int    _MMG5_invmat(double *m,double *mi);
int    _MMG5_invmatg(double m[9],double mi[9]);
double _MMG5_ridSizeInNormalDir(MMG5_pMesh,int,double*,_MMG5_pBezier,double,double);
double _MMG5_ridSizeInTangentDir(MMG5_pMesh, MMG5_pPoint,int,int*,double,double);
extern long _MMG5_safeLL2LCast(long long val);
int    _MMG5_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met);
int    _MMG5_scotchCall(MMG5_pMesh mesh, MMG5_pSol sol);
int    _MMG5_solveDefmetregSys( MMG5_pMesh, double r[3][3], double *, double *,
                                double *, double *, double, double, double);
int    _MMG5_solveDefmetrefSys( MMG5_pMesh,MMG5_pPoint,int*, double r[3][3],
                                double *, double *, double *, double *,
                                double, double, double);
double _MMG5_surftri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
double _MMG5_surftri33_ani(MMG5_pMesh,MMG5_pTria,double*,double*,double*);
double _MMG5_surftri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
extern int    _MMG5_sys33sym(double a[6], double b[3], double r[3]);
int    _MMG5_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol met);
int    _MMG5_interpreg_ani(MMG5_pMesh,MMG5_pSol,MMG5_pTria,char,double,double *mr);
int    _MMG5_interp_iso(double *ma,double *mb,double *mp,double t);
int    _MMG5_intersecmet22(MMG5_pMesh mesh, double *m,double *n,double *mr);
extern int _MMG5_countLocalParamAtTri( MMG5_pMesh,_MMG5_iNode **);
extern int _MMG5_writeLocalParamAtTri( MMG5_pMesh,_MMG5_iNode *,FILE*);
double MMG2_quickarea(double a[2],double b[2],double c[2]);

int MMG5_loadMshMesh_part1(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename,
                           FILE **inm,long *posNodes, long *posElts,
                           long *posNodeData, int *bin, int *iswp,
                           int *nelts);

int MMG5_loadMshMesh_part2(MMG5_pMesh mesh,MMG5_pSol sol,FILE **inm,
                           const long posNodes,const long posElts,
                           const long posNodeData,const int bin,const int iswp,
                           const int nelts);

/* function pointers */
int    (*_MMG5_chkmsh)(MMG5_pMesh,int,int);
int    (*_MMG5_bezierCP)(MMG5_pMesh ,MMG5_Tria *,_MMG5_pBezier ,char );
double (*_MMG5_lenSurfEdg)(MMG5_pMesh mesh,MMG5_pSol sol ,int ,int, char );

#ifdef USE_SCOTCH
int    (*_MMG5_renumbering)(int vertBoxNbr, MMG5_pMesh mesh, MMG5_pSol sol);
#endif

void   _MMG5_Set_commonFunc();

#ifdef __cplusplus
}
#endif

#endif
