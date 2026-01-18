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

#ifndef MMGCOMMON_H
#define MMGCOMMON_H

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
#include "mmg/common/mmgcmakedefines.h"

#if (defined(__APPLE__) && defined(__MACH__))
#include <sys/sysctl.h>
#elif defined(__unix__) || defined(__unix) || defined(unix)
#include <unistd.h>
#elif defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) || defined(__TOS_WIN__) || defined(__WINDOWS__)
#include <windows.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "eigenv_private.h"
#include "libmmgcommon_private.h"

#define MG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

/**
 * Maximum array size when storing list of tria containing a vertex.
 */
#define MMG5_TRIA_LMAX  1024

/**
 * Maximum array size when calling boulep.
 */
#define MMG5_LMAX      10240

/** Maximal number of local parameters */
#define MMG5_LPARMAX   200

/** Check if \a a and \a b have the same sign */
#define MG_SMSGN(a,b)  (((double)(a)*(double)(b) > (0.0)) ? (1) : (0))

/** size of box for renumbering with scotch. */
#define MMG5_BOXSIZE 500

/** Maximal memory used if available memory compitation fail. */
#define MMG5_MEMMAX  800        /**< Default mem if unable to compute memMax */
#define MMG5_BITWIZE_MB_TO_B 20 /**< Bitwise convertion from Mo to O */
#define MMG5_MEMPERCENT 0.5     /**< Percent of RAM used by default */

/* Macro for unset or unititialized mark */
#define MMG5_UNSET -1

/* reference of the boundary that moves in lagrangian mode */
#define MMG5_DISPREF 10

/* million */
#define MMG5_MILLION 1048576

/* numerical accuracy */
#define MMG5_ANGEDG    0.707106781186548   /*0.573576436351046 */
#define MMG5_ANGLIM   -0.999999
#define MMG5_ATHIRD    0.333333333333333

#define MMG5_EPSD      1.e-30
#define MMG5_EPSD2     1.0e-200
#define MMG5_EPS       1.e-06
#define MMG5_EPSOK     1.e-15
#define MMG5_NULKAL    1.e-30

#define MMG5_SQR32     0.866025403784439

#ifndef M_PI
#define M_PI            3.14159265358979323846   /**< pi   */
#define M_PI_2          1.57079632679489661923   /**< pi/2 */
#endif

#define A64TH     0.015625
#define A16TH     0.0625
#define A32TH     0.03125

#define MMG5_MEMMIN 38 /**< minimal memory needed to store the mesh/sol names */

#define MMG5_PATHSEP  '/'

/* Variables for option initialization */
#define MMG5_NONSET_MEM  -1 /**< mem value for unspecified max memory */
#define MMG5_NONSET_HMIN -1 /**< hmin value for unspecified hmin size */
#define MMG5_NONSET_HMAX -1 /**< hmax value for unspecified hmax size */
#define MMG5_NONSET_HSIZ -1 /**< hsiz value for unspecified hsiz map */
#define MMG5_NONSET      -1
#define MMG5_HAUSD     0.01 /**< default value for hausdorff param */
#define MMG5_HGRAD     0.26236426446 /**< default value for gradation (1.3) */
#define MMG5_HGRADREQ  0.83290912294 /**< default value for required gradation (2.3) */
#define MMG5_NOHGRAD     -1 /**< disable gradation */
#define MMG5_LAG         -1 /**< default value for lagrangian option */
#define MMG5_NR          -1 /**< no ridge detection */
#define MMG5_LS         0.0 /**< default level-set to discretize */
#define MMG5_XREG       0.4 /**< default relaxation parameter for coordinate regularization */
#define MMG5_PROCTREE    32 /**< default size of the PROctree */
#define MMG5_OFF          0 /**< 0 */
#define MMG5_ON           1 /**< 1 */
#define MMG5_GAP        0.2 /**< gap value for reallocation */
#define MMG5_HMINCOE  0.001 /**< coefficient to compute the default hmin value */
#define MMG5_HMAXCOE      2 /**< coefficient to compute the default hmax value */
#define MMG5_HMINMAXGAP   5 /**< imposed gap between hmin and hmax if hmax<hmin */
#define MMG5_FEM          1 /**< defaut value for FEM mode */
#define MMG5_FILESTR_LGTH 128 /** Maximal length of a line in input file */

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
#define  MG_OPNBDY    (1 << 7)        /**< 128 open boundary */
#define  MG_OLDPARBDY (1 << 11)       /**< 2048 old parallel boundary */
#define  MG_PARBDYBDY (1 << 12)       /**< 4096 parallel boundary over a boundary */
#define  MG_PARBDY    (1 << 13)       /**< 8192 parallel boundary */
#define  MG_OVERLAP   (1 << 14)       /**< 16384 elements on overlap */
#define  MG_NUL       (1 << 15)       /**< 32768 vertex removed */

/* binary tags for local parameters */
#define  MG_Vert   (1 << 0 )  /**< 1 local parameter applied over vertex */
#define  MG_Tria   (1 << 1 )  /**< 2 local parameter applied over triangle */
#define  MG_Tetra  (1 << 2 )  /**< 4 local parameter applied over tetrahedron */
#define  MG_Edge   (1 << 3 )  /**< 8 local parameter applied over edge */


#define MG_VOK(ppt)      (ppt && ((ppt)->tag < MG_NUL)) /**< Vertex OK */
#define MG_EOK(pt)       (pt && ((pt)->v[0] > 0))       /**< Element OK */

#define MG_SIN(tag)        ((tag & MG_CRN) || (tag & MG_REQ)) /**< Corner or Required */
#define MG_SIN_OR_NOM(tag) ( MG_SIN(tag) || (tag & MG_NOM) ) /**< Corner, Required or non-manifold */
#define MG_RID(tag)        ( ( !( MG_SIN_OR_NOM(tag)) )  && ( tag & MG_GEO ) ) /**< Non-singular ridge point (so ridge metric in aniso mode) */

#define MG_EDG(tag) ((tag & MG_GEO) || (tag & MG_REF)) /**< Edge or Ridge */
#define MG_GEO_OR_NOM(tag) (( tag & MG_GEO ) || ( tag & MG_NOM )) /**< Ridge or non-manifold */
#define MG_EDG_OR_NOM(tag) ( MG_EDG(tag) || (tag & MG_NOM ) ) /**< Edge, ridge or non-manifold */
#define MG_TRUE_BDY(tag) ( (tag & MG_BDY) && !(tag & MG_PARBDY) ) /**< Non parbdy boundary point (true bdy) */



#define MG_SET(flag,bit) ((flag) |= (1 << (bit)))  /**< bit number bit is set to 1 */
#define MG_CLR(flag,bit) ((flag) &= ~(1 << (bit))) /**< bit number bit is set to 0 */
#define MG_GET(flag,bit) ((flag) & (1 << (bit)))   /**< return bit number bit value */

#define MMG5_KA 7 /*!< Key for hash tables. */
#define MMG5_KB 11  /*!< Key for hash tables. */

/* file reading */
#define MMG5_SW 4
#define MMG5_SD 8

/** Reset the customized signals and set the internal counters of points, edges,
 * tria and tetra to the suitable value (needed by users to recover their mesh
 * using the API) */
#define _LIBMMG5_RETURN(mesh,sol,met,val)do      \
  {                                              \
    signal(SIGABRT,SIG_DFL);                     \
    signal(SIGFPE,SIG_DFL);                      \
    signal(SIGILL,SIG_DFL);                      \
    signal(SIGSEGV,SIG_DFL);                     \
    signal(SIGTERM,SIG_DFL);                     \
    signal(SIGINT,SIG_DFL);                      \
    mesh->npi = mesh->np;                        \
    mesh->nti = mesh->nt;                        \
    mesh->nai = mesh->na;                        \
    mesh->nei = mesh->ne;                        \
    mesh->xt  = 0;                               \
    if ( sol ) { sol->npi  = sol->np; }          \
    if ( met ) { met->npi  = met->np; }          \
    return val;                                  \
  }while(0)

/* Macros for memory management */
/** Check if used memory overflow maximal authorized memory.
    Execute the command law if lack of memory. */
#define MMG5_CHK_MEM(mesh,size,string,law) do                          \
  {                                                                     \
    if ( (mesh)->memCur > (mesh)->memMax ) {                            \
      fprintf(stderr,"  ## Error:");                                    \
      fprintf(stderr," unable to allocate %s.\n",string);               \
      fprintf(stderr,"  ## Check the mesh size or ");                   \
      fprintf(stderr,"increase maximal authorized memory with the -m option.\n"); \
      (mesh)->memCur -= (size);                                         \
      law;                                                              \
    }                                                                   \
  }while(0)

static inline
void * mycalloc(size_t c, size_t s) {
  char *ptr;
  ptr = (char *)calloc(c*s+sizeof(size_t),1);
  if (ptr == NULL)
    return NULL;
  else {
    *((size_t*)ptr)=c*s;
    ptr+=sizeof(size_t);
    return (void*)ptr;
  }
}

static inline
void * mymalloc(size_t s) {
  char *ptr;
  ptr = (char *)malloc(s+sizeof(size_t));
  if (ptr == NULL)
    return NULL;
  else {
    *((size_t*)ptr)=s;
    ptr+=sizeof(size_t);
    return (void*)ptr;
  }
}

static inline
void * myrealloc(void * ptr_in, size_t s, size_t oldsize) {
  char *ptr;
  char *ptr_in_c = (char*)ptr_in;

  if ( !ptr_in ) {
    assert ( !oldsize );
    return mymalloc( s );
  }

  ptr_in_c -= sizeof(size_t);
  if (oldsize != *((size_t*)ptr_in_c)) {
    fprintf(stderr, "myrealloc: Error: freed memory mismatch\n");
    assert(0);
  }
  ptr = (char *)realloc(ptr_in_c, s+sizeof(size_t));
  if (ptr == NULL)
    return NULL;
  else {
    *((size_t*)ptr)=s;
    ptr+=sizeof(size_t);
    return (void*)ptr;
  }
}

static inline
size_t myfree(void *ptr) {
  size_t s;
  char * ptr_c = (char*)ptr;

  if ( !ptr ) return 0;

  ptr_c = ptr_c-sizeof(size_t);
  s = *((size_t*)ptr_c);
  free(ptr_c);

  return s;
}

/** Free pointer ptr of mesh structure and compute the new used memory. */
#define MMG5_DEL_MEM(mesh,ptr) do              \
  {                                             \
    size_t size_to_free = myfree(ptr);          \
    (mesh)->memCur -= size_to_free;             \
    ptr = NULL;                                 \
  }while(0)

/** Increment memory counter memCur and check if we don't overflow
    the maximum authorizied memory memMax. */
#define MMG5_ADD_MEM(mesh,size,message,law) do \
  {                                             \
    (mesh)->memCur += (size);                   \
    MMG5_CHK_MEM(mesh,size,message,law);       \
  }while(0)

/** Safe deallocation */
#define MMG5_SAFE_FREE(ptr) do                 \
  {                                             \
    myfree(ptr);                                \
    ptr = NULL;                                 \
  }while(0)

/** Safe allocation with calloc */
#define MMG5_SAFE_CALLOC(ptr,size,type,law) do    \
  {                                                   \
    ptr = (type*)mycalloc(size,sizeof(type));         \
    if ( !ptr ) {                                     \
      perror("  ## Memory problem: calloc");          \
      law;                                            \
    }                                                 \
  }while(0)

/** Safe allocation with malloc */
#define MMG5_SAFE_MALLOC(ptr,size,type,law) do    \
  {                                                   \
    size_t size_to_allocate = (size)*sizeof(type);    \
    ptr = (type*)mymalloc(size_to_allocate);          \
    if ( !ptr ) {                                   \
      perror("  ## Memory problem: malloc");        \
      law;                                          \
    }                                               \
  }while(0)

/** Safe reallocation */
#define MMG5_SAFE_REALLOC(ptr,prevSize,newSize,type,message,law) do    \
  {                                                                     \
    type* tmp;                                                          \
    size_t size_to_allocate = (newSize)*sizeof(type);                   \
                                                                        \
    tmp = (type *)myrealloc((ptr),size_to_allocate,(prevSize)*sizeof(type)); \
    if ( !tmp ) {                                                       \
      MMG5_SAFE_FREE(ptr);                                             \
      perror(" ## Memory problem: realloc");                            \
      law;                                                              \
    }                                                                   \
                                                                        \
    (ptr) = tmp;                                                        \
  }while(0)

/** Check for int32 overflow when trying to reallocate coef*oldSiz+shift array */
#define MMG5_CHK_INT32_OVERFLOW(wantedGap,oldSiz,coef,shift,law) do     \
  {                                                                     \
    /* Check for int32 overflow */                                      \
    if ( sizeof(MMG5_int) == sizeof(int32_t) ) {                        \
      MMG5_int gap_loc = (MMG5_int)((wantedGap) * (oldSiz));            \
      if ( !gap_loc ) gap_loc     = 1;                                  \
                                                                        \
      int32_t max_ne = (INT32_MAX-(shift))/(coef);                      \
      if ( max_ne < (oldSiz)+gap_loc ) {                                \
        /* Detected overflow, target maximal possible size */           \
        gap_loc = max_ne-(oldSiz);                                      \
        if ( gap_loc <=0 ) {                                            \
          /* No possibe realloc without int overflow */                 \
          fprintf(stderr,"  ## Error: %s: %d: Unable to reallocate adja array" \
                  " without int overflow.\n",__func__,__LINE__);        \
          gap_loc = 0;                                                  \
          law;                                                          \
        }                                                               \
        else {                                                          \
          wantedGap = (float)gap_loc/(float)oldSiz;                     \
          printf("wantGap has been modified %15f\n",wantedGap);         \
          wantedGap = (double)gap_loc/(double)oldSiz;                   \
          printf("DwantGap has been modified %15fl\n",wantedGap);       \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }while(0)


/** safe reallocation with memset at 0 for the new values of tab */
#define MMG5_SAFE_RECALLOC(ptr,prevSize,newSize,type,message,law) do \
  {                                                                     \
    type* tmp;                                                          \
    size_t size_to_allocate = (newSize)*sizeof(type);                   \
                                                                        \
    tmp = (type *)myrealloc((ptr),size_to_allocate,(prevSize)*sizeof(type)); \
    if ( !tmp ) {                                                       \
      MMG5_SAFE_FREE(ptr);                                              \
      perror(" ## Memory problem: realloc");                            \
      law;                                                              \
    }                                                                   \
    else {                                                              \
      (ptr) = tmp;                                                      \
      assert(ptr);                                                      \
      if ( newSize > prevSize ) {                                       \
        memset(&((ptr)[prevSize]),0,((newSize)-(prevSize))*sizeof(type)); \
      }                                                                 \
    }                                                                   \
  }while(0)

/** Reallocation of ptr of type type at size (initSize+wantedGap*initSize)
    if possible or at maximum available size if not. Execute the command law
    if reallocation failed. Memset to 0 for the new values of table. */
#define MMG5_TAB_RECALLOC(mesh,ptr,initSize,wantedGap,type,message,law) do \
  {                                                                     \
    MMG5_int    gap;                                                    \
                                                                        \
    assert ( mesh->memCur < mesh->memMax );                             \
                                                                        \
    gap = (MMG5_int)(floor(wantedGap * initSize));                      \
    if ( !gap ) gap     = 1;                                            \
                                                                        \
    if ( mesh->memMax < mesh->memCur + gap*sizeof(type) ) {             \
      gap = (MMG5_int)((mesh->memMax-mesh->memCur)/sizeof(type));       \
      if(gap<1) {                                                       \
        fprintf(stderr,"  ## Error:");                                  \
        fprintf(stderr," unable to allocate %s.\n",message);            \
        fprintf(stderr,"  ## Check the mesh size or ");                 \
        fprintf(stderr,"increase maximal authorized memory with the -m option.\n"); \
        law;                                                            \
      }                                                                 \
    }                                                                   \
                                                                        \
    MMG5_ADD_MEM(mesh,gap*sizeof(type),message,law);                   \
    MMG5_SAFE_RECALLOC((ptr),initSize+1,initSize+gap+1,type,message,law); \
    initSize = initSize+gap;                                            \
  }while(0);

/** Error message when lack of memory */
#define MMG5_INCREASE_MEM_MESSAGE() do                     \
  {                                                         \
    printf("  ## Check the mesh size or increase maximal"); \
    printf(" authorized memory with the -m option.\n");     \
  } while(0)

#define MMG5_SAFELL2LCAST(longlongval) (((longlongval) > (LONG_MAX)) ? 0 : ((long)(longlongval)))
#define MMG5_SAFELL2ICAST(longlongval) (((longlongval) > (INT_MAX)) ? 0 : ((int)(longlongval)))

/** check the return value of fread */
#define MMG_FREAD(ptr,size,count,stream) do             \
  {                                                     \
                                                        \
    if ( count != fread(ptr,size,count,stream) ) {      \
      fputs ( "Reading error", stderr );                \
      return -1;                                        \
    }                                                   \
  } while(0);

/** macro to help to count the number of variadic arguments */
#define CV_VA_NUM_ARGS_HELPER(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...)    N

/** count the number of variadic arguments provided to the macro */
#define CV_VA_NUM_ARGS(...)      CV_VA_NUM_ARGS_HELPER(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

/**
 * check the return value of fscanf
 *
 * \remark don't work without any variadic arg
 * \remark don't work on MSVC because variadic args are not expanded
 */
#ifdef MMG_POSIX
#define MMG_FSCANF(stream,format,...) do                                \
  {                                                                     \
    int io_count   = fscanf(stream,format,__VA_ARGS__);                 \
    int args_count = CV_VA_NUM_ARGS(__VA_ARGS__);                       \
    if ( args_count != io_count ) {                                     \
      fprintf (stderr, "Reading error: fscanf counts %d args while %d provided\n",io_count,args_count ); \
      return -1;                                                        \
    }                                                                   \
  } while(0);
#else
#define MMG_FSCANF(stream,format,...) do                                \
  {                                                                     \
    int io_count   = fscanf(stream,format,__VA_ARGS__);                 \
    int args_count = CV_VA_NUM_ARGS(__VA_ARGS__);                       \
    if ( 0 > io_count ) {                                               \
      fprintf (stderr, "Reading error: fscanf counts %d args\n",io_count); \
      return -1;                                                        \
    }                                                                   \
  } while(0);
#endif

/** Inlined functions for libraries and executables */
#ifdef USE_SCOTCH
/** Warn user that we overflow asked memory during scotch call */
static inline
void MMG5_warnScotch(MMG5_pMesh mesh) {
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
void MMG5_excfun(int sigid) {
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

/**
 * \param fproto function prototype
 *
 * Expand automatically prototype of function pointer in .h/.c files depending
 * on the definition of the MMG_EXTERN and MMG_ASSIGN_NULL preprocessor
 * variables:
 *   - MMG_EXTERN is setted to "extern" in the .h file and empty in the .c one;
 *   - MMG_ASSIGN_NULL is empty in .h file and setted to =NULL in .c one.
 */
#define FUNCTION_POINTER(fproto)\
  MMG_EXTERN fproto MMG_ASSIGN_NULL



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
  static const uint8_t MMG5_inxt2[6] = {1,2,0,1,2}; /*!< next vertex of triangle: {1,2,0} */
  static const uint8_t MMG5_iprv2[3] = {2,0,1}; /*!< previous vertex of triangle: {2,0,1} */

/* Private structures */
/**
 * \struct MMG5_Bezier
 *
 * Store the Bezier definition of a surface triangle.
 *
 * \remark Numbering convention for high order points (b)
 * \verbatim
 *
 *     2                                                       *
 *     |`\                                                     *
 *     |  `\                                                   *
 *     5    `4                                                 *
 *     |      `\                                               *
 *     |        `\                                             *
 *     6          `3                                           *
 *     |            `\                                         *
 *     |              `\                                       *
 *     0 --- 7 --- 8 --- 1                                     *
 *
 * \endverbatim
 *
 */
typedef struct {
  double       b[10][3];/*!< Bezier basis functions */
  double       n[6][3]; /*!< Normals at points */
  double       t[6][3]; /*!< Tangents at points */
  MMG5_pPoint  p[3];    /*!< Triangle vertices */
} MMG5_Bezier;
typedef MMG5_Bezier * MMG5_pBezier;

/**
 * \struct MMG5_iNode
 * \brief Cell for linked list of integer value.
 */
typedef struct MMG5_iNode_s {
  MMG5_int val;
  struct MMG5_iNode_s *nxt;
} MMG5_iNode;

/* Functions declarations */
 void          MMG5_version(MMG5_pMesh,char*);
 extern void MMG5_nsort(int8_t ,double *,int8_t *);
 extern void MMG5_nperm(int8_t n,int8_t shift,int8_t stride,double *val,double *oldval,int8_t *perm);
 extern double MMG5_det3pt1vec(double c0[3],double c1[3],double c2[3],double v[3]);
 extern double MMG5_det4pt(double c0[3],double c1[3],double c2[3],double c3[3]);
 int           MMG5_devangle(double* n1, double *n2, double crit);
 extern double MMG5_orvol(MMG5_pPoint point,MMG5_int *v);
 int           MMG5_Add_inode( MMG5_pMesh mesh, MMG5_iNode **liLi, int val );
 int           MMG5_eigenvmatsym2d(MMG5_pMesh mesh,double m[],double lambda[],double v[][2]);
 int           MMG5_eigenvmatsym3d(MMG5_pMesh mesh,double m[],double lambda[],double v[][3]);
 int           MMG5_eigenvmatnonsym2d(MMG5_pMesh mesh,double m[],double lambda[],double v[][2]);
 int           MMG5_eigenvmatnonsym3d(MMG5_pMesh mesh,double m[],double lambda[],double v[][3]);
 extern void   MMG5_bezierEdge(MMG5_pMesh, MMG5_int, MMG5_int, double*, double*, int8_t,double*);
 int           MMG5_buildridmet(MMG5_pMesh,MMG5_pSol,MMG5_int,double,double,double,double*,double[3][3]);
 extern int    MMG5_buildridmetfic(MMG5_pMesh,double*,double*,double,double,double,double*);
 int           MMG5_buildridmetnor(MMG5_pMesh, MMG5_pSol, MMG5_int,double*, double*,double[3][3]);
 void          MMG5_check_hminhmax(MMG5_pMesh mesh, int8_t sethmin, int8_t sethmax);
 int           MMG5_paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]);
 void          MMG5_transpose3d(double m[3][3]);
 void          MMG5_dotprod(int8_t dim,double *a,double *b,double *result);
 void          MMG5_crossprod3d(double *a,double *b,double *result);
 void          MMG5_mn(double m[6], double n[6], double mn[9] );
 extern int    MMG5_rmtr(double r[3][3],double m[6], double mr[6]);
 int           MMG5_boundingBox(MMG5_pMesh mesh);
 int           MMG5_boulep(MMG5_pMesh mesh,MMG5_int start,int ip,MMG5_int*,MMG5_int *list, MMG5_int *tlist);
 int           MMG5_boulec(MMG5_pMesh, MMG5_int*, MMG5_int,int ip,double *tt);
 int           MMG5_boulen(MMG5_pMesh, MMG5_int*, MMG5_int,int ip,double *nn);
 int           MMG5_bouler(MMG5_pMesh, MMG5_int*, MMG5_int,int ip,MMG5_int *,MMG5_int *,int *, int*, int);
 int           MMG5_boulet(MMG5_pMesh mesh,MMG5_int start,int ip,MMG5_int *list,int8_t s,int8_t *opn);
 double        MMG5_caltri33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt);
 extern double MMG5_caltri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
 extern double MMG5_caltri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
 void          MMG5_defUninitSize(MMG5_pMesh mesh,MMG5_pSol met,int8_t ismet);
 void          MMG5_displayLengthHisto(MMG5_pMesh,MMG5_int,double*,MMG5_int,MMG5_int,double,
                                        MMG5_int,MMG5_int,double,int,double*,MMG5_int*,int8_t);
 void          MMG5_displayLengthHisto_internal( MMG5_int,MMG5_int,MMG5_int,double,
                                                 MMG5_int,MMG5_int,double, MMG5_int,double*,
                                                 MMG5_int*,int8_t,int);
 short         MMG5_dikmov(MMG5_pMesh,MMG5_pSol,short*,short,
                           MMG5_int chkmovmesh(MMG5_pMesh,MMG5_pSol,short,MMG5_int*));
 int           MMG5_minQualCheck ( MMG5_int iel, double minqual, double alpha );
 int           MMG5_elementWeight(MMG5_pMesh,MMG5_pSol,MMG5_pTria,MMG5_pPoint,
                                   MMG5_Bezier*,double r[3][3],double gv[2]);
 void          MMG5_fillDefmetregSys( MMG5_int, MMG5_pPoint, int, MMG5_Bezier,double r[3][3],
                                       double *, double *, double *, double *);
 void          MMG5_Free_ilinkedList( MMG5_pMesh mesh, MMG5_iNode *liLi );
 MMG5_int      MMG5_grad2metSurf(MMG5_pMesh,MMG5_pSol,MMG5_pTria,MMG5_int,MMG5_int);
 int           MMG5_grad2metSurfreq(MMG5_pMesh,MMG5_pSol,MMG5_pTria,MMG5_int,MMG5_int);
 MMG5_int      MMG5_hashFace(MMG5_pMesh,MMG5_Hash*,MMG5_int,MMG5_int,MMG5_int,MMG5_int);
 int           MMG5_hashEdge(MMG5_pMesh mesh,MMG5_Hash *hash,MMG5_int a,MMG5_int b,MMG5_int k);
 int           MMG5_hashUpdate(MMG5_Hash *hash,MMG5_int a,MMG5_int b,MMG5_int k);
 int           MMG5_hashEdgeTag(MMG5_pMesh mesh,MMG5_Hash *hash,MMG5_int a,MMG5_int b,uint16_t k);
 MMG5_int      MMG5_hashGet(MMG5_Hash *hash,MMG5_int a,MMG5_int b);
 int           MMG5_hashNew(MMG5_pMesh mesh, MMG5_Hash *hash,MMG5_int hsiz,MMG5_int hmax);
 int           MMG5_intmetsavedir(MMG5_pMesh mesh, double *m,double *n,double *mr);
 int           MMG5_intridmet(MMG5_pMesh,MMG5_pSol,MMG5_int,MMG5_int,double,double*,double*);
 int           MMG5_mmgIntmet33_ani(double*,double*,double*,double);
 int           MMG5_ismaniball(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int start, int8_t istart);
 int           MMG5_mmgIntextmet(MMG5_pMesh,MMG5_pSol,MMG5_int,double *,double *);
 size_t        MMG5_memSize(void);
 void          MMG5_memOption_memSet(MMG5_pMesh mesh);
 void          MMG5_mmgDefaultValues(MMG5_pMesh mesh);
 int           MMG5_mmgHashTria(MMG5_pMesh mesh, MMG5_int *adja, MMG5_Hash*, int chkISO);
 void          MMG5_mmgInit_parameters(MMG5_pMesh mesh);
 void          MMG5_mmgUsage(char *prog);
 void          MMG5_paramUsage1(void);
 void          MMG5_paramUsage2(void);
 void          MMG5_2d3dUsage(void);
 void          MMG5_lagUsage(void);
 void          MMG5_advancedUsage(void);
 extern int    MMG5_nonUnitNorPts(MMG5_pMesh,MMG5_int,MMG5_int,MMG5_int,double*);
 extern double MMG5_nonorsurf(MMG5_pMesh mesh,MMG5_pTria pt);
 extern int    MMG5_norpts(MMG5_pMesh,MMG5_int,MMG5_int,MMG5_int,double *);
 extern int    MMG5_nortri(MMG5_pMesh mesh,MMG5_pTria pt,double *n);
 void          MMG5_printTria(MMG5_pMesh mesh,char* fileName);
 extern int    MMG5_rotmatrix(double n[3],double r[3][3]);
 int           MMG5_invmat(double *m,double *mi);
 int           MMG5_invmatg(double m[9],double mi[9]);
 int           MMG5_invmat33(double m[3][3],double mi[3][3]);
 int           MMG5_invmat22(double m[2][2],double mi[2][2]);
 int           MMG5_regnor(MMG5_pMesh mesh);
 double        MMG5_ridSizeInNormalDir(MMG5_pMesh,int,double*,MMG5_pBezier,double,double);
 double        MMG5_ridSizeInTangentDir(MMG5_pMesh, MMG5_pPoint,MMG5_int,MMG5_int*,double,double);
 int           MMG5_scale_meshAndSol(MMG5_pMesh,MMG5_pSol,MMG5_pSol,double*);
 int           MMG5_scale_scalarMetric(MMG5_pMesh, MMG5_pSol,double);
 int           MMG5_scotchCall(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol fields,MMG5_int*);
 int           MMG5_check_setted_hminhmax(MMG5_pMesh mesh);
 int           MMG5_solTruncature_iso(MMG5_pMesh mesh, MMG5_pSol met);
 int           MMG5_2dSolTruncature_ani(MMG5_pMesh mesh, MMG5_pSol met);
 int           MMG5_3dSolTruncature_ani(MMG5_pMesh mesh, MMG5_pSol met);
 int           MMG5_truncate_met3d(MMG5_pSol met, MMG5_int ip, double isqhmin, double isqhmax);
 int           MMG5_solveDefmetregSys( MMG5_pMesh, double r[3][3], double *, double *,
                                        double *, double *, double, double, double);
 int           MMG5_solveDefmetrefSys( MMG5_pMesh,MMG5_pPoint,MMG5_int*, double r[3][3],
                                        double *, double *, double *, double *,
                                        double, double, double);
 double        MMG5_surftri_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
 double        MMG5_surftri33_ani(MMG5_pMesh,MMG5_pTria,double*,double*,double*);
 double        MMG5_surftri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
 extern int    MMG5_sys33sym(double a[6], double b[3], double r[3]);
 int           MMG5_interpreg_ani(MMG5_pMesh,MMG5_pSol,MMG5_pTria,int8_t,double,double *mr);
 int           MMG5_interp_iso(double *ma,double *mb,double *mp,double t);
 int           MMG5_intersecmet22(MMG5_pMesh mesh, double *m,double *n,double *mr);
 int           MMG5_countLocalParamAtTri( MMG5_pMesh,MMG5_iNode **);
 int           MMG5_writeLocalParamAtTri( MMG5_pMesh,MMG5_iNode *,FILE*);
 double        MMG2D_quickarea(double a[2],double b[2],double c[2]);
 void          MMG5_build3DMetric(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int ip,double dbuf[6]);
 int           MMG5_loadVtuMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);
 int           MMG5_loadMshMesh_part1(MMG5_pMesh mesh,const char *filename,
                                      FILE **inm,long *posNodes, long *posElts,
                                      long **posNodeData, int *bin, int *iswp,
                                      MMG5_int *nelts,int *nsols);
 int            MMG5_check_readedMesh(MMG5_pMesh mesh,MMG5_int nref);
 int            MMG5_loadMshMesh_part2(MMG5_pMesh mesh,MMG5_pSol *sol,FILE **inm,
                                       const long posNodes,const long posElts,
                                       const long *posNodeData,const int bin,
                                       const int iswp,const MMG5_int nelts,
                                       const int nsols);
int             MMG5_saveMshMesh(MMG5_pMesh,MMG5_pSol*,const char*, int);
int             MMG5_saveDisp(MMG5_pMesh ,MMG5_pSol );
int             MMG5_loadSolHeader(const char*,int,FILE**,int*,int*,int*,MMG5_int*,
                                   int*,int*,int**,long*,int);
int             MMG5_chkMetricType(MMG5_pMesh mesh,int *type,int*, FILE *inm);
int             MMG5_readFloatSol3D(MMG5_pSol,FILE*,int,int,int);
int             MMG5_readDoubleSol3D(MMG5_pSol,FILE*,int,int,MMG5_int);
int             MMG5_saveSolHeader( MMG5_pMesh,const char*,FILE**,int,int*,MMG5_int*,MMG5_int,
                                    int,int,int*,int*,int*);
int             MMG5_saveSolAtTrianglesHeader( MMG5_pMesh,FILE *,int,int,MMG5_int*,int,
                                               int,int*,int*,int*);
int             MMG5_saveSolAtTetrahedraHeader( MMG5_pMesh,FILE *,int,int,MMG5_int*,int,
                                                int,int*,int*,int*);
void            MMG5_writeDoubleSol3D(MMG5_pMesh,MMG5_pSol,FILE*,int,MMG5_int,int);
void            MMG5_printMetStats(MMG5_pMesh mesh,MMG5_pSol met);
void            MMG5_printSolStats(MMG5_pMesh mesh,MMG5_pSol *sol);

int MMG5_defsiz_startingMessage (MMG5_pMesh,MMG5_pSol,const char * funcname );
void MMG5_gradation_info ( MMG5_pMesh );
int MMG5_sum_reqEdgeLengthsAtPoint ( MMG5_pMesh,MMG5_pSol,MMG5_int ip0,MMG5_int ip1 );
int MMG5_compute_meanMetricAtMarkedPoints_iso ( MMG5_pMesh mesh,MMG5_pSol met);
int MMG5_compute_meanMetricAtMarkedPoints_ani ( MMG5_pMesh mesh,MMG5_pSol met);

int  MMG5_reset_metricAtReqEdges_surf ( MMG5_pMesh,MMG5_pSol,int8_t );
void MMG5_mark_pointsOnReqEdge_fromTria ( MMG5_pMesh mesh );
int  MMG5_gradsiz_iso ( MMG5_pMesh mesh,MMG5_pSol met );
int  MMG5_gradsizreq_iso(MMG5_pMesh ,MMG5_pSol );
MMG5_int  MMG5_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met,int *it);
int  MMG5_gradsizreq_ani(MMG5_pMesh mesh,MMG5_pSol met);
int  MMG5_simred2d(MMG5_pMesh,double*,double*,double dm[2],double dn[2],double vp[2][2]);
int  MMG5_simred3d(MMG5_pMesh mesh,double *m,double *n,double dm[3],double dn[3],double vp[3][3]);
extern int  MMG5_updatemet2d_ani(double *m,double *n,double dm[2],double dn[2],double vp[2][2],int8_t ier );
int  MMG5_updatemet3d_ani(double *m,double *n,double dm[3],double dn[3],double vp[3][3],int8_t ier );
void MMG5_gradEigenvreq(double *dm,double *dn,double,int8_t,int8_t *);
int  MMG5_updatemetreq_ani(double *n,double dn[2],double vp[2][2]);
int    MMG5_swapbin(int sbin);
MMG5_int    MMG5_swapbin_int(MMG5_int sbin);
float  MMG5_swapf(float sbin);
double MMG5_swapd(double sbin);
int MMG5_MultiMat_init(MMG5_pMesh);
int MMG5_isLevelSet(MMG5_pMesh,MMG5_int,MMG5_int);
int MMG5_isSplit(MMG5_pMesh ,MMG5_int ,MMG5_int *,MMG5_int *);
int MMG5_isNotSplit(MMG5_pMesh ,MMG5_int);
int MMG5_getStartRef(MMG5_pMesh ,MMG5_int, MMG5_int *);
int MMG5_snpval_ls(MMG5_pMesh mesh,MMG5_pSol sol);
int MMG5_snpval_lssurf(MMG5_pMesh mesh,MMG5_pSol sol);
int MMG5_rmc(MMG5_pMesh ,MMG5_pSol );
int MMG5_resetRef_ls(MMG5_pMesh );
int MMG5_resetRef_lssurf(MMG5_pMesh );
int MMG5_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol);
int MMG5_setref_lssurf(MMG5_pMesh mesh, MMG5_pSol sol);
int MMG5_chkmaniball(MMG5_pMesh mesh, MMG5_int start, int8_t istart);
int MMG5_chkmanimesh(MMG5_pMesh mesh);

/* test functions */
extern double MMG5_test_mat_error( int8_t nelem,double m1[],double m2[] );
int MMG5_test_invmat22();
int MMG5_test_invmat33();
int MMG5_test_eigenvmatsym2d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                             double vpex[][2]);
int MMG5_test_eigenvmatnonsym2d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                                double vpex[][2],double ivpex[][2]);
int MMG5_test_eigenvmatsym3d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                             double vpex[][3]);
int MMG5_test_eigenvmatnonsym3d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                                double vpex[][3],double ivpex[][3]);
int MMG5_test_transpose3d();
int MMG5_test_dotprod();
int MMG5_test_crossprod3d();
int MMG5_test_mn();
extern int MMG5_test_rmtr();
int MMG5_test_rotmatrix();
int MMG5_test_simred2d(MMG5_pMesh mesh,double *mex,double *nex,double *dmex,double *dnex,double vpex[][2]);
int MMG5_test_simred3d(MMG5_pMesh mesh,double *mex,double *nex,double *dmex,double *dnex,double vpex[][3]);
int MMG5_test_updatemet2d_ani();
int MMG5_test_updatemet3d_ani();
int MMG5_test_intersecmet22(MMG5_pMesh mesh);
int MMG5_test_intersecmet33(MMG5_pMesh mesh);

/* tools */
void MMG5_mark_verticesAsUnused ( MMG5_pMesh mesh );
void MMG5_mark_usedVertices ( MMG5_pMesh mesh,void (*delPt)(MMG5_pMesh,MMG5_int) );
void MMG5_keep_subdomainElts ( MMG5_pMesh,int,int (*delElt)(MMG5_pMesh,MMG5_int) );

void   MMG5_Set_commonFunc(void);

#ifdef __cplusplus
}
#endif

#endif
