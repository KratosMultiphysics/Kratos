/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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
#ifndef LIBMMG2D_PRIVATE_H
#define LIBMMG2D_PRIVATE_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>


#include "libmmgcommon_private.h"

#ifdef __cplusplus
extern "C" {
#endif

/* constantes */

#define M_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define M_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define M_LAMBDA  0.34
#define M_MU      0.33

#define MMG2D_EPSD   1.e-10 //e-20??
#define MMG2D_EPSA   1.e-12

#define MMG2D_PRECI  1.
#define MMG2D_SIZE   0.75
#define MMG2D_ALPHA  0.28867513459
#define MMG2D_ALPHAD 3.464101615137755   /* 6.0 / sqrt(3.0)  */
#define MMG2D_BADKAL    0.2
#define MMG2D_NULKAL    1.e-6
#define MMG2D_ANGCORN   -1.e-6
#define MMG2D_SHORTMAX     0x7fff

#define MMG2D_LLONG  2.0
#define MMG2D_LSHRT  0.3
#define MMG2D_LOPTL      1.4
#define MMG2D_LOPTS     0.71

#define MMG2D_NPMAX   50000
#define MMG2D_NEDMAX  100000
#define MMG2D_NEMAX   100000

/** \brief idir[i]: vertices of edge i for a quad */
static const uint8_t MMG2D_idir_q[4][2] = { {0,1},{0,3},{1,2},{2,3} };

/** Free allocated pointers of mesh and sol structure and return value val */
#define MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,val)do               \
  {                                                                 \
    if ( !MMG2D_Free_all(MMG5_ARG_start,                            \
                         MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met, \
                         MMG5_ARG_ppLs,&ls,MMG5_ARG_ppDisp,&disp,   \
                         MMG5_ARG_end) ) {                          \
      return MMG5_LOWFAILURE;                                       \
    }                                                               \
    return val;                                                    \
  }while(0)

/**
 * \param sigid signal number.
 *
 * Signal handling: specify error messages depending from catched signal.
 *
 */
static inline
void MMG2D_excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
  case SIGABRT:
    fprintf(stdout,"  Abnormal stop\n"); break;
  case SIGFPE:
    fprintf(stdout,"  Floating-point exception\n"); break;
  case SIGILL:
    fprintf(stdout,"  Illegal instruction\n"); break;
  case SIGSEGV:
    fprintf(stdout,"  Segmentation fault\n"); break;
  case SIGTERM:
  case SIGINT:
    fprintf(stdout,"  Program killed\n"); break;
  }
  exit(EXIT_FAILURE);
}

static const int MMG2D_iare[3][2] = {{1,2},{2,0},{0,1}};
static const unsigned int MMG2D_idir[5] = {0,1,2,0,1};

/** Reallocation of point table and sol table and creation
    of point ip with coordinates o and tag tag*/
#define MMG2D_POINT_REALLOC(mesh,sol,ip,wantedGap,law,o,tag ) do       \
  {                                                                     \
    MMG5_int klink;                                                     \
                                                                        \
    assert ( mesh && mesh->point );                                     \
    MMG5_TAB_RECALLOC(mesh,mesh->point,mesh->npmax,wantedGap,MMG5_Point, \
                       "larger point table",law);                       \
                                                                        \
    mesh->npnil = mesh->np+1;                                           \
    for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)               \
      mesh->point[klink].tmp  = klink+1;                                \
                                                                        \
    /* solution */                                                      \
    if ( sol ) {                                                        \
      if ( sol->m ) {                                                   \
        MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax-sol->npmax))*sizeof(double), \
                     "larger solution",law);                            \
        MMG5_SAFE_REALLOC(sol->m,sol->size*(sol->npmax+1),              \
                          sol->size*(mesh->npmax+1),                    \
                          double,"larger solution",law);                \
      }                                                                 \
      sol->npmax = mesh->npmax;                                         \
    }                                                                   \
                                                                        \
    /* We try again to add the point */                                 \
    ip = MMG2D_newPt(mesh,o,tag);                                      \
    if ( !ip ) {law;}                                                   \
  }while(0)

/** Reallocation of tria table and creation
    of tria jel */
#define MMG2D_TRIA_REALLOC(mesh,jel,wantedGap,law ) do                  \
  {                                                                     \
   MMG5_int klink,oldSiz;                                               \
                                                                        \
   oldSiz = mesh->ntmax;                                                \
                                                                        \
  int max_factor;                                                       \
  if ( mesh->nquad ) {                                                  \
    /* If mesh contains quads, we need to compute 4*tet_quad to hash quads when
     * packing the mesh so we can't create a mesh larger than INT32_MAX/4 */ \
    max_factor = 4;                                                     \
  }                                                                     \
  else {                                                                \
    /* With only tria, maximal number of tetra is INT32_MAX/3 */        \
    max_factor = 3;                                                     \
  }                                                                     \
  MMG5_CHK_INT32_OVERFLOW(wantedGap,oldSiz,max_factor,max_factor+2,law);\
                                                                        \
  MMG5_TAB_RECALLOC(mesh,mesh->tria,mesh->ntmax,wantedGap,MMG5_Tria,    \
                    "larger tria table",law);                           \
                                                                        \
   mesh->nenil = mesh->nt+1;                                            \
   for (klink=mesh->nenil; klink<mesh->ntmax-1; klink++)                \
     mesh->tria[klink].v[2]  = klink+1;                                 \
                                                                        \
   if ( mesh->adja ) {                                                  \
     /* adja table */                                                   \
     MMG5_ADD_MEM(mesh,3*(mesh->ntmax-oldSiz)*sizeof(MMG5_int),         \
                   "larger adja table",law);                            \
     MMG5_SAFE_RECALLOC(mesh->adja,3*oldSiz+5,3*mesh->ntmax+5,MMG5_int  \
                         ,"larger adja table",law);                     \
   }                                                                    \
                                                                        \
   /* We try again to add the point */                                  \
   jel = MMG2D_newElt(mesh);                                            \
   if ( !jel ) {law;}                                                   \
   }while(0)

/* Prototypes */
/*zaldy*/
MMG5_int MMG2D_newPt(MMG5_pMesh mesh,double c[2],int16_t tag);
void MMG2D_delPt(MMG5_pMesh mesh,MMG5_int ip) ;
void MMG5_delEdge(MMG5_pMesh mesh,MMG5_int iel);
MMG5_int MMG2D_newElt(MMG5_pMesh mesh);
int  MMG2D_delElt(MMG5_pMesh mesh,MMG5_int iel);
int MMG2D_zaldy(MMG5_pMesh mesh);
size_t MMG5_memSize(void);
int MMG2D_memOption(MMG5_pMesh mesh);
int  MMG2D_setMeshSize_alloc(MMG5_pMesh);

int MMG2D_pack(MMG5_pMesh ,MMG5_pSol, MMG5_pSol );
int MMG2D_outqua(MMG5_pMesh ,MMG5_pSol );
int MMG2D_mmg2d1(MMG5_pMesh ,MMG5_pSol );

int  MMG2D_Init_mesh_var( va_list argptr );
int  MMG2D_Free_all_var( va_list argptr );
int  MMG2D_Free_structures_var( va_list argptr );
int  MMG2D_Free_names_var( va_list argptr );

int MMG2D_mmg2d2(MMG5_pMesh , MMG5_pSol);
int MMG2D_mmg2d6(MMG5_pMesh ,MMG5_pSol,MMG5_pSol );
int MMG2D_mmg2d9(MMG5_pMesh ,MMG5_pSol ,MMG5_pSol,MMG5_int** );
int MMG2D_swapdelone(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int8_t ,double ,MMG5_int *);
int MMG5_mmg2dChkmsh(MMG5_pMesh , int, MMG5_int );
int MMG2D_2dMeshCheck(MMG5_pMesh mesh);
int MMG2D_boulep(MMG5_pMesh , MMG5_int , int , MMG5_int * );
int MMG2D_prilen(MMG5_pMesh ,MMG5_pSol );

int MMG2D_coorbary(MMG5_pMesh ,MMG5_pTria ,double c[2],double* ,double* ,double* );
MMG5_int MMG2D_isInTriangle(MMG5_pMesh ,MMG5_int,double c[2]);
int MMG2D_cutEdge(MMG5_pMesh ,MMG5_pTria ,MMG5_pPoint ,MMG5_pPoint );
int MMG2D_cutEdgeTriangle(MMG5_pMesh ,MMG5_int ,MMG5_int ,MMG5_int );
MMG5_int MMG2D_findTria(MMG5_pMesh ,MMG5_int );
int MMG2D_locateEdge(MMG5_pMesh ,MMG5_int ,MMG5_int ,MMG5_int* ,MMG5_int* ) ;
int MMG2D_bdryenforcement(MMG5_pMesh ,MMG5_pSol);
int MMG2D_settagtriangles(MMG5_pMesh ,MMG5_pSol );
MMG5_int MMG2D_findtrianglestate(MMG5_pMesh ,MMG5_int ,MMG5_int ,MMG5_int ,MMG5_int ,MMG5_int ,MMG5_int );

int MMG2D_baseBdry(MMG5_pMesh mesh);

int MMG2D_cavity(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int *);
int MMG2D_delone(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int *,int );
int MMG2D_cenrad_iso(MMG5_pMesh ,double *,double *,double *);

/* Adds Charles */
double MMG2D_caltri_iso_3pt(double *a,double *b,double *c);
int MMG2D_hashTria(MMG5_pMesh );
int MMG2D_hashQuad(MMG5_pMesh mesh);
int MMG2D_cuttri(MMG5_pMesh ,MMG5_pSol,MMG5_pSol );
int MMG2D_split1_sim(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int vx[3]);
int MMG2D_split2_sim(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int vx[3]);
int MMG2D_split3_sim(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int vx[3]);
int MMG2D_split1(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int vx[3]);
int MMG2D_split2(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int vx[3]);
int MMG2D_split3(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int vx[3]);
int MMG2D_splitbar(MMG5_pMesh ,MMG5_int ,MMG5_int );
int MMG2D_assignEdge(MMG5_pMesh );
int MMG2D_bdryEdge(MMG5_pMesh );
int MMG2D_setadj(MMG5_pMesh,int8_t );
int MMG2D_singul(MMG5_pMesh,MMG5_int );
int MMG2D_analys(MMG5_pMesh );
int MMG2D_norver(MMG5_pMesh,MMG5_int );
int MMG2D_regnor(MMG5_pMesh );
int MMG2D_regver(MMG5_pMesh );
int MMG2D_boulen(MMG5_pMesh , MMG5_int ,int8_t ,MMG5_int *,MMG5_int *,double *);
int MMG2D_mmg2d1n(MMG5_pMesh ,MMG5_pSol );
int MMG2D_anatri(MMG5_pMesh ,MMG5_pSol ,int8_t );
int MMG2D_adptri(MMG5_pMesh ,MMG5_pSol );
int MMG2D_defsiz_iso(MMG5_pMesh ,MMG5_pSol );
int MMG2D_defsiz_ani(MMG5_pMesh ,MMG5_pSol );
int MMG2D_defmetbdy_2d(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int8_t );
int MMG2D_defaultmet_2d(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int8_t );
MMG5_int MMG2D_grad2met_ani(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria,MMG5_int,MMG5_int);
int MMG2D_grad2metreq_ani(MMG5_pMesh ,MMG5_pSol ,MMG5_pTria,MMG5_int,MMG5_int);
int MMG2D_gradsiz_ani(MMG5_pMesh ,MMG5_pSol );
int MMG2D_gradsizreq_ani(MMG5_pMesh ,MMG5_pSol );
MMG5_int MMG2D_anaelt(MMG5_pMesh ,MMG5_pSol ,int );
MMG5_int MMG2D_colelt(MMG5_pMesh ,MMG5_pSol ,int );
MMG5_int MMG2D_swpmsh(MMG5_pMesh ,MMG5_pSol ,int );
double MMG2D_lencurv_iso(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int );
double MMG2D_lencurv_ani(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int );
int MMG2D_chkedg(MMG5_pMesh ,MMG5_int );
int MMG2D_bezierCurv(MMG5_pMesh ,MMG5_int ,int8_t ,double ,double *,double *);
int MMG2D_dichoto(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,MMG5_int *);
double MMG2D_quickcal(MMG5_pMesh , MMG5_pTria );
int MMG2D_chkcol(MMG5_pMesh,MMG5_pSol,MMG5_int,int8_t,MMG5_int *,int8_t);
int MMG2D_colver(MMG5_pMesh,int,MMG5_int*);
int MMG2D_colver3(MMG5_pMesh,MMG5_int*);
int MMG2D_colver2(MMG5_pMesh,MMG5_int*);
int MMG2D_boulet(MMG5_pMesh,MMG5_int,int8_t,MMG5_int*);
int MMG2D_bouleendp(MMG5_pMesh,MMG5_int,int8_t,MMG5_int*,MMG5_int*,MMG5_int*);
int MMG2D_savemesh_db(MMG5_pMesh ,char* ,int8_t );
int MMG2D_savemet_db(MMG5_pMesh ,MMG5_pSol ,char* ,int8_t );
int MMG2D_chkswp(MMG5_pMesh , MMG5_pSol ,MMG5_int ,int8_t ,int8_t );
int MMG2D_swapar(MMG5_pMesh ,MMG5_int ,int8_t );
int MMG5_interpmet22(MMG5_pMesh ,double *,double *,double ,double *);
int MMG2D_intmet_iso(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int8_t ,MMG5_int ,double );
int MMG2D_intmet_ani(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int8_t ,MMG5_int ,double );
MMG5_int MMG2D_adpspl(MMG5_pMesh ,MMG5_pSol );
int MMG2D_adpcol(MMG5_pMesh ,MMG5_pSol );
MMG5_int MMG2D_movtri(MMG5_pMesh ,MMG5_pSol ,int ,int8_t );
MMG5_int MMG2D_chkspl(MMG5_pMesh ,MMG5_pSol ,MMG5_int ,int8_t );
int MMG2D_split1b(MMG5_pMesh ,MMG5_int ,int8_t ,MMG5_int );
int MMG2D_movedgpt(MMG5_pMesh ,MMG5_pSol ,int ,MMG5_int *,int8_t );
int MMG2D_movintpt(MMG5_pMesh ,MMG5_pSol ,int ,MMG5_int *,int8_t );
int MMG2D_movintpt_ani(MMG5_pMesh ,MMG5_pSol ,int ,MMG5_int *,int8_t );
int MMG2D_savenor_db(MMG5_pMesh ,char*,int8_t );
int MMG2D_savedisp_db(MMG5_pMesh mesh,MMG5_pSol ,char*,int8_t );
int MMG2D_velextLS(MMG5_pMesh ,MMG5_pSol );

/* tools */
void MMG2D_keep_only1Subdomain ( MMG5_pMesh mesh,int nsd );

/* useful functions to debug */
MMG5_int  MMG2D_indElt(MMG5_pMesh mesh,MMG5_int kel);
MMG5_int  MMG2D_indPt(MMG5_pMesh mesh,MMG5_int kp);

/* Management of local parameters */
int MMG2D_freeLocalPar(MMG5_pMesh );

/* functions pointers */
int    MMG2D_doSol_ani(MMG5_pMesh mesh,MMG5_pSol sol);
int    MMG2D_doSol_iso(MMG5_pMesh mesh,MMG5_pSol sol);
double long_ani(double *ca,double *cb,double *ma,double *mb);
double long_iso(double *ca,double *cb,double *ma,double *mb);
double MMG2D_caltri_ani(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pTria );
double MMG2D_caltri_iso(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pTria );
int    interp_ani(double *,double *,double * ,double );
int    interp_iso(double *,double *,double * ,double );
int    lissmet_iso(MMG5_pMesh mesh,MMG5_pSol sol);
int    lissmet_ani(MMG5_pMesh mesh,MMG5_pSol sol);
int    MMG2D_sum_reqEdgeLengthsAtPoint(MMG5_pMesh,MMG5_pSol,MMG5_pTria,int8_t);
int    MMG2D_set_metricAtPointsOnReqEdges(MMG5_pMesh,MMG5_pSol,int8_t);

#ifdef __cplusplus
}
#endif

#endif
