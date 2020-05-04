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

#ifndef MMGS_H
#define MMGS_H

#include "libmmgs.h"
#include "mmgcommon.h"

#ifdef __cplusplus
extern "C" {
#endif

/* numerical accuracy */
#define MMGS_ALPHAD    3.464101615137755   /* 6.0 / sqrt(3.0)  */

#define MMGS_LOPTL     1.4
#define MMGS_LOPTS     0.71
#define MMGS_LLONG     2.0
#define MMGS_LSHRT     0.3

#define MMGS_LMAX  1024
#define MMGS_BADKAL      2.e-2
#define MMGS_NULKAL      1.e-4

#define MMGS_NPMAX     500000
#define MMGS_NTMAX    1000000
#define MMGS_XPMAX     500000


#define MS_SIN(tag)      ((tag & MG_CRN) || (tag & MG_REQ) || (tag & MG_NOM))


/** Free allocated pointers of mesh and sol structure and return value val */
#define MMGS_RETURN_AND_FREE(mesh,met,val)do                       \
  {                                                                 \
    if ( !MMGS_Free_all(MMG5_ARG_start,                             \
                        MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,  \
                        MMG5_ARG_end) ) {                           \
      return MMG5_LOWFAILURE;                                       \
    }                                                               \
    return val;                                                    \
  }while(0)

/** Reallocation of point table and sol table and creation
    of point ip with coordinates o and tag tag*/
#define MMGS_POINT_REALLOC(mesh,sol,ip,wantedGap,law,o,tag ) do        \
  {                                                                     \
    int klink;                                                          \
                                                                        \
    MMG5_TAB_RECALLOC(mesh,mesh->point,mesh->npmax,wantedGap,MMG5_Point, \
                       "larger point table",law);                       \
                                                                        \
    mesh->npnil = mesh->np+1;                                           \
    for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)               \
      mesh->point[klink].tmp  = klink+1;                                \
                                                                        \
    /* solution */                                                      \
    if ( sol->m ) {                                                     \
      MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax-sol->npmax))*sizeof(double), \
                    "larger solution",law);                             \
      MMG5_SAFE_REALLOC(sol->m,sol->size*(sol->npmax+1),               \
                         sol->size*(mesh->npmax+1),double,              \
                         "larger solution",law);                        \
    }                                                                   \
    sol->npmax = mesh->npmax;                                           \
                                                                        \
    /* We try again to add the point */                                 \
    ip = MMGS_newPt(mesh,o,tag);                                       \
    if ( !ip ) {law;}                                                   \
  }while(0)

/** Reallocation of tria table and creation
    of tria jel */
#define MMGS_TRIA_REALLOC( mesh,jel,wantedGap,law ) do                 \
  {                                                                     \
    int klink,oldSiz;                                                   \
                                                                        \
    oldSiz = mesh->ntmax;                                               \
    MMG5_TAB_RECALLOC(mesh,mesh->tria,mesh->ntmax,wantedGap,MMG5_Tria, \
                       "larger tria table",law);                        \
                                                                        \
    mesh->nenil = mesh->nt+1;                                           \
    for (klink=mesh->nenil; klink<mesh->ntmax-1; klink++)               \
      mesh->tria[klink].v[2]  = klink+1;                                \
                                                                        \
    if ( mesh->adja ) {                                                 \
      /* adja table */                                                  \
      MMG5_ADD_MEM(mesh,3*(mesh->ntmax-oldSiz)*sizeof(int),            \
                    "larger adja table",law);                           \
      MMG5_SAFE_RECALLOC(mesh->adja,3*oldSiz+5,3*mesh->ntmax+5,int     \
                          ,"larger adja table",law);                    \
    }                                                                   \
                                                                        \
    /* We try again to add the point */                                 \
    jel = MMGS_newElt(mesh);                                           \
    if ( !jel ) {law;}                                                  \
  }while(0)

/* prototypes */
int  MMGS_Init_mesh_var( va_list argptr );
int  MMGS_Free_all_var( va_list argptr );
int  MMGS_Free_structures_var( va_list argptr );
int  MMGS_Free_names_var( va_list argptr );

int  MMGS_zaldy(MMG5_pMesh mesh);
int  assignEdge(MMG5_pMesh mesh);
int  MMGS_analys(MMG5_pMesh mesh);
int  MMGS_inqua(MMG5_pMesh,MMG5_pSol);
int  MMGS_outqua(MMG5_pMesh,MMG5_pSol);
int  MMGS_hashTria(MMG5_pMesh );
int  curvpo(MMG5_pMesh ,MMG5_pSol );
int  MMG5_mmgs1(MMG5_pMesh ,MMG5_pSol );
int  MMGS_mmgs2(MMG5_pMesh ,MMG5_pSol );
int  boulet(MMG5_pMesh mesh,int start,int ip,int *list);
int  boulechknm(MMG5_pMesh mesh,int start,int ip,int *list);
int  boulep(MMG5_pMesh mesh,int start,int ip,int *list);
int  bouletrid(MMG5_pMesh mesh,int start,int ip,int *il1,int *l1,int *il2,int *l2,int *ip0,int *ip1);
int  MMGS_newPt(MMG5_pMesh mesh,double c[3],double n[3]);
void MMGS_delPt(MMG5_pMesh mesh,int ip);
int  MMGS_newElt(MMG5_pMesh mesh);
int  MMGS_delElt(MMG5_pMesh mesh,int iel);
int  chkedg(MMG5_pMesh ,int );
int  MMG5_mmgsBezierCP(MMG5_pMesh ,MMG5_Tria*, MMG5_pBezier, char ori);
int  MMGS_bezierInt(MMG5_pBezier ,double *,double *,double *,double *);
int  MMGS_simbulgept(MMG5_pMesh mesh,MMG5_pSol met, int k,int i,int ip);
int  MMGS_split1_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int i, int *vx);
int  MMG5_split2_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  MMGS_split3_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  MMGS_split1(MMG5_pMesh mesh,MMG5_pSol met,int k,int i,int *vx);
int  MMGS_split2(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  MMGS_split3(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  split1b(MMG5_pMesh mesh,int k,char i,int ip);
int  chkcol(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int *list,char typchk);
int  colver(MMG5_pMesh mesh,int *list,int ilist);
int  colver3(MMG5_pMesh mesh,int*list);
int  colver2(MMG5_pMesh mesh,int *ilist);
int  swapar(MMG5_pMesh mesh,int k,int i);
int  chkswp(MMG5_pMesh mesh,MMG5_pSol met,int k,int i,char typchk);
int  swpedg(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist,char typchk);
char typelt(MMG5_pPoint p[3],char *ia);
int  litswp(MMG5_pMesh mesh,int k,char i,double kal);
int  litcol(MMG5_pMesh mesh,int k,char i,double kal);
int  MMG5_mmgsChkmsh(MMG5_pMesh,int,int);
int  paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]);
int  intregmet(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,double s,double mr[6]);
int  MMG5_intridmet(MMG5_pMesh,MMG5_pSol,int,int,double,double*,double*);
int  setref(MMG5_pMesh,int,int,int);
int  delref(MMG5_pMesh);
int  chkmet(MMG5_pMesh,MMG5_pSol);
int  chknor(MMG5_pMesh);
size_t MMG5_memSize(void);
int MMGS_memOption(MMG5_pMesh mesh);
int MMGS_setMeshSize_alloc( MMG5_pMesh mesh );

#ifdef USE_SCOTCH
int MMG5_mmgsRenumbering(int vertBoxNbr, MMG5_pMesh mesh, MMG5_pSol sol);
#endif

/* useful functions to debug */
int  MMGS_indElt(MMG5_pMesh mesh,int kel);
int  MMGS_indPt(MMG5_pMesh mesh,int kp);

/* function pointers */
/* init structures */
void  MMG5_Init_parameters(MMG5_pMesh mesh);
/* iso/aniso computations */
double caleltsig_ani(MMG5_pMesh mesh,MMG5_pSol met,int iel);
double caleltsig_iso(MMG5_pMesh mesh,MMG5_pSol met,int iel);
int    MMGS_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met);
int    MMGS_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met);
void   MMG5_defaultValues(MMG5_pMesh);
int    MMGS_gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met);
int    MMGS_gradsizreq_ani(MMG5_pMesh mesh,MMG5_pSol met);
int    intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    MMGS_intmet33_ani(MMG5_pMesh,MMG5_pSol,int,char,int,double);
int    movridpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movintpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movridpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    MMGS_prilen(MMG5_pMesh mesh,MMG5_pSol met,int);
int    MMGS_set_metricAtPointsOnReqEdges ( MMG5_pMesh,MMG5_pSol );

extern double (*MMG5_calelt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
extern int    (*MMGS_defsiz)(MMG5_pMesh mesh,MMG5_pSol met);
extern int    (*MMGS_gradsiz)(MMG5_pMesh mesh,MMG5_pSol met);
extern int    (*MMGS_gradsizreq)(MMG5_pMesh mesh,MMG5_pSol met);
extern int    (*intmet)(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
extern int    (*movridpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
extern int    (*movintpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);

/**
 * Set common pointer functions between mmgs and mmg3d to the matching mmgs
 * functions.
 */
static inline
void MMGS_Set_commonFunc() {
  MMG5_bezierCP          = MMG5_mmgsBezierCP;
  MMG5_chkmsh            = MMG5_mmgsChkmsh;
  MMG5_indPt             = MMGS_indPt;
  MMG5_indElt            = MMGS_indElt;
  MMG5_grad2met_ani      = MMG5_grad2metSurf;
  MMG5_grad2metreq_ani   = MMG5_grad2metSurfreq;
#ifdef USE_SCOTCH
  MMG5_renumbering       = MMG5_mmgsRenumbering;
#endif
}

#ifdef __cplusplus
}
#endif

#endif
