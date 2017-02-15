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

#ifndef _MMGS_H
#define _MMGS_H

#include "libmmgs.h"
#include "mmgcommon.h"

#ifdef __cplusplus
extern "C" {
#endif

/* numerical accuracy */
#define ALPHAD    3.464101615137755   /* 6.0 / sqrt(3.0)  */

#define LOPTL     1.4
#define LOPTS     0.71
#define LLONG     2.0
#define LSHRT     0.3

#define _MMG5_LMAX  1024
#define BADKAL      2.e-2
#define NULKAL      1.e-4

#define _MMG5_NPMAX     500000
#define _MMG5_NTMAX    1000000
#define _MMG5_XPMAX     500000


#define MS_SIN(tag)      ((tag & MG_CRN) || (tag & MG_REQ) || (tag & MG_NOM))


/** Free allocated pointers of mesh and sol structure and return value val */
#define _MMGS_RETURN_AND_FREE(mesh,met,val)do                 \
  {                                                           \
    MMGS_Free_all(MMG5_ARG_start,                             \
                  MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,  \
                  MMG5_ARG_end);                              \
    return(val);                                              \
  }while(0)

/** Reallocation of point table and sol table and creation
    of point ip with coordinates o and tag tag*/
#define _MMGS_POINT_REALLOC(mesh,sol,ip,wantedGap,law,o,tag ) do        \
  {                                                                     \
    int klink;                                                          \
                                                                        \
    _MMG5_TAB_RECALLOC(mesh,mesh->point,mesh->npmax,wantedGap,MMG5_Point, \
                       "larger point table",law);                       \
                                                                        \
    mesh->npnil = mesh->np+1;                                           \
    for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)               \
      mesh->point[klink].tmp  = klink+1;                                \
                                                                        \
    /* solution */                                                      \
    if ( sol->m ) {                                                     \
      _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax-sol->npmax))*sizeof(double), \
                    "larger solution",law);                             \
      _MMG5_SAFE_REALLOC(sol->m,sol->size*(mesh->npmax+1),double,"larger solution"); \
    }                                                                   \
    sol->npmax = mesh->npmax;                                           \
                                                                        \
    /* We try again to add the point */                                 \
    ip = _MMGS_newPt(mesh,o,tag);                                       \
    if ( !ip ) {law;}                                                   \
  }while(0)

/** Reallocation of tria table and creation
    of tria jel */
#define _MMGS_TRIA_REALLOC( mesh,jel,wantedGap,law ) do                 \
  {                                                                     \
    int klink,oldSiz;                                                   \
                                                                        \
    oldSiz = mesh->ntmax;                                               \
    _MMG5_TAB_RECALLOC(mesh,mesh->tria,mesh->ntmax,wantedGap,MMG5_Tria, \
                       "larger tria table",law);                        \
                                                                        \
    mesh->nenil = mesh->nt+1;                                           \
    for (klink=mesh->nenil; klink<mesh->ntmax-1; klink++)               \
      mesh->tria[klink].v[2]  = klink+1;                                \
                                                                        \
    if ( mesh->adja ) {                                                 \
      /* adja table */                                                  \
      _MMG5_ADD_MEM(mesh,3*(mesh->ntmax-oldSiz)*sizeof(int),            \
                    "larger adja table",law);                           \
      _MMG5_SAFE_RECALLOC(mesh->adja,3*mesh->nt+5,3*mesh->ntmax+5,int   \
                          ,"larger adja table");                        \
    }                                                                   \
                                                                        \
    /* We try again to add the point */                                 \
    jel = _MMGS_newElt(mesh);                                           \
    if ( !jel ) {law;}                                                  \
  }while(0)

/* prototypes */
void _MMGS_Init_mesh_var( va_list argptr );
void _MMGS_Free_all_var( va_list argptr );
void _MMGS_Free_structures_var( va_list argptr );
void _MMGS_Free_names_var( va_list argptr );

int  _MMGS_zaldy(MMG5_pMesh mesh);
int  assignEdge(MMG5_pMesh mesh);
int  _MMGS_analys(MMG5_pMesh mesh);
int  _MMGS_inqua(MMG5_pMesh,MMG5_pSol);
int  _MMGS_outqua(MMG5_pMesh,MMG5_pSol);
int  _MMGS_hashTria(MMG5_pMesh );
int  curvpo(MMG5_pMesh ,MMG5_pSol );
int  _MMG5_mmgs1(MMG5_pMesh ,MMG5_pSol );
int  _MMGS_mmgs2(MMG5_pMesh ,MMG5_pSol );
int  boulet(MMG5_pMesh mesh,int start,int ip,int *list);
int  boulechknm(MMG5_pMesh mesh,int start,int ip,int *list);
int  boulep(MMG5_pMesh mesh,int start,int ip,int *list);
int  bouletrid(MMG5_pMesh mesh,int start,int ip,int *il1,int *l1,int *il2,int *l2,int *ip0,int *ip1);
int  _MMGS_newPt(MMG5_pMesh mesh,double c[3],double n[3]);
void _MMGS_delPt(MMG5_pMesh mesh,int ip);
int  _MMGS_newElt(MMG5_pMesh mesh);
void _MMGS_delElt(MMG5_pMesh mesh,int iel);
int  chkedg(MMG5_pMesh ,int );
int  _MMG5_mmgsBezierCP(MMG5_pMesh ,MMG5_Tria*, _MMG5_pBezier, char ori);
int  _MMGS_bezierInt(_MMG5_pBezier ,double *,double *,double *,double *);
int  _MMGS_simbulgept(MMG5_pMesh mesh,MMG5_pSol met, int k,int i,int ip);
int  _MMGS_split1_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int i, int *vx);
int  _MMG5_split2_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  _MMGS_split3_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  _MMGS_split1(MMG5_pMesh mesh,MMG5_pSol met,int k,int i,int *vx);
int  _MMGS_split2(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  _MMGS_split3(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
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
int  _MMG5_mmgsChkmsh(MMG5_pMesh,int,int);
int  paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]);
int  intregmet(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,double s,double mr[6]);
int  _MMG5_intridmet(MMG5_pMesh,MMG5_pSol,int,int,double,double*,double*);
int  setref(MMG5_pMesh,int,int,int);
int  delref(MMG5_pMesh);
int  chkmet(MMG5_pMesh,MMG5_pSol);
int  chknor(MMG5_pMesh);
long long _MMG5_memSize(void);
void _MMGS_memOption(MMG5_pMesh mesh);

#ifdef USE_SCOTCH
int _MMG5_mmgsRenumbering(int vertBoxNbr, MMG5_pMesh mesh, MMG5_pSol sol);
#endif

/* useful functions to debug */
int  _MMGS_indElt(MMG5_pMesh mesh,int kel);
int  _MMGS_indPt(MMG5_pMesh mesh,int kp);

/* function pointers */
/* init structures */
void  _MMG5_Init_parameters(MMG5_pMesh mesh);
/* iso/aniso computations */
extern double caleltsig_ani(MMG5_pMesh mesh,MMG5_pSol met,int iel);
extern double caleltsig_iso(MMG5_pMesh mesh,MMG5_pSol met,int iel);
int    _MMGS_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met);
int    _MMGS_defsiz_ani(MMG5_pMesh mesh,MMG5_pSol met);
void   _MMG5_defaultValues(MMG5_pMesh);
int    gradsiz_iso(MMG5_pMesh mesh,MMG5_pSol met);
int    gradsiz_ani(MMG5_pMesh mesh,MMG5_pSol met);
int    intmet_iso(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    intmet_ani(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    _MMGS_intmet33_ani(MMG5_pMesh,MMG5_pSol,int,char,int,double);
int    movridpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movintpt_iso(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movridpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    movintpt_ani(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    _MMGS_prilen(MMG5_pMesh mesh,MMG5_pSol met,int);

double (*_MMG5_calelt)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
int    (*_MMG5_defsiz)(MMG5_pMesh mesh,MMG5_pSol met);
int    (*gradsiz)(MMG5_pMesh mesh,MMG5_pSol met);
int    (*intmet)(MMG5_pMesh mesh,MMG5_pSol met,int k,char i,int ip,double s);
int    (*movridpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);
int    (*movintpt)(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ilist);

/**
 * Set common pointer functions between mmgs and mmg3d to the matching mmgs
 * functions.
 */
static inline
void _MMGS_Set_commonFunc() {
  _MMG5_bezierCP          = _MMG5_mmgsBezierCP;
  _MMG5_chkmsh            = _MMG5_mmgsChkmsh;
#ifdef USE_SCOTCH
  _MMG5_renumbering       = _MMG5_mmgsRenumbering;
#endif
}

#ifdef __cplusplus
}
#endif

#endif
