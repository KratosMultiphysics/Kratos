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

#ifndef MMG3D_H
#define MMG3D_H

#include "libmmg3d.h"
#include "libmmgcommon.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Free allocated pointers of mesh and sol structure and return value val */
#define MMG5_RETURN_AND_FREE(mesh,met,disp,val)do                  \
  {                                                                 \
    if ( !MMG3D_Free_all(MMG5_ARG_start,                            \
                         MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met, \
                         MMG5_ARG_ppDisp,&disp,                     \
                         MMG5_ARG_end) ) {                          \
      return MMG5_LOWFAILURE;                                      \
    }                                                               \
    return val;                                                    \
  }while(0)

/** Reallocation of point table and sol table and creation
    of point ip with coordinates o and tag tag*/
#define MMG3D_POINT_REALLOC(mesh,sol,ip,wantedGap,law,o,tag ) do        \
  {                                                                     \
  int klink;                                                            \
  int oldnpmax = mesh->npmax;                                           \
                                                                        \
  MMG5_TAB_RECALLOC(mesh,mesh->point,mesh->npmax,wantedGap,MMG5_Point,  \
                    "larger point table",law);                          \
                                                                        \
  mesh->npnil = mesh->np+1;                                             \
  for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)                 \
    mesh->point[klink].tmp  = klink+1;                                  \
                                                                        \
  /* solution */                                                        \
  if ( sol ) {                                                          \
    if ( sol->m ) {                                                     \
      MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax-sol->npmax))*sizeof(double), \
                   "larger solution",                                   \
                   MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,); \
                   mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point); \
                   mesh->npmax = oldnpmax;                              \
                   mesh->np = mesh->npmax-1;                            \
                   mesh->npnil = 0;                                     \
                   law);                                                \
      MMG5_SAFE_REALLOC(sol->m,sol->size*(sol->npmax+1),                \
                        sol->size*(mesh->npmax+1),                      \
                        double,"larger solution",                       \
                        MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,); \
                        mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point); \
                        mesh->npmax = oldnpmax;                         \
                        mesh->np = mesh->npmax-1;                       \
                        mesh->npnil = 0;                                \
                        law);                                           \
    }                                                                   \
    sol->npmax = mesh->npmax;                                           \
  }                                                                     \
                                                                        \
  /* We try again to add the point */                                   \
  ip = MMG3D_newPt(mesh,o,tag);                                         \
  if ( !ip ) { law; }                                                   \
  }while(0)


/** Reallocation of tetra table and creation
    of tetra jel */
#define MMG3D_TETRA_REALLOC(mesh,jel,wantedGap,law ) do                \
  {                                                                     \
    int klink,oldSiz;                                                   \
                                                                        \
    oldSiz = mesh->nemax;                                               \
    MMG5_TAB_RECALLOC(mesh,mesh->tetra,mesh->nemax,wantedGap,MMG5_Tetra, \
                       "larger tetra table",law);                       \
                                                                        \
    mesh->nenil = mesh->ne+1;                                           \
    for (klink=mesh->nenil; klink<mesh->nemax-1; klink++)               \
      mesh->tetra[klink].v[3]  = klink+1;                               \
                                                                        \
    if ( mesh->adja ) {                                                 \
      /* adja table */                                                  \
      MMG5_ADD_MEM(mesh,4*(mesh->nemax-oldSiz)*sizeof(int),            \
                    "larger adja table",law);                           \
      MMG5_SAFE_RECALLOC(mesh->adja,4*oldSiz+5,4*mesh->nemax+5,int     \
                          ,"larger adja table",law);                    \
    }                                                                   \
                                                                        \
    /* We try again to add the point */                                 \
    jel = MMG3D_newElt(mesh);                                           \
    if ( !jel ) {law;}                                                  \
  }while(0)

/* numerical accuracy */
#define MMG3D_ALPHAD    20.7846096908265 /* 12*sqrt(3) */
#define MMG3D_LLONG     2.5
#define MMG3D_LSHRT     0.3
#define MMG3D_LOPTL     1.3
#define MMG3D_LOPTS     0.6

#define MMG3D_SWAP06       0.0288675 /* 0.6/MMG3D_ALPHAD */
#define MMG3D_SSWAPIMPROVE 1.053
#define MMG3D_LSWAPIMPROVE 1.1

#define MMG3D_BADKAL    0.2
#define MMG3D_MAXKAL     1.


#define MMG3D_NPMAX  1000000 //200000
#define MMG3D_NAMAX   200000 //40000
#define MMG3D_NTMAX  2000000 //400000
#define MMG3D_NEMAX  6000000 //1200000

#define MMG3D_SHORTMAX     0x7fff


/** \brief next vertex of tetra: {1,2,3,0,1,2,3} */
static const unsigned char MMG5_inxt3[7] = { 1,2,3,0,1,2,3 };
/** \brief previous vertex of tetra: {3,0,1,2,3,0,1} */
static const unsigned char MMG5_iprv3[7] = { 3,0,1,2,3,0,1 };
/** \brief idir[i]: vertices of face opposite to vertex i */
static const unsigned char MMG5_idir[4][3] = { {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1} };
/* \var idirinv[i][j]: num of the jth point in the ith face */
static const          char MMG5_idirinv[4][4] = {{-1,0,1,2},{0,-1,2,1},{0,1,-1,2},{0,2,1,-1}};
/** \brief iarf[i]: edges of face opposite to vertex i */
static const unsigned char MMG5_iarf[4][3] = { {5,4,3}, {5,1,2}, {4,2,0}, {3,0,1} };
/** \brief num of the j^th edge in the i^th face */
static const unsigned char MMG5_iarfinv[4][6] = { {-1,-1,-1,2,1,0}, {-1,1,2,-1,-1,0},{2,-1,1,-1,0,-1},{1,2,-1,0,-1,-1}};
/** \brief vertices of extremities of the edges of the tetra */
static const unsigned char MMG5_iare[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} };
/** \brief ifar[i][]: faces sharing the ith edge of the tetra */
static const unsigned char MMG5_ifar[6][2] = { {2,3}, {1,3}, {1,2}, {0,3}, {0,2}, {0,1} };
/** \brief isar[i][]: vertices of extremities of the edge opposite to the ith edge */
static const unsigned char MMG5_isar[6][2] = { {2,3}, {3,1}, {1,2}, {0,3}, {2,0}, {0,1} };
/** \brief arpt[i]: edges passing through vertex i */
static const unsigned char MMG5_arpt[4][3] = { {0,1,2}, {0,4,3}, {1,3,5}, {2,5,4} };

/** \brief idir[i]: vertices of face i for a prism */
static const unsigned char MMG5_idir_pr[5][4] = { {0,1,2,0},{3,5,4,3},{1,4,5,2},{0,2,5,3},{0,3,4,1} };
/** \brief iarf[i]: edges of face i for a prism */
static const unsigned char MMG5_iarf_pr[5][5] = { {0,1,3,0}, {6,8,7,6}, {3,5,8,4}, {5,1,2,7},{0,4,6,2} };

/** Table that associates to each (even) permutation of the 4 vertices of a tetrahedron
 *  the corresponding permutation of its edges.\n Labels :
 *    0  : [0,1,2,3]
 *    1  : [0,2,3,1]
 *    2  : [0,3,1,2]
 *    3  : [1,0,3,2]
 *    4  : [1,2,0,3]
 *    5  : [1,3,2,0]
 *    6  : [2,0,1,3]
 *    7  : [2,1,3,0]
 *    8  : [2,3,0,1]
 *    9  : [3,0,2,1]
 *    10 : [3,1,0,2]
 *    11 : [3,2,1,0]
 *  The edge 0 of the config 1 become the edge 1 of the reference config so permedge[1][0]=1 ...
 */
static const unsigned char MMG5_permedge[12][6] = {
  {0,1,2,3,4,5}, {1,2,0,5,3,4}, {2,0,1,4,5,3}, {0,4,3,2,1,5},
  {3,0,4,1,5,2}, {4,3,0,5,2,1}, {1,3,5,0,2,4}, {3,5,1,4,0,2},
  {5,1,3,2,4,0}, {2,5,4,1,0,3}, {4,2,5,0,3,1}, {5,4,2,3,1,0} };

/**
 * PROctree cell: cellule for point region octree (to speed-up the research of
 * the closest point to another one).
 *
 */
typedef struct MMG3D_PROctree_s
{
  struct MMG3D_PROctree_s* branches; /*!< pointer toward the subtrees of the current PROctree */
  int* v;      /*!< vertex index */
  int  nbVer;  /*!< number of vertices in the sub tree */
  int  depth; /*!< sub tree depth */
} MMG3D_PROctree_s;

/**
 * PROctree global structure (enriched by global variables) for point region
 * octree (to speed-up the research of the closest point to another one).
 */
typedef struct
{
  int nv;  /*!< Max number of points per PROctree cell */
  int nc; /*!< Max number of cells listed per local search in the PROctree (-3)*/
  MMG3D_PROctree_s* q0; /*!<  Pointer toward the first PROctree cell */
} MMG3D_PROctree;
typedef MMG3D_PROctree * MMG3D_pPROctree;


/* PROctree */
void MMG3D_initPROctree_s( MMG3D_PROctree_s* q);
int MMG3D_initPROctree(MMG5_pMesh,MMG3D_pPROctree* q, int nv);
void MMG3D_freePROctree_s(MMG5_pMesh,MMG3D_PROctree_s* q, int nv);
void MMG3D_freePROctree(MMG5_pMesh,MMG3D_PROctree** q);
int MMG3D_isCellIncluded(double* cellCenter, double l, double* zoneCenter, double l0);
void MMG3D_placeInListDouble(double*, double, int, int);
void MMG3D_placeInListPROctree(MMG3D_PROctree_s**, MMG3D_PROctree_s*, int, int);
int MMG3D_seekIndex (double* distList, double dist, int indexMin, int indexMax);
int MMG3D_intersectRect(double *rectin, double *rectinout);
int MMG3D_getListSquareRec(MMG3D_PROctree_s*,double*,double*,
                            MMG3D_PROctree_s***,double*,double*,double, int, int, int*);
int  MMG3D_getListSquare(MMG5_pMesh,double*,MMG3D_PROctree*,double*,MMG3D_PROctree_s***);
int MMG3D_addPROctreeRec(MMG5_pMesh,MMG3D_PROctree_s*,double*, const int, int);
int MMG3D_addPROctree(MMG5_pMesh mesh, MMG3D_PROctree* q, const int no);
int MMG3D_delPROctreeVertex(MMG5_pMesh,MMG3D_PROctree_s* q, int no);
int MMG3D_movePROctree(MMG5_pMesh, MMG3D_pPROctree,int, double*, double*);
void MMG3D_mergeBranchesRec(MMG3D_PROctree_s*, MMG3D_PROctree_s*, int, int , int*);
void MMG3D_mergeBranches(MMG5_pMesh mesh,MMG3D_PROctree_s* q, int dim, int nv);
int MMG3D_delPROctreeRec(MMG5_pMesh,MMG3D_PROctree_s*,double*, const int,const int);
int MMG3D_delPROctree(MMG5_pMesh mesh, MMG3D_pPROctree q, const int no);
void MMG3D_printArbreDepth(MMG3D_PROctree_s* q, int depth, int nv, int dim);
void MMG3D_printArbre(MMG3D_PROctree* q);
void  MMG3D_sizeArbreRec(MMG3D_PROctree_s* q, int nv, int dim, int*,int*);
int*  MMG3D_sizeArbre(MMG3D_PROctree* q, int dim);
int  MMG3D_PROctreein_iso(MMG5_pMesh,MMG5_pSol,MMG3D_pPROctree,int,double);
int  MMG3D_PROctreein_ani(MMG5_pMesh,MMG5_pSol,MMG3D_pPROctree,int,double);

/* prototypes */
int  MMG3D_tetraQual(MMG5_pMesh mesh, MMG5_pSol met,char metRidTyp);
void MMG3D_solTruncature(MMG5_pMesh mesh, MMG5_pSol met);
extern int MMG5_directsurfball(MMG5_pMesh mesh, int ip, int *list, int ilist, double n[3]);

int  MMG3D_Init_mesh_var( va_list argptr );
int  MMG3D_Free_all_var( va_list argptr );
int  MMG3D_Free_structures_var( va_list argptr );
int  MMG3D_Free_names_var( va_list argptr );
int  MMG3D_newPt(MMG5_pMesh mesh,double c[3],int16_t tag);
int  MMG3D_newElt(MMG5_pMesh mesh);
int  MMG3D_delElt(MMG5_pMesh mesh,int iel);
void MMG3D_delPt(MMG5_pMesh mesh,int ip);
int  MMG3D_zaldy(MMG5_pMesh mesh);
void MMG5_freeXTets(MMG5_pMesh mesh);
void MMG5_freeXPrisms(MMG5_pMesh mesh);
void MMG3D_Free_topoTables(MMG5_pMesh mesh);
char MMG5_chkedg(MMG5_pMesh mesh,MMG5_pTria pt,char ori,double,double,int);
int  MMG5_chkBdryTria(MMG5_pMesh mesh);
int  MMG5_mmg3dBezierCP(MMG5_pMesh mesh,MMG5_Tria *pt,MMG5_pBezier pb,char ori);
extern int    MMG5_BezierTgt(double c1[3],double c2[3],double n1[3],double n2[3],double t1[3],double t2[3]);
extern double MMG5_BezierGeod(double c1[3], double c2[3], double t1[3], double t2[3]);
int  MMG3D_bezierInt(MMG5_pBezier pb,double uv[2],double o[3],double no[3],double to[3]);
extern int  MMG5_BezierReg(MMG5_pMesh mesh,int ip0, int ip1, double s, double v[3], double *o, double *no);
extern int  MMG5_BezierRef(MMG5_pMesh mesh,int ip0, int ip1, double s, double *o, double *no, double *to);
extern int  MMG5_BezierEdge(MMG5_pMesh mesh,int ip0, int ip1, double b0[3], double b1[3],char isrid, double v[3]);
extern int  MMG5_BezierRidge(MMG5_pMesh mesh,int ip0, int ip1, double s, double *o, double *no1, double *no2, double *to);
extern int  MMG5_BezierNom(MMG5_pMesh mesh,int ip0,int ip1,double s,double *o,double *no,double *to);
int  MMG5_norface(MMG5_pMesh mesh ,int k, int iface, double v[3]);
int  MMG3D_findEdge(MMG5_pMesh,MMG5_pTetra,int,int,int,int,char*,char* );
int  MMG5_boulernm (MMG5_pMesh mesh, int start, int ip, int *ng, int *nr);
int  MMG5_boulenm(MMG5_pMesh mesh, int start, int ip, int iface, double n[3],double t[3]);
int  MMG5_boulevolp(MMG5_pMesh mesh, int start, int ip, int * list);
int  MMG5_boulesurfvolp(MMG5_pMesh mesh,int start,int ip,int iface,int *listv,
                         int *ilistv,int *lists,int*ilists, int isnm);
int  MMG5_bouletrid(MMG5_pMesh,int,int,int,int *,int *,int *,int *,int *,int *);
int  MMG5_startedgsurfball(MMG5_pMesh mesh,int nump,int numq,int *list,int ilist);
int  MMG5_srcbdy(MMG5_pMesh mesh,int start,int ia);
int  MMG5_coquil(MMG5_pMesh mesh, int start, int ia, int * list);
int  MMG5_coquilface(MMG5_pMesh mesh, int start,char iface,int,int*,int*,int*,int);
int MMG3D_coquilFaceFirstLoop(MMG5_pMesh mesh,int start,int na,int nb,char iface,
                               char ia,int *list,int *ilist,int *it1,int *it2,
                               int *piv,int *adj,char *hasadja,int *nbdy,int silent);
void MMG3D_coquilFaceSecondLoopInit(MMG5_pMesh mesh,int piv,char *iface,char *i,
                                     int *list,int *ilist,int *it1,int *pradj,
                                     int *adj);
void MMG5_coquilFaceErrorMessage(MMG5_pMesh mesh, int k1, int k2);
int16_t MMG5_coquilTravel(MMG5_pMesh,int,int,int*,int*,char*,char*);
int16_t MMG5_openCoquilTravel(MMG5_pMesh,int,int,int*,int*,char*,char*);
int  MMG5_settag(MMG5_pMesh,int,int,int16_t,int);
int  MMG5_deltag(MMG5_pMesh,int,int,int16_t);
int  MMG5_setNmTag(MMG5_pMesh mesh, MMG5_Hash *hash);
int  MMG5_chkcol_int(MMG5_pMesh,MMG5_pSol,int,char,char,int*,int,char);
int  MMG5_chkcol_bdy(MMG5_pMesh,MMG5_pSol,int,char,char,int*,int,int*,int,char);
int  MMG5_chkmanicoll(MMG5_pMesh,int,int,int,int,int,char,char);
int  MMG5_chkmani(MMG5_pMesh mesh);
int  MMG5_colver(MMG5_pMesh,MMG5_pSol,int *,int,char,char);
int  MMG3D_analys(MMG5_pMesh mesh);
void MMG3D_set_reqBoundaries(MMG5_pMesh mesh);
int  MMG5_chkVertexConnectedDomains(MMG5_pMesh mesh);
int  MMG5_norver(MMG5_pMesh mesh);
int  MMG5_setadj(MMG5_pMesh mesh);
int  MMG5_setdhd(MMG5_pMesh mesh);
int  MMG5_singul(MMG5_pMesh mesh);
int  MMG3D_nmgeom(MMG5_pMesh mesh);
int  MMG5_paktet(MMG5_pMesh mesh);
int  MMG5_hashFace(MMG5_pMesh,MMG5_Hash*,int,int,int,int);
int  MMG5_hashGetFace(MMG5_Hash*,int,int,int);
int  MMG3D_hashTria(MMG5_pMesh mesh, MMG5_Hash*);
int   MMG3D_hashPrism(MMG5_pMesh mesh);
int  MMG5_hashPop(MMG5_Hash *hash,int a,int b);
int  MMG5_hPop(MMG5_HGeom *hash,int a,int b,int *ref,int16_t *tag);
int  MMG5_hTag(MMG5_HGeom *hash,int a,int b,int ref,int16_t tag);
int  MMG5_hGet(MMG5_HGeom *hash,int a,int b,int *ref,int16_t *tag);
int  MMG5_hEdge(MMG5_pMesh mesh,MMG5_HGeom *hash,int a,int b,int ref,int16_t tag);
int  MMG5_hNew(MMG5_pMesh mesh,MMG5_HGeom *hash,int hsiz,int hmax);
int  MMG5_hGeom(MMG5_pMesh mesh);
int  MMG5_bdryIso(MMG5_pMesh );
int  MMG5_bdrySet(MMG5_pMesh );
int  MMG5_bdryUpdate(MMG5_pMesh );
int  MMG5_bdryPerm(MMG5_pMesh );
int  MMG5_chkfemtopo(MMG5_pMesh mesh);
int  MMG5_cntbdypt(MMG5_pMesh mesh, int nump);
size_t MMG5_memSize(void);
int  MMG3D_memOption(MMG5_pMesh mesh);
int  MMG3D_memOption_memSet(MMG5_pMesh mesh);
int  MMG3D_memOption_memRepartition(MMG5_pMesh mesh);
int  MMG5_mmg3d1_pattern(MMG5_pMesh ,MMG5_pSol );
int  MMG5_mmg3d1_delone(MMG5_pMesh ,MMG5_pSol );
int  MMG3D_mmg3d2(MMG5_pMesh ,MMG5_pSol );
int  MMG5_mmg3dChkmsh(MMG5_pMesh,int,int);
int   MMG3D_setMeshSize_initData(MMG5_pMesh,int,int,int,int,int,int);
int   MMG3D_setMeshSize_alloc(MMG5_pMesh);
int  MMG3D_split1_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split1(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char metRidTyp);
int  MMG5_split1b(MMG5_pMesh,MMG5_pSol,int*,int,int,int,char,char);
int  MMG5_splitedg(MMG5_pMesh mesh, MMG5_pSol met,int iel, int iar, double crit);
int  MMG3D_split2sf_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split2sf(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG3D_split2_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split2(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG3D_split3_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split3(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG3D_split3cone_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split3cone(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG3D_split3op_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split3op(MMG5_pMesh mesh, MMG5_pSol met, int k, int vx[6],char);
int  MMG3D_split4sf_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split4sf(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG3D_split4op_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split4op(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG3D_split5_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split5(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG3D_split6_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6]);
int  MMG5_split6(MMG5_pMesh mesh,MMG5_pSol met,int k,int vx[6],char);
int  MMG5_split4bar(MMG5_pMesh mesh,MMG5_pSol met,int k,char);
int  MMG3D_simbulgept(MMG5_pMesh mesh,MMG5_pSol met, int *list, int ilist,int);
void MMG5_nsort(int ,double *,char *);
int    MMG3D_optlap(MMG5_pMesh ,MMG5_pSol );
int    MMG5_movintpt_iso(MMG5_pMesh ,MMG5_pSol,MMG3D_pPROctree, int *, int , int);
int    MMG3D_movnormal_iso(MMG5_pMesh ,MMG5_pSol ,int ,int );
int    MMG5_movintptLES_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree,int *,int,int);
int    MMG5_movintpt_ani(MMG5_pMesh ,MMG5_pSol,MMG3D_pPROctree,int *,int ,int);
int    MMG5_movbdyregpt_iso(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree,
                             int*, int, int*, int, int ,int);
int    MMG5_movbdyregpt_ani(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree,
                             int*, int, int*, int, int ,int);
int    MMG5_movbdyrefpt_iso(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree, int*, int,
                             int*, int ,int);
int    MMG5_movbdyrefpt_ani(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree, int*, int,
                             int*, int ,int);
int    MMG5_movbdynompt_iso(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree, int*, int,
                             int*, int ,int);
int    MMG5_movbdynompt_ani(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree, int*, int,
                             int*, int ,int);
int    MMG5_movbdyridpt_iso(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree, int*, int,
                             int*, int ,int);
int    MMG5_movbdyridpt_ani(MMG5_pMesh, MMG5_pSol,MMG3D_pPROctree, int*, int,
                             int*, int ,int);
int    MMG3D_movv_ani(MMG5_pMesh ,MMG5_pSol ,int ,int );
int    MMG3D_movv_iso(MMG5_pMesh ,MMG5_pSol ,int ,int );
int  MMG3D_normalAdjaTri(MMG5_pMesh,int,char,int,double n[3]);
int  MMG5_chkswpbdy(MMG5_pMesh, MMG5_pSol,int*, int, int, int,char);
int  MMG5_swpbdy(MMG5_pMesh,MMG5_pSol,int*,int,int,MMG3D_pPROctree,char);
int  MMG5_swpgen(MMG5_pMesh,MMG5_pSol,int, int, int*,MMG3D_pPROctree,char);
int  MMG5_chkswpgen(MMG5_pMesh,MMG5_pSol,int,int,int*,int*,double,char);
int  MMG3D_swap23(MMG5_pMesh,MMG5_pSol,int,char,int,int,int,int);
int  MMG5_srcface(MMG5_pMesh mesh,int n0,int n1,int n2);
int MMG5_chkptonbdy(MMG5_pMesh,int);
double MMG5_orcal_poi(double a[3],double b[3],double c[3],double d[3]);
int MMG5_countelt(MMG5_pMesh mesh,MMG5_pSol sol, double *weightelt, long *npcible);
/*function for agressive optimization*/
  int MMG3D_opttyp(MMG5_pMesh , MMG5_pSol ,MMG3D_pPROctree ,int);
int MMG3D_swpItem(MMG5_pMesh ,  MMG5_pSol ,MMG3D_pPROctree ,int ,int );
int MMG3D_splitItem(MMG5_pMesh ,  MMG5_pSol ,MMG3D_pPROctree ,int ,int ,double );
int MMG3D_optbdry(MMG5_pMesh ,MMG5_pSol ,MMG3D_pPROctree ,int );
int MMG3D_movetetrapoints(MMG5_pMesh ,MMG5_pSol ,MMG3D_pPROctree ,int ) ;

int MMG5_trydisp(MMG5_pMesh,double *,short);
int MMG5_dichodisp(MMG5_pMesh,double *);
int MMG5_lapantilap(MMG5_pMesh,double *);
int MMG5_ppgdisp(MMG5_pMesh,double *);
int MMG5_denoisbdy(MMG5_pMesh);
int MMG3D_displayQualHisto(int,double,double,double,int,int,int,
                           int his[5],int,int,int);
int MMG3D_displayQualHisto_internal(int,double,double,double,int,int,int,
                                    int his[5],int,int,int);
void MMG3D_computeInqua(MMG5_pMesh,MMG5_pSol,int*,double*,double*,double*,int*,int*,
                        int*,int his[5],int);
int  MMG3D_inqua(MMG5_pMesh mesh,MMG5_pSol met);
void MMG3D_computeOutqua(MMG5_pMesh,MMG5_pSol,int*,double*,double*,double*,int*,int*,
                         int*,int his[5],int*,int);
int  MMG3D_outqua(MMG5_pMesh mesh,MMG5_pSol met);
void MMG3D_computeLESqua(MMG5_pMesh,MMG5_pSol,int*,double*,double*,double*,int*,int*,
                         int*,int his[5],int);
int MMG3D_computePrilen(MMG5_pMesh,MMG5_pSol,double*,double*,double*,int*,int*,int*,
                        int*,int*,int*,char,double**, int [9] );
int  MMG3D_prilen(MMG5_pMesh mesh,MMG5_pSol met,char);
void MMG3D_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met);
void MMG5_defaultValues(MMG5_pMesh);
int  MMG5_intridmet(MMG5_pMesh,MMG5_pSol,int,int,double,double*,double*);
int  MMG5_intregmet(MMG5_pMesh,MMG5_pSol,int,char,double, double*);
int  MMG5_intvolmet(MMG5_pMesh,MMG5_pSol,int,char,double, double*);
int MMG3D_localParamReg(MMG5_pMesh,int,int*,int,int*,int,double*,double*,double*);
int MMG3D_localParamNm(MMG5_pMesh,int,int,int,double*,double*,double*);
int  MMG3D_mark_packedPoints(MMG5_pMesh mesh,int *np,int *nc);
int  MMG3D_pack_tetraAndAdja(MMG5_pMesh mesh);
int  MMG3D_pack_tetra(MMG5_pMesh mesh);
int  MMG3D_pack_prismsAndQuads(MMG5_pMesh mesh);
int  MMG3D_pack_sol(MMG5_pMesh mesh,MMG5_pSol sol);
int  MMG3D_update_eltsVertices(MMG5_pMesh mesh);
int  MMG3D_pack_pointArray(MMG5_pMesh mesh);
int  MMG3D_pack_points(MMG5_pMesh mesh);
void MMG3D_unset_reqBoundaries(MMG5_pMesh mesh);
int  MMG3D_packMesh(MMG5_pMesh,MMG5_pSol,MMG5_pSol);
int  MMG3D_bdryBuild(MMG5_pMesh);

/* useful functions to debug */
int  MMG3D_indElt(MMG5_pMesh mesh,int kel);
int  MMG3D_indPt(MMG5_pMesh mesh,int kp);
void MMG5_printTetra(MMG5_pMesh mesh,char* fileName);


#ifdef USE_SCOTCH
int MMG5_mmg3dRenumbering(int vertBoxNbr, MMG5_pMesh mesh, MMG5_pSol sol);
#endif

int    MMG5_meancur(MMG5_pMesh mesh,int np,double c[3],int ilist,int *list,double h[3]);
double MMG5_surftri(MMG5_pMesh,int,int);
double MMG5_timestepMCF(MMG5_pMesh,double);
int    MMG5_bdyMCF(MMG5_pMesh);
double MMG5_volint(MMG5_pMesh);

/* Lagrangian mode functions */
double MMG5_estavglen(MMG5_pMesh);
int   MMG5_stiffelt(MMG5_pMesh,int,double*,double*);
int  MMG5_mmg3d3(MMG5_pMesh ,MMG5_pSol, MMG5_pSol );
int  MMG5_velextLS(MMG5_pMesh ,MMG5_pSol );
int MMG5_saveDisp(MMG5_pMesh ,MMG5_pSol );

/* Delaunay functions*/
int MMG5_delone(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int *list,int ilist);
  int MMG5_cavity_iso(MMG5_pMesh mesh,MMG5_pSol sol,int iel,int ip,int *list,int lon,double volmin);
  int MMG5_cavity_ani(MMG5_pMesh mesh,MMG5_pSol sol,int iel,int ip,int *list,int lon,double volmin);
int MMG5_cenrad_iso(MMG5_pMesh mesh,double *ct,double *c,double *rad);
int MMG5_cenrad_ani(MMG5_pMesh mesh,double *ct,double *m,double *c,double *rad);

/* mmg3d1.c */
void MMG5_tet2tri(MMG5_pMesh mesh,int k,char ie,MMG5_Tria *ptt);
int  MMG3D_dichoto(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx);
int  MMG3D_dichoto1b(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ret,int);
char MMG5_chkedg(MMG5_pMesh mesh,MMG5_Tria *pt,char ori,double,double,int);
int  MMG3D_splsurfedge( MMG5_pMesh,MMG5_pSol,int,MMG5_pTetra,MMG5_pxTetra,char,
                        char,char,int* );
int  MMG5_anatet(MMG5_pMesh mesh,MMG5_pSol met, char typchk, int patternMode) ;
int  MMG5_movtet(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree PROctree,
                  double clickSurf,double clickVol,int moveVol,int improveSurf,int improveVolSurf,
                  int improveVol,int maxit,int testmark);
int  MMG5_swpmsh(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree PROctree, int);
  int  MMG5_swptet(MMG5_pMesh mesh,MMG5_pSol met,double,double,MMG3D_pPROctree, int,int);

/* pointers */
/* init structures */
void  MMG5_Init_parameters(MMG5_pMesh mesh);
/* iso/aniso computations */
double MMG5_caltet33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt);
extern double MMG5_lenedgCoor_iso(double*, double*, double*, double*);
int    MMG5_intmet_iso(MMG5_pMesh,MMG5_pSol,int,char,int, double);
int    MMG5_intmet_ani(MMG5_pMesh,MMG5_pSol,int,char,int, double);
int    MMG3D_intmet33_ani(MMG5_pMesh,MMG5_pSol,int,char,int, double);
int    MMG5_interp4bar_ani(MMG5_pMesh,MMG5_pSol,int,int,double *);
int    MMG5_interp4bar33_ani(MMG5_pMesh,MMG5_pSol,int,int,double *);
int    MMG5_interp4bar_iso(MMG5_pMesh,MMG5_pSol,int,int,double *);
int    MMG3D_defsiz_iso(MMG5_pMesh,MMG5_pSol );
int    MMG3D_defsiz_ani(MMG5_pMesh ,MMG5_pSol );
int    MMG3D_gradsiz_iso(MMG5_pMesh ,MMG5_pSol );
int    MMG3D_gradsiz_ani(MMG5_pMesh ,MMG5_pSol );
int    MMG3D_gradsizreq_iso(MMG5_pMesh ,MMG5_pSol );
int    MMG3D_gradsizreq_ani(MMG5_pMesh ,MMG5_pSol );
double     MMG5_meansizreg_iso(MMG5_pMesh,MMG5_pSol,int,int*,int,double,double);
int        MMG3D_chk4ridVertices(MMG5_pMesh mesh,MMG5_pTetra pt);
extern int MMG5_moymet(MMG5_pMesh ,MMG5_pSol ,MMG5_pTetra ,double *);
int    MMG3D_set_metricAtPointsOnReqEdges (MMG5_pMesh,MMG5_pSol);
void MMG3D_mark_pointsOnReqEdge_fromTetra (  MMG5_pMesh mesh );

extern double (*MMG5_lenedg)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra );
extern double (*MMG5_lenedgspl)(MMG5_pMesh ,MMG5_pSol ,int, MMG5_pTetra );
extern double (*MMG5_caltet)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt);
extern double (*MMG5_caltri)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt);
extern int    (*MMG3D_defsiz)(MMG5_pMesh ,MMG5_pSol );
extern int    (*MMG3D_gradsiz)(MMG5_pMesh ,MMG5_pSol );
extern int    (*MMG3D_gradsizreq)(MMG5_pMesh ,MMG5_pSol );
extern int    (*MMG5_intmet)(MMG5_pMesh,MMG5_pSol,int,char,int, double);
extern int    (*MMG5_interp4bar)(MMG5_pMesh,MMG5_pSol,int,int,double *);
extern int    (*MMG5_movintpt)(MMG5_pMesh ,MMG5_pSol, MMG3D_pPROctree ,int *, int , int );
extern int    (*MMG5_movbdyregpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int, int ,int);
extern int    (*MMG5_movbdyrefpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int ,int);
extern int    (*MMG5_movbdynompt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int ,int);
extern int    (*MMG5_movbdyridpt)(MMG5_pMesh, MMG5_pSol, MMG3D_pPROctree ,int*, int, int*, int ,int);
extern int    (*MMG5_cavity)(MMG5_pMesh ,MMG5_pSol ,int ,int ,int *,int ,double);
extern int    (*MMG3D_PROctreein)(MMG5_pMesh ,MMG5_pSol ,MMG3D_pPROctree ,int,double );

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Warn user that some tetrahedra of the mesh have been reoriented.
 *
 */
static inline
void MMG5_warnOrientation(MMG5_pMesh mesh) {
  if ( mesh->xt ) {
    if ( mesh->xt != mesh->ne ) {
      fprintf(stderr,"\n  ## Warning: %s: %d tetra on %d reoriented.\n",
              __func__,mesh->xt,mesh->ne);
      fprintf(stderr,"  Your mesh may be non-conform.\n");
    }
    else {
      fprintf(stderr,"\n  ## Warning: %s: all tetra reoriented.\n",__func__);
    }
  }
  mesh->xt = 0;
}

/**
 * Set common pointer functions between mmgs and mmg3d to the matching mmg3d
 * functions.
 */
static inline
void MMG3D_Set_commonFunc() {
  MMG5_bezierCP          = MMG5_mmg3dBezierCP;
  MMG5_chkmsh            = MMG5_mmg3dChkmsh;
  MMG5_indPt             = MMG3D_indPt;
  MMG5_indElt            = MMG3D_indElt;
  MMG5_grad2met_ani      = MMG5_grad2metSurf;
  MMG5_grad2metreq_ani   = MMG5_grad2metSurfreq;
#ifdef USE_SCOTCH
  MMG5_renumbering       = MMG5_mmg3dRenumbering;
#endif
}

#ifdef __cplusplus
}
#endif

#endif
