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

#ifndef PROCTREE_3D_H
#define PROCTREE_3D_H

#include "libmmgtypes.h"

/**
 * PROctree cell: cellule for point region octree (to speed-up the research of
 * the closest point to another one).
 *
 */
typedef struct MMG3D_PROctree_s
{
  struct MMG3D_PROctree_s* branches; /*!< pointer to the subtrees of the current PROctree */
  MMG5_int* v;      /*!< vertex index */
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
int MMG3D_addPROctreeRec(MMG5_pMesh,MMG3D_PROctree_s*,double*, const MMG5_int, int);
int MMG3D_addPROctree(MMG5_pMesh mesh, MMG3D_PROctree* q, const MMG5_int no);
int MMG3D_delPROctreeVertex(MMG5_pMesh,MMG3D_PROctree_s* q, MMG5_int no);
int MMG3D_movePROctree(MMG5_pMesh, MMG3D_pPROctree,MMG5_int, double*, double*);
void MMG3D_mergeBranchesRec(MMG3D_PROctree_s*, MMG3D_PROctree_s*, int, int , int*);
void MMG3D_mergeBranches(MMG5_pMesh mesh,MMG3D_PROctree_s* q, int dim, int nv);
int MMG3D_delPROctreeRec(MMG5_pMesh,MMG3D_PROctree_s*,double*, const MMG5_int,const int);
int MMG3D_delPROctree(MMG5_pMesh mesh, MMG3D_pPROctree q, const int no);
void MMG3D_printArbreDepth(MMG3D_PROctree_s* q, int depth, int nv, int dim);
void MMG3D_printArbre(MMG3D_PROctree* q);
void  MMG3D_sizeArbreRec(MMG3D_PROctree_s* q, int nv, int dim, int*,int*);
int*  MMG3D_sizeArbre(MMG3D_PROctree* q, int dim);
int  MMG3D_PROctreein_iso(MMG5_pMesh,MMG5_pSol,MMG3D_pPROctree,MMG5_int,double);
int  MMG3D_PROctreein_ani(MMG5_pMesh,MMG5_pSol,MMG3D_pPROctree,MMG5_int,double);
int64_t MMG3D_getPROctreeCoordinate(MMG3D_pPROctree q, double* ver, int dim);

#endif
