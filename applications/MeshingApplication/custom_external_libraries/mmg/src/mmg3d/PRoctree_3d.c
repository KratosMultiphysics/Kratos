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

/**
 * \file PROctree_3d.c
 * \brief Tools for local search around coordinates based on PROctree.
 * \author Jean Mercat (Inria/UBordeaux)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * An PROctree of the nodes is created and used for local neighbor search.
 * This helps deciding if a position is too close to other nodes to refine
 * with an insertion of a new node.
 *
 */

#include "libmmgtypes.h"
#include "mmgcommon_private.h"
#include "PRoctree_3d_private.h"
#include <stdio.h>

/**
 * \param q pointer to the PROctree cell
 *
 * Initialisation of the PROctree cell.
 *
 */
void MMG3D_initPROctree_s( MMG3D_PROctree_s* q)
{
  q->nbVer = 0;
  q->depth = 0;
  q->v = NULL;
  q->branches = NULL;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to the global PROctree
 * \param nv maximum number of vertices in each cell before subdivision
 * \return 1 if ok 0 if memory saturated
 *
 * Initialisation of the PROctree cell.
 *
 */
int MMG3D_initPROctree(MMG5_pMesh mesh,MMG3D_pPROctree* q, int nv)
{
  MMG5_int i;

  MMG5_ADD_MEM(mesh,sizeof(MMG3D_PROctree),"PROctree structure",
                return 0);
  MMG5_SAFE_MALLOC(*q,1, MMG3D_PROctree, return 0);


  // set nv to the next power of 2
  nv--;
  nv |= nv >> 1;
  nv |= nv >> 2;
  nv |= nv >> 4;
  nv |= nv >> 8;
  nv |= nv >> 16;
  nv++;
  (*q)->nv = nv;

  // Number maximum of cells listed for the zone search
  (*q)->nc = MG_MAX(2048/nv,16);

  MMG5_ADD_MEM(mesh,sizeof(MMG3D_PROctree_s),"initial PROctree cell",
                return 0);

  MMG5_SAFE_MALLOC((*q)->q0,1, MMG3D_PROctree_s, return 0);
  MMG3D_initPROctree_s((*q)->q0);

  for (i=1;i<=mesh->np; ++i)
  {
    if ( !MG_VOK(&mesh->point[i]) )  continue;
    if (mesh->point[i].tag & MG_BDY) continue;

    if(!MMG3D_addPROctree(mesh, (*q), i))
      return 0;

  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to the PROctree cell
 * \param nv number of vertices in the cell subtree
 *
 * Free the PROctree cell.
 *
 */
void MMG3D_freePROctree_s(MMG5_pMesh mesh,MMG3D_PROctree_s* q, int nv)
{
  int nbBitsInt,depthMax,dim,i,sizBr,nvTemp;

  dim       = mesh->dim;
  sizBr     = 1<<dim;
  nbBitsInt = sizeof(int64_t)*8;
  depthMax  = nbBitsInt/dim - 1;

  if (q->nbVer>nv && q->depth < depthMax )
  {
    for ( i = 0; i<sizBr; i++)
    {
      MMG3D_freePROctree_s(mesh,&(q->branches[i]), nv);
    }
    MMG5_DEL_MEM(mesh,q->branches);
    q->branches = NULL;
  }
  else if (q->nbVer>0)
  {
    if ( q->nbVer<= nv )
    {
      nvTemp = q->nbVer;
      nvTemp--;
      nvTemp |= nvTemp >> 1;
      nvTemp |= nvTemp >> 2;
      nvTemp |= nvTemp >> 4;
      nvTemp |= nvTemp >> 8;
      nvTemp |= nvTemp >> 16;
      nvTemp++;

      MMG5_DEL_MEM(mesh,q->v);
      q->v = NULL;
      q->nbVer = 0;
    }else
    {
      assert(q->v);
      MMG5_DEL_MEM(mesh,q->v);
      q->v = NULL;
      q->nbVer = 0;
    }
  }
}

/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to a pointer to the global PROctree.
 *
 * Free the global PROctree structure.
 *
 */
void MMG3D_freePROctree(MMG5_pMesh mesh,MMG3D_pPROctree *q)
{
  MMG3D_freePROctree_s(mesh,(*q)->q0, (*q)->nv);
  MMG5_DEL_MEM(mesh,(*q)->q0);
  (*q)->q0 = NULL;
  MMG5_DEL_MEM(mesh,*q);
  *q = NULL;
}


/**
 * \param q pointer to the global PROctree.
 * \param ver coordinates of the point.
 * \param dim space dimension (should be 3).
 * \return the integer containing the coordinates
 *
 * Get the integer containing the coordinates
 *
 */
int64_t MMG3D_getPROctreeCoordinate(MMG3D_pPROctree q, double* ver, int dim)
{
  int64_t s    = 1<<20;
  double  prec = 1./(1<<30);
  int place = 0;
  int ix = (int)floor((ver[0]-prec)*s);
  int iy = (int)floor((ver[1]-prec)*s);
  int iz = (int)floor((ver[2]-prec)*s);
  ix = (ix > 0) ? ix:0;
  iy = (iy > 0) ? iy:0;
  iz = (iz > 0) ? iz:0;
  int64_t i=0;
  int j;
  for(j=19; j>=0; j--)
  {
    s=s>>1;
    i += ((ix & s) >> j)<<place;
    place++;
    i += ((iy & s) >> j)<<place;
    place++;
    i += ((iz & s) >> j)<<place;
    place++;
  }
  return i;
}



/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to the global PROctree.
 * \param no index of the moved point.
 * \param newVer new coordinates for the moved point.
 * \param oldVer old coordinates for the moved point.
 * \return 1 if ok 0 if memory saturated
 *
 * Move one point in the PROctree structure. /!\ the vertex of index \a no
 * can have either the new or the old coordinates in the mesh but all
 * other vertice should have the same coordinates as when they were inserted
 * into the PROctree. (ie: one move at a time in the mesh and the PROctree)
 *
 */
int MMG3D_movePROctree(MMG5_pMesh mesh, MMG3D_pPROctree q, MMG5_int no, double* newVer, double* oldVer)
{
  int64_t oldCoor, newCoor;
  double  pt[3];
  int     dim;

  dim = mesh->dim;

  memcpy(&pt, oldVer ,dim*sizeof(double));
  oldCoor = MMG3D_getPROctreeCoordinate(q, oldVer, dim);
  memcpy(&pt, newVer ,dim*sizeof(double));
  newCoor = MMG3D_getPROctreeCoordinate(q, mesh->point[no].c, dim);

  if (newCoor == oldCoor) {
    return 1;
  }
  else // it could be possible to combine delPROctree and addPROctree to keep it local...
  {
    /* delPROctree */
    memcpy(&pt, oldVer ,dim*sizeof(double));
    if (!MMG3D_delPROctreeRec(mesh, q->q0, pt , no, q->nv))
      return 0;

    /* addPROctree */
    memcpy(&pt, newVer ,dim*sizeof(double));
    if(!MMG3D_addPROctreeRec(mesh, q->q0, pt , no, q->nv))
      return 0;
  }
  return 1;
}

/**
 *
 * \param cellCenter 3 coordinates of the center of the PROctree cell to test.
 * \param l size of the cell
 * \param zoneCenter 3 coordinates of the center of the search zone
 * \param radius of the search zone
 * \return wether the cell is included in the search zone.
 *
 */
int MMG3D_isCellIncluded(double* cellCenter, double l, double* zoneCenter, double l0)
{
  double x,y,z,r1;//r2,rmax;

  x = cellCenter[0]-zoneCenter[0];
  y = cellCenter[1]-zoneCenter[1];
  z = cellCenter[2]-zoneCenter[2];

  r1 = sqrt(x*x+y*y+z*z);

  // to avoid sqrt :
  //~ r1 = x*x+y*y+z*z;
  //~ r2 = 3*l*l;
  //~ rmax = MG_MAX(r1,r2);
  //sqrt(3)=1.7320508...
  //~ return (r1+r2+2*rmax<l0*l0);

  return (r1+1.732051*l<l0);
}


/**
 * \param distList list of values.
 * \param dist value to insert in the list.
 * \param index position of the element before the place where \a dist
 * should be inserted.
 * \param size size of the list before insertion.
 *
 * Insert the value \a dist in the list \a distList at position \a index+1.
 * Moves other data so nothing is lost. No memory check performed, this
 * function should be called with coherent parameters.
 *
*/
void MMG3D_placeInListDouble(double* distList, double dist, int index, int size)
{
  memmove(&(distList[index+2]),&(distList[index+1]),(size-(index+1))*sizeof(double));
  distList[index+1] = dist;
}

/**
 * \param qList list of pointer on PROctree.
 * \param q pointer on PROctree to be inserted in the list.
 * \param index position of the element before the place where \a q
 * should be inserted.
 * \param size size of the list before insertion.
 *
 * Insert the pointer \a q in the list \a qList at position \a index+1.
 * Moves other data so nothing is lost. No memory check performed, this
 * function should be called with coherent parameters.
 *
 */
void MMG3D_placeInListPROctree(MMG3D_PROctree_s** qlist, MMG3D_PROctree_s* q, int index, int size)
{
  memmove(&(qlist[index+2]),&(qlist[index+1]),(size-(index+1))*sizeof(MMG3D_PROctree_s*));
  #ifdef DEBUG
  if (index+2+(size-(index+1))>61 || index+1<0)
    fprintf(stderr, "\n  ## Error: %s: index"
            " too large %i > 61\n",__func__, index+2+(size-(index+1));
  #endif
  qlist[index+1] = q;
}

/**
 * \param distList ordered list of value from smallest to largest.
 * \param dist value to be compared to elements in the list.
 * \param indexMin minimum index of the list.
 * \param indexMax maximum index of the list.
 *
 * \return Index of the biggest value of disList that is strictly
 * smaller than dist. Only search in the bounds of indexMin and indexMax.
 *
 */
int MMG3D_seekIndex (double* distList, double dist, int indexMin, int indexMax)
{
  int indexMed;

  if (indexMin > indexMax)
    MMG3D_seekIndex(distList, dist, indexMax, indexMin);
  else if (indexMax - indexMin <2)
  {
    #ifdef DEBUG
      if (indexMin >= 60 || indexMax>=60)
        fprintf(stderr,"\n  ## Error: %s: index should not be that large %i %i.\n",
                __func__,indexMin, indexMax);
    #endif
    if (dist > distList[indexMax])
      return indexMax;
    else
      return indexMin;
  }
  else
  {
    indexMed = (indexMin + indexMax)/2;

    #ifdef DEBUG
      if (indexMed >= 60 || indexMed < 0)
        fprintf(stderr,"\n  ## Error: %s: index should not be that large %i.\n",
                __func__,indexMed);
    #endif

    if (dist > distList[indexMed])
      MMG3D_seekIndex(distList, dist, indexMed, indexMax);
    else
      MMG3D_seekIndex(distList, dist, indexMin, indexMed);
  }
  return 1;
}

/**
 * \param rectin rectangle to intersect, is not modified.
 * \param rectinout rectangle to intersect, is set to the intersection.
 *
 * \return 1 if \a rectinout intersect \a rectin, 0 otherwise (possible because
 * the surface reconstruction may leads to point outside the [0;1]x[0;1]x[0;1]
 * bounding box)
 *
 * Set rectinout to the intersection of the two rectangles.
 *  Rectangles are defined by: the coordinates of the lower left corner
 * of the rectange and the length of the rectangle in each dimension.
*/
int MMG3D_intersectRect(double *rectin, double *rectinout)
{
  double rect1Temp[6], rect2Temp[6];

  rect1Temp[0] = rectin[0];
  rect1Temp[1] = rectin[1];
  rect1Temp[2] = rectin[2];
  rect1Temp[3] = rectin[3]+rectin[0];
  rect1Temp[4] = rectin[4]+rectin[1];
  rect1Temp[5] = rectin[5]+rectin[2];

  rect2Temp[0] = rectinout[0];
  rect2Temp[1] = rectinout[1];
  rect2Temp[2] = rectinout[2];
  rect2Temp[3] = rectinout[3]+rectinout[0];
  rect2Temp[4] = rectinout[4]+rectinout[1];
  rect2Temp[5] = rectinout[5]+rectinout[2];

  rectinout[0] = rect1Temp[0]>rect2Temp[0] ? rect1Temp[0]:rect2Temp[0];
  rectinout[1] = rect1Temp[1]>rect2Temp[1] ? rect1Temp[1]:rect2Temp[1];
  rectinout[2] = rect1Temp[2]>rect2Temp[2] ? rect1Temp[2]:rect2Temp[2];
  rectinout[3] = rect1Temp[3]<rect2Temp[3] ? rect1Temp[3]:rect2Temp[3];
  rectinout[4] = rect1Temp[4]<rect2Temp[4] ? rect1Temp[4]:rect2Temp[4];
  rectinout[5] = rect1Temp[5]<rect2Temp[5] ? rect1Temp[5]:rect2Temp[5];

  rectinout[3] = rectinout[3] - rectinout[0];
  rectinout[4] = rectinout[4] - rectinout[1];
  rectinout[5] = rectinout[5] - rectinout[2];

  if ( rectinout[3]<=0 || rectinout[4]<=0 || rectinout[5]<=0 ) return 0;

  return 1;
}

/**
 * \param q pointer to the PROctree cell.
 * \param center coordinates of the centre of the current subtree.
 * \param rect rectangle that we want to intersect with the subtree. We define
 * it given: the coordinates of one corner of the rectange and the length of
 * the rectangle in each dimension.
 * \param qlist pointer to the list of pointer over the sub PROctrees that
 *  intersect \a rect.
 * \param dist pointer to the list of distances between center of
 * the PROctree cells in qlist and the last 3 elements are the coordinates
 * of the center of the whole recangle.
 * \param ani metric of the point.
 * \param l0 radius of the search zone.
 * \param nc number max of cell in the list +3 (the three last.
 * \param dim dimension =3.
 * \param index number of PROctree cells that intersect \a rect
 *
 * \return 0 if the rectangle doesn't intersect the PROctree (possible due to the
 * surface reconstruction), 1 otherwise.
 *
 * List the number of PROctree cells that intersect the rectangle
 * \a rect. To avoid counting of the cells, a maximum is set.
 *
 */
int MMG3D_getListSquareRec(MMG3D_PROctree_s* q, double* center, double* rect,
                            MMG3D_PROctree_s*** qlist, double* dist, double* ani,
                           double l0, int nc, int dim, int* index)
{
  double recttemp[6];
  double centertemp[3];
  int    recCenter[6];
  double l = 1./(1<<(q->depth+1));
  double distTemp;
  double x,y,z;
  int indexTemp,i,j,k,nBranch;

  // number max of PROctree cells listed for one search
  if ((*index)>nc-4)
    return 1;

  // check if the current cell is included in the search zone, can avoid
  // the loop over the vertices. Never occured in tests unless nc==4, in that
  // case, there is no gain in computing time.
  //~ if (q->nbVer>0 && MMG3D_isCellIncluded(center, l, &(dist[nc-3]), l0))
  //~ {
    //~ (*index)=nc-3;
    //~ fprintf(stdout,"Included cell found\n");
    //~ return;
  //~ }

  assert ( nc>= 3 );
  if (q->branches==NULL && q->v != NULL)
  {
    // the vector dist is of size nc whereas qlist allows nc-3 inputs
    // so the 3 last values can contain the coordinates of the center
    // of the search volume.
    x = dist[nc-3] - center[0];
    y = dist[nc-2] - center[1];
    z = dist[nc-1] - center[2];

    // Should be replaced with distance in metric?
    distTemp = x*x+y*y+z*z;

    // Here the anisotropic distance not tested (not so important, this only
    // reorders the cells)
    //~ distTemp = ani[0]*x*x+ani[3]*y*y+ani[5]*z*z+
                //~ 2*(ani[1]*x*y+ani[2]*x*z+ani[4]*y*z);
    if (*index > 0)
    {
      indexTemp = MMG3D_seekIndex(dist,distTemp,0, *index-1);
      if (indexTemp+1<*index)
      {
        MMG3D_placeInListDouble(dist, distTemp, indexTemp, *index);
        MMG3D_placeInListPROctree((*qlist), q, indexTemp, *index);
      }else
      {
        dist[*index]=distTemp;
        (*qlist)[*index]=q;
      }
    }
    else
    {
      dist[*index]=distTemp;
      (*qlist)[*index]=q;
    }

    (*index)++;

  }else if (q->branches!=NULL)
  {

    // check the position of the search zone in the current cell
    for (i=0;i<3;i++)
    {
      recCenter[i] = (rect[i]>center[i]);
      recCenter[i+3] = ((rect[i]+rect[i+3])>center[i]);
    }

    // three loop describing the 8 branches in binary (k,j,i):(0,0,0),(0,0,1)....(1,1,1)
    for(i=0;i<2;i++)
    {
      for(j=0;j<2;j++)
      {
        for(k=0;k<2;k++)
        {
          // test if that branch intersects the rectangle
          if (((i && recCenter[3]) || (!i && !recCenter[0])) &&
              ((j && recCenter[4]) || (!j && !recCenter[1]))&&
              ((k && recCenter[5]) || (!k && !recCenter[2])))
          {
            // set the branch number
            nBranch = i+2*j+4*k;

            // set recttemp to the cell size of the branch nBranch
            recttemp[0] = center[0]-l*(1-i);
            recttemp[1] = center[1]-l*(1-j);
            recttemp[2] = center[2]-l*(1-k);
            recttemp[3] = recttemp[4] = recttemp[5] = l;
            // intersect the rectangle and the cell and store it in recttemp
            if ( !MMG3D_intersectRect(rect,recttemp) ) return 0;

            // set the new center
            centertemp[0] = center[0]-l/2+i*l;
            centertemp[1] = center[1]-l/2+j*l;
            centertemp[2] = center[2]-l/2+k*l;

            // recursive call in the branch
            if ( !MMG3D_getListSquareRec(&(q->branches[nBranch]),
                                          centertemp, recttemp, qlist, dist,
                                          ani, l0, nc, dim, index) )
              return 0;
          }
        }
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param ani metric to use for the cell ordering from closest to farthest
 * \param q pointer to the global PROctree structure.
 * \param rect rectangle that we want to intersect with the subtree. We define
 * it given: the coordinates of one corner of the rectangle and the length of
 * the rectangle in each dimension.
 * \param qlist pointer to the list of pointer over the sub PROctrees that
 *  intersect \a rect.
 *
 * \return index, the number of subtrees in the list
 * \return  -1 if fail due to lack of memory.
 *
 * List the number of PROctree cells that intersect the rectangle \a rect.
 *
 */
int MMG3D_getListSquare(MMG5_pMesh mesh, double* ani, MMG3D_pPROctree q, double* rect,
                         MMG3D_PROctree_s*** qlist)
{
  double rect2[6], center[3], *dist;
  double l0;
  int    i,index,dim;

  dim = mesh->dim;
  index = 0;

  memcpy(&rect2, rect, sizeof(double)*dim*2);

  //instead of counting exactly the number of cells to be listed, the
  //maximum size is set to nc-3 (so the list dist can have nc-3 values + 3 coordinates of
  //the center of the rectangle)
  assert(q->nc>=3);
  index = q->nc-3;

  MMG5_ADD_MEM(mesh,index*sizeof(MMG3D_PROctree_s*),"PROctree cell",return -1);
  MMG5_SAFE_MALLOC(*qlist,index,MMG3D_PROctree_s*, return -1);

  MMG5_ADD_MEM(mesh,q->nc*sizeof(double),"dist array",return -1);
  MMG5_SAFE_MALLOC(dist,q->nc,double,return -1);

  // Set the center of the zone search
  dist[q->nc-3] = rect[0]+rect[3]/2;
  dist[q->nc-2] = rect[1]+rect[4]/2;
  dist[q->nc-1] = rect[2]+rect[5]/2;

  // Set the radius of the zone search (used for cell inclusion test if activated)
  l0 = rect[3]/2;

  // Initialization of the PROctree cell list
  for (i = 0; i<index; i++)
    (*qlist)[i] = NULL;

  index = 0;

  // Set center of the first PROctree cell
  for (i = 0; i < dim; ++i)
    center[i] = 0.5;

  // Avoid modification of input parameter rect
  memcpy(&rect2, rect, sizeof(double)*dim*2);

  if ( !MMG3D_getListSquareRec(q->q0, center, rect2, qlist, dist, ani, l0,
                                q->nc, dim, &index) ) {
    MMG5_DEL_MEM(mesh,dist);
    return 0;
  }


  if (index>q->nc-4)
  {
    MMG5_DEL_MEM(mesh,dist);
    return 0;
  }

  MMG5_DEL_MEM(mesh,dist);

  return index;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to an PROctree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv maximum number of points in an PROctree cell.
 * \return 1 if ok 0 if memory saturated
 *
 * Add vertex in the suitable quadrant of the PROctree. This function is
 * recursively called until we reach the last one. At each step, the vertex
 * coordinates are scaled such as the quadrant is the [0;1]x[0;1]x[0;1] box.
 *
 */
int MMG3D_addPROctreeRec(MMG5_pMesh mesh, MMG3D_PROctree_s* q, double* ver,
                         const MMG5_int no, int nv)
{
  double   pt[3],quadrant;
  int      dim, nbBitsInt,depthMax,i,j,k;
  int      sizBr;
  int      sizeRealloc;

  nbBitsInt = sizeof(int64_t)*8;
  dim       = mesh->dim;
  depthMax  = nbBitsInt/dim - 1; // maximum depth is to allow integer coordinates
  sizBr     = 1<<dim;

  if ( q->depth < depthMax ) // not at the maximum depth of the tree
  {
    if (q->nbVer < nv)  // not at the maximum number of vertice in the cell
    {

      if(q->nbVer == 0)  // first vertex list allocation
      {
        MMG5_ADD_MEM(mesh,sizeof(MMG5_int),"PROctree vertice table", return 0);
        MMG5_SAFE_MALLOC(q->v,1,MMG5_int,return 0);
      }
      else if(!(q->nbVer & (q->nbVer - 1))) //is a power of 2 -> reallocation of the vertex list
      {
        sizeRealloc = q->nbVer;
        sizeRealloc<<=1;
        MMG5_ADD_MEM(mesh,(sizeRealloc-sizeRealloc/2)*sizeof(MMG5_int),"PROctree realloc",
                      return 0);
        MMG5_SAFE_REALLOC(q->v,q->nbVer,sizeRealloc,MMG5_int,"PROctree",return 0);
      }

      q->v[q->nbVer] = no;
      q->nbVer++;
      return 1;
    }
    else if (q->nbVer == nv && q->branches==NULL)  //vertex list at maximum -> cell subdivision
    {
      /* creation of sub-branch and relocation of vertices in the sub-branches */
      MMG5_ADD_MEM(mesh,sizBr*sizeof(MMG3D_PROctree_s),"PROctree branches",
                    return 0);
      MMG5_SAFE_MALLOC(q->branches,sizBr,MMG3D_PROctree_s,return 0);

      for ( i = 0; i<sizBr; i++)
      {
        MMG3D_initPROctree_s(&(q->branches[i]));
        q->branches[i].depth = q->depth+1;
      }
      q->nbVer++;
      for (i = 0; i<nv; i++)
      {

        memcpy(&pt, mesh->point[q->v[i]].c ,dim*sizeof(double));
        for ( j =0; j < q->depth; j++)
        {
          for (k = 0; k<dim; k++)
          {
            pt[k] -= ((double) (pt[k]>0.5))*0.5;
            pt[k] *= 2;
          }
        }
        if (!MMG3D_addPROctreeRec(mesh, q, pt, q->v[i],nv))
          return 0;
        q->nbVer--;
      }
      if (!MMG3D_addPROctreeRec(mesh, q, ver, no, nv))
        return 0;
      q->nbVer--;
      MMG5_DEL_MEM(mesh,q->v);

    }else // Recursive call in the corresponding sub cell
    {
      quadrant = 0.;
      for ( i = 0; i<dim; i++)
      {
        quadrant += ((double) (ver[i]>0.5))*(1<<i);
        ver[i] -= ((double) (ver[i]>0.5))*0.5;
        ver[i] *= 2;
      }

      q->nbVer++;
      if (!MMG3D_addPROctreeRec(mesh, &(q->branches[(int)quadrant]), ver, no, nv))
        return 0;
    }
  }else // maximum PROctree depth reached
  {
    if (q->nbVer < nv)
    {
      if(q->nbVer == 0) // first allocation
      {
        MMG5_ADD_MEM(mesh,sizeof(MMG5_int),"PROctree vertices table",
                      return 0);
        MMG5_SAFE_MALLOC(q->v,1,MMG5_int,return 0);
      }
      else if(!(q->nbVer & (q->nbVer - 1))) //is a power of 2 -> normal reallocation
      {
        sizeRealloc = q->nbVer;
        sizeRealloc <<= 1;
        MMG5_ADD_MEM(mesh,(sizeRealloc-sizeRealloc/2)*sizeof(MMG5_int),"PROctree realloc",
                      return 0);
        MMG5_SAFE_REALLOC(q->v,q->nbVer,sizeRealloc,MMG5_int,"PROctree",return 0);
      }
    }
    else if (q->nbVer%nv == 0) // special reallocation of the vertex list because it is at maximum depth
    {
      MMG5_ADD_MEM(mesh,nv*sizeof(MMG5_int),"PROctree realloc",
                    return 0);
      MMG5_SAFE_REALLOC(q->v,q->nbVer,q->nbVer+nv,MMG5_int,"PROctree",return 0);
    }

    q->v[q->nbVer] = no;
    q->nbVer++;
  }

  return 1;
}

/**
 * \param pointer to the mesh structure
 * \param q pointer to the global PROctree structure
 * \param no index of the point to add to the PROctree
 *
 * Add the vertex of index \a no to the PROctree.
 *
 */
int MMG3D_addPROctree(MMG5_pMesh mesh, MMG3D_pPROctree q, const MMG5_int no)
{
  double    pt[3];
  int       dim;

  dim = mesh->dim;
  assert(no<=mesh->np);
  memcpy(&pt, mesh->point[no].c ,dim*sizeof(double));
  if (!MMG3D_addPROctreeRec(mesh, q->q0, pt , no, q->nv))
  {
    return 0;
  }

  return 1;
}

/**
 * \param q pointer to a terminal PROctree cell (containing vertex)
 * \param no index of the point to delete from the PROctree
 * \return 1 if ok 0 if memory saturated
 *
 * Delete the vertex of index \a no from the terminal PROctree cell, merge
 * the cells if necessary.
 *
 */
int MMG3D_delPROctreeVertex(MMG5_pMesh mesh, MMG3D_PROctree_s* q, MMG5_int indNo)
{
  MMG5_int* vTemp;
  int i;

  assert(q->v);
  assert(q->nbVer>indNo);
  for(i=0; i<q->nbVer; ++i)
    assert(q->v[i]>0);
  memmove(&q->v[indNo],&q->v[indNo+1], (q->nbVer-indNo-1)*sizeof(MMG5_int));
  --(q->nbVer);
  if (!(q->nbVer & (q->nbVer - 1)) && q->nbVer > 0) // is a power of 2
  {
    MMG5_ADD_MEM(mesh,q->nbVer*sizeof(MMG5_int),"PROctree index",
                  return 0);
    MMG5_SAFE_MALLOC(vTemp,q->nbVer,MMG5_int,return 0);
    memcpy(vTemp, q->v,q->nbVer*sizeof(MMG5_int));
    MMG5_DEL_MEM(mesh,q->v);

    q->v = vTemp;
  }
  return 1;
}

/**
 * \param q0 pointer to an PROctree cell.
 * \param q pointer to an PROctree cell.
 * \param dim dimension of the space (=3).
 * \param nv maximum number of points in an PROctree cell.
 * \param index next index in the array to be filled.
 *
 * Merge sub-branches \a q of \a q0, in their parent \a q0. \a q0 should
 * contain no more than nv vertices.
 *
 */
void MMG3D_mergeBranchesRec(MMG3D_PROctree_s* q0, MMG3D_PROctree_s* q, int dim, int nv, int* index)
{
  int i;

  if (q->v != NULL)
  {

    assert(*index+q->nbVer<=nv);

    memcpy(&(q0->v[*index]), q->v, q->nbVer*sizeof(MMG5_int));
    (*index)+= q->nbVer;
    for(i = 0; i<(*index); ++i)
      assert(q0->v[i]>0);
  }else if (q->branches != NULL)
  {
    for (i = 0; i<(1<<dim); ++i)
      MMG3D_mergeBranchesRec(q0, &(q->branches[i]), dim, nv, index);
  }
}

/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to an PROctree cell.
 * \param dim dimension of the space (=3)
 * \param nv maximum number of points in an PROctree cell.
 *
 * Merge branches that have a parent counting less than nv vertices.
 *
 */
void MMG3D_mergeBranches(MMG5_pMesh mesh,MMG3D_PROctree_s* q, int dim, int nv)
{
  int index;
  int i;

  index = 0;
  assert(q->v);
  assert(q->branches);
  assert(q->nbVer==nv);

  for (i = 0; i<(1<<dim); ++i)
  {
    MMG3D_mergeBranchesRec(q, &(q->branches[i]), dim, nv, &index);
    MMG3D_freePROctree_s(mesh,&(q->branches[i]), nv);
  }
  MMG5_DEL_MEM(mesh,q->branches);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to an PROctree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv maximum number of points in an PROctree cell.
 * \return 1 if ok 0 if memory saturated
 *
 * Delete vertex \a no from the PROctree. This function is recursively
 * called until we reach the terminal PROctree cell containing the vertex
 * \a no. At each step, the vertex coordinates are scaled such as the
 * quadrant is the [0;1]x[0;1]x[0;1] box.
 *
 */
int MMG3D_delPROctreeRec(MMG5_pMesh mesh, MMG3D_PROctree_s* q, double* ver, const MMG5_int no, const int nv)
{
  double quadrant;
  int i;
  int dim = mesh->dim;
  int nbVerTemp;

  if (q->v)
  {
    for ( i = 0; i<q->nbVer; ++i)
    {
      if (q->v[i] == no)
      {
        if (!MMG3D_delPROctreeVertex(mesh, q, i))
          return 0;
        if ( q->nbVer == 0)
        {
          MMG5_DEL_MEM(mesh,q->v);
        }
        break;
      }
    }

  }else if ( q->nbVer == nv+1)
  {
    quadrant = 0.;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }
    --q->nbVer;
    nbVerTemp = q->branches[(int)quadrant].nbVer;

    // warning: calling recursively here is not optimal
    if(!MMG3D_delPROctreeRec(mesh, &(q->branches[(int)quadrant]), ver, no, nv))
      return 0;

    if (nbVerTemp > q->branches[(int)quadrant].nbVer)
    {
      MMG5_ADD_MEM(mesh,nv*sizeof(MMG5_int),"PROctree vertices table",
                    return 0);
      MMG5_SAFE_MALLOC(q->v,nv,MMG5_int,return 0);
      MMG3D_mergeBranches(mesh,q,dim,nv);
    }else
    {
      ++q->nbVer;
    }

  }else if (q->branches != NULL)
  {
    quadrant = 0.;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }

    --q->nbVer;
    nbVerTemp = q->branches[(int)quadrant].nbVer;
    if(!MMG3D_delPROctreeRec(mesh, &(q->branches[(int)quadrant]), ver, no, nv))
      return 0;
    if (nbVerTemp <= q->branches[(int)quadrant].nbVer) // test if deletion worked
    {
      ++q->nbVer;
    }
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param q pointer to the global PROctree.
 * \param no reference of the vertex to be deleted.
 * \return 1 if ok 0 if memory saturated
 *
 * Delete the vertex \a no from the PROctree structure.
 *
 */
int MMG3D_delPROctree(MMG5_pMesh mesh, MMG3D_pPROctree q, const int no)
{
  double pt[3];
  int    dim;

  dim = mesh->dim;

  assert(MG_VOK(&mesh->point[no]));

  memcpy(&pt, mesh->point[no].c ,dim*sizeof(double));
  if(!MMG3D_delPROctreeRec(mesh, q->q0, pt , no, q->nv))
  {
    return 0;
  }
  return 1;
}


/**
 * \param q pointer to an PROctree cell
 * \param depth depth of the subtree
 * \param nv number of vertices in the subtree
 * \param dim dimension in which we work
 *
 * Print the depth \a depth of the subtree of \a q.
 *
 * \warning debug function, not safe
 */
void MMG3D_printArbreDepth(MMG3D_PROctree_s* q, int depth, int nv, int dim)
{
  int i;
  if ( q->depth < depth && q->nbVer > nv)
  {

    for (i = 0; i < (1<<dim); i++)
    {
      MMG3D_printArbreDepth(&(q->branches[i]),depth, nv, dim);
    }
  }else if (q->depth == depth)
  {
    fprintf(stdout,"%i ",q->nbVer);
  }
}

/**
 * \param q pointer to the global PROctree structure
 *
 * Print the PROctree.
 *
 * \warning debug function, not safe
 *
 */
void MMG3D_printArbre(MMG3D_pPROctree q)
{
  int dim;

  dim = 3;
  int i;
  for (i = 0; i<(int)sizeof(int)*8/dim; i++)
  {
    fprintf(stdout,"\n depth %i \n", i);
    MMG3D_printArbreDepth(q->q0, i, q->nv, dim);

  }
  fprintf(stdout,"\n end \n");
}

/**
 * \param q pointer to an PROctree cell
 * \param nv maximum number of vertices in an PROctree leaf
 * \param dim spacial dimension
 *
 * Print the PROctree.
 *
 * \warning debug function, not safe
 *
 */
void MMG3D_printSubArbre(MMG3D_PROctree_s* q, int nv, int dim)
{
  int i;
  for (i = 0; i<(int)sizeof(int)*8/dim; i++)
  {
    fprintf(stdout,"\n depth %i \n", i);
    MMG3D_printArbreDepth(q, i, nv, dim);

  }
  fprintf(stdout,"\n end \n");
}


/**
 * \param q pointer to an PROctree cell
 * \param nv maximum number of vertices in an PROctree leaf
 * \param dim dimension in which we work
 * \param s size of the PROctree
 *
 * Print the memory size of the PROctree.
 *
 * \warning debug function, not safe
 */
void MMG3D_sizeArbreRec(MMG3D_PROctree_s* q, int nv, int dim,int* s1, int* s2)
{
  int i;
  int nVer;
  if (q->branches != NULL)
  {
    for (i= 0; i <(1<<dim); i++)
    {
      MMG3D_sizeArbreRec(&(q->branches[i]),nv,dim, s1, s2);
      (*s1) +=(int)( sizeof(MMG3D_PROctree_s)+(1<<dim)*sizeof(MMG3D_PROctree_s*));
    }
  }else if(q->v != NULL)
  {
    // rounding up to the next higher power of 2
    nVer = q->nbVer;
    nVer--;
    nVer |= nVer >> 1;
    nVer |= nVer >> 2;
    nVer |= nVer >> 4;
    nVer |= nVer >> 8;
    nVer |= nVer >> 16;
    nVer++;
    nVer = (nVer < nv) ? nVer : (int)(((q->nbVer-0.1)/nv+1)*nv);
    (*s2) += nVer*sizeof(MMG5_int);
    (*s1) += sizeof(MMG3D_PROctree_s);
  }else
  {
    (*s1) += sizeof(MMG3D_PROctree_s);
  }
}

/**
 * \param q pointer to the global PROctree structure
 * \param dim dimension in which we work
 *
 * \return the size of the tree or NULL pointer if fail
 *
 * Print the PROctree memory size.
 *
 * \warning debug function, not safe
 *
 */
int* MMG3D_sizeArbre(MMG3D_pPROctree q,int dim)
{
  int *s;
  MMG5_SAFE_MALLOC(s,2, int,return 0);
  s[0] = 0;
  s[1] = 0;
  MMG3D_sizeArbreRec(q->q0, q->nv, dim, &s[0], &s[1]);
  return s;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param PROctree pointer to the PROctree structure.
 * \param ip index of point to check.
 * \param lmax threshold to check minimal distance between points.
 *
 * \return 1 if we can insert \a ip, 0 if we cannot insert the point
 * \return -1 if fail because of memory.
 *
 * Check if the vertex \a ip is not too close from another one (for an isotropic
 * metric).
 *
 */
int MMG3D_PROctreein_iso(MMG5_pMesh mesh,MMG5_pSol sol,MMG3D_pPROctree PROctree,MMG5_int ip,double lmax) {
  MMG5_pPoint      ppt,pp1;
  MMG3D_PROctree_s **lococ;
  double           d2,ux,uy,uz,hpi,hp1,hpi2,methalo[6];
  int              i,j;
  MMG5_int         ip1;
  int              ncells;
  double           ani[6];
  //double          dmax;

  ani[0] = sol->m[ip];
  ani[3] = sol->m[ip];
  ani[5] = sol->m[ip];
  ani[1] = 0;
  ani[2] = 0;
  ani[4] = 0;

  lococ = NULL;
  ppt = &mesh->point[ip];
  // dmax = MG_MAX(0.1,2-lmax);
  //hpi = dmax*sol->m[ip];
  hpi = lmax*sol->m[ip];
  hp1 = hpi*hpi;

  /* methalo is the box that we want to intersect with the PROctree, thus, the limit
   * of the filter. We give: the coordinates of one of the corner of the box and
   * the length of the box in each direction. */
  methalo[0] = ppt->c[0] - hpi;
  methalo[1] = ppt->c[1] - hpi;
  methalo[2] = ppt->c[2] - hpi;
  methalo[3] = methalo[4] = methalo[5] = 2.*hpi;

  ncells = MMG3D_getListSquare(mesh, ani, PROctree, methalo, &lococ);
  if (ncells < 0)
  {
    if ( lococ )
      MMG5_DEL_MEM(mesh,lococ);
    return -1;
  }
  /* Check the PROctree cells */
  for ( i=0; i<ncells; ++i )
  {
    for (j=0; j<lococ[i]->nbVer; ++j)
    {

      ip1  = lococ[i]->v[j];
      pp1  = &mesh->point[ip1];

      hpi2 = lmax * sol->m[ip1];
      //hpi2 = dmax * sol->m[ip1];

      ux = pp1->c[0] - ppt->c[0];
      uy = pp1->c[1] - ppt->c[1];
      uz = pp1->c[2] - ppt->c[2];

      d2 = ux*ux + uy*uy + uz*uz;

      if ( d2 < hp1 || d2 < hpi2*hpi2 )
      {
        MMG5_DEL_MEM(mesh,lococ);
        return 0;
      }
    }
  }
  MMG5_DEL_MEM(mesh,lococ);
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param PROctree pointer to the PROctree structure.
 * \param ip index of point to check.
 * \param lmax threshold to check minimal distance between points.
 *
 * \return 1 if we can insert \a ip, 0 otherwise
 * \return -1 if fail due to lack of memory.
 *
 * Check if the vertex \a ip is not too close from another one (for an
 * anisotropic metric).
 *
 */
int MMG3D_PROctreein_ani(MMG5_pMesh mesh,MMG5_pSol sol,MMG3D_pPROctree PROctree,MMG5_int ip,double lmax) {
  MMG5_pPoint      ppt,pp1;
  MMG3D_PROctree_s **lococ;
  double           d2,ux,uy,uz,methalo[6];
  double           det,dmi, *ma, *mb,m1,m2,m3,dx,dy,dz;
  int              i,j;
  MMG5_int         ip1,iadr;
  int              ncells;
  // double          dmax;

  lococ = NULL;
  ppt = &mesh->point[ip];

  iadr = ip*sol->size;
  ma   = &sol->m[iadr];
  // dmax = MG_MAX(0.1,2-lmax);
  dmi  =(lmax*lmax);
  //dmi  =dmax*dmax;

  // hpi = sol->m[ip]*dmax;

  det = ma[0] * (ma[3]*ma[5] - ma[4]*ma[4])
    - ma[1] * (ma[1]*ma[5] - ma[2]*ma[4])
    + ma[2] * (ma[1]*ma[4] - ma[3]*ma[2]);

  if ( det <= 0. ) return 1;

  det = 1.0 / det;
  m1 = ma[3]*ma[5] - ma[4]*ma[4];
  m2 = ma[0]*ma[5] - ma[2]*ma[2];
  m3 = ma[0]*ma[3] - ma[1]*ma[1];

  if ( m1<=0. || m2<=0. || m3<=0.) return 1;

   dx = lmax * sqrt(m1 * det) ;
   dy = lmax * sqrt(m2 * det) ;
   dz = lmax * sqrt(m3 * det) ;
  /* dx = dmax * sqrt(m1 * det) ; */
  /* dy = dmax * sqrt(m2 * det) ; */
  /* dz = dmax * sqrt(m3 * det) ; */

  /* methalo is the box that we want to intersect with the PROctree, thus, the limit
   * of the filter. We give: the coordinates of one of the corner of the box and
   * the length of the box in each direction. */
  methalo[0] = ppt->c[0] - dx;
  methalo[1] = ppt->c[1] - dy;
  methalo[2] = ppt->c[2] - dz;
  methalo[3] = 2*dx;
  methalo[4] = 2*dy;
  methalo[5] = 2*dz;

  // this function allocates lococ, it has to be deleted after the call
  ncells = MMG3D_getListSquare(mesh,ma,PROctree, methalo, &lococ);
  if (ncells < 0)
  {
    MMG5_DEL_MEM(mesh,lococ);
    return -1;
  }
  /* Check the PROctree cells */
  for ( i=0; i<ncells; ++i )
  {
    for (j=0; j<lococ[i]->nbVer; ++j)
    {
      ip1  = lococ[i]->v[j];
      pp1  = &mesh->point[ip1];

      ux = pp1->c[0] - ppt->c[0];
      uy = pp1->c[1] - ppt->c[1];
      uz = pp1->c[2] - ppt->c[2];

      d2 = ma[0]*ux*ux + ma[3]*uy*uy + ma[5]*uz*uz
        + 2.0*(ma[1]*ux*uy + ma[2]*ux*uz + ma[4]*uy*uz);
      if ( d2 < dmi )
      {
        MMG5_DEL_MEM(mesh,lococ);
        return 0;
      }
      else
      {
        iadr = ip1*sol->size;
        mb   = &sol->m[iadr];
        d2   = mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz
          + 2.0*(mb[1]*ux*uy + mb[2]*ux*uz + mb[4]*uy*uz);
        if ( d2 < dmi ) {
          MMG5_DEL_MEM(mesh,lococ);
          return 0;
        }
      }
    }
  }

  MMG5_DEL_MEM(mesh,lococ);
  return 1;
}
