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

/**
 * \file octree_3d.c
 * \brief Tools for local search around coordinates based on octree.
 * \author Jean Mercat (Inria/UBordeaux)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * An octree of the nodes is created and used for local neighbor search.
 * This helps deciding if a position is too close to other nodes to refine
 * with an insertion of a new node.
 *
 * commande test : ctest -D Experimental -VV avec tunel ssh ouvert
 * ssh : ssh -f -L 2000:vulcain.bordeaux.inria.fr:80 jmercat@vulcain.bordeaux.inria.fr sleep <temps voulu>
 *
 *se connecter : localhost:2000/CDash

 */

#include "mmg3d.h"
#include <stdio.h>

/**
 * \param q pointer toward the octree cell
 *
 * Initialisation of the octree cell.
 *
 */
void _MMG3D_initOctree_s( _MMG3D_octree_s* q)
{
  q->nbVer = 0;
  q->depth = 0;
  q->v = NULL;
  q->branches = NULL;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global octree
 * \param nv maximum number of vertices in each cell before subdivision
 * \return 1 if ok 0 if memory saturated
 *
 * Initialisation of the octree cell.
 *
 */
int _MMG3D_initOctree(MMG5_pMesh mesh,_MMG3D_pOctree* q, int nv)
{
  int i;

  _MMG5_ADD_MEM(mesh,sizeof(_MMG3D_octree),"octree structure",
                return 0);
  _MMG5_SAFE_MALLOC(*q,1, _MMG3D_octree);


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

  _MMG5_ADD_MEM(mesh,sizeof(_MMG3D_octree_s),"initial octree cell",
                return 0);

  _MMG5_SAFE_MALLOC((*q)->q0,1, _MMG3D_octree_s);
  _MMG3D_initOctree_s((*q)->q0);

  for (i=1;i<=mesh->np; ++i)
  {
    if ( !MG_VOK(&mesh->point[i]) )  continue;
    if (mesh->point[i].tag & MG_BDY) continue;

    if(!_MMG3D_addOctree(mesh, (*q), i))
      return 0;

  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the octree cell
 * \param nv number of vertices in the cell subtree
 *
 * Free the octree cell.
 *
 */
void _MMG3D_freeOctree_s(MMG5_pMesh mesh,_MMG3D_octree_s* q, int nv)
{
  int nbBitsInt,depthMax,dim,i,sizTab,sizBr,nvTemp;

  dim       = mesh->dim;
  sizBr     = 1<<dim;
  nbBitsInt = sizeof(int)*8;
  depthMax  = nbBitsInt/dim - 1;

  if (q->nbVer>nv && q->depth < depthMax )
  {
    for ( i = 0; i<sizBr; i++)
    {
      _MMG3D_freeOctree_s(mesh,&(q->branches[i]), nv);
    }
    _MMG5_DEL_MEM(mesh,q->branches,sizBr*sizeof(_MMG3D_octree_s));
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

      _MMG5_DEL_MEM(mesh,q->v,nvTemp*sizeof(int));
      q->v = NULL;
      q->nbVer = 0;
    }else
    {
      if ( q->depth != depthMax )
      {
        sizTab = nv;
      }else
      {
        sizTab = (q->nbVer%nv != 0)? 1 : 0;
        sizTab = nv * ((int)(q->nbVer/nv) + sizTab);
      }
      assert(q->v);
      _MMG5_DEL_MEM(mesh,q->v,sizTab*sizeof(int));
      q->v = NULL;
      q->nbVer = 0;
    }
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward a pointer toward the global octree.
 *
 * Free the global octree structure.
 *
 */
void _MMG3D_freeOctree(MMG5_pMesh mesh,_MMG3D_pOctree *q)
{
  _MMG3D_freeOctree_s(mesh,(*q)->q0, (*q)->nv);
  _MMG5_DEL_MEM(mesh,(*q)->q0,sizeof(_MMG3D_octree_s));
  (*q)->q0 = NULL;
  _MMG5_DEL_MEM(mesh,*q,sizeof(_MMG3D_octree));
  *q = NULL;
}


/**
 * \param q pointer toward the global octree.
 * \param ver coordinates of the point.
 * \param dim space dimension (should be 3).
 * \return the integer containing the coordinates
 *
 * Get the integer containing the coordinates
 *
 */
int _MMG3D_getOctreeCoordinate(_MMG3D_pOctree q, double* ver, int dim)
{
  int s=1<<10;
  double prec = 1./(1<<30);
  int place = 0;
  int ix = floor((ver[0]-prec)*s);
  int iy = floor((ver[1]-prec)*s);
  int iz = floor((ver[2]-prec)*s);
  ix = (ix > 0) ? ix:0;
  iy = (iy > 0) ? iy:0;
  iz = (iz > 0) ? iz:0;
  int i=0;
  int j;
  for(j=9; j>=0; j--)
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
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global octree.
 * \param no index of the moved point.
 * \param newVer new coordinates for the moved point.
 * \param oldVer old coordinates for the moved point.
 * \return 1 if ok 0 if memory saturated
 *
 * Move one point in the octree structure. /!\ the vertex of index \a no
 * can have either the new or the old coordinates in the mesh but all
 * other vertice should have the same coordinates as when they were inserted
 * into the octree. (ie: one move at a time in the mesh and the octree)
 *
 */
int _MMG3D_moveOctree(MMG5_pMesh mesh, _MMG3D_pOctree q, int no, double* newVer, double* oldVer)
{
  int oldCoor, newCoor;
  double pt[3];
  int dim;

  dim = mesh->dim;

  memcpy(&pt, oldVer ,dim*sizeof(double));
  oldCoor = _MMG3D_getOctreeCoordinate(q, oldVer, dim);
  memcpy(&pt, newVer ,dim*sizeof(double));
  newCoor = _MMG3D_getOctreeCoordinate(q, mesh->point[no].c, dim);

  if (newCoor == oldCoor) {
    return 1;
  }
  else // it could be possible to combine delOctree and addOctree to keep it local...
  {
    /* delOctree */
    memcpy(&pt, oldVer ,dim*sizeof(double));
    if (!_MMG3D_delOctreeRec(mesh, q->q0, pt , no, q->nv))
      return 0;

    /* addOctree */
    memcpy(&pt, newVer ,dim*sizeof(double));
    if(!_MMG3D_addOctreeRec(mesh, q->q0, pt , no, q->nv))
      return 0;
  }
  return 1;
}

/**
 *
 * \param cellCenter 3 coordinates of the center of the octree cell to test.
 * \param l size of the cell
 * \param zoneCenter 3 coordinates of the center of the search zone
 * \param radius of the search zone
 * \return wether the cell is included in the search zone.
 *
 */
int _MMG3D_isCellIncluded(double* cellCenter, double l, double* zoneCenter, double l0)
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
void _MMG3D_placeInListDouble(double* distList, double dist, int index, int size)
{
  memmove(&(distList[index+2]),&(distList[index+1]),(size-(index+1))*sizeof(double));
  distList[index+1] = dist;
}

/**
 * \param qList list of pointer on octree.
 * \param q pointer on octree to be inserted in the list.
 * \param index position of the element before the place where \a q
 * should be inserted.
 * \param size size of the list before insertion.
 *
 * Insert the pointer \a q in the list \a qList at position \a index+1.
 * Moves other data so nothing is lost. No memory check performed, this
 * function should be called with coherent parameters.
 *
*/
void _MMG3D_placeInListOctree(_MMG3D_octree_s** qlist, _MMG3D_octree_s* q, int index, int size)
{
  memmove(&(qlist[index+2]),&(qlist[index+1]),(size-(index+1))*sizeof(_MMG3D_octree_s*));
  #ifdef DEBUG
  if (index+2+(size-(index+1)>61 || index+1<0)
    fprintf(stdout, "Error: in placeInListOctree index too large %i > 61\n", index+2+(size-(index+1));
  #endif
  qlist[index+1] = q;
}

/**
 * \param distList ordered list of value from smallest to largest.
 * \param dist value to be compared to elements in the list.
 * \param indexMin minimum index of the list.
 * \param indexMax maximum index of the list.
 *
 * Returns the index of the biggest value of disList that is strictly
 * smaller than dist. Only search in the bounds of indexMin and indexMax.
 *
*/
int _MMG3D_seekIndex (double* distList, double dist, int indexMin, int indexMax)
{
  int indexMed;

  if (indexMin > indexMax)
    _MMG3D_seekIndex(distList, dist, indexMax, indexMin);
  else if (indexMax - indexMin <2)
  {
    #ifdef DEBUG
      if (indexMin >= 60 || indexMax>=60)
        fprintf(stdout,"Error: in seekIndex, index should not be that large %i %i\n",indexMin, indexMax);
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
        fprintf(stdout,"Error: in seekIndex, index should not be that large %i\n",indexMed);
    #endif

    if (dist > distList[indexMed])
      _MMG3D_seekIndex(distList, dist, indexMed, indexMax);
    else
      _MMG3D_seekIndex(distList, dist, indexMin, indexMed);
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
int _MMG3D_intersectRect(double *rectin, double *rectinout)
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

  if ( rectinout[3]<=0 || rectinout[4]<=0 || rectinout[5]<=0 ) return(0);

  return 1;
}

/**
 * \param q pointer toward the octree cell.
 * \param center coordinates of the centre of the current subtree.
 * \param rect rectangle that we want to intersect with the subtree. We define
 * it given: the coordinates of one corner of the rectange and the length of
 * the rectangle in each dimension.
 * \param qlist pointer toward the list of pointer over the sub octrees that
 *  intersect \a rect.
 * \param dist pointer toward the list of distances between center of
 * the octree cells in qlist and the last 3 elements are the coordinates
 * of the center of the whole recangle.
 * \param ani metric of the point.
 * \param l0 radius of the search zone.
 * \param nc number max of cell in the list +3 (the three last.
 * \param dim dimension =3.
 * \param index number of octree cells that intersect \a rect
 *
 * \return 0 if the rectangle doesn't intersect the octree (possible due to the
 * surface reconstruction), 1 otherwise.
 *
 * List the number of octree cells that intersect the rectangle
 * \a rect. To avoid counting of the cells, a maximum is set.
 *
 */
int _MMG3D_getListSquareRec(_MMG3D_octree_s* q, double* center, double* rect,
                            _MMG3D_octree_s*** qlist, double* dist, double* ani, double l0, int nc, int dim, int* index)
{
  double recttemp[6];
  double centertemp[3];
  int    recCenter[6];
  double l = 1./(1<<(q->depth+1));
  double distTemp;
  double x,y,z;
  int indexTemp,i,j,k,nBranch;

  // number max of octree cells listed for one search
  if ((*index)>nc-4)
    return 1;

  // check if the current cell is included in the search zone, can avoid
  // the loop over the vertices. Never occured in tests unless nc==4, in that
  // case, there is no gain in computing time.
  //~ if (q->nbVer>0 && _MMG3D_isCellIncluded(center, l, &(dist[nc-3]), l0))
  //~ {
    //~ (*index)=nc-3;
    //~ fprintf(stdout,"Included cell found\n");
    //~ return;
  //~ }

  if (q->branches==NULL && q->v != NULL)
  {
    // the vector dist is of size nc whereas qlist allows nc-3 inputs
    // so the 3 last values can contain the coordinates of the center
    // of the search volume.
    x = dist[nc-3] - center[0];
    y = dist[nc-2] - center[1];
    z = dist[nc-1] - center[2];
    //#warning should be replaced with distance in metric?
    distTemp = x*x+y*y+z*z;
    //~ #warning anisotropic distance not tested (not so important, this only reorders the cells)
    //~ distTemp = ani[0]*x*x+ani[3]*y*y+ani[5]*z*z+
                //~ 2*(ani[1]*x*y+ani[2]*x*z+ani[4]*y*z);
    if (*index > 0)
    {
      indexTemp = _MMG3D_seekIndex(dist,distTemp,0, *index-1);
      if (indexTemp+1<*index)
      {
        _MMG3D_placeInListDouble(dist, distTemp, indexTemp, *index);
        _MMG3D_placeInListOctree((*qlist), q, indexTemp, *index);
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
            if ( !_MMG3D_intersectRect(rect,recttemp) ) return 0;
            
            // set the new center
            centertemp[0] = center[0]-l/2+i*l;
            centertemp[1] = center[1]-l/2+j*l;
            centertemp[2] = center[2]-l/2+k*l;
            
            // recursive call in the branch
            if ( !_MMG3D_getListSquareRec(&(q->branches[nBranch]),
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
 * \param mesh pointer toward the mesh structure
 * \param ani metric to use for the cell ordering from closest to farthest
 * \param q pointer toward the global octree structure.
 * \param rect rectangle that we want to intersect with the subtree. We define
 * it given: the coordinates of one corner of the rectangle and the length of
 * the rectangle in each dimension.
 * \param qlist pointer toward the list of pointer over the sub octrees that
 *  intersect \a rect.
 *
 * \return index, the number of subtrees in the list, -1 if fail.
 *
 * List the number of octree cells that intersect the rectangle \a rect.
 *
 */
int _MMG3D_getListSquare(MMG5_pMesh mesh, double* ani, _MMG3D_pOctree q, double* rect,
                         _MMG3D_octree_s*** qlist)
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
  index = q->nc-3;

  _MMG5_ADD_MEM(mesh,q->nc*sizeof(_MMG3D_octree_s*),"octree cell",return -1);

  _MMG5_SAFE_MALLOC(*qlist,index,_MMG3D_octree_s*);
  _MMG5_SAFE_MALLOC(dist,index+3,double);

  // Set the center of the zone search
  dist[q->nc-3] = rect[0]+rect[3]/2;
  dist[q->nc-2] = rect[1]+rect[4]/2;
  dist[q->nc-1] = rect[2]+rect[5]/2;

  // Set the radius of the zone search (used for cell inclusion test if activated)
  l0 = rect[3]/2;

  // Initialization of the octree cell list
  for (i = 0; i<index; i++)
    (*qlist)[i] = NULL;

  index = 0;

  // Set center of the first octree cell
  for (i = 0; i < dim; ++i)
    center[i] = 0.5;

  // Avoid modification of input parameter rect
  memcpy(&rect2, rect, sizeof(double)*dim*2);

  if ( !_MMG3D_getListSquareRec(q->q0, center, rect2, qlist, dist, ani, l0,
                                q->nc, dim, &index) )
    return -1;


  if (index>q->nc-4)
  {
    _MMG5_SAFE_FREE(dist);
    return -1;
  }

  _MMG5_SAFE_FREE(dist);

  return index;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an octree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv maximum number of points in an octree cell.
 * \return 1 if ok 0 if memory saturated
 *
 * Add vertex in the suitable quadrant of the octree. This function is
 * recursively called until we reach the last one. At each step, the vertex
 * coordinates are scaled such as the quadrant is the [0;1]x[0;1]x[0;1] box.
 *
 */
int _MMG3D_addOctreeRec(MMG5_pMesh mesh, _MMG3D_octree_s* q, double* ver,
                         const int no, int nv)
{
  double   pt[3];
  int      dim, nbBitsInt,depthMax,i,j,k;
  int      quadrant,sizBr;
  int      sizeRealloc;

  nbBitsInt = sizeof(int)*8;
  dim       = mesh->dim;
  depthMax  = nbBitsInt/dim - 1; // maximum depth is to allow integer coordinates
  sizBr     = 1<<dim;

  if ( q->depth < depthMax ) // not at the maximum depth of the tree
  {
    if (q->nbVer < nv)  // not at the maximum number of vertice in the cell
    {

      if(q->nbVer == 0)  // first vertex list allocation
      {
        _MMG5_ADD_MEM(mesh,sizeof(int),"octree vertice table", return 0);
        _MMG5_SAFE_MALLOC(q->v,1,int);
      }
      else if(!(q->nbVer & (q->nbVer - 1))) //is a power of 2 -> reallocation of the vertex list
      {
        sizeRealloc = q->nbVer;
        sizeRealloc<<=1;
        _MMG5_ADD_MEM(mesh,(sizeRealloc-sizeRealloc/2)*sizeof(int),"octree realloc",
                      return 0);
        _MMG5_SAFE_REALLOC(q->v,sizeRealloc,int,"octree");
      }

      q->v[q->nbVer] = no;
      q->nbVer++;
      return 1;
    }
    else if (q->nbVer == nv && q->branches==NULL)  //vertex list at maximum -> cell subdivision
    {
      /* creation of sub-branch and relocation of vertices in the sub-branches */
      _MMG5_ADD_MEM(mesh,sizBr*sizeof(_MMG3D_octree_s),"octree branches",
                    return 0);
      _MMG5_SAFE_MALLOC(q->branches,sizBr,_MMG3D_octree_s);

      for ( i = 0; i<sizBr; i++)
      {
        _MMG3D_initOctree_s(&(q->branches[i]));
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
        if (!_MMG3D_addOctreeRec(mesh, q, pt, q->v[i],nv))
          return 0;
        q->nbVer--;
      }
      if (!_MMG3D_addOctreeRec(mesh, q, ver, no, nv))
        return 0;
      q->nbVer--;
      _MMG5_DEL_MEM(mesh,q->v,nv*sizeof(int));

    }else // Recursive call in the corresponding sub cell
    {
      quadrant = 0;
      for ( i = 0; i<dim; i++)
      {
        quadrant += ((double) (ver[i]>0.5))*(1<<i);
        ver[i] -= ((double) (ver[i]>0.5))*0.5;
        ver[i] *= 2;
      }

      q->nbVer++;
      if (!_MMG3D_addOctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv))
        return 0;
    }
  }else // maximum octree depth reached
  {
    if (q->nbVer < nv)
    {
      if(q->nbVer == 0) // first allocation
      {
        _MMG5_ADD_MEM(mesh,sizeof(int),"octree vertices table",
                      return 0);
        _MMG5_SAFE_MALLOC(q->v,1,int);
      }
      else if(!(q->nbVer & (q->nbVer - 1))) //is a power of 2 -> normal reallocation
      {
        sizeRealloc = q->nbVer;
        sizeRealloc<<=1;
        _MMG5_ADD_MEM(mesh,(sizeRealloc-sizeRealloc/2)*sizeof(int),"octree realloc",
                      return 0);
        _MMG5_SAFE_REALLOC(q->v,sizeRealloc,int,"octree");
      }
    }
    else if (q->nbVer%nv == 0) // special reallocation of the vertex list because it is at maximum depth
    {
      _MMG5_ADD_MEM(mesh,nv*sizeof(int),"octree realloc",
                    return 0);
      _MMG5_SAFE_REALLOC(q->v,q->nbVer+nv,int,"octree");
    }

    q->v[q->nbVer] = no;
    q->nbVer++;
  }

  return 1;
}

/**
 * \param pointer toward the mesh structure
 * \param q pointer toward the global octree structure
 * \param no index of the point to add to the octree
 *
 * Add the vertex of index \a no to the octree.
 *
 */
int _MMG3D_addOctree(MMG5_pMesh mesh, _MMG3D_pOctree q, const int no)
{
  double pt[3];
  int    dim;

  dim = mesh->dim;
  assert(no<=mesh->np);
  memcpy(&pt, mesh->point[no].c ,dim*sizeof(double));
  if (!_MMG3D_addOctreeRec(mesh, q->q0, pt , no, q->nv))
  {
    return 0;
  }
  memcpy(&pt, mesh->point[no].c ,dim*sizeof(double));

  return 1;
}

/**
 * \param q pointer toward a terminal octree cell (containing vertex)
 * \param no index of the point to delete from the octree
 * \return 1 if ok 0 if memory saturated
 *
 * Delete the vertex of index \a no from the terminal octree cell, merge
 * the cells if necessary.
 *
 */
int _MMG3D_delOctreeVertex(MMG5_pMesh mesh, _MMG3D_octree_s* q, int indNo)
{
  int i;
  int* vTemp;

  assert(q->v);
  assert(q->nbVer>indNo);
  for(i=0; i<q->nbVer; ++i)
    assert(q->v[i]>0);
  memmove(&q->v[indNo],&q->v[indNo+1], (q->nbVer-indNo-1)*sizeof(int));
  --(q->nbVer);
  if (!(q->nbVer & (q->nbVer - 1)) && q->nbVer > 0) // is a power of 2
  {
    _MMG5_ADD_MEM(mesh,q->nbVer*sizeof(int),"octree index",
                  return 0);
    _MMG5_SAFE_MALLOC(vTemp,q->nbVer,int);
    memcpy(vTemp, q->v,q->nbVer*sizeof(int));
    _MMG5_DEL_MEM(mesh,q->v,2*q->nbVer*sizeof(int));

    q->v = vTemp;
  }
  return 1;
}

/**
 * \param q0 pointer toward an octree cell.
 * \param q pointer toward an octree cell.
 * \param dim dimension of the space (=3).
 * \param nv maximum number of points in an octree cell.
 * \param index next index in the array to be filled.
 *
 * Merge sub-branches \a q of \a q0, in their parent \a q0. \a q0 should
 * contain no more than nv vertices.
 *
 */
void _MMG3D_mergeBranchesRec(_MMG3D_octree_s* q0, _MMG3D_octree_s* q, int dim, int nv, int* index)
{
  int i;

  if (q->v != NULL)
  {

    assert(*index+q->nbVer<=nv);

    memcpy(&(q0->v[*index]), q->v, q->nbVer*sizeof(int));
    (*index)+= q->nbVer;
    for(i = 0; i<(*index); ++i)
      assert(q0->v[i]>0);
  }else if (q->branches != NULL)
  {
    for (i = 0; i<(1<<dim); ++i)
      _MMG3D_mergeBranchesRec(q0, &(q->branches[i]), dim, nv, index);
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an octree cell.
 * \param dim dimension of the space (=3)
 * \param nv maximum number of points in an octree cell.
 *
 * Merge branches that have a parent counting less than nv vertices.
 *
 */
void _MMG3D_mergeBranches(MMG5_pMesh mesh,_MMG3D_octree_s* q, int dim, int nv)
{
  int index;
  int i;
  int sizBr;
  sizBr = 1<<dim;
  index = 0;
  assert(q->v);
  assert(q->branches);
  assert(q->nbVer==nv);

  for (i = 0; i<(1<<dim); ++i)
  {
    _MMG3D_mergeBranchesRec(q, &(q->branches[i]), dim, nv, &index);
    _MMG3D_freeOctree_s(mesh,&(q->branches[i]), nv);
  }
  _MMG5_DEL_MEM(mesh,q->branches,sizBr*sizeof(_MMG3D_octree_s));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an octree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv maximum number of points in an octree cell.
 * \return 1 if ok 0 if memory saturated
 *
 * Delete vertex \a no from the octree. This function is recursively
 * called until we reach the terminal octree cell containing the vertex
 * \a no. At each step, the vertex coordinates are scaled such as the
 * quadrant is the [0;1]x[0;1]x[0;1] box.
 *
 */
int _MMG3D_delOctreeRec(MMG5_pMesh mesh, _MMG3D_octree_s* q, double* ver, const int no, const int nv)
{
  int i;
  int quadrant;
  int dim = mesh->dim;
  int nbVerTemp;

  if (q->v)
  {
    for ( i = 0; i<q->nbVer; ++i)
    {
      if (q->v[i] == no)
      {
        if (!_MMG3D_delOctreeVertex(mesh, q, i))
          return 0;
        if ( q->nbVer == 0)
        {
          _MMG5_DEL_MEM(mesh,q->v,sizeof(int));
        }
        break;
      }
    }

  }else if ( q->nbVer == nv+1)
  {
    quadrant = 0;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }
    --q->nbVer;
    nbVerTemp = q->branches[quadrant].nbVer;

    // #warning calling recursively here is not optimal
    if(!_MMG3D_delOctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv))
      return 0;

    if (nbVerTemp > q->branches[quadrant].nbVer)
    {
      _MMG5_ADD_MEM(mesh,nv*sizeof(int),"octree vertices table",
                    return 0);
      _MMG5_SAFE_MALLOC(q->v,nv,int);
      _MMG3D_mergeBranches(mesh,q,dim,nv);
    }else
    {
      ++q->nbVer;
    }

  }else if (q->branches != NULL)
  {
    quadrant = 0;
    for ( i = 0; i<dim; ++i)
    {
      quadrant += ((double) (ver[i]>0.5))*(1<<i);
      ver[i] -= ((double) (ver[i]>0.5))*0.5;
      ver[i] *= 2;
    }

    --q->nbVer;
    nbVerTemp = q->branches[quadrant].nbVer;
    if(!_MMG3D_delOctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv))
      return 0;
    if (nbVerTemp <= q->branches[quadrant].nbVer) // test if deletion worked
    {
      ++q->nbVer;
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global octree.
 * \param no reference of the vertex to be deleted.
 * \return 1 if ok 0 if memory saturated
 *
 * Delete the vertex \a no from the octree structure.
 *
 */
int _MMG3D_delOctree(MMG5_pMesh mesh, _MMG3D_pOctree q, const int no)
{
  double pt[3];
  int    dim;

  dim = mesh->dim;

  assert(MG_VOK(&mesh->point[no]));

  memcpy(&pt, mesh->point[no].c ,dim*sizeof(double));
  if(!_MMG3D_delOctreeRec(mesh, q->q0, pt , no, q->nv))
  {
    return 0;
  }
  return 1;
}


/**
 * \param q pointer toward an octree cell
 * \param depth depth of the subtree
 * \param nv number of vertices in the subtree
 * \param dim dimension in which we work
 *
 * Print the depth \a depth of the subtree of \a q.
 *
 * \warning debug function, not safe
 */
void _MMG3D_printArbreDepth(_MMG3D_octree_s* q, int depth, int nv, int dim)
{
  int i;
  if ( q->depth < depth && q->nbVer > nv)
  {

    for (i = 0; i < (1<<dim); i++)
    {
      _MMG3D_printArbreDepth(&(q->branches[i]),depth, nv, dim);
    }
  }else if (q->depth == depth)
  {
    fprintf(stdout,"%i ",q->nbVer);
  }
}

/**
 * \param q pointer toward the global octree structure
 *
 * Print the octree.
 *
 * \warning debug function, not safe
 *
 */
void _MMG3D_printArbre(_MMG3D_pOctree q)
{
  int dim;

  dim = 3;
  int i;
  for (i = 0; i<sizeof(int)*8/dim; i++)
  {
    fprintf(stdout,"\n profondeur %i \n", i);
    _MMG3D_printArbreDepth(q->q0, i, q->nv, dim);

  }
  fprintf(stdout,"\n fin \n");
}

/**
 * \param q pointer toward an octree cell
 * \param nv maximum number of vertices in an octree leaf
 * \param dim spacial dimension
 *
 * Print the octree.
 *
 * \warning debug function, not safe
 *
 */
void _MMG3D_printSubArbre(_MMG3D_octree_s* q, int nv, int dim)
{
  int i;
  for (i = 0; i<sizeof(int)*8/dim; i++)
  {
    fprintf(stdout,"\n profondeur %i \n", i);
    _MMG3D_printArbreDepth(q, i, nv, dim);

  }
  fprintf(stdout,"\n fin \n");
}


/**
 * \param q pointer toward an octree cell
 * \param nv maximum number of vertices in an octree leaf
 * \param dim dimension in which we work
 * \param s size of the octree
 *
 * Print the memory size of the octree.
 *
 * \warning debug function, not safe
 */
void _MMG3D_sizeArbreRec(_MMG3D_octree_s* q, int nv, int dim,int* s1, int* s2)
{
  int i;
  int nVer;
  if (q->branches != NULL)
  {
    for (i= 0; i <(1<<dim); i++)
    {
      _MMG3D_sizeArbreRec(&(q->branches[i]),nv,dim, s1, s2);
      (*s1) += sizeof(_MMG3D_octree_s)+(1<<dim)*sizeof(_MMG3D_octree_s*);
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
    nVer = (nVer < nv) ? nVer : ((q->nbVer-0.1)/nv+1)*nv;
    (*s2) += nVer*sizeof(int);
    (*s1) += sizeof(_MMG3D_octree_s);
  }else
  {
    (*s1) += sizeof(_MMG3D_octree_s);
  }
}

/**
 * \param q pointer toward the global octree structure
 * \param dim dimension in which we work
 * Print the octree memory size.
 *
 * \warning debug function, not safe
 *
 */
int* _MMG3D_sizeArbre(_MMG3D_pOctree q,int dim)
{
  int *s;
  _MMG5_SAFE_MALLOC(s,2, int);
  s[0] = 0;
  s[1] = 0;
  _MMG3D_sizeArbreRec(q->q0, q->nv, dim, &s[0], &s[1]);
  return s;
}


/**
 * \param q pointer toward an octree cell
 * \param nv number of vertices in the subtree
 * \param dim dimension in which we work
 *
 * Print the memory size of the octree for point stored with a linked list.
 *
 * \warning debug function, not safe
 */
static inline
int _MMG3D_sizeArbreLinkRec(_MMG3D_octree_s* q, int nv, int dim)
{
  int sizeBranches,i;

  sizeBranches = 0;

  if (q->branches != NULL)
  {

    for (i= 0; i <(1<<dim); i++)
    {
      sizeBranches += _MMG3D_sizeArbreLinkRec(&(q->branches[i]), nv, dim)
        +sizeof(_MMG3D_octree_s)+(1<<dim)*sizeof(_MMG3D_octree_s*);
    }
    return sizeBranches;
  }else if(q->v != NULL)
  {
    return sizeof(int)+sizeof(_MMG3D_octree_s);
  }else
  {
    return sizeof(_MMG3D_octree_s);
  }
}

/**
 * \param q pointer toward the global octree
 *
 * Print the memory size of the octree for point stored with a linked list.
 *
 * \warning debug function, not safe
 */
static inline
int _MMG3D_sizeArbreLink(_MMG3D_pOctree q)
{
  int dim;

  dim = 3;
  return _MMG3D_sizeArbreLinkRec(q->q0, q->nv, dim)+q->q0->nbVer*sizeof(int);
}

static inline
int NearNeighborSquare(MMG5_pMesh mesh, double* ani, _MMG3D_pOctree q, int no, double l, int dim)
{
  MMG5_pPoint   ppt,ppt1;
  _MMG3D_octree_s** qlist;
  int ns,nver;
  double rect[2*dim];
  double lmin =10;
  double x,y,z;
  int nmin;
  int i, j;

  ppt = &mesh->point[no];
  rect[0] = ppt->c[0]-l;
  rect[1] = ppt->c[1]-l;
  rect[2] = ppt->c[2]-l;
  rect[3] = 2*l;
  rect[4] = 2*l;
  rect[5] = 2*l;
  qlist = NULL;
  ns = _MMG3D_getListSquare(mesh, ani, q, rect, &qlist);


  for (i = 0; i < ns; i++)
  {
    for (j = 0; j<qlist[i]->nbVer; j++)
    {
      nver = qlist[i]->v[j];
      if(nver != no)
      {
        ppt  = &mesh->point[nver];
        ppt1 = &mesh->point[no];

        x = ppt->c[0] - ppt1->c[0];
        y = ppt->c[1] - ppt1->c[1];
        z = ppt->c[2] - ppt1->c[2];
        x = x*x+y*y+z*z;
        if(lmin>x)
        {
          lmin = x;
          nmin = nver;
        }
      }
    }
  }

  _MMG5_SAFE_FREE(qlist);

  if (sqrt(lmin)<l)
    return nmin;
  else
    return -1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param octree pointer toward the octree structure.
 * \param ip index of point to check.
 *
 * \return 1 if we can insert \a ip, 0 otherwise
 *
 * Check if the vertex \a ip is not too close from another one (for an isotropic
 * metric).
 *
 */
int _MMG3D_octreein_iso(MMG5_pMesh mesh,MMG5_pSol sol,_MMG3D_pOctree octree,int ip,double lmax) {
  MMG5_pPoint     ppt,pp1;
  _MMG3D_octree_s **lococ;
  double          d2,ux,uy,uz,hpi,hp1,hpi2,methalo[6];
  int             ip1,i,j;
  int             ncells;
  double          ani[6];
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

  /* methalo is the box that we want to intersect with the octree, thus, the limit
   * of the filter. We give: the coordinates of one of the corner of the box and
   * the length of the box in each direction. */
  methalo[0] = ppt->c[0] - hpi;
  methalo[1] = ppt->c[1] - hpi;
  methalo[2] = ppt->c[2] - hpi;
  methalo[3] = methalo[4] = methalo[5] = 2.*hpi;

  ncells = _MMG3D_getListSquare(mesh, ani, octree, methalo, &lococ);
  if (ncells < 0)
  {

    _MMG5_DEL_MEM(mesh,lococ,octree->nc*sizeof(_MMG3D_octree_s*));
    return(0);
  }
  /* Check the octree cells */
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
        _MMG5_DEL_MEM(mesh,lococ,octree->nc*sizeof(_MMG3D_octree_s*));
        return(0);
      }
    }
  }
  _MMG5_DEL_MEM(mesh,lococ,octree->nc*sizeof(_MMG3D_octree_s*));
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param octree pointer toward the octree structure.
 * \param ip index of point to check.
 *
 * \return 1 if we can insert \a ip, 0 otherwise
 *
 * Check if the vertex \a ip is not too close from another one (for an
 * anisotropic metric).
 *
 */
int _MMG3D_octreein_ani(MMG5_pMesh mesh,MMG5_pSol sol,_MMG3D_pOctree octree,int ip,double lmax) {
  MMG5_pPoint     ppt,pp1;
  _MMG3D_octree_s **lococ;
  double          d2,ux,uy,uz,methalo[6];
  double          det,dmi, *ma, *mb,m1,m2,m3,dx,dy,dz;
  int             iadr,ip1,i,j;
  int             ncells;
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

  if ( det <= 0. ) return(1);

  det = 1.0 / det;
  m1 = ma[3]*ma[5] - ma[4]*ma[4];
  m2 = ma[0]*ma[5] - ma[2]*ma[2];
  m3 = ma[0]*ma[3] - ma[1]*ma[1];

  if ( m1<=0. || m2<=0. || m3<=0.) return(1);

   dx = lmax * sqrt(m1 * det) ;
   dy = lmax * sqrt(m2 * det) ;
   dz = lmax * sqrt(m3 * det) ;
  /* dx = dmax * sqrt(m1 * det) ; */
  /* dy = dmax * sqrt(m2 * det) ; */
  /* dz = dmax * sqrt(m3 * det) ; */

  /* methalo is the box that we want to intersect with the octree, thus, the limit
   * of the filter. We give: the coordinates of one of the corner of the box and
   * the length of the box in each direction. */
  methalo[0] = ppt->c[0] - dx;
  methalo[1] = ppt->c[1] - dy;
  methalo[2] = ppt->c[2] - dz;
  methalo[3] = 2*dx;
  methalo[4] = 2*dy;
  methalo[5] = 2*dz;

  // this function allocates lococ, it has to be deleted after the call
  ncells = _MMG3D_getListSquare(mesh,ma,octree, methalo, &lococ);
  if (ncells < 0)
  {
    _MMG5_DEL_MEM(mesh,lococ,octree->nc*sizeof(_MMG3D_octree_s*));
    //~ fprintf(stdout,"too many cells\n");
    return(0);
  }
  /* Check the octree cells */
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
        _MMG5_DEL_MEM(mesh,lococ,octree->nc*sizeof(_MMG3D_octree_s*));
        return 0;
      }
      else
      {
        iadr = ip1*sol->size;
        mb   = &sol->m[iadr];
        d2   = mb[0]*ux*ux + mb[3]*uy*uy + mb[5]*uz*uz
          + 2.0*(mb[1]*ux*uy + mb[2]*ux*uz + mb[4]*uy*uz);
        if ( d2 < dmi ) {
          _MMG5_DEL_MEM(mesh,lococ,octree->nc*sizeof(_MMG3D_octree_s*));
          return(0);
        }
      }
    }
  }

  _MMG5_DEL_MEM(mesh,lococ,octree->nc*sizeof(_MMG3D_octree_s*));
  return(1);
}
