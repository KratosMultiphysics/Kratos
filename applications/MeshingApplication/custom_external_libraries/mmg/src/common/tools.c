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
 * \file common/tools.c
 * \brief Various tools for the mmg applications.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgcommon.h"

/**
 * \param mesh pointer toward the mesh stucture.
 * \param ip1 first point of face.
 * \param ip2 second point of face.
 * \param ip3 third point of face.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute non-normalized face normal given three points on the surface.
 *
 */
inline int _MMG5_nonUnitNorPts(MMG5_pMesh mesh,
                                int ip1,int ip2, int ip3,double *n) {
  MMG5_pPoint   p1,p2,p3;
  double        abx,aby,abz,acx,acy,acz;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  p3 = &mesh->point[ip3];

  /* area */
  abx = p2->c[0] - p1->c[0];
  aby = p2->c[1] - p1->c[1];
  abz = p2->c[2] - p1->c[2];

  acx = p3->c[0] - p1->c[0];
  acy = p3->c[1] - p1->c[1];
  acz = p3->c[2] - p1->c[2];

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;

  return(1);
}

/**
 * \param mesh pointer toward the mesh stucture.
 * \param pt triangle for which we compute the surface.
 * \return the computed surface
 *
 * Compute non-oriented surface area of a triangle.
 *
 */
inline double _MMG5_nonorsurf(MMG5_pMesh mesh,MMG5_pTria pt) {
  double   n[3];
  int      ip1, ip2, ip3;

  ip1 = pt->v[0];
  ip2 = pt->v[1];
  ip3 = pt->v[2];

  _MMG5_nonUnitNorPts(mesh,ip1,ip2,ip3,n);

  return(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
}
/**
 * \param mesh pointer toward the mesh stucture.
 * \param ip1 first point of face.
 * \param ip2 second point of face.
 * \param ip3 third point of face.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute normalized face normal given three points on the surface.
 *
 */
inline int _MMG5_norpts(MMG5_pMesh mesh,int ip1,int ip2, int ip3,double *n) {
  double   dd,det;

  _MMG5_nonUnitNorPts(mesh,ip1,ip2,ip3,n);

  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

  if ( det < _MMG5_EPSD2 )  return(0);

  dd = 1.0 / sqrt(det);
  n[0] *= dd;
  n[1] *= dd;
  n[2] *= dd;

  return(1);
}

/**
 * \param mesh pointer toward the mesh stucture.
 * \param pt pointer toward the triangle structure.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute triangle normal.
 *
 */
inline int _MMG5_nortri(MMG5_pMesh mesh,MMG5_pTria pt,double *n) {

  return(_MMG5_norpts(mesh,pt->v[0],pt->v[1],pt->v[2],n));

}

/* Compute product R*M*tR when M is symmetric */
inline int _MMG5_rmtr(double r[3][3],double m[6], double mr[6]){
  double n[3][3];

  n[0][0] = m[0]*r[0][0] + m[1]*r[0][1] + m[2]*r[0][2];
  n[1][0] = m[1]*r[0][0] + m[3]*r[0][1] + m[4]*r[0][2];
  n[2][0] = m[2]*r[0][0] + m[4]*r[0][1] + m[5]*r[0][2];

  n[0][1] = m[0]*r[1][0] + m[1]*r[1][1] + m[2]*r[1][2];
  n[1][1] = m[1]*r[1][0] + m[3]*r[1][1] + m[4]*r[1][2];
  n[2][1] = m[2]*r[1][0] + m[4]*r[1][1] + m[5]*r[1][2];

  n[0][2] = m[0]*r[2][0] + m[1]*r[2][1] + m[2]*r[2][2];
  n[1][2] = m[1]*r[2][0] + m[3]*r[2][1] + m[4]*r[2][2];
  n[2][2] = m[2]*r[2][0] + m[4]*r[2][1] + m[5]*r[2][2];

  mr[0] = r[0][0]*n[0][0] + r[0][1]*n[1][0] + r[0][2]*n[2][0];
  mr[1] = r[0][0]*n[0][1] + r[0][1]*n[1][1] + r[0][2]*n[2][1];
  mr[2] = r[0][0]*n[0][2] + r[0][1]*n[1][2] + r[0][2]*n[2][2];
  mr[3] = r[1][0]*n[0][1] + r[1][1]*n[1][1] + r[1][2]*n[2][1];
  mr[4] = r[1][0]*n[0][2] + r[1][1]*n[1][2] + r[1][2]*n[2][2];
  mr[5] = r[2][0]*n[0][2] + r[2][1]*n[1][2] + r[2][2]*n[2][2];

  return(1);
}

/**
 * \param n pointer toward the vector that we want to send on the third vector
 * of canonical basis.
 * \param r computed rotation matrix.
 *
 * Compute rotation matrix that sends vector \a n to the third vector of
 * canonical basis.
 *
 */
inline int _MMG5_rotmatrix(double n[3],double r[3][3]) {
  double aa,bb,ab,ll,l,cosalpha,sinalpha;

  aa = n[0]*n[0];
  bb = n[1]*n[1];
  ab = n[0]*n[1];
  ll = aa+bb;
  cosalpha = n[2];
  sinalpha = sqrt(1.0- MG_MIN(1.0,cosalpha*cosalpha));

  /* No rotation needed in this case */
  if ( ll < _MMG5_EPS ) {
    if ( n[2] > 0.0 ) {
      r[0][0] = 1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = 1.0;
    }
    else {
      r[0][0] = -1.0 ; r[0][1] = 0.0 ; r[0][2] = 0.0;
      r[1][0] = 0.0 ; r[1][1] = 1.0 ; r[1][2] = 0.0;
      r[2][0] = 0.0 ; r[2][1] = 0.0 ; r[2][2] = -1.0;
    }
  }
  else {
    l = sqrt(ll);

    r[0][0] = (aa*cosalpha + bb)/ll;
    r[0][1] = ab*(cosalpha-1)/ll;
    r[0][2] = -n[0]*sinalpha/l;
    r[1][0] = r[0][1];
    r[1][1] = (bb*cosalpha + aa)/ll;
    r[1][2] = -n[1]*sinalpha/l;
    r[2][0] = n[0]*sinalpha/l;
    r[2][1] = n[1]*sinalpha/l;
    r[2][2] = cosalpha;
  }
  return(1);
}

/**
 * \param m pointer toward a 3x3 symetric matrix
 * \param mi pointer toward the computed 3x3 matrix.
 *
 * Invert \a m (3x3 symetric matrix) and store the result on \a mi
 *
 */
int _MMG5_invmat(double *m,double *mi) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check diagonal matrices */
  vmax = fabs(m[1]);
  maxx = fabs(m[2]);
  if( maxx > vmax ) vmax = maxx;
  maxx = fabs(m[4]);
  if( maxx > vmax ) vmax = maxx;
  if ( vmax < _MMG5_EPS ) {
    mi[0]  = 1./m[0];
    mi[3]  = 1./m[3];
    mi[5]  = 1./m[5];
    mi[1] = mi[2] = mi[4] = 0.0;
    return(1);
  }

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<6; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);
  /* compute sub-dets */
  aa  = m[3]*m[5] - m[4]*m[4];
  bb  = m[4]*m[2] - m[1]*m[5];
  cc  = m[1]*m[4] - m[2]*m[3];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < _MMG5_EPS3 )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[1] = bb*det;
  mi[2] = cc*det;
  mi[3] = (m[0]*m[5] - m[2]*m[2])*det;
  mi[4] = (m[1]*m[2] - m[0]*m[4])*det;
  mi[5] = (m[0]*m[3] - m[1]*m[1])*det;

  return(1);
}

/**
 * \param m initial matrix.
 * \param mi inverted matrix.
 *
 * Invert 3x3 non-symmetric matrix.
 *
 */
int _MMG5_invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);

  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < _MMG5_EPSD )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;

  return(1);
}

/**
 * \param a matrix to invert.
 * \param b last member.
 * \param r vector of unknowns.
 * \return 0 if fail, 1 otherwise.
 *
 * Solve \f$ 3\times 3\f$ symmetric system \f$ A . r = b \f$.
 *
 */
inline int _MMG5_sys33sym(double a[6], double b[3], double r[3]){
  double ia[6],as[6],det,m;
  int    i;

  /* Multiply matrix by a constant coefficient for stability purpose (because of the scaling) */
  m = fabs(a[0]);
  for(i=1;i<6;i++){
    if(fabs(a[i])<m && fabs(a[i]) > 0.){
      m = fabs(a[i]);
    }
  }

  if(m < _MMG5_EPSD) {
    return(0);
  }

  m = 1.0/m;

  for(i=0;i<6;i++){
    as[i] = a[i]*m;
  }

  det = as[0]*(as[3]*as[5]-as[4]*as[4]) - as[1]*(as[1]*as[5]-as[2]*as[4]) \
    + as[2]*(as[1]*as[4]-as[2]*as[3]);

  if(fabs(det) < _MMG5_EPSD)
    return(0);

  det = 1.0/det;

  ia[0] = (as[3]*as[5]-as[4]*as[4]);
  ia[1] = - (as[1]*as[5]-as[2]*as[4]);
  ia[2] = (as[1]*as[4]-as[2]*as[3]);
  ia[3] = (as[0]*as[5]-as[2]*as[2]);
  ia[4] = -(as[0]*as[4]-as[2]*as[1]);
  ia[5] = (as[0]*as[3]-as[1]*as[1]);

  r[0] = ia[0]*b[0] + ia[1]*b[1] + ia[2]*b[2];
  r[1] = ia[1]*b[0] + ia[3]*b[1] + ia[4]*b[2];
  r[2] = ia[2]*b[0] + ia[4]*b[1] + ia[5]*b[2];

  r[0]*=(det*m);
  r[1]*=(det*m);
  r[2]*=(det*m);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param fileName pointer toward the file name.
 *
 * Debug function (not use in clean code): write mesh->tria structure in file.
 *
 */
void _MMG5_printTria(MMG5_pMesh mesh,char* fileName) {
  MMG5_pTria ptt;
  int   k;
  FILE  *inm;

  inm = fopen(fileName,"w");

  fprintf(inm,"----------> %d TRIANGLES <----------\n",mesh->nt);
  for(k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    fprintf(inm,"num %d -> %d %d %d\n",k,ptt->v[0],ptt->v[1],
            ptt->v[2]);
    fprintf(inm,"ref   -> %d\n",ptt->ref);
    fprintf(inm,"tag   -> %d %d %d\n",ptt->tag[0],ptt->tag[1],ptt->tag[2]);
    fprintf(inm,"edg   -> %d %d %d\n",ptt->edg[0],ptt->edg[1],ptt->edg[2]);
    fprintf(inm,"\n");
  }
  fprintf(inm,"---------> END TRIANGLES <--------\n");
  fclose(inm);
}

/**
 * \return the available memory size of the computer.
 *
 * Compute the available memory size of the computer.
 *
 */
long long _MMG5_memSize (void) {
  long long mem;

#if (defined(__APPLE__) && defined(__MACH__))
  size_t size;

  size = sizeof(mem);
  if ( sysctlbyname("hw.memsize",&mem,&size,NULL,0) == -1)
    return(0);

#elif defined(__unix__) || defined(__unix) || defined(unix)
  mem = ((long long)sysconf(_SC_PHYS_PAGES))*
    ((long long)sysconf(_SC_PAGE_SIZE));
#elif defined(_WIN16) || defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) || defined(__TOS_WIN__) || defined(__WINDOWS__)
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  // status.ullTotalPhys is an unsigned long long: we must check that it fits inside a long long
  mem = status.ullTotalPhys & LLONG_MAX;
  if (mem == status.ullTotalPhys) return(mem);
  else return(LLONG_MAX);
#else
  fprintf(stderr,"  ## WARNING: UNKNOWN SYSTEM, RECOVER OF MAXIMAL MEMORY NOT AVAILABLE.\n");
  return(0);
#endif

  return(mem);
}

/**
*Safe cast into a long */
inline
long _MMG5_safeLL2LCast(long long val)
{
  long tmp_l;

  tmp_l  = (long)(val);

  if ( (long long)(tmp_l) != val ) {
        fprintf(stderr,"  ## Error:");
        fprintf(stderr," unable to cast value.n");
        exit(EXIT_FAILURE);
  }
  return(tmp_l);
}



/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param node pointer toward a _MMG5_iNode (cell for linked list)
 * \return 1 if we can alloc the node \a node, 0 otherwise.
 *
 * Node allocation.
 *
 */
inline
int _MMG5_Alloc_inode( MMG5_pMesh mesh, _MMG5_iNode **node ) {

  _MMG5_ADD_MEM(mesh,sizeof(_MMG5_iNode),"boundary reference node",
                return(0););

  _MMG5_SAFE_MALLOC(*node,1,_MMG5_iNode);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param liLi pointer toward the address of the root of the linked list.
 * \param val value to add to the linked list.
 * \return 1 if the node is inserted, 0 if the node is not inserted, -1 if fail.
 *
 * Add a node with value \a val to a sorted linked list with unique entries.
 *
 * \remark as the linked list had unique entries, we don't insert a node if it
 * exists.
 *
 */
inline
int _MMG5_Add_inode( MMG5_pMesh mesh, _MMG5_iNode **liLi, int val ) {
  _MMG5_iNode  *newNode, *cur;

  cur = *liLi;

  /* Travel through the linked list and search if the value val exist or, if
   * not, where to insert it */
  if ( cur ) {
    if ( val < (*liLi)->val ) {
      /* Add a value at the list head */
      if ( !_MMG5_Alloc_inode(mesh,&newNode) ) return(-1);

      newNode->val = val;
      newNode->nxt = (*liLi);

      (*liLi) = newNode;

      return 1;

    }
    else if (val == (*liLi)->val ) return(0);

    while ( cur->nxt && ( val >= (cur->nxt)->val) )
      cur = cur->nxt;

    if ( val == cur->val ) return(0);

    if ( !_MMG5_Alloc_inode(mesh,&newNode) ) return(-1);

    newNode->val = val;
    newNode->nxt = cur->nxt;
    cur->nxt = newNode;
  }
  else {
    if ( !_MMG5_Alloc_inode(mesh,&newNode) ) return(-1);

    newNode->val = val;
    newNode->nxt = NULL;

    *liLi = newNode;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param liLi pointer toward the root of the linked list.
 *
 * Free the memory used by the linked list whose root is \a liLi.
 *
 */
inline
void _MMG5_Free_ilinkedList( MMG5_pMesh mesh, _MMG5_iNode *liLi ) {
  _MMG5_iNode *cur,*nxt;

  cur = liLi;
  while (cur) {
    nxt = cur;
    cur = cur->nxt;

    _MMG5_DEL_MEM(mesh,nxt,sizeof(_MMG5_iNode));
  }
}


/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param node pointer toward a _MMG5_dNode (cell for linked list)
 * \return 1 if we can alloc the node \a node, 0 otherwise.
 *
 * Node allocation.
 *
 */
inline
int _MMG5_Alloc_dnode( MMG5_pMesh mesh, _MMG5_dNode **node ) {

  _MMG5_ADD_MEM(mesh,sizeof(_MMG5_dNode),"node for hausdorff eval",
                return(0););

  _MMG5_SAFE_MALLOC(*node,1,_MMG5_dNode);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param liLi pointer toward the address of the root of the linked list.
 * \param k integer value to add to the linked list.
 * \param val real value to add to the linked list.
 * \return 1 if the node is inserted, 0 if the node is not inserted, -1 if fail.
 *
 * Add a node with integer value \a k and real value \a val to a sorted linked
 * list with unique entries.
 *
 * \remark as the linked list had unique entries, we don't insert a node if it
 * exists.
 *
 */
inline
int _MMG5_Add_dnode( MMG5_pMesh mesh, _MMG5_dNode **liLi, int k, double val ) {
  _MMG5_dNode  *newNode, *cur;

  cur = *liLi;

  /* Travel through the linked list and search if the value val exist or, if
   * not, where to insert it */
  if ( cur ) {
    if ( val < (*liLi)->val ) {
      /* Add a value at the list head */
      if ( !_MMG5_Alloc_dnode(mesh,&newNode) ) return(-1);

      newNode->val = val;
      newNode->nxt = (*liLi);

      (*liLi) = newNode;

      return 1;

    }
    else if (val == (*liLi)->val ) return(0);

    while ( cur->nxt && ( val >= (cur->nxt)->val) )
      cur = cur->nxt;

    if ( val == cur->val ) return(0);

    if ( !_MMG5_Alloc_dnode(mesh,&newNode) ) return(-1);

    newNode->val = val;
    newNode->nxt = cur->nxt;
    cur->nxt = newNode;
  }
  else {
    if ( !_MMG5_Alloc_dnode(mesh,&newNode) ) return(-1);

    newNode->val = val;
    newNode->nxt = NULL;

    *liLi = newNode;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param liLi pointer toward the root of the linked list.
 *
 * Free the memory used by the linked list whose root is \a liLi.
 *
 */
inline
void _MMG5_Free_dlinkedList( MMG5_pMesh mesh, _MMG5_dNode *liLi ) {
  _MMG5_dNode *cur,*nxt;

  cur = liLi;
  while (cur) {
    nxt = cur;
    cur = cur->nxt;

    _MMG5_DEL_MEM(mesh,nxt,sizeof(_MMG5_dNode));
  }
}

/** Compute 3 * 3 determinant : det(c1-c0,c2-c0,v) */
inline double _MMG5_det3pt1vec(double c0[3],double c1[3],double c2[3],double v[3]) {
    double m00,m10,m20,m01,m11,m21,det;

    m00 = c1[0] - c0[0] ; m01 = c2[0] - c0[0];
    m10 = c1[1] - c0[1] ; m11 = c2[1] - c0[1];
    m20 = c1[2] - c0[2] ; m21 = c2[2] - c0[2];
    det = v[0]*(m10*m21 - m20*m11) -v[1]*(m00*m21-m20*m01) + v[2]*(m00*m11-m10*m01);

    return(det);
}

/** Compute 3 * 3 determinant : det(c1-c0,c2-c0,c3-c0) */
inline double _MMG5_det4pt(double c0[3],double c1[3],double c2[3],double c3[3]) {
  double m[3];

  m[0] = c3[0] - c0[0];
  m[1] = c3[1] - c0[1];
  m[2] = c3[2] - c0[2];

  return( _MMG5_det3pt1vec(c0,c1,c2,m) );
}

/**
 * \param point Pointer toward the points array
 * \param v pointer toward the point indices
 *
 * \return the oriented volume of tetra
 *
 * Compute oriented volume of a tetrahedron
 *
 */
inline double _MMG5_orvol(MMG5_pPoint point,int *v) {
    MMG5_pPoint  p0,p1,p2,p3;

    p0 = &point[v[0]];
    p1 = &point[v[1]];
    p2 = &point[v[2]];
    p3 = &point[v[3]];

    return(_MMG5_det4pt(p0->c,p1->c,p2->c,p3->c));
}


/**
 * \param a point coordinates
 * \param b point coor
 * \param c point coor
 *
 * Compute tria area.
 *
 */
double MMG2_quickarea(double a[2],double b[2],double c[2]) {
  double     abx,aby,acx,acy;//,bcx,bcy;
  double     aire;

  abx = b[0] - a[0];
  aby = b[1] - a[1];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  // bcx = c[0] - b[0];
  // bcy = c[1] - b[1];
  //
  /* orientation */
  aire = abx*acy - aby*acx;

  return(aire);
}
