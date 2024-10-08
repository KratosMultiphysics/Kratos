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
 * \file common/tools.c
 * \brief Various tools for the mmg libraries.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgcommon_private.h"

/**
 * \param n array size
 * \param val array of double precision floating points
 * \param perm permutation array
 *
 * naive (increasing) sorting algorithm, for very small tabs ; permutation is
 * stored in perm
 *
 * \remark to use only on very small arrays such as metric tensors (size lower
 * than int8_max).
 */
inline void MMG5_nsort(int8_t n,double *val,int8_t *perm){
    int8_t   i,j;
    int8_t   aux;

    for (i=0; i<n; i++)  perm[i] = i;

    for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
            if ( val[perm[i]] > val[perm[j]] ) {
                aux = perm[i];
                perm[i] = perm[j];
                perm[j] = aux;
            }
        }
    }
}

/**
 * \param n array size
 * \param shift shift to apply when taking array value
 * \param stride stride to apply when taking array value
 * \param val array of double precision floating points
 * \param oldval array to store input values
 * \param perm permutation array
 *
 * Naively permute a small array. Use shift and stride to eventually permute
 * matrix columns.
 *
 * \remark to use only on very small arrays such as metric tensors (size lower
 * than int8_max).
 */
inline void MMG5_nperm(int8_t n,int8_t shift,int8_t stride,double *val,double *oldval,int8_t *perm) {
  int8_t i,k;

  for( i = 0; i < n; i++ )
    oldval[i] = val[shift+i*stride];

  for( i = 0; i < n; i++ ) {
    k = perm[i];
    val[shift+i*stride] = oldval[k];
  }
}

/**
 * \param n1 first normal
 * \param n2 second normal
 * \param crit ridge threshold
 *
 * \return 1 if success, 0 if fail
 *
 * Check if the angle between n1 and n2 is larger than the ridge
 * criterion. If yes, return 1, 0 otherwise (ridge creation).
 *
 */
int MMG5_devangle(double* n1, double *n2, double crit)
{
  double dev;

  dev = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];

  if ( dev < crit ) {
    return 0;
  }

  return 1;
}

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
inline int MMG5_nonUnitNorPts(MMG5_pMesh mesh,
                                MMG5_int ip1,MMG5_int ip2, MMG5_int ip3,double *n) {
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

  return 1;
}

/**
 * \param mesh pointer toward the mesh stucture.
 * \param pt triangle for which we compute the surface.
 * \return the computed surface
 *
 * Compute non-oriented surface area of a triangle.
 *
 */
inline double MMG5_nonorsurf(MMG5_pMesh mesh,MMG5_pTria pt) {
  double    n[3];
  MMG5_int  ip1, ip2, ip3;

  ip1 = pt->v[0];
  ip2 = pt->v[1];
  ip3 = pt->v[2];

  MMG5_nonUnitNorPts(mesh,ip1,ip2,ip3,n);

  return n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
}
/**
 * \param mesh pointer toward the mesh stucture.
 * \param ip1 first point of face.
 * \param ip2 second point of face.
 * \param ip3 third point of face.
 * \param n pointer to store the computed normal.
 * \return 1 if success, 0 otherwise.
 *
 * Compute normalized face normal given three points on the surface.
 *
 */
inline int MMG5_norpts(MMG5_pMesh mesh,MMG5_int ip1,MMG5_int ip2, MMG5_int ip3,double *n) {
  double   dd,det;

  MMG5_nonUnitNorPts(mesh,ip1,ip2,ip3,n);

  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

  if ( det < MMG5_EPSD2 )  return 0;

  dd = 1.0 / sqrt(det);
  n[0] *= dd;
  n[1] *= dd;
  n[2] *= dd;

  return 1;
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
inline int MMG5_nortri(MMG5_pMesh mesh,MMG5_pTria pt,double *n) {

  return MMG5_norpts(mesh,pt->v[0],pt->v[1],pt->v[2],n);

}

/**
 * \param m matrix (stored as double array)
 *
 * Transpose a square matrix in 3D, stored as double array.
 */
void MMG5_transpose3d(double m[3][3]) {
  double swp;
  for( int8_t i = 0; i < 2; i++ ) {
    for( int8_t j = i+1; j < 3; j++ ) {
      swp = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = swp;
    }
  }
}

/**
 * \return 1 if success, 0 if fail.
 *
 * Test the transposition of a 3x3 matrix.
 */
int MMG5_test_transpose3d() {
  double a[3][3] = {{1.,2.,3.},
                    {4.,5.,6.},
                    {7.,8.,9.}};
  double b[3][3] = {{1.,4.,7.},
                    {2.,5.,8.},
                    {3.,6.,9.}};
  double maxerr;

  MMG5_transpose3d(a);

  maxerr = MMG5_test_mat_error(9,(double *)a,(double *)b);
  if( maxerr > MMG5_EPSD ) {
    fprintf(stderr,"  ## Error matrix transposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param dim size of the array
 * \param a first array
 * \param b second array
 * \param result scalar product of the two arrays
 *
 * Compute scalar product of two double precision arrays.
 */
void MMG5_dotprod(int8_t dim,double *a,double *b,double *result) {
  *result = 0.0;
  for( int8_t i = 0; i < dim; i++ )
    *result += a[i]*b[i];
}

/**
 * \return 1 if success, 0 if fail.
 *
 * Test vector scalar product.
 */
int MMG5_test_dotprod() {
  double a[4] = {1.0, 2.0,   3.0,  4.0};
  double b[4] = {1.0, 0.5, 1./3., 0.25};
  double cex = 4.0;
  double cnum,err;

  MMG5_dotprod(4,a,b,&cnum);
  err = fabs(cex-cnum);
  if( err > MMG5_EPSD ) {
    fprintf(stderr,"  ## Error scalar product: in function %s, error %e\n",
      __func__,err);
    return 0;
  }

  return 1;
}

/**
 * \param a first array
 * \param b second array
 * \param result cross product of the two arrays
 *
 * Compute cross product of two double precision arrays in 3D.
 */
void MMG5_crossprod3d(double *a,double *b,double *result) {
  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
}

/**
 * \return 1 if success, 0 if fail.
 *
 * Test vector cross product.
 */
int MMG5_test_crossprod3d() {
  double a[3] = {1., 0., 0.};
  double b[3] = {1./sqrt(2.),1./sqrt(2.),0.};
  double cex[3] = {0.,0.,1./sqrt(2.)};
  double cnum[3],maxerr;

  MMG5_crossprod3d(a,b,cnum);

  maxerr = MMG5_test_mat_error(3,cex,cnum);
  if( maxerr > MMG5_EPSD ) {
    fprintf(stderr,"  ## Error 3D cross product: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param m symmetric matrix
 * \param n symmetric matrix
 * \param mn result
 *
 * Compute product m*n (mn stored by rows for consistency with MMG5_eigenv3d).
 *
 */
void MMG5_mn(double m[6], double n[6], double mn[9] ){

  mn[0] = m[0]*n[0] + m[1]*n[1] + m[2]*n[2];
  mn[1] = m[0]*n[1] + m[1]*n[3] + m[2]*n[4];
  mn[2] = m[0]*n[2] + m[1]*n[4] + m[2]*n[5];
  mn[3] = m[1]*n[0] + m[3]*n[1] + m[4]*n[2];
  mn[4] = m[1]*n[1] + m[3]*n[3] + m[4]*n[4];
  mn[5] = m[1]*n[2] + m[3]*n[4] + m[4]*n[5];
  mn[6] = m[2]*n[0] + m[4]*n[1] + m[5]*n[2];
  mn[7] = m[2]*n[1] + m[4]*n[3] + m[5]*n[4];
  mn[8] = m[2]*n[2] + m[4]*n[4] + m[5]*n[5];

  return;
}

/**
 *
 * Test product of 3x3 symmetric matrices.
 *
 */
int MMG5_test_mn() {
  double m[6] = {1.,2.,3.,4.,5.,6.}; /* Test matrix 1 */
  double n[6] = {2.,3.,4.,5.,6.,7.}; /* Test matrix 2 */
  double mnex[9] = {20., 31., 37.,
                    36., 56., 67.,
                    45., 70., 84.}; /* Exact m*n product */
  double nmex[9] = {20., 36., 45.,
                    31., 56., 70.,
                    37., 67., 84.}; /* Exact n*m product */
  double prodnum[9],maxerr; /* Numerical approximation */

  /** Compute product m*n */
  MMG5_mn(m,n,prodnum);

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(9,mnex,prodnum);
  if( maxerr > MMG5_EPSD ) {
    fprintf(stderr,"  ## Error 3x3 symmetric matrix product m*n: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute product n*m */
  MMG5_mn(n,m,prodnum);

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(9,nmex,prodnum);
  if( maxerr > MMG5_EPSD ) {
    fprintf(stderr,"  ## Error 3x3 symmetric matrix product n*m: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}


/**
 * \param r 3x3 matrix
 * \param m symetric matrix
 * \param mr result
 *
 * \return 1
 *
 * Compute product R*M*tR when M is symmetric
 *
 */
inline int MMG5_rmtr(double r[3][3],double m[6], double mr[6]){
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

  return 1;
}

/**
 *
 * Test computation of product R*M*tR when M is symmetric
 *
 */
inline int MMG5_test_rmtr() {
  double m[6] = {111./2.,-109./2.,  89./2.,111./2.,-91./2.,111./2.}; /* Test matrix */
  double r[3][3] = {{1./sqrt(2.),1./sqrt(2.),         0.},
                    {         0.,1./sqrt(2.),1./sqrt(2.)},
                    {1./sqrt(2.),         0.,1./sqrt(2.)}}; /* Test transformation */
  double outex[6] = {1., 0., 0., 10., 0., 100.}; /* Exact result */
  double outnum[6],maxerr; /* Numerical result */

  /** Compute transformation */
  if( !MMG5_rmtr(r,m,outnum) )
    return 0;

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(6,outex,outnum);
  if( maxerr > 10.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error linear transformation of symmetric matrix: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
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
inline int MMG5_rotmatrix(double n[3],double r[3][3]) {
  double aa,bb,ab,ll,l,cosalpha,sinalpha;

  aa = n[0]*n[0];
  bb = n[1]*n[1];
  ab = n[0]*n[1];
  ll = aa+bb;
  cosalpha = n[2];
  sinalpha = sqrt(1.0- MG_MIN(1.0,cosalpha*cosalpha));

  /* No rotation needed in this case */
  if ( ll < MMG5_EPS ) {
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
  return 1;
}

/**
 *
 * Test computation of the rotation matrix that sends vector \a n to the third
 * vector of canonical basis.
 *
 */
int MMG5_test_rotmatrix() {
  double n[3] = {1./sqrt(1000101.),1000./sqrt(1000101.),10./sqrt(1000101.)}; /* Test unit vector */
  double idex[6] = {1.,0.,0.,1.,0.,1.}; /*Exact identity matrix */
  double ezex [3] = {0.,0.,1.}; /* Exact z-unit vector */
  double R[3][3],idnum[6],eznum[3],maxerr; /* Numerical quantities */

  /** Rodrigues' rotation formula (transposed to give a map from n to [0,0,1]).
   *  Input vector must be a unit vector. */
  if( !MMG5_rotmatrix(n,R) )
    return 0;


  /** Approximate z-unit vector */
  for( int8_t i = 0; i < 3; i++ ) {
    eznum[i] = 0.;
    for( int8_t j = 0; j < 3; j++ )
      eznum[i] += R[i][j]*n[j];
  }

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(3,ezex,eznum);
  if( maxerr > 10.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error vector rotation: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Check orthonormality */
  if( !MMG5_rmtr(R,idex,idnum) )
    return 0;

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(6,idex,idnum);
  if( maxerr > MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error rotation matrix orthonormality: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param m pointer toward a 3x3 symetric matrix
 * \param mi pointer toward the computed 3x3 matrix.
 *
 * Invert \a m (3x3 symetric matrix) and store the result on \a mi
 *
 */
int MMG5_invmat(double *m,double *mi) {
  double  aa,bb,cc,det,vmax,maxx;
  int     k;

  /* check diagonal matrices */
  vmax = fabs(m[1]);
  maxx = fabs(m[2]);
  if( maxx > vmax ) vmax = maxx;
  maxx = fabs(m[4]);
  if( maxx > vmax ) vmax = maxx;
  if ( vmax < MMG5_EPS ) {
    mi[0]  = 1./m[0];
    mi[3]  = 1./m[3];
    mi[5]  = 1./m[5];
    mi[1] = mi[2] = mi[4] = 0.0;
    return 1;
  }

  /* check ill-conditionned matrix */
  vmax = fabs(m[0]);
  for (k=1; k<6; k++) {
    maxx = fabs(m[k]);
    if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return 0;

  /* compute sub-dets */
  aa  = m[3]*m[5] - m[4]*m[4];
  bb  = m[4]*m[2] - m[1]*m[5];
  cc  = m[1]*m[4] - m[2]*m[3];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < MMG5_EPSD2 )  return 0;
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[1] = bb*det;
  mi[2] = cc*det;
  mi[3] = (m[0]*m[5] - m[2]*m[2])*det;
  mi[4] = (m[1]*m[2] - m[0]*m[4])*det;
  mi[5] = (m[0]*m[3] - m[1]*m[1])*det;

  return 1;
}

/**
 * \param m initial matrix.
 * \param mi inverted matrix.
 *
 * Invert 3x3 non-symmetric matrix.
 *
 */
int MMG5_invmatg(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmax,maxx;
  int     k;

  /* check ill-conditionned matrix */
  vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return 0;

  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < MMG5_EPSD )  return 0;
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

  return 1;
}

/**
 * \param m initial matrix.
 * \param mi inverted matrix.
 *
 * Invert 3x3 non-symmetric matrix stored in 2 dimensions
 *
 */
int MMG5_invmat33(double m[3][3],double mi[3][3]) {
  double  aa,bb,cc,det,vmax,maxx;
  int     k,l;

  /* check ill-conditionned matrix */
  vmax = fabs(m[0][0]);
  for (k=0; k<3; k++) {
    for (l=0; l<3; l++) {
      maxx = fabs(m[k][l]);
      if ( maxx > vmax )  vmax = maxx;
    }
  }
  if ( vmax == 0.0 )  return 0;

  /* check diagonal matrices */
  /* lower */
  vmax = fabs(m[1][0]);
  maxx = fabs(m[2][0]);
  if( maxx > vmax ) vmax = maxx;
  maxx = fabs(m[2][1]);
  if( maxx > vmax ) vmax = maxx;
  /* upper */
  maxx = fabs(m[0][1]);
  if( maxx > vmax ) vmax = maxx;
  maxx = fabs(m[0][2]);
  if( maxx > vmax ) vmax = maxx;
  maxx = fabs(m[1][2]);
  if( maxx > vmax ) vmax = maxx;

  if ( vmax < MMG5_EPS ) {
    mi[0][0]  = 1./m[0][0];
    mi[1][1]  = 1./m[1][1];
    mi[2][2]  = 1./m[2][2];
    mi[1][0] = mi[0][1] = mi[2][0] = mi[0][2] = mi[1][2] = mi[2][1] = 0.0;
    return 1;
  }

  /* compute sub-dets */
  aa = m[1][1]*m[2][2] - m[2][1]*m[1][2];
  bb = m[2][1]*m[0][2] - m[0][1]*m[2][2];
  cc = m[0][1]*m[1][2] - m[1][1]*m[0][2];
  det = m[0][0]*aa + m[1][0]*bb + m[2][0]*cc;
  if ( fabs(det) < MMG5_EPSD )  return 0;
  det = 1.0 / det;

  mi[0][0] = aa*det;
  mi[0][1] = bb*det;
  mi[0][2] = cc*det;
  mi[1][0] = (m[2][0]*m[1][2] - m[1][0]*m[2][2])*det;
  mi[1][1] = (m[0][0]*m[2][2] - m[2][0]*m[0][2])*det;
  mi[1][2] = (m[1][0]*m[0][2] - m[0][0]*m[1][2])*det;
  mi[2][0] = (m[1][0]*m[2][1] - m[2][0]*m[1][1])*det;
  mi[2][1] = (m[2][0]*m[0][1] - m[0][0]*m[2][1])*det;
  mi[2][2] = (m[0][0]*m[1][1] - m[1][0]*m[0][1])*det;

  /* Check results */
#ifndef NDEBUG
  double res[3][3];

  res[0][0] = m[0][0] * mi[0][0] + m[0][1] * mi[1][0] + m[0][2] * mi[2][0];
  res[0][1] = m[0][0] * mi[0][1] + m[0][1] * mi[1][1] + m[0][2] * mi[2][1];
  res[0][2] = m[0][0] * mi[0][2] + m[0][1] * mi[1][2] + m[0][2] * mi[2][2];

  res[1][0] = m[1][0] * mi[0][0] + m[1][1] * mi[1][0] + m[1][2] * mi[2][0];
  res[1][1] = m[1][0] * mi[0][1] + m[1][1] * mi[1][1] + m[1][2] * mi[2][1];
  res[1][2] = m[1][0] * mi[0][2] + m[1][1] * mi[1][2] + m[1][2] * mi[2][2];

  res[2][0] = m[2][0] * mi[0][0] + m[2][1] * mi[1][0] + m[2][2] * mi[2][0];
  res[2][1] = m[2][0] * mi[0][1] + m[2][1] * mi[1][1] + m[2][2] * mi[2][1];
  res[2][2] = m[2][0] * mi[0][2] + m[2][1] * mi[1][2] + m[2][2] * mi[2][2];


  assert ( ( fabs(res[0][0]-1.) < MMG5_EPS ) &&
           ( fabs(res[1][1]-1.) < MMG5_EPS ) &&
           ( fabs(res[2][2]-1.) < MMG5_EPS ) &&
           ( fabs(res[0][1]) < MMG5_EPS ) && ( fabs(res[0][2]) < MMG5_EPS ) &&
           ( fabs(res[1][2]) < MMG5_EPS ) &&
           ( fabs(res[1][0]) < MMG5_EPS ) && ( fabs(res[2][0]) < MMG5_EPS ) &&
           ( fabs(res[2][1]) < MMG5_EPS ) && "Matrix inversion" );

#endif

  return 1;
}

/**
 * \param m initial matrix.
 * \param mi inverted matrix.
 *
 * Invert 2x2 non-symmetric matrix stored in 2 dimensions
 *
 */
int MMG5_invmat22(double m[2][2],double mi[2][2]) {
  double det;

  det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
  if ( fabs(det) < MMG5_EPS )  return 0;
  det = 1.0 / det;

  mi[0][0] =  m[1][1]*det;
  mi[0][1] = -m[0][1]*det;
  mi[1][0] = -m[1][0]*det;
  mi[1][1] =  m[0][0]*det;

  return 1;
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
inline int MMG5_sys33sym(double a[6], double b[3], double r[3]){
  double ia[6],as[6],det,m;
  int    i;

  /* If possible : multiply matrix by a constant coefficient for stability
     purpose */
  m = fabs(a[0]);
  m = MG_MIN ( fabs(a[3]), m );
  m = MG_MIN ( fabs(a[5]), m );

  if ( m >= MMG5_EPSD ) {
    /* Normalization */
    m = 1.0/m;
  }
  else {
    /* Unable to normalized */
    m = 1.;
  }

  for(i=0;i<6;i++){
    as[i] = a[i]*m;
  }

  det = as[0]*(as[3]*as[5]-as[4]*as[4]) - as[1]*(as[1]*as[5]-as[2]*as[4]) \
    + as[2]*(as[1]*as[4]-as[2]*as[3]);

  if(fabs(det) < MMG5_EPSD)
    return 0;

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

  r[0] *= (det*m);
  r[1] *= (det*m);
  r[2] *= (det*m);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param fileName pointer toward the file name.
 *
 * Debug function (not use in clean code): write mesh->tria structure in file.
 *
 */
void MMG5_printTria(MMG5_pMesh mesh,char* fileName) {
  MMG5_pTria ptt;
  MMG5_int   k;
  FILE       *inm;

  inm = fopen(fileName,"w");

  fprintf(inm,"----------> %" MMG5_PRId " TRIANGLES <----------\n",mesh->nt);
  for(k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    fprintf(inm,"num %" MMG5_PRId " -> %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",k,ptt->v[0],ptt->v[1],
            ptt->v[2]);
    fprintf(inm,"ref   -> %" MMG5_PRId "\n",ptt->ref);
    fprintf(inm,"tag   -> %d %d %d\n",ptt->tag[0],ptt->tag[1],ptt->tag[2]);
    fprintf(inm,"edg   -> %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",ptt->edg[0],ptt->edg[1],ptt->edg[2]);
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
size_t MMG5_memSize (void) {
  size_t mem;

#if (defined(__APPLE__) && defined(__MACH__))
  size_t size;

  size = sizeof(mem);
  if ( sysctlbyname("hw.memsize",&mem,&size,NULL,0) == -1)
    return 0;

#elif defined(__unix__) || defined(__unix) || defined(unix)
  mem = sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGE_SIZE);

#elif defined(_WIN16) || defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) || defined(__TOS_WIN__) || defined(__WINDOWS__)
  MEMORYSTATUSEX status;

  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  // status.ullTotalPhys is an unsigned long long: check that it fits in a size_t
  mem = status.ullTotalPhys & SIZE_MAX;
  if (mem == status.ullTotalPhys) return mem;
  else return SIZE_MAX;

#else
  fprintf(stderr,"\n  ## Warning: %s: unknown system, recover of maximal memory"
          " not available.\n",__func__);
  return 0;
#endif

  return mem;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * Set the memMax value to its "true" value if memory asked by
 * user. Here the MMG5_MEMPERCENT coef is already applied on memMax.
 *
 */
void MMG5_memOption_memSet(MMG5_pMesh mesh) {

  assert ( mesh->memMax );
  if ( mesh->info.mem <= 0 ) {
    if ( mesh->memMax )
      /* maximal memory = 50% of total physical memory */
      mesh->memMax = (size_t)(MMG5_memSize()*MMG5_MEMPERCENT);
    else {
      /* default value = 800 MB */
      printf("  Maximum memory set to default value: %d MB.\n",MMG5_MEMMAX);
      mesh->memMax = MMG5_MEMMAX*MMG5_MILLION;
    }
  }
  else {
    /* memory asked by user if possible, otherwise total physical memory */
    if ( (size_t)mesh->info.mem*MMG5_MILLION > (mesh->memMax/MMG5_MEMPERCENT) && mesh->memMax ) {
      fprintf(stderr,"\n  ## Warning: %s: asking for %d MB of memory ",
              __func__,mesh->info.mem);
      fprintf(stderr,"when only %zu available.\n",mesh->memMax/MMG5_MILLION);
    }
    else {
      mesh->memMax= (size_t)mesh->info.mem*MMG5_MILLION;
    }
  }
  return;
}


/** Compute 3 * 3 determinant : det(c1-c0,c2-c0,v) */
inline double MMG5_det3pt1vec(double c0[3],double c1[3],double c2[3],double v[3]) {
    double m00,m10,m20,m01,m11,m21,det;

    m00 = c1[0] - c0[0] ; m01 = c2[0] - c0[0];
    m10 = c1[1] - c0[1] ; m11 = c2[1] - c0[1];
    m20 = c1[2] - c0[2] ; m21 = c2[2] - c0[2];
    det = v[0]*(m10*m21 - m20*m11) -v[1]*(m00*m21-m20*m01) + v[2]*(m00*m11-m10*m01);

    return det;
}

/** Compute 3 * 3 determinant : det(c1-c0,c2-c0,c3-c0) */
inline double MMG5_det4pt(double c0[3],double c1[3],double c2[3],double c3[3]) {
  double m[3];

  m[0] = c3[0] - c0[0];
  m[1] = c3[1] - c0[1];
  m[2] = c3[2] - c0[2];

  return  MMG5_det3pt1vec(c0,c1,c2,m) ;
}

/**
 * \param point Pointer toward the points array
 * \param v pointer toward the point indices
 *
 * \return the oriented volume of tetra
 *
 * Compute oriented volume of a tetrahedron (x6)
 *
 */
inline double MMG5_orvol(MMG5_pPoint point,MMG5_int *v) {
    MMG5_pPoint  p0,p1,p2,p3;

    p0 = &point[v[0]];
    p1 = &point[v[1]];
    p2 = &point[v[2]];
    p3 = &point[v[3]];

    return MMG5_det4pt(p0->c,p1->c,p2->c,p3->c);
}


/**
 * \param a point coordinates
 * \param b point coor
 * \param c point coor
 *
 * Compute tria area.
 *
 */
double MMG2D_quickarea(double a[2],double b[2],double c[2]) {
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

  return aire;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Mark all mesh vertices as unused.
 *
 */
void MMG5_mark_verticesAsUnused ( MMG5_pMesh mesh ) {
  MMG5_pPoint ppt;
  MMG5_int    k;

  for ( k=1; k<=mesh->np; k++ ) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;

    /* reset ppt->flag to detect isolated points in keep_subdomainElts */
    ppt->flag = 0;
    ppt->tag |= MG_NUL;
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param delPt function to call to delete point.
 *
 * Mark the mesh vertices that belong to triangles or quadrangles as used (for
 * Mmgs or Mmg2d).
 *
 */
void MMG5_mark_usedVertices ( MMG5_pMesh mesh,void (*delPt)(MMG5_pMesh,MMG5_int)  ) {
  MMG5_pTria  pt;
  MMG5_pQuad  pq;
  MMG5_pPoint ppt;
  int         i;
  MMG5_int    k;

  /* Preserve isolated required points */
  for ( k=1; k<=mesh->np; k++ ) {
    ppt = &mesh->point[k];

    if ( ppt->flag || !(ppt->tag & MG_REQ) ) {
      continue;
    }
    ppt->tag &= ~MG_NUL;
  }

  /* Mark points used by the connectivity */
  for ( k=1; k<=mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[ pt->v[i] ];
      ppt->tag &= ~MG_NUL;
    }
  }

  for ( k=1; k<=mesh->nquad; k++ ) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;

    for (i=0; i<4; i++) {
      ppt = &mesh->point[ pq->v[i] ];
      ppt->tag &= ~MG_NUL;
    }
  }

  /* Finally, clean point array */
  while ( (!MG_VOK(&mesh->point[mesh->np])) && mesh->np ) {
    delPt(mesh,mesh->np);
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param nsd subdomain index.
 * \param delElt function to call to delete elt.
 *
 * Remove triangles that do not belong to subdomain of index \a nsd
 *
 */
void MMG5_keep_subdomainElts ( MMG5_pMesh mesh, int nsd,
                               int (*delElt)(MMG5_pMesh,MMG5_int) ) {
  MMG5_pTria  pt;
  int         i,iv;
  MMG5_int    k,*adja,iadr,iadrv;
  int         nfac = 3; // number of faces per elt

  for ( k=1 ; k <= mesh->nt ; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    /* Mark triangle vertices as seen to be able to detect isolated points */
    mesh->point[pt->v[0]].flag = 1;
    mesh->point[pt->v[1]].flag = 1;
    mesh->point[pt->v[2]].flag = 1;

    if ( pt->ref == nsd ) continue;

    /* Update adjacency relationship: we will delete elt k so k adjacent will
     * not be adjacent to k anymore */
    if ( mesh->adja ) {
      iadr = nfac*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for ( i=0; i<nfac; ++i ) {
        iadrv = adja[i];
        if ( !iadrv ) {
          continue;
        }
        iv = iadrv%nfac;
        iadrv /= nfac;
        mesh->adja[nfac*(iadrv-1)+1+iv] = 0;
      }
    }

    /* Delete element (triangle in 2D and surface, tetra in 3D) */
    delElt(mesh,k);
  }

  return;
}

/**
 * \param nelem number of matrix elements.
 * \param m1 first matrix (single array).
 * \param m2 second matrix (single array).
 *
 * Compute maximum error between two matrices.
 *
 */
inline
double MMG5_test_mat_error( int8_t nelem,double m1[],double m2[] ) {
  double maxerr;
  int8_t k;

  /* Compute max error */
  maxerr = 0;
  for( k = 0; k < nelem; k++ )
    maxerr = MG_MAX(maxerr,fabs(m1[k] - m2[k]));

  return maxerr;
}

/**
 *
 * Test inversion of 2x2 non-symmetric matrix stored in 2 dimensions.
 *
 */
int MMG5_test_invmat22() {
  double A[2][2] = {{4.0,2.0},{1.0,1.0}}; /* Test matrix */
  double iAex[2][2] = {{0.5,-1.0},{-0.5,2.0}}; /* Analytical inverse */
  double iAnum[2][2]; /* Numerical inverse */

  /* Compute matrix inverse */
  if( !MMG5_invmat22(A,iAnum) )
    return 0;

  /* Check error in norm inf */
  double maxerr = MMG5_test_mat_error(4,(double *)iAex,(double *)iAnum);
  if( maxerr > MMG5_EPSD ) {
    fprintf(stderr,"  ## Error matrix inversion: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 *
 * Test inversion of 3x3 non-symmetric matrix stored in 2 dimensions.
 *
 */
int MMG5_test_invmat33() {
  double A[3][3] = {{7.0,2.0,1.0},{0.0,3.0,-1.0},{-3.0,4.0,-2.0}}; /* Test matrix */
  double iAex[3][3] = {{-2.0,8.0,-5.0},{3.0,-11.0,7.0},{9.0,-34.0,21.0}}; /* Analytical inverse */
  double iAnum[3][3]; /* Numerical inverse */

  /* Compute matrix inverse */
  if( !MMG5_invmat33(A,iAnum) )
    return 0;

  /* Check error in norm inf */
  double maxerr = MMG5_test_mat_error(9,(double *)iAex,(double *)iAnum);
  if( maxerr > MMG5_EPSD ) {
    fprintf(stderr,"  ## Error matrix inversion: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}
