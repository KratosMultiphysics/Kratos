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
 * \file common/mettools.c
 * \brief Metric tools for the mmg applications.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"
#include "inlined_functions_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param dim matrix size.
 * \param m matrix array.
 * \param lambda eigenvalues array.
 * \param v array of eigenvectors.
 *
 * Recompose a generic-size symmetric matrix from its eigenvalue
 * decomposition.
 * \warning Eigenvectors in Mmg are stored as matrix rows (the first dimension
 * of the double array spans the number of eigenvectors, the second dimension
 * spans the number of entries of each eigenvector).
 */
static inline
void MMG5_eigenvmat_buildsym(MMG5_pMesh mesh,int8_t dim,double m[],
                             double lambda[],double v[]) {
  int8_t i,j,k,ij;

  /* Storage of a matrix as a one-dimensional array: (dim+1)*dim/2 entries for a
   * symmetric matrix (loop on each symmetric entry only once). */
  ij = 0;

  /* Loop on matrix rows */
  for( i = 0; i < dim; i++ ) {
    /* Loop on the upper-triangular part of the matrix */
    for( j = i; j < dim; j++ ) {
      /* Initialize matrix entry */
      m[ij] = 0.0;
      /* Compute matrix entry as the recomposition of eigenvalues and
       * eigenvectors:
       *
       * M_{ij} = \sum_{k,l} V_{ik} Lambda_{kl} V_{jl} =
       *        = \sum_{k} lambda_{k} V_{ik} V_{jk}
       *
       * Eigenvectors are stored as rows in Mmg (not columns) so their indices
       * have to be exchanged when implementing the above formula. */
      for( k = 0; k < dim; k++ ) {
        m[ij] += lambda[k]*v[k*dim+i]*v[k*dim+j];
      }
      /* Go to next entry */
      ++ij;
    }
  }
  assert( ij == (dim+1)*dim/2 );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m matrix array.
 * \param dim matrix size.
 * \param lambda eigenvalues array.
 * \param v double array of right eigenvectors.
 * \param iv double array of left eigenvectors.
 * \return 1 if success, 0 if failure.
 *
 * Recompose a generic-size non-symmetric matrix from its eigenvalue
 * decomposition.
 * \warning Eigenvectors in Mmg are stored as matrix rows (the first dimension
 * of the double array spans the number of eigenvectors, the second dimension
 * spans the number of entries of each eigenvector).
 */
static inline
void MMG5_eigenvmat_buildnonsym(MMG5_pMesh mesh,int8_t dim,double m[],
                                double lambda[],double v[],double iv[]){
  int8_t i,j,k,ij;

  /* Storage of a matrix as a one-dimensional array: dim^2 entries for a
   * non-symmetric matrix. */
  ij = 0;

  /* Loop on matrix rows */
  for( i = 0; i < dim; i++ ) {
    /* Loop on matrix columns */
    for( j = 0; j < dim; j++ ) {
      /* Initialize matrix entry */
      m[ij] = 0.0;
      /* Compute matrix entry as the recomposition of eigenvalues and
       * eigenvectors:
       *
       * M_{ij} = \sum_{k,l} V_{ik} Lambda_{kl} V^{-1}_{lj} =
       *        = \sum_{k} lambda_{k} V_{ik} V^{-1}_{kj}
       *
       * Eigenvectors are stored as rows in Mmg (not columns) so their indices
       * have to be exchanged when implementing the above formula. */
      for( k = 0; k < dim; k++ ) {
        m[ij] += lambda[k]*v[k*dim+i]*iv[j*dim+k];
      }
      /* Go to next entry */
      ++ij;
    }
  }
  assert( ij == dim*dim );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m matrix array.
 * \param lambda eigenvalues array.
 * \param v double array of eigenvectors.
 * \return 1 if success, 0 if failure.
 *
 * Recompose a 2x2 symmetric matrix from its eigenvalue decomposition.
 * \warning Eigenvectors in Mmg are stored as matrix rows (the first dimension
 * of the double array spans the number of eigenvectors, the second dimension
 * spans the number of entries of each eigenvector).
 */
int MMG5_eigenvmatsym2d(MMG5_pMesh mesh,double m[],double lambda[],double v[][2]) {

  /* Build matrix */
  MMG5_eigenvmat_buildsym(mesh,2,m,lambda,(double *)v);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m matrix array.
 * \param lambda eigenvalues array.
 * \param v double array of eigenvectors.
 * \return 1 if success, 0 if failure.
 *
 * Recompose a 3x3 symmetric matrix from its eigenvalue decomposition.
 * \warning Eigenvectors in Mmg are stored as matrix rows (the first dimension
 * of the double array spans the number of eigenvectors, the second dimension
 * spans the number of entries of each eigenvector).
 */
int MMG5_eigenvmatsym3d(MMG5_pMesh mesh,double m[],double lambda[],double v[][3]) {

  /* Build matrix */
  MMG5_eigenvmat_buildsym(mesh,3,m,lambda,(double *)v);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m matrix array.
 * \param lambda eigenvalues array.
 * \param v double array of eigenvectors.
 * \return 1 if success, 0 if failure.
 *
 * Recompose a 2x2 non-symmetric matrix from its eigenvalue decomposition.
 * \warning Eigenvectors in Mmg are stored as matrix rows (the first dimension
 * of the double array spans the number of eigenvectors, the second dimension
 * spans the number of entries of each eigenvector).
 */
int MMG5_eigenvmatnonsym2d(MMG5_pMesh mesh,double m[],double lambda[],double v[][2]) {
  double iv[2][2];

  /* Compute the inverse of the eigenvectors matrix */
  if( !MMG5_invmat22(v,iv) )
    return 0;

  /* Build matrix */
  MMG5_eigenvmat_buildnonsym(mesh,2,m,lambda,(double *)v,(double *)iv);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m matrix array.
 * \param lambda eigenvalues array.
 * \param v double array of eigenvectors.
 * \return 1 if success, 0 if failure.
 *
 * Recompose a 3x3 non-symmetric matrix from its eigenvalue decomposition.
 * \warning Eigenvectors in Mmg are stored as matrix rows (the first dimension
 * of the double array spans the number of eigenvectors, the second dimension
 * spans the number of entries of each eigenvector).
 */
int MMG5_eigenvmatnonsym3d(MMG5_pMesh mesh,double m[],double lambda[],double v[][3]) {
  double iv[3][3];

  /* Compute the inverse of the eigenvectors matrix */
  if( !MMG5_invmat33(v,iv) )
    return 0;

  /* Build matrix */
  MMG5_eigenvmat_buildnonsym(mesh,3,m,lambda,(double *)v,(double *)iv);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param dim matrix size.
 * \param symmat integer flag (1 if the matrix is symmetric, 0 otherwise).
 * \param m input matrix array.
 * <param mnew output matrix array.
 * \return 1 if success, 0 if failure.
 *
 * Check the recomposition of a matrix from its numerical eigenvalue
 * decomposition.
 */
int MMG5_eigenvmat_check(MMG5_pMesh mesh,int8_t dim,int8_t symmat,double m[],
                         double mnew[]) {

  /* Compute eigendecomposition, recompose matrix from eigendecomposition */
  if( dim == 2 ) {
    double lambda[2],v[2][2];

    if( !MMG5_eigenv2d(symmat,m,lambda,v) )
      return 0;

    if( symmat ) {
      if( !MMG5_eigenvmatsym2d(mesh,mnew,lambda,v) )
        return 0;
    } else {
      if( !MMG5_eigenvmatnonsym2d(mesh,mnew,lambda,v) )
        return 0;
    }

  } else if( dim == 3 ) {
    double lambda[3],v[3][3];

    if( !MMG5_eigenv3d(symmat,m,lambda,v) )
      return 0;

    if( symmat ) {
      if( !MMG5_eigenvmatsym3d(mesh,mnew,lambda,v) )
        return 0;
    } else {
      if( !MMG5_eigenvmatnonsym3d(mesh,mnew,lambda,v) )
        return 0;
    }
  }

  return 1;
}

/**
 * \param dim square matrix size
 * \param lambda eigenvalues array
 * \param vp eigenvectors array
 * \param swap swap array
 * \param perm permutation array
 *
 * Sort and permute eigenvalues (and eigenvectors) in increasing order.
 *
 */
void MMG5_sort_eigenv( int8_t dim,double *lambda,double *vp,
                       double *swap,int8_t *perm ) {
  MMG5_nsort(dim,lambda,perm);
  MMG5_nperm(dim,0,1,lambda,swap,perm);
  for( int8_t i = 0; i < dim; i++ )
    MMG5_nperm(dim,i,dim,vp,swap,perm);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param mex test matrix array.
 * \param lambdaex exact eigenvalues array.
 * \param vpex double array of exact eigenvectors.
 * \return 1 if success, 0 if failure.
 *
 * For a 2x2 symmetric matrix, Test:
 * - the recomposition of the matrix from its exact eigendecomposition;
 * - the computation of the eigenvalues of the matrix;
 * - the computation of the eigenvectors of the matrix;
 * - the recomposition of the matrix from its numerical eigendecomposition.
 *
 */
int MMG5_test_eigenvmatsym2d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                             double vpex[][2]) {
  double mnum[3],lambdanum[2],vpnum[2][2]; /* Numerical quantities */
  double swap[2],maxerr,err;
  int8_t perm[2] = {0,1}; /* eigenvalues permutation array */


  /** Recompose matrix from its eigendecomposition */
  MMG5_eigenvmat_buildsym(mesh,2,mnum,lambdaex,(double *)vpex);

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(3,(double *)mex,(double *)mnum);
  if( maxerr > 10.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute eigendecomposition */
  if( !MMG5_eigenv2d(1,mex,lambdanum,vpnum) )
    return 0;

  /* Naively sort eigenpairs in increasing order */
  MMG5_sort_eigenv(2,lambdanum,(double *)vpnum,swap,perm);

  /* Check eigenvalues error in norm inf */
  maxerr = MMG5_test_mat_error(2,(double *)lambdaex,(double *)lambdanum);
  if( maxerr > 10*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigenvalues: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /* Check eigenvectors error through scalar product */
  maxerr = 0.;
  for( int8_t i = 0; i < 2; i++ ) {
    err = 0.;
    for( int8_t j = 0; j < 2; j++ )
      err += vpex[i][j] * vpnum[i][j];
    err = 1.-fabs(err);
    maxerr = MG_MAX(maxerr,err);
  }
  if( maxerr > MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigenvectors: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute both eigendecomposition and recomposition, and check matrix */
  if( !MMG5_eigenvmat_check(mesh,2,1,mex,mnum) )
    return 0;
  maxerr = MMG5_test_mat_error(3,(double *)mex,(double *)mnum);
  if( maxerr > 10*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigendecomposition and recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param mex test matrix array.
 * \param lambdaex exact eigenvalues array.
 * \param vpex double array of exact right eigenvectors.
 * \param ivpex double array of exact right eigenvectors inverse.
 * \return 1 if success, 0 if failure.
 *
 * For a 2x2 non-symmetric matrix, Test:
 * - the recomposition of the matrix from its exact eigendecomposition;
 * - the computation of the eigenvalues of the matrix;
 * - the computation of the eigenvectors of the matrix;
 * - the recomposition of the matrix from its numerical eigendecomposition.
 *
 */
int MMG5_test_eigenvmatnonsym2d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                                double vpex[][2],double ivpex[][2]) {
  double mnum[4],lambdanum[2],vpnum[2][2]; /* Numerical quantities */
  double swap[2],maxerr,err;
  int8_t perm[2] = {0,1}; /* eigenvalues permutation array */


  /** Recompose matrix from its eigendecomposition */
  MMG5_eigenvmat_buildnonsym(mesh,2,mnum,lambdaex,(double *)vpex,(double *)ivpex);

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(4,(double *)mex,(double *)mnum);
  if( maxerr > 100.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute eigendecomposition */
  if( !MMG5_eigenv2d(0,mex,lambdanum,vpnum) )
    return 0;

  /* Naively sort eigenpairs in increasing order */
  MMG5_sort_eigenv(2,lambdanum,(double *)vpnum,swap,perm);

  /* Check eigenvalues error in norm inf */
  maxerr = MMG5_test_mat_error(2,(double *)lambdaex,(double *)lambdanum);
  if( maxerr > 10*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigenvalues: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /* Check eigenvectors error through scalar product */
  maxerr = 0.;
  for( int8_t i = 0; i < 2; i++ ) {
    err = 0.;
    for( int8_t j = 0; j < 2; j++ )
      err += vpex[i][j] * vpnum[i][j];
    err = 1.-fabs(err);
    maxerr = MG_MAX(maxerr,err);
  }
  if( maxerr > MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigenvectors: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute both eigendecomposition and recomposition, and check matrix */
  if( !MMG5_eigenvmat_check(mesh,2,0,mex,mnum) )
    return 0;
  maxerr = MMG5_test_mat_error(4,(double *)mex,(double *)mnum);
  if( maxerr > 100*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigendecomposition and recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param mex test matrix array.
 * \param lambdaex exact eigenvalues array.
 * \param vpex double array of exact eigenvectors.
 * \return 1 if success, 0 if failure.
 *
 * For a 3x3 symmetric matrix, Test:
 * - the recomposition of the matrix from its exact eigendecomposition;
 * - the computation of the eigenvalues of the matrix;
 * - the computation of the eigenvectors of the matrix;
 *   the recomposition of the matrix from its numerical eigendecomposition.
 *
 */
int MMG5_test_eigenvmatsym3d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                             double vpex[][3]) {
  double mnum[6],lambdanum[3],vpnum[3][3]; /* Numerical quantities */
  double swap[3],maxerr,err;
  int8_t perm[3] = {0,1,2}; /* eigenvalues permutation array */


  /** Recompose matrix from its eigendecomposition */
  MMG5_eigenvmat_buildsym(mesh,3,mnum,lambdaex,(double *)vpex);

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(6,(double *)mex,(double *)mnum);
  if( maxerr > 100.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute eigendecomposition */
  if( !MMG5_eigenv3d(1,mex,lambdanum,vpnum) )
    return 0;

  /* Naively sort eigenpairs in increasing order */
  MMG5_sort_eigenv(3,lambdanum,(double *)vpnum,swap,perm);

  /* Check eigenvalues error in norm inf */
  maxerr = MMG5_test_mat_error(3,(double *)lambdaex,(double *)lambdanum);
  if( maxerr > 10*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigenvalues: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /* Check eigenvectors orthogonality (checking numerical eigenvectors against
   * the exact ones does not make sense in case of multiple eigenvalues) */
  maxerr = 0.;
  for( int8_t i = 0; i < 2; i++ ) {
    for( int8_t j = i+1; j < 3; j++ ) {
      err = 0.;
      for( int8_t k = 0; k < 3; k++ )
        err += vpnum[i][k] * vpnum[j][k];
      err = fabs(err);
      maxerr = MG_MAX(maxerr,err);
    }
  }
  if( maxerr > MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigenvectors orthogonality: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute both eigendecomposition and recomposition, and check matrix */
  if( !MMG5_eigenvmat_check(mesh,3,1,mex,mnum) )
    return 0;
  maxerr = MMG5_test_mat_error(6,(double *)mex,(double *)mnum);
  if( maxerr > 100*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigendecomposition and recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param mex test matrix array.
 * \param lambdaex exact eigenvalues array.
 * \param vpex double array of exact right eigenvectors.
 * \param ivpex double array of exact right eigenvectors inverse.
 * \return 1 if success, 0 if failure.
 *
 * For a 3x3 non-symmetric matrix, Test:
 * - the recomposition of the matrix from its exact eigendecomposition;
 * - the computation of the eigenvalues of the matrix;
 * - the computation of the eigenvectors of the matrix;
 * - the recomposition of the matrix from its numerical eigendecomposition.
 */
int MMG5_test_eigenvmatnonsym3d(MMG5_pMesh mesh,double *mex,double lambdaex[],
                                double vpex[][3],double ivpex[][3]) {
  double mnum[9],lambdanum[3],vpnum[3][3]; /* Numerical quantities */
  double swap[3],maxerr;
  int8_t perm[3] = {0,1,2}; /* eigenvalues permutation array */


  /** Recompose matrix from its eigendecomposition */
  MMG5_eigenvmat_buildnonsym(mesh,3,mnum,lambdaex,(double *)vpex,(double *)ivpex);

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(9,(double *)mex,(double *)mnum);
  if( maxerr > 1000.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }


  /** Compute eigendecomposition */
  if( !MMG5_eigenv3d(0,mex,lambdanum,vpnum) )
    return 0;

  /* Naively sort eigenpairs in increasing order */
  MMG5_sort_eigenv(3,lambdanum,(double *)vpnum,swap,perm);

  /* Check eigenvalues error in norm inf */
  maxerr = MMG5_test_mat_error(3,(double *)lambdaex,(double *)lambdanum);
  if( maxerr > 1000*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigenvalues: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  /* skip eigenvectors check */


  /** Compute both eigendecomposition and recomposition, and check matrix */
  if( !MMG5_eigenvmat_check(mesh,3,0,mex,mnum) )
    return 0;
  maxerr = MMG5_test_mat_error(9,(double *)mex,(double *)mnum);
  if( maxerr > 1000*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error matrix eigendecomposition and recomposition: in function %s, max error %e\n",
      __func__,maxerr);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param t tangent at the ridge point.
 * \param n normal at the ridge point.
 * \param dtan metric size along the tangent direction.
 * \param dv metric size along the \f$t^{}n\f$ direction.
 * \param dn metric size along the normal direction.
 * \param m computed metric at the ridge point.
 * \return 1
 *
 * Build metric tensor at a fictitious ridge point, whose normal and tangent are
 * provided.
 *
 */
inline int
MMG5_buildridmetfic(MMG5_pMesh mesh,double t[3],double n[3],double dtan,
                     double dv,double dn,double m[6]) {
  double u[3],r[3][3];

  u[0] = n[1]*t[2] - n[2]*t[1];
  u[1] = n[2]*t[0] - n[0]*t[2];
  u[2] = n[0]*t[1] - n[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(dtan,dv,dn)*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n[2];

  m[0] = dtan*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1] + dn*r[0][2]*r[0][2];
  m[1] = dtan*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1] + dn*r[0][2]*r[1][2];
  m[2] = dtan*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1] + dn*r[0][2]*r[2][2];
  m[3] = dtan*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1] + dn*r[1][2]*r[1][2];
  m[4] = dtan*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1] + dn*r[1][2]*r[2][2];
  m[5] = dtan*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1] + dn*r[2][2]*r[2][2];

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m pointer toward the first metric to intersect.
 * \param n pointer toward the second metric to intersect.
 * \param mr pointer toward the computed intersected metric.
 * \return 1.
 *
 * Compute the intersected (2 x 2) metric between metrics \a m and \a n,
 * PRESERVING the directions of \a m. Result is stored in \a mr.
 *
 */
int MMG5_intmetsavedir(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  int    i;
  double lambda[2],vp[2][2],siz,isqhmin;

  isqhmin = 1.0 / (mesh->info.hmin * mesh->info.hmin);
  MMG5_eigensym(m,lambda,vp);

  for (i=0; i<2; i++) {
    siz = n[0]*vp[i][0]*vp[i][0] + 2.0*n[1]*vp[i][0]*vp[i][1]
      + n[2]*vp[i][1]*vp[i][1];
    lambda[i] = MG_MAX(lambda[i],siz);
    lambda[i] = MG_MIN(lambda[i],isqhmin);
  }
  mr[0] = lambda[0]*vp[0][0]*vp[0][0] + lambda[1]*vp[1][0]*vp[1][0];
  mr[1] = lambda[0]*vp[0][0]*vp[0][1] + lambda[1]*vp[1][0]*vp[1][1];
  mr[2] = lambda[0]*vp[0][1]*vp[0][1] + lambda[1]*vp[1][1]*vp[1][1];

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param ux distance \f$[p0;p1]\f$ along x axis.
 * \param uy distance \f$[p0;p1]\f$ along y axis.
 * \param uz distance \f$[p0;p1]\f$ along z axis.
 * \param mr computed metric tensor.
 * \param r basis in which the metric is diagona
 *
 * \return 1 if success
 *
 * Build metric tensor at ridge point p0, when computations with respect to p1
 * are to be held. Store the basis vectors in \a r.
 *
 * \remark ALGIANE: a mettre à plat : si p0-p1 est une ridge, on peut
 * reconstruire la mauvaise métrique non? Est-ce qu'il ne vaut mieux pas passer
 * la normale au triangle d'où l'on vient en argument pour le choix de la
 * métrique à reconstruire si c'est possible(quand on vient de grad2metSurfreq)?
 *
 */
int MMG5_buildridmet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int np0,
                     double ux,double uy,double uz,double mr[6],
                     double r[3][3] ) {
  MMG5_pPoint  p0;
  MMG5_pxPoint go;
  double       ps1,ps2,*n1,*n2,*t,*m,dv,dn,u[3];

  p0 = &mesh->point[np0];
  if ( !(MG_GEO & p0->tag) )  return 0;
  m = &met->m[6*np0];
  t = &p0->n[0];

  /* Check that ridge tangent is not null */
  assert ( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] > 0. && "Null tangent");

  go = &mesh->xpoint[p0->xp];

  /* Decide between the two possible configurations */
  n1 = &go->n1[0];
  n2 = &go->n2[0];

  ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
  ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];

  if ( fabs(ps2)<fabs(ps1) ) {
    n1 = &go->n2[0];
    dv = m[2];
    dn = m[4];
  }
  else{
    dv = m[1];
    dn = m[3];
  }

  /* Check that choosed normal is not null */
  /* Remark: a null second normal along ridge point in mmg3d may be the
   * consequence of the computation of the same normal with opposite sign for n1
   * and n2 at a previously inserted ridge point. This append when the ridge
   * delimits a closed angle an we choose the wrong point normal in bezierCP. */
  assert ( n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2] > 0. && "Null normal");

  u[0] = n1[1]*t[2] - n1[2]*t[1];
  u[1] = n1[2]*t[0] - n1[0]*t[2];
  u[2] = n1[0]*t[1] - n1[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) =
   * diag(m[0],dv,dn). Now, compute the metric in the canonical basis. */
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

  mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1] + dn*r[0][2]*r[0][2];
  mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1] + dn*r[0][2]*r[1][2];
  mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1] + dn*r[0][2]*r[2][2];
  mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1] + dn*r[1][2]*r[1][2];
  mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1] + dn*r[1][2]*r[2][2];
  mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1] + dn*r[2][2]*r[2][2];
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param np0 index of edge's extremity.
 * \param nt normal direction at the ridge point.
 * \param mr computed metric tensor.
 * \param r basis in which the metric is diagonal
 *
 * \return 0 if fail, 1 if the metric is build with respect to n1, 2 if it is
 * build with respect to n2.
 *
 * Build metric tensor at ridge point \a p0, when the 'good' normal direction is
 * given by \a nt and store the basis vectors in \a r.
 *
 */
int MMG5_buildridmetnor(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int np0,double nt[3],
                        double mr[6],double r[3][3] ) {
  MMG5_pPoint  p0;
  MMG5_pxPoint go;
  double       ps1,ps2,*n1,*n2,*t,*m,dv,dn,u[3];
  int          ier = 0;

  p0 = &mesh->point[np0];
  if ( !(MG_GEO & p0->tag) )  return 0;
  m = &met->m[6*np0];
  t = &p0->n[0];
  go = &mesh->xpoint[p0->xp];

  /* Decide between the two possible configurations */
  n1 = &go->n1[0];
  n2 = &go->n2[0];

  ps1 = nt[0]*n1[0] + nt[1]*n1[1] + nt[2]*n1[2];
  ps2 = nt[0]*n2[0] + nt[1]*n2[1] + nt[2]*n2[2];

  if ( fabs(ps2) > fabs(ps1) ) {
    n1 = &go->n2[0];
    dv = m[2];
    dn = m[4];
    ier = 2;
  }
  else{
    dv = m[1];
    dn = m[3];
    ier = 1;
  }

  u[0] = n1[1]*t[2] - n1[2]*t[1];
  u[1] = n1[2]*t[0] - n1[0]*t[2];
  u[2] = n1[0]*t[1] - n1[1]*t[0];

  /* If u = n1 ^ t, matrix of the desired metric in (t,u,n1) = diag(m[0],dv,0).
     Now, compute the metric in the canonical basis.*/
  r[0][0] = t[0];  r[0][1] = u[0];  r[0][2] = n1[0];
  r[1][0] = t[1];  r[1][1] = u[1];  r[1][2] = n1[1];
  r[2][0] = t[2];  r[2][1] = u[2];  r[2][2] = n1[2];

  mr[0] = m[0]*r[0][0]*r[0][0] + dv*r[0][1]*r[0][1] + dn*r[0][2]*r[0][2];
  mr[1] = m[0]*r[0][0]*r[1][0] + dv*r[0][1]*r[1][1] + dn*r[0][2]*r[1][2];
  mr[2] = m[0]*r[0][0]*r[2][0] + dv*r[0][1]*r[2][1] + dn*r[0][2]*r[2][2];
  mr[3] = m[0]*r[1][0]*r[1][0] + dv*r[1][1]*r[1][1] + dn*r[1][2]*r[1][2];
  mr[4] = m[0]*r[1][0]*r[2][0] + dv*r[1][1]*r[2][1] + dn*r[1][2]*r[2][2];
  mr[5] = m[0]*r[2][0]*r[2][0] + dv*r[2][1]*r[2][1] + dn*r[2][2]*r[2][2];

  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m pointer toward a \f$(2x2)\f$ metric.
 * \param n pointer toward a \f$(2x2)\f$ metric.
 * \param mr computed \f$(2x2)\f$ metric.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute the intersected (2 x 2) metric from metrics \a m and \a n : take
 * simultaneous reduction, and proceed to truncation in sizes.
 *
 */
int MMG5_intersecmet22(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  double        det,imn[4],lambda[2],vp[2][2],dm[2],dn[2],d0,d1,ip[4];
  double        isqhmin,isqhmax;
  static int8_t mmgWarn0 = 0;
  int           order;

  isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < MMG5_EPS*MMG5_EPS ) {
    if ( !mmgWarn0 ) {
      fprintf(stderr,"\n  ## Warning: %s: null metric det : %E \n",__func__,det);
      mmgWarn0 = 1;
    }
    return 0;
  }
  det = 1.0 / det;

  imn[0] = det * ( m[2]*n[0] - m[1]*n[1]);
  imn[1] = det * ( m[2]*n[1] - m[1]*n[2]);
  imn[2] = det * (-m[1]*n[0] + m[0]*n[1]);
  imn[3] = det * (-m[1]*n[1] + m[0]*n[2]);

  /* Find eigenvalues of imn */
  order = MMG5_eigenv2d(0,imn,lambda,vp);

  if ( !order ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 failing"
              " simultaneous reduction.\n",__func__);
    }
    return 0;
  }

  /* First case : matrices m and n are homothetic : n = lambda0*m */
  if ( order == 2 ) {
    /* Diagonalize m and truncate eigenvalues : trimn, det, etc... are reused */
    if ( fabs(m[1]) < MMG5_EPS ) {
      dm[0]    = m[0];
      dm[1]    = m[2];
      vp[0][0] = 1;
      vp[0][1] = 0;
      vp[1][0] = 0;
      vp[1][1] = 1;
    }
    else {
      MMG5_eigensym(m,dm,vp);
    }
    /* Eigenvalues of the resulting matrix*/
    dn[0] = MG_MAX(dm[0],lambda[0]*dm[0]);
    dn[0] = MG_MIN(isqhmin,MG_MAX(isqhmax,dn[0]));
    dn[1] = MG_MAX(dm[1],lambda[0]*dm[1]);
    dn[1] = MG_MIN(isqhmin,MG_MAX(isqhmax,dn[1]));

    /* Intersected metric = P diag(d0,d1){^t}P, P = (vp0, vp1) stored in columns */
    mr[0] = dn[0]*vp[0][0]*vp[0][0] + dn[1]*vp[1][0]*vp[1][0];
    mr[1] = dn[0]*vp[0][0]*vp[0][1] + dn[1]*vp[1][0]*vp[1][1];
    mr[2] = dn[0]*vp[0][1]*vp[0][1] + dn[1]*vp[1][1]*vp[1][1];

    return 1;
  }

  /* Second case : both eigenvalues of imn are distinct ; theory says qf associated to m and n
     are diagonalizable in basis (vp0, vp1) - the coreduction basis */
  else if( order == 1 ) {

    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp[0][0]*vp[0][0] + 2.0*m[1]*vp[0][0]*vp[0][1] + m[2]*vp[0][1]*vp[0][1];
    dm[1] = m[0]*vp[1][0]*vp[1][0] + 2.0*m[1]*vp[1][0]*vp[1][1] + m[2]*vp[1][1]*vp[1][1];
    dn[0] = n[0]*vp[0][0]*vp[0][0] + 2.0*n[1]*vp[0][0]*vp[0][1] + n[2]*vp[0][1]*vp[0][1];
    dn[1] = n[0]*vp[1][0]*vp[1][0] + 2.0*n[1]*vp[1][0]*vp[1][1] + n[2]*vp[1][1]*vp[1][1];

    /* Diagonal values of the intersected metric */
    d0 = MG_MAX(dm[0],dn[0]);
    d0 = MG_MIN(isqhmin,MG_MAX(d0,isqhmax));

    d1 = MG_MAX(dm[1],dn[1]);
    d1 = MG_MIN(isqhmin,MG_MAX(d1,isqhmax));

    /* Intersected metric = tP^-1 diag(d0,d1)P^-1, P = (vp0, vp1) stored in columns */
    det = vp[0][0]*vp[1][1] - vp[0][1]*vp[1][0];
    if ( fabs(det) < MMG5_EPS )  return 0;
    det = 1.0 / det;

    ip[0] =  vp[1][1]*det;
    ip[1] = -vp[1][0]*det;
    ip[2] = -vp[0][1]*det;
    ip[3] =  vp[0][0]*det;

    mr[0] = d0*ip[0]*ip[0] + d1*ip[2]*ip[2];
    mr[1] = d0*ip[0]*ip[1] + d1*ip[2]*ip[3];
    mr[2] = d0*ip[1]*ip[1] + d1*ip[3]*ip[3];
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param m pointer toward a \f$(3x3)\f$ metric.
 * \param n pointer toward a \f$(3x3)\f$ metric.
 * \param mr computed \f$(3x3)\f$ metric.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute the intersected (3 x 3) metric from metrics \a m and \a n : take
 * simultaneous reduction, and proceed to truncation in sizes.
 *
 */
int MMG5_intersecmet33(MMG5_pMesh mesh, double *m,double *n,double *mr) {
  double  vp[3][3],dm[3],dn[3],d[3],ivp[3][3];
  double  isqhmin,isqhmax;
  int8_t  i;

  isqhmin  = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax  = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  /* Simultaneous reduction */
  if( !MMG5_simred3d(mesh,m,n,dm,dn,vp) )
    return 0;

  /* Diagonal values of the intersected metric */
  for( i = 0; i < 3; i++ ) {
    d[i] = MG_MAX(dm[i],dn[i]);
    d[i] = MG_MIN(isqhmin,MG_MAX(d[i],isqhmax));
  }

  /* Intersected metric = tP^-1 diag(d0,d1,d2)P^-1, P = (vp0, vp1,vp2) stored in
   * columns */

  /* Compute the inverse of the eigenvectors matrix */
  if( !MMG5_invmat33(vp,ivp) )
    return 0;

  /* Recompose matrix */
  MMG5_simredmat(3,mr,d,(double *)ivp);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Test the intersection of (2 x 2) metrics.
 *
 */
int MMG5_test_intersecmet22(MMG5_pMesh mesh) {
  double m[3] = { 508., -504,  502.}; /* Test matrix 1 */
  double n[3] = {4020.,-2020.,1020.}; /* Test matrix 2 */
  double intex[3] = {4500.,-2500.,1500.}; /* Exact intersection */
  double intnum[3],maxerr; /* Numerical quantities */

  /** Compute intersection m^{-1}n */
  if( !MMG5_intersecmet22(mesh,m,n,intnum) )
    return 0;

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(3,intex,intnum);
  if( maxerr > 1000.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error metric intersection: in function %s, line %d, max error %e\n",
      __func__,__LINE__,maxerr);
    return 0;
  }

  /** Iteratively Recompute intersection of the result with (m,n), changing the
   *  matrix to be inverted. */
  for( int8_t it = 0; it < 20; it++ ) {

    /** Compute the intersection with n, invert n */
    if( !MMG5_intersecmet22(mesh,n,intnum,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(3,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }


    /** Compute the intersection with n, invert intnum */
    if( !MMG5_intersecmet22(mesh,intnum,n,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(3,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }


    /** Compute the intersection with n, invert n */
    if( !MMG5_intersecmet22(mesh,m,intnum,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(3,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }


    /** Compute the intersection with m, invert intnum */
    if( !MMG5_intersecmet22(mesh,intnum,m,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(3,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }

  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Test the intersection of (3 x 3) metrics.
 *
 */
int MMG5_test_intersecmet33(MMG5_pMesh mesh) {
  double m[6] = {111./2.,-109./2.,  89./2.,111./2.,-91./2.,111./2.}; /* Test matrix 1 */
  double n[6] = {409./2.,-393./2.,-407./2.,409./2.,391./2.,409./2.}; /* Test matrix 2 */
  double intex[6] = {254.,-246.,-154.,254.,146.,254.}; /* Exact intersection */
  double intnum[6],maxerr; /* Numerical quantities */

  /** Compute intersection m^{-1}n */
  if( !MMG5_intersecmet33(mesh,m,n,intnum) )
    return 0;

  /* Check error in norm inf */
  maxerr = MMG5_test_mat_error(6,intex,intnum);
  if( maxerr > 1000.*MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error metric intersection: in function %s, line %d, max error %e\n",
      __func__,__LINE__,maxerr);
    return 0;
  }

  /** Iteratively Recompute intersection of the result with (m,n), changing the
   *  matrix to be inverted. */
  for( int8_t it = 0; it < 20; it++ ) {

    /** Compute the intersection with n, invert n */
    if( !MMG5_intersecmet33(mesh,n,intnum,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(6,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }


    /** Compute the intersection with n, invert intnum */
    if( !MMG5_intersecmet33(mesh,intnum,n,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(6,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }


    /** Compute the intersection with n, invert n */
    if( !MMG5_intersecmet33(mesh,m,intnum,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(6,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }


    /** Compute the intersection with m, invert intnum */
    if( !MMG5_intersecmet33(mesh,intnum,m,intnum) )
      return 0;

    /* Check error in norm inf */
    maxerr = MMG5_test_mat_error(6,intex,intnum);
    if( maxerr > 1.e5*MMG5_EPSOK ) {
      fprintf(stderr,"  ## Error metric re-intersection: in function %s, line %d, iteration %d, max error %e\n",
        __func__,__LINE__,it,maxerr);
      return 0;
    }

  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param np global index of vertex in which we intersect the metrics.
 * \param me physical metric at point \a np.
 * \param n normal or tangent at point np.
 * \return 0 if fail, 1 otherwise.
 *
 * Intersect the surface metric held in np (supported in tangent plane of \a np)
 * with 3*3 physical metric in \a me. For ridge points, this function fill the
 * \f$ p_0 \rightarrow m[3]\f$ and \f$ p_0 \rightarrow m[4]\f$ fields that contains respectively the
 * specific sizes in the \f$n_1\f$ and \f$n_2\f$ directions.
 *
 */
int MMG5_mmgIntextmet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int np,double me[6],
                       double n[3]) {
  MMG5_pPoint         p0;
  MMG5_pxPoint        go;
  double              hu,isqhmin,isqhmax,dd,alpha1,alpha2,alpha3,u[3];
  double              lambda[3],vp[3][3];
  double              *m,*n1,*n2,*t,r[3][3],mrot[6],mr[3],mtan[3],metan[3];
  int                 order;
  int8_t              i;
  static int8_t       mmgWarn=0, mmgWarn1=0, mmgWarn2=0;

  isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

  p0 = &mesh->point[np];
  m  = &met->m[6*np];

  /** Case of a singular point : since MMG5_defmetsin() computes an isotropic
   * metric, scale it by the physical metric while keeping it isotropic */
  if ( MG_SIN(p0->tag) || (p0->tag & MG_NOM) ) {

    /** 1. Compute minimum size in physical metric.
     * Note that, if we enforce an isotropic metric at singular points, we have
     *  2 direct choices:
     *   - use the maximal size to make the physical metric iso. In a first
     * guess it seems to be a good idea to limit the impact of singular points
     * over the mesh but in practice it leads to collapse too much. See the
     * following pictures of adaptation over a constant aniso metric with mmgs
     * without gradation and with gradation.
     * \image html mmgs-max-on-sin-nograd.png "Max siz on sin points, no gradation" width=30%
     * \image html mmgs-max-on-sin-grad.png "Max siz on sin points, with gradation" width=15%
     *   - use the minimal size to make the physical metric iso. Adaptation
     * over a constant aniso metric with mmgs without gradation
     * shows the creation of a small element at each corner.
     * The computed iso size correspond to the size imposed on the
     * top surface. With gradation, the gradation removes the anisotropy everywhere
     * so the metric at corners has no effect).
     * \image html mmgs-min-on-sin-nograd.png "Max siz on sin points, no gradation" width=30%
     * \image html mmgs-min-on-sin-grad.png "Max siz on sin points, with gradation" width=15%
     *
     */
    order = MMG5_eigenv3d(1,me,lambda,vp);
    if ( !order ) {
      if ( !mmgWarn ) {
        fprintf(stderr,"\n  ## Warning: %s: Unable to diagonalize at least"
                " 1 metric.\n",__func__);
        mmgWarn = 1;
      }
      return 0;
    }

    hu = 0.0;
    for( i = 0; i < 3; i++ ) {
      hu = MG_MAX(hu,lambda[i]);
    }

    /** 2. Truncate physical size with respect to mesh constraints */
    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);

    /** 3. Overwrite diagonal metric if physical metric prescribes smaller sizes */
    if( hu > m[0] )
      m[0] = m[3] = m[5] = hu;

  }
  /** Case of a ridge point : take sizes in 3 directions t,n1,u */
  else if ( p0->tag & MG_GEO ) {
    /* Size prescribed by metric me in direction t */
    t = n;
    hu = me[0]*t[0]*t[0] + me[3]*t[1]*t[1] + me[5]*t[2]*t[2] \
      + 2.0*me[1]*t[0]*t[1] + 2.0*me[2]*t[0]*t[2] + 2.0*me[4]*t[1]*t[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[0] = MG_MAX(m[0],hu);

    /** 1. Size prescribed by metric me in direction u1 = n1 ^ t */
    assert ( p0->xp );
    go = &mesh->xpoint[p0->xp];
    n1 = &go->n1[0];
    n2 = &go->n2[0];

    u[0] = n1[1]*t[2] - n1[2]*t[1];
    u[1] = n1[2]*t[0] - n1[0]*t[2];
    u[2] = n1[0]*t[1] - n1[1]*t[0];
    dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if ( dd > MMG5_EPSD ) {
      dd = 1.0 / sqrt(dd);

      u[0] *= dd;
      u[1] *= dd;
      u[2] *= dd;

      hu = me[0]*u[0]*u[0] + me[3]*u[1]*u[1] + me[5]*u[2]*u[2]          \
        + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

      hu = MG_MIN(isqhmin,hu);
      hu = MG_MAX(isqhmax,hu);
      m[1] = MG_MAX(m[1],hu);
    }
    /** 2. Size prescribed by metric me in direction u2 = n2 ^ t */
    u[0] = n2[1]*t[2] - n2[2]*t[1];
    u[1] = n2[2]*t[0] - n2[0]*t[2];
    u[2] = n2[0]*t[1] - n2[1]*t[0];
    dd = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    if ( dd > MMG5_EPSD ) {
      dd = 1.0 / sqrt(dd);

      u[0] *= dd;
      u[1] *= dd;
      u[2] *= dd;

      hu =     me[0]*u[0]*u[0] +     me[3]*u[1]*u[1] +     me[5]*u[2]*u[2] \
        + 2.0*me[1]*u[0]*u[1] + 2.0*me[2]*u[0]*u[2] + 2.0*me[4]*u[1]*u[2];

      hu = MG_MIN(isqhmin,hu);
      hu = MG_MAX(isqhmax,hu);
      m[2] = MG_MAX(m[2],hu);
    }

    /** 3. Size prescribed by metric me in direction n1 */
    hu = me[0]*n1[0]*n1[0] + me[3]*n1[1]*n1[1] + me[5]*n1[2]*n1[2]
      + 2.0*me[1]*n1[0]*n1[1] + 2.0*me[2]*n1[0]*n1[2] + 2.0*me[4]*n1[1]*n1[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[3] = hu;

    /** 4. Size prescribed by metric me in direction n2 */
    hu = me[0]*n2[0]*n2[0] + me[3]*n2[1]*n2[1] + me[5]*n2[2]*n2[2]
      + 2.0*me[1]*n2[0]*n2[1] + 2.0*me[2]*n2[0]*n2[2] + 2.0*me[4]*n2[1]*n2[2];

    hu = MG_MIN(isqhmin,hu);
    hu = MG_MAX(isqhmax,hu);
    m[4] = hu;
  }
  /** Case of a ref, or regular point : intersect metrics in tangent plane */
  else {
    MMG5_rotmatrix(n,r);

    /** 1. Expression of both metrics in tangent plane */
    MMG5_rmtr(r,m,mrot);
    mtan[0] = mrot[0];
    mtan[1] = mrot[1];
    mtan[2] = mrot[3];


    MMG5_rmtr(r,me,mrot);
    metan[0] = mrot[0];
    metan[1] = mrot[1];
    metan[2] = mrot[3];

    /** 2. Intersection of metrics in the tangent plane */
    if ( !MMG5_intersecmet22(mesh,mtan,metan,mr) ) {
      if ( !mmgWarn1 ) {
        fprintf(stderr,"\n  ## Warning: %s: impossible metric intersection:"
                " surfacic metric skipped.\n",__func__);
        mmgWarn1 = 1;
      }
      m[0] = me[0];
      m[1] = me[1];
      m[2] = me[2];
      m[3] = me[3];
      m[4] = me[4];
      m[5] = me[5];

      return 0;
    }

    /** 3. Back to the canonical basis of \mathbb{R}^3 : me = ^tR*mr*R : mtan and
     * metan are reused */
    mtan[0]  = mr[0]*r[0][0] + mr[1]*r[1][0];
    mtan[1]  = mr[0]*r[0][1] + mr[1]*r[1][1];
    mtan[2]  = mr[0]*r[0][2] + mr[1]*r[1][2];
    metan[0] = mr[1]*r[0][0] + mr[2]*r[1][0];
    metan[1] = mr[1]*r[0][1] + mr[2]*r[1][1];
    metan[2] = mr[1]*r[0][2] + mr[2]*r[1][2];

    alpha1 = r[2][0]*mrot[5];
    alpha2 = r[2][1]*mrot[5];
    alpha3 = r[2][2]*mrot[5];

    m[0] = r[0][0] * mtan[0] + r[1][0] * metan[0] + r[2][0]*alpha1;
    m[1] = r[0][0] * mtan[1] + r[1][0] * metan[1] + r[2][0]*alpha2;
    m[2] = r[0][0] * mtan[2] + r[1][0] * metan[2] + r[2][0]*alpha3;
    m[3] = r[0][1] * mtan[1] + r[1][1] * metan[1] + r[2][1]*alpha2;
    m[4] = r[0][1] * mtan[2] + r[1][1] * metan[2] + r[2][1]*alpha3;
    m[5] = r[0][2] * mtan[2] + r[1][2] * metan[2] + r[2][2]*alpha3;

    /** 4. Truncate the metric in the third direction (because me was not
     * truncated) */
    order = MMG5_eigenv3d(1,m,lambda,vp);
    if ( !order ) {
      if ( !mmgWarn ) {
        fprintf(stderr,"\n  ## Warning: %s: Unable to diagonalize at least"
                " 1 metric.\n",__func__);
        mmgWarn = 1;
      }
      return 0;
    }

    for (i=0; i<3; i++) {
      if( lambda[i]<=0) {
        if ( !mmgWarn2 ) {
          fprintf(stderr,"\n  ## Warning: %s: at least 1 wrong metric "
                  "(eigenvalues : %e %e %e): surfacic metric skipped.\n",
                  __func__,lambda[0],lambda[1],lambda[2]);
          mmgWarn2 = 1;
        }
        m[0] = me[0];
        m[1] = me[1];
        m[2] = me[2];
        m[3] = me[3];
        m[4] = me[4];
        m[5] = me[5];
        return 0;
      }
      lambda[i]=MG_MIN(isqhmin,lambda[i]);
      lambda[i]=MG_MAX(isqhmax,lambda[i]);
    }

    m[0] = vp[0][0]*vp[0][0]*lambda[0] + vp[1][0]*vp[1][0]*lambda[1]
      + vp[2][0]*vp[2][0]*lambda[2];
    m[1] = vp[0][0]*vp[0][1]*lambda[0] + vp[1][0]*vp[1][1]*lambda[1]
      + vp[2][0]*vp[2][1]*lambda[2];
    m[2] = vp[0][0]*vp[0][2]*lambda[0] + vp[1][0]*vp[1][2]*lambda[1]
      + vp[2][0]*vp[2][2]*lambda[2];
    m[3] = vp[0][1]*vp[0][1]*lambda[0] + vp[1][1]*vp[1][1]*lambda[1]
      + vp[2][1]*vp[2][1]*lambda[2];
    m[4] = vp[0][1]*vp[0][2]*lambda[0] + vp[1][1]*vp[1][2]*lambda[1]
      + vp[2][1]*vp[2][2]*lambda[2];
    m[5] = vp[0][2]*vp[0][2]*lambda[0] + vp[1][2]*vp[1][2]*lambda[1]
      + vp[2][2]*vp[2][2]*lambda[2];
  }

  return 1;
}

/**
 * \param c0 table of the coordinates of the starting point.
 * \param n0 normal at the starting point.
 * \param m metric to be transported.
 * \param c1 table of the coordinates of the ending point.
 * \param n1 normal at the ending point.
 * \param mt computed metric.
 * \return 0 if fail, 1 otherwise.
 *
 * Parallel transport of a metric tensor field, attached to point \a c0, with
 * normal \a n0, to point \a c1, with normal \a n1.
 *
 */
int MMG5_paratmet(double c0[3],double n0[3],double m[6],double c1[3],double n1[3],double mt[6]) {
  double  r[3][3],mrot[6],mtan[3],lambda[2],vp[2][2],u[3],ps,ll;

  /* Take the induced metric tensor in the tangent plane by change of basis : R * M * {^t}R*/
  if ( !MMG5_rotmatrix(n0,r) )  return 0;
  MMG5_rmtr(r,m,mrot);
  mtan[0] = mrot[0];
  mtan[1] = mrot[1];
  mtan[2] = mrot[3];

  /* Take eigenvectors of metric tensor in tangent plane */
  MMG5_eigensym(mtan,lambda,vp);

  /* Eigenvector in canonical basis = {t}R*vp[0] */
  u[0] = r[0][0]*vp[0][0] + r[1][0]*vp[0][1];
  u[1] = r[0][1]*vp[0][0] + r[1][1]*vp[0][1];
  u[2] = r[0][2]*vp[0][0] + r[1][2]*vp[0][1];

  /* Projection in the tangent plane of c1 */
  ps = u[0]*n1[0] + u[1]*n1[1] + u[2]*n1[2];
  u[0] -= ps*n1[0];
  u[1] -= ps*n1[1];
  u[2] -= ps*n1[2];
  ll = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
  if ( ll < MMG5_EPSD )  return 0;
  ll = 1.0 / sqrt(ll);
  u[0] *= ll;
  u[1] *= ll;
  u[2] *= ll;

  /* And the transported metric is diag(lambda[0], lambda[1], mrot[5]) in basis
   * (u,n1^u,n1) */
  r[0][0] = u[0];
  r[1][0] = u[1];
  r[2][0] = u[2];

  r[0][1] = n1[1]*u[2] - n1[2]*u[1];
  r[1][1] = n1[2]*u[0] - n1[0]*u[2];
  r[2][1] = n1[0]*u[1] - n1[1]*u[0];

  ll = r[0][1]*r[0][1] + r[1][1]*r[1][1] + r[2][1]*r[2][1];
  if ( ll < MMG5_EPSD )  return 0;
  ll = 1.0 / sqrt(ll);
  r[0][1] *= ll;
  r[1][1] *= ll;
  r[2][1] *= ll;

  r[0][2] = n1[0];
  r[1][2] = n1[1];
  r[2][2] = n1[2];

  /*mt = R * diag(lambda[0], lambda[1], mrot[5])*{^t}R */
  mt[0] = lambda[0]*r[0][0]*r[0][0] + lambda[1]*r[0][1]*r[0][1]
    + mrot[5]*r[0][2]*r[0][2];

  mt[1] = lambda[0]*r[0][0]*r[1][0]
    + lambda[1]*r[0][1]*r[1][1] + mrot[5]*r[0][2]*r[1][2];

  mt[2] = lambda[0]*r[0][0]*r[2][0]
    + lambda[1]*r[0][1]*r[2][1] + mrot[5]*r[0][2]*r[2][2];

  mt[3] = lambda[0]*r[1][0]*r[1][0] + lambda[1]*r[1][1]*r[1][1]
    + mrot[5]*r[1][2]*r[1][2];

  mt[4] = lambda[0]*r[2][0]*r[1][0]
    + lambda[1]*r[2][1]*r[1][1] + mrot[5]*r[2][2]*r[1][2];

  mt[5] = lambda[0]*r[2][0]*r[2][0] + lambda[1]*r[2][1]*r[2][1]
    + mrot[5]*r[2][2]*r[2][2];

  return 1;
}
