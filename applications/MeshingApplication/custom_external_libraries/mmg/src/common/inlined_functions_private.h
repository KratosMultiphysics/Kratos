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
 * \brief inlined Functions
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"

#ifndef _INLINED_FUNC_H
#define _INLINED_FUNC_H

/**
 * \param mesh pointer to the mesh structure.
 * \param np0 index of edge's extremity.
 * \param np1 index of edge's extremity.
 * \param m0 metric at point np0.
 * \param m1 metric at point np1.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of a curve edge according to the prescribed metric, 0 if fail.
 *
 * Compute the curve length of surface edge \f$[np0;np1]\f$ according to the
 * prescribed aniso metrics \a m0 and \a m1.
 *
 * \remark the edge has to be a boundary edge
 */
static inline
double MMG5_lenEdg(MMG5_pMesh mesh,MMG5_int np0,MMG5_int np1,
                    double *m0,double *m1,int8_t isedg) {
  MMG5_pPoint   p0,p1;
  double        gammaprim0[3],gammaprim1[3],t[3],*n1,*n2,ux,uy,uz,ps1,ps2,l0,l1;
  static int8_t mmgWarn=0;

  p0 = &mesh->point[np0];
  p1 = &mesh->point[np1];

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  /* computation of the two tangent vectors to the underlying curve of [i0i1] */
  if ( MG_SIN(p0->tag) || (MG_NOM & p0->tag) ) {
    gammaprim0[0] = ux;
    gammaprim0[1] = uy;
    gammaprim0[2] = uz;
  }
  else if ( isedg ) {
    memcpy(t,p0->n,3*sizeof(double));
    ps1 = ux*t[0] + uy*t[1] + uz*t[2];
    gammaprim0[0] = ps1*t[0];
    gammaprim0[1] = ps1*t[1];
    gammaprim0[2] = ps1*t[2];
  }
  else {
    if ( MG_GEO & p0->tag ) {
      //assert(p0->xp);
      n1 = &mesh->xpoint[p0->xp].n1[0];
      n2 = &mesh->xpoint[p0->xp].n2[0];
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
      ps2 = ux*n2[0] + uy*n2[1] + uz*n2[2];
      if ( fabs(ps2) < fabs(ps1) ) {
        n1  = &mesh->xpoint[p0->xp].n2[0];
        ps1 = ps2;
      }
    }
    else if ( MG_REF & p0->tag || MG_BDY & p0->tag ) {
      // ( MG_BDY  & p0->tag ) => mmg3d
      n1  = &mesh->xpoint[p0->xp].n1[0];
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    }
    else {
      // we come from mmgs because in mmg3d the boundary points are tagged
      // MG_BDY.
      n1  = &(p0->n[0]);
      ps1 = ux*n1[0] + uy*n1[1] + uz*n1[2];
    }
    gammaprim0[0] = ux - ps1*n1[0];
    gammaprim0[1] = uy - ps1*n1[1];
    gammaprim0[2] = uz - ps1*n1[2];
  }

  if ( MG_SIN(p1->tag) || (MG_NOM & p1->tag) ) {
    gammaprim1[0] = -ux;
    gammaprim1[1] = -uy;
    gammaprim1[2] = -uz;
  }
  else if ( isedg ) {
    memcpy(t,p1->n,3*sizeof(double));
    ps1 = -ux*t[0] - uy*t[1] - uz*t[2];
    gammaprim1[0] = ps1*t[0];
    gammaprim1[1] = ps1*t[1];
    gammaprim1[2] = ps1*t[2];
  }
  else {
    if ( MG_GEO & p1->tag ) {
      n1 = &mesh->xpoint[p1->xp].n1[0];
      n2 = &mesh->xpoint[p1->xp].n2[0];
      ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
      ps2 = -ux*n2[0] - uy*n2[1] - uz*n2[2];

      if ( fabs(ps2) < fabs(ps1) ) {
        n1  = &mesh->xpoint[p1->xp].n2[0];
        ps1 = ps2;
      }
    }
    else if ( MG_REF & p1->tag || MG_BDY & p1->tag ) {
      // ( MG_BDY  & p1->tag ) => mmg3d )
      n1  = &mesh->xpoint[p1->xp].n1[0];
      ps1 = - ux*n1[0] - uy*n1[1] - uz*n1[2];
    }
    else {
      // we come from mmgs because in mmg3d the boundary points are tagged
      // MG_BDY.
      n1  = &(p1->n[0]);
      ps1 = -ux*n1[0] - uy*n1[1] - uz*n1[2];
    }
    gammaprim1[0] = - ux - ps1*n1[0];
    gammaprim1[1] = - uy - ps1*n1[1];
    gammaprim1[2] = - uz - ps1*n1[2];
  }

  /* computation of the length of the two tangent vectors in their respective
   * tangent plane */
  /* l_ab = int_a^b sqrt(m_ij d_t x_i(t) d_t x_j(t) ) : evaluated by a 2-point
   * quadrature method. */
  l0 = m0[0]*gammaprim0[0]*gammaprim0[0] + m0[3]*gammaprim0[1]*gammaprim0[1] \
    + m0[5]*gammaprim0[2]*gammaprim0[2] \
    + 2.0*m0[1]*gammaprim0[0]*gammaprim0[1]  + 2.0*m0[2]*gammaprim0[0]*gammaprim0[2] \
    + 2.0*m0[4]*gammaprim0[1]*gammaprim0[2];

  l1 = m1[0]*gammaprim1[0]*gammaprim1[0] + m1[3]*gammaprim1[1]*gammaprim1[1] \
    + m1[5]*gammaprim1[2]*gammaprim1[2] \
    +2.0*m1[1]*gammaprim1[0]*gammaprim1[1]  + 2.0*m1[2]*gammaprim1[0]*gammaprim1[2] \
    + 2.0*m1[4]*gammaprim1[1]*gammaprim1[2];

  if( l0 < 0.) {
    if ( !mmgWarn ) {
      fprintf(stderr,"  ## Warning: %s: at least 1 negative edge length "
              "(%e)\n",__func__,l0);
      mmgWarn = 1;
    }
    return 0.;
  }
  if(l1 < 0.) {
    if ( !mmgWarn ) {
      fprintf(stderr,"  ## Warning: %s: at least 1 negative edge length "
              "(%e)\n",__func__,l1);
      mmgWarn = 1;
    }
    return 0.;
  }
  l0 = 0.5*(sqrt(l0) + sqrt(l1));

  return l0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param np0 index of edge's extremity.
 * \param np1 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of edge according to the prescribed metric, 0 if fail.
 *
 * Compute the curve length of surface edge \f$[i0;i1]\f$ according to the
 * prescribed aniso metric (for special storage of metrics at ridges
 * points). Here the length is computed taking into account the curve nature of
 * the surface edge.
 *
 * \remark the edge has to be a boundary edge
 */
static inline
double MMG5_lenSurfEdg_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int np0,MMG5_int np1,int8_t isedg) {
  MMG5_pPoint   p0,p1;
  double        *m0,*m1,met0[6],met1[6],ux,uy,uz,rbasis[3][3];
  static int8_t mmgWarn = 0;

  p0 = &mesh->point[np0];
  p1 = &mesh->point[np1];

  ux = p1->c[0] - p0->c[0];
  uy = p1->c[1] - p0->c[1];
  uz = p1->c[2] - p0->c[2];

  /* Set metrics */
  if ( MG_SIN(p0->tag) || (MG_NOM & p0->tag)) {
    m0 = &met->m[6*np0];
  }
  else if ( MG_GEO & p0->tag ) {
    /* Note that rbasis isn't used here */
    if ( !MMG5_buildridmet(mesh,met,np0,ux,uy,uz,met0,rbasis) )  {
      if ( !mmgWarn ) {
        fprintf(stderr,"  ## Warning: %s: a- unable to compute at least 1 ridge"
                " metric.\n",__func__);
        mmgWarn = 1;
      }
      return 0.;
    }
    m0 = met0;
  }
  else {
    m0 = &met->m[6*np0];
  }

  if ( MG_SIN(p1->tag) || (MG_NOM & p1->tag)) {
    m1 = &met->m[6*np1];
  }
  else if ( MG_GEO & p1->tag ) {
    /* Note that rbasis isn't used here */
    if ( !MMG5_buildridmet(mesh,met,np1,ux,uy,uz,met1,rbasis) )  {
      if ( !mmgWarn ) {
        fprintf(stderr,"  ## Warning: %s: b- unable to compute at least 1 ridge"
                " metric.\n",__func__);
        mmgWarn = 1;
      }
      return 0.;
    }
    m1 = met1;
  }
  else {
    m1 = &met->m[6*np1];
  }

  return MMG5_lenEdg(mesh,np0,np1,m0,m1,isedg);
}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param np0 index of edge's extremity.
 * \param np1 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise.
 * \return length of edge according to the prescribed metric.
 *
 * Compute the curve length of surface edge \f$[i0;i1]\f$ according to the
 * prescribed aniso metric (for classic storage of metrics at ridges points).
 *
 * \remark the edge has to be a boundary edge
 */
static inline
double MMG5_lenSurfEdg33_ani(MMG5_pMesh mesh,MMG5_pSol met,
                              MMG5_int np0,MMG5_int np1,int8_t isedg) {
  double        *m0,*m1;

  /* Set metrics */
  m0 = &met->m[6*np0];
  m1 = &met->m[6*np1];

  return MMG5_lenEdg(mesh,np0,np1,m0,m1,isedg);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param ip1 index of edge's extremity.
 * \param ip2 index of edge's extremity.
 * \param isedg 1 if the edge is a ridge, 0 otherwise (dummy arg for
 * compatibility with \a lenedg_ani).
 * \return length of edge according to the prescribed metric.
 *
 * Compute the "straight" length of edge \f$[i0;i1]\f$ according to the
 * prescribed iso metric.
 *
 */
static
inline double MMG5_lenSurfEdg_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int ip1,MMG5_int ip2, int8_t isedg) {
  MMG5_pPoint   p1,p2;
  double        h1,h2,l,r,len;

  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];
  h1 = met->m[ip1];
  h2 = met->m[ip2];
  l = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]) \
    + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  l = sqrt(l);
  r = h2 / h1 - 1.0;
  len = fabs(r) < MMG5_EPS ? l / h1 : l / (h2-h1) * log1p(r);

  return len;
}

/**
 * \param dim matrix size.
 * \param m matrix array.
 * \param dm diagonal values array.
 * \param iv array of inverse coreduction basis.
 *
 * Recompose a symmetric matrix from its diagonalization on a simultaneous
 * reduction basis.
 * \warning Eigenvectors in Mmg are stored as matrix rows (the first dimension
 * of the double array spans the number of eigenvectors, the second dimension
 * spans the number of entries of each eigenvector). So the inverse (left
 * eigenvectors) is also stored with transposed indices.
 */
static inline
void MMG5_simredmat(int8_t dim,double *m,double *dm,double *iv) {
  int8_t i,j,k,ij;

  /* Storage of a matrix as a one-dimensional array: dim*(dim+1)/2 entries for
   * a symmetric matrix. */
  ij = 0;

  /* Loop on matrix rows */
  for( i = 0; i < dim; i++ ) {
    /* Loop on the upper-triangular part of the matrix */
    for( j = i; j < dim; j++ ) {
      /* Initialize matrix entry */
      m[ij] = 0.0;
      /* Compute matrix entry as the recomposition of diagonal values after
       * projection on the coreduction basis, using the inverse of the
       * transformation:
       *
       * M_{ij} = \sum_{k,l} V^{-1}_{ki} Lambda_{kl} V^{-1}_{lj} =
       *        = \sum_{k} lambda_{k} V^{-1}_{ki} V^{-1}_{kj}
       *
       * Since the inverse of the transformation is the inverse of an
       * eigenvectors matrix (which is stored in Mmg by columns, and not by
       * rows), the storage of the inverse matrix is also transposed and the
       * indices have to be exchanged when implementing the above formula. */
      for( k = 0; k < dim; k++ ) {
        m[ij] += dm[k]*iv[i*dim+k]*iv[j*dim+k];
      }
      /* Go to the next entry */
      ++ij;
    }
  }
  assert( ij == (dim+1)*dim/2 );
}

#endif
