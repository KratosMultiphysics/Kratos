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
 * \file mmg2d/simred_2d.c
 * \brief Simultaneous reduction function
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh
 * \param m first matrix
 * \param n second matrix
 * \param dm eigenvalues of m in the coreduction basis (to fill)
 * \param dn eigenvalues of n in the coreduction basis (to fill)
 * \param vp coreduction basis (to fill)
 *
 * \return 0 if fail 1 otherwise.
 *
 * Perform simultaneous reduction of matrices \a m and \a n.
 *
 */
int MMG2D_simred(MMG5_pMesh mesh,double *m,double *n,double dm[2],
                 double dn[2],double vp[2][2] ) {

  double       det,dd,sqDelta,trimn,vnorm,lambda[2],imn[4];
  static char  mmgWarn0=0;

  /* Compute imn = M^{-1}N */
  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < MMG5_EPS*MMG5_EPS ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 null metric det : %E \n",
              __func__,det);
    }
    return 0;
  }
  det = 1.0 / det;

  imn[0] = det * ( m[2]*n[0] - m[1]*n[1]);
  imn[1] = det * ( m[2]*n[1] - m[1]*n[2]);
  imn[2] = det * (-m[1]*n[0] + m[0]*n[1]);
  imn[3] = det * (-m[1]*n[1] + m[0]*n[2]);
  dd = imn[0] - imn[3];
  sqDelta = sqrt(fabs(dd*dd + 4.0*imn[1]*imn[2]));
  trimn = imn[0] + imn[3];

  lambda[0] = 0.5 * (trimn - sqDelta);
  if ( lambda[0] < 0.0 ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 metric with a "
              "negative eigenvalue: %f \n",__func__,lambda[0]);
    }
    return 0;
  }

  /* First case : matrices m and n are homothetic: n = lambda0*m */
  if ( sqDelta < MMG5_EPS ) {

    /* Subcase where m is diaonal */
    if ( fabs(m[1]) < MMG5_EPS ) {
      dm[0]   = m[0];
      dm[1]   = m[2];
      vp[0][0] = 1;
      vp[0][1] = 0;
      vp[1][0] = 0;
      vp[1][1] = 1;
    }
    /* Subcase where m is not diagonal; dd,trimn,... are reused */
    else
      MMG5_eigensym(m,dm,vp);

    /* Eigenvalues of metric n */
    dn[0] = lambda[0]*dm[0];
    dn[1] = lambda[0]*dm[1];

  }
  /* Second case: both eigenvalues of imn are distinct ; theory says qf associated to m and n
   are diagonalizable in basis (vp[0], vp[1]) - the coreduction basis */
  else {
    lambda[1] = 0.5 * (trimn + sqDelta);
    assert(lambda[1] >= 0.0);

    vp[0][0] = imn[1];
    vp[0][1] = (lambda[0] - imn[0]);
    vnorm  = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);

    if ( vnorm < MMG5_EPS ) {
      vp[0][0] = (lambda[0] - imn[3]);
      vp[0][1] = imn[2];
      vnorm  = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);
    }

    vnorm   = 1.0 / vnorm;
    vp[0][0] *= vnorm;
    vp[0][1] *= vnorm;

    vp[1][0] = imn[1];
    vp[1][1] = (lambda[1] - imn[0]);
    vnorm  = sqrt(vp[1][0]*vp[1][0] + vp[1][1]*vp[1][1]);

    if ( vnorm < MMG5_EPS ) {
      vp[1][0] = (lambda[1] - imn[3]);
      vp[1][1] = imn[2];
      vnorm  = sqrt(vp[1][0]*vp[1][0] + vp[1][1]*vp[1][1]);
    }

    vnorm   = 1.0 / vnorm;
    vp[1][0] *= vnorm;
    vp[1][1] *= vnorm;

    /* Compute diagonal values in simultaneous reduction basis */
    dm[0] = m[0]*vp[0][0]*vp[0][0] + 2.0*m[1]*vp[0][0]*vp[0][1] + m[2]*vp[0][1]*vp[0][1];
    dm[1] = m[0]*vp[1][0]*vp[1][0] + 2.0*m[1]*vp[1][0]*vp[1][1] + m[2]*vp[1][1]*vp[1][1];
    dn[0] = n[0]*vp[0][0]*vp[0][0] + 2.0*n[1]*vp[0][0]*vp[0][1] + n[2]*vp[0][1]*vp[0][1];
    dn[1] = n[0]*vp[1][0]*vp[1][0] + 2.0*n[1]*vp[1][0]*vp[1][1] + n[2]*vp[1][1]*vp[1][1];
  }

  assert ( dm[0] > MMG5_EPS && dn[0] > MMG5_EPS );
  assert ( dm[1] > MMG5_EPS && dn[1] > MMG5_EPS );

  return 1;
}
