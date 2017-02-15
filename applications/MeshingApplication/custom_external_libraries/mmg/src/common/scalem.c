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
 * \file common/scalem.c
 * \brief Scale and unscale mesh and solution.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 if fail (computed bounding box too small).
 *
 * Compute the mesh bounding box and fill the \a min, \a max and \a delta fields
 * of the \a _MMG5_info structure.
 *
 */
int _MMG5_boundingBox(MMG5_pMesh mesh) {
  MMG5_pPoint    ppt;
  int            k,i;
  double         dd;

  /* compute bounding box */
  for (i=0; i<3; i++) {
    mesh->info.min[i] =  DBL_MAX;
    mesh->info.max[i] = -DBL_MAX;
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for (i=0; i<3; i++) {
      if ( ppt->c[i] > mesh->info.max[i] )  mesh->info.max[i] = ppt->c[i];
      if ( ppt->c[i] < mesh->info.min[i] )  mesh->info.min[i] = ppt->c[i];
    }
    ppt->tmp = 0;
  }
  mesh->info.delta = 0.0;
  for (i=0; i<3; i++) {
    dd = mesh->info.max[i] - mesh->info.min[i];
    if ( dd > mesh->info.delta )  mesh->info.delta = dd;
  }
  if ( mesh->info.delta < _MMG5_EPSD ) {
    fprintf(stderr,"  ## Unable to scale mesh:\n");
    fprintf(stderr,"  ## Check that your mesh contains non-zero points and "
            "valid elements.\n");
    return(0);
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric or solution structure.
 * \return 1 if success, 0 if fail (computed bounding box too small
 * or one af the anisotropic input metric is not valid).
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 */
int _MMG5_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pPoint    ppt;
  double         dd,d1;
  int            k,i,sethmin,sethmax;
  MMG5_pPar      par;
  double         *m;
  double         lambda[3],v[3][3];


  /* compute bounding box */
  if ( ! _MMG5_boundingBox(mesh) ) return(0);

  /* normalize coordinates */
  dd = 1.0 / mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->c[0] = dd * (ppt->c[0] - mesh->info.min[0]);
    ppt->c[1] = dd * (ppt->c[1] - mesh->info.min[1]);
    ppt->c[2] = dd * (ppt->c[2] - mesh->info.min[2]);
  }

  mesh->info.hausd *= dd;
  mesh->info.ls    *= dd;

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hmin  *= dd;
    par->hmax  *= dd;
    par->hausd *= dd;
  }

  /* Check if hmin/hmax have been provided by the user and scale it if yes */
  sethmin = 0;
  sethmax = 0;

  if ( mesh->info.hmin > 0. ) {
    mesh->info.hmin  *= dd;
    sethmin = 1;
  }
  if ( mesh->info.hmax > 0. ) {
    mesh->info.hmax  *= dd;
    sethmax = 1;
  }

  /* Warning: we don't want to compute hmin/hmax from the level-set or the
   * displacement! */
  if ( mesh->info.iso || (mesh->info.lag>-1) || (!met->m && !mesh->info.optim) ) {
    /* Set default values to hmin/hmax from the bounding box if not provided by
     * the user */
    if ( !sethmin )  mesh->info.hmin  = 0.01;

    if ( !sethmax )  mesh->info.hmax  = 2.;

    if ( mesh->info.hmax < mesh->info.hmin ) {
      if ( sethmin && sethmax ) {
        fprintf(stderr,"  ## Error: mismatch parameters:"
                " minimal mesh size larger than maximal one.\n");
        fprintf(stderr,"  Exit program.\n");
        exit(EXIT_FAILURE);
      }
      else if ( sethmin )
        mesh->info.hmax = 100. * mesh->info.hmin;
      else
        mesh->info.hmin = 0.01 * mesh->info.hmax;
    }
    sethmin = 1;
    sethmax = 1;
  }


  /* normalize sizes and if not provided by user, compute hmin/hmax */
  if ( met->m ) {
    if ( met->size == 1 ) {
      for (k=1; k<=mesh->np; k++) {
        met->m[k] *= dd;
        /* Check the metric */
        if ( (!mesh->info.iso) && met->m[k] <= 0) {
          fprintf(stderr,"  ## ERROR: WRONG METRIC AT POINT %d -- \n",k);
          return(0);
        }
      }

      /* compute hmin and hmax parameters if not provided by the user */
      if ( !sethmin ) {
        mesh->info.hmin = FLT_MAX;
        for (k=1; k<=mesh->np; k++)  {
          mesh->info.hmin = MG_MIN(mesh->info.hmin,met->m[k]);
        }
      }
      if ( !sethmax ) {
        mesh->info.hmax = 0.;
        for (k=1; k<=mesh->np; k++)  {
          mesh->info.hmax = MG_MAX(mesh->info.hmax,met->m[k]);
        }
      }

    }
    else if ( met->size==3 ){
      for (k=1; k<=mesh->np; k++) {
        met->m[3*k]   *= dd;
        met->m[3*k+1] *= dd;
        met->m[3*k+2] *= dd;
      }
    } else { //met->size==6
      d1 = 1.0 / (dd*dd);
      for (k=1; k<mesh->np+1; k++) {
        for ( i=0; i<met->size; ++i ) {
          met->m[6*k+i] *= d1;
        }
      }

      /* compute hmin and hmax parameters if not provided by the user and check
       * the input metric */
      if ( !sethmin ) {
        mesh->info.hmin = FLT_MAX;
      }

      if ( !sethmax ) {
        mesh->info.hmax = 0.;
      }

      for (k=1; k<mesh->np+1; k++) {
        m    = &met->m[6*k];

        /* Check the input metric */
        if ( !_MMG5_eigenv(1,m,lambda,v) ) {
          fprintf(stderr,"  ## ERROR: WRONG METRIC AT POINT %d -- \n",k);
          return(0);
        }
        for (i=0; i<3; i++) {
          if(lambda[i]<=0) {
            fprintf(stderr,"  ## ERROR: WRONG METRIC AT POINT %d -- eigenvalue :"
                   " %e %e %e -- det %e\n",k,lambda[0],lambda[1],lambda[2],
                   m[0]*(m[3]*m[5]-m[4]*m[4])-m[1]*(m[1]*m[5]-m[2]*m[4])+
                   m[2]*(m[1]*m[4]-m[2]*m[3]));
            fprintf(stderr,"WRONG METRIC AT POINT %d -- metric %e %e %e %e %e %e\n",
                   k,m[0],m[1],m[2],m[3],m[4],m[5]);
            return(0);
          }
          if ( !sethmin )
            mesh->info.hmin = MG_MIN(mesh->info.hmin,1./sqrt(lambda[i]));
          if ( !sethmax )
            mesh->info.hmax = MG_MAX(mesh->info.hmax,1./sqrt(lambda[i]));
        }
      }
    }
    if ( !sethmin ) {
      mesh->info.hmin *=.1;
      /* Check that user has not given a hmax value lower that the founded
       * hmin. */
      if ( mesh->info.hmin > mesh->info.hmax ) {
        mesh->info.hmin = 0.1*mesh->info.hmax;
      }
    }
    if ( !sethmax ) {
      mesh->info.hmax *=10.;
      /* Check that user has not given a hmin value bigger that the founded
       * hmax. */
      if ( mesh->info.hmax < mesh->info.hmin ) {
        mesh->info.hmax = 10.*mesh->info.hmin;
      }
    }
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric or solution structure.
 * \return 1.
 *
 * Unscale the mesh and the size informations to their initial sizes.
 *
 */
int _MMG5_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pPoint     ppt;
  double          dd;
  int             k,i;
  MMG5_pPar       par;

  /* de-normalize coordinates */
  dd = mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->c[0] = ppt->c[0] * dd + mesh->info.min[0];
    ppt->c[1] = ppt->c[1] * dd + mesh->info.min[1];
    ppt->c[2] = ppt->c[2] * dd + mesh->info.min[2];
  }

  /* unscale paramter values */
  mesh->info.hmin  *= dd;
  mesh->info.hmax  *= dd;
  mesh->info.hausd *= dd;
  mesh->info.ls    *= dd;

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hmin  *= dd;
    par->hmax  *= dd;
    par->hausd *= dd;
  }

  /* unscale sizes */
  if ( met->m ) {
    if ( met->size == 6 ) {
      dd = 1.0 / (dd*dd);
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )  continue;
        for (i=0; i<6; i++)  met->m[6*k+i] *= dd;
      }
    }
    else {
      for (k=1; k<=mesh->np ; k++) {
        ppt = &mesh->point[k];
        if ( MG_VOK(ppt) )  met->m[k] *= dd;
      }
    }
  }

  return(1);
}
