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
 * \file common/scalem.c
 * \brief Scale and unscale mesh and solution.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \return 1 if success, 0 if fail (computed bounding box too small).
 *
 * Compute the mesh bounding box and fill the \a min, \a max and \a delta fields
 * of the \a MMG5_info structure.
 *
 */
int MMG5_boundingBox(MMG5_pMesh mesh) {
  MMG5_pPoint    ppt;
  MMG5_int       k,i;
  double         dd;

  /* compute bounding box */
  for (i=0; i<mesh->dim; i++) {
    mesh->info.min[i] =  DBL_MAX;
    mesh->info.max[i] = -DBL_MAX;
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for (i=0; i<mesh->dim; i++) {
      if ( ppt->c[i] > mesh->info.max[i] )  mesh->info.max[i] = ppt->c[i];
      if ( ppt->c[i] < mesh->info.min[i] )  mesh->info.min[i] = ppt->c[i];
    }
    ppt->tmp = 0;
  }
  mesh->info.delta = 0.0;
  for (i=0; i<mesh->dim; i++) {
    dd = mesh->info.max[i] - mesh->info.min[i];
    if ( dd > mesh->info.delta )  mesh->info.delta = dd;
  }
  if ( mesh->info.delta < MMG5_EPSD ) {
    fprintf(stderr,"\n  ## Error: %s: unable to scale mesh:"
            " Check that your mesh contains non-zero points and "
            "valid elements.\n",__func__);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sethmin 1 if hmin is setted by the user.
 * \param sethmax 1 if hmax is setted by the user.
 *
 * Check the compatibility between the automatically computed hmin/hmax values
 * and the user settings.
 *
 */
void MMG5_check_hminhmax(MMG5_pMesh mesh, int8_t sethmin, int8_t sethmax) {

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

/**
 * \param mesh pointer to the mesh structure.
 * \return 0 if fail, 1 if succeed
 *
 * Check that hmin (resp. hmax) is not user setted if it is negative.
 *
 */
int MMG5_check_setted_hminhmax(MMG5_pMesh mesh) {

  if ( mesh->info.hmin < 0 ) {
    if ( mesh->info.sethmin ) {
      fprintf(stderr,"\n  ## Error: %s: unexpected case (negative user setted"
              " hmin).\n",__func__);
      return 0;
    }
  }

  if ( mesh->info.hmax < 0 ) {
    if ( mesh->info.sethmax ) {
      fprintf(stderr,"\n  ## Error: %s: unexpected case (negative user setted"
              " hmax).\n",__func__);
      return 0;
    }
  }

  return 1;
}

/**
 * \param met pointer to metric.
 * \param ip pointer to global index of point on which metric has to be truncated
 * \param isqhmin inverse square of hmin (min edge size)
 * \param isqhmax inverse square of hmax (max edge size)
 *
 * \return 0 if fail, 1 if succeed
 *
 * Truncation of anisotropic metric at point \a ip by \a hmin and \a hmax.
 *
 */
int MMG5_truncate_met3d(MMG5_pSol met, MMG5_int ip, double isqhmin, double isqhmax) {
  double        v[3][3],lambda[3],*m;
  int           i;
  static int8_t mmgWarn = 0;

  m = &met->m[(MMG5_int)met->size*ip];

  if ( !MMG5_eigenv3d(1,m,lambda,v) ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: Unable to diagonalize at least"
              " 1 metric.\n",__func__);
      mmgWarn = 1;
    }
    return 0;
  }

  for (i=0; i<3; i++) {
    if(lambda[i]<=0) {
      if ( !mmgWarn ) {
        fprintf(stderr,"\n  ## Warning: %s: at least 1 wrong metric "
                "(eigenvalues : %e %e %e).\n",__func__,lambda[0],
                lambda[1],lambda[2]);
        mmgWarn = 1;
      }
      return 0;
    }
    lambda[i]=MG_MIN(isqhmin,lambda[i]);
    lambda[i]=MG_MAX(isqhmax,lambda[i]);
  }

  m[0] = v[0][0]*v[0][0]*lambda[0] + v[1][0]*v[1][0]*lambda[1]
    + v[2][0]*v[2][0]*lambda[2];
  m[1] = v[0][0]*v[0][1]*lambda[0] + v[1][0]*v[1][1]*lambda[1]
    + v[2][0]*v[2][1]*lambda[2];
  m[2] = v[0][0]*v[0][2]*lambda[0] + v[1][0]*v[1][2]*lambda[1]
    + v[2][0]*v[2][2]*lambda[2];
  m[3] = v[0][1]*v[0][1]*lambda[0] + v[1][1]*v[1][1]*lambda[1]
    + v[2][1]*v[2][1]*lambda[2];
  m[4] = v[0][1]*v[0][2]*lambda[0] + v[1][1]*v[1][2]*lambda[1]
    + v[2][1]*v[2][2]*lambda[2];
  m[5] = v[0][2]*v[0][2]*lambda[0] + v[1][2]*v[1][2]*lambda[1]
    + v[2][2]*v[2][2]*lambda[2];

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param dd scaling value.
 *
 * \return 1 if success, 0 if fail.
 *
 * Scale and truncate by hmin/hmax the scalar metric stored in met.  If
 * hmin/hmax are not provided by the user, it is automatically computed from the
 * metric.
 *
 */
int MMG5_scale_scalarMetric(MMG5_pMesh mesh, MMG5_pSol met, double dd) {
  MMG5_int      k;
  int           ier;
  static int8_t mmgWarn0 = 0;

  ++mesh->base;

#ifndef NDEBUG
  for (k=1; k<=mesh->np; k++) {
    assert ( mesh->point[k].flag < mesh->base );
  }
#endif

  for (k=1; k<=mesh->np; k++)  {

    if( !MG_VOK( &mesh->point[k] ) ) {
      continue;
    }

    /* Set point flag to base so MMG5_solTruncature function will use its
     * data to compute the hmin/hmax values if not provided by the user. */
    mesh->point[k].flag = mesh->base;

    /* Check the metric */
    if ( met->m[k] <= 0 ) {
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Error: %s: at least 1 wrong metric.\n",
                __func__);
        return 0;
      }
    }
    /* normalization */
    met->m[k] *= dd;
  }

  ier = MMG5_solTruncature_iso(mesh, met);

  return ier;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param dd scaling value.
 *
 * \return 1 if success, 0 if fail.
 *
 * Scale and truncate by hmin/hmax the scalar metric stored in met.  If
 * hmin/hmax are not provided by the user, it is automatically computed from the
 * metric.
 *
 */
int MMG5_scale_tensorMetric(MMG5_pMesh mesh, MMG5_pSol met, double dd) {
  MMG5_int    k,iadr;
  int         i,ier;

  dd = 1.0 / (dd*dd);

  ++mesh->base;

  for (k=1; k<=mesh->np; k++)  {

    if( !MG_VOK( &mesh->point[k] ) ) {
      continue;
    }

    /* Set point flag to base so MMG5_solTruncature function will use its
     * data to compute the hmin/hmax values if not provided by the user. */
    mesh->point[k].flag = mesh->base;

    iadr = k*met->size;
    for (i=0; i<met->size; i++) {
      met->m[iadr+i] *= dd;
    }
  }

  ier = MMG5_solTruncature_ani(mesh, met);

  return ier;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the solution structure.
 *
 * \return 0 if fail, 1 if succeed.
 *
 * Compute hmin and hmax values from unflagged points (if not setted by the
 * user), check the values and truncate the metric.
 *
 */
int MMG5_solTruncature_iso(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pPoint ppt;
  double      hmin,hmax;
  MMG5_int    k;

  /* Security check: if hmin (resp. hmax) is not setted, it means that sethmin
   * (resp. sethmax) is not setted too */
  if ( !MMG5_check_setted_hminhmax(mesh) ) {
    return 0;
  }

  /* If not provided by the user, compute hmin/hmax from the metric computed by
   * the DoSol function. */
  hmin = FLT_MAX;
  hmax = 0.;
  if ( (!mesh->info.sethmin) || (!mesh->info.sethmax) ) {
    for (k=1; k<=mesh->np; k++)  {
      ppt = &mesh->point[k];
      if ( (!MG_VOK(ppt)) || (ppt->flag < mesh->base) ) continue;
      hmin = MG_MIN(hmin,met->m[k]);
      hmax = MG_MAX(hmax,met->m[k]);
    }
  }

  if ( !mesh->info.sethmin ) {
    mesh->info.hmin = hmin;
  }

  if ( !mesh->info.sethmax ) {
    mesh->info.hmax = hmax;
  }

  MMG5_check_hminhmax(mesh, mesh->info.sethmin, mesh->info.sethmax);

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->flag<mesh->base ) {
      met->m[k] = mesh->info.hmax;
    }
    else {
      met->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,met->m[k]));
    }
  }

  if ( mesh->info.ddebug ) {
    /* print unscaled values for debug purpose */
    fprintf(stdout,"     After truncature computation:   hmin %lf (user setted %d)\n"
            "                                     hmax %lf (user setted %d)\n",
            mesh->info.delta * mesh->info.hmin,mesh->info.sethmin,
            mesh->info.delta * mesh->info.hmax,mesh->info.sethmax);
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the solution structure.
 *
 * \return 0 if fail, 1 if succeed.
 *
 * Compute hmin and hmax values from unflagged points (if not setted by the
 * user), check the values and truncate the 2D metric.
 *
 */
int MMG5_2dSolTruncature_ani(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pPoint ppt;
  MMG5_int    k,iadr;
  double      isqhmin, isqhmax;
  double      lambda[2],vp[2][2];

 /* Security check: if hmin (resp. hmax) is not setted, it means that sethmin
   * (resp. sethmax) is not setted too */
  if ( !MMG5_check_setted_hminhmax(mesh) ) {
    return 0;
  }

  /* If not provided by the user, compute hmin/hmax from the metric computed by
   * the DoSol function. */
  isqhmin = 0.;
  isqhmax = FLT_MAX;
  if ( (!mesh->info.sethmin) || (!mesh->info.sethmax) ) {
    for (k=1; k<=mesh->np; k++)  {
      ppt = &mesh->point[k];
      if ( (!MG_VOK(ppt)) || (ppt->flag < mesh->base) ) continue;
      iadr = met->size*k;

      MMG5_eigensym(met->m+iadr,lambda,vp);

      assert (lambda[0] > 0. && lambda[1] > 0. && "Negative eigenvalue");

      isqhmin = MG_MAX(isqhmin,lambda[0]);
      isqhmin = MG_MAX(isqhmin,lambda[1]);

      isqhmax = MG_MIN(isqhmax,lambda[0]);
      isqhmax = MG_MIN(isqhmax,lambda[1]);
    }
  }

  if ( !mesh->info.sethmin ) {
    mesh->info.hmin = 1./sqrt(isqhmin);
  }

  if ( !mesh->info.sethmax ) {
    mesh->info.hmax = 1./sqrt(isqhmax);
  }

  MMG5_check_hminhmax(mesh, mesh->info.sethmin, mesh->info.sethmax);

  /* vertex size */
  isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

  for (k=1; k<=mesh->np; k++) {
    iadr = 3*k;

    ppt = &mesh->point[k];
    if ( ppt->flag < mesh->base ) {
      met->m[iadr] = met->m[iadr+2] = isqhmax;
      met->m[iadr+1] = 0;
    }
    else {
      MMG5_eigensym(met->m+iadr,lambda,vp);

      lambda[0]=MG_MAX(isqhmax,MG_MIN(isqhmin,lambda[0]));
      lambda[1]=MG_MAX(isqhmax,MG_MIN(isqhmin,lambda[1]));

      met->m[iadr]   = vp[0][0]*vp[0][0]*lambda[0] + vp[1][0]*vp[1][0]*lambda[1];
      met->m[iadr+1] = vp[0][0]*vp[0][1]*lambda[0] + vp[1][0]*vp[1][1]*lambda[1];
      met->m[iadr+2] = vp[0][1]*vp[0][1]*lambda[0] + vp[1][1]*vp[1][1]*lambda[1];
    }
  }

  if ( mesh->info.ddebug ) {
    /* print unscaled values for debug purpose */
    fprintf(stdout,"     After truncature computation:   hmin %lf (user setted %d)\n"
            "                                     hmax %lf (user setted %d)\n",
            mesh->info.delta * mesh->info.hmin,mesh->info.sethmin,
            mesh->info.delta * mesh->info.hmax,mesh->info.sethmax);
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the solution structure.
 *
 * \return 0 if fail, 1 if succeed.
 *
 * Compute hmin and hmax values from unflagged points (if not setted by the
 * user), check the values and truncate the 3D metric.
 *
 */
int MMG5_3dSolTruncature_ani(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_pPoint ppt;
  MMG5_int    k,iadr;
  double      isqhmin, isqhmax, isqhmaxcoe;
  double      lambda[3],vp[3][3];

  /* Security check: if hmin (resp. hmax) is not setted, it means that sethmin
   * (resp. sethmax) is not setted too */
  if ( !MMG5_check_setted_hminhmax(mesh) ) {
    return 0;
  }

  /* If not provided by the user, compute hmin/hmax from the metric computed by
   * the DoSol function. */
  isqhmin = 0.;
  isqhmax = FLT_MAX;
  isqhmaxcoe = 1./(MMG5_HMAXCOE*MMG5_HMAXCOE);
  if ( (!mesh->info.sethmin) || (!mesh->info.sethmax) ) {
    for (k=1; k<=mesh->np; k++)  {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || (ppt->flag < mesh->base) ) {
        continue;
      }
      iadr = met->size*k;

      /* Check metric */
      if (!MMG5_eigenv3d(1,met->m+iadr,lambda,vp) ) {
        fprintf(stdout, " ## Warning: %s: %d: non diagonalizable metric."
                " Impose hmax size at point\n",__func__,__LINE__);
        met->m[iadr+0] = FLT_MIN;
        met->m[iadr+1] = 0;
        met->m[iadr+2] = 0;
        met->m[iadr+3] = FLT_MIN;
        met->m[iadr+4] = 0;
        met->m[iadr+5] = FLT_MIN;
        continue;
      }
      assert ( lambda[0] > 0. && lambda[1] > 0.  && lambda[2] > 0.
               && "Negative eigenvalue");

      /* Take only meaningful values into account: i.e. values at used points
       * and eigenvalues not equal to MMG5_HMAXCOE for Mmgs (value used when
       * length in normal direction can't be computed.) */
      int j;
      for ( j=0; j<3; ++j ) {
        if ( isfinite(lambda[j]) && fabs(lambda[j]-isqhmaxcoe) > MMG5_EPS ) {
          isqhmax = MG_MIN(isqhmax,lambda[j]);
          isqhmin = MG_MAX(isqhmin,lambda[j]);
        }
      }
    }
  }

  if ( !mesh->info.sethmin ) {
    mesh->info.hmin = 1./sqrt(isqhmin);
  }

  if ( !mesh->info.sethmax ) {
    mesh->info.hmax = 1./sqrt(isqhmax);
  }

  /* Check the compatibility between the user settings and the automatically
   * computed values */
  MMG5_check_hminhmax(mesh,mesh->info.sethmin,mesh->info.sethmax);

  /* vertex size */
  isqhmin = 1./(mesh->info.hmin*mesh->info.hmin);
  isqhmax = 1./(mesh->info.hmax*mesh->info.hmax);

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;

    if ( (ppt->flag < mesh->base) || !MMG5_truncate_met3d(met,k,isqhmin,isqhmax) ) {
      /* Fail to diagonalize metric: put hmax */
      iadr = 6*k;
      met->m[iadr]   = isqhmax;
      met->m[iadr+1] = 0.;
      met->m[iadr+2] = 0.;
      met->m[iadr+3] = isqhmax;
      met->m[iadr+4] = 0.;
      met->m[iadr+5] = isqhmax;
  }
  }

  if ( mesh->info.ddebug ) {
    /* print unscaled values for debug purpose */
    fprintf(stdout,"     After truncature computation:   hmin %lf (user setted %d)\n"
            "                                     hmax %lf (user setted %d)\n",
            mesh->info.delta * mesh->info.hmin,mesh->info.sethmin,
            mesh->info.delta * mesh->info.hmax,mesh->info.sethmax);
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to a metric
 * \param sol pointer to a solution structure (level-set or displacement).
 * \param dd pointer to the scaling value (to fill)
 *
 * \return 1 if success, 0 if fail.
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 */
int MMG5_scale_meshAndSol(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,double *dd) {
  MMG5_pPoint    ppt;
  MMG5_pPar      par;
  MMG5_int       k,i;
  int8_t         hsizOrOptim;

  /* sol is a level-set or a displacement so it cannot be an aniso metric */
  if ( sol ) { assert ( sol->type == MMG5_Scalar || sol->type == MMG5_Vector ); }

  /* met is a metric so it cannot be a vector */
  if ( met ) { assert ( met->type != MMG5_Vector ); }

  /* if we are not in iso, isosurf or lagrangian mode, check that sol isn't allocated */
  if ( (!mesh->info.iso) && (!mesh->info.isosurf)  && mesh->info.lag < 0 ) {
    assert ( !(sol && sol->m) );
  }

  /* compute bounding box */
  if ( ! MMG5_boundingBox(mesh) ) return 0;

  /* normalize coordinates */
  *dd = 1.0 / mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for (i=0 ; i<mesh->dim ; i++)
      ppt->c[i] = (*dd) * (ppt->c[i] - mesh->info.min[i]);
  }

  mesh->info.hausd *= (*dd);
  mesh->info.ls    *= (*dd);
  mesh->info.hsiz  *= (*dd);

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hmin  *= (*dd);
    par->hmax  *= (*dd);
    par->hausd *= (*dd);
  }

  /* Security check: if hmin (resp. hmax) is not setted, it means that sethmin
   * (resp. sethmax) is not setted too */
  if ( !MMG5_check_setted_hminhmax(mesh) ) {
    return 0;
  }

  if ( mesh->info.sethmin ) {
    mesh->info.hmin  *= (*dd);
  }
  if ( mesh->info.sethmax ) {
    mesh->info.hmax  *= (*dd);
  }

  hsizOrOptim = ( mesh->info.hsiz > 0. || mesh->info.optim )? 1 : 0;

  /* if we are not in optim or hsiz mode and if we don't have a metric, compute
   * default hmin/hmax */
  if ( (!hsizOrOptim) && (!(met && met->np)) ) {
    /* Set default values to hmin/hmax from the bounding box if not provided by
     * the user */
    if ( !MMG5_Set_defaultTruncatureSizes(mesh,mesh->info.sethmin,mesh->info.sethmax) ) {
      fprintf(stderr,"\n  ## Error: %s: Exit program.\n",__func__);
      return 0;
    }
  }

  if ( sol && sol->np ) {
    for ( k=sol->size; k<sol->size*(mesh->np+1); k++ ) {
      sol->m[k]   *= (*dd);
    }
  }

  return 1;

}
/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param sol pointer to a solution structure (level-set or displacement).
 *
 * \return 1 if success, 0 if fail (computed bounding box too small
 * or one af the anisotropic input metric is not valid).
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 */
int MMG5_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  double  dd;

  if ( !MMG5_scale_meshAndSol(mesh,met,sol,&dd) ) {
    return 0;
  }

  if ( (!met) || (met && !met->np) || (!met->m) ) {
    return 1;
  }

  if ( met->size == 1 ) {
    /* isotropic metric */
    if ( !MMG5_scale_scalarMetric ( mesh, met, dd ) ) {
      return 0;
    }
        }
  else if ( met->size == (mesh->dim-1)*3 ) {
    /* anisotropic metric: met->size=3 in 2D and 6 in 3D */
    if ( !MMG5_scale_tensorMetric ( mesh, met, dd ) ) {
          return 0;
        }
          }
  else {
    fprintf(stderr,"\n  ## Error: %s: unexpected metric size (%d)\n",__func__,met->size);
   }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to a metric.
 * \param sol pointer to a solution structure (level-set or displacement).
 *
 * \return 1.
 *
 * Unscale the mesh and the size informations to their initial sizes.
 *
 */
int MMG5_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  MMG5_pPoint     ppt;
  double          dd;
  int             i;
  MMG5_int        k,iadr;
  MMG5_pPar       par;

  /* de-normalize coordinates */
  dd = mesh->info.delta;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    for ( i=0; i<mesh->dim; ++i ) {
      ppt->c[i] = ppt->c[i] * dd + mesh->info.min[i];
    }
  }

  /* unscale paramter values */
  if ( !mesh->info.sethmin ) {
    mesh->info.hmin = MMG5_NONSET_HMIN;
  }
  else {
    mesh->info.hmin  *= dd;
  }

  if ( !mesh->info.sethmax ) {
    mesh->info.hmax = MMG5_NONSET_HMAX;
  }
  else {
    mesh->info.hmax  *= dd;
  }
  mesh->info.hausd *= dd;
  mesh->info.ls    *= dd;
  mesh->info.hsiz  *= dd;

  /* normalize local parameters */
  for (k=0; k<mesh->info.npar; k++) {
    par = &mesh->info.par[k];
    par->hmin  *= dd;
    par->hmax  *= dd;
    par->hausd *= dd;
  }

  /* de-normalize level-set or displacement */
  if ( sol && sol->np && sol->m ) {
    assert ( sol->size != MMG5_Tensor );
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      iadr = k*sol->size;
      for (i=0; i<sol->size; i++)  sol->m[iadr+i] *= dd;
    }
  }

  /* reset the scaling data to ensure that if we try to unscale again, we will
   * do nothing */
  mesh->info.delta = 1.;
  mesh->info.min[0]= 0.;
  mesh->info.min[1]= 0.;
  mesh->info.min[2]= 0.;

  /* de-normalize metric */
  if ( !(met && met->np && met->m) )  return 1;

  /* unscale sizes */
  switch (met->type) {
  case 1:
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      met->m[k] *= dd;
    }
    break;
  case 3:
    dd = 1.0 / (dd*dd);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      for (i=0; i<met->size; i++)  met->m[met->size*k+i] *= dd;
    }
    break;
  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected metric size (%d)\n",__func__,met->size);
    break;
  }

  return 1;
}
