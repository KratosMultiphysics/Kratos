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
 * \file mmg2d/solmap_2d.c
 * \brief  Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg2d.h"
#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 * \param ani 1 for aniso metric, 0 for iso one.
 *
 * \return 0 if fail, 1 if succeed.
 *
 * Truncate the metric computed by the DoSol function by hmax and hmin values
 * (if setted by the user). Set hmin and hmax if they are not setted.
 *
 */
static inline
int MMG2D_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met, int ani) {
  MMG5_pTria  ptt;
  MMG5_int    k;
  int         i,ier;

  assert ( mesh->info.optim || mesh->info.hsiz > 0. );

  /* Detect the points not used by triangles */
  ++mesh->base;
#ifndef NDEBUG
  for (k=1; k<=mesh->np; k++) {
    assert ( mesh->point[k].flag < mesh->base );
  }
#endif

  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) ) continue;

    for (i=0; i<3; i++) {
      mesh->point[ptt->v[i]].flag = mesh->base;
    }
  }

  /* Compute hmin/hmax on unflagged points and truncate the metric */
  if ( !ani ) {
    ier = MMG5_solTruncature_iso(mesh,met);
  }
  else {
    MMG5_solTruncature_ani = MMG5_2dSolTruncature_ani;
    ier = MMG5_solTruncature_ani(mesh,met);
  }

  return ier;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the sol structure
 * \return 1 if success
 *
 * Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 *
 */
int MMG2D_doSol_iso(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt,pt;
  MMG5_pPoint     p1,p2;
  double          ux,uy,dd;
  MMG5_int        k,ipa,ipb;
  int             i,ib;
  int             MMG_inxtt[5] = {0,1,2,0,1};

  // here we guess that we have less than int32 edges passing through a point
  int             *mark;
  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  /* Memory alloc */
  if ( sol->size!=1 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,sol->size);
    return 0;
  }

  if ( !MMG2D_Set_solSize(mesh,sol,MMG5_Vertex,mesh->np,sol->size) )
    return 0;

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] )  continue;

    for (i=0; i<3; i++) {
      ib  = MMG_inxtt[i+1];
      ipa = ptt->v[i];
      ipb = ptt->v[ib];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      ux  = p1->c[0] - p2->c[0];
      uy  = p1->c[1] - p2->c[1];
      dd  = sqrt(ux*ux + uy*uy);

      sol->m[ipa] += dd;
      mark[ipa]++;
      sol->m[ipb] += dd;
      mark[ipb]++;
    }
  }

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    if ( !mark[k] )  {
      continue;
    }
    sol->m[k] = sol->m[k] / (double)mark[k];
  }
  MMG5_SAFE_FREE(mark);

  /* Computation of hmin/hmax if not provided + size truncature */
  MMG2D_solTruncatureForOptim(mesh,sol,0);

  /* compute quality */
  if ( MMG2D_caltri ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      pt->qual = MMG2D_caltri_iso(mesh,sol,pt);
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the sol structure
 * \return 1 if success
 *
 * Compute anisotropic unit size map using statistical concept of
 * length distribution tensors (formula 5 of \cite COUPEZ20112391).
 *
 */
int MMG2D_doSol_ani(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt,pt;
  MMG5_pPoint     p1,p2;
  double          ux,uy,dd,tensordot[3];
  MMG5_int        k,ipa,ipb,iadr;
  int             i,ib;
  int             MMG_inxtt[5] = {0,1,2,0,1};

  // here we guess that we have less than int32 edges passing through a point
  int             *mark;
  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  /* Memory alloc */
  if ( sol->size!=3 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
              __func__,sol->size);
      return 0;
    }

  if ( !MMG2D_Set_solSize(mesh,sol,MMG5_Vertex,mesh->np,sol->size) )
    return 0;

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] )  continue;

    for (i=0; i<3; i++) {
      ib  = MMG_inxtt[i+1];
      ipa = ptt->v[i];
      ipb = ptt->v[ib];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      ux  = p1->c[0] - p2->c[0];
      uy  = p1->c[1] - p2->c[1];

      tensordot[0] = ux*ux;
      tensordot[1] = ux*uy;
      tensordot[2] = uy*uy;

      iadr = 3*ipa;
      sol->m[iadr]   += tensordot[0];
      sol->m[iadr+1] += tensordot[1];
      sol->m[iadr+2] += tensordot[2];
      mark[ipa]++;

      iadr = 3*ipb;
      sol->m[iadr]   += tensordot[0];
      sol->m[iadr+1] += tensordot[1];
      sol->m[iadr+2] += tensordot[2];
      mark[ipb]++;
    }
  }

  /* Compute metric tensor and hmax if not specified */
  for (k=1; k<=mesh->np; k++) {
    if ( !mark[k] ) {
      continue;
    }

    iadr = 3*k;
    /* Metric = nedges/dim * inv (sum(tensor_dot(edges,edges))).
     * sum(tensor_dot) is stored in sol->m so reuse tensordot to
     * compute M.  */
    dd = 1./(sol->m[iadr]*sol->m[iadr+2] - sol->m[iadr+1]*sol->m[iadr+1]);
    dd *= (double)mark[k]*0.5;

    tensordot[0] = sol->m[iadr+2];
    tensordot[1] = -sol->m[iadr+1];
    tensordot[2] = sol->m[iadr];

    sol->m[iadr]   = dd*tensordot[0];
    sol->m[iadr+1] = dd*tensordot[1];
    sol->m[iadr+2] = dd*tensordot[2];

#ifndef NDEBUG
    /* Check metric */
    double lambda[2],vp[2][2];
    MMG5_eigensym(sol->m+iadr,lambda,vp);

    assert (lambda[0] > 0. && lambda[1] > 0. && "Negative eigenvalue");

    /* Normally the case where the point belongs to only 2 colinear points is
    impossible */
    assert (isfinite(lambda[0]) && isfinite(lambda[1]) && "wrong eigenvalue");
#endif
  }
  MMG5_SAFE_FREE(mark);

  /* Size truncature */
  MMG2D_solTruncatureForOptim(mesh,sol,1);

  /* compute quality */
  if ( MMG2D_caltri ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      pt->qual = MMG2D_caltri_ani(mesh,sol,pt);
    }
  }

  return 1;
}
