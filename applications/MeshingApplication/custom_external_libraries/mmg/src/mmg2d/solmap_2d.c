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

#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the sol structure
 * \return 1 if success
 *
 * Compute isotropic size map according to the mean of the length of the edges
 * passing through a point.
 *
 */
int MMG2D_doSol(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      ptt,pt;
  MMG5_pPoint     p1,p2;
  double          ux,uy,dd;
  int             i,k,ib,iadr,ipa,ipb;
  int             MMG_inxtt[5] = {0,1,2,0,1};

  /* Memory alloc */
  if ( sol->size!=1 && sol->size!=3 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,sol->size);
    return 0;
  }

  if ( !MMG2D_Set_solSize(mesh,sol,MMG5_Vertex,mesh->np,sol->size) )
    return 0;

  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    p1->tagdel = 0;
  }
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !ptt->v[0] )  continue;

    if ( sol->size == 1 ) {
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
        p1->tagdel++;
        sol->m[ipb] += dd;
        p2->tagdel++;
      }
    }
    else if ( sol->size == 3 ) {
      for (i=0; i<3; i++) {
        ib  = MMG_inxtt[i+1];
        ipa = ptt->v[i];
        ipb = ptt->v[ib];
        p1  = &mesh->point[ipa];
        p2  = &mesh->point[ipb];

        ux  = p1->c[0] - p2->c[0];
        uy  = p1->c[1] - p2->c[1];
        dd  = sqrt(ux*ux + uy*uy);

        iadr = 3*ipa;
        sol->m[iadr]   += dd;
        p1->tagdel++;

        iadr = 3*ipb;
        sol->m[iadr]   += dd;
        p2->tagdel++;
      }
    }
    else return 0;
  }

  /* if hmax is not specified, compute it from the metric */
  if ( mesh->info.hmax < 0. ) {
    if ( sol->size == 1 ) {
      dd = 0.;
      for (k=1; k<=mesh->np; k++) {
        p1 = &mesh->point[k];
        if ( !p1->tagdel ) continue;
        dd = MG_MAX(dd,sol->m[k]);
      }
      assert ( dd );
    }
    else if ( sol->size == 3 ) {
      dd = FLT_MAX;
      for (k=1; k<=mesh->np; k++) {
        p1 = &mesh->point[k];
        if ( !p1->tagdel ) continue;
        iadr = 3*k;
        dd = MG_MIN(dd,sol->m[iadr]);
      }
      assert ( dd < FLT_MAX );
      dd = 1./sqrt(dd);
    }
    else {
      fprintf(stderr,"\n  # Error: %s: Unexpected solution size (%d)\n",
              __func__,sol->size);
      return 0;
    }
    mesh->info.hmax = 10.*dd;
  }

  /* vertex size */
  if ( sol->size == 1 ) {
    for (k=1; k<=mesh->np; k++) {
      p1 = &mesh->point[k];
      if ( !p1->tagdel )  {
        sol->m[k] = mesh->info.hmax;
        continue;
      }

      sol->m[k] = sol->m[k] / (double)p1->tagdel;
      p1->tagdel = 0;
    }
  }
  else if ( sol->size == 3 ) {
    for (k=1; k<=mesh->np; k++) {
      p1 = &mesh->point[k];
      iadr = 3*k;

      if ( !p1->tagdel )  {
        sol->m[iadr]   = 1./(mesh->info.hmax*mesh->info.hmax);
        sol->m[iadr+2] = sol->m[iadr];
        continue;
      }
      sol->m[iadr]   = (double)p1->tagdel*(double)p1->tagdel
        / (sol->m[iadr]*sol->m[iadr]);
      sol->m[iadr+2] = sol->m[iadr];
      p1->tagdel = 0;
    }
  }

  /* compute quality */
  if ( MMG2D_caltri ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      pt->qual = MMG2D_caltri(mesh,sol,pt);
    }
  }

  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"   HMAX %f\n",mesh->info.hmax);
  return 1;
}
