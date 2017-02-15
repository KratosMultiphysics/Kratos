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
  int             i,k,ib,ipa,ipb;
  int             MMG_inxtt[5] = {0,1,2,0,1};

  sol->np = mesh->np;
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    p1->tagdel = 0;
  }
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
      p1->tagdel++;
      sol->m[ipb] += dd;
      p2->tagdel++;
    }
  }

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !p1->tagdel )  {
      sol->m[k] = FLT_MAX;
      continue;
    }

    sol->m[k] = sol->m[k] / (double)p1->tagdel;
    p1->tagdel = 0;
  }

  /* compute quality */
  if ( MMG2_caltri_in ) {
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      pt->qual = MMG2_caltri_in(mesh,sol,pt);
    }
  }

  if ( mesh->info.imprim < -4 )
    fprintf(stdout,"     HMIN %f   HMAX %f\n",mesh->info.hmin,mesh->info.hmax);
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 1 if success, 0 if strongly fail.
 *
 * Update the metric definition by taking into accounts the
 * curvature of the external and internal curves present in the mesh.
 *
 **/
int _MMG2D_defBdrySiz(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria       pt;
  MMG5_pPoint      p1,p2;
  double           t1[2],t2[2],b1[2],b2[2],gpp1[2],gpp2[2],pv,M1,M2;
  double           ps1,ps2,ux,uy,ll,li,lm,hmax,hausd,hmin;
  int              k,ip1,ip2,iadr,*adja,adj;
  unsigned char    i,i1,i2;


  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Update isotropic map at boundary edges\n");

  if ( mesh->info.hmax < 0.0 )  mesh->info.hmax = 0.5 * mesh->info.delta;

  hmax = mesh->info.hmax;
  hausd = mesh->info.hausd;
  hmin = mesh->info.hmin;

  /* Only the boundary edges impose a minimal size feature */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    iadr  = 3*(k-1) + 1;
    adja  = &mesh->adja[iadr];

    for (i=0; i<3; i++) {
      adj = adja[i]/3;

      if ( adj && pt->ref==mesh->tria[adj].ref ) continue;

      /* Mark edge as boundary */
      pt->tag[i] |= M_BDRY;

      i1 = _MMG5_inxt2[i];
      i2 = _MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];

      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      ll = ux*ux + uy*uy;
      if ( ll < _MMG5_EPSD ) continue;
      li = 1.0 / sqrt(ll);

      /* Recovery of the two tangent vectors associated to points p1,p2; they
       * need not be oriented in the same fashion */
      if ( M_SIN(p1->tag) /* || (p1->tag & MG_NOM) */ ) {
        t1[0] = li*ux;
        t1[1] = li*uy;
      }
      else {
        t1[0] = p1->n[0];
        t1[1] = p1->n[1];
      }

      if ( M_SIN(p2->tag) /* || (p2->tag & MG_NOM)  */) {
        li = 1.0 / sqrt(ll);
        t2[0] = li*ux;
        t2[1] = li*uy;
      }
      else {
        t2[0] = p2->n[0];
        t2[1] = p2->n[1];
      }

      /* Calculation of the two Bezier coefficients along the curve */
      ps1   = ux*t1[0] + uy*t1[1];
      b1[0] = p1->c[0] + _MMG5_ATHIRD*ps1*t1[0];
      b1[1] = p1->c[1] + _MMG5_ATHIRD*ps1*t1[1];

      ps2   = ux*t2[0]+uy*t2[1];
      b2[0] = p2->c[0] - _MMG5_ATHIRD*ps2*t2[0];
      b2[1] = p2->c[1] - _MMG5_ATHIRD*ps2*t2[1];

      ps1 *= ps1;
      ps2 *= ps2;

      if ( ps1 < _MMG5_EPSD || ps2 < _MMG5_EPSD ) continue;

      /* \gamma^{\prime\prime}(0); \gamma^\prime(0) = ps*t1 by construction */
      gpp1[0] = 6.0*(p1->c[0] - 2.0*b1[0] + b2[0]);
      gpp1[1] = 6.0*(p1->c[1] - 2.0*b1[1] + b2[1]);

      /* Vector product gpp1 ^ t1 */
      pv = gpp1[0]*t1[1] - gpp1[1]*t1[0];
      M1 = fabs(pv)/ps1;

      /* \gamma^{\prime\prime}(1); \gamma^\prime(1) = -ps*t2 by construction */
      gpp2[0] = 6.0*(p2->c[0] - 2.0*b2[0] + b1[0]);
      gpp2[1] = 6.0*(p2->c[1] - 2.0*b2[1] + b1[1]);

      /* Vector product gpp2 ^ t2 */
      pv = gpp2[0]*t2[1] - gpp2[1]*t2[0];
      M2 = fabs(pv)/ps2;

      M1 = MG_MAX(M1,M2);
      if ( M1 < _MMG5_EPSD)
        lm = hmax;
      else {
        lm = 8.0*hausd / M1;
        lm = sqrt(lm);
      }
      met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
      met->m[ip2] = MG_MAX(hmin,MG_MIN(met->m[ip2],lm));
    }
  }

  return(1);
}
