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
 * \file mmg2d/movpt_2d.c
 * \brief Node relocation routines
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"

//extern int8_t ddb;

/**
 * \param mesh pointer to the mesh
 * \param met pointer to the metric structure.
 * \param list pointer to the ball of the point.
 * \param ilist size of the ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 0.9 of the old minimum element quality.
 *
 * \return 0 if we can't move the point (or if we fail), 1 if we can.
 *
 * Relocate boundary vertex whose ball is passed; routine works both in the
 * isotropic and anisotropic case
 *
 */
int MMG2D_movedgpt(MMG5_pMesh mesh,MMG5_pSol met,int ilist,MMG5_int *list, int8_t improve) {
  MMG5_pTria         pt,pt0;
  MMG5_pPoint        p0,p1,p2,ppt;
  double             step,ll1,ll2,o[2],no[2],calold,calnew;
  MMG5_int           k,iel,ip0,ip1,ip2,it1,it2;
  int8_t             i,i1,i2;
  static int8_t      mmgWarn0=0,mmgWarn1=0;

  pt0 = &mesh->tria[0];
  step = 0.1;
  ip1 = ip2 = it1 = it2 = 0;
  calold = calnew = DBL_MAX;

  /* First step: retrieve the two boundary points connected to p0 */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i   = list[k] % 3;
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    pt = &mesh->tria[iel];
    calold = MG_MIN(MMG2D_caltri(mesh,met,pt),calold);

    if ( MG_EDG(pt->tag[i1]) ) {
      if ( ip1 == 0 ) {
        ip1 = pt->v[i2];
        it1 = 3*iel+i1;
      }
      else if ( ip1 != pt->v[i2] ) {
        if ( ip2 == 0 ) {
          ip2 = pt->v[i2];
          it2 = 3*iel+i1;
        }
        else if ( ip2 != pt->v[i2] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the "
                    "intersection of 3 edges. abort.\n",__func__);
          }
          return 0;
        }
      }
    }

    if ( MG_EDG(pt->tag[i2]) ) {
      if ( ip1 == 0 ) {
        ip1 = pt->v[i1];
        it1 = 3*iel+i2;
      }
      else if ( ip1 != pt->v[i1] ) {
        if ( ip2 == 0 ) {
          ip2 = pt->v[i1];
          it2 = 3*iel+i2;
        }
        else if ( ip2 != pt->v[i1] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Warning: %s: at least 1 point at the "
                    "intersection of 3 edges. abort.\n",__func__);
          }
         return 0;
        }
      }
    }
  }

  /* Check that there are exactly two boundary points connected at p0 */
  if ( ip1 == 0 || ip2 == 0 ) {
    if ( !mmgWarn1 ) {
      mmgWarn1 = 1;
      fprintf(stderr,"\n  ## Warning: %s: non singular point at end of edge.\n",
              __func__);
    }
    return 0;
  }

  ip0 = pt->v[i];
  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  /* Calculate length of both edges */
  /* Anisotropic case: ll1 and ll2 = anisotropic edge lengths */
  if ( met->m && met->size == 3 ) {
    ll1 = MMG2D_lencurv_ani(mesh,met,ip0,ip1);
    ll2 = MMG2D_lencurv_ani(mesh,met,ip0,ip2);
  }
  /* In all the remaining cases, ll1 and ll2 = squared edge lengths;
     inconsistency between both cases is not problematic since these values serve only for comparison */
  else {
    ll1 = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0]) + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1]);
    ll2 = (p2->c[0]-p0->c[0])*(p2->c[0]-p0->c[0]) + (p2->c[1]-p0->c[1])*(p2->c[1]-p0->c[1]);
  }

  /* Relocate p0 slightly towards p1 */
  if ( ll1 > ll2 ) {
    iel = it1 / 3;
    i   = it1 % 3;
  }
  /* Relocate p0 slightly towards p2 */
  else {
    iel = it2 / 3;
    i   = it2 % 3;
  }

  i1 = MMG5_inxt2[i];

  pt = &mesh->tria[iel];

  /* step = distance of the relocated position from ip0
     bezierCurv = distance s from inxt2[i] */
  if ( pt->v[i1] == ip0 ) MMG2D_bezierCurv(mesh,iel,i,step,o,no);
  else MMG2D_bezierCurv(mesh,iel,i,1.0-step,o,no);

  /* Evaluate resulting configuration */
  ppt = &mesh->point[0];
  ppt->c[0] = o[0];
  ppt->c[1] = o[1];
  ppt->n[0] = no[0];
  ppt->n[1] = no[1];

  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i   = list[k] % 3;
    pt  = &mesh->tria[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i] = 0;

    calnew = MG_MIN(MMG2D_caltri(mesh,met,pt0),calnew);
  }

  if ( calold < MMG2D_NULKAL && calnew <= calold ) return 0;
  else if ( calnew < MMG2D_NULKAL ) return 0;
  else if ( improve && calnew < 1.02 * calold ) return 0;
  else if ( calnew < 0.3 * calold ) return 0;

  /* Update of the coordinates and normal of the point */
  p0 = &mesh->point[pt->v[i]];
  p0->c[0] = o[0];
  p0->c[1] = o[1];

  p0->n[0] = no[0];
  p0->n[1] = no[1];

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param met pointer to the metric structure.
 * \param list pointer to the ball of the point.
 * \param ilist size of the ball.
 * \param improve force the new minimum element quality to be greater or equal
 * than 0.9 of the old minimum element quality.
 *
 * \return 0 if we can't move the point (or if we fail), 1 if we can.
 *
 * Relocate internal vertex whose ball is passed.
 *
 */
int MMG2D_movintpt(MMG5_pMesh mesh,MMG5_pSol met,int ilist,MMG5_int *list,int8_t improve) {
  MMG5_pTria        pt,pt0;
  MMG5_pPoint       p0,p1,p2,ppt0;
  double            calold,calnew,vol,volbal,b[2];
  MMG5_int          k,iel;
  int8_t            i,i1,i2;

  ppt0 = &mesh->point[0];
  pt0  = &mesh->tria[0];

  volbal = 0.0;
  b[0] = b[1] = 0.0;
  calold = calnew = DBL_MAX;

  /* Calculate the now position for vertex, as well as the quality of the previous configuration */
  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i   = list[k] % 3;
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    pt = &mesh->tria[iel];

    /* Volume of iel */
    p0 = &mesh->point[pt->v[i]];
    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt->v[i2]];
    vol = 0.5* fabs((p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]));

    volbal += vol;

    /* Add coordinates of the centre of mass of iel, weighted by its volume */
    b[0] += MMG5_ATHIRD*vol*(p0->c[0]+p1->c[0]+p2->c[0]);
    b[1] += MMG5_ATHIRD*vol*(p0->c[1]+p1->c[1]+p2->c[1]);

    /* Quality of pt */
    calold = MG_MIN(MMG2D_caltri_iso(mesh,NULL,pt),calold);
  }

  if ( volbal < MMG5_EPSD ) return 0;
  volbal = 1.0 / volbal;
  b[0] *= volbal;
  b[1] *= volbal;

  /* Check the quality of the resulting configuration */
  ppt0->c[0] = b[0];
  ppt0->c[1] = b[1];

  for (k=0; k<ilist; k++) {
    iel = list[k] / 3;
    i   = list[k] % 3;
    pt  = &mesh->tria[iel];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i] = 0;

    calnew = MG_MIN(MMG2D_caltri_iso(mesh,NULL,pt0),calnew);
  }

  if (calold < MMG2D_NULKAL && calnew <= calold) return 0;
  else if (calnew < MMG2D_NULKAL) return 0;
  else if ( improve && calnew < 1.02 * calold ) return 0;
  else if ( calnew < 0.3 * calold ) return 0;

  /* Update of the coordinates of the point */
  pt = &mesh->tria[list[0]/3];
  i  = list[0] % 3;
  p0 = &mesh->point[pt->v[i]];
  p0->c[0] = b[0];
  p0->c[1] = b[1];

  return 1;
}
