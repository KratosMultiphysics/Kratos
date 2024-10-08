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
 * \file mmg2d/isosiz_2d.c
 * \brief Interpolation of metrics
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/
#include "libmmg2d_private.h"
#include "libmmg2d.h"
#include "mmg2dexterns_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param pt tetra to process.
 * \param i index of the edge of the tetra \a pt that we process.
 *
 * \return 1 if success, 0 if fail.
 *
 * Compute the euclidean length of the edge \a i of the tria \a pt, add this
 * length to the metric of the edge extremities and increment the count of times
 * we have processed this extremities.
 *
 */
int MMG2D_sum_reqEdgeLengthsAtPoint(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria pt,int8_t i) {
  MMG5_int         ip0,ip1;

  ip0 = pt->v[MMG5_iprv2[i]];
  ip1 = pt->v[MMG5_inxt2[i]];

  if ( !MMG5_sum_reqEdgeLengthsAtPoint(mesh,met,ip0,ip1) )
    return 0;

  return 1;
}


/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param ismet 1 if user provided metric
 *
 * \return 0 if fail, 1 otherwise
 *
 * Compute the metric at points on required edges as the mean of the lengths of
 * the required eges to which belongs the point. The processeed points are
 * marked with flag 3.
 *
 */
int MMG2D_set_metricAtPointsOnReqEdges ( MMG5_pMesh mesh,MMG5_pSol met, int8_t ismet ) {
  MMG5_pTria pt;
  MMG5_int   k,iadj;
  int        i;

  /* Reset the tria flag */
  for ( k=1; k<=mesh->nt; k++ ) {
    mesh->tria[k].flag = 0;
  }

  /* Reset the input metric at required edges extremities */
  if ( !MMG5_reset_metricAtReqEdges_surf (mesh, met, ismet ) ) {
    return 0;
  }

  /* Process the required edges and add the edge length to the metric of the
   * edge extremities */
  for ( k=1; k<=mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    /* Mark the tria as proceeded */
    pt->flag = 1;
    for ( i=0; i<3; i++ ) {
      if ( (pt->tag[i] & MG_REQ) || (pt->tag[i] & MG_NOSURF) ||
           (pt->tag[i] & MG_PARBDY) ) {
        /* Check if the edge has been proceeded by the neighbour triangle */
        iadj = mesh->adja[3*(k-1)+i+1];
        if ( iadj && mesh->tria[iadj/3].flag ) continue;

        if ( !MMG2D_sum_reqEdgeLengthsAtPoint(mesh,met,pt,i) ) {
          return 0;
        }
      }
    }
  }

  /* Travel the points and compute the metric of the points belonging to
   * required edges as the mean of the required edges length */
  if ( !MMG5_compute_meanMetricAtMarkedPoints ( mesh,met ) ) {
    return 0;
  }

return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 0 if fail, 1 otherwise
 *
 * New version for the definition of a size map; takes into account the
 * curvature of the external and internal curves present in the mesh
 *
 */
int MMG2D_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria  pt;
  MMG5_pPoint p0,p1,p2;
  MMG5_pPar   ppa;
  double      t1[2],t2[2],b1[2],b2[2],gpp1[2],gpp2[2],pv,M1,M2;
  double      ps1,ps2,ux,uy,ll,li,lm,hmax,hausd,hmin,lhmax,lhausd;
  MMG5_int    k,ip,ip1,ip2;
  int         l;
  int8_t      ismet;
  uint8_t     i,i1,i2;

  if ( !MMG5_defsiz_startingMessage (mesh,met,__func__) ) {
    return 0;
  }

  /* Reset the s and flag field : s is used to count the number of req edges to
   * which a req point belongs and flag is used to differenciate the point not
   * yet treated, the points with a required metric and the classical points */
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
    p0->s    = 0;
  }

  /** 1) Size at internal points */
  hmax = mesh->info.hmax;
  hausd = mesh->info.hausd;
  hmin = mesh->info.hmin;

  /* Allocate the structure */
  if ( !met->np ) {
    ismet = 0;

    /* Allocate and store the header informations for each solution */
    if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,MMG5_Scalar) ) {
      return 0;
    }
  }
  else {
    ismet = 1;
    assert ( met->m );
  }

  /** Step 1: Set metric at points belonging to a required edge: compute the
   * metric as the mean of the length of the required eges passing through the
   * point */
  if ( !mesh->info.nosizreq ) {
    if ( !MMG2D_set_metricAtPointsOnReqEdges ( mesh,met,ismet ) ) {
      return 0;
    }
  }

  /** Step 2: size at non required internal points */
  if ( !ismet ) {
    /* Initialize metric with a constant size */
    for ( k=1; k<=mesh->np; k++ ) {
      if ( mesh->point[k].flag ) {
        continue;
      }
      met->m[k] = hmax;
      mesh->point[k].flag = 1;
    }
  }

  /** Step 3: Minimum size feature imposed by the boundary edges */
  for ( k=1; k<=mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for ( i=0; i<3; i++ ) {

      if ( !MG_EDG(pt->tag[i]) ) continue;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];

      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      if ( p1->flag>1 && p2->flag>1 ) continue;

      lhmax = hmax;
      lhausd = hausd;

      /* Retrieve local parameters associated to edge i */
      if ( mesh->info.npar ) {
        for (l=0; l<mesh->info.npar; l++) {
          ppa = &mesh->info.par[l];
          if ( ppa->elt == MMG5_Edg && ppa->ref == pt->edg[i] ) {
            lhmax = ppa->hmax;
            lhausd = ppa->hausd;
            break;
          }
        }
      }


      ux = p2->c[0] - p1->c[0];
      uy = p2->c[1] - p1->c[1];
      ll = ux*ux + uy*uy;
      if ( ll < MMG5_EPSD ) continue;
      li = 1.0 / sqrt(ll);

      /* Recovery of the two tangent vectors associated to points p1,p2; they
       * need not be oriented in the same fashion */
      if ( MG_CRN & p1->tag || (p1->tag & MG_NOM) ) {
        t1[0] = li*ux;
        t1[1] = li*uy;
      }
      else {
        t1[0] = -p1->n[1];
        t1[1] = p1->n[0];
      }

      if ( MG_CRN & p2->tag || (p2->tag & MG_NOM) ) {
        li = 1.0 / sqrt(ll);
        t2[0] = li*ux;
        t2[1] = li*uy;
      }
      else {
        t2[0] = -p2->n[1];
        t2[1] = p2->n[0];
      }

      /* Calculation of the two Bezier coefficients along the curve */
      ps1   = ux*t1[0] + uy*t1[1];
      b1[0] = p1->c[0] + MMG5_ATHIRD*ps1*t1[0];
      b1[1] = p1->c[1] + MMG5_ATHIRD*ps1*t1[1];

      ps2   = ux*t2[0]+uy*t2[1];
      b2[0] = p2->c[0] - MMG5_ATHIRD*ps2*t2[0];
      b2[1] = p2->c[1] - MMG5_ATHIRD*ps2*t2[1];

      ps1 *= ps1;
      ps2 *= ps2;

      if ( ps1 < MMG5_EPSD || ps2 < MMG5_EPSD ) continue;

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
      if ( M1 < MMG5_EPSD )
        lm = lhmax;
      else {
        lm = 8.0*lhausd / M1;
        lm = MG_MIN(lhmax,sqrt(lm));
      }

      if ( p1->flag < 3 ) {
        met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
      }
      if ( p2->flag < 3 ) {
        met->m[ip2] = MG_MAX(hmin,MG_MIN(met->m[ip2],lm));
      }
    }
  }

  /** If local parameters are provided: size truncation on the entire mesh */
  /* Without local parameters information, only the boundary edges impose a
   * minimum size feature */
  if ( mesh->info.npar ) {
    /* Minimum size feature imposed by triangles */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      /* Retrieve local parameters associated to triangle k */
      for (l=0; l<mesh->info.npar; l++) {
        ppa = &mesh->info.par[l];
        if ( ppa->elt == MMG5_Triangle && ppa->ref == pt->ref ) {
          for (i=0; i<3; i++) {
            ip = pt->v[i];
            if ( mesh->point[ip].flag < 3 ) {
              met->m[ip] = MG_MAX(hmin,MG_MIN(met->m[ip],ppa->hmax));
            }
          }
          break;
        }
      }
    }
    /* Minimum size feature imposed by vertices */
    for (k=1; k<=mesh->np; k++) {
      p0 = &mesh->point[k];
      if ( (!MG_VOK(p0)) || p0->flag == 3 ) continue;

      /* Retrieve local parameters associated to vertex k */
      for (l=0; l<mesh->info.npar; l++) {
        ppa = &mesh->info.par[l];
        if ( ppa->elt == MMG5_Vertex && ppa->ref == p0->ref ) {
          met->m[k] = MG_MAX(hmin,MG_MIN(met->m[k],ppa->hmax));
          break;
        }
      }
    }
  }

  return 1;
}
