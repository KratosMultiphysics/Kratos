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
 * \file mmgs/isosiz_s.c
 * \brief Fonctions for isotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmgs_private.h"
#include "libmmgs.h"
#include <math.h>
#include "mmgsexterns_private.h"
#include "mmgexterns_private.h"

#define MAXLEN   1.0e+3

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param hash edge hashtable.
 * \param pt tria to process.
 * \param i index of the edge of the tria \a pt that we process.
 *
 * \return 1 if success, 0 if fail.
 *
 * If the edge \a i of the tria \a pt is seen for the first time, compute its
 * euclidean length, add this length to the metric of the edge extremities and
 * increment the count of times we have processed this extremities.
 *
 */
static inline
int MMGS_sum_reqEdgeLengthsAtPoint(MMG5_pMesh mesh,MMG5_pSol met,MMG5_Hash *hash,
                                  MMG5_pTria pt,int8_t i) {
  MMG5_int         ip0,ip1;

  ip0 = pt->v[MMG5_inxt2[i]];
  ip1 = pt->v[MMG5_iprv2[i]];

  /* Check if the edge is already treated */
  if ( MMG5_hashGet(hash,ip0,ip1) ) return 1;

  /* Mark the edge as treated */
  if ( !MMG5_hashEdge(mesh,hash,ip0,ip1,1) ) return 0;

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
 * Compute the metric at points on trequired adges as the mean of the lengths of
 * the required eges to which belongs the point. The processeed points are
 * marked with flag 3.
 *
 */
int MMGS_set_metricAtPointsOnReqEdges ( MMG5_pMesh mesh,MMG5_pSol met,int8_t ismet ) {
  MMG5_pTria pt;
  MMG5_Hash  hash;
  int        i;
  MMG5_int   k;

  /* Reset the input metric at required edges extremities */
  if ( !MMG5_reset_metricAtReqEdges_surf (mesh, met,ismet ) ) {
    return 0;
  }

  /* Process the required edges and add the edge length to the metric of the
   * edge extremities */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return 0;

  for ( k=1; k<=mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for ( i=0; i<3; i++ ) {
      if ( (pt->tag[i] & MG_REQ) || (pt->tag[i] & MG_NOSURF) ||
           (pt->tag[i] & MG_PARBDY) ) {

        /* Check if the edge has been proceeded by the neighbour triangle */
        if ( !MMGS_sum_reqEdgeLengthsAtPoint(mesh,met,&hash,pt,i) ) {
          MMG5_DEL_MEM(mesh,hash.item);
          return 0;
        }
      }
    }
  }
  MMG5_DEL_MEM(mesh,hash.item);

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
 * \return 1 if success, 0 if fail
 *
 * Define isotropic size map at all vertices of the mesh, associated with
 * geometric approx ; by convention, p0->h stores desired length at point p0
 *
 */
int MMGS_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria  pt;
  MMG5_pPoint p[3],p0;
  MMG5_pPar   par;
  double      n[3][3],t[3][3],nt[3],c1[3],c2[3],*n1,*n2,*t1,*t2;
  double      ps,ps2,ux,uy,uz,ll,l,lm,dd,M1,M2,hausd,hmin,hmax;
  int         j,isloc;
  MMG5_int    k,ip1,ip2;
  int8_t      ismet;
  int8_t      i,i1,i2;

  if ( !MMG5_defsiz_startingMessage (mesh,met,__func__) ) {
    return 0;
  }

  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
    p0->s    = 0;
  }

  /* alloc structure */
  if ( !met->m ) {
    ismet = 0;

    /* Allocate and store the header informations for each solution */
    if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,MMG5_Scalar) ) {
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
    if ( !MMGS_set_metricAtPointsOnReqEdges ( mesh,met,ismet ) ) {
      return 0;
    }
  }

  /** Step 2: size at non required internal points */
  if ( !ismet ) {

    /* init constant size */
    for (k=1; k<=mesh->np; k++) {
      if ( mesh->point[k].flag ) {
        continue;
      }
      met->m[k] = mesh->info.hmax;
      mesh->point[k].flag = 1;
    }
  }

  /** Step 3: Minimum size feature imposed by the hausdorff distance */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    p[0] = &mesh->point[pt->v[0]];
    p[1] = &mesh->point[pt->v[1]];
    p[2] = &mesh->point[pt->v[2]];

    /* normal recovery */
    for (i=0; i<3; i++) {
      if ( MS_SIN(p[i]->tag) ) {
        MMG5_nortri(mesh,pt,n[i]);
      }
      else if ( MG_EDG(p[i]->tag) ) {
        MMG5_nortri(mesh,pt,nt);
        n1  = &mesh->xpoint[p[i]->xp].n1[0];
        n2  = &mesh->xpoint[p[i]->xp].n2[0];
        ps  = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
        ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];
        if ( fabs(ps) > fabs(ps2) )
          memcpy(&n[i],n1,3*sizeof(double));
        else
          memcpy(&n[i],n2,3*sizeof(double));
        memcpy(&t[i],p[i]->n,3*sizeof(double));
      }
      else
        memcpy(&n[i],p[i]->n,3*sizeof(double));
    }

    for (i=0; i<3; i++) {
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];

      /* local parameters */
      hausd = mesh->info.hausd;
      hmin  = mesh->info.hmin;
      hmax  = mesh->info.hmax;
      isloc = 0;
      for (j=0; j<mesh->info.npar; j++) {
        par = &mesh->info.par[j];
        if ( (par->elt == MMG5_Triangle) && (pt->ref == par->ref ) ) {
          if ( !isloc ) {
            hausd = par->hausd;
            hmin  = par->hmin;
            hmax  = par->hmax;
            isloc = 1;
          }
          else {
            hausd = MG_MIN(par->hausd,hausd);
            hmin  = MG_MAX(par->hmin,hmin);
            hmax  = MG_MIN(par->hmax,hmax);
          }
        }
      }

      ux = p[i2]->c[0] - p[i1]->c[0];
      uy = p[i2]->c[1] - p[i1]->c[1];
      uz = p[i2]->c[2] - p[i1]->c[2];
      ll = ux*ux + uy*uy + uz*uz;

      if ( ll < MMG5_EPSD )  continue;

      if ( MG_EDG(pt->tag[i]) ) {
        if ( MS_SIN(p[i1]->tag) ) {
          t[i1][0] = ux;
          t[i1][1] = uy;
          t[i1][2] = uz;
          l = 1.0 / sqrt(ll);
          t[i1][0] *= l;
          t[i1][1] *= l;
          t[i1][2] *= l;
        }
        if ( MS_SIN(p[i2]->tag) ) {
          t[i2][0] = -ux;
          t[i2][1] = -uy;
          t[i2][2] = -uz;
          l = 1.0/sqrt(ll);
          t[i2][0] *= l;
          t[i2][1] *= l;
          t[i2][2] *= l;
        }
        t1 = t[i1];
        t2 = t[i2];

        /* The two Bezier coefficients along curve */
        dd    = (t1[0]*ux + t1[1]*uy + t1[2]*uz)/3.0;
        c1[0] = p[i1]->c[0] + dd * t1[0];
        c1[1] = p[i1]->c[1] + dd * t1[1];
        c1[2] = p[i1]->c[2] + dd * t1[2];

        dd    = -(t2[0]*ux + t2[1]*uy + t2[2]*uz)/3.0;
        c2[0] = p[i2]->c[0] + dd * t2[0];
        c2[1] = p[i2]->c[1] + dd * t2[1];
        c2[2] = p[i2]->c[2] + dd * t2[2];

        M1 = (c2[0]-2.0*c1[0]+p[i1]->c[0])*(c2[0]-2.0*c1[0]+p[i1]->c[0]) \
          + (c2[1]-2.0*c1[1]+p[i1]->c[1])*(c2[1]-2.0*c1[1]+p[i1]->c[1]) \
          + (c2[2]-2.0*c1[2]+p[i1]->c[2])*(c2[2]-2.0*c1[2]+p[i1]->c[2]);

        M2 = (p[i2]->c[0]-2.0*c2[0]+c1[0])*(p[i2]->c[0]-2.0*c2[0]+c1[0]) \
          + (p[i2]->c[1]-2.0*c2[1]+c1[1])*(p[i2]->c[1]-2.0*c2[1]+c1[1])\
          + (p[i2]->c[2]-2.0*c2[2]+c1[2])*(p[i2]->c[2]-2.0*c2[2]+c1[2]);

        M1 = 6.0 * sqrt(M1);
        M2 = 6.0 * sqrt(M2);
        M1 = MG_MAX(M1,M2);

        if ( M1 < MMG5_EPSD )
          lm = MAXLEN;
        else {
          lm = (16.0*ll*hausd) / (3.0*M1);
          lm = sqrt(lm);
        }
        if ( mesh->point[ip1].flag < 3 ) {
          met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
        }
        if ( mesh->point[ip2].flag < 3 ) {
          met->m[ip2] = MG_MAX(hmin,MG_MIN(met->m[ip2],lm));
        }
      }
      else {
        n1 = n[i1];
        n2 = n[i2];

        ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
        c1[0] = (2.0*p[i1]->c[0] + p[i2]->c[0] - ps*n1[0]) / 3.0;
        c1[1] = (2.0*p[i1]->c[1] + p[i2]->c[1] - ps*n1[1]) / 3.0;
        c1[2] = (2.0*p[i1]->c[2] + p[i2]->c[2] - ps*n1[2]) / 3.0;

        ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
        c2[0] = (2.0*p[i2]->c[0] + p[i1]->c[0] - ps*n2[0]) / 3.0;
        c2[1] = (2.0*p[i2]->c[1] + p[i1]->c[1] - ps*n2[1]) / 3.0;
        c2[2] = (2.0*p[i2]->c[2] + p[i1]->c[2] - ps*n2[2]) / 3.0;

        M1 = (c2[0]-2.0*c1[0]+p[i1]->c[0])*(c2[0]-2.0*c1[0]+p[i1]->c[0]) \
          + (c2[1]-2.0*c1[1]+p[i1]->c[1])*(c2[1]-2.0*c1[1]+p[i1]->c[1]) \
          + (c2[2]-2.0*c1[2]+p[i1]->c[2])*(c2[2]-2.0*c1[2]+p[i1]->c[2]);

        M2 = (p[i2]->c[0]-2.0*c2[0]+c1[0])*(p[i2]->c[0]-2.0*c2[0]+c1[0]) \
          + (p[i2]->c[1]-2.0*c2[1]+c1[1])*(p[i2]->c[1]-2.0*c2[1]+c1[1])\
          + (p[i2]->c[2]-2.0*c2[2]+c1[2])*(p[i2]->c[2]-2.0*c2[2]+c1[2]);

        M1 = 6.0 * sqrt(M1);
        M2 = 6.0 * sqrt(M2);
        M1 = MG_MAX(M1,M2);

        if ( M1 < MMG5_EPSD )
          lm = MAXLEN;
        else {
          lm = (16.0*ll*hausd) / (3.0*M1);
          lm = sqrt(lm);
        }
        if ( mesh->point[ip1].flag < 3 ) {
          met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
        }
        if ( mesh->point[ip2].flag < 3 ) {
          met->m[ip2] = MG_MAX(hmin,MG_MIN(met->m[ip2],lm));
        }
      }
    }
  }

  /* take local parameters */
  for (j=0; j<mesh->info.npar; j++) {
    par = &mesh->info.par[j];
    if ( par->elt == MMG5_Triangle ) {
      for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MG_EOK(pt) || pt->ref != par->ref )  continue;
        if ( mesh->point[pt->v[0]].flag < 3 ) {
          met->m[pt->v[0]] = MG_MAX(par->hmin,MG_MIN(met->m[pt->v[0]],par->hmax));
        }
        if ( mesh->point[pt->v[1]].flag < 3 ) {
          met->m[pt->v[1]] = MG_MAX(par->hmin,MG_MIN(met->m[pt->v[1]],par->hmax));
        }
        if ( mesh->point[pt->v[2]].flag < 3 ) {
          met->m[pt->v[2]] = MG_MAX(par->hmin,MG_MIN(met->m[pt->v[2]],par->hmax));
        }
      }
    }
  }
  return 1;
}
