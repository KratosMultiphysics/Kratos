/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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
 * \file common/isosiz.c
 * \brief Fonctions for isotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the meric structure.
 * \param ptt pointer toward the triangle structure.
 * \return The computed area.
 *
 * Compute the area of the surface triangle \a ptt with respect to
 * the isotropic metric \a met.
 *
 */
double MMG5_surftri_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTria ptt) {
  double   *a,*b,*c,abx,aby,abz,acx,acy,acz,det,n[3];

  a = mesh->point[ptt->v[0]].c;
  b = mesh->point[ptt->v[1]].c;
  c = mesh->point[ptt->v[2]].c;

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  n[0] = aby*acz - abz*acy;
  n[1] = abz*acx - abx*acz;
  n[2] = abx*acy - aby*acx;
  det  = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

  return  0.5*sqrt(det) ;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param funcname name of the calling function
 *
 * \return 1 if success, 0 if fail.
 *
 * Print that we enter in the defsiz function in high verbosity level and check
 * the hmax value.
 *
 */
int MMG5_defsiz_startingMessage (MMG5_pMesh mesh,MMG5_pSol met,const char * funcname ) {

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining %stropic map\n",(met->size==1)?"iso":"aniso");

  if ( mesh->info.hmax < 0.0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hmax value.\n",funcname);
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Print gradation values (depending on the verbosity).
 *
 */
void MMG5_gradation_info ( MMG5_pMesh mesh ) {

  if ( mesh->info.imprim > 0 ) {
    if ( mesh->info.hgrad > 0. ) {
      fprintf(stdout,"\n  -- GRADATION : %8f ",
              exp(mesh->info.hgrad) );
      if ( mesh->info.hgradreq > 0. ) {
        fprintf(stdout,"(%8f)",exp(mesh->info.hgradreq));
      }
      fprintf(stdout,"\n");
    }
    else {
      if ( mesh->info.hgradreq > 0. ) {
        fprintf(stdout,"\n  -- GRADATION : DISABLED (%8f)\n",exp(mesh->info.hgradreq));
      }
    }
  }
  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ip0 index of the first edge extremity
 * \param ip1 index of the second edge extremity
 *
 * \return 1 if success, 0 if fail.
 *
 * Compute the euclidean length of the edge \a ip0 \a ip1,
 * add this length to the metric of the edge extremities and
 * increment the count of times we have processed this extremities.
 *
 */
int MMG5_sum_reqEdgeLengthsAtPoint ( MMG5_pMesh mesh,MMG5_pSol met,MMG5_int ip0,MMG5_int ip1 ) {
  MMG5_pPoint p0,p1;
  int    j;
  double      len,dist;

  /* Compute the euclidean edge length */
  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];

  len = 0.;
  for ( j=0; j<mesh->dim; ++j ) {
    dist = p1->c[j]-p0->c[j];
    len += dist*dist;
  }
  len = sqrt(len);

  /* Add the length to the point's metric and increment the number of
   * times the point has been seen */
  met->m[met->size*ip0] += len;
  met->m[met->size*ip1] += len;
  ++p0->s;
  ++p1->s;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Compute the mean metric at mesh points with a non-nul \a s field. At the
 * beginning, the metric of a given point contains the sum of n metrics and the
 * \a s field of the point the number of metrics summed in the point. Set the
 * flag of the processed points to 3.
 *
 */
int MMG5_compute_meanMetricAtMarkedPoints_iso ( MMG5_pMesh mesh,MMG5_pSol met ) {
  MMG5_pPoint p0;
  MMG5_int    k;
  int         mmgWarn = 0;

  for ( k=1; k<=mesh->np; k++ ) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) )  continue;

    if ( !p0->s ) continue;

    met->m[k] /= p0->s;
    p0->flag = 3;

    /* Warn the user that edge size is erased */
    if ( !mmgWarn ) {
      mmgWarn = 1;
      if ( mesh->info.ddebug || (mesh->info.imprim > 4) ) {
        printf("\n  -- SIZEMAP CORRECTION : overwritten of sizes at required vertices\n");
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ismet 1 if user provided metric
 *
 * \return 1 if success, 0 if fail.
 *
 * For a triangle mesh, process the triangles and set to 0 the metrics at points
 * that are at the extremities of a required edge.
 *
 */
int MMG5_reset_metricAtReqEdges_surf ( MMG5_pMesh mesh,MMG5_pSol met,int8_t ismet ) {
  MMG5_pTria  pt;
  int         i,j;
  MMG5_int    k,ip0,ip1,iad0,iad1;

  if ( ismet ) {
    for ( k=1; k<=mesh->nt; k++ ) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for ( i=0; i<3; ++i ) {
        if ( (pt->tag[i] & MG_REQ) || (pt->tag[i] & MG_NOSURF) ||
             (pt->tag[i] & MG_PARBDY) ) {

          ip0 = pt->v[MMG5_iprv2[i]];
          ip1 = pt->v[MMG5_inxt2[i]];
          iad0 = met->size*ip0;
          iad1 = met->size*ip1;

          for ( j=0; j<met->size; ++j ) {

            met->m[ iad0 + j ] = 0.;
            met->m[ iad1 + j ] = 0.;
          }
        }
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Set the s field of the points that belongs to a required edge to 1, set it to
 * 0 otherwise.
 *
 */
void MMG5_mark_pointsOnReqEdge_fromTria (  MMG5_pMesh mesh ) {
  MMG5_pTria  pt;
  MMG5_pPoint ppt;
  MMG5_int    k;
  int8_t      i;

  for ( k=1; k<=mesh->np; k++ ) {
    ppt = &mesh->point[k];
    ppt->s = 0;
  }

  /* Mark the points that belongs to a required edge */
  for ( k=1; k<=mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) { continue; }

    for ( i=0; i<3; ++i ) {
      if ( pt->tag[i] & MG_REQ ) {
        mesh->point[pt->v[MMG5_inxt2[i]]].s = 3*mesh->nt+2;
        mesh->point[pt->v[MMG5_iprv2[i]]].s = 3*mesh->nt+2;
      }
    }
  }
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 0 if fail, 1 otherwise
 *
 * Isotropic mesh gradation routine. The points belonging to a required edge are
 * treated in gradsizreq_iso.
 *
 */
int MMG5_gradsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            hgrad,ll,h1,h2,hn,val;
  int               it,maxit;
  MMG5_int          nup,nu;
  MMG5_int          ip1,ip2,k;
  int8_t            i,j,i1,i2;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  ** Grading mesh\n");
  }

  MMG5_mark_pointsOnReqEdge_fromTria ( mesh );

  for ( k=1; k<=mesh->np; k++ ) {
    mesh->point[k].flag = mesh->base;
  }


  hgrad = mesh->info.hgrad;
  it = 0;
  nup = 0;
  maxit = 100;

  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];
        if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;

        /* Skip points belonging to a required edge */
        if ( p1->s || p2->s ) continue;

        ll = 0.;
        for ( j=0; j<mesh->dim; ++j ) {
          val = p2->c[j]-p1->c[j];
          ll += val * val;
        }
        ll = sqrt(ll);

        h1 = met->m[ip1];
        h2 = met->m[ip2];
        if ( h1 < h2 ) {
          if ( h1 < MMG5_EPSD )  continue;
          hn  = h1 + hgrad*ll;
          if ( h2 > hn ) {
            met->m[ip2] = hn;
            p2->flag    = mesh->base;
            nu++;
          }
        }
        else {
          if ( h2 < MMG5_EPSD )  continue;
          hn = h2 + hgrad*ll;
          if ( h1 > hn ) {
            met->m[ip1] = hn;
            p1->flag    = mesh->base;
            nu++;
          }
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"     gradation: %7"MMG5_PRId" updated, %d iter.\n",nup,it);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return the number of updated metrics.
 *
 * Isotropic mesh gradation routine. The points belonging to a required entity
 * are treated in gradsizreq_iso.
 *
 */
int MMG5_gradsizreq_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            hgrad,ll,h1,h2,hn,ux,uy;
  int               it,maxit,nup,nu;
  MMG5_int          k,ip1,ip2,ipmaster,ipslave;
  uint8_t           i,i1,i2;


  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  ** Grading required points.\n");
  }

  if ( mesh->info.hgrad < 0. ) {
    /** Mark the edges belonging to a required entity */
    MMG5_mark_pointsOnReqEdge_fromTria ( mesh );
  }

  /** Update the sizes and mark the treated points */
  hgrad = mesh->info.hgradreq;
  it = nup = 0;
  maxit = 100;

  do {
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) {
        continue;
      }

      for (i=0; i<3; i++) {
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];

        if ( MMG5_abs ( p1->s - p2->s ) < 2 ) {
          /* No size to propagate */
          continue;
        }
        else if ( p1->s > p2->s ) {
          ipmaster = ip1;
          ipslave  = ip2;
        }
        else {
          assert ( p2->s > p1->s );
          ipmaster = ip2;
          ipslave  = ip1;
        }

        ux = p2->c[0]-p1->c[0];
        uy = p2->c[1]-p1->c[1];
        ll = ux*ux + uy*uy;
        ll = sqrt(ll);

        h1 = met->m[ipmaster];
        h2 = met->m[ipslave];
        if ( h1 < h2 ) {
          if ( h1 < MMG5_EPSD ) {
            continue;
          }
          hn  = h1 + hgrad*ll;
          if ( h2 <= hn ) {
            continue;
          }
        }
        else {
          hn = h1 - hgrad*ll;
          if ( h2 >= hn ) {
            continue;
          }
        }
        met->m[ipslave]           = hn;
        mesh->point[ipslave].s    = mesh->point[ipmaster].s - 1;
        nu++;
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 && nup ) {
    fprintf(stdout,"     gradation (required): %7d updated, %d iter.\n",nup,it);
  }

  return nup;
}
