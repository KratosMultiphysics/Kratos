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
 * \file mmg3d/quality_3d.c
 * \brief Functions to compute elements quality and edge lengths.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg3d.h"
#include "inlined_functions_3d_private.h"
#include "mmg3dexterns_private.h"

extern int8_t ddb;

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the meric structure.
 * \param metRidTyp metric storage (classic or special)
 * \return 1 if success, 0 if fail.
 *
 * Compute the quality of the tetras over the mesh.
 *
 */
int MMG3D_tetraQual(MMG5_pMesh mesh, MMG5_pSol met,int8_t metRidTyp) {
  MMG5_pTetra pt;
  double      minqual;
  MMG5_int    k,iel;

  minqual = 2./MMG3D_ALPHAD;

  /* compute tet quality */
  iel = 1;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) )   continue;

    if ( !metRidTyp && met->size == 6 && met->m ) {
      pt->qual = MMG5_caltet33_ani(mesh,met,pt);
    }
    else if ( !(met && met->m) ) {
      pt->qual = MMG5_caltet_iso(mesh,NULL,pt);
    }
    else {
      pt->qual = MMG5_orcal(mesh,met,k);
    }

    int i=0;
    /* Once metric is stored using 'ridge' convention, a tetra with 4 ridge
     * points has a 0 quality, ignore it for quality checks */
    if ( metRidTyp ) {
      for ( i=0; i<4; ++i ) {
        MMG5_pPoint ppt = &mesh->point[pt->v[i]];
        if ( (MG_SIN(ppt->tag) || MG_NOM & ppt->tag) || !(ppt->tag & MG_GEO) ) {
          break;
        }
      }
    }

    /* Check quality on suitable elements */
    if ( i < 4 && pt->qual < minqual ) {
      minqual = pt->qual;
      iel     = k;
    }
  }

  /* Here the quality is not normalized by alpha, thus we need to
   * normalized it */
  return  MMG5_minQualCheck(iel,minqual,MMG3D_ALPHAD);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the meric structure.
 * \param pt pointer to a tetrahedra.
 * \return The anisotropic quality of the tet or 0.0 if fail.
 *
 * Compute the quality of the tet pt with respect to the anisotropic metric \a
 * met. \f$ Q = V_met(K) / (sum(len(edge_K)^2)^(3/2) \f$ and for a calssic
 * storage of metrics at ridges.
 *
 * \todo test with the square of this measure
 */
inline double MMG5_caltet33_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt) {
  double       cal,abx,aby,abz,acx,acy,acz,adx,ady,adz;
  double       bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double       h1,h2,h3,h4,h5,h6,det,vol,rap,v1,v2,v3,num;
  double       *a,*b,*c,*d;
  double       mm[6];
  int          iad0,iad1,iad2,iad3;
  MMG5_int     ip[4],k;

  ip[0] = pt->v[0];
  ip[1] = pt->v[1];
  ip[2] = pt->v[2];
  ip[3] = pt->v[3];

  iad0  = met->size * ip[0];
  iad1  = met->size * ip[1];
  iad2  = met->size * ip[2];
  iad3  = met->size * ip[3];

  /* average metric */
  for (k=0; k<6; k++)
    mm[k] = 0.25 * (met->m[iad0+k]+met->m[iad1+k]+met->m[iad2+k]+met->m[iad3+k]);

  a = mesh->point[ip[0]].c;
  b = mesh->point[ip[1]].c;
  c = mesh->point[ip[2]].c;
  d = mesh->point[ip[3]].c;

  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];

  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];

  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol <= 0. )  return 0.0;

  det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
    - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
    + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
  if ( det < MMG5_EPSD2 )   {
    return 0.0;
  }
  det = sqrt(det) * vol;

  /* edge lengths */
  h1 = mm[0]*abx*abx + mm[3]*aby*aby + mm[5]*abz*abz
    + 2.0*(mm[1]*abx*aby + mm[2]*abx*abz + mm[4]*aby*abz);
  h2 =  mm[0]*acx*acx + mm[3]*acy*acy + mm[5]*acz*acz
    + 2.0*(mm[1]*acx*acy + mm[2]*acx*acz + mm[4]*acy*acz);
  h3 = mm[0]*adx*adx + mm[3]*ady*ady + mm[5]*adz*adz
    + 2.0*(mm[1]*adx*ady + mm[2]*adx*adz + mm[4]*ady*adz);
  h4 =  mm[0]*bcx*bcx + mm[3]*bcy*bcy + mm[5]*bcz*bcz
    + 2.0*(mm[1]*bcx*bcy + mm[2]*bcx*bcz + mm[4]*bcy*bcz);
  h5 =  mm[0]*bdx*bdx + mm[3]*bdy*bdy + mm[5]*bdz*bdz
    + 2.0*(mm[1]*bdx*bdy + mm[2]*bdx*bdz + mm[4]*bdy*bdz);
  h6 =  mm[0]*cdx*cdx + mm[3]*cdy*cdy + mm[5]*cdz*cdz
    + 2.0*(mm[1]*cdx*cdy + mm[2]*cdx*cdz + mm[4]*cdy*cdz);

  /* quality */
  rap = h1 + h2 + h3 + h4 + h5 + h6;
  num = sqrt(rap) * rap;

  cal = det / num;

  return cal;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param avlen average length (to fill).
 * \param lmin minimal length (to fill).
 * \param lmax max length (to fill).
 * \param ned number of edges (to fill).
 * \param amin (to fill).
 * \param bmin (to fill).
 * \param amax (to fill).
 * \param bmax (to fill).
 * \param nullEdge (to fill).
 * \param metRidTyp (to fill).
 * \param bd_in (to fill).
 * \param hl (to fill).
 *
 * \return 0 if fail, 1 otherwise.
 *
 * Compute the required information to print the length histogram
 *
 */
int MMG3D_computePrilen( MMG5_pMesh mesh, MMG5_pSol met, double* avlen,
                         double* lmin, double* lmax, MMG5_int* ned, MMG5_int* amin,
                         MMG5_int* bmin, MMG5_int* amax,
                         MMG5_int* bmax, MMG5_int* nullEdge, int8_t metRidTyp,
                         double** bd_in, MMG5_int hl[9] )
{
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  MMG5_Hash       hash;
  double          len;
  MMG5_int        k,np,nq,n;
  int8_t          ia,i0,i1,ier,i;
  static double   bd[9]= {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};

  *bd_in = bd;
  memset(hl,0,9*sizeof(MMG5_int));
  *ned = 0;
  *avlen = 0.0;
  *lmax = 0.0;
  *lmin = 1.e30;
  *amin = *amax = *bmin = *bmax = 0;
  *nullEdge = 0;

  if ( (!met) || (!met->m) ) {
    /* the functions that computes the edge length cannot be called without an
     * allocated metric */
    return 0;
  }

  if ( !mesh->ne ) {
    return 0;
  }

  /* Hash all edges in the mesh */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return 0;

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = MMG5_iare[ia][0];
      i1 = MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!MMG5_hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stderr,"  ## Error: %s: function MMG5_hashEdge return 0\n",
                __func__);
        return 0;
      }
    }
  }

  /* Pop edges from hash table, and analyze their length */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    n = 0;
    for(i=0 ; i<4 ; i++) {
      ppt = &mesh->point[pt->v[i]];
      if(!(MG_SIN(ppt->tag) || MG_NOM & ppt->tag) && (ppt->tag & MG_GEO)) continue;
      n++;
    }
    if(!n) {
      continue;
    }
    for(ia=0; ia<6; ia++) {
      i0 = MMG5_iare[ia][0];
      i1 = MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = MMG5_hashPop(&hash,np,nq);
      if( ier ) {
        if ( (!metRidTyp) && met->size==6 && met->m ) {
          // Warning: we may erroneously approximate the length of a curve
          // boundary edge by the length of the straight edge if the "MG_BDY"
          // tag is missing along the edge.
          len = MMG5_lenedg33_ani(mesh,met,ia,pt);
        }
        else
          // Warning: we may erroneously approximate the length of a curve
          // boundary edge by the length of the straight edge if the "MG_BDY"
          // tag is missing along the edge.
          len = MMG5_lenedg(mesh,met,ia,pt);


        if ( !len ) {
          ++(*nullEdge);
        }
        else {
          *avlen += len;
          (*ned)++;

          if( len < (*lmin) ) {
            *lmin = len;
            *amin = np;
            *bmin = nq;
          }

          if ( len > (*lmax) ) {
            *lmax = len;
            *amax = np;
            *bmax = nq;
          }

          /* Locate size of edge among given table */
          for(i=0; i<8; i++) {
            if ( bd[i] <= len && len < bd[i+1] ) {
              hl[i]++;
              break;
            }
          }
          if( i == 8 ) hl[8]++;
        }
      }
    }
  }

  MMG5_DEL_MEM(mesh,hash.item);
  return 1;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage,
 * 1 for special storage.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute sizes of edges of the mesh, and displays histo.
 *
 */
int MMG3D_prilen(MMG5_pMesh mesh, MMG5_pSol met, int8_t metRidTyp) {
  double     avlen, lmin, lmax;
  MMG5_int   ned, nullEdge;
  MMG5_int   amin, bmin, amax, bmax, hl[9];
  double     *bd;

  if (!MMG3D_computePrilen( mesh, met, &avlen, &lmin, &lmax, &ned, &amin,
                            &bmin, &amax, &bmax, &nullEdge, metRidTyp, &bd, hl ) )
    return 0;

  /* Display histogram */
  MMG5_displayLengthHisto(mesh, ned, &avlen, amin, bmin, lmin,
                           amax, bmax, lmax,nullEdge, &bd[0], &hl[0],1);

  return 1;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param ne pointer to the number of used tetra (to fill).
 * \param max pointer to the maximal quality (normalized) to fill.
 * \param avg pointer to the average quality (normalized) to fill.
 * \param min pointer to the minimal quality (normalized) to fill.
 * \param iel pointer to the index of the worst tetra (to fill).
 * \param good number of good elements (to fill).
 * \param med number of elements with a quality greather than 0.5 (to fill).
 * \param his pointer to the mesh histogram (to fill).
 * \param imprim verbosity level
 *
 * Compute the needed quality information in order to print the quality histogram
 * in optimLES mode.
 *
 */
void MMG3D_computeLESqua(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *ne,double *max,double *avg,
                         double *min,MMG5_int *iel,MMG5_int *good,MMG5_int *med,MMG5_int his[5],int imprim) {
  MMG5_pTetra    pt;
  double         rap;
  MMG5_int       k,ok,nex;
  static int8_t  mmgWarn0=0;

  /*compute tet quality*/
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) )   continue;

    pt->qual = MMG5_orcal(mesh,met,k);
  }

  if ( imprim <= 0 )
    return;

  (*min)  = (*avg) = 0.0;
  (*max)  = 1.0;
  (*iel)  = 0;
  (*med) = (*good) = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;
    if ( (!mmgWarn0) && (MMG5_orvol(mesh->point,pt->v) < 0.0) ) {
      mmgWarn0 = 1;
      fprintf(stderr,"  ## Warning: %s: at least 1 negative volume.\n",
              __func__);
    }
    rap = 1. - MMG3D_ALPHAD * pt->qual;
    if ( rap > (*min) ) {
      (*min) = rap;
      (*iel) = ok;
    }
    if ( rap < 0.9 )  (*med)++;
    if ( rap < 0.6 )  (*good)++;
    // if ( rap < MMG3D_BADKAL )  mesh->info.badkal = 1;
    (*avg) += rap;
    (*max)  = MG_MIN((*max),rap);
    if(rap < 0.6)
      his[0] += 1;
    else if(rap < 0.9)
      his[1] += 1;
    else if(rap < 0.93)
      his[2] += 1;
    else if(rap < 0.99)
      his[3] += 1;
    else
      his[4] += 1;
  }

  (*ne) = mesh->ne-nex;

  return;
}

/**
 * \param ne number of used tetra.
 * \param max maximal quality (normalized).
 * \param avg average quality (normalized).
 * \param min minimal quality (normalized).
 * \param iel index of the worst tetra.
 * \param good number of good elements.
 * \param med number of elements with a quality greather than 0.5
 * \param his pointer to the mesh histogram.
 * \param nrid number of tetra with 4 ridge points if we want to warn the user.
 * \param optimLES 1 if we work in optimLES mode, 0 otherwise
 * \param imprim verbosity level
 *
 * \return 0 if the worst element has a nul quality, 1 otherwise.
 *
 * Print the header of the histogram of mesh qualities then call the function
 * that print the histogram for special metric at ridges.
 *
 */
int MMG3D_displayQualHisto(MMG5_int ne,double max,double avg,double min,MMG5_int iel,
                           MMG5_int good,MMG5_int med,MMG5_int his[5],MMG5_int nrid,int optimLES,
                           int imprim) {

  fprintf(stdout,"\n  -- MESH QUALITY");
  if ( optimLES )
    fprintf(stdout," (LES)");
  fprintf(stdout,"  %" MMG5_PRId "\n",ne);

  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%" MMG5_PRId ")\n",
          max,avg / ne,min,iel);

  return ( MMG3D_displayQualHisto_internal(ne,max,avg,min,iel,good,med,his,
                                           nrid,optimLES,imprim) );
}

/**
 * \param ne number of used tetra.
 * \param max maximal quality (normalized).
 * \param avg average quality (normalized).
 * \param min minimal quality (normalized).
 * \param iel index of the worst tetra.
 * \param good number of good elements.
 * \param med number of elements with a quality greather than 0.5
 * \param his pointer to the mesh histogram.
 * \param nrid number of tetra with 4 ridge points if we want to warn the user.
 * \param optimLES 1 if we work in optimLES mode, 0 otherwise
 * \param imprim verbosity level
 *
 * \return 0 if the worst element has a nul quality, 1 otherwise.
 *
 * Print histogram of mesh qualities for special storage of metric at ridges.
 *
 */
int MMG3D_displayQualHisto_internal(MMG5_int ne,double max,double avg,double min,MMG5_int iel,
                                    MMG5_int good,MMG5_int med,MMG5_int his[5],MMG5_int nrid,int optimLES,
                                    int imprim)
{
  const double les_ticks[6] = {0,0.6,0.9,0.93,0.99,1};
  int          i,imax;

  if ( abs(imprim) >= 3 ){
    if ( optimLES ) {
      /* print histo */
      fprintf(stdout,"     HISTOGRAMM:");
      fprintf(stdout,"  %6.2f %% < 0.6\n",100.0*((float)good/(float)ne));
      if ( abs(imprim) > 3 ) {
        fprintf(stdout,"                  %6.2f %% < 0.9\n",100.0*( (float)med/(float)ne));

        assert ( min >= max );
        for ( i=0; i<5; ++i ) {
          if ( max < les_ticks[i+1] && min >= les_ticks[i] ) {
            fprintf(stdout,"     %5.2f < Q < %5.2f   %7"MMG5_PRId"   %6.2f %%\n",
                    les_ticks[i],les_ticks[i+1],his[i],
                    100.*((float)his[i]/(float)ne));
          }
        }
      }
      return 1;
    }
    else {
      /* print histo */
      fprintf(stdout,"     HISTOGRAMM:");
      fprintf(stdout,"  %6.2f %% > 0.12\n",100.0*((float)good/(float)ne));
      if ( abs(imprim) > 3 ) {
        fprintf(stdout,"                  %6.2f %% >  0.5\n",100.0*( (float)med/(float)ne));
        imax = MG_MIN(4,(int)(5.*max));
        for (i=imax; i>=(int)(5*min); i--) {
          fprintf(stdout,"     %5.1f < Q < %5.1f   %7"MMG5_PRId"   %6.2f %%\n",
                  i/5.,i/5.+0.2,his[i],100.*((float)his[i]/(float)ne));
        }
        if ( nrid ) {
          fprintf(stdout,"\n  ## WARNING: %" MMG5_PRId " TETRA WITH 4 RIDGES POINTS:"
                  " UNABLE TO COMPUTE ANISO QUALITY.\n",nrid);
        }
      }
    }
  }

  return MMG5_minQualCheck(iel,min,1.);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param ne pointer to the number of used tetra (to fill).
 * \param max pointer to the maximal quality (normalized) to fill.
 * \param avg pointer to the average quality (normalized) to fill.
 * \param min pointer to the minimal quality (normalized) to fill.
 * \param iel pointer to the index of the worst tetra (to fill).
 * \param good number of good elements (to fill).
 * \param med number of elements with a quality greather than 0.5 (to fill).
 * \param his pointer to the mesh histogram (to fill).
 * \param imprim verbosity level
 *
 * Compute the needed quality information in order to print the quality histogram
 * (for a classic storage of the metric at ridges).
 *
 */
void MMG3D_computeInqua(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *ne,double *max,double *avg,
                        double *min,MMG5_int *iel,MMG5_int *good,MMG5_int *med,MMG5_int his[5],int imprim) {
  MMG5_pTetra   pt;
  double        rap;
  MMG5_int      k,ok,nex;
  int           ir;
  static int8_t mmgWarn0 = 0;

  /*compute tet quality*/
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) )   continue;

    if ( met->m ) {
      if ( met->size == 6) {
        pt->qual = MMG5_caltet33_ani(mesh,met,pt);
      }
      else
        pt->qual = MMG5_orcal(mesh,met,k);
    }
    else // -A option
      pt->qual = MMG5_caltet_iso(mesh,met,pt);
  }
  if ( imprim <= 0 ) return;

  (*min)  = 2.0;
  (*max)  = (*avg) = 0.0;
  (*iel)  = 0;
  (*med) = (*good) = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;
    if ( (!mmgWarn0) && (MMG5_orvol(mesh->point,pt->v) < 0.0) ) {
      mmgWarn0 = 1;
      fprintf(stderr,"  ## Warning: %s: at least 1 negative volume\n",
              __func__);
    }
    rap = MMG3D_ALPHAD * pt->qual;
    if ( rap < (*min) ) {
      (*min) = rap;
      (*iel) = ok;
    }
    if ( rap > 0.5 )  (*med)++;
    if ( rap > 0.12 ) (*good)++;
    if ( rap < MMG3D_BADKAL )  mesh->info.badkal = 1;
    (*avg) += rap;
    (*max)  = MG_MAX((*max),rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

  (*ne) = mesh->ne-nex;

  return;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \return 0 if the worst element has a nul quality, 1 otherwise.
 *
 * Print histogram of mesh qualities for classic storage of metric at ridges.
 *
 */
int MMG3D_inqua(MMG5_pMesh mesh,MMG5_pSol met) {
  double      rapmin,rapmax,rapavg;
  int         k;
  MMG5_int    med,good,iel,ne,his[5];

  ne = iel = good = med = 0;
  for ( k=0; k<5; ++k ) {
    his[k] = 0;
  }
  rapmax = 0.;
  rapmin = 2.;
  rapavg = 2.;

  if( mesh->info.optimLES ) {
    MMG3D_computeLESqua(mesh,met,&ne,&rapmax,&rapavg,&rapmin,&iel,&good,&med,
                        his,mesh->info.imprim);
  }
  else {
    MMG3D_computeInqua(mesh,met,&ne,&rapmax,&rapavg,&rapmin,&iel,&good,&med,
                       his,mesh->info.imprim);
  }

  if ( mesh->info.imprim <= 0 )
    return 1;

  return MMG3D_displayQualHisto(ne,rapmax,rapavg,rapmin,
                                iel,good,med,his,0,mesh->info.optimLES,
                                mesh->info.imprim);
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param ne pointer to the number of used tetra (to fill).
 * \param max pointer to the maximal quality (normalized) to fill.
 * \param avg pointer to the average quality (normalized) to fill.
 * \param min pointer to the minimal quality (normalized) to fill.
 * \param iel pointer to the index of the worst tetra (to fill).
 * \param good number of good elements (to fill).
 * \param med number of elements with a quality greather than 0.5 (to fill).
 * \param his pointer to the mesh histogram (to fill).
 * \param nrid number of tetra with 4 ridge points if we want to warn the user
 *             to fill.
 * \param imprim verbosity level
 *
 * Compute the needed quality information in order to print the quality histogram
 * (for special storage of the metric at ridges).
 *
 */
void MMG3D_computeOutqua(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *ne,double *max,double *avg,
                         double *min,MMG5_int *iel,MMG5_int *good,MMG5_int *med,MMG5_int his[5],
                         MMG5_int *nrid,int imprim) {
  MMG5_pTetra   pt;
  MMG5_pPoint   ppt;
  double        rap;
  int           i,ir,n;
  MMG5_int      k,ok,nex;
  static int8_t mmgWarn0 = 0;

  /*compute tet quality*/
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) )   continue;
    pt->qual = MMG5_orcal(mesh,met,k);
  }

  if ( imprim <= 0 )
    return;

  (*min)  = 2.0;
  (*max)  = (*avg) = 0.0;
  (*iel)  = 0;
  (*med) = (*good) = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = (*nrid) = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;
    if ( (!mmgWarn0) && (MMG5_orvol(mesh->point,pt->v) < 0.0) ) {
      mmgWarn0 = 1;
      fprintf(stderr,"  ## Warning: %s: at least 1 negative volume\n",
              __func__);
    }

    /* Count the number of tets with only ridge metric if special metric storage
     * at ridge. */
    if ( mesh->info.metRidTyp==1 ) {
      n = 0;
      for(i=0 ; i<4 ; i++) {
        ppt = &mesh->point[pt->v[i]];
        if(!(MG_SIN(ppt->tag) || MG_NOM & ppt->tag) && (ppt->tag & MG_GEO)) continue;
        n++;
      }
      if(!n) {
        (*nrid)++;
        continue;
      }
    }

    rap = MMG3D_ALPHAD * pt->qual;
    if ( rap < (*min) ) {
      (*min) = rap;
      (*iel) = ok;
    }
    if ( rap > 0.5 )  (*med)++;
    if ( rap > 0.12 ) (*good)++;
    if ( rap < MMG3D_BADKAL )  mesh->info.badkal = 1;
    (*avg) += rap;
    (*max)  = MG_MAX((*max),rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

  (*ne) = mesh->ne-nex;

  return;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 *
 * \return 0 if the worst element has a nul quality, 1 otherwise.
 *
 * Print histogram of mesh qualities for special storage of metric at ridges.
 *
 */
int MMG3D_outqua(MMG5_pMesh mesh,MMG5_pSol met) {
  double      rapmin,rapmax,rapavg;
  int         k;
  MMG5_int    his[5];
  MMG5_int    med,good,iel,ne,nrid;

  nrid = ne = iel = good = med = 0;
  for ( k=0; k<5; ++k ) {
    his[k] = 0;
  }
  rapmax = 0.;
  rapmin = 2.;
  rapavg = 2.;

  if( mesh->info.optimLES ) {
    MMG3D_computeLESqua(mesh,met,&ne,&rapmax,&rapavg,&rapmin,&iel,&good,&med,
                        his,mesh->info.imprim);
  }
  else {
    MMG3D_computeOutqua(mesh,met,&ne,&rapmax,&rapavg,&rapmin,&iel,&good,&med,
                        his,&nrid,mesh->info.imprim);
  }

  if ( mesh->info.imprim <= 0 )
    return 1;

  return MMG3D_displayQualHisto(ne,rapmax,rapavg,rapmin,
                                iel,good,med,his,nrid,mesh->info.optimLES,
                                mesh->info.imprim);
}

/**
 * \param mesh pointer to the mesh.
 * \param sol, pointer to the sol structure.
 * \param weightelt put weight on elts.
 * \param npcible estimation of the final number of nodes/
 *
 * \return 0 if fail, 1 if success
 *
 * Approximation of the final number of vertex.
 *
 * \warning  call MMG3D_hashTetra(mesh,1) or analysis before using
 *
 */
int MMG5_countelt(MMG5_pMesh mesh,MMG5_pSol sol, double *weightelt, long *npcible) {
  MMG5_pTetra pt;
  double      len;
  int         ia,ipa,ipb,lon,l;
  //int   npbdry;
  int         lenint,loc,nedel,longen;
  double      dned,dnface,dnint/*,dnins*/,w,lenavg,lent[6];
  double      dnpdel,dnadd,leninv,dnaddloc,dnpdelloc;
  int         ddebug=0,ib,nv;
  int64_t     list[MMG3D_LMAX];
  MMG5_int    *pdel,k;
  long        nptot;
  //FILE *inm;

  pdel = (MMG5_int*) calloc(mesh->np+1,sizeof(MMG5_int));
  nptot = (long) mesh->np;

  // svg des poids
  // npbdry = 0;
  // inm = fopen("poid.sol","w");
  // fprintf(inm,"MeshVersionFormatted 2\n Dimension 3 \n SolAtTetrahedra \n %" MMG5_PRId "\n 1 1 \n",mesh->ne);

  // substraction of the half of the number of bdry vertex to avoid the surestimation due of the interface
  // for (k=1; k<=mesh->np; k++) {
  //   if(mesh->point[k].tag & MG_BDY) npbdry++;
  // }
  // nptot -= 0.5*npbdry;

  dnadd = dnpdel = 0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /*longueur moyenne*/
    lenavg = 0;
    nv = 6;
    for(ib=0 ; ib<6 ; ib++) {
      ipa = MMG5_iare[ib][0];
      ipb = MMG5_iare[ib][1];
      lent[ib] = MMG5_lenedg(mesh,sol,ib,pt);
      if ( lent[ib]==0 ) nv--;
      lenavg+=lent[ib];
    }
    if ( nv )
      lenavg /= (double)nv;
    else
      lenavg = 1; // Unable to treat this element

    w = 0;
    if(weightelt)
      weightelt[k] = 0;
    nedel = 0;

    for (ia=0; ia<6; ia++) {
      int8_t isbdy;
      longen = MMG5_coquil(mesh,k,ia,list,&isbdy);
      lon = longen/2;

      if ( lon<=0 ) {
        MMG5_SAFE_FREE(pdel);
        return 0;
      }
      /* if ( isbdy )  { */
       /*    continue; */
      /* } */
      //assert(!(longen%2));
      for (l=1; l<lon; l++)
        if ( list[l] < 6*k )  break;

      if ( l < lon )  {
        loc = 1;
        //continue;
      } else {
        loc = 0;
      }

      dnaddloc = 0;
      //dnpdelloc = 0;

      len = lent[ia];

      if(ddebug) printf("len %e\n",len);
      if(len > 3) {
        loc = 0; //count all edges and divide by lon
        len = lenavg; //if very long edge, take the mean
        lenint = ((int) len);
        if(fabs(lenint -len) > 0.5) lenint++;
        //POURQUOI SURESTIMER ???lenint++;
        //nb de point a inserer sur l'arete : longueur - 1
        dned = lenint - 1;
        //nb de point a l'interieur de la face si toutes les aretes sont coupees le meme nb de fois
        dnface = (lenint+2)*(lenint+1) / 2. - 3 - 3*dned;
        //nb de point a l'interieur du tetra si toutes les aretes sont coupees le meme nb de fois
        dnint = (lenint+3)*(lenint+2)*(lenint+1) / 6. - 4 - 4*dnface - 6*dned;
        //nb de point a inserer pour cette arete de ce tetra : on divise par lon
        //dnins = dned*(1./lon) + (dnface/3. + dnint/6.);//(dnface/12. + dnint/6.);
        //if(!isbdry) {
        //nb points sur l'arete +
        //lon*(2/3 nb point sur la face (ie 1/3 de face et 2 faces adj a l'arete) + 1/6 nb de point interne)
        dnaddloc = dned + lon*(2*dnface/3. + dnint/6.);
        //} else {
        //  dnaddloc = 0.5*(dned + lon*(2*dnface/3. + dnint/6.));
        //}
        dnaddloc *= 1./lon;
        if(!loc) {
          /*on ne compte les points internes que pour les (tres) bons tetras*/
          if( MMG3D_ALPHAD * pt->qual < 0.5) {
            if(MMG3D_ALPHAD * pt->qual >= 1./5.)
              dnaddloc = dned / lon + 2*dnface/3.;
            else
              dnaddloc = dned / lon ;
          }
          //rajout de 30% parce que 1) on vise des longueurs de 0.7 et
          //2) on ne tient pas compte du fait qu'on divise tjs par 2 dans la generation
          if( (MMG3D_ALPHAD*pt->qual <= 1./50.) )
            dnaddloc = 0;
          else  if((MMG3D_ALPHAD*pt->qual <= 1./10.) )
            dnaddloc =  0.2*dnaddloc;
          else if((len > 10) && (MMG3D_ALPHAD*pt->qual >= 1./1.5) ) //on sous-estime uniquement pour les tres bons
            dnaddloc = dnaddloc*0.3 + dnaddloc;
          else if(len < 6 && len>3)
            dnaddloc = 0.7*dnaddloc;


          dnadd += dnaddloc;
        }
      } else if(len > 2.8) {
        //if(!isbdry) {
        dnaddloc = 2.;
        //} else {
        //  dnaddloc = 1;
        //}
        if(!loc){
          // if(!isbdry) {
          dnadd += 2.;
          // } else {
          //   dnadd++;
          // }
        }
        //dnins = 2;
      } else if(len > 1.41) {
        //if(!isbdry)
        dnaddloc = 1;
        if(!loc) {
          //if(!isbdry)
          dnadd += 1.;
        }
        //dnins = 1;
      } else if(len < 0.6) {
        nedel = 1;

        leninv = 1./len;
        if(pt->v[ipa]<pt->v[ipb]) {
          if(!pdel[pt->v[ipa]]) {
            // if(!isbdry) {
            dnpdelloc = (leninv - 1.)/leninv;
            // } else {
            //   dnpdelloc = 0.5*(leninv - 1.)/leninv;
            // }
            if(!loc) {
              dnpdel+=dnpdelloc;
              pdel[pt->v[ipa]]=1;
            }
          } else if(!pdel[pt->v[ipb]]) {
            // if(!isbdry) {
            dnpdelloc = (leninv - 1.)/leninv;
            // } else {
            //   dnpdelloc = 0.5*(leninv - 1.)/leninv;
            // }
            if(!loc) {
              dnpdel +=dnpdelloc;
              pdel[pt->v[ipb]]=1;
            }
          }
        } else {
          if(!pdel[pt->v[ipb]]) {
            // if(!isbdry) {
            dnpdelloc = (leninv - 1.)/leninv;
            // } else {
            //   dnpdelloc = 0.5*(leninv - 1.)/leninv;
            // }
            if(!loc) {
              dnpdel+=dnpdelloc;
              pdel[pt->v[ipb]]=1;
            }
          } else if(!pdel[pt->v[ipa]]) {
            // if(!isbdry) {
            dnpdelloc = (leninv - 1.)/leninv;
            // } else {
            //  dnpdelloc = 0.5*(leninv - 1.)/leninv;
            // }
            if(!loc) {
              dnpdel+=dnpdelloc;
              pdel[pt->v[ipa]]=1;
            }
          }
        }
        //ndel++;
      }

      //pour cette arete de ce tetra :
      //PHASE 1 = dnaddloc + nedel (on compte un si arete trop petite)
      //PHASE 2 = dnaddloc
      if(ddebug) printf("on ajoute %e\n",dnaddloc);
      w += (2*dnaddloc);//1./lon*(2*dnaddloc + dnpdelloc);

    }/*for ia*/
    if(ddebug) printf("on soustrait %d\n",nedel);

    w += 0.5*nedel;

    //si l'elt ne doit pas etre ni splitte ni collapse est ce qu'on le compte ???
    //if(w==0) w+=1;

    //fprintf(inm,"%e\n",w);
    if(weightelt)
      weightelt[k] = 10*w;
  } /*For k*/


  nptot += (long) dnadd - (long) dnpdel;
  *npcible = nptot;
  fprintf(stdout,"  ** ESTIMATION OF THE FINAL NUMBER OF NODES : %ld   \n",nptot);
  if(mesh->info.imprim > 6)
    fprintf(stdout,"  **  %lf ADD DEL %lf\n",dnadd,dnpdel);

  free(pdel);

  //fclose(inm);
  return 1;
}
