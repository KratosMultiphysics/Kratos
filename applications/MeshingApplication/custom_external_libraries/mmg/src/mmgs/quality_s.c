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
 * \file mmgs/quality_s.c
 * \brief Functions to compute elements quality and edge lengths.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmgs_private.h"
#include "inlined_functions_private.h"
#include "mmgsexterns_private.h"
#include "mmgexterns_private.h"

extern int8_t ddb;

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param iel element index
 * \return 0 if fail, -1 if orientation is reversed with regards to orientation
 * of vertices, the computed quality otherwise.
 *
 * Quality function identic to caltri_ani but puts a sign according to deviation
 * to normal to vertices.
 *
 */
double caleltsig_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int iel) {
  MMG5_pTria    pt;
  MMG5_pPoint   pa,pb,pc;
  double        ps1,ps2,abx,aby,abz,acx,acy,acz,dd,rap,anisurf;
  double        n[3],pv[3],l[3],*ncomp,*a,*b,*c;
  MMG5_int      ia,ib,ic;

  pt = &mesh->tria[iel];
  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];

  pa = &mesh->point[ia];
  pb = &mesh->point[ib];
  pc = &mesh->point[ic];

  a = &pa->c[0];
  b = &pb->c[0];
  c = &pc->c[0];

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  pv[0] = aby*acz - abz*acy;
  pv[1] = abz*acx - abx*acz;
  pv[2] = abx*acy - aby*acx;

  dd   = pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2];
  if ( dd < MMG5_EPSD2 )  return 0.0;
  dd = 1.0 / sqrt(dd);

  // If one of the triangle vertex is not REF or GEO, it contains the normal at
  // the C1 surface.
  if ( !MG_EDG(pa->tag) ) {
    memcpy(n,&pa->n[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= dd;
  }
  else if ( !MG_EDG(pb->tag) ) {
    memcpy(n,&pb->n[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= dd;
  }
  else if ( !MG_EDG(pc->tag) ) {
    memcpy(n,&pc->n[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= dd;
  }
  else {
    // We must find the normal at the surface elsewhere. Arbitrary, we take it
    // at point pa.
    memcpy(n,&mesh->xpoint[pa->xp].n1[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= dd;

    if ( (pa->tag & MG_GEO) ) {
      ncomp = &mesh->xpoint[pa->xp].n2[0];
      ps2 = ncomp[0]*pv[0]+ncomp[1]*pv[1]+ncomp[2]*pv[2];
      ps2 *= dd;
      if ( fabs(1.0-fabs(ps1)) > fabs(1.0-fabs(ps2)) ) {
        memcpy(n,ncomp,3*sizeof(double));
        ps1 = ps2;
      }
    }
  }

  /* if orientation is reversed with regards to orientation of vertices */
  if ( ps1 < 0.0 )  return -1.0;

  anisurf = MMG5_surftri_ani(mesh,met,pt);
  if ( anisurf == 0.0 )  return -1.0;

  l[0] = MMG5_lenSurfEdg_ani(mesh,met,ib,ic,( pt->tag[0] & MG_GEO ));
  l[1] = MMG5_lenSurfEdg_ani(mesh,met,ia,ic,( pt->tag[1] & MG_GEO ));
  l[2] = MMG5_lenSurfEdg_ani(mesh,met,ia,ib,( pt->tag[2] & MG_GEO ));

  rap = l[0]*l[0] + l[1]*l[1] + l[2]*l[2];
  if ( rap < MMG5_EPSD2 )  return 0.0;
  return anisurf / rap;
}

/* Same quality function but puts a sign according to deviation to normal to vertices */
double caleltsig_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int iel) {
  MMG5_pTria     pt;
  MMG5_pPoint    pa,pb,pc;
  double         *a,*b,*c,cal,abx,aby,abz,acx,acy,acz,bcx,bcy,bcz,rap;
  double         n[3],*ncomp,pv[3],ps1,ps2,sqcal,invsqcal;
  MMG5_int       ia,ib,ic;

  pt = &mesh->tria[iel];
  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];

  pa = &mesh->point[ia];
  pb = &mesh->point[ib];
  pc = &mesh->point[ic];

  a = &pa->c[0];
  b = &pb->c[0];
  c = &pc->c[0];

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];

  pv[0] = aby*acz - abz*acy;
  pv[1] = abz*acx - abx*acz;
  pv[2] = abx*acy - aby*acx;

  cal   = pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2];
  sqcal = sqrt(cal);

  if ( sqcal < MMG5_EPSD2 )  return 0.0;
  invsqcal = 1.0 / sqcal;

  if ( !MG_EDG(pa->tag) ) {
    memcpy(n,&pa->n[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= invsqcal;
  }
  else if ( !MG_EDG(pb->tag) ) {
    memcpy(n,&pb->n[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= invsqcal;
  }
  else if ( !MG_EDG(pc->tag) ) {
    memcpy(n,&pc->n[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= invsqcal;
  }
  else {
    memcpy(n,&mesh->xpoint[pa->xp].n1[0],3*sizeof(double));
    ps1 = n[0]*pv[0]+n[1]*pv[1]+n[2]*pv[2];
    ps1 *= invsqcal;

    if ( (pa->tag & MG_GEO) ) {
      ncomp = &mesh->xpoint[pa->xp].n2[0];
      ps2 = ncomp[0]*pv[0]+ncomp[1]*pv[1]+ncomp[2]*pv[2];
      ps2 *= invsqcal;
      if ( fabs(1.0-fabs(ps1)) > fabs(1.0-fabs(ps2)) ) {
        memcpy(n,ncomp,3*sizeof(double));
        ps1 = ps2;
      }
    }
  }

  /* if orientation is reversed with regards to orientation of vertex */
  if ( ps1 < 0.0 )  return -1.0;
  if ( cal > MMG5_EPSD2 ) {
    /* qual = 2.*surf / length */
    rap  = abx*abx + aby*aby + abz*abz;
    rap += acx*acx + acy*acy + acz*acz;
    rap += bcx*bcx + bcy*bcy + bcz*bcz;
    if ( rap > MMG5_EPSD2 )
      return sqrt(cal) / rap;
    else
      return 0.0;
  }
  else
    return 0.0;
}


/* coordinates of the center of incircle of p0p1p2 and its 'size' */
inline double incircle(MMG5_pPoint p0,MMG5_pPoint p1,MMG5_pPoint p2,double *o) {
  double   dd,r,rr;

  dd = 1.0 / 3.0;
  o[0] = dd * (p0->c[0] + p1->c[0] + p2->c[0]);
  o[1] = dd * (p0->c[1] + p1->c[1] + p2->c[1]);
  o[2] = dd * (p0->c[2] + p1->c[2] + p2->c[2]);

  rr = sqrt((p0->c[0]-o[0])*(p0->c[0]-o[0]) + (p0->c[1]-o[1])*(p0->c[1]-o[1]) \
            + (p0->c[2]-o[2])*(p0->c[2]-o[2]));

  r = sqrt((p1->c[0]-o[0])*(p1->c[0]-o[0]) + (p1->c[1]-o[1])*(p1->c[1]-o[1]) \
           + (p1->c[2]-o[2])*(p1->c[2]-o[2]));
  rr = MG_MAX(rr,r);

  r = sqrt((p2->c[0]-o[0])*(p2->c[0]-o[0]) + (p2->c[1]-o[1])*(p2->c[1]-o[1]) \
           + (p2->c[2]-o[2])*(p2->c[2]-o[2]));
  rr = MG_MAX(rr,r);

  return rr;
}

inline double diamelt(MMG5_pPoint p0,MMG5_pPoint p1,MMG5_pPoint p2) {
  double  di,dd;

  di = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0])
    + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])
    + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);

  dd = (p2->c[0]-p0->c[0])*(p2->c[0]-p0->c[0])
    + (p2->c[1]-p0->c[1])*(p2->c[1]-p0->c[1])
    + (p2->c[2]-p0->c[2])*(p2->c[2]-p0->c[2]);
  di = MG_MAX(di,dd);

  dd = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0])
    + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1])
    + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
  di = MG_MAX(di,dd);

  return di;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage,
 * 1 for special storage.
 * \return 0 if fail, 1 otherwise.
 *
 * Compute sizes of edges of the mesh, and displays histo.
 *
 */
int MMGS_prilen(MMG5_pMesh mesh, MMG5_pSol met, int metRidTyp) {
  MMG5_pTria      pt;
  MMG5_Hash       hash;
  double          len,avlen,lmin,lmax;
  MMG5_int        ned,hl[9],nullEdge;
  MMG5_int        k,np,nq,amin,bmin,amax,bmax;
  int8_t          ia,i0,i1,i;
  static double   bd[9]= {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};

  memset(hl,0,9*sizeof(MMG5_int));
  ned = 0;
  avlen = 0.0;
  lmax = 0.0;
  lmin = 1.e30;
  amin = amax = bmin = bmax = 0;
  nullEdge = 0;

  /* Hash all edges in the mesh */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return 0;

  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<3; ia++) {
      i0 = MMG5_iprv2[ia];
      i1 = MMG5_inxt2[ia];
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
  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<3; ia++) {
      i0 = MMG5_iprv2[ia];
      i1 = MMG5_inxt2[ia];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash */
      MMG5_hashGet(&hash,np,nq);

      if ( (!metRidTyp) && met->m && met->size>1 ) {
        len = MMG5_lenSurfEdg33_ani(mesh,met,np,nq,(pt->tag[ia] & MG_GEO));
      }
      else
        len = MMG5_lenSurfEdg(mesh,met,np,nq,(pt->tag[ia] & MG_GEO));

      if ( !len ) {
        ++nullEdge;
      }
      else {
        ned ++;
        avlen += len;

        if( len < lmin ) {
          lmin = len;
          amin = np;
          bmin = nq;
        }

        if ( len > lmax ) {
          lmax = len;
          amax = np;
          bmax = nq;
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

  /* Display histogram */
  MMG5_displayLengthHisto(mesh, ned, &avlen, amin, bmin, lmin,
                           amax, bmax, lmax, nullEdge, &bd[0], &hl[0],0);

  MMG5_DEL_MEM(mesh,hash.item);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if the worst element has a nul quality, 1 otherwise.
 *
 * Print histogram of mesh qualities for classical storage of ridges
 * metrics (so before the the MMG5_defsiz function call).
 *
 */
int MMGS_inqua(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  double        rap,rapmin,rapmax,rapavg,med;
  MMG5_int      his[5],k,iel,nex,ok;
  int           i,ir,imax;

  rapmin  = 1.0;
  rapmax  = 0.0;
  rapavg  = med = 0.0;
  iel     = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;

    if ( met->m && (met->size == 6) ) {
      rap = MMGS_ALPHAD * MMG5_caltri33_ani(mesh,met,pt);
    }
    else {
      rap = MMGS_ALPHAD * MMG5_caltri_iso (mesh,NULL,pt);
    }

    if ( rap < rapmin ) {
      rapmin = rap;
      iel    = ok;
    }
    if ( rap > 0.5 )  med++;
    if ( rap < MMGS_BADKAL )  mesh->info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

  fprintf(stdout,"\n  -- MESH QUALITY   %" MMG5_PRId "\n",mesh->nt - nex);
  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%" MMG5_PRId ")\n",
          rapmax,rapavg / (mesh->nt-nex),rapmin,iel);

  if ( mesh->info.imprim >= 3 ){

    /* print histo */
    fprintf(stdout,"     HISTOGRAMM:  %6.2f %% > 0.5\n",100.0*(med/(float)(mesh->nt-nex)));
    imax = MG_MIN(4,(int)(5.*rapmax));
    for (i=imax; i>=(int)(5*rapmin); i--) {
      fprintf(stdout,"     %5.1f < Q < %5.1f   %7"MMG5_PRId"   %6.2f %%\n",
              i/5.,i/5.+0.2,his[i],100.*(his[i]/(float)(mesh->nt-nex)));
    }
  }

  return  MMG5_minQualCheck(iel,rapmin,1.);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if the worst element has a nul quality, 1 otherwise.
 *
 * Print histogram of mesh qualities for special storage of ridges metrics
 * (after the MMG5_defsiz function call).
 *
 */
int MMGS_outqua(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  double        rap,rapmin,rapmax,rapavg,med;
  MMG5_int      nex,his[5],k,iel,ok;
  int           i,ir,imax;

  if ( mesh->info.imprim <= 0 ) return 1;

  rapmin  = 1.0;
  rapmax  = 0.0;
  rapavg  = med = 0.0;
  iel     = 0;

  for (k=0; k<5; k++)  his[k] = 0;

  nex = ok = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;

    rap = MMGS_ALPHAD * MMG5_calelt(mesh,met,pt);

    if ( rap < rapmin ) {
      rapmin = rap;
      iel    = ok;
    }
    if ( rap > 0.5 )  med++;
    if ( rap < MMGS_BADKAL )  mesh->info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

  fprintf(stdout,"\n  -- MESH QUALITY   %" MMG5_PRId "\n",mesh->nt - nex);
  fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%" MMG5_PRId ")\n",
          rapmax,rapavg / (mesh->nt-nex),rapmin,iel);

  if ( mesh->info.imprim >= 3 ){
    /* print histo */
    fprintf(stdout,"     HISTOGRAMM:  %6.2f %% > 0.5\n",100.0*(med/(float)(mesh->nt-nex)));
    imax = MG_MIN(4,(int)(5.*rapmax));
    for (i=imax; i>=(int)(5*rapmin); i--) {
      fprintf(stdout,"     %5.1f < Q < %5.1f   %7"MMG5_PRId"   %6.2f %%\n",
              i/5.,i/5.+0.2,his[i],100.*(his[i]/(float)(mesh->nt-nex)));
    }
  }

  return  MMG5_minQualCheck(iel,rapmin,1.);
}

#define COS145   -0.81915204428899

/* return 0: triangle ok, 1: needle, 2: obtuse; ia: edge problem */
int8_t typelt(MMG5_pPoint p[3],int8_t *ia) {
  double   h1,h2,h3,hmi,hma,ux,uy,uz,vx,vy,vz,wx,wy,wz,dd;

  ux = p[1]->c[0] - p[0]->c[0];
  uy = p[1]->c[1] - p[0]->c[1];
  uz = p[1]->c[2] - p[0]->c[2];
  h1 = ux*ux + uy*uy + uz*uz;

  vx = p[2]->c[0] - p[0]->c[0];
  vy = p[2]->c[1] - p[0]->c[1];
  vz = p[2]->c[2] - p[0]->c[2];
  h2 = vx*vx + vy*vy + vz*vz;
  hmi = h1;
  hma = h2;
  *ia = 2;
  if ( h1 > h2 ) {
    hmi = h2;
    hma = h1;
    *ia = 1;
  }
  wx = p[2]->c[0] - p[1]->c[0];
  wy = p[2]->c[1] - p[1]->c[1];
  wz = p[2]->c[2] - p[1]->c[2];
  h3 = wx*wx + wy*wy + wz*wz;
  if ( h3 < hmi ) {
    hmi = h3;
    *ia = 0;
  }
  else if ( h3 > hma )
    hma = h3;

  /* needle */
  if ( hmi < 0.01 * hma )  return 1;

  /* check obtuse angle */
  dd = (ux*vx + uy*vy + uz*vz) / sqrt(h1*h2);
  if ( dd < COS145 ) {
    *ia = 0;
    return 2;
  }
  dd = (vx*wx + vy*wy + vz*wz) / sqrt(h2*h3);
  if ( dd < COS145 ) {
    *ia = 2;
    return 2;
  }
  dd = -(ux*wx + uy*wy + uz*wz) / sqrt(h1*h3);
  if ( dd < COS145 ) {
    *ia = 1;
    return 2;
  }

  return 0;
}
