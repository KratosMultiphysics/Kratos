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
#include "mmg2d.h"

#define LFILT  0.7


/* create bucket structure and store initial vertices */
pBucket MMG2_newBucket(MMG5_pMesh mesh,int nmax) {
  MMG5_pPoint   ppt;
  pBucket  bucket;
  double   dd;
  int      k,ic,ii,jj;

  /* memory alloc */
  _MMG5_SAFE_CALLOC(bucket,1,Bucket);

  bucket->size = nmax;
  _MMG5_SAFE_CALLOC(bucket->head,nmax*nmax*nmax+1,int);
  _MMG5_SAFE_CALLOC(bucket->link,mesh->npmax+1,int);


  /* insert vertices */
  dd = nmax / (double)PRECI;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !M_VOK(ppt) )  continue;
    ii = M_MAX(0,(int)(dd * ppt->c[0])-1);
    jj = M_MAX(0,(int)(dd * ppt->c[1])-1);
    ic = jj*nmax + ii;

    if ( !bucket->head[ic] )
      bucket->head[ic] = k;
    else {
      bucket->link[k]  = bucket->head[ic];
      bucket->head[ic] = k;
    }
  }

  return(bucket);
}


void MMG2_freeBucket(pBucket bucket) {
  _MMG5_SAFE_FREE(bucket->head);
  _MMG5_SAFE_FREE(bucket->link);
  _MMG5_SAFE_FREE(bucket);
}


/* check and eventually insert vertex */
int buckin_ani(MMG5_pMesh mesh,MMG5_pSol sol,pBucket bucket,int ip) {
  MMG5_pPoint    ppt,pp1;
  double    dd,d2,det,ux,uy,dmi,dx,dy;
  double    *ma,*mb;
  int       i,j,ii,jj,ic,icc,siz,ip1;
  int       iadr,imin,imax,jmin,jmax;

  ppt  = &mesh->point[ip];
  siz  = bucket->size;
  dd   = siz / (double)PRECI;
  iadr = ip*sol->size;
  ma   = &sol->m[iadr];
  dmi  = LFILT*LFILT;

  ii = M_MAX(0,(int)(dd * ppt->c[0])-1);
  jj = M_MAX(0,(int)(dd * ppt->c[1])-1);
  ic = jj*siz + ii;

  /* check current cell */
  if ( bucket->head[ic] ) {
    ip1 = bucket->head[ic];
    pp1 = &mesh->point[ip1];
    ux = pp1->c[0] - ppt->c[0];
    uy = pp1->c[1] - ppt->c[1];
    d2 = ma[0]*ux*ux + ma[2]*uy*uy + 2.0*ma[1]*ux*uy;
    if ( d2 < dmi ) {
      iadr = ip1*sol->size;
      mb = &sol->m[iadr];
      d2 = mb[0]*ux*ux + mb[2]*uy*uy + 2.0*mb[1]*ux*uy;
      if ( d2 < dmi )  return(0);
    }

    while ( bucket->link[ip1] ) {
      ip1 = bucket->link[ip1];
      pp1 = &mesh->point[ip1];
      ux = pp1->c[0] - ppt->c[0];
      uy = pp1->c[1] - ppt->c[1];
      d2 = ma[0]*ux*ux + ma[2]*uy*uy + 2.0*ma[1]*ux*uy;
      if ( d2 < dmi ) {
        iadr = ip1*sol->size;
        mb = &sol->m[iadr];
        d2 = mb[0]*ux*ux + mb[2]*uy*uy + 2.0*mb[1]*ux*uy;
        if ( d2 < dmi )  return(0);
      }
    }
  }

  /* bounding box */
  det = ma[0]*ma[2] - ma[1]*ma[1];
  det = 1.0 / det;
  if ( det < 0.0 )
    return(1);
  else {
    dx = LFILT * sqrt(ma[2] * det) ;
    dy = LFILT * sqrt(ma[1] * det) ;
  }

  imin = (int)(dd * (ppt->c[0]-dx))-1;
  jmin = (int)(dd * (ppt->c[1]-dy))-1;
  imax = (int)(dd * (ppt->c[0]+dx))-1;
  jmax = (int)(dd * (ppt->c[1]+dy))-1;

  imin = M_MAX(0,M_MIN(imin,siz-1));
  imax = M_MIN(siz-1,M_MAX(0,imax));
  jmin = M_MAX(0,M_MIN(jmin,siz-1));
  jmax = M_MIN(siz-1,M_MAX(0,jmax));
  if ( imin == imax && jmin == jmax )  return(1);

  /* explore neighbours */
  for (j=jmin; j<=jmax; j++)
    for (i=imin; i<=imax; i++) {
      icc = j*siz + i;
      ip1 = bucket->head[icc];
      if ( !ip1 )  continue;
      pp1 = &mesh->point[ip1];
      ux = pp1->c[0] - ppt->c[0];
      uy = pp1->c[1] - ppt->c[1];
      d2 = ma[0]*ux*ux + ma[2]*uy*uy + 2.0*ma[1]*ux*uy;
      if ( d2 < dmi ) {
        iadr = ip1*sol->size;
        mb = &sol->m[iadr];
        d2 = mb[0]*ux*ux + mb[2]*uy*uy + 2.0*mb[1]*ux*uy;
        if ( d2 < dmi )  return(0);
      }

      while ( bucket->link[ip1] ) {
        ip1 = bucket->link[ip1];
        pp1 = &mesh->point[ip1];
        ux = pp1->c[0] - ppt->c[0];
        uy = pp1->c[1] - ppt->c[1];
        d2 = ma[0]*ux*ux + ma[2]*uy*uy + 2.0*ma[1]*ux*uy;
        if ( d2 < dmi ) {
          iadr = ip1*sol->size;
          mb = &sol->m[iadr];
          d2 = mb[0]*ux*ux + mb[2]*uy*uy + 2.0*mb[1]*ux*uy;
          if ( d2 < dmi )  return(0);
        }
      }
    }

  return(1);
}


int buckin_iso(MMG5_pMesh mesh,MMG5_pSol sol,pBucket bucket,int ip) {
  MMG5_pPoint    ppt,pp1;
  double    dd,d2,ux,uy,hpi,hp1,hp2;
  int       i,j,ii,jj,ic,icc,siz,ip1;
  int       imin,imax,jmin,jmax;

  ppt = &mesh->point[ip];
  siz = bucket->size;
  dd  = siz / (double)PRECI;
  hpi = LFILT * sol->m[ip];
  hp1 = hpi*hpi;
  ii  = M_MAX(0,(int)(dd * ppt->c[0])-1);
  jj  = M_MAX(0,(int)(dd * ppt->c[1])-1);
  ic  = jj*siz + ii;

  /* check current cell */
  if ( bucket->head[ic] ) {
    ip1 = bucket->head[ic];
    pp1 = &mesh->point[ip1];
    hp2 = LFILT * sol->m[ip1];
    ux = pp1->c[0] - ppt->c[0];
    uy = pp1->c[1] - ppt->c[1];
    d2 = ux*ux + uy*uy;
    if ( d2 < hp1 || d2 < hp2*hp2 )  return(0);

    while ( bucket->link[ip1] ) {
      ip1 = bucket->link[ip1];
      pp1 = &mesh->point[ip1];
      hp2 = LFILT * sol->m[ip1];
      ux = pp1->c[0] - ppt->c[0];
      uy = pp1->c[1] - ppt->c[1];
      d2 = ux*ux + uy*uy;
      if ( d2 < hp1 || d2 < hp2*hp2 )  return(0);
    }
  }

  /* explore neighbors */
  imin = (int)(dd * (ppt->c[0]-hpi))-1;
  jmin = (int)(dd * (ppt->c[1]-hpi))-1;
  imax = (int)(dd * (ppt->c[0]+hpi))-1;
  jmax = (int)(dd * (ppt->c[1]+hpi))-1;

  imin = M_MAX(0,M_MIN(imin,siz-1));
  imax = M_MIN(siz-1,M_MAX(0,imax));
  jmin = M_MAX(0,M_MIN(jmin,siz-1));
  jmax = M_MIN(siz-1,M_MAX(0,jmax));
  if ( imin == imax && jmin == jmax )  return(1);

  for (j=jmin; j<=jmax; j++)
    for (i=imin; i<=imax; i++) {
      icc = j*siz + i;
      ip1 = bucket->head[icc];
      if ( !ip1 )  continue;
      pp1 = &mesh->point[ip1];
      hp2 = LFILT * sol->m[ip1];
      ux = pp1->c[0] - ppt->c[0];
      uy = pp1->c[1] - ppt->c[1];
      d2 = ux*ux + uy*uy;
      if ( d2 < hp1 || d2 < hp2*hp2 )  return(0);

      while ( bucket->link[ip1] ) {
        ip1 = bucket->link[ip1];
        pp1 = &mesh->point[ip1];
        hp2 = LFILT * sol->m[ip1];
        ux = pp1->c[0] - ppt->c[0];
        uy = pp1->c[1] - ppt->c[1];
        d2 = ux*ux + uy*uy;
        if ( d2 < hp1 || d2 < hp2*hp2 )  return(0);
      }
    }

  return(1);
}


int MMG2_addBucket(MMG5_pMesh mesh,pBucket bucket,int ip) {
  MMG5_pPoint    ppt;
  double    dd;
  int       ic,ii,jj,siz;

  ppt = &mesh->point[ip];
  siz = bucket->size;
  dd  = siz / (double)PRECI;
  ii  = M_MAX(0,(int)(dd * ppt->c[0])-1);
  jj  = M_MAX(0,(int)(dd * ppt->c[1])-1);
  ic  = jj*siz + ii;

  /* store new point */
  if ( !bucket->head[ic] ) {
    bucket->head[ic] = ip;
    bucket->link[ip] = 0;
  }
  else {
    bucket->link[ip] = bucket->head[ic];
    bucket->head[ic] = ip;
  }

  return(1);
}


int MMG2_delBucket(MMG5_pMesh mesh,pBucket bucket,int ip) {
  MMG5_pPoint    ppt;
  double    dd;
  int       ic,ii,jj,siz,ip1;

  ppt = &mesh->point[ip];
  siz = bucket->size;
  dd  = siz / (double)PRECI;
  ii  = M_MAX(0,(int)(dd * ppt->c[0])-1);
  jj  = M_MAX(0,(int)(dd * ppt->c[1])-1);
  ic  = jj*siz + ii;

  /* remove vertex from cell */
  if ( bucket->head[ic] ) {
    if ( bucket->head[ic] == ip ) {
      bucket->head[ic] = bucket->link[ip];
      bucket->link[ip] = 0;
    }
    else {
      ip1 = bucket->head[ic];
      while ( ip1 && bucket->link[ip1] != ip ) {
        ip1 = bucket->link[ip1];
      }
      if ( bucket->link[ip1] == ip ) {
        bucket->link[ip1] = bucket->link[ip];
        bucket->link[ip] = 0;
      }
    }
  }

  return(1);
}

