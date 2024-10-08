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
 * \file mmgs/analys_s.c
 * \brief Mesh analysis.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgcommon_private.h"


/**
 * \param mesh pointer toward a MMG5 mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Regularization procedure for derivatives, dual Laplacian
 *
 */
int MMG5_regnor(MMG5_pMesh mesh) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt,p0;
  MMG5_pxPoint  pxp;
  double        *tabl,n[3],*nptr,lm1,lm2,dd,nx,ny,nz,res0,res;
  int           i,it,nit,ilist;
  MMG5_int      k,nn,iel,list[MMG5_LMAX],tlist[MMG5_LMAX],*adja,iad;

  /* assign seed to vertex */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k;
    }
  }

  /* allocate memory for normals */
  MMG5_SAFE_CALLOC(tabl,3*mesh->np+1,double,return 0);

  /* Pointer toward the suitable adjacency array for Mmgs and Mmg3d */
  if ( mesh->adjt ) {
    /* Mmg3d */
    adja = mesh->adjt;
  }
  else {
    /* Mmgs */
    adja = mesh->adja;
  }

  it   = 0;
  nit  = 10;
  res0 = 0.0;
  lm1  = 0.4;
  lm2  = 0.399;
  while ( it++ < nit ) {
    /* step 1: laplacian */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      if ( ppt->tag & MG_CRN || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) continue;

      iel = ppt->s;
      if ( !iel ) continue; // Mmg3d

      pt = &mesh->tria[iel];
      i  = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = MMG5_boulep(mesh,iel,i,adja,list,tlist);

      /* average normal */
      nx = ny = nz = 0.0;
      for (i=1; i<=ilist; i++) {
        p0  = &mesh->point[list[i]];
        if ( p0->tag & MG_CRN || p0->tag & MG_NOM || MG_EDG(p0->tag) ) continue;

        if ( p0->xp ) {
          /* Mmg3d */
          pxp  = &mesh->xpoint[p0->xp];
          nptr = pxp->n1;
        }
        else {
          /* Mmgs */
          nptr = p0->n;
        }

        nx += nptr[0];
        ny += nptr[1];
        nz += nptr[2];
      }

      dd  = nx*nx + ny*ny + nz*nz;
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nx *= dd;
        ny *= dd;
        nz *= dd;
      }

      /* Laplacian */
      if ( ppt->xp ) {
        /* Mmg3d */
        pxp = &mesh->xpoint[ppt->xp];
        nptr = pxp->n1;
      }
      else {
        /* Mmgs */
        nptr = ppt->n;
      }

      iad = 3*(k-1)+1;
      tabl[iad+0] = nptr[0] + lm1 * (nx - nptr[0]);
      tabl[iad+1] = nptr[1] + lm1 * (ny - nptr[1]);
      tabl[iad+2] = nptr[2] + lm1 * (nz - nptr[2]);
    }

    /* step 2: anti-laplacian */
    res = 0;
    nn  = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];

      if ( !MG_VOK(ppt) )  continue;
      if ( ppt->tag & MG_CRN || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) continue;

      iel = ppt->s;
      if ( !iel ) continue; // Mmg3d

      pt = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = MMG5_boulep(mesh,iel,i,adja,list,tlist);

      /* average normal */
      nx = ny = nz = 0.0;
      for (i=1; i<=ilist; i++) {
        iad = 3*(list[i]-1) + 1;
        nx += tabl[iad+0];
        ny += tabl[iad+1];
        nz += tabl[iad+2];
      }
      dd  = nx*nx + ny*ny + nz*nz;
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nx *= dd;
        ny *= dd;
        nz *= dd;
      }

      /* antiLaplacian */
      iad = 3*(k-1)+1;
      n[0] = tabl[iad+0] - lm2 * (nx - tabl[iad+0]);
      n[1] = tabl[iad+1] - lm2 * (ny - tabl[iad+1]);
      n[2] = tabl[iad+2] - lm2 * (nz - tabl[iad+2]);
      nn++;

      if ( ppt->xp ) {
        /* Mmg3d */
        pxp = &mesh->xpoint[ppt->xp];
        nptr = pxp->n1;
      }
      else {
        /* Mmgs */
        nptr = ppt->n;
      }
      res += (nptr[0]-n[0])*(nptr[0]*n[0]) + (nptr[1]-n[1])*(nptr[1]*n[1]) + (nptr[2]-n[2])*(nptr[2]*n[2]);

      nptr[0] = n[0];
      nptr[1] = n[1];
      nptr[2] = n[2];

    }

    if ( it == 1 )  res0 = res;
    if ( res0 > MMG5_EPSD )  res  = res / res0;
    if ( mesh->info.imprim < -1 || mesh->info.ddebug ) {
      fprintf(stdout,"     iter %5d  res %.3E\r",it,res);
      fflush(stdout);
    }
    if ( it > 1 && res < MMG5_EPS )  break;
  }

  /* reset the ppt->s tag */
  for (k=1; k<=mesh->np; ++k) {
    mesh->point[k].s = 0;
  }

  if ( mesh->info.imprim < -1 || mesh->info.ddebug )  fprintf(stdout,"\n");

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     %" MMG5_PRId " normals regularized: %.3e\n",nn,res);

  MMG5_SAFE_FREE(tabl);
  return 1;
}
