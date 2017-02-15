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
 * \file mmg2d/optlap_2d.c
 * \brief Functions to optimize with a laplacian/antilaplacian.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"

/**
 * optimisation by laplacian-antilaplacian
 * \warning UNUSED FUNCTION : change memory allocation to use it
 */

int MMG2_optlap(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria   pt;
  MMG5_pPoint  ppt,ppta,pptb;
  int     k,it,maxtou,i,ipa,ipb,iadr;
  double  omega,mu,*cnew,*cold,*cini,res,res0,*ncount,len;
  int iare[3][2] = {{0,1},{0,2},{1,2}};

  cnew = (double*)malloc(2*(mesh->np+1)*sizeof(double));
  assert(cnew);
  cold = (double*)malloc(2*(mesh->np+1)*sizeof(double));
  assert(cold);
  cini = (double*)malloc(2*(mesh->np+1)*sizeof(double));
  assert(cini);
  ncount = (double*)malloc((mesh->np+1)*sizeof(double));
  assert(ncount);

  omega  = 0.1;
  mu     = 0.;
  res0   = 0;
  maxtou = 20000;
  it     = 1;
  for (k=1 ; k<=mesh->np ; k++) {
    ppt = &mesh->point[k];
    if ( !M_VOK(ppt) ) continue;
    iadr = 2*k;
    cini[iadr + 0] = ppt->c[0];
    cini[iadr + 1] = ppt->c[1];
    cnew[iadr + 0] = 0;
    cnew[iadr + 1] = 0;
    ncount[k] = 0;
  }

  do {
    res    = 0;

    /*laplacian*/
    for (k=1 ; k<=mesh->nt ; k++) {
      pt = &mesh->tria[k];
      if ( !M_EOK(pt) ) continue;

      for (i=0 ; i<3 ; i++) {
        ipa  = pt->v[iare[i][0]];
        ppta = &mesh->point[ipa];

        ipb  = pt->v[iare[i][1]];
        pptb = &mesh->point[ipb];

        len = 1;
        /*sqrt((cini[2*ipa + 0]-cini[2*ipb + 0])*(cini[2*ipa + 0]-cini[2*ipb + 0]) +
          (cini[2*ipa + 1]-cini[2*ipb + 1])*(cini[2*ipa + 1]-cini[2*ipb + 1]));
        */len = 1./len;
        if(ppta->tag & M_BDRY) {
          iadr = 2*ipa;
          cnew[iadr + 0] = 0;
          cnew[iadr + 1] = 0;
          ncount[ipa] += 1;
        } else {
          iadr = 2*ipa;
          cnew[iadr + 0] += pptb->c[0]*len;
          cnew[iadr + 1] += pptb->c[1]*len;
          ncount[ipa]    += len;
        }
        if(pptb->tag & M_BDRY) {
          iadr = 2*ipb;
          cnew[iadr + 0] = 0;
          cnew[iadr + 1] = 0;
          ncount[ipb] += 1;
        } else {
          iadr = 2*ipb;
          cnew[iadr + 0] += ppta->c[0]*len;
          cnew[iadr + 1] += ppta->c[1]*len;
          ncount[ipb]    += len;
        }
      }
    } /* end for k*/

    for (k=1 ; k<=mesh->np ; k++) {
      ppt = &mesh->point[k];
      if ( !M_VOK(ppt) ) continue;
      iadr = 2*k;
      cold[iadr + 0] = ppt->c[0];
      cold[iadr + 1] = ppt->c[1];
      if ( ppt->tag & M_BDRY )  continue;

      cnew[iadr + 0] = ppt->c[0] + omega * (cnew[iadr + 0] / ncount[k] - ppt->c[0]);
      cnew[iadr + 1] = ppt->c[1] + omega * (cnew[iadr + 1] / ncount[k] - ppt->c[1]);
      ppt->c[0] = 0;
      ppt->c[1] = 0;
      ncount[k] = 0;
    }

    /*anti-laplacian*/
    for (k=1 ; k<=mesh->nt ; k++) {
      pt = &mesh->tria[k];
      if ( !M_EOK(pt) ) continue;

      for (i=0 ; i<3 ; i++) {
        ipa  = pt->v[iare[i][0]];
        ppta = &mesh->point[ipa];

        ipb  = pt->v[iare[i][1]];
        pptb = &mesh->point[ipb];

        len = 1;
        /*sqrt((cini[2*ipa + 0]-cini[2*ipb + 0])*(cini[2*ipa + 0]-cini[2*ipb + 0]) +
          (cini[2*ipa + 1]-cini[2*ipb + 1])*(cini[2*ipa + 1]-cini[2*ipb + 1]));
        */len = 1./len;

        if(ppta->tag & M_BDRY) {
          iadr = 2*ipa;
          ncount[ipa] += 1;
        } else {
          iadr = 2*ipa;
          ppta->c[0]  += cnew[iadr + 0]*len;
          ppta->c[1]  += cnew[iadr + 1]*len;
          ncount[ipa] += len;
        }
        if(pptb->tag & M_BDRY) {
          iadr = 2*ipb;
          ncount[ipb] += 1;
        } else {
          iadr = 2*ipb;
          pptb->c[0]  += cnew[iadr + 0]*len;
          pptb->c[1]  += cnew[iadr + 1]*len;
          ncount[ipb] += len;
        }
      }
    } /* end for k*/

    for (k=1 ; k<=mesh->np ; k++) {
      ppt = &mesh->point[k];
      if ( !M_VOK(ppt) ) continue;
      if ( ppt->tag & M_BDRY )  continue;
      iadr = 2*k;
      ppt->c[0] = cnew[iadr + 0] - mu * (ppt->c[0] / ncount[k] - cnew[iadr + 0]);
      ppt->c[1] = cnew[iadr + 1] - mu * (ppt->c[1] / ncount[k] - cnew[iadr + 1]);
      ncount[k] = 0;
      cnew[iadr + 0] = 0;
      cnew[iadr + 1] = 0;
      res+=(ppt->c[0]-cold[iadr + 0])*(ppt->c[0]-cold[iadr + 0])
        + (ppt->c[1]-cold[iadr + 1])*(ppt->c[1]-cold[iadr + 1]);
    }

    if (it==1) res0=res;
    if ( res0 > 1e-10 ) {
      if ( mesh->info.imprim > 4 || mesh->info.ddebug )
        fprintf(stdout,"iteration : %d, residu = %e \r",it,res/res0);
    }

  } while((res0 > 1e-10) && (res/res0 > 1e-10) && (it++ < maxtou));

  fprintf(stdout,"\n");

  free(cold);
  free(cnew);
  free(ncount);
  return(1);
}
