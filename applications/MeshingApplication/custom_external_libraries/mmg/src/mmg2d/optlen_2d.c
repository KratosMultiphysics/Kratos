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
 * \file mmg2d/optlen_2d.c
 * \brief Functions to move the points.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"

#define  HQCOEF    0.9
#define  HCRIT     0.98

int optlen_ani(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppa,ppb;
  pQueue    queue;
  int      *list;
  double    oldc[3],cal,ctg,cx,cy,ux,uy,cpx,cpy,coe,dd,len;
  double    *mb,*mp,*ca,*cb,*qual;
  int       i,j,k,l,iel,lon,nm;
  int       ipa,ipb,nb,nk,npp,iadr,iter,maxtou;
  int nrj;
//  double tmp1,tmp2,tmp3;

  /* queue on quality */
  queue = MMG2_kiuini(mesh,mesh->nt,declic,base - 1);
  assert(queue);

  maxtou = 10;
  nm     = 0;
  npp    = 0;
  nrj = 0;

  _MMG5_SAFE_CALLOC(list,MMG2D_LMAX,int);

  do {
    k = MMG2_kiupop(queue);

    if ( !k )  break;
    npp++;

    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ipa = pt->v[i];
      ppa = &mesh->point[ipa];
      if ( ppa->tag & M_BDRY || ppa->tag & M_REQUIRED || ppa->tag & M_SD)  continue;

      lon   = MMG2_boulep(mesh,k,i,list);
      _MMG5_SAFE_MALLOC(qual,lon+1,double);

      /* optimal point */
      ca   = &ppa->c[0];
      iadr = ipa*sol->size;
      mp   = &sol->m[iadr];
      cx   = 0.0;
      cy   = 0.0;
      nb   = 0;
      cal  = pt->qual;
      for (l=1 ; l<=lon; l++) {
        iel = list[l] / 3;
        nk  = list[l] % 3;
        pt1 = &mesh->tria[iel];
        if ( pt1->qual > cal )  cal = pt1->qual;
        /* tmp1 = 0; */
        /* tmp2 = 0; */
        /* tmp3 = 0; */
        for (j=1; j<3; j++) {
          ipb  = pt1->v[ MMG2_idir[nk + j] ];
          ppb  = &mesh->point[ipb];
          cb   = &ppb->c[0];
          iadr = ipb*sol->size;
          mb   = &sol->m[iadr];

          len = MMG2_length(ca,cb,mp,mb);
          /* optimal point */
          ux  = ppb->c[0] - ppa->c[0];
          uy  = ppb->c[1] - ppa->c[1];
          cx += ppa->c[0] + ux*(1. - 1./len);//ux * len;
          cy += ppa->c[1] + uy*(1. - 1./len);//uy * len;
          nb++;
        }
      }

      //if ( nb < 3 )  continue;
      dd  = 1.0 / (double)nb;
      cpx = cx*dd - ppa->c[0];
      cpy = cy*dd - ppa->c[1];

      /* adjust position */
      coe  = HQCOEF;
      iter = 1;
      if(cal > 10./ALPHA) {
        ctg  = 0.99 * cal;
      } else {
        ctg  = cal * HCRIT;
      }
      memcpy(oldc,ppa->c,3*sizeof(double));
      do {
        ppa->c[0] =/* (1. - coe) * */oldc[0] + coe * cpx;
        ppa->c[1] =/* (1. - coe) * */oldc[1] + coe * cpy;

        for (l=1; l<=lon; l++) {
          iel = list[l] / 3;
          nk  = list[l] % 3;
          pt1 = &mesh->tria[iel];

          cal = MMG2_caltri_in(mesh,sol,pt1);
          if ( cal > ctg )  {
            break;  }
          qual[l] = cal;
        }
        if ( l > lon )  break;
        coe *= 0.5;
      }
      while ( ++iter <= maxtou );
      if ( iter > maxtou ) {
        memcpy(ppa->c,oldc,3*sizeof(double));
        ppa->flag = base - 2;
        nrj++;
        _MMG5_SAFE_FREE(qual);
        continue;
      }

      /* update tria */
      for (l=1; l<=lon; l++) {
        iel = list[l] / 3;
        nk  = list[l] % 3;
        pt1 = &mesh->tria[iel];

        if ( (iel!=k) && (pt1->qual > declic) ) /*k est enleve par le pop*/
          MMG2_kiudel(queue,iel);
        /*else if ( coe > 0.1 )
          kiuput(queue,iel);  */
        pt1->qual = qual[l];
        pt1->flag = base;
        for(i=0; i<2; i++)  mesh->point[pt1->v[i]].flag = base;
      }

      /* interpol metric */
      ppa->flag = base + 1;
      nm++;
      _MMG5_SAFE_FREE(qual);
      break;
    }
  }
  while ( k );
  if ( mesh->info.imprim < - 4 )
    fprintf(stdout,"     %7d PROPOSED  %7d MOVED %d REJ \n",npp,nm,nrj);

  MMG2_kiufree(queue);
  _MMG5_SAFE_FREE(list);

  return(nm);
}

/* optimise using heap */
int optlen_iso(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppa,ppb;
  pQueue    queue;
  int      *list;
  double    oldc[3],cal,ctg,cx,cy,ux,uy,cpx,cpy,coe,dd,len;
  double    hb,hp,*ca,*cb,*qual;
  int       i,j,k,l,iel,lon,nm;
  int       ipa,ipb,nb,nk,npp,iadr,iter,maxtou;
  int nrj;
//  double tmp1,tmp2,tmp3;
  /* queue on quality */
  queue = MMG2_kiuini(mesh,mesh->nt,declic,base - 1);
  assert(queue);

  maxtou = 10;
  nm     = 0;
  npp    = 0;
  nrj = 0;
  _MMG5_SAFE_MALLOC(list,MMG2D_LMAX,int);

  do {
    k = MMG2_kiupop(queue);

    if ( !k )  break;
    npp++;

    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ipa = pt->v[i];
      ppa = &mesh->point[ipa];
      if ( ppa->tag & M_BDRY || ppa->tag & M_REQUIRED || ppa->tag & M_SD)  continue;

      lon   = MMG2_boulep(mesh,k,i,list);
      _MMG5_SAFE_MALLOC(qual,lon+1,double);

      /* optimal point */
      ca   = &ppa->c[0];
      iadr = ipa*sol->size;
      hp   = sol->m[iadr];
      cx   = 0.0;
      cy   = 0.0;
      nb   = 0;
      cal  = pt->qual;
      for (l=1 ; l<=lon; l++) {
        iel = list[l] / 3;
        nk  = list[l] % 3;
        pt1 = &mesh->tria[iel];
        if ( pt1->qual > cal )  cal = pt1->qual;
        /* tmp1 = 0; */
        /* tmp2 = 0; */
        /* tmp3 = 0; */
        for (j=1; j<3; j++) {
          ipb  = pt1->v[ MMG2_idir[nk + j] ];
          ppb  = &mesh->point[ipb];
          cb   = &ppb->c[0];
          iadr = ipb*sol->size;
          hb   = sol->m[iadr];

          len = MMG2_length(ca,cb,&hp,&hb);
          /* optimal point */
          ux  = ppb->c[0] - ppa->c[0];
          uy  = ppb->c[1] - ppa->c[1];
          cx += ppa->c[0] + ux*(1. - 1./len);//ux * len;
          cy += ppa->c[1] + uy*(1. - 1./len);//uy * len;
          nb++;
        }
        /*check amelioration*/
        /*memcpy(oldc,ppa->c,3*sizeof(double));
          ppa->c[0] = oldc[0] + tmp1/3.;
          ppa->c[1] = oldc[1] + tmp2/3.;
          ppa->c[2] = oldc[2] + tmp3/3.;
          if(MMG_caltet(mesh,sol,iel) > pt1->qual) {
          printf("oups %d -- cal of %d ( %d ) %e > %e\n",nb,iel,ipa,pt1->qual,MMG_caltet(mesh,sol,iel));

          //exit(EXIT_FAILURE);
          }
          else {
          //printf("%d -- cal of %d ( %d ) %e > %e\n",nb,iel,ipa,pt1->qual,MMG_caltet(mesh,sol,iel));
          }
          memcpy(ppa->c,oldc,3*sizeof(double));
        */

      }

      //if ( nb < 3 )  continue;
      dd  = 1.0 / (double)nb;
      cpx = cx*dd - ppa->c[0];
      cpy = cy*dd - ppa->c[1];

      /* adjust position */
      coe  = HQCOEF;
      iter = 1;
      if(cal > 10./ALPHA) {
        ctg  = 0.99 * cal;
      } else {
        ctg  = cal * HCRIT;
      }
      memcpy(oldc,ppa->c,3*sizeof(double));
      do {
        ppa->c[0] =/* (1. - coe) * */oldc[0] + coe * cpx;
        ppa->c[1] =/* (1. - coe) * */oldc[1] + coe * cpy;

        for (l=1; l<=lon; l++) {
          iel = list[l] / 3;
          nk  = list[l] % 3;
          pt1 = &mesh->tria[iel];

          cal = MMG2_caltri_in(mesh,sol,pt1);
          if ( cal > ctg )  {

            break;  }
          qual[l] = cal;
        }
        if ( l > lon )  break;
        coe *= 0.5;
      }
      while ( ++iter <= maxtou );
      if ( iter > maxtou ) {
        memcpy(ppa->c,oldc,3*sizeof(double));
        ppa->flag = base - 2;
        nrj++;
        _MMG5_SAFE_FREE(qual);
        continue;
      }

      /* update tria */
      for (l=1; l<=lon; l++) {
        iel = list[l] / 3;
        nk  = list[l] % 3;
        pt1 = &mesh->tria[iel];
        pt1->qual = qual[l];
        pt1->flag = base;
        for(i=0; i<2; i++)  mesh->point[pt1->v[i]].flag = base;

        if ( pt1->qual < declic )
          MMG2_kiudel(queue,iel);
        /*else if ( coe > 0.1 )
          kiuput(queue,iel);  */
      }

      /* interpol metric */
      ppa->flag = base + 1;
      nm++;
      _MMG5_SAFE_FREE(qual);
      break;
    }
  }
  while ( k );
  if ( mesh->info.imprim < - 4 )
    fprintf(stdout,"     %7d PROPOSED  %7d MOVED %d REJ \n",npp,nm,nrj);


  MMG2_kiufree(queue);
  _MMG5_SAFE_FREE(list);
  return(nm);
}




/* optimise using heap : met le point au barycentre */
int optlen_iso_bar(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppa,ppb;
  pQueue    queue;
  int      *list;
  double    oldc[3],cal,ctg,cx,cy,cpx,cpy,coe,dd;
  double    *qual;
  int       i,j,k,l,iel,lon,nm;
  int       ipa,ipb,nb,nk,npp,iter,maxtou;
  int nrj;
//  double tmp1,tmp2,tmp3;
  /* queue on quality */
  queue = MMG2_kiuini(mesh,mesh->nt,declic,base - 1);
  assert(queue);

  maxtou = 15;
  nm     = 0;
  npp    = 0;
  nrj = 0;
  _MMG5_SAFE_MALLOC(list,MMG2D_LMAX,int);

  do {
    k = MMG2_kiupop(queue);

    if ( !k )  break;
    npp++;

    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ipa = pt->v[i];
      ppa = &mesh->point[ipa];
      if ( ppa->tag & M_BDRY || ppa->tag & M_REQUIRED || ppa->tag & M_SD)  continue;

      lon   = MMG2_boulep(mesh,k,i,list);

      _MMG5_SAFE_MALLOC(qual,lon+1,double);

      /* optimal point */
      cx   = 0.0;
      cy   = 0.0;
      nb   = 0;
      cal  = pt->qual;
      for (l=1 ; l<=lon; l++) {
        iel = list[l] / 3;
        nk  = list[l] % 3;
        pt1 = &mesh->tria[iel];
        if ( pt1->qual > cal )  cal = pt1->qual;
        /* tmp1 = 0; */
        /* tmp2 = 0; */
        /* tmp3 = 0; */
        for (j=1; j<3; j++) {
          ipb  = pt1->v[ MMG2_idir[nk + j] ];
          ppb  = &mesh->point[ipb];

          /* optimal point */
          cx += ppb->c[0] ;
          cy += ppa->c[1] ;
          nb++;
        }
        /*check amelioration*/
        /*memcpy(oldc,ppa->c,3*sizeof(double));
          ppa->c[0] = oldc[0] + tmp1/3.;
          ppa->c[1] = oldc[1] + tmp2/3.;
          ppa->c[2] = oldc[2] + tmp3/3.;
          if(MMG_caltet(mesh,sol,iel) > pt1->qual) {
          printf("oups %d -- cal of %d ( %d ) %e > %e\n",nb,iel,ipa,pt1->qual,MMG_caltet(mesh,sol,iel));

          //exit(EXIT_FAILURE);
          }
          else {
          //printf("%d -- cal of %d ( %d ) %e > %e\n",nb,iel,ipa,pt1->qual,MMG_caltet(mesh,sol,iel));
          }
          memcpy(ppa->c,oldc,3*sizeof(double));
        */

      }

      //if ( nb < 3 )  continue;
      dd  = 1.0 / (double)nb;
      cpx = cx*dd ;//- ppa->c[0];
      cpy = cy*dd ;//- ppa->c[1];

      /* adjust position */
      coe  = HQCOEF;
      iter = 1;
      if(cal > 10./ALPHA) {
        ctg  = 0.99 * cal;
      } else {
        ctg  = cal * HCRIT;
      }
      memcpy(oldc,ppa->c,3*sizeof(double));
      do {
        ppa->c[0] = (1. - coe) * oldc[0] + coe * cpx;
        ppa->c[1] = (1. - coe) * oldc[1] + coe * cpy;

        for (l=1; l<=lon; l++) {
          iel = list[l] / 3;
          nk  = list[l] % 3;
          pt1 = &mesh->tria[iel];

          cal = MMG2_caltri_in(mesh,sol,pt1);
          if ( cal > ctg )  {
            break;  }
          qual[l] = cal;
        }
        if ( l > lon )  break;
        coe *= 0.5;
      }
      while ( ++iter <= maxtou );
      if ( iter > maxtou ) {
        memcpy(ppa->c,oldc,3*sizeof(double));
        ppa->flag = base - 2;
        nrj++;
        _MMG5_SAFE_FREE(qual);
        continue;
      }

      /* update tria */
      for (l=1; l<=lon; l++) {
        iel = list[l] / 3;
        nk  = list[l] % 3;
        pt1 = &mesh->tria[iel];
        pt1->qual = qual[l];
        pt1->flag = base;
        for(i=0; i<2; i++)  mesh->point[pt1->v[i]].flag = base;

        if ( pt1->qual < declic )
          MMG2_kiudel(queue,iel);
        /*else if ( coe > 0.1 )
          kiuput(queue,iel);  */
      }

      /* interpol metric */
      ppa->flag = base + 1;
      nm++;
      _MMG5_SAFE_FREE(qual);
      break;
    }
  }
  while ( k );
  if ( mesh->info.imprim < - 4 )
    fprintf(stdout,"     %7d PROPOSED  %7d MOVED %d REJ \n",npp,nm,nrj);

  _MMG5_SAFE_FREE(list);
  MMG2_kiufree(queue);
  return(nm);
}
