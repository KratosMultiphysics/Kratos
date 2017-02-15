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
 * \file mmg2d/split_2d.c
 * \brief Functions for splitting.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg2d.h"

#define QSEUIL 1e4

/*insert ip on edge between k1 and adj1/3 */
int MMG2_split(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int k1,int adj1,double lshort) {
  MMG5_pTria     pt1,pt2,pt3,pt4,ptmp;
  MMG5_pEdge     ped,ped1;
  int       k2,adj2,jel,kel,voy1,voy2,iar1,iar2,iara1,iara2;
  int       *adja,*adja1,*adja2,tmp1,tmp2,piar1,piar2,pvoy1,piara1,piara2,pvoy2;
  int       iadr,tmp,voy,num,newed,num1,num2;
  double    air,cal1,cal2,cal3,cal4,coe,*ca,*cb,*ma,*mb,len;

  coe = QSEUIL/ALPHA;

  voy2  = adj1%3;
  k2    = adj1/3;
  adja2 = &mesh->adja[3*(k2-1) + 1];
  adj2  = adja2[voy2];
  assert(adj2/3==k1);
  voy1  = adj2%3;

  iar1  = MMG2_iare[voy1][0];
  iar2  = MMG2_iare[voy1][1];
  iara1 = MMG2_iare[voy2][0];
  iara2 = MMG2_iare[voy2][1];

  pt1  = &mesh->tria[k1];
  pt2  = &mesh->tria[k2];

  assert(pt2->v[iara1]==pt1->v[iar2]);
  assert(pt2->v[iara2]==pt1->v[iar1]);

  piar1  = pt1->v[iar1];
  piar2  = pt1->v[iar2];
  pvoy1  = pt1->v[voy1];
  piara1 = pt2->v[iara1];
  piara2 = pt2->v[iara2];
  pvoy2  = pt2->v[voy2];

  /*test split ok*/
  /*test area > 0*/
  air = MMG2_quickarea(mesh->point[piar2].c,mesh->point[pvoy1].c,mesh->point[ip].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[piara1].c,mesh->point[ip].c,mesh->point[pvoy2].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[ip].c,mesh->point[pvoy1].c,mesh->point[piar1].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[ip].c,mesh->point[piara2].c,mesh->point[pvoy2].c);
  if(air < EPSA) return(0);

  /*test qual*/
  ptmp = &mesh->tria[0];
  ptmp->v[0] = piar2;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = ip;
  cal1       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal1 > coe) return(0);

  ptmp->v[0] = piara1;
  ptmp->v[1] = ip;
  ptmp->v[2] = pvoy2;
  cal2       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal2 > coe) return(0);

  ptmp->v[0] = ip;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = piar1;
  cal3       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal3 > coe) return(0);

  ptmp->v[0] = ip;
  ptmp->v[1] = piara2;
  ptmp->v[2] = pvoy2;
  cal4       = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal4 > coe) return(0);

  /*test length : pvoy1-ip and pvoy2-ip*/
  ca   = &mesh->point[ip].c[0];
  iadr = ip*sol->size;
  ma   = &sol->m[iadr];

  cb  = &mesh->point[pvoy1].c[0];
  iadr = pvoy1*sol->size;
  mb   = &sol->m[iadr];

  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d %d : %e\n",pvoy1,ip,len);
  if(len < lshort) return(0);

  cb  = &mesh->point[pvoy2].c[0];
  iadr = pvoy2*sol->size;
  mb   = &sol->m[iadr];

  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d %d : %e\n",pvoy2,ip,len);
  if(len < lshort) return(0);

  /*check*/
  cb  = &mesh->point[piara2].c[0];
  iadr = piara2*sol->size;
  mb   = &sol->m[iadr];

  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d (%d)  %d : %e\n",piar1,piara2,ip,len);
  if(len < lshort) return(0);

  cb  = &mesh->point[piar2].c[0];
  iadr = piar2*sol->size;
  mb   = &sol->m[iadr];

  len = MMG2_length(ca,cb,ma,mb);
  //printf("edg %d %d : %e\n",piar2,ip,len);
  if(len < lshort) return(0);

  pt1->v[0] = piar2;
  pt1->v[1] = pvoy1;
  pt1->v[2] = ip;
  pt1->qual = cal1;

  pt2->v[0] = piara1;
  pt2->v[1] = ip;
  pt2->v[2] = pvoy2;
  pt2->qual = cal2;


  jel  = _MMG2D_newElt(mesh);
  if ( !jel ) {
    _MMG5_TRIA_REALLOC(mesh,jel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    pt1  = &mesh->tria[k1];
    pt2  = &mesh->tria[k2];
    ptmp = &mesh->tria[0];
    adja2 =  &mesh->adja[3*(k2-1) + 1];
  }
  kel  = _MMG2D_newElt(mesh);
  if ( !kel ) {
    _MMG5_TRIA_REALLOC(mesh,kel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    pt1  = &mesh->tria[k1];
    pt2  = &mesh->tria[k2];
    ptmp = &mesh->tria[0];
    adja2 =  &mesh->adja[3*(k2-1) + 1];
  }
  pt3  = &mesh->tria[jel];
  pt3->v[0] = ip;
  pt3->v[1] = pvoy1;
  pt3->v[2] = piar1;
  pt3->ref  = pt1->ref;
  pt3->qual = cal3;

  pt4  = &mesh->tria[kel];
  pt4->v[0] = ip;
  pt4->v[1] = piara2;
  pt4->v[2] = pvoy2;
  pt4->ref  = pt2->ref;
  pt4->qual = cal4;

  /*adj*/
  adja1 = &mesh->adja[3*(k1-1) + 1];
  /*printf("adj1 : %d %d %d\n",adja1[0]/3,adja1[1]/3,adja1[2]/3);
    printf("adj2 : %d %d %d\n",adja2[0]/3,adja2[1]/3,adja2[2]/3);
  */
  adja = &mesh->adja[3*(jel-1) + 1];
  adja[0] = adja1[iar2];
  pt3->edg[0] = pt1->edg[iar2];
  tmp = 3*(adja1[iar2]/3 - 1) + 1;
  voy = adja1[iar2]%3;
  if(adja1[iar2]) (&mesh->adja[tmp])[voy] = 3*jel + 0;
  adja[1] = 3*kel + 2;
  num = 0;
  if(pt1->edg[voy1]) {
    /*edge creation*/
    /*split edge piar1 piar2 */
// #warning same tangent
    num = pt1->edg[voy1];
    assert(num);
    ped = &mesh->edge[num];
    newed = _MMG5_newEdge(mesh);
    if ( !newed ) {
      _MMG5_EDGE_REALLOC(mesh,newed,mesh->gap,
                         printf("  ## Error: unable to allocate a new edge.\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      ped = &mesh->edge[num];
    }

    ped1 = &mesh->edge[newed];
    memcpy(ped1,ped,sizeof(MMG5_Edge));
    ped1->a = ip;

    ped->b = ip;

    if(ped->a==piar1) {
      pt3->edg[1] = num;
      pt4->edg[2] = num;
      // printf("on met num pt3  %d %d\n",pt3->v[MMG2_iare[1][0]],pt3->v[MMG2_iare[1][1]]);
      //printf("on met num pt4 %d %d\n",pt4->v[MMG2_iare[2][0]],pt4->v[MMG2_iare[2][1]]);

    } else {
      pt3->edg[1] = newed;
      pt4->edg[2] = newed;
    }
  }
  adja[2] = 3*k1  + 0;
  pt3->edg[2] = 0;

  if(kel) {
    adja = &mesh->adja[3*(kel-1) + 1];
    adja[0] = adja2[iara1];
  }
  pt4->edg[0] = pt2->edg[iara1];


  if(adja2[iara1]) (&mesh->adja[3*(adja2[iara1]/3-1)+1])[adja2[iara1]%3] = 3*kel + 0;
  adja[1] = 3*k2  + 0;
  pt4->edg[1] = 0;
  adja[2] = 3*jel + 1;
  if(pt1->edg[voy1]) assert(pt4->edg[2]);

  tmp1 = adja1[iar1];
  num1 = pt1->edg[iar1];
  tmp2 = adja2[iara2];
  num2 = pt2->edg[iara2];
  adja1[0] = 3*jel + 2;
  pt1->edg[0] = 0;
  adja1[1] = 3*k2  + 2;

  adja1[2] = tmp1;
  pt1->edg[2] = num1;
  if(tmp1)
    (&mesh->adja[3*(tmp1/3-1)+1])[tmp1%3] = 3*k1 + 2;


  adja2[0] = 3*kel + 1;
  pt2->edg[0] = 0;
  adja2[1] = tmp2;
  pt2->edg[1] = num2;
  if(tmp2)
    (&mesh->adja[3*(tmp2/3-1)+1])[tmp2%3] = 3*k2 + 1;
  adja2[2] = 3*k1 + 1;
  if(num) {
    if(ped->a==piar1) {
      pt1->edg[1] = newed;
      pt2->edg[2] = newed;
    } else {
      pt1->edg[1] = num;
      pt2->edg[2] = num;
      //printf("on met num k1 %d %d\n",pt1->v[MMG2_iare[1][0]],pt1->v[MMG2_iare[1][1]]);
      //printf("on met num k2 %d %d\n",pt2->v[MMG2_iare[2][0]],pt2->v[MMG2_iare[2][1]]);

    }
  } else {
    pt1->edg[1] = 0;
    pt2->edg[2] = 0;
  }

  /* if(MMG2D_callbackinsert) */
  /*   MMG2D_callbackinsert((int) ip,(int) k1,(int) k2,(int)jel,(int) kel); */

  return(1);
}

/*insert ip on edge in k1 */
int MMG2_splitbdry(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int k1,int voy1,double *tang) {
  MMG5_pTria     pt1,pt3,ptmp;
  MMG5_pEdge     ped,ped1;
  MMG5_pPoint    ppt;
  int       jel,iar1,iar2,i,num,newed,num1,num2;
  int       *adja,*adja1,tmp1,piar1,piar2,pvoy1,ref1,ref2;
  double    air,cal1,cal2,coe;

  coe = QSEUIL/ALPHA;

  iar1  = MMG2_iare[voy1][0];
  iar2  = MMG2_iare[voy1][1];

  pt1  = &mesh->tria[k1];
  //printf("tr1 %d %d %d -- voy %d : %d %d \n",pt1->v[0],pt1->v[1],pt1->v[2],voy1,iar1,iar2);

  piar1  = pt1->v[iar1];
  piar2  = pt1->v[iar2];
  pvoy1  = pt1->v[voy1];

  /*test split ok*/
  air = MMG2_quickarea(mesh->point[piar2].c,mesh->point[pvoy1].c,mesh->point[ip].c);
  if(air < EPSA) return(0);
  air = MMG2_quickarea(mesh->point[ip].c,mesh->point[pvoy1].c,mesh->point[piar1].c);
  if(air < EPSA) return(0);

  ptmp = &mesh->tria[0];
  ptmp->v[0] = piar2;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = ip;
  cal1 = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal1 > coe) return(0);

  ptmp->v[0] = ip;
  ptmp->v[1] = pvoy1;
  ptmp->v[2] = piar1;
  cal2 = MMG2_caltri_in(mesh,sol,ptmp);
  if(cal2 > coe) return(0);


  pt1->v[0] = piar2;
  pt1->v[1] = pvoy1;
  pt1->v[2] = ip;
  pt1->qual = cal1;

  jel  = _MMG2D_newElt(mesh);
  if ( !jel ) {
    _MMG5_TRIA_REALLOC(mesh,jel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    pt1  = &mesh->tria[k1];
    ptmp = &mesh->tria[0];

  }
  pt3  = &mesh->tria[jel];
  pt3->v[0] = ip;
  pt3->v[1] = pvoy1;
  pt3->v[2] = piar1;
  pt3->ref  = pt1->ref;
  pt3->qual = cal2;

  /*adj*/
  adja1 = &mesh->adja[3*(k1-1) + 1];

  adja = &mesh->adja[3*(jel-1) + 1];
  adja[0] = adja1[iar2];
  if(adja1[iar2])
    (&mesh->adja[3*(adja1[iar2]/3-1)+1])[adja1[iar2]%3] = 3*jel + 0;
  adja[1] = 0;
  adja[2] = 3*k1  + 0;

  tmp1 = adja1[iar1];
  adja1[0] = 3*jel + 2;
  adja1[1] = 0;
  adja1[2] = tmp1;
  if(tmp1)
    (&mesh->adja[3*(tmp1/3-1)+1])[tmp1%3] = 3*k1 + 2;

  /*si dep alors on met la moy des dep dans ip*/
  if(mesh->info.lag >=0) {
    printf("  ## Error: option not implemented: merge option 9\n");
    exit(EXIT_FAILURE);
    /* mesh->disp.mv[2*(ip-1) + 1 + 0] = 0.5*(mesh->disp.mv[2*(piar1-1) + 1 + 0] + mesh->disp.mv[2*(piar2-1) + 1 + 0]); */
    /* mesh->disp.mv[2*(ip-1) + 1 + 1] = 0.5*(mesh->disp.mv[2*(piar1-1) + 1 + 1] + mesh->disp.mv[2*(piar2-1) + 1 + 1]);     */
    /* d1 = mesh->disp.mv[2*(ip-1) + 1 + 0]*mesh->disp.mv[2*(ip-1) + 1 + 0] */
    /*   + mesh->disp.mv[2*(ip-1) + 1 + 1]*mesh->disp.mv[2*(ip-1) + 1 + 1]; */
    /* if ( d1 > 1e-24 )  mesh->point[ip].tag  |= M_MOVE; */
  }

  /*propagation des ref de peau*/
  ref1 = mesh->point[piar1].ref;
  ref2 = mesh->point[piar2].ref;
  if( ref1 || ref2 ) mesh->point[ip].ref = ref1 > ref2 ? ref1 : ref2;

  /*split edge piar1 piar2 if exist */
  num = pt1->edg[voy1];
  assert(num);
  ped = &mesh->edge[num];

  newed = _MMG5_newEdge(mesh);
  if ( !newed ) {
    _MMG5_EDGE_REALLOC(mesh,newed,mesh->gap,
                       printf("  ## Error: unable to allocate a new edge.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
    ped = &mesh->edge[num];
  }

  ped1 = &mesh->edge[newed];
  memcpy(ped1,ped,sizeof(MMG5_Edge));
  ped1->a = ip;
  ppt = &mesh->point[ip];
  for(i=0 ; i<2 ; i++)
    ppt->n[i] = tang[i];

  ped->b = ip;

  num1 = pt1->edg[iar1];
  num2 = pt1->edg[iar2];
  pt1->edg[2] = num1;
  pt3->edg[0] = num2;

  if(ped->a==piar1) {
    pt3->edg[1] = num;
    pt1->edg[1] = newed;
  } else {
    pt3->edg[1] = newed;
    pt1->edg[1] = num;
  }

  pt1->edg[0]  = 0;
  //end add edge

  return(1);
}
