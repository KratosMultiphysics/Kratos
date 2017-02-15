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
#define LLONG1 1.9

/*collapse edge ppb-->ppa*/
int MMG2_colpoi(MMG5_pMesh mesh, MMG5_pSol sol,int iel,int iar,int ia,int ib,double coe) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppa,ppb,pp1,pp2,ppa1,ppb1;
  MMG5_pEdge     ped;
  int       pib,pia,jel,iadr,a1,v1,*adja,iaa,voy,a,a2,v2,adj,adj1,adjj1;
  int     *list,lon,i,kel,num,ed,i1,i2;
  double    declic,*cal,air,coor[2],solu[3],*c1,*c2,*m1,*m2,len;
  double    cbound,capx,capy,cbpx,cbpy,alpha;
  int       nbdry,ref,ip,iadri,*adjai,ibdry[2],it1;

  pt  = &mesh->tria[iel];
  pib = pt->v[ib];
  pia = pt->v[ia];
  ppa = &mesh->point[pia];
  ppb = &mesh->point[pib];

  if(ppb->tag & M_BDRY) return(0);

  iadr = 3*(iel-1) + 1;
  adja = &mesh->adja[iadr];
  jel  = adja[iar]/3;

  list  = (int*)malloc(MMG2D_LMAX*sizeof(int));
  assert(list);

  lon = MMG2_boulep(mesh,iel,ib,list);
  if(!lon) {
    free(list);
    return(0);
  }

  //if vertex between two SD, check geometry criterion
  if(ppb->tag & M_SD) {
    if(!(ppa->tag & M_SD)) {
      free(list);
      return(0);
    }
    /*check pia-pib edge exist otherwise the collapse is forbidden (ex: naca trailing edge)*/
    if(!pt->edg[iar]) {
      free(list);
      return(0);
    }

    /*check geom*/
    cbound = 178.*M_PI/180.;
    nbdry = 0;
    ref = mesh->tria[list[lon]/3].ref;
    for(i=1 ; i<=lon ; i++) {
      kel = list[i]/3;
      ip  = list[i]%3;
      pt1 = &mesh->tria[kel];
      iadri = 3*(kel-1) + 1;
      adjai =  &mesh->adja[iadri];
      if(pt1->ref!=ref) {
        it1 = (i==1)?lon:i-1;
        if(adjai[MMG2_iopp[ip][0]]/3 == list[it1]/3 ) {
          //assert(nbdry<2);
          if(MMG2_iare[MMG2_iopp[ip][0]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][0]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][0]];
          }
        } else {
          assert(adjai[MMG2_iopp[ip][1]]/3 == list[it1]/3);
          if(MMG2_iare[MMG2_iopp[ip][1]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][1]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][0]];
          }
        }
      }
      ref = pt1->ref;
    }
    if(nbdry!=2) {
      _MMG5_SAFE_FREE(list);
      return(0);
    }
    //assert(nbdry==2); //sinon non manifold
    //calcul de l'angle forme par les 3 points
    ppa1  = &mesh->point[ibdry[0]];
    ppb1  = &mesh->point[ibdry[1]];
    capx = ppb->c[0] - ppa1->c[0];
    capy = ppb->c[1] - ppa1->c[1];
    cbpx = ppb->c[0] - ppb1->c[0];
    cbpy = ppb->c[1] - ppb1->c[1];
    alpha = capx*cbpx + capy*cbpy;
    alpha /= sqrt(capx*capx+capy*capy)*sqrt(cbpx*cbpx+cbpy*cbpy);
    alpha = acos(alpha);

    if(alpha < cbound ) {
      free(list);
      return(0);
    }
  }

  cal  = (double*)malloc((lon+1)*sizeof(double));
  assert(cal);

  /*simu colps ppb-->ppa*/
  memcpy(coor, ppb->c,2*sizeof(double));
  memcpy(ppb->c,ppa->c,2*sizeof(double));
  //memcpy(solu,&sol->m[3*pib],sol->size*sizeof(double));
  for(i=0 ; i<sol->size ; i++) {
    solu[i] = sol->m[sol->size*pib + i];
    sol->m[sol->size*pib + i] = sol->m[sol->size*pia + i];
  }

  //memcpy(&sol->m[3*pib],&sol->m[3*pia],sol->size*sizeof(double));


  /*check config*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    air  = MMG2_quickarea(mesh->point[pt1->v[0]].c,mesh->point[pt1->v[1]].c,
                          mesh->point[pt1->v[2]].c);
    if(air < EPSA) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));
      free(cal);
      free(list);
      return(0);
    }
    declic = coe*pt1->qual;
    cal[i] = MMG2_caltri_in(mesh,sol,pt1);
    if (cal[i] > declic) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));
      free(cal);
      free(list);
      return(0);
    }
  }

  /*check lengths*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    /*check first edge containing ib*/
    ed   = (voy+2)%3;
    i1   = pt1->v[MMG2_iare[ed][0]];
    i2   = pt1->v[MMG2_iare[ed][1]];
    pp1  = &mesh->point[i1];
    pp2  = &mesh->point[i2];
    c1   = &pp1->c[0];
    c2   = &pp2->c[0];
    iadr = i1*sol->size;
    m1   = &sol->m[iadr];
    iadr = i2*sol->size;
    m2   = &sol->m[iadr];

    len = MMG2_length(c1,c2,m1,m2);
    if (len > LLONG1) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));
      free(cal);
      free(list);
      return(0);
    }
  }

  /*update tria*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];

    assert(pt1->v[voy]==pib);
    pt1->v[voy] = pia;

    /*edge*/
    num  = pt1->edg[MMG2_iare[voy][0]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;
    }
    num = pt1->edg[ MMG2_iare[voy][1]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;

    }
    pt1->qual = cal[i];
  }

  /*adj of iel*/
  adj  = adja[ib];
  a1   = adja[ia]/3;
  v1   = adja[ia]%3;

  iaa  = MMG2_iare[adja[iar]%3][1];

  adjj1 = MMG2_iare[adja[iar]%3][0];

  /*adj of jel*/
  iadr = 3*(jel-1) + 1;
  adja = &mesh->adja[iadr];
  adj1 = adja[adjj1];
  a2   = adja[iaa]/3;
  v2   = adja[iaa]%3;


  /*adja*/
  iadr = 3*(a1-1) + 1;
  a    = adj/3;
  voy  = adj%3;
  adja = &mesh->adja[iadr];
  adja[v1] = 3*a + voy;

  if(a) {
    iadr = 3*(a-1) + 1;
    adja = &mesh->adja[iadr];
    adja[voy] = 3*a1 + v1;
  }
  /*update edge*/
  num = pt->edg[ib];
  if(num) {
    mesh->tria[a1].edg[v1] = num;
  }

  iadr = 3*(a2-1) + 1;
  a    = adj1/3;
  voy  = adj1%3;
  if(a2) {
    adja = &mesh->adja[iadr];
    adja[v2] = 3*a + voy;
  }
  if(a) {
    iadr = 3*(a-1) + 1;
    adja = &mesh->adja[iadr];
    adja[voy] = 3*a2 + v2;
  }
  /*update edge*/
  pt1  = &mesh->tria[jel];
  num = pt1->edg[adjj1];
  if(a2 && num) {
    if(!((mesh->edge[num].a==mesh->tria[a2].v[MMG2_iare[v2][0]] || mesh->edge[num].a==mesh->tria[a2].v[MMG2_iare[v2][1]])
         && (mesh->edge[num].b==mesh->tria[a2].v[MMG2_iare[v2][0]] ||mesh->edge[num].b==mesh->tria[a2].v[MMG2_iare[v2][1]]))) {
      /* printf("on a un soucis 1\n"); */
      /* printf("pnum %d %d dans %d %d %d\n",mesh->edge[num].a,mesh->edge[num].b,pt1->v[0],pt1->v[1],pt1->v[2]); */
      /* printf("edgea %d %d\n",mesh->tria[a2].v[MMG2_iare[v2][0]],mesh->tria[a2].v[MMG2_iare[v2][1]]); */
      if(mesh->info.imprim > 6)
        fprintf(stdout," ## Warning: bad configuration for collapse\n");
    }
    mesh->tria[a2].edg[v2] = num;
  }
  num = pt1->edg[iaa];
  if(a && num) {
    /* printf("tr %d : %d %d %d et num %d\n",a,mesh->tria[a].edg[0],mesh->tria[a].edg[1],mesh->tria[a].edg[2],num); */
    /* printf("pnum %d %d dans %d %d %d\n",mesh->edge[num].a,mesh->edge[num].b,pt1->v[0],pt1->v[1],pt1->v[2]); */
    /* printf("edgea %d %d\n",mesh->tria[a].v[MMG2_iare[voy][0]],mesh->tria[a].v[MMG2_iare[voy][1]]); */
    if(!((mesh->edge[num].a==mesh->tria[a].v[MMG2_iare[voy][0]] || mesh->edge[num].a==mesh->tria[a].v[MMG2_iare[voy][1]])
         && (mesh->edge[num].b==mesh->tria[a].v[MMG2_iare[voy][0]] || mesh->edge[num].b==mesh->tria[a].v[MMG2_iare[voy][1]]))) {
      /* printf("on a un soucis\n"); */
      /* printf("pnum %d %d dans %d %d %d\n",mesh->edge[num].a,mesh->edge[num].b,pt1->v[0],pt1->v[1],pt1->v[2]); */
      /* printf("edgea %d %d\n",mesh->tria[a].v[MMG2_iare[voy][0]],mesh->tria[a].v[MMG2_iare[voy][1]]); */
      if(mesh->info.imprim > 6) fprintf(stdout," ## Warning: bad configuration for collapse\n");
    }
    mesh->tria[a].edg[voy] = num;
  }

  num = pt->edg[iar];
  if(num) {
    _MMG5_delEdge(mesh,num);
  }


  _MMG2D_delElt(mesh,iel);
  _MMG2D_delElt(mesh,jel);
  memcpy(ppb->c,coor,2*sizeof(double));
  memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));

  free(list);
  free(cal);
  return(1);
}

/**
 * \brief return 1 if the edge does not verify hausd criterion
 */
int MMG2_chkedg(MMG5_pMesh mesh, MMG5_pPoint ppa,MMG5_pPoint ppb) {
  double t0[2],t1[2];
  double ll,l,ux,uy,cosn,ps;
  int    i;

  ux = ppa->c[0] - ppb->c[0];
  uy = ppa->c[1] - ppb->c[1];
  ll = ux*ux + uy*uy;
  l = sqrt(ll);

  if ( l > mesh->info.hmax ) return(1);
  else if ( l < _MMG5_EPSD ) return(0);

  //with tangent, compute control point
  for(i=0 ; i<2 ; i++) {
    t0[i] = ppa->n[i];
    t1[i] = ppb->n[i];
  }
  if(ppa->tag & M_CORNER) {
    for(i=0 ; i<2 ; i++)
      t0[i] = (ppa->c[i] - ppb->c[i])/l;
  }
  if(ppb->tag & M_CORNER) {
    for(i=0 ; i<2 ; i++)
      t1[i] = (ppa->c[i] - ppb->c[i])/l;
  }
  /*check if t0 has the same sens of vect(P0P1)*/
  if(t0[0]/(ppb->c[0]-ppa->c[0]) < 0 || t0[1]/(ppb->c[1]-ppa->c[1])<0) {
    //printf("t0/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t0[i] *= -1;
    }
  }
  /*check if t1 has the opposite sens of vect(P0P1)*/
  if(t1[0]/(ppb->c[0]-ppa->c[0]) > 0 || t1[1]/(ppa->c[1]-ppb->c[1])>0) {
    //printf("t1/pOp1 %e %e\n",t0[0]/(p1->c[0]-p0->c[0]),t0[1]/(p1->c[1]-p0->c[1]));
    for(i=0 ; i<2 ; i++) {
      t1[i] *= -1;
    }
  }
  //compute the distance between mid point and curve (with angle between edge and tang)
  ps = t0[0]*ux + t0[1]*uy;

  ps *= ps;
  cosn = ps/ll ;
  cosn *= fabs(1.0-cosn);
  cosn *= ll;

  if(cosn > 9.*mesh->info.hausd*mesh->info.hausd) return(1);
  //idem for ppb
  ps = -(t1[0]*ux + t1[1]*uy);
  ps *= ps;
  cosn = ps/ll;
  cosn *= fabs(1.0-cosn);
  cosn *=  ll;
  if(cosn > 9.*mesh->info.hausd*mesh->info.hausd) return(1);

  return(0);
}

/**
 * \return -1 if fail, 0 if we can't collapse, 1 otherwise.
 *
 * collapse edge ppb-->ppa and ppbppa is boundary edge
 *
 */
int MMG2_colpoibdry(MMG5_pMesh mesh, MMG5_pSol sol,int iel,int iar,int ia,int ib,double coe) {
  MMG5_pTria     pt,pt1;
  MMG5_pEdge     ped;
  MMG5_pPoint    ppa,ppb,pp1,pp2;//ppa1,ppb1;
  int       pib,pia,jel,iadr,a1,v1,*adja,voy,a,adj;
  int     *list,lon,i,kel,num,i1,i2,ed;
  double    declic,*cal,air,coor[2],solu[3],*c1,*c2,*m1,*m2,len;
  //double    capx,capy,cbpx,cbpy,alpha,cbound;
//  int       iadri,*adjai,nbdry,ip,ibdry;

  pt  = &mesh->tria[iel];
  pib = pt->v[ib];
  pia = pt->v[ia];
  ppa = &mesh->point[pia];
  ppb = &mesh->point[pib];

  if (ppb->tag & M_CORNER) return(0);
  iadr = 3*(iel-1) + 1;
  adja = &mesh->adja[iadr];
  jel  = adja[iar]/3;
  assert(!jel);
  _MMG5_SAFE_MALLOC(list,MMG2D_LMAX,int);

  lon = MMG2_boulep(mesh,iel,ib,list);
  if(!lon) {
    _MMG5_SAFE_FREE(list);
    return(0);
  }

  /*check geom*/
  /* nbdry = 0; */
  /* for(i=1 ; i<=lon ; i++) { */
  /*   kel = list[i]/3; */
  /*   ip  = list[i]%3; */
  /*   pt1 = &mesh->tria[kel]; */
  /*   iadri = 3*(kel-1) + 1; */
  /*   adjai = &mesh->adja[iadri]; */
  /*   if(!adjai[MMG2_iopp[ip][0]]) { */
  /*     assert(nbdry<2); */
  /*     if(MMG2_iare[MMG2_iopp[ip][0]][0]==ip) */
  /*       ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][1]]; */
  /*     else { */
  /*       assert(MMG2_iare[MMG2_iopp[ip][0]][1]==ip) ; */
  /*       ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][0]]; */
  /*     } */
  /*   } */
  /*   if(!adjai[MMG2_iopp[ip][1]]) { */
  /*     assert(nbdry<2); //sinon non manifold */
  /*     if(MMG2_iare[MMG2_iopp[ip][1]][0]==ip) */
  /*       ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][1]]; */
  /*     else { */
  /*       assert(MMG2_iare[MMG2_iopp[ip][1]][1]==ip) ; */
  /*       ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][0]]; */
  /*     } */
  /*   } */
  /* } */
  /* assert(nbdry==2); //sinon non manifold */
/*first check that the two edges verify the hausd criterion*/
  /* ppa1  = &mesh->point[ibdry[0]]; */
  /* ppb1  = &mesh->point[ibdry[1]]; */

  /* if(MMG2_chkedg(mesh,ppb,ppa1))   { */
  /*   _MMG5_SAFE_FREE(list); */
  /*   return(0); */
  /* } */
  /* if(MMG2_chkedg(mesh,ppb,ppb1))  { */
  /*   _MMG5_SAFE_FREE(list); */
  /*   return(0); */
  /* } */

  /* /\*second check that the new edge verify the hausd criteron*\/ */
  /* if(MMG2_chkedg(mesh,ppb1,ppa1))  { */
  /*   _MMG5_SAFE_FREE(list); */
  /*   return(0); */
  /* } */

/* //comment from here */
/*   //calcul de l'angle forme par les 3 points  */
/*   capx = ppb->c[0] - ppa1->c[0];  */
/*   capy = ppb->c[1] - ppa1->c[1];  */
/*   cbpx = ppb->c[0] - ppb1->c[0];  */
/*   cbpy = ppb->c[1] - ppb1->c[1];  */
/*   alpha = capx*cbpx + capy*cbpy; */
/*   alpha /= sqrt(capx*capx+capy*capy)*sqrt(cbpx*cbpx+cbpy*cbpy);   */
/*   alpha = acos(alpha); */
/*   //printf("point %d : %e (= %e)-- %e %e\n",pt->v[j],alpha,alpha*180./M_PI,capx,capy); */

/*   if(alpha < cbound ) { */
/*     free(list); */
/*     return(0); */
/*     //to here */
/*   } else */
  if(lon > 100) {
    _MMG5_SAFE_FREE(list);
    return(0);
  }
  _MMG5_SAFE_CALLOC(cal,lon+1,double);

  /*simu colps ppb-->ppa*/
  memcpy(coor, ppb->c,2*sizeof(double));
  memcpy(ppb->c,ppa->c,2*sizeof(double));
  memcpy(solu,&sol->m[sol->size*pib],sol->size*sizeof(double));
  memcpy(&sol->m[sol->size*pib],&sol->m[sol->size*pia],sol->size*sizeof(double));

  /*check config*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    air  = MMG2_quickarea(mesh->point[pt1->v[0]].c,mesh->point[pt1->v[1]].c,
                          mesh->point[pt1->v[2]].c);
    if(air < EPSA) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));
      _MMG5_SAFE_FREE(cal);
      _MMG5_SAFE_FREE(list);
      return(0);
    }
    declic = coe*pt1->qual;
    cal[i] = MMG2_caltri_in(mesh,sol,pt1);
    if (cal[i] > declic) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));
      _MMG5_SAFE_FREE(cal);
      _MMG5_SAFE_FREE(list);
      return(0);
    }
  }

  /*check lengths*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    if(kel==jel) continue;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    /*check second edge containing ib*/
    ed   = (voy+2)%3;
    i1   = pt1->v[MMG2_iare[ed][0]];
    i2   = pt1->v[MMG2_iare[ed][1]];
    pp1  = &mesh->point[i1];
    pp2  = &mesh->point[i2];
    c1   = &pp1->c[0];
    c2   = &pp2->c[0];
    iadr = i1*sol->size;
    m1   = &sol->m[iadr];
    iadr = i2*sol->size;
    m2   = &sol->m[iadr];

    len = MMG2_length(c1,c2,m1,m2);
    if (len > LLONG1) {
      memcpy(ppb->c,coor,2*sizeof(double));
      memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));
      _MMG5_SAFE_FREE(cal);
      _MMG5_SAFE_FREE(list);
      return(0);
    }
  }
  /*update tria*/
  for(i=2 ; i<=lon ; i++) {
    kel = list[i]/3;
    voy = list[i]%3;
    pt1 = &mesh->tria[kel];
    assert(pt1->v[voy]==pib);
    pt1->v[voy] = pia;

    pt1->qual = cal[i];

    /*edge*/
    num = pt1->edg[ MMG2_iare[voy][0]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;
    }
    num = pt1->edg[ MMG2_iare[voy][1]];
    if(num) {
      ped = &mesh->edge[num];
      if(ped->a==pib) ped->a = pia;
      if(ped->b==pib) ped->b = pia;
    }
  }

  adj  = adja[ib];
  a1   = adja[ia]/3;
  v1   = adja[ia]%3;
  pt1  = &mesh->tria[a1];

  /*adja*/
  a    = adj/3;
  voy  = adj%3;
  if(a1) {
    iadr = 3*(a1-1) + 1;
    adja = &mesh->adja[iadr];
    adja[v1] = 3*a + voy;
  }
  if(a) {
    iadr = 3*(a-1) + 1;
    adja = &mesh->adja[iadr];
    adja[voy] = 3*a1 + v1;
  }

  /*del edge pib pia if exist*/
  num = pt->edg[iar];
  if(!num) {
    /* printf("la edge %d %d iar %d tr %d\n",pia,pib,iar,iel); */
    /* printf("pt %d %d %d\n",pt->v[0],pt->v[1],pt->v[2]); */
    /* printf("pt->ned %d %d %d\n",pt->edg[0],pt->edg[1],pt->edg[2]); */
    fprintf(stdout," ## Error: problem with an edge."
            " Check your data and/or report the bug\n");
    return(-1);
  }
  assert(num);
  _MMG5_delEdge(mesh,num);

  /*check if tr iel has other edge*/
  if(pt->edg[ia]) {
    assert(a);
    ped = &mesh->edge[pt->edg[ia]];
    if(ped->a==pib)
      ped->a = pia;
    else {
      assert(ped->b==pib);
      ped->b=pia;
    }
    mesh->tria[a].edg[voy]=pt->edg[ia];
  }
  if(pt->edg[ib]) {
    assert(a1);
    ped = &mesh->edge[pt->edg[ib]];
    mesh->tria[a1].edg[v1]=pt->edg[ib];
  }
  _MMG2D_delElt(mesh,iel);
  memcpy(ppb->c,coor,2*sizeof(double));
  memcpy(&sol->m[sol->size*pib],solu,sol->size*sizeof(double));

  _MMG5_SAFE_FREE(cal);
  _MMG5_SAFE_FREE(list);

  // if ( !MG2_chkmsh(mesh,0) ) exit(EXIT_FAILURE);

  return(1);
}
