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
 * \file mmg2d/locate_2d.c
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"
#define EPST -1e-18
#define EPSR  1e+18
#define EPSNULL 1e-12
#define EPSNULL2 5e-13



/*renvoie les coor bary de P(c[0],c[1]) dans tria pt ainsi que le det*/
void MMG2_coorbary(MMG5_pMesh mesh,MMG5_pTria pt,double c[2],double* det,double* l1,double* l2) {
  MMG5_pPoint      pt1,pt2,pt3;
  double      a11,a12,a21,a22;
//  double      ax,ay,bx,by,cx,cy;
  double      x,y;

  pt1 = &mesh->point[pt->v[0]];
  pt2 = &mesh->point[pt->v[1]];
  pt3 = &mesh->point[pt->v[2]];
  /*calcul des coor bary dans k*/
  *det = (pt1->c[0]-pt2->c[0])*(pt1->c[1]-pt3->c[1]) -
    (pt1->c[1]-pt2->c[1])*(pt1->c[0]-pt3->c[0]);
  *det = 1./(*det);

  a11 = pt2->c[1]-pt3->c[1];
  a12 = -(pt2->c[0]-pt3->c[0]);
  a21 = -(pt1->c[1]-pt3->c[1]);
  a22 = pt1->c[0]-pt3->c[0];
  *l1 = a11*(c[0]-pt3->c[0]) + a12*(c[1]-pt3->c[1]);
  *l2 = a21*(c[0]-pt3->c[0]) + a22*(c[1]-pt3->c[1]);


  /* ax = pt1->c[0] - c[0]; */
  /* ay = pt1->c[1] - c[1]; */
  /* bx = pt2->c[0] - c[0]; */
  /* by = pt2->c[1] - c[1]; */
  /* cx = pt3->c[0] - c[0]; */
  /* cy = pt3->c[1] - c[1]; */

  x = (pt3->c[1] - pt1->c[1])*(c[0]-pt1->c[0]) -  (pt3->c[0] - pt1->c[0])*(c[1]-pt1->c[1]);
  y = -(pt2->c[1] - pt1->c[1])*(c[0]-pt1->c[0]) +  (pt2->c[0] - pt1->c[0])*(c[1]-pt1->c[1]);
  x*=(*det);
  y*=(*det);
  *l1 = 1-(x+y);
  *l2 = x;
  /* printf("coor %e %e\n",c[0],c[1]); */
  /* printf("coorP2 : %e %e\n",(*det)*((pt3->c[1] - pt1->c[1])*(pt2->c[0]-pt1->c[0])  */
  /*                                   -  (pt3->c[0] - pt1->c[0])*(pt2->c[1]-pt1->c[1])), */
  /*        (*det)*(-(pt2->c[1] - pt1->c[1])*(pt2->c[0]-pt1->c[0]) */
  /*                +  (pt2->c[0] - pt1->c[0])*(pt2->c[1]-pt1->c[1]))); */
  /* printf("coorP3 : %e %e\n",(*det)*((pt3->c[1] - pt1->c[1])*(pt3->c[0]-pt1->c[0])  */
  /*                                   -  (pt3->c[0] - pt1->c[0])*(pt3->c[1]-pt1->c[1])), */
  /*        (*det)*(-(pt2->c[1] - pt1->c[1])*(pt3->c[0]-pt1->c[0]) */
  /*                +  (pt2->c[0] - pt1->c[0])*(pt3->c[1]-pt1->c[1]))); */
  /* printf("det %e -- %e %e -- %e %e\n",*det,x,y,*l1,*l2); */
  /* printf("l3 %e %e\n",y,1. -(*l1+*l2)); */

}

int MMG2_isInTriangle(MMG5_pMesh mesh,int k,double c[2]) {
  MMG5_pTria pt;
  double det,l1,l2,l3;

  pt = &mesh->tria[k];
  if(!M_EOK(pt)) return(0);

  MMG2_coorbary(mesh,pt,&c[0],&det,&l1,&l2);
  l3 = 1. - (l1+l2);
  if(l3>EPST && l1>EPST && l2>EPST) return(k);
  else
    return(0);

}


int MMG2_cutEdge(MMG5_pMesh mesh,MMG5_pTria pt,MMG5_pPoint ppa,MMG5_pPoint ppb,int ivert) {
  double      la[3],lb[3],det;
  int         icompt,i,ireturn;

  MMG2_coorbary(mesh,pt,ppa->c,&det,&la[0],&la[1]);
  la[2] = 1.-(la[0]+la[1]);
  MMG2_coorbary(mesh,pt,ppb->c,&det,&lb[0],&lb[1]);
  lb[2] = 1.-(lb[0]+lb[1]);
  if(mesh->info.ddebug) printf("barya %e %e %e\n",la[0],la[1],la[2]); 
  if(mesh->info.ddebug) printf("baryb %e %e %e\n",lb[0],lb[1],lb[2]); 
  //if(ddebug) exit(EXIT_FAILURE);
  /*if ia ou ib sommets de pt*/
  if(ivert==1) { /*ia sommets*/
    for(i=0 ; i<3 ; i++) {
      if(fabs(la[i]-1.)<1e-12) {
        if(lb[i]<0) return(i+1);
        else return(0);
      }
    }
  } else if(ivert==2) { /*ib sommets*/
    for(i=0 ; i<3 ; i++) {
      if(fabs(lb[i]-1.)<1e-12) {
        if(la[i]<0) return(i+1);
        else return(0);
      }
    }
  }


  icompt = 0;
  for(i=0 ; i<3 ; i++) {
    if((la[i]>=0 && lb[i]<=0) || (la[i]<=0 && lb[i]>=0)) {
      ireturn = i+1;
      icompt++;
    }
  }
  if(mesh->info.ddebug) printf("coor bary %e %e %e\n",la[0],la[1],la[2]);
  if(mesh->info.ddebug)printf("coor bary %e %e %e\n",lb[0],lb[1],lb[2]);

  if(icompt>1) return(ireturn);
  return(0);
}

/*return i>0 if Edge ia-ib intersecte tr k et return -3 if edge ia-ib is in k*/
int MMG2_cutEdgeTriangle(MMG5_pMesh mesh,int k,int ia,int ib) {
  MMG5_pTria   pt;
  MMG5_pPoint  pt1,pt2,pt3,ppa,ppb;
  double  a11,a21,a12,a22,aire1,aire2,aire3,prod1,prod2,prod3;
  int     ibreak,iare,i,ivert;

  ppa = &mesh->point[ia];
  ppb = &mesh->point[ib];

  pt = &mesh->tria[k];
  if(!pt->v[0]) return(0);
  ibreak = 0;
  ivert  = 0;
  if(pt->v[0]==ib || pt->v[1]==ib || pt->v[2]==ib) {
    ibreak = 1;
    ivert  = 2;
  }
  if(pt->v[0]==ia || pt->v[1]==ia || pt->v[2]==ia) ivert = 1;


  pt1 = &mesh->point[pt->v[0]];
  pt2 = &mesh->point[pt->v[1]];
  pt3 = &mesh->point[pt->v[2]];

  /*calcul des aire iaibPi*/
  a11 = ppb->c[0] - ppa->c[0];
  a21 = ppb->c[1] - ppa->c[1];
  a12 = pt1->c[0] - ppa->c[0];
  a22 = pt1->c[1] - ppa->c[1];
  aire1 = a11*a22 - a12*a21;

  a12 = pt2->c[0] - ppa->c[0];
  a22 = pt2->c[1] - ppa->c[1];
  aire2 = a11*a22 - a12*a21;

  a12 = pt3->c[0] - ppa->c[0];
  a22 = pt3->c[1] - ppa->c[1];
  aire3 = a11*a22 - a12*a21;

  prod1 = aire1*aire2;
  prod2 = aire3*aire2;
  prod3 = aire3*aire1;

  if ( prod1 > 0 && ((prod2 < 0 || prod3 < 0))) { /*le tr est coupe par la droite ia-ib*/
    if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
      return(iare);
    }
  }

  if ( prod2 > 0 && ((prod1 < 0 || prod3 < 0))) {
    if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
      return(iare);
    }
  }
  if ( prod3 > 0 && ((prod2 < 0 || prod1 < 0))) {
    if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
      return(iare);
    }
  }
  /*sommet == pt arete*/
  for(i=0 ; i<3 ; i++){
    if(pt->v[i]==ia || ibreak) {
      if((prod1 < 0) || (prod2 < 0) || (prod3 < 0)) {
        if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
          return(iare);
        }
      } else {
        /*check if ia-ib edge de pt*/
        if(ibreak && (pt->v[(i+1)%3]==ia || pt->v[(i+2)%3]==ia)) {
          return(-3);
        }
        else if(pt->v[i]==ia && (pt->v[(i+1)%3]==ib || pt->v[(i+2)%3]==ib)) {
          return(-3);
        }
      }
    }
  }

  return(0);
}


//cherche un des triangles contenant le point ip
int MMG2_findTria(MMG5_pMesh mesh,int ip) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p0,p1,p2;
  int       find,iel,base,iadr,*adja,isign,k;
  double    ax,ay,bx,by,dd,epsra,cx,cy,aire1,aire2,aire3;

  ppt  = &mesh->point[ip];
  ++mesh->base;
  base = ++mesh->base;
  find = 0;
  iel  = 1;
  do {
    pt = &mesh->tria[iel];
    if ( !M_EOK(pt) )  {
      iel++;
      if(iel > mesh->nt) {
        //fprintf(stdout,"exaustif search\n");
        return(0);
      }
      continue;
    }
    if ( pt->base == base )  {
      //printf("gloups  pour %d base %d -- %d\n",iel,base,ip);
      fprintf(stdout,"Warning numerical problem in findTria, please make a bug report\n");
      return(iel);
    }
    /*check not vertex*/
    if (pt->v[0]==ip || pt->v[1]==ip || pt->v[2]==ip) return(iel);
    pt->base = base;
    iadr = 3*(iel-1)+1;
    adja = &mesh->adja[iadr];
    /*calcul des coor bary de ip dans le tr iel*/
    /*det(BM,CM)AM + det(CM,AM)BM + det(AM,BM)CM = 0*/
    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];

    ax = p1->c[0] - p0->c[0];
    ay = p1->c[1] - p0->c[1];
    bx = p2->c[0] - p0->c[0];
    by = p2->c[1] - p0->c[1];
    dd = ax*by - ay*bx;
    isign = dd > 0.0 ? 1 : -1;
    epsra = isign * (EPST*dd);
    /* barycentric */
    bx = p1->c[0] - ppt->c[0];
    by = p1->c[1] - ppt->c[1];
    cx = p2->c[0] - ppt->c[0];
    cy = p2->c[1] - ppt->c[1];

    /* p in 1 */
    aire1 = isign*(bx*cy - by*cx);
    if ( aire1 < epsra && adja[0]) {
      iel = adja[0] / 3;
      continue;
    }

    ax = p0->c[0] - ppt->c[0];
    ay = p0->c[1] - ppt->c[1];
    aire2 = isign*(cx*ay - cy*ax);
    if ( aire2 < epsra && adja[1]) {
      iel = adja[1] / 3;
      continue;
    }

    //aire3 = -epsra*EPSR - (aire1 + aire2);
    /*try to be more robust...*/
    // aire3 = M_MAX(aire3,isign*(bx*ay - by*ax));
    aire3 = (isign*(ax*by - ay*bx));
    if ( aire3 < epsra && adja[2]) {
      iel = adja[2] / 3;
      continue;
    }
    find = 1;
    /*aire1 = M_MAX(aire1,0.0);
      aire2 = M_MAX(aire2,0.0);
      aire3 = M_MAX(aire3,0.0);

      dd = aire1+aire2+aire3;
      if ( dd != 0.0 )  dd = 1.0 / dd;
      cb[0] = aire1 * dd;
      cb[1] = aire2 * dd;
      cb[2] = aire3 * dd;*/

  } while (!find);

  /*exhaustive search*/
  for (k=1 ; k<=mesh->nt ; k++) {
    pt = &mesh->tria[k];
    if(!M_EOK(pt)) continue;
    if (pt->v[0]==ip || pt->v[1]==ip || pt->v[2]==ip) break;
  }
  if(k<=mesh->nt) {
    return(k);
  }
  return(0);
}



/*list tous les tr intersectes par une arete*/
int MMG2_locateEdge(MMG5_pMesh mesh,int ia,int ib,int* kdep,int* list) {
  MMG5_pTria     pt;
  MMG5_pPoint    pt1,pt2,pt3,ppa,ppb,pt4;
  double    a[3],a11,a21,a12,a22,aire1,aire2,aire3,prod1,prod2,prod3;
  double     niaib,npti;
  int       iadr,*adja,k,ibreak,i,ncompt,lon,iare,ivert;
  //int       ktemp;

  k = *kdep;
  ncompt = 0;
  ibreak = 0;
  lon = 0;
  ppa = &mesh->point[ia];
  ppb = &mesh->point[ib];

  pt = &mesh->tria[k];

  ivert = 0;
  if(pt->v[0]==ia || pt->v[1]==ia || pt->v[2]==ia) ivert = 1;

  if ( !ivert ) {

    if ( !(k = MMG2_findTria(mesh,ia) ) ) {
       return 0;
    }
    *kdep = k;
  }

  if ( mesh->info.ddebug || mesh->info.imprim > 6 )
    printf(" Try to enforce edge %d %d\n",ia,ib);

  mesh->base += 2;
  do {
    pt = &mesh->tria[k];

    pt->base = mesh->base;
    iadr = 3*(k-1)+1;
    adja = &mesh->adja[iadr];
    ibreak = 0;
    ncompt++;
    ivert = 0;
    if(pt->v[0]==ib || pt->v[1]==ib || pt->v[2]==ib) {
      ibreak = 1;
      ivert  = 2;
    }
    if(pt->v[0]==ia || pt->v[1]==ia || pt->v[2]==ia) ivert = 1;

    if ( ibreak == 1 && ivert == 1 ) return 4;

    pt1 = &mesh->point[pt->v[0]];
    pt2 = &mesh->point[pt->v[1]];
    pt3 = &mesh->point[pt->v[2]];

    /*calcul des aire ia-ib-Pi*/
    a11 = ppb->c[0] - ppa->c[0];
    a21 = ppb->c[1] - ppa->c[1];
    a12 = pt1->c[0] - ppa->c[0];
    a22 = pt1->c[1] - ppa->c[1];
    aire1 = a11*a22 - a12*a21;

    a12 = pt2->c[0] - ppa->c[0];
    a22 = pt2->c[1] - ppa->c[1];
    aire2 = a11*a22 - a12*a21;

    a12 = pt3->c[0] - ppa->c[0];
    a22 = pt3->c[1] - ppa->c[1];
    aire3 = a11*a22 - a12*a21;

    prod1 = aire1*aire2;
    prod2 = aire3*aire2;
    prod3 = aire3*aire1;

    a[0] = aire1;
    a[1] = aire2;
    a[2] = aire3;

    if ( prod1 > 0 && ((prod2 < 0 || prod3 < 0))) {
      /* ia-ib cut the 0th edge of the tria */
      if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
        pt->base = mesh->base+1;
        list[lon++] = 3*k + iare-1;
      }
      k = adja[0]/3;
      if ((mesh->tria[k].base>=mesh->base) || !k) {
        k = adja[1]/3;
        if (!iare && (mesh->tria[k].base>=mesh->base || !k)) {
          k = adja[2]/3;
          if((mesh->tria[k].base>=mesh->base)) k = 0;
        } else if((mesh->tria[k].base>=mesh->base)) k = 0;
      }
      if(ibreak) break;
      continue;
    }
    if ( prod2 > 0 && ((prod1 < 0 || prod3 < 0))) {
      if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
        pt->base = mesh->base+1;
        list[lon++] = 3*k + iare-1;
      }
      k = adja[1]/3;
      if ((mesh->tria[k].base>=mesh->base) || !k) {
        k = adja[2]/3;
        if (!iare && (!k || mesh->tria[k].base>=mesh->base)) {
          k = adja[0]/3;
          if((mesh->tria[k].base>=mesh->base)) k = 0;
        } else if((mesh->tria[k].base>=mesh->base)) k = 0;
      }
      if(ibreak) break;
      continue;
    }
    if ( prod3 > 0 && ((prod2 < 0 || prod1 < 0))) {
      if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
        pt->base = mesh->base+1;
        list[lon++] = 3*k + iare-1;
      }
      k = adja[2]/3;
      if ((mesh->tria[k].base>=mesh->base) || !k) {
        k = adja[0]/3;
        if (!iare && (!k || mesh->tria[k].base>=mesh->base)) {
          k = adja[1]/3;
          if((mesh->tria[k].base>=mesh->base)) k = 0;
        } else if((mesh->tria[k].base>=mesh->base)) k = 0;
      }
      if(ibreak) break;
      continue;
    }


    /*sommet == pt arete ou pts alignes avec arete*/
    for(i=0 ; i<3 ; i++){
      iare = 0;
      if(pt->v[i]==ia || ibreak) {
        if((prod1 < 0) || (prod2 < 0) || (prod3 < 0)) {
          if ((iare = MMG2_cutEdge(mesh,pt,ppa,ppb,ivert))) {
            pt->base = mesh->base+1;
            list[lon++] = 3*k + iare-1;
          }
          else ibreak = 0;
        } else {
          /*check if ia-ib edge de pt*/
          if(ibreak && (pt->v[(i+1)%3]==ia || pt->v[(i+2)%3]==ia) ){
            pt->base = mesh->base+1;
            list[lon++] = 3*k;
            ibreak = 3;
          }
          else if(pt->v[i]==ia && (pt->v[(i+1)%3]==ib || pt->v[(i+2)%3]==ib)) {
            pt->base = mesh->base+1;
            list[lon++] = 3*k;
            ibreak = 3;
          }
          else if(fabs(prod1)<EPSNULL2 && fabs(prod2)<EPSNULL2 && fabs(prod3)<EPSNULL2) {
            if((a[MMG2_inxt[i]]<0 && a[MMG2_inxt[MMG2_inxt[i]]]>0)
               || (a[MMG2_inxt[i]]>0 && a[MMG2_inxt[MMG2_inxt[i]]]<0)) {
              // printf("l'arete coupe le tr %e %e\n",a[MMG2_inxt[i]], a[MMG2_inxt[MMG2_inxt[i]]]);
              pt->base = mesh->base+1;
              list[lon++] = 3*k;
              ibreak = 3;
            } else if(a[MMG2_inxt[i]]>0 && a[MMG2_inxt[MMG2_inxt[i]]]>0){
              k = adja[MMG2_inxt[MMG2_inxt[i]]]/3;
              ibreak = 1;
              break;
            } else {
              //calcul de ||iaib|| et de ||ptiib|| avec aire(iaibpti)==0
              niaib = sqrt(a11*a11+a21*a21 );
              if(fabs(a[MMG2_inxt[i]])>EPSNULL) {
                pt4 = &mesh->point[pt->v[MMG2_inxt[MMG2_inxt[i]]]];
              } else {
                pt4 = &mesh->point[pt->v[MMG2_inxt[i]]];
              }
              npti = sqrt((ppb->c[0]-pt4->c[0])*(ppb->c[0]-pt4->c[0])+
                          (ppb->c[1]-pt4->c[1])*(ppb->c[1]-pt4->c[1]));
              if(niaib > npti) {
                //on rajoute le triangle
                pt->base = mesh->base+1;
                list[lon++] = 3*k;
                // ibreak = 3;
              }
              k = adja[MMG2_inxt[i]]/3;
              ibreak = 1;
              break;
            }
            /*            k = adja[ktemp]/3;
                          if (!k || mesh->tria[k].base>=mesh->base) {
                          k = adja[(ktemp+1)%3]/3;
                          if(!k || mesh->tria[k].base>=mesh->base) {
                          k = adja[(ktemp+2)%3]/3;
                          if(mesh->tria[k].base>=mesh->base) k=0;
                          }
                          }
                          ibreak=-10;*/
            break;
          } /*end if(fabs(prod1)<EPSNULL2 && fabs(prod2)<EPSNULL2 && fabs(prod3)<EPSNULL2)*/
          else { /*on choisit de passer par l'arete iaPi si aire(iaibPi) >0*/
            assert(pt->v[i]==ia);
            if(a[MMG2_inxt[i]] > 0) {
              iare=MMG2_inxt[MMG2_inxt[i]];
            } else {
              iare=MMG2_inxt[i];
            }
            k = adja[iare]/3;
            if(mesh->info.ddebug) printf("on trouve adj %d (%d\n)\n",iare,k);
            ibreak = 1;
            break;
          }
        } /*end else de if((prod1 < 0) || (prod2 < 0) || (prod3 < 0))*/
        // if(mesh->info.ddebug) printf("pourquoi on passe pas la!!!!!!!!!");
        // if(ddebug) printf("iare %d\n",iare);
        if(iare) {
          //ktemp = k;
          k = adja[i]/3;
          if (!k || mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) {
            k = adja[(i+1)%3]/3;
            if(!k || mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) {
              k = adja[(i+2)%3]/3;
              if(mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) k=0;
            }
          }
        } else {
          //ktemp = k;
          k = adja[(i+1)%3]/3;
          if (!k || mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) {
            k = adja[(i+2)%3]/3;
            if(!k || mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) {
              k = adja[(i+3)%3]/3;
              if(mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) k=0;
            }
          }
        }
        ibreak++;
        break;
      }/*end if ia || ibreak;*/
    }
    //if(ddebug) exit(EXIT_FAILURE);
    if(ibreak==1 || ibreak==-10) continue;
    if(ibreak>1) break;
    /*a-t-on un pts sur l'arete iaib ?*/
    if (fabs(aire1) < EPSNULL || fabs(aire2) < EPSNULL || fabs(aire3) < EPSNULL) {
      fprintf(stdout,"  ## Error: unexpected failure."
              " Check your initial data and/or report the bug\n");
      exit(EXIT_FAILURE);
    }

    k = adja[0]/3;

    if ((mesh->tria[k].base>=mesh->base) || !k || !mesh->tria[k].v[0]) {
      k = adja[1]/3;
      if ((mesh->tria[k].base>=mesh->base) || !k || !mesh->tria[k].v[0]) {
        k = adja[2]/3;
        if((mesh->tria[k].base>=mesh->base) || !k || !mesh->tria[k].v[0]) {
          k=0;
        }
      }
    }

  } while (ncompt < mesh->nt);

  assert(ibreak);

  lon = (ibreak==4)?4:((-1)*lon);

  return(lon);
}
