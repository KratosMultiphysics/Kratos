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
 * \file mmg2d/locate_2d.c
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "libmmg2d_private.h"
#define EPST -1e-18

/* Calculate the barycentric coordinates of point P(c[0],c[1]) in tria pt and the associated determinant */
/* l1 = barycentric coordinate with respect to p1, l2 = barycentric coordinate with respect to p2 */
int MMG2D_coorbary(MMG5_pMesh mesh,MMG5_pTria pt,double c[2],double* det,double* l1,double* l2) {
  MMG5_pPoint      p1,p2,p3;
  double           b2,b3;
  static int8_t    mmgWarn0=0;

  p1 = &mesh->point[pt->v[0]];
  p2 = &mesh->point[pt->v[1]];
  p3 = &mesh->point[pt->v[2]];

  /* Calculate det (p2-p1,p3-p1) */
  *det = (p2->c[0]-p1->c[0])*(p3->c[1]-p1->c[1]) - (p2->c[1]-p1->c[1])*(p3->c[0]-p1->c[0]);
  if ( *det < MMG5_EPSD ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 flat triangle. abort.\n",
              __func__);
    }
    return 0;
  }
  *det = 1.0 / (*det);

  b2 = (c[0]-p1->c[0])*(p3->c[1]-p1->c[1]) - (c[1]-p1->c[1])*(p3->c[0]-p1->c[0]);
  b3 = (p2->c[0]-p1->c[0])*(c[1]-p1->c[1]) -(p2->c[1]-p1->c[1])*(c[0]-p1->c[0]);
  b2 *= (*det);
  b3 *= (*det);
  *l1 = 1.0 - (b2+b3);
  *l2 = b2;

  return 1;
}

/** Check whether c lies in triangle k; return k if so, 0 otherwise */
MMG5_int MMG2D_isInTriangle(MMG5_pMesh mesh,MMG5_int k,double c[2]) {
  MMG5_pTria      pt;
  double          det,l1,l2,l3;
  int8_t          ier;

  pt = &mesh->tria[k];
  if ( !MG_EOK(pt) ) return 0;

  ier = MMG2D_coorbary(mesh,pt,&c[0],&det,&l1,&l2);
  if ( !ier )
    return 0;

  l3 = 1.0 - (l1+l2);
  if ( l3>EPST && l1>EPST && l2>EPST) return k;
  else return 0;
}

/* Check whether edge ppa-ppb crosses triangles pt (in the sense that two edges
   of this triangle are crossed by (ppa-ppb), or only one, and ppa or ppb is a
   vertex of pt; if at least one edge is crossed by ia-ib, return i+1, where i
   is the index of one of the crossed edges */
int MMG2D_cutEdge(MMG5_pMesh mesh,MMG5_pTria pt,MMG5_pPoint ppa,MMG5_pPoint ppb) {
  double      la[3],lb[3],det;
  int         icompt,ireturn;
  int8_t      ier,i;

  ier = MMG2D_coorbary(mesh,pt,ppa->c,&det,&la[0],&la[1]);
  if ( !ier ) return 0;
  la[2] = 1.0-(la[0]+la[1]);

  ier = MMG2D_coorbary(mesh,pt,ppb->c,&det,&lb[0],&lb[1]);
  if ( !ier ) return 0;
  lb[2] = 1.0-(lb[0]+lb[1]);

  /* Check whether ppa or ppb is a vertex of pt */
  for (i=0; i<3; i++) {
    if ( fabs(la[i]-1.0) < 1.0e-12 ) {
      if ( lb[i] < 0.0 ) return i+1;
      else return 0;
    }
    if ( fabs(lb[i]-1.0) < 1.0e-12) {
      if ( la[i] < 0.0 ) return i+1;
      else return 0;
    }
  }

  icompt = 0;
  for (i=0; i<3; i++) {
    if ( ( la[i] >= 0.0 && lb[i] <= 0.0 ) || ( la[i] <= 0.0 && lb[i] >= 0.0 ) ) {
      ireturn = i+1;
      icompt++;
    }
  }

  if ( icompt > 1 ) return ireturn;
  return 0;
}

/* Return i+1>0 if Edge ia-ib intersects triangle k at edge i, 0 if
   it does not intersect k, and return -3 if edge ia-ib is one edge of k*/
int MMG2D_cutEdgeTriangle(MMG5_pMesh mesh,MMG5_int k,MMG5_int ia,MMG5_int ib) {
  MMG5_pTria   pt;
  MMG5_pPoint  p1,p2,p3,ppa,ppb;
  double       a11,a21,a12,a22,area1,area2,area3,prod1,prod2,prod3;
  int          ibreak,iare;
  int8_t       i;

  ppa = &mesh->point[ia];
  ppb = &mesh->point[ib];

  pt = &mesh->tria[k];
  if ( !MG_EOK(pt) ) return 0;
  ibreak = 0;

  if ( pt->v[0] == ib || pt->v[1] == ib || pt->v[2] == ib) ibreak = 1;

  p1 = &mesh->point[pt->v[0]];
  p2 = &mesh->point[pt->v[1]];
  p3 = &mesh->point[pt->v[2]];

  /* Calculation of the areas ia,ib,pi*/
  a11 = ppb->c[0] - ppa->c[0];
  a21 = ppb->c[1] - ppa->c[1];

  a12 = p1->c[0] - ppa->c[0];
  a22 = p1->c[1] - ppa->c[1];
  area1 = a11*a22 - a12*a21;

  a12 = p2->c[0] - ppa->c[0];
  a22 = p2->c[1] - ppa->c[1];
  area2 = a11*a22 - a12*a21;

  a12 = p3->c[0] - ppa->c[0];
  a22 = p3->c[1] - ppa->c[1];
  area3 = a11*a22 - a12*a21;

  prod1 = area1*area2;
  prod2 = area3*area2;
  prod3 = area3*area1;

  /* Both edges p2p3 and p1p3 corresponding to prod2 and prod3 are cut by edge ia,ib */
  if ( prod1 > 0 && ((prod2 < 0 || prod3 < 0))) {
    if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
      return iare;
    }
  }

  /* Both edges corresponding to prod1 and prod3 are cut by edge ia,ib */
  if ( prod2 > 0 && ((prod1 < 0 || prod3 < 0))) {
    if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
      return iare;
    }
  }

  /* Both edges corresponding to prod1 and prod2 are cut by edge ia,ib */
  if ( prod3 > 0 && ((prod2 < 0 || prod1 < 0))) {
    if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
      return iare;
    }
  }

  /* Case where one vertex of pt is ia */
  for(i=0; i<3; i++){
    if ( pt->v[i] == ia || ibreak ) {
      /* One vertex is ia, and the opposite edge is frankly crossed */
      if ( (prod1 < 0) || (prod2 < 0) || (prod3 < 0) ) {
        if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
          return iare;
        }
      }
      else {
        /*check if ia-ib edge de pt*/
        if ( ibreak && ( pt->v[(i+1)%3] == ia || pt->v[(i+2)%3] == ia ) ) {
          return -3;
        }
        else if ( pt->v[i] == ia && ( pt->v[(i+1)%3] == ib || pt->v[(i+2)%3] == ib ) ) {
          return -3;
        }
      }
    }
  }

  return 0;
}

/**
 * \param mesh pointer to mesh
 * \param ip point index
 * \return index of one elt containing k or 0 (if no elt is found)
 *
 * Return the index of one triangle containing k.
 */
MMG5_int MMG2D_findTria(MMG5_pMesh mesh,MMG5_int ip) {
  MMG5_pTria    pt,pt1;
  MMG5_int      iel,base,iadr,*adja,iter,ier;
  int           mvDir[3],jel,i;
  double        l1,l2,l3,det,eps;
  static int8_t mmgWarn0 = 0;

  ++mesh->base;
  base = ++mesh->base;
  iter = 0;
  iel  = 1;
  do {
    mvDir[0] = mvDir[1] = mvDir[2] = 0;
    iter++;
    pt = &mesh->tria[iel];
    if ( !MG_EOK(pt) )  {
      iel++;
      if ( iel > mesh->nt ) return 0;
      continue;
    }

    if ( pt->base == base)  {
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Warning: %s: numerical problem, please make"
                " a bug report.\n",__func__);
      }
      return iel;
    }
    /* Check whether ip is one of the vertices of pt */
    if ( (pt->v[0] == ip) || (pt->v[1] == ip) || (pt->v[2] == ip) ) return iel;

    pt->base = base;
    iadr = 3*(iel-1)+1;
    adja = &mesh->adja[iadr];

    /*compute the barycentric coordinates*/
    ier = MMG2D_coorbary(mesh,pt,mesh->point[ip].c,&det,&l1,&l2);
    if ( !ier ) return 0;

    l3 = 1-l1-l2;
    /*warning -1e12 is too strict*/
    eps = (-1e-9)*MMG2D_quickcal(mesh,pt);

    /*find in which direction we can move*/
    if (l1<eps)
      mvDir[0] = 1;
    if (l2<eps)
      mvDir[1] = 1;
    if (l3<eps)
      mvDir[2] = 1;

    /*if all the barycentric coordinates are positive, we find the triangle*/
    if(!mvDir[0] && !mvDir[1] && !mvDir[2]) {
      break;
    }

    /*first, try to move in an optimal direction*/
    iel = 0;
    for(i=0;i<3; i++) {
      if (!mvDir[i]) continue;
      jel = adja[i]/3;

      if (!jel) continue;
      pt1 = &mesh->tria[jel];

      if(pt1->base == mesh->base) continue;
      iel = jel;
      break;
    }

    /*second, the optimal move is not possible, so try to move in one direction*/
    if(i==3) {
      for(i=0;i<3; i++) {
        if (mvDir[i]) continue;
        jel = adja[i]/3;
        if (!jel) continue;
        pt1 = &mesh->tria[jel];
        if(pt1->base == mesh->base) continue;
        iel = jel;
        break;
      }
    }

    /*we already see all the adj triangle so we can do nothing*/
    if (iel == 0) {
      iter = mesh->nt+1;
    }
  } while ( iter<=mesh->nt );

  return iel;
}

/**
 * \param mesh pointer to the mesh
 * \param ia index of first extremity of the edge
 * \param ib index of second extremity of the edge
 * \param kdep pointer to the index of the first element intersecting the edge
 * \param list pointer to the list of elts intersected by the edge
 *
 * \return 4 if the edge exist in the mesh, 0 if fail, ??? otherwise
 *
 * Calculate the list of all the triangles intersected by edge (ia,ib), starting
 * from kdep = one triangle in the ball of ia; \a list starts at index 0 (and
 * ends at index lon-1). It stores 3*k + iare, where \a k is an intersected tria
 * and \a iare is an intersected edge.
 *
 */
int MMG2D_locateEdge(MMG5_pMesh mesh,MMG5_int ia,MMG5_int ib,MMG5_int* kdep,MMG5_int* list) {
  MMG5_pTria         pt;
  MMG5_pPoint        ppt1,ppt2,ppt3,ppt4,ppa,ppb;
  double             a[3],a11,a21,a12,a22,area1,area2,area3,prod1,prod2,prod3;
  double             niaib,npti;
  MMG5_int           iadr,*adja,k,ibreak,i,ncompt,lon,iare,ivert;
  static int8_t      mmgWarn=0;
  //int       ktemp;

  k = *kdep;
  ncompt = 0;
  lon = 0;
  ppa = &mesh->point[ia];
  ppb = &mesh->point[ib];

  pt = &mesh->tria[k];

  ivert = 0;
  if(pt->v[0]==ia || pt->v[1]==ia || pt->v[2]==ia) ivert = 1;

  if ( !ivert ) {

    if ( !(k = MMG2D_findTria(mesh,ia) ) ) {
      return 0;
    }
    *kdep = k;
  }

  if ( mesh->info.ddebug || mesh->info.imprim > 6 )
    printf(" Try to enforce edge %" MMG5_PRId " %" MMG5_PRId "\n",ia,ib);

  mesh->base += 2;
  do {
    pt = &mesh->tria[k];

    pt->base = mesh->base;
    iadr = 3*(k-1)+1;
    adja = &mesh->adja[iadr];
    ibreak = 0;
    ncompt++;

    /* ib is a vertex of the current triangle */
    if ( pt->v[0] == ib || pt->v[1] == ib || pt->v[2] == ib ) ibreak = 1;

    ppt1 = &mesh->point[pt->v[0]];
    ppt2 = &mesh->point[pt->v[1]];
    ppt3 = &mesh->point[pt->v[2]];

    /* Calculate the 3 areas (ia,ib,pi)*/
    a11 = ppb->c[0] - ppa->c[0];
    a21 = ppb->c[1] - ppa->c[1];
    a12 = ppt1->c[0] - ppa->c[0];
    a22 = ppt1->c[1] - ppa->c[1];
    area1 = a11*a22 - a12*a21;

    a12 = ppt2->c[0] - ppa->c[0];
    a22 = ppt2->c[1] - ppa->c[1];
    area2 = a11*a22 - a12*a21;

    a12 = ppt3->c[0] - ppa->c[0];
    a22 = ppt3->c[1] - ppa->c[1];
    area3 = a11*a22 - a12*a21;

    /** prod1, prod2 and prod3 allow to localize the edge regarding the tria:
        - if all are >0, the edge doesn't cross the tria (same sign for all areas)
        - if 2 are negative and 1 is positive, the edge intersect 2 triangle edges
        - if 2 are 0, the edge pass through 1 of the tria vertices
        - if 3 are 0, the edge coincide with 1 of the tria edges

    */
    prod1 = area1*area2;
    prod2 = area3*area2;
    prod3 = area3*area1;

    a[0] = area1;
    a[1] = area2;
    a[2] = area3;

    /** ia-ib intersect 2 triangle edges => 2 products are <0, 1 is >0 */
    /* Both edges p2p3 and p1p3 in pt have franck intersections with edge (ia - ib) */
    if ( prod1 > 0.0 && ((prod2 < 0.0 || prod3 < 0.0))) { /*le tr est coupe par la droite ia-ib*/
      if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
        pt->base = mesh->base+1;
        list[lon++] = 3*k + iare-1;
      }
      k = adja[0] / 3;
      if ( ( mesh->tria[k].base >= mesh->base ) || !k ) {
        k = adja[1] / 3;
        if ( !iare && (mesh->tria[k].base >= mesh->base || !k ) ) {
          k = adja[2] / 3;
          if ( ( mesh->tria[k].base >= mesh->base ) )
            k = 0;
        }
        else if ( ( mesh->tria[k].base >= mesh->base ) )
          k = 0;
      }
      if ( ibreak ) break;
      continue;
    }

    /* Both edges p1p2 and p1p3 in pt have franck intersections with edge (ia - ib) */
    if ( prod2 > 0.0 && ((prod1 < 0.0 || prod3 < 0.0 ))) {
      if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
        pt->base = mesh->base+1;
        list[lon++] = 3*k + iare-1;
      }
      k = adja[1] / 3;
      if ((mesh->tria[k].base >= mesh->base) || !k) {
        k = adja[2] / 3;
        if ( !iare && (!k || mesh->tria[k].base >= mesh->base) ) {
          k = adja[0] / 3;
          if ( ( mesh->tria[k].base >= mesh->base ) )
            k = 0;
        }
        else if ( ( mesh->tria[k].base >= mesh->base ) )
          k = 0;
      }
      if ( ibreak ) break;
      continue;
    }

    /* Both edges p1p2 and p2p3 in pt have franck intersections with edge (ia - ib) */
    if ( prod3 > 0.0 && ((prod2 < 0.0 || prod1 < 0.0 ))) {
      if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
        pt->base = mesh->base+1;
        list[lon++] = 3*k + iare-1;
      }
      k = adja[2] / 3;
      if ( ( mesh->tria[k].base >= mesh->base ) || !k ) {
        k = adja[0] / 3;
        if ( !iare && (!k || mesh->tria[k].base >= mesh->base) ) {
          k = adja[1] / 3;
          if((mesh->tria[k].base >= mesh->base))
            k = 0;
        }
        else if ( (mesh->tria[k].base >= mesh->base) )
          k = 0;
      }
      if( ibreak ) break;
      continue;
    }

    /* Case where ia or ib is one vertex in pt */
    for(i=0; i<3; i++) {
      iare = 0;
      if ( pt->v[i] == ia || ibreak ) {
        if ( (prod1 < 0.0) || (prod2 < 0.0) || (prod3 < 0.0) ) {
          if ( (iare = MMG2D_cutEdge(mesh,pt,ppa,ppb)) ) {
            pt->base = mesh->base+1;
            list[lon++] = 3*k + iare-1;
          }
          else ibreak = 0;
        }
        else {
          /* Check whether ia-ib is an edge in pt*/
          if ( ibreak && ( pt->v[(i+1)%3] == ia || pt->v[(i+2)%3] == ia) ){
            pt->base = mesh->base+1;
            list[lon++] = 3*k;
            ibreak = 3;
          }
          else if ( pt->v[i] == ia && ( pt->v[(i+1)%3] == ib || pt->v[(i+2)%3] == ib ) ) {
            pt->base = mesh->base+1;
            list[lon++] = 3*k;
            ibreak = 3;
          }
          else if ( fabs(prod1)<MMG5_EPSD2 && fabs(prod2)<MMG5_EPSD2 && fabs(prod3)<MMG5_EPSD2) {
            if ( (a[MMG5_inxt2[i]] < 0.0 && a[MMG5_iprv2[i]] > 0.0 )
                 || (a[MMG5_inxt2[i]] > 0.0 && a[MMG5_iprv2[i]] < 0.0 ) ) {

              pt->base = mesh->base+1;
              list[lon++] = 3*k;
              ibreak = 3;
            }
            else if ( a[MMG5_inxt2[i]] > 0.0 && a[MMG5_iprv2[i]] > 0.0 ){
              k = adja[MMG5_iprv2[i]] / 3;
              //ibreak = 1;
              break;
            }
            else {
              //calcul de ||iaib|| et de ||ptiib|| avec aire(iaibpti)==0
              niaib = sqrt(a11*a11+a21*a21 );
              if ( fabs(a[MMG5_inxt2[i]]) > MMG5_EPSD ) {
                ppt4 = &mesh->point[pt->v[MMG5_iprv2[i]]];
              }
              else {
                ppt4 = &mesh->point[pt->v[MMG5_inxt2[i]]];
              }
              npti = sqrt((ppb->c[0]-ppt4->c[0])*(ppb->c[0]-ppt4->c[0])+
                          (ppb->c[1]-ppt4->c[1])*(ppb->c[1]-ppt4->c[1]));
              if ( niaib > npti ) {
                //on rajoute le triangle
                pt->base = mesh->base+1;
                list[lon++] = 3*k;
                // ibreak = 3;
              }
              k = adja [MMG5_inxt2[i]] / 3;
              //ibreak = 1;
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
          } /*end if(fabs(prod1)<MMG5_EPSD2 && fabs(prod2)<MMG5_EPSD2 && fabs(prod3)<MMG5_EPSD2)*/
          else { /*on choisit de passer par l'arete iaPi si aire(iaibPi) >0*/
            assert ( pt->v[i] == ia );
            if ( a[MMG5_inxt2[i]] > 0.0 )
              iare = MMG5_iprv2[i];
            else
              iare = MMG5_inxt2[i];

            k = adja[iare] / 3;
            ibreak = 1;
            break;
          }
        } /*end else de if((prod1 < 0) || (prod2 < 0) || (prod3 < 0))*/

        if ( iare ) {
          //ktemp = k;
          k = adja[i] / 3;
          if (!k || mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) {
            k = adja[(i+1)%3]/3;
            if(!k || mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) {
              k = adja[(i+2)%3]/3;
              if(mesh->tria[k].base>=mesh->base || !mesh->tria[k].v[0]) k=0;
            }
          }
        }
        else {
          //ktemp = k;
          k = adja[(i+1)%3] / 3;
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

    if ( ibreak == 1 || ibreak == -10 ) continue;
    if ( ibreak > 1 ) break;

    /*a-t-on un pts sur l'arete iaib ?*/
    if (fabs(area1) < MMG5_EPSD || fabs(area2) < MMG5_EPSD || fabs(area3) < MMG5_EPSD) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr,"\n  ## Error: %s: unexpected failure."
                " Check your initial data and/or report the bug. lon:%" MMG5_PRId ". %e %e %e\n"
                " tria %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId ",\n",
                __func__,lon, area1,area2,area3,k,mesh->tria[k].v[0],
                mesh->tria[k].v[1],mesh->tria[k].v[2]);
      }
      return 0;
    }

    //ktemp = k;
    /* printf("adj (base) pour le tri %" MMG5_PRId " : %" MMG5_PRId "(%" MMG5_PRId ") %" MMG5_PRId "(%" MMG5_PRId ") %" MMG5_PRId "(%" MMG5_PRId ")\n",k,adja[0]/3,mesh->tria[adja[0]/3].base>=mesh->base */
    /*        ,adja[1]/3,mesh->tria[adja[1]/3].base>=mesh->base,adja[2]/3,mesh->tria[adja[2]/3].base>=mesh->base); */
    k = adja[0] / 3;
    if ((mesh->tria[k].base>=mesh->base) || !k || !mesh->tria[k].v[0]) {
      k = adja[1] / 3;
      if ((mesh->tria[k].base>=mesh->base) || !k || !mesh->tria[k].v[0]) {
        k = adja[2] / 3;
        if((mesh->tria[k].base>=mesh->base) || !k || !mesh->tria[k].v[0]) {
          k=0;
        }
      }
    }

  }
  while ( ncompt < mesh->nt );

  assert ( ibreak );
  lon = ( ibreak == 4 ) ? 4 : ((-1)*lon);
  return lon;
}
