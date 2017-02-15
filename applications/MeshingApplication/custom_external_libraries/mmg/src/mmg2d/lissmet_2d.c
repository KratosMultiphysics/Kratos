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
 * \file mmg2d/lissmet_2d.c
 * \brief Size gradation functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"

/**
 * Anisotropic gradation (h-gradation procedure). See:
 * http://www.ann.jussieu.fr/frey/publications/ijnme4398.pdf
 */
int lissmet_ani(MMG5_pMesh mesh,MMG5_pSol sol) {
  HashTable      edgeTable;
  Hedge         *pht;
  MMG5_pTria     pt;
  MMG5_pPoint    p1,p2;
  double         hsiz,logh,logs,*ma,*mb,ux,uy,d1,d2,dd,rap,dh;
  double         tail,coef,coef1,coef2,ma1[3],mb1[3],m[3],dd1,dd2;
  int            i,nc,k,itour,maxtou,ncor,a,b,iadr;
  double         SQRT3DIV2=0.8660254037844386;

  hsiz   = mesh->info.hgrad;
  logh   = log(hsiz);
  logs   = 0.001 + logh;
  maxtou = 100;
  ncor   = 0;
  itour  = 0;

  /*alloc hashtable*/
  //#warning revoir le dimensionnement!!!!
  edgeTable.size  = mesh->ntmax;
  edgeTable.nxtmax = 3*mesh->ntmax+1;
  edgeTable.hnxt  = mesh->ntmax;
  _MMG5_SAFE_CALLOC(edgeTable.item,edgeTable.nxtmax,Hedge);

  memset(edgeTable.item,0,edgeTable.nxtmax*sizeof(Hedge));
  for (k=edgeTable.size; k<edgeTable.nxtmax; k++)
    edgeTable.item[k].nxt = k+1;

  /*build edge table*/
  //#warning optimiser!
  for(k=1 ; k<=mesh->nt ; k++) {
    pt = &mesh->tria[k];
    for(i=0 ; i<3 ; i++)
      MMG2_hashEdge(&edgeTable,k,pt->v[MMG2_iare[i][0]],pt->v[MMG2_iare[i][1]]);
  }

  /* reset color */
  for (k=1; k<=mesh->np; k++)
    (&mesh->point[k])->tagdel = mesh->base+1;

  /* analyze mesh edges via hash table */
  do {
    ++mesh->base;
    nc = 0;
    for (k=0; k<edgeTable.size; k++) {
      pht = &edgeTable.item[k];
      /* analyze linked list */
      while ( pht ) {
        if ( !pht->min )  break;
        a  = pht->min;
        b  = pht->max;
        p1 = &mesh->point[a];
        p2 = &mesh->point[b];
        iadr = a*sol->size;
        ma   = &sol->m[iadr];
        iadr = b*sol->size;
        mb   = &sol->m[iadr];

        if ( (p1->tagdel < mesh->base) && (p2->tagdel < mesh->base) ) {
          pht = pht->nxt ? &edgeTable.item[pht->nxt] : 0;
          continue;
        }

        /* compute edge lengths */
        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];

        d1 = ma[0]*ux*ux + ma[2]*uy*uy + 2.0*ma[1]*ux*uy;
        assert(d1 >=0);
        if ( d1 < 0.0 )  d1 = 0.0;
        dd1 = M_MAX(EPSD,sqrt(d1));

        d2 = mb[0]*ux*ux + mb[2]*uy*uy+ 2.0*mb[1]*ux*uy;
        assert(d2 >=0);
        if ( d2 < 0.0 )  d2 = 0.0;
        dd2 = M_MAX(EPSD,sqrt(d2));

        /* swap vertices */
        if ( dd1 > dd2 ) {
          p1   = &mesh->point[b];
          p2   = &mesh->point[a];
          mb   = ma;
          iadr = b*sol->size;
          ma   = &sol->m[iadr];
          dd   = dd1;
          dd1  = dd2;
          dd2  = dd;
        }
        rap = dd2 / dd1;
        dh = rap - 1.0;
        if ( fabs(dh) > EPSD ) {
          // Edge length in the metric
          tail = (dd1+dd2+4*sqrt(0.5*(d1+d2))) / 6.0;
          coef = log(rap) / tail;
          p1->tagdel = mesh->base+1;
          p2->tagdel = mesh->base+1;

          /* update sizes */
          if ( coef > logs ) {
            coef      = exp(tail*logh);
            p1->tagdel = mesh->base;
            p2->tagdel = mesh->base;

            /*coef1 = exp(dd2*logh);*/
            coef1 = 1.0 + logh*dd2;
            coef1 = 1.0 / (coef1*coef1);
            /*coef2 = exp(dd1*logh);*/
            coef2 = 1.0 + logh*dd1;
            coef2 = 1.0 / (coef2*coef2);

            /* metric intersection */
            coef = 1.0 / (coef*coef);
            for (i=0; i<3; i++) {
              ma1[i] = coef * ma[i];
              mb1[i] = coef * mb[i];
            }

            if ( _MMG5_intersecmet22(mesh,ma,mb1,m) ) {
              for (i=0; i<3; i++)  ma[i] = m[i];
            }
            else {
              for (i=0; i<3; i++)  ma[i]  = SQRT3DIV2 * (ma[i]+mb1[i]);
            }
            if ( _MMG5_intersecmet22(mesh,ma1,mb,m) ) {
              for (i=0; i<3; i++)  mb[i] = m[i];
            }
            else {
              for (i=0; i<3; i++)  mb[i] = SQRT3DIV2 * (mb[i]+ma1[i]);
            }
            nc++;
          }
        }
        /* next edge */
        pht = pht->nxt ? &edgeTable.item[pht->nxt] : 0;
      }
    }
    ncor += nc;
  } while ( nc && ++itour < maxtou );
  _MMG5_SAFE_FREE(edgeTable.item);

  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"    gradation: %7d updated, %d iter\n",ncor,itour);
  }

  return(1);
}

int lissmet_iso(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria     ptt;
  MMG5_pPoint    p1,p2;
  double         logh,logs,dd,rap,dh;
  double    tail,coef;
  int   i,nc,k,maxtou;
  int   it,i1,i2;
  double hmax,ax,ay,lograp;

  it     = 0;
  logh   = log(mesh->info.hgrad);
  logs   = 0.01 + logh;
  maxtou = 1000;

  do {
    nc   = 0;
    hmax = 0.;
    for(k=1 ; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      if(!M_EOK(ptt)) continue;
      for(i=0 ; i<3 ; i++) {
        i1 = ptt->v[MMG2_iare[i][0]];
        i2 = ptt->v[MMG2_iare[i][1]];
        /* swap vertices */
        if ( sol->m[i1] > sol->m[i2] ) {
          i1 = ptt->v[MMG2_iare[i][1]];
          i2 = ptt->v[MMG2_iare[i][0]];
        }
        p1 = &mesh->point[i1];
        p2 = &mesh->point[i2];

        /* compute edge size */
        ax = p2->c[0] - p1->c[0];
        ay = p2->c[1] - p1->c[1];

        rap = sol->m[i2] / sol->m[i1];
        dh  = rap - 1.0f;
        dd  = sqrt(ax*ax + ay*ay );
        if ( fabs(dh) > 1e-6 ) {
          lograp = log(rap);
          tail   = dd * dh / (sol->m[i2]*lograp);
          coef   = lograp / tail;

          /* update sizes */
          if ( coef > logs ) {
            sol->m[i2]  = sol->m[i1] * exp(tail*logh);
            hmax = M_MAX(hmax,sol->m[i2]);
            nc++;
          }
        }
      }
    }
  } while(nc && it++ < maxtou);
  return(1);
}
