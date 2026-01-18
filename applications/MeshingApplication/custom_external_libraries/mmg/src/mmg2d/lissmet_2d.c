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
 * \file mmg2d/lissmet_2d.c
 * \brief Size gradation functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \warning unused
 */
#include "libmmg2d_private.h"

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the metric
 *
 * \return 0 if fail, 1 if success
 *
 * Anisotropic gradation (h-gradation procedure). See:
 * \cite borouchaki1998mesh. The Hc-correction method is used (gradation with
 * respect to H-shock measure). Skip edges with a required extremity (treated
 * in lissmetreq_ani).
 */
int lissmet_ani(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_Hash      edgeTable;
  MMG5_hedge     *pht;
  MMG5_pTria     pt;
  MMG5_pPoint    p1,p2;
  double         logh,logs,*ma,*mb,ux,uy,d1,d2,dd,rap,dh;
  double         tail,coef,ma1[3],mb1[3],m[3],dd1,dd2;
  double         SQRT3DIV2=0.8660254037844386;
  int            i,itour,maxtou;
  MMG5_int       ncor,nc,k,iadr,a,b;
  int8_t         ier;
  static int8_t  mmgWarn = 0;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  ** Grading mesh\n");
  }

  MMG5_mark_pointsOnReqEdge_fromTria ( mesh );

  logh   = mesh->info.hgrad;
  logs   = 0.001 + logh;
  maxtou = 100;
  ncor   = 0;
  itour  = 0;

  /* alloc hashtable */
  if ( !MMG5_hashNew(mesh,&edgeTable,mesh->ntmax,3*mesh->ntmax) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate hash table.\n",__func__);
    return 0;
  }

  /* build edge table */
  for(k=1 ; k<=mesh->nt ; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) {
      continue;
    }
    for(i=0 ; i<3 ; i++) {
      a = pt->v[MMG2D_iare[i][0]];
      b = pt->v[MMG2D_iare[i][1]];

      /* Skip edges with a required vertex */
      if ( mesh->point[a].s || mesh->point[b].s ) {
        continue;
      }
      ier = MMG5_hashEdge(mesh,&edgeTable,a,b,k);
      if ( !ier ) {
        if ( !mmgWarn ) {
          mmgWarn = 1;
          fprintf(stderr,"\n  ## Warning: %s: unable to hash at least one edge"
                  " (tria %" MMG5_PRId ", edge %d).\n",__func__,MMG2D_indElt(mesh,k),i);
        }
      }
    }
  }

  /* reset color */
  for (k=1; k<=mesh->np; k++)
    (&mesh->point[k])->tagdel = mesh->base+1;

  /* analyze mesh edges via hash table */
  do {
    ++mesh->base;
    nc = 0;
    for (k=0; k<edgeTable.siz; k++) {
      pht = &edgeTable.item[k];
      /* analyze linked list */
      while ( pht ) {
        if ( !pht->a )  break;
        a  = pht->a;
        b  = pht->b;
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
        dd1 = M_MAX(MMG2D_EPSD,sqrt(d1));

        d2 = mb[0]*ux*ux + mb[2]*uy*uy+ 2.0*mb[1]*ux*uy;
        assert(d2 >=0);
        if ( d2 < 0.0 )  d2 = 0.0;
        dd2 = M_MAX(MMG2D_EPSD,sqrt(d2));

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
        if ( fabs(dh) > MMG2D_EPSD ) {
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

            /* metric intersection */
            coef = 1.0 / (coef*coef);
            for (i=0; i<3; i++) {
              ma1[i] = coef * ma[i];
              mb1[i] = coef * mb[i];
            }

            if ( MMG5_intersecmet22(mesh,ma,mb1,m) ) {
              for (i=0; i<3; i++)  ma[i] = m[i];
            }
            else {
              for (i=0; i<3; i++)  ma[i]  = SQRT3DIV2 * (ma[i]+mb1[i]);
            }
            if ( MMG5_intersecmet22(mesh,ma1,mb,m) ) {
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
  MMG5_SAFE_FREE(edgeTable.item);

  if ( abs(mesh->info.imprim) > 3 && ncor ) {
    fprintf(stdout,"     gradation: %7" MMG5_PRId " updated, %d iter.\n",ncor,itour);
  }

  return 1;
}
