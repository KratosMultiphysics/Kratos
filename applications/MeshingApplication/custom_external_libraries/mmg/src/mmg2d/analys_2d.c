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
 * \file mmg2d/analys_2d.c
 * \brief  Analysis routine for an input mesh without structure
 * passing through a point.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg2d_private.h"

extern int8_t ddb;

/**
 * \param mesh pointer toward the mesh
 * \param init_cc 1 if we need to reinitialized \a cc field of tria because
 * setadj has already been called (isosurf mode)
 *
 * \return 1 if success, 0 if fail
 *
 * Set tags GEO, BDY and REF to triangles and points by traveling the mesh;
 * count number of subdomains or connected components
 *
 */
int MMG2D_setadj(MMG5_pMesh mesh, int8_t init_cc) {
  MMG5_pTria       pt,pt1;
  MMG5_pQuad       pq;
  MMG5_int         *pile,*adja,ipil,k,kk,ncc,ip1,ip2,nr,nref;
  int16_t          tag;
  int8_t           i,ii,i1,i2;

  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  /** Step 1: Tags setting from triangles analysis */
  MMG5_SAFE_MALLOC(pile,mesh->nt+1,MMG5_int,return 0);

  /* Reinitialization of cc field (as we pass in steadj more than once) */
  if ( init_cc ) {
    for ( k=1; k<=mesh->nt; ++k ) {
      mesh->tria[k].cc = 0;
    }
  }

  /* Initialization of the pile */
  ncc = 1;
  pile[1] = 1;
  ipil = 1;
  mesh->tria[1].cc = ncc;
  nr = nref =0;

  while ( ipil > 0 ) {

    while ( ipil > 0 ) {
      k = pile[ipil];
      ipil--;

      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;
      adja = &mesh->adja[3*(k-1)+1];

      for (i=0; i<3; i++) {
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];

        if ( adja[i] ) {
          kk = adja[i] / 3;
          ii = adja[i] % 3;
          pt1 = &mesh->tria[kk];

          if ( (!mesh->info.opnbdy) && (pt->ref==pt1->ref) ) {
            /* Remove BDY/REF and GEO edge tags but required tag must be kept
             * (to preserve required tria info) */
            pt->tag[i] &= ~MG_REF;
            pt->tag[i] &= ~MG_BDY;
            pt->tag[i] &= ~MG_GEO;

            pt1->tag[ii] &= ~MG_REF;
            pt1->tag[ii] &= ~MG_BDY;
            pt1->tag[ii] &= ~MG_GEO;

            pt1->tag[ii] |= pt->tag[i];
            pt->tag[i]   |= pt1->tag[ii];
          }
          else {
            /* Transfert edge tag to adjacent */
            pt1->tag[ii] |= pt->tag[i];
            pt->tag[i] |= pt1->tag[ii];
          }
        }

        /* Transfer tags if i is already an edge (e.g. supplied by the user) */
        if ( MG_EDG(pt->tag[i]) || (pt->tag[i] & MG_REQ) ) {
          mesh->point[ip1].tag |= pt->tag[i];
          mesh->point[ip2].tag |= pt->tag[i];
        }

        /* Case of an external boundary */
        if ( !adja[i] ) {
          tag = ( MG_GEO+MG_BDY );

          if ( ( !(pt->tag[i] & MG_REQ) ) &&
               ( mesh->info.nosurf || (pt->tag[i] & MG_PARBDY ) ) ) {
              pt->tag[i] |= ( tag + MG_REQ + MG_NOSURF );
            }
          else
            pt->tag[i] |= tag;

          if ( ( !(mesh->point[ip1].tag & MG_REQ) ) &&
               ( mesh->info.nosurf || (pt->tag[i] & MG_PARBDY ) ) ) {
            mesh->point[ip1].tag |= ( tag+MG_REQ+MG_NOSURF );
          }
          else {
            mesh->point[ip1].tag |= tag;
          }

          if ( ( !(mesh->point[ip2].tag & MG_REQ) ) &&
               ( mesh->info.nosurf || (pt->tag[i] & MG_PARBDY ) ) ) {
            mesh->point[ip2].tag |= ( tag+MG_REQ+MG_NOSURF );
          }
          else {
            mesh->point[ip2].tag |= tag;
          }
          nr++;
          continue;
        }
        kk = adja[i] / 3;
        ii = adja[i] % 3;
        pt1 = &mesh->tria[kk];

        /* Case of a boundary (except if opnbdy is enabled, the boundary must be
         * between 2 different subdomains) */
        if ( (mesh->info.opnbdy && ( pt->tag[i] || pt1->tag[ii] ) )
             || MMG5_abs(pt1->ref) != MMG5_abs(pt->ref) ) {
          tag = ( MG_REF+MG_BDY );

          if ( !mesh->info.nosurf ) {
            pt->tag[i]   |= tag;
            pt1->tag[ii] |= tag;
            mesh->point[ip1].tag |= tag;
            mesh->point[ip2].tag |= tag;
          }
          else {
            if ( !(pt->tag[i] & MG_REQ) )
              pt->tag[i]   |= ( tag+MG_REQ+MG_NOSURF );
            else
              pt->tag[i]   |= tag;

            if ( !(pt1->tag[ii] & MG_REQ) )
              pt1->tag[ii] |= ( tag+MG_REQ+MG_NOSURF );
            else
              pt1->tag[ii] |= tag;

            if ( !(mesh->point[ip1].tag & MG_REQ ) )
              mesh->point[ip1].tag |= ( tag+MG_REQ+MG_NOSURF );
            else
              mesh->point[ip1].tag |= tag;

            if ( !(mesh->point[ip2].tag & MG_REQ ) )
              mesh->point[ip2].tag |= ( tag+MG_REQ+MG_NOSURF );
            else
              mesh->point[ip2].tag |= tag;
          }
          if ( kk > k ) nref++;
          continue;
        }

        if ( pt1->cc > 0 )  continue;
        ipil++;
        pile[ipil] = kk;
        pt1->cc = ncc;
      }
    }

    /* Find the starting triangle for the next connected component */
    ipil = 0;
    for (kk=1; kk<=mesh->nt; kk++) {
      pt = &mesh->tria[kk];
      if ( pt->v[0] && (pt->cc == 0) ) {
        ipil = 1;
        pile[ipil] = kk;
        ncc++;
        pt->cc = ncc;
        break;
      }
    }
  }
  MMG5_SAFE_FREE(pile);

  /** Step 2: Mark the edges at interface between tria and quads as nosurf and
   * required */
  for ( k=1; k<=mesh->nquad; k++ ) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;

    adja = &mesh->adjq[4*(k-1)+1];
    for (i=0; i<4; i++) {

      kk = adja[i];
      if ( kk >= 0 ) {
        continue;
      }

      /* The edge is at the interface between a quad and a triangle */
      kk = MMG5_abs(kk);
      ii = kk%3;
      pt = &mesh->tria[kk/3];

      if ( !(pt->tag[ii] & MG_REQ) ) {
        pt->tag[ii] |= (MG_REQ + MG_NOSURF + MG_BDY);
      }
      if ( !(pq->tag[i] & MG_REQ) ) {
        pq->tag[i] |= (MG_REQ + MG_NOSURF + MG_BDY);
      }
      assert ( pt->edg[ii] == pq->edg[i] );

      /* Transfer edge tag toward edge vertices */
      i1 = MMG2D_idir_q[i][0];
      i2 = MMG2D_idir_q[i][1];

      if ( !(mesh->point[pq->v[i1]].tag & MG_REQ) ) {
        mesh->point[pq->v[i1]].tag |= (MG_REQ + MG_NOSURF + MG_BDY);
      }
      if ( !(mesh->point[pq->v[i2]].tag & MG_REQ) ) {
        mesh->point[pq->v[i2]].tag |= (MG_REQ + MG_NOSURF + MG_BDY);
      }
    }
  }


  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"     Connected component or subdomains: %" MMG5_PRId " \n",ncc);
    fprintf(stdout,"     Tagged edges: %" MMG5_PRId ",  ridges: %" MMG5_PRId ",  refs: %" MMG5_PRId "\n",nr+nref,nr,nref);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param ref reference of the boundary to analyze (analyze all the boundaries
 * if MMG5_UNSET)
 *
 * \return 1 if success, 0 if fail.
 *
 * Identify singularities in the mesh.
 *
 */
int MMG2D_singul(MMG5_pMesh mesh, MMG5_int ref ) {
  MMG5_pTria          pt;
  MMG5_pPoint         ppt,p1,p2;
  double              ux,uy,uz,vx,vy,vz,dd;
  MMG5_int            list[MMG5_TRIA_LMAX+2],listref[MMG5_TRIA_LMAX+2],k,nm,nre,nc;
  int                 ns,ng,nr;
  int8_t              i;

  nre = nc = nm = 0;

  if ( ref == MMG5_UNSET ) {
    /* reset the ppt->s tag */
    for (k=1; k<=mesh->np; ++k) {
      mesh->point[k].s = 0;
    }
  }
  else {
    /** Mark the points that we don't want to analyze */
    /* First: mark all the points */
    for ( k=1; k<=mesh->np; ++k ) {
      mesh->point[k].s = 1;
    }
    /* Second: if the point belongs to a boundary edge of reference ref, remove
     * the mark */
    for ( k=1; k<=mesh->nt; ++k ) {
      pt = &mesh->tria[k];

      for (i=0; i<3; i++) {
        if ( !MG_EDG(pt->tag[i]) ) continue;
        if ( pt->edg[i] != ref )   continue;

        mesh->point[pt->v[MMG5_iprv2[i]]].s = 0;
        mesh->point[pt->v[MMG5_inxt2[i]]].s = 0;
      }
    }
  }

  /** Singularity identification */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( ! MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];

      if ( ppt->s ) continue;
      ppt->s = 1;

      if ( !MG_VOK(ppt) || MG_CRN & ppt->tag )  continue;

      if ( !MG_EDG(ppt->tag) ) continue;

      ns = MMG5_bouler(mesh,mesh->adja,k,i,list,listref,&ng,&nr,MMG5_TRIA_LMAX);

      if ( !ns )  continue;
      if ( (ng+nr) > 2 ) {
        /* Previous classification may be subject to discussion, and may depend on the user's need */
        ppt->tag |= MG_NOM;
        nm++;
        /* Two ridge curves and one ref curve: non manifold situation */
        /*if ( ng == 2 && nr == 1 ) {
          ppt->tag |= MG_NOM;
          nm++;
          }
          else {
          ppt->tag |= MG_CRN + MG_REQ;
          nre++;
          nc++;
          }*/
      }
      /* One ridge curve and one ref curve meeting at a point */
      else if ( (ng == 1) && (nr == 1) ) {
        ppt->tag |= MG_REQ;
        nre++;
      }
      /* Evanescent ridge */
      else if ( ng == 1 && !nr ){
        ppt->tag |= MG_CRN + MG_REQ;
        nre++;
        nc++;
      }
      /* Evanescent ref curve */
      else if ( nr == 1 && !ng ){
        ppt->tag |= MG_CRN + MG_REQ;
        nre++;
        nc++;
      }
      /* Check ridge angle */
      else {
        assert ( ng == 2 || nr == 2 );
        p1 = &mesh->point[list[1]];
        p2 = &mesh->point[list[2]];
        ux = p1->c[0] - ppt->c[0];
        uy = p1->c[1] - ppt->c[1];
        uz = p1->c[2] - ppt->c[2];
        vx = p2->c[0] - ppt->c[0];
        vy = p2->c[1] - ppt->c[1];
        vz = p2->c[2] - ppt->c[2];
        dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);

        /* If both edges carry different refs, tag vertex as singular */
        if ( listref[1] != listref[2] ) {
          ppt->tag |= MG_REQ;
          nc++;
        }

        /* Check angle */
        if ( fabs(dd) > MMG5_EPSD ) {
          dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
          if ( dd > -mesh->info.dhd ) {
            ppt->tag |= MG_CRN;
            nc++;
          }
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 3 && nc+nre+nm > 0 )
    fprintf(stdout,"     %" MMG5_PRId " corners, %" MMG5_PRId " singular points and %" MMG5_PRId " non manifold points detected\n",nc,nre,nm);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param ref reference of the boundary to analyze (analyze all the boundaries
 * if MMG5_UNSET)
 *
 * \return 1 if success, 0 if fail.
 *
 * Calculate normal vectors at vertices of the mesh.
 *
 */
int MMG2D_norver(MMG5_pMesh mesh, MMG5_int ref) {
  MMG5_pTria       pt,pt1;
  MMG5_pPoint      ppt;
  MMG5_int         k,kk,nn,pleft,pright;
  int8_t           i,ii;

  nn = 0;

  if ( ref == MMG5_UNSET ) {
    /* reset the ppt->s tag */
    for (k=1; k<=mesh->np; ++k) {
      mesh->point[k].s = 0;
    }
  }
  else {
    /** Mark the points that we don't want to analyze */
    /* First: mark all the points */
    for ( k=1; k<=mesh->np; ++k ) {
      mesh->point[k].s = 1;
    }
    /* Second: if the point belongs to a boundary edge or reference ref, remove
     * the mark */
    for ( k=1; k<=mesh->nt; ++k ) {
      pt = &mesh->tria[k];

      for (i=0; i<3; i++) {
        if ( !MG_EDG(pt->tag[i]) ) continue;
        if ( pt->edg[i] != ref )   continue;

        mesh->point[pt->v[MMG5_iprv2[i]]].s = 0;
        mesh->point[pt->v[MMG5_inxt2[i]]].s = 0;
      }
    }
  }

  /* Normal computation */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !MG_EDG(ppt->tag) ) continue;
      if ( ppt->s || (MG_CRN & ppt->tag) || (ppt->tag & MG_NOM) ) continue;

      /* Travel the curve ppt belongs to from left to right until a singularity is met */
      kk = k;
      ii = i;
      do {
        ppt->s = 1;
        if ( !MMG2D_boulen(mesh,kk,ii,&pleft,&pright,ppt->n) ) {
          fprintf(stderr,"\n  ## Error: %s: Impossible to"
                  " calculate normal vector at vertex %" MMG5_PRId ".\n",
                  __func__,MMG2D_indPt(mesh,pt->v[i]));
          return 0;
        }
        nn++;

        kk = pright / 3;
        ii = pright % 3;
        ii = MMG5_iprv2[ii];
        pt1 = &mesh->tria[kk];
        ppt = &mesh->point[pt1->v[ii]];
      }
      while ( !ppt->s && !(MG_CRN  & ppt->tag) && !(ppt->tag & MG_NOM) );

      /* Now travel the curve ppt belongs to from right to left until a singularity is met */
      ppt = &mesh->point[pt->v[i]];
      kk = k;
      ii = i;
      do {
        ppt->s = 1;
        if ( !MMG2D_boulen(mesh,kk,ii,&pleft,&pright,ppt->n) ) {
          fprintf(stderr,"\n  ## Error: %s: Impossible to"
                  " calculate normal vector at vertex %" MMG5_PRId ".\n",
                  __func__,MMG2D_indPt(mesh,pt->v[i]));
          return 0;
        }
        nn++;

        kk = pleft / 3;
        ii = pleft % 3;
        ii = MMG5_inxt2[ii];
        pt1 = &mesh->tria[kk];
        ppt = &mesh->point[pt1->v[ii]];
      }
      while ( !ppt->s && !(MG_CRN & ppt->tag) && !(ppt->tag & MG_NOM) );
    }
  }

  if ( abs(mesh->info.imprim) > 3 && nn > 0 )
    fprintf(stdout,"     %" MMG5_PRId " calculated normal vectors\n",nn);

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 0 if fail, 1 if success
 *
 * Regularize normal vectors at boundary non singular edges with a Laplacian /
 * antilaplacian smoothing
 *
 */
int MMG2D_regnor(MMG5_pMesh mesh) {
  MMG5_pTria            pt;
  MMG5_pPoint           ppt,p1,p2;
  double                *tmp,dd,ps,lm1,lm2,nx,ny,ux,uy,nxt,nyt,res,res0,n[2];
  MMG5_int              k,iel,ip1,ip2,nn,list[MMG5_LMAX];
  int                   it,maxit,ilist;
  int8_t                i,ier;

  it = 0;
  maxit = 10;
  res0 = 0.0;
  nn = 0;
  lm1 = 0.4;
  lm2 = 0.399;

  /* Temporary table for normal vectors */
  MMG5_SAFE_CALLOC(tmp,2*mesh->np+1,double,return 0);

  /* Allocate a seed to each point */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->s = k;
    }
  }

  do {
    /* Step 1: Laplacian */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) continue;
      if ( !MG_EDG(ppt->tag) ) continue;

      nx = ny = 0;
      iel = ppt->s;
      pt  = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k ) i = 1;
      if ( pt->v[2] == k ) i = 2;

      ilist = MMG2D_bouleendp(mesh,iel,i,&ip1,&ip2,list);

      if ( !ilist ) {
        fprintf(stderr,"\n  ## Error: %s: Abort.\n",__func__);
        MMG5_SAFE_FREE(tmp);
        return 0;
      }

      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      /* If p1 is singular, update with the (adequately oriented) normal to the edge ppt-p1 */
      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) {
        ux = p1->c[0] - ppt->c[0];
        uy = p1->c[1] - ppt->c[1];
        dd = ux*ux + uy*uy;
        if ( dd > MMG5_EPSD ) {
          dd = 1.0 / sqrt(dd);
          ux *= dd;
          uy *= dd;
        }
        nxt = -uy;
        nyt = ux;

        ps = nxt*ppt->n[0] + nyt*ppt->n[1];
        if ( ps < 0.0 ) {
          nx += -nxt;
          ny += -nyt;
        }
        else {
          nx += nxt;
          ny += nyt;
        }
      }
      else {
        nx += p1->n[0];
        ny += p1->n[1];
      }

      /* If p2 is singular, update with the (adequately oriented) normal to the edge ppt-p1 */
      if ( MG_SIN(p2->tag) || p2->tag & MG_NOM ) {
        ux = p2->c[0] - ppt->c[0];
        uy = p2->c[1] - ppt->c[1];
        dd = ux*ux + uy*uy;
        if ( dd > MMG5_EPSD ) {
          dd = 1.0 / sqrt(dd);
          ux *= dd;
          uy *= dd;
        }
        nxt = -uy;
        nyt = ux;

        ps = nxt*ppt->n[0] + nyt*ppt->n[1];
        if ( ps < 0.0 ) {
          nx += -nxt;
          ny += -nyt;
        }
        else {
          nx += nxt;
          ny += nyt;
        }
      }
      else {
        nx += p2->n[0];
        ny += p2->n[1];
      }

      /* Average normal normalization */
      dd  = nx*nx + ny*ny;
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nx *= dd;
        ny *= dd;
      }

      /* Laplacian operation */
      tmp[2*(k-1)+1] = ppt->n[0] + lm1 * (nx - ppt->n[0]);
      tmp[2*(k-1)+2] = ppt->n[1] + lm1 * (ny - ppt->n[1]);
    }

    /* Antilaplacian operation */
    res = 0.0;
    for (k=1; k<=mesh->np; k++) {

      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) continue;
      if ( !MG_EDG(ppt->tag) ) continue;

      nx = ny = 0;
      iel = ppt->s;
      pt  = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k ) i = 1;
      if ( pt->v[2] == k ) i = 2;

      ilist = MMG2D_bouleendp(mesh,iel,i,&ip1,&ip2,list);
      if ( !ilist ) {
        fprintf(stderr,"\n  ## Error: %s: Abort.\n",__func__);
        MMG5_SAFE_FREE(tmp);
        return 0;
      }

      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      /* If p1 is singular, update with the (adequately oriented) normal to the edge ppt-p1 */
      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) {
        ux = p1->c[0] - ppt->c[0];
        uy = p1->c[1] - ppt->c[1];
        dd = ux*ux + uy*uy;
        if ( dd > MMG5_EPSD ) {
          dd = 1.0 / sqrt(dd);
          ux *= dd;
          uy *= dd;
        }
        nxt = -uy;
        nyt = ux;

        ps = nxt*ppt->n[0] + nyt*ppt->n[1];
        if ( ps < 0.0 ) {
          nx += -nxt;
          ny += -nyt;
        }
        else {
          nx += nxt;
          ny += nyt;
        }
      }
      else {
        nx += tmp[2*(ip1-1)+1];
        ny += tmp[2*(ip1-1)+2];
      }

      /* If p2 is singular, update with the (adequately oriented) normal to the edge ppt-p1 */
      if ( MG_SIN(p2->tag) || p2->tag & MG_NOM ) {
        ux = p2->c[0] - ppt->c[0];
        uy = p2->c[1] - ppt->c[1];
        dd = ux*ux + uy*uy;
        if ( dd > MMG5_EPSD ) {
          dd = 1.0 / sqrt(dd);
          ux *= dd;
          uy *= dd;
        }
        nxt = -uy;
        nyt = ux;

        ps = nxt*ppt->n[0] + nyt*ppt->n[1];
        if ( ps < 0.0 ) {
          nx += -nxt;
          ny += -nyt;
        }
        else {
          nx += nxt;
          ny += nyt;
        }
      }
      else {
        nx += tmp[2*(ip2-1)+1];
        ny += tmp[2*(ip2-1)+2];
      }

      /* Average normal normalization */
      dd  = nx*nx + ny*ny;
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nx *= dd;
        ny *= dd;
      }

      /* Anti Laplacian operation */
      n[0] = tmp[2*(k-1)+1] - lm2 * (nx - tmp[2*(k-1)+1]);
      n[1] = tmp[2*(k-1)+2] - lm2 * (ny - tmp[2*(k-1)+2]);
      res += (n[0]-ppt->n[0])*(n[0]-ppt->n[0]) + (n[1]-ppt->n[1])*(n[1]-ppt->n[1]);
      ppt->n[0] = n[0];
      ppt->n[1] = n[1];
      nn++;
    }

    /* Normalization */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) continue;
      if ( !MG_EDG(ppt->tag) ) continue;

      dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1];
      if ( dd > MMG5_EPSD ) {
        dd = 1.0 / sqrt(dd);
        ppt->n[0] *= dd;
        ppt->n[1] *= dd;
      }
    }

    if ( it == 0 ) res0 = res;
    if ( res0 > MMG5_EPSD ) res = res / res0;

    if ( mesh->info.imprim < -1 || mesh->info.ddebug ) {
      fprintf(stdout,"     iter %5d  res %.3E",it,res);
      fflush(stdout);
    }
  }
  while ( ++it < maxit && res > MMG5_EPS );


  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     %" MMG5_PRId " normals regularized: %.3e\n",nn,res);

  MMG5_SAFE_FREE(tmp);
  return 1;
}

/**
 * \param mesh pointer towards the mesh
 * \param pt pointer towards current triangle
 * \param k number of current point
 * \param c newly computed coordinates (giving negative area)
 *
 * \return 0 if fail, 1 if success
 *
 * In coordinate regularization, performs a dichotomy between previous point /
 * and newly computed point in the case of negative area
 *
 */
static inline int MMG2D_dichotomy(MMG5_pMesh mesh, MMG5_pTria pt, MMG5_int k, double *c) {

  MMG5_pPoint  ppt;
  double       p[2],co[3][2],vol,to,tp,t;
  int          it,maxit,pos,i,j;

  it = 0;
  maxit = 5;
  to = 0.0;
  tp = 1.0;
  t = 0.5;
  pos = 0;

  ppt = &mesh->point[k];

  for ( i=0 ; i<3 ; i++ ) {
    for ( j=0 ; j<2 ; j++ ) {
      co[i][j] = mesh->point[pt->v[i]].c[j];
    }
  }

  /* initial coordinates of new point */
  p[0] = c[0];
  p[1] = c[1];

  i = 0;
  if ( pt->v[1] == k ) i = 1;
  if ( pt->v[2] == k ) i = 2;

  do {
    co[i][0] = ppt->c[0] + t*(p[0] - ppt->c[0]);
    co[i][1] = ppt->c[1] + t*(p[1] - ppt->c[1]);

    vol = MMG2D_quickarea(co[0],co[1],co[2]);

    if ( vol <= 0.0 ) {
      tp = t;
    }
    else {
      c[0] = co[i][0];
      c[1] = co[i][1];
      to = t;
      pos = 1;
    }

    t = 0.5*(to + tp);
  }
  while ( ++it < maxit );

  if ( pos ) {
    return 1;
  }
  else
    return 0;
}


/**
 * \param mesh pointer toward the mesh
 *
 * \return 0 if fail, 1 if success
 *
 * Regularize vertices coordinates at boundary non singular edges with
 * a Laplacian / antilaplacian smoothing
 *
 */
int MMG2D_regver(MMG5_pMesh mesh) {
  MMG5_pTria            pt;
  MMG5_pPoint           ppt,p1,p2;
  double                *tmp,lm1,lm2,cx,cy,res,res0,c[2],co[3][2],vol;
  MMG5_int              k,kt,iel,ip1,ip2,nn,list[MMG5_LMAX];
  int                   it,maxit,j,ilist,noupdate;
  int8_t                i;

  it = 0;
  maxit=10;
  res0 = 0.0;
  nn = 0;
  lm1 = 0.4;
  lm2 = 0.399;

  /* Temporary table for coordinates */
  MMG5_SAFE_CALLOC(tmp,2*mesh->np+1,double,return 0);

  /* Allocate a seed to each point */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->s = k;
    }
  }

  do {
    /* Step 1: Laplacian */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      tmp[2*(k-1)+1] = ppt->c[0];
      tmp[2*(k-1)+2] = ppt->c[1];

      if ( !MG_VOK(ppt) ) continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) continue;
      if ( !MG_EDG(ppt->tag) ) continue;

      cx = cy = 0;
      iel = ppt->s;
      pt  = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k ) i = 1;
      if ( pt->v[2] == k ) i = 2;

      ilist = MMG2D_bouleendp(mesh,iel,i,&ip1,&ip2,list);

      if ( !ilist ) {
        fprintf(stderr,"\n  ## Error: %s: Abort.\n",__func__);
        MMG5_SAFE_FREE(tmp);
        return 0;
      }

      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      cx += p1->c[0];
      cy += p1->c[1];

      cx += p2->c[0];
      cy += p2->c[1];
      cx *= 0.5;
      cy *= 0.5;

      /* Laplacian operation */
      tmp[2*(k-1)+1] = ppt->c[0] + lm1 * (cx - ppt->c[0]);
      tmp[2*(k-1)+2] = ppt->c[1] + lm1 * (cy - ppt->c[1]);
    }

    /* Antilaplacian operation */
    res = 0.0;
    for (k=1; k<=mesh->np; k++) {

      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) continue;
      if ( !MG_EDG(ppt->tag) ) continue;

      cx = cy = 0;
      iel = ppt->s;
      pt  = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k ) i = 1;
      if ( pt->v[2] == k ) i = 2;

      ilist = MMG2D_bouleendp(mesh,iel,i,&ip1,&ip2,list);

      if ( !ilist ) {
        fprintf(stderr,"\n  ## Error: %s: Abort.\n",__func__);
        MMG5_SAFE_FREE(tmp);
        return 0;
      }

      cx += tmp[2*(ip1-1)+1];
      cy += tmp[2*(ip1-1)+2];

      cx += tmp[2*(ip2-1)+1];
      cy += tmp[2*(ip2-1)+2];
      cx *= 0.5;
      cy *= 0.5;

      /* Anti Laplacian operation */
      c[0] = tmp[2*(k-1)+1] - lm2 * (cx - tmp[2*(k-1)+1]);
      c[1] = tmp[2*(k-1)+2] - lm2 * (cy - tmp[2*(k-1)+2]);

      /* Check if updated point creates a triangle with negative area.
         If it does, performs a dichotomy to find optimal point */
      noupdate = 0;
      for (kt=0; kt<ilist;kt++) {
        pt = &mesh->tria[list[kt]];
        if ( !MG_EOK(pt) ) continue;

        for ( i=0 ; i<3 ; i++ ) {
          for ( j=0 ; j<2 ; j++ ) {
            co[i][j] = mesh->point[pt->v[i]].c[j];
          }
        }

        i = 0;
        if ( pt->v[1] == k ) i = 1;
        if ( pt->v[2] == k ) i = 2;

        for ( j=0 ; j<2 ; j++ ) {
          co[i][j] = c[j];
        }

        vol = MMG2D_quickarea(co[0],co[1],co[2]);

        if ( vol <= 0.0 ) {
          if (!MMG2D_dichotomy(mesh,pt,k,c)) noupdate = 1;
          continue;
        }
      }
      if ( !noupdate ) {
        res += (c[0]-ppt->c[0])*(c[0]-ppt->c[0]) + (c[1]-ppt->c[1])*(c[1]-ppt->c[1]);
        ppt->c[0] = c[0];
        ppt->c[1] = c[1];
        nn++;
      }
    }

    if ( it == 0 ) res0 = res;
    if ( res0 > MMG5_EPSD ) res = res / res0;

    if ( mesh->info.imprim < -1 || mesh->info.ddebug ) {
      fprintf(stdout,"     iter %5d  res %.3E",it,res);
      fflush(stdout);
    }
  }
  while ( ++it < maxit && res > MMG5_EPS );

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     %" MMG5_PRId " coordinates regularized: %.3e\n",nn,res);

  MMG5_SAFE_FREE(tmp);
  return 1;
}

/** preprocessing stage: mesh analysis */
int MMG2D_analys(MMG5_pMesh mesh) {

  /* Transfer the boundary edge references to the triangles, if it has not been
   * already done (option 1) */
  if ( !MMG2D_assignEdge(mesh) ) {
     fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* Creation of tria adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
     fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Creation quadrilaterals adjacency relations in the mesh */
  if ( !MMG2D_hashQuad(mesh) ) {
     fprintf(stderr,"\n  ## Quadrilaterals hashing problem. Exit program.\n");
    return 0;
  }

  /* Set tags to triangles from geometric configuration */
  int8_t init_cc  = 0;
  if ( mesh->info.isosurf ) {
    init_cc = 1;
  }

  /** \warning possible time improvment here: in lssurf mode, the second call of
   * #MMG2D_setadj only adds BDY tag to points that have been inserted by the
   * level-set splitting (the rest of the job of setadj has alredy be done by
   * the first call).  */
  if ( !MMG2D_setadj(mesh,init_cc) ) {
    fprintf(stderr,"\n  ## Problem in function setadj. Exit program.\n");
    return 0;
  }

  /* Identify singularities in the mesh */
  if ( !MMG2D_singul(mesh,MMG5_UNSET) ) {
     fprintf(stderr,"\n  ## Problem in identifying singularities. Exit program.\n");
    return 0;
  }

  /* Regularize vertix vector field with a Laplacian / anti-laplacian smoothing */
  if ( mesh->info.xreg && !MMG2D_regver(mesh) ) {
    fprintf(stderr,"\n  ## Problem in regularizing vertices coordinates. Exit program.\n");
    return 0;
  }

  /* Define normal vectors at vertices on curves */
  if ( !MMG2D_norver(mesh,MMG5_UNSET) ) {
     fprintf(stderr,"\n  ## Problem in calculating normal vectors. Exit program.\n");
    return 0;
  }

  /* Regularize normal vector field with a Laplacian / anti-laplacian smoothing */
  if ( mesh->info.nreg && !MMG2D_regnor(mesh) ) {
      fprintf(stderr,"\n  ## Problem in regularizing normal vectors. Exit program.\n");
      return 0;
  }

  if ( mesh->nquad ) MMG5_DEL_MEM(mesh,mesh->adjq);

  return 1;
}
