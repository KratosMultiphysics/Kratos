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
 * \file mmg2d/mmg2d1.c
 * \brief Mesh adaptation functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"

/* Mesh adaptation routine for the first stages of the algorithm: intertwine splitting
 based on patterns, collapses and swaps.
   typchk = 1 -> adaptation based on edge lengths
   typchk = 2 -> adaptation based on lengths calculated in metric met */
int MMG2D_anatri(MMG5_pMesh mesh,MMG5_pSol met,int8_t typchk) {
  int       it,maxit;
  MMG5_int  ns,nc,nsw,nns,nnc,nnsw;

  nns = nnc = nnsw = 0;
  it = 0;
  maxit = 5;

  /* Main routine; intertwine split, collapse and swaps */
  do {
    if ( typchk == 2 && it == 0 ) {
// #warning Luca: check consistency with 3D
    }
    if ( !mesh->info.noinsert ) {
      /* Memory free */
      MMG5_DEL_MEM(mesh,mesh->adja);
      mesh->adja = 0;

      /* Split long edges according to patterns */
      ns = MMG2D_anaelt(mesh,met,typchk);
      if ( ns < 0 ) {
        fprintf(stderr,"  ## Unable to complete surface mesh. Exit program.\n");
        return 0;
      }

      /* Recreate adjacencies */
      if ( !MMG2D_hashTria(mesh) ) {
        fprintf(stdout,"  ## Hashing problem. Exit program.\n");
        return 0;
      }

      /* Collapse short edges */
      nc = MMG2D_colelt(mesh,met,typchk);
      if ( nc < 0 ) {
        fprintf(stderr,"  ## Unable to collapse mesh. Exiting.\n");
        return 0;
      }

    }
    else {
      ns = 0;
      nc = 0;
    }
    /* Swap edges */
    if ( !mesh->info.noswap ) {
      nsw = MMG2D_swpmsh(mesh,met,typchk);
      if ( nsw < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
    }
    else nsw = 0;

    nns += ns;
    nnc += nc;
    nnsw += nsw;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped\n",ns,nc,nsw);
    if ( it > 3 && MMG5_abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nsw >0 );

  if ( mesh->info.imprim > 0 ) {
    if ( (abs(mesh->info.imprim) < 5 || mesh->info.ddebug ) && nns+nnc > 0 )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped, %d iter.\n",nns,nnc,nnsw,it);
  }
  return 1;
}

/* Travel triangles and split long edges according to patterns */
MMG5_int MMG2D_anaelt(MMG5_pMesh mesh,MMG5_pSol met,int typchk) {
  MMG5_pTria      pt;
  MMG5_pPoint     ppt,p1,p2;
  MMG5_Hash       hash;
  double          len,s,o[2],no[2];
  int             it;
  MMG5_int        ni,npinit,ns,nc,k,nt,ip1,ip2,ip,vx[3];
  int8_t          i,ic,i1,i2,ier;
  static int8_t   mmgWarn0=0;

  s = 0.5;
  ns = 0;
  npinit = mesh->np;

  if ( !MMG5_hashNew(mesh,&hash,mesh->np,3*mesh->np) ) return 0;

  /* Step 1: travel mesh, check edges, and tag those to be split; create the new vertices in hash */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || (pt->ref < 0) ) continue;
    if ( MG_SIN(pt->tag[0]) || MG_SIN(pt->tag[1]) || MG_SIN(pt->tag[2]) )  continue;

    /* Check if pt should be cut */
    pt->flag = 0;

    /* typchk=1 : split on geometric basis and rough size considerations */
    if ( typchk == 1) {
      if ( !MMG2D_chkedg(mesh,k) ) continue;
    }
    /* typchk =2 : split on basis of edge lengths in the metric */
    else if ( typchk ==2 ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);
        if ( len > MMG2D_LLONG ) MG_SET(pt->flag,i);
      }
    }

    /* mesh->info.fem : split edges which are not MG_BDY, but whose vertices are both MG_BDY */
    if ( mesh->info.fem ) {
      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];
        if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) && !(pt->tag[i] & MG_BDY) ) MG_SET(pt->flag,i);
      }
    }

    if ( !pt->flag ) continue;
    ns++;

    /* Travel edges to split and create new points */
    for (i=0; i<3; i++) {
      if ( !MG_GET(pt->flag,i) ) continue;
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      ip = MMG5_hashGet(&hash,ip1,ip2);
      if ( ip > 0 ) continue;

      /* Geometric attributes of the new point */
      ier = MMG2D_bezierCurv(mesh,k,i,s,o,no);
      if ( !ier ) {
        MG_CLR(pt->flag,i);
        continue;
      }
      ip = MMG2D_newPt(mesh,o,pt->tag[i]);
      if ( !ip ) {
        /* reallocation of point table */
        MMG2D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            do {
                              MMG2D_delPt(mesh,mesh->np);
                            } while ( mesh->np>npinit );return -1;,
                            o,pt->tag[i]);
      }
      ppt = &mesh->point[ip];
      if ( MG_EDG(pt->tag[i]) ) {
        ppt->n[0] = no[0];
        ppt->n[1] = no[1];
      }

      /* If there is a metric in the mesh, interpolate it at the new point */
      assert ( met );
      if ( met->m )
        MMG2D_intmet(mesh,met,k,i,ip,s);

      /* Add point to the hashing structure */
      MMG5_hashEdge(mesh,&hash,ip1,ip2,ip);
    }
  }
  if ( !ns ) {
    MMG5_DEL_MEM(mesh,hash.item);
    return ns;
  }

  /* Step 2: Make flags at triangles consistent between themselves (check if adjacent triangle is split) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 ) continue;
    else if ( pt->flag == 7 ) continue;
    nc = 0;

    for (i=0; i<3; i++) {
      i1 = MMG5_iprv2[i];
      i2 = MMG5_inxt2[i];
      if ( !MG_GET(pt->flag,i) && !MG_SIN(pt->tag[i]) ) {
        ip = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( ip > 0 ) {
          MG_SET(pt->flag,i);
          nc++;
        }
      }
    }
    if ( nc > 0 ) ns++;
  }
  if ( mesh->info.ddebug && ns ) {
    fprintf(stdout,"     %" MMG5_PRId " analyzed  %" MMG5_PRId " proposed\n",mesh->nt,ns);
    fflush(stdout);
  }
  /* Step 3: Simulate splitting and delete points leading to an invalid configuration */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  it = 1;
  nc = 0;
  do {
    ni = 0;
    for ( k=1; k<= mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) || pt->ref < 0 ) continue;
      else if ( pt->flag == 0 ) continue;

      vx[0] = vx[1] =vx[2] = 0;
      pt->flag = 0;
      ic = 0;

      for (i=0; i<3; i++) {
        i1 = MMG5_iprv2[i];
        i2 = MMG5_inxt2[i];
        vx[i] = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( vx[i] > 0 ) {
          MG_SET(pt->flag,i);
          if ( mesh->point[vx[i]].flag > 2 )  ic = 1;
        }
      }

      if ( !pt->flag )  continue;
      switch (pt->flag) {
        case 1: case 2: case 4:
          ier = MMG2D_split1_sim(mesh,met,k,vx);
          break;
        case 7:
          ier = MMG2D_split3_sim(mesh,met,k,vx);
          break;
        default:
          ier = MMG2D_split2_sim(mesh,met,k,vx);
          break;
      }
      if ( ier )  continue;

      /* An edge is invalidated in the process */
      ni++;
      if ( ic == 0 && MMG2D_dichoto(mesh,met,k,vx) ) {
        for (i=0; i<3; i++)
          if ( vx[i] > 0 )  mesh->point[vx[i]].flag++;
      }
      /* Relocate point at the center of the edge */
      else {
        for (i=0; i<3; i++) {
          if ( vx[i] > 0 ) {
            p1 = &mesh->point[pt->v[MMG5_iprv2[i]]];
            p2 = &mesh->point[pt->v[MMG5_inxt2[i]]];
            ppt = &mesh->point[vx[i]];
            ppt->c[0] = 0.5 * (p1->c[0] + p2->c[0]);
            ppt->c[1] = 0.5 * (p1->c[1] + p2->c[1]);
          }
        }
      }
    }
    nc += ni;
  }
  while ( ni >0 && ++it <20 );

  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %" MMG5_PRId " corrected,  %" MMG5_PRId " invalid\n",nc,ni);
    fflush(stdout);
  }

  /* step 4: effective splitting */
  ns = 0;
  nt = mesh->nt;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )  continue;
    else if ( pt->flag == 0 )  continue;

    vx[0] = vx[1] = vx[2] = 0;
    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      if ( MG_GET(pt->flag,i) ) {
        vx[i] = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( !vx[i] ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: unable to create point on"
                    " at least 1 edge.\n Exit program.\n",__func__);
          }
          return -1;
        }
      }
    }
    if ( pt->flag == 1 || pt->flag == 2 || pt->flag == 4 ) {
      ier = MMG2D_split1(mesh,met,k,vx);
      ns++;
    }
    else if ( pt->flag == 7 ) {
      ier = MMG2D_split3(mesh,met,k,vx);
      ns++;
    }
    else {
      ier = MMG2D_split2(mesh,met,k,vx);
      ns++;
    }
    if ( !ier ) return -1;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);
  MMG5_DEL_MEM(mesh,hash.item);

  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param vx pointer toward table of edges to split.
 * \return 1.
 *
 * Find acceptable position for splitting.
 *
 */
int MMG2D_dichoto(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,MMG5_int *vx) {
  MMG5_pTria   pt;
  MMG5_pPoint  pa,pb,ps;
  double       o[3][2],p[3][2];
  float        to,tp,t;
  MMG5_int     ia,ib;
  int          ier,it,maxit;
  int8_t       i,i1,i2;

  pt = &mesh->tria[k];

  /* Get point on curve and along segment for edge split */
  for (i=0; i<3; i++) {
    memset(p[i],0,2*sizeof(double));
    memset(o[i],0,2*sizeof(double));
    if ( vx[i] > 0 ) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      ia = pt->v[i1];
      ib = pt->v[i2];
      pa = &mesh->point[ia];
      pb = &mesh->point[ib];
      ps = &mesh->point[vx[i]];
      o[i][0] = 0.5 * (pa->c[0] + pb->c[0]);
      o[i][1] = 0.5 * (pa->c[1] + pb->c[1]);
      p[i][0] = ps->c[0];
      p[i][1] = ps->c[1];
    }
  }
  maxit = 4;
  it = 0;
  tp = 1.0;
  to = 0.0;

  do {
    /* compute new position */
    t = 0.5 * (tp + to);
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
      }
    }
    switch (pt->flag) {
      case 1: case 2: case 4:
        ier = MMG2D_split1_sim(mesh,met,k,vx);
        break;
      case 7:
        ier = MMG2D_split3_sim(mesh,met,k,vx);
        break;
      default:
        ier = MMG2D_split2_sim(mesh,met,k,vx);
        break;
    }
    if ( ier )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );

  /* restore coords of last valid pos. */
  if ( !ier ) {
    t = to;
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
      }
    }
  }
  return 1;
}

/* Travel triangles and collapse short edges */
MMG5_int MMG2D_colelt(MMG5_pMesh mesh,MMG5_pSol met,int typchk) {
  MMG5_pTria   pt;
  MMG5_pPoint  p1,p2;
  double       ux,uy,ll,hmin2;
  MMG5_int     k;
  int          ilist;
  MMG5_int     nc;
  uint8_t      i,i1,i2,open;
  MMG5_int     list[MMG5_TRIA_LMAX+2];

  nc = 0;
  hmin2 = mesh->info.hmin * mesh->info.hmin;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 ) continue;

    /* Travel 3 edges of the triangle and decide whether to collapse p1->p2, based on length criterion */
    pt->flag = 0; // was here before, but probably serves for nothing

    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) continue;

      /* Impossible to collapse a surface point onto a non surface point -- impossible to
       collapse a surface point along a non geometric edge */
      else if ( p1->tag & MG_GEO ) {
        if ( ! (p2->tag & MG_GEO) || !(pt->tag[i] & MG_GEO) ) continue;
      }
      /* Same test for REF points */
      else if ( p1->tag & MG_REF ) {
        if ( ! (p2->tag & MG_GEO || p2->tag & MG_REF) || !(pt->tag[i] & MG_REF) ) continue;
      }

      open = (mesh->adja[3*(k-1)+1+i] == 0) ? 1 : 0;

      /* Check length */
      if ( typchk == 1 ) {
        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];
        ll = ux*ux + uy*uy;
        if ( ll > hmin2 ) continue;
      }
      else {
        ll = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);
        if ( ll > MMG2D_LSHRT ) continue;
      }

      /* Check whether the geometry is preserved */
      ilist = MMG2D_chkcol(mesh,met,k,i,list,typchk);

      if ( ilist > 3 || ( ilist == 3 && open ) ) {
        nc += MMG2D_colver(mesh,ilist,list);
        break;
      }
      else if ( ilist == 3 ) {
        nc += MMG2D_colver3(mesh,list);
        break;
      }
      else if ( ilist == 2 ) {
        nc += MMG2D_colver2(mesh,list);
        break;
      }
    }
  }

  if ( nc > 0 && (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) )
    fprintf(stdout,"     %8" MMG5_PRId " vertices removed\n",nc);

  return nc;
}

/* Travel triangles and swap edges to improve quality */
MMG5_int MMG2D_swpmsh(MMG5_pMesh mesh,MMG5_pSol met,int typchk) {
  MMG5_pTria pt;
  int        it,maxit;
  MMG5_int   ns,nns;
  MMG5_int   k;
  uint8_t    i;

  it = nns = 0;
  maxit = 2;
  mesh->base++;

  do {
    ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) || pt->ref < 0 ) continue;

      for (i=0; i<3; i++) {
        if ( MG_SIN(pt->tag[i]) || MG_EDG(pt->tag[i]) ) continue;
        else if ( MMG2D_chkswp(mesh,met,k,i,typchk) ) {
          ns += MMG2D_swapar(mesh,k,i);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ns > 0 && ++it<maxit );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8" MMG5_PRId " edge swapped\n",nns);
  return nns;
}


/* Mesh adaptation routine for the final stage of the algorithm: intertwine splitting
 based on patterns, collapses, swaps and vertex relocations.*/
int MMG2D_adptri(MMG5_pMesh mesh,MMG5_pSol met) {
  int                  maxit,it;
  MMG5_int             nns,ns,nnc,nc,nnsw,nsw,nnm,nm;

  nns = nnc = nnsw = nnm = it = 0;
  maxit = 5;

  do {

    if ( !mesh->info.noinsert ) {
      ns = MMG2D_adpspl(mesh,met);
      if ( ns < 0 ) {
        fprintf(stderr,"  ## Problem in function adpspl."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }

      nc = MMG2D_adpcol(mesh,met);
      if ( nc < 0 ) {
        fprintf(stderr,"  ## Problem in function adpcol."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else {
      ns = 0;
      nc = 0;
    }
    
    if ( !mesh->info.noswap ) {
      nsw = MMG2D_swpmsh(mesh,met,2);
      if ( nsw < 0 ) {
        fprintf(stderr,"  ## Problem in function swpmsh."
                " Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else
      nsw = 0;

    if ( !mesh->info.nomove ) {
      nm = MMG2D_movtri(mesh,met,1,0);
      if ( nm < 0 ) {
        fprintf(stderr,"  ## Problem in function movtri. "
                "Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else
      nm = 0;

    nns  += ns;
    nnc  += nc;
    nnsw += nsw;
    nnm  += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nsw+nm > 0 )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",ns,nc,nsw,nm);
    if ( ns < 10 && MMG5_abs(nc-ns) < 3 )  break;
    else if ( it > 3 && MMG5_abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && (nc+ns+nsw+nm > 0) );

  /* Last iterations of vertex relocation only */
  if ( !mesh->info.nomove ) {
    nm = MMG2D_movtri(mesh,met,5,1);
    if ( nm < 0 ) {
      fprintf(stderr,"  ## Problem in function movtri. Unable to complete mesh."
              " Exit program.\n");
      return 0;
    }
    nnm += nm;
  }
  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnc > 0 || nns > 0) )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter. \n",nns,nnc,nnsw,nnm,it);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * \return -1 if failed or number of new points.
 *
 * Analysis and splitting routine for edges in the final step of the algorithm;
 * edges are only splitted on a one-by-one basis
 *
 */
MMG5_int MMG2D_adpspl(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria         pt;
  double             lmax,len;
  MMG5_int           ip,ns,k,nt;
  int                ier;
  int8_t             i,i1,i2,imax;

  ns = 0;

  /*loop until nt to avoid the split of new triangle*/
  nt = mesh->nt;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 ) continue;

    imax = -1;
    lmax = -1.0;
    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }

    if ( lmax < MMG2D_LOPTL ) continue;
    else if ( MG_SIN(pt->tag[imax]) ) continue;

    /* Check the feasibility of splitting */
    ip = MMG2D_chkspl(mesh,met,k,imax);

    /* Lack of memory; abort the routine */
    if ( ip < 0 ){
      return ns;
    }
    else if ( ip > 0 ) {
      ier = MMG2D_split1b(mesh,k,imax,ip);

      /* Lack of memory; abort the routine */
      if ( !ier ) {
        MMG2D_delPt(mesh,ip);
        return ns;
      }
      ns += ier;
    }
  }

  return ns;
}

/* Analysis and collapse routine for edges in the final step of the algorithm */
int MMG2D_adpcol(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria        pt;
  MMG5_pPoint       p1,p2;
  double            len;
  MMG5_int          k,nc;
  int               ilist;
  int8_t            i,i1,i2,open;
  MMG5_int          list[MMG5_TRIA_LMAX+2];
  
  nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 ) continue;

    /* Check edge length, and attempt collapse */
    pt->flag = 0;
    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;

      open = ( mesh->adja[3*(k-1)+1+i] == 0 ) ? 1 : 0;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];

      if ( MG_SIN(p1->tag) || p1->tag & MG_NOM ) continue;
      else if ( p1->tag & MG_GEO ) {
        if ( ! (p2->tag & MG_GEO) || !(pt->tag[i] & MG_GEO) ) continue;
      }
      else if ( p1->tag & MG_REF ) {
        if ( ! (p2->tag & MG_GEO || p2->tag & MG_REF) || !(pt->tag[i] & MG_REF) ) continue;
      }

      len = MMG2D_lencurv(mesh,met,pt->v[i1],pt->v[i2]);

      if ( len > MMG2D_LOPTS ) continue;

      ilist = MMG2D_chkcol(mesh,met,k,i,list,2);

      if ( ilist > 3 || ( ilist==3 && open ) ) {
        nc += MMG2D_colver(mesh,ilist,list);
        break;
      }
      else if ( ilist == 3 ) {
        nc += MMG2D_colver3(mesh,list);
        break;
      }
      else if ( ilist == 2 ) {
        nc += MMG2D_colver2(mesh,list);
        break;
      }
    }
  }

  return nc;
}

/* Analyze points to relocate them according to a quality criterion */
MMG5_int MMG2D_movtri(MMG5_pMesh mesh,MMG5_pSol met,int maxit,int8_t improve) {
  MMG5_pTria           pt;
  MMG5_pPoint          p0;
  MMG5_int             nnm,nm,ns,k;
  int                  it,ilist;
  int8_t               i,ier;
  MMG5_int             base,list[MMG5_TRIA_LMAX+2];

  it = nnm = 0;
  base = 0;

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) || pt->ref < 0 ) continue;

      for (i=0; i<3; i++) {
        p0 = &mesh->point[pt->v[i]];
        if ( p0->flag == base || MG_SIN(p0->tag) || p0->tag & MG_NOM ) continue;

        int8_t dummy;
        ilist = MMG5_boulet(mesh,k,i,list,0,&dummy);

        if ( MG_EDG(p0->tag) ) {
          ier = MMG2D_movedgpt(mesh,met,ilist,list,improve);
          if ( ier ) ns++;
        }
        else {
          if ( met->size == 3 && met->m )
            ier = MMG2D_movintpt_ani(mesh,met,ilist,list,improve);
          else
            ier = MMG2D_movintpt(mesh,met,ilist,list,improve);
        }

        if ( ier ) {
          nm++;
          p0->flag = base;
        }
      }
    }
    nnm += nm;
    if ( mesh->info.ddebug )  fprintf(stdout,"     %8" MMG5_PRId " moved, %" MMG5_PRId " geometry\n",nm,ns);
  }
  while ( ++it < maxit && nm > 0 );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm > 0 )
    fprintf(stdout,"     %8" MMG5_PRId " vertices moved, %d iter.\n",nnm,it);

  return nnm;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 1 if success, 0 if strongly fail.
 *
 * Mesh adaptation -- new version of mmg2d1.c
 *
 **/
int MMG2D_mmg2d1n(MMG5_pMesh mesh,MMG5_pSol met) {

  /* Stage 1: creation of a geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !MMG2D_anatri(mesh,met,1) ) {
    fprintf(stderr,"  ## Unable to split mesh-> Exiting.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_ANATRI1 to save adapted mesh at the end of
   * anatri wave */
  if ( getenv("MMG_SAVE_ANATRI1") ) {
    printf("  ## WARNING: EXIT AFTER ANATRI-1."
           " (MMG_SAVE_ANATRI1 env variable is exported).\n");
    return 1;
  }

  /* Stage 2: creation of a computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ){
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");
}
  if ( !MMG2D_defsiz(mesh,met) ) {
    fprintf(stderr,"  ## Metric undefined. Exit program.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_DEFSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_DEFSIZ") ) {
    printf("  ## WARNING: EXIT AFTER DEFSIZ."
           " (MMG_SAVE_DEFSIZ env variable is exported).\n");
    return 1;
  }

  MMG5_gradation_info(mesh);
  if ( mesh->info.hgrad > 0. ) {
    if (!MMG2D_gradsiz(mesh,met) ) {
      fprintf(stderr,"  ## Gradation problem. Exit program.\n");
      return 0;
    }
  }

  if ( mesh->info.hgradreq > 0. ) {
    MMG2D_gradsizreq(mesh,met);
  }
  /* Debug: export variable MMG_SAVE_GRADSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_GRADSIZ") ) {
    printf("  ## WARNING: EXIT AFTER GRADSIZ."
           " (MMG_SAVE_GRADSIZ env variable is exported).\n");
    return 1;
  }

  if ( !MMG2D_anatri(mesh,met,2) ) {
    fprintf(stderr,"  ## Unable to proceed adaptation. Exit program.\n");
    return 0;
  }
  /* Debug: export variable MMG_SAVE_ANATRI1 to save adapted mesh at the end of
   * anatri wave */
  if ( getenv("MMG_SAVE_ANATRI1") ) {
    printf("  ## WARNING: EXIT AFTER ANATRI-2."
           " (MMG_SAVE_ANATRI2 env variable is exported).\n");
    return 1;
  }


  /* Stage 3: fine mesh improvements */
  if ( !MMG2D_adptri(mesh,met) ) {
    fprintf(stderr,"  ## Unable to make fine improvements. Exit program.\n");
    return 0;
  }

  return 1;
}
