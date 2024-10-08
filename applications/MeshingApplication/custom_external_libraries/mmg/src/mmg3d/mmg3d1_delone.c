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
 * \file mmg3d/mmg3d1_delone.c
 * \brief Perform volume and surface mesh adaptation in delaunay mode.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 * Perform volume and surface mesh adaptation in delaunay mode (\a
 * MMG_PATTERN preprocessor flag set to OFF).
 *
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "mmg3dexterns_private.h"

#ifndef MMG_PATTERN

int8_t  ddb;

#define MMG3D_THRES_DEL     1.6
#define MMG3D_LOPTL_DEL     1.41
#define MMG3D_LFILTS_DEL    0.7
#define MMG3D_LFILTL_DEL    0.2

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k index of tetra in which we work.
 * \param imax index in \a k of edge that we consider for split.
 * \param lmax length of edge \a imax.
 * \param lmaxtet length of largest edge of tetra \a k.
 * \param 1 if we want to check tetra with 4 ridge metrics.
 * \param ifilt pointer to store the number of vertices filtered by the PROctree.
 * \param ns pointer toward count of splits (has to be updated)
 * \param warn pointer to store a flag that warn the user in case of
 * reallocation error.
 * \param countMemFailure number of memory errors (to update)
 *
 * \return -2 for low failure (mesh has to be saved).
 *
 * \return -1 for strong failure.
 *
 * \return 0 if edge cannot be splitted and if we want to pass to next loop step
 * (next element or next tetra edge)
 *
 * \return 1 if edge cannot be splitted and we want to try to collapse too long
 * edge.
 *
 * \return 2 if edge has been splitted and we want to treat next element.
 *
 * \return 3 if nothing has been done (no error but no split either).
 *
 * Try to split \a imax if too large.
 *
 */
static inline
int MMG3D_mmg3d1_delone_split(MMG5_pMesh mesh, MMG5_pSol met,
                              MMG3D_pPROctree *PROctree,MMG5_int k,
                              int8_t imax,double lmax,double lmaxtet,
                              int8_t chkRidTet,MMG5_int *ifilt,MMG5_int *ns,
                              int *warn,int8_t *countMemFailure ) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  double        o[3];
  int64_t       list[MMG3D_LMAX+2];
  MMG5_int      ip1,ip2;
  MMG5_int      src,ip;
  int           ilist;
  int8_t        j,i,i1,i2;

  /** Get edge infos */
  pt = &mesh->tetra[k];
  pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

  MMG3D_find_bdyface_from_edge(mesh,pt,imax,&i,&j,&i1,&i2,&ip1,&ip2,&p0,&p1);

  if ( pxt && (pxt->ftag[i] & MG_BDY) ) {

    /** Check edge length */
    if ( lmax < MMG3D_LOPTL_DEL )  {
      /* Edge is small enough: nothing to do */
      return 3;
    }

    /** Edge belongs to a boundary face: try to split using patterns */
    /* Construction of bezier edge */
    double      to[3],no1[3],no2[3];
    MMG5_int    ref;
    int16_t     tag;
    MMG5_pPoint ppt;

    int8_t ier = MMG3D_build_bezierEdge(mesh,k,imax,i,j,pxt,ip1,ip2,p0,p1,
                                        &ref,&tag,o,to,no1,no2,list,&ilist);
    switch (ier) {
    case -1:
      /* Strong failure */
    case 0:
      /* Unable to split edge: pass to next elt */
    case 1:
      /* Unable to split edge: try to collapse shortest edge */

      /* For all previous cases, return ier */
      return ier;
    }
    assert ( ier==2 && "unexpected return value for build_bezierEdge");

    /** b/ Edge splitting */
#ifdef USE_POINTMAP
    src = mesh->point[ip1].src;
#else
    src = 1;
#endif
    ip = MMG3D_newPt(mesh,o,tag,src);
    if ( !ip ){
      /* reallocation of point table */
      MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                          *warn=1;++(*countMemFailure);
                          return 1,
                          o,tag,src);
    }

    ier = 1;
    if ( met && met->m ) {
      ier = MMG5_intmet(mesh,met,k,imax,ip,0.5);
    }
    if ( ier<=0 ) {
      MMG3D_delPt(mesh,ip);
      return 1;
    }

    /* Simulation only if intmet call is successful */
    /* simbulgept needs a valid tangent at ridge point (to build ridge metric in
     * order to comute edge lengths). Thus we need to store the geometric info of
     * point here. */
    ppt = &mesh->point[ip];
    MMG3D_set_geom(mesh,ppt,tag,ref,pxt->ref[i],no1,no2,to);

    ier = MMG3D_simbulgept(mesh,met,list,ilist,ip);

    if ( ier == 2 || ier < 0 ) {
      /* int met failure or sharp angle failure */
      MMG3D_delPt(mesh,ip);
      return 1;
    }
    if ( ier == 0 ) {
      /* very bad quality failure */
      ier = MMG3D_dichoto1b(mesh,met,list,ilist,ip);
    }
    if ( ier == 1 ) {
      ier = MMG5_split1b(mesh,met,list,ilist,ip,1,1,chkRidTet);
    }

    /* if we realloc memory in MMG5_split1b pt and pxt pointers are not valid */
    if ( ier < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to split.\n",__func__);
      MMG3D_delPt(mesh,ip);
      return -1;
    }
    if ( ier == 0 || ier == 2 ) {
      MMG3D_delPt(mesh,ip);
      return 1;
    }

    (*ns)++;

    return 2;
    /* End of case of a bdy face */
  }
  else {
    /** Case of a tetra without xtetra (no boundary faces): split non-bdy
     * edges with Delauney kernel. */
    /* Note that it is possible that non bdy tetra contains a bdy edge, here
     * only non bdy edge are considered */

    int8_t force_splt = 0;
    const int8_t fem_mode = 2; // value of info.fem in case of fem mode

    if ( mesh->info.fem == fem_mode ) {
      /* Force splitting of internal edges connecting bdy points */
      if ( MG_TRUE_BDY(p0->tag) && MG_TRUE_BDY(p1->tag) ) {
        force_splt = 1;
      }
    }

    if ( (!force_splt) && lmax < MMG3D_LOPTL_DEL )  {
      /* Edge is small enough: nothing to do */
      return 3;
    }

    int8_t isbdy;
    ilist = MMG5_coquil(mesh,k,imax,list,&isbdy);

    if ( !ilist ){
      /* Unable to compute edge shell: treat next element */
      return 0;
    }
    else if ( isbdy ) {
      /* Edge is bdy: skip it (we want to treat it from a bdy tetra) */
      return 0;
    }
    else if ( ilist<0 ) {
      return -1;
    }

    o[0] = 0.5*(p0->c[0] + p1->c[0]);
    o[1] = 0.5*(p0->c[1] + p1->c[1]);
    o[2] = 0.5*(p0->c[2] + p1->c[2]);
#ifdef USE_POINTMAP
    src = mesh->point[ip1].src;
#else
    src = 1;
#endif
    ip = MMG3D_newPt(mesh,o,MG_NOTAG,src);

    if ( !ip )  {
      /* reallocation of point table */
      MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                          *warn=1;++(*countMemFailure);
                          return 1,
                          o,MG_NOTAG,src);
    }

    int ier = 1;
    if ( met && met->m ) {
      ier = MMG5_intmet(mesh,met,k,imax,ip,0.5);
    }
    if ( ier<=0 ) {
      MMG3D_delPt(mesh,ip);
      return 1;
    }

    /* Delaunay */
    double lfilt;
    if ( lmaxtet< MMG3D_THRES_DEL ) {
      lfilt = MMG3D_LFILTS_DEL;
    }
    else {
      lfilt = MMG3D_LFILTL_DEL;
    }

    /* No filter for internal edges connecting boundary points: we want to force
     * splitting */
    if ( force_splt ) {
      lfilt = 0;
    }

    ier = 1;
    if ( *PROctree ) {
      ier = MMG3D_PROctreein(mesh,met,*PROctree,ip,lfilt);
    }

    if ( ier == 0 ) {
      /* PROctree allocated and PROctreein refuse the insertion */
      MMG3D_delPt(mesh,ip);
      (*ifilt)++;
      return 1;
    }

    if ( ier < 0 ) {
      /* PROctree allocated but PROctreein fail due to lack of memory */
      MMG3D_freePROctree ( mesh,PROctree );
      MMG3D_delPt(mesh,ip);
      (*ifilt)++;
      return 1;
    }

    int lon = MMG5_cavity(mesh,met,k,ip,list,ilist/2,MMG5_EPSOK);
    if ( lon < 1 ) {
      MMG3D_delPt(mesh,ip);
      return 1;
    }

    int ret = MMG5_delone(mesh,met,ip,list,lon);
    if ( ret > 0 ) {
      if ( *PROctree ) {
        MMG3D_addPROctree(mesh,*PROctree,ip);
      }
      (*ns)++;
      return 2;
    }
    if ( ret == 0 ) {
      MMG3D_delPt(mesh,ip);
      return 1;
    }

    /* Allocation problem ==> savemesh */
    MMG3D_delPt(mesh,ip);
    return -2;
  }

  return 3;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k index of tetra in which we work.
 * \param imin index in \a k of edge that we consider for split.
 * \param lmin length of edge \a imax.
 * \param imax index in \a k of edge that we consider for split.
 * \param lmax length of edge \a imax.
 * \param lmaxtet length of largest edge of tetra \a k.
 * \param 1 if we want to check tetra with 4 ridge metrics.
 * \param ifilt pointer to store the number of vertices filtered by the PROctree.
 * \param ns pointer toward count of splits (has to be updated)
 * \param nc pointer toward count of collapses (has to be updated)
 * \param warn pointer to store a flag that warn the user in case of
 * reallocation error.
 *
 * \return -2 for low failure (mesh has to be saved).
 *
 * \return -1 for strong failure.
 *
 * \return 0 if edge cannot be modified and if we want to pass to next loop step
 * (next element or next tetra edge)
 *
 * \return 1 if edge cannot be modified and we want to treat next elt
 *
 * \return 2 if edge has been modified and we want to treat next element.
 *
 * \return 3 if nothing has been done (no error but no edge modification either).
 *
 * Try to split \a imax edge if too large and to collapse \a imin edge if too
 * small.
 *
 */
static inline
int MMG3D_mmg3d1_delone_splcol(MMG5_pMesh mesh, MMG5_pSol met,
                               MMG3D_pPROctree *PROctree,MMG5_int k,
                               int8_t imin,double lmin,
                               int8_t imax,double lmax,double lmaxtet,
                               int8_t chkRidTet,MMG5_int *ifilt,
                               MMG5_int *ns,MMG5_int *nc,
                               int *warn ) {

  int ier;
  int8_t countMemFailure = 0;

  ier = MMG3D_mmg3d1_delone_split(mesh,met,PROctree,k,imax,lmax,lmaxtet,
                                  chkRidTet,ifilt,ns,warn,&countMemFailure);

  switch ( ier ) {
  case -2:
    /*  Low failure: try to save mesh and exit lib */
  case -1:
    /* Strong failure: exit lib without saving mesh */
  case 0:
    /* Unable to treat too large edge: pass to next edge of element or to next
     * elt */
  case 2:
    /* Edge has been splitted: pass to next element */

    /* For all previous cases, return ier value */
    return ier;
  }

  assert ( (ier==1 || ier==3) && "Check return val of delone_split");

  if ( countMemFailure > 10 ) {
    printf("  ## Error:%s: too much reallocation errors. Exit program.\n",__func__);
    return -1;
  }

  /** If unable to treat edge with ier==1 return value or if edge has
   * not been splitted but slpit_delone has not raised any error: try to
   * collapse short edge. */

  /** 2. Try to merge small edge: if collapse is not possible, pass to
   * next element */
  ier = MMG3D_adpcoledg(mesh,met,PROctree,k,imin,lmin,nc);

  /* Strong failure: ier==-1 */
  /* Unable to treat too small edge: pass to next edge of element: ier==0 */
  /* Edge has been collapsed: pass to next element: ier==2 */
  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param ne number of elements.
 * \param ifilt pointer to store the number of vertices filtered by the PROctree.
 * \param ns pointer to store the number of vertices insertions.
 * \param nc pointer to store the number of collapse.
 * \param warn pointer to store a flag that warn the user in case of
 * reallocation difficulty.
 * \return -1 if fail and we don't save the mesh, 0 if fail but we try to save
 * the mesh, 1 otherwise.
 *
 * \a adpsplcol loop: split edges longer than \ref MMG3D_LOPTL_DEL and
 * collapse edges shorter than \ref MMG3D_LOPTS.
 *
 */
static inline
int MMG5_adpsplcol(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree *PROctree,
                   MMG5_int ne,MMG5_int* ifilt,MMG5_int* ns,MMG5_int* nc,int* warn) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  double        len,lmax;
  MMG5_int      k,base;
  double        lmin;
  double        lmaxtet,lmintet;
  int           ier,imaxtet,imintet;
  int8_t        imin,imax,chkRidTet;
  static int8_t mmgWarn0 = 0;

  base = ++mesh->mark;

  if ( met->size==6 )  chkRidTet=1;
  else chkRidTet=0;

  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt)  || (pt->tag & MG_REQ) )   continue;
    else if ( pt->mark < base-2 )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /** Step 1: find longest and shortest edge  (and try to manage them) */
    imax = -1; lmax = 0.0;
    imin = -1; lmin = DBL_MAX;
    int ii;
    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      len = MMG5_lenedg(mesh,met,ii,pt);

      if ( len > lmax ) {
        lmax = len;
        imax = ii;
      }
      if ( len < lmin ) {
        lmin = len;
        imin = ii;
      }
    }
    /* Check that we have found valid edges */
    if ( imax==-1 ) {
      if ( (mesh->info.ddebug || mesh->info.imprim > 5 ) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  # Warning: %s: all edges of tetra %" MMG5_PRId " are"
                  " boundary and required.\n",
                  __func__,k);
        }
      }
      continue;
    }
    if ( imin==-1 ) {
      if ( (mesh->info.ddebug || mesh->info.imprim > 5 ) ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  # Warning: %s: all edges of tetra %" MMG5_PRId " are"
                  " boundary and required.\n",
                  __func__,k);
        }
      }
      continue;
    }

    ier = MMG3D_mmg3d1_delone_splcol(mesh,met,PROctree,k,imin,lmin,imax,
                                     lmax,lmax,chkRidTet,ifilt,ns,nc,warn);

    switch ( ier ) {
    case -2:
      /* Low failure: try to save mesh and exit lib */
      return 0;
    case -1:
      /* Strong failure: exit lib without saving mesh */
      return -1;
    case 0:
      /* Unable to treat largest/smallest edge: pass to next element */
      continue;
    case 2:
       /* Edge has been collapsed: pass to next element */
      continue;
    }

    assert ( (ier==1 || ier==3) && "Check return val of delone_split");

    /** Step 2: longest and shortest edges are stucked => try the other edges */
    imaxtet = imax;
    imintet = imin;
    lmaxtet = lmax;
    lmintet = lmin;
    assert(lmin);

    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      if ( (ii==imintet) && (lmintet < MMG3D_LOPTS)) continue;
      if ( (ii==imaxtet) && (lmaxtet > MMG3D_LOPTL_DEL) ) continue;

      len = MMG5_lenedg(mesh,met,ii,pt);

      imax = ii;
      lmax = len;
      imin = ii;
      lmin = len;

      /** 1. Try to split too long edge */
      ier = MMG3D_mmg3d1_delone_splcol(mesh,met,PROctree,k,imin,lmin,imax,
                                       lmax,lmaxtet,chkRidTet,ifilt,ns,nc,warn);

      switch ( ier ) {
      case -2:
        /* Low failure: try to save mesh and exit lib */
        return 0;
      case -1:
        /* Strong failure: exit lib without saving mesh */
        return -1;
      case 0:
        /* Unable to treat largest/smallest edge: pass to next edge of element */
        continue;
      case 3:
        /* Edge has not been splitted: pass to next edge */
        continue;
      }
      assert ( ier==2 );
      /* Edge has been splitted: pass to next element */
      break;
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization during insertion phase.
 *
 */
static
int MMG5_optbad(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree PROctree) {
  int           it,maxit;
  MMG5_int      nf,nnf,nnm,nm,nw;
  double        crit;

  /* shape optim */
  it = nnm = nnf = 0;
  maxit = 3;
  crit = 1.053;

  do {
    /* treatment of bad elements*/
    nw = MMG3D_opttyp(mesh,met,PROctree,-1);
    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,PROctree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
                __func__);
        return 0;
      }
      nnf += nf;

      nf += MMG5_swptet(mesh,met,crit,MMG3D_SWAP06,PROctree,2,mesh->mark-1);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      /* move for tria with qual<1., tetra with qual<1, internal move
       * allowed, surface degradation forbidden, volume degradation during the
       * surface move forbidden and volume degradation during volumic move
       * forbidden. Perform 1 iter max (0). */
      nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,0,mesh->mark-1);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",
                __func__);
        return 0;
      }
    }
    else  nm = 0;
    nnm += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8" MMG5_PRId " improved, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nw,nf,nm);
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter. \n",nnf,nnm,it);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param warn set to 1 if we can't insert point due to lack of memory.
 * \return -1 if fail and we dont try to end the remesh process,
 * 0 if fail but we try to end the remesh process and 1 if success.
 *
 * Split edges longer than \ref MMG3D_LOPTL_DEL and collapse edges shorter
 * than \ref MMG3D_LOPTS.
 *
 */
static
int MMG5_adpdel(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree *PROctree, int* warn) {
  int        ier;
  int        it,maxit,noptim;
  MMG5_int   ns,nc,ne,nnm,nm,nnf,nf,nnc,nns,nfilt,ifilt;
  double     maxgap,dd,declic,declicsurf;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = nfilt = 0;
  noptim = 0;
  maxit = 10;
  mesh->gap = maxgap = 0.5;
  declic = 0.5/MMG3D_ALPHAD;
  declicsurf = 1./3.46;

  do {
    if ( !mesh->info.noinsert ) {
      *warn=0;
      ns = nc = 0;
      ifilt = 0;
      ne = mesh->ne;
      ier = MMG5_adpsplcol(mesh,met,PROctree,ne,&ifilt,&ns,&nc,warn);
      if ( ier<=0 ) return -1;
    } /* End conditional loop on mesh->info.noinsert */
    else  ns = nc = ifilt = 0;

    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,*PROctree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
      nnf += nf;
      nf += MMG5_swptet(mesh,met,MMG3D_SSWAPIMPROVE,declic,*PROctree,2,mesh->mark-2);

      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    } else {
      nf = 0;
    }


    if ( !mesh->info.nomove ) {
      /* move for tria with qual<declicsurf, tetra with qual<declic, internal
       * move allowed, surface degradation forbidden, volume degradation during
       * the surface move authorized and volume degradation during volumic move
       * forbidden. Perform 2 iter max (1). */
      nm = MMG5_movtet(mesh,met,*PROctree,declicsurf,declic,1,1,0,1,1,mesh->mark-2);

      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: Unable to improve mesh.\n",__func__);
        return 0;
      }
    }
    else  nm = 0;

    nnm += nm;
    nnc += nc;
    nns += ns;
    nnf += nf;
    nfilt += ifilt;

    /* decrease size of gap for reallocation */

    if ( mesh->gap > maxgap/(double)maxit )
      mesh->gap -= maxgap/(double)maxit;
    else
      mesh->gap -= mesh->gap/(double)maxit;


    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nm+nf > 0)
      fprintf(stdout,"     %8"MMG5_PRId" filtered, %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed,"
              " %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",ifilt,ns,nc,nf,nm);

    /*optimization*/
    dd = MMG5_abs(nc-ns);
    if ( !noptim && (it==5 || ((dd < 5) || (dd < 0.05*MG_MAX(nc,ns)) || !(ns+nc))) ) {
      MMG5_optbad(mesh,met,*PROctree);
      noptim = 1;
    }

    if( it > 5 ) {
      //  if ( ns < 10 && MMG5_abs(nc-ns) < 3 )  break;
      //else if ( it > 3 && MMG5_abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
      dd = MMG5_abs(nc-ns);
      if ( dd < 5 || dd < 0.05*MG_MAX(nc,ns) )   break;
      //else if ( it > 12 && nc >= ns )  break;
    }
  }
  while( ++it < maxit && (noptim || nc+ns > 0) );

  if ( mesh->info.imprim > 0 ) {
    if ( (abs(mesh->info.imprim) < 5) && ( nnc || nns ) ) {
      fprintf(stdout,"     %8"MMG5_PRId" filtered, %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed,"
              " %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter.\n",nfilt,nns,nnc,nnf,nnm, it);
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization for LES computation (improve the element skewness).
 *
 */
static
int MMG5_optetLES(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree PROctree) {
  int       it,maxit;
  MMG5_int  nnf,nf,nw,nm,nnm;
  double    declic;

  it = nnm = nnf = 0;
  maxit = 10;
  declic = 1.01;
  ++mesh->mark;
  do {
    /* treatment of bad elements*/
    if(it < 5) {
      nw = MMG3D_opttyp(mesh,met,PROctree,mesh->mark-2);
    }
    else {
      nw = 0;
    }

    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = MMG5_swptet(mesh,met,declic,MMG3D_SWAP06,PROctree,2,mesh->mark-2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      /* move for tria with qual<1, tetra with qual<1, internal
       * move allowed, surface degradation forbidden, volume degradation during
       * the surface move forbidden and volume degradation during volumic move
       * forbidden. Perform 4 iter max (3). */
      nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,3,mesh->mark-2);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",
          __func__);
        return 0;
      }
    }
    else  nm = 0;
    nnm += nm;

//be careful, this procedure can degrade the worst elt
    if ( !mesh->info.nomove && (it==2)) {
      MMG3D_optlap(mesh,met);
    }

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8" MMG5_PRId " improved, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nw,nf,nm);
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    /* move for tria with qual<declicsurf, tetra with qual<declic, internal
     * move allowed, surface degradation forbidden, volume degradation during
       * the surface move authorized and volume degradation during volumic move
       * forbidden. Perform 4 iter max (3). */
    nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,3,mesh->mark-2);
    if ( nm < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",__func__);
      return 0;
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ) {
    fprintf(stdout,"                                            "
            "                                ");
    fprintf(stdout,"     %8" MMG5_PRId " moved\n",nm);
  }


  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter. \n",nnf,nnm,it);
  }
  return 1;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization using egde swapping and point relocation.
 *
 */
static
int MMG5_optet(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree PROctree) {
  MMG5_pTetra   pt;
  int           it,maxit;
  MMG5_int      nnf,nf,nw,k,nnm,nm;
  double        crit,declic;

  /* shape optim */
  it = nnm = nnf = 0;
  maxit  = 10;
  crit   = MMG3D_SSWAPIMPROVE;
  declic = 0.7/MMG3D_ALPHAD;
  /* mark reinitialization in order to see at least one time each tetra */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    pt->mark = mesh->mark;
  }

  do {
  ++mesh->mark;

    /* treatment of bad elements */
    if(it < 5) {
      nw = MMG3D_opttyp(mesh,met,PROctree,mesh->mark-1);
    }
    else
      nw = 0;
    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,PROctree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
      nnf += nf;

      nf += MMG5_swptet(mesh,met,crit,declic,PROctree,2,mesh->mark-1);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
          __func__);
        return 0;
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      /* move for tria with qual<1., tetra with qual<declic, internal move
       * allowed, surface degradation forbidden, volume degradation during the
       * surface move forbidden and volume degradation during volumic move
       * authorized. Perform 1 iter max (0). */
      nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,declic,1,1,1,1,0,mesh->mark-1);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",__func__);
        return 0;
      }
    }
    else  nm = 0;
    nnm += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8" MMG5_PRId " improved, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nw,nf,nm);
    }

    if ( it > 3 ) {
      if ( !nw && (!nm || !nf) )   break;
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    /* move for tria with qual<1., tetra with qual<1., internal move allowed,
     * surface degradation forbidden, volume degradation during the surface and
     * volume move forbidden. Perform 4 iter max. */
    nm = MMG5_movtet(mesh,met,PROctree,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,3,mesh->mark-2);
    if ( nm < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: Unable to improve mesh.\n",__func__);
      return 0;
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ) {
    fprintf(stdout,"                                            "
            "                                ");
    fprintf(stdout,"     %8" MMG5_PRId " moved\n",nm);
  }

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved, %d iter. \n",nnf,nnm,it);
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param permNodGlob if provided, strore the global permutation of nodes
 * \return 0 if failed, 1 otherwise.
 *
 * Analyze tetrahedra and split long / collapse short, according to
 * prescribed metric.
 *
 */
static
int MMG5_adptet_delone(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree *PROctree,
                       MMG5_int * permNodGlob) {
  MMG5_int  nnf,nf;
  int       warn,ns;

  /** Step 1: few iters of swaps */
  if ( !mesh->info.noswap ) {
    nnf = MMG5_swpmsh(mesh,met,*PROctree,2);
    if ( nnf < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to improve mesh. Exiting.\n",
              __func__);
      return 0;
    }
    nf = MMG5_swptet(mesh,met,MMG3D_SSWAPIMPROVE,MMG3D_SWAP06,*PROctree,2,mesh->mark-2);
    if ( nf < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: Unable to improve mesh. Exiting.\n",
              __func__);
      return 0;
    }
    nnf+=nf;
  } else {
    nnf = nf = 0;
  }

  if ( mesh->info.ddebug ) {
    fprintf(stdout," ------------- Delaunay: INITIAL SWAP %7"MMG5_PRId"\n",nnf);
    MMG3D_outqua(mesh,met);
  }

  /** Step 2: few iters of splits, collapses, swaps and moves */
  warn = 0;

  ns = MMG5_adpdel(mesh,met,PROctree,&warn);

  if ( ns < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to complete mesh. Exit program.\n",
      __func__);
    return 0;
  }

  if ( warn ) {
    fprintf(stderr,"\n  ## Error: %s:",__func__);
    fprintf(stderr," unable to allocate a new point in last call of adpspl.\n");
    fprintf(stderr,"  ## Check the mesh size or ");
    fprintf(stderr,"increase the maximal authorized memory with the -m option.\n");
    fprintf(stderr,"  ## Uncomplete mesh. Exiting\n" );
    return 0;
  }

  /* renumerotation if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

  /** Step 3: Last wave of improvements: few iters of bad elts treatment, swaps
   * and moves */
  if(mesh->info.optimLES) {
    if(!MMG5_optetLES(mesh,met,*PROctree)) return 0;
  }
  else {
    if(!MMG5_optet(mesh,met,*PROctree)) return 0;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param permNodGlob if provided, strore the global permutation of nodes
 * \return 0 if failed, 1 if success.
 *
 * Main adaptation routine.
 *
 */
int MMG5_mmg3d1_delone(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *permNodGlob) {
  MMG3D_pPROctree PROctree = NULL;

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !MMG3D_chkmani(mesh) ) {
    fprintf(stderr,"\n  ## Non orientable implicit surface. Exit program.\n");
    return 0;
  }

  /**--- stage 1: geometric mesh  */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !MMG5_anatet(mesh,met,1,0) ) {
    fprintf(stderr,"\n  ## Unable to split mesh. Exiting.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_ANATET1 to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_ANATET1") ) {
    printf("  ## WARNING: EXIT AFTER ANATET-1."
           " (MMG_SAVE_ANATET1 env variable is exported).\n");
    return 1;
  }

  /**--- stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");


  /* define metric map */
  if ( !MMG3D_defsiz(mesh,met) ) {
    fprintf(stderr,"\n  ## Metric undefined. Exit program.\n");
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
    if ( !MMG3D_gradsiz(mesh,met) ) {
      fprintf(stderr,"\n  ## Gradation problem. Exit program.\n");
      return 0;
    }
  }
  if ( mesh->info.hgradreq > 0. ) {
    MMG3D_gradsizreq(mesh,met);
  }

  /* Debug: export variable MMG_SAVE_GRADSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_GRADSIZ") ) {
    printf("  ## WARNING: EXIT AFTER GRADSIZ."
           " (MMG_SAVE_GRADSIZ env variable is exported).\n");
    return 1;
  }

  /*update quality*/
  if ( !MMG3D_tetraQual(mesh,met,1) ) return 0;

  if ( !MMG5_anatet(mesh,met,2,0) ) {
    fprintf(stderr,"\n  ## Unable to split mesh. Exiting.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_ANATET2 to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_ANATET2") ) {
    printf("  ## WARNING: EXIT AFTER ANATET-2."
           " (MMG_SAVE_ANATET2 env variable is exported).\n");
    return 1;
  }

  /* renumerotation if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) ) {
    return 0;
  }

  if ( mesh->info.PROctree > 0 ) {
    if ( !MMG3D_initPROctree(mesh,&PROctree,mesh->info.PROctree) ) {
      if ( PROctree ) {
        /*free PROctree*/
        MMG3D_freePROctree(mesh,&PROctree);
      }
    }
  }

  if ( !MMG5_adptet_delone(mesh,met,&PROctree,permNodGlob) ) {
    fprintf(stderr,"\n  ## Unable to adapt. Exit program.\n");
    if ( PROctree ) {
      /*free PROctree*/
      MMG3D_freePROctree(mesh,&PROctree);
    }
    return 0;
  }

  /* in test phase: check if no element with 2 bdry faces */
  if ( !MMG5_chkfemtopo(mesh) ) {
    fprintf(stderr,"\n  ## Topology of mesh unsuited for fem computations. Exit program.\n");
    if ( PROctree ) {
      /*free PROctree*/
      MMG3D_freePROctree(mesh,&PROctree);
    }
    return 0;
  }

  int ier = 1;

  if ( mesh->info.iso && !MMG3D_chkmani(mesh) ) {
    fprintf(stderr,"\n  ## Non orientable implicit surface. Exit program.\n");
    ier = 0;
  }

  if ( PROctree ) {
    /*free PROctree*/
    MMG3D_freePROctree(mesh,&PROctree);
  }

  return ier;
}

#endif
