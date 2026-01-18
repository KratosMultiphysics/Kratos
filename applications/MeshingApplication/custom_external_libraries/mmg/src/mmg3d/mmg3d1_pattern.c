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
 * \file mmg3d/mmg3d1_pattern.c
 * \brief Perform volume and surface mesh adaptation with pattern splitting.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 * Perform volume and surface mesh adaptation with pattern splitting
 * (\a MMG_PATTERN preprocessor flag set to ON).
 *
 */

#include "libmmg3d.h"
#include "inlined_functions_3d_private.h"
#include "mmg3dexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param *warn \a warn is set to 1 if we don't have enough memory to complete mesh.
 * \return -1 if failed.
 * \return number of new points.
 *
 * Split edges of length bigger than MMG3D_LOPTL.
 *
 */
static MMG5_int MMG5_adpspl(MMG5_pMesh mesh,MMG5_pSol met, int* warn) {
 MMG5_pTetra   pt;
 MMG5_pxTetra  pxt;
 MMG5_pPoint   p0,p1;
 double        len,lmax,o[3];
 MMG5_int      ns,src,k,ip,ip1,ip2;
 int64_t       list[MMG3D_LMAX+2];
 int           ier,ilist;
 int8_t        imax,j,i,i1,i2;
 int8_t        chkRidTet;
 static int8_t mmgWarn    = 0;

  *warn=0;
  ns = 0;

  if ( met->size==6 )  chkRidTet=1;
  else chkRidTet=0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* find longest edge */
    imax = -1; lmax = 0.0;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      len = MMG5_lenedg(mesh,met,i,pt);

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }
    if ( imax==-1 ) {
      if ( !mmgWarn ) {
        fprintf(stderr,
                "\n  ## Warning: %s: at least 1 tetra with 4 required"
                " or null edges.\n",__func__);
        mmgWarn = 1;
      }
      continue;
    }
    if ( lmax < MMG3D_LOPTL )  continue;

    /* proceed edges according to lengths */
    MMG3D_find_bdyface_from_edge(mesh,pt,imax,&i,&j,&i1,&i2,&ip1,&ip2,&p0,&p1);

    if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
      /* Case of a boundary face */
      if ( !(MG_GET(pxt->ori,i)) ) continue;

      ier = MMG3D_splsurfedge( mesh,met,k,pt,pxt,imax,2,chkRidTet,warn );

      if ( ier==-1 ) { return -1; }
      else if ( !ier ) { continue; }
      else if ( ier==2 ) {
        /* Unable to split due to lack of memory */
        return ns;
      }

      ++ns;
    }
    else {
      /* Case of an internal face */

      /* Skip only boundary edges but try to treat internal edges connecting bdy
       * points */
      int8_t isbdy;
      ilist = MMG5_coquil(mesh,k,imax,list,&isbdy);
      if ( !ilist ) continue;
      else if ( isbdy ) {
        continue;
      }
      else if ( ilist<0 ) {
        return -1;
      }

      o[0] = 0.5*(p0->c[0] + p1->c[0]);
      o[1] = 0.5*(p0->c[1] + p1->c[1]);
      o[2] = 0.5*(p0->c[2] + p1->c[2]);

#ifdef USE_POINTMAP
      src = p0->src;
#else
      src = 1;
#endif
      ip = MMG3D_newPt(mesh,o,MG_NOTAG,src);

      if ( !ip )  {
        /* reallocation of point table */
        MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                             *warn=1;
                             break
                             ,o,MG_NOTAG,src);
      }

      ier = 1;
      if ( met && met->m ) {
        ier = MMG5_intmet(mesh,met,k,imax,ip,0.5);
      }
      if ( !ier ) {
        MMG3D_delPt(mesh,ip);
        return -1;
      }
      else if (ier < 0 ) {
        MMG3D_delPt(mesh,ip);
        continue;
      }

      ier = MMG3D_simbulgept(mesh,met,list,ilist,ip);
      if ( ier == 1 )
        ier = MMG5_split1b(mesh,met,list,ilist,ip,1,1,0);

      if ( ier < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to split.\n",__func__);
        return -1;
      }
      else if ( ier == 0 || ier == 2 ) {
        MMG3D_delPt(mesh,ip);
      }
      else {
       ns++;
      }
    }
  }

  return ns;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \return -1 if failed.
 * \return number of deleted points.
 *
 * Collapse edges of length smaller than MMG3D_LOPTS.
 *
 */
static MMG5_int MMG5_adpcol(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  double        len,lmin;
  MMG5_int      k,nc;
  int           ier;
  int8_t        imin,i;
  static int8_t mmgWarn = 0;

  nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /** Find shortest edge */
    imin = -1; lmin = DBL_MAX;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      len = MMG5_lenedg(mesh,met,i,pt);

      if ( len < lmin ) {
        lmin = len;
        imin = i;
      }
    }
    if ( imin==-1 ) {
      if ( !mmgWarn ) {
        fprintf(stderr,
                "\n  ## Warning: %s: at least 1 tetra with 4 required"
                " or null edges.\n",__func__);
        mmgWarn = 1;
      }
      continue;
    }

    /** Try to collapse this edge */
    ier = MMG3D_adpcoledg(mesh,met,NULL,k,imin,lmin,&nc);
    switch ( ier ) {
    case -1:
      /* chkcol_bdy failure do to error in edge shell computation */
      return -1;
    default:
      assert ( ier==0 || ier==2 || ier==3 );
      /* Pass to next element for other return values:
       *  - 0: for colver failure or impossible collapse
       *  - 2: for successful Collapse
       *  - 3: if edge is large enough or chkcol_[int|bdy] refuse the collapse
       * (due to normal deviation, hausdorff check ...) */
    }
  }

  return nc;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param permNodGlob if provided, strore the global permutation of nodes.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Analyze tetrahedra and split long or collapse short edges according to
 * prescribed metric.
 *
 */
static int MMG5_adptet(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *permNodGlob) {
  int      it1,it,maxit;
  MMG5_int nf,nnf,nnm,nm,nnc,nc,nns,ns;
  int      warn;//,nw;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = warn = 0;
  maxit = 10;
  do {
    if ( !mesh->info.noinsert ) {
      ns = MMG5_adpspl(mesh,met,&warn);
      if ( ns < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to complete mesh."
                " Exit program.\n",__func__);
        return 0;
      }
    }
    else  ns = 0;

    /* renumbering if available and needed */
    if ( it==1 && !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
      return 0;

    if ( !mesh->info.noinsert ) {
      nc = MMG5_adpcol(mesh,met);
      if ( nc < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to complete mesh."
                " Exit program.\n",__func__);
        return 0;
      }
    }
    else  nc = 0;

    if ( !mesh->info.nomove ) {
      nm = MMG5_movtet(mesh,met,NULL,MMG3D_MAXKAL,MMG3D_MAXKAL,1,0,0,0,1,mesh->mark-2);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh."
                " Exiting.\n",__func__);
        return 0;
      }
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,NULL,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh."
                " Exiting.\n",__func__);
        return 0;
      }
      nnf += nf;

      nf = MMG5_swptet(mesh,met,MMG3D_SSWAPIMPROVE,MMG3D_SWAP06,NULL,2,mesh->mark-2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh."
                " Exiting.\n",__func__);
        return 0;
      }
    }
    else  nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    nnm += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",ns,nc,nf,nm);
    if ( ns < 10 && MMG5_abs(nc-ns) < 3 )  break;
    else if ( it > 3 && MMG5_abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && nc+ns > 0 );

  if ( warn ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new point in last"
            " call of MMG5_adpspl.\n",__func__);
    MMG5_INCREASE_MEM_MESSAGE();

    fprintf(stderr,"\n  ## Error: %s: uncomplete mesh."
            " Exiting\n",__func__ );
    return 0;
  }

  /* renumbering if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

  /*shape optim*/
  it1 = it;
  it  = 0;
  maxit = 2;
  do {
/*     /\* treatment of bad elements*\/ */
/*     if( 0 && it < 2) { */
/*       nw = MMG3D_opttyp(mesh,met,NULL); */
/*     } */
/*     else */
/*       nw = 0; */

    if ( !mesh->info.nomove ) {
      nm = MMG5_movtet(mesh,met,NULL,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,0,mesh->mark-2);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",
          __func__);
        return 0;
      }
      nnm += nm;
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
      nf = MMG5_swpmsh(mesh,met,NULL,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: unable to improve mesh."
                " Exiting.\n",__func__);
        return 0;
      }
      nnf += nf;

      nf = MMG5_swptet(mesh,met,MMG3D_SSWAPIMPROVE,MMG3D_SWAP06,NULL,2,mesh->mark-2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Error: %s: Unable to improve mesh."
                " Exiting.\n",__func__);
        return 0;
      }
    }
    else  nf = 0;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && /*nw+*/nf+nm > 0 ){
      fprintf(stdout,"                                            ");
      fprintf(stdout,"%8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nf,nm);
    }
  }
  while( ++it < maxit && /*nw+*/nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = MMG5_movtet(mesh,met,NULL,MMG3D_MAXKAL,MMG3D_MAXKAL,1,1,1,1,3,mesh->mark-2);
    if ( nm < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to improve mesh.\n",
              __func__);
      return 0;
    }
    nnm += nm;
  }
  else  nm = 0;

  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ){
    fprintf(stdout,"                                            ");
    fprintf(stdout,"                  %8" MMG5_PRId " moved\n",nm);
  }

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnc > 0 || nns > 0) )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved,"
              " %d iter. \n",
              nns,nnc,nnf,nnm,it+it1);
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param permNodGlob if provided, strore the global permutation of nodes.
 *
 * \return 0 if failed, 1 if success.
 *
 * Main adaptation routine.
 *
 */
int MMG5_mmg3d1_pattern(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *permNodGlob) {

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !MMG3D_chkmani(mesh) ) {
    fprintf(stderr,"\n  ## Non orientable implicit surface before remeshing. Exit program.\n");
    return 0;
  }

  /**--- stage 1: geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !MMG5_anatet(mesh,met,1,1) ) {
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

  /* renumbering if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

  /**--- Stage 2: computational mesh */
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

  /* Debug: export variable MMG_SAVE_GRADSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_GRADSIZ") ) {
    printf("  ## WARNING: EXIT AFTER GRADSIZ."
           " (MMG_SAVE_GRADSIZ env variable is exported).\n");
    return 1;
  }

  if ( mesh->info.hgrad > 0. ) {
    if ( !MMG3D_gradsiz(mesh,met) ) {
      fprintf(stderr,"\n  ## Gradation problem. Exit program.\n");
      return 0;
    }
  }
  if ( mesh->info.hgradreq > 0. ) {
    MMG3D_gradsizreq(mesh,met);
  }

  /* update quality*/
  if ( !MMG3D_tetraQual(mesh,met,1) ) return 0;

  if ( !MMG5_anatet(mesh,met,2,1) ) {
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

  /* renumbering if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

#ifdef DEBUG
  puts("---------------------------Fin anatet---------------------");
  MMG3D_outqua(mesh,met,mesh->info.optimLES);
#endif
  if ( !MMG5_adptet(mesh,met,permNodGlob) ) {
    fprintf(stderr,"\n  ## Unable to adapt. Exit program.\n");
    return 0;
  }

#ifdef DEBUG
  puts("---------------------Fin adptet-----------------");
  MMG3D_outqua(mesh,met,mesh->info.optimLES);
#endif
  /* in test phase: check if no element with 2 bdry faces */
  if ( !MMG5_chkfemtopo(mesh) ) {
    fprintf(stderr,"\n  ## Topology of mesh unsuited for fem computations. Exit program.\n");
    return 0;
  }

  if ( mesh->info.iso && !MMG3D_chkmani(mesh) ) {
    fprintf(stdout,"\n  ## Warning: %s: Non orientable implicit surface after remeshing.\n",__func__);
    return 1;
  }

  return 1;
}
