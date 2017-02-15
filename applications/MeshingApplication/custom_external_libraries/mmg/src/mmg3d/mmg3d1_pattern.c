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
 * (\a PATTERN preprocessor flag set to ON).
 *
 */

#include "inlined_functions_3d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param *warn \a warn is set to 1 if we don't have enough memory to complete mesh.
 * \return -1 if failed.
 * \return number of new points.
 *
 * Split edges of length bigger than _MMG5_LOPTL.
 *
 */
static int _MMG5_adpspl(MMG5_pMesh mesh,MMG5_pSol met, int* warn) {
  MMG5_pTetra     pt;
  MMG5_pxTetra    pxt;
  MMG5_Tria       ptt;
  MMG5_pPoint     p0,p1,ppt;
  MMG5_pxPoint    pxp;
  double     dd,len,lmax,o[3],to[3],no1[3],no2[3],v[3];
  int        k,ip,ip1,ip2,list[MMG3D_LMAX+2],ilist;
  int        ns,ref,ier;
  int16_t    tag;
  char       imax,j,i,i1,i2,ifa0,ifa1;

  *warn=0;
  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* find longest edge */
    imax = -1; lmax = 0.0;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      len = _MMG5_lenedg(mesh,met,i,pt);

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }
    if ( imax==-1 )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are required or of length null.\n",
              __FILE__,__LINE__,k);
    if ( lmax < _MMG5_LOPTL )  continue;

    /* proceed edges according to lengths */
    ifa0 = _MMG5_ifar[imax][0];
    ifa1 = _MMG5_ifar[imax][1];
    i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
    j  = _MMG5_iarfinv[i][imax];
    i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
    i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    p0  = &mesh->point[ip1];
    p1  = &mesh->point[ip2];

    /* Case of a boundary face */
    if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
      if ( !(MG_GET(pxt->ori,i)) ) continue;
      ref = pxt->edg[_MMG5_iarf[i][j]];
      tag = pxt->tag[_MMG5_iarf[i][j]];
      if ( tag & MG_REQ )  continue;
      tag |= MG_BDY;
      ilist = _MMG5_coquil(mesh,k,imax,list);
      if ( !ilist )  continue;
      else if ( ilist < 0 )
        return(-1);
      if ( tag & MG_NOM ){
        if( !_MMG5_BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
          continue;
        else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
          _MMG5_tet2tri(mesh,k,i,&ptt);
          _MMG5_nortri(mesh,&ptt,no1);
          if ( !MG_GET(pxt->ori,i) ) {
            no1[0] *= -1.0;
            no1[1] *= -1.0;
            no1[2] *= -1.0;
          }
        }
      }
      else if ( tag & MG_GEO ) {
        if ( !_MMG5_BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
          continue;
        if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
          _MMG5_tet2tri(mesh,k,i,&ptt);
          _MMG5_nortri(mesh,&ptt,no1);
          no2[0] = to[1]*no1[2] - to[2]*no1[1];
          no2[1] = to[2]*no1[0] - to[0]*no1[2];
          no2[2] = to[0]*no1[1] - to[1]*no1[0];
          dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
          if ( dd > _MMG5_EPSD2 ) {
            dd = 1.0 / sqrt(dd);
            no2[0] *= dd;
            no2[1] *= dd;
            no2[2] *= dd;
          }
        }
      }
      else if ( tag & MG_REF ) {
        if ( !_MMG5_BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
          continue;
      }
      else {
        if ( !_MMG5_norface(mesh,k,i,v) )  continue;
        if ( !_MMG5_BezierReg(mesh,ip1,ip2,0.5,v,o,no1) )
          continue;
      }
      ip = _MMG3D_newPt(mesh,o,tag);
      if ( !ip ) {
        /* reallocation of point table */
        _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                            *warn=1;
                            break
                            ,o,tag);
      }
      if ( met->m ) {
        ier = _MMG5_intmet(mesh,met,k,imax,ip,0.5);
        if ( !ier ) {
          _MMG3D_delPt(mesh,ip);
          return(-1);
        }
        else if (ier < 0) {
          _MMG3D_delPt(mesh,ip);
          continue;
        }
      }
      ier = _MMG3D_simbulgept(mesh,met,list,ilist,ip);
      if ( !ier ) {
        ier = _MMG3D_dichoto1b(mesh,met,list,ilist,ip);
      }
      ier = _MMG5_split1b(mesh,met,list,ilist,ip,1,1);
      /* if we realloc memory in _MMG5_split1b pt and pxt pointers are not valid */
      pt = &mesh->tetra[k];
      pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

      if ( ier < 0 ) {
        fprintf(stderr," ## Error: unable to split.\n");
        return(-1);
      }
      else if ( !ier ) {
        _MMG3D_delPt(mesh,ip);
        continue;
      }
      ns++;
      ppt = &mesh->point[ip];
      if ( MG_EDG(tag) || (tag & MG_NOM) )
        ppt->ref = ref;
      else
        ppt->ref = pxt->ref[i];

      pxp = &mesh->xpoint[ppt->xp];
      if ( tag & MG_NOM ){
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(ppt->n,to,3*sizeof(double));
      }
      else if ( tag & MG_GEO ) {
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(pxp->n2,no2,3*sizeof(double));
        memcpy(ppt->n,to,3*sizeof(double));
      }
      else if ( tag & MG_REF ) {
        memcpy(pxp->n1,no1,3*sizeof(double));
        memcpy(ppt->n,to,3*sizeof(double));
      }
      else
        memcpy(pxp->n1,no1,3*sizeof(double));
    }

    /* Case of an internal face */
    else {
      if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) continue;
      ilist = _MMG5_coquil(mesh,k,imax,list);
      if ( !ilist ) continue;
      else if ( ilist<0 ) return(-1);
      o[0] = 0.5*(p0->c[0] + p1->c[0]);
      o[1] = 0.5*(p0->c[1] + p1->c[1]);
      o[2] = 0.5*(p0->c[2] + p1->c[2]);

      ip = _MMG3D_newPt(mesh,o,MG_NOTAG);

      if ( !ip )  {
        /* reallocation of point table */
        _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                            *warn=1;
                            break
                            ,o,MG_NOTAG);
      }
      ppt = &mesh->point[ip];
      if ( met->m ) {
        ier = _MMG5_intmet(mesh,met,k,imax,ip,0.5);
        if ( !ier ) {
          _MMG3D_delPt(mesh,ip);
          return(-1);
        }
        else if (ier < 0 ) {
          _MMG3D_delPt(mesh,ip);
          continue;
        }
      }
      ier = _MMG5_split1b(mesh,met,list,ilist,ip,1,1);
      if ( ier < 0 ) {
        fprintf(stderr,"  ## Error: unable to split.\n");
        return(-1);
      }
      else if ( !ier ) {
        _MMG3D_delPt(mesh,ip);
      }
      else {
       ns++;
      }
    }
  }

  return(ns);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return -1 if failed.
 * \return number of deleted points.
 *
 * Collapse edges of length smaller than _MMG5_LOPTS.
 *
 */
static int _MMG5_adpcol(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra     pt;
  MMG5_pxTetra    pxt;
  MMG5_pPoint     p0,p1;
  double     len,lmin;
  int        k,ip,iq,list[MMG3D_LMAX+2],ilist,lists[MMG3D_LMAX+2],ilists,nc;
  int        ier;
  int16_t    tag;
  char       imin,j,i,i1,i2,ifa0,ifa1;

  nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;
    ier = 0;

    /* find shortest edge */
    imin = -1; lmin = DBL_MAX;
    for (i=0; i<6; i++) {
      if ( pt->xt && (pxt->tag[i] & MG_REQ) )  continue;
      len = _MMG5_lenedg(mesh,met,i,pt);

      if ( len < lmin ) {
        lmin = len;
        imin = i;
      }
    }
    if ( imin==-1 )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
              __FILE__,__LINE__,k);
    if ( lmin > _MMG5_LOPTS )  continue;

    // Case of an internal tetra with 4 ridges vertices.
    if ( lmin == 0 ) continue;

    ifa0 = _MMG5_ifar[imin][0];
    ifa1 = _MMG5_ifar[imin][1];
    i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
    j  = _MMG5_iarfinv[i][imin];
    i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
    i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
    ip = pt->v[i1];
    iq = pt->v[i2];
    p0 = &mesh->point[ip];
    p1 = &mesh->point[iq];
    if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;


    /* Case of a boundary face */
    ilist = 0;
    if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
      tag = pxt->tag[_MMG5_iarf[i][j]];
      if ( tag & MG_REQ )  continue;
      tag |= MG_BDY;
      if ( p0->tag > tag )   continue;
      if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;

      if (_MMG5_boulesurfvolp(mesh,k,i1,i,
                              list,&ilist,lists,&ilists,(p0->tag & MG_NOM)) < 0 )
        return(-1);

      ilist = _MMG5_chkcol_bdy(mesh,met,k,i,j,list,ilist,lists,ilists,2);
    }
    /* Case of an internal face */
    else {
      if ( p0->tag & MG_BDY )  continue;
      ilist = _MMG5_boulevolp(mesh,k,i1,list);
      ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,2);
    }
    if ( ilist > 0 ) {
      ier = _MMG5_colver(mesh,met,list,ilist,i2,2);
      if ( ier < 0 )  return(-1);
      else if ( ier ) {
        _MMG3D_delPt(mesh,ier);
        nc++;
      }
    }
    else if (ilist < 0 )  return(-1);
  }

  return(nc);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Analyze tetrahedra and split long or collapse short edges according to
 * prescribed metric.
 *
 */
static int _MMG5_adptet(MMG5_pMesh mesh,MMG5_pSol met) {
  int      it1,it,nnc,nns,nnf,nnm,maxit,nc,ns,nf,nm;
  int      warn;//,nw;
  double   maxgap;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = warn = 0;
  maxit = 10;
  mesh->gap = maxgap = 0.5;
  do {
    if ( !mesh->info.noinsert ) {
      ns = _MMG5_adpspl(mesh,met,&warn);
      if ( ns < 0 ) {
        fprintf(stderr,"  ## Unable to complete mesh. Exit program.\n");
        return(0);
      }
    }
    else  ns = 0;

    /* renumbering if available and needed */
    if ( it==1 && !_MMG5_scotchCall(mesh,met) )
      return(0);

    if ( !mesh->info.noinsert ) {
      nc = _MMG5_adpcol(mesh,met);
      if ( nc < 0 ) {
        fprintf(stderr,"  ## Unable to complete mesh. Exit program.\n");
        return(0);
      }
    }
    else  nc = 0;

    if ( !mesh->info.nomove ) {
      nm = _MMG5_movtet(mesh,met,NULL,1);
      if ( nm < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
      nf = _MMG5_swpmsh(mesh,met,NULL,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = _MMG5_swptet(mesh,met,1.053,NULL,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    nnm += nm;
    /* decrease size of gap for reallocation */
    if ( mesh->gap > maxgap/(double)maxit )
      mesh->gap -= maxgap/(double)maxit;
    else
      mesh->gap -= mesh->gap/(double)maxit;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved\n",ns,nc,nf,nm);
    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && nc+ns > 0 );

  if ( warn ) {
    fprintf(stderr,"  ## Error:");
    fprintf(stderr," unable to allocate a new point in last call"
            " of _MMG5_adpspl.\n");
    _MMG5_INCREASE_MEM_MESSAGE();
    fprintf(stderr,"  ## Uncomplete mesh. Exiting\n" );
    return(0);
  }

  /* renumbering if available */
  if ( !_MMG5_scotchCall(mesh,met) )
    return(0);

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
      nm = _MMG5_movtet(mesh,met,NULL,0);
      if ( nm < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh.\n");
        return(0);
      }
      nnm += nm;
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
      nf = _MMG5_swpmsh(mesh,met,NULL,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = _MMG5_swptet(mesh,met,1.053,NULL,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && /*nw+*/nf+nm > 0 ){
/*       fprintf(stdout,"                         "); */
/*       fprintf(stdout,"%8d improved, %8d swapped, %8d moved\n",nw,nf,nm); */
      fprintf(stdout,"                                            ");
      fprintf(stdout,"%8d swapped, %8d moved\n",nf,nm);
    }
  }
  while( ++it < maxit && /*nw+*/nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = _MMG5_movtet(mesh,met,NULL,3);
    if ( nm < 0 ) {
      fprintf(stderr,"  ## Unable to improve mesh.\n");
      return(0);
    }
    nnm += nm;
  }
  else  nm = 0;

  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ){
    fprintf(stdout,"                                            ");
    fprintf(stdout,"                  %8d moved\n",nm);
  }

  if ( mesh->info.imprim ) {
    if ( abs(mesh->info.imprim) < 5 && (nnc > 0 || nns > 0) )
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped, %8d moved,"
              " %d iter. \n",
              nns,nnc,nnf,nnm,it+it1);
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if failed, 1 if success.
 *
 * Main adaptation routine.
 *
 */
int _MMG5_mmg3d1_pattern(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !_MMG5_chkmani(mesh) ) {
    fprintf(stderr,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  /**--- stage 1: geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !_MMG5_anatet(mesh,met,1,1) ) {
    fprintf(stderr,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }

  /* renumbering if available */
  if ( !_MMG5_scotchCall(mesh,met) )
    return(0);

  /**--- Stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  /* define metric map */
  if ( !_MMG5_defsiz(mesh,met) ) {
    fprintf(stderr,"  ## Metric undefined. Exit program.\n");
    return(0);
  }

  if ( mesh->info.hgrad > 0. ) {
    if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",exp(mesh->info.hgrad));
    if ( !_MMG5_gradsiz(mesh,met) ) {
      fprintf(stderr,"  ## Gradation problem. Exit program.\n");
      return(0);
    }
  }

  if ( !_MMG5_anatet(mesh,met,2,1) ) {
    fprintf(stderr,"  ## Unable to split mesh. Exiting.\n");
    return(0);
  }

  /* renumbering if available */
  if ( !_MMG5_scotchCall(mesh,met) )
    return(0);

#ifdef DEBUG
  puts("---------------------------Fin anatet---------------------");
  _MMG3D_outqua(mesh,met,mesh->info.optimLES);
#endif
  if ( !_MMG5_adptet(mesh,met) ) {
    fprintf(stderr,"  ## Unable to adapt. Exit program.\n");
    return(0);
  }

#ifdef DEBUG
  puts("---------------------Fin adptet-----------------");
  _MMG3D_outqua(mesh,met,mesh->info.optimLES);
#endif
  /* in test phase: check if no element with 2 bdry faces */
  if ( !_MMG5_chkfemtopo(mesh) ) {
    fprintf(stderr,"  ## Topology of mesh unsuited for fem computations. Exit program.\n");
    return(0);
  }

  if ( mesh->info.iso && !_MMG5_chkmani(mesh) ) {
    fprintf(stderr,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  return(1);
}
