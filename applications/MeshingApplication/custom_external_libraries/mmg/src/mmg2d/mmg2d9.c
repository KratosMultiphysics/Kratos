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
 * \file mmg2d/mmg2d9.c
 * \brief Lagrangian meshing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmg2d_private.h"
#define MMG2D_DEGTOL  5.e-1

/* Calculate an estimate of the average (isotropic) length of edges in the mesh */
double MMG2D_estavglen(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    p1,p2;
  MMG5_int       k,na;
  double         len,lent,dna;
  int8_t         i,i1,i2;

  na = 0;
  lent = 0.0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];

      len = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);

      lent += sqrt(len);
      na++;
    }
  }

  dna = (double)na;
  dna = 1.0 / dna;
  lent *= dna;

  return lent;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param disp pointer toward the displacement structure.
 * \param t fraction of displacement to test
 * \param tetIdx to fill with the list of non valid tria if provided.
 *
 * \return 0 if success (movement can be achieved), 1 or the number of invalid
 * tria otherwise.
 *
 * Check if moving mesh with disp for a fraction t yields a
 * valid mesh.
 *
 */
MMG5_int MMG2D_chkmovmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t,MMG5_int *triIdx) {
  MMG5_pTria   pt;
  MMG5_pPoint  ppt;
  double       *v,c[3][2],tau;
  MMG5_int     k,np,idx;
  int8_t       i,j;

  /* Pseudo time-step = fraction of disp to perform */
  tau = (double)t / MMG2D_SHORTMAX;
  idx = 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      v = &disp->m[2*np];
      for (j=0; j<2; j++)
        c[i][j] = ppt->c[j]+tau*v[j];
    }

    //     Other criteria : eg. a rate of degradation, etc... ?
    if( MMG2D_caltri_iso_3pt(c[0],c[1],c[2]) < MMG2D_NULKAL) {
      if ( triIdx ) {
        triIdx[idx++] = k;
      }
      else {
        return 1;
      }
    }
  }

  return idx;
}

/** Perform mesh motion along disp, for a fraction t, and the corresponding updates */
int MMG2D_dispmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t,int itdeg) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  double        *v,tau,ctau,c[3][2],ocal,ncal;
  MMG5_int      k,np;
  int8_t        i,j;

  tau = (double)t /MMG2D_SHORTMAX;
  ctau = 1.0 - tau;

  /* Identify elements which are very distorted in the process */
  for (k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      for (j=0; j<2; j++)
        c[i][j] = ppt->c[j];
    }

    ocal = MMG2D_caltri_iso_3pt(c[0],c[1],c[2]);

    for (i=0; i<3; i++) {
      np = pt->v[i];
      v = &disp->m[2*np];
      for (j=0; j<2; j++)
        c[i][j] += tau*v[j];
    }

    ncal = MMG2D_caltri_iso_3pt(c[0],c[1],c[2]);

    if ( ncal < MMG2D_DEGTOL*ocal )
      pt->cc = itdeg;

  }

  /* Perform physical displacement */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];

    if ( !MG_VOK(ppt) ) continue;
    v = &disp->m[2*k];

    for (i=0; i<2; i++) {
      ppt->c[i] = ppt->c[i] + tau*v[i];
      v[i] *= ctau;
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param disp pointer toward the displacement structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \param *warn \a warn is set to 1 if not enough memory is available to complete mesh.
 * \return -1 if failed.
 * \return number of new points.
 *
 * Split edges of length bigger than MMG5_LOPTL, in the Lagrangian mode.
 * Only affects triangles with cc itdeg
 *
 */
MMG5_int MMG2D_spllag(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met,MMG5_int itdeg,int *warn) {
  MMG5_pTria      pt;
  MMG5_pPoint     p1,p2;
  double          hma2,lmax,len;
  MMG5_int        k,ns,ip,ip1,ip2;
  int8_t          i,i1,i2,imax,ier;
  static int8_t   mmgWarn0=0;

  *warn = 0;
  ns    = 0;
  hma2  = mesh->info.hmax*mesh->info.hmax;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->cc != itdeg ) continue;

    /* Find the longest, internal edge */
    imax = -1;
    lmax = -1.0;

    for (i=0; i<3; i++) {
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      len = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }

    if ( (imax == -1) && (!mmgWarn0) ) {
      mmgWarn0=1;
      fprintf(stderr,"\n  ## Warning: %s: at least 1 tria whose all edges"
              " are required or of length null.\n",__func__);
    }

    if ( lmax < hma2 )  continue;
    else if ( MG_SIN(pt->tag[imax]) ) continue;

    /* Check the feasibility of splitting */
    i1 = MMG5_inxt2[imax];
    i2 = MMG5_iprv2[imax];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];

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

      /* if we realloc memory in split1b pt pointer is not valid aymore. */
      ns += ier;
    }

    /* Interpolate metric, if any */
    if ( met->m )
      met->m[ip] = 0.5*(met->m[ip1]+met->m[ip2]);

    /* Interpolate displacement */
    if ( disp->m ) {
      for (i=0; i<2; i++)
        disp->m[2*ip+i] = 0.5*(disp->m[2*ip1+i]+disp->m[2*ip2+i]);
    }
  }

  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \return -1 if failed.
 * \return number of collapsed points.
 *
 * Attempt to collapse small internal edges in the Lagrangian mode; only affects tria with cc itdeg.
 *
 */
static int MMG2D_coleltlag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  MMG5_pTria     pt;
  MMG5_pPoint    p1,p2;
  double         hmi2,len;
  MMG5_int       nc,k;
  int            ilist;
  MMG5_int       list[MMG5_TRIA_LMAX+2];
  int8_t         i,i1,i2,open;

  nc    = 0;
  hmi2  = mesh->info.hmin*mesh->info.hmin;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->cc != itdeg ) continue;

    for (i=0; i<3; i++) {
      if ( MG_SIN(pt->tag[i]) ) continue;

      open = ( mesh->adja[3*(k-1)+1+i] == 0 ) ? 1 : 0;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];

      if ( MG_SIN(p1->tag) ) continue;
      else if ( p1->tag & MG_GEO ) {
        if ( ! (p2->tag & MG_GEO) || !(pt->tag[i] & MG_GEO) ) continue;
      }

      /* Check length */
      len = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]);
      if ( len > hmi2 )  continue;

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

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param crit coefficient of quality improvment.
 * \param itdeg degraded elements.
 *
 * Internal edge flipping in the Lagrangian mode; only affects trias with cc itdeg
 *
 */
MMG5_int MMG2D_swpmshlag(MMG5_pMesh mesh,MMG5_pSol met,double crit,int itdeg) {
  MMG5_pTria   pt;
  int          it,maxit;
  int8_t       i;
  MMG5_int     k,ns,nns;

  maxit = 2;
  it    = 0;
  nns   = 0;

  do {
    ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;
      if ( pt->cc != itdeg ) continue;

      for (i=0; i<3; i++) {
        /* Prevent swap of a ref or tagged edge */
        if ( MG_SIN(pt->tag[i]) || MG_EDG(pt->tag[i]) ) continue;

        else if ( MMG2D_chkswp(mesh,met,k,i,2) ) {
          ns += MMG2D_swapar(mesh,k,i);
          break;
        }

      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );

  return nns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \return -1 if failed, number of moved points otherwise.
 *
 * Analyze trias with cc = itdeg and move internal points so as to make mesh more uniform.
 *
 */
MMG5_int MMG2D_movtrilag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  MMG5_pTria        pt;
  MMG5_pPoint       p0;
  int               it,maxit,ilist;
  MMG5_int          k,base,list[MMG5_TRIA_LMAX+2],nm,nnm;
  int8_t            i,ier;

  nnm   = 0;
  it    = 0;
  maxit = 5;

  /* Reset point flags */
  base = 1;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  do {
    base++;
    nm = 0;

    for(k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;
      if ( pt->cc != itdeg ) continue;

      for (i=0; i<3; i++) {
        p0 = &mesh->point[pt->v[i]];
        if ( p0->flag == base || MG_SIN(p0->tag) || p0->tag & MG_NOM ) continue;

        int8_t dummy;
        ilist = MMG5_boulet(mesh,k,i,list,0,&dummy);

        if ( MG_EDG(p0->tag) )
          ier = MMG2D_movedgpt(mesh,met,ilist,list,0);
        else
          ier = MMG2D_movintpt(mesh,met,ilist,list,0);

        if ( ier ) {
          nm++;
          p0->flag = base;
        }
      }
    }
    nnm += nm;
  }
  while (++it < maxit && nm > 0 );

  return nnm;
}

/**
 * \param mesh mesh structure
 * \param disp displacement structure
 * \param met metric structure
 * \param invalidTrias array to store the list of invalid tria if we are unable to move
 *
 * \return 0 if fail, 1 if success to move, the opposite of the number of non
 * valid trias if we can't move (-ninvalidTrias).
 *
 * Lagrangian node displacement and meshing.
 * Code for options: info.lag >= 0 -> displacement,
 *                   info.lag > 0  -> displacement+remeshing with swap and moves
 *                   info.lag > 1  -> displacement+remeshing with split+collapse+swap+move
 *
 */
int MMG2D_mmg2d9(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met,MMG5_int **invalidTrias) {
  double             avlen,tau,hmintmp,hmaxtmp;
  int                itmn,itdc,maxitmn,maxitdc,iit,warn;
  MMG5_int           nspl,nnspl,nnnspl,nc,nnc,nnnc,ns,nns,nnns,nm,nnm,nnnm;
  short              t,lastt;
  int8_t             ier;
  MMG5_int           k,ninvalidTrias;

  maxitmn = 10;
  maxitdc = 100;
  t = 0;
  tau = 0.0;
  ninvalidTrias = 0;

  nnnspl = nnnc = nnns = nnnm = lastt = 0;

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** LAGRANGIAN MOTION\n");

  /* Field cc stores information about whether a triangle has been greatly distorted during current step */
  for (k=1; k<=mesh->nt; k++)
    mesh->tria[k].cc = 0;

  /* Estimate of the average, maximum and minimum edge lengths */
  avlen = MMG2D_estavglen(mesh);

  hmintmp = mesh->info.hmin;
  hmaxtmp = mesh->info.hmax;

  mesh->info.hmax = MMG2D_LLONG*avlen;
  mesh->info.hmin = MMG2D_LSHRT*avlen;

  for (itmn=0; itmn<maxitmn; itmn++) {

#ifdef USE_ELAS
    /* Extension of the displacement field */
    if ( !MMG2D_velextLS(mesh,disp) ) {
      fprintf(stderr,"\n  ## Problem in func. MMG2D_velextLS. Exit program.\n");
      return 0;
    }
#else
    fprintf(stderr,"\n  ## Error: %s: you need to compile with the USE_ELAS"
            " CMake's flag set to ON to use the rigidbody movement.\n",__func__);
    return 0;
#endif

    // MMG5_saveDisp(mesh,disp);

    /* Sequence of dichotomy loops to find the largest admissible displacements */
    for (itdc=0; itdc<maxitdc; itdc++) {
      nnspl = nnc = nns = nnm = 0;

      t = MMG5_dikmov(mesh,disp,&lastt,MMG2D_SHORTMAX,MMG2D_chkmovmesh);
      if ( t == 0 ) {
        if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
          fprintf(stderr,"\n   *** Stop: impossible to proceed further\n");
        break;
      }

      ier = MMG2D_dispmesh(mesh,disp,t,itdc);
      if ( !ier ) {
        fprintf(stderr,"\n  ** Impossible motion\n");
        return 0;
      }

      tau = tau + ((double)t /MMG2D_SHORTMAX)*(1.0-tau);
      if ( (abs(mesh->info.imprim) > 3 ) || mesh->info.ddebug )
        printf("   ---> Realized displacement: %f\n",tau);

      /* Local remeshing depending on the option */
      if ( mesh->info.lag > 0 ) {
        for (iit=0; iit<5; iit++) {

          nspl = nc = ns = nm = 0;

          if ( mesh->info.lag > 1 ) {
            if ( !mesh->info.noinsert ) {

              /* Split of points */
              nspl = MMG2D_spllag(mesh,disp,met,itdc,&warn);
              if ( nspl < 0 ) {
                fprintf(stderr,"\n  ## Problem in spllag. Exiting.\n");
                return 0;
              }

              /* Collapse of points */
              nc = MMG2D_coleltlag(mesh,met,itdc);
              if ( nc < 0 ) {
                fprintf(stderr,"\n  ## Problem in coltetlag. Exiting.\n");
                return 0;
              }
            }
          }

          /* Swap of edges in tria that have resulted distorted from the process */
          /* I do not know whether it is safe to put NULL in metric here (a
           * priori ok, since there is no vertex creation or suppression) */
          if ( !mesh->info.noswap ) {
            ns = MMG2D_swpmshlag(mesh,met,1.1,itdc);
            if ( ns < 0 ) {
              fprintf(stderr,"  ## Problem in swapeltlag. Exiting.\n");
              return 0;
            }
          }
          /* Relocate vertices of tria which have been distorted in the displacement process */
          if ( !mesh->info.nomove ) {
            nm = MMG2D_movtrilag(mesh,met,itdc);
            if ( nm < 0 ) {
              fprintf(stderr,"  ## Problem in moveltlag. Exiting.\n");
              return 0;
            }
          }

          if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && (nspl+nc+ns+nm > 0) )
            printf(" %" MMG5_PRId " edges splitted, %" MMG5_PRId " vertices collapsed, %" MMG5_PRId " elements"
                   " swapped, %" MMG5_PRId " vertices moved.\n",nspl,nc,ns,nm);
          nnspl+= nspl;
          nnm  += nm;
          nnc  += nc;
          nns  += ns;
        }
        /* Five iterations of local remeshing have been performed: print final stats */
        if ( abs(mesh->info.imprim) > 3 && abs(mesh->info.imprim) < 5 && (nnspl+nnm+nns+nnc > 0) )
          printf(" %" MMG5_PRId " edges splitted, %" MMG5_PRId " vertices collapsed, %" MMG5_PRId " elements"
                 " swapped, %" MMG5_PRId " vertices moved.\n",nnspl,nnc,nns,nnm);

      }

      nnnspl += nnspl;
      nnnm   += nnm;
      nnnc   += nnc;
      nnns   += nns;

      if ( t == MMG2D_SHORTMAX ) break;
    }
    /* End of dichotomy loop: maximal displacement of the extended velocity
     * field has been performed */
    if ( mesh->info.imprim > 1 && abs(mesh->info.imprim) < 4 ) {
      printf("   ---> Realized displacement: %f\n",tau);
    }

    if ( abs(mesh->info.imprim) > 2 && mesh->info.lag )
      printf(" %" MMG5_PRId " edges splitted, %" MMG5_PRId " vertices collapsed, %" MMG5_PRId " elements"
             " swapped, %" MMG5_PRId " vertices moved.\n",nnnspl,nnnc,nnns,nnnm);

    if ( (t == MMG2D_SHORTMAX) || (t==0 && itdc==0) ) break;
  }

  /* Reinsert standard values for hmin, hmax */
  mesh->info.hmin = hmintmp;
  mesh->info.hmax = hmaxtmp;

  /* If mesh optim with insertion and collapse, perform a new analysis of the
   * MMG2D_DISPREF boundary */
  if ( mesh->info.lag >= 2 ) {
    /* Identify singularities in the mesh */
    if ( !MMG2D_singul(mesh,MMG5_DISPREF) ) {
      fprintf(stderr,"\n  ## Problem in identifying singularities. Exit program.\n");
      return 0;
    }

    /* Define normal vectors at vertices on curves */
    if ( !MMG2D_norver(mesh,MMG5_DISPREF) ) {
      fprintf(stderr,"\n  ## Problem in calculating normal vectors. Exit program.\n");
      return 0;
    }
  }

  if ( tau < MMG5_EPSD2 ) {
    MMG5_SAFE_CALLOC(*invalidTrias,mesh->np,MMG5_int,
                     printf("## Warning: Not enough memory to keep track of"
                            " the invalid triangles.\n");
                     MMG5_DEL_MEM(mesh,disp->m);
                     return 1);
    ninvalidTrias = MMG2D_chkmovmesh(mesh,disp,lastt,*invalidTrias);
    assert ( ninvalidTrias );
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,disp->m);

  if ( ninvalidTrias ) {
    return -ninvalidTrias;
  }
  else {
    return 1;
  }
}
