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

#include "mmg2d.h"
//#include "ls_calls.h"
#define MMG2D_DEGTOL  5.e-1

/* Calculate an estimate of the average (isotropic) length of edges in the mesh */
double MMG2D_estavglen(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    p1,p2;
  int            k,na;
  double         len,lent,dna;
  char           i,i1,i2;

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

/** Compute quality of a triangle from the datum of its 3 vertices */
static
inline double MMG2D_caltri_iso_3pt(double *a,double *b,double *c) {
  double        abx,aby,acx,acy,bcx,bcy,area,h1,h2,h3,hm;

  abx = b[0] - a[0];
  aby = b[1] - a[1];
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];

  /* orientation */
  area = abx*acy - aby*acx;
  if ( area <= 0.0 ) return 0.0;

  /* edge lengths */
  h1 = abx*abx + aby*aby;
  h2 = acx*acx + acy*acy;
  h3 = bcx*bcx + bcy*bcy;

  hm = h1 + h2 + h3;

  if ( hm > MMG2D_EPSD ) {
    return  area / hm;
  }
  else {
    return 0.0;
  }
}

/** Check if moving mesh with disp for a fraction t yields a valid mesh */
int MMG2D_chkmovmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t) {
  MMG5_pTria   pt;
  MMG5_pPoint  ppt;
  double       *v,c[3][2],tau;
  int          k,np;
  char         i,j;

  /* Pseudo time-step = fraction of disp to perform */
  tau = (double)t / MMG2D_SHORTMAX;

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

    if( MMG2D_caltri_iso_3pt(c[0],c[1],c[2]) < MMG2D_NULKAL) return 0;  //     Other criteria : eg. a rate of degradation, etc... ?
  }

  return 1;
}

/* Return the largest fraction t that makes the motion along disp valid */
short MMG2D_dikomv(MMG5_pMesh mesh,MMG5_pSol disp) {
  int     it,maxit;
  short   t,tmin,tmax;
  char    ier;

  maxit = 200;
  it    = 0;

  tmin  = 0;
  tmax  = MMG2D_SHORTMAX;

  /* If full displacement can be achieved */
  if ( MMG2D_chkmovmesh(mesh,disp,tmax) )
    return tmax;

  /* Else, find the largest displacement by dichotomy */
  while( tmin != tmax && it < maxit ) {
    t = (tmin+tmax)/2;

    /* Case that tmax = tmin +1 : check move with tmax */
    if ( t == tmin ) {
      ier = MMG2D_chkmovmesh(mesh,disp,tmax);
      if ( ier )
        return tmax;
      else
        return tmin;
    }

    /* General case: check move with t */
    ier = MMG2D_chkmovmesh(mesh,disp,t);
    if ( ier )
      tmin = t;
    else
      tmax = t;

    it++;
  }

  return tmin;
}

/** Perform mesh motion along disp, for a fraction t, and the corresponding updates */
int MMG2D_dispmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t,int itdeg) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  double        *v,tau,ctau,c[3][2],ocal,ncal;
  int           k,np;
  char          i,j;

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
int MMG2D_spllag(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met,int itdeg,int *warn) {
  MMG5_pTria      pt;
  MMG5_pPoint     p1,p2;
  double          hma2,lmax,len;
  int             k,ns,ip,ip1,ip2;
  char            i,i1,i2,imax,ier;
  static char     mmgWarn0=0;

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
 * Attempt to collapse small internal edges in the Lagrangian mode; only affects tetras with cc itdeg.
 *
 */
static int MMG2D_coleltlag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  MMG5_pTria     pt;
  MMG5_pPoint    p1,p2;
  double         hmi2,len;
  int            nc,k,ilist,list[MMG2D_LONMAX+2];
  char           i,i1,i2,open;

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
int MMG2D_swpmshlag(MMG5_pMesh mesh,MMG5_pSol met,double crit,int itdeg) {
  MMG5_pTria   pt;
  int          k,it,maxit,ns,nns;
  char         i;

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
/** For debugging purposes: save disp */
int MMG2D_saveDisp(MMG5_pMesh mesh,MMG5_pSol disp) {
  FILE        *out;
  int         k;
  char        data[256],*ptr;

  strcpy(data,"disp.sol");
  ptr = strstr(data,".sol");
  if(ptr) *ptr = '\0';
  strcat(data,"disp.sol");

  out = fopen(data,"w");
  printf("save disp\n");
  fprintf(out,"MeshVersionFormatted 1\n\nDimension\n%d\n\n",disp->dim);
  fprintf(out,"SolAtVertices\n%d\n 1 2\n",disp->np);

  /* Print solutions */
  for(k=1; k<= disp->np; k++) {
    fprintf(out,"%f %f\n",disp->m[2*k+0],disp->m[2*k+1]);
  }

  fprintf(out,"\nEnd");
  fclose(out);

  return 1;
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
int MMG2D_movtrilag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  MMG5_pTria        pt;
  MMG5_pPoint       p0;
  int               k,it,base,maxit,nm,nnm,ilist,list[MMG2D_LONMAX+2];
  char              i,ier;

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

        ilist = MMG2D_boulet(mesh,k,i,list);

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
 *
 * \return 0 if fail, 1 if success
 *
 * Lagrangian node displacement and meshing.
 * Code for options: info.lag >= 0 -> displacement,
 *                   info.lag > 0  -> displacement+remeshing with swap and moves
 *                   info.lag > 1  -> displacement+remeshing with split+collapse+swap+move
 *
 */
int MMG2D_mmg2d9(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met) {
  double             avlen,tau,hmintmp,hmaxtmp;
  int                k,itmn,itdc,maxitmn,maxitdc,iit,warn;
  int                nspl,nnspl,nnnspl,nc,nnc,nnnc,ns,nns,nnns,nm,nnm,nnnm;
  short              t;
  char               ier;

  maxitmn = 10;
  maxitdc = 100;
  t = 0;
  tau = 0.0;

  nnnspl = nnnc = nnns = nnnm = 0;

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

  for (itmn=1; itmn<=maxitmn; itmn++) {

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
    //MMG2D_saveDisp(mesh,disp);
    /* Sequence of dichotomy loops to find the largest admissible displacements */
    for (itdc=1; itdc<=maxitdc; itdc++) {
      nnspl = nnc = nns = nnm = 0;

      t = MMG2D_dikomv(mesh,disp);
      if ( t == 0 ) {
        if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
          printf("   *** Stop: impossible to proceed further\n");
        break;
      }

      ier = MMG2D_dispmesh(mesh,disp,t,itdc);
      if ( !ier ) {
        fprintf(stdout,"  ** Impossible motion\n");
        return 0;
      }

      tau = tau + ((double)t /MMG2D_SHORTMAX)*(1.0-tau);
      if ( (abs(mesh->info.imprim) > 3 ) || mesh->info.ddebug )
        printf("   ---> Realized displacement: %f\n",tau);

      /* Local remeshing depending on the option */
      if ( mesh->info.lag > 0 ) {
        for (iit=0; iit<5; iit++) {

          nspl = nc = ns = nm = 0;

          if ( !mesh->info.noinsert ) {

            /* Split of points */
            nspl = MMG2D_spllag(mesh,disp,met,itdc,&warn);
            if ( nspl < 0 ) {
              fprintf(stdout,"  ## Problem in spllag. Exiting.\n");
              return 0;
            }

            /* Collapse of points */
            nc = MMG2D_coleltlag(mesh,met,itdc);
            if ( nc < 0 ) {
              fprintf(stdout,"  ## Problem in coltetlag. Exiting.\n");
              return 0;
            }
          }

          /* Swap of edges in tetra that have resulted distorted from the process */
          /* I do not know whether it is safe to put NULL in metric here (a
           * priori ok, since there is no vertex creation or suppression) */
          if ( !mesh->info.noswap ) {
            ns = MMG2D_swpmshlag(mesh,met,1.1,itdc);
            if ( ns < 0 ) {
              fprintf(stdout,"  ## Problem in swaptetlag. Exiting.\n");
              return 0;
            }
          }
          /* Relocate vertices of tetra which have been distorted in the displacement process */
          if ( !mesh->info.nomove ) {
            nm = MMG2D_movtrilag(mesh,met,itdc);
            if ( nm < 0 ) {
              fprintf(stdout,"  ## Problem in movtetlag. Exiting.\n");
              return 0;
            }
          }

          if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && (nspl+nc+ns+nm > 0) )
            printf(" %d edges splitted, %d vertices collapsed, %d elements"
                   " swapped, %d vertices moved.\n",nspl,nc,ns,nm);
          nnspl+= nspl;
          nnm  += nm;
          nnc  += nc;
          nns  += ns;
        }

        if ( abs(mesh->info.imprim) > 3 && abs(mesh->info.imprim) < 5 && (nnspl+nnm+nns+nnc > 0) )
          printf(" %d edges splitted, %d vertices collapsed, %d elements"
                 " swapped, %d vertices moved.\n",nnspl,nnc,nns,nnm);

      }

      nnnspl += nnspl;
      nnnm   += nnm;
      nnnc   += nnc;
      nnns   += nns;

      if ( t == MMG2D_SHORTMAX ) break;
    }
    if ( mesh->info.imprim > 0 && abs(mesh->info.imprim) < 4 ) {
      printf("   ---> Realized displacement: %f\n",tau);
      if ( abs(mesh->info.imprim) > 2 )
        printf(" %d edges splitted, %d vertices collapsed, %d elements"
               " swapped, %d vertices moved.\n",nnnspl,nnnc,nnns,nnnm);
    }

    if ( t == MMG2D_SHORTMAX ) break;
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

  /* Clean memory */
  MMG5_DEL_MEM(mesh,disp->m);

  return 1;
}
