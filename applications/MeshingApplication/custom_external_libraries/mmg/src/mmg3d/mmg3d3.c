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
 * \file mmg3d/mmg3d3.c
 * \brief Lagrangian meshing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#ifdef USE_ELAS

#include "mmg3d.h"
#include "ls_calls.h"
#include "inlined_functions_3d.h"

#define MMG5_DEGTOL  1.e-1

extern char  ddb;

/** Calculate an estimate of the average (isotropic) length of edges in the mesh */
double MMG5_estavglen(MMG5_pMesh mesh) {
  MMG5_pTetra    pt;
  MMG5_pPoint    p1,p2;
  int       k,na;
  double    len,lent,dna;
  char      i,i1,i2;

  na = 0;
  lent = 0.0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (i=0; i<6; i++) {
      i1 = MMG5_iare[i][0];
      i2 = MMG5_iare[i][1];

      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];

      len = (p2->c[0] - p1->c[0])*(p2->c[0] - p1->c[0])
        + (p2->c[1] - p1->c[1])*(p2->c[1] - p1->c[1])
        + (p2->c[2] - p1->c[2])*(p2->c[2] - p1->c[2]);

      lent += sqrt(len);
      na++;
    }
  }

  dna = (double)na;
  dna = 1.0 / dna;
  lent *= dna;

  return lent;
}

/** Interpolate displacement between v1 and v2 at intermediate position 0<=t<=1 */
static
inline int MMG5_intdispvol(double *v1, double *v2, double *vp, double t) {
  char i;

  for(i=0; i<3; i++)
    vp[i] = (1.0-t)*v1[i] + t*v2[i];

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param disp pointer toward the displacement structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \param *warn \a warn is set to 1 if we don't have enough memory to complete mesh.
 * \return -1 if failed.
 * \return number of new points.
 *
 * Split edges of length bigger than MMG3D_LOPTL, in the Lagrangian mode.
 *
 */
static int MMG5_spllag(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met,int itdeg, int* warn) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_pPoint  p0,p1;
  double       len,lmax,o[3],hma2;
  double      *m1,*m2,*mp;
  int          k,ip,ip1,ip2,list[MMG3D_LMAX+2],ilist,ns,ier,iadr;
  char         imax,i,i1,i2;
  static char  mmgWarn0 = 0;

  *warn=0;
  ns = 0;
  hma2 = mesh->info.hmax*mesh->info.hmax;

  for (k=1; k<=mesh->ne; k++) {
    pt  = &mesh->tetra[k];
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : NULL;

    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;
    if ( pt->mark != itdeg ) continue;

    /* find longest, internal edge */
    imax = -1;
    lmax = 0.0;
    for (i=0; i<6; i++) {
      i1  = MMG5_iare[i][0];
      i2  = MMG5_iare[i][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      p0 = &mesh->point[ip1];
      p1 = &mesh->point[ip2];

      /* Skip the non-internal edges */
      if ( pxt && (pxt->tag[i] & MG_BDY) )  continue;

      len = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0])
        + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])
        + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }
    if ( imax==-1 ) {
      if ( !mmgWarn0 ){
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Warning: %s: all edges of tetra %d are required"
                " or of length null.\n",__func__,k);
      }
      continue;
    }

    if ( lmax < hma2 )  continue;

    /* proceed edges according to lengths */
    i1  = MMG5_iare[imax][0];
    i2  = MMG5_iare[imax][1];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    p0 = &mesh->point[ip1];
    p1 = &mesh->point[ip2];

    /* Deal only with internal faces */
#ifndef NDEBUG
    if ( pxt ) {
      assert( !(pxt->tag[imax] & MG_BDY) ); }
#endif
    ilist = MMG5_coquil(mesh,k,imax,list);

    if ( !ilist ) continue;
    else if ( ilist<0 ) return -1;

    o[0] = 0.5*(p0->c[0] + p1->c[0]);
    o[1] = 0.5*(p0->c[1] + p1->c[1]);
    o[2] = 0.5*(p0->c[2] + p1->c[2]);

    ip = MMG3D_newPt(mesh,o,MG_NOTAG);

    if ( !ip )  {
      /* reallocation of point table */
      MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,*warn=1;break,o,MG_NOTAG);
    }

    /* Interpolation of metric, if any */
    if ( met->m ) {
      if ( !MMG5_intmet(mesh,met,k,imax,ip,0.5) ) {
        MMG3D_delPt(mesh,ip);
        return -1;
      }
    }

    /* Interpolation of displacement */
    if ( disp->m ) {
      iadr = disp->size*ip1;
      m1 = &disp->m[iadr];
      iadr = disp->size*ip2;
      m2 = &disp->m[iadr];
      iadr = disp->size*ip;
      mp = &disp->m[iadr];

      if ( !MMG5_intdispvol(m1,m2,mp,0.5) ) {
        MMG3D_delPt(mesh,ip);
        return -1;
      }
    }

    /* Il y a un check sur la taille des arêtes ici aussi ! */
    ier = MMG5_split1b(mesh,met,list,ilist,ip,1,1,0);
    if ( ier < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to split.\n",__func__);
      return -1;
    }
    else if ( !ier ) {
      MMG3D_delPt(mesh,ip);
    }
    else {
      ns++;
    }
  }

  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param crit coefficient of quality improvment.
 * \param PROctree pointer toward the PROctree structure in delaunay mode and
 * toward the \a NULL pointer otherwise.
 * \param itdeg degraded elements.
 *
 * \return -1 if fail, he number of swap otherwise.
 *
 * Internal edge flipping in the Lagrangian mode; only affects tetra marked with it
 *
 */
int MMG5_swptetlag(MMG5_pMesh mesh,MMG5_pSol met,double crit,MMG3D_pPROctree PROctree,int itdeg) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  int      list[MMG3D_LMAX+2],ilist,k,it,nconf,maxit,ns,nns,ier;
  char     i;

  maxit = 2;
  it = nns = 0;

  do {
    ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
      if ( pt->mark != itdeg ) continue;

      for (i=0; i<6; i++) {
        /* Prevent swap of a ref or tagged edge */
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          if ( pxt->edg[i] || pxt->tag[i] ) continue;
        }

        nconf = MMG5_chkswpgen(mesh,met,k,i,&ilist,list,crit,2);

        if ( nconf<0 ) return -1;
        else if ( nconf ) {
          ier = MMG5_swpgen(mesh,met,nconf,ilist,list,PROctree,2);
          if ( ier > 0 )  ns++;
          else if ( ier < 0 ) return -1;
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
 * Analyze tetrahedra marked with it and move internal points so as to make mesh more uniform.
 *
 */
int MMG5_movtetlag(MMG5_pMesh mesh,MMG5_pSol met, int itdeg) {
  MMG5_pTetra        pt;
  MMG5_pPoint        ppt;
  int           k,ier,nm,nnm,listv[MMG3D_LMAX+2],ilistv,it;
  unsigned char i,base;
  int           maxit;

  base = 1;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  maxit = 2;
  it = nnm = 0;

  do {
    base++;
    nm = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
      if ( pt->mark != itdeg ) continue;

      /* point i */
      for (i=0; i<4; i++) {
        ppt = &mesh->point[pt->v[i]];
        if ( ppt->flag == base )  continue;
        else if ( ppt->tag & MG_BDY ) continue;

        ilistv = MMG5_boulevolp(mesh,k,i,listv);
        if ( !ilistv )  continue;

        ier = MMG5_movintpt_iso(mesh,met, NULL, listv,ilistv,0);

        if ( ier ) {
          nm++;
          ppt->flag = base;

          /* Somehow interpolate displacement ? */
        }
      }
    }
    nnm += nm;
  }
  while( ++it < maxit && nm > 0 );

  return nnm;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param itdeg degraded elements.
 * \return -1 if failed.
 * \return number of collapsed points.
 *
 * Attempt to collapse small internal edges in the Lagrangian mode; only affects tetras marked with it.
 *
 */
static int MMG5_coltetlag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  MMG5_pTetra     pt;
  MMG5_pPoint     p0,p1;
  double     ll,ux,uy,uz,hmi2;
  int        k,nc,list[MMG3D_LMAX+2],ilist,base,nnm;
  int        ier;
  char       i,j,ip,iq,isnm;

  nc = nnm = 0;
  hmi2 = mesh->info.hmin*mesh->info.hmin;

  /* init of point flags, otherwise it can be uninitialized */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  for (k=1; k<=mesh->ne; k++) {
    base = ++mesh->base;
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;
    if ( pt->mark != itdeg ) continue;

    for (i=0; i<4; i++) {
      ier = 0;
      for (j=0; j<3; j++) {
        ip = MMG5_idir[i][MMG5_inxt2[j]];
        iq = MMG5_idir[i][MMG5_iprv2[j]];

        p0 = &mesh->point[pt->v[ip]];
        p1 = &mesh->point[pt->v[iq]];
        if ( p0->flag == base )  continue;
        else if ( p0->tag & MG_BDY ) continue;
        else if ( (p0->tag & MG_REQ) || (p0->tag > p1->tag) )  continue;

        /* check length */
        ux = p1->c[0] - p0->c[0];
        uy = p1->c[1] - p0->c[1];
        uz = p1->c[2] - p0->c[2];
        ll = ux*ux + uy*uy + uz*uz;

        if ( ll > hmi2 )  continue;

        isnm = 0;
        ilist = MMG5_boulevolp(mesh,k,ip,list);
        ilist = MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,2);

        if ( ilist > 0 ) {
          ier = MMG5_colver(mesh,met,list,ilist,iq,2);
          if ( ier < 0 ) return -1;
          else if ( ier ) {
            MMG3D_delPt(mesh,ier);
            break;
          }
        }
        else if ( ilist < 0 ) return -1;
      }
      if ( ier ) {
        p1->flag = base;
        if ( isnm )  nnm++;
        nc++;
        break;
      }
    }
  }

  return nc;
}

/** Check if moving mesh with disp for a fraction t yields a valid mesh */
int MMG5_chkmovmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t) {
  MMG5_pTetra  pt;
  MMG5_pPoint  ppt;
  double       *v,c[4][3],tau;
  int          k,np;
  char         i,j;

  /* Pseudo time-step = fraction of disp to perform */
  tau = (double)t / MMG3D_SHORTMAX;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<4; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      v = &disp->m[3*np];
      for (j=0; j<3; j++)
        c[i][j] = ppt->c[j]+tau*v[j];
    }

    if( MMG5_caltet_iso_4pt(c[0],c[1],c[2],c[3]) < MMG5_EPSOK) return 0;  //     Other criteria : eg. a rate of degradation, etc... ?
  }

  return 1;
}

/** Return the largest fraction t that makes the motion along disp valid */
short MMG5_dikomv(MMG5_pMesh mesh,MMG5_pSol disp) {
  int     it,maxit;
  short   t,tmin,tmax;
  char    ier;

  maxit = 200;
  it = 0;

  tmin = 0;
  tmax = MMG3D_SHORTMAX;

  /* If full displacement can be achieved */
  if ( MMG5_chkmovmesh(mesh,disp,tmax) )
    return tmax;

  /* Else, find the largest displacement by dichotomy */
  while( tmin != tmax && it < maxit ) {
    t = (tmin+tmax)/2;

    /* Case that tmax = tmin +1 : check move with tmax */
    if ( t == tmin ) {
      ier = MMG5_chkmovmesh(mesh,disp,tmax);
      if ( ier )
        return tmax;
      else
        return tmin;
    }

    /* General case: check move with t */
    ier = MMG5_chkmovmesh(mesh,disp,t);
    if ( ier )
      tmin = t;
    else
      tmax = t;

    it++;
  }

  return tmin;
}

/** Perform mesh motion along disp, for a fraction t, and the corresponding updates */
int MMG5_dispmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t,int itdeg) {
  MMG5_pTetra   pt;
  MMG5_pPoint   ppt;
  double        *v,tau,ctau,c[4][3],ocal,ncal;
  int           k,np;
  char          i,j;

  tau = (double)t /MMG3D_SHORTMAX;
  ctau = 1.0 - tau;

  /* Identify elements which are very distorted in the process */
  for (k=1; k<=mesh->ne; k++) {
    pt  = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<4; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      for (j=0; j<3; j++)
        c[i][j] = ppt->c[j];
    }

    ocal = MMG5_caltet_iso_4pt(c[0],c[1],c[2],c[3]);

    for (i=0; i<4; i++) {
      np = pt->v[i];
      v = &disp->m[3*np];
      for (j=0; j<3; j++)
        c[i][j] += tau*v[j];
    }

    ncal = MMG5_caltet_iso_4pt(c[0],c[1],c[2],c[3]);

    if ( ncal < MMG5_DEGTOL*ocal )
      pt->mark = itdeg;

  }

  /* Perform physical displacement */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];

    if ( !MG_VOK(ppt) ) continue;
    v = &disp->m[3*k];

    for (i=0; i<3; i++) {
      ppt->c[i] = ppt->c[i] + tau*v[i];
      v[i] *= ctau;
    }
  }

  return 1;
}

/** For debugging purposes: save disp */
int MMG5_saveDisp(MMG5_pMesh mesh,MMG5_pSol disp) {
  FILE        *out;
  int         k;
  char        data[256],*ptr;

  strcpy(data,disp->namein);
  ptr = strstr(data,".sol");
  *ptr = '\0';
  strcat(data,".o.disp.sol");

  out = fopen(data,"w");

  fprintf(out,"MeshVersionFormatted 1\n\nDimension\n%d\n\n",disp->dim);
  fprintf(out,"SolAtVertices\n%d\n 1 2\n",disp->np);

  /* Print solutions */
  for(k=1; k<= disp->np; k++) {
    fprintf(out,"%f %f %f\n",disp->m[3*k+0],disp->m[3*k+1],disp->m[3*k+2]);
  }

  fprintf(out,"\nEnd");
  fclose(out);

  return 1;
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
int MMG5_mmg3d3(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met) {
  double  avlen,tau;
  int     itdc,itmn,maxitmn,maxitdc,nspl,ns,nm,nc,iit,k,warn;
  int     nns,nnm,nnc,nnspl,nnns,nnnm,nnnc,nnnspl;
  short   t;
  char    ier;

  tau = 0.0;
  maxitmn = 10;
  maxitdc = 100;
  t  = 0;

  //++mesh->info.fem;

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** LAGRANGIAN MOTION\n");

  /* Field mark stores information about whether a tetra has been greatly deformed during current step */
  for (k=1; k<=mesh->ne; k++)
    mesh->tetra[k].mark = 0;

  /* Estimates of the minimum and maximum edge lengths in the mesh */
  avlen = MMG5_estavglen(mesh);
  mesh->info.hmax = MMG3D_LLONG*avlen;
  mesh->info.hmin = MMG3D_LOPTS*avlen;

  //printf("Average length: %f ; proceed with hmin = %f, hmax = %f\n",avlen,mesh->info.hmin,mesh->info.hmax);

  for (itmn=0; itmn<maxitmn; itmn++) {
    nnnspl = nnnc = nnns = nnnm = 0;

    /* Extension of the velocity field */
    if ( !MMG5_velextLS(mesh,disp) ) {
      fprintf(stderr,"\n  ## Problem in func. MMG5_packLS. Exit program.\n");
      return 0;
    }

    //MMG5_saveDisp(mesh,disp);

    /* Dichotomy loop */
    for (itdc=0; itdc<maxitdc; itdc++) {
      nnspl = nnc = nns = nnm = 0;

      t = MMG5_dikomv(mesh,disp);
      if ( t == 0 ) {
        if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
          fprintf(stderr,"\n   *** Stop: impossible to proceed further\n");
        break;
      }

      ier = MMG5_dispmesh(mesh,disp,t,itdc);
      if ( !ier ) {
        fprintf(stderr,"\n  ** Impossible motion\n");
        return 0;
      }

      tau = tau + ((double)t /MMG3D_SHORTMAX)*(1.0-tau);
      if ( (abs(mesh->info.imprim) > 3 ) || mesh->info.ddebug )
        printf("   ---> Realized displacement: %f\n",tau);

      /* Local remeshing depending on the option */
      if ( mesh->info.lag > 0 ) {
        for ( iit=0; iit<5; iit++) {

          nspl = nc = ns = nm = 0;

          if ( mesh->info.lag > 1 ) {
            /* Split of points */
            nspl = MMG5_spllag(mesh,disp,met,itdc,&warn);
            if ( nspl < 0 ) {
              fprintf(stderr,"\n  ## Problem in spllag. Exiting.\n");
              return 0;
            }

            /* Collapse of points */
            nc = MMG5_coltetlag(mesh,met,itdc);
            if ( nc < 0 ) {
              fprintf(stderr,"\n  ## Problem in coltetlag. Exiting.\n");
              return 0;
            }
          }

          /* Swap of edges in tetra that have resulted distorted from the process */
          /* I do not know whether it is safe to put NULL in metric here (a
           * priori ok, since there is no vertex creation or suppression) */
          ns = MMG5_swptetlag(mesh,met,MMG3D_LSWAPIMPROVE,NULL,itdc);
          if ( ns < 0 ) {
            fprintf(stderr,"\n  ## Problem in swaptetlag. Exiting.\n");
            return 0;
          }

          /* Relocate vertices of tetra which have been distorted in the displacement process */
          nm = MMG5_movtetlag(mesh,met,itdc);
          if ( nm < 0 ) {
            fprintf(stderr,"\n  ## Problem in movtetlag. Exiting.\n");
            return 0;
          }

          if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && (nspl+nc+ns+nm > 0) )
            printf(" %d edges splitted, %d vertices collapsed, %d elements"
                   " swapped, %d vertices moved.\n",nspl,nc,ns,nm);
          nnspl+= nspl;
          nnm  += nm;
          nnc  += nc;
          nns  += ns;
        }
        if ( abs(mesh->info.imprim) > 3 && abs(mesh->info.imprim) < 5
             && (nnspl+nnm+nns+nnc > 0) )
          printf(" %d edges splitted, %d vertices collapsed, %d elements"
                 " swapped, %d vertices moved.\n",nnspl,nnc,nns,nnm);
      }

      nnnspl += nnspl;
      nnnm   += nnm;
      nnnc   += nnc;
      nnns   += nns;

      if ( t == MMG3D_SHORTMAX ) break;
    }
    if ( mesh->info.imprim > 0 && abs(mesh->info.imprim) < 4 ) {
      printf("   ---> Realized displacement: %f\n",tau);
      if ( abs(mesh->info.imprim) > 2 )
        printf(" %d edges splitted, %d vertices collapsed, %d elements"
               " swapped, %d vertices moved.\n",nnnspl,nnnc,nnns,nnnm);
    }

    if ( t == MMG3D_SHORTMAX ) break;
  }
  /* Clean memory */
  /* Doing this, memcur of mesh is decreased by size of displacement */
  MMG5_DEL_MEM(mesh,disp->m);

  return 1;
}
#else
/**
 *
 * Hack to avoid to have an empty translation unit (forbidden by ISO C)
 *
 */
static void MMG3D_unused_function(void) {
  return;
}

#endif
