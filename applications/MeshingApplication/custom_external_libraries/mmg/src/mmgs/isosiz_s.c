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
 * \file mmgs/isosiz_s.c
 * \brief Fonctions for isotropic size map computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"
#include <math.h>

#define MAXLEN   1.0e+3

/* Define isotropic size map at all vertices of the mesh, associated with geometric approx ;
   by convention, p0->h stores desired length at point p0 */
int _MMGS_defsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  MMG5_pPoint   p[3];
  MMG5_pPar     par;
  double   n[3][3],t[3][3],nt[3],c1[3],c2[3],*n1,*n2,*t1,*t2;
  double   ps,ps2,ux,uy,uz,ll,l,lm,dd,M1,M2,hausd,hmin,hmax;
  int      k,j,ip1,ip2,isloc;
  char     i,i1,i2;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Defining isotropic map\n");

  if ( mesh->info.hmax < 0.0 )  mesh->info.hmax = 0.5 * mesh->info.delta;

  /* alloc structure */
  if ( !met->m ) {
    met->np    = mesh->np;
    met->npmax = mesh->npmax;
    met->size  = 1;
    _MMG5_ADD_MEM(mesh,(met->npmax+1)*sizeof(double),"solution",return(0));
    _MMG5_SAFE_MALLOC(met->m,mesh->npmax+1,double);

    /* init constant size */
    for (k=1; k<=mesh->np; k++)
      met->m[k] = mesh->info.hmax;
  }

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    p[0] = &mesh->point[pt->v[0]];
    p[1] = &mesh->point[pt->v[1]];
    p[2] = &mesh->point[pt->v[2]];

    /* normal recovery */
    for (i=0; i<3; i++) {
      if ( MS_SIN(p[i]->tag) ) {
        _MMG5_nortri(mesh,pt,n[i]);
      }
      else if ( MG_EDG(p[i]->tag) ) {
        _MMG5_nortri(mesh,pt,nt);
        n1  = &mesh->xpoint[p[i]->xp].n1[0];
        n2  = &mesh->xpoint[p[i]->xp].n2[0];
        ps  = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
        ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];
        if ( fabs(ps) > fabs(ps2) )
          memcpy(&n[i],n1,3*sizeof(double));
        else
          memcpy(&n[i],n2,3*sizeof(double));
        memcpy(&t[i],p[i]->n,3*sizeof(double));
      }
      else
        memcpy(&n[i],p[i]->n,3*sizeof(double));
    }

    for (i=0; i<3; i++) {
      i1  = _MMG5_inxt2[i];
      i2  = _MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];

      /* local parameters */
      hausd = mesh->info.hausd;
      hmin  = mesh->info.hmin;
      hmax  = mesh->info.hmax;
      isloc = 0;
      for (j=0; j<mesh->info.npar; j++) {
        par = &mesh->info.par[j];
        if ( /*( (par->elt == MMG5_Vertex) &&
               ( (p[i1]->ref == par->ref ) || (p[i2]->ref == par->ref) ))
               ||*/ ( (par->elt == MMG5_Triangle) && (pt->ref == par->ref ) ) ) {
          if ( !isloc ) {
            hausd = par->hausd;
            hmin  = par->hmin;
            hmax  = par->hmax;
            isloc = 1;
          }
          else {
            hausd = MG_MIN(par->hausd,hausd);
            hmin  = MG_MAX(par->hmin,hmin);
            hmax  = MG_MIN(par->hmax,hmax);
          }
        }
      }

      ux = p[i2]->c[0] - p[i1]->c[0];
      uy = p[i2]->c[1] - p[i1]->c[1];
      uz = p[i2]->c[2] - p[i1]->c[2];
      ll = ux*ux + uy*uy + uz*uz;

      if ( ll < _MMG5_EPSD )  continue;

      if ( MG_EDG(pt->tag[i]) ) {
        if ( MS_SIN(p[i1]->tag) ) {
          t[i1][0] = ux;
          t[i1][1] = uy;
          t[i1][2] = uz;
          l = 1.0 / sqrt(ll);
          t[i1][0] *= l;
          t[i1][1] *= l;
          t[i1][2] *= l;
        }
        if ( MS_SIN(p[i2]->tag) ) {
          t[i2][0] = -ux;
          t[i2][1] = -uy;
          t[i2][2] = -uz;
          l = 1.0/sqrt(ll);
          t[i2][0] *= l;
          t[i2][1] *= l;
          t[i2][2] *= l;
        }
        t1 = t[i1];
        t2 = t[i2];

        /* The two Bezier coefficients along curve */
        dd    = (t1[0]*ux + t1[1]*uy + t1[2]*uz)/3.0;
        c1[0] = p[i1]->c[0] + dd * t1[0];
        c1[1] = p[i1]->c[1] + dd * t1[1];
        c1[2] = p[i1]->c[2] + dd * t1[2];

        dd    = -(t2[0]*ux + t2[1]*uy + t2[2]*uz)/3.0;
        c2[0] = p[i2]->c[0] + dd * t2[0];
        c2[1] = p[i2]->c[1] + dd * t2[1];
        c2[2] = p[i2]->c[2] + dd * t2[2];

        M1 = (c2[0]-2.0*c1[0]+p[i1]->c[0])*(c2[0]-2.0*c1[0]+p[i1]->c[0]) \
          + (c2[1]-2.0*c1[1]+p[i1]->c[1])*(c2[1]-2.0*c1[1]+p[i1]->c[1]) \
          + (c2[2]-2.0*c1[2]+p[i1]->c[2])*(c2[2]-2.0*c1[2]+p[i1]->c[2]);

        M2 = (p[i2]->c[0]-2.0*c2[0]+c1[0])*(p[i2]->c[0]-2.0*c2[0]+c1[0]) \
          + (p[i2]->c[1]-2.0*c2[1]+c1[1])*(p[i2]->c[1]-2.0*c2[1]+c1[1])\
          + (p[i2]->c[2]-2.0*c2[2]+c1[2])*(p[i2]->c[2]-2.0*c2[2]+c1[2]);

        M1 = 6.0 * sqrt(M1);
        M2 = 6.0 * sqrt(M2);
        M1 = MG_MAX(M1,M2);

        if ( M1 < _MMG5_EPSD )
          lm = MAXLEN;
        else {
          lm = (16.0*ll*hausd) / (3.0*M1);
          lm = sqrt(lm);
        }
        met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
        met->m[ip2] = MG_MAX(hmin,MG_MIN(met->m[ip2],lm));
      }
      else {
        n1 = n[i1];
        n2 = n[i2];

        ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
        c1[0] = (2.0*p[i1]->c[0] + p[i2]->c[0] - ps*n1[0]) / 3.0;
        c1[1] = (2.0*p[i1]->c[1] + p[i2]->c[1] - ps*n1[1]) / 3.0;
        c1[2] = (2.0*p[i1]->c[2] + p[i2]->c[2] - ps*n1[2]) / 3.0;

        ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
        c2[0] = (2.0*p[i2]->c[0] + p[i1]->c[0] - ps*n2[0]) / 3.0;
        c2[1] = (2.0*p[i2]->c[1] + p[i1]->c[1] - ps*n2[1]) / 3.0;
        c2[2] = (2.0*p[i2]->c[2] + p[i1]->c[2] - ps*n2[2]) / 3.0;

        M1 = (c2[0]-2.0*c1[0]+p[i1]->c[0])*(c2[0]-2.0*c1[0]+p[i1]->c[0]) \
          + (c2[1]-2.0*c1[1]+p[i1]->c[1])*(c2[1]-2.0*c1[1]+p[i1]->c[1]) \
          + (c2[2]-2.0*c1[2]+p[i1]->c[2])*(c2[2]-2.0*c1[2]+p[i1]->c[2]);

        M2 = (p[i2]->c[0]-2.0*c2[0]+c1[0])*(p[i2]->c[0]-2.0*c2[0]+c1[0]) \
          + (p[i2]->c[1]-2.0*c2[1]+c1[1])*(p[i2]->c[1]-2.0*c2[1]+c1[1])\
          + (p[i2]->c[2]-2.0*c2[2]+c1[2])*(p[i2]->c[2]-2.0*c2[2]+c1[2]);

        M1 = 6.0 * sqrt(M1);
        M2 = 6.0 * sqrt(M2);
        M1 = MG_MAX(M1,M2);

        if ( M1 < _MMG5_EPSD )
          lm = MAXLEN;
        else {
          lm = (16.0*ll*hausd) / (3.0*M1);
          lm = sqrt(lm);
        }
        met->m[ip1] = MG_MAX(hmin,MG_MIN(met->m[ip1],lm));
        met->m[ip2] = MG_MAX(hmin,MG_MIN(met->m[ip2],lm));
      }
    }
  }

  /* take local parameters */
  for (j=0; j<mesh->info.npar; j++) {
    par = &mesh->info.par[j];
    /* if ( par->elt == MMG5_Vertex ) { */
    /*   for (k=1; k<=mesh->np; k++) { */
    /*     ppt = &mesh->point[k]; */
    /*     if ( !MG_VOK(ppt) || ppt->ref != par->ref )  continue; */
    /*     met->m[k] = MG_MAX(par->hmin,MG_MIN(met->m[k],par->hmax)); */
    /*   } */
    /* } */
    /* else */ if ( par->elt == MMG5_Triangle ) {
      for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MG_EOK(pt) || pt->ref != par->ref )  continue;
        met->m[pt->v[0]] = MG_MAX(par->hmin,MG_MIN(met->m[pt->v[0]],par->hmax));
        met->m[pt->v[1]] = MG_MAX(par->hmin,MG_MIN(met->m[pt->v[1]],par->hmax));
        met->m[pt->v[2]] = MG_MAX(par->hmin,MG_MIN(met->m[pt->v[2]],par->hmax));
      }
    }
  }
  return(1);
}


/** Enforces mesh gradations by truncating size map */
int gradsiz_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  MMG5_pPoint   p1,p2;
  double   ll,hn,h1,h2;
  int      k,nu,nup,it,maxit,ip1,ip2;
  char     i,i1,i2;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** Grading mesh\n");

  mesh->base = 0;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = mesh->base;

  it = nup = 0;
  maxit = 100;
  do {
    mesh->base++;
    nu = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<3; i++) {
        i1  = _MMG5_inxt2[i];
        i2  = _MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        p1 = &mesh->point[ip1];
        p2 = &mesh->point[ip2];
        if ( p1->flag < mesh->base-1 && p2->flag < mesh->base-1 )  continue;

        ll = (p2->c[0]-p1->c[0])*(p2->c[0]-p1->c[0]) + (p2->c[1]-p1->c[1])*(p2->c[1]-p1->c[1]) \
          + (p2->c[2]-p1->c[2])*(p2->c[2]-p1->c[2]);
        ll = sqrt(ll);

        h1 = met->m[ip1];
        h2 = met->m[ip2];
        if ( h1 < h2 ) {
          if ( h1 < _MMG5_EPSD )  continue;
          hn  = h1 + mesh->info.hgrad*ll;
          if ( h2 > hn ) {
            met->m[ip2] = hn;
            p2->flag    = mesh->base;
            nu++;
          }
        }
        else {
          if ( h2 < _MMG5_EPSD )  continue;
          hn = h2 + mesh->info.hgrad*ll;
          if ( h1 > hn ) {
            met->m[ip1] = hn;
            p1->flag    = mesh->base;
            nu++;
          }
        }
      }
    }
    nup += nu;
  }
  while ( ++it < maxit && nu > 0 );

  if ( abs(mesh->info.imprim) > 4 )  fprintf(stdout,"     gradation: %7d updated, %d iter.\n",nup,it);
  return(1);
}

