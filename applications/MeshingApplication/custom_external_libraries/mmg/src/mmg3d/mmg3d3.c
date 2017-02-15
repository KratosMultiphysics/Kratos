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
#define _MMG5_DEGTOL  1.e-1

extern char  ddb;

/** Calculate an estimate of the average (isotropic) length of edges in the mesh */
double _MMG5_estavglen(MMG5_pMesh mesh) {
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
      i1 = _MMG5_iare[i][0];
      i2 = _MMG5_iare[i][1];
      
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
  
  return(lent);
}

/** Interpolate displacement between v1 and v2 at intermediate position 0<=t<=1 */
static
inline int _MMG5_intdispvol(double *v1, double *v2, double *vp, double t) {
  char i;
  
  for(i=0; i<3; i++)
    vp[i] = (1.0-t)*v1[i] + t*v2[i];
  
  return(1);
}

/** compute quality iso of a tetra given by the 4 points a,b,c,d */
static
inline double _MMG5_caltet_iso_4pt(double *a, double *b, double *c, double *d) {
  double     abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz;
  double     vol,v1,v2,v3,rap;
  
  /* volume */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];
  rap = abx*abx + aby*aby + abz*abz;
  
  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];
  rap += acx*acx + acy*acy + acz*acz;
  
  adx = d[0] - a[0];
  ady = d[1] - a[1];
  adz = d[2] - a[2];
  rap += adx*adx + ady*ady + adz*adz;
  
  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;
  if ( vol < _MMG5_EPSD2 )  return(0.0);
  
  bcx = c[0] - b[0];
  bcy = c[1] - b[1];
  bcz = c[2] - b[2];
  rap += bcx*bcx + bcy*bcy + bcz*bcz;
  
  bdx = d[0] - b[0];
  bdy = d[1] - b[1];
  bdz = d[2] - b[2];
  rap += bdx*bdx + bdy*bdy + bdz*bdz;
  
  cdx = d[0] - c[0];
  cdy = d[1] - c[1];
  cdz = d[2] - c[2];
  rap += cdx*cdx + cdy*cdy + cdz*cdz;
  if ( rap < _MMG5_EPSD2 )  return(0.0);
  
  /* quality = vol / len^3/2 */
  rap = rap * sqrt(rap);
  return(vol / rap);
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
 * Split edges of length bigger than _MMG5_LOPTL, in the Lagrangian mode.
 *
 */
static int _MMG5_spllag(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met,int itdeg, int* warn) {
  MMG5_pTetra     pt;
  MMG5_pPoint     p0,p1,ppt;
  double     dd,len,lmax,o[3],hma2;
  double    *m1,*m2,*mp;
  int        k,ip,ip1,ip2,list[MMG3D_LMAX+2],ilist,ns,ier,iadr;
  char       imax,i,i1,i2;
  
  *warn=0;
  ns = 0;
  hma2 = mesh->info.hmax*mesh->info.hmax;
  
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;
    if ( pt->mark != itdeg ) continue;
    
    /* find longest, internal edge */
    imax = -1;
    lmax = 0.0;
    for (i=0; i<6; i++) {
      i1  = _MMG5_iare[i][0];
      i2  = _MMG5_iare[i][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      p0 = &mesh->point[ip1];
      p1 = &mesh->point[ip2];

      if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) )  continue;
      
      len = (p1->c[0]-p0->c[0])*(p1->c[0]-p0->c[0])
        + (p1->c[1]-p0->c[1])*(p1->c[1]-p0->c[1])
        + (p1->c[2]-p0->c[2])*(p1->c[2]-p0->c[2]);
      
      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }
    if ( imax==-1 )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are required or of length null.\n",__FILE__,__LINE__,k);
    
    if ( lmax < hma2 )  continue;
    
    /* proceed edges according to lengths */
    i1  = _MMG5_iare[imax][0];
    i2  = _MMG5_iare[imax][1];
    ip1 = pt->v[i1];
    ip2 = pt->v[i2];
    p0 = &mesh->point[ip1];
    p1 = &mesh->point[ip2];
    
    /* Deal only with internal faces */
    assert( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) );
    ilist = _MMG5_coquil(mesh,k,imax,list);
    if ( !ilist ) continue;
    else if ( ilist<0 ) return(-1);
    o[0] = 0.5*(p0->c[0] + p1->c[0]);
    o[1] = 0.5*(p0->c[1] + p1->c[1]);
    o[2] = 0.5*(p0->c[2] + p1->c[2]);
      
    ip = _MMG3D_newPt(mesh,o,MG_NOTAG);
      
    if ( !ip )  {
      /* reallocation of point table */
      _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,*warn=1;break,o,MG_NOTAG);
    }
    ppt = &mesh->point[ip];
    
    /* Interpolation of metric, if any */
    if ( met->m ) {
      if ( !_MMG5_intmet(mesh,met,k,imax,ip,0.5) ) {
        _MMG3D_delPt(mesh,ip);
        return(-1);
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
      
      if ( !_MMG5_intdispvol(m1,m2,mp,0.5) ) {
        _MMG3D_delPt(mesh,ip);
        return(-1);
      }
    }
    
    /* Il y a un check sur la taille des arêtes ici aussi ! */
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
  
  return(ns);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param crit coefficient of quality improvment.
 * \param octree pointer toward the octree structure in delaunay mode and
 * toward the \a NULL pointer otherwise.
 * \param itdeg degraded elements.
 *
 * Internal edge flipping in the Lagrangian mode; only affects tetra marked with it
 *
 */
int _MMG5_swptetlag(MMG5_pMesh mesh,MMG5_pSol met,double crit,_MMG3D_pOctree octree,int itdeg) {
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

        nconf = _MMG5_chkswpgen(mesh,met,k,i,&ilist,list,crit,2);
        
        if ( nconf ) {
          ier = _MMG5_swpgen(mesh,met,nconf,ilist,list,octree,2);
          if ( ier > 0 )  ns++;
          else if ( ier < 0 ) return(-1);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  return(nns);
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
int _MMG5_movtetlag(MMG5_pMesh mesh,MMG5_pSol met, int itdeg) {
  MMG5_pTetra        pt;
  MMG5_pPoint        ppt;
  int           k,ier,nm,nnm,ns,listv[MMG3D_LMAX+2],ilistv,it;
  unsigned char i,base;
  int           maxit;
  
  base = 1;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;
  
  maxit = 2;
  it = nnm = 0;
  
  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
      if ( pt->mark != itdeg ) continue;
      
      /* point i */
      for (i=0; i<4; i++) {
        ppt = &mesh->point[pt->v[i]];
        if ( ppt->flag == base )  continue;
        else if ( ppt->tag & MG_BDY ) continue;

        ier = 0;
 
        ilistv = _MMG5_boulevolp(mesh,k,i,listv);
        if ( !ilistv )  continue;
        
        ier = _MMG5_movintpt_iso(mesh,met, NULL, listv,ilistv,0);
          
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
  
  return(nnm);
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
static int _MMG5_coltetlag(MMG5_pMesh mesh,MMG5_pSol met,int itdeg) {
  MMG5_pTetra     pt;
  MMG5_pPoint     p0,p1;
  double     ll,ux,uy,uz,hmi2;
  int        k,nc,list[MMG3D_LMAX+2],ilist,base,nnm;
  int16_t    tag;
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
        ip = _MMG5_idir[i][_MMG5_inxt2[j]];
        iq = _MMG5_idir[i][_MMG5_iprv2[j]];
        
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
        ilist = _MMG5_boulevolp(mesh,k,ip,list);
        ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,2);
      
        if ( ilist > 0 ) {
          ier = _MMG5_colver(mesh,met,list,ilist,iq,2);
          if ( ier < 0 ) return(-1);
          else if ( ier ) {
            _MMG3D_delPt(mesh,ier);
            break;
          }
        }
        else if (ilist < 0 ) return(-1);
      }
      if ( ier ) {
        p1->flag = base;
        if ( isnm )  nnm++;
        nc++;
        break;
      }
    }
  }
  
  return(nc);
}

/** Check if moving mesh with disp for a fraction t yields a valid mesh */
int _MMG5_chkmovmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t) {
  MMG5_pTetra  pt;
  MMG5_pPoint  ppt;
  double       *v,c[4][3],tau,cal;
  int          k,np;
  char         i,j;
  
  /* Pseudo time-step = fraction of disp to perform */
  tau = (double)t / _MMG5_SHORTMAX;
  
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
    
    if( _MMG5_caltet_iso_4pt(c[0],c[1],c[2],c[3]) < _MMG5_NULKAL) return(0);  //     Other criteria : eg. a rate of degradation, etc... ?
  }

  return(1);
}

/** Return the largest fraction t that makes the motion along disp valid */
short _MMG5_dikomv(MMG5_pMesh mesh,MMG5_pSol disp) {
  int     it,maxit;
  short   t,tmin,tmax;
  char    ier;
 
  maxit = 200;
  it = 0;
  
  tmin = 0;
  tmax = _MMG5_SHORTMAX;
  
  /* If full displacement can be achieved */
  if ( _MMG5_chkmovmesh(mesh,disp,tmax) )
    return(tmax);

  /* Else, find the largest displacement by dichotomy */
  while( tmin != tmax && it < maxit ) {
    t = (tmin+tmax)/2;
  
    /* Case that tmax = tmin +1 : check move with tmax */
    if ( t == tmin ) {
      ier = _MMG5_chkmovmesh(mesh,disp,tmax);
      if ( ier )
        return(tmax);
      else
        return(tmin);
    }
    
    /* General case: check move with t */
    ier = _MMG5_chkmovmesh(mesh,disp,t);
    if ( ier )
      tmin = t;
    else
      tmax = t;
  
    it++;
  }
  
  return(tmin);
}

/** Perform mesh motion along disp, for a fraction t, and the corresponding updates */
int _MMG5_dispmesh(MMG5_pMesh mesh,MMG5_pSol disp,short t,int itdeg) {
  MMG5_pTetra   pt;
  MMG5_pPoint   ppt;
  double        *v,tau,ctau,c[4][3],ocal,ncal;
  int           k,np;
  char          i,j;
  
  tau = (double)t /_MMG5_SHORTMAX;
  ctau = 1.0 - tau;
  
  /* Identify elements which are very distorted in the process */
  for (k=1; k<=mesh->ne; k++) {
    pt  = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    
    for (i=0; i<4; i++) {
      np = pt->v[i];
      ppt = &mesh->point[np];
      v = &disp->m[3*np];
      for (j=0; j<3; j++)
        c[i][j] = ppt->c[j];
    }
    
    ocal = _MMG5_caltet_iso_4pt(c[0],c[1],c[2],c[3]);
    
    for (i=0; i<4; i++) {
      np = pt->v[i];
      v = &disp->m[3*np];
      for (j=0; j<3; j++)
        c[i][j] += tau*v[j];
    }
    
    ncal = _MMG5_caltet_iso_4pt(c[0],c[1],c[2],c[3]);
    
    if ( ncal < _MMG5_DEGTOL*ocal )
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
  
  return(1);
}

/** For debugging purposes: save disp */
int _MMG5_saveDisp(MMG5_pMesh mesh,MMG5_pSol disp) {
  FILE        *out;
  int         k;
  char        j,data[256],*ptr;
  
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
  
  return(1);
}

/** Lagrangian node displacement and meshing */
int _MMG5_mmg3d3(MMG5_pMesh mesh,MMG5_pSol disp,MMG5_pSol met) {
  double  avlen,tau;
  int     itdc,itmn,maxitmn,maxitdc,nspl,ns,nm,nc,iit,k,warn;
  int     nns,nnm,nnc,nnspl,nnns,nnnm,nnnc,nnnspl;
  short   t;
  char    ier;
  
  tau = 0.0;
  maxitmn = 10;
  maxitdc = 100;
  t  = 0;
  
  //mesh->info.fem = 1;
  
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** LAGRANGIAN MOTION\n");
  
  /* Field mark stores information about whether a tetra has been greatly deformed during current step */
  for (k=1; k<=mesh->ne; k++)
    mesh->tetra[k].mark = 0;

  /* Estimates of the minimum and maximum edge lengths in the mesh */
  avlen = _MMG5_estavglen(mesh);
  mesh->info.hmax = _MMG5_LLONG*avlen;
  mesh->info.hmin = _MMG5_LOPTS*avlen;
  
  //printf("Average length: %f ; proceed with hmin = %f, hmax = %f\n",avlen,mesh->info.hmin,mesh->info.hmax);

  for (itmn=0; itmn<maxitmn; itmn++) {
    nnnspl = nnnc = nnns = nnnm = 0;

    /* Extension of the velocity field */
    if ( !_MMG5_velextLS(mesh,disp) ) {
      fprintf(stderr,"  ## Problem in func. _MMG5_packLS. Exit program.\n");
      return(0);
    }
  
    //_MMG5_saveDisp(mesh,disp);
  
    /* Dichotomy loop */
    for (itdc=0; itdc<maxitdc; itdc++) {
      nnspl = nnc = nns = nnm = 0;

      t = _MMG5_dikomv(mesh,disp);
      if ( t == 0 ) {
        if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
          printf("   *** Stop: impossible to proceed further\n");
        break;
      }
  
      ier = _MMG5_dispmesh(mesh,disp,t,itdc);
      if ( !ier ) {
        fprintf(stderr,"  ** Impossible motion\n");
        return(0);
      }
    
      tau = tau + ((double)t /_MMG5_SHORTMAX)*(1.0-tau);
      if ( (abs(mesh->info.imprim) > 3 ) || mesh->info.ddebug )
        printf("   ---> Realized displacement: %f\n",tau);
    
      /* Local remeshing depending on the option */
      if ( mesh->info.lag > 0 ) {
        for ( iit=0; iit<5; iit++) {
  
          nspl = nc = ns = nm = 0;
    
          if ( mesh->info.lag > 1 ) {
            /* Split of points */
            nspl = _MMG5_spllag(mesh,disp,met,itdc,&warn);
            if ( nspl < 0 ) {
              fprintf(stderr,"  ## Problem in spllag. Exiting.\n");
              return(0);
            }
  
            /* Collapse of points */
            nc = _MMG5_coltetlag(mesh,met,itdc);
            if ( nc < 0 ) {
              fprintf(stderr,"  ## Problem in coltetlag. Exiting.\n");
              return(0);
            }
          }
        
          /* Swap of edges in tetra that have resulted distorted from the process */
          /* I do not know whether it is safe to put NULL in metric here (a
           * priori ok, since there is no vertex creation or suppression) */
          ns = _MMG5_swptetlag(mesh,met,1.1,NULL,itdc);
          if ( ns < 0 ) {
            fprintf(stderr,"  ## Problem in swaptetlag. Exiting.\n");
            return(0);
          }
      
          /* Relocate vertices of tetra which have been distorted in the displacement process */
          nm = _MMG5_movtetlag(mesh,met,itdc);
          if ( nm < 0 ) {
            fprintf(stderr,"  ## Problem in movtetlag. Exiting.\n");
            return(0);
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

      if ( t == _MMG5_SHORTMAX ) break;
    }
    if ( mesh->info.imprim && abs(mesh->info.imprim) < 4 ) {
      printf("   ---> Realized displacement: %f\n",tau);
      if ( abs(mesh->info.imprim) > 2 )
        printf(" %d edges splitted, %d vertices collapsed, %d elements"
               " swapped, %d vertices moved.\n",nnnspl,nnnc,nnns,nnnm);
    }
    
    if ( t == _MMG5_SHORTMAX ) break;
  }
  /* Clean memory */
  /* Doing this, memcur of mesh is decreased by size of displacement */
  _MMG5_DEL_MEM(mesh,disp->m,(disp->size*(disp->npmax+1))*sizeof(double));

  return(1);
}

#endif
