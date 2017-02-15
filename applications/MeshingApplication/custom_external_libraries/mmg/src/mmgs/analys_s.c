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
 * \file mmgs/analys_s.c
 * \brief Mesh analysis.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"

/* topology: set adjacent, detect Moebius, flip faces, count connected comp. */
static int setadj(MMG5_pMesh mesh){
  MMG5_pTria   pt,pt1;
  int          *adja,*adjb,adji1,adji2,*pile,iad,ipil,ip1,ip2,gen;
  int          k,kk,iel,jel,nvf,nf,nr,nt,nre,ncc,ned,ref;
  int16_t      tag;
  char         i,ii,i1,i2,ii1,ii2,voy;

  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  _MMG5_SAFE_MALLOC(pile,mesh->nt+1,int);

  pile[1]  = 1;
  ipil     = 1;
  nvf = nre = nr = nf = nt = ncc = ned = 0;

  while ( ipil > 0 ) {
    ncc++;

    do {
      k  = pile[ipil--];
      pt = &mesh->tria[k];
      pt->flag = 0;
      if ( !MG_EOK(pt) )  continue;

      pt->cc = ncc;
      adja = &mesh->adja[3*(k-1)+1];
      nt++;
      for (i=0; i<3; i++) {
        i1  = _MMG5_inxt2[i];
        i2  = _MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];

        if ( !mesh->point[ip1].tmp )  mesh->point[ip1].tmp = ++nvf;
        if ( !mesh->point[ip2].tmp )  mesh->point[ip2].tmp = ++nvf;

        if ( MG_EDG(pt->tag[i]) ) {
          mesh->point[ip1].tag |= pt->tag[i];
          mesh->point[ip2].tag |= pt->tag[i];
        }

        /* open boundary */
        if ( !adja[i] ) {
          pt->tag[i] |= MG_GEO;
          mesh->point[ip1].tag |= MG_GEO;
          mesh->point[ip2].tag |= MG_GEO;
          nr++;
          ned++;
          continue;
        }

        kk = adja[i] / 3;
        ii = adja[i] % 3;
        if ( kk > k )  ned++;

        /* correct edge tag */
        pt1 = &mesh->tria[kk];
        if ( pt->tag[i] & MG_NOM && !(pt1->tag[ii] & MG_NOM) ) {
          pt1->tag[ii] = pt->tag[i];
          pt1->edg[ii] = pt->edg[i];
          mesh->point[ip1].tag |= MG_NOM;
          mesh->point[ip2].tag |= MG_NOM;
        }
        if ( pt1->tag[ii] & MG_NOM && !(pt->tag[i] & MG_NOM) ) {
          pt->tag[i] = pt1->tag[ii];
          pt->edg[i] = pt1->edg[ii];
          mesh->point[ip1].tag |= MG_NOM;
          mesh->point[ip2].tag |= MG_NOM;
        }
        if ( pt1->cc > 0 )  continue;

        if ( abs(pt1->ref) != abs(pt->ref) ) {
          pt->tag[i]   |= MG_REF;
          pt1->tag[ii] |= MG_REF;
          mesh->point[ip1].tag |= MG_REF;
          mesh->point[ip2].tag |= MG_REF;
          nre++;
        }

        /* store adjacent */
        if ( !pt1->flag ) {
          pt1->flag    = 1;
          if ( !(pt1->tag[ii] & MG_NOM) ) {
            pile[++ipil] = kk;
          }
        }

        /* check orientation */
        ii1 = _MMG5_inxt2[ii];
        ii2 = _MMG5_iprv2[ii];
        if ( pt1->v[ii1] == ip1 ) {
          /* Moebius strip */
          if ( pt1->base < 0 ) {
            pt1->ref      = -abs(pt1->ref);
            pt->tag[i]   |= MG_REF;
            pt1->tag[ii] |= MG_REF;
            nre++;
          }
          /* flip orientation */
          else if ( !(pt->tag[i] & MG_NOM) ) {
            pt1->base   = -pt1->base;
            pt1->v[ii1] = ip2;
            pt1->v[ii2] = ip1;

            /* update adj */
            iad   = 3*(kk-1)+1;
            adjb  = &mesh->adja[iad];
            adji1 = mesh->adja[iad+ii1];
            adji2 = mesh->adja[iad+ii2];
            adjb[ii1] = adji2;
            adjb[ii2] = adji1;

            /* modif tag + ref */
            tag = pt1->tag[ii1];
            pt1->tag[ii1] = pt1->tag[ii2];
            pt1->tag[ii2] = tag;
            ref = pt1->edg[ii1];
            pt1->edg[ii1] = pt1->edg[ii2];
            pt1->edg[ii2] = ref;

            /* modif voyeurs */
            if ( adjb[ii1] ) {
              iel = adjb[ii1] / 3;
              voy = adjb[ii1] % 3;
              mesh->adja[3*(iel-1)+1+voy] = 3*kk + ii1;
            }
            if ( adjb[ii2] ) {
              iel = adjb[ii2] / 3;
              voy = adjb[ii2] % 3;
              mesh->adja[3*(iel-1)+1+voy] = 3*kk + ii2;
            }
            nf++;
          }
        }
      }
    }
    while ( ipil > 0 );

    /* find next triangle */
    ipil = 0;
    for (kk=1; kk<=mesh->nt; kk++) {
      pt = &mesh->tria[kk];
      if ( pt->v[0] && (pt->cc == 0) ) {
        ipil = 1;
        pile[ipil] = kk;
        pt->flag   = 1;
        break;
      }
    }
  }

  /* bilan */
  nr = nre = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      if ( !MG_EDG(pt->tag[i]) )  continue;

      adja = &mesh->adja[3*(k-1)+1];
      jel  = adja[i] / 3;
      if ( !jel || jel > k ) {
        if ( pt->tag[i] & MG_GEO )  nr++;
        if ( pt->tag[i] & MG_REF )  nre++;
      }
    }
  }

  if ( mesh->info.ddebug ) {
    fprintf(stdout,"  a- ridges: %d found.\n",nr);
    fprintf(stdout,"  a- connex: %d connected component(s)\n",ncc);
    fprintf(stdout,"  a- orient: %d flipped\n",nf);
  }
  else if ( abs(mesh->info.imprim) > 4 ) {
    gen = (2 - nvf + ned - nt) / 2;
    fprintf(stdout,"     Connected component: %d,  genus: %d,   reoriented: %d\n",ncc,gen,nf);
    fprintf(stdout,"     Edges: %d,  tagged: %d,  ridges: %d,  refs: %d\n",ned,nr+nre,nr,nre);
  }

  _MMG5_SAFE_FREE(pile);
  return(1);
}

/* Detect non manifold points */
static void nmpoints(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0;
  int        k,np,numt,iel,jel,nmp,*adja;
  char       i0,i1,i,jp;
  
  nmp = 0;
  /* Initialize point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].s = 0;
  
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue; 

    for (i=0; i<3; i++) {
      np = pt->v[i];
      p0 = &mesh->point[np];
      if ( !p0->s )  p0->s = k;
    }  
  }
          
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue; 
    
    for (i=0; i<3; i++) {
      np = pt->v[i];
      p0 = &mesh->point[np];
      numt = p0->s;
      if ( k == numt )  continue;
      jel = k;
      jp  = i;
            
      /* Unfold ball of np, until numt is reached */
      do {
        iel = jel;
        i0  =  jp;
        i1  = _MMG5_inxt2[i0];
        adja = &mesh->adja[3*(iel-1)+1];
        jel = adja[i1] / 3;
        jp  = adja[i1] % 3;
        jp  = _MMG5_inxt2[jp];
      }
      while ( jel && (jel != numt) && (jel !=k) );

      /* Ball has been completely traveled without meeting numt */
      if ( jel == k ) {
        if ( !(p0->tag & MG_CRN) || !(p0->tag & MG_REQ) ) {
          nmp++;
          // p0->tag |= MG_CRN + MG_REQ;
        }
        continue;
      }
      else if ( jel == numt ) 
        continue;
  
      jel = iel;
      jp = i0;
              
      /* At this point, jel =0, i.e. an open boundary has been hit : travel in the opposite sense */
      do {
        iel = jel;
        i0  =  jp;
        i1  = _MMG5_iprv2[i0];
        adja = &mesh->adja[3*(iel-1)+1];
        jel = adja[i1] / 3;
        jp  = adja[i1] % 3;
        jp  = _MMG5_iprv2[jp];
      }
      while ( jel && (jel != numt));  
          
      if ( jel != numt) {
        if ( !(p0->tag & MG_CRN) || !(p0->tag & MG_REQ) ) {
          p0->tag |= MG_CRN + MG_REQ; 
          nmp++;
        }
      }
    }
  }  

  /* reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].s = 0;

  if ( nmp && abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ## %d non manifold points detected\n",nmp);
}

/** improve badly shaped elts for isotropic mesh */
/* static int delbad(MMG5_pMesh mesh) { */
/*   MMG5_pTria    pt; */
/*   MMG5_pPoint   p[3]; */
/*   double   s,kal,declic,ux,uy,uz,vx,vy,vz; */
/*   int     *adja,k,iel,nd,ndd,it; */
/*   char     i,ia,i1,i2,j,typ; */

/*   it = ndd = 0; */
/*   declic = BADKAL / ALPHAD; */

/*   do { */
/*     nd = 0; */
/*     for (k=1; k<=mesh->nt; k++) { */
/*       pt = &mesh->tria[k]; */
/*       if ( !MG_EOK(pt) )  continue; */

/*       kal = _MMG5_calelt(mesh,NULL,pt); */
/*       if ( kal > declic )  continue; */

/*       p[0] = &mesh->point[pt->v[0]]; */
/*       p[1] = &mesh->point[pt->v[1]]; */
/*       p[2] = &mesh->point[pt->v[2]]; */
/*       adja = &mesh->adja[3*(k-1)+1]; */
/*       typ  = typelt(p,&ia); */
      
/*       /\* needle *\/ */
/*       if ( typ == 1 ) { */
/*         if ( litcol(mesh,k,ia,kal) ) { */
/*           nd++; */
/*           continue; */
/*         } */
/*       } */
/*       /\* obtuse *\/ */
/*       else if ( typ == 2 ) { */
/*         /\* delete boundary elt *\/ */
/*         if ( !adja[ia] ) { */
/*           /\* update point coordinates on ridge *\/ */
/*           i1 = _MMG5_inxt2[ia]; */
/*           i2 = _MMG5_iprv2[ia]; */
/*           p[0] = &mesh->point[pt->v[ia]]; */
/*           p[1] = &mesh->point[pt->v[i1]]; */
/*           p[2] = &mesh->point[pt->v[i2]]; */
/*           ux = p[2]->c[0] - p[1]->c[0]; */
/*           uy = p[2]->c[1] - p[1]->c[1]; */
/*           uz = p[2]->c[2] - p[1]->c[2]; */
/*           vx = p[0]->c[0] - p[1]->c[0]; */
/*           vy = p[0]->c[1] - p[1]->c[1]; */
/*           vz = p[0]->c[2] - p[1]->c[2]; */
/*           s  = (ux*vx + uy*vy + uz*vz) / sqrt(ux*ux + uy*uy + uz*uz); */
/*           p[0]->c[0] = vx - s*ux; */
/*           p[0]->c[1] = vy - s*uy; */
/*           p[0]->c[2] = vz - s*uz; */
          
/*           delElt(mesh,k); */
/*           nd++; */
/*           continue; */
/*         } */
/*         if ( litswp(mesh,k,ia,kal) || litcol(mesh,k,ia,kal) ) { */
/*           nd++; */
/*           continue; */
/*         } */
/*       } */

/*       /\* brute force to improve *\/ */
/*       for (i=0; i<3; i++) { */
/*         if ( litswp(mesh,k,i,kal) || litcol(mesh,k,i,kal) ) { */
/*           nd++; */
/*           break; */
/*         } */
/*         else if ( adja[i] ) { */
/*           iel = adja[i] / 3; */
/*           j   = adja[i] % 3; */
/*           if ( litcol(mesh,iel,j,kal) ) { */
/*             nd++; */
/*             break; */
/*           } */
/*         } */
/*       } */
/*     } */
/*     ndd += nd; */
/*     if ( nd && (mesh->info.ddebug || mesh->info.imprim < 0) )  fprintf(stdout,"     %d improved\n",nd); */
/*   } */
/*   while ( nd > 0 && ++it < 5 ); */

/*   if ( abs(mesh->info.imprim) > 4 ) */
/*     fprintf(stdout,"     %d bad elements improved\n",ndd); */

/*   return(1); */
/* } */


/* check for ridges: dihedral angle */
static int setdhd(MMG5_pMesh mesh) {
  MMG5_pTria    pt,pt1;
  double   n1[3],n2[3],dhd;
  int     *adja,k,kk,nr;
  char     i,ii,i1,i2;

  nr = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    /* triangle normal */
    _MMG5_nortri(mesh,pt,n1);
    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_GEO )  continue;
      kk = adja[i] / 3;
      ii = adja[i] % 3;

      /* check angle w. neighbor */
      if ( k < kk ) {
        pt1 = &mesh->tria[kk];
        _MMG5_nortri(mesh,pt1,n2);
        dhd = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
        if ( dhd <= mesh->info.dhd ) {
          pt->tag[i]   |= MG_GEO;
          pt1->tag[ii] |= MG_GEO;
          i1 = _MMG5_inxt2[i];
          i2 = _MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_GEO;
          mesh->point[pt->v[i2]].tag |= MG_GEO;
          nr++;
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 4 && nr > 0 )
    fprintf(stdout,"     %d ridges updated\n",nr);

  return(1);
}

/** check for singularities */
static int _MMG5_singul(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1,p2;
  double         ux,uy,uz,vx,vy,vz,dd;
  int            list[_MMG5_LMAX+2],k,nc,xp,nr,ns,nre;
  char           i;

  nre = nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->s++;
      if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) )  continue;
      else if ( MG_EDG(ppt->tag) ) {
        ns = _MMG5_bouler(mesh,mesh->adja,k,i,list,&xp,&nr, _MMG5_LMAX);

        if ( !ns )  continue;
        if ( (xp+nr) > 2 ) {
          ppt->tag |= MG_CRN + MG_REQ;
          nre++;
          nc++;
        }
        else if ( (xp == 1) && (nr == 1) ) {
          ppt->tag |= MG_REQ;
          nre++;
        }
        else if ( xp == 1 && !nr ){
          ppt->tag |= MG_CRN + MG_REQ;
          nre++;
          nc++;
        }
        else if ( nr == 1 && !xp ){
          ppt->tag |= MG_CRN + MG_REQ;
          nre++;
          nc++;
        }
        /* check ridge angle */
        else {
          p1 = &mesh->point[list[1]];
          p2 = &mesh->point[list[2]];
          ux = p1->c[0] - ppt->c[0];
          uy = p1->c[1] - ppt->c[1];
          uz = p1->c[2] - ppt->c[2];
          vx = p2->c[0] - ppt->c[0];
          vy = p2->c[1] - ppt->c[1];
          vz = p2->c[2] - ppt->c[2];
          dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);
          if ( fabs(dd) > _MMG5_EPSD ) {
            dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
            if ( dd > -mesh->info.dhd ) {
              ppt->tag |= MG_CRN;
              nc++;
            }
          }
        }
      }
    }
  }

  /* check for handle */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->s )  continue;
      nr = boulet(mesh,k,i,list);
      if ( nr != ppt->s ) {
        ppt->tag |= MG_CRN + MG_REQ;
        ppt->s = 0;
        nc++;
      }
    }
  }

  /* reset the ppt->s tag */
  for (k=1; k<=mesh->np; ++k) {
    mesh->point[k].s = 0;
  }

  if ( abs(mesh->info.imprim) > 3 && nre > 0 )
    fprintf(stdout,"     %d corners, %d singular points detected\n",nc,nre);
  return(1);
}


/* compute normals at C1 vertices, for C0: tangents */
static int norver(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  MMG5_pxPoint   go;
  double         n[3],dd;
  int            *adja,k,kk,ier,xp,nn,nt,nf,nnr;
  char           i,ii,i1;

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

  /* 1. process C1 vertices, normals */
  nn = xp = nt = nf = nnr = 0;
  ++mesh->base;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( MS_SIN(ppt->tag) || MG_EDG(ppt->tag) ) {
        if ( mesh->nc1 ) {
          if ( ppt->n[0]*ppt->n[0]+ppt->n[1]*ppt->n[1]+ppt->n[2]*ppt->n[2] > 0 )
            ++nnr;
        }

        if ( MG_EDG(ppt->tag) )  xp++;

        continue;
      }
      else if ( ppt->flag == mesh->base )  continue;
      else if ( mesh->nc1 ) {
        if ( ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2] > 0 )
        continue;
      }

      ier = _MMG5_boulen(mesh,mesh->adja,k,i,ppt->n);
      if ( ier ) {
        ppt->flag = mesh->base;
        nn++;
      }
      else
        nf++;
    }
  }

  /* memory to store normals on both sides of ridges */
  mesh->xpmax = MG_MAX(1.5*xp,_MMG5_XPMAX);
  _MMG5_ADD_MEM(mesh,(mesh->xpmax+1)*sizeof(MMG5_xPoint),"boundary points",return(0));
  _MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xpmax+1,MMG5_xPoint);

  if ( xp ) {
    /* 2. process C0 vertices on curves, tangents */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      adja = &mesh->adja[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = _MMG5_inxt2[i];
        ppt = &mesh->point[pt->v[i]];

        if ( ppt->tag & MG_CRN || ppt->flag == mesh->base )  continue;
        else if ( !MG_EDG(pt->tag[i1]) )  continue;

        ier = _MMG5_boulen(mesh,mesh->adja,k,i,n);
        if ( !ier )  continue;

        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                             "larger xpoint table",
                             mesh->xp--;
                             return(0));
        }
        ppt->xp = mesh->xp;
        go = &mesh->xpoint[mesh->xp];
        memcpy(go->n1,n,3*sizeof(double));

        /* compute n2 along ridge */
        if ( pt->tag[i1] & MG_GEO ) {
          if ( adja[i1] ) {
            kk  = adja[i1] / 3;
            ii  = adja[i1] % 3;
            ii  = _MMG5_inxt2[ii];

            ier = _MMG5_boulen(mesh,mesh->adja,kk,ii,n);
            if ( !ier )  continue;
            memcpy(go->n2,n,3*sizeof(double));

            /* compute tangent as intersection of n1 + n2 */
            ppt->n[0] = go->n1[1]*go->n2[2] - go->n1[2]*go->n2[1];
            ppt->n[1] = go->n1[2]*go->n2[0] - go->n1[0]*go->n2[2];
            ppt->n[2] = go->n1[0]*go->n2[1] - go->n1[1]*go->n2[0];
            ppt->flag = mesh->base;
            dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
            if ( dd > _MMG5_EPSD2 ) {
              dd = 1.0 / sqrt(dd);
              ppt->n[0] *= dd;
              ppt->n[1] *= dd;
              ppt->n[2] *= dd;
            }
            ++nt;
            continue;
          }
        }

        /* compute tgte */
        ier = _MMG5_boulec(mesh,mesh->adja,k,i,ppt->n);
        if ( !ier )  continue;
        dd = go->n1[0]*ppt->n[0] + go->n1[1]*ppt->n[1] + go->n1[2]*ppt->n[2];
        ppt->n[0] -= dd*go->n1[0];
        ppt->n[1] -= dd*go->n1[1];
        ppt->n[2] -= dd*go->n1[2];
        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
        if ( dd > _MMG5_EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          ppt->n[0] *= dd;
          ppt->n[1] *= dd;
          ppt->n[2] *= dd;
          ppt->flag = mesh->base;
          ++nt;
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 4 && nn+nt > 0 ) {
    if ( nnr )
      fprintf(stdout,"     %d input normals ignored\n",nnr);
    fprintf(stdout,"     %d normals,  %d tangents updated  (%d failed)\n",nn,nt,nf);
  }

  return(1);
}

/* regularization procedure for derivatives, dual Laplacian */
static int regnor(MMG5_pMesh mesh) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt,p0;
  double  *tabl,n[3],lm1,lm2,dd,nx,ny,nz,res0,res;
  int      i,k,iad,it,nn,nit,iel,ilist,list[_MMG5_LMAX];

  /* assign seed to vertex */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k;
    }
  }

  /* allocate memory for normals */
  _MMG5_SAFE_CALLOC(tabl,3*mesh->np+1,double);

  it   = 0;
  nit  = 2;
  res0 = 0.0;
  lm1  = 0.4;
  lm2  = 0.399;
  while ( it++ < nit ) {
    /* step 1: laplacian */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || ppt->tag > MG_REF )  continue;

      iel = ppt->s;
      assert(iel);
      pt = &mesh->tria[iel];
      i  = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = boulep(mesh,iel,i,list);  

      /* average normal */
      nx = ny = nz = 0.0;
      for (i=1; i<=ilist; i++) {
        p0  = &mesh->point[list[i]];
        if ( p0->tag > MG_REF )  continue;
        nx += p0->n[0];
        ny += p0->n[1];
        nz += p0->n[2];
      }
      dd  = nx*nx + ny*ny + nz*nz;
      if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nx *= dd;
        ny *= dd;
        nz *= dd;
      }

      /* Laplacian */
      iad = 3*(k-1)+1;
      tabl[iad+0] = ppt->n[0] + lm1 * (nx - ppt->n[0]);
      tabl[iad+1] = ppt->n[1] + lm1 * (ny - ppt->n[1]);
      tabl[iad+2] = ppt->n[2] + lm1 * (nz - ppt->n[2]);
    }

    /* step 2: anti-laplacian */
    res = 0;
    nn  = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || ppt->tag > MG_REF )  continue;

      iel = ppt->s;
      assert(iel);
      pt = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = boulep(mesh,iel,i,list);  

      /* average normal */
      nx = ny = nz = 0.0;
      for (i=1; i<=ilist; i++) {
        iad = 3*(list[i]-1) + 1;
        nx += tabl[iad+0];
        ny += tabl[iad+1];
        nz += tabl[iad+2];
      }
      dd  = nx*nx + ny*ny + nz*nz;
      if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        nx *= dd;
        ny *= dd;
        nz *= dd;
      }

      /* antiLaplacian */
      iad = 3*(k-1)+1;
      n[0] = tabl[iad+0] - lm2 * (nx - tabl[iad+0]);
      n[1] = tabl[iad+1] - lm2 * (ny - tabl[iad+1]);
      n[2] = tabl[iad+2] - lm2 * (nz - tabl[iad+2]);
      nn++;
      res += (ppt->n[0]-n[0])*(ppt->n[0]*n[0]) + (ppt->n[1]-n[1])*(ppt->n[1]*n[1]) + (ppt->n[2]-n[2])*(ppt->n[2]*n[2]); 
    }

    if ( it == 1 )  res0 = res;
    if ( res0 > _MMG5_EPSD )  res  = res / res0;
    if ( mesh->info.imprim < 0 || mesh->info.ddebug ) {
      fprintf(stdout,"     iter %5d  res %.3E\r",it,res); 
      fflush(stdout);
    }
    if ( it > 1 && res < _MMG5_EPS )  break;
  }

  /* reset the ppt->s tag */
  for (k=1; k<=mesh->np; ++k) {
    mesh->point[k].s = 0;
  }

  if ( mesh->info.imprim < 0 || mesh->info.ddebug )  fprintf(stdout,"\n");

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     %d normals regularized: %.3e\n",nn,res);

  _MMG5_SAFE_FREE(tabl);
  return(1);
}


/* preprocessing stage: mesh analysis */
int _MMGS_analys(MMG5_pMesh mesh) {

  /* set tria edges tags */
  if ( !assignEdge(mesh) ) {
    fprintf(stderr,"  ## Analysis problem. Exit program.\n");
    return(0);
  }

  /* create adjacency */
  if ( !_MMGS_hashTria(mesh) ) {
    fprintf(stderr,"  ## Hashing problem. Exit program.\n");
    return(0);
  }

  /* delete badly shaped elts */
  /*if ( mesh->info.badkal && !delbad(mesh) ) {
    fprintf(stderr,"  ## Geometry trouble. Exit program.\n");
    return(0);
    }*/

  /* identify connexity */
  if ( !setadj(mesh) ) {
    fprintf(stderr,"  ## Topology problem. Exit program.\n");
    return(0);
  }

  /* check for nomanifold point */
  nmpoints(mesh);

  /* check for ridges */
  if ( mesh->info.dhd > _MMG5_ANGLIM && !setdhd(mesh) ) {
    fprintf(stderr,"  ## Geometry problem. Exit program.\n");
    return(0);
  }

  /* identify singularities */
  if ( !_MMG5_singul(mesh) ) {
    fprintf(stderr,"  ## Singularity problem. Exit program.\n");
    return(0);
  }

  /* define normals */
  if ( !mesh->xp ) {
    if ( !norver(mesh) ) {
      fprintf(stderr,"  ## Normal problem. Exit program.\n");
      return(0);
    }
    /* regularize normals */
    if ( mesh->info.nreg && !regnor(mesh) ) {
      fprintf(stderr,"  ## Normal regularization problem. Exit program.\n");
      return(0);
    }
  }

  return(1);
}

