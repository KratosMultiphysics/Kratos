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

#include "libmmgs_private.h"
#include "mmgcommon_private.h"

/**
 * \param mesh pointer to the mesh
 *
 * \return 1 if success, 0 if fail
 *
 * topology: set adjacent, detect Moebius, flip faces, count connected comp.
 *
 */
int MMGS_setadj(MMG5_pMesh mesh){
  MMG5_pTria   pt,pt1;
  MMG5_int     *adja,*adjb,adji1,adji2,*pile,iad,ipil,ip1,ip2,gen;
  MMG5_int     k,kk,iel,jel,nvf,nf,nr,nt,nre,nreq,ncc,ned,ref;
  uint16_t     tag;
  int8_t       i,ii,i1,i2,ii1,ii2,voy;

  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  nvf = nf = ncc = ned = 0;

  MMG5_SAFE_MALLOC(pile,mesh->nt+1,MMG5_int,return 0);

  pile[1] = 1;
  ipil    = 1;

  while ( ipil > 0 ) {
    ncc++;

    do {
      k  = pile[ipil--];
      pt = &mesh->tria[k];
      pt->flag = 1;
      if ( !MG_EOK(pt) )  continue;

      pt->cc = ncc;
      adja = &mesh->adja[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];

        if ( !mesh->point[ip1].tmp )  mesh->point[ip1].tmp = ++nvf;
        if ( !mesh->point[ip2].tmp )  mesh->point[ip2].tmp = ++nvf;

        if ( MG_EDG(pt->tag[i]) || pt->tag[i] & MG_REQ ) {
          mesh->point[ip1].tag |= pt->tag[i];
          mesh->point[ip2].tag |= pt->tag[i];
        }

        /* open boundary */
        if ( !adja[i] ) {
          pt->tag[i] |= MG_GEO;
          mesh->point[ip1].tag |= MG_GEO;
          mesh->point[ip2].tag |= MG_GEO;
          ned++;
          continue;
        }

        kk = adja[i] / 3;
        ii = adja[i] % 3;
        if ( kk > k )  ned++;

        /* store adjacent */
        pt1 = &mesh->tria[kk];

        /* correct edge tag */
        if ( (pt->tag[i] & MG_NOM) && !(pt1->tag[ii] & MG_NOM) ) {
          pt1->tag[ii] = pt->tag[i];
          pt1->edg[ii] = pt->edg[i];
          mesh->point[ip1].tag |= MG_NOM;
          mesh->point[ip2].tag |= MG_NOM;
          continue;
        }
        if ( (pt1->tag[ii] & MG_NOM) && !(pt->tag[i] & MG_NOM) ) {
          pt->tag[i] = pt1->tag[ii];
          pt->edg[i] = pt1->edg[ii];
          mesh->point[ip1].tag |= MG_NOM;
          mesh->point[ip2].tag |= MG_NOM;
          continue;
        }

        if ( MMG5_abs(pt1->ref) != MMG5_abs(pt->ref) ) {
          pt->tag[i]   |= MG_REF;
          pt1->tag[ii] |= MG_REF;
          mesh->point[ip1].tag |= MG_REF;
          mesh->point[ip2].tag |= MG_REF;
        }

        /* Do not treat adjacent through a non-manifold edge */
        if ( pt1->tag[ii] & MG_NOM ) {
          continue;
        }

        /* store adjacent */
        if ( !pt1->flag ) {
          pt1->flag    = 1;
            pile[++ipil] = kk;
          }

        /* check orientation */
        ii1 = MMG5_inxt2[ii];
        ii2 = MMG5_iprv2[ii];
        if ( pt1->v[ii1] == ip1 ) {
          /* Moebius strip */
          assert ( pt1->base );
          if ( pt1->base < 0 ) {
            pt1->ref      = -MMG5_abs(pt1->ref);
            /* Add MG_NOM flag because it allows neighbours to have non
             * consistent orientations */
            pt->tag[i]   |= MG_REF + MG_NOM;
            pt1->tag[ii] |= MG_REF + MG_NOM;
          }
          /* flip orientation */
          else {
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
        else {
          /* Mark triangles that have a consistent orientation with their
           * neighbours */
          pt1->base =  -pt1->base;
        }
      }
    }
    while ( ipil > 0 );

    /* find next unmarked triangle */
    ipil = 0;
    for (kk=1; kk<=mesh->nt; kk++) {
      pt = &mesh->tria[kk];
      if ( MG_EOK(pt) && (pt->cc == 0) ) {
        ipil = 1;
        pile[ipil] = kk;
        pt->flag   = 1;
        break;
      }
    }
  }

  /* bilan */
  nr = nre = nreq = nt = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    nt++;
    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( ( !MG_EDG(pt->tag[i]) ) && ( !(pt->tag[i] & MG_REQ) ) )  continue;

      jel  = adja[i] / 3;
      if ( !jel || jel > k ) {
        if ( pt->tag[i] & MG_GEO )  nr++;
        if ( pt->tag[i] & MG_REF )  nre++;
        if ( pt->tag[i] & MG_REQ )  nreq++;
      }
    }
  }

  if ( mesh->info.ddebug ) {
    fprintf(stdout,"  a- ridges: %" MMG5_PRId " found.\n",nr);
    fprintf(stdout,"  a- requir: %" MMG5_PRId " found.\n",nreq);
    fprintf(stdout,"  a- connex: %" MMG5_PRId " connected component(s)\n",ncc);
    fprintf(stdout,"  a- orient: %" MMG5_PRId " flipped\n",nf);
  }
  else if ( abs(mesh->info.imprim) > 3 ) {
    gen = (2 - nvf + ned - nt) / 2;
    fprintf(stdout,"     Connected component: %" MMG5_PRId ",  genus: %" MMG5_PRId ",   reoriented: %" MMG5_PRId "\n",ncc,gen,nf);
    fprintf(stdout,"     Edges: %" MMG5_PRId ",  tagged: %" MMG5_PRId ",  ridges: %" MMG5_PRId ", required: %" MMG5_PRId ", refs: %" MMG5_PRId "\n",
            ned,nr+nre+nreq,nr,nreq,nre);
  }

  MMG5_SAFE_FREE(pile);
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 *
 * \return 1 if succeed, 0 if fail
 *
 * Detect non manifold points
 */
static void nmpoints(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0;
  MMG5_int        k,np,numt,iel,jel,nmp,*adja;
  int8_t          i0,i1,i,jp;
  
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
        i1  = MMG5_inxt2[i0];
        adja = &mesh->adja[3*(iel-1)+1];
        jel = adja[i1] / 3;
        jp  = adja[i1] % 3;
        jp  = MMG5_inxt2[jp];
      }
      while ( jel && (jel != numt) && (jel !=k) );

      /* Ball has been completely traveled without meeting numt */
      if ( jel == k ) {
        if ( !(p0->tag & MG_CRN) || !(p0->tag & MG_REQ) ) {
          nmp++;
          // p0->tag |= MG_CRN + MG_REQ; // Algiane 2022: this line has been
          // commented by commit 1f592c. I think that it is a mistake (some
          // forgotten debug thing). It seems that in any case, nm points are
          // marked as CRN and REQ when checking handles in MMGS_singul
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
        i1  = MMG5_iprv2[i0];
        adja = &mesh->adja[3*(iel-1)+1];
        jel = adja[i1] / 3;
        jp  = adja[i1] % 3;
        jp  = MMG5_iprv2[jp];
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
    fprintf(stdout,"  ## %" MMG5_PRId " non manifold points detected\n",nmp);
}

/** improve badly shaped elts for isotropic mesh */
/* static int delbad(MMG5_pMesh mesh) { */
/*   MMG5_pTria    pt; */
/*   MMG5_pPoint   p[3]; */
/*   double        s,kal,declic,ux,uy,uz,vx,vy,vz; */
/*   MMG5_int      *adja,k,iel,nd,ndd; */
/*   int           it; */
/*   int8_t   i,ia,i1,i2,j,typ; */

/*   it = ndd = 0; */
/*   declic = MMGS_BADKAL / MMGS_ALPHAD; */

/*   do { */
/*     nd = 0; */
/*     for (k=1; k<=mesh->nt; k++) { */
/*       pt = &mesh->tria[k]; */
/*       if ( !MG_EOK(pt) )  continue; */

/*       kal = MMG5_calelt(mesh,NULL,pt); */
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
/*           i1 = MMG5_inxt2[ia]; */
/*           i2 = MMG5_iprv2[ia]; */
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
          
/*           if ( !delElt(mesh,k) )  return 0; */
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
/*     if ( nd && (mesh->info.ddebug || mesh->info.imprim < -1) )  fprintf(stdout,"     %" MMG5_PRId " improved\n",nd); */
/*   } */
/*   while ( nd > 0 && ++it < 5 ); */

/*   if ( abs(mesh->info.imprim) > 4 ) */
/*     fprintf(stdout,"     %" MMG5_PRId " bad elements improved\n",ndd); */

/*   return 1; */
/* } */


/**
 * \param mesh pointer to the mesh structure.
 *
 * \return 1 if succeed, 0 if fail
 *
 * check for ridges: dihedral angle
 */
static int setdhd(MMG5_pMesh mesh) {
  MMG5_pTria    pt,pt1;
  double        n1[3],n2[3],dhd;
  MMG5_int      *adja,k,kk,nr;
  int8_t        i,ii,i1,i2;

  nr = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    /* triangle normal */
    MMG5_nortri(mesh,pt,n1);
    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_GEO )  continue;
      kk = adja[i] / 3;
      ii = adja[i] % 3;

      /* check angle w. neighbor */
      if ( k < kk ) {
        pt1 = &mesh->tria[kk];
        MMG5_nortri(mesh,pt1,n2);
        dhd = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
        if ( dhd <= mesh->info.dhd ) {
          pt->tag[i]   |= MG_GEO;
          pt1->tag[ii] |= MG_GEO;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_GEO;
          mesh->point[pt->v[i2]].tag |= MG_GEO;
          nr++;
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 4 && nr > 0 )
    fprintf(stdout,"     %" MMG5_PRId " ridges updated\n",nr);

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 *
 * \return 1 if succeed, 0 if fail
 *
 * check for singularities
 *
 */
static int MMG5_singul(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1,p2;
  double         ux,uy,uz,vx,vy,vz,dd;
  MMG5_int       list[MMG5_TRIA_LMAX+2],listref[MMG5_TRIA_LMAX+2],k,nc,nre;
  int            xp,nr,ns;
  int8_t         i;

  nre = nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->s++;
      if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) )  continue;
      else if ( MG_EDG(ppt->tag) ) {
        ns = MMG5_bouler(mesh,mesh->adja,k,i,list,listref,&xp,&nr, MMG5_TRIA_LMAX);

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
          if ( fabs(dd) > MMG5_EPSD ) {
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
      int8_t dummy;
      nr = MMG5_boulet(mesh,k,i,list,1,&dummy);
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
    fprintf(stdout,"     %" MMG5_PRId " corners, %" MMG5_PRId " singular points detected\n",nc,nre);
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 *
 * \return 1 if succeed, 0 if fail
 *
 * Compute normals at C1 vertices, for C0: tangents
 * This function allocate the xpoint array. A point will have an xpoint if:
 *   - it is along a reference edge, the xpoint then stores the normal at point while the n field of point containt the tangent at ref edge.
 *   - it is along a ridge, the xpoint then stores both normals at point while ppt->n stores the tangent
 *
 * Corner, required and regular points don't have xpoints.
 *
 */
static int norver(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  MMG5_pxPoint   go;
  double         n[3],dd;
  MMG5_int       *adja,k,kk,ier,xp,nn,nt,nf,nnr;
  int8_t         i,ii,i1;

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

      ier = MMG5_boulen(mesh,mesh->adja,k,i,ppt->n);
      if ( ier ) {
        ppt->flag = mesh->base;
        nn++;
      }
      else
        nf++;
    }
  }

  /* memory to store normals on both sides of ridges */
  mesh->xpmax = MG_MAX(1.5*xp,MMGS_XPMAX);
  /* no need to have more xpoint than point */
  mesh->xpmax = MG_MIN(mesh->npmax,mesh->xpmax);
  MMG5_ADD_MEM(mesh,(mesh->xpmax+1)*sizeof(MMG5_xPoint),"boundary points",return 0);
  MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,return 0);

  if ( xp ) {
    /* 2. process C0 vertices on curves, tangents */
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;

      adja = &mesh->adja[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = MMG5_inxt2[i];
        ppt = &mesh->point[pt->v[i]];

        if ( ppt->tag & MG_CRN || ppt->flag == mesh->base )  continue;
        else if ( !MG_EDG(pt->tag[i1]) )  continue;

        /* As we skip non-manifold point, the edge should be manifold */
        assert ( (!(MG_NOM & pt->tag[i1])) && "Unexpected non-manifold edge" );

        ier = MMG5_boulen(mesh,mesh->adja,k,i,n);
        if ( !ier )  continue;

        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                             "larger xpoint table",
                             mesh->xp--;
                             return 0);
        }
        ppt->xp = mesh->xp;
        go = &mesh->xpoint[mesh->xp];
        memcpy(go->n1,n,3*sizeof(double));

        /* compute n2 along ridge */
        if ( pt->tag[i1] & MG_GEO ) {
          if ( adja[i1] ) {
            kk  = adja[i1] / 3;
            ii  = adja[i1] % 3;
            ii  = MMG5_inxt2[ii];

            ier = MMG5_boulen(mesh,mesh->adja,kk,ii,n);
            if ( !ier )  continue;
            memcpy(go->n2,n,3*sizeof(double));

            /* compute tangent as intersection of n1 + n2 */
            ppt->n[0] = go->n1[1]*go->n2[2] - go->n1[2]*go->n2[1];
            ppt->n[1] = go->n1[2]*go->n2[0] - go->n1[0]*go->n2[2];
            ppt->n[2] = go->n1[0]*go->n2[1] - go->n1[1]*go->n2[0];
            ppt->flag = mesh->base;
            dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
            if ( dd > MMG5_EPSD2 ) {
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
        ier = MMG5_boulec(mesh,mesh->adja,k,i,ppt->n);
        if ( !ier )  continue;
        dd = go->n1[0]*ppt->n[0] + go->n1[1]*ppt->n[1] + go->n1[2]*ppt->n[2];
        ppt->n[0] -= dd*go->n1[0];
        ppt->n[1] -= dd*go->n1[1];
        ppt->n[2] -= dd*go->n1[2];
        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
        if ( dd > MMG5_EPSD2 ) {
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
      fprintf(stdout,"     %" MMG5_PRId " input normals ignored\n",nnr);
    fprintf(stdout,"     %" MMG5_PRId " normals,  %" MMG5_PRId " tangents updated  (%" MMG5_PRId " failed)\n",nn,nt,nf);
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param pt pointer to current triangle
 * \param k number of current point
 * \param c newly computed coordinates (giving negative area)
 * \param n normal of triangle before regularization
 *
 * \return 0 if fail, 1 if success
 *
 * In coordinate regularization, performs a dichotomy between previous point /
 * and newly computed point in the case of negative area
 *
 */
static inline int MMGS_dichotomy(MMG5_pMesh mesh, MMG5_pTria pt, MMG5_int k, double *c, double *n) {

  MMG5_pPoint  ppt;
  double       to,tp,t,nnew[3],p[3],o[3],result;
  int          it,maxit,pos,ier;

  it = 0;
  maxit = 5;
  to = 0.0;
  tp = 1.0;
  t = 0.5;
  pos = 0;

  ppt = &mesh->point[k];

  /* initial coordinates of point before regularization */
  o[0] = ppt->c[0];
  o[1] = ppt->c[1];
  o[2] = ppt->c[2];

  /* initial coordinates of new point */
  p[0] = c[0];
  p[1] = c[1];
  p[2] = c[2];

  do {
    mesh->point[0].c[0] = o[0] + t*(p[0] - o[0]);
    mesh->point[0].c[1] = o[1] + t*(p[1] - o[1]);
    mesh->point[0].c[2] = o[2] + t*(p[2] - o[2]);

    MMG5_nortri(mesh, pt, nnew);
    MMG5_dotprod(3,n,nnew,&result);

    if ( result <= 0.0 ) {
      tp = t;
    }
    else {
      c[0] = mesh->point[0].c[0];
      c[1] = mesh->point[0].c[1];
      c[2] = mesh->point[0].c[2];
      to = t;
      pos = 1;
    }

    t = 0.5*(to + tp);
  }
  while ( ++it < maxit );

  if ( pos ) {
    return 1;
  }
  else
    return 0;
}

/**
 * \param mesh pointer to a MMG5 mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Regularization procedure for vertices coordinates,
 * dual Laplacian for a surface mesh.
 *
 */
int MMGS_regver(MMG5_pMesh mesh) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt,p0;
  MMG5_Tria     tnew;
  double        *tabl,c[3],n[3],nnew[3],*cptr,lm1,lm2,cx,cy,cz,res0,res,result;
  int           i,it,nit,ilist,noupdate,ier;
  MMG5_int      k,kt,nn,iel,list[MMG5_LMAX],tlist[MMG5_LMAX],*adja,iad;

  /* assign seed to vertex */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k;
    }
  }

  /* allocate memory for coordinates */
  MMG5_SAFE_CALLOC(tabl,3*mesh->np+1,double,return 0);

  /* Pointer toward the suitable adjacency array */
  adja = mesh->adja;

  it   = 0;
  nit  = 10;
  res0 = 0.0;
  lm1  = mesh->info.lxreg;
  lm2  = 0.99*lm1;
  while ( it++ < nit ) {
    /* step 1: laplacian */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];

      iad = 3*(k-1)+1;
      tabl[iad+0] = ppt->c[0];
      tabl[iad+1] = ppt->c[1];
      tabl[iad+2] = ppt->c[2];

      if ( !MG_VOK(ppt) )  continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) continue;

      iel = ppt->s;

      pt = &mesh->tria[iel];
      i  = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = MMG5_boulep(mesh,iel,i,adja,list,tlist);

      /* average coordinates */
      cx = cy = cz = 0.0;
      for (i=1; i<=ilist; i++) {
        p0  = &mesh->point[list[i]];

        cptr = p0->c;
        cx += cptr[0];
        cy += cptr[1];
        cz += cptr[2];
      }
      cx /= ilist;
      cy /= ilist;
      cz /= ilist;

      /* Laplacian */
      cptr = ppt->c;
      tabl[iad+0] = cptr[0] + lm1 * (cx - cptr[0]);
      tabl[iad+1] = cptr[1] + lm1 * (cy - cptr[1]);
      tabl[iad+2] = cptr[2] + lm1 * (cz - cptr[2]);
    }

    /* step 2: anti-laplacian */
    res = 0;
    nn  = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];

      if ( !MG_VOK(ppt) )  continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) continue;

      iel = ppt->s;

      pt = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = MMG5_boulep(mesh,iel,i,adja,list,tlist);

      /* average normal */
      cx = cy = cz = 0.0;
      for (i=1; i<=ilist; i++) {
        iad = 3*(list[i]-1) + 1;
        cx += tabl[iad+0];
        cy += tabl[iad+1];
        cz += tabl[iad+2];
      }
      cx /= ilist;
      cy /= ilist;
      cz /= ilist;

      /* antiLaplacian */
      iad = 3*(k-1)+1;
      c[0] = tabl[iad+0] - lm2 * (cx - tabl[iad+0]);
      c[1] = tabl[iad+1] - lm2 * (cy - tabl[iad+1]);
      c[2] = tabl[iad+2] - lm2 * (cz - tabl[iad+2]);

      cptr = ppt->c;

      mesh->point[0].c[0] = c[0];
      mesh->point[0].c[1] = c[1];
      mesh->point[0].c[2] = c[2];

      /* check for negative areas */
      noupdate = 0;
      for (kt = 0 ; kt<ilist ; kt++) {
        pt = &mesh->tria[tlist[kt]];

        if ( !MG_EOK(pt) ) continue;

        MMG5_nortri(mesh, pt, n);

        for (i=0;i<3;i++) {
          tnew.v[i] = pt->v[i];
        }

        i = 0;
        if ( pt->v[1] == k ) i = 1;
        if ( pt->v[2] == k ) i = 2;

        tnew.v[i] = 0;

        MMG5_nortri(mesh, &tnew, nnew);
        MMG5_dotprod(3,n,nnew,&result);
        if ( result < 0.0 ) {
          if (!MMGS_dichotomy(mesh,&tnew,k,c,n))
            noupdate = 1;
          continue;
        }
      }
      if ( !noupdate ) {
        res += (cptr[0]-c[0])*(cptr[0]-c[0]) + (cptr[1]-c[1])*(cptr[1]-c[1]) + (cptr[2]-c[2])*(cptr[2]-c[2]);
        cptr[0] = c[0];
        cptr[1] = c[1];
        cptr[2] = c[2];
        nn++;
      }
    }
      if ( it == 1 )  res0 = res;
      if ( res0 > MMG5_EPSD )  res  = res / res0;
      if ( mesh->info.imprim < -1 || mesh->info.ddebug ) {
        fprintf(stdout,"     iter %5d  res %.3E\r",it,res);
        fflush(stdout);
      }
      if ( it > 1 && res < MMG5_EPS )  break;
    }
  /* reset the ppt->s tag */
  for (k=1; k<=mesh->np; ++k) {
    mesh->point[k].s    = 0;
  }

  if ( mesh->info.imprim < -1 || mesh->info.ddebug )  fprintf(stdout,"\n");

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     %" MMG5_PRId " coordinates regularized: %.3e\n",nn,res);

  MMG5_SAFE_FREE(tabl);
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Remove duplicated triangles.
 *
 */
int MMGS_remDup(MMG5_pMesh mesh) {
  MMG5_Hash     hash;
  MMG5_pTria    ptt;
  MMG5_int      k,jel,dup;

  if ( !mesh->nt ) return 1;

  /* Hash triangles */
  if ( ! MMG5_hashNew(mesh,&hash,(MMG5_int)(0.51*mesh->nt),(MMG5_int)(1.51*mesh->nt)) ) {
    return 0;
  }

  dup = 0;
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    jel = MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k);
    if ( !jel ) {
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
    else if ( jel > 0 ) {
      ++dup;
      /* Remove duplicated face */
      MMGS_delElt(mesh,k);
    }
  }

  if ( abs(mesh->info.imprim) > 5 && dup > 0 ) {
    fprintf(stdout,"  ## ");  fflush(stdout);
    if ( dup > 0 )  fprintf(stdout," %"MMG5_PRId" duplicate removed",dup);
    fprintf(stdout,"\n");
  }

  MMG5_DEL_MEM(mesh,hash.item);

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 *
 * \return 1 if succeed, 0 if fail
 *
 * Preprocessing stage: mesh analysis.
 *
 */
int MMGS_analys_for_norver(MMG5_pMesh mesh) {

  /* create adjacency */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* delete badly shaped elts */
  /*if ( mesh->info.badkal && !delbad(mesh) ) {
    fprintf(stderr,"\n  ## Geometry trouble. Exit program.\n");
    return 0;
    }*/

  /* identify connexity */
  if ( !MMGS_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    return 0;
  }

  /* check for nomanifold point */
  nmpoints(mesh);

  /* check for ridges */
  if ( mesh->info.dhd > MMG5_ANGLIM && !setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    return 0;
  }

  /* identify singularities */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## Singularity problem. Exit program.\n");
    return 0;
  }

  /* define normals */
  if ( !mesh->xp ) {
    if ( !norver(mesh) ) {
      fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
      return 0;
    }
    /* regularize normals */
    if ( mesh->info.nreg && !MMG5_regnor(mesh) ) {
      fprintf(stderr,"\n  ## Normal regularization problem. Exit program.\n");
      return 0;
    }
  }
  return 1;
}


/* preprocessing stage: mesh analysis */
int MMGS_analys(MMG5_pMesh mesh) {

  /* Update tags stored into tria */
  if ( !MMGS_bdryUpdate(mesh) ) {
    fprintf(stderr,"\n  ## Analysis problem. Exit program.\n");
    return 0;
  }

  /* set edges tags and refs to tria */
  if ( !MMGS_assignEdge(mesh) ) {
    fprintf(stderr,"\n  ## Analysis problem. Exit program.\n");
    return 0;
  }

  /* create adjacency */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* delete badly shaped elts */
  /*if ( mesh->info.badkal && !delbad(mesh) ) {
    fprintf(stderr,"\n  ## Geometry trouble. Exit program.\n");
    return 0;
    }*/

  /* identify connexity */
  if ( !MMGS_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    return 0;
  }

  /* check for nomanifold point */
  nmpoints(mesh);

  /* check for ridges */
  if ( mesh->info.dhd > MMG5_ANGLIM && !setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    return 0;
  }

  /* identify singularities */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## Singularity problem. Exit program.\n");
    return 0;
  }

  /* regularize vertices coordinates */
  if ( mesh->info.xreg && !MMGS_regver(mesh) ){
    fprintf(stderr,"\n  ## Coordinates regularization problem. Exit program.\n");
    return 0;
  }

  /* define normals */
  if ( !mesh->xp ) {
    if ( !norver(mesh) ) {
      fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
      return 0;
    }
    /* regularize normals */
    if ( mesh->info.nreg && !MMG5_regnor(mesh) ) {
      fprintf(stderr,"\n  ## Normal regularization problem. Exit program.\n");
      return 0;
    }
  }

  return 1;
}

