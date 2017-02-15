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
 * \file mmg3d/analys_3d.c
 * \brief Mesh analysis.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmg3d.h"

/**
 * \param mesh pointer towarad the mesh structure.
 *
 * Set all boundary triangles to required and add a tag to detect that they are
 * not realy required.
 *
 */
static inline void _MMG5_reqBoundaries(MMG5_pMesh mesh) {
  MMG5_pTria     ptt;
  int            k;

  /* The MG_REQ+MG_NOSURF tag mark the boundary edges that we dont want to touch
   * but that are not really required (-nosurf option) */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !(ptt->tag[0] & MG_REQ) ) {
      ptt->tag[0] |= MG_REQ;
      ptt->tag[0] |= MG_NOSURF;
    }

    if ( !(ptt->tag[1] & MG_REQ) ) {
      ptt->tag[1] |= MG_REQ;
      ptt->tag[1] |= MG_NOSURF;
    }

    if ( !(ptt->tag[2] & MG_REQ) ) {
      ptt->tag[2] |= MG_REQ;
      ptt->tag[2] |= MG_NOSURF;
    }
  }

  return;
}


/**
 * \param mesh pointer towarad the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * topology: set adjacent, detect Moebius, flip faces, count connected comp.
 *
 */
static int _MMG5_setadj(MMG5_pMesh mesh){
  MMG5_pTria   pt,pt1;
  MMG5_pPoint  ppt;
  int    *adja,*adjb,adji1,adji2,*pile,iad,ipil,ip1,ip2,gen;
  int     k,kk,iel,jel,nf,np,nr,nt,nre,nreq,ncc,ned,nvf,edg;
  int16_t tag;
  char    i,ii,i1,i2,ii1,ii2,voy;

  nvf = nf = ncc = ned = 0;
  _MMG5_SAFE_MALLOC(pile,mesh->nt+1,int);

  pile[1] = 1;
  ipil    = 1;
  pt = &mesh->tria[1];
  pt->flag = 1;

  while ( ipil > 0 ) {
    ncc++;
    do {
      k  = pile[ipil--];
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;
      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        i1  = _MMG5_inxt2[i];
        i2  = _MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        if ( !mesh->point[ip1].tmp )  mesh->point[ip1].tmp = ++nvf;
        if ( !mesh->point[ip2].tmp )  mesh->point[ip2].tmp = ++nvf;
        if ( MG_EDG(pt->tag[i]) || pt->tag[i] & MG_REQ ) {
              tag = mesh->point[ip1].tag;
              mesh->point[ip1].tag |= pt->tag[i];
              // Remove the MG_NOSURF tag if the vertex is really required.
              if ( (tag & MG_REQ) && !(tag & MG_NOSURF) ) {
                mesh->point[ip1].tag &= ~MG_NOSURF;
              }
              tag = mesh->point[ip2].tag;
              mesh->point[ip2].tag |= pt->tag[i];
              // Remove the MG_NOSURF tag if the vertex is really required.
              if ( (tag & MG_REQ) && !(tag & MG_NOSURF) ) {
                mesh->point[ip2].tag &= ~MG_NOSURF;
              }
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

        /* non manifold edge */
        if ( pt->tag[i] & MG_NOM ) {
          mesh->point[ip1].tag |= MG_NOM;
          mesh->point[ip2].tag |= MG_NOM;
          continue;
        }

        /* store adjacent */
        pt1 = &mesh->tria[kk];
        if ( abs(pt1->ref) != abs(pt->ref) ) {
          pt->tag[i]   |= MG_REF;
          if ( !(pt->tag[i] & MG_NOM) )  pt1->tag[ii] |= MG_REF;
          mesh->point[ip1].tag |= MG_REF;
          mesh->point[ip2].tag |= MG_REF;
        }

        if ( pt1->flag == 0 ) {
          pt1->flag    = ncc;
          pile[++ipil] = kk;
        }

        /* check orientation */
        ii1 = _MMG5_inxt2[ii];
        ii2 = _MMG5_iprv2[ii];
        if ( pt1->v[ii1] == ip1 ) {
          /* Moebius strip */
          if ( pt1->base < 0 ) {
            fprintf(stderr,"  ## Orientation problem (1).\n");
            return(0);
          }
          /* flip orientation */
          else {
            pt1->base   = -pt1->base;
            pt1->v[ii1] = ip2;
            pt1->v[ii2] = ip1;

            /* update adj */
            iad   = 3*(kk-1)+1;
            adjb  = &mesh->adjt[iad];
            adji1 = mesh->adjt[iad+ii1];
            adji2 = mesh->adjt[iad+ii2];
            adjb[ii1] = adji2;
            adjb[ii2] = adji1;
            tag = pt1->tag[ii1];
            pt1->tag[ii1] = pt1->tag[ii2];
            pt1->tag[ii2] = tag;
            edg = pt1->edg[ii1];
            pt1->edg[ii1] = pt1->edg[ii2];
            pt1->edg[ii2] = edg;

            /* modif voyeurs */
            if ( adjb[ii1] ) {
              iel = adjb[ii1] / 3;
              voy = adjb[ii1] % 3;
              mesh->adjt[3*(iel-1)+1+voy] = 3*kk + ii1;
            }
            if ( adjb[ii2] ) {
              iel = adjb[ii2] / 3;
              voy = adjb[ii2] % 3;
              mesh->adjt[3*(iel-1)+1+voy] = 3*kk + ii2;
            }
            nf++;
          }
        }
      }
    }
    while ( ipil > 0 );

    /* find next unmarked triangle */
    ipil = 0;
    for (kk=1; kk<=mesh->nt; kk++) {
      pt = &mesh->tria[kk];
      if ( MG_EOK(pt) && (pt->flag == 0) ) {
        ipil = 1;
        pile[ipil] = kk;
        pt->flag   = ncc+1;
        break;
      }
    }
  }

  /* bilan */
  np = nr = nre = nreq = nt = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    nt++;
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->tmp ) {
        ppt->tmp = 1;
        np++;
      }
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
    fprintf(stdout,"  a- ridges: %d found.\n",nr);
    fprintf(stdout,"  a- requir: %d found.\n",nreq);
    fprintf(stdout,"  a- connex: %d connected component(s)\n",ncc);
    fprintf(stdout,"  a- orient: %d flipped\n",nf);
  }
  else if ( abs(mesh->info.imprim) > 3 ) {
    gen = (2 - nvf + ned - nt) / 2;
    fprintf(stdout,"     Connected component: %d,  genus: %d,   reoriented: %d\n",ncc,gen,nf);
    fprintf(stdout,"     Edges: %d,  tagged: %d,  ridges: %d, required: %d, refs: %d\n",
            ned,nr+nre+nreq,nr,nreq,nre);
  }
  _MMG5_SAFE_FREE(pile);
  return(1);
}

/** check for ridges: dihedral angle */
static int _MMG5_setdhd(MMG5_pMesh mesh) {
  MMG5_pTria    pt,pt1;
  double        n1[3],n2[3],dhd;
  int          *adja,k,kk,ne,nr;
  char          i,ii,i1,i2;

  ne = nr = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    /* triangle normal */
    _MMG5_nortri(mesh,pt,n1);
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      kk  = adja[i] / 3;
      ii  = adja[i] % 3;
      if ( !kk ) {
        pt->tag[i] |= MG_GEO;
        i1 = _MMG5_inxt2[i];
        i2 = _MMG5_inxt2[i1];
        mesh->point[pt->v[i1]].tag |= MG_GEO;
        mesh->point[pt->v[i2]].tag |= MG_GEO;
        nr++;
      }
      else if ( k < kk ) {
        pt1 = &mesh->tria[kk];
        /* reference curve */
        if ( pt1->ref != pt->ref ) {
          pt->tag[i]   |= MG_REF;
          pt1->tag[ii] |= MG_REF;
          i1 = _MMG5_inxt2[i];
          i2 = _MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_REF;
          mesh->point[pt->v[i2]].tag |= MG_REF;
          ne++;
        }
        /* check angle w. neighbor */
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
  if ( abs(mesh->info.imprim) > 3 && nr > 0 )
    fprintf(stdout,"     %d ridges, %d edges updated\n",nr,ne);

  return(1);
}

/** check for singularities */
static int _MMG5_singul(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1,p2;
  double         ux,uy,uz,vx,vy,vz,dd;
  int            list[MMG3D_LMAX+2],k,nc,xp,nr,ns,nre;
  char           i;

  nre = nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) )
        continue;
      else if ( MG_EDG(ppt->tag) ) {
        ns = _MMG5_bouler(mesh,mesh->adjt,k,i,list,&xp,&nr,MMG3D_LMAX);

        if ( !ns )  continue;
        if ( (xp+nr) > 2 ) {
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( (xp == 1) && (nr == 1) ) {
          ppt->tag |= MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
        }
        else if ( xp == 1 && !nr ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( nr == 1 && !xp ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
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

  if ( abs(mesh->info.imprim) > 3 && nre > 0 )
    fprintf(stdout,"     %d corners, %d singular points detected\n",nc,nre);
  return(1);
}

/** compute normals at C1 vertices, for C0: tangents */
static int _MMG5_norver(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  MMG5_xPoint    *pxp;
  double         n[3],dd;
  int            *adja,k,kk,ng,nn,nt,nf,nnr;
  char           i,ii,i1;

  /* recomputation of normals only if mesh->xpoint has been freed */
  if ( mesh->xpoint ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: no research of boundary points");
      fprintf(stdout," and normals of mesh. ");
      fprintf(stdout,"mesh->xpoint must be freed to enforce analysis.\n");
    }
    return(1);
  }

  /* identify boundary points */
  ++mesh->base;
  mesh->xp = 0;
  nnr      = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->flag == mesh->base )  continue;
      else {
        ++mesh->xp;
        ppt->flag = mesh->base;
        if ( mesh->nc1 ) {
          if ( ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2] > 0 ) {
            if ( ppt->tag & MG_CRN || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) {
              ++nnr;
              continue;
            }
            ppt->xp = -1;
          }
        }
      }
    }
  }

  /* memory to store normals for boundary points */
  mesh->xpmax  = MG_MAX( (long long)(1.5*mesh->xp),mesh->npmax);

  _MMG5_ADD_MEM(mesh,(mesh->xpmax+1)*sizeof(MMG5_xPoint),"boundary points",return(0));
  _MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xpmax+1,MMG5_xPoint);

  /* compute normals + tangents */
  nn = ng = nt = nf = 0;
  mesh->xp = 0;
  ++mesh->base;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->tag & MG_CRN || ppt->tag & MG_NOM || ppt->flag == mesh->base )  continue;

      /* C1 point */
      if ( !MG_EDG(ppt->tag) ) {

        if ( (!mesh->nc1) ||
             ppt->n[0]*ppt->n[0]+ppt->n[1]*ppt->n[1]+ppt->n[2]*ppt->n[2]<=_MMG5_EPSD2 ) {
          if ( !_MMG5_boulen(mesh,mesh->adjt,k,i,ppt->n) ) {
            ++nf;
            continue;
          }
          else ++nn;
        }

        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                             "larger xpoint table",
                             mesh->xp--;
                             return(0));
        }
        ppt->xp = mesh->xp;
        pxp = &mesh->xpoint[ppt->xp];
        memcpy(pxp->n1,ppt->n,3*sizeof(double));
        ppt->n[0] = ppt->n[1] = ppt->n[2] = 0.;
        ppt->flag = mesh->base;

      }

      /* along ridge-curve */
      i1  = _MMG5_inxt2[i];
      if ( !MG_EDG(pt->tag[i1]) )  continue;
      else if ( !_MMG5_boulen(mesh,mesh->adjt,k,i,n) ) {
        ++nf;
        continue;
      }
      ++mesh->xp;
      if(mesh->xp > mesh->xpmax){
        _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                           "larger xpoint table",
                           mesh->xp--;
                           return(0));
      }
      ppt->xp = mesh->xp;
      pxp = &mesh->xpoint[ppt->xp];
      memcpy(pxp->n1,n,3*sizeof(double));

      if ( pt->tag[i1] & MG_GEO && adja[i1] > 0 ) {
        kk = adja[i1] / 3;
        ii = adja[i1] % 3;
        ii = _MMG5_inxt2[ii];
        if ( !_MMG5_boulen(mesh,mesh->adjt,kk,ii,n) ) {
          ++nf;
          continue;
        }
        memcpy(pxp->n2,n,3*sizeof(double));

        /* compute tangent as intersection of n1 + n2 */
        ppt->n[0] = pxp->n1[1]*pxp->n2[2] - pxp->n1[2]*pxp->n2[1];
        ppt->n[1] = pxp->n1[2]*pxp->n2[0] - pxp->n1[0]*pxp->n2[2];
        ppt->n[2] = pxp->n1[0]*pxp->n2[1] - pxp->n1[1]*pxp->n2[0];
        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
        if ( dd > _MMG5_EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          ppt->n[0] *= dd;
          ppt->n[1] *= dd;
          ppt->n[2] *= dd;
        }
        ppt->flag = mesh->base;
        ++nt;
        continue;
      }

      /* compute tgte */
      ppt->flag = mesh->base;
      ++nt;
      if ( !_MMG5_boulec(mesh,mesh->adjt,k,i,ppt->n) ) {
        ++nf;
        continue;
      }
      dd = pxp->n1[0]*ppt->n[0] + pxp->n1[1]*ppt->n[1] + pxp->n1[2]*ppt->n[2];
      ppt->n[0] -= dd*pxp->n1[0];
      ppt->n[1] -= dd*pxp->n1[1];
      ppt->n[2] -= dd*pxp->n1[2];
      dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
      if ( dd > _MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        ppt->n[0] *= dd;
        ppt->n[1] *= dd;
        ppt->n[2] *= dd;
      }
    }
  }
  if ( abs(mesh->info.imprim) > 3 && nn+nt > 0 ) {
    if ( nnr )
      fprintf(stdout,"     %d input normals ignored\n",nnr);
    fprintf(stdout,"     %d normals,  %d tangents updated  (%d failed)\n",nn,nt,nf);
  }
  return(1);
}

/** Define continuous geometric support at non manifold vertices, using volume information */
static void _MMG5_nmgeom(MMG5_pMesh mesh){
  MMG5_pTetra     pt;
  MMG5_pPoint     p0;
  MMG5_pxPoint    pxp;
  int        k,base;
  int        *adja;
  double     n[3],t[3];
  char       i,j,ip,ier;

  base = ++mesh->base;
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      if ( adja[i] ) continue;
      for (j=0; j<3; j++) {
        ip = _MMG5_idir[i][j];
        p0 = &mesh->point[pt->v[ip]];
        if ( p0->flag == base )  continue;
        else if ( !(p0->tag & MG_NOM) )  continue;

        p0->flag = base;
        ier = _MMG5_boulenm(mesh,k,ip,i,n,t);

        if ( !ier ) {
          p0->tag |= MG_REQ;
          p0->tag &= ~MG_NOSURF;
          if ( p0->ref != 0 )
            p0->ref = -abs(p0->ref);
          else
            p0->ref = MG_ISO;
        }
        else {
          if ( !p0->xp ) {
            ++mesh->xp;
            if(mesh->xp > mesh->xpmax){
              _MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,0.2,MMG5_xPoint,
                                 "larger xpoint table",
                                 mesh->xp--;
                                 fprintf(stderr,"  Exit program.\n");
                                 exit(EXIT_FAILURE));
            }
            p0->xp = mesh->xp;
          }
          pxp = &mesh->xpoint[p0->xp];
          memcpy(pxp->n1,n,3*sizeof(double));
          memcpy(p0->n,t,3*sizeof(double));
        }
      }
    }
  }
  /* Mark as required the non-manifold points that do not belong to a surface
   * tetra (a tetra that have a face without adjacent)*/
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !(p0->tag & MG_NOM) || p0->xp ) continue;
    p0->tag |= MG_REQ;
    p0->tag &= ~MG_NOSURF;
  }
}

/** preprocessing stage: mesh analysis */
int _MMG3D_analys(MMG5_pMesh mesh) {
  _MMG5_Hash hash;

  /**--- stage 1: data structures for surface */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** SURFACE ANALYSIS\n");

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"  ## Hashing problem (1). Exit program.\n");
    return(0);
  }

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stdout,"  ## Prism hashing problem. Exit program.\n");
    return(0);
  }
  /* compatibility triangle orientation w/r tetras */
  if ( !_MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"  ## Boundary orientation problem. Exit program.\n");
    return(0);
  }

  /* identify surface mesh */
  if ( !_MMG5_chkBdryTria(mesh) ) {
      fprintf(stderr,"  ## Boundary problem. Exit program.\n");
      return(0);
  }
  _MMG5_freeXTets(mesh);
  _MMG5_freeXPrisms(mesh);

  if ( mesh->info.nosurf ) {
    /* Set surface triangles to required*/
    _MMG5_reqBoundaries(mesh);
  }

  /* create surface adjacency */
  if ( !_MMG3D_hashTria(mesh,&hash) ) {
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    fprintf(stderr,"  ## Hashing problem (2). Exit program.\n");
    return(0);
  }

  /* build hash table for geometric edges */
  if ( !_MMG5_hGeom(mesh) ) {
    fprintf(stderr,"  ## Hashing problem (0). Exit program.\n");
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
    return(0);
  }

  /**--- stage 2: surface analysis */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  /* identify connexity */
  if ( !_MMG5_setadj(mesh) ) {
    fprintf(stderr,"  ## Topology problem. Exit program.\n");
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    return(0);
  }

  /* check for ridges */
  if ( mesh->info.dhd > _MMG5_ANGLIM && !_MMG5_setdhd(mesh) ) {
    fprintf(stderr,"  ## Geometry problem. Exit program.\n");
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    return(0);
  }

  /* identify singularities */
  if ( !_MMG5_singul(mesh) ) {
    fprintf(stderr,"  ## MMG5_Singularity problem. Exit program.\n");
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    return(0);
  }

  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

  /* define (and regularize) normals */
  if ( !_MMG5_norver(mesh) ) {
    fprintf(stderr,"  ## Normal problem. Exit program.\n");
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    return(0);
  }

  /* set bdry entities to tetra */
  if ( !_MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"  ## Boundary problem. Exit program.\n");
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    return(0);
  }

  /* set non-manifold edges sharing non-intersecting multidomains as required */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** UPDATING TOPOLOGY AT NON-MANIFOLD POINTS\n");

  if ( !_MMG5_setNmTag(mesh,&hash) ) {
    fprintf(stderr,"  ## Non-manifold topology problem. Exit program.\n");
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    return(0);
  }

  /* build hash table for geometric edges */
  if ( !mesh->na && !_MMG5_hGeom(mesh) ) {
    fprintf(stderr,"  ## Hashing problem (0). Exit program.\n");
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
    return(0);
  }

  /* Update edges tags and references for xtetras */
  if ( !_MMG5_bdryUpdate(mesh) ) {
    fprintf(stderr,"  ## Boundary problem. Exit program.\n");
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    return(0);
  }

  /* define geometry for non manifold points */
  _MMG5_nmgeom(mesh);

  /* release memory */
  _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
  _MMG5_DEL_MEM(mesh,mesh->adjt,(3*mesh->nt+4)*sizeof(int));
  _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));
  mesh->nt = 0;
  if ( mesh->nprism ) _MMG5_DEL_MEM(mesh,mesh->adjapr,(5*mesh->nprism+6)*sizeof(int));

  return(1);
}
