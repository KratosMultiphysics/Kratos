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
 * \file mmg3d/mmg3d1.c
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

char  ddb;

/**
 * \param mesh pointer toward the mesh structure.
 * \param k tetrahedron index.
 * \param ie face index of tetrahedron.
 * \param ptt pointer toward the output triangle.
 *
 * Set triangle corresponding to face ie of tetra k.
 *
 */
void _MMG5_tet2tri(MMG5_pMesh mesh,int k,char ie,MMG5_Tria *ptt) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  char    i;

  pt = &mesh->tetra[k];
  memset(ptt,0,sizeof(MMG5_Tria));
  ptt->v[0] = pt->v[_MMG5_idir[ie][0]];
  ptt->v[1] = pt->v[_MMG5_idir[ie][1]];
  ptt->v[2] = pt->v[_MMG5_idir[ie][2]];
  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    ptt->ref = pxt->ref[ie];
    for (i=0; i<3; i++) {
      ptt->edg[i] = pxt->edg[_MMG5_iarf[ie][i]];
      ptt->tag[i] = pxt->tag[_MMG5_iarf[ie][i]];
    }
  }
  else {
    for (i=0; i<3; i++) {
      ptt->edg[i] = 0;
      ptt->tag[i] = 0;
    }
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k tetrahedron index.
 * \param vx pointer toward table of edges to split.
 * \return 1.
 *
 * Find acceptable position for splitting.
 *
 */
int _MMG3D_dichoto(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx) {
  MMG5_pTetra  pt;
  MMG5_pPoint  pa,pb,ps;
  double       o[6][3],p[6][3];
  float        to,tp,t;
  int          ia,ib,ier,it,maxit;
  char         i;

  pt = &mesh->tetra[k];
  /* get point on surface and along segment for edge split */
  for (i=0; i<6; i++) {
    memset(p[i],0,3*sizeof(double));
    memset(o[i],0,3*sizeof(double));
    if ( vx[i] > 0 ) {
      ia = pt->v[_MMG5_iare[i][0]];
      ib = pt->v[_MMG5_iare[i][1]];
      pa = &mesh->point[ia];
      pb = &mesh->point[ib];
      ps = &mesh->point[vx[i]];
      o[i][0] = 0.5 * (pa->c[0] + pb->c[0]);
      o[i][1] = 0.5 * (pa->c[1] + pb->c[1]);
      o[i][2] = 0.5 * (pa->c[2] + pb->c[2]);
      p[i][0] = ps->c[0];
      p[i][1] = ps->c[1];
      p[i][2] = ps->c[2];
    }
  }
  maxit = 4;
  it = 0;
  tp = 1.0;
  to = 0.0;
  do {
    /* compute new position */
    t = 0.5 * (tp + to);
    for (i=0; i<6; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
        ps->c[2] = o[i][2] + t*(p[i][2] - o[i][2]);
      }
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32:
      ier = _MMG3D_split1_sim(mesh,met,k,vx);
      break;
    case 11: case 21: case 38: case 56:
      ier = _MMG3D_split3_sim(mesh,met,k,vx);
      break;
    default:
      ier = _MMG5_split2sf_sim(mesh,met,k,vx);
      break;
    }
    if ( ier )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );
  /* restore coords of last valid pos. */
  if ( !ier ) {
    t = to;
    for (i=0; i<6; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
        ps->c[2] = o[i][2] + t*(p[i][2] - o[i][2]);
      }
    }
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param list pointer toward the shell of edge.
 * \param ret double of the number of tetrahedra in the shell.
 * \param ip new point index.
 * \return 1.
 *
 * Find acceptable position for _MMG5_split1b, passing the shell of
 * considered edge, starting from o point.
 *
 */
int _MMG3D_dichoto1b(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ret,int ip) {
  MMG5_pTetra  pt;
  MMG5_pPoint  p0,p1,ppt;
  int          iel,np,nq,it,maxit;
  double       m[3],o[3],tp,to,t;
  char         ia,ier;

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  np = pt->v[_MMG5_iare[ia][0]];
  nq = pt->v[_MMG5_iare[ia][1]];
  p0 = &mesh->point[np];
  p1 = &mesh->point[nq];

  /* initial coordinates for new point */
  ppt = &mesh->point[ip];
  o[0] = ppt->c[0];
  o[1] = ppt->c[1];
  o[2] = ppt->c[2];

  /* midpoint along edge */
  m[0] = 0.5*(p0->c[0] + p1->c[0]);
  m[1] = 0.5*(p0->c[1] + p1->c[1]);
  m[2] = 0.5*(p0->c[2] + p1->c[2]);

  maxit = 4;
  it    = 0;
  ier   = 0;
  tp    = 1.0;
  to    = 0.0;
  do {
    t = 0.5*(to + tp);
    ppt->c[0] = m[0] + t*(o[0]-m[0]);
    ppt->c[1] = m[1] + t*(o[1]-m[1]);
    ppt->c[2] = m[2] + t*(o[2]-m[2]);

    ier = _MMG3D_simbulgept(mesh,met,list,ret,ip);
    if ( ier )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );
  if ( !ier )  t = to;

  ppt->c[0] = m[0] + t*(o[0]-m[0]);
  ppt->c[1] = m[1] + t*(o[1]-m[1]);
  ppt->c[2] = m[2] + t*(o[2]-m[2]);

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param pt pointer toward the triangle.
 * \param ori orientation of the triangle (1 for direct orientation, 0 otherwise).
 * \param hmax maximal edge length.
 * \param hausd maximal hausdorff distance.
 * \param locPar 1 if hmax and hausd are locals parameters.
 * \return 0 if error.
 * \return edges of the triangle pt that need to be split.
 *
 * Find edges of (virtual) triangle pt that need to be split with
 * respect to the Hausdorff criterion.
 *
 */
char _MMG5_chkedg(MMG5_pMesh mesh,MMG5_Tria *pt,char ori, double hmax,
                  double hausd, int locPar) {
  MMG5_pPoint   p[3];
  MMG5_xPoint  *pxp;
//  MMG5_pPar     par;
  double   n[3][3],t[3][3],nt[3],*n1,*n2,t1[3],t2[3];
  double   ps,ps2,ux,uy,uz,ll,il,alpha,dis,hma2;
  int      ia,ib,ic;//l,info;
  char     i,i1,i2;

  ia   = pt->v[0];
  ib   = pt->v[1];
  ic   = pt->v[2];
  p[0] = &mesh->point[ia];
  p[1] = &mesh->point[ib];
  p[2] = &mesh->point[ic];
  pt->flag = 0;

  /* normal recovery */
  for (i=0; i<3; i++) {
    if ( MG_SIN(p[i]->tag) ) {
      _MMG5_nortri(mesh,pt,n[i]);
      if(!ori) {
        n[i][0] *= -1.0;
        n[i][1] *= -1.0;
        n[i][2] *= -1.0;
      }
    }
    else if (p[i]->tag & MG_NOM){
      _MMG5_nortri(mesh,pt,n[i]);
      if(!ori) {
        n[i][0] *= -1.0;
        n[i][1] *= -1.0;
        n[i][2] *= -1.0;
      }
      assert(p[i]->xp);
      pxp = &mesh->xpoint[p[i]->xp];
      memcpy(&t[i],p[i]->n,3*sizeof(double));
    }
    else {
      assert(p[i]->xp);
      pxp = &mesh->xpoint[p[i]->xp];
      if ( MG_EDG(p[i]->tag) ) {
        memcpy(&t[i],p[i]->n,3*sizeof(double));
        _MMG5_nortri(mesh,pt,nt);
        if(!ori) {
          nt[0] *= -1.0;
          nt[1] *= -1.0;
          nt[2] *= -1.0;
        }
        ps  = pxp->n1[0]*nt[0] + pxp->n1[1]*nt[1] + pxp->n1[2]*nt[2];
        ps2 = pxp->n2[0]*nt[0] + pxp->n2[1]*nt[1] + pxp->n2[2]*nt[2];
        if ( fabs(ps) > fabs(ps2) )
          memcpy(&n[i],pxp->n1,3*sizeof(double));
        else
          memcpy(&n[i],pxp->n2,3*sizeof(double));
      }
      else
        memcpy(&n[i],pxp->n1,3*sizeof(double));
    }
  }

  /* analyze edges */
  for (i=0; i<3; i++) {
    i1 = _MMG5_inxt2[i];
    i2 = _MMG5_iprv2[i];

    /* local parameters at vertices */
    /* if ( mesh->info.parTyp & MG_Vert ) { */

    /*   if ( p[i1]->ref == p[i2]->ref ) { */
    /*     for ( l=0; l<mesh->info.npar; ++l) { */
    /*       par = &mesh->info.par[l]; */

    /*       if ( par->elt != MMG5_Vertex || par->ref != p[i1]->ref ) continue; */

    /*       if ( !locPar ) { */
    /*         hausd   = par->hausd; */
    /*         hmax    = par->hmax; */
    /*       } */
    /*       else { */
    /*         hausd = MG_MIN(hausd, par->hausd); */
    /*         hmax = MG_MIN(hmax, par->hmax); */
    /*       } */
    /*       break; */
    /*     } */
    /*   } */
    /*   else { */
    /*     l = 0; */
    /*     info = -1000; */
    /*     do { */
    /*       if ( info >= 0 ) break; */

    /*       par = &mesh->info.par[l]; */
    /*       if ( par->elt != MMG5_Vertex ) continue; */

    /*       if (p[i1]->ref != par->ref && p[i2]->ref != par->ref ) continue; */

    /*       if ( !locPar ) { */
    /*         hausd   = par->hausd; */
    /*         hmin    = par->hmin; */
    /*         hmax    = par->hmax; */
    /*         locPar   = 1; */
    /*       } */
    /*       else { */
    /*         hausd   = MG_MIN(hausd,par->hausd); */
    /*         hmin    = MG_MAX(hmin,par->hmin); */
    /*         hmax    = MG_MIN(hmax,par->hmax); */
    /*       } */

    /*       info    = par->ref; */
    /*     } while ( ++l<mesh->info.npar ); */

    /*     for ( ; l<mesh->info.npar; ++l) { */
    /*       par = &mesh->info.par[l]; */
    /*       if ( par->elt != MMG5_Vertex || par->ref == info ) continue; */
    /*       if (p[i1]->ref != par->ref && p[i2]->ref != par->ref ) continue; */

    /*       hausd   = MG_MIN(hausd,par->hausd); */
    /*       hmax    = MG_MIN(hmax,par->hmax); */
    /*       break; */
    /*     } */
    /*   } */
    /* } */

    hma2 = _MMG5_LLONG*_MMG5_LLONG*hmax*hmax;

    /* check length */
    ux = p[i2]->c[0] - p[i1]->c[0];
    uy = p[i2]->c[1] - p[i1]->c[1];
    uz = p[i2]->c[2] - p[i1]->c[2];
    ll = ux*ux + uy*uy + uz*uz;
    if ( ll < _MMG5_EPSD )  continue;
    else if ( ll > hma2 ) {
      MG_SET(pt->flag,i);
      continue;
    }
    il = 1.0 / sqrt(ll);

    /* Hausdorff w/r tangent direction */
    if ( MG_EDG(pt->tag[i]) || ( pt->tag[i] & MG_NOM )) {
      if ( MG_SIN(p[i1]->tag) ) {
        t1[0] = il * ux;
        t1[1] = il * uy;
        t1[2] = il * uz;
      }
      else {
        if(!((p[i1]->tag & MG_NOM) ||  MG_EDG(p[i1]->tag) ) ) {
          fprintf(stdout,"2. warning geometrical problem\n");
          return(0);
        }
        memcpy(t1,t[i1],3*sizeof(double));
        ps = t1[0]*ux + t1[1]*uy + t1[2]*uz;
        if ( ps < 0.0 ) {
          t1[0] *= -1.0;
          t1[1] *= -1.0;
          t1[2] *= -1.0;
        }
      }
      if ( MG_SIN(p[i2]->tag) ) {
        t2[0] = -il * ux;
        t2[1] = -il * uy;
        t2[2] = -il * uz;
      }
      else {
        if(!((p[i2]->tag & MG_NOM) || MG_EDG(p[i2]->tag) ) ) {
          fprintf(stdout,"2. warning geometrical problem\n");
          return(0);
        }
        memcpy(t2,t[i2],3*sizeof(double));
        ps = - ( t2[0]*ux + t2[1]*uy + t2[2]*uz );
        if ( ps < 0.0 ) {
          t2[0] *= -1.0;
          t2[1] *= -1.0;
          t2[2] *= -1.0;
        }
      }
    }
    else {
      n1 = n[i1];
      n2 = n[i2];
      if ( !_MMG5_BezierTgt(p[i1]->c,p[i2]->c,n1,n2,t1,t2) ) {
        t1[0] = ux * il;
        t1[1] = uy * il;
        t1[2] = uz * il;

        t2[0] = -ux * il;
        t2[1] = -uy * il;
        t2[2] = -uz * il;
      }
    }
    alpha = _MMG5_BezierGeod(p[i1]->c,p[i2]->c,t1,t2);
    ps  = t1[0]*ux + t1[1]*uy + t1[2]*uz;
    ps *= il;
    dis = alpha*alpha*fabs(1.0 - ps*ps);
    if ( dis > hausd*hausd ) {
      MG_SET(pt->flag,i);
      continue;
    }
    ps  = -( t2[0]*ux + t2[1]*uy + t2[2]*uz );
    ps *= il;
    dis = alpha*alpha*fabs(1.0 - ps*ps);

    if ( dis > hausd*hausd ) {
      MG_SET(pt->flag,i);
      continue;
    }
  }
  return(pt->flag);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure (only for delaunay).
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 * \return -1 if failed and swap number otherwise.
 *
 * Search for boundary edges that could be swapped for geometric
 * approximation.
 *
 */
int _MMG5_swpmsh(MMG5_pMesh mesh,MMG5_pSol met,_MMG3D_pOctree octree, int typchk) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  int      k,it,list[MMG3D_LMAX+2],ilist,ret,it1,it2,ns,nns,maxit;
  char     i,j,ia,ier;

  it = nns = 0;
  maxit = 2;
  do {
    ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( (!MG_EOK(pt)) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
      else if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];

      for (i=0; i<4; i++) {
        ier = 0;
        if ( !(pxt->ftag[i] & MG_BDY) ) continue;
        for (j=0; j<3; j++) {
          ia  = _MMG5_iarf[i][j];

          /* No swap of geometric edge */
          if ( MG_EDG(pxt->tag[ia]) || (pxt->tag[ia] & MG_REQ) ||
               (pxt->tag[ia] & MG_NOM) )
            continue;

          ret = _MMG5_coquilface(mesh,k,ia,list,&it1,&it2,0);
          ilist = ret / 2;
          if ( ret < 0 )  return(-1);
          /* CAUTION: trigger collapse with 2 elements */
          if ( ilist <= 1 )  continue;
          ier = _MMG5_chkswpbdy(mesh,met,list,ilist,it1,it2,typchk);
          if ( ier ) {
            ier = _MMG5_swpbdy(mesh,met,list,ret,it1,octree,typchk);
            if ( ier > 0 )  ns++;
            else if ( ier < 0 )  return(-1);
            break;
          }
        }
        if ( ier )  break;
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return(nns);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param crit coefficient of quality improvment.
 * \param octree pointer toward the octree structure in delaunay mode and
 * toward the \a NULL pointer otherwise
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 *
 * Internal edge flipping.
 *
 */
int _MMG5_swptet(MMG5_pMesh mesh,MMG5_pSol met,double crit,
                 _MMG3D_pOctree octree,int typchk) {
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
      if ( pt->qual > 0.0288675 /*0.6/_MMG5_ALPHAD*/ )  continue;

      for (i=0; i<6; i++) {
        /* Prevent swap of a ref or tagged edge */
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          if ( pxt->edg[i] || pxt->tag[i] ) continue;
        }

        nconf = _MMG5_chkswpgen(mesh,met,k,i,&ilist,list,crit,typchk);
        if ( nconf ) {
          ier = _MMG5_swpgen(mesh,met,nconf,ilist,list,octree,typchk);
          if ( ier > 0 )  ns++;
          else if ( ier < 0 ) return(-1);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return(nns);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param maxitin maximum number of iteration.
 * \return -1 if failed, number of moved points otherwise.
 *
 * Analyze tetrahedra and move points so as to make mesh more uniform.
 * In delaunay mode, a negative maxitin means that we don't move internal nodes.
 *
 */
int _MMG5_movtet(MMG5_pMesh mesh,MMG5_pSol met, _MMG3D_pOctree octree,int maxitin) {
  MMG5_pTetra        pt;
  MMG5_pPoint        ppt;
  MMG5_pxTetra       pxt;
  double        *n;
  int           i,k,ier,nm,nnm,ns,lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv,it;
  int           improve;
  unsigned char j,i0,base;
  int           internal,maxit;

  if ( maxitin<0 ) {
    internal = 0;
    maxit = abs(maxitin);
  } else {
    internal=1;
    maxit = maxitin;
  }
  if ( maxit != 1 ) {
    improve   = 1;
  } else {
    improve = 0;
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** OPTIMIZING MESH\n");

  base = 1;
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = base;

  it = nnm = 0;
  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;

      /* point j on face i */
      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
            if ( pxt->tag[_MMG5_iarf[i][j]] & MG_REQ )  continue;
          }
          else  pxt = 0;
          i0  = _MMG5_idir[i][j];
          ppt = &mesh->point[pt->v[i0]];
          if ( ppt->flag == base )  continue;
          else if ( MG_SIN(ppt->tag) )  continue;

          if ( maxit != 1 ) {
            ppt->flag = base;
          }
          ier = 0;
          if ( ppt->tag & MG_BDY ) {
            /* Catch a boundary point by a boundary face */
            if ( !pt->xt || !(MG_BDY & pxt->ftag[i]) )  continue;
            else if( ppt->tag & MG_NOM ){
              if( mesh->adja[4*(k-1)+1+i] ) continue;
              ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,1);
              if( !ier )  continue;
              else if ( ier>0 )
                ier = _MMG5_movbdynompt(mesh,met,octree,listv,ilistv,lists,ilists,improve);
              else
                return(-1);
            }
            else if ( ppt->tag & MG_GEO ) {
              ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
              if ( !ier )  continue;
              else if ( ier>0 )
                ier = _MMG5_movbdyridpt(mesh,met,octree,listv,ilistv,lists,ilists,improve);
              else
                return(-1);
            }
            else if ( ppt->tag & MG_REF ) {
              ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
              if ( !ier )
                continue;
              else if ( ier>0 )
                ier = _MMG5_movbdyrefpt(mesh,met,octree,listv,ilistv,lists,ilists,improve);
              else
                return(-1);
            }
            else {
              ier=_MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
              if ( !ier )
                continue;
              else if ( ier<0 )
                return(-1);

              n = &(mesh->xpoint[ppt->xp].n1[0]);
              // if ( MG_GET(pxt->ori,i) ) {
              /* Useless because if the orientation of the tetra face is
               * compatible with the triangle (MG_GET(ori,i)) we know that we
               * are well orientated. Morever, may introduce numerical errors
               * with wrinkled surfaces. */
                // if ( !_MMG5_directsurfball(mesh, pt->v[i0],lists,ilists,n) )  continue;
              // }
              if ( !MG_GET(pxt->ori,i) ) {
                if ( !_MMG5_directsurfball(mesh,pt->v[i0],lists,ilists,n) )
                  continue;
              }
              ier = _MMG5_movbdyregpt(mesh,met, octree, listv,ilistv,lists,ilists,improve);
              if ( ier )  ns++;
            }
          }
          else if ( internal ) {
            ilistv = _MMG5_boulevolp(mesh,k,i0,listv);
            if ( !ilistv )  continue;
            ier = _MMG5_movintpt(mesh,met,octree,listv,ilistv,improve);
          }
          if ( ier ) {
            nm++;
            if(maxit==1){
              ppt->flag = base;
            }
          }
        }
      }
    }
    nnm += nm;
    if ( mesh->info.ddebug )  fprintf(stdout,"     %8d moved, %d geometry\n",nm,ns);
  }
  while( ++it < maxit && nm > 0 );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm )
    fprintf(stdout,"     %8d vertices moved, %d iter.\n",nnm,it);

  return(nnm);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking permformed for edge length (hmin or LSHORT criterion).
 * \return -1 if failed, number of collapsed points otherwise.
 *
 * Attempt to collapse small edges.
 *
 */
static int _MMG5_coltet(MMG5_pMesh mesh,MMG5_pSol met,char typchk) {
  MMG5_pTetra     pt,ptloc;
  MMG5_pxTetra    pxt;
  MMG5_pPoint     p0,p1;
  MMG5_pPar       par;
  double     ll,ux,uy,uz,hmi2;
  int        k,nc,list[MMG3D_LMAX+2],ilist,ilists,lists[MMG3D_LMAX+2];
  int        base,nnm,l,kk,isloc,ifac1;
  int16_t    tag,isnm;
  char       i,j,ip,iq;
  int        ier;

  nc = nnm = 0;

  /* init point flags */
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    p0->flag = 0;
  }

  for (k=1; k<=mesh->ne; k++) {
    base = ++mesh->base;
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )   continue;

    pxt = pt->xt?  &mesh->xtetra[pt->xt] : 0;

    for (i=0; i<4; i++) {
      ier = 0;

      for (j=0; j<3; j++) {
        if ( pt->xt && (pxt->tag[_MMG5_iarf[i][j]] & MG_REQ) )  continue;
        ip = _MMG5_idir[i][_MMG5_inxt2[j]];
        iq = _MMG5_idir[i][_MMG5_iprv2[j]];

        p0 = &mesh->point[pt->v[ip]];
        p1 = &mesh->point[pt->v[iq]];
        if ( p0->flag == base )  continue;
        else if ( (p0->tag & MG_REQ) || (p0->tag > p1->tag) )  continue;


        /* Ball of point: computed here if needed for the local parameter
         * evaluation, after length check otherwise (because the ball
         * computation is time consuming) */
        if ( mesh->info.npar ) {
          if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
            tag = pxt->tag[_MMG5_iarf[i][j]];
            isnm = (tag & MG_NOM);

            if ( p0->tag > tag ) continue;
            if ( isnm && mesh->adja[4*(k-1)+1+i] )  continue;
            if (_MMG5_boulesurfvolp(mesh,k,ip,i,
                                    list,&ilist,lists,&ilists,p0->tag & MG_NOM) < 0 )
              return(-1);
          }
          else {
            ilist = _MMG5_boulevolp(mesh,k,ip,list);
          }
        }

        /* check length */
        if ( typchk == 1 ) {
          ux = p1->c[0] - p0->c[0];
          uy = p1->c[1] - p0->c[1];
          uz = p1->c[2] - p0->c[2];
          ll = ux*ux + uy*uy + uz*uz;

          /* local parameters*/
          hmi2  = mesh->info.hmin;
          isloc = 0;

          /* Local parameters at vertex */
          // take the min of the local hmin at p0 and hmin at p1

          /* Local parameters at tetra */
          if ( mesh->info.parTyp & MG_Tetra ) {
            l = 0;
            do
            {
              if ( isloc )  break;

              par = &mesh->info.par[l];
              if ( par->elt != MMG5_Tetrahedron )  continue;

              for ( kk=0; kk<ilist; ++kk ) {
                ptloc = &mesh->tetra[list[kk]/4];
                if ( par->ref == ptloc->ref ) {
                  hmi2 = par->hmin;
                  isloc   = 1;
                  break;
                }
              }
            } while ( ++l<mesh->info.npar );

            for ( ; l<mesh->info.npar; ++l ) {
              par = &mesh->info.par[l];
              if ( par->elt != MMG5_Tetrahedron ) continue;

              for ( kk=0; kk<ilist; ++kk ) {
                ptloc = &mesh->tetra[list[kk]/4];
                if ( par->ref == ptloc->ref ) {
                  hmi2 = MG_MAX(hmi2,par->hmin);
                  break;
                }
              }
            }
          }

          /* Local parameters at triangle */
          if ( mesh->info.parTyp & MG_Tria && ( pt->xt && (pxt->ftag[i] & MG_BDY) )) {
            l = 0;
            do
            {
              if ( isloc )  break;

              par = &mesh->info.par[l];
              if ( par->elt != MMG5_Triangle )  continue;

              for ( kk=0; kk<ilists; ++kk ) {
                ptloc = &mesh->tetra[lists[kk]/4];
                ifac1 =  lists[kk] % 4;
                assert(ptloc->xt && (mesh->xtetra[ptloc->xt].ftag[ifac1] & MG_BDY) );

                if ( par->ref == mesh->xtetra[ptloc->xt].ref[ifac1] ) {
                  hmi2  = par->hmin;
                  isloc = 1;
                  break;
                }
              }
            } while ( ++l<mesh->info.npar );

            for ( ; l<mesh->info.npar; ++l ) {
              par = &mesh->info.par[l];
              if ( par->elt != MMG5_Triangle ) continue;

              for ( kk=0; kk<ilists; ++kk ) {
                ptloc = &mesh->tetra[lists[kk]/4];
                ifac1 =  lists[kk] % 4;
                assert(ptloc->xt && (mesh->xtetra[ptloc->xt].ftag[ifac1] & MG_BDY) );

                if ( par->ref == mesh->xtetra[ptloc->xt].ref[ifac1] ) {
                  hmi2  = MG_MAX(hmi2,par->hmin);
                  break;
                }
              }
            }
          }

          hmi2 = hmi2*hmi2;
          if ( ll > hmi2 )  continue;
        }
        else if ( typchk == 2 ) {
          ll = _MMG5_lenedg(mesh,met,_MMG5_iarf[i][j],pt);
          // Case of an internal tetra with 4 ridges vertices.
          if ( ll == 0 ) continue;
          if ( ll > _MMG5_LSHRT )  continue;
        }


        /* Ball computation if not computed before */
        if ( !mesh->info.npar ) {
          if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
            tag = pxt->tag[_MMG5_iarf[i][j]];
            isnm = (tag & MG_NOM);

            if ( p0->tag > tag ) continue;
            if ( isnm && mesh->adja[4*(k-1)+1+i] )  continue;
            if (_MMG5_boulesurfvolp(mesh,k,ip,i,
                                    list,&ilist,lists,&ilists,p0->tag & MG_NOM) < 0 )
              return(-1);
          }
          else {
            ilist = _MMG5_boulevolp(mesh,k,ip,list);
          }
        }

        /* boundary face: collapse ip on iq */
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
          tag = pxt->tag[_MMG5_iarf[i][j]];
          tag |= MG_BDY;

          isnm = ( tag & MG_NOM );
          if ( isnm ) {
            if ( mesh->adja[4*(k-1)+1+i] )  continue;
          }
          if ( p0->tag > tag )  continue;
          ilist = _MMG5_chkcol_bdy(mesh,met,k,i,j,list,ilist,lists,ilists,typchk);
        }
        /* internal face */
        else {
          isnm = 0;
          if ( p0->tag & MG_BDY )  continue;
          ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,typchk);
        }

        if ( ilist > 0 ) {
          ier = _MMG5_colver(mesh,met,list,ilist,iq,typchk);
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
  if ( nc > 0 && (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) )
    fprintf(stdout,"     %8d vertices removed, %8d non manifold,\n",nc,nnm);

  return(nc);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking permformed for edge length (hmax or _MMG5_LLONG criterion).
 * \return -1 if failed.
 * \return number of new points.
 *
 * Analyze volume tetra and split if needed.
 *
 */
static int
_MMG5_anatetv(MMG5_pMesh mesh,MMG5_pSol met,char typchk) {
  MMG5_pTetra   pt;
  MMG5_pPoint   p1,p2;
  MMG5_xTetra  *pxt;
  _MMG5_Hash    hash;
  MMG5_pPar     par;
  double   ll,o[3],ux,uy,uz,hma2;
  int      l,vx[6],k,ip,ip1,ip2,nap,ns,ne,memlack,ier;
  char     i,j,ia;

  /** 1. analysis */
  if ( !_MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return(-1);
  memlack = ns = nap = 0;

  /* Hash all boundary and required edges, and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (i=0; i<6; i++) {
        ip1 = pt->v[_MMG5_iare[i][0]];
        ip2 = pt->v[_MMG5_iare[i][1]];
        ip  = -1;
        if ( !_MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
      }
      continue;
    }

    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] & MG_BDY ) {
        for (j=0; j<3; j++) {
          ip1 = pt->v[_MMG5_idir[i][_MMG5_inxt2[j]]];
          ip2 = pt->v[_MMG5_idir[i][_MMG5_iprv2[j]]];
          ip  = -1;
          if ( !_MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
        }
        break;
      }
    }
  }

  /** 2. Set flags and split internal edges */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = 0;
    for (i=0; i<6; i++) {
      ip  = -1;
      ip1 = pt->v[_MMG5_iare[i][0]];
      ip2 = pt->v[_MMG5_iare[i][1]];
      p1  = &mesh->point[ip1];
      p2  = &mesh->point[ip2];
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( pxt->tag[i] & MG_REQ ) continue;
      }
      else  pxt = 0;
      if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) ) {
        ip = _MMG5_hashGet(&hash,ip1,ip2);
      }
      else {
        if (typchk == 1) {
          ux = p2->c[0] - p1->c[0];
          uy = p2->c[1] - p1->c[1];
          uz = p2->c[2] - p1->c[2];
          ll = ux*ux + uy*uy + uz*uz;

          hma2 = mesh->info.hmax;
          /* Local parameters */
          /* Local param at vertices ip1 and ip2 */
          // take the min of the local hmin at ip1 and hmin at ip2

          /* Local parameters at tetra */
          if ( mesh->info.parTyp & MG_Tetra ) {
            for ( l=0; l<mesh->info.npar; ++l ) {
              par = &mesh->info.par[l];

              if ( par->elt != MMG5_Tetrahedron )  continue;
              if ( par->ref != pt->ref ) continue;

              hma2  = par->hmax;
              break;
            }
          }
          hma2 = _MMG5_LLONG*_MMG5_LLONG*hma2*hma2;
          if ( ll > hma2 )
            ip = _MMG5_hashGet(&hash,ip1,ip2);
        }
        else if ( typchk == 2 ) {
          ll = _MMG5_lenedg(mesh,met,i,pt);
          // Case of an internal tetra with 4 ridges vertices.
          if ( ll == 0 ) continue;
          if ( ll > _MMG5_LLONG )
            ip = _MMG5_hashGet(&hash,ip1,ip2);
        }
      }
      if ( ip < 0 ) continue;
      else if ( !ip ) {
        /* new midpoint */
        o[0] = 0.5 * (p1->c[0]+p2->c[0]);
        o[1] = 0.5 * (p1->c[1]+p2->c[1]);
        o[2] = 0.5 * (p1->c[2]+p2->c[2]);

        ip  = _MMG3D_newPt(mesh,o,0);
        if ( !ip ) {
          /* reallocation of point table */

          _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                              printf("  ## Warning: unable to allocate a new point\n");
                              _MMG5_INCREASE_MEM_MESSAGE();
                              memlack=1;
                              goto split
                              ,o,0);
          p1  = &mesh->point[ip1];
          p2  = &mesh->point[ip2];
        }

        if ( met->m ) {
          if ( typchk == 1 && (met->size>1) )
            ier = _MMG3D_intmet33_ani(mesh,met,k,i,ip,0.5);
          else
            ier = _MMG5_intmet(mesh,met,k,i,ip,0.5);

          if (!ier) {
            // Unable to compute the metric
            return(-1);
          }
          else if ( ier < 0 ) {
            _MMG3D_delPt(mesh,ip);
            continue;
          }
        }

        if ( !_MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
        MG_SET(pt->flag,i);
        nap++;
      }
    }
  }
  if ( !nap )  {
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    return(0);
  }

  /** 3. check and split */
split:
  ns = 0;
  ne = mesh->ne;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    memset(vx,0,6*sizeof(int));
    pt->flag = 0;
    for (ia=0,i=0; i<3; i++) {
      for (j=i+1; j<4; j++,ia++) {
        if ( pt->xt && (mesh->xtetra[pt->xt].tag[ia] & MG_REQ) ) continue;
        vx[ia] = _MMG5_hashGet(&hash,pt->v[i],pt->v[j]);
        if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
      }
    }

    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      _MMG5_split1(mesh,met,k,vx,typchk-1);
      ns++;
      break;
    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      _MMG5_split2sf(mesh,met,k,vx,typchk-1);
      ns++;
      break;

    case 33: case 18: case 12: /* 2 opposite edges split */
      _MMG5_split2(mesh,met,k,vx,typchk-1);
      ns++;
      break;

    case 11: case 21: case 38: case 56: /* 3 edges on the same faces splitted */
      _MMG5_split3(mesh,met,k,vx,typchk-1);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      _MMG5_split3cone(mesh,met,k,vx,typchk-1);
      ns++;
      break;

    case 35: case 19: case 13: case 37: case 22: case 28: case 26:
    case 14: case 49: case 50: case 44: case 41: /* 3 edges on opposite configuration splitted */
      _MMG5_split3op(mesh,met,k,vx,typchk-1);
      ns++;
      break;

    case 23: case 29: case 53: case 60: case 57: case 58:
    case 27: case 15: case 43: case 39: case 54: case 46: /* 4 edges with 3 lying on the same face splitted */
      _MMG5_split4sf(mesh,met,k,vx,typchk-1);
      ns++;
      break;

      /* 4 edges with no 3 lying on the same face splitted */
    case 30: case 45: case 51:
      _MMG5_split4op(mesh,met,k,vx,typchk-1);
      ns++;
      break;

    case 62: case 61: case 59: case 55: case 47: case 31: /* 5 edges split */
      _MMG5_split5(mesh,met,k,vx,typchk-1);
      ns++;
      break;

    case 63: /* 6 edges split */
      _MMG5_split6(mesh,met,k,vx,typchk-1);
      ns++;
      break;
    }
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",nap);

  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  if ( memlack )  return(-1);
  return(nap);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking permformed for edge length (hmax or _MMG5_LLONG criterion).
 * \return -1 if failed.
 * \return number of new points.
 *
 * Analyze tetra and split on geometric criterion.
 *
 */
static int
_MMG5_anatets(MMG5_pMesh mesh,MMG5_pSol met,char typchk) {
  MMG5_pTetra   pt;
  MMG5_pPoint   ppt,p1,p2;
  MMG5_Tria     ptt;
  MMG5_xTetra  *pxt;
  MMG5_xPoint  *pxp;
  _MMG5_Bezier  pb;
  _MMG5_Hash    hash;
  MMG5_pPar     par;
  double        o[3],no[3],to[3],dd,len,hmax,hausd;
  int           vx[6],k,l,ip,ic,it,nap,nc,ni,ne,npinit,ns,ip1,ip2,ier,isloc;
  char          i,j,ia,i1,i2;
  static double uv[3][2] = { {0.5,0.5}, {0.,0.5}, {0.5,0.} };

  /** 1. analysis of boundary elements */
  if ( !_MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) ) return(-1);
  ns = nap = 0;
  npinit=mesh->np;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) || !pt->xt )  continue;

    /* check boundary face cut w/r Hausdorff or hmax */
    pt->flag = 0;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++){
      if ( pxt->ftag[i] & MG_REQ )  continue;
      if ( pxt->ftag[i] & MG_BDY )  break;
    }
    if ( i == 4 )  continue;

    /* virtual triangle */
    _MMG5_tet2tri(mesh,k,i,&ptt);
    if ( typchk == 1 ) {
      if ( !MG_GET(pxt->ori,i) ) continue;

      /* Local parameters for ptt and k */
      hmax  = mesh->info.hmax;
      hausd = mesh->info.hausd;
      isloc = 0;

      if ( mesh->info.parTyp & MG_Tetra ) {
        for ( l=0; l<mesh->info.npar; ++l ) {
          par = &mesh->info.par[l];

          if ( par->elt != MMG5_Tetrahedron )  continue;
          if ( par->ref != pt->ref ) continue;

          hmax = par->hmax;
          hausd = par->hausd;
          isloc = 1;
          break;
        }
      }
      if ( mesh->info.parTyp & MG_Tria ) {
        if ( isloc ) {
          for ( l=0; l<mesh->info.npar; ++l ) {
            par = &mesh->info.par[l];

            if ( par->elt != MMG5_Triangle )  continue;
            if ( par->ref != ptt.ref ) continue;

            hmax = MG_MIN(hmax,par->hmax);
            hausd = MG_MIN(hausd,par->hausd);
            break;
          }
        }
        else {
          for ( l=0; l<mesh->info.npar; ++l ) {
            par = &mesh->info.par[l];

            if ( par->elt != MMG5_Triangle )  continue;
            if ( par->ref != ptt.ref ) continue;

            hmax  = par->hmax;
            hausd = par->hausd;
            isloc = 1;
            break;
          }
        }
      }

      if ( !_MMG5_chkedg(mesh,&ptt,MG_GET(pxt->ori,i),hmax,hausd,isloc) )
        continue;

      /* put back flag on tetra */
      for (j=0; j<3; j++){
        if ( pxt->tag[_MMG5_iarf[i][j]] & MG_REQ )  continue;
        if ( MG_GET(ptt.flag,j) )  MG_SET(pt->flag,_MMG5_iarf[i][j]);
      }
    }
    else if ( typchk == 2 ) {
      for (j=0; j<3; j++) {
        ia = _MMG5_iarf[i][j];
        if ( pxt->tag[ia] & MG_REQ )  continue;
        i1  = _MMG5_iare[ia][0];
        i2  = _MMG5_iare[ia][1];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        len = _MMG5_lenedg(mesh,met,ia,pt);
        // Case of an internal tetra with 4 ridges vertices.
        if ( len == 0 ) continue;
        if ( len > _MMG5_LLONG )  MG_SET(pt->flag,ia);
        /* Treat here the ridges coming from a corner (we can not do that after
         * because the corner don't have xpoints) */
        if ( (mesh->point[ip1].tag & MG_CRN) ||  (mesh->point[ip2].tag & MG_CRN) ) {
          if ( len > _MMG5_LOPTL )  MG_SET(pt->flag,ia);
        }
      }
    }
    if ( !pt->flag )  continue;
    ns++;

    /* geometric support */
    ier = _MMG5_bezierCP(mesh,&ptt,&pb,MG_GET(pxt->ori,i));
    assert(ier);

    /* scan edges in face to split */
    for (j=0; j<3; j++) {
      ia = _MMG5_iarf[i][j];
      if ( !MG_GET(pt->flag,ia) )  continue;
      if ( pxt->tag[ia] & MG_REQ )  continue;
      i1  = _MMG5_iare[ia][0];
      i2  = _MMG5_iare[ia][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      ip  = _MMG5_hashGet(&hash,ip1,ip2);
      if ( ip > 0 && !(ptt.tag[j] & MG_GEO) )  continue;

      ier = _MMG3D_bezierInt(&pb,&uv[j][0],o,no,to);
      /* new point along edge */
      if ( !ip ) {
        ip = _MMG3D_newPt(mesh,o,MG_BDY);
        if ( !ip ) {
          /* reallocation of point table */
          _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                              fprintf(stderr,"  ## Error: unable to allocate a new point.\n");
                              _MMG5_INCREASE_MEM_MESSAGE();
                              do {
                                _MMG3D_delPt(mesh,mesh->np);
                              } while ( mesh->np>npinit );
                              return(-1)
                              ,o,MG_BDY);
          // Now pb->p contain a wrong memory address.
          pb.p[0] = &mesh->point[ptt.v[0]];
          pb.p[1] = &mesh->point[ptt.v[1]];
          pb.p[2] = &mesh->point[ptt.v[2]];
        }
        if ( !_MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return(-1);
        ppt = &mesh->point[ip];
        p1  = &mesh->point[ip1];
        p2  = &mesh->point[ip2];

        if ( met->m ) {
          if ( typchk == 1 && (met->size>1) )
            ier = _MMG3D_intmet33_ani(mesh,met,k,ia,ip,0.5);
          else
            ier = _MMG5_intmet(mesh,met,k,ia,ip,0.5);

          if ( !ier )  return(-1);
          else if ( ier < 0 ) {
            _MMG3D_delPt(mesh,ip);
            continue;
          }
        }

        if ( MG_EDG(ptt.tag[j]) || (ptt.tag[j] & MG_NOM) )
          ppt->ref = ptt.edg[j] ? ptt.edg[j] : ptt.ref;
        else
          ppt->ref = ptt.ref;
        ppt->tag |= ptt.tag[j];
        pxp = &mesh->xpoint[ppt->xp];
        memcpy(pxp->n1,no,3*sizeof(double));
        memcpy(ppt->n,to,3*sizeof(double));
        nap++;
      }
      else if ( MG_EDG(ptt.tag[j]) && !(ptt.tag[j] & MG_NOM) ) {
        ppt = &mesh->point[ip];
        assert(ppt->xp);
        pxp = &mesh->xpoint[ppt->xp];

        dd = no[0]*pxp->n1[0]+no[1]*pxp->n1[1]+no[2]*pxp->n1[2];
        if ( dd > 1.0-_MMG5_EPS ) continue;

        memcpy(pxp->n2,no,3*sizeof(double));
        /* a computation of the tangent with respect to these two normals is possible */
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
      }
    }
  }
  if ( !ns ) {
    _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
    return(ns);
  }

  /** 2. check if split by adjacent; besides, a triangle may have been splitted and not its adjacent
      (thus, the associated n2 may not exist) : update this normal if need be */
  nc = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* update face-edge flag */
    for (i=0; i<4; i++) {
      /* virtual triangle */
      memset(&ptt,0,sizeof(MMG5_Tria));
      if ( pt->xt && pxt->ftag[i] )
        _MMG5_tet2tri(mesh,k,i,&ptt);

      for (j=0; j<3; j++) {
        ia  = _MMG5_iarf[i][j];
        if ( MG_GET(pt->flag,ia) )                continue;
        if ( pt->xt && (pxt->tag[ia] & MG_REQ) )  continue;
        else if ( ptt.tag[j] & MG_REQ )           continue;
        ip1 = pt->v[_MMG5_iare[ia][0]];
        ip2 = pt->v[_MMG5_iare[ia][1]];
        ip  = _MMG5_hashGet(&hash,ip1,ip2);
        if ( ip > 0 ) {

          MG_SET(pt->flag,ia);
          nc++;
          /* ridge on a boundary face */
          if ( !(ptt.tag[j] & MG_GEO) && !(ptt.tag[j] & MG_NOM) )  continue;
          ppt = &mesh->point[ip];
          assert(ppt->xp);
          pxp = &mesh->xpoint[ppt->xp];
          if ( pt->xt )  ier = _MMG5_bezierCP(mesh,&ptt,&pb,MG_GET(pxt->ori,i));
          else  ier = _MMG5_bezierCP(mesh,&ptt,&pb,1);
          ier = _MMG3D_bezierInt(&pb,&uv[j][0],o,no,to);

          dd = no[0]*pxp->n1[0]+no[1]*pxp->n1[1]+no[2]*pxp->n1[2];
          if ( dd > 1.0-_MMG5_EPS ) continue;

          memcpy(pxp->n2,no,3*sizeof(double));
        }
      }
    }
  }
  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %d added\n",nc);
    fflush(stdout);
  }

  /** 3. Simulate splitting and delete points leading to invalid configurations */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  it = 1;
  nc = 0;
  do {
    ni = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) || (pt->tag & MG_REQ) || !pt->flag )  continue;
      memset(vx,0,6*sizeof(int));
      pt->flag = ic = 0;
      for (ia=0,i=0; i<3; i++) {
        for (j=i+1; j<4; j++,ia++) {
          if ( pt->xt && (mesh->xtetra[pt->xt].tag[ia] & MG_REQ) )  continue;
          vx[ia] = _MMG5_hashGet(&hash,pt->v[i],pt->v[j]);
          if ( vx[ia] > 0 ) {
            MG_SET(pt->flag,ia);
            if ( mesh->point[vx[ia]].flag > 2 )  ic = 1;
          }
        }
      }
      if ( !pt->flag )  continue;
      switch (pt->flag) {
      case 1: case 2: case 4: case 8: case 16: case 32:
        ier = _MMG3D_split1_sim(mesh,met,k,vx);
        break;
      case 11: case 21: case 38: case 56:
        ier = _MMG3D_split3_sim(mesh,met,k,vx);
        break;
      default:
        ier = _MMG5_split2sf_sim(mesh,met,k,vx);
        break;
      }
      if ( ier )  continue;

      ni++;
      if ( ic == 0 && _MMG3D_dichoto(mesh,met,k,vx) ) {
        for (ia=0; ia<6; ia++)
          if ( vx[ia] > 0 )  mesh->point[vx[ia]].flag++;
      }
      else {
        for (ia=0,i=0; i<3; i++) {
          for (j=i+1; j<4; j++,ia++) {
            if ( vx[ia] > 0 ) {
              p1 = &mesh->point[pt->v[_MMG5_iare[ia][0]]];
              p2 = &mesh->point[pt->v[_MMG5_iare[ia][1]]];
              ppt = &mesh->point[vx[ia]];
              ppt->c[0] = 0.5 * (p1->c[0] + p2->c[0]);
              ppt->c[1] = 0.5 * (p1->c[1] + p2->c[1]);
              ppt->c[2] = 0.5 * (p1->c[2] + p2->c[2]);
            }
          }
        }
      }
    }
    nc += ni;
  }
  while( ni > 0 && ++it < 20 );
  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %d corrected, %d invalid\n",nc,ni);
    fflush(stdout);
  }

  /** 4. splitting */
  ns = 0;
  ne = mesh->ne;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || !pt->flag || (pt->tag & MG_REQ) )  continue;
    memset(vx,0,6*sizeof(int));
    for (ia=0,i=0; i<3; i++) {
      for (j=i+1; j<4; j++,ia++) {
        if ( MG_GET(pt->flag,ia) )  {
          vx[ia] = _MMG5_hashGet(&hash,pt->v[i],pt->v[j]);
          assert(vx[ia]);
        }
      }
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32:  /* 1 edge split */
      _MMG5_split1(mesh,met,k,vx,typchk-1);
      ns++;
      break;
    case 11: case 21: case 38: case 56: /* 1 face (3 edges) subdivided */
      _MMG5_split3(mesh,met,k,vx,typchk-1);
      ns++;
      break;
    default:
      _MMG5_split2sf(mesh,met,k,vx,typchk-1);
      ns++;
      break;
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"       %7d elements splitted\n",nap);

  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(nap);
}



/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking permformed.
 * \return -1 if failed, number of new points otherwise.
 *
 * Split tetra into 4 when more than 1 boundary face.
 *
 */
static int _MMG5_anatet4(MMG5_pMesh mesh, MMG5_pSol met, char typchk) {
  MMG5_pTetra      pt;
  MMG5_pPoint      ppt;
  MMG5_pxTetra     pxt;
  int         k,ns,ier;
  char        nf,j;

  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    nf = 0;
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      for (j=0; j<4; j++)
        if ( pxt->ftag[j] & MG_BDY )  nf++;
    }
    if ( nf > 1 ) {
      ier  = _MMG5_split4bar(mesh,met,k,typchk-1);
      if ( !ier ) return(-1);
      ns++;
    }
    else {
      nf = 0;
      for (j=0; j<4; j++) {
        ppt = &mesh->point[pt->v[j]];
        if ( ppt->tag & MG_BDY )  nf++;
      }
      if ( nf == 4 ) {
        ier  = _MMG5_split4bar(mesh,met,k,typchk-1);
        if ( !ier ) return(-1);
        ns++;
      }
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d boundary elements splitted\n",ns);
  return(ns);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking for edges length.
 * \param patternMode flag to say if we perform vertex insertion by patterns
 * or by delaunay kernel.
 * \return 0 if fail, number of new points otherwise.
 *
 * Analyze tetrahedra and split if needed.
 *
 */
int _MMG5_anatet(MMG5_pMesh mesh,MMG5_pSol met,char typchk, int patternMode) {
  int     ier,nc,ns,nf,nnc,nns,nnf,it,maxit;

  /* analyze tetras : initial splitting */
  nns = nnc = nnf = it = 0;
  maxit = 5;
  mesh->gap = 0.5;
  do {
    /* memory free */
    _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

    if ( !mesh->info.noinsert ) {

      /* split tetra with more than 2 bdry faces */
      ier = _MMG5_anatet4(mesh,met,typchk);
      if ( ier < 0 )  return(0);
      ns = ier;

      /* analyze surface tetras */
      ier = _MMG5_anatets(mesh,met,typchk);

      if ( ier < 0 ) {
        fprintf(stderr,"  ## Unable to complete surface mesh. Exit program.\n");
        return(0);
      }
      ns += ier;
      if ( patternMode ) {
        /* analyze internal tetras */
        ier = _MMG5_anatetv(mesh,met,typchk);
        if ( ier < 0 ) {
          fprintf(stderr,"  ## Unable to complete volume mesh. Exit program.\n");
          return(0);
        }
        ns += ier;
      }
    }
    else  ns = 0;

    if ( !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"  ## Hashing problem. Exit program.\n");
      return(0);
    }
    if ( typchk == 2 && it == maxit-1 )  mesh->info.fem = 1;

    /* collapse short edges */
    if ( !mesh->info.noinsert ) {
      nc = _MMG5_coltet(mesh,met,typchk);
      if ( nc < 0 ) {
        fprintf(stderr,"  ## Unable to collapse mesh. Exiting.\n");
        return(0);
      }
    }
    else  nc = 0;

    /* attempt to swap */
    if ( !mesh->info.noswap ) {
      nf = _MMG5_swpmsh(mesh,met,NULL,typchk);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = _MMG5_swptet(mesh,met,1.1,NULL,typchk);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nf > 0 ){
#ifndef PATTERN
      fprintf(stdout,"                   ");
#endif
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nf);
    }
    if ( it > 3 && abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nf > 0 );

  if ( mesh->info.imprim ) {
    if ( (abs(mesh->info.imprim) < 5 || mesh->info.ddebug ) && nns+nnc > 0 ) {
#ifndef PATTERN
      fprintf(stdout,"                   ");
#endif
      fprintf(stdout, "     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnf,it);
    }
  }

  return(1);
}
