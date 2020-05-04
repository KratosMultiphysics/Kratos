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

extern char  ddb;

/**
 * \param mesh pointer toward the mesh structure.
 * \param k tetrahedron index.
 * \param ie face index of tetrahedron.
 * \param ptt pointer toward the output triangle.
 *
 * Set triangle corresponding to face ie of tetra k.
 *
 */
void MMG5_tet2tri(MMG5_pMesh mesh,int k,char ie,MMG5_Tria *ptt) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  char    i;

  pt = &mesh->tetra[k];
  memset(ptt,0,sizeof(MMG5_Tria));
  ptt->v[0] = pt->v[MMG5_idir[ie][0]];
  ptt->v[1] = pt->v[MMG5_idir[ie][1]];
  ptt->v[2] = pt->v[MMG5_idir[ie][2]];
  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    ptt->ref = pxt->ref[ie];
    for (i=0; i<3; i++) {
      ptt->edg[i] = pxt->edg[MMG5_iarf[ie][i]];
      ptt->tag[i] = pxt->tag[MMG5_iarf[ie][i]];
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
 * \return 1 if success, 0 if fail.
 *
 * Find acceptable position for splitting.
 *
 */
int MMG3D_dichoto(MMG5_pMesh mesh,MMG5_pSol met,int k,int *vx) {
  MMG5_pTetra  pt;
  MMG5_pPoint  pa,pb,ps;
  double       o[6][3],p[6][3];
  float        to,tp,t;
  int          ia,ib,ier,it,maxit;
  char         i;

  ier = 1;
  pt = &mesh->tetra[k];

  /* get point on surface and along segment for edge split */
  for (i=0; i<6; i++) {
    memset(p[i],0,3*sizeof(double));
    memset(o[i],0,3*sizeof(double));
    if ( vx[i] > 0 ) {
      ia = pt->v[MMG5_iare[i][0]];
      ib = pt->v[MMG5_iare[i][1]];
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
      ier = MMG3D_split1_sim(mesh,met,k,vx);
      break;
    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10:
      ier =MMG3D_split2sf_sim(mesh,met,k,vx);
      break;
    case 33: case 18: case 12:
      ier = MMG3D_split2_sim(mesh,met,k,vx);
      break;
    case 11: case 21: case 38: case 56:
      ier = MMG3D_split3_sim(mesh,met,k,vx);
      break;
    case 7: case 25: case 42: case 52:
      ier =MMG3D_split3cone_sim(mesh,met,k,vx);
      break;
    case 35: case 19: case 13: case 37: case 22: case 28: case 26:
    case 14: case 49: case 50: case 44: case 41:
      ier = MMG3D_split3op_sim(mesh,met,k,vx);
      break;
    case 23: case 29: case 53: case 60: case 57: case 58:
    case 27: case 15: case 43: case 39: case 54: case 46:
      ier = MMG3D_split4sf_sim(mesh,met,k,vx);
      break;
    case 30: case 45: case 51:
      ier = MMG3D_split4op_sim(mesh,met,k,vx);
      break;
    case 62: case 61: case 59: case 55: case 47: case 31:
      ier = MMG3D_split5_sim(mesh,met,k,vx);
      break;
    case 63:
      ier = MMG3D_split6_sim(mesh,met,k,vx);
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

  /* For very ill-shaped elements, we can have no valid position */
  switch (pt->flag) {
  case 1: case 2: case 4: case 8: case 16: case 32:
    ier = MMG3D_split1_sim(mesh,met,k,vx);
    break;
  case 48: case 24: case 40: case 6: case 34: case 36:
  case 20: case 5: case 17: case 9: case 3: case 10:
    ier =MMG3D_split2sf_sim(mesh,met,k,vx);
    break;
  case 33: case 18: case 12:
    ier = MMG3D_split2_sim(mesh,met,k,vx);
    break;
  case 11: case 21: case 38: case 56:
    ier = MMG3D_split3_sim(mesh,met,k,vx);
    break;
  case 7: case 25: case 42: case 52:
    ier =MMG3D_split3cone_sim(mesh,met,k,vx);
    break;
  case 35: case 19: case 13: case 37: case 22: case 28: case 26:
  case 14: case 49: case 50: case 44: case 41:
    ier = MMG3D_split3op_sim(mesh,met,k,vx);
    break;
  case 23: case 29: case 53: case 60: case 57: case 58:
  case 27: case 15: case 43: case 39: case 54: case 46:
    ier = MMG3D_split4sf_sim(mesh,met,k,vx);
    break;
  case 30: case 45: case 51:
    ier = MMG3D_split4op_sim(mesh,met,k,vx);
    break;
  case 62: case 61: case 59: case 55: case 47: case 31:
    ier = MMG3D_split5_sim(mesh,met,k,vx);
    break;
  case 63:
    ier = MMG3D_split6_sim(mesh,met,k,vx);
    break;
  }
  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param list pointer toward the shell of edge.
 * \param ret double of the number of tetrahedra in the shell.
 * \param ip new point index.
 *
 * \return 1 if success
 * \return 0 if fail due to a very bad quality
 * \return 2 if fail due to a sharp angle creation
 *
 * Find acceptable position for MMG5_split1b, passing the shell of
 * considered edge, starting from o point.
 *
 */
int MMG3D_dichoto1b(MMG5_pMesh mesh,MMG5_pSol met,int *list,int ret,int ip) {
  MMG5_pTetra  pt;
  MMG5_pPoint  p0,p1,ppt;
  int          iel,np,nq,it,maxit;
  double       m[3],o[3],tp,to,t;
  char         ia,ier;

  iel = list[0] / 6;
  ia  = list[0] % 6;
  pt  = &mesh->tetra[iel];

  np = pt->v[MMG5_iare[ia][0]];
  nq = pt->v[MMG5_iare[ia][1]];
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
  tp    = 1.0;
  to    = 0.0;
  do {
    t = 0.5*(to + tp);
    ppt->c[0] = m[0] + t*(o[0]-m[0]);
    ppt->c[1] = m[1] + t*(o[1]-m[1]);
    ppt->c[2] = m[2] + t*(o[2]-m[2]);

    ier = MMG3D_simbulgept(mesh,met,list,ret,ip);
    assert ( (!mesh->info.ddebug) || (mesh->info.ddebug && ier != -1) );
    if ( ier == 1 )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );
  if ( !ier )  t = to;

  /* For very ill-shaped elements, we can have novalid position */
  ppt->c[0] = m[0] + t*(o[0]-m[0]);
  ppt->c[1] = m[1] + t*(o[1]-m[1]);
  ppt->c[2] = m[2] + t*(o[2]-m[2]);

  return  MMG3D_simbulgept(mesh,met,list,ret,ip) ;
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
char MMG5_chkedg(MMG5_pMesh mesh,MMG5_Tria *pt,char ori, double hmax,
                  double hausd, int locPar) {
  MMG5_pPoint   p[3];
  MMG5_xPoint  *pxp;
//  MMG5_pPar     par;
  double   n[3][3],t[3][3],nt[3],*n1,*n2,t1[3],t2[3];
  double   ps,ps2,ux,uy,uz,ll,il,alpha,dis,hma2;
  int      ia,ib,ic;//l,info;
  char     i,i1,i2;
  static char mmgWarn0 = 0, mmgWarn1 = 0;

  ia   = pt->v[0];
  ib   = pt->v[1];
  ic   = pt->v[2];
  p[0] = &mesh->point[ia];
  p[1] = &mesh->point[ib];
  p[2] = &mesh->point[ic];
  pt->flag = 0;

  n1 = n2 = NULL;

  /* normal recovery */
  for (i=0; i<3; i++) {
    if ( MG_SIN(p[i]->tag) ) {
      MMG5_nortri(mesh,pt,n[i]);
      if(!ori) {
        n[i][0] *= -1.0;
        n[i][1] *= -1.0;
        n[i][2] *= -1.0;
      }
    }
    else if ( (p[i]->tag & MG_NOM) || (p[i]->tag & MG_OPNBDY) ){
      MMG5_nortri(mesh,pt,n[i]);
      if(!ori) {
        n[i][0] *= -1.0;
        n[i][1] *= -1.0;
        n[i][2] *= -1.0;
      }
      assert(p[i]->xp);
      memcpy(&t[i],p[i]->n,3*sizeof(double));
    }
    else {
      assert(p[i]->xp);
      pxp = &mesh->xpoint[p[i]->xp];
      if ( MG_EDG(p[i]->tag) ) {
        memcpy(&t[i],p[i]->n,3*sizeof(double));
        MMG5_nortri(mesh,pt,nt);
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
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

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

    hma2 = MMG3D_LLONG*MMG3D_LLONG*hmax*hmax;

    /* check length */
    ux = p[i2]->c[0] - p[i1]->c[0];
    uy = p[i2]->c[1] - p[i1]->c[1];
    uz = p[i2]->c[2] - p[i1]->c[2];
    ll = ux*ux + uy*uy + uz*uz;
    if ( ll < MMG5_EPSD )  continue;
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
          if ( !mmgWarn0 ) {
            fprintf(stderr,"\n  ## Warning: %s: a- at least 1 geometrical"
                    " problem\n",__func__);
            mmgWarn0 = 1;
          }
          return 0;
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
          if ( !mmgWarn1 ) {
            fprintf(stderr,"\n  ## Warning: %s: b- at least 1 geometrical"
                    " problem\n",__func__);
            mmgWarn1 = 1;
          }
          return 0;
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
      if ( !MMG5_BezierTgt(p[i1]->c,p[i2]->c,n1,n2,t1,t2) ) {
        t1[0] = ux * il;
        t1[1] = uy * il;
        t1[2] = uz * il;

        t2[0] = -ux * il;
        t2[1] = -uy * il;
        t2[2] = -uz * il;
      }
    }
    alpha = MMG5_BezierGeod(p[i1]->c,p[i2]->c,t1,t2);
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
  return pt->flag;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure (only for delaunay).
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 * \return -1 if failed and swap number otherwise.
 *
 * Search for boundary edges that could be swapped for geometric
 * approximation.
 *
 */
int MMG5_swpmsh(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree PROctree, int typchk) {
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
          ia  = MMG5_iarf[i][j];

          /* No swap of geometric edge */
          if ( MG_EDG(pxt->tag[ia]) || (pxt->tag[ia] & MG_REQ) ||
               (pxt->tag[ia] & MG_NOM) )
            continue;

          ret = MMG5_coquilface(mesh,k,i,ia,list,&it1,&it2,0);
          ilist = ret / 2;
          if ( ret < 0 )  return -1;
          /* CAUTION: trigger collapse with 2 elements */
          if ( ilist <= 1 )  continue;
          ier = MMG5_chkswpbdy(mesh,met,list,ilist,it1,it2,typchk);
          if ( ier <  0 )
            return -1;
          else if ( ier ) {
            ier = MMG5_swpbdy(mesh,met,list,ret,it1,PROctree,typchk);
            if ( ier > 0 )  ns++;
            else if ( ier < 0 )  return -1;
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

  return nns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param crit coefficient of quality improvment.
 * \param PROctree pointer toward the PROctree structure in delaunay mode and
 * toward the \a NULL pointer otherwise
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion)
 * \param testmark all the tets with a mark less than testmark will not be treated.
 *
 * \return -1 if fail, the number of swap otherwise.
 *
 * Internal edge flipping.
 *
 */
int MMG5_swptet(MMG5_pMesh mesh,MMG5_pSol met,double crit,double declic,
                 MMG3D_pPROctree PROctree,int typchk,int testmark) {
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
      else if ( pt->mark < testmark )  continue;

      if ( pt->qual > declic )  continue;

      for (i=0; i<6; i++) {
        /* Prevent swap of a ref or tagged edge */
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
          if ( pxt->edg[i] || pxt->tag[i] ) continue;
        }

        nconf = MMG5_chkswpgen(mesh,met,k,i,&ilist,list,crit,typchk);
        if ( nconf<0 ) return -1;
        else if ( nconf ) {
          ier = MMG5_swpgen(mesh,met,nconf,ilist,list,PROctree,typchk);
          if ( ier > 0 )  ns++;
          else if ( ier < 0 ) return -1;
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ++it < maxit && ns > 0 );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8d edge swapped\n",nns);

  return nns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param clickSurf triangle quality threshold under which we want to move
 * \param clickVol  tetra    quality threshold under which we want to move
 * \param moveVol internal move
 * \param improveSurf forbid surface degradation during the move
 * \param improveVolSurf forbid volume degradation during the surfacic move
 * \param improveVol forbid volume degradation during the move
 * \param maxit maximum number of iteration
 * \param testmark all the tets with a mark less than testmark will not be treated.
 * \return -1 if failed, number of moved points otherwise.
 *
 * Analyze tetrahedra and move points so as to make mesh more uniform.
 *
 */
int MMG5_movtet(MMG5_pMesh mesh,MMG5_pSol met, MMG3D_pPROctree PROctree,
                 double clickSurf,double clickVol,int moveVol, int improveSurf,
                 int improveVolSurf, int improveVol, int maxit,int testmark) {
  MMG5_pTetra        pt;
  MMG5_pPoint        ppt;
  MMG5_pxTetra       pxt;
  MMG5_Tria          tt;
  double        *n,caltri;
  int           i,k,ier,nm,nnm,ns,lists[MMG3D_LMAX+2],listv[MMG3D_LMAX+2],ilists,ilistv,it;
  unsigned char j,i0,base;

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
      else if ( pt->mark < testmark )  continue;

      /* point j on face i */
      for (i=0; i<4; i++) {
        for (j=0; j<3; j++) {
          if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
          }
          else  pxt = 0;

          i0  = MMG5_idir[i][j];
          ppt = &mesh->point[pt->v[i0]];
          if ( ppt->flag == base )  continue;
          else if ( MG_SIN(ppt->tag) )  continue;

          if( pt->xt && (pxt->ftag[i] & MG_BDY)) {
            MMG5_tet2tri(mesh,k,i,&tt);
            caltri = MMG5_caltri(mesh,met,&tt);

            if ( caltri >= clickSurf) {
              j = 3;
              continue;
            }
          }
          if ( maxit != 1 ) {
            ppt->flag = base;
          }
          ier = 0;
          if ( ppt->tag & MG_BDY ) {
            /* Catch a boundary point by a boundary face */
            if ( (!pt->xt) || !(MG_BDY & pxt->ftag[i]) )  continue;
            else if( ppt->tag & MG_NOM ){
              if ( mesh->adja[4*(k-1)+1+i] ) continue;
              ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,1);
              if( !ier )  continue;
              else if ( ier>0 )
                ier = MMG5_movbdynompt(mesh,met,PROctree,listv,ilistv,lists,ilists,improveVolSurf);
              else
                return -1;
            }
            else if ( ppt->tag & MG_GEO ) {
              ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
              if ( !ier )  continue;
              else if ( ier>0 )
                ier = MMG5_movbdyridpt(mesh,met,PROctree,listv,ilistv,lists,ilists,improveVolSurf);
              else
                return -1;
            }
            else if ( ppt->tag & MG_REF ) {
              ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
              if ( !ier )
                continue;
              else if ( ier>0 )
                ier = MMG5_movbdyrefpt(mesh,met,PROctree,listv,ilistv,lists,ilists,improveVolSurf);
              else
                return -1;
            }
            else {
              ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
              if ( !ier )
                continue;
              else if ( ier<0 )
                return -1;

              n = &(mesh->xpoint[ppt->xp].n1[0]);

              /* If the orientation of the tetra face is
               * compatible with the triangle (MG_GET(ori,i)), we know that we
               * are well orientated. */
              if ( !MG_GET(pxt->ori,i) ) {
                if ( !MMG5_directsurfball(mesh,pt->v[i0],lists,ilists,n) )
                  continue;
              }
              ier = MMG5_movbdyregpt(mesh,met,PROctree,listv,ilistv,
                                      lists,ilists,improveSurf,improveVolSurf);
              if (ier < 0 ) return -1;
              else if ( ier )  ns++;
            }
          }
          else if ( moveVol && (pt->qual < clickVol) ) {
            ilistv = MMG5_boulevolp(mesh,k,i0,listv);
            if ( !ilistv )  continue;
            ier = MMG5_movintpt(mesh,met,PROctree,listv,ilistv,improveVol);
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

  return nnm;
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
static int MMG5_coltet(MMG5_pMesh mesh,MMG5_pSol met,char typchk) {
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
        if ( pt->xt && (pxt->tag[MMG5_iarf[i][j]] & MG_REQ) )  continue;
        ip = MMG5_idir[i][MMG5_inxt2[j]];
        iq = MMG5_idir[i][MMG5_iprv2[j]];

        p0 = &mesh->point[pt->v[ip]];
        p1 = &mesh->point[pt->v[iq]];
        if ( p0->flag == base )  continue;
        else if ( (p0->tag & MG_REQ) || (p0->tag > p1->tag) )  continue;


        /* Ball of point: computed here if needed for the local parameter
         * evaluation, after length check otherwise (because the ball
         * computation is time consuming) */
        if ( mesh->info.npar ) {
          if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
            tag = pxt->tag[MMG5_iarf[i][j]];
            isnm = (tag & MG_NOM);

            if ( p0->tag > tag ) continue;
            if ( isnm && mesh->adja[4*(k-1)+1+i] )  continue;
            if (MMG5_boulesurfvolp(mesh,k,ip,i,
                                    list,&ilist,lists,&ilists,p0->tag & MG_NOM) < 0 )
              return -1;
          }
          else {
            ilist = MMG5_boulevolp(mesh,k,ip,list);
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
          if ( ll > hmi2*MMG3D_LSHRT )  continue;
        }
        else if ( typchk == 2 ) {
          ll = MMG5_lenedg(mesh,met,MMG5_iarf[i][j],pt);
          // Case of an internal tetra with 4 ridges vertices.
          if ( ll == 0 ) continue;
          if ( ll > MMG3D_LSHRT )  continue;
        }


        /* Ball computation if not computed before */
        if ( !mesh->info.npar ) {
          if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
            tag = pxt->tag[MMG5_iarf[i][j]];
            isnm = (tag & MG_NOM);

            if ( p0->tag > tag ) continue;
            if ( isnm && mesh->adja[4*(k-1)+1+i] )  continue;
            if (MMG5_boulesurfvolp(mesh,k,ip,i,
                                    list,&ilist,lists,&ilists,p0->tag & MG_NOM) < 0 )
              return -1;
          }
          else {
            ilist = MMG5_boulevolp(mesh,k,ip,list);
          }
        }

        /* boundary face: collapse ip on iq */
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
          tag = pxt->tag[MMG5_iarf[i][j]];
          tag |= MG_BDY;

          isnm = ( tag & MG_NOM );
          if ( isnm ) {
            if ( mesh->adja[4*(k-1)+1+i] )  continue;
          }
          if ( p0->tag > tag )  continue;
          ilist = MMG5_chkcol_bdy(mesh,met,k,i,j,list,ilist,lists,ilists,typchk);
        }
        /* internal face */
        else {
          isnm = 0;
          if ( p0->tag & MG_BDY )  continue;
          ilist = MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,typchk);
        }

        if ( ilist > 0 ) {
          ier = MMG5_colver(mesh,met,list,ilist,iq,typchk);
          if ( ier < 0 ) return -1;
          else if ( ier ) {
            MMG3D_delPt(mesh,ier);
            break;
          }
        }
        else if (ilist < 0 ) return -1;
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

  return nc;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \return 0 if failed, 1 if success
 *
 * Delete the points inserted by pattern if the pattern step fail.
 *
 */
static inline
int MMG3D_delPatternPts(MMG5_pMesh mesh,MMG5_Hash hash)
{
  MMG5_pTetra   pt;
  int           vx[6],k,ia,i,j;

  /* Delete the useless added points */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) ) continue;

    memset(vx,0,6*sizeof(int));
    for (ia=0,i=0; i<3; i++) {
      for (j=i+1; j<4; j++,ia++) {
        if ( pt->xt && (mesh->xtetra[pt->xt].tag[ia] & MG_REQ) ) continue;
        vx[ia] = MMG5_hashGet(&hash,pt->v[i],pt->v[j]);
        if ( vx[ia] > 0 ) {
          MMG3D_delPt(mesh,vx[ia]);
          if ( !MMG5_hashUpdate(&hash,pt->v[i],pt->v[j],-1) ) {
            fprintf(stderr,"\n  ## Error: %s: unable to delete point idx"
                    " along edge %d %d.\n", __func__,
                    MMG3D_indPt(mesh,pt->v[i]),
                    MMG3D_indPt(mesh,pt->v[j]));
            MMG5_DEL_MEM(mesh,hash.item);
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking permformed for edge length (hmax or MMG3D_LLONG criterion).
 * \return -1 if failed.
 * \return number of new points.
 *
 * Analyze volume tetra and split if needed.
 *
 */
static int
MMG5_anatetv(MMG5_pMesh mesh,MMG5_pSol met,char typchk) {
  MMG5_pTetra   pt;
  MMG5_pPoint   p1,p2;
  MMG5_xTetra  *pxt;
  MMG5_Hash    hash;
  MMG5_pPar     par;
  double   ll,o[3],ux,uy,uz,hma2,mincal;
  int      l,vx[6],k,ip,ip1,ip2,nap,ns,ne,memlack,ier;
  char     i,j,ia;

  /** 1. analysis */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return -1;
  memlack = ns = nap = 0;

  /* Hash all boundary and required edges, and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (i=0; i<6; i++) {
        ip1 = pt->v[MMG5_iare[i][0]];
        ip2 = pt->v[MMG5_iare[i][1]];
        ip  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return -1;
      }
      continue;
    }

    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] & MG_BDY ) {
        for (j=0; j<3; j++) {
          ip1 = pt->v[MMG5_idir[i][MMG5_inxt2[j]]];
          ip2 = pt->v[MMG5_idir[i][MMG5_iprv2[j]]];
          ip  = -1;
          if ( !MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return -1;
        }
      }
    }
  }

  /** 2. Set flags and split internal edges */
  mincal = FLT_MAX;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = 0;
    for (i=0; i<6; i++) {
      ip  = -1;
      ip1 = pt->v[MMG5_iare[i][0]];
      ip2 = pt->v[MMG5_iare[i][1]];
      p1  = &mesh->point[ip1];
      p2  = &mesh->point[ip2];
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( pxt->tag[i] & MG_REQ ) continue;
      }
      else  pxt = 0;
      if ( (p1->tag & MG_BDY) && (p2->tag & MG_BDY) ) {
        ip = MMG5_hashGet(&hash,ip1,ip2);
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
          hma2 = MMG3D_LLONG*MMG3D_LLONG*hma2*hma2;
          if ( ll > hma2 ) {
            ip = MMG5_hashGet(&hash,ip1,ip2);
            mincal = MG_MIN(mincal,pt->qual);
          }
        }
        else if ( typchk == 2 ) {
          ll = MMG5_lenedg(mesh,met,i,pt);
          // Case of an internal tetra with 4 ridges vertices.
          if ( ll == 0 ) continue;
          if ( ll > MMG3D_LLONG ) {
            ip = MMG5_hashGet(&hash,ip1,ip2);
            mincal = MG_MIN(mincal,pt->qual);
          }
        }
      }
      if ( ip < 0 ) continue;
      else if ( !ip ) {
        /* new midpoint */
        o[0] = 0.5 * (p1->c[0]+p2->c[0]);
        o[1] = 0.5 * (p1->c[1]+p2->c[1]);
        o[2] = 0.5 * (p1->c[2]+p2->c[2]);

        ip  = MMG3D_newPt(mesh,o,0);
        if ( !ip ) {
          /* reallocation of point table */

          MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                               fprintf(stderr,"\n  ## Warning: %s: unable to"
                                       " allocate a new point\n",__func__);
                               MMG5_INCREASE_MEM_MESSAGE();
                               memlack=1;
                               goto split
                               ,o,0);
        }

        assert ( met );
        if ( met->m ) {
          if ( typchk == 1 && (met->size>1) )
            ier = MMG3D_intmet33_ani(mesh,met,k,i,ip,0.5);
          else
            ier = MMG5_intmet(mesh,met,k,i,ip,0.5);

          if (!ier) {
            // Unable to compute the metric
            return -1;
          }
          else if ( ier < 0 ) {
            MMG3D_delPt(mesh,ip);
            continue;
          }
        }

        if ( !MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return -1;
        MG_SET(pt->flag,i);
        nap++;
      }
    }
  }
  if ( !nap )  {
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /** 3. check and split */
split:
  if ( mincal < MMG5_EPS ) {
    /* Delete the useless added points */
    if ( !MMG3D_delPatternPts(mesh,hash) ) return -1;

    /* Avoid the creation of bad quality elements */
    if ( mesh->info.imprim > 5 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: too bad quality for the worst element."
              " Volumic patterns skipped.\n",__func__);
    }

    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

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
        vx[ia] = MMG5_hashGet(&hash,pt->v[i],pt->v[j]);
        if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
      }
    }

    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      if ( ! MMG5_split1(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;
    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      if ( ! MMG5_split2sf(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 33: case 18: case 12: /* 2 opposite edges split */
      if ( ! MMG5_split2(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 11: case 21: case 38: case 56: /* 3 edges on the same faces splitted */
      if ( ! MMG5_split3(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      if ( ! MMG5_split3cone(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 35: case 19: case 13: case 37: case 22: case 28: case 26:
    case 14: case 49: case 50: case 44: case 41: /* 3 edges on opposite configuration splitted */
      if ( ! MMG5_split3op(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 23: case 29: case 53: case 60: case 57: case 58:
    case 27: case 15: case 43: case 39: case 54: case 46: /* 4 edges with 3 lying on the same face splitted */
      if ( ! MMG5_split4sf(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

      /* 4 edges with no 3 lying on the same face splitted */
    case 30: case 45: case 51:
      if ( ! MMG5_split4op(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 62: case 61: case 59: case 55: case 47: case 31: /* 5 edges split */
      if ( ! MMG5_split5(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 63: /* 6 edges split */
      if ( ! MMG5_split6(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;
    }
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",nap);

  MMG5_DEL_MEM(mesh,hash.item);
  if ( memlack )  return -1;
  return nap;
}

/**
 * \param ppt pointer toward the point that we update
 * \param pxp point toward the \a oot xpoint
 * \param no normal at ppt
 * \return 0 if failed, 1 if success.
 *
 * Starting from a point for which the normal \a pxp->n1 is already stored,
 * store \a no in \a pxp->n2, compute the tangent with respect to this two
 * normals and store it in \a ppt->n.
 *
 */
static inline int
MMG3D_storeGeom(MMG5_pPoint ppt, MMG5_pxPoint pxp, double no[3]) {
  double dd;

  dd = no[0]*pxp->n1[0]+no[1]*pxp->n1[1]+no[2]*pxp->n1[2];
  if ( dd > 1.0-MMG5_EPS ) return 0;

  memcpy(pxp->n2,no,3*sizeof(double));

  /* a computation of the tangent with respect to these two normals is possible */
  ppt->n[0] = pxp->n1[1]*pxp->n2[2] - pxp->n1[2]*pxp->n2[1];
  ppt->n[1] = pxp->n1[2]*pxp->n2[0] - pxp->n1[0]*pxp->n2[2];
  ppt->n[2] = pxp->n1[0]*pxp->n2[1] - pxp->n1[1]*pxp->n2[0];
  dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
  if ( dd > MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    ppt->n[0] *= dd;
    ppt->n[1] *= dd;
    ppt->n[2] *= dd;
  }
  assert ( dd>MMG5_EPSD2 );

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the tetra to split.
 * \param pt tetra to split
 * \param pxt associated xtetra
 * \param imax index of the edge to split to split
 * \param typchk type of check
 * \param chkRidTet check for ridge metric
 * \param *warn \a warn is set to 1 if we don't have enough memory to complete mesh.
 *
 * \return -1 if fail, 0 if we can't split but the upper loop may continue, 1 if
 * the edge is splitted, 2 if we can't split due to lack of memory
 *
 * Split a surface edge using split1b
 *
 */
int MMG3D_splsurfedge( MMG5_pMesh mesh,MMG5_pSol met,int k,
                      MMG5_pTetra pt,MMG5_pxTetra pxt,char imax,char typchk,
                      char chkRidTet,int *warn ) {
  MMG5_Tria    ptt;
  MMG5_pPoint  p0,p1,ppt;
  MMG5_pxPoint pxp;
  double       dd,o[3],to[3],no1[3],no2[3],v[3];
  int          ip,ip1,ip2,list[MMG3D_LMAX+2],ilist;
  int          ref,ier;
  int16_t      tag;
  char         j,i,i1,i2,ifa0,ifa1;

  assert ( pxt = &mesh->xtetra[pt->xt] );

  /* proceed edges according to lengths */
  ifa0 = MMG5_ifar[imax][0];
  ifa1 = MMG5_ifar[imax][1];
  i  = (pxt->ftag[ifa1] & MG_BDY) ? ifa1 : ifa0;
  j  = MMG5_iarfinv[i][imax];
  i1 = MMG5_idir[i][MMG5_inxt2[j]];
  i2 = MMG5_idir[i][MMG5_iprv2[j]];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];
  p0  = &mesh->point[ip1];
  p1  = &mesh->point[ip2];

  ref = pxt->edg[imax];
  tag = pxt->tag[imax];

  if ( tag & MG_REQ ) { return 0; }

  tag |= MG_BDY;

  ilist = MMG5_coquil(mesh,k,imax,list);
  if ( !ilist )  return 0;
  else if ( ilist < 0 ) { return -1; }

  if ( tag & MG_NOM ){
    if( !MMG5_BezierNom(mesh,ip1,ip2,0.5,o,no1,to) ) { return 0; }

    else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
      MMG5_tet2tri(mesh,k,i,&ptt);
      MMG5_nortri(mesh,&ptt,no1);
      if ( !MG_GET(pxt->ori,i) ) {
        no1[0] *= -1.0;
        no1[1] *= -1.0;
        no1[2] *= -1.0;
      }
    }
  }
  else if ( tag & MG_GEO ) {
    if ( !MMG5_BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) ) { return 0; }

    if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
      MMG5_tet2tri(mesh,k,i,&ptt);
      MMG5_nortri(mesh,&ptt,no1);
      no2[0] = to[1]*no1[2] - to[2]*no1[1];
      no2[1] = to[2]*no1[0] - to[0]*no1[2];
      no2[2] = to[0]*no1[1] - to[1]*no1[0];
      dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        no2[0] *= dd;
        no2[1] *= dd;
        no2[2] *= dd;
      }
    }
  }
  else if ( tag & MG_REF ) {
    if ( !MMG5_BezierRef(mesh,ip1,ip2,0.5,o,no1,to) ) { return 0; }
    else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
      MMG5_tet2tri(mesh,k,i,&ptt);
      MMG5_nortri(mesh,&ptt,no1);
      if ( !MG_GET(pxt->ori,i) ) {
        no1[0] *= -1.0;
        no1[1] *= -1.0;
        no1[2] *= -1.0;
      }
    }
  }
  else {
    if ( !MMG5_norface(mesh,k,i,v) ) { return 0; }
    if ( !MMG5_BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) { return 0; }
    else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
      MMG5_tet2tri(mesh,k,i,&ptt);
      MMG5_nortri(mesh,&ptt,no1);
      if ( !MG_GET(pxt->ori,i) ) {
        no1[0] *= -1.0;
        no1[1] *= -1.0;
        no1[2] *= -1.0;
      }
    }
  }

  ip = MMG3D_newPt(mesh,o,tag);
  if ( !ip ) {
    /* reallocation of point table */
    MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                        *warn=1;
                        return 2;
                        ,o,tag);
  }

  if ( met->m ) {
    if ( typchk == 1 && (met->size>1) ) {
      ier = MMG3D_intmet33_ani(mesh,met,k,imax,ip,0.5);
    }
    else {
      ier = MMG5_intmet(mesh,met,k,imax,ip,0.5);
    }
    if ( !ier ) {
      MMG3D_delPt(mesh,ip);
      return -1;
    }
    else if (ier < 0) {
      MMG3D_delPt(mesh,ip);
      return 0;
    }
  }
  ier = MMG3D_simbulgept(mesh,met,list,ilist,ip);
  assert ( (!mesh->info.ddebug) || (mesh->info.ddebug && ier != -1) );
  if ( ier < 0 || ier == 2 ) {
    MMG3D_delPt(mesh,ip);
    return 0;
  }
  else if ( ier == 0 ) {
    ier = MMG3D_dichoto1b(mesh,met,list,ilist,ip);
  }
  if ( ier == 1 ) { ier = MMG5_split1b(mesh,met,list,ilist,ip,1,typchk-1,chkRidTet); }

  /* if we realloc memory in MMG5_split1b pt and pxt pointers are not valid */
  pt = &mesh->tetra[k];
  pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to split.\n",__func__);
    return -1;
  }
  else if ( ier == 0 || ier == 2 ) {
    MMG3D_delPt(mesh,ip);
    return 0;
  }

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
  else {
    memcpy(pxp->n1,no1,3*sizeof(double));
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of tetra thath we check
 * \param pt pointer toward the tetra that we check
 * \param pxt pointer toward the xtetra that we check
 * \param i index of the face in \a k that we check
 * \param ptt pointer toward the virtual triangle build from the face \i of \a k.
 * \param typchk type of checking permformed for edge length (hmax or MMG3D_LLONG criterion).
 *
 * \return 1 if success, 0 if fail.
 *
 * Mark edges to split on geometric criterion (mark stored in \a pt->flag)
 *
 */
static
int MMG3D_chkbdyface(MMG5_pMesh mesh,MMG5_pSol met,int k,MMG5_pTetra pt,
                     MMG5_pxTetra pxt,char i,MMG5_pTria ptt,char typchk ) {

  MMG5_pPar    par;
  double       len,hmax,hausd;
  int          l,ip1,ip2;
  int8_t       isloc;
  char         j,i1,i2,ia;

  if ( typchk == 1 ) {

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
          if ( par->ref != ptt->ref ) continue;

          hmax = MG_MIN(hmax,par->hmax);
          hausd = MG_MIN(hausd,par->hausd);
          break;
        }
      }
      else {
        for ( l=0; l<mesh->info.npar; ++l ) {
          par = &mesh->info.par[l];

          if ( par->elt != MMG5_Triangle )  continue;
          if ( par->ref != ptt->ref ) continue;

          hmax  = par->hmax;
          hausd = par->hausd;
          isloc = 1;
          break;
        }
      }
    }

    if ( !MMG5_chkedg(mesh,ptt,MG_GET(pxt->ori,i),hmax,hausd,isloc) ) {
      return 0;
    }

    /* put back flag on tetra */
    for (j=0; j<3; j++){
      if ( pxt->tag[MMG5_iarf[i][j]] & MG_REQ )  continue;
      if ( MG_GET(ptt->flag,j) )  MG_SET(pt->flag,MMG5_iarf[i][j]);
    }
  }
  else if ( typchk == 2 ) {
    for (j=0; j<3; j++) {
      ia = MMG5_iarf[i][j];
      if ( pxt->tag[ia] & MG_REQ )  continue;
      i1  = MMG5_iare[ia][0];
      i2  = MMG5_iare[ia][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      len = MMG5_lenedg(mesh,met,ia,pt);

      assert( isfinite(len) && (len!=-len) );

      // Case of an internal tetra with 4 ridges vertices.
      if ( len == 0 )  continue;
      if ( len > MMG3D_LLONG )  MG_SET(pt->flag,ia);
      /* Treat here the ridges coming from a corner (we can not do that after
       * because the corner don't have xpoints) */
      if ( (mesh->point[ip1].tag & MG_CRN) ||  (mesh->point[ip2].tag & MG_CRN) ) {
        if ( len > MMG3D_LOPTL )  MG_SET(pt->flag,ia);
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking permformed for edge length (hmax or MMG3D_LLONG criterion).
 * \return -1 if failed.
 * \return number of new points.
 *
 * Split surface edges on geometric criterion.
 *
 */
static int MMG3D_anatets_ani(MMG5_pMesh mesh,MMG5_pSol met,char typchk) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_Tria    ptt;
  double       len,lmax;
  double       ux,uy,uz;
  int          k,ip1,ip2;
  int          ns,ier,warn;
  char         imax,j,i,i1,i2;

  assert ( met->m && met->size==6 );

  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return -1;
  }

  ns = 0;
  warn = 0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) || !pt->xt )  continue;

    /* check boundary face cut w/r Hausdorff or hmax */
    pt->flag = 0;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++){
      if ( pxt->ftag[i] & MG_REQ )     continue;
      if ( !(pxt->ftag[i] & MG_BDY) )  continue;

      /* virtual triangle */
      MMG5_tet2tri(mesh,k,i,&ptt);
      if ( !MG_GET(pxt->ori,i) ) continue;

      if ( !MMG3D_chkbdyface(mesh,met,k,pt,pxt,i,&ptt,typchk) ) { continue; }
    }

    /** Split only the longest edge */
    imax = 6;
    lmax = 0.;
    for ( j=0; j<6; ++j ) {
      if ( !MG_GET(pt->flag,j) ) { continue; }

      i1  = MMG5_iare[j][0];
      i2  = MMG5_iare[j][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];

      ux = mesh->point[ip1].c[0] - mesh->point[ip2].c[0];
      uy = mesh->point[ip1].c[1] - mesh->point[ip2].c[1];
      uz = mesh->point[ip1].c[2] - mesh->point[ip2].c[2];

      len = ux*ux + uy*uy + uz*uz;
      if ( len <= lmax) { continue; }
      lmax = len;
      imax = j;
    }

    pt->flag = 0;
    if ( imax < 6 ) { MG_SET(pt->flag,imax); }

    if ( !pt->flag )  continue;

    ier = MMG3D_splsurfedge( mesh,met,k,pt,pxt,imax,typchk,1,&warn );

    if ( ier==-1 ) { return -1; }
    else if ( !ier ) { continue; }
    else if ( ier==2 ) { break; }

    ++ns;
  }

  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param typchk type of checking permformed for edge length (hmax or MMG3D_LLONG criterion).
 * \return -1 if failed.
 * \return number of new points.
 *
 * Analyze tetra and split on geometric criterion.
 *
 */
static int
MMG3D_anatets_iso(MMG5_pMesh mesh,MMG5_pSol met,char typchk) {
  MMG5_pTetra   pt;
  MMG5_pPoint   ppt,p1,p2;
  MMG5_Tria     ptt,ptt2;
  MMG5_xTetra   *pxt;
  MMG5_xPoint   *pxp;
  MMG5_Bezier   pb,pb2;
  MMG5_Hash     hash;
  double        o[3],no[3],to[3],dd;
  int           vx[6],k,ip,ic,it,nap,nc,ni,ne,ns,ip1,ip2,ier;
  char          i,j,j2,ia,i1,i2,ifac;
  static double uv[3][2] = { {0.5,0.5}, {0.,0.5}, {0.5,0.} };
  static char   mmgWarn = 0, mmgWarn2 = 0;

  /** 1. analysis of boundary elements */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) ) return -1;
  ns = nap = 0;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) || !pt->xt )  continue;

    /* check boundary face cut w/r Hausdorff or hmax */
    pt->flag = 0;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++){
      if ( pxt->ftag[i] & MG_REQ )     continue;
      if ( !(pxt->ftag[i] & MG_BDY) )  continue;

      /* virtual triangle */
      MMG5_tet2tri(mesh,k,i,&ptt);

      if ( typchk == 1 ) {
        if ( !MG_GET(pxt->ori,i) ) { continue; }
      }

      if ( !MMG3D_chkbdyface(mesh,met,k,pt,pxt,i,&ptt,typchk) ) { continue; }

      if ( !pt->flag )  continue;
      ns++;

      /* geometric support */
      ier = MMG5_bezierCP(mesh,&ptt,&pb,MG_GET(pxt->ori,i));
      assert(ier);

      /* scan edges in face to split */
      for (j=0; j<3; j++) {
        ia = MMG5_iarf[i][j];
        if ( !MG_GET(pt->flag,ia) )  continue;
        if ( pxt->tag[ia] & MG_REQ )  continue;
        i1  = MMG5_iare[ia][0];
        i2  = MMG5_iare[ia][1];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        ip  = MMG5_hashGet(&hash,ip1,ip2);
        if ( ip > 0 && !(ptt.tag[j] & MG_GEO) )  continue;

        ier = MMG3D_bezierInt(&pb,&uv[j][0],o,no,to);
        assert(ier);
        /* new point along edge */
        if ( !ip ) {
          ip = MMG3D_newPt(mesh,o,MG_BDY);
          if ( !ip ) {
            /* reallocation of point table */
            MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                                 fprintf(stderr,"\n  ## Error: %s: unable to"
                                         " allocate a new point.\n",__func__);
                                 MMG5_INCREASE_MEM_MESSAGE();
                                 MMG3D_delPatternPts(mesh,hash);return -1;
                                 ,o,MG_BDY);
            // Now pb->p contain a wrong memory address.
            pb.p[0] = &mesh->point[ptt.v[0]];
            pb.p[1] = &mesh->point[ptt.v[1]];
            pb.p[2] = &mesh->point[ptt.v[2]];
          }
          if ( !MMG5_hashEdge(mesh,&hash,ip1,ip2,ip) )  return -1;
          ppt = &mesh->point[ip];

          assert ( met );
          if ( met->m ) {
            if ( typchk == 1 && (met->size>1) )
              ier = MMG3D_intmet33_ani(mesh,met,k,ia,ip,0.5);
            else
              ier = MMG5_intmet(mesh,met,k,ia,ip,0.5);

            if ( !ier ) {
              if ( !mmgWarn ) {
                fprintf(stderr,"\n  ## Error: %s: unable to interpolate at least"
                        " 1 metric.\n",__func__);
                mmgWarn = 1;
              }
              return -1;
            }
            else if ( ier < 0 ) {
              MMG3D_delPt(mesh,ip);
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

          if ( mesh->info.fem<typchk ) {
            if ( MG_EDG(ptt.tag[j]) && !(ptt.tag[j] & MG_NOM) ) {
              /* Update the second normal and the tangent at point ip if the edge
               * is shared by 2 faces (if anatet4 is not called, 1 tetra may have
               * 2 bdry faces) */
              ifac = (MMG5_ifar[ia][0] == i)? MMG5_ifar[ia][1] : MMG5_ifar[ia][0];
              if ( pxt->ftag[ifac] & MG_BDY ) {
                j2   = MMG5_iarfinv[ifac][ia];

                /* Compute tangent and normal with respect to the face ifac */
                /* virtual triangle */
                MMG5_tet2tri(mesh,k,ifac,&ptt2);

                /* geometric support */
                ier = MMG5_bezierCP(mesh,&ptt2,&pb2,MG_GET(pxt->ori,ifac));
                assert(ier);

                ier = MMG3D_bezierInt(&pb2,&uv[j2][0],o,no,to);
                assert(ier);

                if ( !MMG3D_storeGeom(ppt,pxp,no) ) continue;
              }
            }
          }
          nap++;
        }
        else if ( MG_EDG(ptt.tag[j]) && !(ptt.tag[j] & MG_NOM) ) {
          ppt = &mesh->point[ip];
          assert(ppt->xp);
          pxp = &mesh->xpoint[ppt->xp];

          if ( !MMG3D_storeGeom(ppt,pxp,no) ) continue;
        }
      }
    }
  }

  if ( !ns ) {
    MMG5_DEL_MEM(mesh,hash.item);
    return ns;
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
        MMG5_tet2tri(mesh,k,i,&ptt);

      for (j=0; j<3; j++) {
        ia  = MMG5_iarf[i][j];
        if ( MG_GET(pt->flag,ia) )                continue;
        if ( pt->xt && (pxt->tag[ia] & MG_REQ) )  continue;
        else if ( ptt.tag[j] & MG_REQ )           continue;
        ip1 = pt->v[MMG5_iare[ia][0]];
        ip2 = pt->v[MMG5_iare[ia][1]];
        ip  = MMG5_hashGet(&hash,ip1,ip2);
        if ( ip > 0 ) {

          MG_SET(pt->flag,ia);
          nc++;
          /* ridge on a boundary face */
          if ( !(ptt.tag[j] & MG_GEO) && !(ptt.tag[j] & MG_NOM) )  continue;
          ppt = &mesh->point[ip];
          assert(ppt->xp);
          pxp = &mesh->xpoint[ppt->xp];
          if ( pt->xt )  ier = MMG5_bezierCP(mesh,&ptt,&pb,MG_GET(pxt->ori,i));
          else  ier = MMG5_bezierCP(mesh,&ptt,&pb,1);
          assert(ier);

          ier = MMG3D_bezierInt(&pb,&uv[j][0],o,no,to);
          assert(ier);

          dd = no[0]*pxp->n1[0]+no[1]*pxp->n1[1]+no[2]*pxp->n1[2];
          if ( dd > 1.0-MMG5_EPS ) continue;

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
          vx[ia] = MMG5_hashGet(&hash,pt->v[i],pt->v[j]);
          if ( vx[ia] > 0 ) {
            MG_SET(pt->flag,ia);
            if ( mesh->point[vx[ia]].flag > 2 )  ic = 1;
          }
        }
      }
      if ( !pt->flag )  continue;

      switch (pt->flag) {
      case 1: case 2: case 4: case 8: case 16: case 32:
        ier = MMG3D_split1_sim(mesh,met,k,vx);
        break;
      case 48: case 24: case 40: case 6: case 34: case 36:
      case 20: case 5: case 17: case 9: case 3: case 10:
        ier =MMG3D_split2sf_sim(mesh,met,k,vx);
        break;
      case 33: case 18: case 12:
        ier = MMG3D_split2_sim(mesh,met,k,vx);
        break;
      case 11: case 21: case 38: case 56:
        ier = MMG3D_split3_sim(mesh,met,k,vx);
        break;
      case 7: case 25: case 42: case 52:
        ier =MMG3D_split3cone_sim(mesh,met,k,vx);
        break;
      case 35: case 19: case 13: case 37: case 22: case 28: case 26:
      case 14: case 49: case 50: case 44: case 41:
        ier = MMG3D_split3op_sim(mesh,met,k,vx);
        break;
      case 23: case 29: case 53: case 60: case 57: case 58:
      case 27: case 15: case 43: case 39: case 54: case 46:
        ier = MMG3D_split4sf_sim(mesh,met,k,vx);
        break;
      case 30: case 45: case 51:
        ier = MMG3D_split4op_sim(mesh,met,k,vx);
        break;
      case 62: case 61: case 59: case 55: case 47: case 31:
        ier = MMG3D_split5_sim(mesh,met,k,vx);
        break;
      case 63:
        ier = MMG3D_split6_sim(mesh,met,k,vx);
        break;
      }
      if ( ier )  continue;

      ni++;
      if ( ic == 0 && MMG3D_dichoto(mesh,met,k,vx) ) {
        for (ia=0; ia<6; ia++)
          if ( vx[ia] > 0 )  mesh->point[vx[ia]].flag++;
      }
      else {
        if ( it < 20 ) {
          for (ia=0,i=0; i<3; i++) {
            for (j=i+1; j<4; j++,ia++) {
              if ( vx[ia] > 0 ) {
                p1 = &mesh->point[pt->v[MMG5_iare[ia][0]]];
                p2 = &mesh->point[pt->v[MMG5_iare[ia][1]]];
                ppt = &mesh->point[vx[ia]];
                ppt->c[0] = 0.5 * (p1->c[0] + p2->c[0]);
                ppt->c[1] = 0.5 * (p1->c[1] + p2->c[1]);
                ppt->c[2] = 0.5 * (p1->c[2] + p2->c[2]);
              }
            }
          }
        }
        else {
          if ( it==20 && (mesh->info.ddebug || mesh->info.imprim > 5) ) {
            if ( !mmgWarn2 ) {
              fprintf(stderr,"\n  ## Warning: %s: surfacic pattern: unable to find"
                      " a valid split for at least 1 point. Point(s) deletion.",
                      __func__ );
              mmgWarn2 = 1;
            }
          }
          for (ia=0,i=0; i<3; i++) {
            for (j=i+1; j<4; j++,ia++) {
              if ( vx[ia] > 0 ) {
                if ( !MMG5_hashUpdate(&hash,pt->v[MMG5_iare[ia][0]],
                                       pt->v[MMG5_iare[ia][1]],-1) ) {
                  fprintf(stderr,"\n  ## Error: %s: unable to delete point"
                          " idx along edge %d %d.\n",
                          __func__,MMG3D_indPt(mesh,pt->v[MMG5_iare[ia][0]]),
                          MMG3D_indPt(mesh,pt->v[MMG5_iare[ia][1]]));
                  MMG5_DEL_MEM(mesh,hash.item);
                  return -1;
                }
                MMG3D_delPt(mesh,vx[ia]);
              }
            }
          }
        }
      }
    }
    nc += ni;
  }
  while( ni > 0 && ++it < 40 );

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
          vx[ia] = MMG5_hashGet(&hash,pt->v[i],pt->v[j]);
          assert(vx[ia]);
        }
      }
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      if ( ! MMG5_split1(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;
    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      if ( ! MMG5_split2sf(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 33: case 18: case 12: /* 2 opposite edges split */
      if ( ! MMG5_split2(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 11: case 21: case 38: case 56: /* 3 edges on the same faces splitted */
      if ( ! MMG5_split3(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      if ( ! MMG5_split3cone(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 35: case 19: case 13: case 37: case 22: case 28: case 26:
    case 14: case 49: case 50: case 44: case 41: /* 3 edges on opposite configuration splitted */
      if ( ! MMG5_split3op(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 23: case 29: case 53: case 60: case 57: case 58:
    case 27: case 15: case 43: case 39: case 54: case 46: /* 4 edges with 3 lying on the same face splitted */
      if ( ! MMG5_split4sf(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

      /* 4 edges with no 3 lying on the same face splitted */
    case 30: case 45: case 51:
      if ( ! MMG5_split4op(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 62: case 61: case 59: case 55: case 47: case 31: /* 5 edges split */
      if ( ! MMG5_split5(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;

    case 63: /* 6 edges split */
      if ( ! MMG5_split6(mesh,met,k,vx,typchk-1) ) return -1;
      ns++;
      break;
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"       %7d elements splitted\n",nap);

  MMG5_DEL_MEM(mesh,hash.item);
  return nap;
}

static int (*MMG3D_anatets)(MMG5_pMesh mesh,MMG5_pSol met,char typchk);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the tetrahedron with multiple boundary faces (to be swapped)
 * \param metRidTyp metric storage (classic or special)
 * \param ifac face of the tetra \a k that give the best results for the swap23
 * \param conf0 detected configuration for the swap23 of the tetra \a k
 * \param adj neighbour of the tetra k through the face \a ifac (4*k1+ifac1)
 * \param conf1 detected configuration for the swap23 of the tetra \a adj/4
 *
 * \return 0 if failed (too bad quality), the type of operator that creates the
 * best worst quality otherwise (1 if split4bar, 2 if swap23).
 *
 * Simulation of the swap23 and of the split at its barycenter of a tetra when
 * more than 1 boundary face. The quality of the worst created element is
 * computed for both operators and we return the identifier of the operator that
 * give the best results. If the swap23 is choosen, we fill the needed info to
 * perform it (index of the face and tetra that are choosen to swap) and
 * configuration of both tetra.
 *
 */

static int MMG3D_anatet4_sim(MMG5_pMesh mesh,MMG5_pSol met,int k,char metRidTyp,
                             int *ifac,int* conf0,int *adj,int *conf1) {
  MMG5_pTetra          pt,pt1,ptnew;
  MMG5_pxTetra         pxt0,pxt1;
  MMG5_pPoint          ppt,ppt0;
  double               calold0,calold,calnew,calnew0,calnew1,calnew2,calnew3;
  double               worst_split4bar_cal,worst_swap_cal,cb[4];
  int                  loc_conf0,loc_conf1,k1,*adja;
  int                  nbdy,i,j0,j1,np;
  unsigned char        tau0[4],tau1[4];

  pt     = &mesh->tetra[k];
  calold0 = pt->qual;

  assert ( pt->xt );

  pxt0 = &mesh->xtetra[pt->xt];

  ptnew = &mesh->tetra[0];

  /** Step 1: test the split4bar */
  if ( !mesh->info.noinsert ) {
    ppt0 = &mesh->point[0];
    memset(ppt0,0,sizeof(MMG5_Point));

    for (i=0; i<4; i++) {
      ppt   = &mesh->point[pt->v[i]];
      ppt0->c[0] += ppt->c[0];
      ppt0->c[1] += ppt->c[1];
      ppt0->c[2] += ppt->c[2];
    }
    ppt0->c[0] *= 0.25;
    ppt0->c[1] *= 0.25;
    ppt0->c[2] *= 0.25;

    cb[0] = 0.25; cb[1] = 0.25;  cb[2] = 0.25;  cb[3] = 0.25;

    if ( met->m ) {
      if ( !metRidTyp && met->size > 1 )
        MMG5_interp4bar33_ani(mesh,met,k,0,cb);
      else
        MMG5_interp4bar(mesh,met,k,0,cb);
    }

    memcpy(ptnew,pt,sizeof(MMG5_Tetra));
    ptnew->v[0] = 0;
    if ( (!metRidTyp) && met->m && met->size>1 )
      calnew0 = MMG5_caltet33_ani(mesh,met,ptnew);
    else
      calnew0 = MMG5_orcal(mesh,met,0);

    ptnew->v[0] = pt->v[0];
    ptnew->v[1] = 0;
    if ( (!metRidTyp) && met->m && met->size>1 )
      calnew1 = MMG5_caltet33_ani(mesh,met,ptnew);
    else
      calnew1 = MMG5_orcal(mesh,met,0);

    ptnew->v[1] = pt->v[1];
    ptnew->v[2] = 0;
    if ( (!metRidTyp) && met->m && met->size>1 )
      calnew2 = MMG5_caltet33_ani(mesh,met,ptnew);
    else
      calnew2 = MMG5_orcal(mesh,met,0);

    ptnew->v[2] = pt->v[2];
    ptnew->v[3] = 0;
    if ( (!metRidTyp) && met->m && met->size>1 )
      calnew3 = MMG5_caltet33_ani(mesh,met,ptnew);
    else
      calnew3 = MMG5_orcal(mesh,met,0);

    worst_split4bar_cal = MG_MIN(MG_MIN(calnew0,calnew1),MG_MIN(calnew2,calnew3));
  }
  else
    worst_split4bar_cal = 0.;

  /** Step 2: test the swap23 */
  *ifac          = -1;
  *adj           = -1;
  if ( !mesh->info.noswap ) {
    worst_swap_cal = 0.;

    for (j0=0; j0<4; j0++) {
      if ( pxt0->ftag[j0] & MG_BDY ) continue;

      /** Neighbouring element with which we will try to swap */
      adja = &mesh->adja[4*(k-1)+1];
      k1   = adja[j0]/4;
      j1   = adja[j0]%4;

      assert(k1);

      /* Search in which configurations are the tetrahedra (default is case 0-0)
       *
       *           3                    2------------- 0
       *         ,/|`\                  |`\          /|
       *       ,/  |  `\                |  `\       /.|
       *     ,/    '.   `\              '.   `\    / |
       *   ,/       |     `\             |     `\ / .|
       * ,/         |       `\           |       /\.|
       * 0-----------'.--------2          '     /  3
       * `\.         |      ,/            |    / ,/
       *    `\.      |    ,/              |   /,/
       *       `\.   '. ,/                '. ,/
       *          `\. |/                   |/
       *             `1                    1
       */

      /* k may be in configuration 0, 3, 6 or 9. Default is case 0 */
      loc_conf0 = 3*j0;

      switch(loc_conf0) {
      case 3:
        tau0[0] = 1; tau0[1] = 0; tau0[2] = 3; tau0[3] = 2;
        break;
      case 6:
        tau0[0] = 2; tau0[1] = 0; tau0[2] = 1; tau0[3] = 3;
        break;
      case 9:
        tau0[0] = 3; tau0[1] = 0; tau0[2] = 2; tau0[3] = 1;
        break;
      default:
        tau0[0] = 0; tau0[1] = 1; tau0[2] = 2; tau0[3] = 3;
      }

      /* k1 may be in configuration j1, j1+1, j1+2 */
      pt1 = &mesh->tetra[k1];

      if ( pt1->tag & MG_REQ ) continue;

      assert(pt->ref == pt1->ref);
      for ( i=0; i<3; ++i )
        if ( pt->v[MMG5_idir[j0][0]] == pt1->v[MMG5_idir[j1][i]] ) break;

      assert(i<3);
      loc_conf1 = 3*j1+i;

      switch(loc_conf1) {
      case 1:
        tau1[0] = 0; tau1[1] = 2; tau1[2] = 3; tau1[3] = 1;
        break;
      case 2:
        tau1[0] = 0; tau1[1] = 3; tau1[2] = 1; tau1[3] = 2;
        break;
      case 3:
        tau1[0] = 1; tau1[1] = 0; tau1[2] = 3; tau1[3] = 2;
        break;
      case 4:
        tau1[0] = 1; tau1[1] = 3; tau1[2] = 2; tau1[3] = 0;
        break;
      case 5:
        tau1[0] = 1; tau1[1] = 2; tau1[2] = 0; tau1[3] = 3;
        break;
      case 6:
        tau1[0] = 2; tau1[1] = 0; tau1[2] = 1; tau1[3] = 3;
        break;
      case 7:
        tau1[0] = 2; tau1[1] = 1; tau1[2] = 3; tau1[3] = 0;
        break;
      case 8:
        tau1[0] = 2; tau1[1] = 3; tau1[2] = 0; tau1[3] = 1;
        break;
      case 9:
        tau1[0] = 3; tau1[1] = 0; tau1[2] = 2; tau1[3] = 1;
        break;
      case 10:
        tau1[0] = 3; tau1[1] = 2; tau1[2] = 1; tau1[3] = 0;
        break;
      case 11:
        tau1[0] = 3; tau1[1] = 1; tau1[2] = 0; tau1[3] = 2;
        break;
      default:
        tau1[0] = 0; tau1[1] = 1; tau1[2] = 2; tau1[3] = 3;
      }

      /** Do not choose a config that creates a tetra with more than 2 bdries */
      if ( pt1->xt ) {
        pxt1 = &mesh->xtetra[pt1->xt];

        nbdy = 0;
        if ( pxt1->ftag[tau1[1]] ) ++nbdy;
        if ( pxt0->ftag[tau0[1]] ) ++nbdy;
        if ( nbdy > 1 ) continue;

        nbdy = 0;
        if ( pxt1->ftag[tau1[3]] ) ++nbdy;
        if ( pxt0->ftag[tau0[2]] ) ++nbdy;
        if ( nbdy > 1 ) continue;

        nbdy = 0;
        if ( pxt1->ftag[tau1[2]] ) ++nbdy;
        if ( pxt0->ftag[tau0[3]] ) ++nbdy;
        if ( nbdy > 1 ) continue;

      }

      /** Test volume of the 3 created tets */
      calold = MG_MIN(calold0,pt1->qual);

      ptnew = &mesh->tetra[0];
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));
      np    = pt1->v[tau1[0]];

      ptnew->v[tau0[1]] = np;
      if ( (!metRidTyp) && met->m && met->size>1 )
        calnew0 = MMG5_caltet33_ani(mesh,met,ptnew);
      else
        calnew0 = MMG5_orcal(mesh,met,0);
      if ( calnew0 < MMG5_NULKAL ) continue;

      ptnew->v[tau0[1]] = pt->v[tau0[1]];
      ptnew->v[tau0[2]] = np;
      if ( (!metRidTyp) && met->m && met->size>1 )
        calnew1 = MMG5_caltet33_ani(mesh,met,ptnew);
      else
        calnew1 = MMG5_orcal(mesh,met,0);
      if ( calnew1 < MMG5_NULKAL ) continue;

      ptnew->v[tau0[2]] = pt->v[tau0[2]];
      ptnew->v[tau0[3]] = np;
      if ( (!metRidTyp) && met->m && met->size>1 )
        calnew2 = MMG5_caltet33_ani(mesh,met,ptnew);
      else
        calnew2 = MMG5_orcal(mesh,met,0);
      if ( calnew2 < MMG5_NULKAL ) continue;

      calnew = MG_MIN(calnew0,MG_MIN(calnew1,calnew2));

      if ( calold < MMG5_EPSOK ) {
        if ( calnew < calold ) continue;
      }
      else if ( calnew <= MMG5_EPSOK ) continue;

      if ( calnew > worst_swap_cal ) {
        worst_swap_cal = calnew;
        *ifac          = j0;
        *adj           = adja[j0];
        *conf0         = loc_conf0;
        *conf1         = loc_conf1;
      }
    }
  }
  else
    worst_swap_cal = 0.;

  if ( *ifac < 0 ) {
    /* No valid config for the swap23 */
    if ( worst_split4bar_cal < MMG5_EPSOK ) return 0;
    return 1;
  }

  assert ( *adj > 0 );
  if ( worst_swap_cal < worst_split4bar_cal ) {
    if ( worst_split4bar_cal < MMG5_EPSOK ) return 0;
    return  1;
  }

  if ( worst_swap_cal < MMG5_EPSOK ) return 0;

  return 2;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param nf number of swap performed.
 * \param typchk type of checking permformed.
 * \return -1 if failed, number of new points otherwise.
 *
 * Split tetra into 4 when more than 1 boundary face or if 4 boundary vertices.
 *
 */
static int MMG5_anatet4(MMG5_pMesh mesh, MMG5_pSol met,int *nf, char typchk) {
  MMG5_pTetra  pt;
  MMG5_pPoint  ppt;
  MMG5_pxTetra pxt;
  int          k,ns,ier,conf0,conf1,adj,ifac,id_op;
  char         nbdy,j;
#ifndef NDEBUG
  static char mmgWarn=0;
#endif

  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    nbdy = 0;
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      for (j=0; j<4; j++)
        if ( ( pxt->ftag[j] & MG_BDY ) && (!(pxt->ftag[j] & MG_PARBDY)) )  nbdy++;
    }
    if ( nbdy > 1 ) {
      id_op = MMG3D_anatet4_sim(mesh,met,k,typchk-1,&ifac,&conf0,&adj,&conf1);
      if ( !id_op ) {
#ifndef NDEBUG
        if ( (!(mesh->info.noswap && mesh->info.noinsert)) && !mmgWarn ) {
          mmgWarn=1;
          printf("\n  ## Warning: %s: unable to swap or split at least 1 tetra"
                 " with multiple boundary faces.\n",__func__);
        }
#endif
        continue;
      }

      if ( id_op==1 ) {
        /* Split */
        ier  = MMG5_split4bar(mesh,met,k,typchk-1);
        if ( !ier ) return -1;
        ns++;
      }
      else {
        assert ( id_op==2 );

        /* Swap */
        ier = MMG3D_swap23(mesh,met,k,typchk-1,ifac,conf0,adj,conf1);
        if ( ier < 0 ) {
          return -1;
        }
        else if ( ier ) ++(*nf);
      }
    }
    else {
      nbdy = 0;
      for (j=0; j<4; j++) {
        ppt = &mesh->point[pt->v[j]];
        if ( (ppt->tag & MG_BDY) && (!(ppt->tag & MG_PARBDY)) )  nbdy++;
      }
      if ( nbdy == 4 ) {
        if ( !mesh->info.noinsert ) {
          ier  = MMG5_split4bar(mesh,met,k,typchk-1);
          if ( !ier ) return -1;
          ns++;
        }
      }
    }
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     boundary elements: %7d splitted %7d swapped\n",ns,*nf);
  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param nf number of swap performed.
 * \param typchk type of checking permformed.
 * \return -1 if failed, number of new points otherwise.
 *
 * Split tetra into 4 when more than 1 boundary face or if 4 boundary vertices.
 *
 */
static int MMG5_anatet4rid(MMG5_pMesh mesh, MMG5_pSol met,int *nf, char typchk) {
  MMG5_pTetra  pt;
  MMG5_pPoint  ppt;
  int          k,ns,ier;
  char         nrid,j;

  ns = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    nrid = 0;

    for (j=0; j<4; j++) {
      ppt = &mesh->point[pt->v[j]];
      if ( (ppt->tag & MG_GEO) )  nrid++;
    }
    if ( nrid == 4 ) {
      if ( !mesh->info.noinsert ) {
        ier  = MMG5_split4bar(mesh,met,k,typchk-1);
        if ( !ier ) return -1;
        ns++;
      }
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     boundary elements: %7d splitted\n",ns);
  return ns;
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
int MMG5_anatet(MMG5_pMesh mesh,MMG5_pSol met,char typchk, int patternMode) {
  int     ier,nc,ns,nf,nnc,nns,nnf,it,minit,maxit,lastit;

  /* pointer toward the suitable anatets function */
  if ( met->m && met->size==6 ) {
    /* if the aniso metric is not compatible with the geometry, the non
     * conformal surface operators may create spurious ridges */
    MMG3D_anatets = MMG3D_anatets_ani;
  }
  else {
    MMG3D_anatets = MMG3D_anatets_iso;
  }

  /* analyze tetras : initial splitting */
  nns = nnc = nnf = it = 0;
  lastit = 0;
  minit = 3;
  maxit = 6;
  mesh->gap = 0.5;
  do {
    if ( typchk==2 && lastit==1 )  ++mesh->info.fem;


    /* split or swap tetra with more than 2 bdry faces */
    nf = ier = 0;
    if ( mesh->info.fem == typchk ) {
      ier = MMG5_anatet4(mesh,met,&nf,typchk);
      if ( ier < 0 )  return 0;
    }
    else if ( (met->size==6) && (typchk == 1) && lastit ) {
      ier = MMG5_anatet4rid(mesh,met,&nf,typchk);
      if ( ier < 0 )  return 0;
    }
    else ier = 0;
    ns = ier;

    /* memory free */
    if ( mesh->adja )
      MMG5_DEL_MEM(mesh,mesh->adja);

    if ( !mesh->info.noinsert ) {
      /* analyze surface tetras */
      ier = MMG3D_anatets(mesh,met,typchk);

      if ( ier < 0 ) {
        fprintf(stderr,"\n  ## Unable to complete surface mesh. Exit program.\n");
        return 0;
      }
      ns += ier;

      if ( patternMode ) {
        if ( mesh->adja ) { MMG5_DEL_MEM(mesh,mesh->adja); }

        /* analyze internal tetras */
        ier = MMG5_anatetv(mesh,met,typchk);
        if ( ier < 0 ) {
          fprintf(stderr,"\n  ## Unable to complete volume mesh. Exit program.\n");
          return 0;
        }
        ns += ier;
      }
    }

    if ( (!mesh->adja) && (!MMG3D_hashTetra(mesh,1)) ) {
      fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
      return 0;
    }

    /* collapse short edges */
    if ( !mesh->info.noinsert ) {
      nc = MMG5_coltet(mesh,met,typchk);
      if ( nc < 0 ) {
        fprintf(stderr,"\n  ## Unable to collapse mesh. Exiting.\n");
        return 0;
      }
    }
    else  nc = 0;

    /* attempt to swap */
    if ( !mesh->info.noswap ) {
      ier = MMG5_swpmsh(mesh,met,NULL,typchk);
      if ( ier < 0 ) {
        fprintf(stderr,"\n  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
      nf  += ier;

      ier = MMG5_swptet(mesh,met,MMG3D_LSWAPIMPROVE,MMG3D_SWAP06,NULL,typchk,mesh->mark-2);
      if ( ier < 0 ) {
        fprintf(stderr,"\n  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
      nf += ier;
    }
    else nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nf > 0 ){
#ifndef PATTERN
      fprintf(stdout,"                   ");
#endif
      fprintf(stdout,"     %8d splitted, %8d collapsed, %8d swapped\n",ns,nc,nf);
    }

    if ( it > minit-1 && ( !(ns+nc) || (abs(nc-ns) < 0.1 * MG_MAX(nc,ns)) ) ) {
      ++lastit;
      if ( it > minit && lastit>2 ) break;
    }
    else if ( it+2 >= maxit ) {
      /* Last iteration because we reach the maximal number of iter */
      ++lastit;
    }
    else if ( lastit ) {
      /* Avoid the incrementation of mesh->info.fem if we have detected a last
         iteration but anatet4 leads to have nc, nf or ns != 0 so we perform a last
         iter */
      ++lastit;
    }
  }
  while ( ++it < maxit && (ns+nc+nf > 0 || lastit<3) );

  if ( mesh->info.imprim > 0 ) {
    if ( (abs(mesh->info.imprim) < 5 || mesh->info.ddebug ) && nns+nnc > 0 ) {
#ifndef PATTERN
      fprintf(stdout,"                   ");
#endif
      fprintf(stdout, "     %8d splitted, %8d collapsed, %8d swapped, %d iter.\n",nns,nnc,nnf,it);
    }
  }

  return 1;
}
