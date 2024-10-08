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
 * \file mmgs/mmgs1.c
 * \brief Perform surface mesh adaptation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmgs_private.h"
#include "mmgsexterns_private.h"
#include "mmgexterns_private.h"
#include "inlined_functions_private.h"

extern int8_t ddb;

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k element index.
 * \param vx pointer toward table of edges to split.
 * \return 1.
 *
 * Find acceptable position for splitting.
 *
 */
int MMGS_dichoto(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,MMG5_int *vx) {
  MMG5_pTria   pt;
  MMG5_pPoint  pa,pb,ps;
  double       o[3][3],p[3][3];
  float        to,tp,t;
  MMG5_int     ia,ib;
  int          i1,i2,ier,it,maxit;
  int8_t       i,j;

  pt = &mesh->tria[k];
  /* get point on surface and along segment for edge split */
  for (i=0; i<3; i++) {
    memset(p[i],0,3*sizeof(double));
    memset(o[i],0,3*sizeof(double));
    if ( vx[i] > 0 ) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      ia = pt->v[i1];
      ib = pt->v[i2];
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
  j  = 10;
  do {
    /* compute new position */
    t = 0.5 * (tp + to);
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
        ps->c[2] = o[i][2] + t*(p[i][2] - o[i][2]);
        j=i;
      }
    }
    switch (pt->flag) {
    case 1: case 2: case 4:
      ier = MMGS_split1_sim(mesh,met,k,j,vx);
      break;
    case 7:
      ier = MMGS_split3_sim(mesh,met,k,vx);
      break;
    default:
      ier = MMG5_split2_sim(mesh,met,k,vx);
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
    for (i=0; i<3; i++) {
      if ( vx[i] > 0 ) {
        ps = &mesh->point[vx[i]];
        ps->c[0] = o[i][0] + t*(p[i][0] - o[i][0]);
        ps->c[1] = o[i][1] + t*(p[i][1] - o[i][1]);
        ps->c[2] = o[i][2] + t*(p[i][2] - o[i][2]);
      }
    }
  }

  /* For very ill-shaped elements we can have no valid position */
  assert ( j< 10 );
  switch (pt->flag) {
  case 1: case 2: case 4:
    ier = MMGS_split1_sim(mesh,met,k,j,vx);
    break;
  case 7:
    ier = MMGS_split3_sim(mesh,met,k,vx);
    break;
  default:
    ier = MMG5_split2_sim(mesh,met,k,vx);
    break;
  }
  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param iel index of the starting triangle.
 * \param ia local index of the edge to split in \a k.
 * \param ip index of the point that we try to create.
 * \return 1 if success, 0 otherwise.
 *
 * Find acceptable position for MMG5_split1b, starting from point ip.
 *
 */
int MMGS_dichoto1b(MMG5_pMesh mesh, MMG5_pSol met, MMG5_int iel, int ia, MMG5_int ip) {
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1,ppt;
  MMG5_int     np,nq;
  int          it,maxit,i1,i2;
  double       m[3],o[3],tp,to,t;
  int8_t       ier;

  pt  = &mesh->tria[iel];

  i1 = MMG5_inxt2[ia];
  i2 = MMG5_inxt2[i1];
  np = pt->v[i1];
  nq = pt->v[i2];
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

    ier = MMGS_simbulgept(mesh,met,iel,ia,ip);
    if ( ier )
      to = t;
    else
      tp = t;
  }
  while ( ++it < maxit );
  if ( !ier )  t = to;

  /* For very ill-shaped elements, we can have no valid position */
  ppt->c[0] = m[0] + t*(o[0]-m[0]);
  ppt->c[1] = m[1] + t*(o[1]-m[1]);
  ppt->c[2] = m[2] + t*(o[2]-m[2]);

  return  MMGS_simbulgept(mesh,met,iel,ia,ip) ;
}

/* check if edge need to be split and return a binary coding the numbers of the
 * edges of tria iel that should be split according to a hausdorff distance
 * criterion */
int chkedg(MMG5_pMesh mesh,MMG5_int iel) {
  MMG5_pTria    pt;
  MMG5_pPoint   p[3];
  MMG5_pPar     par;
  double        n[3][3],t[3][3],nt[3],c1[3],c2[3],*n1,*n2,t1[3],t2[3];
  double        ps,ps2,cosn,ux,uy,uz,ll,li,dd,hausd,hmax;
  int           l,isloc;
  int8_t        i,i1,i2;
  static int8_t mmgWarn0 = 0, mmgWarn1 = 0;

  pt   = &mesh->tria[iel];
  p[0] = &mesh->point[pt->v[0]];
  p[1] = &mesh->point[pt->v[1]];
  p[2] = &mesh->point[pt->v[2]];

  /* normal recovery */
  for (i=0; i<3; i++) {
    if ( MS_SIN(p[i]->tag) ) {
      MMG5_nortri(mesh,pt,n[i]);
    }
    else if ( MG_EDG(p[i]->tag) ) {
      MMG5_nortri(mesh,pt,nt);
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

  /* analyze edges */
  for (i=0; i<3; i++) {
    if ( pt->tag[i] & MG_REQ) continue;

    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];

    /* local parameters */
    hmax   = mesh->info.hmax;
    hausd  = mesh->info.hausd;
    isloc  = 0;
    for (l=0; l<mesh->info.npar; l++) {
      par = &mesh->info.par[l];
      if ( /*((par->elt == MMG5_Vertex) &&
            ( (p[i1]->ref == par->ref ) || (p[i2]->ref == par->ref) ))
            || */ ((par->elt == MMG5_Triangle) && (pt->ref == par->ref) ) ) {
        if ( !isloc ) {
          hmax  = par->hmax;
          hausd = par->hausd;
          isloc = 1;
        }
        else {
          hausd = MG_MIN(par->hausd,hausd);
          hmax  = MG_MIN(par->hmax,hmax);
        }
      }
    }

    /* check length */
    ux = p[i2]->c[0] - p[i1]->c[0];
    uy = p[i2]->c[1] - p[i1]->c[1];
    uz = p[i2]->c[2] - p[i1]->c[2];
    ll = ux*ux + uy*uy + uz*uz;
    if ( ll > hmax*hmax ) {
      MG_SET(pt->flag,i);
      continue;
    }
    else if ( !MG_EDG(pt->tag[i]) && p[i1]->tag > MG_NOTAG && p[i2]->tag > MG_NOTAG ) {
      MG_SET(pt->flag,i);
      continue;
    }

    /* Hausdorff w/r tangent direction */
    if ( MG_EDG(pt->tag[i]) ) {
      if ( MS_SIN(p[i1]->tag) ) {
        li = 1.0 / sqrt(ll);
        t1[0] = li*ux;
        t1[1] = li*uy;
        t1[2] = li*uz;
      }
      else{
        if(!((p[i1]->tag & MG_NOM) ||  MG_EDG(p[i1]->tag) ) ) {
          //  if(t[i1][0] > 10) {
          if ( !mmgWarn0 ) {
            fprintf(stderr,"\n  ## Warning: %s: a- at least 1 geometrical"
                    " problem\n",__func__);
            mmgWarn0 = 1;
          }
          return 0;
        }
        memcpy(t1,t[i1],3*sizeof(double));
      }

      if ( MS_SIN(p[i2]->tag) ) {
        li = 1.0 / sqrt(ll);
        t2[0] = li*ux;
        t2[1] = li*uy;
        t2[2] = li*uz;
      }
      else{
        if(!((p[i2]->tag & MG_NOM) || MG_EDG(p[i2]->tag) ) ) {
          if ( !mmgWarn1 ) {
            fprintf(stderr,"\n  ## Warning: %s: b- at least 1 geometrical"
                    " problem\n",__func__);
            mmgWarn1 = 1;
          }
          return 0;
        }
        memcpy(t2,t[i2],3*sizeof(double));
      }

      ps = t1[0]*ux + t1[1]*uy + t1[2]*uz;
      ps *= ps;
      cosn = ps/ll ;
      cosn *= (1.0-cosn);
      cosn *= (0.25*ll);
      if ( cosn > hausd*hausd ) {
        MG_SET(pt->flag,i);
        continue;
      }

      ps = -(t2[0]*ux + t2[1]*uy + t2[2]*uz);
      ps *= ps;
      cosn = ps/ll ;
      cosn *= (1.0-cosn);
      cosn *= (0.25*ll);
      if ( cosn > hausd*hausd ) {
        MG_SET(pt->flag,i);
        continue;
      }
    }
    else {
      n1 = n[i1];
      n2 = n[i2];

      ps = ux*n1[0] + uy*n1[1] + uz*n1[2];
      c1[0] = (2.0*p[i1]->c[0] + p[i2]->c[0] - ps*n1[0]) / 3.0 - p[i1]->c[0];
      c1[1] = (2.0*p[i1]->c[1] + p[i2]->c[1] - ps*n1[1]) / 3.0 - p[i1]->c[1];
      c1[2] = (2.0*p[i1]->c[2] + p[i2]->c[2] - ps*n1[2]) / 3.0 - p[i1]->c[2];

      ps = -(ux*n2[0] + uy*n2[1] + uz*n2[2]);
      c2[0] = (2.0*p[i2]->c[0] + p[i1]->c[0] - ps*n2[0]) / 3.0 - p[i2]->c[0];
      c2[1] = (2.0*p[i2]->c[1] + p[i1]->c[1] - ps*n2[1]) / 3.0 - p[i2]->c[1];
      c2[2] = (2.0*p[i2]->c[2] + p[i1]->c[2] - ps*n2[2]) / 3.0 - p[i2]->c[2];

      /* squared cosines */
      ps = c1[0]*ux + c1[1]*uy + c1[2]*uz;
      ps *= ps;
      dd = c1[0]*c1[0] + c1[1]*c1[1] + c1[2]*c1[2];
      cosn  =  ps / (dd*ll);
      cosn *= (1.0-cosn);
      cosn *= (0.25*ll);
      if ( cosn > hausd*hausd ) {
        MG_SET(pt->flag,i);
        continue;
      }

      ps = -c2[0]*ux - c2[1]*uy - c2[2]*uz;
      ps *= ps;
      dd = c2[0]*c2[0]+c2[1]*c2[1]+c2[2]*c2[2];
      cosn  =  ps / (dd*ll);
      cosn *= (1.0-cosn);
      cosn *= (0.25*ll);
      if ( cosn > hausd*hausd ) {
        MG_SET(pt->flag,i);
        continue;
      }
    }
  }

  return pt->flag;
}

/**
 * \param met pointer toward met structure
 * \param typchk type of check to perform: 1 for first stage (adaptation to
 * capture roughly the surface mesh), 2 for second stage of adaptation
 * (rough capture of input metric).
 *
 * Assign functions for computation of edge lengths and tria qualities: they are
 * now used in first stage of adaptation to ensure that, if we enter with a
 * mesh+metric already adapted, we don't delete entierly the work done (massive
 * collapses on planar surfaces for example).
 *
 */
static inline
void MMGS_set_localFunc ( MMG5_pSol met, int8_t typchk,
                          double (**MMGS_lenEdg)(MMG5_pMesh,MMG5_pSol,MMG5_int ,MMG5_int,int8_t),
                          double (**MMGS_caltri)(MMG5_pMesh,MMG5_pSol,MMG5_pTria)) {

  if ( met->m ) {
    if ( (typchk == 1) && (met->size == 6) ) {
      *MMGS_lenEdg = MMG5_lenSurfEdg33_ani;
      *MMGS_caltri = MMG5_caltri33_ani;
    }
    else if ( typchk == 2 ) {
      *MMGS_lenEdg = MMG5_lenSurfEdg;
      *MMGS_caltri = MMG5_calelt;
    }
    else {
      *MMGS_lenEdg = MMG5_lenSurfEdg_iso;
      *MMGS_caltri = MMG5_caltri_iso;
    }
  }

}


static int swpmsh(MMG5_pMesh mesh,MMG5_pSol met,int8_t typchk) {
  MMG5_pTria    pt;
  int           it,maxit;
  MMG5_int      k,ns,nns;
  int8_t        i;

  /* Local function pointer toward the suitable functions to use for edge length
   * and quality computation depending on the adaptation phase */
  double (*MMGS_lenEdg)(MMG5_pMesh mesh,MMG5_pSol sol ,MMG5_int ,MMG5_int, int8_t ) = NULL;
  double (*MMGS_caltri)(MMG5_pMesh mesh,MMG5_pSol sol ,MMG5_pTria pt ) = MMG5_caltri_iso;

  MMGS_set_localFunc ( met, typchk, &MMGS_lenEdg, &MMGS_caltri);

  it = nns = 0;
  maxit = 2;
  mesh->base++;
  do {
    ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

      for (i=0; i<3; i++) {
        if ( MS_SIN(pt->tag[i]) || MG_EDG(pt->tag[i]) )  continue;
        else if ( chkswp(mesh,met,k,i,typchk,
                         MMGS_lenEdg,MMGS_caltri) ) {
          ns += swapar(mesh,k,i);
          break;
        }
      }
    }
    nns += ns;
  }
  while ( ns > 0 && ++it < maxit );
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nns > 0 )
    fprintf(stdout,"     %8" MMG5_PRId " edge swapped\n",nns);

  return nns;
}

/* Analyze triangles and move points to make mesh more uniform */
static int movtri(MMG5_pMesh mesh,MMG5_pSol met,int maxit) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt;
  int           it,ier,ilist;
  MMG5_int      k,base,list[MMG5_TRIA_LMAX+2],nm,ns,nnm;
  int8_t        i;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug )
    fprintf(stdout,"  ** OPTIMIZING MESH\n");

  base = 1;
  for (k=1; k<=mesh->np; k++)  mesh->point[k].flag = base;

  it = nnm = 0;
  do {
    base++;
    nm = ns = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

      for (i=0; i<3; i++) {
        ppt = &mesh->point[pt->v[i]];

        if ( ppt->flag == base || MS_SIN(ppt->tag) || ppt->tag & MG_NOM )
          continue;

        int8_t dummy;
        ilist = MMG5_boulet(mesh,k,i,list,1,&dummy);

        if ( ilist < 1 ) continue;

        if ( MG_EDG(ppt->tag) ) {
          ier = movridpt(mesh,met,list,ilist);
          if ( ier ) {
            ns++;
          }
        }
        else {
          ier = movintpt(mesh,met,list,ilist);
        }

        if ( ier ) {
          nm++;
          ppt->flag = base;
        }
      }
    }
    nnm += nm;
    if ( mesh->info.ddebug )  fprintf(stdout,"     %8" MMG5_PRId " moved, %" MMG5_PRId " geometry\n",nm,ns);
  }
  while ( ++it < maxit && nm > 0);

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && nnm > 0 )
    fprintf(stdout,"     %8" MMG5_PRId " vertices moved, %d iter.\n",nnm,it);

  return nnm;
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
int MMGS_delPatternPts(MMG5_pMesh mesh,MMG5_Hash hash)
{
  MMG5_pTria   pt;
  MMG5_int     vx[3],k;
  int          i,i1,i2;

  /* Delete the useless added points */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )  continue;

    for (i=0; i<3; i++) {
      i1    = MMG5_inxt2[i];
      i2    = MMG5_inxt2[i1];
      vx[i] = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);

      if ( vx[i] > 0 ) {
        MMGS_delPt(mesh,vx[i]);
        if ( !MMG5_hashUpdate(&hash,pt->v[i1],pt->v[i2],0) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to delete point idx"
                  " along edge %" MMG5_PRId " %" MMG5_PRId ".\n", __func__,
                  MMGS_indPt(mesh,pt->v[i1]),
                  MMGS_indPt(mesh,pt->v[i2]));
          MMG5_DEL_MEM(mesh,hash.item);
          return 0;
        }
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 * \param typchk type of check performed depending on the remeshing step
 *
 * \return -1 if fail, the number of split otherwise
 *
 * Analyze triangles and split if needed
 *
 */
static int anaelt(MMG5_pMesh mesh,MMG5_pSol met,int8_t typchk) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt,p1,p2;
  MMG5_Hash     hash;
  MMG5_Bezier   pb;
  MMG5_pxPoint  go;
  double        s,o[3],no[3],to[3],dd,len;
  int           i,ier,it;
  MMG5_int      ns,nc,ni,ic;
  MMG5_int      k,nt,vx[3],j,ip,ip1,ip2;
  int8_t        i1,i2;
  static double uv[3][2] = { {0.5,0.5}, {0.,0.5}, {0.5,0.} };
  static int8_t mmgWarn0=0,mmgWarn1=0,mmgWarn2=0,mmgWarn3=0;

  if ( !MMG5_hashNew(mesh,&hash,mesh->np,3*mesh->np) ) return -1;
  ns = 0;
  s  = 0.5;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) {
      continue;
    }

    pt->flag = 0;

    if ( pt->ref < 0 ) {
      continue;
    }

    /* Required triangle */
    if ( MS_SIN(pt->tag[0]) && MS_SIN(pt->tag[1]) && MS_SIN(pt->tag[2]) )  continue;

    /* check element cut */
    if ( typchk == 1 ) {
      if ( !chkedg(mesh,k) )  continue;
    }
    else if ( typchk == 2 ) {
      for (i=0; i<3; i++) {
        if ( pt->tag[i] & MG_REQ) continue;

        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        len = MMG5_lenSurfEdg(mesh,met,pt->v[i1],pt->v[i2],0);
        if ( !len ) return -1;
        else if ( len > MMGS_LLONG )  {
          MG_SET(pt->flag,i);
        }
      }
      if ( !pt->flag )  continue;
    }
    ns++;

    /* geometric support */
    ier = MMG5_bezierCP(mesh,pt,&pb,1);
    assert(ier);

    /* scan edges to split */
    for (i=0; i<3; i++) {
      if ( !MG_GET(pt->flag,i) )  continue;
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      ip = MMG5_hashGet(&hash,ip1,ip2);

      /* do not compute the point twice except along special edges */
      if ( !MG_EDG(pt->tag[i]) && ip > 0 )  continue;

      /* new point along edge */
      ier = MMGS_bezierInt(&pb,uv[i],o,no,to);
      if ( !ip ) {
        ip = MMGS_newPt(mesh,o,MG_EDG(pt->tag[i]) ? to : no);

        if ( !ip ) {
          /* reallocation of point table */
          MMGS_POINT_REALLOC(mesh,met,ip,mesh->gap,
                              fprintf(stderr,"\n  ## Error: %s: unable to"
                                      " allocate a new point.\n",__func__);
                              MMG5_INCREASE_MEM_MESSAGE();
                              MMGS_delPatternPts( mesh, hash);
                              return -1
                              ,o,MG_EDG(pt->tag[i]) ? to : no);
          // Now pb->p contain a wrong memory address.
          pb.p[0] = &mesh->point[pt->v[0]];
          pb.p[1] = &mesh->point[pt->v[1]];
          pb.p[2] = &mesh->point[pt->v[2]];
        }

        MMG5_hashEdge(mesh,&hash,ip1,ip2,ip);
        p1  = &mesh->point[ip1];
        p2  = &mesh->point[ip2];
        ppt = &mesh->point[ip];

        if ( MG_EDG(pt->tag[i]) ) {
          ++mesh->xp;
          if(mesh->xp > mesh->xpmax){
            /* reallocation of xpoint table */
            MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                               "larger xpoint table",
                               return -1);
          }
          ppt->xp  = mesh->xp;
          ppt->tag = pt->tag[i];
          if ( p1->ref == pt->edg[i] || p2->ref == pt->edg[i] )
            ppt->ref = pt->edg[i];
          ppt->xp  = mesh->xp;
          go = &mesh->xpoint[mesh->xp];
          memset(go->n2,0x00,3*sizeof(double));
          memcpy(go->n1,no,3*sizeof(double));

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
          }
        }
        if ( met->m ) {
          if ( typchk == 1 && (met->size>1) )
            ier = MMGS_intmet33_ani(mesh,met,k,i,ip,s);
          else
            ier = intmet(mesh,met,k,i,ip,s);
        }

        if ( !ier ) {
          if ( !mmgWarn0 ) {
            fprintf(stderr,"\n  ## Error: %s: unable to interpolate at least"
                    " 1 metric.\n",__func__);
            mmgWarn0 = 1;
          }
          return -1;
        }
      }
      else if ( pt->tag[i] & MG_GEO && !(pt->tag[i] & MG_NOM) ) {
        ppt = &mesh->point[ip];
        go  = &mesh->xpoint[ppt->xp];
        memcpy(go->n2,no,3*sizeof(double));

        /* a computation of the tangent with respect to these two normals is possible */
        ppt->n[0] = go->n1[1]*go->n2[2] - go->n1[2]*go->n2[1];
        ppt->n[1] = go->n1[2]*go->n2[0] - go->n1[0]*go->n2[2];
        ppt->n[2] = go->n1[0]*go->n2[1] - go->n1[1]*go->n2[0];
        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
        if ( dd > MMG5_EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          ppt->n[0] *= dd;
          ppt->n[1] *= dd;
          ppt->n[2] *= dd;
        }
        else {
          if ( !mmgWarn1 ) {
            mmgWarn1 = 1;
            fprintf(stderr,"\n  ## Warning: %s: flattened angle around ridge."
                    " Unable to split it.\n",__func__);
          }
          /* Remove the point */
          MMGS_delPt(mesh,ip);
          MMG5_hashUpdate(&hash,ip1,ip2,0);
        }
      }
    }
  }
  if ( !ns ) {
    MMG5_DEL_MEM(mesh,hash.item);
    return ns;
  }

  /* step 2. checking if split by adjacent */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )  continue;
    else if ( pt->flag == 7 )  continue;

    /* geometric support */
    ier = MMG5_bezierCP(mesh,pt,&pb,1);
    assert(ier);
    nc = 0;

    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      if ( !MG_GET(pt->flag,i) ) {
        ip = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( ip > 0 ) {
          MG_SET(pt->flag,i);
          nc++;
          if ( pt->tag[i] & MG_GEO ) {
            /* new point along edge */
            ier = MMGS_bezierInt(&pb,uv[i],o,no,to);
            assert(ier);

            ppt = &mesh->point[ip];
            go  = &mesh->xpoint[ppt->xp];
            memcpy(go->n2,no,3*sizeof(double));

            /* a computation of the tangent with respect to these two normals is possible */
            ppt->n[0] = go->n1[1]*go->n2[2] - go->n1[2]*go->n2[1];
            ppt->n[1] = go->n1[2]*go->n2[0] - go->n1[0]*go->n2[2];
            ppt->n[2] = go->n1[0]*go->n2[1] - go->n1[1]*go->n2[0];
            dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
            if ( dd > MMG5_EPSD2 ) {
              dd = 1.0 / sqrt(dd);
              ppt->n[0] *= dd;
              ppt->n[1] *= dd;
              ppt->n[2] *= dd;
            }
          }
        }
      }
    }
    if ( nc > 0 )  ++ns;
  }
  if ( mesh->info.ddebug && ns ) {
    fprintf(stdout,"     %" MMG5_PRId " analyzed  %" MMG5_PRId " proposed\n",mesh->nt,ns);
    fflush(stdout);
  }

  /** 3. Simulate splitting and delete points leading to invalid configurations */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  it = 1;
  nc = 0;
  j  = 0; // To remove maybe-uninitialized warning on arm
  do {
    ni = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) || pt->ref < 0 ) continue;
      else if ( pt->flag == 0 )  continue;

      vx[0] = vx[1] = vx[2] = 0;
      pt->flag = ic = 0;

      for (i=0; i<3; i++) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];
        vx[i] = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( vx[i] > 0 ) {
          MG_SET(pt->flag,i);
          if ( mesh->point[vx[i]].flag > 2 )  ic = 1;
          j = i;
        }
      }
      if ( !pt->flag )  continue;
      switch (pt->flag) {
      case 1: case 2: case 4:
        ier = MMGS_split1_sim(mesh,met,k,j,vx);
        break;
      case 7:
        ier = MMGS_split3_sim(mesh,met,k,vx);
        break;
      default:
        ier = MMG5_split2_sim(mesh,met,k,vx);
        break;
      }
      if ( ier )  continue;

      ni++;
      if ( ic == 0 && MMGS_dichoto(mesh,met,k,vx) ) {
        for (i=0; i<3; i++)
          if ( vx[i] > 0 )  mesh->point[vx[i]].flag++;
      }
      else {
        if ( it < 20 ) {
          for (i=0; i<3; i++) {
            if ( vx[i] > 0 ) {
              p1 = &mesh->point[pt->v[MMG5_iprv2[i]]];
              p2 = &mesh->point[pt->v[MMG5_inxt2[i]]];
              ppt = &mesh->point[vx[i]];
              ppt->c[0] = 0.5 * (p1->c[0] + p2->c[0]);
              ppt->c[1] = 0.5 * (p1->c[1] + p2->c[1]);
              ppt->c[2] = 0.5 * (p1->c[2] + p2->c[2]);
            }
          }
        }
        else {
          if ( it==20 && (mesh->info.ddebug || mesh->info.imprim > 5) ) {
            if ( !mmgWarn2 ) {
              fprintf(stderr,"\n  ## Warning: %s: surfacic pattern: unable to"
                      " find a valid split for at least 1 point. Point(s)"
                      " deletion.",__func__ );
              mmgWarn2 = 1;
            }
          }
          for (i=0; i<3; i++) {
            if ( vx[i] > 0 ) {
              if ( !MMG5_hashUpdate(&hash,pt->v[MMG5_iprv2[i]],
                                     pt->v[MMG5_inxt2[i]],0) ) {
                fprintf(stderr,"\n  ## Error: %s: unable to delete point"
                        " idx along edge %" MMG5_PRId " %" MMG5_PRId ".\n",
                        __func__,MMGS_indPt(mesh,pt->v[MMG5_iprv2[i]]),
                        MMGS_indPt(mesh,pt->v[MMG5_inxt2[i]]));
                MMG5_DEL_MEM(mesh,hash.item);
                return -1;
              }
              MMGS_delPt(mesh,vx[i]);
            }
          }
        }
      }
    }
    nc += ni;
  }
  while( ni > 0 && ++it < 40 );

  if ( mesh->info.ddebug && nc ) {
    fprintf(stdout,"     %" MMG5_PRId " corrected,  %" MMG5_PRId " invalid\n",nc,ni);
    fflush(stdout);
  }

  /* step 4. splitting */
  ns = 0;
  nt = mesh->nt;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )  continue;
    else if ( pt->flag == 0 )  continue;

    j  = -1;
    vx[0] = vx[1] = vx[2] = 0;
    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      if ( MG_GET(pt->flag,i) ) {
        vx[i] = MMG5_hashGet(&hash,pt->v[i1],pt->v[i2]);
        if ( !vx[i] ) {
          if ( !mmgWarn3 ) {
            mmgWarn3 = 1;
            fprintf(stderr,"\n  ## Error: %s: unable to create point on"
                    " at least 1 edge.\n Exit program.\n",__func__);
          }
          return -1;
        }
        j = i;
      }
    }
    if ( pt->flag == 1 || pt->flag == 2 || pt->flag == 4 ) {
      ier = MMGS_split1(mesh,met,k,j,vx);
      ns++;
    }
    else if ( pt->flag == 7 ) {
      ier = MMGS_split3(mesh,met,k,vx);
      ns++;
    }
    else {
      ier = MMGS_split2(mesh,met,k,vx);
      ns++;
    }
    if ( !ier ) return -1;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);
  MMG5_DEL_MEM(mesh,hash.item);

  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of element to split.
 * \param i index of edge to split.
 * \return -1 if lack of memory, 0 if the edge should not be split and 1
 * if success.
 *
 * Check if splitting edge \a i of element \a k is ok.
 *
 */
MMG5_int chkspl(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int i) {
  MMG5_pTria    pt,pt1;
  MMG5_pPoint   ppt;
  MMG5_pxPoint  go;
  MMG5_Bezier   b;
  double        s,uv[2],o[3],no[3],to[3];
  int           ier;
  MMG5_int      ip,*adja,jel;
  int8_t        i1,i2,j,jj,j2;

  if ( mesh->xp > mesh->xpmax-2 )  return 0;
  pt = &mesh->tria[k];
  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];
  if ( MS_SIN(pt->tag[i1]) || MS_SIN(pt->tag[i2]) )  return 0;
  adja = &mesh->adja[3*(k-1)+1];
  jel  = adja[i] / 3;
  if ( jel ) {
    j   = adja[i] % 3;
    jj  = MMG5_inxt2[j];
    j2  = MMG5_iprv2[j];
    pt1 = &mesh->tria[jel];
    if ( MS_SIN(pt1->tag[jj]) || MS_SIN(pt1->tag[j2]) )  return 0;
  }

  ier = MMG5_bezierCP(mesh,pt,&b,1);
  assert ( ier );

  /* create midedge point */
  uv[0] = 0.5;
  uv[1] = 0.5;
  if (i == 1)         uv[0] = 0.0;
  else if ( i == 2 )  uv[1] = 0.0;

  ier = MMGS_bezierInt(&b,uv,o,no,to);
  assert(ier);
  ip = MMGS_newPt(mesh,o,MG_EDG(pt->tag[i]) ? to : no);
  if ( !ip ) {
    /* reallocation of point table */
    MMGS_POINT_REALLOC(mesh,met,ip,mesh->gap,
                        MMG5_INCREASE_MEM_MESSAGE();
                        return -1
                        ,o,MG_EDG(pt->tag[i]) ? to : no);
  }

  if ( MG_EDG(pt->tag[i]) ) {
    ++mesh->xp;
    ppt = &mesh->point[ip];
    ppt->tag = pt->tag[i];
    ppt->xp  = mesh->xp;
    go = &mesh->xpoint[mesh->xp];
    memcpy(go->n1,no,3*sizeof(double));
  }
  s = 0.5;

  ier = intmet(mesh,met,k,i,ip,s);
  if ( !ier ) return 0;

  return ip;
}

/**
 * \param mesh pointer toward mesh structure
 * \param met pointer toward met structure
 * \param typchk type of check to perform: 1 for first stage (adaptation to
 * capture roughly the surface mesh), 2 for second stage of adaptation
 * (rough capture of input metric).
 *
 * Attempt to collapse small edges
 *
 */
static MMG5_int colelt(MMG5_pMesh mesh,MMG5_pSol met,int8_t typchk) {
  MMG5_pTria    pt;
  MMG5_pPoint   p1,p2;
  MMG5_pPar     par;
  double        ll,ux,uy,uz,hmin;
  MMG5_int      nc,list[MMG5_TRIA_LMAX+2],k;
  int           l,isloc,ier,ilist;
  int8_t        i,i1,i2;

  /* Local function pointer toward the suitable functions to use for edge length
   * and quality computation depending on the adaptation phase */
  double (*MMGS_lenEdg)(MMG5_pMesh mesh,MMG5_pSol sol ,MMG5_int ,MMG5_int, int8_t ) = NULL;
  double (*MMGS_caltri)(MMG5_pMesh mesh,MMG5_pSol sol ,MMG5_pTria pt ) = MMG5_caltri_iso;

  MMGS_set_localFunc ( met, typchk, &MMGS_lenEdg, &MMGS_caltri);

  nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

    /* check edge length */
    pt->flag = 0;

    for (i=0; i<3; i++) {
      if ( MS_SIN(pt->tag[i]) )  continue;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      if ( p1->tag & MG_NOM || p2->tag & MG_NOM )  continue;
      else if ( MS_SIN(p1->tag) )   continue;
      else if ( p1->tag > p2->tag || p1->tag > pt->tag[i] )  continue;

      /* check length */
      if ( typchk == 1 ) {
        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];
        uz = p2->c[2] - p1->c[2];
        ll = ux*ux + uy*uy + uz*uz;

        /* local parameters*/
        hmin  = mesh->info.hmin;
        isloc = 0;
        for (l=0; l<mesh->info.npar; l++) {
          par = &mesh->info.par[l];
          if ( /*((par->elt == MMG5_Vertex) &&
                ( (p1->ref == par->ref ) || (p2->ref == par->ref) ))
                || */((par->elt == MMG5_Triangle) && (pt->ref == par->ref) ) ) {
            if ( !isloc ) {
              hmin  = par->hmin;
              isloc = 1;
            }
            else {
              hmin  = MG_MAX(par->hmin,hmin);
            }
          }
        }
        if ( ll > hmin*hmin )  continue;
      }
      else {
        ll = MMG5_lenSurfEdg(mesh,met,pt->v[i1],pt->v[i2],0);
        if ( !ll ) return -1;
        if ( ll > MMGS_LSHRT )  continue;
      }

      /* check if geometry preserved */
      ilist = chkcol(mesh,met,k,i,list,typchk,MMGS_lenEdg,MMGS_caltri);

      int8_t open = (mesh->adja[3*(k-1)+1+i] == 0);

      if ( ilist+open > 3 ) {
        ier = colver(mesh,list,ilist);
        if ( !ier ) return -1;
        nc += ier;
        break;
      }
      else if ( ilist == 3 ) {
        ier = colver3(mesh,list);
        if ( !ier ) return -1;
        nc += ier;
        break;
      }
      else if ( ilist == 2 ) {
        ier = colver2(mesh,list);
        if ( !ier ) return -1;
        nc += ier;
        break;
      }
    }
  }
  if ( nc > 0 && (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) )
    fprintf(stdout,"     %8" MMG5_PRId " vertices removed\n",nc);

  return nc;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * \return -1 if failed or number of new points.
 *
 * Split edges of length bigger than MMGS_LOPTL.
 *
 */
static MMG5_int adpspl(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  MMG5_pPoint   p1,p2;
  double        len,lmax;
  int           ier;
  MMG5_int      ip,ns,k;
  int8_t        i,i1,i2,imax;

  ns = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

    /* check edge length */
    pt->flag = 0;
    imax = -1;
    lmax = -1.0;
    for (i=0; i<3; i++) {
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      len = MMG5_lenSurfEdg(mesh,met,pt->v[i1],pt->v[i2],0);

      if ( !len ) return -1;

      if ( len > lmax ) {
        lmax = len;
        imax = i;
      }
    }
    if ( lmax < MMGS_LOPTL )  continue;
    else if ( MS_SIN(pt->tag[imax]) )  continue;

    /* check length */
    i1 = MMG5_inxt2[imax];
    i2 = MMG5_iprv2[imax];
    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt->v[i2]];
    if ( p1->tag & MG_NOM || p2->tag & MG_NOM )  continue;

    ip = chkspl(mesh,met,k,imax);
    if ( ip < 0 ) {
      /* Lack of memory, go to collapse step. */
      return (ns);
    }
    else if ( ip > 0 ) {
      ier = MMGS_simbulgept(mesh,met,k,imax,ip);
      if ( !ier ) {
        ier = MMGS_dichoto1b(mesh,met,k,imax,ip);
        if ( !ier ) {
          MMGS_delPt(mesh,ip);
          continue;
        }
      }
      ier = split1b(mesh,k,imax,ip);

      if ( !ier ) {
        /* Lack of memory, go to collapse step. */
        MMGS_delPt(mesh,ip);
        return ns;
      }
      /* if we realloc memory in split1b pt pointer is not valid anymore. */
      ns += ier;
    }
  }
  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return -1 if failed.
 * \return number of deleted points.
 *
 * Collapse edges of length smaller than MMGS_LOPTS.
 *
 */
static MMG5_int adpcol(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt;
  MMG5_pPoint   p1,p2;
  double        len;
  MMG5_int      nc,k,list[MMG5_TRIA_LMAX+2];
  int           ier,ilist;
  int8_t        i,i1,i2;

  nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

    /* check edge length */
    pt->flag = 0;
    for (i=0; i<3; i++) {
      if ( MS_SIN(pt->tag[i]) )  continue;

      /* check length */
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      if ( p1->tag & MG_NOM || p2->tag & MG_NOM )  continue;

      len = MMG5_lenSurfEdg(mesh,met,pt->v[i1],pt->v[i2],0);
      if ( !len ) return -1;
      else if ( len > MMGS_LOPTS )  continue;

      p1 = &mesh->point[pt->v[i1]];
      p2 = &mesh->point[pt->v[i2]];
      if ( MS_SIN(p1->tag) )  continue;
      else if ( p1->tag > p2->tag || p1->tag > pt->tag[i] )  continue;

      /* check if geometry preserved */
      ilist = chkcol(mesh,met,k,i,list,2,MMG5_lenSurfEdg,MMG5_calelt);

      int8_t open = (mesh->adja[3*(k-1)+1+i] == 0);

      if ( ilist+open > 3 ) {
        ier =  colver(mesh,list,ilist);;
        nc +=  ier;
        if ( !ier ) return -1;
        break;
      }
      else if ( ilist == 3 ) {
        ier = colver3(mesh,list);
        nc += ier;
        if ( !ier ) return -1;
        break;
      }
      else if ( ilist == 2 ) {
        ier = colver2(mesh,list);
        nc += ier;
        if ( !ier ) return -1;
        break;
      }
    }
  }
  return nc;
}


/* analyze triangles and split or collapse to match gradation */
static int adptri(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int* permNodGlob) {
  int         it,maxit;
  int64_t     nnc,nns,nnf,nnm;
  MMG5_int    nc,ns,nf,nm;

  /* iterative mesh modifications */
  it = nnc = nns = nnf = nnm = 0;
  maxit = 10;
  do {
    if ( !mesh->info.noinsert ) {
      ns = adpspl(mesh,met);
      if ( ns < 0 ) {
        fprintf(stderr,"\n  ## Unable to complete mesh. Exit program.\n");
        return 0;
      }

      /* renumbering if available and needed */
      if ( it==1 && !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
        return 0;

      nc = adpcol(mesh,met);
      if ( nc < 0 ) {
        fprintf(stderr,"\n  ## Unable to complete mesh. Exit program.\n");
        return 0;
      }
    }
    else {
      ns = 0;
      nc = 0;
    }

    if ( !mesh->info.noswap ) {
      nf = swpmsh(mesh,met,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      nm = movtri(mesh,met,1);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
    }
    else  nm = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    nnm += nm;
    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nf+nm > 0 )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",ns,nc,nf,nm);
    if ( ns < 10 && MMG5_abs(nc-ns) < 3 )  break;
    else if ( it > 3 && MMG5_abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;
  }
  while( ++it < maxit && nc+ns > 0 );

  /* renumbering if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

  /*shape optim*/
  it  = 0;
  maxit = 2;
  do {

    if ( !mesh->info.nomove ) {
      nm = movtri(mesh,met,5);
      if ( nm < 0 ) {
        fprintf(stderr,"\n  ## Unable to improve mesh.\n");
        return 0;
      }
      nnm += nm;
    }
    else  nm = 0;

    if ( !mesh->info.noswap ) {
      nf = swpmsh(mesh,met,2);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
      nnf += nf;
    }
    else  nf = 0;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nf+nm > 0 ){
      fprintf(stdout,"                                            ");
      fprintf(stdout,"%8" MMG5_PRId " swapped, %8" MMG5_PRId " moved\n",nf,nm);
    }
  }
  while( ++it < maxit && nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = movtri(mesh,met,5);
    if ( nm < 0 ) {
      fprintf(stderr,"\n  ## Unable to improve mesh.\n");
      return 0;
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ) {
    fprintf(stdout,"                                                              %8" MMG5_PRId " moved\n",nm);
  }

  if ( mesh->info.imprim > 0 ) {
    if ( abs(mesh->info.imprim) < 5 && (nnc > 0 || nns > 0) )
      fprintf(stdout,"     %8" PRId64 " splitted, %8" PRId64 " collapsed, %8" PRId64 " swapped, %8" PRId64 " moved, %d iter. \n",nns,nnc,nnf,nnm,it);
  }
  return 1;
}

/* analyze tetrahedra and split if needed */
static int anatri(MMG5_pMesh mesh,MMG5_pSol met,int8_t typchk) {
  MMG5_int     nc,ns,nf,nnc,nns,nnf;
  int          it,maxit;

  /* analyze tetras : initial splitting */
  nns = nnc = nnf = it = 0;
  maxit = 5;
  do {
    if ( !mesh->info.noinsert ) {
      /* memory free */
      MMG5_DEL_MEM(mesh,mesh->adja);
      mesh->adja = 0;

      /* analyze surface */
      ns = anaelt(mesh,met,typchk);
      if ( ns < 0 ) {
        fprintf(stderr,"\n  ## Unable to complete surface mesh. Exit program.\n");
        return 0;
      }

      if ( !MMGS_hashTria(mesh) ) {
        fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
        return 0;
      }

      /* collapse short edges */
      nc = colelt(mesh,met,typchk);
      if ( nc < 0 ) {
        fprintf(stderr,"\n  ## Unable to collapse mesh. Exiting.\n");
        return 0;
      }
    }
    else {
      ns = 0;
      nc = 0;
    }

    /* attempt to swap */
    if ( !mesh->info.noswap ) {
      nf = swpmsh(mesh,met,typchk);
      if ( nf < 0 ) {
        fprintf(stderr,"\n  ## Unable to improve mesh. Exiting.\n");
        return 0;
      }
    }
    else nf = 0;

    nnc += nc;
    nns += ns;
    nnf += nf;
    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc > 0 )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped\n",ns,nc,nf);
    if ( it > 3 && MMG5_abs(nc-ns) < 0.1 * MG_MAX(nc,ns) )  break;
  }
  while ( ++it < maxit && ns+nc+nf > 0 );

  if ( mesh->info.imprim > 0 ) {
    if ( (abs(mesh->info.imprim) < 5 || mesh->info.ddebug ) && nns+nnc > 0 )
      fprintf(stdout,"     %8" MMG5_PRId " splitted, %8" MMG5_PRId " collapsed, %8" MMG5_PRId " swapped, %d iter.\n",nns,nnc,nnf,it);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param permNodGlob if provided, strore the global permutation of nodes.
 * \return 0 if failed, 1 if success.
 *
 * Main adaptation routine.
 *
 */
int MMG5_mmgs1(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int *permNodGlob) {

  /* renumbering if available */
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  /*--- stage 1: geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !anatri(mesh,met,1) ) {
    fprintf(stderr,"\n  ## Unable to split mesh-> Exiting.\n");
    return 0;
  }

  /* Debug: export variable MMG_SAVE_ANATRI1 to save adapted mesh at the end of
   * anatri wave */
  if ( getenv("MMG_SAVE_ANATRI1") ) {
    printf("  ## WARNING: EXIT AFTER ANATRI-1."
           " (MMG_SAVE_ANATRI1 env variable is exported).\n");
    return 1;
  }

  /* renumbering if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

  /*--- stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  /* define metric map */
  if ( !MMGS_defsiz(mesh,met) ) {
    fprintf(stderr,"\n  ## Metric undefined. Exit program.\n");
    return 0;
  }

 /* Debug: export variable MMG_SAVE_DEFSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_DEFSIZ") ) {
    printf("  ## WARNING: EXIT AFTER DEFSIZ."
           " (MMG_SAVE_DEFSIZ env variable is exported).\n");
    return 1;
  }

  MMG5_gradation_info(mesh);
  if ( mesh->info.hgrad > 0. ) {
     if (!MMGS_gradsiz(mesh,met) ) {
      fprintf(stderr,"\n  ## Gradation problem. Exit program.\n");
      return 0;
    }
  }

  if ( mesh->info.hgradreq > 0. ) {
    MMGS_gradsizreq(mesh,met);
  }

  /* Debug: export variable MMG_SAVE_GRADSIZ to save adapted mesh at the end of
   * anatet wave */
  if ( getenv("MMG_SAVE_GRADSIZ") ) {
    printf("  ## WARNING: EXIT AFTER GRADSIZ."
           " (MMG_SAVE_GRADSIZ env variable is exported).\n");
    return 1;
  }

  if ( !anatri(mesh,met,2) ) {
    fprintf(stderr,"\n  ## Unable to proceed adaptation. Exit program.\n");
    return 0;
  }
  /* Debug: export variable MMG_SAVE_ANATRI1 to save adapted mesh at the end of
   * anatri wave */
  if ( getenv("MMG_SAVE_ANATRI2") ) {
    printf("  ## WARNING: EXIT AFTER ANATRI-2."
           " (MMG_SAVE_ANATRI2 env variable is exported).\n");
    return 1;
  }

  /* renumbering if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,permNodGlob) )
    return 0;

  if ( !adptri(mesh,met,permNodGlob) ) {
    fprintf(stderr,"\n  ## Unable to adapt. Exit program.\n");
    return 0;
  }

  return 1;
}
