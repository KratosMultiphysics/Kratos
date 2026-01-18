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
 * \file mmgs/colver_s.c
 * \brief Functions for vertices collapsing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmgs_private.h"
#include "mmgexterns_private.h"
#include "inlined_functions_private.h"

/**
 * \param mesh pointer to the mesh
 * \param met pointer to the metric
 * \param k index of the element in wich we collapse
 * \param i index of the edge to collapse
 * \param list pointer to the ball of point
 * \param typchk type of check to perform
 * \param MMGS_lenEdg pointer to the suitable fct to compute edge lengths
 * depending on presence of input metric, metric type (iso/aniso) and \a typchk
 * value (i.e. stage of adaptation)
 * \param MMGS_caltri pointer to the suitable fct to compute tria quality
 * depending on presence of input metric, metric type (iso/aniso) and \a typchk
 * value (i.e. stage of adaptation)
 *
 * \return 0 if we can't move of if we fail, 1 if success
 *
 * check if geometry preserved by collapsing edge i
 *
 */
int chkcol(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i,MMG5_int *list,int8_t typchk,
           double (*MMGS_lenEdg)(MMG5_pMesh,MMG5_pSol,MMG5_int ,MMG5_int,int8_t),
           double (*MMGS_caltri)(MMG5_pMesh,MMG5_pSol,MMG5_pTria)) {
  MMG5_pTria     pt,pt0,pt1,pt2;
  MMG5_pPoint    p1,p2;
  double         len,lon,ps,cosnold,cosnnew,kal,n0old[3],n1old[3],n00old[3];
  double         n0new[3],n1new[3],n00new[3];
  MMG5_int       *adja,jel,kel,ip1,ip2,l,ll;
  int            ilist;
  int8_t         i1,i2,j,jj,j2,lj,open,voy;

  pt0 = &mesh->tria[0];
  pt  = &mesh->tria[k];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

#ifndef NDEBUG
  /* suppression of maybe-uninitialized value on arm */
  lon = 0.;
  n00old[0] = n00old[1] = n00old[2] = 0.;
  n0old[0]  = n0old[1]  = n0old[2]  = 0.;
  n1old[0]  = n1old[1]  = n1old[2]  = 0.;
  n00new[0] = n00new[1] = n00new[2] = 0.;
  n0new[0]  = n0new[1]  = n0new[2]  = 0.;
  n1new[0]  = n1new[1]  = n1new[2]  = 0.;
#endif

  if ( MMGS_lenEdg ) {
    lon = MMGS_lenEdg(mesh,met,ip1,ip2,0);
    if ( !lon ) return 0;
    lon = MG_MIN(lon,MMGS_LSHRT);
    lon = MG_MAX(1.0/lon,MMGS_LLONG);
  }

  /* collect all triangles around vertex i1 */
  ilist = boulechknm(mesh,k,i1,list);
  if ( ilist <= 0 )  return 0;

  /* check for open ball */
  adja = &mesh->adja[3*(k-1)+1];
  open = adja[i] == 0;

  if ( ilist+open > 3 ) {
    /* check references */
    if ( MG_EDG(pt->tag[i2]) ) {
      jel = list[1] / 3;
      pt1 = &mesh->tria[jel];
      if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;
    }

    /* analyze ball */
    assert ( ilist-1+open > 1 );
    for (l=1; l<ilist-1+open; l++) {
      jel = list[l] / 3;
      j   = list[l] % 3;
      jj  = MMG5_inxt2[j];
      j2  = MMG5_iprv2[j];
      pt1 = &mesh->tria[jel];

      /* check length */
      if ( MMGS_lenEdg ) {
        ip1 = pt1->v[j2];
        len = MMGS_lenEdg(mesh,met,ip1,ip2,0);
        if ( len > lon || !len )  return 0;
      }

      /* check normal flipping */
      if ( !MMG5_nortri(mesh,pt1,n1old) )  return 0;
      memcpy(pt0,pt1,sizeof(MMG5_Tria));
      pt0->v[j] = ip2;
      if ( !MMG5_nortri(mesh,pt0,n1new) )  return 0;

      ps = n1new[0]*n1old[0] + n1new[1]*n1old[1]  + n1new[2]*n1old[2];
      if ( ps < 0.0 )  return 0;

      /* keep normals at 1st triangles */
      if ( l == 1 && !open ) {
        memcpy(n00old,n1old,3*sizeof(double));
        memcpy(n00new,n1new,3*sizeof(double));
      }

      /* check normals deviation */
      if ( !(pt1->tag[j2] & MG_GEO) ) {
        if ( l > 1 ) {
          cosnold = n0old[0]*n1old[0] + n0old[1]*n1old[1] + n0old[2]*n1old[2];
          cosnnew = n0new[0]*n1new[0] + n0new[1]*n1new[1] + n0new[2]*n1new[2];
          if ( cosnold < MMG5_ANGEDG ) {
            if ( cosnnew < cosnold )  return 0;
          }
          else if ( cosnnew < MMG5_ANGEDG )  return 0;
        }
      }

      /* check geometric support */
      if ( l == 1 ) {
        pt0->tag[j2] |= pt->tag[i1];
      }
      else if ( l == ilist-2+open ) {
        if ( !open ) {
          ll = list[ilist-1] / 3;
          lj = list[ilist-1] % 3;
          pt0->tag[jj] |= mesh->tria[ll].tag[lj];
        }
        else {
          assert ( list[0]/3 == k );
          pt0->tag[jj] |= pt->tag[i];
        }
      }
      if ( chkedg(mesh,0) )  return 0;

      /* check quality */
      kal = MMGS_ALPHAD*MMGS_caltri(mesh,met,pt0);

      if ( kal < MMGS_NULKAL )  return 0;

      memcpy(n0old,n1old,3*sizeof(double));
      memcpy(n0new,n1new,3*sizeof(double));
    }

    /* check angle between 1st and last triangles */
    if ( !open && !(pt->tag[i] & MG_GEO) ) {
      cosnold = n00old[0]*n1old[0] + n00old[1]*n1old[1] + n00old[2]*n1old[2];
      cosnnew = n00new[0]*n1new[0] + n00new[1]*n1new[1] + n00new[2]*n1new[2];
      if ( cosnold < MMG5_ANGEDG ) {
        if ( cosnnew < cosnold )  return 0;
      }
      else if ( cosnnew < MMG5_ANGEDG )  return 0;

      /* other checks for reference collapse */
      jel = list[ilist-1] / 3;
      j   = list[ilist-1] % 3;
      j   = MMG5_iprv2[j];
      pt  = &mesh->tria[jel];
      if ( MG_EDG(pt->tag[j]) ) {
        jel = list[ilist-2] / 3;
        pt1 = &mesh->tria[jel];
        if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;
      }
    }
  }

  /* specific test: no collapse if any interior edge is EDG */
  else if ( ilist == 3 ) {
    /* Remark: if ilist==3 and i is an open ridge, we pass in the previous
     * test (open+ilist > 3) so here, ip1 is in the middle of the 3
     * triangles */

    p1 = &mesh->point[pt->v[i1]];
    if ( MS_SIN(p1->tag) )  return 0;
    else if ( MG_EDG(pt->tag[i2]) && !MG_EDG(pt->tag[i]) )  return 0;
    else if ( !MG_EDG(pt->tag[i2]) && MG_EDG(pt->tag[i]) )  return 0;
    else if ( MG_EDG(pt->tag[i2]) && MG_EDG(pt->tag[i]) && MG_EDG(pt->tag[i1]) )  return 0;

    /* Check geometric approximation */
    jel = list[1] / 3;
    j   = list[1] % 3;
    jj  = MMG5_inxt2[j];
    j2  = MMG5_iprv2[j];
    pt0 = &mesh->tria[0];
    pt1 = &mesh->tria[jel];
    memcpy(pt0,pt1,sizeof(MMG5_Tria));
    pt0->v[j] = ip2;

    jel = list[2] / 3;
    j   = list[2] % 3;
    pt1 = &mesh->tria[jel];
    pt0->tag[jj] |= pt1->tag[j];
    pt0->tag[j2] |= pt1->tag[MMG5_inxt2[j]];
    if ( chkedg(mesh,0) )  return 0;

    /* check quality */
    kal = MMGS_ALPHAD*MMGS_caltri(mesh,met,pt0);

    if ( kal < MMGS_NULKAL )  return 0;
  }

  /* for specific configurations along open ridge */
  else if ( ilist == 2 ) {
    if ( !open )  return 0;

    jel = list[1] / 3;
    j   = list[1] % 3;

    /* Topological test */
    adja = &mesh->adja[3*(jel-1)+1];
    kel = adja[j] / 3;
    voy = adja[j] % 3;
    pt2 = &mesh->tria[kel];

    if ( pt2->v[voy] == ip2) return 0;

    jj  = MMG5_inxt2[j];
    pt1 = &mesh->tria[jel];
    if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;
    else if ( !(pt1->tag[jj] & MG_GEO) )  return 0;

    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt1->v[j]];
    if ( p2->tag > p1->tag || p2->ref != p1->ref )  return 0;

    /* Check geometric approximation */
    j2  = MMG5_iprv2[j];
    pt0 = &mesh->tria[0];
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i1] = pt1->v[j2];

    if ( chkedg(mesh,0) )  return 0;

    /* check quality */
    kal = MMGS_ALPHAD*MMGS_caltri(mesh,met,pt0);

    if ( kal < MMGS_NULKAL )  return 0;

  }

  return ilist;
}

/* collapse edge i of k, i1->i2 */
int colver(MMG5_pMesh mesh,MMG5_int *list,int ilist) {
  MMG5_pTria    pt,pt1,pt2;
  MMG5_int      *adja,k,iel,jel,kel,ip1,ip2;
  int8_t        i,i1,i2,j,jj,open;

  iel = list[0] / 3;
  i1  = list[0] % 3;
  i   = MMG5_iprv2[i1];
  i2  = MMG5_inxt2[i1];
  pt  = &mesh->tria[iel];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  /* check for open ball */
  adja = &mesh->adja[3*(iel-1)+1];
  open = adja[i] == 0;

  /* update vertex ip1 -> ip2 */
  for (k=1; k<ilist-1+open; k++) {
    jel = list[k] / 3;
    jj  = list[k] % 3;
    pt1 = &mesh->tria[jel];
    pt1->v[jj] = ip2;
    pt1->base  = mesh->base;
  }

  /* update adjacent with 1st elt */
  jel = list[1] / 3;
  jj  = list[1] % 3;
  j   = MMG5_iprv2[jj];
  pt1 = &mesh->tria[jel];
  pt1->tag[j] |= pt->tag[i1];
  pt1->edg[j] = MG_MAX(pt->edg[i1],pt1->edg[j]);
  if ( adja[i1] ) {
    kel = adja[i1] / 3;
    k   = adja[i1] % 3;
    mesh->adja[3*(kel-1)+1+k] = 3*jel + j;
    mesh->adja[3*(jel-1)+1+j] = 3*kel + k;
    pt2 = &mesh->tria[kel];
    pt2->tag[k] |= pt1->tag[j];
    pt2->edg[k] = MG_MAX(pt1->edg[j],pt2->edg[k]);
  }
  else
    mesh->adja[3*(jel-1)+1+j] = 0;

  /* adjacent with last elt */
  if ( !open ) {
    iel = list[ilist-1] / 3;
    i1  = list[ilist-1] % 3;
    pt  = &mesh->tria[iel];

    jel = list[ilist-2] / 3;
    jj  = list[ilist-2] % 3;
    j   = MMG5_inxt2[jj];
    pt1 = &mesh->tria[jel];
    pt1->tag[j] |= pt->tag[i1];
    pt1->edg[j] = MG_MAX(pt->edg[i1],pt1->edg[j]);
    adja = &mesh->adja[3*(iel-1)+1];
    if ( adja[i1] ) {
      kel = adja[i1] / 3;
      k   = adja[i1] % 3;
      mesh->adja[3*(kel-1)+1+k] = 3*jel + j;
      mesh->adja[3*(jel-1)+1+j] = 3*kel + k;
      pt2 = &mesh->tria[kel];
      pt2->tag[k] |= pt1->tag[j];
      pt2->edg[k] = MG_MAX(pt1->edg[j],pt2->edg[k]);
    }
    else
      mesh->adja[3*(jel-1)+1+j] = 0;
  }

  MMGS_delPt(mesh,ip1);
  if ( !MMGS_delElt(mesh,list[0] / 3) ) return 0;
  if ( !open ) {
    if ( !MMGS_delElt(mesh,list[ilist-1] / 3) )  return 0;
  }

  return 1;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param list pointer to the ball of the point to collapse.
 * \return 1 if success, 0 if fail.
 *
 * Collapse edge \f$list[0]\%3\f$ in tet \f$list[0]/3\f$ (\f$ ip->i1\f$ ) for a
 * ball of the collapsed point of size 3: the collapsed point is removed.
 *
 */
int colver3(MMG5_pMesh mesh,MMG5_int* list) {
  MMG5_pTria   pt,pt1,pt2;
  MMG5_int     *adja,iel,jel,kel,mel,ip;
  int8_t       i,i1,j,j1,j2,k,m;

  /* update of new point for triangle list[0] */
  iel = list[0] / 3;
  i   = list[0] % 3;
  i1  = MMG5_inxt2[i];
  pt  = &mesh->tria[iel];
  ip  = pt->v[i];

  jel = list[1] / 3;
  j   = list[1] % 3;
  j1  = MMG5_inxt2[j];
  j2  = MMG5_iprv2[j];
  pt1 = &mesh->tria[jel];

  kel = list[2] / 3;
  k   = list[2] % 3;
  pt2 = &mesh->tria[kel];

  /* update info */
  pt1->v[j]     = pt->v[i1];
  pt1->tag[j1] |= pt2->tag[k];
  pt1->edg[j1]  = MG_MAX(pt1->edg[j1],pt2->edg[k]);
  pt1->tag[j2] |= pt->tag[i];
  pt1->edg[j2]  = MG_MAX(pt1->edg[j2],pt->edg[i]);
  pt1->base     = mesh->base;

  /* update neighbours of new triangle */
  adja = &mesh->adja[3*(jel-1)+1];
  adja[j1] = mesh->adja[3*(kel-1)+1+k];
  adja[j2] = mesh->adja[3*(iel-1)+1+i];

  mel  = adja[j2] / 3;
  if ( mel ) {
    m    = adja[j2] % 3;
    pt   = &mesh->tria[mel];
    pt->tag[m]  = pt1->tag[j2];
    pt->edg[m]  = pt1->edg[j2];
    mesh->adja[3*(mel-1)+1+m] = 3*jel + j2;
  }

  mel = adja[j1] / 3;
  if ( mel ) {
    m    = adja[j1] % 3;
    pt   = &mesh->tria[mel];
    pt->tag[m]  = pt1->tag[j1];
    pt->edg[m]  = pt1->edg[j1];
    mesh->adja[3*(mel-1)+1+m] = 3*jel + j1;
  }

  /* remove vertex + elements */
  MMGS_delPt(mesh,ip);
  if ( !MMGS_delElt(mesh,iel) ) return 0;
  if ( !MMGS_delElt(mesh,kel) ) return 0;

  return 1;
}


/* collapse point along open ridge */
int colver2(MMG5_pMesh mesh,MMG5_int* list) {
  MMG5_pTria   pt,pt1;
  MMG5_int     *adja,iel,jel,kel,ip;
  int8_t       i1,i2,jj,j2,k;

  /* update of new point for triangle list[0] */
  iel = list[0] / 3;
  i1  = list[0] % 3;
  i2  = MMG5_inxt2[i1];
  pt  = &mesh->tria[iel];
  ip  = pt->v[i1];

  jel = list[1] / 3;
  j2  = list[1] % 3;
  jj  = MMG5_iprv2[j2];
  pt1 = &mesh->tria[jel];

  /* update info */
  pt->v[i1] = pt1->v[jj];
  pt->tag[i2] |= pt1->tag[j2];
  pt->edg[i2] = pt1->edg[j2];
  pt->base  = mesh->base;

  /* update neighbours of new triangle */
  adja = &mesh->adja[3*(iel-1)+1];
  adja[i2] = mesh->adja[3*(jel-1)+1+j2];
  adja = &mesh->adja[3*(jel-1)+1];
  kel  = adja[j2] / 3;
  k    = adja[j2] % 3;
  if ( kel )
    mesh->adja[3*(kel-1)+1+k] = 3*iel + i2;

  /* remove vertex + element */
  MMGS_delPt(mesh,ip);
  if ( !MMGS_delElt(mesh,jel) ) return 0;

  return 1;
}

/* collapse edge i of k, i1->i2 */
int litcol(MMG5_pMesh mesh,MMG5_int k,int8_t i,double kali) {
  MMG5_pTria     pt,pt0,pt1;
  MMG5_pPoint    p1,p2;
  double         kal,ps,cosnold,cosnnew;
  double         n0old[3],n0new[3],n1old[3],n1new[3],n00old[3],n00new[3];
  MMG5_int       list[MMG5_TRIA_LMAX+2],jel,ip2,l;
  int            ilist;
  int8_t         i1,i2,j,jj,j2,open;

  pt0 = &mesh->tria[0];
  pt  = &mesh->tria[k];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip2 = pt->v[i2];

#ifndef NDEBUG
  n00old[0] = n00old[1] = n00old[2] = 0.;
  n0old[0]  = n0old[1]  = n0old[2]  = 0.;
  n1old[0]  = n1old[1]  = n1old[2]  = 0.;
  n00new[0] = n00new[1] = n00new[2] = 0.;
  n0new[0]  = n0new[1]  = n0new[2]  = 0.;
  n1new[0]  = n1new[1]  = n1new[2]  = 0.;
#endif

  /* collect all triangles around vertex i1 */
  if ( pt->v[i1] & MG_NOM )  return 0;

  ilist = MMG5_boulet(mesh,k,i1,list,1,&open);

#ifndef NDEBUG
  /* check for open ball */
  int8_t opn;
  MMG5_int *adja = &mesh->adja[3*(k-1)+1];
  opn = adja[i] == 0;
  assert ( opn == open );
#endif

  if ( ilist > 3 ) {
    /* check references */
    jel = list[1] / 3;
    pt1 = &mesh->tria[jel];
    if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;

    /* analyze ball */
    assert ( ilist-1+open > 1 );
    for (l=1; l<ilist-1+open; l++) {
      jel = list[l] / 3;
      j   = list[l] % 3;
      j2  = MMG5_iprv2[j];
      pt1 = &mesh->tria[jel];

      /* check normal flipping */
      if ( !MMG5_nortri(mesh,pt1,n1old) )  return 0;
      memcpy(pt0,pt1,sizeof(MMG5_Tria));
      pt0->v[j] = ip2;
      if ( !MMG5_nortri(mesh,pt0,n1new) )  return 0;
      ps = n1new[0]*n1old[0] + n1new[1]*n1old[1]  + n1new[2]*n1old[2];
      if ( ps < 0.0 )  return 0;

      /* keep normals at 1st triangles */
      if ( l == 1 && !open ) {
        memcpy(n00old,n1old,3*sizeof(double));
        memcpy(n00new,n1new,3*sizeof(double));
      }

      /* check normals deviation */
      if ( !(pt1->tag[j2] & MG_GEO) ) {
        if ( l > 1 ) {
          cosnold = n0old[0]*n1old[0] + n0old[1]*n1old[1] + n0old[2]*n1old[2];
          cosnnew = n0new[0]*n1new[0] + n0new[1]*n1new[1] + n0new[2]*n1new[2];
          if ( cosnold < MMG5_ANGEDG ) {
            if ( cosnnew < MG_MIN(0.0,cosnold) )  return 0;
          }
          else if ( cosnnew < MMG5_ANGEDG )  return 0;
        }

        memcpy(n0old,n1old,3*sizeof(double));
        memcpy(n0new,n1new,3*sizeof(double));
      }
      /* check quality */
      kal = MMGS_ALPHAD*MMG5_caltri_iso(mesh,NULL,pt0);
      if ( kal < MMGS_NULKAL )  return 0;
    }

    /* check angle between 1st and last triangles */
    if ( !open ) {
      cosnold = n00old[0]*n1old[0] + n00old[1]*n1old[1] + n00old[2]*n1old[2];
      cosnnew = n00new[0]*n1new[0] + n00new[1]*n1new[1] + n00new[2]*n1new[2];
      if ( cosnold < MMG5_ANGEDG ) {
        if ( cosnnew < MG_MIN(0.0,cosnold) )  return 0;
      }
      else if ( cosnnew < MMG5_ANGEDG )  return 0;

      /* other reference checks */
      jel = list[ilist-1] / 3;
      pt  = &mesh->tria[jel];
      jel = list[ilist-2] / 3;
      pt1 = &mesh->tria[jel];
      if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;
    }

    return colver(mesh,list,ilist);
  }

  /* specific test: no collapse if any interior edge is EDG */
  else if ( ilist == 3 ) {
    p1 = &mesh->point[pt->v[i1]];
    if ( MS_SIN(p1->tag) )  return 0;
    else if (  MG_EDG(pt->tag[i2]) && !MG_EDG(pt->tag[i]) )  return 0;
    else if ( !MG_EDG(pt->tag[i2]) &&  MG_EDG(pt->tag[i]) )  return 0;
    else if (  MG_EDG(pt->tag[i2]) &&  MG_EDG(pt->tag[i]) && MG_EDG(pt->tag[i1]) )  return 0;

    return colver3(mesh,list);
  }

  /* for specific configurations along open ridge */
  else if ( ilist == 2 ) {
    if ( !open )  return 0;
    jel = list[1] / 3;
    j   = list[1] % 3;
    jj  = MMG5_inxt2[j];
    pt1 = &mesh->tria[jel];
    if ( MMG5_abs(pt->ref) != MMG5_abs(pt1->ref) )  return 0;
    else if ( !(pt1->tag[jj] & MG_GEO) )  return 0;

    p1 = &mesh->point[pt->v[i1]];
    p2 = &mesh->point[pt1->v[jj]];
    if ( p2->tag > p1->tag || p2->ref != p1->ref )  return 0;

    return colver2(mesh,list);
  }

  return 0;
}



