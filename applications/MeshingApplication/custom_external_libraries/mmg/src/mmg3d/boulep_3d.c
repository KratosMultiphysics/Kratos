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
 * \file mmg3d/boulep_3d.c
 * \brief Functions for ball of points computation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"

extern MMG5_Info  info;

/**
 * \brief Given a vertex and a tetrahedron, find all tetrahedra in the ball of
 * this vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param start index of the starting tetrahedra.
 * \param ip local index of the point in the tetrahedra \a start.
 * \param list pointer to the list of the tetra in the volumic ball of
 * \a ip.
 * \return 0 if fail and the number of the tetra in the ball otherwise.
 *
 * Fill the volumic ball (i.e. filled with tetrahedra) of point \a ip in tetra
 * \a start. Results are stored in the form \f$4*kel + jel\f$, kel = number
 * of tetrahedra, jel = local index of p within kel.
 *
 */
int MMG5_boulevolp (MMG5_pMesh mesh, MMG5_int start, int ip, int64_t * list){
  MMG5_pTetra  pt,pt1;
  MMG5_int     base,*adja,nump,k,k1;
  int          ilist,cur;
  int8_t       j,l,i;

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  assert( 0<=ip && ip<4 && "unexpected local index for vertex");
  nump = pt->v[ip];

  /* Store initial tetrahedron */
  pt->flag = base;
  list[0] = 4*start + ip;
  ilist=1;

  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = MMG5_inxt3[i];
      k1 = adja[i];
      if ( !k1 )  continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }
  return ilist;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param pt pointer to the working tetra
 * \param k index of the tetra \a pt.
 * \param na index of the first extermity of the seeking edge.
 * \param nb index of the second extermity of the seeking edge.
 * \param error 1 if we want to print an error message, 0 for a warning.
 * \param mmgWarn static variable to print warning only once (not used if error==1)
 * \param ia pointer to the edge index (to fill).
 *
 * \return 0 if fail, 1 if success.
 *
 * Find the local index of the edge \a ia in the tetra \a pt of index \a k;
 *
 */
int MMG3D_findEdge(MMG5_pMesh mesh,MMG5_pTetra pt,MMG5_int k,MMG5_int na,MMG5_int nb,int error,
                   int8_t *mmgWarn,int8_t *ia) {
  int8_t ipa,ipb;

  /* identification of edge number in tetra k */
  for ((*ia)=0; (*ia)<6; (*ia)++) {
    ipa = MMG5_iare[(*ia)][0];
    ipb = MMG5_iare[(*ia)][1];
    if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
         (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
  }

  /* fail if the edge na-nb is not found */
  if ( (*ia)<6 ) return 1;

  if ( error ) {
    fprintf(stderr,"\n  ## Error: %s: wrong edge's shell: "
            " edge %" MMG5_PRId " %" MMG5_PRId " not found in tetra %" MMG5_PRId ".\n",__func__,
            MMG3D_indPt(mesh,na),
            MMG3D_indPt(mesh,nb),MMG3D_indElt(mesh,k));
    fprintf(stderr,"  Exit program.\n");
  }
  else {
    if ( !(*mmgWarn) ) {
      (*mmgWarn) = 1;
      fprintf(stderr,"\n  ## Warning: %s: at least one wrong edge's"
              " shell.\n",__func__);
    }
  }
  return 0;
}

static inline
void MMG3D_compute_tangent(MMG5_pMesh mesh,int nump,int ip0,int ip1,double t[3]) {
  MMG5_pPoint ppt,p0,p1;
  double      l0,l1,dd;
  int8_t      i;

  ppt = &mesh->point[nump];
  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];

  l0 = (ppt->c[0] - p0->c[0])*(ppt->c[0] - p0->c[0]) \
    + (ppt->c[1] - p0->c[1])*(ppt->c[1] - p0->c[1]) + (ppt->c[2] - p0->c[2])*(ppt->c[2] - p0->c[2]);
  l1 = (ppt->c[0] - p1->c[0])*(ppt->c[0] - p1->c[0]) \
    + (ppt->c[1] - p1->c[1])*(ppt->c[1] - p1->c[1]) + (ppt->c[2] - p1->c[2])*(ppt->c[2] - p1->c[2]);
  l0 = sqrt(l0);
  l1 = sqrt(l1);

  if ( (l0 < MMG5_EPSD2) || (l1 < MMG5_EPSD2) ) {
    for ( i=0; i<3; ++i ) {
      t[i] = p1->c[i] - p0->c[i];
    }
  }
  else if ( l0 < l1 ) {
    dd = l0 / l1;
    for ( i=0; i<3; ++i ) {
      t[i] = dd*(p1->c[i] - ppt->c[i]) + ppt->c[i] - p0->c[i];
    }
  }
  else {
    dd = l1 / l0;
    for ( i=0; i<3; ++i ) {
      t[i] = dd*(p0->c[i] - ppt->c[i]) + ppt->c[i] - p1->c[i];
    }
  }

  return;
}

/**
 * \param mesh pointer to the mesh  structure.
 * \param start tetra index.
 * \param ip point index.
 * \param iface face index.
 * \param n computed normal vector.
 * \param t computed tangent vector.
 * \return 0 if point is singular, 1 otherwise.
 *
 * Define normal and tangent vectors at a non manifold point (\a ip in \a start,
 * supported by face \a iface), enumerating its (outer)surfacic ball.
 *
 */
int MMG5_boulenm(MMG5_pMesh mesh,MMG5_int start,int ip,int iface,
                  double n[3],double t[3]) {
  MMG5_pTetra   pt;
  double        dd,nt[3];
  int           nr,nnm;
  MMG5_int      base,nump,k,*adja,piv,nvstart,aux,na,nb,adj,fstart,ip0,ip1;
  uint16_t      tag;
  int8_t        iopp,ipiv,indb,inda,i,isface;
  int8_t        indedg[4][4] = { {-1,0,1,2}, {0,-1,3,4}, {1,3,-1,5}, {2,4,5,-1} };

  base = ++mesh->base;
  nr  = nnm = 0;
  ip0 = ip1 = 0;

  memset(n,0x00,3*sizeof(double));
  memset(t,0x00,3*sizeof(double));

  pt   = &mesh->tetra[start];
  nump = pt->v[ip];
  k    = start;

  na   = pt->v[ip];
  nb   = pt->v[MMG5_idir[iface][MMG5_inxt2[MMG5_idirinv[iface][ip]]]];
  piv  = pt->v[MMG5_idir[iface][MMG5_iprv2[MMG5_idirinv[iface][ip]]]];

  iopp   = iface;
  fstart = 4*k+iopp;
  do {
    /* computation of normal and tangent at nump */
    if ( MMG5_norface(mesh,k,iopp,nt) ) {
      n[0] += nt[0];
      n[1] += nt[1];
      n[2] += nt[2];
    }

    if ( pt->xt ) {
      for ( inda=0; inda<4; inda++ ){
        if ( pt->v[inda]==na ) break;
      }
      for ( indb=0; indb<4; indb++ ){
        if ( pt->v[indb]==nb ) break;
      }
      assert( (inda < 4) && (indb < 4));
      tag = mesh->xtetra[pt->xt].tag[indedg[inda][indb]];
    }

    else  tag = 0;

    if ( MG_EDG(tag) && !(tag & MG_NOM) )
      nr++;
    else if ( tag & MG_NOM ) {
      nnm++;
      if ( !ip0 )
        ip0 = nb;
      else
        ip1 = nb;
    }

    /* A boundary face has been hit : change travel edge */
    aux     = nb;
    nb      = piv;
    piv     = aux;
    nvstart = k;
    adj     = k;

    /* Now unfold shell of edge (na,nb) starting from k (included) */
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      if ( !MMG3D_findEdge(mesh,pt,k,na,nb,1,NULL,&i) ) return -1;

      /* set sense of travel */
      if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
        adj = adja[ MMG5_ifar[i][0] ] / 4;
        ipiv = MMG5_ifar[i][1];
        iopp = MMG5_ifar[i][0];
        piv = pt->v[ipiv];
      }
      else {
        adj = adja[ MMG5_ifar[i][1] ] / 4;
        ipiv = MMG5_ifar[i][0];
        iopp = MMG5_ifar[i][1];
        piv = pt->v[ipiv];
      }
      isface = (adja[iopp] == 0);
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  if ( (nr > 0 && nnm > 0) || nnm != 2 ) {
    /* We pass here for non-manifold meshes with only edge connection */
    return 0;
  }

  dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( dd > MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    n[0] *= dd;
    n[1] *= dd;
    n[2] *= dd;
  }
  assert( ip0 && ip1 );
  if ( ip0 == ip1 )  return 0;

  MMG3D_compute_tangent(mesh,nump,ip0,ip1,t);

  dd = t[0]*n[0] + t[1]*n[1] + t[2]*n[2];
  t[0] -= dd*n[0];
  t[1] -= dd*n[1];
  t[2] -= dd*n[2];

  dd = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
  if ( dd > MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    t[0] *= dd;
    t[1] *= dd;
    t[2] *= dd;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh  structure.
 * \param start tetra index.
 * \param ip point index.
 * \param t computed tangent vector.
 * \return 0 when more than two NOM points are attached to ip, 1 if sucess.
 *
 * Travel the ball of the internal non manifold point ip in tetra start
 * and calculate the tangent vector to the underlying curve.
 *
 * \remark we are not able to compute tangent along non-manifold points for
 * edge-connected meshes. In this case the point doesn't have xpoint nor tangent
 * or normal.
 */
int MMG5_boulenmInt(MMG5_pMesh mesh,MMG5_int start,int ip,double t[3]) {
  MMG5_pTetra    pt,pt1;
  MMG5_pxTetra   pxt;
  double         dd;
  MMG5_int       base,k,kk,ip0,ip1,nump,na,nb,list[MMG3D_LMAX+2],*adja;
  int            cur,ilist;
  int8_t         i,j,ii,ie;

  base = ++mesh->base;
  ip0 = ip1 = 0;
  cur = ilist = 0;

  /* Store initial tetrahedron */
  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  list[0] = 4*start+ip;
  pt->flag = base;
  ilist++;

  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4;
    pt = &mesh->tetra[k];

    /* If pt bears geometric information, search for endpoints of the NOM curve of ppt */
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      for (j=0; j<3; j++) {
        ie = MMG5_arpt[i][j];
        if ( pxt->tag[ie] & MG_NOM ) {
          na = pt->v[MMG5_iare[ie][0]];
          nb = pt->v[MMG5_iare[ie][1]];
          /* Store nb, if need be */
          if ( na == nump ) {
            if ( ip0 == 0 ) ip0 = nb;
            else if ( ip1 == 0 ) {
              if ( ip0 != nb ) ip1 = nb;
            }
            else {
              if ( ip0 != nb && ip1 != nb ) return 0;
            }
          }
          /* Store na, if need be */
          else {
            if ( ip0 == 0 ) ip0 = na;
            else if ( ip1 == 0 ) {
              if ( ip0 != na ) ip1 = na;
            }
            else {
              if ( ip0 != na && ip1 != na ) return 0;
            }
          }
        }
      }
    }

    /* Pile up tetrahedra in the ball of nump */
    adja = &mesh->adja[4*(k-1)+1];

    for (j=0; j<3; j++) {
      i = MMG5_inxt3[i];
      kk = adja[i] / 4;
      assert ( kk && "point is not an internal nm-point");

      if ( !kk ) {
        /* point is not an internal non manifold point */
        return 0;
      }

      pt1 = &mesh->tetra[kk];
      if ( pt1->flag == base ) continue;

      for (ii=0; ii<4; ii++)
        if ( pt1->v[ii] == nump ) break;
      assert ( ii < 4 );

      list[ilist] = 4*kk+ii;
      pt1->flag = base;
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      ilist++;
    }

    cur++;
  }

  /* At this point, the two points connected to ppt via the NOM curve are ip0 and ip1 */
  MMG3D_compute_tangent(mesh,nump,ip0,ip1,t);

  dd = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
  if ( dd > MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    t[0] *= dd;
    t[1] *= dd;
    t[2] *= dd;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param hash pointer to an allocated hash table.
 * \param start index of the starting tetrahedra.
 * \param ip local index of the point in the tetrahedra \a start.
 * \param ng pointer to the number of ridges.
 * \param nr pointer to the number of reference edges.
 * \param nm pointer to the number of non-manifold edges.
 * \return ns the number of special edges passing through ip, -1 if fail.
 *
 * Count the numer of ridges and reference edges incident to
 * the vertex \a ip when ip is non-manifold.
 *
 */
int MMG5_boulernm(MMG5_pMesh mesh,MMG5_Hash *hash,MMG5_int start,int ip,MMG5_int *ng,MMG5_int *nr,MMG5_int *nm){
  MMG5_pTetra    pt,pt1;
  MMG5_pxTetra   pxt;
  MMG5_hedge    *ph;
  MMG5_int       *adja,nump,k,k1;
  int            ns,ilist,cur;
  MMG5_int       list[MMG3D_LMAX+2],base,ia,ib,a,b,key,jj;
  int8_t         j,l,i;
  uint8_t        ie;

  /* reset the hash table */
  for ( k=0;  k<=hash->max; ++k ) {
    hash->item[k].a = 0;
    hash->item[k].b = 0;
  }

  for ( k=0;  k<=hash->siz; ++k ) {
    hash->item[k].nxt = 0;
  }
  for (k=hash->siz; k<hash->max; k++) {
    hash->item[k].nxt = k+1;
  }

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  nump = pt->v[ip];

  /* Store initial tetrahedron */
  pt->flag = base;
  list[0] = 4*start + ip;
  ilist = 1;

  *ng = *nr = *nm = ns = 0;

  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4; // index of point p in tetra k
    pt = &mesh->tetra[k];

    /* Count the number of ridge of ref edges passing through ip. */
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      for (l=0; l<3; ++l) {
        ie = MMG5_arpt[i][l];

        if ( MG_EDG(pxt->tag[ie]) || (MG_NOM & pxt->tag[ie]) ) {
          /* Seek if we have already seen the edge. If not, hash it and
           * increment ng or nr.*/
          a = pt->v[MMG5_iare[ie][0]];
          b = pt->v[MMG5_iare[ie][1]];
          ia  = MG_MIN(a,b);
          ib  = MG_MAX(a,b);
          key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
          ph  = &hash->item[key];

          if ( ph->a == ia && ph->b == ib )
            continue;
          else if ( ph->a ) {
            while ( ph->nxt && ph->nxt < hash->max ) {
              ph = &hash->item[ph->nxt];
              if ( ph->a == ia && ph->b == ib )  continue;
            }
            ph->nxt   = hash->nxt;
            ph        = &hash->item[hash->nxt];

            if ( hash->nxt >= hash->max-1 ) {
              if ( mesh->info.ddebug )
                fprintf(stderr,"\n  ## Warning: %s: memory alloc problem (edge):"
                        " %" MMG5_PRId "\n",__func__,hash->max);
              MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,MMG5_GAP,MMG5_hedge,
                                 "MMG5_edge",return -1);
              /* ph pointer may be false after realloc */
              ph        = &hash->item[hash->nxt];

              for (jj=ph->nxt; jj<hash->max; jj++)  hash->item[jj].nxt = jj+1;
            }
            hash->nxt = ph->nxt;
          }

          /* insert new edge */
          ph->a = ia;
          ph->b = ib;
          ph->nxt = 0;

          /* Order of following tests impacts the ridge and non-manifold edges
           * count (an edge that has both tags pass only in first test) but
           * should not influence the setting for corners and required tags in
           * setVertexNmTag function) */
          if ( pxt->tag[ie] & MG_GEO ) {
            ++(*ng);
          }
          else if ( pxt->tag[ie] & MG_NOM ) {
            ++(*nm);
          }
          else if ( pxt->tag[ie] & MG_REF ) {
            ++(*nr);
          }
          ++ns;
        }
      }
    }

    /* Continue to travel */
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = MMG5_inxt3[i];
      k1 = adja[i];
      if ( !k1 )  continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }

  return ns;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of the starting tetra.
 * \param ip index in \a start of the looked point.
 * \param iface index in \a start of the starting face.
 * \param listv pointer to the computed volumic ball.
 * \param ilistv pointer to the computed volumic ball size.
 * \param lists pointer to the computed surfacic ball.
 * \param ilists pointer to the computed surfacic ball size.
 * \param isnm 1 if \a ip is non-manifold, 0 otherwise.
 * \return -1 if fail, 1 otherwise.
 *
 * Compute the volumic ball of a SURFACE point \a p, as well as its surfacic
 * ball, starting from tetra \a start, with point \a ip, and face \a if in tetra
 * volumic ball:
 *   - \a listv[k] = 4* tet index + index of point surfacic ball.
 *   - \a lists[k] = 4* tet index + index of boundary face.
 *
 * \warning Don't work for a non-manifold point if \a start has an adjacent
 * through \a iface (for example : a non-manifold subdomain). Thus, if \a ip is
 * non-manifold, must be called only if \a start has no adjacent through iface.
 *
 */
int MMG5_boulesurfvolp(MMG5_pMesh mesh,MMG5_int start,int ip,int iface,
                        int64_t *listv,int *ilistv,MMG5_int *lists,int*ilists, int isnm)
{
  MMG5_pTetra   pt,pt1;
  MMG5_pxTetra  pxt;
  MMG5_int      k,*adja,nump,k1,fstart,piv,na,nb,adj,nvstart,aux,cur,base;
  int8_t        iopp,ipiv,i,j,l,isface;
  static int8_t mmgErr0=0, mmgErr1=0, mmgErr2=0;

  if ( isnm ) assert(!mesh->adja[4*(start-1)+iface+1]);

  base = ++mesh->base;
  *ilists = 0;
  *ilistv = 0;

  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  k = start;

  na = pt->v[ip];
  nb = pt->v[MMG5_idir[iface][MMG5_inxt2[MMG5_idirinv[iface][ip]]]];
  piv = pt->v[MMG5_idir[iface][MMG5_iprv2[MMG5_idirinv[iface][ip]]]];

  iopp = iface;
  fstart = 4*k+iopp;

  do {
    /* A boundary face has been hit : change travel edge */
    lists[(*ilists)] = 4*k+iopp;
    (*ilists)++;

    assert ( mesh->tetra[k].xt && "tetra of surfacic ball has a xtetra (bdy face) ");

    if ( *ilists >= MMG3D_LMAX ) {
      if ( !mmgErr0 ) {
        fprintf(stderr,"\n  ## Warning: %s: problem in surface remesh process."
                " Surface ball of at least 1 point (%" MMG5_PRId ") contains too"
                " many elts.\n"
                "  ##          Try to modify the hausdorff number "
                " or/and the maximum edge size.\n",__func__,
                MMG3D_indPt(mesh,nump));
        mmgErr0 = 1;
      }

      return -1;
    }

    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included)*/
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        listv[(*ilistv)] = 4*k+i;
        (*ilistv)++;
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      if ( !MMG3D_findEdge(mesh,pt,k,na,nb,0,&mmgErr2,&i) ) return -1;

      /* set sense of travel */
      if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
        iopp = MMG5_ifar[i][0];
        ipiv = MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      else {
        ipiv = MMG5_ifar[i][0];
        iopp = MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      if ( isnm ) {
        isface = (adja[iopp] == 0);
      }
      else {
        isface = 0;
        if(pt->xt){
          pxt = &mesh->xtetra[pt->xt];
          isface = (MG_BDY & pxt->ftag[iopp]);
        }
      }
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  /* Now, surfacic ball is complete ; finish travel of volumic ball */
  cur = 0;  // Check numerotation
  while ( cur < (*ilistv) ) {
    k = listv[cur]/4;
    i = listv[cur]%4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = MMG5_inxt3[i];
      k1 = adja[i];
      if ( !k1 )  continue;
      k1/=4;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;

      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);

      /* overflow */
      if ( *ilistv > MMG3D_LMAX-3 ) {
        if ( !mmgErr1 ) {
          fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
                  " Volumic ball of point %" MMG5_PRId " contains too many elts.\n",
                  __func__,MMG3D_indPt(mesh,nump));
          fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
                  " or/and the maximum mesh.\n");
          mmgErr1 = 1;
        }
        return -1;
      }
      listv[(*ilistv)] = 4*k1+j;
      (*ilistv)++;
    }
    cur++;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of the starting tetrahedron.
 * \param ip index in \a start of the desired vertex.
 * \param iface index in \a start of the starting face.
 * \param listv pointer to the computed volumic ball.
 * \param ilistv pointer to the computed volumic ball size.
 * \param lists pointer to the computed surfacic ball.
 * \param ilists pointer to the computed surfacic ball size.
 * \param refmin return the reference of one of the two subdomains in presence
 * \param refplus return the reference of the other subdomain in presence
 * \param isnm is the vertex non-manifold?
 * \return 1 if succesful, a negative value if the ball cannot be computed:
 * -1 if a surface ball had too many elements,
 * -2 if there are more than two references around,
 * -3 if an edge cannot be found, and
 * -4 if a volume ball had too many elements.
 * Among these, -1, -3 and -4 can be taken as a sign that further remeshing
 * is not possible, while -2 just mean the job could not be done.
 *
 * Compute the volumic ball of a SURFACE point \a p, as well as its surfacic
 * ball, starting from tetra \a start, with point \a ip, and face \a if in the
 * volumic ball.
 * \a listv[k] = 4*number of tet + index of point surfacic ball.
 * \a lists[k] = 4*number of tet + index of face.
 *
 * \warning Doesn't work for a non-manifold point if \a start has an adjacent
 * through \a iface (for example: a non-manifold subdomain). Thus, if \a ip is
 * non-manifold, must be called only if \a start has no adjacent through iface.
 *
 */
int MMG5_boulesurfvolpNom(MMG5_pMesh mesh, MMG5_int start, int ip, int iface,
                          int64_t *listv, int *ilistv, MMG5_int *lists, int *ilists,
                          MMG5_int *refmin, MMG5_int *refplus, int isnm)
{
  MMG5_pTetra   pt, pt1;
  MMG5_pxTetra  pxt;
  MMG5_int      k, k1, nump, *adja, piv, na, nb, adj, cur, nvstart, fstart, aux, base;
  int8_t        iopp, ipiv, i, j, l, isface;
  static int8_t mmgErr0=0, mmgErr1=0, mmgErr2=0;

  if ( isnm ) assert(!mesh->adja[4*(start-1)+iface+1]);

  base = ++mesh->base;
  *ilists  = 0;
  *ilistv  = 0;
  *refmin  = -1;
  *refplus = -1;
  
  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  k = start;
  
  na = pt->v[ip];
  nb = pt->v[MMG5_idir[iface][MMG5_inxt2[MMG5_idirinv[iface][ip]]]];
  piv = pt->v[MMG5_idir[iface][MMG5_iprv2[MMG5_idirinv[iface][ip]]]];
  
  iopp = iface;
  fstart = 4*k+iopp;
  
  do {
    /* A boundary face has been hit : change travel edge */
    lists[(*ilists)] = 4*k+iopp;
    (*ilists)++;
    if ( *ilists >= MMG3D_LMAX ) {
      if ( !mmgErr0 ) {
        fprintf(stderr,"\n  ## Warning: %s: problem in surface remesh process."
                " Surface ball of at least 1 point (%" MMG5_PRId ") contains too"
                " many elts.\n"
                "  ##          Try to modify the hausdorff number "
                " or/and the maximum edge size.\n",__func__,
                MMG3D_indPt(mesh,nump));
        mmgErr0 = 1;
      }
      return -1;
    }

    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included)*/
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        listv[(*ilistv)] = 4*k+i;
        (*ilistv)++;

        /* Identify references of both subdomains in presence */
        if ( *refmin == -1 )
          *refmin = pt->ref;
        else {
          if ( *refplus == -1 ) {
            if ( pt->ref != *refmin  ) *refplus = pt->ref;
          }
          else if ( pt->ref != *refmin && pt->ref != *refplus ) return -2;
        }
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      if ( !MMG3D_findEdge(mesh,pt,k,na,nb,0,&mmgErr2,&i) ) return -3;

      /* set sense of travel */
      if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
        iopp = MMG5_ifar[i][0];
        ipiv = MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      else {
        ipiv = MMG5_ifar[i][0];
        iopp = MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      if ( isnm ) {
        isface = (adja[iopp] == 0);
      }
      else {
        isface = 0;
        if(pt->xt){
          pxt = &mesh->xtetra[pt->xt];
          isface = (MG_BDY & pxt->ftag[iopp]);
        }
      }
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  /* Now, surfacic ball is complete ; finish travel of volumic ball */
  cur = 0;  // Check numerotation
  while ( cur < (*ilistv) ) {
    k = listv[cur]/4;
    i = listv[cur]%4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = MMG5_inxt3[i];
      k1 = adja[i];
      if ( !k1 )  continue;
      k1/=4;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;

      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);

      /* overflow */
      if ( *ilistv > MMG3D_LMAX-3 ) {
        if ( !mmgErr1 ) {
          fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
                  " Volumic ball of point %" MMG5_PRId " contains too many elts.\n",
                  __func__,MMG3D_indPt(mesh,nump));
          fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
                  " or/and the maximum mesh.\n");
          mmgErr1 = 1;
        }
        return -4;
      }
      listv[(*ilistv)] = 4*k1+j;
      (*ilistv)++;

      /* Identify references of both subdomains in presence */
      if ( *refmin == -1 )
        *refmin = pt1->ref;
      else {
        if ( *refplus == -1 ) {
          if ( pt1->ref != *refmin  ) *refplus = pt1->ref;
        }
        else if ( pt1->ref != *refmin && pt1->ref != *refplus ) return -2;
      }
    }
    cur++;
  }

  return 1;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param start index of the starting tetrahedron.
 * \param ip index of the looked ridge point.
 * \param iface index in \a start of the starting face.
 * \param il1 pointer to the first ball size.
 * \param l1 pointer to the first computed ball (associated to \a n_1's
 * side).
 * \param il2 pointer to the second ball size.
 * \param l2 pointer to the second computed ball (associated to \a n_2's
 * side).
 * \param ip0 index of the first extremity of the ridge.
 * \param ip1 index of the second extremity of the ridge.
 * \return 0 if fail, 1 otherwise.
 *
 * Computation of the two surface balls of a ridge point: the list \a l1 is
 * associated to the normal of face \a iface. \a ip0 and \a ip1 are the indices
 * of the 2 ending point of the ridge. Both lists are returned enumerated in
 * direct order.
 *
 */
int MMG5_bouletrid(MMG5_pMesh mesh,MMG5_int start,int iface,int ip,int *il1,MMG5_int *l1,
                    int *il2,MMG5_int *l2,MMG5_int *ip0,MMG5_int *ip1)
{
  MMG5_pTetra          pt;
  MMG5_pxTetra         pxt;
  MMG5_pPoint          ppt;
  MMG5_int             k,*adja,*list1,*list2,aux;
  MMG5_int             na, nb, piv,lists[MMG3D_LMAX+2], base;
  MMG5_int             idp, fstart, nvstart, adj;
  int                  ilists, iopp, ipiv,*ilist1,*ilist2;
  int                  ifac,idx,idx2,idx_tmp,i1,isface;
  double               *n1,*n2,nt[3],ps1,ps2;
  int8_t               i;
  static int8_t        mmgErr0=0,mmgErr1=0;

  pt = &mesh->tetra[start];
  if ( !MG_EOK(pt) )  return 0;

#ifndef NDEBUG
  assert(pt->xt);
  pxt = &mesh->xtetra[pt->xt];
  // We must call this function on a well orientated boundary face (to build the
  // direct surfacic ball).
  assert(pxt->ftag[iface] & MG_BDY);
  assert( MG_GET(pxt->ori,iface) );
#endif

  idp = pt->v[ip];
  k = start;

  ppt = &mesh->point[idp];
  assert( ppt->tag & MG_GEO );

  na  = pt->v[ip];
  nb  = pt->v[MMG5_idir[iface][MMG5_inxt2[MMG5_idirinv[iface][ip]]]];
  piv = pt->v[MMG5_idir[iface][MMG5_iprv2[MMG5_idirinv[iface][ip]]]];

  iopp = iface;
  fstart = 4*k+iopp;

  base = ++mesh->base;

  /* Set pointers on lists il1 and il2 to have il1 associated to the normal of
     the face iface.*/
  MMG5_norpts(mesh, pt->v[MMG5_idir[iface][0]],pt->v[MMG5_idir[iface][1]],
               pt->v[MMG5_idir[iface][2]],nt);

  n1 = &(mesh->xpoint[ppt->xp].n1[0]);
  n2 = &(mesh->xpoint[ppt->xp].n2[0]);
  ps1 = n1[0]*nt[0] + n1[1]*nt[1] + n1[2]*nt[2];
  ps2 = n2[0]*nt[0] + n2[1]*nt[1] + n2[2]*nt[2];

  if ( fabs(ps1) < fabs(ps2) ) {
    list1  = l2;
    list2  = l1;
    ilist1 = il2;
    ilist2 = il1;
  }
  else {
    list1  = l1;
    list2  = l2;
    ilist1 = il1;
    ilist2 = il2;
  }
  *ilist1 = 0;
  *ilist2 = 0;

  /* First: fill the surfacic ball. */
  ilists = 0;
  do {
    /* A boundary face has been hit : change travel edge */
    lists[ilists] = 4*k+iopp;
    ilists++;
    if ( ilists >= MMG3D_LMAX ) {
      if ( !mmgErr0 ) {
        fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
                " Volumic ball of point %" MMG5_PRId " contains too many elts.\n",
                __func__,MMG3D_indPt(mesh,idp));
        fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
                " or/and the maximum mesh.\n");
        mmgErr0 = 1;
      }
      return 0;
    }

    aux = nb;
    nb = piv;
    piv = aux;
    nvstart = k;
    adj = k;

    /* Now unfold shell of edge (na,nb) starting from k (included) */
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      pt->flag = base;

      /* identification of edge number in tetra k */
      if ( !MMG3D_findEdge(mesh,pt,k,na,nb,0,&mmgErr1,&i) ) return -1;

      /* set sense of travel */
      if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
        iopp = MMG5_ifar[i][0];
        ipiv = MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      else {
        ipiv = MMG5_ifar[i][0];
        iopp = MMG5_ifar[i][1];
        adj = adja[ iopp ] / 4;
        piv = pt->v[ipiv];
      }
      isface = 0;
      if(pt->xt){
        pxt = &mesh->xtetra[pt->xt];
        isface = (MG_BDY & pxt->ftag[iopp]);
      }
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  // Remark: here the test k!=start is a security bound: theorically it is
  // useless but in case of bad edge tag, it ensure that the loop is not
  // infinite.
  while ( 4*k+iopp != fstart );

  /* Second: travel through the surface ball until meeting a ridge. */
  for (idx=0; idx!=ilists; ++idx) {
    k    = lists[idx]/4;
    ifac = lists[idx]%4;
    pt   = &mesh->tetra[k];
    assert(pt->xt);
    pxt = &mesh->xtetra[pt->xt];

    for ( i=0; i<3; ++i ) {
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i1 = MMG5_inxt2[i];
    if ( pxt->tag[MMG5_iarf[ifac][i1]]  & MG_GEO ) break;
  }
  assert(idx < ilists);
  *ip0 = pt->v[MMG5_idir[ifac][MMG5_iprv2[i]]];

  /* Start from the hit boundary, until another boundary is hit and complete the
   * second ball */
  idx = (idx+1)%ilists;
  for (idx2=idx; idx2!=ilists+idx; ++idx2) {
    idx_tmp = idx2%ilists;
    k       = lists[idx_tmp]/4;
    ifac    = lists[idx_tmp]%4;
    pt      = &mesh->tetra[k];
    assert(pt->xt);
    pxt     = &mesh->xtetra[pt->xt];

    if ( (*ilist2) > MMG3D_LMAX-2 )  return 0;
    list2[(*ilist2)] = 4*k+ifac;
    (*ilist2)++;

    for ( i=0; i<3; ++i ) {
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i1 = MMG5_inxt2[i];
    if ( pxt->tag[MMG5_iarf[ifac][i1]]  & MG_GEO ) break;
  }
  assert(idx2 != ilists+idx);
  *ip1 = pt->v[MMG5_idir[ifac][MMG5_iprv2[i]]];

  /* Start again from the newly hit boundary, until another boundary is hit and
   * complete the first ball */
  idx = (idx2+1)%ilists;
  for (idx2=idx; idx2 != idx+ilists; ++idx2) {
    idx_tmp = idx2%ilists;
    k       = lists[idx_tmp]/4;
    ifac    = lists[idx_tmp]%4;
    pt      = &mesh->tetra[k];
    assert(pt->xt);
    pxt     = &mesh->xtetra[pt->xt];

    if ( (*ilist1) > MMG3D_LMAX-2 )  return 0;
    list1[(*ilist1)] = 4*k+ifac;
    (*ilist1)++;

    for ( i=0; i<3; ++i ) {
      if ( pt->v[MMG5_idir[ifac][i]] == idp ) break;
    }
    assert(i<3);

    i1 = MMG5_inxt2[i];
    if ( pxt->tag[MMG5_iarf[ifac][i1]]  & MG_GEO ) break;
  }
  assert(idx2 != ilists+idx);
  assert(*ip0 == pt->v[MMG5_idir[ifac][MMG5_iprv2[i]]]);

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param start tetra from which we start to travel
 * \param na edge vertex
 * \param nb edge vertex
 * \param tag new edge tag
 * \param edg new edge ref
 * \param piv global index of the pivot to set the sense of travel
 * \param adj index of adjacent tetra for the travel
 *
 * \return -1 if fail, \a start if shell has been completely travelled, 0 otherwise
 *
 * Set tag and ref of the edge \a na \a nb from tetra \a start by traveling
 * its shell in one direction (given by the pivot \a piv).
 *
 */
static inline
int MMG3D_settag_oneDir(MMG5_pMesh  mesh,MMG5_int start, MMG5_int na, MMG5_int nb,
                                     uint16_t tag,int edg, MMG5_int piv,MMG5_int adj) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     *adja;
  uint16_t     taginit;
  int8_t       i;

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,1,NULL,&i) ) {
      return -1;
    }

    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];

      taginit = pxt->tag[i];
      pxt->tag[i] |= tag;
      /* Remove the potential nosurf tag if initially the edge is
       * really required */
      if ( ((taginit & MG_REQ) && !(taginit & MG_NOSURF)) ||
           ((    tag & MG_REQ) && !(    tag & MG_NOSURF)) ) {
        pxt->tag[i] &= ~MG_NOSURF;
      }
      pxt->edg[i]  = MG_MAX(pxt->edg[i],edg);
    }
    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ MMG5_ifar[i][1] ] /4;
      piv = pt->v[ MMG5_ifar[i][0] ];
    }
  }
  return adj;
}

/**
 * \param mesh pointer to the mesh structure
 * \param start tetra from which we start
 * \param ia local index of the edge in \a start
 * \param tag tag to set
 * \param edg edge reference to set
 *
 * \return 1 if success, 0 if fail.
 *
 * Set tag \a tag and ref \a edg of edge \a ia (if need be) in tetra \a start by
 * travelling its shell.
 *
 */
int MMG5_settag(MMG5_pMesh mesh,MMG5_int start,int ia,uint16_t tag,int edg) {
  MMG5_pTetra        pt;
  MMG5_pxTetra       pxt;
  MMG5_int           na,nb,*adja,adj,piv;
  uint16_t           taginit;

  assert( start >= 1 );
  pt = &mesh->tetra[start];
  assert ( MG_EOK(pt) );

  na   = pt->v[ MMG5_iare[ia][0] ];
  nb   = pt->v[ MMG5_iare[ia][1] ];

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][0]] / 4;
  piv = pt->v[MMG5_ifar[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    taginit = pxt->tag[ia];
    pxt->tag[ia] |= tag;
    /* Remove the potential nosurf tag if initially the edge is
     * really required */
    if ( ((taginit & MG_REQ) && !(taginit & MG_NOSURF)) ||
         ((    tag & MG_REQ) && !(    tag & MG_NOSURF)) ) {
      pxt->tag[ia] &= ~MG_NOSURF;
    }
    pxt->edg[ia]  = MG_MAX(pxt->edg[ia],edg);
  }

  adj = MMG3D_settag_oneDir(mesh,start,na,nb,tag,edg,piv,adj);

  /* If all shell has been travelled, stop, else, travel it the other sense */
  if ( adj == start )  return 1;
  else if ( adj < 0 ) return 0;

  assert(!adj);

  pt = &mesh->tetra[start];
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][1]] / 4;
  piv = pt->v[MMG5_ifar[ia][0]];

  adj = MMG3D_settag_oneDir(mesh,start,na,nb,tag,edg,piv,adj);

  if ( adj < 0 ) return 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param start tetra from which we start to travel
 * \param na edge vertex
 * \param nb edge vertex
 * \param tag new edge tag
 * \param piv global index of the pivot to set the sense of travel
 * \param adj index of adjacent tetra for the travel
 *
 * \return -1 if fail, \a start if shell has been completely travelled, 0 otherwise
 *
 * Remove the tag \a tag of edge \a ia in tetra \a start by travelling its
 * shell in one direction (given by the pivot \a piv).
 *
 */
static inline
int MMG3D_deltag_oneDir(MMG5_pMesh  mesh,MMG5_int start, MMG5_int na, MMG5_int nb,
                        uint16_t tag,MMG5_int piv,MMG5_int adj) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     *adja;
  int8_t       i;

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,1,NULL,&i) ) {
      return -1;
    }

    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      if ( (pxt->ftag[MMG5_ifar[i][0]] & MG_BDY) ||
           (pxt->ftag[MMG5_ifar[i][1]] & MG_BDY) ) {
        pxt->tag[i] &= ~tag;
      }
    }
    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ MMG5_ifar[i][0] ] / 4;
      piv = pt->v[ MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ MMG5_ifar[i][1] ] /4;
      piv = pt->v[ MMG5_ifar[i][0] ];
    }
  }
  return adj;
}

/**
 * \param mesh pointer to the mesh structure
 * \param start index of the starting tetra
 * \param ia index of the edge in tetra \a start that we want to modify
 * \param tag tag to remove
 * \return 1 if success, 0 otherwise.
 *
 * Remove the tag \a tag of edge \a ia in tetra \a start by travelling its
 * shell.
 *
 */
int MMG5_deltag(MMG5_pMesh mesh,MMG5_int start,int ia,uint16_t tag) {
  MMG5_pTetra        pt;
  MMG5_pxTetra       pxt;
  MMG5_int           na,nb,*adja,adj,piv;
  int8_t             i;

  assert( start >= 1 );
  pt = &mesh->tetra[start];
  assert ( MG_EOK(pt) );

  na   = pt->v[ MMG5_iare[ia][0] ];
  nb   = pt->v[ MMG5_iare[ia][1] ];

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][0]] / 4;
  piv = pt->v[MMG5_ifar[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( (pxt->ftag[MMG5_ifar[ia][0]] & MG_BDY) ||
         (pxt->ftag[MMG5_ifar[ia][1]] & MG_BDY) ) {
      pxt->tag[ia] &= ~tag;
    }
  }

  adj = MMG3D_deltag_oneDir(mesh,start,na,nb,tag,piv,adj);

  /* If all shell has been travelled, stop, else, travel it the other sense */
  if ( adj == start )  return 1;
  else if ( adj < 0 ) return 0;

  assert(!adj);

  pt = &mesh->tetra[start];
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][1]] / 4;
  piv = pt->v[MMG5_ifar[ia][0]];

  adj = MMG3D_deltag_oneDir(mesh,start,na,nb,tag,piv,adj);

  if ( adj < 0 ) return 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param start index of the starting tetra
 * \param ia index of the edge
 * \param list list of tetra sharing the edge \a ia
 * \param isbdy 1 if edge is bdy, 0 otherwise (note that at interface of
 * 2 domains the edge shell of a bdy edge can be closed)
 *
 * \return 2*ilist if shell is closed, 2*ilist +1 otherwise, 0 if one of the tet
 * of the shell is required, -1 if fail.
 *
 * Find all tets sharing edge ia of tetra start.
 *
 */
int MMG5_coquil(MMG5_pMesh mesh,MMG5_int start,int ia,int64_t*list,int8_t *isbdy) {
  MMG5_pTetra   pt;
  MMG5_int      *adja,piv,na,nb,adj;
  int           ilist;
  int8_t        i;
  static int8_t mmgErr0=0, mmgErr1=0;

  assert ( start >= 1 );
  pt = &mesh->tetra[start];
  assert ( MG_EOK(pt) );

  na   = pt->v[ MMG5_iare[ia][0] ];
  nb   = pt->v[ MMG5_iare[ia][1] ];
  ilist = 0;
  list[ilist] = 6*(int64_t)start+ia;
  ilist++;

  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][0]] / 4; // start travelling by face (ia,0)
  piv = pt->v[MMG5_ifar[ia][1]];
  *isbdy = (pt->xt && (mesh->xtetra[pt->xt].ftag[MMG5_ifar[ia][0]] & MG_BDY))? 1 : 0;

  while ( adj && (adj != start) ) {
    pt = &mesh->tetra[adj];
    if ( pt->tag & MG_REQ )  return 0;

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,0,&mmgErr1,&i) ) return -1;

    list[ilist] = 6*(int64_t)adj +i;
    ilist++;
    /* overflow */
    if ( ilist > MMG3D_LMAX-3 ) {
      if ( !mmgErr0 ) {
        fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
                " Coquil of edge %" MMG5_PRId "-%" MMG5_PRId " contains too many elts.\n",
                __func__,MMG3D_indPt(mesh,na),MMG3D_indPt(mesh,nb));
        fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
                " or/and the maximum mesh.\n");
        mmgErr0 = 1;
      }
      return -1;
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    int8_t travel_fac = MMG5_ifar[i][0];
    if ( pt->v[ travel_fac ] == piv ) {
      adj = adja[ travel_fac ] / 4;
      piv = pt->v[ MMG5_ifar[i][1] ];
    }
    else {
      travel_fac = MMG5_ifar[i][1];
      assert(pt->v[ travel_fac ] == piv );
      adj = adja[ travel_fac ] /4;
      piv = pt->v[ MMG5_ifar[i][0] ];
    }

    /* If we haven't cross yet a bdy face, test traveled triangle. Avoid the
     * test otherwise (to not add a useless access to the xtetra strcuture) */
    if ( !*isbdy ) {
      if ( pt->xt && (mesh->xtetra[pt->xt].ftag[travel_fac] & MG_BDY) ) {
        *isbdy = 1;
      }
    }
  }

  /* At this point, the first travel, in one direction, of the shell is
     complete. Now, analyze why the travel ended. */
  if ( adj == start )  return 2*ilist;
  assert(!adj); // a boundary has been detected

  adj = list[ilist-1] / 6;
  i   = list[ilist-1] % 6;
  ilist = 0;
  *isbdy = 1;

  /* Start back everything from this tetra adj */
  list[ilist] = 6*(int64_t)adj + i;
  ilist++;
  /* overflow */
  if ( ilist > MMG3D_LMAX-3 ) {
    if ( !mmgErr0 ) {
      fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
              " Coquil of edge %" MMG5_PRId "-%" MMG5_PRId " contains too many elts.\n",
              __func__,MMG3D_indPt(mesh,na),MMG3D_indPt(mesh,nb));
      fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
              " or/and the maximum mesh.\n");
      mmgErr0 = 1;
    }
    return -1;
  }

  adja = &mesh->adja[4*(adj-1)+1];
  if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
    adj = adja[ MMG5_ifar[i][0] ];
    piv = pt->v[ MMG5_ifar[i][1] ];
  }
  else {
    adj = adja[ MMG5_ifar[i][1] ];
    piv = pt->v[ MMG5_ifar[i][0] ];
  }

  while ( adj ) {
    adj /= 4;
    pt = &mesh->tetra[adj];
    if ( pt->tag & MG_REQ )  return 0;

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,0,&mmgErr1,&i) ) return -1;

    list[ilist] = 6*(int64_t)adj +i;
    ilist++;
    /* overflow */
    if ( ilist > MMG3D_LMAX-2 ) {
      if ( !mmgErr0 ) {
        fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
                " Coquil of edge %" MMG5_PRId "-%" MMG5_PRId " contains too many elts.\n",
                __func__,MMG3D_indPt(mesh,na),MMG3D_indPt(mesh,nb));
        fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
                " or/and the maximum mesh.\n");
        mmgErr0 = 1;
      }
      return -1;
    }

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
      adj = adja[ MMG5_ifar[i][0] ];
      piv = pt->v[ MMG5_ifar[i][1] ];
    }
    else {
      adj = adja[ MMG5_ifar[i][1] ];
      piv = pt->v[ MMG5_ifar[i][0] ];
    }
  }
  assert(!adj);
  return  2*ilist+1 ;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param start starting tetra.
 * \param ia local edge index in tetra \a start.
 *
 * \return 1 if the edge \a ia in \a start is boundary, 0 otherwise, -1 if fail.
 *
 * Identify whether edge ia in start is a boundary edge by unfolding its shell.
 *
 */
int MMG5_srcbdy(MMG5_pMesh mesh,MMG5_int start,int ia) {
  MMG5_pTetra      pt;
  MMG5_pxTetra     pxt;
  MMG5_int         na,nb,adj,piv,*adja;
  int8_t           iadj,i;

  pt = &mesh->tetra[start];
  na = pt->v[MMG5_iare[ia][0]];
  nb = pt->v[MMG5_iare[ia][1]];

  adja = &mesh->adja[4*(start-1)+1];
  iadj = MMG5_ifar[ia][0];

  if(pt->xt){
    pxt = &mesh->xtetra[pt->xt];
    if( pxt->ftag[iadj] & MG_BDY )
      return 1;
  }

  adj = adja[iadj] / 4;
  piv = pt->v[MMG5_ifar[ia][1]];

  while( adj && ( adj != start ) ) {
    pt = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,1,NULL,&i) ) return -1;

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
      iadj = MMG5_ifar[i][0];
      adj = adja[ iadj ] / 4;
      piv = pt->v[ MMG5_ifar[i][1] ];
    }
    else {
      iadj = MMG5_ifar[i][1];
      adj = adja[ iadj ] /4;
      piv = pt->v[ MMG5_ifar[i][0] ];
    }

    if(pt->xt){
      pxt = &mesh->xtetra[pt->xt];
      if( pxt->ftag[iadj] & MG_BDY )
        return 1;
    }
  }

  return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param k1 should contain a tetra index.
 * \param k2 should contain a tetra index different from k2.
 *
 * Print an error message if MMG5_coquilFace detect a boundary topology problem.
 *
 */
 void MMG5_coquilFaceErrorMessage(MMG5_pMesh mesh, MMG5_int k1, MMG5_int k2) {
  MMG5_pTetra   pt;
  MMG5_int      kel1, kel2;
  static int8_t mmgErr0;

  if ( mmgErr0 ) return;

  mmgErr0 = 1;

  fprintf(stderr,"\n  ## Error: %s: at least 1 problem in surface"
          " remesh process",__func__);
  fprintf(stderr," (potential creation of a lonely boundary face):\n");

  kel1 = MMG3D_indElt(mesh,k1);
  kel2 = MMG3D_indElt(mesh,k2);

  if ( kel1 != 0 ) {
    pt = &mesh->tetra[k1];
    assert ( pt && MG_EOK(pt) );
    fprintf(stderr,"            look at elt %" MMG5_PRId ":",kel1);
    fprintf(stderr," %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId ".\n", MMG3D_indPt(mesh,pt->v[0]),
            MMG3D_indPt(mesh,pt->v[1]),MMG3D_indPt(mesh,pt->v[2]),
            MMG3D_indPt(mesh,pt->v[3]));
    fprintf(stderr,"            adjacent tetras %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
            MMG3D_indElt(mesh,mesh->adja[4*(k1-1)+1]/4),
            MMG3D_indElt(mesh,mesh->adja[4*(k1-1)+2]/4),
            MMG3D_indElt(mesh,mesh->adja[4*(k1-1)+3]/4),
            MMG3D_indElt(mesh,mesh->adja[4*(k1-1)+4]/4));
    fprintf(stderr,"            vertex required? %d %d %d %d\n",
            mesh->point[pt->v[0]].tag & MG_REQ,
            mesh->point[pt->v[1]].tag & MG_REQ,
            mesh->point[pt->v[2]].tag & MG_REQ,
            mesh->point[pt->v[3]].tag & MG_REQ);
  } else if ( kel2 != 0 ) {
    fprintf(stderr,"            look at elt %" MMG5_PRId ":",kel2);
    pt = &mesh->tetra[k2];
    assert ( pt && MG_EOK(pt) );

    fprintf(stderr," %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId ".\n\n",MMG3D_indPt(mesh,pt->v[0]),
            MMG3D_indPt(mesh,pt->v[1]),MMG3D_indPt(mesh,pt->v[2]),
            MMG3D_indPt(mesh,pt->v[3]));
  }
  fprintf(stderr,"\n  ##        Try to modify the hausdorff number,");
  fprintf(stderr," the maximum mesh size or/and the value of angle detection.\n");
  fprintf(stderr," You can also try to run with -noswap option but probably");
  fprintf(stderr," the final mesh will have poor quality.\n\n");
}

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of the starting tetrahedron.
 * \param na global index of the 1st extremity of the edge whose shell is computed
 * \param nb global index of the  2d extremity of the edge whose shell is computed
 * \param iface index of the face from which we come.
 * \param ia index of edge whose shell is computed (in tetra).
 * \param list pointer to the list of tetra in the shell (to fill).
 * \param ilist pointer to the number of tetra in the shell (to fill).
 * \param it1 pointer to the index of the 1st boundary face sharing \a ia
 * \param it2 pointer to the index of the 2d boundary face sharing \a ia
 * (to fill).
 * \param adj pointer to the adjacent to treat in the shell (to update)
 * \param hasadja pointer to 0 if we don't have adja through iface,
 * 1 otherwise (to fill)
 * \param nbdy pointer to the number of boundaries found minus 1 (to update)
 * \param silent if 1, print error message for more than 2 boundary triangles
 * in the shell
 *
 * \return -1 if fail, 1 otherwise
 *
 * Travel in the shell of the edge until meeting the first tetra or reaching a
 * tetra without adjacent. Fill \a it2 and \a list.
 *
 */
int MMG3D_coquilFaceFirstLoop(MMG5_pMesh mesh,MMG5_int start,MMG5_int na,MMG5_int nb,int8_t iface,
                               int8_t ia,int64_t *list,int *ilist,MMG5_int *it1,MMG5_int *it2,
                               MMG5_int *piv,MMG5_int *adj,int8_t *hasadja,int *nbdy,int silent) {

  MMG5_pTetra   pt;
  MMG5_int      pradj,*adja;
  int           pri,ier,ifar_idx;
  int8_t        i;
  static int8_t mmgErr0 = 0;

#ifndef NDEBUG
  MMG5_pxTetra  pxt;
#endif

  pt = &mesh->tetra[start];

  *ilist = 0;

  *it1 = 0;
  *it2 = 0;

  /* Ensure that the first boundary face found is ifac (nedded in multidomain case) */
  ifar_idx = (MMG5_ifar[ia][0]==iface) ? 1 : 0;
  assert ( iface == MMG5_ifar[ia][(ifar_idx+1)%2] );

  (*piv)  = pt->v[MMG5_ifar[ia][ifar_idx]];
  *adj    = start;
  i       = ia;

#ifndef NDEBUG
  pxt = &mesh->xtetra[pt->xt];
  assert ( MG_BDY & pxt->ftag[iface] );
#endif

  (*it1) = 4*start + iface;

  adja       = &mesh->adja[4*(start-1)+1];
  (*hasadja) = (adja[iface] > 0);

  (*nbdy)    = 0;

  do {
    pradj = (*adj);
    pri    = i;

    /* travel through new tetra */
    ier = MMG5_coquilTravel(mesh,na,nb,adj,piv,&iface,&i);

    /* fill the shell */
    list[(*ilist)] = 6*(int64_t)pradj +pri;
    (*ilist)++;

    /* overflow */
    if ( (*ilist) > MMG3D_LMAX-2 ) {
      if ( !mmgErr0 ) {
        fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
                " Coquil of edge %" MMG5_PRId "-%" MMG5_PRId " contains too many elts.\n",
                __func__,MMG3D_indPt(mesh,na),MMG3D_indPt(mesh,nb));
        fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
                " or/and the maximum mesh.\n");
        mmgErr0 = 1;
      }
      return -1;
    }

    if ( ier<0 ) return -1;
    else if ( !ier ) continue;

    if ( !(*it2) ) {
      *it2 = 4*pradj+iface;
      (*nbdy)++;
    }
    else {
      (*nbdy)++;
    }

  } while ( (*adj) && ((*adj) != start) );

  if ( (*adj) != start ) {
    /* The starting boundary face has not been counted (open shell) */
    ++(*nbdy);
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param piv global index of the pivot.
 * \param iface index of the face from which we come.
 * \param i index of edge whose shell is computed (in tetra).
 * \param list pointer to the list of tetra in the shell (to fill).
 * \param ilist pointer to the number of tetra in the shell (to fill).
 * \param it1 pointer to the index of the 1st boundary face sharing \a ia
 * \param pradj pointer to the first tetra of the shell (to fill).
 * \param adj pointer to the adjacent to treat in the shell (to update)
 *
 * Initialize the travel in the shell of the edge in reverse direction than in
 * the \a coquilFaceFirstLoop function.
 *
 */
void MMG3D_coquilFaceSecondLoopInit(MMG5_pMesh mesh,MMG5_int piv,int8_t *iface,
                                     int8_t *ia,int64_t *list,int *ilist,MMG5_int *it1,
                                     MMG5_int *pradj,MMG5_int *adj) {

  MMG5_pTetra   pt;
#ifndef NDEBUG
  MMG5_pxTetra  pxt;
#endif

  assert( !(*adj) );

  (*adj)      = list[(*ilist)-1] / 6;
  (*ia)       = list[(*ilist)-1] % 6;
  (*ilist)    = 0;

  (*pradj) = (*adj);
  pt       = &mesh->tetra[(*adj)];
#ifndef NDEBUG
  assert(pt->xt);
  pxt      = &mesh->xtetra[pt->xt];
#endif

  if ( pt->v[ MMG5_ifar[(*ia)][0] ] == piv ) {
    (*iface) = MMG5_ifar[(*ia)][1];
  }
  else {
    (*iface) = MMG5_ifar[(*ia)][0];
  }

  assert ( pxt->ftag[(*iface)] & MG_BDY );

  *it1 = 4*(*pradj) + (*iface);

}

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of the starting tetrahedron.
 * \param iface index of the boundary face from which we come.
 * \param ia index of edge whose shell is computed (in tetra).
 * \param list pointer to the list of tetra in the shell (to fill).
 * \param it1 pointer to the index of the first boundary face sharing \a ia
 * (to fill).
 * \param it2 pointer to the index of the second boundary face sharing \a ia
 * (to fill).
 * \param silent if 1, print error message for more than 2 boundary triangles
 * in the shell
 * \return -1 if fail, \f$2*ilist\f$ if shell is closed, \f$2*ilist+1\f$
 * otherwise.
 *
 * Find all tets sharing edge \a ia of tetra \a start, and stores boundary faces
 * when met. \f$ it1 \f$ and \f$ it2 = 6*iel + iface\f$, \a iel = index of
 * tetra, \a iface = index of face in tetra.
 *
 * \warning Don't work if \a ia has only one boundary face in its shell.
 */
int MMG5_coquilface(MMG5_pMesh mesh,MMG5_int start,int8_t iface,int ia,int64_t *list,
                     MMG5_int *it1,MMG5_int *it2, int silent) {
  MMG5_pTetra   pt;
  MMG5_int      piv,adj,na,nb,pradj;
  int           ier,nbdy,ilist;
  int8_t        hasadja,i;
  static int8_t mmgErr0=0,mmgErr1=0,mmgWarn0=0;

  pt = &mesh->tetra[start];

  na   = pt->v[ MMG5_iare[ia][0] ];
  nb   = pt->v[ MMG5_iare[ia][1] ];

  /* Travel throug the shell of the edge until reaching a tetra without adjacent
   * or until reaching the starting tetra */
  ier = MMG3D_coquilFaceFirstLoop(mesh,start,na,nb,iface,ia,list,&ilist,it1,it2,
                                   &piv,&adj,&hasadja,&nbdy,silent);

  if ( ier < 0 ) return ier;

  /* At this point, the first travel, in one direction, of the shell is
     complete. Now, analyze why the travel ended. */
  if ( adj == start ) {
    if ( !(*it2) ) {
      if ( !mmgErr0 ) {
        printf("  ## Error: %s: Wrong boundary tags: Only 1 boundary face found in"
               " the shell of the edge\n",__func__);
        mmgErr0 = 1;
      }
      return -1;
    }

    if ( nbdy != 2 ) {
      if ( nbdy < 2 ) {
        MMG5_coquilFaceErrorMessage(mesh, (*it1)/4, (*it2)/4);
        return -1;
      }

      if ( !silent ) {
        if ( !mmgWarn0 ) {
          // Algiane: for a manifold edge 2 cases :
          // 1) the shell is open and we have more than 3 tri sharing the edge
          // (highly non-manifold)
          // 2) we have a non-manifold shape immersed in a domain (3 triangles
          // sharing the edge and a closed shell)
          printf("  ## Warning: %s: you have %d boundary triangles in the closed shell"
                 " of a manifold edge.\n",__func__,nbdy);
          printf("  Problem may occur during remesh process.\n");
          mmgWarn0 = 1;

          /* MMG5_coquilface is called only on edges marked as manifold, check this */
          assert ( pt->xt );
          assert ( !(mesh->xtetra[pt->xt].tag[ia] & MG_NOM) );
        }
      }
    }

    return (2*ilist);
  }

  /* A boundary has been detected : slightly different configuration */
  if ( !hasadja ) return 2*ilist+1;

  /* Start back everything from this tetra adj */
  MMG3D_coquilFaceSecondLoopInit(mesh,piv,&iface,&i,list,&ilist,it1,
                                  &pradj,&adj);

  while ( adj ) {
    pradj = adj;

    ier = MMG5_openCoquilTravel( mesh, na, nb, &adj, &piv, &iface, &i );
    if ( ier<0 ) return -1;

    list[ilist] = 6*(int64_t)pradj +i;
    ilist++;
    /* overflow */
    if ( ilist > MMG3D_LMAX-2 ) {
      if ( !mmgErr1 ) {
        fprintf(stderr,"\n  ## Warning: %s: problem in remesh process."
                " Coquil of edge %" MMG5_PRId "-%" MMG5_PRId " contains too many elts.\n",
                __func__,MMG3D_indPt(mesh,na),MMG3D_indPt(mesh,nb));
        fprintf(stderr,"\n  ##          Try to modify the hausdorff number,"
                " or/and the maximum mesh.\n");
        mmgErr1 = 1;
      }
      return -1;
    }
  }

  assert(!adj);
  *it2 = 4*pradj + iface;

  if ( (!(*it1) || !(*it2)) || ((*it1) == (*it2)) ) {
    MMG5_coquilFaceErrorMessage(mesh, (*it1)/4, (*it2)/4);
    return -1;
  }
  return  2*ilist+1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param na global index of edge extremity.
 * \param nb global index of edge extremity.
 * \param adj starting tetrahedron at the begining and finish tet at the end.
 * \param piv global index of the vertex opposite to the travelling face
 * (updated for the finish tet at the end).
 * \param iface previous traveling face of the tet (suspected to be boundary),
 * updated.
 * \param i local index of the edge \f$[na,nb]\f$ in tet \a adj.
 * \return the tag of the face \a iface of the tetra \a adj, 0 if the tetra
 * is not boundary, -1 if fail.
 *
 * Travel around the edge \f$[na,nb]\f$ from tetra \a adj and through the face
 * \a piv.
 *
 */
int16_t MMG5_coquilTravel(MMG5_pMesh mesh, MMG5_int na, MMG5_int nb, MMG5_int* adj, MMG5_int *piv,
                           int8_t *iface, int8_t *i )
{
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     *adja;
  int16_t      isbdy;

  pt = &mesh->tetra[*adj];
  pxt = &mesh->xtetra[pt->xt];

  /* set new tetra for travel */
  adja = &mesh->adja[4*(*adj-1)+1];
  if ( pt->v[ MMG5_ifar[*i][0] ] == *piv ) {
    *iface = MMG5_ifar[*i][0];
    *adj = adja[ MMG5_ifar[*i][0] ] / 4;
    *piv = pt->v[ MMG5_ifar[*i][1] ];
  }
  else {
    assert(pt->v[ MMG5_ifar[*i][1] ] == *piv );
    *iface = MMG5_ifar[*i][1];
    *adj = adja[ MMG5_ifar[*i][1] ] /4;
    *piv = pt->v[ MMG5_ifar[*i][0] ];
  }
  isbdy = pt->xt ? (pxt->ftag[*iface] & MG_BDY) : 0;

  /* identification of edge number in tetra *adj */
  if ( *adj ) {
    pt = &mesh->tetra[*adj];
    if ( !MMG3D_findEdge(mesh,pt,*adj,na,nb,1,NULL,i) ) return -1;
  }

  return isbdy;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param na global index of edge extremity.
 * \param nb global index of edge extremity.
 * \param adj starting tetrahedron at the begining and finish tet at the end.
 * \param piv global index of the vertex opposite to the travelling face
 * (updated for the finish tet at the end).
 * \param iface traveling face of the tet (suspected to be boundary), updated.
 * \param i local index of the edge \f$[na,nb]\f$ in tet \a adj.
 *
 * \return 1 if success, 0 if fail.
 *
 * Travel around the edge \f$[na,nb]\f$ from tetra \a adj and through the face
 * \a piv. The shell of the edge is open and the tetra \a adj has no neighbour
 * through the face \a iface.
 *
 */
int16_t MMG5_openCoquilTravel(MMG5_pMesh mesh,MMG5_int na,MMG5_int nb,MMG5_int* adj,MMG5_int *piv,
                              int8_t *iface, int8_t *i )
{
  MMG5_pTetra  pt;
  MMG5_int     *adja;

  pt = &mesh->tetra[*adj];

  /* identification of edge number in tetra *adj */
  if ( !MMG3D_findEdge(mesh,pt,*adj,na,nb,1,NULL,i) ) return 0;

  /* set new tetra for travel */
  adja = &mesh->adja[4*(*adj-1)+1];
  if ( pt->v[ MMG5_ifar[*i][0] ] == *piv ) {
    *iface = MMG5_ifar[*i][0];
    *adj = adja[ *iface ] / 4;
    *piv = pt->v[ MMG5_ifar[*i][1] ];
  }
  else {
    assert(pt->v[ MMG5_ifar[*i][1] ] == *piv );
    *iface = MMG5_ifar[*i][1];
    *adj = adja[ *iface ] /4;
    *piv = pt->v[ MMG5_ifar[*i][0] ];
  }

  return 1;
}
