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
 * \file mmg3d/chkmsh_3d.c
 * \brief Check the input mesh validity.
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

#define  MMG5_EPSLOC   1.00005
#define  IEDG(a,b) (((a) > 0) && ((b) > 0)) ? ((a)+(b)) : (((a)+(b))-(1))

extern int8_t ddb;

/**
 * \param mesh pointer toward mesh
 *
 * Test that tetra have positive volumes.
 *
 * \warning Not used.
 */
void MMG5_chkvol(MMG5_pMesh mesh) {
  MMG5_pTetra    pt;
  MMG5_int       k;
#ifdef DEBUG
  int            ier=1;
#endif

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( MMG5_orvol(mesh->point,pt->v) < MMG5_NULKAL ) {
      fprintf(stderr,"\n  ## Warning: %s: tetra %" MMG5_PRId " volume %e\n",__func__,
             k,MMG5_orvol(mesh->point,pt->v));
#ifdef DEBUG
      ier = 0;
#endif
    }
  }
#ifdef DEBUG
  assert(ier);
#endif
}

/**
 * \param mesh pointer toward the mesh
 * \param start tetra from which we start to travel
 * \param na edge vertex
 * \param nb edge vertex
 * \param tag edge tag
 * \param ref edge ref
 * \param piv global index of the pivot to set the sense of travel
 * \param adj index of adjacent tetra for the travel
 *
 * \return -1 if fail, \a start if shell has been completely travelled, 0
 * otherwise (open shell).
 *
 * Test consistency of tag and ref of the edge \a na \a nb from tetra \a start
 * by traveling its shell in one direction (given by the pivot \a piv).
 *
 */
static inline
int MMG3D_chk_shellEdgeTag_oneDir(MMG5_pMesh  mesh,MMG5_int start, MMG5_int na, MMG5_int nb,
                                  int16_t tag,MMG5_int ref, MMG5_int piv,MMG5_int adj) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     *adja;
  int8_t       i;

  while ( adj && (adj != start) ) {
    pt     = &mesh->tetra[adj];

    /* identification of edge number in tetra adj */
    if ( !MMG3D_findEdge(mesh,pt,adj,na,nb,1,NULL,&i) ) return -1;

    /* update edge ref and tag */
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      /* Test only BDY edges */
      if ( pxt->tag[i] & MG_BDY ) {
        assert (pxt->tag[i] == tag && "non consistent tags");
        assert (pxt->edg[i] == ref && "non consistent refs");
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
 * \param mesh pointer toward the mesh
 * \param start tetra from which we start to travel
 * \param ia local index of edge that must be updated
 * \param tag edge tag
 * \param ref edge ref
 * \return 1 if success, 0 if fail.
 *
 * Test consistency of tag and ref of the boundary edge \ia of tetra \a start by
 * traveling its shell.
 *
 */
int MMG3D_chk_shellEdgeTag(MMG5_pMesh  mesh,MMG5_int start, int8_t ia,int16_t tag,MMG5_int ref) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_int     piv,na,nb,adj,*adja;

  pt   = &mesh->tetra[start];

  assert( start >= 1 &&  MG_EOK(pt) && "invalid tetra" );
  assert ( tag & MG_BDY && "Unexpected non boundary tag");

  pxt  = NULL;
  na   = pt->v[MMG5_iare[ia][0]];
  nb   = pt->v[MMG5_iare[ia][1]];

  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    /* Test only BDY edges */
    if ( pxt->tag[ia] & MG_BDY ) {
      assert (pxt->tag[ia] == tag && "non consistent tags"); ;
      assert (pxt->edg[ia] == ref && "non consistent refs"); ;
    }
  }

  /* Travel in one direction */
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][0]] / 4;
  piv = pt->v[MMG5_ifar[ia][1]];

  adj = MMG3D_chk_shellEdgeTag_oneDir(mesh,start,na,nb,tag,ref,piv,adj);

  /* If all shell has been travelled, stop, else, travel it the other sense */
  if ( adj > 0 ) {
    assert ( adj == start );
    return 1;
  }
  else if ( adj < 0 ) return 0;

  assert(!adj);

  pt = &mesh->tetra[start];
  adja = &mesh->adja[4*(start-1)+1];
  adj = adja[MMG5_ifar[ia][1]] / 4;
  piv = pt->v[MMG5_ifar[ia][0]];

  adj = MMG3D_chk_shellEdgeTag_oneDir(mesh,start,na,nb,tag,ref,piv,adj);

  if ( adj < 0 ) return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * Test consistency between the tags in the xtetra of all mesh edges marked as
 * boundaries.
 *
 * \warning Not used.
 */
void MMG3D_chkmeshedgestags(MMG5_pMesh mesh) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_Hash      hash;
  int            i,tag;
  MMG5_int       k,nt,ip1,ip2;

  /* Rough eval of the number of boundary triangles */
  nt = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] & MG_BDY ) {
        ++nt;
      }
    }
  }
  nt = nt/2 + 1;

  /* Travel mesh edges and hash boundary ones */
  MMG5_hashNew(mesh,&hash,nt,3*nt);

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<6; i++) {
      if ( pxt->tag[i] & MG_BDY ) {
        ip1 = pt->v[MMG5_iare[i][0]];
        ip2 = pt->v[MMG5_iare[i][1]];
        tag = MMG5_hashEdgeTag ( mesh,&hash,ip1,ip2,pxt->tag[i]);
        if ( tag != pxt->tag[i] ) {
          fprintf(stderr,"Error: %s: %d: Non consistency at tet %" MMG5_PRId " (%" MMG5_PRId "), edge %d:%" MMG5_PRId "--%" MMG5_PRId "\n ",
                  __func__,__LINE__,k,MMG3D_indElt(mesh,k),i,ip1,ip2);
          assert( tag == pxt->tag[i] && "edge tag error" );
        }
      }
    }
  }
  MMG5_DEL_MEM(mesh,hash.item);
}


/**
 * \param mesh pointer toward the mesh
 * \param ip1 first vertex of edge to test
 * \param ip2 second vertex of edge to test
 * \param tag edge tag
 *
 * Test consistency between the tags of the edge \a ip1 - \a ip2 from all the
 * tetra of the edge shell.
 *
 * \warning Not used.
 */
void MMG3D_chkedgetag(MMG5_pMesh mesh, MMG5_int ip1, MMG5_int ip2, int tag) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_int       k,i1,i2;
  int            i;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (i=0; i<6; i++) {
      i1 = pt->v[MMG5_iare[i][0]];
      i2 = pt->v[MMG5_iare[i][1]];

      if ( ((i1==ip1) && (i2==ip2)) || ((i2==ip1) && (i1==ip2)) ) {
        if ( pxt->tag[i] != tag ) {
          fprintf(stderr,"Error: %s: %d: Non consistency at tet %" MMG5_PRId " (%" MMG5_PRId "), edge %d\n ",
                  __func__,__LINE__,k,MMG3D_indElt(mesh,k),i);
          assert(0);
        }
      }
    }
  }
}

/**
 * \param mesh pointer toward the mesh
 * \param ppt pointer toward unconsistent point
 * \param k tetra index
 * \param i local index of edge in tetra \a k
 * \param ip1 first vertex of edge to test
 * \param ip2 second vertex of edge to test
 *
 * Print error message when an edge vertex has a non-consistent tag.
 *
 */
static inline
void MMG3D_consistency_error_message(MMG5_pMesh mesh,MMG5_pPoint ppt,MMG5_int k,int i,MMG5_int ip1,MMG5_int ip2) {

  assert ( mesh->tetra && "no tetra array");
  MMG5_pTetra   pt = &mesh->tetra[k];

  assert ( pt->xt && "no xtetra");
  MMG5_pxTetra pxt = &mesh->xtetra[pt->xt];

  fprintf(stderr,"Error: %s: %d: Tag error at point %" MMG5_PRId " (%" MMG5_PRId "), "
          "tetra %" MMG5_PRId " (%" MMG5_PRId "), edge %d:%" MMG5_PRId "--%" MMG5_PRId " (%" MMG5_PRId
          "--%" MMG5_PRId ").\n",__func__,__LINE__,
          ip1,MMG3D_indPt(mesh,ip1),k,MMG3D_indElt(mesh,k),i,ip1,ip2,
          MMG3D_indPt(mesh,ip1),MMG3D_indPt(mesh,ip2));
  fprintf(stderr," point tag: %d; edge tag: %d\n",ppt->tag,pxt->tag[i]);

  /** An error has been detected: check the consistency between the tags of
   * tetra edges */
  MMG3D_chkedgetag(mesh,ip1,ip2,pxt->tag[i]);
  assert(0);
}

/**
 * \param mesh
 *
 * Test consistency between points and edges tags. If an error is detected,
 * hash mesh edges to check the consistency between the tags of tetra edges.
 *
 * \warning Not used.
 */
void MMG3D_chkpointtag(MMG5_pMesh mesh) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    p1,p2;
  int            i,i1,i2;
  MMG5_int       k,ip1,ip2;

  /** Check consistency between edge tags and point tags */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];

    for ( i=0; i<6; ++i ) {
      i1  = MMG5_iare[i][0];
      i2  = MMG5_iare[i][1];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      p1 = &mesh->point[ip1];
      p2 = &mesh->point[ip2];

      if ( MG_EDG(pxt->tag[i]) ) {
        if ( !(MG_EDG(p1->tag) || MG_SIN(p1->tag)) ) {
          MMG3D_consistency_error_message(mesh,p1,k,i,ip1,ip2);
        }
        if ( !(MG_EDG(p2->tag) || MG_SIN(p2->tag)) ) {
          MMG3D_consistency_error_message(mesh,p2,k,i,ip1,ip2);
        }
      }

      if ( pxt->tag[i] & MG_NOM ) {
        if ( !MG_SIN_OR_NOM(p1->tag) ) {
          MMG3D_consistency_error_message(mesh,p1,k,i,ip1,ip2);
        }
        if ( !MG_SIN_OR_NOM(p2->tag) ) {
          MMG3D_consistency_error_message(mesh,p2,k,i,ip1,ip2);
        }
      }
    }
  }
}

/**
* \return 0 if fail, 1 otherwise
 *
 * \warning Not used.
 */
int MMG5_chkmshsurf(MMG5_pMesh mesh){
  MMG5_pTria pt;
  MMG5_int   k,k1;
  MMG5_int   *adja,*adja1;
  int8_t     i,voy;

  for (k=1; k<=mesh->nt; k++) {
    pt   = &mesh->tria[k];
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_NOM )  continue;
      k1  = adja[i] / 3;
      voy = adja[i] % 3;

      if(!k1) continue;
      adja1 = &mesh->adjt[3*(k1-1)+1];

      if(adja1[voy] / 3 != k){
        fprintf(stderr,"\n  ## Warning: %s: wrong adjacency relation"
                " for triangles : %" MMG5_PRId " %" MMG5_PRId " \n",__func__,k,k1);
        return 0;
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 otherwise
 *
 * Check the number of boundary faces in each edge shell and the consistency of the edge tag.
 */
static inline
int  MMG3D_chkcoquilface(MMG5_pMesh mesh) {
  MMG5_pTetra pt;
  MMG5_pxTetra pxt;
  MMG5_int k,it1,it2;
  int64_t list[MMG3D_LMAX+2];
  int i,j,ret;
  int8_t ia;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( (!MG_EOK(pt)) || pt->ref < 0 || (pt->tag & MG_REQ) )   continue;
    else if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for (i=0; i<4; i++) {
      if ( !(pxt->ftag[i] & MG_BDY) ) continue;
      for (j=0; j<3; j++) {
        ia  = MMG5_iarf[i][j];

        /* No check for geom edge (I am not sure that it can works) */
        if ( MG_EDG_OR_NOM(pxt->tag[ia]) || (pxt->tag[ia] & MG_REQ) )
        {
          continue;
        }

        ret = MMG5_coquilface(mesh,k,i,ia,list,&it1,&it2,0);
        if ( ret < 0 )  return 0;
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param severe level of performed check (unused)
 * \param base unused argument.
 * \return 0 if fail, 1 if success.
 *
 * Check the mesh validity
 *
 */
int MMG5_mmg3dChkmsh(MMG5_pMesh mesh,int severe,MMG5_int base) {
  MMG5_pTetra     pt,pt1,pt2;
  MMG5_pxTetra    pxt;
  MMG5_int        *adja,*adja1,adj,adj1,k,iadr,iel;
  MMG5_int        a0,a1,a2,b0,b1,b2;
  int             i;
  uint8_t         voy,voy1;
  static int8_t   mmgErr0=0,mmgErr1=0,mmgErr2=0,mmgErr3=0,mmgErr4=0,mmgErr5=0;

  /* Check edge tag consistency (between xtetra) */
  MMG3D_chkmeshedgestags(mesh);

  /* Check point tags consistency with edge tag */
  MMG3D_chkpointtag(mesh);

  if ( !mesh->adja ) return 1;

  /* Check edge tag consistency with number of boundary faces in the edge shell */
  MMG3D_chkcoquilface(mesh);

  for (k=1; k<=mesh->ne; k++) {
    pt1 = &mesh->tetra[k];
    if ( !MG_EOK(pt1) || pt1->ref < 0 )   continue;
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<4; i++) {
      adj = adja[i];

      if ( !adj )  continue;
      adj /= 4;
      voy = adja[i] % 4;

      if ( adj == k ) {
        if ( !mmgErr0 ) {
          fprintf(stderr,"\n  ## Error: %s: 1. at least 1 wrong adjacency %" MMG5_PRId " %" MMG5_PRId "\n",
                  __func__,MMG3D_indElt(mesh,k),MMG3D_indElt(mesh,adj));
          fprintf(stderr,"triangle %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,k),
                  MMG3D_indPt(mesh,pt1->v[0]),MMG3D_indPt(mesh,pt1->v[1]),
                  MMG3D_indPt(mesh,pt1->v[2]),MMG3D_indPt(mesh,pt1->v[3]));
          fprintf(stderr,"adj (%" MMG5_PRId "): %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,k),
                  MMG3D_indElt(mesh,adja[0]/4),MMG3D_indElt(mesh,adja[1]/4),
                  MMG3D_indElt(mesh,adja[2]/4),MMG3D_indElt(mesh,adja[3]/4));
          mmgErr0 = 1;
        }
        return 0;
      }
      pt2 = &mesh->tetra[adj];
      if ( !MG_EOK(pt2) || pt2->ref < 0 ){
        if ( !mmgErr1 ) {
          fprintf(stderr,"\n  ## Error: %s: 4. at least 1 invalid adjacent %" MMG5_PRId " %" MMG5_PRId "\n",
                  __func__,MMG3D_indElt(mesh,adj),MMG3D_indElt(mesh,k));
          fprintf(stderr,"vertices of k   %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,k),
                  MMG3D_indPt(mesh,pt1->v[0]),MMG3D_indPt(mesh,pt1->v[1]),
                  MMG3D_indPt(mesh,pt1->v[2]),MMG3D_indPt(mesh,pt1->v[3]));
          fprintf(stderr,"vertices of adj %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,adj),
                  MMG3D_indPt(mesh,pt2->v[0]),MMG3D_indPt(mesh,pt2->v[1]),
                  MMG3D_indPt(mesh,pt2->v[2]),MMG3D_indPt(mesh,pt2->v[3]));
          fprintf(stderr,"adj(%" MMG5_PRId "): %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,k),
                  MMG3D_indElt(mesh,adja[0]/4),MMG3D_indElt(mesh,adja[1]/4),
                  MMG3D_indElt(mesh,adja[2]/4),MMG3D_indElt(mesh,adja[3]/4));
          mmgErr1 = 1;
        }
        return 0;
      }
      iadr  = (adj-1)*4 + 1;
      adja1 = &mesh->adja[iadr];
      adj1  = adja1[voy] / 4;
      voy1  = adja1[voy] % 4;
      if ( adj1 != k || voy1 != i ) {
        if ( !mmgErr2 ) {
          fprintf(stderr,"\n  ## Error: %s: 2. at least 1 wrong adjacency %" MMG5_PRId " %" MMG5_PRId "\n",
                  __func__,MMG3D_indElt(mesh,k),MMG3D_indElt(mesh,adj1));
          fprintf(stderr,"vertices of %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,k),
                  MMG3D_indPt(mesh,pt1->v[0]),MMG3D_indPt(mesh,pt1->v[1]),
                  MMG3D_indPt(mesh,pt1->v[2]),MMG3D_indPt(mesh,pt1->v[3]));
          fprintf(stderr,"vertices of adj %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,adj),
                  MMG3D_indPt(mesh,pt2->v[0]),MMG3D_indPt(mesh,pt2->v[1]),
                  MMG3D_indPt(mesh,pt2->v[2]),MMG3D_indPt(mesh,pt2->v[3]));
          fprintf(stderr,"adj(%" MMG5_PRId "): %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,k),
                  MMG3D_indElt(mesh,adja[0]/4),MMG3D_indElt(mesh,adja[1]/4),
                  MMG3D_indElt(mesh,adja[2]/4),MMG3D_indElt(mesh,adja[3]/4));
          fprintf(stderr,"adj(%" MMG5_PRId "): %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG3D_indElt(mesh,adj),
                  MMG3D_indElt(mesh,adja1[0]/4),MMG3D_indElt(mesh,adja1[1]/4),
                  MMG3D_indElt(mesh,adja1[2]/4),MMG3D_indElt(mesh,adja1[3]/4));
          mmgErr2 = 1;
        }
        return 0;
      }

      a0 = pt1->v[MMG5_idir[i][0]];
      a1 = pt1->v[MMG5_idir[i][1]];
      a2 = pt1->v[MMG5_idir[i][2]];

      b0 = pt2->v[MMG5_idir[voy][0]];
      b1 = pt2->v[MMG5_idir[voy][1]];
      b2 = pt2->v[MMG5_idir[voy][2]];

      if(!(((a0 == b0)&&(a1 == b1)&&(a2 ==b2))||((a0 == b0)&&(a1 == b2)&&(a2 ==b1))\
           || ((a0 == b1)&&(a1 == b0)&&(a2 ==b2)) || ((a0 == b1)&&(a1 == b2)&&(a2 ==b0)) \
           || ((a0 == b2)&&(a1 == b0)&&(a2 ==b1)) || ((a0 == b2)&&(a1 == b1)&&(a2 ==b0)) )){
        if ( !mmgErr3 ) {
          fprintf(stderr,"\n  ## Warning: %s: Inconsistent faces : tetra %" MMG5_PRId " face %d;"
                  " tetra %" MMG5_PRId " face %i \n",__func__,MMG3D_indElt(mesh,k),i,
                  MMG3D_indElt(mesh,adj),voy);
          fprintf(stderr,"Tet 1 : %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",MMG3D_indPt(mesh,a0),
                  MMG3D_indPt(mesh,a1),MMG3D_indPt(mesh,a2));
          fprintf(stderr,"Tet 2 : %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",MMG3D_indPt(mesh,b0),
                  MMG3D_indPt(mesh,b1),MMG3D_indPt(mesh,b2));
          mmgErr3 = 1;
        }
        return 0;
      }
    }
  }

  /* This test may have to be disactivated : check wether boundary faces (i.e. no neighbour)
     arise only when a BDY face is hit */
  for(k=1;k<=mesh->ne;k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;
    adja = &mesh->adja[4*(k-1)+1];
    for(i=0;i<4;i++){
      if(!adja[i]){
        if(!pt->xt){
          if ( !mmgErr4 ) {
            mmgErr4 = 1;
            fprintf(stderr,"\n  ## Error: %s: Tetra %" MMG5_PRId ": boundary face"
                    " not tagged: %d \n",__func__,MMG3D_indElt(mesh,k),i);
          }
          return 0;
        }
        else{
          pxt = &mesh->xtetra[pt->xt];
          if(!(pxt->ftag[i] & MG_BDY)){
            if ( !mmgErr4 ) {
              mmgErr4 = 1;
              fprintf(stderr,"\n  ## Error: %s: Tetra %" MMG5_PRId ": boundary face"
                      " not tagged : %d \n",__func__,MMG3D_indElt(mesh,k),i);
            }
            return 0;
          }
        }
      }
    }
  }

  /* Case of an implicit surface embedded in mesh : check whether a face separating two
     different subdomains is tagged bdy */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || pt->ref < 0 )   continue;

    adja = &mesh->adja[4*(k-1)+1];
    for(i=0; i<4; i++){
      if(!adja[i]) continue;
      iel = adja[i] / 4;
      pt1 = &mesh->tetra[iel];

      if(pt->ref != pt1->ref){
        if(!pt->xt){
          if ( !mmgErr5 ) {
            mmgErr5 = 1;
            fprintf(stderr,"\n  ## Error: %s: Tetra %" MMG5_PRId " face %d: common"
                    " face is a limit of two subdomains"
                    " and has not xt : %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "  \n",__func__,
                    MMG3D_indElt(mesh,k),i,
                    MMG3D_indPt(mesh,pt->v[MMG5_idir[i][0]]),
                    MMG3D_indPt(mesh,pt->v[MMG5_idir[i][1]]),
                    MMG3D_indPt(mesh,pt->v[MMG5_idir[i][2]]));
          }
          return 0;
        }
        else{
          pxt = &mesh->xtetra[pt->xt];
          if(!(pxt->ftag[i] & MG_BDY)){
            if ( !mmgErr5 ) {
              mmgErr5 = 1;
              fprintf(stderr,"\n  ## Error: %s: Tetra %" MMG5_PRId " %d : common"
                      " face is a limit of two subdomains"
                      " and is not tagged %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " -->%d\n",__func__,
                      MMG3D_indElt(mesh,k),i,
                      MMG3D_indElt(mesh,pt->v[MMG5_idir[i][0]]),
                      MMG3D_indPt(mesh,pt->v[MMG5_idir[i][1]]),
                      MMG3D_indPt(mesh,pt->v[MMG5_idir[i][2]]), pxt->ftag[i]);
            }
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

/**
 * Search boundary faces containing point np.
 *
 * \return 0 if fail, 1 otherwise
 *
 * \warning Not used.
 **/
int MMG5_chkptonbdy(MMG5_pMesh mesh,MMG5_int np){
  MMG5_pTetra      pt;
  MMG5_pxTetra     pxt;
  MMG5_pPoint      p0;
  MMG5_int         k;
  int8_t           i,j,ip;
  static int8_t    mmgWarn0=0,mmgWarn1=0;

  for(k=1;k<=mesh->np;k++)
    mesh->point[k].flag = 0;

  /* Put flag = 1 at each point belonging to a boundary face */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if(!MG_EOK(pt)) continue;
    if(!pt->xt) continue;
    pxt = &mesh->xtetra[pt->xt];
    for(i=0; i<4; i++){
      if(!(pxt->ftag[i] & MG_BDY)) continue;
      for(j=0; j<3; j++){
        ip = MMG5_idir[i][j];
        if(pt->v[ip] == np) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: point %" MMG5_PRId " on face %d of tetra %" MMG5_PRId " :"
                   " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",__func__, MMG3D_indPt(mesh,pt->v[ip]),i,
                   MMG3D_indElt(mesh,k), MMG3D_indPt(mesh,pt->v[0]),
                   MMG3D_indPt(mesh,pt->v[1]),
                   MMG3D_indPt(mesh,pt->v[2]), MMG3D_indPt(mesh,pt->v[3]));
          }
        }
        p0 = &mesh->point[pt->v[ip]];
        p0->flag = 1;
      }
    }
  }


  /* Make sure that all the remaining points are not tagged BDY */
  for(k=1; k<=mesh->np; k++){
    p0 = &mesh->point[k];
    if(!MG_VOK(p0)) continue;
    if(p0->flag) continue;
    if(p0->tag & MG_BDY){
      if ( !mmgWarn1 ) {
        mmgWarn1 = 1;
        fprintf(stderr,"\n  ## Error: %s: point %" MMG5_PRId " tagged bdy while belonging to no BDY face\n",
                __func__,MMG3D_indPt(mesh,k));
      }
      return 0;
    }
  }

  return 1;
}

/**
 *
 * Count how many boundary faces share point nump.
 *
 * \warning Not used.
 */
int MMG5_cntbdypt(MMG5_pMesh mesh, MMG5_int nump){
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_int      k,v0,v1,v2;
  int           nf;
  int8_t        i,j,ip;
  static int8_t mmgWarn0 = 0;

  nf = 0;

  for(k=1; k<=mesh->ne;k++){
    pt = &mesh->tetra[k];
    if(!MG_EOK(pt)) continue;
    if(!pt->xt) continue;
    pxt = &mesh->xtetra[pt->xt];
    for(i=0; i<4; i++){
      if(!(pxt->ftag[i] & MG_BDY)) continue;
      for(j=0; j<3; j++){
        ip = MMG5_idir[i][j];
        if(pt->v[ip] == nump){
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            v0 = pt->v[MMG5_idir[i][0]];
            v1 = pt->v[MMG5_idir[i][0]];
            v2 = pt->v[MMG5_idir[i][0]];

            fprintf(stderr,"\n  ## Error: %s: face %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " in tetra : %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",
                   __func__,MMG3D_indPt(mesh,v0),MMG3D_indPt(mesh,v1),
                   MMG3D_indPt(mesh,v2),
                   MMG3D_indPt(mesh,pt->v[0]),MMG3D_indPt(mesh,pt->v[1]),
                   MMG3D_indPt(mesh,pt->v[2]),MMG3D_indPt(mesh,pt->v[3]));
          }
          nf++;
        }
      }
    }
  }
  return nf;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if fail, 1 otherwise.
 *
 * Count the number of tetras that have several boundary faces, as well as the
 * number of internal edges connecting points of the boundary.
 *
 */
int MMG5_chkfemtopo(MMG5_pMesh mesh) {
  MMG5_pTetra      pt,pt1;
  MMG5_pxTetra     pxt;
  MMG5_pPoint      p0,p1;
  MMG5_int         k,ntet,ned,np,np1,ischk,npchk,iel;
  int64_t          list[MMG3D_LMAX+2];
  int              nf,ilist,l;
  int8_t           i0,j,i,i1,ia,ier;

  ntet = ned = 0;
  for(k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Count elements with at least two boundary faces */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    else if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    nf = 0;
    for (i=0; i<4; i++) {
      if ( pxt->ftag[i] & MG_BDY )  nf++;
    }
    if ( nf >= 2 )  ntet++;
  }
  if ( mesh->info.imprim > 0 && ntet )
    printf("  *** %" MMG5_PRId " tetras with at least 2 boundary faces.\n",ntet);

  /* Count internal edges connecting two points of the boundary */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<4; i++) {
      np = pt->v[i];
      p0 = &mesh->point[np];
      if ( !(p0->tag & MG_BDY) )  continue;

      ischk = p0->flag % 2;
      if ( ischk )  continue;
      p0->flag += 1;

      ilist = MMG5_boulevolp(mesh,k,i,list);
      for (l=0; l<ilist; l++) {
        iel = list[l] / 4;
        i0  = list[l] % 4;
        i1  = i0;

        pt1 = &mesh->tetra[iel];
        for (j=0; j<3; j++) {
          i1  = MMG5_inxt3[i1];
          np1 = pt1->v[i1];
          if ( np1 < np )  continue;
          p1 = &mesh->point[np1];
          if ( !(p1->tag & MG_BDY) )  continue;

          ischk = p1->flag % 2;
          npchk = p1->flag / 2;
          if ( npchk == np )  continue;

          ia = IEDG(i0,i1);
          p1->flag = 2*np + ischk;

          ier = MMG5_srcbdy(mesh,iel,ia);
          if ( ier<0 ) return 0;
          else if ( !ier )  ned++;
        }
      }
    }
  }
  if ( mesh->info.imprim > 0 && ned )
    printf("  *** %" MMG5_PRId " internal edges connecting boundary points.\n",ned);
  return 1;
}

/**
 *
 * Search face n0,n1,n2 in mesh, and get the support tetras, with the
 * corresponding refs.
 *
 * \warning Not used.
 */
int srcface(MMG5_pMesh mesh,MMG5_int n0,MMG5_int n1,MMG5_int n2) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_int      ref,minn,maxn,sn,k,ip0,ip1,ip2,mins,maxs,sum;
  int16_t       tag;
  int8_t        i;
  static int8_t mmgWarn0 = 0;

  minn = MG_MIN(n0,MG_MIN(n1,n2));
  maxn = MG_MAX(n0,MG_MAX(n1,n2));
  sn   = n0 + n1 + n2;
  pxt = 0;

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;

    if( pt->xt ) pxt = &mesh->xtetra[pt->xt];
    for(i=0; i<4; i++) {
      ip0 = pt->v[MMG5_idir[i][0]];
      ip1 = pt->v[MMG5_idir[i][1]];
      ip2 = pt->v[MMG5_idir[i][2]];

      mins = MG_MIN(ip0,MG_MIN(ip1,ip2));
      maxs = MG_MAX(ip0,MG_MAX(ip1,ip2));
      sum  = ip0 + ip1 + ip2;
      tag  = pt->xt ? pxt->ftag[i] : 0;
      ref  = pt->xt ? pxt->ref[i] : 0;

      if( mins == minn && maxs == maxn && sum == sn ) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  ## Error: %s: Face %d in tetra %" MMG5_PRId " with ref %" MMG5_PRId ":"
                  " corresponding ref %" MMG5_PRId " , tag: %d\n",__func__,i,
                  MMG3D_indElt(mesh,k),pt->ref,ref,tag);
        }
      }
    }
  }


  return 1;
}
