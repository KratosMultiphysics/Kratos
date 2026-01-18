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
 * \file mmg3d/optbdry_3d.c
 * \brief Functions for the optimization of very bad elements.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "mmg3dexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k   index of a tetra
 *
 * \return 1 if we move one of the vertices, 0 otherwise.
 *
 * Try to move the vertices of the tetra \a k to improve its quality.
 *
 */
int MMG3D_movetetrapoints(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree PROctree,MMG5_int k) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  /* double        *n; */
  int64_t       listv[MMG3D_LMAX+2];
  MMG5_int      ier/*,lists[MMG3D_LMAX+2]*/,base;
  int           i0,i,j/*,ilists*/,ilistv;
  int           /* improve,*/ internal,nm,/*maxit,*/ns;

  // improve = 1;
  internal = 1;
  nm = ns = 0;
  // maxit = 1;
  base = mesh->base;

  pt = &mesh->tetra[k];

  /* point j on face i */
  for (i=0; i<4; i++) {
    for (j=0; j<3; j++) {
      if ( pt->xt ) {
        pxt = &mesh->xtetra[pt->xt];
        if ( pxt->tag[MMG5_iarf[i][j]] & MG_REQ )  continue;
      }
      else  pxt = 0;
      i0  = MMG5_idir[i][j];
      ppt = &mesh->point[pt->v[i0]];
      if ( ppt->flag == base )  continue;
      else if ( MG_SIN(ppt->tag) )  continue;

      // if ( maxit != 1 )
      ppt->flag = base;

      ier = 0;
      if ( ppt->tag & MG_BDY ) {

        continue;

        /*        /\* Catch a boundary point by a boundary face *\/ */
/*         if ( !pt->xt || !(MG_BDY & pxt->ftag[i]) )  continue; */
/*         else if( ppt->tag & MG_NOM ){ */
/*           if( mesh->adja[4*(k-1)+1+i] ) continue; */
/*           ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,1); */
/*           if( !ier )  continue; */
/*           else if ( ier>0 ) */
/*             ier = MMG5_movbdynompt(mesh,met,PROctree,listv,ilistv,lists,ilists,improve); */
/*           else */
/*             return -1; */
/*         } */
/*         else if ( ppt->tag & MG_GEO ) { */
/*           ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0); */
/*           if ( !ier )  continue; */
/*           else if ( ier>0 ) */
/*             ier = MMG5_movbdyridpt(mesh,met,PROctree,listv,ilistv,lists,ilists,improve); */
/*           else */
/*             return -1; */
/*         } */
/*         else if ( ppt->tag & MG_REF ) { */
/*           ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0); */
/*           if ( !ier ) */
/*             continue; */
/*           else if ( ier>0 ) */
/*             ier = MMG5_movbdyrefpt(mesh,met,PROctree,listv,ilistv,lists,ilists,improve); */
/*           else */
/*             return -1; */
/*         } */
/*         else { */
/*           ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0); */
/*           if ( !ier ) */
/*             continue; */
/*           else if ( ier<0 ) */
/*             return -1; */

/*           n = &(mesh->xpoint[ppt->xp].n1[0]); */
/*           if ( !MG_GET(pxt->ori,i) ) { */
/*             if ( !MMG5_directsurfball(mesh,pt->v[i0],lists,ilists,n) ) */
/*               continue; */
/*           } */
/* //#warning CECILE a modifier pour opttyp */
/*           ier = MMG5_movbdyregpt(mesh,met, PROctree, listv,ilistv,lists,ilists,improve,improve); */
/*           if ( ier )  ns++; */
/*         } */
      }
      else if ( internal ) {
        ilistv = MMG5_boulevolp(mesh,k,i0,listv);
        if ( !ilistv )  continue;
        ier =  MMG3D_movnormal_iso(mesh,met,k,i0);
      }
      if ( ier ) {
        nm++;
        // if(maxit==1)
        ppt->flag = base;
      }
    }
  }
  if(nm)
    return 1;
  else
    return 0;

}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k   index of a tetra
 * \param i   index of point to delete in tetra \a k.
 *
 * \return 1 if success, 0 if we can't delete the point.
 *
 * Try to remove point i of tet k, try the three edges of k containing i.
 *
 */
int MMG3D_coledges(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int i) {
  MMG5_pTetra pt;
  double      len;
  int         ied,iedg,iq,i1,ilistcol;
  int64_t     listcol[MMG3D_LMAX+2];
  int         ier;
  int8_t      iface,ief;

  pt = &mesh->tetra[k];

  if ( MG_SIN(mesh->point[pt->v[i]].tag) ) {
    return 0;
  }

  /*3 possibilities to remove the vertex ib*/
  for(ied = 0 ; ied<3 ;ied++) {
    iedg  = MMG5_arpt[i][ied];
    len =  MMG5_lenedg(mesh,met,iedg,pt);

    if(len > 1.1) continue;
    iface = MMG5_ifar[iedg][0];
    ief   = MMG5_iarfinv[iface][iedg];
    iq    = MMG5_idir[iface][MMG5_iprv2[ief]];
    if(iq==i) {
      iface = MMG5_ifar[iedg][1];
      ief   = MMG5_iarfinv[iface][iedg];
      iq    = MMG5_idir[iface][MMG5_iprv2[ief]];
    }
    i1    = MMG5_idir[iface][MMG5_inxt2[ief]];

    assert( 0<=i1 && i1<4 && "unexpected local index for vertex");
    ilistcol = MMG5_boulevolp(mesh,k,i1,listcol);

    ilistcol = MMG5_chkcol_int(mesh,met,k,iface,ief,listcol,ilistcol,2);
    if ( ilistcol > 0 ) {
      ier = MMG5_colver(mesh,met,listcol,ilistcol,iq,2);
      if ( ilistcol < 0 ) continue;
      if ( ier < 0 ) return -1;
      else if(ier) {
        MMG3D_delPt(mesh,ier);
        return 1;
      }
    }
  }
  return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param k   index of a tetra
 * \param i   index of point to delete in tetra \a k.
 *
 * \return 1 if success, 0 if we can't delete the point.
 *
 * Try to delete the point i of tet k. Try all the edges containing i.
 *
 */
int MMG3D_deletePoint(MMG5_pMesh mesh,  MMG5_pSol met,MMG3D_pPROctree PROctree,
                       MMG5_int k,int i) {
  MMG5_pTetra pt;
  int         il,ilist,ip;
  int64_t     list[MMG3D_LMAX+2];
  MMG5_int    iel;

  pt = &mesh->tetra[k];

  if ( MG_SIN(mesh->point[pt->v[i]].tag) ) {
    return 0;
  }

  assert( 0<=i && i<4 && "unexpected local index for vertex");
  ilist = MMG5_boulevolp(mesh,k,i,list);
  if (ilist > 30 ) return 0;

  for(il = 0 ; il<ilist ; il++) {
    iel = list[il] / 4;
    ip  = list[il] % 4;
    if( MMG3D_coledges(mesh,met,iel,ip) ) {
      return 1;
    }
  }

  return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param PROctree pointer to the PROctree structure.
 * \param k   index of a tetra
 *
 * \return 1 if success, 0 if fail.
 *
 * Try to optimize the tetra k. This tetra has a face on the boundary.
 *
 */
int MMG3D_optbdry(MMG5_pMesh mesh,MMG5_pSol met,MMG3D_pPROctree PROctree,MMG5_int k) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int          ib,i,j;
  int64_t      list[MMG3D_LMAX+2];
  MMG5_int     ipb,it1,it2;
  int          iedg,ier,ilist,ied,ia,ret,imove;

  imove = 0;

  pt = &mesh->tetra[k];
  assert(pt->xt);

  pxt = &mesh->xtetra[pt->xt];

  for(i=0 ; i<4 ; i++)
    if ( pxt->ftag[i] & MG_BDY ) break;

  assert ( i< 4 );
  if ( i== 4 ) return 0;

  ib  = i;
  ipb = pt->v[ib];

  /*check that the vertex is not a boundary one*/
  if ( mesh->point[ipb].tag & MG_BDY ) return 0;

  /* try to move the vertex in order to improve the quality*/
  ier = 0;
  if ( !mesh->info.nomove ) {
    for(j = 0 ; j<3 ; j++) {
      imove = MMG3D_movetetrapoints(mesh,met,PROctree,k);
      ier += imove;
      if(!imove) break;
    }
    if(ier) {
      imove = 1;
    }
  }

  if(!mesh->info.noinsert) {
    /*try to remove the non-bdry vertex*/
    ier = MMG3D_coledges(mesh,met,k,ib);
    if(ier) return 1;

    /* try to remove the non-bdry vertex : with all the edges containing the
     * vertex */
    ier = MMG3D_deletePoint(mesh,met,PROctree,k,i);
    if(ier) return 1;
  }


  /*try to swap the 3 internal edges*/
  if(!mesh->info.noswap) {
    for(ied = 0 ; ied<3 ;ied++) {
      iedg  = MMG5_arpt[i][ied];
      ier = MMG3D_swpItem(mesh,met,PROctree,k,iedg);
      if(ier) {
        return 1;
      }
    }
    /*try to swap the bdry edges*/
    for (j=0; j<3; j++) {
      ia  = MMG5_iarf[i][j];

      /* Mark the edge as boundary in case that the tag is missing */
      pxt->tag[ia] |= MG_BDY;

      /* No swap of geometric edge */
      if ( MG_EDG_OR_NOM(pxt->tag[ia]) || (pxt->tag[ia] & MG_REQ) ) {
        continue;
      }

      ret = MMG5_coquilface(mesh,k,i,ia,list,&it1,&it2,0);
      ilist = ret / 2;
      if ( ret < 0 )  return -1;
      /* CAUTION: trigger collapse with 2 elements */
      if ( ilist <= 1 )  continue;

      /* Here, we work on a boundary edge lying along a boundary face */
      ier = MMG5_chkswpbdy(mesh,met,list,ilist,it1,it2,2);
      if ( ier <  0 )
        return -1;
      else if ( ier ) {
        ier = MMG5_swpbdy(mesh,met,list,ret,it1,PROctree,2);
        if ( ier < 0 )  return -1;
        else if(ier) {
          return 1;
        }
      }
    }
  }


  return imove;

}
