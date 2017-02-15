/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
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
 * \file mmg2d/libmmg2d_tools.c
 * \brief Tools functions for the mmg2d library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/

#include "mmg2d.h"

void MMG2D_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size == 3 ) {
    MMG2_length    = long_ani;
    MMG2_caltri    = caltri_ani;
    MMG2_caltri_in = caltri_ani_in;
    MMG2_buckin    = buckin_ani;
    MMG2_lissmet   = lissmet_ani;
    MMG2_optlen    = optlen_ani;
/*    interp     = interp_ani;
 */
  }
  else {
    MMG2_length     = long_iso;
    MMG2_caltri     = caltri_iso;
    MMG2_caltri_in  = caltri_iso_in;
    MMG2_buckin     = buckin_iso;
    MMG2_lissmet    = lissmet_iso;

    MMG2_optlen     = optlen_iso;
/*    interp     = interp_iso;
 */
  }

  return;
}

int MMG2D_Get_adjaTri(MMG5_pMesh mesh, int kel, int listri[3]) {

  if ( ! mesh->adja ) {
    if (! MMG2_hashel(mesh))
      return(0);
  }

  listri[0] = mesh->adja[3*(kel-1)+1]/3;
  listri[1] = mesh->adja[3*(kel-1)+2]/3;
  listri[2] = mesh->adja[3*(kel-1)+3]/3;

  return(1);
}

inline
int MMG2D_Get_adjaVertices(MMG5_pMesh mesh, int ip, int lispoi[MMG2D_LMAX])
{
  int start;

  if ( !mesh->tria ) return 0;

  start=MMG2_findTria(mesh,ip);
  if ( !start ) return 0;

  return MMG2D_Get_adjaVerticesFast(mesh,ip,start,lispoi);
}

inline
int MMG2D_Get_adjaVerticesFast(MMG5_pMesh mesh, int ip,int start, int lispoi[MMG2D_LMAX])
{
  MMG5_pTria pt;
  int k,prevk,nbpoi,iploc,i,i1,i2,*adja;

  pt   = &mesh->tria[start];

  for ( iploc=0; iploc<3; ++iploc ) {
    if ( pt->v[iploc] == ip ) break;
  }

  assert(iploc!=3);

  k = start;
  i = iploc;
  nbpoi = 0;
  do {
    if ( nbpoi == MMG2D_LMAX ) {
      fprintf(stdout,"  ## Warning: unable to compute adjacent vertices of the"
              " vertex %d:\nthe ball of point contain too many elements.\n",ip);
      return(0);
    }
    i1 = _MMG5_inxt2[i];
    lispoi[nbpoi] = mesh->tria[k].v[i1];
    ++nbpoi;

    adja = &mesh->adja[3*(k-1)+1];
    prevk = k;
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = _MMG5_inxt2[i];
  }
  while ( k && k != start );

  if ( k > 0 ) return(nbpoi);

  /* store the last point of the boundary triangle */
  if ( nbpoi == MMG2D_LMAX ) {
    fprintf(stdout,"  ## Warning: unable to compute adjacent vertices of the"
            " vertex %d:\nthe ball of point contain too many elements.\n",ip);
    return(0);
  }
  i1 = _MMG5_inxt2[i1];
  lispoi[nbpoi] = mesh->tria[prevk].v[i1];
  ++nbpoi;

  /* check if boundary hit */
  k = start;
  i = iploc;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i2 = _MMG5_iprv2[i];
    k  = adja[i2] / 3;
    if ( k == 0 )  break;

    if ( nbpoi == MMG2D_LMAX ) {
      fprintf(stdout,"  ## Warning: unable to compute adjacent vertices of the"
              " vertex %d:\nthe ball of point contain too many elements.\n",ip);
      return(0);
    }
    i  = adja[i2] % 3;
    lispoi[nbpoi] = mesh->tria[k].v[i];
    ++nbpoi;

    i  = _MMG5_iprv2[i];
  }
  while ( k );

  return nbpoi;
}

int MMG2D_Get_triFromEdge(MMG5_pMesh mesh, int ked, int *ktri, int *ied)
{
  int val;

  val = mesh->edge[ked].base;

  if ( !val ) return(0);

  *ktri = val/3;

  *ied = val%3;

  return 1;


}

void MMG2D_Free_triangles(MMG5_pMesh mesh) {

  if ( mesh->adja )
    _MMG5_DEL_MEM(mesh,mesh->adja,(3*mesh->ntmax+5)*sizeof(int));

  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->ntmax+1)*sizeof(MMG5_Tria));

  mesh->nt = 0;
  mesh->nti = 0;
  mesh->nenil = 0;

  return;
}

void MMG2D_Free_edges(MMG5_pMesh mesh) {

  if ( mesh->edge )
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->namax+1)*sizeof(MMG5_Edge));

  if ( mesh->xpoint )
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));

  mesh->na = 0;
  mesh->nai = 0;
  mesh->nanil = 0;

  mesh->xp = 0;

  return;
}

void MMG2D_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol) {

  /* sol */
  if ( sol && sol->m )
    _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));

  return;
}
