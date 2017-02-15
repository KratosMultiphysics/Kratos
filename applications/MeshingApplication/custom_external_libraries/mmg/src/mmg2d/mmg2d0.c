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
 * \file mmg2d/mmg2d0.c
 * \brief Mesh optimization functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"


int MMG2_mmg2d0(MMG5_pMesh mesh,MMG5_pSol sol) {
  int       ns,nm,nsiter,nmiter,nmbar,it,maxtou;
  double    declic;

  /*optim*/
  ns     = 0;
  nm     = 0;
  nmiter = 0;
  nsiter = 0;
  it     = 0;
  maxtou = 100;
  do {
    /*edge flip*/
    if(!mesh->info.noswap) {
      declic = 1.5 / ALPHA;
      nsiter = MMG2_cendel(mesh,sol,declic,-1);
      if ( nsiter && mesh->info.imprim < 0)
        fprintf(stdout,"     %7d SWAPPED\n",nsiter);

      ns+=nsiter;
    }
    /*point relocation*/
    if(!mesh->info.nomove) {
      declic = 1.5 / ALPHA;
      nmiter = MMG2_optlen(mesh,sol,declic,-1);
      if(sol->size==1) nmbar =  optlen_iso_bar(mesh,sol,declic,-1);
      else nmbar=0;
      nm += nmiter+nmbar;
      if ( mesh->info.imprim < 0)
        fprintf(stdout,"     %7d + %7d MOVED \n",nmiter,nmbar);
    }


  } while((nmiter+nsiter > 0) && (++it <= maxtou));
  if ( mesh->info.imprim )
    fprintf(stdout,"     %7d SWAPPED %7d MOVED\n",ns,nm);

  return(1);
}
