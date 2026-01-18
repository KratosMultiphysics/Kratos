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
 * \file mmg2d/tools_2d.c
 * \brief  Various tools.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg2d_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param kel index of the element in the unpacked mesh
 *
 *
 * \return 0 if fail, index of the element in packed numerotation otherwise.
 *
 * find the element index in packed numerotation
 *
 */
MMG5_int MMG2D_indElt(MMG5_pMesh mesh, MMG5_int kel) {
    MMG5_pTria      pt;
    MMG5_int        ne, k;

    ne = 0;
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( MG_EOK(pt) ) {
            ne++;
            if ( k == kel )  return ne;
        }
    }
    return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param kp index of the point in the unpacked mesh
 *
 * \return 0 if fail, index of the point in packed numerotation otherwise.
 *
 * find the point index in packed numerotation
 *
 */
MMG5_int MMG2D_indPt(MMG5_pMesh mesh, MMG5_int kp) {
    MMG5_pPoint ppt;
    MMG5_int    np, k;

    np = 0;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( MG_VOK(ppt) ) {
            np++;
            if ( k == kp )  return np;
        }
    }
    return 0;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param nsd index of subdomain to keep.
 *
 * Keep only subdomain of index \a nsd and remove other subdomains.
 *
 */
void MMG2D_keep_only1Subdomain ( MMG5_pMesh mesh,int nsd ) {

  if ( !nsd ) {
    return;
  }

  if ( mesh->info.imprim > 4 || mesh->info.ddebug ) {
    fprintf(stdout,"\n  -- ONLY KEEP DOMAIN OF REF %d\n",nsd );
  }

  MMG5_mark_verticesAsUnused ( mesh );

  MMG5_keep_subdomainElts ( mesh, nsd, MMG2D_delElt );

  MMG5_mark_usedVertices ( mesh,MMG2D_delPt );

  return;
}
