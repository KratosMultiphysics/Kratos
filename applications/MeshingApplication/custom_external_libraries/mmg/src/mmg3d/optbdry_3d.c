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
 * \file mmg3d/optbdry_3d.c
 * \brief Functions for the optimization of very bad elements.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmg3d.h"

/* /\** */
/*  * \param mesh pointer toward the mesh structure. */
/*  * \param met pointer toward the metric structure. */
/*  * \param k   index of a tetra */
/*  * */
/*  * Try to optimize the tetra k. This tetra has a face on the boundary. */
/*  * */
/*  *\/ */
/* int MMG3D_optbdry(MMG5_pMesh mesh,MMG5_pSol met,int k) { */
/*   MMG5_pTetra pt; */
/*   int    iadr; */
/*   int    ia,ib,i,*adja,ipb; */



/*   pt = &mesh->tetra[k]; */
/*   iadr = 4*(k-1) + 1; */
/*   adja = &mesh->adja[iadr]; */

/*   /\* /\\*essai de collapse du point qui n'est pas sur la peau*\\/ *\/ */
/*   /\* for(i=0 ; i<4 ; i++) if(!adja[i]) break; *\/ */

/*   /\* ib  = i; *\/ */
/*   /\* ipb = pt->v[ib]; *\/ */
/*   /\* if(!mesh->info.noinsert) { *\/ */
/*   /\*   for(i=1 ; i<4 ; i++) { *\/ */
/*   /\*     if(!adja[i]) continue; *\/ */
/*   /\*     for(iedg = 0 ; iedg<3 ;iedg++) { *\/ */
/*   /\*       iq = _MMG5_idir[i][_MMG5_iprv2[iedg]]; *\/ */
/*   /\*       if(mesh->point[ *\/ */
/*   /\*     ia = (ib+i)%4; *\/ */
/*   /\*     ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,2); *\/ */
/*   /\*       if ( ilist > 0 ) { *\/ */
/*   /\*         ier = _MMG5_colver(mesh,met,list,ilist,i2,2); *\/ */
/*   /\*         if ( ilist < 0 ) continue; *\/ */
/*   /\*         if ( ier < 0 ) return(-1); *\/ */
/*   /\*     if(MMG_colpoi2(mesh,sol,k,ia,ib,QDEGRAD)) { *\/ */
/*   /\*       MMG_delPt(mesh,ipb); *\/ */
/*   /\*       break; *\/ */
/*   /\*     } *\/ */
/*   /\*   } *\/ */
/*   /\* } else { *\/ */
/*   /\*   i=4; *\/ */
/*   /\* } *\/ */

/*   /\* /\\*try to move the 4th vertex*\\/ *\/ */
/*   /\* if(i==4) { *\/ */
/*   /\*   //if(k==402140) printf("colpoi refused, try move %d %d %d\n",k,ib,pt->v[ib]); *\/ */
/*   /\*   if(!MMG_movevertexbdry(mesh,sol,k,ib)) return(0); *\/ */
/*   /\*   return(2); *\/ */
/*   /\* } *\/ */

/*   return(1); */

/* } */
