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
 * \file mmg2d/libmmg2df.c
 * \brief Fortran API functions for MMG2D library.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmg2d/liblibmmg2d_private.h file for functions
 * documentation.
 *
 * Define the private Fortran API functions for MMG2D library
 * (incompatible functions with the main binary): adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "libmmg2d.h"
#include "mmgcommon_private.h"

/**
 * See \ref MMG2D_mmg2dlib function in \ref mmg2d/liblibmmg2d_private.h file.
 */
FORTRAN_NAME(MMG2D_MMG2DLIB,mmg2d_mmg2dlib,(MMG5_pMesh *mesh,MMG5_pSol *met
                                            ,int* retval),(mesh,met
                                                           ,retval)){

  *retval = MMG2D_mmg2dlib(*mesh,*met);

  return;
}
/**
 * See \ref MMG2D_mmg2dmesh function in \ref mmg2d/liblibmmg2d_private.h file.
 */
FORTRAN_NAME(MMG2D_MMG2DMESH,mmg2d_mmg2dmesh,(MMG5_pMesh *mesh,MMG5_pSol *met
                                              ,int* retval),(mesh,met
                                                             ,retval)){

  *retval = MMG2D_mmg2dmesh(*mesh,*met);

  return;
}
/**
 * See \ref MMG2D_mmg2dls function in \ref mmg2d/liblibmmg2d_private.h file.
 */
FORTRAN_NAME(MMG2D_MMG2DLS,mmg2d_mmg2dls,(MMG5_pMesh *mesh,MMG5_pSol *sol,
                                          MMG5_pSol *met,int* retval),
             (mesh,sol,met,retval)){

  if ( met ) {
    *retval = MMG2D_mmg2dls(*mesh,*sol,*met);
  }
  else {
    *retval = MMG2D_mmg2dls(*mesh,*sol,NULL);
  }

  return;
}
/**
 * See \ref MMG2D_mmg2dmov function in \ref mmg2d/liblibmmg2d_private.h file.
 */
FORTRAN_NAME(MMG2D_MMG2DMOV,mmg2d_mmg2dmov,(MMG5_pMesh *mesh,MMG5_pSol *met,MMG5_pSol *disp
                                            ,int* retval),(mesh,met,disp
                                                           ,retval)){

  *retval = MMG2D_mmg2dmov(*mesh,*met,*disp);

  return;
}
