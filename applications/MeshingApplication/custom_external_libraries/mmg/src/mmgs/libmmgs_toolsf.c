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
 * \file mmgs/libmmgs_toolsf.c
 * \brief Fortran API functions for MMGS library.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmgs/libmmgs.h file for functions
 * documentation.
 *
 * Define the private Fortran API functions for MMGS library
 * (incompatible functions with the main binary): adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "mmgs.h"

/**
 * See \ref MMGS_setfunc function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_SETFUNC,mmgs_setfunc,
             (MMG5_pMesh *mesh,MMG5_pSol *met),
             (mesh,met)) {
  MMGS_setfunc(*mesh,*met);
  return;
}

/**
 * See \ref MMGS_usage function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_USAGE,mmgs_usage,
             (char *prog,int *strlen),
             (prog,strlen)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,prog,*strlen);
  tmp[*strlen] = '\0';
  MMGS_usage(tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMGS_defaultValues function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_DEFAULTVALUES,mmgs_defaultvalues,
             (MMG5_pMesh *mesh),
             (mesh)) {
  MMGS_defaultValues(*mesh);
  return;
}

/**
 * See \ref MMGS_stockOptions function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_STOCKOPTIONS,mmgs_stockoptions,
               (MMG5_pMesh *mesh, MMG5_Info *info, int* retval),
               (mesh,info,retval)) {
  *retval =  MMGS_stockOptions(*mesh,info);
  return;
}

/**
 * See \ref MMGS_destockOptions function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_DESTOCKOPTIONS,mmgs_destockoptions,
               (MMG5_pMesh *mesh, MMG5_Info *info),
               (mesh,info)) {
  MMGS_destockOptions(*mesh,info);
  return;
}

/**
 * See \ref MMGS_Get_adjaTri function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_GET_ADJATRI,mmgs_get_adjatri,
               (MMG5_pMesh *mesh,int* kel, int* listri, int* retval),
               (mesh,kel,listri,retval)) {
  *retval =  MMGS_Get_adjaTri(*mesh,*kel,listri);
  return;
}

/**
 * See \ref MMGS_Get_adjaVerticesFast function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_GET_ADJAVERTICESFAST,mmgs_get_adjaverticesfast,
               (MMG5_pMesh *mesh,int* ip, int *start, int* lispoi, int* retval),
               (mesh,ip,start,lispoi,retval)) {
  *retval =  MMGS_Get_adjaVerticesFast(*mesh,*ip, *start,lispoi);
  return;
}
