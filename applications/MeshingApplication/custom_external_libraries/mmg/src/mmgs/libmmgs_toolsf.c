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

#include "libmmgs_private.h"
#include "libmmgs.h"

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
 * See \ref MMGS_Get_numberOfNonBdyEdges function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_GET_NUMBEROFNONBDYEDGES,mmgs_get_numberofnonbdyedges,
             (MMG5_pMesh *mesh,MMG5_int* nb_edges, int* retval),
             (mesh,nb_edges,retval)) {
  *retval =  MMGS_Get_numberOfNonBdyEdges(*mesh,nb_edges);
  return;
}

/**
 * See \ref MMGS_Get_nonBdyEdge function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_GET_NONBDYEDGE,mmgs_get_nonbdyedge,
             (MMG5_pMesh *mesh,MMG5_int* e0, MMG5_int* e1,MMG5_int *ref,MMG5_int* idx,int* retval),
             (mesh,e0,e1,ref,idx,retval)) {
  *retval =  MMGS_Get_nonBdyEdge(*mesh,e0,e1,ref,*idx);
  return;
}

/**
 * See \ref MMGS_usage function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_USAGE,mmgs_usage,
             (char *prog,int *strlen,int *retval),
             (prog,strlen,retval)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,prog,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMGS_usage(tmp);
  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMGS_defaultValues function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_DEFAULTVALUES,mmgs_defaultvalues,
             (MMG5_pMesh *mesh, int* retval),
             (mesh,retval)) {
  *retval = MMGS_defaultValues(*mesh);
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
               (MMG5_pMesh *mesh,MMG5_int* kel, MMG5_int* listri, int* retval),
               (mesh,kel,listri,retval)) {
  *retval =  MMGS_Get_adjaTri(*mesh,*kel,listri);
  return;
}

/**
 * See \ref MMGS_Get_adjaVerticesFast function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_GET_ADJAVERTICESFAST,mmgs_get_adjaverticesfast,
               (MMG5_pMesh *mesh,MMG5_int* ip, MMG5_int *start, MMG5_int* lispoi, MMG5_int* retval),
               (mesh,ip,start,lispoi,retval)) {
  *retval =  MMGS_Get_adjaVerticesFast(*mesh,*ip, *start,lispoi);
  return;
}

/**
 * See \ref MMGS_Free_solutions function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_FREE_SOLUTIONS,mmgs_free_solutions,
             (MMG5_pMesh *mesh,MMG5_pSol *sol),(mesh,sol)) {

  MMGS_Free_solutions(*mesh,*sol);

  return;
}

/**
 * See \ref MMGS_doSol function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_DOSOL,mmgs_dosol,
             (MMG5_pMesh *mesh,MMG5_pSol *met,int *retval),
             (mesh,met,retval)) {
  *retval = MMGS_doSol(*mesh,*met);
  return;
}

/**
 * See \ref MMGS_Set_constantSize function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_SET_CONSTANTSIZE,mmgs_set_constantsize,
             (MMG5_pMesh *mesh,MMG5_pSol *met,int *retval),
             (mesh,met,retval)) {
  *retval =  MMGS_Set_constantSize(*mesh,*met);
  return;
}

/**
 * See \ref MMGS_Compute_eigenv function in \ref mmgs/libmmgs.h file.
 */
FORTRAN_NAME(MMGS_COMPUTE_EIGENV,mmgs_compute_eigenv,
             (double m[6],double lambda[3],double vp[3][3],int *retval),
             (m,lambda,vp,retval)) {

  *retval = MMGS_Compute_eigenv(m,lambda,vp);
  return;
}
