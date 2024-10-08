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
 * \file mmg3d/libmmg3d_toolsf.c
 * \brief Fortran API functions for MMG3D library.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \note Please, refer to the \ref mmg3d/libmmg3d.h file for functions
 * documentation.
 *
 * Define the private Fortran API functions for MMG3D library
 * (incompatible functions with the main binary): adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "libmmg3d_private.h"
#include "libmmg3d.h"

/**
 * See \ref MMG3D_setfunc function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SETFUNC,mmg3d_setfunc,
             (MMG5_pMesh *mesh,MMG5_pSol *met),
             (mesh,met)) {
  MMG3D_setfunc(*mesh,*met);
  return;
}

/**
 * See \ref MMG3D_Get_numberOfNonBdyTriangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_NUMBEROFNONBDYTRIANGLES,mmg3d_get_numberofnonbdytriangles,
             (MMG5_pMesh *mesh,MMG5_int* nb_tria, int* retval),
             (mesh,nb_tria,retval)) {
  *retval =  MMG3D_Get_numberOfNonBdyTriangles(*mesh,nb_tria);
  return;
}

/**
 * See \ref MMG3D_Get_nonBdyTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_NONBDYTRIANGLE,mmg3d_get_nonbdytriangle,
             (MMG5_pMesh *mesh,MMG5_int* v0, MMG5_int* v1,MMG5_int* v2,MMG5_int *ref,
              MMG5_int* idx,int* retval),
             (mesh,v0,v1,v2,ref,idx,retval)) {
  *retval =  MMG3D_Get_nonBdyTriangle(*mesh,v0,v1,v2,ref,*idx);
  return;
}


/**
 * See \ref MMG3D_usage function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_USAGE,mmg3d_usage,
             (char *prog, int *strlen),
             (prog,strlen)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,prog,*strlen);
  tmp[*strlen] = '\0';
  MMG3D_usage(tmp);
  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_parsop function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_PARSOP,mmg3d_parsop,
             (MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),
             (mesh,met,retval)) {

  *retval = MMG3D_parsop(*mesh, *met);

  return;
}

/**
 * See \ref MMG3D_defaultValues function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_DEFAULTVALUES,mmg3d_defaultvalues,
             (MMG5_pMesh *mesh),
             (mesh)) {
  MMG3D_defaultValues(*mesh);
  return;
}

/**
 * See \ref MMG3D_stockOptions function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_STOCKOPTIONS,mmg3d_stockoptions,
             (MMG5_pMesh *mesh, MMG5_Info *info, int* retval),
             (mesh,info,retval)) {
  *retval = MMG3D_stockOptions(*mesh,info);
  return;
}

/**
 * See \ref MMG3D_destockOptions function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_DESTOCKOPTIONS,mmg3d_destockoptions,
               (MMG5_pMesh *mesh, MMG5_Info *info),
             (mesh,info)) {
  MMG3D_destockOptions(*mesh,info);
  return;
}

/**
 * See \ref MMG3D_mmg3dcheck function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_MMG3DCHECK,mmg3d_mmg3dcheck,
             (MMG5_pMesh *mesh,MMG5_pSol *met,MMG5_pSol *sol,double *critmin, double *lmin,
              double *lmax, MMG5_int *eltab,int *metRidTyp,int *retval),
             (mesh,met,sol,critmin,lmin,lmax,eltab,metRidTyp,retval)) {
  int8_t tmp = (int8_t)(*metRidTyp);

  *retval = MMG3D_mmg3dcheck(*mesh,*met,*sol,*critmin,*lmin,*lmax,eltab,tmp);

  return;
}


/**
 * See \ref MMG3D_searchqua function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SEARCHQUA,mmg3d_searchqua,
             (MMG5_pMesh *mesh,MMG5_pSol *met,double *critmin,
              MMG5_int *eltab,int8_t *metRidTyp),
             (mesh,met,critmin,eltab,metRidTyp)) {
  MMG3D_searchqua(*mesh,*met,*critmin,eltab,*metRidTyp);
  return;
}

/**
 * See \ref MMG3D_searchlen function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SEARCHLEN,mmg3d_searchlen,
             (MMG5_pMesh *mesh,MMG5_pSol *met, double *lmin,
              double *lmax, MMG5_int *eltab,int8_t *metRidTyp,int *retval),
             (mesh,met,lmin,lmax,eltab,metRidTyp,retval)) {
  *retval = MMG3D_searchlen(*mesh,*met,*lmin,*lmax,eltab,*metRidTyp);
  return;
}

/**
 * See \ref MMG3D_Get_tetFromTria function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TETFROMTRIA,mmg3d_get_tetfromtria,
             (MMG5_pMesh *mesh,MMG5_int *ktri, MMG5_int *ktet, int *iface,int *retval),
             (mesh,ktri,ktet,iface,retval)) {

  *retval = MMG3D_Get_tetFromTria(*mesh,*ktri,ktet,iface);
  return;
}

/**
 * See \ref MMG3D_Get_tetsFromTria function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TETSFROMTRIA,mmg3d_get_tetsfromtria,
             (MMG5_pMesh *mesh,MMG5_int *ktri, MMG5_int ktet[2], int iface[2],int *retval),
             (mesh,ktri,ktet,iface,retval)) {

  *retval = MMG3D_Get_tetsFromTria(*mesh,*ktri,ktet,iface);
  return;
}

/**
 * See \ref MMG3D_Get_adjaTet function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_ADJATET,mmg3d_get_adjatet,
               (MMG5_pMesh *mesh,MMG5_int* kel, MMG5_int listet[4], int* retval),
               (mesh,kel,listet,retval)) {
  *retval =  MMG3D_Get_adjaTet(*mesh,*kel,listet);
  return;
}

/**
 * See \ref MMG3D_Free_solutions function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_FREE_SOLUTIONS,mmg3d_free_solutions,
             (MMG5_pMesh *mesh,MMG5_pSol *sol),(mesh,sol)) {

  MMG3D_Free_solutions(*mesh,*sol);

  return;
}

/**
 * See \ref MMG3D_doSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_DOSOL,mmg3d_dosol,
             (MMG5_pMesh *mesh,MMG5_pSol *met,int *retval),
             (mesh,met,retval)) {
  *retval = MMG3D_doSol(*mesh,*met);
  return;
}

/**
 * See \ref MMG3D_Set_constantSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_CONSTANTSIZE,mmg3d_set_constantsize,
             (MMG5_pMesh *mesh,MMG5_pSol *met,int *retval),
             (mesh,met,retval)) {
  *retval =  MMG3D_Set_constantSize(*mesh,*met);
  return;
}

/**
 * See \ref MMG3D_Compute_eigenv function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_COMPUTE_EIGENV,mmg3d_compute_eigenv,
             (double m[6],double lambda[3],double vp[3][3],int *retval),
             (m,lambda,vp,retval)) {

  *retval = MMG3D_Compute_eigenv(m,lambda,vp);
  return;
}
