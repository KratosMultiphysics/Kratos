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
 * \file mmg3d/API_functionsf_3d.c
 * \brief Fortran API functions for MMG3D library.
 * \author Algiane Froehly  (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete soon...
 * \note Please, refer to the \ref mmg3d/libmmg3d.h file for functions
 * documentation.
 *
 * Define the Fortran API functions for MMG3D library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */

#include "libmmg3d_private.h"
#include "libmmg3d.h"

/**
 * See \ref MMG3D_Init_mesh function in common/libmmgcommon_private.h file.
 */
FORTRAN_VARIADIC ( MMG3D_INIT_MESH, mmg3d_init_mesh,
                   (const int starter, ... ),
                   va_list argptr;
                   int     ier;

                   va_start(argptr, starter);

                   ier = MMG3D_Init_mesh_var(argptr);

                   va_end(argptr);

                   if ( !ier ) exit(EXIT_FAILURE);

                   return;
  )

/**
 * See \ref MMG3D_Init_parameters function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_INIT_PARAMETERS,mmg3d_init_parameters,(MMG5_pMesh *mesh),(mesh)) {
  MMG3D_Init_parameters(*mesh);
  return;
}

/**
 * See \ref MMG3D_Set_inputMeshName function in \ref mmg3d/libmmg3d.h.
 */
FORTRAN_NAME(MMG3D_SET_INPUTMESHNAME, mmg3d_set_inputmeshname,
             (MMG5_pMesh *mesh, char* meshin, int *strlen0, int* retval),
             (mesh,meshin,strlen0,retval)) {
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,meshin,*strlen0);
  tmp[*strlen0] = '\0';
  *retval = MMG3D_Set_inputMeshName(*mesh,tmp);
  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_inputSolName function in \ref mmg3d/libmmg3d.h.
 */
FORTRAN_NAME(MMG3D_SET_INPUTSOLNAME, mmg3d_set_inputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solin, int* strlen0, int* retval),
             (mesh,sol,solin,strlen0,retval)) {

  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,solin,*strlen0);
  tmp[*strlen0] = '\0';
  *retval = MMG3D_Set_inputSolName(*mesh,*sol,tmp);
  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_outputMeshName function in mmg3d/libmmg3d.h.
 */
FORTRAN_NAME(MMG3D_SET_OUTPUTMESHNAME,mmg3d_set_outputmeshname,
             (MMG5_pMesh *mesh, char* meshout, int* strlen0,int* retval),
             (mesh,meshout,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,meshout,*strlen0);
  tmp[*strlen0] = '\0';
  *retval = MMG3D_Set_outputMeshName(*mesh, tmp);
  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_outputSolName function in \ref common/libmmgcommon_private.h file.
 */
FORTRAN_NAME(MMG3D_SET_OUTPUTSOLNAME,mmg3d_set_outputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solout,int* strlen0, int* retval),
             (mesh,sol,solout,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,solout,*strlen0);
  tmp[*strlen0] = '\0';
  *retval = MMG3D_Set_outputSolName(*mesh,*sol,tmp);
  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_SOLSIZE,mmg3d_set_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity,
              MMG5_int* np, int* typSol, int* retval),
             (mesh, sol, typEntity, np, typSol, retval)) {
  *retval = MMG3D_Set_solSize(*mesh,*sol,*typEntity,*np,*typSol);
  return;
}

/**
 * See \ref MMG3D_Set_solAtVerticesSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_SOLSATVERTICESSIZE,mmg3d_set_solsatverticessize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,int* nsols,
              MMG5_int* nentities, int* typSol, int* retval),
             (mesh, sol, nsols, nentities, typSol, retval)) {
  *retval = MMG3D_Set_solsAtVerticesSize(*mesh,sol,*nsols,*nentities,typSol);
  return;
}

/**
 * See \ref MMG3D_Set_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_MESHSIZE,mmg3d_set_meshsize,
             (MMG5_pMesh *mesh, MMG5_int *np, MMG5_int *ne, MMG5_int *nprism,
              MMG5_int *nt, MMG5_int *nquad, MMG5_int *na, int *retval),
             (mesh,np,ne,nprism,nt,nquad,na,retval)) {
  *retval = MMG3D_Set_meshSize(*mesh,*np,*ne,*nprism,*nt,*nquad,*na);
  return;
}

/**
 * See \ref MMG3D_Get_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_SOLSIZE,mmg3d_get_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity, MMG5_int* np, int* typSol, int* retval),
             (mesh,sol,typEntity,np,typSol,retval)) {

  *retval = MMG3D_Get_solSize(*mesh,*sol,typEntity,np,typSol);
  return;
}

/**
 * See \ref MMG3D_Get_solsatverticessize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_SOLSATVERTICESSIZE,mmg3d_get_solsatverticessize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *nsols,
              MMG5_int* nentities, int* typSol, int* retval),
             (mesh,sol,nsols,nentities,typSol,retval)) {

  *retval = MMG3D_Get_solsAtVerticesSize(*mesh,sol,nsols,nentities,typSol);
  return;
}

/**
 * See \ref MMG3D_Get_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_MESHSIZE,mmg3d_get_meshsize,
             (MMG5_pMesh *mesh, MMG5_int* np, MMG5_int* ne,MMG5_int *nprism,
              MMG5_int* nt,MMG5_int *nquad, MMG5_int* na, int* retval),
             (mesh,np,ne,nprism,nt,nquad,na,retval)) {

  *retval = MMG3D_Get_meshSize(*mesh,np,ne,nprism,nt,nquad,na);
  return;
}

/**
 * See \ref MMG3D_Set_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_VERTEX,mmg3d_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
              MMG5_int* pos, int* retval),
             (mesh,c0,c1,c2,ref,pos,retval)) {

  *retval = MMG3D_Set_vertex(*mesh,*c0,*c1,*c2,*ref,*pos);
  return;
}
/**
 * See \ref MMG3D_Get_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_VERTEX,mmg3d_get_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
              int* isCorner, int* isRequired, int* retval),
             (mesh,c0,c1,c2,ref,isCorner,isRequired, retval)) {
  *retval = MMG3D_Get_vertex(*mesh,c0,c1,c2,ref,isCorner,isRequired);
  return;
}
/**
 * See \ref MMG3D_GetByIdx_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GETBYIDX_VERTEX,mmg3d_getbyidx_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
              int* isCorner, int* isRequired, MMG5_int* idx,int* retval),
             (mesh,c0,c1,c2,ref,isCorner,isRequired,idx, retval)) {
  *retval = MMG3D_GetByIdx_vertex(*mesh,c0,c1,c2,ref,isCorner,isRequired,*idx);
  return;
}

/**
 * See \ref MMG3D_Set_vertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_VERTICES,mmg3d_set_vertices,
             (MMG5_pMesh *mesh, double* vertices, MMG5_int* refs, int* retval),
             (mesh,vertices,refs,retval)) {

  *retval = MMG3D_Set_vertices(*mesh,vertices,refs);
  return;
}


/**
 * See \ref MMG3D_Get_vertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_VERTICES,mmg3d_get_vertices,
             (MMG5_pMesh *mesh, double* vertices, MMG5_int* refs,
              int* areCorners, int* areRequired, int* retval),
             (mesh,vertices,refs,areCorners,areRequired, retval)) {
  *retval = MMG3D_Get_vertices(*mesh,vertices,refs,areCorners,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TETRAHEDRON,mmg3d_set_tetrahedron,
             (MMG5_pMesh *mesh, MMG5_int *v0, MMG5_int *v1, MMG5_int *v2, MMG5_int *v3, MMG5_int *ref,
              MMG5_int *pos, int* retval),
             (mesh,v0,v1,v2,v3,ref,pos,retval)){
  *retval = MMG3D_Set_tetrahedron(*mesh,*v0,*v1,*v2,*v3,*ref,*pos);
  return;
}

/**
 * See \ref MMG3D_Set_tetrahedra function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TETRAHEDRA,mmg3d_set_tetrahedra,
             (MMG5_pMesh *mesh, MMG5_int *tetra, MMG5_int *refs, int* retval),
             (mesh,tetra,refs,retval)){
  *retval = MMG3D_Set_tetrahedra(*mesh,tetra,refs);
  return;
}

/**
 * See \ref MMG3D_Get_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TETRAHEDRON,mmg3d_get_tetrahedron,
             (MMG5_pMesh *mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int *v3,
              MMG5_int* ref, int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = MMG3D_Get_tetrahedron(*mesh,v0,v1,v2,v3,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Get_tetrahedra function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TETRAHEDRA,mmg3d_get_tetrahedra,
             (MMG5_pMesh *mesh, MMG5_int* tetra, MMG5_int* refs, int* areRequired,
              int* retval),
             (mesh,tetra,refs,areRequired,retval)) {
  *retval = MMG3D_Get_tetrahedra(*mesh,tetra,refs,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_prism function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_PRISM,mmg3d_set_prism,
             (MMG5_pMesh *mesh, MMG5_int *v0, MMG5_int *v1, MMG5_int *v2, MMG5_int *v3,
              MMG5_int *v4, MMG5_int *v5,MMG5_int *ref,MMG5_int *pos, int* retval),
             (mesh,v0,v1,v2,v3,v4,v5,ref,pos,retval)){
  *retval = MMG3D_Set_prism(*mesh,*v0,*v1,*v2,*v3,*v4,*v5,*ref,*pos);
  return;
}

/**
 * See \ref MMG3D_Set_prisms function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_PRISMS,mmg3d_set_prisms,
             (MMG5_pMesh *mesh, MMG5_int *prisms, MMG5_int *refs, int* retval),
             (mesh,prisms,refs,retval)){
  *retval = MMG3D_Set_prisms(*mesh,prisms,refs);
  return;
}

/**
 * See \ref MMG3D_Get_prism function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_PRISM,mmg3d_get_prism,
             (MMG5_pMesh *mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int *v3,
              MMG5_int *v4, MMG5_int* v5,MMG5_int* ref, int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,v4,v5,ref,isRequired,retval)) {
  *retval = MMG3D_Get_prism(*mesh,v0,v1,v2,v3,v4,v5,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Get_prisms function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_PRISMS,mmg3d_get_prisms,
             (MMG5_pMesh *mesh, MMG5_int* prisms, MMG5_int* refs, int* areRequired,
              int* retval),
             (mesh,prisms,refs,areRequired,retval)) {
  *retval = MMG3D_Get_prisms(*mesh,prisms,refs,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TRIANGLE,mmg3d_set_triangle,
             (MMG5_pMesh *mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref,MMG5_int* pos,
              int* retval),
             (mesh,v0,v1,v2,ref,pos,retval)) {
  *retval = MMG3D_Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
  return;
}

/**
 * See \ref MMG3D_Get_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TRIANGLE,mmg3d_get_triangle,
             (MMG5_pMesh *mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref
              ,int* isRequired, int* retval),
             (mesh,v0,v1,v2,ref,isRequired,retval)) {
  *retval = MMG3D_Get_triangle(*mesh,v0,v1,v2,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Set_triangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TRIANGLES,mmg3d_set_triangles,
             (MMG5_pMesh *mesh, MMG5_int* tria, MMG5_int* refs,
              int* retval),
             (mesh,tria,refs,retval)) {
  *retval = MMG3D_Set_triangles(*mesh, tria, refs);
  return;
}

/**
 * See \ref MMG3D_Get_triangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TRIANGLES,mmg3d_get_triangles,
             (MMG5_pMesh *mesh, MMG5_int* tria, MMG5_int* refs,int* areRequired,
              int* retval),
             (mesh,tria,refs,areRequired,retval)) {
  *retval = MMG3D_Get_triangles(*mesh,tria,refs,areRequired);
  return;
}
/**
 * See \ref MMG3D_Set_quadrilateral function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_QUADRILATERAL,mmg3d_set_quadrilateral,
             (MMG5_pMesh *mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,MMG5_int *v3,
              MMG5_int* ref,MMG5_int* pos,int* retval),
             (mesh,v0,v1,v2,v3,ref,pos,retval)) {
  *retval = MMG3D_Set_quadrilateral(*mesh, *v0, *v1, *v2, *v3,*ref, *pos);
  return;
}

/**
 * See \ref MMG3D_Get_quadrilateral function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_QUADRILATERAL,mmg3d_get_quadrilateral,
             (MMG5_pMesh *mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,MMG5_int *v3,
              MMG5_int* ref,int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = MMG3D_Get_quadrilateral(*mesh,v0,v1,v2,v3,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Set_quadrilaterals function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_QUADRILATERALS,mmg3d_set_quadrilaterals,
             (MMG5_pMesh *mesh, MMG5_int* quads, MMG5_int* refs,
              int* retval),
             (mesh,quads,refs,retval)) {
  *retval = MMG3D_Set_quadrilaterals(*mesh, quads, refs);
  return;
}

/**
 * See \ref MMG3D_Get_quadrilaterals function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_QUADRILATERALS,mmg3d_get_quadrilaterals,
             (MMG5_pMesh *mesh, MMG5_int* quads, MMG5_int* refs,int* areRequired,
              int* retval),
             (mesh,quads,refs,areRequired,retval)) {
  *retval = MMG3D_Get_quadrilaterals(*mesh,quads,refs,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_EDGE,mmg3d_set_edge,
             (MMG5_pMesh *mesh, MMG5_int *v0, MMG5_int *v1, MMG5_int *ref, MMG5_int *pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)){
  *retval = MMG3D_Set_edge(*mesh,*v0,*v1,*ref,*pos);
  return;
}

/**
 * See \ref MMG3D_Get_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_EDGE,mmg3d_get_edge,(MMG5_pMesh *mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref
                                            ,int* isRidge, int* isRequired, int* retval),
             (mesh,e0,e1,ref,isRidge,isRequired,retval)) {
  *retval = MMG3D_Get_edge(*mesh,e0,e1,ref,isRidge,isRequired);
  return;
}
/**
 * See \ref MMG3D_Set_edges function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_EDGES,mmg3d_set_edges,
             (MMG5_pMesh *mesh, MMG5_int *edges, MMG5_int *refs, int* retval),
             (mesh,edges,refs,retval)){
  *retval = MMG3D_Set_edges(*mesh,edges,refs);
  return;
}

/**
 * See \ref MMG3D_Get_edges function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_EDGES,mmg3d_get_edges,(MMG5_pMesh *mesh, MMG5_int* edges,
                                              MMG5_int* refs,int* areRidges,
                                              int* areRequired, int* retval),
             (mesh,edges,refs,areRidges,areRequired,retval)) {
  *retval = MMG3D_Get_edges(*mesh,edges,refs,areRidges,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_corner function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_CORNER,mmg3d_set_corner,(MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG3D_Set_corner(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Unset_corner function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_CORNER,mmg3d_unset_corner,(MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG3D_Unset_corner(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDVERTEX,mmg3d_set_requiredvertex,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG3D_Set_requiredVertex(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Unset_requiredVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_REQUIREDVERTEX,mmg3d_unset_requiredvertex,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG3D_Unset_requiredVertex(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredTetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTETRAHEDRON,mmg3d_set_requiredtetrahedron,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_requiredTetrahedron(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Unset_requiredTetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_REQUIREDTETRAHEDRON,mmg3d_unset_requiredtetrahedron,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Unset_requiredTetrahedron(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredTetrahedra function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTETRAHEDRA,mmg3d_set_requiredtetrahedra,
             (MMG5_pMesh *mesh, MMG5_int *reqIdx, MMG5_int *nreq, int* retval),
             (mesh,reqIdx,nreq,retval)) {
  *retval = MMG3D_Set_requiredTetrahedra(*mesh,reqIdx, *nreq);
  return;
}

/**
 * See \ref MMG3D_Unset_requiredTetrahedra function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_REQUIREDTETRAHEDRA,mmg3d_unset_requiredtetrahedra,
             (MMG5_pMesh *mesh, MMG5_int *reqIdx, MMG5_int *nreq, int* retval),
             (mesh,reqIdx,nreq,retval)) {
  *retval = MMG3D_Unset_requiredTetrahedra(*mesh,reqIdx, *nreq);
  return;
}
/**
 * See \ref MMG3D_Set_requiredTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTRIANGLE,mmg3d_set_requiredtriangle,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_requiredTriangle(*mesh, *k);
  return;
}
/**
 * See \ref MMG3D_Unset_requiredTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_REQUIREDTRIANGLE,mmg3d_unset_requiredtriangle,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Unset_requiredTriangle(*mesh, *k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredTriangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTRIANGLES,mmg3d_set_requiredtriangles,
             (MMG5_pMesh *mesh, MMG5_int *reqIdx, MMG5_int *nreq, int* retval),
             (mesh,reqIdx,nreq,retval)) {
  *retval = MMG3D_Set_requiredTriangles(*mesh, reqIdx, *nreq);
  return;
}

/**
 * See \ref MMG3D_Unset_requiredTriangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_REQUIREDTRIANGLES,mmg3d_unset_requiredtriangles,
             (MMG5_pMesh *mesh, MMG5_int *reqIdx, MMG5_int *nreq, int* retval),
             (mesh,reqIdx,nreq,retval)) {
  *retval = MMG3D_Unset_requiredTriangles(*mesh, reqIdx, *nreq);
  return;
}

/**
 * See \ref MMG3D_Set_parallelTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_PARALLELTRIANGLE,mmg3d_set_paralleltriangle,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_parallelTriangle(*mesh, *k);
  return;
}

/**
 * See \ref MMG3D_Unset_parallelTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_PARALLELTRIANGLE,mmg3d_unset_paralleltriangle,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Unset_parallelTriangle(*mesh, *k);
  return;
}

/**
 * See \ref MMG3D_Set_parallelTriangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_PARALLELTRIANGLES,mmg3d_set_paralleltriangles,
             (MMG5_pMesh *mesh, MMG5_int *parIdx, MMG5_int *npar, int* retval),
             (mesh,parIdx,npar,retval)) {
  *retval = MMG3D_Set_parallelTriangles(*mesh, parIdx, *npar);
  return;
}

/**
 * See \ref MMG3D_Unset_parallelTriangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_PARALLELTRIANGLES,mmg3d_unset_paralleltriangles,
             (MMG5_pMesh *mesh, MMG5_int *parIdx, MMG5_int *npar, int* retval),
             (mesh,parIdx,npar,retval)) {
  *retval = MMG3D_Unset_parallelTriangles(*mesh, parIdx, *npar);
  return;
}

/**
 * See \ref MMG3D_Set_ridge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_RIDGE,mmg3d_set_ridge,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_ridge(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_unset_ridge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_RIDGE,mmg3d_unset_ridge,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Unset_ridge(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredEdge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDEDGE,mmg3d_set_requirededge,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_requiredEdge(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Unset_requiredEdge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_UNSET_REQUIREDEDGE,mmg3d_unset_requirededge,
             (MMG5_pMesh *mesh, MMG5_int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Unset_requiredEdge(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_normalAtVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_NORMALATVERTEX,mmg3d_set_normalatvertex,
             (MMG5_pMesh *mesh, MMG5_int *k, double* n0, double* n1, double* n2,int* retval),
             (mesh,k,n0,n1,n2,retval)) {
  *retval = MMG3D_Set_normalAtVertex(*mesh,*k, *n0, *n1, *n2);
  return;
}
/**
 * See \ref MMG3D_Get_normalAtVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_NORMALATVERTEX,mmg3d_get_normalatvertex,
             (MMG5_pMesh *mesh, MMG5_int *k, double* n0, double* n1, double* n2,int* retval),
             (mesh,k,n0,n1,n2,retval)) {
  *retval = MMG3D_Get_normalAtVertex(*mesh,*k, n0, n1, n2);
  return;
}
/**
 * See \ref MMG3D_Get_tetrahedronQuality function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TETRAHEDRONQUALITY,mmg3d_get_tetrahedronquality,
             (MMG5_pMesh *mesh, MMG5_pSol *met, MMG5_int* k, double* retval),
             (mesh,met,k,retval)) {
  *retval = MMG3D_Get_tetrahedronQuality(*mesh,*met,*k);
  return;
}
/**
 * See \ref MMG3D_Set_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_SCALARSOL,mmg3d_set_scalarsol,
             (MMG5_pSol *met, double *s, MMG5_int *pos, int* retval),
             (met,s,pos,retval)) {
  *retval = MMG3D_Set_scalarSol(*met,*s,*pos);
  return;
}

/**
 * See \ref MMG3D_Get_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_SCALARSOL,mmg3d_get_scalarsol,
             (MMG5_pSol *met, double* s, int* retval),
             (met,s,retval)) {
  *retval = MMG3D_Get_scalarSol(*met,s);
  return;
}
/**
 * See \ref MMG3D_Set_scalarSols function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_SCALARSOLS,mmg3d_set_scalarsols,
             (MMG5_pSol *met, double *s, int* retval),
             (met,s,retval)) {
  *retval = MMG3D_Set_scalarSols(*met,s);
  return;
}

/**
 * See \ref MMG3D_Get_scalarSols function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_SCALARSOLS,mmg3d_get_scalarsols,
             (MMG5_pSol *met, double* s, int* retval),
             (met,s,retval)) {
  *retval = MMG3D_Get_scalarSols(*met,s);
  return;
}

/**
 * See \ref MMG3D_Set_vectorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_VECTORSOL,mmg3d_set_vectorsol,
             (MMG5_pSol *met, double *vx, double *vy, double *vz,
              MMG5_int *pos, int* retval),
             (met,vx,vy,vz,pos,retval)) {
  *retval = MMG3D_Set_vectorSol(*met,*vx,*vy,*vz,*pos);
  return;
}

/**
 * See \ref MMG3D_Get_vectorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_VECTORSOL,mmg3d_get_vectorsol,
             (MMG5_pSol *met, double* vx,double *vy, double *vz, int* retval),
             (met,vx,vy,vz,retval)) {
  *retval = MMG3D_Get_vectorSol(*met,vx,vy,vz);
  return;
}
/**
 * See \ref MMG3D_Set_vectorSols function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_VECTORSOLS,mmg3d_set_vectorsols,
             (MMG5_pSol *met, double *sols, int* retval),
             (met,sols,retval)) {
  *retval = MMG3D_Set_vectorSols(*met,sols);
  return;
}

/**
 * See \ref MMG3D_Get_vectorSols function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_VECTORSOLS,mmg3d_get_vectorsols,
             (MMG5_pSol *met, double* sols, int* retval),
             (met,sols,retval)) {
  *retval = MMG3D_Get_vectorSols(*met,sols);
  return;
}

/**
 * See \ref MMG3D_Set_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TENSORSOL,mmg3d_set_tensorsol,
             (MMG5_pSol *met, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, MMG5_int *pos, int* retval),
             (met,m11,m12,m13,m22,m23,m33,pos,retval)) {
  *retval = MMG3D_Set_tensorSol(*met,*m11,*m12,*m13,*m22,*m23,*m33,*pos);
  return;
}

/**
 * See \ref MMG3D_Get_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TENSORSOL,mmg3d_get_tensorsol,
             (MMG5_pSol *met, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int* retval),
             (met,m11,m12,m13,m22,m23,m33,retval)) {
  *retval = MMG3D_Get_tensorSol(*met,m11,m12,m13,m22,m23,m33);
  return;
}
/**
 * See \ref MMG3D_Set_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TENSORSOLS,mmg3d_set_tensorsols,
             (MMG5_pSol *met, double* sols,int* retval),
             (met,sols,retval)) {
  *retval = MMG3D_Set_tensorSols(*met,sols);
  return;
}

/**
 * See \ref MMG3D_Get_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TENSORSOLS,mmg3d_get_tensorsols,
             (MMG5_pSol *met, double* sols, int* retval),
             (met,sols,retval)) {
  *retval = MMG3D_Get_tensorSols(*met,sols);
  return;
}
/**
 * See \ref MMG3D_Set_ithSol_solsAtVertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_ITHSOL_INSOLSATVERTICES,mmg3d_set_ithsol_insolsatvertices,
             (MMG5_pSol *sol, MMG5_int *i,double *s,MMG5_int *pos, int* retval),
             (sol,i,s,pos,retval)) {
  *retval = MMG3D_Set_ithSol_inSolsAtVertices(*sol,*i,s,*pos);
  return;
}
/**
 * See \ref MMG3D_Get_ithSol_inSolsAtVertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_ITHSOL_INSOLSATVERTICES,mmg3d_get_ithsol_insolsatvertices,
             (MMG5_pSol *sol, MMG5_int* i,double *s,MMG5_int *pos, int* retval),
             (sol,i,s,pos,retval)) {
  *retval = MMG3D_Get_ithSol_inSolsAtVertices(*sol,*i,s,*pos);
  return;
}
/**
 * See \ref MMG3D_Set_ithSols_inSolsAtVertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_ITHSOLS_INSOLSATVERTICES,mmg3d_set_ithsols_insolsatvertices,
             (MMG5_pSol *sol, MMG5_int *i,double *s, int* retval),
             (sol,i,s,retval)) {
  *retval = MMG3D_Set_ithSols_inSolsAtVertices(*sol,*i,s);
  return;
}
/**
 * See \ref MMG3D_Get_ithSols_inSolsAtVertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_ITHSOLS_INSOLSATVERTICES,mmg3d_get_ithsols_insolsatvertices,
             (MMG5_pSol *sol, MMG5_int* i,double *s, int* retval),
             (sol,i,s,retval)) {
  *retval = MMG3D_Get_ithSols_inSolsAtVertices(*sol,*i,s);
  return;
}

/**
 * See \ref MMG3D_Set_handGivenMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_HANDGIVENMESH,mmg3d_set_handgivenmesh,
             (MMG5_pMesh *mesh),
             (mesh)) {
  MMG3D_Set_handGivenMesh(*mesh);
  return;
}

/**
 * See \ref MMG3D_Chk_meshData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_CHK_MESHDATA,mmg3d_chk_meshdata,
             (MMG5_pMesh *mesh,MMG5_pSol *met, int* retval),
             (mesh,met,retval)) {
  *retval = MMG3D_Chk_meshData(*mesh,*met);
  return;
}

/**
 * See \ref MMG3D_Add_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_ADD_TETRAHEDRON,mmg3d_add_tetrahedron,
             (MMG5_pMesh *mesh, MMG5_int *v0, MMG5_int *v1, MMG5_int *v2, MMG5_int *v3, MMG5_int *ref,
              int* retval),
             (mesh,v0,v1,v2,v3,ref,retval)){
  *retval = MMG3D_Add_tetrahedron(*mesh,*v0,*v1,*v2,*v3,*ref);
  return;
}

/**
 * See \ref MMG3D_Add_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_ADD_VERTEX,mmg3d_add_vertex,
             (MMG5_pMesh *mesh, double *c0, double *c1, double *c2, MMG5_int *ref,
              MMG5_int* retval),
             (mesh,c0,c1,c2,ref,retval)){
  *retval = MMG3D_Add_vertex(*mesh,*c0,*c1,*c2,*ref);
  return;
}

/**
 * See \ref MMG3D_Set_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_IPARAMETER,mmg3d_set_iparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, MMG5_int *val, int* retval),
             (mesh,sol,iparam,val,retval)){
  *retval = MMG3D_Set_iparameter(*mesh,*sol,*iparam,*val);
  return;
}

/**
 * See \ref MMG3D_Get_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_IPARAMETER,mmg3d_get_iparameter,
             (MMG5_pMesh *mesh, int *iparam, MMG5_int* retval),
             (mesh,iparam,retval)){
  *retval = MMG3D_Get_iparameter(*mesh,*iparam);
  return;
}

/**
 * See \ref MMG3D_Set_dparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_DPARAMETER,mmg3d_set_dparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *dparam, double *val, int* retval),
             (mesh,sol,dparam,val,retval)){
  *retval = MMG3D_Set_dparameter(*mesh,*sol,*dparam,*val);
  return;
}

/**
 * See \ref MMG3D_Set_localParameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_LOCALPARAMETER,mmg3d_set_localparameter,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, int *typ, MMG5_int *ref,
              double *hmin, double *hmax, double *hausd, int* retval),
             (mesh,sol,typ,ref,hmin, hmax, hausd,retval)){
  *retval = MMG3D_Set_localParameter(*mesh,*sol,*typ,*ref,*hmin,*hmax,*hausd);
  return;
}

/**
 * See \ref MMG3D_Set_multiMat function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_MULTIMAT,mmg3d_set_multimat,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, MMG5_int *ref,int *split,
              MMG5_int* rin,MMG5_int* rex, int* retval),
             (mesh,sol,ref,split,rin,rex,retval)){
  *retval = MMG3D_Set_multiMat(*mesh,*sol,*ref,*split,*rin,*rex);
  return;
}

/**
 * See \ref MMG3D_Set_lsBaseReference function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_LSBASEREFERENCE,mmg3d_set_lsbasereference,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, int *br, int* retval),
             (mesh,sol,br,retval)){
  *retval = MMG3D_Set_lsBaseReference(*mesh,*sol,*br);
  return;
}

/**
 * See \ref MMG3D_Free_allSols function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_FREE_ALLSOLS,mmg3d_free_allsols,
             (MMG5_pMesh *mesh,MMG5_pSol *sol,int* retval),
             (mesh,sol,retval)){

  *retval = MMG3D_Free_allSols(*mesh,sol);

  return;
}
/**
 * See \ref MMG3D_Free_all function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_VARIADIC(MMG3D_FREE_ALL,mmg3d_free_all,
                 (const int starter,...),
                 va_list argptr;
                 int     ier;

                 va_start(argptr, starter);

                 ier = MMG3D_Free_all_var(argptr);

                 va_end(argptr);

                 if ( !ier ) exit(EXIT_FAILURE);
                 return;
  )

/**
 * See \ref MMG3D_Free_structures function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_VARIADIC(MMG3D_FREE_STRUCTURES,mmg3d_free_structures,
                 (const int starter,...),
                 va_list argptr;
                 int     ier;

                 va_start(argptr, starter);

                 ier = MMG3D_Free_structures_var(argptr);

                 va_end(argptr);

                 if ( !ier ) exit(EXIT_FAILURE);

                 return;
  )

/**
 * See \ref MMG3D_Free_names function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_VARIADIC(MMG3D_FREE_NAMES,mmg3d_free_names,
                 (const int starter,...),
                 va_list argptr;
                 int     ier;

                 va_start(argptr, starter);

                 ier = MMG3D_Free_names_var(argptr);

                 va_end(argptr);

                 if ( !ier ) exit(EXIT_FAILURE);
                 return;
  )


/**
 * See \ref MMG3D_loadMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADMESH,mmg3d_loadmesh,
             (MMG5_pMesh *mesh,char* filename, int *strlen0,int* retval),
             (mesh,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadMesh(*mesh,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadVtuMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADVTUMESH,mmg3d_loadvtumesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadVtuMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}
/**
 * See \ref MMG3D_loadVtuMesh_and_allData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADVTUMESH_AND_ALLDATA,mmg3d_loadvtumesh_and_alldata,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadVtuMesh_and_allData(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadVtkMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADVTKMESH,mmg3d_loadvtkmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadVtkMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}
/**
 * See \ref MMG3D_loadvtkMesh_and_allData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADVTKMESH_AND_ALLDATA,mmg3d_loadvtkmesh_and_alldata,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadVtkMesh_and_allData(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}
/**
 * See \ref MMG3D_loadMshMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADMSHMESH,mmg3d_loadmshmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadMshMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadGenericMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADGENERICMESH,mmg3d_loadgenericmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadGenericMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadMshMesh_and_allData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADMSHMESH_AND_ALLDATA,mmg3d_loadmshmesh_and_alldata,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadMshMesh_and_allData(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEMESH,mmg3d_savemesh,
             (MMG5_pMesh *mesh,char* filename, int *strlen0,int* retval),
             (mesh,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveMesh(*mesh,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveVtkMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEVTKMESH,mmg3d_savevtkmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,
              int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveVtkMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}
/**
 * See \ref MMG3D_saveVtkMesh_and_allData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEVTKMESH_AND_ALLDATA,mmg3d_savevtkmesh_and_alldata,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,
              int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveVtkMesh_and_allData(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveVtuMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEVTUMESH,mmg3d_savevtumesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,
              int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveVtuMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}
/**
 * See \ref MMG3D_saveVtuMesh_and_allData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEVTUMESH_AND_ALLDATA,mmg3d_savevtumesh_and_alldata,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,
              int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveVtuMesh_and_allData(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveMshMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEMSHMESH,mmg3d_savemshmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,
              int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveMshMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveMshMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEMSHMESH_AND_ALLDATA,mmg3d_savemshmesh_and_alldata,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,
              int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveMshMesh_and_allData(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveTetgenMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVETETGENMESH,mmg3d_savetetgenmesh,(MMG5_pMesh *mesh,char *meshin,int *strlen0, int* retval),
             (mesh,meshin,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,meshin,*strlen0);
  tmp[*strlen0] = '\0';
  *retval = MMG3D_saveTetgenMesh(*mesh,tmp);
  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveGenericMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEGENERICMESH,mmg3d_savegenericmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0, retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveGenericMesh(*mesh,*sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADSOL,mmg3d_loadsol,
             (MMG5_pMesh *mesh,MMG5_pSol *met,char* filename, int *strlen0,int* retval),
             (mesh,met,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadSol(*mesh,*met,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadAllSols function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADALLSOLS,mmg3d_loadallsols,
             (MMG5_pMesh *mesh,MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_loadAllSols(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVESOL,mmg3d_savesol,
             (MMG5_pMesh *mesh,MMG5_pSol *met,char* filename, int *strlen0,int* retval),
             (mesh,met,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveSol(*mesh,*met,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveAllSols function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEALLSOLS,mmg3d_saveallsols,
             (MMG5_pMesh *mesh,MMG5_pSol *sol,char* filename, int *strlen0,int* retval),
             (mesh,sol,filename,strlen0,retval)){
  char *tmp = NULL;

  MMG5_SAFE_MALLOC(tmp,*strlen0+1,char,return);
  strncpy(tmp,filename,*strlen0);
  tmp[*strlen0] = '\0';

  *retval = MMG3D_saveAllSols(*mesh,sol,tmp);

  MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_switch_metricStorage function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SWITCH_METRIC_STORAGE,mmg3d_swith_metricstorage,
             (MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),
             (mesh,met,retval)){

  *retval = MMG3D_switch_metricStorage(*mesh,*met);

  return;
}

/**
 * See \ref MMG3D_Clean_isoSurf function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_CLEAN_ISOSURF,mmg3d_clean_isosurf,
             (MMG5_pMesh *mesh, int* retval), (mesh, retval)) {
  *retval = MMG3D_Clean_isoSurf(*mesh);
  return;
}
