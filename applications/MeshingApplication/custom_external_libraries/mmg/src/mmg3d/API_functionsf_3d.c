
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

#include "mmg3d.h"

/**
 * See \ref MMG3D_Init_mesh function in common/libmmgcommon.h file.
 */
FORTRAN_VARIADIC ( MMG3D_INIT_MESH, mmg3d_init_mesh,
                 (const int starter, ... ),
                 va_list argptr;

                 va_start(argptr, starter);

                 _MMG3D_Init_mesh_var(argptr);

                 va_end(argptr);

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
 * See \ref MMG3D_Set_inputMeshName function in \ref common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG3D_SET_INPUTMESHNAME, mmg3d_set_inputmeshname,
             (MMG5_pMesh *mesh, char* meshin, int *strlen, int* retval),
             (mesh,meshin,strlen,retval)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG3D_Set_inputMeshName(*mesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_inputSolName function in \ref common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG3D_SET_INPUTSOLNAME, mmg3d_set_inputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solin, int* strlen, int* retval),
             (mesh,sol,solin,strlen,retval)) {

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solin,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG3D_Set_inputSolName(*mesh,*sol,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_outputMeshName function in mmgs/libmmgs.h or
 * mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_OUTPUTMESHNAME,mmg3d_set_outputmeshname,
             (MMG5_pMesh *mesh, char* meshout, int* strlen,int* retval),
             (mesh,meshout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshout,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG3D_Set_outputMeshName(*mesh, tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_outputSolName function in \ref common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG3D_SET_OUTPUTSOLNAME,mmg3d_set_outputsolname,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, char* solout,int* strlen, int* retval),
             (mesh,sol,solout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solout,*strlen);
  tmp[*strlen] = '\0';
  *retval = MMG3D_Set_outputSolName(*mesh,*sol,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_Set_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_SOLSIZE,mmg3d_set_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity,
              int* np, int* typSol, int* retval),
             (mesh, sol, typEntity, np, typSol, retval)) {
  *retval = MMG3D_Set_solSize(*mesh,*sol,*typEntity,*np,*typSol);
  return;
}

/**
 * See \ref MMG3D_Set_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_MESHSIZE,mmg3d_set_meshsize,
             (MMG5_pMesh *mesh, int *np, int *ne, int *nprism,
              int *nt, int *nquad, int *na, int *retval),
             (mesh,np,ne,nprism,nt,nquad,na,retval)) {
  *retval = MMG3D_Set_meshSize(*mesh,*np,*ne,*nprism,*nt,*nquad,*na);
  return;
}

/**
 * See \ref MMG3D_Get_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_SOLSIZE,mmg3d_get_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity, int* np, int* typSol, int* retval),
             (mesh,sol,typEntity,np,typSol,retval)) {

  *retval = MMG3D_Get_solSize(*mesh,*sol,typEntity,np,typSol);
  return;
}

/**
 * See \ref MMG3D_Get_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_MESHSIZE,mmg3d_get_meshsize,
             (MMG5_pMesh *mesh, int* np, int* ne,int *nprism,
              int* nt,int *nquad, int* na, int* retval),
             (mesh,np,ne,nprism,nt,nquad,na,retval)) {

  *retval = MMG3D_Get_meshSize(*mesh,np,ne,nprism,nt,nquad,na);
  return;
}

/**
 * See \ref MMG3D_Set_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_VERTEX,mmg3d_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* pos, int* retval),
             (mesh,c0,c1,c2,ref,pos,retval)) {

  *retval = MMG3D_Set_vertex(*mesh,*c0,*c1,*c2,*ref,*pos);
  return;
}
/**
 * See \ref MMG3D_Get_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_VERTEX,mmg3d_get_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* isCorner, int* isRequired, int* retval),
             (mesh,c0,c1,c2,ref,isCorner,isRequired, retval)) {
  *retval = MMG3D_Get_vertex(*mesh,c0,c1,c2,ref,isCorner,isRequired);
  return;
}

/**
 * See \ref MMG3D_Set_vertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_VERTICES,mmg3d_set_vertices,
             (MMG5_pMesh *mesh, double* vertices, int* refs, int* retval),
             (mesh,vertices,refs,retval)) {

  *retval = MMG3D_Set_vertices(*mesh,vertices,refs);
  return;
}


/**
 * See \ref MMG3D_Get_vertices function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_VERTICES,mmg3d_get_vertices,
             (MMG5_pMesh *mesh, double* vertices, int* refs,
              int* areCorners, int* areRequired, int* retval),
             (mesh,vertices,refs,areCorners,areRequired, retval)) {
  *retval = MMG3D_Get_vertices(*mesh,vertices,refs,areCorners,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TETRAHEDRON,mmg3d_set_tetrahedron,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *v2, int *v3, int *ref,
              int *pos, int* retval),
             (mesh,v0,v1,v2,v3,ref,pos,retval)){
  *retval = MMG3D_Set_tetrahedron(*mesh,*v0,*v1,*v2,*v3,*ref,*pos);
  return;
}

/**
 * See \ref MMG3D_Set_tetrahedra function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TETRAHEDRA,mmg3d_set_tetrahedra,
             (MMG5_pMesh *mesh, int *tetra, int *refs, int* retval),
             (mesh,tetra,refs,retval)){
  *retval = MMG3D_Set_tetrahedra(*mesh,tetra,refs);
  return;
}

/**
 * See \ref MMG3D_Get_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TETRAHEDRON,mmg3d_get_tetrahedron,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int *v3,
              int* ref, int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = MMG3D_Get_tetrahedron(*mesh,v0,v1,v2,v3,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Get_tetrahedra function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TETRAHEDRA,mmg3d_get_tetrahedra,
             (MMG5_pMesh *mesh, int* tetra, int* refs, int* areRequired,
              int* retval),
             (mesh,tetra,refs,areRequired,retval)) {
  *retval = MMG3D_Get_tetrahedra(*mesh,tetra,refs,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_prism function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_PRISM,mmg3d_set_prism,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *v2, int *v3,
              int *v4, int *v5,int *ref,int *pos, int* retval),
             (mesh,v0,v1,v2,v3,v4,v5,ref,pos,retval)){
  *retval = MMG3D_Set_prism(*mesh,*v0,*v1,*v2,*v3,*v4,*v5,*ref,*pos);
  return;
}

/**
 * See \ref MMG3D_Set_prisms function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_PRISMS,mmg3d_set_prisms,
             (MMG5_pMesh *mesh, int *prisms, int *refs, int* retval),
             (mesh,prisms,refs,retval)){
  *retval = MMG3D_Set_prisms(*mesh,prisms,refs);
  return;
}

/**
 * See \ref MMG3D_Get_prism function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_PRISM,mmg3d_get_prism,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int *v3,
              int *v4, int* v5,int* ref, int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,v4,v5,ref,isRequired,retval)) {
  *retval = MMG3D_Get_prism(*mesh,v0,v1,v2,v3,v4,v5,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Get_prisms function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_PRISMS,mmg3d_get_prisms,
             (MMG5_pMesh *mesh, int* prisms, int* refs, int* areRequired,
              int* retval),
             (mesh,prisms,refs,areRequired,retval)) {
  *retval = MMG3D_Get_prisms(*mesh,prisms,refs,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TRIANGLE,mmg3d_set_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (mesh,v0,v1,v2,ref,pos,retval)) {
  *retval = MMG3D_Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
  return;
}

/**
 * See \ref MMG3D_Get_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TRIANGLE,mmg3d_get_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref
              ,int* isRequired, int* retval),
             (mesh,v0,v1,v2,ref,isRequired,retval)) {
  *retval = MMG3D_Get_triangle(*mesh,v0,v1,v2,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Set_triangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_TRIANGLES,mmg3d_set_triangles,
             (MMG5_pMesh *mesh, int* tria, int* refs,
              int* retval),
             (mesh,tria,refs,retval)) {
  *retval = MMG3D_Set_triangles(*mesh, tria, refs);
  return;
}

/**
 * See \ref MMG3D_Get_triangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_TRIANGLES,mmg3d_get_triangles,
             (MMG5_pMesh *mesh, int* tria, int* refs,int* areRequired,
              int* retval),
             (mesh,tria,refs,areRequired,retval)) {
  *retval = MMG3D_Get_triangles(*mesh,tria,refs,areRequired);
  return;
}
/**
 * See \ref MMG3D_Set_quadrilateral function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_QUADRILATERAL,mmg3d_set_quadrilateral,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2,int *v3,
              int* ref,int* pos,int* retval),
             (mesh,v0,v1,v2,v3,ref,pos,retval)) {
  *retval = MMG3D_Set_quadrilateral(*mesh, *v0, *v1, *v2, *v3,*ref, *pos);
  return;
}

/**
 * See \ref MMG3D_Get_quadrilateral function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_QUADRILATERAL,mmg3d_get_quadrilateral,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2,int *v3,
               int* ref,int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = MMG3D_Get_quadrilateral(*mesh,v0,v1,v2,v3,ref,isRequired);
  return;
}
/**
 * See \ref MMG3D_Set_quadrilaterals function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_QUADRILATERALS,mmg3d_set_quadrilaterals,
             (MMG5_pMesh *mesh, int* quads, int* refs,
              int* retval),
             (mesh,quads,refs,retval)) {
  *retval = MMG3D_Set_quadrilaterals(*mesh, quads, refs);
  return;
}

/**
 * See \ref MMG3D_Get_quadrilaterals function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_QUADRILATERALS,mmg3d_get_quadrilaterals,
             (MMG5_pMesh *mesh, int* quads, int* refs,int* areRequired,
              int* retval),
             (mesh,quads,refs,areRequired,retval)) {
  *retval = MMG3D_Get_quadrilaterals(*mesh,quads,refs,areRequired);
  return;
}

/**
 * See \ref MMG3D_Set_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_EDGE,mmg3d_set_edge,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *ref, int *pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)){
  *retval = MMG3D_Set_edge(*mesh,*v0,*v1,*ref,*pos);
  return;
}

/**
 * See \ref MMG3D_Get_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_EDGE,mmg3d_get_edge,(MMG5_pMesh *mesh, int* e0, int* e1, int* ref
                                          ,int* isRidge, int* isRequired, int* retval),
             (mesh,e0,e1,ref,isRidge,isRequired,retval)) {
  *retval = MMG3D_Get_edge(*mesh,e0,e1,ref,isRidge,isRequired);
  return;
}

/**
 * See \ref MMG3D_Set_corner function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_CORNER,mmg3d_set_corner,(MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG3D_Set_corner(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDVERTEX,mmg3d_set_requiredvertex,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG3D_Set_requiredVertex(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredTetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTETRAHEDRON,mmg3d_set_requiredtetrahedron,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_requiredTetrahedron(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_requiredTetrahedra function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTETRAHEDRA,mmg3d_set_requiredtetrahedra,
             (MMG5_pMesh *mesh, int *reqIdx, int *nreq, int* retval),
             (mesh,reqIdx,nreq,retval)) {
  *retval = MMG3D_Set_requiredTetrahedra(*mesh,reqIdx, *nreq);
  return;
}

/**
 * See \ref MMG3D_Set_requiredTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTRIANGLE,mmg3d_set_requiredtriangle,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_requiredTriangle(*mesh, *k);
  return;
}
/**
 * See \ref MMG3D_Set_requiredTriangles function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDTRIANGLES,mmg3d_set_requiredtriangles,
             (MMG5_pMesh *mesh, int *reqIdx, int *nreq, int* retval),
             (mesh,reqIdx,nreq,retval)) {
  *retval = MMG3D_Set_requiredTriangles(*mesh, reqIdx, *nreq);
  return;
}

/**
 * See \ref MMG3D_Set_ridge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_RIDGE,mmg3d_set_ridge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_ridge(*mesh,*k);
  return;
}
/**
 * See \ref MMG3D_Set_requiredEdge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_REQUIREDEDGE,mmg3d_set_requirededge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG3D_Set_requiredEdge(*mesh,*k);
  return;
}

/**
 * See \ref MMG3D_Set_normalAtVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_NORMALATVERTEX,mmg3d_set_normalatvertex,
             (MMG5_pMesh *mesh, int *k, double* n0, double* n1, double* n2,int* retval),
             (mesh,k,n0,n1,n2,retval)) {
  *retval = MMG3D_Set_normalAtVertex(*mesh,*k, *n0, *n1, *n2);
  return;
}
/**
 * See \ref MMG3D_Get_normalAtVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_NORMALATVERTEX,mmg3d_get_normalatvertex,
             (MMG5_pMesh *mesh, int *k, double* n0, double* n1, double* n2,int* retval),
             (mesh,k,n0,n1,n2,retval)) {
  *retval = MMG3D_Get_normalAtVertex(*mesh,*k, n0, n1, n2);
  return;
}

/**
 * See \ref MMG3D_Set_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_SCALARSOL,mmg3d_set_scalarsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
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
              int *pos, int* retval),
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
              double* m22,double *m23, double *m33, int *pos, int* retval),
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
 * See \ref MMG3D_Set_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SET_IPARAMETER,mmg3d_set_iparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, int *val, int* retval),
             (mesh,sol,iparam,val,retval)){
  *retval = MMG3D_Set_iparameter(*mesh,*sol,*iparam,*val);
  return;
}

/**
 * See \ref MMG3D_Get_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_GET_IPARAMETER,mmg3d_get_iparameter,
             (MMG5_pMesh *mesh, int *iparam, int* retval),
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
             (MMG5_pMesh *mesh,MMG5_pSol *sol, int *typ, int *ref,
              double *hmin, double *hmax, double *hausd, int* retval),
             (mesh,sol,typ,ref,hmin, hmax, hausd,retval)){
  *retval = MMG3D_Set_localParameter(*mesh,*sol,*typ,*ref,*hmin,*hmax,*hausd);
  return;
}

/**
 * See \ref MMG5_Free_all function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_VARIADIC(MMG3D_FREE_ALL,mmg3d_free_all,
                 (const int starter,...),
                 va_list argptr;

                 va_start(argptr, starter);

                 _MMG3D_Free_all_var(argptr);

                 va_end(argptr);

                 return;
  )

/**
 * See \ref MMG3D_Free_structures function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_VARIADIC(MMG3D_FREE_STRUCTURES,mmg3d_free_structures,
                 (const int starter,...),
                 va_list argptr;

                 va_start(argptr, starter);

                 _MMG3D_Free_structures_var(argptr);

                 va_end(argptr);

                 return;
  )

/**
 * See \ref MMG3D_Free_names function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_VARIADIC(MMG3D_FREE_NAMES,mmg3d_free_names,
             (const int starter,...),
             va_list argptr;

             va_start(argptr, starter);

             _MMG3D_Free_names_var(argptr);

             va_end(argptr);

             return;
  )


/**
 * See \ref MMG3D_loadMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADMESH,mmg3d_loadmesh,
             (MMG5_pMesh *mesh,char* filename, int *strlen,int* retval),
             (mesh,filename,strlen, retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,filename,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG3D_loadMesh(*mesh,tmp);

  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadMshMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADMSHMESH,mmg3d_loadmshmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen,int* retval),
             (mesh,sol,filename,strlen, retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,filename,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG3D_loadMshMesh(*mesh,*sol,tmp);

  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEMESH,mmg3d_savemesh,
             (MMG5_pMesh *mesh,char* filename, int *strlen,int* retval),
             (mesh,filename,strlen, retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,filename,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG3D_saveMesh(*mesh,tmp);

  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveMshMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVEMSHMESH,mmg3d_savemshmesh,
             (MMG5_pMesh *mesh, MMG5_pSol *sol,char* filename, int *strlen,
              int* retval),
             (mesh,sol,filename,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,filename,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG3D_saveMshMesh(*mesh,*sol,tmp);

  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_loadSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_LOADSOL,mmg3d_loadsol,
             (MMG5_pMesh *mesh,MMG5_pSol *met,char* filename, int *strlen,int* retval),
             (mesh,met,filename,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,filename,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG3D_loadSol(*mesh,*met,tmp);

  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref MMG3D_saveSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG3D_SAVESOL,mmg3d_savesol,
             (MMG5_pMesh *mesh,MMG5_pSol *met,char* filename, int *strlen,int* retval),
             (mesh,met,filename,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,filename,*strlen);
  tmp[*strlen] = '\0';

  *retval = MMG3D_saveSol(*mesh,*met,tmp);

  _MMG5_SAFE_FREE(tmp);

  return;
}


/** Old API ======================================================*/

/**
 * See \ref MMG5_Init_mesh function in common/libmmgcommon.h file.
 */
FORTRAN_NAME(MMG5_INIT_MESH, mmg5_init_mesh,(MMG5_pMesh *mesh, MMG5_pSol *sol
                                             ,MMG5_pSol *disp
               ),(mesh,sol,disp
                 )) {
  MMG5_Init_mesh(mesh,sol,disp);
  return;
}

/**
 * See \ref MMG5_Init_parameters function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_INIT_PARAMETERS,mmg5_init_parameters,(MMG5_pMesh *mesh),(mesh)) {
  MMG5_Init_parameters(*mesh);
  return;
}

/**
 * See \ref MMG5_Set_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SOLSIZE,mmg5_set_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity,
              int* np, int* typSol, int* retval),
             (mesh, sol, typEntity, np, typSol, retval)) {
  *retval = MMG5_Set_solSize(*mesh,*sol,*typEntity,*np,*typSol);
  return;
}

/**
 * See \ref MMG5_Set_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_MESHSIZE,mmg5_set_meshsize,
             (MMG5_pMesh *mesh, int *np, int *ne, int *nt, int *na, int *retval),
             (mesh,np,ne,nt,na,retval)) {
  *retval = MMG5_Set_meshSize(*mesh,*np,*ne,*nt,*na);
  return;
}

/**
 * See \ref MMG5_Get_solSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_SOLSIZE,mmg5_get_solsize,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int* typEntity, int* np, int* typSol, int* retval),
             (mesh,sol,typEntity,np,typSol,retval)) {

  *retval = MMG5_Get_solSize(*mesh,*sol,typEntity,np,typSol);
  return;
}

/**
 * See \ref MMG5_Get_meshSize function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_MESHSIZE,mmg5_get_meshsize,
             (MMG5_pMesh *mesh, int* np, int* ne, int* nt, int* na, int* retval),
             (mesh,np,ne,nt, na,retval)) {

  *retval = MMG5_Get_meshSize(*mesh,np,ne,nt,na);
  return;
}

/**
 * See \ref MMG5_Set_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_VERTEX,mmg5_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* pos, int* retval),
             (mesh,c0,c1,c2,ref,pos,retval)) {

  *retval = MMG5_Set_vertex(*mesh,*c0,*c1,*c2,*ref,*pos);
  return;
}

/**
 * See \ref MMG5_Get_vertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_VERTEX,mmg5_get_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* isCorner, int* isRequired, int* retval),
             (mesh,c0,c1,c2,ref,isCorner,isRequired, retval)) {
  *retval = MMG5_Get_vertex(*mesh,c0,c1,c2,ref,isCorner,isRequired);
  return;
}

/**
 * See \ref MMG5_Set_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_TETRAHEDRON,mmg5_set_tetrahedron,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *v2, int *v3, int *ref,
              int *pos, int* retval),
             (mesh,v0,v1,v2,v3,ref,pos,retval)){
  *retval = MMG5_Set_tetrahedron(*mesh,*v0,*v1,*v2,*v3,*ref,*pos);
  return;
}

/**
 * See \ref MMG5_Get_tetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_TETRAHEDRON,mmg5_get_tetrahedron,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* v3,
              int* ref, int* isRequired, int* retval),
             (mesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = MMG5_Get_tetrahedron(*mesh,v0,v1,v2,v3,ref,isRequired);
  return;
}

/**
 * See \ref MMG5_Set_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_TRIANGLE,mmg5_set_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (mesh,v0,v1,v2,ref,pos,retval)) {
  *retval = MMG5_Set_triangle(*mesh, *v0, *v1, *v2, *ref, *pos);
  return;
}

/**
 * See \ref MMG5_Get_triangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_TRIANGLE,mmg5_get_triangle,
             (MMG5_pMesh *mesh, int* v0, int* v1, int* v2, int* ref
              ,int* isRequired, int* retval),
             (mesh,v0,v1,v2,ref,isRequired,retval)) {
  *retval = MMG5_Get_triangle(*mesh,v0,v1,v2,ref,isRequired);
  return;
}

/**
 * See \ref MMG5_Set_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_EDGE,mmg5_set_edge,
             (MMG5_pMesh *mesh, int *v0, int *v1, int *ref, int *pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)){
  *retval = MMG5_Set_edge(*mesh,*v0,*v1,*ref,*pos);
  return;
}

/**
 * See \ref MMG5_Get_edge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_EDGE,mmg5_get_edge,(MMG5_pMesh *mesh, int* e0, int* e1, int* ref
                                          ,int* isRidge, int* isRequired, int* retval),
             (mesh,e0,e1,ref,isRidge,isRequired,retval)) {
  *retval = MMG5_Get_edge(*mesh,e0,e1,ref,isRidge,isRequired);
  return;
}

/**
 * See \ref MMG5_Set_corner function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_CORNER,mmg5_set_corner,(MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG5_Set_corner(*mesh,*k);
  return;
}

/**
 * See \ref MMG5_Set_requiredVertex function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDVERTEX,mmg5_set_requiredvertex,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval =  MMG5_Set_requiredVertex(*mesh,*k);
  return;
}

/**
 * See \ref MMG5_Set_requiredTetrahedron function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDTETRAHEDRON,mmg5_set_requiredtetrahedron,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG5_Set_requiredTetrahedron(*mesh,*k);
  return;
}

/**
 * See \ref MMG5_Set_requiredTriangle function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDTRIANGLE,mmg5_set_requiredtriangle,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG5_Set_requiredTriangle(*mesh, *k);
  return;
}

/**
 * See \ref MMG5_Set_ridge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_RIDGE,mmg5_set_ridge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG5_Set_ridge(*mesh,*k);
  return;
}

/**
 * See \ref MMG5_Set_requiredEdge function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_REQUIREDEDGE,mmg5_set_requirededge,
             (MMG5_pMesh *mesh, int *k, int* retval),
             (mesh,k,retval)) {
  *retval = MMG5_Set_requiredEdge(*mesh,*k);
  return;
}

/**
 * See \ref MMG5_Set_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_SCALARSOL,mmg5_set_scalarsol,
             (MMG5_pSol *met, double *s, int *pos, int* retval),
             (met,s,pos,retval)) {
  *retval = MMG5_Set_scalarSol(*met,*s,*pos);
  return;
}

/**
 * See \ref MMG5_Get_scalarSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_SCALARSOL,mmg5_get_scalarsol,
             (MMG5_pSol *met, double* s, int* retval),
             (met,s,retval)) {
  *retval = MMG5_Get_scalarSol(*met,s);
  return;
}

/**
 * See \ref MMG5_Set_vectorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_VECTORSOL,mmg5_set_vectorsol,
             (MMG5_pSol *met, double *vx, double *vy, double *vz,
              int *pos, int* retval),
             (met,vx,vy,vz,pos,retval)) {
  *retval = MMG5_Set_vectorSol(*met,*vx,*vy,*vz,*pos);
  return;
}

/**
 * See \ref MMG5_Get_vectorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_VECTORSOL,mmg5_get_vectorsol,
             (MMG5_pSol *met, double* vx,double *vy, double *vz, int* retval),
             (met,vx,vy,vz,retval)) {
  *retval = MMG5_Get_vectorSol(*met,vx,vy,vz);
  return;
}

/**
 * See \ref MMG5_Set_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_TENSORSOL,mmg5_set_tensorsol,
             (MMG5_pSol *met, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int *pos, int* retval),
             (met,m11,m12,m13,m22,m23,m33,pos,retval)) {
  *retval = MMG5_Set_tensorSol(*met,*m11,*m12,*m13,*m22,*m23,*m33,*pos);
  return;
}

/**
 * See \ref MMG5_Get_tensorSol function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_TENSORSOL,mmg5_get_tensorsol,
             (MMG5_pSol *met, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int* retval),
             (met,m11,m12,m13,m22,m23,m33,retval)) {
  *retval = MMG5_Get_tensorSol(*met,m11,m12,m13,m22,m23,m33);
  return;
}

/**
 * See \ref MMG5_Set_handGivenMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_HANDGIVENMESH,mmg5_set_handgivenmesh,
             (MMG5_pMesh *mesh),
             (mesh)) {
  MMG5_Set_handGivenMesh(*mesh);
  return;
}

/**
 * See \ref MMG5_Chk_meshData function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_CHK_MESHDATA,mmg5_chk_meshdata,
             (MMG5_pMesh *mesh,MMG5_pSol *met, int* retval),
             (mesh,met,retval)) {
  *retval = MMG5_Chk_meshData(*mesh,*met);
  return;
}

/**
 * See \ref MMG5_Set_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_IPARAMETER,mmg5_set_iparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *iparam, int *val, int* retval),
             (mesh,sol,iparam,val,retval)){
  *retval = MMG5_Set_iparameter(*mesh,*sol,*iparam,*val);
  return;
}

/**
 * See \ref MMG5_Get_iparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_GET_IPARAMETER,mmg5_get_iparameter,
             (MMG5_pMesh *mesh, int *iparam, int* retval),
             (mesh,iparam,retval)){
  *retval = MMG5_Get_iparameter(*mesh,*iparam);
  return;
}

/**
 * See \ref MMG5_Set_dparameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_DPARAMETER,mmg5_set_dparameter,
             (MMG5_pMesh *mesh, MMG5_pSol *sol, int *dparam, double *val, int* retval),
             (mesh,sol,dparam,val,retval)){
  *retval = MMG5_Set_dparameter(*mesh,*sol,*dparam,*val);
  return;
}

/**
 * See \ref MMG5_Set_localParameter function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SET_LOCALPARAMETER,mmg5_set_localparameter,
             (MMG5_pMesh *mesh,MMG5_pSol *sol, int *typ, int *ref, double *val, int* retval),
             (mesh,sol,typ,ref,val,retval)){
  *retval = MMG5_Set_localParameter(*mesh,*sol,*typ,*ref,*val);
  return;
}


/**
 * See \ref MMG5_Free_structures function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_FREE_STRUCTURES,mmg5_free_structures,(MMG5_pMesh *mesh,MMG5_pSol *met
                                                        ,MMG5_pSol *disp
               ),(mesh,met,disp
                 )){
  MMG5_Free_structures(*mesh,*met,(disp==NULL)?NULL:*disp);
  return;
}

/**
 * See \ref MMG5_Free_names function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_FREE_NAMES,mmg5_free_names,(MMG5_pMesh *mesh,MMG5_pSol *met,
                                              MMG5_pSol *disp),
             (mesh,met,disp))
{
  MMG5_Free_names(*mesh,*met,(disp==NULL)?NULL:*disp);
  return;
}

/**
 * See \ref MMG5_loadMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_LOADMESH,mmg5_loadmesh,(MMG5_pMesh *mesh,int* retval),(mesh, retval)){

  *retval = MMG5_loadMesh(*mesh);

  return;
}
/**
 * See \ref MMG5_saveMesh function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SAVEMESH,mmg5_savemesh,(MMG5_pMesh *mesh,int* retval),(mesh, retval)){

  *retval = MMG5_saveMesh(*mesh);

  return;
}

/**
 * See \ref MMG5_loadMet function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_LOADMET,mmg5_loadmet,(MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),(mesh,met,retval)){

  *retval = MMG5_loadMet(*mesh,*met);

  return;
}

/**
 * See \ref MMG5_saveMet function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_NAME(MMG5_SAVEMET,mmg5_savemet,(MMG5_pMesh *mesh,MMG5_pSol *met,int* retval),(mesh,met,retval)){

  *retval = MMG5_saveMet(*mesh,*met);

  return;
}
