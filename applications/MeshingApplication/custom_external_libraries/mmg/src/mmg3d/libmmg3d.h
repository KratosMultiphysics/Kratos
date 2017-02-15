/* ===========================================================================
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
** ===========================================================================
*/

/**
 * \file mmg3d/libmmg3d.h
 * \brief API headers for the mmg3d library
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning To keep the genheader working, don't break line between the enum
 * name and the opening brace (it creates errors under windows)
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete...
 *
 * \include libexamples/mmg3d/adaptation_example0/example0_a/main.c
 * \include libexamples/mmg3d/example0/adaptation_example0_b/main.c
 * \include libexamples/mmg3d/adaptation_example0_fortran/example0_a/main.F90
 * \include libexamples/mmg3d/adaptation_example0_fortran/example0_b/main.F90
 * \include libexamples/mmg3d/adaptation_example1/main.c
 * \include libexamples/mmg3d/adaptation_example2/main.c
 * \include libexamples/mmg3d/IsosurfDiscretization_example0/main.c
 * \include libexamples/mmg3d/LagrangianMotion_example0/main.c
 */

#ifndef _MMG3DLIB_H
#define _MMG3DLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include "libmmgtypes.h"

/**
 * Maximum array size when storing adjacent points (or ball) of a vertex.
 */
#define MMG3D_LMAX      10240

/**
 * \enum MMG3D_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMG3D_IPARAM asked for integers values ans options prefixed by \a
 * MMG3D_DPARAM asked for real values.
 *
 */
enum MMG3D_Param {
  MMG3D_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  MMG3D_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  MMG3D_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMG3D_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMG3D_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  MMG3D_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  MMG3D_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  MMG3D_IPARAM_optimLES,          /*!< [1/0], Strong mesh optimization for Les computations */
  MMG3D_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  MMG3D_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMG3D_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  MMG3D_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  MMG3D_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  MMG3D_IPARAM_renum,             /*!< [1/0], Turn on/off point relocation with Scotch */
  MMG3D_IPARAM_octree,            /*!< [n], Specify the max number of points per octree cell (DELAUNAY) */
  MMG3D_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  MMG3D_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  MMG3D_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  MMG3D_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMG3D_DPARAM_hgrad,             /*!< [val], Control gradation */
  MMG3D_DPARAM_ls,                /*!< [val], Value of level-set */
  MMG3D_PARAM_size,               /*!< [n], Number of parameters */
};

/*--------------------------- functions header ---------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \param starter dummy argument used to initialize the variadic argument
 * list
 * \param ... variadic arguments that depend to the library function that you
 * want to call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet
 * MMG5_ARG_ppMet, &your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * Here,\a your_mesh is a \a MMG5_pMesh, \a your_metric \a your_level_set and
 * \a your_displacement are \a MMG5_pSol.
 *
 * MMG structures allocation and initialization.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 */
void MMG3D_Init_mesh(const int starter,...);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialize file names to their default values.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_INIT_FILENAMES(mesh,sol)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
void  MMG3D_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_INIT_PARAMETERS(mesh)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
void  MMG3D_Init_parameters(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_INPUTMESHNAME(mesh,meshin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_inputMeshName(MMG5_pMesh mesh,const char* meshin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_OUTPUTMESHNAME(mesh,meshout,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshout\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_INPUTSOLNAME(mesh,sol,solin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol,
                            const char* solin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_OUTPUTSOLNAME(mesh,sol,solout,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solout\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol,
                             const char* solout);

/* init structure sizes */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity type of solutions entities (vertices, triangles...).
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial...).
 * \return 0 if failed, 1 otherwise.
 *
 * Set the solution number, dimension and type.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: typEntity,np,typSol\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity,
                       int np, int typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nprism number of prisms.
 * \param nt number of triangles.
 * \param nquad number of quads.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, prisms, triangles, quadrilaterals and
 * edges of the mesh and allocate the associated tables. If call twice, reset
 * the whole mesh to realloc it at the new size
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_MESHSIZE(mesh,np,ne,nprism,nt,nquad,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER                       :: np,ne,nprism,nt,nquad,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_meshSize(MMG5_pMesh mesh,int np,int ne,int nprism,
                        int nt,int nquad,int na);

/* init structure datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 coordinate of the point along the first dimension.
 * \param c1 coordinate of the point along the second dimension.
 * \param c2 coordinate of the point along the third dimension.
 * \param ref point reference.
 * \param pos position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1,\a c2 and reference \a ref
 * at position \a pos in mesh structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_VERTEX(mesh,c0,c1,c2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: c0,c1,c2\n
 * >     INTEGER, INTENT(IN)           :: ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     double c2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param vertices table of the points coor.
 * The coordinates of the \f$i^{th}\f$ point are stored in
 * vertices[(i-1)*3]\@3.
 * \param refs table of points references.
 * The ref of the \f$i^th\f$ point is stored in refs[i-1].
 * \return 1.
 *
 * Set vertices coordinates and references in mesh structure
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_VERTICES(mesh,vertices,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    REAL(KIND=8), INTENT(IN) :: vertices(*)\n
 * > !    INTEGER,INTENT(IN)       :: refs(*)\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
 int  MMG3D_Set_vertices(MMG5_pMesh mesh, double *vertices,int *refs);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of tetrahedron.
 * \param v1 second vertex of tetrahedron.
 * \param v2 third vertex of tetrahedron.
 * \param v3 fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 * \param pos tetrahedron position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set tetrahedra of vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref at position \a pos in mesh structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TETRAHEDRON(mesh,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_tetrahedron(MMG5_pMesh mesh, int v0, int v1,
                          int v2, int v3, int ref, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param tetra vertices of the tetras of the mesh
 * Vertices of the \f$i^{th}\f$ tetra are stored in tetra[(i-1)*4]\@4.
 * \param refs table of the tetrahedra references.
 * References of the \f$i^{th}\f$ tetra is stored in refs[i-1].
 * \return 0 if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh tetrahedra.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_TETRAHEDRA(mesh,tetra,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    INTEGER, DIMENSION(*), INTENT(IN) :: tetra,refs\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  int  MMG3D_Set_tetrahedra(MMG5_pMesh mesh, int *tetra,
                            int *refs);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of prism.
 * \param v1 second vertex of prism.
 * \param v2 third vertex of prism.
 * \param v3 fourth vertex of prism.
 * \param v4 fifth vertex of prism.
 * \param v5 sixth vertex of prism.
 * \param ref prism reference.
 * \param pos prism position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set prisms of vertices \a v0, \a v1,\a v2,\a v3,\a v4,\a v5 and reference
 * \a ref at position \a pos in mesh structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_PRISM(mesh,v0,v1,v2,v3,v4,v5,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,v4,v5,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_prism(MMG5_pMesh mesh, int v0, int v1,
                     int v2, int v3, int v4, int v5, int ref, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param prisms vertices of the prisms of the mesh
 * Vertices of the \f$i^{th}\f$ prism are stored in prism[(i-1)*6]\@6.
 * \param refs table of the prisms references.
 * References of the \f$i^{th}\f$ prisms is stored in refs[i-1].
 * \return 0 if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh prisms.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_PRISMS(mesh,prisms,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    INTEGER, DIMENSION(*), INTENT(IN) :: prisms,refs\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  int  MMG3D_Set_prisms(MMG5_pMesh mesh, int *prisms,
                        int *refs);

/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TRIANGLE(mesh,v0,v1,v2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_triangle(MMG5_pMesh mesh, int v0, int v1,
                       int v2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param tria pointer toward the table of the tria vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer toward the table of the triangle references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \return 0 if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh triangles.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_TRIANGLES(mesh,tria,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)    :: mesh\n
 * > !    INTEGER,DIMENSION(*), INTENT(IN) :: tria,refs\n
 * > !    INTEGER, INTENT(OUT)             :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  int  MMG3D_Set_triangles(MMG5_pMesh mesh, int *tria, int *refs);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of quadrilateral.
 * \param v1 second vertex of quadrilateral.
 * \param v2 third vertex of quadrilateral.
 * \param v3 fourth vertex of quadrilateral.
 * \param ref quadrilateral reference.
 * \param pos quadrilateral position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set quadrilateral of vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref
 * at position \a pos in mesh structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_QUADRILATERAL(mesh,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_quadrilateral(MMG5_pMesh mesh, int v0, int v1,
                          int v2, int v3, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param quads pointer toward the table of the quads vertices
 * Vertices of the \f$i^{th}\f$ quad are stored in quads[(i-1)*3]\@3.
 * \param refs pointer toward the table of the quadrilateral references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ quad.
 * \return 0 if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh quadrilaterals.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_QUADRILATERALS(mesh,quads,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)    :: mesh\n
 * > !    INTEGER,DIMENSION(*), INTENT(IN) :: quads,refs\n
 * > !    INTEGER, INTENT(OUT)             :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  int  MMG3D_Set_quadrilaterals(MMG5_pMesh mesh, int *quads, int *refs);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first extremity of the edge.
 * \param v1 second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set edges of extremities \a v0, \a v1 and reference \a ref at
 * position \a pos in mesh structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_EDGE(mesh,v0,v1,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_CORNER(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_corner(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set point \a k as required.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDVERTEX(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_requiredVertex(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k element index.
 * \return 1.
 *
 * Set element \a k as required.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTETRAHEDRON(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_requiredTetrahedron(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param reqIdx table of the indices of the required elements.
 * \param nreq number of required elements
 * \return 1.
 *
 * Set the required Tetra.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTETRAHEDRA(mesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER, INTENT(IN)           :: nreq\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int  MMG3D_Set_requiredTetrahedra(MMG5_pMesh mesh, int *reqIdx, int nreq);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as required.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_requiredTriangle(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param reqIdx table of the indices of the required trias.
 * \param nreq number of required trias
 * \return 1.
 *
 * Set the required triangles
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTRIANGLES(mesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER, INTENT(IN)           :: nreq\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int  MMG3D_Set_requiredTriangles(MMG5_pMesh mesh, int *reqIdx, int nreq);

/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set ridge at edge \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_RIDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_ridge(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set edge \a k as required.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_requiredEdge(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k point index
 * \param n0 x componant of the normal at point \a k.
 * \param n1 y componant of the normal at point \a k.
 * \param n2 z componant of the normal at point \a k.
 *
 * \return 1 if success.
 *
 * Set normals (n0,n1,n2) at point \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_NORMALATVERTEX(mesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     REAL(KIND=8), INTENT(IN)      :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_normalAtVertex(MMG5_pMesh mesh, int k, double n0, double n1,
                              double n2) ;

/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_SCALARSOL(met,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: s\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_scalarSol(MMG5_pSol met, double s,int pos);
/**
 * \param met pointer toward the sol structure.
 * \param s table of the scalar solutions values.
 * s[i-1] is the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar solutions at mesh vertices.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_scalarSols(MMG5_pSol met, double *s);
/**
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in solution
 * structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_VECTORSOL(met,vx,vy,vz,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: vx,vy,vz\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz,
                        int pos);
/**
 * \param met pointer toward the sol structure.
 * \param sols table of the vectorial solutions
 * sols[3*(i-1)]\@3 is the solution at vertex i
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial solutions at mesh vertices
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN)      :: sols\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Set_vectorSols(MMG5_pSol met, double *sols);
/**
 * \param met pointer toward the sol structure.
 * \param m11 value of the tensorial solution at position (1,1) in the tensor
 * \param m12 value of the tensorial solution at position (1,2) in the tensor
 * \param m13 value of the tensorial solution at position (1,3) in the tensor
 * \param m22 value of the tensorial solution at position (2,2) in the tensor
 * \param m23 value of the tensorial solution at position (2,3) in the tensor
 * \param m33 value of the tensorial solution at position (3,3) in the tensor
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TENSORSOL(met,m11,m12,m13,m22,m23,m33,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);
/**
 * \param met pointer toward the sol structure.
 * \param sols table of the tensorial solutions.
 * sols[6*(i-1)]\@6 is the solution at vertex i
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: sols\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Set_tensorSols(MMG5_pSol met, double *sols);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * To mark as ended a mesh given without using the API functions
 * (for example, mesh given by mesh->point[i] = 0 ...). Not recommanded.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_HANDGIVENMESH(mesh)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
void MMG3D_Set_handGivenMesh(MMG5_pMesh mesh);

/* check init */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_CHK_MESHDATA(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/** functions to set parameters */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iparam integer parameter to set (see \a MMG3D_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_IPARAMETER(mesh,sol,iparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: iparam,val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam,
                          int val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMG3D_Param structure).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_DPARAMETER(mesh,sol,dparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: dparam\n
 * >     REAL(KIND=8), INTENT(IN)      :: val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam,
                          double val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param hmin minimal edge size.
 * \param hmax maximal edge size.
 * \param hausd value of the Hausdorff number.
 * \return 0 if failed, 1 otherwise.
 *
 * Set local parameters: set the hausdorff value at \a val for all
 * elements of type \a typ and reference \a ref.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_LOCALPARAMETER(mesh,sol,typ,ref,& \n
 * >                                       hmin,hmax,hausd,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: typ,ref\n
 * >     REAL(KIND=8), INTENT(IN)      :: hmin,hmax,hausd\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ,
                              int ref,double hmin,double hmax,double hausd);

/** recover datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param np pointer toward the number of vertices.
 * \param ne pointer toward the number of tetrahedra.
 * \param nprism pointer toward the number of prisms.
 * \param nt pointer toward the number of triangles.
 * \param nquad pointer toward the number of quads.
 * \param na pointer toward the number of edges.
 * \return 1.
 *
 * Get the number of vertices, tetrahedra, prisms, triangles, quadrilaterals and
 * edges of the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_MESHSIZE(mesh,np,ne,nprism,nt,nquad,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER                       :: np,ne,nprism,nt,nquad,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_meshSize(MMG5_pMesh mesh, int* np, int* ne,int *nprism, int* nt,
                        int* nquad, int* na);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity pointer toward the type of entities to which solutions
 * are applied.
 * \param np pointer toward the number of solutions.
 * \param typSol pointer toward the type of the solutions (scalar, vectorial,
 * ...)
 * \return 1.
 *
 * Get the solution number, dimension and type.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER                       :: typEntity,np,typSol\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity,
                       int* np,int* typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 pointer toward the coordinate of the point along the first
 * dimension.
 * \param c1 pointer toward the coordinate of the point along the second
 * dimension.
 * \param c2 pointer toward the coordinate of the point along the third
 * dimension.
 * \param ref pointer to the point reference.
 * \param isCorner pointer toward the flag saying if point is corner.
 * \param isRequired pointer toward the flag saying if point is required.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1,\a c2 and reference \a ref of next
 * vertex of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_VERTEX(mesh,c0,c1,c2,ref,isCorner,isRequired, &\n
 * >                               retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1,c2\n
 * >     INTEGER                       :: ref,isCorner,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2,
                      int* ref,int* isCorner, int* isRequired);

/**
 * \param mesh pointer toward the mesh structure.
 * \param vertices pointer toward the table of the points coordinates.
 * The coordinates of the \f$i^{th}\f$ point are stored in
 * vertices[(i-1)*3]\@3.
 * \param refs pointer to the table of the point references.
 * The ref of the \f$i^th\f$ point is stored in refs[i-1].
 * \param areCorners pointer toward the table of the flags saying if
 * points are corners.
 * areCorners[i-1]=1 if the \f$i^{th}\f$ point is corner.
 * \param areRequired pointer toward the table of flags saying if points
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ point is required.
 * \return 1.
 *
 * Get the coordinates and references of the mesh vertices.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > ! SUBROUTINE MMG3D_GET_VERTICES(mesh,vertices,refs,areCorners,&\n
 * > !                               areRequired,retval)\n
 * > !   MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !   REAL(KIND=8),DIMENSION(*), INTENT(OUT) :: vertices\n
 * > !   INTEGER, DIMENSION(*)                  :: refs,areCorners,areRequired\n
 * > !   INTEGER, INTENT(OUT)          :: retval\n
 * > ! END SUBROUTINE\n
 *
 */
int  MMG3D_Get_vertices(MMG5_pMesh mesh, double* vertices, int* refs,
                        int* areCorners, int* areRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of tetrahedron.
 * \param v1 pointer toward the second vertex of tetrahedron.
 * \param v2 pointer toward the third vertex of tetrahedron.
 * \param v3 pointer toward the fourth vertex of tetrahedron.
 * \param ref pointer toward the tetrahedron reference.
 * \param isRequired pointer toward the flag saying if tetrahedron is
 *  required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref of
 * next tetra of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TETRAHEDRON(mesh,v0,v1,v2,v3,ref,isRequired,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2,v3\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_tetrahedron(MMG5_pMesh mesh, int* v0, int* v1, int* v2,
                           int* v3,int* ref, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param tetra pointer toward the table of the tetrahedra vertices.
 * Vertices of the \f$i^{th}\f$ tetra are stored in tetra[(i-1)*4]\@4.
 * \param refs pointer toward the table of the tetrahedron references.
 * References of the \f$i^{th}\f$ tetra is stored in refs[i-1].
 * \param areRequired pointer toward the table of the flags saying if the
 *  tetrahedra are required. areRequired[i-1]=1 if the \f$i^{th}\f$ tetra
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh tetrahedra.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_TETRAHEDRA(mesh,tetra,refs,areRequired,&\n
 * > !                                   retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: tetra\n
 * > !    INTEGER, DIMENSION(*)         :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  MMG3D_Get_tetrahedra(MMG5_pMesh mesh, int* tetra,int* refs,
                          int* areRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of prism.
 * \param v1 pointer toward the second vertex of prism.
 * \param v2 pointer toward the third vertex of prism.
 * \param v3 pointer toward the fourth vertex of prism.
 * \param v4 pointer toward the fifth vertex of prism.
 * \param v5 pointer toward the sixth vertex of prism.
 * \param ref pointer toward the prism reference.
 * \param isRequired pointer toward the flag saying if prism is
 *  required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0, \a v1, \a v2, \a v3, \a v4, \a v5 and reference \a ref of
 * next prism of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_PRISM(mesh,v0,v1,v2,v3,v4,v5,ref,isRequired,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2,v3,v4,v5\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_prism(MMG5_pMesh mesh, int* v0, int* v1, int* v2,
                     int* v3,int* v4,int* v5,int* ref, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param prisms pointer toward the table of the prisms vertices.
 * Vertices of the \f$i^{th}\f$ prism are stored in prisms[(i-1)*6]\@6.
 * \param refs pointer toward the table of the prism references.
 * References of the \f$i^{th}\f$ prism is stored in refs[i-1].
 * \param areRequired pointer toward the table of the flags saying if the
 *  prisms are required. areRequired[i-1]=1 if the \f$i^{th}\f$ prism
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh prisms.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_PRISMS(mesh,prisms,refs,areRequired,&\n
 * > !                              retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: prisms\n
 * > !    INTEGER, DIMENSION(*)         :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  MMG3D_Get_prisms(MMG5_pMesh mesh, int* prisms,int* refs,
                          int* areRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of triangle.
 * \param v1 pointer toward the second vertex of triangle.
 * \param v2 pointer toward the third vertex of triangle.
 * \param ref pointer toward the triangle reference.
 * \param isRequired pointer toward the flag saying if triangle is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0,\a v1,\a v2 and reference \a ref of next
 * triangle of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TRIANGLE(mesh,v0,v1,v2,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref,
                       int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param tria pointer toward the table of the triangles vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer toward the table of the triangles references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \param areRequired pointer toward table of the flags saying if triangles
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ tria
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh triangles.
 *
 * \remark Fortran interface: (Commentated in order to allow to pass \%val(0)
 * instead of the refs or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_TRIANGLES(mesh,tria,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: tria\n
 * > !    INTEGER, DIMENSION(*)         :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  MMG3D_Get_triangles(MMG5_pMesh mesh, int* tria, int* refs,
                         int* areRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of quadrilateral.
 * \param v1 pointer toward the second vertex of quadrilateral.
 * \param v2 pointer toward the third vertex of quadrilateral.
 * \param v3 pointer toward the fourth vertex of quadrilateral.
 * \param ref pointer toward the quadrilateral reference.
 * \param isRequired pointer toward the flag saying if quadrilateral is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0,\a v1,\a v2,\a v3 and reference \a ref of next
 * quadrilateral of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_QUADRILATERAL(mesh,v0,v1,v2,v3,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2,v3\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_quadrilateral(MMG5_pMesh mesh, int* v0, int* v1, int* v2,int* v3,
                          int* ref, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param quads pointer toward the table of the quadrilaterals vertices
 * Vertices of the \f$i^{th}\f$ quad are stored in tria[(i-1)*4]\@4.
 * \param refs pointer toward the table of the quadrilaterals references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ quad.
 * \param areRequired pointer toward table of the flags saying if quadrilaterals
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ quad
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh quadrilaterals.
 *
 * \remark Fortran interface: (Commentated in order to allow to pass \%val(0)
 * instead of the refs or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_QUADRILATERALS(mesh,quads,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: quads\n
 * > !    INTEGER, DIMENSION(*)         :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  MMG3D_Get_quadrilaterals(MMG5_pMesh mesh, int* quads, int* refs,
                           int* areRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param e0 pointer toward the first extremity of the edge.
 * \param e1 pointer toward the second  extremity of the edge.
 * \param ref pointer toward the edge reference.
 * \param isRidge pointer toward the flag saying if the edge is ridge.
 * \param isRequired pointer toward the flag saying if the edge is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get extremities \a e0, \a e1 and reference \a ref of next edge of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_EDGE(mesh,e0,e1,ref,isRidge,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)          :: e0,e1\n
 * >     INTEGER                       :: ref,isRidge,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k point index
 * \param n0 x componant of the normal at point \a k.
 * \param n1 y componant of the normal at point \a k.
 * \param n2 z componant of the normal at point \a k.
 *
 * \return 1 if success.
 *
 * Get normals (n0,n1,n2) at point \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_NORMALATVERTEX(mesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     REAL(KIND=8)                  :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_normalAtVertex(MMG5_pMesh mesh, int k, double *n0, double *n1,
                              double *n2) ;

/**
 * \param met pointer toward the sol structure.
 * \param s pointer toward the scalar solution value.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solution \a s of next vertex of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_SCALARSOL(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_scalarSol(MMG5_pSol met, double* s);
/**
 * \param met pointer toward the sol structure.
 * \param s table of the scalar solutions at mesh vertices. s[i-1] is
 * the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solutions at mesh vertices.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_Get_scalarSols(MMG5_pSol met, double* s);
/**
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vectorial solution \f$(v_x,v_y,vz)\f$ of next vertex of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_VECTORSOL(met,vx,vy,vz,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: vx,vy,vz\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);
/**
 * \param met pointer toward the sol structure.
 * \param sols table of the solutions at mesh vertices. sols[3*(i-1)]\@3 is
 * the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vectorial solutions at mesh vertices
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Get_vectorSols(MMG5_pSol met, double* sols);
/**
 * \param met pointer toward the sol structure.
 * \param m11 pointer toward the position (1,1) in the solution tensor.
 * \param m12 pointer toward the position (1,2) in the solution tensor.
 * \param m13 pointer toward the position (1,3) in the solution tensor.
 * \param m22 pointer toward the position (2,2) in the solution tensor.
 * \param m23 pointer toward the position (2,3) in the solution tensor.
 * \param m33 pointer toward the position (3,3) in the solution tensor.
 * \return 0 if failed, 1 otherwise.
 *
 * Get tensorial solution of next vertex of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TENSORSOL(met,m11,m12,m13,m22,m23,m33,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33);
/**
 * \param met pointer toward the sol structure.
 * \param sols table of the solutions at mesh vertices.
 * sols[6*(i-1)]\@6 is the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Get tensorial solutions at mesh vertices.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Get_tensorSols(MMG5_pSol met, double *sols);
/**
 * \param mesh pointer toward the mesh structure.
 * \param iparam integer parameter to set (see \a MMG3D_Param structure).
 * \return The value of integer parameter.
 *
 * Get the value of integer parameter \a iparam.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_IPARAMETER(mesh,iparam,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(IN)           :: iparam\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Get_iparameter(MMG5_pMesh mesh, int iparam);

/* input/output functions */
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADMESH(mesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_loadMesh(MMG5_pMesh mesh,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and sol at MSH file format (.msh extension). We read only
 * low-order points, edges, tria, quad, tetra and prisms.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADMSHMESH(mesh,sol,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_loadMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param filename pointer toward the name of file.

 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEMESH(mesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */

int MMG3D_saveMesh(MMG5_pMesh mesh, const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and sol at MSH file format (.msh extension). Write binary file for
 * .mshb extension and ASCII for .msh one.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEMSHMESH(mesh,sol,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int MMG3D_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load metric field.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADSOL(mesh,met,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_loadSol(MMG5_pMesh mesh,MMG5_pSol met, const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVESOL(mesh,met,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_saveSol(MMG5_pMesh mesh,MMG5_pSol met, const char *filename);

/* deallocations */
/**
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend to the library function that you
 * have call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
void MMG3D_Free_all(const int starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend to the library function that you
 * have call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
void MMG3D_Free_structures(const int starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend to the library function that you
 * have call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \a MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
void MMG3D_Free_names(const int starter,...);

/* library */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (metric) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the remesh library.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DLIB(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol met );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (level-set) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the level-set discretization library.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DLS(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_mmg3dls(MMG5_pMesh mesh, MMG5_pSol met );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (output metric) structure.
 * \param disp pointer toward a sol (displacement for the lagrangian motion
 * mode) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the rigidbody movement library.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DMOV(mesh,met,disp,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,disp\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_mmg3dmov(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp );

/** Tools for the library */
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_DEFAULTVALUES(mesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
void MMG3D_defaultValues(MMG5_pMesh mesh);

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Store command line arguments.
 *
 * \remark no matching fortran function.
 *
 */
int  MMG3D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Read local parameters file. This file must have the same name as
 * the mesh with the \a .mmg3d extension or must be named \a
 * DEFAULT.mmg3d.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_PARSOP(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_parsop(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param prog pointer toward the program name.
 *
 * Print help for mmg3d options.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_USAGE(prog,strlen)\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: prog\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >   END SUBROUTINE\n
 *
 */
void  MMG3D_usage(char *prog);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \return 1.
 *
 * Store the info structure in the mesh structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_STOCKOPTIONS(mesh,info,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,info\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 *
 * Recover the info structure stored in the mesh structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_DESTOCKOPTIONS(mesh,info)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,info\n
 * >   END SUBROUTINE\n
 *
 */
void  MMG3D_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/** Checks */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab table of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the _MMG5_defsiz call), 1 for special storage (after this call).
 *
 * Search invalid elements (in term of quality or edge length).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DCHECK(mesh,met,critmin,lmin,lmax,eltab,&\n
 * >                               metridtyp,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)      :: mesh,met\n
 * >     REAL(KIND=8), INTENT(IN)            :: critmin,lmin,lmax\n
 * >     INTEGER,DIMENSION(*), INTENT(OUT)   :: eltab\n
 * >     INTEGER, INTENT(IN)                 :: metridtyp\n
 * >     INTEGER, INTENT(OUT)                :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,double critmin,
                    double lmin, double lmax, int *eltab,char metRidTyp);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param critmin minimum quality for elements.
 * \param eltab pointer toward the table of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the _MMG5_defsiz call), 1 for special storage (after this call).
 *
 * Store elements which have worse quality than \a critmin in \a eltab,
 * \a eltab is allocated and could contain \a mesh->ne elements.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SEARCHQUA(mesh,met,critmin,eltab,metridtyp)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)      :: mesh,met\n
 * >     REAL(KIND=8), INTENT(IN)            :: critmin\n
 * >     INTEGER,DIMENSION(*), INTENT(OUT)   :: eltab\n
 * >     INTEGER, INTENT(IN)                 :: metridtyp\n
 * >   END SUBROUTINE\n
 *
 */
void  MMG3D_searchqua(MMG5_pMesh mesh, MMG5_pSol met, double critmin,
                      int *eltab,char metRidTyp);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab table of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the _MMG5_defsiz call), 1 for special storage (after this call).
 *
 * \return 1 if success, 0 otherwise.
 *
 * Store in \a eltab elements which have edge lengths shorter than \a lmin
 * or longer than \a lmax, \a eltab is allocated and could contain
 * \a mesh->ne elements.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SEARCHLEN(mesh,met,lmin,lmax,eltab,metridtyp,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)      :: mesh,met\n
 * >     REAL(KIND=8), INTENT(IN)            :: lmin,lmax\n
 * >     INTEGER,DIMENSION(*), INTENT(OUT)   :: eltab\n
 * >     INTEGER, INTENT(IN)                 :: metridtyp\n
 * >     INTEGER, INTENT(OUT)                :: retval\n
 * >   END SUBROUTINE\n
 *
 *
 */
int  MMG3D_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin,
                     double lmax,int *eltab,char  metRidTyp);

/** Utils */
/**
 * \brief Return adjacent elements of a tetrahedron.
 * \param mesh pointer toward the mesh structure.
 * \param kel tetrahedron index.
 * \param listet pointer toward the table of the 4 tetra adjacent to \a kel.
 * (the index is 0 if there is no adjacent)
 * \return 1.
 *
 * Find the indices of the 4 adjacent elements of tetrahedron \a
 * kel. \f$listet[i] = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_ADJATET(mesh,kel,listet,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh\n
 * >     INTEGER, INTENT(IN)                :: kel\n
 * >     INTEGER, DIMENSION(4), INTENT(OUT) :: listet\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Get_adjaTet(MMG5_pMesh mesh,int kel, int listet[4]);
/**
 * \param ca pointer toward the coordinates of the first edge's extremity.
 * \param cb pointer toward the coordinates of the second edge's extremity.
 * \param ma pointer toward the metric associated to the first edge's
 * extremity.
 * \param mb pointer toward the metric associated to the second edge's
 * extremity.
 * \return edge length.
 *
 * Compute length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge extremities) according to the size
 * prescription.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LENEDGCOOR(ca,cb,sa,sb,retval)\n
 * >     REAL(KIND=8), INTENT(IN)           :: ca,cb,sa,sb\n
 * >     REAL(KIND=8), INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
double (*MMG3D_lenedgCoor)(double *ca,double *cb,double *sa,double *sb);

/**
 * \param mesh pointer toward the mesh structure.
 * \param pack we pack the mesh at function begining if \f$pack=1\f$.
 * \return 0 if failed, 1 otherwise.
 *
 * Create table of adjacency. Set pack variable to 0 for a compact
 * mesh and to 1 for a mesh that need to be packed.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_HASHTETRA(mesh,pack,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh\n
 * >     INTEGER, INTENT(IN)                :: pack\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  MMG3D_hashTetra(MMG5_pMesh mesh, int pack);

/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the sol structure
 * \return 1 if success
 *
 * Compute isotropic size map according to the mean of the length of the
 * edges passing through a point.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_DOSOL(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_doSol(MMG5_pMesh mesh,MMG5_pSol met);

/** To associate function pointers without calling MMG3D_mmg3dlib */
/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the sol structure (unused).
 *
 * Set function pointers for caltet, lenedg, lenedgCoor defsiz, gradsiz...
 * depending if the readed metric is anisotropic or isotropic
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SETFUNC(mesh,met)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >   END SUBROUTINE\n
 *
 */
void  MMG3D_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \param mesh pointer toward the mesh structure.
 * \param ktri index of the boundary triangle.
 * \param ktet pointer toward the index of the tetra (filled by the function)
 * \param iface pointer toward the index of the face of the tetra \a ktet
 * that correspond to the boundary tria \a ktri.
 * \return 0 if fail, 1 otherwise
 *
 * Fill \a ktet by the index of the tetra to which belong a boundary triangle
 * and \a iface by the index ofthe face of the tetra that correspond to the
 * triangle.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_TETFROMTRIA(mesh,ktri,ktet,iface,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)              :: mesh\n
 * >     INTEGER, INTENT(IN)                      :: ktri\n
 * >     INTEGER, INTENT(OUT)                     :: ktet,iface,retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG3D_Get_tetFromTria(MMG5_pMesh mesh, int ktri, int *ktet, int *iface);


/** Old API ******************************************************** */

/**
 * \enum MMG5_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMG5_IPARAM asked for integers values ans options prefixed by \a
 * MMG5_DPARAM asked for real values.
 *
 */
enum MMG5_Param {
  MMG5_IPARAM_verbose,
  MMG5_IPARAM_mem,
  MMG5_IPARAM_debug,
  MMG5_IPARAM_angle,
  MMG5_IPARAM_iso,
  MMG5_IPARAM_lag,
  MMG5_IPARAM_optim,
  MMG5_IPARAM_noinsert,
  MMG5_IPARAM_noswap,
  MMG5_IPARAM_nomove,
  MMG5_IPARAM_nosurf,
  MMG5_IPARAM_numberOfLocalParam,
  MMG5_IPARAM_renum,
  MMG5_IPARAM_octree,
  MMG5_DPARAM_angleDetection,
  MMG5_DPARAM_hmin,
  MMG5_DPARAM_hmax,
  MMG5_DPARAM_hausd,
  MMG5_DPARAM_hgrad,
  MMG5_DPARAM_ls,
  MMG5_PARAM_size,
};

/*----------------------------- functions header --------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \param mesh adress of a pointer toward a pointer toward the mesh structure
 * \param sol adress of a pointer toward a sol structure
 * (metric or level-set)
 * \param disp adress of a pointer toward a sol structure
 * (displacement for the lagrangian mode).
 *
 * Allocate the mesh and solution structures and initialize it to
 * their default values.
 *
 */
void  MMG5_Init_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol, MMG5_pSol *disp );

/* init structure sizes */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity type of solutions entities (vertices, triangles...).
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial...).
 * \return 0 if failed, 1 otherwise.
 *
 * Set the solution number, dimension and type.
 *
 */
  int  MMG5_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity,int np,
                        int typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of elements (tetrahedra).
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, triangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
 *
 */
int  MMG5_Set_meshSize(MMG5_pMesh mesh, int np, int ne, int nt, int na);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 */
void  MMG5_Init_parameters(MMG5_pMesh mesh);

/* init structure datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 coordinate of the point along the first dimension.
 * \param c1 coordinate of the point along the second dimension.
 * \param c2 coordinate of the point along the third dimension.
 * \param ref point reference.
 * \param pos position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1,\a c2 and reference \a ref
 * at position \a pos in mesh structure
 *
 */
int  MMG5_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     double c2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of tetrahedron.
 * \param v1 second vertex of tetrahedron.
 * \param v2 third vertex of tetrahedron.
 * \param v3 fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 * \param pos tetrahedron position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set tetrahedra of vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref at position \a pos in mesh structure.
 *
 */
int  MMG5_Set_tetrahedron(MMG5_pMesh mesh, int v0, int v1,
                          int v2, int v3, int ref, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure.
 *
 */
int  MMG5_Set_triangle(MMG5_pMesh mesh, int v0, int v1,
                       int v2, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 first extremity of the edge.
 * \param v1 second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set edges of extremities \a v0, \a v1 and reference \a ref at
 * position \a pos in mesh structure
 *
 */
int  MMG5_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref,int pos);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 */
int  MMG5_Set_corner(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set point \a k as required.
 *
 */
int  MMG5_Set_requiredVertex(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k element index.
 * \return 1.
 *
 * Set element \a k as required.
 *
 */
int  MMG5_Set_requiredTetrahedron(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as required.
 *
 */
int  MMG5_Set_requiredTriangle(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set ridge at edge \a k.
 *
 */
int  MMG5_Set_ridge(MMG5_pMesh mesh, int k);
/**
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set edge \a k as required.
 *
 */
int  MMG5_Set_requiredEdge(MMG5_pMesh mesh, int k);

/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 *
 */
int  MMG5_Set_scalarSol(MMG5_pSol met, double s,int pos);
/**
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in solution
 * structure.
 *
 */
int MMG5_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz,
                       int pos);
/**
 * \param met pointer toward the sol structure.
 * \param m11 value of the tensorial solution at position (1,1) in the tensor
 * \param m12 value of the tensorial solution at position (1,2) in the tensor
 * \param m13 value of the tensorial solution at position (1,3) in the tensor
 * \param m22 value of the tensorial solution at position (2,2) in the tensor
 * \param m23 value of the tensorial solution at position (2,3) in the tensor
 * \param m33 value of the tensorial solution at position (3,3) in the tensor
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure.
 *
 */
int MMG5_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);
/**
 * \param mesh pointer toward the mesh structure.
 *
 * To mark as ended a mesh given without using the API functions
 * (for example, mesh given by mesh->point[i] = 0 ...). Not recommanded.
 *
 */
void MMG5_Set_handGivenMesh(MMG5_pMesh mesh);

/* check init */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
 *
 */
int MMG5_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/** functions to set parameters */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iparam integer parameter to set (see \a MMG5_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 */
int  MMG5_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam, int val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param dparam double parameter to set (see \a MMG5_Param structure).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 */
int  MMG5_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam,
                         double val);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param val value of the Hausdorff number.
 * \return 0 if failed, 1 otherwise.
 *
 * Set local parameters: set the hausdorff value at \a val for all
 * elements of type \a typ and reference \a ref.
 *
 */
int  MMG5_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ,
                             int ref, double val);

/** recover datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param np pointer toward the number of vertices.
 * \param ne pointer toward the number of elements (tetrahedra).
 * \param nt pointer toward the number of triangles.
 * \param na pointer toward the number of edges.
 * \return 1.
 *
 * Get the number of vertices, tetrahedra, triangles and edges of the mesh.
 *
 */
int  MMG5_Get_meshSize(MMG5_pMesh mesh, int* np, int* ne, int* nt, int* na);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity pointer toward the type of entities to which solutions
 * are applied.
 * \param np pointer toward the number of solutions.
 * \param typSol pointer toward the type of the solutions (scalar,
 * vectorial...)
 * \return 1.
 *
 * Get the solution number, dimension and type.
 *
 */
int  MMG5_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity,int* np,
                      int* typSol);
/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 pointer toward the coordinate of the point along the first
 * dimension.
 * \param c1 pointer toward the coordinate of the point along the second
 * dimension.
 * \param c2 pointer toward the coordinate of the point along the third
 * dimension.
 * \param ref pointer to the point reference.
 * \param isCorner pointer toward the flag saying if point is corner.
 * \param isRequired pointer toward the flag saying if point is required.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1,\a c2 and reference \a ref of next
 * vertex of mesh.
 *
 */
int  MMG5_Get_vertex(MMG5_pMesh mesh, double* c0,double* c1,double* c2,
                     int* ref,int* isCorner, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of tetrahedron.
 * \param v1 pointer toward the second vertex of tetrahedron.
 * \param v2 pointer toward the third vertex of tetrahedron.
 * \param v3 pointer toward the fourth vertex of tetrahedron.
 * \param ref pointer toward the tetrahedron reference.
 * \param isRequired pointer toward the flag saying if tetrahedron is
 * required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref of
 * next tetra of mesh.
 *
 */
int  MMG5_Get_tetrahedron(MMG5_pMesh mesh, int* v0, int* v1, int* v2,
                          int* v3, int* ref, int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param v0 pointer toward the first vertex of triangle.
 * \param v1 pointer toward the second vertex of triangle.
 * \param v2 pointer toward the third vertex of triangle.
 * \param ref pointer toward the triangle reference.
 * \param isRequired pointer toward the flag saying if triangle is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0,\a v1,\a v2 and reference \a ref of next
 * triangle of mesh.
 *
 */
int  MMG5_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref,
                       int* isRequired);
/**
 * \param mesh pointer toward the mesh structure.
 * \param e0 pointer toward the first extremity of the edge.
 * \param e1 pointer toward the second  extremity of the edge.
 * \param ref pointer toward the edge reference.
 * \param isRidge pointer toward the flag saying if the edge is ridge.
 * \param isRequired pointer toward the flag saying if the edge is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get extremities \a e0, \a e1 and reference \a ref of next edge of mesh.
 *
 */
int  MMG5_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired);
/**
 * \param met pointer toward the sol structure.
 * \param s pointer toward the scalar solution value.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solution \a s of next vertex of mesh.
 *
 */
int  MMG5_Get_scalarSol(MMG5_pSol met, double* s);
/**
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vectorial solution \f$(v_x,v_y,vz)\f$ of next vertex of mesh.
 *
 */
int MMG5_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);
/**
 * \param met pointer toward the sol structure.
 * \param m11 pointer toward the position (1,1) in the solution tensor.
 * \param m12 pointer toward the position (1,2) in the solution tensor.
 * \param m13 pointer toward the position (1,3) in the solution tensor.
 * \param m22 pointer toward the position (2,2) in the solution tensor.
 * \param m23 pointer toward the position (2,3) in the solution tensor.
 * \param m33 pointer toward the position (3,3) in the solution tensor.
 * \return 0 if failed, 1 otherwise.
 *
 * Get tensorial solution of next vertex of mesh.
 *
 */
int MMG5_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33);
/**
 * \param mesh pointer toward the mesh structure.
 * \param iparam integer parameter to set (see \a MMG5_Param structure).
 * \return The value of integer parameter.
 *
 * Get the value of integer parameter \a iparam.
 *
 */
int MMG5_Get_iparameter(MMG5_pMesh mesh, int iparam);

/* input/output functions */
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
 *
 */
int  MMG5_loadMesh(MMG5_pMesh mesh);
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 */
int  MMG5_saveMesh(MMG5_pMesh mesh);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Load metric field.
 *
 */
int  MMG5_loadMet(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 */
int  MMG5_saveMet(MMG5_pMesh mesh, MMG5_pSol met);

/* deallocations */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement for the
 * lagrangian mode).
 *
 * Deallocations before return.
 *
 */
void MMG5_Free_all(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement).
 *
 * Structure deallocations before return.
 *
 */
void MMG5_Free_structures(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp );

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric or solution).
 * \param disp pointer toward a sol structure (displacement).
 *
 * File name deallocations before return.
 *
 */
void MMG5_Free_names(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp);

/* library */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (metric or level-set) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the remesh library.
 *
 */
int  MMG5_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol met );


/* Library tools*/
/**
 * \brief Return adjacent elements of a tetrahedron.
 * \param mesh pointer toward the mesh structure.
 * \param kel tetrahedron index.
 * \param v0 pointer toward the index of the adjacent element of \a kel
 * through its face number 0.
 * \param v1 pointer toward the index of the adjacent element of \a kel
 * through its face number 1.
 * \param v2 pointer toward the index of the adjacent element of \a kel
 * through its face number 2.
 * \param v3 pointer toward the index of the adjacent element of \a kel
 * through its face number 3.
 * \return 1.
 *
 * Find the indices of the 4 adjacent elements of tetrahedron \a
 * kel. \f$v_i = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 */
int MMG5_Get_adjaTet(MMG5_pMesh mesh, int kel, int *v0, int *v1, int *v2,
                     int *v3);
/**
 * \param prog pointer toward the program name.
 *
 * Print help for mmg3d options.
 *
 */
void MMG5_usage(char *prog);
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void MMG5_defaultValues(MMG5_pMesh mesh);

int MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met);

int MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met);
/**
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \return 1.
 *
 * Store the info structure in the mesh structure.
 *
 */
int MMG5_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);

int MMG5_mmg3dcheck(MMG5_pMesh,MMG5_pSol,double, double,double, int *,char);

void MMG5_searchqua(MMG5_pMesh,MMG5_pSol,double, int *, char);

int MMG5_searchlen(MMG5_pMesh,MMG5_pSol, double, double, int *,char);

#ifdef __cplusplus
}
#endif

#endif
