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

/*
 * This file defines the C and Fortran headers of the mmgs API, and
 * their Doxygen documentation.
 *
 * NOTES FOR DEVELOPERS:
 *
 * - The Fortran headers are generated from comment lines that start with '* >'.
 *   They must match the C declarations.
 *
 * - We cannot handle enum types in the Fortran version so enums are replaced
 *   by ints in both versions.
 *
 * - To keep the genheader program working, don't break line between an enum
 *   name and the opening brace (it creates errors under windows)
 *
 * - Since Mmg version 5,
 * -- data structures and parameters that are common between mmg3d, mmg2d
 *    and mmgs use the MMG5_ prefix;
 * -- API functions should have an MMG3D_, MMG2D_, or MMGS_ prefix,
 *    depending on the library; and
 * -- some MMG5_ API functions exists but they are common to the
 *    three libraries.
 *
 */

/**
 * \file mmgs/libmmgs.h
 * \ingroup API
 * \brief API headers and documentation for the mmgs library
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * These are the API functions for the mmgs library. These functions allow to
 * load and save meshes and data defined on meshes; add, extract, or modify mesh
 * data; and to call the library functions that perform remeshing and level-set
 * discretization.
 *
 * Meshes are here defined in terms of vertices and two-dimensional objects:
 * triangles and quadrangles, which live in 3D space. Edges can also be
 * represented. All of these \a entities can have a \a reference: an integer
 * value that can serve as a group identifier. In addition mesh entities can
 * have \a attributes such as "required" or "corner".
 *
 * Data defined on meshes can be for example functions that are meant for
 * level-set discretization and metric tensors that will govern edge
 * lengths. These data can be scalar, vector, or (symmetric) tensor-valued; and
 * there can be more than one data item associated with a mesh entity. These
 * data are often referred to as \a solutions.
 *
 * Two of the functions here are referred to as "programs", because they perform
 * the tasks for which mmgs is meant: meshing and level-set discretization.  The
 * other functions merely serve to load and save data and to perform pre- and
 * post-processing. These programs actually behave much like independent
 * programs: they send diagnostic output to stdout and in rare cases they may
 * call the exit() function.
 *
 */

#ifndef MMGSLIB_H
#define MMGSLIB_H


#ifdef __cplusplus
extern "C" {
#endif

#include "mmg/mmgs/mmgs_export.h"
#include "mmg/common/libmmgtypes.h"

/**
 * Maximum array size when storing adjacent vertices (or ball) of a vertex.
 */
#define MMGS_LMAX      1024

/**
 * \enum MMGS_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * MMGS_IPARAM asked for integers values ans options prefixed by \a
 * MMGS_DPARAM asked for real values.
 *
 */
enum MMGS_Param {
  MMGS_IPARAM_verbose,           /*!< [-1..10], Level of verbosity */
  MMGS_IPARAM_mem,               /*!< [n/-1], Max memory size in MBytes or -1 to keep the default value */
  MMGS_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMGS_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMGS_IPARAM_iso,               /*!< [1/0], Enable level-set discretization */
  MMGS_IPARAM_isosurf,           /*!< [1/0], Enable level-set discretization on the surface part */
  MMGS_IPARAM_isoref,            /*!< [0/n], Iso-surface boundary material reference */
  MMGS_IPARAM_keepRef,           /*!< [1/0], Preserve the initial domain references in level-set mode */
  MMGS_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  MMGS_IPARAM_noinsert,          /*!< [1/0], Avoid/allow vertex insertion */
  MMGS_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMGS_IPARAM_nomove,            /*!< [1/0], Avoid/allow vertex relocation */
  MMGS_IPARAM_nreg,              /*!< [0/1], Disable/enable regularization of normals */
  MMGS_IPARAM_xreg,              /*!< [0/1], Disable/enable regularization by moving vertices */
  MMGS_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  MMGS_IPARAM_numberOfLSBaseReferences, /*!< [n], Number of base references for bubble removal */
  MMGS_IPARAM_numberOfMat,              /*!< [n], Number of material in level-set mode */
  MMGS_IPARAM_numsubdomain,      /*!< [0/n], Save only subdomain n (0==all subdomains) */
  MMGS_IPARAM_renum,             /*!< [1/0], Turn on/off renumbering with Scotch */
  MMGS_IPARAM_anisosize,         /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  MMGS_IPARAM_nosizreq,          /*!< [0/1], Allow/avoid overwritings of sizes at required vertices (advanced usage) */
  MMGS_DPARAM_angleDetection,    /*!< [val], Threshold for angle detection */
  MMGS_DPARAM_hmin,              /*!< [val], Minimal edge length */
  MMGS_DPARAM_hmax,              /*!< [val], Maximal edge length */
  MMGS_DPARAM_hsiz,              /*!< [val], Constant edge length */
  MMGS_DPARAM_hausd,             /*!< [val], Global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMGS_DPARAM_hgrad,             /*!< [val], Gradation */
  MMGS_DPARAM_hgradreq,          /*!< [val], Gradation on required entites (advanced usage) */
  MMGS_DPARAM_ls,                /*!< [val], Function value where the level set is to be discretized */
  MMGS_DPARAM_xreg,              /*!< [val], Relaxation parameter for coordinate regularization (0<val<1) */
  MMGS_DPARAM_rmc,               /*!< [-1/val], Remove small disconnected components in level-set mode */
  MMGS_PARAM_size,               /*!< [n], Number of parameters */
};

/*----------------------------- function headers -----------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \brief Initialize a mesh structure and optionally the associated solution and
 * metric structures.
 *
 * \param starter dummy argument used to initialize the variadic argument list
 * \param ... variadic arguments.
 *
 * For the MMGS_mmgslib function, you need
 * to call the \ref MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the \ref MMGS_mmgsls function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here, \a your_mesh is a \ref MMG5_pMesh, \a your_metric and \a your_level_set
 * are \ref MMG5_pSol.
 *
 * \return 1 on success, 0 on failure
 *
 * \remark No fortran interface, to allow variadic arguments.
 *
 */
LIBMMGS_EXPORT int MMGS_Init_mesh(const int starter,...);

/**
 * \brief Initialize file names to their default values.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_INIT_FILENAMES(mesh,sol)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void  MMGS_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);

/**
 * \brief Initialize the input parameters.
 *
 * \param mesh pointer to the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_INIT_PARAMETERS(mesh)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void  MMGS_Init_parameters(MMG5_pMesh mesh);

/* init file names */
/**
 * \brief Set the name of the input mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_INPUTMESHNAME(mesh,meshin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin);

/**
 * \brief Set the name of the output mesh file.
 *
 * \param mesh pointer to the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_OUTPUTMESHNAME(mesh,meshout,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshout\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout);

/**
 * \brief Set the name of the input solution file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_INPUTSOLNAME(mesh,sol,solin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin);

/**
 * \brief Set the name of the output solution file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param solout name of the output solution file.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_OUTPUTSOLNAME(mesh,sol,solout,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solout\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout);

/**
 * \brief Set the name of the input parameter file.
 *
 * \param mesh pointer to the mesh structure.
 * \param fparamin name of the input parameter file.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_INPUTPARAMNAME(mesh,fparamin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: fparamin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_inputParamName(MMG5_pMesh mesh, const char* fparamin);

/* init structure sizes */
/**
 * \brief Initialize an array of solution fields: set dimension, types and
 * number of fields.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typEntity types of solution entities (vertices, triangles, ...
 *        see \ref MMG5_entities for possible values).
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial, ..., see \ref MMG5_type for possible values)
 * \return 0 on failure, 1 otherwise.
 *
 * To use to initialize an array of solution fields (not used by Mmg itself).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: typEntity,typSol\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: np\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int  MMGS_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, MMG5_int np, int typSol);

/**
 * \brief Initialize an array of solution fields defined at vertices: set
 * dimension, types and number of fields.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an allocatable sol structure.
 * \param nsols number of solutions per entity
 * \param nentities number of entities
 * \param typSol    Array of size nsol listing the type of the solutions
 *                  (scalar, vectorial, ..., see \ref MMG5_type for possible values)
 * \return 0 on failure, 1 otherwise.
 *
 * To use to initialize an array of solution fields (not used by Mmg itself).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_SOLSATVERTICESSIZE(mesh,sol,nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: nsols\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: nentities\n
 * >     INTEGER, INTENT(IN)           :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Set_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol *sol,int nsols,
                                               MMG5_int nentities, int *typSol);

/**
 * \brief Set the number of vertices, triangles and edges of the
 * mesh and allocate the associated tables.
 *
 * \param mesh pointer to the mesh structure.
 * \param np number of vertices.
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 on failure, 1 otherwise.
 *
 * If called again, this function resets the
 * whole mesh to reallocate it at the new size
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_MESHSIZE(mesh,np,nt,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT)            :: np,nt,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_meshSize(MMG5_pMesh mesh, MMG5_int np, MMG5_int nt, MMG5_int na);

/* init structure data */
/**
 * \brief Set the coordinates of a single vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 coordinate of the vertex along the first dimension.
 * \param c1 coordinate of the vertex along the second dimension.
 * \param c2 coordinate of the vertex along the third dimension.
 * \param ref vertex reference.
 * \param pos position of the vertex in the mesh.
 * \return 1.
 *
 * \brief Set vertex coordinates \a c0, \a c1,\a c2 and reference \a ref at
 * position \a pos in the mesh structure (\a pos from 1 to the number of
 * vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_VERTEX(mesh,c0,c1,c2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: c0,c1,c2\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                     double c2, MMG5_int ref, MMG5_int pos);

/**
 * \brief Set the coordinates and references of all vertices in a mesh
 *
 * \param mesh pointer to the mesh structure.
 * \param vertices array of vertex coordinates in the order
 * \f$[x_1, y_1, z_1, x_2, \ldots, z_N]\f$ where \f$N\f$ is the number of vertices.
 * \param refs array of references.
 * The reference of vertex \f$i\f$ is stored in refs[\f$i-1\f$].
 * \return 1.
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs array)
 * > ! SUBROUTINE MMGS_SET_VERTICES(mesh,vertices,refs,retval)\n
 * > !   MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * > !   REAL(KIND=8), DIMENSION(*),INTENT(IN)       :: vertices\n
 * > !   INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN) :: refs\n
 * > !   INTEGER, INTENT(OUT)                        :: retval\n
 * > ! END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int  MMGS_Set_vertices(MMG5_pMesh mesh, double *vertices,MMG5_int *refs);

/**
 * \brief Set the coordinates and reference of a single triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets a triangle with vertices \a v0, \a v1, \a v2 and reference
 * \a ref at position \a pos in the mesh structure (\a pos from 1 to the number
 * of triangles included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_TRIANGLE(mesh,v0,v1,v2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,v2,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_triangle(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                      MMG5_int v2, MMG5_int ref,MMG5_int pos);

/**
 * \brief Set the vertices and references of all triangles in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param tria pointer to an array of vertex numbers
 * Vertices of the \f$i^{th}\f$ triangles are stored in tria[(i-1)*3]\@3.
 * \param refs pointer to an array of triangle references.
 * refs[i-1] is the reference of the \f$i^{th}\f$ triangle.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs array)
 * >  ! SUBROUTINE MMGS_SET_TRIANGLES(mesh,tria,refs,retval)\n
 * >  !   MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * >  !   INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN) :: tria,refs\n
 * >  !   INTEGER, INTENT(OUT)                        :: retval\n
 * >  ! END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Set_triangles(MMG5_pMesh mesh, MMG5_int *tria, MMG5_int *refs);

/**
 * \brief Set the vertices and reference of one edge in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first extremity of the edge.
 * \param v1 second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * Assigns vertices \a v0, \a v1 and reference \a ref to the edge at position \a
 * pos in the mesh structure (\a pos from 1 to the number of edges included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_EDGE(mesh,v0,v1,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_edge(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int ref,MMG5_int pos);

/**
 * \brief Assign the "corner" attribute to a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex number.
 * \return 1.
 *
 * This function sets the corner attribute at vertex \a pos (\a pos from 1 to
 * the number of vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_CORNER(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_corner(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "corner" attribute from a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex number.
 * \return 1.
 *
 * This function removes the corner attribute from vertex \a pos (from 1 to the
 * number of vertices included).
 *
 * \remark Fortran interface
 *
 * >   SUBROUTINE MMGS_UNSET_CORNER(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Unset_corner(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex number.
 * \return 1.
 *
 * This function sets the required attribute at vertex \a k.
 * Vertices with this attribute will not be modified by the remesher.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_REQUIREDVERTEX(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int  MMGS_Set_requiredVertex(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex number.
 * \return 1.
 *
 * This function removes the required attribute from vertex \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_UNSET_REQUIREDVERTEX(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Unset_requiredVertex(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to a triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * This function sets the required attribute at triangle \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_REQUIREDTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_requiredTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * This function removes the required attribute from triangle \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_UNSET_REQUIREDTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int  MMGS_Unset_requiredTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "ridge" attribute to an edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * This function gives the ridge attribute to edge \a k. This influences how
 * this edge is treated by the remesher.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_RIDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_ridge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "ridge" attribute from an edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * This function removes the ridge attribute from edge \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_UNSET_RIDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Unset_ridge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to an edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * This function makes edge \a k a required edge. Required edges will not be
 * modified by the remesher.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_REQUIREDEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_requiredEdge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from an edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * This function removes the "required" attribute from edge \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_UNSET_REQUIREDEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Unset_requiredEdge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Set the vertices and references of all edges in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param edges pointer to an array of edges.
 *              The vertices of edge i should be given in
 *              edges[(i-1)*2] and edges[(i-1)*2+1].
 * \param refs pointer to an array of references.
 *             refs[i-1] is the reference of edge i.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_EDGES(mesh,edges,refs,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: edges(*),refs(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int MMGS_Set_edges(MMG5_pMesh mesh, MMG5_int *edges, MMG5_int* refs);

/**
 * \brief Get vertices, references and attributes of all edges in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param edges pointer to an array of edges.
 *   The vertices of edge i are stored in edges[(i-1)*2] and edges[(i-1)*2+1].
 * \param refs edge references. refs[i-1] is the reference of edge \a i.
 * \param areRidges 1 if the edge is a ridge, 0 otherwise.
 * \param areRequired 1 if the edge is required, 0 otherwise.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_EDGES(mesh,edges,refs,areRidges,areRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: refs(*),edges(*)\n
 * >     INTEGER, INTENT(OUT)          :: areRequired(*),areRidges(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int MMGS_Get_edges(MMG5_pMesh mesh,MMG5_int *edges,MMG5_int* refs,
                                   int *areRidges,int *areRequired);

/**
 * \brief Set the normal orientation at a single vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index
 * \param n0 x componant of the normal at vertex \a k.
 * \param n1 y componant of the normal at vertex \a k.
 * \param n2 z componant of the normal at vertex \a k.
 *
 * \return 1 on success.
 *
 * Set normal (n0,n1,n2) at vertex \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_NORMALATVERTEX(mesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8), INTENT(IN)      :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_normalAtVertex(MMG5_pMesh mesh, MMG5_int k, double n0, double n1, double n2) ;

/**
 * \brief Get the quality measure of a triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of the triangle for which we want to get the quality.
 * \return the computed quality or 0 on failure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_TRIANGLEQUALITY(mesh,met,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,met\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8), INTENT(OUT)     :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  double MMGS_Get_triangleQuality(MMG5_pMesh mesh, MMG5_pSol met, MMG5_int k);

/**
 * \brief Set a single element of a scalar solution structure.
 *
 * \param met pointer to the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets the scalar value \a s at position \a pos in the solution
 * structure (\a pos from 1 to the number of vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_SCALARSOL(met,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: s\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_scalarSol(MMG5_pSol met, double s, MMG5_int pos);

/**
 * \brief Set the values of all elements of a scalar solution structure.
 *
 * \param met pointer to the sol structure.
 * \param s array of scalar solutions values.
 *  s[i-1] is the solution at vertex i.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: s\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_scalarSols(MMG5_pSol met, double *s);

/**
 * \brief Set a single element of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets the vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos
 * in the solution structure (\a pos from 1 to the number of vertices
 * included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_VECTORSOL(met,vx,vy,vz,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: vx,vy,vz\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz, MMG5_int pos);

/**
 * \brief Set all elements of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of vectorial solutions
 * sols[3*(i-1)]\@3 is the solution at vertex i
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets a vector-valued solution at each element of solution
 * structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: sols\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Set_vectorSols(MMG5_pSol met, double *sols);

/**
 * \brief Set a single element of a tensor solution structure.
 *
 * \param met pointer to the sol structure.
 * \param m11 value of the tensorial solution at position (1,1) in the tensor.
 * \param m12 value of the tensorial solution at position (1,2) in the tensor.
 * \param m13 value of the tensorial solution at position (1,3) in the tensor.
 * \param m22 value of the tensorial solution at position (2,2) in the tensor.
 * \param m23 value of the tensorial solution at position (2,3) in the tensor.
 * \param m33 value of the tensorial solution at position (3,3) in the tensor.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets a tensor value at position \a pos in solution
 * structure (\a pos from 1 to the number of vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_TENSORSOL(met,m11,m12,m13,m22,m23,m33,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                                      double m22,double m23, double m33, MMG5_int pos);

/**
 * \brief Set all elements of a tensor solution structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of tensorial solutions.
 * sols[6*(i-1)]\@6 is the solution at vertex i
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: sols\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Set_tensorSols(MMG5_pSol met, double *sols);

/**
 * \brief Set a single element of one out of multiple solution fields that are defined on vertices.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s solution(s) at mesh vertex \a pos.
 * \param pos index of the vertex on which we set the solution.
 *
 * \return 0 on failure, 1 otherwise.
 *
 * Set values of the solution at field \a i of the solution array and at
 * position \pos (\a pos from 1 to the number of vertices included and \a i from 1
 * to the number of solutions). The type of solution is determined from \a sol.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_ITHSOL_INSOLSATVERTICES(sol,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)         :: pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Set_ithSol_inSolsAtVertices(MMG5_pSol sol, int i, double* s, MMG5_int pos);

/**
 * \brief Set all elements of one out of multiple solution fields that are defined on vertices.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s array of solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 on failure, 1 otherwise.
 *
 * Set values of the solution at field \a i of the solution array (\a i from
 * 1 to the number of solutions). The type of solution is determined from \a sol.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_ITHSOLS_INSOLSATVERTICES(sol,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Set_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);

/* check init */
/**
 * \brief Check if the numbers of given entities match with mesh and solution
 * size and check mesh data.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the solution structure.
 * \return 0 on failure, 1 otherwise.
 *
 * This function checks if the numbers of given entities match with the mesh and
 * solution sizes and checks the mesh data. Use of this function is not
 * mandatory.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_CHK_MESHDATA(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/** functions to set parameters */
/**
 * \brief set an integer parameter of the remesher
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure (unused).
 * \param iparam integer parameter to set (see \ref MMGS_Param for a
 *   list of parameters that can be set).
 * \param val value for the parameter.
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets integer parameter \a iparam to value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_IPARAMETER(mesh,sol,iparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     MMG5_DATA_PTR_T                :: sol\n
 * >     INTEGER, INTENT(IN)            :: iparam\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: val\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam, MMG5_int val);

/**
 * \brief set a real-valued parameter of the remesher
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure (unused).
 * \param dparam double parameter to set (see \ref MMGS_Param for a
 *    list of parameters that can be set).
 * \param val value of the parameter.
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets the double parameter \a dparam to value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_DPARAMETER(mesh,sol,dparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     MMG5_DATA_PTR_T               :: sol\n
 * >     INTEGER, INTENT(IN)           :: dparam\n
 * >     REAL(KIND=8), INTENT(IN)      :: val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam, double val);

/**
 * \brief set a local parameter
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param hmin minimal edge length.
 * \param hmax maximal edge length.
 * \param Hausdorff distance.
 * \return 0 on failure, 1 otherwise.
 *
 * Set local parameters: set the hausdorff distance at \a hausd, the minmal edge
 * length at \a hmin and the maximal edge length at \a hmax for all
 * elements of type \a typ and reference \a ref.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_LOCALPARAMETER(mesh,sol,typ,ref,hmin,hmax,hausd,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh,sol\n
 * >     INTEGER, INTENT(IN)            :: typ\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: ref\n
 * >     REAL(KIND=8), INTENT(IN)       :: hmin,hmax,hausd\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ, MMG5_int ref,
                                            double hmin, double hmax, double hausd);

/**
 * \brief Set the reference mapping for the elements of reference
 *  \a ref in level-set discretization mode.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param ref input triangle reference.
 * \param split MMG5_MMAT_NoSplit if the entity must not be split, MMG5_MMAT_Split otherwise
 * \param rmin reference for the negative side after LS discretization
 * \param rplus reference for the positive side after LS discretization
 * \return 0 on failure, 1 otherwise.
 *
 * With this function you can determine which references will be given to the
 * triangles on both sides of the level set, after discretization. Negative and
 * positive here refer to areas where the function is smaller or larger,
 * respectively, than the isovalue of the level set.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_MULTIMAT(mesh,sol,ref,split,rmin,rplus,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref,rmin,rplus\n
 * >     INTEGER, INTENT(IN)           :: split\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Set_multiMat(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int ref,
                                        int split,MMG5_int rmin, MMG5_int rplus);

/**
 * \brief Set a new level-set base reference.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param br new level-set base reference.
 * \return 0 on failure, 1 otherwise.
 *
 * Set a new level-set base reference of ref \a br in level-set discretization
 * mode. Base references are boundary conditions to which an implicit domain can
 * be attached. All implicit volumes that are not attached to listed base
 * references are deleted as spurious volumes by the \a rmc option.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_LSBASEREFERENCE(mesh,sol,br,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: br\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Set_lsBaseReference(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int br);

/** recover data */
/**
 * \brief Get the number of vertices, triangles, and edges of the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param np pointer to the number of vertices.
 * \param nt pointer to the number of triangles.
 * \param na pointer to the number of edges.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_MESHSIZE(mesh,np,nt,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT)            :: np,nt,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_meshSize(MMG5_pMesh mesh, MMG5_int* np, MMG5_int* nt, MMG5_int* na);

/**
 * \brief Get the number of elements, dimension, and type of a solution.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typEntity pointer to the type of entities to which solutions are applied.
 *        (see \ref MMG5_entities for possible values)
 * \param np pointer to the number of elements in the solution.
 * \param typSol pointer to the type of the solution (\ref MMG5_Scalar, \ref MMG5_Vector,
 *    \ref MMG5_Tensor, \ref MMG5_Notype)
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER                       :: typEntity,typSol\n
 * >     INTEGER(MMG5F_INT)            :: np\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, MMG5_int* np,
                                     int* typSol);

/**
 * \brief Get the number of elements, type, and dimensions of several solutions defined on vertices.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an array of sol structures.
 * \param nsols number of solutions per entity
 * \param nentities pointer to the number of entities.
 * \param typSol array of size MMG5_NSOL_MAX to store type of each solution
 * (scalar, vectorial, ..., see \ref MMG5_type for possible values).
 *
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_SOLSATVERTICESSIZE(mesh,sol,nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER                       :: nsols\n
 * >     INTEGER(MMG5F_INT)            :: nentities\n
 * >     INTEGER                       :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Get_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol* sol,int *nsols,
                                                  MMG5_int* nentities,int* typSol);

/**
 * \brief Get the coordinates \a c0, \a c1,\a c2 and reference \a ref of the
 * next vertex of \a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 pointer to the coordinate of the vertex along the first dimension.
 * \param c1 pointer to the coordinate of the vertex along the second dimension.
 * \param c2 pointer to the coordinate of the vertex along the third dimension.
 * \param ref pointer to the vertex reference.
 * \param isCorner pointer to the flag saying if the vertex is corner.
 * \param isRequired pointer to the flag saying if the vertex is required.
 * \return 1.
 *
 * This function retrieves the coordinates \a c0, \a c1,\a c2, reference \a ref,
 * and attributes of the next vertex of a mesh. It is meant to be used in a loop
 * over all vertices. When this function has been called as many times as there
 * are vertices, the internal loop counter will be reset. To obtain data for a
 * specific vertex, the \ref MMGS_GetByIdx_vertex function can be used instead.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_VERTEX(mesh,c0,c1,c2,ref,isCorner,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1,c2\n
 * >     INTEGER(MMG5F_INT)            :: ref\n
 * >     INTEGER                       :: isCorner,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
                                    int* isCorner, int* isRequired);

/**
 * \brief Get the coordinates and reference of a specific vertex in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 pointer to the coordinate of the vertex along the first dimension.
 * \param c1 pointer to the coordinate of the vertex along the second dimension.
 * \param c2 pointer to the coordinate of the vertex along the third dimension.
 * \param ref pointer to the vertex reference.
 * \param isCorner pointer to the flag saying if the vertex is corner.
 * \param isRequired pointer to the flag saying if the vertex is required.
 * \param idx index of vertex to get.
 * \return 1.
 *
 * This function retrieves the coordinates \a c0, \a c1, \a c2 and reference \a ref of
 * vertex \a idx of mesh, as well as its "corner" and "required" attributes.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GETBYIDX_VERTEX(mesh,c0,c1,c2,ref,isCorner,isRequired,idx,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1,c2\n
 * >     INTEGER                       :: isCorner,isRequired\n
 * >     INTEGER(MMG5F_INT)            :: ref,idx\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int  MMGS_GetByIdx_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
                                          int* isCorner, int* isRequired,MMG5_int idx);

/**
 * \brief Get the coordinates, references and attributes of all vertices in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param vertices pointer to the array of coordinates.
 * The coordinates of vertex \a i are stored in
 * vertices[(i-1)*3]\@3.
 * \param refs pointer to the array of vertex references.
 * The ref of vertex \a i is stored in refs[i-1].
 * \param areCorners pointer to the array of flags saying if
 *   vertices are corners.
 * areCorners[i-1]=1 if vertex \a i is corner.
 * \param areRequired pointer to the table of flags saying if vertices
 * are required. areRequired[i-1]=1 if vertex \a i is required.
 * \return 1.
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs,areCorners and areRequired arrays)
 * > !  SUBROUTINE MMGS_GET_VERTICES(mesh,vertices,refs,areCorners,&\n
 * > !                                areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)          :: mesh\n
 * > !    REAL(KIND=8),DIMENSION(*), INTENT(OUT) :: vertices\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)       :: refs\n
 * > !    INTEGER, DIMENSION(*)                  :: areCorners,areRequired\n
 * > !    INTEGER, INTENT(OUT)                   :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_vertices(MMG5_pMesh mesh, double* vertices, MMG5_int* refs,
                                      int* areCorners, int* areRequired);

/**
 * \brief Get the vertices, reference, and required attribute of the next
 * triangle in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the first vertex of the triangle.
 * \param v1 pointer to the second vertex of the triangle.
 * \param v2 pointer to the third vertex of the triangle.
 * \param ref pointer to the triangle reference.
 * \param isRequired pointer to the flag saying if the triangle is required.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the vertices \a v0, \a v1, \a v2, reference \a ref,
 * and required attribute \a isRequired of the next triangle of \a mesh. It is
 * meant to be called in a loop over all triangles. When it has been called as
 * many times as there are triangles, the internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_TRIANGLE(mesh,v0,v1,v2,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: v0,v1,v2,ref\n
 * >     INTEGER                        :: isRequired\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_triangle(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,
                                       MMG5_int* ref, int* isRequired);

/**
 * \brief Get the vertices, references, and required attributes of all triangles
 * in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param tria pointer to an array of vertices
 * Vertices of triangle \a i are stored in tria[(i-1)*3]\@3.
 * \param refs pointer to the array of triangles references.
 * refs[i-1] is the ref of triangle \a i.
 * \param areRequired pointer to an array of flags saying if triangles
 * are required. areRequired[i-1]=1 if triangle \a i
 * is required.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs and areRequired arrays)
 * > !  SUBROUTINE MMGS_GET_TRIANGLES(mesh,tria,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*),INTENT(OUT) :: tria\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)             :: refs\n
 * > !    INTEGER, DIMENSION(*)                        :: areRequired\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_triangles(MMG5_pMesh mesh, MMG5_int* tria, MMG5_int* refs,
                                       int* areRequired);
/**
 * \brief Get the vertices, reference, and attributes of the next edge in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param e0 pointer to the index of the first vertex of the edge.
 * \param e1 pointer to the index of the second  vertex of the edge.
 * \param ref pointer to the edge reference.
 * \param isRidge pointer to the flag saying if the edge is ridge.
 * \param isRequired pointer to the flag saying if the edge is required.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the extremities \a e0, \a e1, reference \a ref, and
 * attributes of the next edge of \a mesh. It is meant to be called in a loop
 * over all edges. When it has been called as many times as there are edges in
 * the mesh, the internal edge counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_EDGE(mesh,e0,e1,ref,isRidge,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: e0,e1\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER                        :: isRidge,isRequired\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_edge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref,
                                  int* isRidge, int* isRequired);

/**
 * \brief Get the normal orientation at an edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex number
 * \param n0 x componant of the normal at vertex \a k.
 * \param n1 y componant of the normal at vertex \a k.
 * \param n2 z componant of the normal at vertex \a k.
 *
 * \return 1 on success.
 *
 * This function retrieves the normal (n0,n1,n2) at vertex \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_NORMALATVERTEX(mesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8)                  :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_normalAtVertex(MMG5_pMesh mesh, MMG5_int k, double *n0, double *n1, double *n2) ;

/**
 * \brief Get the next element of a scalar solution structure.
 *
 * \param met pointer to the sol structure.
 * \param s pointer to the scalar solution value.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the next element \a s of the solution field \a
 * met. It is meant to be called in a loop over all elements. When it has been
 * called as many times as there are elements in the solution, the internal loop
 * counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_SCALARSOL(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_scalarSol(MMG5_pSol met, double* s);

/**
 * \brief Get all elements of a scalar solution structure.
 *
 * \param met pointer to the solution structure.
 * \param s array of scalar solutions at mesh vertices. s[i-1] is
 * the solution at vertex i.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_scalarSols(MMG5_pSol met, double* s);

/**
 * \brief Get the next element of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the next vector-valued element \f$(v_x,v_y,vz)\f$ of
 * a solution field. It is meant to be called in a loop over all elements.  When
 * it has been called as many times as there are elements in the solution, the
 * internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_VECTORSOL(met,vx,vy,vz,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: vx,vy,vz\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);

/**
 * \brief Get all elements of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of solutions at mesh vertices. sols[3*(i-1)]\@3 is
 * the solution at vertex i.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves all elements of a vector-valued solution field.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_vectorSols(MMG5_pSol met, double* sols);

/**
 * \brief Get the next element of a tensor solution structure.
 *
 * \param met pointer to the sol structure.
 * \param m11 pointer to the position (1,1) in the solution tensor.
 * \param m12 pointer to the position (1,2) in the solution tensor.
 * \param m13 pointer to the position (1,3) in the solution tensor.
 * \param m22 pointer to the position (2,2) in the solution tensor.
 * \param m23 pointer to the position (2,3) in the solution tensor.
 * \param m33 pointer to the position (3,3) in the solution tensor.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the next element
 * \f$(m_{11},m_{12},m_{13},m_{22},m_{23},m_{33})\f$ of a tensor-valued solution
 * field.  It is meant to be called in a loop over all vertices. When it has
 * been called as many times as there are elements in the solution, the internal
 * loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_TENSORSOL(met,m11,m12,m13,m22,m23,m33,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_tensorSol(MMG5_pSol met, double *m11, double *m12, double *m13,
                                                     double *m22, double *m23, double *m33);

/**
 * \brief Get all elements of a tensor solution field.
 *
 * \param met pointer to the sol structure.
 * \param sols array of solution values.
 *   sols[6*(i-1)]\@6 is the solution at vertex i.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)           :: met\n
 * >     REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)                    :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_tensorSols(MMG5_pSol met, double *sols);

/**
 * \brief Get one out of several solutions at a specific vertex.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s solution(s) at mesh vertex \a pos. The required size
 *   of this array depends on the type of solution.
 * \param pos index of the vertex on which we get the solution.
 *
 * \return 0 on failure, 1 otherwise.
 *
 * This function retreives the value of field \a i of the solution array at
 * vertex \a pos.  (\a pos from 1 to the number of vertices included and \a i
 * from 1 to the number of solutions). It works for any type of solution; the
 * types are inferred from \a sol.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_ITHSOL_INSOLSATVERTICES(sol,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)         :: pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int  MMGS_Get_ithSol_inSolsAtVertices(MMG5_pSol sol, int i, double* s, MMG5_int pos);

/**
 * \brief Get one out of several solutions at all vertices in the mesh.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s array of solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the values of field \a i of the solution array \a sol
 * (\a i from 1 to \a the number of solutions). It works for any type of solution;
 * the type is inferred from \a sol.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_ITHSOLS_INSOLSATVERTICES(sol,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Get_ithSols_inSolsAtVertices(MMG5_pSol sol, int i, double* s);

/**
 * \brief Get the value of an integer parameter of the remesher.
 *
 * \param mesh pointer to the mesh structure.
 * \param iparam integer parameter to get (see \ref MMGS_Param structure).
 * \return The value of the parameter.
 *
 * This function retrieves the value of integer parameter \a iparam (see \ref
 * MMGS_Param for a list of parameters). It returns the value of the parameter,
 * or zero if the value of \a iparam is not recognized.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_IPARAMETER(mesh,iparam,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: iparam\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_iparameter(MMG5_pMesh mesh, MMG5_int iparam);

/* input/output functions */
/**
 * \brief Load a mesh (in .mesh/.mesb format) from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * This function reads .mesh (ASCII) and .meshb (binary) files. If the name
 * contains ".mesh" the file will be read as an ASCII file and if the name
 * contains .meshb it be read as a binary file. If the file contains neither of
 * these strings the function will first try to open "[filename].meshb"
 * and if this fails it will try "[filename].mesh".
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_loadMesh(MMG5_pMesh mesh, const char* filename);

/**
 * \brief Load a mesh and optionally a solution in VTP (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * This function reads a mesh and optionally one data field in VTK vtp file
 * format (.vtp extension). We read only low-order vertices, edges, triangles
 * and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTPMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtpMesh(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol sol,
                                    const char *filename);

/**
 * \brief Load a mesh and multiple solutions in VTP (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Read a mesh and multiple data fields in VTK vtp file format (.vtp extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTPMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtpMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol,
                                                const char *filename);

/**
 * \brief Load a mesh and possibly data in VTU (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Read a mesh and optionally one data field in VTK vtu file format (.vtu
 * extension). We read only low-order vertices, edges, triangles and
 * quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTUMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtuMesh(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol sol, const char *filename);

/**
 * \brief Load a mesh and multiple solutions in VTU (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Read a mesh and multiple data field in VTK vtu file format (.vtu extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTUMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtuMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol,
                                                const char *filename);

/**
 * \brief Load a mesh and possibly data in VTK format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Read a mesh and optionally one data field in VTK vtk file format (.vtk extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTKMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtkMesh(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol sol,
                                    const char *filename);

/**
 * \brief Load a mesh and multiple solutions in VTK format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Read a mesh and multiple data field in VTK vtk file format (.vtk extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTKMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtkMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol,
                                                const char *filename);

/**
 * \brief Load a mesh and possibly a solution in .msh format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Read a mesh and optionally one data field in MSH file format (.msh extension). We read
 * only low-order vertices, edges, triangles, quadrangles, tetrahedra and prisms.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADMSHMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadMshMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename);

/**
 * \brief Load a mesh and all data from a file in MSH format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to a list of solution structures.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Read a mesh and multiple data in MSH file format (.msh extension). We read only
 * low-order vertices, edges, triangles, quadrangles, tetrahedra and prisms.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADMSHMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadMshMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol, const char *filename);

/**
 * \brief Load a mesh and all data from a file. The format will be guessed from the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADGENERICMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadGenericMesh(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol sol,
                                        const char *filename);

/**
 * \brief Save a mesh in .mesh or .meshb format.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * This function saves a mesh in .mesh or .meshb format (depending on the
 * filename extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_saveMesh(MMG5_pMesh mesh, const char *filename);

/**
 * \brief Write mesh and optionally one data field in MSH  file format (.msh extension).
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * The file is saved in ASCII format for .msh extension, an in binary format for
 * a .mshb extension.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEMSHMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_saveMshMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename);

/**
 * \brief Save a mesh and multiple data fields in MSH format, ascii or binary
 * depending on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function saves a mesh and multiple data fields (that are considered as
 * solutions and not metrics, thus, we do nothing over the ridge vertices) in MSH
 * file format (.msh extension). The file is saved in ASCII format for .msh
 * extension and in binary format for a .mshb extension.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEMSHMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_saveMshMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol, const char *filename);

/**
 * \brief Write mesh and optionally one data field in Vtk file format (.vtk extension).
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEVTKMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_saveVtkMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename);

/**
 * \brief Save a mesh and multiple data fields in VTK format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and a list of data fields in Vtk file format (.vtk extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEVTKMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_saveVtkMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol, const char *filename);

/**
 * \brief Write mesh and optionally one data field vtu Vtk file format (.vtu extension).
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEVTUMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_saveVtuMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename);

/**
 * \brief Write a mesh and multiple data fields in vtu Vtk file format (.vtu extension).
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEVTUMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int MMGS_saveVtuMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol,
                                                  const char *filename);

/**
 * \brief Save a mesh and optionally one data field in VTP format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and optionally one data in polydata Vtk file
 * format (.vtp extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEVTPMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int MMGS_saveVtpMesh(MMG5_pMesh mesh, MMG5_pSol sol,
                                      const char *filename);

/**
 * \brief Save a mesh and multiple data fields in VTP format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and multiple data fields in polydata Vtk file
 * format (.vtp extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEVTPMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int MMGS_saveVtpMesh_and_allData(MMG5_pMesh mesh, MMG5_pSol *sol,
                                                  const char *filename);

/**
 * \brief Save mesh data in a file whose format depends on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEGENERICMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int MMGS_saveGenericMesh(MMG5_pMesh mesh, MMG5_pSol sol,
                                          const char *filename);

/**
 * \brief Load a metric field (or other solution) in medit's .sol format.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * Load metric field. The solution file (in medit file format) must contain
 * only 1 solution: the metric.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADSOL(mesh,met,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_loadSol(MMG5_pMesh mesh, MMG5_pSol met, const char* filename);

/**
 * \brief Load one or more solutions in a solution file in medit file format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solutions array
 * \param filename name of the file to load.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADALLSOLS(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_loadAllSols(MMG5_pMesh mesh, MMG5_pSol *sol,
                                     const char* filename);

/**
 * \brief Write an isotropic or anisotropic metric in medit file format.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVESOL(mesh,met,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_saveSol(MMG5_pMesh mesh, MMG5_pSol met,
                                 const char *filename);

/**
 * \brief Save one or more solutions in a solution file in medit file format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solutions array
 * \param filename name of the solution file.
 * \return 0 or -1 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SAVEALLSOLS(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_saveAllSols(MMG5_pMesh mesh, MMG5_pSol *sol,
                                    const char *filename);

/**
 * \brief Deallocate an array of solution fields
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an array of solution structures (that stores solution fields).
 * \return 1
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_Free_allSols(mesh,sol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Free_allSols(MMG5_pMesh mesh,MMG5_pSol *sol);

/* deallocations */
/**
 * \brief Deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments.
 *
 * For the \ref MMGS_mmgslib function, you need
 * to call the \ref MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the \ref MMGS_mmgsls function, you need
 * to call the \ref MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here, \a your_mesh is a \ref MMG5_pMesh, \a your_metric and \a your_level_set
 * are \ref MMG5_pSol.
 *
 * \return 0 on failure, 1 on success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
LIBMMGS_EXPORT int MMGS_Free_all(const int starter,...);

/**
 * \brief Structure deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments.
 *
 * For the \ref MMGS_mmgslib function, you need
 * to call the \ref MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the \ref MMGS_mmgsls function, you need
 * to call the \ref MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here, \a your_mesh is a pointer to \ref MMG5_pMesh and \a your_metric and
 * \a your_level_set are pointers to \ref MMG5_pSol.
 *
 * \return 0 on failure, 1 on success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
LIBMMGS_EXPORT int MMGS_Free_structures(const int starter,...);

/**
 * \brief Structure deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments.
 *
 * For the \ref MMGS_mmgslib function, you need
 * to call the \ref MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the \ref MMGS_mmgsls function, you need
 * to call the \ref MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here, \a your_mesh is a \ref MMG5_pMesh, \a your_metric and \a your_level_set
 * are \ref MMG5_pSol.
 *
 * \return 0 on failure, 1 on success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
LIBMMGS_EXPORT int MMGS_Free_names(const int starter,...);

/* library */
/**
 * \brief Main "program" for mesh adaptation.
 *
 * \param mesh pointer to the mesh structure.  \param met pointer to the sol
 * (metric) structure.  \return \ref MMG5_SUCCESS on success, \ref
 * MMG5_LOWFAILURE in case there is a failure but a conform mesh can be returned
 * or \ref MMG5_STRONGFAILURE if there is a failure and we can't save the mesh.
 *
 * Main program for the library.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_MMGSLIB(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_mmgslib(MMG5_pMesh mesh, MMG5_pSol met);

/**
 * \brief Main "program" for level-set discretization.
 *
 * \param mesh pointer to the mesh structure.  \param sol pointer to the sol
 * (level-set) structure.  \param met pointer to the sol (metric) structure
 * (optionnal).  \return \ref MMG5_SUCCESS on success, \ref MMG5_LOWFAILURE if
 * there is a a failure but a conform mesh is saved or \ref MMG5_STRONGFAILURE
 * if there is a a failure and we can't save the mesh.
 *
 * Main program for level set discretization library. If a metric \a met is
 * provided, use it to adapt the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_MMGSLS(mesh,sol,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     MMG5_DATA_PTR_T                :: met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_mmgsls(MMG5_pMesh mesh,  MMG5_pSol sol,MMG5_pSol met);

/** To associate function pointers without calling MMGS_mmgslib */
/**
 * \brief Set function pointers for caltet, lenedg, defsiz and gradsiz.
 *
 * \param mesh pointer to the mesh structure (unused).
 * \param met pointer to the sol structure (unused).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SETFUNC(mesh,met)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void  MMGS_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Get the number of non-boundary edges.
 *
 * \param mesh pointer to the mesh structure.
 * \param the number of edges pointer to the number of non boundary edges.
 * \return 0 on failure, 1 otherwise.
 *
 * Get the number of non-boundary edges (for DG methods for example). An edge is
 * boundary if it is located at the interface of 2 domains with different
 * references, if it belongs to one triangle only or if it is a singular edge
 * (ridge or required).
 * Append these edges to the list of edges.
 *
 * \warning reallocate the edge array and append the internal edges. This may
 * modify the behaviour of other functions.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_NUMBEROFNONBDYEDGES(mesh,nb_edges,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: nb_edges\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int MMGS_Get_numberOfNonBdyEdges(MMG5_pMesh mesh, MMG5_int* nb_edges);

/**
 * \brief Get vertices and reference of a non-boundary edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param e0 pointer to the first extremity of the edge (a vertex number).
 * \param e1 pointer to the second  extremity of the edge (a vertex number).
 * \param ref pointer to the edge reference.
 * \param idx index of the non boundary edge to get (between 1 and the number of edges)
 * \return 0 on failure, 1 otherwise.
 *
 * This function returns the vertices \a e0, \a e1 and reference \a ref of the non boundary
 * edge \a idx. An edge is boundary if it is located at
 * the interface of 2 domains with different references, if it belongs to one
 * triangle only or if it is a singular edge (ridge or required).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_NONBDYEDGE(mesh,e0,e1,ref,idx,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: e0,e1\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: idx\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int MMGS_Get_nonBdyEdge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1,
                                         MMG5_int* ref, MMG5_int idx);


/* Tools for the library */
/**
 * \brief Compute an isotropic size map according to the mean of the length of the
 * edges passing through a vertex.
 *
 * \param mesh pointer to the mesh structure
 * \param met pointer to the sol structure
 * \return 1 on success
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_DOSOL(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT extern int (*MMGS_doSol)(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Compute a constant size map.
 *
 * \param mesh pointer to the mesh structure
 * \param met pointer to the sol structure
 * \return 1 on success
 *
 * This function computes a constant size map according to mesh->info.hsiz,
 * mesh->info.hmin and mesh->info.hmax. It updates these 3 values if they are
 * not compatible.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_CONSTANTSIZE(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Print help for mmgs options.
 *
 * \param prog pointer to the program name.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_USAGE(prog,strlen0,retval)\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: prog\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_usage(char *prog);

/**
 * \brief Store command line arguments.
 *
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param sol pointer to a level-set or displacement
 *
 * \return 1.
 *
 * \remark no matching fortran function.
 *
 */
LIBMMGS_EXPORT int  MMGS_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol);

/**
 * \brief Print the default parameter values.
 *
 * \param mesh pointer to the mesh structure.
 * \return 0 on failure, 1 on success.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_DEFAULTVALUES(mesh,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_defaultValues(MMG5_pMesh mesh);

/**
 * \brief Store the info structure in the mesh structure.
 *
 * \param mesh pointer to the mesh structure.
 * \param info pointer to the info structure.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_STOCKOPTIONS(mesh,info,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,info\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/**
 * \brief Recover the info structure stored in the mesh structure.
 *
 * \param mesh pointer to the mesh structure.
 * \param info pointer to the info structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_DESTOCKOPTIONS(mesh,info)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,info\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void MMGS_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/**
 * \brief Return adjacent triangles of a triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param kel triangle index.
 * \param listri pointer to the array of indices of the three adjacent
 * triangles of triangle \a kel (the index is 0 if there is no adjacent triangle).
 * \return 1.
 *
 * Find the indices of the 3 adjacent elements of triangle \a
 * kel. \f$v_i = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_ADJATRI(mesh,kel,listri,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)                :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)                :: kel\n
 * >     INTEGER(MMG5F_INT), DIMENSION(3), INTENT(OUT) :: listri\n
 * >     INTEGER, INTENT(OUT)                          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_adjaTri(MMG5_pMesh mesh, MMG5_int kel, MMG5_int listri[3]);

/**
 * \brief Find adjacent vertices of a triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param ip vertex index.
 * \param start index of a triangle holding \a ip.
 * \param lispoi pointer to an array of size MMGS_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a ip.
 * \return nbpoi the number of adjacent vertices on success, 0 on failure.
 *
 * Find the indices of the adjacent vertices of the vertex \a
 * ip of the triangle \a start.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_ADJAVERTICESFAST(mesh,ip,start,lispoi,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)                        :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)                        :: ip,start\n
 * >     INTEGER(MMG5F_INT), DIMENSION(MMGS_LMAX), INTENT(OUT) :: lispoi\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT)                       :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_adjaVerticesFast(MMG5_pMesh mesh, MMG5_int ip,MMG5_int start, MMG5_int lispoi[MMGS_LMAX]);

/**
 * \brief Compute the real eigenvalues and eigenvectors of a symmetric matrix
 *
 * \param m upper part of a symmetric matrix diagonalizable in |R
 * \param lambda array of eigenvalues
 * \param vp array of eigenvectors
 *
 * \return the order of the eigenvalues
 *
 * Compute the real eigenvalues and eigenvectors of a symmetric matrix m whose
 * upper part is provided (m11, m12, m13, m22, m23, m33 in this order).
 *
 * lambda[0] is the eigenvalue associated to the eigenvector ( v[0][0], v[0,1], v[0,2] )
 * in C and to the eigenvector v(1,:) in fortran
 *
 * lambda[1] is the eigenvalue associated to the eigenvector ( v[1][0], v[1,1], v[1,2] )
 * in C and to the eigenvector v(2,:) in fortran
 *
 * lambda[2] is the eigenvalue associated to the eigenvector ( v[2][0], v[2,1], v[2,2] )
 * in C and to the eigenvector v(3,:) in fortran
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_COMPUTE_EIGENV(m,lambda,vp,retval)\n
 * >     REAL(KIND=8), INTENT(IN)         :: m(*)\n
 * >     REAL(KIND=8), INTENT(OUT)        :: lambda(*),vp(*)\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Compute_eigenv(double m[6],double lambda[3],double vp[3][3]);

/**
 * \brief Free a solution.
 *
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the solution structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_FREE_SOLUTIONS(mesh,sol)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void MMGS_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol);

/**
 * \brief Clean data (triangles and edges) linked to isosurface.
 *
 * \param mesh pointer to mesh structure
 *
 * \return 1 if successful, 0 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_CLEAN_ISOSURF(mesh,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)      :: mesh\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMGS_EXPORT int MMGS_Clean_isoSurf(MMG5_pMesh mesh);

/**
 * \brief Set common function pointers between mmgs and mmg3d to the matching mmgs
 * functions.
 */
LIBMMGS_EXPORT void MMGS_Set_commonFunc(void);

#ifdef __cplusplus
}
#endif

#endif
