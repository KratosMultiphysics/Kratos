/* ===========================================================================
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
** ===========================================================================
*/

/*
 * This file defines the C and Fortran headers of the mmg3d API, and
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
 *
 * \file mmg3d/libmmg3d.h
 * \ingroup API
 * \brief API headers and documentation for the mmg3d library, for volumetric meshes in 3D
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * These are the API functions for the mmg3d library. These functions allow to
 * load and save meshes and data defined on meshes; add, extract, or modify mesh
 * data; and to call the library functions that perform remeshing, level-set
 * discretization, and Lagrangian motion.
 *
 * Meshes are here defined in terms of vertices and three-dimensional objects:
 * tetrahedra and prisms. Optionally lower-dimensional entities can be present:
 * triangles, quadrilaterals and edges. All of these \a entities can have a
 * \a reference: an integer value that can serve as a group identifier. In
 * addition mesh entities can have \a attributes such as "ridge" or "required".
 *
 * Data defined on meshes can be for example functions that are meant for
 * level-set discretization, metric tensors that will govern edge lengths, and
 * vector fields governing lagrangian motion. These data can be scalar, vector,
 * or (symmetric) tensor-valued; and there can be more than one data item
 * associated with a mesh entity. These data are often referred to as \a
 * solutions.
 *
 * Three of the functions here are referred to as "programs", because they
 * perform the tasks for which Mmg is meant: remeshing, level-set discretization
 * and Lagrangian motion. The other functions merely serve to load and save data
 * and to perform pre- and post-processing. These programs actually behave much
 * like independent programs: they send diagnostic output to stdout and in rare
 * cases they may call the exit() function.
 *
 *
 * \htmlonly
 * <h2 class="groupheader">Examples</h2>
 * \endhtmlonly
 *
 * A very simple example code for mesh adaptation with automatic parsing of .mesh files
 * \dontinclude libexamples/mmg3d/adaptation_example0/example0_a/main.c
 * \skipline BEGIN_EXAMPLE
 * \until END_EXAMPLE
 *
 * Mesh adaptation example in which get/set functions are used to provide input to
 * the library and to extract the output mesh.
 * \include libexamples/mmg3d/adaptation_example0/example0_b/main.c
 *
 * Fortran example.
 * \include libexamples/mmg3d/adaptation_example0_fortran/example0_a/main.F90
 *
 * Another Fortran example.
 * \include libexamples/mmg3d/adaptation_example0_fortran/example0_b/main.F90
 *
 * Mesh adaptation example.
 * \include libexamples/mmg3d/adaptation_example1/main.c
 *
 * Another mesh adaptation example.
 * \include libexamples/mmg3d/adaptation_example2/main.c
 *
 * Isosurface discretization example (with metric)
 * \include libexamples/mmg3d/IsosurfDiscretization_lsAndMetric/main.c
 *
 * Lagrangian motion example.
 * \include libexamples/mmg3d/LagrangianMotion_example0/main.c
 */

#ifndef MMG3DLIB_H
#define MMG3DLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mmg/common/libmmgtypes.h"
#include "mmg/mmg3d/mmg3d_export.h"

/**
 * Maximum array size when storing adjacent vertices (or ball) of a vertex.
 */
#define MMG3D_LMAX      10240

/**
 * \enum MMG3D_Param
 * \brief Input parameters for the mmg library.
 *
 * These are the input parameters for the mmg3d library functions. Options prefixed by
 * \a MMG3D_IPARAM require integer values and options prefixed by
 * \a MMG3D_DPARAM require real values. They can be set with the
 * \ref MMG3D_Set_iparameter and \ref MMG3D_Set_dparameter functions,
 * respectively.
 *
 */
enum MMG3D_Param {
  MMG3D_IPARAM_verbose,                   /*!< [-1..10], Level of verbosity */
  MMG3D_IPARAM_mem,                       /*!< [n/-1], Max memory size in MB or keep the default value */
  MMG3D_IPARAM_debug,                     /*!< [1/0], Turn on/off debug mode */
  MMG3D_IPARAM_angle,                     /*!< [1/0], Turn on/off angle detection */
  MMG3D_IPARAM_iso,                       /*!< [1/0], Enable level-set discretization (volume and surfaces) */
  MMG3D_IPARAM_isosurf,                   /*!< [1/0], Enable level-set discretization on the surfaces only */
  MMG3D_IPARAM_nofem,                     /*!< [1/0], Do not attempt to make the mesh suitable for finite-element computations */
  MMG3D_IPARAM_opnbdy,                    /*!< [1/0], Preserve triangles at interface of 2 domains with the same reference */
  MMG3D_IPARAM_lag,                       /*!< [-1/0/1/2], Enable Lagrangian motion */
  MMG3D_IPARAM_optim,                     /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  MMG3D_IPARAM_optimLES,                  /*!< [1/0], Strong mesh optimization for LES computations */
  MMG3D_IPARAM_noinsert,                  /*!< [1/0], Avoid/allow vertex insertion */
  MMG3D_IPARAM_noswap,                    /*!< [1/0], Avoid/allow edge or face flipping */
  MMG3D_IPARAM_nomove,                    /*!< [1/0], Avoid/allow vertex relocation */
  MMG3D_IPARAM_nosurf,                    /*!< [1/0], Avoid/allow surface modifications */
  MMG3D_IPARAM_nreg,                      /*!< [0/1], Enable regularization of normals */
  MMG3D_IPARAM_xreg,                      /*!< [0/1], Enable boundary regularization by moving vertices */
  MMG3D_IPARAM_numberOfLocalParam,        /*!< [n], Number of local parameters (which will be set with \ref MMG3D_Set_localParameter) */
  MMG3D_IPARAM_numberOfLSBaseReferences,  /*!< [n], Number of base references for bubble removal (requires \ref MMG3D_DPARAM_rmc) */
  MMG3D_IPARAM_numberOfMat,               /*!< [n], Number of materials in level-set mode */
  MMG3D_IPARAM_numsubdomain,              /*!< [0/n], Save only the subdomain (reference) n (0==all subdomains) */
  MMG3D_IPARAM_renum,                     /*!< [1/0], Turn on/off renumbering with Scotch */
  MMG3D_IPARAM_anisosize,                 /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  MMG3D_IPARAM_octree,                    /*!< [n], Max number of vertices per PROctree cell (DELAUNAY) */
  MMG3D_IPARAM_nosizreq,                  /*!< [0/1], Allow/avoid overwriting of sizes at required vertices (advanced usage) */
  MMG3D_IPARAM_isoref,                    /*!< [0/n], Isosurface boundary material reference */
  MMG3D_DPARAM_angleDetection,            /*!< [val], Value for angle detection (degrees) */
  MMG3D_DPARAM_hmin,                      /*!< [val], Minimal edge length */
  MMG3D_DPARAM_hmax,                      /*!< [val], Maximal edge length */
  MMG3D_DPARAM_hsiz,                      /*!< [val], Constant edge length */
  MMG3D_DPARAM_hausd,                     /*!< [val], Global Hausdorff distance (on all boundaries in the mesh) */
  MMG3D_DPARAM_hgrad,                     /*!< [val], Gradation */
  MMG3D_DPARAM_hgradreq,                  /*!< [val], Gradation on required entites (advanced usage) */
  MMG3D_DPARAM_ls,                        /*!< [val], Function value where the level set is to be discretized */
  MMG3D_DPARAM_xreg,                      /*!< [val], Relaxation parameter for boundary regularization (0<val<1) */
  MMG3D_DPARAM_rmc,                       /*!< [-1/val], Remove small disconnected components in level-set mode */
  MMG3D_PARAM_size,                       /*!< [n], Number of parameters */
};

/*--------------------------- function headers ---------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \brief Initialize a mesh structure and optionally the associated solution and metric structures.
 *
 * \param starter dummy argument used to initialize the variadic argument
 * list
 * \param ... variadic arguments that depend on the library function that you
 * want to call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet
 * MMG5_ARG_ppMet, &your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * Here, \a your_mesh is a \ref MMG5_pMesh, \a your_metric \a your_level_set and
 * \a your_displacement are \ref MMG5_pSol.
 *
 * \return 1 on success, 0 on failure
 *
 * This function allocates and initializes MMG structures.  All structures of
 * types MMG5_pMesh and MMG5_pSol that will be given as arguments to Mmg
 * functions must be initialized with this function.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 * \warning detected bugs:
 *   - some vertices along open boundaries end up with a normal (while they should not)
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Init_mesh(const int starter,...);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 *
 * \brief Initialize file names to their default values.
 *
 * This function initializes all file names to their default values.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_INIT_FILENAMES(mesh,sol)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT void  MMG3D_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);

/**
 * \brief Initialize parameters to their default values
 *
 * \param mesh pointer to the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_INIT_PARAMETERS(mesh)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT void  MMG3D_Init_parameters(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer to the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * \brief Set the name of input mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_INPUTMESHNAME(mesh,meshin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_inputMeshName(MMG5_pMesh mesh,const char* meshin);

/**
 * \param mesh pointer to the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * \brief Set the name of output mesh file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_OUTPUTMESHNAME(mesh,meshout,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshout\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT  int  MMG3D_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * \brief Set the name of input solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_INPUTSOLNAME(mesh,sol,solin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol,
                                              const char* solin);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 * \brief Set the name of the output solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_OUTPUTSOLNAME(mesh,sol,solout,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solout\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol,
                                              const char* solout);

/**
 * \param mesh pointer to the mesh structure.
 * \param fparamin name of the input parameter file.
 * \return 1.
 *
 * \brief Set the name of the input parameter file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_INPUTPARAMNAME(mesh,fparamin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: fparamin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG3D_EXPORT int  MMG3D_Set_inputParamName(MMG5_pMesh mesh, const char* fparamin);


/* init structure sizes */
/**
 * \brief Initialize a solution field.
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typEntity type of entities on which the solution is defined (vertices, triangles, ...).
 *        See \ref MMG5_entities for the defined entity types. Currently only \ref MMG5_Vertex
 *        is supported.
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial, ...,
 *        see \ref MMG5_type for possible values).
 * \return 0 if failed, 1 otherwise.
 *
 * Initialize a solution field: set dimension, types and number of data.
 * To use to initialize a metric, a level-set or a displacement field.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: np\n
 * >     INTEGER, INTENT(IN)           :: typEntity,typSol\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity,
                         MMG5_int np, int typSol);

/**
 * \brief Initialize an array of solution values defined at vertices
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an allocatable sol structure.
 * \param nsols number of solutions per entity
 * \param nentities number of vertices
 * \param typSol Array of size nsols listing the type of the solutions
 *                 (scalar, vectorial, ..., see \ref MMG5_type for possible values).
 * \return 0 if failed, 1 otherwise.
 *
 * Initialize a solution field defined at vertices: set dimension,
 * types and number of data values. (not used by Mmg itself).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_SOLSATVERTICESSIZE(mesh,sol,nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: nsols\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: nentities\n
 * >     INTEGER, INTENT(IN)           :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Set_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol *sol,int nsols,
                                                   MMG5_int nentities, int *typSol);

/**
 * \brief Set the number of vertices, tetrahedra, prisms, triangles, quadrilaterals,
 * and edges of a mesh.
 * \param mesh pointer to the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nprism number of prisms.
 * \param nt number of triangles.
 * \param nquad number of quads.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * This function sets the number of vertices, tetrahedra, prisms, triangles,
 * quadrilaterals and edges of the mesh and allocates the associated arrays. If
 * called again, it will reset the whole mesh to reallocate it at the new size
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_MESHSIZE(mesh,np,ne,nprism,nt,nquad,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT)            :: np,ne,nprism,nt,nquad,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_meshSize(MMG5_pMesh mesh,MMG5_int np,MMG5_int ne,MMG5_int nprism,
                                          MMG5_int nt,MMG5_int nquad,MMG5_int na);

/* init structure data */
/**
 * \brief Set the coordinates of a single vertex.
 * \param mesh pointer to the mesh structure.
 * \param c0 coordinate of the vertex along the first dimension.
 * \param c1 coordinate of the vertex along the second dimension.
 * \param c2 coordinate of the vertex along the third dimension.
 * \param ref vertex reference.
 * \param pos position of the vertex in the mesh.
 * \return 1.
 *
 * This function sets the coordinates of a vertex \a c0, \a c1,\a c2 and reference \a ref
 * at position \a pos in mesh structure (from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_VERTEX(mesh,c0,c1,c2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: c0,c1,c2\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                                        double c2, MMG5_int ref,MMG5_int pos);

/**
 * \brief Set all vertex coordinates and references in a mesh structure
 *
 * \param mesh pointer to the mesh structure.
 * \param vertices  array of vertex coordinates in the order \f$[x_1, y_1, z_1, x_2, \ldots, z_N]\f$
 *   where \f$N\f$ is the number of vertices in the mesh.
 * \param refs  array of vertex references.
 *   The reference of vertex \f$i\f$ is stored in refs[\f$i-1\f$].
 * \return 1.
 *
 * This function sets the coordinates and references of all vertices in a mesh
 * structure. The number of vertices in the mesh must have been set before.
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_VERTICES(mesh,vertices,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * > !    REAL(KIND=8), INTENT(IN)      :: vertices(*)\n
 * > !    INTEGER(MMG5F_INT),INTENT(IN) :: refs(*)\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_vertices(MMG5_pMesh mesh, double *vertices,MMG5_int *refs);

/**
 * \brief set a single tetrahedron's vertices
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of tetrahedron.
 * \param v1 second vertex of tetrahedron.
 * \param v2 third vertex of tetrahedron.
 * \param v3 fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 * \param pos tetrahedron position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Assign the vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref to the tetrahedron at position \a pos in the mesh structure.
 * \a pos ranges from 1 to nb_tetra included.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TETRAHEDRON(mesh,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,v2,v3,pos\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_tetrahedron(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                            MMG5_int v2, MMG5_int v3, MMG5_int ref, MMG5_int pos);

/**
 * \brief Set the vertices and references of all tetrahedra in a mesh structure
 *
 * \param mesh pointer to the mesh structure.
 * \param tetra vertices of the tetras of the mesh given.
 * The vertices of the \f$i^{th}\f$ tetrahedron are given by
 *   tetra[(i-1)*4] to tetra[(i-1)*4+3] included.
 * \param refs array of the tetrahedra references.
 *   The references of the \f$i^{th}\f$ tetrahedron is given by refs[i-1].
 * \return 0 if failed, 1 otherwise.
 *
 * This function sets the vertices and references of all tetrahedra in a mesh.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_TETRAHEDRA(mesh,tetra,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*), INTENT(IN) :: tetra\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*), INTENT(IN) :: refs\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_tetrahedra(MMG5_pMesh mesh, MMG5_int *tetra,
                                            MMG5_int *refs);

/**
 * \brief Set the vertices and reference of a single prism in a mesh.
 *
 * \param mesh pointer to the mesh structure.
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
 * This function sets the six vertices \a v0, \a v1,\a v2,\a v3,\a v4,\a v5 and reference
 * \a ref for the prism at position \a pos in the mesh structure (from 1 to nb_prisms included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_PRISM(mesh,v0,v1,v2,v3,v4,v5,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,v2,v3,v4,v5,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_prism(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                      MMG5_int v2, MMG5_int v3, MMG5_int v4, MMG5_int v5, MMG5_int ref, MMG5_int pos);

/**
 * \brief Set the vertices and references of all prisms in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param prisms vertices of the prisms of the mesh
 *   Vertices of the \f$i^{th}\f$ prism are stored in prism[(i-1)*6]
 *   to prism[(i-1)*6+5] included.
 * \param refs array of the prisms references.
 *   The references of the \f$i^{th}\f$ prism is stored in refs[i-1].
 * \return 0 if failed, 1 otherwise.
 *
 * This function sets the vertices and references of all prisms in a mesh.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_PRISMS(mesh,prisms,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*), INTENT(IN) :: prisms\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*), INTENT(IN) :: refs\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_prisms(MMG5_pMesh mesh, MMG5_int *prisms,
                                        MMG5_int *refs);

/**
 * \brief Set the vertices and reference of a single triangle in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * This function defines a triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure (from 1 to nb_triangle included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TRIANGLE(mesh,v0,v1,v2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)         :: v0,v1,v2,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_triangle(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                          MMG5_int v2, MMG5_int ref,MMG5_int pos);

/**
 * \brief Set the vertices and references of all triangles in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param tria pointer to the array of the tria vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer to the array of the triangle references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \return 0 if failed, 1 otherwise.
 *
 * This function sets the vertices and references of all triangles in a mesh.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_TRIANGLES(mesh,tria,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * > !    INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN) :: tria\n
 * > !    INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN) :: refs\n
 * > !    INTEGER, INTENT(OUT)                        :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_triangles(MMG5_pMesh mesh, MMG5_int *tria, MMG5_int *refs);

/**
 * \brief Set the vertices and reference of a single quadrilateral in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of quadrilateral.
 * \param v1 second vertex of quadrilateral.
 * \param v2 third vertex of quadrilateral.
 * \param v3 fourth vertex of quadrilateral.
 * \param ref quadrilateral reference.
 * \param pos quadrilateral position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set a quadrilateral of vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref
 * at position \a pos in mesh structure (from 1 to nb_quadrangles included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_QUADRILATERAL(mesh,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_quadrilateral(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                               MMG5_int v2, MMG5_int v3, MMG5_int ref,MMG5_int pos);

/**
 * \brief Set the vertices and references of all quadrilaterals in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param quads pointer to the array of the quads vertices
 * Vertices of the \f$i^{th}\f$ quadra are stored in quads[(i-1)*3]\@3.
 * \param refs pointer to the array of the quadrilateral references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ quadra.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG3D_SET_QUADRILATERALS(mesh,quads,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * > !    INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN) :: quads\n
 * > !    INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN) :: refs\n
 * > !    INTEGER, INTENT(OUT)                        :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_quadrilaterals(MMG5_pMesh mesh, MMG5_int *quads, MMG5_int *refs);

/**
 * \brief Set the vertices and reference of a single edge in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first extremity of the edge.
 * \param v1 second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set edges of extremities \a v0, \a v1 and reference \a ref at
 * position \a pos in mesh structure  (from 1 to nb_edges included)
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_EDGE(mesh,v0,v1,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_edge(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int ref,MMG5_int pos);

/**
 * \brief Assign the "corner" attribute to a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set the "corner" attribute at vertex \a k. This affects how the vertex is
 * treated during remeshing.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_CORNER(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_corner(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "corner" attribute from a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Remove corner attribute from vertex \a k (from 1 to the number of vertices included).
 *
 * \remark Fortran interface
 *
 * >   SUBROUTINE MMG3D_UNSET_CORNER(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_corner(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set vertex \a k as required (\a k from 1 to nb_vertices included). This
 * prevents the remesher from moving the vertex.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDVERTEX(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_requiredVertex(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove required attribute from a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * This function removes the required attribute from vertex \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_REQUIREDVERTEX(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_requiredVertex(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to a tetrahedron.
 *
 * \param mesh pointer to the mesh structure.
 * \param k element index.
 * \return 1.
 *
 * Set element \a k as required (\a k from 1 to nb_tetra included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTETRAHEDRON(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_requiredTetrahedron(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from a tetrahedron.
 *
 * \param mesh pointer to the mesh structure.
 * \param k element index.
 * \return 1.
 *
 * Remove required attribute from element \a k (\a k from 1 to
 * nb_tetra included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_REQUIREDTETRAHEDRON(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_requiredTetrahedron(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to multiple tetrahedra.
 *
 * \param mesh pointer to the mesh structure.
 * \param reqIdx array of the indices of the required elements.
 * \param nreq number of required elements
 * \return 1.
 *
 * Determine which tetrahedra have the "required" attribute.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTETRAHEDRA(mesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)              :: nreq\n
 * >     INTEGER, INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_requiredTetrahedra(MMG5_pMesh mesh, MMG5_int *reqIdx, MMG5_int nreq);

/**
 * \brief Remove the "required" attribute from multiple tetrahedra.
 *
 * \param mesh pointer to the mesh structure.
 * \param reqIdx array of the indices of the required elements.
 * \param nreq number of required elements
 * \return 1.
 *
 * Remove required attribute from a list of Tetra whose indices are contained in
 * array \a reqIdx.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_REQUIREDTETRAHEDRA(mesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)              :: nreq\n
 * >     INTEGER, INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_requiredTetrahedra(MMG5_pMesh mesh, MMG5_int *reqIdx, MMG5_int nreq);

/**
 * \brief Assign the "required" attribute to a single triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as required (\a k from 1 to nb_tria included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_requiredTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from a single triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Remove required attribute from triangle \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_REQUIREDTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_requiredTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to multiple triangles.
 *
 * \param mesh pointer to the mesh structure.
 * \param reqIdx array of the indices of the required triangles.
 * \param nreq number of required triangles
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDTRIANGLES(mesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)              :: nreq\n
 * >     INTEGER, INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_requiredTriangles(MMG5_pMesh mesh, MMG5_int *reqIdx, MMG5_int nreq);

/**
 * \brief Remove the "required" attribute from multiple triangles.
 *
 * \param mesh pointer to the mesh structure.
 * \param reqIdx array of the indices of the required trias.
 * \param nreq number of required trias
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_REQUIREDTRIANGLES(mesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)              :: nreq\n
 * >     INTEGER, INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_requiredTriangles(MMG5_pMesh mesh, MMG5_int *reqIdx, MMG5_int nreq);

/**
 * \brief Assign the "parallel" attribute to a single triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as parallel (triangle at the interface between two
 * processors, ie, this triangle is required). (\a k from 1 to nb_tria included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_PARALLELTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_parallelTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "parallel" attribute from a single triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Remove parallel attribute from triangle \a k (i.e. triangles aren't preserved
 * anymore). (\a k from 1 to nb_tria included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_PARALLELTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_parallelTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "parallel" attribute to multiple triangles.
 *
 * \param mesh pointer to the mesh structure
 * \param parIdx array of indices of triangles
 * \param npar number of triangles
 * \return 1.
 *
 * Set the parallel triangles (triangles at the interface between processors, i.e.
 * these triangles are required).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_PARALLELTRIANGLES(mesh,parIdx,npar,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*),INTENT(IN) :: parIdx\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)              :: npar\n
 * >     INTEGER, INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_parallelTriangles(MMG5_pMesh mesh, MMG5_int *parIdx, MMG5_int npar);

/**
 * \brief Remove the "parallel" attribute from multiple triangles.
 *
 * \param mesh pointer to the mesh structure.
 * \param parIdx array of the indices of triangles
 * \param npar number of triangles
 * \return 1.
 *
 * Remove parallel attributes from triangles (i.e.
 * triangles aren't preserved anymore).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_PARALLELTRIANGLES(mesh,parIdx,npar,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*),INTENT(IN) :: parIdx\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)              :: npar\n
 * >     INTEGER, INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_parallelTriangles(MMG5_pMesh mesh, MMG5_int *parIdx, MMG5_int npar);

/**
 * \brief Assign the "ridge" attribute to a single edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set the "ridge" attribute at edge \a k. This affects how the remesher treats
 * the edge and the adjacent triangles or tetrahedra.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_RIDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_ridge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "ridge" attribute from a single edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Remove ridge attribute at edge \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_RIDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_ridge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to a single edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set edge \a k as required.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_REQUIREDEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_requiredEdge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from a single edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Remove required attribute from edge \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_UNSET_REQUIREDEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Unset_requiredEdge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Set the normal orientation at a single vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index
 * \param n0 x componant of the normal at vertex \a k.
 * \param n1 y componant of the normal at vertex \a k.
 * \param n2 z componant of the normal at vertex \a k.
 *
 * \return 1 if success.
 *
 * Set normal (n0,n1,n2) at vertex \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_NORMALATVERTEX(mesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8), INTENT(IN)      :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_normalAtVertex(MMG5_pMesh mesh, MMG5_int k, double n0, double n1,
                                                double n2) ;

/**
 * \brief Set a single element of a scalar solution structure.
 *
 * \param met pointer to the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 * (pos from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_SCALARSOL(met,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: s\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_scalarSol(MMG5_pSol met, double s,MMG5_int pos);

/**
 * \brief Set the values of all elements of a scalar solution structure.
 *
 * \param met pointer to the sol structure.
 * \param s array of the scalar solutions values.
 *   s[i-1] is the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar solutions at mesh vertices.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: s\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_scalarSols(MMG5_pSol met, double *s);

/**
 * \brief Set a single element of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in solution
 * structure. (pos from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_VECTORSOL(met,vx,vy,vz,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: vx,vy,vz\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz,
                                          MMG5_int pos);

/**
 * \brief Set all elements of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of the vectorial solutions
 * sols[3*(i-1)]\@3 is the solution at vertex i
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial solutions at mesh vertices
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: sols\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Set_vectorSols(MMG5_pSol met, double *sols);

/**
 * \brief Set a single element of a tensor solution structure.
 *
 * \param met pointer to the sol structure.
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
 * structure. (pos from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TENSORSOL(met,m11,m12,m13,m22,m23,m33,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                                         double m22,double m23, double m33, MMG5_int pos);

/**
 * \brief Set all elements of a tensor solution structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of the tensorial solutions.
 * sols[6*(i-1)]\@6 is the solution at vertex i
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values by array.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: sols\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Set_tensorSols(MMG5_pSol met, double *sols);

/**
 * \brief Set a single element of one out of multiple solution fields that are defined on vertices.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s solution(s) at mesh vertex \a pos.
 * \param pos index of the vertex on which we set the solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Set values of the solution at the ith field of the solution array and at
 * position \a pos (\a pos from 1 to \a nb_vertices included and \a i from 1 to \a
 * nb_sols).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_ITHSOL_INSOLSATVERTICES(sol,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)         :: pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);

/**
 * \brief Set all elements of one out of multiple solution fields that are defined on vertices.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s array of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Set values of the ith field of the solution array by array (\a i from 1 to \a
 * nb_sols).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_ITHSOLS_INSOLSATVERTICES(sol,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);

/**
 * \brief Finish providing mesh data without using the API functions.
 *
 * \param mesh pointer to the mesh structure.
 *
 * To mark as ended a mesh given without using the API functions (for example,
 * mesh given by mesh->point[i] = 0 ...). This function performs verifications,
 * e.g. to make sure that all tetrahedra are consistently oriented.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_HANDGIVENMESH(mesh)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT void MMG3D_Set_handGivenMesh(MMG5_pMesh mesh);

/* check init */
/**
 * \brief Check if the number of given entities match with mesh and sol size
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
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
 LIBMMG3D_EXPORT int MMG3D_Chk_meshData(MMG5_pMesh mesh, MMG5_pSol met);

/* functions to set parameters
 *
 * NOTE iparam and dparam are int rather than enum MMG3D_Param because
 * genheader cannot handle enums in function arguments; i.e. the Fortran
 * API will break.
 */

/**
 * \brief set an integer parameter of the remesher
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure (unused).
 * \param iparam integer parameter to set (see \ref MMG3D_Param for a
 *               list of parameters that can be set).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * This function sets the integer parameter \a iparam to value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_IPARAMETER(mesh,sol,iparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     MMG5_DATA_PTR_T                :: sol\n
 * >     INTEGER, INTENT(IN)            :: iparam\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: val\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_iparameter(MMG5_pMesh mesh,MMG5_pSol sol, int iparam,
                                           MMG5_int val);

/**
 * \brief set a real-valued parameter of the remesher
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure (unused).
 * \param dparam double parameter to set (see \ref MMG3D_Param for a
 *               list of parameters that can be set).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * This function sets the double parameter \a dparam to value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_DPARAMETER(mesh,sol,dparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     MMG5_DATA_PTR_T               :: sol\n
 * >     INTEGER, INTENT(IN)           :: dparam\n
 * >     REAL(KIND=8), INTENT(IN)      :: val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_dparameter(MMG5_pMesh mesh,MMG5_pSol sol, int dparam,
                                           double val);
/**
 * \brief set a local parameter
 *
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the sol structure
 * \param typ type of entity (triangle, edge,...)
 * \param ref reference of the entity
 * \param hmin minimal edge size
 * \param hmax maximal edge size
 * \param hausd Hausdorff distance
 * \return 0 if failed, 1 otherwise.
 *
 * Set local parameters: set the Hausdorff distance, minimum edge length, and
 * maximum edge length for all entities of type \a typ and reference \a ref.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_LOCALPARAMETER(mesh,sol,typ,ref,& \n
 * >                                       hmin,hmax,hausd,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: typ\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref\n
 * >     REAL(KIND=8), INTENT(IN)      :: hmin,hmax,hausd\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ,
                                                MMG5_int ref,double hmin,double hmax,double hausd);

/**
 * \brief Set the reference mapping for the elements of reference
 *  \a ref in level-set discretization mode.
 *
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the sol structure
 * \param ref input tetrahedron reference
 * \param split MMG5_MMAT_NoSplit if the entity must not be splitted, MMG5_MMAT_Split otherwise
 * \param rmin reference for the negative side after LS discretization
 * \param rplus reference for the positive side after LS discretization
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_MULTIMAT(mesh,sol,ref,split,rmin,rplus,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: split\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref,rmin,rplus\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 * With this function you can determine which references will be given to the
 * tetrahedra on both sides of the level set, after discretization. Negative and
 * positive here refer to volumes where the function is smaller or larger,
 * respectively, than the isovalue of the level set.
 */
 LIBMMG3D_EXPORT int  MMG3D_Set_multiMat(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int ref,int split,
                                         MMG5_int rmin, MMG5_int rplus);

/**
 * \brief Set a new level-set base reference.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param br new level-set base reference.
 * \return 0 if failed, 1 otherwise.
 *
 * Set a new level-set base reference of ref \a br in LS discretization
 * mode. Base references are boundary conditions to which implicit domains can
 * be attached. All implicit volumes that are not attached to listed base
 * references are deleted as spurious volumes by the \a rmc option.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_LSBASEREFERENCE(mesh,sol,br,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER(MMG5F_INT)            :: br\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG3D_EXPORT int  MMG3D_Set_lsBaseReference(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int br);


/** recover data */
/**
 * \brief Get the number of vertices, tetrahedra, prisms, triangles,
 * quadrilaterals and edges of the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param np pointer to the number of vertices.
 * \param ne pointer to the number of tetrahedra.
 * \param nprism pointer to the number of prisms.
 * \param nt pointer to the number of triangles.
 * \param nquad pointer to the number of quads.
 * \param na pointer to the number of edges.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_MESHSIZE(mesh,np,ne,nprism,nt,nquad,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT)            :: np,ne,nprism,nt,nquad,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_meshSize(MMG5_pMesh mesh, MMG5_int* np, MMG5_int* ne,MMG5_int *nprism, MMG5_int* nt,
                                          MMG5_int* nquad, MMG5_int* na);

/**
 * \brief Get the number of elements, dimension, and type of a solution structure.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typEntity pointer to the type of entities to which solutions
 *        are applied (see \ref MMG5_entities for possible values)
 * \param np pointer to the number of solutions.
 * \param typSol pointer to the type of the solutions
 *        (scalar, vectorial, ..., see \ref MMG5_type for possible values)
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER                       :: typEntity,typSol\n
 * >     INTEGER(MMG5F_INT)            :: np\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity,
                                         MMG5_int* np,int* typSol);

/**
 * \brief Get the number of elements, type, and dimensions of several solutions defined on vertices.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an array of sol structure.
 * \param nsols pointer to the number of solutions per entity.
 * \param nentities pointer to the number of solutions.
 * \param typSol array of size MMG5_NSOLS_MAX to store type of each solution
 *        (scalar, vectorial, ..., see \ref MMG5_type for possible values).
 *
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_SOLSATVERTICESSIZE(mesh,sol,nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER                       :: nsols\n
 * >     INTEGER(MMG5F_INT)            :: nentities\n
 * >     INTEGER                       :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol* sol,int *nsols,
                                                    MMG5_int* nentities,int* typSol);

/**
 * \brief Get the coordinates \a c0, \a c1,\a c2 and reference \a ref of the
 * next vertex of \a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 pointer to the coordinate of the vertex along the first
 * dimension.
 * \param c1 pointer to the coordinate of the vertex along the second
 * dimension.
 * \param c2 pointer to the coordinate of the vertex along the third
 * dimension.
 * \param ref pointer to the vertex reference.
 * \param isCorner pointer to the flag saying if vertex is corner.
 * \param isRequired pointer to the flag saying if vertex is required.
 * \return 1.
 *
 * This function retrieves the coordinates \a c0, \a c1,\a c2 and reference \a
 * ref of the next vertex of a mesh. It is meant to be used in a loop over all
 * vertices. When this function has been called as many times as there are
 * vertices, the internal loop counter will be reset. To obtain data for a
 * specific vertex, the \ref MMG3D_GetByIdx_vertex function can be used instead.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_VERTEX(mesh,c0,c1,c2,ref,isCorner,isRequired, &\n
 * >                               retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1,c2\n
 * >     INTEGER(MMG5F_INT)            :: ref\n
 * >     INTEGER                       :: isCorner,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2,
                                        MMG5_int* ref,int* isCorner, int* isRequired);

/**
 * \brief Get the coordinates and reference of a specific vertex in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 pointer to the coordinate of the vertex along the first dimension.
 * \param c1 pointer to the coordinate of the vertex along the second dimension.
 * \param c2 pointer to the coordinate of the vertex along the third dimension.
 * \param ref pointer to the vertex reference.
 * \param isCorner pointer to the flag saying if vertex is corner.
 * \param isRequired pointer to the flag saying if vertex is required.
 * \param idx index of vertex to get.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1, \a c2 and reference \a ref of
 * vertex \a idx of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GETBYIDX_VERTEX(mesh,c0,c1,c2,ref,isCorner,isRequired,idx,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1,c2\n
 * >     INTEGER(MMG5F_INT)            :: ref,idx\n
 * >     INTEGER                       :: isCorner,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_GetByIdx_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
                                            int* isCorner, int* isRequired,MMG5_int idx);

/**
 * \brief Get the coordinates and references of all vertices in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param vertices pointer to the array of coordinates.
 * The coordinates of the \f$i^{th}\f$ vertex are stored in
 * vertices[(i-1)*3]\@3.
 * \param refs pointer to the array of the vertex references.
 * The ref of the \f$i^th\f$ vertex is stored in refs[i-1].
 * \param areCorners pointer to the array of the flags saying if
 * vertices are corners.
 * areCorners[i-1]=1 if the \f$i^{th}\f$ vertex is corner.
 * \param areRequired pointer to the array of flags saying if vertices
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ vertex is required.
 * \return 1.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > ! SUBROUTINE MMG3D_GET_VERTICES(mesh,vertices,refs,areCorners,&\n
 * > !                               areRequired,retval)\n
 * > !   MMG5_DATA_PTR_T,INTENT(INOUT)          :: mesh\n
 * > !   REAL(KIND=8),DIMENSION(*), INTENT(OUT) :: vertices\n
 * > !   INTEGER(MMG5F_INT), DIMENSION(*)       :: refs\n
 * > !   INTEGER, DIMENSION(*)                  :: areCorners,areRequired\n
 * > !   INTEGER, INTENT(OUT)                   :: retval\n
 * > ! END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_vertices(MMG5_pMesh mesh, double* vertices, MMG5_int* refs,
                                          int* areCorners, int* areRequired);

/**
 * \brief Get the vertices and reference of the next tetrahedron in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the first vertex of tetrahedron.
 * \param v1 pointer to the second vertex of tetrahedron.
 * \param v2 pointer to the third vertex of tetrahedron.
 * \param v3 pointer to the fourth vertex of tetrahedron.
 * \param ref pointer to the tetrahedron reference.
 * \param isRequired pointer to the flag saying if tetrahedron is
 *  required.
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the vertices \a v0, \a v1, \a v2, \a v3 and reference
 * \a ref of the next tetrahedron of \a mesh. It is meant to be called in a loop
 * over all tetrahedra. When it has been called as many times as there are
 * tetrahedra, the internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TETRAHEDRON(mesh,v0,v1,v2,v3,ref,isRequired,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: v0,v1,v2,v3\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER                        :: isRequired\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_tetrahedron(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,
                                             MMG5_int* v3,MMG5_int* ref, int* isRequired);

/**
 * \brief Get the vertices and reference of all tetrahedra in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param tetra pointer to the array where the vertices are to be stored.
 * Vertices of the \f$i^{th}\f$ tetra are stored in tetra[(i-1)*4] to tetra[(i-1)*4+3]
 * \param refs pointer to the array of the tetrahedron references.
 * References of the \f$i^{th}\f$ tetra is stored in refs[i-1].
 * \param areRequired pointer to the array of the flags saying if the
 *  tetrahedra are required. areRequired[i-1]=1 if the \f$i^{th}\f$ tetrahedron
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_TETRAHEDRA(mesh,tetra,refs,areRequired,&\n
 * > !                                   retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*),INTENT(OUT) :: tetra\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)             :: refs\n
 * > !    INTEGER, DIMENSION(*)                        :: areRequired\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_tetrahedra(MMG5_pMesh mesh, MMG5_int* tetra,MMG5_int* refs,
                                            int* areRequired);

/**
 * \brief Get the vertices and reference of the next prism in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the first vertex of prism.
 * \param v1 pointer to the second vertex of prism.
 * \param v2 pointer to the third vertex of prism.
 * \param v3 pointer to the fourth vertex of prism.
 * \param v4 pointer to the fifth vertex of prism.
 * \param v5 pointer to the sixth vertex of prism.
 * \param ref pointer to the prism reference.
 * \param isRequired pointer to the flag saying if prism is
 *  required.
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the vertices \a v0, \a v1, \a v2, \a v3, \a v4, \a v5
 * and reference \a ref of the next prism of \a mesh. It is meant to be called
 * in a loop over all prisms. When it has been called as many times as there are
 * prisms, the internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_PRISM(mesh,v0,v1,v2,v3,v4,v5,ref,isRequired,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)   :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT) :: v0,v1,v2,v3,v4,v5\n
 * >     INTEGER(MMG5F_INT)              :: ref\n
 * >     INTEGER                         :: isRequired\n
 * >     INTEGER, INTENT(OUT)            :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_prism(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,
                                       MMG5_int* v3,MMG5_int* v4,MMG5_int* v5,MMG5_int* ref, int* isRequired);

/**
 * \brief Get the vertices and references of all prisms in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param prisms pointer to the array where the vertices are to be stored
 * Vertices of the \f$i^{th}\f$ prism are stored in prisms[(i-1)*6] to prisms[(i-1)*6+5].
 * \param refs pointer to the array of the prism references.
 * The reference of the \f$i^{th}\f$ prism is stored in refs[i-1].
 * \param areRequired pointer to the array of the flags saying if the
 *  prisms are required. areRequired[i-1]=1 if the \f$i^{th}\f$ prism
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_PRISMS(mesh,prisms,refs,areRequired,&\n
 * > !                              retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*),INTENT(OUT) :: prisms\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)             :: refs\n
 * > !    INTEGER, DIMENSION(*)                        :: areRequired\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_prisms(MMG5_pMesh mesh, MMG5_int* prisms,MMG5_int* refs,
                                        int* areRequired);

/**
 * \brief Get the vertices and reference of the next triangle in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the first vertex of triangle.
 * \param v1 pointer to the second vertex of triangle.
 * \param v2 pointer to the third vertex of triangle.
 * \param ref pointer to the triangle reference.
 * \param isRequired pointer to the flag saying if triangle is required.
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the vertices \a v0, \a v1, \a v2, and reference \a
 * ref of the next triangle of \a mesh. It is meant to be called in a loop over
 * all triangles. When it has been called as many times as there are triangles,
 * the internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TRIANGLE(mesh,v0,v1,v2,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)   :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT) :: v0,v1,v2\n
 * >     INTEGER(MMG5F_INT)              :: ref\n
 * >     INTEGER                         :: isRequired\n
 * >     INTEGER, INTENT(OUT)            :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_triangle(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref,
                                          int* isRequired);

/**
 * \brief Get the vertices and references of all triangles in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param tria pointer to the array where the vertices are to be stored
 * Vertices of the \f$i^{th}\f$ triangle are stored in tria[(i-1)*3] to tria[(i-1)*3+2].
 * \param refs pointer to the array where the references are to be stored.
 * refs[i-1] is the reference of the \f$i^{th}\f$ triangle.
 * \param areRequired pointer to array of the flags saying if triangles
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ tria
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface: (Commentated in order to allow to pass \%val(0)
 * instead of the refs or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_TRIANGLES(mesh,tria,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*),INTENT(OUT) :: tria\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)             :: refs\n
 * > !    INTEGER, DIMENSION(*)                        :: areRequired\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_triangles(MMG5_pMesh mesh, MMG5_int* tria, MMG5_int* refs,
                                           int* areRequired);

/**
 * \brief Get the vertices and reference of the next quadrilateral of the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the first vertex of quadrilateral.
 * \param v1 pointer to the second vertex of quadrilateral.
 * \param v2 pointer to the third vertex of quadrilateral.
 * \param v3 pointer to the fourth vertex of quadrilateral.
 * \param ref pointer to the quadrilateral reference.
 * \param isRequired pointer to the flag saying if quadrilateral is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get the vertices \a v0,\a v1,\a v2,\a v3 and reference \a ref of the next
 * quadrilateral of mesh. This function is meant to be called in a loop over all
 * quadrilaterals. When it has been called as many times as there are
 * quadrilaterals, the internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_QUADRILATERAL(mesh,v0,v1,v2,v3,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)   :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT) :: v0,v1,v2,v3,ref\n
 * >     INTEGER                         :: isRequired\n
 * >     INTEGER, INTENT(OUT)            :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_quadrilateral(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,MMG5_int* v3,
                                               MMG5_int* ref, int* isRequired);

/**
 * \brief Get the vertices and references of all quadrilaterals of the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param quads pointer to the array where the vertices will be stored.
 * Vertices of the \f$i^{th}\f$ quadrilateral are stored in quads[(i-1)*4] to quads[(i-1)*4+3].
 * \param refs pointer to the array of the quadrilaterals references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ quadrilateral.
 * \param areRequired pointer to array of the flags saying if quadrilaterals
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ quadrilateral
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface: (Commentated in order to allow to pass \%val(0)
 * instead of the refs or areRequired arrays)
 *
 * > !  SUBROUTINE MMG3D_GET_QUADRILATERALS(mesh,quads,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*),INTENT(OUT) :: quads\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)             :: refs\n
 * > !    INTEGER, DIMENSION(*)                        :: areRequired\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_quadrilaterals(MMG5_pMesh mesh, MMG5_int* quads, MMG5_int* refs,
                                                int* areRequired);

/**
 * \brief Get the vertices and reference of the next edge in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param e0 pointer to the first extremity of the edge.
 * \param e1 pointer to the second  extremity of the edge.
 * \param ref pointer to the edge reference.
 * \param isRidge pointer to the flag saying if the edge is ridge.
 * \param isRequired pointer to the flag saying if the edge is required.
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the extremities \a e0, \a e1 and reference \a ref of
 * next edge of \a mesh. It is meant to be called in a loop over all edges. When
 * it has been called as many times as there are edges in the mesh, the internal
 * edge counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_EDGE(mesh,e0,e1,ref,isRidge,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)   :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT) :: e0,e1\n
 * >     INTEGER(MMG5F_INT)              :: ref\n
 * >     INTEGER                         :: isRidge,isRequired\n
 * >     INTEGER, INTENT(OUT)            :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_edge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref,
                                      int* isRidge, int* isRequired);

/**
 * \brief Set the vertices and references of all edges in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param edges pointer to the array of edges.
 * The vertices of the \f$i^{th}\f$ edge should be given in edges[(i-1)*2] and edges[(i-1)*2+1].
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_EDGES(mesh,edges,refs,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: edges(*)\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: refs(*)\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Set_edges(MMG5_pMesh mesh, MMG5_int *edges, MMG5_int* refs);

/**
 * \brief Get the vertices and references of all edges in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param edges pointer to the array of edges.
 * The vertices of the \f$i^{th}\f$ edge are stored in edges[(i-1)*2] and edges[(i-1)*2+1].
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \param areRidges 1 if the edge is a ridge, 0 otherwise.
 * \param areRequired 1 if the edge is required, 0 otherwise.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_EDGES(mesh,edges,refs,areRidges,areRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: edges(*)\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: refs(*)\n
 * >     INTEGER, INTENT(OUT)           :: areRequired(*),areRidges(*)\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Get_edges(MMG5_pMesh mesh,MMG5_int *edges,MMG5_int* refs,
                                      int *areRidges,int *areRequired);

/**
 * \brief Get the normal orientation at a single mesh vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index
 * \param n0 x componant of the normal at vertex \a k.
 * \param n1 y componant of the normal at vertex \a k.
 * \param n2 z componant of the normal at vertex \a k.
 *
 * \return 1 if success.
 *
 * This function retrieves the normal (n0,n1,n2) at vertex \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_NORMALATVERTEX(mesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8)                  :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Get_normalAtVertex(MMG5_pMesh mesh, MMG5_int k, double *n0, double *n1,
                                               double *n2) ;

/**
 * \brief Get the quality measure of a single tetrahedron in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure (may be NULL for an isotropic metric).
 * \param k index of the tetrahedron for which we want to get the quality (from 1 to
 *        the number of tetrahedra included)
 * \return the computed quality or 0 in case of failure.
 *
 * This function returns the quality measure of tetrahedron \a k. Quality values
 * range from 0 (degenerate) to 1 (best attainable). The function returns 0
 * if the tetrahedron is flat or has a negative volume, and also if \a k is out
 * of range. In the latter case it will also print a diagnostic message to
 * standard output.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TETRAHEDRONQUALITY(mesh,met,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,met\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8), INTENT(OUT)     :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT double MMG3D_Get_tetrahedronQuality(MMG5_pMesh mesh, MMG5_pSol met, MMG5_int k);

/**
 * \brief Get the next element of a scalar solution structure defined at vertices.
 *
 * \param met pointer to the sol structure.
 * \param s pointer to the scalar solution value.
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the solution \a s of the next vertex of \a mesh. It
 * is meant to be called in a loop over all vertices. When it has been called as
 * many times as there are vertices in the mesh, the internal loop counter will
 * be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_SCALARSOL(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Get_scalarSol(MMG5_pSol met, double* s);

/**
 * \brief Get all elements of a scalar solution structure defined at vertices.
 *
 * \param met pointer to the sol structure.
 * \param s array of the scalar solutions at mesh vertices.
 * The solution at vertex i will be stored in s[i-1].
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int  MMG3D_Get_scalarSols(MMG5_pSol met, double* s);

/**
 * \brief Get the next element of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the next vector-valued element \f$(v_x,v_y,vz)\f$ of
 * the solution. It is meant to be called in a loop over all elements.  When it
 * has been called as many times as there are elements in the solution, the
 * internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_VECTORSOL(met,vx,vy,vz,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: vx,vy,vz\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz);

/**
 * \brief Get all elements of a vector solution structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of the solutions at mesh vertices. sols[3*(i-1)]\@3 is
 * the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vectorial solutions at mesh vertices
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Get_vectorSols(MMG5_pSol met, double* sols);

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
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the next element
 * \f$(m_{11},m_{12},m_{13},m_{22},m_{23},m_{33})\f$ of a tensor-valued solution
 * field.  It is meant to be called in a loop over all vertices. When it has
 * been called as many times as there are elements in the solution, the internal
 * loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TENSORSOL(met,m11,m12,m13,m22,m23,m33,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                                         double *m22,double *m23, double *m33);

/**
 * \brief Get all elements of a tensor solution field.
 *
 * \param met pointer to the sol structure.
 * \param sols array of the solutions at mesh vertices.
 * The solution at vertex \a i will be stored in sols[6*(i-1)] to sols[6*(i-1)+5].
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)           :: met\n
 * >     REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)                    :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Get_tensorSols(MMG5_pSol met, double *sols);

/**
 * \brief Get one out of several solutions at a specific vertex.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s solution(s) at mesh vertex \a pos. The required size
 *   of this array depends on the type of solution.
 * \param pos index of the vertex on which we get the solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Get values of the ith field of the solution array at vertex \a pos.
 * (pos from 1 to nb_vertices included and \a i from 1 to \a nb_sols).
 * The type of solution is inferred from \a sol.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_ITHSOL_INSOLSATVERTICES(sol,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)         :: pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);

/**
 * \brief Get one out of several solutions at all vertices in the mesh.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s array of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * This function retrieves the values of the solution at the ith field of the
 * solution array (\a i from 1 to \a nb_sols).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_ITHSOLS_INSOLSATVERTICES(sol,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Get_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);

/**
 * \brief Get the value of an integer parameter of the remesher.
 *
 * \param mesh pointer to the mesh structure.
 * \param iparam integer parameter to get (see \ref MMG3D_Param for a
 *               list of parameters that can be set).
 * \return The value of integer parameter.
 *
 * Get the value of integer parameter \a iparam.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_IPARAMETER(mesh,iparam,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: iparam\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_Get_iparameter(MMG5_pMesh mesh, MMG5_int iparam);

/**
 * \brief Add a tetrahedron to the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of tetrahedron.
 * \param v1 second vertex of tetrahedron.
 * \param v2 third vertex of tetrahedron.
 * \param v3 fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 *
 * This function adds a tetrahedron with vertices \a v0, \a v1, \a v2, \a v3 and reference
 * \a ref at the first available position of the mesh.
 *
 * \return 0 if unable to create the tetrahedron, the unit-offset index of the new tet if it
 * has strictly positive volume, a negative index if it has a zero or negative volume.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_ADD_TETRAHEDRON(mesh,v0,v1,v2,v3,ref,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,v2,v3,ref\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_Add_tetrahedron(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                             MMG5_int v2, MMG5_int v3, MMG5_int ref);

/**
 * \brief Add a vertex to the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 x coor of the new vertex
 * \param c1 y coor of the new vertex
 * \param c2 z coor of the new vertex
 * \param ref vertex reference.
 *
 * \return 0 if unable to create the vertex, the index of the new vertex
 * otherwise.
 *
 * This function adds a vertex with coordinates \a c0 \a c1 \a c2 and reference
 * \a ref at the first available position of the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_ADD_VERTEX(mesh,c0,c1,c2,ref,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)   :: mesh\n
 * >     REAL(KIND=8), INTENT(IN)        :: c0,c1,c2\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)  :: ref\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT) :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT MMG5_int  MMG3D_Add_vertex(MMG5_pMesh mesh, double c0, double c1,
                                             double c2, MMG5_int ref);

/* input/output functions */
/**
 * \brief Load a mesh (in .mesh/.mesb format) from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * Read mesh data.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_loadMesh(MMG5_pMesh mesh,const char *filename);

/**
 * \brief Load a mesh and possibly a solution in .msh format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * This function reads a mesh and 0 or 1 data fields in MSH file format (.msh
 * extension). We read only low-order vertices, edges, triangles, quadrangles,
 * tetrahedra and prisms.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADMSHMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_loadMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a mesh and possibly a solution in VTU (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * This function reads a mesh and 0 or 1 data field in VTU (VTK) file format (.vtu
 * extension). We read only low-order vertices, edges, tria, quadra, tetra and
 * prisms. Point and cell references must be stored in PointData or CellData
 * whose names contain the "medit:ref" keyword.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADVTUMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_loadVtuMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a mesh and multiple solutions in VTU (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * This functionreads a mesh and a list of data in VTU file format (.vtu extension). We read
 * only low-order vertices, edges, tria, quadra, tetra and prisms. Point and cell
 * references must be stored in PointData or CellData whose names contains the
 * "medit:ref" keyword.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADVTUMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG3D_EXPORT int MMG3D_loadVtuMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Load a mesh and possibly a solution from a file in VTK format.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * This function reads a mesh and 0 or 1 data fields in VTK file format (.vtu extension). We read
 * only low-order vertices, edges, tria, quadra, tetra and prisms. Point and cell
 * references must be stored in PointData or CellData whose names contain the
 * "medit:ref" keyword.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADVTKMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_loadVtkMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a mesh and multiple solutions from a file in VTK format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * Read mesh and a list of data in VTK file format (.vtu extension). We read
 * only low-order vertices, edges, tria, quadra, tetra and prisms. Point and cell
 * references must be stored in PointData or CellData whose names contains the
 * "medit:ref" keyword.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADVTKMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_loadVtkMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Load a mesh and all data from a file in MSH format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to a list of solution structures.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * Read mesh and a list of data in MSH file format (.msh extension). We read only
 * low-order vertices, edges, tria, quadra, tetra and prisms.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADMSHMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT  int MMG3D_loadMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Read mesh data in a format determined by the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient memory, file
 * format...), 1 if success.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADGENERICMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT int MMG3D_loadGenericMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,const char *filename);

/**
 * \brief Save a mesh in .mesh/.meshb format.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename pointer to the name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveMesh(MMG5_pMesh mesh, const char *filename);

/**
 * \brief Save a mesh in MSH format, ascii or binary depending on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data in MSH file format (.msh extension). Write binary
 * file for .mshb extension and ASCII for .msh one.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEMSHMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Save a mesh and data in MSH format, ascii or binary depending on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields (that are considered as solutions and
 * not metrics, thus, we do nothing over the ridge vertices) in MSH file format
 * (.msh extension).  Save file in ASCII format for .msh extension, in binary
 * format for .mshb one.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEMSHMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Save a mesh and optionally one solution in VTK format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data in Vtk file format (.vtk extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEVTKMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveVtkMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Save a mesh and multiple data fields in VTK format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields in Vtk file format (.vtk extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEVTKMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveVtkMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Save a mesh and optionally one data field in VTU format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data in vtu Vtk file format (.vtu extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEVTUMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveVtuMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Save a mesh and multiple data fields in VTU format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields in vtu Vtk file format (.vtu extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEVTUMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveVtuMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Save data in Tetgen's Triangle format.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the readed file.
 * \return 0 or -1 if fail, 1 otherwise.
 *
 * Save mesh data in Triangle (or equivalent to Tetgen in 3D) file format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVETETGENMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveTetgenMesh(MMG5_pMesh ,const char *);

/**
 * \brief Save mesh data in a file whose format depends on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEGENERICMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveGenericMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a metric field (or other solution).
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient
 * memory, file format...), 1 if successful.
 *
 * Load metric field. The solution file must contains only 1 solution: the
 * metric
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADSOL(mesh,met,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_loadSol(MMG5_pMesh mesh,MMG5_pSol met, const char *filename);

/**
 * \brief Load one or more solutions in a solution file in medit file format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solutions array
 * \param filename name of file.
 *
 * \return 0 if file is not found, -1 if fail for another reason (insufficient
 * memory, file format...), 1 if successful.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LOADALLSOLS(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_loadAllSols(MMG5_pMesh mesh,MMG5_pSol *sol, const char* filename);

/**
 * \brief Write isotropic or anisotropic metric.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVESOL(mesh,met,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveSol(MMG5_pMesh mesh,MMG5_pSol met, const char *filename);

/**
 * \brief Save 1 or more solutions in medit solution file format
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solutions array
 * \param filename name of the solution file.
 * \return 0 or -1 if fail, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SAVEALLSOLS(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_saveAllSols(MMG5_pMesh  mesh,MMG5_pSol *sol ,const char *filename);

/**
 * \brief Deallocate an array of solution fields
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an array of solution structure (that stores solution fields).
 * \return 1
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_Free_allSols(mesh,sol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Free_allSols(MMG5_pMesh mesh,MMG5_pSol *sol);

/* deallocations */
/**
 * \brief Deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend on the library function that you
 * have call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * \return 1 if success, 0 if fail
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Free_all(const int starter,...);

/**
 * \brief Structure deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend on the library function that you
 * have call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the \ref MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * \return 0 if fail, 1 if success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Free_structures(const int starter,...);

/**
 * \brief Structure deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend on the library function that you
 * have call.
 *
 * For the MMG3D_mmg3dlib function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the MMG3D_mmg3dls function, you need
 * to call the \ref MMG3D_Init_mesh function with the following arguments :
 * MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the \ref MMG3D_mmg3dmov function, you must call
 * : MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * \return 0 if fail, 1 if success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Free_names(const int starter,...);

/* library */
/**
 * \brief Main "program" for the mesh adaptation library.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol (metric) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DLIB(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_mmg3dlib(MMG5_pMesh mesh, MMG5_pSol met );

/**
 * \brief Main "program" for the level-set discretization library.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol (level-set) structure.
 * \param met pointer to a sol structure (metric), optional.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the level-set discretization library. If a metric \a met is
 * provided, use it to adapt the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DLS(mesh,sol,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     MMG5_DATA_PTR_T                :: met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_mmg3dls(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSol met );

/**
 * \brief Main program for the rigid-body movement library.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol (output metric) structure.
 * \param disp pointer to a sol (displacement for the lagrangian motion
 * mode) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DMOV(mesh,met,disp,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,disp\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_mmg3dmov(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol disp );

/** Tools for the library */
/**
 * \brief Print the default parameters values.
 *
 * \param mesh pointer to the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_DEFAULTVALUES(mesh,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_defaultValues(MMG5_pMesh mesh);

/**
 * \brief Store command-line arguments.
 *
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer to the mesh structure.
 * \param met pointer to a metric
 * \param sol pointer to a level-set or displacement
 * \return 1 if we want to run Mmg after, 0 if not or if fail.
 *
 * \remark no matching fortran function.
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol);

/**
 * \brief Read a file containing Local parameters (.mmg3d extension)
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
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
  LIBMMG3D_EXPORT int  MMG3D_parsop(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Print help for mmg3d options.
 *
 * \param prog pointer to the program name.
 * \param return 1 if success, 0 if fail.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_USAGE(prog,strlen0,retval)\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: prog\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_usage(char *prog);

/**
 * \brief Store the info structure in the mesh structure.
 *
 * \param mesh pointer to the mesh structure.
 * \param info pointer to the info structure.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_STOCKOPTIONS(mesh,info,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,info\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_stockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/**
 * \brief Recover the info structure stored in the mesh structure.
 *
 * \param mesh pointer to the mesh structure.
 * \param info pointer to the info structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_DESTOCKOPTIONS(mesh,info)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,info\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT void  MMG3D_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/** Checks */
/**
 * \brief Search invalid elements (in term of quality or edge length) in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure (metric).
 * \param sol pointer to the sol structure (ls or displacement).
 * \param critmin minimum quality for elements.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab array of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the MMG5_defsiz call), 1 for special storage (after this call).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_MMG3DCHECK(mesh,met,sol,critmin,lmin,lmax,eltab,&\n
 * >                               metridtyp,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)               :: mesh,met,sol\n
 * >     REAL(KIND=8), INTENT(IN)                     :: critmin,lmin,lmax\n
 * >     INTEGER(MMG5F_INT),DIMENSION(*), INTENT(OUT) :: eltab\n
 * >     INTEGER, INTENT(IN)                          :: metridtyp\n
 * >     INTEGER, INTENT(OUT)                         :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,double critmin,
                                       double lmin, double lmax, MMG5_int *eltab,int8_t metRidTyp);

/**
 * \brief List bad elements.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param critmin minimum quality for elements.
 * \param eltab pointer to the array of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the MMG5_defsiz call), 1 for special storage (after this call).
 *
 * Store elements which have worse quality than \a critmin in \a eltab,
 * \a eltab is allocated and could contain \a mesh->ne elements.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SEARCHQUA(mesh,met,critmin,eltab,metridtyp)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)               :: mesh,met\n
 * >     REAL(KIND=8), INTENT(IN)                     :: critmin\n
 * >     INTEGER(MMG5F_INT),DIMENSION(*), INTENT(OUT) :: eltab\n
 * >     INTEGER, INTENT(IN)                          :: metridtyp\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT void  MMG3D_searchqua(MMG5_pMesh mesh, MMG5_pSol met, double critmin,
                                        MMG5_int *eltab,int8_t metRidTyp);

/**
 * \brief List edges that are too short or too long.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param lmin minimum edge length.
 * \param lmax maximum ede length.
 * \param eltab array of invalid elements.
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage
 * (before the MMG5_defsiz call), 1 for special storage (after this call).
 *
 * \return 1 if success, 0 otherwise.
 *
 * Store in \a eltab elements which have edge lengths shorter than \a lmin
 * or longer than \a lmax, \a eltab is allocated and could contain
 * \a mesh->ne elements.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SEARCHLEN(mesh,met,lmin,lmax,eltab,metridtyp,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)               :: mesh,met\n
 * >     REAL(KIND=8), INTENT(IN)                     :: lmin,lmax\n
 * >     INTEGER(MMG5F_INT),DIMENSION(*), INTENT(OUT) :: eltab\n
 * >     INTEGER, INTENT(IN)                          :: metridtyp\n
 * >     INTEGER, INTENT(OUT)                         :: retval\n
 * >   END SUBROUTINE\n
 *
 *
 */
  LIBMMG3D_EXPORT int  MMG3D_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin,
                                       double lmax,MMG5_int *eltab,int8_t  metRidTyp);

/** Utils */
/**
 * \brief Return adjacent elements of a tetrahedron.
 *
 * \param mesh pointer to the mesh structure.
 * \param kel tetrahedron index.
 * \param listet pointer to the array of the 4 tetra adjacent to \a kel.
 * (the index is 0 if there is no adjacent)
 * \return 1.
 *
 * Find the indices of the 4 adjacent elements of tetrahedron \a
 * kel. \f$listet[i] = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_ADJATET(mesh,kel,listet,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)                :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)                :: kel\n
 * >     INTEGER(MMG5F_INT), DIMENSION(4), INTENT(OUT) :: listet\n
 * >     INTEGER, INTENT(OUT)                          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Get_adjaTet(MMG5_pMesh mesh,MMG5_int kel, MMG5_int listet[4]);

/**
 * \brief Compute the length of an edge according to the size prescription.
 *
 * \param ca pointer to the coordinates of the first edge's extremity.
 * \param cb pointer to the coordinates of the second edge's extremity.
 * \param ma pointer to the metric associated to the first edge's
 * extremity.
 * \param mb pointer to the metric associated to the second edge's
 * extremity.
 * \return edge length.
 *
 * Compute the length of edge \f$[ca,cb]\f$ (with \a ca and \a cb
 * coordinates of edge endpoints) according to the size
 * prescription.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_LENEDGCOOR(ca,cb,sa,sb,retval)\n
 * >     REAL(KIND=8), INTENT(IN)           :: ca,cb,sa,sb\n
 * >     REAL(KIND=8), INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT extern  double (*MMG3D_lenedgCoor)(double *ca,double *cb,double *sa,double *sb);

/**
 * \brief Create array of adjacency.
 *
 * \param mesh pointer to the mesh structure.
 * \param pack we pack the mesh at function begining if \f$pack=1\f$.
 * \return 0 if failed, 1 otherwise.
 *
 * Create array of adjacency. Set pack variable to 0 for a compact
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
  LIBMMG3D_EXPORT int  MMG3D_hashTetra(MMG5_pMesh mesh, int pack);

/**
 * \brief Compute isotropic size map according to the mean of the length of the
 * edges passing through a vertex.
 *
 * \param mesh pointer to the mesh structure
 * \param met pointer to the sol structure
 * \return 1 if success
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_DOSOL(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT extern int (*MMG3D_doSol)(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Compute a constant size map according to the hsiz, hmin and hmax parameters.
 *
 * \param mesh pointer to the mesh structure
 * \param met pointer to the sol structure
 * \return 1 if success
 *
 * This function computes a constant size map according to mesh->info.hsiz,
 * mesh->info.hmin and mesh->info.hmax. It updates these 3 values if not
 * compatible.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SET_CONSTANTSIZE(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Swap the m22 and m23 values of the metric.
 *
 * \param mesh pointer to the mesh structure
 * \param met pointer to the sol structure
 * \return 1 if success
 *
 * Switch the m22 and m23 value of the metric to allow to pass from the API
 * storage to the medit storage.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SWITCH_METRICSTORAGE(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_switch_metricStorage(MMG5_pMesh mesh, MMG5_pSol met);


/** To associate function pointers without calling MMG3D_mmg3dlib */
/**
 * \brief Set function pointers for caltet, lenedg, lenedgCoor defsiz, gradsiz...
 * depending if the metric that was read is anisotropic or isotropic
 *
 * \param mesh pointer to the mesh structure (unused).
 * \param met pointer to the sol structure (unused).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_SETFUNC(mesh,met)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG3D_EXPORT void  MMG3D_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Get the number of non-boundary triangles.
 *
 * \param mesh pointer to the mesh structure.
 * \param nb_tria pointer to the number of non-boundary triangles.
 * \return 0 if failed, 1 otherwise.
 *
 * Get the number of non-boundary triangles (for DG methods for example).
 * A triangle is
 * boundary if it is located at the interface of 2 domains with different
 * references or if it belongs to one tetra only.
 * Append these triangles to the list of triangles.
 *
 * \warning reallocates the triangle array and appends the internal triangles.
 * This may modify the behaviour of other functions.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_NUMBEROFNONBDYTRIANGLESS(mesh,nb_tria,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)   :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT) :: nb_tria\n
 * >     INTEGER, INTENT(OUT)            :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Get_numberOfNonBdyTriangles(MMG5_pMesh mesh, MMG5_int* nb_tria);

/**
 * \brief Get vertices and reference of a non-boundary triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the firts vertex of the triangle
 * \param v1 pointer to the second vertex of the triangle.
 * \param v2 pointer to the third vertex of the triangle.
 * \param ref pointer to the triangle reference.
 * \param idx index of the non-boundary triangle to get (between 1 and nb_tria)
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and reference \a ref of the idx^th non-boundary
 * triangle (for DG methods for example). A tria is boundary if it is located at
 * the interface of 2 domains with different references or if it belongs to one
 * tetra only.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_NONBDYTRIANGLE(mesh,v0,v1,v2,ref,idx,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: v0,v1,v2\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: idx\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Get_nonBdyTriangle(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref, MMG5_int idx);

/**
 * \brief Get a tetrahedron given one of its triangles and the index by which it
 * refers to this triangle (DEPRECATED).
 *
 * \param mesh pointer to the mesh structure.
 * \param ktri index of the boundary triangle.
 * \param ktet pointer to an integer that will contains the tetra index.
 * \param iface pointer to the triangle in \a ktet.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Fill \a ktet by the index of a tetrahedron to which belongs a boundary triangle
 * and \a iface by the index of the triangle in the tetra.
 *
 * \warning will be deprecated in release 5.5
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TETFROMTRIA(mesh,ktri,ktet,iface,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)              :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)           :: ktri\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT)          :: ktet\n
 * >     INTEGER, INTENT(OUT)                     :: iface\n
 * >     INTEGER, INTENT(OUT)                     :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Get_tetFromTria(MMG5_pMesh mesh, MMG5_int ktri, MMG5_int *ktet, int *iface);

/**
 * \brief Get two tetrahedra given a triangle and face indices.
 *
 * \param mesh pointer to the mesh structure.
 * \param ktri index of the boundary triangle.
 * \param ktet array of size 2 that will contain the indices of the tetra
 * (filled by the function).
 * \param iface pointer to an array of size 2 that will contains the indices
 * of the faces of the tetras \a ktet[i] that corresponds to the boundary tria
 * \a ktri.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Fill \a ktet by the indices of the tetrahedra that have a boundary triangle
 * and \a iface by the indices of the faces of the tetras that correspond to the
 * triangle. Fill ktet[1] and iface[1] by 0 if the triangle belongs to 1 tetrahedron only.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_GET_TETSFROMTRIA(mesh,ktri,ktet,iface,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)                  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)               :: ktri\n
 * >     INTEGER(MMG5F_INT), DIMENSION(2), INTENT(OUT):: ktet\n
 * >     INTEGER, DIMENSION(2), INTENT(OUT)           :: iface\n
 * >     INTEGER, INTENT(OUT)                         :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Get_tetsFromTria(MMG5_pMesh mesh, MMG5_int ktri, MMG5_int ktet[2], int iface[2]);

/**
 * \brief Compute the real eigenvalues and eigenvectors of a symmetric matrix
 *
 * \param m upper part of a symmetric matrix diagonalizable in |R
 * \param lambda array of the metric eigenvalues
 * \param vp array of the metric eigenvectors
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
 * >   SUBROUTINE MMG3D_COMPUTE_EIGENV(m,lambda,vp,retval)\n
 * >     REAL(KIND=8), INTENT(IN)         :: m(*)\n
 * >     REAL(KIND=8), INTENT(OUT)        :: lambda(*),vp(*)\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Compute_eigenv(double m[6],double lambda[3],double vp[3][3]);

/**
 * \brief Clean data (triangles and edges) linked to isosurface.
 *
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if successful, 0 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_CLEAN_ISOSURF(mesh,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)      :: mesh\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT int MMG3D_Clean_isoSurf(MMG5_pMesh mesh);

/**
 * \brief Free the solution structure of a given mesh.
 *
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the solution structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG3D_FREE_SOLUTIONS(mesh,sol)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG3D_EXPORT void MMG3D_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol);


 /**
  * Set common pointer functions between mmgs and mmg3d to the matching mmg3d
  * functions.
  */
  LIBMMG3D_EXPORT void MMG3D_Set_commonFunc(void);
#ifdef __cplusplus
}
#endif

#endif
