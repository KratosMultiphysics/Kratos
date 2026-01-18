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
 * This file defines the C and Fortran headers of the mmg2d API, and
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
 * \file mmg2d/libmmg2d.h
 * \ingroup API
 * \brief API headers and documentation for the mmg2d library
 * \author Cecile Dobrzynski and Algiane Froehly (Bx INP/Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * These are the API functions for the mmg2d library. These functions allow to
 * load and save meshes and data defined on meshes; add, extract, or modify mesh
 * data; and to call the library functions that perform meshing, remeshing,
 * level-set discretization, and Lagrangian motion.
 *
 * Meshes are here defined in terms of vertices and two-dimensional objects:
 * triangles and quadrangles. Edges can also be represented. All of these \a
 * entities can have a \a reference: an integer value that can serve as a group
 * identifier. In addition mesh entities can have \a attributes such as
 * "required" or "corner".
 *
 * Data defined on meshes can be for example functions that are meant for
 * level-set discretization, metric tensors that will govern edge lengths, and
 * vector fields governing lagrangian motion. These data can be scalar, vector,
 * or (symmetric) tensor-valued; and there can be more than one data item
 * associated with a mesh entity. These data are often referred to as \a
 * solutions.
 *
 * Four of the functions here are referred to as "programs", because they
 * perform the tasks for which Mmg is meant: (re)meshing, level-set discretization
 * and Lagrangian motion. The other functions merely serve to load and save data
 * and to perform pre- and post-processing. These programs actually behave much
 * like independent programs: they send diagnostic output to stdout and in rare
 * cases they may call the exit() function.
 *
 */

#ifndef MMG2DLIB_H
#define MMG2DLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mmg/common/libmmgtypes.h"
#include "mmg/mmg2d/mmg2d_export.h"

/**
 * Maximum array size when storing adjacent vertices (or ball) of a vertex.
 */
#define MMG2D_LMAX   1024

/**
 * \enum MMG2D_Param
 * \brief Input parameters for mmg library.
 *
 * These are the input parameters for mmg2d library functions. Options prefixed by \a
 * MMG2D_IPARAM require integer values and options prefixed by \a
 * MMG2D_DPARAM require real values. They can be set with the
 * \ref MMG2D_Set_iparameter and \ref MMG2D_Set_dparameter functions,
 * respectively.
 *
 */
  enum MMG2D_Param {
    MMG2D_IPARAM_verbose,           /*!< [-1..10], Level of verbosity */
    MMG2D_IPARAM_mem,               /*!< [n/-1], Max memory size in Mbytes or keep the default value */
    MMG2D_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
    MMG2D_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
    MMG2D_IPARAM_iso,               /*!< [1/0], Enable level-set discretization */
    MMG2D_IPARAM_isosurf,           /*!< [1/0], Enable level-set discretization on the surface part only */
    MMG2D_IPARAM_opnbdy,            /*!< [1/0], Preserve edges at the interface of 2 domains with same reference */
    MMG2D_IPARAM_lag,               /*!< [-1/0/1/2], Enable Lagrangian motion */
    MMG2D_IPARAM_3dMedit,           /*!< [0/1/2], Read/write 2D mesh in 3D (Medit only). out if val=1 in/out if val=2 */
    MMG2D_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
    MMG2D_IPARAM_noinsert,          /*!< [1/0], Avoid/allow vertex insertion */
    MMG2D_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
    MMG2D_IPARAM_nomove,            /*!< [1/0], Avoid/allow vertex relocation */
    MMG2D_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
    MMG2D_IPARAM_nreg,              /*!< [0/1], Enable normal regularization */
    MMG2D_IPARAM_xreg,              /*!< [0/1], Enable regularization by moving vertices */
    MMG2D_IPARAM_numsubdomain,      /*!< [0/n], Save only the subdomain n (0==all subdomains) */
    MMG2D_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
    MMG2D_IPARAM_numberOfLSBaseReferences,   /*!< [n], Number of base references for bubble removal */
    MMG2D_IPARAM_numberOfMat,                /*!< [n], Number of materials in level-set mode */
    MMG2D_IPARAM_anisosize,                 /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
    MMG2D_IPARAM_nosizreq,          /*!< [0/1], Allow/avoid overwriting of sizes at required vertices (advanced usage) */
    MMG2D_DPARAM_angleDetection,    /*!< [val], Threshold for angle detection */
    MMG2D_DPARAM_hmin,              /*!< [val], Minimal edge length */
    MMG2D_DPARAM_hmax,              /*!< [val], Maximal edge length */
    MMG2D_DPARAM_hsiz,              /*!< [val], Constant edge length */
    MMG2D_DPARAM_hausd,             /*!< [val], Global Hausdorff distance (on all the boundary surfaces of the mesh) */
    MMG2D_DPARAM_hgrad,             /*!< [val], Global gradation */
    MMG2D_DPARAM_hgradreq,          /*!< [val], Global gradation on required entites (advanced usage) */
    MMG2D_DPARAM_ls,                /*!< [val], Function value where the level set is to be discretized */
    MMG2D_DPARAM_xreg,              /*!< [val], Relaxation parameter for coordinate regularization (0<val<1) */
    MMG2D_DPARAM_rmc,               /*!< [-1/val], Remove small disconnected components in level-set mode */
    MMG2D_IPARAM_nofem,             /*!< [1/0], Do not attempt to make the mesh suitable for finite-element computations */
    MMG2D_IPARAM_isoref,            /*!< [0/n], Iso-surface boundary material reference */
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
 * \return 1 on success, 0 on failure
 *
 * For the MMG2D_mmgslib function, you need
 * to call the \ref MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the \ref MMG2D_mmgsls function, you need
 * to call the \a MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here, \a your_mesh is a \ref MMG5_pMesh, \a your_metric and \a your_level_set
 * are \ref MMG5_pSol.
 *
 * \remark No fortran interface, to allow variadic arguments.
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Init_mesh(const int starter,...);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 *
 * \brief Initialize file names to their default values.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_INIT_FILENAMES(mesh,sol)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT void  MMG2D_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);

/**
 * \param mesh pointer to the mesh structure.
 *
 * \brief Initialize the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_INIT_PARAMETERS(mesh)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT void  MMG2D_Init_parameters(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer to the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * \brief Set the name of the input mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_INPUTMESHNAME(mesh,meshin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin);

/**
 * \param mesh pointer to the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * \brief Set the name of the output mesh file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_OUTPUTMESHNAME(mesh,meshout,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshout\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * \brief Set the name of the input solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_INPUTSOLNAME(mesh,sol,solin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param solout name of the output solution file.
 * \return 0 on failure, 1 otherwise.
 *
 *  \brief Set the name of the output solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_OUTPUTSOLNAME(mesh,sol,solout,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solout\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout);

/**
 * \param mesh pointer to the mesh structure.
 * \param fparamin name of the input parameter file.
 * \return 1.
 *
 * \brief Set the name of the input parameter file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_INPUTPARAMNAME(mesh,fparamin,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: fparamin\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_inputParamName(MMG5_pMesh mesh, const char* fparamin);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure (unused).
 * \param iparam integer parameter to set (see \a MMG2D_Param structure).
 * \param val value for the parameter.
 * \return 0 on failure, 1 otherwise.
 *
 * \brief Set integer parameter \a iparam to value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_IPARAMETER(mesh,sol,iparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     MMG5_DATA_PTR_T               :: sol\n
 * >     INTEGER, INTENT(IN)           :: iparam\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, MMG5_int val);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param dparam double parameter to set (see \a MMG2D_Param structure).
 * \param val value of the parameter.
 * \return 0 on failure, 1 otherwise.
 *
 * \brief Set double parameter \a dparam to value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_DPARAMETER(mesh,sol,dparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     MMG5_DATA_PTR_T               :: sol\n
 * >     INTEGER, INTENT(IN)           :: dparam\n
 * >     REAL(KIND=8), INTENT(IN)      :: val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val);

/**
 * \brief Set local parameters.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param hmin minimal edge size.
 * \param hmax maximal edge size.
 * \param hausd value of the Hausdorff number.
 * \return 0 on failure, 1 otherwise.
 *
 * Set local parameters: set the Hausdorff distance, minimal desired edge
 * length, and maximal desired edge length for all entities of type \a typ and
 * reference \a ref.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_LOCALPARAMETER(mesh,sol,typ,ref,& \n
 * >                                       hmin,hmax,hausd,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: typ\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref\n
 * >     REAL(KIND=8), INTENT(IN)      :: hmin,hmax,hausd\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_localParameter(MMG5_pMesh mesh, MMG5_pSol sol, int typ,
                                                MMG5_int ref,double hmin,double hmax,double hausd);

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param ref input tetra reference.
 * \param split \ref MMG5_MMAT_NoSplit if the entity must not be split, \ref MMG5_MMAT_Split otherwise
 * \param rmin reference for the negative side after LS discretization
 * \param rplus reference for the positive side after LS discretization
 * \return 0 on failure, 1 otherwise.
 *
 * \brief Set the reference mapping for the elements of ref \a ref in LS discretization mode.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_MULTIMAT(mesh,sol,ref,split,rmin,rplus,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref,rmin,rplus\n
 * >     INTEGER, INTENT(IN)           :: split\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 * With this function you can determine which references will be given to the
 * triangles on both sides of the level set, after discretization. Negative and
 * positive here refer to areas where the function is smaller or larger,
 * respectively, than the isovalue of the level set.
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_multiMat(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int ref,int split,
                                          MMG5_int rmin, MMG5_int rplus);

/**
 * \brief Set level-set base reference.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param br new level-set base reference.
 * \return 0 on failure, 1 otherwise.
 *
 * Set a new level-set base reference of ref \a br in ls discretization
 * mode. Based references are boundary conditions to which implicit domain can
 * be attached. All implicit volumes that are not attached to listed base
 * references are deleted as spurious volumes by the \a rmc option.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_LSBASEREFERENCE(mesh,sol,br,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: br\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG2D_EXPORT int  MMG2D_Set_lsBaseReference(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int br);


/* init structure datas */
/**
 * \brief Set the numbers of entities in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param np number of vertices.
 * \param nt number of triangles.
 * \param nquad number of quads.
 * \param na number of edges.
 * \return 0 on failure, 1 otherwise.
 *
 * Set the number of vertices, triangles, quadrangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_MESHSIZE(mesh,np,nt,nquad,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT)            :: np,nt,nquad,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_meshSize(MMG5_pMesh mesh, MMG5_int np, MMG5_int nt, MMG5_int nquad, MMG5_int na);

/**
 * \brief Set the size and type of a solution field.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typEntity type of solutions entities (vertices, triangles, ..., 
 *        see \ref MMG5_entities for possible values).
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial, ..., see \ref MMG5_type for possible values).
 * \return 0 on failure, 1 otherwise.
 *
 * Initialize an array of solution field: set dimension, types and number of
 * data.
 * To use to initialize an array of solution fields (not used by Mmg itself).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: typEntity\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: np\n
 * >     INTEGER, INTENT(IN)           :: typSol\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity,
                                        MMG5_int np, int typSol);

/**
 * \brief Initialize an array of solutions field defined at vertices
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an allocatable sol structure.
 * \param nsols number of solutions per entity
 * \param nentities number of entities
 * \param typSol    Array of size nsol listing the type of the solutions
 *                  (scalar, vectorial, ..., see \ref MMG5_type for possible values).
 * \return 0 on failure, 1 otherwise.
 *
 * Initialize an array of solutions field defined at vertices: set dimension,
 * types and number of data.
 * To use to initialize an array of solution fields (not used by Mmg itself).
 *
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_SOLSATVERTICESSIZE(mesh,sol,nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(IN)           :: nsols\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: nentities\n
 * >     INTEGER, INTENT(IN)           :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol *sol,int nsols,
                                                   MMG5_int nentities, int *typSol);

/**
 * \brief Set the coordinates and reference of a single vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 coordinate of the vertex along the first dimension.
 * \param c1 coordinate of the vertex along the second dimension.
 * \param ref vertex reference.
 * \param pos position of the vertex in the mesh.
 * \return 1 on success, 0 in case of failure
 *
 * Set vertex of coordinates \a c0, \a c1 and reference \a ref
 * at position \a pos in mesh structure (from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_VERTEX(mesh,c0,c1,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: c0,c1\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                                        MMG5_int ref,MMG5_int pos);

/**
 * \brief Set the coordinates and references of all vertices in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param vertices array of vertex coordinates in the order \f$[x_1, y_1, x_2, \ldots, y_N]\f$
 * \param refs array of vertex references \f$[r_1, r_2, \ldots, r_N]\f$
 * \return 1.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs array)
 *
 * > !  SUBROUTINE MMG2D_SET_VERTICES(mesh,vertices,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)              :: mesh\n
 * > !    REAL(KIND=8), DIMENSION(*),INTENT(IN)      :: vertices\n
 * > !    INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN):: refs\n
 * > !    INTEGER, INTENT(OUT)                       :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_vertices(MMG5_pMesh mesh, double *vertices,MMG5_int *refs);

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
 * \remark Fortran interface
 *
 * >   SUBROUTINE MMG2D_SET_CORNER(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_corner(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "corner" attribute from a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Remove corner attribute from vertex \a k (from 1 to the number of vertices
 * included).
 *
 * \remark Fortran interface
 *
 * >   SUBROUTINE MMG2D_UNSET_CORNER(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Unset_corner(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Assign the "required" attribute to a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set vertex \a k as required (\a k from 1 to the number of vertices
 * included). This prevents the remesher from moving the vertex.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_REQUIREDVERTEX(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_requiredVertex(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from a vertex.
 *
 * \param mesh pointer to the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * This function removes the "required" attribute from vertex \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_UNSET_REQUIREDVERTEX(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Unset_requiredVertex(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Set the vertices and reference of a single triangle in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * This function defines a triangle with vertices \a v0, \a v1, \a v2 and
 * reference \a ref at position \a pos in mesh structure (from 1 to the number
 * of triangles included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_TRIANGLE(mesh,v0,v1,v2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,v2,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_triangle(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                         MMG5_int v2, MMG5_int ref, MMG5_int pos);

/**
 * \brief Set the vertices and references of all triangles in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param tria pointer to the array of the triangle's vertices.
 * The vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer to the array of the triangle references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \return 0 on failure, 1 otherwise.
 *
 * This function sets the vertices and references of all triangles in a mesh.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs array)
 * > !  SUBROUTINE MMG2D_SET_TRIANGLES(mesh,tria,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)               :: mesh\n
 * > !    INTEGER(MMG5F_INT),DIMENSION(*), INTENT(IN) :: tria,refs\n
 * > !    INTEGER, INTENT(OUT)                        :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_triangles(MMG5_pMesh mesh, MMG5_int *tria, MMG5_int *refs);

/**
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * \brief Give triangle \a k the "required" attribute.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_REQUIREDTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_requiredTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \param mesh pointer to the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * \brief Remove the "required" attribute from triangle \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_UNSET_REQUIREDTRIANGLE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Unset_requiredTriangle(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Set the vertices and reference of a single quadrangle in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of quadrangle.
 * \param v1 second vertex of quadrangle.
 * \param v2 third vertex of quadrangle.
 * \param v3 fourth vertex of quadrangle.
 * \param ref quadrangle reference.
 * \param pos quadrangle position in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * Define a quadrangle with vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref at position \a pos in mesh structure (from 1 to nb_quad included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_QUADRILATERAL(mesh,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_quadrilateral(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
                                               MMG5_int v2, MMG5_int v3, MMG5_int ref, MMG5_int pos);

/**
 * \brief Set the vertices and references of all quadrangles in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param quadra vertices of the quadrangles of the mesh
 * Vertices of the \f$i^{th}\f$ quadrangle are stored in quadra[(i-1)*4]\@4.
 * \param refs array of references.
 * The reference of the \f$i^{th}\f$ quadrangle is stored in refs[i-1].
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE MMG2D_SET_QUADRILATERALS(mesh,quadra,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*), INTENT(IN) :: quadra,refs\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_quadrilaterals(MMG5_pMesh mesh, MMG5_int *quadra,
                                                MMG5_int *refs);

/**
 * \brief Define a single edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 first vertex of edge.
 * \param v1 second vertex of edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * Define an edge with vertices \a v0, \a v1 and reference \a ref at position \a
 * pos in the mesh structure (\a pos from 1 to the number of edges included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_EDGE(mesh,v0,v1,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: v0,v1\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos,ref\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_edge(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int ref, MMG5_int pos);

/**
 * \brief Set the vertices and references of all edges in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param edges pointer to the array of edges.
 * Vertices of the \f$i^{th}\f$ edge are stored in edge[(i-1)*2]\@2.
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_EDGES(mesh,edges,refs,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: edges(*),refs(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_edges(MMG5_pMesh mesh, MMG5_int *edges, MMG5_int* refs);

/**
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * \brief Give edge \a k the "required" attribute.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_REQUIREDEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_requiredEdge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Remove the "required" attribute from edge \a k.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_UNSET_REQUIREDEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Unset_requiredEdge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Give edge \a k the "parallel" attribute.
 *
 * \param mesh pointer to the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_PARALLELEDGE(mesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_parallelEdge(MMG5_pMesh mesh, MMG5_int k);

/**
 * \brief Set a single value of a sol structure.
 *
 * \param met pointer to the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure.
 * (pos from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_SCALARSOL(met,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: s\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_scalarSol(MMG5_pSol met, double s, MMG5_int pos);

/**
 * \brief Set all values of a scalar sol structure.
 *
 * \param met pointer to the sol structure.
 * \param s array of scalar solutions values.
 * s[i-1] is the solution at vertex i.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: s\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_scalarSols(MMG5_pSol met, double *s);

/**
 * \brief Set a single vector value in a sol structure.
 *
 * \param met pointer to the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param pos position of the solution in the mesh (begins at 1).
 * \return 0 on failure, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y)\f$ at position \a pos in solution
 * structure. ( pos from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_VECTORSOL(met,vx,vy,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: vx,vy\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_vectorSol(MMG5_pSol met, double vx,double vy,
                                          MMG5_int pos);

/**
 * \brief Set all values in a vector sol structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of vectorial solutions
 * sols[2*(i-1)]\@2 is the solution at vertex i
 * \return 0 on failure, 1 otherwise.
 *
 * Set vectorial solutions at mesh vertices
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)        :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN):: sols\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_vectorSols(MMG5_pSol met, double *sols);

/**
 * \brief Set a single element of a tensor sol structure.
 *
 * \param met pointer to the sol structure.
 * \param m11 value at position (1,1) in the solution tensor.
 * \param m12 value at position (1,2) in the solution tensor.
 * \param m22 value at position (2,2) in the solution tensor.
 * \param pos position of the solution in the mesh.
 * \return 0 on failure, 1 otherwise.
 *
 * Set tensor value \a s at position \a pos in solution structure
 * (\a pos from 1 to nb_vertices included).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_TENSORSOL(met,m11,m12,m22,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(IN)      :: m11,m12,m22\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_tensorSol(MMG5_pSol met, double m11, double m12, double m22,
                                          MMG5_int pos);

/**
 * \brief Set all elements of a tensor sol structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of tensorial solutions.
 * sols[3*(i-1)]\@3 is the solution at vertex i
 * \return 0 on failure, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: met\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: sols\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_tensorSols(MMG5_pSol met, double *sols);

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
 * Set values of the solution at the ith field of the solution array.
 * (\a pos from 1 to \a nb_vertices included and \a i from 1 to \a nb_sols).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_ITHSOL_INSOLSATVERTICES(sol,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)         :: pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Set_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);

/**
 * \brief Set all elements of one out of multiple solution fields that are defined on vertices.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s array of solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[2*(k-1)]\@2 for a vectorial solution
 * and s[3*(k-1)]\@3 for a tensor solution.
 *
 * \return 0 on failure, 1 otherwise.
 *
 * Set scalar values of the solution at the ith field of the solution array.
 * (\a i from 1 to \a nb_sols)
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_ITHSOLS_INSOLSATVERTICES(sol,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);

/** recover datas */
/**
 * \param mesh pointer to the mesh structure.
 * \param np pointer to the number of vertices.
 * \param nt pointer to the number of triangles.
 * \param nquad pointer to the number of quads.
 * \param na pointer to the number of edges.
 * \return 1.
 *
 * \brief Get the number of vertices, triangles and edges of the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_MESHSIZE(mesh,np,nt,nquad,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER(MMG5F_INT)            :: np,nt,nquad,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_meshSize(MMG5_pMesh mesh, MMG5_int* np, MMG5_int* nt, MMG5_int* nquad, MMG5_int* na);

/**
 * \brief Get the number of solutions, their dimension and their type.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \param typEntity pointer to the type of entities to which solutions are applied
 *        (see \ref MMG5_entities for possible values).
 * \param np pointer to the number of solutions.
 * \param typSol pointer to the type of the solutions
 *               (scalar, vectorial, ..., see \ref MMG5_type for possible values)
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_SOLSIZE(mesh,sol,typEntity,np,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER                       :: typEntity,typSol\n
 * >     INTEGER(MMG5F_INT)            :: np\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG2D_EXPORT int  MMG2D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, MMG5_int* np,
                                        int* typSol);

/**
 * \brief Get the number of elements and dimension of a solution defined on vertices.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an array of sol structure.
 * \param nentities pointer to the number of entities.
 * \param typSol array of size MMG5_NSOL_MAX to store type of each solution
 *        (scalar, vectorial, ..., see \ref MMG5_type for possible values).
 *
 * \return 1.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_SOLSATVERTICESSIZE(mesh,sol,nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER                       :: nsols\n
 * >     INTEGER(MMG5F_INT)            :: nentities\n
 * >     INTEGER                       :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol* sol,int *nsols,
                                                    MMG5_int* nentities,int* typSol);

/**
 * \brief Get the coordinates and reference ref of the next vertex of a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 pointer to the coordinate of the vertex along the first dimension.
 * \param c1 pointer to the coordinate of the vertex along the second dimension.
 * \param ref pointer to the vertex reference.
 * \param isCorner pointer to the flag saying if vertex is corner.
 * \param isRequired pointer to the flag saying if vertex is required.
 * \return 1.
 *
 * This function retrieves the coordinates \a c0 and \a c1, and reference \a
 * ref of the next vertex of a mesh. It is meant to be used in a loop over all
 * vertices. When this function has been called as many times as there are
 * vertices, the internal loop counter will be reset. To obtain data for a
 * specific vertex, the \ref MMG2D_GetByIdx_vertex function can be used instead.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_VERTEX(mesh,c0,c1,ref,isCorner,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1\n
 * >     INTEGER(MMG5F_INT)            :: ref\n
 * >     INTEGER                       :: isCorner,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, MMG5_int* ref,
                                        int* isCorner, int* isRequired);

/**
 * \brief Get the coordinates and reference of a specific vertex in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param c0 pointer to the coordinate of the vertex along the first dimension.
 * \param c1 pointer to the coordinate of the vertex along the second dimension.
 * \param ref pointer to the vertex reference.
 * \param isCorner pointer to the flag saying if vertex is corner.
 * \param isRequired pointer to the flag saying if vertex is required.
 * \param idx index of vertex to get.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1 and reference \a ref of
 * vertex \a idx of mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GETBYIDX_VERTEX(mesh,c0,c1,ref,isCorner,isRequired,idx,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1\n
 * >     INTEGER                       :: isCorner,isRequired\n
 * >     INTEGER(MMG5F_INT)            :: ref,idx\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_GetByIdx_vertex(MMG5_pMesh mesh, double* c0, double* c1, MMG5_int* ref,
                                             int* isCorner, int* isRequired,MMG5_int idx);

/**
 * \brief Get the coordinates and references of all vertices in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param vertices pointer to the array of vertex coordinates.
 * The coordinates of the \f$i^{th}\f$ vertex are stored in
 * vertices[(i-1)*2]\@2.
 * \param refs pointer to the array of references.
 * The ref of the \f$i^th\f$ vertex is stored in refs[i-1].
 * \param areCorners pointer to the array of flags saying if
 * vertices are corners.
 * areCorners[i-1]=1 if the \f$i^{th}\f$ vertex is corner.
 * \param areRequired pointer to the array of flags saying if vertices
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ vertex is required.
 * \return 1.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners and areRequired arrays)
 * > !  SUBROUTINE MMG2D_GET_VERTICES(mesh,vertices,refs,areCorners,&\n
 * > !                                areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)          :: mesh\n
 * > !    REAL(KIND=8),DIMENSION(*), INTENT(OUT) :: vertices\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)       :: refs\n
 * > !    INTEGER, DIMENSION(*)                  :: areCorners,areRequired\n
 * > !    INTEGER, INTENT(OUT)                   :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_vertices(MMG5_pMesh mesh, double* vertices, MMG5_int* refs,
                                          int* areCorners, int* areRequired);

/**
 * \brief Get the vertices and reference of the next triangle in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the first vertex of triangle.
 * \param v1 pointer to the second vertex of triangle.
 * \param v2 pointer to the third vertex of triangle.
 * \param ref pointer to the triangle reference.
 * \param isRequired pointer to the flag saying if triangle is required.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the vertices \a v0, \a v1, \a v2, and reference \a
 * ref of the next triangle of \a mesh. It is meant to be called in a loop over
 * all triangles. When it has been called as many times as there are triangles,
 * the internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_TRIANGLE(mesh,v0,v1,v2,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: v0,v1,v2\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER                        :: isRequired\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_triangle(MMG5_pMesh mesh, MMG5_int* v0,
                                         MMG5_int* v1, MMG5_int* v2, MMG5_int* ref
                                         ,int* isRequired);

/**
 * \brief Get the vertices and references of all triangles in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param tria pointer to the array of triangles vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer to the array of triangles references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \param areRequired pointer to array of flags saying if triangles
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ tria
 * is required.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs and areRequired arrays)
 * > !  SUBROUTINE MMG2D_GET_TRIANGLES(mesh,tria,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*),INTENT(OUT) :: tria\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)             :: refs\n
 * > !    INTEGER, DIMENSION(*)                        :: areRequired\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_triangles(MMG5_pMesh mesh, MMG5_int* tria, MMG5_int* refs,
                                           int* areRequired);

/**
 * \brief Get the vertices and reference of the next quadrangle of the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param v0 pointer to the first vertex of quadrangle.
 * \param v1 pointer to the second vertex of quadrangle.
 * \param v2 pointer to the third vertex of quadrangle.
 * \param v3 pointer to the fourth vertex of quadrangle.
 * \param ref pointer to the quadrangle reference.
 * \param isRequired pointer to the flag saying if quadrangle is
 *  required.
 * \return 0 on failure, 1 otherwise.
 *
 * Get the vertices \a v0,\a v1,\a v2,\a v3 and reference \a ref of the next
 * quadrangle of mesh. This function is meant to be called in a loop over all
 * quadrangles. When it has been called as many times as there are
 * quadrangles, the internal loop counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_QUADRILATERAL(mesh,v0,v1,v2,v3,ref,isRequired,&\n
 * >                                      retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: v0,v1,v2,v3\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER                        :: isRequired\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_quadrilateral(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2,
                                               MMG5_int* v3,MMG5_int* ref, int* isRequired);

/**
 * \brief Get the vertices and references of all quadrangles of the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param quadra pointer to the array of quadrangles vertices.
 * Vertices of the \f$i^{th}\f$ quadrangle are stored in quadra[(i-1)*4]\@4.
 * \param refs pointer to the array of quadrlaterals references.
 * References of the \f$i^{th}\f$ quad is stored in refs[i-1].
 * \param areRequired pointer to the array of flags saying if the
 *  quadrangles are required. areRequired[i-1]=1 if the \f$i^{th}\f$ quad
 * is required.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > !  SUBROUTINE MMG2D_GET_QUADRILATERALS(mesh,quadra,refs,areRequired,&\n
 * > !                                      retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)                :: mesh\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*),INTENT(OUT) :: quadra\n
 * > !    INTEGER(MMG5F_INT), DIMENSION(*)             :: refs\n
 * > !    INTEGER, DIMENSION(*)                        :: areRequired\n
 * > !    INTEGER, INTENT(OUT)                         :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
 LIBMMG2D_EXPORT int  MMG2D_Get_quadrilaterals(MMG5_pMesh mesh, MMG5_int* quadra,MMG5_int* refs,
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
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves the extremities \a e0, \a e1 and reference \a ref of
 * next edge of \a mesh. It is meant to be called in a loop over all edges. When
 * it has been called as many times as there are edges in the mesh, the internal
 * edge counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_EDGE(mesh,e0,e1,ref,isRidge,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: e0,e1\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER                        :: isRidge,isRequired\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_edge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref,
                                     int* isRidge, int* isRequired);

/**
 * \brief Get the vertices and references of all edges in a mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param edges pointer to the array of edges.
 * Vertices of the \f$i^{th}\f$ edge are stored in edge[(i-1)*2]\@2.
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \param areRidges 1 if the edge is a ridge, 0 otherwise.
 * \param areRequired 1 if the edge is required, 0 otherwise.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areRidges or areRequired arrays)
 *
 * \remark Fortran interface:
 * > !  SUBROUTINE MMG2D_GET_EDGES(mesh,edges,refs,areRidges,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * > !    INTEGER(MMG5F_INT), INTENT(IN) :: edges(*)\n
 * > !    INTEGER(MMG5F_INT), INTENT(OUT):: refs(*)\n
 * > !    INTEGER, INTENT(OUT)           :: areRequired(*),areRidges(*)\n
 * > !    INTEGER, INTENT(OUT)           :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_edges(MMG5_pMesh mesh,MMG5_int *edges,MMG5_int* refs,
                                      int *areRidges,int *areRequired);

/**
 * \brief Get the quality measure of a single triangle in the mesh.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of the triangle for which we want to get the quality.
 * \return the computed quality, or 0 on failure
 *
 * Get quality of triangle \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_TRIANGLEQUALITY(mesh,met,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,met\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8), INTENT(OUT)     :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT double MMG2D_Get_triangleQuality(MMG5_pMesh mesh,MMG5_pSol met, MMG5_int k);

/**
 * \param met pointer to the sol structure.
 * \param s pointer to the scalar solution value.
 * \return 0 on failure, 1 otherwise.
 *
 * \brief Get the scalar solution value \a s of next element of a solution.
 *
 * This function is meant to be called in a loop over all elements. When it has
 * been called as many times as there are elements in the solution, the internal
 * counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_SCALARSOL(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_scalarSol(MMG5_pSol met, double* s);

/**
 * \brief Get all elements of a scalar sol structure.
 *
 * \param met pointer to the sol structure.
 * \param s array of scalar solutions at mesh vertices. s[i-1] is
 * the solution at vertex i.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_SCALARSOLS(met,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_scalarSols(MMG5_pSol met, double* s);

/**
 * \brief Get the next element of a vector sol structure.
 *
 * \param met pointer to the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \return 0 on failure, 1 otherwise.
 *
 * This function retrieves vectorial solution \f$(v_x,v_y)\f$ of the next
 * element of \a met.  It is meant to be called in a loop over all
 * elements. When it has been called as many times as there are elements in the
 * solution, the internal counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_VECTORSOL(met,vx,vy,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: vx,vy\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy);

/**
 * \brief Get all elements of a vector sol structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array to store the data in the order \f$[x_1, y_1, x_2, \ldots, y_N]\f$
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_VECTORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: met\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_vectorSols(MMG5_pSol met, double* sols);

/**
 * \brief Get the next element of a tensor sol structure.
 *
 * \param met pointer to the sol structure.
 * \param m11 pointer to the position (1,1) in the solution tensor.
 * \param m12 pointer to the position (1,2) in the solution tensor.
 * \param m22 pointer to the position (2,2) in the solution tensor.
 * \return 0 on failure, 1 otherwise.
 *
 * This function is meant to be called in a loop over all elements. When it has
 * been called as many times as there are elements in the solution, the internal
 * counter will be reset.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_TENSORSOL(met,m11,m12,m22,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: m11,m12,m22\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12,double *m22);

/**
 * \brief Get all elements of a tensor sol structure.
 *
 * \param met pointer to the sol structure.
 * \param sols array of solutions at mesh vertices.
 * sols[3*(i-1)]\@3 is the solution at vertex i.
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)           :: met\n
 * >     REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)                    :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_tensorSols(MMG5_pSol met, double *sols);

/**
 * \brief Get one out of several scalar solutions at a specific vertex.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s solution(s) at mesh vertex \a pos.
 * \param pos index of the vertex on which we get the solution.
 *
 * \return 0 on failure, 1 otherwise.
 *
 * Get values of the ith field of the solution array at vertex \a pos.  (pos
 * from 1 to the number of vertices included and \a i from 1 to the number of
 * solutions).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_ITHSOL_INSOLSATVERTICES(sol,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)         :: pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);

/**
 * \brief Get one out of several scalar solutions at all vertices in the mesh.
 *
 * \param sol pointer to the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s array of solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[2*(k-1)]\@2 for a vectorial solution
 * and s[3*(k-1)]\@3 for a tensor solution.
 *
 * \return 0 on failure, 1 otherwise.
 *
 * Get values of the solution at the ith field of the solution array.
 * (\a i from 1 to \a nb_sols)
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_ITHSOLS_INSOLSATVERTICES(sol,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: sol\n
 * >     INTEGER, INTENT(IN)                    :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_Get_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);


/**
 * \brief Check if the number of given entities match with mesh and sol size
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \return 0 on failure, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_Chk_meshData(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Deallocate an array of solution fields
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to an array of solution structure (that stores solution fields).
 * \return 1
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_Free_allSols(mesh,sol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Free_allSols(MMG5_pMesh mesh,MMG5_pSol *sol);

/* deallocations */
/**
 * \brief Deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend to the library function that you
 * have call.
 *
 * For the \ref MMG2D_mmg2dlib function, you need
 * to call the \ref MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the \ref MMG2D_mmg2dls function, you need
 * to call the \ref MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the \ref MMG2D_mmg2dmov function, you must call
 * : MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 * MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * \return 0 on failure, 1 on success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
LIBMMG2D_EXPORT int MMG2D_Free_all(const int starter,...);

/**
 * \brief Structure deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend to the library function that you
 * have call.
 *
 * For the \ref MMG2D_mmg2dlib function, you need
 * to call the \ref MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the \ref MMG2D_mmg2dls function, you need
 * to call the \ref MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the \ref MMG2D_mmg2dmov function, you must call
 * : MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * \return 0 on failure, 1 on success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
 LIBMMG2D_EXPORT int MMG2D_Free_structures(const int starter,...);

/**
 * \brief Structure deallocations before return.
 *
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend to the library function that you
 * have call.
 *
 * For the MMG2D_mmg2dlib function, you need
 * to call the \a MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&your_metric,MMG5_ARG_end).
 *
 * For the MMG2D_mmg2dls function, you need
 * to call the \a MMG2D_Init_mesh function with the following arguments :
 * MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * For the MMG2D_mmg2dmov function, you must call
 * : MMG2D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh,
 *  MMG5_ARG_ppMet,&empty_metric,MMG5_ARG_ppDisp, &your_displacement,
 * MMG5_ARG_end).
 *
 * \return 0 on failure, 1 otherwise
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Free_names(const int starter,...);

/**
 * \brief Load a mesh (in .mesh/.mesb format) from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the readed file.
 *
 * \return 0 if the file is not found, -1 in case of failure for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadMesh(MMG5_pMesh mesh,const char * filename);

/**
 * \brief Load a mesh and possibly a solution in VTP (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * Read a mesh and 0 or 1 data fields in VTK vtp file format (.vtp extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADVTPMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadVtpMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a mesh and multiple solutions in VTP (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * Read a mesh and a list of data fields in VTK vtp file format (.vtp extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADVTPMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadVtpMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Load a mesh and possibly data in VTU (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * Read a mesh and 0 or 1 data fields in VTK vtu file format (.vtu extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADVTUMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG2D_EXPORT int MMG2D_loadVtuMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a mesh and multiple solutions in VTU (VTK) format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * Read a mesh and a list of data fields in VTK vtu file format (.vtu extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADVTUMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG2D_EXPORT int MMG2D_loadVtuMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Load a mesh and possibly data in VTK format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * Read mesh and 0 or 1 data fields in VTK file format (.vtk extension). We
 * read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADVTKMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG2D_EXPORT int MMG2D_loadVtkMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a mesh and multiple solutions in VTK format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * This function reads a mesh and a list of data fields in VTK file format (.vtk
 * extension). We read only low-order vertices, edges, triangles and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADVTKMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 LIBMMG2D_EXPORT int MMG2D_loadVtkMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Load a mesh and possibly a solution in .msh format from file.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (lack of
 * memory, file format...), 1 on success.
 *
 * This function reads a mesh and 0 or 1 data fields in MSH file format (.msh
 * extension). We read only low-order vertices, edges, triangles, and quadrangles.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADMSHMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Load a mesh and all data from a file in MSH format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to a list of solution structures.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (lack of
 * memory, file format...), 1 on success.
 *
 * This function reads a mesh and all data fields from a file in MSH file format
 * (.msh extension). We read only low-order vertices, edges, triangles,
 * quadrangles, tetrahedra and prisms.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADMSHMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

  /* FIXME: why is it called medit format and is this really specific for metrics? */
/**
 * \brief Load a metric field (or other solution) in medit's .sol format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure..
 * \param filename name of the solution file.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * This function loads a metric field. The file in medit format must contain 1
 * solution: the metric.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADSOL(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadSol(MMG5_pMesh mesh,MMG5_pSol sol,const char * filename);

/**
 * \brief Read mesh data in a format determined by the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure or the NULL pointer.
 * \param sol pointer to the level-set structure or the NULL pointer.
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason (insufficient memory, file
 * format...), 1 on success.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADGENERICMESH(mesh,met,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadGenericMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,const char *filename);

/**
 * \brief Load one or more solutions in a solution file in medit file format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solutions array
 * \param filename name of the file to load.
 *
 * \return 0 if the file is not found, -1 if failing for another reason
 * (insufficient memory, file format...), 1 on success.
 *
 * Load 1 or more solutions in a solution file in medit file format
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_LOADALLSOLS(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_loadAllSols(MMG5_pMesh mesh,MMG5_pSol *sol, const char* filename);

  /* FIXME: why is this here, neither implemented nor documented? */
  LIBMMG2D_EXPORT int MMG2D_loadVect(MMG5_pMesh ,char *);

/**
 * \brief Save a mesh in .mesh/.meshb format.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the readed file.
 * \return 0 or -1 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveMesh(MMG5_pMesh ,const char *);

/**
 * \brief Save a mesh and optionally one data field in MSH format, ascii or
 * binary depending on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and optionally one data field in MSH file format
 * (.msh extension). It uses ASCII format for .msh extension, binary format for
 * .msb extension.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEMSHMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Save a mesh and multiple data fields in MSH format, ascii or binary
 * depending on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and a list of data fields in MSH  file format (.msh extension).
 * It uses ASCII format for .msh extension, binary format for .mshb extension.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEMSHMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Save a mesh and optionally one solution in VTK format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and 0 or 1 data fields in Vtk file format (.vtk extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEVTKMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveVtkMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

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
 * >   SUBROUTINE MMG2D_SAVEVTKMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveVtkMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Save a mesh and optionally one data field in VTU format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and 0 or 1 data fields in vtu Vtk file format (.vtu extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEVTUMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveVtuMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Save a mesh and multiple data fields in VTU format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and a list of data fields in vtu Vtk file format (.vtu extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEVTUMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveVtuMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Save a mesh and optionally one data field in VTP format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and optionally one data field in polydata Vtk
 * file format (.vtp extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEVTPMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveVtpMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Save a mesh and multiple data fields in VTP format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the file to write.
 * \return 0 on failure, 1 otherwise.
 *
 * This function writes a mesh and a list of data fields in polydata Vtk file
 * format (.vtp extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEVTPMESH_AND_ALLDATA(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveVtpMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \brief Save data in Tetgen's Triangle format.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the file to write
 * \return 0 or -1 on failure, 1 otherwise.
 *
 * This function saves mesh data in Triangle (or equivalent to Tetgen in 2D)
 * file format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVETETGENMESH(mesh,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveTetgenMesh(MMG5_pMesh ,const char *);

/**
 * \brief Save mesh data in a file whose format depends on the filename extension.
 *
 * \param mesh pointer to the mesh structure.
 * \param filename name of the file to write
 * \return 0 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEGENERICMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveGenericMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \brief Save metric field in medit solution file format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solution structure.
 * \param filename name of the solution file to write.
 * \return 0 or -1 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVESOL(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveSol(MMG5_pMesh  mesh,MMG5_pSol sol ,const char *filename);

/**
 * \brief Save one or more solutions in a solution file in medit file format.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the solutions array
 * \param filename name of the solution file.
 * \return 0 or -1 on failure, 1 otherwise.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SAVEALLSOLS(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_saveAllSols(MMG5_pMesh  mesh,MMG5_pSol *sol ,const char *filename);

  /* FIXME: why is this here? */
  LIBMMG2D_EXPORT int MMG2D_saveVect(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename,double lambda);

/**
 * \brief Main "program" for the mesh adaptation library.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to a sol structure (metric).
 * \return \ref MMG5_SUCCESS if successful, \ref MMG5_LOWFAILURE in case there is a  failure
 * but a conform mesh is returned and \ref MMG5_STRONGFAILURE if there is a failure and we
 * can't save the mesh.
 *
 * This function adapts a given mesh, trying to improve the quality, under the
 * given metric and parameters.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_MMG2DLIB(mesh,sol,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_mmg2dlib(MMG5_pMesh mesh,MMG5_pSol sol);

/**
 * \brief Main "program" for the mesh generation library.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to a sol structure (metric).
 * \return \ref MMG5_SUCCESS if successful, \ref MMG5_LOWFAILURE if there is a failure
 * but a conform mesh is returned and \ref MMG5_STRONGFAILURE if there is a failure and we
 * can't save the mesh.
 *
 * FIXME: This function creates a triangular mesh from a given polygon, right?
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_MMG2DMESH(mesh,sol,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_mmg2dmesh(MMG5_pMesh mesh,MMG5_pSol sol);

/**
 * \brief Main "program" for the level-set discretization library.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to a sol structure (level-set function).
 * \param met pointer to a sol structure (metric).
 * \return \ref MMG5_SUCCESS if successful, \ref MMG5_LOWFAILURE if there is a failure
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if there is a failure and we
 * can't save the mesh.
 *
 * This is the main program for the level-set discretization library. If a
 * metric \a met is provided, it is used to adapt the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_MMG2DLS(mesh,sol,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     MMG5_DATA_PTR_T                :: met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_mmg2dls(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) ;

/**
 * \brief Main "program" for the rigid-body movement library.
 *
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to a sol structure (displacement).
 * \param disp pointer to a sol (displacement for the lagrangian motion
 * mode) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if there is a failure
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if there is a failure and we
 * can't save the mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_MMG2DMOV(mesh,sol,disp,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol,disp\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_mmg2dmov(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol disp);

/* Tools for the library */

/**
 * \brief Print the default parameters values.
 *
 * \param mesh pointer to the mesh structure.
 * \return 0 on failure, 1 on success.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_DEFAULTVALUES(mesh,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_defaultValues(MMG5_pMesh mesh);

/**
 * \brief Store command line arguments.
 *
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer to the mesh structure.
 * \param met pointer to a metric
 * \param sol pointer to a level-set or displacement function
 * \return 1 if we want to run Mmg after, 0 if not or in case of failure
 *
 * \remark no matching fortran function.
 *
 */
  LIBMMG2D_EXPORT int MMG2D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol);

/**
 * \brief Read a file containing Local parameters (.mmg2d extension)
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \return 1.
 *
 * This function reads a local parameters file. This file must have the same
 * name as the mesh with the \a .mmg2d extension or must be named \a
 * DEFAULT.mmg2d.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_PARSOP(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_parsop(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Print help for mmg2d options.
 *
 * \param prog pointer to the program name.
 * \param return 1 on success, 0 on failure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_USAGE(prog,strlen0,retval)\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: prog\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int  MMG2D_usage(char *prog);

/**
 * \brief Compute unit tensor according to the lengths of the
 * edges passing through a vertex.
 *
 * \param mesh pointer to the mesh structure
 * \param met pointer to the sol structure
 * \return 1 on success
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_DOSOL(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT extern int (*MMG2D_doSol)(MMG5_pMesh mesh ,MMG5_pSol met );

/**
 * \brief Compute a constant size map according to the hsiz, hmin and hmax parameters.
 *
 * \param mesh pointer to the mesh structure
 * \param met pointer to the sol structure
 * \return 1 on success
 *
 * This function computes a constant size map according to mesh->info.hsiz,
 * mesh->info.hmin and mesh->info.hmax. It updates these 3 values if not
 * compatible.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SET_CONSTANTSIZE(mesh,met,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh,met\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \brief Set function pointers for length, caltri... depending if case is iso or aniso
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to a sol structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SETFUNC(mesh,met)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN) :: mesh,met\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT void MMG2D_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

  /* FIXME: is this description correct? */
/**
 * \brief Get the number of non-boundary edges.
 *
 * \param mesh pointer to the mesh structure.
 * \param nb_edges pointer to the number of non boundary edges.
 * \return 0 on failure, 1 otherwise.
 *
 * This function extracts the number of non boundary edges (for DG methods for
 * example). An edge is boundary if it is located at the interface of two
 * domains with different references, if it belongs to one triangle only or if
 * it is a singular edge (ridge or required).
 *
 * Append these edges to the list of edges.
 *
 * \warning reallocate the edge array and append the internal edges. This may
 * modify the behaviour of other functions.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_NUMBEROFNONBDYEDGES(mesh,nb_edges,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: nb_edges\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_numberOfNonBdyEdges(MMG5_pMesh mesh, MMG5_int* nb_edges);

/**
 * \brief Get vertices and reference of a non-boundary edge.
 *
 * \param mesh pointer to the mesh structure.
 * \param e0 pointer to the first extremity of the edge.
 * \param e1 pointer to the second  extremity of the edge.
 * \param ref pointer to the edge reference.
 * \param idx index of the non boundary edge to get (between 1 and nb_edges)
 * \return 0 on failure, 1 otherwise.
 *
 * This function returns the extremities \a e0, \a e1 and reference \a ref of
 * the idx^th non boundary edge (for DG methods for example). An edge is
 * boundary if it is located at the interface of 2 domains with different
 * references, if it belongs to one triangle only or if it is a singular edge
 * (ridge or required).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_NONBDYEDGE(mesh,e0,e1,ref,idx,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: e0,e1\n
 * >     INTEGER(MMG5F_INT)             :: ref\n
 * >     INTEGER(MMG5F_INT), INTENT(IN) :: idx\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_nonBdyEdge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref, MMG5_int idx);

/**
 * \brief Return adjacent elements of a triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param kel triangle index.
 * \param listri pointer to the array of indices of the three adjacent
 * triangles of the elt \a kel (the index is 0 if there is no adjacent).
 * \return 1.
 *
 * Find the indices of the 3 adjacent elements of triangle \a
 * kel. \f$v_i = 0\f$ if the \f$i^{th}\f$ face has no adjacent element
 * (so we are on a boundary face).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_ADJATRI(mesh,kel,listri,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: mesh\n
 * >     INTEGER, INTENT(IN)                :: kel\n
 * >     INTEGER, DIMENSION(3), INTENT(OUT) :: listri\n
 * >     INTEGER, INTENT(OUT)               :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_adjaTri(MMG5_pMesh mesh, MMG5_int kel, MMG5_int listri[3]);

/**
 * \brief Return adjacent vertices of a triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param ip vertex index.
 * \param lispoi pointer to an array of size MMG2D_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a ip.
 * \return nbpoi the number of adjacent vertices if success, 0 on failure.
 *
 * Find the indices of the adjacent vertices of the vertex \a
 * ip.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_ADJAVERTICES(mesh,ip,lispoi,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)                         :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)                         :: ip\n
 * >     INTEGER(MMG5F_INT), DIMENSION(MMG2D_LMAX), INTENT(OUT) :: lispoi\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT MMG5_int MMG2D_Get_adjaVertices(MMG5_pMesh mesh, MMG5_int ip, MMG5_int lispoi[MMG2D_LMAX]);

/**
 * \brief Return adjacent vertices of a triangle.
 *
 * \param mesh pointer to the mesh structure.
 * \param ip vertex index.
 * \param start index of a triangle holding \a ip.
 * \param lispoi pointer to an array of size MMG2D_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a ip.
 * \return nbpoi the number of adjacent vertices if success, 0 on failure.
 *
 * Find the indices of the adjacent vertices of the vertex \a
 * ip of the triangle \a start.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_ADJAVERTICESFAST(mesh,ip,start,lispoi,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)                         :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)                         :: ip,start\n
 * >     INTEGER(MMG5F_INT), DIMENSION(MMG2D_LMAX), INTENT(OUT) :: lispoi\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT)                        :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT MMG5_int MMG2D_Get_adjaVerticesFast(MMG5_pMesh mesh, MMG5_int ip,MMG5_int start,
                                                      MMG5_int lispoi[MMG2D_LMAX]);

/**
 * \brief Find a triangle given an adjacent triangle and an edge number.
 *
 * \param mesh pointer to the mesh structure.
 * \param ked index of the boundary edge.
 * \param ktri pointer to the index of the tri (filled by the function).
 * \param ied pointer to the index of the edge of the triangle \a ktri that
 * correspond to the boundary edge \a ked.
 * \return 0 on failure, 1 otherwise
 *
 * Fill \a ktri by the index of the triangle to which belong a boundary edge
 * and \a ied by the index of the edge of the triangle that correspond to the
 * edge.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_TRIFROMEDGE(mesh,ked,ktri,ied,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)        :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)     :: ked\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT)    :: ktri\n
 * >     INTEGER, INTENT(OUT)               :: retval,ied\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_triFromEdge(MMG5_pMesh mesh, MMG5_int ked, MMG5_int *ktri, int *ied);

/**
 * \brief Find two triangles given the edge that they share.
 *
 * \param mesh pointer to the mesh structure.
 * \param ked index of the boundary edge.
 * \param ktri pointer to an array of size 2 to fill by the indices of the
 * triangles that share the edge \a ked (filled by the function).
 * \param ied pointer to an array of size two to fill by the indices of the
 * edge in each triangle.
 *
 * \return 0 on failure, 1 otherwise
 *
 * Fill \a ktri by the indices of the triangles to which belong a boundary edge
 * and \a ied by the indices of the matching edge in each triangle. If \a ked
 * belongs to one triangle only, ktri[1] = ied[1] = 0.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_GET_TRISFROMEDGE(mesh,ked,ktri,ied,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)                  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(IN)               :: ked\n
 * >     INTEGER(MMG5F_INT), DIMENSION(2),INTENT(OUT) :: ktri\n
 * >     INTEGER, INTENT(OUT)                         :: retval,ied\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Get_trisFromEdge(MMG5_pMesh mesh, MMG5_int ked, MMG5_int ktri[2],int ied[2]);

/**
 * \brief Compute the real eigenvalues and eigenvectors of a symmetric matrix
 *
 * \param m upper part of a symMetric matrix diagonalizable in |R
 * \param lambda array of the metric eigenvalues
 * \param vp array of the metric eigenvectors
 *
 * \return the order of the eigenvalues
 *
 * This function computes the real eigenvalues and eigenvectors of a symmetric matrix m
 * whose upper part is provided (m11, m12, m22, in this order).
 *
 * lambda[0] is the eigenvalue associated to the eigenvector ( v[0][0], v[0,1] )
 * in C and to the eigenvector v(1,:) in fortran
 *
 * lambda[1] is the eigenvalue associated to the eigenvector ( v[1][0], v[1,1] )
 * in C and to the eigenvector v(2,:) in fortran
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_COMPUTE_EIGENV(m,lambda,vp,retval)\n
 * >     REAL(KIND=8), INTENT(IN)         :: m(*)\n
 * >     REAL(KIND=8), INTENT(OUT)        :: lambda(*),vp(*)\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_Compute_eigenv(double m[3],double lambda[2],double vp[2][2]);

/**
 * \brief Reset the vertex tags.
 *
 * \param mesh pointer to the mesh structure
 *
 * This function resets the tags of all vertices. Be careful: all the tags are deleted.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_RESET_VERTICESTAGS(mesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT void MMG2D_Reset_verticestags(MMG5_pMesh mesh);

/**
 * \brief Free the mesh elements (and the adjacency information).
 *
 * \param mesh pointer to the mesh structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_FREE_TRIANGLES(mesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT void MMG2D_Free_triangles(MMG5_pMesh mesh);

/**
 * \brief Free the mesh edges (and the associated xpoints).
 *
 * \param mesh pointer to the mesh structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_FREE_EDGES(mesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT void MMG2D_Free_edges(MMG5_pMesh mesh);

/**
 * \brief Free the solution.
 *
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the solution structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_FREE_SOLUTIONS(mesh,sol)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT void MMG2D_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol);


/**
 * \brief Set common function pointers between mmgs and mmg2d to the matching mmg2d
 * functions.
 */
  LIBMMG2D_EXPORT void MMG2D_Set_commonFunc(void);

/**
 * \brief Normalize the mesh and size information.
 *
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param ls pointer to a solution structure (level-set or displacement).
 *
 * \return 1 on success, 0 in case of failure (computed bounding box too small
 * or one af the anisotropic input metric is not valid).
 *
 * This function scales the mesh and the size information between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG2D_SCALEMESH(mesh,met,ls,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,ls\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG2D_EXPORT int MMG2D_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol ls);

#ifdef __cplusplus
}
#endif

#endif
