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
 * \file mmgs/libmmgs.h
 * \brief API headers for the mmgs library
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning To keep the genheader working, don't break line between the enum
 * name and the opening brace (it creates errors under windows)
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
 * Maximum array size when storing adjacent points (or ball) of a vertex.
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
  MMGS_IPARAM_verbose,           /*!< [-1..10], Tune level of verbosity */
  MMGS_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  MMGS_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  MMGS_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  MMGS_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  MMGS_IPARAM_isosurf,           /*!< [1/0], Level-set meshing on the surface part */
  MMGS_IPARAM_isoref,            /*!< [0/n], Iso-surface boundary material reference */
  MMGS_IPARAM_keepRef,           /*!< [1/0], Preserve the initial domain references in level-set mode */
  MMGS_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  MMGS_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  MMGS_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  MMGS_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  MMGS_IPARAM_nreg,              /*!< [0/1], Disabled/enabled normal regularization */
  MMGS_IPARAM_xreg,              /*!< [0/1], Disabled/enabled coordinates regularization */
  MMGS_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  MMGS_IPARAM_numberOfLSBaseReferences, /*!< [n], Number of base references for bubble removal */
  MMGS_IPARAM_numberOfMat,              /*!< [n], Number of material in ls mode */
  MMGS_IPARAM_numsubdomain,      /*!< [0/n], Save the subdomain nb (0==all subdomain) */
  MMGS_IPARAM_renum,             /*!< [1/0], Turn on/off point relocation with Scotch */
  MMGS_IPARAM_anisosize,         /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  MMGS_IPARAM_nosizreq,          /*!< [0/1], Allow/avoid overwritten of sizes at required points (advanced usage) */
  MMGS_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  MMGS_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  MMGS_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  MMGS_DPARAM_hsiz,              /*!< [val], Constant mesh size */
  MMGS_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  MMGS_DPARAM_hgrad,             /*!< [val], Control gradation */
  MMGS_DPARAM_hgradreq,          /*!< [val], Control gradation on required entites (advanced usage) */
  MMGS_DPARAM_ls,                /*!< [val], Value of level-set */
  MMGS_DPARAM_rmc,               /*!< [-1/val], Remove small connex componants in level-set mode */
  MMGS_PARAM_size,               /*!< [n], Number of parameters */
};

/*----------------------------- functions header -----------------------------*/
/* Initialization functions */
/* init structures */
/**
 * \param starter dummy argument used to initialize the variadic argument list
 * \param ... variadic arguments.
 *
 * For the MMGS_mmgslib function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the MMGS_mmgsls function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here,\a your_mesh is a \a MMG5_pMesh, \a your_metric and \a your_level_set
 * are \a MMG5_pSol.
 *
 * \return 1 if success, 0 if fail
 *
 * MMG structures allocation and initialization.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 */
LIBMMGS_EXPORT int MMGS_Init_mesh(const int starter,...);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialize file names to their default values.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_INIT_FILENAMES(mesh,sol)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void  MMGS_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
/**
 * \param mesh pointer toward the mesh structure.
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
 * \param mesh pointer toward the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
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
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
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

/* init structure sizes */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity type of solutions entities (vertices, triangles...).
 * \param np number of solutions.
 * \param typSol type of solution (scalar, vectorial...).
 * \return 0 if failed, 1 otherwise.
 *
 * Initialize an array of solutions field: set dimension, types and number of
 * data.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward an allocatable sol structure.
 * \param nsols number of solutions per entity
 * \param nentities number of entities
 * \param typSol    Array of size nsol listing the type of the solutions
 *                  (scalar, vectorial...).
 * \return 0 if failed, 1 otherwise.
 *
 * Initialize an array of solutions field defined at vertices: set dimension,
 * types and number of data.
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
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param nt number of triangles.
 * \param na number of edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices, triangles and edges of the
 * mesh and allocate the associated tables. If call twice, reset the
 * whole mesh to realloc it at the new size
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
 * at position \a pos in mesh structure (\a pos from 1 to nb_vertices included).
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
                     double c2, MMG5_int ref,MMG5_int pos);
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
 * \param mesh pointer toward the mesh structure.
 * \param v0 first vertex of triangle.
 * \param v1 second vertex of triangle.
 * \param v2 third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure.(\a pos from 1 to nb_tria included).
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
 * \param mesh pointer toward the mesh structure.
 * \param tria pointer toward the table of the tria vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer toward the table of the triangle references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \return 0 if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh triangles.
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
 * \param mesh pointer toward the mesh structure.
 * \param v0 first extremity of the edge.
 * \param v1 second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set edges of extremities \a v0, \a v1 and reference \a ref at
 * position \a pos in mesh structure (\a pos from 1 to nb_edges included).
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
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set corner at point \a pos (\a pos from 1 to nb_vertices included).
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
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Remove corner attribute at point \a pos (from 1 to nb_vertices included).
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
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Set point \a k as required.
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
 * \param mesh pointer toward the mesh structure.
 * \param k vertex index.
 * \return 1.
 *
 * Remove required attribute from point \a k.
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
 * \param mesh pointer toward the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Set triangle \a k as required.
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
 * \param mesh pointer toward the mesh structure.
 * \param k triangle index.
 * \return 1.
 *
 * Remove required attribute from triangle \a k.
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
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set ridge at edge \a k.
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
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Remove ridge attribute at edge \a k.
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
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Set edge \a k as required.
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
 * \param mesh pointer toward the mesh structure.
 * \param k edge index.
 * \return 1.
 *
 * Remove required attribute from edge \a k.
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
 * \param mesh pointer toward the mesh structure.
 * \param edges pointer toward the array of edges.
 * Vertices of the \f$i^{th}\f$ edge are stored in edge[(i-1)*2]\@2.
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \return 0 if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh edges.
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
 * \param mesh pointer toward the mesh structure.
 * \param edges pointer toward the array of edges.
 * Vertices of the \f$i^{th}\f$ edge are stored in edge[(i-1)*2]\@2.
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \param areRidges 1 if the edge is a ridge, 0 otherwise.
 * \param areRequired 1 if the edge is required, 0 otherwise.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh edges.
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
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param k index of the triangle for which we want to get the quality.
 * \return the computed quality or 0. if fail.
 *
 * Get quality of tria \a k.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_GET_TRIANGLEQUALITY(mesh,met,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,met\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: k\n
 * >     REAL(KIND=8), INTENT(OUT)     :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  double MMGS_Get_triangleQuality(MMG5_pMesh mesh,MMG5_pSol met, MMG5_int k);

/**
 * \param met pointer toward the sol structure.
 * \param s solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure
 * (\a pos from 1 to nb_vertices included).
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
LIBMMGS_EXPORT int  MMGS_Set_scalarSol(MMG5_pSol met, double s,MMG5_int pos);
/**
 * \param met pointer toward the sol structure.
 * \param s table of the scalar solutions values.
 * s[i-1] is the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Set scalar solutions at mesh vertices.
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
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in solution
 * structure.
 * (\a pos from 1 to nb_vertices included).
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
 * \param met pointer toward the sol structure.
 * \param sols table of the vectorial solutions
 * sols[3*(i-1)]\@3 is the solution at vertex i
 * \return 0 if failed, 1 otherwise.
 *
 * Set vectorial solutions at mesh vertices
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
 * \param met pointer toward the sol structure.
 * \param m11 value of the tensorial solution at position (1,1) in the tensor.
 * \param m12 value of the tensorial solution at position (1,2) in the tensor.
 * \param m13 value of the tensorial solution at position (1,3) in the tensor.
 * \param m22 value of the tensorial solution at position (2,2) in the tensor.
 * \param m23 value of the tensorial solution at position (2,3) in the tensor.
 * \param m33 value of the tensorial solution at position (3,3) in the tensor.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure.
 * (\a pos from 1 to nb_vertices included).
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
 * \param met pointer toward the sol structure.
 * \param sols table of the tensorial solutions.
 * sols[6*(i-1)]\@6 is the solution at vertex i
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values by array.
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
 * \param sol pointer toward the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s solution(s) at mesh vertex \a pos.
 * \param pos index of the vertex on which we set the solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Set values of the solution at the ith field of the solution array and at
 * position \pos (\a pos from 1 to nb_vertices included and i from 1 to nb_sols).
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
  LIBMMGS_EXPORT int  MMGS_Set_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);
/**
 * \param sol pointer toward the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s table of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Set values of the solution at the ith field of the solution array (\a i from
 * 1 to nb_sols).
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
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Check if the number of given entities match with mesh and sol size
 * (not mandatory) and check mesh datas.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure (unused).
 * \param iparam integer parameter to set (see \a MMGS_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure (unused).
 * \param dparam double parameter to set (see \a MMGS_Param structure).
 * \param val value of the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param hmin minimal edge size.
 * \param hmax maximal edge size.
 * \param hausd value of the Hausdorff number.
 * \return 0 if failed, 1 otherwise.
 *
 * Set local parameters: set the hausdorff value at \a hausd, the minmal edge
 * size value at \a hmin and the maximal edge size value at \a hmax for all
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param ref input tetra reference.
 * \param split MMG5_MMAT_NoSplit if the entity must not be splitted, MMG5_MMAT_Split otherwise
 * \param rin internal reference after ls discretization
 * \param rex external reference after ls discretization
 * \return 0 if failed, 1 otherwise.
 *
 * Set the reference mapping for the elements of ref \a ref in ls discretization mode.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SET_MULTIMAT(mesh,sol,ref,split,rin,rex,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >     INTEGER(MMG5F_INT), INTENT(IN):: ref,rin,rex\n
 * >     INTEGER, INTENT(IN)           :: split\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMGS_EXPORT int  MMGS_Set_multiMat(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int ref,
                                        int split,MMG5_int rin, MMG5_int rex);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param br new level-set base reference.
 * \return 0 if failed, 1 otherwise.
 *
 * Set a new level-set base reference of ref \a br in ls discretization
 * mode. Based references are boundary conditions to which implicit domain can
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

/** recover datas */
/**
 * \param mesh pointer toward the mesh structure.
 * \param np pointer toward the number of vertices.
 * \param nt pointer toward the number of triangles.
 * \param na pointer toward the number of edges.
 * \return 1.
 *
 * Get the number of vertices, triangles and edges of the mesh.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param typEntity pointer toward the type of entities to which solutions are applied.
 * \param np pointer toward the number of solutions.
 * \param typSol pointer toward the type of the solutions (scalar, vectorial...)
 * \return 1.
 *
 * Get the solution number, dimension and type.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward an array of sol structure.
 * \param nsols number of solutions per entity
 * \param nentities pointer toward the number of entities.
 * \param typSol array of size MMG5_NSOL_MAX to store type of each solution
 * (scalar, vector..).
 *
 * \return 1.
 *
 * Get the solution number, dimension and type.
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
 * \param mesh pointer toward the mesh structure.
 * \param c0 pointer toward the coordinate of the point along the first dimension.
 * \param c1 pointer toward the coordinate of the point along the second dimension.
 * \param c2 pointer toward the coordinate of the point along the third dimension.
 * \param ref pointer to the point reference.
 * \param isCorner pointer toward the flag saying if point is corner.
 * \param isRequired pointer toward the flag saying if point is required.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1,\a c2 and reference \a ref of next
 * vertex of mesh.
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
 * \param mesh pointer toward the mesh structure.
 * \param c0 pointer toward the coordinate of the point along the first dimension.
 * \param c1 pointer toward the coordinate of the point along the second dimension.
 * \param c2 pointer toward the coordinate of the point along the third dimension.
 * \param ref pointer to the point reference.
 * \param isCorner pointer toward the flag saying if point is corner.
 * \param isRequired pointer toward the flag saying if point is required.
 * \param idx index of point to get.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1, \a c2 and reference \a ref of
 * vertex \a idx of mesh.
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
 * >   SUBROUTINE MMGS_GET_TRIANGLE(mesh,v0,v1,v2,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)  :: mesh\n
 * >     INTEGER(MMG5F_INT), INTENT(OUT):: v0,v1,v2,ref\n
 * >     INTEGER                        :: isRequired\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int  MMGS_Get_triangle(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref,
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
 * \param met pointer toward the sol structure.
 * \param s pointer toward the scalar solution value.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solution \a s of next vertex of mesh.
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
 * \param met pointer toward the sol structure.
 * \param s table of the scalar solutions at mesh vertices. s[i-1] is
 * the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Get solutions at mesh vertices.
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
 * \param met pointer toward the sol structure.
 * \param vx x value of the vectorial solution.
 * \param vy y value of the vectorial solution.
 * \param vz z value of the vectorial solution.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vectorial solution \f$(v_x,v_y,vz)\f$ of next vertex of mesh.
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
 * \param met pointer toward the sol structure.
 * \param sols table of the solutions at mesh vertices. sols[3*(i-1)]\@3 is
 * the solution at vertex i.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vectorial solutions at mesh vertices
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
 * >   SUBROUTINE MMGS_GET_TENSORSOL(met,m11,m12,m13,m22,m23,m33,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: met\n
 * >     REAL(KIND=8), INTENT(OUT)     :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
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
 * >   SUBROUTINE MMGS_GET_TENSORSOLS(met,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)           :: met\n
 * >     REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: sols\n
 * >     INTEGER, INTENT(OUT)                    :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_Get_tensorSols(MMG5_pSol met, double *sols);
/**
 * \param sol pointer toward the array of solutions
 * \param i position of the solution field that we want to set.
 * \param s solution(s) at mesh vertex \a pos.
 * \param pos index of the vertex on which we get the solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Get values of the ith field of the solution array at vertex \a pos.
 * (\a pos from 1 to nb_vertices included and \a i from 1 to nb_sols).
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
 LIBMMGS_EXPORT int  MMGS_Get_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos);
/**
 * \param sol pointer toward the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s table of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Get values of the solution at the ith field of the solution array.
 * (\a i from 1 to \a nb_sols).
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
  LIBMMGS_EXPORT int  MMGS_Get_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double* s);
/**
 * \param mesh pointer toward the mesh structure.
 * \param iparam integer parameter to set (see \a MMGS_Param structure).
 * \return The value of integer parameter.
 *
 * Get the value of integer parameter \a iparam.
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
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and 0 or 1 data field at VTK vtp file format (.vtp extension). We
 * read only low-order points, edges, tria and quad.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTPMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtpMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and a list of data fields at VTK vtp file format (.vtp extension). We
 * read only low-order points, edges, tria and quad.
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
LIBMMGS_EXPORT int MMGS_loadVtpMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and 0 or 1 data field at VTK vtu file format (.vtu extension). We
 * read only low-order points, edges, tria and quad.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTUMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtuMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and a list of data field at VTK vtu file format (.vtu extension). We
 * read only low-order points, edges, tria and quad.
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
LIBMMGS_EXPORT int MMGS_loadVtuMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and 0 or 1 data field at VTK vtk file format (.vtk extension). We
 * read only low-order points, edges, tria and quad.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADVTKMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadVtkMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and a list of data field at VTK vtk file format (.vtk extension). We
 * read only low-order points, edges, tria and quad.
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
LIBMMGS_EXPORT int MMGS_loadVtkMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and 0 or 1 data field at MSH file format (.msh extension). We read
 * only low-order points, edges, tria, quad, tetra and prisms.
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
LIBMMGS_EXPORT int MMGS_loadMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a list of solution structures.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh and a list of data at MSH file format (.msh extension). We read only
 * low-order points, edges, tria, quadra, tetra and prisms.
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
LIBMMGS_EXPORT int MMGS_loadMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Read mesh data.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_LOADGENERICMESH(mesh,sol,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen0\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT int MMGS_loadGenericMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data field at MSH  file format (.msh extension).
 * Save file at ASCII format for .msh extension, at binary format for .mshb one.
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
LIBMMGS_EXPORT int MMGS_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields (that are considered as solutions and
 * not metrics, thus, we do nothing over the ridge points) at MSH file format
 * (.msh extension).  Save file at ASCII format for .msh extension, at binary
 * format for .mshb one.
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
LIBMMGS_EXPORT int MMGS_saveMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data at Vtk file format (.vtk extension).
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
LIBMMGS_EXPORT int MMGS_saveVtkMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields at Vtk file format (.vtk extension).
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
LIBMMGS_EXPORT int MMGS_saveVtkMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data at vtu Vtk file format (.vtu extension).
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
LIBMMGS_EXPORT int MMGS_saveVtuMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields at vtu Vtk file format (.vtu extension).
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
  LIBMMGS_EXPORT int MMGS_saveVtuMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data at polydata Vtk file format (.vtp extension).
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
  LIBMMGS_EXPORT int MMGS_saveVtpMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields at polydata Vtk file format (.vtp extension).
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
  LIBMMGS_EXPORT int MMGS_saveVtpMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data in a file whose format depends on the filename extension.
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
  LIBMMGS_EXPORT int MMGS_saveGenericMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Load metric field. The solution file (at medit file format) must contains
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
  LIBMMGS_EXPORT int  MMGS_loadSol(MMG5_pMesh mesh,MMG5_pSol met, const char* filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solutions array
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Load 1 or more solutions in a solution file at medit file format.
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
LIBMMGS_EXPORT int  MMGS_loadAllSols(MMG5_pMesh mesh,MMG5_pSol *sol, const char* filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric at medit file format.
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
LIBMMGS_EXPORT int  MMGS_saveSol(MMG5_pMesh mesh, MMG5_pSol met, const char *filename);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solutions array
 * \param filename name of the solution file.
 * \return 0 or -1 if fail, 1 otherwise.
 *
 * Save 1 or more solutions in a solution file at medit file format.
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
LIBMMGS_EXPORT int MMGS_saveAllSols(MMG5_pMesh  mesh,MMG5_pSol *sol ,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward an array of solution structure (that stores solution fields).
 * \return 1
 *
 * Deallocation of an array of solution fields
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
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments.
 *
 * For the MMGS_mmgslib function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the MMGS_mmgsls function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here,\a your_mesh is a \a MMG5_pMesh, \a your_metric and \a your_level_set
 * are \a MMG5_pSol.
 *
 * \return 0 if fail, 1 if success
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
LIBMMGS_EXPORT int MMGS_Free_all(const int starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments.
 *
 * For the MMGS_mmgslib function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the MMGS_mmgsls function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here,\a your_mesh is a \a MMG5_pMesh, \a your_metric and \a your_level_set
 * are \a MMG5_pSol.
 *
 * Here, \a your_mesh is a pointer toward \a MMG5_pMesh and \a your_metric and
 * \a your_level_set a pointer toward \a MMG5_pSol.
 *
 * \return 0 if fail, 1 if success
 *
 * Structure deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
LIBMMGS_EXPORT int MMGS_Free_structures(const int starter,...);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments.
 *
 * For the MMGS_mmgslib function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet,
 * &your_metric,MMG5_ARG_end).
 *
 * For the MMGS_mmgsls function, you need
 * to call the \a MMGS_Init_mesh function with the following arguments :
 * MMGS_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppLs,
 * &your_level_set,MMG5_ARG_end).
 *
 * Here,\a your_mesh is a \a MMG5_pMesh, \a your_metric and \a your_level_set
 * are \a MMG5_pSol.
 *
 * \return 0 if fail, 1 if success
 *
 * Structure deallocations before return.
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
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol (metric) structure.
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
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
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol (level-set) structure.
 * \param met pointer toward the sol (metric) structure (optionnal).
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if fail but a
 * conform mesh is saved or \ref MMG5_STRONGFAILURE if fail and we can't save
 * the mesh.
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
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the sol structure (unused).
 *
 * Set function pointers for caltet, lenedg, defsiz and gradsiz.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_SETFUNC(mesh,met)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void  MMGS_setfunc(MMG5_pMesh mesh,MMG5_pSol met);

/**
 * \param mesh pointer toward the mesh structure.
 * \param nb_edges pointer toward the number of non boundary edges.
 * \return 0 if failed, 1 otherwise.
 *
 * Get the number of non boundary edges (for DG methods for example). An edge is
 * boundary if it is located at the interface of 2 domains with different
 * references, if it belongs to one triangle only or if it is a singular edge
 * (ridge or required).
 * Append these edges to the list of edge.
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
 * \param mesh pointer toward the mesh structure.
 * \param e0 pointer toward the first extremity of the edge.
 * \param e1 pointer toward the second  extremity of the edge.
 * \param ref pointer toward the edge reference.
 * \param idx index of the non boundary edge to get (between 1 and nb_edges)
 * \return 0 if failed, 1 otherwise.
 *
 * Get extremities \a e0, \a e1 and reference \a ref of the idx^th non boundary
 * edge (for DG methods for example). An edge is boundary if it is located at
 * the interface of 2 domains witch different references, if it belongs to one
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
  LIBMMGS_EXPORT int MMGS_Get_nonBdyEdge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref, MMG5_int idx);


/* Tools for the library */
/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the sol structure
 * \return 1 if success
 *
 * Compute isotropic size map according to the mean of the length of the
 * edges passing through a point.
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
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the sol structure
 * \return 1 if success
 *
 * Compute constant size map according to mesh->info.hsiz, mesh->info.hmin and
 * mesh->info.hmax. Update this 3 value if not compatible.
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
 * \param prog pointer toward the program name.
 *
 * Print help for mmgs options.
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
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param sol pointer toward a level-set or displacement
 *
 * \return 1.
 *
 * Store command line arguments.
 *
 * \remark no matching fortran function.
 *
 */
LIBMMGS_EXPORT int  MMGS_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol);
/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
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
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 * \return 1.
 *
 * Store the info structure in the mesh structure.
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
 * \param mesh pointer toward the mesh structure.
 * \param info pointer toward the info structure.
 *
 * Recover the info structure stored in the mesh structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_DESTOCKOPTIONS(mesh,info)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,info\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void MMGS_destockOptions(MMG5_pMesh mesh, MMG5_Info *info);

/**
 * \brief Return adjacent elements of a triangle.
 * \param mesh pointer toward the mesh structure.
 * \param kel triangle index.
 * \param listri pointer toward the table of the indices of the three adjacent
 * triangles of the elt \a kel (the index is 0 if there is no adjacent).
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
 * \brief Return adjacent elements of a triangle.
 * \param mesh pointer toward the mesh structure.
 * \param ip vertex index.
 * \param start index of a triangle holding \a ip.
 * \param lispoi pointer toward an array of size MMGS_LMAX that will contain
 * the indices of adjacent vertices to the vertex \a ip.
 * \return nbpoi the number of adjacent points if success, 0 if fail.
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
 * \param m upper part of a symetric matric diagonalizable in |R
 * \param lambda array of the metric eigenvalues
 * \param vp array of the metric eigenvectors
 *
 * \return the order of the eigenvalues
 *
 * Compute the real eigenvalues and eigenvectors of a symetric matrice m whose
 * upper part is provided (m11, m12, m13, m22, m23, m33 in this order).
 * lambda[0] is the eigenvalue associated to the eigenvector ( v[0][0], v[0,1], v[0,2] )
 * in C and to the eigenvector v(1,:) in fortran
 * lambda[1] is the eigenvalue associated to the eigenvector ( v[1][0], v[1,1], v[1,2] )
 * in C and to the eigenvector v(2,:) in fortran
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
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the solution structure
 *
 * Free the solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMGS_FREE_SOLUTIONS(mesh,sol)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMGS_EXPORT void MMGS_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol);

/**
 * \param mesh pointer toward mesh sructure
 *
 * \return 1 if successful, 0 otherwise.
 *
 * Clean data (triangles and edges) linked to isosurface.
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
 * Set common pointer functions between mmgs and mmg3d to the matching mmgs
 * functions.
 */
LIBMMGS_EXPORT void MMGS_Set_commonFunc(void);

#ifdef __cplusplus
}
#endif

#endif
