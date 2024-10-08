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
 * \brief API header for the common part of the MMG libraries.
 * \author Algiane Froehly  (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning To keep the genheader working, don't break line between the enum
 * name and the opening brace (it creates errors under windows)
 */

#ifndef MMGLIBCOMMON_H
#define MMGLIBCOMMON_H

#include <stdarg.h>

#include "libmmgtypes.h"

#include "chrono_private.h"

#include "mmg_core_export_private.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MMG5_VOLFRAC     1.e-5

/*----------------------------- functions header -----------------------------*/
/* Initialization functions */
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialize file names to their default values.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_INIT_FILENAMES(mesh,sol)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh,sol\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG_CORE_EXPORT void  MMG5_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_INIT_PARAMETERS(mesh)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG_CORE_EXPORT void  (MMG5_Init_parameters)(MMG5_pMesh mesh);

/* init file names */
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SET_INPUTMESHNAME(mesh,meshin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG_CORE_EXPORT int  MMG5_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SET_OUTPUTMESHNAME(mesh,meshout,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshout\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG_CORE_EXPORT int  MMG5_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SET_INPUTSOLNAME(mesh,sol,solin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG_CORE_EXPORT int  MMG5_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SET_OUTPUTSOLNAME(mesh,sol,solout,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solout\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG_CORE_EXPORT int MMG5_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param ls pointer toward a solution structure (level-set or displacement).
 *
 * \return 1 if success, 0 if fail (computed bounding box too small
 * or one af the anisotropic input metric is not valid).
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SCALEMESH(mesh,met,ls,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,ls\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG_CORE_EXPORT int MMG5_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol ls);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a metric.
 * \param sol pointer toward a solution structure (level-set or displacement).
 *
 * \return 1.
 *
 * Unscale the mesh and the size informations to their initial sizes.
 *
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_UNSCALEMESH(mesh,met,ls,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met,ls\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  LIBMMG_CORE_EXPORT int MMG5_unscaleMesh(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol ls);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param hsiz wanted edge size
 *
 * fill the metric field with the size \a hsiz
 *
 * \Remark not for extern users.
 *
 */
LIBMMG_CORE_EXPORT void MMG5_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met,double hsiz);

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
 */
LIBMMG_CORE_EXPORT int  MMG5_Set_multiMat(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int ref,int split,
                                          MMG5_int rin, MMG5_int rex);


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param br new level-set base reference.
 * \return 0 if failed, 1 otherwise.
 *
 * Set a new level-set base reference of ref \a br in ls discretization
 * mode. Base references are boundary conditions to which implicit domain can
 * be attached. All implicit volumes that are not attached to listed based
 * references are deleted as spurious volumes by the \a rmc option.
 *
 */
LIBMMG_CORE_EXPORT int  MMG5_Set_lsBaseReference(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_int br);

/* deallocations */
LIBMMG_CORE_EXPORT void MMG5_Free_structures(MMG5_pMesh mesh,MMG5_pSol sol);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * File name deallocations before return.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SETMMGFREE_NAMES(mesh,met)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >   END SUBROUTINE\n
 *
 */
LIBMMG_CORE_EXPORT void MMG5_mmgFree_names(MMG5_pMesh mesh, MMG5_pSol met);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sethmin 1 if hmin is already setted (>0.)
 * \param sethmax 1 if hmax is already setted (>0.)
 *
 * \return 1 if success, 0 if we detect mismatch parameters
 *
 * Set default values for hmin and hmax  from the bounding box.
 *
 * \Remark not for extern users.
 *
 */
LIBMMG_CORE_EXPORT extern int MMG5_Set_defaultTruncatureSizes(MMG5_pMesh mesh,int8_t sethmin,int8_t sethmax);

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric.
 * \param hsiz computed constant size to impose.
 *
 * \return 1 if success, 0 if fail
 *
 * Compute the constant size to impose according to hmin and hmax and store it in \a hsiz.
 * Fill hmin and hamx if they are not setted by the user.
 *
 */
LIBMMG_CORE_EXPORT int MMG5_Compute_constantSize(MMG5_pMesh mesh,MMG5_pSol met,double *hsize);

/**
 * \param tag input entity tag
 *
 * \return the list of the flags contained in \a tag
 *
 * Print the name associated to the \a typ value in the \a MMG5_type enum.
 *
 * \warning for debug purpose, no thread safe.
 */
const char* MMG5_Get_tagName(int tag);

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward an array of solution structure (that stores solution fields).
 * \return 1
 *
 * Deallocation of an array of solution fields
 *
 */
LIBMMG_CORE_EXPORT int MMG5_Free_allSols(MMG5_pMesh mesh,MMG5_pSol *sol);

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 *
 * \return 1 if success, 0 if fail.
 *
 * Save node list at .node file format (Tetgen/Triangle).
 */
LIBMMG_CORE_EXPORT int MMG5_saveNode(MMG5_pMesh mesh,const char *filename);

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \param ext file extension (.poly or .edge)
 *
 * \return 1 if success, 0 if fail.
 *
 * Save edge list at .edge file format (Tetgen/Triangle).
 */
LIBMMG_CORE_EXPORT int MMG5_saveEdge(MMG5_pMesh mesh,const char *filename,const char *ext);

/* Useful tools to manage C strings */
/**
 * \param path string containing a filename and its path
 *
 * \return a pointer toward the allocated string that contains the file basename.
 *
 * Extract basename from a path (allocate a string to store it).
 *
 */
LIBMMG_CORE_EXPORT char *MMG5_Get_basename(char *path);

/**
 * \param ptr pointer toward the file extension (dot included)
 * \param fmt default file format.
 *
 * \return and index associated to the file format detected from the extension.
 *
 * Get the wanted file format from the mesh extension. If \a fmt is provided, it
 * is used as default file format (\a ptr==NULL), otherwise, the default file
 * format is the medit one.
 *
 */
LIBMMG_CORE_EXPORT int MMG5_Get_format( char *ptr, int fmt );

/**
 * \param filename string containing a filename
 *
 * \return pointer toward the filename extension or toward the end of the string
 * if no extension have been founded
 *
 * Get the extension of the filename string. Do not consider '.o' as an extension.
 *
 */
LIBMMG_CORE_EXPORT char *MMG5_Get_filenameExt( char *filename );

/**
 * \param path string containing a filename and its path
 *
 * \return a pointer toward the path allocated here
 *
 * Remove filename from a path and return the path in a newly allocated string.
 *
 */
LIBMMG_CORE_EXPORT char *MMG5_Get_path(char *path);

/**
 * \param path path from which we want to remove the extension.
 *
 * \return allocated string or NULL if the allocation fail.
 *
 * Allocate a new string and copy \a path without extension in it.
 *
 */
LIBMMG_CORE_EXPORT char *MMG5_Remove_ext (char* path,char *ext);

/* Enum utilities: print enum fields under a string form */
/**
 * \param fmt file format.
 *
 * \return The name of the file format in a string.
 *
 * Print the name of the file format associated to \a fmt.
 *
 */
LIBMMG_CORE_EXPORT const char* MMG5_Get_formatName(enum MMG5_Format fmt);

/**
 * \param ent MMG5_entities enum
 *
 * \return the name of the enum field
 *
 * Print the name associated to the \a ent value in the \a MMG5_entities enum.
 *
 */
LIBMMG_CORE_EXPORT const char* MMG5_Get_entitiesName(enum MMG5_entities ent);

/**
 * \param typ MMG5_type enum
 *
 * \return the name of the enum field
 *
 * Print the name associated to the \a typ value in the \a MMG5_type enum.
 *
 */
LIBMMG_CORE_EXPORT const char* MMG5_Get_typeName(enum MMG5_type typ);


/**
 * \param mesh pointer toward mesh
 *
 * \return 1 if successful, 0 otherwise
 *
 * Clean non-ridge edges belonging to isosurface.
 *
 */
LIBMMG_CORE_EXPORT int MMG5_Clean_isoEdges(MMG5_pMesh mesh);

#ifdef __cplusplus
}
#endif

#endif
