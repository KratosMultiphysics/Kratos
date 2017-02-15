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
 * \file common/libmmgcommon.h
 * \brief API header for the common part of the MMG libraries.
 * \author Algiane Froehly  (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning To keep the genheader working, don't break line between the enum
 * name and the opening brace (it creates errors under windows)
 */

#ifndef _MMGLIBCOMMON_H
#define _MMGLIBCOMMON_H

#include <stdarg.h>

#include "chrono.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "libmmgtypes.h"

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
void  MMG5_Init_fileNames(MMG5_pMesh mesh, MMG5_pSol sol);
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
void  (_MMG5_Init_parameters)(MMG5_pMesh mesh);

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
int  MMG5_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin);
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
int  MMG5_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout);
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
int  MMG5_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin);
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
int  MMG5_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout);

/* deallocations */
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 *
 * File name deallocations before return.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SET_MMGFREE_NAMES(mesh,met)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,met\n
 * >   END SUBROUTINE\n
 *
 */
void MMG5_mmgFree_names(MMG5_pMesh mesh, MMG5_pSol met);
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and sol at MSH  file format (.msh extension).
 * Write binary file for .mshb extension.and ASCII for .msh one.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE MMG5_SAVEMSHMESH(mesh,sol,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: mesh,sol\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int MMG5_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename);

#ifdef __cplusplus
}
#endif

#endif
