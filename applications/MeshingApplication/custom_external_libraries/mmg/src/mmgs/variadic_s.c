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
 * \file mmgs/variadic_s.c
 * \brief C variadic functions definitions for MMGS library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmgs/libmmgs.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * variadic functions definitions for MMGS library.
 *
 */

#include "mmgs.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Allocate the mesh and solutions structures at \a MMGS format.
 *
 */
static inline
void _MMGS_Alloc_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol) {

  /* mesh allocation */
  if ( *mesh )  _MMG5_SAFE_FREE(*mesh);
  _MMG5_SAFE_CALLOC(*mesh,1,MMG5_Mesh);

  /* sol allocation */
  if ( !sol ) {
    fprintf(stderr,"  ## Error: an allocatable solution structure of type \"MMG5_pSol\""
           " is needed.\n");
    fprintf(stderr,"            Exit program.\n");
    exit(EXIT_FAILURE);
  }

  if ( *sol )  _MMG5_DEL_MEM(*mesh,*sol,sizeof(MMG5_Sol));
  _MMG5_SAFE_CALLOC(*sol,1,MMG5_Sol);

  return;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * Initialization of mesh and solution structures to their default
 * values (default names, versions, dimensions...).
 *
 */
static inline
void _MMGS_Init_woalloc_mesh(MMG5_pMesh mesh, MMG5_pSol sol ) {

  _MMGS_Set_commonFunc();

  (mesh)->dim  = 3;
  (mesh)->ver  = 2;
  (sol)->dim   = 3;
  (sol)->ver   = 2;
  (sol)->size  = 1;

  /* Default parameters values */
  MMGS_Init_parameters(mesh);

  /* Default vaules for file names */
  MMGS_Init_fileNames(mesh,sol);

  return;
}

/**
 * \param argptr list of the mmg structures that must be initialized. Each
 * structure must follow one of the \a MMG5_ARG* preprocessor variable that allow
 * to identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * Internal function for structure allocations (taking a va_list argument).
 *
 */
void _MMGS_Init_mesh_var( va_list argptr ) {
  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  sol = NULL;


  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case MMG5_ARG_ppMet: case MMG5_ARG_ppLs:
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"  ## Error: MMGS_Init_mesh:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs.\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMGS_Init_mesh:\n"
            " you need to initialize the mesh structure that"
            " will contain your mesh.\n");
    exit(EXIT_FAILURE);
  }

  if ( !sol ) {
    fprintf(stderr,"  ## Error: MMGS_Init_mesh:\n"
            " you need to initialize a solution structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or the"
            " MMG5_ARG_ppLs preprocessor variable) that will contain the output"
            " mesh metric informations, and the input one, if provided.\n.");
    exit(EXIT_FAILURE);
  }

  /* allocations */
  _MMGS_Alloc_mesh(mesh,sol);

  /* initialisations */
  _MMGS_Init_woalloc_mesh(*mesh,*sol);

  return;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the \a MMG5_ARG preprocessor variable that allow to
 * identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * Internal function for deallocations before return (taking a va_list as
 * argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 */
void _MMGS_Free_all_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  sol = NULL;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet): case(MMG5_ARG_ppLs):
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"  ## Error: MMGS_Free_all:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet or "
              "MMG5_ARG_ppLs.\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMGS_Free_all:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n");
    exit(EXIT_FAILURE);
  }

  if ( !sol ) {
    fprintf(stderr,"  ## Error: MMGS_Free_all:\n"
            " you need to provide your metric structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or"
            " the MMG5_ARG_ppLs preprocessor variable)"
            " to allow to free the associated memory.\n");
  }


  MMGS_Free_structures(MMG5_ARG_start,
                        MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                        MMG5_ARG_end);

  _MMG5_SAFE_FREE(*mesh);

  if ( sol )
    _MMG5_SAFE_FREE(*sol);

  return;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the \a MMG5_ARG* preprocessor variable that allow
 * to identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * Internal function for structures deallocations before return (taking a
 * va_list as argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void _MMGS_Free_structures_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  sol = NULL;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet): case(MMG5_ARG_ppLs):
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"  ## Error: MMGS_Free_structures:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet or"
              " MMG5_ARG_ppLs.\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMGS_Free_structures:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n");
    exit(EXIT_FAILURE);
  }

  MMGS_Free_names(MMG5_ARG_start,
                   MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                   MMG5_ARG_end);

 /* mesh */
  assert(mesh && *mesh);
  if ( (*mesh)->point )
    _MMG5_DEL_MEM((*mesh),(*mesh)->point,((*mesh)->npmax+1)*sizeof(MMG5_Point));

  if ( (*mesh)->edge )
    _MMG5_DEL_MEM((*mesh),(*mesh)->edge,((*mesh)->na+1)*sizeof(MMG5_Edge));

  if ( (*mesh)->adja )
    _MMG5_DEL_MEM((*mesh),(*mesh)->adja,(3*(*mesh)->ntmax+5)*sizeof(int));

  if ( (*mesh)->xpoint )
    _MMG5_DEL_MEM((*mesh),(*mesh)->xpoint,((*mesh)->xpmax+1)*sizeof(MMG5_xPoint));

  if ( (*mesh)->tria )
    _MMG5_DEL_MEM((*mesh),(*mesh)->tria,((*mesh)->ntmax+1)*sizeof(MMG5_Tria));

  /* sol */
  if ( sol && (*sol) && (*sol)->m )
    _MMG5_DEL_MEM((*mesh),(*sol)->m,((*sol)->size*((*sol)->npmax+1))*sizeof(double));

  /* (*mesh)->info */
  if ( (*mesh)->info.npar && (*mesh)->info.par )
    _MMG5_DEL_MEM((*mesh),(*mesh)->info.par,(*mesh)->info.npar*sizeof(MMG5_Par));

  if ( (*mesh)->info.imprim>5 || (*mesh)->info.ddebug )
    printf("  MEMORY USED AT END (bytes) %ld\n",_MMG5_safeLL2LCast((*mesh)->memCur));

  return;
}

/**
 * \param argptr list of the mmg structures for whose we want to deallocate the
 * name. Each structure must follow one of the \a MMG5_ARG preprocessor variable
 * that allow to identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh
 * structure (that will contain the mesh and identified by the MMG5_ARG_ppMesh
 * keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * Internal function for name deallocations before return (taking a va_list as
 * argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void _MMGS_Free_names_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  sol = NULL;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet): case(MMG5_ARG_ppLs):
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"  ## Error: MMGS_Free_names:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet "
              " or MMG5_ARG_ppLs\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMGS_Free_names:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n");
    exit(EXIT_FAILURE);
  }
  if ( !sol ) {
    fprintf(stderr,"  ## Error: MMGS_Free_names:\n"
            " you need to provide your metric structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or the "
            " MMG5_ARG_ppLs preprocessor variable)"
            " to allow to free the associated memory.\n");
  }

  /* mesh & met */
  MMG5_mmgFree_names(*mesh,*sol);

  return;
}
