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

#include "libmmgs_private.h"
#include "libmmgs.h"
#include "mmgsexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the sol structure.
 * \param ls pointer to the sol structure.
 *
 * \return 0 if fail, 1 if success
 *
 * Allocate the mesh and solutions structures at \a MMGS format.
 *
 */
static inline
int MMGS_Alloc_mesh(MMG5_pMesh *mesh, MMG5_pSol *met, MMG5_pSol *ls) {

  /* mesh allocation */
  if ( *mesh )  MMG5_SAFE_FREE(*mesh);
  MMG5_SAFE_CALLOC(*mesh,1,MMG5_Mesh,return 0);

 /* metric allocation */
  if ( met ) {
    if ( *met )  MMG5_DEL_MEM(*mesh,*met);
    MMG5_SAFE_CALLOC(*met,1,MMG5_Sol,return 0);
  }

  /* level-set allocation in ls mode */
  if ( ls ) {
    if ( *ls )
      MMG5_DEL_MEM(*mesh,*ls);
    MMG5_SAFE_CALLOC(*ls,1,MMG5_Sol,return 0);
  }

  return 1;
}
/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to a sol structure (metric).
 * \param ls pointer to the level-set (in ls-mode).
 *
 * Initialization of mesh and solution structures to their default
 * values (default names, versions, dimensions...).
 *
 */
static inline
void MMGS_Init_woalloc_mesh(MMG5_pMesh mesh, MMG5_pSol *met,MMG5_pSol *ls ) {

  MMGS_Set_commonFunc();

  (mesh)->dim   = 3;
  (mesh)->ver   = 2;
  (mesh)->nsols = 0;

  if ( met && *met ) {
    (*met)->dim  = 3;
    (*met)->ver  = 2;
    (*met)->size = 1;
    (*met)->type = 1;
  }

  if ( ls && *ls ) {
    (*ls)->dim  = 3;
    (*ls)->ver  = 2;
    (*ls)->size = 1;
    (*ls)->type = 1;
  }

  /* Default parameters values */
  MMGS_Init_parameters(mesh);

  /* Default vaules for file names */
  if ( met ) {
    MMGS_Init_fileNames(mesh,*met);
  }
  else {
    MMGS_Init_fileNames(mesh,NULL);
  }

  if ( ls && *ls ) {
    MMGS_Set_inputSolName(mesh,*ls,"");
    MMGS_Set_outputSolName(mesh,*ls,"");
  }

  return;
}

/**
 * \param argptr list of the mmg structures that must be initialized. Each
 * structure must follow one of the \a MMG5_ARG* preprocessor variable that allow
 * to identify it.
 *
 * \a argptr contains at least a pointer to a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer to a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 * \return 0 if fail, 1 if success
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * Internal function for structure allocations (taking a va_list argument).
 *
 */
int MMGS_Init_mesh_var( va_list argptr ) {
  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*ls;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  mesh = NULL;
  sol = ls = NULL;


  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case MMG5_ARG_ppMet:
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppLs):
      ls = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMGS_Init_mesh:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs.\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMGS_Init_mesh:\n"
            " you need to initialize the mesh structure that"
            " will contain your mesh.\n",__func__);
    return 0;
  }

  /* allocations */
  if ( !MMGS_Alloc_mesh(mesh,sol,ls) )  return 0;

  /* initialisations */
  MMGS_Init_woalloc_mesh(*mesh,sol,ls);

  return 1;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the \a MMG5_ARG preprocessor variable that allow to
 * identify it.
 *
 * \a argptr contains at least a pointer to a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer to a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * \return 0 if fail, 1 if success
 *
 * Internal function for deallocations before return (taking a va_list as
 * argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 */
int MMGS_Free_all_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*sols,*ls;
  int            typArg;
  int            meshCount,metCount,lsCount,fieldsCount;

  meshCount = metCount = lsCount = fieldsCount = 0;
  sol = sols = ls = NULL;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet):
      ++metCount;
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppLs):
      ++lsCount;
      ls = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppSols):
      ++fieldsCount;
      sols = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMGS_Free_all:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet or "
              "MMG5_ARG_ppLs.\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMGS_Free_all:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  if ( metCount > 1 || lsCount > 1 || fieldsCount > 1 ) {
    fprintf(stdout,"\n  ## Warning: %s: MMGS_Free_all:\n"
            " This function can free only one structure of each type.\n"
            " Probable memory leak.\n",
            __func__);
  }

  if ( !MMGS_Free_structures(MMG5_ARG_start,
                             MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                             MMG5_ARG_ppLs, ls,MMG5_ARG_ppSols, sols,
                             MMG5_ARG_end) )
    return 0;

  if ( sol )
    MMG5_SAFE_FREE(*sol);

  if ( sols ) {
    MMG5_DEL_MEM(*mesh,*sols);
  }

  MMG5_SAFE_FREE(*mesh);

  return 1;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the \a MMG5_ARG* preprocessor variable that allow
 * to identify it.
 *
 * \a argptr contains at least a pointer to a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer to a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * \return 0 if fail, 1 if success
 *
 * Internal function for structures deallocations before return (taking a
 * va_list as argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
int MMGS_Free_structures_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*ls,*sols;
  int            typArg,i;
  int            meshCount;

  meshCount = 0;
  sol = ls = sols = NULL;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet):
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppLs):
      ls = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppSols):
      sols = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMGS_Free_structures:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet or"
              " MMG5_ARG_ppLs.\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMGS_Free_structures:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  MMGS_Free_names(MMG5_ARG_start,
                  MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                  MMG5_ARG_ppLs, ls, MMG5_ARG_ppSols, sols,
                  MMG5_ARG_end);

  /* mesh */
  assert(mesh && *mesh);

  /* ls */
  if ( ls && (*ls) && (*ls)->m ) {
    MMG5_DEL_MEM((*mesh),(*ls)->m);
  }

  /* met */
  if ( sol && (*sol) && (*sol)->m ) {
    MMG5_DEL_MEM((*mesh),(*sol)->m);
  }

  /* field */
  if ( sols && (*mesh)->nsols ) {
    for ( i=0; i<(*mesh)->nsols; ++i ) {
      MMG5_DEL_MEM((*mesh),(*sols)[i].m);
    }
  }

  MMG5_Free_structures(*mesh,NULL);

  return 1;
}

/**
 * \param argptr list of the mmg structures for whose we want to deallocate the
 * name. Each structure must follow one of the \a MMG5_ARG preprocessor variable
 * that allow to identify it.
 *
 * \a argptr contains at least a pointer to a \a MMG5_pMesh
 * structure (that will contain the mesh and identified by the MMG5_ARG_ppMesh
 * keyword).
 *
 *  To call the \a MMGS_mmgslib function, you must also provide
 * a pointer to a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMGS_mmgsls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 * \return 0 if fail, 1 if success
 *
 * Internal function for name deallocations before return (taking a va_list as
 * argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
int MMGS_Free_names_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      psl,*sol,*ls,*sols;
  int            typArg,i;
  int            meshCount;

  meshCount = 0;
  sol = ls = sols = NULL;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet):
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppLs):
      ls = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppSols):
      sols = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMGS_Free_names:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet "
              " or MMG5_ARG_ppLs\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMGS_Free_names:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  /* mesh & met */
  assert(mesh && *mesh );

  if ( sol ) {
    MMG5_mmgFree_names(*mesh,*sol);
  }
  else {
    MMG5_mmgFree_names(*mesh,NULL);
  }
  /* ls */
  if ( ls && *ls ) {
    if ( (*ls)->namein ) {
      MMG5_DEL_MEM(*mesh,(*ls)->namein);
    }

    if ( (*ls)->nameout ) {
      MMG5_DEL_MEM(*mesh,(*ls)->nameout);
    }
  }
  /* Fields */
  if ( sols ) {
    for ( i=0; i<(*mesh)->nsols; ++i ) {
      psl = (*sols) + i;
      if ( psl->namein ) {
        MMG5_DEL_MEM(*mesh,psl->namein);
      }
      if ( psl->nameout ) {
        MMG5_DEL_MEM(*mesh,psl->nameout);
      }
    }
  }


  return 1;
}
