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
 * \file mmg2d/variadic_2d.c
 * \brief C variadic functions definitions for MMG2D library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmg2d/libmmg2d.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * variadic functions definitions for MMG2D library.
 *
 */

#include "mmg2d.h"
#include "mmg2dexterns.c"

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 *
 * \return 1 if success, 0 if fail
 *
 * Allocate the mesh and solutions structures at \a MMG2D format.
 *
 */
static inline
int MMG2D_Alloc_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol, MMG5_pSol *disp) {

  /* mesh allocation */
  if ( *mesh )  MMG5_SAFE_FREE(*mesh);
  MMG5_SAFE_CALLOC(*mesh,1,MMG5_Mesh,return 0);

  /* sol allocation */
  if ( !sol ) {
    fprintf(stderr,"\n  ## Error: %s: an allocatable solution"
            " structure of type \"MMG5_pSol\"" " is needed.\n",__func__);
    fprintf(stderr,"            Exit program.\n");
    return 0;
  }

  if ( *sol )  MMG5_DEL_MEM(*mesh,*sol);
  MMG5_SAFE_CALLOC(*sol,1,MMG5_Sol,return 0);

  /* Displacement allocation */
  if ( disp ) {
    if ( *disp )
      MMG5_DEL_MEM(*mesh,*disp);
    MMG5_SAFE_CALLOC(*disp,1,MMG5_Sol,return 0);
  }

  return 1;
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a sol structure (metric or level-set).
 * \param disp pointer toward a sol structure (displacement).
 *
 * Initialization of mesh and solution structures to their default
 * values (default names, versions, dimensions...).
 *
 */
static inline
void MMG2D_Init_woalloc_mesh(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSol disp) {

  MMG2D_Set_commonFunc();

  (mesh)->dim   = 2;
  (mesh)->ver   = 2;
  (mesh)->nsols = 0;

  (sol)->dim  = 2;
  (sol)->ver  = 2;
  (sol)->size = 1;

  if ( disp ) {
    (disp)->dim  = 2;
    (disp)->ver  = 2;
    (disp)->size = 2;
  }

  /* Default parameters values */
  MMG2D_Init_parameters(mesh);

  /* Default vaules for file names */
  MMG2D_Init_fileNames(mesh,sol);

  if ( disp ) {
    MMG2D_Set_inputSolName(mesh,disp,"");
    MMG2D_Set_outputSolName(mesh,disp,"");
  }

  return;
}

/**
 * \param argptr list of the mmg structures that must be initialized. Each
 * structure must follow one of the \a MMG5_ARG* preprocessor variable that allow to identify
 * it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword)
 *
 *  To call the \a MMG2D_mmg2dlib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMG2D_mmg2dls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 *  To call the \a MMG2D_mmg2dmov library, you must also provide a
 * pointer toward a \a MMG5_pSol structure storing the displacement (and
 * identified by the MMG5_ARG_ppDisp keyword).
 *
 * \return 0 if fail, 1 otherwise
 *
 * Internal function for structure allocations (taking a va_list argument).
 *
 */
int MMG2D_Init_mesh_var( va_list argptr ) {
  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  mesh = NULL;
  disp = sol = NULL;


  while ( (typArg = va_arg(argptr,int          )) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet): case MMG5_ARG_ppLs:
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG2D_Init_mesh:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              "  MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG2D_Init_mesh:\n"
            " you need to initialize the mesh structure that"
            " will contain your mesh.\n",__func__);
    return 0;
  }

  if ( !sol ) {
    fprintf(stderr,"\n  ## Error: %s: MMG2D_Init_mesh:\n"
            " you need to initialize a solution structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or"
            " MMG5_ARG_ppLs preprocessor variable)"
            " that will contain the output mesh metric"
            " informations, and the input one, if provided.\n.",__func__);
    return 0;
  }

  /* allocations */
  if ( !MMG2D_Alloc_mesh(mesh,sol,disp) )  return 0;

  /* initialisations */
  if ( disp )
    MMG2D_Init_woalloc_mesh(*mesh,*sol,*disp);
  else
    MMG2D_Init_woalloc_mesh(*mesh,*sol,NULL);

  return 1;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the \a MMG5_ARG preprocessor variable that allow to
 * identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword)
 *
 *  To call the \a MMG2D_mmg2dlib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMG2D_mmg2dls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 *  To call the \a MMG2D_mmg2dmov library, you must also provide a
 * pointer toward a \a MMG5_pSol structure storing the displacement (and
 * identified by the MMG5_ARG_ppDisp keyword).
 *
 * \return 0 if fail, 1 if success
 *
 * Internal function for deallocations before return (taking a va_list as
 * argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 */
int MMG2D_Free_all_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      psl,*sol,*disp,*sols;
  int            typArg;
  int            meshCount,ier,i;

  meshCount = 0;
  disp = sol = sols = NULL;

  while ( (typArg = va_arg(argptr,int          )) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet):case(MMG5_ARG_ppLs):
      sol = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppSols):
      sols = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG2D_Free_all:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh or MMG5_ARG_ppMet\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG2D_Free_all:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  if ( !disp )
    ier = MMG2D_Free_structures(MMG5_ARG_start,
                                MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                                MMG5_ARG_end);
  else
    ier = MMG2D_Free_structures(MMG5_ARG_start,
                                MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                                MMG5_ARG_ppDisp, disp,
                                MMG5_ARG_end);

  if ( sol )
    MMG5_SAFE_FREE(*sol);

  if ( disp )
    MMG5_SAFE_FREE(*disp);

  if ( sols ) {
    for ( i=0; i<(*mesh)->nsols; ++i ) {
      psl = (*sols) + i;
      if ( psl->m ) {
        MMG5_DEL_MEM(*mesh,psl->m);
      }
    }
    MMG5_DEL_MEM(*mesh,*sols);
  }

  MMG5_SAFE_FREE(*mesh);

  return ier;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the \a MMG5_ARG* preprocessor variable that allow
 * to identify it.  \a argptr contains at least a pointer toward a \a MMG5_pMesh
 * structure (that will contain the mesh and identified by the MMG5_ARG_ppMesh
 * keyword) and a pointer toward a \a MMG5_pSol structure (that will contain the
 * ouput metric (and the input one, if provided) and identified by the
 * MMG5_ARG_ppMet keyword).
 *
 *  To call the \a MMG2D_mmg2dls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 *  To call the \a MMG2D_mmg2dmov library, you must also provide a
 * pointer toward a \a MMG5_pSol structure storing the displacement (and
 * identified by the MMG5_ARG_ppDisp keyword).
 *
 * \return 1 if success, 0 if fail
 *
 * Internal function for structures deallocations before return (taking a
 * va_list as argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
int MMG2D_Free_structures_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  mesh = NULL;
  disp = sol = NULL;

  while ( (typArg = va_arg(argptr,int          )) != MMG5_ARG_end )
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
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG2D_Free_structures:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh or MMG5_ARG_ppMet\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG2D_Free_structures:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  if ( !disp ) {
    if ( !MMG2D_Free_names(MMG5_ARG_start,
                           MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                           MMG5_ARG_end) )
      return 0;
  }
  else {
    if ( !MMG2D_Free_names(MMG5_ARG_start,
                           MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                           MMG5_ARG_ppDisp, disp,
                           MMG5_ARG_end) )
      return 0;
  }

  /* mesh */
  assert(mesh && *mesh);

  if ( (*mesh)->edge )
    MMG5_DEL_MEM((*mesh),(*mesh)->edge);

  if ( (*mesh)->adja )
    MMG5_DEL_MEM((*mesh),(*mesh)->adja);

  if ( (*mesh)->tria )
    MMG5_DEL_MEM((*mesh),(*mesh)->tria);

  /* disp */
  if ( disp && (*disp) && (*disp)->m )
    MMG5_DEL_MEM((*mesh),(*disp)->m);

  if ( sol ) {
    MMG5_Free_structures(*mesh,*sol);
  }
  else {
    MMG5_Free_structures(*mesh,NULL);
  }

  return 1;
}

/**
 * \param argptr list of the mmg structures for whose we want to deallocate the
 * name. Each structure must follow one of the \a MMG5_ARG* preprocessor variable
 * that allow to identify it.  \a argptr contains at least a pointer toward a \a
 * MMG5_pMesh structure (that will contain the mesh and identified by the
 * MMG5_ARG_ppMesh keyword) and a pointer toward a \a MMG5_pSol structure (that
 * will contain the ouput metric (and the input one, if provided) and identified
 * by the MMG5_ARG_ppMet keyword).
 *
 * Internal function for name deallocations before return (taking a va_list as
 * argument).
 *
 * \return 0 if fail, 1 if success
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
int MMG2D_Free_names_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  disp = sol = NULL;

  while ( (typArg = va_arg(argptr,int          )) != MMG5_ARG_end )
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
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG2D_Free_names:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following"
              " preprocessor variable: MMG5_ARG_ppMesh or MMG5_ARG_ppMet\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG2D_Free_names:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  /* mesh & met */
  if ( sol ) {
    MMG5_mmgFree_names(*mesh,*sol);
  }
  else {
    MMG5_mmgFree_names(*mesh,NULL);
  }

  /* disp */
  if ( disp && *disp ) {
    if ( (*disp)->namein ) {
      MMG5_DEL_MEM(*mesh,(*disp)->namein);
    }

    if ( (*disp)->nameout ) {
      MMG5_DEL_MEM(*mesh,(*disp)->nameout);
    }
  }

  return 1;
}
