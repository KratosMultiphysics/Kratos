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
 * \file mmg3d/variadic_3d.c
 * \brief C variadic functions definitions for MMG3D library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete soon...
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmg3d/libmmg3d.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * variadic functions definitions for MMG3D library.
 *
 */

#include "libmmg3d_private.h"
#include "libmmg3d.h"
#include "mmg3dexterns_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric).
 * \param ls pointer toward the level-set (in ls-mode).
 * \param disp pointer toward a sol structure (displacement).
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate the mesh and solutions structures at \a MMG3D format.
 *
 */
static inline
int MMG3D_Alloc_mesh(MMG5_pMesh *mesh, MMG5_pSol *met, MMG5_pSol *ls, MMG5_pSol *disp
  ) {

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
 * \param met pointer toward a sol structure (metric).
 * \param ls pointer toward the level-set (in ls-mode).
 * \param disp pointer toward a sol structure (displacement).
 *
 * Initialization of mesh and solution structures to their default
 * values (default names, versions, dimensions...).
 *
 */
static inline
void MMG3D_Init_woalloc_mesh(MMG5_pMesh mesh, MMG5_pSol *met,MMG5_pSol *ls, MMG5_pSol *disp
  ) {

  MMG3D_Set_commonFunc();

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

  if ( disp && *disp ) {
    (*disp)->dim  = 3;
    (*disp)->ver  = 2;
    (*disp)->size = 2;
    (*disp)->type = 2;
  }

  /* Default parameters values */
  MMG3D_Init_parameters(mesh);

  /* Default vaules for file names */
  if ( met ) {
    MMG3D_Init_fileNames(mesh,*met);
  }
  else {
    MMG3D_Init_fileNames(mesh,NULL);
  }

  if ( ls && *ls ) {
    MMG3D_Set_inputSolName(mesh,*ls,"");
    MMG3D_Set_outputSolName(mesh,*ls,"");
  }

  if ( disp && *disp ) {
    MMG3D_Set_inputSolName(mesh,*disp,"");
    MMG3D_Set_outputSolName(mesh,*disp,"");
  }

  return;
}

/**
 * \param argptr list of the mmg structures that must be initialized. Each
 * structure must follow one of the MMG5_ARG* preprocessor variable that allow
 * to identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword)
 *
 *  To call the \a MMG3D_mmg3dlib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMG3D_mmg3dls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 *  To call the \a MMG3D_mmg3dmov library, you must also provide a
 * pointer toward a \a MMG5_pSol structure storing the displacement (and
 * identified by the MMG5_ARG_ppDisp keyword).
 *
 * \return 1 if success, 0 if fail
 *
 * Internal function for structure allocations (taking a va_list argument).
 *
 */
int MMG3D_Init_mesh_var( va_list argptr ) {
  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp,*ls;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  mesh = NULL;
  disp = sol = ls = NULL;


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
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG3D_Init_mesh:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one"
              " of the MMG5_ARG* preprocessor variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              "  MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG3D_Init_mesh:\n"
            " you need to initialize the mesh structure that"
            " will contain your mesh.\n",__func__);
    return 0;
  }

  /* allocations */
  if ( !MMG3D_Alloc_mesh(mesh,sol,ls,disp) ) return 0;

  /* initialisations */
  MMG3D_Init_woalloc_mesh(*mesh,sol,ls,disp);

  return 1;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the MMG5_ARG* preprocessor variable that allow
 * to identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword)
 *
 *  To call the \a MMG3D_mmg3dlib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMG3D_mmg3dls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 *  To call the \a MMG3D_mmg3dmov library, you must also provide a
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
int MMG3D_Free_all_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp,*sols,*ls;
  int            typArg;
  int            meshCount,metCount,lsCount,dispCount,fieldsCount;

  meshCount = metCount = lsCount = dispCount = fieldsCount = 0;
  mesh = NULL;
  disp = sol = sols = ls = NULL;

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
    case(MMG5_ARG_ppDisp):
      ++dispCount;
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppSols):
      ++fieldsCount;
      sols = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG3D_Free_all:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following preprocessor"
              " variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG3D_Free_all:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  if ( metCount > 1 || lsCount > 1 || dispCount > 1 || fieldsCount > 1 ) {
    fprintf(stdout,"\n  ## Warning: %s: MMG3D_Free_all:\n"
            " This function can free only one structure of each type.\n"
            " Probable memory leak.\n",
            __func__);
  }

  if ( !MMG3D_Free_structures(MMG5_ARG_start,
                              MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                              MMG5_ARG_ppLs, ls,MMG5_ARG_ppDisp, disp,
                              MMG5_ARG_ppSols, sols,
                              MMG5_ARG_end) ) {
    return 0;
  }

  if ( sol ) {
    MMG5_SAFE_FREE(*sol);
  }

  if ( disp ) {
    MMG5_SAFE_FREE(*disp);
  }

  if ( ls ) {
    MMG5_SAFE_FREE(*ls);
  }

  if ( sols ) {
    MMG5_DEL_MEM(*mesh,*sols);
  }

  MMG5_SAFE_FREE(*mesh);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a solution / level-set.
 * \param sol pointer toward a displacement.
 *
 * Free mesh arrays.
 *
 */
void MMG3D_Free_arrays(MMG5_pMesh *mesh,MMG5_pSol *sol,MMG5_pSol *ls,
                       MMG5_pSol *disp,MMG5_pSol *field){
  int i;
  if ( (*mesh)->tetra )
    MMG5_DEL_MEM((*mesh),(*mesh)->tetra);

  if ( (*mesh)->prism )
    MMG5_DEL_MEM((*mesh),(*mesh)->prism);

  if ( (*mesh)->edge )
    MMG5_DEL_MEM((*mesh),(*mesh)->edge);

  if ( (*mesh)->adjt )
    MMG5_DEL_MEM(*mesh,(*mesh)->adjt);

  if ( (*mesh)->adja )
    MMG5_DEL_MEM((*mesh),(*mesh)->adja);

  if ( (*mesh)->adjapr )
    MMG5_DEL_MEM((*mesh),(*mesh)->adjapr);

  if ( (*mesh)->htab.geom )
    MMG5_DEL_MEM((*mesh),(*mesh)->htab.geom);

  if ( (*mesh)->tria )
    MMG5_DEL_MEM((*mesh),(*mesh)->tria);

  if ( (*mesh)->quadra )
    MMG5_DEL_MEM((*mesh),(*mesh)->quadra);

  if ( (*mesh)->xtetra )
    MMG5_DEL_MEM((*mesh),(*mesh)->xtetra);

  if ( (*mesh)->xprism )
    MMG5_DEL_MEM((*mesh),(*mesh)->xprism);

  /* disp */
  if ( disp && (*disp) && (*disp)->m )
    MMG5_DEL_MEM((*mesh),(*disp)->m);

  /* ls */
  if ( ls && (*ls) && (*ls)->m )
    MMG5_DEL_MEM((*mesh),(*ls)->m);

  /* field */
  if ( field && (*mesh)->nsols ) {
    for ( i=0; i<(*mesh)->nsols; ++i ) {
      MMG5_DEL_MEM((*mesh),(*field)[i].m);
    }
  }

  if ( sol ) {
    MMG5_Free_structures(*mesh,*sol);
  }
  else {
    MMG5_Free_structures(*mesh,NULL);
  }
  return;
}

/**
 * \param argptr list of the mmg structures that must be deallocated. Each
 * structure must follow one of the MMG5_ARG* preprocessor variable that allow
 * to identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword)
 *
 *  To call the \a MMG3D_mmg3dlib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMG3D_mmg3dls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 *  To call the \a MMG3D_mmg3dmov library, you must also provide a
 * pointer toward a \a MMG5_pSol structure storing the displacement (and
 * identified by the MMG5_ARG_ppDisp keyword).
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
int MMG3D_Free_structures_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*ls,*disp,*sols;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  mesh = NULL;
  disp = sol = ls = sols = NULL;

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
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppSols):
      sols = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG3D_Free_structures:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following preprocessor"
              " variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG3D_Free_structures:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  if ( !MMG3D_Free_names(MMG5_ARG_start,
                         MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                         MMG5_ARG_ppLs, ls,
                         MMG5_ARG_ppDisp, disp,
                         MMG5_ARG_ppSols, sols,
                         MMG5_ARG_end) ) {
    return 0;
  }

  /* mesh */
  assert(mesh && *mesh);


  MMG3D_Free_arrays(mesh,sol,ls,disp,sols);

  return 1;
}

/**
 * \param argptr list of the mmg structures for whose we want to deallocate the
 * name. Each structure must follow one of the \a MMG5_ARG* preprocessor
 * variable that allow to identify it.
 *
 * \a argptr contains at least a pointer toward a \a MMG5_pMesh structure
 * (that will contain the mesh and identified by the MMG5_ARG_ppMesh keyword)
 *
 *  To call the \a MMG3D_mmg3dlib function, you must also provide
 * a pointer toward a \a MMG5_pSol structure (that will contain the ouput
 * metric (and the input one, if provided) and identified by the MMG5_ARG_ppMet
 * keyword).
 *
 *  To call the \a MMG3D_mmg3dls function, you must also provide a pointer
 * toward a \a MMG5_pSol structure (that will contain the level-set function and
 * identified by the MMG5_ARG_ppLs keyword).
 *
 *  To call the \a MMG3D_mmg3dmov library, you must also provide a
 * pointer toward a \a MMG5_pSol structure storing the displacement (and
 * identified by the MMG5_ARG_ppDisp keyword).
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
int MMG3D_Free_names_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      psl,*sol,*disp,*ls,*sols;
  int            typArg,i;
  int            meshCount;

  meshCount = 0;
  disp = sol = ls = sols = NULL;

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
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    case(MMG5_ARG_ppSols):
      sols = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s: MMG3D_Free_names:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one of the following preprocessor"
              " variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s: MMG3D_Free_names:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n",__func__);
    return 0;
  }

  /* mesh & met */
  if (!sol)
    MMG5_mmgFree_names(*mesh,NULL);
  else
    MMG5_mmgFree_names(*mesh,*sol);

  /* disp */
  if ( disp && *disp ) {
    if ( (*disp)->namein ) {
      MMG5_DEL_MEM(*mesh,(*disp)->namein);
    }

    if ( (*disp)->nameout ) {
      MMG5_DEL_MEM(*mesh,(*disp)->nameout);
    }
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
