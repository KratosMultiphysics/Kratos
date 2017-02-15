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

#include "mmg3d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward a sol structure (metric or level-set).
 * \param disp pointer toward a sol structure (displacement).
 *
 * Allocate the mesh and solutions structures at \a MMG3D format.
 *
 */
static inline
void _MMG3D_Alloc_mesh(MMG5_pMesh *mesh, MMG5_pSol *sol, MMG5_pSol *disp
  ) {

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

  /* displacement allocation */
  if ( disp ) {
    if ( *disp )
      _MMG5_DEL_MEM(*mesh,*disp,sizeof(MMG5_Sol));
    _MMG5_SAFE_CALLOC(*disp,1,MMG5_Sol);
  }

  return;
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
void _MMG3D_Init_woalloc_mesh(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSol disp
  ) {

  _MMG3D_Set_commonFunc();

  (mesh)->dim  = 3;
  (mesh)->ver  = 2;
  (sol)->dim   = 3;
  (sol)->ver   = 2;
  (sol)->size  = 1;
  if ( disp ) {
    (disp)->dim  = 3;
    (disp)->ver  = 2;
    (disp)->size = 2;
  }

  /* Default parameters values */
  MMG3D_Init_parameters(mesh);

  /* Default vaules for file names */
  MMG3D_Init_fileNames(mesh,sol);

  if ( disp ) {
    MMG3D_Set_inputSolName(mesh,disp,"");
    MMG3D_Set_outputSolName(mesh,disp,"");
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
 * Internal function for structure allocations (taking a va_list argument).
 *
 */
void _MMG3D_Init_mesh_var( va_list argptr ) {
  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  disp = sol = NULL;


  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
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
      fprintf(stderr,"  ## Error: MMG3D_Init_mesh:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one"
              " of the MMG5_ARG* preprocessor variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              "  MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMG3D_Init_mesh:\n"
            " you need to initialize the mesh structure that"
            " will contain your mesh.\n");
    exit(EXIT_FAILURE);
  }

  if ( !sol ) {
    fprintf(stderr,"  ## Error: MMG3D_Init_mesh:\n"
            " you need to initialize a solution structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or"
            " MMG5_ARG_ppLs preprocessor variable) that will contain the output"
            " mesh metric"
            " informations, and the input one, if provided.\n.");
    exit(EXIT_FAILURE);
  }

  /* allocations */
  _MMG3D_Alloc_mesh(mesh,sol,disp);

  /* initialisations */
  if ( disp )
    _MMG3D_Init_woalloc_mesh(*mesh,*sol,*disp);
  else
    _MMG3D_Init_woalloc_mesh(*mesh,*sol,NULL);

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
 * Internal function for deallocations before return (taking a va_list as
 * argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 */
void _MMG3D_Free_all_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  disp = sol = NULL;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
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
    default:
      fprintf(stderr,"  ## Error: MMG3D_Free_all:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one of the following preprocessor"
              " variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMG3D_Free_all:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n");
    exit(EXIT_FAILURE);
  }

  if ( !sol ) {
    fprintf(stderr,"  ## Error: MMG3D_Free_all:\n"
            " you need to provide your metric structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or"
            " MMG5_ARG_ppLs preprocessor variable)"
            " to allow to free the associated memory.\n");
  }


  if ( !disp )
    MMG3D_Free_structures(MMG5_ARG_start,
                          MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                          MMG5_ARG_end);
  else
    MMG3D_Free_structures(MMG5_ARG_start,
                          MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                          MMG5_ARG_ppDisp, disp,
                          MMG5_ARG_end);

  _MMG5_SAFE_FREE(*mesh);

  if ( sol )
    _MMG5_SAFE_FREE(*sol);

  if ( disp )
    _MMG5_SAFE_FREE(*disp);

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
 * Internal function for structures deallocations before return (taking a
 * va_list as argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void _MMG3D_Free_structures_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  disp = sol = NULL;

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
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"  ## Error: MMG3D_Free_structures:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one of the following preprocessor"
              " variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMG3D_Free_structures:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n");
    exit(EXIT_FAILURE);
  }

  if ( !disp )
    MMG3D_Free_names(MMG5_ARG_start,
                     MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                     MMG5_ARG_end);
  else
    MMG3D_Free_names(MMG5_ARG_start,
                     MMG5_ARG_ppMesh, mesh, MMG5_ARG_ppMet, sol,
                     MMG5_ARG_ppDisp, disp,
                     MMG5_ARG_end);

 /* mesh */
  assert(mesh && *mesh);
  if ( (*mesh)->point )
    _MMG5_DEL_MEM((*mesh),(*mesh)->point,((*mesh)->npmax+1)*sizeof(MMG5_Point));

  if ( (*mesh)->tetra )
    _MMG5_DEL_MEM((*mesh),(*mesh)->tetra,((*mesh)->nemax+1)*sizeof(MMG5_Tetra));

  if ( (*mesh)->prism )
    _MMG5_DEL_MEM((*mesh),(*mesh)->prism,((*mesh)->nprism+1)*sizeof(MMG5_Prism));

  if ( (*mesh)->edge )
    _MMG5_DEL_MEM((*mesh),(*mesh)->edge,((*mesh)->na+1)*sizeof(MMG5_Edge));

  if ( (*mesh)->adjt )
    _MMG5_DEL_MEM(*mesh,(*mesh)->adjt,(3*(*mesh)->nt+4)*sizeof(int));

  if ( (*mesh)->adja )
    _MMG5_DEL_MEM((*mesh),(*mesh)->adja,(4*(*mesh)->nemax+5)*sizeof(int));

  if ( (*mesh)->adjapr )
    _MMG5_DEL_MEM((*mesh),(*mesh)->adjapr,(5*(*mesh)->nprism+6)*sizeof(int));

  if ( (*mesh)->xpoint )
    _MMG5_DEL_MEM((*mesh),(*mesh)->xpoint,((*mesh)->xpmax+1)*sizeof(MMG5_xPoint));

  if ( (*mesh)->htab.geom )
    _MMG5_DEL_MEM((*mesh),(*mesh)->htab.geom,((*mesh)->htab.max+1)*sizeof(MMG5_hgeom));

  if ( (*mesh)->tria )
    _MMG5_DEL_MEM((*mesh),(*mesh)->tria,((*mesh)->nt+1)*sizeof(MMG5_Tria));

 if ( (*mesh)->quad )
    _MMG5_DEL_MEM((*mesh),(*mesh)->quad,((*mesh)->nquad+1)*sizeof(MMG5_Quad));

  if ( (*mesh)->xtetra )
    _MMG5_DEL_MEM((*mesh),(*mesh)->xtetra,((*mesh)->xtmax+1)*sizeof(MMG5_xTetra));

  if ( (*mesh)->xprism )
    _MMG5_DEL_MEM((*mesh),(*mesh)->xprism,((*mesh)->xpr+1)*sizeof(MMG5_xPrism));

  /* sol */
  if ( sol && (*sol) && (*sol)->m )
    _MMG5_DEL_MEM((*mesh),(*sol)->m,((*sol)->size*((*sol)->npmax+1))*sizeof(double));

  /* disp */
  if ( disp && (*disp) && (*disp)->m )
    _MMG5_DEL_MEM((*mesh),(*disp)->m,((*disp)->size*((*disp)->npmax+1))*sizeof(double));

  /* (*mesh)->info */
  if ( (*mesh)->info.npar && (*mesh)->info.par )
    _MMG5_DEL_MEM((*mesh),(*mesh)->info.par,(*mesh)->info.npar*sizeof(MMG5_Par));

  if ( (*mesh)->info.imprim>5 || (*mesh)->info.ddebug )
    printf("  MEMORY USED AT END (bytes) %ld\n",_MMG5_safeLL2LCast((*mesh)->memCur));

  return;
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
 * Internal function for name deallocations before return (taking a va_list as
 * argument).
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 */
void _MMG3D_Free_names_var(va_list argptr)
{

  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg;
  int            meshCount;

  meshCount = 0;
  disp = sol = NULL;

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
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      break;
    default:
      fprintf(stderr,"  ## Error: MMG3D_Free_names:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one of the following preprocessor"
              " variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              " MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      exit(EXIT_FAILURE);
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"  ## Error: MMG3D_Free_names:\n"
            " you need to provide your mesh structure"
            " to allow to free the associated memory.\n");
    exit(EXIT_FAILURE);
  }
  if ( !sol ) {
    fprintf(stderr,"  ## Error: MMG3D_Free_names:\n"
            " you need to provide your metric structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or"
            " MMG5_ARG_ppLs preprocessor variable)"
            " to allow to free the associated memory.\n");
  }

  /* mesh & met */
  MMG5_mmgFree_names(*mesh,*sol);

  /* disp */
  if ( disp && *disp ) {
    if ( (*disp)->namein ) {
      _MMG5_DEL_MEM(*mesh,(*disp)->namein,(strlen((*disp)->namein)+1)*sizeof(char));
    }

    if ( (*disp)->nameout ) {
      _MMG5_DEL_MEM(*mesh,(*disp)->nameout,(strlen((*disp)->nameout)+1)*sizeof(char));
    }
  }

  return;
}
