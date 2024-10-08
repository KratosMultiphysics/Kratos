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

#include "libmmg2d.h"
#include "libmmg2d_private.h"

mytime   MMG5_ctim[TIMEMAX];

/**
 * Print elapsed time at end of process.
 */
static void MMG5_endcod(void) {
  char   stim[32];

  chrono(OFF,&MMG5_ctim[0]);
  printim(MMG5_ctim[0].gdif,stim);
  fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdyRefs pointer toward the list of the boundary references.
 * \return npar, the number of local parameters at edges if success,
 * 0 otherwise.
 *
 * Count the local default values at edges and fill the list of the boundary
 * references.
 *
 */
static inline
int MMG2D_countLocalParamAtEdg( MMG5_pMesh mesh,MMG5_iNode **bdyRefs) {
  int         npar,ier;
  MMG5_int    k;

  /** Count the number of different boundary references and list it */
  (*bdyRefs) = NULL;

  k = mesh->na? mesh->edge[1].ref : 0;

  /* Try to alloc the first node */
  ier = MMG5_Add_inode( mesh, bdyRefs, k );
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate the first boundary"
           " reference node.\n",__func__);
    return 0;
  }
  else {
    assert(ier);
    npar = 1;
  }

  for ( k=1; k<=mesh->na; ++k ) {
    ier = MMG5_Add_inode( mesh, bdyRefs, mesh->edge[k].ref );

    if ( ier < 0 ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to list the edge references.\n"
              "              Uncomplete parameters file.\n",__func__);
      break;
    }
    else if ( ier ) ++npar;
  }

  return npar;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \param out pointer toward the file in which to write.
 * \return 1 if success, 0 otherwise.
 *
 * Write the local default values at edges in the parameter file.
 *
 */
static inline
int MMG2D_writeLocalParamAtEdg( MMG5_pMesh mesh, MMG5_iNode *bdryRefs,
                                 FILE *out ) {
  MMG5_iNode *cur;

  cur = bdryRefs;
  while( cur ) {
    fprintf(out,"%"MMG5_PRId" Edge %e %e %e \n",cur->val,
            mesh->info.hmin, mesh->info.hmax,mesh->info.hausd);
    cur = cur->nxt;
  }

  MMG5_Free_ilinkedList(mesh,bdryRefs);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 otherwise.
 *
 * Write a DEFAULT.mmg2d file containing the default values of parameters that
 * can be locally defined.
 *
 */
static inline
int MMG2D_writeLocalParam( MMG5_pMesh mesh ) {
  MMG5_iNode  *edgRefs,*triRefs;
  int          nparEdg,nparTri;
  char         *ptr,data[MMG5_FILESTR_LGTH];
  FILE         *out;

  /** Save the local parameters file */
  strcpy(data,mesh->namein);

  ptr = MMG5_Get_filenameExt(data);

  if ( ptr ) *ptr = '\0';
  strcat(data,".mmg2d");

  if ( !(out = fopen(data,"wb")) ) {
    fprintf(stderr,"\n  ** UNABLE TO OPEN %s.\n",data);
    return 0;
  }

  fprintf(stdout,"\n  %%%% %s OPENED\n",data);

  nparEdg = MMG2D_countLocalParamAtEdg( mesh, &edgRefs);
  if ( !nparEdg ) {
    fclose(out);
    return 0;
  }

  nparTri = MMG5_countLocalParamAtTri( mesh, &triRefs);
  if ( !nparTri ) {
    fclose(out);
    return 0;
  }

  fprintf(out,"parameters\n %d\n",nparTri+nparEdg);

  /** Write local param at triangles */
  if (! MMG2D_writeLocalParamAtEdg(mesh,edgRefs,out) ) {
    fclose(out);
    return 0;
  }

  /** Write local param at tetra */
  if (! MMG5_writeLocalParamAtTri(mesh,triRefs,out) ) {
    fclose(out);
    return 0;
  }

  fclose(out);
  fprintf(stdout,"  -- WRITING COMPLETED\n");

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric).
 * \param sol pointer toward a sol structure (metric).
 *
 * \return \ref MMG5_SUCCESS if success, \ref MMG5_LOWFAILURE if failed
 * but a conform mesh is saved and \ref MMG5_STRONGFAILURE if failed and we
 * can't save the mesh.
 *
 * Program to save the local default parameter file: read the mesh and metric
 * (needed to compite the hmax/hmin parameters), scale the mesh and compute the
 * hmax/hmin param, unscale the mesh and write the default parameter file.
 *
 */
static inline
int MMG2D_defaultOption(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  mytime    ctim[TIMEMAX];
  double    hsiz;
  char      stim[32];

  signal(SIGABRT,MMG2D_excfun);
  signal(SIGFPE,MMG2D_excfun);
  signal(SIGILL,MMG2D_excfun);
  signal(SIGSEGV,MMG2D_excfun);
  signal(SIGTERM,MMG2D_excfun);
  signal(SIGINT,MMG2D_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  if ( mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: "
            "unable to save of a local parameter file with"
            " the default parameters values because local parameters"
            " are provided.\n",__func__);
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
  }


  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  if ( met && met->np && (met->np != mesh->np) ) {
    fprintf(stderr,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,met->m);
    met->np = 0;
  }

  if ( sol && sol->np && (sol->np != mesh->np) ) {
    fprintf(stderr,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,sol->m);
    sol->np = 0;
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG2D_setfunc(mesh,met);
  MMG2D_Set_commonFunc();

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- DEFAULT PARAMETERS COMPUTATION\n");
  }

  /* scaling mesh and hmin/hmax computation*/
  if ( !MMG5_scaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  /* specific meshing + update hmin/hmax */
  if ( mesh->info.optim ) {
    if ( !MMG2D_doSol(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,sol) )
        _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }
  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) ) {
     if ( !MMG5_unscaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
     _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  /* unscaling mesh */
  if ( !MMG5_unscaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  /* Save the local parameters file */
  mesh->mark = 0;
  if ( !MMG2D_writeLocalParam(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to save the local parameters file.\n"
            "            Exit program.\n",__func__);
     _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
  }

  _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);
}

int main(int argc,char *argv[]) {
  MMG5_pMesh    mesh;
  MMG5_pSol     sol,met,disp,ls;
  int           ier,ierSave,fmtin,fmtout;
  char          stim[32],*ptr;

  /* Select line buffering even if the output is not a terminal and force stderr
   * and stdout to print in the same order as the events */
  setvbuf(stdout, NULL, _IOLBF, 1024);
  setvbuf(stderr, NULL, _IOLBF, 1024);

  /* Version info */
#ifndef MMG_COMPARABLE_OUTPUT
  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MMG_VERSION_RELEASE,MMG_RELEASE_DATE);
  fprintf(stdout,"     %s\n",MMG_COPYRIGHT);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
#endif

  /* Print timer at exit */
  atexit(MMG5_endcod);

  MMG2D_Set_commonFunc();
  tminit(MMG5_ctim,TIMEMAX);
  chrono(ON,&MMG5_ctim[0]);

  /* assign default values */
  mesh = NULL;
  met  = NULL;
  ls   = NULL;
  disp = NULL;

  if ( !MMG2D_Init_mesh(MMG5_ARG_start,
                        MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                        MMG5_ARG_ppLs,&ls,
                        MMG5_ARG_ppDisp,&disp,
                        MMG5_ARG_end) )
    return MMG5_STRONGFAILURE;

  /* reset default values for file names */
  if ( !MMG2D_Free_names(MMG5_ARG_start,
                         MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                         MMG5_ARG_ppLs,&ls,
                         MMG5_ARG_ppDisp,&disp,
                         MMG5_ARG_end) )
    return MMG5_STRONGFAILURE;

  /* Set default metric size */
  if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Scalar) )
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);

  /* Read command line */
  if ( !MMG2D_parsar(argc,argv,mesh,met,ls) )  return MMG5_STRONGFAILURE;

  /* load data */
  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);

  /* For each mode: pointer over the solution structure to load */
  if ( mesh->info.lag >= 0 ) {
    sol = disp;
  }
  else if ( mesh->info.iso || mesh->info.isosurf ) {
    sol = ls;
  }
  else {
    sol = met;
  }

  /* read mesh/sol files */
  ptr   = MMG5_Get_filenameExt(mesh->namein);
  fmtin = MMG5_Get_format(ptr,MMG5_FMT_MeditASCII);

  switch ( fmtin ) {
  case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    ier = MMG2D_loadMshMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtp ):
    ier = MMG2D_loadVtpMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtu ):
    ier = MMG2D_loadVtuMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtk ):
    ier = MMG2D_loadVtkMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_MeditASCII ): case ( MMG5_FMT_MeditBinary ):
    ier = MMG2D_loadMesh(mesh,mesh->namein);
    if ( ier <  1 ) { break; }

    /* Read displacement in lag mode */
    if ( mesh->info.lag >= 0 ) {
      /* In Lagrangian mode, the name of the displacement file has been parsed in ls */
      if ( !MMG2D_Set_inputSolName(mesh,disp,ls->namein) ) {
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
      MMG5_DEL_MEM(mesh,ls->namein);
    }

    if ( mesh->info.lag >= 0 || mesh->info.iso || mesh->info.isosurf ) {
      /* displacement or isovalue are mandatory */
      if (  MMG2D_loadSol(mesh,sol,sol->namein) < 1 ) {
        /* displacement or isovalue are mandatory */
        fprintf(stdout,"  ## ERROR: UNABLE TO LOAD SOLUTION.\n");
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
    }
    else {
      /* Facultative metric */
      if ( MMG2D_loadSol(mesh,met,met->namein) == -1 ) {
        fprintf(stdout,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
    }
    /* In iso mode: read metric if any */
    if ( ( mesh->info.iso || mesh->info.isosurf ) && met->namein ) {
      if (  MMG2D_loadSol(mesh,met,met->namein) < 1 ) {
        fprintf(stdout,"  ## ERROR: UNABLE TO LOAD METRIC.\n");
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
    }
    break;
  default:
    fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTED.\n",MMG5_Get_formatName(fmtin) );
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
  }

  if ( ier < 1) {
    if ( ier==0 ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",mesh->namein);
      fprintf(stderr,"  ** UNABLE TO OPEN INPUT FILE.\n");
    }
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
  }

  /* Check input data */
  if ( mesh->info.lag >= 0 ) {
    if ( met->namein ) {
      fprintf(stdout,"  ## WARNING: MESH ADAPTATION UNAVAILABLE IN"
              " LAGRANGIAN MODE. METRIC IGNORED.\n");
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
    }
  }
  else if ( mesh->info.iso || mesh->info.isosurf ) {
    if ( ls->m == NULL ) {
      fprintf(stderr,"\n  ## ERROR: NO ISOVALUE DATA.\n");
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
    }
  }

  /* Read parameter file */
  if ( !MMG2D_parsop(mesh,met) )
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);

  chrono(OFF,&MMG5_ctim[1]);
  if ( mesh->info.imprim >= 0 ) {
    printim(MMG5_ctim[1].gdif,stim);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
  }

  if ( mesh->mark ) {
    /* Save a local parameters file containing the default parameters */
    ier = MMG2D_defaultOption(mesh,met,disp);
    MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,ier);
  }
  /* Lagrangian mode */
  else if ( mesh->info.lag > -1 ) {
    ier = MMG2D_mmg2dmov(mesh,met,disp);
  }
  /* Level Set mode */
  else if ( mesh->info.iso || mesh->info.isosurf ) {
    ier = MMG2D_mmg2dls(mesh,ls,met);
  }
  /* Mesh generation mode */
  else if ( !mesh->nt ) {
    ier = MMG2D_mmg2dmesh(mesh,met);
  }
  /* Remeshing mode */
  else {
    if ( met && ls && met->namein && ls->namein ) {
      fprintf(stdout,"\n  ## ERROR: IMPOSSIBLE TO PROVIDE BOTH A METRIC"
              " AND A SOLUTION IN ADAPTATION MODE.\n");
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
    }
    
    ier = MMG2D_mmg2dlib(mesh,met);
  }

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&MMG5_ctim[1]);
    if ( mesh->info.imprim > 0 )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);

    ptr    = MMG5_Get_filenameExt(mesh->nameout);
    fmtout = MMG5_Get_format(ptr,fmtin);

    switch ( fmtout ) {
    case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
      ierSave = MMG2D_saveMshMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtp ):
      ierSave = MMG2D_saveVtpMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtu ):
      ierSave = MMG2D_saveVtuMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtk ):
      ierSave = MMG2D_saveVtkMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_Tetgen ):
      ierSave = MMG2D_saveTetgenMesh(mesh,mesh->nameout);
      /* This format dont allow to save a solution: use a .sol file */
      if ( !ierSave ) {
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
      if ( met && met->np ) {
        ierSave = MMG2D_saveSol(mesh,met,mesh->nameout);
      }
      break;
    default:
      ierSave = MMG2D_saveMesh(mesh,mesh->nameout);
      if ( !ierSave ) {
        MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      }
      if ( met && met->np ) {
        ierSave = MMG2D_saveSol(mesh,met,mesh->nameout);
      }
      break;
    }

    if ( !ierSave )
      MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);

    chrono(OFF,&MMG5_ctim[1]);
    if ( mesh->info.imprim > 0 ) fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* free mem */
  MMG2D_RETURN_AND_FREE(mesh,met,ls,disp,ier);
}
