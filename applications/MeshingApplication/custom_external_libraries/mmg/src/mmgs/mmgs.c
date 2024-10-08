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
 * \file mmgs/mmgs.c
 * \brief Main file for MMGS executable: perform surface mesh adaptation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmgs.h"
#include "libmmgs_private.h"
#include <math.h>

mytime         MMG5_ctim[TIMEMAX];


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
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Read local parameters file. This file must have the same name as
 * the mesh with the \a .mmgs extension or must be named \a
 * DEFAULT.mmgs.
 *
 */
static int MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  float      fp1,fp2,hausd;
  int        i,j,ret,npar,nbr,split;
  MMG5_int   ref,rin,rex,br;
  char       *ptr,buf[256],data[256];
  FILE       *in;
  fpos_t     position;

  /* check for parameter file */
  strcpy(data,mesh->namein);

  ptr = MMG5_Get_filenameExt(data);

  if ( ptr )  *ptr = '\0';
  strcat(data,".mmgs");

  in = fopen(data,"rb");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.mmgs");
    in = fopen(data,"rb");
    if ( !in ) {
      return 1;
    }
  }
  if ( mesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n  %%%% %s OPENED\n",data);
  }

  /* read parameters */
  mesh->info.npar = 0;
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%255s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; (size_t)i<strlen(data); i++) data[i] = tolower(data[i]);

    /* Read user-defined references for the LS mode */
    if ( !strcmp(data,"lsreferences") ) {
      ret = fscanf(in,"%d",&npar);
      if ( !ret ) {
        fprintf(stderr,"  %%%% Wrong format for lsreferences: %d\n",npar);
        return 0;
      }

      if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_numberOfMat,npar) ) {
        return 0;
      }
      for (i=0; i<mesh->info.nmat; i++) {
        MMG_FSCANF(in,"%" MMG5_PRId "",&ref);
        fgetpos(in,&position);
        MMG_FSCANF(in,"%255s",data);
        split = MMG5_MMAT_NoSplit;
        rin = rex = ref;
        if ( strcmp(data,"nosplit") ) {
          fsetpos(in,&position);
          split = MMG5_MMAT_Split;
          MMG_FSCANF(in,"%" MMG5_PRId "",&rin);
          MMG_FSCANF(in,"%" MMG5_PRId "",&rex);
        }
        if ( !MMGS_Set_multiMat(mesh,met,ref,split,rin,rex) ) {
          return 0;
        }
      }
    }
    /* Read user-defined local parameters and store them in the structure info->par */
    else if ( !strcmp(data,"parameters") ) {
      MMG_FSCANF(in,"%d",&npar);

      if ( !ret ) {
        fprintf(stderr,"  %%%% Wrong format for parameters: %d\n",npar);
        return 0;
      }
      else if ( npar > MMG5_LPARMAX ) {
        fprintf(stderr,"  %%%% Too many local parameters %d. Abort\n",npar);
        return 0;
      }

      /* Allocate memory and fill the info->par table (adding one, corresponding to the command line data) */
      if ( npar ) {
        if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_numberOfLocalParam,npar) )
          return 0;

        for (i=0; i<mesh->info.npar; i++) {
          MMG_FSCANF(in,"%" MMG5_PRId " %255s ",&ref,buf);
          ret = fscanf(in,"%f %f %f",&fp1,&fp2,&hausd);

          if ( !ret ) {
            fprintf(stderr,"  %%%% Wrong format: %s\n",buf);
            return 0;
          }

          for (j=0; (size_t)j<strlen(buf); j++)  buf[j] = tolower(buf[j]);

          if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") ) {
            if ( !MMGS_Set_localParameter(mesh,met,MMG5_Triangle,ref,fp1,fp2,hausd) ) {
              return 0;
            }
          }
          else {
            fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
            return 0;
          }
        }
      }
    }
    /* Read user-defined references where connected components should stay attached in ls mode */
    else if ( !strcmp(data,"lsbasereferences") ) {
      MMG_FSCANF(in,"%d",&nbr);
      if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_numberOfLSBaseReferences,nbr) )
        return 0;

      for (i=0; i<mesh->info.nbr; i++) {
        MMG_FSCANF(in,"%" MMG5_PRId "",&br);
        if ( !MMGS_Set_lsBaseReference(mesh,met,br) ) {
          return 0;
        }
      }
    }
    else {
      fprintf(stderr,"  %%%% Wrong format: %s\n",data);
      return 0;
    }
  }

  fclose(in);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 otherwise.
 *
 * Write a DEFAULT.mmg3d file containing the default values of parameters that
 * can be locally defined.
 *
 */
static inline
int MMGS_writeLocalParam( MMG5_pMesh mesh ) {
  MMG5_iNode  *triRefs;
  int          npar;
  char         *ptr,data[MMG5_FILESTR_LGTH];
  FILE         *out;

  strcpy(data,mesh->namein);

  ptr = MMG5_Get_filenameExt(data);

  if ( ptr ) *ptr = '\0';

  strcat(data,".mmgs");

  /** Save the local parameters file */
  if ( !(out = fopen(data,"wb")) ) {
    fprintf(stderr,"\n  ** UNABLE TO OPEN %s.\n",data);
    return 0;
  }

  fprintf(stdout,"\n  %%%% %s OPENED\n",data);


  npar = MMG5_countLocalParamAtTri( mesh, &triRefs);

  if ( !npar ) {
    fclose(out);
    return 0;
  }

  fprintf(out,"parameters\n %d\n",npar);

  if ( !MMG5_writeLocalParamAtTri(mesh, triRefs, out) ) {
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
int MMGS_defaultOption(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  mytime    ctim[TIMEMAX];
  double    hsiz;
  char      stim[32];

  MMGS_Set_commonFunc();

  signal(SIGABRT,MMG5_excfun);
  signal(SIGFPE,MMG5_excfun);
  signal(SIGILL,MMG5_excfun);
  signal(SIGSEGV,MMG5_excfun);
  signal(SIGTERM,MMG5_excfun);
  signal(SIGINT,MMG5_excfun);

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
  MMGS_setfunc(mesh,met);

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- DEFAULT PARAMETERS COMPUTATION\n");
  }

  /* scaling mesh and hmin/hmax computation*/
  if ( !MMG5_scaleMesh(mesh,met,sol) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  /* Specific meshing + hmin/hmax update */
  if ( mesh->info.optim ) {
    if ( !MMGS_doSol(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,sol) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
        _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) ) {
     if ( !MMG5_unscaleMesh(mesh,met,sol) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
       _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  /* unscaling mesh */
  if ( !MMG5_unscaleMesh(mesh,met,sol) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  /* Save the local parameters file */
  mesh->mark = 0;
  if ( !MMGS_writeLocalParam(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to save the local parameters file.\n"
            "            Exit program.\n",__func__);
       _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
  }

    _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);
}


int main(int argc,char *argv[]) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol,met,ls;
  int        ier,ierSave,fmtin,fmtout;
  char       stim[32],*ptr;

  /* Select line buffering even if the output is not a terminal and force stderr
   * and stdout to print in the same order as the events */
  setvbuf(stdout, NULL, _IOLBF, 1024);
  setvbuf(stderr, NULL, _IOLBF, 1024);

  /* Version info */
#ifndef MMG_COMPARABLE_OUTPUT
  fprintf(stdout,"  -- MMGS, Release %s (%s) \n",MMG_VERSION_RELEASE,MMG_RELEASE_DATE);
  fprintf(stdout,"     %s\n",MMG_COPYRIGHT);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
#endif

  MMGS_Set_commonFunc();

  /* Print timer at exit */
  atexit(MMG5_endcod);

  tminit(MMG5_ctim,TIMEMAX);
  chrono(ON,&MMG5_ctim[0]);

  /* assign default values */
  mesh = NULL;
  met  = NULL;
  ls   = NULL;

  MMGS_Init_mesh(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                 MMG5_ARG_ppLs,&ls,
                 MMG5_ARG_end);

  /* reset default values for file names */
  MMGS_Free_names(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                  MMG5_ARG_ppLs,&ls,
                  MMG5_ARG_end);

  /* command line */
  if ( !MMGS_parsar(argc,argv,mesh,met,ls) )  return MMG5_STRONGFAILURE;

  /* load data */
  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);

  /* For each mode: pointer over the solution structure to load */
  if ( mesh->info.iso ) {
    sol = ls;
  }
  else {
    sol = met;
  }

  /* read mesh file */
  ptr   = MMG5_Get_filenameExt(mesh->namein);
  fmtin = MMG5_Get_format(ptr,MMG5_FMT_MeditASCII);

  switch ( fmtin ) {

  case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    ier = MMGS_loadMshMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtp ):
    ier = MMGS_loadVtpMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtu ):
    ier = MMGS_loadVtuMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_VtkVtk ):
    ier = MMGS_loadVtkMesh(mesh,sol,mesh->namein);
    break;

  case ( MMG5_FMT_MeditASCII ): case ( MMG5_FMT_MeditBinary ):
    ier = MMGS_loadMesh(mesh,mesh->namein);
    if ( ier <  1 ) { break; }

    /* read level-set in iso mode */
    if ( mesh->info.iso || mesh->info.isosurf ) {
      if ( MMGS_loadSol(mesh,ls,ls->namein) < 1 ) {
        fprintf(stderr,"\n  ## ERROR: UNABLE TO LOAD LEVEL-SET.\n");
        MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
      }
      if ( met->namein ) {
        if ( MMGS_loadSol(mesh,met,met->namein) < 1 ) {
          fprintf(stdout,"  ## ERROR: UNABLE TO LOAD METRIC.\n");
          MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
        }
      }
    }
    else {
      /* read metric if any */
      if ( MMGS_loadSol(mesh,met,met->namein) == -1 ) {
        fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
      }
    }
    break;

  default:
    fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTD.\n",MMG5_Get_formatName(fmtin) );
    MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
  }

  if ( ier<1 ) {
    if ( ier==0 ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",mesh->namein);
      fprintf(stderr,"  ** UNABLE TO OPEN INPUT FILE.\n");
    }
    MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
  }

  /* Check input data */
  if ( mesh->info.iso || mesh->info.isosurf ) {
    if ( ls->m==NULL ) {
      fprintf(stderr,"\n  ## ERROR: NO ISOVALUE DATA.\n");
      MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
    }
  }

  /* Read parameter file */
  if ( !MMG5_parsop(mesh,met) )
    MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_LOWFAILURE);

  chrono(OFF,&MMG5_ctim[1]);

  if ( mesh->info.imprim >= 0 ) {
    printim(MMG5_ctim[1].gdif,stim);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
  }

  if ( mesh->mark ) {
    /* Save a local parameters file containing the default parameters */
    ier = MMGS_defaultOption(mesh,met,ls);
    MMGS_RETURN_AND_FREE(mesh,met,ls,ier);
  }
  else if ( mesh->info.iso || mesh->info.isosurf ) {
    ier = MMGS_mmgsls(mesh,ls,met);
  }
  else {
    if ( met && ls && met->namein && ls->namein ) {
      fprintf(stdout,"\n  ## ERROR: IMPOSSIBLE TO PROVIDE BOTH A METRIC"
              " AND A SOLUTION IN ADAPTATION MODE.\n");
      MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
    }

    ier = MMGS_mmgslib(mesh,met);
  }

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&MMG5_ctim[1]);
    if ( mesh->info.imprim > 0 )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);

    ptr    = MMG5_Get_filenameExt(mesh->nameout);
    fmtout = MMG5_Get_format(ptr,fmtin);

    switch ( fmtout ) {
    case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
      ierSave = MMGS_saveMshMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtp ):
      ierSave = MMGS_saveVtpMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtu ):
      ierSave = MMGS_saveVtuMesh(mesh,met,mesh->nameout);
      break;
    case ( MMG5_FMT_VtkVtk ):
      ierSave = MMGS_saveVtkMesh(mesh,met,mesh->nameout);
      break;
    default:
      ierSave = MMGS_saveMesh(mesh,mesh->nameout);
      if ( !ierSave ) {
        MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);
      }
      if ( met && met->np ) {
        ierSave = MMGS_saveSol(mesh,met,met->nameout);
      }
      break;
    }

    if ( !ierSave )
      MMGS_RETURN_AND_FREE(mesh,met,ls,MMG5_STRONGFAILURE);

    chrono(OFF,&MMG5_ctim[1]);
    if ( mesh->info.imprim > 0 )  fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* release memory */
  /* free mem */
  MMGS_RETURN_AND_FREE(mesh,met,ls,ier);

  return 0;
}
