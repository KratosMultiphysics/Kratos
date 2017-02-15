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
 * \file mmgs/mmgs.c
 * \brief Main file for MMGS executable: perform surface mesh adaptation.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgs.h"

#include <math.h>

mytime         MMG5_ctim[TIMEMAX];


/**
 * Print elapsed time at end of process.
 */
static void _MMG5_endcod() {
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
static int _MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  float      fp1,fp2,hausd;
  int        ref,i,j,ret,npar;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  strcpy(data,mesh->namein);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mmgs");
  in = fopen(data,"rb");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.mmgs");
    in = fopen(data,"rb");
    if ( !in )  return(1);
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* read parameters */
  mesh->info.npar = 0;
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"parameters") ) {
      fscanf(in,"%d",&npar);
      if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_numberOfLocalParam,npar) )
        exit(EXIT_FAILURE);

      for (i=0; i<mesh->info.npar; i++) {
        fscanf(in,"%d %s ",&ref,buf);
        ret = fscanf(in,"%f %f %f",&fp1,&fp2,&hausd);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);

        if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") ) {
          if ( !MMGS_Set_localParameter(mesh,met,MMG5_Triangle,ref,fp1,fp2,hausd) ) {
            exit(EXIT_FAILURE);
          }
        }
        /* else if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") ) { */
        /*   if ( !MMGS_Set_localParameter(mesh,met,MMG5_Vertex,ref,fp1,fp2,hausd) ) { */
        /*     exit(EXIT_FAILURE); */
        /*   } */
        /* } */
        else {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
      }
    }
  }
  fclose(in);
  return(1);
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
int _MMGS_writeLocalParam( MMG5_pMesh mesh ) {
  _MMG5_iNode  *triRefs;
  int          npar;
  char         *ptr,data[128];
  FILE         *out;

  strcpy(data,mesh->namein);
  ptr = strstr(data,".mesh");
  if ( ptr ) *ptr = '\0';
  strcat(data,".mmgs");

  /** Save the local parameters file */
  if ( !(out = fopen(data,"wb")) ) {
    fprintf(stderr,"\n  ** UNABLE TO OPEN %s.\n",data);
    return(0);
  }

  fprintf(stdout,"\n  %%%% %s OPENED\n",data);


  npar = _MMG5_countLocalParamAtTri( mesh, &triRefs);

  if ( !npar ) return 0;

  fprintf(out,"parameters\n %d\n",npar);

  if ( !_MMG5_writeLocalParamAtTri(mesh, triRefs, out) ) return 0;

  fclose(out);
  fprintf(stdout,"  -- WRITING COMPLETED\n");

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward a sol structure (metric).
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
int _MMGS_defaultOption(MMG5_pMesh mesh,MMG5_pSol met) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  _MMGS_Set_commonFunc();

  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  if ( mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: "
            "unable to save of a local parameter file with"
            " the default parameters values because local parameters"
            " are provided.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
  }


  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stderr,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMGS_setfunc(mesh,met);

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMGS: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- DEFAULT PARAMETERS COMPUTATION\n");
  }

  /* scaling mesh and hmin/hmax computation*/
  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* unscaling mesh */
  if ( !_MMG5_unscaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* Save the local parameters file */
  mesh->mark = 0;
  if ( !_MMGS_writeLocalParam(mesh) ) {
    fprintf(stderr,"  ## Error: Unable to save the local parameters file.\n"
            "            Exit program.\n");
     _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);
  }

  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}


int main(int argc,char *argv[]) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier,ierSave,msh;
  char       stim[32];

  fprintf(stdout,"  -- MMGS, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  _MMGS_Set_commonFunc();

  /* trap exceptions */
  atexit(_MMG5_endcod);

  tminit(MMG5_ctim,TIMEMAX);
  chrono(ON,&MMG5_ctim[0]);

  /* assign default values */
  mesh = NULL;
  met  = NULL;

  MMGS_Init_mesh(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                 MMG5_ARG_end);

  /* reset default values for file names */
  MMGS_Free_names(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                  MMG5_ARG_end);

  /* command line */
  if ( !MMGS_parsar(argc,argv,mesh,met) )  return(MMG5_STRONGFAILURE);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);

  /* read mesh file */
  msh = 0;
  ier = MMGS_loadMesh(mesh,mesh->namein);
  if ( !ier ) {
    ier = MMGS_loadMshMesh(mesh,met,mesh->namein);
    msh = 1;
  }

  if ( !msh ) {
    ier = MMGS_loadSol(mesh,met,met->namein);
    if ( ier==-1 ) {
      fprintf(stderr,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      _MMGS_RETURN_AND_FREE(mesh,met,MMG5_STRONGFAILURE);
    }
  }

  if ( !_MMG5_parsop(mesh,met) )
    _MMGS_RETURN_AND_FREE(mesh,met,MMG5_LOWFAILURE);

  chrono(OFF,&MMG5_ctim[1]);
  printim(MMG5_ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  if ( mesh->mark ) {
    /* Save a local parameters file containing the default parameters */
    ier = _MMGS_defaultOption(mesh,met);
    _MMGS_RETURN_AND_FREE(mesh,met,ier);
  }
  else if ( mesh->info.iso ) {
    ier = MMGS_mmgsls(mesh,met);
  }
  else {
    ier = MMGS_mmgslib(mesh,met);
  }

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&MMG5_ctim[1]);
    if ( mesh->info.imprim )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);

    if ( !strcmp(&mesh->nameout[strlen(mesh->nameout)-5],".mesh") ||
         !strcmp(&mesh->nameout[strlen(mesh->nameout)-6],".meshb") )
      msh = 0;

    else if (!strcmp(&mesh->nameout[strlen(mesh->nameout)-4],".msh") ||
             !strcmp(&mesh->nameout[strlen(mesh->nameout)-5],".mshb") )
      msh = 1;

    if ( !msh )
      ierSave = MMGS_saveMesh(mesh,mesh->nameout);
    else
      ierSave = MMGS_saveMshMesh(mesh,met,mesh->nameout);

    if ( !ierSave )
      _MMGS_RETURN_AND_FREE(mesh,met,MMG5_STRONGFAILURE);

    if ( !msh && !MMGS_saveSol(mesh,met,met->nameout) )
      _MMGS_RETURN_AND_FREE(mesh,met,MMG5_STRONGFAILURE);

    chrono(OFF,&MMG5_ctim[1]);
    if ( mesh->info.imprim )  fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* release memory */
  /* free mem */
  _MMGS_RETURN_AND_FREE(mesh,met,ier);

  return(0);
}
