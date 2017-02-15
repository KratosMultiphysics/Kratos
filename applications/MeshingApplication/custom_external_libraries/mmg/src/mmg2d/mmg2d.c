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
#include "mmg2d.h"

mytime   MMG5_ctim[TIMEMAX];

static void usage(char *name) {
  _MMG5_mmgUsage(name);

  fprintf(stdout,"-lag [0/1/2]    Lagrangian mesh displacement according to mode 0/1/2\n");
  fprintf(stdout,"-mov filedep    (with -lag option)\n");
  fprintf(stdout,"-nsd val        only if no given triangle, save the subdomain nb (0==all subdomain)\n");
  fprintf(stdout,"-msh val        read and write to gmsh visu if val = 1 (out) if val=2 (in and out)\n");
  fprintf(stdout,"-degrad Qw Qdeg (with -lag option) : threshold for optimization\n");

  /* fprintf(stdout,"-per          obsolete : to deal with periodic mesh on a square\n");*/

  fprintf(stdout,"\n");
  /*

    fprintf(stdout,"-optim       mesh optimization\n");
  */
  fprintf(stdout,"-nosurf       no surfacic modifications\n");

  fprintf(stdout,"-noinsert     no insertion/suppression point\n");
  fprintf(stdout,"-noswap       no edge flipping\n");
  fprintf(stdout,"-nomove       no point relocation\n");
  fprintf(stdout,"-bucket val   Specify the size of bucket per dimension \n");
  fprintf(stdout,"\n\n");

  exit(EXIT_FAILURE);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param qdegrad optimization thresholds.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
static inline void _MMG5_defaultValues(MMG5_pMesh mesh, double qdegrad[2]) {
  mesh->info.dhd = cos((135.*M_PI)/180.);

  _MMG5_mmgDefaultValues(mesh);

  fprintf(stdout,"Optimization threshold "
          "   (-degrad) : %e %e\n",qdegrad[0],qdegrad[1]);
  fprintf(stdout,"Bucket size per dimension (-bucket) : %d\n",
          mesh->info.octree);
  fprintf(stdout,"\n\n");

  exit(EXIT_FAILURE);
}

int parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,double *qdegrad) {
  int     i;
  // char   *ptr;
  char    namein[128];

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      _MMG5_defaultValues(mesh,qdegrad);
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        usage(argv[0]);
        break;
      case 'a':
        if ( !strcmp(argv[i],"-ar") && ++i < argc )
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_angleDetection,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        break;
      case 'b':
        if ( !strcmp(argv[i],"-bucket") && ++i < argc )
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_bucket,
                                     atoi(argv[i])) )
            exit(EXIT_FAILURE);
        break;
      case 'd':  /* debug */
        if ( !strcmp(argv[i],"-degrad") ) {
          ++i;
          qdegrad[0] = atof(argv[i++])/ALPHA;
          qdegrad[1] = atof(argv[i]);
        } else {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_debug,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hmin,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hmax,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hausd,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_hgrad,
                                     atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else
          usage(argv[0]);
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG2D_Set_inputMeshName(mesh, argv[i]) )
              exit(EXIT_FAILURE);

            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,5) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            usage(argv[0]);
          }
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-lag") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_lag,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            usage(argv[0]);
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            usage(argv[0]);
            i--;
          }
        }
        else if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_iso,1) )
            exit(EXIT_FAILURE);
          /* if ( ++i < argc && isdigit(argv[i][0]) ) { */
          /*   if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_ls,atof(argv[i])) ) */
          /*     exit(EXIT_FAILURE); */
          /* } */
          /* else if ( i == argc ) { */
          /*   fprintf(stderr,"Missing argument option %c%c\n",argv[i-1][1],argv[i-1][2]); */
          /*   usage(argv[0]); */
          /* } */
          /* else i--; */
        }
        break;
      case 'm':  /* memory */
        if (!strcmp(argv[i],"-m") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_mem,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            usage(argv[0]);
          }
        } else if(!strcmp(argv[i],"-msh") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_msh,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            usage(argv[0]);
          }
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_angle,0) )
            exit(EXIT_FAILURE);
        } else if ( !strcmp(argv[i],"-nsd") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_numsubdomain,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            usage(argv[0]);
          }
        } else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_noswap,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_noinsert,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nomove,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nosurf") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nosurf,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMG2D_Set_outputMeshName(mesh,argv[i]) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            usage(argv[0]);
          }
        }
        break;
      case 'p':
        if ( !strcmp(argv[i],"-per") ) {
          fprintf(stdout,"WARNING OBSOLETE OPTION\n");
          mesh->info.renum = -10;
        }
        break;
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG2D_Set_inputSolName(mesh,met,argv[i]) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            usage(argv[0]);
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        usage(argv[0]);
      }

    }

    else {
      if ( mesh->namein == NULL ) {
        if ( !MMG2D_Set_inputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
        if ( mesh->info.imprim == -99 )  {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,5) )
            exit(EXIT_FAILURE);
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG2D_Set_outputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else {
        fprintf(stdout,"  Argument %s ignored\n",argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,i) )
      exit(EXIT_FAILURE);
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",namein);
    if ( !MMG2D_Set_inputMeshName(mesh,namein) )
      exit(EXIT_FAILURE);
  }
  if ( mesh->nameout == NULL ) {
    if ( !MMG2D_Set_outputMeshName(mesh,"") )
      exit(EXIT_FAILURE);
  }
  if ( met->namein == NULL ) {
    if ( !MMG2D_Set_inputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  if ( met->nameout == NULL ) {
    if ( !MMG2D_Set_outputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  /* if ( mesh->namedep == NULL ) { */
  /*   mesh->namedep = (char *)calloc(128,sizeof(char)); */
  /*   assert(mesh->namedep); */
  /*   strcpy(mesh->namedep,mesh->namein); */
  /*   ptr = strstr(mesh->namedep,".mesh"); */
  /*   if ( ptr ) *ptr = '\0'; */
  /* } */

  return(1);
}


static void endcod() {
  // double   ttot,ttim[TIMEMAX];
  // int      k,call[TIMEMAX];

  //chrono(OFF,&ctim[0]);
//#warning message endcod : comment for merge
  /* for (k=0; k<TIMEMAX; k++) { */
  /*   call[k] = ctim[k].call; */
  /*   ttim[k] = ctim[k].call ? gttime(ctim[k]) : 0.0; */
  /* } */
  /* ttot    = ttim[1]+ttim[2]+ttim[3]+ttim[4]; */
  /* ttim[0] = M_MAX(ttim[0],ttot); */

  /* fprintf(stdout,"\n  -- CPU REQUIREMENTS\n"); */
  /* fprintf(stdout,"  in/out %8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[1]/ttim[0],call[1],ttim[1]/(float)call[1]); */
  /* fprintf(stdout,"  analysis %8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[2]/ttim[0],call[2],ttim[2]/(float)call[2]); */
  /* fprintf(stdout,"  remeshing  %8.2f %%    %3d. calls,   %7.2f sec/call\n", */
  /*         100.*ttim[3]/ttim[0],call[3],ttim[3]/(float)call[3]); */
  /* fprintf(stdout,"  total     %8.2f %%    %3d. calls,  %7.2f sec/call\n", */
  /*         100.*ttot/ttim[0],call[0],ttot/(float)call[0]); */

  /* fprintf(stdout,"\n   ELAPSED TIME  %.2f SEC.  (%.2f)\n",ttim[0],ttot); */
}

int main(int argc,char *argv[]) {
  MMG5_pMesh    mesh;
  MMG5_pSol     sol;
  double        qdegrad[2];
  int           ier,ierSave,msh;
  char          stim[32];

  /* interrupts */
  atexit(endcod);

  _MMG2D_Set_commonFunc();
  tminit(MMG5_ctim,TIMEMAX);
  chrono(ON,&MMG5_ctim[0]);

  /* assign default values */
  mesh = NULL;
  sol  = NULL;

  MMG2D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                  MMG5_ARG_end);

  /* reset default values for file names */
  MMG2D_Free_names(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                   MMG5_ARG_end);

  qdegrad[0] = 10./ALPHA;
  qdegrad[1] = 1.3;

//  sol.type = 1;

  if ( !parsar(argc,argv,mesh,sol,qdegrad) )  return(MMG5_STRONGFAILURE);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&MMG5_ctim[1]);

  /* read mesh file */
  msh = 0;
  ier = MMG2D_loadMesh(mesh,mesh->namein);
  if ( !ier ) {
    ier = MMG2D_loadMshMesh(mesh,sol,mesh->namein);
    msh = 1;
  }

  if ( ier < 1)

    _MMG2D_RETURN_AND_FREE(mesh,sol,MMG5_STRONGFAILURE);
  /* Set default metric size */
  if ( !msh && !MMG2D_Set_solSize(mesh,sol,MMG5_Vertex,0,MMG5_Scalar) )
    _MMG2D_RETURN_AND_FREE(mesh,sol,MMG5_STRONGFAILURE);

  /* read displacement if any */
  if ( mesh->info.lag >= 0 ) {

    if ( !msh ) {
      ier = MMG2D_loadSol(mesh,sol,sol->namein);
      if ( ier < 1 ) {
        fprintf(stdout,"  ## ERROR: UNABLE TO LOAD DISPLACEMENT.\n");
        _MMG2D_RETURN_AND_FREE(mesh,sol,MMG5_STRONGFAILURE);
      }
    }
    if ( sol->size != 2 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE.\n");
      _MMG2D_RETURN_AND_FREE(mesh,sol,MMG5_STRONGFAILURE);
    }
  }
  /* read metric if any */
  else if ( !msh ) {
    ier = MMG2D_loadSol(mesh,sol,sol->namein);
    if ( ier == -1 ) {
      if ( (sol->size != MMG5_Scalar) && (sol->size != MMG5_Tensor) ) {
        fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        _MMG2D_RETURN_AND_FREE(mesh,sol,MMG5_STRONGFAILURE);
      }
    }
    if ( mesh->info.iso && !ier ) {
      fprintf(stdout,"  ## ERROR: NO ISOVALUE DATA.\n");
      _MMG2D_RETURN_AND_FREE(mesh,sol,MMG5_STRONGFAILURE);
    }
  }

  chrono(OFF,&MMG5_ctim[1]);
  printim(MMG5_ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  if ( mesh->info.lag > -1 ) {
    ier = MMG2D_mmg2dmov(mesh,sol);
  }
  else if ( mesh->info.iso ) {
    ier = MMG2D_mmg2dls(mesh,sol);
  }
  else if ( !mesh->nt ) {
    ier = MMG2D_mmg2dmesh(mesh,sol);
  }
  else {
    ier = MMG2D_mmg2dlib(mesh,sol);
  }

  if ( !strcmp(&mesh->nameout[strlen(mesh->nameout)-5],".mesh") ||
       !strcmp(&mesh->nameout[strlen(mesh->nameout)-6],".meshb") )
    msh = 0;
  else if (!strcmp(&mesh->nameout[strlen(mesh->nameout)-4],".msh") ||
           !strcmp(&mesh->nameout[strlen(mesh->nameout)-5],".mshb") )
    msh = 1;

  if ( !msh )
    ierSave = MMG2D_saveMesh(mesh,mesh->nameout);
  else
    ierSave = MMG2D_saveMshMesh(mesh,sol,mesh->nameout);

  if ( !ierSave )
    _MMG2D_RETURN_AND_FREE(mesh,sol,MMG5_STRONGFAILURE);

  if( !msh &&  sol->np )
    MMG2D_saveSol(mesh,sol,mesh->nameout);

  chrono(OFF,&MMG5_ctim[1]);
  if ( mesh->info.imprim ) fprintf(stdout,"  -- WRITING COMPLETED\n");

  /* free mem */
  chrono(OFF,&MMG5_ctim[0]);
  printim(MMG5_ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG2D: ELAPSED TIME  %s\n",stim);

  _MMG2D_RETURN_AND_FREE(mesh,sol,ier);
}
