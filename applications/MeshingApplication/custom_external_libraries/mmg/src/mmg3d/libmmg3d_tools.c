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
 * \file mmg3d/libmmg3d_tools.c
 * \brief Tools functions for the mmg3d library
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/

#include "inlined_functions_3d.h"

void MMG3D_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size == 1 || ( met->size == 3 && mesh->info.lag >= 0 ) ) {
    if ( mesh->info.optimLES ) {
      MMG5_caltet          = MMG3D_caltetLES_iso;
      MMG5_movintpt        = MMG5_movintpt_iso;
    }
    else {
      MMG5_caltet          = MMG5_caltet_iso;
      MMG5_movintpt        = MMG5_movintpt_iso;
    }
    MMG5_caltri          = MMG5_caltri_iso;
    MMG5_lenedg          = MMG5_lenedg_iso;
    MMG3D_lenedgCoor     = MMG5_lenedgCoor_iso;
    MMG5_lenSurfEdg      = MMG5_lenSurfEdg_iso;
    MMG5_intmet          = MMG5_intmet_iso;
    MMG5_lenedgspl       = MMG5_lenedg_iso;
    MMG5_movbdyregpt     = MMG5_movbdyregpt_iso;
    MMG5_movbdyrefpt     = MMG5_movbdyrefpt_iso;
    MMG5_movbdynompt     = MMG5_movbdynompt_iso;
    MMG5_movbdyridpt     = MMG5_movbdyridpt_iso;
    MMG5_interp4bar      = MMG5_interp4bar_iso;
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_iso;
    MMG3D_defsiz         = MMG3D_defsiz_iso;
    MMG3D_gradsiz        = MMG3D_gradsiz_iso;
    MMG3D_gradsizreq     = MMG3D_gradsizreq_iso;

#ifndef PATTERN
    MMG5_cavity          = MMG5_cavity_iso;
    MMG3D_PROctreein     = MMG3D_PROctreein_iso;
#endif
  }
  else if ( met->size == 6 ) {
    if ( !met->m && !mesh->info.optim && mesh->info.hsiz<=0. ) {
      MMG5_caltet          = MMG5_caltet_iso;
      MMG5_caltri          = MMG5_caltri_iso;
      MMG5_lenedg         = MMG5_lenedg_iso;
      MMG3D_lenedgCoor     = MMG5_lenedgCoor_iso;
      MMG5_lenSurfEdg     = MMG5_lenSurfEdg_iso;
    }
    else {
      MMG5_caltet         = MMG5_caltet_ani;
      MMG5_caltri         = MMG5_caltri_ani;
      MMG5_lenedg         = MMG5_lenedg_ani;
      MMG3D_lenedgCoor     = MMG5_lenedgCoor_ani;
      MMG5_lenSurfEdg     = MMG5_lenSurfEdg_ani;
    }
    MMG5_intmet         = MMG5_intmet_ani;
    MMG5_lenedgspl      = MMG5_lenedg_ani;
    MMG5_movintpt       = MMG5_movintpt_ani;
    MMG5_movbdyregpt     = MMG5_movbdyregpt_ani;
    MMG5_movbdyrefpt     = MMG5_movbdyrefpt_ani;
    MMG5_movbdynompt     = MMG5_movbdynompt_ani;
    MMG5_movbdyridpt     = MMG5_movbdyridpt_ani;
    MMG5_interp4bar      = MMG5_interp4bar_ani;
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_ani;
    MMG3D_defsiz         = MMG3D_defsiz_ani;
    MMG3D_gradsiz        = MMG3D_gradsiz_ani;
    MMG3D_gradsizreq     = MMG3D_gradsizreq_ani;
#ifndef PATTERN
    MMG5_cavity         = MMG5_cavity_ani;
    MMG3D_PROctreein      = MMG3D_PROctreein_ani;
#endif
  }
}

int MMG3D_Get_adjaTet(MMG5_pMesh mesh, int kel, int listet[4]) {
  int idx;

  if ( ! mesh->adja ) {
    if (! MMG3D_hashTetra(mesh, 0))
      return 0;
  }

  idx = 4*(kel-1);
  listet[0] = mesh->adja[idx+1]/4;
  listet[1] = mesh->adja[idx+2]/4;
  listet[2] = mesh->adja[idx+3]/4;
  listet[3] = mesh->adja[idx+4]/4;

  return 1;
}

int MMG3D_usage(char *prog) {

  MMG5_mmgUsage(prog);

  fprintf(stdout,"-A           enable anisotropy (without metric file).\n");
  fprintf(stdout,"-opnbdy      preserve input triangles at the interface of"
          " two domains of the same reference.\n");

#ifdef USE_ELAS
  fprintf(stdout,"-lag [n] Lagrangian mesh displacement according to mode [0/1/2]\n");
  fprintf(stdout,"             0: displacement\n");
  fprintf(stdout,"             1: displacement + remeshing (swap and move)\n");
  fprintf(stdout,"             2: displacement + remeshing (split, collapse,"
          " swap and move)\n");
#endif
#ifndef PATTERN
  fprintf(stdout,"-octree val  Specify the max number of points per octree cell \n");
#endif
#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]      Turn on or off the renumbering using SCOTCH [1/0] \n");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"-nofem       do not force Mmg to create a finite element mesh \n");
  fprintf(stdout,"-optim       mesh optimization\n");
  fprintf(stdout,"-optimLES    strong mesh optimization for LES computations\n");
  fprintf(stdout,"-noinsert    no point insertion/deletion \n");
  fprintf(stdout,"-noswap      no edge or face flipping\n");
  fprintf(stdout,"-nomove      no point relocation\n");
  fprintf(stdout,"-nosurf      no surface modifications\n");
  fprintf(stdout,"\n\n");

  return 1;
}

int MMG3D_defaultValues(MMG5_pMesh mesh) {

  MMG5_mmgDefaultValues(mesh);

#ifndef PATTERN
  fprintf(stdout,"Max number of point per octree cell (-octree) : %d\n",
          mesh->info.PROctree);
#endif
#ifdef USE_SCOTCH
  fprintf(stdout,"SCOTCH renumbering                  : enabled\n");
#else
  fprintf(stdout,"SCOTCH renumbering                  : disabled\n");
#endif
  fprintf(stdout,"\n\n");

  return 1;
}

int MMG3D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  int     i;
  char    namein[128];

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      if ( !MMG3D_defaultValues(mesh) ) return 0;
      return 0;
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        MMG3D_usage(argv[0]);
        return 0;

      case 'a':
        if ( !strcmp(argv[i],"-ar") && ++i < argc )
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_angleDetection,
                                    atof(argv[i])) )
            return 0;
        break;
      case 'A': /* anisotropy */
        if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          return 0;
        break;
      case 'd':
        if ( !strcmp(argv[i],"-default") ) {
          mesh->mark=1;
        }
        else {
          /* debug */
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_debug,1) ) {
            return 0;
          }
        }
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hmin,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hmax,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hsiz") && ++i < argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hsiz,
                                     atof(argv[i])) )
            return 0;

        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hausd,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hgradreq") && ++i <= argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hgradreq,
                                    atof(argv[i])) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hgrad,
                                    atof(argv[i])) )
            return 0;
        }
        else {
          MMG3D_usage(argv[0]);
          return 0;
        }
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG3D_Set_inputMeshName(mesh, argv[i]) )
              return 0;

          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            MMG3D_usage(argv[0]);
            return 0;
          }
        }
        break;

      case 'l':
        if ( !strcmp(argv[i],"-lag") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_lag,atoi(argv[i])) )
              return 0;
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG3D_usage(argv[0]);
            return 0;
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG3D_usage(argv[0]);
            i--;
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_iso,1) )
            return 0;
          if ( ++i < argc && (isdigit(argv[i][0]) ||
                              (argv[i][0]=='-' && isdigit(argv[i][1])) ) ) {
            if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_ls,atof(argv[i])) )
              return 0;
          }
          else i--;
        }
        break;
      case 'm':  /* memory */
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_mem,atoi(argv[i])) )
            return 0;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMG3D_usage(argv[0]);
          return 0;
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nofem") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nofem,1) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_angle,0) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_noswap,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_noinsert,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nomove,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nosurf") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nosurf,1) )
            return 0;
        }
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMG3D_Set_outputMeshName(mesh,argv[i]) )
              return 0;
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMG3D_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-opnbdy") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_opnbdy,1) )
            return 0;
        }
#ifndef PATTERN
        else if ( !strcmp(argv[i],"-octree") && ++i < argc ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_octree,
                                     atoi(argv[i])) )
            return 0;
        }
#endif
        else if( !strcmp(argv[i],"-optimLES") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_optimLES,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-optim") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_optim,1) )
            return 0;
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_renum,atoi(argv[i])) )
                return 0;
            }
            else {
              fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
              MMG3D_usage(argv[0]);
              return 0;
            }
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG3D_usage(argv[0]);
            return 0;
          }
        }
        break;
#endif
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG3D_Set_inputSolName(mesh,met,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMG3D_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( isdigit(argv[i][0]) ||
               (argv[i][0]=='-' && isdigit(argv[i][1])) ) {
            if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,atoi(argv[i])) )
              return 0;
          }
          else {
            i--;
            fprintf(stderr,"Missing argument option %s\n",argv[i]);
          }
        }
        else {
          fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
          MMG3D_usage(argv[0]);
          return 0;
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        MMG3D_usage(argv[0]);
        return 0;
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        if ( !MMG3D_Set_inputMeshName(mesh,argv[i]) )
          return 0;
        if ( mesh->info.imprim == -99 ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,5) )
            return 0;
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG3D_Set_outputMeshName(mesh,argv[i]) )
          return 0;
      }
      else {
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        MMG3D_usage(argv[0]);
        return 0;
      }
    }
    i++;
  }

  /* check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,i) )
      return 0;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%127s",namein);
    if ( !MMG3D_Set_inputMeshName(mesh,namein) )
      return 0;
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMG3D_Set_outputMeshName(mesh,"") )
      return 0;
  }

  if ( met->namein == NULL ) {
    if ( !MMG3D_Set_inputSolName(mesh,met,"") )
      return 0;
  }
  if ( met->nameout == NULL ) {
    if ( !MMG3D_Set_outputSolName(mesh,met,"") )
      return 0;
  }

  return 1;
}

int MMG3D_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  float       fp1,fp2,hausd;
  int         ref,i,j,ret,npar;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  strcpy(data,mesh->namein);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".mmg3d");
  in = fopen(data,"rb");
  if ( !in ) {
    strcat(data,".mmg3d5");
    in = fopen(data,"rb");
    if ( !in ) {
      sprintf(data,"%s","DEFAULT.mmg3d");
      in = fopen(data,"rb");
      if ( !in ) {
        sprintf(data,"%s","DEFAULT.mmg3d5");
        in = fopen(data,"rb");
        if ( !in ) {
          return 1;
        }
      }
    }
  }
  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"\n  %%%% %s OPENED\n",data);

  /* read parameters */
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%255s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"parameters") ) {
      fscanf(in,"%d",&npar);
      if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_numberOfLocalParam,npar) )
        return 0;

      for (i=0; i<mesh->info.npar; i++) {
        ret = fscanf(in,"%d %255s ",&ref,buf);
        if ( ret )
          ret = fscanf(in,"%f %f %f",&fp1,&fp2,&hausd);

        if ( !ret ) {
          fprintf(stderr,"  %%%% Wrong format: %s\n",buf);
          return 0;
        }

        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);

        if ( (!strcmp(buf,"triangles") || !strcmp(buf,"triangle")) ) {
          if ( !MMG3D_Set_localParameter(mesh,met,MMG5_Triangle,ref,fp1,fp2,hausd) ) {
            return 0;
          }
        }
        /* else if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") ) { */
        /*   if ( !MMG3D_Set_localParameter(mesh,met,MMG5_Vertex,ref,fp1,fp2,hausd) ) { */
        /*     return 0; */
        /*   } */
        /* } */
        else if ( !strcmp(buf,"tetrahedra") || !strcmp(buf,"tetrahedron") ) {
          if ( !MMG3D_Set_localParameter(mesh,met,MMG5_Tetrahedron,ref,fp1,fp2,hausd) ) {
            return 0;
          }
        }
        else {
          fprintf(stderr,"  %%%% Wrong format: %s\n",buf);
          return 0;
        }
      }
    }
  }
  fclose(in);
  return 1;
}

int MMG3D_stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(&mesh->info,info,sizeof(MMG5_Info));
  MMG3D_memOption(mesh);
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
      return 0;
    } else if(mesh->info.mem < 39)
      return 0;
  }
  return 1;
}

void MMG3D_destockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(info,&mesh->info,sizeof(MMG5_Info));
  return;
}

int MMG3D_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,double critmin, double lmin,
                    double lmax, int *eltab,char metRidTyp) {

  mytime    ctim[TIMEMAX];
  int       ier;
  char      stim[32];

  MMG3D_Set_commonFunc();

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  MMG3D_Free_topoTables(mesh);

  signal(SIGABRT,MMG5_excfun);
  signal(SIGFPE,MMG5_excfun);
  signal(SIGILL,MMG5_excfun);
  signal(SIGSEGV,MMG5_excfun);
  signal(SIGTERM,MMG5_excfun);
  signal(SIGINT,MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag > -1 ) {
    fprintf(stderr,"\n  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n"
            "            You must call the MMG3D_mmg3dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stderr,"\n  ## Error: level-set discretisation unavailable"
            " (MMG3D_IPARAM_iso):\n"
            "          You must call the MMG3D_mmg3dmov function to use this option.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  MMG5_warnOrientation(mesh);



  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,met->m);
    met->np = 0;
  }
  else if ( met->size!=1 && met->size!=6 ) {
    fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG3D_setfunc(mesh,met);

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* scaling mesh */
  if ( !MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);


  MMG3D_searchqua(mesh,met,critmin,eltab,metRidTyp);
  ier = MMG3D_searchlen(mesh,met,lmin,lmax,eltab,metRidTyp);
  if ( !ier )
    _LIBMMG5_RETURN(mesh,met,MMG5_LOWFAILURE);

  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}

void MMG3D_searchqua(MMG5_pMesh mesh,MMG5_pSol met,double critmin, int *eltab,
                    char metRidTyp) {
  MMG5_pTetra   pt;
  double   rap;
  int      k;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if( !MG_EOK(pt) )
      continue;

    if ( (!metRidTyp) && met->m && met->size>1 ) {
      rap = MMG3D_ALPHAD * MMG5_caltet33_ani(mesh,met,pt);
    }
    else {
      rap = MMG3D_ALPHAD * MMG5_caltet(mesh,met,pt);
    }

    if ( rap == 0.0 || rap < critmin ) {
      eltab[k] = 1;
    }
  }
  return;
}

int MMG3D_Get_tetFromTria(MMG5_pMesh mesh, int ktri, int *ktet, int *iface) {
  int val;

  val = mesh->tria[ktri].cc;

  if ( !val ) {
    fprintf(stderr,"  ## Error: %s: the main fonction of the Mmg library must be"
            " called before this function.\n",__func__);
    return 0;
  }

  *ktet = val/4;

  *iface = val%4;

  return 1;
}


int MMG3D_Get_tetsFromTria(MMG5_pMesh mesh, int ktri, int ktet[2], int iface[2])
{
  int ier;
  int itet;
#ifndef NDEBUG
  int ia0,ib0,ic0,ia1,ib1,ic1;
#endif

  ktet[0]  =  ktet[1] = 0;
  iface[0] = iface[1] = 0;

  ier = MMG3D_Get_tetFromTria(mesh,ktri,ktet,iface);

  if ( !ier ) return 0;

  if ( !mesh->adja ) {
    if (!MMG3D_hashTetra(mesh, 0) )
      return 0;
  }

  itet = mesh->adja[4*(*ktet-1) + *iface + 1 ];

  if ( itet ) {
    ktet[1]  = itet/4;
    iface[1] = itet%4;

#ifndef NDEBUG
    ia0 = mesh->tetra[ktet[0]].v[MMG5_idir[iface[0]][0]];
    ib0 = mesh->tetra[ktet[0]].v[MMG5_idir[iface[0]][1]];
    ic0 = mesh->tetra[ktet[0]].v[MMG5_idir[iface[0]][2]];

    ia1 = mesh->tetra[ktet[1]].v[MMG5_idir[iface[1]][0]];
    ib1 = mesh->tetra[ktet[1]].v[MMG5_idir[iface[1]][1]];
    ic1 = mesh->tetra[ktet[1]].v[MMG5_idir[iface[1]][2]];

    assert ( ( (ia0 == ia1) && ((ib0 == ic1) && (ic0 == ib1 )) ) ||
             ( (ia0 == ib1) && ((ib0 == ia1) && (ic0 == ic1 )) ) ||
             ( (ia0 == ic1) && ((ib0 == ib1) && (ic0 == ia1 )) ) );
#endif
  }

  return 1;
}


int MMG3D_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin,
                   double lmax, int *eltab,char metRidTyp) {
  MMG5_pTetra          pt;
 MMG5_Hash           hash;
  double          len;
  int             k,np,nq;
  char            ia,i0,i1,ier;

  /* Hash all edges in the mesh */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return 0;

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = MMG5_iare[ia][0];
      i1 = MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!MMG5_hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stderr,"\n  ## Error: %s: function MMG5_hashEdge return 0\n",
                __func__);
        return 0;
      }
    }
  }

  /* Pop edges from hash table, and analyze their length */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = MMG5_iare[ia][0];
      i1 = MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = MMG5_hashPop(&hash,np,nq);
      if( ier ) {
        if ( (!metRidTyp) && met->m && met->size>1 ) {
          len = MMG5_lenedg(mesh,met,ia,pt);
        }
        else {
          len = MMG5_lenedg33_ani(mesh,met,ia,pt);
        }

        if( (len < lmin) || (len > lmax) ) {
          eltab[k] = 1;
          break;
        }
      }
    }
  }
  MMG5_DEL_MEM(mesh,hash.item);
  return 1;
}

int MMG3D_doSol(MMG5_pMesh mesh,MMG5_pSol met) {
    MMG5_pTetra  pt;
    MMG5_pPoint  p1,p2;
    double       ux,uy,uz,dd;
    int          i,k,iadr,ia,ib,ipa,ipb,type;
    int         *mark;

    MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

    /* Memory alloc */
    if ( met->size==1 ) type=1;
    else if ( met->size==6 ) type = 3;
    else {
      fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
              __func__,met->size);
      return 0;
    }

    if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
      return 0;

    /* edges */
    dd = 0.;
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !MG_EOK(pt) )  continue;

        if ( met->size == 1 ) {
          for (i=0; i<6; i++) {
            ia  = MMG5_iare[i][0];
            ib  = MMG5_iare[i][1];
            ipa = pt->v[ia];
            ipb = pt->v[ib];
            p1  = &mesh->point[ipa];
            p2  = &mesh->point[ipb];

            ux  = p1->c[0] - p2->c[0];
            uy  = p1->c[1] - p2->c[1];
            uz  = p1->c[2] - p2->c[2];
            dd  = sqrt(ux*ux + uy*uy + uz*uz);

            met->m[ipa] += dd;
            mark[ipa]++;
            met->m[ipb] += dd;
            mark[ipb]++;
          }
        }
        else if ( met->size == 6 ) {
          for (i=0; i<6; i++) {
            ia  = MMG5_iare[i][0];
            ib  = MMG5_iare[i][1];
            ipa = pt->v[ia];
            ipb = pt->v[ib];
            p1  = &mesh->point[ipa];
            p2  = &mesh->point[ipb];

            ux  = p1->c[0] - p2->c[0];
            uy  = p1->c[1] - p2->c[1];
            uz  = p1->c[2] - p2->c[2];
            dd  = sqrt(ux*ux + uy*uy + uz*uz);

            iadr = 6*ipa;
            met->m[iadr]   += dd;
            mark[ipa]++;

            iadr = 6*ipb;
            met->m[iadr] += dd;
            mark[ipb]++;
          }
        }
        else {
          MMG5_SAFE_FREE(mark);
          return 0;
        }
    }

    /* if hmax is not specified, compute it from the metric */
    if ( mesh->info.hmax < 0. ) {
      if ( met->size == 1 ) {
        dd = 0.;
        for (k=1; k<=mesh->np; k++) {
          if ( !mark[k] ) continue;
          dd = MG_MAX(dd,met->m[k]);
        }
        assert ( dd );
      }
      else if ( met->size == 6 ) {
        dd = FLT_MAX;
        for (k=1; k<=mesh->np; k++) {
          if ( !mark[k] ) continue;
          iadr = 6*k;
          dd = MG_MIN(dd,met->m[iadr]);
        }
        assert ( dd < FLT_MAX );
        dd = 1./sqrt(dd);
      }
      mesh->info.hmax = 10.*dd;
    }


    /* vertex size */
    if ( met->size == 1 ) {
      for (k=1; k<=mesh->np; k++) {
        if ( !mark[k] ) {
          met->m[k] = mesh->info.hmax;
          continue;
        }
        else
          met->m[k] = met->m[k] / (double)mark[k];
      }
    }
    else if ( met->size == 6 ) {
      for (k=1; k<=mesh->np; k++) {
        iadr = 6*k;
        if ( !mark[k] ) {
          met->m[iadr]   = 1./(mesh->info.hmax*mesh->info.hmax);
          met->m[iadr+3] = met->m[iadr];
          met->m[iadr+5] = met->m[iadr];
          continue;
        }
        else {
          met->m[iadr]   = (double)mark[k]*(double)mark[k]/(met->m[iadr]*met->m[iadr]);
          met->m[iadr+3] = met->m[iadr];
          met->m[iadr+5] = met->m[iadr];
        }
      }
    }

    MMG5_SAFE_FREE(mark);
    return 1;
}

int MMG3D_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met) {
  double      hsiz;
  int         type;

  /* Memory alloc */
  if ( met->size==1 ) type=1;
  else if ( met->size==6 ) type = 3;
  else {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
  }
  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) )
    return 0;

  mesh->info.hsiz = hsiz;

  MMG5_Set_constantSize(mesh,met,hsiz);

  return 1;
}

int MMG3D_switch_metricStorage(MMG5_pMesh mesh, MMG5_pSol met) {
  int    k;
  double tmp;

  if ( met->size!=6 ) { return 1; }

  for ( k=1; k<=met->np; ++k ) {

    tmp = met->m[6*k+2];

    met->m[6*k+2] = met->m[6*k+3];
    met->m[6*k+3] = tmp;

  }
  return 1;
}
