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
      _MMG5_caltet          = _MMG3D_caltetLES_iso;
      _MMG5_movintpt        = _MMG5_movintptLES_iso;
    }
    else {
      _MMG5_caltet          = _MMG5_caltet_iso;
      _MMG5_movintpt        = _MMG5_movintpt_iso;
    }
    _MMG5_caltri          = _MMG5_caltri_iso;
    _MMG5_lenedg          = _MMG5_lenedg_iso;
    MMG3D_lenedgCoor      = _MMG5_lenedgCoor_iso;
    _MMG5_lenSurfEdg      = _MMG5_lenSurfEdg_iso;
    _MMG5_intmet          = _MMG5_intmet_iso;
    _MMG5_lenedgspl       = _MMG5_lenedg_iso;
    _MMG5_movbdyregpt     = _MMG5_movbdyregpt_iso;
    _MMG5_movbdyrefpt     = _MMG5_movbdyrefpt_iso;
    _MMG5_movbdynompt     = _MMG5_movbdynompt_iso;
    _MMG5_movbdyridpt     = _MMG5_movbdyridpt_iso;
    _MMG5_interp4bar      = _MMG5_interp4bar_iso;
    _MMG5_defsiz          = _MMG3D_defsiz_iso;
    _MMG5_gradsiz         = _MMG5_gradsiz_iso;
#ifndef PATTERN
    _MMG5_cavity          = _MMG5_cavity_iso;
    _MMG3D_octreein       = _MMG3D_octreein_iso;
#endif
  }
  else if ( met->size == 6 ) {
    if ( !met->m ) {
      _MMG5_caltet          = _MMG5_caltet_iso;
      _MMG5_caltri          = _MMG5_caltri_iso;
      _MMG5_lenedg         = _MMG5_lenedg_iso;
      MMG3D_lenedgCoor     = _MMG5_lenedgCoor_iso;
      _MMG5_lenSurfEdg     = _MMG5_lenSurfEdg_iso;
    }
    else {
      _MMG5_caltet         = _MMG5_caltet_ani;
      _MMG5_caltri         = _MMG5_caltri_ani;
      _MMG5_lenedg         = _MMG5_lenedg_ani;
      MMG3D_lenedgCoor     = _MMG5_lenedgCoor_ani;
      _MMG5_lenSurfEdg     = _MMG5_lenSurfEdg_ani;
    }
    _MMG5_intmet         = _MMG5_intmet_ani;
    _MMG5_lenedgspl      = _MMG5_lenedg_ani;
    _MMG5_movintpt       = _MMG5_movintpt_ani;
   _MMG5_movbdyregpt     = _MMG5_movbdyregpt_ani;
   _MMG5_movbdyrefpt     = _MMG5_movbdyrefpt_ani;
   _MMG5_movbdynompt     = _MMG5_movbdynompt_ani;
   _MMG5_movbdyridpt     = _MMG5_movbdyridpt_ani;
    _MMG5_interp4bar     = _MMG5_interp4bar_ani;
    _MMG5_defsiz         = _MMG3D_defsiz_ani;
    _MMG5_gradsiz        = _MMG5_gradsiz_ani;
#ifndef PATTERN
    _MMG5_cavity         = _MMG5_cavity_ani;
    _MMG3D_octreein      = _MMG3D_octreein_ani;
#endif
  }
}

int MMG3D_Get_adjaTet(MMG5_pMesh mesh, int kel, int listet[4]) {
  int idx;

  if ( ! mesh->adja ) {
    if (! MMG3D_hashTetra(mesh, 0))
      return(0);
  }

  idx = 4*(kel-1);
  listet[0] = mesh->adja[idx+1]/4;
  listet[1] = mesh->adja[idx+2]/4;
  listet[2] = mesh->adja[idx+3]/4;
  listet[3] = mesh->adja[idx+4]/4;

  return(1);
}

void MMG3D_usage(char *prog) {

  _MMG5_mmgUsage(prog);
  fprintf(stdout,"-A           enable anisotropy (without metric file).\n");

#ifdef USE_ELAS
  fprintf(stdout,"-lag [0/1/2] Lagrangian mesh displacement according to mode 0/1/2\n");
#endif
#ifndef PATTERN
  fprintf(stdout,"-octree val  Specify the max number of points per octree cell \n");
#endif
#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]      Turn on or off the renumbering using SCOTCH [1/0] \n");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"-optim       mesh optimization\n");
  fprintf(stdout,"-optimLES    strong mesh optimization for LES computations\n");
  fprintf(stdout,"-noinsert    no point insertion/deletion \n");
  fprintf(stdout,"-noswap      no edge or face flipping\n");
  fprintf(stdout,"-nomove      no point relocation\n");
  fprintf(stdout,"-nosurf      no surface modifications\n");
  fprintf(stdout,"\n\n");

  exit(EXIT_FAILURE);
}

void MMG3D_defaultValues(MMG5_pMesh mesh) {

  _MMG5_mmgDefaultValues(mesh);

#ifndef PATTERN
  fprintf(stdout,"Max number of point per octree cell (-octree) : %d\n",
          mesh->info.octree);
#endif
#ifdef USE_SCOTCH
  fprintf(stdout,"SCOTCH renumbering                  : enabled\n");
#else
  fprintf(stdout,"SCOTCH renumbering                  : disabled\n");
#endif
  fprintf(stdout,"\n\n");

  exit(EXIT_FAILURE);
}

int MMG3D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  int     i;
  char    namein[128];

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      MMG3D_defaultValues(mesh);
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        MMG3D_usage(argv[0]);
        break;

      case 'a':
        if ( !strcmp(argv[i],"-ar") && ++i < argc )
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_angleDetection,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        break;
      case 'A': /* anisotropy */
        if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          exit(EXIT_FAILURE);
        break;
      case 'd':
        if ( !strcmp(argv[i],"-default") ) {
          mesh->mark=1;
        }
        else {
          /* debug */
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_debug,1) ) {
            exit(EXIT_FAILURE);
          }
        }
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hmin,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hmax,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hausd,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_hgrad,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else
          MMG3D_usage(argv[0]);
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG3D_Set_inputMeshName(mesh, argv[i]) )
              exit(EXIT_FAILURE);

          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            MMG3D_usage(argv[0]);
          }
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-lag") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_lag,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else if ( i == argc ) {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG3D_usage(argv[0]);
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG3D_usage(argv[0]);
            i--;
          }
        }
        else if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_iso,1) )
            exit(EXIT_FAILURE);
          if ( ++i < argc && (isdigit(argv[i][0]) ||
                              (argv[i][0]=='-' && isdigit(argv[i][1])) ) ) {
            if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_ls,atof(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else i--;
        }
        break;
      case 'm':  /* memory */
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_mem,atoi(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMG3D_usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_angle,0) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_noswap,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_noinsert,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nomove,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nosurf") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nosurf,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMG3D_Set_outputMeshName(mesh,argv[i]) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMG3D_usage(argv[0]);
          }
        }
#ifndef PATTERN
        else if ( !strcmp(argv[i],"-octree") && ++i < argc ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_octree,
                                     atoi(argv[i])) )
            exit(EXIT_FAILURE);
        }
#endif
        else if( !strcmp(argv[i],"-optimLES") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_optimLES,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-optim") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_optim,1) )
            exit(EXIT_FAILURE);
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_renum,atoi(argv[i])) )
                exit(EXIT_FAILURE);
            }
            else {
              fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
              MMG3D_usage(argv[0]);
            }
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMG3D_usage(argv[0]);
          }
        }
        break;
#endif
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG3D_Set_inputSolName(mesh,met,argv[i]) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMG3D_usage(argv[0]);
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMG3D_usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        MMG3D_usage(argv[0]);
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        if ( !MMG3D_Set_inputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
        if ( mesh->info.imprim == -99 ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,5) )
            exit(EXIT_FAILURE);
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG3D_Set_outputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else {
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        MMG3D_usage(argv[0]);
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
      exit(EXIT_FAILURE);
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",namein);
    if ( !MMG3D_Set_inputMeshName(mesh,namein) )
      exit(EXIT_FAILURE);
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMG3D_Set_outputMeshName(mesh,"") )
      exit(EXIT_FAILURE);
  }

  if ( met->namein == NULL ) {
    if ( !MMG3D_Set_inputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  if ( met->nameout == NULL ) {
    if ( !MMG3D_Set_outputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }

  return(1);
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
          return(1);
        }
      }
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /* read parameters */
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"parameters") ) {
      fscanf(in,"%d",&npar);
      if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_numberOfLocalParam,npar) )
        exit(EXIT_FAILURE);

      for (i=0; i<mesh->info.npar; i++) {
        ret = fscanf(in,"%d %s ",&ref,buf);
        ret = fscanf(in,"%f %f %f",&fp1,&fp2,&hausd);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);

        if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") ) {
          if ( !MMG3D_Set_localParameter(mesh,met,MMG5_Triangle,ref,fp1,fp2,hausd) ) {
            exit(EXIT_FAILURE);
          }
        }
        /* else if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") ) { */
        /*   if ( !MMG3D_Set_localParameter(mesh,met,MMG5_Vertex,ref,fp1,fp2,hausd) ) { */
        /*     exit(EXIT_FAILURE); */
        /*   } */
        /* } */
        else if ( !strcmp(buf,"tetrahedra") || !strcmp(buf,"tetrahedron") ) {
          if ( !MMG3D_Set_localParameter(mesh,met,MMG5_Tetrahedron,ref,fp1,fp2,hausd) ) {
            exit(EXIT_FAILURE);
          }
        }
        else {
          fprintf(stderr,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
      }
    }
  }
  fclose(in);
  return(1);
}

int MMG3D_stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(&mesh->info,info,sizeof(MMG5_Info));
  _MMG3D_memOption(mesh);
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
      return(0);
    } else if(mesh->info.mem < 39)
      return(0);
  }
  return(1);
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

  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- MMG3d, Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  _MMG3D_Set_commonFunc();

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  _MMG3D_Free_topoTables(mesh);

  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag > -1 ) {
    fprintf(stderr,"  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n"
            "            You must call the MMG3D_mmg3dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stderr,"  ## Error: level-set discretisation unavailable"
            " (MMG3D_IPARAM_iso):\n"
            "          You must call the MMG3D_mmg3dmov function to use this option.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  _MMG5_warnOrientation(mesh);



  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 && met->size!=6 ) {
    fprintf(stderr,"  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG3D_setfunc(mesh,met);

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* scaling mesh */
  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);


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
      rap = _MMG5_ALPHAD * _MMG5_caltet33_ani(mesh,met,pt);
    }
    else {
      rap = _MMG5_ALPHAD * _MMG5_caltet(mesh,met,pt);
    }

    if ( rap == 0.0 || rap < critmin ) {
      eltab[k] = 1;
    }
  }
  return;
}

int MMG3D_Get_tetFromTria(MMG5_pMesh mesh, int ktri, int *ktet, int *iface)
{
  int val;

  val = mesh->tria[ktri].cc;

  if ( !val ) return(0);

  *ktet = val/4;

  *iface = val%4;

  return 1;
}

int MMG3D_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin,
                   double lmax, int *eltab,char metRidTyp) {
  MMG5_pTetra          pt;
 _MMG5_Hash           hash;
  double          len;
  int             k,np,nq;
  char            ia,i0,i1,ier;

  /* Hash all edges in the mesh */
  if ( !_MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return(0);

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = _MMG5_iare[ia][0];
      i1 = _MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!_MMG5_hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stderr,"%s:%d: Error: function _MMG5_hashEdge return 0\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
      }
    }
  }

  /* Pop edges from hash table, and analyze their length */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = _MMG5_iare[ia][0];
      i1 = _MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = _MMG5_hashPop(&hash,np,nq);
      if( ier ) {
        if ( (!metRidTyp) && met->m && met->size>1 ) {
          len = _MMG5_lenedg(mesh,met,ia,pt);
        }
        else {
          len = _MMG5_lenedg33_ani(mesh,met,ia,pt);
        }

        if( (len < lmin) || (len > lmax) ) {
          eltab[k] = 1;
          break;
        }
      }
    }
  }
  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(1);
}

int MMG3D_doSol(MMG5_pMesh mesh,MMG5_pSol met) {
    MMG5_pTetra     pt;
    MMG5_pPoint     p1,p2;
    double     ux,uy,uz,dd;
    int        i,k,ia,ib,ipa,ipb;
    int       *mark;

    _MMG5_SAFE_CALLOC(mark,mesh->np+1,int);

    /* Memory alloc */
    met->np     = mesh->np;
    met->npmax  = mesh->npmax;
    met->size   = 1;
    met->dim    = mesh->dim;

    _MMG5_ADD_MEM(mesh,(met->size*(met->npmax+1))*sizeof(double),"solution",return(0));
    _MMG5_SAFE_CALLOC(met->m,met->size*(met->npmax+1),double);

    /* edges */
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( !pt->v[0] )  continue;

        for (i=0; i<6; i++) {
            ia  = _MMG5_iare[i][0];
            ib  = _MMG5_iare[i][1];
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

    /* vertex size */
    for (k=1; k<=mesh->np; k++) {
        p1 = &mesh->point[k];
        if ( !mark[k] ) {
            met->m[k] = FLT_MAX;
            continue;
        }
        else
            met->m[k] = met->m[k] / (double)mark[k];
    }
    _MMG5_SAFE_FREE(mark);
    return(1);
}

/** Old API °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°*/
int MMG5_Get_adjaTet(MMG5_pMesh mesh, int kel, int *v0, int *v1, int *v2, int *v3) {
  printf("  ##  MMG5_Get_adjaTet: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );
  int listet[4],ier;

  ier = MMG3D_Get_adjaTet(mesh,kel,listet);
  if (ier!=1) return ier;

  *v0=listet[0];
  *v1=listet[1];
  *v2=listet[2];
  *v3=listet[3];
  return 1;
}

void MMG5_usage(char *prog) {
  printf("  ## MMG5_usage: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  MMG3D_usage(prog);
  return;
}

void MMG5_defaultValues(MMG5_pMesh mesh) {
  printf("  ## MMG5_defaultValues: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  MMG3D_defaultValues(mesh);
  return;
}

int MMG5_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  printf("  ## MMG5_parsar: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  return(MMG3D_parsar(argc,argv, mesh, met) );
}

int MMG5_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  printf("  ## MMG5_parsop: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  return(MMG3D_parsop(mesh,met));
}

int MMG5_stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {
  printf("  ## MMG5_stockOptions: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  return(MMG3D_stockOptions(mesh,info) );
}

int MMG5_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,double critmin, double lmin,
                    double lmax, int *eltab,char metRidTyp) {
  printf("  ## MMG5_mmg3dcheck: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  return(MMG3D_mmg3dcheck( mesh, met, critmin, lmin,lmax, eltab, metRidTyp));
}

void MMG5_searchqua(MMG5_pMesh mesh,MMG5_pSol met,double critmin, int *eltab,
                    char metRidTyp) {
  printf("  ## MMG5_searchqua: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  MMG3D_searchqua( mesh, met,critmin,eltab, metRidTyp);
  return;
}

int MMG5_searchlen(MMG5_pMesh mesh, MMG5_pSol met, double lmin,
                   double lmax, int *eltab,char metRidTyp) {
  printf("  ## MMG5_searchlen: "
         "MMG5_ API is deprecated (replaced by the MMG3D_ one) and will"
        " be removed soon\n." );

  return(MMG3D_searchlen(mesh, met,  lmin,lmax,eltab,metRidTyp));
}
