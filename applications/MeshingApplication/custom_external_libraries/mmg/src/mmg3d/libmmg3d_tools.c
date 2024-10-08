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

#include "libmmg3d.h"
#include "mmgcommon_private.h"
#include "inlined_functions_3d_private.h"
#include "mmgversion.h"
#include "mmg3dexterns_private.h"
#include "mmgexterns_private.h"


void MMG3D_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( mesh->info.ani || (met && met->size == 6) ) {
    /* Force data consistency: if aniso metric is provided, met->size==6 and
     * info.ani==0; with -A option, met->size==1 and info.ani==1 */
    met->size = 6;
    mesh->info.ani = 1;

    if ( (!met->m) && (!mesh->info.optim) && mesh->info.hsiz<=0. ) {
      MMG5_caltet          = MMG5_caltet_iso;
      MMG5_caltri          = MMG5_caltri_iso;
      MMG3D_doSol          = MMG3D_doSol_iso;
      MMG5_lenedg          = MMG5_lenedg_iso;
      MMG3D_lenedgCoor     = MMG5_lenedgCoor_iso;
      MMG5_lenSurfEdg      = MMG5_lenSurfEdg_iso;
    }
    else {
      MMG5_caltet          = MMG5_caltet_ani;
      MMG5_caltri          = MMG5_caltri_ani;
      MMG3D_doSol          = MMG3D_doSol_ani;
      MMG5_lenedg          = MMG5_lenedg_ani;
      MMG3D_lenedgCoor     = MMG5_lenedgCoor_ani;
      MMG5_lenSurfEdg      = MMG5_lenSurfEdg_ani;
    }
    MMG5_intmet          = MMG5_intmet_ani;
    MMG5_lenedgspl       = MMG5_lenedg_ani;
    MMG5_movintpt        = MMG5_movintpt_ani;
    MMG5_movbdyregpt     = MMG5_movbdyregpt_ani;
    MMG5_movbdyrefpt     = MMG5_movbdyrefpt_ani;
    MMG5_movbdynompt     = MMG5_movbdynompt_ani;
    MMG5_movbdyridpt     = MMG5_movbdyridpt_ani;
    MMG5_interp4bar      = MMG5_interp4bar_ani;
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_ani;
    MMG3D_defsiz         = MMG3D_defsiz_ani;
    MMG3D_gradsiz        = MMG3D_gradsiz_ani;
    MMG3D_gradsizreq     = MMG3D_gradsizreq_ani;
#ifndef MMG_PATTERN
    MMG5_cavity          = MMG5_cavity_ani;
    MMG3D_PROctreein     = MMG3D_PROctreein_ani;
#endif
  }
  else {
    if ( mesh->info.optimLES ) {
      MMG5_caltet          = MMG3D_caltetLES_iso;
      MMG5_movintpt        = MMG5_movintpt_iso;
    }
    else {
      MMG5_caltet          = MMG5_caltet_iso;
      MMG5_movintpt        = MMG5_movintpt_iso;
    }
    MMG5_caltri          = MMG5_caltri_iso;
    MMG3D_doSol          = MMG3D_doSol_iso;
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

#ifndef MMG_PATTERN
    MMG5_cavity          = MMG5_cavity_iso;
    MMG3D_PROctreein     = MMG3D_PROctreein_iso;
#endif
  }
}

int MMG3D_Get_adjaTet(MMG5_pMesh mesh, MMG5_int kel, MMG5_int listet[4]) {
  MMG5_int idx;

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

  /* Common generic options, file options and mode options */
  MMG5_mmgUsage(prog);

  /* Lagrangian option (only for mmg2d/3d) */
  MMG5_lagUsage();

  /* Common parameters (first section) */
  MMG5_paramUsage1( );

  /* Parameters shared by mmg2d and 3d only*/
  MMG5_2d3dUsage();

#ifndef MMG_PATTERN
  fprintf(stdout,"-octree val  specify the max number of points per octree cell \n");
#endif
#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]      turn on or off the renumbering using SCOTCH [1/0] \n");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"-nofem       do not force Mmg to create a finite element mesh \n");
  fprintf(stdout,"-nosurf      no surface modifications\n");

  fprintf(stdout,"\n");

  /* Common parameters (second section) */
  MMG5_paramUsage2();

  fprintf(stdout,"-optimLES    enable skewness improvement (for LES computations)\n");

  /* Common options for advanced users */
  MMG5_advancedUsage();

  fprintf(stdout,"\n\n");

  return 1;
}

int MMG3D_defaultValues(MMG5_pMesh mesh) {

  MMG5_mmgDefaultValues(mesh);

#ifndef MMG_PATTERN
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

// In ls mode : metric must be provided using -met option (-sol or default is the ls).
// In adp mode : -sol or -met or default allow to store the metric.
int MMG3D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  MMG5_pSol tmp = NULL;
  int     i;
  char    namein[MMG5_FILESTR_LGTH];

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
        else if ( !strcmp(argv[i],"-isoref") && ++i <= argc ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_isoref,
                                     atoi(argv[i])) )
            return 0;
        }
        else {
          MMG3D_usage(argv[0]);
          return 0;
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
        else if ( !strcmp(argv[i],"-lssurf") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_isosurf,1) )
            return 0;
          if ( ++i < argc && (isdigit(argv[i][0]) ||
                              (argv[i][0]=='-' && isdigit(argv[i][1])) ) ) {
            if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_ls,atof(argv[i])) )
              return 0;
          }
          else i--;
        }
        break;
      case 'm':
        if ( !strcmp(argv[i],"-met") ) {
          if ( !met ) {
            fprintf(stderr,"No metric structure allocated for %c%c%c option\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            return 0;
          }
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
        else if ( !strcmp(argv[i],"-m") ) {
          /* memory */
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_mem,atoi(argv[i])) )
            return 0;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMG3D_usage(argv[0]);
          return 0;
        }
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nofem") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nofem,1) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nreg") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nreg,1) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_angle,0) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nsd") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_numsubdomain,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
            MMG3D_usage(argv[0]);
            return 0;
          }
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
        else if( !strcmp(argv[i],"-nosizreq") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nosizreq,1) ) {
            return 0;
          }
        }
        break;
      case 'o':
        if ( (!strcmp(argv[i],"-out")) || (!strcmp(argv[i],"-o")) ) {
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
#ifndef MMG_PATTERN
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
      case 'r':
        if ( !strcmp(argv[i],"-rmc") ) {
          if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_rmc,0) )
            return 0;
          if ( ++i < argc && (isdigit(argv[i][0]) ) ) {
            if ( !MMG3D_Set_dparameter(mesh,met,MMG3D_DPARAM_rmc,atof(argv[i])) )
              return 0;
          }
          else i--;
        }
#ifdef USE_SCOTCH
        else if ( !strcmp(argv[i],"-rn") ) {
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
#endif
        else {
          fprintf(stderr,"Unrecognized option %s\n",argv[i]);
          MMG3D_usage(argv[0]);
          return 0;
        }
        break;
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          /* For retrocompatibility, store the metric if no sol structure available */
          tmp = sol ? sol : met;
          assert(tmp);
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG3D_Set_inputSolName(mesh,tmp,argv[i]) )
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
      case 'x':
        if ( !strcmp(argv[i],"-xreg") ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_xreg,1) )
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
    MMG_FSCANF(stdin,"%d",&i);
    if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,i) )
      return 0;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%127s",namein);
    if ( !MMG3D_Set_inputMeshName(mesh,namein) )
      return 0;
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMG3D_Set_outputMeshName(mesh,"") )
      return 0;
  }

  /* adp mode: if the metric name has been stored in sol, move it in met */
  if ( met->namein==NULL && sol && sol->namein && !(mesh->info.iso || mesh->info.isosurf || mesh->info.lag>=0) ) {
    if ( !MMG3D_Set_inputSolName(mesh,met,sol->namein) )
      return 0;
    MMG5_DEL_MEM(mesh,sol->namein);
  }

  /* default : store solution (resp. displacement) name in iso
   * (resp. lagrangian) mode, metric name otherwise */
  tmp = ( mesh->info.iso || mesh->info.isosurf || mesh->info.lag >=0 ) ? sol : met;
  assert ( tmp );
  if ( tmp->namein == NULL ) {
    if ( !MMG3D_Set_inputSolName(mesh,tmp,"") ) { return 0; }
  }
  if ( met->nameout == NULL ) {
    if ( !MMG3D_Set_outputSolName(mesh,met,"") )
      return 0;
  }

  return 1;
}

int MMG3D_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  float      fp1,fp2,hausd;
  int        i,j,ret,npar,nbr,split;
  MMG5_int   ref,rin,rex,br;
  char       *ptr,buf[256],data[256];
  FILE       *in;
  fpos_t     position;

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

    /* Read user-defined references for the LS mode */
    if ( !strcmp(data,"lsreferences") ) {
      ret = fscanf(in,"%d",&npar);
      if ( !ret ) {
        fprintf(stderr,"  %%%% Wrong format for lsreferences: %d\n",npar);
        return 0;
      }

      if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_numberOfMat,npar) ) {
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
        if ( !MMG3D_Set_multiMat(mesh,met,ref,split,rin,rex) ) {
          return 0;
        }
      }
    }
    /* Read user-defined local parameters and store them in the structure info->par */
    else if ( !strcmp(data,"parameters") ) {
      ret = fscanf(in,"%d",&npar);

      if ( !ret ) {
        fprintf(stderr,"  %%%% Wrong format for parameters: %d\n",npar);
        return 0;
      }

      if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_numberOfLocalParam,npar) ) {
        return 0;
      }

      for (i=0; i<mesh->info.npar; i++) {
        ret = fscanf(in,"%" MMG5_PRId " %255s ",&ref,buf);
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
    else if ( !strcmp(data,"lsbasereferences") ) {
      MMG_FSCANF(in,"%d",&nbr);
      if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_numberOfLSBaseReferences,nbr) )
        return 0;
      
      for (i=0; i<mesh->info.nbr; i++) {
        MMG_FSCANF(in,"%" MMG5_PRId "",&br);
        if ( !MMG3D_Set_lsBaseReference(mesh,met,br) ) {
          return 0;
        }
      }
    }
  }
  fclose(in);
  return 1;
}

int MMG3D_freeLocalPar(MMG5_pMesh mesh) {

  free(mesh->info.par);
  mesh->info.npar = 0;

  return 1;
}

int MMG3D_Get_numberOfNonBdyTriangles(MMG5_pMesh mesh, MMG5_int* nb_tria) {
  MMG5_pTetra pt,pt1;
  MMG5_pPrism pp;
  MMG5_pTria  ptt;
  MMG5_Hash   hash;
  int         i;
  MMG5_int    ref,*adja,j,k,iel;

  *nb_tria = 0;
  memset ( &hash, 0x0, sizeof(MMG5_Hash));

  if ( !mesh->tetra ) {
    /* No triangle at all */
    return 1;
  }

  /** First step: Mesh analysis to detect the tetra/prisms boundary faces and to
   * store the info in the xtetra/xprisms structures */
  if ( !mesh->adja ) {
    /* create tetra adjacency */
    if ( !MMG3D_hashTetra( mesh,0 ) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to create "
              "adjacency table.\n",__func__);
      return 0;
    }
  }

  if ( !mesh->adjapr ) {
    /* create prism adjacency */
    if ( !MMG3D_hashPrism(mesh) ) {
      fprintf(stderr,"\n  ## Error: %s: Prism hashing problem.\n",__func__);
      return 0;
    }
  }

  /* If mesh->xtetra is filled, we assume that the surface analysis is
   * complete */
  if ( !mesh->xtetra ) {
    /* compatibility triangle orientation w/r tetras */
    if ( !MMG5_bdryPerm(mesh) ) {
      fprintf(stderr,"\n  ## Error: %s: Boundary orientation problem.\n",__func__);
      return 0;
    }
    /* identify surface mesh */
    if ( !MMG5_chkBdryTria(mesh) ) {
      fprintf(stderr,"\n  ## Error: %s: Boundary problem.\n",__func__);
      return 0;
    }
    MMG5_freeXTets(mesh);
    MMG5_freeXPrisms(mesh);

    /* create surface adjacency */
    if ( !MMG3D_hashTria(mesh,&hash) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      fprintf(stderr,"\n  ## Error: %s: Hashing problem.\n",__func__);
      return 0;
    }

    /* identify connexity and flip orientation of faces if needed */
    if ( !MMG5_setadj(mesh) ) {
      fprintf(stderr,"\n  ## Error: %s: Topology problem.\n",__func__);
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }

    /* set bdry entities to tetra and fill the orientation field */
    if ( !MMG5_bdrySet(mesh) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      fprintf(stderr,"\n  ## Error: %s: Boundary problem.\n",__func__);
      return 0;
    }
    MMG5_DEL_MEM(mesh,hash.item);
  }

  /** Second step: Count the number of non boundary faces */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[4*(k-1)+1];

    for ( i=0; i<4; i++ ) {
      iel = adja[i] / 4;
      assert ( iel != k );

      pt1 = &mesh->tetra[iel];

      if ( (!iel) || (pt->ref != pt1->ref) ||
           (mesh->info.opnbdy && pt->xt &&
            (mesh->xtetra[pt->xt].ftag[i] & MG_BDY) ) ) {
        /* Do not treat boundary faces */
        continue;
      }
      if ( k < iel ) {
        /* Treat face from the tetra with lowest index */
        ++(*nb_tria);
      }
    }
  }
  for ( k=1; k<=mesh->nprism; k++ ) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) ) continue;

    adja = &mesh->adjapr[5*(k-1)+1];

    for ( i=0; i<2; i++ ) {
      iel = adja[i] / 5;

      if ( iel<0 ) {
        ref = mesh->tetra[MMG5_abs(iel)].ref;
      } else {
        ref = mesh->prism[iel].ref;
      }

      if ( (!iel) || (pp->ref != ref) ||
           (mesh->info.opnbdy && pp->xpr &&
            (mesh->xprism[pp->xpr].ftag[i] & MG_BDY) ) ) {
        /* Do not treat boundary faces */
        continue;
      }
      if ( k < iel ) {
        /* Treat face from the element with lowest index */
        ++(*nb_tria);
      }
    }
  }

  if ( !(*nb_tria) ) {
    return 1;
  }

  /** Third step: Append the non boundary edges to the boundary edges array */
  if ( mesh->nt ) {
    MMG5_ADD_MEM(mesh,(*nb_tria)*sizeof(MMG5_Tria),"non boundary triangles",
                 printf("  Exit program.\n");
                 MMG5_DEL_MEM(mesh,hash.item);
                 return 0);
    MMG5_SAFE_RECALLOC(mesh->tria,(mesh->nt+1),(mesh->nt+(*nb_tria)+1),
                       MMG5_Tria,"non bdy tria arrray",return 0);
  }
  else {
    MMG5_ADD_MEM(mesh,((*nb_tria)+1)*sizeof(MMG5_Tria),"non boundary triangles",
                 printf("  Exit program.\n");
                 MMG5_DEL_MEM(mesh,hash.item);
                 return 0);
    MMG5_SAFE_RECALLOC(mesh->tria,0,((*nb_tria)+1),
                       MMG5_Tria,"non bdy tria arrray",return 0);
  }

  j = mesh->nt+1;
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[4*(k-1)+1];

    for ( i=0; i<4; i++ ) {
      iel = adja[i] / 4;
      assert ( iel != k );

      pt1 = &mesh->tetra[iel];

      if ( (!iel) || (pt->ref != pt1->ref) ||
           (mesh->info.opnbdy && pt->xt &&
            (mesh->xtetra[pt->xt].ftag[i] & MG_BDY)) ) {
        /* Do not treat boundary faces */
        continue;
      }
      if ( k < iel ) {
        /* Treat edge from the triangle with lowest index */
        ptt = &mesh->tria[j++];
        assert ( ptt );
        ptt->v[0]   = pt->v[MMG5_idir[i][0]];
        ptt->v[1]   = pt->v[MMG5_idir[i][1]];
        ptt->v[2]   = pt->v[MMG5_idir[i][2]];
        ptt->ref    = mesh->xtetra[pt->xt].ref[i];
      }
    }
  }

  for ( k=1; k<=mesh->nprism; k++ ) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) ) continue;

    adja = &mesh->adjapr[5*(k-1)+1];

    for ( i=0; i<2; i++ ) {
      iel = adja[i] / 5;

      if ( iel<0 ) {
        ref = mesh->tetra[MMG5_abs(iel)].ref;
      } else {
        ref = mesh->prism[iel].ref;
      }
      if ( (!iel) || (pp->ref != ref) ||
           (mesh->info.opnbdy && pp->xpr &&
            (mesh->xprism[pp->xpr].ftag[i] & MG_BDY)) ) {
        /* Do not treat boundary faces */
        continue;
      }
      if ( k < iel ) {
        /* Treat edge from the triangle with lowest index */
        ptt = &mesh->tria[j++];
        assert ( ptt );
        ptt->v[0]   = pp->v[MMG5_idir_pr[i][0]];
        ptt->v[1]   = pp->v[MMG5_idir_pr[i][1]];
        ptt->v[2]   = pp->v[MMG5_idir_pr[i][2]];
        ptt->ref    = mesh->xprism[pp->xpr].ref[i];
      }
    }
  }

  return 1;
}

int MMG3D_Get_nonBdyTriangle(MMG5_pMesh mesh,MMG5_int* v0,MMG5_int* v1,MMG5_int* v2,
                             MMG5_int* ref,MMG5_int idx) {
  MMG5_pTria ptt;
  size_t     nt_tot=0;
  char       *ptr_c = (char*)mesh->tria;

  if ( !mesh->tria ) {
    fprintf(stderr,"\n  ## Error: %s: triangle array is not allocated.\n"
            " Please, call the MMG3D_Get_numberOfNonBdyTriangles function"
            " before the %s one.\n",
            __func__,__func__);
    return 0;
  }

  ptr_c = ptr_c-sizeof(size_t);
  nt_tot = *((size_t*)ptr_c);

  if ( mesh->nt==(MMG5_int)nt_tot ) {
    fprintf(stderr,"\n  ## Error: %s: no internal triangle.\n"
            " Please, call the MMG3D_Get_numberOfNonBdyTriangles function"
            " before the %s one and check that the number of internal"
            " triangles is non null.\n",
            __func__,__func__);
    return 0;
  }

  if ( mesh->nt+idx > (MMG5_int)nt_tot ) {
    fprintf(stderr,"\n  ## Error: %s: Can't get the internal triangle of index %" MMG5_PRId "."
            " Index must be between 1 and %zu.\n",
            __func__,idx,nt_tot-(size_t)mesh->nt);
    return 0;
  }

  ptt = &mesh->tria[mesh->nt+idx];

  *v0  = ptt->v[0];
  *v1  = ptt->v[1];
  *v2  = ptt->v[2];

  if ( ref != NULL ) {
    *ref = ptt->ref;
  }

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

int MMG3D_mmg3dcheck(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol,double critmin, double lmin,
                    double lmax, MMG5_int *eltab,int8_t metRidTyp) {

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
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stderr,"\n  ## Error: level-set discretisation unavailable"
            " (MMG3D_IPARAM_iso):\n"
            "          You must call the MMG3D_mmg3dmov function to use this option.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  MMG5_warnOrientation(mesh);

  if ( met ) {
  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,met->m);
    met->np = 0;
  }
  else if ( met->size!=1 && met->size!=6 ) {
    fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }
  if ( sol ) {
    if ( sol->np && (sol->np != mesh->np) ) {
      fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
      MMG5_DEL_MEM(mesh,sol->m);
      sol->np = 0;
    }
    else if ( sol->size!=1 && sol->size!=6 ) {
      fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMG3D_setfunc(mesh,met);

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  %s\n   MODULE MMG3D: IMB-LJLL : %s (%s)\n  %s\n",
            MG_STR,MMG_VERSION_RELEASE,MMG_RELEASE_DATE,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* scaling mesh */
  if ( !MMG5_scaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);


  MMG3D_searchqua(mesh,met,critmin,eltab,metRidTyp);
  ier = MMG3D_searchlen(mesh,met,lmin,lmax,eltab,metRidTyp);
  if ( !ier )
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);

  _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);
}

void MMG3D_searchqua(MMG5_pMesh mesh,MMG5_pSol met,double critmin, MMG5_int *eltab,
                    int8_t metRidTyp) {
  MMG5_pTetra   pt;
  double        rap;
  MMG5_int      k;

  assert ( met );

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

int MMG3D_Get_tetFromTria(MMG5_pMesh mesh, MMG5_int ktri, MMG5_int *ktet, int *iface) {
  MMG5_int val;

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


int MMG3D_Get_tetsFromTria(MMG5_pMesh mesh, MMG5_int ktri, MMG5_int ktet[2], int iface[2])
{
  int      ier;
  MMG5_int itet;
#ifndef NDEBUG
  MMG5_int ia0,ib0,ic0,ia1,ib1,ic1;
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
                    double lmax, MMG5_int *eltab,int8_t metRidTyp) {
  MMG5_pTetra pt;
  MMG5_Hash   hash;
  double      len;
  MMG5_int    k,np,nq;
  int8_t      ia,i0,i1,ier;

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

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 * \param ani 1 for aniso metric, 0 for iso one
 *
 * \return 0 if fail, 1 if succeed.
 *
 * Truncate the metric computed by the DoSol function by hmax and hmin values
 * (if setted by the user). Set hmin and hmax if they are not setted.
 *
 */
static inline
int MMG3D_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met,int ani) {
  MMG5_pTetra pt;
  MMG5_int    k;
  int         i,ier;

  assert ( mesh->info.optim );

  /* Detect the point not used by the mesh */
  ++mesh->base;

#ifndef NDEBUG
  for (k=1; k<=mesh->np; k++) {
    assert ( mesh->point[k].flag < mesh->base );
  }
#endif

  /* For mmg3d, detect points used by tetra */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<4; i++) {
      mesh->point[pt->v[i]].flag = mesh->base;
    }
  }

  /* Compute hmin/hmax on unflagged points and truncate the metric */
  if ( !ani ) {
    ier = MMG5_solTruncature_iso(mesh,met);
  }
  else {
    MMG5_solTruncature_ani = MMG5_3dSolTruncature_ani;
    ier = MMG5_3dSolTruncature_ani(mesh,met);
  }

  return ier;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 1 if succeed, 0 if fail
 *
 * Compute isotropic size map according to the mean of the length of the
 * edges passing through a point.
 *
 */
int MMG3D_doSol_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra  pt;
  MMG5_pPoint  p1,p2;
  double       ux,uy,uz,dd;
  MMG5_int     k,ipa,ipb;
  int          i,ia,ib,type;
  // we guess that we have less than INT32_MAX edges passing through a point
  int         *mark;

  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  /* Memory alloc */
  if ( met->size!=1 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
          }

  type=1;
  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

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

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    if ( !mark[k] ) {
      continue;
    }
    else
      met->m[k] = met->m[k] / (double)mark[k];
  }

  MMG5_SAFE_FREE(mark);

  /* Metric truncation */
  MMG3D_solTruncatureForOptim(mesh,met,0);

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param met pointer toward the metric
 *
 * \return 1 if succeed, 0 if fail
 *
 * Compute anisotropic size map according to the mean of the length of the
 * edges passing through a point.
 *
 */
int MMG3D_doSol_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTetra  pt;
  MMG5_pPoint  p1,p2;
  double       u[3],dd,tensordot[6];
  MMG5_int     k,iadr,ipa,ipb;
  int          i,j,type;
  // we guess that we have less than INT32_MAX edges passing through a point
  int          *mark;

  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  /* Memory alloc */
  if ( met->size!=6 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
  }

  type = 3;
  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  /* Travel the tetra edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<6; i++) {
      ipa = pt->v[MMG5_iare[i][0]];
      ipb = pt->v[MMG5_iare[i][1]];
      p1  = &mesh->point[ipa];
      p2  = &mesh->point[ipb];

      u[0]  = p1->c[0] - p2->c[0];
      u[1]  = p1->c[1] - p2->c[1];
      u[2]  = p1->c[2] - p2->c[2];

      tensordot[0] = u[0]*u[0];
      tensordot[1] = u[0]*u[1];
      tensordot[2] = u[0]*u[2];
      tensordot[3] = u[1]*u[1];
      tensordot[4] = u[1]*u[2];
      tensordot[5] = u[2]*u[2];

      iadr = 6*ipa;
      for ( j=0; j<6; ++j ) {
        met->m[iadr+j]   += tensordot[j];
      }
      mark[ipa]++;

      iadr = 6*ipb;
      for ( j=0; j<6; ++j ) {
        met->m[iadr+j]   += tensordot[j];
      }
      mark[ipb]++;
    }
  }

  for (k=1; k<=mesh->np; k++) {
    if ( !mark[k] ) {
      continue;
    }

    /* Metric = nedges/dim * inv (sum(tensor_dot(edges,edges))).
     * sum(tensor_dot) is stored in sol->m so reuse tensordot to
     * compute M.  */
    iadr = 6*k;
    if ( !MMG5_invmat(met->m+iadr,tensordot) ) {
      /* Non invertible matrix: impose FLT_MIN, it will be truncated by hmax
       * later */
      fprintf(stdout, " ## Warning: %s: %d: non invertible matrix."
             " Impose hmax size at point\n",__func__,__LINE__);
      met->m[iadr+0] = FLT_MIN;
      met->m[iadr+1] = 0;
      met->m[iadr+2] = 0;
      met->m[iadr+3] = FLT_MIN;
      met->m[iadr+4] = 0;
      met->m[iadr+5] = FLT_MIN;
      continue;
    }

    dd = (double)mark[k]/3.;

    for ( j=0; j<6; ++j ) {
      met->m[iadr+j] = dd*tensordot[j];
    }

#ifndef NDEBUG
    /* Check metric */
    double lambda[3],vp[3][3];
    if (!MMG5_eigenv3d(1,met->m+iadr,lambda,vp) ) {
      fprintf(stdout, " ## Warning: %s: %d: non diagonalizable metric.",
              __func__,__LINE__);
    }

    assert ( lambda[0] > 0. && lambda[1] > 0.  && lambda[2] > 0.
             && "Negative eigenvalue" );
    assert ( isfinite(lambda[0]) && isfinite(lambda[1]) && isfinite(lambda[2])
             && "Infinite eigenvalue" );
#endif
  }

  MMG5_SAFE_FREE(mark);

  MMG3D_solTruncatureForOptim(mesh,met,1);

  return 1;
}

int MMG3D_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met) {
  double      hsiz;
  int         type;

  /* Set solution size */
  if ( mesh->info.ani ) {
    met->size = 6;
    type = 3;
  }
  else {
    met->size = 1;
    type = 1;
  }

  /* Memory alloc */
  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) )
    return 0;

  mesh->info.hsiz = hsiz;

  MMG5_Set_constantSize(mesh,met,hsiz);

  return 1;
}

int MMG3D_switch_metricStorage(MMG5_pMesh mesh, MMG5_pSol met) {
  MMG5_int    k;
  double      tmp;

  if ( met->size!=6 ) { return 1; }

  for ( k=1; k<=met->np; ++k ) {

    tmp = met->m[6*k+2];

    met->m[6*k+2] = met->m[6*k+3];
    met->m[6*k+3] = tmp;

  }
  return 1;
}

int MMG3D_Compute_eigenv(double m[6],double lambda[3],double vp[3][3]) {

  return  MMG5_eigenv3d(1,m,lambda,vp);

}

void MMG3D_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol) {

  /* sol */
  if ( !sol ) return;

  if ( sol->m )
    MMG5_DEL_MEM(mesh,sol->m);

  if ( sol->namein ) {
    MMG5_DEL_MEM(mesh,sol->namein);
  }

  if ( sol->nameout ) {
    MMG5_DEL_MEM(mesh,sol->nameout);
  }

  memset ( sol, 0x0, sizeof(MMG5_Sol) );

  /* Reset state to a scalar status */
  sol->dim  = 3;
  sol->ver  = 2;
  sol->size = 1;
  sol->type = 1;

  return;
}

int MMG3D_Clean_isoSurf(MMG5_pMesh mesh) {
  MMG5_int   k,nref;

  nref = 0;
  /** Step 1: a. deletion of triangles that belong to isosurf */
  if ( mesh->tria ) {

    MMG5_int nt = mesh->nt;
    k  = 1;
    do {
      MMG5_pTria ptt = &mesh->tria[k];

      if ( !MG_EOK(ptt) ) {
        continue;
      }

      if ( MMG5_abs(ptt->ref) == mesh->info.isoref ) {
        /* Current tria will be suppressed: search last non isosurf tria to fill
         * empty position */
        MMG5_pTria ptt1 = &mesh->tria[mesh->nt];
        assert( ptt1 );

        while ( ((!MG_EOK(ptt1)) || (MMG5_abs(ptt1->ref) == mesh->info.isoref))
                && k < mesh->nt ) {
          --mesh->nt;
          ptt1 = &mesh->tria[mesh->nt];

        }
        if ( ptt != ptt1 ) {
          /* We don't find any tria to keep after index k */
          memcpy(ptt,ptt1,sizeof(MMG5_Tria));
          --mesh->nt;
        }
      }
      /* Initially negative refs were used to mark isosurface: keep following
       * piece of code for retrocompatibility */
      if ( ptt->ref < 0 ) {
        ptt->ref = -ptt->ref;
        ++nref;
      }
    }
    while ( ++k < mesh->nt );

    /* At the end of the loop, either k==mesh->nt, either k==mesh->nt+1 (because
     * tria at idx mesh->nt was iso or unused and element mesh->nt+1 has been
     * copied into k) */
    assert ( (k==mesh->nt) || (k==mesh->nt+1) );

    /* Check if last element is iso */
    MMG5_pTria ptt = &mesh->tria[mesh->nt];
    if ( (!MG_EOK(ptt)) || (MMG5_abs(ptt->ref) == mesh->info.isoref) ) {
      --mesh->nt;
    }

    if ( mesh->info.imprim > 4 ) {
      fprintf(stdout,"     Deleted iso triangles: %" MMG5_PRId "\n",nt-mesh->nt);
    }

    if( !mesh->nt ) {
      MMG5_DEL_MEM(mesh,mesh->tria);
    }
    else if ( mesh->nt < nt ) {
      MMG5_ADD_MEM(mesh,(mesh->nt-nt)*sizeof(MMG5_Tria),"triangles",
                   fprintf(stderr,"  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->tria,nt+1,(mesh->nt+1),MMG5_Tria,
                         "triangles",return 0);
    }
  }

  /** Step 2: deletion of edges that belong to isosurf */
  return MMG5_Clean_isoEdges(mesh);
}
