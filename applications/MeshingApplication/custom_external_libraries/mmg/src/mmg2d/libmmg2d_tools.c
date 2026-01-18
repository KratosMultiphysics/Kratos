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
 * \file mmg2d/libmmg2d_tools.c
 * \brief Tools functions for the mmg2d library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 **/

#include "libmmg2d.h"
#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"
#include "mmgexterns_private.h"

void MMG2D_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( mesh->info.ani || (met && met->size==3 ) ) {
    /* Force data consistency: if aniso metric is provided, met->size==3 and
     * info.ani==0; with -A option, met->size==1 and info.ani==1 */
    met->size = 3;
    mesh->info.ani = 1;

    /* Set pointers */
    MMG2D_lencurv    = MMG2D_lencurv_ani;
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_ani;
    MMG2D_defsiz     = MMG2D_defsiz_ani;
    MMG2D_gradsiz    = lissmet_ani;
    MMG2D_gradsizreq = MMG5_gradsizreq_ani;
    MMG2D_caltri     = MMG2D_caltri_ani;
    MMG2D_intmet     = MMG2D_intmet_ani;
    MMG2D_doSol      = MMG2D_doSol_ani;
  }
  else {
    MMG2D_lencurv    = MMG2D_lencurv_iso;
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_iso;
    MMG2D_defsiz     = MMG2D_defsiz_iso;
    MMG2D_gradsiz    = MMG5_gradsiz_iso;
    MMG2D_gradsizreq = MMG5_gradsizreq_iso;
    MMG2D_caltri     = MMG2D_caltri_iso;
    MMG2D_intmet     = MMG2D_intmet_iso;
    MMG2D_doSol      = MMG2D_doSol_iso;
  }
  return;
}

int MMG2D_usage(char *name) {

  /* Common generic options, file options and mode options */
  MMG5_mmgUsage(name);

  /* Lagrangian option (only for mmg2d/3d) */
  MMG5_lagUsage();

  /* Common parameters (first section) */
  MMG5_paramUsage1();

  /* Parameters shared by mmg2d and 3d only*/
  MMG5_2d3dUsage();

  /* Specific parameters */
  fprintf(stdout,"-3dMedit val read and write for gmsh visu: output only if val=1, input and output if val=2, input if val=3\n");
  fprintf(stdout,"\n");

  fprintf(stdout,"-nofem       do not force Mmg to create a finite element mesh \n");
  fprintf(stdout,"-nosurf      no surface modifications\n");

  /* Common parameters (second section) */
  MMG5_paramUsage2();

  /* Common options for advanced users */
  MMG5_advancedUsage();

  fprintf(stdout,"\n\n");

  return 1;
}

// In ls mode : metric must be provided using -met option (-sol or default is the ls).
// In adp mode : -sol or -met or default allow to store the metric.
int MMG2D_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  MMG5_pSol tmp = NULL;
  double  val;
  int     i;
  char    namein[MMG5_FILESTR_LGTH],*endptr;
  int     param;

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      MMG2D_defaultValues(mesh);
      return 0;
    }
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) ) {
      MMG2D_usage(argv[0]);
      return 0;
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case 'a':
        if ( !strcmp(argv[i],"-ar") ) {
          if ( i >= argc -1 ) {
            fprintf(stderr,"\nMissing argument option %s\n",argv[i]);
            return 0;
          }
          else {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_angleDetection,val))
                return 0;
            }
            else {
              /* argument is not a number */
              fprintf(stderr,"\nMissing argument option %s\n",argv[i]);
              return 0;
            }
          }
        }
        break;
      case 'A': /* anisotropy */
        if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          return 0;
        break;
      case 'd':
        if ( !strcmp(argv[i],"-default") ) {
          mesh->mark=1;
        } else {  /* debug */
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_debug,1) )
            return 0;
        }
        break;
      case 'f':
        if ( !strcmp(argv[i],"-f") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG2D_Set_inputParamName(mesh,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing filename for %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'h':
        param = MMG5_UNSET;
        if ( i >= argc -1 ) {
          fprintf(stderr,"\nMissing argument option %s\n",argv[i]);
          return 0;
        }
        else {
          if ( !strcmp(argv[i],"-hmin") ) {
            param = MMG2D_DPARAM_hmin;
          }
          else if ( !strcmp(argv[i],"-hmax") ) {
            param = MMG2D_DPARAM_hmax;
          }
          else if ( !strcmp(argv[i],"-hsiz") ) {
            param = MMG2D_DPARAM_hsiz;
          }
          else if ( !strcmp(argv[i],"-hausd") ) {
            param = MMG2D_DPARAM_hausd;
          }
          else if ( !strcmp(argv[i],"-hgradreq") ) {
            param = MMG2D_DPARAM_hgradreq;
          }
          else if ( !strcmp(argv[i],"-hgrad") ) {
            param = MMG2D_DPARAM_hgrad;
          }
          else {
            /* Arg unknown by Mmg: arg starts with -h but is not known */
            MMG2D_usage(argv[0]);
            return 0;
          }

          assert ( param != MMG5_UNSET );

          val = strtof(argv[i+1],&endptr);
          if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
            ++i;
            if ( !MMG2D_Set_dparameter(mesh,met,param,val) ){
              return 0;
            }
          } else {
            /* argument is not a number */
            fprintf(stderr,"\nMissing argument option %s\n",argv[i]);
            return 0;
          }
        }

        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG2D_Set_inputMeshName(mesh, argv[i]) )
              return 0;

            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,5) )
              return 0;
          }else{
            fprintf(stderr,"\nMissing filname for %s\n",argv[i]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-isoref") && ++i <= argc ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_isoref,
                                     atoi(argv[i])) )
            return 0;
        }
        else {
          MMG2D_usage(argv[0]);
          return 0;
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-lag") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_lag,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing or unexpected argument option %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_iso,1) )
            return 0;

          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_ls,val))
                return 0;
            }
          }
        }
        else if ( !strcmp(argv[i],"-lssurf") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_isosurf,1) )
            return 0;

          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_ls,val))
                return 0;
            }
          }
        }
        break;
      case 'm':  /* memory */
        if ( !strcmp(argv[i],"-met") ) {
          if ( !met ) {
            fprintf(stderr,"\nNo metric structure allocated for %s option\n",
                    argv[i]);
            return 0;
          }
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG2D_Set_inputSolName(mesh,met,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing filname for %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        else if (!strcmp(argv[i],"-m") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_mem,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing argument option %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nofem") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nofem,1) )
            return 0;
        }
        if ( !strcmp(argv[i],"-nreg") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nreg,1) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_angle,0) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nsd") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_numsubdomain,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing argument option %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        } else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_noswap,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_noinsert,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nomove,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nosurf") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nosurf,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nosizreq") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_nosizreq,1) ) {
            return 0;
          }
        }
        break;
      case 'o':
        if ( (!strcmp(argv[i],"-out")) || (!strcmp(argv[i],"-o")) ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMG2D_Set_outputMeshName(mesh,argv[i]) )
              return 0;
          }else{
            fprintf(stderr,"\nMissing filname for %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-opnbdy") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_opnbdy,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-optim") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_optim,1) )
            return 0;
        }
        break;
      case 'r':
        if ( !strcmp(argv[i],"-rmc") ) {
          if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_rmc,0) )
            return 0;
          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_rmc,val))
                return 0;
            }
          }
        }
        break;
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          /* For retrocompatibility, store the metric if no sol structure available */
          tmp = sol ? sol : met;

          assert(tmp);
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMG2D_Set_inputSolName(mesh,tmp,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing filname for %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,atoi(argv[i])) )
              return 0;
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"\nMissing argument option for %s\n",argv[i-1]);
          MMG2D_usage(argv[0]);
          return 0;
        }
        break;
      case 'x':
        if ( !strcmp(argv[i],"-xreg") ) {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_xreg,1) )
            return 0;
          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMG2D_Set_dparameter(mesh,met,MMG2D_DPARAM_xreg,val))
                return 0;
            }
          }
        }
        break;
      case '3':
        if(!strcmp(argv[i],"-3dMedit") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_3dMedit,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing argument option for %s\n",argv[i-1]);
            MMG2D_usage(argv[0]);
            return 0;
          }
        }
        break;
      default:
        fprintf(stderr,"\nUnrecognized option %s\n",argv[i]);
        MMG2D_usage(argv[0]);
        return 0;
      }

    }

    else {
      if ( mesh->namein == NULL ) {
        if ( !MMG2D_Set_inputMeshName(mesh,argv[i]) )
          return 0;
        if ( mesh->info.imprim == -99 )  {
          if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,5) )
            return 0;
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG2D_Set_outputMeshName(mesh,argv[i]) )
          return 0;
      }
      else {
        fprintf(stdout,"\nArgument %s ignored\n",argv[i]);
        MMG2D_usage(argv[0]);
        return 0;
      }
    }
    i++;
  }

  /** check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%d",&i);
    if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_verbose,i) )
      return 0;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"\n  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%127s",namein);
    if ( !MMG2D_Set_inputMeshName(mesh,namein) )
      return 0;
  }
  if ( mesh->nameout == NULL ) {
    if ( !MMG2D_Set_outputMeshName(mesh,"") )
      return 0;
  }

  /* adp mode: if the metric name has been stored in sol, move it in met */
  if ( met->namein==NULL && sol && sol->namein && !(mesh->info.iso || mesh->info.isosurf || mesh->info.lag>=0) ) {
    if ( !MMG2D_Set_inputSolName(mesh,met,sol->namein) )
      return 0;
    MMG5_DEL_MEM(mesh,sol->namein);
  }

  /* default : store solution (resp. displacement) name in iso
   * (resp. lagrangian) mode, metric name otherwise */
  tmp = ( mesh->info.iso || mesh->info.isosurf || mesh->info.lag >=0 ) ? sol : met;
  assert ( tmp );
  if ( tmp->namein == NULL ) {
    if ( !MMG2D_Set_inputSolName(mesh,tmp,"") )
      return 0;
  }
  if ( met->nameout == NULL ) {
    if ( !MMG2D_Set_outputSolName(mesh,met,"") )
      return 0;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
int MMG2D_defaultValues(MMG5_pMesh mesh) {

  MMG5_mmgDefaultValues(mesh);

  fprintf(stdout,"\n\n");

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param met pointer to the metric
 *
 * \return 1 if success, 0 if fail
 *
 * Read parameter file DEFAULT.mmg2d
 *
 */
int MMG2D_parsop(MMG5_pMesh mesh,MMG5_pSol met) {
  int        ret,i,j,npar,nbr,split;
  MMG5_int   ref,rin,rex,br;
  float      fp1,fp2,fp3;
  char       *ptr,data[256];
  FILE       *in;
  fpos_t     position;

  /* Check for parameter file */
  if (mesh->info.fparam) {
    strcpy(data,mesh->info.fparam);
  }
  else {
    strcpy(data,mesh->namein);
  }

  ptr = MMG5_Get_filenameExt(data);

  if ( ptr ) *ptr = '\0';
  strcat(data,".mmg2d");

  in = fopen(data,"rb");

  if ( !in ) {
    if ( !mesh->info.fparam ) {
      sprintf(data,"%s","DEFAULT.mmg2d");
      in = fopen(data,"rb");
      if ( !in )
        return 1;
    }
    else if (mesh->info.fparam ) {
      fprintf(stderr,"  ** In %s: %s file NOT FOUND. \n",__func__,data);
      fprintf(stdout,"  ## ERROR: UNABLE TO LOAD PARAMETER FILE.\n");
      return 0;
    }
  }
  if ( mesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n  %%%% %s OPENED\n",data);
  }

  /* Read parameters */
  while ( !feof(in) ) {
    ret = fscanf(in,"%255s",data);
    if ( !ret || feof(in) ) break;
    for (i=0; (size_t)i<strlen(data); i++) data[i] = tolower(data[i]);

    /* Read user-defined references for the LS mode */
    if ( !strcmp(data,"lsreferences") ) {
      ret = fscanf(in,"%d",&npar);
      if ( !ret ) {
        fprintf(stderr,"  %%%% Wrong format for lsreferences: %d\n",npar);
        return 0;
      }

      if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_numberOfMat,npar) ) {
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
        if ( !MMG2D_Set_multiMat(mesh,met,ref,split,rin,rex) ) {
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
      else if ( npar > MMG5_LPARMAX ) {
        fprintf(stderr,"  %%%% Too many local parameters %d. Abort\n",npar);
        return 0;
      }

      /* Allocate memory and fill the info->par table (adding one, corresponding to the command line data) */
      if ( npar ) {
        if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_numberOfLocalParam,npar) )
          return 0;

        for (i=0; i<mesh->info.npar; i++) {
          ret = fscanf(in,"%" MMG5_PRId " %255s",&ref,data);
          if ( ret ) ret = fscanf(in,"%f %f %f",&fp1,&fp2,&fp3);

          if ( !ret ) {
            fprintf(stderr,"  %%%% Wrong format: %s\n",data);
            return (0);
          }

          for (j=0; (size_t)j<strlen(data); j++) data[j] = tolower(data[j]);
          if ( !strcmp(data,"triangles") || !strcmp(data,"triangle") ) {
            if ( !MMG2D_Set_localParameter(mesh,met,MMG5_Triangle,ref,fp1,fp2,fp3) ) {
              return 0;
            }
          }
          else if ( !strcmp(data,"edges") || !strcmp(data,"edge") ) {
            if ( !MMG2D_Set_localParameter(mesh,met,MMG5_Edg,ref,fp1,fp2,fp3) ) {
              return 0;
            }
          }
          else {
            fprintf(stderr,"  %%%% Wrong format: %s\n",data);
            return 0;
          }
        }
      }
    }
    /* Read user-defined references where connected components should stay attached in ls mode */
    else if ( !strcmp(data,"lsbasereferences") ) {
      MMG_FSCANF(in,"%d",&nbr);
      if ( !MMG2D_Set_iparameter(mesh,met,MMG2D_IPARAM_numberOfLSBaseReferences,nbr) )
        return 0;

      for (i=0; i<mesh->info.nbr; i++) {
        MMG_FSCANF(in,"%" MMG5_PRId "",&br);
        if ( !MMG2D_Set_lsBaseReference(mesh,met,br) ) {
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

/* Free the structure dedicated to the management of multiple local parameters */
int MMG2D_freeLocalPar(MMG5_pMesh mesh) {

  free(mesh->info.par);
  mesh->info.npar = 0;

  return 1;
}

int MMG2D_Get_numberOfNonBdyEdges(MMG5_pMesh mesh, MMG5_int* nb_edges) {
  MMG5_pTria pt,pt1;
  MMG5_pEdge ped;
  MMG5_int   *adja,k,i,j,i1,i2,iel;

  *nb_edges = 0;
  if ( mesh->tria ) {
    /* Create the triangle adjacency if needed */
    if ( !mesh->adja ) {
      if ( !MMG2D_hashTria( mesh ) ) {
        fprintf(stderr,"\n  ## Error: %s: unable to create "
                "adjacency table.\n",__func__);
        return 0;
      }
    }

    /* Count the number of non boundary edges */
    for ( k=1; k<=mesh->nt; k++ ) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      adja = &mesh->adja[3*(k-1)+1];

      for ( i=0; i<3; i++ ) {
        iel = adja[i] / 3;
        assert ( iel != k );

        pt1 = &mesh->tria[iel];

        if ( (!iel) || (pt->ref != pt1->ref) ||
             ((pt->ref==pt1->ref) && MG_SIN(pt->tag[i])) ||
             (mesh->info.opnbdy && pt->tag[i]) ) {
          /* Do not treat boundary edges */
          continue;
        }
        if ( k < iel ) {
          /* Treat edge from the triangle with lowest index */
          ++(*nb_edges);
        }
      }
    }

    /* Append the non boundary edges to the boundary edges array */
    if ( mesh->namax ) {
      MMG5_ADD_MEM(mesh,(*nb_edges)*sizeof(MMG5_Edge),"non boundary edges",
                   printf("  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->edge,(mesh->namax+1),(mesh->namax+(*nb_edges)+1),
                         MMG5_Edge,"non bdy edges arrray",return 0);
    }
    else {
      MMG5_ADD_MEM(mesh,((*nb_edges)+1)*sizeof(MMG5_Edge),"non boundary edges",
                   printf("  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->edge,0,((*nb_edges)+1),
                         MMG5_Edge,"non bdy edges arrray",return 0);
    }

    j = mesh->namax+1;
    for ( k=1; k<=mesh->nt; k++ ) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      adja = &mesh->adja[3*(k-1)+1];

      for ( i=0; i<3; i++ ) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        iel = adja[i] / 3;
        assert ( iel != k );

        pt1 = &mesh->tria[iel];

        if ( (!iel) || (pt->ref != pt1->ref) ||
             ((pt->ref==pt1->ref) && MG_SIN(pt->tag[i])) ||
             (mesh->info.opnbdy && pt->tag[i]) ) {
          /* Do not treat boundary edges */
          continue;
        }
        if ( k < iel ) {
          /* Treat edge from the triangle with lowest index */
          ped = &mesh->edge[j++];
          assert ( ped );
          ped->a   = pt->v[i1];
          ped->b   = pt->v[i2];
          ped->ref = pt->edg[i];
        }
      }
    }
  }
  return 1;
}

int MMG2D_Get_nonBdyEdge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref, MMG5_int idx) {
  MMG5_pEdge ped;
  size_t     na_tot=0;
  char       *ptr_c = (char*)mesh->edge;

  if ( !mesh->edge ) {
    fprintf(stderr,"\n  ## Error: %s: edge array is not allocated.\n"
            " Please, call the MMG2D_Get_numberOfNonBdyEdges function"
            " before the %s one.\n",
            __func__,__func__);
    return 0;
  }

  ptr_c = ptr_c-sizeof(size_t);
  na_tot = *((size_t*)ptr_c);

  if ( mesh->namax==(MMG5_int)na_tot ) {
    fprintf(stderr,"\n  ## Error: %s: no internal edge.\n"
            " Please, call the MMG2D_Get_numberOfNonBdyEdges function"
            " before the %s one and check that the number of internal"
            " edges is non null.\n",
            __func__,__func__);
    return 0;
  }

  if ( mesh->namax+idx > (MMG5_int)na_tot ) {
    fprintf(stderr,"\n  ## Error: %s: Can't get the internal edge of index %" MMG5_PRId "."
            " Index must be between 1 and %"MMG5_PRId".\n",
            __func__,idx,(MMG5_int)na_tot-mesh->namax);
    return 0;
  }

  ped = &mesh->edge[mesh->namax+idx];

  *e0  = ped->a;
  *e1  = ped->b;

  if ( ref != NULL ) {
    *ref = mesh->edge[mesh->namax+idx].ref;
  }

  return 1;
}

int MMG2D_Get_adjaTri(MMG5_pMesh mesh, MMG5_int kel, MMG5_int listri[3]) {

  if ( ! mesh->adja ) {
    if (! MMG2D_hashTria(mesh))
      return 0;
  }

  listri[0] = mesh->adja[3*(kel-1)+1]/3;
  listri[1] = mesh->adja[3*(kel-1)+2]/3;
  listri[2] = mesh->adja[3*(kel-1)+3]/3;

  return 1;
}

MMG5_int MMG2D_Get_adjaVertices(MMG5_pMesh mesh, MMG5_int ip, MMG5_int lispoi[MMG2D_LMAX])
{
  MMG5_int start;

  if ( !mesh->tria ) return 0;

  start=MMG2D_findTria(mesh,ip);
  if ( !start ) return 0;

  return MMG2D_Get_adjaVerticesFast(mesh,ip,start,lispoi);
}

MMG5_int MMG2D_Get_adjaVerticesFast(MMG5_pMesh mesh, MMG5_int ip,MMG5_int start, MMG5_int lispoi[MMG2D_LMAX])
{
  MMG5_pTria pt;
  int        iploc,i,i1,i2;
  MMG5_int   prevk,k,*adja,nbpoi;

  pt   = &mesh->tria[start];

  for ( iploc=0; iploc<3; ++iploc ) {
    if ( pt->v[iploc] == ip ) break;
  }

  assert(iploc!=3);

  k = start;
  i = iploc;
  nbpoi = 0;
  do {
    if ( nbpoi == MMG2D_LMAX ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent"
              " vertices of the vertex %" MMG5_PRId ":\nthe ball of point contain too many"
              " elements.\n",__func__,ip);
      return 0;
    }
    i1 = MMG5_inxt2[i];
    lispoi[nbpoi] = mesh->tria[k].v[i1];
    ++nbpoi;

    adja = &mesh->adja[3*(k-1)+1];
    prevk = k;
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];
  }
  while ( k && k != start );

  if ( k > 0 ) return nbpoi;

  /* store the last point of the boundary triangle */
  if ( nbpoi == MMG2D_LMAX ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent vertices of the"
            " vertex %" MMG5_PRId ":\nthe ball of point contain too many elements.\n",
            __func__,ip);
    return 0;
  }
  i1 = MMG5_inxt2[i1];
  lispoi[nbpoi] = mesh->tria[prevk].v[i1];
  ++nbpoi;

  /* check if boundary hit */
  k = start;
  i = iploc;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i2 = MMG5_iprv2[i];
    k  = adja[i2] / 3;
    if ( k == 0 )  break;

    if ( nbpoi == MMG2D_LMAX ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to compute adjacent vertices of the"
              " vertex %" MMG5_PRId ":\nthe ball of point contain too many elements.\n",
              __func__,ip);
      return 0;
    }
    i  = adja[i2] % 3;
    lispoi[nbpoi] = mesh->tria[k].v[i];
    ++nbpoi;

    i  = MMG5_iprv2[i];
  }
  while ( k );

  return nbpoi;
}

int MMG2D_Get_triFromEdge(MMG5_pMesh mesh, MMG5_int ked, MMG5_int *ktri, int *ied)
{
  MMG5_int val;

  val = mesh->edge[ked].base;

  if ( !val ) {
    fprintf(stderr,"  ## Error: %s: the main fonction of the Mmg library must be"
            " called before this function.\n",__func__);
    return 0;
  }

  *ktri = val/3;

  *ied = val%3;

  return 1;
}

int MMG2D_Get_trisFromEdge(MMG5_pMesh mesh, MMG5_int ked, MMG5_int ktri[2], int ied[2])
{
  int ier;
  MMG5_int itri;
#ifndef NDEBUG
  MMG5_int ia0,ib0,ia1,ib1;
#endif

  ktri[0]  =  ktri[1] = 0;
  ied[0]   =  ied[1]  = 0;

  ier = MMG2D_Get_triFromEdge(mesh, ked, ktri, ied);

  if ( !ier ) return 0;

  if ( !mesh->adja ) {
    if (!MMG2D_hashTria(mesh) )
      return 0;
  }

  itri = mesh->adja[3*(*ktri-1) + *ied + 1 ];

  if ( itri ) {
    ktri[1]  = itri/3;
    ied[1]   = itri%3;

#ifndef NDEBUG
    ia0 = mesh->tria[ktri[0]].v[MMG5_inxt2[ied[0]]];
    ib0 = mesh->tria[ktri[0]].v[MMG5_iprv2[ied[0]]];

    ia1 = mesh->tria[ktri[1]].v[MMG5_inxt2[ied[1]]];
    ib1 = mesh->tria[ktri[1]].v[MMG5_iprv2[ied[1]]];

    assert ( ( (ia0 == ia1) && (ib0 == ib1) ) ||
             ( (ia0 == ib1) && (ib0 == ia1) ) );
#endif
  }

  return 1;
}

int MMG2D_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met) {
  double      hsiz;

  /* Set solution size */
  if ( mesh->info.ani ) {
    met->size = 3;
  }
  else {
    met->size = 1;
  }

  /* Memory alloc */
  if ( !MMG2D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,met->size) )
    return 0;

  if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) )
    return 0;

  mesh->info.hsiz = hsiz;

  MMG5_Set_constantSize(mesh,met,hsiz);

  return 1;
}

int MMG2D_Compute_eigenv(double m[3],double lambda[2],double vp[2][2]) {

  return  MMG5_eigensym(m,lambda,vp);

}


void MMG2D_Reset_verticestags(MMG5_pMesh mesh) {
  MMG5_int k;

  for ( k=1; k<=mesh->np;  ++k ) {
    mesh->point[k].tag = 0;
  }

}

void MMG2D_Free_triangles(MMG5_pMesh mesh) {

  if ( mesh->adja )
    MMG5_DEL_MEM(mesh,mesh->adja);

  if ( mesh->tria )
    MMG5_DEL_MEM(mesh,mesh->tria);

  mesh->nt = 0;
  mesh->nti = 0;
  mesh->nenil = 0;

  return;
}

void MMG2D_Free_edges(MMG5_pMesh mesh) {

  if ( mesh->edge )
    MMG5_DEL_MEM(mesh,mesh->edge);

  if ( mesh->xpoint )
    MMG5_DEL_MEM(mesh,mesh->xpoint);

  mesh->na = 0;
  mesh->nai = 0;
  mesh->nanil = 0;

  mesh->xp = 0;

  return;
}

void MMG2D_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol) {

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
  sol->dim  = 2;
  sol->ver  = 2;
  sol->size = 1;
  sol->type = 1;

  return;
}
