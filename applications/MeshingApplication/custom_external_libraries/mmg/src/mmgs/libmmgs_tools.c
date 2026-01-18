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
 * \file mmgs/libmmgs_tools.c
 * \brief Tools functions for the mmgs library
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmgs.h"
#include "libmmgs_private.h"
#include "inlined_functions_private.h"
#include "mmgsexterns_private.h"
#include "mmgexterns_private.h"


void MMGS_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( (!mesh->info.ani) && ((!met) || (met->size < 6)) ) {
    MMG5_calelt      = MMG5_caltri_iso;
    MMGS_doSol       = MMGS_doSol_iso;
    MMG5_lenSurfEdg  = MMG5_lenSurfEdg_iso;
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_iso;
    MMGS_defsiz      = MMGS_defsiz_iso;
    MMGS_gradsiz     = MMG5_gradsiz_iso;
    MMGS_gradsizreq  = MMG5_gradsizreq_iso;
    intmet           = intmet_iso;
    movintpt         = movintpt_iso;
    movridpt         = movridpt_iso;
  }
  else {
    /* Force data consistency: if aniso metric is provided, met->size==6 and
     * info.ani==0; with -A option, met->size==1 and info.ani==1 */
    met->size = 6;
    mesh->info.ani = 1;

    /* Set function pointers */
    if ( (!met->m) && (!mesh->info.optim) && mesh->info.hsiz<=0. ) {
      MMG5_calelt     = MMG5_caltri_iso;
      MMGS_doSol      = MMGS_doSol_iso;
      MMG5_lenSurfEdg = MMG5_lenSurfEdg_iso;
    }
    else {
      MMG5_calelt     = MMG5_caltri_ani;
      MMGS_doSol      = MMGS_doSol_ani;
      MMG5_lenSurfEdg = MMG5_lenSurfEdg_ani;
    }
    MMG5_compute_meanMetricAtMarkedPoints = MMG5_compute_meanMetricAtMarkedPoints_ani;
    MMGS_defsiz      = MMGS_defsiz_ani;
    MMGS_gradsiz     = MMGS_gradsiz_ani;
    MMGS_gradsizreq  = MMG5_gradsizreq_ani;
    intmet        = intmet_ani;
    movintpt      = movintpt_ani;
    movridpt      = movridpt_ani;
  }
}

int MMGS_usage(char *prog) {

  /* Common generic options, file options and mode options */
  MMG5_mmgUsage(prog);

  /* Common parameters (first section) */
  MMG5_paramUsage1( );

  /* Specific options */
  fprintf(stdout,"-keep-ref    preserve initial domain references in level-set mode.\n");

#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]      Turn on or off the renumbering using SCOTCH [0/1] \n");
#endif

  fprintf(stdout,"\n");

  /* Common parameters (second section) */
  MMG5_paramUsage2();

  /* Common options for advanced users */
  MMG5_advancedUsage();

  fprintf(stdout,"\n\n");

  return 1;
}

int MMGS_defaultValues(MMG5_pMesh mesh) {

  MMG5_mmgDefaultValues(mesh);
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
int MMGS_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol sol) {
  MMG5_pSol tmp = NULL;
  double val;
  int    i,param;
  char   namein[MMG5_FILESTR_LGTH],*endptr;

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      MMGS_defaultValues(mesh);
      return 0;
    }
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) ) {
      MMGS_usage(argv[0]);
      return 0;
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case 'a': /* ridge angle */
        if ( !strcmp(argv[i],"-ar") ) {
          if ( i >= argc -1 ) {
            fprintf(stderr,"\nMissing argument option %s\n",argv[i]);
            return 0;
          }
          else {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_angleDetection,val))
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
        if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          return 0;
        break;
      case 'f':
        if ( !strcmp(argv[i],"-f") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMGS_Set_inputParamName(mesh,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing filename for %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
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
            param = MMGS_DPARAM_hmin;
          }
          else if ( !strcmp(argv[i],"-hmax") ) {
            param = MMGS_DPARAM_hmax;
          }
          else if ( !strcmp(argv[i],"-hsiz") ) {
            param = MMGS_DPARAM_hsiz;
          }
          else if ( !strcmp(argv[i],"-hausd") ) {
            param = MMGS_DPARAM_hausd;
          }
          else if ( !strcmp(argv[i],"-hgradreq") ) {
            param = MMGS_DPARAM_hgradreq;
          }
          else if ( !strcmp(argv[i],"-hgrad") ) {
            param = MMGS_DPARAM_hgrad;
          }
          else {
            /* Arg unknown by Mmg: arg starts with -h but is not known */
            MMGS_usage(argv[0]);
            return 0;
          }

          assert ( param != MMG5_UNSET );

          val = strtof(argv[i+1],&endptr);
          if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
            ++i;
            if ( !MMGS_Set_dparameter(mesh,met,param,val) ){
              return 0;
            }
          } else {
            /* argument is not a number */
            fprintf(stderr,"\nMissing argument option %s\n",argv[i]);
            return 0;
          }
        }

        break;
      case 'd':
        if ( !strcmp(argv[i],"-default") ) {
          mesh->mark=1;
        }
        else {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_debug,1) )
            return 0;
        }
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMGS_Set_inputMeshName(mesh, argv[i]) )
              return 0;

            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
              return 0;
          }else{
            fprintf(stderr,"\nMissing filname for %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-isoref") && ++i <= argc ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_isoref,
                                    atoi(argv[i])) )
            return 0;
        }
        else {
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
      case 'k':
        if ( !strcmp(argv[i],"-keep-ref") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_keepRef,1) )
            return 0;
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_iso,1) )
            return 0;

          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_ls,val))
                return 0;
            }
          }
        }
        else if ( !strcmp(argv[i],"-lssurf") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_isosurf,1) )
            return 0;

          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_ls,val))
                return 0;
            }
          }
        }
        break;
      case 'm':
        if ( !strcmp(argv[i],"-met") ) {
          if ( !met ) {
            fprintf(stderr,"\nNo metric structure allocated for %s option\n",
                    argv[i-1]);
            return 0;
          }
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMGS_Set_inputSolName(mesh,met,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing filname for %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if  (!strcmp(argv[i],"-m") ) {
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_mem,atoi(argv[i])) )
            return 0;
        }
        else {
          fprintf(stderr,"\nMissing argument option %s\n",argv[i-1]);
          MMGS_usage(argv[0]);
          return 0;
        }
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_angle,0) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nsd") ) {
          if ( ++i < argc && isdigit(argv[i][0]) ) {
            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_numsubdomain,atoi(argv[i])) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing argument option %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noswap,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noinsert,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nomove,1) )
            return 0;
        }
        else if ( !strcmp(argv[i],"-nreg") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nreg,1) )
            return 0;
        }
        else if( !strcmp(argv[i],"-nosizreq") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nosizreq,1) ) {
            return 0;
          }
        }
        break;
      case 'o':
        if ( (!strcmp(argv[i],"-out")) || (!strcmp(argv[i],"-o")) ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
              return 0;
          }else{
            fprintf(stderr,"\nMissing filname for %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        else if( !strcmp(argv[i],"-optim") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_optim,1) )
            return 0;
        }
        else {
          fprintf(stderr,"\nUnrecognized option %s\n",argv[i]);
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
      case 'r':
        if ( !strcmp(argv[i],"-rmc") ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_rmc,0) )
            return 0;
          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_rmc,val))
                return 0;
            }
          }
        }
#ifdef USE_SCOTCH
        else if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_renum,atoi(argv[i])) )
                return 0;
            }
            else {
              fprintf(stderr,"\nMissing argument option %s\n",argv[i-1]);
              MMGS_usage(argv[0]);
              return 0;
            }
          }
          else {
            fprintf(stderr,"\nMissing argument option %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
#endif
        else {
          fprintf(stderr,"\nUnrecognized option %s\n",argv[i]);
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          /* For retrocompatibility, store the metric if no sol structure available */
          tmp = sol ? sol : met;
          assert(tmp);
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMGS_Set_inputSolName(mesh,tmp,argv[i]) )
              return 0;
          }
          else {
            fprintf(stderr,"\nMissing filname for %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,atoi(argv[i])) )
              return 0;
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"\nMissing argument option %s\n",argv[i-1]);
          MMGS_usage(argv[0]);
          return 0;
        }
        break;
      case 'x':
        if ( !strcmp(argv[i],"-xreg") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_xreg,1) )
            return 0;
          if ( i < argc -1 ) {
            val = strtof(argv[i+1],&endptr);
            if ( endptr == &(argv[i+1][strlen(argv[i+1])]) ) {
              ++i;
              if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_xreg,val))
                return 0;
            }
          }
        }
        break;
      default:
        fprintf(stderr,"\nUnrecognized option %s\n",argv[i]);
        MMGS_usage(argv[0]);
        return 0;
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        if ( !MMGS_Set_inputMeshName(mesh,argv[i]) )
          return 0;
        if ( mesh->info.imprim == -99 ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
            return 0;
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
          return 0;
      }
      else {
        fprintf(stdout,"\nArgument %s ignored\n",argv[i]);
        MMGS_usage(argv[0]);
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
    if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,i) )
      return 0;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"\n  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    MMG_FSCANF(stdin,"%127s",namein);
    if ( !MMGS_Set_inputMeshName(mesh,namein) )
      return 0;
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMGS_Set_outputMeshName(mesh,"") )
      return 0;
  }

  /* adp mode: if the metric name has been stored in sol, move it in met */
  if ( met->namein==NULL && sol && sol->namein &&
       !(mesh->info.iso || mesh->info.isosurf || mesh->info.lag>=0) ) {
    if ( !MMGS_Set_inputSolName(mesh,met,sol->namein) )
      return 0;
    MMG5_DEL_MEM(mesh,sol->namein);
  }

  /* default : store solution name in iso mode, metric name otherwise */
  tmp = ( mesh->info.iso || mesh->info.isosurf || mesh->info.lag >=0 ) ? sol : met;
  assert ( tmp );
  if ( tmp->namein == NULL ) {
    if ( !MMGS_Set_inputSolName(mesh,tmp,"") ) { return 0; }
  }
  if ( met->nameout == NULL ) {
    if ( !MMGS_Set_outputSolName(mesh,met,"") )
      return 0;
  }
  return 1;
}

int MMGS_freeLocalPar(MMG5_pMesh mesh) {

  free(mesh->info.par);
  mesh->info.npar = 0;

  return 1;
}

int MMGS_stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(&mesh->info,info,sizeof(MMG5_Info));
  MMGS_memOption(mesh);
  if( mesh->info.mem > 0) {
    if ( mesh->npmax < mesh->np || mesh->ntmax < mesh->nt ) {
      return 0;
    } else if(mesh->info.mem < 39)
      return 0;
  }
  return 1;
}

void MMGS_destockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(info,&mesh->info,sizeof(MMG5_Info));
  return;
}

int MMGS_Get_numberOfNonBdyEdges(MMG5_pMesh mesh, MMG5_int* nb_edges) {
  MMG5_pTria pt,pt1;
  MMG5_pEdge ped;
  MMG5_int   *adja,k,j,iel;
  int        i,i1,i2;

  *nb_edges = 0;
  if ( mesh->tria ) {
    /* Create the triangle adjacency if needed */
    if ( !mesh->adja ) {
      if ( !MMGS_hashTria(mesh) ) {
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
        /* Do not treat boundary edges */
        if ( MG_EDG(pt->tag[i]) ) continue;

        iel = adja[i] / 3;
        assert ( iel != k );

        pt1 = &mesh->tria[iel];

        if ( (!iel) || (pt->ref != pt1->ref) ) {
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
    if ( mesh->na ) {
      MMG5_ADD_MEM(mesh,(*nb_edges)*sizeof(MMG5_Edge),"non boundary edges",
                   printf("  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->edge,(mesh->na+1),(mesh->na+(*nb_edges)+1),
                         MMG5_Edge,"non bdy edges arrray",return 0);
    }
    else {
      MMG5_ADD_MEM(mesh,((*nb_edges)+1)*sizeof(MMG5_Edge),"non boundary edges",
                   printf("  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->edge,0,(*nb_edges)+1,
                         MMG5_Edge,"non bdy edges arrray",return 0);
    }

    j = mesh->na+1;
    for ( k=1; k<=mesh->nt; k++ ) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      adja = &mesh->adja[3*(k-1)+1];

      for ( i=0; i<3; i++ ) {
        /* Do not treat boundary edges */
        if ( MG_EDG(pt->tag[i]) ) continue;

        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        iel = adja[i] / 3;
        assert ( iel != k );

        pt1 = &mesh->tria[iel];

        if ( (!iel) || (pt->ref != pt1->ref) ) {
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

int MMGS_Get_nonBdyEdge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref, MMG5_int idx) {
  MMG5_pEdge ped;
  size_t     na_tot=0;
  char       *ptr_c = (char*)mesh->edge;

  if ( !mesh->edge ) {
    fprintf(stderr,"\n  ## Error: %s: edge array is not allocated.\n"
            " Please, call the MMGS_Get_numberOfNonBdyEdges function"
            " before the %s one.\n",
            __func__,__func__);
    return 0;
  }

  ptr_c = ptr_c-sizeof(size_t);
  na_tot = *((size_t*)ptr_c);

  if ( (size_t)mesh->namax==na_tot ) {
    fprintf(stderr,"\n  ## Error: %s: no internal edge.\n"
            " Please, call the MMGS_Get_numberOfNonBdyEdges function"
            " before the %s one and check that the number of internal"
            " edges is non null.\n",
            __func__,__func__);
  }

  if ( (size_t)mesh->namax+idx > na_tot ) {
    fprintf(stderr,"\n  ## Error: %s: Can't get the internal edge of index %" MMG5_PRId "."
            " Index must be between 1 and %zu.\n",
            __func__,idx,na_tot-(size_t)mesh->namax);
  }

  ped = &mesh->edge[mesh->na+idx];

  *e0  = ped->a;
  *e1  = ped->b;

  if ( ref != NULL ) {
    *ref = mesh->edge[mesh->na+idx].ref;
  }

  return 1;
}

int MMGS_Get_adjaTri(MMG5_pMesh mesh, MMG5_int kel, MMG5_int listri[3]) {

  if ( ! mesh->adja ) {
    if (! MMGS_hashTria(mesh))
      return 0;
  }

  listri[0] = mesh->adja[3*(kel-1)+1]/3;
  listri[1] = mesh->adja[3*(kel-1)+2]/3;
  listri[2] = mesh->adja[3*(kel-1)+3]/3;

  return 1;
}

int MMGS_Get_adjaVerticesFast(MMG5_pMesh mesh, MMG5_int ip,MMG5_int start, MMG5_int lispoi[MMGS_LMAX])
{
  MMG5_pTria pt;
  MMG5_int   k,prevk,*adja;
  int        i,i1,i2,iploc,nbpoi;

  pt   = &mesh->tria[start];

  for ( iploc=0; iploc<3; ++iploc ) {
    if ( pt->v[iploc] == ip ) break;
  }

  assert(iploc!=3);

  k = start;
  i = iploc;
  nbpoi = 0;
  do {
    if ( nbpoi == MMGS_LMAX ) {
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
  if ( nbpoi == MMGS_LMAX ) {
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

    if ( nbpoi == MMGS_LMAX ) {
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

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the solution structure.
 * \param ani 1 for aniso metric, 0 for iso one
 *
 * \return 0 if fail, 1 if succeed.
 *
 * Truncate the metric computed by the DoSol function by hmax and hmin values
 * (if setted by the user). Set hmin and hmax if they are not setted.
 *
 */
static inline
int MMGS_solTruncatureForOptim(MMG5_pMesh mesh, MMG5_pSol met,int ani) {
  MMG5_int k;
  int      i,ier;

  assert ( mesh->info.optim );

  /* Detect points not used by the mesh */
  ++mesh->base;

#ifndef NDEBUG
  for (k=1; k<=mesh->np; k++) {
    assert ( mesh->point[k].flag < mesh->base );
  }
#endif

  for (k=1; k<=mesh->nt; k++) {
    MMG5_pTria ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) ) continue;

    for (i=0; i<3; i++) {
      mesh->point[ptt->v[i]].flag = mesh->base;
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
 * \param mesh pointer to the mesh
 * \param met pointer to the metric
 *
 * \return 1 if succeed, 0 if fail
 *
 * \remark need the normal at vertices.
 *
 * Compute isotropic size map according to the mean of the length of the
 * edges passing through a point.
 *
 */
int MMGS_doSol_iso(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria   ptt;
  MMG5_pPoint  p1,p2;
  double       ux,uy,uz,dd;
  MMG5_int     k,ipa,ipb;
  int          i,type;
  //We guess that we have less than INT32_MAX edges
  // passing through each point
  int          *mark;

  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  /* Memory alloc */
  if ( met->size!=1 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
  }

  type = 1;
  if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  /* Travel the triangles edges and add the edge contribution to edges
   * extermities */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) )  continue;

    for (i=0; i<3; i++) {
      ipa = ptt->v[i];
      ipb = ptt->v[MMG5_inxt2[i]];
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

  MMGS_solTruncatureForOptim(mesh,met,0);

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param k index of starting triangle
 * \param i local index of point \a p1 in \a k
 * \param p1 point on which we want to compute the 3D unit tensor
 * \param m pointer to store computed metric
 *
 * \return 1 if succeed, 0 if the tensor of covariance is non invertible.
 *
 * Compute the 3D unit tensor at point \a p1, the vertex number \a i of tetra \a
 * k.
 *
 * \remark If all edges starting from \a p1 belong to the same plane, the tensor
 * of covariance is not invertible and \a m can't be computed.
 *
 */
static inline
int MMGS_unitTensor_3D( MMG5_pMesh mesh,MMG5_int k,int i,MMG5_pPoint p1,double *m) {
  MMG5_int    list[MMGS_LMAX+2];
  int         ilist,j;
  int8_t      open;

  /** Step 1: compute ball of point */
  ilist = MMG5_boulet(mesh,k,i,list,1,&open);
  if ( ilist < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to compute ball of point.\n",
            __func__);
    return 0;
  }

  /** Step 2: compute unit tensor */
  if ( (!MG_SIN(p1->tag)) && (MG_GEO & p1->tag) && open ) {
    /* MG_GEO point along an open boundary: we want to compute the 2D unit
     * tensor */
    return 0;
  }

  for ( j=0; j<6; ++j ) {
    m[j] = 0.;
  }

  for ( j=0; j<ilist; ++j ) {
    /* Compute euclidean edge length */
    double u[3];
    MMG5_int iel = list[j] / 3;
    int8_t   i1  = list[j] % 3;
    int8_t   i2  = MMG5_inxt2[i1];

    MMG5_pTria ptt = &mesh->tria[iel];
    MMG5_pPoint p2 = &mesh->point[ptt->v[i2]];

    u[0]  = p1->c[0] - p2->c[0];
    u[1]  = p1->c[1] - p2->c[1];
    u[2]  = p1->c[2] - p2->c[2];

    m[0] += u[0]*u[0];
    m[1] += u[0]*u[1];
    m[2] += u[0]*u[2];
    m[3] += u[1]*u[1];
    m[4] += u[1]*u[2];
    m[5] += u[2]*u[2];
  }

  double tensordot[6];
  int ier = MMG5_invmat(m,tensordot);

  /* Invmat make the assumtion that input matrix is invertible (always the case
   * normally). Here, if all edges belongs to the same plane, it is not the case
   * so we have to check the values of the resulting matrix even if invmat
   * succeed */
  if ( ier ) {
    for ( j=0; j<6; ++j ) {
      if ( !isfinite(tensordot[j]) ) {
        ier = 0;
        break;
      }
    }
  }

  if ( ier ) {
    double lambda[3],vp[3][3];
    if ( !MMG5_eigenv3d(1,tensordot,lambda,vp) ) {
      ier = 0;
    }
    else if ( (!isfinite(lambda[0])) || (!isfinite(lambda[1])) || (!isfinite(lambda[2])) ) {
      ier = 0;
    }
    else if ( lambda[0]<=0 || lambda[1]<=0 || lambda[2]<=0  ) {
      ier = 0;
    }
  }

  /* Non invertible matrix: put FLT_MIN waiting for a treatment later */
  if ( !ier ) {
    m[0] = FLT_MIN;
    m[1] = 0;
    m[2] = 0;
    m[3] = FLT_MIN;
    m[4] = 0;
    m[5] = FLT_MIN;
    return 0;
  }

  double dd = (double)ilist/3.;
  for ( j=0; j<6; ++j ) {
    m[j] = dd*tensordot[j];
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param p0 starting point
 * \param k index of starting element
 * \param i local index of \a p0 in \a k
 * \param ilist computed number of tria in the ball of \a p0
 * \param r rotation that send the normal at p0 onto the z vector
 * \param lipoint rotated ball of point \a p0
 * \param n normal at point \a p0
 *
 * \return 1 if success, 0 otherwise.
 *
 * Compute the rotation matrix that sends the tangent plane at \a p0 onto z=0
 * and apply this rotation to the opened ball of \a p0.
 *
 */
static inline
int MMGS_surfopenballRotation(MMG5_pMesh mesh,MMG5_pPoint p0,MMG5_int k, int i,
                              int ilist,double r[3][3],double *lispoi,double n[3]) {
  MMG5_pTria  ptt;
  MMG5_pPoint p2;
  double      ux,uy,uz;
  MMG5_int    *adja,kold;
  int         jold,i1;

  if ( !MMG5_rotmatrix(n,r) ) {
    return 0;
  }

  /* Enumeration of the ball starting from a boundary */
  MMG5_int iel = k;
  int j        = i;
  do {
    adja  = &mesh->adja[3*(iel-1)+1];
    i1    = MMG5_iprv2[j];
    kold  = iel;
    jold  = j;
    iel   = adja[i1] / 3;
    j     = (int)adja[i1] % 3;
    j     = MMG5_iprv2[j];
  }
  // Remark: here the test k!=start is a security bound: theorically it is
  // useless but in case of bad edge tag, it ensure that the loop is not
  // infinite.
  while ( iel && iel != k );

  iel  = kold;
  j    = jold;
  int idx = 0;
  do {
    ptt       = &mesh->tria[iel];
    // here, p1 == &mesh->point[ptt->v[j]]

    adja      = &mesh->adja[3*(iel-1)+1];
    i1        = MMG5_inxt2[j];
    p2 = &mesh->point[ptt->v[i1]];

    ux = p2->c[0] - p0->c[0];
    uy = p2->c[1] - p0->c[1];
    uz = p2->c[2] - p0->c[2];

    lispoi[3*idx+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    lispoi[3*idx+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
    lispoi[3*idx+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

    iel       = adja[i1] / 3;
    j         = (int)adja[i1] % 3;
    j         = MMG5_inxt2[j];

    ++idx;
  }
  while ( iel );

  /* last point : the half-ball is open : ilist tria, and ilist +1 points ;
     lists are enumerated in direct order */
  i1 = MMG5_inxt2[i1];
  p2 = &mesh->point[ptt->v[i1]];

  ux = p2->c[0] - p0->c[0];
  uy = p2->c[1] - p0->c[1];
  uz = p2->c[2] - p0->c[2];

  lispoi[3*idx+1] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
  lispoi[3*idx+2] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;
  lispoi[3*idx+3] =  r[2][0]*ux + r[2][1]*uy + r[2][2]*uz;

  assert ( idx == ilist );

  /* Check all projections over tangent plane. */
  for (j=0; j<ilist; j++) {
    double area = lispoi[3*j+1]*lispoi[3*(j+1)+2] - lispoi[3*j+2]*lispoi[3*(j+1)+1];
    if ( area <= 0.0 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute ball rotation.\n",
              __func__);
      return 0;
    }
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param k index of starting triangle
 * \param i local index of point \a p1 in \a k
 * \param p1 point on which we want to compute the 3D unit tensor
 * \param m pointer to store computed metric
 * \param isqhmax squared inverse of MMG5_HMAXCOE
 *
 * \return 1 if succeed, 0 if fail.
 *
 * Compute the 2D unit tensor at point \a p1, the vertex number \a i of tetra \a
 * k : edges of the ball of point are projected onto the tangent plane.
 *
 */
static inline
int MMGS_unitTensor_2D ( MMG5_pMesh mesh,MMG5_int k,int i,MMG5_pPoint p1,
                         double *m,double isqhmax) {
  MMG5_pTria  ptt;
  double      r[3][3],lispoi[3*MMGS_LMAX+1],b0[3],b1[3],b2[3],dd;
  MMG5_int    list[MMGS_LMAX+2];
  int         ilist,j;
  int8_t      opn;

  /** Step 1: compute ball of point */

  /* Possible improvement: if we have called MMGS_unitTensor_3D previously,
   * boulet is already computed */
  ilist = MMG5_boulet(mesh,k,i,list,1,&opn);
  if ( ilist < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to compute ball of point.\n",
            __func__);
    return 0;
  }


  /** Step 2: Store or compute the suitable normal depending on the type of
   * point we are processing. */
  /* If point \a p1 is corner or required we compute the normal
  at point as the * mean of the normals at triangles of the ball. If the point
  is ridge (can be * non-manifold), we compute it as the mean of n1 and n2. In
  the other cases, * the normal is stored in p1->n */
  double nn[3] = {0.,0.,0.};
  if ( MG_SIN(p1->tag) ) {
    /* Corner or required point: compute the normal at point as the mean of
     * normals at triangles */
    for ( j=0; j<ilist; ++j ) {
      ptt = &mesh->tria[list[j]/3];
      double n[3];

      MMG5_nortri(mesh,ptt,n);
      nn[0] += n[0];  nn[1] += n[1];  nn[2] += n[2];
    }
    /* normalization */
    dd = nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2];
    if ( dd > MMG5_EPSD2 ) {
      dd = 1.0 / sqrt(dd);
      nn[0] *= dd;
      nn[1] *= dd;
      nn[2] *= dd;
    }
  }
  else if ( MG_GEO & p1->tag ) {
    if ( MG_NOM & p1->tag ) {
      fprintf(stderr,"\n  ## Error: %s: we should not pass here with a "
              "non-manifold point: it sould always be posible to compute the 3D"
              " unit tensor on such points.\n",
              __func__);
      return 0;
    }
    /* Manifold ridge point: compute the normal as the mean of the 2 computed
     * normals */
    assert ( p1->xp );
    nn[0] = mesh->xpoint[p1->xp].n1[0] + mesh->xpoint[p1->xp].n2[0];
    nn[1] = mesh->xpoint[p1->xp].n1[1] + mesh->xpoint[p1->xp].n2[1];
    nn[2] = mesh->xpoint[p1->xp].n1[2] + mesh->xpoint[p1->xp].n2[2];

    /* normalization */
    dd = nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2];
    if ( dd > MMG5_EPSD2 ) {
      dd = 1.0 / sqrt(dd);
      nn[0] *= dd;
      nn[1] *= dd;
      nn[2] *= dd;
    }
  }
  else if ( MG_REF & p1->tag ) {
    /* Reference point: normal is stored in the xpoint */
    assert ( p1->xp );
    nn[0] = mesh->xpoint[p1->xp].n1[0];
    nn[1] = mesh->xpoint[p1->xp].n1[1];
    nn[2] = mesh->xpoint[p1->xp].n1[2];
  }
  else {
    /* Regular point: normal is stored in the field n of the point */
    nn[0] = p1->n[0]; nn[1] = p1->n[1]; nn[2] = p1->n[2];
  }

  /** Step 3: Rotation of the ball of p0 so lispoi will contain all the points
     of the ball of p0, rotated so that t_{p_0}S = [z = 0] */
  if ( opn ) {
    if ( !MMGS_surfopenballRotation(mesh,p1,k,i,ilist,r,lispoi,nn)  ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute opened ball rotation.\n",
              __func__);
      return 0;
    }
  }
  else {
    if ( !MMGS_surfballRotation(mesh,p1,list,ilist,r,lispoi,nn)  ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute ball rotation.\n",
              __func__);
      return 0;
    }
  }

  /** Step 4: computation of unit tensor */
  double tensordot[3];
  tensordot[0] = 0.;
  tensordot[1] = 0.;
  tensordot[2] = 0.;

  for ( j=0; j<ilist; ++j ) {
    /* Compute the 2D metric tensor from the projection of the vectors on
     * the tangent plane */
    tensordot[0] += lispoi[3*j+1]*lispoi[3*j+1];
    tensordot[1] += lispoi[3*j+1]*lispoi[3*j+2];
    tensordot[2] += lispoi[3*j+2]*lispoi[3*j+2];
  }

  /* Metric = nedges/dim * inv (sum(tensor_dot(edges,edges)))  */
  dd = 1./(tensordot[0]*tensordot[2] - tensordot[1]*tensordot[1]);
  dd *= (double)ilist/2.;

  double tmp = tensordot[0];

  tensordot[0] = dd*tensordot[2];
  tensordot[1] = -dd*tensordot[1];
  tensordot[2] = dd*tmp;

#ifndef NDEBUG
  /* Check 2D metric */
  double lambda[2],vp[2][2];
  MMG5_eigensym(tensordot,lambda,vp);

  assert (lambda[0] > 0. && lambda[1] > 0. && "Negative eigenvalue");

  /* Normally the case where the point belongs to only 2 colinear points is
     impossible */
  assert (isfinite(lambda[0]) && isfinite(lambda[1]) && "wrong eigenvalue");
#endif

  /* At this point, tensordot (with 0 replaced by isqhmax in the z
     direction) is the desired metric, except it is expressed in the rotated
     canonical basis, that is tensordot = R * metric in cb * ^t R, so metric
     in cb = ^tR*intm*R */
  // intm = intm[0]  intm[1]    0
  //        intm[1]  intm[2]    0
  //           0       0     isqhmax

  /* b0 and b1 serve now for nothing : let them be the lines of matrix intm*R */

  // Remark: here, we put a 'fake' value along the normal direction but we can't
  // use a nan or inf value because it impacts the computation of the rotation
  // for the other directions
  b0[0] = tensordot[0]*r[0][0] + tensordot[1]*r[1][0];
  b0[1] = tensordot[0]*r[0][1] + tensordot[1]*r[1][1];
  b0[2] = tensordot[0]*r[0][2] + tensordot[1]*r[1][2];
  b1[0] = tensordot[1]*r[0][0] + tensordot[2]*r[1][0];
  b1[1] = tensordot[1]*r[0][1] + tensordot[2]*r[1][1];
  b1[2] = tensordot[1]*r[0][2] + tensordot[2]*r[1][2];

  b2[0] = isqhmax*r[2][0];
  b2[1] = isqhmax*r[2][1];
  b2[2] = isqhmax*r[2][2];

  m[0] = r[0][0] * b0[0] + r[1][0] * b1[0] + r[2][0] * b2[0];
  m[1] = r[0][0] * b0[1] + r[1][0] * b1[1] + r[2][0] * b2[1];
  m[2] = r[0][0] * b0[2] + r[1][0] * b1[2] + r[2][0] * b2[2];

  m[3] = r[0][1] * b0[1] + r[1][1] * b1[1] + r[2][1] * b2[1];
  m[4] = r[0][1] * b0[2] + r[1][1] * b1[2] + r[2][1] * b2[2];

  m[5] = r[0][2] * b0[2] + r[1][2] * b1[2] + r[2][2] * b2[2];

  return 1;
}



/**
 * \param mesh pointer to the mesh
 * \param met pointer to the metric
 *
 * \return 1 if succeed, 0 if fail
 *
 * \remark need the normal at vertices.
 * \remark hmax/hmin may be not setted here so we can't use it.
 *
 * Compute the unit metric tensor at mesh vertices (from edges passing through
 * the vertices).
 *
 */
int MMGS_doSol_ani(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria   ptt;
  MMG5_pPoint  p1;
  double       *m,tensordot[6];
  MMG5_int     k,iadr,ip;
  int          i,j,type;

  /* Memory alloc */
  if ( met->size!=6 ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected size of metric: %d.\n",
            __func__,met->size);
    return 0;
  }

  type = 3;
  if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate metric.\n",
            __func__);
    return 0;
  }

  /** Increment base marker to detect points already processed */
  ++mesh->base;

#ifndef NDEBUG
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    assert ( p1->flag < mesh->base );
  }
#endif

  double isqhmax = 1./(MMG5_HMAXCOE*MMG5_HMAXCOE);

  /** Step 1: treat non-manifold points: travel triangle edges and add edge
   * contribution to extremities (we have to do that because ball of
   * non-manifold points can't be computed) */

  /* mark array will be used only for non-manifold points */
  int     *mark;
  MMG5_SAFE_CALLOC(mark,mesh->np+1,int,return 0);

  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) )  continue;

    for (i=0; i<3; i++) {
      MMG5_int ipa = ptt->v[MMG5_iprv2[i]];
      MMG5_int ipb = ptt->v[MMG5_inxt2[i]];
      p1  = &mesh->point[ipa];
      MMG5_pPoint p2  = &mesh->point[ipb];

      if ( (!(MG_NOM & p1->tag)) && !(MG_NOM & p2->tag) ) {
        continue;
      }

      double u[3];
      u[0]  = p1->c[0] - p2->c[0];
      u[1]  = p1->c[1] - p2->c[1];
      u[2]  = p1->c[2] - p2->c[2];

      tensordot[0] = u[0]*u[0];
      tensordot[1] = u[0]*u[1];
      tensordot[2] = u[0]*u[2];
      tensordot[3] = u[1]*u[1];
      tensordot[4] = u[1]*u[2];
      tensordot[5] = u[2]*u[2];

      if ( MG_NOM & p1->tag ) {
        iadr = 6*ipa;
        for ( j=0; j<6; ++j ) {
          met->m[iadr+j]   += tensordot[j];
        }
        mark[ipa]++;
      }

      if ( MG_NOM & p2->tag ) {
        iadr = 6*ipb;
        for ( j=0; j<6; ++j ) {
          met->m[iadr+j]   += tensordot[j];
        }
        mark[ipb]++;
      }
    }
  }

  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( !mark[k] ) {
      continue;
    }

    p1->flag = mesh->base;

    /* Metric = nedges/dim * inv (sum(tensor_dot(edges,edges))).
     * sum(tensor_dot) is stored in sol->m so reuse tensordot to
     * compute M.  */
    iadr = 6*k;
    if ( !MMG5_invmat(met->m+iadr,tensordot) ) {
      /* Non invertible matrix: impose isqhmax, it will be truncated by hmax
       * later */
      fprintf(stdout, " ## Warning: %s: %d: non invertible matrix."
             " Impose hmax size at point\n",__func__,__LINE__);
      met->m[iadr+0] = isqhmax;
      met->m[iadr+2] = 0;
      met->m[iadr+3] = isqhmax;
      met->m[iadr+4] = 0;
      met->m[iadr+5] = isqhmax;
      continue;
    }

    double dd = (double)mark[k]/3.;

    for ( j=0; j<6; ++j ) {
      met->m[iadr+j] = dd*tensordot[j];
    }
  }

  MMG5_SAFE_FREE(mark);


  /* Step 2: travel the triangles to treat the other type of points: On corner
   * or ridge point, try to compute 3D unit tensor, on other points, projection
   * of the bal of point onto the tangent plane and computation of the 2D unit
   * tensor. */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) )  continue;

    for (i=0; i<3; i++) {
      ip = ptt->v[i];

      p1   = &mesh->point[ip];
      if ( p1->flag == mesh->base ) {
        continue;
      }

      iadr = 6*ip;
      m    = &met->m[iadr];

      if ( MG_CRN & p1->tag ) {
        /** Corner point (no normal at vertex): if the corner defines an angle
         * we can compute the 3D unit tensor (non singular), in the other cases
         * (corner along a flat surface, for example because it is at the
         * intersection of 3 specific edges or because it is provided by the
         * user) we will project the edges onto the tangent plane and compute
         * the 2D unit tensor (as for required and regular points) */
        int ier = MMGS_unitTensor_3D( mesh, k, i, p1, m);

        if ( !ier ) {
          /* Non invertible matrix: compute the 2D unit tensor */
          ier = MMGS_unitTensor_2D( mesh, k, i, p1, m, isqhmax);
        }

        if ( !ier ) {
          fprintf(stderr,"\n  ## Error: %s: unable to compute anisotropic unit"
                  " tensor at corner point %"MMG5_PRId".\n",__func__,MMGS_indPt(mesh,ip));
          return 0;
        }
        p1->flag = mesh->base;
      }
      else if ( MG_REQ & p1->tag ) {
        /** Required point (no normal at vertex): projection of edges onto the
         * tangent plane and computation of the 2D unit tensor */
        if ( ! MMGS_unitTensor_2D( mesh, k, i, p1, m, isqhmax) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to compute anisotropic unit"
                  " tensor at required point %"MMG5_PRId".\n",__func__,MMGS_indPt(mesh,ip));
          return 0;
        }
        p1->flag = mesh->base;
      }
      else if ( p1->tag & MG_GEO && !(p1->tag & MG_NOM) ) {
        /** Ridge point (2 normals): normally we can compute the 3D unit
         * tensor. If the angle is too flat, the computation fails and we try to
         * compute the 2D unit tensor on the tangent plane of the point. */
        int ier = MMGS_unitTensor_3D( mesh, k, i, p1, m);

        if ( !ier ) {
          /* Non invertible matrix: compute the 2D unit tensor */
          ier = MMGS_unitTensor_2D( mesh, k, i, p1, m, isqhmax);
        }

        if ( !ier ) {
          fprintf(stderr,"\n  ## Error: %s: unable to compute anisotropic unit"
                  " tensor at ridge point %"MMG5_PRId".\n",__func__,
                  MMGS_indPt(mesh,ip));
        }
        p1->flag = mesh->base;
      }
      else {
        /** Regular or reference point: projection of edges onto the tangent
         * plane and computation of the 2D unit metric tensor */

        assert ( !(MG_NOM & p1->tag) ); // non-manifold points should have been already treated

        if ( ! MMGS_unitTensor_2D( mesh, k, i, p1, m, isqhmax) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to compute anisotropic unit"
                  " tensor at required point %"MMG5_PRId".\n",__func__,
                  MMGS_indPt(mesh,ip));
          return 0;
        }
        p1->flag = mesh->base;
      }
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
    }
  }

  /** Step 3: check computed metric */
#ifndef NDEBUG
  for (k=1; k<=mesh->np; k++) {
    p1 = &mesh->point[k];
    if ( p1->flag != mesh->base ) {
      continue;
    }
    iadr = 6*k;

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
  }
#endif

  /** Step 4: computation of hmin/hmax (if needed) and metric truncature */
  MMGS_solTruncatureForOptim(mesh,met,1);

  return 1;
}

int MMGS_Set_constantSize(MMG5_pMesh mesh,MMG5_pSol met) {
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
  if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type) )
    return 0;

  if ( !MMG5_Compute_constantSize(mesh,met,&hsiz) )
    return 0;

  mesh->info.hsiz = hsiz;

  MMG5_Set_constantSize(mesh,met,hsiz);

  return 1;
}

int MMGS_Compute_eigenv(double m[6],double lambda[3],double vp[3][3]) {

  return  MMG5_eigenv3d(1,m,lambda,vp);

}

void MMGS_Free_solutions(MMG5_pMesh mesh,MMG5_pSol sol) {

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

int MMGS_Clean_isoSurf(MMG5_pMesh mesh) {

  return MMG5_Clean_isoEdges(mesh);
}
