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

#include "mmgs.h"
#include "inlined_functions.h"

void MMGS_setfunc(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( met->size < 6 ) {
    _MMG5_calelt      = _MMG5_caltri_iso;
    _MMG5_lenSurfEdg  = _MMG5_lenSurfEdg_iso;
    _MMG5_defsiz      = _MMGS_defsiz_iso;
    gradsiz           = gradsiz_iso;
    intmet            = intmet_iso;
    movintpt          = movintpt_iso;
    movridpt          = movridpt_iso;
  }
  else {
    if ( !met->m ) {
      _MMG5_calelt     = _MMG5_caltri_iso;
      _MMG5_lenSurfEdg = _MMG5_lenSurfEdg_iso;
    }
    else {
      _MMG5_calelt     = _MMG5_caltri_ani;
      _MMG5_lenSurfEdg = _MMG5_lenSurfEdg_ani;
    }
    _MMG5_defsiz  = _MMGS_defsiz_ani;
    gradsiz       = gradsiz_ani;
    intmet        = intmet_ani;
    movintpt      = movintpt_ani;
    movridpt      = movridpt_ani;
  }
}

void MMGS_usage(char *prog) {
  _MMG5_mmgUsage(prog);

  fprintf(stdout,"-A           enable anisotropy (without metric file).\n");
  fprintf(stdout,"-keep-ref    preserve initial domain references in level-set mode.\n");
  fprintf(stdout,"-nreg        normal regul.\n");
#ifdef USE_SCOTCH
  fprintf(stdout,"-rn [n]      Turn on or off the renumbering using SCOTCH [0/1] \n");
#endif
  fprintf(stdout,"\n\n");

  exit(EXIT_FAILURE);
}

void MMGS_defaultValues(MMG5_pMesh mesh) {

  _MMG5_mmgDefaultValues(mesh);
#ifdef USE_SCOTCH
  fprintf(stdout,"SCOTCH renumbering                  : enabled\n");
#else
  fprintf(stdout,"SCOTCH renumbering                  : disabled\n");
#endif
  fprintf(stdout,"\n\n");

  exit(EXIT_FAILURE);
}

int MMGS_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  int    i;
  char   namein[128];

  /* First step: search if user want to see the default parameters values. */
  for ( i=1; i< argc; ++i ) {
    if ( !strcmp(argv[i],"-val") ) {
      MMGS_defaultValues(mesh);
    }
  }

  /* Second step: read all other arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        MMGS_usage(argv[0]);
        break;
      case 'a': /* ridge angle */
        if ( !strcmp(argv[i],"-ar") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_angleDetection,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'A': /* anisotropy */
        if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Tensor) )
          exit(EXIT_FAILURE);
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hmin,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hmax,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hausd") && ++i <= argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hausd,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-hgrad") && ++i <= argc ) {
          if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_hgrad,
                                    atof(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else
          MMGS_usage(argv[0]);
        break;
      case 'd':
        if ( !strcmp(argv[i],"-default") ) {
          mesh->mark=1;
        }
        else {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_debug,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMGS_Set_inputMeshName(mesh, argv[i]) )
              exit(EXIT_FAILURE);

            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            MMGS_usage(argv[0]);
          }
        }
        break;
      case 'k':
        if ( !strcmp(argv[i],"-keep-ref") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_keepRef,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'l':
        if ( !strcmp(argv[i],"-ls") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_iso,1) )
            exit(EXIT_FAILURE);
          if ( ++i < argc && (isdigit(argv[i][0]) ||
                              (argv[i][0]=='-' && isdigit(argv[i][1])) ) ) {
            if ( !MMGS_Set_dparameter(mesh,met,MMGS_DPARAM_ls,atof(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else i--;
        }
        break;
      case 'm':
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_mem,atoi(argv[i])) )
            exit(EXIT_FAILURE);
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMGS_usage(argv[0]);
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-nr") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_angle,0) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-noswap") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noswap,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-noinsert") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_noinsert,1) )
            exit(EXIT_FAILURE);
        }
        else if( !strcmp(argv[i],"-nomove") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nomove,1) )
            exit(EXIT_FAILURE);
        }
        else if ( !strcmp(argv[i],"-nreg") ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_nreg,1) )
            exit(EXIT_FAILURE);
        }
        break;
      case 'o':
        if ( !strcmp(argv[i],"-out") ) {
          if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
            if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
              exit(EXIT_FAILURE);
          }else{
            fprintf(stderr,"Missing filname for %c%c%c\n",
                    argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMGS_usage(argv[0]);
          }
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_renum,atoi(argv[i])) )
                exit(EXIT_FAILURE);
            }
            else {
              fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
              MMGS_usage(argv[0]);
            }
          }
          else {
            fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
            MMGS_usage(argv[0]);
          }
        }
        break;
#endif
      case 's':
        if ( !strcmp(argv[i],"-sol") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( !MMGS_Set_inputSolName(mesh,met,argv[i]) )
              exit(EXIT_FAILURE);
          }
          else {
            fprintf(stderr,"Missing filname for %c%c%c\n",argv[i-1][1],argv[i-1][2],argv[i-1][3]);
            MMGS_usage(argv[0]);
          }
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( argv[i][0] == '-' || isdigit(argv[i][0]) ) {
            if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,atoi(argv[i])) )
              exit(EXIT_FAILURE);
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          MMGS_usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        MMGS_usage(argv[0]);
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        if ( !MMGS_Set_inputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
        if ( mesh->info.imprim == -99 ) {
          if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,5) )
            exit(EXIT_FAILURE);
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMGS_Set_outputMeshName(mesh,argv[i]) )
          exit(EXIT_FAILURE);
      }
      else {
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        MMGS_usage(argv[0]);
      }
    }
    i++;
  }

  /* check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    if ( !MMGS_Set_iparameter(mesh,met,MMGS_IPARAM_verbose,i) )
      exit(EXIT_FAILURE);
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%s",namein);
    if ( !MMGS_Set_inputMeshName(mesh,namein) )
      exit(EXIT_FAILURE);
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMGS_Set_outputMeshName(mesh,"") )
      exit(EXIT_FAILURE);
  }

  if ( met->namein == NULL ) {
    if ( !MMGS_Set_inputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  if ( met->nameout == NULL ) {
    if ( !MMGS_Set_outputSolName(mesh,met,"") )
      exit(EXIT_FAILURE);
  }
  return(1);
}

int MMGS_stockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(&mesh->info,info,sizeof(MMG5_Info));
  _MMGS_memOption(mesh);
  if( mesh->info.mem > 0) {
    if ( mesh->npmax < mesh->np || mesh->ntmax < mesh->nt ) {
      return(0);
    } else if(mesh->info.mem < 39)
      return(0);
  }
  return(1);
}

void MMGS_destockOptions(MMG5_pMesh mesh, MMG5_Info *info) {

  memcpy(info,&mesh->info,sizeof(MMG5_Info));
  return;
}

int MMGS_Get_adjaTri(MMG5_pMesh mesh, int kel, int listri[3]) {

  if ( ! mesh->adja ) {
    if (! _MMGS_hashTria(mesh))
      return(0);
  }

  listri[0] = mesh->adja[3*(kel-1)+1]/3;
  listri[1] = mesh->adja[3*(kel-1)+2]/3;
  listri[2] = mesh->adja[3*(kel-1)+3]/3;

  return(1);
}

inline
int MMGS_Get_adjaVerticesFast(MMG5_pMesh mesh, int ip,int start, int lispoi[MMGS_LMAX])
{
  MMG5_pTria pt;
  int k,prevk,nbpoi,iploc,i,i1,i2,*adja;

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
      fprintf(stderr,"  ## Warning: unable to compute adjacent vertices of the"
              " vertex %d:\nthe ball of point contain too many elements.\n",ip);
      return(0);
    }
    i1 = _MMG5_inxt2[i];
    lispoi[nbpoi] = mesh->tria[k].v[i1];
    ++nbpoi;

    adja = &mesh->adja[3*(k-1)+1];
    prevk = k;
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = _MMG5_inxt2[i];
  }
  while ( k && k != start );

  if ( k > 0 ) return(nbpoi);

  /* store the last point of the boundary triangle */
  if ( nbpoi == MMGS_LMAX ) {
    fprintf(stderr,"  ## Warning: unable to compute adjacent vertices of the"
            " vertex %d:\nthe ball of point contain too many elements.\n",ip);
    return(0);
  }
  i1 = _MMG5_inxt2[i1];
  lispoi[nbpoi] = mesh->tria[prevk].v[i1];
  ++nbpoi;

  /* check if boundary hit */
  k = start;
  i = iploc;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i2 = _MMG5_iprv2[i];
    k  = adja[i2] / 3;
    if ( k == 0 )  break;

    if ( nbpoi == MMGS_LMAX ) {
      fprintf(stderr,"  ## Warning: unable to compute adjacent vertices of the"
              " vertex %d:\nthe ball of point contain too many elements.\n",ip);
      return(0);
    }
    i  = adja[i2] % 3;
    lispoi[nbpoi] = mesh->tria[k].v[i];
    ++nbpoi;

    i  = _MMG5_iprv2[i];
  }
  while ( k );

  return nbpoi;
}
