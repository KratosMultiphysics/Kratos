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

/**
 * Pack the mesh \a mesh and its associated metric \a met and return \a val.
 */
#define _MMG2D_RETURN_AND_PACK(mesh,met,val)do                          \
  {                                                                     \
  if ( !MMG2_tassage(mesh,met) ) {                                      \
    mesh->npi = mesh->np;                                               \
    mesh->nti = mesh->nt;                                               \
    mesh->nai = mesh->na;                                               \
    mesh->nei = mesh->ne;                                               \
    met->npi  = met->np;                                                \
    return(MMG5_LOWFAILURE);                                            \
  }                                                                     \
    _LIBMMG5_RETURN(mesh,met,val);                                      \
    }while(0)


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the solution structure.
 *
 * Truncate a scalar metric by hmax and hmin values.
 *
 */
static inline
void _MMG2D_scalarSolTruncature(MMG5_pMesh mesh, MMG5_pSol met) {
  int         k;

  /* vertex size */
  for (k=1; k<=mesh->np; k++) {
    met->m[k] = MG_MIN(mesh->info.hmax,MG_MAX(mesh->info.hmin,met->m[k]));
  }
  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \return 0 if memory problem (uncomplete mesh), 1 otherwise.
 *
 * Pack the mesh and metric and create explicitly all the mesh structures
 * (edges).
 *
 * \warning edges are not packed.
 */
static inline
int MMG2_tassage(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pEdge         ped;
  MMG5_pTria         pt,ptnew;
  MMG5_pPoint        ppt,pptnew;
  _MMG5_Hash         hash;
  int                np,nt,k,nbl,isol,isolnew,i,memWarn,num;
  int                iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;

  /* compact vertices */
  np=0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & MG_NUL ) {
      ppt->tmp = 0;
      continue;
    }
    ppt->tmp = ++np;
  }


  /* compact edges */
  memWarn = 0; num = 0;

  if ( (!mesh->na) && mesh->adja ) {
    //printf("NO EDGES\n");
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if (!pt->v[0]) continue;
      for (i=0 ; i<3 ; i++) {
        if (!(&mesh->adja[3*(k-1)+1])[i]) ++mesh->na;
      }
    }

    if ( mesh->na ) {
      assert(!mesh->edge);
      _MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"final edges",
                    memWarn=1);

      if ( memWarn ) {
        if ( mesh->info.ddebug )
          printf("  -- Attempt to allocate a smallest edge table...\n");
        mesh->namax = mesh->na;
        memWarn = 0;
        _MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"final edges",
                      printf("  ## Warning: uncomplete mesh.\n");
                      memWarn=1);
      }

      if ( memWarn ) {
        mesh->na = 0;
        goto triangles;
      }

      /* We have enough memory to allocate the edge table */
      _MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge);

      for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if (!pt->v[0]) continue;
        for (i=0 ; i<3 ; i++) {
          if ((&mesh->adja[3*(k-1)+1])[i]) continue;
          ++num;
          ped = &mesh->edge[num];
          ped->a = pt->v[MMG2_iare[i][0]];
          ped->b = pt->v[MMG2_iare[i][1]];
          ped->base = 3*k+i;
          ped->ref  = M_MIN(mesh->point[pt->v[MMG2_iare[i][0]]].ref,
                            mesh->point[pt->v[MMG2_iare[i][1]]].ref);
        }
      }
    }
  }

  nbl = 0;

  for (k=1; k<=mesh->na; k++) {
    ped  = &mesh->edge[k];

    if ( !ped->a )
      continue;
    else if ( !ped->b ) {
      ped->a = 0;
      continue;
    }

    ped->a = mesh->point[ped->a].tmp;
    ped->b = mesh->point[ped->b].tmp;

    /* impossible to do that without update triangle....*/
    /* if(k!=nbl) { */
    /*   pednew = &mesh->edge[nbl]; */
    /*   memcpy(pednew,ped,sizeof(MMG5_Edge)); */
    /*   memset(ped,0,sizeof(MMG5_Tria)); */
    /* } */
    /*nbl++;*/
  }

  ped = &mesh->edge[mesh->na];

  while ( !ped->a ) {
    ped = &mesh->edge[--mesh->na];
  }

  /* mesh->na = nbl;*/

  /* compact triangle */
triangles:
  nt  = 0;
  nbl = 1;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] )  {
      continue;
    }
    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    nt++;
    if(k!=nbl) {
      ptnew = &mesh->tria[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tria));
      //and the adjacency
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadrnew = 3*(nbl-1) + 1;
      adjanew = &mesh->adja[iadrnew];
      for(i=0 ; i<3 ; i++) {
        adjanew[i] = adja[i];
        if(!adja[i]) continue;
        iadrv = 3*(adja[i]/3-1) +1;
        adjav = &mesh->adja[iadrv];
        voy = i;
        adjav[adja[i]%3] = 3*nbl + voy;

        adja[i] = 0;
      }
      memset(pt,0,sizeof(MMG5_Tria));
    }
    nbl++;
  }
  mesh->nt = nt;

  /* Travel through the tria and hash the boundary edges in order to recover
   * from which tria comes a boundary edges */
  hash.item = NULL;
  if ( _MMG5_hashNew(mesh,&hash,mesh->nt,3*mesh->nt) ) {

    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if (!pt->v[0]) continue;
      for (i=0 ; i<3 ; i++) {
        if ( !_MMG5_hashEdge(mesh,&hash,pt->v[_MMG5_inxt2[i]],pt->v[_MMG5_iprv2[i]],3*k+i) ) {
          fprintf(stdout,"  ## Warning: unable hash boundary edges.\n");
          break;
        }
      }
    }

  }

  for (k=1; k<=mesh->na; k++) {
    ped  = &mesh->edge[k];
    if(!ped->a || !ped->b) continue;

    ped->base = _MMG5_hashGet(&hash,ped->a,ped->b);
  }

  if ( hash.item )  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));


  /* compact metric */
  if ( sol->m ) {
    nbl = 1;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( ppt->tag & MG_NUL )  continue;
      isol    = k * sol->size;
      isolnew = nbl* sol->size;

      for (i=0; i<sol->size; i++)
        sol->m[isolnew + i] = sol->m[isol + i];
      ++nbl;
    }
  }

  /*compact vertices*/
  np  = 0;
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( ppt->tag & MG_NUL )  continue;

    /* remove the required tag at vertices in -nosurf mode */
    if ( mesh->info.nosurf && (ppt->tag & M_NOSURF) ) {
      ppt->tag &= ~M_REQUIRED;
    }

    if(k!=nbl) {
      pptnew = &mesh->point[nbl];
      memcpy(pptnew,ppt,sizeof(MMG5_Point));
      ppt->tag   &= ~MG_NUL;
      assert(ppt->tmp == nbl);
    }
    np++;
    if(k != nbl) {
      ppt = &mesh->point[k];
      memset(ppt,0,sizeof(MMG5_Point));
      ppt->tag    = MG_NUL;
    }
    nbl++;
  }
  mesh->np = np;
  if ( sol->m ) sol->np  = np;

  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  if(mesh->np < mesh->npmax - 3) {
    mesh->npnil = mesh->np + 1;
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;
  } else {
    mesh->npnil = 0;
  }

  /*to do only if the edges are packed*/
  /* if(mesh->na < mesh->namax - 3) { */
  /*   mesh->nanil = mesh->na + 1; */
  /*   for (k=mesh->nanil; k<mesh->namax-1; k++) */
  /*     mesh->edge[k].b = k+1; */
  /* } else { */
  /*   mesh->nanil = 0; */
  /* } */

  if(mesh->nt < mesh->ntmax - 3) {
    mesh->nenil = mesh->nt + 1;
    for (k=mesh->nenil; k<mesh->ntmax-1; k++)
      mesh->tria[k].v[2] = k+1;
  } else {
    mesh->nenil = 0;
  }

  if ( memWarn ) return 0;

  return(1);
}

/*
  opt[0] = option
  opt[1] = ddebug
  opt[2] = noswap
  opt[3] = noinsert
  opt[4] = nomove
  opt[5] = imprim
  opt[6] = nr

  optdbl[0] = hgrad
  optdbl[1] =ar
*/

int MMG2D_mmg2dlib(MMG5_pMesh mesh,MMG5_pSol sol)
//,void (*titi)(int ,int ,int,int,int)
{
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /*uncomment to callback*/
  //MMG2D_callbackinsert = titi;
  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( !mesh->nt ) {
    fprintf(stdout,"\n  ## ERROR: NO TRIANGLES IN THE MESH. \n");
    fprintf(stdout,"          To generate a mesh from boundaries call the"
            " MMG2D_mmg2dmesh function\n.");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stdout,"  ## Error: level-set discretisation unavailable"
            " (MMG2D_IPARAM_iso):\n"
            "          You must call the MMG2D_mmg2dls function to use this option.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  } else if ( mesh->info.lag >= 0 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG2D_IPARAM_lag):\n"
            "            You must call the MMG2D_mmg2dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG2DLIB: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  sol->np = mesh->np;
  sol->ver  = mesh->ver;

  if ( !sol->m ) {
    /* mem alloc */
    _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax+1))*sizeof(double),
                  "initial solution",return(0));
    _MMG5_SAFE_CALLOC(sol->m,sol->size*(mesh->npmax+1),double);
    sol->np = 0;
  } else   if ( sol->np && (sol->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER : %d != %d\n",sol->np,mesh->np);
    //exit(EXIT_FAILURE);
  }
  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* /\* default values *\/ */
  /* mesh->info.imprim = opt[5]; */
  /* mesh->info.mem    = 0; */
  /* mesh->info.ddebug = opt[1]; */
  /* mesh->info.iso = 0; */
  /* mesh->info.lag = -1; */
  /* mesh->info.hmin = -1; */
  /* mesh->info.hmax = -1; */
  /* mesh->info.hausd = 0.01; */
  /* switch(opt[0]) { */
  /* case 0: */
  /* case 1: */
  /* case 2: */
  /*   break; */
  /* case 6: */
  /*   mesh->info.iso = 1; */
  /*   break; */
  /* case 9: */
  /* case 99: */
  /*   mesh->info.lag = 0; */
  /*   break; */
  /* default: */
  /*   fprintf(stdout,"option not recognized %d\n",opt[0]); */
  /*   exit(EXIT_FAILURE); */
  /* } */
  /* mesh->info.noswap = opt[2]; */
  /* mesh->info.nomove = opt[4]; */
  /* mesh->info.noinsert = opt[3]; */
  /* mesh->info.hgrad  = optdbl[0]; */
  /* if(opt[6]) */
  /*   mesh->info.dhd  = -1; */
  /* else */
  /*   mesh->info.dhd  = 180.-optdbl[1]; */
  /* /\*this options are not used inside library version*\/ */
  /* //qdegrad[0] = 10./ALPHA; */
  /* //qdegrad[1] = 1.3; */
  /* mesh->info.renum = 0; */

  /* sol->type = 1; */

  MMG2D_setfunc(mesh,sol);
  _MMG2D_Set_commonFunc();

  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          MG_STR,MG_VER,MG_REL,sol->size == 1 ? "ISO" : "ANISO",MG_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* analysis */
  chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");

  if ( !MMG2_scaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( !sol->np ) {
    if ( !MMG2D_doSol(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_scalarSolTruncature(mesh,sol);
  }

  if ( (mesh)->adja )
     _MMG5_DEL_MEM((mesh),(mesh)->adja,(3*(mesh)->ntmax+5)*sizeof(int));

  if ( !MMG2_hashel(mesh) )
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /*geom : corner detection*/
  if ( mesh->info.dhd>0 )
    if( !MMG2_evalgeom(mesh) ) _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Update the metric definition by taking into accounts the
     curvature of the external and internal curves present in the mesh */
  if ( sol->size==1 && !_MMG2D_defBdrySiz(mesh,sol) ) return(0);

  /*mesh gradation*/
  if( mesh->info.hgrad > 0 ) {
    if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",mesh->info.hgrad);
    MMG2_lissmet(mesh,sol);
  }
  MMG2_outqua(mesh,sol);

  if ( abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* remeshing */
  chrono(ON,&ctim[3]);

  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : MESH ADAPTATION\n");
  if ( (!mesh->info.noinsert) && !MMG2_mmg2d1(mesh,sol) ) {
    if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  /* optimisation */
  chrono(ON,&ctim[4]);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 3 : MESH OPTIMISATION\n");
  if ( !MMG2_mmg2d0(mesh,sol) ) {
    if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[4]));
  printim(ctim[4].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG2D: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  MMG2_outqua(mesh,sol);

  if ( abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if (!MMG2_tassage(mesh,sol) ) _LIBMMG5_RETURN(mesh,sol,MMG5_LOWFAILURE);

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG2DLIB: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,sol,MMG5_SUCCESS);

}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail (lack of memory), 1 otherwise.
 *
 * Clean the mesh structure when we just call the MMG2D_Free_Triangles and
 * MMG2D_Free_Edges functions between 2 call of the MMG2D_mmg2dmesh function:
 *   - Allocate the tria and edge structures if needed;
 *   - Reset the tags at vertices.
 *
 */
static inline
int _MMG2D_restart(MMG5_pMesh mesh){
  int k;

 /** If needed, reallocate the missing structures */
  if ( !mesh->tria ) {
    /* If we call the library more than one time and if we free the triangles
     * using the MMG2D_Free_triangles function we need to reallocate it */
    _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),
                  "initial triangles",return(0));
    _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria);
    mesh->nenil = mesh->nt + 1;
    for ( k=mesh->nenil; k<mesh->ntmax-1; k++) {
      mesh->tria[k].v[2] = k+1;
    }
  }
  if ( !mesh->edge ) {
    /* If we call the library more than one time and if we free the triangles
     * using the MMG2D_Free_triangles function we need to reallocate it */
    _MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),
                  "initial edges",return(0));
    _MMG5_SAFE_CALLOC(mesh->edge,mesh->namax+1,MMG5_Edge);
    mesh->nanil = mesh->na + 1;
    for ( k=mesh->nanil; k<mesh->namax-1; k++) {
      mesh->edge[k].b = k+1;
    }
  }

  for ( k=1; k<=mesh->np;  ++k ) {
    mesh->point[k].tag = 0;
  }
  return 1;
}

int MMG2D_mmg2dmesh(MMG5_pMesh mesh,MMG5_pSol sol) {
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /*uncomment for callback*/
  //MMG2D_callbackinsert = titi;

  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->nt ) {
    fprintf(stdout,"  ## Error: your mesh contains already triangles.\n"
            " The mesh generation option is unavailable.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso ) {
    fprintf(stdout,"  ## Error: level-set discretisation unavailable"
            " (MMG2D_IPARAM_iso):\n"
            "          You must call the MMG2D_mmg2dls function to use this option.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  } else if ( mesh->info.lag >= 0 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG2D_IPARAM_lag):\n"
            "            You must call the MMG2D_mmg2dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG2DMESH: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  sol->ver  = mesh->ver;
  if ( !sol->m ) {
    /* mem alloc */
    _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax+1))*sizeof(double),
                  "initial solution",return(0));
    _MMG5_SAFE_CALLOC(sol->m,sol->size*(mesh->npmax+1),double);
    sol->np = 0;
  } else   if ( sol->np && (sol->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER : %d != %d\n",sol->np,mesh->np);
    //exit(EXIT_FAILURE);
  }
  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  MMG2D_setfunc(mesh,sol);
  _MMG2D_Set_commonFunc();

  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          MG_STR,MG_VER,MG_REL,sol->size == 1 ? "ISO" : "ANISO",MG_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* analysis */
  chrono(ON,&ctim[2]);

  if ( !_MMG2D_restart(mesh) ) {
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  };

  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( !sol->np ) {
    if ( !MMG2D_doSol(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_scalarSolTruncature(mesh,sol);
  }

  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  /*geom : corner detection*/
  if ( mesh->info.dhd>0 )
    if( !MMG2_evalgeom(mesh) ) _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* specific meshing */
  chrono(ON,&ctim[3]);

  /* memory alloc */
  _MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int);

  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : MESH GENERATION\n");

  if ( !MMG2_mmg2d2(mesh,sol) )  {
    if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);

  /* remeshing */
  chrono(ON,&ctim[4]);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 3 : MESH IMPROVEMENT\n");
  }

  /* analysis */
  /*geom : corner detection*/
  if ( mesh->info.dhd>0 )
    if( !MMG2_evalgeom(mesh) ) _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Update the metric definition by taking into accounts the
     curvature of the external and internal curves present in the mesh */
  if ( sol->size==1 && !_MMG2D_defBdrySiz(mesh,sol) ) return(0);

  /*mesh gradation*/
  if( mesh->nt && mesh->info.hgrad > 0 ) {
    if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",mesh->info.hgrad);
    MMG2_lissmet(mesh,sol);
  }
  if ( mesh->nt )  MMG2_outqua(mesh,sol);

  if ( mesh->nt && abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }
  if ( (!mesh->info.noinsert) && !MMG2_mmg2d1(mesh,sol) ) {
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  /* optimisation */
  chrono(ON,&ctim[5]);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 4 : MESH OPTIMISATION\n");
  if ( !MMG2_mmg2d0(mesh,sol) ) {
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[5]));
  printim(ctim[5].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 4 COMPLETED.     %s\n",stim);
  }

  chrono(OFF,&(ctim[4]));
  printim(ctim[4].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMGS: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  MMG2_outqua(mesh,sol);

  if ( abs(mesh->info.imprim) > 1 )  {
    MMG2_prilen(mesh,sol);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if (!MMG2_tassage(mesh,sol) ) _LIBMMG5_RETURN(mesh,sol,MMG5_LOWFAILURE);

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG2DMESH: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,sol,MMG5_SUCCESS);

}

int MMG2D_mmg2dls(MMG5_pMesh mesh,MMG5_pSol sol)
{
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  /*uncomment for callback*/
  //MMG2D_callbackinsert = titi;

  /* interrupts */
  signal(SIGABRT,_MMG2_excfun);
  signal(SIGFPE,_MMG2_excfun);
  signal(SIGILL,_MMG2_excfun);
  signal(SIGSEGV,_MMG2_excfun);
  signal(SIGTERM,_MMG2_excfun);
  signal(SIGINT,_MMG2_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag >= 0 ) {
    fprintf(stdout,"  ## Error: lagrangian mode unavailable (MMG2D_IPARAM_lag):\n"
            "            You must call the MMG2D_mmg2dmov function to move a rigidbody.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMG2DLS: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  sol->ver  = mesh->ver;

  if ( !mesh->nt ) {
    fprintf(stdout,"\n  ## ERROR: NO TRIANGLES IN THE MESH \n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  }
  else if ( !sol->m ) {
    fprintf(stdout,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  } else   if ( sol->size != 1 ) {
    fprintf(stdout,"  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  } else if ( sol->np && (sol->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));
    sol->np = 0;
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);


  MMG2D_setfunc(mesh,sol);
  _MMG2D_Set_commonFunc();

  fprintf(stdout,"\n  %s\n   MODULE MMG2D-IMB/LJLL : %s (%s) %s\n  %s\n",
          MG_STR,MG_VER,MG_REL,sol->size == 1 ? "ISO" : "ANISO",MG_STR);
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM NUMBER OF POINTS    (NPMAX) : %8d\n",mesh->npmax);
    fprintf(stdout,"  MAXIMUM NUMBER OF TRIANGLES (NTMAX) : %8d\n",mesh->ntmax);
  }

  /* analysis */
  chrono(ON,&ctim[2]);
  if ( mesh->info.imprim )   fprintf(stdout,"\n  -- PHASE 1 : DATA ANALYSIS\n");
  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** SETTING ADJACENCIES\n");
  if ( !MMG2_scaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  if ( mesh->nt && !MMG2_hashel(mesh) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  if ( mesh->info.ddebug && !_MMG5_chkmsh(mesh,1,0) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  /*geom : corner detection*/
  if ( mesh->info.dhd>0 )
    if( !MMG2_evalgeom(mesh) ) _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  MMG2_outqua(mesh,sol);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  chrono(ON,&ctim[3]);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 2 : LEVEL-SET DISCRETIZATION\n");
  }

  /* specific meshing */
  if (! MMG2_mmg2d6(mesh,sol) ) _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);


  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);


  /* mesh improvement */
  chrono(ON,&ctim[4]);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  -- PHASE 3 : MESH IMPROVEMENT\n");
  }

  /* analysis */
  /* metric creation */
  _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax+1))*sizeof(double),
                "initial solution",return(0));
  _MMG5_SAFE_CALLOC(sol->m,sol->size*(mesh->npmax+1),double);
  sol->np = 0;
  if ( !MMG2D_doSol(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
  _MMG2D_scalarSolTruncature(mesh,sol);

  /*geom : corner detection*/
  if ( mesh->info.dhd>0 )
    if( !MMG2_evalgeom(mesh) ) _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  /* Update the metric definition by taking into accounts the
     curvature of the external and internal curves present in the mesh */
  if ( sol->size==1 && !_MMG2D_defBdrySiz(mesh,sol) ) return(0);

  /*mesh gradation*/
  if( mesh->info.hgrad > 0 ) {
    if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",
                                       mesh->info.hgrad);
    MMG2_lissmet(mesh,sol);
  }

  if ( (!mesh->info.noinsert) && !MMG2_mmg2d1(mesh,sol) ) {
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  /* optimisation */
  chrono(ON,&ctim[5]);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 4 : MESH OPTIMISATION\n");
  if ( !MMG2_mmg2d0(mesh,sol) ) {
    _MMG2D_RETURN_AND_PACK(mesh,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[5]));
  printim(ctim[5].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 4 COMPLETED.     %s\n",stim);
  }

  chrono(OFF,&(ctim[4]));
  printim(ctim[4].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 5 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMGS: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  if ( !MMG2_unscaleMesh(mesh,sol) )  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);

  MMG2_outqua(mesh,sol);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if (!MMG2_tassage(mesh,sol) ) _LIBMMG5_RETURN(mesh,sol,MMG5_LOWFAILURE);

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMG2DLIB: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,sol,MMG5_SUCCESS);

}

int MMG2D_mmg2dmov(MMG5_pMesh mesh,MMG5_pSol sol)
// MMG5_pSol,met)
{
  //mytime    ctim[TIMEMAX];
  //char      stim[32];

  fprintf(stdout,"  -- MMG2D, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  fprintf(stdout,"\n  ## ERROR: OPTIONALITY NOT YET IMPLEMENTED. \n");

  _LIBMMG5_RETURN(mesh,sol,MMG5_STRONGFAILURE);
}
