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
 * \file mmgs/libmmgs.c
 * \brief API functions for MMGS library.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \todo documentation doxygen
 *
 * Private API functions for MMGS library: incompatible functions
 * with the main binary.
 *
 */

#include "mmgs.h"

/**
 * Pack the mesh \a mesh and its associated metric \a met and return \a val.
 */
#define _MMGS_RETURN_AND_PACK(mesh,met,val)do                           \
  {                                                                     \
    if ( !_MMGS_packMesh(mesh,met) )  {                                 \
      mesh->npi = mesh->np;                                             \
      mesh->nti = mesh->nt;                                             \
      mesh->nai = mesh->na;                                             \
      mesh->nei = mesh->ne;                                             \
      met->npi  = met->np;                                              \
      return(MMG5_LOWFAILURE);                                          \
    }                                                                   \
    _LIBMMG5_RETURN(mesh,met,val);                                      \
  }while(0)

/** Free adja, xtetra and xpoint tables */
static inline
void _MMGS_Free_topoTables(MMG5_pMesh mesh) {
  int k;

  mesh->xp = 0;
  if ( mesh->adja )
    _MMG5_DEL_MEM(mesh,mesh->adja,(3*mesh->ntmax+5)*sizeof(int));

  if ( mesh->xpoint )
    _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));

  for(k=1; k <=mesh->np; k++) {
    mesh->point[k].xp = 0;
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param met pointer toward the solution (metric) structure.
 *
 * Pack the sparse mesh and create edges before getting
 * out of library
 *
 */
static inline
int _MMGS_packMesh(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pTria    pt,ptnew;
  MMG5_pPoint   ppt,pptnew;
  int           np,nc,nr, k,nt,nbl,imet,imetnew,i,na,jel;
  int           iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;
  char          i1,i2;

  /* compact vertices */
  np = nc = nr = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++np;
    if ( ppt->tag & MG_CRN )  nc++;
    ppt->ref = abs(ppt->ref);
  }

  /* compact triangles */
  nt  = 0;
  na  = 0;
  nbl = 1;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    nt++;
    if ( k!=nbl ) {
      ptnew = &mesh->tria[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tria));

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
      }
    }
    nbl++;

    /* Count the edges */
    for(i=0 ; i<3 ; i++) {
      if ( !MG_EDG(pt->tag[i]) ) continue;

      adja = &mesh->adja[3*(k-1)+1];
      jel  = adja[i] / 3;

      if ( jel && jel <= k ) continue;
      ++na;
    }

  }
  mesh->nt = nt;

  /* compact metric */
  nbl = 1;
  if ( met && met->m ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      imet    = k   * met->size;
      imetnew = nbl * met->size;

      for (i=0; i<met->size; i++)
        met->m[imetnew + i] = met->m[imet + i];
      ++nbl;
    }
  }

  /*compact vertices*/
  np  = 0;
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    np++;
    if ( k!=nbl ) {
      pptnew = &mesh->point[nbl];
      memmove(pptnew,ppt,sizeof(MMG5_Point));
      memset(ppt,0,sizeof(MMG5_Point));
      ppt->tag    = MG_NUL;
    }
    nbl++;
  }
  mesh->np = np;
  if ( met && met->m )
    met->np  = np;

  /* memory alloc */
  mesh->na = 0;
  if ( mesh->edge ) {
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));
    _MMG5_SAFE_FREE(mesh->edge);
  }

  if ( na ) {
    _MMG5_ADD_MEM(mesh,(na+1)*sizeof(MMG5_Edge),"final edges",
                  na = 0;
                  printf("  ## Warning: uncomplete mesh\n")
      );
  }

  if ( na ) {
    _MMG5_SAFE_CALLOC(mesh->edge,na+1,MMG5_Edge);
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( MG_EOK(pt) ) {
        for (i=0; i<3; i++) {
          if ( !MG_EDG(pt->tag[i]) )  continue;

          adja = &mesh->adja[3*(k-1)+1];
          jel  = adja[i] / 3;
          if ( jel && jel <= k ) continue;
          i1 = _MMG5_inxt2[i];
          i2 = _MMG5_inxt2[i1];
          mesh->na++;
          mesh->edge[mesh->na].a    = mesh->point[pt->v[i1]].tmp;
          mesh->edge[mesh->na].b    = mesh->point[pt->v[i2]].tmp;
          mesh->edge[mesh->na].ref  = pt->edg[i];
          mesh->edge[mesh->na].tag |= pt->tag[i];
          if ( pt->tag[i] & MG_GEO )  nr++;
        }
      }
    }
  }

  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  mesh->npnil = mesh->np + 1;
  for(k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  mesh->nenil = mesh->nt + 1;
  for(k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  /* to could save the mesh, the adjacency have to be correct */
  if ( mesh->info.ddebug && (!_MMG5_chkmsh(mesh,1,1) ) ) {
    fprintf(stderr,"  ##  Problem. Invalid mesh.\n");
    return(0);
  }

  if ( mesh->info.imprim ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d   RIDGES  %8d\n",mesh->na,nr);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);
  }
  return(1);
}

int MMGS_mmgsls(MMG5_pMesh mesh,MMG5_pSol met)
{
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMGS, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  _MMGS_Set_commonFunc();

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  _MMGS_Free_topoTables(mesh);

  /* trap exceptions */
  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMGSLS: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stderr,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
    met->np = 0;
  }
  else if ( met->size!=1 ) {
    fprintf(stderr,"  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  MMGS_setfunc(mesh,met);

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMGSLS: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMGS_inqua(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }
  }

  /* specific meshing */
  if ( !met->np ) {
    fprintf(stderr,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }
  if ( !_MMGS_mmgs2(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* mesh analysis */
  if ( !_MMGS_analys(mesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
      fprintf(stdout,"\n  -- PHASE 2 : ISOTROPIC MESHING\n");
  }

  if ( !_MMG5_mmgs1(mesh,met) ) {
    if ( (!mesh->adja) && !_MMGS_hashTria(mesh) ) {
      fprintf(stderr,"  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMGSLS: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* save file */
  if (!_MMGS_outqua(mesh,met) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  if ( !_MMGS_packMesh(mesh,met) )     _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMGSLS: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}

int MMGS_mmgslib(MMG5_pMesh mesh,MMG5_pSol met)
{
  mytime    ctim[TIMEMAX];
  char      stim[32];

  fprintf(stdout,"  -- MMGS, Release %s (%s) \n",MG_VER,MG_REL);
  fprintf(stdout,"     %s\n",MG_CPY);
  fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

  _MMGS_Set_commonFunc();

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  _MMGS_Free_topoTables(mesh);

  /* trap exceptions */
  signal(SIGABRT,_MMG5_excfun);
  signal(SIGFPE,_MMG5_excfun);
  signal(SIGILL,_MMG5_excfun);
  signal(SIGSEGV,_MMG5_excfun);
  signal(SIGTERM,_MMG5_excfun);
  signal(SIGINT,_MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  if ( mesh->info.iso ) {
    fprintf(stderr,"  ## Error: level-set discretisation unavailable"
            " (MMGS_IPARAM_iso):\n"
            "          You must call the MMGS_mmgsls function to use this option.\n");
    _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  _MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim ) fprintf(stdout,"\n  -- MMGS: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

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
  MMGS_setfunc(mesh,met);

  if ( mesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE MMGS: IMB-LJLL : %s (%s)\n  %s\n",MG_STR,MG_VER,MG_REL,MG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  if ( !_MMG5_scaleMesh(mesh,met) ) _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  /* mesh analysis */
  if ( !_MMGS_analys(mesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  if ( abs(mesh->info.imprim) > 0 ) {
    if ( !_MMGS_inqua(mesh,met) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
      _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
    }
  }

  if ( mesh->info.imprim > 1 && !mesh->info.iso && met->m )
    _MMGS_prilen(mesh,met,0);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim ) {
      fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");
  }

  if ( !_MMG5_mmgs1(mesh,met) ) {
    if ( (!mesh->adja) && !_MMGS_hashTria(mesh) ) {
      fprintf(stderr,"  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    }
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMGS: IMB-LJLL \n  %s\n",MG_STR,MG_STR);
  }

  /* save file */
  if (!_MMGS_outqua(mesh,met) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
    _MMGS_RETURN_AND_PACK(mesh,met,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 1 )
    _MMGS_prilen(mesh,met,1);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if ( !_MMG5_unscaleMesh(mesh,met) )  _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);

  if ( !_MMGS_packMesh(mesh,met) )     _LIBMMG5_RETURN(mesh,met,MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim )
    fprintf(stdout,"\n   MMGSLIB: ELAPSED TIME  %s\n",stim);
  _LIBMMG5_RETURN(mesh,met,MMG5_SUCCESS);
}
