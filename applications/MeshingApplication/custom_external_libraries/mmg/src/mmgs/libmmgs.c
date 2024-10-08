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

#include "libmmgs.h"
#include "libmmgs_private.h"
#include "mmgsexterns_private.h"
#include "mmgexterns_private.h"

/**
 * Pack the mesh \a mesh and its associated metric \a met and return \a val.
 */
#define MMGS_RETURN_AND_PACK(mesh,met,sol,val)do  \
  {                                               \
    if ( !MMGS_packMesh(mesh,met,sol) )  {        \
      mesh->npi = mesh->np;                       \
      mesh->nti = mesh->nt;                       \
      mesh->nai = mesh->na;                       \
      mesh->nei = mesh->ne;                       \
      met->npi  = met->np;                        \
      if ( met ) { met->npi  = met->np; }         \
      if ( sol ) { sol->npi  = sol->np; }         \
      return MMG5_LOWFAILURE;                     \
    }                                             \
    _LIBMMG5_RETURN(mesh,met,sol,val);            \
  }while(0)

/** Free adja, xtetra and xpoint tables */
static inline
void MMGS_Free_topoTables(MMG5_pMesh mesh) {
  MMG5_int k;

  mesh->xp = 0;
  if ( mesh->adja )
    MMG5_DEL_MEM(mesh,mesh->adja);

  if ( mesh->xpoint )
    MMG5_DEL_MEM(mesh,mesh->xpoint);

  for(k=1; k <=mesh->np; k++) {
    mesh->point[k].xp = 0;
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param sol pointer toward a solution structure.
 * \param met pointer toward the solution (metric) structure.
 *
 * Pack the sparse mesh and create edges before getting
 * out of library
 *
 */
static inline
int MMGS_packMesh(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  MMG5_pTria    pt,ptnew;
  MMG5_pPoint   ppt,pptnew;
  MMG5_int      iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;
  MMG5_int      k,nt,np,na,jel,nc,nr,nbl;
  int           imet,imetnew,i;
  int8_t        i1,i2;

  /* Remove non wanted subdomains if needed */
  MMGS_keep_only1Subdomain ( mesh, mesh->info.nsd );

  /* compact vertices */
  np = nc = nr = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++np;
    if ( ppt->tag & MG_CRN )  nc++;
    ppt->ref = MMG5_abs(ppt->ref);
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

      assert ( jel != k );
      if ( jel ) {
        /* If we have an adjacent, treat the edge either from the tria with
         * higher ref or, if both tria have the same references, from the tria
         * with higher index */
        if ( mesh->tria[jel].ref < mesh->tria[k].ref ) {
          continue;
        }
        else if ( (mesh->tria[jel].ref == mesh->tria[k].ref) && (k < jel)  ) {
          continue;
        }
      }
      ++na;
    }

  }
  mesh->nt = nt;

  /* compact solutions (metric, level-set, displacement...) */
  if ( met && met->m ) {
    nbl = 1;
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

  if ( sol && sol->m ) {
    nbl = 1;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      imet    = k   * sol->size;
      imetnew = nbl * sol->size;

      for (i=0; i<sol->size; i++)
        sol->m[imetnew + i] = sol->m[imet + i];
      ++nbl;
    }
  }

  /* compact vertices */
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

  if ( sol && sol->m )
    sol->np  = np;

  /* memory alloc */
  mesh->na = 0;
  if ( mesh->edge ) {
    MMG5_DEL_MEM(mesh,mesh->edge);
    MMG5_SAFE_FREE(mesh->edge);
  }

  if ( na ) {
    MMG5_ADD_MEM(mesh,(na+1)*sizeof(MMG5_Edge),"final edges",
                 na = 0;
                 printf("  ## Warning: uncomplete mesh\n")
      );
  }

  if ( na ) {
    MMG5_SAFE_CALLOC(mesh->edge,na+1,MMG5_Edge,return 0);
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( MG_EOK(pt) ) {
        for (i=0; i<3; i++) {
          if ( !MG_EDG(pt->tag[i]) )  continue;

          adja = &mesh->adja[3*(k-1)+1];
          jel  = adja[i] / 3;
          if ( jel ) {
            /* If we have an adjacent, treat the edge either from the tria with
             * higher ref or, if both tria have the same references, from the tria
             * with higher index */
            if ( mesh->tria[jel].ref < mesh->tria[k].ref ) {
              continue;
            }
            else if ( (mesh->tria[jel].ref == mesh->tria[k].ref) && (k < jel)  ) {
              continue;
            }
          }

          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
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
  if ( mesh->info.ddebug && (!MMG5_chkmsh(mesh,1,1) ) ) {
    fprintf(stderr,"\n  ##  Warning: %s: invalid mesh.\n",__func__);
    return 0;
  }

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8" MMG5_PRId "   CORNERS %8" MMG5_PRId "\n",mesh->np,nc);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8" MMG5_PRId "\n",mesh->nt);

    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8" MMG5_PRId "   RIDGES  %8" MMG5_PRId "\n",mesh->na,nr);
  }
  return 1;
}

int MMGS_mmgsls(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol umet)
{
  MMG5_pSol met=NULL;
  mytime    ctim[TIMEMAX];
  int8_t    mettofree = 0;
  char      stim[32];

  /** In debug mode, check that all structures are allocated */
  assert ( mesh );
  assert ( sol );
  assert ( mesh->point );
  assert ( mesh->tria );

  MMG5_version(mesh,"S");

  if ( (!mesh->info.iso) && (!mesh->info.isosurf) ) {
    fprintf(stdout,"\n  ## WARNING: ISO MODE NOT PROVIDED: ENABLING ISOVALUE DISCRETIZATION MODE (-ls) \n");
    mesh->info.iso = 1;
  }

  if ( !umet ) {
    /* User doesn't provide the metric (library mode only), allocate our own one */
    MMG5_SAFE_CALLOC(met,1,MMG5_Sol,_LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE));
    mettofree = 1;
  }
  else {
    /* Using appli we always pass here */
    met = umet;
  }

  MMGS_Set_commonFunc();

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  MMGS_Free_topoTables(mesh);

  /* trap exceptions */
  signal(SIGABRT,MMG5_excfun);
  signal(SIGFPE,MMG5_excfun);
  signal(SIGILL,MMG5_excfun);
  signal(SIGSEGV,MMG5_excfun);
  signal(SIGTERM,MMG5_excfun);
  signal(SIGINT,MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

#ifdef USE_SCOTCH
  MMG5_warnScotch(mesh);
#endif

  /* Check options */
  /* specific meshing */
  if ( met && met->np ) {
    if ( mesh->info.optim ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }

    if ( mesh->info.hsiz>0. ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  if ( mesh->info.optim &&  mesh->info.hsiz>0. ) {
    printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ AND OPTIM OPTIONS CAN NOT BE USED"
           " TOGETHER.\n");
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.optim ) {
    printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION IS NOT AVAILABLE IN"
           " SURFACIC LEVEL_SET DISCRETIZATION MODE.\n");
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- MMGSLS: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));

  if ( sol->np && (sol->np != mesh->np) ) {
    fprintf(stderr,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    if ( mettofree ) {  MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,sol,met,MMG5_STRONGFAILURE);
  }
  else if ( sol->size!=1 ) {
    fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
    if ( mettofree ) {  MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  if ( met && met->np &&  (met->np != mesh->np) ) {
    fprintf(stdout,"\n  ## WARNING: WRONG METRIC NUMBER. IGNORED\n");
    if ( mettofree ) {  MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,sol,met,MMG5_STRONGFAILURE);
  }

  /* clean old isosurface */
  if ( !MMGS_Clean_isoSurf(mesh) ) {
    fprintf(stderr,"\n  ## Unable to clean old isosurface.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 1 : ISOSURFACE DISCRETIZATION\n");
  }

  /* scaling mesh */
  if ( !MMG5_scaleMesh(mesh,met,sol) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  MMGS_setfunc(mesh,met);

  if ( mesh->info.imprim > 0 || mesh->info.imprim < -1 ) {
    if ( !MMGS_inqua(mesh,met) ) {
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      if ( !MMG5_unscaleMesh(mesh,met,sol) ) {
        _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      }
      MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  if ( !sol->np ) {
    fprintf(stderr,"\n  ## ERROR: A VALID SOLUTION FILE IS NEEDED \n");
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  /* Specific meshing: compute optim option here because after isovalue
   * discretization mesh elements have too bad qualities */
  if ( mesh->info.optim ) {
    /* Mean metric computation */
    if ( !MMGS_doSol(mesh,met) ) {
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      if ( !MMG5_unscaleMesh(mesh,met,sol) ) {
        _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE); }
      MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  /* Discretization of the mesh->info.ls isovalue of sol in the mesh */
  if ( !MMGS_mmgs2(mesh,sol,met) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 2 : ANALYSIS\n");
  }

  /* Specific meshing: compute hsiz on mesh after isovalue discretization */
  if ( mesh->info.hsiz > 0. ) {
    if ( !MMGS_Set_constantSize(mesh,met) ) {
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      if ( !MMG5_unscaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  /* mesh analysis */
  if ( !MMGS_analys(mesh) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( !MMG5_unscaleMesh(mesh,met,sol) )
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[4]));
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 3 : MESH IMPROVEMENT\n");
  }

  if ( !MMG5_mmgs1(mesh,met,NULL) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( (!mesh->adja) && !MMGS_hashTria(mesh) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
    if ( !MMG5_unscaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[4]));
  printim(ctim[4].gdif,stim);
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
  }

  /* save file */
  if (!MMGS_outqua(mesh,met) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( !MMG5_unscaleMesh(mesh,met,sol) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim > 0 )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if ( !MMG5_unscaleMesh(mesh,met,sol) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  if ( !MMGS_packMesh(mesh,met,sol) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n   MMGSLS: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMGS\n  %s\n\n",MG_STR,MG_STR);
  }

  if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
  _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);
}

int MMGS_mmgslib(MMG5_pMesh mesh,MMG5_pSol met)
{
  MMG5_pSol sol=NULL;
  mytime    ctim[TIMEMAX];
  char      stim[32];

  /** In debug mode, check that all structures are allocated */
  assert ( mesh );
  assert ( met );
  assert ( mesh->point );
  assert ( mesh->tria );

  MMG5_version(mesh,"S");

  MMGS_Set_commonFunc();

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  MMGS_Free_topoTables(mesh);

  /* trap exceptions */
  signal(SIGABRT,MMG5_excfun);
  signal(SIGFPE,MMG5_excfun);
  signal(SIGILL,MMG5_excfun);
  signal(SIGSEGV,MMG5_excfun);
  signal(SIGTERM,MMG5_excfun);
  signal(SIGINT,MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  if ( mesh->info.iso || mesh->info.isosurf ) {
    fprintf(stderr,"\n  ## ERROR: LEVEL-SET DISCRETISATION UNAVAILABLE"
            " (MMGS_IPARAM_iso or MMGS_IPARAM_isosurf):\n"
            "          YOU MUST CALL THE MMGS_MMGSLS FUNCTION TO USE THIS OPTION.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- MMGS: INPUT DATA\n");

  /* Check input */
  chrono(ON,&(ctim[1]));

  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,met->m);
    met->np = 0;
  }
  else if ( met->size!=1 && met->size!=6 ) {
    fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  /* specific meshing */

  if ( met->np ) {
    if ( mesh->info.optim ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }

    if ( mesh->info.hsiz>0. ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }
  if ( mesh->info.optim &&  mesh->info.hsiz>0. ) {
    printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ AND OPTIM OPTIONS CAN NOT BE USED"
           " TOGETHER.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* scaling mesh */
  if ( !MMG5_scaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  MMGS_setfunc(mesh,met);

  /* specific meshing */
  if ( mesh->info.hsiz > 0. ) {
    if ( !MMGS_Set_constantSize(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  /* mesh analysis */
  if ( !MMGS_analys(mesh) ) {
    if ( !MMG5_unscaleMesh(mesh,met,NULL) )
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  /* specific meshing: optim mode needs normal at vertices */
  if ( mesh->info.optim ) {
    if ( !MMGS_doSol(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  if ( mesh->info.imprim > 0 ||  mesh->info.imprim < -1 ) {
    if ( !MMGS_inqua(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,NULL) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  if ( mesh->info.imprim > 1 && met->m ) {
    MMGS_prilen(mesh,met,0);
  }

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");
  }

  if ( !MMG5_mmgs1(mesh,met,NULL) ) {
    if ( (!mesh->adja) && !MMGS_hashTria(mesh) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
    if ( !MMG5_unscaleMesh(mesh,met,NULL) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }

  /* save file */
  if (!MMGS_outqua(mesh,met) ) {
    if ( !MMG5_unscaleMesh(mesh,met,NULL) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMGS_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 4 && met->m )
    MMGS_prilen(mesh,met,1);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim > 0 )  fprintf(stdout,"\n  -- MESH PACKED UP\n");

  if ( !MMG5_unscaleMesh(mesh,met,NULL) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  if ( !MMGS_packMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n   MMGSLIB: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMGS\n  %s\n\n",MG_STR,MG_STR);
  }

  _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);


}

/**
* Set common pointer functions between mmgs and mmg3d to the matching mmgs
* functions.
*/
void MMGS_Set_commonFunc(void) {
    MMG5_bezierCP = MMG5_mmgsBezierCP;
    MMG5_chkmsh = MMG5_mmgsChkmsh;
    MMG5_indPt = MMGS_indPt;
    MMG5_indElt = MMGS_indElt;
    MMG5_grad2met_ani = MMG5_grad2metSurf;
    MMG5_grad2metreq_ani = MMG5_grad2metSurfreq;
    MMG5_solTruncature_ani = MMG5_3dSolTruncature_ani;

#ifdef USE_SCOTCH
    MMG5_renumbering = MMG5_mmgsRenumbering;
#endif
}
