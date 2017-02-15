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
 * \file mmg2d/zaldy_2d.c
 * \brief Memory management.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */
#include "mmg2d.h"


/* get new point address */
int _MMG2D_newPt(MMG5_pMesh mesh,double c[2],int16_t tag) {
  MMG5_pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return(0);

  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,2*sizeof(double));
  ppt->tag   &= ~MG_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;
  ppt->xp     = 0;
  //ppt->fla   = mesh->flag;

  return(curpt);
}


void _MMG2D_delPt(MMG5_pMesh mesh,int ip) {
  MMG5_pPoint   ppt;
  MMG5_pxPoint  pxp;

  ppt = &mesh->point[ip];
  if ( ppt->xp ) {
    pxp = &mesh->xpoint[ppt->xp];
    memset(pxp,0,sizeof(MMG5_xPoint));
  }

  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;

  mesh->npnil = ip;
  if ( ip == mesh->np )  mesh->np--;
}

/* get new elt address */
int _MMG5_newEdge(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nanil ) {
    return(0);
  }
  curiel = mesh->nanil;
  if ( mesh->nanil > mesh->na )  mesh->na = mesh->nanil;
  mesh->nanil = mesh->edge[curiel].b;
  mesh->edge[curiel].b = 0;

  return(curiel);
}


void _MMG5_delEdge(MMG5_pMesh mesh,int iel) {
  MMG5_pEdge    pt;

  pt = &mesh->edge[iel];
  if ( !pt->a ) {
    fprintf(stdout,"  ## INVALID EDGE.\n");
    return;
  }
  memset(pt,0,sizeof(MMG5_Edge));
  pt->b = mesh->nanil;
  mesh->nanil = iel;
  if ( iel == mesh->na )  mesh->na--;
}

/* get new elt address */
int _MMG2D_newElt(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nenil ) {
    return(0);
  }
  curiel = mesh->nenil;
  if ( mesh->nenil > mesh->nt )  mesh->nt = mesh->nenil;
  mesh->nenil = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0;
  mesh->tria[curiel].ref = 0;
  mesh->tria[curiel].base = 0;
  mesh->tria[curiel].edg[0] = 0;
  mesh->tria[curiel].edg[1] = 0;
  mesh->tria[curiel].edg[2] = 0;


  return(curiel);
}


void _MMG2D_delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTria    pt;
  int      iadr;

  pt = &mesh->tria[iel];
  if ( !M_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT.\n");
    return;
  }
  memset(pt,0,sizeof(MMG5_Tria));
  pt->v[2] = mesh->nenil;
  pt->qual = 0.0;
  iadr = (iel-1)*3 + 1;
  if ( mesh->adja )
    memset(&mesh->adja[iadr],0,3*sizeof(int));

  mesh->nenil = iel;
  if ( iel == mesh->nt )  mesh->nt--;
}


/* check if n elets available */
int _MMG5_getnElt(MMG5_pMesh mesh,int n) {
  int     curiel;

  if ( !mesh->nenil )  return(0);
  curiel = mesh->nenil;
  do {
    curiel = mesh->tria[curiel].v[2];
  }
  while (--n);

  return(n == 0);
}

/** memory repartition for the -m option */
void _MMG2D_memOption(MMG5_pMesh mesh) {
  long long  million = 1048576L;
  int        ctri,npask,bytes,memtmp;

  mesh->memMax = _MMG5_memSize();

  mesh->npmax = MG_MAX(1.5*mesh->np,_MMG2D_NPMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMG2D_NEMAX);
  mesh->namax = M_MAX(1.5*mesh->na,_MMG2D_NEDMAX);
  mesh->xpmax  = M_MAX(0.1*mesh->xp,0.1*_MMG2D_NPMAX);

  if ( mesh->info.mem <= 0 ) {
    if ( mesh->memMax && (mesh->memMax >2000*million))
      /* maximal memory = 2Go */
      mesh->memMax = 2000*million;
    else {
      /* default value = 800 Mo */
      printf(" ## Maximum memory set to default value: %d Mo.\n",_MMG5_MEMMAX);
      mesh->memMax = _MMG5_MEMMAX*million;
    }
  }
  else {
    /* memory asked by user if possible, otherwise total physical memory */
    if ( (long long)(mesh->info.mem)*million > mesh->memMax && mesh->memMax ) {
      fprintf(stdout,"  ## Warning: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"when only %ld available.\n",_MMG5_safeLL2LCast((long long)(mesh->memMax/million)));
    }
    else {
      mesh->memMax= (long long)(mesh->info.mem)*million;
    }

    /* if asked memory is lower than default _MMG2D_NPMAX/_MMG2D_NTMAX we take lower values */
    ctri = 2;

    /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
     * point+tria+tets+adja+adjt+sol+item *
     * warning: we exceed memory in saveMesh when we call _MMG5_hNew */
    bytes = sizeof(MMG5_Point) +  0.1*sizeof(MMG5_xPoint) +
      2*sizeof(MMG5_Tria) + 3*sizeof(int)
      + sizeof(MMG5_Sol) /*+ sizeof(Displ)*/
      + sizeof(int) + 5*sizeof(int);

    /*init allocation need 38Mo*/
    npask = (int)((double)(mesh->info.mem-38) / bytes * (int)million);
    mesh->npmax = MG_MIN(npask,mesh->npmax);
    mesh->ntmax = MG_MIN(ctri*npask,mesh->ntmax);
    mesh->namax = MG_MIN(ctri*npask,mesh->namax);
    mesh->xpmax = MG_MIN(0.1*npask,0.1*mesh->xp);

    /*check if the memory asked is enough to load the mesh*/
    if(mesh->np &&
       (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->namax < mesh->na) ){
      memtmp = (int)(mesh->np * bytes /(int)million + 38);
      memtmp = MG_MAX(memtmp, (int)(mesh->nt * bytes /(ctri* (int)million) + 38));
      memtmp = MG_MAX(memtmp, (int)(mesh->na * bytes /(ctri*(int)million) + 38));
      mesh->memMax = (long long) memtmp+1;
      fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
              memtmp+1);
    }
    if(mesh->info.mem < 39) {
      mesh->memMax = (long long) 39;
      fprintf(stdout,"  ## ERROR: asking for %d Mo of memory ",mesh->info.mem);
      fprintf(stdout,"is not enough to load mesh. You need to ask %d Mo minimum\n",
              39);
    }
  }

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (Mo)    %ld\n",
            _MMG5_safeLL2LCast((long long)(mesh->memMax/million)));

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  _MMG2D_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  _MMG2D_NTMAX    %d\n",mesh->ntmax);
    fprintf(stdout,"  _MMG2D_NAMAX    %d\n",mesh->namax);
  }

  return;
}

/* allocate main structure */
int MMG2D_zaldy(MMG5_pMesh mesh) {
  int     k;

  _MMG2D_memOption(mesh);

  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point);

  if ( mesh->xp ) {
    _MMG5_ADD_MEM(mesh,(mesh->xpmax+1)*sizeof(MMG5_xPoint),"initial xpoint",return(0));
    _MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xpmax+1,MMG5_xPoint);
    memset(&mesh->xpoint[0],0,sizeof(MMG5_xPoint));
  }
  _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",return(0));
  _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria);
  memset(&mesh->tria[0],0,sizeof(MMG5_Tria));

  _MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"initial edges",return(0));
  _MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge);

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;
  mesh->nanil = mesh->na + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    /* Set tangent field of point to 0 */
    mesh->point[k].n[0] = 0;
    mesh->point[k].n[1] = 0;
    mesh->point[k].n[2] = 0;
    /* link */
    mesh->point[k].tmp  = k+1;
  }

  for (k=mesh->nanil; k<mesh->namax-1; k++)
    mesh->edge[k].b = k+1;

  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return(1);
}

