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


/* Create a new vertex in the mesh, and return its number */
int MMG2D_newPt(MMG5_pMesh mesh,double c[2],int16_t tag) {
  MMG5_pPoint  ppt;
  int     curpt;

  if ( !mesh->npnil )  return 0;

  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,2*sizeof(double));
  ppt->tag   &= ~MG_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;
  ppt->tag = tag;

  return curpt;
}

/* Delete a point in the mesh and update the garbage collector accordingly */
void MMG2D_delPt(MMG5_pMesh mesh,int ip) {
  MMG5_pPoint   ppt;

  ppt = &mesh->point[ip];

  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;

  mesh->npnil = ip;
  if ( ip == mesh->np )  mesh->np--;
}

void MMG5_delEdge(MMG5_pMesh mesh,int iel) {
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

/* Create a new triangle in the mesh and return its address */
int MMG2D_newElt(MMG5_pMesh mesh) {
  int     curiel;

  if ( !mesh->nenil ) {
    return 0;
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

  return curiel;
}

/* Delete a triangle in the mesh and update the garbage collector accordingly */
int MMG2D_delElt(MMG5_pMesh mesh,int iel) {
  MMG5_pTria    pt;
  int      iadr;

  pt = &mesh->tria[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stdout,"  ## INVALID ELEMENT.\n");
    return 0;
  }
  memset(pt,0,sizeof(MMG5_Tria));
  pt->v[2] = mesh->nenil;
  pt->qual = 0.0;
  iadr = (iel-1)*3 + 1;
  if ( mesh->adja )
    memset(&mesh->adja[iadr],0,3*sizeof(int));

  mesh->nenil = iel;
  if ( iel == mesh->nt )  mesh->nt--;
  return 1;
}


/* check if n elets available */
int MMG5_getnElt(MMG5_pMesh mesh,int n) {
  int     curiel;

  if ( !mesh->nenil )  return 0;
  curiel = mesh->nenil;
  do {
    curiel = mesh->tria[curiel].v[2];
  }
  while (--n);

  return n == 0;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * Set the memMax value to its "true" value (50% of the RAM or memory asked by
 * user) and perform memory repartition for the -m option.  If -m is not given,
 * memMax is the detected RAM. If -m is provided, check the user option and set
 * memMax to the available RAM if the user ask for too much memory. Last,
 * perform the memory repartition between the mmg arrays with respect to the
 * memMax value.
 *
 * \remark Here, mesh->npmax/ntmax must be setted.
 *
 */
static inline
int MMG2D_memOption_memSet(MMG5_pMesh mesh) {
  size_t   usedMem,avMem,reservedMem,npadd;
  int      ctri,bytes;

  MMG5_memOption_memSet(mesh);

  /* init allocation need MMG5_MEMMIN B */
  reservedMem = MMG5_MEMMIN;

  /* Compute the needed initial memory */
  usedMem = reservedMem + (mesh->np+1)*sizeof(MMG5_Point)
    + (mesh->nt+1)*sizeof(MMG5_Tria) + (3*mesh->nt+1)*sizeof(int)
    + (mesh->na+1)*sizeof(MMG5_Edge) + (mesh->np+1)*sizeof(double);

  if ( usedMem > mesh->memMax  ) {
    fprintf(stderr,"\n  ## Error: %s: %zu MB of memory ",__func__,
            mesh->memMax/MMG5_MILLION);
    fprintf(stderr,"is not enough to load mesh. You need to ask %zu MB minimum\n",
            usedMem/MMG5_MILLION+1);
    return 0;
  }

  ctri = 2;

  /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
   * point+tria+edges+adjt+ aniso sol */
  bytes = sizeof(MMG5_Point) +
    2*sizeof(MMG5_Tria) + 3*2*sizeof(int)
    + 0.2*sizeof(MMG5_Edge) + 3*sizeof(double);

  avMem = mesh->memMax-usedMem;

  /* If npadd is exactly the maximum memory available, we will use all the
   * memory and the analysis step will fail. As arrays may be reallocated, we
   * can have smaller values for npmax and ntmax (npadd/2). */
  npadd = avMem/(double)(2*bytes);
  mesh->npmax = MG_MIN(mesh->npmax,mesh->np+npadd);
  mesh->ntmax = MG_MIN(mesh->ntmax,ctri*npadd+mesh->nt);
  mesh->namax = MG_MIN(mesh->namax,ctri*npadd+mesh->na);

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (MB)    %zu\n",
            mesh->memMax/MMG5_MILLION);
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MMG2D_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  MMG2D_NTMAX    %d\n",mesh->ntmax);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition for the -m option
 *
 */
int MMG2D_memOption(MMG5_pMesh mesh) {

  mesh->memMax = MMG5_memSize();

  mesh->npmax = MG_MAX(1.5*mesh->np,MMG2D_NPMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,MMG2D_NEMAX);
  mesh->namax = mesh->na;

  return  MMG2D_memOption_memSet(mesh);
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh.
 *
 */
int MMG2D_setMeshSize_alloc( MMG5_pMesh mesh ) {
  int k;

  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                printf("  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,return 0);

  MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",return 0);
  MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria,return 0);
  memset(&mesh->tria[0],0,sizeof(MMG5_Tria));

  mesh->namax = mesh->na;
  if ( mesh->na ) {
    MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"initial edges",return 0);
    MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge,return 0);
  }

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;
  mesh->nanil = 0;

  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    /* Set tangent field of point to 0 */
    mesh->point[k].n[0] = 0;
    mesh->point[k].n[1] = 0;
    mesh->point[k].n[2] = 0;
    /* link */
    mesh->point[k].tmp  = k+1;
  }

  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * allocate main structure
 *
 */
int MMG2D_zaldy(MMG5_pMesh mesh) {

  if ( !MMG2D_memOption(mesh) )  return 0;

  return  MMG2D_setMeshSize_alloc(mesh);
}
