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
 * \file mmgs/zaldy_s.c
 * \brief Memory management.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmgs_private.h"

/* get new point address */
MMG5_int MMGS_newPt(MMG5_pMesh mesh,double c[3],double n[3]) {
  MMG5_pPoint  ppt;
  MMG5_int     curpt;

  if ( !mesh->npnil )  return 0;

  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,3*sizeof(double));
  if ( n )
    memcpy(ppt->n,n,3*sizeof(double));
  ppt->tag   &= ~MG_NUL;
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;

  return curpt;
}

void MMGS_delPt(MMG5_pMesh mesh,MMG5_int ip) {
  MMG5_pPoint   ppt;

  ppt = &mesh->point[ip];
  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;
  mesh->npnil = ip;
  if ( ip == mesh->np ) {
    while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;
  }
}

MMG5_int MMGS_newElt(MMG5_pMesh mesh) {
  MMG5_int     curiel;

  if ( !mesh->nenil )  return 0;
  curiel = mesh->nenil;

  if ( mesh->nenil > mesh->nt )  mesh->nt = mesh->nenil;
  mesh->nenil = mesh->tria[curiel].v[2];
  mesh->tria[curiel].v[2] = 0;

  return curiel;
}

/**
 * \param mesh pointer toward the mesh
 * \param iel index of the element to delete
 *
 * \return 1 if success, 0 if fail
 *
 * Delete the element \a iel
 *
 */
int MMGS_delElt(MMG5_pMesh mesh,MMG5_int iel) {
  MMG5_pTria    pt;

  pt = &mesh->tria[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stderr,"\n  ## INVALID ELEMENT %" MMG5_PRId ".\n",iel);
    return 0;
  }
  memset(pt,0,sizeof(MMG5_Tria));
  pt->v[2] = mesh->nenil;
  if ( mesh->adja )
    memset(&mesh->adja[3*(iel-1)+1],0,3*sizeof(MMG5_int));
  mesh->nenil = iel;
  if ( iel == mesh->nt ) {
    while ( !MG_EOK((&mesh->tria[mesh->nt])) )  mesh->nt--;
  }
  return 1;
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
 * \remark Here, mesh->npmax/nemax/ntmax must be setted.
 *
 */
static inline
int MMGS_memOption_memSet(MMG5_pMesh mesh) {
  size_t     usedMem,avMem,reservedMem,npadd;
  int        bytes;

  MMG5_memOption_memSet(mesh);

  /* init allocation need MMG5_MEMMIN B */
  reservedMem = MMG5_MEMMIN;

  /* Compute the needed initial memory */
  usedMem = reservedMem + (mesh->np+1)*sizeof(MMG5_Point)
    + (mesh->nt+1)*sizeof(MMG5_Tria) + (3*mesh->nt+1)*sizeof(MMG5_int)
    + (mesh->np+1)*sizeof(double);

  if ( usedMem > mesh->memMax  ) {
    fprintf(stderr,"\n  ## Error: %s: %zu MB of memory ",__func__,mesh->memMax/MMG5_MILLION);
    fprintf(stderr,"is not enough to load mesh. You need to ask %zu MB minimum\n",
            usedMem/MMG5_MILLION+1);
    return 0;
  }

  /* point+xpoint+tria+adja+aniso sol */
  bytes = sizeof(MMG5_Point) + sizeof(MMG5_xPoint) +
    2*sizeof(MMG5_Tria) + 3*sizeof(MMG5_int) + 6*sizeof(double);

  avMem = mesh->memMax-usedMem;

  /* If npadd is exactly the maximum memory available, we will use all the
   * memory and the analysis step will fail. As arrays may be reallocated, we
   * can have smaller values for xpmax and ntmax (npadd/2). */
  npadd = avMem/(2*bytes);
  mesh->npmax = MG_MIN(mesh->npmax,mesh->np+npadd);
  mesh->ntmax = MG_MIN(mesh->ntmax,2*npadd+mesh->nt);

  if ( sizeof(MMG5_int) == sizeof(int32_t) ) {
    /** Check that we will not overflow int32_max when allocating adja array */
    /* maximal number of triangles, taking the
     * computation of adjacency relationships into account */
    int coef = 3;
    int32_t int32_ntmax = (INT32_MAX-(coef+1))/coef;

    if ( int32_ntmax < mesh->ntmax ) {
      if ( int32_ntmax <= mesh->nt ) {
        /* No possible allocation without int32 overflow */
        fprintf(stderr,"\n  ## Error: %s: with %" MMG5_PRId " triangles Mmg will overflow"
                " the 32-bit integer.\n",__func__,mesh->nt);
        fprintf(stderr,"Please, configure Mmg with MMG5_INT=int64_t argument.\n");
        return 0;
      }
      else {
        /* Correction of maximal number of triangles */
        mesh->ntmax = int32_ntmax;
      }
    }
  }

  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (MB)    %zu\n",
            mesh->memMax/MMG5_MILLION);
  }

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MMG2D_NPMAX    %" MMG5_PRId "\n",mesh->npmax);
    fprintf(stdout,"  MMG2D_NTMAX    %" MMG5_PRId "\n",mesh->ntmax);
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
int MMGS_memOption(MMG5_pMesh mesh) {

  mesh->memMax = MMG5_memSize();

  mesh->npmax = MG_MAX(1.5*mesh->np,MMGS_NPMAX);
  mesh->ntmax = MG_MAX(1.5*mesh->nt,MMGS_NTMAX);

  return  MMGS_memOption_memSet(mesh);
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh.
 *
 */
int MMGS_setMeshSize_alloc( MMG5_pMesh mesh ) {
  MMG5_int k;

  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,return 0);
  MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria,return 0);


  mesh->namax = mesh->na;
  if ( mesh->na ) {
    MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",return 0);
    MMG5_SAFE_CALLOC(mesh->edge,(mesh->na+1),MMG5_Edge,return 0);
  }

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++)
    mesh->point[k].tmp  = k+1;

  for (k=mesh->nenil; k<mesh->ntmax-1; k++)
    mesh->tria[k].v[2] = k+1;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 1 if success, 0 if fail
 *
 * allocate main structure
 *
 */
int MMGS_zaldy(MMG5_pMesh mesh) {

  if ( !MMGS_memOption(mesh) )  return 0;

  return  MMGS_setMeshSize_alloc(mesh);
}
