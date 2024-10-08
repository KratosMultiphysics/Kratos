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
 * \file mmg3d/zaldy_3d.c
 * \brief Memory management.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmg3d_private.h"

/** get new point address */
MMG5_int MMG3D_newPt(MMG5_pMesh mesh,double c[3],int16_t tag,MMG5_int src) {

  MMG5_pPoint  ppt;
  MMG5_int     curpt;

  if ( !mesh->npnil )  return 0;
  curpt = mesh->npnil;
  if ( mesh->npnil > mesh->np )  mesh->np = mesh->npnil;
  ppt   = &mesh->point[curpt];
  memcpy(ppt->c,c,3*sizeof(double));
  mesh->npnil = ppt->tmp;
  ppt->tmp    = 0;

  ppt->ref = 0;
  ppt->xp = 0;
  ppt->flag = 0;
  /* point on geometry */
  if ( tag & MG_BDY ) {
    mesh->xp++;
    if(mesh->xp > mesh->xpmax){
      /* reallocation of xpoint table */
      MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                         "larger xpoint table",
                         return 0);
    }
    ppt->xp  = mesh->xp;
  }
  assert(tag < 24704);
  assert(tag >= 0);
  ppt->n[0]   = 0;
  ppt->n[1]   = 0;
  ppt->n[2]   = 0;
  ppt->tag    = tag;
  ppt->tagdel = 0;
#ifdef USE_POINTMAP
  assert( src );
  ppt->src = src;
#endif
  return curpt;
}

void MMG3D_delPt(MMG5_pMesh mesh,MMG5_int ip) {
  MMG5_pPoint   ppt;
  MMG5_xPoint   *pxp;

  ppt = &mesh->point[ip];
  if ( ppt->xp ) {
    pxp = &mesh->xpoint[ppt->xp];
    memset(pxp,0,sizeof(MMG5_xPoint));
  }
  memset(ppt,0,sizeof(MMG5_Point));
  ppt->tag    = MG_NUL;
  ppt->tmp    = mesh->npnil;
  mesh->npnil = ip;
  if ( ip == mesh->np ) {
    while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;
  }
}

/** get new elt address */
MMG5_int MMG3D_newElt(MMG5_pMesh mesh) {
  MMG5_int     curiel;

  if ( !mesh->nenil )  return 0;
  curiel = mesh->nenil;

  if ( mesh->nenil > mesh->ne )  mesh->ne = mesh->nenil;
  mesh->nenil = mesh->tetra[curiel].v[3];
  mesh->tetra[curiel].v[3] = 0;
  mesh->tetra[curiel].mark=mesh->mark;

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
int MMG3D_delElt(MMG5_pMesh mesh,MMG5_int iel) {
  MMG5_pTetra   pt;
  MMG5_int      iadr;

  pt = &mesh->tetra[iel];
  if ( !MG_EOK(pt) ) {
    fprintf(stderr,"\n  ## INVALID ELEMENT %" MMG5_PRId ".\n",iel);
    return 0;
  }
  memset(pt,0,sizeof(MMG5_Tetra));
  pt->v[3] = mesh->nenil;
  iadr = 4*(iel-1) + 1;
  if ( mesh->adja )
    memset(&mesh->adja[iadr],0,4*sizeof(MMG5_int));
  mesh->nenil = iel;
  if ( iel == mesh->ne ) {
    while ( !MG_EOK((&mesh->tetra[mesh->ne])) )  mesh->ne--;
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
int MMG3D_memOption_memSet(MMG5_pMesh mesh) {

  MMG5_memOption_memSet(mesh);

  return  MMG3D_memOption_memRepartition( mesh );
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition for the memMax amout of memory available.
 *
 */
int MMG3D_memOption_memRepartition(MMG5_pMesh mesh) {
  size_t     usedMem,avMem,reservedMem,npadd;
  int        ctri;
  MMG5_int   bytes;

  /* init allocation need MMG5_MEMMIN B */
  reservedMem = MMG5_MEMMIN +
    (mesh->nprism*sizeof(MMG5_Prism)+mesh->nquad*sizeof(MMG5_Quad));

  /* Compute the needed initial memory */
  usedMem = reservedMem + (mesh->np+1)*sizeof(MMG5_Point)
    + (mesh->nt+1)*sizeof(MMG5_Tria) + (mesh->ne+1)*sizeof(MMG5_Tetra)
    + (3*mesh->nt+1)*sizeof(MMG5_int)   + (4*mesh->ne+1)*sizeof(MMG5_int)
    + (mesh->np+1)*sizeof(double);

  if ( usedMem > mesh->memMax  ) {
    fprintf(stderr,"\n  ## Error: %s: %zu MB of memory ",__func__,mesh->memMax/MMG5_MILLION);
    fprintf(stderr,"is not enough to load mesh. You need to ask %zu MB minimum\n",
            usedMem/MMG5_MILLION+1);
    return 0;
  }

  ctri = 2;
  /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
   * point+tria+tets+adja+adjt+aniso sol+item */
  bytes = sizeof(MMG5_Point) + sizeof(MMG5_xPoint) +
    6*sizeof(MMG5_Tetra) + ctri*sizeof(MMG5_xTetra) +
    4*6*sizeof(MMG5_int) + ctri*3*sizeof(MMG5_int) +
    4*sizeof(MMG5_hedge)+6*sizeof(double);

#ifdef USE_SCOTCH
  /* bytes = bytes + vertTab + edgeTab + PermVrtTab *
   * + vertOldTab + sortPartTab - adja */
  bytes = bytes + 3*6*sizeof(MMG5_int);
#endif

  avMem = mesh->memMax-usedMem;

  /* If npadd is exactly the maximum memory available, we will use all the
   * memory and the analysis step will fail. As arrays may be reallocated, we
   * can have smaller values for npmax,ntmax and nemax (npadd/2). */
  npadd = avMem/(2*bytes);
  mesh->npmax = MG_MIN(mesh->npmax,mesh->np+npadd);
  mesh->ntmax = MG_MIN(mesh->ntmax,ctri*npadd+mesh->nt);
  mesh->nemax = MG_MIN(mesh->nemax,6*npadd+mesh->ne);

  if ( sizeof(MMG5_int) == sizeof(int32_t) ) {
    /** Check that we will not overflow int32_max when allocating adja array */

    int coef;
    if ( mesh->nprism ) {
      coef = 5;
    }
    else {
      coef = 4;
    }

    /* maximal number of triangles, taking the
     * computation of adjacency relationships into account */
    int32_t int32_nemax = (INT32_MAX-(coef+1))/coef;

    if ( int32_nemax < mesh->nemax ) {
      if ( int32_nemax <= mesh->ne ) {
        /* No possible allocation without int32 overflow */
        fprintf(stderr,"\n  ## Error: %s: with %" MMG5_PRId " tetrahedra Mmg will overflow"
                " the 32-bit integer.\n",__func__,mesh->ne);
        fprintf(stderr,"Please, configure Mmg with MMG5_INT=int64_t argument.\n");
        return 0;
      }
      else {
        /* Correction of maximal number of tetrahedra */
        mesh->nemax = int32_nemax;
      }
    }
  }

  /* check if the memory asked is enough to load the mesh*/
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (MB)    %zu\n",
            mesh->memMax/MMG5_MILLION);
  }
  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MMG3D_NPMAX    %" MMG5_PRId "\n",mesh->npmax);
    fprintf(stdout,"  MMG3D_NTMAX    %" MMG5_PRId "\n",mesh->ntmax);
    fprintf(stdout,"  MMG3D_NEMAX    %" MMG5_PRId "\n",mesh->nemax);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition for the -m option.
 *
 */
int MMG3D_memOption(MMG5_pMesh mesh) {

  mesh->npmax = MG_MAX((int)(1.5*mesh->np),MMG3D_NPMAX);
  mesh->nemax = MG_MAX((int)(1.5*mesh->ne),MMG3D_NEMAX);
  mesh->ntmax = MG_MAX((int)(1.5*mesh->nt),MMG3D_NTMAX);

  return  MMG3D_memOption_memSet(mesh);
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh.
 *
 */
int MMG3D_setMeshSize_alloc( MMG5_pMesh mesh ) {
  MMG5_int k;

  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,return 0);

  MMG5_ADD_MEM(mesh,(mesh->nemax+1)*sizeof(MMG5_Tetra),"initial tetrahedra",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->tetra,mesh->nemax+1,MMG5_Tetra,return 0);

  if ( mesh->nprism ) {
    MMG5_ADD_MEM(mesh,(mesh->nprism+1)*sizeof(MMG5_Prism),"initial prisms",return 0);
    MMG5_SAFE_CALLOC(mesh->prism,(mesh->nprism+1),MMG5_Prism,return 0);
  }

  if ( mesh->nt ) {
    MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_Tria),"initial triangles",return 0);
    MMG5_SAFE_CALLOC(mesh->tria,mesh->nt+1,MMG5_Tria,return 0);
    memset(&mesh->tria[0],0,sizeof(MMG5_Tria));
  }

  if ( mesh->nquad ) {
    MMG5_ADD_MEM(mesh,(mesh->nquad+1)*sizeof(MMG5_Quad),"initial quadrilaterals",return 0);
    MMG5_SAFE_CALLOC(mesh->quadra,(mesh->nquad+1),MMG5_Quad,return 0);
  }

  mesh->namax = mesh->na;
  if ( mesh->na ) {
    MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",return 0);
    MMG5_SAFE_CALLOC(mesh->edge,(mesh->na+1),MMG5_Edge,return 0);
  }

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->ne + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    /* link */
    mesh->point[k].tmp  = k+1;
  }

  for (k=mesh->nenil; k<mesh->nemax-1; k++)
    mesh->tetra[k].v[3] = k+1;

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
int MMG3D_zaldy(MMG5_pMesh mesh) {

  if ( !MMG3D_memOption(mesh) )  return 0;

  return  MMG3D_setMeshSize_alloc(mesh);
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Free xtetra structure.
 *
 */
void MMG5_freeXTets(MMG5_pMesh mesh) {
  MMG5_pTetra pt;
  MMG5_int    k;

  for (k=1; k<=mesh->ne; k++) {
    pt     = &mesh->tetra[k];
    pt->xt = 0;
  }
  if ( mesh->xtetra )
    MMG5_DEL_MEM(mesh,mesh->xtetra);
  mesh->xt = 0;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Free xprism structure.
 *
 */
void MMG5_freeXPrisms(MMG5_pMesh mesh) {
  MMG5_pPrism pp;
  MMG5_int    k;

  for (k=1; k<=mesh->nprism; k++) {
    pp      = &mesh->prism[k];
    pp->xpr = 0;
  }
  if ( mesh->xprism )
    MMG5_DEL_MEM(mesh,mesh->xprism);
  mesh->xpr = 0;
}
