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
 * \file mmg3d/libmmg3d.c
 * \brief Most of the API functions of the MMG3D library.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete soon...
 * \todo documentation doxygen
 *
 * Private API functions for MMG3D library: incompatible functions
 * with the main binary.
 *
 */

#include "libmmg3d.h"
#include "inlined_functions_3d_private.h"
#include "mmg3dexterns_private.h"
#include "mmgexterns_private.h"

/**
 * Pack the mesh \a mesh and its associated metric \a met and/or solution \a sol
 * and return \a val.
 */
#define MMG5_RETURN_AND_PACK(mesh,met,sol,val)do  \
  {                                               \
    if ( !MMG3D_packMesh(mesh,met,sol) )  {       \
      mesh->npi = mesh->np;                       \
      mesh->nti = mesh->nt;                       \
      mesh->nai = mesh->na;                       \
      mesh->nei = mesh->ne;                       \
      mesh->xt   = 0;                             \
      if ( met ) { met->npi  = met->np; }         \
      if ( sol ) { sol->npi  = sol->np; }         \
      return MMG5_STRONGFAILURE;                  \
    }                                             \
    _LIBMMG5_RETURN(mesh,met,sol,val);            \
  }while(0)

/** Free adja, xtetra and xpoint tables */
void MMG3D_Free_topoTables(MMG5_pMesh mesh) {
  MMG5_int k;

  mesh->xp = 0;
  if ( mesh->adja )
    MMG5_DEL_MEM(mesh,mesh->adja);

  MMG5_freeXTets(mesh);

  if ( mesh->adjapr )
    MMG5_DEL_MEM(mesh,mesh->adjapr);

  MMG5_freeXPrisms(mesh);

  if ( mesh->xpoint )
    MMG5_DEL_MEM(mesh,mesh->xpoint);

  for(k=1; k <=mesh->np; k++) {
    mesh->point[k].xp = 0;
  }

  return;
}

/**
 * \param mesh pointer to the mesh structure (unused).
 * \return -1 if fail, the number of detected ridges otherwise
 *
 * Create the boundary entities of the mesh (triangles and edges).
 *
 * \warning mesh must be packed and hashed
 */
MMG5_int MMG3D_bdryBuild(MMG5_pMesh mesh) {
  MMG5_pTetra pt;
  MMG5_hgeom  *ph;
  MMG5_int    k,nr;
  int         i;


  /* rebuild triangles*/
  if ( mesh->tria )
    MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to rebuild triangles\n",__func__);
    return -1;
  }

  /* build hash table for edges */
  if ( mesh->htab.geom )
    MMG5_DEL_MEM(mesh,mesh->htab.geom);

  if ( mesh->edge )
    MMG5_DEL_MEM(mesh,mesh->edge);
  mesh->na = 0;

  nr = 0;
  /* in the worst case (all edges are marked), we will have around 1 edge per *
   * triangle (we count edges only one time) */
  if ( MMG5_hNew(mesh,&mesh->htab,mesh->nt,3*(mesh->nt)) ) {
    for (k=1; k<=mesh->ne; k++) {
      pt   = &mesh->tetra[k];
      if ( MG_EOK(pt) &&  pt->xt ) {
        for (i=0; i<6; i++) {
          if ( mesh->xtetra[pt->xt].edg[i] ||
               ( mesh->xtetra[pt->xt].tag[i] & MG_REQ ||
                 MG_EDG(mesh->xtetra[pt->xt].tag[i])) )
            if ( !MMG5_hEdge(mesh,&mesh->htab,pt->v[MMG5_iare[i][0]],pt->v[MMG5_iare[i][1]],
                             mesh->xtetra[pt->xt].edg[i],mesh->xtetra[pt->xt].tag[i]))
              return -1;
        }
      }
    }

    /* edges + ridges + required edges */
    for (k=0; k<=mesh->htab.max; k++) {
      ph = &mesh->htab.geom[k];
      if ( !(ph->a) )  continue;
      mesh->na++;
    }
    if ( mesh->na ) {
      MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"edges",
                   mesh->na = 0;
                   printf("  ## Warning: uncomplete mesh\n"));
    }

    if ( mesh->na ) {
      MMG5_SAFE_CALLOC(mesh->edge,mesh->na+1,MMG5_Edge,return -1);

      mesh->na = 0;
      for (k=0; k<=mesh->htab.max; k++) {
        ph = &mesh->htab.geom[k];
        if ( !ph->a )  continue;
        mesh->na++;
        mesh->edge[mesh->na ].a  = ph->a;
        mesh->edge[mesh->na ].b  = ph->b;
        mesh->edge[mesh->na].tag = ( ph->tag | MG_REF ) ;
        mesh->edge[mesh->na].ref = ph->ref;

        if ( MG_GEO & ph->tag ) nr++;
      }
    }
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
  }
  else
    mesh->memCur -= (long long)((3*mesh->nt+2)*sizeof(MMG5_hgeom));

  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  if ( mesh->info.imprim > 0 ) {
    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8" MMG5_PRId "   RIDGES  %8" MMG5_PRId "\n",mesh->na,nr);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES  %8" MMG5_PRId "\n",mesh->nt);
  }

  return nr;
}

#ifndef USE_SCOTCH
/**
 * \param mesh pointer to the mesh structure (unused).
 * \param np pointer to the number of used points
 * \param nc pointer to the number of packed corners
 * \return 1 if success, 0 if fail.
 *
 * Pack the mesh points and store the packed point index in tmp.
 * Don't preserve numbering order.
 *
 */

int MMG3D_mark_packedPoints(MMG5_pMesh mesh,MMG5_int *np,MMG5_int *nc) {
  MMG5_pPoint   ppt;
  MMG5_int      k;

  for ( k=1; k<=mesh->np; ++k ) {
    mesh->point[k].tmp = 0;
  }

  (*nc) = 0;
  k     = 1;
  (*np) = mesh->np+1;

  do {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) {
      --(*np);
      ppt = &mesh->point[*np];
      assert ( ppt );

      /* Search the last used point */
      while ( !MG_VOK(ppt) && k < *np ) {
        --(*np);
        ppt = &mesh->point[*np];
      }
    }

#ifndef NDEBUG
    if ( !MG_VOK(ppt) ) {
      assert ( k==*np );
    }
    if ( k==*np ) {
      assert ( !MG_VOK(ppt) );
    }
#endif
    if ( k==*np ) {
      break;
    }

    ppt->tmp = k;

    if ( ppt->tag & MG_NOSURF ) {
      ppt->tag &= ~MG_NOSURF;
      ppt->tag &= ~MG_REQ;
    }

    if ( ppt->tag & MG_CRN )  (*nc)++;

    ppt->ref = MMG5_abs(ppt->ref);
  }
  while ( ++k < (*np) );

  assert ( k==*np );
  --(*np);

#ifndef NDEBUG
  for ( k=1; k<=mesh->np; ++k  ) {
    if ( MG_VOK(&mesh->point[k]) ) {
      /* Check that all the used points have a tmp field  */
      assert(mesh->point[k].tmp);
    }
  }
#endif

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Pack the sparse tetra and the associated adja array.
 * Don't preserve numbering order.
 *
 */
int MMG3D_pack_tetraAndAdja(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  MMG5_int      k,iadr,iadr1,iadrv,*adjav,*adja,*adja1;
  int           i,voy;

  if ( !mesh->ne ) {
    return 1;
  }

  k = 1;
  do {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) {
      pt1 = &mesh->tetra[mesh->ne];
      assert( pt && pt1 && MG_EOK(pt1) );
      memcpy(pt,pt1,sizeof(MMG5_Tetra));

      /* treat adja array */
      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadr1 = 4*(mesh->ne-1) + 1;
      adja1 = &mesh->adja[iadr1];
      for(i=0 ; i<4 ; i++) {
        adja[i] = adja1[i];
        if(!adja1[i]) continue;
        iadrv = 4*(adja1[i]/4-1) + 1;
        adjav = &mesh->adja[iadrv];
        voy = i;
        adjav[adja1[i]%4] = 4*k + voy;
      }

      if ( !MMG3D_delElt(mesh,mesh->ne) )  return 0;
    }
  }
  while ( ++k < mesh->ne );

  /* Recreate nil chain */
  if ( mesh->ne >= mesh->nemax-1 )
    mesh->nenil = 0;
  else
    mesh->nenil = mesh->ne + 1;

  mesh->nei = mesh->ne;

  if ( mesh->nenil )
    for(k=mesh->nenil; k<mesh->nemax-1; k++)
      mesh->tetra[k].v[3] = k+1;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Pack the sparse tetra. Doesn't pack the adjacency array.
 * Don't preserve numbering order.
 *
 */
int MMG3D_pack_tetra(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,pt1;
  MMG5_int      k;

  if ( !mesh->ne ) {
    return 1;
  }

  if ( mesh->tetra ) {
    k = 1;
    do {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) {
        pt1 = &mesh->tetra[mesh->ne];
        assert( pt && pt1 && MG_EOK(pt1) );
        memcpy(pt,pt1,sizeof(MMG5_Tetra));
        if ( !MMG3D_delElt(mesh,mesh->ne) )  return 0;
      }
    }
    while ( ++k < mesh->ne );

    /* Recreate nil chain */
    if ( mesh->ne >= mesh->nemax-1 )
      mesh->nenil = 0;
    else
      mesh->nenil = mesh->ne + 1;

    if ( mesh->nenil )
      for(k=mesh->nenil; k<mesh->nemax-1; k++)
        mesh->tetra[k].v[0] = 0;
  }
  mesh->nei = mesh->ne;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Pack prisms and quads. Don't preserve numbering order.
 *
 */
int MMG3D_pack_prismsAndQuads(MMG5_pMesh mesh) {
  MMG5_pPrism   pp,pp1;
  MMG5_pQuad    pq,pq1;
  MMG5_int      k;

  if ( mesh->prism ) {
    k = 1;
    do {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) ) {
        pp1 = &mesh->prism[mesh->nprism];
        assert( pp && pp1 && MG_EOK(pp1) );
        memcpy(pp,pp1,sizeof(MMG5_Prism));
        --mesh->nprism;
      }
    }
    while ( ++k < mesh->nprism );
  }

  if ( mesh->quadra ) {
    k = 1;
    do {
      pq = &mesh->quadra[k];
      if ( !MG_EOK(pq) ) {
        pq1 = &mesh->quadra[mesh->nquad];
        assert( pq && pq1 && MG_EOK(pq1) );
        memcpy(pq,pq1,sizeof(MMG5_Quad));
        --mesh->nquad;
      }
    }
    while ( ++k < mesh->nquad );
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure (unused).
 * \param met pointer to the solution (metric or level-set) structure.
 * \return 1 if success, 0 if fail.
 *
 * Pack a sparse solution structure. Don't preserve numbering order.
 *
 */
int MMG3D_pack_sol(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pPoint   ppt,ppt1;
  MMG5_int      np,k,isol,isol1;
  int           i;

  if ( sol && sol->m ) {
    k  = 1;
    np = mesh->np;
    /* Search the last used point */
    while ( (!MG_VOK(mesh->point+np)) && np ) {
      --np;
    }

    do {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) {
        /* Copy last used point into first used one */
        isol    = k  * sol->size;
        isol1   = np * sol->size;

        assert( sol->m+isol && sol->m+isol1 && MG_VOK(&mesh->point[np]) );

        for (i=0; i<sol->size; i++) {
          sol->m[isol + i] = sol->m[isol1 + i];
        }

        /* Search the last used point */
        do {
          --np;
          ppt1 = &mesh->point[np];
        }
        while ( !MG_VOK(ppt1) && k < np );
      }
    }
    while ( ++k < np );

    sol->npi = sol->np = np;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure (unused).
 * \return the number of corners if success, -1 otherwise
 *
 * Pack a sparse point array. Don't preserve numbering order.
 *
 */
int MMG3D_pack_pointArray(MMG5_pMesh mesh) {
  MMG5_pPoint   ppt,ppt1;
  MMG5_int      k;

  k = 1;

  while ( (!MG_VOK(&mesh->point[mesh->np]))  && mesh->np  ) {
    MMG3D_delPt(mesh,mesh->np);
  }

  do {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) {
      ppt1 = &mesh->point[mesh->np];
      assert( ppt && ppt1 && MG_VOK(ppt1) );
      memcpy(ppt,ppt1,sizeof(MMG5_Point));

      /* delete point without deleting xpoint */
      memset(ppt1,0,sizeof(MMG5_Point));
      ppt1->tag    = MG_NUL;

      while ( !MG_VOK((&mesh->point[mesh->np])) )  mesh->np--;

    }

    /* Copy the normal stored in the xpoint into ppt->n. */
    if ( ppt->tag & MG_BDY &&
         !(ppt->tag & MG_CRN || ppt->tag & MG_NOM || MG_EDG(ppt->tag)) ) {

      if ( ppt->xp && mesh->xpoint ) {
        memcpy(ppt->n,mesh->xpoint[ppt->xp].n1,3*sizeof(double));
        ++mesh->nc1;
      }
    }
  }
  while ( ++k < mesh->np );

  /* Recreate nil chain */
  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  if ( mesh->np >= mesh->npmax-1 )
    mesh->npnil = 0;
  else
    mesh->npnil = mesh->np + 1;

  if ( mesh->npnil )
    for(k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;

  /* Set mesh->npi to suitable value */
  mesh->npi = mesh->np;

  return 1;
}

#else

/**
 * \param mesh pointer to the mesh structure (unused).
 * \param np pointer to the number of packed points
 * \param nc pointer to the number of packed corners
 * \return 1 if success, 0 if fail.
 *
 * Count the number of packed points and store the packed point index in tmp.
 * Preserve numbering order.
 *
 */

int MMG3D_mark_packedPoints(MMG5_pMesh mesh,MMG5_int *np,MMG5_int *nc) {
  MMG5_pPoint   ppt;
  MMG5_int      k;

  (*np) = (*nc) = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++(*np);

    if ( ppt->tag & MG_NOSURF ) {
      ppt->tag &= ~MG_NOSURF;
      ppt->tag &= ~MG_REQ;
    }

    if ( ppt->tag & MG_CRN )  (*nc)++;

    ppt->ref = MMG5_abs(ppt->ref);
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Pack the sparse tetra and the associated adja array.
 * Preserve numbering order.
 *
 */
int MMG3D_pack_tetraAndAdja(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,ptnew;
  MMG5_int      ne,nbl,k,iadr,iadrnew,iadrv,*adjav,*adja,*adjanew;
  int           i,voy;

  ne  = 0;
  nbl = 1;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    ne++;
    if ( k!=nbl ) {
      ptnew = &mesh->tetra[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));

      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadrnew = 4*(nbl-1) + 1;
      adjanew = &mesh->adja[iadrnew];
      for(i=0 ; i<4 ; i++) {
        adjanew[i] = adja[i];
        if(!adja[i]) continue;
        iadrv = 4*(adja[i]/4-1) +1;
        adjav = &mesh->adja[iadrv];
        voy = i;
        adjav[adja[i]%4] = 4*nbl + voy;
      }
    }
    nbl++;
  }
  mesh->nei = mesh->ne = ne;

  /* Recreate nil chain */
  if ( mesh->ne >= mesh->nemax-1 )
    mesh->nenil = 0;
  else
    mesh->nenil = mesh->ne + 1;

  if ( mesh->nenil )
    for(k=mesh->nenil; k<mesh->nemax-1; k++)
      mesh->tetra[k].v[3] = k+1;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Pack the sparse tetra. Doesn't pack the adjacency array.
 * Preserve numbering order.
 *
 */
int MMG3D_pack_tetra(MMG5_pMesh mesh) {
  MMG5_pTetra   pt,ptnew;
  MMG5_int      ne,nbl,k;

  ne  = 0;
  nbl = 1;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    ne++;
    if ( k!=nbl ) {
      ptnew = &mesh->tetra[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));
    }
    nbl++;
  }
  mesh->nei = mesh->ne = ne;

  /* Recreate nil chain */
  if ( mesh->ne >= mesh->nemax-1 )
    mesh->nenil = 0;
  else
    mesh->nenil = mesh->ne + 1;

  if ( mesh->nenil )
    for(k=mesh->nenil; k<mesh->nemax-1; k++)
      mesh->tetra[k].v[0] = 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Pack prisms and quads. Preserve numbering order.
 *
 */
int MMG3D_pack_prismsAndQuads(MMG5_pMesh mesh) {
  MMG5_pPrism   pp,ppnew;
  MMG5_pQuad    pq,pqnew;
  MMG5_int      k,ne,nbl;

  ne  = 0;
  nbl = 1;
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) )  continue;

    ++ne;
    if ( k!=nbl ) {
      ppnew = &mesh->prism[nbl];
      memcpy(ppnew,pp,sizeof(MMG5_Prism));
    }
    ++nbl;
  }
  mesh->nprism = ne;

  ne  = 0;
  nbl = 1;
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;

    ++ne;
    if ( k!=nbl ) {
      pqnew = &mesh->quadra[nbl];
      memcpy(pqnew,pq,sizeof(MMG5_Quad));
    }
    ++nbl;
  }
  mesh->nquad = ne;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure (unused).
 * \param met pointer to the solution (metric or level-set) structure.
 * \return 1 if success, 0 if fail.
 *
 * Pack a sparse solution structure. Preserve numbering order.
 *
 */
int MMG3D_pack_sol(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pPoint   ppt;
  MMG5_int      k,isol,isolnew,np,nbl;
  int           i;

  np  = 0;
  nbl = 1;
  if ( sol && sol->m ) {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;

      ++np;

      if ( k!= nbl ) {
        isol    = k   * sol->size;
        isolnew = nbl * sol->size;

        for (i=0; i<sol->size; i++)
          sol->m[isolnew + i] = sol->m[isol + i];
      }
      ++nbl;
    }
    sol->npi = sol->np = np;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure (unused).
 * \return the number of corners if success, -1 otherwise
 *
 * Pack a sparse point array. Preserve numbering order.
 *
 */
int MMG3D_pack_pointArray(MMG5_pMesh mesh) {
  MMG5_pPoint   ppt,pptnew;
  MMG5_int      k,np,nbl;

  nbl       = 1;
  mesh->nc1 = 0;
  np        = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;

    if ( ppt->tag & MG_BDY &&
         !(ppt->tag & MG_CRN || ppt->tag & MG_NOM || MG_EDG(ppt->tag)) ) {

      if ( ppt->xp && mesh->xpoint ) {
        memcpy(ppt->n,mesh->xpoint[ppt->xp].n1,3*sizeof(double));
        ++mesh->nc1;
      }
    }

    np++;
    if ( k!=nbl ) {
      pptnew = &mesh->point[nbl];
      memmove(pptnew,ppt,sizeof(MMG5_Point));
      memset(ppt,0,sizeof(MMG5_Point));
      ppt->tag    = MG_NUL;
    }

    nbl++;
  }
  mesh->npi = mesh->np = np;

  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  if ( mesh->np >= mesh->npmax-1 )
    mesh->npnil = 0;
  else
    mesh->npnil = mesh->np + 1;

  if ( mesh->npnil )
    for(k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;

  return 1;
}

#endif


/**
 * \param mesh pointer to the mesh structure (unused).
 * \return 1 if success, 0 otherwise
 *
 * Update the element vertices indices with the pack point index stored in the
 * tmp field of the points.
 *
 */
int MMG3D_update_eltsVertices(MMG5_pMesh mesh) {
  MMG5_pTetra   pt;
  MMG5_pPrism   pp;
  MMG5_pQuad    pq;
  MMG5_int      k;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    pt->v[3] = mesh->point[pt->v[3]].tmp;
  }
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) )  continue;

    pp->v[0] = mesh->point[pp->v[0]].tmp;
    pp->v[1] = mesh->point[pp->v[1]].tmp;
    pp->v[2] = mesh->point[pp->v[2]].tmp;
    pp->v[3] = mesh->point[pp->v[3]].tmp;
    pp->v[4] = mesh->point[pp->v[4]].tmp;
    pp->v[5] = mesh->point[pp->v[5]].tmp;
  }
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;

    pq->v[0] = mesh->point[pq->v[0]].tmp;
    pq->v[1] = mesh->point[pq->v[1]].tmp;
    pq->v[2] = mesh->point[pq->v[2]].tmp;
    pq->v[3] = mesh->point[pq->v[3]].tmp;
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure (unused).
 * \return the number of corners if success, -1 otherwise
 *
 * Pack a sparse point array and update the element vertices according to their
 * new indices.
 *
 */
MMG5_int MMG3D_pack_points(MMG5_pMesh mesh) {
  MMG5_int np, nc;

  /** Store in tmp the pack index of each point and count the corner*/
  if ( !MMG3D_mark_packedPoints(mesh,&np,&nc) ) return -1;

  /** Update the element vertices indices */
  if ( !MMG3D_update_eltsVertices(mesh) ) return -1;

  /** Pack the point array */
  if ( MMG3D_pack_pointArray(mesh)<0 ) return -1;

  return nc;
}

/**
 * \param mesh pointer towarad the mesh structure.
 *
 * Set all boundary triangles to required and add a tag to detect that they are
 * not realy required.
 *
 */
void MMG3D_unset_reqBoundaries(MMG5_pMesh mesh) {
  MMG5_pTetra pt;
  MMG5_int    k;
  int         i;

  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if ( MG_EOK(pt) &&  pt->xt ) {

      for (i=0; i<6; i++) {
        if ( mesh->xtetra[pt->xt].tag[i] & MG_NOSURF ) {
          mesh->xtetra[pt->xt].tag[i] &= ~MG_REQ;
          mesh->xtetra[pt->xt].tag[i] &= ~MG_NOSURF;
        }
      }
    }
  }
  return;
}

/**
 * \param mesh pointer to the mesh structure (unused).
 * \param sol pointer to a solution structure.
 * \param met pointer to a solution structure.
 * \return 1 if success, 0 if fail or if we are unable to build triangles.
 *
 * Pack the sparse mesh and create triangles and edges before getting
 * out of library
 *
 */
int MMG3D_packMesh(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  MMG5_int           nc,nr;

  /* Remove non wanted subdomains if needed */
  MMG3D_keep_only1Subdomain ( mesh, mesh->info.nsd );

  /* compact vertices */
  if ( !mesh->point ) {
    fprintf(stderr, "\n  ## Error: %s: points array not allocated.\n",
            __func__);
    return 0;
  }
  if ( !mesh->tetra ) {
    fprintf(stderr, "\n  ## Error: %s: tetra array not allocated.\n",
            __func__);
    return 0;
  }

  /* compact tetrahedra */
  if ( mesh->adja ) {
    if ( !MMG3D_pack_tetraAndAdja(mesh) ) return 0;
  }
  else {
    if ( !MMG3D_pack_tetra(mesh) ) return 0;
  }

  /* update prisms and quads vertex indices */
  if ( !MMG3D_pack_prismsAndQuads(mesh) ) return 0;

  /* compact solutions (metric, level-set, displacement...) */
  if ( met && met->m )
    if ( !MMG3D_pack_sol(mesh,met) ) return 0;

  if ( sol && sol->m )
    if ( !MMG3D_pack_sol(mesh,sol) ) return 0;

  /*compact vertices*/
  nc = MMG3D_pack_points(mesh);
  if ( nc<0 ) return 0;

  if ( met  && met->m  ) assert(met->np ==mesh->np);
  if ( sol && sol->m ) assert(sol->np == mesh->np);

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: prism hashing problem. Exit program.\n",
            __func__);
    return 0;
  }

  /* Remove the MG_REQ tags added by the nosurf option */
  MMG3D_unset_reqBoundaries(mesh);

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8" MMG5_PRId "   CORNERS %8" MMG5_PRId "\n",mesh->np,nc);
    fprintf(stdout,"     NUMBER OF TETRAHEDRA %8" MMG5_PRId "\n",mesh->ne);
  }

  nr = MMG3D_bdryBuild(mesh);
  if ( nr < 0 ) return 0;

  /* to could save the mesh, the adjacency have to be correct */
  if ( mesh->info.ddebug && (!MMG5_chkmsh(mesh,1,1) ) ) {
    fprintf(stderr,"\n  ##  Warning: %s: invalid mesh.\n",__func__);
    return 0;
  }

  return 1;
}

int MMG3D_mmg3dlib(MMG5_pMesh mesh,MMG5_pSol met) {
  MMG5_pSol sol=NULL; // unused
  mytime    ctim[TIMEMAX];
  char      stim[32];

  /** In debug mode, check that all structures are allocated */
  assert ( mesh );
  assert ( met );
  assert ( mesh->point );
  assert ( mesh->tetra );

  MMG5_version(mesh,"3D");

  MMG3D_Set_commonFunc();


  MMG5_warnOrientation(mesh);

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  MMG3D_Free_topoTables(mesh);

  signal(SIGABRT,MMG5_excfun);
  signal(SIGFPE,MMG5_excfun);
  signal(SIGILL,MMG5_excfun);
  signal(SIGSEGV,MMG5_excfun);
  signal(SIGTERM,MMG5_excfun);
  signal(SIGINT,MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag > -1 ) {
    fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MODE UNAVAILABLE (MMG3D_IPARAM_lag):\n"
            "            YOU MUST CALL THE MMG3D_MMG3DMOV FUNCTION TO MOVE A RIGIDBODY.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.iso || mesh->info.isosurf ) {
    fprintf(stderr,"\n  ## ERROR: LEVEL-SET DISCRETISATION UNAVAILABLE"
            " (MMG3D_IPARAM_iso or MMG3D_IARAM_isosurf ):\n"
            "          YOU MUST CALL THE MMG3D_MMG3DMOV FUNCTION TO USE THIS OPTION.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.optimLES && met->size==6 ) {
    fprintf(stdout,"\n  ## ERROR: STRONG MESH OPTIMIZATION FOR LES METHODS"
            " UNAVAILABLE (MMG3D_IPARAM_optimLES) WITH AN ANISOTROPIC METRIC.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- MMG3DLIB: INPUT DATA\n");

  chrono(ON,&(ctim[1]));

  /* check input */
  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
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

#ifdef USE_SCOTCH
  MMG5_warnScotch(mesh);
#endif

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  --  INPUT DATA COMPLETED.     %s\n",stim);

  /* analysis */
  chrono(ON,&(ctim[2]));
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  /* reset fem value to user setting (needed for multiple library call) */
  mesh->info.fem = mesh->info.setfem;

  /* scaling mesh */
  if ( !MMG5_scaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);

  MMG3D_setfunc(mesh,met);

  /* specific meshing */
  if ( mesh->info.optim ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG3D_Set_constantSize(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  if ( !MMG3D_tetraQual(mesh,met,0) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);

  if ( mesh->info.imprim > 0  ||  mesh->info.imprim < -1 ) {
    if ( !MMG3D_inqua(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  /* mesh analysis */
  if ( !MMG3D_analys(mesh) ) {
    if ( !MMG5_unscaleMesh(mesh,met,NULL) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 1 && met->m ) MMG3D_prilen(mesh,met,0);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");
  }

  /* renumerotation if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,NULL) )
  {
    if ( !MMG5_unscaleMesh(mesh,met,NULL) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

#ifdef MMG_PATTERN
  if ( !MMG5_mmg3d1_pattern(mesh,met,NULL) ) {
    if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
    if ( !MMG5_unscaleMesh(mesh,met,NULL) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }
#else
  if ( !MMG5_mmg3d1_delone(mesh,met,NULL) ) {
    if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
    if ( !MMG5_unscaleMesh(mesh,met,NULL) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }
#endif

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }

  /* last renum to give back a good numbering to the user */
  if ( !MMG5_scotchCall(mesh,met,NULL,NULL) )
  {
    if ( !MMG5_unscaleMesh(mesh,met,NULL) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  /* save file */
  if ( !MMG3D_outqua(mesh,met) ) {
    if ( !MMG5_unscaleMesh(mesh,met,NULL) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 4 && met->m )
    MMG3D_prilen(mesh,met,1);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim > 0 )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !MMG5_unscaleMesh(mesh,met,NULL) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  if ( !MMG3D_packMesh(mesh,met,sol) )       _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3D\n  %s\n\n",MG_STR,MG_STR);
  }
  _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);
}

/* Level set function discretization and remeshing mode */
int MMG3D_mmg3dls(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol umet) {
  MMG5_pSol met=NULL;
  mytime    ctim[TIMEMAX];
  char      stim[32];
  int8_t    mettofree = 0;

  /** In debug mode, check that all structures are allocated */
  assert ( mesh );
  assert ( sol );
  assert ( mesh->point );
  assert ( mesh->tetra );

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

  MMG3D_Set_commonFunc();

  signal(SIGABRT,MMG5_excfun);
  signal(SIGFPE,MMG5_excfun);
  signal(SIGILL,MMG5_excfun);
  signal(SIGSEGV,MMG5_excfun);
  signal(SIGTERM,MMG5_excfun);
  signal(SIGINT,MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.lag > -1 ) {
    fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MODE UNAVAILABLE (MMG3D_IPARAM_lag):\n"
            "            YOU MUST CALL THE MMG3D_MMG3DMOV FUNCTION TO MOVE A RIGIDBODY.\n");
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.optimLES ) {
    fprintf(stdout,"\n  ## ERROR: STRONG MESH OPTIMIZATION FOR LES METHODS"
            " UNAVAILABLE (MMG3D_IPARAM_optimLES) IN ISOSURFACE"
            " DISCRETIZATION MODE.\n");
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  /* Specific meshing */
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

#ifdef USE_SCOTCH
  MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- MMG3DLS: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  MMG5_warnOrientation(mesh);

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  MMG3D_Free_topoTables(mesh);

  if ( sol->np && (sol->np != mesh->np) ) {
    if ( mettofree ) {  MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,sol,met,MMG5_STRONGFAILURE);
  }
  else if ( sol->size!=1 ) {
    fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE.\n");
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  if ( met && met->np &&  (met->np != mesh->np) ) {
    fprintf(stdout,"\n  ## WARNING: WRONG METRIC NUMBER. IGNORED\n");
    if ( mettofree ) {  MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,sol,met,MMG5_STRONGFAILURE);
  }

  /* clean old isosurface */
  if ( !MMG3D_Clean_isoSurf(mesh) ) {
    fprintf(stderr,"\n  ## Unable to clean old isosurface.\n");
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  -- INPUT DATA COMPLETED.     %s\n",stim);

  chrono(ON,&(ctim[2]));

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 1 : ISOSURFACE DISCRETIZATION\n");
  }

  /* reset fem value to user setting (needed for multiple library call) */
  mesh->info.fem = mesh->info.setfem;

  /* scaling mesh */
  if ( !MMG5_scaleMesh(mesh,met,sol) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }

  MMG3D_setfunc(mesh,met);

  if ( !MMG3D_tetraQual(mesh,met,0) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 0 ||  mesh->info.imprim < -1 ) {
    if ( !MMG3D_inqua(mesh,met) ) {
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      if ( !MMG5_unscaleMesh(mesh,met,sol) ) {
        _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      }
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_LOWFAILURE);
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
    if ( !MMG3D_doSol(mesh,met) ) {
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      if ( !MMG5_unscaleMesh(mesh,met,sol) ) {
        _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE); }
      MMG5_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
    }
  }

  /* Discretization of the mesh->info.ls isovalue of sol in the mesh */
  if ( !MMG3D_mmg3d2(mesh,sol,met) ) {
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
    if ( !MMG3D_Set_constantSize(mesh,met) ) {
      if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
      if ( !MMG5_unscaleMesh(mesh,met,sol) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
  }

  /* mesh analysis */
  if ( !MMG3D_analys(mesh) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( !MMG5_unscaleMesh(mesh,met,sol) )  _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,sol,met,MMG5_LOWFAILURE);
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

  /* renumerotation if available */
  if ( !MMG5_scotchCall(mesh,met,NULL,NULL) )
  {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( !MMG5_unscaleMesh(mesh,met,sol) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,sol,met,MMG5_LOWFAILURE);
  }

#ifdef MMG_PATTERN
  if ( !MMG5_mmg3d1_pattern(mesh,met,NULL) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
    if ( !MMG5_unscaleMesh(mesh,met,sol) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,sol,met,MMG5_LOWFAILURE);
  }
#else
  if ( !MMG5_mmg3d1_pattern(mesh,met,NULL) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
      _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    }
    if ( !MMG5_unscaleMesh(mesh,met,sol) )    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,sol,met,MMG5_LOWFAILURE);
  }
#endif

  chrono(OFF,&(ctim[4]));
  printim(ctim[4].gdif,stim);
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
  }

  /* last renum to give back a good numbering to the user */
  if ( !MMG5_scotchCall(mesh,met,NULL,NULL) )
  {
    if ( !MMG5_unscaleMesh(mesh,met,NULL) ) _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,sol,MMG5_LOWFAILURE);
  }

  /* save file */
  if ( !MMG3D_outqua(mesh,met) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    if ( !MMG5_unscaleMesh(mesh,met,sol) )   _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,sol,met,MMG5_LOWFAILURE);
  }

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim > 0 )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !MMG5_unscaleMesh(mesh,met,sol) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  if ( !MMG3D_packMesh(mesh,sol,met) ) {
    if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
    _LIBMMG5_RETURN(mesh,met,sol,MMG5_STRONGFAILURE);
  }
  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n   MMG3DLS: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3D\n  %s\n\n",MG_STR,MG_STR);
  }
  if ( mettofree ) { MMG5_DEL_MEM(mesh,met->m);MMG5_SAFE_FREE (met); }
  _LIBMMG5_RETURN(mesh,met,sol,MMG5_SUCCESS);
}


int MMG3D_mmg3dmov(MMG5_pMesh mesh,MMG5_pSol met, MMG5_pSol disp) {
  mytime    ctim[TIMEMAX];
  char      stim[32];
  MMG5_int  *invalidTets;
  MMG5_int  k;
  MMG5_int  ier;

  /** In debug mode, check that all structures are allocated */
  assert ( mesh );
  assert ( met );
  assert ( mesh->point );
  assert ( mesh->tetra );

  MMG5_version(mesh,"3D");

  MMG3D_Set_commonFunc();

  signal(SIGABRT,MMG5_excfun);
  signal(SIGFPE,MMG5_excfun);
  signal(SIGILL,MMG5_excfun);
  signal(SIGSEGV,MMG5_excfun);
  signal(SIGTERM,MMG5_excfun);
  signal(SIGINT,MMG5_excfun);

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Check options */
  if ( mesh->info.iso || mesh->info.isosurf ) {
    fprintf(stderr,"\n  ## ERROR: LEVEL-SET DISCRETISATION UNAVAILABLE"
            " (MMG3D_IPARAM_iso || MMG3D_IPARAM_isosurf ):\n"
            "          YOU MUST CALL THE MMG3D_mmg3dls FUNCTION TO USE THIS OPTION.\n");
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }
  else if ( mesh->info.optimLES ) {
    fprintf(stdout,"\n  ## ERROR: STRONG MESH OPTIMIZATION FOR LES METHODS"
            " UNAVAILABLE (MMG3D_IPARAM_optimLES) IN LAGRANGIAN MODE.\n");
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }
  if ( mesh->info.optim ) {
    printf("\n  ## ERROR: OPTIM OPTION UNAVAILABLE IN LAGRANGIAN MODE\n");
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }
  if ( mesh->info.hsiz>0. ) {
    printf("\n  ## ERROR: HSIZ OPTION UNAVAILABLE IN LAGRANGIAN MODE\n");
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }

#ifdef USE_SCOTCH
  MMG5_warnScotch(mesh);
#endif

  if ( mesh->info.imprim > 0 ) fprintf(stdout,"\n  -- MMG3DMOV: INPUT DATA\n");
  /* load data */
  chrono(ON,&(ctim[1]));
  MMG5_warnOrientation(mesh);

  /** Free topologic tables (adja, xpoint, xtetra) resulting from a previous
   * run */
  MMG3D_Free_topoTables(mesh);

  if ( mesh->info.lag == -1 ) {
    if ( mesh->info.imprim > 0 )
      fprintf(stdout,"\n  ## Warning: displacement mode for the rigidbody"
              " movement is not set.\n"
              "               Lagrangian displacement computed according"
              " to mode 1.\n");
    mesh->info.lag = 1;
  }

#ifndef USE_ELAS
  fprintf(stderr,"\n  ## ERROR: YOU NEED TO COMPILE WITH THE USE_ELAS"
          " CMake's FLAG SET TO ON TO USE THE RIGIDBODY MOVEMENT LIBRARY.\n");
  _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
#endif

  if ( !disp ) {
    fprintf(stderr,"\n  ## ERROR: IN LAGRANGIAN MODE, A STRUCTURE OF TYPE"
            " \"MMG5_pSoL\" IS NEEDED TO STORE THE DISPLACEMENT FIELD.\n"
            "            THIS STRUCTURE MUST BE DIFFERENT FROM THE ONE USED"
            " TO STORE THE METRIC.\n");
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }
  if (disp->np && (disp->np != mesh->np) ) {
    fprintf(stdout,"\n  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,disp->m);
    disp->np = 0;
  }
  else if (disp->size!=3) {
    fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MOTION OPTION NEED A VECTORIAL DISPLACEMENT\n");
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
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

  /* reset fem value to user setting (needed for multiple library call) */
  mesh->info.fem = mesh->info.setfem;

  /* scaling mesh */
  if ( !MMG5_scaleMesh(mesh,met,disp) )   _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);

  MMG3D_setfunc(mesh,met);

  if ( !MMG3D_tetraQual(mesh,met,0) )   _LIBMMG5_RETURN(mesh,met,disp,MMG5_LOWFAILURE);

  if ( mesh->info.imprim > 0  ||  mesh->info.imprim < -1 ) {
    if ( !MMG3D_inqua(mesh,met) ) {
      _LIBMMG5_RETURN(mesh,met,disp,MMG5_LOWFAILURE);
    }
  }

  /* mesh analysis */
  if ( !MMG3D_analys(mesh) ) {
    if ( !MMG5_unscaleMesh(mesh,met,disp) )    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 4 && met->m ) {
    MMG3D_prilen(mesh,met,0);
  }

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);

  /* mesh adaptation */
  chrono(ON,&(ctim[3]));
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 2 : LAGRANGIAN MOTION\n");
  }

  /* renumerotation if available */
  if ( !MMG5_scotchCall(mesh,met,disp,NULL) )
  {
    if ( !MMG5_unscaleMesh(mesh,met,disp) )    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  /* Lagrangian mode */
  invalidTets = NULL;
  ier = MMG5_mmg3d3(mesh,disp,met,&invalidTets);
  if ( !ier ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }
  else if ( ier < 0 ) {
    printf("\n  ## Warning: Unable to perform any movement "
           "(%" MMG5_PRId " intersecting tetrahedra).\n",-ier);
    if ( mesh->info.imprim > 1 ) {
      printf("     List of invalid tets: ");
      for ( k=0; k<-ier; ++k ) {
        printf("%" MMG5_PRId " ",MMG3D_indElt(mesh,invalidTets[k]));
      }
      printf("\n\n");
    }
    MMG5_SAFE_FREE(invalidTets);
  }
  disp->npi = disp->np;

  if ( (ier > 0) && mesh->info.optim ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      if ( !MMG5_unscaleMesh(mesh,met,disp) )    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
      MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
  }

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }

  /* End with a classical remeshing stage, provided mesh->info.lag >= 1 */
  if ( (ier > 0) && (mesh->info.lag >= 1) ) {
    chrono(ON,&(ctim[4]));
    if ( mesh->info.imprim > 0 ) {
      fprintf(stdout,"\n  -- PHASE 3 : MESH IMPROVEMENT\n");
    }

    /* renumerotation if available */
    if ( !MMG5_scotchCall(mesh,met,disp,NULL) )
    {
      if ( !MMG5_unscaleMesh(mesh,met,disp) )    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
      MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }

#ifdef MMG_PATTERN
    if ( !MMG5_mmg3d1_pattern(mesh,met,NULL) ) {
      if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
        fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
        _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
      }
      if ( !MMG5_unscaleMesh(mesh,met,disp) )    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
      MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
#else
    if ( !MMG5_mmg3d1_delone(mesh,met,NULL) ) {
      if ( !(mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
        fprintf(stderr,"\n  ## Hashing problem. Invalid mesh.\n");
        _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
      }
      if ( !MMG5_unscaleMesh(mesh,met,disp) )    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
      MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
    }
#endif

    chrono(OFF,&(ctim[4]));
    printim(ctim[4].gdif,stim);
    if ( mesh->info.imprim > 0 ) {
      fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
    }
  }

  /* last renum to give back a good numbering to the user */
  if ( !MMG5_scotchCall(mesh,met,disp,NULL) )
  {
    if ( !MMG5_unscaleMesh(mesh,met,disp) ) _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
    MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  /* save file */
  if ( !MMG3D_outqua(mesh,met) ) {
    if ( !MMG5_unscaleMesh(mesh,met,disp) ) {
      disp->npi = disp->np;
      _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
    }
    MMG5_RETURN_AND_PACK(mesh,met,disp,MMG5_LOWFAILURE);
  }

  if ( mesh->info.imprim > 1 && met->m )
    MMG3D_prilen(mesh,met,1);

  chrono(ON,&(ctim[1]));
  if ( mesh->info.imprim > 0 )  fprintf(stdout,"\n  -- MESH PACKED UP\n");
  if ( !MMG5_unscaleMesh(mesh,met,disp) ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }
  if ( !MMG3D_packMesh(mesh,met,disp) ) {
    disp->npi = disp->np;
    _LIBMMG5_RETURN(mesh,met,disp,MMG5_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[1]));

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( mesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n   MMG3DMOV: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE MMG3D\n  %s\n\n",MG_STR,MG_STR);
  }
  disp->npi = disp->np;
  _LIBMMG5_RETURN(mesh,met,disp,MMG5_SUCCESS);
}

/**
* Set common pointer functions between mmgs and mmg3d to the matching mmg3d
* functions.
*/
void MMG3D_Set_commonFunc(void) {
    MMG5_bezierCP = MMG5_mmg3dBezierCP;
    MMG5_chkmsh = MMG5_mmg3dChkmsh;
    MMG5_indPt = MMG3D_indPt;
    MMG5_indElt = MMG3D_indElt;
    MMG5_grad2met_ani = MMG5_grad2metSurf;
    MMG5_grad2metreq_ani = MMG5_grad2metSurfreq;
    MMG5_solTruncature_ani = MMG5_3dSolTruncature_ani;

#ifdef USE_SCOTCH
    MMG5_renumbering = MMG5_mmg3dRenumbering;
#endif
}
