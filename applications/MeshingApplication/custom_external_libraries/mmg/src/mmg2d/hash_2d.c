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
#include "mmg2d.h"

#define KTA     7
#define KTB    11

int MMG2D_hashNew(HashTable *hash,int hsize,int hmax) {
  int   k;

  hash->size  = hsize;
  hash->nxtmax =hmax+1;
  hash->hnxt  = hsize;
  MMG5_SAFE_CALLOC(hash->item,hash->nxtmax,Hedge,return 0);

  for (k=hash->size; k<hash->nxtmax; k++)
    hash->item[k].nxt = k+1;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \return 1 if success, 0 if fail
 *
 * Create adjacency relations between the triangles in the mesh
 *
 */
int MMG2D_hashTria(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1;
  int      *hcode,*link,inival,hsize,iadr;
  unsigned char   i,ii,i1,i2;
  unsigned int    key;

  if ( mesh->adja )  return 1;
  if ( !mesh->nt )  return 0;

  /* memory alloc */
  MMG5_SAFE_CALLOC(hcode,mesh->nt+1,int,return 0);

  /* memory alloc */
  MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                return 0;);
  MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int,return 0);

  link  = mesh->adja;
  hsize = mesh->nt;

  /* init */
  inival = INT_MAX;
  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !pt->v[0] )  continue;
    for (i=0; i<3; i++) {
      i1 = MMG2D_idir[i+1];
      i2 = MMG2D_idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }

      /* compute key */
      key = KTA*mins + KTB*maxs;
      key = key % hsize + 1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = MMG2D_idir[i+1];
    i2 = MMG2D_idir[i+2];
    pt = &mesh->tria[k];

    mins = M_MIN(pt->v[i1],pt->v[i2]);
    maxs = M_MAX(pt->v[i1],pt->v[i2]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;

    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = MMG2D_idir[ii+1];
      i2 = MMG2D_idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }

      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  MMG5_SAFE_FREE(hcode);

  return 1;
}

/*hash edge :
  return 1 if edge exist in the table*/
int MMG2D_hashEdge(pHashTable edgeTable,int iel,int ia, int ib) {
  int         key,mins,maxs;
  Hedge      *ha;
  static char mmgErr = 0;

  /* compute key */
  if ( ia < ib ) {
    mins = ia;
    maxs = ib;
  }
  else {
    mins = ib;
    maxs = ia;
  }

  key = KTA*mins + KTB*maxs;
  key = key % edgeTable->size;
  ha  = &edgeTable->item[key];
  if ( ha->min ) {
    /* edge exist*/
    if ( ha->min == mins && ha->max == maxs ) {
      return ha->iel;
    }
    else {
      while ( ha->nxt && ha->nxt < edgeTable->nxtmax ) {
        ha = &edgeTable->item[ha->nxt];
        if ( ha->min == mins && ha->max == maxs )
          return ha->iel;
      }
      ha->nxt = edgeTable->hnxt;
      ha      = &edgeTable->item[edgeTable->hnxt];
      ++edgeTable->hnxt;
      if ( edgeTable->hnxt == edgeTable->nxtmax ) {
        if ( !mmgErr ) {
          mmgErr = 1;
          fprintf(stderr,"\n  ## Error: %s: memory alloc problem (edge): %d.\n",
                  __func__,edgeTable->nxtmax);
        }
        return 0;
      }
    }
  }
  /* insert */
  ha->min = mins;
  ha->max = maxs;
  ha->iel = iel;
  ha->nxt = 0;

  return 0;

}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 0 if fail, 1 otherwise
 *
 * Transfer some input edge data to the corresponding triangles fields
 *
 */
int MMG2D_assignEdge(MMG5_pMesh mesh) {
  MMG5_Hash      hash;
  MMG5_pTria      pt;
  MMG5_pEdge      pa;
  int             k,ia;
  char            i,i1,i2;

  if ( !mesh->na ) return 1;

  /* Temporarily allocate a hash structure for storing edges */
  hash.siz = mesh->na;
  hash.max = 3*mesh->na+1;

  MMG5_ADD_MEM(mesh,(hash.max+1)*sizeof(MMG5_hedge),"hash table",return 0);
  MMG5_SAFE_CALLOC(hash.item,hash.max+1,MMG5_hedge,return 0);

  hash.nxt = mesh->na;

  for (k=mesh->na; k<hash.max; k++)
    hash.item[k].nxt = k+1;

  /* hash mesh edges */
  for (k=1; k<=mesh->na; k++)
    MMG5_hashEdge(mesh,&hash,mesh->edge[k].a,mesh->edge[k].b,k);

  /* set references to triangles */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      ia = MMG5_hashGet(&hash,pt->v[i],pt->v[i1]);
      if ( ia ) {
        i2 = MMG5_inxt2[i1];
        pa = &mesh->edge[ia];
        pt->edg[i2] = pa->ref;
        pt->tag[i2] |= pa->tag;
      }
    }
  }

  /* Delete the hash for edges */
  MMG5_DEL_MEM(mesh,hash.item);
  MMG5_DEL_MEM(mesh,mesh->edge);
  mesh->na = 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 1 if success, 0 if fail
 *
 * Create the edges in the mesh from the information stored in the triangles, or
 * by identifying the different components of the mesh.
 *
 * \remark Possible extension needed to take into account constrained edges
 *
 */
int MMG2D_bdryEdge(MMG5_pMesh mesh) {
  MMG5_pTria      pt,pt1;
  MMG5_pEdge      pa;
  MMG5_pPoint     p0;
  int             k,*adja,natmp,iel;
  char            i,i1,i2;

  natmp = 0;
  mesh->na = 0;

  /* First step: Count number of boundary edges */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      iel = adja[i] / 3;
      pt1 = &mesh->tria[iel];

      if ( iel && pt->ref <= pt1->ref ) continue;
      natmp++;
    }
  }

  /* Second step: Create edge mesh and store the corresponding edges */
  MMG5_ADD_MEM(mesh,(natmp+1)*sizeof(MMG5_Edge),"edges",return 0);
  MMG5_SAFE_CALLOC(mesh->edge,natmp+1,MMG5_Edge,return 0);

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      iel = adja[i] / 3;
      pt1 = &mesh->tria[iel];

      if ( iel && pt->ref <= pt1->ref ) continue;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];

      mesh->na++;
      pa = &mesh->edge[mesh->na];
      pa->a = pt->v[i1];
      pa->b = pt->v[i2];

      pa->tag = pt->tag[i];
      pa->ref = pt->edg[i];

    }
  }

  /* Set point tags */
  for (k=1; k<=mesh->na; k++) {
    pa = &mesh->edge[k];
    p0 = &mesh->point[pa->a];
    p0->tag |= MG_BDY;

    p0 = &mesh->point[pa->b];
    p0->tag |= MG_BDY;
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \return 0 if memory problem (uncomplete mesh), 1 otherwise.
 *
 * Pack the mesh and metric and create explicitly all the mesh structures
 * (edges).
 *
 */
int MMG2D_pack(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria         pt,ptnew,pt1;
  MMG5_pEdge         ped;
  MMG5_pPoint        ppt,pptnew;
  int                np,ned,nt,k,iel,nbl,isol,isolnew,memWarn,nc;
  int                iadr,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;
  char               i,i1,i2;

  /* Recreate adjacencies if need be */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Warning: %s: hashing problem. Exit program.\n",
            __func__);
    return 0;
  }

  /* Pack vertex indices */
  np = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;
    ppt->tmp = ++np;

    if ( ppt->tag & MG_CRN )  nc++;

    if ( ppt->tag & MG_NOSURF ) {
      ppt->tag &= ~MG_NOSURF;
      ppt->tag &= ~MG_REQ;
    }
  }

  /* Count the number of edges in the mesh */
  memWarn = 0;
  ned = 0;

  if ( mesh->edge ) {
    fprintf(stderr,"\n  ## Warning: %s: unexpected edge table..."
            " Ignored data.\n",__func__);
    MMG5_DEL_MEM(mesh,mesh->edge);
    mesh->na = 0;
  }

  mesh->na = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      iel = adja[i] / 3;

      if ( pt->tag[i] & MG_NOSURF ) {
        pt->tag[i] &= ~MG_REQ;
        pt->tag[i] &= ~MG_NOSURF;
      }

      pt1 = &mesh->tria[iel];
      if ( (!iel) || (pt->ref > pt1->ref) ) {
        ++mesh->na;
      }
      else if ( (pt->ref==pt1->ref) && MG_SIN(pt->tag[i]) ) {
        ++mesh->na;
      }
    }
  }

  /* Pack edges */
  mesh->namax = mesh->na+1;
  if ( mesh->na ) {

    MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"final edges", memWarn=1);

    if ( memWarn ) {
      if ( mesh->info.ddebug )
        printf("  -- Attempt to allocate a smallest edge table...\n");
      mesh->namax = mesh->na;
      memWarn = 0;
      MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"final edges",
                    fprintf(stderr,"\n  ## Warning: %s: uncomplete mesh.\n",
                            __func__);
                    memWarn=1);
    }

    if ( memWarn )
      mesh->na = 0;
    else {
      /* We have enough memory to allocate the edge table */
      MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge, return 0);

      nt = 0;
      for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MG_EOK(pt) ) continue;
        ++nt;
        adja = &mesh->adja[3*(k-1)+1];

        for (i=0; i<3; i++) {
          i1 = MMG5_inxt2[i];
          i2 = MMG5_iprv2[i];
          iel = adja[i] / 3;
          pt1 = &mesh->tria[iel];
          if ( !iel || (pt->ref > pt1->ref) ||
               ((pt->ref==pt1->ref) && MG_SIN(pt->tag[i])) ) {
            ++ned;
            ped = &mesh->edge[ned];
            ped->a = pt->v[i1];
            ped->b = pt->v[i2];
            /* the base field is used to be able to recover the tria (and its face)
             * from which comes a boundary edge */
            ped->base = 3*nt+i;
            ped->ref = pt->edg[i];
            ped->tag = pt->tag[i];
          }
        }
      }
    }
  }

  for (k=1; k<=mesh->na; k++) {
    ped  = &mesh->edge[k];
    if ( !ped->a ) continue;
    ped->a = mesh->point[ped->a].tmp;
    ped->b = mesh->point[ped->b].tmp;
  }

  /* Pack triangles */
  nt  = 0;
  nbl = 1;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    nt++;

    if ( k != nbl ) {
      ptnew = &mesh->tria[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tria));

      /* Update the adjacency */
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      iadrnew = 3*(nbl-1) + 1;
      adjanew = &mesh->adja[iadrnew];

      for(i=0; i<3; i++) {
        adjanew[i] = adja[i];
        if ( !adja[i] ) continue;
        iadrv = 3*(adja[i]/3-1)+1;
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

  /* Pack metric map */
  if ( sol && sol->m ) {
    nbl = 1;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      isol    = k * sol->size;
      isolnew = nbl * sol->size;

      for (i=0; i<sol->size; i++)
        sol->m[isolnew + i] = sol->m[isol + i];
      ++nbl;
    }
  }

  /* Pack vertices*/
  np  = 0;
  nbl = 1;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;

    if ( k != nbl ) {
      pptnew = &mesh->point[nbl];
      memcpy(pptnew,ppt,sizeof(MMG5_Point));
      ppt->tag   = 0;
      assert ( ppt->tmp == nbl );
    }
    np++;
    if ( k != nbl ) {
      ppt = &mesh->point[k];
      memset(ppt,0,sizeof(MMG5_Point));
      ppt->tag    = 0;
    }
    nbl++;
  }
  mesh->np = np;
  if ( sol && sol->m ) sol->np  = np;

  /* Reset ppt->tmp field */
  for(k=1 ; k<=mesh->np ; k++)
    mesh->point[k].tmp = 0;

  if(mesh->np < mesh->npmax - 3) {
    mesh->npnil = mesh->np + 1;
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;
  }
  else {
    mesh->npnil = 0;
  }

  /* Reset garbage collector */
  if ( mesh->nt < mesh->ntmax - 3 ) {
    mesh->nenil = mesh->nt + 1;
    for (k=mesh->nenil; k<mesh->ntmax-1; k++)
      mesh->tria[k].v[2] = k+1;
  }
  else {
    mesh->nenil = 0;
  }

  if ( mesh->info.imprim > 0 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8d\n",mesh->nt);

    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8d\n",mesh->na);
  }

  if ( memWarn ) return 0;

  return 1;
}
