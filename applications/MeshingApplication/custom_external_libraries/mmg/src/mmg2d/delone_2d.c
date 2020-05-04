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

#define  MMG2D_AREAMIN      1e-15 //1e-20 failed : creation of too bad element

#define KTA     7
#define KTB    11

/* Cavity correction for quality */
static int MMG2D_correction_iso(MMG5_pMesh mesh,int ip,int *list,int ilist,int nedep) {
  MMG5_pTria      pt;
  MMG5_pPoint     ppt,p1,p2;
  double           dd,ux,uy,vx,vy;
  int             *adja,i,ipil,iel,lon,iadr,adj,ib,ic,base,ncor,nei[3];

  ppt  = &mesh->point[ip];
  if ( !MG_VOK(ppt) )  return ilist;
  base = mesh->base;
  lon  = ilist;
  do {
    ipil = lon-1;
    ncor = 0;

    while ( ipil >= 0 ) {
      iel  = list[ipil];
      iadr = (iel-1)*3 + 1;
      adja = &mesh->adja[iadr];
      nei[0]  = adja[0] /3;
      nei[1]  = adja[1] /3;
      nei[2]  = adja[2] /3;
      pt   = &mesh->tria[iel];

      for (i=0; i<3; i++) {
        adj = nei[i];
        /* Consider only the external faces of the cavity */
        if ( adj && mesh->tria[adj].base == base )  continue;

        ib = pt->v[ MMG5_inxt2[i] ];
        ic = pt->v[ MMG5_iprv2[i] ];

        p1 = &mesh->point[ib];
        p2 = &mesh->point[ic];

        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];

        vx = ppt->c[0] - p1->c[0];
        vy = ppt->c[1] - p1->c[1];

        /* area PBC */
        dd =  ux*vy - uy*vx;
        if ( dd < MMG2D_AREAMIN )  break;

      }

      /* Remove triangle iel from the cavity if it leads to a degenerate triangle after insertion of ppt */
      if ( i < 3 /*||  pt->tag & MG_REQ*/ ) {
        /* remove iel from list */
        pt->base = base-1;
        list[ipil] = list[--lon];
        ncor = 1;
        break;
      }
      else
        ipil--;
    }
  }
  while ( ncor > 0 && lon >= nedep );
  return lon;
}

/* Hashing routine for maintaining adjacencies during Delaunization; hash mesh edge v[0],v[1] (face i of iel) */
int MMG2D_hashEdgeDelone(MMG5_pMesh mesh,HashTable *hash,int iel,int i,int *v) {
  int             *adja,iadr,jel,j,key,mins,maxs;
  Hedge           *ha;
  static char     mmgWarn0=0;

  /* Compute key */
  if ( v[0] < v[1] ) {
    mins = v[0];
    maxs = v[1];
  }
  else {
    mins = v[1];
    maxs = v[0];
  }

  key = KTA*mins + KTB*maxs;
  key = key % hash->size;
  ha  = &hash->item[key];

  if ( ha->min ) {
    /* identical face */
    if ( ha->min == mins && ha->max == maxs ) {
      iadr = (iel-1)*3 + 1;
      adja = &mesh->adja[iadr];
      adja[i] = ha->iel;

      jel  = ha->iel /3;
      j    = ha->iel % 3;
      iadr = (jel-1)*3 + 1;
      adja = &mesh->adja[iadr];
      adja[j] = iel*3 + i;
      return 1;
    }
    else {
      while ( ha->nxt && ha->nxt < hash->nxtmax ) {
        ha = &hash->item[ha->nxt];
        if ( ha->min == mins && ha->max == maxs ) {
          iadr = (iel-1)*3 + 1;
          adja = &mesh->adja[iadr];
          adja[i] = ha->iel;

          jel  = ha->iel /3;
          j    = ha->iel % 3;
          iadr = 3*(jel-1) + 1;
          adja = &mesh->adja[iadr];
          adja[j] = 3*iel+i;
          return 1;
        }
      }
    }
    ha->nxt   = hash->hnxt;
    ha        = &hash->item[hash->hnxt];
    ha->min     = mins;
    ha->max     = maxs;
    ha->iel     = iel*3 + i;
    hash->hnxt = ha->nxt;
    ha->nxt   = 0;

    if ( hash->hnxt >= hash->nxtmax ) {
      if(mesh->info.imprim > 6) {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  ## Warning: %s: overflow.\n",__func__);
        }
      }
      return 0;
    }
    return 1;
  }

  /* If ha->man does not exist, insert it in the hash table */
  ha->min = mins;
  ha->max = maxs;
  ha->iel = 3*iel+i;
  ha->nxt = 0;

  return 1;
}

/**  Create the cavity point ip, starting from triangle list[0];
     Return a negative value for ilist if one of the triangles of the cavity is required */
int MMG2D_cavity(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int *list) {
  MMG5_pTria      pt,pt1,ptc;
  MMG5_pPoint     ppt;
  double          c[2],crit,dd,eps,rad,ct[6];
  int             *adja,*adjb,adj,adi,voy,i,j,ilist,ipil,jel,iadr,base,nei[3],l,tref; //isreq;
  static char     mmgWarn0=0;

  ppt = &mesh->point[ip];
  base  = ++mesh->base;
  //isreq = 0;
  tref = mesh->tria[list[0]].ref;
  mesh->tria[list[0]].base = base;

  /* Pile up cavity by adjacency */
  eps   = 1. + MMG5_EPSOK;
  ilist = 1;
  ipil  = 0;

  do {
    jel  = list[ipil];
    iadr = (jel-1)*3 + 1;
    adja = &mesh->adja[iadr];
    nei[0]  = adja[0];
    nei[1]  = adja[1];
    nei[2]  = adja[2];
    ptc  = &mesh->tria[jel];

    for (i=0; i<3; i++) {
      adj = nei[i] /3;
      voy = nei[i] % 3;
      if ( !adj )  continue;
      pt  = &mesh->tria[adj];

      /* Case where the triangle has already been piled, or a boundary face is hit */
      if ( pt->base == base || pt->ref != ptc->ref )  continue;

      /* Store the 6 coordinates of the vertices of pt */
      for (j=0,l=0; j<3; j++,l+=2) {
        memcpy(&ct[l],mesh->point[pt->v[j]].c,2*sizeof(double));
      }

      if ( !MMG2D_cenrad_iso(mesh,ct,c,&rad) )  continue;
      crit = eps * rad;

      /* Delaunay criterion */
      dd = (ppt->c[0] - c[0]) * (ppt->c[0] - c[0]) + (ppt->c[1] - c[1]) * (ppt->c[1] - c[1]);
      if ( dd > crit )  continue;

      /* Lost face(s); I don't understand this test: the algorithm is supposed to stop when changing references */
      iadr = (adj-1)*3 + 1;
      adjb = &mesh->adja[iadr];

      for (j=0; j<3; j++) {
        if ( j == voy )  continue;
        adi = adjb[j] /3;
        if ( !adi )  continue;
        pt1 = &mesh->tria[adi];
        if ( pt1->base == base && adi != jel ) {
          if ( pt1->ref != tref ) {
            break;
          }
        }
      }
      /* store tria */
      if ( j == 3 ) {
        //if ( pt->tag & MG_REQ ) isreq = 1;
        pt->base = base;
        list[ilist++] = adj;
      }
      else {
        if ( !mmgWarn0 ) {
          mmgWarn0 = 1;
          fprintf(stderr,"\n  ## Error: %s: we pass here at least one time but one"
                  " should never go through here.\n",__func__);
        }
      }
    }
    if ( ilist > MMG2D_LONMAX - 3 ) return -1;

    ++ipil;
  }
  while ( ipil < ilist );

  ilist = MMG2D_correction_iso(mesh,ip,list,ilist,1);
  //if ( isreq ) ilist = -fabs(ilist);
  return ilist;
}

/* Insertion in point ip in the cavity described by list */
int MMG2D_delone(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int *list,int ilist) {
  MMG5_pTria      pt,pt1;
  MMG5_pPoint     ppt;
  int             *adja,*adjb,i,j,k,iel,jel,old,v[2],iadr,base,size,nei[3],iadrold;
  int             tref,ielnum[3*MMG2D_LONMAX+1];
  short           i1;
  char            alert;
  HashTable       hedg;
  static char     mmgWarn0=0,mmgWarn1=0;

  /* Reset tagdel field */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].tagdel = 0;

  /* Triangles in the cavity are those s.t. pt->base == base */
  base = mesh->base;
  /* Count the number of external faces in the cavity, and tag the corresponding vertices */
  size = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tria[old];
    iadr = (old-1)*3 + 1;
    adja = &mesh->adja[iadr];
    nei[0]  = adja[0]/3 ;
    nei[1]  = adja[1]/3 ;
    nei[2]  = adja[2]/3 ;
    for (i=0; i<3; i++) {
      jel = nei[i];
      if ( (!jel) || (mesh->tria[jel].base != base) ) {
        for (j=0; j<2; j++) {
          i1  = MMG2D_iare[i][j];
          ppt = &mesh->point[ pt1->v[i1] ];
          ppt->tagdel = 1;
        }
        size++;
      }
    }
  }

  /* Check for an isolated vertex (the cavity should ne contain any internal vertex) */
  alert = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tria[old];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      if ( !ppt->tagdel ) {
        alert = 1;
      }
    }
  }
  /* Reset tagdel field */
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tria[old];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      ppt->tagdel = 0;
    }
  }
  if ( alert )  return 0;

  /* Hash table parameters */
  if ( size >= 3*MMG2D_LONMAX )  return 0;
  if ( !MMG2D_hashNew(&hedg,size,3*size) ) { /*3*size is enough */
    fprintf(stderr,"\n  ## Warning: %s: unable to complete mesh.\n",__func__);
    return -1;
  }

  /* Allocate memory for "size" new triangles */
  ielnum[0] = size;
  for (k=1; k<=size; k++) {
    ielnum[k] = MMG2D_newElt(mesh);
    if ( !ielnum[k] ) {
      MMG2D_TRIA_REALLOC(mesh,ielnum[k],mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                 " a new element.\n",__func__);
                          MMG5_INCREASE_MEM_MESSAGE();
                          printf("  Exit program.\n");return -1);
    }
  }

  size = 1;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt   = &mesh->tria[old];

    iadrold = (old-1)*3 + 1;
    adja = &mesh->adja[iadrold];
    nei[0]  = adja[0];
    nei[1]  = adja[1];
    nei[2]  = adja[2];

    for (i=0; i<3; i++) {
      jel = nei[i] / 3;
      j   = nei[i] % 3;

      /* Catch the associated external face */
      if ( (!jel) || (mesh->tria[jel].base != base) ) {
        assert ( size <= ielnum[0] );
        iel = ielnum[size++];
        assert(iel);

        pt1 = &mesh->tria[iel];
        memcpy(pt1,pt,sizeof(MMG5_Tria));
        pt1->v[i] = ip;
        pt1->qual = MMG2D_caltri_iso(mesh,sol,pt1);
        pt1->ref = pt->ref;

        if ( (!mmgWarn0) && (pt1->qual < MMG2D_AREAMIN) ) {
          mmgWarn0 = 1;
          fprintf(stderr,"  ## Warning: %s: creation of a very bad element.\n",
                  __func__);
        }

        /* Update adjacency via the external face */
        iadr = (iel-1)*3 + 1;
        adjb = &mesh->adja[iadr];
        adjb[i] = adja[i];

        if ( jel ) {
          iadr = (jel-1)*3 + 1;
          adjb = &mesh->adja[iadr];
          adjb[j] = iel*3 + i;
        }
        /* Update adjacency via the internal faces */
        for (j=0; j<3; j++) {
          if ( j != i ) {
            v[0] = pt1->v[ MMG5_inxt2[j] ];
            v[1] = pt1->v[ MMG5_iprv2[j] ];
            MMG2D_hashEdgeDelone(mesh,&hedg,iel,j,v);
          }
        }
      }
    }
  }

  /* Remove the old triangles */
  tref = mesh->tria[list[0]].ref;
  for (k=0; k<ilist; k++) {
    if ( (!mmgWarn1) && (tref != mesh->tria[list[k]].ref) ) {
      mmgWarn1 = 1;
      fprintf(stderr,"\n  ## Warning: %s: sud-domain ignored.\n",__func__);
    }
    MMG2D_delElt(mesh,list[k]);
  }

  //ppt = &mesh->point[ip];
  //  ppt->flag = mesh->flag;
  MMG5_SAFE_FREE(hedg.item);
  return 1;
}
