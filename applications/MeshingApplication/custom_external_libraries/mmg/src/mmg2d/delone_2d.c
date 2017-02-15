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

#define  _MMG2_EPSRAD       1.00005
#define  _MMG2_AREAMIN       1e-30

#define KTA     7
#define KTB    11

/* cavity correction for quality */
static int
_MMG2_correction_iso(MMG5_pMesh mesh,int ip,int *list,int ilist,int nedep) {
  MMG5_pPoint     ppt,p1,p2;
  MMG5_pTria      pt;
  double           dd,ux,uy,vx,vy;
  int             *adja,i,ipil,iel,lon,iadr,adj,ib,ic,base,ncor;
  int              vois[3];

  ppt  = &mesh->point[ip];
  if ( !M_VOK(ppt) )  return(ilist);
  base = mesh->base;
  lon  = ilist;
  do {
    ipil = lon-1;
    ncor = 0;

    while ( ipil >= 0 ) {
      iel  = list[ipil];
      iadr = (iel-1)*3 + 1;
      adja = &mesh->adja[iadr];
      vois[0]  = adja[0] /3;
      vois[1]  = adja[1] /3;
      vois[2]  = adja[2] /3;
      pt   = &mesh->tria[iel];
      for (i=0; i<3; i++) {
        adj = vois[i];
        if ( adj && mesh->tria[adj].base == base )  continue;

        ib = pt->v[ MMG2_iare[i][0] ];
        ic = pt->v[ MMG2_iare[i][1] ];

        p1 = &mesh->point[ib];
        p2 = &mesh->point[ic];

        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];

        vx = ppt->c[0] - p1->c[0];
        vy = ppt->c[1] - p1->c[1];

        /* area PBC */
        dd =  ux*vy - uy*vx;
        if ( dd < _MMG2_AREAMIN )  break;

      }
      if ( i < 3 /*||  pt->tag & MG_REQ*/ ) {
        if ( ipil < nedep )
        {/*printf("on veut tout retirer ? %d %d -- %d\n",ipil,nedep,iel);*/return(0);   }
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
  return(lon);
}

/* hash mesh edge v[0],v[1] (face i of iel) */
int _MMG2_hashEdgeDelone(MMG5_pMesh mesh,HashTable *hash,int iel,int i,int *v) {
  int             *adja,iadr,jel,j,key,mins,maxs;
  Hedge     *ha;

  /* compute key */
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
      return(1);
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
          iadr = (jel-1)*3 + 1;
          adja = &mesh->adja[iadr];
          adja[j] = iel*3 + i;
          return(1);
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
      if(mesh->info.imprim > 6) fprintf(stdout," ## Warning: overflow\n");
      return(0);
    }
    return(1);
  }

  /* insert */
  ha->min = mins;
  ha->max = maxs;
  ha->iel = iel*3 + i;
  ha->nxt = 0;

  return(1);
}

/** Return a negative value for ilist if one of the tet of the cavity is required */
int _MMG2_cavity(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int *list) {
  MMG5_pPoint     ppt;
  MMG5_pTria      pt,pt1,ptc;
  double     c[2],crit,dd,eps,ray,ct[6];
  int        *adja,*adjb,adj,adi,voy,i,j,ilist,ipil,jel,iadr,base;
  int         vois[3],l;
  int         tref;//,isreq;

  ppt = &mesh->point[ip];
  base  = ++mesh->base;

  //isreq = 0;

  tref = list[0];
  mesh->tria[list[0]].base = base;

  /* grow cavity by adjacency */
  eps   = _MMG2_EPSRAD*_MMG2_EPSRAD;
  ilist = 1;
  ipil  = 0;

  do {
    jel  = list[ipil];
    iadr = (jel-1)*3 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0];
    vois[1]  = adja[1];
    vois[2]  = adja[2];
    ptc  = &mesh->tria[jel];

    for (i=0; i<3; i++) {
      adj = vois[i] /3;
      voy = vois[i] % 3;
      if ( !adj )  continue;
      pt  = &mesh->tria[adj];
      /* boundary face */

      if ( pt->base == base || pt->ref != ptc->ref )  continue;
      for (j=0,l=0; j<3; j++,l+=2) {
        memcpy(&ct[l],mesh->point[pt->v[j]].c,2*sizeof(double));
      }

      if ( !_MMG2_cenrad_iso(mesh,ct,c,&ray) )  continue;
      crit = eps * ray;

      /* Delaunay criterion */
      dd = (ppt->c[0] - c[0]) * (ppt->c[0] - c[0]) \
        + (ppt->c[1] - c[1]) * (ppt->c[1] - c[1]);
      if ( dd > crit )  continue;

      /* lost face(s) */
      iadr = (adj-1)*3 + 1;
      adjb = &mesh->adja[iadr];

      for (j=0; j<3; j++) {
        if ( j == voy )  continue;
        adi = adjb[j] /3;
        if ( !adi )  continue;
        pt1 = &mesh->tria[adi];
        if ( pt1->base == base && adi != jel ) {
          if ( !adi || pt1->ref != tref )  break;
        }
      }
      /* store tria */
      if ( j == 3 ) {
        //if ( pt->tag & M_REQUIRED ) isreq = 1;
        pt->base = base;
        list[ilist++] = adj;
      }
    }
    if ( ilist > MMG2_LONMAX - 3 ) return(-1);

    ++ipil;
  }
  while ( ipil < ilist );

  ilist = _MMG2_correction_iso(mesh,ip,list,ilist,1);

  //if ( isreq ) ilist = -fabs(ilist);

  return(ilist);
}
/* cavity -> ball */
int _MMG2_delone(MMG5_pMesh mesh,MMG5_pSol sol,int ip,int *list,int ilist) {
  MMG5_pPoint     ppt;
  MMG5_pTria      pt,pt1;
  int        *adja,*adjb,i,j,k,iel,jel,old,v[2],iadr,base,size;
  int         vois[3],iadrold;
  short       i1;
  char        alert;
  int         tref,ielnum[3*MMG2_LONMAX+1];
  HashTable   hedg;

  base = mesh->base;
  /* external faces */
  size = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tria[old];
    iadr = (old-1)*3 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0]/3 ;
    vois[1]  = adja[1]/3 ;
    vois[2]  = adja[2]/3 ;
    for (i=0; i<3; i++) {
      jel = vois[i];
      if ( !jel || mesh->tria[jel].base != base ) {
        for (j=0; j<2; j++) {
          i1  = MMG2_iare[i][j];
          ppt = &mesh->point[ pt1->v[i1] ];
          ppt->tagdel |= M_MOVE;
        }
        size++;
      }
    }
  }
  /* check isolated vertex */
  alert = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tria[old];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      if ( !(ppt->tagdel & M_MOVE) )  alert = 1;
    }
  }
  /* reset tag */
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tria[old];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      ppt->tagdel &= ~M_MOVE;
    }
  }
  if ( alert )  {return(0);}
  /* hash table params */
  if ( size > 3*MMG2_LONMAX )  return(0);
  if ( !MMG2_hashNew(&hedg,size,3*size) ) { /*3*size suffit */
    fprintf(stdout,"  ## Unable to complete mesh.\n");
    return(-1);
  }

  /*tria allocation : we create "size" tria*/
  ielnum[0] = size;
  for (k=1 ; k<=size ; k++) {
    ielnum[k] = _MMG2D_newElt(mesh);
    if(!ielnum[k]) {
      _MMG5_TRIA_REALLOC(mesh,ielnum[k],mesh->gap,
                         printf("  ## Error: unable to allocate a new element.\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      pt1  = &mesh->tria[old];
    }
  }

  size = 1;
  for (k=0; k<ilist; k++) {
    old  = list[k];

    iadrold = (old-1)*3 + 1;
    adja = &mesh->adja[iadrold];
    vois[0]  = adja[0];
    vois[1]  = adja[1];
    vois[2]  = adja[2];

    pt   = &mesh->tria[old];

    for (i=0; i<3; i++) {
      jel = vois[i] /3;
      j   = vois[i] % 3;

      /* external face */
      if ( !jel || (mesh->tria[jel].base != base) ) {
        iel = ielnum[size++];
        assert(iel);

        pt1 = &mesh->tria[iel];
        memcpy(pt1,pt,sizeof(MMG5_Tria));
        pt1->v[i] = ip;
        pt1->qual = MMG2_caltri_in(mesh,sol,pt1);
        pt1->ref = mesh->tria[old].ref;

        if ( pt1->qual < 1e-10 ) {
          fprintf(stdout,"  ## Warning: creation of a very bad element.\n");
        }

        iadr = (iel-1)*3 + 1;
        adjb = &mesh->adja[iadr];
        adjb[i] = adja[i];


        if ( jel ) {
          iadr = (jel-1)*3 + 1;
          adjb = &mesh->adja[iadr];
          adjb[j] = iel*3 + i;
        }
        /* internal edges (p1,p2) */
        for (j=0; j<3; j++) {
          if ( j != i ) {
            v[0] = pt1->v[ MMG2_iare[j][0] ];
            v[1] = pt1->v[ MMG2_iare[j][1] ];

            _MMG2_hashEdgeDelone(mesh,&hedg,iel,j,v);
          }
        }
      }
    }
  }

  /* remove old tria */
  tref = mesh->tria[list[0]].ref;
  for (k=0; k<ilist; k++) {
    if(tref!=mesh->tria[list[k]].ref) {
      fprintf(stdout,"  ## Warning: sud-domain ignored\n");
    }
    _MMG2D_delElt(mesh,list[k]);
  }

  //ppt = &mesh->point[ip];
  //  ppt->flag = mesh->flag;
  _MMG5_SAFE_FREE(hedg.item);
  return(1);
}
