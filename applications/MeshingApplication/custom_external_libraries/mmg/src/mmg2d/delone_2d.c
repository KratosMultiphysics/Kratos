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
#include "libmmg2d_private.h"

#define  MMG2D_AREAMIN      1e-15 //1e-20 failed : creation of too bad element

#define KTA     7
#define KTB    11

/* Cavity correction for quality */
static int MMG2D_correction_iso(MMG5_pMesh mesh,MMG5_int ip,MMG5_int *list,int ilist,int nedep) {
  MMG5_pTria      pt;
  MMG5_pPoint     ppt,p1,p2;
  double          dd,ux,uy,vx,vy;
  MMG5_int        *adja,iel,iadr,adj,ib,ic,ncor,nei[3],base;
  int             lon,i,ipil;

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

/**
 * \param mesh pointer toward the mesh structure
 * \param hash pointer toward the hash table structure
 * \param iel index of triangle
 * \param i index of the face of the element
 *
 * \return 0 if fail, 1 if success
 *
 * Update of the adjacency relationships adjacencies: hash mesh edge
 * \f$v[0],v[1]\f$ (face \a i of \a iel) and fill the adjacency arrays if the
 * face has already been seen.
 *
 */
static int MMG2D_hashEdgeDelone(MMG5_pMesh mesh,MMG5_Hash *hash,MMG5_int iel,int i) {
  MMG5_pTria      pt;
  MMG5_int        *adja,iadr,jel,ip1,ip2;
  int             j;
  int16_t         i1,i2;

  pt  = &mesh->tria[iel];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  ip1 = pt->v[i1];
  ip2 = pt->v[i2];

  /* Search if the face is already hashed */
  jel = MMG5_hashGet(hash,ip1,ip2);

  if ( !jel ) {
    /* if not, hash the edge and store its unique key */
    jel = MMG5_hashEdge(mesh,hash,ip1,ip2,3*iel+i);
    if ( !jel ) {
      printf("  # Error: %s: Unable to add edge %" MMG5_PRId " %" MMG5_PRId " within the hash table\n",
             __func__,MMG2D_indPt(mesh,ip1),MMG2D_indPt(mesh,ip2));
      return 0;
    }
  }
  else {
    /* otherwise, update the adjacency array */
    iadr = (iel-1)*3 + 1;
    adja = &mesh->adja[iadr];
    adja[i] = jel;

    j    = jel % 3;
    jel /= 3;
    iadr = (jel-1)*3 + 1;
    adja = &mesh->adja[iadr];
    adja[j] = iel*3 + i;
  }

  return 1;
}

/**  Create the cavity point ip, starting from triangle list[0];
     Return a negative value for ilist if one of the triangles of the cavity is required */
int MMG2D_cavity(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int ip,MMG5_int *list) {
  MMG5_pTria      pt,pt1,ptc;
  MMG5_pPoint     ppt;
  double          c[2],crit,dd,eps,rad,ct[6];
  MMG5_int        tref,*adja,*adjb,adj,adi,jel,iadr,nei[3],l,base; //isreq;
  int             voy,ilist,ipil;
  int             i,j;
  static int8_t   mmgWarn0=0;

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
    if ( ilist > MMG5_TRIA_LMAX - 3 ) return -1;

    ++ipil;
  }
  while ( ipil < ilist );

  ilist = MMG2D_correction_iso(mesh,ip,list,ilist,1);
  //if ( isreq ) ilist = -fabs(ilist);
  return ilist;
}

/**
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the solution (metric) structure
 * \param ip index of point to insert
 * \param list Cavity of the point \a ip.
 * \param ilist number of trias in the cavity of \a ip.
 *
 * \return 0 if the point can't be inserted, -1 if lack of memory, 1 if success.
 *
 *  Insertion in point ip in the cavity described by list.
 *
 */
int MMG2D_delone(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int ip,MMG5_int *list,int ilist) {
  MMG5_pTria      pt,pt1;
  MMG5_pPoint     ppt;
  MMG5_int        base,*adja,*adjb,iel,jel,old,iadr,size,nei[3],iadrold;
  int             i,j,k;
  MMG5_int        ielnum[3*MMG5_TRIA_LMAX+1],tref;
  int8_t          ier;
  short           i1;
  int8_t          alert;
  MMG5_Hash       hedg;
  static int8_t   mmgWarn0=0,mmgWarn1=0;

  /* Reset tagdel field */
  for (k=1; k<ilist; k++)
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
  if ( size >= 3*MMG5_TRIA_LMAX )  return 0;
  if ( !MMG5_hashNew(mesh,&hedg,size,3*size) ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to allocate hash table.\n",__func__);
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
            ier = MMG2D_hashEdgeDelone(mesh,&hedg,iel,j);
            if ( !ier ) {
              fprintf(stderr,"  ## Warning: %s: unable to update adjacency"
                      " relationship (elt %" MMG5_PRId ", edge %d).\n",
                      __func__,MMG2D_indElt(mesh,iel),j);
              return -1;
            }
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
