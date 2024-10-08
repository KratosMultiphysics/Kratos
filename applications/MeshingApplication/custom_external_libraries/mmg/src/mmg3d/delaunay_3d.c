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
 * \file mmg3d/delaunay_3d.c
 * \brief Functions for mesh modifications in Delaunay mode.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \remark Delaunay mode only (\a MMG_PATTERN flag set to \a OFF).
 * \todo doxygen documentation.
 */

#include "inlined_functions_3d_private.h"

#ifndef MMG_PATTERN

#define MMG3D_EPSRAD       1.00005
/* For Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2 test case:
   pbs with MMG3D_EPSCON=5e-4 (MMG3D does not insert enough vertex...)
*/
#define MMG3D_EPSCON       1e-5 //5.0e-4
#define MMG3D_LONMAX       4096

// uncomment to debug
// int MMG_cas;
// extern int MMG_npuiss,MMG_nvol,MMG_npres;

#define MMG3D_KTA     7
#define MMG3D_KTB    11
#define MMG3D_KTC    13

/* hash mesh edge v[0],v[1] (face i of iel) */
int MMG5_hashEdgeDelone(MMG5_pMesh mesh,MMG5_Hash *hash,MMG5_int iel,int i,MMG5_int *v) {
  MMG5_int       key,*adja,iadr,jel,mins,maxs;
  int            j;
  MMG5_hedge     *ha;

  /* compute key */
  if ( v[0] < v[1] ) {
    mins = v[0];
    maxs = v[1];
  }
  else {
    mins = v[1];
    maxs = v[0];
  }
  key = (MMG3D_KTA*(int64_t)mins + MMG3D_KTB*(int64_t)maxs)%hash->siz;
  ha  = &hash->item[key];

  if ( ha->a ) {
    /* identical face */
    if ( ha->a == mins && ha->b == maxs ) {
      iadr = (iel-1)*4 + 1;
      adja = &mesh->adja[iadr];
      adja[i] = ha->k;

      jel  = ha->k >> 2;
      j    = ha->k % 4;
      iadr = (jel-1)*4 + 1;
      adja = &mesh->adja[iadr];
      adja[j] = iel*4 + (MMG5_int)i;
      return 1;
    }
    else {
      while ( ha->nxt && ha->nxt < hash->max ) {
        ha = &hash->item[ha->nxt];
        if ( ha->a == mins && ha->b == maxs ) {
          iadr = (iel-1)*4 + 1;
          adja = &mesh->adja[iadr];
          adja[i] = ha->k;

          jel  = ha->k >> 2;
          j    = ha->k % 4;
          iadr = (jel-1)*4 + 1;
          adja = &mesh->adja[iadr];
          adja[j] = iel*4 + (MMG5_int)i;
          return 1;
        }
      }
    }
    ha->nxt   = hash->nxt;
    ha        = &hash->item[hash->nxt];
    ha->a     = mins;
    ha->b     = maxs;
    ha->k     = iel*4 + (MMG5_int)i;
    hash->nxt = ha->nxt;
    ha->nxt   = 0;

    if ( hash->nxt >= hash->max ) {
      MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,MMG5_GAP,MMG5_hedge,"face",
                         return 0;);
      for (j=hash->nxt; j<hash->max; j++)  hash->item[j].nxt = j+1;
    }
    return 1;
  }

  /* insert */
  ha->a = mins;
  ha->b = maxs;
  ha->k = iel*4 + (MMG5_int)i;
  ha->nxt = 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure.
 * \param ip index of the point to insert.
 * \param list pointer toward the list of the tetra in the cavity (computed by
 * \ref MMG5_cavity).
 * \param ilist number of tetra inside the cavity.
 * \return 1 if sucess, 0 or -1 if fail.
 *
 * Insertion of the vertex \a ip. The cavity of \a ip become its ball.
 *
 */
int MMG5_delone(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int ip,int64_t *list,int ilist) {
  MMG5_pPoint   ppt;
  MMG5_pTetra   pt,pt1;
  MMG5_xTetra   xt;
  MMG5_pxTetra  pxt0;
  MMG5_int      base,*adja,*adjb,iel,jel,old,v[3],iadr;
  int           i,j,k,l,m,size;
  MMG5_int      vois[4],iadrold,ielnum[3*MMG3D_LONMAX+1];
  short         i1;
  char          alert;
  int           isused = 0,ixt;
  MMG5_Hash     hedg;
#ifndef NDEBUG
  MMG5_int      tref;
#endif

  base = mesh->base;
  /* external faces */
  size = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tetra[old];
    iadr = (old-1)*4 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0] >> 2;
    vois[1]  = adja[1] >> 2;
    vois[2]  = adja[2] >> 2;
    vois[3]  = adja[3] >> 2;
    for (i=0; i<4; i++) {
      jel = vois[i];
      if ( (!jel) || mesh->tetra[jel].flag != base ) {
        for (j=0; j<3; j++) {
          i1  = MMG5_idir[i][j];
          ppt = &mesh->point[ pt1->v[i1] ];
          ppt->tagdel |= MG_NOM;
        }
        size++;
      }
    }
  }
  /* check isolated vertex */
  alert = 0;
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tetra[old];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      if ( !(ppt->tagdel & MG_NOM) )  alert = 1;
    }
  }
  /* reset tag */
  for (k=0; k<ilist; k++) {
    old  = list[k];
    pt1  = &mesh->tetra[old];
    for (i=0; i<4; i++) {
      ppt = &mesh->point[ pt1->v[i] ];
      ppt->tagdel &= ~MG_NOM;
    }
  }
  if ( alert )  {return 0;}
  /* hash table params */
  if ( size > 3*MMG3D_LONMAX )  return 0;
  if ( !MMG5_hashNew(mesh,&hedg,size,3*size) ) { /*3*size suffit */
    fprintf(stderr,"\n  ## Error: %s: unable to complete mesh.\n",__func__);
    return -1;
  }

  /*tetra allocation : we create "size" tetra*/
  ielnum[0] = size;
  for (k=1 ; k<=size ; k++) {
    ielnum[k] = MMG3D_newElt(mesh);

    if ( !ielnum[k] ) {
      MMG3D_TETRA_REALLOC(mesh,ielnum[k],mesh->gap,
                          fprintf(stderr,"\n  ## Error: %s: unable to allocate a"
                                  " new element.\n",__func__);
                          return -1);
    }
  }

  size = 1;
  for (k=0; k<ilist; k++) {
    old  = list[k];

    iadrold = (old-1)*4 + 1;
    adja = &mesh->adja[iadrold];
    vois[0]  = adja[0];
    vois[1]  = adja[1];
    vois[2]  = adja[2];
    vois[3]  = adja[3];

    pt   = &mesh->tetra[old];
    if(pt->xt) {
      pxt0 = &mesh->xtetra[pt->xt];
      memcpy(&xt,pxt0,sizeof(MMG5_xTetra));
      isused=0;
      ixt = 1;
    } else {
      ixt = 0;
    }

    for (i=0; i<4; i++) {
      jel = vois[i] /4;
      j   = vois[i] % 4;

      /* external face */
      if ( !jel || (mesh->tetra[jel].flag != base) ) {
        iel = ielnum[size++];
        assert(iel);

        pt1 = &mesh->tetra[iel];
        memcpy(pt1,pt,sizeof(MMG5_Tetra));
        pt1->v[i] = ip;
        pt1->qual = MMG5_orcal(mesh,sol,iel);
        pt1->ref = mesh->tetra[old].ref;
        pt1->mark = mesh->mark;
        iadr = (iel-1)*4 + 1;
        adjb = &mesh->adja[iadr];
        adjb[i] = adja[i];

        if(ixt) {
          if( xt.ref[i] || xt.ftag[i]) {
            if(!isused) {
              pt1->xt = pt->xt;
              pt->xt = 0;
              pxt0 = &mesh->xtetra[pt1->xt];
              memset(pxt0,0,sizeof(MMG5_xTetra));
              pxt0->ref[i]   = xt.ref[i] ; pxt0->ftag[i]  = xt.ftag[i];
              pxt0->edg[MMG5_iarf[i][0]] = xt.edg[MMG5_iarf[i][0]];
              pxt0->edg[MMG5_iarf[i][1]] = xt.edg[MMG5_iarf[i][1]];
              pxt0->edg[MMG5_iarf[i][2]] = xt.edg[MMG5_iarf[i][2]];
              pxt0->tag[MMG5_iarf[i][0]] = xt.tag[MMG5_iarf[i][0]];
              pxt0->tag[MMG5_iarf[i][1]] = xt.tag[MMG5_iarf[i][1]];
              pxt0->tag[MMG5_iarf[i][2]] = xt.tag[MMG5_iarf[i][2]];
              pxt0->ori = xt.ori;
              isused=1;
            } else {
              mesh->xt++;
              if ( mesh->xt > mesh->xtmax ) {
                MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                                   "larger xtetra table",
                                   mesh->xt--;
                                   fprintf(stderr,"  Exit program.\n"); return -1;);
              }
              pt1->xt = mesh->xt;
              pxt0 = &mesh->xtetra[pt1->xt];
              pxt0->ref[i]   = xt.ref[i] ; pxt0->ftag[i]  = xt.ftag[i];
              pxt0->edg[MMG5_iarf[i][0]] = xt.edg[MMG5_iarf[i][0]];
              pxt0->edg[MMG5_iarf[i][1]] = xt.edg[MMG5_iarf[i][1]];
              pxt0->edg[MMG5_iarf[i][2]] = xt.edg[MMG5_iarf[i][2]];
              pxt0->tag[MMG5_iarf[i][0]] = xt.tag[MMG5_iarf[i][0]];
              pxt0->tag[MMG5_iarf[i][1]] = xt.tag[MMG5_iarf[i][1]];
              pxt0->tag[MMG5_iarf[i][2]] = xt.tag[MMG5_iarf[i][2]];
              pxt0->ori = xt.ori;
            }
          }
          else {
            pt1->xt = 0;
          }
        }

        if ( jel ) {
          iadr = (jel-1)*4 + 1;
          adjb = &mesh->adja[iadr];
          adjb[j] = iel*4 + i;
        }

        /* internal faces (p1,p2,ip) */
        for (j=0; j<4; j++) {
          if ( j != i ) {
            m = 0;
            for (l=0; l<3; l++)
              if ( pt1->v[ MMG5_idir[j][l] ] != ip ) {
                v[m] = pt1->v[ MMG5_idir[j][l] ];
                m++;
              }
            MMG5_hashEdgeDelone(mesh,&hedg,iel,j,v);
          }
        }
      }
    }
  }

  /* remove old tetra */
#ifndef NDEBUG
  tref = mesh->tetra[list[0]].ref;
#endif
  for (k=0; k<ilist; k++) {
    assert(tref==mesh->tetra[list[k]].ref);
    if ( !MMG3D_delElt(mesh,list[k]) ) return -1;
  }

  // ppt = &mesh->point[ip];
  // ppt->flag = mesh->flag;
  MMG5_DEL_MEM(mesh,hedg.item);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the met structure
 * \param ip index of the point to insert
 * \param list poiner toward the cavity of the point
 * \param ilist number of elts in the cavity
 * \param nedep ???
 * \param volmin minimal authorized volume
 *
 * \return 0 if fail (unused point), the size of the corrected cavity otherwise
 *
 * Cavity correction for quality (aniso).
 *
 */
static int MMG5_correction_ani(MMG5_pMesh mesh,MMG5_pSol met,int ip,int64_t* list,
                                int ilist,int nedep,double volmin) {
  MMG5_pPoint   ppt,p1,p2,p3;
  MMG5_pTetra   pt;
  double        dd,det,nn,eps,eps2,ux,uy,uz,vx,vy,vz,v1,v2,v3;
  double        *ma,*mb,*mc,*md,mm[6],h1,h2,h3;
  MMG5_int      *adja,i,j,iel,iadr,adj,ib,ic,id;
  MMG5_int      base,vois[4];
  int           ipil,lon,ncor;

  ppt  = &mesh->point[ip];
  if ( ppt->tag & MG_NUL )  return ilist;
  base = mesh->base;
  lon  = ilist;
  eps  = MMG3D_EPSCON;
  eps2 = eps*eps;

  /* average metric */
  memset(mm,0,6*sizeof(double));
  ma   = &met->m[6*ip];

  do {
    ipil = lon-1;
    ncor = 0;

    while ( ipil >= 0 ) {
      iel  = list[ipil];
      iadr = (iel-1)*4 + 1;
      adja = &mesh->adja[iadr];
      vois[0]  = adja[0] >> 2;
      vois[1]  = adja[1] >> 2;
      vois[2]  = adja[2] >> 2;
      vois[3]  = adja[3] >> 2;
      pt   = &mesh->tetra[iel];

      // MMG_cas=0; // uncomment to debug
      for (i=0; i<4; i++) {
        adj = vois[i];
        // MMG_cas = 0;// uncomment to debug
        if ( adj && mesh->tetra[adj].flag == base)  continue;

        ib = pt->v[ MMG5_idir[i][0] ];
        ic = pt->v[ MMG5_idir[i][1] ];
        id = pt->v[ MMG5_idir[i][2] ];

        p1 = &mesh->point[ib];
        p2 = &mesh->point[ic];
        p3 = &mesh->point[id];

        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];
        uz = p2->c[2] - p1->c[2];

        vx = p3->c[0] - p1->c[0];
        vy = p3->c[1] - p1->c[1];
        vz = p3->c[2] - p1->c[2];

        /* volume PABC */
        v1 = uz*vy - uy*vz;
        v2 = ux*vz - uz*vx;
        v3 = uy*vx - ux*vy;
        dd = v1*(ppt->c[0]-p1->c[0]) + v2*(ppt->c[1]-p1->c[1]) \
          + v3*(ppt->c[2]-p1->c[2]);
        // MMG_cas=1; // uncomment to debug

        /*test sur le volume avec un eps local*/
        h1 = ux*ux + uy*uy + uz*uz;
        h2 = vx*vx + vy*vy + vz*vz;
        h3 = (p2->c[0] - p3->c[0])*(p2->c[0] - p3->c[0]) + (p2->c[1] - p3->c[1])*(p2->c[1] - p3->c[1])
          + (p2->c[2] - p3->c[2])*(p2->c[2] - p3->c[2]);
        if ( dd < volmin*sqrt(h1*h2*h3) )  break;

        /* average metric */
        mb   = &met->m[6*ib];
        mc   = &met->m[6*ic];
        md   = &met->m[6*id];
        for (j=0; j<6; j++)
          mm[j] = 0.25 * (ma[j]+mb[j]+mc[j]+md[j]);

        det = mm[0] * ( mm[3]*mm[5] - mm[4]*mm[4]) \
          - mm[1] * ( mm[1]*mm[5] - mm[2]*mm[4]) \
          + mm[2] * ( mm[1]*mm[4] - mm[2]*mm[3]);
        if ( det < MMG5_EPSOK )  break;

        /* point close to face */
        // MMG_cas=2; // uncomment to debug
        nn = mm[0]*v1*v1 + mm[3]*v2*v2 + mm[5]*v3*v3 \
          + 2.0*(mm[1]*v1*v2 + mm[2]*v1*v3 + mm[4]*v2*v3);

        if ( det*dd*dd < nn * eps2 )  break;
        // MMG_cas=0; // uncomment to debug
      }
      if ( i < 4 || pt->tag & MG_REQ ) {
        if ( ipil <= nedep )   {
          return 0;
        }
        /* remove iel from list */
        pt->flag = base-1;
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
 * \param met pointer toward the met structure
 * \param ip index of the point to insert
 * \param list poiner toward the cavity of the point
 * \param ilist number of elts in the cavity
 * \param nedep ???
 * \param volmin minimal authorized volume
 *
 * \return 0 if fail (unused point), the size of the corrected cavity otherwise
 *
 * Cavity correction for quality (iso).
 *
 */
static int
MMG5_correction_iso(MMG5_pMesh mesh,int ip,int64_t *list,int ilist,int nedep,double volmin) {
  MMG5_pPoint ppt,p1,p2,p3;
  MMG5_pTetra      pt;
  double           dd,nn,eps,eps2,ux,uy,uz,vx,vy,vz,v1,v2,v3;
  MMG5_int         *adja,i,iel,iadr,adj,ib,ic,id;
  MMG5_int         base,vois[4];
  int              ncor,ipil,lon;

  ppt  = &mesh->point[ip];
  if ( ppt->tag & MG_NUL )  return ilist;
  base = mesh->base;
  lon  = ilist;
  eps  = MMG3D_EPSCON;
  eps2 = eps*eps;
  do {
    ipil = lon-1;
    ncor = 0;

    while ( ipil >= 0 ) {
      iel  = list[ipil];
      iadr = (iel-1)*4 + 1;
      adja = &mesh->adja[iadr];
      vois[0]  = adja[0] >> 2;
      vois[1]  = adja[1] >> 2;
      vois[2]  = adja[2] >> 2;
      vois[3]  = adja[3] >> 2;
      pt   = &mesh->tetra[iel];
      // MMG_cas=0; // uncomment to debug
      for (i=0; i<4; i++) {
        adj = vois[i];
        // MMG_cas = 0; // uncomment to debug
        if ( adj && mesh->tetra[adj].flag == base )  continue;

        ib = pt->v[ MMG5_idir[i][0] ];
        ic = pt->v[ MMG5_idir[i][1] ];
        id = pt->v[ MMG5_idir[i][2] ];

        p1 = &mesh->point[ib];
        p2 = &mesh->point[ic];
        p3 = &mesh->point[id];

        ux = p2->c[0] - p1->c[0];
        uy = p2->c[1] - p1->c[1];
        uz = p2->c[2] - p1->c[2];

        vx = p3->c[0] - p1->c[0];
        vy = p3->c[1] - p1->c[1];
        vz = p3->c[2] - p1->c[2];

        /* volume PABC */
        v1 = uz*vy - uy*vz;
        v2 = ux*vz - uz*vx;
        v3 = uy*vx - ux*vy;
        dd = v1*(ppt->c[0]-p1->c[0]) + v2*(ppt->c[1]-p1->c[1]) \
          + v3*(ppt->c[2]-p1->c[2]);
        // MMG_cas=1; // uncomment to debug

        if ( dd < volmin )  break;

        /* point close to face */
        nn = (v1*v1 + v2*v2 + v3*v3);

        // MMG_cas=2; // uncomment to debug
        if ( dd*dd < nn * eps2 )  break;
        // MMG_cas=0; //uncomment to debug
      }
      if ( i < 4 ||  pt->tag & MG_REQ ) {
        if ( ipil <= nedep )  {
          return 0;
        }

        /* remove iel from list */
        pt->flag = base-1;
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
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param iel tetra index.
 * \param ip point local index in \a iel.
 * \param list pointer toward the list of tetra in the shell of edge where
 * ip will be inserted.
 * \param lon number of tetra in the list.
 * \return ilist number of tetra inside the cavity or -ilist if one of the tet
 * of the cavity is required.
 *
 * Mark elements in cavity and update the list of tetra in the cavity.
 *
 */
int MMG5_cavity_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int iel,int ip,int64_t* list,int lon,double volmin) {
  MMG5_pPoint    ppt;
  MMG5_pTetra    pt,pt1;
  double         c[3],eps,dd,ray,ux,uy,uz,crit;
  double         *mj,*mp,ct[12];
  MMG5_int       *adja,*adjb,adj,adi,ia,jel,iadr;
  MMG5_int       base,vois[4],tref;
  int            k,voy,ilist,ipil,i,j,isreq,l;

  if ( lon < 1 )  return 0;
  ppt = &mesh->point[ip];
  if ( ppt->tag & MG_NUL )  return 0;
  base  = ++mesh->base;

  isreq = 0;

  tref = mesh->tetra[list[0]/6].ref;
  for (k=0; k<lon; k++) {
    mesh->tetra[list[k]/6].flag = base;

    if ( !mesh->info.opnbdy ) {
      if ( tref != mesh->tetra[list[k]/6].ref ) {
        return 0;
      }
    }
    else {
      pt = &mesh->tetra[list[k]/6];
      if ( pt->xt ) {
        l  = list[k]%6;
        if ( (mesh->xtetra[pt->xt].ftag[MMG5_ifar[l][0]] & MG_BDY) ||
             (mesh->xtetra[pt->xt].ftag[MMG5_ifar[l][1]] & MG_BDY) )
          return 0;
      }
    }
  }
  for (k=0; k<lon; k++)
    list[k] = list[k] / 6;

  /* grow cavity by adjacency */
  eps   = MMG3D_EPSRAD * MMG3D_EPSRAD;
  ilist = lon;
  ipil  = 0;
  iadr  = ip*6;
  mp    = &met->m[iadr];

  do {
    jel  = list[ipil];
    iadr = (jel-1)*4 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0];
    vois[1]  = adja[1];
    vois[2]  = adja[2];
    vois[3]  = adja[3];

    for (i=0; i<4; i++) {
      adj = vois[i];

      if ( !adj )  continue;

      adj >>= 2;
      voy = vois[i] % 4;
      pt  = &mesh->tetra[adj];

      /* boundary face */
      if ( pt->flag == base )  continue;
      if ( pt->xt && (mesh->xtetra[pt->xt].ftag[voy] & MG_BDY) ) continue;

      for (j=0,l=0; j<4; j++,l+=3) {
        memcpy(&ct[l],mesh->point[pt->v[j]].c,3*sizeof(double));
      }


      /* Delaunay kernel */
      if ( !MMG5_cenrad_ani(mesh,ct,mp,c,&ray) )  continue;

      ux = ppt->c[0] - c[0];
      uy = ppt->c[1] - c[1];
      uz = ppt->c[2] - c[2];
      dd =      mp[0]*ux*ux + mp[3]*uy*uy + mp[5]*uz*uz \
        + 2.0*(mp[1]*ux*uy + mp[2]*ux*uz + mp[4]*uy*uz);
      crit = eps * ray;
      if ( dd > crit )  continue;

      /* mixed metrics */
      crit = sqrt(dd/ray);
      for (j=0; j<4; j++) {
        ia   = pt->v[j];
        iadr = 6*ia;
        mj   = &met->m[iadr];
        if ( !MMG5_cenrad_ani(mesh,ct,mj,c,&ray) )  continue;
        ux = ppt->c[0] - c[0];
        uy = ppt->c[1] - c[1];
        uz = ppt->c[2] - c[2];
        dd =      mj[0]*ux*ux + mj[3]*uy*uy + mj[5]*uz*uz \
          + 2.0*(mj[1]*ux*uy + mj[2]*ux*uz + mj[4]*uy*uz);
        crit += sqrt(dd/ray);
      }
      crit *= MMG3D_EPSRAD;
      if ( crit > 5.0 ) continue;

      /* lost face(s) */
      iadr = (adj-1)*4 + 1;
      adjb = &mesh->adja[iadr];

      for (j=0; j<4; j++) {
        if ( j == voy )  continue;
        adi = adjb[j];

        if ( !adi )  continue;

        adi >>= 2;
        assert(adi !=jel);

        pt1 = &mesh->tetra[adi];
        if ( pt1->flag == base ) {
          if ( pt1->xt && (mesh->xtetra[pt1->xt].ftag[adjb[j]%4] & MG_BDY) ) break;
        }
      }
      /* store tetra */
      if ( j == 4 ) {
        if ( pt->tag & MG_REQ ) isreq = 1;
        pt->flag = base;
        list[ilist++] = adj;
      }
    }
    if ( ilist > MMG3D_LONMAX - 3 )  return -1;
    ++ipil;
  }
  while ( ipil < ilist );

  /* global overflow */
  ilist = MMG5_correction_ani(mesh,met,ip,list,ilist,lon,volmin);

  if ( isreq ) ilist = -abs(ilist);

  // uncomment to debug
  /* if(MMG_cas==1) MMG_nvol++; */
  /* else if(MMG_cas==2 || MMG_cas>20) { */
  /*   MMG_npuiss++; */
  /*   if(MMG_cas>20) MMG_npres++; */
  /* } */

  return ilist;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \param iel tetra index.
 * \param ip point local index in \a iel.
 * \param list pointer toward the list of tetra in the shell of edge where
 * ip will be inserted.
 * \param lon number of tetra in the list.
 * \return ilist number of tetra inside the cavity or -ilist if one of the tet
 * of the cavity is required.
 *
 * Mark elements in cavity and update the list of tetra in the cavity.
 *
 */
int MMG5_cavity_iso(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int iel,int ip,int64_t *list,int lon,double volmin) {
  MMG5_pPoint      ppt;
  MMG5_pTetra      pt,pt1;
  double           c[3],crit,dd,eps,ray,ct[12];
  MMG5_int         *adja,*adjb,adj,adi,jel,iadr;
  MMG5_int         tref,base,vois[4],l;
  int              isreq;
  int              voy,i,j,k,ipil,ilist;

  if ( lon < 1 )  return 0;
  ppt = &mesh->point[ip];
  if ( ppt->tag & MG_NUL )  return 0;
  base  = ++mesh->base;

  isreq = 0;

  tref = mesh->tetra[list[0]/6].ref;
  for (k=0; k<lon; k++) {
    mesh->tetra[list[k]/6].flag = base;

    if ( !mesh->info.opnbdy ) {
      if ( tref != mesh->tetra[list[k]/6].ref ) {
        return 0;
      }
    }
    else {
      pt = &mesh->tetra[list[k]/6];
      if ( pt->xt ) {
        l  = list[k]%6;
        if ( (mesh->xtetra[pt->xt].ftag[MMG5_ifar[l][0]] & MG_BDY) ||
             (mesh->xtetra[pt->xt].ftag[MMG5_ifar[l][1]] & MG_BDY) )
          return 0;
      }
    }
  }

  for (k=0; k<lon; k++)
    list[k] = list[k] / 6;

  /* grow cavity by adjacency */
  eps   = MMG3D_EPSRAD*MMG3D_EPSRAD;
  ilist = lon;
  ipil  = 0;

  do {
    jel  = list[ipil];
    iadr = (jel-1)*4 + 1;
    adja = &mesh->adja[iadr];
    vois[0]  = adja[0];
    vois[1]  = adja[1];
    vois[2]  = adja[2];
    vois[3]  = adja[3];

    for (i=0; i<4; i++) {
      adj = vois[i];
      if ( !adj )  continue;

      adj >>= 2;
      voy = vois[i] % 4;
      pt  = &mesh->tetra[adj];
      /* boundary face */
      if ( pt->flag == base )  continue;
      if ( pt->xt && (mesh->xtetra[pt->xt].ftag[voy] & MG_BDY) ) continue;

      for (j=0,l=0; j<4; j++,l+=3) {
        memcpy(&ct[l],mesh->point[pt->v[j]].c,3*sizeof(double));
      }

      if ( !MMG5_cenrad_iso(mesh,ct,c,&ray) )  continue;
      crit = eps * ray;

      /* Delaunay criterion */
      dd = (ppt->c[0] - c[0]) * (ppt->c[0] - c[0]) \
        + (ppt->c[1] - c[1]) * (ppt->c[1] - c[1]) \
        + (ppt->c[2] - c[2]) * (ppt->c[2] - c[2]);
      if ( dd > crit )  continue;

      /* lost face(s) */
      iadr = (adj-1)*4 + 1;
      adjb = &mesh->adja[iadr];

      for (j=0; j<4; j++) {
        if ( j == voy )  continue;
        adi = adjb[j];
        if ( !adi )  continue;

        adi >>= 2;
        assert(adi !=jel);

        pt1 = &mesh->tetra[adi];
        if ( pt1->flag == base )
          if ( pt1->xt && (mesh->xtetra[pt1->xt].ftag[adjb[j]%4] & MG_BDY) ) break;
      }
      /* store tetra */
      if ( j == 4 ) {
        if ( pt->tag & MG_REQ ) isreq = 1;
        pt->flag = base;
        list[ilist++] = adj;
      }
    }
    if ( ilist > MMG3D_LONMAX - 3 ) return -1;

    ++ipil;
  }
  while ( ipil < ilist );

  /* global overflow: obsolete avec la reallocation */
  ilist = MMG5_correction_iso(mesh,ip,list,ilist,lon,volmin);

  if ( isreq ) ilist = -abs(ilist);

  // uncomment to debug
  /* if(MMG_cas==1) MMG_nvol++; */
  /* else if(MMG_cas==2 || MMG_cas>20) { */
  /*   MMG_npuiss++; */
  /*   if(MMG_cas>20) MMG_npres++; */
  /* } */
  return ilist;
}

#endif
