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

/**
 * \param mesh pointer toward the mesh structure.
 * \param severe level of performed check
 * \param base unused argument.
 * \return 0 if fail, 1 if success.
 *
 * Check the mesh validity
 *
 */
int MMG5_mmg2dChkmsh(MMG5_pMesh mesh, int severe,int base) {
  MMG5_pPoint    ppt;
  MMG5_pTria     pt1,pt2;
  MMG5_pEdge     ped;
  int           *adja,*adja1,adj,adj1,k,i,iadr;
  int            kk,l,nk,j,ip,lon,len;
  int           *list;
  unsigned char  voy,voy1;
  static char    mmgErr0=0,mmgErr1=0,mmgErr2=0,mmgErr3=0,mmgErr4=0,mmgErr5=0;
  static char    mmgErr6=0;

  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if ( !MG_EOK(pt1) )  continue;
    iadr = (k-1)*3 + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<3; i++) {
      adj = adja[i] / 3;
      voy = adja[i] % 3;
      if ( !adj )  continue;

      if ( adj == k ) {
        if ( !mmgErr0 ) {
          mmgErr0 = 1;
          fprintf(stderr,"\n  ## Error: %s: 1. at least 1 wrong"
                  " adjacency %d %d\n",__func__,MMG2D_indElt(mesh,k),
                  MMG2D_indElt(mesh,adj));
          fprintf(stderr,"vertices of %d: %d %d %d \n",MMG2D_indElt(mesh,k),
                  MMG2D_indPt(mesh,pt1->v[0]),MMG2D_indPt(mesh,pt1->v[1]),
                  MMG2D_indPt(mesh,pt1->v[2]));
          fprintf(stderr,"adj of %d: %d %d %d \n",
                  k,adja[0]/3,adja[1]/3,adja[2]/3);
        }
        return 0;
      }
      pt2 = &mesh->tria[adj];
      if ( !MG_EOK(pt2) ) {
        if ( !mmgErr1 ) {
          mmgErr1 = 1;
          fprintf(stderr,"\n  ## Error: %s: 4. at least 1 invalid"
                  " adjacent %d %d\n",__func__,MMG2D_indElt(mesh,adj),
                  MMG2D_indElt(mesh,k));
          fprintf(stderr,"vertices of %d: %d %d %d\n",
                  MMG2D_indElt(mesh,k),MMG2D_indPt(mesh,pt1->v[0]),
                  MMG2D_indPt(mesh,pt1->v[1]),MMG2D_indPt(mesh,pt1->v[2]));
          fprintf(stderr,"vertices adj %d: %d %d %d \n",
                  adj,pt2->v[0],pt2->v[1],pt2->v[2]);
          fprintf(stderr,"adj of %d: %d %d %d\n",MMG2D_indElt(mesh,k),
                  MMG2D_indElt(mesh,adja[0]/3),MMG2D_indElt(mesh,adja[1]/3),
                  MMG2D_indElt(mesh,adja[2]/3));
        }
        return 0;
      }
      iadr  = (adj-1)*3 + 1;
      adja1 = &mesh->adja[iadr];
      adj1  = adja1[voy] / 3;
      voy1  = adja1[voy] % 3;
      if ( adj1 != k || voy1 != i ) {
        if ( !mmgErr2 ) {
          mmgErr2 = 1;
          fprintf(stderr,"\n  ## Error: %s: 2. at least 1 wrong adjacency"
                  " %d %d\n",__func__,MMG2D_indElt(mesh,k),MMG2D_indElt(mesh,adj1));
          fprintf(stderr,"vertices of %d: %d %d %d \n",MMG2D_indElt(mesh,k),
                  MMG2D_indPt(mesh,pt1->v[0]),MMG2D_indPt(mesh,pt1->v[1]),
                  MMG2D_indPt(mesh,pt1->v[2]));
          fprintf(stderr,"adj(k) %d: %d %d %d \n",MMG2D_indElt(mesh,adj),
                  MMG2D_indPt(mesh,pt2->v[0]),MMG2D_indPt(mesh,pt2->v[1]),
                  MMG2D_indPt(mesh,pt2->v[2]));
          fprintf(stderr,"adj(%d): %d %d %d\n",
                  MMG2D_indElt(mesh,k),MMG2D_indElt(mesh,adja[0]/3),
                  MMG2D_indElt(mesh,adja[1]/3),MMG2D_indElt(mesh,adja[2]/3));
          fprintf(stderr,"adj(%d): %d %d %d %d\n",
                  MMG2D_indElt(mesh,adj),MMG2D_indElt(mesh,adja1[0]/3),
                  MMG2D_indElt(mesh,adja1[1]/3),MMG2D_indElt(mesh,adja1[2]/3),
                  MMG2D_indElt(mesh,adja1[3]/3));
        }
        return 0;
      }

      /*chk edge*/
      if(pt1->edg[i]) {
        ped = &mesh->edge[pt1->edg[i]];
        if ( !mmgErr3 ) {
          mmgErr3 = 1;
          if(!(((ped->a==pt1->v[MMG2D_iare[i][0]]) || (ped->a==pt1->v[MMG2D_iare[i][1]]))
               || ((ped->b==pt1->v[MMG2D_iare[i][0]]) || (ped->b==pt1->v[MMG2D_iare[i][1]])))) {
            fprintf(stderr,"\n  ## Error: %s: 3. at least 1 wrong edge in triangle %d\n",
                    __func__,MMG2D_indElt(mesh,k));
            fprintf(stderr,"vertices of %d: %d %d %d \n",MMG2D_indElt(mesh,k),
                    MMG2D_indPt(mesh,pt1->v[0]),MMG2D_indPt(mesh,pt1->v[1]),
                    MMG2D_indPt(mesh,pt1->v[2]));
          }
          fprintf(stderr,"edge %d : %d %d\n",i,ped->a,ped->b);
          return 0;
        }
      }

    }
  }

  if ( !severe )  return 1;

  MMG5_SAFE_CALLOC(list,MMG2D_LMAX,int,return 0);

  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if ( !MG_EOK(pt1) )  continue;
    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<3; i++) {
      adj = (adja[i]-1) / 3 + 1;
      if ( !adj )  continue;

      ip  = pt1->v[i];
      ppt = &mesh->point[ip];
      if ( !MG_VOK(ppt) ) {
        if ( !mmgErr4 ) {
          mmgErr4 = 1;
          fprintf(stderr,"\n  ## Error: %s: 6. at least 1 unused vertex %d  %d\n",__func__,
                  MMG2D_indElt(mesh,k),MMG2D_indPt(mesh,ip));
          fprintf(stderr,"%d %d %d\n",MMG2D_indPt(mesh,pt1->v[0]),
                  MMG2D_indPt(mesh,pt1->v[1]),MMG2D_indPt(mesh,pt1->v[2]));
        }
        MMG5_SAFE_FREE(list);
        return 0;
      }
      lon = MMG2D_boulep(mesh,k,i,list);
      for (l=1; l<=lon; l++) {
        kk  = list[l] / 3;
        nk  = list[l] % 3;
        pt2 = &mesh->tria[kk];
        if ( pt2->v[nk] != ip ) {
          if ( !mmgErr5 ) {
            mmgErr5 = 1;
            fprintf(stderr,"\n  ## Error: %s: 5. at least 1 wrong ball %d, %d\n",
                    __func__,MMG2D_indPt(mesh,ip),MMG2D_indPt(mesh,pt2->v[nk]));
          }
          MMG5_SAFE_FREE(list);
          return 0;
        }
      }
      if ( lon < 1 )  continue;
      len = 0;
      for (kk=1; kk<=mesh->nt; kk++) {
        pt2 = &mesh->tria[kk];
        if ( !pt2->v[0] )  continue;
        for (j=0; j<3; j++)
          if ( pt2->v[j] == ip ) {
            len++;
            break;
          }
      }
      if ( len != lon ) {
        if ( !mmgErr6 ) {
          mmgErr6 = 1;
          fprintf(stderr,"\n  ## Error: %s: 7. at least 1 incorrect ball %d: %d %d\n",
                  __func__,MMG2D_indPt(mesh,pt1->v[i]),lon,len);
        }
        MMG5_SAFE_FREE(list);
        return 0;
      }
    }
  }
  MMG5_SAFE_FREE(list);
  return 1;
}

/* Check of adjacency relations and edge tags */
int MMG2D_chkmsh(MMG5_pMesh mesh) {
  MMG5_pTria        pt,pt1;
  MMG5_pPoint       p1,p2;
  int               *adja,*adjaj,k,jel;
  char              i,i1,i2,j;
  static char       mmgErr0=0,mmgErr1=0,mmgErr2=0,mmgErr3=0,mmgErr4=0;
  static char       mmgErr6=0,mmgErr5=0;

  /* Check adjacencies */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];

    for (i=0; i<3; i++) {
      jel = adja[i] / 3;
      j   = adja[i] % 3;

      if ( !jel ) {
        if ( !(pt->tag[i] & MG_GEO ) ) {
          if ( !mmgErr0 ) {
            mmgErr0 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong  edge tag"
                    " (edge %d in tria %d.)\n",
                    __func__,i,MMG2D_indElt(mesh,k));
          }
          return 0;
        }
      }
      else {
        adjaj = &mesh->adja[3*(jel-1)+1];
        if ( adjaj[j] / 3 != k ) {
         if ( !mmgErr1 ) {
            mmgErr1 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong adjacency"
                    " (%d %d).\n",__func__,MMG2D_indElt(mesh,k),
                    MMG2D_indElt(mesh,jel));
         }
          return 0;
        }
      }
    }
  }

  /* Check consistency between tags of edges and vertices */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];
      if ( pt->tag[i] & MG_GEO ) {
        if ( !(mesh->point[pt->v[i1]].tag & MG_GEO) && !( MG_SIN(mesh->point[pt->v[i1]].tag) )) {
         if ( !mmgErr2 ) {
            mmgErr2 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %d: edge %d, vertex %d)\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i1]));
         }
          return 0;
        }
        if ( !(mesh->point[pt->v[i2]].tag & MG_GEO) && !( MG_SIN(mesh->point[pt->v[i2]].tag) )) {
          if ( !mmgErr2 ) {
            mmgErr2 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %d: edge %d, vertex %d)\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i2]));
          }
          return 0;
        }
      }

      if ( pt->tag[i] & MG_REF ) {
        if ( !(mesh->point[pt->v[i1]].tag & MG_REF) && !( MG_SIN(mesh->point[pt->v[i1]].tag)) ) {
          if ( !mmgErr3 ) {
            mmgErr3 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %d: edge %d, vertex %d)\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i1]));
          }
          return 0;
        }
        if ( !(mesh->point[pt->v[i2]].tag & MG_REF) && !( MG_SIN(mesh->point[pt->v[i2]].tag)) ) {
          if ( !mmgErr3 ) {
            mmgErr3 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %d: edge %d, vertex %d)\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i2]));
          }
          return 0;
        }
      }

    }
  }

  /* Check consistency between edge tags and triangle refs */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      jel = adja[i] / 3;

      if ( ( pt->tag[i] & MG_GEO ) && jel ) {
        if ( !mmgErr4 ) {
          mmgErr4 = 1;
          fprintf(stderr,"\n  ## Error: %s: at least 1 edge tagged boundary"
                  " while it has a neighbour (%d %d).\n",__func__,
                  MMG2D_indPt(mesh,pt->v[i1]),MMG2D_indPt(mesh,pt->v[i2]));
        }
        return 0;
      }

      if ( pt->tag[i] & MG_REF ) {
        pt1 = &mesh->tria[jel];
        if ( pt->ref == pt1->ref ) {
          if ( !mmgErr5 ) {
            mmgErr5 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge tagged ref while"
                    " both corresponding triangles have same ref (%d %d).\n",
                    __func__,
                    MMG2D_indPt(mesh,pt->v[i1]), MMG2D_indPt(mesh,pt->v[i2]));
          }
          {
            fprintf(stderr,"Saving mesh...\n");
            if ( !MMG2D_hashTria(mesh) ) {
              fprintf(stdout,"  ## Hashing problem. Exit program.\n");
              return 0;
            }

            MMG2D_bdryEdge(mesh);
            MMG2D_savemesh_db(mesh,mesh->nameout,0);
            return 0;
          }

          return 0;
        }
      }
    }
  }

  /* Check consistency between REF, GEO and BDY tags between edges and points */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_GEO || pt->tag[i] & MG_REF ) {
        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];

        if ( !(pt->tag[i] & MG_BDY) ) {
          if ( !mmgErr6 ) {
            mmgErr6 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge (%d %d) tagged"
                    " %d, but not MG_BDY\n",__func__,MMG2D_indPt(mesh,pt->v[i1]),
                    MMG2D_indPt(mesh,pt->v[i2]),pt->tag[i]);
          }
          return 0;
        }

        p1 = &mesh->point[pt->v[i1]];
        p2 = &mesh->point[pt->v[i2]];

        if ( !(p1->tag & MG_BDY) ) {
          if ( !mmgErr6 ) {
            mmgErr6 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge (%d %d) tagged %d,"
                    " with a point (%d) not tagged BDY.\n",__func__,
                    MMG2D_indPt(mesh,pt->v[i1]),MMG2D_indPt(mesh,pt->v[i2]),
                    pt->tag[i],MMG2D_indPt(mesh,pt->v[i1]));
          }
          return 0;
        }

        if ( !(p2->tag & MG_BDY) ) {
          if ( !mmgErr6 ) {
            mmgErr6 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge (%d %d) tagged %d,"
                    " with a point (%d) not tagged BDY.\n",__func__,
                    MMG2D_indPt(mesh,pt->v[i1]),MMG2D_indPt(mesh,pt->v[i2]),
                    pt->tag[i],MMG2D_indPt(mesh,pt->v[i2]));
          }
          return 0;
        }
      }
    }
  }

  return 1;
}

/* Check orientation of elements in the mesh */
int MMG2D_chkor(MMG5_pMesh mesh) {
  MMG5_pTria        pt;
  MMG5_pPoint       p0,p1,p2;
  double            det;
  int               k;

  for (k=1; k<=mesh->np; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    p0 = &mesh->point[pt->v[0]];
    p1 = &mesh->point[pt->v[1]];
    p2 = &mesh->point[pt->v[2]];

    det = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);

    if ( det <= 0.0 ) return 0;
  }

  return 1;
}
