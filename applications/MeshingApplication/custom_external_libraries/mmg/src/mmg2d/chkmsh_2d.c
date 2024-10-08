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

#include "libmmg2d.h"
#include "libmmg2d_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param severe level of performed check
 * \param base 1 if we want to test opnbdy edge tags (consistent only after analysis)
 * \return 0 if fail, 1 if success.
 *
 * Check the mesh validity
 *
 */
int MMG5_mmg2dChkmsh(MMG5_pMesh mesh, int severe,MMG5_int base) {
  MMG5_pPoint    ppt;
  MMG5_pTria     pt,pt1,pt2;
  MMG5_int       *adja,*adja1,adj,adj1,k,iadr;
  int            i,j,lon,len;
  MMG5_int       kk,l,nk,ip;
  MMG5_int       *list;
  uint8_t        voy,voy1;
  static int8_t  mmgErr0=0,mmgErr1=0,mmgErr2=0,mmgErr3=0,mmgErr4=0,mmgErr5=0;
  static int8_t  mmgErr6=0,mmgErr7=0,mmgErr8=0,mmgErr9=0,mmgErr10=0,mmgErr11=0;
  static int8_t  mmgErr12=0,mmgErr13=0,mmgErr14=0,mmgErr15=0;

  /** Check adjacency relationships */
  for (k=1; k<=mesh->nt; k++) {
    pt1 = &mesh->tria[k];
    if ( !MG_EOK(pt1) )  continue;
    iadr = (k-1)*3 + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<3; i++) {
      adj = adja[i] / 3;
      voy = adja[i] % 3;
      if ( !adj ) {
        continue;
      }

      /* Check that curent elt is not adjacent to itself */
      if ( adj == k ) {
        if ( !mmgErr0 ) {
          mmgErr0 = 1;
          fprintf(stderr,"\n  ## Error: %s: 1. at least 1 wrong"
                  " adjacency %" MMG5_PRId " %" MMG5_PRId "\n",__func__,MMG2D_indElt(mesh,k),
                  MMG2D_indElt(mesh,adj));
          fprintf(stderr,"vertices of %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",MMG2D_indElt(mesh,k),
                  MMG2D_indPt(mesh,pt1->v[0]),MMG2D_indPt(mesh,pt1->v[1]),
                  MMG2D_indPt(mesh,pt1->v[2]));
          fprintf(stderr,"adj of %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",
                  k,adja[0]/3,adja[1]/3,adja[2]/3);
        }
        return 0;
      }

      /* Check existence of adja elt */
      pt2 = &mesh->tria[adj];
      if ( !MG_EOK(pt2) ) {
        if ( !mmgErr1 ) {
          mmgErr1 = 1;
          fprintf(stderr,"\n  ## Error: %s: 4. at least 1 invalid"
                  " adjacent %" MMG5_PRId " %" MMG5_PRId "\n",__func__,MMG2D_indElt(mesh,adj),
                  MMG2D_indElt(mesh,k));
          fprintf(stderr,"vertices of %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
                  MMG2D_indElt(mesh,k),MMG2D_indPt(mesh,pt1->v[0]),
                  MMG2D_indPt(mesh,pt1->v[1]),MMG2D_indPt(mesh,pt1->v[2]));
          fprintf(stderr,"vertices adj %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",
                  adj,pt2->v[0],pt2->v[1],pt2->v[2]);
          fprintf(stderr,"adj of %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",MMG2D_indElt(mesh,k),
                  MMG2D_indElt(mesh,adja[0]/3),MMG2D_indElt(mesh,adja[1]/3),
                  MMG2D_indElt(mesh,adja[2]/3));
        }
        return 0;
      }

      /* Check consistency of adjacency array */
      iadr  = (adj-1)*3 + 1;
      adja1 = &mesh->adja[iadr];
      adj1  = adja1[voy] / 3;
      voy1  = adja1[voy] % 3;
      if ( adj1 != k || voy1 != i ) {
        if ( !mmgErr2 ) {
          mmgErr2 = 1;
          fprintf(stderr,"\n  ## Error: %s: 2. at least 1 inconsistency in adja array"
                  " %" MMG5_PRId " %" MMG5_PRId "\n",__func__,MMG2D_indElt(mesh,k),MMG2D_indElt(mesh,adj1));
          fprintf(stderr,"vertices of %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",MMG2D_indElt(mesh,k),
                  MMG2D_indPt(mesh,pt1->v[0]),MMG2D_indPt(mesh,pt1->v[1]),
                  MMG2D_indPt(mesh,pt1->v[2]));
          fprintf(stderr,"adj(k) %" MMG5_PRId ": %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " \n",MMG2D_indElt(mesh,adj),
                  MMG2D_indPt(mesh,pt2->v[0]),MMG2D_indPt(mesh,pt2->v[1]),
                  MMG2D_indPt(mesh,pt2->v[2]));
          fprintf(stderr,"adj(%" MMG5_PRId "): %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
                  MMG2D_indElt(mesh,k),MMG2D_indElt(mesh,adja[0]/3),
                  MMG2D_indElt(mesh,adja[1]/3),MMG2D_indElt(mesh,adja[2]/3));
          fprintf(stderr,"adj(%" MMG5_PRId "): %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
                  MMG2D_indElt(mesh,adj),MMG2D_indElt(mesh,adja1[0]/3),
                  MMG2D_indElt(mesh,adja1[1]/3),MMG2D_indElt(mesh,adja1[2]/3),
                  MMG2D_indElt(mesh,adja1[3]/3));
        }
        return 0;
      }
    }
  }

  /** Check consistency between tags of edges and vertices */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      int i1 = MMG5_inxt2[i];
      int i2 = MMG5_iprv2[i];
      /* Check that a GEO edg has GEO or singular vertices (CRN or REQ) */
      if ( pt->tag[i] & MG_GEO ) {
        if ( !(mesh->point[pt->v[i1]].tag & MG_GEO) && !( MG_SIN(mesh->point[pt->v[i1]].tag) )) {
         if ( !mmgErr3 ) {
            mmgErr3 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %" MMG5_PRId ": edge %d, vertex %" MMG5_PRId ")\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i1]));
         }
          return 0;
        }
        if ( !(mesh->point[pt->v[i2]].tag & MG_GEO) && !( MG_SIN(mesh->point[pt->v[i2]].tag) )) {
          if ( !mmgErr4 ) {
            mmgErr4 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %" MMG5_PRId ": edge %d, vertex %" MMG5_PRId ")\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i2]));
          }
          return 0;
        }
      }

      /* Check that a REF edg has REF or singular vertices (CRN or REQ) */
      if ( pt->tag[i] & MG_REF ) {
        if ( !(mesh->point[pt->v[i1]].tag & MG_REF) && !( MG_SIN(mesh->point[pt->v[i1]].tag)) ) {
          if ( !mmgErr5 ) {
            mmgErr5 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %" MMG5_PRId ": edge %d, vertex %" MMG5_PRId ")\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i1]));
          }
          return 0;
        }
        if ( !(mesh->point[pt->v[i2]].tag & MG_REF) && !( MG_SIN(mesh->point[pt->v[i2]].tag)) ) {
          if ( !mmgErr6 ) {
            mmgErr6 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 tag inconsistency"
                    " (triangle %" MMG5_PRId ": edge %d, vertex %" MMG5_PRId ")\n",__func__,
                    MMG2D_indElt(mesh,k),i,MMG2D_indPt(mesh,pt->v[i2]));
          }
          return 0;
        }
      }

    }
  }

  /** Check consistency between edge tags and triangle refs */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      int i1 = MMG5_inxt2[i];
      int i2 = MMG5_iprv2[i];

      MMG5_int jel = adja[i] / 3;

      if ( base ) {
        if ( !jel ) {
          /* Check that a boundary edge has a bdy tag */
          if ( !(pt->tag[i] & MG_GEO ) ) {
            if ( !mmgErr7 ) {
              mmgErr7 = 1;
              fprintf(stderr,"\n  ## Error: %s: at least 1 wrong  edge tag"
                      " (edge %d in tria %" MMG5_PRId ".)\n",
                      __func__,i,MMG2D_indElt(mesh,k));
            }
            return 0;
          }
        }
        else if ( pt->tag[i] & MG_GEO ) {
          /* Check that an internal edge (even between 2 domains) is not MG_GEO
           * (used to mark opn bdy) */
          if ( !mmgErr8 ) {
            mmgErr8 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge tagged boundary"
                    " while it has a neighbour (%" MMG5_PRId " %" MMG5_PRId ").\n",__func__,
                    MMG2D_indPt(mesh,pt->v[i1]),MMG2D_indPt(mesh,pt->v[i2]));
          }
          return 0;
        }
      }

      if ( !mesh->info.opnbdy ) {
        /* Check that we don't have REF tag between tria of same refs */
        if ( pt->tag[i] & MG_REF ) {
          pt1 = &mesh->tria[jel];
          if ( pt->ref == pt1->ref ) {
            if ( !mmgErr9 ) {
              mmgErr9 = 1;
              fprintf(stderr,"\n  ## Error: %s: at least 1 edge tagged ref while"
                      " both corresponding triangles have same ref (%" MMG5_PRId
                      " %" MMG5_PRId ").\n",__func__,
                      MMG2D_indPt(mesh,pt->v[i1]), MMG2D_indPt(mesh,pt->v[i2]));
            }
            return 0;
          }
        }
      }
    }
  }

  /** Check consistency between REF, GEO and BDY tags between edges and points */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      if ( pt->tag[i] & MG_GEO || pt->tag[i] & MG_REF ) {
        int i1 = MMG5_inxt2[i];
        int i2 = MMG5_iprv2[i];

        /* Check that a GEO or REF edge is BDY too */
        if ( !(pt->tag[i] & MG_BDY) ) {
          if ( !mmgErr10 ) {
            mmgErr10 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge (%" MMG5_PRId " %" MMG5_PRId ") tagged"
                    " %d, but not MG_BDY\n",__func__,MMG2D_indPt(mesh,pt->v[i1]),
                    MMG2D_indPt(mesh,pt->v[i2]),pt->tag[i]);
          }
          return 0;
        }

        /* Check that a bdy edge has bdy vertices */
        MMG5_pPoint p1 = &mesh->point[pt->v[i1]];
        if ( !(p1->tag & MG_BDY) ) {
          if ( !mmgErr11 ) {
            mmgErr11 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge (%" MMG5_PRId " %" MMG5_PRId ") tagged %d,"
                    " with a point (%" MMG5_PRId ") not tagged BDY.\n",__func__,
                    MMG2D_indPt(mesh,pt->v[i1]),MMG2D_indPt(mesh,pt->v[i2]),
                    pt->tag[i],MMG2D_indPt(mesh,pt->v[i1]));
          }
          return 0;
        }

        MMG5_pPoint p2 = &mesh->point[pt->v[i2]];
        if ( !(p2->tag & MG_BDY) ) {
          if ( !mmgErr12 ) {
            mmgErr12 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 edge (%" MMG5_PRId " %" MMG5_PRId ") tagged %d,"
                    " with a point (%" MMG5_PRId ") not tagged BDY.\n",__func__,
                    MMG2D_indPt(mesh,pt->v[i1]),MMG2D_indPt(mesh,pt->v[i2]),
                    pt->tag[i],MMG2D_indPt(mesh,pt->v[i2]));
          }
          return 0;
        }
      }
    }
  }

  if ( !severe )  return 1;


  MMG5_SAFE_CALLOC(list,MMG2D_LMAX,MMG5_int,return 0);

  /** Checks on vertices */
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

      /* Check that a used elt has used vertices */
      if ( !MG_VOK(ppt) ) {
        if ( !mmgErr13 ) {
          mmgErr13 = 1;
          fprintf(stderr,"\n  ## Error: %s: at least 1 unused vertex %" MMG5_PRId
                  "  %" MMG5_PRId "\n",__func__,
                  MMG2D_indElt(mesh,k),MMG2D_indPt(mesh,ip));
          fprintf(stderr,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
                  MMG2D_indPt(mesh,pt1->v[0]),
                  MMG2D_indPt(mesh,pt1->v[1]),MMG2D_indPt(mesh,pt1->v[2]));
        }
        MMG5_SAFE_FREE(list);
        return 0;
      }

      /* Check the validity of the ball of point */
      lon = MMG2D_boulep(mesh,k,i,list);
      for (l=1; l<=lon; l++) {
        kk  = list[l] / 3;
        nk  = list[l] % 3;
        pt2 = &mesh->tria[kk];
        if ( pt2->v[nk] != ip ) {
          if ( !mmgErr14 ) {
            mmgErr14 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong ball %" MMG5_PRId
                    ", %" MMG5_PRId "\n",
                    __func__,MMG2D_indPt(mesh,ip),MMG2D_indPt(mesh,pt2->v[nk]));
          }
          MMG5_SAFE_FREE(list);
          return 0;
        }
      }
      if ( lon < 1 )  continue;

      /* Check that the number of valid triangles containing the vertex ip is
       * equal to the length of the ball of ip. */
      len = 0;
      for (kk=1; kk<=mesh->nt; kk++) {
        pt2 = &mesh->tria[kk];
        if ( !MG_EOK(pt2) )  continue;
        for (j=0; j<3; j++)
          if ( pt2->v[j] == ip ) {
            len++;
            break;
          }
      }
      if ( len != lon ) {
        if ( !mmgErr15 ) {
          mmgErr15 = 1;
          fprintf(stderr,"\n  ## Error: %s: at least 1 incorrect ball %"
                  MMG5_PRId ": %d %d\n",
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
