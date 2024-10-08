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


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Check if all edges exist in the mesh and if not force them.
 *
 */
int MMG2D_bdryenforcement(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      pt,pt1;
  MMG5_pEdge      ped;
  MMG5_pPoint     ppt;
  MMG5_int        list[MMG5_TRIA_LMAX],list2[3];
  MMG5_int        k,l,kk,nex,kdep,lon,iel;
  int             iare,ied;
  MMG5_int        ia,ib,ilon,rnd,idep,*adja,ir,adj;
  int8_t          i,i1,i2,j;
//  int       iadr2,*adja2,ndel,iadr,ped0,ped1;
  static int8_t   mmgWarn0=0,mmgWarn1=0,mmgWarn2=0,mmgWarn3=0;
  static int8_t   mmgWarn4=0,mmgWarn5=0,mmgWarn6=0,mmgWarn7=0;
  static int8_t   mmgWarn8=0;

  nex = 0;

  /* Associate seed to each point in the mesh */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) ppt->s = 0;
  }

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++)
      mesh->point[pt->v[i]].s = k;
  }

  /* Check for the existing edges in the mesh */
  /* No command for identifying whether an edge is valid or not? */
  for(k=1; k<=mesh->na; k++) {
    ped = &mesh->edge[k];
    if ( !ped->a ) continue;
    if ( ped->base < 0 ) {
      nex++;
      continue;
    }

    ppt = &mesh->point[ped->a];
    kdep = ppt->s;
    pt = &mesh->tria[kdep];

    if ( pt->v[0] == ped->a )
      j=0;
    else if ( pt->v[1] == ped->a )
      j=1;
    else
      j=2;

    int8_t dummy;
    lon = MMG5_boulet(mesh,kdep,j,list,0,&dummy);

    if ( lon <= 0 ) {
      if ( !mmgWarn0 ) {
        mmgWarn0=1;
        fprintf(stderr,"\n  ## Error: %s: at least 1 wrong ball "
                "(point %" MMG5_PRId " of triangle %" MMG5_PRId ").\n",__func__,
                MMG2D_indPt(mesh, mesh->tria[kdep].v[j]),MMG2D_indElt(mesh,kdep));
      }
      return 0;
    }

    for (l=0; l<lon; l++) {
      iel = list[l] / 3;
      pt  = &mesh->tria[iel];
      i = list[l] % 3;
      i1 = MMG5_inxt2[i];
      i2 = MMG5_iprv2[i];

      if ( pt->v[i1] == ped->b ) {
        ped->base = -1;
        nex++;
        break;
      }
      else if ( pt->v[i2] == ped->b ) {
        ped->base = -1;
        nex++;
        break;
      }
    }

    if ( l >= lon ) {
      if ( (mesh->info.imprim > 5) && (!mmgWarn1) ) {
        mmgWarn1 = 1;
        fprintf(stderr,"\n  ## Warning: %s: at least 1 missing edge (%" MMG5_PRId " %" MMG5_PRId ").\n",
                __func__,MMG2D_indPt(mesh,ped->a),MMG2D_indPt(mesh,ped->b));
      }
      ped->base = kdep;
    }
  }

  /** Now treat the missing edges */
  if ( nex != mesh->na ) {
    if(mesh->info.imprim > 5)
      printf(" ** number of missing edges : %" MMG5_PRId "\n",mesh->na-nex);

    for (kk=1; kk<=mesh->na; kk++) {
      ped = &mesh->edge[kk];
      if ( !ped->a || ped->base < 0 ) continue;
      ia = ped->a;
      ib = ped->b;
      kdep = ped->base;

      if(mesh->info.ddebug)
        printf("\n  -- edge enforcement %" MMG5_PRId " %" MMG5_PRId "\n",ia,ib);

      /* List of the triangles intersected by the edge */
      if ( !(lon=MMG2D_locateEdge(mesh,ia,ib,&kdep,list)) ) {
        if ( mesh->info.ddebug && (!mmgWarn2) ) {
          fprintf(stderr,"\n  ## Error: %s: at least 1 edge not found.\n",
                  __func__);
          mmgWarn2=1;
        }
        return 0;
      }

      /* Failure */
      if ( !( lon < 0 || lon == 4 ) ) {
        if ( mesh->info.ddebug && (!mmgWarn3) ) {
          mmgWarn3=1;
          fprintf(stderr,"\n ## Error: %s: Unable to force at least"
                  " 1 edge (%" MMG5_PRId " %" MMG5_PRId " -- %" MMG5_PRId ").\n",__func__,MMG2D_indPt(mesh,ia),
                  MMG2D_indPt(mesh,ib),lon);
        }
        return 0;
      }

      /* Insertion of the edge */
      lon = -lon;

      /* Considered edge actually exists */
      if ( lon == 4 ) {
        if ( mesh->info.ddebug && (!mmgWarn4) ) {
          mmgWarn4=1;
          fprintf(stderr,"\n  ## Warning: %s: existing edge.\n",__func__);
        }
      }
      if( (!mmgWarn5) && (lon>MMG5_TRIA_LMAX) ) {
        mmgWarn5 = 1;
        fprintf(stderr,"\n  ## Error: %s: at least 1 edge intersecting too many"
               " triangles (%" MMG5_PRId ")\n",__func__,lon);
      }
      if ( lon<2 ) {
        if ( mesh->info.ddebug && (!mmgWarn6) ) {
          mmgWarn6 = 1;
          fprintf(stderr,"\n  ## Warning: %s: few edges... %" MMG5_PRId "\n",__func__,lon);
        }
      }

      /* Randomly swap edges in the list, while... */
      ilon = lon;
      srand(time(NULL));

      while ( ilon > 0 ) {
        rnd = ( rand() % lon );
        k   = list[rnd] / 3;

        /* Check for triangle k in the pipe until one triangle with base+1 is met (indicating that it is
           crossed by the considered edge) */
        l = 0;
        while ( l++ < lon ) {
          pt = &mesh->tria[k];
          if ( pt->base == mesh->base+1 ) break;
          k = list[(++rnd)%lon] / 3;
        }

        assert ( l <= lon );
        idep = list[rnd] % 3;
        // if(mesh->info.ddebug) printf("i= %" MMG5_PRId " < %" MMG5_PRId " ? on demarre avec %" MMG5_PRId "\n",i,lon+1,k);
        adja = &mesh->adja[3*(k-1)+1];

        for (i=0; i<3; i++) {
          ir = (idep+i) % 3;

          /* Check the adjacent triangle in the pipe */
          adj = adja[ir] / 3;
          pt1 = &mesh->tria[adj];
          if ( pt1->base != (mesh->base+1) ) continue;

          /* Swap edge ir in triangle k, corresponding to a situation where both triangles are to base+1 */
          if ( !MMG2D_swapdelone(mesh,sol,k,ir,1e+4,list2) ) {
            if ( mesh->info.ddebug && (!mmgWarn7) ) {
              mmgWarn7 = 1;
              fprintf(stderr,"\n  ## Warning: %s: unable to swap at least 1"
                      " edge.\n",__func__);
            }
            continue;
          }

          /* Is new triangle intersected by ia-ib ?? */
          for (ied=1; ied<3; ied++) {
            iare = MMG2D_cutEdgeTriangle(mesh,list2[ied],ia,ib);
            if ( !iare ) {
              /* tr not in pipe */
              ilon--;
              if ( mesh->info.ddebug && (!mmgWarn8) ) {
                mmgWarn8 = 1;
                fprintf(stderr,"\n  ## Warning: %s: at least 1 triangle (%" MMG5_PRId ")"
                        " not intersected ==> %" MMG5_PRId "\n",__func__,
                        MMG2D_indElt(mesh,list2[ied]),ilon);
              }
              mesh->tria[list2[ied]].base = mesh->base;
            }
            else if ( iare < 0 ) {
              /* ia-ib is one edge of k */
              mesh->tria[list2[ied]].base = mesh->base;
              ilon -= 2;
            }
            else {
              if ( mesh->info.ddebug )
                printf("  ** tr intersected %" MMG5_PRId " \n",list2[ied]);
              mesh->tria[list2[ied]].base = mesh->base+1;
            }
          }
          break;
        }
     }
    }
  }

  /* Reset ->s field of vertices */
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) ppt->s = 0;
  }

/*   //check if there are more bdry edges.. and delete tr      */
/* #warning a optimiser en mettant un pile et un while */
/*   //printf("tr 1 : %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",mesh->tria[1].v[0],mesh->tria[1].v[1],mesh->tria[1].v[2]);  */
/*   do { */
/*     ndel=0; */
/*     for(k=1 ; k<=mesh->nt ; k++) { */
/*       pt = &mesh->tria[k]; */
/*       if(!pt->v[0]) continue;  */
/*       iadr = 3*(k-1) + 1;                       */
/*       adja = &mesh->adja[iadr];  */

/*       for(i=0 ; i<3 ; i++)  { */
/*  if(adja[i]) continue; */
/*         ped0 = pt->v[MMG2D_iare[i][0]]; */
/*         ped1 = pt->v[MMG2D_iare[i][1]]; */
/*         for(j=1 ; j<=mesh->na ; j++) { */
/*           ped = &mesh->edge[j]; */
/*           if((ped->a == ped0 && ped->b==ped1) || (ped->b == ped0 && ped->a==ped1)) break;   */
/*         } */
/*         if(j<=mesh->na)  */
/*           continue; */
/*  else  */
/*    break; */
/*       } */
/*       if(i==3) continue;  */
  /*       fprintf(stdout,"WARNING BDRY EDGES MISSING : DO YOU HAVE DUPLICATED VERTEX ? %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n", */
/*        pt->v[0],pt->v[1],pt->v[2]); */
/*       for(i=0 ; i<3 ; i++) {         */
/*  if(!adja[i]) continue; */
/*  iadr2 = 3*(adja[i]/3-1) + 1; */
/*  adja2 = &mesh->adja[iadr2]; */
/*  adja2[adja[i]%3] = 0; */
/*       } */
/*       MMG2D_delElt(mesh,k);  */
/*       ndel++; */

/*     }           */
/*   } while(ndel); */

  return 1;
}
