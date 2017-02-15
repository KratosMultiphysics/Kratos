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


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Check if all edges exist and if not force them.
 *
 */
int MMG2_bdryenforcement(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria      pt,pt1;
  MMG5_pEdge      ped;
  int       i,j,k,kk,nex,*list,kdep,lon,voy,iel,iare,ied;
  int       ia,ib,ilon,rnd,idep,*adja,ir,adj,list2[3],iter;//,iadr,ped0,ped1;
//  int       iadr2,*adja2,ndel;

  _MMG5_SAFE_CALLOC(list,MMG2D_LMAX,int);

  nex = 0;

  /** Liste edges to enforce */
  for(i=1 ; i<=mesh->na ; i++) {
    /* Edge to create */
    ped = &mesh->edge[i];
    if(!ped->a) continue;
    if(ped->base < 0) {
      nex++;
      continue;
    }
    /* Find first edge extremity */
    kdep = MMG2_findTria(mesh,ped->a);
    assert(kdep);

    if(mesh->tria[kdep].v[0]==ped->a)
      j=0;
    else if(mesh->tria[kdep].v[1]==ped->a)
      j=1;
    else
      j=2;

    lon = MMG2_boulep(mesh,kdep,j,list);
    if(lon>MMG2D_LMAX) {
      printf("  ## Error: too many triangles (%d) around the vertex %d\n",
             lon,ped->a);
      exit(EXIT_FAILURE);
    } else if(!lon) {
      printf(" ## Error: problem with point %d of triangle %d\n",
             mesh->tria[kdep].v[j],kdep);
      _MMG5_SAFE_FREE(list);
      return(0);
    }
    assert(lon);

    /* Find second edge extremity */
    for (j=1 ; j<=lon ; j++) {
      iel = list[j];
      pt = &mesh->tria[iel/3];
      voy = iel%3;
      if (pt->v[MMG2_iare[(voy+1)%3][0]]==ped->b) {
        ped->base = -1;
        nex++;
        break;
      }
      else if (pt->v[MMG2_iare[(voy+2)%3][1]]==ped->b) {
        ped->base = -1;
        nex++;
        break;
      }
    }


    if(j>lon) {
      /* Missing edge */
      if(mesh->info.imprim > 5) printf("  ** missing edge %d %d \n",
                                       ped->a,ped->b);
      ped->base = kdep;
    }
  }

  /** Enforce edges */
  if(nex!=mesh->na) {
    if(mesh->info.imprim > 4)
      printf(" ** number of missing edges : %d\n",mesh->na-nex);

    for(kk=1 ; kk<=mesh->na ; kk++) {
      ped = &mesh->edge[kk];
    if( !ped->a || ped->base < 0) continue;
      ia = ped->a;
      ib = ped->b;
      kdep = ped->base;
      if(mesh->info.ddebug)
        printf("\n  -- edge enforcement %d %d\n",ia,ib);

      if(!(lon=MMG2_locateEdge(mesh,ia,ib,&kdep,list))) {
        if(mesh->info.ddebug)
          printf("  ## Error: edge not found\n");
        _MMG5_SAFE_FREE(list);
        return(0);
      }
      if(!(lon<0 || lon==4)) {
        if(mesh->info.ddebug)
          printf("  ** Unable to force edge: %d %d -- %d\n",ia,ib,lon);
        exit(EXIT_FAILURE);
      }
      /*edge exist*/
      if(lon==4) {
        if(mesh->info.ddebug) printf("  ** Existing edge\n");
        //exit(EXIT_FAILURE);
      }
      if(lon>MMG2D_LMAX) {
        printf(" ## Error: too many triangles (%d)\n",lon);
        exit(EXIT_FAILURE);
      }
      if(lon<2) {
        if(mesh->info.ddebug) printf("  ** few edges... %d\n",lon);
        //exit(EXIT_FAILURE);
      }
      lon = -lon;
      ilon = lon;

      /*retournement d'arêtes aleatoirement dans la liste, tant que */
      srand(time(NULL));
      iter=0;
      while (ilon>0 && iter++<2*lon) {
        rnd = (rand()%lon);

        k = list[rnd]/3;
        if(mesh->info.ddebug) {
          printf("  ** Random edge swap\n");
          /* for(i=0 ; i<lon ; i++) { */
          /*   if((mesh->tria[list[i]/3]).base == mesh->base+1) */
          /*     printf("%de tr : %d base %d == %d\n",i+1,list[i]/3,mesh->tria[list[i]/3].base,mesh->base+1); */
          /*   else */
          /*     printf("pas base %de tr %d \n",i+1,list[i]/3); */

          /*   printf("vertex %d %d %d\n",mesh->tria[list[i]/3].v[0],mesh->tria[list[i]/3].v[1], */
          /*          mesh->tria[list[i]/3].v[2]); */
          /* } */
        }

        /*check k in Pipe*/
        i=0;
        k = list[rnd]/3;
        while(i++<lon) {
          pt = &mesh->tria[k];
          if(pt->base == mesh->base+1) break;
          k = list[(++rnd)%lon]/3;
        }
        assert(i<=lon);
        idep = list[rnd]%3;
        // if(mesh->info.ddebug) printf("i= %d < %d ? on demarre avec %d\n",i,lon+1,k);
        adja = &mesh->adja[3*(k-1)+1];
        for(i=0 ; i<3 ; i++) {
          ir = (idep+i)%3;
          /*check adj in Pipe*/
          adj = adja[ir]/3;
          voy = adja[ir]%3;
          pt1 = &mesh->tria[adj];
          if (pt1->base != (mesh->base+1)) {
            continue;
          }
          /************************/
          /********swap***********/
          /************************/
          if(!MMG2_swapar(mesh,sol,k,ir,1e+4,list2)) {
            if(mesh->info.ddebug) printf("  ## Warning: unable to swap\n");
            continue;
          }
          if(mesh->info.ddebug) printf("  ** Successful swap\n");
          /*new tr intersecté par ia-ib ??*/
          for(ied=1 ; ied<3 ; ied++) {
            iare = MMG2_cutEdgeTriangle(mesh,list2[ied],ia,ib);
            if(!iare) { /*tr not in pipe*/
              ilon--;
              if(mesh->info.ddebug)
                printf("  ## Warning: tr %d not intersected ==> %d\n",list2[ied],ilon);
              mesh->tria[list2[ied]].base = mesh->base;
            } else if(iare < 0) {
              mesh->tria[list2[ied]].base = mesh->base;
              ilon -= 2;
            } else {
              if(mesh->info.ddebug) printf("  ** tr intersected %d \n",list2[ied]);
              mesh->tria[list2[ied]].base = mesh->base+1;
            }
          }
          break;
        }
     }
  }/*end k --> mesh->na*/
  }

  _MMG5_SAFE_FREE(list);


/*   //check if there are more bdry edges.. and delete tr      */
/* #warning a optimiser en mettant un pile et un while */
/*   //printf("tr 1 : %d %d %d\n",mesh->tria[1].v[0],mesh->tria[1].v[1],mesh->tria[1].v[2]);  */
/*   do { */
/*     ndel=0; */
/*     for(k=1 ; k<=mesh->nt ; k++) { */
/*       pt = &mesh->tria[k]; */
/*       if(!pt->v[0]) continue;  */
/*       iadr = 3*(k-1) + 1;                       */
/*       adja = &mesh->adja[iadr];  */

/*       for(i=0 ; i<3 ; i++)  { */
/*  if(adja[i]) continue; */
/*         ped0 = pt->v[MMG2_iare[i][0]]; */
/*         ped1 = pt->v[MMG2_iare[i][1]]; */
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
  /*       fprintf(stdout,"WARNING BDRY EDGES MISSING : DO YOU HAVE DUPLICATED VERTEX ? %d %d %d\n", */
/*        pt->v[0],pt->v[1],pt->v[2]); */
/*       for(i=0 ; i<3 ; i++) {         */
/*  if(!adja[i]) continue; */
/*  iadr2 = 3*(adja[i]/3-1) + 1; */
/*  adja2 = &mesh->adja[iadr2]; */
/*  adja2[adja[i]%3] = 0; */
/*       } */
/*       MMG2_delElt(mesh,k);  */
/*       ndel++; */

/*     }           */
/*   } while(ndel); */

  return(1);
}
