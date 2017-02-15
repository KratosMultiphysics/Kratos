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


/* read mesh data */
int MMG2_evalgeom(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppt,ppa,ppb;
  double    capx,capy,cbpx,cbpy,alpha,rbound;
  int       *list,k,j,lon,iadr,*adja,nbdry,ibdry[2],ip,i,iel;
  int       ref,nc;

  nc = 0;
  rbound = mesh->info.dhd*M_PI/180.;
  /*corners detection*/
  _MMG5_SAFE_MALLOC(list,MMG2D_LMAX,int);

  ++mesh->base;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if(!M_EOK(pt)) continue;
    for(j=0 ; j<3 ; j++) {
      ppt = &mesh->point[pt->v[j]];
      if( (!(ppt->tag & M_BDRY)) && !(ppt->tag & M_SD)) continue;
      if(ppt->tagdel == mesh->base) continue;

      lon = MMG2_boulep(mesh,k,j,list);
      assert(lon);
      /*bdry triangles*/
      nbdry = 0;
      ref = mesh->tria[list[1]/3].ref;
      for(i=1 ; i<=lon ; i++) {
        iel = list[i]/3;
        ip  = list[i]%3;
        pt1 = &mesh->tria[iel];
        if(pt1->ref!=ref) continue;
        iadr = 3*(iel-1) + 1;
        adja = &mesh->adja[iadr];
        if(!adja[MMG2_iopp[ip][0]] || (mesh->tria[adja[MMG2_iopp[ip][0]]/3].ref != ref)) {

          if(nbdry>=2) fprintf(stdout,"NON MANIFOLD MESH\n");
          if(MMG2_iare[MMG2_iopp[ip][0]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][0]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][0]][0]];
          }
          //printf("1) tr %d : %d -- edge %d = %d %d\n",k,iel,MMG2_iopp[ip][0],MMG2_iare[MMG2_iopp[ip][0]][0],
          //                  MMG2_iare[MMG2_iopp[ip][0]][1]);
        }
        if(!adja[MMG2_iopp[ip][1]] || (mesh->tria[adja[MMG2_iopp[ip][1]]/3].ref != ref)) {
          if(nbdry>=2) fprintf(stdout,"NON MANIFOLD MESH\n");
          if(MMG2_iare[MMG2_iopp[ip][1]][0]==ip)
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][1]];
          else {
            assert(MMG2_iare[MMG2_iopp[ip][1]][1]==ip) ;
            ibdry[nbdry++] = pt1->v[MMG2_iare[MMG2_iopp[ip][1]][0]];
          }
          //printf("2) tr %d : %d -- edge %d = %d %d\n",k,iel);
        }
      }
      if(nbdry!=2) {
        fprintf(stdout,"NON MANIFOLD DOMAIN (NO CORNERS DETECTION) %d -- vertex %d\n",nbdry,k);
        continue;
      }
      //calcul de l'angle forme par les 3 points
      ppa  = &mesh->point[ibdry[0]];
      ppb  = &mesh->point[ibdry[1]];
      capx = ppt->c[0] - ppa->c[0];
      capy = ppt->c[1] - ppa->c[1];
      cbpx = ppt->c[0] - ppb->c[0];
      cbpy = ppt->c[1] - ppb->c[1];
      alpha = capx*cbpx + capy*cbpy;
      alpha /= sqrt(capx*capx+capy*capy)*sqrt(cbpx*cbpx+cbpy*cbpy);
      alpha = acos(alpha);
      //printf("point %d : %e (= %e)-- %e %e\n",pt->v[j],alpha,alpha*180./M_PI,capx,capy);
      if(alpha < rbound && (!(ppt->tag & M_CORNER)) ) {
        ++nc;
        ppt->tag |= M_CORNER;
      }

      ppt->tagdel++;
    }
  }

  _MMG5_SAFE_FREE(list);

  if ( abs(mesh->info.imprim) > 3 && nc )
    fprintf(stdout,"     %d corners detected\n",nc);
  return(1);
}
