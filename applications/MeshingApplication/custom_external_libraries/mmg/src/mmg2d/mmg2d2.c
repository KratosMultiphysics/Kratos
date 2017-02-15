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
/**
 * \file mmg2d/mmg2d2.c
 * \brief Mesh generation functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "mmg2d.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \return 1 if success, 0 if fail.
 *
 * Remove the bounding box triangles.
 *
 */
int MMG2_removeBBtriangles(MMG5_pMesh mesh) {
  MMG5_pTria   pt;
  int     ip1,ip2,ip3,ip4,k,iadr,*adja,iadr2,*adja2;
  int     i,nd;

  /*BB vertex*/
  ip1=(mesh->np-3);
  ip2=(mesh->np-2);
  ip3=(mesh->np-1);
  ip4=(mesh->np);

  nd = 0;
  for(k=1 ; k<=mesh->nt ; k++) {
    pt  = &mesh->tria[k];
    if(!pt->v[0]) continue;
    if(pt->base<0) {
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for(i=0 ; i<3 ; i++) {
        if(!adja[i]) continue;
        iadr2 = 3*(adja[i]/3-1) + 1;
        adja2 = &mesh->adja[iadr2];
        adja2[adja[i]%3] = 0;
      }
      _MMG2D_delElt(mesh,k);
      continue;
    } else if(!pt->base) {
      printf("  ## Warning: undetermined triangle %d : %d %d %d\n",
             k,pt->v[0],pt->v[1],pt->v[2]);
      nd++;
    }
  }

  if(!nd) {
    _MMG2D_delPt(mesh,ip1);
    _MMG2D_delPt(mesh,ip2);
    _MMG2D_delPt(mesh,ip3);
    _MMG2D_delPt(mesh,ip4);
  } else {
    fprintf(stdout,"  ## Error: procedure failed : %d indetermined triangles\n",nd);
    return(0);
  }
  return(1);
}

/*tag des triangles : in = base ; out = -base ; indetermine = 0*/
int MMG2_settagtriangles(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria  pt;
  int    base,nd,iter,maxiter,k;
  int    ip1,ip2,ip3,ip4;

  /*BB vertex*/
  ip1=(mesh->np-3);
  ip2=(mesh->np-2);
  ip3=(mesh->np-1);
  ip4=(mesh->np);

  base = ++mesh->base;
  iter    = 0;
  maxiter = 3;
  do {
    nd = 0;
    for(k=1 ; k<=mesh->nt ; k++) {
      pt = &mesh->tria[k];
      if ( !M_EOK(pt) )  continue;
      if(!MMG2_findtrianglestate(mesh,k,ip1,ip2,ip3,ip4,base)) nd++ ;
    }
    if(mesh->info.ddebug) printf(" ** how many undetermined triangles ? %d\n",nd);
  } while (nd && ++iter<maxiter);

  return(1);
}

/**
 * \param mesh poitner toward the mesh structure
 * \param sol pointer toward the solution structure
 * \return  0 if fail.
 *
 * Insertion of the list of points inside the mesh
 *
 */
int MMG2_insertpointdelone(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pPoint  ppt;
  //  MMG5_pTria   pt;
  int     list[MMG2_LONMAX],lon;
  int     k,i;
  // int     kk, iadr,*adja,

  for(k=1 ; k<=mesh->np - 4 ; k++) {

    ppt = &mesh->point[k];

    /*recherche du triangle contenant le point : lel*/
    list[0] = MMG2_findTria(mesh,k);
    if(!list[0]) {
      if ( mesh->info.ddebug )
        printf(" ** exhaustive search of point location.\n");

      for(i = 1 ; i<= mesh->nt ; i++) {
        list[0] = MMG2_isInTriangle(mesh,i,&ppt->c[0]);
        if(list[0]) break;
      }
      if(i>mesh->nt) {
        fprintf(stdout,"  ## Error: no triangle found for vertex %d\n",k);
        return(0);
      }
    }

    lon = _MMG2_cavity(mesh,sol,k,list);

    if ( lon < 1 ) {
      fprintf(stdout,"  ## Error: unable to insert vertex %d\n",k);
      return(0);

    } else {
      _MMG2_delone(mesh,sol,k,list,lon);
    }
  }

  return(1);
}

/*insertion of the list of points inside the mesh*/
/*return 0 if pbs occur*/
int MMG2_insertpoint(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria   pt,pt1;
  MMG5_pPoint  ppt;
  double  declic;
  int     k,nsiter,lel,mel,nel,ia,ib,ic,aext0,aext1,aext2;
  int     iadr,*adja,nflat,ie1,ie2,ie3,voy,text,atext1,atext2;
  double  EPSDD = 1e+6;

  for(k=1 ; k<=mesh->np - 4 ; k++) {
    if(1){
      if(mesh->info.ddebug) printf("------------------------------------try swap to increase mesh quality\n");
      declic = 1.1 / ALPHA;
      nsiter = 1;
      while (nsiter) {
        nsiter = MMG2_cendel(mesh,sol,declic,-1);
        if ( nsiter && mesh->info.imprim < 0)
          fprintf(stdout,"     %7d SWAPPED\n",nsiter);
      }
    }
    ppt = &mesh->point[k];
    if(ppt->tmp==1) continue;
    /*recherche du triangle contenant le point : lel*/
    lel = MMG2_findTria(mesh,k);
    assert(lel);

    nflat = 0; // to avoid bad triangle  1, 2, 4, 3, 5, 7
    pt  = &mesh->tria[lel];
    iadr = 3*(lel-1) + 1;
    adja = &mesh->adja[iadr];
    ia  = pt->v[0];
    ib  = pt->v[1];
    ic  = pt->v[2];
    aext0 = adja[0];
    aext1 = adja[1];
    aext2 = adja[2];

    /*creation de trois triangles*/
    mel = _MMG2D_newElt(mesh);
    if ( !mel ) {
      _MMG5_TRIA_REALLOC(mesh,mel,mesh->gap,
                         printf("  ## Error: unable to allocate a new element.\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      pt  = &mesh->tria[lel];
    }
    nel = _MMG2D_newElt(mesh);
    if ( !nel ) {
      _MMG5_TRIA_REALLOC(mesh,mel,mesh->gap,
                         printf("  ## Error: unable to allocate a new element.\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         printf("  Exit program.\n");
                         exit(EXIT_FAILURE));
      pt  = &mesh->tria[lel];
    }
    pt->v[2] = k;
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
    if(pt->qual > EPSDD) {
      nflat = 1;
    } else {
      adja[0] = 3*nel + 1;
      adja[1] = 3*mel + 0;
      if(aext2)
        (&mesh->adja[3*(aext2/3-1) + 1])[aext2%3] = 3*lel + 2;
      //printf("adj of %d : %d %d %d -- %d(%d) : %d\n",lel,nel,mel,aext2/3,aext2/3,aext2%2,lel);
    }
    pt  = &mesh->tria[nel];
    pt->base = mesh->base;
    iadr = 3*(nel-1) + 1;
    adja = &mesh->adja[iadr];
    pt->v[0] = ib;
    pt->v[1] = ic;
    pt->v[2] = k;
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
    if(pt->qual > EPSDD) {
      nflat += 2;
    } else {
      adja[0] = 3*mel + 1;
      adja[1] = 3*lel + 0;
      adja[2] = aext0;
      if(aext0)
        (&mesh->adja[3*(aext0/3-1) + 1])[aext0%3] = 3*nel + 2;
      //printf("adj of %d : %d %d %d -- %d(%d) : %d\n",nel,mel,lel,aext0/3,aext0/3,aext0%3,nel);
    }
    pt  = &mesh->tria[mel];
    pt->base = mesh->base;
    iadr = 3*(mel-1) + 1;
    adja = &mesh->adja[iadr];
    pt->v[0] = ic;
    pt->v[1] = ia;
    pt->v[2] = k;
    pt->qual = MMG2_caltri_in(mesh,sol,pt);
    if(pt->qual > EPSDD) {
      nflat += 4;
      /*on coupe le tr aext1 en 2*/
      text = aext1/3;
      voy  = aext1%3;
      pt1  = &mesh->tria[text];
      ie1  = pt1->v[voy];
      ie2  = pt1->v[MMG2_inxt[voy]];
      ie3  = pt1->v[MMG2_inxt[voy+1]];
      pt1->v[0] = ie1;
      pt1->v[1] = ie2;
      pt1->v[2] = k;
      pt1->qual = MMG2_caltri_in(mesh,sol,pt1);
      pt->v[0] = ie3;
      pt->v[1] = ie1;
      pt->v[2] = k;
      pt->qual = MMG2_caltri_in(mesh,sol,pt);
      /*maj des adj*/
      /*text*/
      iadr = 3*(text-1) + 1;
      adja = &mesh->adja[iadr];
      atext1 = adja[MMG2_inxt[voy]];
      atext2 = adja[MMG2_inxt[voy+1]];
      adja[0] = 3*lel + 1;
      adja[1] = 3*mel + 0;
      adja[2] = atext2;
      if(atext2)
        (&mesh->adja[3*(atext2/3-1) + 1])[atext2%3] = 3*text + 2;
      //if(ddebug) printf("adj of %d : %d %d %d \n",text,adja[0]/3,adja[1]/3,adja[2]/3);
      /*mel*/
      iadr = 3*(mel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[0] = 3*text + 1;
      adja[1] = 3*nel + 0;
      adja[2] = atext1;
      if(atext1)
        (&mesh->adja[3*(atext1/3-1) + 1])[atext1%3] = 3*mel + 2;
      if(nflat==6 /*|| (*adj)==7*/) {//CECILE 14/04/14:*adj ne veut rien dire!!!
        adja[1] = 3*(aext0/3);
        //pas de suite sinon on perd de l'info (&mesh->adja[3*(aext0/3-1) + 1])[0] = 3*mel + 1;

      } else {
        /*nel*/
        iadr = 3*(nel-1) + 1;
        adja = &mesh->adja[iadr];
        adja[0] = 3*mel + 1;
      }

      /*lel*/
      iadr = 3*(lel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[1] = 3*text + 0;

    } else {
      adja[0] = 3*lel + 1;
      adja[1] = 3*nel + 0;
      adja[2] = aext1;
      if(aext1)
        (&mesh->adja[3*(aext1/3-1) + 1])[aext1%3] = 3*mel + 2;
      //printf("adj of %d : %d %d %d -- %d(%d) : %d\n",mel,lel,nel,aext1/3,aext1/3,aext1%3,mel);

    }
    if(nflat==1 || nflat==3 || nflat==7 || nflat==5){
      /*on coupe le tr aext2 en 2*/
      text = aext2/3;
      voy  = aext2%3;
      pt1  = &mesh->tria[text];
      ie1  = pt1->v[voy];
      ie2  = pt1->v[MMG2_inxt[voy]];
      ie3  = pt1->v[MMG2_inxt[voy+1]];
      pt1->v[0] = ie1;
      pt1->v[1] = ie2;
      pt1->v[2] = k;
      pt1->qual = MMG2_caltri_in(mesh,sol,pt1);
      pt = &mesh->tria[lel];
      pt->v[0] = ie3;
      pt->v[1] = ie1;
      pt->v[2] = k;
      pt->qual = MMG2_caltri_in(mesh,sol,pt);
      /*maj des adj*/
      /*text*/
      iadr = 3*(text-1) + 1;
      adja = &mesh->adja[iadr];
      atext1 = adja[MMG2_inxt[voy]];
      atext2 = adja[MMG2_inxt[voy+1]];
      adja[0] = 3*nel + 1;
      adja[1] = 3*lel + 0;
      adja[2] = atext2;
      if(atext2)
        (&mesh->adja[3*(atext2/3-1) + 1])[atext2%3] = 3*text + 2;

      /*lel*/
      iadr = 3*(lel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[0] = 3*text + 1;
      adja[1] = 3*mel + 0;
      adja[2] = atext1;
      if(atext1)
        (&mesh->adja[3*(atext1/3-1) + 1])[atext1%3] = 3*lel + 2;
      if(nflat==5 || nflat==7) {
        adja[1] = 3*(aext1/3);
        if(aext1)
          (&mesh->adja[3*(aext1/3-1) + 1])[0] = 3*lel + 1;
      }

      /*nel*/
      iadr = 3*(nel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[1] = 3*text + 0;

      /*mel*/
      iadr = 3*(mel-1) + 1;
      adja = &mesh->adja[iadr];
      if(!(nflat==5 || nflat==7))
        adja[0] = 3*lel + 1;

    }
    if(nflat==2 || nflat==3 || nflat==7 || nflat == 6){
      /*on coupe le tr aext0 en 2*/
      text = aext0/3;
      voy  = aext0%3;
      pt1  = &mesh->tria[text];
      ie1  = pt1->v[voy];
      ie2  = pt1->v[MMG2_inxt[voy]];
      ie3  = pt1->v[MMG2_inxt[voy+1]];
      pt1->v[0] = ie1;
      pt1->v[1] = ie2;
      pt1->v[2] = k;
      pt1->qual = MMG2_caltri_in(mesh,sol,pt1);
      pt = &mesh->tria[nel];
      pt->v[0] = ie3;
      pt->v[1] = ie1;
      pt->v[2] = k;
      pt->qual = MMG2_caltri_in(mesh,sol,pt);
      /*maj des adj*/
      /*text*/
      iadr = 3*(text-1) + 1;
      adja = &mesh->adja[iadr];
      atext1 = adja[MMG2_inxt[voy]];
      atext2 = adja[MMG2_inxt[voy+1]];
      adja[0] = 3*mel + 1;
      adja[1] = 3*nel + 0;
      adja[2] = atext2;
      if(atext2)
        (&mesh->adja[3*(atext2/3-1) + 1])[atext2%3] = 3*text + 2;
      /*nel*/
      iadr = 3*(nel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[0] = 3*text + 1;
      adja[1] = 3*lel + 0;
      adja[2] = atext1;
      if(atext1)
        (&mesh->adja[3*(atext1/3-1) + 1])[atext1%3] = 3*nel + 2;
      if (nflat==3 || nflat==7) {
        adja[1] = 3*(aext2/3);
        //(&mesh->adja[3*(aext2/3-1) + 1])[aext2%3] = 3*nel + 1;
      }

      /*mel*/
      iadr = 3*(mel-1) + 1;
      adja = &mesh->adja[iadr];
      adja[1] = 3*text + 0;
    }
    if(mesh->info.ddebug) {
      if ( !_MMG5_chkmsh(mesh,0,0) ) exit(EXIT_FAILURE);
    }

  }
  return(1);
}

/*put different ref on different SD*/
int MMG2_markSD(MMG5_pMesh mesh) {
  MMG5_pTria   pt,pt1;
  MMG5_pEdge   ped;
  MMG5_pPoint  ppt;
  int     k,i,j,iadr,*adja,ped0,ped1,*list,ipil,ncurc,nref;
  int     kinit,nt,nsd,ip1,ip2,ip3,ip4,ned,iel,voy;

  if ( !MMG2_hashel(mesh) )  return(0);
  for(k=1 ; k<=mesh->nt ; k++) mesh->tria[k].flag = mesh->mark;
  _MMG5_SAFE_CALLOC(list,mesh->nt,int);
  kinit = 0;
  nref  = 0;
  ip1=(mesh->np);
  for(k=1 ; k<=mesh->nt ; k++) {
    if ( !mesh->tria[k].v[0] ) continue;
    pt = &mesh->tria[k];
    pt->flag = mesh->mark;
    pt->ref  = 0;
    list[k-1] = 0;
    if((!kinit) && (pt->v[0]==ip1 || pt->v[1]==ip1 || pt->v[2]==ip1)) kinit = k;
  }
  do {
    nref++;
    list[0] = kinit;
    ipil = 0;
    ncurc = 0;
    do {
      k = list[ipil];
      pt = &mesh->tria[k];
      pt->ref = nref;
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for(i=0 ; i<3 ; i++)  {
        pt1 = &mesh->tria[adja[i]/3];

        if(pt1->ref==nref) continue;
        ped0 = pt->v[MMG2_iare[i][0]];
        ped1 = pt->v[MMG2_iare[i][1]];
//#warning optimize th
        for(j=1 ; j<=mesh->na ; j++) {
          ped = &mesh->edge[j];
          if((ped->a == ped0 && ped->b==ped1) || (ped->b == ped0 && ped->a==ped1)) break;
        }
        if(j<=mesh->na) continue;

        pt1->ref = nref;
        if(adja[i])
          list[++ncurc] = adja[i]/3;
      }
      ++ipil ;
    } while (ipil <= ncurc) ;
    kinit = 0;
    for(k=1 ; k<=mesh->nt ; k++) {
      if ( !mesh->tria[k].v[0] ) continue;
      mesh->tria[k].flag = mesh->mark;
      list[k-1] = 0;
      if(!kinit && !(mesh->tria[k].ref)) kinit = k;
    }
  } while (kinit);

  fprintf(stdout," %8d SUB-DOMAINS\n",nref-1); //because we have BB triangles

  /*remove BB triangles*/
  /*BB vertex*/
  ip1=(mesh->np-3);
  ip2=(mesh->np-2);
  ip3=(mesh->np-1);
  ip4=(mesh->np);

  if(nref!=1) {
    nt = mesh->nt;
    for(k=1 ; k<=nt ; k++) {
      if ( !mesh->tria[k].v[0] ) continue;
      pt = &mesh->tria[k];
      for(i=0 ; i<3 ; i++)
        mesh->point[pt->v[i]].tag = MG_NUL;
      if(pt->ref != 1) continue;
      /*update adjacencies*/
      iadr = 3*(k-1)+1;
      adja = &mesh->adja[iadr];
      for(i=0 ; i<3 ; i++) {
        if(!adja[i]) continue;
        iel = adja[i]/3;
        voy = adja[i]%3;
        (&mesh->adja[3*(iel-1)+1])[voy] = 0;
      }
      _MMG2D_delElt(mesh,k);

    }
  } else { /*remove all the triangle containing one of the BB vertex*/
    nt = mesh->nt;
    for(k=1 ; k<=nt ; k++) {
      if ( !mesh->tria[k].v[0] ) continue;
      pt = &mesh->tria[k];
      for(i=0 ; i<3 ; i++)
        mesh->point[pt->v[i]].tag = MG_NUL;
      if(!(pt->v[0]==ip1 || pt->v[1]==ip1 || pt->v[2]==ip1 ||
           pt->v[0]==ip2 || pt->v[1]==ip2 || pt->v[2]==ip2 ||
           pt->v[0]==ip3 || pt->v[1]==ip3 || pt->v[2]==ip3 ||
           pt->v[0]==ip4 || pt->v[1]==ip4 || pt->v[2]==ip4 )) continue;
      /*update adjacencies*/
      iadr = 3*(k-1)+1;
      adja = &mesh->adja[iadr];
      for(i=0 ; i<3 ; i++) {
        if(!adja[i]) continue;
        iel = adja[i]/3;
        voy = adja[i]%3;
        (&mesh->adja[3*(iel-1)+1])[voy] = 0;
      }
      _MMG2D_delElt(mesh,k);
    }
  }


  _MMG2D_delPt(mesh,ip1);
  _MMG2D_delPt(mesh,ip2);
  _MMG2D_delPt(mesh,ip3);
  _MMG2D_delPt(mesh,ip4);

  if(mesh->info.renum) {
    nsd = mesh->info.renum;
    nt = mesh->nt;
    for(k=1 ; k<=nt ; k++) {
      if ( !mesh->tria[k].v[0] ) continue;
      pt = &mesh->tria[k];
      pt->ref--;
      if(mesh->tria[k].ref == nsd) continue;
      _MMG2D_delElt(mesh,k);
    }
  }

  /*remove vertex*/
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ pt->v[i] ];
      ppt->tag &= ~MG_NUL;
    }
  }
  /*remove edge*/
  ned = mesh->na;
  for (k=1; k<=ned; k++) {
    ped = &mesh->edge[k];
    if ( !ped->a )  continue;
    ppt = &mesh->point[ ped->a ];
    if(!M_VOK(ppt)) {
      _MMG5_delEdge(mesh,k);
      continue;
    }
    ppt = &mesh->point[ ped->b ];
    if(!M_VOK(ppt)) {
      _MMG5_delEdge(mesh,k);
      continue;
    }
  }

  _MMG5_SAFE_FREE(list);
  return(1);
}

/*cherche si le  tr est dedans ou dehors ou indetermine (donc edge a forcer ou cas non convexe)*/
/*return <0 value if triangle outside ; > 0 if triangle inside*/
int MMG2_findtrianglestate(MMG5_pMesh mesh,int k,int ip1,int ip2,int ip3,int ip4,int base) {
  MMG5_pTria     pt,pt1;
  MMG5_pEdge     ped;
  int       i,j,nb,ped0,ped1;
  int       numed[3],nbed,iadr,*adja,adj;

  pt = &mesh->tria[k];
  /*Count number of BB points*/
  nb = 0;
  for(i=0 ; i<3 ; i++) {
    if(pt->v[i]==ip1 || pt->v[i]==ip2 || pt->v[i]==ip3 || pt->v[i]==ip4) {
      nb++;
    }
  }
  /*triangle to be deleted*/
  if(nb==3 || nb==2) {
    pt->base = -base;
    pt->ref  = 3;
    return(-base);
  }

  /*we have to check the status of neighbour*/
  /*number of bdry edges*/
//#warning to optimize
  nbed = 0;
  for(i=0 ; i<3 ; i++) {
    ped0 = pt->v[MMG2_iare[i][0]];
    ped1 = pt->v[MMG2_iare[i][1]];
    numed[i] = 0;
    for(j=1 ; j<=mesh->na ; j++) {
      ped = &mesh->edge[j];
      if((ped->a == ped0 && ped->b==ped1) || (ped->b == ped0 && ped->a==ped1)){
        ped->base = -1;
        nbed++;
        numed[i] = j;
        break;
      }
    }
  }
  if(nbed<2) {
    if(nb) {
      pt->base = -base;
      pt->ref = 3;
      return(-base);
    } else {
//#warning check if it is true with 1 bdry edge
      pt->base = base;
      return(base);
    }
  } else if(nbed==2) {
    if(nb) {
      pt->base = -base;
      pt->ref = 3;
      return(-base);
    } else {
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for(i=0 ; i<3 ; i++) {
        adj = adja[i]/3;
        pt1 = &mesh->tria[adj];
        if(!pt1->base) continue;
        if(abs(pt1->base)<mesh->base) continue;
        if(numed[i]) { //wrong if the edge is an internal boundary
          pt->base = -pt1->base;
          return(pt->base);
        } else {
          pt->base = pt1->base;
          return(pt->base);
        }
      }
      if(i==3) {
        pt->base = 0;
        return(0);
      }
    }
  } else {
    assert(!nb);
    fprintf(stdout,"  ## Warning: connex component with only one triangle\n");
    pt->base = base;
    return(base);
  }


  /*le tr est indetermine*/
  pt->base = 0;

  return(0);
}

/**
 *
 * \return <0 value if triangle outside ; > 0 if triangle inside
 *
 * Seek if triangle is inside, outside or not determinated (thus, the edge must
 * be enforced or the case is non convexe.
 *
 */
int MMG2_findpos(MMG5_pMesh mesh,MMG5_pTria pt,int ip1,int ip2,int ip3,int ip4,int base) {
  MMG5_pEdge     ped;
  int            i,j,nb,ped0,ped1;
  unsigned int   MMG2_iare2[3][2] = { {0,1}, {0,2}, {1,2} };

  /*au moins 1 point de la BB ??*/
  nb = 0;
  for(i=0 ; i<3 ; i++) {
    if(pt->v[i]==ip1 || pt->v[i]==ip2 || pt->v[i]==ip3 || pt->v[i]==ip4) {
      nb++;
    }
  }
  if(nb==3 || nb==2) {
    pt->base = -base;
    pt->ref  = 3;
    return(-base);
  } else if(!nb){
    pt->base = base;
    //pt->ref  = base;
    return(base);
  }
  /*1 pt de la BB => on ne sait pas si le tr est dedans ou dehors*/
  /*contient un edge ?*/
  for(i=0 ; i<3 ; i++) {
    ped0 = pt->v[MMG2_iare2[i][0]];
    ped1 = pt->v[MMG2_iare2[i][1]];
    for(j=1 ; j<=mesh->na ; j++) {
      ped = &mesh->edge[j];
      if((ped->a == ped0 && ped->b==ped1) || (ped->b == ped0 && ped->a==ped1)){
        ped->base = -1;
        break;
      }
    }
    if(j<=mesh->na) break;
  }
  if(i<3) { /*le tr est dehors*/
    //printf("on a trouve edge : %d %d\n",ped0,ped1);
    pt->base = -base;
    pt->ref  = 3;
    return(base);
  }

  /*le tr est indetermine*/
  pt->base = 0;

  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure.
 * \return 0 if fail, 1 if success.
 *
 * Mesh triangulation.
 *
 **/
int MMG2_mmg2d2(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria     pt;
  //MMG5_pEdge     ped;
  MMG5_pPoint    ppt,ppt2;
  double    c[2],dd;
  int       j,k,kk,ip1,ip2,ip3,ip4,jel,kel,nt,iadr,*adja;
  int       *numper;

  mesh->base = 0;
  /*if existed triangle : delete them*/
  if(mesh->nt) {
    nt = mesh->nt;
    for(k=1 ; k<=nt ; k++) {
      _MMG2D_delElt(mesh,k);
      iadr = 3*(k-1) + 1;
      (&mesh->adja[iadr])[0] = 0;
      (&mesh->adja[iadr])[1] = 0;
      (&mesh->adja[iadr])[2] = 0;
    }
  }
  if(mesh->info.renum==-10) {
    /*traitement des points periodiques*/
    _MMG5_SAFE_CALLOC(numper,mesh->np+1,int);
    for(k=1 ; k<=mesh->np ; k++) {
      ppt = &mesh->point[k];
      for(kk=k ; kk<=mesh->np ; kk++) {
        if(k==kk) continue;
        ppt2 = &mesh->point[kk];
        dd = (ppt->c[0]-ppt2->c[0])*(ppt->c[0]-ppt2->c[0])+(ppt->c[1]-ppt2->c[1])*(ppt->c[1]-ppt2->c[1]);
        if(dd<1e-6) {
          //printf("point img %d %d\n",k,kk);
          ppt2->tmp = 1;
          if(!numper[k]) {
            numper[k] = kk;
          } else if(numper[k]!=kk){
            j = numper[k];
            // printf("j = %d %d\n",j,numper[j]) ;
            while(numper[j] && numper[j]!=kk) {
              j = numper[j];
              //printf("j = %d %d\n",j,numper[j]) ;
            }
            if(numper[j]!=kk) numper[j] = kk;
          }
        }
      }
    }
  }
  /*add bounding box vertex*/
  c[0] = -0.5;//mesh->info.min[0] - 1.;
  c[1] = -0.5;// mesh->info.min[1] - 1.;
  ip1 = _MMG2D_newPt(mesh,c,0);
  if ( !ip1 ) {
    /* reallocation of point table */
    _MMG2D_POINT_REALLOC(mesh,sol,ip1,mesh->gap,
                         printf("  ## Error: unable to allocate a new point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         return(-1)
                         ,c,0);
  }

  c[0] = -0.5;//mesh->info.min[0] - 1.;
  c[1] =  PRECI / mesh->info.delta *(mesh->info.max[1]-mesh->info.min[1])
    + 0.5;//mesh->info.max[1] + 1.;
  ip2 = _MMG2D_newPt(mesh,c,0);
  if ( !ip2 ) {
    /* reallocation of point table */
    _MMG2D_POINT_REALLOC(mesh,sol,ip2,mesh->gap,
                         printf("  ## Error: unable to allocate a new point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         return(-1)
                         ,c,0);
  }

  c[0] =  PRECI / mesh->info.delta *(mesh->info.max[0]-mesh->info.min[0])
    + 0.5;//mesh->info.max[0] + 1.;
  c[1] = -0.5;//mesh->info.min[1] - 1.;
  ip3 = _MMG2D_newPt(mesh,c,0);
  if ( !ip3 ) {
    /* reallocation of point table */
    _MMG2D_POINT_REALLOC(mesh,sol,ip3,mesh->gap,
                         printf("  ## Error: unable to allocate a new point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         return(-1)
                         ,c,0);
  }

  c[0] =  PRECI / mesh->info.delta *(mesh->info.max[0]-mesh->info.min[0])
    + 0.5;//mesh->info.max[0] + 1.;
  c[1] = PRECI / mesh->info.delta *(mesh->info.max[1]-mesh->info.min[1])
    + 0.5;//mesh->info.max[1] + 1.;
  ip4 = _MMG2D_newPt(mesh,c,0);
  if ( !ip4 ) {
    /* reallocation of point table */
    _MMG2D_POINT_REALLOC(mesh,sol,ip4,mesh->gap,
                         printf("  ## Error: unable to allocate a new point\n");
                         _MMG5_INCREASE_MEM_MESSAGE();
                         return(-1)
                         ,c,0);
  }

  assert(ip1==(mesh->np-3));
  assert(ip2==(mesh->np-2));
  assert(ip3==(mesh->np-1));
  assert(ip4==(mesh->np));
  /*creation des deux premiers triangles + adjacence*/
  jel  = _MMG2D_newElt(mesh);
  if ( !jel ) {
    _MMG5_TRIA_REALLOC(mesh,jel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
  }
  pt   = &mesh->tria[jel];
  pt->v[0] = ip1;
  pt->v[1] = ip4;
  pt->v[2] = ip2;
  pt->base = mesh->base;

  kel  = _MMG2D_newElt(mesh);
  if ( !kel ) {
    _MMG5_TRIA_REALLOC(mesh,kel,mesh->gap,
                       printf("  ## Error: unable to allocate a new element.\n");
                       _MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");
                       exit(EXIT_FAILURE));
  }
  pt   = &mesh->tria[kel];
  pt->v[0] = ip1;
  pt->v[1] = ip3;
  pt->v[2] = ip4;
  pt->base = mesh->base;

  iadr = 3*(jel-1) + 1;
  adja = &mesh->adja[iadr];
  adja[2] = 3*kel + 1;

  iadr = 3*(kel-1) + 1;
  adja = &mesh->adja[iadr];
  adja[1] = 3*jel + 2;

  /*vertex insertion*/
  if(!MMG2_insertpointdelone(mesh,sol)) return(0);
  fprintf(stdout,"  -- END OF INSERTION PHASE\n");

  /*bdry enforcement*/
  if(!MMG2_bdryenforcement(mesh,sol)) {
    printf("  ## Error: unable to enforce the boundaries.\n");
    return(0);
  }

  if(mesh->info.ddebug)
    if ( !_MMG5_chkmsh(mesh,1,0) ) exit(EXIT_FAILURE);

  /*mark SD and remove BB*/
  if(mesh->na)
    MMG2_markSD(mesh);
  else {
    /*tag des triangles : in = base ; out = -base ; indetermine = 0*/
    if(!MMG2_settagtriangles(mesh,sol)) return(0);
    if(!MMG2_removeBBtriangles(mesh)) return(0);
  }

  /*mark vertex bdry and compute tangent*/
  MMG2_baseBdry(mesh);

  return(1);
}
