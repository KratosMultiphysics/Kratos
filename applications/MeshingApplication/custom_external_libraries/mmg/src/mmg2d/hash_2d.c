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

#define KTA     7
#define KTB    11

/*hash edge :
  return 1 if edge exist in the table*/
int MMG2_FindEdge(pHashTable edgeTable,int ia, int ib) {
  int       key,mins,maxs;
  Hedge     *ha;

  /* compute key */
  if ( ia < ib ) {
    mins = ia;
    maxs = ib;
  }
  else {
    mins = ib;
    maxs = ia;
  }
  key = KTA*mins + KTB*maxs;
  key = key % edgeTable->size;
  ha  = &edgeTable->item[key];
  if ( ha->min ) {
    /* edge exist*/
    if ( ha->min == mins && ha->max == maxs ) {
      return(ha->iel);
    }
    else {
      while ( ha->nxt && ha->nxt < edgeTable->nxtmax ) {
        ha = &edgeTable->item[ha->nxt];
        if ( ha->min == mins && ha->max == maxs ) {
          return(ha->iel);
        }
      }
      ha->nxt = edgeTable->hnxt;
      ha      = &edgeTable->item[edgeTable->hnxt];
      ++edgeTable->hnxt;
      if ( edgeTable->hnxt == edgeTable->nxtmax ) {
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",edgeTable->nxtmax);
        assert(0);
        return(0);
      }
    }
  }
 
  return(0);

}

int MMG2_hashNew(HashTable *hash,int hsize,int hmax) {
  int   k;

  hash->size  = hsize;
  hash->nxtmax =hmax+1;
  hash->hnxt  = hsize;
  _MMG5_SAFE_CALLOC(hash->item,hash->nxtmax,Hedge);

  for (k=hash->size; k<hash->nxtmax; k++)
    hash->item[k].nxt = k+1;

  return(1);
}
int MMG2_hashel(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  int       k,kk,pp,l,ll,mins,mins1,maxs,maxs1;
  int      *hcode,*link,inival,hsize,iadr;
  unsigned char  *hvoy,i,ii,i1,i2;
  unsigned int    key;

  if ( mesh->adja )  return(1);
  if ( !mesh->nt )  return(0);

  /* memory alloc */
  _MMG5_SAFE_CALLOC(hcode,mesh->nt+1,int);

  /* memory alloc */
  _MMG5_ADD_MEM(mesh,(3*mesh->ntmax+5)*sizeof(int),"adjacency table",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->adja,3*mesh->ntmax+5,int);

  link  = mesh->adja;
  hsize = mesh->nt;
  hvoy  = (unsigned char*)hcode;

  /* init */
  inival = 2147483647;
  for (k=0; k<=mesh->nt; k++)
    hcode[k] = -inival;

  /* build hash table */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !pt->v[0] )  continue;
    for (i=0; i<3; i++) {
      i1 = MMG2_idir[i+1];
      i2 = MMG2_idir[i+2];
      if ( pt->v[i1] < pt->v[i2] ) {
        mins = pt->v[i1];
        maxs = pt->v[i2];
      }
      else {
        mins = pt->v[i2];
        maxs = pt->v[i1];
      }

      /* compute key */
      key = KTA*mins + KTB*maxs;
      key = key % hsize + 1;

      /* insert */
      iadr = 3*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=3*mesh->nt; l>0; l--) {
    if ( link[l] >= 0 )  continue;
    k = (l-1) / 3 + 1;
    i = (l-1) % 3;
    i1 = MMG2_idir[i+1];
    i2 = MMG2_idir[i+2];
    pt = &mesh->tria[k];

    mins = M_MIN(pt->v[i1],pt->v[i2]);
    maxs = M_MAX(pt->v[i1],pt->v[i2]);

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    hvoy[l] = 0;
    while ( ll != inival ) {
      kk = (ll-1) / 3 + 1;
      ii = (ll-1) % 3;
      i1 = MMG2_idir[ii+1];
      i2 = MMG2_idir[ii+2];
      pt1  = &mesh->tria[kk];
      if ( pt1->v[i1] < pt1->v[i2] ) {
        mins1 = pt1->v[i1];
        maxs1 = pt1->v[i2];
      }
      else {
        mins1 = pt1->v[i2];
        maxs1 = pt1->v[i1];
      }

      if ( mins1 == mins  && maxs1 == maxs ) {
        /* adjacent found */
        if ( pp != 0 )  link[pp] = link[ll];
        link[l] = 3*kk + ii;
        link[ll]= 3*k + i;
        break;
      }
      pp = ll;
      ll = -link[ll];
    }
  }
  _MMG5_SAFE_FREE(hcode);

  MMG2_baseBdry(mesh);
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param ip1,ip2,ip3 integer
 * \param t  : computed tangent at vertex ip2
 *
 * Compute the tangent at vertex ip2
 *
 */
int MMG2_computetangent(MMG5_pMesh mesh,int ip1,int ip2,int ip3,double *t) {
  double c2[2],c1[2],c3[2],dd1,dd2;
  double dd,n1[2],n2[2],n[2];
//  double t1[2],t2[2],theta,kappa;
  int    i;

  for(i=0 ; i<2 ; i++) {
    c1[i] =  mesh->point[ip1].c[i];
    c2[i] =  mesh->point[ip2].c[i];
    c3[i] =  mesh->point[ip3].c[i];
  }
  if(mesh->point[ip2].tag & M_CORNER) {
    t[0] = c2[0]-c1[0];
    t[1] = c2[1]-c1[1];
  }


  /*mean normal*/
  n1[0] = - (c2[1] - c1[1]);
  n1[1] =  (c2[0] - c1[0]);
  dd1 = 1;//sqrt(n1[0]*n1[0]+n1[1]*n1[1]);
  n2[0] = - (c3[1] - c2[1]);
  n2[1] =  (c3[0] - c2[0]);
  dd2 = 1;//sqrt(n2[0]*n2[0]+n2[1]*n2[1]);

  dd = 0;
  for(i=0 ; i<2 ; i++) {
    n[i] = 0.5*(n1[i]/dd1+n2[i]/dd2);
    dd += n[i]*n[i];
  }
  dd = sqrt(dd);
  n[0] /= dd;
  n[1]/=dd;
  t[0] = n[1];
  t[1] = -n[0];
  /* t1[0] = c2[0] - c1[0]; */
  /* t1[1] = c2[1] - c1[1]; */
  /* t2[0] = c3[0] - c2[0]; */
  /* t2[1] = c3[1] - c2[1]; */
  /* theta = (t1[0]*t2[0]+t1[1]*t2[1])/(sqrt(t1[0]*t1[0]+t1[1]*t1[1])*sqrt(t2[0]*t2[0]+t2[1]*t2[1])); */
  /* //printf("theta %e %e == %e\n",theta,acos(theta),M_Pi/2.); */
  /* theta = acos(theta); */
  /* //sign */
  /* dd = t1[0]*t2[1]-t1[1]*t2[0]; */
  /* if(dd<0) {/\*puts("eheheheheheheheh 222");*\/theta *= -1.;} */
  /* dd = sqrt(t1[0]*t1[0]+t1[1]*t1[1]) +sqrt(t2[0]*t2[0]+t2[1]*t2[1]); */
  /* //printf("arc length %e\n",dd); */
  /* kappa = theta / (0.5*dd); */
  /* printf("pour %d on trouve kappa %e\n",ip2,kappa); */

  /* if(fabs(kappa)<1e-8){ */
  /*   printf("STRAIGHT EDGE\n"); */
  /*   return(1); */
  /* }  */

  /* r = 1./kappa; */
  /* if(r>0) { */
  /*   n[0] = - (c2[1] - c1[1]); */
  /*   n[1] =  (c2[0] - c1[0]); */
  /* } else { */
  /*   n[0] =  (c2[1] - c1[1]); */
  /*   n[1] = - (c2[0] - c1[0]); */
  /*   r *= -1; */
  /* } */
  /* dd = sqrt(n[0]*n[0]+n[1]*n[1]); */
  /* n[0] /= dd; */
  /* n[1] /= dd; */
  /* //can be < 0 !! ==> straight edge */
  /* if((r*r - 0.25*dd*dd) < 0) { */
  /*   printf("STRAIGHT EDGE 2\n"); */
  /*   return(1); */
  /* } */
  /* center[0] = 0.5*(c1[0]+c2[0]) + (sqrt(r*r - 0.25*dd*dd))*(n[0]); */
  /* center[1] = 0.5*(c1[1]+c2[1]) + (sqrt(r*r - 0.25*dd*dd))*(n[1]); */

  return(1);
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param ped pointer toward an Edge.
 * \param k integer
 * \param ik
 * \param i2
 * \param nv
 * \return 1
 *
 * Compute the tangent.
 *
 */
int MMG2_tangent(MMG5_pMesh mesh,MMG5_pEdge ped,int k,int ik,int i2,int nv) {
  MMG5_pTria pt,ptiel;
  MMG5_pEdge pediel;
  MMG5_pPoint pta,ptb,ptaiel,ptbiel;

  int iadj,ip,i,iadr,*adja,iel,voyiel,vi1,vi2,num,ipnv;

  pt = &mesh->tria[k];
  iadr = 3*(k-1)+1;
  adja = &mesh->adja[iadr];
  iadj = i2;
  ip = pt->v[ik];
  pta = &mesh->point[ped->a];
  ptb = &mesh->point[ped->b];

  if(!adja[iadj]) {
    if(nv) {
      MMG2_computetangent(mesh,pt->v[i2],ped->b,ip,ptb->n);
    }
    else {
      MMG2_computetangent(mesh,pt->v[i2],ped->a,ip,pta->n);
    }
    num = pt->edg[i2];
    assert(num);
    pediel = &mesh->edge[num];
    ptaiel = &mesh->point[pediel->a];
    ptbiel = &mesh->point[pediel->b];
    ipnv = (nv==0)?ped->a:ped->b;
    if(pediel->a==ipnv){
      for(i=0 ; i<2 ; i++) {
        if(nv)
          ptaiel->n[i] = ptb->n[i];
        else
          ptaiel->n[i] = pta->n[i];
      }
    } else {
      assert(pediel->b==ipnv);
      for(i=0 ; i<2 ; i++) {
        if(nv)
          ptbiel->n[i] = ptb->n[i];
        else
          ptbiel->n[i] = pta->n[i];
      }
    }
    return(1);

  }
  do {
    iel = adja[iadj]/3;
    ptiel = &mesh->tria[iel];
    voyiel = adja[iadj]%3;
    /*find ptiel*/
    iadr = 3*(iel-1) + 1;
    adja = &mesh->adja[iadr];

    vi1 = MMG2_idir[voyiel+1];
    vi2 = MMG2_idir[voyiel+2];
    if(ptiel->v[vi2]==ip) {
      num = ptiel->edg[vi2];
      iadj = vi2;
    } else {
      num = ptiel->edg[vi1];
      iadj = vi1;
    }
    ip = ptiel->v[voyiel];
  } while (!num && iel);
  if(nv) {
    MMG2_computetangent(mesh,pt->v[i2],ped->b,ptiel->v[voyiel],ptb->n);
  }
  else {
    MMG2_computetangent(mesh,pt->v[i2],ped->a,ptiel->v[voyiel],pta->n);
  }

  return(1);
}

/* base boundary vertices and compute tangents*/
int MMG2_baseBdry(MMG5_pMesh mesh) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppt;
  MMG5_pEdge     ped;
  int      *adja,adj,iadr,k,i,ip,ned,num,i1,i2;
  HashTable edgeT;

  /*edge treatment*/
  edgeT.size  = mesh->namax;
  edgeT.nxtmax = 3*mesh->namax+1;
  edgeT.hnxt  = mesh->namax;
  _MMG5_SAFE_CALLOC(edgeT.item,edgeT.nxtmax,Hedge);
  
  memset(edgeT.item,0,edgeT.nxtmax*sizeof(Hedge));
  for (k=edgeT.size; k<edgeT.nxtmax; k++)
    edgeT.item[k].nxt = k+1;
  ned = 0;
  if(mesh->na) {
    ned = mesh->na;
 
    for(k=1 ; k<=mesh->na ; k++) {
      ped = &mesh->edge[k];
      if(!ped->a) continue;
      MMG2_hashEdge(&edgeT,k,ped->a,ped->b);
    }
  }

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;

    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[iadr];

    for (i=0; i<3; i++) {
      adj = adja[i]/3;
      pt1 = &mesh->tria[adj];
      pt->edg[i] = 0;
      if ( !adj  ) {
        num = 0;
        if(ned) {
          num = MMG2_FindEdge(&edgeT,pt->v[MMG2_iopp[i][0]],pt->v[MMG2_iopp[i][1]]);
        }
        ip  = pt->v[MMG2_iopp[i][0]];
        mesh->point[ip].tag |= M_BDRY;
        if ( mesh->info.nosurf && ( !( mesh->point[ip].tag & M_REQUIRED) ) ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }

        ip  = pt->v[MMG2_iopp[i][1]];
        mesh->point[ip].tag |= M_BDRY;
        if ( mesh->info.nosurf && ( !( mesh->point[ip].tag & M_REQUIRED) ) ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }

        if(num) {
          pt->edg[i] = num;
        } else {
          num = _MMG5_newEdge(mesh);
          if ( !num ) {
            _MMG5_EDGE_REALLOC(mesh,num,mesh->gap,
                               printf("  ## Error: unable to allocate a new edge.\n");
                               _MMG5_INCREASE_MEM_MESSAGE();
                               printf("  Exit program.\n");
                               exit(EXIT_FAILURE));
          }
          pt->edg[i] = num;
          ped = &mesh->edge[num];
          ped->a = pt->v[MMG2_iopp[i][0]];
          ped->b = pt->v[MMG2_iopp[i][1]];
        }
      } else if(pt->ref != pt1->ref) {
        num = 0;
        if(ned) {
          num = MMG2_FindEdge(&edgeT,pt->v[MMG2_iopp[i][0]],pt->v[MMG2_iopp[i][1]]);
        }
        ip  = pt->v[MMG2_iopp[i][0]];
        mesh->point[ip].tag |= M_SD;
        if ( mesh->info.nosurf && ( !( mesh->point[ip].tag & M_REQUIRED) ) ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }

        ip  = pt->v[MMG2_iopp[i][1]];
        mesh->point[ip].tag |= M_SD;
        if ( mesh->info.nosurf  && ( !( mesh->point[ip].tag & M_REQUIRED) )  ) {
          mesh->point[ip].tag |= M_REQUIRED;
          mesh->point[ip].tag |= M_NOSURF;
        }

        if(num) {
          pt->edg[i] = num;
        } else {
          num = _MMG5_newEdge(mesh);
          if ( !num ) {
            _MMG5_EDGE_REALLOC(mesh,num,mesh->gap,
                               printf("  ## Error: unable to allocate a new edge.\n");
                               _MMG5_INCREASE_MEM_MESSAGE();
                               printf("  Exit program.\n");
                               exit(EXIT_FAILURE));
          }
          pt->edg[i] = num;
          ped = &mesh->edge[num];
          ped->a = pt->v[MMG2_iopp[i][0]];
          ped->b = pt->v[MMG2_iopp[i][1]];
          ned++;
          MMG2_hashEdge(&edgeT,num,ped->a,ped->b);
        }
      }
    }
  }

  if(edgeT.item)  _MMG5_SAFE_FREE(edgeT.item);

  /*compute tangents*/
  mesh->base++;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !M_EOK(pt) )  continue;
    iadr = 3*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for (i=0; i<3; i++) {
      num = pt->edg[i];
      if(!num) continue;
      ped = &mesh->edge[num];

      i1 = MMG2_idir[i+1];
      i2 = MMG2_idir[i+2];
      ppt = &mesh->point[ped->a];
      /*tangent on ped->a*/
      if(ppt->flag != mesh->base) {
        if(pt->v[i1]==ped->a) {
          assert(pt->v[i2]==ped->b);
          MMG2_tangent(mesh,ped,k,i,i2,0);
          ppt->flag = mesh->base;
        } else {
          assert(pt->v[i2]==ped->a);
          MMG2_tangent(mesh,ped,k,i,i1,0);
          ppt->flag = mesh->base;
        }

      } else {
        //printf("on a deja calcule pour %d %e %e\n",ped->a,ped->t0[0],ped->t0[1]);
      }



      /* ppt = &mesh->point[ped->b]; */
      /* /\*tangent on ped->b*\/ */
      /* if(ppt->flag != mesh->flag) { */
      /*  if(pt->v[i1]==ped->b) { */
      /*    MMG2_tangent(mesh,ped,k,i,i2,1); */
      /*  } else { */
      /*    MMG2_tangent(mesh,ped,k,i,i1,1);     */
      /*  } */
      /*  ppt->flag = mesh->flag; */

      /* } else { */
      /*  //printf("on a deja calcule pour %d %e %e\n",ped->b,ped->t1[0],ped->t1[1]); */
      /* } */

    }
  }
  return(1);
}


/*hash edge :
  return 1 if edge exist in the table*/
int MMG2_hashEdge(pHashTable edgeTable,int iel,int ia, int ib) {
  int       key,mins,maxs;
  Hedge     *ha;

  /* compute key */
  if ( ia < ib ) {
    mins = ia;
    maxs = ib;
  }
  else {
    mins = ib;
    maxs = ia;
  }

  key = KTA*mins + KTB*maxs;
  key = key % edgeTable->size;
  ha  = &edgeTable->item[key];
  if ( ha->min ) {
    /* edge exist*/
    if ( ha->min == mins && ha->max == maxs ) {
      return(ha->iel);
    }
    else {
      while ( ha->nxt && ha->nxt < edgeTable->nxtmax ) {
        ha = &edgeTable->item[ha->nxt];
        if ( ha->min == mins && ha->max == maxs )
          return(ha->iel);
      }
      ha->nxt = edgeTable->hnxt;
      ha      = &edgeTable->item[edgeTable->hnxt];
      ++edgeTable->hnxt;
      if ( edgeTable->hnxt == edgeTable->nxtmax ) {
        fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",edgeTable->nxtmax);
        assert(0);
        return(0);
      }
    }
  }
  /* insert */
  ha->min = mins;
  ha->max = maxs;
  ha->iel = iel;
  ha->nxt = 0;

  return(0);

}
