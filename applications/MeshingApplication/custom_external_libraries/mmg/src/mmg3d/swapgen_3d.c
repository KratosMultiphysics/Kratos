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
 * \file mmg3d/swapgen_3d.c
 * \brief Functions for swapping process inside the mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "inlined_functions_3d.h"

/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the metric structure.
 * \param start tetrahedra in which the swap should be performed
 * \param ia edge that we want to swap
 * \param ilist pointer to store the size of the shell of the edge
 * \param list pointer to store the shell of the edge
 * \param crit improvment coefficient
 * \return 0 if fail, the index of point corresponding to the swapped
 * configuration otherwise (\f$4*k+i\f$).
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 *
 * Check whether swap of edge \a ia in \a start should be performed, and
 * return \f$4*k+i\f$ the index of point corresponding to the swapped
 * configuration. The shell of edge is built during the process.
 *
 */
int _MMG5_chkswpgen(MMG5_pMesh mesh,MMG5_pSol met,int start,int ia,
                    int *ilist,int *list,double crit,char typchk) {
  MMG5_pTetra    pt,pt0;
  MMG5_pPoint    p0;
  double    calold,calnew,caltmp;
  int       na,nb,np,adj,piv,npol,refdom,k,l,iel;
  int       *adja,pol[MMG3D_LMAX+2];
  char      i,ipa,ipb,ip,ier;

  pt  = &mesh->tetra[start];
  refdom = pt->ref;

  pt0 = &mesh->tetra[0];
  na  = pt->v[_MMG5_iare[ia][0]];
  nb  = pt->v[_MMG5_iare[ia][1]];
  calold = pt->qual;

  /* Store shell of ia in list, and associated pseudo polygon in pol */
  (*ilist) = 0;
  npol = 0;
  list[(*ilist)] = 6*start+ia;
  (*ilist)++;
  adja = &mesh->adja[4*(start-1)+1];
  adj  = adja[_MMG5_ifar[ia][0]];      // start travelling by face (ia,0)
  piv  = pt->v[_MMG5_ifar[ia][1]];
  pol[npol] = 4*start + _MMG5_ifar[ia][1];
  npol++;

  while ( adj ) {
    adj /= 4;
    if ( adj ==start ) break;

    pt = &mesh->tetra[adj];
    if ( pt->tag & MG_REQ ) return(0);

    /* Edge is on a boundary between two different domains */
    if ( pt->ref != refdom )  return(0);
    calold = MG_MIN(calold, pt->qual);
    /* identification of edge number in tetra adj */
    for (i=0; i<6; i++) {
      ipa = _MMG5_iare[i][0];
      ipb = _MMG5_iare[i][1];
      if ( (pt->v[ipa] == na && pt->v[ipb] == nb) ||
           (pt->v[ipa] == nb && pt->v[ipb] == na))  break;
    }
    assert(i<6);
    list[(*ilist)] = 6*adj +i;
    (*ilist)++;
    /* overflow */
    if ( (*ilist) > MMG3D_LMAX-3 )  return(0);

    /* set new triangle for travel */
    adja = &mesh->adja[4*(adj-1)+1];
    if ( pt->v[ _MMG5_ifar[i][0] ] == piv ) {
      pol[npol] = 4*adj + _MMG5_ifar[i][1];
      npol++;
      adj = adja[ _MMG5_ifar[i][0] ];
      piv = pt->v[ _MMG5_ifar[i][1] ];
    }
    else {
      assert(pt->v[ _MMG5_ifar[i][1] ] == piv);
      pol[npol] = 4*adj + _MMG5_ifar[i][0];
      npol++;
      adj = adja[ _MMG5_ifar[i][1] ];
      piv = pt->v[ _MMG5_ifar[i][0] ];
    }
  }
  //CECILE : je vois pas pourquoi ca ameliore de faire ce test
  //plus rapide mais du coup on elimine des swap...
  //4/01/14 commentaire
  // if ( calold*_MMG5_ALPHAD > 0.5 )  return(0);

  /* Prevent swap of an external boundary edge */
  if ( !adj )  return(0);

  assert(npol == (*ilist)); // du coup, apres on pourra virer npol

  /* Find a configuration that enhances the worst quality within the shell */
  for (k=0; k<npol; k++) {
    iel = pol[k] / 4;
    ip  = pol[k] % 4;
    np  = mesh->tetra[iel].v[ip];
    calnew = 1.0;
    ier = 1;

    if ( mesh->info.fem ) {
      p0 = &mesh->point[np];
      if ( p0->tag & MG_BDY ) {
        for (l=0; l<npol;l++) {
          if ( k < npol-1 ) {
            if ( l == k || l == k+1 )  continue;
          }
          else {
            if ( l == npol-1 || l == 0 )  continue;
          }
          iel = pol[l] / 4;
          ip  = pol[l] % 4;
          pt = &mesh->tetra[iel];
          p0 = &mesh->point[pt->v[ip]];
          if ( p0->tag & MG_BDY ) {
            ier = 0;
            break;
          }
        }
      }
      if ( !ier )  continue;
      ier = 1;
    }

    for (l=0; l<(*ilist); l++) {
      /* Do not consider tets of the shell of collapsed edge */
      if ( k < npol-1 ) {
        if ( l == k || l == k+1 )  continue;
      }
      else {
        if ( l == npol-1 || l == 0 )  continue;
      }
      iel = list[l] / 6;
      i   = list[l] % 6;
      pt  = &mesh->tetra[iel];

      /* Prevent from creating a tetra with 4 ridges vertices */
      if ( mesh->point[np].tag & MG_GEO ) {
        if ( ( mesh->point[pt->v[_MMG5_ifar[i][0]]].tag & MG_GEO ) &&
             ( mesh->point[pt->v[_MMG5_ifar[i][1]]].tag & MG_GEO ) ) {
          if ( ( mesh->point[_MMG5_iare[i][0]].tag & MG_GEO ) ||
               ( mesh->point[_MMG5_iare[i][1]].tag & MG_GEO ) )
            return(0);
        }
      }

      /* First tetra obtained from iel */
      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[_MMG5_iare[i][0]] = np;


      if ( typchk==1 && met->size > 1 && met->m )
        caltmp = _MMG5_caltet33_ani(mesh,met,pt0);
      else
        caltmp = _MMG5_orcal(mesh,met,0);

      calnew = MG_MIN(calnew,caltmp);

      ier = (calnew > crit*calold);
      if ( !ier )  break;

      /* Second tetra obtained from iel */
      memcpy(pt0,pt,sizeof(MMG5_Tetra));
      pt0->v[_MMG5_iare[i][1]] = np;

      if ( typchk==1 && met->size > 1 && met->m )
        caltmp = _MMG5_caltet33_ani(mesh,met,pt0);
      else
        caltmp = _MMG5_orcal(mesh,met,0);

      calnew = MG_MIN(calnew,caltmp);

      ier = (calnew > crit*calold);
      if ( !ier )  break;
    }
    if ( ier )  return(pol[k]);
  }
  return(0);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param nconf configuration.
 * \param ilist number of tetrahedra in the shell of the edge that we want
 *  to swap.
 * \param list pointer toward the shell of the edge that we want to swap.
 * \param octree pointer toward the octree structure in Delaunay mode,
 * NULL pointer in pattern mode.
 * \param typchk type of checking permformed for edge length (hmin or LSHORT
 * criterion).
 * \return -1 if lack of memory, 0 if fail to swap, 1 otherwise.
 *
 * Perform swap of edge whose shell is passed according to configuration nconf.
 *
 */
int _MMG5_swpgen(MMG5_pMesh mesh,MMG5_pSol met,int nconf,int ilist,int *list,
                 _MMG3D_pOctree octree, char typchk) {
  MMG5_pTetra    pt;
  MMG5_pPoint    p0,p1;
  int       iel,na,nb,np,nball,ret,start;
  double    m[3];
  char      ia,ip,iq;
  int       ier;

  iel = list[0] / 6;
  ia  = list[0] % 6;

  pt = &mesh->tetra[iel];
  na = pt->v[_MMG5_iare[ia][0]];
  nb = pt->v[_MMG5_iare[ia][1]];
  p0 = &mesh->point[na];
  p1 = &mesh->point[nb];

  /* Temporarily create midpoint at swapped edge */
  m[0] = 0.5*(p0->c[0] + p1->c[0]);
  m[1] = 0.5*(p0->c[1] + p1->c[1]);
  m[2] = 0.5*(p0->c[2] + p1->c[2]);

  np  = _MMG3D_newPt(mesh,m,0);
  if(!np){
    _MMG5_POINT_REALLOC(mesh,met,np,mesh->gap,
                        printf("  ## Error: unable to allocate a new point\n");
                        _MMG5_INCREASE_MEM_MESSAGE();
                        return(-1)
                        ,m,0);
  }
  if ( met->m ) {
    if ( typchk == 1 && (met->size>1) ) {
      if ( _MMG3D_intmet33_ani(mesh,met,iel,ia,np,0.5)<=0 )  return(0);
    }
    else {
      if ( _MMG5_intmet(mesh,met,iel,ia,np,0.5)<=0 ) return(0);
    }
  }

  /** First step : split of edge (na,nb) */
  ret = 2*ilist + 0;
  ier = _MMG5_split1b(mesh,met,list,ret,np,0,typchk-1);
  /* pointer adress may change if we need to realloc memory during split */
  pt = &mesh->tetra[iel];

  if ( ier < 0 ) {
    fprintf(stdout,"  ## Warning: unable to swap internal edge.\n");
    return(-1);
  }
  else if ( !ier )  return(0);

  /** Second step : collapse of np towards enhancing configuration */
  start = nconf / 4;
  iq = nconf % 4;

  pt = &mesh->tetra[start];
  for (ip=0; ip<4; ip++) {
    if ( pt->v[ip] == np )  break;
  }
  assert(ip<4);

  memset(list,0,(MMG3D_LMAX+2)*sizeof(int));
  nball = _MMG5_boulevolp(mesh,start,ip,list);

  ier = _MMG5_colver(mesh,met,list,nball,iq,typchk);
  if ( ier < 0 ) {
    fprintf(stdout,"  ## Warning: unable to swap internal edge.\n");
    return(-1);
  }
  else if ( ier ) _MMG3D_delPt(mesh,ier);

  return(1);
}
