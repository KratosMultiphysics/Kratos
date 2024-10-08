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
 * \file mmg3d/opttyp_3d.c
 * \brief Functions for the optimization of very bad elements.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg3d.h"
#include "inlined_functions_3d_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param iel element index.
 * \param item bad entity.
 * \return -1 if fail, element type otherwise.
 *
 * Identify the type of element.\n
 * Type:
 *    - 0: element is ok.
 *    - 1: empty volume but 4 valid faces (sliver)
 *    - 2: empty volume and a vertex close to the opposite face.
 * 4 valid faces (chinese hat).
 *    - 3: 3 valid faces, 1 obtuse face (aileron)
 *    - 4: 2 valid faces, 2 acute faces (=> 1 small edge)
 *    - 5: 1 valid face, 3 small edges
 *    - 6: 2 faces with big edges, 2 faces with small edges
 *    - 7: 4 faces with big edges
 *
 * Element caracteristics by type:
 *    - 0: 0 obtuse face, 0 acute face, 0 big edge, 0 small edge.
 *    - 1: 0 obtuse face, 0 acute face, 0 big edge, 0 small edge.
 *    - 2: 0 obtuse face, 0 acute face, 0 big edge, 0 small edge.
 *    - 3: 1 obtuse face, 0 acute face, 1 big edge, 0 small edge.
 *    - 4: 0 obtuse face, 2 acute face, 0 big edge, 1 small edge.
 *    - 5: 0 obtuse face, 3 acute face, 0 big edge, 3 small edge.
 *    - 6: 2 obtuse face, 2 acute face, 1 big edge, 1 small edge.
 *    - 7: 0 obtuse face, 4 acute face, 0 big edge, 2 small edge.
 *
 */
static int MMG3D_typelt(MMG5_pMesh mesh,MMG5_int iel,int *item) {
  MMG5_pTetra    pt;
  MMG5_pPoint    pa,pb,pc,pd;
  double         abx,aby,abz,acx,acy,acz,adx,ady,adz,v1,v2,v3,vol;
  double         bcx,bcy,bcz,bdx,bdy,bdz,cdx,cdy,cdz,h[6],volchk,ssmall;
  double         s[4],dd,rapmin,rapmax,surmin,surmax;
  int            i,k,isur,isurmax,isurmin,iarmax,iarmin;
  MMG5_int       ia,ib,ic,id;
  int            nobtus,naigu;
  short          i0,i1,i2;
  double         EPSVOL = 0.001;
  double         RAPMAX = 0.25;

  pt = &mesh->tetra[iel];
  if ( !pt->v[0] )  return -1;

  ia = pt->v[0];
  ib = pt->v[1];
  ic = pt->v[2];
  id = pt->v[3];
  pa = &mesh->point[ia];
  pb = &mesh->point[ib];
  pc = &mesh->point[ic];
  pd = &mesh->point[id];

  /* volume */
  abx = pb->c[0] - pa->c[0];
  aby = pb->c[1] - pa->c[1];
  abz = pb->c[2] - pa->c[2];

  acx = pc->c[0] - pa->c[0];
  acy = pc->c[1] - pa->c[1];
  acz = pc->c[2] - pa->c[2];

  adx = pd->c[0] - pa->c[0];
  ady = pd->c[1] - pa->c[1];
  adz = pd->c[2] - pa->c[2];

  v1  = acy*adz - acz*ady;
  v2  = acz*adx - acx*adz;
  v3  = acx*ady - acy*adx;
  vol = abx * v1 + aby * v2 + abz * v3;

  /* max edge */
  h[0] = abx*abx + aby*aby + abz*abz;
  h[1] = acx*acx + acy*acy + acz*acz;
  h[2] = adx*adx + ady*ady + adz*adz;

  bcx = pc->c[0] - pb->c[0];
  bcy = pc->c[1] - pb->c[1];
  bcz = pc->c[2] - pb->c[2];

  bdx = pd->c[0] - pb->c[0];
  bdy = pd->c[1] - pb->c[1];
  bdz = pd->c[2] - pb->c[2];

  cdx = pd->c[0] - pc->c[0];
  cdy = pd->c[1] - pc->c[1];
  cdz = pd->c[2] - pc->c[2];

  h[3] = bcx*bcx + bcy*bcy + bcz*bcz;
  h[4] = bdx*bdx + bdy*bdy + bdz*bdz;
  h[5] = cdx*cdx + cdy*cdy + cdz*cdz;

  /* face areas */
  dd = cdy*bdz - cdz*bdy;
  s[0] = dd * dd;
  dd = cdz*bdx - cdx*bdz;
  s[0] = s[0] + dd * dd;
  dd = cdx*bdy - cdy*bdx;
  s[0] = s[0] + dd * dd;
  s[0] = sqrt(s[0]);

  s[1] = sqrt(v1*v1 + v2*v2 + v3*v3);

  dd = bdy*adz - bdz*ady;
  s[2] = dd * dd;
  dd = bdz*adx - bdx*adz;
  s[2] = s[2] + dd * dd;
  dd = bdx*ady - bdy*adx;
  s[2] = s[2] + dd * dd;
  s[2] = sqrt(s[2]);

  dd = aby*acz - abz*acy;
  s[3] = dd * dd;
  dd = abz*acx - abx*acz;
  s[3] = s[3] + dd * dd;
  dd = abx*acy - aby*acx;
  s[3] = s[3] + dd * dd;
  s[3] = sqrt(s[3]);

  /* classification */
  rapmin = h[0];
  rapmax = h[0];
  iarmin = 0;
  iarmax = 0;
  for (i=1; i<6; i++) {
    if ( h[i] < rapmin ) {
      rapmin = h[i];
      iarmin = i;
    }
    else if ( h[i] > rapmax ) {
      rapmax = h[i];
      iarmax = i;
    }
  }
  rapmin = sqrt(rapmin);
  rapmax = sqrt(rapmax);
  volchk = EPSVOL * rapmin*rapmin*rapmin;

  /* small volume: types 1,2,3,4 */
  item[0] = item[1] = -1;
  if ( vol < volchk ) {
    // puts("volume nul : type 1,2,3,4");
    ssmall = 0.4 * (s[0]+s[1]+s[2]+s[3]);
    isur   = 0;
    for (i=0; i<4; i++)
      isur += s[i] > ssmall;

    /* types 2,3 */
    item[0] = iarmax;
    item[1] = MMG5_isar[iarmax][0];
    if ( isur == 1 ) {
      surmin   = s[0];
      isurmin = 0;
      surmax   = s[0];
      isurmax = 0;
      for (i=1; i<4; i++) {
        if ( s[i] < surmin ) {
          surmin  = s[i];
          isurmin = i;
        }
        else if ( s[i] > surmax ) {
          surmax  = s[i];
          isurmax = i;
        }
      }
      dd = surmin / surmax;
      if ( dd < RAPMAX ) {
        item[1] = MMG5_isar[iarmax][0];
        return 3;
      }
      else {
        item[0] = isurmax;
        item[1] = isurmin;
        return 2;
      }
    }

    /* type 1 */
    isur = 0;
    if ( s[0]+s[1] > ssmall )  isur = 1;
    if ( s[0]+s[2] > ssmall )  isur++;
    if ( s[0]+s[3] > ssmall )  isur++;

    if ( isur > 2 ) {
      dd = rapmin / rapmax;
      item[0] = 0; //iarmin;
      // Changed by 0 because we overflow idir
      item[1] = 0; //MMG5_idir[iarmin][0];
      if ( dd < 0.01 )  return 4;
      if ( s[0]+s[1] > ssmall ) {
        item[0] = 0;
        return 1;
      }
      if ( s[0]+s[2] > ssmall ) {
        item[0] = 1;
        return 1;
      }
      if ( s[0]+s[3] > ssmall ) {
        item[0] = 2;
        return 1;
      }
    }

    //puts("default");
    item[0] = 0;
    return 1;
  }/* end chkvol */

  dd = rapmin / rapmax;
  /* types 3,6,7 */
  if ( dd < RAPMAX ) {
    /* One edge is 4 time smaller than another one */

    for (i=0; i<6; i++)  h[i] = sqrt(h[i]);

    nobtus = 0;
    for (k=0; k<4; k++) {
      for (i=0; i<3; i++) {
        i0 = MMG5_idir[k][i];
        i1 = MMG5_idir[k][MMG5_inxt2[i]];
        i2 = MMG5_idir[k][MMG5_inxt2[i+1]];
        if ( h[i0]+h[i1] < 1.2*h[i2] ) {
          /* Obtuse face */
          nobtus++;
          item[0] = i2;
          item[1] = MMG5_idir[k][MMG5_inxt2[i+1]];
        }
      }
    }

    switch(nobtus){
    case 0 :
      break;
    case 1:
      item[0] = iarmax;
      item[1] = MMG5_isar[iarmax][0];
      return 3;
    case 2:
      item[0] = iarmin;
      item[1] = iarmax;
      return 6;
    default:
      item[0] = iarmin;
      item[1] = iarmax;
      return 7;
    }
  }

  /* type 4,5,7 */
  else if ( dd < 0.7*RAPMAX ) {
    naigu = 0;
    for (k=0; k<4; k++) {
      for (i=0; i<3; i++) {
        i0 = MMG5_idir[k][i];
        i1 = MMG5_idir[k][MMG5_inxt2[i]];
        i2 = MMG5_idir[k][MMG5_inxt2[i+1]];
        if ( h[i0]+h[i1] > 1.5*h[i2] )  naigu++;
      }
    }
    switch(naigu){
    case 0 :
      break;
    case 1:
      break;
    case 2:
      item[0] = iarmin;
      return 4;
    case 3:
      /* Item to define */
      return 5;
    default:
      item[0] = iarmin;
      item[1] = iarmax;
      return 7;
    }
  }
  item[0] = 0;
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k elt index.
 * \param iar index of edge to not try to swap.
 * \return -1 if fail, 0 if we don't swap anything, 1 otherwise.
 *
 * Try to swap edge \a iar of tetra \a k.
 *
 */
int MMG3D_swpItem(MMG5_pMesh mesh,  MMG5_pSol met,MMG3D_pPROctree PROctree,MMG5_int k,int iar) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  int64_t       list[MMG3D_LMAX+2];
  MMG5_int      nconf;
  int           lon,ier;
  double        OCRIT = 1.01;

  ier = 0;
  pt = &mesh->tetra[k];
  /* Prevent swap of a ref or tagged edge */
  if ( pt->xt ) {
    pxt = &mesh->xtetra[pt->xt];
    if ( pxt->edg[iar] || pxt->tag[iar] ) return 0;
  }

  nconf = MMG5_chkswpgen(mesh,met,k,iar,&lon,list,OCRIT,2);
  if ( nconf ) {
    ier = MMG5_swpgen(mesh,met,nconf,lon,list,PROctree,2);
    if ( ier < 0 ) return -1;
    else
      return ier;
  }

  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k elt index.
 * \param iar index of edge to not try to swap.
 * \return -1 if fail, 0 if we don't swap anything, 1 otherwise.
 *
 * Try to swap all edges of tetra \a k except of the edge number \a iar.
 *
 */
static inline
int MMG3D_swpalmostall(MMG5_pMesh mesh,  MMG5_pSol met,MMG3D_pPROctree PROctree,
                        MMG5_int k,int iar) {
  int           i,ier;

  ier = 0;
  for(i=0 ; i<6 ; i++) {
    if(i==iar) continue;
    ier = MMG3D_swpItem(mesh,met,PROctree,k,i);
    if ( ier < 0 ) return -1;
    else if(ier)
      return ier;
  }
  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k elt index.
 * \param iar index of edge to split.
 * \param OCRIT quality threshold.
 * \return 1 if success, 0 otherwise
 *
 * Try to split edge number \a iar of tetra \a k
 *
 */
int MMG3D_splitItem(MMG5_pMesh mesh,  MMG5_pSol met,MMG3D_pPROctree PROctree,
                     MMG5_int k,int iar,double OCRIT) {
  MMG5_pTetra   pt;
  double        len;
  double        LLONG2 = 0.1;
  int           j;
  MMG5_int      ier;

  ier = 0;
  pt = &mesh->tetra[k];

  if ( mesh->info.noinsert ) return 0;

  len = MMG5_lenedg(mesh,met,iar,pt);
  if(len > LLONG2) {
    ier = MMG5_splitedg(mesh,met,k,iar,OCRIT);
  }

  if ( ier && !mesh->info.nomove ) {
    /*if there is a reallocation inside splitedg, the pointer is unvalid*/
    pt = &mesh->tetra[k];
    for(j=0 ; j<4 ; j++) {
      if(pt->v[j] == ier) break;
    }
    assert(j<4);
    if(met->size!=1)
      ier = MMG3D_movv_ani(mesh,met,k,j);
    else
      ier = MMG3D_movv_iso(mesh,met,k,j);

  }
  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param k elt index.
 * \param iar index of edge to not split.
 * \return 1 if success, 0 otherwise
 *
 * Try to split evry edge of tetra \a k except of edge number \a iar.
 *
 */
static inline
int MMG3D_splitalmostall(MMG5_pMesh mesh,  MMG5_pSol met,MMG3D_pPROctree PROctree,
                          MMG5_int k,int iar) {
  int           i,ier;
  double        OCRIT=1.01;

  ier = 0;

  /* if(MMG5_orvolnorm(mesh,k) < 5.e-9) { */
  /*   OCRIT *= 0.5; */
  /* } else */
  /*   OCRIT *= 0.75; */

  for(i=0 ; i<6 ; i++) {
    if(i==iar) continue;
    ier = MMG3D_splitItem(mesh,met,PROctree,k,i,OCRIT);
    if(ier) return ier;
  }

  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param PROctree pointer toward the PROctree structure.
 * \param testmark all the tets with a mark less than testmark will not be treated.
 * \return 0 if fail, number of improved elts otherwise.
 *
 * Travel across the mesh to detect element with very bad quality (less than
 * 0.2) and try to improve them by every means.
 *
 */
MMG5_int MMG3D_opttyp(MMG5_pMesh mesh, MMG5_pSol met,MMG3D_pPROctree PROctree,MMG5_int testmark) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  double         crit;
  int            ityp,item[2];
  MMG5_int       k,ntot,ne,nd,cs[10],ds[10];
  int            ier,i,npeau;
  int            it,maxit;
//  double         OCRIT = 1.01;
  MMG5_int       nbdy,nbdy2,base ;

  ntot = 0;
  crit = 0.2 / MMG3D_ALPHAD;
  base = testmark;

  it = 0;
  maxit = 10;
  do {
    ne = mesh->ne;
    nd = 0;
    nbdy = nbdy2 = 0;
    memset(cs,0,10*sizeof(MMG5_int));
    memset(ds,0,10*sizeof(MMG5_int));

    for (k=1 ; k<=ne ; k++) {
      pt = &mesh->tetra[k];
      if(!MG_EOK(pt)  || (pt->tag & MG_REQ) ) continue;
      else if ( pt->mark < base )  continue;

      if(pt->qual > crit) continue;

      ityp = MMG3D_typelt(mesh,k,item);
      cs[ityp]++;

      /*tet with bdry faces*/
      pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

      /*optim bdry tetra*/
      npeau = 0;
      if ( pxt ) {
        for(i=0 ; i<4 ; i++) {
          if ( pxt->ftag[i] & MG_BDY ) npeau++;
        }

        if(npeau>1) {
          nbdy++;
          continue;
        } else if ( npeau ) {
          nbdy2++;

          ier = MMG3D_optbdry(mesh,met,PROctree,k);
          if(ier) {
            nd++;
            ds[ityp]++;
            continue;
          }
        }
      }

      switch(ityp) {

      case 1:  /* sliver */
      case 3:  /* fin */
      case 6:  /* no good face: move away closest vertices */
      case 7:
      default:

        if(mesh->info.noswap) break;

        ier = MMG3D_swpItem(mesh,met,PROctree,k,item[0]);

        if(ier > 0) {
          nd++;
          ds[ityp]++;
          break;
        } else if(!ier) {
          /* second try to split the biggest edge */
          if(!mesh->info.noinsert) {
            /* if(MMG5_orvolnorm(mesh,k) < 5.e-9) { */
            /*   OCRIT *= 0.5; */
            /* } else */
            /*   OCRIT *= 0.75; */
            ier = MMG3D_splitItem(mesh,met,PROctree,k,item[0],1.01);

            if(ier) {
              nd++;
              ds[ityp]++;
              break;
            }
          } /* end noinsert */

          ier = MMG3D_swpalmostall(mesh,met,PROctree,k,item[0]);

          if(ier > 0) {
            nd++;
            ds[ityp]++;
            break;
          }

          if ( !mesh->info.noinsert ) {
            ier = MMG3D_splitalmostall(mesh,met,PROctree,k,item[0]);

            if(ier > 0) {
              nd++;
              ds[ityp]++;
              break;
            }
          }
        }
        if ( !mesh->info.nomove ) {
          for(i=0 ; i<4 ; i++) {
            if ( ((met->size!=1) && MMG3D_movv_ani(mesh,met,k,i)) ||
                 ((met->size==1) && MMG3D_movv_iso(mesh,met,k,i)) ) {
              nd++;
              ds[ityp]++;
              break;
            }
          }
        }
        break;
      case 2: /*chapeau*/
        if ( !mesh->info.nomove ) {
          if ( ( (met->size!=1) && MMG3D_movv_ani(mesh,met,k,item[0])) ||
               ((met->size==1) && MMG3D_movv_iso(mesh,met,k,item[0])) ) {
            nd++;
            ds[ityp]++;
          } else {
            for(i=0 ; i<4 ; i++) {
              if(item[0]==i) continue;
              if( ((met->size!=1) && MMG3D_movv_ani(mesh,met,k,i)) ||
                  ((met->size==1) && MMG3D_movv_iso(mesh,met,k,i)) ) {
                nd++;
                ds[ityp]++;
                break;
              }
            }
          }
        }
        break;
      } /* end switch */
    } /* end for k */

    /* printf("bdry : %" MMG5_PRId " %" MMG5_PRId "\n",nbdy,nbdy2); */
    /*  for (k=0; k<=7; k++) */
    /*    if ( cs[k] ) */
    /*    printf("  optim [%" MMG5_PRId "]      = %5d   %5d  %6.2f %%\n",k,cs[k],ds[k],100.0*ds[k]/cs[k]); */

    ntot += nd;
    if(base==-1) base = mesh->mark-1;
  } while (nd && it++<maxit);

  return ntot;
}
