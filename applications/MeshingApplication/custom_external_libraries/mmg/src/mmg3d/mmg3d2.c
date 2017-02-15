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
 * \file mmg3d/mmg3d2.c
 * \brief Create implicit surface in mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmg3d.h"

extern char  ddb;

/**
 * \remark Not used.
 *
 * solve 3*3 non symmetric system Ar = b
 *
 */
static inline int
_MMG5_invsl(double A[3][3],double b[3],double r[3]) {
  double detA;

  detA = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2]) \
    - A[0][1]*(A[1][0]*A[2][2] - A[2][0]*A[1][2]) \
    + A[0][2]*(A[1][0]*A[2][1] - A[2][0]*A[1][1]);
  if ( detA < _MMG5_EPSD )  return(0);
  detA = 1.0 / detA;

  r[0] =  b[0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2]) \
    - A[0][1]*(b[1]*A[2][2] - b[2]*A[1][2]) \
    + A[0][2]*(b[1]*A[2][1] - b[2]*A[1][1]);

  r[1] = A[0][0]*(b[1]*A[2][2] - b[2]*A[1][2]) \
    - b[0]*(A[1][0]*A[2][2] - A[2][0]*A[1][2]) \
    + A[0][2]*(A[1][0]*b[2] - A[2][0]*b[1]);

  r[2] = A[0][0]*(A[1][1]*b[2] - A[2][1]*b[1]) \
    - A[0][1]*(A[1][0]*b[2] - A[2][0]*b[1]) \
    + b[0]*(A[1][0]*A[2][1] - A[2][0]*A[1][1]);

  r[0] *= detA;
  r[1] *= detA;
  r[2] *= detA;

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \param k index of the starting tetra.
 * \param indp local index (inside the tria \a k) of the vertex that we check.
 * \return 1 if success, 0 if fail
 *
 * Check whether snapping the value of vertex \a indp to 0 exactly
 * leads to a non manifold situation.
 *
 * \warning: we assume that the triangle \a start has vertex \a istart
 * with value 0 and the other two with changing values.
 *
 */

static int
_MMG5_ismaniball(MMG5_pMesh mesh,MMG5_pSol sol,int k,int indp) {
  MMG5_pTetra   pt,pt1;
  double   v,v0,v1,v2;
  int      *adja,list[MMG3D_LMAX+1],bdy[MMG3D_LMAX+1],ibdy,np,ilist,base,cur,iel,jel,res,l;
  char     i,i0,i1,i2,j0,j1,j2,j,ip,nzeros,nopp,nsame;

  pt = &mesh->tetra[k];
  np = pt->v[indp];
  if ( fabs(sol->m[np]-mesh->info.ls) > _MMG5_EPSD2 )  return(1);

  memset(bdy,0,(MMG3D_LMAX+1)*sizeof(int));

  memset(list,0,(MMG3D_LMAX+1)*sizeof(int));

  /* Sign of a starting point in ball of np */
  for (j=0; j<3; j++) {
    ip = _MMG5_idir[indp][j];
    if ( sol->m[pt->v[ip]]-mesh->info.ls != 0.0 )  break;
  }
  if ( j == 3) {
    fprintf(stderr,"  *** Problem in function _MMG5_ismaniball : tetra %d : 4 null values",k);
    exit(EXIT_FAILURE);
  }

  v = sol->m[pt->v[ip]]-mesh->info.ls;
  base = ++mesh->base;
  pt->flag = base;
  ilist = 0;
  list[ilist] = 4*k+indp;
  ilist++;

  /* travel list and pile up, by adjacency, faces of ball of np while they have at least
     a vertex with same sign as v */
  res = cur = 0;
  while ( cur < ilist ) {
    iel = list[cur] / 4;
    i   = list[cur] % 4;
    pt  = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];

    /* Store a face for starting back enumeration with the opposite sign */
    if ( !res ) {
      for (j=0; j<3; j++) {
        i1 = _MMG5_idir[i][j];
        v1 = sol->m[pt->v[i1]]-mesh->info.ls;
        if ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) {
          res = 4*iel + i;
          break;
        }
      }
    }

    /* Pile up faces sharing a vertex with same sign as v */
    for (j=0; j<3; j++) {
      i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
      i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
      v1 = sol->m[pt->v[i1]]-mesh->info.ls;
      v2 = sol->m[pt->v[i2]]-mesh->info.ls;

      if ( ( ( v1 != 0.0 ) && MG_SMSGN(v,v1) ) ||
           ( ( v2 != 0.0 ) && MG_SMSGN(v,v2) ) ) {
        jel = adja[_MMG5_idir[i][j]];
        if( !jel ) continue;

        jel /= 4;
        pt1 = &mesh->tetra[jel];

        if ( pt1->flag == base )  continue;
        for (ip=0; ip<4; ip++) {
          if ( pt1->v[ip] == np )  break;
        }
        assert( ip < 4 );
        pt1->flag   = base;
        list[ilist] = 4*jel + ip;
        ilist++;
        assert(ilist < MMG3D_LMAX);
      }
    }
    cur++;
  }
  /* 0 value has been snapped accidentally */
  if ( !res ) {
    return(0);
  }


  /* Fill in list bdy, corresponding to the support tetras of the boundary to be created */
  ibdy = 0;
  for(l=0; l<ilist; l++) {
    iel = list[l] / 4;
    i   = list[l] % 4;
    pt  = &mesh->tetra[iel];

    nzeros = nsame = nopp = 0;

    i0 = _MMG5_idir[i][0];
    i1 = _MMG5_idir[i][1];
    i2 = _MMG5_idir[i][2];

    v0 = sol->m[pt->v[i0]]-mesh->info.ls;
    v1 = sol->m[pt->v[i1]]-mesh->info.ls;
    v2 = sol->m[pt->v[i2]]-mesh->info.ls;

    if ( v0 == 0.0 )
      nzeros++;
    else if ( MG_SMSGN(v,v0) )
      nsame++;
    else
      nopp++;

    if ( v1 == 0.0 )
      nzeros++;
    else if ( MG_SMSGN(v,v1) )
      nsame++;
    else
      nopp++;

    if ( v2 == 0.0 )
      nzeros++;
    else if ( MG_SMSGN(v,v2) )
      nsame++;
    else
      nopp++;

    if ( ( nzeros == 2 && nsame == 1 ) || ( nsame >= 1 && nopp >= 1 ) )  {
      bdy[ibdy] = list[l];
      ibdy++;
    }
  }

  /* Reset the current part of the ball, and start back the process with the other sign */
  iel = res / 4;
  pt = &mesh->tetra[iel];
  base = ++mesh->base;
  pt->flag = base;

  memset(list,0,(MMG3D_LMAX+1)*sizeof(int));
  ilist = cur = 0;
  list[ilist] = res;
  ilist++;
  while ( cur < ilist ) {
    iel = list[cur] / 4;
    i   = list[cur] % 4;
    pt  = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];

    /* Pile up faces sharing a vertex with opposite sign to v */
    for (j=0; j<3; j++) {
      i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
      i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
      v1 = sol->m[pt->v[i1]]-mesh->info.ls;
      v2 = sol->m[pt->v[i2]]-mesh->info.ls;

      if ( v1 == 0.0 && v2 == 0.0 ) {
        jel = adja[_MMG5_idir[i][j]];
        if( !jel ) continue;
        jel /=4 ;
        pt1 = &mesh->tetra[jel];
        pt1->flag = base;
      }

      else if ( ( ( v1 != 0.0 ) && (!MG_SMSGN(v,v1)) ) || ( ( v2 != 0.0 ) && (!MG_SMSGN(v,v2)) ) ) {
        jel = adja[_MMG5_idir[i][j]];
        if( !jel ) continue;
        jel /= 4;
        pt1 = &mesh->tetra[jel];

        j0 = _MMG5_idir[i][0];
        j1 = _MMG5_idir[i][1];
        j2 = _MMG5_idir[i][2];

        v0 = sol->m[pt1->v[j0]]-mesh->info.ls;
        v1 = sol->m[pt1->v[j1]]-mesh->info.ls;
        v2 = sol->m[pt1->v[j2]]-mesh->info.ls;

        nzeros = nsame = nopp = 0;

        if ( v0 == 0.0 )
          nzeros++;
        else if ( MG_SMSGN(v,v0) )
          nsame++;
        else
          nopp++;

        if ( v1 == 0.0 )
          nzeros++;
        else if ( MG_SMSGN(v,v1) )
          nsame++;
        else
          nopp++;

        if ( v2 == 0.0 )
          nzeros++;
        else if ( MG_SMSGN(v,v2) )
          nsame++;
        else
          nopp++;

        if ( ( nzeros == 2 && nsame == 1 ) || ( nsame >= 1 && nopp >= 1 ) )  {
          if ( pt1->flag < base - 1 ) return(0);
        }

        if ( pt1->flag == base ) continue;
        for (ip=0; ip<4; ip++) {
          if ( pt1->v[ip] == np )  break;
        }
        assert( ip < 4 );
        pt1->flag   = base;
        list[ilist] = 4*jel + ip;
        ilist++;
        assert(ilist < MMG3D_LMAX);
      }
    }
    cur++;
  }

  /* Now, all elements of bdy should have been marked by a flag base + 1 */
  for (l=0; l<ibdy; l++) {
    iel = bdy[l] / 4;
    pt = &mesh->tetra[iel];
    if ( pt->flag != base ) return(0);
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set function.
 * \param tmp saving of the level-set values before the snap.
 * \return 1 if success, 0 if fail.
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
static int _MMG3D_snpval_ls(MMG5_pMesh mesh,MMG5_pSol sol,double *tmp) {
  MMG5_pTetra   pt;
  MMG5_pPoint   p0;
  int      k,nc,ns,ip;
  char     i;

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"  ## Hashing problem (1). Exit program.\n");
    return(0);
  }

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]-mesh->info.ls) < _MMG5_EPS ) {
      if ( mesh->info.ddebug )  fprintf(stdout,"  Snapping value %d ; previous value : %E\n",k,fabs(sol->m[k]));
      tmp[k] = ( fabs(sol->m[k]-mesh->info.ls) < _MMG5_EPSD ) ? (mesh->info.ls-100.0*_MMG5_EPS) : sol->m[k];
      p0->flag = 1;
      sol->m[k] = mesh->info.ls;
      ns++;
    }
  }

  /* Check snapping did not lead to a nonmanifold situation */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<4; i++) {
      ip = pt->v[i];
      p0 = &mesh->point[ip];
      if ( p0->flag ) {
        if ( !_MMG5_ismaniball(mesh,sol,k,i) ) {
          sol->m[ip] = tmp[ip];
          nc++;
        }
        p0->flag = 0;
        tmp[ip]  = mesh->info.ls;
      }
    }
  }
  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+nc > 0 )
    fprintf(stdout,"     %8d points snapped, %d corrected\n",ns,nc);

  /* memory free */
  _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \return 1 if success, 0 otherwise.
 *
 * Proceed to discretization of the implicit function carried by sol into mesh,
 * once values of sol have been snapped/checked
 *
 */
static int _MMG3D_cuttet_ls(MMG5_pMesh mesh, MMG5_pSol sol/*,double *tmp*/){
  MMG5_pTetra   pt;
  MMG5_pPoint   p0,p1;
  _MMG5_Hash     hash;
  double   c[3],v0,v1,s;
  int      vx[6],nb,k,ip0,ip1,np,ns,ne;
  char     ia;
  /* Commented because unused */
  /*MMG5_pPoint  p[4];*/
  /*double   *grad,A[3][3],b[3],*g0,*g1,area,a,d,dd,s1,s2;*/
  /*int       ip[4],ng*/
  /*char    i,ier;*/

  /* reset point flags and h */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[_MMG5_iare[ia][0]];
      ip1 = pt->v[_MMG5_iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      v0  = sol->m[ip0]-mesh->info.ls;
      v1  = sol->m[ip1]-mesh->info.ls;
      if ( fabs(v0) > _MMG5_EPSD2 && fabs(v1) > _MMG5_EPSD2 && v0*v1 < 0.0 ) {
        if ( !p0->flag ) {
          p0->flag = nb;
          nb++;
        }
        if ( !p1->flag ) {
          p1->flag = nb;
          nb++;
        }
      }
    }
  }
  if ( ! nb )  return(1);

  /* Store gradients of level set function at those points */
  /* Commented because unused */
  /* grad = (double*)calloc(3*nb+1,sizeof(double)); */
  /* assert(grad); */

  /* for (k=1; k<=mesh->ne; k++) { */
  /*   pt = &mesh->tetra[k]; */
  /*   ia = 0; */
  /*   for (i=0; i<4; i++) { */
  /*     ip[i] = pt->v[i]; */
  /*     p[i]  = &mesh->point[ip[i]]; */
  /*     if ( p[i]->flag == 0 )  ia++; */
  /*   } */
  /*   if ( ia == 4 )  continue; */

  /*   A[0][0] = p[1]->c[0] - p[0]->c[0];  A[0][1] = p[1]->c[1] - p[0]->c[1];  A[0][2] = p[1]->c[2] - p[0]->c[2]; */
  /*   A[1][0] = p[2]->c[0] - p[0]->c[0];  A[1][1] = p[2]->c[1] - p[0]->c[1];  A[1][2] = p[2]->c[2] - p[0]->c[2]; */
  /*   A[2][0] = p[3]->c[0] - p[0]->c[0];  A[2][1] = p[3]->c[1] - p[0]->c[1];  A[2][2] = p[3]->c[2] - p[0]->c[2]; */

  /*   b[0] = sol->m[ip[1]] - sol->m[ip[0]]; */
  /*   b[1] = sol->m[ip[2]] - sol->m[ip[0]]; */
  /*   b[2] = sol->m[ip[3]] - sol->m[ip[0]]; */

  /*   area = _MMG5_det4pt(p[0]->c,p[1]->c,p[2]->c,p[3]->c); */
  /*   ier  = _MMG5_invsl(A,b,c); */
  /*   if ( !ier )  continue; */

  /*   for (i=0; i<4; i++) { */
  /*     if ( p[i]->flag ) { */
  /*       ng = p[i]->flag; */
  /*       tmp[ip[i]] += fabs(area); */
  /*       grad[3*(ng-1)+1] += (area*c[0]); */
  /*       grad[3*(ng-1)+2] += (area*c[1]); */
  /*       grad[3*(ng-1)+3] += (area*c[2]); */
  /*     } */
  /*   } */
  /* } */
  /* for (k=1; k<=mesh->np; k++) { */
  /*   p0 = &mesh->point[k]; */
  /*   if ( p0->flag ) { */
  /*     area = MG_MAX(_MMG5_EPSD2,tmp[k]); */
  /*     area = 1.0 / area; */
  /*     ng   = p0->flag; */
  /*     grad[3*(ng-1)+1] *= area; */
  /*     grad[3*(ng-1)+2] *= area; */
  /*     grad[3*(ng-1)+3] *= area; */
  /*   } */
  /* } */

  /* Create intersection points at 0 isovalue and set flags to tetras */
  if ( !_MMG5_hashNew(mesh,&hash,nb,7*nb) ) return(0);
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[_MMG5_iare[ia][0]];
      ip1 = pt->v[_MMG5_iare[ia][1]];
      np  = _MMG5_hashGet(&hash,ip0,ip1);
      if ( np )  continue;

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0]-mesh->info.ls;
      v1 = sol->m[ip1]-mesh->info.ls;
      if ( fabs(v0) < _MMG5_EPSD2 || fabs(v1) < _MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      /* g0 = &grad[3*(p0->flag -1)+1]; */
      /* g1 = &grad[3*(p1->flag -1)+1]; */
      /* a = 0.5 * ((g1[0]-g0[0])*(p1->c[0]-p0->c[0]) + (g1[1]-g0[1])*(p1->c[1]-p0->c[1]) \ */
      /*            + (g1[2]-g0[2])*(p1->c[2]-p0->c[2])); */
      /* d  = v1 - v0 - a; */
      /* dd = d*d - 4.0*a*v0; */
      /* dd = MG_MAX(_MMG5_EPSD2,dd); */
      /* dd = sqrt(dd); */
      /* if ( fabs(a) < _MMG5_EPSD2 ) */
      /*   s = v0 / (v0-v1); */
      /* else { */
      /*   s1 = 0.5*( dd -d) / a; */
      /*   s2 = 0.5*(-dd -d) / a; */
      /*   if ( s1 > 0.0 && s1 < 1.0 ) */
      /*     s = s1; */
      /*   else if (s2 > 0.0 && s2 < 1.0) */
      /*     s = s2; */
      /*   else */
      /*     s = MG_MIN(fabs(s1),fabs(s1-1.0)) < MG_MIN(fabs(s2),fabs(s2-1.0)) ? s1 : s2 ; */
      /* } */
      // IMPORTANT A REGARDER
      s = v0 / (v0-v1);

      s = MG_MAX(MG_MIN(s,1.0-_MMG5_EPS),_MMG5_EPS);
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

      np = _MMG3D_newPt(mesh,c,0);
      if ( !np ) {
        _MMG5_POINT_REALLOC(mesh,sol,np,0.2,
                            fprintf(stderr,"  ## Error: unable to allocate a new point\n");
                            _MMG5_INCREASE_MEM_MESSAGE();
                            return(0)
                            ,c,0);
      }
      sol->m[np] = mesh->info.ls;
      _MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }

  /* Proceed to splitting, according to flags to tets */
  ne = mesh->ne;
  ns = 0;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ) )  continue;
    pt->flag = 0;
    memset(vx,0,6*sizeof(int));
    for (ia=0; ia<6; ia++) {
      vx[ia] = _MMG5_hashGet(&hash,pt->v[_MMG5_iare[ia][0]],pt->v[_MMG5_iare[ia][1]]);
      if ( vx[ia] )  MG_SET(pt->flag,ia);
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      _MMG5_split1(mesh,sol,k,vx,1);
      ns++;
      break;

    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      _MMG5_split2sf(mesh,sol,k,vx,1);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      _MMG5_split3cone(mesh,sol,k,vx,1);
      ns++;
      break;

    case 30: case 45: case 51:
      _MMG5_split4op(mesh,sol,k,vx,1);
      ns++;
      break;

    default :
      assert(pt->flag == 0);
      break;
    }
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",ns);

  _MMG5_DEL_MEM(mesh,hash.item,(hash.max+1)*sizeof(_MMG5_hedge));
  return(ns);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \return 1.
 *
 * Set references to tets according to the sign of the level set function.
 *
 */
static int _MMG3D_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTetra   pt;
  double        v;
  int      k,ip;
  char     nmns,npls,nz,i;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    nmns = npls = nz = 0;
    for (i=0; i<4; i++) {
      ip = pt->v[i];
      v  = sol->m[ip]-mesh->info.ls;
      if ( v > 0.0 )
        npls++;
      else if ( v < 0.0 )
        nmns++;
      else
        nz ++;
    }
    assert(nz < 4);
    if ( npls ) {
      assert(!nmns);
      pt->ref = MG_PLUS;
    }
    else {
      assert(nmns);
      pt->ref = MG_MINUS;
    }
  }
  return(1);
}

/** Check whether implicit surface is orientable in ball of point ip in tet iel ;
    Beware : may return 0 when implicit boundary is tangent to outer boundary */
int _MMG5_chkmaniball(MMG5_pMesh mesh, int start, char ip){
  MMG5_pTetra    pt,pt1;
  int       ref,base,ilist,nump,k,cur,k1,nref;
  int       *adja,list[MMG3D_LMAX+2];
  char      i,l,j;

  base = ++mesh->base;
  ilist = 0;

  pt = &mesh->tetra[start];
  nump = pt->v[ip];
  ref = pt->ref;

  /* Store initial tetrahedron */
  pt->flag = base;
  list[ilist] = 4*start+ip;
  ilist++;

  /* explore list, and find all tets in ball of p belonging to the component ref */
  cur = 0;
  while( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4;
    pt = &mesh->tetra[k];

    adja = &mesh->adja[4*(k-1)+1];
    for(l=0; l<3; l++){
      i = _MMG5_inxt3[i];

      /* Travel only through non boundary faces. */
      k1 = adja[i];
      if(!k1) continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];

      if( pt1 ->ref != ref ) continue;

      if( pt1->flag == base ) continue;
      pt1->flag = base;

      for(j=0; j<4 ; j++){
        if(pt1->v[j] == nump)
          break;
      }
      assert(j<4);

      /* overflow */
      assert ( ilist <= MMG3D_LMAX-3 );
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }

  /* Number of caught tets with ref ptstart->ref*/
  nref = ilist;

  /* Complete ball of point */
  cur = 0;
  while(cur < ilist){
    k = list[cur] / 4;
    i = list[cur] % 4;
    pt = &mesh->tetra[k];

    adja = &mesh->adja[4*(k-1)+1];
    for(l=0; l<3; l++){
      i = _MMG5_inxt3[i];

      k1 = adja[i];
      if ( !k1 ) continue;
      k1/=4;

      pt1 = &mesh->tetra[k1];
      if(pt1->flag == base) continue;
      pt1->flag = base;

      for(j=0; j<4 ; j++){
        if(pt1->v[j] == nump)
          break;
      }
      assert(j<4);

      /* overflow */
      assert ( ilist <= MMG3D_LMAX-3 );
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }

  /* Elements from nref to ilist-1 must not have ref ptstart->ref */
  for(cur=nref; cur<ilist; cur++) {
    k = list[cur] / 4;
    pt = &mesh->tetra[k];
    if( pt->ref == ref ) {
      fprintf(stderr,"   *** Topological problem:");
      fprintf(stderr," non manifold surface at point %d \n",nump);
      return(0);
    }
  }

  return(1);
}

/** Check whether implicit surface enclosed in volume is orientable */
int _MMG5_chkmani(MMG5_pMesh mesh){
  MMG5_pTetra    pt,pt1;
  int       k,iel,ref;
  int       *adja;
  char      i,j,ip,cnt;

  for(k=1; k<=mesh->np; k++){
    mesh->point[k].flag = 0;
  }

  /** First test : check whether a tetra has 4 boundary faces */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )   continue;
    adja = &mesh->adja[4*(k-1)+1];

    ref = pt->ref;
    cnt = 0;
    for(i=0; i<4; i++) {
      if( !adja[i] ) {
        cnt++;
      }
      else {
        pt1 = &mesh->tetra[adja[i]/4];
        if ( pt1->ref != ref ) cnt++;
      }
    }
    if ( cnt == 4 ) {
      fprintf(stdout,"Tetra %d : 4 boundary faces \n",k);
      //return(0);
    }
  }

  /** Second test : Check whether configuration is manifold in each ball */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ))   continue;
    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      if(!adja[i]) continue;
      iel = adja[i] / 4;
      pt1 = &mesh->tetra[iel];
      if(pt1->ref == pt->ref) continue;

      for(j=0; j<3; j++){
        ip = _MMG5_idir[i][j];

        if(!_MMG5_chkmaniball(mesh,k,ip))
          return(0);
      }
    }
  }

  if ( mesh->info.imprim || mesh->info.ddebug )
    fprintf(stdout,"  *** Manifold implicit surface.\n");
  return(1);
}

/** Check whether implicit surface enclosed in volume is orientable (perform an additionnal
    test w.r.t. _MMG5_chkmani) */
int _MMG5_chkmani2(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTetra    pt,pt1;
  int       k,iel;
  int       *adja;
  char      i,j,ip,cnt;

  for(k=1; k<=mesh->np; k++){
    mesh->point[k].flag = 0;
  }

  /** First test : assure no tetra has its 4 vertices on implicit boundary */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ))   continue;

    cnt = 0;
    for(j=0; j<4; j++) {
      if( sol->m[pt->v[j]]-mesh->info.ls == 0.0 ) cnt++;
    }
    if(cnt == 4) {
      fprintf(stderr,"Problem in tetra %d : 4 vertices on implicit boundary",k);
      exit(EXIT_FAILURE);
    }
  }

  /** Second test : check whether configuration is manifold in each ball */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ))   continue;
    adja = &mesh->adja[4*(k-1)+1];

    for(i=0; i<4; i++){
      if(!adja[i]) continue;
      iel = adja[i] / 4;
      pt1 = &mesh->tetra[iel];
      if(pt1->ref == pt->ref) continue;

      for(j=0; j<3; j++){
        ip = _MMG5_idir[i][j];

        if(!_MMG5_chkmaniball(mesh,k,ip)){
          fprintf(stderr,"Non orientable implicit surface : ball of point %d\n",pt->v[ip]);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if ( mesh->info.ddebug )  fprintf(stdout,"  *** Manifold implicit surface.\n");
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param k index of element in which we collapse.
 * \param iface face through wich we perform the collapse
 * \param iedg edge to collapse
 * \param ndepmin index of an elt with ref MG_MINUS and outside the shell of edge.
 * \param ndepplus ndex of an elt with ref MG_PLUS and outside the shell of edge.
 * \param isminp 1 if we have found a tetra with ref MG_MINUS
 * \param isplp 1 if we have found a tetra with ref MG_PLUS
 * \return 0 if we create a non manifold situation, 1 otherwise
 *
 * Check whether collapse of point np to nq does not create a non manifold
 * situation at nq ndepmin, ndepplus = tetra of ref minus, plus in ball of np,
 * not in shell of (np,nq).
 *
 */
int _MMG5_chkmanicoll(MMG5_pMesh mesh,int k,int iface,int iedg,int ndepmin,int ndepplus,char isminp,char isplp) {
  MMG5_pTetra    pt,pt1;
  int       nump,numq,ilist,ref,cur,stor,iel,jel,base,ndepmq,ndeppq;
  int       list[MMG3D_LMAX+2],*adja,*adja1;
  char      i,j,ip,jp,iq,jq,voy,indp,indq,isminq,isplq,ismin,ispl;

  ilist = 0;
  ndepmq = ndeppq = 0;
  isplq = isminq = 0;

  pt    = &mesh->tetra[k];
  ip    = _MMG5_idir[iface][_MMG5_inxt2[iedg]];
  iq    = _MMG5_idir[iface][_MMG5_iprv2[iedg]];
  nump  = pt->v[ip];
  numq  = pt->v[iq];

  /* Case when nump does not have any interior (resp. ext.) tetra which will not
     disappear : search for start in ball of q */
  if ( !ndepmin || !ndepplus ) {
    base  = ++mesh->base;

    pt = &mesh->tetra[k];
    for(j=0; j<4; j++) {
      if ( pt->v[j] == numq ) break;
    }
    assert( j < 4 );
    list[ilist] = 4*k+j;
    ilist++;
    assert( ilist < MMG3D_LMAX+2 );
    pt->flag = base;

    if ( pt->ref == MG_MINUS ) isminq = 1;
    else if ( pt->ref == MG_PLUS ) isplq = 1;

    cur = 0;
    while( cur < ilist ) {
      iel = list[cur] / 4;
      i = list[cur] % 4;
      adja = &mesh->adja[4*(iel-1)+1];

      for (j=0; j<3; j++) {
        i = _MMG5_inxt3[i];
        jel = adja[i];
        if ( !jel ) continue;

        jel /= 4;
        pt1 = &mesh->tetra[jel];

        if ( pt1->ref == MG_MINUS ) isminq = 1;
        else if ( pt1->ref == MG_PLUS ) isplq = 1;

        if ( pt1->flag == base ) continue;
        pt1->flag = base;
        for(iq=0; iq<4; iq++)
          if ( pt1->v[iq] == numq ) break;
        assert( iq < 4 );
        /* overflow */
        assert ( ilist < MMG3D_LMAX+2 );

        list[ilist] = 4*jel+iq;
        ilist++;

        /* check if jel is an available starting tetra for further enumeration */
        if ( !ndeppq && pt1->ref == MG_PLUS ) {
          for(ip=0; ip<4; ip++)
            if ( pt1->v[ip] == nump ) break;
          if( ip == 4 ) ndeppq = jel;
        }
        if( !ndepmq && pt1->ref == MG_MINUS ) {
          for(ip=0; ip<4; ip++)
            if ( pt1->v[ip] == nump ) break;
          if( ip == 4 ) ndepmq = jel;
        }
      }
      cur++;
    }

    memset(list,0,(MMG3D_LMAX+2)*sizeof(int));
    ilist = 0;
  }

  ispl = ( isplp || isplq ) ? 1 : 0;
  ismin = ( isminp || isminq ) ? 1 : 0;

  /** First step : pile up tetras of future ball of nq, crossing
      through the shell of (np,nq), as long as they have same ref as ndepmin
      list[l] <= 0 if element of ball of np, >= 0, if element of ball of nq */
  base  = ++mesh->base;

  if( ndepmin ) {
    pt = &mesh->tetra[ndepmin];
    ref = pt->ref;

    for(j=0; j<4; j++) {
      if ( pt->v[j] == nump ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = - (4*ndepmin+j);
    ilist++;
  }
  else if ( ndepmq ) {
    pt = &mesh->tetra[ndepmq];
    ref = pt->ref;

    for(j=0; j<4; j++) {
      if ( pt->v[j] == numq ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = 4*ndepmq+j;
    ilist++;
  }
  else {
    if ( ismin && ispl )
      return(0);
    else
      return(1);
  }


  cur = 0;
  while ( cur < ilist ) {
    stor = list[cur];
    /* Element belongs to the ball of np */
    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;
      for (i=0; i<3; i++) {
        jp = _MMG5_inxt3[jp];
        jel = adja[jp];
        if ( !jel ) continue;

        jel /= 4;
        voy = adja[jp] % 4;

        pt1 = &mesh->tetra[jel];
        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == numq ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          jel = adja1[j];
          if (!jel ) continue;

          jel /= 4;
          pt1 = &mesh->tetra[jel];

          if ( pt1->ref != ref) continue;   // ICI, il ne faut pas autoriser à passer si on a à nouveau un tet de la coquille (avant de marquer)
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == nump ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;

          for(j=0; j<4; j++)
            if ( pt1->v[j] == numq ) break;
          assert( j< 4);

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if ( pt1->v[j] == nump ) break;
          assert( j< 4 );

          list[ilist] = - (4*jel+j);
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
      }
    }
    /* Element belongs to the ball of nq */
    else {
      iel = stor / 4;
      iq  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;
      for (i=0; i<3; i++) {
        jq = _MMG5_inxt3[jq];
        jel = adja[jq];
        if ( !jel ) continue;

        jel /= 4;
        voy = adja[jq] % 4;

        pt1 = &mesh->tetra[jel];
        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == nump ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4);

          jel = adja1[j];
          if (!jel ) continue;

          jel /= 4;

          pt1 = &mesh->tetra[jel];
          if ( pt1->ref != ref) continue;
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == numq ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          list[ilist] = -(4*jel+j);
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4);

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
      }
    }
    cur++;
  }

  assert( cur == ilist );

  /** Second step : same process, starting with a tetra of different reference, in the ball of np */
  if( ndepplus ) {
    pt = &mesh->tetra[ndepplus];
    for(j=0; j<4; j++) {
      if ( pt->v[j] == nump ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = - (4*ndepplus+j);
    ilist++;
    ref = pt->ref;
  }
  else if ( ndeppq ) {
    pt = &mesh->tetra[ndeppq];
    for(j=0; j<4; j++) {
      if ( pt->v[j] == numq ) break;
    }
    assert( j < 4 );

    pt->flag = base;
    list[ilist] = 4*ndeppq+j;
    ilist++;
    ref = pt->ref;
  }
  else {
    if ( ismin && ispl )
      return(0);
    else
      return(1);
  }

  while ( cur < ilist ) {
    stor = list[cur];
    /* Element belongs to the ball of np */
    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;

      for (i=0; i<3; i++) {
        jp = _MMG5_inxt3[jp];
        jel = adja[jp];

        if ( !jel ) continue;

        jel /=4;
        voy = adja[jp] % 4;

        pt1 = &mesh->tetra[jel];
        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == numq ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          jel = adja1[j];
          if (!jel ) continue;

          jel /= 4;
          pt1 = &mesh->tetra[jel];
          if ( pt1->ref != ref) continue;
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == nump ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;
          for(j=0; j<4; j++)
            if ( pt1->v[j] == numq ) break;
          assert( j< 4 );

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          list[ilist] = - (4*jel+j);
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
      }
    }
    /* Element belongs to the ball of nq */
    else {
      iel = stor / 4;
      iq  = stor % 4;

      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;

      for (i=0; i<3; i++) {
        jq = _MMG5_inxt3[jq];
        jel = adja[jq];

        if ( !jel ) continue;

        jel /= 4;
        voy = adja[jq] % 4;

        pt1 = &mesh->tetra[jel];

        if ( pt1->ref != ref ) continue;

        /* Current tetra is neighbour of a tetra of the shell of (np,nq) */
        if( pt1->v[voy] == nump ) {
          adja1 = &mesh->adja[4*(jel-1)+1];
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4);

          jel = adja1[j];
          if (!jel ) continue;

          jel /= 4;
          pt1 = &mesh->tetra[jel];
          if ( pt1->ref != ref) continue;
          if ( pt1->flag == base ) continue;

          /* New tetra to be added must not be itself an element of the shell */
          for(j=0; j<4; j++) {
            if ( pt1->v[j] == numq ) break;
          }
          if ( j<4 ) continue;

          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == nump ) break;
          assert( j< 4);

          list[ilist] = -(4*jel+j);
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
        else {
          if ( pt1->flag == base ) continue;
          pt1->flag = base;
          for(j=0; j<4; j++)
            if (pt1->v[j] == numq ) break;
          assert( j< 4 );

          list[ilist] = 4*jel+j;
          ilist++;
          assert( ilist < MMG3D_LMAX+1 );
        }
      }
    }
    cur++;
  }
  assert( cur == ilist );

  /* At this point, all elements of ball np \cup ball nq \setminus shell have been tagged
     unless the future ball of nq, ending up from collapse is non manifold */
  cur = 0;
  while ( cur < ilist ) {
    stor = list[cur];

    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;
      for(i=0; i<3; i++) {
        jp = _MMG5_inxt3[jp];
        jel = adja[jp];

        if ( !jel ) continue;

        jel /= 4;
        pt1 = &mesh->tetra[jel];
        if (pt1->flag == base ) continue;
        pt1->flag = base;

        indp = -1;
        indq = -1;
        for(j=0; j<4; j++) {
          if ( pt1->v[j] == nump )
            indp = j;
          else if ( pt1->v[j] == numq )
            indq = j;
        }
        assert( indp >= 0 && indp < 4 );

        /* Only tets of the shell of (np,nq) can be added, unless future ball is non manifold */
        if ( indq == -1 ) {
          fprintf(stdout,"  ## Warning: we should rarely passed here. ");
          fprintf(stdout,"tetra %d =  %d %d %d %d, ref = %d\n",
                  jel,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3],pt1->ref);
          return(0);
        }

        list[ilist] = -(4*jel+indp);
        ilist++;
        assert( ilist < MMG3D_LMAX +1 );
      }
    }
    else {
      iel = stor / 4;
      iq  = stor % 4;
      pt = &mesh->tetra[iel];
      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;
      for(i=0; i<3; i++) {
        jq = _MMG5_inxt3[jq];
        jel = adja[jq];

        if ( !jel ) continue;

        jel /= 4;
        pt1 = &mesh->tetra[jel];
        if (pt1->flag == base ) continue;
        pt1->flag = base;

        indp = -1;
        indq = -1;

        for(j=0; j<4; j++) {
          if ( pt1->v[j] == nump )
            indp = j;
          else if ( pt1->v[j] == numq )
            indq = j;
        }
        assert( indq >= 0 && indq < 4 );

        /* Only tets of the shell of (np,nq) can be added, unless future ball is non manifold */
        if ( indp == -1 ) {
          fprintf(stdout,"  ## Warning: we should rarely passed here. ");
          fprintf(stdout,"tetra %d =  %d %d %d %d, ref = %d\n",
                  jel,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3],pt1->ref);
          return(0);
        }

        list[ilist] = 4*jel+indq;
        ilist++;
        assert( ilist < MMG3D_LMAX +1 );
      }
    }
    cur++;
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int _MMG3D_mmg3d2(MMG5_pMesh mesh,MMG5_pSol sol) {
  double   *tmp;

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION\n");

  if ( mesh->nprism || mesh->nquad ) {
    fprintf(stderr,"  ## Error: Isosurface extraction not available with hybrid"
            " meshes. Exit program.\n");
    return(0);
  }

  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                fprintf(stderr,"  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double);

  /* Snap values of level set function if need be, then discretize it */
  if ( !_MMG3D_snpval_ls(mesh,sol,tmp) ) {
    fprintf(stderr,"  ## Problem with implicit function. Exit program.\n");
    return(0);
  }
  _MMG5_DEL_MEM(mesh,tmp,(mesh->npmax+1)*sizeof(double));

  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"  ## Hashing problem. Exit program.\n");
    return(0);
  }

  /* compatibility triangle orientation w/r tetras */
  if ( !_MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"  ## Boundary orientation problem. Exit program.\n");
    return(0);
  }

  if ( !_MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"  ## Boundary problem. Exit program.\n");
    return(0);
  }

  /* build hash table for initial edges */
  if ( !_MMG5_hGeom(mesh) ) {
    fprintf(stderr,"  ## Hashing problem (0). Exit program.\n");
    return(0);
  }

  if ( !_MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"  ## Problem in setting boundary. Exit program.\n");
    return(0);
  }

  if ( !_MMG3D_cuttet_ls(mesh,sol/*,tmp*/) ) {
    fprintf(stderr,"  ## Problem in discretizing implicit function. Exit program.\n");
    return(0);
  }

  _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
  _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));

  mesh->nt = 0;

  if ( !_MMG3D_setref_ls(mesh,sol) ) {
    fprintf(stderr,"  ## Problem in setting references. Exit program.\n");
    return(0);
  }

  /* Clean memory */
  _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));

  return(1);
}
