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

#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "mmg3dexterns_private.h"

extern int8_t ddb;

/**
 * \param ppt array of points containing the tetra vertices
 * \param i0 local index of the vertex that has a sign different to the other vertices.
 * \param part_opp 0 if we want to compute the area containing the vertex \a i0,
 *   1 if we want the area that do not contains \a i0.
 *
 * \return the computed area (multiplied by 6) if sucess or 0.0 if fail.
 *
 * Compute the area (x6) defined by the level-set inside the tetra with vertices
 * \a ppt. This tetra must be splitted by the level-set such has it has exactly
 * 1 vertex (the vertex \a i0) with sign opposite to the other vertices. If \a
 * part_opp == 0, we compte the area that contains \a i0, otherwise we compute
 * the complementary area.
 *
 */
static inline
double MMG3D_vfrac_1vertex(MMG5_pPoint ppt[4],int8_t i0,double v[4],int8_t part_opp) {
  double      vfrac,lam,area,o1[3],o2[3],o3[3];
  int8_t      i1,i2,i3;

  i1 = MMG5_idir[i0][0];
  i2 = MMG5_idir[i0][1];
  i3 = MMG5_idir[i0][2];

  lam = v[i0] / (v[i0]-v[i1]);
  o1[0] = ppt[i0]->c[0] + lam*(ppt[i1]->c[0]-ppt[i0]->c[0]);
  o1[1] = ppt[i0]->c[1] + lam*(ppt[i1]->c[1]-ppt[i0]->c[1]);
  o1[2] = ppt[i0]->c[2] + lam*(ppt[i1]->c[2]-ppt[i0]->c[2]);

  lam = v[i0] / (v[i0]-v[i2]);
  o2[0] = ppt[i0]->c[0] + lam*(ppt[i2]->c[0]-ppt[i0]->c[0]);
  o2[1] = ppt[i0]->c[1] + lam*(ppt[i2]->c[1]-ppt[i0]->c[1]);
  o2[2] = ppt[i0]->c[2] + lam*(ppt[i2]->c[2]-ppt[i0]->c[2]);

  lam = v[i0] / (v[i0]-v[i3]);
  o3[0] = ppt[i0]->c[0] + lam*(ppt[i3]->c[0]-ppt[i0]->c[0]);
  o3[1] = ppt[i0]->c[1] + lam*(ppt[i3]->c[1]-ppt[i0]->c[1]);
  o3[2] = ppt[i0]->c[2] + lam*(ppt[i3]->c[2]-ppt[i0]->c[2]);

  vfrac = fabs (MMG5_det4pt(ppt[i0]->c,o1,o2,o3));

  if ( !part_opp ) {
    return vfrac;
  }
  else {
    area = MMG5_det4pt(ppt[0]->c,ppt[1]->c,ppt[2]->c,ppt[3]->c);
    vfrac  = fabs(area) - vfrac;
    return vfrac;
  }

  /* Should not pass here */
  return 0.0;
}


/**
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the ls function
 * \param k index of the triangle
 * \return volfrac
 *
 * Calculate the area of the positive (if pm == 1) or negative (if pm == -1)
 * subdomain inside tetra k defined by the ls function in sol
 *
 **/
double MMG3D_vfrac(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int k,int pm) {
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt[4];
  double         v[4],vfm,vfp,lam,eps,o[18];
  int            nplus,nminus,nzero;
  MMG5_int       flag,ip[4];
  int8_t         cfg,ia;
  int8_t         i,i0,i1,imin1,imin2,iplus1,iplus2,iz;
  uint8_t        tau[4];
  const uint8_t  *taued;

  eps = MMG5_EPS*MMG5_EPS;
  pt = &mesh->tetra[k];

  ip[0] = pt->v[0];
  ip[1] = pt->v[1];
  ip[2] = pt->v[2];
  ip[3] = pt->v[3];

  ppt[0] = &mesh->point[ip[0]];
  ppt[1] = &mesh->point[ip[1]];
  ppt[2] = &mesh->point[ip[2]];
  ppt[3] = &mesh->point[ip[3]];

  v[0] = sol->m[ip[0]];
  v[1] = sol->m[ip[1]];
  v[2] = sol->m[ip[2]];
  v[3] = sol->m[ip[3]];

  /* Identify number of zero, positive and negative vertices, and corresponding
   * indices */
  nplus = nminus = nzero = 0;
  imin1 = imin2 = iplus1 = iplus2 = iz = -1;

  for (i=0; i<4; i++) {
    if ( fabs(v[i]) < eps ) {
      nzero++;
      if ( iz < 0 ) iz = i;
    }
    else if ( v[i] >= eps ) {
      nplus++;
      if ( iplus1 < 0 ) iplus1 = i;
      else if ( iplus2 < 0 ) iplus2 = i;
    }
    else {
      nminus++;
      if ( imin1 < 0 ) imin1 = i;
      else if ( imin2 < 0 ) imin2 = i;
    }
  }

  /* Degenerate case */
  if ( nzero == 4 ) return 0.0;

  /* Whole tetra is positive */
  if ( nminus == 0 ) {
    vfp = MMG5_det4pt(ppt[0]->c,ppt[1]->c,ppt[2]->c,ppt[3]->c);
    vfp = fabs(vfp);
    if ( pm == 1 ) return vfp;
    else           return 0.0;
  }

  /* Whole tetra is negative */
  if ( nplus == 0 ) {
    vfm = MMG5_det4pt(ppt[0]->c,ppt[1]->c,ppt[2]->c,ppt[3]->c);
    vfm = fabs(vfm);
    if ( pm == -1 ) return vfm;
    else            return 0.0;
  }

  /* Exactly one vertex is negative */
  if ( nminus == 1 ) {
    return MMG3D_vfrac_1vertex(ppt,imin1,v,pm!=-1);
  }

  /* Exactly one vertex is positive */
  if ( nplus == 1 ) {
    return MMG3D_vfrac_1vertex(ppt,iplus1,v,pm!=1);
  }

  flag = 0;
  for ( ia=0; ia<18; ++ia ) {
    o[ia] = 0.0;
  }

  /* We have exactly 2 negative vertices and 2 positive ones */
  assert ( nplus==2 && nminus==2 );

  /* Config detection */
  for ( ia=0; ia<6; ++ia ) {
    i0 = MMG5_iare[ia][0];
    i1 = MMG5_iare[ia][1];

    if ( fabs(v[i0]) < MMG5_EPSD2 || fabs(v[i1]) < MMG5_EPSD2 )  continue;
    else if ( MG_SMSGN(v[i0],v[i1]) )  continue;

    MG_SET(flag,ia);

    /* Computation of the intersection between edges and isovalue */
    lam = v[i0] / (v[i0]-v[i1]);
    o[3*ia  ] = ppt[i0]->c[0] + lam*(ppt[i1]->c[0]-ppt[i0]->c[0]);
    o[3*ia+1] = ppt[i0]->c[1] + lam*(ppt[i1]->c[1]-ppt[i0]->c[1]);
    o[3*ia+2] = ppt[i0]->c[2] + lam*(ppt[i1]->c[2]-ppt[i0]->c[2]);
  }

  assert ( flag==30 || flag==45 || flag==51 );

  /* Set permutation of vertices : reference configuration 30 */
  tau[0] = 0 ; tau[1] = 1 ; tau[2] = 2 ; tau[3] = 3;
  taued = &MMG5_permedge[0][0];

  switch(flag){
  case 45:
    tau[0] = 1 ; tau[1] = 3 ; tau[2] = 2 ; tau[3] = 0;
    taued = &MMG5_permedge[5][0];
    break;

  case 51:
    tau[0] = 1 ; tau[1] = 2 ; tau[2] = 0 ; tau[3] = 3;
    taued = &MMG5_permedge[4][0];
    break;
  }

  if ( pm < 0 ) {
    if ( v[tau[0]] < 0.0 ) {
      /* compute the area that contains tau[0] and tau[1] */
      cfg = 0;
    }
    else {
      /* compute the area that contains tau[2] and tau[3] */
      cfg = 2;
    }
  }
  else {
    assert ( pm > 0 );
    if ( v[tau[0]] < 0.0 ) {
      /* compute the area that contains tau[0] and tau[1] */
      cfg = 2;
    }
    else {
      cfg = 0;
    }
  }
  assert ( cfg == 0 || cfg == 2 );

  /* Computation of the area depending on the detected config */
  if ( cfg == 0 ) {
    vfp  = fabs ( MMG5_det4pt(ppt[tau[0]]->c,ppt[tau[1]]->c,&o[3*taued[3]],&o[3*taued[4]]) );
    vfp += fabs ( MMG5_det4pt(ppt[tau[0]]->c,&o[3*taued[4]],&o[3*taued[3]],&o[3*taued[2]]) );
    vfp += fabs ( MMG5_det4pt(ppt[tau[0]]->c,&o[3*taued[3]],&o[3*taued[1]],&o[3*taued[2]]) );
#ifdef NDEBUG
    return vfp;
#endif
  }
  else if ( cfg == 2 ) {
    vfm  = fabs ( MMG5_det4pt(&o[3*taued[2]],&o[3*taued[4]],ppt[tau[2]]->c,ppt[tau[3]]->c) );
    vfm += fabs ( MMG5_det4pt(&o[3*taued[2]],&o[3*taued[3]],ppt[tau[2]]->c,&o[3*taued[4]]) );
    vfm += fabs ( MMG5_det4pt(&o[3*taued[1]],&o[3*taued[3]],ppt[tau[2]]->c,&o[3*taued[2]]) );
#ifdef NDEBUG
    return vfm;
#endif
  }

#ifndef NDEBUG
  /** Checks for debug mode */
  /* Compute the complementary area */
  if ( cfg == 0 ) {
    vfm  = fabs ( MMG5_det4pt(&o[3*taued[2]],&o[3*taued[4]],ppt[tau[2]]->c,ppt[tau[3]]->c) );
    vfm += fabs ( MMG5_det4pt(&o[3*taued[2]],&o[3*taued[3]],ppt[tau[2]]->c,&o[3*taued[4]]) );
    vfm += fabs ( MMG5_det4pt(&o[3*taued[1]],&o[3*taued[3]],ppt[tau[2]]->c,&o[3*taued[2]]) );
  }
  else if ( cfg == 2 ) {
    vfp  = fabs ( MMG5_det4pt(ppt[tau[0]]->c,ppt[tau[1]]->c,&o[3*taued[3]],&o[3*taued[4]]) );
    vfp += fabs ( MMG5_det4pt(ppt[tau[0]]->c,&o[3*taued[4]],&o[3*taued[3]],&o[3*taued[2]]) );
    vfp += fabs ( MMG5_det4pt(ppt[tau[0]]->c,&o[3*taued[3]],&o[3*taued[1]],&o[3*taued[2]]) );
  }

  /* vfp and vfm have been computed: check that the sum of both area is the
   * area of the whole tetra */
  double vf;
  vf = fabs ( MMG5_det4pt(ppt[0]->c,ppt[1]->c,ppt[2]->c,ppt[3]->c) );

  assert ( fabs(vf-(vfp+vfm)) < MMG5_EPSOK );

  if ( cfg==0 ) {
    return vfp;
  }
  else {
    return vfm;
  }
#endif

  /* Should not pass here */
  return 0.0;
}

/**
 * \param mesh pointer to the mesh.
 *
 * Reset mesh->info.isoref vertex and tetra references to 0.
 *
 * \warning to improve: for now, entities linked to the old ls (corners, required
 * points, normals/tangents, triangles and edges) are deleted in loadMesh. It
 * would be better to analyze wich entities must be kept and which ones must be
 * deleted depending on the split/nosplit info.
 */
int MMG3D_resetRef_ls(MMG5_pMesh mesh) {
  MMG5_pTetra     pt;
  MMG5_pPoint     p0;
  MMG5_int        k,ref;
  int8_t          i;

  /* Travel edges and reset tags at edges extremities */
  /* Reset ref and tags at ISO points */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<4; i++) {
      p0 = &mesh->point[pt->v[i]];
      /* Reset triangles */

      /* Reset vertices */
      if ( p0->ref == mesh->info.isoref ) {
        p0->ref = 0;
        /* Reset tags */
        p0->tag &= ~MG_CRN;
        p0->tag &= ~MG_REQ;
      }
    }
  }

  /* Reset the tetra references to their initial distribution */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];

    if ( !MG_EOK(pt) ) continue;

    /* If no material map is provided, reference is resetted to 0, otherwise,
     * reference has to exist in the material map because it is not possible to
     * decide for the user if a ref that is not listed has to be preserved or
     * resetted to 0 */
    if( !MMG5_getStartRef(mesh,pt->ref,&ref) ) return 0;
    pt->ref = ref;
  }

  return 1;
}


/**
 * \remark Not used.
 *
 * solve 3*3 non symmetric system Ar = b
 *
 */
static inline int
MMG5_invsl(double A[3][3],double b[3],double r[3]) {
  double detA;

  detA = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2]) \
    - A[0][1]*(A[1][0]*A[2][2] - A[2][0]*A[1][2]) \
    + A[0][2]*(A[1][0]*A[2][1] - A[2][0]*A[1][1]);
  if ( detA < MMG5_EPSD )  return 0;
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

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set values.
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

int MMG3D_ismaniball(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int k,int indp) {
  MMG5_pTetra   pt,pt1;
  double        v,v0,v1,v2;
  int           ibdy,ilist,cur,l;
  MMG5_int      *adja,list[MMG3D_LMAX+1],bdy[MMG3D_LMAX+1],jel,np,iel,res,base;
  int8_t        i,i0,i1,i2,j0,j1,j2,j,ip,nzeros,nopp,nsame;
  static int8_t mmgWarn0 = 0;

  pt = &mesh->tetra[k];
  np = pt->v[indp];
  if ( fabs(sol->m[np]) > MMG5_EPSD2 )  return 1;

  memset(bdy,0,(MMG3D_LMAX+1)*sizeof(MMG5_int));
  memset(list,0,(MMG3D_LMAX+1)*sizeof(MMG5_int));

  /* Sign of a starting point in ball of np */
  for (j=0; j<3; j++) {
    ip = MMG5_idir[indp][j];
    if ( sol->m[pt->v[ip]] != 0.0 )  break;
  }
  if ( j == 3 ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s:  at least 1 tetra with 4 null"
              " values.\n",__func__);
    }
    return 0;
  }

  v = sol->m[pt->v[ip]];
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
        i1 = MMG5_idir[i][j];
        v1 = sol->m[pt->v[i1]];
        if ( ( v1 != 0.0 ) && !MG_SMSGN(v,v1) ) {
          res = 4*iel + i;
          break;
        }
      }
    }

    /* Pile up faces sharing a vertex with same sign as v */
    for (j=0; j<3; j++) {
      i1 = MMG5_idir[i][MMG5_inxt2[j]];
      i2 = MMG5_idir[i][MMG5_iprv2[j]];
      v1 = sol->m[pt->v[i1]];
      v2 = sol->m[pt->v[i2]];

      if ( ( ( v1 != 0.0 ) && MG_SMSGN(v,v1) ) ||
           ( ( v2 != 0.0 ) && MG_SMSGN(v,v2) ) ) {
        jel = adja[MMG5_idir[i][j]];
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

  /* Fill in list bdy, corresponding to the support tetras of the boundary to be created */
  ibdy = 0;
  for(l=0; l<ilist; l++) {
    iel = list[l] / 4;
    i   = list[l] % 4;
    pt  = &mesh->tetra[iel];

    nzeros = nsame = nopp = 0;

    i0 = MMG5_idir[i][0];
    i1 = MMG5_idir[i][1];
    i2 = MMG5_idir[i][2];

    v0 = sol->m[pt->v[i0]];
    v1 = sol->m[pt->v[i1]];
    v2 = sol->m[pt->v[i2]];

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

    /* If no starting face with one vertex with opposite sign to v has been found,
     the only possibility for an admissible config is that adjacent to a face with 3 values equal to 0 has such vertex;
        v0,v1 are reused */
    if ( !res && nzeros == 2 && nsame == 1 ) {
      for (j=0; j<3; j++) {
        i0 = MMG5_idir[i][j];
        v0 = sol->m[pt->v[i0]];
        if ( v0 != 0.0 && MG_SMSGN(v,v0) ) break;
      }

      adja = &mesh->adja[4*(iel-1)+1];
      jel = adja[i0] / 4;
      j0 = adja[i0] % 4;
      pt1 = &mesh->tetra[jel];
      v1 = sol->m[pt1->v[j0]];
      if ( v1 != 0.0 && !MG_SMSGN(v,v1) ) {
        for (j=0; j<4; j++)
          if ( pt1->v[j] == np ) break;
        res = 4*jel+j;
      }
    }

    if ( ( nzeros == 2 && nsame == 1 ) || ( nsame >= 1 && nopp >= 1 ) )  {
      bdy[ibdy] = list[l];
      ibdy++;
    }
  }

  /* Invalid configuration has been created */
  if ( !res )
    return 0;

  /* Reset the current part of the ball, and start back the process with the other sign */
  iel = res / 4;
  pt = &mesh->tetra[iel];
  base = ++mesh->base;
  pt->flag = base;

  memset(list,0,(MMG3D_LMAX+1)*sizeof(MMG5_int));
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
      i1 = MMG5_idir[i][MMG5_inxt2[j]];
      i2 = MMG5_idir[i][MMG5_iprv2[j]];
      v1 = sol->m[pt->v[i1]];
      v2 = sol->m[pt->v[i2]];

      if ( v1 == 0.0 && v2 == 0.0 ) {
        jel = adja[MMG5_idir[i][j]];
        if( !jel ) continue;
        jel /=4 ;
        pt1 = &mesh->tetra[jel];
        pt1->flag = base;
      }

      else if ( ( ( v1 != 0.0 ) && (!MG_SMSGN(v,v1)) ) || ( ( v2 != 0.0 ) && (!MG_SMSGN(v,v2)) ) ) {
        jel = adja[MMG5_idir[i][j]];
        if( !jel ) continue;
        jel /= 4;
        pt1 = &mesh->tetra[jel];

        j0 = MMG5_idir[i][0];
        j1 = MMG5_idir[i][1];
        j2 = MMG5_idir[i][2];

        v0 = sol->m[pt1->v[j0]];
        v1 = sol->m[pt1->v[j1]];
        v2 = sol->m[pt1->v[j2]];

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
          if ( pt1->flag < base - 1 ) return 0;
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
    if ( pt->flag != base ) return 0;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set function.
 * \return 1 if success, 0 if fail.
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
int MMG3D_snpval_ls(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTetra   pt;
  MMG5_pPoint   p0;
  double        *tmp;
  MMG5_int      k,nc,ns,ip,ncg;
  int8_t        i;

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Error: %s: hashing problem (1). Exit program.\n",
      __func__);
    return 0;
  }

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double,return 0);

  /* Include tetras with very poor quality that are connected to the negative part */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] ) continue;
    if ( pt->qual < MMG5_EPS ) {
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        if ( sol->m[ip] < 1000.0*MMG5_EPS ) break;
      }
      if ( i < 4 ) {
        for (i=0; i<4; i++) {
          ip = pt->v[i];
          sol->m[ip] = -1000.0*MMG5_EPS;
        }
      }
    }
  }

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]) < MMG5_EPS ) {
      if ( mesh->info.ddebug )
        fprintf(stderr,"  ## Warning: %s: snapping value at vertex %" MMG5_PRId "; "
                "previous value: %E.\n",__func__,k,fabs(sol->m[k]));

      tmp[k] = ( fabs(sol->m[k]) < MMG5_EPSD ) ?
        (-100.0*MMG5_EPS) : sol->m[k];
      p0->flag = 1;
      sol->m[k] = 0;
      ns++;
    }
  }

  ncg = 0;
  do {
    nc = 0;
    /* Check snapping did not lead to a nonmanifold situation */
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        p0 = &mesh->point[ip];
        if ( p0->flag == 1 ) {
          if ( !MMG3D_ismaniball(mesh,sol,k,i) ) {
            if ( tmp[ip] < 0.0 )
              sol->m[ip] = -100.0*MMG5_EPS;
            else
              sol->m[ip] = +100.0*MMG5_EPS;

            p0->flag = 0;
            nc++;
          }
        }
      }
    }
    ncg += nc;
  }
  while ( nc );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+ncg > 0 )
    fprintf(stdout,"     %8" MMG5_PRId " points snapped, %" MMG5_PRId " corrected\n",ns,ncg);

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;
  
  /* memory free */
  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,tmp);

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the level-set
 *
 * \return 1 if success, 0 otherwise
 *
 * Removal of small parasitic components (bubbles of material, etc) with volume
 * less than mesh->info.rmc (default VOLFRAC) * volume of the mesh.
 *
 */
int MMG3D_rmc(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTetra    pt,pt1,pt2;
  MMG5_pxTetra   pxt;
  double         volc,voltot,v0,v1,v2,v3;
  int            l,cur,ipile;
  MMG5_int       ncp,ncm,base,k,kk,ll,ip0,ip1,ip2,ip3,*adja,*pile;
  int8_t         i,j,i1,onbr;

  ncp = 0;
  ncm = 0;

  /* Erase tetra flags */
  for (k=1; k<=mesh->ne; k++) mesh->tetra[k].flag = 0;

  /* Calculate volume of the total mesh (x6 to avoid useless division)*/
  voltot = 0.0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    voltot += fabs ( MMG5_orvol(mesh->point,pt->v) );
  }

  /* Memory allocation for pile */
  MMG5_ADD_MEM(mesh,(mesh->ne+1)*sizeof(MMG5_int),"temporary table",
               printf("  Exit program.\n");
               return 0);
  MMG5_SAFE_CALLOC(pile,mesh->ne+1,MMG5_int,return 0);

  /* Investigate only positive connected components */
  base = ++mesh->base;
  
  for (k=1; k<=mesh->ne; k++) {
    ipile = 0;
    volc  = 0.0;
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->flag == base ) continue;

    /* Checks signs of the LS function at the 4 vertices of pt */
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];
    ip3 = pt->v[3];

    v0 = sol->m[ip0];
    v1 = sol->m[ip1];
    v2 = sol->m[ip2];
    v3 = sol->m[ip3];

    if ( v0 <= 0.0 && v1 <= 0.0 && v2 <= 0.0 && v3 <= 0.0 ) continue;

    /* Add tetra to pile if one vertex is > 0 */
    pt->flag = base;
    pile[ipile] = k;
    ipile++;
    if ( ipile > mesh->ne ) {
      fprintf(stderr,"\n  ## Problem in length of pile; function rmc.\n"
              " Check that the level-set intersect the mesh.\n"
              " Exit program.\n");

      return 0;
    }

    /* Pile up all the positive connected component attached to the first tetra */
    cur = 0;
    do {
      kk = pile[cur];
      pt1 = &mesh->tetra[kk];

      /* Add local volume fraction of the positive subdomain to volc */
      volc += MMG3D_vfrac(mesh,sol,kk,1);

      /* Add adjacent tetra to kk via positive vertices to the pile, if need be */
      adja = &mesh->adja[4*(kk-1)+1];
      for (i=0; i<4; i++) {
        ip0 = pt1->v[i];
        if ( sol->m[ip0] <= 0.0 ) continue;

        for ( i1=0; i1<3; ++i1 ) {
          ll = adja[MMG5_idir[i][i1]] / 4;
          if ( !ll ) continue;

          pt2 = &mesh->tetra[ll];
          if ( MG_EOK(pt2) && pt2->flag != base ) {
            pt2->flag   = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->ne ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }
      }
    }
    while ( ++cur < ipile );

    /* Remove connected component if its volume is too small */
    if ( volc < mesh->info.rmc * voltot ) {
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tetra[pile[l]];
        for (i=0; i<4; i++) {
          ip0 = pt1->v[i];
          if ( sol->m[ip0] > 0.0 ) {
            sol->m[ip0] = - 100*MMG5_EPS;
          }
        }
      }
      ncp++;
    }
  }
  
  /* Investigate only negative connected components */
  base = ++mesh->base;

  for (k=1; k<=mesh->ne; k++) {
    ipile = 0;
    volc  = 0.0;
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->flag == base ) continue;

    /* Checks signs of the LS function at the 4 vertices of pt */
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];
    ip3 = pt->v[3];

    v0 = sol->m[ip0];
    v1 = sol->m[ip1];
    v2 = sol->m[ip2];
    v3 = sol->m[ip3];

    if ( v0 >= 0.0 && v1 >= 0.0 && v2 >= 0.0 && v3 >= 0.0 ) continue;

    /* Add tetra to pile if one vertex is > 0 */
    pt->flag = base;
    pile[ipile] = k;
    ipile++;
    if ( ipile > mesh->ne ) {
      fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
      return 0;
    }

    /* Pile up all the negative connected component attached to the first tetra */
    cur = 0;
    do {
      kk = pile[cur];
      pt1 = &mesh->tetra[kk];

      /* Add local volume fraction of the negative subdomain to volc */
      volc += MMG3D_vfrac(mesh,sol,kk,-1);

      /* Add adjacent tetra to kk via negative vertices to the pile, if need be */
      adja = &mesh->adja[4*(kk-1)+1];
      for (i=0; i<4; i++) {
        ip0 = pt1->v[i];
        if ( sol->m[ip0] >= 0.0 ) continue;

        for ( i1=0; i1<3; ++i1 ) {
          ll = adja[MMG5_idir[i][i1]] / 4;
          if ( !ll ) continue;

          pt2 = &mesh->tetra[ll];
          if ( MG_EOK(pt2) && pt2->flag != base ) {
            pt2->flag   = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->ne ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }
      }
    }
    while ( ++cur < ipile );

    /* Remove connected component if its volume is too small */
    if ( volc < mesh->info.rmc * voltot ) {
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tetra[pile[l]];
        for (i=0; i<4; i++) {
          ip0 = pt1->v[i];
          if ( sol->m[ip0] < 0.0 ) sol->m[ip0] = 100*MMG5_EPS;
        }
      }
      ncm++;
    }
    
    /* Remove connected component if it is not attached to one base reference */
    if ( mesh->info.nbr ) {
      onbr = 0;
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tetra[pile[l]];
        if ( pt1->xt ) {
          pxt = &mesh->xtetra[pt1->xt];
          for (i=0; i<4; i++) {
            if ( MMG5_isbr(mesh,pxt->ref[i]) ) {
              for (j=0; j<3; j++) {
                ip0 = pt1->v[MMG5_idir[i][j]];
                if ( sol->m[ip0] < 0.0 )  {
                  onbr = 1;
                  break;
                }
              }
            }
          }
        }
        if ( onbr ) break;
      }

      if ( !onbr ) {
        for (l=0; l<ipile; l++) {
          pt1 = &mesh->tetra[pile[l]];
          for (i=0; i<4; i++) {
            ip0 = pt1->v[i];
            if ( sol->m[ip0] < 0.0 ) sol->m[ip0] = 100*MMG5_EPS;
          }
        }
        ncm++;
      }
    }
  }

  /* Erase tetra flags */
  for (k=1; k<=mesh->ne; k++) mesh->tetra[k].flag = 0;

  /* Release memory */
  MMG5_DEL_MEM(mesh,pile);

  if ( mesh->info.imprim > 0 || mesh->info.ddebug ) {
    printf("\n  *** Removed %" MMG5_PRId " positive parasitic bubbles and %" MMG5_PRId " negative parasitic bubbles\n",ncp,ncm);
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set values.
 * \param met pointer to a metric (non-mandatory).
 * \return 1 if success, 0 otherwise.
 *
 * Proceed to discretization of the implicit function carried by sol into mesh,
 * once values of sol have been snapped/checked
 *
 */
int MMG3D_cuttet_ls(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met){
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  MMG5_Hash     hash;
  double        c[3],v0,v1,s;
  int           ier;
  MMG5_int      vx[6],k,ip0,ip1,np,nb,ns,ne,src,refext,refint;
  int8_t        ia,j,npneg;
  static int8_t mmgWarn = 0;

  /* reset point flags and h */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      v0  = sol->m[ip0];
      v1  = sol->m[ip1];
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        if ( !p0->flag ) {
          p0->flag = ++nb;
        }
        if ( !p1->flag ) {
          p1->flag = ++nb;
        }
      }
    }
  }
  if ( ! nb )  return 1;

  /* Create intersection points at 0 isovalue and set flags to tetras */
  if ( !MMG5_hashNew(mesh,&hash,nb,7*nb) ) return 0;
  /* Hash all boundary and required edges, and put ip = -1 in hash structure */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (ia=0; ia<6; ia++) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
      continue;
    }

    if ( !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for (ia=0; ia<4; ia++) {
      if ( !(pxt->ftag[ia] & MG_BDY) ) continue;

      for (j=0; j<3; j++) {
        if ( !(pxt->tag[ MMG5_iarf[ia][j] ] & MG_REQ) ) continue;

        ip0 = pt->v[MMG5_idir[ia][MMG5_inxt2[j]]];
        ip1 = pt->v[MMG5_idir[ia][MMG5_iprv2[j]]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
    }
  }


  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      np  = MMG5_hashGet(&hash,ip0,ip1);

      if ( np>0 )  continue;

      if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0];
      v1 = sol->m[ip1];
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      npneg = (np<0);

      s = v0 / (v0-v1);

      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

#ifdef USE_POINTMAP
      src = p0->src;
#else
      src = 1;
#endif
      np = MMG3D_newPt(mesh,c,0,src);
      if ( !np ) {
        MMG5_int oldnpmax = mesh->npmax;
        MMG3D_POINT_REALLOC(mesh,sol,np,MMG5_GAP,
                             fprintf(stderr,"\n  ## Error: %s: unable to"
                                     " allocate a new point\n",__func__);
                             MMG5_INCREASE_MEM_MESSAGE();
                             return 0
                             ,c,0,src);
        if( met ) {
          if( met->m ) {
            MMG5_ADD_MEM(mesh,(met->size*(mesh->npmax-met->npmax))*sizeof(double),
                         "larger solution",
                         MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                         mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                         mesh->npmax = oldnpmax;
                         mesh->np = mesh->npmax-1;
                         mesh->npnil = 0;
                         return 0);
            MMG5_SAFE_REALLOC(met->m,met->size*(met->npmax+1),
                              met->size*(mesh->npmax+1),
                              double,"larger solution",
                              MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                              mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                              mesh->npmax = oldnpmax;
                              mesh->np = mesh->npmax-1;
                              mesh->npnil = 0;
                              return 0);
          }
          met->npmax = mesh->npmax;
        }
      }
      sol->m[np] = 0;
      /* If user provide a metric, interpolate it at the new point */
      if ( met && met->m ) {
        if ( met->size > 1 ) {
          ier = MMG3D_intmet33_ani(mesh,met,k,ia,np,s);
        }
        else {
          ier = MMG5_intmet_iso(mesh,met,k,ia,np,s);
        }
        if ( ier <= 0 ) {
          // Unable to compute the metric
          fprintf(stderr,"\n  ## Error: %s: unable to"
                  " interpolate the metric during the level-set"
                  " discretization\n",__func__);
          return 0;
        }
      }

      if ( npneg ) {
        /* We split a required edge */
        if ( !mmgWarn ) {
          mmgWarn = 1;
          fprintf(stderr,"  ## Warning: %s: the level-set intersect at least"
                  " one required entity. Required entity ignored.\n\n",__func__);
        }
        MMG5_hashUpdate(&hash,ip0,ip1,np);
      }
      else
        MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }

  /* Proceed to splitting, according to flags to tets */
  ne  = mesh->ne;
  ns  = 0;
  ier = 1;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = 0;
    memset(vx,0,6*sizeof(MMG5_int));
    for (ia=0; ia<6; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
      if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
    }
    switch (pt->flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
      ier = MMG5_split1(mesh,met,k,vx,1);
      ns++;
      break;

    case 48: case 24: case 40: case 6: case 34: case 36:
    case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
      ier = MMG5_split2sf(mesh,met,k,vx,1);
      ns++;
      break;

    case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
      ier = MMG5_split3cone(mesh,met,k,vx,1);
      ns++;
      break;

    case 30: case 45: case 51:
      ier = MMG5_split4op(mesh,met,k,vx,1);
      ns++;
      break;

    default :
      assert(pt->flag == 0);
      break;
    }
    if ( !ier ) return 0;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;
}


/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set values.
 * \return 1.
 *
 * Set references to tets according to the sign of the level set function.
 *
 */
int MMG3D_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTetra   pt;
  double        v;
  int           ier;
  MMG5_int      ref,refint,refext,k,ip;
  int8_t        nmns,npls,nz,i;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    ref = pt->ref;

    nmns = npls = nz = 0;
    for (i=0; i<4; i++) {
      ip = pt->v[i];
      v  = sol->m[ip];
      if ( v > 0.0 )
        npls++;
      else if ( v < 0.0 )
        nmns++;
      else
        nz ++;
    }
    /* Remark: this test is not consistent with the test of the lssurf option
     * because it autorizes the level-set to be superposed with the surface */
    assert(nz < 4);
    ier = MMG5_isSplit(mesh,ref,&refint,&refext);

    if ( npls ) {
      if ( ier ) {
        assert(!nmns);
        pt->ref = refext;
      }
    }
    else {
      if ( ier ) {
        assert(nmns);
        pt->ref = refint;
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \return 1 if success, 0 if the xtetra array can't be reallocated.
 *
 * Update the xtetra array to store the new bdy faces created by the isosurface
 * discretization.
 *
 */
int MMG3D_update_xtetra ( MMG5_pMesh mesh ) {
  MMG5_pTetra   pt,pt1,ptmax,ptmin;
  MMG5_pxTetra  pxt;
  int           i,j,imax,imin;
  MMG5_int      *adja,k,jel;

  if ( (!mesh->info.iso) || (!mesh->info.opnbdy) ) {
    /* In non opnbdy mode, info stored in xtetra is not used */
    /* In non ls mode, xtetra are alread updated */
    return 1;
  }

  /* Opnbdy mode uses data stored in xtetra but in iso mode, the new triangles
   * created by the ls discretization haven't been stored inside the xtetra */
  if ( !mesh->xtetra ) {
    fprintf(stderr,"\n  ## Error: %s: the xtetra array must be allocated.\n",
      __func__);
    return 0;
  }
  if ( !mesh->adja ) {
    fprintf(stderr,"\n  ## Error: %s: the ajda array must be allocated.\n",
      __func__);
    return 0;
  }


  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    adja = &mesh->adja[4*(k-1)+1];

    for (i=0; i<4; i++) {
      if ( !adja[i] ) {
        /* Face is already stored */
        continue;
      }

      jel = adja[i]/4;
      pt1 = &mesh->tetra[jel];

      if ( pt->ref == pt1->ref ) {
        /* Potential opnbdy face is already stored */
        continue;
      }

      j = adja[i]%4;
      /* Detection of the tetra of higher ref */
      if ( pt->ref > pt1->ref ) {
        ptmax = pt;
        imax  = i;
        ptmin = pt1;
        imin  = j;
      }
      else {
        ptmax = pt1;
        imax  = j;
        ptmin = pt;
        imin  = i;
      }

      /* Update the xtetra array for both tetra */
      /* Tetra ptmax */
      if ( !ptmax->xt ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                            "larger xtetra table",
                            mesh->xt--;
                            fprintf(stderr,"  Exit program.\n");return 0;);
        }
        ptmax->xt = mesh->xt;
      }

      pxt = &mesh->xtetra[ptmax->xt];
      pxt->ref[imax]   = mesh->info.isoref;
      pxt->ftag[imax] |= MG_BDY;
      MG_SET(pxt->ori,imax);

      /* Tetra ptmin */
      if ( !ptmin->xt ) {
        mesh->xt++;
        if ( mesh->xt > mesh->xtmax ) {
          MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,MMG5_GAP,MMG5_xTetra,
                            "larger xtetra table",
                            mesh->xt--;
                            fprintf(stderr,"  Exit program.\n");return 0;);
        }
        ptmin->xt = mesh->xt;
      }

      pxt = &mesh->xtetra[ptmin->xt];
      pxt->ref[imin]   = mesh->info.isoref;
      pxt->ftag[imin] |= MG_BDY;
      MG_CLR(pxt->ori,imin);
    }
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param start index of the starting tetra
 * \param ip point index
 *
 * \return 1 if success, 0 if fail
 *
 * Check whether implicit surface is orientable in ball of point ip in tet iel ;
 * Beware : may return 0 when implicit boundary is tangent to outer boundary
 *
 */
int MMG3D_chkmaniball(MMG5_pMesh mesh, MMG5_int start, int8_t ip){
  MMG5_pTetra    pt,pt1;
  int            ilist,cur,nref;
  MMG5_int       base,ref,*adja,list[MMG3D_LMAX+2],k,k1,nump;
  int8_t         i,l,j;

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

    adja = &mesh->adja[4*(k-1)+1];
    for(l=0; l<3; l++){
      i = MMG5_inxt3[i];

      /* Travel only through non boundary faces. */
      k1 = adja[i];
      if(!k1) continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      if( MMG5_isNotSplit(mesh,pt1->ref) ) continue;

      if( pt1->ref != ref ) continue;

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

    adja = &mesh->adja[4*(k-1)+1];
    for(l=0; l<3; l++){
      i = MMG5_inxt3[i];

      k1 = adja[i];
      if ( !k1 ) continue;
      k1/=4;

      pt1 = &mesh->tetra[k1];
      if( MMG5_isNotSplit(mesh,pt1->ref) ) continue;

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
      fprintf(stderr,"   *** Topological problem\n");
      fprintf(stderr,"       non manifold surface at point %" MMG5_PRId " %" MMG5_PRId "\n",nump, MMG3D_indPt(mesh,nump));
      fprintf(stderr,"       non manifold surface at tet %" MMG5_PRId " (ip %d)\n", MMG3D_indElt(mesh,start),ip);
      fprintf(stderr,"       nref (color %d) %" MMG5_PRId "\n",nref,ref);
      return 0;
    }
  }

  return 1;
}

/** Check whether implicit surface enclosed in volume is orientable */
int MMG3D_chkmani(MMG5_pMesh mesh){
  MMG5_pTetra   pt,pt1;
  MMG5_int      ref;
  MMG5_int      iel,k,*adja;
  int8_t        i,j,ip,cnt;
  static int8_t mmgWarn0 = 0;

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
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Warning: %s: at least 1 tetra with 4 boundary"
                " faces.\n",__func__);
      }
      //return 0;
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
      if( !MMG5_isLevelSet(mesh,pt1->ref,pt->ref) ) continue;

      for(j=0; j<3; j++){
        ip = MMG5_idir[i][j];

        /* If the starting point is MG_PARBDY: this is not a non-manifold topology */
        /*    - True  for centralized input in parmmg */
        /*    - Wrong for distributed input in parmmg: TODO */
        if ( !(mesh->point[pt->v[ip]].tag & MG_PARBDY)) {
          if(!MMG3D_chkmaniball(mesh,k,ip))
            return 0;
        }
      }
    }
  }

  if ( mesh->info.imprim > 0 || mesh->info.ddebug )
    fprintf(stdout,"  *** Manifold implicit surface.\n");
  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the metric
 *
 * \return 1 if success, 0 otherwise.
 *
 * Check whether implicit surface enclosed in volume is orientable (perform an
 * additionnal test w.r.t. MMG3D_chkmani)
 *
 */
int MMG3D_chkmani2(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTetra    pt,pt1;
  MMG5_int       k,iel;
  MMG5_int       *adja;
  int8_t         i,j,ip,cnt;

  for(k=1; k<=mesh->np; k++){
    mesh->point[k].flag = 0;
  }

  /** First test : assure no tetra has its 4 vertices on implicit boundary */
  for(k=1; k<=mesh->ne; k++){
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || (pt->tag & MG_REQ))   continue;

    cnt = 0;
    for(j=0; j<4; j++) {
      if( sol->m[pt->v[j]] == 0.0 ) cnt++;
    }
    if(cnt == 4) {
      fprintf(stderr,"\n  ## Error: %s: tetra %" MMG5_PRId ": 4 vertices on implicit boundary.\n",
              __func__,k);
      return 0;
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
        ip = MMG5_idir[i][j];

        /* If the starting point is MG_PARBDY: this is not a non-manifold topology */
        /*    - True  for centralized input in parmmg */
        /*    - Wrong for distributed input in parmmg: TODO */
        if ( !(mesh->point[pt->v[ip]].tag & MG_PARBDY)) {
          if(!MMG3D_chkmaniball(mesh,k,ip)){
            fprintf(stderr,"\n  ## Error: %s: non orientable implicit surface:"
                    " ball of point %" MMG5_PRId ".\n",__func__,pt->v[ip]);
            return 0;
          }
        }
      }
    }
  }
  if ( mesh->info.ddebug )  fprintf(stdout,"  *** Manifold implicit surface.\n");
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param k index of element in which we collapse.
 * \param iface face through wich we perform the collapse
 * \param iedg edge to collapse
 * \param ndepmin index of an elt with ref refmin and outside the shell of edge.
 * \param ndepplus ndex of an elt with ref refplus and outside the shell of edge.
 * \param refmin reference of one of the two subdomains in presence
 * \param refplus reference of the other subdomain in presence
 * \param isminp 1 if we have found a tetra with ref refmin
 * \param isplp 1 if we have found a tetra with ref refplus
 * \return 0 if we create a non manifold situation, 1 otherwise
 *
 * Check whether collapse of point np to nq does not create a non manifold
 * situation at nq ndepmin, ndepplus = tetra of ref minus, plus in ball of np,
 * not in shell of (np,nq).
 *
 */
int MMG3D_chkmanicoll(MMG5_pMesh mesh,MMG5_int k,int iface,int iedg,MMG5_int ndepmin,MMG5_int ndepplus,MMG5_int refmin,MMG5_int refplus,int8_t isminp,int8_t isplp) {
  MMG5_pTetra    pt,pt1;
  int            ilist,cur;
  MMG5_int       stor;
  MMG5_int       ref,nump,numq,list[MMG3D_LMAX+2],*adja,*adja1,iel,jel,ndepmq,ndeppq,base;
  int8_t         i,j,ip,jp,iq,jq,voy,indp,indq,isminq,isplq,ismin,ispl;

  ilist = 0;
  ndepmq = ndeppq = 0;
  isplq = isminq = 0;

  pt    = &mesh->tetra[k];
  ip    = MMG5_idir[iface][MMG5_inxt2[iedg]];
  iq    = MMG5_idir[iface][MMG5_iprv2[iedg]];
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

    if ( pt->ref == refmin ) isminq = 1;
    else if ( pt->ref == refplus ) isplq = 1;

    cur = 0;
    while( cur < ilist ) {
      iel = list[cur] / 4;
      i = list[cur] % 4;
      adja = &mesh->adja[4*(iel-1)+1];

      for (j=0; j<3; j++) {
        i = MMG5_inxt3[i];
        jel = adja[i];
        if ( !jel ) continue;

        jel /= 4;
        pt1 = &mesh->tetra[jel];

        if ( pt1->ref == refmin ) isminq = 1;
        else if ( pt1->ref == refplus ) isplq = 1;

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
        if ( !ndeppq && pt1->ref == refplus ) {
          for(ip=0; ip<4; ip++)
            if ( pt1->v[ip] == nump ) break;
          if( ip == 4 ) ndeppq = jel;
        }
        if( !ndepmq && pt1->ref == refmin ) {
          for(ip=0; ip<4; ip++)
            if ( pt1->v[ip] == nump ) break;
          if( ip == 4 ) ndepmq = jel;
        }
      }
      cur++;
    }

    memset(list,0,(MMG3D_LMAX+2)*sizeof(MMG5_int));
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
      return 0;
    else
      return 1;
  }


  cur = 0;
  while ( cur < ilist ) {
    stor = list[cur];
    /* Element belongs to the ball of np */
    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;

      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;
      for (i=0; i<3; i++) {
        jp = MMG5_inxt3[jp];
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

      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;
      for (i=0; i<3; i++) {
        jq = MMG5_inxt3[jq];
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
      return 0;
    else
      return 1;
  }

  while ( cur < ilist ) {
    stor = list[cur];
    /* Element belongs to the ball of np */
    if ( stor <= 0 ) {
      stor *= -1;
      iel = stor / 4;
      ip  = stor % 4;

      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;

      for (i=0; i<3; i++) {
        jp = MMG5_inxt3[jp];
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

      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;

      for (i=0; i<3; i++) {
        jq = MMG5_inxt3[jq];
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
      adja = &mesh->adja[4*(iel-1)+1];

      jp = ip;
      for(i=0; i<3; i++) {
        jp = MMG5_inxt3[jp];
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

        /* Only tets of the shell of (np,nq) can be added, unless future ball is
         * non manifold */
        if ( indq == -1 ) {
          if ( mesh->info.ddebug ) {
          fprintf(stderr,"\n  ## Warning: %s: we should rarely passed here. "
                  "tetra %" MMG5_PRId " =  %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId ", ref = %" MMG5_PRId ".",__func__,
                  jel,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3],pt1->ref);
          }
          return 0;
        }

        list[ilist] = -(4*jel+indp);
        ilist++;
        assert( ilist < MMG3D_LMAX +1 );
      }
    }
    else {
      iel = stor / 4;
      iq  = stor % 4;
      adja = &mesh->adja[4*(iel-1)+1];

      jq = iq;
      for(i=0; i<3; i++) {
        jq = MMG5_inxt3[jq];
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
          if ( mesh->info.ddebug ) {
          fprintf(stderr,"\n  ## Warning: %s: we should rarely passed here. "
                  "tetra %" MMG5_PRId " =  %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId ", ref = %" MMG5_PRId "\n",__func__,
                  jel,pt1->v[0],pt1->v[1],pt1->v[2],pt1->v[3],pt1->ref);
          }
          return 0;
        }

        list[ilist] = 4*jel+indq;
        ilist++;
        assert( ilist < MMG3D_LMAX +1 );
      }
    }
    cur++;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set.
 * \param met pointer to  a metric (optionnal).
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int MMG3D_mmg3d2(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  char str[16]="";

  /* Set function pointers */
  if ( mesh->info.isosurf ) {
    strcat(str,"(BOUNDARY PART)");

    MMG3D_snpval   = MMG3D_snpval_lssurf;
    MMG3D_resetRef = MMG3D_resetRef_lssurf;
    MMG3D_cuttet   = MMG3D_cuttet_lssurf;
    MMG3D_setref   = MMG3D_setref_lssurf;
  }
  else {
    MMG3D_snpval   = MMG3D_snpval_ls;
    MMG3D_resetRef = MMG3D_resetRef_ls;
    MMG3D_cuttet   = MMG3D_cuttet_ls;
    MMG3D_setref   = MMG3D_setref_ls;
  }

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION %s\n",str);

  if ( mesh->nprism || mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }

  /* Work only with the 0 level set */
  MMG5_int k;
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /* Snap values of level set function if need be */
  if ( !MMG3D_snpval(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    return 0;
  }

  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Compatibility triangle orientation w/r tetras */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }

  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  /* Build hash table for initial edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    return 0;
  }

  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* Reset the mesh->info.isoref field everywhere it appears */
  if ( !MMG3D_resetRef(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    return 0;
  }

  /* Removal of small parasitic components */
  if ( mesh->info.iso ) {
    if ( mesh->info.rmc > 0. && !MMG3D_rmc(mesh,sol) ) {
      fprintf(stderr,"\n  ## Error in removing small parasitic components."
              " Exit program.\n");
      return 0;
    }
  }
  else {
    /* RMC : on verra */
    if ( mesh->info.rmc > 0 ) {
      fprintf(stdout,"\n  ## Warning: rmc option not implemented for boundary"
              " isosurface extraction.\n");
    }
  }

#ifdef USE_POINTMAP
  /* Initialize source point with input index */
  MMG5_int ip;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].src = ip;
#endif

  if ( !MMG3D_cuttet(mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
    return 0;
  }

  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);

  mesh->nt = 0;

  if ( !MMG3D_setref(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Clean old bdy analysis */
  for ( MMG5_int k=1; k<=mesh->np; ++k ) {
    if ( mesh->point[k].tag & MG_BDY ) {
      mesh->point[k].tag &= ~MG_BDY;
    }
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);

  return 1;
}
