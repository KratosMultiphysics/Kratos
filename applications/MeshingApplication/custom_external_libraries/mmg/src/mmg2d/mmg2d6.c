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
 * \file mmg2d/mmg2d6.c
 * \brief Isosurface discretization.
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
 * \param sol pointer toward the sol structure.
 * \return 1 if success.
 *
 * Isosurface discretization
 *
 **/

/* Identify whether a triangle with reference ref should be split, and the labels of the resulting triangles */
int MMG2D_isSplit(MMG5_pMesh mesh,int ref,int *refint,int *refext) {
  MMG5_pMat    pm;
  int          k;

  /* Check in the info->mat table if reference ref is supplied by the user */
  for (k=0; k<mesh->info.nmat; k++) {
    pm = &mesh->info.mat[k];
    if ( pm->ref == ref ) {
      if ( !pm->dospl ) return 0;
      else {
        *refint = pm->rin;
        *refext = pm->rex;
        return 1;
      }
    }
  }

  /* Default case: split with references MG_MINUS, MG_PLUS */
  *refint = MG_MINUS;
  *refext = MG_PLUS;
  return 1;

}

/* Retrieve the initial domain reference associated to the (split) reference ref */
int MMG2D_getIniRef(MMG5_pMesh mesh,int ref) {
  MMG5_pMat     pm;
  int           k;

  for (k=0; k<mesh->info.nmat; k++) {
    pm = &mesh->info.mat[k];
    if ( pm->ref == ref && !pm->dospl ) return pm->ref;
    if ( ref == pm->rin || ref == pm->rex ) return pm->ref;
  }
  return ref;
}

/* Reset MG_ISO vertex and edge references to 0 */
int MMG2D_resetRef(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0;
  int             k,ref;
  char            i;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      p0 = &mesh->point[pt->v[i]];
      if ( pt->edg[i] == MG_ISO ) pt->edg[i] = 0;
      if ( p0->ref == MG_ISO ) p0->ref = 0;
    }
  }

  /* Reset the triangle references to their initial distribution */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;
    ref = MMG2D_getIniRef(mesh,pt->ref);
    pt->ref = ref;
  }

  return 1;
}

/* Check whether snapping the value of vertex i of k to 0 exactly leads to a non manifold situation
 assumption: the triangle k has vertex i with value 0 and the other two with changing values */
int MMG2D_ismaniball(MMG5_pMesh mesh, MMG5_pSol sol, int start, char istart) {
  MMG5_pTria       pt;
  double           v1, v2;
  int              *adja,k,ip1,ip2,end1;
  char             i,i1,smsgn;
  static char      mmgWarn=0;

  k = start;
  i = MMG5_inxt2[istart];

  /* First loop: stop if an external boundary, or a change in signs (or a 0) is met
     recall that MG_SMGSGN(a,b) = 1 provided a*b >0 */
  do{
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i] / 3;
    i1 = adja[i] % 3;
    i = MMG5_iprv2[i1];

    if ( k==0 ) break;

    pt = &mesh->tria[k];
    ip1 = pt->v[i1];
    ip2 = pt->v[i];

    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    smsgn = MG_SMSGN(v1,v2) ? 1 : 0;
  }
  while ( smsgn );

  end1 = k;
  k = start;
  i = MMG5_iprv2[istart];

  /* Second loop: same travel in the opposite sense */
  do{
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i] / 3;
    i1 = adja[i] % 3;
    i = MMG5_inxt2[i1];

    if ( k==0 ) break;

    pt = &mesh->tria[k];
    ip1 = pt->v[i1];
    ip2 = pt->v[i];

    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    smsgn = MG_SMSGN(v1,v2) ? 1 : 0;
  }
  while ( smsgn );

  /* If first stop was due to an external boundary, the second one must too;
     else, the final triangle for the first travel must be that of the second one */
  if ( k != end1 ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: unsnap at least 1 point "
              "(point %d in tri %d).\n",__func__,MMG2D_indElt(mesh,start),
              MMG2D_indPt(mesh,mesh->tria[start].v[istart]));
    }
    return 0;
  }
  return 1;
}

/* Snap values of sol very close to 0 to 0 exactly (to avoid very small triangles in cutting) */
int MMG2D_snapval(MMG5_pMesh mesh, MMG5_pSol sol, double *tmp) {
  MMG5_pTria       pt,pt1;
  MMG5_pPoint      p0;
  double           v1,v2;
  int              k,kk,iel,ns,nc,ip,ip1,ip2,npl,nmn,ilist,list[MMG2D_LONMAX+2];
  char             i,j,j1,j2;

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]) < MMG5_EPS ) {
      tmp[k] =  - 100.0*MMG5_EPS;
      p0->flag = 1;
      sol->m[k] = 0.0;
      ns++;
    }
  }

  /* Check that the snapping process has not led to a nonmanifold situation */
  /* TO DO: Check first that every snapped point corresponds to a change in signs */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ip1 = pt->v[MMG5_inxt2[i]];
      ip2 = pt->v[MMG5_iprv2[i]];

      p0 = &mesh->point[ip];
      v1 = sol->m[ip1];
      v2 = sol->m[ip2];

      /* Catch a snapped point by a triangle where there is a sign change */
      if ( p0->flag && !(MG_SMSGN(v1,v2)) ) {
        if ( !MMG2D_ismaniball(mesh,sol,k,i) ) {
          sol->m[ip] = tmp[ip];
          nc++;
        }
        p0->flag = 0;
      }
    }
  }

  /* Check that the ls function does not show isolated spots with 0 values (without sign changes) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      if ( fabs(sol->m[ip]) >= MMG5_EPS ) continue;
      npl = nmn = 0;
      ilist = MMG2D_boulet(mesh,k,i,list);
      for(kk=0; kk<ilist; kk++) {
        iel = list[kk] / 3;
        j = list[kk] % 3;
        j1 = MMG5_inxt2[j];
        j2 = MMG5_iprv2[i];
        pt1 = &mesh->tria[iel];
        ip1 = pt1->v[j1];
        ip2 = pt1->v[j2];
        if ( sol->m[ip1] >= MMG5_EPS ) npl = 1;
        else if ( sol->m[ip1] <= -MMG5_EPS ) nmn = 1;

        if ( sol->m[ip2] >= MMG5_EPS ) npl = 1;
        else if ( sol->m[ip2] <= -MMG5_EPS ) nmn = 1;
      }

      if ( npl == 1 && nmn == 0 )
        sol->m[ip] = 100.0*MMG5_EPS;
      else if ( npl == 0 && nmn == 1 )
        sol->m[ip] = 100.0*MMG5_EPS;
    }

  }


  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+nc > 0 )
    fprintf(stdout,"     %8d points snapped, %d corrected\n",ns,nc);

  return 1;
}

/* Check whether the ball of vertex i in tria start is manifold;
 by assumption, i inxt[i] is one edge of the implicit boundary */
int MMG2D_chkmaniball(MMG5_pMesh mesh, int start, char istart) {
  MMG5_pTria         pt;
  int                *adja,k,refstart;
  char               i,i1;

  pt = &mesh->tria[start];
  k = start;
  i = istart;
  refstart = pt->ref;

  /* First travel, while another part of the implicit boundary is not met */
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  while ( k && ( mesh->tria[k].ref == refstart ) );

  /* Case where a boundary is hit: travel in the other sense from start, and make sure
   that a boundary is hit too */
  if ( k == 0 ) {
    k = start;
    i = istart;

    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_iprv2[i];
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_iprv2[i];

    /* Check of the way the point is caught (the left-hand edge is not an external edge) */
    assert ( k );

    do {
      adja = &mesh->adja[3*(k-1)+1];
      i1 = MMG5_iprv2[i];

      k = adja[i1] / 3;
      i = adja[i1] % 3;
      i = MMG5_iprv2[i];
    }
    while ( k && ( mesh->tria[k].ref != refstart ) );

    if ( k == 0 ) return 1;
    else          return 0;

  }

  /* General case: go on travelling until another implicit boundary is met */
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_inxt2[i];
  }
  while ( k && ( mesh->tria[k].ref != refstart ) );

  /* At least 3 boundary segments meeting at p */
  if ( k != start )
    return 0;

  return 1;
}

/* Check whether the resulting two subdomains occupying mesh are manifold */
int MMG2D_chkmanimesh(MMG5_pMesh mesh) {
  MMG5_pTria      pt,pt1;
  int             *adja,k,cnt,iel;
  char            i,i1;
  static char     mmgWarn=0;

  /* First check: check whether one triangle in the mesh has 3 boundary faces */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];
    cnt = 0;
    for (i=0; i<3; i++) {
      iel = adja[i] / 3;

      if (!iel ) {
        cnt++;
        continue;
      }
      else {
        pt1 = &mesh->tria[iel];
        if ( pt1->ref != pt->ref ) cnt++;
      }
    }
    if( cnt == 3 ) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr,"\n  ## Warning: %s: at least 1 triangle with 3 boundary"
                " edges.\n",__func__);
      }
      /* return 0; */
    }
  }

  /* Second check: check whether the configuration is manifold in the ball of each point;
     each vertex on the implicit boundary is caught in such a way that i1 inxt[i1] is one edge of the implicit
     boundary */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      adja = &mesh->adja[3*(k-1)+1];
      iel = adja[i] / 3;

      if (! iel ) continue;
      pt1 = &mesh->tria[iel];
      if ( pt->ref == pt1->ref ) continue;

      i1 = MMG5_inxt2[i];
      if ( !MMG2D_chkmaniball(mesh,k,i1) )
        return 0;
    }
  }

  if ( mesh->info.imprim > 0 || mesh->info.ddebug )
    fprintf(stdout,"  *** Manifold implicit surface.\n");
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the level-set
 *
 * \return 1 if success, 0 otherwise
 *
 * Effective discretization of the 0 level set encoded in sol in the mesh
 *
 */
int MMG2D_cuttri_ls(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash   hash;
  double       v0,v1,s,c[2];
  int          k,ip0,ip1,nb,np,nt,ns,refint,refext,vx[3];
  char         i,i0,i1,ier;

  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Evaluate the number of intersected edges by the 0 level set */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];

      if ( p0->flag && p1->flag ) continue;

      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        nb++;
        if ( !p0->flag ) p0->flag = nb;
        if ( !p1->flag ) p1->flag = nb;
      }
    }
  }
  if ( !nb ) return 1;

  /* Create the intersection points between the edges in the mesh and the 0 level set */
  if ( !MMG5_hashNew(mesh,&hash,nb,2*nb) ) return 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];

      np = MMG5_hashGet(&hash,ip0,ip1);
      if ( np ) continue;

      if ( !MMG2D_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      /* Intersection point between edge p0p1 and the 0 level set */
      s = v0/(v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);

      np = MMG2D_newPt(mesh,c,0);
      if ( !np ) {
        fprintf(stderr,"\n  ## Error: %s: Insufficient memory; abort\n",
          __func__);
        return 0;
      }
      sol->m[np] = 0.0;
      MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }

  /* Proceed to splitting by calling patterns */
  nt  = mesh->nt;
  ns  = 0;
  ier = 1;
  for (k=1; k<=nt; k++) {

    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    pt->flag = 0;

    for (i=0; i<3; i++) {
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      vx[i] = MMG5_hashGet(&hash,ip0,ip1);

      if ( vx[i] ) MG_SET(pt->flag,i);
    }

    switch( pt->flag ) {
      /* 1 edge split -> 0-+ */
      case 1: case 2: case 4:
        ier = MMG2D_split1(mesh,sol,k,vx);
        ns++;
        break;

      /* 2 edge split -> +-- or -++ */
      case 3: case 5: case 6:
        ier = MMG2D_split2(mesh,sol,k,vx);
        ns++;
        break;

      default:
        assert(pt->flag==0);
        break;
    }
    if ( !ier ) return 0;
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;

}

/* Set references to the new triangles */
int MMG2D_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria    pt;
  double        v,v1;
  int           k,ip,ip1,ier,ref,refint,refext;
  char          i,i1,i2,nmn,npl,nz;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    ref = pt->ref;
    nmn = npl = nz = 0;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      v = sol->m[ip];

      if ( v > 0.0 )
        npl++;
      else if ( v < 0.0 )
        nmn++;
      else
        nz++;
    }

    assert(nz < 3);
    ier = MMG2D_isSplit(mesh,ref,&refint,&refext);

    if ( npl ) {
      if ( ier ) {
        assert ( !nmn );
        pt->ref = refext;
      }
    }
    else {
      if ( ier ) {
        assert ( !npl );
        pt->ref = refint;
      }
    }

    /* Set MG_ISO ref at ls edges and at the points of these edges */
    if ( nz == 2 ) {
      for (i=0; i<3; i++) {
        ip  = pt->v[MMG5_inxt2[i]];
        ip1 = pt->v[MMG5_iprv2[i]];
        v   = sol->m[ip];
        v1  = sol->m[ip1];
        if ( v == 0.0 && v1 == 0.0) {
          pt->edg[i]  = MG_ISO;
          pt->tag[i] |= MG_REF;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].ref = MG_ISO;
          mesh->point[pt->v[i2]].ref = MG_ISO;
        }
      }
    }

  }

  return 1;
}

/* Main function of the -ls mode */
int MMG2D_mmg2d6(MMG5_pMesh mesh, MMG5_pSol sol) {
  double *tmp;
  int k;

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION\n");

  /* Work only with the 0 level set */
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /* Allocate memory for tmp */
  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                printf("  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double,return 0);

  /* Snap values of the level set function which are very close to 0 to 0 exactly */
  if ( !MMG2D_snapval(mesh,sol,tmp) ) {
    fprintf(stderr,"\n  ## Wrong input implicit function. Exit program.\n");
    return 0;
  }

  MMG5_DEL_MEM(mesh,tmp);

  /* Creation of adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* No need to keep adjacencies from now on */
  MMG5_DEL_MEM(mesh,mesh->adja);

  /* Transfer the boundary edge references to the triangles */
  if ( !MMG2D_assignEdge(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* Reset the MG_ISO field everywhere it appears */
  if ( !MMG2D_resetRef(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    return 0;
  }

  /* Effective splitting of the crossed triangles */
  if ( !MMG2D_cuttri_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in cutting triangles. Exit program.\n");
    return 0;
  }

  /* Set references on the interior / exterior triangles*/
  if ( !MMG2D_setref_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Creation of adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Check that the resulting mesh is manifold */
  if ( !MMG2D_chkmanimesh(mesh) ) {
    fprintf(stderr,"\n  ## No manifold resulting situation. Exit program.\n");
    return 0;
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);
  sol->np = 0;

  if ( mesh->info.mat )
    MMG5_SAFE_FREE( mesh->info.mat );

  return 1;
}
