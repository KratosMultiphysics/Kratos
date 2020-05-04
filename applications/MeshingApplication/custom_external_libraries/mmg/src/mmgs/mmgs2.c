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
 * \file mmgs/mmgs2.c
 * \brief Create implicit surface in mesh.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \param start index of the starting tria
 * \param istart local index (inside the tria \a start) of the vertex that we check.
 * \return 1 if success, 0 if fail
 *
 * Check whether snapping the value of vertex \a istart of \a start to 0 exactly
 * leads to a non manifold situation.
 *
 * \warning: we assume that the triangle \a start has vertex \a istart
 * with value 0 and the other two with changing values.
 *
 */
static int
MMGS_ismaniball(MMG5_pMesh mesh, MMG5_pSol sol, int start, char istart) {
  MMG5_pTria       pt;
  double           v1, v2;
  int              *adja,k,ip1,ip2,end1;
  char             i,i1,smsgn;
  static char      mmgWarn=0;

  k = start;
  i = MMG5_inxt2[istart];

  /* First loop: stop if an external boundary, or a change in signs (or a 0) is
     met recall that MG_SMGSGN(a,b) = 1 provided a*b >0 */
  do{
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i] / 3;
    i1 = adja[i] % 3;
    i = MMG5_iprv2[i1];

    if ( k==0 ) break;

    pt = &mesh->tria[k];
    ip1 = pt->v[i1];
    ip2 = pt->v[i];

    v1 = sol->m[ip1]-mesh->info.ls;
    v2 = sol->m[ip2]-mesh->info.ls;

    smsgn = MG_SMSGN(v1,v2) ? 1 : 0;
  }
  while ( smsgn && k!=start );

  assert(k!=start); //unexpected case

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

    v1 = sol->m[ip1]-mesh->info.ls;
    v2 = sol->m[ip2]-mesh->info.ls;

    smsgn = MG_SMSGN(v1,v2) ? 1 : 0;
  }
  while ( smsgn && k!=start );

  assert(k!=start); //unexpected case

  /* If first stop was due to an external boundary, the second one must too
     (end1 = k = 0);
     else, the final triangle for the first travel must be that
     of the second one */
  if ( k != end1 ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: triangle %d, point %d; unsnap \n",
              __func__,start,mesh->tria[start].v[istart]);
    }
    return 0;
  }
  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set function.
 * \return 1 if success, 0 if fail.
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
static
int MMGS_snpval_ls(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria    pt;
  MMG5_pPoint   p0;
  double        *tmp,v1,v2;
  int           k,nc,ns,ip,ip1,ip2;
  char          i;

  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double,return 0);

  /* create tetra adjacency */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: hashing problem (1). Exit program.\n",
            __func__);
    return 0;
  }

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]-mesh->info.ls) < MMG5_EPS ) {
      if ( mesh->info.ddebug )
        fprintf(stderr,"  ## Warning: %s: snapping value %d; "
                "previous value : %E\n",__func__,k,fabs(sol->m[k]));

      tmp[k] = ( fabs(sol->m[k]-mesh->info.ls) < MMG5_EPSD ) ?
        (mesh->info.ls-100.0*MMG5_EPS) : sol->m[k];
      p0->flag = 1;
      sol->m[k] = mesh->info.ls;
      ns++;
    }
  }

  /* Check that the snapping process has not led to a nonmanifold situation */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ip1 = pt->v[MMG5_inxt2[i]];
      ip2 = pt->v[MMG5_iprv2[i]];

      p0 = &mesh->point[ip];
      v1 = sol->m[ip1]-mesh->info.ls;
      v2 = sol->m[ip2]-mesh->info.ls;
      if ( p0->flag && !(MG_SMSGN(v1,v2)) ) {
        if ( !MMGS_ismaniball(mesh,sol,k,i) ) {
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
  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,tmp);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of starting tria.
 * \param istart local index of point that we check (in tria \a start)
 * \return 1 if the ball is manifold, 0 otherwise.
 *
 * Check whether the ball of vertex i in tria start is manifold;
 *
 * \warning i inxt[i] is one edge of the implicit boundary.
 *
 */
int MMGS_chkmaniball(MMG5_pMesh mesh, int start, char istart) {
  MMG5_pTria         pt;
  int                *adja,k;
  char               i,i1;

  pt = &mesh->tria[start];
  k = start;
  i = istart;

  i1 = MMG5_iprv2[i];
  assert( MG_EDG(pt->tag[i1]) && (pt->edg[i1]==MG_ISO) );

  /* First travel, while another part of the implicit boundary is not met */
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;

    if ( !k || mesh->tria[k].edg[i]==MG_ISO ) break;

    i = MMG5_inxt2[i];
  }
  while ( k!=start );

  assert(k!=start); //unexpected case

  /* Case where a boundary is hit: travel in the other sense from start, and make sure
   that a boundary is hit too */
  if ( k == 0 ) {
    k = start;
    i = istart;

    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_iprv2[i];
    k = adja[i1] / 3;

    /* Check of the way the point is caught (the left-hand edge is not an external edge) */
    assert ( k );

    i = adja[i1] % 3;
    i = MMG5_iprv2[i];

    do {
      adja = &mesh->adja[3*(k-1)+1];
      i1 = MMG5_iprv2[i];

      k = adja[i1] / 3;
      i = adja[i1] % 3;

      if ( (!k) || mesh->tria[k].edg[i]==MG_ISO ) break;

      i = MMG5_iprv2[i];
    }
    while ( k!=start );

    assert(k!=start); //unexpected case

    return !k;
  }

  /* General case: go on travelling until another implicit boundary is met */
  i = MMG5_inxt2[i];
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;

    if ( (!k) || mesh->tria[k].edg[i]==MG_ISO ) break;

    i = MMG5_inxt2[i];
  }
  while ( k!=start );

  /* At least 3 boundary segments meeting at p */
  if ( k != start )
    return 0;

  return 1;
}

/**
 * \param mesh pointer toward the mesh.
 * \return 1 if the mesh is manifold, 0 otherwise.
 *
 * Check whether the resulting two subdomains occupying mesh are manifold.
 *
 */
static
int MMGS_chkmanimesh(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  int             *adja,k,cnt,iel;
  char            i,i1;
  static char     mmgWarn0 = 0;


  /* First check: check whether one triangle in the mesh has 3 boundary faces */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    cnt = 0;
    adja = &mesh->adja[3*(k-1)+1];
    for (i=0; i<3; i++) {
      iel = adja[i] / 3;

      if (!iel ) {
        cnt++;
        continue;
      }
      else {
        if ( pt->edg[i] == MG_ISO ) cnt++;
      }
    }
    if( cnt == 3 ) {
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Warning: %s: at least 1 triangle with 3 boundary"
                " edges.\n",__func__);
      }
      return 0;
    }
  }

  /* Second check: check whether the configuration is manifold in the ball of
     each point; each vertex on the implicit boundary is caught in such a way
     that i1 inxt[i1] is one edge of the implicit boundary */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      adja = &mesh->adja[3*(k-1)+1];
      iel = adja[i] / 3;

      if ( (!iel) || (pt->edg[i] != MG_ISO) ) continue;

      i1 = MMG5_inxt2[i];
      if ( !MMGS_chkmaniball(mesh,k,i1) )
        return 0;
    }
  }

  if ( mesh->info.imprim > 0 || mesh->info.ddebug )
    fprintf(stdout,"  *** Manifold implicit surface.\n");

  return 1;
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
static int MMGS_cuttri_ls(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash   hash;
  double       c[3],v0,v1,s;
  int          vx[3],nb,k,ip0,ip1,np,ns,nt,ier;
  char         ia;
  /* reset point flags and h */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* compute the number nb of intersection points on edges */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    for (ia=0; ia<3; ia++) {
      ip0 = pt->v[MMG5_inxt2[ia]];
      ip1 = pt->v[MMG5_iprv2[ia]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      v0  = sol->m[ip0]-mesh->info.ls;
      v1  = sol->m[ip1]-mesh->info.ls;
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
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
  if ( ! nb )  return 1;

  /* Create intersection points at 0 isovalue and set flags to trias */
  if ( !MMG5_hashNew(mesh,&hash,nb,3*nb) ) return 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (ia=0; ia<3; ia++) {
      ip0 = pt->v[MMG5_inxt2[ia]];
      ip1 = pt->v[MMG5_iprv2[ia]];
      np  = MMG5_hashGet(&hash,ip0,ip1);
      if ( np )  continue;

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0]-mesh->info.ls;
      v1 = sol->m[ip1]-mesh->info.ls;
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      s = v0 / (v0-v1);

      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

      np = MMGS_newPt(mesh,c,NULL);
      if ( !np ) {
        MMGS_POINT_REALLOC(mesh,sol,np,0.2,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            return 0
                            ,c,NULL);
      }
      sol->m[np] = mesh->info.ls;
      MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
    }
  }

  /* Proceed to splitting, according to flags to tris */
  nt  = mesh->nt;
  ns  = 0;
  ier = 1;
  for (k=1; k<=nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = 0;
    memset(vx,0,3*sizeof(int));
    for (ia=0; ia<3; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_inxt2[ia]],pt->v[MMG5_iprv2[ia]]);
      if ( vx[ia] ) {
        MG_SET(pt->flag,ia);
      }
    }
    switch (pt->flag) {
    case 1: /* 1 edge split */
      ier = MMGS_split1(mesh,sol,k,0,vx);
      ns++;
      break;

    case 2: /* 1 edge split */
      ier = MMGS_split1(mesh,sol,k,1,vx);
      ns++;
      break;

    case 4: /* 1 edge split */
      ier = MMGS_split1(mesh,sol,k,2,vx);
      ns++;
      break;

    case 3: case 5: case 6: /* 2 edges split */
      ier = MMGS_split2(mesh,sol,k,vx);
      ns++;
      break;

    case 7: /* 3 edges splitted */
      ier =MMGS_split3(mesh,sol,k,vx);
      ns++;
      break;

    default :
      assert(pt->flag == 0);
      break;
    }
    if ( !ier ) return 0;
  }
  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7d splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \return 1.
 *
 * Set references to tris according to the sign of the level set function.
 *
 */
static int MMGS_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTria   pt;
  double       v,v1;
  int          k,ip,ip1;
  char         nmns,npls,nz,i;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    nmns = npls = nz = 0;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      v = sol->m[ip]-mesh->info.ls;
      if ( v > 0.0 )
        npls++;
      else if ( v < 0.0 )
        nmns++;
      else
        nz ++;
    }
    assert(nz < 3);

    if ( mesh->info.iso != 2 ) {
      /* Keep the initial triangle references of the mesh */
      if ( npls ) {
        assert(!nmns);
        pt->ref = MG_PLUS;
      }
      else {
        assert(nmns);
        pt->ref = MG_MINUS;
      }
    }

    // Set MG_ISO ref at ls edges
    if ( nz == 2 ) {
      for (i=0; i<3; i++) {
        ip  = pt->v[MMG5_inxt2[i]];
        ip1 = pt->v[MMG5_iprv2[i]];
        v   = sol->m[ip] -mesh->info.ls;
        v1  = sol->m[ip1]-mesh->info.ls;
        if ( v == 0.0 && v1 == 0.0) {
          pt->edg[i]  = MG_ISO;
          pt->tag[i] |= MG_REF;
        }
      }
    }

  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solution structure
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int MMGS_mmgs2(MMG5_pMesh mesh,MMG5_pSol sol) {

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION\n");

  /* Snap values of level set function if need be, then discretize it */
  if ( !MMGS_snpval_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    return 0;
  }

  /* Check the initial mesh adjacencies */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  MMG5_DEL_MEM(mesh,mesh->adja);

  if ( !MMGS_cuttri_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
    return 0;
  }

  if ( !MMGS_setref_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Creation of adjacency relations in the mesh */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* Check that the resulting mesh is manifold */
  if ( !MMGS_chkmanimesh(mesh) ) {
    fprintf(stderr,"\n  ## No manifold resulting situation. Exit program.\n");
    return 0;
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);

  return 1;
}
