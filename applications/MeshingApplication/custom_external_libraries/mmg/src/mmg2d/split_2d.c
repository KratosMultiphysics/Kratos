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
 * \file mmg2d/split_2d.c
 * \brief Functions for splitting.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmg2d_private.h"
#include "mmg2dexterns_private.h"

extern uint8_t ddb;

/**
 * \param mesh pointer to the mesh
 * \param met pointer to the metric
 * \param k triangle index
 * \param i local index of the edge to split
 *
 * \return 1 if we can split, 0 if not, -1 if fail.
 *
 * Check whether splitting of edge i in tria k is possible and return the newly created point;
 * possibly perform a dichotomy to find the latest valid position for the point.
 *
 */
MMG5_int MMG2D_chkspl(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int8_t i) {
  MMG5_pTria           pt,pt1,pt0;
  MMG5_pPoint          p1,p2,ppt;
  double               mid[2],o[2],no[2],calnew,caltmp,tp,to,t,calseuil;
  MMG5_int             ip,jel,*adja,npinit;
  int                  it,maxit;
  const double         s = 0.5;
  int8_t               i1,i2,j,j1,j2,ier,isv;
  assert ( met );

  calseuil = 1e-4 / MMG2D_ALPHAD;
  npinit = mesh->np;


  pt  = &mesh->tria[k];
  pt0 = &mesh->tria[0];
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];

  p1 = &mesh->point[pt->v[i1]];
  p2 = &mesh->point[pt->v[i2]];

  adja = &mesh->adja[3*(k-1)+1];

  jel  = adja[i] / 3;
  j    = adja[i] % 3;
  j1   = MMG5_inxt2[j];
  j2   = MMG5_iprv2[j];

  /* Midpoint of edge i */
  mid[0] = s*(p1->c[0]+p2->c[0]);
  mid[1] = s*(p1->c[1]+p2->c[1]);

  /* If the splitted edge is not geometric, the new point is simply its midpoint */
  if ( !MG_EDG(pt->tag[i]) ) {
    ip = MMG2D_newPt(mesh,mid,0);
    if ( !ip ) {
      /* reallocation of point table */
      MMG2D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                           printf("  ## Error: unable to allocate a new point.\n");
                           MMG5_INCREASE_MEM_MESSAGE();
                           do {
                             MMG2D_delPt(mesh,mesh->np);
                           } while ( mesh->np>npinit );return -1;,
                           mid,pt->tag[i]);

    }
    /* If there is a metric in the mesh, interpolate it at the new point */
    if ( met->m )
      MMG2D_intmet(mesh,met,k,i,ip,s);

    ppt = &mesh->point[ip];
    if ( pt->tag[i] ) ppt->tag = pt->tag[i];
    if ( pt->edg[i] ) ppt->ref = pt->edg[i];

    /* Check quality of the four new elements */
    calnew = DBL_MAX;
    memcpy(pt0,pt,sizeof(MMG5_Tria));
    pt0->v[i2] = ip;

    caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
    calnew = MG_MIN(calnew,caltmp);

    pt0->v[i1] = ip; pt0->v[i2] = pt->v[i2];
    caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
    calnew = MG_MIN(calnew,caltmp);

    if ( jel ) {
      pt1 = &mesh->tria[jel];
      memcpy(pt0,pt1,sizeof(MMG5_Tria));
      pt0->v[j1] = ip;
      caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
      calnew = MG_MIN(calnew,caltmp);

      pt0->v[j1] = pt1->v[j1] ; pt0->v[j2] = ip;
      caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
      calnew = MG_MIN(calnew,caltmp);
    }

    /* Delete point and abort splitting if one of the created triangles
       has very bad quality.
       MMG5_EPSOK is not sufficient :
       we were created very bad element and were not able to delete them */
    if ( (calnew < calseuil)  ) {
      MMG2D_delPt(mesh,ip);
      return 0;
    }
  }
  /* Otherwise, the new point is inserted on the underlying curve to the edge;
     a dichotomy is applied to find the largest distance to the edge that yields an admissible configuration */
  else {
    ier = MMG2D_bezierCurv(mesh,k,i,s,o,no);
    if ( !ier ) return 0;

    ip  = MMG2D_newPt(mesh,o,pt->tag[i]);
    if ( !ip ) {
      /* reallocation of point table */
      MMG2D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                           printf("  ## Error: unable to allocate a new point.\n");
                           MMG5_INCREASE_MEM_MESSAGE();
                           do {
                             MMG2D_delPt(mesh,mesh->np);
                           } while ( mesh->np>npinit ); return -1;,
                           o,pt->tag[i]);
    }
    if ( met->m )
      MMG2D_intmet(mesh,met,k,i,ip,s);

    ppt = &mesh->point[ip];
    if ( pt->tag[i] ) ppt->tag = pt->tag[i];
    if ( pt->edg[i] ) ppt->ref = pt->edg[i];

    ppt->n[0] = no[0];
    ppt->n[1] = no[1];

    isv   = 0;
    it    = 0;
    maxit = 5;
    tp    = 1.0;
    t     = 1.0;
    to    = 0.0;

    do {
      ppt->c[0] = mid[0] + t*(o[0] - mid[0]);
      ppt->c[1] = mid[1] + t*(o[1] - mid[1]);

      /* Check quality of the four new elements */
      calnew = DBL_MAX;
      memcpy(pt0,pt,sizeof(MMG5_Tria));
      pt0->v[i2] = ip;
      caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
      calnew = MG_MIN(calnew,caltmp);

      pt0->v[i1] = ip; pt0->v[i2] = pt->v[i2];
      caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
      calnew = MG_MIN(calnew,caltmp);

      if ( jel ) {
        pt1 = &mesh->tria[jel];
        memcpy(pt0,pt1,sizeof(MMG5_Tria));
        pt0->v[j1] = ip;
        caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
        calnew = MG_MIN(calnew,caltmp);

        pt0->v[j1] = pt1->v[j1] ; pt0->v[j2] = ip;
        caltmp = MMG2D_ALPHAD*MMG2D_caltri(mesh,met,pt0);
        calnew = MG_MIN(calnew,caltmp);
      }

      ier = ( calnew > MMG5_EPSOK );
      if ( ier ) {
        isv = 1;
        to = t;
        if ( t == tp ) break;
      }
      else
        tp = t;

      /* If no admissible position has been found, do the last iteration with the midpoint m */
      if ( (it == maxit-2) && !isv )
        t = 0.0;
      else
        t = 0.5*(to+tp);
    }
    while ( ++it < maxit );

    /* One satisfying position has been found: to */
    if ( isv ) {
      ppt->c[0] = mid[0] + to*(o[0] - mid[0]);
      ppt->c[1] = mid[1] + to*(o[1] - mid[1]);
    }
    /* No satisfying position has been found */
    else {
      MMG2D_delPt(mesh,ip);
      return 0;
    }
  }

  return ip;
}

/**
 * \parma mesh pointer to the mesh
 * \param k index of the tria to split
 * \param i local index of the edge to split
 * \param ip global index of the new point
 *
 * \return 1 if success, 0 if fail
 *
 * Effective splitting of edge i in tria k: point ip is introduced and the
 * adjacency structure in the mesh is preserved
 *
 */
int MMG2D_split1b(MMG5_pMesh mesh,MMG5_int k,int8_t i,MMG5_int ip) {
  MMG5_pTria         pt,pt1;
  MMG5_int           *adja,iel,jel,kel,mel;
  int8_t             i1,i2,m,j,j1,j2;

  iel = MMG2D_newElt(mesh);
  if ( !iel ) {
    MMG2D_TRIA_REALLOC(mesh,iel,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
  }

  pt = &mesh->tria[k];
  pt->flag = 0;
  pt->base = mesh->base;

  i1 = MMG5_inxt2[i];
  i2 = MMG5_iprv2[i];

  adja = &mesh->adja[3*(k-1)+1];
  jel  = adja[i] / 3;
  j    = adja[i] % 3;

  pt1 = &mesh->tria[iel];
  memcpy(pt1,pt,sizeof(MMG5_Tria));
  memcpy(&mesh->adja[3*(iel-1)+1],&mesh->adja[3*(k-1)+1],3*sizeof(MMG5_int));

  /* Update both triangles */
  pt->v[i2]  = ip;
  pt1->v[i1] = ip;

  pt->tag[i1] = MG_NOTAG;
  pt->edg[i1] = 0;

  pt1->tag[i2] = MG_NOTAG;
  pt1->edg[i2] = 0;

  /* Update adjacencies */
  mel = adja[i1] / 3;
  m   = adja[i1] % 3;
  mesh->adja[3*(k-1)+1+i1] = 3*iel+i2;
  mesh->adja[3*(iel-1)+1+i2] = 3*k+i1;
  if ( mel )
    mesh->adja[3*(mel-1)+1+m] = 3*iel+i1;

  if ( jel ) {
    kel = MMG2D_newElt(mesh);
    if ( !kel ) {
      MMG2D_TRIA_REALLOC(mesh,kel,mesh->gap,
                          printf("  ## Error: unable to allocate a new element.\n");
                          MMG5_INCREASE_MEM_MESSAGE();
                          printf("  Exit program.\n");return 0);
    }

    pt  = &mesh->tria[jel];
    pt1 = &mesh->tria[kel];
    j1 = MMG5_inxt2[j];
    j2 = MMG5_iprv2[j];

    pt->flag = 0;
    pt->base = mesh->base;

    memcpy(pt1,pt,sizeof(MMG5_Tria));
    memcpy(&mesh->adja[3*(kel-1)+1],&mesh->adja[3*(jel-1)+1],3*sizeof(MMG5_int));

    /* Update triangles */
    pt->v[j1]    = ip;
    pt1->v[j2]   = ip;
    pt->tag[j2]  = MG_NOTAG;
    pt->edg[j2]  = 0;
    pt1->tag[j1] = MG_NOTAG;
    pt1->edg[j1] = 0;

    /* Update adjacencies */
    adja = &mesh->adja[3*(jel-1)+1];
    mel  = adja[j2] / 3;
    m    = adja[j2] % 3;
    mesh->adja[3*(jel-1)+1+j2] = 3*kel+j1;
    mesh->adja[3*(kel-1)+1+j1] = 3*jel+j2;
    if ( mel )
      mesh->adja[3*(mel-1)+1+m] = 3*kel+j2;

    mesh->adja[3*(iel-1)+1+i] = 3*kel+j;
    mesh->adja[3*(kel-1)+1+j] = 3*iel+i;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the metric
 * \param k triangle index
 * \param vx list of new point indices for each edge
 *
 * \return 0 if fail, 1 if success
 *
 * Simulate the split of one edge in triangle k
 *
 */
int MMG2D_split1_sim(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int k, MMG5_int vx[3]) {
  MMG5_pTria  pt,pt0;
  double      cal;
  uint8_t     tau[3];

  pt = &mesh->tria[k];
  pt0 = &mesh->tria[0];
  memcpy(pt0,pt,sizeof(MMG5_Tria));

  /* Set permutation from the reference configuration (case 1: edge 0 is splitted) to the actual one */
  tau[0] = 0; tau[1] = 1; tau[2] = 2;

  switch ( pt->flag ) {
  case 2:
    tau[0] = 1; tau[1] = 2; tau[2] = 0;
    break;

  case 4:
    tau[0] = 2; tau[1] = 0; tau[2] = 1;
    break;
  }

  pt0->v[tau[2]] = vx[tau[0]];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  pt0->v[tau[2]] = pt->v[tau[2]];
  pt0->v[tau[1]] = vx[tau[0]];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the metric
 * \param k triangle index
 * \param vx list of new point indices for each edge
 *
 * \return 0 if fail, 1 if success
 *
 * Split 1 edge of triangle k
 *
 */
int MMG2D_split1(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int k, MMG5_int vx[3]) {
  MMG5_pTria   pt,pt1;
  MMG5_pPoint  p0;
  MMG5_int     iel;
  uint8_t      tau[3];

  pt = &mesh->tria[k];

  /* Set permutation from the reference configuration (case 1: edge 0 is splitted) to the actual one */
  tau[0] = 0; tau[1] = 1; tau[2] = 2;

  switch ( pt->flag ) {
  case 2:
    tau[0] = 1; tau[1] = 2; tau[2] = 0;
    break;

  case 4:
    tau[0] = 2; tau[1] = 0; tau[2] = 1;
    break;
  }

  pt->flag = 0;

  /* Update of point references */
  p0 = &mesh->point[vx[tau[0]]];

  if ( pt->edg[tau[0]] > 0 )
    p0->ref = pt->edg[tau[0]];

  iel = MMG2D_newElt(mesh);
  if ( !iel ) {
    MMG2D_TRIA_REALLOC(mesh,iel,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
    pt = &mesh->tria[k];
  }
  pt1 = &mesh->tria[iel];
  memcpy(pt1,pt,sizeof(MMG5_Tria));

  /* Generic formulation for the split of one edge */
  /* Update of vertices */
  pt->v[tau[2]] = vx[tau[0]];
  pt1->v[tau[1]] = vx[tau[0]];

  /* Update of edge references and tags*/
  pt->tag[tau[1]] = MG_NOTAG;
  pt->edg[tau[1]] = 0;

  pt1->tag[tau[2]] = MG_NOTAG;
  pt1->edg[tau[2]] = 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the metric
 * \param k triangle index
 * \param vx list of new point indices for each edge
 *
 * \return 0 if fail, 1 if success
 *
 * Simulate the split of two edges in triangle k
 *
 */
int MMG2D_split2_sim(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int k, MMG5_int vx[3]) {
  MMG5_pTria  pt,pt0;
  double      cal;
  uint8_t     tau[3];

  pt = &mesh->tria[k];
  pt0 = &mesh->tria[0];
  memcpy(pt0,pt,sizeof(MMG5_Tria));

  /* Set permutation from the reference configuration (case 6: edges 1,2 are splitted) to the actual one */
  tau[0] = 0; tau[1] = 1; tau[2] = 2;

  switch ( pt->flag ) {
  case 5:
    tau[0] = 1; tau[1] = 2; tau[2] = 0;
    break;

  case 3:
    tau[0] = 2; tau[1] = 0; tau[2] = 1;
    break;
  }

  pt0->v[tau[1]] = vx[tau[2]] ; pt0->v[tau[2]] = vx[tau[1]];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  pt0->v[tau[1]] = pt->v[tau[1]] ; pt0->v[tau[2]] = pt->v[tau[2]];
  pt0->v[tau[0]] = vx[tau[2]];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  pt0->v[tau[0]] = vx[tau[1]] ; pt0->v[tau[1]] = vx[tau[2]];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the metric
 * \param k triangle index
 * \param vx list of new point indices for each edge
 *
 * \return 0 if fail, 1 if success
 *
 * Split 2 edges of triangle k
 *
 */
int MMG2D_split2(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int k, MMG5_int vx[3]) {
  MMG5_pTria  pt,pt1,pt2;
  MMG5_pPoint p1,p2;
  MMG5_int    iel,jel;
  uint8_t     tau[3];

  pt = &mesh->tria[k];

  /* Set permutation from the reference configuration (case 6: edges 1,2 are splitted) to the actual one */
  tau[0] = 0; tau[1] = 1; tau[2] = 2;

  switch ( pt->flag ) {
  case 5:
    tau[0] = 1; tau[1] = 2; tau[2] = 0;
    break;

  case 3:
    tau[0] = 2; tau[1] = 0; tau[2] = 1;
    break;
  }

  pt->flag = 0;

  /* Update of point references */
  p1 = &mesh->point[vx[tau[1]]];
  p2 = &mesh->point[vx[tau[2]]];

  if ( pt->edg[tau[1]] > 0 )
    p1->ref = pt->edg[tau[1]];

  if ( pt->edg[tau[2]] > 0 )
    p2->ref = pt->edg[tau[2]];

  iel = MMG2D_newElt(mesh);
  if ( !iel ) {
    MMG2D_TRIA_REALLOC(mesh,iel,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
    pt = &mesh->tria[k];
  }

  jel = MMG2D_newElt(mesh);
  if ( !jel ) {
    MMG2D_TRIA_REALLOC(mesh,jel,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
    pt = &mesh->tria[k];
  }

  pt1 = &mesh->tria[iel];
  pt2 = &mesh->tria[jel];
  memcpy(pt1,pt,sizeof(MMG5_Tria));
  memcpy(pt2,pt,sizeof(MMG5_Tria));


  /* Generic formulation for the split of two edges */
  /* Update of vertices */
  pt->v[tau[1]] = vx[tau[2]] ; pt->v[tau[2]] = vx[tau[1]];
  pt1->v[tau[0]] = vx[tau[2]];
  pt2->v[tau[0]] = vx[tau[1]]; pt2->v[tau[1]] = vx[tau[2]];

  /* Update of edge references and tags*/
  pt->tag[tau[0]] = MG_NOTAG;
  pt->edg[tau[0]] = 0;

  pt1->tag[tau[1]] = MG_NOTAG;
  pt1->edg[tau[1]] = 0;

  pt2->tag[tau[0]] = MG_NOTAG;   pt2->tag[tau[2]] = MG_NOTAG;
  pt2->edg[tau[0]] = MG_NOTAG;   pt2->edg[tau[2]] = MG_NOTAG;

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the metric
 * \param k triangle index
 * \param vx list of new point indices for each edge
 *
 * \return 0 if fail, 1 if success
 *
 * Simulate the split of three edges in triangle k
 *
 */
int MMG2D_split3_sim(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int k, MMG5_int vx[3]) {
  MMG5_pTria         pt,pt0;
  double             cal;

  pt = &mesh->tria[k];
  pt0 = &mesh->tria[0];
  memcpy(pt0,pt,sizeof(MMG5_Tria));

  pt0->v[1] = vx[2] ; pt0->v[2] = vx[1];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  pt0->v[0] = vx[2] ; pt0->v[1] = pt->v[1]; pt0->v[2] = vx[0];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  pt0->v[0] = vx[1] ; pt0->v[1] = vx[0]; pt0->v[2] = pt->v[2];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  pt0->v[1] = vx[2]; pt0->v[2] = vx[0];
  cal = MMG2D_quickcal(mesh,pt0);
  if ( cal < MMG5_EPSD )  return 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the metric
 * \param k triangle index
 * \param vx list of new point indices for each edge
 *
 * \return 0 if fail, 1 if success
 *
 * Split the three edges of triangle k
 *
 */
int MMG2D_split3(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int k, MMG5_int vx[3]) {
  MMG5_pTria          pt,pt1,pt2,pt3;
  MMG5_pPoint         p0,p1,p2;
  MMG5_int            iel,jel,kel;

  pt = &mesh->tria[k];
  pt->flag = 0;

  /* Update of point references */
  p0 = &mesh->point[vx[0]];
  p1 = &mesh->point[vx[1]];
  p2 = &mesh->point[vx[2]];

  if ( pt->edg[0] > 0 )
    p0->ref = pt->edg[0];

  if ( pt->edg[1] > 0 )
    p1->ref = pt->edg[1];

  if ( pt->edg[2] > 0 )
    p2->ref = pt->edg[2];

  iel = MMG2D_newElt(mesh);
  if ( !iel ) {
    MMG2D_TRIA_REALLOC(mesh,iel,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);

    pt = &mesh->tria[k];
  }

  jel = MMG2D_newElt(mesh);

  if ( !jel ) {
    MMG2D_TRIA_REALLOC(mesh,jel,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
    pt = &mesh->tria[k];
  }

  kel = MMG2D_newElt(mesh);

  if ( !kel ) {
    MMG2D_TRIA_REALLOC(mesh,kel,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
    pt = &mesh->tria[k];
  }

  pt1 = &mesh->tria[iel];
  pt2 = &mesh->tria[jel];
  pt3 = &mesh->tria[kel];
  memcpy(pt1,pt,sizeof(MMG5_Tria));
  memcpy(pt2,pt,sizeof(MMG5_Tria));
  memcpy(pt3,pt,sizeof(MMG5_Tria));

  /* Update of vertices */
  pt->v[1] = vx[2] ; pt->v[2] = vx[1];
  pt1->v[0] = vx[2] ; pt1->v[2] = vx[0];
  pt2->v[0] = vx[1]; pt2->v[1] = vx[0];
  pt3->v[0] = vx[1] ; pt3->v[1] = vx[2] ; pt3->v[2] = vx[0];

  /* Update of tags and references */
  pt->tag[0] = MG_NOTAG;
  pt->edg[0] = 0;

  pt1->tag[1] = MG_NOTAG;
  pt1->edg[1] = 0;

  pt2->tag[2] = MG_NOTAG;
  pt2->edg[2] = 0;

  pt3->tag[0] = pt3->tag[1] = pt3->tag[2] = MG_NOTAG;
  pt3->edg[0] = pt3->edg[1] = pt3->edg[2] = 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param k index of the tria to split
 * \param ip global index of the new point
 *
 * \return 1 if success, 0 if fail
 *
 * Insert the point ip inside the tria k
 *
 */
int MMG2D_splitbar(MMG5_pMesh mesh,MMG5_int k,MMG5_int ip) {
  MMG5_pTria         pt,pt0,pt1,pt2;
  MMG5_pPoint        p0,p1,p2,ppt;
  MMG5_int           *adja,iel1,iel2,jel0,jel2;
  MMG5_int           ip0,ip1,ip2;
  int8_t             j2,j0;
  double             cal,calseuil;

  pt  = &mesh->tria[k];
  pt0 = &mesh->tria[0];
  ppt = &mesh->point[ip];
  ip0 = pt->v[0];
  p0  = &mesh->point[ip0];
  ip1 = pt->v[1];
  p1  = &mesh->point[ip1];
  ip2 = pt->v[2];
  p2 = &mesh->point[ip2];

  calseuil = MMG5_EPSOK ;

  /* Check quality of the three new elements */
  cal = MMG2D_quickarea(ppt->c,p1->c,p2->c);
  if ( (cal < calseuil)  ) {
     return 0;
  }

  cal = MMG2D_quickarea(p0->c,ppt->c,p2->c);
  if ( (cal < calseuil)  ) {
      return 0;
  }
  pt0->v[0] = ip0;
  cal = MMG2D_quickarea(p0->c,p1->c,ppt->c);
  if ( (cal < calseuil)  ) {
      return 0;
  }

  iel1 = MMG2D_newElt(mesh);
  if ( !iel1 ) {
    MMG2D_TRIA_REALLOC(mesh,iel1,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
  }
  iel2 = MMG2D_newElt(mesh);
  if ( !iel2 ) {
    MMG2D_TRIA_REALLOC(mesh,iel2,mesh->gap,
                        printf("  ## Error: unable to allocate a new element.\n");
                        MMG5_INCREASE_MEM_MESSAGE();
                        printf("  Exit program.\n");return 0);
  }

  pt->flag = 0;
  pt->base = mesh->base;

  adja = &mesh->adja[3*(k-1)+1];
  jel0  = adja[0] / 3;
  j0    = adja[0] % 3;
#ifndef NDEBUG
  int8_t jel1 = adja[1] / 3;
  int8_t j1   = adja[1] % 3;
#endif
  jel2  = adja[2] / 3;
  j2    = adja[2] % 3;

  pt1 = &mesh->tria[iel1];
  memcpy(pt1,pt,sizeof(MMG5_Tria));
  memcpy(&mesh->adja[3*(iel1-1)+1],&mesh->adja[3*(k-1)+1],3*sizeof(MMG5_int));
  pt2 = &mesh->tria[iel2];
  memcpy(pt2,pt,sizeof(MMG5_Tria));
  memcpy(&mesh->adja[3*(iel2-1)+1],&mesh->adja[3*(k-1)+1],3*sizeof(MMG5_int));

  /* Update the three triangles */
  pt->v[1]  = ip;
  pt1->v[2] = ip;
  pt2->v[0] = ip;

  pt->tag[1] = MG_NOTAG;
  pt->edg[1] = 0;

  pt1->tag[2] = MG_NOTAG;
  pt1->edg[2] = 0;

  pt2->tag[0] = MG_NOTAG;
  pt2->edg[0] = 0;

  /* Update external adjacencies */
#ifndef NDEBUG
  assert(mesh->adja[3*(k-1)+1+1] == 3*jel1+j1);
  if ( jel1 ) {
    assert(mesh->adja[3*(jel1-1)+1+j1] == 3*k+1);
  }
#endif

  mesh->adja[3*(iel1-1)+1+2] = 3*jel2+j2;
  if ( jel2 )
    mesh->adja[3*(jel2-1)+1+j2] = 3*iel1+2;

  mesh->adja[3*(iel2-1)+1+0] = 3*jel0+j0;
  if ( jel0 )
    mesh->adja[3*(jel0-1)+1+j0] = 3*iel2+0;

  /*update internal adjacencies*/
  mesh->adja[3*(k-1)+1+0] = 3*iel2+1;
  mesh->adja[3*(iel2-1)+1+1] = 3*k+0;

  mesh->adja[3*(k-1)+1+2] = 3*iel1+1;
  mesh->adja[3*(iel1-1)+1+1] = 3*k+2;

  mesh->adja[3*(iel1-1)+1+0] = 3*iel2+2;
  mesh->adja[3*(iel2-1)+1+2] = 3*iel1+0;

  return 1;
}
