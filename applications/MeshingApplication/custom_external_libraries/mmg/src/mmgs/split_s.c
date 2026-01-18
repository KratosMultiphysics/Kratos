/* ============================================================================
**  This file is part of the mmg software package for the tetrahedra
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
 * \file mmgs/split_s.c
 * \brief Functions to create new points.
 * \author Charles Dapogny (UPMC
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux
 * \author Pascal Frey (UPMC
 * \author Algiane Froehly (Inria/UBordeaux
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "libmmgs_private.h"
#include "mmgsexterns_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of element to split.
 * \param i index of edge to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \param k index of element to split.
 * \return 0 if split leads to invalid element, else 1.
 *
 * Simulate the splitting of element \a k along edge \a i. Check that the new
 * triangles are not empty (otherwise we can create a 0 surface triangle).
 *
 */
int MMGS_split1_sim(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int i,MMG5_int *vx) {
  MMG5_pTria      pt,pt0;
  double          n[3],nref[3],vnew,vold;
  int             is;

  pt  = &mesh->tria[k];
  MMG5_nonUnitNorPts(mesh, pt->v[0], pt->v[1],pt->v[2],nref);

  vold = nref[0]*nref[0] + nref[1]*nref[1] + nref[2]*nref[2];
  if ( vold < MMG5_EPSOK ) return 0;

  pt0 = &mesh->tria[0];

  /* Test volume of the two created triangles */
  memcpy(pt0,pt,sizeof(MMG5_Tria));

  is         = MMG5_iprv2[i];
  pt0->v[is] = vx[i];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  pt0->v[is] = pt->v[is];
  is         = MMG5_inxt2[i];
  pt0->v[is] = vx[i];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of element to split.
 * \param i index of edge to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \return 1 if success, 0 if fail.
 *
 * Split element \a k along edge \a i.
 *
 */
int MMGS_split1(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,int i,MMG5_int *vx) {
  MMG5_pTria      pt,pt1;
  MMG5_pPoint     ppt;
  MMG5_int        iel;
  int8_t          i1,i2;

  iel = MMGS_newElt(mesh);
  if ( !iel ) {
    MMGS_TRIA_REALLOC(mesh,iel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       fprintf(stderr,"  Exit program.\n");
                       return 0);
  }

  pt  = &mesh->tria[k];
  pt->flag = 0;
  pt1 = &mesh->tria[iel];
  pt1 = memcpy(pt1,pt,sizeof(MMG5_Tria));

  i1 = MMG5_inxt2[i];
  i2 = MMG5_inxt2[i1];


  if ( pt->edg[i] > 0 ) {
    ppt = &mesh->point[vx[i]];
    ppt->ref = pt->edg[i];
  }

  pt->v[i2]   = pt1->v[i1]   = vx[i];
  pt->tag[i1] = pt1->tag[i2] = MG_NOTAG;
  pt->edg[i1] = pt1->edg[i2] = 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of the starting triangle.
 * \param i local index of the edge to split in \a k.
 * \param ip index of the point that we try to create.
 *
 * \return 0 if final position is invalid or if computation of bezier patch
 * fails, 1 if all checks are ok.
 *
 * Simulate the creation of the point \a ip, to be inserted at an
 * edge. Check that the new triangles are not empty (otherwise we can create a 0
 * surface triangle).
 *
 * \remark Don't work for non-manifold edge.
 */
int MMGS_simbulgept(MMG5_pMesh mesh,MMG5_pSol met, MMG5_int k,int i,MMG5_int ip) {
  MMG5_pTria     pt,pt0;
  MMG5_pPoint    ppt,ppt0;
  double         cal;
  MMG5_int       kadja;
  int            iadja,is;
  static int     mmgErr0 = 0, mmgErr1 = 0;

  pt0  = &mesh->tria[0];
  ppt0 = &mesh->point[0];
  ppt  = &mesh->point[ip];

  /* MMG5_calelt function needs the normal(s) at point in aniso mode so we can't
   * copy only the point coordinates */
  memcpy(ppt0 ,ppt ,sizeof(MMG5_Point));
  ppt0->tag = mesh->point[ip].tag;

  memcpy(&met->m[0],&met->m[met->size*ip], met->size*sizeof(double));

  /* Copy tria to split in tria 0 for simu purpose */
  pt = &mesh->tria[k];
  memcpy(pt0,pt,sizeof(MMG5_Tria));
  is         = MMG5_iprv2[i];
  pt0->v[is] = 0;

  /* For now, if ip is a ridge point, only 1 normal has been computed: allocate
   * a new xpoint and update the second normal */
  MMG5_int *adja = &mesh->adja[3*(k-1)+1];
  kadja = adja[i] / 3;
  iadja = adja[i] % 3;

  /* update normal n2 if need be */
  int8_t compute_n2 = kadja && (pt0->tag[i] & MG_GEO) ;
  if ( compute_n2 ) {
    MMG5_Bezier b;
    int ier = MMG5_bezierCP(mesh,&mesh->tria[kadja],&b,1);
    if ( !ier ) {
      if( !mmgErr0 ) {
        mmgErr0 = 1;
        fprintf(stderr,"\n  ## Warning: %s: function MMG5_bezierCP return 0.\n",
                __func__);
      }
      assert(0);
    }
    double uv[2],o[3],no[3],to[3];
    uv[0] = 0.5;
    uv[1] = 0.5;
    if ( iadja == 1 )       uv[0] = 0.0;
    else if ( iadja == 2 )  uv[1] = 0.0;

    ier = MMGS_bezierInt(&b,uv,o,no,to);
    if ( !ier ) {
      if( !mmgErr1 ) {
        mmgErr1 = 1;
        fprintf(stderr,"  ## Warning: %s: function MMGS_bezierInt return 0.\n",
                __func__);
      }
      assert(0);
    }
    MMG5_int nxp = mesh->xp + 1;
    if ( nxp > mesh->xpmax ) {
      MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                        "larger xpoint table",
                        return 0);
    }
    ppt0->xp = nxp;
    MMG5_pxPoint go  = &mesh->xpoint[ppt0->xp];
    memcpy(go->n2,no,3*sizeof(double));
    assert ( ppt->xp && "missing xpoint at ridge point" );
    MMG5_pxPoint pxp = &mesh->xpoint[ppt->xp];
    memcpy(go->n1,pxp->n1,3*sizeof(double));
  }

  // Check the validity of the two triangles created from k.
  cal        = MMG5_calelt(mesh,met,pt0);
  if ( cal < MMG5_EPSOK )  return 0;

  pt0->v[is] = pt->v[is];
  is         = MMG5_inxt2[i];
  pt0->v[is] = 0;
  cal        = MMG5_calelt(mesh,met,pt0);
  if ( cal < MMG5_EPSOK )  return 0;

  // Check the validity of the two triangles created from the triangle adjacent
  // to k by edge i.
  if ( !kadja ) return 1;

  pt = &mesh->tria[kadja];
  memcpy(pt0,pt,sizeof(MMG5_Tria));
  is         = MMG5_iprv2[iadja];
  pt0->v[is] = 0;
  cal        = MMG5_calelt(mesh,met,pt0);
  if ( cal < MMG5_EPSOK )  return 0;

  pt0->v[is] = pt->v[is];
  is         = MMG5_inxt2[iadja];
  pt0->v[is] = 0;
  cal        = MMG5_calelt(mesh,met,pt0);
  if ( cal < MMG5_EPSOK )  return 0;

  /* If point is ridge: copy computed second normal into xpoint */
  if ( compute_n2 ) {
    MMG5_pxPoint go  = &mesh->xpoint[ppt0->xp];
    MMG5_pxPoint pxp = &mesh->xpoint[ppt->xp];
    memcpy(pxp->n2,go->n2,3*sizeof(double));
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param k index of element to split.
 * \param i index of edge to split.
 * \param ip index of the new point.
 * \return 0 if lack of memory, 1 otherwise.
 *
 * Split element \a k along edge \a i, inserting point \a ip and updating
 * the adjacency relations.
 *
 * \remark do not call this function in non-manifold case
 */
int split1b(MMG5_pMesh mesh,MMG5_int k,int8_t i,MMG5_int ip) {
  MMG5_pTria     pt,pt1;
  MMG5_pPoint    ppt;
  MMG5_int       *adja,iel,jel,kel,mel;
  int8_t         i1,i2,j,j1,j2,m;

  iel = MMGS_newElt(mesh);
  if ( !iel )  {
    MMGS_TRIA_REALLOC(mesh,iel,mesh->gap,
                       MMG5_INCREASE_MEM_MESSAGE();
                       return 0);
  }
  pt = &mesh->tria[k];
  pt->flag = 0;
  pt->base = mesh->base;

  pt1 = &mesh->tria[iel];
  memcpy(pt1,pt,sizeof(MMG5_Tria));
  memcpy(&mesh->adja[3*(iel-1)+1],&mesh->adja[3*(k-1)+1],3*sizeof(MMG5_int));

  ppt = &mesh->point[ip];
  if ( pt->edg[i] )  ppt->ref = pt->edg[i];
  if ( pt->tag[i] )  ppt->tag = pt->tag[i];

  adja = &mesh->adja[3*(k-1)+1];
  jel = adja[i] / 3;
  j   = adja[i] % 3;

  /* update two triangles */
  i1  = MMG5_inxt2[i];
  i2  = MMG5_iprv2[i];
  pt->v[i2]   = ip;
  pt->tag[i1] = MG_NOTAG;
  pt->edg[i1] = 0;
  pt1->v[i1]   = ip;
  pt1->tag[i2] = MG_NOTAG;
  pt1->edg[i2] = 0;

  /* update adjacency relations */
  mel = adja[i1] / 3;
  m   = adja[i1] % 3;
  mesh->adja[3*(k-1)+1+i1]   = 3*iel+i2;
  mesh->adja[3*(iel-1)+1+i2] = 3*k+i1;
  mesh->adja[3*(iel-1)+1+i1] = 3*mel+m;
  if(mel)
    mesh->adja[3*(mel-1)+1+m]  = 3*iel+i1;

  if ( jel ) {
    kel = MMGS_newElt(mesh);
    if ( !kel )  {
      MMGS_TRIA_REALLOC(mesh,kel,mesh->gap,
                         MMG5_INCREASE_MEM_MESSAGE();
                         if ( !MMGS_delElt(mesh,iel) )  return 0;
                         return 0);
    }
    pt  = &mesh->tria[jel];
    pt1 = &mesh->tria[kel];
    pt->flag = 0;
    pt->base = mesh->base;
    memcpy(pt1,pt,sizeof(MMG5_Tria));
    memcpy(&mesh->adja[3*(kel-1)+1],&mesh->adja[3*(jel-1)+1],3*sizeof(MMG5_int));

    j1 = MMG5_inxt2[j];
    j2 = MMG5_iprv2[j];
    pt->v[j1]    = ip;
    pt->tag[j2]  = MG_NOTAG;
    pt->edg[j2]  = 0;
    pt1->v[j2]   = ip;
    pt1->tag[j1] = MG_NOTAG;
    pt1->edg[j1] = 0;

    /* update adjacency */
    adja = &mesh->adja[3*(jel-1)+1];
    mel  = adja[j2] / 3;
    m    = adja[j2] % 3;
    mesh->adja[3*(jel-1)+1+j2] = 3*kel+j1;
    mesh->adja[3*(kel-1)+1+j1] = 3*jel+j2;
    mesh->adja[3*(kel-1)+1+j2] = 3*mel+m;
    if(mel)
      mesh->adja[3*(mel-1)+1+m]  = 3*kel+j2;

    mesh->adja[3*(iel-1)+1+i]  = 3*kel+j;
    mesh->adja[3*(kel-1)+1+j]  = 3*iel+i;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \return 0 if split leads to invalid element, else 1.
 *
 * Simulate the splitting of element \a k along the 2 edges \a i1 and \a i2.
 * Check that the new triangles are not empty (otherwise we can create a 0
 * surface triangle).
 *
 */
int MMG5_split2_sim(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,MMG5_int *vx) {
  MMG5_pTria    pt,pt0;
  double        n[3],nref[3],vold,vnew;
  int           i1,i2,i;

  pt  = &mesh->tria[k];
  MMG5_nonUnitNorPts(mesh, pt->v[0], pt->v[1],pt->v[2],nref);

  vold = nref[0]*nref[0] + nref[1]*nref[1] + nref[2]*nref[2];
  if ( vold < MMG5_EPSOK ) return 0;

  pt0 = &mesh->tria[0];

  memcpy(pt0,pt,sizeof(MMG5_Tria));

  i = 0;
  if ( !vx[0] )  i = 1;
  else if ( !vx[1] )  i = 2;
  i1 = MMG5_inxt2[i];
  i2 = MMG5_inxt2[i1];

  /* Check the quality of the 3 new triangles */
  /* Tri 1 */
  pt0->v[i2] = vx[i];
  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  /* Tri 2 */
  pt0->v[i1] = vx[i];
  pt0->v[i2] = vx[i1];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  /* Tri 3 */
  pt0->v[i2] = pt->v[i2];
  pt0->v[i1] = vx[i];
  pt0->v[i]  = vx[i1];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \return 1 if success, 0 if fail.
 *
 * Split element \a k along the 2 edges \a i1 and \a i2.
 *
 */
int MMGS_split2(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,MMG5_int *vx) {
  MMG5_pTria    pt,pt1,pt2;
  MMG5_pPoint   p3,p4;
  MMG5_int      iel,jel;
  int8_t        i,i1,i2;

  /* create 2 elements */
  iel = MMGS_newElt(mesh);
  if ( !iel ) {
    MMGS_TRIA_REALLOC(mesh,iel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       fprintf(stderr,"  Exit program.\n");
                       return 0);
  }
  jel = MMGS_newElt(mesh);
  if ( !jel ) {
    MMGS_TRIA_REALLOC(mesh,jel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       fprintf(stderr,"  Exit program.\n");
                       return 0);
  }

  pt  = &mesh->tria[k];
  pt->flag = 0;
  pt1 = &mesh->tria[iel];
  pt2 = &mesh->tria[jel];
  pt1 = memcpy(pt1,pt,sizeof(MMG5_Tria));
  pt2 = memcpy(pt2,pt,sizeof(MMG5_Tria));

  i = 0;
  if ( !vx[0] )  i = 1;
  else if ( !vx[1] )  i = 2;
  i1 = MMG5_inxt2[i];
  i2 = MMG5_inxt2[i1];

  p3 = &mesh->point[vx[i]];
  p4 = &mesh->point[vx[i1]];

  /* update refs */
  if ( pt->edg[i] > 0 )   p3->ref = pt->edg[i];
  if ( pt->edg[i1] > 0 )  p4->ref = pt->edg[i1];

  pt->v[i1] = pt1->v[i2] = pt2->v[i1] = vx[i];
  pt->v[i2] = pt2->v[i]  = vx[i1];

  pt->tag[i] = pt->tag[i2] = pt1->tag[i1] = pt2->tag[i2] = MG_NOTAG;
  pt->edg[i] = pt->edg[i2] = pt1->edg[i1] = pt2->edg[i2] = 0;

  /* alternate configs */
  /* pt->v[i2]  = pt1->v[i]  = pt2->v[i] = vx[i1]; */
  /* pt1->v[i2] = pt2->v[i1] = vx[i]; */

  /* pt->tag[i] = pt1->tag[i1] = pt1->tag[i2] = pt2->tag[i2] = MG_NOTAG; */
  /* pt->edg[i] = pt1->edg[i1] = pt1->edg[i2] = pt2->edg[i2] = 0; */

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \return 0 if split leads to invalid element, else 1.
 *
 * Simulate the splitting of element \a k along the 3 edges. Check that the new
 * triangles are not empty (otherwise we can create a 0 surface triangle).
 *
 */
int MMGS_split3_sim(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,MMG5_int *vx) {
  MMG5_pTria    pt,pt0;
  double        n[3],nref[3],vnew,vold;

  pt   = &mesh->tria[k];
  MMG5_nonUnitNorPts(mesh, pt->v[0], pt->v[1],pt->v[2],nref);

  vold = nref[0]*nref[0] + nref[1]*nref[1] + nref[2]*nref[2];
  if ( vold < MMG5_EPSOK ) return 0;

  pt0  = &mesh->tria[0];

  memcpy(pt0,pt,sizeof(MMG5_Tria));

  /* Check the 4 new triangles */
  /* Tri 1 */
  pt0->v[1]  = vx[2];
  pt0->v[2]  = vx[1];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  /* Tri 2 */
  pt0->v[1]  = pt->v[1];
  pt0->v[0]  = vx[2];
  pt0->v[2]  = vx[0];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  /* Tri 3 */
  pt0->v[2]  = pt->v[2];
  pt0->v[0]  = vx[1];
  pt0->v[1]  = vx[0];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }


  /* Tri 4 */
  pt0->v[0]  = vx[2];
  pt0->v[1]  = vx[0];
  pt0->v[2]  = vx[1];

  MMG5_nonUnitNorPts(mesh, pt0->v[0], pt0->v[1],pt0->v[2],n);

  vnew = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
  if ( vnew < MMG5_EPSOK )  return 0;

  /* Check if we create a tri with wrong orientation */
  if ( nref[0]*n[0]+nref[1]*n[1]+nref[2]*n[2] < 0 ) {
    return 0;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param met pointer to the metric structure.
 * \param k index of element to split.
 * \param vx \f$vx[i]\f$ is the index of the point to add on the edge \a i.
 * \return 1 if success, 0 if fail.
 *
 * Split element \a k along the 3 edges
 *
 */
int MMGS_split3(MMG5_pMesh mesh,MMG5_pSol met,MMG5_int k,MMG5_int *vx) {
  MMG5_pTria    pt,pt1,pt2,pt3;
  MMG5_pPoint   p3,p4,p5;
  MMG5_int      iel,jel,kel;

  /* create 3 elements */
  iel = MMGS_newElt(mesh);
  if ( !iel ) {
    MMGS_TRIA_REALLOC(mesh,iel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       fprintf(stderr,"  Exit program.\n");
                       return 0);
  }
  jel = MMGS_newElt(mesh);
  if ( !jel ) {
    MMGS_TRIA_REALLOC(mesh,jel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       fprintf(stderr,"  Exit program.\n");
                       return 0);
  }
  kel = MMGS_newElt(mesh);
  if ( !kel ) {
    MMGS_TRIA_REALLOC(mesh,kel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       fprintf(stderr,"  Exit program.\n");
                       return 0);
  }

  pt  = &mesh->tria[k];
  pt->flag = 0;
  pt1 = &mesh->tria[iel];
  pt2 = &mesh->tria[jel];
  pt3 = &mesh->tria[kel];
  pt1 = memcpy(pt1,pt,sizeof(MMG5_Tria));
  pt2 = memcpy(pt2,pt,sizeof(MMG5_Tria));
  pt3 = memcpy(pt3,pt,sizeof(MMG5_Tria));

  p3 = &mesh->point[vx[0]];
  p4 = &mesh->point[vx[1]];
  p5 = &mesh->point[vx[2]];

  /* update refs */
  if ( pt->edg[0] > 0 )  p3->ref = pt->edg[0];
  if ( pt->edg[1] > 0 )  p4->ref = pt->edg[1];
  if ( pt->edg[2] > 0 )  p5->ref = pt->edg[2];

  /* update topo */
  pt->v[1]  = pt1->v[0] = pt3->v[0] = vx[2];
  pt->v[2]  = pt2->v[0] = pt3->v[2] = vx[1];
  pt1->v[2] = pt2->v[1] = pt3->v[1] = vx[0];

  pt->tag[0]  = pt1->tag[1] = pt2->tag[2] = MG_NOTAG;
  pt->edg[0]  = pt1->edg[1] = pt2->edg[2] = 0;

  pt3->tag[0] = pt3->tag[1] = pt3->tag[2] = MG_NOTAG;
  pt3->edg[0] = pt3->edg[1] = pt3->edg[2] = 0;

  return 1;
}
