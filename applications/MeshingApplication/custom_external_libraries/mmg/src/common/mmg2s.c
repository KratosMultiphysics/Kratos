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
#include "mmgcommon_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set
 *
 * \return 1 if success, 0 if fail
 *
 * Snap values of sol very close to 0 to 0 exactly in the case of surface ls splitting
 * (to avoid very small triangles in cutting)
 */
int MMG5_snpval_lssurf(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTria       pt;
  MMG5_pPoint      p0,p1;
  double           v0,v1,*tmp;
  int              k,ns,nc,ip0,ip1;
  int8_t           i,i0,i1;

  /* Allocate memory for tmp */
  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                printf("  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double,return 0);

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]) < MMG5_EPS ) {
      tmp[k] =  sol->m[k];
      p0->flag = 1;
      sol->m[k] = 0.0;
      ns++;
    }
  }

  /* Unsnap values that have been put to 0, entailing a component reduced to a point */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      if ( !(pt->tag[i] & MG_REF) ) continue;
      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      ip0 = pt->v[i0];
      ip1 = pt->v[i1];

      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      p0 = &mesh->point[pt->v[i0]];
      p1 = &mesh->point[pt->v[i1]];

      if ( fabs(v0) < MMG5_EPS && fabs(v1) < MMG5_EPS ) {
        if ( p0->flag ) {
          if ( tmp[ip0] < 0.0 )
            sol->m[ip0] = -100.0*MMG5_EPS;
          else
            sol->m[ip0] = 100.0*MMG5_EPS;
          nc++;
          p0->flag = 0;
        }
        else if ( p1->flag ) {
          if ( tmp[ip1] < 0.0 )
            sol->m[ip1] = -100.0*MMG5_EPS;
          else
            sol->m[ip1] = 100.0*MMG5_EPS;
          nc++;
          p1->flag = 0;
        }
      }
    }
  }

  MMG5_DEL_MEM ( mesh, tmp );

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+nc > 0 )
    fprintf(stdout,"     %8d points snapped, %d corrected\n",ns,nc);

  return 1;
}

/**
 * \param mesh pointer to the mesh
 *
 * Reset mesh->info.isoref vertex references to 0.
 *
 */
int MMG5_resetRef_lssurf(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0,p1;
  MMG5_int        k,ref;
  int8_t          i,i0,i1;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      if ( !(pt->tag[i] & MG_REF) ) continue;

      if( !MMG5_getStartRef(mesh,pt->edg[i],&ref) ) return 0;
      pt->edg[i] = ref;

      i0 = MMG5_inxt2[i];
      i1 = MMG5_inxt2[i0];

      p0 = &mesh->point[pt->v[i0]];
      p1 = &mesh->point[pt->v[i1]];

      if ( p0->ref == mesh->info.isoref ) p0->ref = 0;
      if ( p1->ref == mesh->info.isoref ) p1->ref = 0;
    }
  }

  return 1;
}

/* Set references to the new triangles */
int MMG5_setref_lssurf(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria     pt;
  double         v,v1;
  MMG5_int       k,ip1,ref,refint,refext;
  int8_t         ier,i,i1,j,nmn,npl,nz;

  /* Travel all surface edges (via triangles) */
  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      if ( !(pt->tag[i] & MG_REF) ) continue;
      ref = pt->edg[i];
      nmn = npl = nz = 0;

      i1 = i;
      for (j=0; j<2; j++) {
        i1  = MMG5_inxt2[i1];
        ip1 = pt->v[i1];
        v1  = sol->m[ip1];

        if ( v1 > 0.0 )
          npl++;
        else if ( v1 < 0.0 )
          nmn++;
        else
          nz++;
      }

      assert(nz < 2);
      ier = MMG5_isSplit(mesh,ref,&refint,&refext);
      if ( npl ) {
        if (ier ) {
          assert ( !nmn );
          pt->edg[i] = refext;
        }
      }
      else {
        if ( ier ) {
          assert ( !npl );
          pt->edg[i] = refint;
        }
      }
    }
  }

  /* Set MG_ISO ref at vertices on the ls */
  for (k=1; k<=mesh->np; k++) {
    v = sol->m[k];
    if ( v == 0.0 ) {
      mesh->point[k].ref = MG_ISO;
    }
  }

  return 1;
}
