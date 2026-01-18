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
 */

#include "libmmgs_private.h"
#include "mmgexterns_private.h"

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
static int MMGS_cuttri(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;
  double       c[3],v0,v1,s;
  MMG5_int     vx[3],k,np,refint,refext;
  MMG5_int     ip0,ip1,ns,nt,ier,nb;
  int8_t       i;

  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Evaluate the number of intersected edges by the 0 level set */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      /* If only surface edges are discretized, skip non boundary entities: as
       * mmgs doesn't add MG_BDY tags, we check if an edge is bdy from the
       * MG_REF tag */
      if ( mesh->info.isosurf && !(pt->tag[i] & MG_REF) ) continue;

      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];

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

  /* Create the intersection points between the edges in the mesh and the 0
   * level set */
  if ( !MMG5_hashNew(mesh,&hash,nb,2*nb) ) return 0;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      if ( mesh->info.isosurf && !(pt->tag[i] & MG_REF) ) {
        continue;
      }

      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];

      np = MMG5_hashGet(&hash,ip0,ip1);
      if ( np ) continue;

      /* Look either at the triangle ref or at the boundary one */
      MMG5_int ref;
      if ( mesh->info.isosurf ) {
        ref = pt->edg[i];
      }
      else {
        ref = pt->ref;
      }

      /* If user asks to keep input refs, ignore multi-mat mode */
      if ( mesh->info.iso!=2 ) {
        if ( !MMG5_isSplit(mesh,ref,&refint,&refext) ) continue;
      }

      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];

      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
      else if ( MG_SMSGN(v0,v1) )  continue;
      else if ( !p0->flag || !p1->flag )  continue;

      /* Intersection point between edge p0p1 and the 0 level set */
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
      sol->m[np] = 0.0;

      /* If user provide a metric, interpolate it at the new point */
      if ( met && met->m ) {
        if ( met->size > 1 ) {
          ier = MMGS_intmet33_ani(mesh,met,k,i,np,s);
        }
        else {
          ier = intmet_iso(mesh,met,k,i,np,s);
        }
        if ( ier <= 0 ) {
          // Unable to compute the metric
          fprintf(stderr,"\n  ## Error: %s: unable to"
                  " interpolate the metric during the level-set"
                  " discretization\n",__func__);
          return 0;
        }
      }

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
    memset(vx,0,3*sizeof(MMG5_int));
    for (i=0; i<3; i++) {
      vx[i] = MMG5_hashGet(&hash,pt->v[MMG5_inxt2[i]],pt->v[MMG5_iprv2[i]]);
      if ( vx[i] ) {
        MG_SET(pt->flag,i);
      }
    }
    switch (pt->flag) {
    case 1: /* 1 edge split */
      ier = MMGS_split1(mesh,met,k,0,vx);
      ns++;
      break;

    case 2: /* 1 edge split */
      ier = MMGS_split1(mesh,met,k,1,vx);
      ns++;
      break;

    case 4: /* 1 edge split */
      ier = MMGS_split1(mesh,met,k,2,vx);
      ns++;
      break;

    case 3: case 5: case 6: /* 2 edges split */
      ier = MMGS_split2(mesh,met,k,vx);
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

  /* reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;

}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set
 * \param met pointer to a metric (optionnal)
 *
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int MMGS_mmgs2(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  char str[16]="";
  MMG5_int k;

  assert ( (mesh->info.iso || mesh->info.isosurf) && "level-set discretization mode not specified" );

  assert ( (!(mesh->info.iso && mesh->info.isosurf)) && "unable to determine level-set discretization mode" );

  /* Set function pointers */
  if ( mesh->info.isosurf ) {
    strcat(str,"(BOUNDARY PART)");

    MMG5_snpval   = MMG5_snpval_lssurf;
    MMG5_resetRef = MMG5_resetRef_lssurf;
    MMG5_setref   = MMG5_setref_lssurf;
  }
  else {
    MMG5_snpval   = MMG5_snpval_ls;
    MMG5_resetRef = MMG5_resetRef_ls;
    MMG5_setref   = MMG5_setref_ls;
  }

  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION %s\n",str);
  }

  /* Work only with the 0 level set */
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /* Transfer the boundary edge references to the triangles */
  if ( !MMGS_assignEdge(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* Snap values of level set function if need be, then discretize it */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Error: %s: hashing problem (1). Exit program.\n",
            __func__);
    return 0;
  }
  /* identify connexity and reorient triangles (needed for manifoldness checks) */
  if ( !MMGS_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    return 0;
  }

  /* Snap values of the level set function which are very close to 0 to 0 exactly */
  if ( !MMG5_snpval(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    return 0;
  }

  if ( mesh->info.iso ) {
    /* Removal of small parasitic components */
    if ( mesh->info.rmc > 0. && !MMG5_rmc(mesh,sol) ) {
      fprintf(stderr,"\n  ## Error in removing small parasitic components. Exit program.\n");
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

  MMG5_DEL_MEM(mesh,mesh->adja);

  if ( mesh->info.iso != 2 ) {
    /* Reset the mesh->info.isoref field everywhere it appears */
    if ( !MMG5_resetRef(mesh) ) {
      fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
      return 0;
    }
  }

  if ( !MMGS_cuttri(mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
    return 0;
  }

  if ( !MMG5_setref(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Creation of adjacency relations in the mesh */
  if ( !MMGS_hashTria(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  if ( mesh->info.iso ) {
    /* Check that the resulting mesh is manifold */
    if ( !MMG5_chkmanimesh(mesh) ) {
      fprintf(stderr,"\n  ## No manifold resulting situation. Exit program.\n");
      return 0;
    }
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);
  sol->np = 0;

  MMG5_DEL_MEM( mesh,mesh->info.mat );

  return 1;
}
