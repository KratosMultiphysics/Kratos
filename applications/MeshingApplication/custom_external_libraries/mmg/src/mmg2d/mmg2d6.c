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
#include "libmmg2d_private.h"
#include "mmgexterns_private.h"
#include "mmg2dexterns_private.h"

/**
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the level-set
 * \param met pointer toward a metric (non-mandatory)
 *
 * \return 1 if success, 0 otherwise
 *
 * Effective discretization of the 0 level set encoded in sol in the mesh.
 * Only the boundary part of the domain is discretized if mesh->info.isosurf is 1.
 *
 */
int MMG2D_cuttri(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSol met){
  MMG5_pTria   pt;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;
  double       v0,v1,s,c[2];
  MMG5_int     k,ip0,ip1,nb,np,nt,ns,refint,refext,vx[3];
  int8_t       i,ier;

  /* Reset flag field for points */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Estimate the number of intersected edges or surface edges by the 0 level set */
  nb = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {

      /* If only surface edges are discretized, skip non boundary entities */
      if ( mesh->info.isosurf && !(pt->tag[i] & MG_BDY) ) continue;

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

  /* Create the intersection points between the edges or surface edges
   * in the mesh and the 0 level set */
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

      if ( !MMG5_isSplit(mesh,ref,&refint,&refext) ) continue;

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

      np = MMG2D_newPt(mesh,c,0);
      if ( !np ) {
       /* reallocation of point table */
        MMG2D_POINT_REALLOC(mesh,met,np,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            return 0;,
                            c,0);
      }
      sol->m[np] = 0.0;
      /* If there is a metric in the mesh, interpolate it at the new point */
      if ( met && met->m )
        MMG2D_intmet(mesh,met,k,i,np,s);

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
      ip0 = pt->v[MMG5_inxt2[i]];
      ip1 = pt->v[MMG5_iprv2[i]];

      vx[i] = MMG5_hashGet(&hash,ip0,ip1);

      if ( vx[i] ) MG_SET(pt->flag,i);
    }

    switch( pt->flag ) {
      /* 1 edge split -> 0-+ */
      case 1: case 2: case 4:
        ier = MMG2D_split1(mesh,met,k,vx);
        ns++;
        break;

      /* 2 edge split -> +-- or -++ */
      case 3: case 5: case 6:
        ier = MMG2D_split2(mesh,met,k,vx);
        ns++;
        break;

      default:
        assert(pt->flag==0);
        break;
    }
    if ( !ier ) return 0;
  }

  if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);

  MMG5_DEL_MEM(mesh,hash.item);
  return ns;

}

/* Main function of the -ls mode */
int MMG2D_mmg2d6(MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met) {
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

  if ( mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }

  /* Work only with the 0 level set */
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /* Transfer the boundary edge references to the triangles */
  if ( !MMG2D_assignEdge(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* Processing of bdy infos for isosurface boundary discretization */
  if ( mesh->info.isosurf ) {
    /* Creation of tria adjacency relations in the mesh */
    if ( !MMG2D_hashTria(mesh) ) {
       fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
      return 0;
    }

    /* Set tags to triangles from geometric configuration */
    if ( !MMG2D_setadj(mesh,0) ) {
      fprintf(stderr,"\n  ## Problem in function setadj. Exit program.\n");
      return 0;
    }
  }

  /* Snap values of the level set function which are very close to 0 to 0 exactly */
  if ( !MMG5_snpval(mesh,sol) ) {
    fprintf(stderr,"\n  ## Wrong input implicit function. Exit program.\n");
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

  /* No need to keep adjacencies from now on */
  MMG5_DEL_MEM(mesh,mesh->adja);

  /* Reset the mesh->info.isoref field everywhere it appears */
  if ( !MMG5_resetRef(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    return 0;
  }

  /* Effective splitting of the crossed triangles */
  if ( !MMG2D_cuttri(mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in cutting triangles. Exit program.\n");
    return 0;
  }

  /* Set references on the interior / exterior triangles*/
  if ( !MMG5_setref(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Creation of adjacency relations in the mesh */
  if ( !MMG2D_hashTria(mesh) ) {
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
