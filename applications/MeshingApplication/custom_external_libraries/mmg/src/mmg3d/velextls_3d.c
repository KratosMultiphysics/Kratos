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
 * \file mmg3d/velextls_3d.c
 * \brief Tools for interfacing mmg with LS, for extension of the displacement field.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */
#ifdef USE_ELAS

#include "libmmg3d_private.h"
#include "ls_calls.h"

#define MMG5_DEGTOL    0.75
#define _LS_LAMBDA      10.0e5
#define _LS_MU          8.2e5

/**
 * \param mesh pointer to the mesh
 * \param disp pointer to the displacement
 * \param lsst pointer to the elastic structure (mesh + sol + info)
 * \param npfin pointer to the final number of points in the packed mesh
 * for the elasticity call
 *
 * \return invperm array of the permutation (vertices) from the submesh (given
 * to the elastic library) toward the global one.
 *
 * Create submesh for solving the linear elasticity velocity extension problem.
 * invperm stores the permutation [ new pt nb -> old pt nb ] (for unpacking purposes)
 * Fill npf = number of vertices in the packed mesh.
 *
 */
MMG5_int* MMG5_packLS(MMG5_pMesh mesh,MMG5_pSol disp,LSst *lsst,MMG5_int *npfin) {
  MMG5_pTetra    pt,pt1;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    p0;
  double         u[3];
  int            n,nlay,ilist,ilisto,ilistck;
  MMG5_int       k,ip,npf,ntf,iel,jel,*perm,*invperm,*adja,*list,vper[4];
  int            refdirh,refdirnh;
  int8_t         i,j,jface;

  nlay = 20;
  refdirh = 0;
  refdirnh = 1;
  npf = 0;
  ntf = 0;
  u[0] = u[1] = u[2] = 0.0;
  *npfin = 0;

  MMG5_ADD_MEM(mesh,(mesh->ne+1)*sizeof(MMG5_int),"element list",return NULL);
  MMG5_SAFE_CALLOC(list,mesh->ne+1,MMG5_int,return NULL);

  MMG5_ADD_MEM(mesh,(mesh->np+1)*sizeof(MMG5_int),"point permutation",return NULL);
  MMG5_SAFE_CALLOC(perm,mesh->np+1,MMG5_int,return NULL);

  ilist = ilisto = ilistck = 0;

  for (k=1; k<=mesh->ne; k++)
    mesh->tetra[k].mark = 0;      // A faire fusionner avec celui de mmg3d3.c

  /* Step 1: pile all the tetras containing a triangle with ref DISPREF */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for(i=0; i<4; i++) {
      if ( (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == MMG5_DISPREF) ) {
        ilist++;
        list[ilist] = k;
        MG_SET(pt->mark,0);

        for(j=0; j<4; j++) {
          ip = pt->v[j];
          if ( !perm[ip] ) {
            npf++;
            perm[ip] = npf;
          }
        }
        break;
      }
    }
  }
  if ( !npf ) {
    fprintf(stderr,
            "\n  ## Error: %s: no triangle with reference %d in the mesh.\n"
            "              Nothing to move.\n",__func__,MMG5_DISPREF);
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return NULL;
  }

  /* Step 2: create a layer around these tetras */
  for(n=0; n<nlay; n++) {
    ilistck = ilisto;
    ilisto = ilist;

    for(k=ilistck+1; k<=ilisto; k++) {
      iel = list[k];
      adja = &mesh->adja[4*(iel-1)+1];

      for(i=0; i<4; i++) {
        jel = adja[i] / 4;
        if ( !jel ) continue;
        pt1 = &mesh->tetra[jel];
        if ( MG_EOK(pt1) && (!MG_GET(pt1->mark,0) ) ) {
          ilist++;
          assert( ilist <= mesh->ne );
          MG_SET(pt1->mark,0);
          list[ilist] = jel;

          for(j=0; j<4; j++) {
            ip = pt1->v[j];
            if ( !perm[ip] ) {
              npf++;
              perm[ip] = npf;
            }
          }
        }
      }
    }
  }

  /* Creation of the inverse permutation table */
  MMG5_ADD_MEM ( mesh,(npf+1)*sizeof(MMG5_int),"permutation table",
                  MMG5_DEL_MEM ( mesh,list );
                  MMG5_DEL_MEM ( mesh,perm );
                  return NULL );
  MMG5_SAFE_CALLOC ( invperm,(npf+1),MMG5_int,
                      MMG5_DEL_MEM ( mesh,list );
                      MMG5_DEL_MEM ( mesh,perm );
                      return NULL );

  /* Step 3: count of the surface triangles in the new mesh
     Code for pt->mark : if MG_GET(pt->mark,0) = in the list
     if !MG_GET(pt->mark,0) = not in the list
     if MG_GET(pt->mark, i+1) : face i has a boundary triangle which has already been counted */
  for(k=1; k<=ilist; k++) {
    iel = list[k];
    pt = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];
    if (pt->xt) pxt = &mesh->xtetra[pt->xt];

    for(i=0; i<4; i++) {
      jel = adja[i] / 4;
      jface = adja[i] % 4;

      /* Face i carries a non homogeneous Dirichlet BC */
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == MMG5_DISPREF) ) {
        /* If this triangle has not been taken into account */
        if ( MG_GET(pt->mark,i+1) ) continue;

        ntf++;
        if ( !jel ) continue;
        pt1 = &mesh->tetra[jel];
        MG_SET(pt1->mark,jface+1);
      }
      /* iel has no neighbour through face i within the list */
      else if ( !jel || ( !mesh->tetra[jel].mark ) ) ntf++;
    }
  }

  /* Step 4: creation of the mesh for elasticity */
  if ( !LS_mesh(lsst,npf,0,ntf,ilist) ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn LS_mesh. Exiting.\n",
            __func__);
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    MMG5_DEL_MEM ( mesh,invperm );
    return NULL;
  }

  /* Set verbosity and debug info */
  LS_setPar(lsst, (mesh->info.imprim > 0), 0);

  /* Step 5: fill the LS mesh */
  /* Add vertices */
  for(k=1; k<=mesh->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;

    p0 = &mesh->point[k];
    invperm[ip] = k;
    if ( !LS_addVer(lsst,ip,p0->c,p0->ref) ) {
      fprintf(stderr,"\n  ## Error: %s: problem in fn LS_addVer. Exiting.\n",
              __func__);
      MMG5_DEL_MEM ( mesh,list );
      MMG5_DEL_MEM ( mesh,perm );
      MMG5_DEL_MEM ( mesh,invperm );
      return NULL;
    }
  }

  /* Add tetrahedra */
  for(k=1; k<=ilist; k++) {
    iel = list[k];
    pt = &mesh->tetra[iel];

    for(i=0; i<4; i++)
      vper[i] = perm[pt->v[i]];

    if (!LS_addTet(lsst,(int)k,(int*)vper,0) ) {
      fprintf(stderr,"\n  ## Error: %s: problem in fn LS_addTet. Exiting.\n",
              __func__);
      MMG5_DEL_MEM ( mesh,list );
      MMG5_DEL_MEM ( mesh,perm );
      MMG5_DEL_MEM ( mesh,invperm );
      return NULL;
    }
  }

  /* Add surface triangles */
  ntf = 0;
  for(k=1; k<=ilist; k++) {
    iel = list[k];
    pt = &mesh->tetra[iel];
    adja = &mesh->adja[4*(iel-1)+1];
    if (pt->xt) pxt = &mesh->xtetra[pt->xt];

    for(i=0; i<4; i++) {
      jel = adja[i] / 4;

      /* Face i carries a non homogeneous Dirichlet BC */
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == MMG5_DISPREF) ) {
        /* If this triangle has not been taken into account */
        if ( MG_GET(pt->mark,i+1) ) continue;
        ntf++;
        for (j=0; j<3; j++)
          vper[j] = perm[pt->v[MMG5_idir[i][j]]];

        if ( !LS_addTri(lsst,(int)ntf,(int*)vper,refdirnh) ) {
          fprintf(stderr,"\n  ## Error: %s: problem in fn LS_addTri. Exiting.\n",
                  __func__);
          MMG5_DEL_MEM ( mesh,list );
          MMG5_DEL_MEM ( mesh,perm );
          MMG5_DEL_MEM ( mesh,invperm );
          return NULL;
        }
      }
      /* iel has no neighbour through face i within list */
      else if ( !jel || (!mesh->tetra[jel].mark) ) {
        ntf++;
        for (j=0; j<3; j++)
          vper[j] = perm[pt->v[MMG5_idir[i][j]]];

        if ( !LS_addTri(lsst,(int)ntf,(int*)vper,refdirh) ) {
          fprintf(stderr,"\n  ## Error: %s: problem in fn LS_addTri. Exiting.\n",
                  __func__);
          MMG5_DEL_MEM ( mesh,list );
          MMG5_DEL_MEM ( mesh,perm );
          MMG5_DEL_MEM ( mesh,invperm );
          return NULL;
        }
      }
    }
  }

  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && (ilist+npf+ntf > 0) )
    printf("Number of packed tetra %d, points %" MMG5_PRId
           ", triangles %" MMG5_PRId "\n",ilist,npf,ntf);

  /* Add boundary conditions */
  if ( !LS_setBC(lsst,Dirichlet,refdirnh,'f',LS_tri,NULL) ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn LS_set BC. Exiting.\n",
            __func__);
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    MMG5_DEL_MEM ( mesh,invperm );
    return NULL;
  }

  if ( !LS_setBC(lsst,Dirichlet,refdirh,'v',LS_tri,u) ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn LS_set BC. Exiting.\n",
            __func__);
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    MMG5_DEL_MEM ( mesh,invperm );
    return NULL;
  }

  /* Add materials */
  if ( !LS_setLame(lsst,0,_LS_LAMBDA,_LS_MU) ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn LS_setLame. Exiting.\n",
            __func__);
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    MMG5_DEL_MEM ( mesh,invperm );
    return NULL;
  }

  /* Transfer displacement */
  if ( !LS_newSol(lsst) ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn LS_newSol. Exiting.\n",
            __func__);
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    MMG5_DEL_MEM ( mesh,invperm );
    return NULL;
  }

  for(k=1; k<=mesh->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;

    if ( !LS_addSol(lsst,ip,&disp->m[3*k]) ) {
      fprintf(stderr,"\n  ## Error: %s: problem in fn LS_addSol. Exiting.\n",
              __func__);
      MMG5_DEL_MEM ( mesh,list );
      MMG5_DEL_MEM ( mesh,perm );
      MMG5_DEL_MEM ( mesh,invperm );
      return NULL;
    }
  }

  *npfin = npf;
  MMG5_DEL_MEM ( mesh,list );
  MMG5_DEL_MEM ( mesh,perm );

  return invperm;
}

/**
 * \param mesh pointer to the mesh
 * \param disp pointer to the displacement
 * \param lsst pointer to the elastic structure (mesh + sol + info)
 * \param npf pointer to the number of points in the submesh
 * \param invperm array of the permutation from the submesh toward
 * the global one
 *
 * \return 1 if success, 0 otherwise
 *
 * Transfer solution from the submesh to the global mesh
 *
 */
int MMG5_unpackLS(MMG5_pMesh mesh,MMG5_pSol disp,LSst *lsst,MMG5_int npf,MMG5_int *invperm) {
  double     *u;
  MMG5_int   ip,k;
  int8_t     i;

  u = LS_getSol(lsst);

  for(k=1; k<=mesh->np; k++) {
    for(i=0; i<3; i++)
      disp->m[3*k+i] = 0.0;
  }

  for(k=1; k<=npf; k++) {
    ip = invperm[k];

    for(i=0; i<3; i++)
      disp->m[3*ip+i] = u[3*(k-1)+i];
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh.
 * \param disp pointer to the displacement.
 *
 * \return 0 if fail, 1 if success.
 *
 * Extension of the displacement at the nodes of triangles tagged MMG5_DISPREF
 *
 */
int MMG5_velextLS(MMG5_pMesh mesh,MMG5_pSol disp) {
  LSst        *lsst;
  MMG5_int    npf,*invperm;

  /* LibElas is not compatible with int64: Check for int32 overflow */
  if ( mesh->np > INT_MAX || mesh->ne > INT_MAX || sizeof(MMG5_int) == 8 ) {
    fprintf(stderr,"\n  ## Error: %s: impossible to call elasticity library"
            " with int64 integers.\n",__func__);
    return 0;
  }

  /* Creation of the data structure for the submesh */
  lsst    = LS_init(mesh->dim,mesh->ver,P1,1);
  invperm = MMG5_packLS(mesh,disp,lsst,&npf);

  if ( !npf ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn MMG5_packLS. Exiting.\n",
            __func__);
    return 0;
  }

  /* Resolution of the elasticity system on the submesh */
  if ( !LS_elastic(lsst) ) {
    fprintf(stderr,"\n  ## Error: %s: Problem in fn elasti1. Exiting.\n",
            __func__);
    return 0;
  }

  /* Update of the displacement */
  if ( !MMG5_unpackLS(mesh,disp,lsst,npf,invperm) ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn MMG5_unpackLS. Exiting.\n",
            __func__);
    return 0;
  }

  /* Free memory */
  MMG5_DEL_MEM ( mesh, invperm );

  if ( !LS_stop(lsst) ) {
    fprintf(stderr,"\n  ## Error: %s: problem in fn LS_stop. Exiting.\n",
            __func__);
    return 0;
  }

  return 1;
}
#else
/**
 *
 * Hack to avoid to have an empty translation unit (forbidden by ISO C)
 *
 */
#ifdef _WIN32
void MMG3D_unused_function(void) {
  return;
}
#else
void __attribute__((unused)) MMG3D_unused_function(void) {
  return;
}
#endif

#endif
