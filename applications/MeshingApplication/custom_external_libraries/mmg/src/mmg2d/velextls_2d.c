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
 * \file mmg2d/mmg2d9.c
 * \brief Velocity extension for Lagrangian meshing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#ifdef USE_ELAS

#define _LS_LAMBDA      10.0e5
#define _LS_MU          8.2e5

#include "libmmg2d_private.h"
#include "ls_calls.h"

/** Create submesh for solving the linear elasticity velocity extension problem.
    invperm stores the permutation [ new pt nb -> old pt nb ] (for unpacking purposes)
    Return: npf = number of vertices in the packed mesh.
*/
MMG5_int* MMG2D_packLS(MMG5_pMesh mesh,MMG5_pSol disp,LSst *lsst,MMG5_int *npfin) {
  MMG5_pTria      pt,pt1;
  MMG5_pPoint     p0;
  double          u[2];
  MMG5_int        k,iel,jel,n,npf,nef,ip,nlay,refdirh,refdirnh,ilist,ilisto,ilistck;
  MMG5_int        vper[3],*perm,*list,*adja,*invperm;
  int8_t          i,j,jedg;

  nlay       = 20;
  refdirh    = 0;
  refdirnh   = 1;
  npf        = 0;
  nef        = 0;
  u[0]       = u[1] = 0.0;
  ilist      = ilisto = ilistck = 0;
  MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_int),"element list",return NULL);
  MMG5_SAFE_CALLOC(list,mesh->nt+1,MMG5_int,return NULL);

  MMG5_ADD_MEM(mesh,(mesh->np+1)*sizeof(MMG5_int),"point permutation",return NULL);
  MMG5_SAFE_CALLOC(perm,mesh->np+1,MMG5_int,return NULL);


  /* Reset flag field at triangles */
  for(k=1; k<=mesh->nt; k++)
    mesh->tria[k].flag = 0;

  /* Step 1: pile up all the triangles with one edge with ref DISPREF, and get the corresponding points */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      if ( pt->edg[i] == MMG5_DISPREF ) {
        ilist++;
        list[ilist] = k;
        MG_SET(pt->flag,0);

        for (j=0; j<3; j++) {
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

  /* Step 2: Create a hull of nlay layers around these triangles */
  for (n=0; n<nlay; n++) {
    ilistck = ilisto;
    ilisto  = ilist;

    for (k=ilistck+1; k<=ilisto; k++) {
      iel   = list[k];
      adja  = &mesh->adja[3*(iel-1)+1];

      for (i=0; i<3; i++) {
        jel = adja[i] / 3;
        if ( !jel ) continue;
        pt1 = &mesh->tria[jel];

        if ( MG_EOK(pt1) && ( !MG_GET(pt1->flag,0) ) ) {
          ilist++;
          assert ( ilist <= mesh->nt );
          MG_SET(pt1->flag,0);
          list[ilist] = jel;

          for (j=0; j<3; j++) {
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
  if ( !npf ) {
    fprintf(stderr,
            "\n  ## Error: %s: no triangle with reference %d in the mesh.\n"
            "              Nothing to move.\n",__func__,MMG5_DISPREF);
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return NULL;
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
     Code for pt->flag : if pt->flag = 1 -> in the list
     if pt->flag = 0 -> not in the list
     if MG_GET(pt->mark, i+1) : edge i has already been counted */
  for (k=1; k<=ilist; k++) {
    iel  = list[k];
    pt   = &mesh->tria[iel];
    adja = &mesh->adja[3*(iel-1)+1];

    for (i=0; i<3; i++) {
      jel   = adja[i] / 3;
      jedg  = adja[i] % 3;

      /* Face i carries a non homogeneous Dirichlet BC */
      if ( pt->edg[i] == MMG5_DISPREF ) {

        /* If this triangle has not been taken into account */
        if ( MG_GET(pt->flag,i+1) ) continue;

        nef++;
        if ( !jel ) continue;
        pt1 = &mesh->tria[jel];
        MG_SET(pt1->flag,jedg+1);
      }
      /* iel has no neighbour through face i within the list */
      else if ( !jel || !MG_GET(mesh->tria[jel].flag,0) ) nef++;
    }
  }

  /* Step 4: creation of the mesh for elasticity */
  if ( !LS_mesh(lsst,npf,nef,ilist,0) ) {
    fprintf(stdout,"  ## Problem in function LS_mesh. Exiting.\n");
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return 0;
  }

  /* Set verbosity and debug info */
  LS_setPar(lsst, (mesh->info.imprim > 0),0);

  /* Step 5: fill the LS mesh */
  /* Add vertices, and fill in table invperm on the fly */
  for (k=1; k<=mesh->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;

    p0 = &mesh->point[k];
    invperm[ip] = k;

    if ( !LS_addVer(lsst,ip,p0->c,p0->ref) ) {
      fprintf(stdout,"  ## Problem in fn LS_addVer. Exiting.\n");
      MMG5_DEL_MEM ( mesh,list );
      MMG5_DEL_MEM ( mesh,perm );
      return 0;
    }
  }

  /* Add triangles */
  for (k=1; k<=ilist; k++) {
    iel = list[k];
    pt = &mesh->tria[iel];

    for (i=0; i<3; i++)
      vper[i] = perm[pt->v[i]];

    if (!LS_addTri(lsst,(int)k,(int*)vper,0) ) {
      fprintf(stdout,"  ## Problem in fn LS_addTet. Exiting.\n");
      MMG5_DEL_MEM ( mesh,list );
      MMG5_DEL_MEM ( mesh,perm );
      return 0;
    }
  }

  /* Add boundary edges */
  nef = 0;
  for (k=1; k<=ilist; k++) {
    iel    = list[k];
    pt     = &mesh->tria[iel];
    adja   = &mesh->adja[3*(iel-1)+1];

    for (i=0; i<3; i++) {
      jel  = adja[i] / 3;

      /* Case where face i carries a non homogeneous Dirichlet BC */
      if ( pt->edg[i] == MMG5_DISPREF ) {
        /* If this triangle has not been taken into account */
        if ( MG_GET(pt->flag,i+1) ) continue;
        nef ++;

        vper[0] = perm[pt->v[MMG5_inxt2[i]]];
        vper[1] = perm[pt->v[MMG5_iprv2[i]]];

        if ( !LS_addEdg(lsst,(int)nef,(int*)vper,refdirnh) ) {
          fprintf(stdout,"  ## Problem in fn LS_addEdg. Exiting.\n");
          MMG5_DEL_MEM ( mesh,list );
          MMG5_DEL_MEM ( mesh,perm );
          return 0;
        }
      }
      /* iel has no neighbour through face i within the list */
      else if ( !jel || !MG_GET(mesh->tria[jel].flag,0) ) {
        nef++;
        vper[0] = perm[pt->v[MMG5_inxt2[i]]];
        vper[1] = perm[pt->v[MMG5_iprv2[i]]];

        if ( !LS_addEdg(lsst,(int)nef,(int*)vper,refdirh) ) {
          fprintf(stdout,"  ## Problem in fn LS_addEdg. Exiting.\n");
          MMG5_DEL_MEM ( mesh,list );
          MMG5_DEL_MEM ( mesh,perm );
          return 0;
        }
      }
    }
  }

  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && (ilist+npf+nef > 0) )
    printf("Number of packed triangles %" MMG5_PRId ", points %" MMG5_PRId ", edges %" MMG5_PRId "\n",ilist,npf,nef);

  /* Add boundary conditions */
  /*if ( !LS_setBC(lsst,Dirichlet,refdirnh,'v',LS_edg,v) ) {
    fprintf(stdout,"  ## Problem in fn LS_set BC. Exiting.\n");
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return 0;
    }*/
  if ( !LS_setBC(lsst,Dirichlet,refdirnh,'f',LS_edg,NULL) ) {
    fprintf(stdout,"  ## Problem in fn LS_set BC. Exiting.\n");
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return 0;
  }

  if ( !LS_setBC(lsst,Dirichlet,refdirh,'v',LS_edg,u) ) {
    fprintf(stdout,"  ## Problem in fn LS_set BC. Exiting.\n");
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return 0;
  }

  /* Add materials */
  if ( !LS_setLame(lsst,0,_LS_LAMBDA,_LS_MU) ) {
    fprintf(stdout,"  ## Problem in fn LS_setLame. Exiting.\n");
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return 0;
  }

  /* Transfer displacement */
  if ( !LS_newSol(lsst) ) {
    fprintf(stdout,"  ## Problem in fn LS_CreaSol. Exiting.\n");
    MMG5_DEL_MEM ( mesh,list );
    MMG5_DEL_MEM ( mesh,perm );
    return 0;
  }

  for(k=1; k<=mesh->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;
    if ( !LS_addSol(lsst,ip,&disp->m[2*k]) ) {
      fprintf(stdout,"  ## Problem in fn LS_addSol. Exiting.\n");
      return 0;
    }
  }

  *npfin = npf;
  MMG5_DEL_MEM ( mesh,list );
  MMG5_DEL_MEM ( mesh,perm );
  return invperm;
}

/** Transfer solution from the submesh to the global mesh */
int MMG2D_unpackLS(MMG5_pMesh mesh,MMG5_pSol disp,LSst *lsst,MMG5_int npf,MMG5_int *invperm) {
  double      *u;
  MMG5_int    k,ip;
  int8_t      i;

  u = LS_getSol(lsst);

  for (k=1; k<=mesh->np; k++) {
    for (i=0; i<2; i++)
      disp->m[2*k+i] = 0.0;
  }

  for (k=1; k<=npf; k++) {
    ip = invperm[k];
    for (i=0; i<2; i++)
      disp->m[2*ip+i] = u[2*(k-1)+i];
  }

  return 1;
}

/** Extension of the displacement at the nodes of edges tagged MMG5_DISPREF */
int MMG2D_velextLS(MMG5_pMesh mesh,MMG5_pSol disp) {
  LSst       *lsst;
  MMG5_int   npf,*invperm;

  /* Creation of the data structure for storing the submesh */
  lsst = LS_init(mesh->dim,mesh->ver,P1,1);

  /* Creation of the submesh */
  invperm = MMG2D_packLS(mesh,disp,lsst,&npf);

  if ( !npf ) {
    fprintf(stdout,"  ## Problem in fn MMG2D_packLS. Exiting.\n");
    return 0;
  }

  /* Resolution of the elasticity system on the submesh */
  if ( !LS_elastic(lsst) ) {
    fprintf(stdout,"  ## Problem in function elasti1. Exiting.\n");
    return 0;
  }

  /* Update of the displacement */
  if ( !MMG2D_unpackLS(mesh,disp,lsst,npf,invperm) ) {
    fprintf(stdout,"  ## Problem in fn MMG2D_unpackLS. Exiting.\n");
    return 0;
  }

  /* Memory release */
  MMG5_DEL_MEM ( mesh,invperm );

  /* Release of the Lsst structure */
  if ( !LS_stop(lsst) ) {
    fprintf(stdout,"  ## Problem in fn LS_stop. Exiting.\n");
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
void MMG2D_unused_function(void) {
  return;
}
#else
void __attribute__((unused)) MMG2D_unused_function(void) {
  return;
}
#endif

#endif
