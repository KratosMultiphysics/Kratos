/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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

#include "mmg3d.h"
#include "ls_calls.h"
#define _MMG5_DEGTOL    0.75
#define _MMG5_DISPREF   0
#define _LS_LAMBDA      10.0e5
#define _LS_MU          8.2e5

/** Create submesh for solving the linear elasticity velocity extension problem.
 invperm stores the permutation [ new pt nb -> old pt nb ] (for unpacking purposes)
 Return: npf = number of vertices in the packed mesh.
 */
int* _MMG5_packLS(MMG5_pMesh mesh,MMG5_pSol disp,LSst *lsst,int *npfin) {
  MMG5_pTetra    pt,pt1;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    p0;
  double         u[3];
  int            k,n,ip,iel,jel,nlay,npf,ntf,ilist,ilisto,ilistck,vper[4],*list,*perm,*invperm,*adja;
  int            refdirh,refdirnh;
  char           i,j,jface;
  
  nlay = 20;
  refdirh = 0;
  refdirnh = 1;
  npf = 0;
  ntf = 0;
  u[0] = u[1] = u[2] = 0.0;
  list = (int*)calloc(mesh->ne+1,sizeof(int));
  perm = (int*)calloc(mesh->np+1,sizeof(int));
  ilist = ilisto = ilistck = 0;
  
  for (k=1; k<=mesh->ne; k++)
    mesh->tetra[k].mark = 0;      // A faire fusionner avec celui de mmg3d3.c
  
  /* Step 1: pile all the tetras containing a triangle with ref DISPREF */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) || !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];
    
    for(i=0; i<4; i++) {
      if ( (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == _MMG5_DISPREF) ) {
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
  
  /* Step 2: create a layer around these tetras */
  for(n=0; n<nlay; n++) {
    ilistck = ilisto;
    ilisto = ilist;
    
    for(k=ilistck+1; k<=ilisto; k++) {
      iel = list[k];
      pt = &mesh->tetra[iel];
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
  invperm = (int*)calloc(npf+1,sizeof(int));
  
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
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == _MMG5_DISPREF) ) {
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
    fprintf(stderr,"  ## Problem in fn LS_mesh. Exiting.\n");
    return(0);
  }

  /* Set verbosity and debug info */
  if ( !mesh->info.imprim )
    LS_setPar(lsst,0,0);
  else
    LS_setPar(lsst,1,0);
  
  /* Step 5: fill the LS mesh */
  /* Add vertices */
  for(k=1; k<=mesh->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;
    
    p0 = &mesh->point[k];
    invperm[ip] = k;
    if ( !LS_addVer(lsst,ip,p0->c,p0->ref) ) {
      fprintf(stderr,"  ## Problem in fn LS_addVer. Exiting.\n");
      return(0);
    }
  }
  
  /* Add tetrahedra */
  for(k=1; k<=ilist; k++) {
    iel = list[k];
    pt = &mesh->tetra[iel];
    
    for(i=0; i<4; i++)
      vper[i] = perm[pt->v[i]];
    
    if (!LS_addTet(lsst,k,vper,0) ) {
      fprintf(stderr,"  ## Problem in fn LS_addTet. Exiting.\n");
      return(0);
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
      jface = adja[i] % 4;
      
      /* Face i carries a non homogeneous Dirichlet BC */
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) && (pxt->ref[i] == _MMG5_DISPREF) ) {
        /* If this triangle has not been taken into account */
        if ( MG_GET(pt->mark,i+1) ) continue;
        ntf++;
        for (j=0; j<3; j++)
          vper[j] = perm[pt->v[_MMG5_idir[i][j]]];
        
        if ( !LS_addTri(lsst,ntf,vper,refdirnh) ) {
          fprintf(stderr,"  ## Problem in fn LS_addTri. Exiting.\n");
          return(0);
        }
      }
      /* iel has no neighbour through face i within list */
      else if ( !jel || (!mesh->tetra[jel].mark) ) {
        ntf++;
        for (j=0; j<3; j++)
          vper[j] = perm[pt->v[_MMG5_idir[i][j]]];
        
        if ( !LS_addTri(lsst,ntf,vper,refdirh) ) {
          fprintf(stderr,"  ## Problem in fn LS_addTri. Exiting.\n");
          return(0);
        }
      }
    }
  }

  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && (ilist+npf+ntf > 0) )
       printf("Number of packed tetra %d, points %d, triangles %d\n",ilist,npf,ntf);
  
  /* Add boundary conditions */
  if ( !LS_setBC(lsst,Dirichlet,refdirnh,'f',LS_tri,NULL) ) {
    fprintf(stderr,"  ## Problem in fn LS_set BC. Exiting.\n");
    return(0);
  }
  
  if ( !LS_setBC(lsst,Dirichlet,refdirh,'v',LS_tri,u) ) {
    fprintf(stderr,"  ## Problem in fn LS_set BC. Exiting.\n");
    return(0);
  }
  
  /* Add materials */
  if ( !LS_setLame(lsst,0,_LS_LAMBDA,_LS_MU) ) {
    fprintf(stderr,"  ## Problem in fn LS_setLame. Exiting.\n");
    return(0);
  }
  
  /* Transfer displacement */
  if ( !LS_newSol(lsst) ) {
    fprintf(stderr,"  ## Problem in fn LS_newSol. Exiting.\n");
    return(0);
  }
  
  for(k=1; k<=mesh->np; k++) {
    ip = perm[k];
    if ( !ip ) continue;
    
    if ( !LS_addSol(lsst,ip,&disp->m[3*k]) ) {
      fprintf(stderr,"  ## Problem in fn LS_addSol. Exiting.\n");
      return(0);
    }
  }
  
  *npfin = npf;
  free(list);
  free(perm);
  return(invperm);
}

/** Transfer solution from the submesh to the global mesh */
int _MMG5_unpackLS(MMG5_pMesh mesh,MMG5_pSol disp,LSst *lsst,int npf,int *invperm) {
  double     *u;
  int        k,ip;
  char       i;
  
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

  return(1);
}

/** Extension of the displacement at the nodes of triangles tagged _MMG5_DISPREF */
int _MMG5_velextLS(MMG5_pMesh mesh,MMG5_pSol disp) {
  LSst        *lsst;
  int         npf,*invperm;
  
  /* Creation of the data structure for the submesh */
  lsst = LS_init(mesh->dim,mesh->ver,P1,1);
  invperm = _MMG5_packLS(mesh,disp,lsst,&npf);
  
  if ( !npf ) {
    fprintf(stderr,"  ## Problem in fn MMG5_packLS. Exiting.\n");
    return(0);
  }
  
  /* Resolution of the elasticity system on the submesh */
  if ( !LS_elastic(lsst) ) {
    fprintf(stderr,"  ## Problem in fn elasti1. Exiting.\n");
    return(0);
  }
  
  /* Update of the displacement */
  if ( !_MMG5_unpackLS(mesh,disp,lsst,npf,invperm) ) {
    fprintf(stderr,"  ## Problem in fn _MMG5_unpackLS. Exiting.\n");
    return(0);
  }
  
  /* Free memory */
  free(invperm);
  
  
  if ( !LS_stop(lsst) ) {
    fprintf(stderr,"  ## Problem in fn LS_stop. Exiting.\n");
    return(0);
  }
  
  return(1);
}

#endif
