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
 * Compare the preservation of parallel triangles of too files when the input
 * triangle index is store in the triangle references (so each triangle of one
 * file has a unique reference and among the 2 files, two maching triangles have
 * to have the same ref).
 *
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/** Include the mmg3d library hader file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "libmmg3d_private.h"
#include "mmg/mmg3d/libmmg3d.h"

#define  TRIA_PARBDY(tag)    ((tag[0] & MG_PARBDY) && (tag[1] & MG_PARBDY) && (tag[2] & MG_PARBDY) )

int main(int argc,char *argv[]) {
  MMG5_pMesh      mesh1,mesh2;
  MMG5_pSol       sol1,sol2;
  MMG5_HGeom      hash;
  char            *file1, *file2;

  fprintf(stdout,"  -- TEST MMG3DLIB \n");

  if ( argc != 3 ) {
    printf(" Usage: %s filein fileout \n",argv[0]);
    return(1);
  }

  /* Name and path of the mesh file */
  file1 = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( file1 == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(file1,argv[1]);

  file2 = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( file2 == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(file2,argv[2]);

  /** Read meshes */
  mesh1 = NULL;
  mesh2 = NULL;
  sol1  = NULL;
  sol2  = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh1,
                  MMG5_ARG_ppMet,&sol1,
                  MMG5_ARG_end);

  if ( MMG3D_loadMesh(mesh1,file1) != 1 )  {
    fprintf(stderr,"Error: %s: %d: File not found %s\n.",__func__,__LINE__,file1);
    exit(EXIT_FAILURE);
  }


  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mesh2,
                  MMG5_ARG_ppMet,&sol2,
                  MMG5_ARG_end);

  if ( MMG3D_loadMesh(mesh2,file2) != 1 )  {
    fprintf(stderr,"Error: %s: %d: File not found %s\n.",__func__,__LINE__,file2);
    exit(EXIT_FAILURE);
  }

  /** Hash parallel triangle using their reference as key */
  MMG5_int k;
  MMG5_int maxref = 0;
  MMG5_int npar1   = 0;
  for ( k=1; k<=mesh1->nt; ++k ) {
    MMG5_pTria pt = &mesh1->tria[k];
    if ( (!MG_EOK(pt)) || !(TRIA_PARBDY(pt->tag)) ) {
      continue;
    }
    ++npar1;
    maxref = (maxref < pt->ref) ? pt->ref : maxref;
  }

  /* Check data consistency with mesh2 */
  MMG5_int npar2   = 0;
  for ( k=1; k<=mesh2->nt; ++k ) {
    MMG5_pTria pt = &mesh2->tria[k];
    if ( (!MG_EOK(pt)) || !(TRIA_PARBDY(pt->tag)) ) {
      continue;
    }
    ++npar2;
    if ( pt->ref > maxref ) {
      fprintf(stderr,"Error: %s: %d:"
              " maximal reference of parallel boundaries mimatch:"
              " max ref in mesh1 is %"MMG5_PRId" while ref %"MMG5_PRId" "
              "founded in mesh2.\n",
              __func__,__LINE__,maxref,pt->ref);
      exit(EXIT_FAILURE);
    }
  }

  if ( npar2 != npar1 ) {
    fprintf(stderr,"Error: %s: %d:"
            " number of parallel boundaries mimatch: %"MMG5_PRId" versus %"MMG5_PRId".\n",
            __func__,__LINE__,npar1,npar2);
    exit(EXIT_FAILURE);
  }

  hash.siz = mesh1->nt + 1;
  hash.max = mesh1->nt + 1;
  hash.nxt = 0;
  MMG5_SAFE_CALLOC(hash.geom,(hash.max+1),MMG5_hgeom,
                   perror("  ## Memory problem: calloc");exit(EXIT_FAILURE));

  for ( k=1; k<=mesh1->nt; ++k ) {
    MMG5_pTria pt = &mesh1->tria[k];
    if ( (!MG_EOK(pt)) || !TRIA_PARBDY(pt->tag) ) {
      continue;
    }
    hash.geom[pt->ref].a = k;
  }

  /** Travel meshes and compare parallel triangles with same references: they
   * should have vertices with same references and same coordinates (but we can
   * have different indices) */
  for ( k=1; k<=mesh2->nt; ++k ) {
    MMG5_pTria pt2 = &mesh2->tria[k];
    if ( (!MG_EOK(pt2)) || !(TRIA_PARBDY(pt2->tag)) ) {
      continue;
    }

    MMG5_int k1 = hash.geom[pt2->ref].a;
    if ( !k1 ) {
      fprintf(stderr,"Error: %s: %d:"
            " face of parallel ref %"MMG5_PRId" not found in hashtable.\n",
              __func__,__LINE__,pt2->ref);
      exit(EXIT_FAILURE);
    }
    MMG5_pTria pt1 = &mesh1->tria[k1];
    if ( (!MG_EOK(pt1)) || !(TRIA_PARBDY(pt1->tag)) ) {
      fprintf(stderr,"Error: %s: %d:"
            " parallel tria %"MMG5_PRId" in mesh2 (ref %"MMG5_PRId") is not"
              " valid in mesh1 (tria %"MMG5_PRId").\n",
              __func__,__LINE__,k,pt2->ref,k1);
      exit(EXIT_FAILURE);
    }

    /* Coordinates comparison */
    MMG5_pPoint ppt1,ppt2;

    /* Find matching vertices */
    ppt1 = &mesh1->point[pt1->v[0]];
    int shift;
    for ( shift=0; shift<3; ++shift) {
      ppt2 = &mesh2->point[pt2->v[shift]];
      if ( ppt2->ref == ppt1->ref ) {
        break;
      }
    }
    if ( shift==3 ) {
      fprintf(stderr,"Error: %s: %d:"
              " ref of vertices mismatch.\n",
              __func__,__LINE__);
      fprintf(stderr," Vertices of triangle %"MMG5_PRId" in mesh1 have refs"
              " %"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId"\n",
              k1,mesh1->point[pt1->v[0]].ref,mesh1->point[pt1->v[1]].ref,
              mesh1->point[pt1->v[2]].ref);
      fprintf(stderr," Vertices of triangle %"MMG5_PRId" in mesh2 have refs"
              " %"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId"\n",
              k,mesh2->point[pt2->v[0]].ref,mesh2->point[pt2->v[1]].ref,
              mesh2->point[pt2->v[2]].ref);
      exit(EXIT_FAILURE);
    }

    int i;
    for (i=0; i<3; ++i) {
      ppt1 = &mesh1->point[pt1->v[i]];

      int i2 = (i+shift)%3;
      ppt2 = &mesh2->point[pt2->v[i2]];

      if ( ppt1->ref != ppt2->ref ) {
        fprintf(stderr,"Error: %s: %d:"
                " ref of vertices %d and %d mismatch.\n",
                __func__,__LINE__,i,i2);
        fprintf(stderr," Vertices of triangle %"MMG5_PRId" in mesh1 have "
                "refs %"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId"\n",
                k1,mesh1->point[pt1->v[0]].ref,mesh1->point[pt1->v[1]].ref,
                mesh1->point[pt1->v[2]].ref);
        fprintf(stderr," Vertices of triangle %"MMG5_PRId" in mesh2 have "
                "refs %"MMG5_PRId" %"MMG5_PRId" %"MMG5_PRId"\n",
                k,mesh2->point[pt2->v[0]].ref,mesh2->point[pt2->v[1]].ref,
                mesh2->point[pt2->v[2]].ref);
        exit(EXIT_FAILURE);
      }

      int j;
      for (j=0; j<3; ++j) {
        if ( fabs(ppt1->c[j]-ppt2->c[j]) > 1e-5 ) {
          fprintf(stderr,"Error: %s: %d:"
                  " Elts %"MMG5_PRId" and %"MMG5_PRId": coor %d of"
                  " vertices %d and %d (ref %"MMG5_PRId")"
                  " differs: %15lf versus %15lf.\n",__func__,__LINE__,
                  k1,k,j,i,i2,ppt1->ref,ppt1->c[j],ppt2->c[j]);
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  MMG5_SAFE_FREE(hash.geom);

  /** 3) Free the MMG3D5 structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh1,
                 MMG5_ARG_ppMet,&sol1,
                 MMG5_ARG_end);
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh2,
                 MMG5_ARG_ppMet,&sol2,
                 MMG5_ARG_end);

  free(file1);
  file1 = NULL;

  free(file2);
  file2 = NULL;

  return(0);
}
