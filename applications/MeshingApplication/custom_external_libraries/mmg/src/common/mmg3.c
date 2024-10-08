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
 * \brief common functions for lagrangian meshing.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "mmgcommon_private.h"

/**
 * \param mesh pointer toward the mesh structure
 * \param disp pointer toward the displacement field
 * \param lastt 0 if a movement is possible, pointer toward the last tested fraction otherwise
 * \param shortmax maximal parameter t (MMG2D_SHORTMAX or MMG3D_SHORTMAX)
 * \param chkmovmesh function that has to be called to check motion validity
 *
 * \return the largest fraction t allowing a valid motion
 *
 * Generic function to compute the the largest fraction t that makes the motion
 * along disp valid.
 *
 */
short MMG5_dikmov ( MMG5_pMesh mesh,MMG5_pSol disp,short *lastt,short shortmax,
                    MMG5_int chkmovmesh(MMG5_pMesh,MMG5_pSol,short,MMG5_int*) ) {
  int     it,maxit;
  short   t,tmin,tmax;
  int8_t  ier;

  maxit = 200;
  it = 0;

  tmin = 0;
  tmax = shortmax;

  *lastt = 0;

  /* If full displacement can be achieved */
  if ( !chkmovmesh(mesh,disp,tmax,NULL) )
    return tmax;

  /* Else, find the largest displacement by dichotomy */
  while( tmin != tmax && it < maxit ) {
    t = (tmin+tmax)/2;

    /* Case that tmax = tmin +1 : check move with tmax */
    if ( t == tmin ) {
      ier = chkmovmesh(mesh,disp,tmax,NULL);
      if ( !ier ) {
        return tmax;
      }
      else {
        if ( tmin==0 ) {
          *lastt = tmax;
        }
        return tmin;
      }
    }

    /* General case: check move with t */
    ier = chkmovmesh(mesh,disp,t,NULL);
    if ( !ier ) {
      tmin = t;
    }
    else
      tmax = t;

    it++;
  }

  if ( tmin==0 ) {
    *lastt=t;
  }

  return tmin;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param disp pointer toward the displacement field
 *
 * \return 1 if success, 0 if fail.
 *
 * For debugging purposes: save displacement field.
 *
 */
int MMG5_saveDisp(MMG5_pMesh mesh,MMG5_pSol disp) {
  FILE        *out;
  MMG5_int    k;
  char        data[256],*ptr;

  strcpy(data,disp->namein);
  ptr = strstr(data,".sol");
  *ptr = '\0';
  strcat(data,".o.disp.sol");

  out = fopen(data,"w");

  fprintf(out,"MeshVersionFormatted 1\n\nDimension\n%d\n\n",disp->dim);
  fprintf(out,"SolAtVertices\n%"MMG5_PRId"\n 1 2\n",disp->np);

  /* Print solutions */
  for(k=1; k<= disp->np; k++) {
    int i;
    for ( i=0; i<mesh->dim; ++i ) {
      fprintf(out," %f",disp->m[mesh->dim*k+i]);
    }
    fprintf(out,"\n");
  }

  fprintf(out,"\nEnd");
  fclose(out);

  return 1;
}
