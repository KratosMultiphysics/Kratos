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
 * Example of use of the mmg3d library (mmg3d application)
 *
 * \author Charles Dapogny (LJLL, UPMC)
 * \author Cécile Dobrzynski (Inria / IMB, Université de Bordeaux)
 * \author Pascal Frey (LJLL, UPMC)
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <ctype.h>

/** Include the mmg3d library hader file */
// if the header file is in the "include" directory
// #include "libmmg3d.h"
// if the header file is in "include/mmg/mmg3d"
#include "mmg/mmg3d/libmmg3d.h"

mytime    ctim[TIMEMAX];

#define RETURN_AND_FREE(mesh,met,val)do    \
    {                                                 \
      MMG5_Free_all(mesh,met);                        \
      return(val);                                    \
    }while(0)

static inline void endcod() {
  char    stim[32];

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3DLIB: ELAPSED TIME  %s\n",stim);
}

int main(int argc,char *argv[]) {
  MMG5_pMesh      mesh;
  MMG5_pSol       met,disp;
  int             ier,dummy,typSol,np;
  char            stim[32];

  atexit(endcod);

  tminit(ctim,TIMEMAX);
  chrono(ON,&ctim[0]);

  /* assign default values */
  mesh = NULL;
  met  = NULL;
  disp = NULL;

  MMG5_Init_mesh(&mesh,&met,&disp);

  /* reset default values for file names */
  MMG5_Free_names(mesh,met,disp);

  /* command line */
  if ( !MMG5_parsar(argc,argv,mesh,met) )  return(1);

  /* load data */
  fprintf(stdout,"\n  -- INPUT DATA\n");
  chrono(ON,&ctim[1]);
  /* read mesh file */
  if ( MMG5_loadMesh(mesh)<1 ) {
    MMG5_Free_all(mesh,met,disp );
    return(MMG5_STRONGFAILURE);
  }
  /* Set default metric size */
  if ( !MMG5_Set_solSize(mesh,met,MMG5_Vertex,0,MMG5_Scalar) ) {
    MMG5_Free_all(mesh,met,disp);
    return(MMG5_STRONGFAILURE);
  }

  /* read displacement if any */
  if ( MMG5_Get_iparameter(mesh, MMG5_IPARAM_lag) > -1 ) {
    if ( !MMG5_Set_solSize(mesh,disp,MMG5_Vertex,0,MMG5_Vector) ) {
      MMG5_Free_all(mesh,met,disp);
      return(MMG5_STRONGFAILURE);
    }

    if ( !MMG5_Set_inputSolName(mesh,disp,met->namein) ) {
      MMG5_Free_all(mesh,met,disp);
      return(MMG5_STRONGFAILURE);
    }
    ier = MMG5_loadMet(mesh,disp);
    if ( ier == 0 ) {
      fprintf(stdout,"  ## ERROR: NO DISPLACEMENT FOUND.\n");
      MMG5_Free_all(mesh,met,disp);
      return(MMG5_STRONGFAILURE);
    }
    else if ( ier == -1 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      MMG5_Free_all(mesh,met,disp);
      return(MMG5_STRONGFAILURE);
    }
  }
  /* read metric if any */
  else {
    ier = MMG5_loadMet(mesh,met);
    if ( ier == -1 ) {
      fprintf(stdout,"  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      MMG5_Free_all(mesh,met,disp);
      return(MMG5_STRONGFAILURE);
    }
    else {
      MMG5_Get_solSize( mesh, met, &dummy, &dummy, &typSol);
      if ( (typSol != MMG5_Scalar) && (typSol != MMG5_Tensor) ) {
        fprintf(stdout,"  ## ERROR: WRONG DATA TYPE.\n");
        MMG5_Free_all(mesh,met,disp);
        return(MMG5_STRONGFAILURE);
      }
    }
    if ( MMG5_Get_iparameter(mesh, MMG5_IPARAM_iso) && !ier ) {
      fprintf(stdout,"  ## ERROR: NO ISOVALUE DATA.\n");
      MMG5_Free_all(mesh,met,disp);
      return(MMG5_STRONGFAILURE);
    }
    if ( !MMG5_parsop(mesh,met) )
      return(MMG5_LOWFAILURE);
  }

  chrono(OFF,&ctim[1]);
  printim(ctim[1].gdif,stim);
  fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);

  if ( MMG3D_Get_iparameter(mesh, MMG3D_IPARAM_lag) > -1 ) {
    ier = MMG3D_mmg3dmov(mesh,met,disp);
  }
  else {
    ier = MMG5_mmg3dlib(mesh,met);
  }

  if ( ier != MMG5_STRONGFAILURE ) {
    chrono(ON,&ctim[1]);
    if ( mesh->info.imprim )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);
    if ( !MMG5_saveMesh(mesh) )         {
      MMG5_Free_all(mesh,met,disp );
      return(MMG5_STRONGFAILURE);
    }
    if ( !MMG5_saveMet(mesh,met) )     {
      MMG5_Free_all(mesh,met,disp);
      return(MMG5_STRONGFAILURE);
    }
    chrono(OFF,&ctim[1]);
    if ( mesh->info.imprim )
      fprintf(stdout,"  -- WRITING COMPLETED\n");
  }

  /* free mem */
  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  fprintf(stdout,"\n   MMG3D: ELAPSED TIME  %s\n",stim);
  MMG5_Free_all(mesh,met,disp);
  return(ier);
}
