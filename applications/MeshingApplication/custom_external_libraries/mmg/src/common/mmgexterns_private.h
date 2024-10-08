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

#ifndef MMGEXTERNS_H
#define MMGEXTERNS_H

#include "mmgcommon_private.h"

#ifndef MMG_EXTERN
#define MMG_EXTERN extern
#define MMG_ASSIGN_NULL
#endif

FUNCTION_POINTER ( int  (*MMG5_chkmsh)(MMG5_pMesh,int,MMG5_int) );
FUNCTION_POINTER ( int  (*MMG5_bezierCP)(MMG5_pMesh ,MMG5_Tria *,MMG5_pBezier ,int8_t ) );
FUNCTION_POINTER ( double (*MMG5_lenSurfEdg)(MMG5_pMesh mesh,MMG5_pSol sol ,MMG5_int ,MMG5_int, int8_t ) );
FUNCTION_POINTER ( MMG5_int  (*MMG5_indElt)(MMG5_pMesh mesh,MMG5_int kel) );
FUNCTION_POINTER ( MMG5_int  (*MMG5_indPt)(MMG5_pMesh mesh,MMG5_int kp) );
FUNCTION_POINTER ( MMG5_int  (*MMG5_grad2met_ani)(MMG5_pMesh,MMG5_pSol,MMG5_pTria,MMG5_int,MMG5_int) );
FUNCTION_POINTER ( int (*MMG5_grad2metreq_ani)(MMG5_pMesh,MMG5_pSol,MMG5_pTria,MMG5_int,MMG5_int) );
FUNCTION_POINTER ( int (*MMG5_compute_meanMetricAtMarkedPoints)( MMG5_pMesh,MMG5_pSol) );
FUNCTION_POINTER ( int (*MMG5_solTruncature_ani)(MMG5_pMesh mesh, MMG5_pSol met) );
FUNCTION_POINTER ( int (*MMG5_resetRef)(MMG5_pMesh) );
FUNCTION_POINTER ( int (*MMG5_setref)(MMG5_pMesh,MMG5_pSol) );
FUNCTION_POINTER ( int (*MMG5_snpval)(MMG5_pMesh,MMG5_pSol));

#ifdef USE_SCOTCH
FUNCTION_POINTER ( int  (*MMG5_renumbering)(int,MMG5_pMesh,MMG5_pSol,MMG5_pSol,MMG5_int*) );
#endif

#undef MMG_EXTERN
#undef MMG_ASSIGN_NULL

#endif
