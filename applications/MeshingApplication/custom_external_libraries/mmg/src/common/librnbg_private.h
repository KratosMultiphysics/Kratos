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
 * \brief header file for the librnbg.c librnbg_s.c librnbg_3d.c files
 * \author Cedric Lachat  (Inria/UBordeaux)
 * \version 5
 * \date 2013
 * \copyright GNU Lesser General Public License.
 */

#ifdef USE_SCOTCH

#ifndef __RENUM__
#define __RENUM__

#include <scotch.h>

#define HASHPRIME 37

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define SCOTCH_5 (!strcmp(TOSTRING(SCOTCH_VERSION),"5.0") ||            \
                  !strcmp(TOSTRING(SCOTCH_VERSION),"5.1") || !strcmp(TOSTRING(SCOTCH_VERSION),"5"))

#define SCOTCH_6 !strcmp(TOSTRING(SCOTCH_VERSION),"6")

#define SCOTCH_7 !strcmp(TOSTRING(SCOTCH_VERSION),"7")

#define CHECK_SCOTCH(t,m,e) if(0!=t){perror(m);return e;}

int    _SCOTCHintSort2asc1(SCOTCH_Num * sortPartTb, MMG5_int vertNbr);
int    MMG5_kPartBoxCompute(SCOTCH_Graph*, MMG5_int, MMG5_int, SCOTCH_Num*,MMG5_pMesh);
void   MMG5_swapNod(MMG5_pMesh,MMG5_pPoint, double*, MMG5_pSol,MMG5_int*, MMG5_int, MMG5_int, int);

#endif /* __RENUM__ */
#endif
