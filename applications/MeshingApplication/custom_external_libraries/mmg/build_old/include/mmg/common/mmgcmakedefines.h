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

#ifndef MMGCMAKEDEFINE_H
#define MMGCMAKEDEFINE_H

/* inttypes.h is needed to handle prints of MMG5_int using PRId32 and PRId64 macros */
#include <inttypes.h>

#define MMG_POSIX
#define MMG_GNU

#define MMG5_int int32_t /*!< Integer type for C */

#define MMG5_PRId PRId32 /*!< Printing format for MMG5_int type */

#define MMG5_SWAPBIN(x) MMG5_swapbin(x) /*!< MMG5_swapbin function for MMG5_int */

#define MMG5_abs(x) abs(x) /*!< Abs function for MMG5_int */

/* #undef USE_POINTMAP */
#define MMG_DYN_LIB 0

#endif
