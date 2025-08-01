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
#ifndef _WIN32
#include "mmg/common/git_log_mmg.h"
#endif

#ifndef MMGVERSION_H
#define MMGVERSION_H

#define MMG_VERSION_RELEASE "5.7.0"
#define MMG_VERSION_MAJOR 5
#define MMG_VERSION_MINOR 7
#define MMG_VERSION_PATCH 0
#define MMG_RELEASE_DATE "Dec. 13, 2022"

#define MMG_COPYRIGHT   "Copyright (c) Bdx INP/CNRS/Inria/UPMC, 2004-"

#define MMG_VERSION_EQ(MAJOR,MINOR) \
((MMG_VERSION_MAJOR == (MAJOR)) && (MMG_VERSION_MINOR == (MINOR)))

#define MMG_VERSION_ MMG_VERSION_EQ

#define MMG_VERSION_LT(MAJOR,MINOR)                                  \
  (MMG_VERSION_MAJOR < (MAJOR) || (MMG_VERSION_MAJOR == (MAJOR) &&   \
                                   (MMG_VERSION_MINOR < (MINOR) )))

#define MMG_VERSION_LE(MAJOR,MINOR) \
  (MMG_VERSION_LT(MAJOR,MINOR) || MMG_VERSION_EQ(MAJOR,MINOR))

#define MMG_VERSION_GT(MAJOR,MINOR) (0 == MMG_VERSION_LE(MAJOR,MINOR))

#define MMG_VERSION_GE(MAJOR,MINOR) (0 == MMG_VERSION_LT(MAJOR,MINOR))

#endif
