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

/*
 * WARNING: To keep the genheader program working, don't break line between an enum
 * name and the opening brace (it creates errors under windows)
 */

/* The following comment block defines the Doxygen group with label API, which
 * includes
 *   mmg3d/libmmg3d.h
 *   mmg2d/libmmg2.h
 *   mmgs/libmmg2d.h
 * and others
 * with the aim to generate a central index into the API documentation.
 */

/**
 * \defgroup API Application Programming Interface
 *
 * Mmg comes with three application programming interfaces (APIs), one
 * corresponding to each of the command-line programs. In total there are over
 * 360 API functions. You can find their declarations and descriptions in the
 * files listed above.
 *
 * Examples showing how to use these functions can be found under the "Related
 * Pages" tab.
 *
 */

/**
 * \file common/libmmg.h
 * \brief Wrapper for include files for the mmg library.
 * \author Algiane Froehly (Inria / IMB, Université de Bordeaux)
 * \author Mark Potse (IMB, Université de Bordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#ifndef MMGLIB_H
#define MMGLIB_H

#include "mmg/mmg2d/libmmg2d.h"
#include "mmg/mmgs/libmmgs.h"
#include "mmg/mmg3d/libmmg3d.h"

#endif
