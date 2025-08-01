/* ============================================================================ =
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright(c) Bx INP / CNRS / Inria / UBordeaux / UPMC, 2004 -
**
**  mmg is free software : you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE.See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg(in
**  files COPYING.LESSER and COPYING).If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** ============================================================================ =
*/

#ifndef MMG_EXPORT_H
#define MMG_EXPORT_H

#include "mmg/common/mmgcmakedefines.h"

#if(MMG_DYN_LIB)
#  if defined(_WIN32) || defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) || defined(__TOS_WIN__) || defined(__WINDOWS__)
#    define MMG_DECL_EXPORT     __declspec(dllexport)
#    define MMG_DECL_IMPORT     __declspec(dllimport)
#    define MMG_DECL_HIDDEN
#  elif(__GNUC__ >= 4)
#    define MMG_DECL_EXPORT     __attribute__((visibility("default")))
#    define MMG_DECL_IMPORT     __attribute__((visibility("default")))
#    define MMG_DECL_HIDDEN     __attribute__((visibility("hidden")))
#  endif
#endif

#if !defined(MMG_DECL_EXPORT)
#  define MMG_DECL_EXPORT
#  define MMG_DECL_IMPORT
#  define MMG_DECL_HIDDEN
#endif

#endif
