## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================

## =============================================================================
##
## Installation and copy of common headers
##
## =============================================================================

SET(MMGCOMMON_SHRT_INCLUDE    mmg/common )
SET(MMGCOMMON_INCLUDE         ${PROJECT_BINARY_DIR}/include/${MMGCOMMON_SHRT_INCLUDE} )


############################################################################
#####
#####         Headers
#####
############################################################################

# common header files needed for mmg<X> libraries
#
# Remark: header installation would need to be cleaned, for now, to allow
# independent build of each project and because mmgs and mmg2d have been added
# to mmg3d without rethinking the install architecture, the header files that
# are common between codes are copied in all include directories (mmg/,
# mmg/mmg3d/, mmg/mmgs/, mmg/mmg2d/).  they are also copied in build directory
# to enable library call without installation.
SET( mmgcommon_headers
  ${MMGCOMMON_SOURCE_DIR}/mmg_export.h
  ${MMGCOMMON_SOURCE_DIR}/libmmgtypes.h
  ${MMGCOMMON_BINARY_DIR}/mmgcmakedefines.h
  ${MMGCOMMON_BINARY_DIR}/mmgcmakedefinesf.h
  ${MMGCOMMON_BINARY_DIR}/mmgversion.h
  )

IF ( PERL_FOUND )
  LIST ( APPEND mmgcommon_headers  ${MMGCOMMON_BINARY_DIR}/libmmgtypesf.h )
ENDIF ( )

IF ( MMG_INSTALL_PRIVATE_HEADERS )
  LIST ( APPEND mmgcommon_headers
    ${MMGCOMMON_SOURCE_DIR}/libmmgcommon_private.h
    ${MMGCOMMON_SOURCE_DIR}/mmgcommon_private.h
    ${MMGCOMMON_SOURCE_DIR}/mmg_core_export_private.h
    ${MMGCOMMON_SOURCE_DIR}/inlined_functions_private.h
    ${MMGCOMMON_SOURCE_DIR}/mmgexterns_private.h
    ${MMGCOMMON_SOURCE_DIR}/librnbg_private.h
    ${MMGCOMMON_SOURCE_DIR}/chrono_private.h
    ${MMGCOMMON_SOURCE_DIR}/eigenv_private.h
    ${MMGCOMMON_SOURCE_DIR}/vtkparser.hpp
    )
ENDIF()

IF (NOT WIN32 OR MINGW)
  LIST(APPEND mmgcommon_headers  ${MMGCOMMON_BINARY_DIR}/git_log_mmg.h )
ENDIF()

# Install header files in /usr/local or equivalent
INSTALL(FILES ${mmgcommon_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/common COMPONENT headers)

# Copy header files in project directory at build step
COPY_1_HEADER_AND_CREATE_TARGET(
  ${MMGCOMMON_SOURCE_DIR} libmmgtypes ${MMGCOMMON_INCLUDE} "mmgcommon")

COPY_1_HEADER_AND_CREATE_TARGET(
  ${MMGCOMMON_BINARY_DIR} mmgcmakedefines ${MMGCOMMON_INCLUDE} "mmgcommon")

COPY_1_HEADER_AND_CREATE_TARGET(
  ${MMGCOMMON_SOURCE_DIR} mmg_export
  ${MMGCOMMON_INCLUDE} "mmgcommon")

# Fortran generation
COPY_HEADER (
  ${MMGCOMMON_BINARY_DIR} mmgcmakedefinesf.h
  ${MMGCOMMON_INCLUDE} mmgcmakedefinesf.h
  mmgmmgcommon_mmgcmakedefines_header copymmgcommon_mmgcmakedefinesf )

COPY_1_HEADER_AND_CREATE_TARGET(
  ${MMGCOMMON_BINARY_DIR} mmgversion ${MMGCOMMON_INCLUDE} "mmgcommon")

SET ( tgt_list copymmgcommon_libmmgtypes copymmgcommon_mmgcmakedefines
  copymmgcommon_mmg_export copymmgcommon_mmgcmakedefinesf
  copymmgcommon_mmgversion )

if (PERL_FOUND)
  COPY_HEADER (
    ${MMGCOMMON_BINARY_DIR} libmmgtypesf.h
    ${MMGCOMMON_INCLUDE} libmmgtypesf.h
    mmg_fortran_header copymmgcommon_libmmgtypesf )

  LIST ( APPEND tgt_list copymmgcommon_libmmgtypesf )
endif()

IF (NOT WIN32 OR MINGW)
  COPY_HEADER (
    ${MMGCOMMON_BINARY_DIR} git_log_mmg.h
    ${MMGCOMMON_INCLUDE} git_log_mmg.h
    GenerateGitHashMmg copymmgcommon_mmggithash )

  LIST ( APPEND tgt_list copymmgcommon_mmggithash)
ENDIF ()

IF ( MMG_INSTALL_PRIVATE_HEADERS )
  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} libmmgcommon_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} mmgcommon_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} mmg_core_export_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} inlined_functions_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} mmgexterns_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} librnbg_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} chrono_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMGCOMMON_SOURCE_DIR} eigenv_private ${MMGCOMMON_INCLUDE} "mmgcommon")

  ADD_CUSTOM_TARGET(mmgcommon_vtkparser_header ALL
    DEPENDS
    ${MMGCOMMON_SOURCE_DIR}/vtkparser.hpp )

  COPY_HEADER (
    ${MMGCOMMON_SOURCE_DIR} vtkparser.hpp
    ${MMGCOMMON_INCLUDE} vtkparser.hpp
    mmgcommon_vtkparser_header copymmgcommon_vtkparser )


  LIST ( APPEND tgt_list copymmgcommon_libmmgcommon_private
    copymmgcommon_mmgcommon_private copymmgcommon_mmg_core_export_private
    copymmgcommon_inlined_functions_private copymmgcommon_mmgexterns_private
    copymmgcommon_librnbg_private copymmgcommon_chrono_private
    copymmgcommon_eigenv_private copymmgcommon_vtkparser)

ENDIF()

ADD_CUSTOM_TARGET (copy_mmgcommon_headers ALL
    DEPENDS ${tgt_list} )
